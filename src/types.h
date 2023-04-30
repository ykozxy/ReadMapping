#ifndef READMAPPING_TYPES_H
#define READMAPPING_TYPES_H

#include <utility>
#include <vector>
#include <string>
#include <iostream>
#include <set>

#include "logging.h"

#define SUBSTR_MIN_LENGTH 15  // Substrings less than this length will not be considered
#define SUBSTR_COMBINE_LENGTH 8  // Substrings less than this length will be combined
#define SUBSTR_THRESHOLD 40  // Substrings less than this length will only be split two-fold
#define SUBSTR_DIVISIONS {3, 4, 5, 6}  // The number of divisions to make when splitting the read
#define MAX_GAP 25  // Maximum indices to search backward/forward in local search per position
#define LOCAL_MATCH_SCORE_TOLERANCE 4  // Paths with difference with top score less than this will be considered in backtrace in local match
#define GLOBAL_MATCH_MAX_LENGTH 100000  // Maximum length of a global match
#define NUM_THREADS 8  // Number of threads to use for parallelization
#define TIMEOUT 5  // Timeout in seconds for each read
#define READMAP_PADDING_THRESHOLD 10  // Maximum padding adds to the left of a read when splitting it


/* Enum class of Insertion/Deletion (Indel) */
enum class IndelType {
    INSERTION,
    DELETION
};

/* Represents an indel */
struct Indel {
    IndelType type;
    int position;
    int length;

    Indel(IndelType type, int position, int length) {
        this->type = type;
        this->position = position;
        this->length = length;
    }

    bool operator==(const Indel &other) const {
        return position == other.position;
    }

    bool operator<(const Indel &other) const {
        return position < other.position;
    }

    [[nodiscard]] bool equals(const Indel &other) const {
        return type == other.type && position == other.position && length == other.length;
    }

    [[nodiscard]] std::string to_string() const {
        std::string result = "Indel{";
        if (type == IndelType::INSERTION) {
            result += "I";
        } else {
            result += "D";
        }
        result += "@" + std::to_string(position) + ", l=" + std::to_string(length) + "}";
        return result;
    }
};

/* Represents a mismatch */
struct Mismatch {
    int position;
    char new_base;

    Mismatch(int position, char new_base) {
        this->position = position;
        this->new_base = new_base;
    }

    bool operator==(const Mismatch &other) const {
        return position == other.position;
    }

    bool operator<(const Mismatch &other) const {
        return position < other.position;
    }

    [[nodiscard]] bool equals(const Mismatch &other) const {
        return position == other.position && new_base == other.new_base;
    }
};

/* Enum class of Mutation (Mismatch and Indel) */
enum class MutationType {
    INSERTION,
    DELETION,
    MISMATCH
};

/* Represents a mutation */
struct Mutation {
    MutationType type;
    int position;
    char new_base{};

    Mutation(MutationType type, int position) {
        if (type == MutationType::MISMATCH || type == MutationType::INSERTION) {
            throw std::invalid_argument("Mismatch or insertion must have a new base");
        }
        this->type = type;
        this->position = position;
    }

    Mutation(MutationType type, int position, char new_base) {
        this->type = type;
        this->position = position;
        this->new_base = new_base;
    }

    bool operator==(const Mutation &other) const {
        return position == other.position;
    }

    bool operator<(const Mutation &other) const {
        return position < other.position;
    }

    [[nodiscard]] bool equals(const Mutation &other) const {
        if (type != other.type) {
            return false;
        }
        switch (type) {
            case MutationType::INSERTION:
            case MutationType::MISMATCH:
                return position == other.position && new_base == other.new_base;
            case MutationType::DELETION:
                return position == other.position;
        }
    }

    [[nodiscard]] std::string to_string() const {
        std::string result = "Mutation{";
        if (type == MutationType::INSERTION) {
            result += "I";
        } else if (type == MutationType::DELETION) {
            result += "D";
        } else {
            result += "M";
        }
        result += "@" + std::to_string(position);
        if (type == MutationType::MISMATCH || type == MutationType::INSERTION) {
            result += ", new=" + std::string(1, new_base);
        }
        result += "}";
        return result;
    }
};

/* Represents a mapping of a read sequence to a genome */
struct Mapping {
    int position;
    std::set<Mismatch> mismatches;
    std::multiset<Indel> indels;

    Mapping() : Mapping(-1) {}

    explicit Mapping(int position) {
        this->position = position;
    }

    void merge(Mapping *other, Logger *logger) {
        if (other->position < 0) {
            logger->log_msg(LogLevel::WARN, "Mapping::merge", "Cannot merge mapping with negative position.");
            return;
        }

        if (position < 0) {
            position = other->position;
        } else {
            position = std::min(position, other->position);
        }

        mismatches.insert(other->mismatches.begin(), other->mismatches.end());
        for (auto indel: other->indels) {
            add_indel(indel, logger);
        }
    }

    /**
     * Returns the number of mismatches and indels
     * @return number of mismatches and indels
     */
    [[nodiscard]] int num_errors() const {
        int score = (int) mismatches.size();
        for (auto indel: indels) {
            score += indel.length;
        }
        return score;
    }

    /**
     * Returns the adjustment to the length of the reference sequence due to indels.
     * @return
     */
    [[nodiscard]] int length_adjust() const {
        int adjust = 0;
        for (auto indel: indels) {
            if (indel.type == IndelType::INSERTION) {
                adjust -= indel.length;
            } else {
                adjust += indel.length;
            }
        }
        return adjust;
    }

    /**
     * Adds an indel to the mapping
     * @param indel indel to add
     */
    void add_indel(Indel indel, Logger *logger) {
        // look for indel at the same pos:
        Indel *indel_ptr = nullptr;
        for (const auto &i: indels) {
            if (i.position == indel.position) {
                // if indel is the same type, add lengths
                if (i.type == indel.type) {
                    indel_ptr = const_cast<Indel *>(&i);
                    indel_ptr->length += indel.length;
                } else {
                    // if indel is different type, throw error
//                    throw std::runtime_error("Indels at the same position have different type");
                    logger->log_msg(LogLevel::ERROR, "Mapping::add_indel",
                                    "Indels at same position have different type");
                }
                return;
            }
        }
        // if no indel at same pos, add indel
        indels.insert(indel);
    }

    /**
     * Visualize the mapping by providing a reference and read sequence.
     * @param reference - The reference sequence
     * @param read      - The read sequence
     */
    void visulize(std::string reference, std::string read) {
        int MAX_PADDING = 10;

        std::string modified_reference = reference.substr(0, position);
        for (int i = position; i < (int) reference.length(); i++) {
            for (auto indel: indels) {
                if (indel.position == i) {
                    if (indel.type == IndelType::INSERTION) {
                        for (int j = 0; j < indel.length; j++) {
                            modified_reference += "-";
                        }
                    }
                }
            }
            modified_reference += reference[i];
        }

        std::string modified_read;
        int del_length = 0;
        for (int i = 0; i < (int) read.length(); i++) {
            for (auto indel: indels) {
                if (indel.position - position - del_length == i) {
                    if (indel.type == IndelType::DELETION) {
                        for (int j = 0; j < indel.length; j++) {
                            modified_read += "-";
                        }
                        del_length += indel.length;
                    }
                }
            }
            modified_read += read[i];
        }

        std::string read_padding, reference_padding;
        for (int i = 0; i < position; i++) {
            read_padding += " ";
        }
        if (position > MAX_PADDING) {
            read_padding = "   " + read_padding.substr(position - MAX_PADDING);
            reference_padding = "...";
            modified_reference = modified_reference.substr(position - MAX_PADDING);
        }
        if ((int) (modified_reference.length() - modified_read.length()) > MAX_PADDING * 2) {
            modified_reference = modified_reference.substr(0, MAX_PADDING * 2 + modified_read.size()) + "...";
        }
        std::cout << "Ref : " << reference_padding << modified_reference << std::endl;
        std::cout << "Read: " << read_padding << modified_read << std::endl;
    }

    [[nodiscard]] std::string to_string() const {
        // format: Mapping{pos=position, mis={mismatch1, mismatch2}, indels={indel1, indel2, ...}}
        std::string result = "Mapping{pos=" + std::to_string(position) + ", mis={";
        for (auto mismatch: mismatches) {
            result += std::string(1, mismatch.new_base) + "@" + std::to_string(mismatch.position) + ", ";
        }
        if (!mismatches.empty())
            result = result.substr(0, result.length() - 2);
        result += "}, indels={";
        for (auto indel: indels) {
            result += indel.to_string() + ", ";
        }
        if (!indels.empty())
            result = result.substr(0, result.length() - 2);
        result += "}}";
        return result;
    }
};

/* Represents a read sequence */
struct ReadSequence {
    std::string read; // The read sequence
    Mapping *mapping; // The mapping of the read sequence to a genome

    explicit ReadSequence(std::string read) {
        this->read = std::move(read);
        this->mapping = nullptr;
    }

    ~ReadSequence() {
        delete mapping;
    }

    [[nodiscard]] std::string to_string() const {
        return "ReadSequence{read=" + read + ", mapping=" + mapping->to_string() + "}";
    }
};

/* Represents a pair of read sequences */
struct ReadPair {
    ReadSequence *left;
    ReadSequence *right;

    ReadPair(ReadSequence *left, ReadSequence *right) {
        this->left = left;
        this->right = right;
    }

    ~ReadPair() {
        delete left;
        delete right;
    }

    [[nodiscard]] std::string to_string() const {
        return "ReadPair{left=" + left->read + ", right=" + right->read + "}";
    }
};

#endif //READMAPPING_TYPES_H

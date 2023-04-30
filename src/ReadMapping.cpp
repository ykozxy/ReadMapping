#include "ReadMapping.h"

#include <iostream>
#include <fstream>
#include <future>

#include "logging.h"

using namespace std;

ReadMapping::ReadMapping(const string &reference_file, const string &reads_file, Logger *logger) : cache(nullptr) {
    this->logger = logger;
    read_reference(reference_file);
    read_reads(reads_file);
    logger->log_msg(LogLevel::INFO, "ReadMapping", "Building index...");
    build_index();
}

ReadMapping::~ReadMapping() {
    for (auto *read: reads) {
        delete read;
    }

    for (auto p: all_mutations) {
        delete p.second;
    }
}

void ReadMapping::do_mapping() {
    logger->create_progress_bar((int) reads.size());

//    map_worker(0, (int) reads.size());
//    logger->finish_progress_bar();
//    return;

    vector<future<void>> futures;

    // Calculate number reads per worker
    int read_per_worker = (int) reads.size() / NUM_THREADS;
    vector<int> num_reads;
    num_reads.reserve(NUM_THREADS);
    for (int i = 0; i < NUM_THREADS; i++) {
        num_reads.push_back(read_per_worker);
    }

    // Distribute the remaining reads
    int remaining_reads = (int) reads.size() % NUM_THREADS;
    for (int i = 0; i < remaining_reads; i++) {
        num_reads[i]++;
    }

    // Create workers
    int start = 0;
    for (int i = 0; i < NUM_THREADS; i++) {
        futures.push_back(async(launch::async, &ReadMapping::map_worker, this, start, start + num_reads[i]));
        start += num_reads[i];
    }

    // Join workers
    for (auto &f: futures) {
        f.wait();
    }
    logger->finish_progress_bar();
}

void ReadMapping::map_worker(int start, int end) {
    for (int i = start; i < end; i++) {
        try {
            auto as = async(launch::async, &ReadMapping::map_pair, this, reads[i]);

            // Wait for mapping to finish
            auto status = as.wait_for(chrono::seconds(TIMEOUT));
            switch (status) {
                case future_status::ready:
                    break;
                case future_status::timeout:
                    logger->log_msg(LogLevel::ERROR, "map_worker", "i=" + to_string(i) + ", Mapping timed out.");
                    continue;
                case future_status::deferred:
                    logger->log_msg(LogLevel::ERROR, "map_worker", "i=" + to_string(i) + ", Mapping deferred.");
                    continue;
            }

            logger->update_progress_bar();
        } catch (exception &e) {
            logger->log_msg(LogLevel::ERROR, "map_worker", "i=" + to_string(i) + ", " + e.what());
            continue;
        }
    }
}

void ReadMapping::map_pair(ReadPair *read_pair) {
    vector<Mapping *> mappings_left;
    for (int i = 0; i <= READMAP_PADDING_THRESHOLD; i++) {
        map_read(read_pair->left, mappings_left, i);
        if (!mappings_left.empty()) {
            break;
        }
    }

//    map_read(read_pair->left, mappings_left);
    if (mappings_left.empty()) {
        logger->log_msg(LogLevel::WARN, "map_pair", "Left pair (" + read_pair->left->read + ") has no mapping.");
//        return;
    }

    vector<Mapping *> mappings_right;
    for (int i = 0; i <= READMAP_PADDING_THRESHOLD; i++) {
        map_read(read_pair->right, mappings_right, i);
        if (!mappings_right.empty()) {
            break;
        }
    }

//    map_read(read_pair->right, mappings_right);
    if (mappings_right.empty()) {
        logger->log_msg(LogLevel::WARN, "map_pair", "Right pair (" + read_pair->right->read + ") has no mapping.");
        if (!mappings_left.empty()) {
            read_pair->left->mapping = mappings_left[0];
        }
        return;
    }
    if (mappings_left.empty()) {
        read_pair->right->mapping = mappings_right[0];
        return;
    }

    // Find one combination of mappings with the smallest positive distance
    int min_distance = INT32_MAX;
    Mapping *mapping_left = nullptr;
    Mapping *mapping_right = nullptr;
    for (auto *mapping_l: mappings_left) {
        for (auto *mapping_r: mappings_right) {
            int distance = mapping_r->position - (mapping_l->position + mapping_l->length_adjust());
            if (distance > 0 && distance < min_distance) {
                min_distance = distance;
                mapping_left = mapping_l;
                mapping_right = mapping_r;
            }
        }
    }

    if (mapping_left == nullptr || mapping_right == nullptr) {
        logger->log_msg(LogLevel::ERROR, "map_pair", "No combination of mappings found.");
        logger->log_msg(LogLevel::ERROR, "map_pair", "\t left: " + read_pair->left->read);
        logger->log_msg(LogLevel::ERROR, "map_pair", "\t right: " + read_pair->right->read);
        return;
    }

    read_pair->left->mapping = mapping_left;
    read_pair->right->mapping = mapping_right;

//    lock_guard<mutex> lock(log_mutex);
//    read_pair->left->mapping->visulize(m_string, read_pair->left->read);
//    read_pair->right->mapping->visulize(m_string, read_pair->right->read);
//    cout << endl;
}

void ReadMapping::map_read(ReadSequence *read, std::vector<Mapping *> &mappings, int padding) const {
    /* If read is too short, abort */
    if (read->read.length() <= SUBSTR_MIN_LENGTH) {
        return;
    }

    /* For short reads, divide it into longer substrings */
    vector<int> divisions = SUBSTR_DIVISIONS;
    if (read->read.length() - padding <= SUBSTR_THRESHOLD) {
        divisions = {2};
    }

    // from longer to shorter substrings
    for (auto div: divisions) {
        string read_copy = read->read.substr(padding, read->read.length() - padding);
        size_t length = read_copy.length() / div;

        // split read into substrings
        vector<string> substrings;
        size_t cumulate_length = 0;
        while (cumulate_length < read_copy.length()) {
            substrings.emplace_back(
                    read_copy.substr(cumulate_length, min(length, (int) read_copy.length() - cumulate_length))
            );
            cumulate_length += length;
        }

        // Add the padding to the first substring
        substrings[0] = read->read.substr(0, padding) + substrings[0];

        // if the last substring is too short, combine it with the previous substring
        if ((int) substrings.back().length() < SUBSTR_COMBINE_LENGTH) {
            substrings[substrings.size() - 2] += substrings.back();
            substrings.pop_back();
        }

        // find all exact matches
        vector<vector<int>> exact_matches;
        for (const auto &substring: substrings) {
            vector<int> res;
            if ((res = cache->find(substring)).empty()) {
                exact_matches.emplace_back();
            } else {
                exact_matches.emplace_back(res);
            }
        }

        // sort all exact matches
        for (auto &match: exact_matches) {
            sort(match.begin(), match.end());
        }

        // find all possible combinations of match positions
        vector<vector<int>> combinations;
        get_increasing_sequences(exact_matches, combinations);

        // Error check: if no combination found, skip
        if (combinations.empty()) {
            logger->log_msg(LogLevel::DEBUG, "map_read", "No combination found for read " + read->read);
            continue;
        }

        // Skip if no success match
        int success_match_count = 0;
        for (auto &p: combinations[0]) {
            if (p >= 0) {
                success_match_count++;
            }
        }
        if (success_match_count == 0) {
            logger->log_msg(LogLevel::DEBUG, "map_read",
                            "No success match (div=" + to_string(div) + ") for read " + read->read);
            continue;
        }
        if (success_match_count < ceil((double) combinations.size() / 2.0)) {
            logger->log_msg(LogLevel::DEBUG, "map_read",
                            "Less than half success match (div=" + to_string(div) + ") for read " + read->read);
        }

        // check each combination
        for (auto &combination: combinations) {
            // error check
            if (combination.empty()) {
                logger->log_msg(LogLevel::WARN, "map_read", "Empty combination found for read " + read->read);
                continue;
            }

            // Sanity check: if the gap between two matches is too large, skip
            int max_gap = 0;
            int last_pos = -1;
            for (auto &p: combination) {
                if (p >= 0) {
                    if (last_pos >= 0) {
                        max_gap = max(max_gap, p - last_pos);
                    }
                    last_pos = p;
                }
            }
            if (max_gap > MAX_GAP) {
                logger->log_msg(LogLevel::DEBUG, "map_read",
                                "Too large gap (div=" + to_string(div) + ") for read " + read->read);
                continue;
            }

            // merging mismatches
            vector<string> merged_substrings;
            for (size_t i = 0; i < substrings.size(); i++) {
                if (i == 0) {
                    merged_substrings.emplace_back(substrings[i]);
                } else if (combination[i] < 0 && combination[i - 1] < 0) {
                    merged_substrings.back() += substrings[i];
                    combination.erase(combination.begin() + (int) i);
                    substrings.erase(substrings.begin() + (int) i);
                    i--;
                } else {
                    merged_substrings.emplace_back(substrings[i]);
                }
            }

            bool uncertain;
            auto m = match_sequence(combination, merged_substrings, read, uncertain);
            if (uncertain || m == nullptr) {
                if (m == nullptr) {
                    logger->log_msg(LogLevel::DEBUG, "map_read", "Quick match failed for read " + read->read);
                } else {
                    logger->log_msg(LogLevel::DEBUG, "map_read", "Quick match is uncertain about read " + read->read);

                }

                auto m2 = match_sequence_accurate(combination, merged_substrings, read);
                if (m == nullptr) {
                    m = m2;
                } else if (m2 != nullptr) {
                    // Compare m with m2
                    m = m->num_errors() < m2->num_errors() ? m : m2;
                }

//                if (m2 == nullptr || m2->num_errors() >= 8) {
//                    logger->log_msg(LogLevel::WARN, "map_read", "Quick and accurate match failed for read " + read->read);
//                    auto m3 = global_match(read);
//
//                    if (m3 != nullptr && m3->num_errors() < 8) {
//                        m = m3;
//                    }
//                }

                if (m2 != nullptr && m2->num_errors() >= 8) {
                    logger->log_msg(LogLevel::DEBUG, "map_read",
                                    "Too many mismatch (div=" + to_string(div) + ") for read " + read->read);
                    continue;
                }
            }

            if (m == nullptr) {
                continue;
            }
            mappings.emplace_back(m);
        }

        // Stop once we have found a mapping
        if (!mappings.empty()) {
            return;
        }
    }

    if (mappings.empty()) {
        // Use global match
        // FIXME: this is extremely slow
//        auto m = global_match(read);
//        if (m != nullptr && m->num_errors() < 8) {
//            mappings.emplace_back(m);
//        }
    }
}

Mapping *
ReadMapping::match_sequence(const vector<int> &sequence,
                            const vector<string> &substrings,
                            const ReadSequence *read,
                            bool &uncertain) const {
    uncertain = false;
    vector<Mapping *> all_maps;

    // handle first substring mismatch
    if (sequence[0] < 0) {
        int search_left_bound = max(0,
                                    sequence[1] - average_read_length - (int) substrings[0].length());
        int search_right_bound = sequence[1];
        auto best_map = local_match(
                reference_dna,
                substrings[0],
                search_left_bound,
                search_right_bound, logger);
        if (best_map == nullptr) {
            logger->log_msg(LogLevel::WARN, "map_sequence", "No mapping for segment 0 of read " + read->read);
            return nullptr;
        }
        all_maps.emplace_back(best_map);
    } else {
        all_maps.emplace_back(new Mapping(sequence[0]));
    }

    // handle middle substring mismatch
    for (size_t i = 1; i < substrings.size() - 1; i++) {
        if (sequence[i] < 0) {
            // mismatch found!
            // left bound is the end of the previous substring
            int search_left_bound = (int) substrings[i - 1].length() + sequence[i - 1];
            // right bound is the start of the next substring
            int search_right_bound = sequence[i + 1];
            auto best_map = local_match(
                    reference_dna,
                    substrings[i],
                    search_left_bound,
                    search_right_bound, logger);
            if (best_map == nullptr) {
                logger->log_msg(LogLevel::WARN, "map_sequence",
                                "No mapping for segment " + to_string(i) + " of read " + read->read);
                continue;
            }
            all_maps.emplace_back(best_map);
        } else {
            // exact match!
            // push a new mapping without mismatch or indel
            all_maps.push_back(new Mapping(sequence[i]));
        }
    }

    // handle last substring mismatch
    if (sequence.back() < 0) {
        int search_left_bound = sequence[sequence.size() - 2] +
                                (int) substrings[substrings.size() - 2].length();
        int search_right_bound = min((int) reference_dna.length(),
                                     search_left_bound + average_read_length);
        auto best_map = local_match(
                reference_dna,
                substrings.back(),
                search_left_bound,
                search_right_bound, logger);
        if (best_map == nullptr) {
            logger->log_msg(LogLevel::WARN, "map_sequence",
                            "No mapping for segment " + to_string(sequence.size() - 1) + " of read " + read->read);
            return nullptr;
        }
        all_maps.emplace_back(best_map);
    } else {
        all_maps.emplace_back(new Mapping(sequence.back()));
    }

    // Merge all mappings
    auto *merged_map = new Mapping();
    int end_pos = all_maps[0]->position;
    for (size_t i = 0; i < all_maps.size(); ++i) {
        auto map = all_maps[i];

        // Sanity check
        if (end_pos != map->position) {
            logger->log_msg(LogLevel::DEBUG, "map_sequence",
                            "Merged mapping is not continuous for read " + read->read + ", attempting to fix");

            // Try to fix the problem
            int indel_length;
            if (map->position > end_pos) {
                // Add deletions
                indel_length = map->position - end_pos;
                map->add_indel({IndelType::DELETION, end_pos + 2, indel_length}, logger);
            } else {
                // Add insertions
                indel_length = end_pos - map->position;
                map->add_indel({IndelType::INSERTION, end_pos + 2, indel_length}, logger);
            }

            // Prevent too-long indels
            if (indel_length >= 5) {
                logger->log_msg(LogLevel::DEBUG, "map_sequence",
                                "Inter-map indel length too long for read " + read->read + ", skipping.");
                return nullptr;
            }

            uncertain = true;
        }

        // Check if we have a mismatch / indel at the start or end of substring
        // FIXME: check uncertain condition position
        int new_end_pos = end_pos + (int) substrings[i].length() + map->length_adjust();
        if (!uncertain) {
            for (auto indel: map->indels) {
                if (indel.position <= end_pos + 1 || indel.position >= new_end_pos - 1) {
                    logger->log_msg(LogLevel::DEBUG, "map_sequence",
                                    "Indels close to substr boundaries for read " + read->read);
                    uncertain = true;
                    break;
                }
            }
            for (auto mismatch: map->mismatches) {
                if (mismatch.position <= end_pos + 2 || mismatch.position >= new_end_pos - 2) {
                    logger->log_msg(LogLevel::DEBUG, "map_sequence",
                                    "Mismatches close to substr boundaries for read " + read->read);
                    uncertain = true;
                    break;
                }
            }
        }

        merged_map->merge(map, logger);

        end_pos = new_end_pos;
    }

    return merged_map;
}

Mapping *
ReadMapping::match_sequence_accurate(const vector<int> &sequence,
                                     const vector<string> &substrings,
                                     const ReadSequence *read) const {
    // Find the left-most and right-most range
    int left_bound = -1;
    int right_bound = -1;
    for (int i: sequence) {
        if (i >= 0) {
            left_bound = i;
            break;

        }
    }
    for (int i = (int) sequence.size() - 1; i >= 0; i--) {
        if (sequence[i] >= 0) {
            right_bound = sequence[i] + (int) substrings[i].length();
            break;
        }
    }

    left_bound = max(0, left_bound - average_read_length);
    right_bound = min((int) reference_dna.length() - 1, right_bound + average_read_length);

    return local_match(reference_dna, read->read, left_bound, right_bound, logger, true);
}


Mapping *
ReadMapping::global_match(const ReadSequence *read) const {
    if (reference_dna.length() <= GLOBAL_MATCH_MAX_LENGTH) {
        return local_match(reference_dna, read->read, 0, (int) reference_dna.length(), logger);
    }

    // Split the reference DNA into chunks
    mutex m;
    Mapping *best_mapping = nullptr;

    auto worker = [&](int left_bound, int right_bound) {
        auto *mapping = local_match(reference_dna, read->read, left_bound, right_bound, logger);
        lock_guard<mutex> lock(m);
        if (best_mapping == nullptr) {
            best_mapping = mapping;
        } else if (mapping != nullptr) {
            best_mapping = best_mapping->num_errors() < mapping->num_errors() ? best_mapping : mapping;
        }
    };
    vector<future<void>> futures;

    for (int i = 0; i < (int) reference_dna.length(); i += GLOBAL_MATCH_MAX_LENGTH - average_read_length * 2) {
        int left_bound = i;
        int right_bound = min(i + GLOBAL_MATCH_MAX_LENGTH, (int) reference_dna.length());
        futures.emplace_back(async(launch::async, worker, left_bound, right_bound));
    }

    for (auto &f: futures) {
        f.wait();
    }
    return best_mapping;
}


void ReadMapping::calculate_mutations(bool verbose) {
    all_mutations.clear();

    // Construct map for start and end positions of read sequences
    map<unsigned int, vector<ReadSequence *>> start_pos;
    map<unsigned int, vector<ReadSequence *>> end_pos;

    for (auto *read: reads) {
        for (auto r: {read->left, read->right}) {
            auto mapping = r->mapping;
            if (mapping == nullptr) {
                continue;
            }

            auto start = mapping->position;
            auto end = start + mapping->length_adjust() + r->read.length();
            start_pos[start].push_back(r);
            end_pos[end].push_back(r);
        }
    }

    // Iterate over all positions in the reference DNA
    vector<ReadSequence *> active_sequence;
    logger->create_progress_bar((int) reference_dna.size());
    long double coverage = 0;
    for (int i = 0; i < (int) reference_dna.length(); i++) {
        logger->update_progress_bar();

        // Remove mappings that end at this position
        if (end_pos.find(i) != end_pos.end()) {
            for (auto *read: end_pos[i]) {
                active_sequence.erase(remove(active_sequence.begin(), active_sequence.end(), read),
                                      active_sequence.end());
            }
        }

        // Add new sequences to active sequence
        if (start_pos.find(i) != start_pos.end()) {
            for (auto *sequence: start_pos[i]) {
                active_sequence.push_back(sequence);
            }
        }
        coverage = (coverage * i + (double) active_sequence.size()) / (i + 1);

        // Calculate mutations at this position
        vector<pair<Mutation *, int>> mutation_count;
        for (auto read: active_sequence) {
            int read_pos = 0;
            int ref_pos = 0;
            vector<Mutation *> mut;

            Mutation *m = nullptr;

            // Check indels at this position
            for (auto indel: read->mapping->indels) {
                if (indel.type == IndelType::DELETION) {
                    read_pos += indel.length;
                } else {
                    ref_pos += indel.length;
                }

                if (indel.type == IndelType::INSERTION &&
                    i + 1 >= indel.position &&
                    i + 1 < indel.position + indel.length) {
                    try {
                        m = new Mutation(
                                MutationType::INSERTION,
                                i,  // FIXME: in output, the insert is AFTER the position
                                read->read.at(i - read->mapping->position + read_pos + 1));
                        mut.push_back(m);
                    } catch (out_of_range &e) {
                        logger->log_msg(LogLevel::ERROR, "calculate_mutations",
                                        "Out of range error for read " + read->read + " at position " + to_string(i));
                        break;
                    }
                    if (verbose) cout << m->to_string() << endl;
                    break;
                } else if (indel.type == IndelType::DELETION &&
                           i >= indel.position &&
                           i < indel.position + indel.length) {
                    m = new Mutation(MutationType::DELETION, i, reference_dna.at(i + ref_pos));
                    mut.push_back(m);
                    if (verbose) cout << m->to_string() << endl;
                    break;
                }
            }

            // Check mismatches at this position
            for (auto mismatch: read->mapping->mismatches) {
                if (i == mismatch.position) {
                    auto new_m = new Mutation(MutationType::MISMATCH, i, mismatch.new_base);
                    mut.push_back(new_m);
                    if (verbose) std::cout << new_m->to_string() << std::endl;

                    if (m != nullptr) {
                        logger->log_msg(LogLevel::WARN, "calculate_mutations",
                                        "Multiple mutations at position " + to_string(i) + " in read " + read->read);
                        break;
                    }
                    m = new_m;
                    break;
                }
            }

            if (m == nullptr) {
                continue;
            }

            // Check if mutation already exists
            for (auto t: mut) {
                bool found = false;
                for (auto &p: mutation_count) {
                    if (p.first->equals(*t)) {
                        p.second++;
                        delete t;
                        t = nullptr;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    mutation_count.emplace_back(t, 1);
                }
            }
        }

        if (mutation_count.empty()) {
            continue;
        }

        // Count mutations
        int sum_count = 0;
        int max_count = 0;
        Mutation *max_mutation = nullptr;
        for (auto &p: mutation_count) {
            sum_count += p.second;
            if (p.second > max_count) {
                max_count = p.second;
                max_mutation = p.first;
            }
        }

        // Get majority mutation
        // FIXME: adjust this bound to see if setting a min value (3) produces better results
        if (max_count > ((int) active_sequence.size()) / 2) {
            // Majority of reads have the same mutation at this position -> insert into map
            all_mutations[i] = max_mutation;
        } else {
            // No majority -> delete mutation
            logger->log_msg(LogLevel::DEBUG, "calculate_mutations",
                            "No majority mutation at position " + to_string(i));
            delete max_mutation;
            continue;
        }

        // Delete all other mutations
        for (auto &p: mutation_count) {
            if (p.first != max_mutation) {
                delete p.first;
            }
        }
    }
    logger->finish_progress_bar();

//    coverage /= (int) reference_dna.length();
    logger->log_msg(LogLevel::INFO, "calculate_mutations", "Average coverage: " + to_string(coverage));
}

void ReadMapping::output_mutations(const string &output_file) {
    ofstream out(output_file);
    if (!out.is_open()) {
        logger->log_msg(LogLevel::ERROR, "output_mutations", "Could not open output file " + output_file);
        return;
    }

    for (auto &p: all_mutations) {
        auto pos = p.first;
        auto mut = p.second;

        out << ">";
        switch (mut->type) {
            case MutationType::INSERTION:
                out << "I";
                break;
            case MutationType::DELETION:
                out << "D";
                break;
            case MutationType::MISMATCH:
                out << "S";
                break;
        }
        out << mut->position << " ";

        switch (mut->type) {
            case MutationType::INSERTION:
                out << mut->new_base;
                break;
            case MutationType::DELETION:
                out << reference_dna.at(pos);
                break;
            case MutationType::MISMATCH:
                out << reference_dna.at(pos) << " " << mut->new_base;
                break;
        }
        out << endl;
    }

    out.close();
}

const string &ReadMapping::get_reference_dna() const {
    return reference_dna;
}

const vector<ReadPair *> &ReadMapping::get_reads() const {
    return reads;
}

void ReadMapping::build_index() {
    average_read_length = 0;
    for (auto *read: reads) {
        average_read_length += (int) read->left->read.length();
        average_read_length += (int) read->right->read.length();
    }
    average_read_length /= (int) reads.size() * 2;

    cache = new Cache(reference_dna);
    cache->build_suffix_array();
}

void ReadMapping::read_reference(const string &input_file) {
    string input_dna;
    ifstream input_file_stream(input_file);
    if (input_file_stream.is_open()) {
        // Skip first line
        string line;
        getline(input_file_stream, line);

        // Read the rest of the file
        while (getline(input_file_stream, line)) {
            input_dna += line;
        }
    } else {
        logger->log_msg(LogLevel::ERROR, "read_reference", "Unable to open file " + input_file);
        perror("ReadMapping::read_reference");
        ::exit(1);
    }
    reference_dna = input_dna;
}

void ReadMapping::read_reads(const string &input_file) {
    ifstream input_file_stream(input_file);
    if (input_file_stream.is_open()) {
        string name;

        ReadSequence *last = nullptr;
        while (getline(input_file_stream, name)) {
            string line;
            getline(input_file_stream, line);

            if (last == nullptr) {
                last = new ReadSequence(line);
            } else {
                auto *new_read = new ReadSequence(line);
                auto *read_pair = new ReadPair(last, new_read);
                reads.push_back(read_pair);
                last = nullptr;
            }
        }
    } else {
        logger->log_msg(LogLevel::ERROR, "read_reads", "Unable to open file " + input_file);
        perror("ReadMapping::read_reads");
        ::exit(1);
    }
}


#include "algorithms.h"
#include <iostream>

#include "logging.h"

using namespace std;

using coord = pair<int, int>;

//void get_increasing_sequences(const vector<vector<int>> &sequences, vector<vector<int>> &result, int pos) {
//    if (pos == 0) {
//        // Initialize result vector
//        if (sequences[0].empty()) {
//            result.emplace_back(1, -1);
//        } else {
//            for (auto p: sequences[0]) {
//                result.emplace_back(1, p);
//            }
//        }
//
//        get_increasing_sequences(sequences, result, pos + 1);
//        return;
//    }
//
//    if (pos == (int) sequences.size()) {
//        // End case: all positions have been processed
//        return;
//    }
//
//    size_t orig_size = result.size();
//    for (size_t i = 0; i < orig_size; i++) {
//        auto &combination = result[i];
//
//        // Calculate the last greatest value in the combination
//        int last_greatest_val = -1;
//        for (int j = pos - 1; j >= 0; j--) {
//            if (combination[j] >= 0) {
//                last_greatest_val = combination[j];
//                break;
//            }
//        }
//
//        if (sequences[pos].empty()) {
//            vector<int> combination_copy = combination;
//            combination_copy.emplace_back(-1);
//            result.emplace_back(combination_copy);
//            continue;
//        }
//
//        for (auto p: sequences[pos]) {
//            if (last_greatest_val < 0 || last_greatest_val < p) {
//                vector<int> combination_copy = combination;
//                combination_copy.emplace_back(p);
//                result.emplace_back(combination_copy);
//            }
//        }
//    }
//
//    result.erase(result.begin(), result.begin() + (int) orig_size);
//    get_increasing_sequences(sequences, result, pos + 1);
//}

// Find the index of the first element in v that is greater than or equal to val
int binary_search(const vector<int> &v, int val) {
    if (v.back() < val) {
        return v.size();
    }

    int l = 0, r = v.size() - 1;
    while (l < r) {
        int mid = (l + r) / 2;
        if (v[mid] < val) {
            l = mid + 1;
        } else {
            r = mid;
        }
    }

    return l;
}

void generateCombinationsHelper(const vector<vector<int>> &sequences,
                                vector<vector<int>> &result,
                                vector<int> &current,
                                int index, int last_index = -1) {
    if (index == sequences.size()) {
        result.push_back(current);
        return;
    }

    if (sequences[index].empty()) {
        current[index] = -1;
        generateCombinationsHelper(sequences, result, current, index + 1, last_index);
    } else {
        int start_pos = binary_search(sequences[index],
                                      (index == 0 || last_index == -1) ? 0 : current[last_index] + 1);
        for (int i = start_pos; i < sequences[index].size(); ++i) {
            if (last_index != -1 && current[last_index] + MAX_GAP < sequences[index][i]) {
                break;
            }

            current[index] = sequences[index][i];
            generateCombinationsHelper(sequences, result, current, index + 1, index);
        }
    }
}


void get_increasing_sequences(const vector<vector<int>> &sequences, vector<vector<int>> &result, int pos) {
    vector<int> current(sequences.size());
    generateCombinationsHelper(sequences, result, current, 0);
}


Mapping *
local_match(const std::string &reference, const std::string &read, int left_bound, int right_bound, Logger *logger,
            int max_gap,
            bool verbose) {
    int n = (int) right_bound - left_bound + 1;
    int m = (int) read.size();

    if (n < m) {
        logger->log_msg(LogLevel::DEBUG, "local_match",
                        "Reference (" + reference.substr(left_bound, n) + ") is shorter than read (" + read + ").");
//        return nullptr;
    }

    vector<vector<int>> score(n + 1, vector<int>(m + 1, 0));
    vector<vector<coord>> backtrack(n + 1, vector<coord>(m + 1, {-1, -1}));

    // Similarity score
    auto sim = [](char a, char b) {
        return a == b ? 2 : -2;
    };

    // Gap panelty score
    auto gap = [](int i) {
        // Use affine panelty here
//        return i + 10;
        return i + 1;
    };

    // dp process
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int match = score[i - 1][j - 1] + sim(reference[left_bound + i - 1], read[j - 1]);
            coord match_pos = {i - 1, j - 1};
            int del_score = -1;
            coord del_pos;
            int ins = -1;
            coord ins_pos;

            // Find the best score for deletion
            for (int k = 1; k <= min(i, max_gap); k++) {
                int tmp = score[i - k][j] - gap(k);
                if (tmp > del_score) {
                    del_score = tmp;
                    del_pos = {i - k, j};
                }
            }

            // Find the best score for insertion
            for (int k = 1; k <= min(j, max_gap); k++) {
                int tmp = score[i][j - k] - gap(k);
                if (tmp > ins) {
                    ins = tmp;
                    ins_pos = {i, j - k};
                }
            }

            // Find the best score
            int best_score = max(match, max(del_score, ins));
            if (best_score < 0) {
                best_score = 0;
            }

            // caveat: didn't consider multiple paths with same score
            if (best_score == match) {
                score[i][j] = match;
                backtrack[i][j] = match_pos;
            } else if (best_score == del_score) {
                score[i][j] = del_score;
                backtrack[i][j] = del_pos;
            } else {
                score[i][j] = ins;
                backtrack[i][j] = ins_pos;
            }
        }
    }

    // Find the best scores
    int best_score = -1;
    vector<coord> best_pos;
    for (int i = 0; i <= n; i++) {
        if (score[i][m] > best_score) {
            best_score = score[i][m];
            best_pos.clear();
            best_pos.emplace_back(i, m);
        } else if (score[i][m] == best_score) {
            best_pos.emplace_back(i, m);
        }
    }
    for (int i = 0; i <= n; i++) {
        if (score[i][m] != best_score && score[i][m] >= best_score - LOCAL_MATCH_SCORE_TOLERANCE) {
            best_pos.emplace_back(i, m);
        }
    }

    // Backtrack to find the best alignment
    vector<Mapping *> best_alignments;
    for (auto p: best_pos) {
        auto *map = new Mapping();
        coord last;
        while (p.first > 0 && p.second > 0) {
            if (verbose) cout << "(" << p.first << "," << p.second << ") ";
            last = p;

            coord prev = backtrack[p.first][p.second];
            // FIXME: wrong deletion and insertion positions
            if (prev.second == p.second) {
                // Deletion
                map->add_indel({IndelType::DELETION, left_bound + prev.first, abs(p.first - prev.first)}, logger);
                if (verbose) cout << "D" << abs(p.first - prev.first) << " -> ";
            } else if (prev.first == p.first) {
                // Insertion
                map->add_indel({IndelType::INSERTION, left_bound + p.first, abs(p.second - prev.second)}, logger);
                if (verbose) cout << "I" << abs(p.second - prev.second) << " -> ";
            } else if (reference[left_bound + p.first - 1] != read[p.second - 1]) {
                // Mismatch
                map->mismatches.insert({left_bound + p.first - 1, read[p.second - 1]});
                if (verbose) cout << "M" << " -> ";
            } else {
                // Match
                if (verbose) cout << "C" << " -> ";
            }

            p = prev;
        }
        if (verbose) cout << endl;

        // Read's front before reference's front
        if (last.second - 1 > 0) {
            map->add_indel({IndelType::INSERTION, left_bound, last.second - 1}, logger);
            if (verbose) cout << "I" << last.second - 1 << endl;
        }

        map->position = left_bound + last.first - 1;
        best_alignments.emplace_back(map);
        if (verbose) cout << map->to_string() << endl << endl;
    }

    sort(best_alignments.begin(), best_alignments.end(), [](const Mapping *a, const Mapping *b) {
        assert(a != nullptr && b != nullptr);
        int a_err = a->num_errors();
        int b_err = b->num_errors();
        if (a_err == b_err) {
            return a->indels.size() < b->indels.size();
        } else {
            return a_err < b_err;
        }
    });

    if (best_alignments.empty()) {
        logger->log_msg(LogLevel::DEBUG, "local_match", "No alignment found for read " + read);
        return nullptr;
    }

    return best_alignments[0];
}

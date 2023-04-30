#ifndef READMAPPING_ALGORITHMS_H
#define READMAPPING_ALGORITHMS_H

#include <string>
#include <vector>
#include <future>
#include <functional>

#include "types.h"

int binary_search(const std::vector<int>& v, int val);

/**
 * Given a vector of vectors of positions, returns all possible combinations of positions
 * such that the positions are strictly increasing in each combination.
 * @param sequences - The vector of vectors of positions
 * @param result    - The result vector of vectors of positions
 * @param pos       - The current position in the sequences vector
 */
void get_increasing_sequences(const std::vector<std::vector<int>> &sequences,
                              std::vector<std::vector<int>> &result,
                              int pos = 0);

/**
 * Match a read to a reference sequence using the Smith-Waterman algorithm.
 * Reference: https://www.wikiwand.com/en/Smith%E2%80%93Waterman_algorithm
 * @param reference   - The reference sequence
 * @param read        - The read sequence
 * @param left_bound  - The left bound of the search
 * @param right_bound - The right bound of the search
 * @param logger      - The logger object
 * @param max_gap     - The maximum gap allowed between two matches
 * @param verbose     - Whether to print debug information
 * @return            - The optimal mapping of the read to the reference
 */
Mapping *local_match(const std::string &reference, const std::string &read, int left_bound, int right_bound,
                     Logger *logger,
                     int max_gap = MAX_GAP, bool verbose = false);

#endif //READMAPPING_ALGORITHMS_H

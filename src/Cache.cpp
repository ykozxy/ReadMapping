#include "Cache.h"

#include <algorithm>
#include <functional>

#include "logging.h"

using namespace std;

#define L_TYPE true
#define S_TYPE false


/**
 * Build LS type array of input string s
 * @return  - type array: true if a position is L type, false if S type.
 */
VB Cache::build_typemap() {
    VB res(m_string.size() + 1);

    // Last pos is S type
    res[m_string.size()] = S_TYPE;

    if (m_string.empty()) return res;

    // Second to last pos is L type
    res[m_string.size() - 1] = L_TYPE;

    for (int i = (int) m_string.size() - 2; i >= 0; i--) {
        if (m_string[i] > m_string[i + 1]) {
            res[i] = L_TYPE;
        } else if (m_string[i] == m_string[i + 1] && res[i + 1]) {
            res[i] = L_TYPE;
        } else {
            res[i] = S_TYPE;
        }
    }

    return res;
}


/**
 * Build bucket sizes for input string s
 * @param alphabet_size - size of alphabet
 * @return              - vector of bucket sizes, in order A, C, G, T
 */
VI Cache::get_bucket_sizes(int alphabet_size) {
    VI res(alphabet_size, 0);
    for (auto c: m_string) {
        res[(int) c]++;
    }
    return res;
}


/**
 * Given a vector of bucket sizes, return a vector of bucket heads
 * @param bucket_sizes - vector of bucket sizes
 * @return             - vector of bucket heads
 */
VI Cache::get_bucket_heads(const VI &bucket_sizes) {
    int offset = 1;
    VI res;
    res.reserve(bucket_sizes.size());
    for (auto s: bucket_sizes) {
        res.push_back(offset);
        offset += s;
    }
    return res;
}


/**
 * Given a vector of bucket sizes, return a vector of bucket tails
 * @param bucket_sizes - vector of bucket sizes
 * @return             - vector of bucket tails
 */
VI Cache::get_bucket_tails(const VI &bucket_sizes) {
    int offset = 1;
    VI res;
    res.reserve(bucket_sizes.size());
    for (auto s: bucket_sizes) {
        offset += s;
        res.push_back(offset - 1);
    }
    return res;
}


/**
 * Check if a position is LMS
 * @param type_map - type array
 * @param i        - position
 * @return         - true if LMS, false otherwise
 */
bool Cache::is_LMS(const VB &type_map, size_t i) {
    if (i == 0) return false;
    if (type_map[i] == S_TYPE && type_map[i - 1] == L_TYPE) return true;
    return false;
}


/**
 * Check if two LMS substrings are equal
 * @param type_map - L/S type array
 * @param i        - start position of first LMS substring
 * @param j        - start position of second LMS substring
 * @return         - true if equal, false otherwise
 */
bool Cache::compare_LMS_substrings(const VB &type_map, size_t i, size_t j) {
    if (i == m_string.size() || j == m_string.size()) return false;

    int idx = 0;
    while (true) {
        bool i_LMS = is_LMS(type_map, i + idx);
        bool j_LMS = is_LMS(type_map, j + idx);

        if (idx > 0 && i_LMS && j_LMS) {
            return true;
        } else if (i_LMS != j_LMS || (m_string[i + idx] != m_string[j + idx])) {
            return false;
        }

        idx++;
    }
}


/**
 * Guess the order of LMS substrings
 * @param bucket_sizes - vector of bucket sizes
 * @param type_map     - type array
 * @return             - vector of LMS indices
 */
VI Cache::guess_LMS_sort(const VI &bucket_sizes, const VB &type_map) {
    VI res(m_string.size() + 1, -1);
    auto bucket_tails = get_bucket_tails(bucket_sizes);

    // Do bucket sort on LMS indices
    for (size_t i = 0; i < m_string.size(); i++) {
        if (is_LMS(type_map, i)) {
            res[bucket_tails[m_string[i]]--] = (int) i;
        }
    }

    res[0] = (int) m_string.size();
    return res;
}


/**
 * Induce sort L-type suffixes.
 * @param arr          - suffix array to be sorted
 * @param bucket_sizes - vector of bucket sizes
 * @param type_map     - L/S type array
 */
void Cache::induce_sort_L(VI &arr, const VI &bucket_sizes, const VB &type_map) {
    auto bucket_heads = get_bucket_heads(bucket_sizes);
    for (int i: arr) {
        if (i == -1) continue;

        int j = i - 1;
        if (j >= 0 && type_map[j] == L_TYPE) {
            // Only consider L-type suffixes
            arr[bucket_heads[m_string[j]]++] = j;
        }
    }
}


/**
 * Induce sort S-type suffixes.
 * @param arr          - suffix array to be sorted
 * @param bucket_sizes - vector of bucket sizes
 * @param type_map     - L/S type array
 */
void Cache::induce_sort_S(VI &arr, const VI &bucket_sizes, const VB &type_map) {
    auto bucket_tails = get_bucket_tails(bucket_sizes);
    for (int i = (int) arr.size() - 1; i >= 0; i--) {
//        if (arr[i] == -1) continue;

        int j = arr[i] - 1;
        if (j >= 0 && type_map[j] == S_TYPE) {
            // Only consider S-type suffixes
            arr[bucket_tails[m_string[j]]--] = j;
        }
    }
}


/**
 * Construct summarize string of LMS substrings
 * @param arr            - suffix array
 * @param type_map       - L/S type array
 * @param summary_str    - output summarize string
 * @param summary_offset - output summarize offset
 * @return               - alphabet size of summarize string
 */
int Cache::summerize_suffix_array(const VI &arr,
                                  const VB &type_map,
                                  VI &summary_str,
                                  VI &summary_offset) {
    VI lms_names(m_string.size() + 1, -1);

    int cur_name = 0;
    lms_names[arr[0]] = cur_name;
    int last_lms_offset = arr[0];

    for (int i = 1; i < (int) arr.size(); i++) {
        int offset = arr[i];
        if (is_LMS(type_map, offset)) {
            if (!compare_LMS_substrings(type_map, last_lms_offset, offset)) {
                cur_name++;
            }
            last_lms_offset = offset;
            lms_names[offset] = cur_name;
        }
    }

    // Construct summary string
    summary_str.clear();
    summary_offset.clear();
    for (int i = 0; i < (int) lms_names.size(); i++) {
        if (lms_names[i] != -1) {
            summary_str.push_back(lms_names[i]);
            summary_offset.push_back(i);
        }
    }

    return cur_name + 1;
}


VI Cache::construct_summery_suffix_array(const VI &summary_str, int alphabet_size) {
    VI summary_suffix_array;
    if (alphabet_size == summary_str.size()) {
        // Summary string is unique
        summary_suffix_array = VI(summary_str.size() + 1, -1);
        summary_suffix_array[0] = (int) summary_str.size();
        for (int i = 0; i < (int) summary_str.size(); i++) {
            summary_suffix_array[summary_str[i] + 1] = i;
        }
    } else {
        // Summary string is not unique
        Cache cache(summary_str);
        cache.build_suffix_array(alphabet_size);
        summary_suffix_array = cache.get_suffix_array();
    }

    return summary_suffix_array;
}

VI Cache::accurate_LMS_sort(const VI &bucket_sizes, const VI &summary_suffix_array, const VI &summary_offset) {
    VI res(m_string.size() + 1, -1);

    auto bucket_tails = get_bucket_tails(bucket_sizes);
    for (int i = (int) summary_suffix_array.size() - 1; i > 1; i--) {
        int offset = summary_offset[summary_suffix_array[i]];

        res[bucket_tails[m_string[offset]]--] = offset;
    }

    res[0] = (int) m_string.size();
    return res;
}


void Cache::build_suffix_array(int alphabet_size) {
    auto type_map = build_typemap();
    auto bucket_sizes = get_bucket_sizes(alphabet_size);

    auto guess = guess_LMS_sort(bucket_sizes, type_map);
    induce_sort_L(guess, bucket_sizes, type_map);
    induce_sort_S(guess, bucket_sizes, type_map);

    VI summary_str, summary_offset;
    int summary_alphabet_size = summerize_suffix_array(guess, type_map, summary_str, summary_offset);

    auto summary_suffix_array = construct_summery_suffix_array(summary_str, summary_alphabet_size);

    auto res = accurate_LMS_sort(bucket_sizes, summary_suffix_array, summary_offset);
    induce_sort_L(res, bucket_sizes, type_map);
    induce_sort_S(res, bucket_sizes, type_map);

    suffix_array = res;

    if (m_string.back() != -1) {
        // Push an extra -1 to the end of the string
        m_string.push_back(-1);
    }
}

Cache::Cache(const string &reference) {
    for (char c: reference) {
        m_string.push_back(c);
    }
}

Cache::Cache(const std::vector<int> &reference) {
    m_string = reference;
}

Cache::~Cache() = default;

vector<int> Cache::get_suffix_array() const {
    return suffix_array;
}

int Cache::binary_search(const string &pattern, bool left_bound) {
    int start = 0;
    int end = (int) suffix_array.size() - 1;
    int mid;
    while (start <= end) {
        mid = start + (end - start) / 2;
        int cmp = compare_string(pattern, suffix_array[mid], pattern.size());
        if (cmp == 0) {
            if (left_bound) {
                if (mid == 0 || compare_string(pattern, suffix_array[mid - 1], pattern.size()) != 0) {
                    return mid;
                } else {
                    end = mid - 1;
                }
            } else {
                if (mid == (int) suffix_array.size() - 1 ||
                    compare_string(pattern, suffix_array[mid + 1], pattern.size()) != 0) {
                    return mid;
                } else {
                    start = mid + 1;
                }
            }
        } else if (cmp < 0) {
            start = mid + 1;
        } else {
            end = mid - 1;
        }
    }
    return -1;
}

VI Cache::find(const string &pattern) {
    VI occurrences;

    if (pattern.empty()) {
        return occurrences;
    }

    // Find the leftmost occurrence of the pattern
    int left_bound = binary_search(pattern, true);
    if (left_bound == -1) {
        return occurrences;
    }

    // Find the rightmost occurrence of the pattern
    int right_bound = binary_search(pattern, false);
    if (right_bound == -1) {
        return occurrences;
    }

    // Iterate through the suffix array within the bounds and add the occurrences
    for (int i = left_bound; i <= right_bound; i++) {
        occurrences.push_back(suffix_array[i]);
    }

    return occurrences;
}


int Cache::compare_string(const std::string &pattern, int start, int size) {
    for (int i = 0; i < size; i++) {
        if ((int) pattern[i] < m_string[start + i]) {
            return 1;
        } else if ((int) pattern[i] > m_string[start + i]) {
            return -1;
        }
    }
    return 0;
}

#ifndef READMAPPING_CACHE_H
#define READMAPPING_CACHE_H

#include <string>
#include <vector>

#include "types.h"

using VI = std::vector<int>;
using VB = std::vector<bool>;

/* A suffix tree */
class Cache {
public:
    explicit Cache(const std::string &reference);

    explicit Cache(const std::vector<int> &reference);

    ~Cache();

    void build_suffix_array(int alphabet_size = 256);

    std::vector<int> find(const std::string &pattern);

    [[nodiscard]] std::vector<int> get_suffix_array() const;

private:
    std::vector<int> suffix_array;
    std::vector<int> m_string;
    std::string reference;

    int binary_search(const std::string &pattern, bool left_bound);

    int compare_string(const std::string &pattern, int start, int size);

    /* SAIS algorithm functions */
    VB build_typemap();

    VI get_bucket_sizes(int alphabet_size);

    static VI get_bucket_heads(const VI &bucket_sizes);

    static VI get_bucket_tails(const VI &bucket_sizes);

    static bool is_LMS(const VB &type_map, size_t i);

    bool compare_LMS_substrings(const VB &type_map, size_t i, size_t j);

    VI guess_LMS_sort(const VI &bucket_sizes, const VB &type_map);

    void induce_sort_L(VI &arr, const VI &bucket_sizes, const VB &type_map);

    void induce_sort_S(VI &arr, const VI &bucket_sizes, const VB &type_map);

    int summerize_suffix_array(const VI &arr,
                               const VB &type_map,
                               VI &summary_str,
                               VI &summary_offset);

    static VI construct_summery_suffix_array(const VI &summary_str, int alphabet_size);

    VI accurate_LMS_sort(const VI &bucket_sizes, const VI &summary_suffix_array, const VI &summary_offset);
};


#endif //READMAPPING_CACHE_H

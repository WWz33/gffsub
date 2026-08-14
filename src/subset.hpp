#ifndef SUBSET_HPP
#define SUBSET_HPP

#include "annotation.hpp"
#include "filter.hpp"
#include "region.hpp"
#include "selector_filter.hpp"

#include <cstddef>
#include <optional>
#include <string>
#include <unordered_set>
#include <vector>

namespace gffsub {

struct SubsetParams {
    std::optional<Region> region;
    std::string bed_file;
    std::string seqid_filter;       // ^-prefixed for exclude
    std::string source_filter;
    std::optional<std::optional<double>> score_filter;
    std::optional<char> strand_filter;
    std::optional<char> phase_filter;
    std::string feature;
    bool longest = false;
    size_t threads = 1;
    std::vector<GrepFilter> grep_filters;
    std::vector<ExprNode> include_exprs;
    std::vector<ExprNode> exclude_exprs;
    bool invert_grep = false;
};

// Apply region, bed, seqid, source, score, strand, phase, grep, expr, and
// feature filters in sequence. Modifies data in place.
void subset(GffData& data, IdIndex& idx, const SubsetParams& params);

}  // namespace gffsub

#endif  // SUBSET_HPP

#ifndef GFFSUB_SELECTOR_FILTER_HPP
#define GFFSUB_SELECTOR_FILTER_HPP

#include "annotation.hpp"
#include "expr_parser.hpp"

#include <optional>
#include <regex>
#include <string>
#include <string_view>
#include <vector>

namespace gffsub {

struct GrepFilter {
    std::string field;
    std::string pattern;
    bool use_regex = false;
    bool ignore_case = false;
    std::optional<std::regex> compiled;
};

void filter_by_grep(GffData& data, const std::vector<GrepFilter>& filters, bool invert);
void filter_by_expr(GffData& data, const std::vector<ExprNode>& filters, bool include);

}  // namespace gffsub

#endif  // GFFSUB_SELECTOR_FILTER_HPP

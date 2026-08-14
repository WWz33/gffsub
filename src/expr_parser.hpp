#ifndef GFFSUB_EXPR_PARSER_HPP
#define GFFSUB_EXPR_PARSER_HPP

#include <optional>
#include <regex>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace gffsub {

enum class ExprOp {
    Equal,
    NotEqual,
    Regex,
    NotRegex,
    Less,
    LessEqual,
    Greater,
    GreaterEqual
};

struct ExprFilter {
    std::string field;
    ExprOp op;
    std::string value;
    bool ignore_case = false;
    std::optional<std::regex> compiled;
};

struct ExprNode {
    enum class Kind {
        Predicate,
        Not,
        And,
        Or
    };

    Kind kind = Kind::Predicate;
    ExprFilter predicate;
    std::vector<ExprNode> children;
};

std::optional<std::pair<std::string, std::string>> parse_field_pattern(std::string_view value);
bool parse_expr_filters(std::string_view expr, std::vector<ExprNode>& out, std::string& error);
bool compile_filter_regexes(std::vector<struct GrepFilter>& grep_filters,
                            std::vector<ExprNode>& include_filters,
                            std::vector<ExprNode>& exclude_filters,
                            bool ignore_case,
                            std::string& error);

}  // namespace gffsub

#endif  // GFFSUB_EXPR_PARSER_HPP

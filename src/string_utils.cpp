#include "string_utils.hpp"

#include <cctype>

namespace gffsub {

std::string to_lower(std::string_view sv) {
    std::string s;
    s.reserve(sv.size());
    for (char ch : sv) {
        s.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(ch))));
    }
    return s;
}

std::vector<std::string> split_line(const std::string& line, char delimiter) {
    std::vector<std::string> cols;
    cols.reserve(delimiter == '\t' ? 9 : 4);
    size_t start = 0;
    while (true) {
        auto pos = line.find(delimiter, start);
        if (pos == std::string::npos) {
            cols.emplace_back(line.substr(start));
            break;
        }
        cols.emplace_back(line.substr(start, pos - start));
        start = pos + 1;
    }
    return cols;
}

std::unordered_set<std::string> parse_list(std::string_view sv, bool& exclude) {
    exclude = false;
    if (!sv.empty() && sv.front() == '^') {
        exclude = true;
        sv.remove_prefix(1);
    }
    std::unordered_set<std::string> result;
    size_t pos = 0;
    while (pos < sv.size()) {
        auto comma = sv.find(',', pos);
        if (comma == std::string_view::npos) comma = sv.size();
        if (comma > pos) result.emplace(sv.substr(pos, comma - pos));
        pos = comma + 1;
    }
    return result;
}

}  // namespace gffsub

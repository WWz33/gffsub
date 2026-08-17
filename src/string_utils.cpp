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

}  // namespace gffsub

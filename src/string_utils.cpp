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

// URL-decode a GFF3 attribute value per spec: %XX -> byte.
std::string url_decode(std::string_view input) {
    std::string out;
    out.reserve(input.size());
    for (size_t i = 0; i < input.size(); ++i) {
        if (input[i] == '%' && i + 2 < input.size()) {
            auto hex_val = [](char c) -> int {
                if (c >= '0' && c <= '9') return c - '0';
                if (c >= 'a' && c <= 'f') return c - 'a' + 10;
                if (c >= 'A' && c <= 'F') return c - 'A' + 10;
                return -1;
            };
            const int hi = hex_val(input[i + 1]);
            const int lo = hex_val(input[i + 2]);
            if (hi >= 0 && lo >= 0) {
                out.push_back(static_cast<char>(hi * 16 + lo));
                i += 2;
            } else {
                out.push_back('%');
            }
        } else {
            out.push_back(input[i]);
        }
    }
    return out;
}

}  // namespace gffsub

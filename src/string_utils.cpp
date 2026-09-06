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

std::vector<std::string_view> split_line(std::string_view line, char delimiter) {
    std::vector<std::string_view> cols;
    cols.reserve(delimiter == '\t' ? 9 : 4);
    size_t start = 0;
    while (true) {
        auto pos = line.find(delimiter, start);
        if (pos == std::string_view::npos) {
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

int natural_compare(std::string_view a, std::string_view b) {
    size_t i = 0, j = 0;
    while (i < a.size() && j < b.size()) {
        if (std::isdigit(static_cast<unsigned char>(a[i])) &&
            std::isdigit(static_cast<unsigned char>(b[j]))) {
            // Compare digit runs numerically. Strip leading zeros, then
            // longer run wins, then lexicographic on equal-length digits.
            size_t ia = i, jb = j;
            while (ia < a.size() && std::isdigit(static_cast<unsigned char>(a[ia]))) ++ia;
            while (jb < b.size() && std::isdigit(static_cast<unsigned char>(b[jb]))) ++jb;
            std::string_view da = a.substr(i, ia - i), db = b.substr(j, jb - j);
            const auto nza = da.find_first_not_of('0');
            const auto nzb = db.find_first_not_of('0');
            da.remove_prefix(nza == std::string_view::npos ? da.size() : nza);
            db.remove_prefix(nzb == std::string_view::npos ? db.size() : nzb);
            if (da.size() != db.size()) return da.size() < db.size() ? -1 : 1;
            const int c = da.compare(db);
            if (c != 0) return c;
            i = ia; j = jb;
            continue;
        }
        if (a[i] != b[j]) return a[i] < b[j] ? -1 : 1;
        ++i; ++j;
    }
    if (i < a.size()) return 1;
    if (j < b.size()) return -1;
    return 0;
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

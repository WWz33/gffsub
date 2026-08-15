#include "parser.hpp"

#include <string>
#include <string_view>

namespace gffsub {

namespace {

// URL-decode a GFF3 attribute value per spec: %XX → byte.
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

}  // namespace

std::unordered_map<std::string, std::vector<std::string>> parse_attributes(std::string_view attrs) {
    std::unordered_map<std::string, std::vector<std::string>> parsed;
    size_t pos = 0;
    while (pos < attrs.size()) {
        // Find the end of the current tag=value pair (next ';').
        size_t pair_end = attrs.find(';', pos);
        if (pair_end == std::string_view::npos) {
            pair_end = attrs.size();
        }
        const auto pair = attrs.substr(pos, pair_end - pos);

        // Skip empty fragments (e.g. from "key=val;;key2=val2" or leading/trailing ';').
        if (!pair.empty()) {
            const size_t eq = pair.find('=');
            if (eq != std::string_view::npos && eq != 0 && eq + 1 <= pair.size()) {
                const std::string key{pair.substr(0, eq)};
                std::string_view value = pair.substr(eq + 1);
                if (!key.empty()) {
                    // Split the RAW value on ',' before URL-decoding each part:
                    // ',' is the value-list separator and a literal comma inside
                    // a single value must be escaped as %2C (GFF3 spec).
                    size_t part_start = 0;
                    while (part_start <= value.size()) {
                        size_t part_end = value.find(',', part_start);
                        if (part_end == std::string_view::npos) {
                            part_end = value.size();
                        }
                        const std::string part = url_decode(value.substr(part_start, part_end - part_start));
                        if (!part.empty()) {
                            parsed[key].push_back(part);
                        }
                        if (part_end == value.size()) {
                            break;
                        }
                        part_start = part_end + 1;
                    }
                }
            }
        }

        pos = (pair_end < attrs.size()) ? pair_end + 1 : attrs.size();
    }
    return parsed;
}

}  // namespace gffsub

#include "parser.hpp"
#include "string_utils.hpp"

#include <string>
#include <string_view>

namespace gffsub {

namespace {

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

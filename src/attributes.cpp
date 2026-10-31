#include "parser.hpp"
#include "gtf_parser.hpp"
#include "string_utils.hpp"

#include <algorithm>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

namespace gffsub {

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
                // Trim spaces around the key: emitters write `ID=x; Parent=y`
                // (spec allows unescaped spaces in column 9).
                const std::string key{trim_view(pair.substr(0, eq))};
                std::string_view value = pair.substr(eq + 1);
                if (!key.empty()) {
                    // GFF3 reserves the comma list for five tags; for any other
                    // tag a comma is part of the value. Splitting everything
                    // would turn `inference=ab initio,similar to` into two
                    // values and duplicate the tag on GTF output.
                    static const std::unordered_set<std::string> kMultiValue = {
                        "Parent", "Alias", "Note", "Dbxref", "Ontology_term"
                    };
                    if (kMultiValue.count(key) == 0) {
                        const std::string single = url_decode(value);
                        if (!single.empty()) {
                            parsed[key].push_back(single);
                        }
                        pos = (pair_end < attrs.size()) ? pair_end + 1 : attrs.size();
                        continue;
                    }
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

namespace {

// Split a raw Parent value on commas, trimming each part and keeping order.
void append_parent_ids(std::vector<std::string>& out, std::string_view raw) {
    size_t start = 0;
    while (start <= raw.size()) {
        const size_t comma = raw.find(',', start);
        const size_t end = (comma == std::string_view::npos) ? raw.size() : comma;
        std::string part{trim_view(raw.substr(start, end - start))};
        if (!part.empty() &&
            std::find(out.begin(), out.end(), part) == out.end()) {
            out.push_back(std::move(part));
        }
        if (comma == std::string_view::npos) break;
        start = comma + 1;
    }
}

}  // namespace

std::vector<std::string> parent_ids(const GffRecord& rec) {
    std::vector<std::string> out;
    bool from_parent_attr = false;
    if (rec.src_fmt == InputFormat::GTF) {
        // GTF uses `key "value";`; parse_attributes cannot read it. A quoted
        // comma list is split here so the index, --drop-orphans, --longest,
        // and the GTF writer all see the same parents. Most GTF rows have no
        // Parent attribute, so skip the parse when the text is absent.
        if (rec.attr_raw.find("Parent") != std::string_view::npos) {
            for (const auto& [key, value] : parse_gtf_attributes(rec.attr_raw)) {
                if (key == "Parent") {
                    from_parent_attr = true;
                    append_parent_ids(out, value);
                }
            }
        }
    } else {
        const auto attrs = parse_attributes(rec.attr_raw);
        const auto it = attrs.find("Parent");
        if (it != attrs.end()) {
            from_parent_attr = true;
            // parse_attributes already split raw commas; a %2C comma is part
            // of one value, so do not split again.
            for (const auto& value : it->second) {
                if (std::find(out.begin(), out.end(), value) == out.end()) {
                    out.push_back(value);
                }
            }
        }
    }
    const bool parents_from_gtf_attr =
        from_parent_attr && rec.src_fmt == InputFormat::GTF && !out.empty();
    if (rec.parent_id && !parents_from_gtf_attr &&
        std::find(out.begin(), out.end(), *rec.parent_id) == out.end()) {
        // GTF without a Parent attribute: the parser synthesizes parent_id
        // from transcript_id/gene_id. GFF3: first Parent value.
        append_parent_ids(out, *rec.parent_id);
    }
    return out;
}

}  // namespace gffsub

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

namespace {

// GFF3 reserves the comma list for these five tags (gff3.md:94). For any
// other tag a comma is part of the value: splitting `gene_id=a,b` would
// turn one value into two.
bool is_gff3_multi_value(std::string_view key) {
    return key == "Parent" || key == "Alias" || key == "Note" ||
           key == "Dbxref" || key == "Ontology_term";
}

// GTF identity keys. GTF has no URL escaping, but a value may carry the same
// %2C escape a GFF3 value does, so these are decoded before the identity is
// stored. Descriptive attributes pass through verbatim.
bool is_identity_key(std::string_view key) {
    return key == "ID" || key == "Parent" ||
           key == "gene_id" || key == "transcript_id";
}

void append_split_parts(std::vector<std::pair<std::string, std::string>>& out,
                        const std::string& key, std::string_view raw,
                        bool trim_parts) {
    size_t start = 0;
    while (start <= raw.size()) {
        const size_t comma = raw.find(',', start);
        const size_t end = (comma == std::string_view::npos) ? raw.size() : comma;
        std::string_view part = raw.substr(start, end - start);
        if (trim_parts) part = trim_view(part);
        std::string decoded = url_decode(part);
        if (!decoded.empty()) {
            out.emplace_back(key, std::move(decoded));
        }
        if (comma == std::string_view::npos) break;
        start = comma + 1;
    }
}

}  // namespace

std::vector<std::pair<std::string, std::string>> decode_attributes(
    std::string_view column9, InputFormat format) {
    std::vector<std::pair<std::string, std::string>> out;
    if (format == InputFormat::GTF) {
        for (auto& [key, value] : parse_gtf_attributes(column9)) {
            if (key == "Parent") {
                append_split_parts(out, key, value, /*trim_parts=*/true);
            } else if (is_identity_key(key)) {
                std::string decoded = url_decode(value);
                if (!decoded.empty()) {
                    out.emplace_back(std::move(key), std::move(decoded));
                }
            } else {
                out.emplace_back(std::move(key), std::move(value));
            }
        }
        return out;
    }

    // GFF3 (BED records carry no column 9).
    size_t pos = 0;
    while (pos < column9.size()) {
        size_t pair_end = column9.find(';', pos);
        if (pair_end == std::string_view::npos) {
            pair_end = column9.size();
        }
        const auto fragment = column9.substr(pos, pair_end - pos);
        pos = (pair_end < column9.size()) ? pair_end + 1 : column9.size();
        if (fragment.empty()) continue;

        const size_t eq = fragment.find('=');
        if (eq == std::string_view::npos || eq == 0) continue;
        // Trim around the key: emitters write `ID=x; Parent=y`.
        std::string key{trim_view(fragment.substr(0, eq))};
        if (key.empty()) continue;
        const auto value = fragment.substr(eq + 1);

        if (is_gff3_multi_value(key)) {
            // Split the RAW value on ',' before decoding: a literal comma
            // inside one value must be escaped as %2C (GFF3 spec). Parent
            // parts are identity, so a space after the comma is trimmed.
            append_split_parts(out, key, value, key == "Parent");
        } else {
            std::string decoded = url_decode(value);
            if (!decoded.empty()) {
                out.emplace_back(std::move(key), std::move(decoded));
            }
        }
    }
    return out;
}

std::unordered_map<std::string, std::vector<std::string>> parse_attributes(
    std::string_view attrs) {
    std::unordered_map<std::string, std::vector<std::string>> parsed;
    for (auto& [key, value] : decode_attributes(attrs, InputFormat::GFF3)) {
        parsed[key].push_back(std::move(value));
    }
    return parsed;
}

std::vector<std::string> attr_values(const GffRecord& rec, std::string_view key) {
    std::vector<std::string> out;
    for (auto& [k, value] : decode_attributes(rec.attr_raw, rec.src_fmt)) {
        if (k == key) {
            out.push_back(std::move(value));
        }
    }
    // Synthesized fields: BED records have no column 9; the GTF parser fills
    // ID/parent_id/gene_id/transcript_id from the same column. Append the
    // cached value when the decoded list did not already carry it.
    const std::optional<std::string>* cached = nullptr;
    if (key == "ID") cached = &rec.id;
    else if (key == "Parent") cached = &rec.parent_id;
    else if (key == "gene_id") cached = &rec.gene_id;
    else if (key == "transcript_id") cached = &rec.transcript_id;
    if (cached && *cached &&
        std::find(out.begin(), out.end(), **cached) == out.end()) {
        out.push_back(**cached);
    }
    return out;
}

std::vector<std::string> parent_ids(const GffRecord& rec) {
    std::vector<std::string> out;
    for (auto& id : attr_values(rec, "Parent")) {
        if (std::find(out.begin(), out.end(), id) == out.end()) {
            out.push_back(std::move(id));
        }
    }
    if (rec.src_fmt == InputFormat::GTF &&
        rec.feat_class == FeatureClass::Transcript && rec.id) {
        // A transcript is never its own parent. A non-standard Parent
        // attribute that lists sibling transcripts would otherwise link the
        // transcript into a cycle.
        out.erase(std::remove(out.begin(), out.end(), *rec.id), out.end());
    }
    return out;
}

}  // namespace gffsub

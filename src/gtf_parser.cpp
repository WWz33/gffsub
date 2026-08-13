#include "gtf_parser.hpp"

#include <optional>
#include <string>

namespace gffsub {
namespace {

std::optional<std::string> extract_quoted_value(const std::string& attrs, const std::string& key) {
    // GTF attributes are "; "-delimited: key "value";
    // Match the key as a whole token, not as a substring of a longer name.
    size_t pos = 0;
    while (pos < attrs.size()) {
        // Find the next occurrence of the key.
        size_t hit = attrs.find(key, pos);
        if (hit == std::string::npos) {
            return std::nullopt;
        }
        // Left boundary: must be at start, or preceded by ';' (optionally spaces).
        bool left_ok = (hit == 0);
        if (!left_ok) {
            size_t p = hit;
            while (p > 0 && (attrs[p - 1] == ' ' || attrs[p - 1] == '\t')) {
                --p;
            }
            left_ok = (p > 0 && attrs[p - 1] == ';');
        }
        if (!left_ok) {
            pos = hit + 1;
            continue;
        }
        // Right boundary: after the key, skip spaces, expect '"'.
        size_t q = hit + key.size();
        while (q < attrs.size() && (attrs[q] == ' ' || attrs[q] == '\t')) {
            ++q;
        }
        if (q >= attrs.size() || attrs[q] != '"') {
            pos = hit + 1;
            continue;
        }
        size_t q1 = q;
        size_t q2 = attrs.find('"', q1 + 1);
        if (q2 == std::string::npos) {
            return std::nullopt;
        }
        return attrs.substr(q1 + 1, q2 - q1 - 1);
    }
    return std::nullopt;
}

}  // namespace

void apply_gtf_attributes(GffRecord& rec) {
    if (!rec.gene_id) {
        rec.gene_id = extract_quoted_value(rec.attr_raw, "gene_id");
    }
    if (!rec.transcript_id) {
        rec.transcript_id = extract_quoted_value(rec.attr_raw, "transcript_id");
    }
}

}  // namespace gffsub

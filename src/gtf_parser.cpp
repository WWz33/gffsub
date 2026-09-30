#include "gtf_parser.hpp"

#include <optional>
#include <string>
#include <vector>
#include "string_utils.hpp"

namespace gffsub {

std::optional<std::string> extract_quoted_value(std::string_view attrs, std::string_view key) {
    // Whole-token match over the parsed pairs: avoids matching a key that
    // only appears as a substring of a longer one (`ref_gene_id`) or inside
    // another attribute's value, and accepts bare (unquoted) values.
    for (const auto& [k, v] : parse_gtf_attributes(attrs)) {
        if (k == key) {
            return v;
        }
    }
    return std::nullopt;
}

void apply_gtf_attributes(GffRecord& rec) {
    // One pass over column 9; both keys come from the same parse.
    if (!rec.gene_id || !rec.transcript_id) {
        for (const auto& [key, value] : parse_gtf_attributes(rec.attr_raw)) {
            if (!rec.gene_id && key == "gene_id") {
                rec.gene_id = value;
            } else if (!rec.transcript_id && key == "transcript_id") {
                rec.transcript_id = value;
            }
            if (rec.gene_id && rec.transcript_id) break;
        }
    }
}

namespace {

std::string gtf_unescape(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (size_t j = 0; j < s.size(); ++j) {
        if (s[j] == '\\' && j + 1 < s.size() && s[j + 1] == '"') {
            out.push_back('"');
            ++j;
        } else if (s[j] == '\\' && j + 1 < s.size() && s[j + 1] == '\\') {
            out.push_back('\\');
            ++j;
        } else {
            out.push_back(s[j]);
        }
    }
    return out;
}

// Emit one `;`-separated fragment as a key/value pair. Quoted values keep
// `;` verbatim (handled by the caller's quote tracking); bare values
// (non-standard GTF emitted by some tools) are accepted as-is.
void emit_gtf_fragment(std::string_view frag,
                       std::vector<std::pair<std::string, std::string>>& out) {
    frag = trim_view(frag);
    if (frag.empty()) return;

    const auto q1 = frag.find('"');
    if (q1 == std::string_view::npos) {
        // Bare `key value` form.
        const auto sp = frag.find_first_of(" \t");
        if (sp == std::string_view::npos) return;  // key without a value
        const auto key = trim_view(frag.substr(0, sp));
        const auto value = trim_view(frag.substr(sp));
        if (key.empty() || value.empty()) return;
        out.emplace_back(std::string{key}, std::string{value});
        return;
    }

    const auto key = trim_view(frag.substr(0, q1));
    if (key.empty()) return;
    // Find the closing quote, honoring backslash escapes.
    size_t q2 = q1 + 1;
    while (q2 < frag.size()) {
        if (frag[q2] == '\\' && q2 + 1 < frag.size()) {
            q2 += 2;
            continue;
        }
        if (frag[q2] == '"') break;
        ++q2;
    }
    if (q2 >= frag.size()) return;  // unclosed quote: drop the fragment
    out.emplace_back(std::string{key}, gtf_unescape(std::string{frag.substr(q1 + 1, q2 - q1 - 1)}));
}

}  // namespace

std::vector<std::pair<std::string, std::string>> parse_gtf_attributes(std::string_view attrs) {
    std::vector<std::pair<std::string, std::string>> result;
    // Split on ';' outside double quotes so a semicolon inside a quoted value
    // (`note "a;b"`) is preserved; AGAT's parser does the same.
    size_t frag_start = 0;
    bool in_quote = false;
    for (size_t i = 0; i < attrs.size(); ++i) {
        const char c = attrs[i];
        if (c == '\\' && in_quote && i + 1 < attrs.size()) {
            ++i;
            continue;
        }
        if (c == '"') {
            in_quote = !in_quote;
        } else if (c == ';' && !in_quote) {
            emit_gtf_fragment(attrs.substr(frag_start, i - frag_start), result);
            frag_start = i + 1;
        }
    }
    emit_gtf_fragment(attrs.substr(frag_start), result);
    return result;
}

std::string gtf_attrs_to_gff3(const GffRecord& rec) {
    static const auto url_escape = [](const std::string& s) {
        std::string out;
        out.reserve(s.size());
        for (const unsigned char ch : s) {
            if (ch == ',') { out += "%2C"; }
            else if (ch == ';') { out += "%3B"; }
            else if (ch == '=') { out += "%3D"; }
            else if (ch == '&') { out += "%26"; }
            else if (ch == '%') { out += "%25"; }
            else if (ch == '"') { out += "%22"; }
            else if (ch == '\t') { out += "%09"; }
            else if (ch == '\n') { out += "%0A"; }
            else if (ch == '\r') { out += "%0D"; }
            else if (ch < 0x20 || ch == 0x7F) {
                static const char hex[] = "0123456789ABCDEF";
                out += '%';
                out += hex[ch >> 4];
                out += hex[ch & 0x0F];
            } else {
                out.push_back(static_cast<char>(ch));
            }
        }
        return out;
    };

    std::vector<std::string> parts;
    if (rec.id) {
        parts.push_back("ID=" + url_escape(*rec.id));
    }
    if (rec.parent_id) {
        parts.push_back("Parent=" + url_escape(*rec.parent_id));
    }

    // Convert the remaining `key "value";` pairs to key=value.
    // Skip gene_id and transcript_id — they were synthesized into ID=/Parent=
    // above and must not be re-emitted (AGAT replaces them with Parent=).
    for (const auto& [key, value] : parse_gtf_attributes(rec.attr_raw)) {
        if (key == "gene_id" || key == "transcript_id") {
            continue;
        }
        parts.push_back(url_escape(key) + "=" + url_escape(value));
    }

    if (parts.empty()) {
        return ".";
    }
    std::string out;
    for (size_t i = 0; i < parts.size(); ++i) {
        if (i > 0) {
            out += ';';
        }
        out += parts[i];
    }
    return out;
}

}  // namespace gffsub

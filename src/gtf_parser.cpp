#include "gtf_parser.hpp"
#include "parser.hpp"

#include <optional>
#include <string>
#include <vector>
#include "string_utils.hpp"

namespace gffsub {

namespace {

std::string gtf_unescape(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (size_t j = 0; j < s.size(); ++j) {
        if (s[j] != '\\' || j + 1 >= s.size()) {
            out.push_back(s[j]);
            continue;
        }
        // GXF2 defines C-style escapes for non-printing characters, and the
        // GTF writer emits \t, \n and \r, so decoding them keeps a round trip
        // value-for-value.
        switch (s[j + 1]) {
            case '"':  out.push_back('"');  ++j; break;
            case '\\': out.push_back('\\'); ++j; break;
            case 't':  out.push_back('\t'); ++j; break;
            case 'n':  out.push_back('\n'); ++j; break;
            case 'r':  out.push_back('\r'); ++j; break;
            default:   out.push_back('\\'); break;  // unknown escape kept as-is
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
        // Bare form: `key value` (non-standard GTF from some tools) or
        // GFF3-style `key=value`. Cut at whichever separator comes first so
        // `ID=e1` does not fall through the "key without a value" drop.
        const auto sp = frag.find_first_of(" \t");
        const auto eq = frag.find('=');
        size_t cut = std::string_view::npos;
        if (sp != std::string_view::npos) cut = sp;
        if (eq != std::string_view::npos && (cut == std::string_view::npos || eq < cut)) cut = eq;
        if (cut == std::string_view::npos) return;  // key without a value
        const auto key = trim_view(frag.substr(0, cut));
        const auto value = trim_view(frag.substr(cut + (cut == eq ? 1 : 0)));
        if (key.empty() || value.empty()) return;
        out.emplace_back(std::string{key}, std::string{value});
        return;
    }

    std::string key{trim_view(frag.substr(0, q1))};
    // Non-standard `key="value"` (GFF3 separator before the quote): the `=`
    // belongs to the separator, not to the tag name.
    while (!key.empty() && (key.back() == '=' || key.back() == ' ' || key.back() == '\t')) {
        key.pop_back();
    }
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
    // GTF2.2: a comment starts at `#` and continues to the end of the line,
    // and may follow the attributes. `#` inside a quoted value is data.
    {
        bool quoted = false;
        for (size_t i = 0; i < attrs.size(); ++i) {
            const char c = attrs[i];
            if (c == '\\' && quoted && i + 1 < attrs.size()) {
                ++i;
                continue;
            }
            if (c == '"') {
                quoted = !quoted;
            } else if (c == '#' && !quoted) {
                attrs = attrs.substr(0, i);
                break;
            }
        }
    }
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
    {
        // Build Parent from the full parent list, not just the scalar first
        // value: a GTF `Parent "t1,t2"` row must keep both links in GFF3.
        const auto parents = parent_ids(rec);
        if (!parents.empty()) {
            std::string parent_part = "Parent=";
            for (size_t i = 0; i < parents.size(); ++i) {
                if (i > 0) parent_part += ',';
                parent_part += url_escape(parents[i]);
            }
            parts.push_back(std::move(parent_part));
        }
    }

    // Convert the remaining `key "value";` pairs to key=value.
    // transcript_id becomes Parent= and is never re-emitted. gene_id becomes
    // ID= on a gene row and Parent= on a transcript row, but on an exon/CDS
    // row it is the only gene reference, so it stays as a gene_id attribute.
    // A bare ID/Parent fragment produced the synthesized pair above and is
    // skipped.
    for (const auto& [key, value] : decode_attributes(rec.attr_raw, InputFormat::GTF)) {
        if (key == "transcript_id" || key == "ID" || key == "Parent") {
            continue;
        }
        if (key == "gene_id" &&
            (rec.feat_class == FeatureClass::Gene ||
             rec.feat_class == FeatureClass::Transcript)) {
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

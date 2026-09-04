#include "gtf_parser.hpp"

#include <optional>
#include <string>
#include <vector>

namespace gffsub {

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
        // Find closing quote, handling escaped quotes (\" -> literal quote in value)
        size_t q2 = q1 + 1;
        while (q2 < attrs.size()) {
            if (attrs[q2] == '\\' && q2 + 1 < attrs.size()) {
                q2 += 2;  // skip escaped char
                continue;
            }
            if (attrs[q2] == '"') break;
            ++q2;
        }
        if (q2 >= attrs.size()) {
            return std::nullopt;
        }
        // Unescape: remove backslash before quotes and backslashes
        std::string value = attrs.substr(q1 + 1, q2 - q1 - 1);
        std::string unescaped;
        unescaped.reserve(value.size());
        for (size_t j = 0; j < value.size(); ++j) {
            if (value[j] == '\\' && j + 1 < value.size() && value[j + 1] == '"') {
                unescaped.push_back('"');
                ++j;
            } else if (value[j] == '\\' && j + 1 < value.size() && value[j + 1] == '\\') {
                unescaped.push_back('\\');
                ++j;
            } else {
                unescaped.push_back(value[j]);
            }
        }
        return unescaped;
    }
    return std::nullopt;
}

void apply_gtf_attributes(GffRecord& rec) {
    if (!rec.gene_id) {
        rec.gene_id = extract_quoted_value(rec.attr_raw, "gene_id");
    }
    if (!rec.transcript_id) {
        rec.transcript_id = extract_quoted_value(rec.attr_raw, "transcript_id");
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

}  // namespace

std::vector<std::pair<std::string, std::string>> parse_gtf_attributes(const std::string& attrs) {
    std::vector<std::pair<std::string, std::string>> result;
    size_t pos = 0;
    while (pos < attrs.size()) {
        size_t end = attrs.find(';', pos);
        if (end == std::string::npos) {
            end = attrs.size();
        }
        std::string frag = attrs.substr(pos, end - pos);
        pos = (end < attrs.size()) ? end + 1 : attrs.size();

        const auto first = frag.find_first_not_of(" \t");
        if (first == std::string::npos) {
            continue;
        }
        const auto last = frag.find_last_not_of(" \t");
        frag = frag.substr(first, last - first + 1);

        const auto q1 = frag.find('"');
        if (q1 == std::string::npos || q1 == 0) {
            continue;
        }
        const auto q2 = frag.rfind('"');
        if (q2 == q1) {
            continue;
        }
        std::string key = frag.substr(0, q1);
        const auto klast = key.find_last_not_of(" \t");
        if (klast == std::string::npos) {
            continue;
        }
        key = key.substr(0, klast + 1);
        const std::string raw_value = frag.substr(q1 + 1, q2 - q1 - 1);
        if (key.empty() || raw_value.empty()) {
            continue;
        }
        result.emplace_back(std::move(key), gtf_unescape(raw_value));
    }
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

#include "annotation.hpp"
#include "feature_types.hpp"
#include "parser.hpp"
#include "record.hpp"
#include "region.hpp"
#include "gtf_parser.hpp"
#include "input_source.hpp"
#include "string_utils.hpp"
#include <iostream>
#include <stdexcept>

namespace gffsub {

namespace {

// A filtering tool should not reject an annotation because a handful of lines
// are malformed, but it should say so. Counted per parse and reported once.
struct ParseIssues {
    size_t short_lines = 0;       // fewer than 9 columns (or 3 for BED)
    size_t extra_columns = 0;     // more than 9 columns (or 12 for BED)
    size_t bad_coordinates = 0;   // non-integer, non-positive, or start > end
    size_t cds_without_phase = 0; // CDS whose phase is not 0, 1 or 2
    size_t bad_strand = 0;        // strand that is not +, -, ., or ?
    size_t bad_phase = 0;         // phase that is not 0, 1, 2, or .
};
thread_local ParseIssues g_issues;

void reset_issues() { g_issues = ParseIssues{}; }

void report_issues() {
    if (g_issues.short_lines) {
        std::cerr << "Warning: " << g_issues.short_lines
                  << " line(s) skipped: too few columns\n";
    }
    if (g_issues.extra_columns) {
        std::cerr << "Warning: " << g_issues.extra_columns
                  << " line(s) skipped: too many columns\n";
    }
    if (g_issues.bad_coordinates) {
        std::cerr << "Warning: " << g_issues.bad_coordinates
                  << " line(s) skipped: invalid coordinates\n";
    }
    if (g_issues.bad_strand) {
        std::cerr << "Warning: " << g_issues.bad_strand
                  << " line(s) skipped: invalid strand (expected +, -, ., or ?)\n";
    }
    if (g_issues.bad_phase) {
        std::cerr << "Warning: " << g_issues.bad_phase
                  << " line(s) skipped: invalid phase (expected 0, 1, 2, or .)\n";
    }
    if (g_issues.cds_without_phase) {
        std::cerr << "Warning: " << g_issues.cds_without_phase
                  << " CDS line(s) without phase 0/1/2\n";
    }
}

}  // namespace

static std::optional<std::string> extract_attr_value(std::string_view attrs, std::string_view key) {
    size_t pos = 0;
    while (pos < attrs.size()) {
        // Find the end of the current tag=value pair (next ';').
        size_t pair_end = attrs.find(';', pos);
        if (pair_end == std::string_view::npos) {
            pair_end = attrs.size();
        }
        const auto pair = attrs.substr(pos, pair_end - pos);
        pos = (pair_end < attrs.size()) ? pair_end + 1 : attrs.size();

        // Skip empty fragments (e.g. from ";;").
        if (pair.empty()) {
            continue;
        }

        const size_t eq = pair.find('=');
        if (eq == std::string_view::npos || eq == 0) {
            continue;
        }
        // Trim around the key so `ID=x; Parent=y` resolves: the attribute
        // indexer and --out-attrs trim the same way (attributes.cpp).
        if (trim_view(pair.substr(0, eq)) == key) {
            const auto value = pair.substr(eq + 1);
            if (value.empty()) {
                return std::nullopt;
            }
            // Split the RAW value on ',' before decoding: a literal comma in a
            // single value must stay escaped as %2C (GFF3 spec). Single-value
            // fields take the first part; the index builds full lists from
            // parse_attributes.
            const size_t comma = value.find(',');
            const auto first_part = (comma == std::string_view::npos)
                                        ? value
                                        : value.substr(0, comma);
            if (first_part.empty()) {
                return std::nullopt;
            }
            return url_decode(first_part);
        }
    }
    return std::nullopt;
}

int parse_file(const std::string& filename, GffData& data, InputFormat format) {
    // Bytes and gzip decoding live in input_source; record fields are views
    // into what it loaded.
    if (load_input(filename, data, format) != 0) return -1;
    // Reserve records capacity: average GFF3 line ~130 bytes. Cap to avoid
    // bad_alloc on comment-heavy files.
    size_t hint = data.content().size() / 130;
    if (hint > 1u << 20) hint = 1u << 20;
    data.reserve(hint);
    return parse_content(data, format);
}

int parse_content(GffData& data, InputFormat format) {
    std::string_view content{data.content()};
    // Strip a UTF-8 byte-order mark: it is not part of the first seqid and
    // would otherwise glue to it, or hide a ## directive from the parser.
    if (content.size() >= 3 &&
        static_cast<unsigned char>(content[0]) == 0xEF &&
        static_cast<unsigned char>(content[1]) == 0xBB &&
        static_cast<unsigned char>(content[2]) == 0xBF) {
        content.remove_prefix(3);
    }
    bool in_fasta = false;
    size_t pos = 0;
    reset_issues();
    data.source_format = format;

    while (pos < content.size()) {
        size_t eol = content.find('\n', pos);
        if (eol == std::string_view::npos) eol = content.size();
        std::string_view line = content.substr(pos, eol - pos);
        pos = eol + 1;

        if (in_fasta) continue;

        // Strip CR from CRLF line endings before any column parsing.
        if (!line.empty() && line.back() == '\r') line.remove_suffix(1);

        if (line.rfind("##FASTA", 0) == 0) { in_fasta = true; continue; }
        if (line.empty()) continue;
        // GFF3 spec: a line beginning with `>` opens an implied ##FASTA section.
        if (line[0] == '>') { in_fasta = true; continue; }
        if (line[0] == '#') {
            // Capture ## directives (##gff-version, ##sequence-region, etc.)
            // for re-emission on output. Single-line # comments are skipped.
            if (line.size() > 1 && line[1] == '#') {
                data.directives.emplace_back(line);
            }
            continue;
        }

        GffRecord rec;
        rec.line_idx = static_cast<int>(data.size());
        rec.kept = true;
        rec.src_fmt = format;

        if (format == InputFormat::GFF3 || format == InputFormat::GTF) {
            auto cols = split_line(line, '\t');
            if (cols.size() < 9) { ++g_issues.short_lines; continue; }
            // GFF3 and GTF both define exactly nine columns. A tenth field has
            // no place in the format and would be dropped on output, so the
            // row is reported instead of silently losing data.
            if (cols.size() > 9) { ++g_issues.extra_columns; continue; }

            try {
                rec.seqid = cols[0];
                rec.source = cols[1];
                rec.type = cols[2];
                rec.feat_class = classify_type(rec.type);
                {
                    // GFF3 requires positive integer coordinates, start <= end;
                    // reject trailing garbage. Invalid lines are skipped.
                    size_t pos = 0;
                    rec.start = std::stoll(std::string{cols[3]}, &pos);
                    if (pos != cols[3].size()) throw std::invalid_argument{"start"};
                    pos = 0;
                    rec.end = std::stoll(std::string{cols[4]}, &pos);
                    if (pos != cols[4].size()) throw std::invalid_argument{"end"};
                    if (rec.start < 1 || rec.end < 1 || rec.start > rec.end) {
                        throw std::invalid_argument{"coordinates"};
                    }
                }
                rec.score_raw = cols[5];
                if (cols[5] != ".") {
                    try {
                        size_t spos = 0;
                        const double s = std::stod(std::string{cols[5]}, &spos);
                        if (spos == cols[5].size()) rec.score = s;
                    } catch (...) { rec.score = std::nullopt; }
                }
                // GFF3 spec: strand is +, -, ., or ?. Phase is 0, 1, 2, or .; it is
                    // required on every CDS. A multi-character field is not
                    // a valid value, and truncating it to one character would
                    // invent data.
                    const auto valid_strand = [](std::string_view s) {
                        return s.size() == 1 &&
                               (s[0] == '+' || s[0] == '-' || s[0] == '.' || s[0] == '?');
                    };
                    const auto valid_phase = [](std::string_view s) {
                        return s.size() == 1 &&
                               (s[0] == '0' || s[0] == '1' || s[0] == '2' || s[0] == '.');
                    };
                    if (!valid_strand(cols[6])) { ++g_issues.bad_strand; continue; }
                    if (!valid_phase(cols[7])) { ++g_issues.bad_phase; continue; }
                    rec.strand = cols[6][0];
                    rec.phase = cols[7][0];
                    // GFF3 spec: phase is required for every CDS and is 0, 1 or 2.
                    if (rec.feat_class == FeatureClass::CDS &&
                        rec.phase != '0' && rec.phase != '1' && rec.phase != '2') {
                        ++g_issues.cds_without_phase;
                    }
                rec.attr_raw = cols[8];

                // The `;`/`=` scan is GFF3-only: it is not quote-aware, so on
                // GTF a quoted value containing "; Parent=..." would be read
                // as a real Parent key. GTF col9 goes through the quote-aware
                // GTF parser instead, which accepts both `key "value";` and
                // bare `key=value` fragments (mixed-format files).
                if (format == InputFormat::GTF) {
                    rec.gene_id = extract_quoted_value(rec.attr_raw, "gene_id");
                    rec.transcript_id = extract_quoted_value(rec.attr_raw, "transcript_id");
                    rec.id = extract_quoted_value(rec.attr_raw, "ID");
                    rec.parent_id = extract_quoted_value(rec.attr_raw, "Parent");
                    // A quoted GTF Parent may carry a comma list. The scalar
                    // parent_id holds one ID, the way the GFF3 scalar holds the
                    // first Parent value; parent_ids() expands the full list.
                    if (rec.parent_id) {
                        const auto comma = rec.parent_id->find(',');
                        if (comma != std::string::npos) {
                            rec.parent_id = trim_copy(rec.parent_id->substr(0, comma));
                        }
                    }
                } else {
                    rec.id = extract_attr_value(rec.attr_raw, "ID");
                    rec.parent_id = extract_attr_value(rec.attr_raw, "Parent");
                    if (rec.attr_raw.find("gene_id=") != std::string::npos) {
                        rec.gene_id = extract_attr_value(rec.attr_raw, "gene_id");
                    }
                    if (rec.attr_raw.find("transcript_id=") != std::string::npos) {
                        rec.transcript_id = extract_attr_value(rec.attr_raw, "transcript_id");
                    }
                    // Multi-parent Parent=tx1,tx2: extract_attr_value already keeps
                    // only the first raw comma part; the index builds the full list
                    // from parse_attributes.
                }

                if (format == InputFormat::GTF) {
                    // GTF has no ID=/Parent= attributes; synthesize from
                    // gene_id/transcript_id so the index can build parent/child links.
                    // Guard against empty-string optionals (e.g. Ensembl transcript_id "").
                    if (!rec.id || rec.id->empty()) {
                        if (rec.feat_class == FeatureClass::Gene) {
                            if (rec.gene_id && !rec.gene_id->empty()) {
                                rec.id = rec.gene_id;
                            } else {
                                rec.id = std::nullopt;
                            }
                        } else if (rec.feat_class == FeatureClass::Transcript) {
                            if (rec.transcript_id && !rec.transcript_id->empty()) {
                                rec.id = rec.transcript_id;
                            } else {
                                rec.id = std::nullopt;
                            }
                        } else {
                            // exon/CDS/etc: no ID — linked via parent_id only
                            rec.id = std::nullopt;
                        }
                    }
                    // A GTF transcript's hierarchy parent is its gene. A
                    // non-standard Parent attribute on a transcript row does
                    // not replace that link: doing so hides the gene from
                    // --longest and from the parent index. gene_id is the only
                    // gene reference GTF defines.
                    if (rec.feat_class == FeatureClass::Gene) {
                        rec.parent_id = std::nullopt;
                    } else if (rec.feat_class == FeatureClass::Transcript) {
                        rec.parent_id = (rec.gene_id && !rec.gene_id->empty())
                                            ? rec.gene_id
                                            : std::nullopt;
                    } else if (!rec.parent_id || rec.parent_id->empty()) {
                        if (rec.transcript_id && !rec.transcript_id->empty()) {
                            rec.parent_id = rec.transcript_id;
                        } else {
                            rec.parent_id = std::nullopt;
                        }
                    }
                }
            } catch (const std::exception&) {
                ++g_issues.bad_coordinates;
                continue;
            }
        } else if (format == InputFormat::BED) {
            auto cols = split_line(line, '\t');
            if (cols.size() < 3) { ++g_issues.short_lines; continue; }
            if (cols.size() > 12) { ++g_issues.extra_columns; continue; }

            try {
                rec.seqid = cols[0];
                {
                    size_t pos = 0;
                    const int64_t bed_start = std::stoll(std::string{cols[1]}, &pos);
                    if (pos != cols[1].size()) throw std::invalid_argument{"start"};
                    pos = 0;
                    rec.end = std::stoll(std::string{cols[2]}, &pos);
                    if (pos != cols[2].size()) throw std::invalid_argument{"end"};
                    if (bed_start < 0 || rec.end < bed_start) {
                        throw std::invalid_argument{"coordinates"};
                    }
                    // Zero-length BED interval [s,s) is an insertion site;
                    // map it to the one-base 1-based position s+1 so that
                    // end == start (GFF3 zero-length convention).
                    if (rec.end == bed_start) rec.end = bed_start + 1;
                    rec.start = bed_start + 1;  // BED 0-based half-open -> 1-based inclusive
                }
                rec.source = "gffsub";  // string literal: static storage, no allocation
                rec.type = kRegionType;
                rec.feat_class = FeatureClass::Region;
                rec.score_raw = (cols.size() > 4) ? cols[4] : ".";
                if (cols.size() > 4 && cols[4] != ".") {
                    try {
                        size_t spos = 0;
                        const double s = std::stod(std::string{cols[4]}, &spos);
                        if (spos == cols[4].size()) rec.score = s;
                    } catch (...) { rec.score = std::nullopt; }
                }
                rec.strand = (cols.size() > 5 && (cols[5] == "+" || cols[5] == "-" || cols[5] == "."))
                                 ? cols[5][0]
                                 : '.';
                if (cols.size() > 5 && !cols[5].empty() &&
                    cols[5] != "+" && cols[5] != "-" && cols[5] != ".") {
                    ++g_issues.bad_strand;
                }
                rec.phase = '.';
                rec.id = (cols.size() > 3 && !cols[3].empty()) ? std::optional<std::string>(cols[3]) : std::nullopt;
                rec.parent_id = std::nullopt;
                rec.gene_id = std::nullopt;
                rec.transcript_id = std::nullopt;
                rec.attr_raw = {};  // synthesized as "ID=..." at print_gff3 time
            } catch (const std::exception&) {
                ++g_issues.bad_coordinates;
                continue;
            }
        }

        data.append(std::move(rec));
    }

    report_issues();
    return 0;
}

int parse_stdin(GffData& data, InputFormat& format_out) {
    // Transparent gzip: `zcat file.gff3.gz | gffsub -` gets the same result
    // as passing the .gz path directly.
    if (load_stream(std::cin, data) != 0) return -1;
    const std::string_view content{data.content()};
    format_out = infer_format_from_content(content);
    size_t hint = content.size() / 130;
    if (hint > 1u << 20) hint = 1u << 20;
    data.reserve(hint);
    return parse_content(data, format_out);
}

}  // namespace gffsub

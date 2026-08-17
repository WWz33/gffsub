#include "annotation.hpp"
#include "feature_types.hpp"
#include "parser.hpp"
#include "record.hpp"
#include "region.hpp"
#include "gtf_parser.hpp"
#include <fstream>

namespace gffsub {

static std::string url_decode(std::string_view input) {
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
        const auto found_key = pair.substr(0, eq);
        if (found_key == key) {
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

static std::vector<std::string> split_line(const std::string& line, char delimiter) {
    std::vector<std::string> cols;
    cols.reserve(delimiter == '\t' ? 9 : 4);
    size_t start = 0;
    while (true) {
        auto pos = line.find(delimiter, start);
        if (pos == std::string::npos) {
            cols.emplace_back(line.substr(start));
            break;
        }
        cols.emplace_back(line.substr(start, pos - start));
        start = pos + 1;
    }
    return cols;
}

int parse_file(const std::string& filename, GffData& data, InputFormat format) {
    std::ifstream file(filename);
    if (!file.is_open()) return -1;

    // Reserve records capacity from file size to avoid ~20 reallocations
    // when pushing 1M+ records. Average GFF3 line is ~130 bytes. Only run
    // the seek optimization on seekable regular files; on a FIFO/pipe the
    // seeks fail and would leave failbit set, silently zeroing the output.
    {
        const auto pos = file.tellg();
        file.seekg(0, std::ios::end);
        const auto end = file.tellg();
        if (file && end > pos) {
            // Cap to avoid bad_alloc on sparse / comment-heavy files where the
            // byte/record ratio is far lower than 130.
            size_t hint = static_cast<size_t>(end - pos) / 130;
            if (hint > 1u << 20) hint = 1u << 20;
            data.reserve(hint);
            file.seekg(pos);  // restore (only valid when seeks succeeded)
        } else {
            file.clear();  // non-seekable input (FIFO/pipe): drop failbit, stay at pos 0
        }
    }

    std::string line;
    bool in_fasta = false;

    while (std::getline(file, line)) {
        if (in_fasta) continue;

        // Strip CR from CRLF line endings before any column parsing.
        if (!line.empty() && line.back() == '\r') line.pop_back();

        if (line.rfind("##FASTA", 0) == 0) { in_fasta = true; continue; }
        if (line.empty()) continue;
        if (line[0] == '#') {
            // Capture ## directives (##gff-version, ##sequence-region, etc.)
            // for re-emission on output. Single-line # comments are skipped.
            if (line.size() > 1 && line[1] == '#') {
                data.directives.push_back(line);
            }
            continue;
        }

        GffRecord rec;
        rec.line_idx = static_cast<int>(data.size());
        rec.kept = true;
        rec.src_fmt = format;

        if (format == InputFormat::GFF3 || format == InputFormat::GTF) {
            auto cols = split_line(line, '\t');
            if (cols.size() < 9) continue;

            try {
                rec.seqid = std::move(cols[0]);
                rec.source = std::move(cols[1]);
                rec.type = std::move(cols[2]);
                rec.feat_class = classify_type(rec.type);
                {
                    // GFF3 requires positive integer coordinates, start <= end;
                    // reject trailing garbage. Invalid lines are skipped.
                    size_t pos = 0;
                    rec.start = std::stoll(cols[3], &pos);
                    if (pos != cols[3].size()) throw std::invalid_argument{"start"};
                    pos = 0;
                    rec.end = std::stoll(cols[4], &pos);
                    if (pos != cols[4].size()) throw std::invalid_argument{"end"};
                    if (rec.start < 1 || rec.end < 1 || rec.start > rec.end) {
                        throw std::invalid_argument{"coordinates"};
                    }
                }
                rec.score_raw = cols[5];
                if (cols[5] != ".") {
                    try {
                        size_t spos = 0;
                        const double s = std::stod(cols[5], &spos);
                        if (spos == cols[5].size()) rec.score = s;
                    } catch (...) { rec.score = std::nullopt; }
                }
                rec.strand = cols[6].empty() ? '.' : cols[6][0];
                rec.phase = cols[7].empty() ? '.' : cols[7][0];
                rec.attr_raw = std::move(cols[8]);

                // Extract ID and Parent from col9. gene_id and transcript_id
                // are GTF conventions but also appear in some GFF3 files, so
                // they are only scanned when the key is present in col9 (quick
                // substring check avoids a full col9 walk when absent).
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

                if (format == InputFormat::GTF) {
                    apply_gtf_attributes(rec);
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
                    if (!rec.parent_id || rec.parent_id->empty()) {
                        if (rec.feat_class == FeatureClass::Gene) {
                            rec.parent_id = std::nullopt;
                        } else if (rec.feat_class == FeatureClass::Transcript &&
                                   rec.gene_id && !rec.gene_id->empty()) {
                            rec.parent_id = rec.gene_id;
                        } else if (rec.transcript_id && !rec.transcript_id->empty()) {
                            rec.parent_id = rec.transcript_id;
                        } else {
                            rec.parent_id = std::nullopt;
                        }
                    }
                }
            } catch (const std::exception&) {
                continue;
            }
        } else if (format == InputFormat::BED) {
            auto cols = split_line(line, '\t');
            if (cols.size() < 3) continue;

            try {
                rec.seqid = cols[0];
                {
                    size_t pos = 0;
                    const int64_t bed_start = std::stoll(cols[1], &pos);
                    if (pos != cols[1].size()) throw std::invalid_argument{"start"};
                    pos = 0;
                    rec.end = std::stoll(cols[2], &pos);
                    if (pos != cols[2].size()) throw std::invalid_argument{"end"};
                    if (bed_start < 0 || rec.end <= bed_start) {
                        throw std::invalid_argument{"coordinates"};
                    }
                    rec.start = bed_start + 1;  // BED 0-based half-open -> 1-based inclusive
                }
                rec.source = "gffsub";
                rec.type = std::string(kRegionType);
                rec.feat_class = FeatureClass::Region;
                rec.score_raw = (cols.size() > 4) ? cols[4] : ".";
                if (cols.size() > 4 && cols[4] != ".") {
                    try {
                        size_t spos = 0;
                        const double s = std::stod(cols[4], &spos);
                        if (spos == cols[4].size()) rec.score = s;
                    } catch (...) { rec.score = std::nullopt; }
                }
                rec.strand = (cols.size() > 5 && !cols[5].empty()) ? cols[5][0] : '.';
                rec.phase = '.';
                rec.id = (cols.size() > 3 && !cols[3].empty()) ? std::optional<std::string>(cols[3]) : std::nullopt;
                rec.parent_id = std::nullopt;
                rec.gene_id = std::nullopt;
                rec.transcript_id = std::nullopt;
                rec.attr_raw = rec.id ? "ID=" + *rec.id : "";
            } catch (const std::exception&) {
                continue;
            }
        }

        data.append(std::move(rec));
    }

    return 0;
}

namespace {

// Detect the input format by sniffing the first non-comment feature line.
// Heuristics mirror AGAT's select_gff_format: column-9 character shape wins
// over the filename extension, so a GFF3 file renamed .gtf is still parsed as
// GFF3.
//
//   GFF3: >=8 cols, col9 contains both '=' and ';'
//   GTF:  >=8 cols, col9 contains '"' (quoted attribute values)
//   BED:  3..12 cols, no col9, cols[1] and cols[2] are integers
//   default: GFF3
InputFormat sniff_format(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) return InputFormat::GFF3;
    // On non-seekable input (FIFO/pipe) reading here would consume bytes the
    // later parse_file() re-open cannot recover; fall back to GFF3, the
    // project's primary format.
    {
        const auto pos = f.tellg();
        f.seekg(0, std::ios::end);
        if (!f) return InputFormat::GFF3;
        f.seekg(pos);
    }
    std::string line;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty() || line[0] == '#') continue;
        if (line.rfind("##FASTA", 0) == 0) break;

        const auto cols = split_line(line, '\t');
        if (cols.size() >= 9) {
            const auto& a = cols[8];
            if (a.find('=') != std::string::npos &&
                a.find(';') != std::string::npos) {
                return InputFormat::GFF3;
            }
            if (a.find('"') != std::string::npos) {
                return InputFormat::GTF;
            }
            // Col9 present but neither GFF3 nor GTF shape; fall back to GFF3
            // (lenient) rather than guessing BED from a 9-column line.
            return InputFormat::GFF3;
        }
        if (cols.size() >= 3 && cols.size() <= 12) {
            // Tentative BED: require integer start/end.
            const auto is_int = [](const std::string& s) {
                if (s.empty()) return false;
                size_t pos = 0;
                try {
                    std::stoll(s, &pos);
                } catch (...) {
                    return false;
                }
                return pos == s.size();
            };
            if (is_int(cols[1]) && is_int(cols[2])) {
                return InputFormat::BED;
            }
        }
        // Not a recognizable feature line; keep scanning.
    }
    return InputFormat::GFF3;
}

}  // namespace

InputFormat infer_input_format(const std::string& path) {
    return sniff_format(path);
}

}  // namespace gffsub

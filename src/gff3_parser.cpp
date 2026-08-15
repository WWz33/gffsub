#include "annotation.hpp"
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

int parse_file(const std::string& filename, GffData& data, IdIndex& idx, InputFormat format) {
    std::ifstream file(filename);
    if (!file.is_open()) return -1;

    std::string line;
    bool in_fasta = false;

    while (std::getline(file, line)) {
        if (in_fasta) continue;

        // Strip CR from CRLF line endings before any column parsing.
        if (!line.empty() && line.back() == '\r') line.pop_back();

        if (line.rfind("##FASTA", 0) == 0) { in_fasta = true; continue; }
        if (line.empty() || line[0] == '#') continue;

        GffRecord rec;
        rec.line_idx = static_cast<int>(data.size());
        rec.kept = true;
        rec.src_fmt = format;

        if (format == InputFormat::GFF3 || format == InputFormat::GTF) {
            auto cols = split_line(line, '\t');
            if (cols.size() < 9) continue;

            try {
                rec.seqid = cols[0];
                rec.source = cols[1];
                rec.type = cols[2];
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
                rec.attr_raw = cols[8];

                rec.id = extract_attr_value(cols[8], "ID");
                rec.parent_id = extract_attr_value(cols[8], "Parent");
                rec.gene_id = extract_attr_value(cols[8], "gene_id");
                rec.transcript_id = extract_attr_value(cols[8], "transcript_id");
                // Multi-parent Parent=tx1,tx2: extract_attr_value already keeps
                // only the first raw comma part; the index builds the full list
                // from parse_attributes.

                if (format == InputFormat::GTF) {
                    apply_gtf_attributes(rec);
                    // GTF has no ID=/Parent= attributes; synthesize from
                    // gene_id/transcript_id so the index can build parent/child links.
                    // Guard against empty-string optionals (e.g. Ensembl transcript_id "").
                    if (!rec.id || rec.id->empty()) {
                        if (rec.type == "gene") {
                            // gene: ID = gene_id
                            if (rec.gene_id && !rec.gene_id->empty()) {
                                rec.id = rec.gene_id;
                            } else {
                                rec.id = std::nullopt;
                            }
                        } else if (rec.type == "transcript" || rec.type == "mRNA") {
                            // transcript: ID = transcript_id
                            if (rec.transcript_id && !rec.transcript_id->empty()) {
                                rec.id = rec.transcript_id;
                            } else {
                                rec.id = std::nullopt;
                            }
                        } else {
                            // exon/CDS/etc: no synthesized ID — linked via parent_id only
                            rec.id = std::nullopt;
                        }
                    }
                    if (!rec.parent_id || rec.parent_id->empty()) {
                        if (rec.type == "gene") {
                            rec.parent_id = std::nullopt;
                        } else if ((rec.type == "transcript" || rec.type == "mRNA") &&
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
                rec.type = "region";
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

        if (rec.id) {
            idx.add(*rec.id, rec.line_idx);
        }

        data.append(rec);
    }

    return 0;
}

}  // namespace gffsub

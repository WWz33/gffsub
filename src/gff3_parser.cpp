#include "gff3.hpp"
#include "gtf_parser.hpp"
#include <fstream>

namespace gffsub {

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
        if (eq == std::string_view::npos || eq == 0 || eq + 1 >= pair.size()) {
            continue;
        }
        const auto found_key = pair.substr(0, eq);
        if (found_key == key) {
            return std::string(pair.substr(eq + 1));
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
    int line_num = 0;
    bool in_fasta = false;

    while (std::getline(file, line)) {
        line_num++;
        if (in_fasta) continue;

        if (line.rfind("##FASTA", 0) == 0) { in_fasta = true; continue; }
        if (line.empty() || line[0] == '#') continue;

        GffRecord rec;
        rec.line_idx = static_cast<int>(data.size());
        rec.kept = true;

        if (format == InputFormat::GFF3 || format == InputFormat::GTF) {
            auto cols = split_line(line, '\t');
            if (cols.size() < 9) continue;

            try {
                rec.seqid = cols[0];
                rec.source = cols[1];
                rec.type = cols[2];
                rec.start = std::stoll(cols[3]);
                rec.end = std::stoll(cols[4]);
                rec.score = (cols[5] == ".") ? std::nullopt : std::optional(std::stod(cols[5]));
                rec.score_raw = cols[5];
                rec.strand = cols[6].empty() ? '.' : cols[6][0];
                rec.phase = cols[7].empty() ? '.' : cols[7][0];
                rec.attr_raw = cols[8];

                rec.id = extract_attr_value(cols[8], "ID");
                rec.parent_id = extract_attr_value(cols[8], "Parent");
                rec.gene_id = extract_attr_value(cols[8], "gene_id");
                rec.transcript_id = extract_attr_value(cols[8], "transcript_id");

                if (format == InputFormat::GTF) {
                    apply_gtf_attributes(rec);
                }
            } catch (const std::exception&) {
                return -1;
            }
        } else if (format == InputFormat::BED) {
            auto cols = split_line(line, '\t');
            if (cols.size() < 3) continue;

            try {
                rec.seqid = cols[0];
                rec.start = std::stoll(cols[1]) + 1;
                rec.end = std::stoll(cols[2]);
                rec.source = "gffsub";
                rec.type = "region";
                rec.score = cols.size() > 4 ? std::optional(std::stod(cols[4])) : std::nullopt;
                rec.score_raw = (cols.size() > 4) ? cols[4] : ".";
                rec.strand = cols.size() > 5 ? (cols[5][0]) : '.';
                rec.phase = '.';
                rec.id = cols.size() > 3 ? std::optional(cols[3]) : std::nullopt;
                rec.parent_id = std::nullopt;
                rec.gene_id = std::nullopt;
                rec.transcript_id = std::nullopt;
                rec.attr_raw = rec.id ? "ID=" + *rec.id : "";
            } catch (const std::exception&) {
                return -1;
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

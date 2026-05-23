#include "gff3.hpp"
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <getopt.h>

using namespace gffsub;

static void usage(const char* prog) {
    std::cerr << "Usage: " << prog << " <input.gff3> [options]\n"
        << "\n"
        << "Selector Options:\n"
        << "  --id ID\n"
        << "      Extract a feature by exact GFF3 ID. May be repeated.\n"
        << "      Default GFF3 output matches: " << prog << " query <input.gff3> --id ID\n"
        << "  --id-list FILE\n"
        << "      Extract exact feature IDs listed one per non-empty line.\n"
        << "      Default GFF3 output matches: " << prog << " query <input.gff3> --id-list FILE\n"
        << "  --name NAME\n"
        << "      Extract a gene by ID/Name/gene_id/locus_tag/Alias/Dbxref.\n"
        << "      Default GFF3 output matches: " << prog << " query <input.gff3> --name NAME\n"
        << "  --attr KEY=VALUE\n"
        << "      Extract features by an exact GFF3 attribute value. May be repeated.\n"
        << "      Default GFF3 output matches: " << prog << " query <input.gff3> --attr KEY=VALUE\n"
        << "  --include-children\n"
        << "      Include descendants of records matched by --id, --id-list, --name, or --attr.\n"
        << "  --attrs KEYS\n"
        << "      Output selected attributes as TSV/JSON fields. Implies query summary output.\n"
        << "  --summary-format FMT\n"
        << "      Output query summary instead of GFF3. Choices: tsv, json.\n"
        << "  --upstream N, --downstream N, --strand-aware\n"
        << "      Extract a window around --id. Output matches the window command.\n"
        << "  --qc\n"
        << "      Run annotation QC. Output matches the qc command.\n"
        << "\n"
        << "Input/Region Options:\n"
        << "  -r, --region CHR:START-END\n"
        << "      Extract features overlapping the specified genomic region.\n"
        << "      Coordinates are 1-based and inclusive (GFF format).\n"
        << "      Example: -r chr1:1000000-2000000\n"
        << "\n"
        << "  -b, --bed FILE\n"
        << "      Extract features using genomic regions from a BED file.\n"
        << "      BED files use 0-based half-open coordinates, automatically\n"
        << "      converted to 1-based for internal processing.\n"
        << "\n"
        << "Feature Filter Options:\n"
        << "  -f, --feature TYPE, --type TYPE\n"
        << "      Filter features by type (3rd column in GFF/GTF).\n"
        << "      Examples: gene, mRNA, exon, CDS, transcript\n"
        << "\n"
        << "  -L, --longest\n"
        << "      Keep only the longest transcript isoform for each gene.\n"
        << "      Per-gene comparison (AGAT logic): if gene has CDS isoforms,\n"
        << "      only compare by CDS length; otherwise compare by exon length.\n"
        << "\n"
        << "  -@, --threads N\n"
        << "      Number of threads for parallel processing (default: 1).\n"
        << "      Currently used with --longest for multi-chromosome parallelization.\n"
        << "\n"
        << "Output Options:\n"
        << "  -t, --output-format FMT\n"
        << "      Output format. Choices: gff3, gtf2, gtf3, bed\n"
        << "      gff3  - GFF3 format (##gff-version 3)\n"
        << "      gtf2  - GTF2 format (##gtf-version 2)\n"
        << "      gtf3  - GTF3/Ensembl format (##gtf-version 2.2.1)\n"
        << "      bed   - BED format (0-based half-open coordinates)\n"
        << "      Default: gff3\n"
        << "\n"
        << "  -o, --output FILE\n"
        << "      Output file path. If not specified, writes to stdout.\n"
        << "\n"
        << "  -h, --help\n"
        << "      Display this help message.\n"
        << "\n"
        << "Examples:\n"
        << "  " << prog << " annotation.gff3 --id GeneA\n"
        << "  " << prog << " annotation.gff3 --id-list genes.txt\n"
        << "  " << prog << " annotation.gff3 --id GeneA --include-children\n"
        << "  " << prog << " annotation.gff3 --id GeneA --summary-format tsv\n"
        << "  " << prog << " annotation.gff3 --id GeneA --upstream 2000 --downstream 500\n"
        << "  " << prog << " annotation.gff3 --qc\n"
        << "  " << prog << " annotation.gff3 --name ABC1\n"
        << "  " << prog << " annotation.gff3 --attr biotype=protein_coding\n"
        << "  " << prog << " annotation.gff3 -r chr1:1-100000 -f gene\n"
        << "  " << prog << " annotation.gff3 --bed regions.bed -f exon\n"
        << "  " << prog << " annotation.gff3 --longest\n"
        << "  " << prog << " annotation.gff3 --longest -@ 6\n"
        << "  " << prog << " annotation.gff3 -r chr1:1-100000 -t gtf3 -o out.gtf\n";
}

static void query_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " query <input.gff3> [options]\n"
        << "\n"
        << "Query Options:\n"
        << "  --id ID                 Query a feature by ID.\n"
        << "  --name NAME             Query a gene by Name/Alias/gene_id/locus_tag/ID.\n"
        << "  --id-list FILE          Query one feature ID per line.\n"
        << "  --region CHR:START-END  Query features overlapping a 1-based inclusive region.\n"
        << "  --type TYPE             Restrict query output by feature type.\n"
        << "  --attr KEY=VALUE        Query features by an exact GFF3 attribute value.\n"
        << "  --attrs KEYS            Output selected attributes as extra TSV/JSON fields.\n"
        << "  --include-children      Include descendants of matched IDs.\n"
        << "  --summary-format FMT    Output query summary instead of GFF3. Choices: tsv, json.\n"
        << "  -h, --help              Display this help message.\n";
}

static void window_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " window <input.gff3> [options]\n"
        << "\n"
        << "Window Options:\n"
        << "  --id ID                 Target feature ID or gene lookup key.\n"
        << "  --upstream N            Bases to add upstream of the target (default: 0).\n"
        << "  --downstream N          Bases to add downstream of the target (default: 0).\n"
        << "  --strand-aware          Interpret upstream/downstream by feature strand.\n"
        << "  -h, --help              Display this help message.\n";
}

static void qc_usage(const char* prog) {
    std::cerr << "Usage: " << prog << " qc <input.gff3>\n"
        << "\n"
        << "QC checks:\n"
        << "  duplicate_id      Repeated ID attributes.\n"
        << "  invalid_range     start greater than end.\n"
        << "  missing_parent    Parent points to an absent ID.\n"
        << "  child_outside_parent  Child coordinates outside parent coordinates.\n";
}

struct SummaryRow {
    std::string query_id;
    std::string matched_id;
    std::string matched_by;
    std::string seqid;
    int64_t start = 0;
    int64_t end = 0;
    char strand = '.';
    std::string type;
    std::string parent_id;
    size_t child_count = 0;
    size_t transcript_count = 0;
    size_t exon_count = 0;
    int64_t cds_length = 0;
    std::string status;
    std::vector<std::string> attrs;
};

static bool append_unique(GffData& out, std::unordered_set<int>& seen, const GffRecord& rec) {
    if (!seen.insert(rec.line_idx).second) {
        return false;
    }
    out.append(rec);
    return true;
}

static std::string record_id(const GffRecord& rec) {
    if (rec.id) return *rec.id;
    if (rec.gene_id) return *rec.gene_id;
    if (rec.transcript_id) return *rec.transcript_id;
    return "";
}

static void add_feature_counts(SummaryRow& row, const std::vector<GffRecord>& records) {
    for (const auto& rec : records) {
        if (rec.type == "mRNA" || rec.type == "transcript") {
            ++row.transcript_count;
        } else if (rec.type == "exon") {
            ++row.exon_count;
        } else if (rec.type == "CDS") {
            row.cds_length += rec.end - rec.start + 1;
        }
    }
}

static SummaryRow make_summary_row(const AnnotationIndex& index,
                                   const std::string& query_id,
                                   const std::string& matched_by,
                                   const GffRecord& rec) {
    SummaryRow row;
    row.query_id = query_id;
    row.matched_id = record_id(rec);
    row.matched_by = matched_by;
    row.seqid = rec.seqid;
    row.start = rec.start;
    row.end = rec.end;
    row.strand = rec.strand;
    row.type = rec.type;
    row.parent_id = rec.parent_id.value_or("");
    row.status = "found";

    if (rec.id) {
        const auto children = index.children_of(*rec.id);
        row.child_count = children.size();
        const auto model = index.gene_model(*rec.id);
        if (model) {
            add_feature_counts(row, model->records);
        } else {
            add_feature_counts(row, children);
        }
    }

    return row;
}

static SummaryRow make_not_found_row(const std::string& query_id, const std::string& matched_by) {
    SummaryRow row;
    row.query_id = query_id;
    row.matched_by = matched_by;
    row.status = "not_found";
    return row;
}

static bool contains_record(const std::vector<GffRecord>& records, int line_idx) {
    for (const auto& rec : records) {
        if (rec.line_idx == line_idx) {
            return true;
        }
    }
    return false;
}

static std::string infer_gene_match_key(const AnnotationIndex& index, const std::string& query, const GffRecord& rec) {
    if (rec.id && *rec.id == query) {
        return "ID";
    }
    if (rec.gene_id && *rec.gene_id == query) {
        return "gene_id";
    }
    for (const char* key : {"Name", "locus_tag", "Alias", "Dbxref"}) {
        if (contains_record(index.with_attribute(key, query), rec.line_idx)) {
            return key;
        }
    }
    return "name";
}

static std::string json_escape(const std::string& value) {
    std::ostringstream out;
    for (const char ch : value) {
        switch (ch) {
            case '\\': out << "\\\\"; break;
            case '"': out << "\\\""; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            default: out << ch; break;
        }
    }
    return out.str();
}

static std::string trim_copy(std::string_view value) {
    size_t start = 0;
    while (start < value.size() && std::isspace(static_cast<unsigned char>(value[start]))) {
        ++start;
    }

    size_t end = value.size();
    while (end > start && std::isspace(static_cast<unsigned char>(value[end - 1]))) {
        --end;
    }

    return std::string{value.substr(start, end - start)};
}

static std::vector<std::string> split_attr_keys(std::string_view keys) {
    std::vector<std::string> result;
    size_t pos = 0;
    while (pos <= keys.size()) {
        const size_t comma = keys.find(',', pos);
        const size_t end = (comma == std::string_view::npos) ? keys.size() : comma;
        const auto key = trim_copy(keys.substr(pos, end - pos));
        if (!key.empty()) {
            result.push_back(key);
        }
        if (comma == std::string_view::npos) {
            break;
        }
        pos = comma + 1;
    }
    return result;
}

static std::string join_values(const std::vector<std::string>& values) {
    std::ostringstream out;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i > 0) {
            out << ',';
        }
        out << values[i];
    }
    return out.str();
}

static std::vector<std::string> extract_selected_attrs(const std::string& attrs,
                                                       const std::vector<std::string>& keys) {
    const auto parsed = parse_attributes(attrs);
    std::vector<std::string> values;
    values.reserve(keys.size());
    for (const auto& key : keys) {
        const auto it = parsed.find(key);
        if (it == parsed.end()) {
            values.emplace_back();
        } else {
            values.push_back(join_values(it->second));
        }
    }
    return values;
}

static void print_qc_row(std::ostream& out,
                         const char* severity,
                         const char* code,
                         int line_idx,
                         const std::string& id,
                         const std::string& message) {
    out << severity << '\t'
        << code << '\t'
        << line_idx << '\t'
        << id << '\t'
        << message << '\n';
}

static void print_summary_tsv(std::ostream& out,
                              const std::vector<SummaryRow>& rows,
                              const std::vector<std::string>& selected_attrs) {
    out << "query_id\tmatched_id\tmatched_by\tseqid\tstart\tend\tstrand\ttype\tparent_id\t"
        << "child_count\ttranscript_count\texon_count\tcds_length\tstatus";
    for (const auto& key : selected_attrs) {
        out << '\t' << key;
    }
    out << '\n';
    for (const auto& row : rows) {
        out << row.query_id << '\t'
            << row.matched_id << '\t'
            << row.matched_by << '\t'
            << row.seqid << '\t'
            << row.start << '\t'
            << row.end << '\t'
            << row.strand << '\t'
            << row.type << '\t'
            << row.parent_id << '\t'
            << row.child_count << '\t'
            << row.transcript_count << '\t'
            << row.exon_count << '\t'
            << row.cds_length << '\t'
            << row.status;
        for (size_t i = 0; i < selected_attrs.size(); ++i) {
            out << '\t';
            if (i < row.attrs.size()) {
                out << row.attrs[i];
            }
        }
        out << '\n';
    }
}

static void print_summary_json(std::ostream& out,
                               const std::vector<SummaryRow>& rows,
                               const std::vector<std::string>& selected_attrs) {
    out << "[\n";
    for (size_t i = 0; i < rows.size(); ++i) {
        const auto& row = rows[i];
        out << "  {"
            << "\"query_id\":\"" << json_escape(row.query_id) << "\","
            << "\"matched_id\":\"" << json_escape(row.matched_id) << "\","
            << "\"matched_by\":\"" << json_escape(row.matched_by) << "\","
            << "\"seqid\":\"" << json_escape(row.seqid) << "\","
            << "\"start\":" << row.start << ','
            << "\"end\":" << row.end << ','
            << "\"strand\":\"" << row.strand << "\","
            << "\"type\":\"" << json_escape(row.type) << "\","
            << "\"parent_id\":\"" << json_escape(row.parent_id) << "\","
            << "\"child_count\":" << row.child_count << ','
            << "\"transcript_count\":" << row.transcript_count << ','
            << "\"exon_count\":" << row.exon_count << ','
            << "\"cds_length\":" << row.cds_length << ','
            << "\"status\":\"" << row.status << "\"";
        if (!selected_attrs.empty()) {
            out << ",\"attrs\":{";
            for (size_t j = 0; j < selected_attrs.size(); ++j) {
                out << "\"" << json_escape(selected_attrs[j]) << "\":\"";
                if (j < row.attrs.size()) {
                    out << json_escape(row.attrs[j]);
                }
                out << "\"";
                if (j + 1 < selected_attrs.size()) {
                    out << ',';
                }
            }
            out << "}";
        }
        out << "}";
        if (i + 1 < rows.size()) {
            out << ',';
        }
        out << '\n';
    }
    out << "]\n";
}

static int run_query(int argc, char* argv[], const char* prog) {
    if (argc == 2 && (std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")) {
        query_usage(prog);
        return 0;
    }

    if (argc < 2) {
        query_usage(prog);
        return 1;
    }

    const std::string input_file = argv[1];
    std::vector<std::string> ids;
    std::string name;
    std::string id_list_file;
    std::string region_str;
    std::string feature_type;
    std::string summary_format;
    std::vector<std::string> selected_attrs;
    std::vector<std::pair<std::string, std::string>> attr_filters;
    bool include_children = false;

    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        auto require_value = [&](const char* option) -> std::optional<std::string> {
            if (i + 1 >= argc) {
                std::cerr << "Error: " << option << " requires a value\n";
                return std::nullopt;
            }
            ++i;
            return std::string{argv[i]};
        };

        if (arg == "--id") {
            auto value = require_value("--id");
            if (!value) return 1;
            ids.push_back(*value);
        } else if (arg == "--name") {
            auto value = require_value("--name");
            if (!value) return 1;
            name = *value;
        } else if (arg == "--id-list") {
            auto value = require_value("--id-list");
            if (!value) return 1;
            id_list_file = *value;
        } else if (arg == "--region") {
            auto value = require_value("--region");
            if (!value) return 1;
            region_str = *value;
        } else if (arg == "--type") {
            auto value = require_value("--type");
            if (!value) return 1;
            feature_type = *value;
        } else if (arg == "--attr") {
            auto value = require_value("--attr");
            if (!value) return 1;
            const auto equal_pos = value->find('=');
            if (equal_pos == std::string::npos || equal_pos == 0 || equal_pos + 1 == value->size()) {
                std::cerr << "Error: --attr expects KEY=VALUE\n";
                return 1;
            }
            attr_filters.emplace_back(value->substr(0, equal_pos), value->substr(equal_pos + 1));
        } else if (arg == "--attrs") {
            auto value = require_value("--attrs");
            if (!value) return 1;
            const auto keys = split_attr_keys(*value);
            if (keys.empty()) {
                std::cerr << "Error: --attrs expects a comma-separated list of keys\n";
                return 1;
            }
            selected_attrs.insert(selected_attrs.end(), keys.begin(), keys.end());
        } else if (arg == "--summary-format") {
            auto value = require_value("--summary-format");
            if (!value) return 1;
            summary_format = *value;
            if (summary_format != "tsv" && summary_format != "json") {
                std::cerr << "Error: --summary-format expects tsv or json\n";
                return 1;
            }
        } else if (arg == "--include-children") {
            include_children = true;
        } else if (arg == "-h" || arg == "--help") {
            query_usage(prog);
            return 0;
        } else {
            std::cerr << "Error: unknown query option " << arg << '\n';
            query_usage(prog);
            return 1;
        }
    }

    if (!id_list_file.empty()) {
        std::ifstream in{id_list_file};
        if (!in.is_open()) {
            std::cerr << "Error: cannot open " << id_list_file << '\n';
            return 1;
        }
        std::string line;
        while (std::getline(in, line)) {
            if (!line.empty()) {
                ids.push_back(line);
            }
        }
    }

    gffsub::AnnotationIndex index = gffsub::AnnotationIndex::from_gff3(input_file);
    GffData result;
    std::unordered_set<int> seen;
    std::vector<SummaryRow> summary_rows;
    const bool emit_summary = !summary_format.empty() || !selected_attrs.empty();

    auto add_summary = [&](const std::string& query_id, const std::string& matched_by, const GffRecord& rec) {
        if (emit_summary) {
            auto row = make_summary_row(index, query_id, matched_by, rec);
            if (!selected_attrs.empty()) {
                row.attrs = extract_selected_attrs(rec.attr_raw, selected_attrs);
            }
            summary_rows.push_back(std::move(row));
        }
    };

    auto add_match = [&](const GffRecord& rec, const std::string& query_id, const std::string& matched_by) {
        if (feature_type.empty() || rec.type == feature_type) {
            if (append_unique(result, seen, rec)) {
                add_summary(query_id, matched_by, rec);
            }
        }
        if (include_children && rec.id) {
            for (const auto& child : index.descendants_of(*rec.id)) {
                if (feature_type.empty() || child.type == feature_type) {
                    if (append_unique(result, seen, child)) {
                        add_summary(query_id, "child", child);
                    }
                }
            }
        }
    };

    for (const auto& id : ids) {
        const auto rec = index.find_by_id(id);
        if (rec) {
            add_match(*rec, id, "ID");
        } else if (emit_summary) {
            auto row = make_not_found_row(id, "ID");
            row.attrs.assign(selected_attrs.size(), "");
            summary_rows.push_back(std::move(row));
        }
    }

    if (!name.empty()) {
        const auto rec = index.find_gene(name);
        if (rec) {
            add_match(*rec, name, infer_gene_match_key(index, name, *rec));
        } else if (emit_summary) {
            auto row = make_not_found_row(name, "name");
            row.attrs.assign(selected_attrs.size(), "");
            summary_rows.push_back(std::move(row));
        }
    }

    if (!region_str.empty()) {
        const auto region = parse_region(region_str);
        if (!region) {
            std::cerr << "Error: invalid region format " << region_str << '\n';
            return 1;
        }
        for (const auto& rec : index.overlap(region->seqid, region->start, region->end)) {
            add_match(rec, region_str, "region");
        }
    }

    for (const auto& [key, value] : attr_filters) {
        bool matched = false;
        for (const auto& rec : index.with_attribute(key, value)) {
            matched = true;
            add_match(rec, key + "=" + value, key);
        }
        if (!matched && emit_summary) {
            auto row = make_not_found_row(key + "=" + value, key);
            row.attrs.assign(selected_attrs.size(), "");
            summary_rows.push_back(std::move(row));
        }
    }

    std::sort(result.records.begin(), result.records.end(),
              [](const GffRecord& lhs, const GffRecord& rhs) {
                  return lhs.line_idx < rhs.line_idx;
              });
    if (summary_format == "json") {
        print_summary_json(std::cout, summary_rows, selected_attrs);
    } else if (emit_summary) {
        print_summary_tsv(std::cout, summary_rows, selected_attrs);
    } else {
        print_gff3(std::cout, result);
    }
    return 0;
}

static int run_window(int argc, char* argv[], const char* prog) {
    if (argc < 2) {
        window_usage(prog);
        return 1;
    }

    const std::string input_file = argv[1];
    std::string id;
    int64_t upstream = 0;
    int64_t downstream = 0;
    bool strand_aware = false;

    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        auto require_value = [&](const char* option) -> std::optional<std::string> {
            if (i + 1 >= argc) {
                std::cerr << "Error: " << option << " requires a value\n";
                return std::nullopt;
            }
            ++i;
            return std::string{argv[i]};
        };

        if (arg == "--id") {
            auto value = require_value("--id");
            if (!value) return 1;
            if (!id.empty()) {
                std::cerr << "Error: window requires exactly one --id\n";
                return 1;
            }
            id = *value;
        } else if (arg == "--upstream") {
            auto value = require_value("--upstream");
            if (!value) return 1;
            upstream = std::stoll(*value);
            if (upstream < 0) {
                std::cerr << "Error: --upstream must be non-negative\n";
                return 1;
            }
        } else if (arg == "--downstream") {
            auto value = require_value("--downstream");
            if (!value) return 1;
            downstream = std::stoll(*value);
            if (downstream < 0) {
                std::cerr << "Error: --downstream must be non-negative\n";
                return 1;
            }
        } else if (arg == "--strand-aware") {
            strand_aware = true;
        } else if (arg == "-h" || arg == "--help") {
            window_usage(prog);
            return 0;
        } else {
            std::cerr << "Error: unknown window option " << arg << '\n';
            window_usage(prog);
            return 1;
        }
    }

    if (id.empty()) {
        std::cerr << "Error: window requires --id\n";
        return 1;
    }

    const auto index = gffsub::AnnotationIndex::from_gff3(input_file);
    auto target = index.find_by_id(id);
    if (!target) {
        target = index.find_gene(id);
    }
    if (!target) {
        std::cerr << "Error: cannot find " << id << '\n';
        return 1;
    }

    const auto region = window_region(*target, upstream, downstream, strand_aware);
    GffData result;
    for (const auto& rec : index.overlap(region.seqid, region.start, region.end)) {
        result.append(rec);
    }
    print_gff3(std::cout, result);
    return 0;
}

static int run_qc(int argc, char* argv[], const char* prog) {
    if (argc != 2) {
        qc_usage(prog);
        return 1;
    }

    GffData data;
    IdIndex idx;
    const std::string input_file = argv[1];
    if (parse_file(input_file, data, idx, InputFormat::GFF3) != 0) {
        std::cerr << "Error: cannot parse " << input_file << '\n';
        return 1;
    }

    std::unordered_map<std::string, const GffRecord*> by_id;
    std::unordered_map<std::string, int> id_counts;
    for (const auto& rec : data.records) {
        if (rec.id) {
            ++id_counts[*rec.id];
            by_id.emplace(*rec.id, &rec);
        }
    }

    std::cout << "severity\tcode\tline_idx\tid\tmessage\n";

    for (const auto& [id, count] : id_counts) {
        if (count > 1) {
            print_qc_row(std::cout, "error", "duplicate_id", -1, id, "ID appears more than once");
        }
    }

    for (const auto& rec : data.records) {
        const std::string id = record_id(rec);
        if (rec.start > rec.end) {
            print_qc_row(std::cout, "error", "invalid_range", rec.line_idx, id, "start is greater than end");
        }

        const auto attrs = parse_attributes(rec.attr_raw);
        const auto parent_it = attrs.find("Parent");
        if (parent_it == attrs.end()) {
            continue;
        }
        for (const auto& parent_id : parent_it->second) {
            const auto parent_record_it = by_id.find(parent_id);
            if (parent_record_it == by_id.end()) {
                print_qc_row(std::cout, "error", "missing_parent", rec.line_idx, id,
                             "Parent " + parent_id + " was not found");
                continue;
            }

            const auto& parent = *parent_record_it->second;
            if (rec.seqid != parent.seqid || rec.start < parent.start || rec.end > parent.end) {
                print_qc_row(std::cout, "warning", "child_outside_parent", rec.line_idx, id,
                             "child is outside Parent " + parent_id);
            }
        }
    }

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc > 1 && std::string(argv[1]) == "query") {
        return run_query(argc - 1, argv + 1, argv[0]);
    }
    if (argc > 1 && std::string(argv[1]) == "window") {
        return run_window(argc - 1, argv + 1, argv[0]);
    }
    if (argc > 1 && std::string(argv[1]) == "qc") {
        return run_qc(argc - 1, argv + 1, argv[0]);
    }

    std::vector<std::string> ids;
    std::string id_list_file;
    std::string name;
    std::vector<std::pair<std::string, std::string>> attr_filters;
    bool include_children = false;
    std::vector<std::string> selected_attrs;
    std::string summary_format;
    std::string upstream_arg;
    std::string downstream_arg;
    bool strand_aware = false;
    bool do_qc = false;
    std::string region_str;
    std::string bed_file;
    std::string feature;
    bool do_longest = false;
    size_t num_threads = 6;
    bool threads_set = false;
    std::string output_format = "gff3";
    std::string output_file;

    enum {
        OPT_ID = 1000,
        OPT_ID_LIST,
        OPT_NAME,
        OPT_ATTR,
        OPT_ATTRS,
        OPT_SUMMARY_FORMAT,
        OPT_UPSTREAM,
        OPT_DOWNSTREAM,
        OPT_STRAND_AWARE,
        OPT_QC
    };
    static struct option long_options[] = {
        {"id",            required_argument, nullptr, OPT_ID},
        {"id-list",       required_argument, nullptr, OPT_ID_LIST},
        {"name",          required_argument, nullptr, OPT_NAME},
        {"attr",          required_argument, nullptr, OPT_ATTR},
        {"attrs",         required_argument, nullptr, OPT_ATTRS},
        {"summary-format", required_argument, nullptr, OPT_SUMMARY_FORMAT},
        {"upstream",      required_argument, nullptr, OPT_UPSTREAM},
        {"downstream",    required_argument, nullptr, OPT_DOWNSTREAM},
        {"strand-aware",  no_argument,       nullptr, OPT_STRAND_AWARE},
        {"qc",            no_argument,       nullptr, OPT_QC},
        {"include-children", no_argument,     nullptr, 'C'},
        {"region",        required_argument, nullptr, 'r'},
        {"bed",           required_argument, nullptr, 'b'},
        {"feature",       required_argument, nullptr, 'f'},
        {"type",          required_argument, nullptr, 'f'},
        {"longest",       no_argument,       nullptr, 'L'},
        {"threads",       required_argument, nullptr, '@'},
        {"output-format", required_argument, nullptr, 't'},
        {"output",        required_argument, nullptr, 'o'},
        {"help",          no_argument,       nullptr, 'h'},
        {nullptr,        0,                 nullptr, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "r:b:f:L@:t:o:h", long_options, nullptr)) != -1) {
        switch (opt) {
            case OPT_ID: ids.emplace_back(optarg); break;
            case OPT_ID_LIST: id_list_file = optarg; break;
            case OPT_NAME: name = optarg; break;
            case OPT_ATTR: {
                const std::string value{optarg};
                const auto equal_pos = value.find('=');
                if (equal_pos == std::string::npos || equal_pos == 0 || equal_pos + 1 == value.size()) {
                    std::cerr << "Error: --attr expects KEY=VALUE\n";
                    return 1;
                }
                attr_filters.emplace_back(value.substr(0, equal_pos), value.substr(equal_pos + 1));
                break;
            }
            case OPT_ATTRS: {
                const auto keys = split_attr_keys(optarg);
                if (keys.empty()) {
                    std::cerr << "Error: --attrs expects a comma-separated list of keys\n";
                    return 1;
                }
                selected_attrs.insert(selected_attrs.end(), keys.begin(), keys.end());
                break;
            }
            case OPT_SUMMARY_FORMAT:
                summary_format = optarg;
                if (summary_format != "tsv" && summary_format != "json") {
                    std::cerr << "Error: --summary-format expects tsv or json\n";
                    return 1;
                }
                break;
            case 'C': include_children = true; break;
            case OPT_UPSTREAM: upstream_arg = optarg; break;
            case OPT_DOWNSTREAM: downstream_arg = optarg; break;
            case OPT_STRAND_AWARE: strand_aware = true; break;
            case OPT_QC: do_qc = true; break;
            case 'r': region_str = optarg; break;
            case 'b': bed_file = optarg; break;
            case 'f': feature = optarg; break;
            case 'L': do_longest = true; break;
            case '@': {
                threads_set = true;
                size_t t = std::stoul(optarg);
                if (t == 0) t = 1;
                if (t > 256) t = 256; // cap to prevent over-subscription
                num_threads = t;
                break;
            }
            case 't': output_format = optarg; break;
            case 'o': output_file = optarg; break;
            case 'h': usage(argv[0]); return 0;
            default: usage(argv[0]); return 1;
        }
    }

    if (optind >= argc) {
        usage(argv[0]);
        return 1;
    }

    // Validate output format
    OutputFormat fmt = OutputFormat::GFF3;
    if (output_format == "gtf") {
        fmt = OutputFormat::GTF2;
    } else if (output_format == "gtf2") {
        fmt = OutputFormat::GTF2;
    } else if (output_format == "gtf3") {
        fmt = OutputFormat::GTF3;
    } else if (output_format == "bed") {
        fmt = OutputFormat::BED;
    } else if (output_format == "gff3") {
        fmt = OutputFormat::GFF3;
    } else {
        std::cerr << "Error: unknown output format " << output_format << '\n';
        std::cerr << "Supported formats: gff3, gtf2, gtf3, bed\n";
        return 1;
    }

    std::string input_file = argv[optind];

    if (do_qc) {
        if (!ids.empty() || !id_list_file.empty() || !name.empty() || !attr_filters.empty() || include_children ||
            !selected_attrs.empty() || !summary_format.empty() || !upstream_arg.empty() || !downstream_arg.empty() ||
            strand_aware || !region_str.empty() || !bed_file.empty() || !feature.empty() || do_longest ||
            output_format != "gff3" || !output_file.empty()) {
            std::cerr << "Error: --qc only supports the input file\n";
            return 1;
        }

        std::vector<std::string> qc_args{"qc", input_file};
        std::vector<char*> qc_argv;
        qc_argv.reserve(qc_args.size());
        for (auto& arg : qc_args) {
            qc_argv.push_back(arg.data());
        }
        return run_qc(static_cast<int>(qc_argv.size()), qc_argv.data(), argv[0]);
    }

    if (!upstream_arg.empty() || !downstream_arg.empty() || strand_aware) {
        if (ids.size() != 1) {
            std::cerr << "Error: window shortcut requires exactly one --id\n";
            return 1;
        }
        if (!id_list_file.empty() || !name.empty() || !attr_filters.empty() || include_children ||
            !selected_attrs.empty() || !summary_format.empty() || !region_str.empty() || !bed_file.empty() ||
            !feature.empty() || do_longest || output_format != "gff3" || !output_file.empty()) {
            std::cerr << "Error: window shortcut only supports --id, --upstream, --downstream, and --strand-aware\n";
            return 1;
        }

        std::vector<std::string> window_args{"window", input_file};
        for (const auto& id : ids) {
            window_args.push_back("--id");
            window_args.push_back(id);
        }
        if (!upstream_arg.empty()) {
            window_args.push_back("--upstream");
            window_args.push_back(upstream_arg);
        }
        if (!downstream_arg.empty()) {
            window_args.push_back("--downstream");
            window_args.push_back(downstream_arg);
        }
        if (strand_aware) {
            window_args.push_back("--strand-aware");
        }

        std::vector<char*> window_argv;
        window_argv.reserve(window_args.size());
        for (auto& arg : window_args) {
            window_argv.push_back(arg.data());
        }
        return run_window(static_cast<int>(window_argv.size()), window_argv.data(), argv[0]);
    }

    if (!summary_format.empty() || !selected_attrs.empty()) {
        if (!bed_file.empty() || do_longest || threads_set || output_format != "gff3" || !output_file.empty()) {
            std::cerr << "Error: --summary-format/--attrs only supports query-style selectors; "
                      << "do not combine with --bed, --longest, --threads, --output-format, or --output\n";
            return 1;
        }

        std::vector<std::string> query_args{"query", input_file};
        for (const auto& id : ids) {
            query_args.push_back("--id");
            query_args.push_back(id);
        }
        if (!id_list_file.empty()) {
            query_args.push_back("--id-list");
            query_args.push_back(id_list_file);
        }
        if (!name.empty()) {
            query_args.push_back("--name");
            query_args.push_back(name);
        }
        if (!region_str.empty()) {
            query_args.push_back("--region");
            query_args.push_back(region_str);
        }
        if (!feature.empty()) {
            query_args.push_back("--type");
            query_args.push_back(feature);
        }
        for (const auto& [key, value] : attr_filters) {
            query_args.push_back("--attr");
            query_args.push_back(key + "=" + value);
        }
        if (include_children) {
            query_args.push_back("--include-children");
        }
        if (!summary_format.empty()) {
            query_args.push_back("--summary-format");
            query_args.push_back(summary_format);
        }
        if (!selected_attrs.empty()) {
            query_args.push_back("--attrs");
            query_args.push_back(join_values(selected_attrs));
        }

        std::vector<char*> query_argv;
        query_argv.reserve(query_args.size());
        for (auto& arg : query_args) {
            query_argv.push_back(arg.data());
        }
        return run_query(static_cast<int>(query_argv.size()), query_argv.data(), argv[0]);
    }

    if (!id_list_file.empty()) {
        std::ifstream in{id_list_file};
        if (!in.is_open()) {
            std::cerr << "Error: cannot open " << id_list_file << '\n';
            return 1;
        }
        std::string line;
        while (std::getline(in, line)) {
            if (!line.empty()) {
                ids.push_back(line);
            }
        }
    }

    GffData data;
    IdIndex idx;

    // Parse input file
    if (parse_file(input_file, data, idx, InputFormat::GFF3) != 0) {
        std::cerr << "Error: cannot parse " << input_file << '\n';
        return 1;
    }

    if (!ids.empty() || !name.empty() || !attr_filters.empty()) {
        const auto index = gffsub::AnnotationIndex::from_gff3(input_file);
        std::unordered_set<int> selected_lines;
        auto add_selected = [&](const GffRecord& rec) {
            if (selected_lines.insert(rec.line_idx).second && include_children && rec.id) {
                for (const auto& child : index.descendants_of(*rec.id)) {
                    selected_lines.insert(child.line_idx);
                }
            }
        };
        for (const auto& id : ids) {
            const auto rec = index.find_by_id(id);
            if (rec) {
                add_selected(*rec);
            }
        }
        if (!name.empty()) {
            const auto rec = index.find_gene(name);
            if (rec) {
                add_selected(*rec);
            }
        }
        for (const auto& [key, value] : attr_filters) {
            for (const auto& rec : index.with_attribute(key, value)) {
                add_selected(rec);
            }
        }
        for (auto& rec : data.records) {
            rec.kept = selected_lines.count(rec.line_idx) > 0;
        }
    }

    Region region{"", 0, 0};
    std::optional<Region> parsed_region;

    // Apply region filters
    if (!region_str.empty()) {
        parsed_region = parse_region(region_str);
        if (!parsed_region) {
            std::cerr << "Error: invalid region format " << region_str << '\n';
            return 1;
        }
        region = *parsed_region;
        filter_by_region(data, region);
    }

    if (!bed_file.empty()) {
        filter_by_regions_from_file(data, bed_file);
    }

    // Apply feature filters
    if (do_longest) {
        filter_longest_isoform(data, idx, feature, num_threads);
    } else if (!feature.empty()) {
        filter_by_feature(data, feature);
    }

    // Output
    std::ofstream out_file;
    std::ostream* out = &std::cout;
    if (!output_file.empty()) {
        out_file.open(output_file);
        if (!out_file.is_open()) {
            std::cerr << "Error: cannot open " << output_file << '\n';
            return 1;
        }
        out = &out_file;
    }

    switch (fmt) {
        case OutputFormat::GFF3: print_gff3(*out, data); break;
        case OutputFormat::GTF2: print_gtf(*out, data, fmt); break;
        case OutputFormat::GTF3: print_gtf3(*out, data); break;
        case OutputFormat::BED:  print_bed(*out, data); break;
    }

    return 0;
}

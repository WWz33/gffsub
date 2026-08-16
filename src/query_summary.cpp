#include "query_summary.hpp"

#include "annotation.hpp"
#include "gtf_parser.hpp"
#include "parser.hpp"

#include <cstdio>
#include <sstream>

namespace gffsub {
namespace {

void add_feature_counts(SummaryRow& row, const std::vector<GffRecord>& records) {
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

bool contains_record(const std::vector<GffRecord>& records, int line_idx) {
    for (const auto& rec : records) {
        if (rec.line_idx == line_idx) {
            return true;
        }
    }
    return false;
}

std::string json_escape(const std::string& value) {
    std::ostringstream out;
    for (const unsigned char uc : value) {
        const char ch = static_cast<char>(uc);
        switch (ch) {
            case '\\': out << "\\\\"; break;
            case '"': out << "\\\""; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            case '\b': out << "\\b"; break;
            case '\f': out << "\\f"; break;
            default:
                if (uc < 0x20) {
                    // RFC 8259: control chars must be \uXXXX escaped
                    char buf[8];
                    std::snprintf(buf, sizeof(buf), "\\u%04x", uc);
                    out << buf;
                } else {
                    out << ch;
                }
                break;
        }
    }
    return out.str();
}

std::string join_values(const std::vector<std::string>& values) {
    std::ostringstream out;
    for (size_t i = 0; i < values.size(); ++i) {
        if (i > 0) {
            out << ',';
        }
        out << values[i];
    }
    return out.str();
}

}  // namespace

std::string record_id(const GffRecord& rec) {
    if (rec.id) return *rec.id;
    if (rec.gene_id) return *rec.gene_id;
    if (rec.transcript_id) return *rec.transcript_id;
    return "";
}

SummaryRow make_summary_row(const AnnotationIndex& index,
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
        row.child_count = index.children_of(*rec.id).size();
        // All counts describe the matched feature's subtree: for a gene this
        // covers the whole model, for a transcript its own exons/CDS.
        add_feature_counts(row, index.descendants_of(*rec.id));
    }

    return row;
}

SummaryRow make_not_found_row(const std::string& query_id, const std::string& matched_by) {
    SummaryRow row;
    row.query_id = query_id;
    row.matched_by = matched_by;
    row.status = "not_found";
    return row;
}

std::string infer_gene_match_key(const AnnotationIndex& index, const std::string& query, const GffRecord& rec) {
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

std::vector<std::string> extract_output_attrs(const GffRecord& rec,
                                              const std::vector<std::string>& keys) {
    const auto parsed = parse_attributes(rec.attr_raw);
    std::vector<std::string> values;
    values.reserve(keys.size());
    for (const auto& key : keys) {
        if (key == "gene_id" && rec.gene_id) {
            values.push_back(*rec.gene_id);
            continue;
        }
        if (key == "transcript_id" && rec.transcript_id) {
            values.push_back(*rec.transcript_id);
            continue;
        }
        const auto it = parsed.find(key);
        if (it != parsed.end()) {
            values.push_back(join_values(it->second));
            continue;
        }
        // GTF fallback: `key "value";` attributes that parse_attributes cannot parse.
        const auto gtf = extract_quoted_value(rec.attr_raw, key);
        if (gtf) {
            values.push_back(*gtf);
            continue;
        }
        values.emplace_back();
    }
    return values;
}

void print_summary_tsv(std::ostream& out,
                       const std::vector<SummaryRow>& rows,
                       const std::vector<std::string>& output_attrs) {
    // TSV values cannot contain literal tabs/newlines/CR; replace them with
    // C-style escapes so decoded attribute values (e.g. %09 -> tab) don't
    // break column structure.
    auto tsv_escape = [](const std::string& s) {
        std::string out;
        out.reserve(s.size());
        for (char ch : s) {
            if (ch == '\t') out += "\\t";
            else if (ch == '\n') out += "\\n";
            else if (ch == '\r') out += "\\r";
            else out.push_back(ch);
        }
        return out;
    };
    out << "query_id\tmatched_id\tmatched_by\tseqid\tstart\tend\tstrand\ttype\tparent_id\t"
        << "child_count\ttranscript_count\texon_count\tcds_length\tstatus";
    for (const auto& key : output_attrs) {
        out << '\t' << key;
    }
    out << '\n';
    for (const auto& row : rows) {
        out << tsv_escape(row.query_id) << '\t'
            << tsv_escape(row.matched_id) << '\t'
            << tsv_escape(row.matched_by) << '\t'
            << tsv_escape(row.seqid) << '\t'
            << row.start << '\t'
            << row.end << '\t'
            << row.strand << '\t'
            << tsv_escape(row.type) << '\t'
            << tsv_escape(row.parent_id) << '\t'
            << row.child_count << '\t'
            << row.transcript_count << '\t'
            << row.exon_count << '\t'
            << row.cds_length << '\t'
            << row.status;
        for (size_t i = 0; i < output_attrs.size(); ++i) {
            out << '\t';
            if (i < row.attrs.size()) {
                out << tsv_escape(row.attrs[i]);
            }
        }
        out << '\n';
    }
}

void print_summary_json(std::ostream& out,
                        const std::vector<SummaryRow>& rows,
                        const std::vector<std::string>& output_attrs) {
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
            << "\"status\":\"" << json_escape(row.status) << "\"";
        if (!output_attrs.empty()) {
            out << ",\"attrs\":{";
            for (size_t j = 0; j < output_attrs.size(); ++j) {
                out << "\"" << json_escape(output_attrs[j]) << "\":\"";
                if (j < row.attrs.size()) {
                    out << json_escape(row.attrs[j]);
                }
                out << "\"";
                if (j + 1 < output_attrs.size()) {
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

}  // namespace gffsub

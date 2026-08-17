#include "query_summary.hpp"

#include "gtf_parser.hpp"

#include <algorithm>
#include <unordered_set>

namespace gffsub {

namespace {

std::string record_id_impl(const GffRecord& rec) {
    if (rec.id) return *rec.id;
    if (rec.gene_id) return *rec.gene_id;
    if (rec.transcript_id) return *rec.transcript_id;
    return "";
}

void add_feature_counts(SummaryRow& row, const std::vector<GffRecord>& records) {
    for (const auto& rec : records) {
        if (rec.feat_class == FeatureClass::Transcript) {
            ++row.transcript_count;
        } else if (rec.feat_class == FeatureClass::Exon) {
            ++row.exon_count;
        } else if (rec.feat_class == FeatureClass::CDS) {
            row.cds_length += rec.end - rec.start + 1;
        }
    }
}

std::string tsv_escape(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char ch : s) {
        if (ch == '\t') out += "\\t";
        else if (ch == '\n') out += "\\n";
        else if (ch == '\r') out += "\\r";
        else out.push_back(ch);
    }
    return out;
}

}  // namespace

std::string record_id(const GffRecord& rec) {
    return record_id_impl(rec);
}

SummaryRow make_summary_row(const AnnotationIndex& index, const GffRecord& rec) {
    SummaryRow row;
    row.seqid = rec.seqid;
    row.start = rec.start;
    row.end = rec.end;
    row.strand = rec.strand;
    row.type = rec.type;
    row.length = rec.end - rec.start + 1;

    if (rec.id) {
        row.child_count = index.children_of(*rec.id).size();
        add_feature_counts(row, index.descendants_of(*rec.id));
    }

    return row;
}

void print_summary(std::ostream& out, const std::vector<SummaryRow>& rows) {
    out << "seqid\tstart\tend\tstrand\ttype\tlength\tchild_count\ttranscript_count\texon_count\tcds_length\n";

    std::unordered_set<std::string> seqids;
    for (const auto& row : rows) {
        seqids.insert(row.seqid);
    }
    const bool has_all = seqids.size() > 1;

    int64_t all_length = 0;
    size_t all_child_count = 0;
    size_t all_transcript_count = 0;
    size_t all_exon_count = 0;
    int64_t all_cds_length = 0;

    for (const auto& row : rows) {
        out << tsv_escape(row.seqid) << '\t'
            << row.start << '\t'
            << row.end << '\t'
            << row.strand << '\t'
            << tsv_escape(row.type) << '\t'
            << row.length << '\t'
            << row.child_count << '\t'
            << row.transcript_count << '\t'
            << row.exon_count << '\t'
            << row.cds_length << '\n';

        all_length += row.length;
        all_child_count += row.child_count;
        all_transcript_count += row.transcript_count;
        all_exon_count += row.exon_count;
        all_cds_length += row.cds_length;
    }

    if (has_all) {
        out << "all\tNA\tNA\tNA\tNA\t"
            << all_length << '\t'
            << all_child_count << '\t'
            << all_transcript_count << '\t'
            << all_exon_count << '\t'
            << all_cds_length << '\n';
    }
}

}  // namespace gffsub

#include "output.hpp"
#include "gtf_parser.hpp"
#include "record.hpp"
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

namespace gffsub {

void print_gff3(std::ostream& out, const GffData& data) {
    out << "##gff-version 3\n";
    for (const auto& rec : data) {
        if (!rec.kept) continue;
        const std::string& score_str = rec.score_raw.empty() ? std::string{"."} : rec.score_raw;
        // GTF attr_raw is `key "value";` which is invalid GFF3 column 9;
        // rewrite it as tag=value with synthesized ID=/Parent=.
        const std::string& raw_col9 = (rec.src_fmt == InputFormat::GTF) ? gtf_attrs_to_gff3(rec) : rec.attr_raw;
        const std::string& col9 = raw_col9.empty() ? std::string{"."} : raw_col9;
        out << rec.seqid << '\t' << rec.source << '\t' << rec.type << '\t'
            << rec.start << '\t' << rec.end << '\t' << score_str << '\t'
            << rec.strand << '\t' << rec.phase << '\t' << col9 << '\n';
    }
}

static std::string build_gtf_attrs(const std::string& gene_id_val, const std::string& transcript_id_val, bool is_gene) {
    // GTF2.2: gene_id required on every line. transcript_id only on non-gene features.
    // Escape quotes and backslashes in values (GTF attributes are quoted strings).
    auto gtf_escape = [](const std::string& s) {
        std::string out;
        out.reserve(s.size());
        for (const char ch : s) {
            if (ch == '"' || ch == '\\') {
                out.push_back('\\');
            }
            out.push_back(ch);
        }
        return out;
    };
    std::string result;
    result.reserve(64);
    result = "gene_id \"";
    result += gtf_escape(gene_id_val);
    result += "\";";
    if (!is_gene) {
        result += " transcript_id \"";
        result += gtf_escape(transcript_id_val);
        result += "\";";
    }
    return result;
}

void print_gtf(std::ostream& out, const GffData& data, OutputFormat fmt) {
    // GTF header per AGAT spec
    if (fmt == OutputFormat::GTF3) {
        out << "##gtf-version 2.2.1\n";
    } else {
        out << "##gtf-version 2\n";
    }

    // Build mappings from ALL records (not just kept) so parent mRNAs
    // filtered out by subset still resolve gene_id for surviving children.
    std::unordered_map<std::string, std::string> mRNA_to_gene;
    std::unordered_set<std::string> gene_ids;
    for (const auto& rec : data) {
        if ((rec.type == "mRNA" || rec.type == "transcript") && rec.parent_id && rec.id) {
            mRNA_to_gene[*rec.id] = *rec.parent_id;
        }
        if (rec.type == "gene" && rec.id) {
            gene_ids.insert(*rec.id);
        }
    }

    // GTF3 feature types (static to avoid recreation on each call)
    static const std::unordered_set<std::string> gtf3_types = {
        "gene", "transcript", "exon", "CDS", "start_codon", "stop_codon",
        "five_prime_utr", "three_prime_utr", "Selenocysteine", "mRNA"
    };

    for (const auto& rec : data) {
        if (!rec.kept) continue;

        // Filter by GTF3 feature types if using GTF3
        if (fmt == OutputFormat::GTF3) {
            if (gtf3_types.count(rec.type) == 0) continue;
        }

        const std::string& score_str = rec.score_raw.empty() ? std::string{"."} : rec.score_raw;

        std::string gene_id_val;
        std::string transcript_id_val;

        if (rec.type == "gene") {
            gene_id_val = rec.id ? *rec.id : (rec.gene_id ? *rec.gene_id : "");
        } else if (rec.type == "mRNA" || rec.type == "transcript") {
            if (rec.parent_id) {
                gene_id_val = *rec.parent_id;
            } else if (rec.gene_id) {
                gene_id_val = *rec.gene_id;
            }
            transcript_id_val = rec.id ? *rec.id : (rec.transcript_id ? *rec.transcript_id : "");
        } else {
            // Child features
            if (rec.parent_id && mRNA_to_gene.count(*rec.parent_id)) {
                gene_id_val = mRNA_to_gene[*rec.parent_id];
                transcript_id_val = *rec.parent_id;
            } else if (rec.parent_id && gene_ids.count(*rec.parent_id)) {
                // Parent is a gene (e.g. TF_binding_site): no transcript to
                // reference — leave transcript_id empty per the inter/inter_CNS
                // convention instead of fabricating one from the gene ID.
                gene_id_val = *rec.parent_id;
            } else if (rec.parent_id) {
                gene_id_val = *rec.parent_id;
                transcript_id_val = *rec.parent_id;
            } else if (rec.gene_id) {
                gene_id_val = *rec.gene_id;
            }
        }

        // GTF2.2 requires gene_id on every feature line. When it cannot be
        // resolved, emit an empty value (gene_id "";) per the inter/inter_CNS
        // convention rather than dropping the feature silently.
        std::string gtf_type = rec.type;
        if (rec.type == "mRNA") {
            gtf_type = (fmt == OutputFormat::GTF3) ? "transcript" : "mRNA";
        }

        std::string attrs = build_gtf_attrs(gene_id_val, transcript_id_val, rec.type == "gene");

        out << rec.seqid << '\t' << rec.source << '\t' << gtf_type << '\t'
            << rec.start << '\t' << rec.end << '\t' << score_str << '\t'
            << rec.strand << '\t' << rec.phase << '\t' << attrs << '\n';
    }
}

void print_gtf3(std::ostream& out, const GffData& data) {
    print_gtf(out, data, OutputFormat::GTF3);
}

void print_bed(std::ostream& out, const GffData& data) {
    // BED is 0-based half-open, GFF is 1-based inclusive
    // BED start = GFF start - 1, BED end = GFF end
    for (const auto& rec : data) {
        if (!rec.kept) continue;

        std::string name = rec.id ? *rec.id : rec.type;
        const std::string& score_str = (rec.score_raw.empty() || rec.score_raw == ".") ? std::string{"0"} : rec.score_raw;

        out << rec.seqid << '\t'
            << (rec.start - 1) << '\t'
            << rec.end << '\t'
            << name << '\t'
            << score_str << '\t'
            << rec.strand << '\n';
    }
}

}  // namespace gffsub

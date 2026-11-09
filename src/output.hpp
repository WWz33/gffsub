#ifndef GFFSUB_OUTPUT_HPP
#define GFFSUB_OUTPUT_HPP

#include "record.hpp"

#include <cstdint>
#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

namespace gffsub {

// Gene/transcript labels and stop_codon intervals resolved from the whole
// annotation. GTF output needs records the caller may not print: the gene
// behind a kept exon, the transcript label, the stop codon that trims the
// terminal CDS. Build it once from the full GffData, not from a subset.
struct GtfStopCodon {
    int64_t start = 0;
    int64_t end = 0;
    std::string seqid;
    char strand = '.';
};

struct GtfModel {
    std::unordered_map<std::string, std::string> tx_label;   // transcript ID -> emitted transcript_id
    std::unordered_map<std::string, std::string> gene_label; // gene ID -> emitted gene_id
    std::unordered_map<std::string, std::string> tx_parent;  // transcript ID -> gene ID
    std::unordered_map<std::string, std::vector<GtfStopCodon>> stop_codons;
};

GtfModel build_gtf_model(const std::vector<GffRecord>& records);
inline GtfModel build_gtf_model(const GffData& data) {
    return build_gtf_model(data.records);
}

// out_attrs: non-empty whitelist of attribute tags to emit in column 9.
// ID and Parent are always kept (required for a self-consistent feature
// tree); "." passes through untouched.
void print_gff3(std::ostream& out, const GffData& data,
                const std::vector<std::string>& out_attrs = {});
// The model must be built from the same annotation the records came from;
// pass the full GffData even when only a subset of records is printed.
void print_gtf(std::ostream& out, const GffData& data, const GtfModel& model,
               OutputFormat fmt, const std::vector<std::string>& out_attrs = {});
void print_bed(std::ostream& out, const GffData& data);

}  // namespace gffsub

#endif  // GFFSUB_OUTPUT_HPP

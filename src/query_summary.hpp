#ifndef QUERY_SUMMARY_HPP
#define QUERY_SUMMARY_HPP

#include "annotation.hpp"
#include "record.hpp"

#include <cstdint>
#include <ostream>
#include <string>
#include <vector>

namespace gffsub {

struct SummaryRow {
    std::string seqid;
    int64_t start = 0;
    int64_t end = 0;
    char strand = '.';
    std::string type;
    int64_t length = 0;
    size_t child_count = 0;
    size_t transcript_count = 0;
    size_t exon_count = 0;
    int64_t cds_length = 0;
};

// Resolve the ID of a record: ID, gene_id, transcript_id, or empty.
std::string record_id(const GffRecord& rec);

// Build a summary row for a record, counting descendants from the index.
SummaryRow make_summary_row(const AnnotationIndex& index, const GffRecord& rec);

// Print summary rows as TSV. Appends an "all" row summing numeric fields
// when more than one distinct seqid is present.
void print_summary(std::ostream& out, const std::vector<SummaryRow>& rows);

}  // namespace gffsub

#endif  // QUERY_SUMMARY_HPP

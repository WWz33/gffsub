#ifndef QUERY_SUMMARY_HPP
#define QUERY_SUMMARY_HPP

#include "record.hpp"

#include <ostream>
#include <string>
#include <vector>

namespace gffsub {

// Resolve the ID of a record: ID, gene_id, transcript_id, or empty.
std::string record_id(const GffRecord& rec);

// Print an aggregated summary of records as TSV: one row per
// (seqid, type) with feature count and total bases. When more than one
// distinct seqid is present, per-type "all" rows are appended.
void print_summary(std::ostream& out, const std::vector<GffRecord>& records);

}  // namespace gffsub

#endif  // QUERY_SUMMARY_HPP

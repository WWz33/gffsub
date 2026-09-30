#ifndef GFFSUB_OUTPUT_HPP
#define GFFSUB_OUTPUT_HPP

#include "record.hpp"

#include <ostream>
#include <string>
#include <vector>

namespace gffsub {

// out_attrs: non-empty whitelist of attribute tags to emit in column 9.
// ID and Parent are always kept (required for a self-consistent feature
// tree); "." passes through untouched.
void print_gff3(std::ostream& out, const GffData& data,
                const std::vector<std::string>& out_attrs = {});
void print_gtf3(std::ostream& out, const GffData& data,
                const std::vector<std::string>& out_attrs = {});
void print_gtf(std::ostream& out, const GffData& data, OutputFormat fmt,
               const std::vector<std::string>& out_attrs = {});
void print_bed(std::ostream& out, const GffData& data);

}  // namespace gffsub

#endif  // GFFSUB_OUTPUT_HPP

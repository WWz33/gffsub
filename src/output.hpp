#ifndef GFFSUB_OUTPUT_HPP
#define GFFSUB_OUTPUT_HPP

#include "record.hpp"

#include <ostream>

namespace gffsub {

void print_gff3(std::ostream& out, const GffData& data);
void print_gtf3(std::ostream& out, const GffData& data);
void print_gtf(std::ostream& out, const GffData& data, OutputFormat fmt);
void print_bed(std::ostream& out, const GffData& data);

}  // namespace gffsub

#endif  // GFFSUB_OUTPUT_HPP

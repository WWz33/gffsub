#ifndef GFFSUB_CLI_USAGE_HPP
#define GFFSUB_CLI_USAGE_HPP

#include <ostream>

namespace gffsub {

// Explicit -h/--help writes to stdout; a usage error writes to stderr.
void usage(std::ostream& out, const char* prog);
void query_usage(std::ostream& out, const char* prog);
void window_usage(std::ostream& out, const char* prog);

}  // namespace gffsub

#endif  // GFFSUB_CLI_USAGE_HPP

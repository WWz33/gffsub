#ifndef GFFSUB_INPUT_SOURCE_HPP
#define GFFSUB_INPUT_SOURCE_HPP

#include "feature_types.hpp"

#include <iosfwd>
#include <string>
#include <string_view>

namespace gffsub {

class GffData;

// Read the whole input into data.buffer; a regular file is mmap'd read-only
// instead of copied when possible. gzip input (RFC 1952) is inflated
// transparently, including members concatenated by `cat a.gz b.gz`; truncated
// streams fail. Returns 0 on success, -1 on open or decode error.
//
// Non-regular inputs (FIFO, device) cannot be sniffed before reading without
// consuming bytes, so the format is resolved from the loaded content and
// written back through `format`.
int load_input(const std::string& path, GffData& data, InputFormat& format);

// Read an already-open stream (stdin, FIFO, device) into data.buffer and
// inflate it when it starts with the gzip magic. The stream is read once, so
// callers must not reopen it afterwards.
int load_stream(std::istream& in, GffData& data);

// Detect the input format by sniffing the first non-comment feature line.
InputFormat infer_input_format(const std::string& path);
InputFormat infer_format_from_content(std::string_view content);

}  // namespace gffsub

#endif  // GFFSUB_INPUT_SOURCE_HPP

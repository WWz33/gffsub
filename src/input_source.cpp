#include "input_source.hpp"

#include "record.hpp"
#include "string_utils.hpp"

#include <algorithm>
#include <cctype>
#include <fcntl.h>
#include <fstream>
#include <istream>
#include <optional>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

namespace gffsub {

namespace {

// gzip magic bytes (RFC 1952).
bool is_gzip_magic(const unsigned char* magic) {
    return magic[0] == 0x1f && magic[1] == 0x8b;
}

// Inflate a whole gzip file; false on open/decode error. Handles gzip
// members concatenated by `cat a.gz b.gz` (gzread does this natively) and
// detects truncation: a cut-off member latches an error that gzread reports
// as end-of-file, so the error code is checked after the loop.
bool inflate_file(const std::string& path, std::string& out) {
    gzFile gz = gzopen(path.c_str(), "rb");
    if (gz == nullptr) return false;
    std::string result;
    char buf[1u << 20];
    int n = 0;
    while ((n = gzread(gz, buf, sizeof buf)) > 0) {
        result.append(buf, static_cast<size_t>(n));
    }
    int errnum = 0;
    (void)gzerror(gz, &errnum);
    gzclose(gz);
    if (n < 0 || errnum != Z_OK) return false;
    out = std::move(result);
    return true;
}

// Inflate a gzip stream already held in memory (stdin and FIFO paths).
// Feeds the input in chunks (avail_in is a 32-bit count, so a >4GiB stream
// would otherwise be truncated) and continues across concatenated gzip
// members. Truncated input fails: running out of input before Z_STREAM_END
// surfaces as Z_BUF_ERROR with no pending input.
bool inflate_memory(std::string_view in, std::string& out) {
    z_stream strm{};
    if (inflateInit2(&strm, 15 + 16) != Z_OK) return false;  // 16 = gzip wrapper
    std::string result;
    char buf[1u << 20];
    size_t consumed = 0;
    for (;;) {
        if (strm.avail_in == 0) {
            if (consumed >= in.size()) {
                // Out of input before Z_STREAM_END: the member is truncated.
                inflateEnd(&strm);
                return false;
            }
            const size_t chunk = std::min<size_t>(in.size() - consumed, 1u << 30);
            strm.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(in.data() + consumed));
            strm.avail_in = static_cast<uInt>(chunk);
            consumed += chunk;
        }
        strm.next_out = reinterpret_cast<Bytef*>(buf);
        strm.avail_out = sizeof buf;
        const int ret = inflate(&strm, Z_NO_FLUSH);
        if (ret == Z_OK || ret == Z_STREAM_END) {
            result.append(buf, sizeof(buf) - strm.avail_out);
        }
        if (ret == Z_STREAM_END) {
            if (strm.avail_in == 0 && consumed >= in.size()) break;
            // Concatenated members: rewind for the next gzip header.
            if (inflateReset2(&strm, 15 + 16) != Z_OK) {
                inflateEnd(&strm);
                return false;
            }
            continue;
        }
        if (ret != Z_OK) {
            inflateEnd(&strm);
            return false;
        }
    }
    inflateEnd(&strm);
    out = std::move(result);
    return true;
}

// istreambuf_iterator reads byte-by-byte through virtual calls; chunked reads
// are ~10x faster on 100MB+ streams.
std::string read_stream_chunked(std::istream& in) {
    std::string out;
    constexpr size_t kChunk = 1u << 20;
    while (in) {
        const auto old = out.size();
        out.resize(old + kChunk);
        in.read(out.data() + old, kChunk);
        out.resize(old + static_cast<size_t>(in.gcount()));
    }
    return out;
}

// Detect the input format by sniffing the first non-comment feature line.
// Heuristics mirror AGAT's select_gff_format: column-9 character shape wins
// over the filename extension, so a GFF3 file renamed .gtf is still parsed as
// GFF3.
//
//   GTF:  >=9 cols, col9 contains '"' — quoted values are GTF's defining
//         feature and are forbidden raw in GFF3 (spec: encode as %22), so
//         this is checked before '=' to keep GTF values that contain '='
//         (legal in GTF) from being misread as GFF3.
//   GFF3: >=9 cols, col9 contains '='
//   BED:  3..12 cols, no col9, cols[1] and cols[2] are integers
//   default: GFF3
// Decide the format from one feature line. Returns nullopt when the line is
// not recognizable and the caller should keep scanning.
std::optional<InputFormat> sniff_line(std::string_view line) {
    const auto cols = split_line(line, '\t');
    const auto is_int = [](std::string_view s) {
        if (s.empty()) return false;
        size_t pos = 0;
        try {
            std::stoll(std::string{s}, &pos);
        } catch (...) {
            return false;
        }
        return pos == s.size();
    };
    if (cols.size() >= 9) {
        const auto& a = cols[8];
        // Column 9 shape decides. The first fragment separates the two grammars
        // even when a value carries the other syntax's punctuation:
        // GFF3 `ID=g1;Name="A"` (quote inside a value) and GTF
        // `note "a=b"` (`=` inside a value) must not flip the verdict.
        const auto semi = a.find(';');
        const auto first = a.substr(0, semi == std::string_view::npos ? a.size() : semi);
        const auto eq = first.find('=');
        const auto quote = first.find('"');
        if (quote != std::string_view::npos && (eq == std::string_view::npos || quote < eq)) {
            return InputFormat::GTF;  // `tag "value"`
        }
        if (eq != std::string_view::npos) {
            return InputFormat::GFF3;  // `tag=value`
        }
        // Unquoted GTF attributes: `key value; ...` (no quotes, no '=').
        // Shape check keeps BED9 (itemRgb, no semicolon/space) out of this.
        if (a.find(';') != std::string_view::npos) {
            const auto sp = a.find_first_of(" \t");
            if (sp != std::string_view::npos && sp > 0) {
                bool ident = true;
                for (size_t i = 0; i < sp; ++i) {
                    const unsigned char c = static_cast<unsigned char>(a[i]);
                    if (!(std::isalnum(c) || c == '_')) {
                        ident = false;
                        break;
                    }
                }
                if (ident) return InputFormat::GTF;
            }
        }
        // Col9 present but neither GFF3 nor GTF shape. BED9/BED12 lines
        // (itemRgb or a plain col9) also land here: columns 2-3 are integer
        // coordinates (0-based start/end), while a GFF3 source column (col 2)
        // is a label, so two leading integers indicate BED.
        if (is_int(cols[1]) && is_int(cols[2])) {
            return InputFormat::BED;
        }
        return InputFormat::GFF3;
    }
    if (cols.size() >= 3 && cols.size() <= 12) {
        // Tentative BED: require integer start/end.
        if (is_int(cols[1]) && is_int(cols[2])) {
            return InputFormat::BED;
        }
    }
    return std::nullopt;
}

InputFormat sniff_format(const std::string& path) {
    // Only regular files are sniffed in a separate open. Opening a FIFO here
    // would connect to the writer and closing it again would kill the writer
    // (EPIPE) before load_input can read; non-regular inputs are sniffed from
    // the loaded content inside load_input instead.
    struct stat pst{};
    if (::stat(path.c_str(), &pst) != 0 || !S_ISREG(pst.st_mode)) {
        return InputFormat::GFF3;
    }
    std::ifstream f(path, std::ios::binary);
    if (!f.is_open()) return InputFormat::GFF3;
    {
        unsigned char magic[2] = {0, 0};
        f.read(reinterpret_cast<char*>(magic), 2);
        if (f.gcount() == 2 && is_gzip_magic(magic)) {
            f.close();
            gzFile gz = gzopen(path.c_str(), "rb");
            if (gz == nullptr) return InputFormat::GFF3;
            // Read incrementally until one feature line is seen, so a gzipped
            // file with a long header of comments/directives still sniffs
            // from real content (the plain-file path has no byte cap). Line
            // handling mirrors the plain path: strip CR, skip blank and '#'.
            std::string line;
            char c = 0;
            while (gzread(gz, &c, 1) == 1) {
                if (c != '\n') {
                    line.push_back(c);
                    continue;
                }
                if (!line.empty() && line.back() == '\r') line.pop_back();
                if (!line.empty()) {
                    if (line.rfind("##FASTA", 0) == 0) {
                        line.clear();
                        break;
                    }
                    if (line[0] != '#') {
                        if (const auto fmt = sniff_line(line)) {
                            gzclose(gz);
                            return *fmt;
                        }
                    }
                }
                line.clear();
            }
            // Last line without a trailing newline.
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (!line.empty() && line[0] != '#' && line.rfind("##FASTA", 0) != 0) {
                if (const auto fmt = sniff_line(line)) {
                    gzclose(gz);
                    return *fmt;
                }
            }
            gzclose(gz);
            return InputFormat::GFF3;
        }
        f.clear();
        f.seekg(0);
    }
    // On non-seekable input (FIFO/pipe) reading here would consume bytes the
    // later load_input() re-open cannot recover; fall back to GFF3, the
    // project's primary format.
    {
        const auto pos = f.tellg();
        f.seekg(0, std::ios::end);
        if (!f) return InputFormat::GFF3;
        f.seekg(pos);
    }
    std::string line;
    while (std::getline(f, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        if (line.rfind("##FASTA", 0) == 0) break;
        if (line[0] == '#') continue;
        if (const auto fmt = sniff_line(line)) return *fmt;
    }
    return InputFormat::GFF3;
}

}  // namespace

int load_input(const std::string& filename, GffData& data, InputFormat& format) {
    // Probe for gzip only on regular files: reading a FIFO/pipe here would
    // consume bytes the later streaming read cannot recover.
    struct stat pst{};
    const bool regular = (::stat(filename.c_str(), &pst) == 0) && S_ISREG(pst.st_mode);
    std::ifstream probe(filename, std::ios::binary);
    if (!probe.is_open()) return -1;
    unsigned char magic[2] = {0, 0};
    bool gzip_input = false;
    if (regular) {
        probe.read(reinterpret_cast<char*>(magic), 2);
        gzip_input = probe.gcount() == 2 && is_gzip_magic(magic);
        probe.close();
    }

    if (gzip_input) {
        if (!inflate_file(filename, data.buffer)) return -1;
    } else if (regular) {
        int fd = ::open(filename.c_str(), O_RDONLY);
        if (fd < 0) return -1;
        struct stat st{};
        if (::fstat(fd, &st) == 0 && S_ISREG(st.st_mode) && st.st_size > 0) {
            void* map = ::mmap(nullptr, static_cast<size_t>(st.st_size),
                               PROT_READ, MAP_PRIVATE, fd, 0);
            if (map != MAP_FAILED) {
                data.set_mapping(map, static_cast<size_t>(st.st_size));
            }
        }
        ::close(fd);
        if (data.content().empty()) {
            // Failed mapping: read normally.
            std::ifstream file(filename, std::ios::binary);
            if (!file.is_open()) return -1;
            data.buffer = read_stream_chunked(file);
        }
    } else {
        // FIFO/device: a second open() would block waiting for a new writer,
        // so read the already-open stream once. Decompress and sniff from the
        // loaded content, matching the stdin path.
        if (load_stream(probe, data) != 0) return -1;
        probe.close();
        format = infer_format_from_content(data.content());
    }
    return 0;
}

int load_stream(std::istream& in, GffData& data) {
    data.buffer = read_stream_chunked(in);
    // Transparent gzip: `zcat file.gff3.gz | gffsub -` gets the same result
    // as passing the .gz path directly.
    if (data.buffer.size() >= 2) {
        const unsigned char magic[2] = {
            static_cast<unsigned char>(data.buffer[0]),
            static_cast<unsigned char>(data.buffer[1]),
        };
        if (is_gzip_magic(magic)) {
            std::string inflated;
            if (!inflate_memory(data.buffer, inflated)) return -1;
            data.buffer = std::move(inflated);
        }
    }
    return 0;
}

InputFormat infer_input_format(const std::string& path) {
    return sniff_format(path);
}

InputFormat infer_format_from_content(std::string_view content) {
    size_t pos = 0;
    while (pos < content.size()) {
        size_t eol = content.find('\n', pos);
        if (eol == std::string_view::npos) eol = content.size();
        std::string_view line = content.substr(pos, eol - pos);
        pos = eol + 1;
        if (!line.empty() && line.back() == '\r') line.remove_suffix(1);
        if (line.empty()) continue;
        if (line.rfind("##FASTA", 0) == 0) break;
        if (line[0] == '#') continue;
        if (const auto fmt = sniff_line(line)) return *fmt;
    }
    return InputFormat::GFF3;
}

}  // namespace gffsub

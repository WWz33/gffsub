#include "annotation.hpp"
#include "feature_types.hpp"
#include "parser.hpp"
#include "record.hpp"
#include "region.hpp"
#include "gtf_parser.hpp"
#include "string_utils.hpp"
#include <cctype>
#include <algorithm>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

namespace gffsub {

namespace {

// Defined below next to the other stream helpers.
std::string read_stream_chunked(std::istream& in);

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

}  // namespace

static std::optional<std::string> extract_attr_value(std::string_view attrs, std::string_view key) {
    size_t pos = 0;
    while (pos < attrs.size()) {
        // Find the end of the current tag=value pair (next ';').
        size_t pair_end = attrs.find(';', pos);
        if (pair_end == std::string_view::npos) {
            pair_end = attrs.size();
        }
        const auto pair = attrs.substr(pos, pair_end - pos);
        pos = (pair_end < attrs.size()) ? pair_end + 1 : attrs.size();

        // Skip empty fragments (e.g. from ";;").
        if (pair.empty()) {
            continue;
        }

        const size_t eq = pair.find('=');
        if (eq == std::string_view::npos || eq == 0) {
            continue;
        }
        // Trim around the key so `ID=x; Parent=y` resolves: the attribute
        // indexer and --out-attrs trim the same way (attributes.cpp).
        if (trim_view(pair.substr(0, eq)) == key) {
            const auto value = pair.substr(eq + 1);
            if (value.empty()) {
                return std::nullopt;
            }
            // Split the RAW value on ',' before decoding: a literal comma in a
            // single value must stay escaped as %2C (GFF3 spec). Single-value
            // fields take the first part; the index builds full lists from
            // parse_attributes.
            const size_t comma = value.find(',');
            const auto first_part = (comma == std::string_view::npos)
                                        ? value
                                        : value.substr(0, comma);
            if (first_part.empty()) {
                return std::nullopt;
            }
            return url_decode(first_part);
        }
    }
    return std::nullopt;
}

int parse_file(const std::string& filename, GffData& data, InputFormat format) {
    // Read the whole input into data.buffer; record fields are string_views
    // into it. stdin ("-") goes through parse_stdin instead: the format must
    // be sniffed from the buffer after the single read.
    //
    // gzip input (magic 1f 8b) is inflated transparently; mmap-able regular
    // files are mapped read-only to skip the read + zero-fill copy.
    {
        // Probe for gzip only on regular files: reading a FIFO/pipe here
        // would consume bytes the later streaming read cannot recover.
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
            // FIFO/device: a second open() would block waiting for a new
            // writer, so read the already-open stream once. Decompress and
            // sniff from the loaded content, matching the stdin path.
            data.buffer = read_stream_chunked(probe);
            probe.close();
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
            format = infer_format_from_content(data.content());
        }
    }
    // Reserve records capacity: average GFF3 line ~130 bytes. Cap to avoid
    // bad_alloc on comment-heavy files.
    size_t hint = data.content().size() / 130;
    if (hint > 1u << 20) hint = 1u << 20;
    data.reserve(hint);
    return parse_content(data, format);
}

int parse_content(GffData& data, InputFormat format) {
    const std::string_view content{data.content()};
    bool in_fasta = false;
    size_t pos = 0;

    while (pos < content.size()) {
        size_t eol = content.find('\n', pos);
        if (eol == std::string_view::npos) eol = content.size();
        std::string_view line = content.substr(pos, eol - pos);
        pos = eol + 1;

        if (in_fasta) continue;

        // Strip CR from CRLF line endings before any column parsing.
        if (!line.empty() && line.back() == '\r') line.remove_suffix(1);

        if (line.rfind("##FASTA", 0) == 0) { in_fasta = true; continue; }
        if (line.empty()) continue;
        if (line[0] == '#') {
            // Capture ## directives (##gff-version, ##sequence-region, etc.)
            // for re-emission on output. Single-line # comments are skipped.
            if (line.size() > 1 && line[1] == '#') {
                data.directives.emplace_back(line);
            }
            continue;
        }

        GffRecord rec;
        rec.line_idx = static_cast<int>(data.size());
        rec.kept = true;
        rec.src_fmt = format;

        if (format == InputFormat::GFF3 || format == InputFormat::GTF) {
            auto cols = split_line(line, '\t');
            if (cols.size() < 9) continue;

            try {
                rec.seqid = cols[0];
                rec.source = cols[1];
                rec.type = cols[2];
                rec.feat_class = classify_type(rec.type);
                {
                    // GFF3 requires positive integer coordinates, start <= end;
                    // reject trailing garbage. Invalid lines are skipped.
                    size_t pos = 0;
                    rec.start = std::stoll(std::string{cols[3]}, &pos);
                    if (pos != cols[3].size()) throw std::invalid_argument{"start"};
                    pos = 0;
                    rec.end = std::stoll(std::string{cols[4]}, &pos);
                    if (pos != cols[4].size()) throw std::invalid_argument{"end"};
                    if (rec.start < 1 || rec.end < 1 || rec.start > rec.end) {
                        throw std::invalid_argument{"coordinates"};
                    }
                }
                rec.score_raw = cols[5];
                if (cols[5] != ".") {
                    try {
                        size_t spos = 0;
                        const double s = std::stod(std::string{cols[5]}, &spos);
                        if (spos == cols[5].size()) rec.score = s;
                    } catch (...) { rec.score = std::nullopt; }
                }
                rec.strand = cols[6].empty() ? '.' : cols[6][0];
                rec.phase = cols[7].empty() ? '.' : cols[7][0];
                rec.attr_raw = cols[8];

                // The `;`/`=` scan is GFF3-only: it is not quote-aware, so on
                // GTF a quoted value containing "; Parent=..." would be read
                // as a real Parent key. GTF col9 goes through the quote-aware
                // GTF parser instead, which accepts both `key "value";` and
                // bare `key=value` fragments (mixed-format files).
                if (format == InputFormat::GTF) {
                    rec.gene_id = extract_quoted_value(rec.attr_raw, "gene_id");
                    rec.transcript_id = extract_quoted_value(rec.attr_raw, "transcript_id");
                    rec.id = extract_quoted_value(rec.attr_raw, "ID");
                    rec.parent_id = extract_quoted_value(rec.attr_raw, "Parent");
                } else {
                    rec.id = extract_attr_value(rec.attr_raw, "ID");
                    rec.parent_id = extract_attr_value(rec.attr_raw, "Parent");
                    if (rec.attr_raw.find("gene_id=") != std::string::npos) {
                        rec.gene_id = extract_attr_value(rec.attr_raw, "gene_id");
                    }
                    if (rec.attr_raw.find("transcript_id=") != std::string::npos) {
                        rec.transcript_id = extract_attr_value(rec.attr_raw, "transcript_id");
                    }
                    // Multi-parent Parent=tx1,tx2: extract_attr_value already keeps
                    // only the first raw comma part; the index builds the full list
                    // from parse_attributes.
                }

                if (format == InputFormat::GTF) {
                    apply_gtf_attributes(rec);
                    // GTF has no ID=/Parent= attributes; synthesize from
                    // gene_id/transcript_id so the index can build parent/child links.
                    // Guard against empty-string optionals (e.g. Ensembl transcript_id "").
                    if (!rec.id || rec.id->empty()) {
                        if (rec.feat_class == FeatureClass::Gene) {
                            if (rec.gene_id && !rec.gene_id->empty()) {
                                rec.id = rec.gene_id;
                            } else {
                                rec.id = std::nullopt;
                            }
                        } else if (rec.feat_class == FeatureClass::Transcript) {
                            if (rec.transcript_id && !rec.transcript_id->empty()) {
                                rec.id = rec.transcript_id;
                            } else {
                                rec.id = std::nullopt;
                            }
                        } else {
                            // exon/CDS/etc: no ID — linked via parent_id only
                            rec.id = std::nullopt;
                        }
                    }
                    if (!rec.parent_id || rec.parent_id->empty()) {
                        if (rec.feat_class == FeatureClass::Gene) {
                            rec.parent_id = std::nullopt;
                        } else if (rec.feat_class == FeatureClass::Transcript &&
                                   rec.gene_id && !rec.gene_id->empty()) {
                            rec.parent_id = rec.gene_id;
                        } else if (rec.transcript_id && !rec.transcript_id->empty()) {
                            rec.parent_id = rec.transcript_id;
                        } else {
                            rec.parent_id = std::nullopt;
                        }
                    }
                }
            } catch (const std::exception&) {
                continue;
            }
        } else if (format == InputFormat::BED) {
            auto cols = split_line(line, '\t');
            if (cols.size() < 3) continue;

            try {
                rec.seqid = cols[0];
                {
                    size_t pos = 0;
                    const int64_t bed_start = std::stoll(std::string{cols[1]}, &pos);
                    if (pos != cols[1].size()) throw std::invalid_argument{"start"};
                    pos = 0;
                    rec.end = std::stoll(std::string{cols[2]}, &pos);
                    if (pos != cols[2].size()) throw std::invalid_argument{"end"};
                    if (bed_start < 0 || rec.end < bed_start) {
                        throw std::invalid_argument{"coordinates"};
                    }
                    // Zero-length BED interval [s,s) is an insertion site;
                    // map it to the one-base 1-based position s+1 so that
                    // end == start (GFF3 zero-length convention).
                    if (rec.end == bed_start) rec.end = bed_start + 1;
                    rec.start = bed_start + 1;  // BED 0-based half-open -> 1-based inclusive
                }
                rec.source = "gffsub";  // string literal: static storage, no allocation
                rec.type = kRegionType;
                rec.feat_class = FeatureClass::Region;
                rec.score_raw = (cols.size() > 4) ? cols[4] : ".";
                if (cols.size() > 4 && cols[4] != ".") {
                    try {
                        size_t spos = 0;
                        const double s = std::stod(std::string{cols[4]}, &spos);
                        if (spos == cols[4].size()) rec.score = s;
                    } catch (...) { rec.score = std::nullopt; }
                }
                rec.strand = (cols.size() > 5 && !cols[5].empty()) ? cols[5][0] : '.';
                rec.phase = '.';
                rec.id = (cols.size() > 3 && !cols[3].empty()) ? std::optional<std::string>(cols[3]) : std::nullopt;
                rec.parent_id = std::nullopt;
                rec.gene_id = std::nullopt;
                rec.transcript_id = std::nullopt;
                rec.attr_raw = {};  // synthesized as "ID=..." at print_gff3 time
            } catch (const std::exception&) {
                continue;
            }
        }

        data.append(std::move(rec));
    }

    return 0;
}

namespace {

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
        if (a.find('"') != std::string_view::npos) {
            return InputFormat::GTF;
        }
        if (a.find('=') != std::string_view::npos) {
            return InputFormat::GFF3;
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
    // (EPIPE) before parse_file can read; non-regular inputs are sniffed from
    // the loaded content inside parse_file instead.
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
    // later parse_file() re-open cannot recover; fall back to GFF3, the
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

namespace {

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

}  // namespace

int parse_stdin(GffData& data, InputFormat& format_out) {
    data.buffer = read_stream_chunked(std::cin);
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
    const std::string_view content{data.content()};
    format_out = infer_format_from_content(content);
    size_t hint = content.size() / 130;
    if (hint > 1u << 20) hint = 1u << 20;
    data.reserve(hint);
    return parse_content(data, format_out);
}

}  // namespace gffsub

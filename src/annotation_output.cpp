#include "feature_types.hpp"
#include "output.hpp"
#include "gtf_parser.hpp"
#include "parser.hpp"
#include "record.hpp"
#include <unordered_map>
#include <unordered_set>

namespace gffsub {

namespace {

bool is_space(char c) { return c == ' ' || c == '\t'; }

// Keep only the listed tags (plus ID/Parent) in a GFF3 column-9 string.
// Segments are re-emitted verbatim so existing URL encoding is preserved.
std::string project_col9(std::string_view col9, const std::vector<std::string>& keep_tags) {
    if (keep_tags.empty() || col9.empty() || col9 == ".") {
        return std::string{col9};
    }
    std::unordered_set<std::string> keep{keep_tags.begin(), keep_tags.end()};
    keep.insert("ID");
    keep.insert("Parent");
    std::string out;
    size_t pos = 0;
    while (pos <= col9.size()) {
        size_t end = col9.find(';', pos);
        if (end == std::string_view::npos) end = col9.size();
        std::string_view pair = col9.substr(pos, end - pos);
        const auto eq = pair.find('=');
        if (eq != std::string_view::npos && eq > 0) {
            size_t kstart = 0;
            while (kstart < eq && is_space(pair[kstart])) ++kstart;
            size_t kend = eq;
            while (kend > kstart && is_space(pair[kend - 1])) --kend;
            if (keep.count(std::string{pair.substr(kstart, kend - kstart)}) > 0) {
                if (!out.empty()) out += ';';
                out += pair;
            }
        }
        if (end == col9.size()) break;
        pos = end + 1;
    }
    if (out.empty()) out = ".";
    return out;
}

// A sub-subset must not declare landmarks it no longer contains: drop
// ##sequence-region directives whose seqid is absent from the kept records.
bool directive_survives(const std::string& d,
                        const std::unordered_set<std::string>& kept_seqids) {
    constexpr std::string_view kPrefix = "##sequence-region";
    if (d.rfind(kPrefix, 0) != 0) return true;
    // Require the directive name to end here (whitespace or EOL) so an
    // application-specific ##sequence-regionfoo is not misread.
    if (d.size() > kPrefix.size() &&
        d[kPrefix.size()] != ' ' && d[kPrefix.size()] != '\t') {
        return true;
    }
    size_t p = d.find_first_not_of(" \t", kPrefix.size());
    if (p == std::string::npos) return true;
    size_t q = d.find_first_of(" \t", p);
    const std::string seqid = d.substr(p, q == std::string::npos ? q : q - p);
    return kept_seqids.count(seqid) > 0;
}

std::unordered_set<std::string> collect_kept_seqids(const GffData& data) {
    std::unordered_set<std::string> seqids;
    for (const auto& rec : data) {
        if (rec.kept) seqids.insert(std::string{rec.seqid});
    }
    return seqids;
}

}  // namespace

void print_gff3(std::ostream& out, const GffData& data,
                const std::vector<std::string>& out_attrs) {
    // ##gff-version must be the topmost line (GFF3 spec). Take it from the
    // input when present, otherwise synthesize; ##gtf-version is a GTF
    // directive and must not leak into GFF3 output.
    std::string_view version_line = "##gff-version 3";
    for (const auto& d : data.directives) {
        if (d.rfind("##gff-version", 0) == 0) {
            version_line = d;
            break;
        }
    }
    out << version_line << '\n';
    const auto kept_seqids = collect_kept_seqids(data);
    for (const auto& d : data.directives) {
        if (d.rfind("##gff-version", 0) == 0) continue;
        if (d.rfind("##gtf-version", 0) == 0) continue;
        if (!directive_survives(d, kept_seqids)) continue;
        out << d << '\n';
    }
    for (const auto& rec : data) {
        if (!rec.kept) continue;
        const std::string_view score_str = rec.score_raw.empty() ? std::string_view{"."} : rec.score_raw;
        std::string col9_storage;  // only used when col9 must be synthesized
        std::string_view col9;
        if (rec.src_fmt == InputFormat::GTF) {
            // GTF attr_raw is `key "value";` which is invalid GFF3 column 9;
            // rewrite it as tag=value with synthesized ID=/Parent=.
            col9_storage = gtf_attrs_to_gff3(rec);
            col9 = col9_storage;
        } else if (rec.attr_raw.empty() && rec.id) {
            // BED-sourced record: attr_raw was left empty at parse time.
            col9_storage = "ID=" + *rec.id;
            col9 = col9_storage;
        } else {
            col9 = rec.attr_raw;
        }
        if (col9.empty()) col9 = ".";
        std::string col9_projected;
        if (!out_attrs.empty()) {
            col9_projected = project_col9(col9, out_attrs);
            col9 = col9_projected;
        }
        out << rec.seqid << '\t' << rec.source << '\t' << rec.type << '\t'
            << rec.start << '\t' << rec.end << '\t' << score_str << '\t'
            << rec.strand << '\t' << rec.phase << '\t' << col9 << '\n';
    }
}

static std::string build_gtf_attrs(const std::string& gene_id_val, const std::string& transcript_id_val,
                                   bool is_gene, const GffRecord& rec,
                                   const std::vector<std::string>& out_attrs) {
    // GTF2.2: gene_id required on every line, transcript_id on non-gene
    // features; other attributes may follow (spec: "Any other attributes or
    // comments must appear after these two"). Preserve them like AGAT does
    // ("attribute conserved: All"): URL-decode GFF3 values, quote, escape.
    auto gtf_escape = [](const std::string& s) {
        std::string out;
        out.reserve(s.size());
        for (const char ch : s) {
            if (ch == '"') { out += "\\\""; continue; }
            if (ch == '\\') { out += "\\\\"; continue; }
            if (ch == '\n') { out += "\\n"; continue; }
            if (ch == '\t') { out += "\\t"; continue; }
            if (ch == '\r') { out += "\\r"; continue; }
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
    // Pass through the remaining attributes. Skip keys already emitted above.
    // GTF-source col9 is `key "value";`, which parse_attributes (tag=value)
    // cannot read — it would drop or corrupt every attribute. Use the GTF
    // pair parser for GTF-source records.
    std::unordered_set<std::string> keep;
    if (!out_attrs.empty()) keep.insert(out_attrs.begin(), out_attrs.end());
    const auto wanted = [&](std::string_view key) {
        if (key == "gene_id" || key == "transcript_id") return false;
        if (keep.empty()) return true;
        return keep.count(std::string{key}) > 0;
    };
    if (rec.src_fmt == InputFormat::GTF) {
        for (const auto& [key, value] : parse_gtf_attributes(rec.attr_raw)) {
            if (!wanted(key)) continue;
            result += ' ';
            result += key;
            result += " \"";
            result += gtf_escape(value);
            result += "\";";
        }
    } else {
        for (const auto& [key, values] : parse_attributes(rec.attr_raw)) {
            if (!wanted(key) || key == "Parent" || key == "ID") continue;
            for (const auto& v : values) {
                result += ' ';
                result += key;
                result += " \"";
                result += gtf_escape(v);
                result += "\";";
            }
        }
    }
    return result;
}

void print_gtf(std::ostream& out, const GffData& data, OutputFormat fmt,
               const std::vector<std::string>& out_attrs) {
    // GTF header per AGAT spec
    if (fmt == OutputFormat::GTF3) {
        out << "##gtf-version 2.2.1\n";
    } else {
        out << "##gtf-version 2\n";
    }
    // Preserve ## directives from input (##sequence-region, ##species, etc.)
    // but skip ##gff-version / ##gtf-version — the header line above already
    // declares the output format version.
    const auto kept_seqids = collect_kept_seqids(data);
    for (const auto& d : data.directives) {
        if (d.rfind("##gff-version", 0) == 0) continue;
        if (d.rfind("##gtf-version", 0) == 0) continue;
        if (!directive_survives(d, kept_seqids)) continue;
        out << d << '\n';
    }

    // Build mappings from ALL records (not just kept) so parent mRNAs
    // filtered out by subset still resolve gene_id for surviving children.
    std::unordered_map<std::string, std::string> mRNA_to_gene;
    std::unordered_set<std::string> gene_ids;
    for (const auto& rec : data) {
        if ((rec.feat_class == FeatureClass::Transcript) && rec.parent_id && rec.id) {
            mRNA_to_gene[*rec.id] = *rec.parent_id;
        }
        if (rec.feat_class == FeatureClass::Gene && rec.id) {
            gene_ids.insert(*rec.id);
        }
    }

    for (const auto& rec : data) {
        if (!rec.kept) continue;

        // Normalize the type label first, then filter — so that aliases
        // like "5'-utr" or "five_prime_UTR" pass the GTF3 whitelist after
        // being normalized to "five_prime_utr".
        std::string gtf_type{gtf_type_label(rec.type, fmt)};

        if (fmt == OutputFormat::GTF3) {
            if (!gtf_type_emittable(gtf_type, fmt)) continue;
        }

        const std::string score_str = rec.score_raw.empty() ? std::string{"."} : std::string{rec.score_raw};

        std::string gene_id_val;
        std::string transcript_id_val;

        if (rec.feat_class == FeatureClass::Gene) {
            if (rec.id) {
                gene_id_val = *rec.id;
            } else if (rec.gene_id) {
                gene_id_val = *rec.gene_id;
            }
        } else if (rec.feat_class == FeatureClass::Transcript) {
            if (rec.parent_id) {
                gene_id_val = *rec.parent_id;
            } else if (rec.gene_id) {
                gene_id_val = *rec.gene_id;
            }
            if (rec.id) {
                transcript_id_val = *rec.id;
            } else if (rec.transcript_id) {
                transcript_id_val = *rec.transcript_id;
            }
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
                // Parent unresolvable to a gene/transcript record in this file
                // (flat GTF with only exon/CDS lines): trust the source's
                // gene_id attribute over the parent ID.
                gene_id_val = rec.gene_id ? *rec.gene_id : *rec.parent_id;
                transcript_id_val = *rec.parent_id;
            } else if (rec.gene_id) {
                gene_id_val = *rec.gene_id;
            }
        }

        // GTF2.2 requires gene_id on every feature line. When it cannot be
        // resolved, emit an empty value (gene_id "";) per the inter/inter_CNS
        // convention rather than dropping the feature silently.
        std::string attrs = build_gtf_attrs(gene_id_val, transcript_id_val,
                                            rec.feat_class == FeatureClass::Gene, rec,
                                            out_attrs);

        out << rec.seqid << '\t' << rec.source << '\t' << gtf_type << '\t'
            << rec.start << '\t' << rec.end << '\t' << score_str << '\t'
            << rec.strand << '\t' << rec.phase << '\t' << attrs << '\n';
    }
}

void print_gtf3(std::ostream& out, const GffData& data,
                const std::vector<std::string>& out_attrs) {
    print_gtf(out, data, OutputFormat::GTF3, out_attrs);
}

void print_bed(std::ostream& out, const GffData& data) {
    // BED is 0-based half-open, GFF is 1-based inclusive
    // BED start = GFF start - 1, BED end = GFF end
    for (const auto& rec : data) {
        if (!rec.kept) continue;

        std::string name = rec.id ? *rec.id : std::string{rec.type};
        const std::string score_str = (rec.score_raw.empty() || rec.score_raw == ".") ? std::string{"0"} : std::string{rec.score_raw};

        out << rec.seqid << '\t'
            << (rec.start - 1) << '\t'
            << rec.end << '\t'
            << name << '\t'
            << score_str << '\t'
            << rec.strand << '\n';
    }
}

}  // namespace gffsub

#include "feature_types.hpp"
#include "output.hpp"
#include "gtf_parser.hpp"
#include "parser.hpp"
#include "record.hpp"
#include "string_utils.hpp"
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

namespace gffsub {

namespace {

// Loud enough to notice on a broken file, quiet enough for a big one.
void warn_duplicate_transcript_id(std::string_view id) {
    constexpr int kLimit = 5;
    static int warned = 0;
    if (warned > kLimit) return;
    if (warned == kLimit) {
        std::cerr << "Warning: further duplicate-transcript-ID warnings suppressed\n";
    } else {
        std::cerr << "Warning: transcript ID '" << id
                  << "' is defined more than once with different parents; GTF"
                     " gene_id assignment may be wrong\n";
    }
    ++warned;
}

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
            // Same trim as attributes.cpp so a key matches whichever parser
            // produced it (space/tab/CR around the key).
            const std::string key{trim_view(pair.substr(0, eq))};
            if (keep.count(key) > 0) {
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

// Format header plus the input's ## directives, minus the version lines the
// header already declares and minus landmarks the kept rows no longer cover.
void emit_prologue(std::ostream& out, const GffData& data, std::string_view version_line) {
    if (!version_line.empty()) {
        out << version_line << '\n';
    }
    const auto kept_seqids = collect_kept_seqids(data);
    for (const auto& d : data.directives) {
        if (d.rfind("##gff-version", 0) == 0) continue;
        if (d.rfind("##gtf-version", 0) == 0) continue;
        if (!directive_survives(d, kept_seqids)) continue;
        out << d << '\n';
    }
}

}  // namespace

void print_gff3(std::ostream& out, const GffData& data,
                const std::vector<std::string>& out_attrs) {
    // ##gff-version must be the topmost line and names version 3 (GFF3 spec).
    // An input line claiming an older version is not reproduced: the rows
    // below are GFF3 whatever the source said.
    emit_prologue(out, data, "##gff-version 3");
    for (const auto& rec : data) {
        if (!rec.kept) continue;
        const std::string_view score_str = rec.score_raw.empty() ? std::string_view{"."} : rec.score_raw;
        // GTF spellings such as 5UTR are not SO terms; GFF3 column 3 must be
        // a Sequence Ontology name.
        const std::string_view type_label = gff3_type_label(rec.type);
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
        out << rec.seqid << '\t' << rec.source << '\t' << type_label << '\t'
            << rec.start << '\t' << rec.end << '\t' << score_str << '\t'
            << rec.strand << '\t' << rec.phase << '\t' << col9 << '\n';
    }
}

static std::string build_gtf_attrs(const std::string& gene_id_val, const std::string& transcript_id_val,
                                   bool is_gene, const GffRecord& rec,
                                   const std::vector<std::string>& out_attrs,
                                   bool drop_parent = false,
                                   std::string_view id_override = {}) {
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
    // One decoder handles both source formats.
    //
    // drop_parent drops Parent from a child row that expansion split under
    // several parents. Each such row names its own transcript, and repeating
    // the original parent list makes a re-read expand the row again (two
    // parents double the child row count per pass). A single-parent row keeps
    // Parent: it names the same parent transcript_id already carries.
    std::unordered_set<std::string> keep;
    if (!out_attrs.empty()) keep.insert(out_attrs.begin(), out_attrs.end());
    const auto wanted = [&](std::string_view key) {
        if (key == "gene_id" || key == "transcript_id") return false;
        // ID and Parent carry the feature tree, so --out-attrs never drops
        // them. Both output formats follow the same rule.
        if (key == "ID") return true;
        if (key == "Parent") return !drop_parent;
        if (keep.empty()) return true;
        return keep.count(std::string{key}) > 0;
    };
    for (const auto& [key, value] : decode_attributes(rec.attr_raw, rec.src_fmt)) {
        if (!wanted(key)) continue;
        const std::string_view out_value =
            (key == "ID" && !id_override.empty()) ? id_override : std::string_view{value};
        result += ' ';
        result += key;
        result += " \"";
        result += gtf_escape(std::string{out_value});
        result += "\";";
    }
    return result;
}

GtfModel build_gtf_model(const std::vector<GffRecord>& records) {
    GtfModel model;
    // Build mappings from ALL records (not just kept) so parent mRNAs
    // filtered out by subset still resolve gene_id for surviving children.
    // A transcript ID defined twice with different parents is invalid input
    // (spec: IDs are unique); the last definition wins for every child, so
    // report it rather than hand out a silently wrong gene_id.
    //
    // The emitted identifier prefers the record's gene_id / transcript_id
    // attribute over its ID/Parent. Ensembl GFF3 labels genes `gene:gene-X`
    // and carries the canonical `gene_id=gene-X` alongside, and GTF consumers
    // expect the unprefixed form.
    // stop_codon intervals per transcript, used to trim the terminal CDS
    // (GTF2.2: "the stop codon is not included in the CDS for the terminal
    // exon"). Keyed by the Parent value as written on the stop_codon line.
    // seqid and strand travel with the interval: a stop codon on another
    // chromosome or the opposite strand does not trim this CDS.
    for (const auto& rec : records) {
        if (rec.feat_class == FeatureClass::Transcript && rec.id) {
            if (rec.parent_id) {
                const auto existing = model.tx_parent.find(*rec.id);
                if (existing != model.tx_parent.end() && existing->second != *rec.parent_id) {
                    warn_duplicate_transcript_id(*rec.id);
                }
                model.tx_parent[*rec.id] = *rec.parent_id;
            }
            model.tx_label[*rec.id] = rec.transcript_id ? *rec.transcript_id : *rec.id;
        }
        if (rec.feat_class == FeatureClass::Gene && rec.id) {
            model.gene_label[*rec.id] = rec.gene_id ? *rec.gene_id : *rec.id;
        }
        if (rec.feat_class == FeatureClass::Unknown &&
            to_lower(rec.type) == "stop_codon") {
            // Key by every Parent: a shared stop codon trims each parent's
            // terminal CDS, the same way a shared exon expands to every parent.
            for (const auto& pid : parent_ids(rec)) {
                model.stop_codons[pid].push_back(GtfStopCodon{
                    rec.start, rec.end, std::string{rec.seqid}, rec.strand});
            }
        }
    }
    return model;
}

void print_gtf(std::ostream& out, const GffData& data, const GtfModel& model,
               OutputFormat fmt, const std::vector<std::string>& out_attrs) {
    // Header names the GTF flavor whose vocabulary the rows below use.
    std::string_view version_line;
    switch (fmt) {
        case OutputFormat::GTF3:  version_line = "##gtf-version 3"; break;
        case OutputFormat::GTF25: version_line = "##gtf-version 2.5"; break;
        case OutputFormat::GTF22: version_line = "##gtf-version 2.2"; break;
        default: break;
    }
    emit_prologue(out, data, version_line);

    // Resolve the emitted gene_id / transcript_id for one parent link. The
    // parent chain decides, so the IDs match the parent transcript/gene rows.
    // A child that also carries its own transcript_id keeps that value only
    // when the row is not an expanded multi-parent row: on expansion the
    // current parent names the row, otherwise every row would name the same
    // transcript.
    const auto resolve_child = [&](const std::string& pid, const GffRecord& rec,
                                   bool expanded, std::string& gene_id_val, std::string& transcript_id_val) {
        if (rec.gene_id) gene_id_val = *rec.gene_id;
        if (rec.transcript_id && !expanded) transcript_id_val = *rec.transcript_id;
        const auto tx = model.tx_label.find(pid);
        if (tx != model.tx_label.end()) {
            if (transcript_id_val.empty()) transcript_id_val = tx->second;
            if (gene_id_val.empty()) {
                const auto tp = model.tx_parent.find(pid);
                if (tp != model.tx_parent.end()) {
                    const auto g = model.gene_label.find(tp->second);
                    gene_id_val = (g != model.gene_label.end()) ? g->second : tp->second;
                }
            }
        } else if (model.gene_label.count(pid) > 0) {
            // Parent is a gene (e.g. TF_binding_site): no transcript to
            // reference — leave transcript_id empty per the inter/inter_CNS
            // convention instead of inventing one.
            if (gene_id_val.empty()) {
                const auto g = model.gene_label.find(pid);
                gene_id_val = (g != model.gene_label.end()) ? g->second : pid;
            }
        } else {
            // Parent unresolvable to a gene/transcript record in this file
            // (flat GTF with only exon/CDS lines): trust the parent ID as the
            // transcript.
            if (transcript_id_val.empty()) transcript_id_val = pid;
            if (gene_id_val.empty()) gene_id_val = pid;
        }
        // A transcript with no gene row leaves gene_id unresolved. GTF2.2
        // requires the attribute, so use the transcript identity rather than
        // emitting `gene_id "";`.
        if (gene_id_val.empty()) {
            gene_id_val = transcript_id_val.empty() ? pid : transcript_id_val;
        }
    };

    for (const auto& rec : data) {
        if (!rec.kept) continue;

        // Normalize the type label first, then filter — so that aliases
        // like "5'-utr" or "five_prime_UTR" pass the GTF3 whitelist after
        // being normalized to "five_prime_utr".
        std::string gtf_type{gtf_type_label(rec.type, fmt)};

        if (!gtf_type_emittable(gtf_type, fmt)) continue;

        const std::string score_str = rec.score_raw.empty() ? std::string{"."} : std::string{rec.score_raw};

        // One output row. parent_id selects the stop_codon used to trim the
        // terminal CDS (empty for gene and transcript rows). drop_parent is set
        // when expansion split one source row under several parents;
        // id_override gives each expanded row a unique ID.
        const auto emit_row = [&](const std::string& gene_id_val, const std::string& transcript_id_val,
                                  const std::string& parent_id, bool drop_parent,
                                  std::string_view id_override = {}) {
            // GTF2.2: gene_id required on every feature line. When it cannot
            // be resolved, emit an empty value (gene_id "";) per the
            // inter/inter_CNS convention rather than dropping the feature.
            std::string attrs = build_gtf_attrs(gene_id_val, transcript_id_val,
                                                rec.feat_class == FeatureClass::Gene, rec,
                                                out_attrs, drop_parent, id_override);

            // GTF excludes the stop codon from the terminal CDS. Trim only
            // when the stop_codon sits hard against the CDS 3' end of the same
            // transcript, so partially overlapping or shifted annotations pass
            // through untouched.
            int64_t start = rec.start;
            int64_t end = rec.end;
            if (rec.feat_class == FeatureClass::CDS && !parent_id.empty()) {
                const auto it = model.stop_codons.find(parent_id);
                if (it != model.stop_codons.end()) {
                    for (const auto& stop : it->second) {
                        // A stop codon on another seqid or strand is a
                        // different feature with a colliding Parent value.
                        if (stop.seqid != rec.seqid || stop.strand != rec.strand) {
                            continue;
                        }
                        if (rec.strand == '+' && stop.end == rec.end && stop.start > rec.start) {
                            end = stop.start - 1;
                            break;
                        }
                        if (rec.strand == '-' && stop.start == rec.start && stop.end < rec.end) {
                            start = stop.end + 1;
                            break;
                        }
                    }
                }
            }

            out << rec.seqid << '\t' << rec.source << '\t' << gtf_type << '\t'
                << start << '\t' << end << '\t' << score_str << '\t'
                << rec.strand << '\t' << rec.phase << '\t' << attrs << '\n';
        };

        if (rec.feat_class == FeatureClass::Gene) {
            std::string gene_id_val;
            if (rec.id) {
                const auto it = model.gene_label.find(*rec.id);
                gene_id_val = (it != model.gene_label.end()) ? it->second : *rec.id;
            } else if (rec.gene_id) {
                gene_id_val = *rec.gene_id;
            }
            emit_row(gene_id_val, std::string{}, std::string{}, true);
            continue;
        }

        if (rec.feat_class == FeatureClass::Transcript) {
            std::string gene_id_val;
            std::string transcript_id_val;
            if (rec.id) {
                const auto it = model.tx_label.find(*rec.id);
                transcript_id_val = (it != model.tx_label.end()) ? it->second : *rec.id;
            } else if (rec.transcript_id) {
                transcript_id_val = *rec.transcript_id;
            }
            if (rec.parent_id) {
                const auto it = model.gene_label.find(*rec.parent_id);
                gene_id_val = (it != model.gene_label.end()) ? it->second : *rec.parent_id;
            } else if (rec.gene_id) {
                gene_id_val = *rec.gene_id;
            }
            // GTF2.2 requires gene_id on every line. A transcript with no
            // gene row and no gene_id falls back to its own transcript_id
            // instead of emitting the invalid `gene_id "";`.
            if (gene_id_val.empty()) gene_id_val = transcript_id_val;
            emit_row(gene_id_val, transcript_id_val, std::string{}, false);
            continue;
        }

        // GTF has no multi-parent feature: a shared exon becomes one row per
        // parent transcript. parent_ids handles Parent=a,b and a quoted GTF
        // Parent "a,b" the same way.
        const auto parents = parent_ids(rec);

        if (parents.empty()) {
            std::string gene_id_val;
            std::string transcript_id_val;
            if (rec.gene_id) gene_id_val = *rec.gene_id;
            if (rec.transcript_id) transcript_id_val = *rec.transcript_id;
            // A top-level feature with no parent still needs both mandatory
            // GTF attributes; the record's own ID is the only identity left.
            if (transcript_id_val.empty() && rec.id) transcript_id_val = *rec.id;
            if (gene_id_val.empty()) gene_id_val = transcript_id_val;
            emit_row(gene_id_val, transcript_id_val, std::string{}, false);
            continue;
        }

        const bool expanded = parents.size() > 1;
        for (const auto& pid : parents) {
            std::string gene_id_val;
            std::string transcript_id_val;
            resolve_child(pid, rec, expanded, gene_id_val, transcript_id_val);
            // One row per parent is a distinct feature: AGAT gives each
            // expanded row a unique ID (agat_how_does_it_work.md:26).
            std::string expanded_id;
            if (expanded && rec.id) {
                expanded_id = *rec.id + "." + pid;
            }
            emit_row(gene_id_val, transcript_id_val, pid, expanded, expanded_id);
        }
    }
}

void print_bed(std::ostream& out, const GffData& data) {
    // BED is 0-based half-open, GFF is 1-based inclusive
    // BED start = GFF start - 1, BED end = GFF end
    for (const auto& rec : data) {
        if (!rec.kept) continue;

        // A parsed record always has start >= 1; a library-built record can
        // start at 0, where start - 1 would print -1. Clamp to the first base
        // instead of writing an invalid coordinate.
        const int64_t bed_start = rec.start > 1 ? rec.start - 1 : 0;

        std::string name = rec.id ? *rec.id : std::string{rec.type};
        // A GFF3 ID may decode %09/%0A into a real tab or newline. BED has
        // no escape syntax, so a control character would add a column or a
        // line. Replace them with a space.
        for (char& ch : name) {
            if (ch == '\t' || ch == '\n' || ch == '\r') ch = ' ';
        }
        // BED score is an integer from 0 to 1000. Out-of-range, fractional,
        // unparseable, and non-finite scores have no BED representation, so
        // they fall back to 0.
        std::string score_str = "0";
        if (rec.score && std::isfinite(*rec.score)) {
            const double rounded = *rec.score < 0 ? 0 : (*rec.score > 1000 ? 1000 : *rec.score);
            score_str = std::to_string(static_cast<int64_t>(rounded));
        }

        out << rec.seqid << '\t'
            << bed_start << '\t'
            << rec.end << '\t'
            << name << '\t'
            << score_str << '\t'
            << rec.strand << '\n';
    }
}

}  // namespace gffsub

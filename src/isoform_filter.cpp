#include "filter.hpp"
#include "parser.hpp"
#include "record.hpp"
#include <algorithm>
#include <future>
#include <map>
#include <unordered_map>

namespace gffsub {

void filter_longest_isoform(GffData& data, std::string_view longest_type_sv, size_t num_threads) {
    std::string isoform_type{longest_type_sv};
    if (isoform_type.empty()) {
        // Auto-detect: pick the most frequent transcript-class type among
        // kept records (mRNA, transcript, ncRNA, tRNA, ...). Ties prefer
        // mRNA, then transcript, then lexicographic order for determinism.
        std::map<std::string, size_t> counts;
        for (const auto& rec : data.records) {
            if (!rec.kept) continue;
            if (rec.feat_class == FeatureClass::Transcript) ++counts[rec.type];
        }
        if (!counts.empty()) {
            const auto best = std::max_element(
                counts.begin(), counts.end(),
                [](const auto& a, const auto& b) {
                    if (a.second != b.second) return a.second < b.second;
                    const auto rank = [](const std::string& t) {
                        if (t == "mRNA") return 0;
                        if (t == "transcript") return 1;
                        return 2;
                    };
                    const int ra = rank(a.first);
                    const int rb = rank(b.first);
                    if (ra != rb) return ra > rb;
                    return a.first > b.first;
                });
            isoform_type = best->first;
        }
    }

    // Build gene -> [isoform indices] index once.
    // Use parse_attributes to get ALL Parent values so a multi-parent isoform
    // (Parent=geneA,geneB) is registered under every gene it belongs to.
    std::unordered_map<std::string, std::vector<int>> gene_to_isoforms;
    for (int i = 0; i < static_cast<int>(data.records.size()); ++i) {
        const auto& rec = data.records[i];
        if (!rec.kept) continue;
        if (rec.type != isoform_type) continue;
        const auto attrs = parse_attributes(rec.attr_raw);
        const auto parent_it = attrs.find("Parent");
        if (parent_it != attrs.end()) {
            for (const auto& parent_id : parent_it->second) {
                auto& vec = gene_to_isoforms[parent_id];
                if (std::find(vec.begin(), vec.end(), i) == vec.end()) {
                    vec.push_back(i);
                }
            }
        } else if (rec.parent_id) {
            gene_to_isoforms[*rec.parent_id].push_back(i);
        }
    }

    // Build isoform -> [child indices] index once, using the full Parent list
    // (parse_attributes splits Parent=a,b) so shared children link to all parents.
    std::unordered_map<std::string, std::vector<int>> isoform_to_children;
    for (int i = 0; i < static_cast<int>(data.records.size()); ++i) {
        const auto& rec = data.records[i];
        if (!rec.kept) continue;
        const auto attrs = parse_attributes(rec.attr_raw);
        const auto parent_it = attrs.find("Parent");
        if (parent_it != attrs.end()) {
            for (const auto& parent_id : parent_it->second) {
                auto& vec = isoform_to_children[parent_id];
                if (std::find(vec.begin(), vec.end(), i) == vec.end()) {
                    vec.push_back(i);
                }
            }
        } else if (rec.parent_id) {
            isoform_to_children[*rec.parent_id].push_back(i);
        }
    }

    // Group genes by chromosome
    std::unordered_map<std::string, std::vector<int>> chrom_to_gene_idx;
    for (int i = 0; i < static_cast<int>(data.records.size()); ++i) {
        const auto& rec = data.records[i];
        if (!rec.kept) continue;
        if (rec.feat_class == FeatureClass::Gene && rec.id) {
            chrom_to_gene_idx[rec.seqid].push_back(i);
        }
    }

    auto process_chromosome = [&](const std::string& chrom, const std::vector<int>& gene_indices) {
        (void)chrom; // reserved for future debugging
        for (int gene_idx : gene_indices) {
            const auto& gene = data.records[gene_idx];
            if (!gene.id) continue;

            // Find isoforms for this gene using index
            auto isoform_it = gene_to_isoforms.find(*gene.id);
            if (isoform_it == gene_to_isoforms.end()) {
                // gene with no isoform children: keep only if it has any
                // non-isoform children (e.g. TF_binding_site). A gene with
                // no children at all is dropped.
                auto child_it = isoform_to_children.find(*gene.id);
                if (child_it == isoform_to_children.end()) {
                    data.records[gene_idx].kept = false;
                }
                continue;
            }
            if (isoform_it->second.size() <= 1) continue;

            const auto& isoform_indices = isoform_it->second;

            // Per-gene check: does ANY isoform have CDS?
            bool gene_has_cds = false;
            for (int iso_idx : isoform_indices) {
                const auto& iso = data.records[iso_idx];
                if (!iso.id) continue;
                auto child_it = isoform_to_children.find(*iso.id);
                if (child_it != isoform_to_children.end()) {
                    for (int child_idx : child_it->second) {
                        if (data.records[child_idx].feat_class == FeatureClass::CDS) {
                            gene_has_cds = true;
                            break;
                        }
                    }
                }
                if (gene_has_cds) break;
            }

            // Find longest isoform
            int longest_idx = -1;
            int64_t max_len = -1;

            for (int iso_idx : isoform_indices) {
                const auto& iso = data.records[iso_idx];
                if (!iso.id) continue;

                auto child_it = isoform_to_children.find(*iso.id);
                if (child_it == isoform_to_children.end()) continue;

                int64_t len = 0;
                bool found = false;

                if (gene_has_cds) {
                    // Distinct CDS IDs under one transcript are alternative
                    // protein variants (spec EDEN example: cds00003/cds00004);
                    // score by the LONGEST variant, do not sum across variants.
                    // Lines sharing an ID are one discontinuous CDS (summed).
                    std::unordered_map<std::string, int64_t> variant_len;
                    for (int child_idx : child_it->second) {
                        const auto& child = data.records[child_idx];
                        if (child.feat_class == FeatureClass::CDS) {
                            const std::string key = child.id ? *child.id : std::string{};
                            variant_len[key] += child.end - child.start + 1;
                        }
                    }
                    if (variant_len.empty()) continue; // isoform without CDS is skipped
                    for (const auto& [key, vlen] : variant_len) {
                        (void)key;
                        if (vlen > len) len = vlen;
                    }
                } else {
                    for (int child_idx : child_it->second) {
                        const auto& child = data.records[child_idx];
                        if (child.feat_class == FeatureClass::Exon) {
                            len += child.end - child.start + 1;
                            found = true;
                        }
                    }
                    // No exon children: cannot be "longest" by exon span,
                    // mirroring gffread's covlen (sum of exon lengths).
                    if (!found) continue;
                }

                if (len > max_len) {
                    max_len = len;
                    longest_idx = iso_idx;
                }
            }

            // Mark longest as kept, others as not kept
            if (longest_idx >= 0) {
                for (int iso_idx : isoform_indices) {
                    data.records[iso_idx].kept = (iso_idx == longest_idx);
                }
                // Drop children of all isoforms, then re-keep children of the
                // longest. Two-pass so a child shared between the longest and a
                // dropped isoform (multi-parent Parent=a,b) stays kept.
                for (int iso_idx : isoform_indices) {
                    const auto& iso = data.records[iso_idx];
                    if (!iso.id) continue;
                    auto child_it = isoform_to_children.find(*iso.id);
                    if (child_it != isoform_to_children.end()) {
                        for (int child_idx : child_it->second) {
                            data.records[child_idx].kept = false;
                        }
                    }
                }
                const auto& longest = data.records[longest_idx];
                if (longest.id) {
                    auto child_it = isoform_to_children.find(*longest.id);
                    if (child_it != isoform_to_children.end()) {
                        for (int child_idx : child_it->second) {
                            data.records[child_idx].kept = true;
                        }
                    }
                }
            }
        }
    };

    if (num_threads <= 1) {
        for (auto& [chrom, gene_indices] : chrom_to_gene_idx) {
            process_chromosome(chrom, gene_indices);
        }
    } else {
        std::vector<std::future<void>> futures;
        for (auto& kv : chrom_to_gene_idx) {
            const std::string& chrom = kv.first;
            std::vector<int> gene_indices = kv.second;
            auto f = std::async(std::launch::async, [&, chrom, gene_indices]() {
                process_chromosome(chrom, gene_indices);
            });
            futures.push_back(std::move(f));
        }
        for (auto& f : futures) {
            f.get();
        }
    }
}

}  // namespace gffsub

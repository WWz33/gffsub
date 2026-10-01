#include "filter.hpp"
#include "parser.hpp"
#include "record.hpp"
#include <algorithm>
#include <future>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <unordered_set>

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
            if (rec.feat_class == FeatureClass::Transcript) ++counts[std::string{rec.type}];
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

    // Genes are the Parent of the competing isoforms. They come from a gene
    // (L1) record, or, when the file has none (flat GTF, GFF3 without gene
    // rows), from the parent ID attribute itself. Both are planned the same
    // way; only an actual gene record can be dropped as childless.
    struct ChromGroup {
        std::vector<int> gene_indices;
        std::vector<std::string> virtual_gene_ids;
    };
    std::unordered_map<std::string, ChromGroup> chrom_to_genes;
    std::unordered_set<std::string> genes_with_record;
    for (int i = 0; i < static_cast<int>(data.records.size()); ++i) {
        const auto& rec = data.records[i];
        if (!rec.kept) continue;
        if (rec.feat_class == FeatureClass::Gene && rec.id) {
            genes_with_record.insert(*rec.id);
            chrom_to_genes[std::string{rec.seqid}].gene_indices.push_back(i);
        }
    }
    for (const auto& [gene_id, iso_indices] : gene_to_isoforms) {
        if (iso_indices.empty() || genes_with_record.count(gene_id) > 0) continue;
        // Bucket by the seqid of the first isoform so the per-chromosome
        // batching below applies to recordless genes too.
        chrom_to_genes[std::string{data.records[iso_indices.front()].seqid}]
            .virtual_gene_ids.push_back(gene_id);
    }

    // Two-phase design: phase 1 computes each gene's decision read-only
    // (parallel-safe), phase 2 applies every write single-threaded. Genes on
    // different chromosomes can share children via multi-parent transcripts
    // (e.g. trans-splicing), so workers must never write record state.
    struct GenePlan {
        int gene_idx = -1;
        bool drop_gene = false;
        int longest_idx = -1;
        std::vector<int> isoform_indices;
    };

    auto plan_gene = [&](const std::string& gene_id, int gene_idx) -> std::optional<GenePlan> {
        auto isoform_it = gene_to_isoforms.find(gene_id);
        if (isoform_it == gene_to_isoforms.end()) {
            // gene with no isoform children: keep only if it has any
            // non-isoform children (e.g. TF_binding_site). A gene with
            // no children at all is dropped.
            auto child_it = isoform_to_children.find(gene_id);
            if (child_it == isoform_to_children.end()) {
                GenePlan plan;
                plan.gene_idx = gene_idx;
                plan.drop_gene = true;
                return plan;
            }
            return std::nullopt;
        }
        if (isoform_it->second.size() <= 1) return std::nullopt;

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
            // Segment lengths are int64 and summed without clamping, so a
            // pathological record near the coordinate limit must saturate
            // instead of wrapping negative and beating every real isoform.
            constexpr int64_t kLenMax = std::numeric_limits<int64_t>::max();
            const auto add_segment = [&](int64_t seg) {
                len = (seg > kLenMax - len) ? kLenMax : len + seg;
            };

            if (gene_has_cds) {
                // CDS length is the sum of every CDS segment under the
                // transcript. Segments either share an ID (one discontinuous
                // CDS) or carry distinct IDs (e.g. a CDS split by a
                // translational frameshift, as in the GFF3 spec example);
                // both count towards the transcript's coding length.
                for (int child_idx : child_it->second) {
                    const auto& child = data.records[child_idx];
                    if (child.feat_class == FeatureClass::CDS) {
                        add_segment(child.end - child.start + 1);
                    }
                }
                if (len == 0) continue; // isoform without CDS is skipped
            } else {
                for (int child_idx : child_it->second) {
                    const auto& child = data.records[child_idx];
                    if (child.feat_class == FeatureClass::Exon) {
                        add_segment(child.end - child.start + 1);
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

        GenePlan plan;
        plan.gene_idx = gene_idx;
        plan.longest_idx = longest_idx;
        plan.isoform_indices = isoform_indices;
        return plan;
    };

    auto plan_chromosome = [&](const ChromGroup& group) {
        std::vector<GenePlan> plans;
        for (int gene_idx : group.gene_indices) {
            const auto& gene = data.records[gene_idx];
            if (!gene.id) continue;
            if (auto plan = plan_gene(*gene.id, gene_idx)) {
                plans.push_back(std::move(*plan));
            }
        }
        for (const auto& gene_id : group.virtual_gene_ids) {
            if (auto plan = plan_gene(gene_id, -1)) {
                plans.push_back(std::move(*plan));
            }
        }
        return plans;
    };

    std::vector<GenePlan> all_plans;
    if (num_threads <= 1) {
        for (auto& [chrom, group] : chrom_to_genes) {
            (void)chrom;
            auto plans = plan_chromosome(group);
            all_plans.insert(all_plans.end(),
                             std::make_move_iterator(plans.begin()),
                             std::make_move_iterator(plans.end()));
        }
    } else {
        // Cap in-flight threads at num_threads; scaffold-heavy files can
        // have thousands of chromosomes, one future each would exhaust
        // thread limits. Batch: launch up to num_threads, then wait.
        std::vector<std::future<std::vector<GenePlan>>> futures;
        for (auto& kv : chrom_to_genes) {
            ChromGroup group = kv.second;
            futures.push_back(std::async(std::launch::async, [&, group]() {
                return plan_chromosome(group);
            }));
            if (futures.size() >= num_threads) {
                for (auto& f : futures) {
                    auto plans = f.get();
                    all_plans.insert(all_plans.end(),
                                     std::make_move_iterator(plans.begin()),
                                     std::make_move_iterator(plans.end()));
                }
                futures.clear();
            }
        }
        for (auto& f : futures) {
            auto plans = f.get();
            all_plans.insert(all_plans.end(),
                             std::make_move_iterator(plans.begin()),
                             std::make_move_iterator(plans.end()));
        }
    }

    // Phase 2: apply writes single-threaded. Global drop pass then global
    // re-keep pass. Only the re-keep pass writes true, so a record claimed by
    // any gene's longest isoform stays kept regardless of plan order; this
    // covers both multi-parent children and a multi-parent isoform that is
    // one gene's winner and another gene's loser.
    for (const auto& plan : all_plans) {
        if (plan.drop_gene) {
            data.records[plan.gene_idx].kept = false;
        }
    }
    for (const auto& plan : all_plans) {
        if (plan.longest_idx < 0) continue;
        for (int iso_idx : plan.isoform_indices) {
            data.records[iso_idx].kept = false;
            const auto& iso = data.records[iso_idx];
            if (!iso.id) continue;
            auto child_it = isoform_to_children.find(*iso.id);
            if (child_it != isoform_to_children.end()) {
                for (int child_idx : child_it->second) {
                    data.records[child_idx].kept = false;
                }
            }
        }
    }
    for (const auto& plan : all_plans) {
        if (plan.longest_idx < 0) continue;
        data.records[plan.longest_idx].kept = true;
        const auto& longest = data.records[plan.longest_idx];
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

}  // namespace gffsub

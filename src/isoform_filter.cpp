#include "filter.hpp"
#include "parser.hpp"
#include "record.hpp"
#include <algorithm>
#include <future>
#include <unordered_map>

namespace gffsub {

void filter_longest_isoform(GffData& data, IdIndex& /*idx*/, std::string_view feature_type, size_t num_threads) {
    std::string isoform_type{feature_type};
    if (isoform_type.empty()) {
        // Auto-detect: GFF3 uses "mRNA", GTF uses "transcript".
        // Pick whichever exists in the data, preferring mRNA.
        bool has_mrna = false;
        bool has_transcript = false;
        for (const auto& rec : data.records) {
            if (!rec.kept) continue;
            if (rec.type == "mRNA") has_mrna = true;
            else if (rec.type == "transcript") has_transcript = true;
        }
        isoform_type = has_mrna ? "mRNA" : (has_transcript ? "transcript" : "mRNA");
    }

    // Build gene -> [isoform indices] index once
    std::unordered_map<std::string, std::vector<int>> gene_to_isoforms;
    for (int i = 0; i < static_cast<int>(data.records.size()); ++i) {
        const auto& rec = data.records[i];
        if (!rec.kept) continue;
        if (rec.type == isoform_type && rec.parent_id) {
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
        if (rec.type == "gene" && rec.id) {
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
                // gene with no isoform children is dropped
                data.records[gene_idx].kept = false;
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
                        if (data.records[child_idx].type == "CDS") {
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
                    for (int child_idx : child_it->second) {
                        const auto& child = data.records[child_idx];
                        if (child.type == "CDS") {
                            len += child.end - child.start + 1;
                            found = true;
                        }
                    }
                    if (!found) continue; // isoform without CDS is skipped
                } else {
                    for (int child_idx : child_it->second) {
                        const auto& child = data.records[child_idx];
                        if (child.type == "exon") {
                            len += child.end - child.start + 1;
                            found = true;
                        }
                    }
                    if (!found) {
                        len = iso.end - iso.start + 1;
                    }
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

    if (!feature_type.empty()) {
        filter_by_feature(data, feature_type);
    }
}

}  // namespace gffsub

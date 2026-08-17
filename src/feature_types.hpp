#ifndef GFFSUB_FEATURE_TYPES_HPP
#define GFFSUB_FEATURE_TYPES_HPP

#include <cstdint>
#include <string_view>

namespace gffsub {

enum class OutputFormat { GFF3, GTF2, GTF3, BED };
enum class InputFormat { GFF3, GTF, BED };

// Coarse biological class of a GFF/GTF feature type, derived once at parse
// time by classify_type(). Covers the SO feature space via suffix matching
// (gffread-style), so miRNA/tRNA/pseudogene/protein_coding_gene etc. are
// recognized without a static type list.
enum class FeatureClass : uint8_t {
    Unknown = 0,
    Gene,         // gene, pseudogene, protein_coding_gene, *_gene
    Transcript,   // mRNA, transcript, miRNA, tRNA, ncRNA, *_transcript
    Exon,         // exon, pseudogenic_exon
    CDS,          // CDS, pseudogenic_cds
    UTR,          // *utr (five_prime_utr, three_prime_utr, 3utr, 5utr, UTR)
    Match,        // cDNA_match, protein_match, *_match (alignment features)
    Region        // BED-derived region rows
};

// Type string used for BED rows (BED has no type column).
inline constexpr std::string_view kRegionType = "region";

// Classify a feature type string into FeatureClass. Suffix-based, case-
// insensitive. Order matters: match/utr/exon/transcript checked before gene
// so e.g. gene_segment (ends "segment") and transcript_region (ends "region")
// fall through to Unknown as intended. AGAT feature_levels.yaml level3
// values (e.g. sig_peptide: exon) describe the PARENT type, not the self
// type, so they do not influence this classification. SO accessions
// (SO:0000234) are not recognized.
FeatureClass classify_type(std::string_view type);

// True if the type is an AGAT feature_levels.yaml "spread" feature: one
// biological feature that may span multiple disjunct lines sharing an ID
// (GFF3 spec: multiline features must have an ID).
bool is_spread_feature(std::string_view type);

// Resolve the GTF column-3 label for a feature type. GTF3 rewrites "mRNA" to
// "transcript" and normalizes UTR spellings to five_prime_utr /
// three_prime_utr. GTF2 normalizes UTR to 5UTR / 3UTR. Every other type
// passes through unchanged.
std::string_view gtf_type_label(std::string_view type, OutputFormat fmt);

// Should this type survive into the given output format? GTF3: strict
// whitelist (case-insensitive, so GFF3 "five_prime_UTR"/"CDS" match). GTF2:
// passthrough. GFF3/BED: always true.
bool gtf_type_emittable(std::string_view type, OutputFormat fmt);

}  // namespace gffsub

#endif  // GFFSUB_FEATURE_TYPES_HPP

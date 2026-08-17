// Unit tests for feature_types module: classify_type, gtf_type_label,
// gtf_type_emittable, is_spread_feature.
#include "feature_types.hpp"

#include <cassert>
#include <cstdio>
#include <string>
#include <string_view>

using namespace gffsub;

static int failures = 0;

#define CHECK(expr) \
    do { if (!(expr)) { \
        std::fprintf(stderr, "FAIL %s:%d: %s\n", __FILE__, __LINE__, #expr); \
        ++failures; \
    } } while (0)

static void check_eq(FeatureClass got, FeatureClass want, const char* label) {
    if (got != want) {
        std::fprintf(stderr, "FAIL classify_type(%s): got %d, want %d\n", label, static_cast<int>(got), static_cast<int>(want));
        ++failures;
    }
}

// ---------------------------------------------------------------------------
// classify_type
// ---------------------------------------------------------------------------

static void test_classify_gene() {
    check_eq(classify_type("gene"), FeatureClass::Gene, "gene");
    check_eq(classify_type("pseudogene"), FeatureClass::Gene, "pseudogene");
    check_eq(classify_type("protein_coding_gene"), FeatureClass::Gene, "protein_coding_gene");
    check_eq(classify_type("mirna_gene"), FeatureClass::Gene, "mirna_gene");
    check_eq(classify_type("ncrna_gene"), FeatureClass::Gene, "ncrna_gene");
    check_eq(classify_type("transposable_element_gene"), FeatureClass::Gene, "transposable_element_gene");
    check_eq(classify_type("GENE"), FeatureClass::Gene, "GENE");
    check_eq(classify_type("Pseudogene"), FeatureClass::Gene, "Pseudogene");
}

static void test_classify_transcript() {
    check_eq(classify_type("mRNA"), FeatureClass::Transcript, "mRNA");
    check_eq(classify_type("transcript"), FeatureClass::Transcript, "transcript");
    check_eq(classify_type("miRNA"), FeatureClass::Transcript, "miRNA");
    check_eq(classify_type("tRNA"), FeatureClass::Transcript, "tRNA");
    check_eq(classify_type("ncRNA"), FeatureClass::Transcript, "ncRNA");
    check_eq(classify_type("lnc_RNA"), FeatureClass::Transcript, "lnc_RNA");
    check_eq(classify_type("rRNA"), FeatureClass::Transcript, "rRNA");
    check_eq(classify_type("snoRNA"), FeatureClass::Transcript, "snoRNA");
    check_eq(classify_type("pseudogenic_transcript"), FeatureClass::Transcript, "pseudogenic_transcript");
    check_eq(classify_type("primary_transcript"), FeatureClass::Transcript, "primary_transcript");
    check_eq(classify_type("MRNA"), FeatureClass::Transcript, "MRNA");
    check_eq(classify_type("MiRNA"), FeatureClass::Transcript, "MiRNA");
}

static void test_classify_exon() {
    check_eq(classify_type("exon"), FeatureClass::Exon, "exon");
    check_eq(classify_type("pseudogenic_exon"), FeatureClass::Exon, "pseudogenic_exon");
    check_eq(classify_type("EXON"), FeatureClass::Exon, "EXON");
}

static void test_classify_cds() {
    check_eq(classify_type("CDS"), FeatureClass::CDS, "CDS");
    check_eq(classify_type("pseudogenic_cds"), FeatureClass::CDS, "pseudogenic_cds");
    check_eq(classify_type("cds"), FeatureClass::CDS, "cds");
}

static void test_classify_utr() {
    check_eq(classify_type("five_prime_UTR"), FeatureClass::UTR, "five_prime_UTR");
    check_eq(classify_type("three_prime_utr"), FeatureClass::UTR, "three_prime_utr");
    check_eq(classify_type("5UTR"), FeatureClass::UTR, "5UTR");
    check_eq(classify_type("3utr"), FeatureClass::UTR, "3utr");
    check_eq(classify_type("UTR"), FeatureClass::UTR, "UTR");
    check_eq(classify_type("5'-utr"), FeatureClass::UTR, "5'-utr");
}

static void test_classify_match() {
    check_eq(classify_type("cDNA_match"), FeatureClass::Match, "cDNA_match");
    check_eq(classify_type("protein_match"), FeatureClass::Match, "protein_match");
    check_eq(classify_type("match"), FeatureClass::Match, "match");
    check_eq(classify_type("est_match"), FeatureClass::Match, "est_match");
}

static void test_classify_order_sensitive() {
    // gene_segment ends "segment", not "gene" -> Unknown (not Gene)
    check_eq(classify_type("gene_segment"), FeatureClass::Unknown, "gene_segment");
    check_eq(classify_type("v_gene_segment"), FeatureClass::Unknown, "v_gene_segment");
    // transcript_region ends "region", not "transcript" -> Unknown
    check_eq(classify_type("transcript_region"), FeatureClass::Unknown, "transcript_region");
    // match checked before utr/exon — no conflict in practice, but verify
    check_eq(classify_type("cdna_match"), FeatureClass::Match, "cdna_match");
}

static void test_classify_unknown() {
    check_eq(classify_type(""), FeatureClass::Unknown, "empty");
    check_eq(classify_type("region"), FeatureClass::Unknown, "region");
    check_eq(classify_type("promoter"), FeatureClass::Unknown, "promoter");
    check_eq(classify_type("enhancer"), FeatureClass::Unknown, "enhancer");
    check_eq(classify_type("intron"), FeatureClass::Unknown, "intron");
    check_eq(classify_type("start_codon"), FeatureClass::Unknown, "start_codon");
    check_eq(classify_type("stop_codon"), FeatureClass::Unknown, "stop_codon");
    check_eq(classify_type("sig_peptide"), FeatureClass::Unknown, "sig_peptide");
}

// ---------------------------------------------------------------------------
// gtf_type_label
// ---------------------------------------------------------------------------

static void test_gtf_type_label_mrna() {
    CHECK(gtf_type_label("mRNA", OutputFormat::GTF3) == "transcript");
    CHECK(gtf_type_label("mRNA", OutputFormat::GTF2) == "mRNA");
    CHECK(gtf_type_label("mRNA", OutputFormat::GFF3) == "mRNA");
}

static void test_gtf_type_label_passthrough() {
    CHECK(gtf_type_label("exon", OutputFormat::GTF3) == "exon");
    CHECK(gtf_type_label("CDS", OutputFormat::GTF2) == "CDS");
    CHECK(gtf_type_label("gene", OutputFormat::GTF3) == "gene");
    CHECK(gtf_type_label("miRNA", OutputFormat::GTF3) == "miRNA");
}

static void test_gtf_type_label_utr() {
    // GTF3 normalizes UTR spellings
    CHECK(gtf_type_label("five_prime_UTR", OutputFormat::GTF3) == "five_prime_utr");
    CHECK(gtf_type_label("three_prime_UTR", OutputFormat::GTF3) == "three_prime_utr");
    CHECK(gtf_type_label("5UTR", OutputFormat::GTF3) == "five_prime_utr");
    CHECK(gtf_type_label("3utr", OutputFormat::GTF3) == "three_prime_utr");
    // GTF2 uses 5UTR/3UTR form
    CHECK(gtf_type_label("five_prime_UTR", OutputFormat::GTF2) == "5UTR");
    CHECK(gtf_type_label("three_prime_utr", OutputFormat::GTF2) == "3UTR");
    CHECK(gtf_type_label("5UTR", OutputFormat::GTF2) == "5UTR");
    // Plain UTR stays as-is (can't determine 5' or 3')
    CHECK(gtf_type_label("UTR", OutputFormat::GTF3) == "UTR");
}

// ---------------------------------------------------------------------------
// gtf_type_emittable
// ---------------------------------------------------------------------------

static void test_gtf_type_emittable() {
    // GTF3 whitelist (case-insensitive)
    CHECK(gtf_type_emittable("gene", OutputFormat::GTF3));
    CHECK(gtf_type_emittable("transcript", OutputFormat::GTF3));
    CHECK(gtf_type_emittable("exon", OutputFormat::GTF3));
    CHECK(gtf_type_emittable("CDS", OutputFormat::GTF3));
    CHECK(gtf_type_emittable("start_codon", OutputFormat::GTF3));
    CHECK(gtf_type_emittable("stop_codon", OutputFormat::GTF3));
    CHECK(gtf_type_emittable("mRNA", OutputFormat::GTF3));
    CHECK(gtf_type_emittable("five_prime_UTR", OutputFormat::GTF3));  // case-insensitive bug fix
    CHECK(gtf_type_emittable("three_prime_utr", OutputFormat::GTF3));
    CHECK(gtf_type_emittable("Selenocysteine", OutputFormat::GTF3));
    CHECK(gtf_type_emittable("inter", OutputFormat::GTF3));
    CHECK(gtf_type_emittable("inter_CNS", OutputFormat::GTF3));
    CHECK(gtf_type_emittable("intron_CNS", OutputFormat::GTF3));
    // Non-whitelisted types filtered
    CHECK(!gtf_type_emittable("miRNA", OutputFormat::GTF3));
    CHECK(!gtf_type_emittable("promoter", OutputFormat::GTF3));
    CHECK(!gtf_type_emittable("pseudogene", OutputFormat::GTF3));
    // GTF2/GFF3/BED: passthrough
    CHECK(gtf_type_emittable("miRNA", OutputFormat::GTF2));
    CHECK(gtf_type_emittable("anything", OutputFormat::GFF3));
    CHECK(gtf_type_emittable("anything", OutputFormat::BED));
}

// ---------------------------------------------------------------------------
// is_spread_feature
// ---------------------------------------------------------------------------

static void test_is_spread_feature() {
    CHECK(is_spread_feature("CDS"));
    CHECK(is_spread_feature("pseudogenic_cds"));
    CHECK(is_spread_feature("five_prime_utr"));
    CHECK(is_spread_feature("three_prime_UTR"));
    CHECK(is_spread_feature("utr"));
    CHECK(is_spread_feature("3utr"));
    CHECK(is_spread_feature("start_codon"));
    CHECK(is_spread_feature("stop_codon"));
    CHECK(is_spread_feature("uorf"));
    // Non-spread
    CHECK(!is_spread_feature("exon"));
    CHECK(!is_spread_feature("gene"));
    CHECK(!is_spread_feature("mRNA"));
}

// ---------------------------------------------------------------------------

int main() {
    test_classify_gene();
    test_classify_transcript();
    test_classify_exon();
    test_classify_cds();
    test_classify_utr();
    test_classify_match();
    test_classify_order_sensitive();
    test_classify_unknown();

    test_gtf_type_label_mrna();
    test_gtf_type_label_passthrough();
    test_gtf_type_label_utr();

    test_gtf_type_emittable();
    test_is_spread_feature();

    if (failures == 0) {
        std::printf("feature_types_smoke OK\n");
        return 0;
    }
    std::fprintf(stderr, "feature_types_smoke FAILED: %d assertions\n", failures);
    return 1;
}

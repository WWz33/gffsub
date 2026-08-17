#include "feature_types.hpp"
#include "string_utils.hpp"

#include <string>
#include <unordered_set>

namespace gffsub {

namespace {

bool ends_with(std::string_view s, std::string_view suffix) {
    return s.size() >= suffix.size() &&
           s.compare(s.size() - suffix.size(), suffix.size(), suffix) == 0;
}

}  // namespace

FeatureClass classify_type(std::string_view type) {
    const std::string lc = to_lower(type);
    if (ends_with(lc, "match"))                                     return FeatureClass::Match;
    if (lc.find("utr") != std::string::npos)                       return FeatureClass::UTR;
    if (ends_with(lc, "exon"))                                      return FeatureClass::Exon;
    if (ends_with(lc, "rna") || ends_with(lc, "transcript"))       return FeatureClass::Transcript;
    if (ends_with(lc, "gene"))                                     return FeatureClass::Gene;
    if (ends_with(lc, "cds"))                                      return FeatureClass::CDS;
    return FeatureClass::Unknown;
}

bool is_spread_feature(std::string_view type) {
    const std::string lc = to_lower(type);
    if (ends_with(lc, "cds"))  return true;
    if (ends_with(lc, "utr"))  return true;
    if (lc == "start_codon" || lc == "stop_codon" || lc == "uorf") return true;
    return false;
}

std::string_view gtf_type_label(std::string_view type, OutputFormat fmt) {
    if (fmt != OutputFormat::GTF2 && fmt != OutputFormat::GTF3) {
        return type;
    }
    const std::string lc = to_lower(type);

    // mRNA -> transcript (GTF3 only; GTF2 keeps "mRNA")
    if (fmt == OutputFormat::GTF3 && lc == "mrna") {
        return "transcript";
    }

    // UTR spelling normalization per GTF version
    // GFF3 uses five_prime_UTR / three_prime_UTR (mixed case)
    if (lc == "five_prime_utr" || lc == "5utr" || lc == "5'-utr") {
        return (fmt == OutputFormat::GTF3) ? "five_prime_utr" : "5UTR";
    }
    if (lc == "three_prime_utr" || lc == "3utr" || lc == "3'-utr") {
        return (fmt == OutputFormat::GTF3) ? "three_prime_utr" : "3UTR";
    }

    return type;
}

bool gtf_type_emittable(std::string_view type, OutputFormat fmt) {
    if (fmt != OutputFormat::GTF3) {
        return true;
    }
    // GTF2.2 spec: required CDS/start_codon/stop_codon, optional
    // 5UTR/3UTR/exon/inter/inter_CNS/intron_CNS. GTF3 (Ensembl) extends with
    // gene/transcript/five_prime_utr/three_prime_utr/Selenocysteine/mRNA.
    // Case-insensitive so GFF3 "five_prime_UTR"/"CDS" match.
    static const std::unordered_set<std::string> gtf3_whitelist = {
        "gene", "transcript", "exon", "cds", "start_codon", "stop_codon",
        "five_prime_utr", "three_prime_utr", "selenocysteine", "mrna",
        "5utr", "3utr", "inter", "inter_cns", "intron_cns"
    };
    return gtf3_whitelist.count(to_lower(type)) > 0;
}

}  // namespace gffsub

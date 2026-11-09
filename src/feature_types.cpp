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
    if (fmt != OutputFormat::GTF3 && fmt != OutputFormat::GTF25 && fmt != OutputFormat::GTF22) {
        return type;
    }
    const std::string lc = to_lower(type);

    // mRNA -> transcript. GTF2.2 has no transcript feature at all, so the
    // rename is moot there: the whitelist drops the row.
    if (lc == "mrna" && (fmt == OutputFormat::GTF3 || fmt == OutputFormat::GTF25)) {
        return "transcript";
    }

    // UTR spelling per GTF version
    // GFF3 uses five_prime_UTR / three_prime_UTR (mixed case)
    if (lc == "five_prime_utr" || lc == "5utr" || lc == "5'-utr") {
        if (fmt == OutputFormat::GTF3)  return "five_prime_utr";
        if (fmt == OutputFormat::GTF25) return "UTR";
        return "5UTR";
    }
    if (lc == "three_prime_utr" || lc == "3utr" || lc == "3'-utr") {
        if (fmt == OutputFormat::GTF3)  return "three_prime_utr";
        if (fmt == OutputFormat::GTF25) return "UTR";
        return "3UTR";
    }

    // GTF2.2 requires exact capitalization for the accepted types. A GFF3
    // file that spells them differently is normalized here.
    if (lc == "gene")           return "gene";
    if (lc == "exon")           return "exon";
    if (lc == "transcript")     return "transcript";
    if (lc == "cds")            return "CDS";
    if (lc == "start_codon")    return "start_codon";
    if (lc == "stop_codon")     return "stop_codon";
    if (lc == "selenocysteine") return "Selenocysteine";

    return type;
}

std::string_view gff3_type_label(std::string_view type) {
    const std::string lc = to_lower(type);
    if (lc == "5utr" || lc == "5'-utr" || lc == "five_prime_utr") {
        return "five_prime_UTR";
    }
    if (lc == "3utr" || lc == "3'-utr" || lc == "three_prime_utr") {
        return "three_prime_UTR";
    }
    return type;
}

bool gtf_type_emittable(std::string_view type, OutputFormat fmt) {
    // GTF2.2 spec: required CDS/start_codon/stop_codon, optional
    // 5UTR/3UTR/exon/inter/inter_CNS/intron_CNS. GTF2.5 (Ensembl release 75)
    // adds gene/transcript/Selenocysteine and collapses the two UTRs into
    // "UTR". GTF3 (release 82) splits them back into five_prime_utr /
    // three_prime_utr. Case-insensitive so GFF3 "five_prime_UTR"/"CDS" match.
    static const std::unordered_set<std::string> gtf3_whitelist = {
        "gene", "transcript", "exon", "cds", "start_codon", "stop_codon",
        "five_prime_utr", "three_prime_utr", "selenocysteine"
    };
    static const std::unordered_set<std::string> gtf25_whitelist = {
        "gene", "transcript", "exon", "cds", "start_codon", "stop_codon",
        "utr", "selenocysteine"
    };
    static const std::unordered_set<std::string> gtf22_whitelist = {
        "cds", "start_codon", "stop_codon", "exon", "5utr", "3utr",
        "inter", "inter_cns", "intron_cns"
    };

    const std::string lc = to_lower(type);
    switch (fmt) {
        case OutputFormat::GTF3:  return gtf3_whitelist.count(lc) > 0;
        case OutputFormat::GTF25: return gtf25_whitelist.count(lc) > 0;
        case OutputFormat::GTF22: return gtf22_whitelist.count(lc) > 0;
        default:                  return true;
    }
}

}  // namespace gffsub

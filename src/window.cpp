#include "window.hpp"

#include "region.hpp"

#include <optional>
#include <stdexcept>

namespace gffsub {

GffData window(const AnnotationIndex& index, const WindowParams& params) {
    // For discontinuous features (multi-line CDS sharing one ID), compute the
    // span across all segments so the window covers the whole feature.
    std::optional<GffRecord> target;
    auto all_segments = index.find_all_by_id(params.id);
    if (!all_segments.empty()) {
        target = all_segments[0];
        for (const auto& seg : all_segments) {
            if (seg.start < target->start) target->start = seg.start;
            if (seg.end > target->end) target->end = seg.end;
        }
    }
    if (!target) {
        target = index.find_gene(params.id);
    }
    if (!target) {
        throw std::runtime_error("cannot find " + params.id);
    }

    const auto region = window_region(*target, params.upstream, params.downstream, params.strand_aware);
    GffData result;
    for (const auto& rec : index.overlap(region.seqid, region.start, region.end)) {
        result.append(rec);
    }
    return result;
}

}  // namespace gffsub

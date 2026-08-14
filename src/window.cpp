#include "window.hpp"

#include "region.hpp"

#include <optional>
#include <stdexcept>

namespace gffsub {

GffData window(const AnnotationIndex& index, const WindowParams& params) {
    auto target = index.find_by_id(params.id);
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

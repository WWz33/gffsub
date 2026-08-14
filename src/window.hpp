#ifndef WINDOW_HPP
#define WINDOW_HPP

#include "annotation.hpp"
#include "record.hpp"

#include <cstdint>
#include <string>

namespace gffsub {

struct WindowParams {
    std::string id;
    int64_t upstream = 0;
    int64_t downstream = 0;
    bool strand_aware = false;
};

// Find a feature by ID (falls back to gene name lookup), build a window
// around it, and return all overlapping records.
GffData window(const AnnotationIndex& index, const WindowParams& params);

}  // namespace gffsub

#endif  // WINDOW_HPP

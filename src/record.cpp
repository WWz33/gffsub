// GffData storage-layer definitions split out of record.hpp so the header
// does not pull mmap/unistd into every TU that includes it, and library
// consumers that touch GffData without gff3_parser.o still link.
#include "record.hpp"

#include <sys/mman.h>

namespace gffsub {

void GffData::release_mapping() {
    if (mapped_data_ != nullptr) {
        ::munmap(mapped_data_, mapped_size_);
        mapped_data_ = nullptr;
        mapped_size_ = 0;
    }
}

}  // namespace gffsub

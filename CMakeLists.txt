cmake_minimum_required(VERSION 3.10)
project(gffsub VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(gffsub
    src/gffsub.cpp
    src/gff3_parser.cpp
    src/gff3_filter.cpp
    src/gff3_output.cpp
)

target_include_directories(gffsub PRIVATE
    src
)

target_compile_options(gffsub PRIVATE
    -Wall -O2
)

# Installation
install(TARGETS gffsub DESTINATION bin)

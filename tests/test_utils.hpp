// Shared test helpers for gffsub smoke tests. Header-only so each test
// binary stays self-contained with no extra link step.
#ifndef GFFSUB_TEST_UTILS_HPP
#define GFFSUB_TEST_UTILS_HPP

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace test_utils {

inline std::string read_file(const std::string& path) {
    std::ifstream in{path};
    std::ostringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}

inline int run_command(const std::string& command) {
    const int status = std::system(command.c_str());
    if (status != 0) {
        std::cerr << "command failed: " << command << '\n';
    }
    return status;
}

// Expect a command to FAIL (non-zero exit). Returns 0 on success of the
// expectation, 1 if the command unexpectedly succeeded.
inline int expect_command_failure(const std::string& command) {
    const int status = std::system(command.c_str());
    if (status == 0) {
        std::cerr << "command unexpectedly succeeded: " << command << '\n';
        return 1;
    }
    return 0;
}

inline bool contains(const std::string& text, const std::string& needle) {
    return text.find(needle) != std::string::npos;
}

inline int require_contains(const std::string& path, const std::string& needle) {
    const auto text = read_file(path);
    if (!contains(text, needle)) {
        std::cerr << "missing '" << needle << "' in " << path << '\n';
        return 1;
    }
    return 0;
}

inline int require_not_contains(const std::string& path, const std::string& needle) {
    const auto text = read_file(path);
    if (contains(text, needle)) {
        std::cerr << "unexpected '" << needle << "' in " << path << '\n';
        return 1;
    }
    return 0;
}

inline int compare_files(const std::string& lhs_path, const std::string& rhs_path) {
    const auto lhs = read_file(lhs_path);
    const auto rhs = read_file(rhs_path);
    if (lhs != rhs) {
        std::cerr << "output mismatch: " << lhs_path << " vs " << rhs_path << '\n';
        return 1;
    }
    return 0;
}

inline int require_exit_one_with_error(const std::string& command,
                                       const std::string& err_path,
                                       const std::string& message) {
    if (expect_command_failure(command) != 0 ||
        require_contains(err_path, "Error:") != 0 ||
        require_contains(err_path, message) != 0) {
        return 1;
    }
    return 0;
}

}  // namespace test_utils

#endif  // GFFSUB_TEST_UTILS_HPP

// corefine_faults_hello — Phase 0 build-environment smoke binary
// for PLAN_cgal_corefine.md.
//
// Verifies that:
//   1. CGAL ≥ 5.4 (tested 5.6.1) is available as headers,
//   2. CGAL::Polygon_mesh_processing::corefinement and
//      CGAL::Surface_mesh include cleanly,
//   3. Boost headers and GMP/MPFR (via CGAL::CGAL) link cleanly,
//   4. C++17 (`std::string_view`, `std::filesystem`) is available.
//
// No CLI arguments.  Exits 0 on success and prints a single line
// `CGAL version: <X.Y>` to stdout for the Phase-0 acceptance test.
//
// Phase 1 will replace this binary with the real `corefine_faults`
// tool that drives `CGAL::Polygon_mesh_processing::corefine`.

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string_view>

#include <CGAL/version.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

namespace fs = std::filesystem;

int main(int argc, char** argv) {
    (void)argc;
    (void)argv;

    // C++17 features used elsewhere in the tool — exercise them here so
    // the environment is verified end-to-end at Phase 0.
    constexpr std::string_view kPrefix = "CGAL version: ";
    const fs::path cwd = fs::current_path();

    // Instantiate the kernel + Surface_mesh + a corefinement-related
    // typedef so the headers are actually compiled (not just included).
    using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
    using Mesh = CGAL::Surface_mesh<Kernel::Point_3>;
    [[maybe_unused]] Mesh m;

    // CGAL 5.x exposes the version both as a string and as a numeric
    // macro; print the string form for human readability.
    std::cout << kPrefix << CGAL_VERSION_STR << '\n';

    // Optional second line, kept short so test grep can match the
    // first line by prefix without ambiguity.
    std::cerr << "[corefine_faults_hello] cwd=" << cwd.string() << '\n';

    return EXIT_SUCCESS;
}

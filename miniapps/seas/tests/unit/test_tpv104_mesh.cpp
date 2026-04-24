// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Unit tests for TPV104 mesh files (§4.10 Step 10 gates T_TPV104_MESH_0..3).
//
// T_TPV104_MESH_0: byte-identical geometry invariant vs TPV102 reference.
// T_TPV104_MESH_1: vertex / element / physical-group counts within 5 %
//                  at matched resolution.  DEFERRED — requires .msh build.
// T_TPV104_MESH_2: wave-speed arrival on the coarse mesh.  DEFERRED —
//                  requires the Step 9 driver + solver to be available.
// T_TPV104_MESH_3: physical-group parity with driver.  DEFERRED — same.

#include "test_macros.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;

static std::string locate(const std::string &tail)
{
   const std::vector<std::string> prefixes = {"", "../", "../../"};
   for (const auto &p : prefixes)
   {
      std::ifstream f(p + tail);
      if (f.is_open()) { return p + tail; }
   }
   return "";
}

static std::string slurp(const std::string &path)
{
   std::ifstream f(path);
   std::ostringstream oss;
   oss << f.rdbuf();
   return oss.str();
}

// R3-005 (review round 3): sentinel-anchor tail extraction.
// Locate the first line that begins with `"lc = "` and return everything
// from that line onward.  This is immune to header-comment-block length
// changes and handles CRLF line endings by stripping trailing `\r`.
// A hardcoded line-number cutoff (e.g. "lines >= 30") would silently
// desync if a header were edited on only one of the two files being
// diffed.
static std::string TailFromCode(const std::string &path)
{
   std::ifstream in(path);
   std::string line;
   std::ostringstream oss;
   bool in_code = false;
   while (std::getline(in, line))
   {
      if (!line.empty() && line.back() == '\r') { line.pop_back(); }
      if (!in_code && line.rfind("lc = ", 0) == 0) { in_code = true; }
      if (in_code) { oss << line << "\n"; }
   }
   return oss.str();
}

// ----------------------------------------------------------------------------
// T_TPV104_MESH_0 — byte-identical geometry invariant.
// Diff tpv104_200m.geo against tpv102_200m.geo starting from line 30.
// Lines 1-29 are the header comment block (per plan §4.10 Step 10).
// ----------------------------------------------------------------------------
void TestGeometryByteIdentical()
{
   std::cout << "\n[T_TPV104_MESH_0] byte-identical geometry vs TPV102\n";

   const std::string tpv104_path = locate("tpv104/mesh/tpv104_200m.geo");
   const std::string tpv102_path = locate("tpv102/mesh/tpv102_200m.geo");

   TEST_ASSERT(!tpv104_path.empty(), "locate tpv104/mesh/tpv104_200m.geo");
   TEST_ASSERT(!tpv102_path.empty(), "locate tpv102/mesh/tpv102_200m.geo");
   if (tpv104_path.empty() || tpv102_path.empty()) { return; }

   // R3-005: align both files on the `"lc = "` sentinel instead of a
   // hardcoded line-number cutoff (which was brittle against
   // header-length drift between TPV102 and TPV104).
   const std::string body_104 = TailFromCode(tpv104_path);
   const std::string body_102 = TailFromCode(tpv102_path);

   TEST_ASSERT(!body_104.empty(), "TPV104 geometry body non-empty");
   TEST_ASSERT(!body_102.empty(), "TPV102 geometry body non-empty");
   TEST_ASSERT(body_104 == body_102,
               "lines 30+ byte-identical TPV104 ↔ TPV102");

   if (body_104 != body_102)
   {
      // Print the first divergent line for diagnostics.
      size_t n = std::min(body_104.size(), body_102.size());
      for (size_t i = 0; i < n; ++i)
      {
         if (body_104[i] != body_102[i])
         {
            std::cerr << "  first diff at byte " << i << "\n";
            std::cerr << "    104: "
                      << body_104.substr(i, 60) << "\n"
                      << "    102: "
                      << body_102.substr(i, 60) << "\n";
            break;
         }
      }
   }
}

// ----------------------------------------------------------------------------
// T_TPV104_MESH_0b — 1000m variant has the same geometry body.
// ----------------------------------------------------------------------------
void TestGeometryByteIdentical1000m()
{
   std::cout << "\n[T_TPV104_MESH_0b] 1000m byte-identical vs TPV102\n";

   const std::string tpv104_path = locate("tpv104/mesh/tpv104_1000m.geo");
   const std::string tpv102_path = locate("tpv102/mesh/tpv102_1000m.geo");

   TEST_ASSERT(!tpv104_path.empty(), "locate tpv104/mesh/tpv104_1000m.geo");
   TEST_ASSERT(!tpv102_path.empty(), "locate tpv102/mesh/tpv102_1000m.geo");
   if (tpv104_path.empty() || tpv102_path.empty()) { return; }

   // R3-005: same sentinel-anchor approach — no hardcoded cutoff.
   TEST_ASSERT(TailFromCode(tpv104_path) == TailFromCode(tpv102_path),
               "1000m geometry bodies byte-identical TPV104 ↔ TPV102 "
               "(lc-anchored, CRLF-tolerant)");
}

// ----------------------------------------------------------------------------
// T_TPV104_MESH_0c — 500m variant exists and has the expected lc_fault
// value (= 500).  The 500m mesh is derived from the 200m template per
// plan §4.10 Step 10; it is NOT expected to byte-match TPV102 (there is
// no TPV102 500m counterpart).
// ----------------------------------------------------------------------------
void TestMesh500mPresentAndScaled()
{
   std::cout << "\n[T_TPV104_MESH_0c] 500m variant exists with lc_fault=500\n";

   const std::string path = locate("tpv104/mesh/tpv104_500m.geo");
   TEST_ASSERT(!path.empty(), "locate tpv104/mesh/tpv104_500m.geo");
   if (path.empty()) { return; }

   const std::string body = slurp(path);

   // R3-006 (review round 3): tolerate whitespace variation around
   // `=` and `;`.  Gmsh treats  `lc_fault=500;` and `lc_fault  = 500 ;`
   // identically, so a literal-string match is unnecessarily brittle.
   const size_t name_pos = body.find("lc_fault");
   TEST_ASSERT(name_pos != std::string::npos,
               "500m variant declares lc_fault");
   const size_t eq_pos = body.find('=', name_pos);
   const size_t semi_pos = (eq_pos != std::string::npos)
                           ? body.find(';', eq_pos)
                           : std::string::npos;
   TEST_ASSERT(eq_pos != std::string::npos
               && semi_pos != std::string::npos,
               "`lc_fault = <value>;` assignment found");
   std::string val = body.substr(eq_pos + 1, semi_pos - eq_pos - 1);
   val.erase(std::remove_if(val.begin(), val.end(),
                            [](unsigned char c) { return std::isspace(c); }),
             val.end());
   TEST_ASSERT(val == "500",
               "500m variant sets lc_fault to 500 (whitespace-tolerant)");

   TEST_ASSERT(body.find("TPV104") != std::string::npos,
               "500m header references TPV104");
}

// ----------------------------------------------------------------------------
// T_TPV104_MESH_1 — vertex/element counts.  DEFERRED: requires gmsh
// build to produce .msh files; gmsh invocation runs in the pythonenv
// conda env, not in the C++ unit-test harness.  The gate is tracked in
// the Phase-2 acceptance matrix under Step 10.
// ----------------------------------------------------------------------------
void TestMeshCountsDeferred()
{
   std::cout << "\n[T_TPV104_MESH_1] vertex/element counts — DEFERRED "
             << "(requires gmsh build; tracked by Phase-2 acceptance "
             << "matrix)\n";
}

// ----------------------------------------------------------------------------
// T_TPV104_MESH_2 — wave-speed arrival check.  DEFERRED until Step 9
// driver is available; the test requires running the driver with a
// delta-pulse source and measuring P-wave arrival at a station.
// ----------------------------------------------------------------------------
void TestWaveArrivalDeferred()
{
   std::cout << "\n[T_TPV104_MESH_2] wave-speed arrival — DEFERRED "
             << "(requires Step 9 driver)\n";
}

// ----------------------------------------------------------------------------
// T_TPV104_MESH_3 — physical-group parity.  DEFERRED — same reason.
// ----------------------------------------------------------------------------
void TestPhysicalGroupParityDeferred()
{
   std::cout << "\n[T_TPV104_MESH_3] physical-group parity — DEFERRED "
             << "(requires Step 9 driver)\n";
}

int main(int argc, char *argv[])
{
   TestGeometryByteIdentical();
   TestGeometryByteIdentical1000m();
   TestMesh500mPresentAndScaled();
   TestMeshCountsDeferred();
   TestWaveArrivalDeferred();
   TestPhysicalGroupParityDeferred();

   TEST_PRINT_RESULTS();
   return (num_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

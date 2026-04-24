// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Unit tests for SEAS_DIAG_TPV104_STATE probe instrumentation
// (§4.10 Step 11 gates T_TPV104_PROBE_1..3).

#include "test_macros.hpp"

#include <algorithm>
#include <cstdlib>
#include <dirent.h>
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

// Count occurrences of a needle in a haystack.
static int count_occurrences(const std::string &hay, const std::string &needle)
{
   int n = 0;
   size_t pos = 0;
   while ((pos = hay.find(needle, pos)) != std::string::npos)
   {
      ++n;
      pos += needle.size();
   }
   return n;
}

// ----------------------------------------------------------------------------
// T_TPV104_PROBE_1 — probe output header format.
// Every probe emission in the iterator must be guarded by
// `#ifdef SEAS_DIAG_TPV104_STATE ... #endif`, use `GetProbeFile(...)`,
// and write lines prefixed by `# probe=` / `# t qp_id ...`.
// ----------------------------------------------------------------------------
void TestProbeFormat()
{
   std::cout << "\n[T_TPV104_PROBE_1] probe file header format\n";

   const std::string path = locate("dynamic/tpv104_substep_iterator.cpp");
   TEST_ASSERT(!path.empty(), "locate iterator source");
   if (path.empty()) { return; }

   const std::string body = slurp(path);

   // Header line prefix `# probe=<name>  code=MFEM  rank=<r> ...` must
   // be emitted by `GetProbeFile` exactly once per file open.
   TEST_ASSERT(body.find("# probe=") != std::string::npos,
               "iterator emits `# probe=<name>` header per §5.1");
   TEST_ASSERT(body.find("code=MFEM") != std::string::npos,
               "header labels the emitter as MFEM");

   // Every Probe-N block must be inside `#ifdef SEAS_DIAG_TPV104_STATE`.
   const int n_ifdef = count_occurrences(body, "#ifdef SEAS_DIAG_TPV104_STATE");
   const int n_endif = count_occurrences(body, "#endif  // SEAS_DIAG_TPV104_STATE");
   TEST_ASSERT(n_ifdef >= 1,
               "iterator has ≥ 1 SEAS_DIAG_TPV104_STATE guard");
   TEST_ASSERT(n_ifdef == n_endif,
               "every SEAS_DIAG_TPV104_STATE #ifdef has a matching #endif");

   // Every named probe referenced by the tools (Probes 1, 2, 3, 4, 5)
   // must be emitted by the iterator under some code path.  The probe
   // names match the pytest fixtures in tpv104/scripts/tests/.
   const std::vector<std::string> required_probes = {
      "trial_traction",   // Probe 1
      "state_evolution",  // Probe 2
      "friction_coeff",   // Probe 3
      "slip_rate",        // Probe 4
      "corrected_imposed" // Probe 5
   };
   for (const auto &name : required_probes)
   {
      const std::string key = "GetProbeFile(\"" + name + "\")";
      TEST_ASSERT(body.find(key) != std::string::npos,
                  ("iterator emits probe '" + name + "'").c_str());
   }
}

// ----------------------------------------------------------------------------
// T_TPV104_PROBE_2 — disabled build emits no probe output.
// Sanity: the iterator source does NOT reference any probe output
// facility outside the SEAS_DIAG_TPV104_STATE guards (grep asserts
// this by counting total GetProbeFile refs vs inside-guard refs).
//
// An alternative runtime check (build a no-diag iterator, run it,
// confirm no tpv104_probe_*.txt files appear) requires running the
// driver — this static check avoids that dependency.
// ----------------------------------------------------------------------------
void TestDisabledBuildNoProbe()
{
   std::cout << "\n[T_TPV104_PROBE_2] disabled build emits no probe output\n";

   const std::string path = locate("dynamic/tpv104_substep_iterator.cpp");
   TEST_ASSERT(!path.empty(), "locate iterator source");
   if (path.empty()) { return; }

   const std::string body = slurp(path);

   // Walk the body and count `GetProbeFile(...)` references that occur
   // OUTSIDE any `#ifdef SEAS_DIAG_TPV104_STATE` ... `#endif` region.
   // (Function-definition site is OK because the function body itself
   // is wrapped in the guard at a higher level.)
   int depth = 0;
   int unguarded_hits = 0;
   const std::string open_tag  = "#ifdef SEAS_DIAG_TPV104_STATE";
   const std::string close_tag = "#endif  // SEAS_DIAG_TPV104_STATE";
   size_t pos = 0;
   while (pos < body.size())
   {
      const size_t next_open  = body.find(open_tag,  pos);
      const size_t next_close = body.find(close_tag, pos);
      const size_t next_probe = body.find("GetProbeFile(", pos);
      const size_t next = std::min({next_open, next_close, next_probe});
      if (next == std::string::npos) { break; }
      if (next == next_open)  { ++depth; pos = next + open_tag.size(); }
      else if (next == next_close) { --depth; pos = next + close_tag.size(); }
      else  // next_probe
      {
         if (depth == 0)
         {
            // Skip the function-definition line — its `GetProbeFile`
            // occurrence is inside an `#ifdef SEAS_DIAG_TPV104_STATE`
            // block that wraps the whole helper; the depth counter
            // does handle this correctly.  If we observe this at
            // depth 0, it is a real leak.
            ++unguarded_hits;
         }
         pos = next_probe + std::string("GetProbeFile(").size();
      }
   }
   TEST_ASSERT(unguarded_hits == 0,
               "no unguarded GetProbeFile() call (disabled build compiles clean)");
}

// ----------------------------------------------------------------------------
// T_TPV104_PROBE_3 — SEAS_DIAG_TPV104_STATE does not affect BP5.
// Static grep: no BP5 source file references the flag.
// ----------------------------------------------------------------------------
void TestBP5UnaffectedByDiag()
{
   std::cout << "\n[T_TPV104_PROBE_3] BP5 unaffected by SEAS_DIAG_TPV104_STATE\n";

   // Enumerate files under miniapps/seas/bp5/ and
   // miniapps/seas/drivers/seas_bp1* / seas_bp2* / seas_bp5*.
   const std::vector<std::string> bp5_source_roots = {
      locate("bp5"),                         // bp5/ directory (config, etc.)
      locate("drivers"),                     // driver directory (bp1/2/5 drivers live here)
      locate("solver"),                      // SEASOperator used by BP5
      locate("domain"),                      // elasticity / BCs
      locate("fault"),                       // fault basis
      locate("friction/dieterich_ruina.hpp") // extreme-care friction solver
   };

   int total_grep_hits = 0;
   int dirs_scanned = 0;
   for (const auto &root : bp5_source_roots)
   {
      if (root.empty()) { continue; }
      // Treat root either as a file or a directory.  Files: grep the
      // content directly.  Directories: listdir, then grep each entry
      // recursively one level down (sufficient for our layout).
      std::ifstream f_test(root);
      if (!f_test.is_open()) {/* maybe a directory */}

      DIR *d = opendir(root.c_str());
      if (d == nullptr)
      {
         // Treat as a single file.
         const std::string body = slurp(root);
         total_grep_hits += count_occurrences(body,
                                              "SEAS_DIAG_TPV104_STATE");
         continue;
      }
      ++dirs_scanned;
      struct dirent *entry;
      while ((entry = readdir(d)) != nullptr)
      {
         const std::string name(entry->d_name);
         if (name == "." || name == "..") { continue; }
         // Only scan text-ish source files (defensive against binary
         // files, .msh files with huge bodies, etc.).
         const bool is_source =
            name.size() > 4 && (
               name.substr(name.size() - 4) == ".cpp" ||
               name.substr(name.size() - 4) == ".hpp" ||
               name.substr(name.size() - 2) == ".h"
            );
         if (!is_source) { continue; }
         // R5-006 follow-up: `drivers/tpv104_driver.cpp` legitimately
         // calls `Tpv104SubStepIterator::CloseAllProbeFiles()` and has
         // a comment mentioning `SEAS_DIAG_TPV104_STATE`.  The intent
         // of this grep is "BP5 / shared drivers don't reference it" —
         // so skip any file whose name starts with `tpv104_`.
         if (name.rfind("tpv104_", 0) == 0) { continue; }
         const std::string path = root + "/" + name;
         const std::string body = slurp(path);
         const int n = count_occurrences(body, "SEAS_DIAG_TPV104_STATE");
         total_grep_hits += n;
         if (n > 0)
         {
            std::cerr << "  FAIL: " << path << " references "
                      << "SEAS_DIAG_TPV104_STATE (" << n << " times)\n";
         }
      }
      closedir(d);
   }
   TEST_ASSERT(total_grep_hits == 0,
               "BP5 / drivers / solver / domain / fault / dieterich_ruina "
               "sources do NOT reference SEAS_DIAG_TPV104_STATE");
   std::cout << "  directories scanned: " << dirs_scanned << "\n";
}

int main(int argc, char *argv[])
{
   TestProbeFormat();
   TestDisabledBuildNoProbe();
   TestBP5UnaffectedByDiag();

   TEST_PRINT_RESULTS();
   return (num_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

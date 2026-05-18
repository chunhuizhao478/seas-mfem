// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// Phase 1 of miniapps/seas/io/PLAN_paraview_compaction_2026-04-28.md.
//
// Bit-exactness test (modulo ASCII setprecision(10) truncation):
//
//   1. Build a tiny inline tetrahedral mesh with 2 fault faces.
//   2. Run `WriteFaultSurfaceVTU` twice on identical input data:
//        (a) `SetLegacyAsciiVTU(true)`  -> ASCII per-rank VTU
//        (b) default                    -> single binary VTU on rank 0
//   3. Parse both VTUs' CellData via `vtu::ParseFaultVTUCellData`,
//      which transparently handles ASCII and base64-binary formats.
//   4. Assert per-cell max-abs-diff < 1e-9 across all 12 standard fields.
//
// Also verifies file count: the binary path must produce exactly
//   FaultSurface/{fault_surface_c0.vtu, fault_surface.pvd}
// — no per-rank files, no PVTU.
//
// Phase 1 acceptance criterion: "the new VTU is byte-readable by
// ParaView 5.10 and renders identical geometry/fields to the old
// ASCII output (visual diff, plus per-field max-abs-diff < 1e-9
// verified by `test_fault_surface_vtu_binary.cpp`)."

#include "mfem.hpp"
#include "../../io/paraview_output.hpp"
#include "../../io/fault_vtu_binary.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>      // ::rmdir
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

#define TEST_NEAR(val, exp, tol, msg) do { \
   num_tests++; \
   double v_ = (val), e_ = (exp), t_ = (tol); \
   if (std::abs(v_ - e_) <= t_) { num_passed++; \
      std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
      << " (got " << v_ << ", expected " << e_ \
      << ", |diff|=" << std::abs(v_-e_) << ", tol=" << t_ << ")\n"; } \
} while (0)

// ---------------------------------------------------------------------------
// Build a small mesh with 1 interior fault face (one tetrahedral split into
// two halves separated by a triangle interior face).  Mirrors the structure
// used by the other fault-surface tests.
// ---------------------------------------------------------------------------
static Mesh BuildTwoTetMesh()
{
   // Two tets sharing a triangular interior face at z=0.
   //   A = (0,0,-1), B = (1,0, 0), C = (0,1, 0), D = (0,0, 0): tet 1
   //   B, C, D, E = (0,0, 1):                                   tet 2
   // Interior face = triangle (B, C, D).
   Mesh m(/*dim=*/3, /*nv=*/5, /*nelem=*/2, /*nbdr=*/0, /*sdim=*/3);
   m.AddVertex(0., 0., -1.);  // 0 A
   m.AddVertex(1., 0.,  0.);  // 1 B
   m.AddVertex(0., 1.,  0.);  // 2 C
   m.AddVertex(0., 0.,  0.);  // 3 D
   m.AddVertex(0., 0.,  1.);  // 4 E
   const int t1[4] = {0,1,2,3};
   const int t2[4] = {1,2,3,4};
   m.AddTet(t1, /*attr=*/1);
   m.AddTet(t2, /*attr=*/2);
   m.FinalizeTetMesh(/*generate_edges=*/0, /*refine=*/0, /*fix_orient=*/true);
   return m;
}

static int CountFilesIn(const std::string &dir)
{
   DIR *d = opendir(dir.c_str());
   if (!d) { return -1; }
   int n = 0;
   struct dirent *e;
   while ((e = readdir(d)) != nullptr)
   {
      const std::string name(e->d_name);
      if (name == "." || name == "..") { continue; }
      ++n;
   }
   closedir(d);
   return n;
}

static void RemoveDir(const std::string &dir)
{
   DIR *d = opendir(dir.c_str());
   if (!d) { return; }
   struct dirent *e;
   while ((e = readdir(d)) != nullptr)
   {
      const std::string name(e->d_name);
      if (name == "." || name == "..") { continue; }
      const std::string path = dir + "/" + name;
      ::remove(path.c_str());
   }
   closedir(d);
   ::rmdir(dir.c_str());
}

int main(int argc, char *argv[])
{
   (void)argc; (void)argv;
   std::cout << "=== test_fault_surface_vtu_binary ===\n";

   Mesh mesh = BuildTwoTetMesh();
   const int nbf = 1;
   const int n_dofs = 1;        // 1 fault face × 1 nbf

   Array<int> fault_faces, empty_shared;
   for (int f = 0; f < mesh.GetNumFaces(); ++f)
   {
      if (mesh.GetFaceInformation(f).IsInterior())
      { fault_faces.Append(f); break; }
   }
   TEST_ASSERT(fault_faces.Size() == 1, "found 1 interior face");

   // Synthesize per-DOF data.  Use identifiable values so cross-format
   // comparison is unambiguous.
   Vector slip(2 * n_dofs);          slip = 0.123456789012345;
   Vector slip_rate(2 * n_dofs);     slip_rate = 1.234567890123456e-3;
   Vector traction(2 * n_dofs);      traction = -1.5e7;
   Vector state(n_dofs);             state = 0.456789012345678;
   Vector normal_stress(n_dofs);     normal_stress = 5.0e7;
   Vector a(n_dofs);                 a = 0.0123;
   Vector Dc(n_dofs);                Dc = 0.14;
   Vector x2(n_dofs);                x2 = 1.0;
   Vector x3(n_dofs);                x3 = -3.0;

   // -- Path A: legacy ASCII per-rank writer ---------------------------------
   const std::string prefixA = "/tmp/test_fault_vtu_binary_ASCII";
   RemoveDir(prefixA + "/FaultSurface");
   ::mkdir(prefixA.c_str(), 0755);
   {
      ParaViewOutput<Mesh> pv(prefixA, mesh, /*order=*/1);
      pv.InitFaultOutputBP5(fault_faces, empty_shared, nbf);
      pv.SetLegacyAsciiVTU(true);
      pv.WriteFaultSurfaceVTU(prefixA, /*cycle=*/0, /*time=*/0.0,
                              /*rank=*/0, /*nranks=*/1,
                              slip, slip_rate, traction, state, normal_stress,
                              a, Dc, x2, x3);
   }
   const std::string vtuA = prefixA + "/FaultSurface/fault_surface_r0_c0.vtu";
   TEST_ASSERT(std::ifstream(vtuA).good(), "ASCII path emitted r0_c0.vtu");

   // -- Path B: Phase 1 single binary VTU on rank 0 --------------------------
   const std::string prefixB = "/tmp/test_fault_vtu_binary_BIN";
   RemoveDir(prefixB + "/FaultSurface");
   ::mkdir(prefixB.c_str(), 0755);
   {
      ParaViewOutput<Mesh> pv(prefixB, mesh, /*order=*/1);
      pv.InitFaultOutputBP5(fault_faces, empty_shared, nbf);
      // R-002 (REVIEW.md 2026-04-28): on HDF5-enabled builds the
      // default is now Hdf5; this test asserts the binary VTU file
      // layout, so force the Vtu backend.  legacy_ascii_ default =
      // false; emit binary.
      pv.SetFaultOutputMode(
         ParaViewOutput<Mesh>::FaultOutputMode::Vtu);
      pv.WriteFaultSurfaceVTU(prefixB, /*cycle=*/0, /*time=*/0.0,
                              /*rank=*/0, /*nranks=*/1,
                              slip, slip_rate, traction, state, normal_stress,
                              a, Dc, x2, x3);
   }
   const std::string vtuB = prefixB + "/FaultSurface/fault_surface_c0.vtu";
   TEST_ASSERT(std::ifstream(vtuB).good(), "BIN path emitted c0.vtu");

   // File-count check: binary path produces exactly 2 files (the VTU and the
   // PVD).  Legacy ASCII path produces 3 files (r0_c0.vtu, c0.pvtu, pvd).
   const int n_bin = CountFilesIn(prefixB + "/FaultSurface");
   const int n_asc = CountFilesIn(prefixA + "/FaultSurface");
   TEST_ASSERT(n_bin == 2,
               "binary writer file count == 2 (got " + std::to_string(n_bin) + ")");
   TEST_ASSERT(n_asc == 3,
               "legacy ASCII writer file count == 3 (got " +
               std::to_string(n_asc) + ")");

   // -- Compare CellData per-field --------------------------------------------
   const std::vector<std::string> fields = {
      "slip_dip", "slip_strike", "slip_rate_dip", "slip_rate_strike",
      "traction_dip", "traction_strike", "state_variable", "normal_stress",
      "param_a", "param_Dc", "fault_x2", "fault_x3"
   };
   for (const auto &fname : fields)
   {
      auto vA = vtu::ParseFaultVTUCellData(vtuA, fname);
      auto vB = vtu::ParseFaultVTUCellData(vtuB, fname);
      TEST_ASSERT(vA.size() == vB.size() && !vA.empty(),
                  "field " + fname + ": both writers emit same array length");
      double max_abs = 0.0;
      const std::size_t n = std::min(vA.size(), vB.size());
      for (std::size_t i = 0; i < n; ++i)
      {
         max_abs = std::max(max_abs, std::abs(vA[i] - vB[i]));
      }
      TEST_NEAR(max_abs, 0.0, 1e-9,
                "field " + fname + " ASCII vs binary max-abs-diff");
   }

   // Cleanup
   RemoveDir(prefixA + "/FaultSurface");
   ::rmdir(prefixA.c_str());
   RemoveDir(prefixB + "/FaultSurface");
   ::rmdir(prefixB.c_str());

   std::cout << "\n=== Summary: " << num_passed << " / " << num_tests
             << " passed; " << num_failed << " failed ===\n";
   return num_failed > 0 ? 1 : 0;
}

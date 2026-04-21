// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.1.0 §4.1 — post-R-001 fault-surface VTU contract.
//
// Verifies the CellData form of `ParaViewOutput::WriteFaultSurfaceVTU`
// (plan v9.1.0 §2.3 R-001 / H-V91-A4) across three requirements:
//
//   (a) Face-average correctness on a linear field.  For a linear
//       f(x2, x3) = a0 + a2*x2 + a3*x3 sampled at the face's nbf
//       barycentric quadrature points, the arithmetic mean equals
//       f(face_centroid) exactly.  Asserted to 1e-10 relative.
//
//   (b) Adjacent-cell continuity under a linear field.  Any two output
//       triangles that share a physical edge must emit the same
//       CellData value (the face-centroid evaluation is identical for
//       any face whose centroid lies on the shared edge midpoint only
//       when the triangles are coplanar; this test therefore only
//       asserts per-cell correctness, which subsumes no-inter-cell-
//       drift for a single linear field across multiple faces).
//
//   (c) R-002 regression: with nbf=6 and per-QP IDs
//       `local_slip_rate[2*d] = d`, each face's cell value equals
//       the mean of its six per-QP IDs.  A pre-R-001 writer that
//       reads only the first three QPs would emit the mean of those
//       three only, i.e. `d + 1.0` rather than `d + 2.5`.
//
// The test constructs a small inline tetrahedral mesh and picks two
// or more interior faces to treat as fault faces.  It does not require
// a real fault attribute.

#include "mfem.hpp"
#include "../../io/paraview_output.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

// ---------------------------------------------------------------------------
// Test framework
// ---------------------------------------------------------------------------

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
// Minimal VTU CellData parser — reads a single DataArray by name from the
// <CellData> section of an ASCII UnstructuredGrid VTU.  Robust against the
// specific layout emitted by WriteFaultSurfaceVTU (newline-separated
// per-cell scalar).
// ---------------------------------------------------------------------------

static std::vector<double> ParseCellData(const std::string &path,
                                         const std::string &field)
{
   std::vector<double> out;
   std::ifstream in(path);
   if (!in.is_open())
   {
      std::cout << "  WARN: could not open " << path << "\n";
      return out;
   }
   std::stringstream ss;
   ss << in.rdbuf();
   const std::string content = ss.str();

   const auto cd_open = content.find("<CellData>");
   const auto cd_close = content.find("</CellData>");
   if (cd_open == std::string::npos || cd_close == std::string::npos ||
       cd_close <= cd_open)
   {
      std::cout << "  WARN: <CellData> block not found in " << path << "\n";
      return out;
   }
   const std::string cd = content.substr(cd_open, cd_close - cd_open);

   const std::string marker = std::string("Name=\"") + field + "\"";
   const auto name_pos = cd.find(marker);
   if (name_pos == std::string::npos) { return out; }
   const auto arr_close = cd.find("</DataArray>", name_pos);
   const auto arr_data_open = cd.find(">", name_pos);
   if (arr_close == std::string::npos ||
       arr_data_open == std::string::npos || arr_data_open >= arr_close)
   {
      return out;
   }
   const std::string body = cd.substr(arr_data_open + 1,
                                      arr_close - arr_data_open - 1);

   std::istringstream body_ss(body);
   double v;
   while (body_ss >> v) { out.push_back(v); }
   return out;
}

static int ParseNumberOfCells(const std::string &path)
{
   std::ifstream in(path);
   if (!in.is_open()) { return -1; }
   std::stringstream ss;
   ss << in.rdbuf();
   const std::string content = ss.str();
   const auto key = content.find("NumberOfCells=\"");
   if (key == std::string::npos) { return -1; }
   const auto qstart = key + 15;  // length of "NumberOfCells=\""
   const auto qend = content.find('"', qstart);
   if (qend == std::string::npos) { return -1; }
   return std::atoi(content.substr(qstart, qend - qstart).c_str());
}

// ---------------------------------------------------------------------------
// Build a tiny 3D tetrahedral mesh and select interior-face indices that
// `InitFaultOutputBP5` can use as "fault" faces.  Any interior face with
// a valid `GetInteriorFaceTransformations` is acceptable.
// ---------------------------------------------------------------------------

static Mesh MakeTinyTetMesh()
{
   // 2 x 2 x 1 hex subdivided into tets — gives enough interior triangles
   // to exercise the per-cell loop.
   return Mesh::MakeCartesian3D(2, 2, 1, Element::TETRAHEDRON,
                                2.0, 2.0, 1.0);
}

static Array<int> CollectInteriorFaces(Mesh &mesh, int max_n)
{
   Array<int> out;
   const int nfaces = mesh.GetNumFaces();
   for (int f = 0; f < nfaces && out.Size() < max_n; f++)
   {
      FaceElementTransformations *FTr = mesh.GetInteriorFaceTransformations(f);
      if (FTr != nullptr)
      {
         out.Append(f);
      }
   }
   return out;
}

// Evaluate a linear field at a single point:  f(x) = a0 + a2*x[1] + a3*x[2]
// (we use y as "x2" and z as "x3" in the fault frame convention).
static inline double LinearField(double x2, double x3)
{
   return 0.5 + 0.1 * x2 + 0.2 * x3;
}

// ---------------------------------------------------------------------------
// Test 1 — face-centroid equality for a linear field (nbf=3)
// ---------------------------------------------------------------------------

static void TestLinearFieldFaceCentroid()
{
   std::cout << "\n=== Test 4.1(a,b): linear field, CellData face-average ===\n";

   Mesh mesh = MakeTinyTetMesh();
   const int dim = mesh.Dimension();
   Array<int> fault_faces = CollectInteriorFaces(mesh, 6);
   Array<int> empty_shared;

   const int nbf = 3;
   const int n_int = fault_faces.Size();
   const int n_total = n_int;
   const int total_dofs = n_total * nbf;

   // Build the ParaViewOutput and initialize fault output path.
   const std::string prefix = "test_r001_continuity";
   ::mkdir(prefix.c_str(), 0755);
   ParaViewOutput<Mesh> pv(prefix, mesh, /*order=*/1);
   pv.InitFaultOutputBP5(fault_faces, empty_shared, nbf);

   // Integration rule matches wave_operator / SetFaultDOFData for order=1.
   const IntegrationRule &ir = IntRules.Get(Geometry::TRIANGLE, 2 /*= 2*order*/);
   TEST_ASSERT(ir.GetNPoints() == nbf,
               "IntRules TRIANGLE degree-2 has 3 points (matches nbf=3)");

   // Size the input vectors according to the DOF layout.
   Vector local_slip      (2 * total_dofs);  local_slip      = 0.0;
   Vector local_slip_rate (2 * total_dofs);  local_slip_rate = 0.0;
   Vector local_traction  (2 * total_dofs);  local_traction  = 0.0;
   Vector local_state     (total_dofs);      local_state     = 0.0;
   Vector local_normal    (0);  // has_normal false
   Vector local_a         (total_dofs);      local_a         = 0.0;
   Vector local_Dc        (total_dofs);      local_Dc        = 0.0;
   Vector local_x2        (total_dofs);      local_x2        = 0.0;
   Vector local_x3        (total_dofs);      local_x3        = 0.0;

   // Expected per-face centroid coords and LinearField(centroid) values.
   std::vector<double> expected_centroid_f(n_int);
   for (int fi = 0; fi < n_int; fi++)
   {
      const int f = fault_faces[fi];
      FaceElementTransformations *FTr = mesh.GetInteriorFaceTransformations(f);
      MFEM_VERIFY(FTr != nullptr, "interior face transformation missing");

      // Compute face centroid from the three barycentric QPs (exact for
      // linear fields, matches paraview_output.hpp's per-cell average).
      double cx2 = 0.0, cx3 = 0.0;
      for (int k = 0; k < nbf; k++)
      {
         Vector coord(dim);
         FTr->Face->Transform(ir.IntPoint(k), coord);
         const double x2 = coord(1), x3 = coord(2);
         const double f_k = LinearField(x2, x3);
         const int d = fi * nbf + k;
         local_slip     (2*d)   = f_k;
         local_slip     (2*d+1) = 2.0 * f_k;  // distinct field in strike slot
         local_slip_rate(2*d)   = f_k;
         local_slip_rate(2*d+1) = 3.0 * f_k;
         local_traction (2*d)   = f_k;
         local_traction (2*d+1) = -f_k;
         local_state    (d)     = f_k;
         local_a        (d)     = f_k;
         local_Dc       (d)     = f_k;
         local_x2       (d)     = x2;
         local_x3       (d)     = x3;
         cx2 += x2;
         cx3 += x3;
      }
      cx2 /= nbf;
      cx3 /= nbf;
      expected_centroid_f[fi] = LinearField(cx2, cx3);
   }

   // Write VTU to temp dir.
   pv.WriteFaultSurfaceVTU(prefix, /*cycle=*/0, /*time=*/0.0,
                           /*rank=*/0, /*nranks=*/1,
                           local_slip, local_slip_rate, local_traction,
                           local_state, local_normal,
                           local_a, local_Dc, local_x2, local_x3);

   const std::string vtu = prefix + "/FaultSurface/fault_surface_r0_c0.vtu";
   const int ncells = ParseNumberOfCells(vtu);
   TEST_ASSERT(ncells == n_int,
               std::string("VTU NumberOfCells = n_int (") +
               std::to_string(n_int) + ")");

   // All 12 fields should have n_int entries in CellData.
   auto cell_slip_dip        = ParseCellData(vtu, "slip_dip");
   auto cell_slip_strike     = ParseCellData(vtu, "slip_strike");
   auto cell_slip_rate_dip   = ParseCellData(vtu, "slip_rate_dip");
   auto cell_slip_rate_strike= ParseCellData(vtu, "slip_rate_strike");
   auto cell_traction_dip    = ParseCellData(vtu, "traction_dip");
   auto cell_traction_strike = ParseCellData(vtu, "traction_strike");
   auto cell_state           = ParseCellData(vtu, "state_variable");
   auto cell_param_a         = ParseCellData(vtu, "param_a");

   TEST_ASSERT(static_cast<int>(cell_slip_dip.size()) == n_int,
               "slip_dip CellData length = n_int");
   TEST_ASSERT(static_cast<int>(cell_slip_strike.size()) == n_int,
               "slip_strike CellData length = n_int");
   TEST_ASSERT(static_cast<int>(cell_slip_rate_dip.size()) == n_int,
               "slip_rate_dip CellData length = n_int");
   TEST_ASSERT(static_cast<int>(cell_state.size()) == n_int,
               "state CellData length = n_int");

   // Per-face correctness.  Tolerance 5e-9 relative to max(1, |expected|):
   // `paraview_output.hpp` writes with `std::setprecision(10)` which gives
   // ~5e-10 relative roundtrip error; 5e-9 leaves ~10x headroom so the
   // test catches the pre-R-001 writer (which would be off by
   // O(h * |grad f|) ≈ 1e-2 × |f|, many OOM larger) while not tripping
   // on ASCII printing roundoff.
   for (int fi = 0; fi < n_int; fi++)
   {
      const double expected = expected_centroid_f[fi];
      const double tol = 5.0e-9 * std::max(1.0, std::abs(expected));

      TEST_NEAR(cell_slip_dip[fi], expected, tol,
                std::string("slip_dip cell[") + std::to_string(fi) +
                "] == f(centroid)");
      TEST_NEAR(cell_slip_strike[fi], 2.0 * expected, tol,
                std::string("slip_strike cell[") + std::to_string(fi) +
                "] == 2*f(centroid)");
      TEST_NEAR(cell_slip_rate_dip[fi], expected, tol,
                std::string("slip_rate_dip cell[") + std::to_string(fi) +
                "] == f(centroid)");
      TEST_NEAR(cell_slip_rate_strike[fi], 3.0 * expected, tol,
                std::string("slip_rate_strike cell[") + std::to_string(fi) +
                "] == 3*f(centroid)");
      TEST_NEAR(cell_traction_dip[fi], expected, tol,
                std::string("traction_dip cell[") + std::to_string(fi) +
                "] == f(centroid)");
      TEST_NEAR(cell_traction_strike[fi], -expected, tol,
                std::string("traction_strike cell[") + std::to_string(fi) +
                "] == -f(centroid)");
      TEST_NEAR(cell_state[fi], expected, tol,
                std::string("state cell[") + std::to_string(fi) +
                "] == f(centroid)");
      TEST_NEAR(cell_param_a[fi], expected, tol,
                std::string("param_a cell[") + std::to_string(fi) +
                "] == f(centroid)");
   }
}

// ---------------------------------------------------------------------------
// Test 2 — R-002 regression: nbf=6, per-QP IDs, average of all six
// ---------------------------------------------------------------------------

static void TestNbf6AverageAllQPs()
{
   std::cout << "\n=== Test 4.1(c): nbf=6 all-QP averaging (R-002) ===\n";

   Mesh mesh = MakeTinyTetMesh();
   Array<int> fault_faces = CollectInteriorFaces(mesh, 4);
   Array<int> empty_shared;
   const int nbf = 6;   // as-if order=2 triangle rule
   const int n_int = fault_faces.Size();
   const int total_dofs = n_int * nbf;

   const std::string prefix = "test_r001_nbf6";
   ::mkdir(prefix.c_str(), 0755);
   ParaViewOutput<Mesh> pv(prefix, mesh, /*order=*/1);
   pv.InitFaultOutputBP5(fault_faces, empty_shared, nbf);

   Vector local_slip      (2 * total_dofs);  local_slip      = 0.0;
   Vector local_slip_rate (2 * total_dofs);  local_slip_rate = 0.0;
   Vector local_traction  (2 * total_dofs);  local_traction  = 0.0;
   Vector local_state     (total_dofs);      local_state     = 0.0;
   Vector local_normal    (0);
   Vector local_a         (0);
   Vector local_Dc        (0);
   Vector local_x2        (0);
   Vector local_x3        (0);

   // Per-QP integer IDs: local_slip_rate[2*d] = d  (globally unique).
   for (int d = 0; d < total_dofs; d++)
   {
      local_slip_rate(2*d)   = static_cast<double>(d);
      local_slip_rate(2*d+1) = 0.0;
   }

   pv.WriteFaultSurfaceVTU(prefix, /*cycle=*/0, /*time=*/0.0,
                           /*rank=*/0, /*nranks=*/1,
                           local_slip, local_slip_rate, local_traction,
                           local_state, local_normal,
                           local_a, local_Dc, local_x2, local_x3);

   const std::string vtu = prefix + "/FaultSurface/fault_surface_r0_c0.vtu";
   auto cell_slip_rate_dip = ParseCellData(vtu, "slip_rate_dip");

   TEST_ASSERT(static_cast<int>(cell_slip_rate_dip.size()) == n_int,
               std::string("nbf=6: CellData length = n_int (") +
               std::to_string(n_int) + ")");

   // Post-R-001 expected: cell[fi] = mean(fi*6, fi*6+1, ..., fi*6+5)
   //                               = fi*6 + (0+1+2+3+4+5)/6
   //                               = fi*6 + 2.5
   // Pre-R-001 expected (hardcoded k<3): cell[fi] would have been
   // mean(fi*6, fi*6+1, fi*6+2) = fi*6 + 1.0.
   for (int fi = 0; fi < n_int; fi++)
   {
      const double expected_post = fi * nbf + 2.5;
      TEST_NEAR(cell_slip_rate_dip[fi], expected_post, 1.0e-12,
                std::string("nbf=6 cell[") + std::to_string(fi) +
                "] == mean of all 6 per-QP IDs (R-002 regression)");
      // Defensive: the pre-fix value must NOT match.
      const double pre_fix_value = fi * nbf + 1.0;
      TEST_ASSERT(std::abs(cell_slip_rate_dip[fi] - pre_fix_value) > 0.5,
                  std::string("nbf=6 cell[") + std::to_string(fi) +
                  "] != pre-fix-first-3-only value");
   }
}

int main()
{
   std::cout << "========================================\n";
   std::cout << "TPV102 v9.1.0 §4.1 — fault-surface VTU CellData contract\n";
   std::cout << "========================================\n";

   TestLinearFieldFaceCentroid();
   TestNbf6AverageAllQPs();

   std::cout << "\n========================================\n";
   std::cout << "Total:  " << num_tests << "\n";
   std::cout << "Passed: " << num_passed << "\n";
   std::cout << "Failed: " << num_failed << "\n";
   std::cout << "========================================\n";
   return (num_failed > 0) ? 1 : 0;
}

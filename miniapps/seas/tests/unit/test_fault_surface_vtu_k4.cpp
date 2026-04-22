// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.2.0 §19 — R-V92-E02 stage-4 VTU diagnostic.
//
// Verifies that `WriteFaultSurfaceVTU` correctly handles the optional
// `_k4` stage-4 buffers:
//
//   (a) When all three `_k4` vectors are EMPTY (default), the output
//       VTU contains only the original 12 CellData fields — no
//       `<field>_k4` fields appear.  This is the backwards-compat path.
//
//   (b) When all three `_k4` vectors are POPULATED, the output VTU
//       contains 5 additional CellData fields:
//         slip_rate_dip_k4, slip_rate_strike_k4,
//         traction_dip_k4,  traction_strike_k4,
//         normal_stress_k4.
//       Each `_k4` field's per-cell value equals the face-averaged
//       `_k4` input value (using the same k<nbf averaging as the
//       corresponding non-k4 field).
//
//   (c) With k4 inputs deliberately DIFFERENT from the averaged
//       inputs, the emitted `_k4` CellData values reflect the k4
//       input — i.e. the diagnostic actually distinguishes the two
//       buffers.  This is what addresses R-V92-E02: if RK4-averaging
//       corrupts the end-of-pipeline DOFData, the `<field>_k4 -
//       <field>` delta is nonzero in the production 12 s run and
//       H-V92-K is confirmed.

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

// Parse a named CellData array from the emitted VTU.  Returns an empty
// vector if the field is absent.
static std::vector<double> ParseCellData(const std::string &path,
                                         const std::string &field)
{
   std::vector<double> out;
   std::ifstream in(path);
   if (!in.is_open()) { return out; }
   std::stringstream ss;
   ss << in.rdbuf();
   const std::string content = ss.str();
   const auto cd_open  = content.find("<CellData>");
   const auto cd_close = content.find("</CellData>");
   if (cd_open == std::string::npos || cd_close == std::string::npos) { return out; }
   const std::string cd = content.substr(cd_open, cd_close - cd_open);
   const std::string marker = std::string("Name=\"") + field + "\"";
   const auto name_pos = cd.find(marker);
   if (name_pos == std::string::npos) { return out; }
   const auto arr_close = cd.find("</DataArray>", name_pos);
   const auto arr_open  = cd.find(">", name_pos);
   if (arr_close == std::string::npos || arr_open == std::string::npos ||
       arr_open >= arr_close) { return out; }
   const std::string body = cd.substr(arr_open + 1, arr_close - arr_open - 1);
   std::istringstream body_ss(body);
   double v;
   while (body_ss >> v) { out.push_back(v); }
   return out;
}

static bool CellDataContainsField(const std::string &path,
                                  const std::string &field)
{
   std::ifstream in(path);
   if (!in.is_open()) { return false; }
   std::stringstream ss;
   ss << in.rdbuf();
   const std::string content = ss.str();
   const auto cd_open  = content.find("<CellData>");
   const auto cd_close = content.find("</CellData>");
   if (cd_open == std::string::npos || cd_close == std::string::npos)
   {
      return false;
   }
   const std::string cd = content.substr(cd_open, cd_close - cd_open);
   return cd.find(std::string("Name=\"") + field + "\"") != std::string::npos;
}

static Mesh MakeTinyTetMesh()
{
   return Mesh::MakeCartesian3D(2, 2, 1, Element::TETRAHEDRON,
                                2.0, 2.0, 1.0);
}

static Array<int> CollectInteriorFaces(Mesh &mesh, int max_n)
{
   Array<int> out;
   const int nfaces = mesh.GetNumFaces();
   for (int f = 0; f < nfaces && out.Size() < max_n; f++)
   {
      if (mesh.GetInteriorFaceTransformations(f) != nullptr) { out.Append(f); }
   }
   return out;
}

// ---------------------------------------------------------------------------
// Test (a) — backwards compat: k4 vectors EMPTY ⇒ no _k4 fields in VTU.
// ---------------------------------------------------------------------------
static void TestBackwardsCompatNoK4()
{
   std::cout << "\n=== R-V92-E02 (a): empty k4 ⇒ no _k4 fields ===\n";

   Mesh mesh = MakeTinyTetMesh();
   Array<int> fault_faces = CollectInteriorFaces(mesh, 4);
   Array<int> empty_shared;
   const int nbf = 3;
   const int n_int = fault_faces.Size();
   const int total_dofs = n_int * nbf;

   const std::string prefix = "test_r_v92_e02_no_k4";
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

   // Default constructed (empty) k4 vectors — default-argument path.
   pv.WriteFaultSurfaceVTU(prefix, /*cycle=*/0, /*time=*/0.0,
                           /*rank=*/0, /*nranks=*/1,
                           local_slip, local_slip_rate, local_traction,
                           local_state, local_normal,
                           local_a, local_Dc, local_x2, local_x3);

   const std::string vtu = prefix + "/FaultSurface/fault_surface_r0_c0.vtu";
   TEST_ASSERT(!CellDataContainsField(vtu, "normal_stress_k4"),
               "no normal_stress_k4 when k4 vectors are empty");
   TEST_ASSERT(!CellDataContainsField(vtu, "traction_dip_k4"),
               "no traction_dip_k4 when k4 vectors are empty");
   TEST_ASSERT(!CellDataContainsField(vtu, "traction_strike_k4"),
               "no traction_strike_k4 when k4 vectors are empty");
   TEST_ASSERT(!CellDataContainsField(vtu, "slip_rate_dip_k4"),
               "no slip_rate_dip_k4 when k4 vectors are empty");
   TEST_ASSERT(!CellDataContainsField(vtu, "slip_rate_strike_k4"),
               "no slip_rate_strike_k4 when k4 vectors are empty");
   TEST_ASSERT(CellDataContainsField(vtu, "normal_stress"),
               "original normal_stress field still present (backwards compat)");
}

// ---------------------------------------------------------------------------
// Test (b) + (c) — k4 vectors populated with values distinct from the
// averaged inputs; verify both sets of fields appear and carry the
// expected per-face-averaged values.
// ---------------------------------------------------------------------------
static void TestK4FieldsPopulatedAndDistinct()
{
   std::cout << "\n=== R-V92-E02 (b+c): populated k4 ⇒ _k4 fields carry distinct values ===\n";

   Mesh mesh = MakeTinyTetMesh();
   Array<int> fault_faces = CollectInteriorFaces(mesh, 4);
   Array<int> empty_shared;
   const int nbf = 3;
   const int n_int = fault_faces.Size();
   const int total_dofs = n_int * nbf;

   const std::string prefix = "test_r_v92_e02_with_k4";
   ::mkdir(prefix.c_str(), 0755);
   ParaViewOutput<Mesh> pv(prefix, mesh, /*order=*/1);
   pv.InitFaultOutputBP5(fault_faces, empty_shared, nbf);

   // Averaged inputs: constant.  sigma_n_avg = 120e6 (equilibrium).
   // V = 0 (pre-rupture).  traction_avg = 0.
   Vector local_slip      (2 * total_dofs);  local_slip      = 0.0;
   Vector local_slip_rate (2 * total_dofs);  local_slip_rate = 0.0;
   Vector local_traction  (2 * total_dofs);  local_traction  = 0.0;
   Vector local_state     (total_dofs);      local_state     = 0.0;
   Vector local_normal    (total_dofs);      local_normal    = 120e6;
   Vector local_a         (0);
   Vector local_Dc        (0);
   Vector local_x2        (0);
   Vector local_x3        (0);

   // k4 inputs: sigma_n_k4 varies per DOF in a triangle wave — the
   // face-averaged value on face fi should be the arithmetic mean of
   // the three per-DOF values = 70e6 + fi * 1e6 for face fi.  This
   // proves the per-face averaging loop handles the k4 buffers with
   // the same stride / bounds as the averaged buffers.
   Vector local_slip_rate_k4 (2 * total_dofs);  local_slip_rate_k4 = 0.0;
   Vector local_traction_k4  (2 * total_dofs);  local_traction_k4  = 0.0;
   Vector local_normal_k4    (total_dofs);      local_normal_k4    = 0.0;

   for (int fi = 0; fi < n_int; fi++)
   {
      const int base = fi * nbf;
      // Normal stress k4: {A-1, A, A+1}, mean = A.  Choose A = 70e6 + fi * 1e6.
      const double A = 70e6 + static_cast<double>(fi) * 1e6;
      local_normal_k4(base + 0) = A - 1.0e6;
      local_normal_k4(base + 1) = A;
      local_normal_k4(base + 2) = A + 1.0e6;
      // slip_rate_dip_k4: face-averaged value = fi (2*d index is slip_dip comp).
      // All three DOFs get the same value fi for simplicity.
      for (int k = 0; k < nbf; k++)
      {
         const int d = base + k;
         local_slip_rate_k4(2*d)     = static_cast<double>(fi);
         local_slip_rate_k4(2*d + 1) = 5.0 + static_cast<double>(fi);
         local_traction_k4 (2*d)     = 2e6 * static_cast<double>(fi + 1);
         local_traction_k4 (2*d + 1) = -3e6 * static_cast<double>(fi + 1);
      }
   }

   pv.WriteFaultSurfaceVTU(prefix, /*cycle=*/0, /*time=*/0.0,
                           /*rank=*/0, /*nranks=*/1,
                           local_slip, local_slip_rate, local_traction,
                           local_state, local_normal,
                           local_a, local_Dc, local_x2, local_x3,
                           local_slip_rate_k4, local_traction_k4,
                           local_normal_k4);

   const std::string vtu = prefix + "/FaultSurface/fault_surface_r0_c0.vtu";

   // All five _k4 fields must exist.
   TEST_ASSERT(CellDataContainsField(vtu, "normal_stress_k4"),
               "normal_stress_k4 field present");
   TEST_ASSERT(CellDataContainsField(vtu, "traction_dip_k4"),
               "traction_dip_k4 field present");
   TEST_ASSERT(CellDataContainsField(vtu, "traction_strike_k4"),
               "traction_strike_k4 field present");
   TEST_ASSERT(CellDataContainsField(vtu, "slip_rate_dip_k4"),
               "slip_rate_dip_k4 field present");
   TEST_ASSERT(CellDataContainsField(vtu, "slip_rate_strike_k4"),
               "slip_rate_strike_k4 field present");

   // Values.
   auto cell_sn_k4  = ParseCellData(vtu, "normal_stress_k4");
   auto cell_sn     = ParseCellData(vtu, "normal_stress");
   auto cell_srd_k4 = ParseCellData(vtu, "slip_rate_dip_k4");
   auto cell_srs_k4 = ParseCellData(vtu, "slip_rate_strike_k4");
   auto cell_td_k4  = ParseCellData(vtu, "traction_dip_k4");
   auto cell_ts_k4  = ParseCellData(vtu, "traction_strike_k4");

   TEST_ASSERT(static_cast<int>(cell_sn_k4.size()) == n_int,
               "normal_stress_k4 CellData length = n_int");
   TEST_ASSERT(static_cast<int>(cell_sn.size()) == n_int,
               "normal_stress CellData length = n_int");

   for (int fi = 0; fi < n_int; fi++)
   {
      const double expected_sn_k4 = 70e6 + static_cast<double>(fi) * 1e6;
      const double expected_sn_avg = 120e6;
      const double expected_srd_k4 = static_cast<double>(fi);
      const double expected_srs_k4 = 5.0 + static_cast<double>(fi);
      const double expected_td_k4  = 2e6 * static_cast<double>(fi + 1);
      const double expected_ts_k4  = -3e6 * static_cast<double>(fi + 1);

      TEST_NEAR(cell_sn_k4[fi], expected_sn_k4, 1.0,
                std::string("normal_stress_k4 cell[") + std::to_string(fi) +
                "] == 70 + fi MPa");
      TEST_NEAR(cell_sn[fi], expected_sn_avg, 1.0,
                std::string("normal_stress cell[") + std::to_string(fi) +
                "] == 120 MPa (unchanged by k4 inputs)");
      TEST_NEAR(cell_srd_k4[fi], expected_srd_k4, 1e-10,
                std::string("slip_rate_dip_k4 cell[") + std::to_string(fi) +
                "] == fi");
      TEST_NEAR(cell_srs_k4[fi], expected_srs_k4, 1e-10,
                std::string("slip_rate_strike_k4 cell[") + std::to_string(fi) +
                "] == 5 + fi");
      TEST_NEAR(cell_td_k4[fi], expected_td_k4, 1.0,
                std::string("traction_dip_k4 cell[") + std::to_string(fi) +
                "] == 2 MPa * (fi+1)");
      TEST_NEAR(cell_ts_k4[fi], expected_ts_k4, 1.0,
                std::string("traction_strike_k4 cell[") + std::to_string(fi) +
                "] == -3 MPa * (fi+1)");

      // The discriminator: normal_stress_k4 != normal_stress (by construction).
      const double delta = cell_sn_k4[fi] - cell_sn[fi];
      TEST_ASSERT(std::abs(delta) > 1e6,
                  std::string("delta normal_stress_k4 - normal_stress "
                              "for cell[") + std::to_string(fi) +
                  "] is nonzero (diagnostic is live)");
   }
}

// ---------------------------------------------------------------------------
// Test (d) — only two of three _k4 vectors populated ⇒ _k4 emission OFF.
// Protects against a subtle caller bug where the driver forgets to size
// one of the three buffers; the writer falls back to backwards-compat.
// ---------------------------------------------------------------------------
static void TestPartialK4DisablesEmission()
{
   std::cout << "\n=== R-V92-E02 (d): partial k4 inputs ⇒ no _k4 fields ===\n";

   Mesh mesh = MakeTinyTetMesh();
   Array<int> fault_faces = CollectInteriorFaces(mesh, 4);
   Array<int> empty_shared;
   const int nbf = 3;
   const int n_int = fault_faces.Size();
   const int total_dofs = n_int * nbf;

   const std::string prefix = "test_r_v92_e02_partial";
   ::mkdir(prefix.c_str(), 0755);
   ParaViewOutput<Mesh> pv(prefix, mesh, /*order=*/1);
   pv.InitFaultOutputBP5(fault_faces, empty_shared, nbf);

   Vector local_slip      (2 * total_dofs);  local_slip      = 0.0;
   Vector local_slip_rate (2 * total_dofs);  local_slip_rate = 0.0;
   Vector local_traction  (2 * total_dofs);  local_traction  = 0.0;
   Vector local_state     (total_dofs);      local_state     = 0.0;
   Vector local_normal    (total_dofs);      local_normal    = 120e6;
   Vector local_a         (0);
   Vector local_Dc        (0);
   Vector local_x2        (0);
   Vector local_x3        (0);

   // Only normal_stress_k4 populated; other two empty.
   Vector local_normal_k4 (total_dofs);       local_normal_k4 = 100e6;
   Vector local_slip_rate_k4;   // empty
   Vector local_traction_k4;    // empty

   pv.WriteFaultSurfaceVTU(prefix, /*cycle=*/0, /*time=*/0.0,
                           /*rank=*/0, /*nranks=*/1,
                           local_slip, local_slip_rate, local_traction,
                           local_state, local_normal,
                           local_a, local_Dc, local_x2, local_x3,
                           local_slip_rate_k4, local_traction_k4,
                           local_normal_k4);

   const std::string vtu = prefix + "/FaultSurface/fault_surface_r0_c0.vtu";
   TEST_ASSERT(!CellDataContainsField(vtu, "normal_stress_k4"),
               "partial k4 inputs ⇒ NO normal_stress_k4 (fail-safe)");
   TEST_ASSERT(!CellDataContainsField(vtu, "traction_dip_k4"),
               "partial k4 inputs ⇒ NO traction_dip_k4");
}

int main()
{
   std::cout << "========================================\n";
   std::cout << "TPV102 v9.2.0 §19 — R-V92-E02 stage-4 VTU diagnostic\n";
   std::cout << "========================================\n";

   TestBackwardsCompatNoK4();
   TestK4FieldsPopulatedAndDistinct();
   TestPartialK4DisablesEmission();

   std::cout << "\n========================================\n";
   std::cout << "Total:  " << num_tests << "\n";
   std::cout << "Passed: " << num_passed << "\n";
   std::cout << "Failed: " << num_failed << "\n";
   std::cout << "========================================\n";
   return (num_failed > 0) ? 1 : 0;
}

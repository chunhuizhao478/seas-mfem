// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// Unit test: --paraview-fields / SetFaultVTUFields filter on the
// fault-surface VTU writer.
//
// Contract (io/paraview_output.hpp::SetFaultVTUFields):
//
//   (a) EMPTY filter (default, no SetFaultVTUFields call) ⇒ ALL 12
//       standard fields emitted in the per-rank VTU CellData block
//       AND in the PVTU index.  Plus the 5 `_k4` fields when
//       stage-4 buffers are supplied.  Byte-identical to the
//       pre-filter writer behaviour.
//
//   (b) Non-empty filter ⇒ only the named CellData arrays appear in
//       the per-rank VTU AND the PVTU.  Fields not in the filter
//       are silently omitted.  Points / Cells / connectivity are
//       unaffected.
//
//   (c) Filter entry that does NOT match any known field is
//       silently skipped — the resulting VTU just has fewer fields
//       than the filter size.  (This matches the CLI contract:
//       typos fail quietly, producing a smaller-than-expected VTU,
//       rather than aborting.)
//
//   (d) Filter still works correctly in combination with the
//       existing `_k4` stage-4 diagnostic: requesting
//       `normal_stress_k4` alone emits only that field (the
//       non-k4 `normal_stress` is filtered out, even though the
//       writer receives data for both).

#include "mfem.hpp"
#include "../../io/paraview_output.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <set>
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

// ---------------------------------------------------------------------------
// Helpers — copy/adapted from test_fault_surface_vtu_k4.cpp.
// ---------------------------------------------------------------------------

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

static int CellDataFieldCount(const std::string &path)
{
   std::ifstream in(path);
   if (!in.is_open()) { return -1; }
   std::stringstream ss;
   ss << in.rdbuf();
   const std::string content = ss.str();
   const auto cd_open  = content.find("<CellData>");
   const auto cd_close = content.find("</CellData>");
   if (cd_open == std::string::npos || cd_close == std::string::npos)
   {
      return -1;
   }
   const std::string cd = content.substr(cd_open, cd_close - cd_open);
   int count = 0;
   size_t pos = 0;
   while ((pos = cd.find("<DataArray", pos)) != std::string::npos)
   {
      ++count;
      pos += 10;
   }
   return count;
}

// Same count helper but operating on PVTU's <PCellData> / <PDataArray>.
static int PVTUFieldCount(const std::string &path)
{
   std::ifstream in(path);
   if (!in.is_open()) { return -1; }
   std::stringstream ss;
   ss << in.rdbuf();
   const std::string content = ss.str();
   const auto cd_open  = content.find("<PCellData>");
   const auto cd_close = content.find("</PCellData>");
   if (cd_open == std::string::npos || cd_close == std::string::npos)
   {
      return -1;
   }
   const std::string cd = content.substr(cd_open, cd_close - cd_open);
   int count = 0;
   size_t pos = 0;
   while ((pos = cd.find("<PDataArray", pos)) != std::string::npos)
   {
      ++count;
      pos += 11;
   }
   return count;
}

static bool PVTUContainsField(const std::string &path,
                              const std::string &field)
{
   std::ifstream in(path);
   if (!in.is_open()) { return false; }
   std::stringstream ss;
   ss << in.rdbuf();
   const std::string content = ss.str();
   return content.find(std::string("Name=\"") + field + "\"") != std::string::npos;
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
      if (mesh.GetInteriorFaceTransformations(f) != nullptr)
      {
         out.Append(f);
      }
   }
   return out;
}

// Shared fixture — build a ParaViewOutput on a tiny tet mesh and
// prepare the local vectors WriteFaultSurfaceVTU expects.
struct Fixture
{
   Mesh mesh;
   Array<int> fault_faces;
   Array<int> empty_shared;
   int nbf;
   int n_int;
   int total_dofs;
   std::unique_ptr<ParaViewOutput<Mesh>> pv;
   std::string prefix;

   Vector local_slip, local_slip_rate, local_traction, local_state;
   Vector local_normal, local_a, local_Dc, local_x2, local_x3;
   Vector local_slip_rate_k4, local_traction_k4, local_normal_k4;

   // R-006: prepend /tmp/ so test artefacts do not pollute the working
   // tree.  Caller passes a bare scenario name; the fixture absolutises it.
   Fixture(const std::string &p)
      : mesh(MakeTinyTetMesh()),
        nbf(3),
        prefix("/tmp/" + p)
   {
      fault_faces = CollectInteriorFaces(mesh, 4);
      n_int = fault_faces.Size();
      total_dofs = n_int * nbf;

      ::mkdir(prefix.c_str(), 0755);
      pv = std::make_unique<ParaViewOutput<Mesh>>(prefix, mesh, 1);
      pv->InitFaultOutputBP5(fault_faces, empty_shared, nbf);
      // The fixture below parses both VTU and PVTU with regex.  Phase 1's
      // single binary VTU writer drops the per-rank PVTU entirely, so opt
      // into the legacy ASCII back end to keep the existing assertions
      // intact (the filter contract itself is identical between writers).
      pv->SetLegacyAsciiVTU(true);

      local_slip     .SetSize(2 * total_dofs);  local_slip      = 0.1;
      local_slip_rate.SetSize(2 * total_dofs);  local_slip_rate = 0.2;
      local_traction .SetSize(2 * total_dofs);  local_traction  = 0.3;
      local_state    .SetSize(total_dofs);      local_state     = 0.4;
      local_normal   .SetSize(total_dofs);      local_normal    = 120e6;
      local_a        .SetSize(total_dofs);      local_a         = 0.012;
      local_Dc       .SetSize(total_dofs);      local_Dc        = 0.14;
      local_x2       .SetSize(total_dofs);      local_x2        = 5.0;
      local_x3       .SetSize(total_dofs);      local_x3        = -7.5;

      local_slip_rate_k4.SetSize(2 * total_dofs); local_slip_rate_k4 = 0.25;
      local_traction_k4 .SetSize(2 * total_dofs); local_traction_k4  = 0.35;
      local_normal_k4   .SetSize(total_dofs);     local_normal_k4    = 118e6;
   }

   void Write(int cycle = 0, bool with_k4 = false)
   {
      if (with_k4)
      {
         pv->WriteFaultSurfaceVTU(prefix, cycle, 0.0,
                                  /*rank=*/0, /*nranks=*/1,
                                  local_slip, local_slip_rate, local_traction,
                                  local_state, local_normal,
                                  local_a, local_Dc, local_x2, local_x3,
                                  local_slip_rate_k4, local_traction_k4,
                                  local_normal_k4);
      }
      else
      {
         pv->WriteFaultSurfaceVTU(prefix, cycle, 0.0,
                                  /*rank=*/0, /*nranks=*/1,
                                  local_slip, local_slip_rate, local_traction,
                                  local_state, local_normal,
                                  local_a, local_Dc, local_x2, local_x3);
      }
   }

   std::string VTUPath(int cycle = 0) const
   {
      return prefix + "/FaultSurface/fault_surface_r0_c"
                    + std::to_string(cycle) + ".vtu";
   }
   std::string PVTUPath(int cycle = 0) const
   {
      return prefix + "/FaultSurface/fault_surface_c"
                    + std::to_string(cycle) + ".pvtu";
   }
};

// ---------------------------------------------------------------------------
// Test (a) — no filter set ⇒ all 12 fields emitted; PVTU lists all 12.
// Backward-compatibility gate.
// ---------------------------------------------------------------------------
static void TestEmptyFilterEmitsAllFields()
{
   std::cout << "\n=== (a) empty filter emits all 12 fields ===\n";
   Fixture f("test_pv_field_filter_empty");
   // Explicit: do NOT call SetFaultVTUFields — default state.
   f.Write();

   TEST_ASSERT(CellDataFieldCount(f.VTUPath()) == 12,
               "VTU CellData has 12 fields (no filter)");
   TEST_ASSERT(PVTUFieldCount(f.PVTUPath()) == 12,
               "PVTU PCellData has 12 fields (no filter)");
   for (const char *name : {"slip_dip", "slip_strike", "slip_rate_dip",
                            "slip_rate_strike", "traction_dip",
                            "traction_strike", "state_variable",
                            "normal_stress", "param_a", "param_Dc",
                            "fault_x2", "fault_x3"})
   {
      TEST_ASSERT(CellDataContainsField(f.VTUPath(), name),
                  std::string("VTU contains field: ") + name);
   }
}

// ---------------------------------------------------------------------------
// Test (b) — filter = {slip_rate_strike} ⇒ exactly one field in VTU
// and PVTU; other named fields absent.  This is the BP5 500-yr
// sbatch usage pattern.
// ---------------------------------------------------------------------------
static void TestSingleFieldFilter()
{
   std::cout << "\n=== (b) single-field filter (slip_rate_strike) ===\n";
   Fixture f("test_pv_field_filter_single");
   f.pv->SetFaultVTUFields({"slip_rate_strike"});
   f.Write();

   TEST_ASSERT(CellDataFieldCount(f.VTUPath()) == 1,
               "VTU CellData has exactly 1 field");
   TEST_ASSERT(CellDataContainsField(f.VTUPath(), "slip_rate_strike"),
               "VTU contains slip_rate_strike");
   for (const char *name : {"slip_dip", "slip_strike", "slip_rate_dip",
                            "traction_dip", "traction_strike",
                            "state_variable", "normal_stress",
                            "param_a", "param_Dc", "fault_x2", "fault_x3"})
   {
      TEST_ASSERT(!CellDataContainsField(f.VTUPath(), name),
                  std::string("VTU does NOT contain ") + name);
   }
   TEST_ASSERT(PVTUFieldCount(f.PVTUPath()) == 1,
               "PVTU PCellData has exactly 1 field");
   TEST_ASSERT(PVTUContainsField(f.PVTUPath(), "slip_rate_strike"),
               "PVTU lists slip_rate_strike");
   TEST_ASSERT(!PVTUContainsField(f.PVTUPath(), "traction_strike"),
               "PVTU does NOT list traction_strike");
}

// ---------------------------------------------------------------------------
// Test (c) — filter with a few requested fields + one unknown name.
// Unknown name is silently skipped; known names emitted.
// ---------------------------------------------------------------------------
static void TestMultiFieldAndUnknown()
{
   std::cout << "\n=== (c) multi-field filter + unknown name skip ===\n";
   Fixture f("test_pv_field_filter_multi");
   f.pv->SetFaultVTUFields({"slip_rate_strike", "normal_stress",
                            "DOES_NOT_EXIST"});
   f.Write();

   TEST_ASSERT(CellDataFieldCount(f.VTUPath()) == 2,
               "VTU CellData has exactly 2 fields (unknown name skipped)");
   TEST_ASSERT(CellDataContainsField(f.VTUPath(), "slip_rate_strike"),
               "VTU contains slip_rate_strike");
   TEST_ASSERT(CellDataContainsField(f.VTUPath(), "normal_stress"),
               "VTU contains normal_stress");
   TEST_ASSERT(!CellDataContainsField(f.VTUPath(), "DOES_NOT_EXIST"),
               "VTU does NOT contain DOES_NOT_EXIST");
   TEST_ASSERT(PVTUFieldCount(f.PVTUPath()) == 2,
               "PVTU has exactly 2 fields");
}

// ---------------------------------------------------------------------------
// Test (d) — filter + _k4 buffers.  Requesting only `normal_stress_k4`
// excludes the non-k4 `normal_stress` even though both data arrays
// are supplied.
// ---------------------------------------------------------------------------
static void TestFilterWithK4Buffers()
{
   std::cout << "\n=== (d) filter + stage-4 buffers ===\n";
   Fixture f("test_pv_field_filter_k4");
   f.pv->SetFaultVTUFields({"normal_stress_k4"});
   f.Write(/*cycle=*/0, /*with_k4=*/true);

   TEST_ASSERT(CellDataFieldCount(f.VTUPath()) == 1,
               "VTU has exactly 1 field (normal_stress_k4 only)");
   TEST_ASSERT(CellDataContainsField(f.VTUPath(), "normal_stress_k4"),
               "VTU contains normal_stress_k4");
   TEST_ASSERT(!CellDataContainsField(f.VTUPath(), "normal_stress"),
               "VTU does NOT contain non-k4 normal_stress");
   TEST_ASSERT(!CellDataContainsField(f.VTUPath(), "slip_rate_strike_k4"),
               "VTU does NOT contain unlisted slip_rate_strike_k4");
   TEST_ASSERT(PVTUFieldCount(f.PVTUPath()) == 1,
               "PVTU has exactly 1 field");
}

// ---------------------------------------------------------------------------
// Test (e) — clearing the filter (replacing with empty set) restores
// full emission behaviour on a subsequent write.
// ---------------------------------------------------------------------------
static void TestClearingFilterRestoresAll()
{
   std::cout << "\n=== (e) clearing filter restores all-field emission ===\n";
   Fixture f("test_pv_field_filter_clear");
   f.pv->SetFaultVTUFields({"slip_rate_strike"});
   f.Write(/*cycle=*/0);

   TEST_ASSERT(CellDataFieldCount(f.VTUPath(0)) == 1,
               "cycle 0 filtered (1 field)");

   // Clear to empty set ⇒ default emit-all behaviour.
   f.pv->SetFaultVTUFields({});
   f.Write(/*cycle=*/1);

   TEST_ASSERT(CellDataFieldCount(f.VTUPath(1)) == 12,
               "cycle 1 un-filtered (12 fields)");
}

int main()
{
   std::cout << "========================================\n";
   std::cout << "Fault-surface VTU --paraview-fields filter tests\n";
   std::cout << "========================================\n";

   TestEmptyFilterEmitsAllFields();
   TestSingleFieldFilter();
   TestMultiFieldAndUnknown();
   TestFilterWithK4Buffers();
   TestClearingFilterRestoresAll();

   std::cout << "\n========================================\n";
   std::cout << "Total:  " << num_tests  << "\n";
   std::cout << "Passed: " << num_passed << "\n";
   std::cout << "Failed: " << num_failed << "\n";
   std::cout << "========================================\n";
   return (num_failed > 0) ? 1 : 0;
}

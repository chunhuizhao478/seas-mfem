// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 pepper-bug unit-test plan (2026-04-22) — Phase 4 P4-T2:
// fault DOF-offset invariants (serial).
//
// ============================================================================
// Motivation
// ============================================================================
// WaveOperator maintains two maps:
//   fault_face_dof_offset_[f]   — offset into fault_dof_data_ for
//                                 interior fault face f.
//   shared_fault_dof_offset_[sf] — offset for shared fault face sf.
// The union of { [offset, offset + nqp_per_face) } entries must be a
// DISJOINT contiguous range covering [0, fault_dof_data_.size()).
//
// A subtle bug in SetFaultDOFData (wave_operator.hpp:278-292) could
// produce overlapping or non-contiguous ranges, silently corrupting
// the DOFData lookup at every fault QP.  This test gates that
// invariant on a serial 2-tet fixture (trivially 1 interior, 0 shared)
// and on a 4-tet fixture (multiple interior fault faces).
//
// ============================================================================
// Usage:  ./seas_test_fault_dof_offset_invariants   (serial, <1s)

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <set>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST(cond, msg) do { \
   num_tests++; \
   if (cond) { \
      num_passed++; \
      std::cout << "  PASSED: " << msg << "\n"; \
   } else { \
      num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; \
   } \
} while (0)

// 2-tet mesh with one interior fault face.
static Mesh BuildTwoTetFaultMesh()
{
   Mesh mesh(3, 5, 2, 0);
   real_t verts[5][3] = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0},
      {0.0,  1.0, 0.0}, {0.0, -1.0, 0.0},
   };
   for (int v = 0; v < 5; v++) { mesh.AddVertex(verts[v]); }
   mesh.AddTet(0, 1, 2, 4, 1);
   mesh.AddTet(0, 1, 2, 3, 1);
   mesh.FinalizeTopology();
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         real_t cy = 0;
         for (int v = 0; v < fv.Size(); v++)
         { cy += mesh.GetVertex(fv[v])[1]; }
         cy /= fv.Size();
         if (std::abs(cy) < 1e-10)
         { mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3); }
         continue;
      }
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 1);
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// 8-element hex mesh split via MakeCartesian3D; y=0 layer becomes fault.
// Produces 4 interior fault faces (each hex face on y=0).
static Mesh BuildEightHexMultiFaultMesh()
{
   const real_t L = 1.0;
   Mesh mesh = Mesh::MakeCartesian3D(2, 2, 2, Element::HEXAHEDRON,
                                     L, 2.0*L, L);
   for (int v = 0; v < mesh.GetNV(); v++)
   {
      real_t *x = mesh.GetVertex(v);
      x[1] -= L;   // shift so y ∈ [-L, L]
   }
   const real_t tol = 1e-8;
   for (int be = 0; be < mesh.GetNBE(); be++)
   { mesh.SetBdrAttribute(be, 1); }
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(f);
      if (!ftr) { continue; }
      const IntegrationPoint &ip = Geometries.GetCenter(ftr->GetGeometryType());
      ftr->Face->SetIntPoint(&ip);
      Vector c(3); ftr->Face->Transform(ip, c);
      if (std::abs(c(1)) > tol) { continue; }
      Array<int> verts; mesh.GetFaceVertices(f, verts);
      if (verts.Size() == 4)
      { mesh.AddBdrQuad(verts[0], verts[1], verts[2], verts[3], 3); }
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

static void RunFixture(const std::string &label, Mesh &mesh)
{
   std::cout << "\n-- Fixture " << label << " --\n";

   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {};

   const int order = 1;
   WaveOperator<Mesh> wave(mesh, order,
                           TPV102Params::lambda,
                           TPV102Params::mu,
                           TPV102Params::rho, bc);

   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   const Array<int> &shr_faces = wave.GetFaultSharedFaces();
   const int n_int = int_faces.Size();
   const int n_shr = shr_faces.Size();
   TEST(n_shr == 0, label + ": serial build has 0 shared fault faces");
   std::cout << "    interior fault faces = " << n_int << "\n";

   if (n_int == 0) { return; }

   // Determine nqp_per_face from the first interior fault face.
   auto *ftr0 = mesh.GetInteriorFaceTransformations(int_faces[0]);
   const int nqp =
      IntRules.Get(ftr0->GetGeometryType(), 2*order).GetNPoints();

   // Build fault_coords + DOFData so WaveOperator populates the offset maps.
   std::vector<Vector> fault_coords;
   for (int i = 0; i < n_int; i++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir =
         IntRules.Get(ftr->GetGeometryType(), 2*order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   }
   std::vector<DOFData> dof_data;
   const int num_fault = n_int * nqp;
   InitializeFaultDOFs(dof_data, num_fault, fault_coords);

   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp);

   // Invariant 1: fault_face_dof_offset_ has exactly n_int entries.
   const auto &ifmap = wave.GetFaultFaceDofOffset();
   const auto &sfmap = wave.GetSharedFaultDofOffset();
   TEST((int)ifmap.size() == n_int,
        label + ": fault_face_dof_offset_ size == n_int");
   TEST((int)sfmap.size() == 0,
        label + ": shared_fault_dof_offset_ empty (serial)");

   // Invariant 2: every mesh face in int_faces has an offset entry.
   for (int i = 0; i < n_int; i++)
   {
      auto it = ifmap.find(int_faces[i]);
      TEST(it != ifmap.end(),
           label + ": face " + std::to_string(int_faces[i]) +
           " has offset entry");
   }

   // Invariant 3: union of [offset, offset+nqp) ranges is disjoint AND
   // covers exactly [0, n_int * nqp).
   std::set<int> seen;
   bool all_in_range = true;
   bool no_dups = true;
   for (const auto &p : ifmap)
   {
      const int off = p.second;
      for (int q = 0; q < nqp; q++)
      {
         const int idx = off + q;
         if (idx < 0 || idx >= num_fault) { all_in_range = false; }
         if (!seen.insert(idx).second)    { no_dups = false; }
      }
   }
   TEST(all_in_range, label + ": every DOF offset + q in [0, num_fault)");
   TEST(no_dups,      label + ": no DOF offset duplicates (disjoint)");
   TEST((int)seen.size() == num_fault,
        label + ": union covers exactly [0, num_fault)");

   // Invariant 4: SetFaultDOFData is idempotent (no silent corruption
   // on a second call).
   wave.SetFaultDOFData(&dof_data, nqp);
   const auto &ifmap2 = wave.GetFaultFaceDofOffset();
   bool matches = (ifmap.size() == ifmap2.size());
   if (matches)
   {
      for (const auto &p : ifmap)
      {
         auto it = ifmap2.find(p.first);
         if (it == ifmap2.end() || it->second != p.second)
         { matches = false; break; }
      }
   }
   TEST(matches,
        label + ": SetFaultDOFData is idempotent (second call preserves map)");
}

int main()
{
   std::cout << "\n=== TPV102 Pepper-Bug Phase 4 P4-T2: "
             << "Fault DOF-Offset Invariants ===\n";

   {
      Mesh m1 = BuildTwoTetFaultMesh();
      RunFixture("2-tet (1 interior fault face)", m1);
   }
   {
      Mesh m2 = BuildEightHexMultiFaultMesh();
      RunFixture("8-hex (4 interior fault faces)", m2);
   }

   std::cout << "\n========================================\n";
   std::cout << "  Phase 4 P4-T2 results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}

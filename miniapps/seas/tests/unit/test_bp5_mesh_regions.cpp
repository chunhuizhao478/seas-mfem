// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

// Unit tests for BP5 mesh regions and parameter mapping:
//   - Boundary attributes (1-6)
//   - Fault region parameters (a, Dc, tau_pre)
//   - Fault face → element mapping validity

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../fault/fault_geometry.hpp"
#include "../../fault/rate_state_fault.hpp"
#include "../../friction/dieterich_ruina.hpp"
#include "../../friction/state_evolution.hpp"
#include "../../config/bp5_params.hpp"
#include "test_macros.hpp"

#include <iostream>
#include <cmath>
#include <memory>
#include <set>
#include <map>

using namespace mfem;
using namespace mfem::seas;

// Helper: Create 3D hex mesh
static Mesh Create3DMesh(int nx, int ny, int nz,
                          real_t Lx, real_t Ly, real_t Lz)
{
   Mesh mesh = Mesh::MakeCartesian3D(2 * nx, 2 * ny, nz,
                                      Element::HEXAHEDRON,
                                      2.0 * Lx, 2.0 * Ly, Lz);

   for (int i = 0; i < mesh.GetNV(); i++)
   {
      real_t *v = mesh.GetVertex(i);
      v[0] -= Lx;
      v[1] -= Ly;
   }

   for (int be = 0; be < mesh.GetNBE(); be++)
   {
      ElementTransformation *T = mesh.GetBdrElementTransformation(be);
      const IntegrationPoint &ip = Geometries.GetCenter(T->GetGeometryType());
      T->SetIntPoint(&ip);
      Vector center(3);
      T->Transform(ip, center);

      real_t tol = 1e-6;
      if (std::abs(center(0) - (-Lx)) < tol)      { mesh.SetBdrAttribute(be, 1); }
      else if (std::abs(center(0) - Lx) < tol)     { mesh.SetBdrAttribute(be, 2); }
      else if (std::abs(center(1) - Ly) < tol)     { mesh.SetBdrAttribute(be, 3); }
      else if (std::abs(center(1) - (-Ly)) < tol)  { mesh.SetBdrAttribute(be, 4); }
      else if (std::abs(center(2) - 0.0) < tol)    { mesh.SetBdrAttribute(be, 5); }
      else if (std::abs(center(2) - Lz) < tol)     { mesh.SetBdrAttribute(be, 6); }
   }

   mesh.SetAttributes();
   return mesh;
}

// =============================================================================
// Test 1: Boundary Attributes
// =============================================================================
void TestBoundaryAttributes()
{
   std::cout << "\n--- Test: Boundary Attributes ---\n";

   real_t Lx = 50e3, Ly = 60e3, Lz = 40e3;
   Mesh mesh = Create3DMesh(2, 2, 1, Lx, Ly, Lz);

   // Count faces per boundary attribute
   std::map<int, int> attr_count;
   for (int i = 0; i < mesh.GetNBE(); i++)
   {
      attr_count[mesh.GetBdrAttribute(i)]++;
   }

   // Should have all 6 boundary attributes
   TEST_ASSERT(attr_count.count(1) > 0, "Boundary attr 1 (x=-Lx) exists");
   TEST_ASSERT(attr_count.count(2) > 0, "Boundary attr 2 (x=+Lx) exists");
   TEST_ASSERT(attr_count.count(3) > 0, "Boundary attr 3 (y=+Ly) exists");
   TEST_ASSERT(attr_count.count(4) > 0, "Boundary attr 4 (y=-Ly) exists");
   TEST_ASSERT(attr_count.count(5) > 0, "Boundary attr 5 (z=0) exists");
   TEST_ASSERT(attr_count.count(6) > 0, "Boundary attr 6 (z=Lz) exists");

   int total = 0;
   for (auto &p : attr_count)
   {
      std::cout << "  Attr " << p.first << ": " << p.second << " faces\n";
      total += p.second;
   }
   TEST_ASSERT(total == mesh.GetNBE(), "All boundary elements have attributes");
}

// =============================================================================
// Test 2: Fault Region Parameters
// =============================================================================
void TestFaultRegionParameters()
{
   std::cout << "\n--- Test: Fault Region Parameters ---\n";

   BP5Params params;
   real_t Lx = 50e3, Ly = 60e3, Lz = 40e3;
   Mesh mesh = Create3DMesh(2, 4, 2, Lx, Ly, Lz);

   ElasticityDomainOperator<Mesh> domain_op(
      mesh, 1, params.lambda(), params.mu(),
      params.Vp, params.Wf, params.lf, DGMethod::BR2);

   int nf = domain_op.GetNumFaultDOFs();
   if (nf == 0)
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   FaultGeometry<Mesh> fault_geom(domain_op, params);

   const Vector &a_vals = fault_geom.GetAValues();
   const Vector &dc_vals = fault_geom.GetDcValues();
   const Vector &tau_pre = fault_geom.GetTauPre();

   TEST_ASSERT(a_vals.Size() == nf, "a values size = nf");
   TEST_ASSERT(dc_vals.Size() == nf, "Dc values size = nf");
   TEST_ASSERT(tau_pre.Size() == 2 * nf, "tau_pre size = 2*nf");

   // Check a values are in expected range
   TEST_ASSERT(a_vals.Min() >= params.a0 - 1e-10,
               "a_min >= a0");
   TEST_ASSERT(a_vals.Max() <= params.amax + 1e-10,
               "a_max <= amax");

   // Check Dc values are positive
   TEST_ASSERT(dc_vals.Min() > 0.0, "All Dc > 0");

   // Check tau_pre is finite and positive in magnitude
   bool all_finite = true;
   bool any_positive = false;
   for (int i = 0; i < nf; i++)
   {
      real_t td = tau_pre(2 * i);
      real_t ts = tau_pre(2 * i + 1);
      if (!std::isfinite(td) || !std::isfinite(ts)) { all_finite = false; }
      real_t mag = std::sqrt(td * td + ts * ts);
      if (mag > 0.0) { any_positive = true; }
   }
   TEST_ASSERT(all_finite, "All tau_pre components finite");
   TEST_ASSERT(any_positive, "At least one tau_pre > 0");

   std::cout << "  nf=" << nf
             << ", a_range=[" << a_vals.Min() << ", " << a_vals.Max() << "]"
             << ", Dc_range=[" << dc_vals.Min() << ", " << dc_vals.Max() << "]"
             << "\n";
}

// =============================================================================
// Test 3: Fault Face → Element Mapping
// =============================================================================
void TestFaultFaceElementMapping()
{
   std::cout << "\n--- Test: Fault Face → Element Mapping ---\n";

   BP5Params params;
   real_t Lx = 50e3, Ly = 60e3, Lz = 40e3;
   Mesh mesh = Create3DMesh(2, 2, 1, Lx, Ly, Lz);

   ElasticityDomainOperator<Mesh> domain_op(
      mesh, 1, params.lambda(), params.mu(),
      params.Vp, params.Wf, params.lf, DGMethod::BR2);

   int nf = domain_op.GetNumFaultDOFs();
   if (nf == 0)
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   const Array<int> &fault_int_faces = domain_op.GetFaultInteriorFaces();
   int nf_int = fault_int_faces.Size();

   TEST_ASSERT(nf_int > 0, "Has interior fault faces");

   int NE = mesh.GetNE();
   bool all_valid = true;
   bool all_distinct = true;
   bool all_on_fault = true;

   for (int i = 0; i < nf_int; i++)
   {
      FaceElementTransformations *FTr =
         mesh.GetInteriorFaceTransformations(fault_int_faces[i]);
      int e1 = FTr->Elem1No;
      int e2 = FTr->Elem2No;

      if (e1 < 0 || e1 >= NE || e2 < 0 || e2 >= NE) { all_valid = false; }
      if (e1 == e2) { all_distinct = false; }

      // Check face centroid at x ≈ 0
      const IntegrationPoint &ip =
         Geometries.GetCenter(FTr->GetGeometryType());
      FTr->SetIntPoint(&ip);
      Vector fc(3);
      FTr->Transform(ip, fc);
      if (std::abs(fc(0)) > 1.0) { all_on_fault = false; }
   }

   TEST_ASSERT(all_valid, "All element indices valid");
   TEST_ASSERT(all_distinct, "elem1 != elem2");
   TEST_ASSERT(all_on_fault, "All face centroids at x ≈ 0");

   std::cout << "  nf_int=" << nf_int << ", NE=" << NE << "\n";
}

// =============================================================================
// Main
// =============================================================================
int main()
{
   std::cout << "========================================\n";
   std::cout << "  BP5 Mesh Regions Unit Tests\n";
   std::cout << "========================================\n";

   TestBoundaryAttributes();
   TestFaultRegionParameters();
   TestFaultFaceElementMapping();

   TEST_PRINT_RESULTS();
   return num_failed;
}

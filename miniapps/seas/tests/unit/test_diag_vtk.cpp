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

// Unit tests for diagnostic VTK output:
//   - Fault face → element mapping
//   - State-based field values after initialization
//   - L2 p=0 field correctness

#include "mfem.hpp"
#include "../../solver/seas_operator.hpp"
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
#include <vector>
#include <set>

using namespace mfem;
using namespace mfem::seas;

// Helper: Create 3D hex mesh (same as test_bp5_integration.cpp)
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

// Shared fixture
struct DiagVTKFixture
{
   BP5Params params;
   std::unique_ptr<Mesh> mesh;
   std::unique_ptr<ElasticityDomainOperator<Mesh>> domain_op;
   std::unique_ptr<FaultGeometry<Mesh>> fault_geom;
   std::unique_ptr<DieterichRuinaFriction> friction;
   std::unique_ptr<AgingLawPsi> evolution;
   std::unique_ptr<RateStateFaultOperator<Mesh, 2>> fault_op;
   std::unique_ptr<BP5SEASOp> seas_op;
   int nf = 0;

   bool Setup()
   {
      real_t Lx = 50e3, Ly = 60e3, Lz = 40e3;
      mesh = std::make_unique<Mesh>(Create3DMesh(1, 1, 1, Lx, Ly, Lz));

      domain_op = std::make_unique<ElasticityDomainOperator<Mesh>>(
         *mesh, 1, params.lambda(), params.mu(),
         params.Vp, params.Wf, params.lf, DGMethod::BR2);

      nf = domain_op->GetNumFaultDOFs();
      if (nf == 0) { return false; }

      fault_geom = std::make_unique<FaultGeometry<Mesh>>(*domain_op, params);

      DieterichRuinaFriction::Constants fc;
      fc.V0 = params.V0;
      fc.f0 = params.f0;
      fc.b = params.b;
      fc.Dc = params.L0;
      friction = std::make_unique<DieterichRuinaFriction>(fc);

      evolution = std::make_unique<AgingLawPsi>(params.b, params.V0, params.f0);

      fault_op = std::make_unique<RateStateFaultOperator<Mesh, 2>>(
         fault_geom.get(), friction.get(), evolution.get(), params);

      seas_op = std::make_unique<BP5SEASOp>(domain_op.get(), fault_op.get());

      return true;
   }
};

// =============================================================================
// Test 1: Fault Face → Element Mapping
// =============================================================================
void TestFaultFaceElementMapping()
{
   std::cout << "\n--- Test: Fault Face → Element Mapping ---\n";

   DiagVTKFixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   const Array<int> &fault_int_faces = fix.domain_op->GetFaultInteriorFaces();
   int nf_int = fault_int_faces.Size();

   TEST_ASSERT(nf_int > 0, "Has interior fault faces");

   // Build mapping
   std::vector<int> elem1(nf_int), elem2(nf_int);
   int NE = fix.mesh->GetNE();
   bool all_valid = true;
   bool all_distinct = true;

   for (int i = 0; i < nf_int; i++)
   {
      FaceElementTransformations *FTr =
         fix.mesh->GetInteriorFaceTransformations(fault_int_faces[i]);
      elem1[i] = FTr->Elem1No;
      elem2[i] = FTr->Elem2No;

      if (elem1[i] < 0 || elem1[i] >= NE) { all_valid = false; }
      if (elem2[i] < 0 || elem2[i] >= NE) { all_valid = false; }
      if (elem1[i] == elem2[i]) { all_distinct = false; }
   }

   TEST_ASSERT(all_valid, "All element indices valid");
   TEST_ASSERT(all_distinct, "elem1 != elem2 for all faces");

   // Verify face centroids at x ≈ 0
   bool all_on_fault = true;
   for (int i = 0; i < nf_int; i++)
   {
      // Get face centroid
      FaceElementTransformations *FTr =
         fix.mesh->GetInteriorFaceTransformations(fault_int_faces[i]);
      const IntegrationPoint &ip =
         Geometries.GetCenter(FTr->GetGeometryType());
      FTr->SetIntPoint(&ip);
      Vector face_center(3);
      FTr->Transform(ip, face_center);
      if (std::abs(face_center(0)) > 1.0) { all_on_fault = false; }
   }

   TEST_ASSERT(all_on_fault, "All fault face centroids at x ≈ 0");

   // Verify L2 p=0 field has correct size
   L2_FECollection l2_fec(0, 3);
   FiniteElementSpace l2_fes(fix.mesh.get(), &l2_fec);
   TEST_ASSERT(l2_fes.GetNDofs() == NE,
               "L2 p=0 has one DOF per element");

   // Map test value to fault-adjacent elements
   GridFunction test_field(&l2_fes);
   test_field = 0.0;
   for (int i = 0; i < nf_int; i++)
   {
      test_field(elem1[i]) = 1.0;
      test_field(elem2[i]) = 1.0;
   }

   int n_nonzero = 0;
   for (int i = 0; i < NE; i++)
   {
      if (test_field(i) > 0.5) { n_nonzero++; }
   }

   // Each fault face marks 2 elements (may overlap)
   TEST_ASSERT(n_nonzero > 0, "At least some elements marked as fault-adjacent");
   TEST_ASSERT(n_nonzero <= 2 * nf_int,
               "At most 2*nf_int elements marked");

   std::cout << "  nf_int=" << nf_int << ", fault-adjacent elements="
             << n_nonzero << "\n";
}

// =============================================================================
// Test 2: State-Based Field Values After Init
// =============================================================================
void TestDiagVTKFieldValues()
{
   std::cout << "\n--- Test: Diagnostic VTK Field Values ---\n";

   DiagVTKFixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   // Initialize
   Vector state(fix.fault_op->StateSize());
   fix.seas_op->SetInitialCondition(state);

   int N = fix.nf;

   // Extract slip
   Vector slip;
   fix.fault_op->GetSlip(state, slip);
   TEST_ASSERT(slip.Size() == 2 * N, "Slip vector size = 2*N");

   // Initial slip should be zero
   real_t slip_norm = slip.Norml2();
   TEST_NEAR(slip_norm, 0.0, 1e-15, "Initial slip is zero");

   // Extract theta (which is psi in psi-space)
   Vector theta;
   fix.fault_op->GetTheta(state, theta);
   TEST_ASSERT(theta.Size() == N, "Theta vector size = N");
   TEST_ASSERT(theta.Min() > 0.0, "All theta/psi > 0");

   // Slip rate
   const Vector &V_rate = fix.fault_op->GetSlipRate();
   TEST_ASSERT(V_rate.Size() == 2 * N, "Slip rate vector size = 2*N");

   // V_magnitude should be > 0
   bool all_V_positive = true;
   for (int i = 0; i < N; i++)
   {
      real_t vm = std::sqrt(V_rate(2*i)*V_rate(2*i) +
                            V_rate(2*i+1)*V_rate(2*i+1));
      if (vm <= 0.0) { all_V_positive = false; }
   }
   TEST_ASSERT(all_V_positive, "All V_magnitude > 0");

   // Traction
   const Vector &traction = fix.seas_op->GetTraction();
   TEST_ASSERT(traction.Size() == 2 * N, "Traction vector size = 2*N");

   // Total stress = tau_pre + traction should be finite and reasonable
   const Vector &tau_pre = fix.fault_geom->GetTauPre();
   bool all_tau_finite = true;
   for (int i = 0; i < 2 * N; i++)
   {
      real_t total = tau_pre(i) + traction(i);
      if (!std::isfinite(total)) { all_tau_finite = false; }
   }
   TEST_ASSERT(all_tau_finite, "All total stress components finite");

   // Check total stress magnitude is physically reasonable (5-25 MPa)
   bool any_reasonable = false;
   for (int i = 0; i < N; i++)
   {
      real_t td = tau_pre(2*i) + traction(2*i);
      real_t ts = tau_pre(2*i+1) + traction(2*i+1);
      real_t tm = std::sqrt(td*td + ts*ts);
      if (tm > 5e6 && tm < 30e6) { any_reasonable = true; }
   }
   TEST_ASSERT(any_reasonable,
               "At least one stress magnitude in 5-30 MPa range");

   // a values from FaultGeometry
   const Vector &a_vals = fix.fault_geom->GetAValues();
   TEST_ASSERT(a_vals.Size() == N, "a values vector size = N");
   TEST_ASSERT(a_vals.Min() > 0.0, "All a values > 0");

   // Dc values
   const Vector &dc_vals = fix.fault_geom->GetDcValues();
   TEST_ASSERT(dc_vals.Size() == N, "Dc values vector size = N");
   TEST_ASSERT(dc_vals.Min() > 0.0, "All Dc values > 0");

   std::cout << "  N=" << N << ", slip_norm=" << slip_norm
             << ", theta_min=" << theta.Min()
             << ", a_min=" << a_vals.Min()
             << ", dc_min=" << dc_vals.Min() << "\n";
}

// =============================================================================
// Test 3: L2 p=0 Field Mapping Consistency
// =============================================================================
void TestL2FieldMapping()
{
   std::cout << "\n--- Test: L2 p=0 Field Mapping Consistency ---\n";

   DiagVTKFixture fix;
   if (!fix.Setup())
   {
      std::cout << "  (Skipped: no fault faces found)\n";
      return;
   }

   // Initialize
   Vector state(fix.fault_op->StateSize());
   fix.seas_op->SetInitialCondition(state);

   const Array<int> &fault_int_faces = fix.domain_op->GetFaultInteriorFaces();
   int nf_int = fault_int_faces.Size();

   // Build mapping
   std::vector<int> elem1(nf_int), elem2(nf_int);
   for (int i = 0; i < nf_int; i++)
   {
      FaceElementTransformations *FTr =
         fix.mesh->GetInteriorFaceTransformations(fault_int_faces[i]);
      elem1[i] = FTr->Elem1No;
      elem2[i] = FTr->Elem2No;
   }

   // Map a_vals to L2 p=0 field and verify consistency
   const Vector &a_vals = fix.fault_geom->GetAValues();

   L2_FECollection l2_fec(0, 3);
   FiniteElementSpace l2_fes(fix.mesh.get(), &l2_fec);
   GridFunction a_field(&l2_fes);
   a_field = 0.0;

   for (int i = 0; i < nf_int; i++)
   {
      a_field(elem1[i]) = a_vals(i);
      a_field(elem2[i]) = a_vals(i);
   }

   // Check that mapped elements match FaultGeometry values
   bool all_match = true;
   for (int i = 0; i < nf_int; i++)
   {
      if (std::abs(a_field(elem1[i]) - a_vals(i)) > 1e-15) { all_match = false; }
      if (std::abs(a_field(elem2[i]) - a_vals(i)) > 1e-15) { all_match = false; }
   }
   TEST_ASSERT(all_match, "L2 field matches FaultGeometry a values");

   // Non-fault elements should be zero
   std::set<int> fault_elems;
   for (int i = 0; i < nf_int; i++)
   {
      fault_elems.insert(elem1[i]);
      fault_elems.insert(elem2[i]);
   }

   bool non_fault_zero = true;
   for (int i = 0; i < fix.mesh->GetNE(); i++)
   {
      if (fault_elems.count(i) == 0 && std::abs(a_field(i)) > 1e-15)
      {
         non_fault_zero = false;
      }
   }
   TEST_ASSERT(non_fault_zero, "Non-fault elements have a = 0");

   std::cout << "  nf_int=" << nf_int << ", fault_elems=" << fault_elems.size()
             << "\n";
}

// =============================================================================
// Main
// =============================================================================
int main()
{
   std::cout << "========================================\n";
   std::cout << "  Diagnostic VTK Unit Tests\n";
   std::cout << "========================================\n";

   TestFaultFaceElementMapping();
   TestDiagVTKFieldValues();
   TestL2FieldMapping();

   TEST_PRINT_RESULTS();
   return num_failed;
}

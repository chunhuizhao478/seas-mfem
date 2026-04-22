// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 pepper-bug unit-test plan (2026-04-22) — Test 5 (NEW):
// Shared-fault Pelties-7 ANALYTIC probe.
//
// ============================================================================
// Motivation
// ============================================================================
// `test_shared_fault_dof_data_consistency` Phase D asserts that both
// ranks produce bit-identical DOFData on shared fault QPs — but does
// NOT check whether the value matches the analytic Pelties 7 prediction.
// If both ranks compute the SAME WRONG value (e.g., a systemic sign
// bug or wrong canonical-frame orientation), Phase D passes while the
// output is still wrong.
//
// This test closes that gap.  On a 4-rank 4-km tet fixture with shared
// fault faces, it drives Q with an antisymmetric velocity perturbation
// (VX = sgn(y) · V_test) that makes the Pelties 7c trial traction
// CLOSED-FORM: tau_2_trial = Zs · V_test (same formula as
// test_interior_fault_flux_path T2 but on shared faces).
//
// Two assertions per shared fault QP:
//   I1. |tau2_corr - tau2_0 - tau_2_analytic| < 10% · |tau_2_analytic|
//   I2. tau1_corr (dip) stays ~0 within numerical noise (no leakage
//       between strike and dip channels — the pepper-plot signature).
//
// ============================================================================
// Fixture
// ============================================================================
// Reuses the 4-km Cartesian tet box from
// test_shared_fault_dof_data_consistency (n = 2·2·2 hex → ~192 tets),
// partitioned to 4 MPI ranks so ParMETIS reliably cuts across the
// fault at y = 0.  Fault attr = 3, absorbing attr = 5, natural attr = 1
// (top at z = 0).
//
// ============================================================================
// Failure interpretation
// ============================================================================
// I1 FAIL  ⇒ Shared-fault branch does NOT match Pelties 7c.  Directly
//            confirms H-V92-P variant where BOTH ranks compute the same
//            wrong value.  Contrast with Phase D (A=B) — this is A=truth.
// I2 FAIL  ⇒ Strike→dip leakage on the shared branch — canonical-frame
//            bug.  If interior P3-T3 passes but this fails, the bug is
//            in the shared-branch canonical frame construction only.
//
// ============================================================================
// Usage:  mpirun -np 4 ./seas_test_shared_fault_flux_path_analytic

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../dynamic/tpv102_setup_total.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

// 4-km Cartesian tet box, fault at y=0 (attr 3), top at z=0 (attr 1),
// absorbing elsewhere (attr 5).  Matches test_shared_fault_dof_data_consistency.
static Mesh BuildFixtureMesh()
{
   const int nx = 2, ny = 2, nz = 2;
   const double Lx = 2.0e3, Ly = 2.0e3, Lz = 4.0e3;
   Mesh mesh = Mesh::MakeCartesian3D(2*nx, 2*ny, nz,
                                     Element::TETRAHEDRON,
                                     2.0*Lx, 2.0*Ly, Lz);
   for (int v = 0; v < mesh.GetNV(); v++)
   {
      real_t *x = mesh.GetVertex(v);
      x[0] -= Lx; x[1] -= Ly; x[2] -= Lz;
   }
   const real_t tol = 1e-6;
   for (int be = 0; be < mesh.GetNBE(); be++)
   {
      auto *Tr = mesh.GetBdrElementTransformation(be);
      const IntegrationPoint &ip = Geometries.GetCenter(Tr->GetGeometryType());
      Tr->SetIntPoint(&ip);
      Vector c(3); Tr->Transform(ip, c);
      if (std::abs(c(2)) < tol) { mesh.SetBdrAttribute(be, 1); }
      else                      { mesh.SetBdrAttribute(be, 5); }
   }
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(f);
      if (!ftr) { continue; }
      const IntegrationPoint &ip = Geometries.GetCenter(ftr->GetGeometryType());
      ftr->Face->SetIntPoint(&ip);
      Vector c(3); ftr->Face->Transform(ip, c);
      if (std::abs(c(1)) > tol) { continue; }
      Array<int> verts; mesh.GetFaceVertices(f, verts);
      if (verts.Size() == 3)
      { mesh.AddBdrTriangle(verts[0], verts[1], verts[2], 3); }
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

static void BuildFaultCoords(ParMesh &pmesh,
                              const Array<int> &int_faces,
                              const Array<int> &shr_faces,
                              int order,
                              std::vector<Vector> &fault_coords)
{
   fault_coords.clear();
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr = pmesh.GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir =
         IntRules.Get(ftr->GetGeometryType(), 2*order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         ftr->SetAllIntPoints(&ir.IntPoint(q));
         Vector phys(3); ftr->Face->Transform(ir.IntPoint(q), phys);
         fault_coords.push_back(phys);
      }
   }
   for (int i = 0; i < shr_faces.Size(); i++)
   {
      auto *ftr = pmesh.GetSharedFaceTransformations(shr_faces[i]);
      const IntegrationRule &ir =
         IntRules.Get(ftr->GetGeometryType(), 2*order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         ftr->SetAllIntPoints(&ir.IntPoint(q));
         Vector phys(3); ftr->Face->Transform(ir.IntPoint(q), phys);
         fault_coords.push_back(phys);
      }
   }
}

int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);
   int rank = 0, nprocs = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   if (rank == 0)
   {
      std::cout << "\n=== TPV102 Pepper-Bug Test 5: Shared-Fault Pelties-7 "
                << "ANALYTIC Probe ===\n  MPI ranks: " << nprocs << "\n";
   }

   if (nprocs < 4)
   {
      if (rank == 0)
      { std::cout << "  SKIPPED: requires np >= 4 (ParMETIS cuts)\n"; }
      MPI_Finalize();
      return 77;
   }

   Mesh serial_mesh = BuildFixtureMesh();
   ParMesh pmesh(MPI_COMM_WORLD, serial_mesh);

   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {5};

   const int order = 1;
   WaveOperator<ParMesh> wave(pmesh, order,
                              TPV102Params::lambda,
                              TPV102Params::mu,
                              TPV102Params::rho, bc);

   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   const Array<int> &shr_faces = wave.GetFaultSharedFaces();
   const int n_int = int_faces.Size();
   const int n_shr = shr_faces.Size();

   int total_shared = 0;
   MPI_Allreduce(&n_shr, &total_shared, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if (rank == 0)
   {
      std::cout << "  total shared fault faces across ranks: "
                << total_shared << "\n";
   }
   if (total_shared == 0)
   {
      if (rank == 0) { std::cout << "  SKIPPED: no shared fault faces\n"; }
      MPI_Finalize();
      return 77;
   }

   int nqp_per_face = 0;
   if (n_shr > 0)
   {
      auto *ftr = pmesh.GetSharedFaceTransformations(shr_faces[0]);
      nqp_per_face = IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
   }
   else if (n_int > 0)
   {
      auto *ftr = pmesh.GetInteriorFaceTransformations(int_faces[0]);
      nqp_per_face = IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
   }
   int g_nqp = 0;
   MPI_Allreduce(&nqp_per_face, &g_nqp, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   nqp_per_face = g_nqp;

   const int num_fault_total = (n_int + n_shr) * nqp_per_face;

   std::vector<Vector> fault_coords;
   BuildFaultCoords(pmesh, int_faces, shr_faces, order, fault_coords);

   std::vector<DOFData> dof_data;
   if (num_fault_total > 0)
   {
      InitializeFaultDOFs(dof_data, num_fault_total, fault_coords);
      ZeroDOFDataPreStressTotal(dof_data, num_fault_total);
   }

   FaultFaceFlux fault_flux(TPV102Params::rho, TPV102Params::cp,
                            TPV102Params::cs);
   wave.SetFaultFlux(&fault_flux);
   wave.SetFaultDOFData(&dof_data, nqp_per_face);

   real_t bulk_bg[NUM_STATE] = {0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;
   wave.SetAbsorbingBackground(bulk_bg);

   const int size = wave.Height();
   const auto &fes_const = wave.GetFESpace();
   auto &fes = const_cast<ParFiniteElementSpace &>(fes_const);
   const int ndof_total = fes.GetNDofs();

   // Build Q with antisymmetric VX: +V_test on y>0 DOFs, -V_test on y<0 DOFs.
   // Mimics the anti-symmetric pattern from test_interior_fault_flux_path T2
   // but on the ParMesh.  Expected Pelties 7c at the fault QPs:
   //   v_t2⁺ = +V_test,  v_t2⁻ = -V_test,  v_t2⁻ - v_t2⁺ = -2V_test
   //   tau_2_trial = eta_s · (-2V_test) = -Zs · V_test
   const real_t V_test = 1.0e-5;
   const real_t Zs     = TPV102Params::rho * TPV102Params::cs;
   const real_t expected_dtau2 = -Zs * V_test;

   Vector Q(size);
   InitializeStateTotal(Q, ndof_total,
                        TPV102Params::sigma_n, TPV102Params::tau_ini);

   // Direct per-element antisymmetric VX assignment (mirrors the pattern
   // in test_interior_fault_flux_path T2).  Each element's VX DOFs get
   // +V_test if centroid y>0, -V_test if centroid y<0.  No DOF in this
   // mesh straddles y=0 (the 4-ny layer split is at y ∈ {-2000, -1000, 0,
   // 1000, 2000} m), so assignment is unambiguous.
   const int ndof_per_elem = fes.GetFE(0)->GetDof();
   real_t local_vx_plus = 0, local_vx_minus = 0;
   int local_plus_count = 0, local_minus_count = 0;
   for (int e = 0; e < pmesh.GetNE(); e++)
   {
      real_t cy = 0;
      Array<int> ev; pmesh.GetElementVertices(e, ev);
      for (int v = 0; v < ev.Size(); v++)
      { cy += pmesh.GetVertex(ev[v])[1]; }
      cy /= ev.Size();
      const real_t val = (cy > 0.0) ? +V_test : -V_test;
      if (cy > 0) { local_vx_plus = val; local_plus_count++; }
      else        { local_vx_minus = val; local_minus_count++; }
      Array<int> edofs; fes.GetElementDofs(e, edofs);
      for (int j = 0; j < ndof_per_elem; j++)
      {
         Q(VX * ndof_total + edofs[j]) += val;
      }
   }
   // Sanity print — confirm both signs are placed on this rank.
   int gp = 0, gm = 0;
   MPI_Reduce(&local_plus_count, &gp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(&local_minus_count, &gm, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
   if (rank == 0)
   {
      std::cout << "  global +side elems = " << gp
                << "  -side elems = " << gm << "\n";
      std::cout << "  per-rank 0 vx_plus,vx_minus = "
                << local_vx_plus << ", " << local_vx_minus << "\n";
   }

   Vector k(size);
   wave.Mult(Q, k);

   // Collect local tau2_corr, tau1_corr, sigma_n_corr on SHARED fault QPs.
   // Layout: dof_data = [interior QPs ..., shared QPs ...]
   int n_local_fail_tau2 = 0, n_local_fail_tau1 = 0;
   real_t worst_tau2_rel = 0.0;
   real_t worst_tau1_abs = 0.0;
   real_t local_avg_dtau2 = 0.0;
   int n_local_shared = 0;
   for (int i = 0; i < n_shr; i++)
   {
      for (int q = 0; q < nqp_per_face; q++)
      {
         const int dof_idx = (n_int + i) * nqp_per_face + q;
         if (dof_idx < 0 || dof_idx >= num_fault_total) { continue; }
         n_local_shared++;
         const DOFData &fd = dof_data[dof_idx];

         const real_t actual_dtau2 = fd.tau2_corr - TPV102Params::tau_ini;
         local_avg_dtau2 += actual_dtau2;
         const real_t err_tau2 = std::abs(actual_dtau2 - expected_dtau2);
         const real_t rel_tau2 = err_tau2 / std::abs(expected_dtau2);
         if (rel_tau2 > worst_tau2_rel) { worst_tau2_rel = rel_tau2; }
         if (rel_tau2 > 0.05)  // 5% tolerance (Pelties 7c)
         { n_local_fail_tau2++; }

         const real_t actual_tau1 = fd.tau1_corr;   // tau1_0 = 0 for TPV102
         const real_t abs_tau1 = std::abs(actual_tau1);
         if (abs_tau1 > worst_tau1_abs) { worst_tau1_abs = abs_tau1; }
         // tau1_corr should stay within friction-correction scale
         // 0.1 · eta_s · V_test ≈ 4.6 Pa
         if (abs_tau1 > 0.1 * 0.5 * Zs * V_test)
         { n_local_fail_tau1++; }
      }
   }

   int g_fail_tau2 = 0, g_fail_tau1 = 0, g_shared = 0;
   real_t g_sum_dtau2 = 0.0;
   MPI_Reduce(&local_avg_dtau2, &g_sum_dtau2, 1, MPI_DOUBLE, MPI_SUM, 0,
              MPI_COMM_WORLD);
   MPI_Reduce(&n_local_fail_tau2, &g_fail_tau2, 1, MPI_INT, MPI_SUM, 0,
              MPI_COMM_WORLD);
   MPI_Reduce(&n_local_fail_tau1, &g_fail_tau1, 1, MPI_INT, MPI_SUM, 0,
              MPI_COMM_WORLD);
   MPI_Reduce(&n_local_shared, &g_shared, 1, MPI_INT, MPI_SUM, 0,
              MPI_COMM_WORLD);
   real_t g_worst_tau2_rel = 0, g_worst_tau1_abs = 0;
   MPI_Reduce(&worst_tau2_rel, &g_worst_tau2_rel, 1, MPI_DOUBLE, MPI_MAX, 0,
              MPI_COMM_WORLD);
   MPI_Reduce(&worst_tau1_abs, &g_worst_tau1_abs, 1, MPI_DOUBLE, MPI_MAX, 0,
              MPI_COMM_WORLD);

   int exit_code = 0;
   if (rank == 0)
   {
      const real_t g_avg_dtau2 = (g_shared > 0) ? g_sum_dtau2 / g_shared : 0.0;
      std::cout << "  expected Δτ_2   = " << std::scientific
                << std::setprecision(3) << expected_dtau2 << " Pa\n";
      std::cout << "  actual avg Δτ_2 = " << g_avg_dtau2 << " Pa\n";
      std::cout << "  local-shared QPs (all ranks): " << g_shared << "\n";
      std::cout << "  worst rel err τ_2 = " << g_worst_tau2_rel << "\n";
      std::cout << "  worst |τ_1_corr|  = " << g_worst_tau1_abs << " Pa\n";

      std::cout << "\n  I1 (Pelties 7c shared-fault match):\n";
      if (g_fail_tau2 == 0)
      {
         std::cout << "    PASSED: 0 / " << g_shared
                   << " shared QPs outside 5% tol\n";
      }
      else
      {
         std::cout << "    FAILED: " << g_fail_tau2 << " / " << g_shared
                   << " shared QPs outside 5% tol (worst rel "
                   << g_worst_tau2_rel << ")\n";
         exit_code = 1;
      }
      std::cout << "  I2 (strike→dip leakage check):\n";
      if (g_fail_tau1 == 0)
      {
         std::cout << "    PASSED: 0 / " << g_shared
                   << " shared QPs with |τ_1_corr| > 0.1·η_s·V\n";
      }
      else
      {
         std::cout << "    FAILED: " << g_fail_tau1 << " / " << g_shared
                   << " shared QPs with |τ_1_corr| > 0.1·η_s·V "
                   << "(worst " << g_worst_tau1_abs << " Pa)\n";
         exit_code = 1;
      }
      std::cout << "\n========================================\n";
      std::cout << "  Test 5 shared-fault analytic: "
                << ((exit_code == 0) ? "PASS" : "FAIL") << "\n";
      std::cout << "========================================\n";
   }

   MPI_Bcast(&exit_code, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Finalize();
   return exit_code;
}

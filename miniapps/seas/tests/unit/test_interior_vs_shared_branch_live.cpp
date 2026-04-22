// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 pepper-bug unit-test plan (2026-04-22) — Phase 1 P1-T1:
// LIVE interior-vs-shared branch cross-check on the SAME physical face.
//
// ============================================================================
// Motivation
// ============================================================================
// `test_fault_flux_interior_vs_shared_branch_equivalence` is a tautology
// post-v9.0.0 Pelties-9 fix (its own header says "REPLACE this test
// with a live cross-check").  This file is the replacement.
//
// Strategy:
//   1. On rank 0, build a serial 2-tet mesh with interior fault at y=0.
//      Run wave.Mult(Q, k) — the INTERIOR fault branch fires.
//      Record k[c, dof_0_of_tet0] for c ∈ {SXX..VZ} into k_serial[9].
//   2. On all ranks, build the SAME serial mesh, construct a ParMesh via
//      explicit partitioning (rank 0 gets tet 0, rank 1 gets tet 1) so
//      the former interior fault face becomes a SHARED fault face.
//      Run wave.Mult on the ParMesh — the SHARED fault branch fires.
//      Record k[c, dof_0_of_rank0_tet] into k_parallel[9] on rank 0.
//   3. Assert |k_serial[c] - k_parallel[c]| ≤ 1e-10 × max(|k_serial[c]|, 1.0)
//      for every c.  Any failure confirms an interior-vs-shared
//      branch divergence on the same physical face — pepper signature.
//
// The test also demonstrates NON-tautology: with env var
// SEAS_TEST_TAMPER=1, it flips shared_fault_elem1_on_plus_[0] via a
// test-only mutator (added under SEAS_TEST_INTERNAL) BEFORE the
// parallel Mult.  This MUST fail the assertion.
//
// ============================================================================
// Usage:  mpirun -np 2 ./seas_test_interior_vs_shared_branch_live
//         (rank 0 runs BOTH the serial comparison and the parallel twin;
//          rank 1 only participates in the parallel half)

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../dynamic/tpv102_setup_total.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

static Mesh BuildTwoTetFaultMesh()
{
   Mesh mesh(3, 5, 2, 0);
   real_t verts[5][3] = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0},
      {0.0,  1.0, 0.0},   // V3 on +y side
      {0.0, -1.0, 0.0},   // V4 on -y side
   };
   for (int v = 0; v < 5; v++) { mesh.AddVertex(verts[v]); }
   mesh.AddTet(0, 1, 2, 4, 1);   // tet 0: − side
   mesh.AddTet(0, 1, 2, 3, 1);   // tet 1: + side
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

// Fill fault_coords / DOFData from the wave operator's fault-face lists.
template <typename MeshT>
static int SetupFault(WaveOperator<MeshT> &wave, MeshT &mesh, int order,
                      std::vector<DOFData> &dof_data, FaultFaceFlux &ff,
                      std::vector<Vector> &fault_coords)
{
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   const Array<int> &shr_faces = wave.GetFaultSharedFaces();

   int nqp = 0;
   if (int_faces.Size() > 0)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[0]);
      MFEM_VERIFY(ftr, "interior fault face FTR null");
      nqp = IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
   }
#ifdef MFEM_USE_MPI
   if constexpr (std::is_same_v<MeshT, ParMesh>)
   {
      if (nqp == 0 && shr_faces.Size() > 0)
      {
         auto *ftr = mesh.GetSharedFaceTransformations(shr_faces[0]);
         MFEM_VERIFY(ftr, "shared fault face FTR null");
         nqp = IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
      }
   }
#endif

   fault_coords.clear();
   for (int i = 0; i < int_faces.Size(); i++)
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
#ifdef MFEM_USE_MPI
   if constexpr (std::is_same_v<MeshT, ParMesh>)
   {
      for (int i = 0; i < shr_faces.Size(); i++)
      {
         auto *ftr = mesh.GetSharedFaceTransformations(shr_faces[i]);
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
   }
#endif

   const int n_fault =
      (int_faces.Size() + shr_faces.Size()) * nqp;
   if (n_fault > 0)
   {
      InitializeFaultDOFs(dof_data, n_fault, fault_coords);
      ZeroDOFDataPreStressTotal(dof_data, n_fault);
   }
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp);
   return n_fault;
}

// Set Q with asymmetric VX perturbation on the − side only (tet 0 or
// whatever element has centroid y < 0 on this (par)mesh).
template <typename MeshT>
static void SetOneSidedQ(MeshT &mesh, Vector &Q, int comp,
                         int ndof_total, const FiniteElementSpace &fes,
                         real_t val, int side_y_sign)
{
   for (int e = 0; e < mesh.GetNE(); e++)
   {
      real_t cy = 0;
      Array<int> ev; mesh.GetElementVertices(e, ev);
      for (int v = 0; v < ev.Size(); v++)
      { cy += mesh.GetVertex(ev[v])[1]; }
      cy /= ev.Size();
      const bool this_side = (side_y_sign < 0) ? (cy < 0.0) : (cy > 0.0);
      if (!this_side) { continue; }
      Array<int> edofs; fes.GetElementDofs(e, edofs);
      for (int j = 0; j < edofs.Size(); j++)
      {
         Q(comp * ndof_total + edofs[j]) += val;
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
      std::cout << "\n=== TPV102 Pepper-Bug Phase 1 P1-T1: "
                << "LIVE interior-vs-shared branch cross-check ===\n"
                << "  MPI ranks: " << nprocs << "\n";
   }

   if (nprocs < 2)
   {
      if (rank == 0)
      { std::cout << "  SKIPPED: requires np >= 2\n"; }
      MPI_Finalize();
      return 77;
   }

   const bool tamper =
      (std::getenv("SEAS_TEST_TAMPER") != nullptr) &&
      (std::string(std::getenv("SEAS_TEST_TAMPER")) == "1");
   if (rank == 0 && tamper)
   {
      std::cout << "  *** SEAS_TEST_TAMPER=1 set — will corrupt "
                << "shared_fault_elem1_on_plus_[0] to verify test is "
                << "non-tautological.  Expected: FAIL.\n";
   }

   const int order = 1;
   const real_t V_test = 1.0e-5;

   // Record k at the fault-adjacent VX DOF in the serial run (rank 0).
   real_t k_serial[NUM_STATE] = {0};
   if (rank == 0)
   {
      std::cout << "\n-- Serial (interior-fault branch) --\n";
      Mesh mesh = BuildTwoTetFaultMesh();
      BoundaryConfig bc;
      bc.natural_attrs   = {1};
      bc.fault_attr      = 3;
      bc.absorbing_attrs = {};
      WaveOperator<Mesh> wave(mesh, order,
                              TPV102Params::lambda,
                              TPV102Params::mu,
                              TPV102Params::rho, bc);
      real_t bulk_bg[NUM_STATE] = {0};
      bulk_bg[SYY] =  TPV102Params::sigma_n;
      bulk_bg[SXY] = -TPV102Params::tau_ini;
      wave.SetAbsorbingBackground(bulk_bg);

      FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
      std::vector<DOFData> dof_data;
      std::vector<Vector> fault_coords;
      SetupFault(wave, mesh, order, dof_data, ff, fault_coords);

      const auto &fes = wave.GetFESpace();
      const int ndof_total = fes.GetNDofs();
      const int size = wave.Height();

      Vector Q(size);
      InitializeStateTotal(Q, ndof_total,
                           TPV102Params::sigma_n, TPV102Params::tau_ini);
      SetOneSidedQ(mesh, Q, VX, ndof_total, fes, +V_test, /*side=*/-1);

      Vector k(size);
      wave.Mult(Q, k);

      // Record k[c, elem=0-tet-on-minus-side, dof=0] for each component.
      // Find the element whose centroid is y<0 (the − side).
      int e_minus = -1;
      for (int e = 0; e < mesh.GetNE(); e++)
      {
         real_t cy = 0;
         Array<int> ev; mesh.GetElementVertices(e, ev);
         for (int v = 0; v < ev.Size(); v++)
         { cy += mesh.GetVertex(ev[v])[1]; }
         cy /= ev.Size();
         if (cy < 0) { e_minus = e; break; }
      }
      MFEM_VERIFY(e_minus >= 0, "serial: no element with cy<0");
      Array<int> edofs; fes.GetElementDofs(e_minus, edofs);
      const int dof0 = edofs[0];
      for (int c = 0; c < NUM_STATE; c++)
      {
         k_serial[c] = k(c * ndof_total + dof0);
      }
      std::cout << "    serial e_minus=" << e_minus << " dof0=" << dof0 << "\n";
      for (int c = 0; c < NUM_STATE; c++)
      {
         std::cout << "    k_serial[" << c << "] = "
                   << std::scientific << std::setprecision(6)
                   << k_serial[c] << "\n";
      }
   }

   MPI_Bcast(k_serial, NUM_STATE, MPI_DOUBLE, 0, MPI_COMM_WORLD);

   // Parallel: ParMesh with one-tet-per-rank (explicit partitioning).
   if (rank == 0)
   {
      std::cout << "\n-- Parallel (shared-fault branch, np=2, explicit partition) --\n";
   }

   Mesh serial_mesh = BuildTwoTetFaultMesh();
   // Explicit partitioning: 2 tets → tet 0 to rank 0, tet 1 to rank 1.
   // Only meaningful when nprocs=2; for nprocs > 2 we skip because the
   // partitioning array has size == serial_mesh.GetNE() == 2.
   if (nprocs != 2)
   {
      if (rank == 0)
      {
         std::cout << "  SKIPPED: parallel branch requires np == 2 "
                   << "(fixture is 2 elements)\n";
      }
      MPI_Finalize();
      return 77;
   }
   std::vector<int> partitioning(serial_mesh.GetNE());
   // tet 0 (has V4 at y=-1 → − side) → rank 0
   // tet 1 (has V3 at y=+1 → + side) → rank 1
   for (int e = 0; e < serial_mesh.GetNE(); e++)
   {
      real_t cy = 0;
      Array<int> ev; serial_mesh.GetElementVertices(e, ev);
      for (int v = 0; v < ev.Size(); v++)
      { cy += serial_mesh.GetVertex(ev[v])[1]; }
      cy /= ev.Size();
      partitioning[e] = (cy < 0) ? 0 : 1;
   }
   ParMesh pmesh(MPI_COMM_WORLD, serial_mesh, partitioning.data());

   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {};
   WaveOperator<ParMesh> wave(pmesh, order,
                              TPV102Params::lambda,
                              TPV102Params::mu,
                              TPV102Params::rho, bc);
   real_t bulk_bg[NUM_STATE] = {0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;
   wave.SetAbsorbingBackground(bulk_bg);

#ifdef SEAS_TEST_INTERNAL
   if (tamper)
   {
      // Only tamper on the rank that owns at least one shared fault face.
      if (wave.GetFaultSharedFaces().Size() > 0)
      {
         std::cout << "  [rank " << rank
                   << "] TAMPERING: flipping shared_fault_elem1_on_plus_[0]\n";
         wave.TamperSharedFaultElem1OnPlus(0);
      }
   }
#else
#  error "This test requires -DSEAS_TEST_INTERNAL"
#endif

   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   std::vector<DOFData> dof_data;
   std::vector<Vector> fault_coords;
   SetupFault(wave, pmesh, order, dof_data, ff, fault_coords);

   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   const int size = wave.Height();

   Vector Q(size);
   InitializeStateTotal(Q, ndof_total,
                        TPV102Params::sigma_n, TPV102Params::tau_ini);
   SetOneSidedQ(pmesh, Q, VX, ndof_total, fes, +V_test, /*side=*/-1);

   Vector k(size);
   wave.Mult(Q, k);

   // Rank 0 should own exactly one element (the − side tet).
   real_t k_parallel[NUM_STATE] = {0};
   if (rank == 0)
   {
      const int ne_local = pmesh.GetNE();
      std::cout << "    rank 0 ne_local=" << ne_local << "\n";
      if (ne_local >= 1)
      {
         // Find the − side element on rank 0 (should be the only one).
         int e_minus = -1;
         for (int e = 0; e < ne_local; e++)
         {
            real_t cy = 0;
            Array<int> ev; pmesh.GetElementVertices(e, ev);
            for (int v = 0; v < ev.Size(); v++)
            { cy += pmesh.GetVertex(ev[v])[1]; }
            cy /= ev.Size();
            if (cy < 0) { e_minus = e; break; }
         }
         MFEM_VERIFY(e_minus >= 0, "parallel rank 0: no element with cy<0");
         Array<int> edofs; fes.GetElementDofs(e_minus, edofs);
         const int dof0 = edofs[0];
         for (int c = 0; c < NUM_STATE; c++)
         {
            k_parallel[c] = k(c * ndof_total + dof0);
         }
         std::cout << "    parallel e_minus=" << e_minus
                   << " dof0=" << dof0 << "\n";
         for (int c = 0; c < NUM_STATE; c++)
         {
            std::cout << "    k_parallel[" << c << "] = "
                      << std::scientific << std::setprecision(6)
                      << k_parallel[c] << "\n";
         }
      }
   }

   // Compare on rank 0.
   int exit_code = 0;
   if (rank == 0)
   {
      std::cout << "\n-- Comparison --\n";
      int n_fail = 0;
      real_t worst_rel = 0.0;
      for (int c = 0; c < NUM_STATE; c++)
      {
         const real_t scale =
            std::max({std::abs(k_serial[c]), std::abs(k_parallel[c]),
                      real_t(1.0)});
         const real_t d = std::abs(k_serial[c] - k_parallel[c]);
         const real_t r = d / scale;
         if (r > worst_rel) { worst_rel = r; }
         if (r > 1.0e-10)
         {
            n_fail++;
            std::cout << "    FAIL c=" << c
                      << "  serial=" << std::scientific
                      << std::setprecision(6) << k_serial[c]
                      << "  parallel=" << k_parallel[c]
                      << "  rel=" << r << "\n";
         }
      }
      std::cout << "  worst relative |k_serial - k_parallel| = "
                << worst_rel << " over all 9 components\n";
      std::cout << "  components failing rel>1e-10 : " << n_fail << "\n";
      const bool expected_fail = tamper;
      if (expected_fail)
      {
         if (n_fail > 0)
         {
            std::cout << "  PASSED (tamper mode): test DETECTED the "
                      << "injected bug\n";
         }
         else
         {
            std::cout << "  FAILED (tamper mode): test is TAUTOLOGICAL — "
                      << "tamper did not break anything\n";
            exit_code = 1;
         }
      }
      else
      {
         if (n_fail == 0)
         {
            std::cout << "  PASSED: interior branch matches shared "
                      << "branch on the same physical face\n";
         }
         else
         {
            std::cout << "  FAILED: interior-vs-shared branch divergence "
                      << "— PEPPER BUG SIGNATURE\n";
            exit_code = 1;
         }
      }
      std::cout << "========================================\n";
   }

   MPI_Bcast(&exit_code, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Finalize();
   return exit_code;
}

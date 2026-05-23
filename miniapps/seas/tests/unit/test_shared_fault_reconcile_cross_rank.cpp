// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Phase 1 of PLAN_shared_fault_reconcile_fix_2026-05-23.md — the deterministic
// local RED->GREEN oracle for the shared-fault cross-rank reconcile.
//
// ============================================================================
// What this test proves
// ============================================================================
// Phase 0 (Frontera job 7747304) localised the SAFS cross-rank divergence to
// SOURCE (1) INTERPOLATION: on a shared fault face the two ranks feed the
// friction solve inputs that differ at ~1e-14 (the same physical QP reached
// through `shape1.Loc1` on the owner vs `shape2.Loc2` on the ghost).  At the
// LSW slip-onset kink that 1e-14 tips one rank to slip while the other locks
// (worst_rel -> 1).  Rate-and-state stays bounded (~1e-14) only because its
// solve is smooth — the SAME inconsistency is present, just not amplified.
//
// The clean 2-tet planar fixture is naturally bit-identical across ranks
// (shape1==shape2 by symmetry), so it cannot make the seed on its own.  We
// therefore INJECT a deterministic 1-ULP cross-rank difference on ONE rank via
// the SEAS_TEST_INTERNAL hook in FaultFaceFlux::ComputeTrialTraction, then run
// the np=2 shared-fault path and check WaveOperator's own R-101 consistency
// verifier (the exact Frontera check) every step.
//
// Oracle (per friction law):
//   - BEFORE the Phase-2 reconcile: worst_rel != 0 every-step max  ->  RED
//        (LSW -> O(1) once the kink is reached; rate-state ~1e-14 — both nonzero,
//         proving the defect is METHOD-INVARIANT, not LSW-specific).
//   - AFTER the Phase-2 reconcile: worst_rel == 0 (bit-identical across ranks)
//        for BOTH laws  ->  GREEN.
//
// The assertion is `max_worst_rel == 0` (exact cross-rank bit-identity), which
// is the post-fix contract from the plan (physically-exact AND cross-rank
// bit-identical).  It is RED today because the injection makes the two ranks
// differ; it goes GREEN when the reconcile copies the owner's state to the peer.
//
// Usage:  mpirun -np 2 ./seas_test_shared_fault_reconcile_cross_rank
//         Optional: SEAS_RECONCILE_NSTEPS=<N> (default 200).
//         Requires MFEM built with MPI and -DSEAS_TEST_INTERNAL.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../dynamic/tpv102_setup_total.hpp"
#include "../../dynamic/tpv205_setup.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../config/tpv205_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_TRUE(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; \
      std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

namespace {

constexpr real_t kL = 1000.0;
constexpr real_t kDt = 1.0e-4;
constexpr int    kOrder = 1;

// Planar 2-tet fault (y=0 plane); identical to the rake-sweep / multistep tests.
Mesh BuildTwoTetFaultMesh()
{
   Mesh mesh(3, 5, 2, 0);
   real_t verts[5][3] = {
      {0.0, 0.0, 0.0}, {kL, 0.0, 0.0}, {0.0, 0.0, kL},
      {0.0,  kL, 0.0},   // V3 on +y side
      {0.0, -kL, 0.0},   // V4 on -y side
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
         for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
         cy /= fv.Size();
         if (std::abs(cy) < 1e-8)
         { mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3); }  // fault attr 3
         continue;
      }
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 1);         // natural attr 1
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// Gather the physical coords of every fault QP (interior + shared) in the same
// order InitializeFaultDOFs expects.  Mirrors the rake-sweep helper.
template <typename MeshT>
int CollectFaultQPCoords(WaveOperator<MeshT> &wave, MeshT &mesh, int order,
                         std::vector<Vector> &coords)
{
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   const Array<int> &shr_faces = wave.GetFaultSharedFaces();
   int nqp = 0;
   if (int_faces.Size() > 0)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[0]);
      MFEM_VERIFY(ftr, "interior fault FTR null");
      nqp = IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
   }
#ifdef MFEM_USE_MPI
   if constexpr (std::is_same_v<MeshT, ParMesh>)
   {
      if (nqp == 0 && shr_faces.Size() > 0)
      {
         auto *ftr = mesh.GetSharedFaceTransformations(shr_faces[0]);
         MFEM_VERIFY(ftr, "shared fault FTR null");
         nqp = IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
      }
   }
#endif
   coords.clear();
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         coords.push_back(phys);
      }
   }
#ifdef MFEM_USE_MPI
   if constexpr (std::is_same_v<MeshT, ParMesh>)
   {
      for (int i = 0; i < shr_faces.Size(); i++)
      {
         auto *ftr = mesh.GetSharedFaceTransformations(shr_faces[i]);
         const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                                  2*order);
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            ftr->SetAllIntPoints(&ip);
            Vector phys(3); ftr->Face->Transform(ip, phys);
            coords.push_back(phys);
         }
      }
   }
#endif
   return (int_faces.Size() + shr_faces.Size()) * nqp;
}

ParMesh MakePartitioned2Tet()
{
   Mesh serial_mesh = BuildTwoTetFaultMesh();
   std::vector<int> part(serial_mesh.GetNE(), 0);
   for (int e = 0; e < serial_mesh.GetNE(); e++)
   {
      real_t cy = 0; Array<int> ev; serial_mesh.GetElementVertices(e, ev);
      for (int v = 0; v < ev.Size(); v++) { cy += serial_mesh.GetVertex(ev[v])[1]; }
      cy /= ev.Size();
      part[e] = (cy < 0.0) ? 0 : 1;   // −side→rank0, +side→rank1
   }
   return ParMesh(MPI_COMM_WORLD, serial_mesh, part.data());
}

// Run the np=2 shared-fault path for one friction law with the 1-ULP seed
// injected on rank 0; return the every-step max of the R-101 worst_rel and
// whether the fault ever slipped (onset sanity).
struct LegResult { double max_wr = 0.0; double max_slip_rate = 0.0; int onset_steps = 0; };

LegResult RunLeg(FaultFrictionLaw law, int rank, int nsteps,
                 bool disable_reconcile = false)
{
   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr = 3;
   bc.absorbing_attrs = {};

   ParMesh pmesh = MakePartitioned2Tet();

   const bool lsw = (law == FaultFrictionLaw::LSW);
   const real_t lambda = lsw ? TPV205Params::lambda : TPV102Params::lambda;
   const real_t mu     = lsw ? TPV205Params::mu     : TPV102Params::mu;
   const real_t rho    = lsw ? TPV205Params::rho    : TPV102Params::rho;
   const real_t cp     = lsw ? TPV205Params::cp     : TPV102Params::cp;
   const real_t cs     = lsw ? TPV205Params::cs     : TPV102Params::cs;

   WaveOperator<ParMesh> wave(pmesh, kOrder, lambda, mu, rho, bc);
   wave.SetFaultFrictionLaw(law);

   FaultFaceFlux ff(rho, cp, cs);
   wave.SetFaultFlux(&ff);

   std::vector<Vector> coords;
   const int n_fault = CollectFaultQPCoords(wave, pmesh, kOrder, coords);
   int nqp = 0;
   {
      const Array<int> &ifc = wave.GetFaultInteriorFaces();
      const Array<int> &sfc = wave.GetFaultSharedFaces();
      const int nf = ifc.Size() + sfc.Size();
      if (nf > 0 && n_fault > 0) { nqp = n_fault / nf; }
   }

   std::vector<DOFData> dof_data;
   Vector Q(wave.Height()), Q_new(wave.Height());
   Q = 0.0;
   const int ndof_total = wave.GetFESpace().GetNDofs();

   if (lsw)
   {
      // TPV205 (LSW) — fluctuation-Q with a DETERMINISTIC prestress knife's
      // edge: put the shear entirely in the (exact, roundoff-free) prestress
      // so |tau| == tau_str EXACTLY, with bulk Q = 0 so the unperturbed trial
      // is exactly 0.  Then the unperturbed rank is exactly LOCKED (V=0) and
      // the seed-injected rank crosses the kink (V>0) => worst_rel -> 1.
      // mu_d == mu_s (no weakening) keeps tau_str fixed so the edge is stable.
      // mu_s * sigma_n = 0.5 * 120e6 = 60e6 is exact in binary.
      real_t bulk_bg[NUM_STATE] = {0};
      wave.SetAbsorbingBackground(bulk_bg);
      if (n_fault > 0)
      {
         InitializeFaultDOFs_TPV205(dof_data, n_fault, coords);
         const real_t mu_s = 0.5, sig_n = 120.0e6;
         for (int i = 0; i < n_fault; i++)
         {
            dof_data[i].lsw_mu_s = mu_s;
            dof_data[i].lsw_mu_d = mu_s;          // no weakening: tau_str fixed
            dof_data[i].lsw_d_c  = 0.40;
            dof_data[i].sigma_n0 = sig_n;
            dof_data[i].tau1_0   = 0.0;
            dof_data[i].tau2_0   = mu_s * sig_n;  // exactly at the kink
            dof_data[i].sigma_n_nuc = 0.0;
            dof_data[i].tau1_nuc = 0.0;
            dof_data[i].tau2_nuc = 0.0;
         }
      }
   }
   else
   {
      // TPV102 (rate-and-state) — total-Q (prestress in bulk Q).
      real_t bulk_bg[NUM_STATE] = {0};
      bulk_bg[SYY] =  TPV102Params::sigma_n;
      bulk_bg[SXY] = -TPV102Params::tau_ini;
      wave.SetAbsorbingBackground(bulk_bg);
      if (n_fault > 0)
      {
         InitializeFaultDOFs(dof_data, n_fault, coords);
         ZeroDOFDataPreStressTotal(dof_data, n_fault);
         for (int i = 0; i < n_fault; i++)
         { dof_data[i].tau2_nuc = TPV102Params::nuc_dtau; }  // persistent driver
      }
      InitializeStateTotal(Q, ndof_total, TPV102Params::sigma_n,
                           TPV102Params::tau_ini);
   }

   wave.SetFaultDOFData(&dof_data, nqp);

   // Inject the cross-rank strike-traction difference on rank 0 ONLY.
   //  - rate-state: the FAITHFUL seed scale (~1e-6 Pa = ~1.6e-14 of the 26 MPa
   //    trial — the real shared-face interpolation gap).  Smooth friction maps
   //    it to a bounded ~1e-14 worst_rel (defect present, not amplified).
   //  - LSW: injected at the AMPLIFIED scale (~10 MPa) so the injected rank
   //    crosses the slip-onset kink into PHYSICAL slip (~2 m/s) while the
   //    unperturbed rank stays exactly locked (V=0) — directly reproducing the
   //    SAFS failure mode (one slips / one locks => worst_rel -> 1).  On the
   //    real mesh this O(1) split EMERGES as slip-weakening amplifies the
   //    1e-14 seed over the driven nucleation (at onset the ranks' |tau|
   //    already differ by ~7 MPa, job 7747304); a static 2-tet has no
   //    nucleation pump to grow the seed, so we inject at the amplified scale.
   //    The reconcile must zero the split at EITHER scale.
   const real_t inject_pa = lsw ? 1.0e7 : 1.0e-6;
#ifdef SEAS_TEST_INTERNAL
   FaultFaceFlux::s_seas_test_tau2_trial_perturb_pa = (rank == 0) ? inject_pa : 0.0;
   // Set IDENTICALLY on all ranks (the reconcile's MPI exchange is collective).
   FaultFaceFlux::s_seas_test_disable_reconcile = disable_reconcile;
#endif
   (void)inject_pa;
   (void)disable_reconcile;

   LegResult res;
   for (int step = 0; step < nsteps; step++)
   {
      wave.AdvanceADER(Q, kDt, /*ader_order=*/2, Q_new);
      Q.Swap(Q_new);

      double wr = 0.0; int wf = -1;
      wave.VerifySharedFaultDOFDataConsistency(1.0e-10, &wr, &wf,
                                               /*abort_on_fail=*/false);
      res.max_wr = std::max(res.max_wr, wr);

      for (const auto &d : dof_data)
      {
         res.max_slip_rate = std::max(res.max_slip_rate,
                                      std::abs((double)d.slip_rate));
         if (std::abs((double)d.slip_rate) > 0.0) { res.onset_steps++; break; }
      }
   }
   return res;
}

} // anonymous

int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);
   int rank = 0, nprocs = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   int nsteps = 200;
   if (const char *e = std::getenv("SEAS_RECONCILE_NSTEPS"))
   { int v = std::atoi(e); if (v > 0) { nsteps = v; } }

   if (rank == 0)
   {
      std::cout << "\n=== shared-fault cross-rank reconcile oracle (Phase 1) ==="
                << "\n  np=" << nprocs << ", planar 2-tet fault, dt="
                << std::scientific << std::setprecision(2) << kDt
                << ", N=" << nsteps << " ADER-2 steps"
                << "\n  injecting a 1-ULP strike-traction seed on rank 0 "
                   "(SEAS_TEST_INTERNAL)\n";
   }
   if (nprocs != 2)
   {
      if (rank == 0) { std::cout << "  SKIPPED: requires np == 2\n"; }
      MPI_Finalize();
      return 77;
   }
#ifndef SEAS_TEST_INTERNAL
   if (rank == 0)
   {
      std::cout << "  SKIPPED: built without -DSEAS_TEST_INTERNAL "
                   "(injection hook compiled out)\n";
   }
   MPI_Finalize();
   return 77;
#else
   struct { FaultFrictionLaw law; const char *name; } legs[] = {
      {FaultFrictionLaw::RateAndState, "rate-and-state (TPV102)"},
      {FaultFrictionLaw::LSW,          "linear slip-weakening (TPV205)"},
   };

   for (const auto &leg : legs)
   {
      LegResult r = RunLeg(leg.law, rank, nsteps);

      // worst_rel is already an MPI_MAX inside the verifier (rank-0 holds the
      // global value); reduce our running max to be safe.
      double g_wr = r.max_wr;
      MPI_Reduce(&r.max_wr, &g_wr, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

      if (rank == 0)
      {
         std::cout << "\n  [" << leg.name << "]"
                   << "  max worst_rel = " << std::scientific
                   << std::setprecision(6) << g_wr
                   << "  (slipped: " << (r.onset_steps > 0 ? "yes" : "no")
                   << ", max|V|=" << r.max_slip_rate << ")\n";
         // GREEN (post Phase-2 reconcile): cross-rank bit-identical.
         TEST_TRUE(g_wr == 0.0,
                   std::string("cross-rank DOFData bit-identical despite "
                               "injected strike-traction seed — ") + leg.name);
      }
   }

   // Phase-3 NEGATIVE leg: disable the reconcile -> the R-101 guard must TRIP
   // (worst_rel > tol), proving the guard is not a no-op and the reconcile is
   // what makes the ranks consistent.  LSW (the O(1) split) is the sharp case.
   {
      LegResult r = RunLeg(FaultFrictionLaw::LSW, rank, nsteps,
                           /*disable_reconcile=*/true);
      double g_wr = r.max_wr;
      MPI_Reduce(&r.max_wr, &g_wr, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
      if (rank == 0)
      {
         std::cout << "\n  [NEGATIVE: LSW, reconcile DISABLED]"
                   << "  max worst_rel = " << std::scientific
                   << std::setprecision(6) << g_wr << "\n";
         TEST_TRUE(g_wr > 1.0e-10,
                   "R-101 guard TRIPS with the reconcile disabled (proves the "
                   "guard catches a real cross-rank desync)");
      }
   }

   if (rank == 0)
   {
      std::cout << "\n========================================\n"
                << "  Results: " << num_passed << " passed, "
                << num_failed << " failed out of " << num_tests << " tests\n"
                << "========================================\n";
   }
   int exit_code = (num_failed == 0) ? 0 : 1;
   MPI_Bcast(&exit_code, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Finalize();
   return exit_code;
#endif
}

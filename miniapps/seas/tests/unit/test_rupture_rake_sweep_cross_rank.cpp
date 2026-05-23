// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// GENERAL rake-sweep reproducer for the weak-channel tangential-traction
// instability (cross-rank V1/V2 divergence).  Phase 1 of
// PLAN_dip_channel_v1_cross_rank_stability_2026-05-22.md.
//
// ============================================================================
// Motivation
// ============================================================================
// The SAFS run trips R-101 on a cross-rank divergence of V1 (dip slip-rate) at
// a shared fault QP.  The hypothesis (H-B) is that this is NOT dip-specific:
// the WEAK tangential channel (the slip component small relative to |τ|) is
// ill-conditioned (its rake fraction τk/|τ| amplifies a 1-ULP wiggle) AND has
// a closed-loop gain > 1, so it amplifies the unavoidable cross-rank FP seed
// (the self/nbr swap in the shared-face rotation) to O(1).  At the SAFS rake
// the dip channel is weak; under a dip-dominated load the STRIKE channel would
// be the unstable one.
//
// This test drives the SAME 2-tet planar TPV102 fault (so the tilt/side effect
// is removed and only the rake matters) with an OBLIQUE nucleation rake φ:
//
//     tau1_nuc = Δτ·sin φ   (dip)        tau2_nuc = Δτ·cos φ   (strike)
//
//   φ =  0°  → pure strike  (dip   = weak channel)  == the proven-stable TPV path
//   φ = 90°  → pure dip     (strike = weak channel)
//   φ = 45°  → balanced
//
// and MONITORS the cross-rank consistency each step via the non-aborting mode
// of WaveOperator::VerifySharedFaultDOFDataConsistency (rank0-vs-rank1 — the
// exact Frontera R-101 check).  It also tracks max|V1|,max|V2| and a serial
// boundedness leg.
//
// Interpretation (the headline):
//   - φ=0 (pure strike) cross-rank-stable AND φ>0 grows / exceeds 1e-10
//     AND the growth is NOT symmetric about φ=45°  ⇒  NON-COVARIANT
//     weak-channel instability CONFIRMED (H-B).  Phase 3 fix required.
//   - everything bit-consistent for all φ  ⇒  already covariant (re-examine).
//   - per-rank τk_0/τk_nuc differ at step 0  ⇒  consistency gap (H-A/H-C).
//
// This test is RED before the Phase-3 fix (by design) and GREEN after.
//
// Usage:  mpirun -np 2 ./seas_test_rupture_rake_sweep_cross_rank
//         Optional: SEAS_RAKE_NSTEPS=<N> to override the step count.
//         Requires MFEM built with MPI.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../dynamic/tpv102_setup_total.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_LE(v, tol, msg) do { \
   num_tests++; \
   double vv = (v), tt = (tol); \
   if (vv <= tt) { num_passed++; \
      std::cout << "  PASSED: " << msg << "  (got " << std::scientific \
                << std::setprecision(3) << vv << ", tol " << tt << ")\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << "  (got " << std::scientific << std::setprecision(3) \
                << vv << ", expected <= " << tt << ")\n"; } \
} while (0)

namespace {

constexpr real_t kL = 1000.0;
constexpr real_t kDt = 1.0e-4;
constexpr int    kOrder = 1;

// Planar 2-tet TPV102 fault (y=0 plane); identical to the multistep test.
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

// Walk the wave operator's fault faces, build DOFData, and set an OBLIQUE
// persistent nucleation at rake φ.  Mirrors the multistep SetupFault but
// splits Δτ between the dip (tau1_nuc) and strike (tau2_nuc) channels.
template <typename MeshT>
int SetupFaultOblique(WaveOperator<MeshT> &wave, MeshT &mesh, int order,
                      std::vector<DOFData> &dof_data, FaultFaceFlux &ff,
                      std::vector<Vector> &fault_coords, real_t rake_deg)
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

   fault_coords.clear();
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*order);
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
         const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                                  2*order);
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

   const int n_fault = (int_faces.Size() + shr_faces.Size()) * nqp;
   if (n_fault > 0)
   {
      InitializeFaultDOFs(dof_data, n_fault, fault_coords);
      ZeroDOFDataPreStressTotal(dof_data, n_fault);
      const real_t r    = rake_deg * M_PI / 180.0;
      const real_t dtau = TPV102Params::nuc_dtau;
      for (int i = 0; i < n_fault; i++)
      {
         dof_data[i].tau1_nuc = dtau * std::sin(r);  // dip
         dof_data[i].tau2_nuc = dtau * std::cos(r);  // strike
      }
   }
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp);
   return n_fault;
}

void LocalMaxV(const std::vector<DOFData> &d, real_t &mv1, real_t &mv2,
               real_t &mvabs)
{
   mv1 = mv2 = mvabs = 0.0;
   for (const auto &x : d)
   {
      mv1   = std::max(mv1,   std::abs(x.V1));
      mv2   = std::max(mv2,   std::abs(x.V2));
      mvabs = std::max(mvabs, std::abs(x.slip_rate));
   }
}

} // anonymous

int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);
   int rank = 0, nprocs = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   int nsteps = 400;
   if (const char *e = std::getenv("SEAS_RAKE_NSTEPS"))
   {
      int v = std::atoi(e);
      if (v > 0) { nsteps = v; }
   }

   if (rank == 0)
   {
      std::cout << "\n=== GENERAL rake-sweep weak-channel cross-rank "
                << "reproducer ===\n"
                << "  np = " << nprocs << ", planar 2-tet TPV102 fault, dt = "
                << std::scientific << std::setprecision(2) << kDt
                << " s, N = " << nsteps << " ADER-2 steps\n"
                << "  oblique nucleation: tau1_nuc=Δτ·sinφ (dip), "
                << "tau2_nuc=Δτ·cosφ (strike); Δτ=" << TPV102Params::nuc_dtau
                << " Pa\n"
                << "  cross-rank check = VerifySharedFaultDOFDataConsistency "
                << "(non-aborting); φ=0 is the proven-stable TPV anchor\n";
   }
   if (nprocs != 2)
   {
      if (rank == 0) { std::cout << "  SKIPPED: requires np == 2\n"; }
      MPI_Finalize();
      return 77;
   }

   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr = 3;
   bc.absorbing_attrs = {};
   real_t bulk_bg[NUM_STATE] = {0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;

   const real_t rakes[] = {0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0};
   const int n_rake = static_cast<int>(sizeof(rakes)/sizeof(rakes[0]));

   const int chk[3] = {std::min(50, nsteps-1), std::min(150, nsteps-1),
                       nsteps-1};

   if (rank == 0)
   {
      std::cout << "\n" << std::setw(7) << "rake"
                << std::setw(8) << "weakch"
                << std::setw(13) << "maxV1"
                << std::setw(13) << "maxV2"
                << std::setw(13) << ("wr@" + std::to_string(chk[0]))
                << std::setw(13) << ("wr@" + std::to_string(chk[1]))
                << std::setw(13) << ("wr@" + std::to_string(chk[2])) << "\n";
   }

   real_t pure_strike_wr = 0.0;     // rake = 0 (TPV anchor) final worst_rel
   real_t worst_oblique_wr = 0.0;   // max over rake > 0 of final worst_rel
   real_t any_blowup = 0.0;         // max |V| over the whole sweep (NaN guard)

   for (int ri = 0; ri < n_rake; ri++)
   {
      const real_t phi = rakes[ri];

      Mesh serial_mesh = BuildTwoTetFaultMesh();
      std::vector<int> part(serial_mesh.GetNE(), 0);
      for (int e = 0; e < serial_mesh.GetNE(); e++)
      {
         real_t cy = 0; Array<int> ev; serial_mesh.GetElementVertices(e, ev);
         for (int v = 0; v < ev.Size(); v++)
         { cy += serial_mesh.GetVertex(ev[v])[1]; }
         cy /= ev.Size();
         part[e] = (cy < 0.0) ? 0 : 1;   // −side→rank0, +side→rank1
      }
      ParMesh pmesh(MPI_COMM_WORLD, serial_mesh, part.data());

      WaveOperator<ParMesh> wave_p(pmesh, kOrder, TPV102Params::lambda,
                                   TPV102Params::mu, TPV102Params::rho, bc);
      wave_p.SetAbsorbingBackground(bulk_bg);
      FaultFaceFlux ff_p(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
      std::vector<DOFData> dof_p;
      std::vector<Vector> fc_p;
      SetupFaultOblique(wave_p, pmesh, kOrder, dof_p, ff_p, fc_p, phi);

      const int ndof_total_p = wave_p.GetFESpace().GetNDofs();
      Vector Q_p(wave_p.Height()), Q_p_new(wave_p.Height());
      InitializeStateTotal(Q_p, ndof_total_p, TPV102Params::sigma_n,
                           TPV102Params::tau_ini);

      std::vector<double> wr_hist(nsteps, 0.0);
      real_t maxV1 = 0.0, maxV2 = 0.0;
      for (int step = 0; step < nsteps; step++)
      {
         wave_p.AdvanceADER(Q_p, kDt, /*ader_order=*/2, Q_p_new);
         Q_p.Swap(Q_p_new);

         double wr = 0.0; int wf = -1;
         wave_p.VerifySharedFaultDOFDataConsistency(1.0e-10, &wr, &wf,
                                                    /*abort_on_fail=*/false);
         wr_hist[step] = wr;

         real_t mv1, mv2, mvabs; LocalMaxV(dof_p, mv1, mv2, mvabs);
         maxV1 = std::max(maxV1, mv1);
         maxV2 = std::max(maxV2, mv2);
      }

      // Reduce per-rank maxima to rank 0 for reporting.
      real_t gV1 = 0.0, gV2 = 0.0;
      MPI_Reduce(&maxV1, &gV1, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, 0,
                 MPI_COMM_WORLD);
      MPI_Reduce(&maxV2, &gV2, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, 0,
                 MPI_COMM_WORLD);

      const double wr_final = wr_hist[chk[2]];
      if (phi == 0.0) { pure_strike_wr = wr_final; }
      else { worst_oblique_wr = std::max(worst_oblique_wr, wr_final); }
      any_blowup = std::max({any_blowup, gV1, gV2});

      if (rank == 0)
      {
         const char *weak = (phi < 45.0) ? "dip(V1)"
                          : (phi > 45.0) ? "strike(V2)" : "balanced";
         std::cout << std::setw(6) << std::fixed << std::setprecision(0) << phi
                   << std::setw(11) << weak
                   << std::setw(13) << std::scientific << std::setprecision(3)
                   << gV1 << std::setw(13) << gV2
                   << std::setw(13) << wr_hist[chk[0]]
                   << std::setw(13) << wr_hist[chk[1]]
                   << std::setw(13) << wr_hist[chk[2]] << "\n";
      }
   }

   if (rank == 0)
   {
      std::cout << "\n  pure-strike (φ=0) final worst_rel = "
                << std::scientific << std::setprecision(6) << pure_strike_wr
                << "\n  worst oblique  (φ>0) final worst_rel = "
                << worst_oblique_wr << "\n";

      const bool stable_anchor   = (pure_strike_wr <= 1.0e-10);
      const bool oblique_diverges = (worst_oblique_wr > 1.0e-10);
      std::cout << "\n  VERDICT: ";
      if (std::isnan(any_blowup) || any_blowup > 1.0e6)
      {
         std::cout << "blow-up (|V| non-physical) — instability is "
                      "destructive, not just a consistency trip.\n";
      }
      else if (stable_anchor && oblique_diverges)
      {
         std::cout << "NON-COVARIANT weak-channel instability CONFIRMED (H-B): "
                      "φ=0 stable, φ>0 diverges cross-rank.  Phase-3 covariant "
                      "fix required.\n";
      }
      else if (!stable_anchor)
      {
         std::cout << "even pure strike diverges — re-examine the cross-rank "
                      "seed / consistency (H-A/H-C) before H-B.\n";
      }
      else
      {
         std::cout << "all rakes cross-rank-consistent — scheme already "
                      "covariant on this fixture (re-examine reproduction).\n";
      }

      // Anchor: pure strike (the proven-stable TPV regime) MUST stay
      // cross-rank-consistent.  (Sanity — should always pass.)
      TEST_LE(pure_strike_wr, 1.0e-10,
              "pure-strike rake=0 cross-rank consistent (TPV anchor)");
      // CAPTURED BUG (RED until the Phase-3 covariant fix): every rake — not
      // just pure strike — must be cross-rank-consistent.
      TEST_LE(worst_oblique_wr, 1.0e-10,
              "oblique-rake cross-rank V1/V2 consistent at ALL rakes "
              "(RED until Phase-3 covariant fix)");

      std::cout << "\n========================================\n"
                << "  Results: " << num_passed << " passed, "
                << num_failed << " failed out of " << num_tests << " tests\n"
                << "========================================\n";
   }

   int exit_code = (num_failed == 0) ? 0 : 1;
   MPI_Bcast(&exit_code, 1, MPI_INT, 0, MPI_COMM_WORLD);
   MPI_Finalize();
   return exit_code;
}

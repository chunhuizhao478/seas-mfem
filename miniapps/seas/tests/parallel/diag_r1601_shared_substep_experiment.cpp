// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// diag_r1601_shared_substep_experiment.cpp — DIAGNOSTIC HARNESS (not a test;
// registered in no aggregate).  Root-cause experiment for the R-1601
// shared-fault substep divergence — see
// debug_document/tpv104_debug_document/R1601_root_cause_2026-07-10.md.
//
// Mirrors AdvanceADERWithSubStep_Spatial (drivers/spatial_dyn_driver.cpp:399)
// for N macro steps on the 2-tet shared-fault fixture and prints
// partition-invariant per-step metrics:
//
//   STEP k  maxQ=<global max|Q|>  l2=<global ||Q||_2>  slip2=...  V2=...
//
// Run matrix (compare stdout across runs).  NOTE: since Phase 2 landed the
// unified dispatch IS production — the former SEAS_DIAG_SHARED_SUBSTEP_CONSUME
// / _DUMP env gates were retired along with the R-1601 fallback, so np=2 with
// no env vars now exercises the unified path and must match np=1:
//   mpirun -np 1 ./seas_diag_r1601_experiment    -> ground truth
//   mpirun -np 2 ./seas_diag_r1601_experiment    -> unified (must match np=1)

#include "mfem.hpp"

#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/friction_substep_iterator.hpp"
#include "../../dynamic/friction_solver.hpp"
#include "../../friction/state_evolution.hpp"
#include "../../domain/boundary_config.hpp"
#include "shared_fault_fixtures.hpp"

#include <cmath>
#include <cstdio>
#include <string>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <vector>

using namespace mfem;
using namespace mfem::seas;
using mfem::seas::test_fixtures::BuildTwoTetSharedFault;
using mfem::seas::test_fixtures::PartitionByYSign;

namespace
{

constexpr real_t kRho = 2670.0;
constexpr real_t kCp  = 6000.0;
constexpr real_t kCs  = 3464.0;
const     real_t kZp  = kRho * kCp;
const     real_t kZs  = kRho * kCs;
const     real_t kMu  = kRho * kCs * kCs;
const     real_t kLam = kRho * kCp * kCp - 2.0 * kMu;

DOFData MakeFaultDof()
{
   DOFData d;
   d.Zp_plus  = kZp;  d.Zp_minus = kZp;
   d.Zs_plus  = kZs;  d.Zs_minus = kZs;
   d.eta_p    = 0.5 * kZp;
   d.eta_s    = 0.5 * kZs;
   d.sigma_n0 = 120.0e6;
   d.tau1_0   = 0.0;
   // Super-critical: 90 MPa > mu_s*sigma_n = 81.24 MPa, so the fault slips
   // from step 1 and any sign/routing error can feed back through Q.
   d.tau2_0   = 90.0e6;
   d.lsw_mu_s = 0.677;
   d.lsw_mu_d = 0.525;
   d.lsw_d_c  = 0.40;
   d.slip1    = 0.0;
   d.slip2    = 0.0;
   d.slip_rate_substep_max = 0.0;
   d.sigma_n_substep_min   = std::numeric_limits<real_t>::max();
   return d;
}

// Rate-and-state (aging, TPV102-style) DOFData — the HISTORICAL R-1601 law.
DOFData MakeAgingDof()
{
   DOFData d;
   d.Zp_plus  = kZp;  d.Zp_minus = kZp;
   d.Zs_plus  = kZs;  d.Zs_minus = kZs;
   d.eta_p    = 0.5 * kZp;
   d.eta_s    = 0.5 * kZs;
   d.sigma_n0 = 50.0e6;
   d.tau1_0   = 0.0;
   d.tau2_0   = 29.2e6;
   d.a   = 0.010;
   d.b   = 0.015;
   d.Dc  = 2.0;
   d.psi = 0.60;
   d.slip_rate_substep_max = 0.0;
   d.sigma_n_substep_min   = std::numeric_limits<real_t>::max();
   return d;
}

}  // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   const int rank   = Mpi::WorldRank();
   const int nprocs = Mpi::WorldSize();
   MFEM_VERIFY(nprocs == 1 || nprocs == 2,
               "experiment fixture is defined for np in {1,2}; got "
               << nprocs);

   // Fixture select: default = 2-tet (pure shared at np=2).
   // SEAS_DIAG_FIXTURE=hex4 -> 4-hex two-fault mesh with the MIXED partition
   // (rank 0 owns interior AND shared fault QPs), which distinguishes
   // absolute vs rebased indexing into the substep buffer (R-1304).
   const char *fix_env = std::getenv("SEAS_DIAG_FIXTURE");
   const std::string fix = fix_env ? fix_env : "2tet";
   ParMesh pmesh =
      (fix == "hex4")
         ? mfem::seas::test_fixtures::PartitionMixedInteriorShared(
              mfem::seas::test_fixtures::BuildFourHexTwoFault())
      : (fix == "tet2x2")
         ? mfem::seas::test_fixtures::PartitionOneInteriorOneShared(
              mfem::seas::test_fixtures::BuildTwoTetPairsTwoFaults())
         : PartitionByYSign(BuildTwoTetSharedFault());
   if (rank == 0) { std::printf("# fixture=%s\n", fix.c_str()); }

   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {};

   const int order = 1;
   const char *law_env = std::getenv("SEAS_DIAG_LAW");
   const bool use_rs = (law_env != nullptr && std::string(law_env) == "rs");
   if (rank == 0) { std::printf("# law=%s\n", use_rs ? "rate-state-aging" : "lsw"); }

   WaveOperator<ParMesh> wave(pmesh, order, kLam, kMu, kRho, bc);
   FaultFaceFlux flux(kRho, kCp, kCs);
   wave.SetFaultFlux(&flux);
   wave.SetFaultFrictionLaw(use_rs ? FaultFrictionLaw::RateAndState
                                   : FaultFrictionLaw::LSW);
   real_t Q_bg[NUM_STATE] = {0};
   wave.SetAbsorbingBackground(Q_bg);

   const int n_total = wave.GetNumTotalFaultQPs();
   const int nbf     = wave.GetNbfPerFace();
   if (rank == 0)
   {
      std::printf("# np=%d  interiorQPs=%d sharedQPs=%d total=%d\n",
                  nprocs, wave.GetNumLocalFaultQPs(),
                  wave.GetNumSharedFaultQPs(), n_total);
   }

   // Persistent fault state, exactly like the driver.
   std::vector<DOFData> dof_data(static_cast<size_t>(n_total));
   for (auto &d : dof_data) { d = use_rs ? MakeAgingDof() : MakeFaultDof(); }
   wave.SetFaultDOFData(&dof_data, nbf);

   std::vector<Vector> fault_coords(static_cast<size_t>(n_total), Vector(3));
   for (auto &v : fault_coords) { v = 0.0; }

   std::unique_ptr<IFrictionIterator> iterator_owner;
   if (use_rs)
   {
      iterator_owner = std::make_unique<RateStateAgingIterator>(
         flux, AgingLawPsi(/*b=*/0.015, /*V0=*/1.0e-6, /*f0=*/0.6),
         FrictionSolver::Method::Brent);
   }
   else
   {
      iterator_owner = std::make_unique<LinearSlipWeakeningIterator>(flux);
   }
   IFrictionIterator &iterator = *iterator_owner;
   const int    O  = 2;
   const real_t dt = 1.0e-6;   // ~5x below the 2-tet CFL limit (h_min~0.2m, cp=6km/s)
   {
      // Driver-style init; AdvanceADERWithSubStep rescales each step anyway.
      std::vector<real_t> deltaT(O, dt / O), weights(O, 1.0 / O);
      iterator.SetSubSteps(deltaT, weights);
   }

   const int ndof = wave.GetScalarNDof();
   Vector Q(NUM_STATE * ndof);
   Q = 0.0;
   for (int i = 0; i < ndof; ++i)
   {
      Q(SXY * ndof + i) = 1.0e6;
      Q(VX  * ndof + i) = 1.0e-3;
   }

   const auto nuc_noop = [](real_t, real_t) {};
   const int  NSTEPS   = 12;
   real_t t = 0.0;

   for (int step = 0; step < NSTEPS; ++step)
   {
      // ---- mirror of AdvanceADERWithSubStep_Spatial ---------------------
      const std::vector<real_t> cfg_dT = iterator.GetDeltaT();
      const std::vector<real_t> cfg_w  = iterator.GetTimeWeights();
      const real_t sum_dT =
         std::accumulate(cfg_dT.begin(), cfg_dT.end(), static_cast<real_t>(0));
      std::vector<real_t> dT(O);
      for (int o = 0; o < O; ++o) { dT[o] = cfg_dT[o] * (dt / sum_dT); }
      iterator.SetSubSteps(dT, cfg_w);

      std::vector<real_t> tau_nodes(O);
      real_t acc = 0.0;
      for (int o = 0; o < O; ++o)
      {
         tau_nodes[o] = acc + 0.5 * dT[o];
         acc += dT[o];
      }

      std::vector<Vector> Q_per_node;
      wave.ComputeADERSubStepStates(Q, dt, O, tau_nodes, Q_per_node);

      std::vector<std::vector<real_t>> Qp(O), Qm(O);
      for (int o = 0; o < O; ++o)
      { wave.EvaluateBulkAtFaultQPsCanonical(Q_per_node[o], Qp[o], Qm[o]); }

      const size_t n_words =
         static_cast<size_t>(NUM_STATE) * static_cast<size_t>(n_total);
      std::vector<real_t> Ip(n_words, 0.0), Im(n_words, 0.0);
      if (n_total > 0)
      {
         iterator.Advance(dof_data, fault_coords, Qp, Qm, dt, t,
                          Ip.data(), Im.data(), nuc_noop);
      }

      wave.SetSubStepFaultImposedStates(
         n_total > 0 ? Ip.data() : nullptr,
         n_total > 0 ? Im.data() : nullptr, n_total);
      Vector Q_new(Q.Size());
      wave.SetTime(t);
      wave.AdvanceADER(Q, dt, O, Q_new);
      wave.ResetSubStepFaultImposedStates();
      Q = Q_new;
      t += dt;
      // -------------------------------------------------------------------

      real_t loc_max = 0.0, loc_ss = 0.0;
      for (int i = 0; i < Q.Size(); ++i)
      {
         loc_max = std::max(loc_max, std::abs(Q(i)));
         loc_ss += Q(i) * Q(i);
      }
      real_t glob_max = 0.0, glob_ss = 0.0;
      MPI_Allreduce(&loc_max, &glob_max, 1, MPITypeMap<real_t>::mpi_type,
                    MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&loc_ss, &glob_ss, 1, MPITypeMap<real_t>::mpi_type,
                    MPI_SUM, MPI_COMM_WORLD);

      // Fault probe: QP 0 on the lowest rank that owns fault QPs.
      real_t probe[3] = {0, 0, 0};
      int have = (n_total > 0) ? rank : std::numeric_limits<int>::max();
      int owner = 0;
      MPI_Allreduce(&have, &owner, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      if (rank == owner && n_total > 0)
      {
         probe[0] = dof_data[0].slip2;
         probe[1] = dof_data[0].V2;
         probe[2] = dof_data[0].tau2_corr;
      }
      MPI_Bcast(probe, 3, MPITypeMap<real_t>::mpi_type, owner, MPI_COMM_WORLD);

      if (rank == 0)
      {
         std::printf("STEP %2d  maxQ=%.9e  l2=%.9e  slip2=%.9e  V2=%.9e  "
                     "tau2c=%.9e\n",
                     step, glob_max, std::sqrt(glob_ss),
                     probe[0], probe[1], probe[2]);
         std::fflush(stdout);
      }
      if (!(glob_max < 1e20))
      {
         if (rank == 0)
         { std::printf("BLOWUP at step %d (maxQ=%.3e) — R-1601 reproduced\n",
                       step, glob_max); }
         break;
      }
   }
   return 0;
}

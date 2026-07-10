// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_shared_fault_substep_parity_np2.cpp — Phase 0(b) of
// document/PLAN_unify_interior_shared_fault_substep_2026-07-09.md
// (review finding R-001).
//
// THE ORACLE for the unified substep dispatch (Phase 2 landed: the R-1601
// discard fallback is retired; both face classes consume the buffer).
//
// On the ADER + sub-step path the friction iterator computes a per-QP,
// time-integrated imposed state `I_imp` for EVERY fault QP (interior AND
// shared) and publishes it to the wave operator via
// `SetSubStepFaultImposedStates`.  The ADER corrector then:
//
//   * INTERIOR fault QPs  -> CONSUME that buffer (interior branch gate)
//   * SHARED   fault QPs  -> CONSUME it too (Phase 2 unified dispatch; the
//                            historical R-1601 discard fallback is retired —
//                            justification: debug_document/tpv104_debug_document/
//                            R1601_root_cause_2026-07-10.md)
//
// HOW THIS TEST DETECTS IT, without a friction solve or cross-partition QP
// matching: install a SENTINEL buffer and ask whether the corrector noticed.
//
//   Q_no_buf   = AdvanceADER(Q)                      with no buffer installed
//   Q_with_buf = AdvanceADER(Q)                      with a sentinel buffer
//
//   np=1 (fault face INTERIOR): Q_with_buf != Q_no_buf  -> buffer consumed.
//   np=2 (fault face SHARED)  : Q_with_buf != Q_no_buf  -> buffer consumed.
//
// Both rank counts must observe the sentinel (diff > 0): consumption is the
// unconditional production contract on interior AND shared fault QPs.
//
// DOFData is re-seeded before every AdvanceADER call.  Without that, the inline
// shared solve's WriteBackState (and the boss-broadcast reconcile) would mutate
// the fault state between the two runs, and Q_new would differ for a reason
// having nothing to do with the buffer.

#include "mfem.hpp"

#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../domain/boundary_config.hpp"
#include "shared_fault_fixtures.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

using namespace mfem;
using namespace mfem::seas;
using mfem::seas::test_fixtures::BuildTwoTetSharedFault;
using mfem::seas::test_fixtures::PartitionByYSign;

namespace
{

int g_rank = 0, g_nprocs = 1;
int g_num_tests = 0, g_num_passed = 0, g_num_failed = 0;

#define PARITY_ASSERT(cond, msg)                                              \
   do {                                                                       \
      ++g_num_tests;                                                          \
      if (cond) { ++g_num_passed;                                             \
         if (g_rank == 0) { std::cout << "  PASSED: " << msg << "\n"; } }     \
      else { ++g_num_failed;                                                  \
         std::cerr << "  FAILED [rank " << g_rank << ", line " << __LINE__    \
                   << "]: " << msg << "\n"; }                                 \
   } while (0)

constexpr real_t kRho = 2670.0;
constexpr real_t kCp  = 6000.0;
constexpr real_t kCs  = 3464.0;
const     real_t kZp  = kRho * kCp;
const     real_t kZs  = kRho * kCs;
const     real_t kMu  = kRho * kCs * kCs;
const     real_t kLam = kRho * kCp * kCp - 2.0 * kMu;

// LSW fault state.  `EvaluateADER_LSW` aborts unless at least one of
// lsw_mu_s / lsw_mu_d / lsw_d_c is non-zero, so seed all of them.
DOFData MakeFaultDof()
{
   DOFData d;
   d.Zp_plus  = kZp;  d.Zp_minus = kZp;
   d.Zs_plus  = kZs;  d.Zs_minus = kZs;
   d.eta_p    = 0.5 * kZp;
   d.eta_s    = 0.5 * kZs;
   d.sigma_n0 = 120.0e6;
   d.tau1_0   = 0.0;
   d.tau2_0   = 70.0e6;
   d.lsw_mu_s = 0.677;
   d.lsw_mu_d = 0.525;
   d.lsw_d_c  = 0.40;
   d.slip1    = 0.0;
   d.slip2    = 0.0;
   d.slip_rate_substep_max = 0.0;
   d.sigma_n_substep_min   = std::numeric_limits<real_t>::max();
   return d;
}

real_t MaxAbsDiff(const Vector &a, const Vector &b)
{
   MFEM_VERIFY(a.Size() == b.Size(), "MaxAbsDiff: size mismatch");
   real_t m = 0.0;
   for (int i = 0; i < a.Size(); ++i)
   { m = std::max(m, std::abs(a(i) - b(i))); }
   return m;
}

}  // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank   = Mpi::WorldRank();
   g_nprocs = Mpi::WorldSize();

   MFEM_VERIFY(g_nprocs == 1 || g_nprocs == 2,
               "fixture is defined for np=1 (interior fault face) or np=2 "
               "(shared fault face); got np=" << g_nprocs);

   if (g_rank == 0)
   {
      std::cout << "\n=== Phase 0(b): does the ADER corrector consume the "
                   "sub-step buffer? (np=" << g_nprocs << ") ===\n";
   }

   ParMesh pmesh = PartitionByYSign(BuildTwoTetSharedFault());

   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {};

   const int order = 1;
   WaveOperator<ParMesh> wave(pmesh, order, kLam, kMu, kRho, bc);
   FaultFaceFlux flux(kRho, kCp, kCs);
   wave.SetFaultFlux(&flux);
   wave.SetFaultFrictionLaw(FaultFrictionLaw::LSW);
   // AdvanceADER's corrector requires the absorbing-background pointer to have
   // been set on EVERY rank (R-1505 collective guard, wave_operator.inl:3963);
   // a rest state (Q_bg = 0) is correct for this fluctuation-Q fixture.
   real_t Q_bg[NUM_STATE] = {0};
   wave.SetAbsorbingBackground(Q_bg);

   const int n_int   = wave.GetNumLocalFaultQPs();
   const int n_shr   = wave.GetNumSharedFaultQPs();
   const int n_total = wave.GetNumTotalFaultQPs();
   const int nbf     = wave.GetNbfPerFace();

   if (g_rank == 0)
   {
      std::cout << "  fault QPs: interior=" << n_int << " shared=" << n_shr
                << " total=" << n_total << " (nbf/face=" << nbf << ")\n";
   }
   PARITY_ASSERT(n_total > 0, "fixture has fault QPs on this rank");
   if (g_nprocs == 1)
   { PARITY_ASSERT(n_int > 0 && n_shr == 0, "np=1: fault face is INTERIOR"); }
   else
   { PARITY_ASSERT(n_shr > 0 && n_int == 0, "np=2: fault face is SHARED"); }

#ifdef SEAS_TEST_INTERNAL
   // ------------------------------------------------------------------
   // NEGATIVE CONTROL for the Phase-3 reconcile assertion (review
   // REVIEW_phase2_3_unify_2026-07-10.md R-001).  With
   // SEAS_TEST_RECONCILE_NEGATIVE=1 at np=2, rank 1 tampers its
   // shared-fault routing bit, which swaps that rank's +/- inputs on the
   // INLINE (no-buffer) path and makes its DOFData/I_imp payload disagree
   // with rank 0's.  The reconcile assertion MUST abort (job exit != 0);
   // reaching the exit(1) below means the guard is dead code.
   //
   // NOTE: the tamper is only visible to the INLINE path — with the
   // sentinel buffer installed, both ranks' payloads are identical by
   // construction (seeded DOFData + the same sentinel), so the buffer leg
   // cannot serve as the negative control.
   // ------------------------------------------------------------------
   if (g_nprocs == 2 &&
       std::getenv("SEAS_TEST_RECONCILE_NEGATIVE") != nullptr)
   {
      if (g_rank == 1) { wave.TamperSharedFaultElem1OnPlus(0); }
      std::vector<DOFData> dof_neg(static_cast<size_t>(n_total));
      for (auto &d : dof_neg)
      {
         DOFData seed;
         seed.Zp_plus  = kZp;  seed.Zp_minus = kZp;
         seed.Zs_plus  = kZs;  seed.Zs_minus = kZs;
         seed.eta_p    = 0.5 * kZp;
         seed.eta_s    = 0.5 * kZs;
         seed.sigma_n0 = 120.0e6;
         seed.tau2_0   = 90.0e6;   // super-critical => nonzero solve outputs
         seed.lsw_mu_s = 0.677;
         seed.lsw_mu_d = 0.525;
         seed.lsw_d_c  = 0.40;
         d = seed;
      }
      wave.SetFaultDOFData(&dof_neg, nbf);
      wave.ResetSubStepFaultImposedStates();   // INLINE path (see NOTE above)
      Vector Q_neg(NUM_STATE * wave.GetScalarNDof());
      Q_neg = 0.0;
      // The field must be ASYMMETRIC across the fault: with a uniform Q,
      // Q_plus == Q_minus and the tampered +/- swap exchanges two identical
      // vectors — no observable mismatch.  Each rank owns one element of the
      // 2-tet fixture, so a rank-dependent amplitude makes the two fault
      // sides genuinely different.
      const real_t vx_side = 1.0e-3 * static_cast<real_t>(1 + g_rank);
      for (int i = 0; i < wave.GetScalarNDof(); ++i)
      {
         Q_neg(SXY * wave.GetScalarNDof() + i) = 1.0e6;
         Q_neg(VX  * wave.GetScalarNDof() + i) = vx_side;
      }
      Vector Q_neg_out(Q_neg.Size());
      wave.SetTime(0.0);
      wave.AdvanceADER(Q_neg, 1.0e-4, 2, Q_neg_out);   // must MFEM_ABORT
      if (g_rank == 0)
      {
         std::cerr << "NEGATIVE CONTROL FAILED: the reconcile assertion did "
                      "not fire on a manufactured cross-rank mismatch.\n";
      }
      return 1;
   }
#endif  // SEAS_TEST_INTERNAL

   // Non-trivial initial state: a uniform shear-stress + velocity field, so the
   // fault flux is not identically zero.
   const int ndof = wave.GetScalarNDof();
   Vector Q(NUM_STATE * ndof);
   Q = 0.0;
   for (int i = 0; i < ndof; ++i)
   {
      Q(SXY * ndof + i) = 1.0e6;
      Q(VX  * ndof + i) = 1.0e-3;
   }

   const real_t dt = 1.0e-4;
   const int    O  = 2;

   // A sentinel that the inline solve could never produce.
   std::vector<real_t> sentinel(static_cast<size_t>(NUM_STATE) * n_total, 1.0e-3);

   // Re-seed DOFData before EVERY AdvanceADER: the inline shared solve calls
   // WriteBackState and the boss-broadcast mutates it, so a stale copy would
   // make the two runs differ for the wrong reason.
   auto run = [&](bool install_buffer)
   {
      std::vector<DOFData> dof(static_cast<size_t>(n_total));
      for (auto &d : dof) { d = MakeFaultDof(); }
      wave.SetFaultDOFData(&dof, nbf);

      if (install_buffer && n_total > 0)
      {
         wave.SetSubStepFaultImposedStates(sentinel.data(), sentinel.data(),
                                           n_total);
      }
      else
      {
         wave.ResetSubStepFaultImposedStates();
      }

      Vector Q_new(Q.Size());
      Q_new = 0.0;
      wave.SetTime(0.0);
      wave.AdvanceADER(Q, dt, O, Q_new);
      wave.ResetSubStepFaultImposedStates();
      return Q_new;
   };

   const Vector Q_no_buf   = run(/*install_buffer=*/false);
   const Vector Q_with_buf = run(/*install_buffer=*/true);
   const real_t diff       = MaxAbsDiff(Q_no_buf, Q_with_buf);

   // Reduce so a rank that owns no fault QPs cannot mask a difference.
   real_t diff_global = 0.0;
   MPI_Allreduce(&diff, &diff_global, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX,
                 MPI_COMM_WORLD);
   if (g_rank == 0)
   {
      std::cout << "  max|Q_with_buf - Q_no_buf| = " << diff_global << "\n";
   }

   // Phase 2 of PLAN_unify_interior_shared_fault_substep (the R-1601 fallback
   // is retired; see debug_document/tpv104_debug_document/
   // R1601_root_cause_2026-07-10.md): interior AND shared fault QPs consume
   // the sub-step buffer identically.  The former RED oracle is now the
   // unconditional GREEN contract on both rank counts.
   if (g_nprocs == 1)
   {
      PARITY_ASSERT(diff_global > 0.0,
                    "np=1: INTERIOR fault QPs CONSUME the sub-step buffer "
                    "(installing a sentinel changes Q_new)");
   }
   else
   {
      PARITY_ASSERT(diff_global > 0.0,
                    "np=2: SHARED fault QPs CONSUME the sub-step buffer, "
                    "exactly as interior QPs do (unified substep dispatch; "
                    "supersedes the R-1601 discard fallback)");
   }

   int failed_global = 0, passed_global = 0, tests_global = 0;
   MPI_Allreduce(&g_num_failed, &failed_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&g_num_passed, &passed_global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&g_num_tests,  &tests_global,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if (g_rank == 0)
   {
      std::cout << "\n========================================\n";
      std::cout << "  Phase 0(b) substep-buffer parity: " << passed_global
                << " / " << tests_global << " passed, " << failed_global
                << " failed\n";
      std::cout << "========================================\n";
   }
   return (failed_global == 0) ? 0 : 1;
}

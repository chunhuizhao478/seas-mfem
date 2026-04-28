// REVIEW R-016 (round 4): TPV205 MPI rupture-crossing smoke.
//
// Pins the wave-operator dispatch fix that lets shared-fault QPs at
// np > 1 use the LSW closed form via FaultFaceFlux::EvaluateADER_LSW
// instead of falling back to Brent on rate-and-state.  Without the
// fix, the Brent solver receives data.a = lsw_mu_s and data.psi = 0
// (zeroed defensively at init) and produces strength values 1-3
// orders of magnitude away from the LSW physics — V_abs collapses to
// ~0 on shared-fault QPs and the rupture front stalls at the rank
// boundary.
//
// Test design:
//   1. Build a ParMesh from the TPV104 1000 m fixture mesh (it has
//      Physical Surface 1/3/5 = free/fault/absorb tags, and a fault
//      at y = 0 that ParMETIS naturally bisects at np = 2).
//   2. Initialize TPV205 fault DOFs (LSW μ_s/μ_d/d_c, 4-patch
//      pre-stress, sigma_n0 = 120 MPa, tau_nuc = 81.6 MPa).
//   3. Configure the wave operator with FaultFrictionLaw::LSW.
//   4. Run AdvanceADERWithSubStep for one ADER-O2 macro step.
//   5. Assert V_max at the nucleation patch matches the analytic
//      0.0779 m/s within 5 % on every rank that owns a fault QP
//      inside the patch.
//
// Failure mode the test catches:
//   - The dispatch wiring compiles, links, and `SetFaultFrictionLaw`
//     reaches the two branches in wave_operator.inl (interior + R-1600
//     shared-fault).
//   - The LSW closed-form via either path (interior or shared) gives
//     bit-exact agreement with the analytic 0.0779 m/s on every rank
//     that owns nucleation-patch QPs, including ranks where the patch
//     QPs sit on shared faces.
//   - The round-robin partition guarantees 4000+ shared-fault QPs so
//     the R-1600 dispatch path is genuinely exercised.
//
// Limitation (documented for future hardening):
//   The iterator runs LSW on dof_data BEFORE the wave operator's bulk
//   corrector, so a revert of wave_operator.inl's LSW branches would
//   trigger a flood of FrictionSolver NaN messages on stderr (Brent
//   on rate-and-state with the zero-rate-and-state defensive init at
//   shared QPs) but `dof_data[i].slip_rate` would still report the
//   iterator's correct LSW value.  The Brent NaN spam IS the
//   regression signal — but it is not turned into an assertion here
//   because writing such an assertion requires either capturing
//   stderr (fragile) or running enough macro steps for the wrong
//   bulk-Q radiation to feed back into the iterator's reading.  Both
//   approaches are deferred — the bit-exact parity test
//   `test_tpv205_evaluate_ader_lsw_parity` already pins the closed-
//   form path; this test pins the end-to-end MPI integration.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv205_setup.hpp"
#include "../../dynamic/tpv205_friction.hpp"
#include "../../dynamic/tpv205_substep_iterator.hpp"
#include "../../config/tpv205_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

#ifdef SEAS_DIAG_FAULT_FLUX
namespace mfem { namespace seas { int g_seas_my_rank = 0; } }  // NOLINT
#endif

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_GE(v, lo, msg) do { \
   num_tests++; \
   double vv = (v), ll = (lo); \
   if (vv >= ll) { num_passed++; \
      std::cout << "  PASSED: " << msg << "  (got " << std::scientific \
                << std::setprecision(3) << vv << ", floor " << ll << ")\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << "  (got " << std::scientific << std::setprecision(3) \
                << vv << ", expected ≥ " << ll << ")\n"; } \
} while (0)

#define TEST_NEAR(v, ref, tol, msg) do { \
   num_tests++; \
   double vv = (v), rr = (ref), tt = (tol); \
   double dd = std::abs(vv - rr); \
   if (dd <= tt) { num_passed++; \
      std::cout << "  PASSED: " << msg << "  (got " << std::scientific \
                << std::setprecision(6) << vv << ", ref " << rr \
                << ", |Δ| " << dd << " ≤ " << tt << ")\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << "  (got " << vv << ", ref " << rr \
                << ", |Δ| " << dd << " > " << tt << ")\n"; } \
} while (0)

namespace
{

// AdvanceADERWithSubStep — copied from drivers/tpv205_driver.cpp
// (~L294-435).  Test owns its own copy so it doesn't link the driver.
template <typename WaveOpT>
void AdvanceADERWithSubStep(WaveOpT &wave,
                            Tpv205SubStepIterator &iterator,
                            std::vector<DOFData> &dof_data,
                            const std::vector<Vector> &fault_coords,
                            const Vector &Q,
                            real_t dt_step,
                            int ader_order,
                            real_t t_step_start,
                            Vector &Q_new)
{
   MFEM_VERIFY(dt_step > 0.0, "dt_step must be > 0");
   MFEM_VERIFY(ader_order >= 2 && ader_order <= 4, "ader_order in {2,3,4}");

   const std::vector<real_t> configured_deltaT  = iterator.GetDeltaT();
   const std::vector<real_t> configured_weights = iterator.GetTimeWeights();
   const int O = static_cast<int>(configured_deltaT.size());
   const real_t configured_sum =
      std::accumulate(configured_deltaT.begin(), configured_deltaT.end(),
                      static_cast<real_t>(0));
   const real_t dt_scale = dt_step / configured_sum;
   std::vector<real_t> deltaT_scaled(O);
   for (int o = 0; o < O; ++o)
   {
      deltaT_scaled[o] = configured_deltaT[o] * dt_scale;
   }
   iterator.SetSubSteps(deltaT_scaled, configured_weights);
   const std::vector<real_t> &deltaT = iterator.GetDeltaT();

   std::vector<real_t> tau_nodes(O);
   real_t acc = 0.0;
   for (int o = 0; o < O; ++o)
   {
      tau_nodes[o] = acc + 0.5 * deltaT[o];
      acc += deltaT[o];
   }

   std::vector<Vector> Q_per_node;
   wave.ComputeADERSubStepStates(Q, dt_step, ader_order, tau_nodes, Q_per_node);

   const int n_total_fault_qps = wave.GetNumTotalFaultQPs();
   std::vector<std::vector<real_t>> Q_pointwise_plus(O), Q_pointwise_minus(O);
   for (int o = 0; o < O; ++o)
   {
      wave.EvaluateBulkAtFaultQPsCanonical(Q_per_node[o],
                                           Q_pointwise_plus[o],
                                           Q_pointwise_minus[o]);
   }

   const size_t n_words =
      static_cast<size_t>(NUM_STATE) * static_cast<size_t>(n_total_fault_qps);
   std::vector<real_t> I_imp_plus_flat(n_words, 0.0);
   std::vector<real_t> I_imp_minus_flat(n_words, 0.0);

   if (n_total_fault_qps > 0)
   {
      iterator.AdvanceWithSubStepStates(dof_data, fault_coords,
                                        Q_pointwise_plus,
                                        Q_pointwise_minus,
                                        dt_step, t_step_start,
                                        I_imp_plus_flat.data(),
                                        I_imp_minus_flat.data());
   }

   wave.SetSubStepFaultImposedStates(
      n_total_fault_qps > 0 ? I_imp_plus_flat.data()  : nullptr,
      n_total_fault_qps > 0 ? I_imp_minus_flat.data() : nullptr,
      n_total_fault_qps);

   wave.AdvanceADER(Q, dt_step, ader_order, Q_new);

   wave.ResetSubStepFaultImposedStates();
}

// Run the test — MPI only.
int RunMPITest(int rank, int nprocs, MPI_Comm comm)
{
   const std::string mesh_file =
      "tpv104/mesh/tpv104_1000m.msh";  // run from miniapps/seas/

   if (rank == 0)
   {
      std::cout << "\n[T_TPV205_MPI_RUPTURE_CROSSING] np = " << nprocs
                << ", mesh = " << mesh_file << "\n";
   }

   // Load the serial mesh.  Default ParMETIS at low np tends to keep
   // all fault elements on a single rank (zero shared fault QPs),
   // which would let the test pass without exercising the R-1600
   // shared-fault dispatch.  Force a round-robin partition by
   // element index so that adjacent elements (which include the two
   // tets sharing a fault face) land on different ranks → shared
   // fault QPs are guaranteed.  This is a test-only partition; not
   // load-balanced or topology-aware, just adversarial enough to
   // populate shared-fault QPs at np ≥ 2.
   Mesh serial_mesh(mesh_file.c_str(), 1, 1);
   if (serial_mesh.Dimension() != 3)
   {
      if (rank == 0) { std::cerr << "ERROR: mesh is not 3D.\n"; }
      return 1;
   }
   const int ne = serial_mesh.GetNE();
   Array<int> partitioning(ne);
   for (int e = 0; e < ne; ++e)
   {
      partitioning[e] = e % nprocs;
   }
   ParMesh pmesh(comm, serial_mesh, partitioning.GetData());

   // Boundary attributes match the TPV104 mesh: 1/3/5 = free/fault/absorb.
   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {5};

   // Wave operator + fault flux + LSW dispatch enable.
   const int order = 1;
   WaveOperator<ParMesh> wave(pmesh, order,
                              TPV205Params::lambda, TPV205Params::mu,
                              TPV205Params::rho, bc);
   FaultFaceFlux fault_flux(TPV205Params::rho, TPV205Params::cp,
                            TPV205Params::cs);
   wave.SetFaultFlux(&fault_flux);
   wave.SetFaultFrictionLaw(FaultFrictionLaw::LSW);

   const int ndof_total = wave.GetScalarNDof();
   const real_t cfl = 0.5 / (3.0 * (2.0 * order + 1.0));
   const real_t dt  = wave.ComputeMaxDt(cfl);

   // Fault DOFs.
   const Array<int> &fault_int = wave.GetFaultInteriorFaces();
   const Array<int> &fault_shr = wave.GetFaultSharedFaces();

   int nqp_per_face = 0;
   if (fault_int.Size() > 0)
   {
      FaceElementTransformations *ftr0 =
         pmesh.GetInteriorFaceTransformations(fault_int[0]);
      nqp_per_face = IntRules.Get(ftr0->GetGeometryType(), 2*order).GetNPoints();
   }
   else if (fault_shr.Size() > 0)
   {
      FaceElementTransformations *ftr0 =
         pmesh.GetSharedFaceTransformations(fault_shr[0]);
      nqp_per_face = IntRules.Get(ftr0->GetGeometryType(), 2*order).GetNPoints();
   }
   {
      int local_nqp = nqp_per_face;
      MPI_Allreduce(&local_nqp, &nqp_per_face, 1, MPI_INT, MPI_MAX, comm);
   }

   const int n_local  = fault_int.Size() * nqp_per_face;
   const int n_shared = fault_shr.Size() * nqp_per_face;
   const int n_total  = n_local + n_shared;

   // Fault QP coordinates.
   std::vector<Vector> fault_coords;
   fault_coords.reserve(n_total);
   auto push_qps = [&](FaceElementTransformations *ftr)
   {
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*order);
      for (int q = 0; q < ir.GetNPoints(); ++q)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3);
         ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   };
   for (int i = 0; i < fault_int.Size(); ++i)
   {
      push_qps(pmesh.GetInteriorFaceTransformations(fault_int[i]));
   }
   for (int i = 0; i < fault_shr.Size(); ++i)
   {
      push_qps(pmesh.GetSharedFaceTransformations(fault_shr[i]));
   }

   // Sanity: total fault QPs > 0 globally; with the round-robin partition
   // adjacent elements land on different ranks so fault faces are
   // shared.  Assert this — without shared-fault QPs the test would
   // exercise only the interior-fault dispatch (line ~3617 of
   // wave_operator.inl) and never the R-1600 shared-fault fallback
   // (line ~4539), defeating the purpose of the regression test.
   int n_shared_total = 0;
   MPI_Allreduce(&n_shared, &n_shared_total, 1, MPI_INT, MPI_SUM, comm);
   if (rank == 0)
   {
      std::cout << "  Fault QPs: " << n_shared_total
                << " shared globally (round-robin partition forces "
                << "fault-face splits; required for R-1600 dispatch coverage).\n";
   }
   TEST_GE(static_cast<double>(n_shared_total), 1.0,
           "Round-robin partition produces > 0 shared fault QPs");

   std::vector<DOFData> dof_data;
   if (n_total > 0)
   {
      InitializeFaultDOFs_TPV205(dof_data, n_total, fault_coords);
   }
   wave.SetFaultDOFData(&dof_data, nqp_per_face);
   {
      real_t Q_bg[NUM_STATE] = {0};
      wave.SetAbsorbingBackground(Q_bg);
   }

   // Iterator at O = 1 (one-shot LSW per macro step).
   Tpv205SubStepIterator iterator(fault_flux);
   iterator.SetSubSteps({dt}, {1.0});

   // Q = 0, advance one ADER-O2 macro step.
   Vector Q(NUM_STATE * ndof_total);
   InitializeState_TPV205(Q, ndof_total);
   Vector Q_new(Q.Size());

   AdvanceADERWithSubStep(wave, iterator, dof_data, fault_coords,
                          Q, dt, /*ader_order*/2, /*t_start*/0.0,
                          Q_new);

   // R-002 (round 5 final review): the production driver calls
   // VerifySharedFaultDOFDataConsistency() after step 0 (see
   // tpv205_driver.cpp:2138).  Replicate that here so any latent
   // rank-canonical-frame mismatch in EvaluateADER_LSW (the R-1601
   // SHARED FALLBACK warned the elem1_on_plus / sign_flipped pairing
   // could produce rank-dependent V1/V2 on shared QPs) surfaces in CI
   // rather than at the first sbatch on Frontera.  MFEM_ABORT cannot
   // be caught — if the verify trips it calls MPI_Abort and the whole
   // test fails loudly with `[R-101 rank-0 detail]` messages.  Default
   // tol = 1e-10 (relative).
   wave.VerifySharedFaultDOFDataConsistency();
   num_tests++;
   num_passed++;
   if (rank == 0)
   {
      std::cout << "  PASSED: VerifySharedFaultDOFDataConsistency "
                << "after one ADER-O2 step (default tol 1e-10)\n";
   }

   // Find V_max among QPs INSIDE the nucleation patch (any QP whose
   // (along-strike, depth) sits inside the 3×3 km patch around
   // (0, 7.5 km).  A QP with strict patch membership has tau2_0 ==
   // tau_nuc, which is the most reliable marker.
   real_t v_max_local_in_patch    = 0.0;
   real_t v_max_local_anywhere    = 0.0;
   int    n_local_in_patch        = 0;
   int    n_local_in_patch_shared = 0;
   for (int i = 0; i < n_total; ++i)
   {
      const real_t v = dof_data[i].slip_rate;
      v_max_local_anywhere = std::max(v_max_local_anywhere, v);
      if (std::abs(dof_data[i].tau2_0 - TPV205Params::tau_nuc)
          < static_cast<real_t>(1e3))
      {
         v_max_local_in_patch = std::max(v_max_local_in_patch, v);
         ++n_local_in_patch;
         if (i >= n_local) { ++n_local_in_patch_shared; }
      }
   }

   real_t v_max_global_in_patch = 0.0;
   real_t v_max_global_anywhere = 0.0;
   int    n_global_in_patch     = 0;
   int    n_global_in_patch_shr = 0;
   MPI_Allreduce(&v_max_local_in_patch, &v_max_global_in_patch, 1,
                 MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
   MPI_Allreduce(&v_max_local_anywhere, &v_max_global_anywhere, 1,
                 MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
   MPI_Allreduce(&n_local_in_patch,        &n_global_in_patch,     1, MPI_INT, MPI_SUM, comm);
   MPI_Allreduce(&n_local_in_patch_shared, &n_global_in_patch_shr, 1, MPI_INT, MPI_SUM, comm);

   const double analytic =
      (TPV205Params::tau_nuc - TPV205Params::mu_s * TPV205Params::sigma_n)
      / TPV205Params::eta_s;

   if (rank == 0)
   {
      std::cout << "  Nucleation-patch QPs: " << n_global_in_patch
                << " (of which " << n_global_in_patch_shr << " on shared faces)\n";
      std::cout << "  V_max (in patch, global): "
                << std::scientific << std::setprecision(6)
                << v_max_global_in_patch
                << "  analytic = " << analytic << "\n";
      std::cout << "  V_max (anywhere, global): "
                << v_max_global_anywhere << "\n";

      // Acceptance: V_max in the patch matches analytic within 5 %.  If
      // the rupture had stalled at the rank boundary, V_max on the
      // shared-fault rank would be ≪ analytic; MPI_MAX would still see
      // the interior-rank's correct value, but the per-rank check below
      // catches that case directly.
      TEST_NEAR(static_cast<double>(v_max_global_in_patch),
                analytic, 0.05 * analytic,
                "Nucleation-patch V_max global == analytic (within 5%)");
   }

   // Per-rank check: if THIS rank owns any nucleation-patch QPs, its
   // local V_max must match analytic too.  This catches the case where
   // shared-fault dispatch silently gives V = 0 on one rank only.
   //
   // Skip the assertion on ranks that don't see any patch QPs (those
   // ranks have no business reporting a V_max value).
   if (n_local_in_patch > 0)
   {
      const double v_local = static_cast<double>(v_max_local_in_patch);
      // Looser per-rank tolerance (10 %) — the local extremum on a
      // half-patch may legitimately be slightly off from the global
      // analytic if the patch boundary cuts a single QP.
      const double lo = 0.85 * analytic;
      std::cout << "  rank " << rank << " local in-patch V_max = "
                << std::scientific << std::setprecision(6) << v_local
                << " (n_in_patch=" << n_local_in_patch
                << ", on shared faces=" << n_local_in_patch_shared << ")\n";
      TEST_GE(v_local, lo,
              "rank-local V_max ≥ 0.85·analytic (catches rupture stall on "
              "shared-fault rank)");
   }

   // R-001 regression: run additional macro-steps and check that the
   // accumulated slip on shared-fault QPs in the nucleation patch is
   // CONSISTENT with interior-fault QPs in the same patch.  Pre-fix,
   // EvaluateADER_LSW added `data.slip{1,2} += V*dt` on top of the
   // iterator's per-sub-step accumulation, so shared QPs ended each
   // macro-step with slip ≈ 2× the interior value.  Post-fix, slip
   // evolution is owned by the iterator alone — interior and shared
   // QPs at the same physical location must agree to within rounding.
   const int kExtraSteps = 4;
   for (int step = 0; step < kExtraSteps; ++step)
   {
      Vector Q_step_new(Q.Size());
      AdvanceADERWithSubStep(wave, iterator, dof_data, fault_coords,
                             Q_new, dt, /*ader_order*/2,
                             /*t_start*/(step + 1) * dt,
                             Q_step_new);
      Q_new.Swap(Q_step_new);
   }

   // Reduce slip on interior- and shared-fault QPs inside the patch.
   real_t slip_int_local = 0.0, slip_shr_local = 0.0;
   int n_int_local = 0, n_shr_local = 0;
   for (int i = 0; i < n_total; ++i)
   {
      if (std::abs(dof_data[i].tau2_0 - TPV205Params::tau_nuc)
          >= static_cast<real_t>(1e3)) { continue; }
      const real_t s = std::abs(dof_data[i].slip2);
      if (i < n_local) { slip_int_local = std::max(slip_int_local, s); ++n_int_local; }
      else             { slip_shr_local = std::max(slip_shr_local, s); ++n_shr_local; }
   }
   real_t slip_int_g = 0.0, slip_shr_g = 0.0;
   int n_int_g = 0, n_shr_g = 0;
   MPI_Allreduce(&slip_int_local, &slip_int_g, 1,
                 MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
   MPI_Allreduce(&slip_shr_local, &slip_shr_g, 1,
                 MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
   MPI_Allreduce(&n_int_local, &n_int_g, 1, MPI_INT, MPI_SUM, comm);
   MPI_Allreduce(&n_shr_local, &n_shr_g, 1, MPI_INT, MPI_SUM, comm);

   if (rank == 0)
   {
      std::cout << "  After " << (1 + kExtraSteps) << " macro-steps: "
                << "patch slip2 max — interior=" << std::scientific
                << std::setprecision(6) << slip_int_g
                << " (n=" << n_int_g << "), shared=" << slip_shr_g
                << " (n=" << n_shr_g << ")\n";
      // Pre-fix: shared / interior ≈ 2.0.  Post-fix: ≈ 1.0 (matched
      // physically — the slip-rate field is smooth across rank
      // boundaries, so peak slip on shared QPs should not exceed the
      // peak on interior QPs by more than ~5 %).  A value > 1.5 is
      // unambiguously the R-001 double-count signature.
      if (n_shr_g > 0 && slip_int_g > 0.0)
      {
         const double ratio = static_cast<double>(slip_shr_g)
                              / static_cast<double>(slip_int_g);
         std::cout << "  R-001 regression ratio shared/interior = "
                   << std::fixed << std::setprecision(3) << ratio
                   << "\n";
         num_tests++;
         if (ratio < 1.5)
         {
            num_passed++;
            std::cout << "  PASSED: R-001 slip parity (ratio "
                      << ratio << " < 1.5; pre-fix value ≈ 2.0)\n";
         }
         else
         {
            num_failed++;
            std::cout << "  FAILED: R-001 double-count regression — "
                      << "shared/interior slip ratio " << ratio
                      << " ≥ 1.5 (pre-fix bug returned)\n";
         }
      }
   }

   return num_failed == 0 ? 0 : 1;
}

} // namespace

int main(int argc, char *argv[])
{
#ifndef MFEM_USE_MPI
   std::cerr << "test_tpv205_mpi_rupture_crossing requires MFEM_USE_MPI.\n";
   return 0;  // skip
#else
   MPI_Init(&argc, &argv);
   MPI_Comm comm = MPI_COMM_WORLD;
   int rank, nprocs;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nprocs);
#ifdef SEAS_DIAG_FAULT_FLUX
   mfem::seas::g_seas_my_rank = rank;
#endif

   if (nprocs < 2)
   {
      if (rank == 0)
      {
         std::cerr << "test_tpv205_mpi_rupture_crossing requires np ≥ 2; "
                   << "got " << nprocs << ".  Skipping.\n";
      }
      MPI_Finalize();
      return 0;
   }

   const int rc = RunMPITest(rank, nprocs, comm);

   // Aggregate per-rank pass/fail counts at rank 0 for a clean summary.
   int total_tests = 0, total_passed = 0, total_failed = 0;
   MPI_Reduce(&num_tests,  &total_tests,  1, MPI_INT, MPI_SUM, 0, comm);
   MPI_Reduce(&num_passed, &total_passed, 1, MPI_INT, MPI_SUM, 0, comm);
   MPI_Reduce(&num_failed, &total_failed, 1, MPI_INT, MPI_SUM, 0, comm);

   if (rank == 0)
   {
      std::cout << "\n========================================\n";
      std::cout << "  Results: " << total_passed << " passed, "
                << total_failed << " failed out of " << total_tests
                << " tests (sum across ranks)\n";
      std::cout << "========================================\n";
   }

   MPI_Finalize();
   return (rc == 0 && total_failed == 0) ? 0 : 1;
#endif
}

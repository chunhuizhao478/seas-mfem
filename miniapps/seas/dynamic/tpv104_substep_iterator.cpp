// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// Implementation of TPV104 sub-step fault iterator (§4.10 Step 7).

#include "tpv104_substep_iterator.hpp"
// friction_coeff_stable.hpp reached transitively via tpv104_friction_solver.hpp.
#include "../friction/friction_coeff_stable.hpp"

#include <cmath>
#include <cstring>
#include <numeric>
#include <stdexcept>
#include <string>

#ifdef SEAS_DIAG_TPV104_STATE
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <ios>
#include <memory>
#include <sstream>
#include <unordered_map>
#endif  // SEAS_DIAG_TPV104_STATE

namespace mfem
{
namespace seas
{

#ifdef SEAS_DIAG_TPV104_STATE
// Probe writer — one file per probe per rank (TPV104 Phase-3 §5.1 format).
// Files are opened lazily on first emit and closed either by an
// explicit `CloseAllProbeFiles` (R5-006 — call before MPI_Finalize) or
// at process exit via static-destructor order.
namespace
{
// Accessor returns a reference to the static registry so
// `CloseAllProbeFiles` can clear it without breaking the anonymous-
// namespace isolation.  Magic-statics make this thread-safe under C++11+.
using ProbeFileMap = std::unordered_map<std::string,
                                        std::unique_ptr<std::ofstream>>;
ProbeFileMap &GetProbeFileRegistry()
{
   static ProbeFileMap files;
   return files;
}

std::ofstream &GetProbeFile(const char *probe_name)
{
   ProbeFileMap &files = GetProbeFileRegistry();
   auto it = files.find(probe_name);
   if (it == files.end())
   {
      const char *dir_env = std::getenv("SEAS_DIAG_TPV104_DIR");
      std::string dir = dir_env ? dir_env : ".";

      // R4-005 (review round 4): derive the per-process rank suffix in
      // the following priority order:
      //   (1) SEAS_DIAG_TPV104_RANK env var (explicit caller override);
      //   (2) MPI_Comm_rank(MPI_COMM_WORLD) under MFEM_USE_MPI if MPI
      //       has been initialised — avoids all four ranks clobbering
      //       the same file under mpirun -np 4 without per-process env
      //       injection;
      //   (3) "0" fallback for serial / pre-MPI_Init paths.
      // Same priority for the `nprocs` banner field.
      std::string rank = "0";
      std::string nprocs = "1";
      if (const char *rank_env = std::getenv("SEAS_DIAG_TPV104_RANK"))
      {
         rank = rank_env;
      }
#ifdef MFEM_USE_MPI
      else
      {
         int mpi_inited = 0;
         MPI_Initialized(&mpi_inited);
         if (mpi_inited)
         {
            int r = 0, n = 1;
            MPI_Comm_rank(MPI_COMM_WORLD, &r);
            MPI_Comm_size(MPI_COMM_WORLD, &n);
            rank   = std::to_string(r);
            nprocs = std::to_string(n);
         }
      }
#endif
      if (const char *np_env = std::getenv("SEAS_DIAG_TPV104_NPROCS"))
      {
         nprocs = np_env;
      }

      std::string path = dir + "/tpv104_probe_" + std::string(probe_name)
                         + "_rank" + rank + ".txt";
      auto f = std::make_unique<std::ofstream>(path);
      if (f->is_open())
      {
         (*f) << "# probe=" << probe_name
              << "  code=MFEM  rank=" << rank
              << "  nprocs=" << nprocs << "\n";
         // R4-004 disclosure: flag channels whose MFEM emission
         // differs from the reference runtime in per-sub-step cadence.
         if (probe_name == std::string("trial_traction"))
         {
            (*f) << "# cadence_note: MFEM emits O identical rows per "
                    "macro-step (Q̄ = I/dt_macro used for every "
                    "sub-step); reference emits O distinct rows. "
                    "Step-13 probe_diff must coarsen by sub-step "
                    "averaging or mask before comparing.\n";
         }
      }
      it = files.emplace(probe_name, std::move(f)).first;
   }
   return *(it->second);
}
}  // anonymous namespace
#endif  // SEAS_DIAG_TPV104_STATE

Tpv104SubStepIterator::Tpv104SubStepIterator(FaultFaceFlux &flux,
                                             const SlipLawSRWPsi &state_evo)
   : flux_(flux), state_evo_(state_evo)
{
}

void Tpv104SubStepIterator::CloseAllProbeFiles()
{
#ifdef SEAS_DIAG_TPV104_STATE
   ProbeFileMap &files = GetProbeFileRegistry();
   for (auto &kv : files)
   {
      if (kv.second && kv.second->is_open())
      {
         kv.second->flush();
         kv.second->close();
      }
   }
   files.clear();
#endif  // SEAS_DIAG_TPV104_STATE
}

void Tpv104SubStepIterator::SetSubSteps(std::vector<real_t> deltaT,
                                        std::vector<real_t> time_weights)
{
   if (deltaT.empty())
   {
      throw std::runtime_error(
         "Tpv104SubStepIterator::SetSubSteps: deltaT must have at least "
         "one entry.");
   }
   if (deltaT.size() != time_weights.size())
   {
      throw std::runtime_error(
         "Tpv104SubStepIterator::SetSubSteps: deltaT and time_weights "
         "must have equal size; got deltaT.size() = "
         + std::to_string(deltaT.size())
         + ", time_weights.size() = "
         + std::to_string(time_weights.size()));
   }
   for (size_t o = 0; o < deltaT.size(); ++o)
   {
      if (!std::isfinite(deltaT[o]) || deltaT[o] <= 0.0)
      {
         throw std::runtime_error(
            "Tpv104SubStepIterator::SetSubSteps: deltaT[" + std::to_string(o)
            + "] must be finite and positive; got " + std::to_string(deltaT[o]));
      }
      if (!std::isfinite(time_weights[o]))
      {
         throw std::runtime_error(
            "Tpv104SubStepIterator::SetSubSteps: time_weights["
            + std::to_string(o) + "] must be finite; got "
            + std::to_string(time_weights[o]));
      }
   }
   const real_t wsum = std::accumulate(time_weights.begin(),
                                       time_weights.end(),
                                       static_cast<real_t>(0));
   if (std::abs(wsum - 1.0) > 1e-12)
   {
      throw std::runtime_error(
         "Tpv104SubStepIterator::SetSubSteps: time_weights must sum to "
         "1 within 1e-12; got sum = " + std::to_string(wsum));
   }

   deltaT_       = std::move(deltaT);
   time_weights_ = std::move(time_weights);
}

// R4-002: default-argument value lives in the header (Brent).  The
// definition omits the default per standard practice but must match
// the header's declared default.
void Tpv104SubStepIterator::Advance(std::vector<DOFData> &dof_data,
                                    const std::vector<Vector> &fault_coords,
                                    const std::vector<real_t> &V_w,
                                    const real_t *I_plus_flat,
                                    const real_t *I_minus_flat,
                                    real_t dt_macro,
                                    real_t t_macro_start,
                                    real_t *I_imp_plus_flat,
                                    real_t *I_imp_minus_flat,
                                    FrictionSolver::Method method)
{
   if (deltaT_.empty())
   {
      throw std::runtime_error(
         "Tpv104SubStepIterator::Advance: SetSubSteps has not been "
         "called; cannot iterate without a configured quadrature.");
   }
   if (!std::isfinite(dt_macro) || dt_macro <= 0.0)
   {
      throw std::runtime_error(
         "Tpv104SubStepIterator::Advance: dt_macro must be finite and "
         "positive; got " + std::to_string(dt_macro));
   }
   if (!std::isfinite(t_macro_start))
   {
      throw std::runtime_error(
         "Tpv104SubStepIterator::Advance: t_macro_start must be finite; "
         "got " + std::to_string(t_macro_start));
   }
   if (I_plus_flat == nullptr || I_minus_flat == nullptr
       || I_imp_plus_flat == nullptr || I_imp_minus_flat == nullptr)
   {
      throw std::runtime_error(
         "Tpv104SubStepIterator::Advance: all four I_* pointers must be "
         "non-null.");
   }
   const int n = static_cast<int>(dof_data.size());
   if (static_cast<int>(fault_coords.size()) != n
       || static_cast<int>(V_w.size()) != n)
   {
      throw std::runtime_error(
         "Tpv104SubStepIterator::Advance: dof_data, fault_coords, and "
         "V_w must have equal size; got "
         + std::to_string(dof_data.size()) + ", "
         + std::to_string(fault_coords.size()) + ", "
         + std::to_string(V_w.size()));
   }

   // Verify Σ deltaT == dt_macro.  Accept small accumulated drift.
   // R5-007 (review round 5): scale the tolerance by O so the check
   // is not artificially tight under very high-order quadratures.
   // Floor at 1e-12 keeps the gate useful at low O.  `eps ≈ 2.2e-16`
   // per add, so `O * eps` with a constant safety margin is the
   // natural drift bound.
   const real_t dtsum = std::accumulate(deltaT_.begin(), deltaT_.end(),
                                        static_cast<real_t>(0));
   const real_t rel = std::abs(dtsum - dt_macro) / std::max(dt_macro,
                                                            1e-300);
   const int    O_size = static_cast<int>(deltaT_.size());
   const real_t sum_tol = std::max<real_t>(1e-12,
                                           static_cast<real_t>(10.0) * O_size
                                             * std::numeric_limits<real_t>::epsilon());
   if (rel > sum_tol)
   {
      throw std::runtime_error(
         "Tpv104SubStepIterator::Advance: Σ deltaT[o] must equal "
         "dt_macro within " + std::to_string(sum_tol)
         + " rel; got Σ = " + std::to_string(dtsum)
         + ", dt_macro = " + std::to_string(dt_macro)
         + ", O = " + std::to_string(O_size));
   }

   // Zero the output accumulators.
   const size_t nwords = static_cast<size_t>(NUM_STATE) * n;
   std::memset(I_imp_plus_flat,  0, nwords * sizeof(real_t));
   std::memset(I_imp_minus_flat, 0, nwords * sizeof(real_t));

   const real_t inv_dt_macro = 1.0 / dt_macro;

   // Cache per-QP time-averaged bulk state once — Q_avg is constant
   // across sub-steps under our interpretation of the predictor (§4.10
   // Step 7 ambiguity resolution).  Stack-allocated per-DOF working
   // buffers are tiny (9 real_t each); a single heap buffer of size
   // NUM_STATE keeps the per-sub-step loop allocation-free.
   real_t Q_avg_plus[NUM_STATE];
   real_t Q_avg_minus[NUM_STATE];
   real_t Q_imp_plus[NUM_STATE];
   real_t Q_imp_minus[NUM_STATE];

   // Track elapsed sub-step time so nucleation endpoints are absolute.
   real_t t_sub_cursor = t_macro_start;

   const int O = static_cast<int>(deltaT_.size());
   for (int o = 0; o < O; ++o)
   {
      const real_t dt_sub       = deltaT_[o];
      const real_t weight       = time_weights_[o];
      const real_t t_sub_end    = t_sub_cursor + dt_sub;
      const real_t accum_scale  = weight * dt_macro;

      // 1. Inject one nucleation increment at the sub-step endpoint.
      //    The accumulator writes only `tau2_nuc` (R3-001 strike-slip
      //    invariant); the increment ΔS is monotonically non-negative
      //    (R3-004) and finite-dt-asserted (R3-007).
      ApplyNucleationIncremental_TPV104(dof_data, fault_coords,
                                        t_sub_end, dt_sub);

      // 2. Per-QP friction pipeline + ψ update + imposed-state accumulator.
      const bool last_sub_step = (o == O - 1);
      for (int i = 0; i < n; ++i)
      {
         DOFData      &d  = dof_data[i];
         const real_t *Ip = I_plus_flat  + static_cast<ptrdiff_t>(i) * NUM_STATE;
         const real_t *Im = I_minus_flat + static_cast<ptrdiff_t>(i) * NUM_STATE;

         // Time-average predictor -> pointwise Q̄ for friction solve.
         for (int c = 0; c < NUM_STATE; ++c)
         {
            Q_avg_plus[c]  = Ip[c] * inv_dt_macro;
            Q_avg_minus[c] = Im[c] * inv_dt_macro;
         }

         // Friction stage chain.  R5-003 (plan §4.10.X) closed by
         // expanding the `FrictionSolver::Method` enum with
         // `NewtonRaphsonStable`, so all four solvers (Brent, legacy
         // Newton, stable-asinh Newton, Hybrid) dispatch through the
         // same `ComputeStageState(..., method)` entry point.
         // `fault_face_flux.cpp` is untouched ([C2]).
         EvalStageState s;
         flux_.ComputeStageState(d, Q_avg_plus, Q_avg_minus, s, method);

#ifdef SEAS_DIAG_TPV104_STATE
         // Probe 1 (Stage A) — trial traction pre-friction.
         // Fields: t, qp_id, sigma_n_trial, tau1_trial, tau2_trial  [Pa].
         //
         // R4-004 deviation disclosure (review round 4): the MFEM
         // iterator evaluates the trial traction at Q̄ = I_±/dt_macro
         // for every sub-step o (see §4.10 Step 7 ambiguity resolution
         // in the header docstring).  The reference runtime evaluates
         // at per-sub-step qInterpolated[o].  Probe-1 rows emitted by
         // MFEM are therefore IDENTICAL across the O sub-steps of one
         // macro-step, while the reference emits O distinct rows.
         // The Step-13 Python probe-diff must mask (or coarsen by
         // sub-step-averaging) this channel when comparing across
         // codes until the iterator routes the per-sub-step ADER
         // predictor through the trial traction call (a Step-9
         // driver-level extension that requires editing the ADER
         // predictor interface — currently blocked by [C2]).
         {
            std::ofstream &f = GetProbeFile("trial_traction");
            if (f.is_open())
            {
               f << std::scientific << std::setprecision(16)
                 << (t_sub_cursor + dt_sub) << " " << i << " "
                 << s.sigma_n_trial << " "
                 << s.tau1_trial    << " "
                 << s.tau2_trial    << "\n";
            }
         }
         // Probe 3 (Stage C) — friction coefficient.
         // Fields: t, qp_id, V, psi, a, V0, mu.
         // R4-003 / R4-006 (review round 4): emit μ via the same
         // stable-asinh formula exposed for the Step-5 Newton solver,
         // NOT a naive inline `a·asinh(V·exp(ψ/a)/(2V₀))` that overflows
         // at ψ/a > ~700.  The probe now matches any downstream
         // FrictionCoefficientStable caller bit-for-bit.  When the
         // friction solver produced V_abs = 0 (e.g. rest state), μ is
         // defined as 0 by convention (no slip, no friction force) —
         // NaN propagation via asinh(0·exp(ψ/a)) is avoided either
         // way; the explicit guard makes the probe row unambiguous.
         {
            std::ofstream &f = GetProbeFile("friction_coeff");
            if (f.is_open())
            {
               const real_t V0_scalar = state_evo_.GetV0();
               const real_t mu = (s.V_abs > 0.0)
                  ? friction_stable::FrictionCoefficientStable(
                       s.V_abs, d.psi, d.a, V0_scalar)
                  : 0.0;
               f << std::scientific << std::setprecision(16)
                 << (t_sub_cursor + dt_sub) << " " << i << " "
                 << s.V_abs << " " << d.psi << " "
                 << d.a << " " << V0_scalar << " "
                 << mu << "\n";
            }
         }
         // Probe 4 (Stage D) — slip-rate magnitude.
         // Fields: t, qp_id, Theta, psi, sigma_n, eta_s, a, V_abs.
         {
            std::ofstream &f = GetProbeFile("slip_rate");
            if (f.is_open())
            {
               f << std::scientific << std::setprecision(16)
                 << (t_sub_cursor + dt_sub) << " " << i << " "
                 << s.Theta << " " << d.psi << " "
                 << s.sigma_n_total << " " << d.eta_s << " "
                 << d.a << " " << s.V_abs << "\n";
            }
         }
#endif  // SEAS_DIAG_TPV104_STATE

         // ψ update per sub-step (§3.12 directive — FVW analytic step).
         // Uses per-QP V_w[i] and per-QP d.a (NOT the base-virtual
         // default V_w_default_ — R-001 silent-fallthrough guard).
#ifdef SEAS_DIAG_TPV104_STATE
         const real_t psi_in_probe = d.psi;
#endif  // SEAS_DIAG_TPV104_STATE

         // R4-001 (review round 4): accumulate per-sub-step slip.
         // Each sub-step contributes  V · dt_sub  to the total slip in
         // both fault-tangent components.  Without this, `d.slip1` and
         // `d.slip2` remain at their init values (0 after
         // `InitializeFaultDOFs_TPV104`), and the Step-9 station writer
         // silently emits columns of zeros for `slip1 / slip2`.
         // WriteBackState below captures V1/V2 from the FINAL sub-step
         // only; relying on an external post-hoc  slip += V·dt_macro
         // would use that biased last-iterate V (not a time-average).
         // Accumulation here uses the sub-step V1/V2 that drove the
         // ψ integration, so slip and ψ evolve on the same time grid.
         d.slip1 += s.V1 * dt_sub;
         d.slip2 += s.V2 * dt_sub;

         d.psi = UpdateStateAnalyticSlipLawSRW(d.psi, s.V_abs,
                                               d.Dc, dt_sub,
                                               V_w[i], d.a,
                                               state_evo_.GetB(),
                                               state_evo_.GetV0(),
                                               state_evo_.GetF0(),
                                               state_evo_.GetMuW());
#ifdef SEAS_DIAG_TPV104_STATE
         // Probe 2 (Stage B) — state evolution.
         // Fields: t, qp_id, psi_in, V, L, dt, V_w, a, b, V0, f0, muW, psi_out.
         {
            std::ofstream &f = GetProbeFile("state_evolution");
            if (f.is_open())
            {
               f << std::scientific << std::setprecision(16)
                 << (t_sub_cursor + dt_sub) << " " << i << " "
                 << psi_in_probe << " " << s.V_abs << " "
                 << d.Dc << " " << dt_sub << " "
                 << V_w[i] << " " << d.a << " "
                 << state_evo_.GetB()   << " "
                 << state_evo_.GetV0()  << " "
                 << state_evo_.GetF0()  << " "
                 << state_evo_.GetMuW() << " "
                 << d.psi << "\n";
            }
         }
#endif  // SEAS_DIAG_TPV104_STATE

         // Build per-sub-step imposed state.  Pure function on (d, s,
         // Q_avg_plus, Q_avg_minus).
         flux_.BuildImposedState(d, s, Q_avg_plus, Q_avg_minus,
                                 Q_imp_plus, Q_imp_minus);

#ifdef SEAS_DIAG_TPV104_STATE
         // Probe 5 (Stage E) — corrected traction + imposed state.
         // Fields: t, qp_id, sigma_n_corr, tau1_corr, tau2_corr,
         //         Qimp_plus_vn, Qimp_plus_vt1, Qimp_plus_vt2,
         //         Qimp_minus_vn, Qimp_minus_vt1, Qimp_minus_vt2.
         {
            std::ofstream &f = GetProbeFile("corrected_imposed");
            if (f.is_open())
            {
               f << std::scientific << std::setprecision(16)
                 << (t_sub_cursor + dt_sub) << " " << i << " "
                 << s.sigma_n_corr << " "
                 << s.tau1_corr    << " "
                 << s.tau2_corr    << " "
                 << Q_imp_plus[VX]  << " "
                 << Q_imp_plus[VY]  << " "
                 << Q_imp_plus[VZ]  << " "
                 << Q_imp_minus[VX] << " "
                 << Q_imp_minus[VY] << " "
                 << Q_imp_minus[VZ] << "\n";
            }
         }
#endif  // SEAS_DIAG_TPV104_STATE

         // Accumulate time-weighted imposed state in time-integrated
         // form: I_imp_± += timeWeights[o] · dt_macro · Q_imp_±^{(o)}.
         real_t *Iout_p = I_imp_plus_flat
                          + static_cast<ptrdiff_t>(i) * NUM_STATE;
         real_t *Iout_m = I_imp_minus_flat
                          + static_cast<ptrdiff_t>(i) * NUM_STATE;
         for (int c = 0; c < NUM_STATE; ++c)
         {
            Iout_p[c] += accum_scale * Q_imp_plus[c];
            Iout_m[c] += accum_scale * Q_imp_minus[c];
         }

         // After the final sub-step: write slip_rate, V1, V2, tau*_corr,
         // sigma_n_corr onto DOFData so the probe format and downstream
         // station writers see the macro-step-terminal values.
         if (last_sub_step)
         {
            flux_.WriteBackState(d, s);
         }
      }

      t_sub_cursor = t_sub_end;
   }
}

// ---------------------------------------------------------------------------
// AdvanceWithSubStepStates — SeisSol-equivalent per-sub-step Q variant.
//
// Differs from Advance() at exactly one place: the per-sub-step working
// state Q_avg_plus / Q_avg_minus is taken from the per-sub-step input
// arrays instead of being computed as I/dt_macro.  Everything downstream
// (nucleation, ComputeStageState, ψ update, slip accumulation, imposed
// state, accumulator, WriteBackState) is bit-identical to Advance().
//
// Closes the R4-004 cadence-deviation note: trial traction is now
// evaluated at distinct per-sub-step Q values (matching SeisSol's
// `qInterpolated[o]` → `precomputeStressFromQInterpolated[o]`), not at a
// single time-averaged Q̄ replicated O times.
// ---------------------------------------------------------------------------
void Tpv104SubStepIterator::AdvanceWithSubStepStates(
   std::vector<DOFData> &dof_data,
   const std::vector<Vector> &fault_coords,
   const std::vector<real_t> &V_w,
   const std::vector<std::vector<real_t>> &Q_pointwise_plus_per_substep,
   const std::vector<std::vector<real_t>> &Q_pointwise_minus_per_substep,
   real_t dt_macro,
   real_t t_macro_start,
   real_t *I_imp_plus_flat,
   real_t *I_imp_minus_flat,
   FrictionSolver::Method method)
{
   if (deltaT_.empty())
   {
      throw std::runtime_error(
         "Tpv104SubStepIterator::AdvanceWithSubStepStates: SetSubSteps "
         "has not been called; cannot iterate without a configured "
         "quadrature.");
   }
   if (!std::isfinite(dt_macro) || dt_macro <= 0.0)
   {
      throw std::runtime_error(
         "Tpv104SubStepIterator::AdvanceWithSubStepStates: dt_macro "
         "must be finite and positive; got " + std::to_string(dt_macro));
   }
   if (!std::isfinite(t_macro_start))
   {
      throw std::runtime_error(
         "Tpv104SubStepIterator::AdvanceWithSubStepStates: t_macro_start "
         "must be finite; got " + std::to_string(t_macro_start));
   }
   if (I_imp_plus_flat == nullptr || I_imp_minus_flat == nullptr)
   {
      throw std::runtime_error(
         "Tpv104SubStepIterator::AdvanceWithSubStepStates: I_imp_*_flat "
         "must both be non-null.");
   }

   const int O = static_cast<int>(deltaT_.size());
   if (static_cast<int>(Q_pointwise_plus_per_substep.size()) != O ||
       static_cast<int>(Q_pointwise_minus_per_substep.size()) != O)
   {
      throw std::runtime_error(
         "Tpv104SubStepIterator::AdvanceWithSubStepStates: Q_pointwise_*"
         "_per_substep must each have size O = " + std::to_string(O)
         + "; got plus = "
         + std::to_string(Q_pointwise_plus_per_substep.size())
         + ", minus = "
         + std::to_string(Q_pointwise_minus_per_substep.size()));
   }

   const int n = static_cast<int>(dof_data.size());
   if (static_cast<int>(fault_coords.size()) != n
       || static_cast<int>(V_w.size()) != n)
   {
      throw std::runtime_error(
         "Tpv104SubStepIterator::AdvanceWithSubStepStates: dof_data, "
         "fault_coords, V_w must have equal size; got "
         + std::to_string(dof_data.size()) + ", "
         + std::to_string(fault_coords.size()) + ", "
         + std::to_string(V_w.size()));
   }

   // Verify each per-sub-step Q has length NUM_STATE * n.
   const size_t expected_words =
      static_cast<size_t>(NUM_STATE) * static_cast<size_t>(n);
   for (int o = 0; o < O; o++)
   {
      if (Q_pointwise_plus_per_substep[o].size() != expected_words ||
          Q_pointwise_minus_per_substep[o].size() != expected_words)
      {
         throw std::runtime_error(
            "Tpv104SubStepIterator::AdvanceWithSubStepStates: "
            "Q_pointwise_*[" + std::to_string(o) + "] must have size "
            "NUM_STATE * n = " + std::to_string(expected_words)
            + "; got plus = "
            + std::to_string(Q_pointwise_plus_per_substep[o].size())
            + ", minus = "
            + std::to_string(Q_pointwise_minus_per_substep[o].size()));
      }
   }

   // Verify Σ deltaT == dt_macro (same tolerance policy as Advance).
   const real_t dtsum = std::accumulate(deltaT_.begin(), deltaT_.end(),
                                        static_cast<real_t>(0));
   const real_t rel = std::abs(dtsum - dt_macro) / std::max(dt_macro,
                                                            1e-300);
   const real_t sum_tol = std::max<real_t>(
      1e-12,
      static_cast<real_t>(10.0) * O
         * std::numeric_limits<real_t>::epsilon());
   if (rel > sum_tol)
   {
      throw std::runtime_error(
         "Tpv104SubStepIterator::AdvanceWithSubStepStates: Σ deltaT "
         "must equal dt_macro within " + std::to_string(sum_tol)
         + " rel; got Σ = " + std::to_string(dtsum)
         + ", dt_macro = " + std::to_string(dt_macro));
   }

   // Zero the output accumulators.
   const size_t nwords = expected_words;
   std::memset(I_imp_plus_flat,  0, nwords * sizeof(real_t));
   std::memset(I_imp_minus_flat, 0, nwords * sizeof(real_t));

   real_t Q_imp_plus[NUM_STATE];
   real_t Q_imp_minus[NUM_STATE];

   real_t t_sub_cursor = t_macro_start;

   for (int o = 0; o < O; ++o)
   {
      const real_t dt_sub      = deltaT_[o];
      const real_t weight      = time_weights_[o];
      const real_t t_sub_end   = t_sub_cursor + dt_sub;
      const real_t accum_scale = weight * dt_macro;

      // 1. Per-sub-step nucleation increment (same call as Advance).
      ApplyNucleationIncremental_TPV104(dof_data, fault_coords,
                                        t_sub_end, dt_sub);

      // 2. Per-QP friction pipeline + ψ + accumulator.  Sub-step Q comes
      //    from the predictor's per-sub-step pointwise output, not from
      //    I/dt_macro.
      const bool last_sub_step = (o == O - 1);
      const real_t *Qp_o = Q_pointwise_plus_per_substep[o].data();
      const real_t *Qm_o = Q_pointwise_minus_per_substep[o].data();

      for (int i = 0; i < n; ++i)
      {
         DOFData &d = dof_data[i];

         // Per-QP working states Q̃_± at sub-step time τ_o.
         const real_t *Q_tilde_plus  =
            Qp_o + static_cast<ptrdiff_t>(i) * NUM_STATE;
         const real_t *Q_tilde_minus =
            Qm_o + static_cast<ptrdiff_t>(i) * NUM_STATE;


         EvalStageState s;
         flux_.ComputeStageState(d, Q_tilde_plus, Q_tilde_minus, s, method);

         // R4-001: per-sub-step slip accumulation (matches Advance).
         d.slip1 += s.V1 * dt_sub;
         d.slip2 += s.V2 * dt_sub;

         d.psi = UpdateStateAnalyticSlipLawSRW(d.psi, s.V_abs,
                                               d.Dc, dt_sub,
                                               V_w[i], d.a,
                                               state_evo_.GetB(),
                                               state_evo_.GetV0(),
                                               state_evo_.GetF0(),
                                               state_evo_.GetMuW());

         flux_.BuildImposedState(d, s, Q_tilde_plus, Q_tilde_minus,
                                 Q_imp_plus, Q_imp_minus);

         real_t *Iout_p = I_imp_plus_flat
                          + static_cast<ptrdiff_t>(i) * NUM_STATE;
         real_t *Iout_m = I_imp_minus_flat
                          + static_cast<ptrdiff_t>(i) * NUM_STATE;
         for (int c = 0; c < NUM_STATE; ++c)
         {
            Iout_p[c] += accum_scale * Q_imp_plus[c];
            Iout_m[c] += accum_scale * Q_imp_minus[c];
         }

         if (last_sub_step)
         {
            flux_.WriteBackState(d, s);
         }
      }

      t_sub_cursor = t_sub_end;
   }
}

} // namespace seas
} // namespace mfem

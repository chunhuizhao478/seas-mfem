// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// SCEC TPV102 benchmark driver — 3D dynamic rupture on a vertical
// strike-slip fault with the regularised rate-and-state ageing law
// (SCEC TPV101/102 spec, ageing law in ψ-space).
//
// Fully parallel MPI driver using ParMesh and WaveOperator<ParMesh>.
// Mirrors drivers/tpv104_driver.cpp's ADER time-stepping flow with the
// only-change-is-the-friction-law substitutions:
//   - TPV102Params material / friction constants (matches PDF spec).
//   - InitializeFaultDOFs (per-QP a(x,z) via SCEC boxcar, ψ_ini from
//     equilibrium inversion).
//   - No V_w side-channel — TPV102 has no weakening velocity.
//   - ApplyNucleationIncremental_TPV102 (cumulative per-sub-step
//     accumulator into tau2_nuc).
//   - UpdateStateAnalytic (the exact ageing-law analytic ψ update from
//     friction/state_evolution.hpp; SCEC Eq. (2): dθ/dt = 1 − Vθ/L).
//   - AgingLawPsi (not SlipLawSRWPsi) feeds the substep iterator.
//   - TPV102StationWriter / TPV102SurfaceStationWriter from the existing
//     dynamic/tpv102_setup.hpp (9 fault stations + 6 free-surface).
//
// Usage:
//   ibrun ./seas_tpv102_driver \
//         --mesh tpv102/mesh/tpv102_200m.msh \
//         --tfinal 12.0 --ader-order 2 \
//         --friction-solver newton-stable \
//         --output-dir tpv102/results

#include "mfem.hpp"
#include "../dynamic/wave_state.hpp"
#include "../dynamic/wave_operator.hpp"
#include "../dynamic/fault_face_flux.hpp"
#include "../dynamic/friction_solver.hpp"
#include "../dynamic/tpv102_setup.hpp"
#include "../dynamic/tpv102_nucleation.hpp"
#include "../dynamic/tpv102_substep_iterator.hpp"
#include "../config/tpv102_params.hpp"
#include "../domain/boundary_config.hpp"
#include "../friction/state_evolution.hpp"
#include "../dynamic/seas_diag_rank.hpp"
#include "../dynamic/fault_locality_partition.hpp"
#include "../io/paraview_output.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <limits>
#include <climits>
#include <memory>
#include <string>
#include <vector>
#include <sys/stat.h>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

#ifdef SEAS_DIAG_FAULT_FLUX
// Global rank cache declared in dynamic/seas_diag_rank.hpp — one
// definition for the whole link.  Same pattern as tpv102_driver.cpp.
namespace mfem { namespace seas { int g_seas_my_rank = 0; } }  // NOLINT
#endif

// --------------------------------------------------------------------------
// Small CLI parsing helpers (same convention as tpv102_driver.cpp so
// launch scripts can share argument shapes).
// --------------------------------------------------------------------------
static std::string GetStringArg(int argc, char *argv[], const char *flag,
                                const std::string &default_val)
{
   for (int i = 1; i < argc - 1; ++i)
   {
      if (std::string(argv[i]) == flag) { return argv[i + 1]; }
   }
   return default_val;
}

static real_t GetRealArg(int argc, char *argv[], const char *flag,
                         real_t default_val)
{
   const std::string val = GetStringArg(argc, argv, flag, "");
   if (val.empty()) { return default_val; }
   return std::stod(val);
}

static int GetIntArg(int argc, char *argv[], const char *flag,
                     int default_val)
{
   const std::string val = GetStringArg(argc, argv, flag, "");
   if (val.empty()) { return default_val; }
   return std::stoi(val);
}

static bool HasFlag(int argc, char *argv[], const char *flag)
{
   for (int i = 1; i < argc; ++i)
   {
      if (std::string(argv[i]) == flag) { return true; }
   }
   return false;
}

// Map the --friction-solver CLI name to FrictionSolver::Method (plan §4.10.X).
//
// R7-001/R7-005 note: on the current driver path this value is kept only
// for future iterator wiring.  The production time loop runs Brent via
// wave.AdvanceADER -> FaultFaceFlux::EvaluateADER (fluctuation-Q
// dispatch — wave_operator.inl:3614); the returned Method is not
// routed through that call.
// Nevertheless, MapSolver is kept strict so that (a) `newton` (bare) is
// the canonical shorthand for the stable-asinh variant and (b) unknown
// strings abort loudly rather than silently defaulting — both become
// load-bearing the moment the iterator is wired into the time loop.
static FrictionSolver::Method MapSolver(const std::string &s)
{
   if (s == "newton-stable" || s == "newton")
   {
      return FrictionSolver::Method::NewtonRaphsonStable;
   }
   if (s == "brent")        { return FrictionSolver::Method::Brent; }
   if (s == "newton-legacy"){ return FrictionSolver::Method::NewtonRaphson; }
   if (s == "hybrid")       { return FrictionSolver::Method::HybridNRBisection; }
   MFEM_ABORT("--friction-solver: unknown value '" << s
              << "'.  Accepted: newton-stable | newton | brent | "
              << "newton-legacy | hybrid.");
}

static std::string SolverBanner(const std::string &s)
{
   if (s == "newton-stable" || s == "newton")
   {
      return "Newton-Raphson (stable-asinh, plan §4.10 Step 5)";
   }
   if (s == "brent")        { return "Brent (log10-V, Tandem-verified; legacy dispatch)"; }
   if (s == "newton-legacy"){ return "Newton-Raphson (legacy, MFEM-native μ)"; }
   if (s == "hybrid")       { return "Hybrid NR+Bisection (legacy MFEM μ)"; }
   return std::string("unknown (") + s + ")";
}

// R8-001: route the dispatched-classification tags through an enum +
// `switch` rather than hardcoded constants.  Today — under R7-001
// option (b) — GetDispatched*(cli) unconditionally returns the
// Brent / OneShot / SlipSRW enum value regardless of CLI, because the
// driver's time loop ignores the CLI solver / iterator flags.  When
// R7-001 option (a) wires the iterator in, updating the Get*
// functions to branch on `cli` is the single required edit; TagOf() /
// BannerOf() stay correct because they're derived from the enum.
// Forgetting to extend those switches for a new enum value is a
// compiler warning (-Wswitch), preventing hardcoded-tag drift.
enum class DispatchedSolver { Brent, NewtonRaphsonStable, NewtonRaphsonLegacy, Hybrid };
enum class DispatchedIterator { OneShot, SubStep };
enum class DispatchedLaw { SlipSRW, Aging };

static DispatchedSolver GetDispatchedSolver(const std::string &/*friction_solver_cli*/)
{
   // R7-001 option (b): wave.AdvanceADER -> EvaluateADER hard-codes
   // the default Method::Brent argument; the CLI value does not reach
   // the solver dispatch.  Replace with a CLI-aware switch when the
   // iterator is wired (option a).
   return DispatchedSolver::Brent;
}
static DispatchedIterator GetDispatchedIterator(const std::string &fault_iterator_cli)
{
   // R-602/R-603 (round-7 implementation): the CLI value now controls
   // dispatch.  "substep" routes the time loop through the
   // Tpv102SubStepIterator + per-sub-step ADER predictor path; any other
   // value (including "one-shot") keeps the legacy single-shot
   // wave.AdvanceADER call.  Default unchanged: one-shot.
   if (fault_iterator_cli == "substep")
   {
      return DispatchedIterator::SubStep;
   }
   return DispatchedIterator::OneShot;
}
static DispatchedLaw GetDispatchedLaw(const std::string &/*fric_law_cli*/)
{
   return DispatchedLaw::Aging;
}

static const char *TagOf(DispatchedSolver s)
{
   switch (s)
   {
      case DispatchedSolver::Brent:                return "brent";
      case DispatchedSolver::NewtonRaphsonStable:  return "newton-stable";
      case DispatchedSolver::NewtonRaphsonLegacy:  return "newton-legacy";
      case DispatchedSolver::Hybrid:               return "hybrid";
   }
   return "unknown";
}
static const char *TagOf(DispatchedIterator i)
{
   switch (i)
   {
      case DispatchedIterator::OneShot: return "oneshot";
      case DispatchedIterator::SubStep: return "substep";
   }
   return "unknown";
}
static const char *TagOf(DispatchedLaw l)
{
   switch (l)
   {
      case DispatchedLaw::SlipSRW: return "slip-srw";
      case DispatchedLaw::Aging:   return "aging";
   }
   return "unknown";
}

static std::string BannerOf(DispatchedSolver s)
{
   switch (s)
   {
      case DispatchedSolver::Brent:
         return "Brent (hard-coded via EvaluateADER fluctuation-Q "
                "dispatch; --friction-solver flag IGNORED)";
      case DispatchedSolver::NewtonRaphsonStable:
         return "Newton-Raphson (stable-asinh, plan §4.10 Step 5)";
      case DispatchedSolver::NewtonRaphsonLegacy:
         return "Newton-Raphson (legacy, MFEM-native μ)";
      case DispatchedSolver::Hybrid:
         return "Hybrid NR+Bisection (legacy MFEM μ)";
   }
   return "unknown";
}
static std::string BannerOf(DispatchedIterator i)
{
   switch (i)
   {
      case DispatchedIterator::OneShot:
         return "one-shot (default; legacy wave.AdvanceADER dispatch)";
      case DispatchedIterator::SubStep:
         return "sub-step (Tpv102SubStepIterator + per-sub-step ADER "
                "predictor — round-7 R-602/R-603, opt-in via "
                "--fault-iterator substep)";
   }
   return "unknown";
}
static std::string BannerOf(DispatchedLaw l)
{
   switch (l)
   {
      case DispatchedLaw::SlipSRW:
         return "slip-SRW (ψ-space, macro-step analytic "
                "UpdateStateAnalyticSlipLawSRW)";
      case DispatchedLaw::Aging:
         return "aging (ψ-space, macro-step analytic UpdateStateAnalytic; "
                "SCEC TPV101/102 dθ/dt = 1 − Vθ/L)";
   }
   return "unknown";
}

// ---------------------------------------------------------------------------
// AdvanceADERWithSubStep — round-7 implementation of the SeisSol-equivalent
// per-sub-step ADER dispatch.  Composes:
//
//   1. wave.ComputeADERSubStepStates(Q, dt, order, tau_nodes) → Q_per_node[o]
//   2. wave.EvaluateBulkAtFaultQPsCanonical(Q_per_node[o]) for each o →
//        Q_pointwise_plus_per_substep, Q_pointwise_minus_per_substep
//   3. iterator.AdvanceWithSubStepStates(...) →
//        accumulated I_imp_plus_flat, I_imp_minus_flat, plus DOFData updates
//   4. wave.SetSubStepFaultImposedStates(...) — the fault branch of the
//        upcoming AdvanceADER call will consume these.
//   5. wave.AdvanceADER(...) — bulk corrector runs as usual but the fault
//        branch substitutes the iterator's I_imp instead of running
//        EvaluateADER inline.
//   6. wave.ResetSubStepFaultImposedStates() — restore default for safety.
//
// Sub-step quadrature: `iterator.GetDeltaT()` and `iterator.GetTimeWeights()`
// must be configured by SetSubSteps before this is called.  The driver
// passes the cumulative-prefix nodes
//   tau_nodes[o] = Σ_{o'<=o} deltaT[o']
// (the SUB-STEP END NODES on [0, dt]).  This matches the cadence the
// existing iterator uses for nucleation endpoints (`t_sub_end`) and ψ
// updates (`dt_sub` per sub-step), so trial traction at sub-step `o`
// is now consistent with the rest of that sub-step's bookkeeping.
//
// Returns 0 on success.  Aborts via MFEM_ABORT on contract violations
// (predictor / iterator size mismatches).
// ---------------------------------------------------------------------------
template <typename WaveOpT>
static void AdvanceADERWithSubStep(
   WaveOpT &wave,
   mfem::seas::Tpv102SubStepIterator &iterator,
   std::vector<mfem::seas::DOFData> &dof_data,
   const std::vector<mfem::Vector> &fault_coords,
   const mfem::Vector &Q,
   mfem::real_t dt_step,
   int ader_order,
   mfem::real_t t_step_start,
   mfem::seas::FrictionSolver::Method method,
   mfem::Vector &Q_new)
{
   using mfem::real_t;
   using mfem::Vector;

   MFEM_VERIFY(dt_step > 0.0,
               "AdvanceADERWithSubStep: dt_step must be > 0, got "
               << dt_step);
   MFEM_VERIFY(ader_order >= 2 && ader_order <= 4,
               "AdvanceADERWithSubStep: ader_order must be in {2,3,4}, "
               "got " << ader_order);

   // R-1002: read configured deltaT/weights ONCE up front; then rescale
   // deltaT per call so the iterator's Σ deltaT == dt_step verify holds
   // even when the time loop's dt_step varies (e.g., the final time
   // step where dt_step = tfinal - t < auto-CFL dt).  The configured
   // ratios deltaT[o]/Σ deltaT are preserved; weights stay unchanged.
   const std::vector<real_t> configured_deltaT = iterator.GetDeltaT();
   const std::vector<real_t> configured_weights = iterator.GetTimeWeights();
   const int O = static_cast<int>(configured_deltaT.size());
   MFEM_VERIFY(O >= 1,
               "AdvanceADERWithSubStep: iterator has empty deltaT; "
               "SetSubSteps must be called before dispatching this path.");
   const mfem::real_t configured_sum =
      std::accumulate(configured_deltaT.begin(), configured_deltaT.end(),
                      static_cast<mfem::real_t>(0));
   MFEM_VERIFY(configured_sum > 0.0,
               "AdvanceADERWithSubStep: configured Σ deltaT = "
               << configured_sum << " ≤ 0");

   const mfem::real_t dt_scale = dt_step / configured_sum;
   std::vector<mfem::real_t> deltaT_scaled(O);
   for (int o = 0; o < O; o++)
   {
      deltaT_scaled[o] = configured_deltaT[o] * dt_scale;
   }
   iterator.SetSubSteps(deltaT_scaled, configured_weights);
   const std::vector<real_t> &deltaT = iterator.GetDeltaT();   // = deltaT_scaled

   // R-1001: SUB-STEP MIDPOINT nodes on [0, dt_step].  tau_nodes[o] is
   // the midpoint of sub-step o relative to the macro-step start.  For
   // ADER-2 predictor (Q linear in τ), the midpoint Q(τ_o) equals the
   // sub-step's time-average — restoring the T_TPV102_SSI_3 contract
   // (substep at O=1 with deltaT={dt}, weights={1.0} == one-shot at
   // O=1).  Pre-fix: cumulative-end nodes (τ_0=dt) gave Q(dt) instead
   // of Q̄=I/dt, breaking the SSI_3 bit-equivalence.  For O ≥ 3,
   // midpoint rule is O(dt²)-accurate; full SeisSol parity (Gauss-
   // Lobatto) is R-1005, separate.
   std::vector<real_t> tau_nodes(O);
   real_t acc = 0.0;
   for (int o = 0; o < O; o++)
   {
      tau_nodes[o] = acc + 0.5 * deltaT[o];
      acc += deltaT[o];
   }

   // Predictor: per-sub-step pointwise Q in the bulk.
   std::vector<Vector> Q_per_node;
   wave.ComputeADERSubStepStates(Q, dt_step, ader_order, tau_nodes,
                                 Q_per_node);
   MFEM_VERIFY(static_cast<int>(Q_per_node.size()) == O,
               "AdvanceADERWithSubStep: ComputeADERSubStepStates returned "
               << Q_per_node.size() << " nodes, expected " << O);

   // Per-sub-step Q at fault QPs in canonical frame.
   //
   // R-1003: sized to GetNumTotalFaultQPs() (interior + shared) so the
   // shared-fault iterator slice has storage.  The iterator's verify at
   // tpv102_substep_iterator.cpp:573-589 expects each Q_pointwise_*[o]
   // and the I_imp_*_flat output to be sized NUM_STATE * dof_data.size(),
   // and dof_data is sized num_fault_total = local + shared (driver L1015).
   // At np=1, shared == 0 so total == local and the buffer is byte-identical
   // to the pre-R-1003 sizing.
   const int n_total_fault_qps = wave.GetNumTotalFaultQPs();
   std::vector<std::vector<real_t>> Q_pointwise_plus(O), Q_pointwise_minus(O);
   for (int o = 0; o < O; o++)
   {
      wave.EvaluateBulkAtFaultQPsCanonical(Q_per_node[o],
                                           Q_pointwise_plus[o],
                                           Q_pointwise_minus[o]);
   }

   // Iterator: per-sub-step friction + ψ + slip + accumulator.
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
                                        I_imp_minus_flat.data(),
                                        method);
   }

   // Hand the iterator's output to the wave op so the upcoming
   // AdvanceADER's fault branch substitutes it for inline EvaluateADER.
   wave.SetSubStepFaultImposedStates(
      n_total_fault_qps > 0 ? I_imp_plus_flat.data()  : nullptr,
      n_total_fault_qps > 0 ? I_imp_minus_flat.data() : nullptr,
      n_total_fault_qps);

   // Bulk corrector: runs unchanged for non-fault faces; fault branch
   // consumes the side-channel imposed states.
   wave.AdvanceADER(Q, dt_step, ader_order, Q_new);

   // Restore default behavior for any subsequent direct AdvanceADER call.
   wave.ResetSubStepFaultImposedStates();
}

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
   MPI_Init(&argc, &argv);
   MPI_Comm comm = MPI_COMM_WORLD;
   int rank, nprocs;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nprocs);
#else
   int rank = 0, nprocs = 1;
#endif

#ifdef SEAS_DIAG_FAULT_FLUX
   mfem::seas::g_seas_my_rank = rank;
#endif

   if (rank == 0)
   {
#ifdef SEAS_DIAG_FAULT_FLUX
      const char *diag_fault_flux = "SEAS_DIAG_FAULT_FLUX = ON";
#else
      const char *diag_fault_flux = "SEAS_DIAG_FAULT_FLUX = OFF";
#endif
      // SEAS_DIAG_TPV102_STATE banner intentionally omitted: the TPV102
      // substep iterator does not implement the per-sub-step probe
      // writers (no parity with tpv104_substep_iterator.cpp's
      // GetProbeFile path).  Advertising ON/OFF here would be a false
      // contract.  Re-add only when the iterator actually emits the
      // probe traces.
      std::fprintf(stderr, "[BUILD] %s\n", diag_fault_flux);
      std::ofstream binfo("build_info.txt");
      if (binfo.is_open())
      {
         binfo << "[BUILD] " << diag_fault_flux << "\n";
         binfo.close();
      }
   }

   // -----------------------------------------------------------------------
   // Parse command-line arguments
   // -----------------------------------------------------------------------
   std::string mesh_file = GetStringArg(argc, argv, "--mesh",
                                        "tpv102/mesh/tpv102_200m.msh");
   real_t mesh_scale     = GetRealArg(argc, argv, "--mesh-scale", 1.0);
   int order             = GetIntArg(argc, argv, "--order", 1);
   std::string bc_mode   = GetStringArg(argc, argv, "--bc-mode", "absorbing");
   real_t tfinal         = GetRealArg(argc, argv, "--tfinal", TPV102Params::t_final);
   std::string output_dir    = GetStringArg(argc, argv, "--output-dir", "tpv102/results");
   std::string output_prefix = GetStringArg(argc, argv, "--output-prefix", "tpv102");
   real_t cfl_factor     = GetRealArg(argc, argv, "--cfl", 0.5);
   int bc_free           = GetIntArg(argc, argv, "--bc-free", 1);
   int bc_fault          = GetIntArg(argc, argv, "--bc-fault", 3);
   int bc_absorb         = GetIntArg(argc, argv, "--bc-absorb", 5);
   int ader_order        = GetIntArg(argc, argv, "--ader-order", 2);
   real_t dt_override    = GetRealArg(argc, argv, "--dt", 0.0);
   real_t output_dt      = GetRealArg(argc, argv, "--output-dt", 0.01);
   bool disable_nucleation = HasFlag(argc, argv, "--disable-nucleation");
   bool debug_qnorm      = HasFlag(argc, argv, "--debug-qnorm");
   bool dry_run          = HasFlag(argc, argv, "--dry-run");
   bool verify_dispatch  = HasFlag(argc, argv, "--verify-dispatch");
   bool fault_locality_part = HasFlag(argc, argv, "--partition-fault-locality");
   std::string partition_file = GetStringArg(argc, argv, "--partition-file", "");

   // ParaView output controls — mirror tpv102_driver.cpp + BP5 conventions:
   //   --paraview              : enable PVD/VTU output, interval matches --output-dt
   //   --paraview-every N      : write every N steps
   //   --paraview-dt X         : write every X seconds (overrides step interval)
   //   --paraview-bulk-dt X    : enable a SECOND collection in ParaView_bulk/
   //                             with velocity + sigma_yy/sigma_xy/sigma_xz at
   //                             coarser cadence (typical: 0.05 s)
   //   --pv-low-order          : linear tets only (~40x smaller volume output)
   //   --no-domain-pv          : suppress fault-schedule volume save (fault-
   //                             surface PVD/VTU still written; bulk collection
   //                             unaffected)
   bool use_paraview = HasFlag(argc, argv, "--paraview");
   bool pv_low_order = HasFlag(argc, argv, "--pv-low-order");
   bool pv_no_domain = HasFlag(argc, argv, "--no-domain-pv");
   int  paraview_step_interval = GetIntArg(argc, argv, "--paraview-every", 0);
   real_t paraview_dt_flag     = GetRealArg(argc, argv, "--paraview-dt", 0.0);
   real_t paraview_bulk_dt     = GetRealArg(argc, argv, "--paraview-bulk-dt", 0.0);
   if (paraview_step_interval > 0 || paraview_dt_flag > 0.0
       || paraview_bulk_dt > 0.0)
   {
      use_paraview = true;
   }
   // `--dry-run` is a shortcut for "no mesh, no time-stepping,
   // just print banner + verify wiring compiles/runs".  Used by
   // test_tpv102_smoke.cpp and the banner-check sbatch on Frontera.
   // `--verify-dispatch` (R7-004) additionally emits machine-readable
   // [dispatch] lines describing what the production time loop ACTUALLY
   // runs — the tri-consistency check between banner claims and runtime
   // dispatch.  Under R7-001 option (b) the dispatch is always
   // Brent/one-shot/slip-SRW, independent of --friction-solver /
   // --fault-iterator / --fric-law values.

   // TPV102-specific CLI (plan §4.10 Step 9).
   // R7-001/R7-003/R7-006: these flags are accepted so smoke tests and
   // sbatch scripts can exercise banner parity, but on the current
   // driver path (one-shot wave.AdvanceADER) they have NO effect on the
   // dispatched solver / iterator / friction law.  Runtime is:
   //   - friction solver: Brent (hard-coded via EvaluateADER)
   //   - fault iterator : one-shot (Tpv102SubStepIterator not wired
   //                      — requires exposing per-sub-step I± from
   //                      wave_operator.inl; that file is on the
   //                      extreme-care no-touch list)
   //   - friction law   : slip-SRW via the free function
   //                      UpdateStateAnalyticSlipLawSRW (per-macro-step)
   // The banner below mirrors this disclosure exactly.
   std::string friction_solver =
      GetStringArg(argc, argv, "--friction-solver", "newton-stable");
   std::string fric_law =
      GetStringArg(argc, argv, "--fric-law", "aging");
   // Default = "one-shot": legacy wave.AdvanceADER dispatch (production
   // path, byte-identical to pre-R-602 behavior).  Pass --fault-iterator
   // substep to opt in to the per-sub-step Tpv102SubStepIterator path.
   // (Pre-R-602 the default was "substep" but DispatchedIterator unconditionally
   //  returned OneShot; the string was banner-only.  Now that
   //  GetDispatchedIterator routes the string, the default has to flip
   //  to keep production behavior unchanged.)
   std::string fault_iterator =
      GetStringArg(argc, argv, "--fault-iterator", "one-shot");

   // R-1601: validate `--fault-iterator` value at parse time.  Without
   // this, typos (e.g., `sub-step` with a hyphen, `SUBSTEP` mixed case)
   // silently fall through to one-shot AND bypass the R-1503 guard
   // (`fault_iterator == "substep"` exact-match), defeating both the
   // user's intended dispatch routing AND the safety net.  Mirror the
   // loud-abort pattern used by `--mixed-flux` and `--friction-solver`.
   if (fault_iterator != "one-shot" && fault_iterator != "substep")
   {
      MFEM_ABORT("--fault-iterator: unknown value '" << fault_iterator
                 << "'.  Accepted: one-shot | substep.");
   }

   // Round-11 Mixed-Flux dispatch (Zhang et al. 2023, MIXED_FLUX_PLAN.md).
   // Default = "none": upwind everywhere, byte-identical to pre-Mixed-Flux
   // behavior.  Accepted values:
   //   - "none"           upwind everywhere (default)
   //   - "adjacent"       central flux on faces adjacent to fault (Mixed-Flux 2)
   //   - "all-continuous" central on every interior non-fault face (Mixed-Flux 1)
   std::string mixed_flux_str =
      GetStringArg(argc, argv, "--mixed-flux", "none");
   MixedFluxMode mixed_flux_mode = MixedFluxMode::None;
   if      (mixed_flux_str == "none")           { mixed_flux_mode = MixedFluxMode::None; }
   else if (mixed_flux_str == "adjacent")       { mixed_flux_mode = MixedFluxMode::Adjacent; }
   else if (mixed_flux_str == "all-continuous") { mixed_flux_mode = MixedFluxMode::AllContinuous; }
   else
   {
      MFEM_ABORT("--mixed-flux: unknown value '" << mixed_flux_str
                 << "'.  Accepted: none | adjacent | all-continuous.");
   }

   // R-1105 / Phase 6 §2: driver-level fast-fail mutual-exclusion guard
   // for `--use-precomputed-face-fluxes` × `--mixed-flux != none`.  The
   // wave-operator-level guard at SetMixedFluxMode (wave_operator.inl
   // L1160-1167) catches the combination, but only AFTER mesh
   // construction and ParMesh distribution — wasting minutes of work
   // on a doomed run at production scale (e.g., 1.5M tets at np=128).
   // This driver-level check fires at CLI-parse time, before the mesh
   // is read.  TPV102 currently doesn't expose
   // `--use-precomputed-face-fluxes` (the precomputed-flux path is a
   // TPV102 opt-in); the guard is defensive — kicks in the moment a
   // future driver enhancement adds the flag.
   if (HasFlag(argc, argv, "--use-precomputed-face-fluxes") &&
       mixed_flux_mode != MixedFluxMode::None)
   {
      if (rank == 0)
      {
         std::cerr
            << "[FATAL] --use-precomputed-face-fluxes is mutually "
            "exclusive with --mixed-flux != none.  The precomputed-flux "
            "path bakes upwind dispatch into its tables; mixed-flux "
            "would be silently ignored.  Pick one or the other.\n";
      }
#ifdef MFEM_USE_MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
#else
      std::abort();
#endif
   }

   // R-1204: --mixed-flux adjacent has no test coverage at --ader-order > 2.
   // Plan §Risk R5 flagged this; until a higher-order MPI gate lands, abort
   // on the unvalidated combination.  Override via
   //   SEAS_FORCE_MIXED_FLUX_ADER_O_GT2=1
   // for experimental runs.  This guard is at CLI-parse time, before mesh
   // construction.
   if (mixed_flux_mode != MixedFluxMode::None && ader_order > 2)
   {
      const char *force = std::getenv("SEAS_FORCE_MIXED_FLUX_ADER_O_GT2");
      if (!(force && force[0] == '1'))
      {
         if (rank == 0)
         {
            std::cerr
               << "[FATAL] --mixed-flux " << mixed_flux_str
               << " --ader-order " << ader_order
               << ": untested combination (R-1204, plan §Risk R5).  "
               "Mixed-flux dispatch was validated at ADER-O2 only.  "
               "Set SEAS_FORCE_MIXED_FLUX_ADER_O_GT2=1 to override "
               "for experimental runs.\n";
         }
#ifdef MFEM_USE_MPI
         MPI_Abort(MPI_COMM_WORLD, 1);
#else
         std::abort();
#endif
      }
      else if (rank == 0)
      {
         // R-1407: override exercised — emit a loud warning so a user
         // who set the env var in their shell rc can see they're
         // running unverified code paths.  No silent bypass.
         std::cerr
            << "[WARNING] SEAS_FORCE_MIXED_FLUX_ADER_O_GT2=1 — running "
            "--mixed-flux " << mixed_flux_str
            << " --ader-order " << ader_order
            << " is UNTESTED.  Numerical correctness is NOT guaranteed.  "
            "If a regression is observed, cite R-1204/R-1407 in the "
            "report.  Disable the override (`unset "
            "SEAS_FORCE_MIXED_FLUX_ADER_O_GT2`) for production runs.\n";
      }
   }

   // R-1503 + R-1605 / Plan §Risk R5: substep iterator + mixed-flux has
   // NO test coverage at ANY rank count.  The substep dispatch path
   // interleaves per-substep fault-state setting with the predictor /
   // corrector; a subtle ordering bug interacting with the mixed-flux
   // central dispatch could produce wrong rupture-front velocities.
   // The guard fires on np >= 1 (R-1605: np=1 is also uncovered — the
   // local-reproducer configuration a developer would use; see plan
   // §Risk R5).  Override via
   //   SEAS_FORCE_MIXED_FLUX_SUBSTEP_MPI=1
   // for experimental runs (loud warning when used).
   if (mixed_flux_mode != MixedFluxMode::None &&
       fault_iterator == "substep")
   {
      const char *force =
         std::getenv("SEAS_FORCE_MIXED_FLUX_SUBSTEP_MPI");
      if (!(force && force[0] == '1'))
      {
         if (rank == 0)
         {
            std::cerr
               << "[FATAL] --mixed-flux " << mixed_flux_str
               << " --fault-iterator substep at np=" << nprocs
               << ": untested combination (R-1503/R-1605, Plan §Risk "
               "R5).  No test exercises substep iterator + mixed-flux "
               "dispatch at ANY np (np=1 included; the local-reproducer "
               "configuration is also uncovered).  Set "
               "SEAS_FORCE_MIXED_FLUX_SUBSTEP_MPI=1 to override for "
               "experimental runs.\n";
         }
#ifdef MFEM_USE_MPI
         MPI_Abort(MPI_COMM_WORLD, 1);
#else
         std::abort();
#endif
      }
      else if (rank == 0)
      {
         std::cerr
            << "[WARNING] SEAS_FORCE_MIXED_FLUX_SUBSTEP_MPI=1 — running "
            "--fault-iterator substep + --mixed-flux " << mixed_flux_str
            << " at np=" << nprocs << " is UNTESTED at any rank count.  "
            "Numerical correctness is NOT guaranteed.  Cite R-1503/R-1605 "
            "if a regression is observed.\n";
      }
   }

   // ader-order is accepted verbatim; wave.AdvanceADER clamps/validates
   // internally.  Allowing it through avoids false warnings when the
   // smoke test exercises --ader-order 5 for banner-text verification.
   if (ader_order < 1) { ader_order = 2; }

   if (rank == 0)
   {
      std::cout << "========================================\n";
      std::cout << "SCEC TPV102 Dynamic Rupture Simulation\n";
      std::cout << "========================================\n";
      std::cout << "Ranks: " << nprocs << "\n";
      std::cout << "Mesh: " << mesh_file << "\n";
      std::cout << "Scale: " << mesh_scale << "\n";
      std::cout << "Order: " << order << "\n";
      std::cout << "BC: " << bc_mode << "\n";
      std::cout << "Tfinal: " << tfinal << " s\n";
      std::cout << "Output: " << output_dir << "/" << output_prefix << "\n";
      std::cout << "CFL: " << cfl_factor << "\n";
      std::cout << "BC attrs: free=" << bc_free
                << ", fault=" << bc_fault
                << ", absorb=" << bc_absorb << "\n";
      // R7-001 option (b): the banner describes what the driver
      // ACTUALLY runs, not what the CLI requested.  The CLI values
      // (friction_solver, fault_iterator) are echoed on a separate
      // "CLI parsed (banner-only)" line so log parsers can still see
      // them, but the four load-bearing lines are the dispatch truth
      // derived from GetDispatched* / BannerOf (R8-001).
      const DispatchedSolver   actual_solver = GetDispatchedSolver(friction_solver);
      const DispatchedIterator actual_iter   = GetDispatchedIterator(fault_iterator);
      const DispatchedLaw      actual_law    = GetDispatchedLaw(fric_law);

      std::cout << "Time integrator: ADER-O" << ader_order
                << " (one-shot via wave.AdvanceADER)\n";
      std::cout << "Fault iterator: " << BannerOf(actual_iter) << "\n";
      std::cout << "Friction solver: " << BannerOf(actual_solver) << "\n";
      std::cout << "Friction law: " << BannerOf(actual_law) << "\n";
      // Round-11 Mixed-Flux banner (Zhang et al. 2023).
      const char *mixed_flux_banner =
         (mixed_flux_mode == MixedFluxMode::None)
            ? "none (upwind everywhere, default)"
       : (mixed_flux_mode == MixedFluxMode::Adjacent)
            ? "adjacent (Mixed-Flux 2 per Zhang et al. 2023, central on "
              "fault-adjacent non-fault interior faces)"
       : (mixed_flux_mode == MixedFluxMode::AllContinuous)
            ? "all-continuous (Mixed-Flux 1, central on every interior "
              "non-fault face)"
            : "?";
      std::cout << "Mixed flux: " << mixed_flux_banner << "\n";
      std::cout << "Nucleation: "
                << (disable_nucleation ? "DISABLED"
                                       : "enabled (TPV102, macro-step "
                                         "incremental telescoping)")
                << "\n";
      std::cout << "CLI parsed (banner-only, not dispatched): "
                << "friction_solver=" << friction_solver
                << ", fault_iterator=" << fault_iterator
                << ", fric_law=" << fric_law
                << ", mixed_flux=" << mixed_flux_str << "\n";
      std::cout << "========================================\n\n";
      mkdir(output_dir.c_str(), 0755);
   }

   // R8-004: per-rank machine-readable tri-consistency lines.  Emit
   // OUTSIDE the rank==0 guard so multi-rank Frontera logs can self-
   // verify that every rank dispatches the same solver.  Each line is
   // prefixed with `rank=<N>` so interleaved output is still parseable.
   if (verify_dispatch)
   {
      const DispatchedSolver   actual_solver = GetDispatchedSolver(friction_solver);
      const DispatchedIterator actual_iter   = GetDispatchedIterator(fault_iterator);
      const DispatchedLaw      actual_law    = GetDispatchedLaw(fric_law);
      std::cout << "[dispatch] rank=" << rank
                << " friction_solver_actual=" << TagOf(actual_solver) << "\n";
      std::cout << "[dispatch] rank=" << rank
                << " fault_iterator_actual=" << TagOf(actual_iter) << "\n";
      std::cout << "[dispatch] rank=" << rank
                << " friction_law_actual=" << TagOf(actual_law) << "\n";
   }

#ifdef MFEM_USE_MPI
   MPI_Barrier(comm);
#endif

   // R7-005: validate --friction-solver eagerly so typos abort before
   // any simulation work (including --dry-run).
   //
   // R8-002: preserve the Method value in a named local (`method`)
   // rather than discarding MapSolver's return.  Under R7-001 option
   // (b) the value is not routed through the solver dispatch, so it
   // is cast to void here.  When R7-001 option (a) lands, the
   // required edit is a single site: remove `(void)method;` and pass
   // `method` into `Tpv102SubStepIterator::Advance`.  Keeping the
   // named local ensures grep / IDE reference finds the linkage point.
   // R-602/R-603 (round-7): `method` is now passed into
   // AdvanceADERWithSubStep on the substep dispatch path.  On the
   // legacy one-shot path it remains unused (Brent is hard-coded inside
   // EvaluateADER); the previous (void)method silencer is removed.
   const FrictionSolver::Method method = MapSolver(friction_solver);

   // --dry-run: no mesh, no simulation.  Print a canonical end-of-run
   // line that test_tpv102_smoke.cpp scrapes ("[dry-run] OK.").  Used
   // as the cheapest possible Frontera startup verification.
   if (dry_run)
   {
      if (rank == 0)
      {
         std::cout << "[tpv102_driver] --dry-run: banner printed, "
                   << "no mesh / no simulation.\n";
         std::cout << "[dry-run] OK.\n";
      }
#ifdef MFEM_USE_MPI
      MPI_Finalize();
#endif
      return 0;
   }

   // -----------------------------------------------------------------------
   // 1. Load serial mesh, partition to ParMesh
   // -----------------------------------------------------------------------
   Mesh serial_mesh(mesh_file.c_str(), 1, 1);
   MFEM_VERIFY(serial_mesh.Dimension() == 3, "TPV102 requires 3D mesh");
   if (mesh_scale != 1.0)
   {
      serial_mesh.SetCurvature(1, false, 3, Ordering::byVDIM);
      Vector &nodes = *serial_mesh.GetNodes();
      nodes *= mesh_scale;
   }

#ifdef MFEM_USE_MPI
   // G1 fault-locality partition (TPV102 dynamic only).  See
   // dynamic/fault_locality_partition.hpp for the algorithm and the
   // 2026-04-25_pm debug doc Section 13 for the motivation: ParMETIS
   // splits the mesh across the fault plane at np ≥ 4, producing a
   // partition-induced mirror-symmetry break that the rupture amplifies
   // to mm-scale slip_dip pollution.  G1 forces every fault face's two
   // adjacent elements to be co-resident on the same rank — generalizes
   // beyond y=0 mirror to arbitrary fault geometries.
   //
   // BP5 path NOT touched per user directive 2026-04-25.
   //
   // Lifetime note: Array<int> custom_partitioning is hoisted OUTSIDE
   // the if-block to guarantee its data outlives the ParMesh ctor, in
   // case MFEM stores the pointer rather than copying.  Job 7677864
   // showed bit-exact baseline output despite --partition-file being
   // passed, suggesting the partition was being silently dropped.
   std::unique_ptr<ParMesh> pmesh_ptr;
   Array<int> custom_partitioning;
   if (!partition_file.empty())
   {
      // Load explicit partitioning array from a sidecar file (typically
      // produced by tpv102/mesh/build_symmirror_mesh.py --emit-partition).
      const bool ok = seas::LoadPartitioningFromFile(partition_file, nprocs,
                                                     serial_mesh.GetNE(),
                                                     custom_partitioning);
      MFEM_VERIFY(ok, "Failed to load partition file: " << partition_file
                  << " (np_expected=" << nprocs
                  << ", ne_expected=" << serial_mesh.GetNE() << ")");
      if (rank == 0)
      {
         // Verify partition is actually distributed across all ranks (not
         // accidentally all-zeros or all-one-rank).  Compute count of
         // unique ranks and per-rank element count from the loaded array.
         std::vector<int> rank_count(nprocs, 0);
         for (int e = 0; e < custom_partitioning.Size(); e++)
         {
            const int r = custom_partitioning[e];
            if (r >= 0 && r < nprocs) { rank_count[r]++; }
         }
         int n_used = 0, min_e = INT_MAX, max_e = 0;
         for (int r = 0; r < nprocs; r++)
         {
            if (rank_count[r] > 0)
            {
               n_used++;
               if (rank_count[r] < min_e) { min_e = rank_count[r]; }
               if (rank_count[r] > max_e) { max_e = rank_count[r]; }
            }
         }
         std::cout << "[partition] loaded from " << partition_file
                   << " (np=" << nprocs << ", ne=" << custom_partitioning.Size()
                   << ", ranks_used=" << n_used << "/" << nprocs
                   << ", elems_per_rank=" << min_e << ".." << max_e
                   << ", first_5=[" << custom_partitioning[0] << ","
                   << custom_partitioning[1] << ","
                   << custom_partitioning[2] << ","
                   << custom_partitioning[3] << ","
                   << custom_partitioning[4] << "])" << std::endl;
      }
      pmesh_ptr.reset(new ParMesh(comm, serial_mesh,
                                  custom_partitioning.GetData()));
      // Verify ParMesh actually used our partition: each rank's local
      // element count must equal rank_count[my_rank].  ABORT if not —
      // job 7677864 silently produced bit-exact baseline output despite
      // --partition-file being passed; we want loud failure if MFEM
      // ignores the partition.
      const int local_ne = pmesh_ptr->GetNE();
      int expected_ne = 0;
      for (int e = 0; e < custom_partitioning.Size(); e++)
      {
         if (custom_partitioning[e] == rank) { expected_ne++; }
      }
      int local_honored = (local_ne == expected_ne) ? 1 : 0;
      int all_honored = 0;
      MPI_Allreduce(&local_honored, &all_honored, 1, MPI_INT, MPI_MIN, comm);
      int sum_local = 0, sum_expected = 0;
      MPI_Allreduce(&local_ne, &sum_local, 1, MPI_INT, MPI_SUM, comm);
      MPI_Allreduce(&expected_ne, &sum_expected, 1, MPI_INT, MPI_SUM, comm);
      if (rank == 0)
      {
         std::cout << "[partition] ParMesh local NE: rank0_actual="
                   << local_ne << " rank0_expected=" << expected_ne
                   << "  global_actual=" << sum_local
                   << " global_expected=" << sum_expected
                   << "  all_honored=" << (all_honored ? "YES" : "NO")
                   << std::endl;
      }
      if (!all_honored)
      {
         if (rank == 0)
         {
            std::cerr << "[partition] FATAL: at least one rank's local NE "
                      << "does not match the partition file.  MFEM is not "
                      << "applying our custom partitioning array.  Aborting "
                      << "to avoid silent fallback to ParMETIS-default."
                      << std::endl;
         }
         MPI_Abort(comm, 73);
      }
   }
   else if (fault_locality_part)
   {
      // Identify fault faces in the serial mesh by bdr_attr == bc_fault.
      // bc_fault default is 3 per TPV102 driver; same convention used by
      // wave_operator.inl when populating fault_interior_faces_.
      const int bc_fault_attr = GetIntArg(argc, argv, "--bc-fault", 3);
      Array<int> fault_faces;
      seas::FindFaultFaceIndices(serial_mesh, bc_fault_attr, fault_faces);
      Array<int> partitioning;
      int n_relocated = 0;
      seas::BuildFaultLocalityPartitioning(serial_mesh, fault_faces, nprocs,
                                            partitioning, &n_relocated);
      const int violations = seas::VerifyFaultLocality(serial_mesh,
                                                       fault_faces,
                                                       partitioning);
      if (rank == 0)
      {
         std::cout << "[partition] fault-locality partitioning enabled: "
                   << fault_faces.Size() << " fault faces, "
                   << n_relocated << " elements relocated, "
                   << violations << " violations" << std::endl;
      }
      MFEM_VERIFY(violations == 0,
                  "Fault-locality partition has " << violations
                  << " violations — partitioning logic bug.");
      pmesh_ptr.reset(new ParMesh(comm, serial_mesh, partitioning.GetData()));
   }
   else
   {
      pmesh_ptr.reset(new ParMesh(comm, serial_mesh));
   }
   ParMesh &pmesh = *pmesh_ptr;
   using MeshT = ParMesh;
#else
   Mesh &pmesh = serial_mesh;
   using MeshT = Mesh;
#endif

   int ne_local = pmesh.GetNE();
   int ne_global = ne_local;
#ifdef MFEM_USE_MPI
   MPI_Allreduce(&ne_local, &ne_global, 1, MPI_INT, MPI_SUM, comm);
#endif
   if (rank == 0)
   {
      std::cout << "Mesh: " << ne_global << " elements, "
                << nprocs << " ranks\n";
   }

   // H_α diagnostic: dump per-rank local-element-order vs physical
   // y-coordinate so we can detect when MFEM's re-ordering breaks the
   // y-mirror-pairing within a single rank's local domain.  Gated by
   // env var SEAS_DIAG_PARMESH_ORDER=1.
   if (std::getenv("SEAS_DIAG_PARMESH_ORDER") != nullptr)
   {
      // Each rank writes its local element centroids to a sidecar file:
      //   /tmp/parmesh_order_rank<R>.txt with lines:
      //     local_id  cx  cy  cz
      // Then offline we compare across ranks for the y-mirror property.
      char fname[256];
      std::snprintf(fname, sizeof(fname),
                    "/tmp/parmesh_order_rank%d_np%d.txt", rank, nprocs);
      std::ofstream ofs(fname);
      ofs.precision(15);
      for (int e = 0; e < pmesh.GetNE(); e++)
      {
         ElementTransformation *Tr = pmesh.GetElementTransformation(e);
         IntegrationPoint ip;
         ip.x = ip.y = ip.z = 0.25;   // tet barycentre in reference
         Vector phys(3);
         Tr->Transform(ip, phys);
         ofs << e << " " << phys(0) << " " << phys(1)
             << " " << phys(2) << "\n";
      }
      ofs.close();
      if (rank == 0)
      {
         std::cout << "[diag-order] Wrote per-rank local element "
                   << "ordering to /tmp/parmesh_order_rank*_np"
                   << nprocs << ".txt" << std::endl;
      }
   }

   // -----------------------------------------------------------------------
   // 2. Boundary conditions
   // -----------------------------------------------------------------------
   BoundaryConfig bc;
   bc.natural_attrs   = {bc_free};
   bc.fault_attr      = bc_fault;
   bc.absorbing_attrs = {bc_absorb};

   {
      int n_free = 0, n_fault = 0, n_absorb = 0;
      for (int b = 0; b < pmesh.GetNBE(); b++)
      {
         int attr = pmesh.GetBdrAttribute(b);
         if (bc.natural_attrs.count(attr))   { ++n_free; }
         if (attr == bc.fault_attr)          { ++n_fault; }
         if (bc.absorbing_attrs.count(attr)) { ++n_absorb; }
      }
      int n_free_g = n_free, n_fault_g = n_fault, n_absorb_g = n_absorb;
#ifdef MFEM_USE_MPI
      MPI_Allreduce(&n_free,   &n_free_g,   1, MPI_INT, MPI_SUM, comm);
      MPI_Allreduce(&n_fault,  &n_fault_g,  1, MPI_INT, MPI_SUM, comm);
      MPI_Allreduce(&n_absorb, &n_absorb_g, 1, MPI_INT, MPI_SUM, comm);
#endif
      if (rank == 0)
      {
         std::cout << "BC faces — free: " << n_free_g
                   << ", fault: " << n_fault_g
                   << ", absorb: " << n_absorb_g << "\n";
      }
      MFEM_VERIFY(n_fault_g > 0,
                  "No fault faces with attr=" << bc.fault_attr
                  << " found.  Check --bc-fault or mesh Physical Surface tags.");
   }

   // -----------------------------------------------------------------------
   // 3. WaveOperator
   // -----------------------------------------------------------------------
   WaveOperator<MeshT> wave(pmesh, order,
                            TPV102Params::lambda, TPV102Params::mu,
                            TPV102Params::rho, bc);

   int ndof_total  = wave.GetScalarNDof();
   int ndof_global = ndof_total;
#ifdef MFEM_USE_MPI
   MPI_Allreduce(&ndof_total, &ndof_global, 1, MPI_INT, MPI_SUM, comm);
#endif
   if (rank == 0)
   {
      std::cout << "DOFs per component (global): " << ndof_global << "\n";
      std::cout << "Total DOFs (global): " << NUM_STATE * ndof_global << "\n";
   }

   // -----------------------------------------------------------------------
   // 4. Time step
   // -----------------------------------------------------------------------
   real_t cfl    = cfl_factor / (3.0 * (2.0 * order + 1.0));
   real_t dt_cfl = wave.ComputeMaxDt(cfl);
   real_t dt     = (dt_override > 0.0) ? dt_override : dt_cfl;
   int nsteps    = (tfinal > 0.0) ? static_cast<int>(std::ceil(tfinal / dt)) : 0;
   if (rank == 0)
   {
      std::cout << "CFL: " << cfl << ", dt_cfl = " << dt_cfl << " s\n";
      std::cout << "dt: " << dt << " s, steps: " << nsteps << "\n\n";
   }

   // -----------------------------------------------------------------------
   // 5. Fault DOF data (local partition) — TPV102-specific init
   // -----------------------------------------------------------------------
   L2_FECollection fec(order, 3, BasisType::GaussLobatto);
#ifdef MFEM_USE_MPI
   ParFiniteElementSpace fes(&pmesh, &fec);
#else
   FiniteElementSpace fes(&pmesh, &fec);
#endif

   const Array<int> &fault_int_faces = wave.GetFaultInteriorFaces();
   const Array<int> &fault_shr_faces = wave.GetFaultSharedFaces();

   int nqp_per_face = 0;
   if (fault_int_faces.Size() > 0)
   {
      FaceElementTransformations *ftr0 =
         pmesh.GetInteriorFaceTransformations(fault_int_faces[0]);
      MFEM_VERIFY(ftr0, "fault interior face has null transformation");
      nqp_per_face = IntRules.Get(ftr0->GetGeometryType(), 2 * order).GetNPoints();
   }
#ifdef MFEM_USE_MPI
   else if (fault_shr_faces.Size() > 0)
   {
      FaceElementTransformations *ftr0 =
         pmesh.GetSharedFaceTransformations(fault_shr_faces[0]);
      MFEM_VERIFY(ftr0, "fault shared face has null transformation");
      nqp_per_face = IntRules.Get(ftr0->GetGeometryType(), 2 * order).GetNPoints();
   }
   {
      int local_nqp = nqp_per_face;
      MPI_Allreduce(&local_nqp, &nqp_per_face, 1, MPI_INT, MPI_MAX, comm);
   }
#endif

   const int num_fault_local  = fault_int_faces.Size() * nqp_per_face;
   const int num_shared_fault = fault_shr_faces.Size() * nqp_per_face;
   const int num_fault_total  = num_fault_local + num_shared_fault;

   std::vector<Vector> fault_coords;
   fault_coords.reserve(num_fault_total);

   auto push_qps = [&](FaceElementTransformations *ftr)
   {
      const IntegrationRule &ir = IntRules.Get(
         ftr->GetGeometryType(), 2 * order);
      for (int q = 0; q < ir.GetNPoints(); ++q)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3);
         ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   };
   for (int i = 0; i < fault_int_faces.Size(); ++i)
   {
      push_qps(pmesh.GetInteriorFaceTransformations(fault_int_faces[i]));
   }
#ifdef MFEM_USE_MPI
   for (int i = 0; i < fault_shr_faces.Size(); ++i)
   {
      push_qps(pmesh.GetSharedFaceTransformations(fault_shr_faces[i]));
   }
#endif

   int num_fault_global = num_fault_total;
#ifdef MFEM_USE_MPI
   MPI_Allreduce(&num_fault_total, &num_fault_global, 1, MPI_INT, MPI_SUM, comm);
#endif
   if (rank == 0)
   {
      std::cout << "Fault QPs (global): " << num_fault_global
                << " (local: " << num_fault_local
                << ", shared: " << num_shared_fault << ")\n";
   }

#ifdef SEAS_DIAG_TPV102_FAULT_BASIS
   // D1 instrumentation companion: dump fault QP coordinates indexed
   // by dof_idx (= the same index keyed by [diag-flip] in
   // wave_operator.inl).  Join via:
   //   awk '/diag-coords/{c[$2]=$0} /diag-flip/{f[$2]=$0} END{for(k in f)print f[k]" "c[k]}'
   for (int i = 0; i < num_fault_total; ++i)
   {
      const Vector &xqp = fault_coords[i];
      std::fprintf(stderr,
         "[diag-coords] dof=%d (x,y,z)=(%+9.1f,%+9.1f,%+9.1f) rank=%d\n",
         i, xqp(0), xqp(1), xqp(2), rank);
   }
#endif

   std::vector<DOFData> dof_data;
   if (num_fault_total > 0)
   {
      // TPV102 fluctuation-Q init: Q = 0, pre-stress lives in DOFData.
      InitializeFaultDOFs(dof_data, num_fault_total, fault_coords);
   }

   // Resolve the rank that owns the hypocenter QP (closest local fault QP
   // to (hypo_along_strike, -hypo_down_dip) in x/z), then tag that DOF
   // with diag_print = true so the C-1 EVAL / C-1n NORMAL probes in
   // dynamic/fault_face_flux.cpp emit one line per Evaluate call instead
   // of flooding stderr from every QP.  Mirrors the TPV102 pattern at
   // tpv102_driver.cpp:540-588.
   int hypo_rank = 0;
   int hypo_dof_local = -1;
   {
      real_t local_min_dist2 = std::numeric_limits<real_t>::max();
      for (int i = 0; i < num_fault_total; i++)
      {
         real_t dx = fault_coords[i](0) - TPV102Params::hypo_along_strike;
         real_t dz = std::abs(fault_coords[i](2)) - TPV102Params::hypo_down_dip;
         real_t d2 = dx*dx + dz*dz;
         if (d2 < local_min_dist2)
         {
            local_min_dist2 = d2;
            hypo_dof_local = i;
         }
      }
#ifdef MFEM_USE_MPI
      struct MinDist { double d; int r; };
      static_assert(offsetof(MinDist, d) == 0,
                    "MinDist.d must be at offset 0 for MPI_DOUBLE_INT");
      static_assert(offsetof(MinDist, r) == sizeof(double),
                    "MinDist.r must follow d with no padding for MPI_DOUBLE_INT");
      MinDist in{static_cast<double>(local_min_dist2), rank};
      MinDist out{};
      MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, comm);
      hypo_rank = out.r;
#endif
   }

#ifdef SEAS_DIAG_FAULT_FLUX
   // Only the rank that won the MINLOC sets diag_print = true on its
   // closest hypocenter DOF; at most one DOF is flagged globally.
   // Probe block: dynamic/fault_face_flux.cpp:259-303 (C-1 EVAL +
   // C-1n NORMAL).  Single-rank stderr printf, no MPI calls — no
   // deadlock at any rank count.
   //
   // C-2 BULK-PROBE EXTENSION: also identify the fault face containing
   // the diag DOF and the two adjacent tets, plus the 6 non-fault
   // interior faces of those tets.  Pushed to WaveOperator via setters
   // so the C-2A (Mult per-call) and C-2B (per non-fault face) probes
   // can fire without MPI on this single rank.  Other ranks have the
   // setters left at default -1 → silent.
   int diag_face_idx_for_dof = -1;
   int diag_elem_plus  = -1, diag_elem_minus = -1;
   int diag_face_dof_plus = -1, diag_face_dof_minus = -1;
   std::vector<int> diag_nonfault_faces;
   if (rank == hypo_rank && hypo_dof_local >= 0 && num_fault_total > 0)
   {
      dof_data[hypo_dof_local].diag_print = true;
      const Vector &hpos = fault_coords[hypo_dof_local];
      std::fprintf(stderr,
         "[diag] rank %d tagging hypo DOF %d at (%.1f, %.1f, %.1f)\n",
         rank, hypo_dof_local, hpos(0), hpos(1), hpos(2));

      // Map hypo_dof_local back to (interior fault face index, q-of-face).
      // Layout: dof_data[i] lives on interior face i / nqp_per_face,
      // QP index = i % nqp_per_face.  Shared-fault DOFs follow at
      // [num_fault_local, num_fault_total) — currently not used as
      // diag target (the hypo nucleus lives on interior faces).
      if (hypo_dof_local < num_fault_local)
      {
         const int fi = hypo_dof_local / nqp_per_face;     // index into fault_int_faces
         const int qi = hypo_dof_local % nqp_per_face;
         if (fi >= 0 && fi < fault_int_faces.Size())
         {
            diag_face_idx_for_dof = fault_int_faces[fi];
            FaceElementTransformations *ftr =
               pmesh.GetInteriorFaceTransformations(diag_face_idx_for_dof);
            if (ftr)
            {
               const int e1 = ftr->Elem1No;
               const int e2 = ftr->Elem2No;
               // Determine canonical + / - via FaultBasis sign_flipped at
               // this QP (matches wave_operator.inl:1226 logic:
               // elem1_on_plus = !sign_flipped).
               bool elem1_on_plus = true;
               const FaultBasis *fb = wave.GetFaultBasis();
               if (fb)
               {
                  const int fb_idx =
                     wave.LookupInteriorFaultBasisIndex(diag_face_idx_for_dof);
                  if (fb_idx >= 0 && fb_idx < fb->NumFaces())
                  {
                     const FaultBasisData &bd = fb->GetBasis(fb_idx);
                     if (qi < static_cast<int>(bd.qp_data.size()))
                     {
                        elem1_on_plus = !bd.qp_data[qi].sign_flipped;
                     }
                  }
               }
               diag_elem_plus  = elem1_on_plus ? e1 : e2;
               diag_elem_minus = elem1_on_plus ? e2 : e1;

               // Find the local DOF index on each tet that is closest
               // to the hypocenter QP physical position AND lies ON
               // THE FAULT FACE (|y - hpos.y| < tol).  Without the
               // on-face restriction, the apex DOF (y ≈ ±141 m for a
               // 200 m mesh) can win the closest-by-Euclidean-distance
               // contest when hpos is near a fault-face edge, putting
               // C-2A at qualitatively different positions on E+ vs
               // E- and producing a spurious asymmetry signature.
               // Returns the on-face DOF closest to hpos in (x,z); if
               // no DOF is within tol of hpos.y, falls back to the
               // closest by full distance and reports the y offset.
               auto closest_face_dof_idx = [&](int e,
                                               real_t &out_y_off) -> int
               {
                  const FiniteElement *fe = wave.GetFESpace().GetFE(e);
                  ElementTransformation *Tr =
                     wave.GetFESpace().GetElementTransformation(e);
                  const IntegrationRule &nodes = fe->GetNodes();
                  const real_t y_tol = 1e-3;  // 1 mm — much smaller than 200 m mesh
                  int best_on_face = -1;
                  real_t best_d2_on_face =
                     std::numeric_limits<real_t>::max();
                  int best_any = -1;
                  real_t best_d2_any =
                     std::numeric_limits<real_t>::max();
                  real_t best_y_any = 0.0;
                  for (int k = 0; k < nodes.GetNPoints(); k++)
                  {
                     Vector phys(3);
                     Tr->Transform(nodes.IntPoint(k), phys);
                     const real_t dx = phys(0) - hpos(0);
                     const real_t dy = phys(1) - hpos(1);
                     const real_t dz = phys(2) - hpos(2);
                     const real_t d2_full = dx*dx + dy*dy + dz*dz;
                     const real_t d2_xz   = dx*dx + dz*dz;
                     if (std::abs(dy) < y_tol && d2_xz < best_d2_on_face)
                     {
                        best_d2_on_face = d2_xz; best_on_face = k;
                     }
                     if (d2_full < best_d2_any)
                     {
                        best_d2_any = d2_full; best_any = k;
                        best_y_any = phys(1);
                     }
                  }
                  if (best_on_face >= 0)
                  {
                     out_y_off = 0.0;
                     return best_on_face;
                  }
                  out_y_off = best_y_any - hpos(1);
                  return best_any;
               };
               real_t y_off_p = 0.0, y_off_m = 0.0;
               diag_face_dof_plus  =
                  closest_face_dof_idx(diag_elem_plus,  y_off_p);
               diag_face_dof_minus =
                  closest_face_dof_idx(diag_elem_minus, y_off_m);

               // Print physical coordinates of the two diag DOFs so
               // the C-2A/B/C analysis can verify they are at mirror
               // positions before drawing conclusions about bulk
               // asymmetry.
               auto dof_phys = [&](int e, int k) -> Vector
               {
                  const FiniteElement *fe = wave.GetFESpace().GetFE(e);
                  ElementTransformation *Tr =
                     wave.GetFESpace().GetElementTransformation(e);
                  Vector phys(3);
                  Tr->Transform(fe->GetNodes().IntPoint(k), phys);
                  return phys;
               };
               Vector p_plus  = dof_phys(diag_elem_plus,  diag_face_dof_plus);
               Vector p_minus = dof_phys(diag_elem_minus, diag_face_dof_minus);
               std::fprintf(stderr,
                  "[diag-c2-pos] rank=%d  DOF+ at (%+.3e,%+.3e,%+.3e) "
                  "y_off=%+.3e   DOF- at (%+.3e,%+.3e,%+.3e) y_off=%+.3e   "
                  "dx=%+.3e dz=%+.3e\n",
                  rank,
                  p_plus(0),  p_plus(1),  p_plus(2),  y_off_p,
                  p_minus(0), p_minus(1), p_minus(2), y_off_m,
                  p_plus(0) - p_minus(0), p_plus(2) - p_minus(2));

               // Collect the non-fault interior faces of the two diag
               // tets.  For tets, GetElementFaces returns 4 faces; we
               // skip the fault face itself and any boundary faces
               // (which are still printed elsewhere).
               Array<int> faces_p, ori_p, faces_m, ori_m;
               pmesh.GetElementFaces(diag_elem_plus,  faces_p, ori_p);
               pmesh.GetElementFaces(diag_elem_minus, faces_m, ori_m);
               auto add_nonfault = [&](const Array<int> &fs)
               {
                  for (int k = 0; k < fs.Size(); k++)
                  {
                     if (fs[k] == diag_face_idx_for_dof) { continue; }
                     diag_nonfault_faces.push_back(fs[k]);
                  }
               };
               add_nonfault(faces_p);
               add_nonfault(faces_m);

               std::fprintf(stderr,
                  "[diag-c2] rank %d face=%d e+=%d e-=%d "
                  "face_dof+=%d face_dof-=%d nonfault_faces=[",
                  rank, diag_face_idx_for_dof,
                  diag_elem_plus, diag_elem_minus,
                  diag_face_dof_plus, diag_face_dof_minus);
               for (size_t k = 0; k < diag_nonfault_faces.size(); k++)
               {
                  std::fprintf(stderr, "%s%d",
                               k == 0 ? "" : ",", diag_nonfault_faces[k]);
               }
               std::fprintf(stderr, "]\n");
            }
         }
      }
   }
   wave.SetDiagBulkElems(diag_elem_plus, diag_elem_minus);
   wave.SetDiagBulkFaceDofs(diag_face_dof_plus, diag_face_dof_minus);
   wave.SetDiagNonFaultFaces(diag_nonfault_faces);
#endif

   FaultFaceFlux fault_flux(TPV102Params::rho, TPV102Params::cp,
                            TPV102Params::cs);
   wave.SetFaultFlux(&fault_flux);
   wave.SetFaultDOFData(&dof_data, nqp_per_face);

   // Zero Q_bg — fluctuation-Q dispatch.
   {
      real_t Q_bg[NUM_STATE] = {0};
      wave.SetAbsorbingBackground(Q_bg);
   }

   // Round-11 Mixed-Flux dispatch wiring (R-1205 required call order):
   //   1. WaveOperator ctor (already done)
   //   2. wave.SetFaultFlux         (already done)
   //   3. wave.SetFaultDOFData      (already done)
   //   4. wave.SetAbsorbingBackground (already done)
   //   5. wave.SetMixedFluxMode     <-- HERE; setter cross-checks
   //                                    bc_.fault_attr > 0 for Adjacent,
   //                                    and aborts if precomputed-flux
   //                                    is also enabled (R-1203).
   wave.SetMixedFluxMode(mixed_flux_mode);
   // R-1202: report the GLOBAL central-set size and per-rank min/max
   // for load-balance diagnostic.  The Round-2 banner reported rank-0's
   // local size as if it were global, masking 127× of the count at
   // np=128 production.  All ranks participate in the reduction; only
   // rank 0 prints.
   if (mixed_flux_mode != MixedFluxMode::None)
   {
      const long long local_size_ll =
         static_cast<long long>(wave.GetCentralFluxFaceSet().size());
#ifdef MFEM_USE_MPI
      // R-1409: Allreduce so EVERY rank can self-verify its own log
      // (per-rank log-grep parity with R8-004's tri-consistency lines).
      // Cost is negligible (3 × long_long per call, called once per run).
      // R-1602: use the WaveOperator's ParMesh communicator
      // (`pmesh.GetComm()`), not `MPI_COMM_WORLD`, for symmetry with
      // SetMixedFluxMode's consensus check (which uses `pmesh.GetComm()`
      // at wave_operator.inl:1199-1213).  Currently both communicators
      // coincide because TPV102 builds its ParMesh on MPI_COMM_WORLD,
      // but a future multi-region driver running on a sub-comm would
      // deadlock or diverge if the two collectives use mismatched
      // communicators.
      MPI_Comm wave_comm = pmesh.GetComm();
      long long global_sum = 0, local_max = 0, local_min = 0;
      MPI_Allreduce(&local_size_ll, &global_sum, 1, MPI_LONG_LONG, MPI_SUM,
                    wave_comm);
      MPI_Allreduce(&local_size_ll, &local_max, 1, MPI_LONG_LONG, MPI_MAX,
                    wave_comm);
      MPI_Allreduce(&local_size_ll, &local_min, 1, MPI_LONG_LONG, MPI_MIN,
                    wave_comm);
      if (rank == 0)
      {
         std::cout << "[mixed-flux] mode=" << mixed_flux_str
                   << "  |central_set|_global=" << global_sum
                   << "  per-rank min=" << local_min
                   << " max=" << local_max
                   << "  (Zhang et al. 2023 mixed-flux dispatch)\n";
         if (local_min > 0 && local_max > 4 * local_min)
         {
            std::cout << "[mixed-flux] WARNING: rank load imbalance "
                      << static_cast<double>(local_max) /
                         static_cast<double>(local_min)
                      << "× (max/min); consider "
                      << "--partition-fault-locality\n";
         }
         else if (local_min == 0 && local_max > 0)
         {
            std::cout << "[mixed-flux] WARNING: at least one rank has "
                      << "zero central-set entries while another has "
                      << local_max << "; partition is severely "
                      << "fault-asymmetric.\n";
         }
      }
      // R-1409: machine-readable per-rank line.  Every rank prints once;
      // log parsers can grep "[mixed-flux] rank=N" to verify dispatch
      // engagement and global counts visible from any single log file.
      std::cout << "[mixed-flux] rank=" << rank
                << "  local_set_size=" << local_size_ll
                << "  global_sum=" << global_sum
                << "  global_min=" << local_min
                << "  global_max=" << local_max << "\n";
#else
      std::cout << "[mixed-flux] mode=" << mixed_flux_str
                << "  |central_set|=" << local_size_ll
                << " (serial)  (Zhang et al. 2023 mixed-flux dispatch)\n";
#endif
   }

   // TPV102 uses the regularised rate-and-state ageing law (SCEC TPV101/102
   // §"Friction Law" Eq. (2)).  AgingLawPsi from friction/state_evolution.hpp
   // implements dψ/dt in ψ-space; the iterator queries (b, V0, f0) accessors
   // for the analytic ψ update inside the per-sub-step loop.
   mfem::seas::AgingLawPsi state_evo(
      TPV102Params::b, TPV102Params::V0, TPV102Params::f0);

   // R-602/R-603 substep iterator (default OFF; opt-in via
   // --fault-iterator substep).  SetSubSteps configures the ADER-O
   // quadrature: at order O, equal-width sub-steps with equal weights
   // 1/O is the simplest valid quadrature on [0, dt] satisfying
   //   Σ deltaT[o] == dt_macro,  Σ time_weights[o] == 1.
   // The iterator's own per-call argument validation enforces this.
   mfem::seas::Tpv102SubStepIterator substep_iterator(fault_flux, state_evo);
   {
      const int O = std::max(1, ader_order);
      std::vector<real_t> deltaT(O, 1.0 / static_cast<real_t>(O));
      std::vector<real_t> weights(O, 1.0 / static_cast<real_t>(O));
      // The iterator interprets deltaT in absolute (physical) time units,
      // so seed it with the auto-CFL `dt` here.  The Σ deltaT==dt_macro
      // check inside Advance/AdvanceWithSubStepStates is RELATIVE; the
      // helper AdvanceADERWithSubStep rescales the configured deltaT to
      // the actual dt_step at each call (see the dt_scale loop above),
      // so the final macro-step (where dt_step = tfinal - t < dt) is
      // handled correctly.
      for (int o = 0; o < O; o++) { deltaT[o] = dt / static_cast<real_t>(O); }
      substep_iterator.SetSubSteps(deltaT, weights);
   }

   const bool use_substep_iterator =
      (GetDispatchedIterator(fault_iterator) == DispatchedIterator::SubStep);
   if (rank == 0 && use_substep_iterator)
   {
      std::cout << "[tpv102_driver] --fault-iterator substep ACTIVE: "
                << "ADER-O" << ader_order
                << " per-sub-step Q via ComputeADERSubStepStates + "
                << "Tpv102SubStepIterator::AdvanceWithSubStepStates.\n";
   }
   // R-1003 (LANDED): the shared-fault ADER branch
   // (`ComputeADERSharedFaceFluxRHS`) now consults `substep_I_imp_*_flat_`
   // under the same absolute-index gate as the interior branch.  Driver
   // sizing (this file, `AdvanceADERWithSubStep` helper above) and
   // `WaveOperator::EvaluateBulkAtFaultQPsCanonical` (wave_operator.inl)
   // cover the full `[0, GetNumTotalFaultQPs())` index space, including
   // `[GetNumLocalFaultQPs(), GetNumTotalFaultQPs())` for shared-fault QPs.
   // The pre-R-1003 abort that rejected np>1 with --fault-iterator substep
   // is therefore retired.  See miniapps/seas/debug_document/
   // tpv102_debug_document/SUBSTEP_ITERATOR_MPI_REVIEW.md for the
   // pre-landing review and the merge-blocking MPI parity tests
   // (`test_tpv102_substep_iterator_mpi.cpp` and friends).
   if (use_substep_iterator && nprocs > 1 && rank == 0)
   {
      std::cout << "[tpv102_driver] R-1003 path: --fault-iterator substep "
                << "active under MPI (nprocs=" << nprocs << ").  Shared-"
                << "fault QPs flow through the substep side-channel.\n";
   }
   // -----------------------------------------------------------------------
   // 6. Initialize Q = 0 (fluctuation-Q).
   // -----------------------------------------------------------------------
   Vector Q(NUM_STATE * ndof_total);
   InitializeState(Q, ndof_total);

   MFEM_VERIFY(wave.GetAbsorbingBackground() != nullptr,
               "tpv102_driver: SetAbsorbingBackground(Q_bg=0) not called.");

   // R-403 PROBE: end-to-end y-mirror test of the wave operator.
   // Triggered by SEAS_DIAG_R403=1.  Bypasses the time loop entirely:
   //   - Q = 0 (already)
   //   - tau*_nuc = 0 (no nucleation injected; Q stays zero in bulk)
   //   - call wave.Mult(Q, dQdt) once
   //   - dump dQdt per-DOF with physical position to a text file
   //   - exit
   // Post-processing (separate Python) verifies that for every
   // (+y, -y) DOF mirror pair, the channel-signed antisymmetric-mirror
   // relation holds at FP-bit precision.  This is a STRICT SUPERSET of
   // probing CalcOrtho normals: a PASS rules out every operator-side
   // mirror leak (CalcOrtho, Loc1.Transform, CalcShape, face-flux
   // accumulation order, Elem1/Elem2 assignment, etc.).
   {
      const char *r403 = std::getenv("SEAS_DIAG_R403");
      if (r403 && r403[0] != '\0' &&
          !(r403[0] == '0' && r403[1] == '\0'))
      {
         Vector dQdt(NUM_STATE * ndof_total);
         dQdt = 0.0;
         wave.Mult(Q, dQdt);

         std::ostringstream pathss;
         pathss << output_dir << "/r403_rhs_rank" << rank << ".txt";
         std::ofstream fs(pathss.str());
         if (fs.is_open())
         {
            fs << std::scientific << std::setprecision(17);
            fs << "# columns: dof_global_idx x y z channel value\n";
            const int ne_local = pmesh.GetNE();
            const int ndof_per_el = wave.GetFESpace().GetFE(0)->GetDof();
            for (int e = 0; e < ne_local; e++)
            {
               const FiniteElement *fe = wave.GetFESpace().GetFE(e);
               ElementTransformation *Tr =
                  wave.GetFESpace().GetElementTransformation(e);
               const IntegrationRule &nodes = fe->GetNodes();
               for (int i = 0; i < ndof_per_el; i++)
               {
                  const IntegrationPoint &ip = nodes.IntPoint(i);
                  Tr->SetIntPoint(&ip);
                  Vector x(3);
                  Tr->Transform(ip, x);
                  const int dof_global = e * ndof_per_el + i;
                  for (int c = 0; c < NUM_STATE; c++)
                  {
                     fs << dof_global << " "
                        << x(0) << " " << x(1) << " " << x(2) << " "
                        << c << " "
                        << dQdt[c * ndof_total + dof_global] << "\n";
                  }
               }
            }
            fs.close();
         }
         if (rank == 0)
         {
            std::cout << "[R-403] dumped per-DOF rhs to "
                      << output_dir << "/r403_rhs_rank<R>.txt\n"
                      << "[R-403] exiting after probe.\n";
         }
#ifdef MFEM_USE_MPI
         MPI_Finalize();
#endif
         return 0;
      }
   }

   if (rank == 0)
   {
      std::cout << "Background: τ_strike = "
                << TPV102Params::tau_ini / 1e6 << " MPa, "
                << "σ_n = " << TPV102Params::sigma_n / 1e6 << " MPa, "
                << "θ_ini(a_vw=" << TPV102Params::a_vw
                << ") ≈ 1.606e9 s (= 50.9 yr, SCEC spec), V_ini = "
                << TPV102Params::V_ini << " m/s\n\n";
   }

   // -----------------------------------------------------------------------
   // 7. Station output
   // -----------------------------------------------------------------------
   auto stations = DefaultStations();
   TPV102StationWriter station_writer;
#ifdef MFEM_USE_MPI
   station_writer.Open(output_dir, output_prefix, stations,
                       fault_coords, num_fault_local, comm);
#else
   station_writer.Open(output_dir, output_prefix, stations,
                       fault_coords, num_fault_local);
#endif
   station_writer.WriteStep(0.0, dof_data);

   auto surface_stations = DefaultSurfaceStations();
   TPV102SurfaceStationWriter surface_writer;
   surface_writer.Open(output_dir, output_prefix, surface_stations,
                       pmesh, fes);
   surface_writer.WriteStep(0.0, Q);

   // -----------------------------------------------------------------------
   // 7b. ParaView output (mirrors tpv102_driver.cpp / BP5 seas::ParaViewOutput
   //     pattern).  Two collections:
   //       - pv_out (output_dir/ParaView): velocity + mpi_rank volume +
   //         fault-surface PVD/VTU (slip, slip_rate, traction dip+strike,
   //         psi, sigma_n, plus static a, Dc, x2, x3) at the fault schedule.
   //       - pv_bulk_out (output_dir/ParaView_bulk): velocity + sigma_yy +
   //         sigma_xy + sigma_xz + mpi_rank at --paraview-bulk-dt cadence.
   //
   // R-801 / BP5 component convention enforced project-wide: comp 0 = dip,
   // comp 1 = strike.  TPV102 is pure strike-slip so the strike channel
   // carries the rupture; dip stays near zero.
   //
   // Under R7-001 option (b) ADER one-shot, there is no "stage-4" snapshot
   // distinct from the time-averaged DOFData — write the same values into
   // both the averaged and the _k4 buffers so the ParaView Calculator
   // delta `<field>_k4 - <field>` reads zero everywhere (consistent with
   // what tpv102_driver.cpp does on the ADER path).
   // -----------------------------------------------------------------------
   using PvFES = typename seas::GFType<MeshT>::FESType;
   using PvGF  = typename seas::GFType<MeshT>::type;

   std::unique_ptr<seas::ParaViewOutput<MeshT>> pv_out;
   std::unique_ptr<L2_FECollection> pv_vel_fec, pv_rank_fec;
   std::unique_ptr<PvFES> pv_vel_fes, pv_rank_fes;
   std::unique_ptr<PvGF>  pv_vel_gf,  pv_rank_gf;

   std::unique_ptr<seas::ParaViewOutput<MeshT>> pv_bulk_out;
   std::unique_ptr<L2_FECollection> pv_bulk_sigma_fec;
   std::unique_ptr<PvFES> pv_bulk_sigma_fes;
   std::unique_ptr<PvGF>  pv_bulk_syy_gf, pv_bulk_sxy_gf, pv_bulk_sxz_gf;

   Vector pv_local_slip, pv_local_slip_rate, pv_local_traction;
   Vector pv_local_state, pv_local_normal_stress;
   Vector pv_local_a, pv_local_Dc, pv_local_x2, pv_local_x3;
   Vector pv_local_slip_rate_k4, pv_local_traction_k4, pv_local_normal_stress_k4;

   if (use_paraview)
   {
      if (rank == 0) { mkdir((output_dir + "/ParaView").c_str(), 0755); }
#ifdef MFEM_USE_MPI
      MPI_Barrier(comm);
#endif
      pv_out = std::make_unique<seas::ParaViewOutput<MeshT>>(
         output_dir + "/ParaView", pmesh, order);

      if (pv_low_order)
      {
         pv_out->SetHighOrderOutput(false);
         pv_out->SetLevelsOfDetail(1);
      }

      // Velocity: 3-component vector L2, byNODES so memcpy from Q's
      // [VX..VZ] block is a single contiguous copy.
      pv_vel_fec = std::make_unique<L2_FECollection>(order, 3, BasisType::GaussLobatto);
      pv_vel_fes = std::make_unique<PvFES>(&pmesh, pv_vel_fec.get(),
                                           3, Ordering::byNODES);
      pv_vel_gf  = std::make_unique<PvGF>(pv_vel_fes.get());
      *pv_vel_gf = 0.0;
      pv_out->RegisterDomainField("velocity", pv_vel_gf.get());

      // MPI rank: L2 p=0 (one value per element).
      pv_rank_fec = std::make_unique<L2_FECollection>(0, 3);
      pv_rank_fes = std::make_unique<PvFES>(&pmesh, pv_rank_fec.get());
      pv_rank_gf  = std::make_unique<PvGF>(pv_rank_fes.get());
      *pv_rank_gf = static_cast<real_t>(rank);
      pv_out->RegisterDomainField("mpi_rank", pv_rank_gf.get());

      // Fault L2-p0 fields keyed off the wave operator's canonical face
      // lists (interior + shared) so ParaView indexing matches DOFData.
      pv_out->InitFaultOutputBP5(fault_int_faces, fault_shr_faces,
                                 nqp_per_face);

      pv_local_slip.SetSize(2 * num_fault_total);
      pv_local_slip_rate.SetSize(2 * num_fault_total);
      pv_local_traction.SetSize(2 * num_fault_total);
      pv_local_state.SetSize(num_fault_total);
      pv_local_normal_stress.SetSize(num_fault_total);
      pv_local_slip_rate_k4.SetSize(2 * num_fault_total);
      pv_local_traction_k4.SetSize(2 * num_fault_total);
      pv_local_normal_stress_k4.SetSize(num_fault_total);

      pv_local_a.SetSize(num_fault_total);
      pv_local_Dc.SetSize(num_fault_total);
      pv_local_x2.SetSize(num_fault_total);
      pv_local_x3.SetSize(num_fault_total);
      for (int i = 0; i < num_fault_total; i++)
      {
         pv_local_a(i)  = dof_data[i].a;
         pv_local_Dc(i) = dof_data[i].Dc;
         pv_local_x2(i) = fault_coords[i](0);
         pv_local_x3(i) = fault_coords[i](2);
      }
      pv_out->SetFaultParamsBP5(pv_local_a, pv_local_Dc,
                                pv_local_x2, pv_local_x3);

      // Schedule (CLI > step-interval > output_dt step interval):
      const int output_interval_for_pv =
         std::max(1, static_cast<int>(output_dt / dt));
      if (paraview_step_interval > 0)
      {
         pv_out->output_every_n_steps = paraview_step_interval;
      }
      else if (paraview_dt_flag > 0.0)
      {
         pv_out->fixed_dt = paraview_dt_flag;
      }
      else
      {
         pv_out->output_every_n_steps = output_interval_for_pv;
      }

      if (rank == 0)
      {
         std::cout << "ParaView output: ON (prefix="
                   << output_dir << "/ParaView)\n";
         if (pv_no_domain)
         {
            std::cout << "  Mode: fault-surface PVD only (--no-domain-pv)\n";
         }
         if (paraview_step_interval > 0)
         {
            std::cout << "  Interval: every " << paraview_step_interval
                      << " steps (--paraview-every)\n";
         }
         else if (paraview_dt_flag > 0.0)
         {
            std::cout << "  Interval: every " << paraview_dt_flag
                      << " s (--paraview-dt)\n";
         }
         else
         {
            std::cout << "  Interval: every " << output_interval_for_pv
                      << " steps (matches --output-dt=" << output_dt
                      << " s)\n";
         }
      }

      if (paraview_bulk_dt > 0.0)
      {
         pv_bulk_out = std::make_unique<seas::ParaViewOutput<MeshT>>(
            output_dir + "/ParaView_bulk", pmesh, order);
         if (pv_low_order)
         {
            pv_bulk_out->SetHighOrderOutput(false);
            pv_bulk_out->SetLevelsOfDetail(1);
         }
         pv_bulk_out->RegisterDomainField("velocity", pv_vel_gf.get());
         pv_bulk_out->RegisterDomainField("mpi_rank", pv_rank_gf.get());

         pv_bulk_sigma_fec = std::make_unique<L2_FECollection>(
            order, 3, BasisType::GaussLobatto);
         pv_bulk_sigma_fes = std::make_unique<PvFES>(&pmesh,
                                                     pv_bulk_sigma_fec.get());
         pv_bulk_syy_gf = std::make_unique<PvGF>(pv_bulk_sigma_fes.get());
         pv_bulk_sxy_gf = std::make_unique<PvGF>(pv_bulk_sigma_fes.get());
         pv_bulk_sxz_gf = std::make_unique<PvGF>(pv_bulk_sigma_fes.get());
         *pv_bulk_syy_gf = 0.0;
         *pv_bulk_sxy_gf = 0.0;
         *pv_bulk_sxz_gf = 0.0;
         pv_bulk_out->RegisterDomainField("sigma_yy", pv_bulk_syy_gf.get());
         pv_bulk_out->RegisterDomainField("sigma_xy", pv_bulk_sxy_gf.get());
         pv_bulk_out->RegisterDomainField("sigma_xz", pv_bulk_sxz_gf.get());

         pv_bulk_out->fixed_dt = paraview_bulk_dt;

         if (rank == 0)
         {
            std::cout << "  Bulk collection: ON (prefix="
                      << output_dir << "/ParaView_bulk, every "
                      << paraview_bulk_dt
                      << " s; fields: velocity, sigma_yy, sigma_xy, "
                      << "sigma_xz, mpi_rank)\n";
         }
      }
   }

   // paraview_write: MPI-collective snapshot writer.  V_max must already be
   // globally reduced; step_num/time identify the frame.
   auto paraview_write = [&](int step_num, real_t time, real_t V_max)
   {
      if (!pv_out) { return; }

      const bool fault_wants =
         pv_out->PeekShouldWrite(step_num, time, V_max);
      const bool bulk_wants = pv_bulk_out &&
         pv_bulk_out->PeekShouldWrite(step_num, time, V_max);
      if (!fault_wants && !bulk_wants) { return; }

      if (!pv_no_domain || bulk_wants)
      {
         std::memcpy(pv_vel_gf->GetData(),
                     Q.GetData() + VX * ndof_total,
                     3 * ndof_total * sizeof(real_t));
      }

      if (bulk_wants)
      {
         std::memcpy(pv_bulk_syy_gf->GetData(),
                     Q.GetData() + SYY * ndof_total,
                     ndof_total * sizeof(real_t));
         std::memcpy(pv_bulk_sxy_gf->GetData(),
                     Q.GetData() + SXY * ndof_total,
                     ndof_total * sizeof(real_t));
         std::memcpy(pv_bulk_sxz_gf->GetData(),
                     Q.GetData() + SXZ * ndof_total,
                     ndof_total * sizeof(real_t));
         pv_bulk_out->ForceSave(step_num, time);
      }

      if (!fault_wants) { return; }

      for (int i = 0; i < num_fault_total; i++)
      {
         const DOFData &d = dof_data[i];
         pv_local_slip(2*i + 0)      = d.slip1;       // dip
         pv_local_slip(2*i + 1)      = d.slip2;       // strike
         pv_local_slip_rate(2*i + 0) = d.V1;
         pv_local_slip_rate(2*i + 1) = d.V2;
         pv_local_traction(2*i + 0)  = d.tau1_corr;
         pv_local_traction(2*i + 1)  = d.tau2_corr;
         pv_local_state(i)           = d.psi;
         pv_local_normal_stress(i)   = d.sigma_n_corr;

         // ADER one-shot: no stage-4 distinct from averaged DOFData (R7-001
         // option b).  Mirror tpv102 ADER path: write the same values into
         // _k4 buffers so the VTU delta reads exactly zero.
         pv_local_slip_rate_k4(2*i + 0) = d.V1;
         pv_local_slip_rate_k4(2*i + 1) = d.V2;
         pv_local_traction_k4(2*i + 0)  = d.tau1_corr;
         pv_local_traction_k4(2*i + 1)  = d.tau2_corr;
         pv_local_normal_stress_k4(i)   = d.sigma_n_corr;
      }

      if (pv_no_domain)
      {
         pv_out->CommitSchedule(time);
      }
      else
      {
         pv_out->UpdateFaultFieldsBP5(pv_local_slip, pv_local_slip_rate,
                                      pv_local_traction, pv_local_state,
                                      pv_local_normal_stress);
         pv_out->ForceSave(step_num, time);
         // ForceSave only advances last_write_time_; the V_max-adaptive
         // schedule additionally needs current_regime_ / last_v_max_
         // advanced for the next PeekShouldWrite to use the correct
         // regime interval (paraview_output.hpp:985-990).
         pv_out->CommitSchedule(time, V_max);
      }

      pv_out->WriteFaultSurfaceVTU(
         output_dir, step_num, time, rank, nprocs,
         pv_local_slip, pv_local_slip_rate, pv_local_traction,
         pv_local_state, pv_local_normal_stress,
         pv_local_a, pv_local_Dc, pv_local_x2, pv_local_x3,
         pv_local_slip_rate_k4, pv_local_traction_k4,
         pv_local_normal_stress_k4);
   };

   // Initial snapshot at t=0 (V_max = V_ini since fault is quasi-static).
   paraview_write(0, 0.0, TPV102Params::V_ini);

   // If tfinal == 0: init-only run, stations at t=0 already written.
   // Skip the time loop and go straight to summary.  This is the
   // Phase-2 P2_D gate — the init sbatch verifies ψ_ini at every
   // station from the single t=0 row.
   if (nsteps == 0)
   {
      if (rank == 0)
      {
         std::cout << "[tpv102_driver] tfinal = 0 — init-only run. "
                   << "Station t=0 row written, exiting.\n";
      }
      station_writer.Flush();
      station_writer.Close();
      surface_writer.Flush();
      surface_writer.Close();
#ifdef MFEM_USE_MPI
      MPI_Finalize();
#endif
      return 0;
   }

   // -----------------------------------------------------------------------
   // 8. ADER time stepping (plan §4.10 Step 9: ADER-only, no RK4 path)
   // -----------------------------------------------------------------------
   int output_interval = std::max(1, static_cast<int>(output_dt / dt));

   real_t t = 0.0;
   real_t V_max_global = 0.0;
   Vector Q_new(Q.Size());
   std::vector<real_t> psi_n(num_fault_total);

   if (rank == 0)
   {
      std::cout << "Starting ADER-O(" << ader_order << ") time loop...\n";
   }

   for (int step = 0; step < nsteps; ++step)
   {
      real_t dt_step = std::min(dt, tfinal - t);
      if (dt_step <= 0.0) { break; }

      // Save ψ at step start for the analytic update below.
      for (int i = 0; i < num_fault_total; ++i)
      {
         psi_n[i] = dof_data[i].psi;
      }

      // Nucleation (§3.9 directive): cumulative per-macro-step increment.
      // ApplyNucleationIncremental_TPV102 adds  Δτ · smoothStepIncrement
      // to `tau2_nuc` at every call, so over [0, T_nuc] the channel
      // telescopes to the full perturbation.
      //
      // R-1008 (round-10) NUCLEATION DOUBLE-COUNT FIX:
      //   The substep path's Tpv102SubStepIterator::AdvanceWithSubStepStates
      //   (and legacy Advance) already calls ApplyNucleationIncremental_TPV102
      //   ONCE PER SUB-STEP internally (line 295 / its AdvanceWithSubStepStates
      //   sibling).  Σ_o ΔS over the sub-steps telescopes to the same
      //   macro-step increment ΔS(t+dt) − ΔS(t).  If the driver ALSO
      //   calls the accumulator here, tau2_nuc is incremented TWICE per
      //   macro-step → 2× the spec Δτ₀ at full ramp → 2-million× rupture
      //   over-acceleration (terminal velocity by t=0.5 s instead of
      //   ~t=1 s).  Skip the driver-level call when the iterator owns
      //   nucleation cadence; one-shot path keeps the driver-level call.
      if (!disable_nucleation && num_fault_total > 0 && !use_substep_iterator)
      {
         ApplyNucleationIncremental_TPV102(dof_data, fault_coords,
                                           t + dt_step, dt_step);
      }

      // ADER predictor-corrector.  Default (use_substep_iterator==false):
      // one-shot AdvanceADER runs the bulk wave update + the fault-face
      // Riemann solve (through FaultFaceFlux::EvaluateADER); the solve
      // reads DOFData.psi and DOFData.tau*_nuc as set above, and writes
      // V1/V2/slip_rate/tau*_corr/sigma_n_corr back onto DOFData.
      //
      // Sub-step path (use_substep_iterator==true, --fault-iterator
      // substep):  AdvanceADERWithSubStep composes ComputeADERSubStepStates
      // (per-sub-step pointwise Q via Taylor expansion),
      // EvaluateBulkAtFaultQPsCanonical (per-sub-step canonical Q at fault
      // QPs), Tpv102SubStepIterator::AdvanceWithSubStepStates (per-sub-step
      // friction + ψ + slip + accumulated I_imp), and AdvanceADER (with
      // the iterator's I_imp installed via SetSubStepFaultImposedStates so
      // the fault branch consumes them in lieu of inline EvaluateADER).
      // At O=1 the two paths are bit-identical (T_TPV102_SSI_3 contract).
      if (use_substep_iterator)
      {
         AdvanceADERWithSubStep(wave, substep_iterator, dof_data,
                                fault_coords, Q, dt_step,
                                ader_order, /*t_step_start=*/t,
                                method, Q_new);
      }
      else
      {
         wave.AdvanceADER(Q, dt_step, ader_order, Q_new);
      }
      Q.Swap(Q_new);
      t += dt_step;

#ifdef SEAS_DIAG_FAULT_FLUX
      // C-2C BULK-DELTA: per-macro-step delta of bulk Q at the diag
      // DOFs.  Captures the NET change Q^{n+1} - Q^n at the two diag
      // tets' face DOFs after one full ADER step (predictor + corrector
      // including fault Riemann + bulk wave op).  Lets the post-run
      // analyzer correlate per-step asymmetry growth against C-1n.
      if (diag_elem_plus >= 0 && diag_elem_minus >= 0 &&
          diag_face_dof_plus >= 0 && diag_face_dof_minus >= 0)
      {
         static std::vector<real_t> prev_Qp(NUM_STATE, 0.0);
         static std::vector<real_t> prev_Qm(NUM_STATE, 0.0);
         static bool have_prev = false;

         const real_t *Qd = Q.GetData();
         const int dof_off_p = diag_elem_plus  * wave.GetNDof();
         const int dof_off_m = diag_elem_minus * wave.GetNDof();
         real_t Qp[NUM_STATE], Qm[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            Qp[c] = Qd[c*ndof_total + dof_off_p + diag_face_dof_plus];
            Qm[c] = Qd[c*ndof_total + dof_off_m + diag_face_dof_minus];
         }
         real_t dQp[NUM_STATE], dQm[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            dQp[c] = have_prev ? (Qp[c] - prev_Qp[c]) : 0.0;
            dQm[c] = have_prev ? (Qm[c] - prev_Qm[c]) : 0.0;
            prev_Qp[c] = Qp[c];
            prev_Qm[c] = Qm[c];
         }
         have_prev = true;
         std::fprintf(stderr,
            "[C-2C BULK-DELTA] rank=%d t=%.4e  "
            "Q+_SXX=%+.4e Q-_SXX=%+.4e diff_SXX=%+.4e  "
            "dQ+_SXX=%+.4e dQ-_SXX=%+.4e d_diff_SXX=%+.4e  "
            "Q+_SYY=%+.4e Q-_SYY=%+.4e diff_SYY=%+.4e  "
            "Q+_SZZ=%+.4e Q-_SZZ=%+.4e diff_SZZ=%+.4e\n",
            g_seas_my_rank, t,
            Qp[SXX], Qm[SXX], Qp[SXX] - Qm[SXX],
            dQp[SXX], dQm[SXX], dQp[SXX] - dQm[SXX],
            Qp[SYY], Qm[SYY], Qp[SYY] - Qm[SYY],
            Qp[SZZ], Qm[SZZ], Qp[SZZ] - Qm[SZZ]);
      }
#endif

      // ψ update + slip accumulation.  TPV102 uses the regularised
      // rate-and-state ageing law (SCEC TPV101/102 Eq. (2): dθ/dt = 1 −
      // Vθ/L).  `UpdateStateAnalytic` from friction/state_evolution.hpp
      // is the exact ageing-law analytic update for constant V over dt:
      //   θ(t+dt) = θ·exp(−Vdt/L) + (L/V)·(1 − exp(−Vdt/L))
      // converted back to ψ-space.
      //
      // The iterator's per-sub-step path already accumulates slip1/slip2
      // and writes ψ per sub-step (R-1009 slip double-count avoidance).
      // On the one-shot path the driver owns both the slip integration
      // and the macro-step ψ update from psi_n[i].
      for (int i = 0; i < num_fault_total; ++i)
      {
         dof_data[i].psi = UpdateStateAnalytic(
            psi_n[i],
            dof_data[i].slip_rate,
            dof_data[i].Dc,
            dt_step,
            TPV102Params::f0,
            TPV102Params::b,
            TPV102Params::V0);
         if (!use_substep_iterator)
         {
            dof_data[i].slip1 += dof_data[i].V1 * dt_step;
            dof_data[i].slip2 += dof_data[i].V2 * dt_step;
         }
      }

      // V_max tracking.
      real_t V_max_local = 0.0;
      for (int i = 0; i < num_fault_total; ++i)
      {
         V_max_local = std::max(V_max_local, dof_data[i].slip_rate);
      }
#ifdef MFEM_USE_MPI
      real_t V_max_step;
      MPI_Allreduce(&V_max_local, &V_max_step, 1,
                    MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
#else
      real_t V_max_step = V_max_local;
#endif
      V_max_global = std::max(V_max_global, V_max_step);

      // Output cadence.
      if (step % output_interval == 0 || step == nsteps - 1)
      {
         station_writer.WriteStep(t, dof_data);
         surface_writer.WriteStep(t, Q);
         if (rank == 0)
         {
            std::cout << "Step " << step << "/" << nsteps
                      << ", t = " << t << " s"
                      << ", V_max = " << V_max_step << " m/s\n";
         }
      }

      // ParaView snapshot — has its own internal schedule (matches
      // --paraview-dt / --paraview-every / output_dt), so we always
      // call and let pv_out / pv_bulk_out's PeekShouldWrite gate the I/O.
      paraview_write(step + 1, t, V_max_step);

      if (step == 0)
      {
         wave.VerifySharedFaultDOFDataConsistency();
      }

      // NaN tripwire.
      {
         real_t local_nan = std::isnan(Q.Norml2()) ? 1.0 : 0.0;
         real_t global_nan = local_nan;
#ifdef MFEM_USE_MPI
         MPI_Allreduce(&local_nan, &global_nan, 1,
                       MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
#endif
         if (global_nan > 0.0)
         {
            std::cerr << "ERROR: NaN detected at step " << step
                      << ", t = " << t << " s (rank " << rank << ")\n";
#ifdef MFEM_USE_MPI
            MPI_Finalize();
#endif
            return 1;
         }
      }

      if (debug_qnorm && step % output_interval == 0 && rank == 0)
      {
         std::cout << "  [qnorm] ||Q||_2 = " << Q.Norml2() << "\n";
      }
   }

   // -----------------------------------------------------------------------
   // 9. Summary
   // -----------------------------------------------------------------------
   station_writer.Flush();
   station_writer.Close();
   surface_writer.Flush();
   surface_writer.Close();

   // Diagnostic dump: end-of-run per-fault-QP state, gated by env var.
   // Writes <output_dir>/fault_qp_dump_rank<R>.txt, one row per local
   // fault QP.  No MPI calls; each rank writes only its own file.
   {
      const char *dump_env = std::getenv("SEAS_DIAG_DUMP_FAULT_QPS");
      if (dump_env && dump_env[0] != '\0' &&
          !(dump_env[0] == '0' && dump_env[1] == '\0'))
      {
         std::ostringstream pathss;
         pathss << output_dir << "/fault_qp_dump_rank" << rank << ".txt";
         std::ofstream fs(pathss.str());
         if (fs.is_open())
         {
            fs << std::scientific << std::setprecision(17);
            fs << "# columns: dof_idx x y z V1 V2 slip1 slip2 "
                  "tau1_corr tau2_corr sigma_n_corr psi\n";
            const int n = static_cast<int>(dof_data.size());
            for (int i = 0; i < n; i++)
            {
               const Vector &c = fault_coords[i];
               const DOFData &d = dof_data[i];
               fs << i << " " << c(0) << " " << c(1) << " " << c(2)
                  << " " << d.V1 << " " << d.V2
                  << " " << d.slip1 << " " << d.slip2
                  << " " << d.tau1_corr << " " << d.tau2_corr
                  << " " << d.sigma_n_corr << " " << d.psi << "\n";
            }
            fs.close();
            if (rank == 0)
            {
               std::cout << "[diag-dump] wrote " << n
                         << " fault QPs to fault_qp_dump_rank<R>.txt\n";
            }
         }
      }
   }

   if (rank == 0)
   {
      std::cout << "\n========================================\n";
      std::cout << "TPV102 run complete.\n";
      std::cout << "  tfinal  = " << tfinal << " s\n";
      std::cout << "  steps   = " << nsteps << "\n";
      std::cout << "  V_max   = " << V_max_global << " m/s\n";
      std::cout << "  output  = " << output_dir << "/\n";
      std::cout << "========================================\n";
   }

#ifdef MFEM_USE_MPI
   MPI_Finalize();
#endif
   return 0;
}

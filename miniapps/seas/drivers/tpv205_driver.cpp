// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// SCEC TPV5 / TPV205 benchmark driver — 3D dynamic rupture on a vertical
// right-lateral strike-slip fault in a half-space, with linear slip-
// weakening (LSW) friction and four pre-stress patches (background
// 70 MPa, nucleation 81.6 MPa, left 78 MPa, right 62 MPa).
//
// Fully parallel MPI driver using ParMesh and WaveOperator<ParMesh>.
// Mirrors drivers/tpv104_driver.cpp's ADER time-stepping flow with the
// LSW-specific substitutions:
//   - TPV205Params material + LSW friction constants.
//   - InitializeFaultDOFs_TPV205 (per-QP μ_s/μ_d/d_c via patch +
//     strength-barrier lookup; per-QP τ_strike(x, z) via patch lookup).
//   - No state-evolution instance (LSW has no ψ).
//   - No nucleation accumulator (TPV205 nucleates by static patch
//     pre-stress; the τ_strike field is at its full patch value from
//     t = 0 with no time-varying perturbation).
//   - Tpv205SubStepIterator does the LSW closed-form solve internally
//     using FaultFaceFlux::ComputeTrialTraction + LSW + BuildImposedState.
//   - All ADER dispatches route through AdvanceADERWithSubStep (the
//     legacy wave.AdvanceADER fault branch only supports rate-and-state
//     friction via Brent — not applicable for LSW).
//   - TPV205StationWriter / TPV205SurfaceStationWriter (SCEC trace
//     column order; final column is μ_eff(δ) instead of ψ).
//
// Usage:
//   ibrun ./seas_tpv205_driver \
//         --mesh tpv205/mesh/tpv2053d_200m.msh \
//         --tfinal 12.0 --ader-order 2 \
//         --friction-solver newton-stable \
//         --output-dir tpv205/results

#include "mfem.hpp"
#include "../dynamic/wave_state.hpp"
#include "../dynamic/wave_operator.hpp"
#include "../dynamic/fault_face_flux.hpp"
#include "../dynamic/friction_solver.hpp"
#include "../dynamic/tpv205_setup.hpp"
#include "../dynamic/tpv205_friction.hpp"
#include "../dynamic/tpv205_substep_iterator.hpp"
#include "../config/tpv205_params.hpp"
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
// R7-001/R7-005 note: TPV205 production runs LSW closed-form via
// Tpv205SubStepIterator + AdvanceADERWithSubStep; the wave operator's
// fault dispatch routes through FaultFaceFlux::EvaluateADER_LSW on
// shared-fault QPs (R-1601 fallback).  The --friction-solver CLI is
// accepted for launch-script symmetry with TPV102/TPV104 but is
// intentionally ignored (LSW has no root finder).
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
enum class DispatchedSolver { Brent, NewtonRaphsonStable, NewtonRaphsonLegacy, Hybrid, LSWClosedForm };
enum class DispatchedIterator { OneShot, SubStep };
enum class DispatchedLaw { SlipSRW, Aging, LSW };

static DispatchedSolver GetDispatchedSolver(const std::string &/*friction_solver_cli*/)
{
   // R-004: TPV205 always uses the LSW closed-form solve via
   // Tpv205SubStepIterator + AdvanceADERWithSubStep (no Newton, no
   // Brent, no rate-and-state EvaluateADER).  The --friction-solver
   // CLI is accepted for launch-script symmetry with TPV102/TPV104
   // but intentionally ignored — banner reflects the actual dispatch.
   return DispatchedSolver::LSWClosedForm;
}
static DispatchedIterator GetDispatchedIterator(const std::string &fault_iterator_cli)
{
   // R-602/R-603 (round-7 implementation): the CLI value now controls
   // dispatch.  "substep" routes the time loop through the
   // Tpv205SubStepIterator + per-sub-step ADER predictor path; any other
   // value (including "one-shot") keeps the legacy single-shot
   // wave.AdvanceADER call.  Default unchanged: one-shot.
   if (fault_iterator_cli == "substep")
   {
      return DispatchedIterator::SubStep;
   }
   return DispatchedIterator::OneShot;
}
static DispatchedLaw GetDispatchedLaw(const std::string &fric_law_cli)
{
   // R-411: previously the CLI value was silently ignored, so
   // `--fric-law aging` on TPV205 mapped to LSW with no warning.  Hard-
   // reject any value other than the empty default or "lsw" so typos
   // surface at parse time.  TPV205 only implements LSW closed-form
   // (SCEC TPV5 §7-11); the flag is still accepted (and routed to LSW)
   // for launch-script symmetry with TPV102/TPV104.
   if (!fric_law_cli.empty() && fric_law_cli != "lsw")
   {
      MFEM_ABORT("TPV205 only supports --fric-law lsw (linear slip-"
                 "weakening, SCEC TPV5 §7-11 closed-form).  Got: '"
                 << fric_law_cli << "'.  Use --fric-law lsw or omit "
                 "the flag.");
   }
   return DispatchedLaw::LSW;
}

static const char *TagOf(DispatchedSolver s)
{
   switch (s)
   {
      case DispatchedSolver::Brent:                return "brent";
      case DispatchedSolver::NewtonRaphsonStable:  return "newton-stable";
      case DispatchedSolver::NewtonRaphsonLegacy:  return "newton-legacy";
      case DispatchedSolver::Hybrid:               return "hybrid";
      case DispatchedSolver::LSWClosedForm:        return "lsw-closed-form";
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
      case DispatchedLaw::LSW:     return "lsw";
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
      case DispatchedSolver::LSWClosedForm:
         return "LSW closed-form (no root-finder; "
                "--friction-solver flag IGNORED)";
   }
   return "unknown";
}
static std::string BannerOf(DispatchedIterator i)
{
   switch (i)
   {
      case DispatchedIterator::OneShot:
         return "one-shot (iterator with O = 1; LSW closed-form per macro-step)";
      case DispatchedIterator::SubStep:
         return "sub-step (Tpv205SubStepIterator + per-sub-step ADER "
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
         return "aging (ψ-space, forward-Euler)";
      case DispatchedLaw::LSW:
         return "lsw (linear slip-weakening; closed-form per sub-step, "
                "SCEC TPV5 §7-11)";
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
   mfem::seas::Tpv205SubStepIterator &iterator,
   std::vector<mfem::seas::DOFData> &dof_data,
   const std::vector<mfem::Vector> &fault_coords,
   const mfem::Vector &Q,
   mfem::real_t dt_step,
   int ader_order,
   mfem::real_t t_step_start,
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
   // sub-step's time-average — restoring the T_TPV205_SSI_3 contract
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
   // tpv205_substep_iterator.cpp:573-589 expects each Q_pointwise_*[o]
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
                                        I_imp_minus_flat.data());
   }

   // R-015: RAII guard so that an exception thrown inside
   // `wave.AdvanceADER` (or any subsequent driver code that runs before
   // the explicit Reset) cannot leave stale `I_imp_*_flat` pointers
   // dangling in the wave operator.  The flat buffers are stack-vectors
   // local to this function — once we return, those pointers become
   // invalid and the next call into `wave.AdvanceADER` would read freed
   // memory.  The guard's destructor calls Reset whether the scope
   // exits normally or via exception.
   struct SubStepFaultImposedGuard
   {
      WaveOpT &wave_;
      explicit SubStepFaultImposedGuard(WaveOpT &w) : wave_(w) {}
      ~SubStepFaultImposedGuard() { wave_.ResetSubStepFaultImposedStates(); }
   };

   // Hand the iterator's output to the wave op so the upcoming
   // AdvanceADER's fault branch substitutes it for inline EvaluateADER.
   wave.SetSubStepFaultImposedStates(
      n_total_fault_qps > 0 ? I_imp_plus_flat.data()  : nullptr,
      n_total_fault_qps > 0 ? I_imp_minus_flat.data() : nullptr,
      n_total_fault_qps);
   SubStepFaultImposedGuard imposed_guard(wave);

   // Bulk corrector: runs unchanged for non-fault faces; fault branch
   // consumes the side-channel imposed states.  `imposed_guard`'s
   // destructor calls ResetSubStepFaultImposedStates when this scope
   // exits, regardless of normal return or exception path.
   wave.AdvanceADER(Q, dt_step, ader_order, Q_new);
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
#ifdef SEAS_DIAG_TPV205_STATE
      const char *diag_tpv205 = "SEAS_DIAG_TPV205_STATE = ON";
#else
      const char *diag_tpv205 = "SEAS_DIAG_TPV205_STATE = OFF";
#endif
      std::fprintf(stderr, "[BUILD] %s\n", diag_fault_flux);
      std::fprintf(stderr, "[BUILD] %s\n", diag_tpv205);
      std::ofstream binfo("build_info.txt");
      if (binfo.is_open())
      {
         binfo << "[BUILD] " << diag_fault_flux << "\n";
         binfo << "[BUILD] " << diag_tpv205 << "\n";
         binfo.close();
      }
   }

   // -----------------------------------------------------------------------
   // Parse command-line arguments
   // -----------------------------------------------------------------------
   // R-012: default matches the .geo-produced filename
   // (tpv205/mesh/tpv2053d_200m.geo → tpv2053d_200m.msh, with the `3d`
   // infix).  Pre-fix the default referenced a `tpv205_200m.msh` path
   // that gmsh never produces; first invocation aborted on file-not-
   // found.  README usage line at the top of this file is the source
   // of truth.
   std::string mesh_file = GetStringArg(argc, argv, "--mesh",
                                        "tpv205/mesh/tpv2053d_200m.msh");
   real_t mesh_scale     = GetRealArg(argc, argv, "--mesh-scale", 1.0);
   int order             = GetIntArg(argc, argv, "--order", 1);
   std::string bc_mode   = GetStringArg(argc, argv, "--bc-mode", "absorbing");
   real_t tfinal         = GetRealArg(argc, argv, "--tfinal", TPV205Params::t_final);
   std::string output_dir    = GetStringArg(argc, argv, "--output-dir", "tpv205/results");
   std::string output_prefix = GetStringArg(argc, argv, "--output-prefix", "tpv205");
   real_t cfl_factor     = GetRealArg(argc, argv, "--cfl", 0.5);
   // TPV205 mesh (tpv2053d_*.geo) emits Physical Surface 101/103/105.
   int bc_free           = GetIntArg(argc, argv, "--bc-free",
                                     TPV205Params::bc_free_default);
   int bc_fault          = GetIntArg(argc, argv, "--bc-fault",
                                     TPV205Params::bc_fault_default);
   int bc_absorb         = GetIntArg(argc, argv, "--bc-absorb",
                                     TPV205Params::bc_absorb_default);
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
   //   --no-volume-pv / --no-domain-pv (deprecated alias) :
   //                             suppress fault-schedule volume save.
   //   --volume-pv-dt X        : volume cadence override (s).  Re-enables
   //                             the volume save even with --no-volume-pv.
   //   --paraview-fault-{vtu,hdf5,legacy-ascii} : back-end selection.
   //   --paraview-{fault,bulk}-{zfp-tol,deflate-level} : VTKHDF chunk filter.
   //   --paraview-max-snapshots N / --paraview-{co,nucleation,inter}seismic-dt X :
   //                             snapshot cap and per-regime cadences.
   //   See tpv102_driver.cpp for the canonical doc-block.
   bool use_paraview = HasFlag(argc, argv, "--paraview");
   bool pv_low_order = HasFlag(argc, argv, "--pv-low-order");
   // Phase 4: see tpv102_driver.cpp for the canonical doc-block on
   // --no-volume-pv / --volume-pv-dt.
   bool pv_no_domain = HasFlag(argc, argv, "--no-domain-pv")
                       || HasFlag(argc, argv, "--no-volume-pv");
   real_t volume_pv_dt = GetRealArg(argc, argv, "--volume-pv-dt", 0.0);
   int  paraview_step_interval = GetIntArg(argc, argv, "--paraview-every", 0);
   real_t paraview_dt_flag     = GetRealArg(argc, argv, "--paraview-dt", 0.0);
   real_t paraview_bulk_dt     = GetRealArg(argc, argv, "--paraview-bulk-dt", 0.0);
   if (paraview_step_interval > 0 || paraview_dt_flag > 0.0
       || paraview_bulk_dt > 0.0)
   {
      use_paraview = true;
   }

   // R-101 / PLAN_paraview_compaction_2026-04-28 §Phase 2b / 2d.3 / 3
   // CLI flags.  See tpv102_driver.cpp for the canonical doc-block.
   const bool   paraview_force_vtu       = HasFlag(argc, argv, "--paraview-fault-vtu");
   const bool   paraview_force_hdf5      = HasFlag(argc, argv, "--paraview-fault-hdf5");
   // R-305: legacy per-rank ASCII back end (debugging only).
   const bool   paraview_legacy_ascii    = HasFlag(argc, argv, "--paraview-fault-legacy-ascii");
   const real_t paraview_fault_zfp_tol   = GetRealArg(argc, argv, "--paraview-fault-zfp-tol",   0.0);
   const int    paraview_fault_deflate   = GetIntArg (argc, argv, "--paraview-fault-deflate-level", -1);
   // Phase 2d.3: bulk-side flags (see tpv102_driver.cpp for canonical comment).
   const real_t paraview_bulk_zfp_tol    = GetRealArg(argc, argv, "--paraview-bulk-zfp-tol",   0.0);
   const int    paraview_bulk_deflate    = GetIntArg (argc, argv, "--paraview-bulk-deflate-level", -1);
   // Phase 6.3 / 6.3a: volume back-end + primary-collection compression.
   const bool   paraview_volume_force_vtu  = HasFlag(argc, argv, "--paraview-volume-vtu");
   const bool   paraview_volume_force_hdf5 = HasFlag(argc, argv, "--paraview-volume-hdf5");
   const real_t paraview_volume_zfp_tol    = GetRealArg(argc, argv, "--paraview-volume-zfp-tol",   0.0);
   const int    paraview_volume_deflate    = GetIntArg (argc, argv, "--paraview-volume-deflate-level", -1);
   const int    paraview_max_snapshots   = GetIntArg (argc, argv, "--paraview-max-snapshots", 0);
   const real_t paraview_coseismic_dt    = GetRealArg(argc, argv, "--paraview-coseismic-dt",   -1.0);
   const real_t paraview_nucleation_dt   = GetRealArg(argc, argv, "--paraview-nucleation-dt",  -1.0);
   const real_t paraview_interseismic_dt = GetRealArg(argc, argv, "--paraview-interseismic-dt",-1.0);
   if (paraview_force_vtu || paraview_force_hdf5
       || paraview_fault_zfp_tol > 0.0
       || paraview_fault_deflate >= 0
       || paraview_bulk_zfp_tol > 0.0
       || paraview_bulk_deflate >= 0
       || paraview_volume_force_vtu || paraview_volume_force_hdf5
       || paraview_volume_zfp_tol > 0.0
       || paraview_volume_deflate >= 0
       || paraview_max_snapshots > 0
       || paraview_coseismic_dt > 0.0
       || paraview_nucleation_dt > 0.0
       || paraview_interseismic_dt > 0.0)
   {
      use_paraview = true;
   }
   if (paraview_force_vtu && paraview_force_hdf5)
   {
      MFEM_ABORT("--paraview-fault-vtu and --paraview-fault-hdf5 are "
                 "mutually exclusive.");
   }
   if (paraview_legacy_ascii && paraview_force_hdf5)
   {
      MFEM_ABORT("--paraview-fault-legacy-ascii implies the binary VTU "
                 "back end and is incompatible with --paraview-fault-hdf5.");
   }
#ifndef MFEM_USE_HDF5
   if (paraview_force_hdf5)
   {
      MFEM_ABORT("--paraview-fault-hdf5 requires the seas-mfem build to "
                 "define MFEM_USE_HDF5=YES; current build has it disabled.");
   }
#endif
#ifndef MFEM_USE_H5Z_ZFP
   if (paraview_fault_zfp_tol > 0.0)
   {
      MFEM_ABORT("--paraview-fault-zfp-tol requires the seas-mfem build "
                 "to define MFEM_USE_H5Z_ZFP=YES; current build has it "
                 "disabled.");
   }
   if (paraview_bulk_zfp_tol > 0.0)
   {
      MFEM_ABORT("--paraview-bulk-zfp-tol requires the seas-mfem build "
                 "to define MFEM_USE_H5Z_ZFP=YES; current build has it "
                 "disabled.");
   }
#endif
#ifndef MFEM_USE_HDF5
   // R-306 / plan §Phase 2d.3: deflate-level flags require HDF5.
   if (paraview_fault_deflate >= 0)
   {
      MFEM_ABORT("--paraview-fault-deflate-level requires the seas-mfem "
                 "build to define MFEM_USE_HDF5=YES; current build has "
                 "it disabled.");
   }
   if (paraview_bulk_deflate >= 0)
   {
      MFEM_ABORT("--paraview-bulk-deflate-level requires the seas-mfem "
                 "build to define MFEM_USE_HDF5=YES; current build has "
                 "it disabled.");
   }
#endif
   if (paraview_fault_zfp_tol > 0.0 && paraview_fault_deflate >= 0)
   {
      MFEM_ABORT("--paraview-fault-zfp-tol and --paraview-fault-deflate-level "
                 "are mutually exclusive — choose ZFP-accuracy OR deflate, "
                 "not both.");
   }
   if (paraview_bulk_zfp_tol > 0.0 && paraview_bulk_deflate >= 0)
   {
      MFEM_ABORT("--paraview-bulk-zfp-tol and --paraview-bulk-deflate-level "
                 "are mutually exclusive — choose ZFP-accuracy OR deflate, "
                 "not both.");
   }
   // Phase 6.3 / 6.3a: volume back-end + volume compression validation.
   if (paraview_volume_force_vtu && paraview_volume_force_hdf5)
   {
      MFEM_ABORT("--paraview-volume-vtu and --paraview-volume-hdf5 are "
                 "mutually exclusive.");
   }
#ifndef MFEM_USE_HDF5
   if (paraview_volume_force_hdf5)
   {
      MFEM_ABORT("--paraview-volume-hdf5 requires the seas-mfem build to "
                 "define MFEM_USE_HDF5=YES; current build has it disabled.");
   }
   if (paraview_volume_zfp_tol > 0.0)
   {
      MFEM_ABORT("--paraview-volume-zfp-tol requires the seas-mfem build "
                 "to define MFEM_USE_HDF5=YES; current build has it disabled.");
   }
   if (paraview_volume_deflate >= 0)
   {
      MFEM_ABORT("--paraview-volume-deflate-level requires the seas-mfem "
                 "build to define MFEM_USE_HDF5=YES; current build has it "
                 "disabled.");
   }
#endif
#ifndef MFEM_USE_H5Z_ZFP
   if (paraview_volume_zfp_tol > 0.0)
   {
      MFEM_ABORT("--paraview-volume-zfp-tol requires MFEM_USE_H5Z_ZFP=YES; "
                 "current build has it disabled.");
   }
#endif
   if (paraview_volume_zfp_tol > 0.0 && paraview_volume_deflate >= 0)
   {
      MFEM_ABORT("--paraview-volume-zfp-tol and --paraview-volume-deflate-level "
                 "are mutually exclusive — choose ZFP-accuracy OR deflate, "
                 "not both.");
   }
   if (paraview_bulk_dt <= 0.0
       && (paraview_bulk_zfp_tol > 0.0 || paraview_bulk_deflate >= 0))
   {
      mfem::out
         << "warning: --paraview-bulk-zfp-tol / --paraview-bulk-deflate-level "
            "set but --paraview-bulk-dt not provided; the secondary bulk "
            "collection is disabled, the flag has no effect.\n";
   }
   // `--dry-run` is a shortcut for "no mesh, no time-stepping,
   // just print banner + verify wiring compiles/runs".  Used by
   // test_tpv205_smoke.cpp and the banner-check sbatch on Frontera.
   // `--verify-dispatch` (R7-004) additionally emits machine-readable
   // [dispatch] lines describing what the production time loop ACTUALLY
   // runs — the tri-consistency check between banner claims and runtime
   // dispatch.  Under R7-001 option (b) the dispatch is always
   // Brent/one-shot/slip-SRW, independent of --friction-solver /
   // --fault-iterator / --fric-law values.

   // TPV205-specific CLI (plan §4.10 Step 9).
   // R7-001/R7-003/R7-006: these flags are accepted so smoke tests and
   // sbatch scripts can exercise banner parity, but on the current
   // driver path (one-shot wave.AdvanceADER) they have NO effect on the
   // dispatched solver / iterator / friction law.  Runtime is:
   //   - friction solver: LSW closed-form (Tpv205SubStepIterator +
   //                      AdvanceADERWithSubStep; shared-fault QPs at
   //                      np > 1 dispatch through
   //                      FaultFaceFlux::EvaluateADER_LSW per the
   //                      R-1601 fallback in wave_operator.inl)
   //   - fault iterator : substep at O = ader_order, or one-shot at
   //                      O = 1 (both route through AdvanceADERWithSubStep
   //                      — TPV205 has no legacy wave.AdvanceADER fault
   //                      branch since LSW is not in the rate-and-state
   //                      EvaluateADER closure)
   //   - friction law   : LSW (linear slip-weakening; SCEC TPV5 §7-11)
   // The banner below mirrors this disclosure exactly.
   std::string friction_solver =
      GetStringArg(argc, argv, "--friction-solver", "newton-stable");
   std::string fric_law =
      GetStringArg(argc, argv, "--fric-law", "lsw");
   // Default = "one-shot": legacy wave.AdvanceADER dispatch (production
   // path, byte-identical to pre-R-602 behavior).  Pass --fault-iterator
   // substep to opt in to the per-sub-step Tpv205SubStepIterator path.
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
   // is read.  TPV205 currently doesn't expose
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
      std::cout << "SCEC TPV205 Dynamic Rupture Simulation\n";
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
                << " (Tpv205SubStepIterator + AdvanceADERWithSubStep)\n";
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
      std::cout << "Nucleation: static patch pre-stress at t = 0 "
                << "(no time-varying perturbation; SCEC TPV5 §7-9)\n";
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

   // Validate --friction-solver eagerly so typos abort before any
   // simulation work (including --dry-run).  TPV205 LSW does not use
   // the rate-and-state Newton/Brent solver, but the CLI flag is
   // accepted for symmetry with TPV102/TPV104 launch scripts; the
   // returned value is intentionally unused.
   (void) MapSolver(friction_solver);

   // --dry-run: no mesh, no simulation.  Print a canonical end-of-run
   // line that test_tpv205_smoke.cpp scrapes ("[dry-run] OK.").  Used
   // as the cheapest possible Frontera startup verification.
   if (dry_run)
   {
      if (rank == 0)
      {
         std::cout << "[tpv205_driver] --dry-run: banner printed, "
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
   MFEM_VERIFY(serial_mesh.Dimension() == 3, "TPV205 requires 3D mesh");
   if (mesh_scale != 1.0)
   {
      serial_mesh.SetCurvature(1, false, 3, Ordering::byVDIM);
      Vector &nodes = *serial_mesh.GetNodes();
      nodes *= mesh_scale;
   }

#ifdef MFEM_USE_MPI
   // G1 fault-locality partition (TPV205 dynamic only).  See
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
      // produced by tpv205/mesh/build_symmirror_mesh.py --emit-partition).
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
      // bc_fault default is 3 per TPV205 driver; same convention used by
      // wave_operator.inl when populating fault_interior_faces_.
      const int bc_fault_attr = GetIntArg(argc, argv, "--bc-fault",
                                          TPV205Params::bc_fault_default);
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
                            TPV205Params::lambda, TPV205Params::mu,
                            TPV205Params::rho, bc);

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
   // 5. Fault DOF data (local partition) — TPV205-specific init
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

#ifdef SEAS_DIAG_TPV205_FAULT_BASIS
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
      // TPV205 fluctuation-Q init: Q = 0, patch-dependent pre-stress
      // lives in DOFData (tau2_0).  No V_w side-channel — TPV205 uses
      // linear slip-weakening, not the FVW slip law.
      InitializeFaultDOFs_TPV205(dof_data, num_fault_total, fault_coords);

      // R-013: --disable-nucleation overrides the patch-dependent τ_strike
      // with the uniform background value (70 MPa) everywhere on the
      // fault.  This kills the spontaneous-rupture trigger (the 81.6 MPa
      // nucleation patch can no longer exceed μ_s·σ_n = 81.24 MPa) and
      // is useful for verifying the strength-barrier semantics or
      // running a "what if there were no patches" diagnostic.  μ_s,
      // μ_d, d_c, and the strength-barrier mask are unaffected.
      if (disable_nucleation)
      {
         // R-002: also overwrite tau2_corr so the t=0 station row reports
         // the uniform background value the override establishes.  Pre-fix
         // tau2_corr was left at the patch value (set by
         // InitializeFaultDOFs_TPV205) so the first row of any patch-
         // overlapping station file leaked the original 81.6 / 78.0 / 62.0
         // MPa even though the user requested --disable-nucleation.  From
         // step 1 onward WriteBackState rebuilds tau2_corr from the new
         // tau2_0 = 70 MPa.
         for (auto &d : dof_data)
         {
            d.tau2_0    = TPV205Params::tau_back;
            d.tau2_corr = TPV205Params::tau_back;
         }
      }
   }
   // Print the disable-nucleation banner on rank 0 unconditionally
   // (rank 0 may have zero local fault QPs at higher rank counts; the
   // banner reflects the run-wide state).
   if (rank == 0 && disable_nucleation)
   {
      std::cout << "[tpv205_driver] --disable-nucleation: "
                << "overriding all 4 patches with the uniform "
                << "background τ_strike = " << TPV205Params::tau_back/1e6
                << " MPa.  Spontaneous rupture suppressed.\n";
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
         real_t dx = fault_coords[i](0) - TPV205Params::hypo_along_strike;
         real_t dz = std::abs(fault_coords[i](2)) - TPV205Params::hypo_down_dip;
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

   FaultFaceFlux fault_flux(TPV205Params::rho, TPV205Params::cp,
                            TPV205Params::cs);
   wave.SetFaultFlux(&fault_flux);
   // R-016: tell the wave operator to dispatch the LSW closed-form on
   // the fault branch (interior + R-1600 shared-fault fallback).
   // Without this, shared-fault QPs at np > 1 fall through to
   // fault_flux_->EvaluateADER (Brent on rate-and-state) which produces
   // wrong physics for LSW and stalls the rupture front at MPI rank
   // boundaries.  Must be called BEFORE the time loop.
   wave.SetFaultFrictionLaw(FaultFrictionLaw::LSW);
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
      // coincide because TPV205 builds its ParMesh on MPI_COMM_WORLD,
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

   // TPV205 uses linear slip-weakening (LSW), which has no state
   // variable — there is no AgingLawPsi / SlipLawSRWPsi instance.  The
   // sub-step iterator binds only to the FaultFaceFlux for its pure
   // helpers (`ComputeTrialTraction`, `BuildImposedState`,
   // `WriteBackState`) and does the LSW closed-form solve internally.
   mfem::seas::Tpv205SubStepIterator substep_iterator(fault_flux);

   // TPV205 ALWAYS routes through the substep iterator — the legacy
   // wave.AdvanceADER fault branch only supports rate-and-state friction
   // (Brent solver), not LSW.  --fault-iterator one-shot is honored as
   // O = 1 (single sub-step covering the whole macro-step); --fault-iterator
   // substep uses O = ader_order.  In both cases the dispatch path is
   // AdvanceADERWithSubStep (per-sub-step pointwise Q + LSW closed-form).
   const bool substep_quadrature =
      (GetDispatchedIterator(fault_iterator) == DispatchedIterator::SubStep);
   if (substep_quadrature)
   {
      // Re-configure with O = ader_order sub-steps (default).
      const int O = std::max(1, ader_order);
      std::vector<real_t> deltaT(O, 0.0), weights(O, 1.0 / static_cast<real_t>(O));
      for (int o = 0; o < O; o++) { deltaT[o] = dt / static_cast<real_t>(O); }
      substep_iterator.SetSubSteps(deltaT, weights);
   }
   else
   {
      // O = 1: single sub-step per macro-step (one-shot semantics).
      std::vector<real_t> deltaT(1, dt);
      std::vector<real_t> weights(1, 1.0);
      substep_iterator.SetSubSteps(deltaT, weights);
   }
   if (rank == 0)
   {
      std::cout << "[tpv205_driver] iterator quadrature: "
                << "O = " << substep_iterator.NumSubSteps()
                << " (ALL TPV205 runs use the iterator — LSW is not "
                << "supported in the legacy wave.AdvanceADER fault "
                << "branch).\n";
   }
   if (nprocs > 1 && rank == 0)
   {
      std::cout << "[tpv205_driver] MPI dispatch (nprocs=" << nprocs
                << "): shared-fault QPs use FaultFrictionLaw::LSW "
                << "dispatch in wave_operator (REVIEW R-016 fix; "
                << "shared-fault closed-form LSW solve via "
                << "FaultFaceFlux::EvaluateADER_LSW).\n";
   }
   // -----------------------------------------------------------------------
   // 6. Initialize Q = 0 (fluctuation-Q).
   // -----------------------------------------------------------------------
   Vector Q(NUM_STATE * ndof_total);
   InitializeState_TPV205(Q, ndof_total);

   MFEM_VERIFY(wave.GetAbsorbingBackground() != nullptr,
               "tpv205_driver: SetAbsorbingBackground(Q_bg=0) not called.");

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
                << TPV205Params::tau_back / 1e6 << " MPa, "
                << "σ_n = " << TPV205Params::sigma_n / 1e6 << " MPa\n";
      std::cout << "Patches:    nucleation 81.6 MPa @ (0, 7.5 km), "
                << "left 78.0 MPa @ (-7.5, 7.5 km), "
                << "right 62.0 MPa @ (+7.5, 7.5 km)  "
                << "(side 3 km, depth 7.5 km)\n";
      std::cout << "LSW:        μ_s = " << TPV205Params::mu_s
                << " (rupture area), " << TPV205Params::mu_s_barrier
                << " (strength barrier)  μ_d = " << TPV205Params::mu_d
                << "  d_c = " << TPV205Params::d_c << " m\n";
      std::cout << "Rupture area: |x| < " << TPV205Params::rupture_along_strike_half/1e3
                << " km, 0 ≤ depth ≤ " << TPV205Params::rupture_depth/1e3
                << " km, V_ini = " << TPV205Params::V_ini << " m/s\n\n";
   }

   // -----------------------------------------------------------------------
   // 7. Station output
   // -----------------------------------------------------------------------
   auto stations = DefaultStations_TPV205();
   TPV205StationWriter station_writer;
#ifdef MFEM_USE_MPI
   station_writer.Open(output_dir, output_prefix, stations,
                       fault_coords, num_fault_local, comm);
#else
   station_writer.Open(output_dir, output_prefix, stations,
                       fault_coords, num_fault_local);
#endif
   station_writer.WriteStep(0.0, dof_data);

   auto surface_stations = DefaultSurfaceStations_TPV205();
   if (rank == 0 && surface_stations.empty())
   {
      std::cout << "[tpv205_driver] surface stations: NONE configured "
                << "(R-006).  DefaultSurfaceStations_TPV205 returns an "
                << "empty list until the SCEC TPV5 §III.2 receiver "
                << "coordinates are transcribed; no off-fault traces "
                << "will be written this run.\n";
   }
   TPV205SurfaceStationWriter surface_writer;
   surface_writer.Open(output_dir, output_prefix, surface_stations,
                       pmesh, fes);
   surface_writer.WriteStep(0.0, Q);

   // -----------------------------------------------------------------------
   // 7b. ParaView output (mirrors tpv102_driver.cpp / BP5 seas::ParaViewOutput
   //     pattern).  Two collections:
   //       - pv_out (output_dir/ParaView): velocity + mpi_rank volume +
   //         fault-surface PVD/VTU (slip, slip_rate, traction dip+strike,
   //         mu_eff(δ), sigma_n, plus static μ_s, d_c, x2, x3) at the
   //         fault schedule.
   //       - pv_bulk_out (output_dir/ParaView_bulk): velocity + sigma_yy +
   //         sigma_xy + sigma_xz + mpi_rank at --paraview-bulk-dt cadence.
   //
   // R-801 / BP5 component convention enforced project-wide: comp 0 = dip,
   // comp 1 = strike.  TPV205 is pure strike-slip so the strike channel
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
      // Phase 6.3: select the volume back end.
      auto volume_mode =
         seas::ParaViewOutput<MeshT>::DefaultVolumeOutputMode();
      if (paraview_volume_force_vtu)
      { volume_mode = seas::ParaViewOutput<MeshT>::VolumeOutputMode::Vtu; }
      if (paraview_volume_force_hdf5)
      { volume_mode = seas::ParaViewOutput<MeshT>::VolumeOutputMode::Hdf5; }
      pv_out = std::make_unique<seas::ParaViewOutput<MeshT>>(
         output_dir + "/ParaView", pmesh, order,
         /*collection_name=*/"volume", volume_mode);

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
         // R-016: ParaView "a" / "Dc" channels carry LSW μ_s / d_c
         // (mirrors the station-writer convention for trace columns).
         // TPV104's fault output uses d.a = direct-effect parameter;
         // TPV205 reuses the same channel slot for LSW μ_s but reads
         // it from the dedicated lsw_mu_s field, not the repurposed
         // d.a slot (which is now zeroed at init for TPV205).
         pv_local_a(i)  = dof_data[i].lsw_mu_s;
         pv_local_Dc(i) = dof_data[i].lsw_d_c;
         pv_local_x2(i) = fault_coords[i](0);
         // R-014: write x3 as POSITIVE down-dip depth so ParaView
         // filters and the kStationsTPV205 station labels (which use
         // positive depth, e.g. `x2_0_x3_7.5`) reference the same
         // coordinate.  Mesh z is negative below the surface.
         pv_local_x3(i) = std::abs(fault_coords[i](2));
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

      // R-101 wiring: Phase 2b / 2d.3 / 3 CLI overrides.  See
      // tpv102_driver.cpp for the canonical doc-block.
      if (paraview_force_vtu)
      {
         pv_out->SetFaultOutputMode(
            seas::ParaViewOutput<MeshT>::FaultOutputMode::Vtu);
      }
      if (paraview_force_hdf5)
      {
         pv_out->SetFaultOutputMode(
            seas::ParaViewOutput<MeshT>::FaultOutputMode::Hdf5);
      }
      if (paraview_legacy_ascii)
      {
         pv_out->SetFaultOutputMode(
            seas::ParaViewOutput<MeshT>::FaultOutputMode::Vtu);
         pv_out->SetLegacyAsciiVTU(true);
      }
#ifdef MFEM_USE_HDF5
      if (paraview_fault_zfp_tol > 0.0)
      {
         pv_out->SetFaultHDFCompression(
            mfem::ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy,
            paraview_fault_zfp_tol);
      }
      else if (paraview_fault_deflate >= 0)
      {
         pv_out->SetFaultHDFCompression(
            mfem::ParaViewHDFDataCollection::HDFCompression::Deflate,
            static_cast<double>(paraview_fault_deflate));
      }
      // Phase 6.3a: --paraview-volume-* controls primary `pv_out`.
      // R-310: --paraview-bulk-* is RE-ROUTED to `pv_bulk_out` below.
      if (paraview_volume_zfp_tol > 0.0)
      {
         pv_out->SetVolumeHDFCompression(
            mfem::ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy,
            paraview_volume_zfp_tol);
      }
      else if (paraview_volume_deflate >= 0)
      {
         pv_out->SetVolumeHDFCompression(
            mfem::ParaViewHDFDataCollection::HDFCompression::Deflate,
            static_cast<double>(paraview_volume_deflate));
      }
#endif
      if (paraview_max_snapshots > 0)
      {
         pv_out->GetSchedule().max_total_snapshots = paraview_max_snapshots;
      }
      if (paraview_coseismic_dt    > 0.0)
      { pv_out->GetSchedule().dt_coseismic    = paraview_coseismic_dt; }
      if (paraview_nucleation_dt   > 0.0)
      { pv_out->GetSchedule().dt_nucleation   = paraview_nucleation_dt; }
      if (paraview_interseismic_dt > 0.0)
      { pv_out->GetSchedule().dt_interseismic = paraview_interseismic_dt; }
      pv_out->GetSchedule().Validate();
      pv_out->SetTotalRunTime(tfinal);

      // Phase 4: volume-PV decouple (see tpv102_driver.cpp for canonical comment).
      const bool volume_save_enabled = (volume_pv_dt > 0.0) || !pv_no_domain;
      pv_out->SetVolumeSaveEnabled(volume_save_enabled);
      if (volume_pv_dt > 0.0) { pv_out->SetVolumePVDt(volume_pv_dt); }

      if (rank == 0)
      {
         std::cout << "ParaView output: ON (prefix="
                   << output_dir << "/ParaView)\n";
         if (pv_no_domain && volume_pv_dt <= 0.0)
         {
            std::cout << "  Mode: fault-surface PVD only "
                         "(--no-volume-pv / --no-domain-pv)\n";
         }
         if (volume_pv_dt > 0.0)
         {
            std::cout << "  Volume cadence: every " << volume_pv_dt
                      << " s (--volume-pv-dt)\n";
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
         // Phase 6.4: secondary collection inherits the primary's
         // back end; distinct collection name "wave_bulk".
         pv_bulk_out = std::make_unique<seas::ParaViewOutput<MeshT>>(
            output_dir + "/ParaView_bulk", pmesh, order,
            /*collection_name=*/"wave_bulk", volume_mode);
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

         // Phase 6.4 R-310: --paraview-bulk-* now applies to the
         // secondary wavefield collection.
#ifdef MFEM_USE_HDF5
         if (paraview_bulk_zfp_tol > 0.0)
         {
            pv_bulk_out->SetVolumeHDFCompression(
               mfem::ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy,
               paraview_bulk_zfp_tol);
         }
         else if (paraview_bulk_deflate >= 0)
         {
            pv_bulk_out->SetVolumeHDFCompression(
               mfem::ParaViewHDFDataCollection::HDFCompression::Deflate,
               static_cast<double>(paraview_bulk_deflate));
         }
#endif

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

      // Phase 4: see tpv102_driver.cpp for canonical comment.
      const bool volume_active = pv_out->GetVolumeSaveEnabled();
      if (volume_active || bulk_wants)
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
         // R-005 / R-016: report μ_eff(δ) so the ParaView "state" channel
         // actually weakens during rupture and matches the station-writer
         // column 9.  Reads the LSW-native fields directly (no repurposing).
         {
            const real_t delta = std::sqrt(d.slip1 * d.slip1
                                           + d.slip2 * d.slip2);
            pv_local_state(i) = mfem::seas::LSWFrictionCoefficient_TPV205(
                                   delta,
                                   d.lsw_mu_s, d.lsw_mu_d, d.lsw_d_c);
         }
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

      // Phase 4: see tpv102_driver.cpp for canonical comment.
      if (pv_out->GetVolumeSaveEnabled())
      {
         pv_out->UpdateFaultFieldsBP5(pv_local_slip, pv_local_slip_rate,
                                      pv_local_traction, pv_local_state,
                                      pv_local_normal_stress);
         pv_out->ForceSave(step_num, time);
         // R-003: ForceSave only advances last_write_time_; the V_max-
         // adaptive schedule additionally needs current_regime_ /
         // last_v_max_ advanced for the next PeekShouldWrite to use the
         // correct regime interval (paraview_output.hpp:985-990).
         // Mirrors the TPV102/TPV104 R-006 fix.
         pv_out->CommitSchedule(time, V_max);
      }
      else
      {
         pv_out->CommitSchedule(time);
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
   paraview_write(0, 0.0, TPV205Params::V_ini);

   // If tfinal == 0: init-only run, stations at t=0 already written.
   // Skip the time loop and go straight to summary.  This is the
   // Phase-2 P2_D gate — the init sbatch verifies ψ_ini at every
   // station from the single t=0 row.
   if (nsteps == 0)
   {
      if (rank == 0)
      {
         std::cout << "[tpv205_driver] tfinal = 0 — init-only run. "
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

   if (rank == 0)
   {
      std::cout << "Starting ADER-O(" << ader_order << ") time loop...\n";
   }

   for (int step = 0; step < nsteps; ++step)
   {
      real_t dt_step = std::min(dt, tfinal - t);
      if (dt_step <= 0.0) { break; }

      // TPV205 nucleates by static patch pre-stress (set in
      // InitializeFaultDOFs_TPV205), so there is no per-step nucleation
      // increment — the τ_strike(x, z) field is already at its full
      // patch value from t = 0.  `disable_nucleation` is honored by
      // skipping the patch lookup at init time (driver path below).
      //
      // The substep iterator dispatch is mandatory for TPV205 interior
      // QPs: LSW friction is not implemented in the rate-and-state
      // EvaluateADER closure that the legacy `wave.AdvanceADER` fault
      // branch uses.  All TPV205 runs route through
      // `AdvanceADERWithSubStep`, which uses the iterator's pure
      // ComputeTrialTraction + LSW closed-form + BuildImposedState
      // pipeline and then hands the I_imp side-channel to the wave
      // operator's bulk corrector.  Shared-fault QPs at np > 1
      // additionally re-dispatch through FaultFaceFlux::EvaluateADER_LSW
      // via the R-1601 SHARED FALLBACK in wave_operator.inl.
      AdvanceADERWithSubStep(wave, substep_iterator, dof_data,
                             fault_coords, Q, dt_step,
                             ader_order, /*t_step_start=*/t,
                             Q_new);
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

      // LSW has no state variable — no ψ update.  Slip is accumulated
      // PER SUB-STEP inside Tpv205SubStepIterator (see step 5 of
      // StepOneQP_), so the driver does NOT re-integrate slip here.

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
            // R-001 (round 5): final column is the LSW effective
            // friction coefficient μ_eff(δ), not psi (TPV205 has no
            // state variable; d.psi is defensively zeroed at init and
            // never written).
            fs << "# columns: dof_idx x y z V1 V2 slip1 slip2 "
                  "tau1_corr tau2_corr sigma_n_corr mu_eff_delta\n";
            const int n = static_cast<int>(dof_data.size());
            for (int i = 0; i < n; i++)
            {
               const Vector &c = fault_coords[i];
               const DOFData &d = dof_data[i];
               const real_t delta = std::sqrt(d.slip1 * d.slip1
                                              + d.slip2 * d.slip2);
               const real_t mu_eff =
                  mfem::seas::LSWFrictionCoefficient_TPV205(
                     delta, d.lsw_mu_s, d.lsw_mu_d, d.lsw_d_c);
               fs << i << " " << c(0) << " " << c(1) << " " << c(2)
                  << " " << d.V1 << " " << d.V2
                  << " " << d.slip1 << " " << d.slip2
                  << " " << d.tau1_corr << " " << d.tau2_corr
                  << " " << d.sigma_n_corr << " " << mu_eff << "\n";
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
      std::cout << "TPV205 run complete.\n";
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

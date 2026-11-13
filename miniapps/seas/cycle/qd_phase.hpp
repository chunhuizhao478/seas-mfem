// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// qd_phase.hpp — fully dynamic BP5 cycle (plan rev 3, Section 3.3 d):
// QDPhase, adapter 1, the quasi-dynamic regime.
//
// It builds the QD stack of seas_bp5_full (tests/verification/
// bp5_verification_full.cpp:1549-1652, 1773-1802, 2481-2550) from deck
// values instead of BP5Params' analytic fields:
//
//   QDPhase(pmesh, mpi, cfg):
//     dom  = ElasticityDomainOperator(pmesh, p, lambda, mu, Vp, IP, solver,
//                                     FarField)       # Dirichlet {5}, natural {1}
//     seed = BP5Params from deck values (sigma_n, b, f0, V0, V_bg, Vp, rho, cs)
//     geom = FaultGeometry(dom, seed, mpi, /*compute_bp5_params=*/false)   # H1a
//     Table() = owned fault DOFs: XYZ, bulk element, attribute
//   Setup(par, pre):
//     geom.SetRateStatePerDOF(par.rs, (0, 1))                 # H1a: a, Dc, eta
//     geom.SetInitialFieldsPerDOF(pre.tau0, par.V_init, pre.sigma0)   # H1b
//     fault = RateStateFaultOperator(geom, friction, aging, seed)   # caches them
//     fault.SetScecPsiInit(ic == "scec_fd")
//     op    = SEASQuasiDynamicOperator(dom, fault), elastic sigma_n on
//     [cycle.qd].operators = "matrix": op.SetQDOperators(true)  (QD speed item 1)
//     [cycle.qd.green]: op.SetQDGreen(config, check_every)      (items 3, 4)
//     op.SetInitialCondition(state)
//     a = "amax" with ic = "scec_fd" (BP3-FD): psi(state) = pre.psi0, Solve(0)
//   InitTS(dt0): PETSc TS (rk 5dp from the options file), adaptive, dt0;
//     refuses a TS that is not RK with the FSAL property (seas_bp5_full's
//     check, QD speed-up v2 item A): the monitor reads the caches of the
//     last stage, which is the accepted state only for an FSAL tableau
//   Run(t, dt, t_final, max_steps, monitor): TSSolve; `monitor` runs at the
//     start step and after each accepted step, as seas_bp5_full's monitor,
//     and `false` stops the run at that step.  It runs from TS's post-step
//     hook, not from a TSMonitor: TSSolve calls the monitors before TSStep,
//     which resets the converged reason, so a stop set there was lost and
//     the run went on (found by the P3 smoke test; it also affected P1's
//     NaN stop).  Run calls it for the start step before TSSolve.
//   Solve(t): one right-hand side at (t, state) for the switch to the
//     dynamic regime (plan 6.2); Displacement, Slip, SlipRate, Psi, SigmaN
//     then hand the static state to the bridge (cycle/state_bridge.hpp).
//     With the Green's function a stage solves no displacement:
//     Displacement() solves it for the last stage when it is stale
//     (RefreshDisplacement), so ask for it only where it is used.
//   Velocity(v): the quasi-static velocity du/dt of the last evaluation.
//     K u = B s + D g(t) with g(t) = t g_hat (checked once), so du/dt
//     solves K du/dt = B V + D g_hat: the displacement at t = 1 with the
//     slip rate V as slip.  One solve (the face loops or the stage
//     operators, as the stages); the state and the caches are untouched.
//   Reenter(s, psi, t, step, dt0): the switch back (plan 6.3): the state
//     from the bridge, one Solve(t) (static u and V from the balance with
//     eta), then a new TS at t with first step dt0 whose step counter
//     continues at `step`.  The new TS starts like a fresh QD run from that
//     state (test_cycle_qd_restart_after_switch): the right-hand side is a
//     function of (t, state) alone (no warm start), so the history of this
//     object does not enter.  The monitor skips the re-entry step, as a
//     restart skips its checkpoint step: the rows at that time are the
//     dynamic regime's last, written before the switch.
//
// The kernel reads sigma_n, b, f0 and V0 as scalars outside SAFS mode
// (fault/rate_state_fault.hpp: sigma_n_bp5_, AgingLawPsi, psi_init), so
// Setup refuses per-DOF values that are not globally uniform and equal to
// the seed.  a, Dc and eta are per DOF.
//
// [stress] kind = "constant_tensor" (Chunhui Zhao, 2026-10-29):
// ProjectTensor gives the tensor's traction on the owned fault DOFs with the
// dynamic operator's own projection (domain/fault_stress_projection.hpp) for
// BuildTensorPrestress; Setup then takes sigma0 per DOF (no uniform
// sigma_n), switches the kernel to its per-DOF mode (SetSAFSMode: tau_pre and
// sigma_n read per DOF from FaultGeometry, which H1b filled), and, with ic =
// "scec_fd", writes BuildTensorPrestress's psi0 into the state as for a =
// "amax" (with "gold" the kernel inverts the same per-DOF balance itself).
// SigmaN adds the elastic normal traction to the DOF's own sigma0.
//
// Restart: a checkpoint written at step k, and a run restarted from it, both
// re-evaluate the first RK stage of step k + 1 (TSRestartStep) instead of
// reusing the FSAL stage, so the restarted run repeats the original run's
// steps bit for bit (test_cycle_restart).  A new TS evaluates its first
// stage anyway, so without the restart at checkpoint steps a restarted run
// would match the uninterrupted one only to round-off (seas_bp5_full: about
// 1e-11 relative after 25 steps).  The price: every checkpoint step moves the
// run by the same round-off against a run without checkpoints.
//
// [cycle.relaxation] (Chunhui Zhao, 2026-11-02; cycle/normal_stress_relaxation.hpp):
// Setup makes the relaxation, hands its reduction to the kernel
// (SetSigmaNReduction) and fills it at every stage time (the operator's
// friction stage hook); the post-step hook accepts each step (before the
// monitor, so the monitor, the switch and a checkpoint see r_{n+1}); Solve
// sets dsig from its evaluation; Reenter restarts the clock with r kept.
// Off (the default): no hook, no reduction, every bit unchanged.
//
// [cycle.opening] mode = "shear_free" (Chunhui Zhao, 2026-11-04: "start
// version A"; fault opening design rev 1.1, cycle/fault_opening.hpp): Setup
// makes the open set (OpenSetQD) on the Green's function, hands its mask to
// the kernel (SetFrozenNodes: open DOFs V = 0, dpsi/dt = 0) and to the
// relaxation (SetFrozen: r frozen), and its stage solve to the operator's
// Green's-function stage hook (tau = 0 at open DOFs).  SettleOpenSet_
// applies the rule at an evaluated state (write the solved open slip into
// the state, flip the DOFs whose sigma0 + T_n crossed 0, evaluate again,
// until nothing flips): after each accepted step (post-step hook, after the
// relaxation accepted the step, whose dsig is then reset to the new
// evaluation's, and before the monitor; a change restarts the RK stage,
// TSRestartStep, and the next step starts from [time].dt_initial, because
// the static step was not error-controlled: T5 on D9f, 0a2c611a, stepped
// 8,488 s into such a state and aborted), in Reenter with the dynamic
// regime's flags as the start, and at the start of a run (SettleOpening).
// Rank 0 logs every flip and every settle with flips ("[cycle.opening]
// flip", "settle", "next dt").  Off (the default): no mask, no hook, every
// bit unchanged.
//
// mode = "contact" (version B; fault opening plan rev 1, 2026-11-05,
// cycle/fault_opening.hpp file comment): the same machinery with the open
// unknowns grown to (delta, g) pairs — the stage solve imposes tau = 0 AND
// sigma = 0 on the open set, the rule opens at sigma < -sigma_tol and
// closes at solved g < 0 (hysteresis, D-B7), the flip log gains g, the
// stage operators are built with slip_comps = 3 for the unit g columns
// (QDUnitColumn; the dense G stays tangential, D-B6), the kernel sees the
// tapered floor behind [cycle.opening].floor_taper (plan eq. 6.2,
// RateStateFaultOperator::SetFloorTaper), and a friction solve without a
// finite V rejects the TS step instead of aborting (D-B5,
// SetNonFiniteVReject + TSSetFunctionDomainError in InitTS).
//
// [cycle.qd] exact_retry (Chunhui Zhao, 2026-11-07, the card "Safe retry
// only"; cycle/exact_retry.hpp): TSPreStep keeps a copy of the state each
// step starts from and TSPreStage copies it back at the first stage of a
// retry, so a rejected step never hands its retry PETSc's rollback, which
// turns a last-stage Inf rate into NaN in the state (Vista 1127480 at
// 12,222 yr: ten retries from psi = NaN, DIVERGED_STEP_REJECTED); the driver
// then scans the whole state for NaN after each step.  Off (the default):
// no hook, every bit unchanged.
//
// [cycle.qd] integrator = "esdirk" (option I; Chunhui Zhao, 2026-11-09:
// "let's proceed with implementation, take all decisions p1-p6 by
// recommended items"; cycle/qd_implicit.hpp, phase I3): the fault state is
// advanced by QDImplicitFault on its own TS (PETSc arkimex, fully implicit)
// and its own vector U = (in-plane slip, psi, gap) per owned DOF, built
// from the planar dense Green's function (PlanarProblem).  Everything else
// keeps reading this phase as before:
//   after each accepted step (PostStepImplicit_) U is written into State()
//   (slip at the in-plane slot, psi) and into the open set (flags: the
//   branch of the gap row, s = sigma - gamma g < 0; gap: g there, 0
//   elsewhere, and the same cleaned g back into U, so U is a function of
//   State() and the gap bit for bit and a restart repeats the run), then
//   one right-hand-side evaluation of the kernel at that state fills the
//   operator's caches (traction, sigma_n, slip rate, max V) for the monitor,
//   the outputs, the switch and the checkpoint.  In that evaluation the
//   open set only adds the carried gap's tractions (OpenSetQD::StageCarried:
//   no open-set solve, the stage solves fixed slip and gap), and no DOF is
//   frozen: at an open DOF sigma = 0, the tapered law has no strength, and
//   the kernel's friction solve is V = -tau / eta, the radiation row (P3).
//   MaxV() leaves the open DOFs out, as the DP5 contact path does (there
//   they are frozen): the switch and the hand-back guard see closed DOFs.
//   Start, switch back: the consistent start (plan 2.10; ConsistentStart:
//   the gap rows at fixed slip and psi) takes the place of the settle.
//   No relaxation (refused), no post-step settle, no reset of the step to
//   dt_initial after a flip (the flip happened inside an error-controlled
//   stage solve).
//   The checkpoint is unchanged (version 7 already carries the flags and
//   the gap); a checkpoint step restarts the stage and drops the LU factor
//   in the running run as in the restarted one.
// Off (the default "dp5"): none of this is built, every bit unchanged.
//
// This file and dyn_phase.hpp are the only cycle files that include solver
// headers (rule M4).

#ifndef MFEM_SEAS_CYCLE_QD_PHASE_HPP
#define MFEM_SEAS_CYCLE_QD_PHASE_HPP

#include "mfem.hpp"
#include "cycle_config.hpp"
#include "cycle_frames.hpp"
#include "loading.hpp"
#include "finite_width.hpp"
#include "fault_fields.hpp"
#include "prestress.hpp"
#include "fd_cycle_output.hpp"
#include "bench_output.hpp"
#include "normal_stress_relaxation.hpp"
#include "fault_opening.hpp"
#include "exact_retry.hpp"
#include "qd_implicit.hpp"
#include "qd_timers.hpp"

#include "../solver/seas_operator.hpp"
#include "../domain/elasticity_operator.hpp"
#include "../domain/fault_stress_projection.hpp"
#include "../fault/fault_geometry.hpp"
#include "../fault/rate_state_fault.hpp"
#include "../friction/dieterich_ruina.hpp"
#include "../friction/state_evolution.hpp"
#include "../config/bp5_params.hpp"
#include "../common/mpi_context.hpp"

#ifdef MFEM_USE_PETSC
#include "petsc.h"
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iomanip>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{
namespace cycle
{

class QDPhase
{
public:
   using DomainOp = ElasticityDomainOperator<ParMesh>;
   using FaultOp = RateStateFaultOperator<ParMesh, 2>;

   /// Called at the start step and after each accepted step with the step
   /// number and time.  Return false to stop the run after this call.
   using Monitor = std::function<bool(int step, real_t t)>;

   QDPhase(ParMesh &pmesh, MPIContext &mpi, const CycleConfig &cfg)
      : pmesh_(pmesh), mpi_(mpi), cfg_(cfg)
   {
      const spatial::SpatialFrictionConfig &s = cfg.shared;
      const spatial::RateStateBlock &rsb = *s.rate_state;
      const real_t lambda = s.material_fallback.lambda;
      const real_t mu = s.material_fallback.mu;
      const real_t rho = s.material_fallback.rho;

      // Seed: every scalar the kernel reads outside SAFS mode, from the deck.
      seed_.rho = rho;
      seed_.cs = std::sqrt(mu / rho);
      seed_.sigma_n = cfg.sigma_n;
      seed_.b = rsb.b_default;
      seed_.f0 = rsb.f_0_default;
      seed_.V0 = rsb.V_0_default;
      seed_.L0 = rsb.Dc_default;
      seed_.V_init = cfg.prestress.V_background;
      seed_.Vp = cfg.plate_rate;

      const SolverType st = (cfg.qd.solver == "mumps_blr") ? SolverType::MUMPS_BLR
                                                           : SolverType::MUMPS;
      // Legacy ctor: Dirichlet {5} (walls and creep plane), natural {1},
      // fault 3 (BCMode::FarField), exactly as seas_bp5_full builds it, with
      // the cycle's one fault orientation system (phase A): the operator
      // builds the FaultFrameField; the dynamic regime reads it.  Loading:
      // [cycle.loading] (cycle/loading.hpp; default the BP5 function).
      // The prescribed-slip faces must lie in the loading's band (loading.hpp).
      CheckCreepInBand(pmesh, {5}, cfg.loading);
      // [cycle.pml] (PML plan rev2, P1): the band, before the operator, so
      // the loading takes the D3 wall speed factor
      if (cfg.pml.On()) { SetupPmlBand_(pmesh); }
      FaultFrameSource frames;
      frames.config = std::make_shared<const FaultFrameConfig>(CycleFrameConfig(cfg));
      dom_ = std::make_unique<DomainOp>(pmesh, s.mesh.order, lambda, mu, cfg.plate_rate,
                                        seed_.Wf, seed_.lf, DGMethod::IP, st,
                                        BCMode::FarField, BasisType::GaussLobatto, frames,
                                        MakeCycleLoading(cfg, pml_wall_factor_));
      if (pml_band_ && cfg.pml.Cfs()) { dom_->SetPmlBand(pml_band_); }
      if (st == SolverType::MUMPS_BLR) { dom_->SetBLRTol(cfg.qd.blr_tol); }
      if (cfg.penalty_factor != 1.0) { dom_->SetPenaltyFactor(cfg.penalty_factor); }
      // [cycle.finite_width]: the 2.5D term, before the first K is used
      if (cfg.finite_width.On()) { SetupFiniteWidth_(); }

      geom_ = std::make_unique<FaultGeometry<ParMesh>>(*dom_, seed_, &mpi_,
                                                       /*compute_bp5_params=*/false);
      BuildTable_();
   }

   /// [cycle.pml]: the band (null when off) and the D3 wall speed factor
   /// (1 when off), for the dynamic phase (P2) to use the same ones.
   std::shared_ptr<const PmlBand> PmlBandPtr() const { return pml_band_; }
   real_t PmlWallSpeedFactor() const { return pml_wall_factor_; }

   /// The owned fault DOFs in the operator's owned order.
   const FaultDOFTable &Table() const { return table_; }
   /// [cycle.finite_width]: U_pl1, the substrate displacement at t = 1 on
   /// the operator's vdofs (cycle/finite_width.hpp), for DynPhase (one
   /// u_pl in both regimes); null when the term is off.
   const Vector *SubstrateShape() const
   { return cfg_.finite_width.On() ? &U_pl1_ : nullptr; }
   /// [stress] kind = "constant_tensor": the tensor's traction on the owned
   /// fault DOFs as kernel components (tau_1, tau_2, sigma_n), 3 per DOF in
   /// Table() order, by the dynamic operator's projection (file comment).
   /// Collective.
   void ProjectTensor(const std::array<real_t, 9> &S, Vector &K_own)
   {
      Vector K;
      ProjectStressToFaultKernel(*dom_, [S](const real_t[3], real_t s[9])
      { for (int i = 0; i < 9; i++) { s[i] = S[i]; } }, K);
      dom_->RestrictToOwnedFault(K, K_own, 3);
      MFEM_VERIFY(K_own.Size() == 3 * table_.Size(), "QDPhase::ProjectTensor: "
                  << K_own.Size() << " values for " << table_.Size() << " owned fault DOFs");
   }
   int NumGlobalFaultDOFs() const { return geom_->NumGlobalFaultDOFs(); }
   const BP5Params &Seed() const { return seed_; }
   DomainOp &Domain() { return *dom_; }
   FaultGeometry<ParMesh> &Geometry() { return *geom_; }

   /// Hooks H1a and H1b, then the kernel, the operator and the initial
   /// condition (one domain solve: assembly and factorisation).
   void Setup(const FaultParams &par, const Prestress &pre)
   {
      const int N = table_.Size();
      MFEM_VERIFY(par.rs.a.Size() == N && pre.tau0.Size() == 2 * N,
                  "QDPhase::Setup: parameters are not on this rank's owned fault DOFs");
      // The kernel's scalars (see the file comment); the tensor's sigma_n
      // is per DOF.
      const bool tensor = cfg_.TensorPrestress();
      if (!tensor) { AssertUniform_(par.rs.sigma_n_eff, seed_.sigma_n, "sigma_n"); }
      AssertUniform_(par.rs.b, seed_.b, "b");
      AssertUniform_(par.rs.f_0, seed_.f0, "f_0");
      AssertUniform_(par.rs.V_0, seed_.V0, "V_0");
      if (!tensor) { AssertUniform_(pre.sigma0, seed_.sigma_n, "sigma0"); }
      eta_ = GlobalUniform_(par.rs.eta, "eta");
      MFEM_VERIFY(pre.sigma0.Size() == N, "QDPhase::Setup: " << pre.sigma0.Size()
                  << " sigma0 values for " << N << " owned fault DOFs");
      sigma0_ = pre.sigma0;

      Vector init_vel_dir(2);
      init_vel_dir(0) = 0.0;   // dip
      init_vel_dir(1) = 1.0;   // strike
      geom_->SetRateStatePerDOF(par.rs, init_vel_dir);                  // H1a
      geom_->SetInitialFieldsPerDOF(pre.tau0, par.V_init, pre.sigma0);  // H1b

      DieterichRuinaFriction::Constants fc;
      fc.V0 = seed_.V0;
      fc.f0 = seed_.f0;
      fc.b = seed_.b;
      fc.Dc = seed_.L0;   // scalar fallback; the kernel uses per-DOF Dc
      friction_ = std::make_unique<DieterichRuinaFriction>(fc);
      aging_ = std::make_unique<AgingLawPsi>(seed_.b, seed_.V0, seed_.f0);
      fault_ = std::make_unique<FaultOp>(geom_.get(), friction_.get(), aging_.get(),
                                         seed_, &mpi_);
      fault_->SetScecPsiInit(!cfg_.IsGold());
      // The tensor: tau_pre and sigma_n per DOF, from H1b's copies.
      if (tensor) { fault_->SetSAFSMode(true, &geom_->GetTauPre(), &geom_->sigma_n_per_dof()); }
      fault_->SetSigmaNFloor(cfg_.sigma_n_floor);   // [cycle].sigma_n_floor_pa (0 = off)
      if (cfg_.opening.Contact())
      {
         // [cycle.opening] mode = "contact": the tapered floor (plan eq.
         // 6.2, behind [cycle.opening].floor_taper) and the non-finite step
         // rejection (D-B5; the TS domain check is set in InitTS)
         fault_->SetFloorTaper(cfg_.opening.floor_taper);
         fault_->SetNonFiniteVReject(true);
      }

      op_ = std::make_unique<PBP5SEASOp>(dom_.get(), fault_.get(), &mpi_);
      op_->SetElasticSigmaN(true);
      if (cfg_.Relaxation())
      {
         // [cycle.relaxation]: r = 0 at t = 0, so the start is unchanged
         relax_ = std::make_unique<NormalStressRelaxation>(N, cfg_.relax_t);
         fault_->SetSigmaNReduction(&relax_->Reduction());
         NormalStressRelaxation *rx = relax_.get();
         op_->SetFrictionStageHook([rx](real_t t) { rx->Fill(t); });
      }
      // QD speed items 1, 3, 4 (seas_bp5_full --qd-operators, --qd-green*):
      // after SetElasticSigmaN, which they build N_u, N_delta for.
      matrix_ = (cfg_.qd.operators == "matrix");
      if (matrix_)
      {
         // [cycle.opening] mode = "contact": the stage operators carry the
         // opening columns (slip_comps = 3) for the unit g columns of the
         // open-set solve (QDUnitColumn); the stage path itself and the
         // dense G stay tangential (D-B6)
         op_->SetQDOperators(true, cfg_.opening.Contact() ? 3 : 2);
         dom_->GetQDStageOperators()->PrintSummary(mpi_.GetComm(), std::cout);
      }
      if (cfg_.qd.green.enabled)
      {
         const QDGreenDeck &g = cfg_.qd.green;
         QDGreenConfig gc;
         gc.block = g.block;
         gc.cache_path = g.cache;
         gc.tol = g.tol;
         gc.compress = g.compress;
         gc.leaf = g.leaf;
         gc.eta = g.eta;
         gc.hcache_path = g.hcache;
         gc.compress_mem_gib = g.compress_mem_gib;
         op_->SetQDGreen(gc, g.check_every);
      }
      if (cfg_.opening.On()) { SetupOpening_(); }   // [cycle.opening]
#ifdef MFEM_USE_PETSC
      if (cfg_.qd.Esdirk()) { SetupImplicit_(par); }   // [cycle.qd] integrator = "esdirk"
#else
      MFEM_VERIFY(!cfg_.qd.Esdirk(), "[cycle.qd].integrator = \"esdirk\" needs PETSc");
#endif

      // R-003 of seas_bp5_full: allocate at least one entry so PETSc's
      // PlaceMemory/ResetMemory never sees an empty Memory on a rank without
      // fault DOFs; SetSize keeps the allocation when shrinking.
      const int actual = fault_->StateSize();
      state_.SetSize(std::max(actual, 1));
      state_.SetSize(actual);
      op_->SetInitialCondition(state_);
      // BP3-FD (a = "amax", ic = "scec_fd"; cycle/prestress.hpp): the
      // kernel's scec psi0 is BP5's uniform f0 + b ln(V0 / V_bg) and its
      // gold inversion subtracts eta |V_init|, so neither is BP3-FD's eq. 23.
      // The state takes BuildPrestress's psi0, and one evaluation at t = 0
      // leaves the slip rate and the traction consistent with it.
      if ((cfg_.prestress.AMax() || tensor) && !cfg_.IsGold())
      {
         for (int i = 0; i < N; i++)
         {
            state_(FaultOp::StatePerNode * i + FaultOp::PsiIndex) = pre.psi0(i);
         }
         Solve(0.0);
      }
   }

   Vector &State() { return state_; }
   /// Slip (kernel pairs, 2 N) and psi (N) into State(); nothing is solved.
   void SetFaultState(const Vector &slip_k, const Vector &psi)
   {
      const int N = fault_->NumNodes();
      MFEM_VERIFY(slip_k.Size() == 2 * N && psi.Size() == N,
                  "QDPhase::SetFaultState: " << slip_k.Size() << " slip and " << psi.Size()
                  << " psi values for " << N << " owned fault DOFs");
      for (int i = 0; i < N; i++)
      {
         state_(FaultOp::StatePerNode * i + FaultOp::SlipIndex) = slip_k(2 * i);
         state_(FaultOp::StatePerNode * i + FaultOp::SlipIndex + 1) = slip_k(2 * i + 1);
         state_(FaultOp::StatePerNode * i + FaultOp::PsiIndex) = psi(i);
      }
   }
   FaultOp &Fault() { return *fault_; }

   /// The opt-in sigma_n' floor ([cycle].sigma_n_floor_pa) since the last
   /// ResetFloorCount, over all ranks: the most RHS evaluations with a floored
   /// DOF on one rank, the most DOFs floored in one evaluation on one rank, and
   /// the lowest true sigma_n' below the floor with its position (min_sigma =
   /// inf: the floor never acted).  Collective.
   struct FloorInfo
   {
      long calls = 0;
      int max_nodes = 0;
      double min_sigma = std::numeric_limits<double>::infinity();
      double x[3] = {0.0, 0.0, 0.0};
   };
   FloorInfo FloorCount() const
   {
      const auto &f = fault_->GetFloorStats();
      FloorInfo r;
      long c = f.calls;
      int m = f.max_nodes;
      MPI_Allreduce(&c, &r.calls, 1, MPI_LONG, MPI_MAX, mpi_.GetComm());
      MPI_Allreduce(&m, &r.max_nodes, 1, MPI_INT, MPI_MAX, mpi_.GetComm());
      struct { double v; int rank; } in{static_cast<double>(f.min_sigma), mpi_.Rank()}, out;
      MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, mpi_.GetComm());
      r.min_sigma = out.v;
      if (std::isfinite(out.v))
      {
         if (mpi_.Rank() == out.rank)
         {
            for (int k = 0; k < 3; k++) { r.x[k] = table_.coords(3 * f.min_node + k); }
         }
         MPI_Bcast(r.x, 3, MPI_DOUBLE, out.rank, mpi_.GetComm());
      }
      return r;
   }
   void ResetFloorCount() { fault_->ResetFloorStats(); }
   PBP5SEASOp &Operator() { return *op_; }
   /// The global max |V| of the last evaluation.  Option I: over the closed
   /// DOFs (an open DOF slides at -tau / eta, which is not a rupture: as in
   /// the DP5 contact path, open DOFs stay out of max V and so out of the
   /// switch and the hand-back guard).  Collective.
   real_t MaxV() const
   {
#ifdef MFEM_USE_PETSC
      if (imp_)
      {
         const Vector &V = fault_->GetSlipRate();
         V.HostRead();
         real_t m = 0.0;
         for (int j = 0; j < fault_->NumNodes(); j++)
         {
            if (!open_->IsOpen(j)) { m = std::max(m, std::hypot(V(2 * j), V(2 * j + 1))); }
         }
         return mpi_.GlobalMax(m);
      }
#endif
      return op_->GetMaxSlipRate();
   }
   real_t Eta() const { return eta_; }

   /// Plan 6.2, step 1: one right-hand-side evaluation at (t, State()) (one
   /// solve), which leaves the static displacement, the traction, the
   /// normal traction and the slip rate consistent with the state, with the
   /// walls and the creep plane at Vp t.  Returns the global max V.
   real_t Solve(real_t t)
   {
      op_->SetTime(t);
      Vector rate(state_.Size());
      op_->Mult(state_, rate);
      if (relax_)
      {
         // dsig at t_n: every caller solves at the relaxation's time (the
         // start, the switch at an accepted step, the re-entry)
         MFEM_VERIFY(t == relax_->Time(), "QDPhase::Solve: t = " << t << " s is not the "
                     "relaxation's time " << relax_->Time() << " s");
         relax_->SetDelta(RelaxDelta_());
      }
      return MaxV();
   }

   /// [cycle.relaxation], or null when off.
   const NormalStressRelaxation *Relaxation() const { return relax_.get(); }

   /// [cycle.opening], or null when off.
   const OpenSetQD *Opening() const { return open_.get(); }
   /// [cycle.opening]: the open flags per owned DOF (1 open).
   const std::vector<char> &OpenFlags() const
   {
      MFEM_VERIFY(open_, "QDPhase::OpenFlags: [cycle.opening] is off");
      return open_->Flags();
   }
   /// [cycle.opening]: set the open flags (a restart: before the first
   /// stage).  Collective.
   void SetOpenFlags(const std::vector<char> &f)
   {
      MFEM_VERIFY(open_, "QDPhase::SetOpenFlags: [cycle.opening] is off");
      open_->SetFlags(f);
   }
   /// Contact: the gap per owned DOF [m] (0 at closed DOFs).
   void Gap(Vector &g) const
   {
      MFEM_VERIFY(open_ && open_->Contact(), "QDPhase::Gap: mode contact only");
      const std::vector<real_t> &gv = open_->Gap();
      g.SetSize(static_cast<int>(gv.size()));
      for (size_t j = 0; j < gv.size(); j++) { g(static_cast<int>(j)) = gv[j]; }
   }
   /// Contact: replace the gap (a restart: after SetOpenFlags, before the
   /// first stage; g != 0 only at open DOFs).  Collective.
   void SetGap(const Vector &g)
   {
      MFEM_VERIFY(open_ && open_->Contact(), "QDPhase::SetGap: mode contact only");
      std::vector<real_t> gv(static_cast<size_t>(g.Size()));
      for (int j = 0; j < g.Size(); j++) { gv[static_cast<size_t>(j)] = g(j); }
      open_->SetGap(gv);
   }
   /// Contact: the largest gap over all ranks [m].  Collective.
   real_t MaxGapGlobal() const
   {
      MFEM_VERIFY(open_ && open_->Contact(), "QDPhase::MaxGapGlobal: mode contact only");
      real_t m = 0.0;
      for (real_t g : open_->Gap()) { m = std::max(m, g); }
      return mpi_.GlobalMax(m);
   }
   /// [cycle.opening]: one Solve(t) and the rule at (t, State()) until
   /// nothing flips (the start of a run).  Returns the DOFs that flipped.
   /// Collective.
   long SettleOpening(real_t t)
   {
      MFEM_VERIFY(open_, "QDPhase::SettleOpening: [cycle.opening] is off");
#ifdef MFEM_USE_PETSC
      if (imp_)
      {
         const long flips = ConsistentStart_(t, "start at step", 0);
         Solve(t);
         return flips;
      }
#endif
      Solve(t);
      return SettleOpenSet_(t, state_, "start at step", 0);
   }
   /// [cycle.opening]: DOFs that flipped in the last Reenter's settle.
   long ReenterFlips() const { return reenter_flips_; }
   /// [cycle.opening]: largest |tau_pre + T| of the last evaluation at the
   /// open DOFs [Pa] over all ranks (0 when none: the stage solve's residual).
   real_t OpenShearResidual() const
   {
      real_t m = 0.0;
      if (open_)
      {
         const Vector &T = op_->GetTraction();
         const Vector &tp = fault_->TauPre();
         for (int j = 0; j < open_->NumOwned(); j++)
         {
            if (!open_->IsOpen(j)) { continue; }
            m = std::max(m, std::hypot(tp(2 * j) + T(2 * j), tp(2 * j + 1) + T(2 * j + 1)));
         }
      }
      return mpi_.GlobalMax(m);
   }
   /// A checkpoint's relaxation state (restart), before the first stage.
   void SetRelaxationState(real_t t_n, const Vector &r, const Vector &dsig)
   {
      MFEM_VERIFY(relax_, "QDPhase::SetRelaxationState: [cycle.relaxation] is off");
      relax_->SetState(t_n, r, dsig);
   }
   /// [cycle.relaxation] over all ranks: r range and the lowest sigma_n'
   /// the law sees before the floor (sigma0 + dsig - r at t_n) with its
   /// position.  Collective; every rank must call it (off: zeros, inf).
   struct RelaxInfo
   {
      double r_min = 0.0, r_max = 0.0;
      double min_sigma = std::numeric_limits<double>::infinity();
      double x[3] = {0.0, 0.0, 0.0};
   };
   RelaxInfo RelaxationSummary() const
   {
      RelaxInfo info;
      double lo = std::numeric_limits<double>::infinity(), hi = -lo, smin = lo;
      int imin = -1;
      if (relax_)
      {
         const Vector &r = relax_->R(), &d = relax_->Delta();
         for (int i = 0; i < r.Size(); i++)
         {
            lo = std::min(lo, static_cast<double>(r(i)));
            hi = std::max(hi, static_cast<double>(r(i)));
            const double sl = sigma0_(i) + d(i) - r(i);
            if (sl < smin) { smin = sl; imin = i; }
         }
      }
      MPI_Allreduce(&lo, &info.r_min, 1, MPI_DOUBLE, MPI_MIN, mpi_.GetComm());
      MPI_Allreduce(&hi, &info.r_max, 1, MPI_DOUBLE, MPI_MAX, mpi_.GetComm());
      struct { double v; int rank; } in{smin, mpi_.Rank()}, out;
      MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, mpi_.GetComm());
      info.min_sigma = out.v;
      if (std::isfinite(out.v))
      {
         if (mpi_.Rank() == out.rank)
         {
            for (int k = 0; k < 3; k++) { info.x[k] = table_.coords(3 * imin + k); }
         }
         MPI_Bcast(info.x, 3, MPI_DOUBLE, out.rank, mpi_.GetComm());
      }
      if (!relax_) { info.r_min = info.r_max = 0.0; }
      return info;
   }
   /// What the last evaluation left, for the bridge (owned fault DOFs; slip
   /// and slip rate as kernel (dip, strike) pairs, 2 N values).  The
   /// displacement is solved here when the Green's function left it stale
   /// (collective then; every rank must ask).  Contact solves on every call,
   /// with the gap as the third slip component (B5): the 2 N slip must not
   /// meet the widened stage operators, and the gap's elastic field belongs
   /// in the volume solve.
   const Vector &Displacement() const
   {
      if (open_ && open_->Contact())
      {
         Gap(gap_refresh_);
         op_->RefreshDisplacement(gap_refresh_);
      }
      else
      {
         op_->RefreshDisplacement();
      }
      return op_->GetDisplacement();
   }
   bool GreenActive() const { return op_->QDGreenActive(); }
   bool MatrixStages() const { return matrix_; }

   /// The quasi-static velocity du/dt of the last evaluation (see the file
   /// comment): a vector of the displacement space.  One solve; collective.
   void Velocity(Vector &v)
   {
      if (!loading_checked_) { CheckLinearLoading_(); }
      const Vector &V = fault_->GetSlipRate();
      if (!vel_gf_) { vel_gf_ = std::make_unique<ParGridFunction>(&dom_->GetFESpace()); }
      if (matrix_)
      {
         Vector T(V.Size());
         if (open_ && open_->Contact())
         {
            // contact (B6): the stage operators take 3 components per DOF,
            // so the tangential rates ride with a zero third component.  The
            // open set's gap rate is not tracked, so its du/dt share is left
            // out of this diagnostic velocity.
            const int n = V.Size() / 2;
            Vector V3(3 * n);
            V.HostRead();
            V3.HostWrite();
            for (int j = 0; j < n; j++)
            {
               V3(3 * j) = V(2 * j);
               V3(3 * j + 1) = V(2 * j + 1);
               V3(3 * j + 2) = 0.0;
            }
            dom_->SolveStage(1.0, V3, *vel_gf_, T, nullptr);
         }
         else
         {
            dom_->SolveStage(1.0, V, *vel_gf_, T, nullptr);
         }
      }
      else
      {
         Vector local;
         dom_->ExpandOwnedToLocalFault(V, local, dom_->NumSlipComponents());
         dom_->Solve(1.0, local, *vel_gf_);
      }
      v = *vel_gf_;
   }

   /// |V| on the owned fault DOFs from the kernel's slip rate of the last
   /// evaluation (the one max V comes from).
   void SlipRateMagnitude(Vector &s) const
   {
      const Vector &V = fault_->GetSlipRate();
      const int N = fault_->NumNodes();
      s.SetSize(N);
      for (int i = 0; i < N; i++) { s(i) = std::hypot(V(2 * i), V(2 * i + 1)); }
   }
   void Slip(Vector &s) const { fault_->GetSlip(state_, s); }
   const Vector &SlipRate() const { return fault_->GetSlipRate(); }
   void Psi(Vector &psi) const
   {
      const int N = fault_->NumNodes();
      psi.SetSize(N);
      for (int i = 0; i < N; i++) { psi(i) = state_(FaultOp::StatePerNode * i + FaultOp::PsiIndex); }
   }
   /// sigma_n_eff of the last friction solve (N values): sigma0 plus the
   /// elastic normal traction when that solve used it, as the kernel forms
   /// it (fault/rate_state_fault.hpp, ComputeRHS).
   void SigmaN(Vector &sn) const
   {
      const int N = fault_->NumNodes();
      const Vector *nt = op_->GetFrictionNormalTraction();
      sn.SetSize(N);
      for (int i = 0; i < N; i++) { sn(i) = nt ? sigma0_(i) + (*nt)(i) : sigma0_(i); }
   }
   /// The kernel's traction of the last evaluation (2 N, kernel frame) and
   /// the elastic normal traction (N, compression positive; empty unless
   /// elastic sigma_n was on in that evaluation).
   const Vector &Traction() const { return op_->GetTraction(); }
   const Vector &NormalTraction() const { return op_->GetNormalTraction(); }

   /// Station fields from the accepted state, as ParallelBP5BenchmarkOutput::
   /// Write forms them: V re-solved from (traction, state) and the normal
   /// traction the solver used, both from the last right-hand-side
   /// evaluation (station sigma_n fix ef70517 of gpu-qd-bp5).  sigma_n is
   /// SigmaN (the snapshots, and the station column of station_sigma_n).
   void FillStationFields(StationFields &f) const
   {
      const Vector &trac = op_->GetTraction();
      Vector slip, theta, V;
      fault_->GetSlip(state_, slip);
      fault_->GetTheta(state_, theta);
      fault_->RecomputeSlipRate(trac, state_, V, op_->GetFrictionNormalTraction());
      const Vector &tau_pre = geom_->GetTauPre();
      const int N = fault_->NumNodes();
      for (Vector *v : {&f.slip_dip, &f.slip_strike, &f.V_dip, &f.V_strike,
                        &f.tau_a_dip, &f.tau_a_strike, &f.tau_b_dip, &f.tau_b_strike})
      {
         v->SetSize(N);
      }
      for (int i = 0; i < N; i++)
      {
         f.slip_dip(i) = slip(2 * i);
         f.slip_strike(i) = slip(2 * i + 1);
         f.V_dip(i) = V(2 * i);
         f.V_strike(i) = V(2 * i + 1);
         f.tau_a_dip(i) = tau_pre(2 * i);
         f.tau_a_strike(i) = tau_pre(2 * i + 1);
         f.tau_b_dip(i) = trac(2 * i);
         f.tau_b_strike(i) = trac(2 * i + 1);
      }
      f.theta = theta;
      SigmaN(f.sigma_n);
      f.eta = eta_;
      // [cycle.relaxation]: the reduction the last evaluation used (r at
      // the accepted step, where the monitor reads it)
      if (relax_) { f.sigma_red = relax_->Reduction(); }
      else { f.sigma_red.SetSize(0); }
      OpenToStation_(f);
   }

   /// Stations on this rank's owned fault DOFs, in the plane frame of the
   /// deck: 2-D coordinates (x - o).e1 and (x - o).e2 for stations and DOFs,
   /// then LocateOnFaces (fd_cycle_output.hpp: face-basis match, else the
   /// nearest DOF, as seas_bp5_full's Probe2DInterpolator, but the face
   /// chosen where several meet is the same on any partition).
   std::vector<StationLocation> LocateStations(const OutputConfig &oc) const
   {
      std::vector<std::array<real_t, 2>> pts;
      for (const StationSpec &s : oc.stations)
      {
         pts.push_back({Project_(s.x.data(), oc.plane_origin, oc.plane_axis1),
                        Project_(s.x.data(), oc.plane_origin, oc.plane_axis2)});
      }
      return LocatePlanePoints(oc, pts);
   }

   /// Points given in plane coordinates (x2, x3), located as the stations
   /// (on a 2-D mesh the fault faces are segments; x3 is the distance off
   /// the fault line, 0 for a point on it).
   std::vector<StationLocation> LocatePlanePoints(
      const OutputConfig &oc, const std::vector<std::array<real_t, 2>> &pts) const
   {
      const int N = table_.Size();
      Vector x2(N), x3(N);
      for (int i = 0; i < N; i++)
      {
         const real_t *x = table_.coords.GetData() + 3 * i;
         x2(i) = Project_(x, oc.plane_origin, oc.plane_axis1);
         x3(i) = Project_(x, oc.plane_origin, oc.plane_axis2);
      }
      return LocateOnFaces(x2, x3, dom_->GetNbfPerFace(), pts, dom_->Dimension() - 1);
   }

   /// The fault data of BenchmarkOutput (bench_output.hpp) on the owned
   /// DOFs: positions, plane coordinates, the face quadrature of the
   /// traction (FaceQuadrature: basis at the points, w_q |J_q| and the
   /// points' plane coordinates per owned face) and the profile points.
   /// Each owned DOF must sit at its face's node (checked), so the owned
   /// DOFs are grouped by face as the face basis numbers them.
   BenchFaultGeometry BenchGeometry(const CycleConfig &c) const
   {
      const OutputConfig &oc = c.output;
      BenchFaultGeometry g;
      const int N = table_.Size();
      g.X = table_.coords;
      g.x2.SetSize(N);
      g.x3.SetSize(N);
      for (int i = 0; i < N; i++)
      {
         const real_t *x = table_.coords.GetData() + 3 * i;
         g.x2(i) = Project_(x, oc.plane_origin, oc.plane_axis1);
         g.x3(i) = Project_(x, oc.plane_origin, oc.plane_axis2);
      }
      const FaceQuadrature *fq = dom_->GetFaceQuadrature();
      MFEM_VERIFY(fq, "QDPhase::BenchGeometry: the domain has no fault face quadrature");
      const int nbf = dom_->GetNbfPerFace();
      const DenseMatrix &e_q = fq->BasisAtQuadPoints();
      const IntegrationRule &ir = fq->GetQuadRule();
      const IntegrationRule &nodes = fq->GetNodalRule();
      const int nq = ir.GetNPoints();
      MFEM_VERIFY(e_q.Height() == nbf && e_q.Width() == nq && nodes.GetNPoints() == nbf,
                  "QDPhase::BenchGeometry: face basis " << e_q.Height() << " x " << e_q.Width()
                  << ", " << nodes.GetNPoints() << " nodes, " << nbf << " DOFs per face");
      g.nbf = nbf;
      g.e_q = e_q;
      // Per local fault face (interior, then shared): w_q |J_q|, x2_q, x3_q,
      // and the face's nodes; copied to each of its DOFs and restricted to
      // owned like any fault field.
      const Array<int> &fint = dom_->GetFaultInteriorFaces();
      const Array<int> &fsh = dom_->GetFaultSharedFaces();
      const int nf = fint.Size() + fsh.Size();
      const int per = 3 * nq + 3;
      Vector loc(per * nbf * nf);
      // A 2-D mesh: the face is a segment, the normal and the points have 2
      // components, padded with z = 0 as the DOF table's positions.
      const int sdim = pmesh_.SpaceDimension();
      Vector nor(sdim), x(sdim);
      real_t xp[3] = {0.0, 0.0, 0.0};
      for (int f = 0; f < nf; f++)
      {
         FaceElementTransformations *T = (f < fint.Size())
            ? pmesh_.GetInteriorFaceTransformations(fint[f])
            : pmesh_.GetSharedFaceTransformations(fsh[f - fint.Size()]);
         MFEM_VERIFY(T, "QDPhase::BenchGeometry: fault face " << f << " has no transformation");
         std::vector<real_t> fd(3 * nq);
         for (int q = 0; q < nq; q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            T->SetAllIntPoints(&ip);
            CalcOrtho(T->Jacobian(), nor);
            T->Transform(ip, x);
            for (int d = 0; d < sdim; d++) { xp[d] = x(d); }
            fd[q] = ip.weight * nor.Norml2();
            fd[nq + q] = Project_(xp, oc.plane_origin, oc.plane_axis1);
            fd[2 * nq + q] = Project_(xp, oc.plane_origin, oc.plane_axis2);
         }
         for (int k = 0; k < nbf; k++)
         {
            real_t *r = &loc(per * (f * nbf + k));
            for (int i = 0; i < 3 * nq; i++) { r[i] = fd[i]; }
            const IntegrationPoint &np = nodes.IntPoint(k);
            T->SetAllIntPoints(&np);
            T->Transform(np, x);
            for (int d = 0; d < 3; d++) { r[3 * nq + d] = (d < sdim) ? x(d) : 0.0; }
         }
      }
      Vector own;
      dom_->RestrictToOwnedFault(loc, own, per);
      MFEM_VERIFY(own.Size() == per * N && N % nbf == 0, "QDPhase::BenchGeometry: "
                  << own.Size() / per << " owned values for " << N << " owned DOFs");
      const int nfo = N / nbf;
      g.wJ.SetSize(nfo * nq);
      g.q2.SetSize(nfo * nq);
      g.q3.SetSize(nfo * nq);
      real_t scale = 1.0, worst = 0.0;
      for (int i = 0; i < 3 * N; i++) { scale = std::max(scale, std::abs(table_.coords(i))); }
      for (int f = 0; f < nfo; f++)
      {
         const real_t *r0 = &own(per * f * nbf);
         for (int q = 0; q < nq; q++)
         {
            g.wJ(f * nq + q) = r0[q];
            g.q2(f * nq + q) = r0[nq + q];
            g.q3(f * nq + q) = r0[2 * nq + q];
         }
         for (int k = 0; k < nbf; k++)
         {
            const int j = f * nbf + k;
            const real_t *r = &own(per * j);
            for (int i = 0; i < 3 * nq; i++)
            {
               MFEM_VERIFY(r[i] == r0[i], "QDPhase::BenchGeometry: owned DOFs " << f * nbf
                           << " and " << j << " are not on one face (DOFs not grouped by face)");
            }
            for (int d = 0; d < 3; d++)
            {
               worst = std::max(worst, std::abs(r[3 * nq + d] - table_.coords(3 * j + d)));
            }
         }
      }
      worst = mpi_.GlobalMax(worst);
      MFEM_VERIFY(worst <= 1e-9 * scale, "QDPhase::BenchGeometry: an owned fault DOF is "
                  << worst << " m from its face node (face basis order != DOF order)");
      std::vector<std::array<real_t, 2>> pts;
      for (const ProfileSpec &p : c.bench.profiles)
      {
         for (int i = 0; i < p.points; i++)
         {
            const real_t s = static_cast<real_t>(i) / (p.points - 1);
            pts.push_back({p.from[0] + s * (p.to[0] - p.from[0]),
                           p.from[1] + s * (p.to[1] - p.from[1])});
         }
      }
      g.profile_locs = LocatePlanePoints(oc, pts);
      return g;
   }

#ifdef MFEM_USE_PETSC
   /// PETSc TS as seas_bp5_full --petsc-ts sets it up.  The options file
   /// (cfg.qd.petsc_options) must already be loaded by MFEMInitializePetsc.
   void InitTS(real_t dt0)
   {
      if (imp_) { InitImplicitTS_(dt0); return; }
      ode_ = std::make_unique<PetscODESolver>(mpi_.GetComm(), "");
      ode_->SetAbsTol(cfg_.qd.ts_atol);
      ode_->SetRelTol(cfg_.qd.ts_rtol);
      ode_->Init(*op_, PetscODESolver::ODE_SOLVER_GENERAL);
      TS ts = *ode_;
      TSAdapt adapt;
      MFEM_VERIFY(TSGetAdapt(ts, &adapt) == PETSC_SUCCESS, "TSGetAdapt failed");
      // PetscODESolver turns adaptivity off; RK45 error control needs it.
      MFEM_VERIFY(TSAdaptSetType(adapt, TSADAPTBASIC) == PETSC_SUCCESS,
                  "TSAdaptSetType(BASIC) failed");
      MFEM_VERIFY(TSSetFromOptions(ts) == PETSC_SUCCESS, "TSSetFromOptions failed");
      {
         // seas_bp5_full's check (QD speed-up v2 item A): the monitor reads
         // the operator's caches from the last stage of each step, the
         // accepted state only for an RK tableau with the FSAL property.
         TSType type = nullptr;
         MFEM_VERIFY(TSGetType(ts, &type) == PETSC_SUCCESS, "TSGetType failed");
         PetscBool is_rk = PETSC_FALSE, fsal = PETSC_FALSE;
         MFEM_VERIFY(PetscStrcmp(type, TSRK, &is_rk) == PETSC_SUCCESS, "PetscStrcmp failed");
         TSRKType rk = nullptr;
         if (is_rk)
         {
            MFEM_VERIFY(TSRKGetType(ts, &rk) == PETSC_SUCCESS, "TSRKGetType failed");
            MFEM_VERIFY(TSRKGetTableau(ts, nullptr, nullptr, nullptr, nullptr, nullptr,
                                       nullptr, nullptr, &fsal) == PETSC_SUCCESS,
                        "TSRKGetTableau failed");
         }
         MFEM_VERIFY(is_rk && fsal, "QDPhase: PETSc TS " << type << (is_rk ? " " : "")
                     << (is_rk ? rk : "") << " from " << cfg_.qd.petsc_options << ": the "
                     "monitor needs an RK tableau with the FSAL property (3bs, 5dp, 5bs, "
                     "6vr), because it reads the operator's caches from the last stage");
      }
      MFEM_VERIFY(TSSetTimeStep(ts, dt0) == PETSC_SUCCESS, "TSSetTimeStep failed");
      MFEM_VERIFY(TSSetApplicationContext(ts, this) == PETSC_SUCCESS,
                  "TSSetApplicationContext failed");
      MFEM_VERIFY(TSSetPostStep(ts, &QDPhase::PostStepCallback_) == PETSC_SUCCESS,
                  "TSSetPostStep failed");
      if (open_ && open_->Contact())
      {
         // [cycle.opening] mode = "contact" (D-B5): a non-finite stage state
         // (a friction solve that could not finish sent a non-finite rate
         // into the next stage, SetNonFiniteVReject) rejects the step, so
         // the controller shrinks dt and retries; the abort remains where
         // the adapt hits its minimum step (PETSc).
         MFEM_VERIFY(TSSetFunctionDomainError(ts, &QDPhase::StageStateFinite_)
                     == PETSC_SUCCESS, "TSSetFunctionDomainError failed");
      }
      if (cfg_.qd.exact_retry)
      {
         // [cycle.qd] exact_retry (cycle/exact_retry.hpp): the retry of a
         // rejected step starts from a copy of the step's start state, not
         // from PETSc's rollback.  Off (the default): no hook, every bit as
         // before.
         MFEM_VERIFY(TSSetPreStep(ts, &QDPhase::PreStepSave_) == PETSC_SUCCESS,
                     "TSSetPreStep failed");
         MFEM_VERIFY(TSSetPreStage(ts, &QDPhase::PreStageRestore_) == PETSC_SUCCESS,
                     "TSSetPreStage failed");
      }
   }

   /// [cycle.qd] exact_retry: TSPreStep, copy the state the step starts from.
   static PetscErrorCode PreStepSave_(TS ts)
   {
      PetscFunctionBeginUser;
      void *ctx = nullptr;
      PetscCall(TSGetApplicationContext(ts, &ctx));
      PetscCall(static_cast<QDPhase *>(ctx)->retry_.Save(ts));
      PetscFunctionReturn(PETSC_SUCCESS);
   }

   /// [cycle.qd] exact_retry: TSPreStage, at the first stage of a retry put
   /// the copy back over PETSc's rollback.  Rank 0 logs a retry whose
   /// rollback had left non-finite entries (the case that ended Vista
   /// 1127480).  Collective: every rank sees the same rejections.
   static PetscErrorCode PreStageRestore_(TS ts, PetscReal t)
   {
      PetscFunctionBeginUser;
      void *ctx = nullptr;
      PetscCall(TSGetApplicationContext(ts, &ctx));
      QDPhase *self = static_cast<QDPhase *>(ctx);
      bool restored = false;
      long bad = 0;
      PetscCall(self->retry_.RestoreIfRetry(ts, restored, &bad));
      if (!restored) { PetscFunctionReturn(PETSC_SUCCESS); }
      long bad_all = bad;
      MPI_Comm comm = PetscObjectComm(reinterpret_cast<PetscObject>(ts));
      PetscCallMPI(MPI_Allreduce(&bad, &bad_all, 1, MPI_LONG, MPI_SUM, comm));
      if (bad_all > 0)
      {
         self->retry_nonfinite_++;
         if (self->mpi_.Rank() == 0)
         {
            PetscInt step = 0;
            PetscReal dt = 0.0;
            PetscCall(TSGetStepNumber(ts, &step));
            PetscCall(TSGetTimeStep(ts, &dt));
            std::cout << "  [cycle.qd] exact retry: the rollback of a rejected try of step "
                      << step << " left " << bad_all << " non-finite state entries; "
                      << "the retry (dt = " << static_cast<double>(dt) << " s, t = "
                      << static_cast<double>(t) << " s) starts from the saved copy"
                      << std::endl;
         }
      }
      PetscFunctionReturn(PETSC_SUCCESS);
   }

   /// [cycle.opening] mode = "contact" (D-B5): PETSc calls this on every
   /// stage state; accept = false rejects the step.  The answer must agree
   /// on every rank (PETSc), hence the reduction.
   static PetscErrorCode StageStateFinite_(TS ts, PetscReal t, Vec Y,
                                           PetscBool *accept)
   {
      PetscFunctionBeginUser;
      PetscInt n = 0;
      const PetscScalar *a = nullptr;
      PetscCall(VecGetLocalSize(Y, &n));
      PetscCall(VecGetArrayRead(Y, &a));
      PetscInt ok = 1;
      for (PetscInt i = 0; i < n; i++)
      {
         if (!std::isfinite(static_cast<double>(PetscRealPart(a[i]))))
         {
            ok = 0;
            break;
         }
      }
      if (!ok)
      {
         // diagnosis (2026-11-07): name the entries (print only)
         void *ctx = nullptr;
         PetscCall(TSGetApplicationContext(ts, &ctx));
         if (ctx)
         {
            static_cast<const QDPhase *>(ctx)->ReportNonFinite_(
               a, n, "stage state", static_cast<real_t>(t));
         }
      }
      PetscCall(VecRestoreArrayRead(Y, &a));
      MPI_Comm comm = PetscObjectComm(reinterpret_cast<PetscObject>(ts));
      PetscCall(MPIU_Allreduce(MPI_IN_PLACE, &ok, 1, MPIU_INT, MPI_MIN, comm));
      *accept = ok ? PETSC_TRUE : PETSC_FALSE;
      if (!ok)
      {
         PetscMPIInt rank = 0;
         MPI_Comm_rank(comm, &rank);
         if (rank == 0)
         {
            std::cout << "  [cycle.opening] non-finite stage state at t = "
                      << static_cast<double>(t)
                      << " s: the step is rejected and dt shrinks (D-B5)\n";
            PetscReal dt_try = 0.0;
            PetscCall(TSGetTimeStep(ts, &dt_try));
            std::cout << "  [cycle.opening] diag: the rejected step tried dt = "
                      << static_cast<double>(dt_try) << " s" << std::endl;
         }
      }
      PetscFunctionReturn(PETSC_SUCCESS);
   }

   /// [cycle.opening] mode = "contact", diagnosis of D-B5 (2026-11-07): on
   /// this rank, print up to 5 non-finite entries of a TS vector (the
   /// fault state: per owned DOF the dip slip, the strike slip, psi) with
   /// the DOF's coordinates.  At most 200 lines per rank per run; print
   /// only.
   void ReportNonFinite_(const PetscScalar *a, PetscInt n, const char *what,
                         real_t t) const
   {
      const int per = FaultOp::StatePerNode;
      int shown = 0;
      for (PetscInt k = 0; k < n && shown < 5 && nonfinite_lines_ < 200; k++)
      {
         const double v = static_cast<double>(PetscRealPart(a[k]));
         if (std::isfinite(v)) { continue; }
         const int i = static_cast<int>(k) / per;
         const int c = static_cast<int>(k) % per;
         std::cout << std::scientific << std::setprecision(9)
                   << "  [cycle.opening] diag non-finite " << what << " at t = "
                   << t << " s: rank " << mpi_.Rank() << " DOF " << i << " x "
                   << geom_->GetCoordsX2()(i) << " y " << geom_->GetCoordsX3()(i)
                   << " component " << c
                   << (c == FaultOp::PsiIndex ? " (psi)" : " (slip)")
                   << " value " << v << std::endl;
         std::cout << std::defaultfloat;
         shown++;
         nonfinite_lines_++;
      }
   }

   /// Plan 6.3, the switch back: SetFaultState, Solve(t), and a new TS at t
   /// (first step dt0, step counter `step`).  Returns the global max V of
   /// the solve.  Contact (B4): `gap` is the dynamic regime's gap per owned
   /// DOF (g = d . n through the bridge); the carried set starts the first
   /// solve at it, and the settle's rule takes over.  Collective.
   real_t Reenter(const Vector &slip_k, const Vector &psi, real_t t, int step, real_t dt0,
                  const std::vector<char> *open_flags = nullptr,
                  const Vector *gap = nullptr)
   {
      SetFaultState(slip_k, psi);
      if (relax_) { relax_->Restart(t); }   // r kept from the switch; the clock at t
      if (open_)
      {
         // [cycle.opening]: the dynamic regime's open DOFs start the set;
         // the solve releases their shear (a static step), then the rule
         MFEM_VERIFY(open_flags, "QDPhase::Reenter: [cycle.opening] needs the dynamic "
                     "regime's open flags");
         open_->SetFlags(*open_flags);
         if (open_->Contact())
         {
            MFEM_VERIFY(gap, "QDPhase::Reenter: mode contact needs the dynamic "
                        "regime's gap");
            // CloseDOF's projection returns a closing DOF's gap to zero only
            // to rounding, and the bridge's dot products add their own, so
            // residue far below any physical opening is dropped; more than
            // that at a closed DOF (or a negative gap) is a real mismatch
            // between the flags and the state, and OpenSetQD::SetGap aborts.
            const real_t residue = 1e-9;   // m
            Vector gk(*gap);
            for (int j = 0; j < gk.Size(); j++)
            {
               if (std::abs(gk(j)) <= residue &&
                   (!(*open_flags)[static_cast<size_t>(j)] || gk(j) < 0.0))
               {
                  gk(j) = 0.0;
               }
            }
            SetGap(gk);
         }
      }
      real_t Vmax = 0.0;
      if (imp_)
      {
         // option I: the gap rows at the handed-back slip and psi (plan
         // 2.10), then the caches of that state
         reenter_flips_ = ConsistentStart_(t, "switch back at step", step);
         Vmax = Solve(t);
      }
      else
      {
         Vmax = Solve(t);
         if (open_)
         {
            reenter_flips_ = SettleOpenSet_(t, state_, "switch back at step", step);
            Vmax = op_->GetMaxSlipRate();
         }
      }
      InitTS(dt0);
      TS ts = Ts_();
      MFEM_VERIFY(TSSetStepNumber(ts, step) == PETSC_SUCCESS, "TSSetStepNumber failed");
      skip_monitor_step_ = step;   // that time's rows are the dynamic regime's
      return Vmax;
   }

   /// Continue from a checkpoint: the TS step counter and a fresh first stage.
   void RestoreTS(int step)
   {
      TS ts = Ts_();
      MFEM_VERIFY(TSSetStepNumber(ts, step) == PETSC_SUCCESS, "TSSetStepNumber failed");
      u_stale_ = true;   // option I: U from the restored State() and gap (Run)
      RestartStage();
      skip_monitor_step_ = step;   // its output is already in the files
   }

   /// The next step re-evaluates its first RK stage (see the file comment).
   void RestartStage()
   {
      TS ts = Ts_();
      MFEM_VERIFY(TSRestartStep(ts) == PETSC_SUCCESS, "TSRestartStep failed");
      // option I: a restarted run has no LU factor either
      if (imp_) { imp_->DropFactor(); }
   }

   /// Steps until t_final, max_steps, or `mon` returns false.  t and dt
   /// return the final time and the next step size.
   void Run(real_t &t, real_t &dt, real_t t_final, int max_steps, const Monitor &mon)
   {
      monitor_ = mon;
      stopped_ = false;
      // The start step (TSSolve's first TSMonitor call in seas_bp5_full).
      TS ts = Ts_();
      PetscInt s0 = 0;
      MFEM_VERIFY(TSGetStepNumber(ts, &s0) == PETSC_SUCCESS, "TSGetStepNumber failed");
      MFEM_VERIFY(TSSetTime(ts, t) == PETSC_SUCCESS && TSSetTimeStep(ts, dt) == PETSC_SUCCESS,
                  "TSSetTime/TSSetTimeStep failed");
      if (imp_)
      {
         // option I: U from State() and the gap when either was set from
         // outside (a new TS, a restart), then PetscODESolver::Run's calls
         if (u_stale_) { LoadU_(); }
         if (imp_phase_new_)
         {
            imp_mark_t0_ = imp_t_last_ = t;
            imp_phase_new_ = false;
         }
         const double wall0 = MPI_Wtime();
         if (CallMonitor_(static_cast<int>(s0), t))
         {
            MFEM_VERIFY(TSSetMaxSteps(ts, max_steps) == PETSC_SUCCESS &&
                        TSSetMaxTime(ts, t_final) == PETSC_SUCCESS &&
                        TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP) == PETSC_SUCCESS,
                        "TSSetMaxSteps/TSSetMaxTime failed");
            MFEM_VERIFY(TSSolve(ts, U_imp_) == PETSC_SUCCESS, "TSSolve failed");
            PetscReal pt = 0.0;
            MFEM_VERIFY(TSGetTime(ts, &pt) == PETSC_SUCCESS, "TSGetTime failed");
            t = static_cast<real_t>(pt);
            MFEM_VERIFY(TSGetTimeStep(ts, &pt) == PETSC_SUCCESS, "TSGetTimeStep failed");
            dt = static_cast<real_t>(pt);
         }
         imp_wall_ += MPI_Wtime() - wall0;
         monitor_ = nullptr;
         return;
      }
      if (CallMonitor_(static_cast<int>(s0), t))
      {
         ode_->SetMaxIter(max_steps);
         ode_->Run(state_, t, dt, t_final);
      }
      monitor_ = nullptr;
   }

   bool StoppedByMonitor() const { return stopped_; }

   real_t NextDt() const
   {
      PetscReal dt;
      TSGetTimeStep(Ts_(), &dt);
      return dt;
   }
   int StepNumber() const
   {
      PetscInt n;
      TSGetStepNumber(Ts_(), &n);
      return static_cast<int>(n);
   }
   int Rejections() const
   {
      PetscInt n;
      TSGetStepRejections(Ts_(), &n);
      return static_cast<int>(n);
   }
   /// [cycle.qd] exact_retry: retries started from the saved copy, and
   /// those whose PETSc rollback had left non-finite entries (this phase's
   /// lifetime, all TS objects).
   long ExactRetries() const { return retry_.Restores(); }
   long ExactRetriesNonFinite() const { return retry_nonfinite_; }
   TS GetTS() const { return Ts_(); }
   /// [cycle.qd] integrator = "esdirk": the implicit fault system (its
   /// counts), or null with "dp5".
   const QDImplicitFault *Implicit() const { return imp_.get(); }
   /// Option I, diagnostics (phase I5): what the implicit solve did since
   /// the last call (or the start): accepted steps, simulated time, the
   /// solver's counts, and the wall time of the QD loop (the monitor's
   /// output included).  The counts agree on every rank; the wall time is
   /// this rank's.
   struct ImplicitPhaseInfo
   {
      long steps = 0;
      real_t t0 = 0.0, t1 = 0.0;    ///< [s]
      long newton_its = 0, factors = 0, failed_solves = 0, columns = 0;
      long max_open = 0;            ///< most open points in one solve (whole run)
      double wall = 0.0;            ///< [s]
   };
   ImplicitPhaseInfo ImplicitPhase()
   {
      ImplicitPhaseInfo r;
      if (!imp_) { return r; }
      const QDImplicitFault::Counts &c = imp_->GetCounts();
      r.steps = imp_steps_ - imp_mark_steps_;
      r.t0 = imp_mark_t0_;
      r.t1 = imp_t_last_;
      r.newton_its = c.newton_its - imp_mark_.newton_its;
      r.factors = c.factors - imp_mark_.factors;
      r.failed_solves = c.newton_fails - imp_mark_.newton_fails;
      r.columns = c.columns - imp_mark_.columns;
      r.max_open = c.max_open;
      r.wall = imp_wall_;
      imp_mark_ = c;
      imp_mark_steps_ = imp_steps_;
      imp_wall_ = 0.0;
      imp_phase_new_ = true;
      return r;
   }
   /// Option I: failed stage solves of this phase's current TS (each one a
   /// rejected step).
   long ImplicitSNESFailures() const
   {
      PetscInt n = 0;
      if (imp_ && ts_imp_) { TSGetSNESFailures(ts_imp_, &n); }
      return static_cast<long>(n);
   }
#endif

   ~QDPhase()
   {
#ifdef MFEM_USE_PETSC
      if (ts_imp_) { TSDestroy(&ts_imp_); }
      if (U_imp_) { VecDestroy(&U_imp_); }
#endif
   }

private:
   /// [cycle.opening]: the open set on the Green's function (file comment).
   void SetupOpening_()
   {
      MFEM_VERIFY(op_->QDGreenActive(), "[cycle.opening] needs [cycle.qd.green] "
                  "(decision D6)");
      const QDGreenOperator *G = dom_->GetQDGreen();
      MFEM_VERIFY(G && pmesh_.Dimension() == 2 && G->Planar(), "[cycle.opening]: a 2-D "
                  "mesh with the planar Green's function only (design rev 1.1, Section 6)");
      DomainOp *dom = dom_.get();
      if (cfg_.opening.Contact())
      {
         // mode = "contact" (fault opening plan rev 1, Section 5): the open
         // unknowns are (delta, g) pairs; the delta columns stay the Green's
         // function's, the g columns are pure solve-path responses of a unit
         // opening through the stage operators' opening columns, fetched
         // once per DOF at its first opening (D-B6)
         MFEM_VERIFY(dom_->GetQDStageOperators() &&
                     dom_->GetQDStageOperators()->slip_comps == 3,
                     "[cycle.opening] mode = \"contact\": the stage operators "
                     "carry no opening columns (Setup builds them with "
                     "slip_comps = 3)");
         open_ = std::make_unique<OpenSetQD>(
            mpi_.GetComm(), fault_->NumNodes(), G->PlanarSlot(),
            [dom](const Vector &e, Vector &T, Vector &Tn) { dom->GreenStage(0.0, e, T, &Tn); },
            fault_->TauPre(),
            [dom](const Vector &e, Vector &T, Vector &Tn) { dom->QDUnitColumn(e, T, &Tn); },
            sigma0_, cfg_.opening.sigma_tol);
      }
      else
      {
         open_ = std::make_unique<OpenSetQD>(
            mpi_.GetComm(), fault_->NumNodes(), G->PlanarSlot(),
            [dom](const Vector &e, Vector &T, Vector &Tn) { dom->GreenStage(0.0, e, T, &Tn); },
            fault_->TauPre());
      }
      OpenSetQD *os = open_.get();
      if (cfg_.qd.Esdirk())
      {
         // option I (file comment): the stage solves fixed slip and gap, so
         // an evaluation adds the carried gap's tractions only, and no DOF
         // is frozen (an open DOF slides at V = -tau / eta, P3)
         op_->SetGreenStageHook([os](real_t, Vector &, Vector &T, Vector *Tn)
         { os->StageCarried(T, Tn); });
         return;
      }
      fault_->SetFrozenNodes(&open_->Mask());
      if (relax_) { relax_->SetFrozen(&open_->Mask()); }
      op_->SetGreenStageHook([os](real_t, Vector &slip, Vector &T, Vector *Tn)
      { os->Stage(slip, T, Tn); });
   }

#ifdef MFEM_USE_PETSC
   /// The TS in use: option I's own, or PetscODESolver's.
   TS Ts_() const { return imp_ ? ts_imp_ : static_cast<TS>(*ode_); }

   /// [cycle.qd] integrator = "esdirk": the implicit fault system on the
   /// planar dense Green's function (cycle/qd_implicit.hpp, PlanarProblem;
   /// the deck checks are CycleConfig's, the Green's function's are
   /// PlanarProblem's).  Collective.
   void SetupImplicit_(const FaultParams &par)
   {
      MFEM_VERIFY(open_ && open_->Contact() && op_->QDGreenActive(), "[cycle.qd].integrator "
                  "= \"esdirk\" needs [cycle.opening] mode = \"contact\" and "
                  "[cycle.qd.green]");
      MFEM_VERIFY(pmesh_.Dimension() == 2, "[cycle.qd].integrator = \"esdirk\": a 2-D "
                  "mesh only (option I, decision P5)");
      MFEM_VERIFY(!relax_, "[cycle.qd].integrator = \"esdirk\" with the relaxation");
      const QDGreenOperator *G = dom_->GetQDGreen();
      MFEM_VERIFY(G, "[cycle.qd].integrator = \"esdirk\": no Green's function");
      PlanarFaultInputs in;
      in.tau_pre = fault_->TauPre();
      in.sigma0 = sigma0_;
      in.a = par.rs.a;
      in.eta = par.rs.eta;
      in.Dc = par.rs.Dc;
      in.b = seed_.b;
      in.f0 = seed_.f0;
      in.V0 = seed_.V0;
      in.sigma_floor = cfg_.sigma_n_floor;
      in.floor_taper = cfg_.opening.floor_taper;
      QDImplicitSettings set;
      set.esdirk_type = cfg_.qd.esdirk_type;
      set.newton_tol = cfg_.qd.newton_tol;
      set.newton_max_it = cfg_.qd.newton_max_it;
      set.jacobian_reuse = cfg_.qd.jacobian_reuse;
      set.newton_local = cfg_.qd.newton_local;
      set.dense_lu = cfg_.qd.dense_lu;
      set.stage_guess = cfg_.qd.stage_guess;
      set.atol = cfg_.qd.ts_atol;
      set.rtol = cfg_.qd.ts_rtol;
      DomainOp *dom = dom_.get();
      imp_ = std::make_unique<QDImplicitFault>(
         mpi_.GetComm(),
         PlanarProblem(mpi_.GetComm(), *G, in,
                       [dom](const Vector &e, Vector &T, Vector &Tn)
      { dom->QDUnitColumn(e, T, &Tn); }), set);
      MFEM_VERIFY(imp_->NumOwned() == fault_->NumNodes(), "option I: " << imp_->NumOwned()
                  << " points for " << fault_->NumNodes() << " owned fault DOFs");
      if (mpi_.Rank() == 0)
      {
         const double mib = 3.0 * 8.0 * imp_->NumGlobal() * static_cast<double>(imp_->NumGlobal())
                            / (1024.0 * 1024.0);
         std::cout << "  [cycle.qd] integrator esdirk (option I): arkimex " << set.esdirk_type
                   << " fully implicit on " << imp_->NumGlobal() << " fault points, dense LU "
                   << (set.dense_lu == "shared" ? "shared by all ranks (ScaLAPACK), G "
                       : "") << "on every rank (about " << std::setprecision(4) << mib << " MiB per rank "
                   "for G and the factor), newton_tol " << set.newton_tol << ", newton_max_it "
                   << set.newton_max_it << ", jacobian_reuse " << set.jacobian_reuse
                   << ", newton_local " << (set.newton_local ? "on" : "off")
                   << (set.stage_guess == "extrapolate"
                       ? ", stage_guess extrapolate (Newton starts from the last "
                         "step's dense output)" : "") << "\n";
      }
   }

   /// Option I: a new TS at the current State() and gap (first step dt0).
   /// The options file's step-controller options apply (-ts_adapt_*); its
   /// -ts_type is replaced by arkimex.
   void InitImplicitTS_(real_t dt0)
   {
      MPI_Comm comm = mpi_.GetComm();
      if (ts_imp_) { MFEM_VERIFY(TSDestroy(&ts_imp_) == PETSC_SUCCESS, "TSDestroy failed"); }
      TS ts = nullptr;
      MFEM_VERIFY(TSCreate(comm, &ts) == PETSC_SUCCESS, "TSCreate failed");
      ts_imp_ = ts;
      TSAdapt adapt;
      MFEM_VERIFY(TSGetAdapt(ts, &adapt) == PETSC_SUCCESS, "TSGetAdapt failed");
      MFEM_VERIFY(TSAdaptSetType(adapt, TSADAPTBASIC) == PETSC_SUCCESS,
                  "TSAdaptSetType(BASIC) failed");
      MFEM_VERIFY(TSSetFromOptions(ts) == PETSC_SUCCESS, "TSSetFromOptions failed");
      MFEM_VERIFY(imp_->SetupTS(ts) == PETSC_SUCCESS, "QDImplicitFault::SetupTS failed");
      MFEM_VERIFY(imp_->SetTolerances(ts) == PETSC_SUCCESS,
                  "QDImplicitFault::SetTolerances failed");
      // a failed stage solve or error test shrinks the step; the run ends
      // where the controller reaches its smallest step (PETSc)
      MFEM_VERIFY(TSSetMaxStepRejections(ts, -1) == PETSC_SUCCESS,
                  "TSSetMaxStepRejections failed");
      MFEM_VERIFY(TSSetTimeStep(ts, dt0) == PETSC_SUCCESS, "TSSetTimeStep failed");
      MFEM_VERIFY(TSSetApplicationContext(ts, this) == PETSC_SUCCESS,
                  "TSSetApplicationContext failed");
      MFEM_VERIFY(TSSetPostStep(ts, &QDPhase::PostStepImplicit_) == PETSC_SUCCESS,
                  "TSSetPostStep failed");
      if (cfg_.qd.exact_retry)
      {
         MFEM_VERIFY(TSSetPreStep(ts, &QDPhase::PreStepSave_) == PETSC_SUCCESS,
                     "TSSetPreStep failed");
         MFEM_VERIFY(TSSetPreStage(ts, &QDPhase::PreStageRestore_) == PETSC_SUCCESS,
                     "TSSetPreStage failed");
      }
      if (!U_imp_)
      {
         MFEM_VERIFY(VecCreateMPI(comm, QDImplicitFault::kPer * imp_->NumOwned(),
                                  PETSC_DETERMINE, &U_imp_) == PETSC_SUCCESS,
                     "VecCreateMPI failed");
      }
      MFEM_VERIFY(TSSetSolution(ts, U_imp_) == PETSC_SUCCESS, "TSSetSolution failed");
      imp_->DropFactor();
      u_stale_ = true;
   }

   /// Option I: per owned DOF (in-plane slip, psi, gap) from State() and
   /// the open set's gap into @a u (3 N).
   void StateToU_(real_t *u) const
   {
      const int N = fault_->NumNodes();
      const int spn = FaultOp::StatePerNode, slot = open_->Slot();
      const std::vector<real_t> &g = open_->Gap();
      state_.HostRead();
      for (int o = 0; o < N; o++)
      {
         u[3 * o] = state_(spn * o + FaultOp::SlipIndex + slot);
         u[3 * o + 1] = state_(spn * o + FaultOp::PsiIndex);
         u[3 * o + 2] = g[static_cast<size_t>(o)];
      }
   }

   /// Option I: the gap columns of every point with g != 0 in @a u (a
   /// restart, a switch back), in global order.  Collective.
   void EnsureColumns_(const real_t *u)
   {
      const int N = fault_->NumNodes(), NG = imp_->NumGlobal();
      std::vector<int> has(static_cast<size_t>(NG), 0);
      for (int o = 0; o < N; o++)
      {
         if (u[3 * o + 2] != 0.0) { has[static_cast<size_t>(imp_->FirstOwned() + o)] = 1; }
      }
      MPI_Allreduce(MPI_IN_PLACE, has.data(), NG, MPI_INT, MPI_MAX, mpi_.GetComm());
      for (int l = 0; l < NG; l++) { if (has[static_cast<size_t>(l)]) { imp_->EnsureColumn(l); } }
   }

   /// Option I: U from State() and the gap (a new TS, a restart).  Collective.
   void LoadU_()
   {
      PetscScalar *a = nullptr;
      MFEM_VERIFY(VecGetArray(U_imp_, &a) == PETSC_SUCCESS, "VecGetArray failed");
      StateToU_(reinterpret_cast<real_t *>(a));
      EnsureColumns_(reinterpret_cast<real_t *>(a));
      MFEM_VERIFY(VecRestoreArray(U_imp_, &a) == PETSC_SUCCESS, "VecRestoreArray failed");
      u_stale_ = false;
   }

   /// Option I: an accepted or consistent U (3 N, owned) at time t into
   /// State() (slip at the in-plane slot, psi) and the open set (flags: the
   /// branch of the gap row; gap: max(g, 0) there, 0 elsewhere), and the
   /// same cleaned gap back into @a u.  Returns the DOFs whose flag
   /// changed, over all ranks.  Collective.
   long UToState_(real_t t, real_t *u) { return UToState_(t, u, u, nullptr); }

   /// The same from a read-only @a u.  @a u_out: where the cleaned gap is
   /// written (null: nowhere).  @a dirty (when not null): the entries, over
   /// all ranks, whose gap in @a u differs from the cleaned one.
   long UToState_(real_t t, const real_t *u, real_t *u_out, long *dirty)
   {
      const int N = fault_->NumNodes();
      const int spn = FaultOp::StatePerNode, slot = open_->Slot();
      imp_->Evaluate(t, u, imp_V_, imp_tau_, imp_sig_, imp_open_);
      std::vector<char> flags(static_cast<size_t>(N));
      std::vector<real_t> gap(static_cast<size_t>(N));
      const std::vector<char> &old = open_->Flags();
      long changed = 0, ndirty = 0;
      state_.HostReadWrite();
      for (int o = 0; o < N; o++)
      {
         const size_t k = static_cast<size_t>(o);
         const bool is_open = imp_open_[k] != 0 && u[3 * o + 2] > 0.0;
         flags[k] = is_open ? 1 : 0;
         gap[k] = is_open ? u[3 * o + 2] : 0.0;
         if (u[3 * o + 2] != gap[k]) { ndirty++; }
         if (u_out) { u_out[3 * o + 2] = gap[k]; }
         if ((old[k] != 0) != is_open) { changed++; }
         state_(spn * o + FaultOp::SlipIndex + slot) = u[3 * o];
         state_(spn * o + FaultOp::PsiIndex) = u[3 * o + 1];
      }
      MPI_Allreduce(MPI_IN_PLACE, &changed, 1, MPI_LONG, MPI_SUM, mpi_.GetComm());
      if (dirty)
      {
         MPI_Allreduce(MPI_IN_PLACE, &ndirty, 1, MPI_LONG, MPI_SUM, mpi_.GetComm());
         *dirty = ndirty;
      }
      if (changed > 0) { open_->SetFlags(flags); }
      open_->SetGap(gap);
      return changed;
   }

   /// Option I, the consistent start (plan 2.10) at (t, State()): the gap
   /// rows min(sigma, gamma g) = 0 at fixed slip and psi, from the gap the
   /// open set holds; the result into the open set.  Rank 0 logs a change.
   /// Returns the DOFs whose flag changed.  Collective.
   long ConsistentStart_(real_t t, const char *where, long step)
   {
      const int N = fault_->NumNodes();
      std::vector<real_t> u(static_cast<size_t>(std::max(3 * N, 1)));
      StateToU_(u.data());
      EnsureColumns_(u.data());
      real_t g_before = 0.0;
      for (int o = 0; o < N; o++) { g_before = std::max(g_before, u[3 * o + 2]); }
      g_before = mpi_.GlobalMax(g_before);
      const long n_open = imp_->ConsistentStart(t, u.data());
      const long changed = UToState_(t, u.data());
      const real_t g_after = MaxGapGlobal();
      if (mpi_.Rank() == 0 && (changed > 0 || n_open > 0))
      {
         std::cout << "  [cycle.qd] esdirk consistent start (" << where << " " << step
                   << ", t = " << t << " s): " << n_open << " open, " << changed
                   << " flags changed; largest gap " << g_before << " -> " << g_after
                   << " m\n";
      }
      u_stale_ = true;
      return changed;
   }

   /// Option I: after each accepted step (file comment).
   static PetscErrorCode PostStepImplicit_(TS ts)
   {
      PetscFunctionBeginUser;
      void *ctx = nullptr;
      PetscCall(TSGetApplicationContext(ts, &ctx));
      QDPhase *self = static_cast<QDPhase *>(ctx);
      PetscInt step;
      PetscReal t;
      PetscCall(TSGetStepNumber(ts, &step));
      PetscCall(TSGetTime(ts, &t));
      Vec U;
      PetscCall(TSGetSolution(ts, &U));
      PetscScalar *a = nullptr;
      {
         QDTimer tm("QDPostState");
         if (self->cfg_.qd.stage_guess == "extrapolate")
         {
            // PETSc restarts the step (TSRestartStep: no stage history, so no
            // extrapolated first guess) whenever PostStep touched the
            // solution vector.  Read U, and write the cleaned gap back only
            // on the steps where the cleaning changes an entry (then on
            // every rank: PETSc's test is per rank and must agree).
            const PetscScalar *ar = nullptr;
            long dirty = 0;
            PetscCall(VecGetArrayRead(U, &ar));
            self->UToState_(static_cast<real_t>(t), reinterpret_cast<const real_t *>(ar),
                            nullptr, &dirty);
            PetscCall(VecRestoreArrayRead(U, &ar));
            if (dirty > 0)
            {
               PetscCall(VecGetArray(U, &a));
               const std::vector<real_t> &g = self->open_->Gap();
               for (size_t o = 0; o < g.size(); o++) { a[3 * o + 2] = g[o]; }
               PetscCall(VecRestoreArray(U, &a));
            }
         }
         else
         {
            PetscCall(VecGetArray(U, &a));
            self->UToState_(static_cast<real_t>(t), reinterpret_cast<real_t *>(a));
            PetscCall(VecRestoreArray(U, &a));
         }
      }
      self->imp_steps_++;
      self->imp_t_last_ = static_cast<real_t>(t);
      {
         // the kernel's caches at the accepted state, for the monitor
         QDTimer tm("QDPostKernel");
         self->op_->SetTime(static_cast<real_t>(t));
         self->imp_rate_.SetSize(self->state_.Size());
         self->op_->Mult(self->state_, self->imp_rate_);
      }
      if (!self->CallMonitor_(static_cast<int>(step), static_cast<real_t>(t)))
      {
         PetscCall(TSSetConvergedReason(ts, TS_CONVERGED_USER));
      }
      PetscFunctionReturn(PETSC_SUCCESS);
   }
#endif

   /// [cycle.opening]: the rule at an evaluated state x at t (the last
   /// evaluation is at (t, x)): write the solved open slip into x, flip the
   /// DOFs whose true sigma (SigmaN: sigma0 + T_n, no r, no floor) crossed
   /// 0, evaluate again, until nothing flips.  After a change dsig of the
   /// relaxation is the new evaluation's.  Rank 0 logs each flip (DOF,
   /// position, sigma_n', the in-plane shear tau_pre + T it held) and, for a
   /// settle with flips, the passes, the largest static slip |Delta| of the
   /// open-set solves and max V before and after (`where`, `step`: the log's
   /// context).  Returns the DOFs flipped.  Collective.
   long SettleOpenSet_(real_t t, Vector &x, const char *where, long step)
   {
      const int N = fault_->NumNodes();
      const int spn = FaultOp::StatePerNode, slot = open_->Slot();
      const real_t V_before = op_->GetMaxSlipRate();   // global: collective
      real_t delta_max = 0.0;
      long flipped = 0;
      int passes = 0;
      for (int pass = 0;; pass++)
      {
         {
            const Vector &s = op_->GetStageSlip();
            s.HostRead();
            x.HostReadWrite();
            for (int j = 0; j < N; j++)
            {
               if (open_->IsOpen(j)) { x(spn * j + FaultOp::SlipIndex + slot) = s(2 * j + slot); }
            }
         }
         Vector sn;
         SigmaN(sn);
         Vector tau(N);   // in-plane shear of this evaluation (kernel frame)
         {
            const Vector &T = op_->GetTraction();
            const Vector &tp = geom_->GetTauPre();
            T.HostRead();
            tp.HostRead();
            for (int j = 0; j < N; j++) { tau(j) = tp(2 * j + slot) + T(2 * j + slot); }
         }
         const long n = open_->Update(sn);
         if (n == 0) { break; }
         flipped += n;
         passes = pass + 1;
         LogFlips_(t, where, step, pass + 1, n, sn, tau, x);
         MFEM_VERIFY(pass < cfg_.opening.max_passes, "[cycle.opening]: the open set "
                     "did not settle at t = " << t << " s: " << n << " DOFs still flipped "
                     "after " << cfg_.opening.max_passes << " passes ([cycle.opening]."
                     "max_passes); " << open_->NumOpenGlobal() << " open now");
         op_->SetTime(t);
         Vector rate(x.Size());
         op_->Mult(x, rate);
         delta_max = std::max(delta_max, open_->LastStepMax());
      }
      if (flipped > 0)
      {
         const real_t V_after = op_->GetMaxSlipRate();   // collective
         if (mpi_.Rank() == 0)
         {
            std::cout << "  [cycle.opening] settle (" << where << " " << step << ", t = " << t
                      << " s): " << flipped << " flips in " << passes << " passes, "
                      << open_->NumOpenGlobal() << " open; largest static slip |Delta| "
                      << delta_max << " m; max V " << V_before << " -> " << V_after
                      << " m/s\n";
         }
      }
      if (flipped > 0 && relax_) { relax_->SetDelta(RelaxDelta_()); }
      return flipped;
   }

   /// [cycle.opening]: rank 0 prints the DOFs the last OpenSetQD::Update
   /// flipped (at most 8 per pass, then a count), with the state x's psi
   /// there (frozen while open).  Collective.
   void LogFlips_(real_t t, const char *where, long step, int pass, long n,
                  const Vector &sn, const Vector &tau, const Vector &x) const
   {
      constexpr int kFields = 8;   // global DOF, x, y, sigma_n', tau, now open, psi, g
      x.HostRead();
      const std::vector<int> &fl = open_->LastFlipped();
      const std::vector<real_t> &fg = open_->LastFlippedGaps();   // contact; else empty
      std::vector<double> mine(kFields * fl.size());
      for (size_t k = 0; k < fl.size(); k++)
      {
         const int j = fl[k];
         double *r = mine.data() + kFields * k;
         r[0] = static_cast<double>(open_->GlobalIndex(j));
         r[1] = table_.coords(3 * j);
         r[2] = table_.coords(3 * j + 1);
         r[3] = sn(j);
         r[4] = tau(j);
         r[5] = open_->IsOpen(j) ? 1.0 : 0.0;
         r[6] = x(FaultOp::StatePerNode * j + FaultOp::PsiIndex);
         r[7] = (k < fg.size()) ? static_cast<double>(fg[k]) : 0.0;
      }
      const MPI_Comm comm = mpi_.GetComm();
      const int nranks = mpi_.Size();
      int cnt = static_cast<int>(mine.size());
      std::vector<int> counts(nranks), displs(nranks);
      MPI_Gather(&cnt, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);
      int total = 0;
      if (mpi_.Rank() == 0)
      {
         for (int r = 0; r < nranks; r++) { displs[r] = total; total += counts[r]; }
      }
      std::vector<double> all(std::max(total, 1));
      MPI_Gatherv(mine.data(), cnt, MPI_DOUBLE, all.data(), counts.data(), displs.data(),
                  MPI_DOUBLE, 0, comm);
      if (mpi_.Rank() != 0) { return; }
      const int m = total / kFields;
      for (int k = 0; k < std::min(m, 8); k++)
      {
         const double *r = all.data() + kFields * k;
         std::cout << "  [cycle.opening] flip (" << where << " " << step << ", t = " << t
                   << " s, pass " << pass << "): DOF " << static_cast<long long>(r[0])
                   << " at (" << r[1] / 1e3 << ", " << r[2] / 1e3 << ") km "
                   << (r[5] != 0.0 ? "closed -> open" : "open -> closed") << ", sigma_n' "
                   << r[3] / 1e6 << " MPa, shear " << r[4] / 1e6 << " MPa, psi " << r[6];
         // contact (plan 5.5): the released shear is the printed shear of a
         // closed -> open flip; the gap column is the g the solve wanted
         if (open_->Contact()) { std::cout << ", g " << r[7] << " m"; }
         std::cout << "\n";
      }
      if (n > 8) { std::cout << "  [cycle.opening] ... and " << n - 8 << " more flips\n"; }
   }

   /// [cycle.opening]: the open flags into the station fields (empty when off).
   void OpenToStation_(StationFields &f) const
   {
      if (!open_) { f.open.SetSize(0); f.gap.SetSize(0); return; }
      f.open = open_->Mask();
      if (open_->Contact())
      {
         const std::vector<real_t> &gv = open_->Gap();
         f.gap.SetSize(static_cast<int>(gv.size()));
         for (size_t j = 0; j < gv.size(); j++) { f.gap(static_cast<int>(j)) = gv[j]; }
      }
      else { f.gap.SetSize(0); }
   }

   void BuildTable_()
   {
      // Owned XYZ from the geometry (GetFaultDOFCoords3D restricted to owned).
      table_.coords = geom_->fault_dof_coords_3d();
      // Bulk element of each local fault DOF: the Elem1 side of its fault
      // face, in the operator's local order (interior faces, then shared
      // faces, nbf DOFs per face), restricted to owned like any fault field.
      const int nbf = dom_->GetNbfPerFace();
      const int nloc = dom_->GetNumFaultDOFs();
      Vector loc_elem(nloc);
      loc_elem = -1.0;
      const Array<int> &fint = dom_->GetFaultInteriorFaces();
      const Array<int> &fsh = dom_->GetFaultSharedFaces();
      for (int i = 0; i < fint.Size(); i++)
      {
         FaceElementTransformations *T = pmesh_.GetInteriorFaceTransformations(fint[i]);
         MFEM_VERIFY(T != nullptr, "QDPhase: fault face " << fint[i] << " has no transformation");
         for (int k = 0; k < nbf; k++) { loc_elem(i * nbf + k) = T->Elem1No; }
      }
      for (int i = 0; i < fsh.Size(); i++)
      {
         FaceElementTransformations *T = pmesh_.GetSharedFaceTransformations(fsh[i]);
         MFEM_VERIFY(T != nullptr, "QDPhase: shared fault face " << fsh[i]
                     << " has no transformation");
         const int f = fint.Size() + i;
         for (int k = 0; k < nbf; k++) { loc_elem(f * nbf + k) = T->Elem1No; }
      }
      MFEM_VERIFY((fint.Size() + fsh.Size()) * nbf == nloc,
                  "QDPhase: " << fint.Size() + fsh.Size() << " fault faces x " << nbf
                  << " DOFs != " << nloc << " local fault DOFs");
      Vector own_elem;
      dom_->RestrictToOwnedFault(loc_elem, own_elem);
      const int N = own_elem.Size();
      MFEM_VERIFY(table_.coords.Size() == 3 * N,
                  "QDPhase: " << table_.coords.Size() / 3 << " owned coordinates vs "
                  << N << " owned DOFs");
      table_.elem.SetSize(N);
      table_.attr.SetSize(N);
      for (int i = 0; i < N; i++)
      {
         const int e = static_cast<int>(own_elem(i));
         MFEM_VERIFY(e >= 0 && e < pmesh_.GetNE(), "QDPhase: owned fault DOF " << i
                     << " has no local bulk element");
         table_.elem[i] = e;
         table_.attr[i] = pmesh_.GetAttribute(e);
      }
   }

   /// [cycle.pml] (PML plan rev2, P1, Section 7.2): finds the walls on the
   /// mesh, refuses a fault or prescribed-slip face with a vertex inside a
   /// band, and sets the D3 wall speed factor.  Collective.
   void SetupPmlBand_(ParMesh &pmesh)
   {
      const PmlConfig &P = cfg_.pml;
      PmlBandParams prm;
      prm.thickness = P.thickness;
      prm.power = P.power;
      prm.R0 = P.R0;
      // the sponge has no stretch: f_alpha only feeds PmlBand's alpha,
      // which nothing reads then
      prm.f_alpha = P.Sponge() ? 1.0 : P.f_alpha;
      prm.c_ref = P.c_ref;
      auto band = std::make_shared<PmlBand>(pmesh, P.wall_attrs, prm, pmesh.GetComm());
      // every interior boundary element (fault faces, prescribed-slip faces)
      // must lie outside the bands (Section 3.7)
      long long inside = 0;
      Array<int> v;
      for (int be = 0; be < pmesh.GetNBE(); be++)
      {
         if (pmesh.GetFaceInformation(pmesh.GetBdrElementFaceIndex(be)).IsBoundary())
         { continue; }
         pmesh.GetBdrElementVertices(be, v);
         for (int k = 0; k < v.Size(); k++)
         {
            if (band->InBand(pmesh.GetVertex(v[k]))) { inside++; break; }
         }
      }
      MPI_Allreduce(MPI_IN_PLACE, &inside, 1, MPI_LONG_LONG, MPI_SUM, pmesh.GetComm());
      MFEM_VERIFY(inside == 0, "[cycle.pml]: " << inside << " interior faces (fault or "
                  "prescribed slip) have a vertex inside a band; thin the band or "
                  "move the walls");
      // [cycle.finite_width]: the walls move with the substrate (factor 1)
      if (P.Sponge() || P.wall_speed_factor == 1.0 || cfg_.finite_width.On())
      { pml_wall_factor_ = 1.0; }
      else if (P.AutoWallSpeed())
      {
         const auto &m = cfg_.loading.plane_normal;
         pml_wall_factor_ = band->WallSpeedFactor(m.data());
      }
      else { pml_wall_factor_ = P.wall_speed_factor; }
      if (mpi_.Rank() == 0 && P.Sponge())
      {
         std::cout << "[cycle.pml] sponge (dynamic phase only; QD unchanged): "
                   << band->Walls().size() << " walls, L = " << P.thickness << " m, N = "
                   << P.power << ", R0 = " << P.R0 << ", c = " << P.c_ref << " m/s; d0 = "
                   << band->D0() << " 1/s\n";
      }
      else if (mpi_.Rank() == 0)
      {
         std::cout << "[cycle.pml] QD band: " << band->Walls().size() << " walls, L = "
                   << P.thickness << " m, N = " << P.power << ", R0 = " << P.R0
                   << ", f_alpha = " << P.f_alpha << " Hz, c = " << P.c_ref
                   << " m/s; d0 = " << band->D0() << " 1/s, s0 at the wall = "
                   << band->S0Wall() << ", Delta L = " << band->DeltaL()
                   << " m; wall speed factor (D3) = " << pml_wall_factor_ << "\n";
         for (const PmlWall &w : band->Walls())
         {
            std::cout << "  wall attr " << w.attr << ": n = (" << w.n[0] << ", " << w.n[1]
                      << ", " << w.n[2] << "), n . x = " << w.h << " m\n";
         }
      }
      pml_band_ = band;
   }

   /// [cycle.finite_width] (cycle/finite_width.hpp): U_pl1 by the deck's
   /// substrate, then ElasticityDomainOperator::SetFiniteWidth (QD time is
   /// physical time: no offset).  "locked_reference" solves K u = b(t = 1)
   /// with zero slip and beta = 0 first; SetFiniteWidth drops that K and
   /// its factorization.  Collective.
   void SetupFiniteWidth_()
   {
      const FiniteWidthConfig &F = cfg_.finite_width;
      MFEM_VERIFY(pmesh_.Dimension() == 2, "[cycle.finite_width] needs a 2-D map-view "
                  "mesh; this one is " << pmesh_.Dimension() << "-D");
      ParFiniteElementSpace &fes = dom_->GetFESpace();
      if (F.substrate == "loading_plane")
      {
         InterpolateSubstrate(fes, MakePlateSubstrate(cfg_.loading, cfg_.plate_rate), U_pl1_);
      }
      else
      {
         MFEM_VERIFY(F.substrate == "locked_reference", "QDPhase: substrate \""
                     << F.substrate << "\"");
         Vector zero(2 * dom_->GetNumFaultDOFs());
         zero = 0.0;
         ParGridFunction u(&fes);
         dom_->Solve(1.0, zero, u);
         U_pl1_ = u;
      }
      dom_->SetFiniteWidth(F.BetaOverMu(), U_pl1_, nullptr);
   }

   /// Velocity's premise: the Dirichlet loading is g(t) = t g_hat
   /// (the --qd-green build's check, here on the assembled right-hand side).
   void CheckLinearLoading_()
   {
      Vector r0, r1, r2;
      dom_->AssembleDirichletOnlyRHS(r0, 0.0);
      dom_->AssembleDirichletOnlyRHS(r1, 1.0);
      dom_->AssembleDirichletOnlyRHS(r2, 2.0);
      real_t loc[3] = {r0.Normlinf(), 0.0, r1.Normlinf()}, glob[3];
      for (int i = 0; i < r1.Size(); i++) { loc[1] = std::max(loc[1], std::abs(r2(i) - 2.0 * r1(i))); }
      MPI_Allreduce(loc, glob, 3, MPITypeMap<real_t>::mpi_type, MPI_MAX, mpi_.GetComm());
      MFEM_VERIFY(glob[0] == 0.0 && glob[1] <= 1e-12 * glob[2],
                  "QDPhase::Velocity needs Dirichlet loading g(t) = t g_hat: max|b(0)| = "
                  << glob[0] << ", max|b(2) - 2 b(1)| = " << glob[1] << " (max|b(1)| = "
                  << glob[2] << ")");
      loading_checked_ = true;
   }

   static real_t Project_(const real_t *x, const std::array<real_t, 3> &o,
                          const std::array<real_t, 3> &e)
   {
      return (x[0] - o[0]) * e[0] + (x[1] - o[1]) * e[1] + (x[2] - o[2]) * e[2];
   }

   /// The global min and max of v; refuses a spread beyond round-off.
   real_t GlobalUniform_(const Vector &v, const char *name) const
   {
      real_t lo = std::numeric_limits<real_t>::max();
      real_t hi = std::numeric_limits<real_t>::lowest();
      for (int i = 0; i < v.Size(); i++) { lo = std::min(lo, v(i)); hi = std::max(hi, v(i)); }
      lo = mpi_.GlobalMin(lo);
      hi = mpi_.GlobalMax(hi);
      MFEM_VERIFY(hi >= lo, "QDPhase: no fault DOF on any rank");
      MFEM_VERIFY(hi - lo <= 1e-12 * (std::abs(0.5 * (lo + hi)) + 1.0),
                  "QDPhase: per-DOF " << name << " varies over the fault (" << lo << " .. "
                  << hi << "); the QD kernel reads it as one scalar outside SAFS mode");
      return lo;
   }

   void AssertUniform_(const Vector &v, real_t want, const char *name) const
   {
      const real_t got = GlobalUniform_(v, name);
      MFEM_VERIFY(got == want, "QDPhase: per-DOF " << name << " = " << got
                  << " differs from the kernel's scalar " << want);
   }

#ifdef MFEM_USE_PETSC
   /// [cycle.relaxation]: dsig of the last evaluation, the elastic normal
   /// traction the kernel used (the accepted state's after a step: FSAL).
   const Vector &RelaxDelta_() const
   {
      const Vector *nt = op_->GetFrictionNormalTraction();
      MFEM_VERIFY(nt, "QDPhase: [cycle.relaxation] needs the elastic normal traction "
                  "of the last evaluation");
      return *nt;
   }

   /// The monitor at one step; false: stop there.  The skipped step (a
   /// restart's or a re-entry's) counts as `true`.
   bool CallMonitor_(int step, real_t t)
   {
      if (step == skip_monitor_step_)
      {
         skip_monitor_step_ = -1;
         return true;
      }
      QDTimer tm("QDMonitor");
      if (monitor_ && !monitor_(step, t))
      {
         stopped_ = true;
         return false;
      }
      return true;
   }

   /// After each accepted step (TSSolve, before its loop test, so the
   /// reason set here ends the run at this step).
   static PetscErrorCode PostStepCallback_(TS ts)
   {
      PetscFunctionBeginUser;
      void *ctx = nullptr;
      PetscCall(TSGetApplicationContext(ts, &ctx));
      QDPhase *self = static_cast<QDPhase *>(ctx);
      PetscInt step;
      PetscReal t;
      PetscCall(TSGetStepNumber(ts, &step));
      PetscCall(TSGetTime(ts, &t));
      if (self->open_ && self->open_->Contact())
      {
         // D-B5 checks the stage states only: a non-finite rate of the
         // last stage reaches the accepted state through the update (a
         // zero weight times NaN), and every retry of the next step then
         // fails.  Diagnosis (2026-11-07): name such entries (print only).
         Vec U;
         PetscCall(TSGetSolution(ts, &U));
         PetscInt n = 0;
         const PetscScalar *a = nullptr;
         PetscCall(VecGetLocalSize(U, &n));
         PetscCall(VecGetArrayRead(U, &a));
         self->ReportNonFinite_(a, n, "accepted state", static_cast<real_t>(t));
         PetscCall(VecRestoreArrayRead(U, &a));
      }
      if (self->relax_)
      {
         // [cycle.relaxation]: accept the step (the file comment)
         self->relax_->Accept(static_cast<real_t>(t), self->RelaxDelta_());
      }
      if (self->open_)
      {
         // [cycle.opening]: the solved open slip into the TS solution, then
         // the rule; a change is a static step, so the next step restarts
         Vec U;
         PetscCall(TSGetSolution(ts, &U));
         PetscInt n = 0;
         PetscCall(VecGetLocalSize(U, &n));
         PetscScalar *a = nullptr;
         PetscCall(VecGetArray(U, &a));
         long flips = 0;
         {
            Vector x(reinterpret_cast<real_t *>(a), static_cast<int>(n));
            flips = self->SettleOpenSet_(static_cast<real_t>(t), x, "after step",
                                         static_cast<long>(step));
         }
         PetscCall(VecRestoreArray(U, &a));
         if (flips > 0)
         {
            // The flip changed the state by a static step that no error
            // control saw (V can jump by orders of magnitude), so the next
            // step starts again from [time].dt_initial, as after a switch back.
            PetscReal dt_old = 0.0;
            PetscCall(TSGetTimeStep(ts, &dt_old));
            PetscCall(TSRestartStep(ts));
            PetscCall(TSSetTimeStep(ts, self->cfg_.dt_init));
            if (self->mpi_.Rank() == 0)
            {
               std::cout << "  [cycle.opening] next dt " << dt_old << " -> " << self->cfg_.dt_init
                         << " s ([time].dt_initial) after the flips at step " << step << "\n";
            }
         }
      }
      if (!self->CallMonitor_(static_cast<int>(step), static_cast<real_t>(t)))
      {
         PetscCall(TSSetConvergedReason(ts, TS_CONVERGED_USER));
      }
      PetscFunctionReturn(PETSC_SUCCESS);
   }

   std::unique_ptr<PetscODESolver> ode_;
   /// [cycle.qd] integrator = "esdirk" (option I): the implicit fault
   /// system, its TS and its state vector (3 per owned DOF); null with "dp5".
   std::unique_ptr<QDImplicitFault> imp_;
   TS ts_imp_ = nullptr;
   Vec U_imp_ = nullptr;
   bool u_stale_ = true;   ///< U must be rebuilt from State() and the gap
   std::vector<real_t> imp_V_, imp_tau_, imp_sig_;   ///< UToState_ scratch
   std::vector<char> imp_open_;
   Vector imp_rate_;                                 ///< PostStepImplicit_ scratch
   /// ImplicitPhase(): the marks of the last call.
   QDImplicitFault::Counts imp_mark_;
   long imp_steps_ = 0, imp_mark_steps_ = 0;
   real_t imp_mark_t0_ = 0.0, imp_t_last_ = 0.0;
   double imp_wall_ = 0.0;
   bool imp_phase_new_ = true;
   /// [cycle.qd] exact_retry (declared after ode_, so destroyed before it).
   ExactRetry retry_;
   long retry_nonfinite_ = 0;
   /// ReportNonFinite_: lines printed on this rank (capped).
   mutable long nonfinite_lines_ = 0;
#endif

   ParMesh &pmesh_;
   MPIContext &mpi_;
   const CycleConfig &cfg_;
   BP5Params seed_;
   std::unique_ptr<DomainOp> dom_;
   std::shared_ptr<const PmlBand> pml_band_;   ///< [cycle.pml] (null: off)
   real_t pml_wall_factor_ = 1.0;              ///< D3 wall speed factor
   Vector U_pl1_;   ///< [cycle.finite_width] substrate shape (SubstrateShape)
   std::unique_ptr<FaultGeometry<ParMesh>> geom_;
   FaultDOFTable table_;
   std::unique_ptr<DieterichRuinaFriction> friction_;
   std::unique_ptr<AgingLawPsi> aging_;
   std::unique_ptr<FaultOp> fault_;
   std::unique_ptr<PBP5SEASOp> op_;
   Vector state_;
   Vector sigma0_;   ///< sigma0 per owned fault DOF (Setup), for SigmaN
   std::unique_ptr<NormalStressRelaxation> relax_;   ///< [cycle.relaxation] or null
   std::unique_ptr<OpenSetQD> open_;                ///< [cycle.opening] or null
   mutable Vector gap_refresh_;                     ///< Displacement(), contact
   long reenter_flips_ = 0;
   real_t eta_ = 0.0;
   bool matrix_ = false;
   bool loading_checked_ = false;
   std::unique_ptr<ParGridFunction> vel_gf_;
   Monitor monitor_;
   bool stopped_ = false;
   int skip_monitor_step_ = -1;
};

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_QD_PHASE_HPP

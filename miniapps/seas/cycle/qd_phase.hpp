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
//     op.SetInitialCondition(state)
//   InitTS(dt0): PETSc TS (rk 5dp from the options file), adaptive, dt0
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
// Restart: a checkpoint written at step k, and a run restarted from it, both
// re-evaluate the first RK stage of step k + 1 (TSRestartStep) instead of
// reusing the FSAL stage, so the restarted run repeats the original run's
// steps bit for bit (test_cycle_restart).  A new TS evaluates its first
// stage anyway, so without the restart at checkpoint steps a restarted run
// would match the uninterrupted one only to round-off (seas_bp5_full: about
// 1e-11 relative after 25 steps).  The price: every checkpoint step moves the
// run by the same round-off against a run without checkpoints.
//
// This file and dyn_phase.hpp are the only cycle files that include solver
// headers (rule M4).

#ifndef MFEM_SEAS_CYCLE_QD_PHASE_HPP
#define MFEM_SEAS_CYCLE_QD_PHASE_HPP

#include "mfem.hpp"
#include "cycle_config.hpp"
#include "fault_fields.hpp"
#include "prestress.hpp"
#include "fd_cycle_output.hpp"

#include "../solver/seas_operator.hpp"
#include "../domain/elasticity_operator.hpp"
#include "../fault/fault_geometry.hpp"
#include "../fault/rate_state_fault.hpp"
#include "../friction/dieterich_ruina.hpp"
#include "../friction/state_evolution.hpp"
#include "../config/bp5_params.hpp"
#include "../io/bp5_benchmark_output.hpp"   // Probe2DInterpolator
#include "../common/mpi_context.hpp"

#ifdef MFEM_USE_PETSC
#include "petsc.h"
#endif

#include <algorithm>
#include <cmath>
#include <functional>
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
      // fault 3 (BCMode::FarField), exactly as seas_bp5_full builds it.
      dom_ = std::make_unique<DomainOp>(pmesh, s.mesh.order, lambda, mu, cfg.plate_rate,
                                        seed_.Wf, seed_.lf, DGMethod::IP, st,
                                        BCMode::FarField, BasisType::GaussLobatto);
      if (st == SolverType::MUMPS_BLR) { dom_->SetBLRTol(cfg.qd.blr_tol); }
      if (cfg.penalty_factor != 1.0) { dom_->SetPenaltyFactor(cfg.penalty_factor); }

      geom_ = std::make_unique<FaultGeometry<ParMesh>>(*dom_, seed_, &mpi_,
                                                       /*compute_bp5_params=*/false);
      BuildTable_();
   }

   /// The owned fault DOFs in the operator's owned order.
   const FaultDOFTable &Table() const { return table_; }
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
      // The kernel's scalars (see the file comment).
      AssertUniform_(par.rs.sigma_n_eff, seed_.sigma_n, "sigma_n");
      AssertUniform_(par.rs.b, seed_.b, "b");
      AssertUniform_(par.rs.f_0, seed_.f0, "f_0");
      AssertUniform_(par.rs.V_0, seed_.V0, "V_0");
      AssertUniform_(pre.sigma0, seed_.sigma_n, "sigma0");
      eta_ = GlobalUniform_(par.rs.eta, "eta");

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

      op_ = std::make_unique<PBP5SEASOp>(dom_.get(), fault_.get(), &mpi_);
      op_->SetElasticSigmaN(true);

      // R-003 of seas_bp5_full: allocate at least one entry so PETSc's
      // PlaceMemory/ResetMemory never sees an empty Memory on a rank without
      // fault DOFs; SetSize keeps the allocation when shrinking.
      const int actual = fault_->StateSize();
      state_.SetSize(std::max(actual, 1));
      state_.SetSize(actual);
      op_->SetInitialCondition(state_);
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
   PBP5SEASOp &Operator() { return *op_; }
   real_t MaxV() const { return op_->GetMaxSlipRate(); }
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
      return op_->GetMaxSlipRate();
   }
   /// What the last evaluation left, for the bridge (owned fault DOFs; slip
   /// and slip rate as kernel (dip, strike) pairs, 2 N values).
   const Vector &Displacement() const { return op_->GetDisplacement(); }
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
      for (int i = 0; i < N; i++) { sn(i) = nt ? seed_.sigma_n + (*nt)(i) : seed_.sigma_n; }
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
   /// SigmaN (for the snapshots; the station rows do not use it).
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
   }

   /// Stations on this rank's owned fault DOFs, in the plane frame of the
   /// deck: 2-D coordinates (x - o).e1 and (x - o).e2 for stations and DOFs,
   /// then seas_bp5_full's Probe2DInterpolator (face-basis match, else the
   /// nearest DOF).
   std::vector<StationLocation> LocateStations(const OutputConfig &oc) const
   {
      const int N = table_.Size();
      Vector x2(N), x3(N);
      for (int i = 0; i < N; i++)
      {
         const real_t *x = table_.coords.GetData() + 3 * i;
         x2(i) = Project_(x, oc.plane_origin, oc.plane_axis1);
         x3(i) = Project_(x, oc.plane_origin, oc.plane_axis2);
      }
      std::vector<Probe2DInterpolator::Station> st;
      for (const StationSpec &s : oc.stations)
      {
         st.push_back({s.name, Project_(s.x.data(), oc.plane_origin, oc.plane_axis1),
                       Project_(s.x.data(), oc.plane_origin, oc.plane_axis2)});
      }
      Probe2DInterpolator interp(x2, x3, st, dom_->GetNbfPerFace(),
                                 BasisType::GaussLobatto, false);
      std::vector<StationLocation> out(st.size());
      for (size_t s = 0; s < st.size(); s++)
      {
         StationLocation &L = out[s];
         const int si = static_cast<int>(s);
         L.exact = interp.HasExactMatch(si);
         L.distance = interp.GetMatchDistance(si);
         if (N == 0 || (interp.GetNearestDOF(si) < 0 && !L.exact))
         {
            L.distance = std::numeric_limits<real_t>::max();
         }
         if (L.exact)
         {
            L.start = interp.GetFaceStart(si);
            L.weights = interp.GetExactWeights(si);
         }
         L.nearest = (N > 0) ? interp.GetNearestDOF(si) : -1;
      }
      return out;
   }

#ifdef MFEM_USE_PETSC
   /// PETSc TS as seas_bp5_full --petsc-ts sets it up.  The options file
   /// (cfg.qd.petsc_options) must already be loaded by MFEMInitializePetsc.
   void InitTS(real_t dt0)
   {
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
      MFEM_VERIFY(TSSetTimeStep(ts, dt0) == PETSC_SUCCESS, "TSSetTimeStep failed");
      MFEM_VERIFY(TSSetApplicationContext(ts, this) == PETSC_SUCCESS,
                  "TSSetApplicationContext failed");
      MFEM_VERIFY(TSSetPostStep(ts, &QDPhase::PostStepCallback_) == PETSC_SUCCESS,
                  "TSSetPostStep failed");
   }

   /// Plan 6.3, the switch back: SetFaultState, Solve(t), and a new TS at t
   /// (first step dt0, step counter `step`).  Returns the global max V of
   /// the solve.  Collective.
   real_t Reenter(const Vector &slip_k, const Vector &psi, real_t t, int step, real_t dt0)
   {
      SetFaultState(slip_k, psi);
      const real_t Vmax = Solve(t);
      InitTS(dt0);
      TS ts = *ode_;
      MFEM_VERIFY(TSSetStepNumber(ts, step) == PETSC_SUCCESS, "TSSetStepNumber failed");
      skip_monitor_step_ = step;   // that time's rows are the dynamic regime's
      return Vmax;
   }

   /// Continue from a checkpoint: the TS step counter and a fresh first stage.
   void RestoreTS(int step)
   {
      TS ts = *ode_;
      MFEM_VERIFY(TSSetStepNumber(ts, step) == PETSC_SUCCESS, "TSSetStepNumber failed");
      RestartStage();
      skip_monitor_step_ = step;   // its output is already in the files
   }

   /// The next step re-evaluates its first RK stage (see the file comment).
   void RestartStage()
   {
      TS ts = *ode_;
      MFEM_VERIFY(TSRestartStep(ts) == PETSC_SUCCESS, "TSRestartStep failed");
   }

   /// Steps until t_final, max_steps, or `mon` returns false.  t and dt
   /// return the final time and the next step size.
   void Run(real_t &t, real_t &dt, real_t t_final, int max_steps, const Monitor &mon)
   {
      monitor_ = mon;
      stopped_ = false;
      // The start step (TSSolve's first TSMonitor call in seas_bp5_full).
      TS ts = *ode_;
      PetscInt s0 = 0;
      MFEM_VERIFY(TSGetStepNumber(ts, &s0) == PETSC_SUCCESS, "TSGetStepNumber failed");
      MFEM_VERIFY(TSSetTime(ts, t) == PETSC_SUCCESS && TSSetTimeStep(ts, dt) == PETSC_SUCCESS,
                  "TSSetTime/TSSetTimeStep failed");
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
      TSGetTimeStep(*ode_, &dt);
      return dt;
   }
   int StepNumber() const
   {
      PetscInt n;
      TSGetStepNumber(*ode_, &n);
      return static_cast<int>(n);
   }
   int Rejections() const
   {
      PetscInt n;
      TSGetStepRejections(*ode_, &n);
      return static_cast<int>(n);
   }
   TS GetTS() const { return *ode_; }
#endif

private:
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
   /// The monitor at one step; false: stop there.  The skipped step (a
   /// restart's or a re-entry's) counts as `true`.
   bool CallMonitor_(int step, real_t t)
   {
      if (step == skip_monitor_step_)
      {
         skip_monitor_step_ = -1;
         return true;
      }
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
      if (!self->CallMonitor_(static_cast<int>(step), static_cast<real_t>(t)))
      {
         PetscCall(TSSetConvergedReason(ts, TS_CONVERGED_USER));
      }
      PetscFunctionReturn(PETSC_SUCCESS);
   }

   std::unique_ptr<PetscODESolver> ode_;
#endif

   ParMesh &pmesh_;
   MPIContext &mpi_;
   const CycleConfig &cfg_;
   BP5Params seed_;
   std::unique_ptr<DomainOp> dom_;
   std::unique_ptr<FaultGeometry<ParMesh>> geom_;
   FaultDOFTable table_;
   std::unique_ptr<DieterichRuinaFriction> friction_;
   std::unique_ptr<AgingLawPsi> aging_;
   std::unique_ptr<FaultOp> fault_;
   std::unique_ptr<PBP5SEASOp> op_;
   Vector state_;
   real_t eta_ = 0.0;
   Monitor monitor_;
   bool stopped_ = false;
   int skip_monitor_step_ = -1;
};

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_QD_PHASE_HPP

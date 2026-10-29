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
#include "cycle_frames.hpp"
#include "loading.hpp"
#include "fault_fields.hpp"
#include "prestress.hpp"
#include "fd_cycle_output.hpp"
#include "bench_output.hpp"

#include "../solver/seas_operator.hpp"
#include "../domain/elasticity_operator.hpp"
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
      // fault 3 (BCMode::FarField), exactly as seas_bp5_full builds it, with
      // the cycle's one fault orientation system (phase A): the operator
      // builds the FaultFrameField; the dynamic regime reads it.  Loading:
      // [cycle.loading] (cycle/loading.hpp; default the BP5 function).
      // The prescribed-slip faces must lie in the loading's band (loading.hpp).
      CheckCreepInBand(pmesh, {5}, cfg.loading);
      FaultFrameSource frames;
      frames.config = std::make_shared<const FaultFrameConfig>(CycleFrameConfig(cfg));
      dom_ = std::make_unique<DomainOp>(pmesh, s.mesh.order, lambda, mu, cfg.plate_rate,
                                        seed_.Wf, seed_.lf, DGMethod::IP, st,
                                        BCMode::FarField, BasisType::GaussLobatto, frames,
                                        MakeCycleLoading(cfg));
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
      // QD speed items 1, 3, 4 (seas_bp5_full --qd-operators, --qd-green*):
      // after SetElasticSigmaN, which they build N_u, N_delta for.
      matrix_ = (cfg_.qd.operators == "matrix");
      if (matrix_)
      {
         op_->SetQDOperators(true);
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
      if (cfg_.prestress.AMax() && !cfg_.IsGold())
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
   /// and slip rate as kernel (dip, strike) pairs, 2 N values).  The
   /// displacement is solved here when the Green's function left it stale
   /// (collective then; every rank must ask).
   const Vector &Displacement() const
   {
      op_->RefreshDisplacement();
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
         dom_->SolveStage(1.0, V, *vel_gf_, T, nullptr);
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

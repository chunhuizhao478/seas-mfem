// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dyn_phase.hpp — fully dynamic BP5 cycle (plan rev 3, Sections 3.3 e and
// 6.1-6.2): DynPhase, adapter 2, the dynamic regime.
//
// It builds the Nitsche operator W on the QD regime's mesh, boundaries,
// material and penalty, as drivers/nitsche_dyn_driver.cpp builds it
// (:700-930), with the cycle's fault data instead of a benchmark's:
//
//   DynPhase(pmesh, cfg, table, par, pre):
//     bc  = {Dirichlet {5}, natural {1}, fault 3}: the QD legacy
//           constructor's lists (BCMode::FarField, stage 1 of plan 6.4)
//     bc.default_dirichlet_func = (x, tau) -> MakeBP5DirichletFunc(Vp)(x, t0 + tau)
//     fs  = closure RateStateClosure (aging law, the deck's b, f0, V0),
//           orientation and frame from [fault_geometry],
//           prestress zero (H2 below), parameters (a, Dc) = par at the DOF
//     W   = NitscheWaveOperator(pmesh, p, material, bc, {penalty}, &fs)
//     the owned fault DOFs of W must be the QD table's, point for point
//     frames = (eps_f, n_k, t1_k, t2_k) per owned DOF        (state_bridge.hpp)
//     W.SetBackgroundTraction(frames.TractionToGlobal(pre.tau0, pre.sigma0))  # H2
//     dt  = ComputeMaxDt(cfl C(p)); W.SetTimeStep(dt)  (once: dt is fixed)
//     [numerics] kv_eta > 0: W.SetKelvinVoigt(...)    (v_D = 0 on Dirichlet faces)
//     [friction.prakash_clifton]: W.SetPrakashClifton(v_star, length)
//   Enter(u_s, d, V, psi, sigma_n, t_s):                     (plan 6.2)
//     t0 = t_s; tau = 0; W.SetTime(0)        # the loads read t_s + tau (R10)
//     W.SetFaultState(d, V, psi); Prakash-Clifton: sigma* = sigma_n
//     U = u_s; U' = 0; a fresh CentralDifferenceSolver, Init(W)
//   Step(): one central-difference step; W.ImplicitSolve advances the fault
//
// The wall and creep loads are evaluated at t0 + tau, never at tau alone:
// W's clock stays the stepper clock tau, so a large t0 does not round the
// step (plan J3).  The Kelvin-Voigt term sees v_D = 0 on the Dirichlet
// faces (no bc_velocity, the Nitsche driver's production setting); the
// true wall and creep velocities are Vp/2 and Vp, about 1e-9 m/s.
//
// Station fields in the dynamic regime use the kernel components of the
// end-of-step traction X* + gamma V (the friction strength along -V, as the
// Nitsche driver's stations; gamma includes sigma_F^eta with Kelvin-Voigt)
// and eta = 0, so the station row -(tau + eta V) has the QD rows' sign.
//
// This file and qd_phase.hpp are the only cycle files that include solver
// headers (rule M4).

#ifndef MFEM_SEAS_CYCLE_DYN_PHASE_HPP
#define MFEM_SEAS_CYCLE_DYN_PHASE_HPP

#include "mfem.hpp"
#include "cycle_config.hpp"
#include "fault_fields.hpp"
#include "prestress.hpp"
#include "state_bridge.hpp"
#include "fd_cycle_output.hpp"

#include "../domain/nitsche_wave_operator.hpp"
#include "../domain/boundary_config.hpp"
#include "../domain/domain_config.hpp"
#include "../solver/nitsche_cfl.hpp"
#include "../friction/rate_state_closure.hpp"
#include "../friction/dieterich_ruina.hpp"
#include "../fault/fault_frame.hpp"
#include "../fault/fault_orientation.hpp"
#include "../dynamic/heterogeneous_material.hpp"

#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <string>

namespace mfem
{
namespace seas
{
namespace cycle
{

class DynPhase
{
public:
   /// The measured stable constant C(p) of the explicit step, p = 1..4
   /// (Nitsche P3; the same values as drivers/nitsche_dyn_driver.cpp:91).
   static constexpr real_t kCeig[5] = {0.0, 0.113, 0.212, 0.282, 0.332};

   /// `table`, `par` and `pre` are the QD regime's, on the same owned fault
   /// DOFs.  `dt_scale` multiplies the CFL step (tests only).
   DynPhase(ParMesh &pmesh, const CycleConfig &cfg, const FaultDOFTable &table,
            const FaultParams &par, const Prestress &pre, real_t dt_scale = 1.0)
      : pmesh_(pmesh), cfg_(cfg), t0_(std::make_shared<real_t>(0.0))
   {
      const spatial::SpatialFrictionConfig &s = cfg.shared;
      const spatial::RateStateBlock &rsb = *s.rate_state;
      const DynamicConfig &dc = cfg.dyn;
      const int p = s.mesh.order;
      MFEM_VERIFY(p >= 1 && p <= 4, "DynPhase: order " << p << " is outside 1..4");
      MFEM_VERIFY(rsb.state_evolution == spatial::StateEvolutionKind::AgingLaw,
                  "DynPhase: the QD regime runs the aging law; the deck asks for another");
      material_ = std::make_unique<MaterialField>(MaterialField::MakeConstant(
         s.material_fallback.lambda, s.material_fallback.mu, s.material_fallback.rho));

      // Boundaries: the QD legacy constructor's lists, loads read at t0 + tau.
      bc_.dirichlet_attrs = {5};
      bc_.natural_attrs = {1};
      bc_.fault_attr = cfg.fault_attr;
      {
         const DirichletFunc load = MakeBP5DirichletFunc(cfg.plate_rate);
         std::shared_ptr<const real_t> t0 = t0_;
         bc_.default_dirichlet_func = [load, t0](const Vector &x, real_t tau, Vector &u)
         { load(x, *t0 + tau, u); };
      }

      // Fault setup: closure, orientation, frame, parameters.
      fr_ = std::make_unique<DieterichRuinaFriction>(DieterichRuinaFriction::Constants{
         rsb.V_0_default, rsb.f_0_default, rsb.b_default, rsb.Dc_default});
      closure_ = std::make_unique<RateStateClosure>(*fr_, RateStateClosure::StateLaw::Aging);
      NitscheFaultSetup fs;
      fs.closure = closure_.get();
      FaultOrientationConfig &oc = fs.orientation;
      for (int i = 0; i < 3; i++)
      {
         oc.ref_normal[i] = dc.ref_normal[i];
         oc.up[i] = dc.up[i];
         fs.up[i] = dc.up[i];
      }
      oc.pin = FaultOrientationConfig::Pin::Reference;
      if (dc.pin_min_margin > 0.0) { oc.c_min = dc.pin_min_margin; }
      oc.fault_attrs = {cfg.fault_attr};
      fs.sigma_min = dc.sigma_min;
      // The prestress arrives through H2 after setup; the projection sees 0.
      fs.prestress = [](const real_t[3], real_t sig[9])
      { for (int i = 0; i < 9; i++) { sig[i] = 0.0; } };
      const int N = table.Size();
      MFEM_VERIFY(par.rs.a.Size() == N && par.rs.Dc.Size() == N && pre.tau0.Size() == 2 * N
                  && pre.sigma0.Size() == N,
                  "DynPhase: parameters are not on this rank's owned fault DOFs");
      std::map<std::array<real_t, 3>, int> where;
      for (int i = 0; i < N; i++)
      {
         const real_t *x = table.coords.GetData() + 3 * i;
         where[{x[0], x[1], x[2]}] = i;
      }
      const Vector a = par.rs.a, Dc = par.rs.Dc;
      fs.parameters = [where, a, Dc](const real_t x[3], FaultPointState &st)
      {
         const auto it = where.find({x[0], x[1], x[2]});
         MFEM_VERIFY(it != where.end(), "DynPhase: the dynamic operator's fault DOF at ("
                     << x[0] << ", " << x[1] << ", " << x[2] << ") is not an owned "
                     "fault DOF of the QD regime on this rank");
         st.a = a(it->second);
         st.Dc = Dc(it->second);
      };

      DomainConfig dcfg;
      dcfg.penalty_factor = cfg.penalty_factor;
      W_ = std::make_unique<NitscheWaveOperator>(pmesh, p, *material_, bc_, dcfg, &fs);
      MFEM_VERIFY(W_->Orientation().GetReport().junction_edges == 0,
                  "DynPhase: the fault has " << W_->Orientation().GetReport().junction_edges
                  << " junction edges ([fault_geometry] junction_policy = \"refuse\")");
      MFEM_VERIFY(W_->NumOwnedFaultDOFs() == N, "DynPhase: " << W_->NumOwnedFaultDOFs()
                  << " owned fault DOFs in the dynamic operator, " << N << " in the QD regime");
      for (int i = 0; i < 3 * N; i++)
      {
         MFEM_VERIFY(W_->FaultPositions()(i) == table.coords(i),
                     "DynPhase: owned fault DOF " << i / 3 << " sits at a different point "
                     "in the two operators (the bridge needs the same owned order)");
      }

      BuildFrames_();
      Vector S0;
      frames_.TractionToGlobal(pre.tau0, pre.sigma0, S0);
      W_->SetBackgroundTraction(S0);   // hook H2

      dt_cfl_ = ComputeMaxDt(W_->Sample(), pmesh, p, dc.cfl * kCeig[p]);
      dt_ = dt_scale * dt_cfl_;
      W_->SetTimeStep(dt_);

      if (dc.kv_eta > 0.0)
      {
         W_->SetKelvinVoigt(dc.kv_eta,
                            dc.kv_scale == "node" ? NitscheWaveOperator::KVScale::Node
                                                  : NitscheWaveOperator::KVScale::Dt,
                            dc.kv_layer_width,
                            dc.kv_length == "edge" ? NitscheWaveOperator::KVLength::MeanEdge
                                                   : NitscheWaveOperator::KVLength::Inscribed);
      }
      if (dc.pc) { W_->SetPrakashClifton(dc.pc_v_star, dc.pc_length); }

      const int n = W_->Height();
      U_.SetSize(n);
      Ud_.SetSize(n);
      U_ = 0.0;
      Ud_ = 0.0;
   }

   /// Plan 6.2: hand the static QD state at t_s over.  u_s: the QD
   /// displacement; d, V: global three-vectors (3 N); psi, sigma_n: N.
   void Enter(const Vector &u_s, const Vector &d, const Vector &V, const Vector &psi,
              const Vector &sigma_n, real_t t_s)
   {
      MFEM_VERIFY(u_s.Size() == W_->Height(), "DynPhase::Enter: displacement size "
                  << u_s.Size() << " != " << W_->Height());
      *t0_ = t_s;
      tau_ = 0.0;
      W_->SetTime(0.0);
      MFEM_VERIFY(W_->GetTime() == 0.0, "DynPhase::Enter: the operator's clock is "
                  << W_->GetTime() << ", not 0 (the loads would read t_s + clock)");
      W_->SetFaultState(d, V, psi);
      if (W_->HasPrakashClifton()) { W_->SetPrakashCliftonState(sigma_n); }
      U_ = u_s;
      Ud_ = 0.0;
      ode_ = std::make_unique<CentralDifferenceSolver>();
      ode_->Init(*W_);
      steps_ = 0;
   }

   /// One central-difference step: tau += dt.
   void Step()
   {
      MFEM_VERIFY(ode_, "DynPhase::Step before Enter");
      real_t h = dt_;
      ode_->Step(U_, Ud_, tau_, h);
      steps_++;
   }

   real_t Tau() const { return tau_; }
   real_t T0() const { return *t0_; }
   real_t Time() const { return *t0_ + tau_; }
   real_t Dt() const { return dt_; }
   real_t DtCFL() const { return dt_cfl_; }
   long Steps() const { return steps_; }
   /// Global max |V| of the last step (0 before the first).
   real_t MaxV() const { return W_->LastStep().max_V; }

   /// The last step's global fault diagnostics, in the cycle's own struct
   /// so the driver does not read the operator's record (rule M4).
   struct StepInfo
   {
      real_t max_V = 0.0;
      real_t max_V_dt_over_Dc = 0.0;   ///< the operator aborts above 2
      real_t min_sigma_eff = 0.0, max_sigma_eff = 0.0;
      real_t min_sigma_pc = 0.0, max_sigma_pc = 0.0;   ///< 0 without Prakash-Clifton
   };
   StepInfo LastStepInfo() const
   {
      const NitscheStepRecord &r = W_->LastStep();
      StepInfo s;
      s.max_V = r.max_V;
      s.max_V_dt_over_Dc = r.max_V_dt_over_Dc;
      s.min_sigma_eff = r.min_sigma_eff;
      s.max_sigma_eff = r.max_sigma_eff;
      s.min_sigma_pc = r.min_sigma_pc;
      s.max_sigma_pc = r.max_sigma_pc;
      return s;
   }

   /// Owned fault DOFs: positions (3 N), |V| (N), slip as kernel pairs (2 N),
   /// psi (N).
   int NumOwnedFaultDOFs() const { return W_->NumOwnedFaultDOFs(); }
   const Vector &FaultPositions() const { return W_->FaultPositions(); }
   void SlipRateMagnitude(Vector &v) const
   {
      const Vector &V = W_->SlipRate();
      const int N = W_->NumOwnedFaultDOFs();
      v.SetSize(N);
      for (int j = 0; j < N; j++)
      {
         v(j) = std::sqrt(V(3 * j) * V(3 * j) + V(3 * j + 1) * V(3 * j + 1)
                          + V(3 * j + 2) * V(3 * j + 2));
      }
   }
   /// Slip as global three-vectors (3 N), Bridge::ToQD's input.
   const Vector &Slip() const { return W_->Slip(); }
   void SlipKernel(Vector &sk) const { frames_.ToKernel(W_->Slip(), sk); }
   const Vector &Psi() const { return W_->State(); }

   bool HasKelvinVoigt() const { return W_->HasKelvinVoigt(); }
   real_t KelvinVoigtBetaMin() const { return W_->KelvinVoigtBetaMin(); }
   real_t KelvinVoigtBetaMax() const { return W_->KelvinVoigtBetaMax(); }
   long long KelvinVoigtElements() const { return W_->KelvinVoigtLayerElements(); }
   bool HasPrakashClifton() const { return W_->HasPrakashClifton(); }
   int NumSharedFaultFaces() const { return W_->NumSharedFaultFacesGlobal(); }

   NitscheWaveOperator &Op() { return *W_; }
   const NitscheWaveOperator &Op() const { return *W_; }
   const FrameMap &Frames() const { return frames_; }
   const Vector &U() const { return U_; }
   const Vector &Ud() const { return Ud_; }
   const DieterichRuinaFriction &Friction() const { return *fr_; }
   const BoundaryConfig &Boundaries() const { return bc_; }

   /// Station fields on owned fault DOFs (kernel components; see the file
   /// comment).  Before the first step the traction is T_h(U, d) plus the
   /// background at tau.
   void FillStationFields(StationFields &f)
   {
      const int N = W_->NumOwnedFaultDOFs();
      Vector sk, vk, T, tau_k, sn;
      frames_.ToKernel(W_->Slip(), sk);
      frames_.ToKernel(W_->SlipRate(), vk);
      const Vector &Xs = W_->LastShear();
      // steps_ > 0: after a re-Enter (P4) the operator still holds the last
      // phase's shear until this phase's first step.
      if (steps_ > 0 && Xs.Size() == 2 * N)
      {
         T.SetSize(3 * N);
         const Vector &V = W_->SlipRate();
         for (int j = 0; j < N; j++)
         {
            const real_t *dv = &W_->FaultFrames()(9 * j + 3), *sv = &W_->FaultFrames()(9 * j + 6);
            const real_t g = W_->ClosureGamma(j);
            const real_t X1 = Xs(2 * j) + g * fault_frame::Dot(&V(3 * j), dv);
            const real_t X2 = Xs(2 * j + 1) + g * fault_frame::Dot(&V(3 * j), sv);
            fault_frame::Embed(X1, X2, dv, sv, &T(3 * j));
         }
      }
      else
      {
         Vector S;
         W_->TrialTraction(U_, W_->Slip(), T);
         W_->BackgroundTraction(tau_, S);
         T += S;
      }
      frames_.TractionToKernel(T, tau_k, sn);
      for (Vector *v : {&f.slip_dip, &f.slip_strike, &f.V_dip, &f.V_strike,
                        &f.tau_a_dip, &f.tau_a_strike, &f.tau_b_dip, &f.tau_b_strike,
                        &f.theta})
      {
         v->SetSize(N);
      }
      for (int j = 0; j < N; j++)
      {
         f.slip_dip(j) = sk(2 * j);
         f.slip_strike(j) = sk(2 * j + 1);
         f.V_dip(j) = vk(2 * j);
         f.V_strike(j) = vk(2 * j + 1);
         f.tau_a_dip(j) = tau_k(2 * j);
         f.tau_a_strike(j) = tau_k(2 * j + 1);
         f.tau_b_dip(j) = 0.0;
         f.tau_b_strike(j) = 0.0;
         f.theta(j) = fr_->PsiToTheta(W_->State()(j), W_->FaultDc()(j));
      }
      f.eta = 0.0;
   }

   /// Plan R8: max over local fault DOFs of |eps_f t_k (this copy) - eps_f
   /// t_k (owner's copy)|, global.  0 up to round-off means the QD kernel's
   /// copy of kernel components and the operator's copy of global vectors
   /// put the same slip on every face.
   real_t GhostFrameMismatch() const { return ghost_mismatch_; }
   /// max over local fault DOFs of |eps_f n_k - n| with n the oriented normal.
   real_t NormalMismatch() const { return normal_mismatch_; }

private:
   void BuildFrames_()
   {
      ElasticityDomainOperator<ParMesh> &E = W_->Elasticity();
      DenseMatrix Bk;
      E.GetFaultDOFBasis(Bk);
      const int nd = E.GetNumFaultDOFs();
      const int nbf = E.GetNbfPerFace();
      MFEM_VERIFY(Bk.Height() == 9 && Bk.Width() == nd, "DynPhase: kernel basis is "
                  << Bk.Height() << " x " << Bk.Width() << ", expected 9 x " << nd);
      Vector loc(FrameMap::kPerDOF * nd), gl(6 * nd);
      real_t nmis = 0.0;
      for (int k = 0; k < nd; k++)
      {
         const real_t nk[3] = {Bk(0, k), Bk(1, k), Bk(2, k)};
         const real_t *n = W_->Orientation().Normal(k / nbf);
         const real_t eps = fault_frame::EpsilonF(n, nk);
         real_t *f = &loc(FrameMap::kPerDOF * k);
         f[0] = eps;
         for (int c = 0; c < 3; c++)
         {
            f[1 + c] = Bk(c, k);
            f[4 + c] = Bk(3 + c, k);
            f[7 + c] = Bk(6 + c, k);
            gl(6 * k + c) = eps * Bk(3 + c, k);
            gl(6 * k + 3 + c) = eps * Bk(6 + c, k);
            nmis = std::max(nmis, std::abs(eps * nk[c] - n[c]));
         }
      }
      Vector own;
      E.RestrictToOwnedFault(loc, own, FrameMap::kPerDOF);
      frames_ = FrameMap(own);
      // R8: expand the owner's eps_f t_k back to every local copy.
      Vector gown, gback;
      E.RestrictToOwnedFault(gl, gown, 6);
      E.ExpandOwnedToLocalFault(gown, gback, 6);
      real_t gmis = 0.0;
      for (int i = 0; i < gl.Size(); i++) { gmis = std::max(gmis, std::abs(gback(i) - gl(i))); }
      real_t red[2] = {gmis, nmis};
      MPI_Allreduce(MPI_IN_PLACE, red, 2, MPITypeMap<real_t>::mpi_type, MPI_MAX,
                    pmesh_.GetComm());
      ghost_mismatch_ = red[0];
      normal_mismatch_ = red[1];
   }

   ParMesh &pmesh_;
   const CycleConfig &cfg_;
   std::shared_ptr<real_t> t0_;   ///< held by the wall-load wrapper too
   std::unique_ptr<MaterialField> material_;
   BoundaryConfig bc_;
   std::unique_ptr<DieterichRuinaFriction> fr_;
   std::unique_ptr<RateStateClosure> closure_;
   std::unique_ptr<NitscheWaveOperator> W_;
   std::unique_ptr<CentralDifferenceSolver> ode_;
   FrameMap frames_;
   Vector U_, Ud_;
   real_t tau_ = 0.0;
   real_t dt_ = 0.0, dt_cfl_ = 0.0;
   long steps_ = 0;
   real_t ghost_mismatch_ = 0.0, normal_mismatch_ = 0.0;
};

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_DYN_PHASE_HPP

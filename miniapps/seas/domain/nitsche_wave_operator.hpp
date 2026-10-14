// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// nitsche_wave_operator.hpp — Nitsche dynamic rupture, plan item N-01.
//
// P3 content: the fault-free explicit elastodynamic operator
//
//    M a = F(t) - K u,    a = d2u/dt2,
//
// P4: optional first-order absorbing (dashpot) boundary on
// BoundaryConfig::absorbing_attrs, adding a velocity-only damping C:
//
//    M a + C v + K u = F(t).
//
// With central difference (gamma = 1/2) the step solves
// (M + dt/2 C) A^{n+1} = F(t_{n+1}) - K U^{n+1} - C U*  (U* = predicted
// velocity), so that M A^{n+1} + C Udot^{n+1} = F - K U^{n+1} holds exactly
// with Udot^{n+1} = U* + dt/2 A^{n+1}.  Without absorbing faces nothing
// changes (bit-identical to P3).
//
// as an MFEM SecondOrderTimeDependentOperator, to be advanced by MFEM's
// CentralDifferenceSolver (Newmark beta = 0, gamma = 1/2; see
// solver/nitsche_cfl.hpp).  K is the IP stiffness of ElasticityDomainOperator
// built through the Coefficient constructor (C-09) from the per-element
// NitscheMaterialSample; M^{-1} is NitscheMass (C-02) from the same sample.
//
// P6 (build sequence rev 5.2; plan rev 8.2 "The Nitsche pipeline"): the
// fault.  When the mesh carries fault faces the operator needs a
// NitscheFaultSetup and the body of ImplicitSolve is boxes 3-13 of the
// time-indexed algorithm:
//
//    3  d*     = d^n + dt/2 V^n                      half-step slip
//    4  psi*   = psi^n + dt r(|V^n|, psi^n)           predictor
//    5  d*_k   = ToKernel(d*)                         dictionary, per DOF
//    6  T_k    = ComputeTractionDiagnostics(U^{n+1}, d*_k)   (unchanged call)
//    7  RestrictToOwnedFault                          (unchanged call)
//    8  T      = FromKernel(T_k);  S = S_0 + G(t_{n+1}) S_nuc  (same projection)
//       X_tan  = P_tan(T + S), sigma_eff = -(T + S).n;  range check
//    9  V^{n+1}: closure with gamma_F = 1/2 sigma_F dt and psi*
//   10  psi^{n+1} = psi^n + dt/2 [r(|V^n|, psi^n) + r(|V^{n+1}|, psi*)]
//   11  d^{n+1} = d* + dt/2 V^{n+1}
//   12  G d^{n+1} = AssembleSlipOnlyRHS(ToKernel(d^{n+1})), F(t_{n+1}), K U^{n+1}
//   13  (M + dt/2 C) A^{n+1} = F - K U^{n+1} + G d^{n+1} - C Udot*
//
// The fault state (d, V global three-vectors, psi a scalar per DOF; spec A3)
// lives on the OWNED fault DOFs; ExpandOwnedToLocalFault /
// RestrictToOwnedFault bracket the kernel calls as in solver/seas_operator.hpp.
// The dictionary crossing (eps_f, the kernel triad) is applied per LOCAL DOF,
// because on a shared face each rank's kernel triad follows its own element.
// The traction is the numerical flux {sigma n} - sigma_F([[u]] - delta) of
// ComputeTractionDiagnostics; do not replace it by the plain stress average
// (Uphoff et al. section 5.1: unstable).
//
// Without fault faces boxes 3-11 are absent and every bit is P3/P4's.
//
// NitscheMaterialSample lives in nitsche_material_sample.hpp (moved there in
// P3 to break an include cycle with nitsche_mass.hpp); it is re-exported here.
//
// Scope (Chunhui Zhao, 2026-10-13): the Nitsche path uses tetrahedral meshes only.

#ifndef MFEM_SEAS_NITSCHE_WAVE_OPERATOR_HPP
#define MFEM_SEAS_NITSCHE_WAVE_OPERATOR_HPP

#include "mfem.hpp"
#include "nitsche_material_sample.hpp"
#include "nitsche_mass.hpp"
#include "elasticity_operator.hpp"
#include "../integrator/dg_absorbing_bc_integrator.hpp"
#include "../fault/fault_orientation.hpp"
#include "../fault/fault_frame.hpp"
#include "../friction/fault_closure.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>

namespace mfem
{
namespace seas
{

#ifdef MFEM_USE_MPI

/// Everything the fault path needs besides the mesh and the material.
/// Stresses are global tensors in the code frame (spec I1, I2; plan N-07),
/// tension positive, row-major 3x3.
struct NitscheFaultSetup
{
   const IFaultClosure *closure = nullptr;    ///< N-03; must outlive the operator
   FaultOrientationConfig orientation;        ///< N-08
   real_t up[3] = {0.0, 0.0, 1.0};            ///< N-09 frame (F1)
   /// Prestress sigma_0(x).
   std::function<void(const real_t x[3], real_t sigma[9])> prestress;
   /// Nucleation dsigma(x, t) = nucleation_spatial(x) nucleation_temporal(t)
   /// (TPV102: dtau F(r) (e_X e_Y + e_Y e_X) times G(t)).  Empty = none.
   std::function<void(const real_t x[3], real_t sigma[9])> nucleation_spatial;
   std::function<real_t(real_t t)> nucleation_temporal;
   /// Per-point friction parameters (a, Dc) at a fault DOF position.
   std::function<void(const real_t x[3], FaultPointState &s)> parameters;
   /// Permitted range sigma_n_eff >= sigma_min ([pore_pressure]
   /// min_sigma_n_pa, default 0; spatial/code/spatial_friction.cpp:957-960).
   real_t sigma_min = 0.0;
};

/// Per-step fault diagnostics (global over ranks).
struct NitscheStepRecord
{
   long step = 0;
   real_t t = 0.0;
   real_t min_sigma_eff = 0.0;
   real_t max_sigma_eff = 0.0;
   real_t max_V = 0.0;
   real_t max_V_dt_over_Dc = 0.0;   ///< predictor stability: abort above 2
};

/// Explicit elastodynamics d2u/dt2 = M^{-1} (F(t) - K u - C v + G d), with
/// the rate-and-state fault advanced inside ImplicitSolve when present.
class NitscheWaveOperator : public SecondOrderTimeDependentOperator
{
public:
   /// Facewise-constant tolerance of the scalar gamma (plan, "The scalar
   /// gamma"): a floating-point tolerance for an exact identity.
   static constexpr real_t kPenaltySpreadTol = 1e-10;

   /// `material` is sampled once per element; K and M share that sample.
   /// A mesh with fault faces (bc.fault_attr) requires `fault`; a mesh
   /// without them must not be given one.
   NitscheWaveOperator(ParMesh &mesh, int order, const MaterialField &material,
                       const BoundaryConfig &bc,
                       const DomainConfig &config = {},
                       const NitscheFaultSetup *fault = nullptr)
      : SecondOrderTimeDependentOperator(0, 0.0),
        mesh_(mesh), order_(order)
   {
      sample_ = std::make_unique<NitscheMaterialSample>(mesh, material);

      // Scalar fallbacks for the operator's scalar-only accessors: global
      // element means of the sample (not used by the IP assembly).
      real_t sums[3] = {0.0, 0.0, 0.0};
      for (int i = 0; i < sample_->LambdaGF().Size(); ++i)
      {
         sums[0] += sample_->LambdaGF()(i);
         sums[1] += sample_->MuGF()(i);
      }
      sums[2] = mesh.GetNE();
      MPI_Allreduce(MPI_IN_PLACE, sums, 3, MPITypeMap<real_t>::mpi_type,
                    MPI_SUM, mesh.GetComm());
      const real_t lam_fb = sums[0] / sums[2];
      const real_t mu_fb  = sums[1] / sums[2];

      elast_ = std::make_unique<ElasticityDomainOperator<ParMesh>>(
         mesh, order, sample_->Lambda(), sample_->Mu(), lam_fb, mu_fb,
         /*Vp=*/0.0, /*Wf=*/0.0, /*lf=*/0.0, bc, DGMethod::IP,
         SolverType::CG_AMG, config);

      int nfd = elast_->GetNumFaultDOFs();
      MPI_Allreduce(MPI_IN_PLACE, &nfd, 1, MPI_INT, MPI_SUM, mesh.GetComm());
      has_fault_ = (nfd > 0);
      MFEM_VERIFY(!has_fault_ || fault != nullptr,
                  "NitscheWaveOperator: the mesh carries " << nfd
                  << " fault DOFs (attr " << bc.fault_attr << ") but no "
                  "NitscheFaultSetup was given.");
      MFEM_VERIFY(has_fault_ || fault == nullptr,
                  "NitscheWaveOperator: a NitscheFaultSetup was given but the "
                  "mesh carries no fault faces (attr " << bc.fault_attr << ").");

      RhoPerElement(*sample_, rho_e_);
      mass_ = std::make_unique<NitscheMass>(mesh, order, rho_e_);
      BuildDampingBlocks_();

      const int n = elast_->GetFESpace().GetVSize();
      MFEM_VERIFY(mass_->VSize() == n,
                  "NitscheWaveOperator: mass size " << mass_->VSize()
                  << " != stiffness vector size " << n);
      height = width = n;
      Kx_.SetSize(n);
      F_.SetSize(n);

      if (has_fault_) { SetupFault_(*fault); }
   }

   /// Time step used to check the Newmark coefficients in ImplicitSolve.
   /// With absorbing faces it also builds (M + dt/2 C)^{-1}.
   void SetTimeStep(real_t dt)
   {
      dt_ = dt;
      if (!damping_blocks_.empty())
      {
         mass_damped_ = std::make_unique<NitscheMass>(
            mesh_, order_, rho_e_, &damping_blocks_, 0.5 * dt);
      }
   }
   real_t GetTimeStep() const { return dt_; }

   /// First-order form (Operator::Mult) is not used: second-order solvers
   /// call the three-argument Mult and ImplicitSolve below.
   void Mult(const Vector &, Vector &) const override
   {
      MFEM_ABORT("NitscheWaveOperator is second order; use "
                 "Mult(x, dxdt, y) / ImplicitSolve with a "
                 "SecondOrderODESolver.");
   }

   /// y = M^{-1}(F(t) - K x + G d - C dxdt) at the STORED fault state and
   /// t = GetTime().  Called once by NewmarkSolver for the initial
   /// acceleration A^0; advances nothing (T-21(a)).
   void Mult(const Vector &x, const Vector &dxdt, Vector &y) const override
   {
      Acceleration_(x, &dxdt, has_fault_ ? &d_ : nullptr, *mass_, y);
   }

   /// Central difference (beta = 0, gamma = 1/2): fac0 = beta dt^2 = 0 and
   /// fac1 = gamma dt = dt/2.  x = U^{n+1}, dxdt = Udot*, GetTime() = t_{n+1}
   /// (linalg/operator.hpp:808-810).  With a fault, boxes 3-13 advance
   /// (d, V, psi) in place; exactly one call per step (T-21).
   void ImplicitSolve(const real_t fac0, const real_t fac1,
                      const Vector &x, const Vector &dxdt, Vector &k) override
   {
      MFEM_VERIFY(fac0 == 0.0,
                  "NitscheWaveOperator::ImplicitSolve: fac0 = " << fac0
                  << " != 0; only the central-difference scheme "
                  "(Newmark beta = 0) is supported.");
      MFEM_VERIFY(fac1 == 0.5 * dt_,
                  "NitscheWaveOperator::ImplicitSolve: fac1 = " << fac1
                  << " != dt/2 = " << 0.5 * dt_
                  << "; call SetTimeStep(dt) with the solver's dt "
                  "(Newmark gamma = 1/2).");
      const real_t t = GetTime();
      MFEM_VERIFY(t != last_solve_time_,
                  "NitscheWaveOperator::ImplicitSolve called twice at t = "
                  << t << "; exactly one call per time step is allowed.");
      last_solve_time_ = t;
      if (has_fault_) { AdvanceFault_(x, t); }
      if (damping_blocks_.empty())
      {
         Acceleration_(x, nullptr, has_fault_ ? &d_ : nullptr, *mass_, k);
      }
      else
      {
         MFEM_VERIFY(mass_damped_ && mass_damped_->HalfDt() == fac1,
                     "NitscheWaveOperator::ImplicitSolve: the damped mass "
                     "(M + dt/2 C) was built for dt/2 = "
                     << (mass_damped_ ? mass_damped_->HalfDt() : 0.0)
                     << " but fac1 = " << fac1 << "; call SetTimeStep(dt).");
         Acceleration_(x, &dxdt, has_fault_ ? &d_ : nullptr, *mass_damped_, k);
      }
      if (has_fault_) { last_accel_ = k; }
   }

   const NitscheMaterialSample &Sample() const { return *sample_; }
   ParMesh &GetMesh() const { return mesh_; }
   int GetOrder() const { return order_; }
   ElasticityDomainOperator<ParMesh> &Elasticity() { return *elast_; }
   const NitscheMass &Mass() const { return *mass_; }
   /// (M + dt/2 C)^{-1}; null without absorbing faces or before SetTimeStep.
   const NitscheMass *DampedMass() const { return mass_damped_.get(); }
   /// Per-element dashpot blocks C_e (element-local byNODES, 3 nd x 3 nd).
   const std::map<int, DenseMatrix> &DampingBlocks() const { return damping_blocks_; }
   bool HasDamping() const { return !damping_blocks_.empty(); }

   /// Cv = C v (zero without absorbing faces).
   void ApplyDamping(const Vector &v, Vector &Cv) const
   {
      Cv.SetSize(v.Size());
      Cv = 0.0;
      if (damping_blocks_.empty()) { return; }
      const int n = mass_->NumScalarDofs();
      Vector x, y;
      for (const auto &kv : damping_blocks_)
      {
         const Array<int> &dofs = mass_->ElementScalarDofs(kv.first);
         const int nd = dofs.Size();
         x.SetSize(3 * nd);
         y.SetSize(3 * nd);
         for (int c = 0; c < 3; ++c)
            for (int k = 0; k < nd; ++k) { x(c * nd + k) = v(c * n + dofs[k]); }
         kv.second.Mult(x, y);
         for (int c = 0; c < 3; ++c)
            for (int k = 0; k < nd; ++k) { Cv(c * n + dofs[k]) += y(c * nd + k); }
      }
   }

   // ------------------------------------------------------------------
   // Fault (P6).  All per-DOF vectors are on the OWNED fault DOFs: 3
   // components (global x, y, z) for d, V, positions, tractions; 9 for the
   // frame (n, d, s); 1 for psi, a, Dc, sigma_F.
   // ------------------------------------------------------------------
   bool HasFault() const { return has_fault_; }
   int NumOwnedFaultDOFs() const { return no_; }

   /// Initial (or handed-over) fault state.  Sizes 3 no, 3 no, no.
   void SetFaultState(const Vector &d, const Vector &V, const Vector &psi)
   {
      MFEM_VERIFY(has_fault_, "SetFaultState: the operator has no fault");
      MFEM_VERIFY(d.Size() == 3 * no_ && V.Size() == 3 * no_ && psi.Size() == no_,
                  "SetFaultState: sizes " << d.Size() << ", " << V.Size()
                  << ", " << psi.Size() << " != " << 3 * no_ << ", "
                  << 3 * no_ << ", " << no_);
      d_ = d;
      V_ = V;
      psi_ = psi;
   }
   const Vector &Slip() const { return d_; }          ///< d,   3 no
   const Vector &SlipRate() const { return V_; }      ///< V,   3 no
   const Vector &State() const { return psi_; }       ///< psi, no
   const Vector &FaultPositions() const { return x_own_; }   ///< 3 no
   const Vector &FaultFrames() const { return frame_own_; }  ///< 9 no: n, d, s
   const Vector &FaultSigmaF() const { return sigF_own_; }   ///< no
   const Vector &FaultA() const { return a_own_; }
   const Vector &FaultDc() const { return Dc_own_; }
   const FaultOrientation &Orientation() const { return *orient_; }
   const NitscheStepRecord &LastStep() const { return record_; }
   /// Box 8/9 values of the last advance on owned DOFs: X_tan in the frame
   /// (dip, strike; 2 no) and sigma_n_eff (no).
   const Vector &LastShear() const { return Xlast_; }
   const Vector &LastSigmaEff() const { return siglast_; }
   /// Box 6-8 trial traction T_h(U^{n+1}, d*) of the last advance (3 no).
   const Vector &LastTrialTraction() const { return Town_; }
   /// A^{n+1} returned by the last ImplicitSolve (energy diagnostics).
   const Vector &LastAcceleration() const { return last_accel_; }
   /// Projected prestress (3 no) and nucleation spatial traction (3 no);
   /// the nucleation enters as NucleationFactor(t) * NucleationTraction().
   const Vector &PrestressTraction() const { return S0_own_; }
   const Vector &NucleationTraction() const { return Snuc_own_; }
   real_t NucleationFactor(real_t t) const { return nuc_temporal_ ? nuc_temporal_(t) : 0.0; }
   const IFaultClosure &Closure() const { return *closure_; }
   real_t SigmaMin() const { return sigma_min_; }
   /// Global number of shared (cross-rank) fault faces, each counted once.
   int NumSharedFaultFacesGlobal() const { return n_shared_global_; }
   /// Fault-face penalty mean and relative spread, in local face order.
   const Vector &FacePenalty() const { return sigF_face_; }
   const Vector &FacePenaltySpread() const { return spread_face_; }

   /// Projected prestress plus nucleation traction at time t (box 8).
   void BackgroundTraction(real_t t, Vector &S) const
   {
      S = S0_own_;
      if (nuc_temporal_)
      {
         const real_t g = nuc_temporal_(t);
         if (g != 0.0) { S.Add(g, Snuc_own_); }
      }
   }

   /// Boxes 5-8 (traction part): T_h(x, d) through the dictionary, as global
   /// three-vectors on owned DOFs.  `d` is on owned DOFs.
   void TrialTraction(const Vector &x, const Vector &d, Vector &T)
   {
      ToKernelLocal_(d, slip_k_);
      u_gf_->Vector::operator=(x);
      elast_->ComputeTractionDiagnostics(*u_gf_, slip_k_, Tk_, Tks_, Tkc_, jr_k_,
                                         &snk_);
      Tloc_.SetSize(3 * nd_);
      for (int k = 0; k < nd_; k++)
      {
         real_t nk[3], t1k[3], t2k[3];
         KernelTriad_(k, nk, t1k, t2k);
         fault_frame::FromKernel(Tk_(2 * k), Tk_(2 * k + 1), snk_(k), eps_loc_[k],
                                 nk, t1k, t2k, &Tloc_(3 * k));
      }
      elast_->RestrictToOwnedFault(Tloc_, T, 3);
   }

   /// Box 12, slip lifting: G d for d on owned DOFs.
   void SlipLoad(const Vector &d, Vector &Gd) const
   {
      ToKernelLocal_(d, slip_k_);
      elast_->AssembleSlipOnlyRHS(Gd, slip_k_);
   }

   /// Box 9 at one owned DOF j: X = (T + S) in the frame, the range check,
   /// the closure with gamma_F = 1/2 sigma_F dt and state psi.  Returns the
   /// global slip rate in V[3], sigma_eff and X_tan in the frame.
   void CloseDOF(int j, const real_t TS[3], real_t psi, real_t V[3],
                 real_t &sigma_eff, real_t X[2]) const
   {
      const real_t *n = &frame_own_(9 * j), *dv = &frame_own_(9 * j + 3),
                    *sv = &frame_own_(9 * j + 6);
      fault_frame::ToLocal(TS, 0.0, n, dv, sv, X[0], X[1], sigma_eff);
      if (!(sigma_eff >= sigma_min_))
      {
         MFEM_ABORT("NitscheWaveOperator: sigma_n_eff left its range on rank "
                    << Mpi::WorldRank() << ", owned fault DOF " << j
                    << " at (" << x_own_(3 * j) << ", " << x_own_(3 * j + 1)
                    << ", " << x_own_(3 * j + 2) << "), step " << steps_
                    << ", t = " << GetTime() << ": sigma_n_eff = " << sigma_eff
                    << " < sigma_min = " << sigma_min_);
      }
      FaultPointState st;
      st.psi = psi;
      st.a = a_own_(j);
      st.Dc = Dc_own_(j);
      real_t Vl[2];
      closure_->SolveSlipRate(X, 0.5 * sigF_own_(j) * dt_, sigma_eff, st, Vl);
      fault_frame::Embed(Vl[0], Vl[1], dv, sv, V);
   }

   /// r(|V|, psi) at owned DOF j (box 4, 10).
   real_t StateRate(int j, real_t Vabs, real_t psi) const
   {
      FaultPointState st;
      st.psi = psi;
      st.a = a_own_(j);
      st.Dc = Dc_own_(j);
      return closure_->StateRate(Vabs, st);
   }

   /// Boxes 12-13 for a given owned slip d (exposed for the T-21(c) hand
   /// loop): k = (M + fac C)^{-1}(F(GetTime()) - K x + G d - C v).
   void AccelerationWithSlip(const Vector &x, const Vector &v, const Vector &d,
                             Vector &k) const
   {
      if (damping_blocks_.empty()) { Acceleration_(x, nullptr, &d, *mass_, k); }
      else { Acceleration_(x, &v, &d, *mass_damped_, k); }
   }

private:
   /// Shared by Mult and ImplicitSolve so both give identical bits:
   /// a = Minv (F(t) - K x + G d - C v) with t = GetTime().  Without
   /// absorbing faces the C v term is skipped entirely, without a fault the
   /// G d term (P3/P4 bits).
   void Acceleration_(const Vector &x, const Vector *v, const Vector *d,
                      const NitscheMass &Minv, Vector &a) const
   {
      elast_->ApplyStiffness(x, Kx_);
      elast_->AssembleDirichletOnlyRHS(F_, GetTime());
      F_ -= Kx_;
      if (d)
      {
         SlipLoad(*d, Gd_);
         F_ += Gd_;
      }
      if (!damping_blocks_.empty() && v)
      {
         ApplyDamping(*v, Cv_);
         F_ -= Cv_;
      }
      Minv.Apply(F_, a);
   }

   /// Boxes 3-11 at t = t_{n+1}, x = U^{n+1}.
   void AdvanceFault_(const Vector &x, real_t t)
   {
      ++steps_;
      const real_t dt = dt_;
      // 3, 4
      dstar_.SetSize(3 * no_);
      psistar_.SetSize(no_);
      r0_.SetSize(no_);
      for (int j = 0; j < no_; j++)
      {
         const real_t *V = &V_(3 * j);
         for (int c = 0; c < 3; c++) { dstar_(3 * j + c) = d_(3 * j + c) + 0.5 * dt * V[c]; }
         const real_t Vabs = std::sqrt(V[0] * V[0] + V[1] * V[1] + V[2] * V[2]);
         r0_(j) = StateRate(j, Vabs, psi_(j));
         psistar_(j) = psi_(j) + dt * r0_(j);
      }
      // 5-8
      TrialTraction(x, dstar_, Town_);
      BackgroundTraction(t, S_);
      // 9-11
      real_t min_sig = std::numeric_limits<real_t>::infinity();
      real_t max_sig = -std::numeric_limits<real_t>::infinity();
      real_t max_V = 0.0, max_ratio = 0.0;
      Xlast_.SetSize(2 * no_);
      siglast_.SetSize(no_);
      for (int j = 0; j < no_; j++)
      {
         real_t TS[3], Vn[3], X[2], sig;
         for (int c = 0; c < 3; c++) { TS[c] = Town_(3 * j + c) + S_(3 * j + c); }
         CloseDOF(j, TS, psistar_(j), Vn, sig, X);
         const real_t Vabs = std::sqrt(Vn[0] * Vn[0] + Vn[1] * Vn[1] + Vn[2] * Vn[2]);
         // C-07 anti-parallel: V . X_tan < 0 wherever |V| > 0.
         const real_t *dv = &frame_own_(9 * j + 3), *sv = &frame_own_(9 * j + 6);
         const real_t VX = X[0] * fault_frame::Dot(Vn, dv) + X[1] * fault_frame::Dot(Vn, sv);
         MFEM_VERIFY(Vabs == 0.0 || VX < 0.0,
                     "NitscheWaveOperator: slip rate not anti-parallel to the "
                     "shear traction on rank " << Mpi::WorldRank() << ", owned "
                     "fault DOF " << j << ", step " << steps_ << ": V.X = " << VX);
         const real_t ratio = Vabs * dt / Dc_own_(j);
         MFEM_VERIFY(ratio <= 2.0,
                     "NitscheWaveOperator: predictor stability |V| dt / Dc = "
                     << ratio << " > 2 on rank " << Mpi::WorldRank()
                     << ", owned fault DOF " << j << " at (" << x_own_(3 * j)
                     << ", " << x_own_(3 * j + 1) << ", " << x_own_(3 * j + 2)
                     << "), step " << steps_ << ", t = " << t);
         psi_(j) = psi_(j) + 0.5 * dt * (r0_(j) + StateRate(j, Vabs, psistar_(j)));
         for (int c = 0; c < 3; c++)
         {
            V_(3 * j + c) = Vn[c];
            d_(3 * j + c) = dstar_(3 * j + c) + 0.5 * dt * Vn[c];
         }
         Xlast_(2 * j) = X[0];
         Xlast_(2 * j + 1) = X[1];
         siglast_(j) = sig;
         min_sig = std::min(min_sig, sig);
         max_sig = std::max(max_sig, sig);
         max_V = std::max(max_V, Vabs);
         max_ratio = std::max(max_ratio, ratio);
      }
      real_t red[4] = {-min_sig, max_sig, max_V, max_ratio};
      MPI_Allreduce(MPI_IN_PLACE, red, 4, MPITypeMap<real_t>::mpi_type, MPI_MAX,
                    mesh_.GetComm());
      record_.step = steps_;
      record_.t = t;
      record_.min_sigma_eff = -red[0];
      record_.max_sigma_eff = red[1];
      record_.max_V = red[2];
      record_.max_V_dt_over_Dc = red[3];
   }

   /// Kernel triad (n_k, t1_k, t2_k) of local fault DOF k.
   void KernelTriad_(int k, real_t nk[3], real_t t1k[3], real_t t2k[3]) const
   {
      for (int c = 0; c < 3; c++)
      {
         nk[c] = Bk_(c, k);
         t1k[c] = Bk_(3 + c, k);
         t2k[c] = Bk_(6 + c, k);
      }
   }

   /// Owned global three-vectors -> local kernel (dip, strike) pairs.
   void ToKernelLocal_(const Vector &d, Vector &slip_k) const
   {
      elast_->ExpandOwnedToLocalFault(d, dloc_, 3);
      slip_k.SetSize(2 * nd_);
      for (int k = 0; k < nd_; k++)
      {
         real_t nk[3], t1k[3], t2k[3];
         KernelTriad_(k, nk, t1k, t2k);
         fault_frame::ToKernel(&dloc_(3 * k), eps_loc_[k], t1k, t2k,
                               slip_k(2 * k), slip_k(2 * k + 1));
      }
   }

   /// Face transformation of local fault face fi (interior, then shared).
   FaceElementTransformations *FaultFaceTr_(int fi) const
   {
      const int n_int = elast_->GetFaultInteriorFaces().Size();
      FaceElementTransformations *Tr =
         (fi < n_int) ? mesh_.GetInteriorFaceTransformations(
                           elast_->GetFaultInteriorFaces()[fi])
                      : mesh_.GetSharedFaceTransformations(
                           elast_->GetFaultSharedFaces()[fi - n_int]);
      MFEM_VERIFY(Tr != nullptr, "NitscheWaveOperator: no transformation for "
                  "local fault face " << fi);
      return Tr;
   }

   /// The prestress / nucleation traction through the SAME road as the
   /// kernel's traction: sigma(x_q) n_q with n_q the face's unit normal
   /// CalcOrtho / nl_q (the normal ComputeTractionAtQuadPoints uses), at the
   /// same quadrature points, projected with ProjectTractionToFaultDOFs
   /// (same rule, basis, nl_q weights and triad), then FromKernel ->
   /// global three-vectors on owned DOFs, i.e. sigma nhat.  (Projecting
   /// sigma nhat directly and then applying FromKernel would apply eps_f
   /// twice: the first P6 run caught exactly that as sigma_n_eff = -120 MPa.)
   void ProjectStress_(const std::function<void(const real_t[3], real_t[9])> &sig,
                       Vector &S_own) const
   {
      const FaceQuadrature *fq = elast_->GetFaceQuadrature();
      const DenseMatrix &e_q = fq->BasisAtQuadPoints();
      const FaultBasis *fb = elast_->GetFaultBasis();
      Vector S_loc(3 * nd_);
      S_loc = 0.0;
      Vector nor(3), xq(3), Tq, nlq, out;
      for (int fi = 0; fi < nloc_faces_; fi++)
      {
         FaceElementTransformations *Tr = FaultFaceTr_(fi);
         const IntegrationRule &ir =
            IntRules.Get(Tr->GetGeometryType(), 2 * order_ + 1);
         const int nq = ir.GetNPoints();
         MFEM_VERIFY(e_q.Width() == nq && e_q.Height() == nbf_,
                     "NitscheWaveOperator: face basis at quadrature points is "
                     << e_q.Height() << " x " << e_q.Width() << ", expected "
                     << nbf_ << " x " << nq);
         Tq.SetSize(3 * nq);
         nlq.SetSize(nq);
         for (int q = 0; q < nq; q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            Tr->SetAllIntPoints(&ip);
            CalcOrtho(Tr->Jacobian(), nor);
            const real_t nl = nor.Norml2();
            nlq(q) = nl;
            const real_t nq3[3] = {nor(0) / nl, nor(1) / nl, nor(2) / nl};
            Tr->Transform(ip, xq);
            real_t s[9];
            sig(xq.GetData(), s);
            for (int a = 0; a < 3; a++)
            {
               Tq(a * nq + q) = s[3 * a] * nq3[0] + s[3 * a + 1] * nq3[1]
                                + s[3 * a + 2] * nq3[2];
            }
         }
         const FaultBasisData &b = fb->GetBasis(fi);
         real_t bv[3][3];
         for (int c = 0; c < 3; c++)
         {
            bv[0][c] = b.normal[c];
            bv[1][c] = b.tangent1[c];
            bv[2][c] = b.tangent2[c];
         }
         DGElasticityIPCombinedIntegrator::ProjectTractionToFaultDOFs(
            3, 3, Tq, nlq, ir, nbf_, e_q, bv, out,
            b.qp_data.empty() ? nullptr : &b.qp_data);
         for (int kk = 0; kk < nbf_; kk++)
         {
            const int k = fi * nbf_ + kk;
            real_t nk[3], t1k[3], t2k[3];
            KernelTriad_(k, nk, t1k, t2k);
            // sigma_n = -(projected normal component), as ComputeTractionImpl.
            fault_frame::FromKernel(out(nbf_ + kk), out(2 * nbf_ + kk), -out(kk),
                                    eps_loc_[k], nk, t1k, t2k, &S_loc(3 * k));
         }
      }
      elast_->RestrictToOwnedFault(S_loc, S_own, 3);
   }

   void SetupFault_(const NitscheFaultSetup &fs)
   {
      MFEM_VERIFY(fs.closure != nullptr, "NitscheFaultSetup: closure is null");
      MFEM_VERIFY(static_cast<bool>(fs.prestress), "NitscheFaultSetup: no prestress");
      MFEM_VERIFY(static_cast<bool>(fs.parameters), "NitscheFaultSetup: no parameters");
      MFEM_VERIFY(static_cast<bool>(fs.nucleation_spatial)
                  == static_cast<bool>(fs.nucleation_temporal),
                  "NitscheFaultSetup: nucleation needs both its spatial and "
                  "its temporal factor");
      closure_ = fs.closure;
      sigma_min_ = fs.sigma_min;
      nuc_temporal_ = fs.nucleation_temporal;
      nd_ = elast_->GetNumFaultDOFs();
      no_ = elast_->GetNumOwnedFaultDOFs();
      nbf_ = elast_->GetNbfPerFace();
      const Array<int> &fint = elast_->GetFaultInteriorFaces();
      const Array<int> &fsh = elast_->GetFaultSharedFaces();
      nloc_faces_ = fint.Size() + fsh.Size();
      MFEM_VERIFY(nd_ == nloc_faces_ * nbf_, "NitscheWaveOperator: " << nd_
                  << " local fault DOFs != " << nloc_faces_ << " faces x " << nbf_);
      int nsh = fsh.Size();
      MPI_Allreduce(MPI_IN_PLACE, &nsh, 1, MPI_INT, MPI_SUM, mesh_.GetComm());
      n_shared_global_ = nsh / 2;

      // N-08 orientation; its constructor refuses an inadmissible fault.
      orient_ = std::make_unique<FaultOrientation>(mesh_, fint, fsh, fs.orientation);

      // Kernel triad and eps_f per local DOF (dictionary, plan Part IV).
      elast_->GetFaultDOFBasis(Bk_);
      eps_loc_.assign(nd_, 1.0);
      for (int k = 0; k < nd_; k++)
      {
         real_t nk[3], t1k[3], t2k[3];
         KernelTriad_(k, nk, t1k, t2k);
         eps_loc_[k] = fault_frame::EpsilonF(orient_->Normal(k / nbf_), nk);
      }

      // N-09 frame per local face -> per owned DOF.
      Vector frame_loc(9 * nd_);
      for (int fi = 0; fi < nloc_faces_; fi++)
      {
         fault_frame::FaceVertices fv;
         fv.nv = 3;
         orient_->FaceVertexData(fi, fv.key, fv.x);
         fault_frame::Frame fr;
         const bool ok = fault_frame::ComputeFrame(
            orient_->Normal(fi), fs.up, fault_frame::DegeneracyThreshold(1.0), fv, fr);
         MFEM_VERIFY(ok, "NitscheWaveOperator: degenerate frame (F6(b)) on "
                     "local fault face " << fi);
         for (int kk = 0; kk < nbf_; kk++)
         {
            real_t *f = &frame_loc(9 * (fi * nbf_ + kk));
            for (int c = 0; c < 3; c++) { f[c] = fr.n[c]; f[3 + c] = fr.d[c]; f[6 + c] = fr.s[c]; }
         }
      }
      elast_->RestrictToOwnedFault(frame_loc, frame_own_, 9);

      // C-04 sigma_F per face; refuse a face whose penalty is not constant.
      elast_->GetFaultFacePenalties(sigF_face_, spread_face_);
      for (int fi = 0; fi < nloc_faces_; fi++)
      {
         if (spread_face_(fi) > kPenaltySpreadTol)
         {
            long long key[3];
            real_t xv[3][3];
            orient_->FaceVertexData(fi, key, xv);
            MFEM_ABORT("NitscheWaveOperator: the penalty sigma_F is not "
                       "constant on fault face " << fi << " of rank "
                       << Mpi::WorldRank() << " (centroid ("
                       << (xv[0][0] + xv[1][0] + xv[2][0]) / 3 << ", "
                       << (xv[0][1] + xv[1][1] + xv[2][1]) / 3 << ", "
                       << (xv[0][2] + xv[1][2] + xv[2][2]) / 3
                       << ")): relative spread " << spread_face_(fi) << " > "
                       << kPenaltySpreadTol << "; the scalar gamma_F would not "
                       "be the kernel's derivative (plan, \"The scalar gamma\").");
         }
      }
      Vector sig_loc(nd_);
      for (int k = 0; k < nd_; k++) { sig_loc(k) = sigF_face_(k / nbf_); }
      elast_->RestrictToOwnedFault(sig_loc, sigF_own_, 1);

      // Positions and friction parameters per owned DOF.
      Vector x_loc;
      elast_->GetFaultDOFCoords3D(x_loc);
      elast_->RestrictToOwnedFault(x_loc, x_own_, 3);
      a_own_.SetSize(no_);
      Dc_own_.SetSize(no_);
      for (int j = 0; j < no_; j++)
      {
         FaultPointState st;
         fs.parameters(&x_own_(3 * j), st);
         MFEM_VERIFY(st.a > 0.0 && st.Dc > 0.0, "NitscheFaultSetup: parameters "
                     "give a = " << st.a << ", Dc = " << st.Dc << " at owned "
                     "fault DOF " << j);
         a_own_(j) = st.a;
         Dc_own_(j) = st.Dc;
      }

      // Projected prestress and nucleation (box 8).
      ProjectStress_(fs.prestress, S0_own_);
      if (fs.nucleation_spatial) { ProjectStress_(fs.nucleation_spatial, Snuc_own_); }
      else { Snuc_own_.SetSize(3 * no_); Snuc_own_ = 0.0; }

      d_.SetSize(3 * no_);   d_ = 0.0;
      V_.SetSize(3 * no_);   V_ = 0.0;
      psi_.SetSize(no_);     psi_ = 0.0;
      u_gf_ = std::make_unique<ParGridFunction>(&elast_->GetFESpace());
      Gd_.SetSize(height);
   }

   /// C_e for every element with an absorbing face (N-06), material from the
   /// sample (not the operator's fallback scalars).
   void BuildDampingBlocks_()
   {
      damping_blocks_.clear();
      const auto &fes = elast_->GetFESpace();
      const ParFiniteElementSpace *mfes = sample_->RhoGF().ParFESpace();
      Array<int> mdofs;
      for (int be : elast_->GetAbsorbingBdrElements())
      {
         FaceElementTransformations *Tr = mesh_.GetBdrFaceTransformations(be);
         MFEM_VERIFY(Tr != nullptr, "NitscheWaveOperator: no face "
                     "transformation for absorbing boundary element " << be);
         const int e = Tr->Elem1No;
         const FiniteElement &fe = *fes.GetFE(e);
         mfes->GetElementDofs(e, mdofs);
         const real_t lam = sample_->LambdaGF()(mdofs[0]);
         const real_t mu  = sample_->MuGF()(mdofs[0]);
         const real_t rho = sample_->RhoGF()(mdofs[0]);
         const real_t cp = std::sqrt((lam + 2.0 * mu) / rho);
         const real_t cs = std::sqrt(mu / rho);
         const int nd = fe.GetDof();
         auto it = damping_blocks_.find(e);
         if (it == damping_blocks_.end())
         {
            it = damping_blocks_.emplace(e, DenseMatrix(3 * nd)).first;
            it->second = 0.0;
         }
         DGAbsorbingBCIntegrator::AddFaceBlock(fe, *Tr, rho, cp, cs, it->second);
      }
   }

   ParMesh &mesh_;
   int order_;
   std::unique_ptr<NitscheMaterialSample> sample_;
   std::unique_ptr<ElasticityDomainOperator<ParMesh>> elast_;
   std::unique_ptr<NitscheMass> mass_;
   std::unique_ptr<NitscheMass> mass_damped_;
   std::map<int, DenseMatrix> damping_blocks_;
   Vector rho_e_;
   real_t dt_ = 0.0;
   real_t last_solve_time_ = std::numeric_limits<real_t>::quiet_NaN();
   mutable Vector Kx_, F_, Cv_;

   // Fault (P6)
   bool has_fault_ = false;
   const IFaultClosure *closure_ = nullptr;
   real_t sigma_min_ = 0.0;
   std::function<real_t(real_t)> nuc_temporal_;
   std::unique_ptr<FaultOrientation> orient_;
   int nd_ = 0, no_ = 0, nbf_ = 1, nloc_faces_ = 0, n_shared_global_ = 0;
   DenseMatrix Bk_;                    // kernel triad per local DOF (9 x nd)
   std::vector<real_t> eps_loc_;       // eps_f per local DOF
   Vector frame_own_, x_own_, sigF_own_, a_own_, Dc_own_;
   Vector sigF_face_, spread_face_;
   Vector S0_own_, Snuc_own_;
   Vector d_, V_, psi_;
   long steps_ = 0;
   NitscheStepRecord record_;
   std::unique_ptr<ParGridFunction> u_gf_;
   Vector dstar_, psistar_, r0_, Town_, S_, Xlast_, siglast_;
   Vector Tk_, Tks_, Tkc_, jr_k_, snk_, Tloc_;
   mutable Vector dloc_, slip_k_, Gd_;
   Vector last_accel_;
};

#endif // MFEM_USE_MPI

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_NITSCHE_WAVE_OPERATOR_HPP

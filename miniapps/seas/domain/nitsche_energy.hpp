// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// nitsche_energy.hpp — Nitsche dynamic rupture, plan "The quantity velocity
// Verlet conserves" and T-15: the time-discrete energy balance of the
// Nitsche operator with the fault.
//
//    E_h = 1/2 Udot.M Udot + E_el(U, d) - dt^2/8 A.M A
//    E_el = 1/2 U.K U - U.G d + 1/2 int_F sigma_F |d|^2
//         = (SIPG volume + skeleton) - int_F {sigma n}.([[u]] - d)
//           + 1/2 int_F sigma_F |[[u]] - d|^2
//
// and the stepwise balance E_h^{n+1} - E_h^n = dt Pbar^{n+1/2} + rho^n with
// Pbar the trapezoidal average of
//
//    P_ext  = Udot.F                            (external load)
//    P_abs  = -Udot.C Udot                      (dashpot, <= 0)
//    P_pre  = -int_F S_0 . V                    (prestress)
//    P_nuc  = -G(t) int_F S_nuc . V             (nucleation)
//    P_fric =  int_F (T_h + S) . V              (friction, <= 0: E4)
//    P_kv   = -(v.K_eta v - 2 v.G_eta V + int_F sigma_F^eta |V|^2)
//                                               (Kelvin-Voigt layer, <= 0)
//
// With the Kelvin-Voigt layer the friction power uses the total traction
// T_h + T_eta(v, V) + S that the closure sees (P_fric stays <= 0 node by
// node), and the remaining damping power is P_kv = -2 E_el^eta(v, V).
//
// since dE_el/dd = T_h (the kernel's numerical-flux traction), the fault
// power is int_F T_h . V = P_fric + P_pre + P_nuc.  No eta_s V^2 term (E6).
// Fault integrals use the face mass M_F = sum_q w_q nl_q phi_i phi_j of the
// traction projection, on owned faces only.  Dot products of DOF fields are
// frame-invariant, so the kernel's own components are used for the parts.
//
// An observer: it reads the operator through accessors and changes nothing
// in the step.  Scope (Chunhui Zhao, 2026-10-14).

#ifndef MFEM_SEAS_NITSCHE_ENERGY_HPP
#define MFEM_SEAS_NITSCHE_ENERGY_HPP

#include "mfem.hpp"
#include "nitsche_wave_operator.hpp"

#include <memory>
#include <vector>

namespace mfem
{
namespace seas
{

#ifdef MFEM_USE_MPI

struct NitscheEnergyRecord
{
   real_t t = 0.0;
   real_t E_kin = 0.0, E_el = 0.0, E_mod = 0.0, E_h = 0.0;
   real_t E_sipg = 0.0, E_cons = 0.0, E_pen = 0.0;   ///< parts of E_el
   real_t P_ext = 0.0, P_abs = 0.0, P_pre = 0.0, P_nuc = 0.0, P_fric = 0.0;
   real_t P_kv = 0.0;   ///< Kelvin-Voigt layer (0 without it)
   /// E4 at the nodes: max over fault DOFs of (T_h + S).V (<= 0 up to
   /// round-off).  The integral P_fric of the two face interpolants is not
   /// sign-definite at p >= 2 (the P2 vertex functions integrate to zero),
   /// so E4 is a nodal statement.
   real_t fric_node_max = 0.0;
   real_t fric_node_scale = 0.0;   ///< max |T_h + S| |V| over DOFs
   real_t P() const { return P_ext + P_abs + P_pre + P_nuc + P_fric + P_kv; }
};

class NitscheEnergyMonitor
{
public:
   /// drop_consistency: T-15 negative control (E_el without the fault
   /// consistency term).
   explicit NitscheEnergyMonitor(NitscheWaveOperator &op, bool drop_consistency = false)
      : op_(op), E_(op.Elasticity()), drop_cons_(drop_consistency)
   {
      MFEM_VERIFY(op.HasFault(), "NitscheEnergyMonitor: the operator has no fault");
      comm_ = op.GetMesh().GetComm();
      mass_ = std::make_unique<ParBilinearForm>(&E_.GetFESpace());
      mass_->AddDomainIntegrator(new VectorMassIntegrator(
         const_cast<NitscheMaterialSample &>(op.Sample()).Rho()));
      mass_->Assemble();
      mass_->Finalize();
      BuildFaceMass_();
   }

   /// Record at the current state.  U, Udot, A at time t; the fault state
   /// is the operator's.  The fault traction is evaluated at (U, d).
   NitscheEnergyRecord Evaluate(const Vector &U, const Vector &Ud, const Vector &A,
                                real_t t)
   {
      NitscheEnergyRecord r;
      r.t = t;
      const real_t dt = op_.GetTimeStep();
      Vector MUd(Ud.Size()), MA(A.Size()), KU, Gd, CUd, F;
      mass_->SpMat().Mult(Ud, MUd);
      mass_->SpMat().Mult(A, MA);
      E_.ApplyStiffness(U, KU);
      op_.SlipLoad(op_.Slip(), Gd);
      op_.ApplyDamping(Ud, CUd);
      E_.AssembleDirichletOnlyRHS(F, t);
      const Vector &d = op_.Slip(), &V = op_.SlipRate();
      r.E_kin = 0.5 * Dot_(Ud, MUd);
      const real_t pen_dd = 0.5 * FaultDot_(d, d, true);
      r.E_el = 0.5 * Dot_(U, KU) - Dot_(U, Gd) + pen_dd;
      r.E_mod = dt * dt / 8.0 * Dot_(A, MA);
      // parts at quadrature level through the kernel's decomposition
      Parts_(U, d, r.E_cons, r.E_pen);
      r.E_sipg = r.E_el - r.E_cons - r.E_pen;
      if (drop_cons_) { r.E_el -= r.E_cons; }
      r.E_h = r.E_kin + r.E_el - r.E_mod;
      r.P_ext = Dot_(Ud, F);
      r.P_abs = -Dot_(Ud, CUd);
      // fault power
      Vector T, S0 = op_.PrestressTraction(), Sn = op_.NucleationTraction();
      op_.TrialTraction(U, d, T);
      const real_t g = op_.NucleationFactor(t);
      Vector TS(T);
      TS += S0;
      if (g != 0.0) { TS.Add(g, Sn); }
      if (op_.HasKelvinVoigt())
      {
         Vector KvU, GvV, Tkv;
         op_.ApplyKVStiffness(Ud, KvU);
         op_.KVSlipLoad(V, GvV);
         op_.KVTrialTraction(Ud, V, Tkv);
         TS += Tkv;
         r.P_kv = -(Dot_(Ud, KvU) - 2.0 * Dot_(Ud, GvV)
                    + FaultDot_(V, V, true, &op_.FacePenaltyKV()));
      }
      r.P_fric = FaultDot_(TS, V, false);
      real_t nm[2] = {-1e300, 0.0};
      for (int j = 0; j < op_.NumOwnedFaultDOFs(); j++)
      {
         real_t dp = 0.0, a2 = 0.0, b2 = 0.0;
         for (int c = 0; c < 3; c++)
         {
            dp += TS(3 * j + c) * V(3 * j + c);
            a2 += TS(3 * j + c) * TS(3 * j + c);
            b2 += V(3 * j + c) * V(3 * j + c);
         }
         nm[0] = std::max(nm[0], dp);
         nm[1] = std::max(nm[1], std::sqrt(a2 * b2));
      }
      MPI_Allreduce(MPI_IN_PLACE, nm, 2, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm_);
      r.fric_node_max = nm[0];
      r.fric_node_scale = nm[1];
      r.P_pre = -FaultDot_(S0, V, false);
      r.P_nuc = -g * FaultDot_(Sn, V, false);
      return r;
   }

   /// Balance bookkeeping: call Start once, then Step after every step.
   void Start(const NitscheEnergyRecord &r0) { prev_ = r0; sum_rho_ = 0.0; sum_abs_ = 0.0; n_ = 0; }
   real_t Step(const NitscheEnergyRecord &r)
   {
      const real_t dt = r.t - prev_.t;
      const real_t rho = (r.E_h - prev_.E_h) - 0.5 * dt * (r.P() + prev_.P());
      sum_rho_ += rho;
      sum_abs_ += std::abs(rho);
      min_fric_ = std::min({min_fric_, r.P_fric, prev_.P_fric});
      max_fric_ = std::max({max_fric_, r.P_fric, prev_.P_fric});
      node_fric_max_ = std::max(node_fric_max_, r.fric_node_max);
      node_fric_scale_ = std::max(node_fric_scale_, r.fric_node_scale);
      prev_ = r;
      n_++;
      return rho;
   }
   real_t SumRho() const { return sum_rho_; }
   real_t SumAbsRho() const { return sum_abs_; }
   real_t MaxFrictionPower() const { return max_fric_; }
   real_t MinFrictionPower() const { return min_fric_; }
   /// max over the run and the fault DOFs of (T_h + S).V, and of |T_h+S||V|
   real_t NodalFrictionMax() const { return node_fric_max_; }
   real_t NodalFrictionScale() const { return node_fric_scale_; }

private:
   real_t Dot_(const Vector &a, const Vector &b) const
   {
      return InnerProduct(comm_, a, b);
   }

   /// int_F a . b over owned fault faces, a and b on owned DOFs (3 comps).
   /// sigma = true weights by sigma_F of the face.
   real_t FaultDot_(const Vector &a, const Vector &b, bool sigma,
                    const Vector *face_weight = nullptr) const
   {
      Vector al, bl;
      E_.ExpandOwnedToLocalFault(a, al, 3);
      E_.ExpandOwnedToLocalFault(b, bl, 3);
      const Array<int> &own = E_.GetOwnedFaultFaceMap();
      const int nbf = E_.GetNbfPerFace();
      real_t s = 0.0;
      for (int of = 0; of < own.Size(); of++)
      {
         const int lf = own[of];
         const DenseMatrix &M = MF_[lf];
         const real_t w = sigma ? (face_weight ? (*face_weight)(lf) : op_.FacePenalty()(lf))
                                : 1.0;
         real_t fs = 0.0;
         for (int i = 0; i < nbf; i++)
            for (int j = 0; j < nbf; j++)
               for (int c = 0; c < 3; c++)
               {
                  fs += al(3 * (lf * nbf + i) + c) * M(i, j) * bl(3 * (lf * nbf + j) + c);
               }
         s += w * fs;
      }
      MPI_Allreduce(MPI_IN_PLACE, &s, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm_);
      return s;
   }

   /// Fault consistency and penalty energies from the kernel's
   /// decomposition at (U, d): corr = -sigma_F R, R = [[u]] - d, so
   /// -int {sigma n}.R = (1/sigma_F) int Tst.corr and
   /// 1/2 int sigma_F |R|^2 = 1/(2 sigma_F) int |corr|^2 (both exact: R and
   /// corr lie in the face space, so their L2 projections are themselves).
   void Parts_(const Vector &U, const Vector &d, real_t &E_cons, real_t &E_pen)
   {
      const int nd = E_.GetNumFaultDOFs(), nbf = E_.GetNbfPerFace();
      // the operator's own local kernel slip for d
      Vector dl;
      E_.ExpandOwnedToLocalFault(d, dl, 3);
      DenseMatrix B;
      E_.GetFaultDOFBasis(B);
      Vector slip(2 * nd);
      for (int k = 0; k < nd; k++)
      {
         real_t nk[3], t1k[3], t2k[3];
         for (int c = 0; c < 3; c++) { nk[c] = B(c, k); t1k[c] = B(3 + c, k); t2k[c] = B(6 + c, k); }
         const real_t e = fault_frame::EpsilonF(op_.Orientation().Normal(k / nbf), nk);
         fault_frame::ToKernel(&dl(3 * k), e, t1k, t2k, slip(2 * k), slip(2 * k + 1));
      }
      ParGridFunction u(&E_.GetFESpace());
      u = U;
      Vector T, Ts, Tc, R, sn, sns, snc;
      E_.ComputeTractionDiagnostics(u, slip, T, Ts, Tc, R, &sn, &sns, &snc);
      const Array<int> &own = E_.GetOwnedFaultFaceMap();
      real_t s[2] = {0.0, 0.0};
      for (int of = 0; of < own.Size(); of++)
      {
         const int lf = own[of];
         const DenseMatrix &M = MF_[lf];
         const real_t sF = op_.FacePenalty()(lf);
         real_t c1 = 0.0, c2 = 0.0;
         for (int i = 0; i < nbf; i++)
            for (int j = 0; j < nbf; j++)
            {
               const int a = lf * nbf + i, b = lf * nbf + j;
               const real_t m = M(i, j);
               c1 += m * (Ts(2 * a) * Tc(2 * b) + Ts(2 * a + 1) * Tc(2 * b + 1) + sns(a) * snc(b));
               c2 += m * (Tc(2 * a) * Tc(2 * b) + Tc(2 * a + 1) * Tc(2 * b + 1) + snc(a) * snc(b));
            }
         s[0] += c1 / sF;
         s[1] += 0.5 * c2 / sF;
      }
      MPI_Allreduce(MPI_IN_PLACE, s, 2, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm_);
      E_cons = s[0];
      E_pen = s[1];
   }

   /// M_F per local fault face on the traction rule (nl-weighted).
   void BuildFaceMass_()
   {
      const int nbf = E_.GetNbfPerFace();
      const DenseMatrix &e_q = E_.GetFaceQuadrature()->BasisAtQuadPoints();
      ParMesh &pm = op_.GetMesh();
      const Array<int> &fint = E_.GetFaultInteriorFaces();
      const Array<int> &fsh = E_.GetFaultSharedFaces();
      const int nloc = fint.Size() + fsh.Size();
      MF_.assign(nloc, DenseMatrix(nbf));
      Vector nor(3);
      for (int f = 0; f < nloc; f++)
      {
         FaceElementTransformations *Tr = (f < fint.Size())
            ? pm.GetInteriorFaceTransformations(fint[f])
            : pm.GetSharedFaceTransformations(fsh[f - fint.Size()]);
         const IntegrationRule &ir = IntRules.Get(Tr->GetGeometryType(),
                                                  2 * op_.GetOrder() + 1);
         DenseMatrix &M = MF_[f];
         M = 0.0;
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            Tr->SetAllIntPoints(&ir.IntPoint(q));
            CalcOrtho(Tr->Jacobian(), nor);
            const real_t wn = ir.IntPoint(q).weight * nor.Norml2();
            for (int i = 0; i < nbf; i++)
               for (int j = 0; j < nbf; j++) { M(i, j) += wn * e_q(i, q) * e_q(j, q); }
         }
      }
   }

   NitscheWaveOperator &op_;
   ElasticityDomainOperator<ParMesh> &E_;
   bool drop_cons_;
   MPI_Comm comm_;
   std::unique_ptr<ParBilinearForm> mass_;
   std::vector<DenseMatrix> MF_;
   NitscheEnergyRecord prev_;
   real_t sum_rho_ = 0.0, sum_abs_ = 0.0;
   real_t max_fric_ = -1e300, min_fric_ = 1e300;
   real_t node_fric_max_ = -1e300, node_fric_scale_ = 0.0;
   long n_ = 0;
};

#endif // MFEM_USE_MPI

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_NITSCHE_ENERGY_HPP

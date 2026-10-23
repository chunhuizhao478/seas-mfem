// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// poro_mass.hpp — poroelastic option, step P2.4 (plan rev 3).
//
// PoroMass applies the inverse of the matrix of the acceleration solve (9.2)
// of derivation rev 1.2,
//
//    a = (M + dt D Gamma)^{-1} F,        x = [u; w],  F = [F_u; F_w],
//
// with the mass and drag of (8.2): on element e
//
//    M|_e = [rho  rho_f] (x) M_e (x) I_3,   D|_e = [0 0] (x) M_e (x) I_3 + C_e,
//           [rho_f   m ]                           [0 b]
//
// M_e the scalar DG mass of NitscheMass and C_e the elastic operator's
// absorbing dashpot (u rows only; NitscheWaveOperator::DampingBlocks()).
// Gamma is the velocity coefficient of (9.5): its w row is [0, gamma_w], with
// gamma_w per element (FittedGamma, step P2.5); its u row is [1/2, g_uw],
// g_uw = 0 (Gamma_diag, "diagonal") or (rho_f / rho)(1/2 - gamma_w) (Gamma_mom,
// "momentum", option F3).  So with m' = m + dt gamma_w b
//
//    (M + dt D Gamma)|_e = [ rho M_e + dt/2 C_e    rho_f M_e + dt g_uw C_e ]
//                          [ rho_f M_e             m' M_e                  ]
//
// (each block (x) I_3, C_e already 3 nd x 3 nd).  PoroMass eliminates w, the
// Schur form of plan rev 3 section 5.1 (step P2.4).  With r = rho_f / m' and
// y_w = M_e^{-1} F_w:
//
//    S_u = (rho - r rho_f) M_e + dt (1/2 - g_uw r) C_e    (the Schur complement)
//    a_u = S_u^{-1} (F_u - r F_w - (dt g_uw / m') C_e y_w)
//    a_w = y_w / m' - r a_u
//
// Away from absorbing faces (C_e = 0) this is (8.3): a_u = W_uu M_e^{-1}
// (F_u - r F_w), W_uu = 1 / (rho - rho_f^2 / m'), the same for both Gamma
// (Proposition 9.3).  S_u has the form of NitscheMass's own block, rho M_e +
// half_dt C_e, with rho_eff = rho - r rho_f and half_dt_eff = dt (1/2 - g_uw r)
// per element, so PoroMass composes one NitscheMass built from those: its
// Apply gives a_u on every element (the scalar path, or the dense inverse of
// S_u where C_e is present), and its ElementMassInverse(e) is M_e^{-1}.  With
// rho_f = 0, r = 0 and g_uw = 0, so rho_eff = rho, half_dt_eff = dt * 0.5 and
// the u rows are today's NitscheMass(rho, C, dt/2) bit for bit (test U2-11).
// rho_eff > 0 because rho m - rho_f^2 > 0 and m' >= m; half_dt_eff >= dt/2
// because g_uw <= 0, so S_u stays SPD.
//
// Input per local element: the Biot constants (rho, rho_f, m, b; the sample's,
// or any with rho m - rho_f^2 > 0, rho_f = 0 included) and gamma_w in
// [1/2, 1].  dt and Gamma are fixed at construction: the stepper builds a new
// PoroMass when dt changes.  Layout: x = [u; w], each in the elastic
// operator's vector DG layout (byNODES, n = 3 ndofs).  Element-local, no
// communication; the constructor is collective (NitscheMass builds a
// ParFiniteElementSpace).
//
// Scope (Chunhui Zhao, 2026-10-13): tetrahedral meshes only.

#ifndef MFEM_SEAS_PORO_MASS_HPP
#define MFEM_SEAS_PORO_MASS_HPP

#include "mfem.hpp"
#include "../../domain/nitsche_mass.hpp"
#include "poro_material.hpp"

#include <map>
#include <memory>
#include <vector>

namespace mfem
{
namespace seas
{
namespace poro
{

/// The u row of Gamma (9.5): [multiphysics.poroelastic] velocity_coupling.
enum class VelocityCoupling
{
   kMomentum,   ///< Gamma_mom (F3), "momentum"
   kDiagonal    ///< Gamma_diag, "diagonal"
};

/// g_uw of (9.5): the (u, w) entry of Gamma.
inline real_t GammaUW(VelocityCoupling coupling, real_t rho, real_t rho_f, real_t gamma_w)
{
   return coupling == VelocityCoupling::kMomentum ? rho_f / rho * (0.5 - gamma_w) : 0.0;
}

#ifdef MFEM_USE_MPI

/// The sample's constants of every local element, for PoroMass.
inline std::vector<BiotConstants> ElementConstants(const PoroMaterialSample &sample)
{
   const int ne = sample.AlphaGF().ParFESpace()->GetNE();
   std::vector<BiotConstants> c(ne);
   for (int e = 0; e < ne; e++) { c[e] = sample.Constants(e); }
   return c;
}

class PoroMass
{
public:
   /// `c` and `gamma_w` per local element of `mesh`; `damping_blocks` the
   /// elastic operator's dashpot blocks C_e (NitscheWaveOperator::
   /// DampingBlocks(); kept, not owned), or null without absorbing faces.
   PoroMass(ParMesh &mesh, int order, const std::vector<BiotConstants> &c,
            real_t dt, const Vector &gamma_w, VelocityCoupling coupling,
            const std::map<int, DenseMatrix> *damping_blocks = nullptr)
      : dt_(dt), coupling_(coupling), gamma_w_(gamma_w)
   {
      const int ne = mesh.GetNE();
      MFEM_VERIFY(dt > 0.0, "PoroMass: dt (" << dt << ") must be > 0.");
      MFEM_VERIFY(static_cast<int>(c.size()) == ne && gamma_w.Size() == ne,
                  "PoroMass: " << c.size() << " constants and " << gamma_w.Size()
                  << " values of gamma_w for " << ne << " local elements.");
      Vector rho_eff(ne), half_dt_eff(ne);
      r_.SetSize(ne);
      inv_mp_.SetSize(ne);
      s_.SetSize(ne);
      for (int e = 0; e < ne; e++)
      {
         const BiotConstants &k = c[e];
         const real_t g = gamma_w(e);
         MFEM_VERIFY(g >= 0.5 && g <= 1.0, "PoroMass: gamma_w (" << g
                     << ") at element " << e << " is not in [1/2, 1] (derivation "
                     "(9.4)).");
         MFEM_VERIFY(k.rho > 0.0 && k.m > 0.0 && k.rho_f >= 0.0 && k.b >= 0.0
                     && k.rho * k.m - k.rho_f * k.rho_f > 0.0,
                     "PoroMass: element " << e << " has rho " << k.rho << ", rho_f "
                     << k.rho_f << ", m " << k.m << ", b " << k.b
                     << "; the mass needs rho, m > 0, rho_f, b >= 0 and "
                     "rho m - rho_f^2 > 0.");
         const real_t mp = k.m + dt * g * k.b;
         const real_t g_uw = GammaUW(coupling, k.rho, k.rho_f, g);
         r_(e) = k.rho_f / mp;
         inv_mp_(e) = 1.0 / mp;
         s_(e) = dt * g_uw / mp;
         rho_eff(e) = k.rho - r_(e) * k.rho_f;
         half_dt_eff(e) = dt * (0.5 - g_uw * r_(e));
      }
      static const std::map<int, DenseMatrix> kNoBlocks;
      blocks_ = damping_blocks ? damping_blocks : &kNoBlocks;
      inner_ = std::make_unique<NitscheMass>(mesh, order, rho_eff, blocks_, 0.0,
                                             &half_dt_eff);
      n_ = inner_->VSize();
      g_.SetSize(n_);
      yw_.SetSize(n_);
   }

   PoroMass(const PoroMass &) = delete;
   PoroMass &operator=(const PoroMass &) = delete;

   /// a = (M + dt D Gamma)^{-1} F, F and a of size 2 n.  Operation order per
   /// element: y_w = M_e^{-1} F_w; g = F_u - r F_w, minus s C_e y_w (s =
   /// dt g_uw / m') where C_e is present and s != 0; a_u = NitscheMass(rho_eff,
   /// half_dt_eff).Apply(g); a_w = y_w (1/m') - r a_u.
   void Apply(const Vector &F, Vector &a) const
   {
      MFEM_VERIFY(F.Size() == 2 * n_, "PoroMass::Apply: input size " << F.Size()
                  << " != 2 n = " << 2 * n_);
      a.SetSize(2 * n_);
      const int ndofs = inner_->NumScalarDofs();
      const int ne = r_.Size();
      Vector x, y, z, cz;
      for (int e = 0; e < ne; e++)
      {
         const Array<int> &dofs = inner_->ElementScalarDofs(e);
         const DenseMatrix &Minv = inner_->ElementMassInverse(e);
         const int nd = dofs.Size();
         const real_t r = r_(e);
         x.SetSize(nd);
         y.SetSize(nd);
         for (int c = 0; c < 3; c++)
         {
            const int off = c * ndofs;
            for (int k = 0; k < nd; k++) { x(k) = F(n_ + off + dofs[k]); }
            Minv.Mult(x, y);
            for (int k = 0; k < nd; k++)
            {
               const int i = off + dofs[k];
               yw_(i) = y(k);
               g_(i) = F(i) - r * F(n_ + i);
            }
         }
         if (s_(e) != 0.0)
         {
            auto it = blocks_->find(e);
            if (it != blocks_->end())
            {
               z.SetSize(3 * nd);
               cz.SetSize(3 * nd);
               for (int c = 0; c < 3; c++)
                  for (int k = 0; k < nd; k++) { z(c * nd + k) = yw_(c * ndofs + dofs[k]); }
               it->second.Mult(z, cz);
               for (int c = 0; c < 3; c++)
                  for (int k = 0; k < nd; k++) { g_(c * ndofs + dofs[k]) -= s_(e) * cz(c * nd + k); }
            }
         }
      }
      Vector au;
      au.MakeRef(a, 0, n_);
      inner_->Apply(g_, au);
      for (int e = 0; e < ne; e++)
      {
         const Array<int> &dofs = inner_->ElementScalarDofs(e);
         const real_t r = r_(e), inv_mp = inv_mp_(e);
         for (int c = 0; c < 3; c++)
         {
            const int off = c * ndofs;
            for (int k = 0; k < dofs.Size(); k++)
            {
               const int i = off + dofs[k];
               a(n_ + i) = yw_(i) * inv_mp - r * a(i);
            }
         }
      }
   }

   /// Size of one field (u or w); F and a have 2 n.
   int NumU() const { return n_; }
   real_t Dt() const { return dt_; }
   VelocityCoupling Coupling() const { return coupling_; }
   const Vector &GammaW() const { return gamma_w_; }
   /// The NitscheMass of S_u (rho_eff, half_dt_eff); its ElementMassInverse
   /// is M_e^{-1}.
   const NitscheMass &SchurMass() const { return *inner_; }

private:
   real_t dt_;
   VelocityCoupling coupling_;
   Vector gamma_w_;
   Vector r_, inv_mp_, s_;    ///< rho_f / m', 1 / m', dt g_uw / m' per element
   const std::map<int, DenseMatrix> *blocks_ = nullptr;
   std::unique_ptr<NitscheMass> inner_;
   int n_ = 0;
   mutable Vector g_, yw_;
};

#endif // MFEM_USE_MPI

} // namespace poro
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_PORO_MASS_HPP

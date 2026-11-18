// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// width_memory.hpp — the 2.5D width memory of viscoelastic attenuation plan
// rev6 (section 5.6, Result 6, eqs. (2.41b), (2.42b), (5.4), (5.5)).
//
// In the 2.5D finite-width model the volume term -beta (u - u_pl), beta =
// b mu_R, is the shear stress of the rock outside the slab.  With the outside
// rock as viscoelastic as the slab (assumption A6), the term gets its own
// memory, one displacement vector per mechanism:
//
//    eta_l' + omega_l eta_l = u' - u_pl'                                 (5.4)
//    f_eta = - sum_l b Y^S_l mu_U M_1 eta_l                              (5.5)
//
// M_1 the mass matrix with unit density (b mu_R M_1 is the M_beta of
// ElasticityDomainOperator::SetFiniteWidth).  eta lives on the displacement
// DOFs (it is a displacement), element by element; no face terms.  The
// step is the exact exponential one of the strain memory (5.3), applied to
// the increment of w = u - u_pl, u_pl(t) = (t_off + t) U_pl1:
//
//    eta_l^{n+1} = exp(-omega_l dt) eta_l^n + phi(omega_l dt) (w^{n+1} - w^n)
//
// The coefficients b Y^S_l mu_U are per element, from the fit at the
// element's material (as the strain memory's lambda_l, mu_l).  Elements
// excluded from the rheology (the PML band, decision D4) are excluded here
// too: their eta stays 0.  Reset(U, t) gives the entry state eta = 0,
// w^n = U - u_pl(t) (2.42d).
//
// Local time stepping (plan section 5.8, phase V6): the memory has no face
// terms, so each element steps at its own cluster's ticks with its dt_c
// (AdvanceLts) and its force goes into its own cluster's (AddLoadCluster).
//
// Scope (Chunhui Zhao, 2026-11-18).

#ifndef MFEM_SEAS_WIDTH_MEMORY_HPP
#define MFEM_SEAS_WIDTH_MEMORY_HPP

#include "mfem.hpp"
#include "viscoelastic.hpp"
#include "../domain/nitsche_material_sample.hpp"

#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <vector>

namespace mfem
{
namespace seas
{

#ifdef MFEM_USE_MPI

class WidthMemory
{
public:
   /// ufes: the displacement space (2-D, byNODES); beta_over_mu: b of
   /// [cycle.finite_width] (> 0); U_pl1: the substrate displacement at
   /// physical time 1 on the vdofs of ufes; t_offset: the physical time of
   /// the operator's t = 0 (read at every call; null = 0).  The fit is the
   /// model's, at each element's material.
   WidthMemory(ParFiniteElementSpace &ufes, const NitscheMaterialSample &sample,
               const ViscoelasticModel &model, real_t beta_over_mu, const Vector &U_pl1,
               std::shared_ptr<const real_t> t_offset)
      : ufes_(ufes), L_(model.L()), omega_(model.Omega()), U_pl1_(U_pl1),
        t_off_(std::move(t_offset))
   {
      ParMesh &pm = *ufes.GetParMesh();
      MFEM_VERIFY(std::isfinite(beta_over_mu) && beta_over_mu > 0.0,
                  "WidthMemory: b = " << beta_over_mu << " must be > 0 (without the 2.5D "
                  "term there is no width memory; plan section 5.6)");
      MFEM_VERIFY(pm.Dimension() == 2 && ufes.GetVDim() == 2
                  && ufes.GetOrdering() == Ordering::byNODES,
                  "WidthMemory: a 2-D displacement space with 2 components, byNODES");
      MFEM_VERIFY(U_pl1.Size() == ufes.GetVSize(), "WidthMemory: U_pl1 has " << U_pl1.Size()
                  << " entries, the space " << ufes.GetVSize());
      ne_ = pm.GetNE();
      const int n = ufes.GetVSize();
      eta_.assign(static_cast<std::size_t>(L_) * n, 0.0);
      w_old_.SetSize(n);
      w_old_ = 0.0;
      coef_.assign(static_cast<std::size_t>(ne_) * L_, 0.0);
      excluded_.assign(ne_, 0);
      Me_.resize(ne_);
      vdofs_.resize(ne_);
      MassIntegrator mi;
      for (int e = 0; e < ne_; e++)
      {
         ElementTransformation *T = pm.GetElementTransformation(e);
         mi.AssembleElementMatrix(*ufes.GetFE(e), *T, Me_[e]);
         ufes.GetElementVDofs(e, vdofs_[e]);
         const IntegrationPoint &c = Geometries.GetCenter(pm.GetElementBaseGeometry(e));
         const real_t lam = sample.LambdaGF().GetValue(e, c);
         const real_t mu = sample.MuGF().GetValue(e, c);
         const real_t rho = sample.RhoGF().GetValue(e, c);
         const GMBFit f = model.Fit(lam, mu, rho);
         MFEM_VERIFY(f.error.empty(), "WidthMemory: element " << e << ": " << f.error);
         for (int l = 0; l < L_; l++)
         {
            coef_[static_cast<std::size_t>(e) * L_ + l] = beta_over_mu * f.YS[l] * f.mu_U;
         }
      }
   }

   int L() const { return L_; }
   /// b Y^S_l mu_U of element e.
   real_t Coef(int e, int l) const { return coef_[static_cast<std::size_t>(e) * L_ + l]; }
   const Vector &SubstrateShape() const { return U_pl1_; }
   real_t TimeOffset() const { return t_off_ ? *t_off_ : 0.0; }

   /// Elements without memory (1 = excluded), as RheologyCoupler::SetExcluded.
   void SetExcluded(const std::vector<char> &mask)
   {
      MFEM_VERIFY(static_cast<int>(mask.size()) == ne_, "WidthMemory::SetExcluded: "
                  << mask.size() << " flags for " << ne_ << " elements");
      excluded_ = mask;
   }

   /// Entry state (2.42d) at the operator's time t.
   void Reset(const Vector &U, real_t t)
   {
      std::fill(eta_.begin(), eta_.end(), 0.0);
      W_(U, t, w_old_);
      last_t_ = t;
   }

   /// Step to t = t_{n+1} with U = U^{n+1}.
   void Advance(const Vector &U, real_t t, real_t dt)
   {
      MFEM_VERIFY(t > last_t_, "WidthMemory::Advance: t = " << t
                  << " does not follow the last step's t = " << last_t_);
      last_t_ = t;
      if (dt != dt_cached_) { Factors_(dt); }
      W_(U, t, w_new_);
      const int n = ufes_.GetVSize();
      for (int e = 0; e < ne_; e++)
      {
         if (excluded_[e]) { continue; }
         const Array<int> &vd = vdofs_[e];
         for (int a = 0; a < vd.Size(); a++)
         {
            const int k = vd[a];
            const real_t dw = w_new_(k) - w_old_(k);
            for (int l = 0; l < L_; l++)
            {
               real_t &h = eta_[static_cast<std::size_t>(l) * n + k];
               h = a_[l] * h + b_[l] * dw;
            }
         }
      }
      w_old_ = w_new_;
   }

   /// Local time stepping: clusters of the local elements.
   void SetLts(const std::vector<int> &cl)
   {
      MFEM_VERIFY(static_cast<int>(cl.size()) == ne_, "WidthMemory::SetLts: " << cl.size()
                  << " clusters for " << ne_ << " elements");
      cl_ = cl;
   }

   /// One tick (clusters 0..top end at t, U = u(t)): the elements of those
   /// clusters step with dt(c) on their own DOFs.
   void AdvanceLts(int top, const Vector &U, real_t t, const std::function<real_t(int)> &dt)
   {
      MFEM_VERIFY(!cl_.empty(), "WidthMemory::AdvanceLts: SetLts first");
      W_(U, t, w_new_);
      const int n = ufes_.GetVSize();
      real_t ah[16], bh[16];
      MFEM_VERIFY(L_ <= 16, "WidthMemory: L = " << L_ << " > 16");
      int c_last = -1;
      for (int e = 0; e < ne_; e++)
      {
         const int c = cl_[e];
         if (c > top) { continue; }
         const Array<int> &vd = vdofs_[e];
         if (!excluded_[e])
         {
            if (c != c_last)
            {
               const real_t h = dt(c);
               for (int l = 0; l < L_; l++)
               {
                  const real_t x = omega_[l] * h;
                  ah[l] = std::exp(-x);
                  bh[l] = (x > 0.0) ? -std::expm1(-x) / x : 1.0;
               }
               c_last = c;
            }
            for (int a = 0; a < vd.Size(); a++)
            {
               const int k = vd[a];
               const real_t dw = w_new_(k) - w_old_(k);
               for (int l = 0; l < L_; l++)
               {
                  real_t &h = eta_[static_cast<std::size_t>(l) * n + k];
                  h = ah[l] * h + bh[l] * dw;
               }
            }
         }
         for (int a = 0; a < vd.Size(); a++) { w_old_(vd[a]) = w_new_(vd[a]); }
      }
   }

   /// F += f_eta of the elements of cluster c.
   void AddLoadCluster(int c, Vector &F) const
   {
      MFEM_VERIFY(!cl_.empty(), "WidthMemory::AddLoadCluster: SetLts first");
      AddLoad_(F, c);
   }

   /// F += f_eta = - sum_l b Y^S_l mu_U M_1 eta_l (5.5).
   void AddLoad(Vector &F) const { AddLoad_(F, -1); }

private:
   void AddLoad_(Vector &F, int only) const
   {
      const int n = ufes_.GetVSize();
      Vector x, y;
      for (int e = 0; e < ne_; e++)
      {
         if (excluded_[e] || (only >= 0 && cl_[e] != only)) { continue; }
         const DenseMatrix &M = Me_[e];
         const int nd = M.Height();
         const Array<int> &vd = vdofs_[e];
         x.SetSize(nd);
         y.SetSize(nd);
         for (int c = 0; c < 2; c++)
         {
            x = 0.0;
            for (int l = 0; l < L_; l++)
            {
               const real_t cl = Coef(e, l);
               const real_t *h = &eta_[static_cast<std::size_t>(l) * n];
               for (int j = 0; j < nd; j++) { x(j) += cl * h[vd[c * nd + j]]; }
            }
            M.Mult(x, y);
            for (int j = 0; j < nd; j++) { F(vd[c * nd + j]) -= y(j); }
         }
      }
   }

public:

   /// 1/2 sum_l eta_l . (b Y^S_l mu_U M_1) eta_l and the dissipation rate
   /// sum_l omega_l eta_l . (b Y^S_l mu_U M_1) eta_l, this rank.
   void Energy(real_t &E, real_t &D) const
   {
      const int n = ufes_.GetVSize();
      E = D = 0.0;
      Vector x, y;
      for (int e = 0; e < ne_; e++)
      {
         const DenseMatrix &M = Me_[e];
         const int nd = M.Height();
         const Array<int> &vd = vdofs_[e];
         x.SetSize(nd);
         y.SetSize(nd);
         for (int l = 0; l < L_; l++)
         {
            const real_t *h = &eta_[static_cast<std::size_t>(l) * n];
            for (int c = 0; c < 2; c++)
            {
               for (int j = 0; j < nd; j++) { x(j) = h[vd[c * nd + j]]; }
               M.Mult(x, y);
               const real_t q = Coef(e, l) * (x * y);
               E += 0.5 * q;
               D += omega_[l] * q;
            }
         }
      }
   }

   /// Largest |eta_l| component on this rank (the hand-back report).
   real_t MaxEta() const
   {
      real_t m = 0.0;
      for (real_t h : eta_) { m = std::max(m, std::abs(h)); }
      return m;
   }
   /// eta_l (vdofs); the writable form is for tests that start from a
   /// state other than the entry state (after Reset).
   const real_t *Eta(int l) const { return &eta_[static_cast<std::size_t>(l) * ufes_.GetVSize()]; }
   real_t *Eta(int l) { return &eta_[static_cast<std::size_t>(l) * ufes_.GetVSize()]; }

private:
   void W_(const Vector &U, real_t t, Vector &w) const
   {
      MFEM_VERIFY(U.Size() == U_pl1_.Size(), "WidthMemory: U has " << U.Size()
                  << " entries, U_pl1 " << U_pl1_.Size());
      w = U;
      w.Add(-(TimeOffset() + t), U_pl1_);
   }
   void Factors_(real_t dt)
   {
      a_.resize(L_);
      b_.resize(L_);
      for (int l = 0; l < L_; l++)
      {
         const real_t x = omega_[l] * dt;
         a_[l] = std::exp(-x);
         b_[l] = (x > 0.0) ? -std::expm1(-x) / x : 1.0;
      }
      dt_cached_ = dt;
   }

   ParFiniteElementSpace &ufes_;
   int L_ = 0, ne_ = 0;
   std::vector<real_t> omega_;
   Vector U_pl1_;
   std::shared_ptr<const real_t> t_off_;
   std::vector<real_t> eta_;              ///< L blocks of vsize
   Vector w_old_, w_new_;
   std::vector<real_t> coef_;             ///< b Y^S_l mu_U per element
   std::vector<char> excluded_;
   std::vector<DenseMatrix> Me_;          ///< unit-density scalar mass per element
   std::vector<Array<int>> vdofs_;
   std::vector<real_t> a_, b_;
   std::vector<int> cl_;                  ///< LTS clusters (empty: global stepping)
   real_t dt_cached_ = -1.0;
   real_t last_t_ = -std::numeric_limits<real_t>::infinity();
};

#endif // MFEM_USE_MPI

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_WIDTH_MEMORY_HPP

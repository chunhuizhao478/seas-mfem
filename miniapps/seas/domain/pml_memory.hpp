// domain/pml_memory.hpp - the dynamic PML's memory terms (PML design plan
// rev2, Sections 3.1-3.7; implementation step P2, Section 7.3).
//
// The dynamic phase starts at the QD state u_s: U = u_s + w (3.1).  With
// K_p the SIPG matrix with C and the shared penalty kappa_gamma P (Section
// 3.5), the internal force is (3.12)
//   F_int(U) = K_p U + b_band + B(w, Phi)
//   b_band   = (K0 - K_p) u_s, with the wall data at t_s: the band terms of
//              the zero-frequency tensor C0 minus those of C, on u_s
//   B        = the same band terms with the time-domain band tensor C_pml
//              acting on w (and on the wall jump w - (g(t) - g(t_s))),
//              minus those of C: the memory part
// and the mass side adds C_d w' + E(w, chi) (3.9)-(3.10).  In the inner
// domain everything here is zero.
//
// Band terms ("Delta sigma" = C_pml * G - C : G for an input G):
//   volume       + int  Delta sigma(grad w) : grad v
//   consistency  - int  {Delta sigma(grad w)} n . [v]          (per side material)
//   symmetry     + eps/2 int {grad v} : Delta sigma([w] x n)    (eps = -1, SIPG;
//                eps on a wall face, with [w] = w - (g(t) - g(t_s)))
// The penalty is not filtered: it is all in K_p.
//
// C_pml * G at one point, three stages (Section 3.3):
//   1. G~ = G - sum_a phi_a x n_a,      phi_a' + (alpha + d_a) phi_a = d_a G n_a
//   2. sigma~ = lam tr(G~) I + mu (G~ + G~^T)
//   3. sigma = sum_b [ (S / s_b) * (sigma~ n_b) ] x n_b over an orthonormal
//      frame {n_b} holding the active normals; S / s_b is the cascade of
//      filter A (y = z + m, m' + alpha m = d_a z) over the active a != b.
// Mass side (3.7)-(3.9): y_0 = w, y_a = y_(a-1) + chi_a,
//   chi_a' + alpha chi_a = d_a y_(a-1);  rho y_J'' = rho (w'' + sum_a d_a w' + E_J).
// Every filter steps by (3.16):
//   m^(n+1) = e m^n + (d / q)(1 - e)(z^n + z^(n+1)) / 2,  e = exp(-q dt),
// q = alpha (A) or alpha + d (B), so a constant input settles to the
// zero-frequency gain exactly: b_band + B -> (K0 - K_p) U at rest
// (Proposition 3 holds for the discrete operators).
//
// 2.5D finite width (SetFiniteWidth; 2.5D with the PML, piece F-a): the term
// beta (u - u_pl) of rho u'' = div sigma - beta (u - u_pl) gets the mass-side
// stretch S.  K_p + M_beta and the substrate load already hold the plain
// term, so the band adds, at its volume points,
//   beta [ (S0 - 1) (u_s - u_pl(t)) + (y_J - w) ],   y_J - w = sum_a chi_a,
// S0 = prod_a (1 + d_a / alpha) the zero-frequency value of S and u_pl(t) =
// (t_off + t) U_pl1.  At rest y_J = S0 w, so the total is beta S0 (U - u_pl):
// the QD band's term (PmlBandBetaMassIntegrator).  No new state.
//
// Faces: interior and shared faces with a quadrature point inside a band,
// and Dirichlet wall faces inside a band.  Each side of a face keeps its own
// two pipelines (its material); a rank adds only its own side's test terms
// on a shared face.  Natural faces carry nothing (sigma_pml n = 0 weakly).
//
// Quadrature: volume and mass side 2p + N (the band integrator of P1, so
// K_p + B at rest equals the QD K0 to rounding); faces 2p + 1.

#ifndef MFEM_SEAS_PML_MEMORY_HPP
#define MFEM_SEAS_PML_MEMORY_HPP

#include "mfem.hpp"
#include "pml_band.hpp"
#include "elasticity_operator.hpp"
#include "../common/rank_costs.hpp"
#include "../integrator/pml_elasticity_integrators.hpp"

#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <vector>

namespace mfem
{
namespace seas
{

// The per-point kernels of the device form (UCERF 3-D plan rev5, P1.10;
// PmlDynamicTerms::SetDevice): StepPipe_, Stage12_, Hooke_, StepMass_ and
// MassE_ for D = 3 on flat arrays, the same operations in the same order.
// n: the point's frame (3 x 3, row b = n_b); c: its filter coefficients,
// kCoef reals: d[3], eB[3], gB[3], eA[3], gA[3].
namespace pmldev
{
constexpr int kCoef = 16;
constexpr int kMaxNd = 10;   ///< nodes per element (p <= 2 tetrahedra)

MFEM_HOST_DEVICE inline void Hooke(real_t lam, real_t mu, const real_t H[3][3], real_t s[3][3])
{
   real_t tr = 0.0;
   for (int k = 0; k < 3; k++) { tr += H[k][k]; }
   for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
      { s[i][j] = mu * (H[i][j] + H[j][i]) + (i == j ? lam * tr : 0.0); }
}

MFEM_HOST_DEVICE inline void Stage12(int na, const real_t *n, real_t lam, real_t mu,
                                     const real_t *st, real_t s2[3][3])
{
   const real_t *Gp = st, *phi = st + 9;
   real_t Gt[3][3];
   for (int i = 0; i < 3; i++)
      for (int l = 0; l < 3; l++)
      {
         real_t v = Gp[i * 3 + l];
         for (int k = 0; k < na; k++) { v -= phi[k * 3 + i] * n[3 * k + l]; }
         Gt[i][l] = v;
      }
   Hooke(lam, mu, Gt, s2);
}

MFEM_HOST_DEVICE inline void StepPipe(int na, const real_t *n, const real_t *c, real_t lam,
                                      real_t mu, const real_t G[3][3], real_t *st,
                                      real_t ds[3][3])
{
   const real_t *eB = c + 3, *gB = c + 6, *eA = c + 9, *gA = c + 12;
   real_t *Gp = st;
   real_t *phi = st + 9;
   for (int k = 0; k < na; k++)
   {
      const real_t *nk = n + 3 * k;
      for (int i = 0; i < 3; i++)
      {
         real_t zn = 0.0, zo = 0.0;
         for (int l = 0; l < 3; l++) { zn += G[i][l] * nk[l]; zo += Gp[i * 3 + l] * nk[l]; }
         phi[k * 3 + i] = eB[k] * phi[k * 3 + i] + gB[k] * (zo + zn);
      }
   }
   for (int i = 0; i < 9; i++) { Gp[i] = G[i / 3][i % 3]; }
   real_t st2[3][3];
   Stage12(na, n, lam, mu, st, st2);
   real_t *f = phi + na * 3;
   real_t sig[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
   for (int b = 0; b < 3; b++)
   {
      const real_t *nb = n + 3 * b;
      real_t z[3];
      for (int i = 0; i < 3; i++)
      {
         z[i] = 0.0;
         for (int l = 0; l < 3; l++) { z[i] += st2[i][l] * nb[l]; }
      }
      for (int k = 0; k < na; k++)
      {
         if (k == b) { continue; }
         real_t *m = f, *zo = f + 3;
         for (int i = 0; i < 3; i++)
         {
            m[i] = eA[k] * m[i] + gA[k] * (zo[i] + z[i]);
            zo[i] = z[i];
            z[i] += m[i];
         }
         f += 6;
      }
      for (int i = 0; i < 3; i++)
         for (int j = 0; j < 3; j++) { sig[i][j] += z[i] * nb[j]; }
   }
   real_t Gn[3][3], s0[3][3];
   for (int i = 0; i < 3; i++)
      for (int l = 0; l < 3; l++) { Gn[i][l] = Gp[i * 3 + l]; }
   Hooke(lam, mu, Gn, s0);
   for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) { ds[i][j] = sig[i][j] - s0[i][j]; }
}

MFEM_HOST_DEVICE inline void StepMass(int na, const real_t *c, const real_t *w, real_t *st)
{
   const real_t *eA = c + 9, *gA = c + 12;
   real_t y[3];
   for (int i = 0; i < 3; i++) { y[i] = w[i]; }
   real_t *f = st + 3;
   for (int k = 0; k < na; k++)
   {
      real_t *chi = f, *zo = f + 3;
      for (int i = 0; i < 3; i++)
      {
         chi[i] = eA[k] * chi[i] + gA[k] * (zo[i] + y[i]);
         zo[i] = y[i];
         y[i] += chi[i];
      }
      f += 6;
   }
   for (int i = 0; i < 3; i++) { st[i] = w[i]; }
}

MFEM_HOST_DEVICE inline void MassE(int na, const real_t *c, real_t al, const real_t *st,
                                   real_t E[3])
{
   real_t y[3], Yd[3];
   for (int i = 0; i < 3; i++) { y[i] = st[i]; Yd[i] = 0.0; E[i] = 0.0; }
   const real_t *f = st + 3;
   for (int k = 0; k < na; k++)
   {
      const real_t *chi = f;
      for (int i = 0; i < 3; i++)
      {
         const real_t dy = c[k] * y[i] - al * chi[i];
         E[i] += c[k] * Yd[i] - al * dy;
         Yd[i] += dy;
         y[i] += chi[i];
      }
      f += 6;
   }
}
}  // namespace pmldev

class PmlDynamicTerms
{
public:
   /// lam_e, mu_e, rho_e: element values (local elements) of the material;
   /// lam_nbr, mu_nbr: face-neighbour elements (shared faces), may be empty
   /// on one rank.
   PmlDynamicTerms(ElasticityDomainOperator<ParMesh> &elast,
                   std::shared_ptr<const PmlBand> band,
                   const Vector &lam_e, const Vector &mu_e, const Vector &rho_e,
                   const Vector &lam_nbr, const Vector &mu_nbr)
      : elast_(elast), band_(std::move(band)), D_(band_->Dim())
   {
      MFEM_VERIFY(elast_.GetMesh().Dimension() == D_, "PmlDynamicTerms: dimension");
      MFEM_VERIFY(D_ == 2 || D_ == 3, "PmlDynamicTerms: 2-D or 3-D, not " << D_);
      eps_ = elast_.GetEpsilon();
      Build_(lam_e, mu_e, rho_e, lam_nbr, mu_nbr);
   }

   /// The filter coefficients for this step (once, before Enter).
   void SetTimeStep(real_t dt)
   {
      MFEM_VERIFY(dt > 0.0, "PmlDynamicTerms::SetTimeStep: dt = " << dt);
      dt_ = dt;
      for (Coef &c : coefs_) { c.Set(band_->Alpha(), dt); }
      dev_coef_dirty_ = true;
   }
   real_t TimeStep() const { return dt_; }

   /// C_d of (3.10): rho sum_a d_a (N_i . N_j) on the band elements, added to
   /// `blocks` (element-local byNODES, D nd x D nd).
   void AddDampingBlocks(std::map<int, DenseMatrix> &blocks) const
   {
      const FiniteElementSpace &fes = elast_.GetFESpace();
      for (const VolElem &ve : vol_)
      {
         const int nd = fes.GetFE(ve.e)->GetDof();
         auto it = blocks.find(ve.e);
         if (it == blocks.end())
         {
            it = blocks.emplace(ve.e, DenseMatrix(D_ * nd)).first;
            it->second = 0.0;
         }
         DenseMatrix &B = it->second;
         for (int q = 0; q < ve.nq; q++)
         {
            const VolPt &p = vpt_[ve.p0 + q];
            const real_t *phi = &vshape_[p.shape];
            real_t dsum = 0.0;
            for (int k = 0; k < p.nact; k++) { dsum += coefs_[p.coef].d[k]; }
            const real_t s = p.wdet * ve.rho * dsum;
            for (int k = 0; k < nd; k++)
               for (int l = 0; l < nd; l++)
               {
                  const real_t v = s * phi[k] * phi[l];
                  for (int c = 0; c < D_; c++) { B(c * nd + k, c * nd + l) += v; }
               }
         }
      }
   }

   /// The largest penalty factor kappa_gamma (3.14) over the band's volume
   /// and face points, all ranks (collective).  The SIPG penalty scales by it,
   /// so the explicit step's stable dt scales by about 1 / sqrt(kappa_max)
   /// (measured 2-D: test_pml_dyn CFL line).
   real_t KappaMax() const
   {
      real_t k = 1.0;
      for (const BandPt &bp : bpt_) { k = std::max(k, bp.pp.kappa); }
      const ParMesh &pm = static_cast<const ParMesh &>(elast_.GetMesh());
      MPI_Allreduce(MPI_IN_PLACE, &k, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, pm.GetComm());
      return k;
   }

   /// k(e): the largest kappa_gamma over the band's volume and face points of
   /// local element e (both local sides of a face), 1 where the band adds
   /// none.  The SIPG penalty of a face scales by kappa at its points, so
   /// the stable step of e alone scales by about 1 / sqrt(k(e)).
   void ElementKappa(Vector &k) const
   {
      k.SetSize(elast_.GetMesh().GetNE());
      k = 1.0;
      for (const VolElem &ve : vol_)
      {
         for (int q = 0; q < ve.nq; q++)
         {
            k(ve.e) = std::max(k(ve.e), bpt_[vpt_[ve.p0 + q].band].pp.kappa);
         }
      }
      for (const Face &f : faces_)
      {
         for (int q = 0; q < f.nq; q++)
         {
            const real_t kq = bpt_[fpt_[f.p0 + q].band].pp.kappa;
            k(f.e1) = std::max(k(f.e1), kq);
            if (f.kind == 0) { k(f.e2) = std::max(k(f.e2), kq); }
         }
      }
   }

   /// t(e) = 1 for the local elements whose rows the band terms reach
   /// (AddForce) or whose values they read: the elements with a band volume
   /// point and the local sides of every band face.  Local time stepping
   /// (NitscheLtsRun) keeps all of them in one cluster.
   void TouchedElements(std::vector<char> &t) const
   {
      t.assign(elast_.GetMesh().GetNE(), 0);
      for (const VolElem &ve : vol_) { t[ve.e] = 1; }
      for (const Face &f : faces_)
      {
         t[f.e1] = 1;
         if (f.kind == 0) { t[f.e2] = 1; }
      }
   }

   int NumBandElements() const { return static_cast<int>(vol_.size()); }
   int NumBandFaces() const { return static_cast<int>(faces_.size()); }
   size_t MemoryReals() const { return state_.size(); }

   /// The 2.5D term in the band (header): beta = beta_over_mu mu_e, U_pl1 the
   /// substrate displacement at physical time 1 (vdofs), t_offset the
   /// physical time of the operator's t = 0 (null = 0), as
   /// ElasticityDomainOperator::SetFiniteWidth.  Before Enter.  Collective.
   void SetFiniteWidth(real_t beta_over_mu, const Vector &U_pl1,
                       std::shared_ptr<const real_t> t_offset)
   {
      MFEM_VERIFY(!entered_, "PmlDynamicTerms::SetFiniteWidth after Enter");
      MFEM_VERIFY(std::isfinite(beta_over_mu) && beta_over_mu > 0.0,
                  "PmlDynamicTerms::SetFiniteWidth: beta / mu = " << beta_over_mu);
      beta_over_mu_ = beta_over_mu;
      beta_t_off_ = std::move(t_offset);
      Gather_(U_pl1);
      BetaStatic_(fb1_);
   }
   bool HasFiniteWidth() const { return beta_over_mu_ > 0.0; }

   /// Entry at the switch (Section 3.1, Proposition 2): memory to zero,
   /// u_s kept, b_band = band terms of (C0 - C) on u_s with the wall jump
   /// u_s - g(t_s); t: the operator's time at entry.  Collective (face
   /// neighbour exchange).
   void Enter(const Vector &u_s, real_t t)
   {
      us_ = u_s;
      std::fill(state_.begin(), state_.end(), 0.0);
      if (dev_built_) { d_state_ = 0.0; }   // the device form's memory
      // g(t_s) at the wall face points
      for (Face &f : faces_)
      {
         if (!f.wall) { continue; }
         for (int q = 0; q < f.nq; q++)
         {
            FacePt &fp = fpt_[f.p0 + q];
            Vector x(fp.x, D_), g(D_);
            elast_.EvalDirichlet(f.attr, x, t, g);
            for (int c = 0; c < D_; c++) { fp.gs[c] = g(c); }
         }
      }
      // b_band: the static band terms on u_s
      b_band_.SetSize(u_s.Size());
      b_band_ = 0.0;
      Gather_(u_s);
      Terms_(b_band_, /*static_mode=*/true, t, /*jump_from_gs=*/true);
      if (beta_over_mu_ > 0.0) { BetaStatic_(fbs_); }   // beta (S0 - 1) u_s
      t_last_ = t;
      fcache_.SetSize(u_s.Size());
      fcache_ = 0.0;
      // the device form's memory is zero here, and so is its B + E
      if (!dev_) { ForceFromState_(fcache_); }
      entered_ = true;
   }

   /// One step: w = x - u_s at the step's new time t; every filter by (3.16).
   /// Collective.
   void Advance(const Vector &x, real_t t)
   {
      MFEM_VERIFY(entered_, "PmlDynamicTerms::Advance before Enter");
      const RankCosts::Scope cost(RankCosts::Pml);
      if (dev_)
      {
         UpdateAllDevice_(x);
         t_last_ = t;
         return;
      }
      w_ = x;
      w_ -= us_;
      Gather_(w_);
      if (D_ == 2) { UpdateAll_<2>(t); }
      else { UpdateAll_<3>(t); }
      t_last_ = t;
   }

   /// F_int += b_band + B(w, Phi) + E(w, chi), from the stored state (the
   /// B + E vector is formed once per step, in Advance, and at Enter).
   void AddForce(Vector &F) const
   {
      MFEM_VERIFY(entered_, "PmlDynamicTerms::AddForce before Enter");
      MFEM_VERIFY(fcache_.Size() == F.Size(), "PmlDynamicTerms::AddForce: size");
      const RankCosts::Scope cost(RankCosts::Pml);
      if (dev_) { b_band_.UseDevice(true); }
      F += b_band_;
      F += fcache_;   // B + E of the stored state (UpdateAll_ / Enter)
      if (beta_over_mu_ > 0.0)
      {
         // 2.5D: beta (S0 - 1) (u_s - u_pl(t)); the chi part is in fcache_
         F += fbs_;
         F.Add(-((beta_t_off_ ? *beta_t_off_ : 0.0) + t_last_), fb1_);
      }
   }

   const Vector &BandLoad() const { return b_band_; }

   /// The device form (UCERF 3-D plan rev5, P1.10): Advance steps every
   /// pipeline in kernels on the mfem device and leaves B + E there; Enter
   /// and the static terms stay on the host (once per phase).  3-D, one
   /// element type, no band wall (Dirichlet) faces and no 2.5D term yet.
   /// The arithmetic per point is the host's; the rows of an element are
   /// summed in the host's order (its volume points, then its faces in
   /// face order).  Collective on the first Advance.
   void SetDevice(bool on)
   {
      MFEM_VERIFY(!on || Device::IsEnabled(), "PmlDynamicTerms::SetDevice: no mfem device");
      MFEM_VERIFY(!on || D_ == 3, "PmlDynamicTerms::SetDevice: 3-D only");
      dev_ = on;
   }
   bool OnDevice() const { return dev_; }

   // The functions holding the kernels are public (nvcc extended lambdas).
   void UpdateAllDevice_(const Vector &x)
   {
      if (!dev_built_) { BuildDevice_(); }
      if (dev_coef_dirty_) { PackCoefs_(); }
      const int n = us_.Size(), nd = dnd_, nv = 3 * nd;
      // w = x - u_s, with the face-neighbour values
      {
         const real_t *xx = x.Read(true), *us = us_.Read(true);
         real_t *g = gf_->Write(true);
         mfem::forall(n, [=] MFEM_HOST_DEVICE (int i) { g[i] = xx[i] - us[i]; });
      }
      {
         const RankCosts::Scope wait(RankCosts::WPml);
         gf_->ExchangeFaceNbrData();
      }
      fcache_.SetSize(n);
      fcache_.UseDevice(true);
      fcache_ = 0.0;
      const real_t al = band_->Alpha();
      const real_t *w = gf_->Read(true);
      const real_t *co = d_coef_.Read(true), *fr = d_frame_.Read(true);
      const int *bi = d_bi_.Read(true), *vd = d_vd_.Read(true);
      real_t *st = d_state_.ReadWrite(true);
      real_t *fc = fcache_.ReadWrite(true);
      const int nve = static_cast<int>(vol_.size());
      if (nve > 0)
      {
         const int *ie = d_vie_.Read(true), *ip = d_vip_.Read(true);
         const real_t *re = d_vre_.Read(true), *rp = d_vrp_.Read(true);
         const real_t *vsh = d_vsh_.Read(true), *vdsh = d_vdsh_.Read(true);
         mfem::forall(nve, [=] MFEM_HOST_DEVICE (int e)
         {
            const int p0 = ie[3 * e], nq = ie[3 * e + 1];
            const int *v = vd + ie[3 * e + 2];
            const real_t lam = re[3 * e], mu = re[3 * e + 1], rho = re[3 * e + 2];
            real_t ue[3 * pmldev::kMaxNd], r[3 * pmldev::kMaxNd];
            for (int i = 0; i < nv; i++) { ue[i] = w[v[i]]; r[i] = 0.0; }
            for (int q = 0; q < nq; q++)
            {
               const int *pi = ip + 5 * (p0 + q);
               const int b = pi[0], na = bi[2 * b];
               if (na == 0) { continue; }
               const real_t *c = co + pmldev::kCoef * bi[2 * b + 1], *nf = fr + 9 * b;
               const real_t *sh = vsh + pi[1], *dsh = vdsh + pi[2];
               const real_t wdet = rp[p0 + q];
               real_t G[3][3], wv[3], ds[3][3], E[3];
               for (int i = 0; i < 3; i++)
                  for (int l = 0; l < 3; l++)
                  {
                     real_t s = 0.0;
                     for (int k = 0; k < nd; k++) { s += ue[i * nd + k] * dsh[k * 3 + l]; }
                     G[i][l] = s;
                  }
               pmldev::StepPipe(na, nf, c, lam, mu, G, st + pi[3], ds);
               for (int i = 0; i < 3; i++)
               {
                  real_t s = 0.0;
                  for (int k = 0; k < nd; k++) { s += ue[i * nd + k] * sh[k]; }
                  wv[i] = s;
               }
               pmldev::StepMass(na, c, wv, st + pi[4]);
               for (int k = 0; k < nd; k++)
                  for (int p = 0; p < 3; p++)
                  {
                     real_t s = 0.0;
                     for (int m = 0; m < 3; m++) { s += dsh[k * 3 + m] * ds[p][m]; }
                     r[p * nd + k] += wdet * s;
                  }
               pmldev::MassE(na, c, al, st + pi[4], E);
               for (int k = 0; k < nd; k++)
                  for (int cc = 0; cc < 3; cc++) { r[cc * nd + k] += wdet * rho * E[cc] * sh[k]; }
            }
            for (int i = 0; i < nv; i++) { fc[v[i]] = r[i]; }
         });
      }
      const int nf = static_cast<int>(faces_.size());
      if (nf > 0)
      {
         const int *ie = d_fie_.Read(true), *ip = d_fip_.Read(true);
         const int *nvd = (d_nvd_.Size() > 0) ? d_nvd_.Read(true) : nullptr;
         const real_t *wn = (gf_->FaceNbrData().Size() > 0) ? gf_->FaceNbrData().Read(true)
                                                            : nullptr;
         const real_t *re = d_fre_.Read(true), *rp = d_frp_.Read(true);
         const real_t *fsh = d_fsh_.Read(true), *fdsh = d_fdsh_.Read(true);
         const real_t eps = eps_;
         real_t *fb = d_fbuf_.Write(true);
         mfem::forall(nf, [=] MFEM_HOST_DEVICE (int f)
         {
            const int kind = ie[5 * f], nq = ie[5 * f + 1], p0 = ie[5 * f + 2];
            const int *v1 = vd + ie[5 * f + 3];
            real_t u[2][3 * pmldev::kMaxNd], r[2][3 * pmldev::kMaxNd];
            for (int i = 0; i < nv; i++) { u[0][i] = w[v1[i]]; r[0][i] = 0.0; r[1][i] = 0.0; }
            if (kind == 0)
            {
               const int *v2 = vd + ie[5 * f + 4];
               for (int i = 0; i < nv; i++) { u[1][i] = w[v2[i]]; }
            }
            else
            {
               const int *v2 = nvd + ie[5 * f + 4];
               for (int i = 0; i < nv; i++) { u[1][i] = wn[v2[i]]; }
            }
            const real_t lm[2][2] = {{re[4 * f], re[4 * f + 1]}, {re[4 * f + 2], re[4 * f + 3]}};
            const int nlocal = (kind == 1) ? 1 : 2;
            for (int q = 0; q < nq; q++)
            {
               const int *pi = ip + 9 * (p0 + q);
               const real_t *pr = rp + 4 * (p0 + q);
               const real_t wq = pr[0], *nor = pr + 1;
               const int b = pi[0], na = bi[2 * b];
               const real_t *c = co + pmldev::kCoef * bi[2 * b + 1], *nfm = fr + 9 * b;
               real_t a[2][3], jmp[3];
               for (int s = 0; s < 2; s++)
               {
                  const real_t *sh = fsh + pi[1 + s];
                  for (int i = 0; i < 3; i++)
                  {
                     real_t t = 0.0;
                     for (int k = 0; k < nd; k++) { t += u[s][i * nd + k] * sh[k]; }
                     a[s][i] = t;
                  }
               }
               for (int cc = 0; cc < 3; cc++) { jmp[cc] = a[0][cc] - a[1][cc]; }
               real_t J[3][3];
               for (int i = 0; i < 3; i++)
                  for (int l = 0; l < 3; l++) { J[i][l] = jmp[i] * nor[l]; }
               real_t dcs[2][3][3], dss[2][3][3];
               for (int s = 0; s < 2; s++)
               {
                  const real_t *dsh = fdsh + pi[3 + s];
                  real_t Gs[3][3];
                  for (int i = 0; i < 3; i++)
                     for (int l = 0; l < 3; l++)
                     {
                        real_t t = 0.0;
                        for (int k = 0; k < nd; k++) { t += u[s][i * nd + k] * dsh[k * 3 + l]; }
                        Gs[i][l] = t;
                     }
                  for (int i = 0; i < 3; i++)
                     for (int j = 0; j < 3; j++) { dcs[s][i][j] = 0.0; dss[s][i][j] = 0.0; }
                  if (na == 0) { continue; }
                  pmldev::StepPipe(na, nfm, c, lm[s][0], lm[s][1], Gs, st + pi[5 + s], dcs[s]);
                  pmldev::StepPipe(na, nfm, c, lm[s][0], lm[s][1], J, st + pi[7 + s], dss[s]);
               }
               real_t avg_n[3];
               for (int i = 0; i < 3; i++)
               {
                  real_t t = 0.0;
                  for (int l = 0; l < 3; l++) { t += 0.5 * (dcs[0][i][l] + dcs[1][i][l]) * nor[l]; }
                  avg_n[i] = t;
               }
               const real_t csym = 0.5 * eps;
               for (int s = 0; s < nlocal; s++)
               {
                  const real_t *sh = fsh + pi[1 + s], *dsh = fdsh + pi[3 + s];
                  const real_t sgn = (s == 0) ? -1.0 : 1.0;
                  for (int k = 0; k < nd; k++)
                     for (int p = 0; p < 3; p++)
                     {
                        real_t t = sgn * sh[k] * avg_n[p];
                        real_t sy = 0.0;
                        for (int m = 0; m < 3; m++) { sy += dsh[k * 3 + m] * dss[s][p][m]; }
                        t += csym * sy;
                        r[s][p * nd + k] += wq * t;
                     }
               }
            }
            real_t *o = fb + static_cast<std::size_t>(2 * nv) * f;
            for (int i = 0; i < nv; i++) { o[i] = r[0][i]; o[nv + i] = r[1][i]; }
         });
         // each element adds its faces' terms, in face order
         const int nt = d_te_.Size();
         const int *te = d_te_.Read(true), *tp = d_tptr_.Read(true), *tf = d_tfs_.Read(true);
         const real_t *fbr = d_fbuf_.Read(true);
         mfem::forall(nt, [=] MFEM_HOST_DEVICE (int t)
         {
            const int *v = vd + te[t];
            for (int j = tp[t]; j < tp[t + 1]; j++)
            {
               const real_t *o = fbr + static_cast<std::size_t>(nv) * tf[j];
               for (int i = 0; i < nv; i++) { fc[v[i]] += o[i]; }
            }
         });
      }
   }

   /// B_0 + K_p acting on a field u (static mode, wall jump u - g(t)): for
   /// tests, K0 u = K_p u + StaticBand(u) to rounding.  Collective.
   void StaticBand(const Vector &u, real_t t, Vector &out)
   {
      out.SetSize(u.Size());
      out = 0.0;
      Gather_(u);
      Terms_(out, true, t, false);
      if (beta_over_mu_ > 0.0)
      {
         Vector fb;
         BetaStatic_(fb);   // 2.5D: + beta (S0 - 1) u
         out += fb;
      }
   }

private:
   // ---------------- per-point filter coefficients
   struct Coef
   {
      int nact = 0;
      real_t d[3] = {0, 0, 0};
      real_t eB[3], gB[3];   // filter B: q = alpha + d
      real_t eA[3], gA[3];   // filter A: q = alpha
      void Set(real_t alpha, real_t dt)
      {
         for (int k = 0; k < nact; k++)
         {
            const real_t qB = alpha + d[k], qA = alpha;
            eB[k] = std::exp(-qB * dt);
            gB[k] = d[k] / qB * (1.0 - eB[k]) * 0.5;
            eA[k] = std::exp(-qA * dt);
            gA[k] = d[k] / qA * (1.0 - eA[k]) * 0.5;
         }
      }
   };
   // the band at a point: active walls, frame, coefficients
   struct BandPt
   {
      int nact = 0;
      real_t n[3][3];      // frame: n[b] for b < D; the first nact are active normals
      int coef = -1;
      PmlPoint pp;         // zero-frequency data (static mode)
   };
   struct VolPt
   {
      int nact = 0;
      int coef = -1;
      int band = -1;       // index into bpt_
      real_t wdet = 0.0;   // weight * det J
      int shape = 0;       // offset into vshape_ (nd values)
      int dshape = 0;      // offset into vdshape_ (nd * D values, physical)
      size_t st_g = 0;     // pipeline state (stress)
      size_t st_m = 0;     // mass-side state
   };
   struct VolElem
   {
      int e = -1, nq = 0, p0 = 0;
      int nd = 0, vd = 0;  // dofs per component; offset into vdofs_ (D nd entries)
      real_t lam = 0, mu = 0, rho = 0;
   };
   struct FacePt
   {
      int band = -1;
      real_t w = 0.0;      // quadrature weight (the face measure is in nor)
      real_t nor[3] = {0, 0, 0};
      real_t x[3] = {0, 0, 0};
      real_t gs[3] = {0, 0, 0};
      int shape[2] = {0, 0}, dshape[2] = {0, 0};
      size_t st_c[2] = {0, 0}, st_s[2] = {0, 0};
   };
   struct Face
   {
      int kind = 0;        // 0 interior, 1 shared, 2 wall
      int e1 = -1, e2 = -1;   // e2: local element (interior) or face-neighbour index (shared)
      int nd1 = 0, nd2 = 0;
      bool wall = false;
      int attr = 0;
      int nq = 0, p0 = 0;
      int vd1 = 0, vd2 = 0;   // offsets into vdofs_ (e1, local e2) or nbr_vdofs_ (shared e2)
      real_t lam[2] = {0, 0}, mu[2] = {0, 0};
   };

   ElasticityDomainOperator<ParMesh> &elast_;
   std::shared_ptr<const PmlBand> band_;
   int D_;
   real_t eps_ = -1.0;
   real_t dt_ = 0.0;
   bool entered_ = false;
   real_t t_last_ = 0.0;

   std::vector<Coef> coefs_;
   std::vector<BandPt> bpt_;
   std::vector<VolPt> vpt_;
   std::vector<VolElem> vol_;
   std::vector<FacePt> fpt_;
   std::vector<Face> faces_;
   std::vector<real_t> vshape_, vdshape_, fshape_, fdshape_;
   std::vector<real_t> state_;
   std::vector<int> vdofs_, nbr_vdofs_;   // element vdofs, cached by Build_

   Vector us_, w_, b_band_;
   // 2.5D finite width (SetFiniteWidth; 0 = off): int beta (S0 - 1) f . v
   // on the band points for f = u_s (fbs_) and f = U_pl1 (fb1_)
   real_t beta_over_mu_ = 0.0;
   std::shared_ptr<const real_t> beta_t_off_;
   Vector fbs_, fb1_;
   Vector fcache_;   // B(w, Phi) + E(w, chi) of the stored state
   // the device form (SetDevice, BuildDevice_): the tables above as flat arrays
   bool dev_ = false, dev_built_ = false, dev_coef_dirty_ = true;
   int dnd_ = 0;
   Vector d_coef_, d_frame_, d_vre_, d_vrp_, d_vsh_, d_vdsh_, d_fre_, d_frp_, d_fsh_, d_fdsh_;
   Vector d_state_, d_fbuf_;
   Array<int> d_bi_, d_vie_, d_vip_, d_fie_, d_fip_, d_vd_, d_nvd_, d_te_, d_tptr_, d_tfs_;
   mutable Vector ue_, u1_, u2_, r_, r1_, r2_;   // scratch (per element / face)
   std::unique_ptr<ParGridFunction> gf_;   // values with face-neighbour data

   // ---------------- the pipeline at one point
   // state layout for an input G (stress pipeline): Gprev[D*D], phi[nact][D],
   // then for b = 0..D-1, for each active a != b: m[D], z[D].
   int PipeSize_(int nact) const
   {
      int f3 = 0;
      for (int b = 0; b < D_; b++) { f3 += (b < nact) ? nact - 1 : nact; }
      return D_ * D_ + nact * D_ + f3 * 2 * D_;
   }
   // mass side: wprev[D], then per active a: chi[D], z[D]
   int MassSize_(int nact) const { return D_ + nact * 2 * D_; }

   size_t Alloc_(int n)
   {
      const size_t off = state_.size();
      state_.resize(off + n, 0.0);
      return off;
   }

   // The per-point kernels take the dimension D as a template argument (the
   // loops unroll; the same operations in the same order as with a run-time
   // D, so the same bits); the ...R_ forms dispatch on D_ for the cold paths.
   template <int D>
   void Hooke_(real_t lam, real_t mu, const real_t H[3][3], real_t s[3][3]) const
   {
      real_t tr = 0.0;
      for (int k = 0; k < D; k++) { tr += H[k][k]; }
      for (int i = 0; i < D; i++)
         for (int j = 0; j < D; j++)
         { s[i][j] = mu * (H[i][j] + H[j][i]) + (i == j ? lam * tr : 0.0); }
   }

   void HookeR_(real_t lam, real_t mu, const real_t H[3][3], real_t s[3][3]) const
   {
      if (D_ == 2) { Hooke_<2>(lam, mu, H, s); }
      else { Hooke_<3>(lam, mu, H, s); }
   }

   /// Steps the stress pipeline at one point with the new input G (3.16).
   /// With ds, also Delta sigma of the new state (PipeDelta_'s value, the
   /// same operations: the cascade outputs y_b are the z it forms here).
   template <int D>
   void StepPipe_(const BandPt &bp, real_t lam, real_t mu, const real_t G[3][3],
                  real_t *st, real_t (*ds)[3] = nullptr) const
   {
      const int na = bp.nact;
      const Coef &c = coefs_[bp.coef];
      real_t *Gp = st;
      real_t *phi = st + D * D;
      // stage 1: z_a = G n_a ; z_a^n from Gprev
      for (int k = 0; k < na; k++)
      {
         const real_t *n = bp.n[k];
         for (int i = 0; i < D; i++)
         {
            real_t zn = 0.0, zo = 0.0;
            for (int l = 0; l < D; l++) { zn += G[i][l] * n[l]; zo += Gp[i * D + l] * n[l]; }
            phi[k * D + i] = c.eB[k] * phi[k * D + i] + c.gB[k] * (zo + zn);
         }
      }
      for (int i = 0; i < D * D; i++) { Gp[i] = G[i / D][i % D]; }
      // stage 2 and 3 with the updated phi: the cascades' new inputs
      real_t st2[3][3];
      Stage12_<D>(bp, lam, mu, st, st2);
      real_t *f = phi + na * D;
      real_t sig[3][3] = {{0}};
      for (int b = 0; b < D; b++)
      {
         real_t z[3];
         for (int i = 0; i < D; i++)
         {
            z[i] = 0.0;
            for (int l = 0; l < D; l++) { z[i] += st2[i][l] * bp.n[b][l]; }
         }
         for (int k = 0; k < na; k++)
         {
            if (k == b) { continue; }
            real_t *m = f, *zo = f + D;
            for (int i = 0; i < D; i++)
            {
               m[i] = c.eA[k] * m[i] + c.gA[k] * (zo[i] + z[i]);
               zo[i] = z[i];
               z[i] += m[i];          // this filter's output feeds the next
            }
            f += 2 * D;
         }
         if (ds)
         {
            for (int i = 0; i < D; i++)
               for (int j = 0; j < D; j++) { sig[i][j] += z[i] * bp.n[b][j]; }
         }
      }
      if (ds)
      {
         real_t Gn[3][3], s0[3][3];
         for (int i = 0; i < D; i++)
            for (int l = 0; l < D; l++) { Gn[i][l] = Gp[i * D + l]; }
         Hooke_<D>(lam, mu, Gn, s0);
         for (int i = 0; i < D; i++)
            for (int j = 0; j < D; j++) { ds[i][j] = sig[i][j] - s0[i][j]; }
      }
   }

   /// G~ from the stored state, then sigma~ (stages 1 and 2).
   template <int D>
   void Stage12_(const BandPt &bp, real_t lam, real_t mu, const real_t *st,
                 real_t s2[3][3]) const
   {
      const real_t *Gp = st, *phi = st + D * D;
      real_t Gt[3][3];
      for (int i = 0; i < D; i++)
         for (int l = 0; l < D; l++)
         {
            real_t v = Gp[i * D + l];
            for (int k = 0; k < bp.nact; k++) { v -= phi[k * D + i] * bp.n[k][l]; }
            Gt[i][l] = v;
         }
      Hooke_<D>(lam, mu, Gt, s2);
   }
   void Stage12R_(const BandPt &bp, real_t lam, real_t mu, const real_t *st,
                  real_t s2[3][3]) const
   {
      if (D_ == 2) { Stage12_<2>(bp, lam, mu, st, s2); }
      else { Stage12_<3>(bp, lam, mu, st, s2); }
   }

   /// Delta sigma = sigma_pml - C : Gprev from the stored state.
   void PipeDelta_(const BandPt &bp, real_t lam, real_t mu, const real_t *st,
                   real_t ds[3][3]) const
   {
      const int D = D_, na = bp.nact;
      real_t s2[3][3];
      Stage12R_(bp, lam, mu, st, s2);
      real_t sig[3][3] = {{0}};
      const real_t *f = st + D * D + na * D;
      for (int b = 0; b < D; b++)
      {
         real_t y[3];
         for (int i = 0; i < D; i++)
         {
            y[i] = 0.0;
            for (int l = 0; l < D; l++) { y[i] += s2[i][l] * bp.n[b][l]; }
         }
         for (int k = 0; k < na; k++)
         {
            if (k == b) { continue; }
            for (int i = 0; i < D; i++) { y[i] += f[i]; }   // y = z + m (cascade)
            f += 2 * D;
         }
         for (int i = 0; i < D; i++)
            for (int j = 0; j < D; j++) { sig[i][j] += y[i] * bp.n[b][j]; }
      }
      real_t G[3][3], s0[3][3];
      for (int i = 0; i < D; i++)
         for (int l = 0; l < D; l++) { G[i][l] = st[i * D + l]; }
      HookeR_(lam, mu, G, s0);
      for (int i = 0; i < D; i++)
         for (int j = 0; j < D; j++) { ds[i][j] = sig[i][j] - s0[i][j]; }
   }

   /// Zero-frequency Delta sigma: C0 : G - C : G.
   void StaticDelta_(const BandPt &bp, real_t lam, real_t mu, const real_t G[3][3],
                     real_t ds[3][3]) const
   {
      real_t s0[3][3], s[3][3];
      band_->ApplyC0(bp.pp, lam, mu, G, s0);
      HookeR_(lam, mu, G, s);
      for (int i = 0; i < D_; i++)
         for (int j = 0; j < D_; j++) { ds[i][j] = s0[i][j] - s[i][j]; }
   }

   /// Mass side: steps chi with the new w at the point.
   template <int D>
   void StepMass_(const BandPt &bp, const real_t *w, real_t *st) const
   {
      const Coef &c = coefs_[bp.coef];
      real_t *wp = st;
      real_t y[3];
      for (int i = 0; i < D; i++) { y[i] = w[i]; }
      real_t *f = st + D;
      for (int k = 0; k < bp.nact; k++)
      {
         real_t *chi = f, *zo = f + D;
         for (int i = 0; i < D; i++)
         {
            chi[i] = c.eA[k] * chi[i] + c.gA[k] * (zo[i] + y[i]);
            zo[i] = y[i];
            y[i] += chi[i];
         }
         f += 2 * D;
      }
      for (int i = 0; i < D; i++) { wp[i] = w[i]; }
   }

   /// E_J (3.9) from the stored mass state:
   ///   Yd_a = Yd_(a-1) + d_a y_(a-1) - alpha chi_a   (Yd = y' - w')
   ///   E_a  = E_(a-1) + d_a Yd_(a-1) - alpha (d_a y_(a-1) - alpha chi_a)
   template <int D>
   void MassE_(const BandPt &bp, const real_t *st, real_t E[3]) const
   {
      const Coef &c = coefs_[bp.coef];
      const real_t al = band_->Alpha();
      real_t y[3], Yd[3];
      for (int i = 0; i < D; i++) { y[i] = st[i]; Yd[i] = 0.0; E[i] = 0.0; }
      const real_t *f = st + D;
      for (int k = 0; k < bp.nact; k++)
      {
         const real_t *chi = f;
         for (int i = 0; i < D; i++)
         {
            const real_t dy = c.d[k] * y[i] - al * chi[i];   // chi_a'
            E[i] += c.d[k] * Yd[i] - al * dy;
            Yd[i] += dy;
            y[i] += chi[i];
         }
         f += 2 * D;
      }
   }

   /// y_J - w = sum_a chi_a from the stored mass state (2.5D term).
   void MassChi_(const BandPt &bp, const real_t *st, real_t X[3]) const
   {
      const int D = D_;
      for (int i = 0; i < D; i++) { X[i] = 0.0; }
      const real_t *f = st + D;
      for (int k = 0; k < bp.nact; k++)
      {
         for (int i = 0; i < D; i++) { X[i] += f[i]; }
         f += 2 * D;
      }
   }

   /// out = int beta (S0 - 1) u . v over the band's volume points, u the
   /// gathered field (2.5D term, static part).
   void BetaStatic_(Vector &out)
   {
      const int D = D_;
      Vector &ue = ue_, &r = r_;
      real_t uv[3];
      out.SetSize(elast_.GetFESpace().GetVSize());
      out = 0.0;
      for (const VolElem &ve : vol_)
      {
         const int nd = ve.nd;
         OwnValues_(ve.vd, D * nd, ue);
         r.SetSize(D * nd);
         r = 0.0;
         for (int q = 0; q < ve.nq; q++)
         {
            const VolPt &p = vpt_[ve.p0 + q];
            if (p.nact == 0) { continue; }
            const real_t *sh = &vshape_[p.shape];
            ValAt_(D, nd, ue, sh, uv);
            const real_t s = p.wdet * beta_over_mu_ * ve.mu * (bpt_[p.band].pp.S0 - 1.0);
            for (int k = 0; k < nd; k++)
               for (int c = 0; c < D; c++) { r(c * nd + k) += s * uv[c] * sh[k]; }
         }
         AddElem_(ve.vd, D * nd, r, out);
      }
   }

   // ---------------- setup
   int NewBandPt_(const real_t *x)
   {
      BandPt bp;
      bp.pp = band_->At(x);
      bp.nact = bp.pp.nact;
      const auto &W = band_->Walls();
      for (int k = 0; k < bp.nact; k++)
         for (int i = 0; i < 3; i++) { bp.n[k][i] = (i < D_) ? W[bp.pp.wall[k]].n[i] : 0.0; }
      // complete the frame (Gram-Schmidt against the coordinate axes)
      int nb = bp.nact;
      for (int ax = 0; ax < D_ && nb < D_; ax++)
      {
         real_t v[3] = {0, 0, 0};
         v[ax] = 1.0;
         for (int b = 0; b < nb; b++)
         {
            real_t dot = 0.0;
            for (int i = 0; i < D_; i++) { dot += v[i] * bp.n[b][i]; }
            for (int i = 0; i < D_; i++) { v[i] -= dot * bp.n[b][i]; }
         }
         real_t nn = 0.0;
         for (int i = 0; i < D_; i++) { nn += v[i] * v[i]; }
         if (nn < 1e-6) { continue; }
         for (int i = 0; i < 3; i++) { bp.n[nb][i] = (i < D_) ? v[i] / std::sqrt(nn) : 0.0; }
         nb++;
      }
      MFEM_VERIFY(nb == D_, "PmlDynamicTerms: frame completion failed");
      Coef c;
      c.nact = bp.nact;
      for (int k = 0; k < bp.nact; k++) { c.d[k] = bp.pp.d[k]; }
      coefs_.push_back(c);
      bp.coef = static_cast<int>(coefs_.size()) - 1;
      bpt_.push_back(bp);
      return static_cast<int>(bpt_.size()) - 1;
   }

   void Build_(const Vector &lam_e, const Vector &mu_e, const Vector &rho_e,
               const Vector &lam_nbr, const Vector &mu_nbr)
   {
      ParMesh &mesh = const_cast<ParMesh &>(elast_.GetMesh());
      auto &fes = const_cast<ParFiniteElementSpace &>(elast_.GetFESpace());
      const PmlBand &band = *band_;
      const int N = band.Params().power;
      Vector x(D_);
      // volume
      for (int e = 0; e < mesh.GetNE(); e++)
      {
         if (!PmlBandElement(mesh, e, band)) { continue; }
         const FiniteElement &fe = *fes.GetFE(e);
         ElementTransformation &T = *mesh.GetElementTransformation(e);
         const IntegrationRule &ir = IntRules.Get(fe.GetGeomType(), 2 * fe.GetOrder() + N);
         VolElem ve;
         ve.e = e;
         ve.lam = lam_e(e);
         ve.mu = mu_e(e);
         ve.rho = rho_e(e);
         ve.p0 = static_cast<int>(vpt_.size());
         const int nd = fe.GetDof();
         Vector sh(nd);
         DenseMatrix ds(nd, D_);
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            T.SetIntPoint(&ip);
            T.Transform(ip, x);
            VolPt p;
            p.band = NewBandPt_(x.GetData());
            p.nact = bpt_[p.band].nact;
            p.coef = bpt_[p.band].coef;
            p.wdet = ip.weight * T.Weight();
            fe.CalcShape(ip, sh);
            fe.CalcPhysDShape(T, ds);
            p.shape = static_cast<int>(vshape_.size());
            for (int k = 0; k < nd; k++) { vshape_.push_back(sh(k)); }
            p.dshape = static_cast<int>(vdshape_.size());
            for (int k = 0; k < nd; k++)
               for (int j = 0; j < D_; j++) { vdshape_.push_back(ds(k, j)); }
            p.st_g = Alloc_(PipeSize_(p.nact));
            p.st_m = Alloc_(MassSize_(p.nact));
            vpt_.push_back(p);
         }
         ve.nq = ir.GetNPoints();
         vol_.push_back(ve);
      }
      // faces
      mesh.ExchangeFaceNbrData();
      fes.ExchangeFaceNbrData();
      auto add_face = [&](FaceElementTransformations &Tr, int kind, int attr,
                          real_t l1, real_t m1, real_t l2, real_t m2)
      {
         const FiniteElement &fe1 = *fes.GetFE(Tr.Elem1No);
         const FiniteElement *fe2 = nullptr;
         if (kind == 0) { fe2 = fes.GetFE(Tr.Elem2No); }
         if (kind == 1) { fe2 = fes.GetFaceNbrFE(Tr.Elem2No - mesh.GetNE()); }
         const int order = 2 * std::max(fe1.GetOrder(), fe2 ? fe2->GetOrder() : 0) + 1;
         const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), order);
         // any point inside a band?
         bool any = false;
         for (int q = 0; q < ir.GetNPoints() && !any; q++)
         {
            Tr.Transform(ir.IntPoint(q), x);
            any = band.InBand(x.GetData());
         }
         if (!any) { return; }
         Face f;
         f.kind = kind;
         f.wall = (kind == 2);
         f.attr = attr;
         f.e1 = Tr.Elem1No;
         f.e2 = (kind == 0) ? Tr.Elem2No : (kind == 1 ? Tr.Elem2No - mesh.GetNE() : -1);
         f.nd1 = fe1.GetDof();
         f.nd2 = fe2 ? fe2->GetDof() : 0;
         f.lam[0] = l1; f.mu[0] = m1; f.lam[1] = l2; f.mu[1] = m2;
         f.p0 = static_cast<int>(fpt_.size());
         f.nq = ir.GetNPoints();
         const int nsides = fe2 ? 2 : 1;
         Vector nor(D_), sh;
         DenseMatrix ds;
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            Tr.SetAllIntPoints(&ip);
            Tr.Transform(ip, x);
            FacePt fp;
            fp.band = NewBandPt_(x.GetData());
            fp.w = ip.weight;
            CalcOrtho(Tr.Jacobian(), nor);
            for (int i = 0; i < D_; i++) { fp.nor[i] = nor(i); fp.x[i] = x(i); }
            for (int s = 0; s < nsides; s++)
            {
               const FiniteElement &fe = (s == 0) ? fe1 : *fe2;
               ElementTransformation &ET = (s == 0) ? *Tr.Elem1 : *Tr.Elem2;
               const IntegrationPoint &eip = (s == 0) ? Tr.GetElement1IntPoint()
                                                      : Tr.GetElement2IntPoint();
               const int nd = fe.GetDof();
               sh.SetSize(nd);
               ds.SetSize(nd, D_);
               fe.CalcShape(eip, sh);
               fe.CalcPhysDShape(ET, ds);
               fp.shape[s] = static_cast<int>(fshape_.size());
               for (int k = 0; k < nd; k++) { fshape_.push_back(sh(k)); }
               fp.dshape[s] = static_cast<int>(fdshape_.size());
               for (int k = 0; k < nd; k++)
                  for (int j = 0; j < D_; j++) { fdshape_.push_back(ds(k, j)); }
               const int na = bpt_[fp.band].nact;
               fp.st_c[s] = Alloc_(PipeSize_(na));
               fp.st_s[s] = Alloc_(PipeSize_(na));
            }
            fpt_.push_back(fp);
         }
         faces_.push_back(f);
      };
      for (int fi = 0; fi < mesh.GetNumFaces(); fi++)
      {
         FaceElementTransformations *Tr = mesh.GetInteriorFaceTransformations(fi);
         if (!Tr) { continue; }
         add_face(*Tr, 0, 0, lam_e(Tr->Elem1No), mu_e(Tr->Elem1No),
                  lam_e(Tr->Elem2No), mu_e(Tr->Elem2No));
      }
      for (int sf = 0; sf < mesh.GetNSharedFaces(); sf++)
      {
         FaceElementTransformations *Tr = mesh.GetSharedFaceTransformations(sf, true);
         if (!Tr) { continue; }
         const int nb = Tr->Elem2No - mesh.GetNE();
         add_face(*Tr, 1, 0, lam_e(Tr->Elem1No), mu_e(Tr->Elem1No), lam_nbr(nb), mu_nbr(nb));
      }
      const Array<int> &dm = elast_.DirichletBdrMarker();
      for (int be = 0; be < mesh.GetNBE(); be++)
      {
         const int attr = mesh.GetBdrAttribute(be);
         if (dm.Size() == 0 || dm[attr - 1] == 0) { continue; }
         int face = -1, info = 0;
         mesh.GetBdrElementFace(be, &face, &info);
         if (mesh.FaceIsInterior(face)) { continue; }
         FaceElementTransformations *Tr = mesh.GetBdrFaceTransformations(be);
         if (!Tr) { continue; }
         add_face(*Tr, 2, attr, lam_e(Tr->Elem1No), mu_e(Tr->Elem1No), 0.0, 0.0);
      }
      gf_ = std::make_unique<ParGridFunction>(&fes);
      // element vdofs, once (the per-step loops read them from here)
      Array<int> vd;
      auto own = [&](int e)
      {
         fes.GetElementVDofs(e, vd);
         const int off = static_cast<int>(vdofs_.size());
         for (int i = 0; i < vd.Size(); i++) { vdofs_.push_back(vd[i]); }
         return off;
      };
      for (VolElem &ve : vol_)
      {
         ve.nd = fes.GetFE(ve.e)->GetDof();
         ve.vd = own(ve.e);
      }
      for (Face &f : faces_)
      {
         f.vd1 = own(f.e1);
         if (f.kind == 0) { f.vd2 = own(f.e2); }
         if (f.kind == 1)
         {
            fes.GetFaceNbrElementVDofs(f.e2, vd);
            f.vd2 = static_cast<int>(nbr_vdofs_.size());
            for (int i = 0; i < vd.Size(); i++) { nbr_vdofs_.push_back(vd[i]); }
         }
      }
   }

   // ---------------- values at the points
   void Gather_(const Vector &u)
   {
      gf_->SetFromTrueDofs(u);   // DG: true dofs = local vdofs
      {
         const RankCosts::Scope wait(RankCosts::WPml);
         gf_->ExchangeFaceNbrData();
      }
      // the host loops read both through raw pointers (no-ops without a device)
      gf_->HostRead();
      gf_->FaceNbrData().HostRead();
   }

   // ---------------- the device form (SetDevice)
   void PackCoefs_()
   {
      d_coef_.SetSize(pmldev::kCoef * static_cast<int>(coefs_.size()));
      real_t *h = d_coef_.HostWrite();
      for (std::size_t i = 0; i < coefs_.size(); i++)
      {
         const Coef &c = coefs_[i];
         real_t *o = h + pmldev::kCoef * i;
         for (int k = 0; k < 3; k++)
         {
            const bool a = k < c.nact;
            o[k] = a ? c.d[k] : 0.0;
            o[3 + k] = a ? c.eB[k] : 0.0;
            o[6 + k] = a ? c.gB[k] : 0.0;
            o[9 + k] = a ? c.eA[k] : 0.0;
            o[12 + k] = a ? c.gA[k] : 0.0;
         }
         o[15] = 0.0;
      }
      dev_coef_dirty_ = false;
   }

   void BuildDevice_()
   {
      MFEM_VERIFY(beta_over_mu_ == 0.0, "PmlDynamicTerms: the device form has no 2.5D term");
      MFEM_VERIFY(state_.size() < static_cast<std::size_t>(std::numeric_limits<int>::max()),
                  "PmlDynamicTerms: the device form holds the memory of one rank in one "
                  "vector, at most 2^31 reals; this rank has " << state_.size()
                  << " (use more ranks)");
      const FiniteElementSpace &fes = elast_.GetFESpace();
      dnd_ = (fes.GetNE() > 0) ? fes.GetFE(0)->GetDof() : 0;
      const int nd = dnd_, nv = 3 * nd;
      MFEM_VERIFY(nd <= pmldev::kMaxNd, "PmlDynamicTerms: the device form covers elements of "
                  "at most " << pmldev::kMaxNd << " nodes, not " << nd);
      auto ints = [](const std::vector<int> &h, Array<int> &d)
      {
         d.SetSize(static_cast<int>(h.size()));
         int *p = d.HostWrite();
         for (std::size_t i = 0; i < h.size(); i++) { p[i] = h[i]; }
      };
      auto reals = [](const std::vector<real_t> &h, Vector &d)
      {
         MFEM_VERIFY(h.size() < static_cast<std::size_t>(std::numeric_limits<int>::max()),
                     "PmlDynamicTerms: a device table of " << h.size() << " reals");
         d.SetSize(static_cast<int>(h.size()));
         real_t *p = d.HostWrite();
         for (std::size_t i = 0; i < h.size(); i++) { p[i] = h[i]; }
      };
      for (int v : vdofs_) { MFEM_VERIFY(v >= 0, "PmlDynamicTerms: negative vdof (not DG)"); }
      for (int v : nbr_vdofs_) { MFEM_VERIFY(v >= 0, "PmlDynamicTerms: negative neighbour vdof"); }
      ints(vdofs_, d_vd_);
      ints(nbr_vdofs_, d_nvd_);
      reals(vshape_, d_vsh_);
      reals(vdshape_, d_vdsh_);
      reals(fshape_, d_fsh_);
      reals(fdshape_, d_fdsh_);
      std::vector<int> hi;
      std::vector<real_t> hr;
      for (const BandPt &b : bpt_)
      {
         hi.push_back(b.nact);
         hi.push_back(b.coef);
         for (int k = 0; k < 3; k++)
            for (int i = 0; i < 3; i++) { hr.push_back(b.n[k][i]); }
      }
      ints(hi, d_bi_);
      reals(hr, d_frame_);
      hi.clear();
      hr.clear();
      std::vector<int> pi;
      std::vector<real_t> pr;
      for (const VolElem &ve : vol_)
      {
         MFEM_VERIFY(ve.nd == nd, "PmlDynamicTerms: the device form needs one element type");
         hi.push_back(ve.p0);
         hi.push_back(ve.nq);
         hi.push_back(ve.vd);
         hr.push_back(ve.lam);
         hr.push_back(ve.mu);
         hr.push_back(ve.rho);
         for (int q = 0; q < ve.nq; q++)
         {
            const VolPt &p = vpt_[ve.p0 + q];
            pi.push_back(p.band);
            pi.push_back(p.shape);
            pi.push_back(p.dshape);
            pi.push_back(static_cast<int>(p.st_g));
            pi.push_back(static_cast<int>(p.st_m));
            pr.push_back(p.wdet);
         }
      }
      MFEM_VERIFY(pr.size() == vpt_.size(), "PmlDynamicTerms: volume points out of order");
      ints(hi, d_vie_);
      reals(hr, d_vre_);
      ints(pi, d_vip_);
      reals(pr, d_vrp_);
      hi.clear();
      hr.clear();
      pi.clear();
      pr.clear();
      // element -> (face, side) in face order; an element's vdofs offset
      std::map<int, std::pair<int, std::vector<int>>> touch;
      for (std::size_t fi = 0; fi < faces_.size(); fi++)
      {
         const Face &f = faces_[fi];
         MFEM_VERIFY(f.kind != 2, "PmlDynamicTerms: the device form has no band wall "
                     "(Dirichlet) faces yet (UCERF 3-D plan rev5, P1.11); attribute " << f.attr);
         MFEM_VERIFY(f.nd1 == nd && f.nd2 == nd, "PmlDynamicTerms: one element type");
         MFEM_VERIFY(f.p0 == static_cast<int>(pr.size() / 4), "PmlDynamicTerms: face points "
                     "out of order");
         hi.push_back(f.kind);
         hi.push_back(f.nq);
         hi.push_back(f.p0);
         hi.push_back(f.vd1);
         hi.push_back(f.vd2);
         hr.push_back(f.lam[0]);
         hr.push_back(f.mu[0]);
         hr.push_back(f.lam[1]);
         hr.push_back(f.mu[1]);
         for (int q = 0; q < f.nq; q++)
         {
            const FacePt &p = fpt_[f.p0 + q];
            pi.push_back(p.band);
            pi.push_back(p.shape[0]);
            pi.push_back(p.shape[1]);
            pi.push_back(p.dshape[0]);
            pi.push_back(p.dshape[1]);
            pi.push_back(static_cast<int>(p.st_c[0]));
            pi.push_back(static_cast<int>(p.st_c[1]));
            pi.push_back(static_cast<int>(p.st_s[0]));
            pi.push_back(static_cast<int>(p.st_s[1]));
            pr.push_back(p.w);
            for (int i = 0; i < 3; i++) { pr.push_back(p.nor[i]); }
         }
         auto &t1 = touch[f.e1];
         t1.first = f.vd1;
         t1.second.push_back(2 * static_cast<int>(fi));
         if (f.kind == 0)
         {
            auto &t2 = touch[f.e2];
            t2.first = f.vd2;
            t2.second.push_back(2 * static_cast<int>(fi) + 1);
         }
      }
      ints(hi, d_fie_);
      reals(hr, d_fre_);
      ints(pi, d_fip_);
      reals(pr, d_frp_);
      std::vector<int> te, tptr(1, 0), tfs;
      for (const auto &kv : touch)
      {
         te.push_back(kv.second.first);
         for (int fs : kv.second.second) { tfs.push_back(fs); }
         tptr.push_back(static_cast<int>(tfs.size()));
      }
      ints(te, d_te_);
      ints(tptr, d_tptr_);
      ints(tfs, d_tfs_);
      MFEM_VERIFY(static_cast<std::size_t>(2 * nv) * faces_.size()
                  < static_cast<std::size_t>(std::numeric_limits<int>::max()),
                  "PmlDynamicTerms: the device form's face buffer is too large for one rank");
      d_fbuf_.SetSize(2 * nv * static_cast<int>(faces_.size()));
      d_fbuf_.UseDevice(true);
      d_state_.SetSize(static_cast<int>(state_.size()));
      d_state_.UseDevice(true);
      {
         // the memory so far (zero right after Enter)
         real_t *h = d_state_.HostWrite();
         for (std::size_t i = 0; i < state_.size(); i++) { h[i] = state_[i]; }
      }
      us_.UseDevice(true);
      gf_->UseDevice(true);
      dev_built_ = true;
   }

   // element values (byNODES) from the cached vdofs (as GetSubVector:
   // a negative index -1-j reads -data[j]), and the element-vector add
   static void Take_(const std::vector<int> &map, int off, int n, const real_t *d, Vector &ue)
   {
      ue.SetSize(n);
      for (int i = 0; i < n; i++)
      {
         const int j = map[off + i];
         ue(i) = (j >= 0) ? d[j] : -d[-1 - j];
      }
   }
   void OwnValues_(int off, int n, Vector &ue) const { Take_(vdofs_, off, n, gf_->GetData(), ue); }
   void NbrValues_(int off, int n, Vector &ue) const
   {
      Take_(nbr_vdofs_, off, n, const_cast<ParGridFunction &>(*gf_).FaceNbrData().GetData(), ue);
   }
   void AddElem_(int off, int n, const Vector &r, Vector &out) const
   {
      Array<int> vd(const_cast<int *>(vdofs_.data()) + off, n);
      out.AddElementVector(vd, r);
   }

   static void GradAt_(int D, int nd, const Vector &ue, const real_t *ds, real_t G[3][3])
   {
      for (int i = 0; i < D; i++)
         for (int l = 0; l < D; l++)
         {
            real_t v = 0.0;
            for (int k = 0; k < nd; k++) { v += ue(i * nd + k) * ds[k * D + l]; }
            G[i][l] = v;
         }
   }
   static void ValAt_(int D, int nd, const Vector &ue, const real_t *sh, real_t u[3])
   {
      for (int i = 0; i < D; i++)
      {
         real_t v = 0.0;
         for (int k = 0; k < nd; k++) { v += ue(i * nd + k) * sh[k]; }
         u[i] = v;
      }
   }

   /// static_mode: Delta sigma at zero frequency on the gathered field into
   /// `out` (the wall jump: u - g(t_s) when jump_from_gs, else u - g(t)).
   void Terms_(Vector &out, bool static_mode, real_t t, bool jump_from_gs)
   {
      MFEM_VERIFY(static_mode, "PmlDynamicTerms::Terms_: static mode only");
      const int D = D_;
      Vector &ue = ue_, &r = r_;
      real_t G[3][3], ds[3][3];
      for (const VolElem &ve : vol_)
      {
         const int nd = ve.nd;
         OwnValues_(ve.vd, D * nd, ue);
         r.SetSize(D * nd);
         r = 0.0;
         for (int q = 0; q < ve.nq; q++)
         {
            const VolPt &p = vpt_[ve.p0 + q];
            if (p.nact == 0) { continue; }
            GradAt_(D, nd, ue, &vdshape_[p.dshape], G);
            StaticDelta_(bpt_[p.band], ve.lam, ve.mu, G, ds);
            AddVol_(D, nd, &vdshape_[p.dshape], p.wdet, ds, r);
         }
         AddElem_(ve.vd, D * nd, r, out);
      }
      for (const Face &f : faces_)
      {
         auto delta = [&](const FacePt &fp, int s, const real_t Gs[3][3], const real_t J[3][3],
                          real_t dc[3][3], real_t dsym[3][3])
         {
            const BandPt &bp = bpt_[fp.band];
            if (bp.nact == 0)
            {
               for (int i = 0; i < 3; i++)
                  for (int j = 0; j < 3; j++) { dc[i][j] = dsym[i][j] = 0.0; }
               return;
            }
            StaticDelta_(bp, f.lam[s], f.mu[s], Gs, dc);
            StaticDelta_(bp, f.lam[s], f.mu[s], J, dsym);
         };
         if (D == 2) { FaceInputs_<2>(f, t, jump_from_gs, delta, out); }
         else { FaceInputs_<3>(f, t, jump_from_gs, delta, out); }
      }
   }

   static void AddVol_(int D, int nd, const real_t *dsh, real_t wdet, const real_t ds[3][3],
                       Vector &r)
   {
      for (int k = 0; k < nd; k++)
         for (int p = 0; p < D; p++)
         {
            real_t v = 0.0;
            for (int m = 0; m < D; m++) { v += dsh[k * D + m] * ds[p][m]; }
            r(p * nd + k) += wdet * v;
         }
   }

   /// For every face point: the inputs of each side (grad w_s, J = [w] x nor),
   /// `delta` gives (Delta sigma of the consistency input, of the symmetry
   /// input) for side s, and the test terms of the local sides go to `out`.
   template <int D, typename DeltaF>
   void FaceInputs_(const Face &f, real_t t, bool jump_from_gs, DeltaF delta, Vector &out)
   {
      Vector &u1 = u1_, &u2 = u2_, &r1 = r1_, &r2 = r2_;
      OwnValues_(f.vd1, D * f.nd1, u1);
      if (f.kind == 0) { OwnValues_(f.vd2, D * f.nd2, u2); }
      if (f.kind == 1) { NbrValues_(f.vd2, D * f.nd2, u2); }
      const int nsides = (f.kind == 2) ? 1 : 2;
      r1.SetSize(D * f.nd1);
      r2.SetSize(D * std::max(f.nd2, 1));
      r1 = 0.0;
      r2 = 0.0;
      for (int q = 0; q < f.nq; q++)
      {
         const FacePt &fp = fpt_[f.p0 + q];
         real_t a1[3], a2[3] = {0, 0, 0}, jmp[3];
         ValAt_(D, f.nd1, u1, &fshape_[fp.shape[0]], a1);
         if (nsides == 2) { ValAt_(D, f.nd2, u2, &fshape_[fp.shape[1]], a2); }
         else
         {
            // wall: [w] = w - (g(t) - g(t_s)) (dynamic) or u - g (static)
            Vector xv(const_cast<real_t *>(fp.x), D), g(D);
            if (jump_from_gs) { for (int c = 0; c < D; c++) { a2[c] = fp.gs[c]; } }
            else
            {
               elast_.EvalDirichlet(f.attr, xv, t, g);
               for (int c = 0; c < D; c++) { a2[c] = g(c); }
               if (dyn_jump_)
               {
                  for (int c = 0; c < D; c++) { a2[c] -= fp.gs[c]; }
               }
            }
         }
         for (int c = 0; c < D; c++) { jmp[c] = a1[c] - a2[c]; }
         real_t J[3][3] = {{0}};
         for (int i = 0; i < D; i++)
            for (int l = 0; l < D; l++) { J[i][l] = jmp[i] * fp.nor[l]; }
         real_t dcs[2][3][3], dss[2][3][3];
         for (int s = 0; s < nsides; s++)
         {
            real_t Gs[3][3];
            GradAt_(D, s == 0 ? f.nd1 : f.nd2, s == 0 ? u1 : u2, &fdshape_[fp.dshape[s]], Gs);
            delta(fp, s, Gs, J, dcs[s], dss[s]);
         }
         // consistency: -{Delta sigma} nor . [v];  symmetry: eps c {grad v} : dsym
         real_t avg_n[3];
         for (int i = 0; i < D; i++)
         {
            real_t v = 0.0;
            for (int l = 0; l < D; l++)
            {
               v += (nsides == 2 ? 0.5 * (dcs[0][i][l] + dcs[1][i][l]) : dcs[0][i][l])
                    * fp.nor[l];
            }
            avg_n[i] = v;
         }
         const real_t csym = (nsides == 2) ? 0.5 * eps_ : eps_;
         const int nlocal = (f.kind == 1) ? 1 : nsides;
         for (int s = 0; s < nlocal; s++)
         {
            const int nd = (s == 0) ? f.nd1 : f.nd2;
            const real_t *sh = &fshape_[fp.shape[s]];
            const real_t *dsh = &fdshape_[fp.dshape[s]];
            const real_t sgn = (s == 0) ? -1.0 : 1.0;   // [v] = v1 - v2
            Vector &r = (s == 0) ? r1 : r2;
            for (int k = 0; k < nd; k++)
               for (int p = 0; p < D; p++)
               {
                  real_t v = sgn * sh[k] * avg_n[p];
                  real_t sy = 0.0;
                  for (int m = 0; m < D; m++) { sy += dsh[k * D + m] * dss[s][p][m]; }
                  v += csym * sy;
                  r(p * nd + k) += fp.w * v;
               }
         }
      }
      AddElem_(f.vd1, D * f.nd1, r1, out);
      if (f.kind == 0) { AddElem_(f.vd2, D * f.nd2, r2, out); }
   }

   bool dyn_jump_ = false;   // FaceInputs_: subtract g(t_s) from g(t)

   /// Every pipeline steps with the gathered w at time t, and fcache_ =
   /// B + E of the new state (the work of ForceFromState_, in the same
   /// order: each volume point's force right after its step, then each
   /// face's force after the face's step; bitwise the same vector).
   template <int D>
   void UpdateAll_(real_t t)
   {
      Vector &ue = ue_, &r = r_;
      real_t G[3][3], wv[3], ds[3][3], E[3];
      fcache_.SetSize(us_.Size());
      fcache_ = 0.0;
      for (const VolElem &ve : vol_)
      {
         const int nd = ve.nd;
         OwnValues_(ve.vd, D * nd, ue);
         r.SetSize(D * nd);
         r = 0.0;
         for (int q = 0; q < ve.nq; q++)
         {
            const VolPt &p = vpt_[ve.p0 + q];
            if (p.nact == 0) { continue; }
            const BandPt &bp = bpt_[p.band];
            GradAt_(D, nd, ue, &vdshape_[p.dshape], G);
            StepPipe_<D>(bp, ve.lam, ve.mu, G, &state_[p.st_g], ds);
            ValAt_(D, nd, ue, &vshape_[p.shape], wv);
            StepMass_<D>(bp, wv, &state_[p.st_m]);
            AddVol_(D, nd, &vdshape_[p.dshape], p.wdet, ds, r);
            MassE_<D>(bp, &state_[p.st_m], E);
            const real_t *sh = &vshape_[p.shape];
            for (int k = 0; k < nd; k++)
               for (int c = 0; c < D; c++) { r(c * nd + k) += p.wdet * ve.rho * E[c] * sh[k]; }
            if (beta_over_mu_ > 0.0) { AddBetaChi_(ve, p, nd, r); }
         }
         AddElem_(ve.vd, D * nd, r, fcache_);
      }
      // faces: step both pipelines of each side; their Delta sigma of the
      // new state go straight into the face's test terms (ForceFace_'s
      // assembly, the same loop in FaceInputs_)
      dyn_jump_ = true;
      for (const Face &f : faces_)
      {
         FaceInputs_<D>(f, t, false, [&](const FacePt &fp, int s, const real_t Gs[3][3],
                                            const real_t J[3][3], real_t dc[3][3],
                                            real_t dsym[3][3])
         {
            const BandPt &bp = bpt_[fp.band];
            for (int i = 0; i < 3; i++)
               for (int j = 0; j < 3; j++) { dc[i][j] = dsym[i][j] = 0.0; }
            if (bp.nact == 0) { return; }
            StepPipe_<D>(bp, f.lam[s], f.mu[s], Gs, &state_[fp.st_c[s]], dc);
            StepPipe_<D>(bp, f.lam[s], f.mu[s], J, &state_[fp.st_s[s]], dsym);
         }, fcache_);
      }
      dyn_jump_ = false;
   }

   /// The force of one volume point from its stored state into r.
   void VolForce_(const VolElem &ve, const VolPt &p, int nd, real_t ds[3][3], real_t E[3],
                  Vector &r) const
   {
      const int D = D_;
      PipeDelta_(bpt_[p.band], ve.lam, ve.mu, &state_[p.st_g], ds);
      AddVol_(D, nd, &vdshape_[p.dshape], p.wdet, ds, r);
      if (D_ == 2) { MassE_<2>(bpt_[p.band], &state_[p.st_m], E); }
      else { MassE_<3>(bpt_[p.band], &state_[p.st_m], E); }
      const real_t *sh = &vshape_[p.shape];
      for (int k = 0; k < nd; k++)
         for (int c = 0; c < D; c++) { r(c * nd + k) += p.wdet * ve.rho * E[c] * sh[k]; }
      if (beta_over_mu_ > 0.0) { AddBetaChi_(ve, p, nd, r); }
   }

   /// 2.5D term, memory part: r += int beta (y_J - w) . v at one point.
   void AddBetaChi_(const VolElem &ve, const VolPt &p, int nd, Vector &r) const
   {
      real_t X[3];
      MassChi_(bpt_[p.band], &state_[p.st_m], X);
      const real_t s = p.wdet * beta_over_mu_ * ve.mu;
      const real_t *sh = &vshape_[p.shape];
      for (int k = 0; k < nd; k++)
         for (int c = 0; c < D_; c++) { r(c * nd + k) += s * X[c] * sh[k]; }
   }

   /// B + E from the stored state into out (at Enter; UpdateAll_ forms the
   /// same vector each step).
   void ForceFromState_(Vector &out)
   {
      const int D = D_;
      Vector &r = r_;
      real_t ds[3][3], E[3];
      for (const VolElem &ve : vol_)
      {
         const int nd = ve.nd;
         r.SetSize(D * nd);
         r = 0.0;
         for (int q = 0; q < ve.nq; q++)
         {
            const VolPt &p = vpt_[ve.p0 + q];
            if (p.nact == 0) { continue; }
            VolForce_(ve, p, nd, ds, E, r);
         }
         AddElem_(ve.vd, D * nd, r, out);
      }
      // faces: Delta sigma from the stored pipelines; the inputs are not
      // needed (zero field), only the test functions
      for (const Face &f : faces_)
      {
         ForceFace_(f, out);
      }
   }

   void ForceFace_(const Face &f, Vector &out) const
   {
      const int D = D_;
      const int nsides = (f.kind == 2) ? 1 : 2;
      Vector &r1 = r1_, &r2 = r2_;
      r1.SetSize(D * f.nd1);
      r2.SetSize(D * std::max(f.nd2, 1));
      r1 = 0.0;
      r2 = 0.0;
      for (int q = 0; q < f.nq; q++)
      {
         const FacePt &fp = fpt_[f.p0 + q];
         const BandPt &bp = bpt_[fp.band];
         if (bp.nact == 0) { continue; }
         real_t dcs[2][3][3] = {{{0}}}, dss[2][3][3] = {{{0}}};
         for (int s = 0; s < nsides; s++)
         {
            PipeDelta_(bp, f.lam[s], f.mu[s], &state_[fp.st_c[s]], dcs[s]);
            PipeDelta_(bp, f.lam[s], f.mu[s], &state_[fp.st_s[s]], dss[s]);
         }
         real_t avg_n[3];
         for (int i = 0; i < D; i++)
         {
            real_t v = 0.0;
            for (int l = 0; l < D; l++)
            {
               v += (nsides == 2 ? 0.5 * (dcs[0][i][l] + dcs[1][i][l]) : dcs[0][i][l])
                    * fp.nor[l];
            }
            avg_n[i] = v;
         }
         const real_t csym = (nsides == 2) ? 0.5 * eps_ : eps_;
         const int nlocal = (f.kind == 1) ? 1 : nsides;
         for (int s = 0; s < nlocal; s++)
         {
            const int nd = (s == 0) ? f.nd1 : f.nd2;
            const real_t *sh = &fshape_[fp.shape[s]];
            const real_t *dsh = &fdshape_[fp.dshape[s]];
            const real_t sgn = (s == 0) ? -1.0 : 1.0;
            Vector &r = (s == 0) ? r1 : r2;
            for (int k = 0; k < nd; k++)
               for (int p = 0; p < D; p++)
               {
                  real_t v = sgn * sh[k] * avg_n[p];
                  real_t sy = 0.0;
                  for (int m = 0; m < D; m++) { sy += dsh[k * D + m] * dss[s][p][m]; }
                  v += csym * sy;
                  r(p * nd + k) += fp.w * v;
               }
         }
      }
      AddElem_(f.vd1, D * f.nd1, r1, out);
      if (f.kind == 0) { AddElem_(f.vd2, D * f.nd2, r2, out); }
   }
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_PML_MEMORY_HPP

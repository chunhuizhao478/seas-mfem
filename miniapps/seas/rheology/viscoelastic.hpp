// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// viscoelastic.hpp — viscoelastic attenuation: a generalized Maxwell body
// (GMB) with L relaxation mechanisms (viscoelastic attenuation plan rev6,
// sections 2 and 5.3-5.4; the rheology plan rev 1 reserved this file as
// phase R7).
//
// Stress convention (plan decision D1, section 2.7): the deck's moduli are
// the RELAXED moduli lambda_R, mu_R, so C_0 = C_R and the stiffness K is the
// elastic one.  The model adds (plan eq. 2.11 / 5.1)
//
//    Delta sigma = sum_l C_l : xi_l,   C_l : xi = lambda_l tr(xi) I + 2 mu_l xi,
//
// with the spring strain xi_l of mechanism l obeying (2.8)
//
//    d xi_l / dt + omega_l xi_l = d eps / dt,      xi_l = 0 at entry (2.42d),
//
// stepped exactly for a strain linear over the step (2.28):
//
//    xi_l^{n+1} = exp(-omega_l dt) xi_l^n + phi(omega_l dt) (eps^{n+1} - eps^n),
//    phi(x) = (1 - exp(-x)) / x.
//
// FitGMB (section 5.3) gives omega_l, Y^P_l, Y^S_l from Q_P, Q_S over the band
// [f_c / sqrt(R), f_c sqrt(R)] (2L - 1 log-spaced points, least squares,
// (2.19)-(2.20a)), the unrelaxed moduli (2.21) and lambda_l = Y^lam_l lambda_U,
// mu_l = Y^mu_l mu_U with Y^lam of (2.13).  With moduli "reference" (SeisSol's
// convention) the deck's moduli are those at f_c and (2.23) gives the
// unrelaxed ones; FitGMB then also returns the relaxed moduli.
//
// Q per element: either fixed (Q_s, Q_p) or the relation (2.17a)
// Q_S = q_s_over_vs V_S (V_S in m/s), Q_P = q_p_over_q_s Q_S, from the
// node's C_0 and density.  The per-node coefficients lambda_l, mu_l are
// part of the node state (set by InitState), so a heterogeneous mesh needs
// no lookup per step.
//
// Scope (Chunhui Zhao, 2026-10-20; plan rev6, 2026-11-17).

#ifndef MFEM_SEAS_VISCOELASTIC_HPP
#define MFEM_SEAS_VISCOELASTIC_HPP

#include "mfem.hpp"
#include "rheology_model.hpp"

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{

/// Plan phase that implements this model (named in older refusals).
constexpr const char *kViscoelasticPhase = "R7";

/// How the deck's moduli are read (plan section 2.7, decision D1).
enum class GMBModuli { Relaxed, Reference };

/// The band and the number of mechanisms (one per run).
struct GMBBand
{
   int L = 3;               ///< mechanisms
   real_t f_c = 0.1;        ///< central frequency [Hz]
   real_t R = 100.0;        ///< f_max / f_min
};

/// Result of FitGMB for one material.
struct GMBFit
{
   std::vector<real_t> omega;       ///< relaxation frequencies omega_l [rad/s]
   std::vector<real_t> YP, YS;      ///< anelastic coefficients (2.12)
   std::vector<real_t> Ylam, Ymu;   ///< (2.13)
   std::vector<real_t> lam_l, mu_l; ///< lambda_l = Ylam_l lam_U, mu_l = Ymu_l mu_U
   real_t lam_U = 0.0, mu_U = 0.0;  ///< unrelaxed moduli
   real_t lam_R = 0.0, mu_R = 0.0;  ///< relaxed moduli (input, or (2.23) with Reference)
   std::string error;               ///< "" when the fit is acceptable
};

namespace gmb
{

/// The 2L - 1 fit frequencies omega_k (2.19a), increasing; L = 1: omega_c.
inline std::vector<real_t> FitFrequencies(const GMBBand &b)
{
   const int K = 2 * b.L - 1;
   std::vector<real_t> w(K);
   const real_t wc = 2.0 * M_PI * b.f_c;
   if (K == 1) { w[0] = wc; return w; }
   const real_t wmin = wc / std::sqrt(b.R);
   for (int k = 0; k < K; k++)
   {
      w[k] = wmin * std::pow(b.R, static_cast<real_t>(k) / (K - 1));
   }
   return w;
}

/// The relaxation frequencies omega_l = omega_{2l-1} (2.19a).
inline std::vector<real_t> RelaxationFrequencies(const GMBBand &b)
{
   const std::vector<real_t> wk = FitFrequencies(b);
   std::vector<real_t> wl(b.L);
   for (int l = 0; l < b.L; l++) { wl[l] = wk[2 * l]; }
   return wl;
}

/// Least squares min |A y - r| (A: m x n row-major, m >= n) by Householder
/// QR.  Returns false when A is rank deficient (|R_jj| <= 1e-14 max |R_ii|).
inline bool LeastSquares(int m, int n, std::vector<real_t> A, std::vector<real_t> r,
                         std::vector<real_t> &y)
{
   for (int j = 0; j < n; j++)
   {
      real_t nrm = 0.0;
      for (int i = j; i < m; i++) { nrm += A[i * n + j] * A[i * n + j]; }
      nrm = std::sqrt(nrm);
      if (nrm == 0.0) { return false; }
      const real_t alpha = (A[j * n + j] > 0.0) ? -nrm : nrm;
      std::vector<real_t> v(m, 0.0);
      for (int i = j; i < m; i++) { v[i] = A[i * n + j]; }
      v[j] -= alpha;
      real_t vv = 0.0;
      for (int i = j; i < m; i++) { vv += v[i] * v[i]; }
      if (vv == 0.0) { continue; }
      for (int c = j; c < n; c++)
      {
         real_t s = 0.0;
         for (int i = j; i < m; i++) { s += v[i] * A[i * n + c]; }
         s = 2.0 * s / vv;
         for (int i = j; i < m; i++) { A[i * n + c] -= s * v[i]; }
      }
      real_t s = 0.0;
      for (int i = j; i < m; i++) { s += v[i] * r[i]; }
      s = 2.0 * s / vv;
      for (int i = j; i < m; i++) { r[i] -= s * v[i]; }
   }
   real_t dmax = 0.0;
   for (int j = 0; j < n; j++) { dmax = std::max(dmax, std::fabs(A[j * n + j])); }
   y.assign(n, 0.0);
   for (int j = n - 1; j >= 0; j--)
   {
      if (!(std::fabs(A[j * n + j]) > 1e-14 * dmax)) { return false; }
      real_t s = r[j];
      for (int c = j + 1; c < n; c++) { s -= A[j * n + c] * y[c]; }
      y[j] = s / A[j * n + j];
   }
   return true;
}

/// Y_l for one constant Q over the band: (2.20), (2.20a).
inline bool FitY(const GMBBand &b, real_t Q, std::vector<real_t> &Y)
{
   const std::vector<real_t> wk = FitFrequencies(b);
   const std::vector<real_t> wl = RelaxationFrequencies(b);
   const int K = static_cast<int>(wk.size()), L = b.L;
   std::vector<real_t> A(static_cast<std::size_t>(K) * L), r(K, 1.0 / Q);
   for (int k = 0; k < K; k++)
   {
      for (int l = 0; l < L; l++)
      {
         A[k * L + l] = (wl[l] * wk[k] + wl[l] * wl[l] / Q) / (wl[l] * wl[l] + wk[k] * wk[k]);
      }
   }
   return LeastSquares(K, L, A, r, Y);
}

/// psi_1, psi_2 of (2.22a) at angular frequency w: M(w) / M_U = psi_1 + i psi_2.
inline void Psi(const std::vector<real_t> &wl, const std::vector<real_t> &Y, real_t w,
                real_t &psi1, real_t &psi2)
{
   psi1 = 1.0;
   psi2 = 0.0;
   for (std::size_t l = 0; l < wl.size(); l++)
   {
      const real_t d = wl[l] * wl[l] + w * w;
      psi1 -= Y[l] * wl[l] * wl[l] / d;
      psi2 += Y[l] * wl[l] * w / d;
   }
}

/// Q^{-1}(w) = Im M / Re M of the fitted body (2.16a).
inline real_t InverseQ(const std::vector<real_t> &wl, const std::vector<real_t> &Y, real_t w)
{
   real_t p1, p2;
   Psi(wl, Y, w, p1, p2);
   return p2 / p1;
}

} // namespace gmb

/// The fit of plan section 5.3.  lam_in, mu_in: the deck's moduli of the
/// material (relaxed with GMBModuli::Relaxed, at f_c with Reference).
inline GMBFit FitGMB(real_t Q_P, real_t Q_S, const GMBBand &band, real_t lam_in, real_t mu_in,
                     GMBModuli moduli = GMBModuli::Relaxed)
{
   GMBFit f;
   std::ostringstream err;
   if (!(band.L >= 1 && band.L <= 10))
   {
      err << "FitGMB: mechanisms = " << band.L << " (1 to 10)";
      f.error = err.str();
      return f;
   }
   if (!(band.f_c > 0.0) || !(band.R >= 1.0) || (band.L > 1 && !(band.R > 1.0)))
   {
      err << "FitGMB: f_c = " << band.f_c << " Hz, R = " << band.R
          << " (need f_c > 0 and R > 1 for more than one mechanism)";
      f.error = err.str();
      return f;
   }
   if (!(Q_P > 1.0) || !(Q_S > 1.0) || !std::isfinite(Q_P) || !std::isfinite(Q_S))
   {
      err << "FitGMB: Q_P = " << Q_P << ", Q_S = " << Q_S << " (need finite Q > 1)";
      f.error = err.str();
      return f;
   }
   if (!(mu_in > 0.0) || !(lam_in + 2.0 * mu_in > 0.0))
   {
      err << "FitGMB: lambda = " << lam_in << ", mu = " << mu_in << " Pa";
      f.error = err.str();
      return f;
   }
   const int L = band.L;
   f.omega = gmb::RelaxationFrequencies(band);
   if (!gmb::FitY(band, Q_P, f.YP) || !gmb::FitY(band, Q_S, f.YS))
   {
      f.error = "FitGMB: the least-squares system is rank deficient (band too narrow for L)";
      return f;
   }
   real_t sP = 0.0, sS = 0.0;
   for (int l = 0; l < L; l++)
   {
      if (!(f.YP[l] >= 0.0) || !(f.YS[l] >= 0.0))
      {
         err << "FitGMB: mechanism " << l + 1 << " has Y^P = " << f.YP[l] << ", Y^S = "
             << f.YS[l] << " < 0 (section 2.8: a negative spring); widen the band or "
             "change L";
         f.error = err.str();
         return f;
      }
      sP += f.YP[l];
      sS += f.YS[l];
   }
   if (!(sP < 1.0) || !(sS < 1.0))
   {
      err << "FitGMB: sum Y^P = " << sP << ", sum Y^S = " << sS
          << " (each must be < 1; Q is too low for this band)";
      f.error = err.str();
      return f;
   }
   real_t P_U, mu_U;
   if (moduli == GMBModuli::Relaxed)
   {
      f.lam_R = lam_in;
      f.mu_R = mu_in;
      P_U = (lam_in + 2.0 * mu_in) / (1.0 - sP);   // (2.21)
      mu_U = mu_in / (1.0 - sS);
   }
   else
   {
      // (2.22)-(2.23): M_U = M_ref (r + psi_1) / (2 r^2) at omega_c
      const real_t wc = 2.0 * M_PI * band.f_c;
      auto unrelaxed = [&](real_t Mref, const std::vector<real_t> &Y)
      {
         real_t p1, p2;
         gmb::Psi(f.omega, Y, wc, p1, p2);
         const real_t r = std::sqrt(p1 * p1 + p2 * p2);
         return Mref * (r + p1) / (2.0 * r * r);
      };
      P_U = unrelaxed(lam_in + 2.0 * mu_in, f.YP);
      mu_U = unrelaxed(mu_in, f.YS);
      f.mu_R = mu_U * (1.0 - sS);
      f.lam_R = P_U * (1.0 - sP) - 2.0 * f.mu_R;
   }
   f.mu_U = mu_U;
   f.lam_U = P_U - 2.0 * mu_U;
   f.Ylam.resize(L);
   f.Ymu.resize(L);
   f.lam_l.resize(L);
   f.mu_l.resize(L);
   real_t skap = 0.0;
   for (int l = 0; l < L; l++)
   {
      f.Ymu[l] = f.YS[l];
      // (2.12b): lambda_l = Y^P_l (lam_U + 2 mu_U) - 2 Y^S_l mu_U
      f.lam_l[l] = f.YP[l] * P_U - 2.0 * f.YS[l] * mu_U;
      f.mu_l[l] = f.YS[l] * mu_U;
      f.Ylam[l] = (f.lam_U != 0.0) ? f.lam_l[l] / f.lam_U : 0.0;   // (2.13)
      const real_t kap_l = f.lam_l[l] + 2.0 / 3.0 * f.mu_l[l];
      if (!(kap_l >= 0.0))
      {
         err << "FitGMB: mechanism " << l + 1 << " has a negative bulk part kappa_l = "
             << kap_l << " Pa (section 2.8): Q_P / Q_S too low for this lambda / mu";
         f.error = err.str();
         return f;
      }
      skap += kap_l;
   }
   const real_t kap_U = f.lam_U + 2.0 / 3.0 * f.mu_U;
   if (!(skap < kap_U))
   {
      err << "FitGMB: sum kappa_l = " << skap << " >= kappa_U = " << kap_U
          << " Pa (section 2.8: the relaxed bulk modulus would be <= 0)";
      f.error = err.str();
      return f;
   }
   return f;
}

/// How Q is given (deck: fixed values or the relation (2.17a)).
struct GMBQuality
{
   bool from_vs = true;          ///< Q_S = q_s_over_vs V_S, Q_P = q_p_over_q_s Q_S
   real_t q_s_over_vs = 0.05;    ///< [s/m] (V_S in m/s)
   real_t q_p_over_q_s = 1.5;
   real_t Q_s = 0.0, Q_p = 0.0;  ///< fixed values (from_vs = false)

   void At(real_t lam, real_t mu, real_t rho, real_t &QP, real_t &QS) const
   {
      if (from_vs)
      {
         QS = q_s_over_vs * std::sqrt(mu / rho);
         QP = q_p_over_q_s * QS;
      }
      else
      {
         QS = Q_s;
         QP = Q_p;
      }
      (void)lam;
   }
};

/// The pointwise GMB (plan section 5.4).  State per node, in this order:
///    xi_1 .. xi_L (6 each), eps^n (6), lambda_1 .. lambda_L, mu_1 .. mu_L
/// so NumStateVars = 6 L + 6 + 2 L.  The coefficients are written by
/// InitState from the node's C_0 (= C_R) and density.
class ViscoelasticModel : public IRheologyModel
{
public:
   ViscoelasticModel(const GMBBand &band, const GMBQuality &q,
                     GMBModuli moduli = GMBModuli::Relaxed)
      : band_(band), q_(q), moduli_(moduli), omega_(gmb::RelaxationFrequencies(band))
   {
      MFEM_VERIFY(moduli == GMBModuli::Relaxed,
                  "ViscoelasticModel: moduli = \"reference\" needs the relaxed moduli of "
                  "FitGMB in C_0 and K (plan section 2.7); not wired in a run yet");
      MFEM_VERIFY(band.L >= 1, "ViscoelasticModel: L = " << band.L);
   }

   const char *Name() const override { return "viscoelastic"; }
   int NumStateVars() const override { return 8 * band_.L + 6; }
   bool IsLinearInStrain() const override { return true; }
   int L() const { return band_.L; }
   const std::vector<real_t> &Omega() const { return omega_; }
   const GMBBand &Band() const { return band_; }

   /// The fit at one material (for diagnostics and the step factor).
   GMBFit Fit(real_t lam, real_t mu, real_t rho) const
   {
      real_t QP, QS;
      q_.At(lam, mu, rho, QP, QS);
      return FitGMB(QP, QS, band_, lam, mu, moduli_);
   }

   /// xi = 0, eps^n = 0 (SetEntryStrain sets it); coefficients from the
   /// fit at the node's material.
   void InitState(const RheologyPointContext &ctx, real_t *q) const override
   {
      MFEM_VERIFY(ctx.rho0 > 0.0, "ViscoelasticModel: the node has no density (rho0 = "
                  << ctx.rho0 << ")");
      const GMBFit f = Fit(ctx.lambda0, ctx.mu0, ctx.rho0);
      MFEM_VERIFY(f.error.empty(), f.error << " (node at " << ctx.x[0] << ", " << ctx.x[1]
                  << ", " << ctx.x[2] << ")");
      const int L = band_.L;
      for (int a = 0; a < 6 * L; a++) { q[a] = 0.0; }
      for (int a = 0; a < 6; a++) { q[6 * L + a] = 0.0; }
      for (int l = 0; l < L; l++)
      {
         q[6 * L + 6 + l] = f.lam_l[l];
         q[7 * L + 6 + l] = f.mu_l[l];
      }
   }

   /// eps^n of a node state = the entry strain (2.42d), so the first step's
   /// increment is that of the dynamic motion only.
   void SetEntryStrain(const RheologyPointContext &, const real_t eps[6],
                       real_t *q) const override
   {
      for (int a = 0; a < 6; a++) { q[6 * band_.L + a] = eps[a]; }
   }

   /// One step (5.3) and (2.29).  The factors exp(-omega_l dt) and
   /// phi(omega_l dt) are cached for the last dt.
   bool Update(const RheologyPointContext &, const real_t eps[6], real_t dt, real_t *q,
               real_t dsig[6]) const override
   {
      const int L = band_.L;
      if (dt != dt_cached_) { Factors_(dt); }
      real_t *eold = q + 6 * L;
      real_t de[6];
      for (int a = 0; a < 6; a++) { de[a] = eps[a] - eold[a]; }
      for (int l = 0; l < L; l++)
      {
         real_t *xi = q + 6 * l;
         const real_t al = a_[l], bl = b_[l];
         for (int a = 0; a < 6; a++) { xi[a] = al * xi[a] + bl * de[a]; }
         const real_t lam = q[6 * L + 6 + l], mu = q[7 * L + 6 + l];
         const real_t ltr = lam * (xi[0] + xi[1] + xi[2]);
         dsig[0] += ltr + 2.0 * mu * xi[0];
         dsig[1] += ltr + 2.0 * mu * xi[1];
         dsig[2] += ltr + 2.0 * mu * xi[2];
         dsig[3] += 2.0 * mu * xi[3];
         dsig[4] += 2.0 * mu * xi[4];
         dsig[5] += 2.0 * mu * xi[5];
      }
      for (int a = 0; a < 6; a++) { eold[a] = eps[a]; }
      return true;
   }

   /// C_U / C_R of the P modulus, 1 / (1 - sum Y^P) (see (2.31)): above 1;
   /// the solver cuts its step by 1 / sqrt(this) (plan section 5.9).
   real_t StiffnessBound(const RheologyPointContext &ctx, const real_t *q) const override
   {
      const int L = band_.L;
      real_t sl = 0.0;
      for (int l = 0; l < L; l++) { sl += q[6 * L + 6 + l] + 2.0 * q[7 * L + 6 + l]; }
      const real_t PR = ctx.lambda0 + 2.0 * ctx.mu0;
      return (PR + sl) / PR;
   }

   int NumDiagnostics() const override { return 2; }
   const char *DiagnosticName(int i) const override
   {
      return i == 0 ? "visco_dsig_vm" : "visco_dissipation";
   }
   /// 0: von Mises norm of Delta sigma; 1: dissipation rate D of (2.26),
   /// sum_l omega_l xi_l : C_l : xi_l [W/m^3].
   void Diagnostics(const RheologyPointContext &, const real_t *q, real_t *out) const override
   {
      const int L = band_.L;
      real_t ds[6] = {0, 0, 0, 0, 0, 0};
      real_t D = 0.0;
      for (int l = 0; l < L; l++)
      {
         const real_t *xi = q + 6 * l;
         const real_t lam = q[6 * L + 6 + l], mu = q[7 * L + 6 + l];
         real_t s[6];
         sym6::ApplyIsotropic(lam, mu, xi, s);
         for (int a = 0; a < 6; a++) { ds[a] += s[a]; }
         D += omega_[l] * sym6::DoubleDot(xi, s);
      }
      const real_t m = sym6::Trace(ds) / 3.0;
      real_t dev[6] = {ds[0] - m, ds[1] - m, ds[2] - m, ds[3], ds[4], ds[5]};
      out[0] = std::sqrt(1.5 * sym6::DoubleDot(dev, dev));
      out[1] = D;
   }

private:
   void Factors_(real_t dt) const
   {
      const int L = band_.L;
      a_.resize(L);
      b_.resize(L);
      for (int l = 0; l < L; l++)
      {
         const real_t x = omega_[l] * dt;
         a_[l] = std::exp(-x);
         b_[l] = (x > 0.0) ? -std::expm1(-x) / x : 1.0;
      }
      dt_cached_ = dt;
   }

   GMBBand band_;
   GMBQuality q_;
   GMBModuli moduli_;
   std::vector<real_t> omega_;
   mutable std::vector<real_t> a_, b_;
   mutable real_t dt_cached_ = -1.0;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_VISCOELASTIC_HPP

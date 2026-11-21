// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// strain_split.hpp — degraded small-strain isotropic elasticity with a
// tension/compression split, at one quadrature point.
//
// Port of SmallDeformationIsotropicElasticity::computeStress{No,Spectral,VolDev}
// Decomposition (AD, spe_hftc) and NDSmallDeformationIsotropicElasticity
// (non-AD, cmame_revision); both evaluate the same formulas:
//
//   spectral   lambda = K - 2G/LIBMESH_DIM with LIBMESH_DIM = 3, also in 2D:
//              the 2D models are plane strain with eps_zz = 0, so callers pass
//              the strain padded to 3x3.
//              eps+ = sum_i <eps_i>_+ n_i (x) n_i
//              sigma+ = lambda <tr eps>_+ I + 2G eps+
//              sigma  = g sigma+ + (sigma0 - sigma+),  sigma0 = K tr eps I + 2G dev eps
//              psi+ = 1/2 lambda <tr eps>_+^2 + G eps+ : eps+
//              psi- = (1/2 lambda (tr eps)^2 + G eps : eps) - psi+
//   voldev     sigma- = K <tr eps>_- I,  sigma = g (sigma0 - sigma-) + sigma-
//              psi- = 1/2 K <tr eps>_-^2,  psi+ = psi0 - psi-
//   none       sigma = g sigma0,  psi+ = 1/2 sigma0 : eps,  psi- = 0
//
// sigma is continuous across repeated eigenvalues. SplitTangentApply gives
// the consistent tangent (Newton in the quasi-static phase R).
//
// Storage: 3x3 column-major arrays (MFEM DenseMatrix layout).

#ifndef MFEM_PFF_STRAIN_SPLIT_HPP
#define MFEM_PFF_STRAIN_SPLIT_HPP

#include "mfem.hpp"
#include "../../../linalg/kernels.hpp" // mfem::kernels::CalcEigenvalues<3>
#include "../common/pff_types.hpp"

namespace mfem
{
namespace pff
{

/// Result of the split at one point.
struct SplitResult
{
   real_t sigma[9];        ///< degraded Cauchy stress, 3x3 column-major
   real_t psi_active;      ///< psi_e^+ (undegraded)
   real_t psi_inactive;    ///< psi_e^-
   real_t psi;             ///< g psi+ + psi-
};

MFEM_HOST_DEVICE inline real_t MacaulayPos(real_t x) { return x > 0.0 ? x : 0.0; }
MFEM_HOST_DEVICE inline real_t MacaulayNeg(real_t x) { return x < 0.0 ? x : 0.0; }

/// Positive spectral part eps+ = sum_i <eps_i>_+ n_i n_i of a symmetric 3x3.
MFEM_HOST_DEVICE inline void SpectralPositivePart(const real_t *eps, real_t *eps_pos)
{
   real_t lam[3], vec[9];
   kernels::CalcEigenvalues<3>(eps, lam, vec);
   for (int i = 0; i < 9; i++) { eps_pos[i] = 0.0; }
   for (int k = 0; k < 3; k++)
   {
      const real_t lp = MacaulayPos(lam[k]);
      if (lp == 0.0) { continue; }
      const real_t *n = vec + 3*k;
      for (int j = 0; j < 3; j++)
         for (int i = 0; i < 3; i++)
         {
            eps_pos[i + 3*j] += lp * n[i] * n[j];
         }
   }
}

/// Degraded stress and energies for strain eps (3x3, symmetric).
MFEM_HOST_DEVICE inline void ComputeSplit(Decomposition dec, const real_t *eps,
                                          real_t K, real_t G, real_t g,
                                          SplitResult &r)
{
   const real_t tr = eps[0] + eps[4] + eps[8];
   real_t sigma0[9];
   for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
      {
         const real_t dij = (i == j) ? 1.0 : 0.0;
         // K tr I + 2G dev(eps) = (K - 2G/3) tr I + 2G eps
         sigma0[i + 3*j] = K * tr * dij + 2.0 * G * (eps[i + 3*j] - dij * tr / 3.0);
      }
   real_t eps_eps = 0.0;
   for (int i = 0; i < 9; i++) { eps_eps += eps[i] * eps[i]; }

   if (dec == Decomposition::None)
   {
      real_t s0e = 0.0;
      for (int i = 0; i < 9; i++) { s0e += sigma0[i] * eps[i]; r.sigma[i] = g * sigma0[i]; }
      r.psi_active = 0.5 * s0e;
      r.psi_inactive = 0.0;
   }
   else if (dec == Decomposition::Spectral)
   {
      const real_t lambda = K - 2.0 * G / 3.0;
      const real_t trp = MacaulayPos(tr);
      real_t ep[9];
      SpectralPositivePart(eps, ep);
      real_t ep_ep = 0.0;
      for (int j = 0; j < 3; j++)
         for (int i = 0; i < 3; i++)
         {
            const int ij = i + 3*j;
            const real_t sp = lambda * trp * ((i == j) ? 1.0 : 0.0) + 2.0 * G * ep[ij];
            r.sigma[ij] = g * sp + (sigma0[ij] - sp);
            ep_ep += ep[ij] * ep[ij];
         }
      const real_t psi0 = 0.5 * lambda * tr * tr + G * eps_eps;
      r.psi_active = 0.5 * lambda * trp * trp + G * ep_ep;
      r.psi_inactive = psi0 - r.psi_active;
   }
   else // VolDev
   {
      const real_t trn = MacaulayNeg(tr);
      real_t dev_dev = 0.0;
      for (int j = 0; j < 3; j++)
         for (int i = 0; i < 3; i++)
         {
            const int ij = i + 3*j;
            const real_t dij = (i == j) ? 1.0 : 0.0;
            const real_t dev = eps[ij] - dij * tr / 3.0;
            dev_dev += dev * dev;
            const real_t sn = K * trn * dij;
            r.sigma[ij] = g * (sigma0[ij] - sn) + sn;
         }
      const real_t psi0 = 0.5 * K * tr * tr + G * dev_dev;
      r.psi_inactive = 0.5 * K * trn * trn;
      r.psi_active = psi0 - r.psi_inactive;
   }
   r.psi = g * r.psi_active + r.psi_inactive;
}

/// Directional derivative of the degraded stress: dsig = (d sigma / d eps) : deps
/// for a symmetric strain increment deps (3x3), with ComputeSplit's formulas.
///   spectral: sigma = sigma0 - (1 - g) sigma+,
///     d sigma+ = lambda H(tr eps) tr(deps) I + 2G d eps+,
///     d eps+ = sum_a H(l_a) (n_a.deps.n_a) M_a
///            + 1/2 sum_{a != b} theta_ab (n_a.deps.n_b) (n_a n_b + n_b n_a),
///     theta_ab = (<l_a>+ - <l_b>+)/(l_a - l_b), and H(l_a) when l_a = l_b
///     (the limit; H(0) = 0). MOOSE forms the same derivative by AD.
///   voldev: sigma = g (sigma0 - sigma-) + sigma-, d sigma- = K H(-tr) tr(deps) I.
///   none: g C : deps.
MFEM_HOST_DEVICE inline void SplitTangentApply(Decomposition dec, const real_t *eps,
                                               real_t K, real_t G, real_t g,
                                               const real_t *deps, real_t *dsig)
{
   const real_t lambda = K - 2.0 * G / 3.0;
   const real_t tr = eps[0] + eps[4] + eps[8];
   const real_t dtr = deps[0] + deps[4] + deps[8];
   real_t dsig0[9];
   for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
      {
         dsig0[i + 3*j] = lambda * dtr * ((i == j) ? 1.0 : 0.0) + 2.0 * G * deps[i + 3*j];
      }
   if (dec == Decomposition::None)
   {
      for (int i = 0; i < 9; i++) { dsig[i] = g * dsig0[i]; }
      return;
   }
   if (dec == Decomposition::VolDev)
   {
      const real_t kn = (tr < 0.0) ? K * dtr : 0.0;
      for (int j = 0; j < 3; j++)
         for (int i = 0; i < 3; i++)
         {
            const real_t dn = (i == j) ? kn : 0.0;
            dsig[i + 3*j] = g * (dsig0[i + 3*j] - dn) + dn;
         }
      return;
   }
   real_t lam[3], vec[9];
   kernels::CalcEigenvalues<3>(eps, lam, vec);
   // q_ab = n_a . deps . n_b
   real_t q[9];
   for (int a = 0; a < 3; a++)
      for (int b = 0; b < 3; b++)
      {
         real_t v = 0.0;
         for (int k = 0; k < 3; k++)
            for (int l = 0; l < 3; l++) { v += vec[k + 3*a] * deps[k + 3*l] * vec[l + 3*b]; }
         q[a + 3*b] = v;
      }
   const real_t scale = fabs(lam[0]) + fabs(lam[1]) + fabs(lam[2]);
   real_t dep[9];
   for (int i = 0; i < 9; i++) { dep[i] = 0.0; }
   for (int a = 0; a < 3; a++)
      for (int b = 0; b < 3; b++)
      {
         real_t c;
         if (a == b) { c = (lam[a] > 0.0) ? q[a + 3*a] : 0.0; }
         else
         {
            const real_t dl = lam[a] - lam[b];
            const real_t th = (fabs(dl) > 1e-12 * scale && dl != 0.0)
                              ? (MacaulayPos(lam[a]) - MacaulayPos(lam[b])) / dl
                              : ((lam[a] > 0.0) ? 1.0 : 0.0);
            c = 0.5 * th * q[a + 3*b];
         }
         if (c == 0.0) { continue; }
         const real_t *na = vec + 3*a, *nb = vec + 3*b;
         for (int j = 0; j < 3; j++)
            for (int i = 0; i < 3; i++)
            {
               dep[i + 3*j] += (a == b) ? c * na[i] * na[j] : c * (na[i] * nb[j] + nb[i] * na[j]);
            }
      }
   const real_t ht = (tr > 0.0) ? dtr : 0.0;
   for (int j = 0; j < 3; j++)
      for (int i = 0; i < 3; i++)
      {
         const int ij = i + 3*j;
         const real_t dsp = lambda * ht * ((i == j) ? 1.0 : 0.0) + 2.0 * G * dep[ij];
         dsig[ij] = dsig0[ij] - (1.0 - g) * dsp;
      }
}

/// Pad a dim x dim symmetric strain (column-major) to 3x3 (plane strain in 2D).
MFEM_HOST_DEVICE inline void PadStrain3(int dim, const real_t *e, real_t *e3)
{
   for (int i = 0; i < 9; i++) { e3[i] = 0.0; }
   for (int j = 0; j < dim; j++)
      for (int i = 0; i < dim; i++) { e3[i + 3*j] = e[i + dim*j]; }
}

} // namespace pff
} // namespace mfem

#endif

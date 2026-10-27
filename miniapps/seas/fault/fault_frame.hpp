// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// fault_frame.hpp — Nitsche dynamic rupture, plan item N-09 (fault
// coordinate spec rev 8.3, section 4) and the kernel dictionary of plan
// Part IV.  Pure functions: no mesh ownership, nothing stored.
//
//   F1  s = (up x n) / |up x n|,  d = s x n.
//   F2  tangent1 = d (dip), tangent2 = s (strike)  (Tandem ordering).
//   F3  (d, s, n) right-handed and orthonormal.
//   F5  F1 is used where |up x n| >= eps_deg = c kappa_F eps_mach / eps_ang.
//   F6  below eps_deg: edge fallback by the lexicographically smallest
//       sorted vertex-key pair, directed small -> large key, projected into
//       the tangent plane, |p_F| >= |e_F|/2 guard, next edge otherwise.
//   F7  the caller counts faces with Frame::fallback.
//   F8  (transport) is not implemented; not needed for TPV102.
//
//   Crossings (spec A3, C3):
//     ToLocal(T, sigma0_n, n, d, s) -> X_i = (T + sigma0_n n).t_i,
//                                      sigma_n_eff = -(T + sigma0_n n).n
//     Embed(V1, V2, d, s)           -> V = V1 d + V2 s
//
//   Dictionary (plan Part IV): the kernel stores (n_k, t1_k, t2_k) =
//   eps_f (n, d, s) with eps_f = sign(n . n_k) on a planar fault.
//     ToKernel:   d_i^k = eps_f delta . t_i^k
//     FromKernel: T = eps_f (tau1^k t1^k + tau2^k t2^k - sigma_n n_k)
//   Not wired into the Nitsche operator until P6.
//
// Scope (Chunhui Zhao, 2026-10-14): tetrahedral meshes, planar fault faces.

#ifndef MFEM_SEAS_FAULT_FRAME_HPP
#define MFEM_SEAS_FAULT_FRAME_HPP

#include "mfem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace mfem
{
namespace seas
{
namespace fault_frame
{

inline void Cross(const real_t a[3], const real_t b[3], real_t c[3])
{
   c[0] = a[1] * b[2] - a[2] * b[1];
   c[1] = a[2] * b[0] - a[0] * b[2];
   c[2] = a[0] * b[1] - a[1] * b[0];
}

inline real_t Dot(const real_t a[3], const real_t b[3])
{
   return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/// The constitutive frame at a point.
struct Frame
{
   real_t n[3] = {0, 0, 0};   ///< oriented unit normal (input)
   real_t d[3] = {0, 0, 0};   ///< tangent1 = dip
   real_t s[3] = {0, 0, 0};   ///< tangent2 = strike
   bool fallback = false;     ///< F6 used (components are frame, not dip/strike)
};

/// F5 threshold eps_deg = c kappa_F eps_mach / eps_ang.  With c = 45,
/// kappa_F = 1 and eps_ang = 1e-6 this is ~1e-8, the spec's value for
/// shape-regular faces in double precision (delta ~ 1e-14).
inline real_t DegeneracyThreshold(real_t kappa_F, real_t eps_ang = 1e-6,
                                  real_t c = 45.0)
{
   return c * kappa_F * std::numeric_limits<real_t>::epsilon() / eps_ang;
}

/// Face data needed by the F6 fallback: the face's vertex keys
/// (partition-invariant) and coordinates.  nv = 3 for a triangle.
struct FaceVertices
{
   int nv = 3;
   long long key[4] = {0, 0, 0, 0};
   real_t x[4][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
};

/// F1-F3 with the F5 threshold and the F6 fallback.  `n` must be the
/// oriented unit normal.  Returns false (degenerate geometry, F6(b)) when
/// every edge fails the guard; the caller refuses the face.
inline bool ComputeFrame(const real_t n[3], const real_t up[3], real_t eps_deg,
                         const FaceVertices &fv, Frame &fr)
{
   for (int i = 0; i < 3; i++) { fr.n[i] = n[i]; }
   real_t s[3];
   Cross(up, n, s);
   const real_t sl = std::sqrt(Dot(s, s));
   if (sl >= eps_deg)
   {
      for (int i = 0; i < 3; i++) { fr.s[i] = s[i] / sl; }
      Cross(fr.s, n, fr.d);
      fr.fallback = false;
      return true;
   }
   // F6(a): edges sorted by their sorted key pair, directed small -> large.
   std::vector<std::pair<std::pair<long long, long long>, std::pair<int, int>>> edges;
   for (int i = 0; i < fv.nv; i++)
   {
      const int j = (i + 1) % fv.nv;
      int a = i, b = j;
      if (fv.key[b] < fv.key[a]) { std::swap(a, b); }
      edges.push_back({{fv.key[a], fv.key[b]}, {a, b}});
   }
   std::sort(edges.begin(), edges.end());
   for (const auto &e : edges)
   {
      const int a = e.second.first, b = e.second.second;
      real_t ev[3], p[3];
      for (int i = 0; i < 3; i++) { ev[i] = fv.x[b][i] - fv.x[a][i]; }
      const real_t en = Dot(ev, n);
      for (int i = 0; i < 3; i++) { p[i] = ev[i] - en * n[i]; }
      const real_t pl = std::sqrt(Dot(p, p)), el = std::sqrt(Dot(ev, ev));
      if (pl >= 0.5 * el)                                     // F6(b) guard
      {
         for (int i = 0; i < 3; i++) { fr.s[i] = p[i] / pl; }  // F6(c)
         Cross(fr.s, n, fr.d);
         fr.fallback = true;
         return true;
      }
   }
   return false;
}

/// Global -> local before the constitutive evaluation (A3, C3).
inline void ToLocal(const real_t T[3], real_t sigma0_n, const real_t n[3],
                    const real_t d[3], const real_t s[3],
                    real_t &X1, real_t &X2, real_t &sigma_n_eff)
{
   real_t Tt[3];
   for (int i = 0; i < 3; i++) { Tt[i] = T[i] + sigma0_n * n[i]; }
   X1 = Dot(Tt, d);
   X2 = Dot(Tt, s);
   sigma_n_eff = -Dot(Tt, n);
}

/// Local -> global after the constitutive evaluation (A3).
inline void Embed(real_t V1, real_t V2, const real_t d[3], const real_t s[3],
                  real_t V[3])
{
   for (int i = 0; i < 3; i++) { V[i] = V1 * d[i] + V2 * s[i]; }
}

// --------------------------------------------------------------------------
// Dictionary (plan Part IV).
// --------------------------------------------------------------------------

/// eps_f = sign(n . n_k) for the oriented normal n and the kernel's stored
/// normal n_k of the same face.  A comparison of two normals of one face
/// (harmless under spec P1).  Aborts if the two are not (anti)parallel.
inline real_t EpsilonF(const real_t n[3], const real_t n_k[3])
{
   const real_t c = Dot(n, n_k);
   MFEM_VERIFY(std::abs(c) > 0.5,
               "fault_frame::EpsilonF: oriented normal and kernel normal are "
               "not parallel (n.n_k = " << c << ")");
   return (c > 0.0) ? 1.0 : -1.0;
}

/// Into the kernel: d_i^k = eps_f delta . t_i^k.
inline void ToKernel(const real_t delta[3], real_t eps_f, const real_t t1k[3],
                     const real_t t2k[3], real_t &d1k, real_t &d2k)
{
   d1k = eps_f * Dot(delta, t1k);
   d2k = eps_f * Dot(delta, t2k);
}

/// Out of the kernel: T = eps_f (tau1^k t1^k + tau2^k t2^k - sigma_n n_k).
inline void FromKernel(real_t tau1k, real_t tau2k, real_t sigma_n, real_t eps_f,
                       const real_t nk[3], const real_t t1k[3],
                       const real_t t2k[3], real_t T[3])
{
   for (int i = 0; i < 3; i++)
   {
      T[i] = eps_f * (tau1k * t1k[i] + tau2k * t2k[i] - sigma_n * nk[i]);
   }
}

} // namespace fault_frame
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FAULT_FRAME_HPP

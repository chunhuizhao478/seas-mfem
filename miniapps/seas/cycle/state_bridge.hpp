// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// state_bridge.hpp — fully dynamic BP5 cycle (plan rev 3, Sections 3.3 f
// and 11.4): the map between the two regimes' fault layouts.
//
// The QD kernel keeps two slip components per owned fault DOF in the kernel
// frame (t1 = dip, t2 = strike of the kernel's triad n_k, t1_k, t2_k); the
// Nitsche operator keeps the global vector delta.  With eps_f = sign(n . n_k)
// for the oriented normal n (fault/fault_frame.hpp:160-186):
//
//    s_k = eps_f delta . t_k,  k = 1, 2,      delta = eps_f (s_1 t_1 + s_2 t_2)
//
// exact up to round-off, because t_1, t_2 are orthonormal, delta has no
// normal component and eps_f^2 = 1.  The same map takes the slip rate, and,
// with -sigma_n n_k added, a traction:
//
//    T = eps_f (tau_1 t_1 + tau_2 t_2 - sigma_n n_k)       (FromKernel)
//
// psi is the same variable in both and is copied unchanged.
//
// FrameMap holds (eps_f, n_k, t1_k, t2_k) of each OWNED fault DOF as the
// dynamic operator applies it: the owner's local copy, the copy whose
// traction RestrictToOwnedFault keeps.  DynPhase builds it from the
// operator (rule M4: this file includes no solver header); the conversions
// call the dictionary functions themselves, so their doubles are the
// operator's.
//
// Ghost copies (plan R8): the QD kernel copies kernel components to a
// DOF's other local copies, the Nitsche operator copies the global vector
// and applies each copy's own triad.  The two agree when eps_f t_k of every
// copy equals the owner's; DynPhase::GhostFrameMismatch measures that and
// test_cycle_frame_map checks it on 1, 2, 4 and 8 ranks.

#ifndef MFEM_SEAS_CYCLE_STATE_BRIDGE_HPP
#define MFEM_SEAS_CYCLE_STATE_BRIDGE_HPP

#include "mfem.hpp"
#include "../fault/fault_frame.hpp"

namespace mfem
{
namespace seas
{
namespace cycle
{

class FrameMap
{
public:
   static constexpr int kPerDOF = 10;   ///< eps_f, n_k, t1_k, t2_k

   FrameMap() = default;
   /// `packed`: kPerDOF values per owned fault DOF.
   explicit FrameMap(const Vector &packed) : f_(packed)
   {
      MFEM_VERIFY(f_.Size() % kPerDOF == 0, "FrameMap: " << f_.Size()
                  << " values is not a multiple of " << kPerDOF);
   }

   int Size() const { return f_.Size() / kPerDOF; }
   real_t Eps(int j) const { return f_(kPerDOF * j); }
   const real_t *Nk(int j) const { return f_.GetData() + kPerDOF * j + 1; }
   const real_t *T1k(int j) const { return f_.GetData() + kPerDOF * j + 4; }
   const real_t *T2k(int j) const { return f_.GetData() + kPerDOF * j + 7; }
   const Vector &Packed() const { return f_; }

   /// Kernel (dip, strike) pairs (2 N) -> global three-vectors (3 N).
   void ToGlobal(const Vector &k, Vector &g) const
   {
      const int N = Size();
      MFEM_VERIFY(k.Size() == 2 * N, "FrameMap::ToGlobal: size " << k.Size()
                  << " != 2 x " << N);
      g.SetSize(3 * N);
      for (int j = 0; j < N; j++)
      {
         fault_frame::FromKernel(k(2 * j), k(2 * j + 1), 0.0, Eps(j), Nk(j), T1k(j), T2k(j),
                                 &g(3 * j));
      }
   }

   /// Global three-vectors (3 N) -> kernel (dip, strike) pairs (2 N).
   void ToKernel(const Vector &g, Vector &k) const
   {
      const int N = Size();
      MFEM_VERIFY(g.Size() == 3 * N, "FrameMap::ToKernel: size " << g.Size()
                  << " != 3 x " << N);
      k.SetSize(2 * N);
      for (int j = 0; j < N; j++)
      {
         fault_frame::ToKernel(&g(3 * j), Eps(j), T1k(j), T2k(j), k(2 * j), k(2 * j + 1));
      }
   }

   /// Kernel shear (2 N) and sigma_n (N, compression positive) -> the global
   /// traction (3 N, tension positive).
   void TractionToGlobal(const Vector &tau, const Vector &sigma_n, Vector &T) const
   {
      const int N = Size();
      MFEM_VERIFY(tau.Size() == 2 * N && sigma_n.Size() == N,
                  "FrameMap::TractionToGlobal: sizes " << tau.Size() << ", "
                  << sigma_n.Size() << " for " << N << " DOFs");
      T.SetSize(3 * N);
      for (int j = 0; j < N; j++)
      {
         fault_frame::FromKernel(tau(2 * j), tau(2 * j + 1), sigma_n(j), Eps(j), Nk(j),
                                 T1k(j), T2k(j), &T(3 * j));
      }
   }

   /// The inverse: tau_k = eps_f T . t_k, sigma_n = -eps_f T . n_k.
   void TractionToKernel(const Vector &T, Vector &tau, Vector &sigma_n) const
   {
      const int N = Size();
      MFEM_VERIFY(T.Size() == 3 * N, "FrameMap::TractionToKernel: size " << T.Size()
                  << " != 3 x " << N);
      tau.SetSize(2 * N);
      sigma_n.SetSize(N);
      for (int j = 0; j < N; j++)
      {
         fault_frame::ToKernel(&T(3 * j), Eps(j), T1k(j), T2k(j), tau(2 * j), tau(2 * j + 1));
         sigma_n(j) = -Eps(j) * fault_frame::Dot(&T(3 * j), Nk(j));
      }
   }

private:
   Vector f_;
};

/// Bridge::ToNitsche and Bridge::ToQD (plan 3.3 f) on plain vectors: the QD
/// side as kernel pairs (slip, slip rate: 2 N) and psi (N), the Nitsche side
/// as global three-vectors (3 N) and psi (N).
namespace bridge
{

inline void ToNitsche(const FrameMap &f, const Vector &slip_k, const Vector &V_k,
                      const Vector &psi_qd, Vector &d, Vector &V, Vector &psi)
{
   MFEM_VERIFY(psi_qd.Size() == f.Size(), "bridge::ToNitsche: " << psi_qd.Size()
               << " psi values for " << f.Size() << " DOFs");
   f.ToGlobal(slip_k, d);
   f.ToGlobal(V_k, V);
   psi = psi_qd;
}

inline void ToQD(const FrameMap &f, const Vector &d, const Vector &psi_dyn,
                 Vector &slip_k, Vector &psi)
{
   MFEM_VERIFY(psi_dyn.Size() == f.Size(), "bridge::ToQD: " << psi_dyn.Size()
               << " psi values for " << f.Size() << " DOFs");
   f.ToKernel(d, slip_k);
   psi = psi_dyn;
}

} // namespace bridge

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_STATE_BRIDGE_HPP

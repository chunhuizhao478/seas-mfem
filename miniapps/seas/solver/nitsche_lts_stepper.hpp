// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// nitsche_lts_stepper.hpp — Nitsche profiling plan item 7 (clustered local
// time stepping), phases P1 and P2 of the design lts-kernels-design-rev1
// (§3.4, §3.7): the operator split by element owner and the multirate
// leapfrog, with Kelvin-Voigt damping and absorbing faces per cluster (P2).
// Serial (one rank): no fault closure, no data F(t).  Not used by any run
// yet.
//
// NitscheLtsSplit: K = sum_e G_e with
//
//    G_e = vol_e + sum_{interior F of e} a_F^(e) + sum_{Dirichlet F of e} a_F,
//
// a_F^(e) the half of the face form that belongs to side e
// (DGElasticityIPCombinedIntegrator::SetSide; design §7.3), and cluster c
// owns K^(c) = sum_{e in c} G_e.  The integrators, coefficients, DOF
// numbering and face loops are those of
// ElasticityDomainOperator::AssembleStiffness (BilinearForm::Assemble), so
// sum_c K^(c) equals the assembled K up to the order of the sums.  Built on
// the Kelvin-Voigt moduli beta_e lambda_e, beta_e mu_e it splits K_eta the
// same way (design §3.4 step 2); since every term of G_e uses side e's
// moduli, G_e of K_eta is beta_e times G_e of K.  AssembleAbsorbingClusters
// gives cluster c the dashpot blocks C_e of its elements (C^(c)).
//
// NitscheMultirateLeapfrog: one sync interval [t, t + 2^L dt0], L = Nc - 1,
// of the kick-drift-kick loop of design §3.4 step 3, with dt_c = 2^c dt0 and
//
//    a_c = B_k^{-1} r_c,   r_c = -(K^(c) u + K_eta^(c) v + C^(c) v),
//
// at a tick where clusters 0 .. k end, B_k = M + 1/2 dt_{c_e} C_e on the
// elements with c_e <= k and M elsewhere (NitscheMass with the damping
// blocks of those elements and half_dt_per_elem; B_k = M without absorbing
// faces):
//
//    for m = 0 .. 2^L - 1:
//       for c with m % 2^c == 0:        v += dt_c/2 a_c       (start kicks)
//       u += dt0 v,  t += dt0                                  (drift)
//       E = {c : (m + 1) % 2^c == 0}                           (clusters ending)
//       for c in E:                     a_c from (u, v)       (all at one u, v)
//       for c in E:                     v += dt_c/2 a_c       (end kicks)
//
// The damping terms see v before the end kicks: for a cluster's own DOFs that
// is v at the middle of its step, the half-step velocity Udot* the global
// step damps with (design §3.4 step 5).  At a tick where element e's own
// cluster ends, B_e (v_new - v) = sum of the kicks gives M_e (v_new - v) =
// kicks - 1/2 dt_{c_e} C_e v_new, the global step's implicit dashpot, and the
// next start kicks reuse those a_c, so the dashpot impulse is the trapezoid
// rule on C v at e's own step boundaries.  At a tick where only faster
// clusters end, their kicks on e's rows are elastic impulses and take M_e:
// with B_e there the kick would carry an extra 1/2 dt_{c_e} C_e dv, an O(dt)
// global error (found by the order test L-10).  Init(u, v) computes
// every a_c at the start with M^{-1}, as NewmarkSolver's first Mult does.  At
// the end of the interval every cluster has ended a step, so v is
// synchronised.  With one cluster this is central difference in
// kick-drift-kick form, equal to CentralDifferenceSolver on
// NitscheWaveOperator up to round-off.

#ifndef MFEM_SEAS_NITSCHE_LTS_STEPPER_HPP
#define MFEM_SEAS_NITSCHE_LTS_STEPPER_HPP

#include "mfem.hpp"
#include "../domain/nitsche_mass.hpp"
#include "../integrator/dg_elasticity_ip_combined_integrator.hpp"

#include <cmath>
#include <map>
#include <memory>
#include <utility>
#include <vector>

namespace mfem
{
namespace seas
{

#ifdef MFEM_USE_MPI

class NitscheLtsSplit
{
public:
   /// fes: the stiffness space (ElasticityDomainOperator::GetFESpace);
   /// lambda, mu: its coefficients; dirichlet_marker: the boundary marker of
   /// its Dirichlet faces (empty: none); penalty_factor: the deck's.
   NitscheLtsSplit(ParFiniteElementSpace &fes, Coefficient &lambda, Coefficient &mu,
                   const Array<int> &dirichlet_marker, real_t penalty_factor)
      : fes_(fes), mesh_(*fes.GetParMesh()), vol_(lambda, mu),
        face_(lambda, mu, 3, -1.0, penalty_factor), dir_(lambda, mu, 3, -1.0, penalty_factor)
   {
      MFEM_VERIFY(mesh_.GetNRanks() == 1, "NitscheLtsSplit: serial (one rank) only in P1");
      dirichlet_marker.Copy(dirichlet_marker_);
      const int ne = mesh_.GetNE();
      elem_faces_.assign(ne, {});
      for (int f = 0; f < mesh_.GetNumFaces(); f++)
      {
         int e1, e2;
         mesh_.GetFaceElements(f, &e1, &e2);
         if (e1 >= 0 && e2 >= 0)
         {
            elem_faces_[e1].push_back({f, 1});
            elem_faces_[e2].push_back({f, 2});
         }
      }
      for (int be = 0; be < mesh_.GetNBE(); be++)
      {
         const int a = mesh_.GetBdrAttribute(be);
         if (a > dirichlet_marker_.Size() || !dirichlet_marker_[a - 1]) { continue; }
         int f, o;
         mesh_.GetBdrElementFace(be, &f, &o);
         int e1, e2;
         mesh_.GetFaceElements(f, &e1, &e2);
         if (e2 >= 0) { continue; }   // BilinearForm skips interior faces too
         elem_faces_[e1].push_back({be, 0});
      }
   }

   /// K^(c) for every cluster c (cluster[e] in [0, nc)), as SparseMatrix on
   /// the local vector DOFs, assembled with skip_zeros = 0 as K is.
   void AssembleClusters(const std::vector<int> &cluster, int nc)
   {
      const int n = fes_.GetVSize();
      K_.clear();
      for (int c = 0; c < nc; c++) { K_.push_back(std::make_unique<SparseMatrix>(n, n)); }
      Array<int> vd, vd2;
      DenseMatrix elmat;
      for (int e = 0; e < mesh_.GetNE(); e++)
      {
         fes_.GetElementVDofs(e, vd);
         vol_.AssembleElementMatrix(*fes_.GetFE(e), *fes_.GetElementTransformation(e), elmat);
         K_[cluster[e]]->AddSubMatrix(vd, vd, elmat, 0);
      }
      for (int f = 0; f < mesh_.GetNumFaces(); f++)
      {
         FaceElementTransformations *tr = mesh_.GetInteriorFaceTransformations(f);
         if (!tr) { continue; }
         fes_.GetElementVDofs(tr->Elem1No, vd);
         fes_.GetElementVDofs(tr->Elem2No, vd2);
         vd.Append(vd2);
         for (int s = 1; s <= 2; s++)
         {
            face_.SetSide(s);
            face_.AssembleFaceMatrix(*fes_.GetFE(tr->Elem1No), *fes_.GetFE(tr->Elem2No),
                                     *tr, elmat);
            K_[cluster[s == 1 ? tr->Elem1No : tr->Elem2No]]->AddSubMatrix(vd, vd, elmat, 0);
         }
      }
      for (int e = 0; e < mesh_.GetNE(); e++)
      {
         for (const auto &fs : elem_faces_[e])
         {
            if (fs.second != 0) { continue; }
            FaceElementTransformations *tr = mesh_.GetBdrFaceTransformations(fs.first);
            fes_.GetElementVDofs(e, vd);
            dir_.AssembleFaceMatrix(*fes_.GetFE(e), *fes_.GetFE(e), *tr, elmat);
            K_[cluster[e]]->AddSubMatrix(vd, vd, elmat, 0);
         }
      }
      for (auto &k : K_) { k->Finalize(0); }
   }

   /// C^(c) for every cluster: the dashpot blocks C_e (element-local
   /// byNODES, NitscheWaveOperator::DampingBlocks) of its elements.
   void AssembleAbsorbingClusters(const std::map<int, DenseMatrix> &blocks,
                                  const std::vector<int> &cluster, int nc)
   {
      const int n = fes_.GetVSize();
      C_.clear();
      for (int c = 0; c < nc; c++) { C_.push_back(std::make_unique<SparseMatrix>(n, n)); }
      Array<int> vd;
      for (const auto &b : blocks)
      {
         fes_.GetElementVDofs(b.first, vd);   // byNODES: component-major, as C_e
         C_[cluster[b.first]]->AddSubMatrix(vd, vd, b.second, 0);
      }
      for (auto &k : C_) { k->Finalize(0); }
   }

   int NumClusters() const { return int(K_.size()); }
   const SparseMatrix &K(int c) const { return *K_[c]; }
   const SparseMatrix &Absorbing(int c) const { return *C_[c]; }

   /// G_e on its patch: vdofs = the vector DOFs of e, then of each interior-
   /// face neighbour in face order; G the dense patch matrix.
   void ElementPatch(int e, Array<int> &vdofs, DenseMatrix &G)
   {
      Array<int> vd, vd2, nb;
      fes_.GetElementVDofs(e, vdofs);
      nb.Append(e);
      for (const auto &fs : elem_faces_[e])
      {
         if (fs.second == 0) { continue; }
         int e1, e2;
         mesh_.GetFaceElements(fs.first, &e1, &e2);
         const int o = (fs.second == 1) ? e2 : e1;
         fes_.GetElementVDofs(o, vd);
         vdofs.Append(vd);
      }
      const int m = vdofs.Size();
      G.SetSize(m);
      G = 0.0;
      // local index of each patch vdof (a patch has at most 5 elements)
      auto local = [&](const Array<int> &v, Array<int> &loc)
      {
         loc.SetSize(v.Size());
         for (int i = 0; i < v.Size(); i++)
         {
            loc[i] = -1;
            for (int j = 0; j < m; j++) { if (vdofs[j] == v[i]) { loc[i] = j; break; } }
            MFEM_VERIFY(loc[i] >= 0, "ElementPatch: vdof outside the patch");
         }
      };
      auto add = [&](const Array<int> &v, const DenseMatrix &a)
      {
         Array<int> loc;
         local(v, loc);
         for (int i = 0; i < v.Size(); i++)
         {
            for (int j = 0; j < v.Size(); j++) { G(loc[i], loc[j]) += a(i, j); }
         }
      };
      DenseMatrix elmat;
      fes_.GetElementVDofs(e, vd);
      vol_.AssembleElementMatrix(*fes_.GetFE(e), *fes_.GetElementTransformation(e), elmat);
      add(vd, elmat);
      for (const auto &fs : elem_faces_[e])
      {
         if (fs.second == 0)
         {
            FaceElementTransformations *tr = mesh_.GetBdrFaceTransformations(fs.first);
            fes_.GetElementVDofs(e, vd);
            dir_.AssembleFaceMatrix(*fes_.GetFE(e), *fes_.GetFE(e), *tr, elmat);
            add(vd, elmat);
            continue;
         }
         FaceElementTransformations *tr = mesh_.GetInteriorFaceTransformations(fs.first);
         fes_.GetElementVDofs(tr->Elem1No, vd);
         fes_.GetElementVDofs(tr->Elem2No, vd2);
         vd.Append(vd2);
         face_.SetSide(fs.second);
         face_.AssembleFaceMatrix(*fes_.GetFE(tr->Elem1No), *fes_.GetFE(tr->Elem2No), *tr,
                                  elmat);
         add(vd, elmat);
      }
   }

private:
   ParFiniteElementSpace &fes_;
   ParMesh &mesh_;
   ElasticityIntegrator vol_;
   DGElasticityIPCombinedIntegrator face_, dir_;
   Array<int> dirichlet_marker_;
   /// per element: (interior face, side 1|2) or (boundary element, 0)
   std::vector<std::vector<std::pair<int, int>>> elem_faces_;
   std::vector<std::unique_ptr<SparseMatrix>> K_, C_;
};

class NitscheMultirateLeapfrog
{
public:
   /// K[c] = K^(c), c = 0 .. Nc-1 (cluster c steps at 2^c dt0); Kv[c] =
   /// K_eta^(c) and C[c] = C^(c), each empty for none; Binv[k] = B_k^{-1} for
   /// k = 0 .. Nc-1 (header), empty for M^{-1} (Minv) at every tick.
   NitscheMultirateLeapfrog(std::vector<const Operator *> K, const NitscheMass &Minv,
                            real_t dt0, std::vector<const Operator *> Kv = {},
                            std::vector<const Operator *> C = {},
                            std::vector<const NitscheMass *> Binv = {})
      : K_(std::move(K)), Kv_(std::move(Kv)), C_(std::move(C)), Minv_(Minv),
        Binv_(std::move(Binv)), dt0_(dt0)
   {
      MFEM_VERIFY(!K_.empty() && K_.size() <= 30, "NitscheMultirateLeapfrog: "
                  << K_.size() << " clusters");
      MFEM_VERIFY((Kv_.empty() || Kv_.size() == K_.size())
                  && (C_.empty() || C_.size() == K_.size()),
                  "NitscheMultirateLeapfrog: " << K_.size() << " clusters, "
                  << Kv_.size() << " damping and " << C_.size() << " absorbing operators");
      MFEM_VERIFY(C_.empty() == Binv_.empty() && (Binv_.empty() || Binv_.size() == K_.size()),
                  "NitscheMultirateLeapfrog: absorbing faces need one B_k^{-1} per cluster "
                  "(" << Binv_.size() << " for " << K_.size() << " clusters), and B_k need them");
      L_ = int(K_.size()) - 1;
      a_.resize(K_.size());
      for (auto &a : a_) { a.SetSize(Minv_.VSize()); }
      r_.SetSize(Minv_.VSize());
   }

   real_t SyncInterval() const { return std::ldexp(dt0_, L_); }
   real_t Dt0() const { return dt0_; }

   /// a_c at the start state (u, v), for every cluster, with M^{-1}.
   void Init(const Vector &u, const Vector &v)
   {
      for (int c = 0; c <= L_; c++) { Force_(c, u, v, Minv_); }
   }

   /// One sync interval from (u, v) at t; v is synchronised at the end.
   void Step(Vector &u, Vector &v, real_t &t) { Step_(u, v, t, false); }

   /// The same interval for the recursion's own state: from v+ = v after the
   /// start kicks at t, to v+ at t + SyncInterval() (the next start kicks
   /// applied).  With damping a_c depends on the velocity before the end
   /// kicks, so (u, v) at sync points does not close the recursion, (u, v+)
   /// does; its map has the stepping's spectrum (P2 eigenvalue scan).
   void StepKicked(Vector &u, Vector &v, real_t &t) { Step_(u, v, t, true); }

private:
   void Step_(Vector &u, Vector &v, real_t &t, bool kicked)
   {
      const int nt = 1 << L_;
      const real_t t0 = t;
      for (int m = 0; m < nt; m++)
      {
         for (int c = 0; c <= L_; c++)
         {
            if (m % (1 << c) == 0 && !(kicked && m == 0))
            {
               v.Add(0.5 * std::ldexp(dt0_, c), a_[c]);
            }
         }
         u.Add(dt0_, v);
         t = t0 + (m + 1) * dt0_;
         int top = 0;   // clusters 0 .. top end a step at t
         while (top < L_ && (m + 1) % (1 << (top + 1)) == 0) { top++; }
         const NitscheMass &mass = Binv_.empty() ? Minv_ : *Binv_[top];
         for (int c = 0; c <= top; c++) { Force_(c, u, v, mass); }
         for (int c = 0; c <= top; c++) { v.Add(0.5 * std::ldexp(dt0_, c), a_[c]); }
      }
      if (kicked)
      {
         for (int c = 0; c <= L_; c++) { v.Add(0.5 * std::ldexp(dt0_, c), a_[c]); }
      }
   }

   void Force_(int c, const Vector &u, const Vector &v, const NitscheMass &mass)
   {
      K_[c]->Mult(u, r_);
      if (!Kv_.empty()) { Kv_[c]->AddMult(v, r_); }
      if (!C_.empty()) { C_[c]->AddMult(v, r_); }
      r_.Neg();
      mass.Apply(r_, a_[c]);
   }

   std::vector<const Operator *> K_, Kv_, C_;
   const NitscheMass &Minv_;
   std::vector<const NitscheMass *> Binv_;
   real_t dt0_;
   int L_ = 0;
   std::vector<Vector> a_;
   Vector r_;
};

#endif // MFEM_USE_MPI

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_NITSCHE_LTS_STEPPER_HPP

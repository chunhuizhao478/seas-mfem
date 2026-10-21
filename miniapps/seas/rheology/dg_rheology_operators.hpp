// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dg_rheology_operators.hpp — the three operators that couple the rheology
// field to the DG elasticity solver (rheology plan rev 1, sections 4.2 and
// 6.5, Appendix A).
//
//   S  strain at the rheology nodes:   eps_i = sym(grad u_h)(x_i)
//   P  force of the stress correction: (P dsig)_w = - sum_K int_K eps(w) : dsig_h
//                                              + sum_F int_F {dsig_h n} . [w]
//   Q  its traction on the fault:      {dsig_h n} at the fault quadrature
//                                      points, projected by the operator
//                                      like the prestress (FaceAverage)
//
// The faces of P are those that carry the SIPG consistency term in K:
// interior faces (fault faces included), shared faces and Dirichlet faces
// (one-sided there); free and absorbing faces have none.  n is the unit
// normal outward from element 1 and [w] = w_1 - w_2, as in
// DGElasticityIPCombinedIntegrator::AssembleFaceMatrix (consistency term
// -{sigma(u) n} . [w], blocks (0,0) c0 = -1/2 and (1,0) c0 = +1/2), so that
// P (C_0 : eps(v)) = -K v for a continuous v (test O4) and a uniform dsig
// gives only the traction on the free and absorbing faces (Proposition A.1,
// test O2).
//
// Tables are per reference element (gradients at the rheology nodes, the
// volume rule of order 2p - 2) and per element (inverse Jacobian, |J|: the
// elements must be affine); per face the shape values at the 2p + 1 face
// rule, the order of K's face terms.  P runs on the active elements and on
// the faces that touch one, so a run in which nothing yields does no P work.
//
// Scope (Chunhui Zhao, 2026-10-20).

#ifndef MFEM_SEAS_DG_RHEOLOGY_OPERATORS_HPP
#define MFEM_SEAS_DG_RHEOLOGY_OPERATORS_HPP

#include "mfem.hpp"
#include "rheology_field.hpp"

#include <set>
#include <vector>

namespace mfem
{
namespace seas
{

#ifdef MFEM_USE_MPI

class DGRheologyOperators
{
public:
   /// ufes: the displacement space (DG, 3 components, byNODES, as
   /// ElasticityDomainOperator::SetupFESpace); dirichlet: the boundary
   /// attributes with the SIPG boundary integrator.  Collective (shared
   /// face transformations).
   DGRheologyOperators(ParFiniteElementSpace &ufes, const RheologyField &rf,
                       const std::set<int> &dirichlet)
   {
      ParMesh &pm = *ufes.GetParMesh();
      MFEM_VERIFY(ufes.GetVDim() == 3 && ufes.GetOrdering() == Ordering::byNODES,
                  "DGRheologyOperators: the displacement space must have 3 components, "
                  "byNODES");
      ne_ = pm.GetNE();
      MFEM_VERIFY(ne_ == rf.NumElements(), "DGRheologyOperators: field and space differ");
      MFEM_VERIFY(ufes.GetTrueVSize() == ufes.GetVSize(),
                  "DGRheologyOperators: the displacement space must be discontinuous");
      const FiniteElement &fer = rf.NodeElement();
      nr_ = rf.NodesPerElement();
      const FiniteElement *feu = (ne_ > 0) ? ufes.GetFE(0)
                                 : ufes.FEColl()->FiniteElementForGeometry(Geometry::TETRAHEDRON);
      nd_ = feu->GetDof();
      p_ = feu->GetOrder();
      MFEM_VERIFY(fer.GetOrder() == p_ - 1, "DGRheologyOperators: rheology space has order "
                  << fer.GetOrder() << ", expected " << p_ - 1);

      // reference gradients at the rheology nodes
      DenseMatrix ds(nd_, 3);
      Gn_.resize(static_cast<std::size_t>(nr_) * nd_ * 3);
      for (int i = 0; i < nr_; i++)
      {
         feu->CalcDShape(fer.GetNodes().IntPoint(i), ds);
         for (int j = 0; j < nd_; j++)
         {
            for (int k = 0; k < 3; k++) { Gn_[(static_cast<std::size_t>(i) * nd_ + j) * 3 + k] = ds(j, k); }
         }
      }
      // volume rule 2p - 2: grad(phi) psi has degree 2p - 2 on an affine tet
      const IntegrationRule &irv = IntRules.Get(Geometry::TETRAHEDRON, 2 * p_ - 2);
      nqv_ = irv.GetNPoints();
      wv_.resize(nqv_);
      Pv_.resize(static_cast<std::size_t>(nqv_) * nr_);
      Gv_.resize(static_cast<std::size_t>(nqv_) * nd_ * 3);
      Vector sh(nr_);
      for (int q = 0; q < nqv_; q++)
      {
         const IntegrationPoint &ip = irv.IntPoint(q);
         wv_[q] = ip.weight;
         fer.CalcShape(ip, sh);
         for (int i = 0; i < nr_; i++) { Pv_[static_cast<std::size_t>(q) * nr_ + i] = sh(i); }
         feu->CalcDShape(ip, ds);
         for (int j = 0; j < nd_; j++)
         {
            for (int k = 0; k < 3; k++) { Gv_[(static_cast<std::size_t>(q) * nd_ + j) * 3 + k] = ds(j, k); }
         }
      }
      // per element: inverse Jacobian, |J|, displacement vdofs
      invJ_.resize(static_cast<std::size_t>(ne_) * 9);
      detJ_.resize(ne_);
      vdofs_.resize(static_cast<std::size_t>(ne_) * 3 * nd_);
      Array<int> vd;
      for (int e = 0; e < ne_; e++)
      {
         ElementTransformation *T = pm.GetElementTransformation(e);
         MFEM_VERIFY(T->OrderJ() == 0, "DGRheologyOperators: element " << e
                     << " is not affine; the rheology nodes need straight tetrahedra");
         T->SetIntPoint(&Geometries.GetCenter(Geometry::TETRAHEDRON));
         const DenseMatrix &Ji = T->InverseJacobian();
         for (int a = 0; a < 3; a++)
         {
            for (int b = 0; b < 3; b++) { invJ_[static_cast<std::size_t>(e) * 9 + 3 * a + b] = Ji(a, b); }
         }
         detJ_[e] = T->Weight();
         ufes.GetElementVDofs(e, vd);
         MFEM_VERIFY(vd.Size() == 3 * nd_, "DGRheologyOperators: element " << e << " has "
                     << vd.Size() << " vdofs");
         for (int a = 0; a < 3 * nd_; a++)
         {
            MFEM_VERIFY(vd[a] >= 0, "DGRheologyOperators: negative vdof");
            vdofs_[static_cast<std::size_t>(e) * 3 * nd_ + a] = vd[a];
         }
      }
      // faces with a consistency term
      for (int f = 0; f < pm.GetNumFaces(); f++)
      {
         FaceElementTransformations *Tr = pm.GetInteriorFaceTransformations(f);
         if (Tr) { AddFace_(*Tr, Face::Interior, *feu, fer); }
      }
      for (int sf = 0; sf < pm.GetNSharedFaces(); sf++)
      {
         FaceElementTransformations *Tr = pm.GetSharedFaceTransformations(sf);
         MFEM_VERIFY(Tr, "DGRheologyOperators: no transformation for shared face " << sf);
         AddFace_(*Tr, Face::Shared, *feu, fer);
      }
      for (int be = 0; be < pm.GetNBE(); be++)
      {
         if (!dirichlet.count(pm.GetBdrAttribute(be))) { continue; }
         FaceElementTransformations *Tr = pm.GetBdrFaceTransformations(be);
         if (Tr && Tr->Elem2No < 0) { AddFace_(*Tr, Face::Boundary, *feu, fer); }
      }
   }

   /// Faces of one kind: 0 interior (fault faces included), 1 shared, 2
   /// Dirichlet boundary.
   int NumFaces(int kind) const
   {
      int n = 0;
      for (const FaceRec &r : faces_) { n += (static_cast<int>(r.kind) == kind) ? 1 : 0; }
      return n;
   }

   /// Faces of one kind that touch an active element (the faces AddFaces
   /// works on), after RheologyCoupler::Refresh.
   int NumActiveFaces(const RheologyField &rf, int kind) const
   {
      int n = 0;
      for (const FaceRec &r : faces_)
      {
         if (static_cast<int>(r.kind) != kind) { continue; }
         n += (rf.Active(r.e1) || (r.kind == Face::Interior && rf.Active(r.e2))
               || (r.kind == Face::Shared && rf.NeighbourActive(r.e2))) ? 1 : 0;
      }
      return n;
   }

   /// S: the strain of U (local displacement vector) at the nr nodes of
   /// element e, eps[6 i + c] in the [xx, yy, zz, xy, yz, xz] order.
   void Strain(int e, const Vector &U, real_t *eps) const
   {
      const int *vd = &vdofs_[static_cast<std::size_t>(e) * 3 * nd_];
      const real_t *Ji = &invJ_[static_cast<std::size_t>(e) * 9];
      for (int i = 0; i < nr_; i++)
      {
         // H_ref(c, k) = sum_j u_c,j dphi_j/dxi_k ; H = H_ref J^{-1}
         real_t Hr[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
         const real_t *G = &Gn_[static_cast<std::size_t>(i) * nd_ * 3];
         for (int c = 0; c < 3; c++)
         {
            for (int j = 0; j < nd_; j++)
            {
               const real_t u = U(vd[c * nd_ + j]);
               Hr[3 * c] += u * G[3 * j];
               Hr[3 * c + 1] += u * G[3 * j + 1];
               Hr[3 * c + 2] += u * G[3 * j + 2];
            }
         }
         real_t H[9];
         for (int c = 0; c < 3; c++)
         {
            for (int d = 0; d < 3; d++)
            {
               H[3 * c + d] = Hr[3 * c] * Ji[d] + Hr[3 * c + 1] * Ji[3 + d] + Hr[3 * c + 2] * Ji[6 + d];
            }
         }
         sym6::FromRowMajor(H, eps + 6 * i);
      }
   }

   /// F += the volume part of P dsig on element e:
   /// F_(c, j) -= int_e sum_d dphi_j/dx_d dsig_cd.
   void AddVolume(int e, const RheologyField &rf, Vector &F) const
   {
      const int *vd = &vdofs_[static_cast<std::size_t>(e) * 3 * nd_];
      const real_t *Ji = &invJ_[static_cast<std::size_t>(e) * 9];
      for (int q = 0; q < nqv_; q++)
      {
         // dsig_h at the point, row-major
         real_t s6[6] = {0, 0, 0, 0, 0, 0};
         for (int i = 0; i < nr_; i++)
         {
            const real_t w = Pv_[static_cast<std::size_t>(q) * nr_ + i];
            const real_t *ds = rf.Dsig(e, i);
            for (int c = 0; c < 6; c++) { s6[c] += w * ds[c]; }
         }
         real_t S[9];
         sym6::ToRowMajor(s6, S);
         // A(c, k) = w |J| sum_d S(c, d) Ji(k, d)
         const real_t wq = wv_[q] * detJ_[e];
         real_t A[9];
         for (int c = 0; c < 3; c++)
         {
            for (int k = 0; k < 3; k++)
            {
               A[3 * c + k] = wq * (S[3 * c] * Ji[3 * k] + S[3 * c + 1] * Ji[3 * k + 1]
                                    + S[3 * c + 2] * Ji[3 * k + 2]);
            }
         }
         const real_t *G = &Gv_[static_cast<std::size_t>(q) * nd_ * 3];
         for (int c = 0; c < 3; c++)
         {
            for (int j = 0; j < nd_; j++)
            {
               F(vd[c * nd_ + j]) -= G[3 * j] * A[3 * c] + G[3 * j + 1] * A[3 * c + 1]
                                     + G[3 * j + 2] * A[3 * c + 2];
            }
         }
      }
   }

   /// F += the face part of P dsig on every face that touches an active
   /// element: element 1 gets +int phi_1 {dsig n}, a local element 2
   /// -int phi_2 {dsig n} (the rank of a shared face's element 2 adds its own
   /// side with its own outward normal).  After rf.ExchangeNeighbours().
   void AddFaces(const RheologyField &rf, Vector &F) const
   {
      for (const FaceRec &r : faces_)
      {
         const bool on = rf.Active(r.e1)
                         || (r.kind == Face::Interior && rf.Active(r.e2))
                         || (r.kind == Face::Shared && rf.NeighbourActive(r.e2));
         if (!on) { continue; }
         const int *vd1 = &vdofs_[static_cast<std::size_t>(r.e1) * 3 * nd_];
         const int *vd2 = (r.kind == Face::Interior)
                          ? &vdofs_[static_cast<std::size_t>(r.e2) * 3 * nd_] : nullptr;
         const real_t *d = &fdata_[r.off];
         const std::size_t stride = FaceStride_(r.kind);
         for (int q = 0; q < r.nq; q++, d += stride)
         {
            const real_t *wn = d;
            const real_t *phi1 = d + 3;
            const real_t *phi2 = phi1 + nd_;
            const real_t *psi1 = phi1 + nd_ * (r.kind == Face::Interior ? 2 : 1);
            const real_t *psi2 = psi1 + nr_;
            real_t s6[6];
            SideStress_(rf, r, psi1, psi2, s6);
            real_t S[9];
            sym6::ToRowMajor(s6, S);
            real_t t[3];
            for (int c = 0; c < 3; c++) { t[c] = S[3 * c] * wn[0] + S[3 * c + 1] * wn[1] + S[3 * c + 2] * wn[2]; }
            for (int c = 0; c < 3; c++)
            {
               for (int j = 0; j < nd_; j++) { F(vd1[c * nd_ + j]) += phi1[j] * t[c]; }
            }
            if (vd2)
            {
               for (int c = 0; c < 3; c++)
               {
                  for (int j = 0; j < nd_; j++) { F(vd2[c * nd_ + j]) -= phi2[j] * t[c]; }
               }
            }
         }
      }
   }

   /// Q at one point of a fault face: s = {dsig_h}, row-major, with the
   /// integration points of Tr already set (as ProjectStress_ does).
   /// Element 2 of a shared face is face-neighbour element Elem2No - NE.
   void FaceAverage(const RheologyField &rf, FaceElementTransformations &Tr, real_t s[9]) const
   {
      const FiniteElement &fer = rf.NodeElement();
      Vector psi(nr_);
      real_t s6[6] = {0, 0, 0, 0, 0, 0};
      MFEM_VERIFY(Tr.Elem2No >= 0, "DGRheologyOperators::FaceAverage: not an interior face");
      fer.CalcShape(Tr.GetElement1IntPoint(), psi);
      for (int i = 0; i < nr_; i++)
      {
         const real_t *ds = rf.Dsig(Tr.Elem1No, i);
         for (int c = 0; c < 6; c++) { s6[c] += 0.5 * psi(i) * ds[c]; }
      }
      fer.CalcShape(Tr.GetElement2IntPoint(), psi);
      const int e2 = Tr.Elem2No;
      for (int i = 0; i < nr_; i++)
      {
         for (int c = 0; c < 6; c++)
         {
            const real_t v = (e2 < ne_) ? rf.Dsig(e2, i)[c] : rf.NeighbourDsig(e2 - ne_, i, c);
            s6[c] += 0.5 * psi(i) * v;
         }
      }
      sym6::ToRowMajor(s6, s);
   }

private:
   enum class Face { Interior = 0, Shared = 1, Boundary = 2 };

   struct FaceRec
   {
      Face kind;
      int e1, e2;   ///< e2: local element, face-neighbour element, or -1
      int nq;
      std::size_t off;
   };

   /// Per point: w n (3), phi_1 (nd), [phi_2 (nd), interior], [psi_1, psi_2
   /// (nr each), not boundary] -- boundary: psi_1 only.
   std::size_t FaceStride_(Face k) const
   {
      switch (k)
      {
         case Face::Interior: return 3 + 2 * nd_ + 2 * nr_;
         case Face::Shared: return 3 + nd_ + 2 * nr_;
         case Face::Boundary: return 3 + nd_ + nr_;
      }
      return 0;
   }

   /// {dsig_h} at a face point (one-sided on a boundary face), 6 components.
   void SideStress_(const RheologyField &rf, const FaceRec &r, const real_t *psi1,
                    const real_t *psi2, real_t s6[6]) const
   {
      const real_t h = (r.kind == Face::Boundary) ? 1.0 : 0.5;
      for (int c = 0; c < 6; c++) { s6[c] = 0.0; }
      for (int i = 0; i < nr_; i++)
      {
         const real_t *ds = rf.Dsig(r.e1, i);
         for (int c = 0; c < 6; c++) { s6[c] += h * psi1[i] * ds[c]; }
      }
      if (r.kind == Face::Interior)
      {
         for (int i = 0; i < nr_; i++)
         {
            const real_t *ds = rf.Dsig(r.e2, i);
            for (int c = 0; c < 6; c++) { s6[c] += h * psi2[i] * ds[c]; }
         }
      }
      else if (r.kind == Face::Shared)
      {
         for (int i = 0; i < nr_; i++)
         {
            for (int c = 0; c < 6; c++) { s6[c] += h * psi2[i] * rf.NeighbourDsig(r.e2, i, c); }
         }
      }
   }

   void AddFace_(FaceElementTransformations &Tr, Face kind, const FiniteElement &feu,
                 const FiniteElement &fer)
   {
      FaceRec r;
      r.kind = kind;
      r.e1 = Tr.Elem1No;
      r.e2 = (kind == Face::Interior) ? Tr.Elem2No
             : (kind == Face::Shared) ? Tr.Elem2No - ne_ : -1;
      MFEM_VERIFY(r.e1 >= 0 && r.e1 < ne_, "DGRheologyOperators: face element 1 = " << r.e1);
      MFEM_VERIFY(kind != Face::Interior || (r.e2 >= 0 && r.e2 < ne_),
                  "DGRheologyOperators: interior face element 2 = " << r.e2);
      MFEM_VERIFY(kind != Face::Shared || r.e2 >= 0,
                  "DGRheologyOperators: shared face neighbour " << r.e2);
      const IntegrationRule &ir = IntRules.Get(Tr.GetGeometryType(), 2 * p_ + 1);
      r.nq = ir.GetNPoints();
      r.off = fdata_.size();
      const std::size_t stride = FaceStride_(kind);
      fdata_.resize(r.off + stride * r.nq);
      Vector nor(3), phi(nd_), psi(nr_);
      for (int q = 0; q < r.nq; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Tr.SetAllIntPoints(&ip);
         real_t *d = &fdata_[r.off + stride * q];
         CalcOrtho(Tr.Jacobian(), nor);   // outward from element 1, |nor| = face |J|
         for (int c = 0; c < 3; c++) { d[c] = ip.weight * nor(c); }
         real_t *o = d + 3;
         feu.CalcShape(Tr.GetElement1IntPoint(), phi);
         for (int j = 0; j < nd_; j++) { o[j] = phi(j); }
         o += nd_;
         if (kind == Face::Interior)
         {
            feu.CalcShape(Tr.GetElement2IntPoint(), phi);
            for (int j = 0; j < nd_; j++) { o[j] = phi(j); }
            o += nd_;
         }
         fer.CalcShape(Tr.GetElement1IntPoint(), psi);
         for (int i = 0; i < nr_; i++) { o[i] = psi(i); }
         o += nr_;
         if (kind != Face::Boundary)
         {
            fer.CalcShape(Tr.GetElement2IntPoint(), psi);
            for (int i = 0; i < nr_; i++) { o[i] = psi(i); }
         }
      }
      faces_.push_back(r);
   }

   int ne_ = 0, nr_ = 0, nd_ = 0, p_ = 0, nqv_ = 0;
   std::vector<real_t> Gn_;      ///< [i][j][k] dphi_j/dxi_k at rheology node i
   std::vector<real_t> wv_;      ///< volume weights
   std::vector<real_t> Pv_;      ///< [q][i] psi_i at volume point q
   std::vector<real_t> Gv_;      ///< [q][j][k] dphi_j/dxi_k at volume point q
   std::vector<real_t> invJ_;    ///< [e] row-major J^{-1}
   std::vector<real_t> detJ_;
   std::vector<int> vdofs_;      ///< [e][c nd + j]
   std::vector<FaceRec> faces_;
   std::vector<real_t> fdata_;
};

#endif // MFEM_USE_MPI

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_DG_RHEOLOGY_OPERATORS_HPP

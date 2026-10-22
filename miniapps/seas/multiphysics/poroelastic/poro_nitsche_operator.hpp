// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// poro_nitsche_operator.hpp — poroelastic option, step P2.2 (plan rev 3).
//
// PoroNitscheOperator is the bulk (u, w) operator of derivation rev 1.2,
// eq. (8.1): the stiffness
//
//        K = [ K_E^u  Z_uw ]      K_E^u = today's elastic SIPG stiffness built
//            [ Z_wu   Z_ww ]      from the undrained sample (lambda_u, mu, rho)
//
// and the right-hand side of the acceleration solve (9.2),
//
//   F_u = Force(u, v*_u, d) + F_V,u(t)          - Z_uw w
//   F_w = F_D,w(t) + F_V,w(t) + F_init - Z_wu u - Z_ww w - b M_e v*_w
//
// Force(u, v*_u, d) is the elastic operator's own force (hook H1): its
// Dirichlet load, - K_E^u u, the slip load G^u d, the dashpot - C v*_u and the
// Kelvin-Voigt terms, all with lambda_u (option F2, derivation Proposition
// 6.2).  So with w = 0 the u rows are the elastic force bit for bit (U2-4).
// The operator owns that NitscheWaveOperator; the fault advance, the time
// step and the mass are the callers' (steps P2.4, P2.5).
//
// The Z blocks are DGPoroDivIntegrator (dg_poro_div_integrator.hpp), B_c of
// (6.3), assembled once into HypreParMatrix:
//   Z_wu   c = alpha M, volume and interior faces (the permeable fault faces
//          included, derivation section 7.1), and the boundary rows of
//          section 6.6: kTrialNormal on Dirichlet-u faces, kTestNormal on
//          no-flow faces, kBothNormal on faces that are both;
//   Z_uw   = Z_wu^T;
//   Z_ww   c = M with the fluid penalty kappa_w of (6.5), volume and interior
//          faces, and kBothNormal with the one-sided penalty on no-flow faces.
// Free-surface, absorbing and Dirichlet faces are drained (p = 0): no fluid
// terms (section 6.6).
//
// Loads, each assembled with MFEM's ParLinearForm loops:
//   F_D,w  the Dirichlet load of the kTrialNormal rows, -int alpha M (div q)
//          (g.n) on Dirichlet-u faces, with the elastic operator's own
//          DirichletFunc per attribute, at t = GetTime() in every call;
//   F_V    (6.6) for an injected volume V(x, t) (SetInjectedVolume, or
//          SetSource with a gaussian PoroSource): c = alpha M on the solid
//          rows, with the Dirichlet faces; c = M on the fluid rows, with the
//          no-flow faces; a well PoroSource adds its point load instead
//          (poro_source.hpp), at the source's own time, which must be the
//          operator's;
//   F_init (6.6) for a starting pressure (SetInitialPressure), c = 1 on the
//          fluid rows, with the no-flow faces; assembled once.
// The drag b M_e v*_w is element-local: M_e is the scalar DG mass of
// NitscheMass (same rule IntRules.Get(geom, 2p), same weights), so the
// PoroMass of step P2.4 inverts exactly this matrix.
//
// Layout: x = [u; w], each in the elastic operator's vector DG layout
// (byNODES, size n = Elastic().Height()).
//
// Refused at construction: a Dirichlet or no-flow attribute on a face that is
// not a true boundary face (the elastic operator's interior "skeleton"
// Dirichlet faces of the quasi-dynamic BP5 loading have no fluid rows here),
// a no-flow attribute on the fault, and kappa_w <= 0.
//
// Scope (Chunhui Zhao, 2026-10-13): tetrahedral meshes only.

#ifndef MFEM_SEAS_PORO_NITSCHE_OPERATOR_HPP
#define MFEM_SEAS_PORO_NITSCHE_OPERATOR_HPP

#include "mfem.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "dg_poro_div_integrator.hpp"
#include "poro_material.hpp"
#include "poro_source.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <set>
#include <vector>

namespace mfem
{
namespace seas
{
namespace poro
{

#ifdef MFEM_USE_MPI

/// Fluid walls that no deck sets yet (derivation section 6.6): w.n = 0 on
/// these boundary attributes, for the quasi-static tests V-Q1 to V-Q3 of
/// phase P7.  Every other boundary face is drained (p = 0).
struct PoroBoundaryConfig
{
   std::set<int> no_flow_attrs;
};

class PoroNitscheOperator
{
public:
   /// `sample` is the Biot sample of `mesh` (kept, not owned); the elastic
   /// operator is built from its undrained material.  kappa_w is the fluid
   /// penalty of (6.5) ([multiphysics.poroelastic] fluid_penalty).
   PoroNitscheOperator(ParMesh &mesh, int order, const PoroMaterialSample &sample,
                       const BoundaryConfig &bc, real_t kappa_w,
                       const PoroBoundaryConfig &pbc = {},
                       const DomainConfig &config = {},
                       const NitscheFaultSetup *fault = nullptr)
      : mesh_(mesh), sample_(sample), bc_(bc), pbc_(pbc), kappa_w_(kappa_w)
   {
      MFEM_VERIFY(kappa_w > 0.0, "PoroNitscheOperator: the fluid penalty kappa_w ("
                  << kappa_w << ") must be > 0.");
      MFEM_VERIFY(sample.AlphaGF().ParFESpace()->GetParMesh() == &mesh,
                  "PoroNitscheOperator: the Biot sample belongs to another mesh.");
      CheckBoundary_();

      elastic_ = std::make_unique<NitscheWaveOperator>(
                    mesh, order, sample.UndrainedMaterial(), bc, config, fault);
      n_ = elastic_->Height();
      ParFiniteElementSpace &fes = Space_();

      // alpha M per element, as one value (the product of the sample's two).
      l2_fec_ = std::make_unique<L2_FECollection>(0, mesh.Dimension());
      l2_fes_ = std::make_unique<ParFiniteElementSpace>(&mesh, l2_fec_.get());
      alphaM_gf_ = std::make_unique<ParGridFunction>(l2_fes_.get());
      Array<int> dofs;
      for (int e = 0; e < mesh.GetNE(); e++)
      {
         l2_fes_->GetElementDofs(e, dofs);
         const BiotConstants &c = sample.Constants(e);
         (*alphaM_gf_)(dofs[0]) = c.alpha * c.M;
      }
      alphaM_gf_->ExchangeFaceNbrData();
      c_aM_ = std::make_unique<GridFunctionCoefficient>(alphaM_gf_.get());
      c_M_ = std::make_unique<GridFunctionCoefficient>(&sample.BiotModulusGF());
      c_one_ = std::make_unique<ConstantCoefficient>(1.0);

      // Boundary markers.
      const int nattr = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
      auto marker = [nattr](auto pick)
      {
         Array<int> m(nattr);
         m = 0;
         for (int a = 1; a <= nattr; a++) { m[a - 1] = pick(a) ? 1 : 0; }
         return m;
      };
      const auto &D = bc.dirichlet_attrs;
      const auto &NF = pbc.no_flow_attrs;
      dir_marker_ = marker([&](int a) { return D.count(a) > 0; });
      nf_marker_ = marker([&](int a) { return NF.count(a) > 0; });
      Array<int> dir_only = marker([&](int a) { return D.count(a) && !NF.count(a); });
      Array<int> nf_only = marker([&](int a) { return NF.count(a) && !D.count(a); });
      Array<int> both = marker([&](int a) { return D.count(a) && NF.count(a); });
      auto any = [](const Array<int> &m) { for (int v : m) { if (v) { return true; } } return false; };
      has_dir_ = any(dir_marker_);
      has_nf_ = any(nf_marker_);

      // Z_wu (c = alpha M) and Z_uw = Z_wu^T.
      {
         ParBilinearForm z(&fes);
         z.AddDomainIntegrator(new DGPoroDivIntegrator(*c_aM_));
         z.AddInteriorFaceIntegrator(new DGPoroDivIntegrator(*c_aM_));
         if (any(dir_only))
         {
            z.AddBdrFaceIntegrator(new DGPoroDivIntegrator(
               *c_aM_, 0.0, DGPoroDivIntegrator::kTrialNormal), dir_only);
         }
         if (any(nf_only))
         {
            z.AddBdrFaceIntegrator(new DGPoroDivIntegrator(
               *c_aM_, 0.0, DGPoroDivIntegrator::kTestNormal), nf_only);
         }
         if (any(both))
         {
            z.AddBdrFaceIntegrator(new DGPoroDivIntegrator(
               *c_aM_, 0.0, DGPoroDivIntegrator::kBothNormal), both);
         }
         z.Assemble(0);
         z.Finalize(0);
         Zwu_.reset(z.ParallelAssemble());
      }
      Zuw_.reset(Zwu_->Transpose());

      // Z_ww (c = M, with the fluid penalty).
      {
         ParBilinearForm z(&fes);
         z.AddDomainIntegrator(new DGPoroDivIntegrator(*c_M_, kappa_w));
         z.AddInteriorFaceIntegrator(new DGPoroDivIntegrator(*c_M_, kappa_w));
         if (has_nf_)
         {
            z.AddBdrFaceIntegrator(new DGPoroDivIntegrator(
               *c_M_, kappa_w, DGPoroDivIntegrator::kBothNormal), nf_marker_);
         }
         z.Assemble(0);
         z.Finalize(0);
         Zww_.reset(z.ParallelAssemble());
      }

      // The source integrators' coefficients (their boundary flags are not
      // used by AssembleSourceRHS / AssembleDirichletRHS).
      B_aM_ = std::make_unique<DGPoroDivIntegrator>(*c_aM_);
      B_M_ = std::make_unique<DGPoroDivIntegrator>(*c_M_);
      B_one_ = std::make_unique<DGPoroDivIntegrator>(*c_one_);

      // F_D,w: one VectorFunctionCoefficient per Dirichlet attribute, with
      // the function the elastic operator uses there.
      if (has_dir_)
      {
         fdw_ = std::make_unique<ParLinearForm>(&fes);
         for (int a : D)
         {
            auto it = bc.dirichlet_funcs.find(a);
            const DirichletFunc &f = (it != bc.dirichlet_funcs.end())
                                     ? it->second : bc.default_dirichlet_func;
            MFEM_VERIFY(f, "PoroNitscheOperator: no DirichletFunc for boundary "
                        "attribute " << a << "; set default_dirichlet_func or a "
                        "per-attribute function in BoundaryConfig.");
            g_.push_back(std::make_unique<VectorFunctionCoefficient>(3, f));
            g_marker_.push_back(std::make_unique<Array<int>>(
                                   marker([a](int b) { return b == a; })));
            fdw_->AddBdrFaceIntegrator(new DGPoroDivDirichletIntegrator(*B_aM_, *g_.back()),
                                       *g_marker_.back());
         }
      }

      // Drag: b_e and the scalar mass M_e of NitscheMass.
      const NitscheMass &mass = elastic_->Mass();
      b_.SetSize(mesh.GetNE());
      me_.resize(mesh.GetNE());
      for (int e = 0; e < mesh.GetNE(); e++)
      {
         b_(e) = sample.Constants(e).b;
         const FiniteElement *fe = fes.GetFE(e);
         ElementTransformation *T = mesh.GetElementTransformation(e);
         const int nd = fe->GetDof();
         MFEM_VERIFY(nd == mass.ElementScalarDofs(e).Size(),
                     "PoroNitscheOperator: element " << e << " has " << nd
                     << " DOFs in the stiffness space and "
                     << mass.ElementScalarDofs(e).Size() << " in the mass.");
         DenseMatrix &M = me_[e];
         M.SetSize(nd);
         M = 0.0;
         Vector shape(nd);
         const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(), 2 * fe->GetOrder());
         for (int j = 0; j < ir.GetNPoints(); ++j)
         {
            const IntegrationPoint &ip = ir.IntPoint(j);
            T->SetIntPoint(&ip);
            fe->CalcShape(ip, shape);
            const real_t w = ip.weight * T->Weight();
            for (int k = 0; k < nd; ++k)
            {
               for (int l = 0; l < nd; ++l) { M(k, l) += w * shape(k) * shape(l); }
            }
         }
      }
      tmp_.SetSize(n_);
   }

   /// The elastic operator of the solid block (undrained sample).
   NitscheWaveOperator &Elastic() { return *elastic_; }
   const NitscheWaveOperator &Elastic() const { return *elastic_; }
   const PoroMaterialSample &Sample() const { return sample_; }
   ParMesh &GetMesh() const { return mesh_; }
   /// Size of one field (u or w); x = [u; w] has 2 n.
   int NumU() const { return n_; }
   int Height() const { return 2 * n_; }
   real_t FluidPenalty() const { return kappa_w_; }

   const HypreParMatrix &Zwu() const { return *Zwu_; }
   const HypreParMatrix &Zuw() const { return *Zuw_; }
   const HypreParMatrix &Zww() const { return *Zww_; }
   /// The scalar DG mass M_e of local element e (as NitscheMass's) and b_e.
   const DenseMatrix &ElementMass(int e) const { return me_[e]; }
   real_t Drag(int e) const { return b_(e); }

   void SetTime(real_t t) { elastic_->SetTime(t); }
   real_t GetTime() const { return elastic_->GetTime(); }

   /// The injected volume V(x, t) of derivation (4.3), or null for none.
   /// Kept, not owned; F_V is assembled at GetTime() in every Force call.
   void SetInjectedVolume(Coefficient *V)
   {
      V_ = V;
      fvu_.reset();
      fvw_.reset();
      if (!V) { return; }
      ParFiniteElementSpace &fes = Space_();
      fvu_ = std::make_unique<ParLinearForm>(&fes);
      fvu_->AddDomainIntegrator(new DGPoroDivSourceIntegrator(*B_aM_, *V));
      fvu_->AddInteriorFaceIntegrator(new DGPoroDivSourceIntegrator(*B_aM_, *V));
      if (has_dir_)
      {
         fvu_->AddBdrFaceIntegrator(new DGPoroDivSourceIntegrator(*B_aM_, *V), dir_marker_);
      }
      fvw_ = std::make_unique<ParLinearForm>(&fes);
      fvw_->AddDomainIntegrator(new DGPoroDivSourceIntegrator(*B_M_, *V));
      fvw_->AddInteriorFaceIntegrator(new DGPoroDivSourceIntegrator(*B_M_, *V));
      if (has_nf_)
      {
         fvw_->AddBdrFaceIntegrator(new DGPoroDivSourceIntegrator(*B_M_, *V), nf_marker_);
      }
   }

   /// The fluid source of step P2.3, or null for none: a gaussian through
   /// SetInjectedVolume, a well through its point load.  Kept, not owned; the
   /// caller advances it, and every Force call checks that its Time() is
   /// GetTime() (to 1e-12 relative).
   void SetSource(PoroSource *src)
   {
      src_ = src;
      SetInjectedVolume(src && !src->IsWell() ? &src->InjectedVolume() : nullptr);
   }

   /// The starting pressure p_init(x) (derivation section 4.2): F_init of
   /// (6.6), fluid rows, assembled once here.
   void SetInitialPressure(Coefficient &p_init)
   {
      ParLinearForm f(&Space_());
      f.AddDomainIntegrator(new DGPoroDivSourceIntegrator(*B_one_, p_init));
      f.AddInteriorFaceIntegrator(new DGPoroDivSourceIntegrator(*B_one_, p_init));
      if (has_nf_)
      {
         f.AddBdrFaceIntegrator(new DGPoroDivSourceIntegrator(*B_one_, p_init), nf_marker_);
      }
      f.Assemble();
      finit_ = f;
      has_init_ = true;
   }

   /// Kx = [K_E^u u + Z_uw w; Z_wu u + Z_ww w]: the stiffness, no loads.
   void ApplyStiffness(const Vector &x, Vector &Kx) const
   {
      MFEM_VERIFY(x.Size() == 2 * n_, "PoroNitscheOperator::ApplyStiffness: size "
                  << x.Size() << " != 2 n = " << 2 * n_);
      Kx.SetSize(2 * n_);
      Vector xu, xw, Ku, Kw;
      Halves_(x, xu, xw);
      Ku.MakeRef(Kx, 0, n_);
      Kw.MakeRef(Kx, n_, n_);
      elastic_->Elasticity().ApplyStiffness(xu, Ku);
      Zuw_->Mult(xw, tmp_);
      Ku += tmp_;
      Zwu_->Mult(xu, Kw);
      Zww_->Mult(xw, tmp_);
      Kw += tmp_;
   }

   /// F of (9.2) at t = GetTime(), before the mass: the formulas of the
   /// header.  v = [v*_u; v*_w] or null (no velocity terms); d on the owned
   /// fault DOFs or null, as for NitscheWaveOperator::Force.
   void Force(const Vector &x, const Vector *v, const Vector *d, Vector &F) const
   {
      MFEM_VERIFY(x.Size() == 2 * n_, "PoroNitscheOperator::Force: x has size "
                  << x.Size() << " != 2 n = " << 2 * n_);
      MFEM_VERIFY(!v || v->Size() == 2 * n_, "PoroNitscheOperator::Force: v has "
                  "size " << v->Size() << " != 2 n = " << 2 * n_);
      const real_t t = GetTime();
      F.SetSize(2 * n_);
      Vector xu, xw, vu, vw, Fu, Fw;
      Halves_(x, xu, xw);
      if (v) { Halves_(*v, vu, vw); }
      Fu.MakeRef(F, 0, n_);
      Fw.MakeRef(F, n_, n_);
      if (src_)
      {
         MFEM_VERIFY(std::abs(src_->Time() - t) <= 1e-12 * std::max(real_t(1.0), std::abs(t)),
                     "PoroNitscheOperator::Force: the source is at t = " << src_->Time()
                     << " and the operator at t = " << t << "; advance the source "
                     "first (step 2b of plan section 5.2).");
      }
      if (V_)
      {
         V_->SetTime(t);
         fvu_->Assemble();
         fvw_->Assemble();
      }

      // u rows: the elastic force with lambda_u, then (6.6) and (6.3).
      elastic_->Force(xu, v ? &vu : nullptr, d, Fu);
      if (V_) { Fu += *fvu_; }
      Zuw_->Mult(xw, tmp_);
      Fu -= tmp_;

      // w rows.
      if (fdw_)
      {
         for (auto &g : g_) { g->SetTime(t); }
         fdw_->Assemble();
         Fw = *fdw_;
      }
      else { Fw = 0.0; }
      if (V_) { Fw += *fvw_; }
      if (has_init_) { Fw += finit_; }
      Zwu_->Mult(xu, tmp_);
      Fw -= tmp_;
      Zww_->Mult(xw, tmp_);
      Fw -= tmp_;
      if (v)
      {
         ApplyDrag(vw, tmp_);
         Fw -= tmp_;
      }
      if (src_ && src_->IsWell()) { src_->AddWellLoad(F); }
   }

   /// Dv = b_e M_e v per element and component (the w rows of D v, 8.2).
   void ApplyDrag(const Vector &vw, Vector &Dv) const
   {
      MFEM_VERIFY(vw.Size() == n_, "PoroNitscheOperator::ApplyDrag: size "
                  << vw.Size() << " != n = " << n_);
      Dv.SetSize(n_);
      const NitscheMass &mass = elastic_->Mass();
      const int ns = mass.NumScalarDofs();
      Vector xe, ye;
      for (size_t e = 0; e < me_.size(); e++)
      {
         const Array<int> &dofs = mass.ElementScalarDofs(static_cast<int>(e));
         const int nd = dofs.Size();
         xe.SetSize(nd);
         ye.SetSize(nd);
         for (int c = 0; c < 3; c++)
         {
            const int off = c * ns;
            for (int k = 0; k < nd; k++) { xe(k) = vw(off + dofs[k]); }
            me_[e].Mult(xe, ye);
            for (int k = 0; k < nd; k++) { Dv(off + dofs[k]) = b_(e) * ye(k); }
         }
      }
   }

private:
   ParFiniteElementSpace &Space_() { return elastic_->Elasticity().GetFESpace(); }

   void Halves_(const Vector &x, Vector &a, Vector &b) const
   {
      a.MakeRef(const_cast<Vector &>(x), 0, n_);
      b.MakeRef(const_cast<Vector &>(x), n_, n_);
   }

   /// Dirichlet and no-flow attributes must sit on true boundary faces, and
   /// no-flow must not be the fault.
   void CheckBoundary_() const
   {
      for (int a : pbc_.no_flow_attrs)
      {
         MFEM_VERIFY(a != bc_.fault_attr, "PoroNitscheOperator: no-flow attribute "
                     << a << " is the fault attribute; a sealed fault arrives in "
                     "step P6 (derivation section 7.5).");
         MFEM_VERIFY(a >= 1 && mesh_.bdr_attributes.Size() && a <= mesh_.bdr_attributes.Max(),
                     "PoroNitscheOperator: no-flow attribute " << a
                     << " is not a boundary attribute of the mesh.");
      }
      for (int be = 0; be < mesh_.GetNBE(); be++)
      {
         const int a = mesh_.GetBdrAttribute(be);
         const bool dir = bc_.dirichlet_attrs.count(a) > 0;
         const bool nf = pbc_.no_flow_attrs.count(a) > 0;
         if (!dir && !nf) { continue; }
         const int f = mesh_.GetBdrElementFaceIndex(be);
         MFEM_VERIFY(mesh_.GetFaceInformation(f).IsBoundary(),
                     "PoroNitscheOperator: boundary element " << be << " with "
                     << (dir ? "Dirichlet" : "no-flow") << " attribute " << a
                     << " lies on an interior face; the poroelastic option "
                     "supports " << (dir ? "Dirichlet" : "no-flow")
                     << " conditions on the outer boundary only.");
      }
   }

   ParMesh &mesh_;
   const PoroMaterialSample &sample_;
   BoundaryConfig bc_;
   PoroBoundaryConfig pbc_;
   real_t kappa_w_;
   int n_ = 0;

   std::unique_ptr<NitscheWaveOperator> elastic_;
   std::unique_ptr<L2_FECollection> l2_fec_;
   std::unique_ptr<ParFiniteElementSpace> l2_fes_;
   std::unique_ptr<ParGridFunction> alphaM_gf_;
   std::unique_ptr<GridFunctionCoefficient> c_aM_, c_M_;
   std::unique_ptr<ConstantCoefficient> c_one_;
   Array<int> dir_marker_, nf_marker_;
   bool has_dir_ = false, has_nf_ = false;

   std::unique_ptr<HypreParMatrix> Zwu_, Zuw_, Zww_;
   std::unique_ptr<DGPoroDivIntegrator> B_aM_, B_M_, B_one_;

   std::vector<std::unique_ptr<VectorFunctionCoefficient>> g_;
   std::vector<std::unique_ptr<Array<int>>> g_marker_;
   std::unique_ptr<ParLinearForm> fdw_;

   Coefficient *V_ = nullptr;
   PoroSource *src_ = nullptr;
   std::unique_ptr<ParLinearForm> fvu_, fvw_;
   Vector finit_;
   bool has_init_ = false;

   Vector b_;
   std::vector<DenseMatrix> me_;
   mutable Vector tmp_;
};

#endif // MFEM_USE_MPI

} // namespace poro
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_PORO_NITSCHE_OPERATOR_HPP

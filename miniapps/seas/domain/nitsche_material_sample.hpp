// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// nitsche_material_sample.hpp — Nitsche dynamic rupture.
//
// The per-element material sample that feeds the spatially varying
// ElasticityDomainOperator constructor (C-09), the mass (C-02) and the CFL
// bound (N-02).  A MaterialField is sampled once at every element centroid,
// exactly as BimaterialWaveOperator::BuildGodunovFluxPool_ does
// (dynamic/bimaterial_wave_operator.inl:224-248: Geometries.GetCenter,
// rho > 0 and lambda + 2 mu > 0 checks), and stored as three piecewise-
// constant (order-0 L2) ParGridFunctions with face-neighbour data
// exchanged, wrapped as GridFunctionCoefficients.
//
// (Moved unchanged from nitsche_wave_operator.hpp in P3 to break the
// nitsche_mass.hpp <-> nitsche_wave_operator.hpp include cycle.)
//
// Scope (Chunhui Zhao, 2026-10-13): the Nitsche path uses tetrahedral meshes
// only; triangles in 2D (2D in-plane design rev 3, phase C).

#ifndef MFEM_SEAS_NITSCHE_MATERIAL_SAMPLE_HPP
#define MFEM_SEAS_NITSCHE_MATERIAL_SAMPLE_HPP

#include "mfem.hpp"
#include "../dynamic/heterogeneous_material.hpp"

#include <memory>

namespace mfem
{
namespace seas
{

#ifdef MFEM_USE_MPI

/// One element's (lambda, mu, rho) by MaterialField::sampling (centroid,
/// stiffest quadrature point or SeisSol cell average over IntRules(geom,
/// 2 order)); Mode::Constant and Mode::Coefficient.  Shared by
/// NitscheMaterialSample and DynPhase::PlanLts (serial mesh).
inline void SampleElementMaterial(const MaterialField &material, ElementTransformation &Tr,
                                  int e, Geometry::Type gtype, int order,
                                  real_t &lam, real_t &mu, real_t &rho)
{
   ElementTransformation *T = &Tr;
   const IntegrationPoint &ip = Geometries.GetCenter(gtype);
   T->SetIntPoint(&ip);
   if (material.mode != MaterialField::Mode::Coefficient
            || material.sampling == MaterialField::Sampling::Centroid)
   {
      material.EvalAt(e, *T, ip, lam, mu, rho);
   }
   else
   {
      // heterogeneity plan rev4 (D11): reduce the quadrature points of
      // IntRules(geom, 2 order) to one value per element
      const IntegrationRule &ir = IntRules.Get(gtype, 2 * order);
      if (material.sampling == MaterialField::Sampling::StiffestPoint)
      {
         real_t best = -1.0;
         for (int q = 0; q < ir.GetNPoints(); ++q)
         {
            const IntegrationPoint &qp = ir.IntPoint(q);
            T->SetIntPoint(&qp);
            real_t l, m, r;
            material.EvalAt(e, *T, qp, l, m, r);
            const real_t c2 = (l + 2.0 * m) / r;
            if (c2 > best) { best = c2; lam = l; mu = m; rho = r; }
         }
      }
      else   // CellAverage: SeisSol ParameterDB.cpp computeAveragedMaterial
      {
         real_t wsum = 0.0, rho_m = 0.0, mu_inv = 0.0, ve = 0.0;
         real_t l0 = 0.0, m0 = 0.0, r0 = 0.0;
         bool uniform = true;   // equal points: return them exactly (gate U4)
         for (int q = 0; q < ir.GetNPoints(); ++q)
         {
            const IntegrationPoint &qp = ir.IntPoint(q);
            T->SetIntPoint(&qp);
            real_t l, m, r;
            material.EvalAt(e, *T, qp, l, m, r);
            if (q == 0) { l0 = l; m0 = m; r0 = r; }
            else { uniform = uniform && l == l0 && m == m0 && r == r0; }
            MFEM_VERIFY(m > 0.0, "NitscheMaterialSample: CellAverage needs mu > 0 (element "
                        << e << ")");
            const real_t w = qp.weight;
            wsum += w;
            rho_m += w * r;
            mu_inv += w / m;
            ve += w * l / (2.0 * m * (3.0 * l + 2.0 * m));
         }
         if (uniform) { lam = l0; mu = m0; rho = r0; }
         else
         {
            rho_m /= wsum; mu_inv /= wsum; ve /= wsum;
            rho = rho_m;
            mu = 1.0 / mu_inv;
            lam = 4.0 * mu * mu * ve / (1.0 - 6.0 * mu * ve);
         }
      }
   }
}

/// Per-element (lambda, mu, rho) sample of a MaterialField on a ParMesh.
class NitscheMaterialSample
{
public:
   /// @a order: the DG order, which sets the rule IntRules(geom, 2 order) of
   /// MaterialField::Sampling::StiffestPoint and ::CellAverage (the mass
   /// matrix's rule); unused for the centroid sample.
   NitscheMaterialSample(ParMesh &pmesh, const MaterialField &material, int order = 1)
      : fec_(0, pmesh.Dimension()),
        fes_(&pmesh, &fec_),
        lambda_gf_(&fes_), mu_gf_(&fes_), rho_gf_(&fes_),
        lambda_coef_(&lambda_gf_), mu_coef_(&mu_gf_), rho_coef_(&rho_gf_)
   {
      MFEM_VERIFY(material.mode == MaterialField::Mode::Constant
                  || material.mode == MaterialField::Mode::Coefficient
                  || material.mode == MaterialField::Mode::GridFunction,
                  "NitscheMaterialSample: unrecognised material.mode ("
                  << static_cast<int>(material.mode) << ").");

      const int ne = pmesh.GetNE();
      const int dim = pmesh.Dimension();
      MFEM_VERIFY(dim == 2 || dim == 3, "NitscheMaterialSample: mesh dimension "
                  << dim << "; 2 or 3 only");
      const Geometry::Type simplex = (dim == 3) ? Geometry::TETRAHEDRON
                                                : Geometry::TRIANGLE;
      for (int e = 0; e < ne; ++e)
      {
         MFEM_VERIFY(pmesh.GetElementBaseGeometry(e) == simplex,
                     "NitscheMaterialSample: element " << e << " is a "
                     << Geometry::Name[pmesh.GetElementBaseGeometry(e)]
                     << "; the Nitsche path supports tetrahedra (3D) and "
                     "triangles (2D) only.");
      }

      Array<int> dofs;
      for (int e = 0; e < ne; ++e)
      {
         ElementTransformation *T = pmesh.GetElementTransformation(e);
         const Geometry::Type   gtype = pmesh.GetElementBaseGeometry(e);
         const IntegrationPoint &ip   = Geometries.GetCenter(gtype);
         real_t lam, mu, rho;
         if (material.mode == MaterialField::Mode::GridFunction)
         {
            MFEM_VERIFY(material.lambda_gf && material.mu_gf && material.rho_gf,
                        "NitscheMaterialSample: GridFunction mode requires "
                        "three GFs.");
            lam = material.lambda_gf->GetValue(e, ip);
            mu  = material.mu_gf->GetValue(e, ip);
            rho = material.rho_gf->GetValue(e, ip);
         }
         else { SampleElementMaterial(material, *T, e, gtype, order, lam, mu, rho); }
         MFEM_VERIFY(rho > 0.0,
                     "NitscheMaterialSample: rho (" << rho
                     << ") must be > 0 at element " << e);
         MFEM_VERIFY(lam + 2.0 * mu > 0.0,
                     "NitscheMaterialSample: (lambda + 2*mu) ("
                     << (lam + 2.0 * mu) << ") must be > 0 at element "
                     << e << " (lambda=" << lam << ", mu=" << mu << ").");

         fes_.GetElementDofs(e, dofs);
         MFEM_ASSERT(dofs.Size() == 1, "order-0 L2 has one DOF per element");
         lambda_gf_(dofs[0]) = lam;
         mu_gf_(dofs[0])     = mu;
         rho_gf_(dofs[0])    = rho;
      }

      lambda_gf_.ExchangeFaceNbrData();
      mu_gf_.ExchangeFaceNbrData();
      rho_gf_.ExchangeFaceNbrData();
   }

   Coefficient &Lambda() { return lambda_coef_; }
   Coefficient &Mu() { return mu_coef_; }
   Coefficient &Rho() { return rho_coef_; }

   const ParGridFunction &LambdaGF() const { return lambda_gf_; }
   const ParGridFunction &MuGF() const { return mu_gf_; }
   const ParGridFunction &RhoGF() const { return rho_gf_; }

private:
   L2_FECollection fec_;
   ParFiniteElementSpace fes_;
   ParGridFunction lambda_gf_, mu_gf_, rho_gf_;
   GridFunctionCoefficient lambda_coef_, mu_coef_, rho_coef_;
};

#endif // MFEM_USE_MPI

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_NITSCHE_MATERIAL_SAMPLE_HPP

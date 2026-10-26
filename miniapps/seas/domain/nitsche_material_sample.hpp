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

/// Per-element (lambda, mu, rho) sample of a MaterialField on a ParMesh.
class NitscheMaterialSample
{
public:
   NitscheMaterialSample(ParMesh &pmesh, const MaterialField &material)
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
         else
         {
            material.EvalAt(e, *T, ip, lam, mu, rho);
         }
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

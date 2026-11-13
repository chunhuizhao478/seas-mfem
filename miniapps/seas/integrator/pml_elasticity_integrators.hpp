// integrator/pml_elasticity_integrators.hpp - the QD band volume integrator.
//
// PML design plan rev2, Section 4.3 (P1, Section 7.2): on band elements the
// volume term is  int C0 grad u : grad v  with the full gradient (C0 has the
// major symmetry but not the minor one, (2.22)).
//
// With g~ = Lambda0^-1 grad(phi) (Lambda0 is symmetric) the entry for test
// (k, p) and trial (l, u) is, from (4.2),
//   S0 [ lam g~_k,p g~_l,u + mu ( delta_pu g~_k . g~_l + g~_k,u g~_l,p ) ]
// which is the plain isotropic elasticity entry with the stretched shape
// gradients and the weight S0 (the map of Proposition 1).
//
// Elements with no vertex inside a band are handed to MFEM's
// ElasticityIntegrator unchanged, so the matrix outside the bands is the
// same bit for bit (test T0).  The band face terms are in
// DGElasticityIPCombinedIntegrator::SetPmlBand.

#ifndef MFEM_SEAS_PML_ELASTICITY_INTEGRATORS_HPP
#define MFEM_SEAS_PML_ELASTICITY_INTEGRATORS_HPP

#include "mfem.hpp"
#include "../domain/pml_band.hpp"

namespace mfem
{
namespace seas
{

/// True when some vertex of element `e` lies strictly inside a band.
inline bool PmlBandElement(const Mesh &mesh, int e, const PmlBand &band)
{
   Array<int> v;
   mesh.GetElementVertices(e, v);
   for (int k = 0; k < v.Size(); k++)
   {
      if (band.InBand(mesh.GetVertex(v[k]))) { return true; }
   }
   return false;
}

class PmlBandElasticityIntegrator : public BilinearFormIntegrator
{
public:
   PmlBandElasticityIntegrator(Coefficient &lambda, Coefficient &mu,
                               const PmlBand &band)
      : lambda_(lambda), mu_(mu), band_(band), plain_(lambda, mu) {}

   void AssembleElementMatrix(const FiniteElement &el, ElementTransformation &Tr,
                              DenseMatrix &elmat) override
   {
      const Mesh *mesh = Tr.mesh;
      MFEM_VERIFY(mesh, "PmlBandElasticityIntegrator: no mesh on the transformation");
      if (!PmlBandElement(*mesh, Tr.ElementNo, band_))
      {
         plain_.AssembleElementMatrix(el, Tr, elmat);   // outside the bands
         return;
      }
      const int D = el.GetDim();
      const int nd = el.GetDof();
      MFEM_VERIFY(D == band_.Dim(), "PmlBandElasticityIntegrator: dimension");
      elmat.SetSize(D * nd);
      elmat = 0.0;

      // order: today's 2p (affine) plus the profile's power, so the
      // quadratic profile is integrated as well as the rest
      const int order = 2 * el.GetOrder() + band_.Params().power;
      const IntegrationRule &ir = IntRules.Get(el.GetGeomType(), order);

      DenseMatrix dshape(nd, D), gt(nd, D);
      Vector x(D);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Tr.SetIntPoint(&ip);
         el.CalcPhysDShape(Tr, dshape);
         Tr.Transform(ip, x);
         const PmlPoint p = band_.At(x);
         const real_t lam = lambda_.Eval(Tr, ip);
         const real_t mu = mu_.Eval(Tr, ip);
         const real_t w = ip.weight * Tr.Weight() * p.S0;
         for (int k = 0; k < nd; k++)          // g~ = Lambda0^-1 grad phi
            for (int i = 0; i < D; i++)
            {
               real_t v = 0.0;
               for (int j = 0; j < D; j++) { v += p.Linv[i][j] * dshape(k, j); }
               gt(k, i) = v;
            }
         for (int k = 0; k < nd; k++)
            for (int l = 0; l < nd; l++)
            {
               real_t gg = 0.0;
               for (int i = 0; i < D; i++) { gg += gt(k, i) * gt(l, i); }
               for (int pp = 0; pp < D; pp++)
                  for (int u = 0; u < D; u++)
                  {
                     real_t v = lam * gt(k, pp) * gt(l, u) + mu * gt(k, u) * gt(l, pp);
                     if (pp == u) { v += mu * gg; }
                     elmat(pp * nd + k, u * nd + l) += w * v;
                  }
            }
      }
   }

private:
   Coefficient &lambda_;
   Coefficient &mu_;
   const PmlBand &band_;
   ElasticityIntegrator plain_;
};

/// The 2.5D finite-width mass term inside the QD band (2.5D with the PML,
/// piece F-a): int beta S0 u . v on band elements, S0 the zero-frequency
/// value of the mass-side stretch S = prod_a s_a (the term beta (u - u_pl)
/// gets the stretch of the inertia term).  The rule is the band rule 2p + N
/// of PmlBandElasticityIntegrator and of the dynamic band's volume points,
/// so the dynamic M_beta + band part at rest equals this matrix to
/// rounding.  Elements with no vertex inside a band are handed to MFEM's
/// VectorMassIntegrator unchanged.
class PmlBandBetaMassIntegrator : public BilinearFormIntegrator
{
public:
   PmlBandBetaMassIntegrator(Coefficient &beta, const PmlBand &band, int vdim)
      : beta_(beta), band_(band), vdim_(vdim), plain_(beta)
   { plain_.SetVDim(vdim); }

   void AssembleElementMatrix(const FiniteElement &el, ElementTransformation &Tr,
                              DenseMatrix &elmat) override
   {
      const Mesh *mesh = Tr.mesh;
      MFEM_VERIFY(mesh, "PmlBandBetaMassIntegrator: no mesh on the transformation");
      if (!PmlBandElement(*mesh, Tr.ElementNo, band_))
      {
         plain_.AssembleElementMatrix(el, Tr, elmat);   // outside the bands
         return;
      }
      const int D = el.GetDim();
      const int nd = el.GetDof();
      MFEM_VERIFY(D == band_.Dim() && vdim_ == D, "PmlBandBetaMassIntegrator: dimension");
      elmat.SetSize(D * nd);
      elmat = 0.0;
      const int order = 2 * el.GetOrder() + band_.Params().power;
      const IntegrationRule &ir = IntRules.Get(el.GetGeomType(), order);
      Vector shape(nd), x(D);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         Tr.SetIntPoint(&ip);
         el.CalcPhysShape(Tr, shape);
         Tr.Transform(ip, x);
         const PmlPoint p = band_.At(x);
         const real_t w = ip.weight * Tr.Weight() * beta_.Eval(Tr, ip) * p.S0;
         for (int k = 0; k < nd; k++)
            for (int l = 0; l < nd; l++)
            {
               const real_t v = w * shape(k) * shape(l);
               for (int c = 0; c < D; c++) { elmat(c * nd + k, c * nd + l) += v; }
            }
      }
   }

private:
   Coefficient &beta_;
   const PmlBand &band_;
   int vdim_;
   VectorMassIntegrator plain_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_PML_ELASTICITY_INTEGRATORS_HPP

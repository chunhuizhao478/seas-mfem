// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// poro_material.hpp — the Biot material of the poroelastic option
// (poroelastic plan rev 3, step P1.1; poro-derivation-rev1.2 §4.1, eq. (4.4)
// and (8.4)).
//
// The ten Biot inputs, under SeisSol's names (derivation §4.1): the drained
// frame moduli lambda, mu come from the elastic MaterialField (the deck's
// [material_constant_fallback] or a material file), the other eight from
// [multiphysics.poroelastic] (BiotFluidInputs; the tests also give them per
// element).  The bulk density rho is the elastic material's own rho, which
// must equal (1 - phi) rho_s + phi rho_f to 1e-12 relative (plan §4.2), so the
// solid block keeps the exact density of the elastic run.
//
//   derived constants, per element (eq. (4.4); SeisSol Helper.h forms)
//     K_d = lambda + 2 mu / 3          alpha = 1 - K_d / K_s
//     M   = K_s / (1 - phi - K_d / K_s + phi K_s / K_f)
//     lambda_u = lambda + alpha^2 M    K_u = K_d + alpha^2 M    B = alpha M / K_u
//     m = T rho_f / phi                b = eta / kappa
//     tau = (rho m - rho_f^2) / (b rho)                           eq. (8.4)
//
// PoroMaterialSample samples every element at its centroid, as
// NitscheMaterialSample does, and stores the constants as order-0 L2
// ParGridFunctions with face-neighbour data exchanged.  A second constructor
// (step P2.5) takes the fluid constants (rho_f, alpha, M, m, b) per element
// directly, for media the eight inputs cannot describe: alpha = 0, rho_f = 0
// or b = 0, which tests U2-15, U2-16 and U2-18 need; lambda_u, K_u, B and tau
// then follow from them (CompleteBiotConstants).  UndrainedSample() is a
// NitscheMaterialSample with (lambda_u, mu, rho): the elastic operator, its
// penalty, traction kernel, slip load and Kelvin-Voigt layer are built from it
// (derivation Proposition 6.2, option F2).
//
// Scope (Chunhui Zhao, 2026-10-13): tetrahedral meshes only.

#ifndef MFEM_SEAS_PORO_MATERIAL_HPP
#define MFEM_SEAS_PORO_MATERIAL_HPP

#include "mfem.hpp"
#include "../../domain/nitsche_material_sample.hpp"
#include "poro_config.hpp"

#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{
namespace poro
{

/// The eight Biot inputs of [multiphysics.poroelastic] (SI units).
struct BiotFluidInputs
{
   real_t porosity = 0.0;       ///< phi
   real_t permeability = 0.0;   ///< kappa [m^2]
   real_t viscosity = 0.0;      ///< eta [Pa s]
   real_t tortuosity = 0.0;     ///< T
   real_t rho_fluid = 0.0;      ///< rho_f [kg/m^3]
   real_t rho_solid = 0.0;      ///< rho_s, of the grains [kg/m^3]
   real_t bulk_fluid = 0.0;     ///< K_f [Pa]
   real_t bulk_solid = 0.0;     ///< K_s, of the grains [Pa]
};

/// The eight inputs from the deck table.
inline BiotFluidInputs FluidInputsFromConfig(const PoroConfig &cfg)
{
   BiotFluidInputs f;
   f.porosity = cfg.porosity;
   f.permeability = cfg.permeability_m2;
   f.viscosity = cfg.viscosity_pa_s;
   f.tortuosity = cfg.tortuosity;
   f.rho_fluid = cfg.rho_fluid;
   f.rho_solid = cfg.rho_solid;
   f.bulk_fluid = cfg.bulk_fluid_pa;
   f.bulk_solid = cfg.bulk_solid_pa;
   return f;
}

/// (1 - phi) rho_s + phi rho_f, in the order of PoroConfig::BulkDensity.
inline real_t BulkDensity(const BiotFluidInputs &f)
{
   return (1.0 - f.porosity) * f.rho_solid + f.porosity * f.rho_fluid;
}

/// The constants of one element: the elastic inputs and eq. (4.4), (8.4).
struct BiotConstants
{
   real_t lambda = 0.0, mu = 0.0, rho = 0.0;   ///< drained frame, bulk density
   real_t rho_f = 0.0;                          ///< fluid density
   real_t Kd = 0.0;                             ///< drained bulk modulus
   real_t alpha = 0.0;                          ///< Biot coefficient
   real_t M = 0.0;                              ///< Biot modulus
   real_t lambda_u = 0.0;                       ///< undrained Lame modulus
   real_t Ku = 0.0;                             ///< undrained bulk modulus
   real_t B = 0.0;                              ///< Skempton coefficient
   real_t m = 0.0;                              ///< apparent fluid mass T rho_f / phi
   real_t b = 0.0;                              ///< drag eta / kappa
   real_t tau = 0.0;                            ///< relaxation time, eq. (8.4)
};

/// Eq. (4.4) and (8.4) from the ten inputs, with rho the bulk density.
inline BiotConstants DeriveBiotConstants(real_t lambda, real_t mu, real_t rho,
                                         const BiotFluidInputs &f)
{
   BiotConstants c;
   c.lambda = lambda;
   c.mu = mu;
   c.rho = rho;
   c.rho_f = f.rho_fluid;
   const real_t phi = f.porosity, Ks = f.bulk_solid, Kf = f.bulk_fluid;
   c.Kd = lambda + 2.0 / 3.0 * mu;
   c.alpha = 1.0 - c.Kd / Ks;
   c.M = Ks / (1.0 - phi - c.Kd / Ks + phi * Ks / Kf);
   c.lambda_u = lambda + c.alpha * c.alpha * c.M;
   c.Ku = c.Kd + c.alpha * c.alpha * c.M;
   c.B = c.alpha * c.M / c.Ku;
   c.m = f.tortuosity * f.rho_fluid / phi;
   c.b = f.viscosity / f.permeability;
   c.tau = (rho * c.m - f.rho_fluid * f.rho_fluid) / (c.b * rho);
   return c;
}

/// Empty when the inputs make a valid Biot medium, otherwise the reason.
inline std::string CheckBiotInputs(real_t lambda, real_t mu, real_t rho,
                                   const BiotFluidInputs &f)
{
   std::ostringstream err;
   if (!(f.porosity > 0.0 && f.porosity < 1.0))
   {
      err << "porosity " << f.porosity << " is not in (0, 1)";
   }
   else if (!(f.permeability > 0.0 && f.viscosity > 0.0 && f.rho_fluid > 0.0
              && f.rho_solid > 0.0 && f.bulk_fluid > 0.0 && f.bulk_solid > 0.0))
   {
      err << "permeability, viscosity, rho_fluid, rho_solid, bulk_fluid and "
          << "bulk_solid must be > 0";
   }
   else if (!(f.tortuosity >= 1.0))
   {
      err << "tortuosity " << f.tortuosity << " is below 1";
   }
   else if (!(mu > 0.0 && lambda + 2.0 / 3.0 * mu > 0.0))
   {
      err << "the drained frame needs mu > 0 and K_d = lambda + 2 mu / 3 > 0 "
          << "(lambda " << lambda << ", mu " << mu << ")";
   }
   else if (!(std::abs(rho - BulkDensity(f)) <= 1e-12 * std::abs(rho)))
   {
      err.precision(17);
      err << "rho " << rho << " differs from the bulk density (1 - phi) rho_s"
          << " + phi rho_f = " << BulkDensity(f) << " by more than 1e-12 relative";
   }
   else
   {
      const BiotConstants c = DeriveBiotConstants(lambda, mu, rho, f);
      if (!(c.alpha > 0.0 && c.alpha <= 1.0))
      {
         err << "Biot coefficient alpha = 1 - K_d / K_s = " << c.alpha
             << " is not in (0, 1] (K_d " << c.Kd << ", K_s " << f.bulk_solid << ")";
      }
      else if (!(c.M > 0.0))
      {
         err << "Biot modulus M = " << c.M << " is not > 0";
      }
      else if (!(rho * c.m - f.rho_fluid * f.rho_fluid > 0.0))
      {
         err << "rho m - rho_f^2 = " << rho * c.m - f.rho_fluid * f.rho_fluid
             << " is not > 0";
      }
   }
   return err.str();
}

/// The fluid constants of one element given directly (step P2.5), for media
/// the eight inputs cannot describe (alpha = 0, rho_f = 0, b = 0).
struct BiotFluidConstants
{
   real_t rho_f = 0.0;   ///< fluid density, >= 0
   real_t alpha = 0.0;   ///< Biot coefficient, in [0, 1]
   real_t M = 0.0;       ///< Biot modulus, > 0
   real_t m = 0.0;       ///< apparent fluid mass, > 0
   real_t b = 0.0;       ///< drag, >= 0
};

/// The constants of one element from the elastic inputs and the given fluid
/// constants: K_d, lambda_u, K_u, B by eq. (4.4) and tau by (8.4), in the
/// operation order of DeriveBiotConstants; tau = +inf when b = 0.
inline BiotConstants CompleteBiotConstants(real_t lambda, real_t mu, real_t rho,
                                           const BiotFluidConstants &k)
{
   BiotConstants c;
   c.lambda = lambda;
   c.mu = mu;
   c.rho = rho;
   c.rho_f = k.rho_f;
   c.Kd = lambda + 2.0 / 3.0 * mu;
   c.alpha = k.alpha;
   c.M = k.M;
   c.lambda_u = lambda + c.alpha * c.alpha * c.M;
   c.Ku = c.Kd + c.alpha * c.alpha * c.M;
   c.B = c.alpha * c.M / c.Ku;
   c.m = k.m;
   c.b = k.b;
   c.tau = c.b > 0.0 ? (rho * c.m - c.rho_f * c.rho_f) / (c.b * rho)
           : std::numeric_limits<real_t>::infinity();
   return c;
}

/// Empty when the elastic inputs and the given fluid constants make a valid
/// mass, drag and stiffness, otherwise the reason.
inline std::string CheckBiotConstants(real_t lambda, real_t mu, real_t rho,
                                      const BiotFluidConstants &k)
{
   std::ostringstream err;
   if (!(mu > 0.0 && lambda + 2.0 / 3.0 * mu > 0.0 && rho > 0.0))
   {
      err << "the frame needs mu > 0, K_d = lambda + 2 mu / 3 > 0 and rho > 0 (lambda "
          << lambda << ", mu " << mu << ", rho " << rho << ")";
   }
   else if (!(k.alpha >= 0.0 && k.alpha <= 1.0))
   {
      err << "Biot coefficient alpha = " << k.alpha << " is not in [0, 1]";
   }
   else if (!(k.M > 0.0 && k.m > 0.0))
   {
      err << "Biot modulus M = " << k.M << " and fluid mass m = " << k.m
          << " must be > 0";
   }
   else if (!(k.rho_f >= 0.0 && k.b >= 0.0))
   {
      err << "rho_f = " << k.rho_f << " and b = " << k.b << " must be >= 0";
   }
   else if (!(rho * k.m - k.rho_f * k.rho_f > 0.0))
   {
      err << "rho m - rho_f^2 = " << rho * k.m - k.rho_f * k.rho_f << " is not > 0";
   }
   return err.str();
}

#ifdef MFEM_USE_MPI

/// Per-element Biot constants of a ParMesh (eq. (4.4), (8.4)), and the
/// undrained sample for the elastic block.
class PoroMaterialSample
{
public:
   /// The eight fluid inputs of local element e, whose centroid is x.
   using FluidField = std::function<BiotFluidInputs(int e, const Vector &x)>;

   /// Uniform fluid inputs from the deck table.
   PoroMaterialSample(ParMesh &pmesh, const MaterialField &material,
                      const PoroConfig &cfg)
      : PoroMaterialSample(pmesh, material, UniformField(FluidInputsFromConfig(cfg)))
   { }

   PoroMaterialSample(ParMesh &pmesh, const MaterialField &material,
                      const FluidField &fluid)
      : PoroMaterialSample(pmesh, material, ElementFn(
                              [&fluid](int e, const Vector &x, real_t lam, real_t mu,
                                       real_t rho, BiotConstants &c)
   {
      const BiotFluidInputs f = fluid(e, x);
      const std::string err = CheckBiotInputs(lam, mu, rho, f);
      if (err.empty()) { c = DeriveBiotConstants(lam, mu, rho, f); }
      return err;
   }), ByElement{})
   { }

   /// The fluid constants of local element e, whose centroid is x (step P2.5).
   using ConstantsField = std::function<BiotFluidConstants(int e, const Vector &x)>;

   /// The fluid constants given directly (CheckBiotConstants,
   /// CompleteBiotConstants); the elastic inputs come from `material` as in
   /// the other constructors.
   PoroMaterialSample(ParMesh &pmesh, const MaterialField &material,
                      const ConstantsField &constants)
      : PoroMaterialSample(pmesh, material, ElementFn(
                              [&constants](int e, const Vector &x, real_t lam, real_t mu,
                                           real_t rho, BiotConstants &c)
   {
      const BiotFluidConstants k = constants(e, x);
      const std::string err = CheckBiotConstants(lam, mu, rho, k);
      if (err.empty()) { c = CompleteBiotConstants(lam, mu, rho, k); }
      return err;
   }), ByElement{})
   { }

private:
   /// The constants of element e (centroid x) from its elastic inputs, or the
   /// reason they are refused.
   using ElementFn = std::function<std::string(int e, const Vector &x, real_t lam,
                                               real_t mu, real_t rho, BiotConstants &c)>;
   struct ByElement { };

   PoroMaterialSample(ParMesh &pmesh, const MaterialField &material,
                      const ElementFn &element, ByElement)
      : drained_(pmesh, material),
        fec_(0, pmesh.Dimension()),
        fes_(&pmesh, &fec_),
        rho_f_gf_(&fes_), m_gf_(&fes_), b_gf_(&fes_), alpha_gf_(&fes_),
        M_gf_(&fes_), tau_gf_(&fes_),
        lambda_u_gf_(std::make_shared<ParGridFunction>(&fes_)),
        mu_gf_(std::make_shared<ParGridFunction>(&fes_)),
        rho_gf_(std::make_shared<ParGridFunction>(&fes_))
   {
      const ParFiniteElementSpace *dfes = drained_.RhoGF().ParFESpace();
      const int ne = pmesh.GetNE();
      constants_.resize(ne);
      Array<int> dofs, ddofs;
      Vector x(pmesh.SpaceDimension());
      for (int e = 0; e < ne; ++e)
      {
         // The elastic inputs exactly as the elastic run samples them.
         dfes->GetElementDofs(e, ddofs);
         const real_t lam = drained_.LambdaGF()(ddofs[0]);
         const real_t mu  = drained_.MuGF()(ddofs[0]);
         const real_t rho = drained_.RhoGF()(ddofs[0]);

         ElementTransformation *T = pmesh.GetElementTransformation(e);
         T->Transform(Geometries.GetCenter(pmesh.GetElementBaseGeometry(e)), x);
         BiotConstants c;
         const std::string err = element(e, x, lam, mu, rho, c);
         MFEM_VERIFY(err.empty(), "PoroMaterialSample: element " << e
                     << " (centroid " << x(0) << ", " << x(1) << ", " << x(2)
                     << "): " << err);
         constants_[e] = c;

         fes_.GetElementDofs(e, dofs);
         MFEM_ASSERT(dofs.Size() == 1, "order-0 L2 has one DOF per element");
         const int k = dofs[0];
         rho_f_gf_(k) = c.rho_f;
         m_gf_(k) = c.m;
         b_gf_(k) = c.b;
         alpha_gf_(k) = c.alpha;
         M_gf_(k) = c.M;
         tau_gf_(k) = c.tau;
         (*lambda_u_gf_)(k) = c.lambda_u;
         (*mu_gf_)(k) = mu;
         (*rho_gf_)(k) = rho;
      }
      for (ParGridFunction *gf : {&rho_f_gf_, &m_gf_, &b_gf_, &alpha_gf_, &M_gf_,
                                  &tau_gf_, lambda_u_gf_.get(), mu_gf_.get(),
                                  rho_gf_.get()})
      {
         gf->ExchangeFaceNbrData();
      }
      // The order-0 shape function is 1, so the undrained sample reads these
      // values back unchanged (test U1-4).
      undrained_ = std::make_unique<NitscheMaterialSample>(
                      pmesh, MaterialField::MakeGridFunction(rho_gf_, lambda_u_gf_, mu_gf_));
   }

public:
   /// (lambda, mu, rho) as the elastic run samples them.
   const NitscheMaterialSample &DrainedSample() const { return drained_; }
   /// (lambda_u, mu, rho) for the elastic block (Proposition 6.2).
   NitscheMaterialSample &UndrainedSample() { return *undrained_; }
   const NitscheMaterialSample &UndrainedSample() const { return *undrained_; }
   /// (lambda_u, mu, rho) as the order-0 MaterialField UndrainedSample() reads,
   /// for building the elastic operator of the solid block (step P2.2): the
   /// operator's own sample reads back the same values.
   MaterialField UndrainedMaterial() const
   {
      return MaterialField::MakeGridFunction(rho_gf_, lambda_u_gf_, mu_gf_);
   }

   /// All constants of local element e.
   const BiotConstants &Constants(int e) const { return constants_[e]; }

   const ParGridFunction &RhoFluidGF() const { return rho_f_gf_; }
   const ParGridFunction &FluidMassGF() const { return m_gf_; }
   const ParGridFunction &DragGF() const { return b_gf_; }
   const ParGridFunction &AlphaGF() const { return alpha_gf_; }
   const ParGridFunction &BiotModulusGF() const { return M_gf_; }
   const ParGridFunction &TauGF() const { return tau_gf_; }
   const ParGridFunction &LambdaUGF() const { return *lambda_u_gf_; }

   /// The same inputs on every element.
   static FluidField UniformField(const BiotFluidInputs &f)
   {
      return [f](int, const Vector &) { return f; };
   }

private:
   NitscheMaterialSample drained_;
   L2_FECollection fec_;
   ParFiniteElementSpace fes_;
   ParGridFunction rho_f_gf_, m_gf_, b_gf_, alpha_gf_, M_gf_, tau_gf_;
   std::shared_ptr<ParGridFunction> lambda_u_gf_, mu_gf_, rho_gf_;
   std::vector<BiotConstants> constants_;
   std::unique_ptr<NitscheMaterialSample> undrained_;
};

#endif // MFEM_USE_MPI

} // namespace poro
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_PORO_MATERIAL_HPP

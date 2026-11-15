// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Unit tests of the damaged poroelastic properties (multiphysics/poroelastic/poro_damage.hpp)
// with the cmame_revision 3d_hydromech/em_5J parameters.

#include "mfem.hpp"
#include "../../multiphysics/poroelastic/poro_damage.hpp"
#include "../test_framework.hpp"

using namespace mfem;
using namespace mfem::pff;

static PoroDamageParams CmameParams()
{
   PoroDamageParams P;
   const real_t E = 50e9, nu = 0.3;
   P.K = E / 3.0 / (1.0 - 2.0 * nu);
   P.K_s = 50e9;
   P.K_f = 2.24e9;
   P.phi0 = 0.008;
   P.phi_lb = 0.008;
   P.phi_ub = 0.999;
   P.g_min = 1e-6;
   P.viscosity = 1e-3;
   P.perm_model = PermeabilityModel::NormalStrain;
   P.k0 = 5e-19;
   P.perm_exponent = 10;
   P.d_threshold = 0.5;
   P.fc = 1.0;
   P.anisotropic = true;
   return P;
}

static void TestBiotAndPorosity()
{
   std::cout << "\n[Biot coefficient, porosity, Biot modulus]\n";
   const PoroDamageParams P = CmameParams();
   TEST_NEAR(BiotCoefficient(P, 0.0), 1.0 - P.K / P.K_s, 1e-15, "alpha(0) = 1 - K/K_s");
   TEST_NEAR(BiotCoefficient(P, 1.0), 1.0 - 1e-6 * P.K / P.K_s, 1e-15, "alpha(1) = 1 - g_min K/K_s");
   TEST_NEAR(BiotCoefficient(P, 1.5), BiotCoefficient(P, 1.0), 0.0, "d clamped to [0,1]");
   TEST_NEAR(Porosity(P, 0.0), 0.008, 1e-16, "phi(0) = phi0");
   TEST_NEAR(Porosity(P, 1.0), 0.999, 1e-16, "phi(1) clamped to phi_ub");
   TEST_NEAR(Porosity(P, 0.5), 0.008 + 0.992 * 0.75, 1e-15, "phi(0.5) = phi0 + (1-phi0)(1-(1-d)^2)");
   PoroDamageParams Ps = P; Ps.porosity_model = PorosityModel::Strain;
   TEST_NEAR(Porosity(Ps, 0.9, 1e-3), 0.009, 1e-16, "strain model: phi0 + eps_1");
   TEST_NEAR(Porosity(Ps, 0.9, -1e-3), 0.008, 1e-16, "strain model: clamped at phi_lb");
   const real_t a0 = BiotCoefficient(P, 0.0), phi = Porosity(P, 0.0);
   TEST_REL(BiotModulus(P, a0, phi), 1.0 / (phi / P.K_f + (a0 - phi) / P.K_s), 1e-14,
            "1/M = phi/K_f + (alpha - phi)/K_s");
}

static void TestHistory()
{
   std::cout << "\n[enhanced history energy]\n";
   const PoroDamageParams P = CmameParams();
   const real_t cp = HistoryPressureCoefficient(P);
   TEST_REL(cp, (0.008 - 1.0) * (1.0/2.24e9 - 1.0/50e9) - P.K / (50e9*50e9), 1e-14,
            "c_p = (phi0-1)(1/K_f - 1/K_s) - K/K_s^2");
   TEST_ASSERT(cp < 0.0, "c_p < 0: pore pressure lowers the driving energy");
   const real_t p0 = 0.0965e6;   // initial pore pressure of the inputs
   TEST_NEAR(EnhancedHistory(P, 1e3, p0, 0.0), std::max(0.0, 2e3 + cp*p0*p0), 1e-9, "H = <2 psi+ + c_p p^2>_+");
   TEST_NEAR(EnhancedHistory(P, 0.0, 1e8, 0.0), 0.0, 0.0, "Macaulay bracket: H >= 0");
   TEST_NEAR(EnhancedHistory(P, 1.0, 0.0, 5.0), 5.0, 0.0, "H never decreases");
}

static void TestPermeability()
{
   std::cout << "\n[normal-strain permeability (Heider 2021)]\n";
   const PoroDamageParams P = CmameParams();
   real_t k[9];
   const real_t grad[3] = {3.0, 0.0, 0.0};
   real_t eps[9] = {2e-3, 0, 0, 0, -1e-4, 0, 0, 0, 0};
   const real_t h = 5e-5;
   EffectivePermeability(P, 0.4, grad, eps, h, k);
   TEST_NEAR(k[0], P.k0, 0.0, "d < d_th: k = k0 I");
   EffectivePermeability(P, 0.8, grad, eps, h, k);
   const real_t w = h * (1.0 + 2e-3), kw = w * w / 12.0, wt = std::pow(0.8, 10);
   TEST_REL(k[0], P.k0, 1e-14, "normal direction (n = e_x): k_xx = k0");
   TEST_REL(k[4], P.k0 + wt * kw, 1e-12, "in-plane: k_yy = k0 + d^n (h(1+eps_nn))^2/12");
   TEST_REL(k[8], P.k0 + wt * kw, 1e-12, "in-plane: k_zz likewise");
   TEST_NEAR(k[3], 0.0, 0.0, "off-diagonal zero for n = e_x");
   const real_t zero[3] = {0.0, 0.0, 0.0};
   EffectivePermeability(P, 0.8, zero, eps, h, k);
   TEST_NEAR(k[4], P.k0, 0.0, "no damage gradient and no regularisation: k = k0 I");
   PoroDamageParams Pd = P; Pd.perm_model = PermeabilityModel::DarcyPoiseuille; Pd.wc = 1e-4;
   EffectivePermeability(Pd, 0.5, grad, eps, h, k);
   TEST_REL(k[0], Pd.k0 + std::pow(0.5, 10) * (std::pow(0.5e-4, 2) / 12.0 - Pd.k0), 1e-12,
            "darcy_poiseuille: k0 + d^n (w^2/12 - k0), w = d w_c");
}

int main()
{
   TestBiotAndPorosity();
   TestHistory();
   TestPermeability();
   return TestSummary("test_poro_damage");
}

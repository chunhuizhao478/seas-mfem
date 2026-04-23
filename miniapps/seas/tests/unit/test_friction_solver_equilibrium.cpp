// TPV102 friction-solver equilibrium test.  At the TPV102 background
// equilibrium (Theta = tau_ini, sigma_n = sigma_n, psi = ComputeInitialPsi(a)),
// the solver MUST return V = V_ini.  The pepper-bug investigation found
// the solver returning V = 0.116 m/s instead — orders of magnitude wrong.
//
// This test isolates the friction solver from all DG / wave-operator
// machinery to confirm the solver itself is the bug.

#include "mfem.hpp"
#include "../../dynamic/friction_solver.hpp"
#include "../../config/tpv102_params.hpp"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>

using namespace mfem;
using namespace mfem::seas;

int main()
{
   std::cout << "\n=== TPV102 friction-solver equilibrium test ===\n";

   const real_t a = TPV102Params::a_vw;
   const real_t Dc = TPV102Params::Dc;
   const real_t sigma_n = TPV102Params::sigma_n;
   const real_t tau_ini = TPV102Params::tau_ini;
   const real_t V_ini = TPV102Params::V_ini;
   const real_t eta_s = TPV102Params::eta_s;
   const real_t psi_eq = ComputeInitialPsi(a);

   std::cout << "  Inputs:\n"
             << "    a       = " << std::scientific << std::setprecision(6)
             << a << "\n"
             << "    Dc      = " << Dc << "\n"
             << "    sigma_n = " << sigma_n << " Pa\n"
             << "    tau_ini = " << tau_ini << " Pa\n"
             << "    V_ini   = " << V_ini << " m/s\n"
             << "    eta_s   = " << eta_s << " Pa·s/m\n"
             << "    psi_eq  = " << psi_eq << "  (from ComputeInitialPsi)\n";

   // Solve the friction equation with TPV102-equilibrium inputs.
   FrictionSolver solver;
   const real_t V_back = solver.Solve(tau_ini, psi_eq, sigma_n, eta_s, a,
                                      FrictionSolver::Method::Brent);
   std::cout << "\n  solver.Solve(Theta=tau_ini, psi=psi_eq, sigma_n, "
             << "eta_s, a) = "
             << std::scientific << std::setprecision(6) << V_back << " m/s\n"
             << "  expected   = " << V_ini << " m/s\n";

   // The solver should return V = V_ini at equilibrium.  If it returns
   // something else (e.g. 0.116 m/s as observed in the wave operator),
   // the friction solver is broken at this configuration.
   const real_t rel_err = std::abs(V_back - V_ini) / std::max(V_ini, 1e-30);
   std::cout << "  relative error = " << rel_err << "\n";

   if (rel_err < 1e-3)
   {
      std::cout << "  PASS: friction solver returns V_ini at equilibrium\n";
   }
   else
   {
      std::cout << "  FAIL: friction solver does NOT return V_ini\n";
   }

   // Reproduce the observed psi from the wave operator (0.7355) to see
   // if a slightly-off psi causes the wrong V.
   const real_t psi_observed = 0.7355;
   const real_t V_observed_psi = solver.Solve(tau_ini, psi_observed, sigma_n,
                                              eta_s, a,
                                              FrictionSolver::Method::Brent);
   std::cout << "\n  Re-solve with psi = " << psi_observed
             << " (the value printed by [C-1 EVAL-TOTAL]):\n"
             << "    V_back = " << std::scientific << std::setprecision(6)
             << V_observed_psi << " m/s  (expected ~V_ini = " << V_ini
             << ")\n";

   return 0;
}

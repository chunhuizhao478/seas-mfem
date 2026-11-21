// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dashpot_reflection.cpp — C1 gate 3: a plane P (or S) wave at normal
// incidence on a Lysmer dashpot face reflects less than 2 % of its amplitude.
//
//   ./dashpot_reflection [-w p|s] [-nx 400] [-ny 4]
//
// Strip [0, L] x [0, H] of triangles, acrylic (E 2.8 GPa, nu 0.35, rho 1190).
// The other component is held at 0 everywhere, so the wave is exactly plane
// (P: u_y = 0, speed c_p; S: u_x = 0, speed c_s). Dashpots on x = L only; the
// wave starts at x0 = L/3 travelling in +x: u = A exp(-((x - x0)/w)^2),
// v = -c du/dx. After it has crossed x = L the energy left in the strip is
// the reflected energy; the amplitude ratio is sqrt(E_left / E_0).
// Explicit central difference as in pff_driver (cfl 0.4, zeta = 0).

#include "mfem.hpp"
#include "../../domain/explicit_elasticity_operator.hpp"
#include "../../solver/central_difference.hpp"
#include "../../solver/cfl.hpp"

#include <cmath>
#include <iostream>

using namespace mfem;
using namespace mfem::pff;

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   const char *wave = "p";
   int nx = 400, ny = 4;
   OptionsParser args(argc, argv);
   args.AddOption(&wave, "-w", "--wave", "p or s.");
   args.AddOption(&nx, "-nx", "--nx", "Elements along x.");
   args.AddOption(&ny, "-ny", "--ny", "Elements along y.");
   args.Parse();
   if (!args.Good()) { args.PrintUsage(mfem::out); return 1; }
   const bool P = std::string(wave) == "p";

   const real_t L = 0.1, H = L * ny / nx;
   Mesh serial = Mesh::MakeCartesian2D(nx, ny, Element::TRIANGLE, true, L, H);
   const int max_attr = serial.bdr_attributes.Max();
   ParMesh mesh(MPI_COMM_WORLD, serial);
   serial.Clear();
   H1_FECollection h1(1, 2);
   ParFiniteElementSpace V(&mesh, &h1, 2, Ordering::byNODES), S(&mesh, &h1);

   ExplicitElasticityParams p;
   const real_t E = 2.8e9, nu = 0.35;
   p.K = E / 3.0 / (1.0 - 2.0 * nu);
   p.G = E / 2.0 / (1.0 + nu);
   p.density = 1190.0;
   p.decomposition = Decomposition::Spectral;
   p.pressure_marker.SetSize(max_attr); p.pressure_marker = 0;
   p.absorbing_marker.SetSize(max_attr); p.absorbing_marker = 0;
   p.absorbing_marker[2 - 1] = 1;   // MakeCartesian2D: attribute 2 is x = L
   PulseLoadParams pp;
   PulseLoad pulse(pp);

   const int nt = S.GetTrueVSize(), n = V.GetTrueVSize();
   Array<int> fixed;
   const int other = P ? 1 : 0;
   for (int i = 0; i < nt; i++) { fixed.Append(other * nt + i); }
   ExplicitElasticityOperator op(V, S, p, pulse, fixed);
   const real_t c = P ? p.Cp() : p.Cs();

   const real_t x0 = L / 3.0, w = L / 40.0, A = 1e-6;
   const int comp = P ? 0 : 1;
   ParGridFunction ug(&S), vg(&S);
   FunctionCoefficient uf([&](const Vector &X) { const real_t s = (X(0) - x0) / w; return A * std::exp(-s * s); });
   FunctionCoefficient vf([&](const Vector &X)
   {
      const real_t s = (X(0) - x0) / w;
      return -c * A * std::exp(-s * s) * (-2.0 * s / w);
   });
   ug.ProjectCoefficient(uf);
   vg.ProjectCoefficient(vf);
   Vector ut(nt), vt(nt);
   ug.GetTrueDofs(ut);
   vg.GetTrueDofs(vt);

   CentralDifferenceState st;
   st.u.SetSize(n); st.u = 0.0;
   st.v_half.SetSize(n); st.v_half = 0.0;
   for (int i = 0; i < nt; i++) { st.u(comp * nt + i) = ut(i); st.v_half(comp * nt + i) = vt(i); }
   Vector D(nt);
   D = 0.0;
   st.F_prev.SetSize(n);
   op.InternalForce(st.u, D, st.F_prev);   // F^{-1} := F^0 (no KV at the first step)

   const real_t omega = std::sqrt(op.MaxEigenvalue(60));
   const real_t dt = 0.4 * CriticalStep(omega, 0.0);
   CentralDifference cd(op);
   CentralDifferenceStepInfo info;

   real_t psi0 = 0.0;
   Vector F0(n);
   op.InternalForce(st.u, D, F0, &psi0);
   const real_t E0 = op.KineticEnergy(st.v_half) + psi0;

   // Travel 2 (L - x0) / c: the incident pulse is fully through x = L.
   const real_t T = 2.0 * (L - x0) / c;
   const int steps = static_cast<int>(std::ceil(T / dt));
   for (int k = 0; k < steps; k++) { cd.Step(st, D, 0.0, 0.0, dt, info); }
   real_t psi = 0.0;
   Vector F(n);
   op.InternalForce(st.u, D, F, &psi);
   const real_t E1 = op.KineticEnergy(info.v_n) + psi;
   const real_t ratio = std::sqrt(std::max(E1, 0.0) / E0);
   if (Mpi::Root())
   {
      mfem::out << "dashpot_reflection: " << (P ? "P" : "S") << " wave, c = " << c
                << " m/s, h = " << L / nx << ", pulse width " << w << ", dt = " << dt
                << ", " << steps << " steps\n  E0 = " << E0 << ", E_left = " << E1
                << ", amplitude ratio sqrt(E_left/E0) = " << ratio << "\n"
                << (ratio < 0.02 ? "PASS" : "FAIL") << " (gate < 2 %)\n";
   }
   return ratio < 0.02 ? 0 : 1;
}

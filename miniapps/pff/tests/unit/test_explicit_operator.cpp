// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_explicit_operator.cpp — checks of the explicit elasticity pieces
// (kernels/split_force_kernels.hpp, domain/explicit_elasticity_operator.hpp)
// on small triangle and tetrahedron meshes, np = 1:
//   - F_int at d = 0 is linear and symmetric (x . F(y) = y . F(x)), and
//     int psi = 1/2 u . F(u);
//   - F_int equals the assembled SplitElasticityIntegrator (the qs path) for
//     the spectral split with a nonzero damage field (independent code path);
//   - lumped mass sums to rho |Omega|;
//   - the dashpot blocks on a face x = const sum to rho c_p |face| (xx) and
//     rho c_s |face| (yy), with no xy coupling (FarmsNonReflectDashpotBC);
//   - the pressure vector of a constant Peak on one face is -P n |face|.
//   - the central difference with zeta = 0 and no dashpots conserves the
//     discrete energy H = 1/2 v^{n-1/2} M v^{n-1/2} + 1/2 u^n K u^{n-1}.
//   - mass-proportional damping (quench = "mass", C = alpha M): a rigid
//     translation decays by the closed-form factor of the half-step
//     velocity update, v^{N-1/2} = v0 (1 - a/2) ((1 - a/2)/(1 + a/2))^{N-1},
//     a = alpha dt, close to exp(-alpha t), and KE + mass_damping_work = KE0.
//   - SplitElasticityIntegrator::AssembleElementGrad (spectral, d != 0) =
//     central differences of AssembleElementVector (the qs Newton tangent,
//     RACCOON mode II), on a quad and a triangle.

#include "mfem.hpp"
#include "../test_framework.hpp"
#include "../../domain/explicit_elasticity_operator.hpp"
#include "../../solver/central_difference.hpp"
#include "../../solver/cfl.hpp"
#include "../../integrator/split_elasticity_integrator.hpp"

using namespace mfem;
using namespace mfem::pff;

static ExplicitElasticityParams Params(int max_attr)
{
   ExplicitElasticityParams p;
   const real_t E = 2.8e9, nu = 0.35;
   p.K = E / 3.0 / (1.0 - 2.0 * nu);
   p.G = E / 2.0 / (1.0 + nu);
   p.density = 1190.0;
   p.decomposition = Decomposition::Spectral;
   p.pressure_marker.SetSize(max_attr);
   p.pressure_marker = 0;
   p.absorbing_marker.SetSize(max_attr);
   p.absorbing_marker = 0;
   return p;
}

static void Check(Mesh &serial, const char *label)
{
   std::cout << "\n[" << label << "]\n";
   const int dim = serial.Dimension();
   const int max_attr = serial.bdr_attributes.Max();
   ParMesh mesh(MPI_COMM_WORLD, serial);
   H1_FECollection h1(1, dim);
   ParFiniteElementSpace V(&mesh, &h1, dim, Ordering::byNODES), S(&mesh, &h1);

   // Cartesian meshes: attribute 2 is x = xmax (right) in 2D, attribute 3 in 3D.
   const int right = (dim == 2) ? 2 : 3;
   ExplicitElasticityParams p = Params(max_attr);
   p.absorbing_marker[right - 1] = 1;
   p.pressure_marker[right - 1] = 1;
   PulseLoadParams pp;
   pp.peak_pressure = 1e6;
   PulseLoad pulse(pp);
   Array<int> fixed;
   ExplicitElasticityOperator op(V, S, p, pulse, fixed);
   const int n = V.GetTrueVSize(), nt = S.GetTrueVSize();

   // Element volumes are positive.
   real_t vol = 0.0;
   for (int e = 0; e < mesh.GetNE(); e++) { vol += mesh.GetElementVolume(e); }
   TEST_NEAR(op.LumpedMass().Sum(), p.density * vol, 1e-9 * p.density * vol,
             "lumped mass sums to rho |Omega|");

   Vector x(n), y(n), Fx(n), Fy(n), D0(nt);
   D0 = 0.0;
   x.Randomize(1); y.Randomize(2);
   x *= 1e-4; y *= 1e-4;
   real_t psi = 0.0;
   op.InternalForce(x, D0, Fx, &psi);
   op.InternalForce(y, D0, Fy);
   const real_t a = x * Fy, b = y * Fx;
   TEST_REL(a, b, 1e-12, "K symmetric at d = 0: x.K y = y.K x");
   TEST_REL(psi, 0.5 * (x * Fx), 1e-12, "int psi = 1/2 u.K u at d = 0");
   Vector xy(n), Fxy(n);
   add(x, 2.0, y, xy);
   op.InternalForce(xy, D0, Fxy);
   Vector lin(n);
   add(Fx, 2.0, Fy, lin);
   lin -= Fxy;
   TEST_ASSERT(lin.Normlinf() <= 1e-10 * Fxy.Normlinf(), "F_int linear at d = 0");

   // Against the qs integrator with damage (spectral split, nonlinear).
   ParGridFunction dg(&S);
   FunctionCoefficient dfun([](const Vector &X) { return 0.4 + 0.3 * std::sin(7.0 * X(0)) * X(1); });
   dg.ProjectCoefficient(dfun);
   Vector D(nt);
   dg.GetTrueDofs(D);
   Vector Fk(n), Fq(n);
   op.InternalForce(x, D, Fk);
   ParNonlinearForm N(&V);
   auto *integ = new SplitElasticityIntegrator(p.K, p.G, Decomposition::Spectral, p.deg, dg);
   integ->SetIntRule(&IntRules.Get(mesh.GetTypicalElementGeometry(), 2));
   N.AddDomainIntegrator(integ);
   N.Mult(x, Fq);
   Fq -= Fk;
   TEST_ASSERT(Fq.Normlinf() <= 1e-10 * Fk.Normlinf(),
               "kernel F_int = SplitElasticityIntegrator (spectral, d != 0)");

   // Dashpot on the right face: unit-speed x then y velocity.
   const real_t face = (dim == 2) ? 1.0 : 1.0;   // unit square / cube faces
   Vector v(n), f(n);
   v = 0.0;
   for (int i = 0; i < nt; i++) { v(i) = 1.0; }
   op.DashpotForce(v, f);
   real_t fx = 0.0, fy = 0.0;
   for (int i = 0; i < nt; i++) { fx += f(i); fy += f(nt + i); }
   TEST_REL(fx, p.density * p.Cp() * face, 1e-12, "dashpot: sum C_xx = rho c_p |face|");
   TEST_ASSERT(std::abs(fy) <= 1e-12 * fx, "dashpot: no xy coupling on a flat x face");
   v = 0.0;
   for (int i = 0; i < nt; i++) { v(nt + i) = 1.0; }
   op.DashpotForce(v, f);
   fy = 0.0;
   for (int i = 0; i < nt; i++) { fy += f(nt + i); }
   TEST_REL(fy, p.density * p.Cs() * face, 1e-12, "dashpot: sum C_yy = rho c_s |face|");

   // Pressure: -P n |face| with n = +x on the right face.
   const Vector &bv = op.PressureVector();
   real_t bx = 0.0, by = 0.0;
   for (int i = 0; i < nt; i++) { bx += bv(i); by += bv(nt + i); }
   TEST_REL(bx, -1e6 * face, 1e-12, "pressure vector = -P n |face| (outward n)");
   TEST_ASSERT(std::abs(by) <= 1e-12 * std::abs(bx), "pressure vector has no tangential part");

   // Discrete energy of the central difference (no dashpots, zeta = 0).
   ExplicitElasticityParams p0 = Params(max_attr);
   p0.pressure_marker[right - 1] = 1;
   ExplicitElasticityOperator op0(V, S, p0, pulse, fixed);
   CentralDifference cd(op0);
   CentralDifferenceState st;
   st.u.SetSize(n); st.u = 0.0;
   st.v_half.SetSize(n); st.v_half.Randomize(3);
   st.F_prev.SetSize(n); st.F_prev = 0.0;
   const real_t dt = 0.4 * CriticalStep(std::sqrt(op0.MaxEigenvalue(60)), 0.0);
   CentralDifferenceStepInfo info;
   // First step from u^0 = 0 with v^0 (v_half holds v^0); H measured from step 1 on.
   cd.Step(st, D0, 0.0, 0.0, dt, info);
   auto H = [&]() { return op0.KineticEnergy(st.v_half) + 0.5 * (st.u * st.F_prev); };
   const real_t H1 = H();
   real_t worst = 0.0;
   for (int k = 0; k < 2000; k++)
   {
      cd.Step(st, D0, 0.0, 0.0, dt, info);
      worst = std::max(worst, std::abs(H() - H1));
   }
   std::cout << "  max |H - H1| / H1 over 2000 steps: " << worst / H1 << "\n";
   TEST_ASSERT(worst <= 1e-10 * H1, "central difference conserves H over 2000 steps");

   // Mass-proportional damping on a rigid translation (F_int = 0, no dashpots).
   {
      const real_t alpha = 2e5, dtm = 2e-7;   // alpha dt = 0.04
      CentralDifferenceState sm;
      sm.u.SetSize(n); sm.u = 0.0;
      sm.v_half.SetSize(n); sm.v_half = 0.0;
      for (int i = 0; i < nt; i++) { sm.v_half(i) = 1.0; sm.v_half(nt + i) = -0.5; }
      sm.F_prev.SetSize(n); sm.F_prev = 0.0;
      const real_t ke0 = op0.KineticEnergy(sm.v_half);
      const int N = 250;
      real_t Wm = 0.0, worst_w = 0.0;
      Vector fm_prev(n), u_prev(sm.u);
      fm_prev = 0.0;
      for (int k = 0; k < N; k++)
      {
         cd.Step(sm, D0, 0.0, 0.0, dtm, info, alpha);
         // trapezoid in u from u^{k-1} to u^k (info at t^k), as in pff_driver
         Vector uk(sm.u);
         uk.Add(-dtm, sm.v_half);
         for (int i = 0; i < n; i++) { Wm += 0.5 * (fm_prev(i) + info.f_mass(i)) * (uk(i) - u_prev(i)); }
         u_prev = uk;
         fm_prev = info.f_mass;
         worst_w = std::max(worst_w, std::abs(op0.KineticEnergy(info.v_n) + Wm - ke0) / ke0);
      }
      const real_t a = alpha * dtm;
      const real_t expect = (1.0 - 0.5 * a) * std::pow((1.0 - 0.5 * a) / (1.0 + 0.5 * a), N - 1);
      TEST_REL(sm.v_half(0), expect, 1e-12, "mass damping: discrete decay factor (rigid x)");
      TEST_REL(sm.v_half(nt + 1), -0.5 * expect, 1e-12, "mass damping: same factor in y");
      // v^{N-1/2} at t = (N - 1/2) dt vs exp(-alpha t): O((alpha dt)^2) per unit alpha t
      TEST_REL(sm.v_half(0), std::exp(-alpha * (N - 0.5) * dtm), 2e-3,
               "mass damping: close to exp(-alpha t)");
      std::cout << "  mass damping: max |KE + W_mass - KE0| / KE0 = " << worst_w << "\n";
      TEST_ASSERT(worst_w <= 1e-3, "mass damping: KE + mass_damping_work = KE0");
   }
}

static void CheckSpectralGrad(Element::Type type, const char *label)
{
   Mesh m = Mesh::MakeCartesian2D(2, 2, type, true, 1.0, 1.0);
   H1_FECollection h1(1, 2);
   FiniteElementSpace V(&m, &h1, 2, Ordering::byNODES), S(&m, &h1);
   GridFunction dg(&S);
   FunctionCoefficient dfun([](const Vector &X) { return 0.3 + 0.4 * X(0) * X(1); });
   dg.ProjectCoefficient(dfun);
   ExplicitElasticityParams p = Params(1);
   SplitElasticityIntegrator integ(p.K, p.G, Decomposition::Spectral, p.deg, dg);
   real_t worst = 0.0;
   for (int e = 0; e < m.GetNE(); e++)
   {
      const FiniteElement &fe = *V.GetFE(e);
      ElementTransformation &T = *V.GetElementTransformation(e);
      const int n = fe.GetDof() * 2;
      Vector u(n), up(n), um(n), fp, fm;
      u.Randomize(10 + e);
      u -= 0.5;
      u *= 1e-3;
      DenseMatrix J;
      integ.AssembleElementGrad(fe, T, u, J);
      const real_t h = 1e-9;
      for (int k = 0; k < n; k++)
      {
         up = u; um = u; up(k) += h; um(k) -= h;
         integ.AssembleElementVector(fe, T, up, fp);
         integ.AssembleElementVector(fe, T, um, fm);
         for (int i = 0; i < n; i++)
         {
            const real_t fd = (fp(i) - fm(i)) / (2.0 * h);
            worst = std::max(worst, std::abs(fd - J(i, k)) / (J.MaxMaxNorm() + 1e-30));
         }
      }
   }
   std::cout << "\n[" << label << "] spectral element tangent vs FD: " << worst << "\n";
   TEST_ASSERT(worst < 1e-6, "spectral AssembleElementGrad = FD of AssembleElementVector");
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   Mesh tri = Mesh::MakeCartesian2D(6, 5, Element::TRIANGLE, true, 1.0, 1.0);
   Check(tri, "triangles");
   Mesh tet = Mesh::MakeCartesian3D(3, 3, 3, Element::TETRAHEDRON, 1.0, 1.0, 1.0);
   Check(tet, "tetrahedra");
   CheckSpectralGrad(Element::QUADRILATERAL, "quads");
   CheckSpectralGrad(Element::TRIANGLE, "triangles");
   return TestSummary("test_explicit_operator");
}

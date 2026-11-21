// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_damage_step.cpp — the pieces of the explicit driver's damage step
// (drivers/pff_driver.cpp, phase C2), np = 1:
//   - QsDamageOperator::Residual(x, D) = A(D) x - b(D) (the matrix-free
//     residual of the skip test equals the assembled one), AT1 and AT2;
//   - closed form at d = 0 and uniform D (ADPFFSource, farms@cdms
//     pulsepower/spe_hftc/code_2d/fracture.i, p = 2, eta = 1e-6):
//       sum_i r_i = (Gc/(c0 l) alpha'(0) - 2 (1 - eta) D) |Omega|,
//     AT1 (alpha' = 1): r > 0 (no damage, the step is skipped) iff
//     D < D_c = Gc/(2 c0 l (1-eta)); AT2 (alpha'(0) = 0): r < 0 for any D > 0;
//   - node-interior seeding (PLAN.md parity details: nodes on the box
//     boundary end at 0, interior nodes at the box value).

#include "mfem.hpp"
#include "../test_framework.hpp"
#include "../../domain/qs_damage_operator.hpp"
#include "../../phasefield/initial_damage.hpp"

using namespace mfem;
using namespace mfem::pff;

static PhaseFieldModel Model(CrackModel m)
{
   PhaseFieldModel pf;
   pf.crack.model = m;
   pf.Gc = 2720.57;   // spe_hftc: 8 l ft^2/(3E), ft = 69 MPa, l = 0.6 mm, E = 2.8 GPa
   pf.l = 6e-4;
   pf.deg.p = 2.0;
   pf.deg.eta = 1e-6;
   return pf;
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   Mesh serial = Mesh::MakeCartesian2D(8, 8, Element::TRIANGLE, true, 1e-2, 1e-2);
   ParMesh mesh(MPI_COMM_WORLD, serial);
   H1_FECollection h1(1, 2);
   L2_FECollection l2(0, 2);
   ParFiniteElementSpace S(&mesh, &h1), L(&mesh, &l2);
   const IntegrationRule &ir = IntRules.Get(Geometry::TRIANGLE, 2);
   const int n = S.GetTrueVSize();

   for (CrackModel m : {CrackModel::AT1, CrackModel::AT2})
   {
      const PhaseFieldModel pf = Model(m);
      const char *name = m == CrackModel::AT1 ? "AT1" : "AT2";
      QsDamageOperator op(S, pf, &ir);
      ParGridFunction D(&L);
      FunctionCoefficient Df([](const Vector &x) { return 1e6 * (0.5 + 80.0 * x(0) + 30.0 * x(1)); });
      D.ProjectCoefficient(Df);
      Vector x(n), r(n), Ax(n);
      x.Randomize(7);
      op.Assemble(D);
      op.A().Mult(x, Ax);
      Ax -= op.b();
      op.Residual(x, D, r);
      r -= Ax;
      std::cout << name << ": |Residual - (A x - b)|_inf / |A x - b|_inf = "
                << r.Normlinf() / Ax.Normlinf() << "\n";
      TEST_ASSERT(r.Normlinf() <= 1e-12 * Ax.Normlinf(),
                  "matrix-free damage residual = A x - b");

      // Closed form at d = 0 and uniform D.
      const real_t area = 1e-4;
      const real_t Dc = pf.Gc / (2.0 * pf.crack.c0() * pf.l * (1.0 - pf.deg.eta));
      for (real_t f : {0.5, 2.0})
      {
         D = f * Dc;
         Vector zero(n);
         zero = 0.0;
         op.Residual(zero, D, r);
         const bool at1 = m == CrackModel::AT1;
         const real_t expect = ((at1 ? pf.Gc / (pf.crack.c0() * pf.l) : 0.0)
                                - 2.0 * (1.0 - pf.deg.eta) * f * Dc) * area;
         TEST_REL(r.Sum(), expect, 1e-12,
                  "sum r = (Gc/(c0 l) alpha'(0) - 2(1-eta) D)|Omega| at d = 0");
         TEST_ASSERT((at1 && f < 1.0) ? r.Min() > 0.0 : r.Max() < 0.0,
                     at1 ? "AT1: r > 0 everywhere below the onset D_c, < 0 above"
                     : "AT2: r < 0 everywhere for D > 0 (no threshold)");
      }
   }

   // Seeding: box [0, 5 mm] x [0, 5 mm] on the 8 x 8 grid of 1.25 mm cells.
   std::vector<DamageBox> boxes(1);
   boxes[0].lo = {{-1.0, -1.0, 0.0}};
   boxes[0].hi = {{5e-3, 5e-3, 0.0}};
   boxes[0].value = 0.9;
   Vector seed;
   SeedInitialDamage(S, boxes, seed);
   ParGridFunction xg(&S), yg(&S);
   FunctionCoefficient xf([](const Vector &p) { return p(0); }), yf([](const Vector &p) { return p(1); });
   xg.ProjectCoefficient(xf);
   yg.ProjectCoefficient(yf);
   Vector X(n), Y(n);
   xg.GetTrueDofs(X);
   yg.GetTrueDofs(Y);
   int wrong = 0, seeded = 0;
   for (int i = 0; i < n; i++)
   {
      const bool interior = X(i) < 5e-3 - 1e-9 && Y(i) < 5e-3 - 1e-9;
      if (seed(i) != (interior ? 0.9 : 0.0)) { wrong++; }
      if (seed(i) > 0.0) { seeded++; }
   }
   std::cout << "seeded nodes " << seeded << " (expected 16)\n";
   TEST_ASSERT(wrong == 0 && seeded == 16,
               "seed 0.9 on box-interior nodes, 0 on the box boundary and outside");
   return TestSummary("test_damage_step");
}

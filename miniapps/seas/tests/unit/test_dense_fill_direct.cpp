// Direct fill of the dense stiffness blocks (DenseBlockStiffness::InitDirect,
// Chunhui Zhao 2026-11-18): the blocks filled straight from the element and
// face matrices (ElasticityDomainOperator::VisitStiffnessMatrices) must be the
// blocks read from the assembled K and K_eta, and give the same products.
//
//   mpirun -np N ./seas_test_dense_fill_direct        (N = 1, 2, 3)
//
// Mesh: 4 x 3 x 3 boxes cut into tets (unit cube scaled to km), one Dirichlet
// face (its boundary-face matrices enter K), heterogeneous lambda and mu, a
// second operator with other coefficients in the place of K_eta; orders 1 to 4.

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../domain/dense_block_stiffness.hpp"

#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

using namespace mfem;
using namespace mfem::seas;

namespace
{
int g_rank = 0, g_fail = 0, g_pass = 0;

void Check(bool ok, const std::string &msg)
{
   if (g_rank != 0) { return; }
   (ok ? g_pass : g_fail)++;
   std::cout << (ok ? "  PASSED  " : "  FAILED  ") << msg << "\n";
}

std::string Sci(double x)
{
   std::ostringstream o;
   o << std::scientific << std::setprecision(3) << x;
   return o.str();
}

real_t MaxAll(real_t v)
{
   MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, MPI_COMM_WORLD);
   return v;
}
real_t MinAll(real_t v)
{
   MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPITypeMap<real_t>::mpi_type, MPI_MIN, MPI_COMM_WORLD);
   return v;
}
} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   int pmax = 4;
   for (int i = 1; i + 1 < argc; i++) { if (std::string(argv[i]) == "--max-order") { pmax = std::atoi(argv[i + 1]); } }

   for (int p = 1; p <= pmax; p++)
   {
      Mesh smesh = Mesh::MakeCartesian3D(4, 3, 3, Element::TETRAHEDRON, 4e3, 3e3, 3e3);
      ParMesh pm(MPI_COMM_WORLD, smesh);
      if (g_rank == 0)
      {
         std::cout << "order " << p << ", " << smesh.GetNE() << " tets, " << Mpi::WorldSize() << " rank(s)\n";
      }
      BoundaryConfig bc;
      bc.dirichlet_attrs = {1};
      bc.natural_attrs = {6};
      bc.absorbing_attrs = {2, 3, 4, 5};
      bc.fault_attr = 7;   // none in this mesh
      bc.default_dirichlet_func = [](const Vector &, real_t, Vector &u) { u.SetSize(3); u = 0.0; };
      DomainConfig cfg;
      FunctionCoefficient lam([](const Vector &x) { return 3.0e10 * (1.0 + 0.3 * std::sin(x(0) / 900.0) + 1e-4 * x(2)); });
      FunctionCoefficient mu([](const Vector &x) { return 2.5e10 * (1.0 + 0.2 * std::cos(x(1) / 700.0) + 5e-5 * x(0)); });
      FunctionCoefficient elam([](const Vector &x) { return 4.0e6 * (1.0 + 0.4 * std::cos(x(2) / 500.0)); });
      FunctionCoefficient emu([](const Vector &x) { return 3.0e6 * (1.0 + 0.1 * std::sin(x(0) / 300.0 + x(1) / 800.0)); });
      ElasticityDomainOperator<ParMesh> elast(pm, p, lam, mu, 3.0e10, 2.5e10, 0.0, 0.0, 0.0, bc,
                                              DGMethod::IP, SolverType::OperatorOnly, cfg);
      ElasticityDomainOperator<ParMesh> kv(pm, p, elam, emu, 4.0e6, 3.0e6, 0.0, 0.0, 0.0, bc,
                                           DGMethod::IP, SolverType::OperatorOnly, cfg);
      auto *pfes = dynamic_cast<ParFiniteElementSpace *>(&elast.GetFESpace());
      const int n = pfes->GetVSize();

      using Sink = DenseBlockStiffness::DirectSink;
      std::function<bool(Sink &)> vk = [&elast](Sink &s) { return elast.VisitStiffnessMatrices(s); };
      std::function<bool(Sink &)> ve = [&kv](Sink &s) { return kv.VisitStiffnessMatrices(s); };
      // direct fill FIRST, so it cannot lean on anything the assembly set up
      DenseBlockStiffness e1, e2;
      int okd[2] = {e1.InitDirect(*pfes, vk, std::function<bool(Sink &)>()) ? 1 : 0,
                    e2.InitDirect(*pfes, vk, ve) ? 1 : 0};
      MPI_Allreduce(MPI_IN_PLACE, okd, 2, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      Check(okd[0] == 1, "p" + std::to_string(p) + ": InitDirect fills the blocks of K");
      Check(okd[1] == 1, "p" + std::to_string(p) + ": InitDirect fills the blocks of K and K_eta");

      DenseBlockStiffness d1, d2;
      int ok[2] = {d1.Init(*pfes, elast.ParStiffness(), nullptr) ? 1 : 0,
                   d2.Init(*pfes, elast.ParStiffness(), &kv.ParStiffness()) ? 1 : 0};
      MPI_Allreduce(MPI_IN_PLACE, ok, 2, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      Check(ok[0] == 1 && ok[1] == 1, "p" + std::to_string(p) + ": Init reads the assembled K and K_eta");
      if (!(okd[0] && okd[1] && ok[0] && ok[1])) { continue; }

      const real_t lay = MinAll(std::min(e1.MaxBlockDifference(d1), e2.MaxBlockDifference(d2)));
      Check(lay >= 0.0, "p" + std::to_string(p) + ": same block layout as the assembled fill");
      const real_t b1 = MaxAll(e1.MaxBlockDifference(d1)), b2 = MaxAll(e2.MaxBlockDifference(d2));
      Check(b1 == 0.0, "p" + std::to_string(p) + ": blocks of K bitwise equal to the assembled ones (max diff " + Sci(b1) + ")");
      Check(b2 == 0.0, "p" + std::to_string(p) + ": blocks of K and K_eta bitwise equal (max diff " + Sci(b2) + ")");
      const real_t a0 = MaxAll(e2.Asymmetry(0)) / MaxAll(e2.MaxEntry(0));
      const real_t a1 = MaxAll(e2.Asymmetry(1)) / MaxAll(e2.MaxEntry(1));
      Check(a0 < 1e-12 && a1 < 1e-12, "p" + std::to_string(p) + ": element and face matrices symmetric (rel "
            + Sci(a0) + ", " + Sci(a1) + ")");
      const real_t m0 = MaxAll(std::abs(e2.MaxEntry(0) - d2.MaxEntry(0))), m1 = MaxAll(std::abs(e2.MaxEntry(1) - d2.MaxEntry(1)));
      Check(m0 == 0.0 && m1 == 0.0, "p" + std::to_string(p) + ": same largest entries");

      Vector u(n), v(n), y(n), z(n), Ku(n), Ev(n);
      std::mt19937_64 rng(4321 + g_rank + 17 * p);
      std::uniform_real_distribution<real_t> U(-1.0, 1.0);
      for (int i = 0; i < n; i++) { u(i) = 1e-5 * U(rng); v(i) = 1e-5 * U(rng); }
      d1.Mult(u, y);
      e1.Mult(u, z);
      z -= y;
      const real_t r1 = MaxAll(z.Normlinf());
      Check(r1 == 0.0, "p" + std::to_string(p) + ": K u bitwise equal, halo through the space's exchange (max diff " + Sci(r1) + ")");
      d2.Mult(u, v, y);
      e2.Mult(u, v, z);
      z -= y;
      const real_t r2 = MaxAll(z.Normlinf());
      Check(r2 == 0.0, "p" + std::to_string(p) + ": K u + K_eta v bitwise equal (max diff " + Sci(r2) + ")");
      // against the assembled operators themselves
      elast.ApplyStiffness(u, Ku);
      kv.ApplyStiffness(v, Ev);
      Ku += Ev;
      e2.Mult(u, v, z);
      z -= Ku;
      const real_t r3 = MaxAll(z.Normlinf()) / MaxAll(Ku.Normlinf());
      Check(r3 < 1e-12, "p" + std::to_string(p) + ": K u + K_eta v equals the two ApplyStiffness calls (rel diff " + Sci(r3) + ")");
   }
   if (g_rank == 0) { std::cout << g_pass << " passed, " << g_fail << " failed\n"; }
   int fail = g_fail;
   MPI_Bcast(&fail, 1, MPI_INT, 0, MPI_COMM_WORLD);
   return fail ? 1 : 0;
}

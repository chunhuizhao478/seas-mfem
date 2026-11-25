// test_qd_pcg.cpp — UCERF 3-D plan rev5, Q1.1: the PCG + BoomerAMG back end
// of the QD stage solve (solver/qd_linear_solver.hpp, SolverType::QD_PCG).
//
//   mpirun -np N ./seas_test_qd_pcg [--mesh FILE] [--order P] [--nrhs K]
//          [--tol T] [--precs "amg aux_h1"] [--strength S]
//          [--smoother default|l1_jacobi]
//   mpirun -np 1 ./seas_test_qd_pcg --refuse N        (must abort, N = 1..11)
//
// U-Q1 (default run): on the T-plane mesh (ucerf3d/mesh/T-plane.msh; walls 5
//   Dirichlet with the BP5 loading, top 1 natural, bottom 2 and the creep
//   faces 6 left unlisted, fault 3), K random stage right-hand sides
//   b = B delta + D g(t), half of them with t = 0 (slip only), solved by
//   MUMPS and by PCG with each preconditioner of --precs (default both,
//   "amg" and "aux_h1").  Pass: displacement and fault traction (shear and
//   normal) within 100 x tol of MUMPS (relative, max norm) (bound set 11-25:
//   at the same true residual amg misses MUMPS by up to 31 x tol in
//   displacement and aux_h1 by up to 56 x tol in fault traction), and every PCG
//   solve converged.  MUMPS runs first and is freed before the AMG set-up,
//   so the two factorizations are never held together.  The iteration
//   count and seconds of every solve are printed.
// U-Q2 (--refuse N): the deck keys and the operator refuse a bad set-up:
//   1 pcg_rel_tol with solver = "mumps"      2 pcg_rel_tol = 0
//   3 pcg_max_it = 0                         4 pcg_amg_strength = 1
//   5 pcg_smoother = "jacobi"                6 solver = "cg"
//   7 solver = "pcg" with [cycle.qd.green]   8 a QD_PCG operator with BR2
//     (K not symmetric positive definite)
//   9 a PCG solve that cannot converge in pcg_max_it = 2 iterations
//  10 pcg_preconditioner = "ilu"           11 pcg_amg_strength with "aux_h1"
//   The make target checks each exit code and message.
//
// Exit code 0 when every check passes.

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../constitutive/linear_elastic.hpp"
#include "../../config/bp5_params.hpp"
#include "../../solver/qd_linear_solver.hpp"
#include "../../cycle/cycle_config.hpp"

#include <unistd.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace
{

MPI_Comm comm = MPI_COMM_WORLD;
int myrank = 0;
int failures = 0;

real_t GlobalMaxAbs(const Vector &v)
{
   real_t a = 0.0;
   for (int i = 0; i < v.Size(); i++) { a = std::max(a, std::abs(v(i))); }
   MPI_Allreduce(MPI_IN_PLACE, &a, 1, MFEM_MPI_REAL_T, MPI_MAX, comm);
   return a;
}

real_t RelDiff(const Vector &ref, const Vector &x)
{
   Vector d(ref);
   d -= x;
   const real_t r = GlobalMaxAbs(ref);
   return r > 0.0 ? GlobalMaxAbs(d) / r : GlobalMaxAbs(d);
}

void Check(const std::string &name, real_t err, real_t bound)
{
   const bool ok = (err <= bound);
   if (!ok) { failures++; }
   if (myrank == 0)
   {
      std::printf("  %-52s %10.3e  <= %10.3e  %s\n", name.c_str(), err, bound,
                  ok ? "PASS" : "FAIL");
   }
}

std::unique_ptr<ElasticityDomainOperator<ParMesh>> MakeDomain(
   ParMesh &pmesh, int order, const LinearElastic &model, const BP5Params &params,
   DGMethod method, SolverType st)
{
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.dirichlet_attrs = {5};
   bc.natural_attrs = {1};
   bc.default_dirichlet_func = MakeBP5DirichletFunc(params.Vp);
   return std::make_unique<ElasticityDomainOperator<ParMesh>>(
      pmesh, order, model, params.Vp, params.Wf, params.lf, bc, method, st);
}

// ---------------------------------------------------------------- U-Q2
const char *kDeck = "config/fd_cycle/bp5_fd_first_event.toml";

std::string ReadFile(const std::string &path)
{
   std::ifstream f(path);
   MFEM_VERIFY(f.good(), "cannot read " << path);
   std::stringstream s;
   s << f.rdbuf();
   return s.str();
}

/// The deck with `lines` added at the top of its [cycle.qd] table (or a new
/// table when it has none).
std::string WithQD(const std::string &deck, const std::string &lines)
{
   const std::string h = "[cycle.qd]\n";
   const size_t p = deck.find(h);
   if (p == std::string::npos) { return deck + "\n" + h + lines + "\n"; }
   std::string out = deck;
   out.insert(p + h.size(), lines + "\n");
   return out;
}

/// Remove `key = ...` lines of the [cycle.qd] table (so a test can set the key).
std::string DropQDKey(const std::string &deck, const std::string &key)
{
   std::istringstream in(deck);
   std::ostringstream out;
   std::string line, table;
   while (std::getline(in, line))
   {
      if (!line.empty() && line[0] == '[') { table = line; }
      const size_t k = line.find_first_not_of(" \t");
      if (table == "[cycle.qd]" && k != std::string::npos
          && line.compare(k, key.size(), key) == 0
          && line.find('=', k) != std::string::npos
          && line.substr(k, line.find('=', k) - k).find_first_not_of(" \t", key.size())
             == std::string::npos)
      {
         continue;
      }
      out << line << "\n";
   }
   return out.str();
}

cycle::CycleConfig LoadText(const std::string &tag, const std::string &text)
{
   const std::string dir = "/tmp";
   const std::string path = dir + "/seas_test_qd_pcg_" + tag + "_"
                            + std::to_string(getpid()) + ".toml";
   {
      std::ofstream f(path);
      f << text;
   }
   cycle::CycleConfig c = cycle::CycleConfig::Load(path);
   std::remove(path.c_str());
   return c;
}

int Refuse(int n, const std::string &mesh_file)
{
   if ((n >= 1 && n <= 7) || n == 10 || n == 11)
   {
      const std::string deck = DropQDKey(ReadFile(kDeck), "solver");
      std::string q;
      switch (n)
      {
         case 1: q = "solver = \"mumps\"\npcg_rel_tol = 1e-8"; break;
         case 2: q = "solver = \"pcg\"\npcg_rel_tol = 0.0"; break;
         case 3: q = "solver = \"pcg\"\npcg_max_it = 0"; break;
         case 4: q = "solver = \"pcg\"\npcg_amg_strength = 1.0"; break;
         case 5: q = "solver = \"pcg\"\npcg_smoother = \"jacobi\""; break;
         case 6: q = "solver = \"cg\""; break;
         case 7: q = "solver = \"pcg\"\noperators = \"matrix\""; break;
         case 10: q = "solver = \"pcg\"\npcg_preconditioner = \"ilu\""; break;
         case 11: q = "solver = \"pcg\"\npcg_preconditioner = \"aux_h1\"\n"
                      "pcg_amg_strength = 0.5"; break;
      }
      std::string text = WithQD(n == 7 ? DropQDKey(deck, "operators") : deck, q);
      if (n == 7) { text += "\n[cycle.qd.green]\n"; }
      if (myrank == 0) { std::cout << "--refuse " << n << ": deck [cycle.qd] " << q << "\n"; }
      LoadText("refuse" + std::to_string(n), text);   // must abort
      if (myrank == 0) { std::cout << "--refuse " << n << ": the deck was accepted\n"; }
      return 0;
   }
   MFEM_VERIFY(n == 8 || n == 9, "--refuse N: N = 1..11");
   Mesh smesh(mesh_file.c_str(), 1, 1);
   ParMesh pmesh(comm, smesh);
   smesh.Clear();
   BP5Params params;
   LinearElastic model(params.lambda(), params.mu());
   auto dom = MakeDomain(pmesh, 1, model, params,
                         n == 8 ? DGMethod::BR2 : DGMethod::IP, SolverType::QD_PCG);
   QDPCGConfig pc;
   if (n == 9) { pc.max_it = 2; }
   dom->SetQDPCGConfig(pc);
   ParGridFunction u(dynamic_cast<ParFiniteElementSpace *>(&dom->GetFESpace()));
   Vector slip(2 * dom->GetNumFaultDOFs());
   slip = 1.0;
   if (myrank == 0) { std::cout << "--refuse " << n << ": a QD_PCG solve\n"; }
   dom->Solve(0.0, slip, u);   // must abort
   if (myrank == 0) { std::cout << "--refuse " << n << ": the solve was accepted\n"; }
   return 0;
}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   myrank = Mpi::WorldRank();

   std::string mesh_file = "ucerf3d/mesh/T-plane.msh";
   int order = 1, nrhs = 10, refuse = 0, seed = 2026;
   real_t tol = 1e-10, strength = 0.8;
   std::string smoother = "default", precs = "amg aux_h1", reference = "mumps";
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "Gmsh mesh (T-plane).");
   args.AddOption(&order, "-o", "--order", "Polynomial order.");
   args.AddOption(&nrhs, "-n", "--nrhs", "Number of right-hand sides.");
   args.AddOption(&tol, "-t", "--tol", "PCG relative tolerance.");
   args.AddOption(&precs, "-p", "--precs", "Preconditioners to test.");
   args.AddOption(&strength, "-s", "--strength", "AMG strength threshold.");
   args.AddOption(&smoother, "-x", "--smoother", "default or l1_jacobi.");
   args.AddOption(&seed, "-r", "--seed", "Random seed.");
   args.AddOption(&refuse, "-R", "--refuse", "Refusal case 1..9 (must abort).");
   args.AddOption(&reference, "-ref", "--reference",
                  "mumps (U-Q1) or none (timing only, for meshes whose MUMPS "
                  "factor does not fit: checks |b - K u|/|b| <= 10 x tol).");
   args.Parse();
   if (!args.Good())
   {
      if (myrank == 0) { args.PrintUsage(std::cout); }
      return 1;
   }
   if (refuse) { return Refuse(refuse, mesh_file); }
   MFEM_VERIFY(reference == "mumps" || reference == "none",
               "--reference must be mumps or none");
   const bool ref = (reference == "mumps");

   Mesh smesh(mesh_file.c_str(), 1, 1);
   ParMesh pmesh(comm, smesh);
   smesh.Clear();
   BP5Params params;
   LinearElastic model(params.lambda(), params.mu());
   auto dom = MakeDomain(pmesh, order, model, params, DGMethod::IP,
                         SolverType::OperatorOnly);
   auto *dg = dynamic_cast<ParFiniteElementSpace *>(&dom->GetFESpace());
   dom->BuildQDStageOperators(true);
   HypreParMatrix &K = dom->ParStiffness();
   const int nown = dom->GetNumOwnedFaultDOFs();
   const long long N = K.GetGlobalNumRows();
   if (myrank == 0)
   {
      std::printf("%s: %s, p = %d, %d ranks, %lld DOFs, "
                  "%d right-hand sides, tol %.1e, preconditioners %s, AMG strength "
                  "%.2f, smoother %s\n",
                  ref ? "U-Q1 PCG against MUMPS" : "PCG timing (no reference)",
                  mesh_file.c_str(), order, Mpi::WorldSize(),
                  N, nrhs, tol, precs.c_str(), strength, smoother.c_str());
   }

   std::mt19937 gen(seed + myrank);
   std::uniform_real_distribution<double> U(-1.0, 1.0);
   std::vector<Vector> slip(nrhs), b(nrhs), xm(nrhs);
   std::vector<real_t> t(nrhs);
   for (int k = 0; k < nrhs; k++)
   {
      slip[k].SetSize(2 * nown);
      for (int i = 0; i < slip[k].Size(); i++) { slip[k](i) = U(gen); }
      t[k] = (k % 2) ? 3.15e9 * (k + 1) : 0.0;   // up to 1000 yr of loading
      dom->QDStageRHS(t[k], slip[k], b[k]);
   }

   // ---- MUMPS, then freed
   if (ref)
   {
      QDLinearSolverSetup setup;
      auto mumps = MakeQDLinearSolver(QDLinearSolverKind::MUMPS, K, 0, 1e-12, &setup);
      if (myrank == 0) { std::cout << QDLinearSolverSummary(setup) << "\n"; }
      for (int k = 0; k < nrhs; k++)
      {
         xm[k].SetSize(b[k].Size());
         xm[k] = 0.0;
         mumps->Mult(b[k], xm[k]);
      }
   }

   // ---- PCG through the factory, each preconditioner
   MFEM_VERIFY(smoother == "default" || smoother == "l1_jacobi", "--smoother");
   std::istringstream ps(precs);
   std::string prec;
   while (ps >> prec)
   {
      QDPCGConfig pc;
      pc.rel_tol = tol;
      pc.preconditioner = prec;
      pc.amg_strength = strength;
      pc.l1_jacobi = (smoother == "l1_jacobi");
      QDLinearSolverSetup setup;
      auto solver = MakeQDLinearSolver(QDLinearSolverKind::PCG, K, 0, 1e-12,
                                       &setup, &pc, dg);
      if (myrank == 0) { std::cout << QDLinearSolverSummary(setup) << "\n"; }
      const auto *pcg = dynamic_cast<const QDPCGSolver *>(solver.get());
      MFEM_VERIFY(pcg, "the factory did not return a QDPCGSolver");

      real_t worst_u = 0.0, worst_T = 0.0, worst_N = 0.0, worst_r = 0.0;
      std::vector<double> times;
      for (int k = 0; k < nrhs; k++)
      {
         Vector x(b[k].Size());
         x = 0.0;
         MPI_Barrier(comm);
         const double t0 = MPI_Wtime();
         solver->Mult(b[k], x);   // aborts when it does not converge
         const double sec = MPI_Wtime() - t0;
         times.push_back(sec);
         // true residual |b - K x| / |b| (CG stops on the preconditioned norm)
         Vector r(b[k].Size());
         K.Mult(x, r);
         r -= b[k];
         const real_t rr = std::sqrt(InnerProduct(comm, r, r) / InnerProduct(comm, b[k], b[k]));
         worst_r = std::max(worst_r, rr);
         real_t eu = -1.0, eT = -1.0, eN = -1.0;
         if (ref)
         {
            Vector Tm, Tc, Nm, Nc;
            dom->QDStageTraction(xm[k], slip[k], Tm, &Nm);
            dom->QDStageTraction(x, slip[k], Tc, &Nc);
            eu = RelDiff(xm[k], x);
            eT = RelDiff(Tm, Tc);
            // normal traction relative to the shear scale (T_n may be ~0 for slip)
            Vector dN(Nm);
            dN -= Nc;
            eN = GlobalMaxAbs(dN) / std::max(GlobalMaxAbs(Nm), GlobalMaxAbs(Tm));
         }
         worst_u = std::max(worst_u, eu);
         worst_T = std::max(worst_T, eT);
         worst_N = std::max(worst_N, eN);
         if (myrank == 0)
         {
            std::printf("  %-6s rhs %2d  t = %9.3e s  iterations %4d  rel. residual %.2e (true %.2e)  "
                        "%.3f s | u %.2e  T %.2e  T_n %.2e\n", prec.c_str(), k, t[k],
                        pcg->LastIterations(), pcg->LastRelNorm(), rr, sec, eu, eT, eN);
         }
      }
      if (myrank == 0)
      {
         std::vector<double> ts(times);
         std::sort(ts.begin(), ts.end());
         std::printf("  %s: mean iterations %.1f over %lld solves, median %.3f s per "
                     "solve (rank 0)\n", prec.c_str(),
                     double(pcg->TotalIterations()) / double(pcg->NumSolves()),
                     pcg->NumSolves(), ts.empty() ? 0.0 : ts[ts.size() / 2]);
      }
      Check(prec + " true residual |b - K u|/|b|, worst", worst_r, 10.0 * tol);
      if (!ref) { continue; }
      Check("U-Q1 " + prec + " displacement, max over right-hand sides", worst_u,
            100.0 * tol);
      Check("U-Q1 " + prec + " fault shear traction", worst_T, 100.0 * tol);
      Check("U-Q1 " + prec + " fault normal traction (relative to shear)", worst_N,
            100.0 * tol);
   }

   if (myrank == 0)
   {
      std::cout << (failures ? "SOME CHECKS FAILED" : "ALL CHECKS PASSED")
                << " (" << failures << " failures)\n";
   }
   return failures ? 1 : 0;
}

// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// GPU plan G0: ElasticityDomainOperator::AssembleMatrixOnly() (seas_bp5_full
// --exit-after-assembly) must describe the same matrix the solver path
// factorizes.  Two operators on one ParMesh write matrix statistics: one
// through the first Solve() (assembly + MUMPS factorization), one through
// AssembleMatrixOnly() (assembly only).  Every count, norm and per-rank
// figure must agree; the matrix-only file skips the symmetry check (-1), and
// its per-rank nonzero counts must add up to the global count.
//
// Usage:  mpirun -np 1 ./seas_test_g0_matrix_only
//         mpirun -np 2 ./seas_test_g0_matrix_only

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../config/bp5_params.hpp"
#include "test_macros.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

using namespace mfem;
using namespace mfem::seas;

namespace
{

// Tet box [-Lx, Lx] x [-Ly, Ly] x [-Lz, 0], top/bottom natural (1), sides
// Dirichlet (5), fault elements (3) on the interior faces at y = 0.
Mesh MakeBoxWithFault(int nx, int ny, int nz, real_t Lx, real_t Ly, real_t Lz)
{
   Mesh mesh = Mesh::MakeCartesian3D(2 * nx, 2 * ny, nz, Element::TETRAHEDRON,
                                     2.0 * Lx, 2.0 * Ly, Lz);
   for (int i = 0; i < mesh.GetNV(); i++)
   {
      real_t *v = mesh.GetVertex(i);
      v[0] -= Lx; v[1] -= Ly; v[2] -= Lz;
   }
   for (int be = 0; be < mesh.GetNBE(); be++)
   {
      ElementTransformation *T = mesh.GetBdrElementTransformation(be);
      const IntegrationPoint &ip = Geometries.GetCenter(T->GetGeometryType());
      T->SetIntPoint(&ip);
      Vector c(3);
      T->Transform(ip, c);
      const bool top_bottom = std::abs(c(2)) < 1e-6 || std::abs(c(2) + Lz) < 1e-6;
      mesh.SetBdrAttribute(be, top_bottom ? 1 : 5);
   }
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      auto *FTr = mesh.GetInteriorFaceTransformations(f);
      if (!FTr) { continue; }
      const IntegrationPoint &ip = Geometries.GetCenter(FTr->GetGeometryType());
      FTr->Face->SetIntPoint(&ip);
      Vector c(3);
      FTr->Face->Transform(ip, c);
      if (std::abs(c(1)) > 1e-6) { continue; }
      Array<int> v;
      mesh.GetFaceVertices(f, v);
      mesh.AddBdrTriangle(v[0], v[1], v[2], 3);
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// The numbers in a matrix-stats JSON file, keyed "a" or "a.b" (one level of
// nesting, as WriteMatrixStatsJson writes).  Booleans become 0/1; strings
// are skipped.
std::map<std::string, double> ReadStats(const std::string &path)
{
   std::map<std::string, double> out;
   std::ifstream f(path);
   std::stringstream ss;
   ss << f.rdbuf();
   const std::string s = ss.str();
   std::string outer;
   size_t i = 0;
   int depth = 0;
   while (i < s.size())
   {
      const char c = s[i];
      if (c == '{') { depth++; i++; continue; }
      if (c == '}') { depth--; if (depth <= 1) { outer.clear(); } i++; continue; }
      if (c != '"') { i++; continue; }
      const size_t e = s.find('"', i + 1);
      const std::string key = s.substr(i + 1, e - i - 1);
      i = e + 1;
      while (i < s.size() && (s[i] == ' ' || s[i] == ':')) { i++; }
      const std::string full = (depth >= 2 && !outer.empty()) ? outer + "." + key : key;
      if (s[i] == '{') { outer = key; continue; }
      if (s[i] == '"')
      {
         i = s.find('"', i + 1) + 1;  // string value
         continue;
      }
      if (s.compare(i, 4, "true") == 0) { out[full] = 1.0; i += 4; continue; }
      if (s.compare(i, 5, "false") == 0) { out[full] = 0.0; i += 5; continue; }
      size_t used = 0;
      try { out[full] = std::stod(s.substr(i), &used); }
      catch (...) { used = 1; }
      i += used;
   }
   return out;
}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   const int rank = Mpi::WorldRank();
   const int nranks = Mpi::WorldSize();
   if (rank == 0)
   {
      std::cout << "========================================\n"
                << "  G0 matrix-only assembly test (" << nranks << " ranks)\n"
                << "========================================\n";
   }

   const real_t Lx = 100.0e3, Ly = 60.0e3, Lz = 50.0e3;
   Mesh serial = MakeBoxWithFault(2, 2, 2, Lx, Ly, Lz);
   ParMesh pmesh(MPI_COMM_WORLD, serial);
   BP5Params params;

   const long long global_ne = pmesh.GetGlobalNE();  // collective
   const std::string full_path = "g0_matrix_only_test_full.json";
   const std::string only_path = "g0_matrix_only_test_only.json";

   for (int order = 1; order <= 2; order++)
   {
      if (rank == 0) { std::cout << "\n--- order " << order << " ---\n"; }
      {
         ElasticityDomainOperator<ParMesh> op(
            pmesh, order, params.lambda(), params.mu(), params.Vp, Lz,
            2.0 * Ly, DGMethod::IP, SolverType::MUMPS);
         op.SetMatrixStatsConfig(true, full_path);
         Vector slip(op.GetNumFaultDOFs() * op.NumSlipComponents());
         slip = 0.0;
         ParGridFunction u(&op.GetFESpace());
         op.Solve(0.0, slip, u);
      }
      {
         ElasticityDomainOperator<ParMesh> op(
            pmesh, order, params.lambda(), params.mu(), params.Vp, Lz,
            2.0 * Ly, DGMethod::IP, SolverType::MUMPS);
         op.SetMatrixStatsConfig(true, only_path);
         op.AssembleMatrixOnly();
      }
      MPI_Barrier(MPI_COMM_WORLD);

      if (rank == 0)
      {
         const auto full = ReadStats(full_path);
         const auto only = ReadStats(only_path);
         TEST_ASSERT(!full.empty() && !only.empty(),
                     "both statistics files written");

         // Everything except the symmetry check must be identical.
         int compared = 0, differ = 0;
         for (const auto &kv : full)
         {
            if (kv.first.rfind("symmetry", 0) == 0) { continue; }
            auto it = only.find(kv.first);
            compared++;
            if (it == only.end() || it->second != kv.second)
            {
               differ++;
               std::cout << "  differs: " << kv.first << " " << kv.second
                         << " vs " << (it == only.end() ? -999.0 : it->second)
                         << "\n";
            }
         }
         std::cout << "  " << compared << " values compared\n";
         TEST_ASSERT(compared > 20 && differ == 0,
                     "matrix-only statistics equal the solver path's");
         TEST_ASSERT(only.count("symmetry.relative_to_fnorm") &&
                     only.at("symmetry.relative_to_fnorm") == -1.0 &&
                     full.at("symmetry.relative_to_fnorm") >= 0.0,
                     "symmetry check skipped only in matrix-only mode");

         const double nnz = only.at("global_nnz");
         TEST_ASSERT(only.at("n_ranks") == nranks, "rank count recorded");
         TEST_ASSERT(std::abs(only.at("local_nnz.mean") * nranks - nnz)
                     <= 1e-9 * nnz,
                     "per-rank nonzeros add up to the global count");
         TEST_ASSERT(only.at("local_nnz.min") <= only.at("local_nnz.mean") &&
                     only.at("local_nnz.mean") <= only.at("local_nnz.max") &&
                     only.at("local_nnz.max_fits_int32") == 1.0,
                     "per-rank nonzero range is ordered and fits int32");
         const int block = (order == 1 ? 4 : 10) * 3;
         TEST_ASSERT(only.at("dg.block_size") == block,
                     "DG block size is (dofs per tet) x 3");
         TEST_ASSERT(only.at("global_n") == block * global_ne,
                     "global size is block size x elements");
         std::remove(full_path.c_str());
         std::remove(only_path.c_str());
      }
      MPI_Barrier(MPI_COMM_WORLD);
   }

   int failed = 0;
   if (rank == 0)
   {
      TEST_PRINT_RESULTS();
      failed = num_failed;
   }
   MPI_Bcast(&failed, 1, MPI_INT, 0, MPI_COMM_WORLD);
   return failed;
}

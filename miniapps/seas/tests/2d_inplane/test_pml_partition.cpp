// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/2d_inplane/test_pml_partition.cpp — cycle/pml_partition.hpp
// ([cycle.pml] partition_weight: the mesh cut that shares the band's work).
//
// On a 2-D box of triangles (40 x 40 squares of 5 km, split in two; the
// band 50 km thick on all four walls, attributes 1-4), with P ranks:
//   W1  partition_weight = 0, or [cycle.pml] off: no cut is returned (the
//       caller keeps MFEM's default cut).
//   W2  partition_weight = 4: every element has a rank in [0, P), every rank
//       has elements, and every rank holds the same cut (it is broadcast).
//   W3  the weighted work (1 per element outside, 1 + 4 inside a band) of
//       the busiest rank is within 10 % of the mean; with MFEM's default
//       cut of the same mesh it is printed for comparison (not checked).
//   W4  the band elements found = the count from the geometry (centroid
//       within 50 km of a wall).
// Usage (from miniapps/seas): mpirun -np P ./seas_test_pml_partition
// (P = 1: W1 only, a cut needs more than one rank).

#include "mfem.hpp"
#include "../../cycle/pml_partition.hpp"
#include <iostream>

using namespace mfem;
using namespace mfem::seas;

namespace
{
int g_fail = 0, g_pass = 0;
void Check(bool ok, const std::string &what, int rank)
{
   int all = ok ? 1 : 0;
   MPI_Allreduce(MPI_IN_PLACE, &all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
   if (all) { g_pass++; } else { g_fail++; }
   if (rank == 0) { std::cout << (all ? "  PASSED: " : "  FAILED: ") << what << std::endl; }
}
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   const int rank = Mpi::WorldRank(), P = Mpi::WorldSize();
   const real_t Lbox = 200e3, Lband = 50e3, w = 4.0;
   const int n = 40;
   Mesh smesh = Mesh::MakeCartesian2D(n, n, Element::TRIANGLE, true, Lbox, Lbox);
   const int ne = smesh.GetNE();

   cycle::CycleConfig cfg;
   cfg.pml.on = true;
   cfg.pml.kind = "cfs";
   cfg.pml.wall_attrs = {1, 2, 3, 4};
   cfg.pml.thickness = Lband;
   cfg.pml.power = 2;
   cfg.pml.R0 = 1e-3;
   cfg.pml.f_alpha = 0.01;
   cfg.pml.c_ref = 6000.0;
   if (rank == 0) { std::cout << "test_pml_partition: " << P << " rank(s), " << ne << " triangles\n"; }

   // W1
   cfg.pml.partition_weight = 0.0;
   Check(cycle::PlanPmlPartition(smesh, cfg, MPI_COMM_WORLD).part.empty(),
         "W1 partition_weight = 0: no cut (the default cut stays)", rank);
   cfg.pml.partition_weight = w;
   cfg.pml.on = false;
   Check(cycle::PlanPmlPartition(smesh, cfg, MPI_COMM_WORLD).part.empty(),
         "W1 [cycle.pml] off: no cut", rank);
   cfg.pml.on = true;

   if (P > 1)
   {
      const cycle::PmlPartition R = cycle::PlanPmlPartition(smesh, cfg, MPI_COMM_WORLD);
      if (rank == 0) { std::cout << "  " << R.report << "\n"; }
      // W2
      bool ok = static_cast<int>(R.part.size()) == ne;
      std::vector<long> cnt(P, 0);
      long sum = 0;
      for (int e = 0; ok && e < ne; e++)
      {
         ok = R.part[e] >= 0 && R.part[e] < P;
         if (ok) { cnt[R.part[e]]++; sum += static_cast<long>(R.part[e]) * (e % 97 + 1); }
      }
      for (int r = 0; ok && r < P; r++) { ok = cnt[r] > 0; }
      long smin = sum, smax = sum;
      MPI_Allreduce(MPI_IN_PLACE, &smin, 1, MPI_LONG, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &smax, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
      Check(ok && smin == smax, "W2 every element has a rank, every rank has elements, "
            "the same cut on every rank", rank);
      // W3, W4
      std::vector<double> work(P, 0.0), work0(P, 0.0);
      int *def = smesh.GeneratePartitioning(P);
      long inside = 0;
      Vector c(2);
      double total = 0.0;
      for (int e = 0; e < ne; e++)
      {
         smesh.GetElementCenter(e, c);
         const bool in = c(0) < Lband || c(0) > Lbox - Lband || c(1) < Lband || c(1) > Lbox - Lband;
         const double cost = in ? 1.0 + w : 1.0;
         inside += in ? 1 : 0;
         if (ok) { work[R.part[e]] += cost; }
         work0[def[e]] += cost;
         total += cost;
      }
      delete [] def;
      double wmax = 0.0, w0max = 0.0;
      for (int r = 0; r < P; r++) { wmax = std::max(wmax, work[r]); w0max = std::max(w0max, work0[r]); }
      const double imb = wmax / (total / P), imb0 = w0max / (total / P);
      std::ostringstream s;
      s << "W3 weighted work max/mean = " << imb << " (limit 1.10); MFEM's default cut: " << imb0;
      Check(ok && imb < 1.10, s.str(), rank);
      const std::string key = " " + std::to_string(inside) + " of " + std::to_string(ne) + " elements inside";
      Check(rank != 0 || R.report.find(key) != std::string::npos,
            "W4 band elements: the report has" + key + " a band (the geometry's count)", rank);
   }
   if (rank == 0)
   {
      std::cout << "test_pml_partition: " << g_pass << " PASSED, " << g_fail << " FAILED" << std::endl;
   }
   return g_fail ? 1 : 0;
}

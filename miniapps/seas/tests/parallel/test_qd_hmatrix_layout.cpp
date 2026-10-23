// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

// test_qd_hmatrix_layout.cpp — QD speed item C (qd-c-design-rev1 §5.1): the
// bytes layout of the compressed Green's function (qd_hmatrix_operator.hpp)
// gives bitwise the same product as the count layout.  No mesh: G is a
// smooth 1/r kernel on a jittered planar point cloud.
//
//   mpirun -np N ./seas_test_qd_hmatrix_layout [--nx NX] [--ny NY]
//          [--leaf L] [--eps EPS] [--file PATH]
//
//   1. Compress G in the count layout and write the H file; the product
//      agrees with the dense G to within 10·EPS (relative, max norm).
//   2. Read the file back with fake nodes of 1, 2 and 4 ranks and the real
//      nodes, in the count layout and in the bytes layout with piece
//      divisors 4 and 1000 (1000 cuts every V^T into single rows).  For
//      three slip vectors, y gathered in global order must equal the
//      product of step 1 byte for byte (memcmp), and T_D its rows.
//   3. In the bytes layout every rank holds exactly the bytes PlanLayout
//      planned for it, and divisor 1000 makes more V^T pieces than the
//      count layout.
//
// Exit code 0 when every check passes.

#include "mfem.hpp"
#include "../../domain/qd_hmatrix_operator.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace
{

int g_fail = 0;

void Check(bool ok, const std::string &what, int rank)
{
   if (!ok) { g_fail = 1; }
   if (rank == 0)
   {
      std::cout << (ok ? "  PASS  " : "  FAIL  ") << what << std::endl;
   }
}

/// Even split of nf fault DOFs over P ranks (QDGreenOperator's layout).
std::vector<long long> CountLayout(long long nf, int P)
{
   std::vector<long long> lo(P + 1);
   for (int r = 0; r <= P; r++) { lo[r] = nf * r / P; }
   return lo;
}

/// @a y (this rank's rows, layout @a lo) gathered on rank 0 in order g.
std::vector<real_t> GatherRows(MPI_Comm comm, const std::vector<real_t> &y,
                               const std::vector<long long> &lo, int ncomp)
{
   int rank, P;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &P);
   std::vector<int> cnt(P), dsp(P);
   for (int r = 0; r < P; r++)
   {
      cnt[r] = static_cast<int>(ncomp * (lo[r + 1] - lo[r]));
      dsp[r] = static_cast<int>(ncomp * lo[r]);
   }
   std::vector<real_t> all(rank == 0 ? ncomp * lo[P] : 0);
   MPI_Gatherv(y.data(), cnt[rank], MPITypeMap<real_t>::mpi_type, all.data(),
               cnt.data(), dsp.data(), MPITypeMap<real_t>::mpi_type, 0, comm);
   return all;
}

bool SameBytes(const std::vector<real_t> &a, const std::vector<real_t> &b)
{
   return a.size() == b.size() &&
          (a.empty() || std::memcmp(a.data(), b.data(),
                                    a.size() * sizeof(real_t)) == 0);
}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   MPI_Comm comm = MPI_COMM_WORLD;
   const int rank = Mpi::WorldRank(), P = Mpi::WorldSize();
   int nx = 40, ny = 30, leaf = 16;
   double eps = 1e-8;
   std::string path = "qd_hmatrix_layout_test.bin";
   for (int i = 1; i < argc; i++)
   {
      const std::string a = argv[i];
      if (a == "--nx" && i + 1 < argc) { nx = std::atoi(argv[++i]); }
      else if (a == "--ny" && i + 1 < argc) { ny = std::atoi(argv[++i]); }
      else if (a == "--leaf" && i + 1 < argc) { leaf = std::atoi(argv[++i]); }
      else if (a == "--eps" && i + 1 < argc) { eps = std::atof(argv[++i]); }
      else if (a == "--file" && i + 1 < argc) { path = argv[++i]; }
   }
   const int ncomp = 3;
   const long long nf = static_cast<long long>(nx) * ny;
   const double eta = 2.0, mem_node = 1e12;
   if (rank == 0)
   {
      std::cout << "test_qd_hmatrix_layout: " << P << " ranks, N_f " << nf
                << ", leaf " << leaf << ", eps " << eps << std::endl;
   }

   // Points: a jittered nx × ny grid in the plane z = 0.3·x (spacing 1).
   // Keys: the point, then a nearby "face centroid".
   std::vector<real_t> X(3 * nf), keys(6 * nf);
   for (long long g = 0; g < nf; g++)
   {
      const int i = static_cast<int>(g % nx), j = static_cast<int>(g / nx);
      const double x = i + 0.3 * std::sin(1.7 * g + 0.4);
      const double y = j + 0.3 * std::cos(2.3 * g + 0.1);
      X[3 * g] = x;
      X[3 * g + 1] = y;
      X[3 * g + 2] = 0.3 * x;
      for (int d = 0; d < 3; d++)
      {
         keys[6 * g + d] = X[3 * g + d];
         keys[6 * g + 3 + d] = X[3 * g + d] + 0.1 * (d + 1);
      }
   }
   // G(ncomp·g + c, 2·h + e): smooth in |x_g - x_h|, not symmetric
   auto Gval = [&](long long g, int c, long long h, int e)
   {
      double r2 = 0.25, dx = 0.0;
      for (int d = 0; d < 3; d++)
      {
         const double t = X[3 * g + d] - X[3 * h + d];
         r2 += t * t;
         if (d == c % 2) { dx = t; }
      }
      return (1.0 + 0.1 * c + 0.05 * e) / std::sqrt(r2) +
             (0.2 + 0.03 * c * e) * dx / r2;
   };
   const std::vector<long long> lo0 = CountLayout(nf, P);
   const int nb0 = static_cast<int>(lo0[rank + 1] - lo0[rank]);
   std::vector<real_t> rows(static_cast<size_t>(nb0) * ncomp * 2 * nf),
       td(static_cast<size_t>(nb0) * ncomp);
   for (int i = 0; i < nb0; i++)
   {
      const long long g = lo0[rank] + i;
      for (int c = 0; c < ncomp; c++)
      {
         real_t *row = rows.data() + (static_cast<size_t>(i) * ncomp + c) * 2
                       * nf;
         for (long long h = 0; h < nf; h++)
         {
            for (int e = 0; e < 2; e++) { row[2 * h + e] = Gval(g, c, h, e); }
         }
         td[ncomp * i + c] = std::sin(0.01 * g + c);
      }
   }
   // Three slip vectors, the same on every rank
   std::vector<std::vector<real_t>> slips(3, std::vector<real_t>(2 * nf));
   for (long long q = 0; q < 2 * nf; q++)
   {
      slips[0][q] = std::sin(0.37 * q + 0.1) + 0.5 * std::cos(0.011 * q);
      slips[1][q] = (q % 7 == 0) ? 1.0 : 0.0;
      slips[2][q] = std::exp(-1e-3 * q) * ((q % 2) ? -1.0 : 1.0);
   }

   // 1. Compress in the count layout and write the file
   unsetenv("SEAS_QD_HMAT_NODE_SIZE");
   unsetenv("SEAS_QD_HMAT_PIECE_DIV");
   std::vector<std::vector<real_t>> yref(3);
   const std::vector<real_t> tdref = GatherRows(comm, td, lo0, ncomp);
   {
      QDHMatrix h(comm, nf, ncomp, lo0, keys, leaf, eta);
      h.Compress(rows.data(), eps, mem_node, true);
      h.WriteFile(path, "test_qd_hmatrix_layout", td.data());
      const std::string sum = h.Summary("compressed");
      if (rank == 0) { std::cout << sum; }
      double err = 0.0, ymax = 0.0;
      std::vector<real_t> y(static_cast<size_t>(nb0) * ncomp);
      for (int s = 0; s < 3; s++)
      {
         h.Mult(slips[s].data(), y.data());
         yref[s] = GatherRows(comm, y, lo0, ncomp);
         for (size_t i = 0; i < y.size(); i++)
         {
            const real_t *row = rows.data() + i * 2 * nf;
            real_t d = 0.0;
            for (long long q = 0; q < 2 * nf; q++)
            {
               d += row[q] * slips[s][q];
            }
            err = std::max(err, std::abs(y[i] - d));
            ymax = std::max(ymax, std::abs(d));
         }
      }
      double loc[2] = {err, ymax}, glob[2];
      MPI_Allreduce(loc, glob, 2, MPI_DOUBLE, MPI_MAX, comm);
      std::ostringstream o;
      o << "compressed product vs dense G: max rel. difference "
        << glob[0] / glob[1] << " <= 10 eps";
      Check(glob[0] <= 10.0 * eps * glob[1], o.str(), rank);
   }
   rows.clear();
   rows.shrink_to_fit();

   // 2. and 3. Read back: fake nodes x layouts
   const char *node_sizes[] = {"", "1", "2", "4"};
   for (const char *ns : node_sizes)
   {
      if (ns[0]) { setenv("SEAS_QD_HMAT_NODE_SIZE", ns, 1); }
      else { unsetenv("SEAS_QD_HMAT_NODE_SIZE"); }
      long long count_pieces = 0;
      for (int div : {0, 4, 1000})
      {
         if (div > 0)
         {
            setenv("SEAS_QD_HMAT_PIECE_DIV", std::to_string(div).c_str(), 1);
         }
         QDHMatrix h(comm, nf, ncomp, lo0, keys, leaf, eta);
         std::string info;
         const bool found = h.ReadHeader(path, eps, info);
         if (!found)
         {
            Check(false, "H file found", rank);
            break;
         }
         if (div > 0) { h.PlanLayout(); }
         const std::vector<long long> lo = h.RowLo();
         const int nb = static_cast<int>(lo[rank + 1] - lo[rank]);
         std::vector<real_t> tdr(static_cast<size_t>(nb) * ncomp);
         h.ReadData(path, mem_node, tdr.data());
         std::ostringstream tag;
         tag << "nodes " << (ns[0] ? ns : "real") << ", layout " << h.Layout();
         if (div > 0) { tag << " (piece divisor " << div << ")"; }
         bool same = SameBytes(GatherRows(comm, tdr, lo, ncomp), tdref);
         std::vector<real_t> y(static_cast<size_t>(nb) * ncomp);
         for (int s = 0; s < 3; s++)
         {
            h.Mult(slips[s].data(), y.data());
            same = SameBytes(GatherRows(comm, y, lo, ncomp), yref[s]) && same;
         }
         int ok = (rank != 0 || same) ? 1 : 0, all_ok = 0;
         MPI_Allreduce(&ok, &all_ok, 1, MPI_INT, MPI_MIN, comm);
         Check(all_ok == 1, tag.str() + ": y and T_D bitwise equal to the "
               "compressed product", rank);
         // Bytes held per rank: max / mean, and planned = held
         const double held = h.HeldBytes();
         double hmax = 0.0, hsum = 0.0;
         MPI_Allreduce(&held, &hmax, 1, MPI_DOUBLE, MPI_MAX, comm);
         MPI_Allreduce(&held, &hsum, 1, MPI_DOUBLE, MPI_SUM, comm);
         if (rank == 0)
         {
            std::printf("        largest rank %.0f B, mean %.0f B (%.3fx), "
                        "%lld V^T pieces\n", hmax, hsum / P, hmax * P / hsum,
                        h.NumPieces());
         }
         if (div == 0) { count_pieces = h.NumPieces(); }
         else
         {
            int eq = h.PlannedBytes() == held ? 1 : 0, all_eq = 0;
            MPI_Allreduce(&eq, &all_eq, 1, MPI_INT, MPI_MIN, comm);
            Check(all_eq == 1, tag.str() + ": bytes held = bytes planned "
                  "on every rank", rank);
         }
         if (div == 1000)
         {
            Check(h.NumPieces() > count_pieces, tag.str() + ": V^T split "
                  "into more pieces than the count layout", rank);
         }
         const std::string sum = h.Summary("read");
         if (rank == 0 && div == 4 && !ns[0]) { std::cout << sum; }
      }
      unsetenv("SEAS_QD_HMAT_PIECE_DIV");
   }
   MPI_Barrier(comm);
   if (rank == 0) { std::remove(path.c_str()); }
   int fail_all = 0;
   MPI_Allreduce(&g_fail, &fail_all, 1, MPI_INT, MPI_MAX, comm);
   if (rank == 0)
   {
      std::cout << (fail_all ? "FAILED" : "ALL PASSED") << std::endl;
   }
   return fail_all ? 1 : 0;
}

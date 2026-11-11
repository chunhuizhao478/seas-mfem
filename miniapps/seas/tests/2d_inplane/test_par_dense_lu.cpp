// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/2d_inplane/test_par_dense_lu.cpp — cycle/par_dense_lu.hpp
// (ParDenseLU, the dense LU shared by all ranks) against the serial LU the
// implicit QD fault solve uses on every rank (mfem::LUFactors), and the
// time of both.
//
//   L1  sizes 1, 7, 64, 65, 200, 692 (the sjf1 fault): a dense matrix
//       shaped like the Newton matrix (a smooth kernel that decays with
//       |i - j| plus a diagonal), a row that needs pivoting, 3 right-hand
//       sides: ParDenseLU's solution equals the serial one to 1e-10
//       (relative, max norm) and the residual |A x - b| / |b| is below
//       1e-8 (the scaled row makes x about 1e6 |b|, so 1e-10 is round-off),
//       on every rank.
//   L2  a singular matrix (two equal rows): Factor returns false on every
//       rank.
//   L3  time (-n N, default 692; -reps R, default 3): the serial factor and
//       solve on every rank at once (as the implicit solve does now) and
//       ParDenseLU's, best of R.  Printed, not checked.  UCERF: -n 3290.
//
// Usage (from miniapps/seas):
//   mpirun -np P ./seas_test_par_dense_lu [-n N] [-reps R] [-only-time]

#include "mfem.hpp"
#include "../../cycle/par_dense_lu.hpp"
#include <cmath>
#include <iostream>
#include <vector>

using namespace mfem;

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

/// The test matrix: a kernel like the fault's (decays with distance, both
/// signs) plus a diagonal that changes along the fault; row 2 is scaled so
/// that partial pivoting has to swap.
struct TestMatrix
{
   int M;
   bool singular;
   double operator()(int i, int j) const
   {
      int ii = i;
      if (singular && i == M - 1) { ii = 0; }   // the last row repeats the first
      const double d = std::abs(double(ii - j));
      double v = 1.0 / (1.0 + d * d) - 0.3 / (1.0 + 0.1 * d) + 0.01 * std::sin(0.37 * ii + 0.11 * j);
      if (ii == j) { v += 0.5 + 2.0 * (1.0 + std::cos(0.05 * ii)); }
      if (ii == 2 && M > 4) { v *= 1e-6; }
      return v;
   }
};

double Rhs(int i, int k) { return std::sin(0.3 * i + 1.7 * k) + 0.1 * k; }

void Compare(int M, int rank)
{
   TestMatrix A{M, false};
   // serial LU (column-major), as QDImplicitFault::Factor_
   std::vector<double> a(size_t(M) * M);
   std::vector<int> piv(M);
   for (int j = 0; j < M; j++) { for (int i = 0; i < M; i++) { a[size_t(j) * M + i] = A(i, j); } }
   LUFactors lu(a.data(), piv.data());
   const bool ok_s = lu.Factor(M);
   seas::ParDenseLU plu(MPI_COMM_WORLD);
   const bool ok_p = plu.Factor(M, A);
   double worst = 0.0, worst_res = 0.0;
   for (int k = 0; k < 3; k++)
   {
      std::vector<double> b(M), xs(M), xp(M);
      for (int i = 0; i < M; i++) { b[i] = Rhs(i, k); xs[i] = b[i]; }
      lu.Solve(M, 1, xs.data());
      plu.Solve(b.data(), xp.data());
      double nx = 0.0, nd = 0.0, nb = 0.0, nr = 0.0;
      for (int i = 0; i < M; i++)
      {
         nx = std::max(nx, std::abs(xs[i]));
         nd = std::max(nd, std::abs(xs[i] - xp[i]));
         double r = -b[i];
         for (int j = 0; j < M; j++) { r += A(i, j) * xp[j]; }
         nr = std::max(nr, std::abs(r));
         nb = std::max(nb, std::abs(b[i]));
      }
      worst = std::max(worst, nd / nx);
      worst_res = std::max(worst_res, nr / nb);
   }
   std::ostringstream s;
   s << "L1 M = " << M << ": factors ok, |x_par - x_serial| / |x| = " << worst
     << " (tolerance 1e-10), |A x - b| / |b| = " << worst_res << " (tolerance 1e-8)";
   Check(ok_s && ok_p && worst < 1e-10 && worst_res < 1e-8, s.str(), rank);
}

void Singular(int rank)
{
   const int M = 100;
   TestMatrix A{M, true};
   seas::ParDenseLU plu(MPI_COMM_WORLD);
   Check(!plu.Factor(M, A), "L2 a singular matrix (M = 100): Factor returns false", rank);
}

void Time(int M, int reps, int rank, int P)
{
   TestMatrix A{M, false};
   std::vector<double> b(M), x(M);
   for (int i = 0; i < M; i++) { b[i] = Rhs(i, 0); }
   double tf_s = 1e300, ts_s = 1e300, tf_p = 1e300, ts_p = 1e300, ta_s = 1e300;
   for (int r = 0; r < reps; r++)
   {
      // serial on every rank at once (the implicit solve today)
      std::vector<double> a(size_t(M) * M);
      std::vector<int> piv(M);
      MPI_Barrier(MPI_COMM_WORLD);
      double t0 = MPI_Wtime();
      for (int j = 0; j < M; j++) { for (int i = 0; i < M; i++) { a[size_t(j) * M + i] = A(i, j); } }
      MPI_Barrier(MPI_COMM_WORLD);
      double t1 = MPI_Wtime();
      LUFactors lu(a.data(), piv.data());
      lu.Factor(M);
      MPI_Barrier(MPI_COMM_WORLD);
      double t2 = MPI_Wtime();
      x = b;
      lu.Solve(M, 1, x.data());
      MPI_Barrier(MPI_COMM_WORLD);
      double t3 = MPI_Wtime();
      ta_s = std::min(ta_s, t1 - t0); tf_s = std::min(tf_s, t2 - t1); ts_s = std::min(ts_s, t3 - t2);
      seas::ParDenseLU plu(MPI_COMM_WORLD);
      MPI_Barrier(MPI_COMM_WORLD);
      t0 = MPI_Wtime();
      plu.Factor(M, A);
      MPI_Barrier(MPI_COMM_WORLD);
      t1 = MPI_Wtime();
      plu.Solve(b.data(), x.data());
      MPI_Barrier(MPI_COMM_WORLD);
      t2 = MPI_Wtime();
      tf_p = std::min(tf_p, t1 - t0); ts_p = std::min(ts_p, t2 - t1);
      if (rank == 0 && r == 0)
      {
         std::cout << "  L3 process grid " << plu.GridRows() << " x " << plu.GridCols()
                   << ", block " << seas::ParDenseLU::kBlock << std::endl;
      }
   }
   if (rank == 0)
   {
      std::cout << "  L3 M = " << M << ", " << P << " ranks, best of " << reps << ":\n"
                << "     serial on every rank: assembly " << ta_s << " s, factor " << tf_s
                << " s, solve " << ts_s << " s\n"
                << "     shared (ScaLAPACK):   assembly + factor " << tf_p << " s, solve "
                << ts_p << " s\n"
                << "     factor ratio " << (ta_s + tf_s) / tf_p << ", solve ratio "
                << ts_s / ts_p << std::endl;
   }
}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   const int rank = Mpi::WorldRank(), P = Mpi::WorldSize();
   int n = 692, reps = 3;
   bool only_time = false;
   OptionsParser args(argc, argv);
   args.AddOption(&n, "-n", "--size", "Matrix size of the timing L3.");
   args.AddOption(&reps, "-reps", "--reps", "Repetitions of the timing L3.");
   args.AddOption(&only_time, "-only-time", "--only-time", "-all", "--all", "Run L3 only.");
   args.Parse();
   if (!args.Good()) { if (rank == 0) { args.PrintUsage(std::cout); } return 1; }
   if (rank == 0)
   {
      std::cout << "test_par_dense_lu: " << P << " rank(s), ScaLAPACK "
                << (seas::ParDenseLU::Available() ? "available" : "NOT available") << std::endl;
   }
   if (!seas::ParDenseLU::Available())
   {
      if (rank == 0) { std::cout << "  SKIPPED: this build has no ScaLAPACK" << std::endl; }
      return 0;
   }
   if (!only_time)
   {
      for (int M : {1, 7, 64, 65, 200, 692}) { Compare(M, rank); }
      Singular(rank);
   }
   Time(n, reps, rank, P);
   if (rank == 0)
   {
      std::cout << "test_par_dense_lu: " << g_pass << " PASSED, " << g_fail << " FAILED"
                << std::endl;
   }
   return g_fail ? 1 : 0;
}

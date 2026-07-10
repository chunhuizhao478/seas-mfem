// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_station_tiebreak_np2.cpp — regression guard for the deterministic
// station tie-break (np4_attractor_root_cause_2026-07-10.md).
//
// The historical station->QP assignment broke exact-distance ties by
// rank-LOCAL scan order + lowest-candidate-RANK in the MPI owner
// reduction; different partitions therefore reported different
// (individually correct) physical points — the "np>=4 attractor".
//
// This test constructs a synthetic EXACT two-QP tie with distinct
// coordinates, distributed across np=2 so that the LOWEST rank holds the
// lexicographically LARGER point:
//
//   QP A (lex smaller): (x, y, z) = (-128, 0, -7628),  slip2 = 111
//   QP B             : (x, y, z) = (+128, 0, -7372),  slip2 = 222
//   station "tie" at (along_strike, down_dip) = (0, 7500):
//     A: dx = -128, dz = |z|-7500 = +128   -> dist2 = 2*128^2
//     B: dx = +128, dz = -128              -> dist2 = 2*128^2  (BIT-equal:
//     128 and the z values are exactly representable)
//
// Phase 1 (rank 0 holds B, rank 1 holds A): pre-fix, rank 0 wins the
//   MPI_MIN rank tie-break and the trace reports B (222).  Post-fix the
//   owner is the rank holding the lexicographic winner A -> 111.
// Phase 2 (swapped layout): the reported value must be IDENTICAL (111) —
//   partition-invariance of the selection.
// Phase 3 (no-tie control, station "notie" near B): the unique nearest QP
//   B (222) must be reported in both layouts — the historical pick is
//   preserved for non-tied stations.
// Phase 4 (serial-scan contract): FindNearestDOF_TPV104 with BOTH tied QPs
//   local, B enumerated FIRST, must return A — documents the intentional
//   one-time change of tied picks at np=1 (the pre-fix scan returned the
//   first-enumerated member of the tie).
//
// Runs meaningfully only at np=2 (aborts otherwise).

#include "mfem.hpp"

#include "../../dynamic/tpv104_stations.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace
{

int g_rank = 0;
int g_failed = 0;

void Check(bool ok, const std::string &what)
{
   if (g_rank == 0)
   {
      std::cout << (ok ? "  [PASS] " : "  [FAIL] ") << what << std::endl;
   }
   if (!ok) { g_failed++; }
}

// Coordinates chosen exactly representable so the tie is BIT-exact.
Vector MakeCoord(real_t x, real_t y, real_t z)
{
   Vector v(3);
   v(0) = x; v(1) = y; v(2) = z;
   return v;
}

const real_t kAX = -128.0, kAZ = -7628.0, kASlip = 111.0;   // lex smaller
const real_t kBX = 128.0,  kBZ = -7372.0, kBSlip = 222.0;
const real_t kFarSlip = 999.0;

/// Build this rank's QP set.  `rank_has_A`: this rank holds tie-member A
/// (plus one far QP); the other rank holds B (plus one far QP).
void BuildRankQPs(bool rank_has_A,
                  std::vector<Vector> &coords, std::vector<DOFData> &dofs)
{
   coords.clear(); dofs.clear();
   DOFData far_d; far_d.slip2 = kFarSlip;
   // A far QP FIRST in local order on every rank, so the tie member is
   // never trivially at local index 0.
   coords.push_back(MakeCoord(1500.0, 0.0, -500.0));
   dofs.push_back(far_d);
   DOFData d;
   if (rank_has_A)
   {
      coords.push_back(MakeCoord(kAX, 0.0, kAZ));
      d.slip2 = kASlip;
   }
   else
   {
      coords.push_back(MakeCoord(kBX, 0.0, kBZ));
      d.slip2 = kBSlip;
   }
   dofs.push_back(d);
}

/// Open + write one step through the production writer with EXPLICIT
/// per-rank coords/dofs, then read back the h-slip (column 2) of the
/// requested station's first data row on rank 0.  Returns NaN on any read
/// failure.
real_t RunWriterScenarioCoords(const std::string &out_dir,
                               const std::string &prefix,
                               const std::vector<TPV104Station> &stations,
                               const std::string &read_station_name,
                               const std::vector<Vector> &coords,
                               const std::vector<DOFData> &dofs)
{
   TPV104StationWriter writer;
   writer.Open(out_dir, prefix, stations, coords,
               static_cast<int>(coords.size()), MPI_COMM_WORLD);
   writer.WriteStep(0.0, dofs);
   writer.Flush();
   writer.Close();
   MPI_Barrier(MPI_COMM_WORLD);

   real_t val = std::numeric_limits<real_t>::quiet_NaN();
   if (g_rank == 0)
   {
      const std::string fname = out_dir + "/" + prefix + "_station_"
                                + read_station_name + ".dat";
      std::ifstream in(fname);
      std::string line;
      while (std::getline(in, line))
      {
         if (line.empty() || line[0] == '#') { continue; }
         std::istringstream iss(line);
         real_t t, hslip;
         if (iss >> t >> hslip) { val = hslip; }
         break;
      }
   }
   MPI_Bcast(&val, 1, MPITypeMap<real_t>::mpi_type, 0, MPI_COMM_WORLD);
   return val;
}

/// Open + write one step through the production writer, then read back the
/// h-slip (column 2) of the requested station's first data row on rank 0.
/// Returns NaN on any read failure.
real_t RunWriterScenario(const std::string &out_dir,
                         const std::string &prefix,
                         const std::vector<TPV104Station> &stations,
                         const std::string &read_station_name,
                         bool this_rank_has_A)
{
   std::vector<Vector> coords;
   std::vector<DOFData> dofs;
   BuildRankQPs(this_rank_has_A, coords, dofs);

   TPV104StationWriter writer;
   writer.Open(out_dir, prefix, stations, coords,
               static_cast<int>(coords.size()), MPI_COMM_WORLD);
   writer.WriteStep(0.0, dofs);
   writer.Flush();
   writer.Close();
   MPI_Barrier(MPI_COMM_WORLD);

   real_t val = std::numeric_limits<real_t>::quiet_NaN();
   if (g_rank == 0)
   {
      const std::string fname = out_dir + "/" + prefix + "_station_"
                                + read_station_name + ".dat";
      std::ifstream in(fname);
      std::string line;
      while (std::getline(in, line))
      {
         if (line.empty() || line[0] == '#') { continue; }
         std::istringstream iss(line);
         real_t t, hslip;
         if (iss >> t >> hslip) { val = hslip; }
         break;
      }
   }
   // Everyone learns the value so Check() counts consistently.
   MPI_Bcast(&val, 1, MPITypeMap<real_t>::mpi_type, 0, MPI_COMM_WORLD);
   return val;
}

} // namespace

int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);
   int nprocs = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
   if (nprocs != 2)
   {
      if (g_rank == 0)
      {
         std::cerr << "test_station_tiebreak_np2 requires exactly np=2 "
                      "(got " << nprocs << ")" << std::endl;
      }
      MPI_Finalize();
      return 1;
   }
   if (g_rank == 0)
   {
      std::cout << "=== test_station_tiebreak_np2 ===" << std::endl;
   }

   const std::string out_dir = "output_test_station_tiebreak";
   if (g_rank == 0) { ::mkdir(out_dir.c_str(), 0755); }
   MPI_Barrier(MPI_COMM_WORLD);

   const std::vector<TPV104Station> stations = {
      { 0.0,   7500.0, "tie"   },
      { 100.0, 7400.0, "notie" },
   };

   // Phase 1: rank 0 holds B (lex LARGER), rank 1 holds A (lex smaller).
   // The lexicographic winner A lives on the HIGHER rank — the pre-fix
   // lowest-rank tie-break would report B (222).
   const real_t v1 = RunWriterScenario(out_dir, "p1", stations, "tie",
                                       /*this_rank_has_A=*/g_rank == 1);
   Check(v1 == kASlip,
         "phase 1: tied station reports the lexicographic winner A "
         "(got " + std::to_string(v1) + ", want 111 — 222 means the "
         "pre-fix lowest-rank pick regressed)");

   // Phase 2: swapped layout — selection must be partition-invariant.
   const real_t v2 = RunWriterScenario(out_dir, "p2", stations, "tie",
                                       /*this_rank_has_A=*/g_rank == 0);
   Check(v2 == kASlip,
         "phase 2: swapped layout reports the SAME point A (got "
         + std::to_string(v2) + ")");
   Check(v1 == v2, "phase 1 == phase 2 (partition invariance)");

   // Phase 3: no-tie control — unique nearest QP is B in both layouts.
   const real_t v3a = RunWriterScenario(out_dir, "p3a", stations, "notie",
                                        /*this_rank_has_A=*/g_rank == 1);
   const real_t v3b = RunWriterScenario(out_dir, "p3b", stations, "notie",
                                        /*this_rank_has_A=*/g_rank == 0);
   Check(v3a == kBSlip && v3b == kBSlip,
         "phase 3: non-tied station keeps the historical unique pick B "
         "in both layouts (got " + std::to_string(v3a) + " / "
         + std::to_string(v3b) + ")");

   // Phase 4: serial-scan contract.  Both tie members local, B FIRST in
   // enumeration order: the pre-fix scan returned B (first minimal); the
   // lexicographic rule must return A.
   {
      std::vector<Vector> coords;
      coords.push_back(MakeCoord(kBX, 0.0, kBZ));       // index 0 = B
      coords.push_back(MakeCoord(kAX, 0.0, kAZ));       // index 1 = A
      coords.push_back(MakeCoord(1500.0, 0.0, -500.0)); // far
      const TPV104Station tie_st = { 0.0, 7500.0, "tie" };
      const int pick = FindNearestDOF_TPV104(tie_st, coords, 3);
      Check(pick == 1,
            "phase 4: serial scan picks lexicographic winner A at a tie "
            "even when B is enumerated first (got index "
            + std::to_string(pick) + ")");
      const TPV104Station notie_st = { 100.0, 7400.0, "notie" };
      const int pick2 = FindNearestDOF_TPV104(notie_st, coords, 3);
      Check(pick2 == 0,
            "phase 4: serial scan keeps the unique nearest pick B for a "
            "non-tied station (got index " + std::to_string(pick2) + ")");
   }

   // Phase 5 (R-201, REVIEW_station_tiebreak_2026-07-10.md): the osculating
   // band.  Three QPs by distance from station "tie" (G = |(-128, +128)|):
   //   A  at d = G            (x = -128,  marker 111)
   //   P  at d = G + 0.8*tol  (x = -127,  marker 333)   [inside the window]
   //   q* at d = G + 1.5*tol  (x = -130,  marker 444, lex-SMALLEST)
   //                                                    [outside the window]
   // Pre-R-201 (local-best anchoring), a rank holding {P, q*} would admit
   // q* (0.7*tol above ITS best P) and q* would win the reduction — but
   // only when co-resident with P: the pick depended on the layout.  With
   // the global-min-anchored re-scan, q* is outside the window everywhere
   // and A must win in BOTH layouts.
   {
      const real_t G   = std::sqrt(static_cast<real_t>(2.0) * 128.0 * 128.0);
      const real_t tol = StationTieTol(G);
      auto qp_at = [&](real_t x, real_t d_target, real_t marker,
                       std::vector<Vector> &coords, std::vector<DOFData> &dofs)
      {
         const real_t dz = std::sqrt(d_target * d_target - x * x);
         coords.push_back(MakeCoord(x, 0.0, -(7500.0 + dz)));
         DOFData d; d.slip2 = marker; dofs.push_back(d);
      };
      const std::vector<TPV104Station> st5 = { { 0.0, 7500.0, "tie" } };

      // Layout A: rank 0 = {P, q*, far}; rank 1 = {A, far}.
      std::vector<Vector> cA; std::vector<DOFData> dA;
      if (g_rank == 0)
      {
         qp_at(-127.0, G + 0.8 * tol, 333.0, cA, dA);
         qp_at(-130.0, G + 1.5 * tol, 444.0, cA, dA);
      }
      else
      {
         qp_at(-128.0, G, kASlip, cA, dA);
      }
      cA.push_back(MakeCoord(1500.0, 0.0, -500.0));
      { DOFData fd; fd.slip2 = kFarSlip; dA.push_back(fd); }
      const real_t v5a = RunWriterScenarioCoords(out_dir, "p5a", st5, "tie",
                                                 cA, dA);

      // Layout B: rank 0 = {A, P, far}; rank 1 = {q*, far}.
      std::vector<Vector> cB; std::vector<DOFData> dB;
      if (g_rank == 0)
      {
         qp_at(-128.0, G, kASlip, cB, dB);
         qp_at(-127.0, G + 0.8 * tol, 333.0, cB, dB);
      }
      else
      {
         qp_at(-130.0, G + 1.5 * tol, 444.0, cB, dB);
      }
      cB.push_back(MakeCoord(1500.0, 0.0, -500.0));
      { DOFData fd; fd.slip2 = kFarSlip; dB.push_back(fd); }
      const real_t v5b = RunWriterScenarioCoords(out_dir, "p5b", st5, "tie",
                                                 cB, dB);

      Check(v5a == kASlip && v5b == kASlip,
            "phase 5: osculating-band q* (1.5*tol beyond the global min) is "
            "excluded in BOTH layouts; the true tie winner A is reported "
            "(got " + std::to_string(v5a) + " / " + std::to_string(v5b)
            + "; 444 means local-best anchoring regressed)");
      Check(v5a == v5b,
            "phase 5: layout A == layout B (partition invariance in the "
            "osculating band)");
   }

   int failed_global = 0;
   MPI_Allreduce(&g_failed, &failed_global, 1, MPI_INT, MPI_SUM,
                 MPI_COMM_WORLD);
   if (g_rank == 0)
   {
      std::cout << (failed_global == 0
                       ? "ALL station tie-break checks PASSED"
                       : "station tie-break checks FAILED")
                << std::endl;
   }
   MPI_Finalize();
   return failed_global == 0 ? 0 : 1;
}

// R-104 (REVIEW.md 2026-05-16 round 2) regression test — MPI-only.
//
// This binary exercises the rank-0 gate around the cap-exhausted
// warning in `AdaptiveSchedule::RecomputeIntervalForCap`.  The gate
// is the R-001 morning-round fix that wraps `mfem::out << ...` in an
// `MPI_Comm_rank == 0` check so the warning fires once per ParaView
// instance instead of N times (where N = MPI rank count).
//
// We can't fold this into `test_paraview_schedule_cap` because that
// test's existing R-106 sub-test calls `pv.Save()` on a serial
// `mfem::Mesh`, which on an MPI-initialised build tries to use
// `MFEM_COMM_WORLD` for rank-suffixed file naming and deadlocks two
// ranks on the implicit `MPI_Finalize` at program exit.  This binary
// stays serial-Save-free so it can safely initialise MPI.
//
// Run with:
//   mpirun --oversubscribe -np 2 ./seas_test_paraview_rank0_warning_gate
//
// Exit code 0 = PASS (rank 0 printed exactly one warning, rank 1
// stayed silent).  Exit code 1 = FAIL.
// Skipped (exit 0) with a message on np=1 — the gate has nothing to
// verify on a single rank.

#include "mfem.hpp"
#include "../../io/paraview_output.hpp"

#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>

int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);

   int rank = 0, nranks = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nranks);

   if (nranks < 2)
   {
      if (rank == 0)
      {
         std::cout << "SKIP: needs -np >= 2 (got np=" << nranks
                   << "); the rank-0 gate has nothing to verify on a "
                      "single rank.\n";
      }
      MPI_Finalize();
      return 0;
   }

   // Per-rank capture of mfem::out so we can count which ranks
   // produced the warning.
   std::ostringstream buf;
   std::streambuf *orig = mfem::out.rdbuf(buf.rdbuf());

   mfem::Mesh smesh = mfem::Mesh::MakeCartesian2D(1, 1,
                                                  mfem::Element::TRIANGLE);
   mfem::seas::ParaViewOutput<mfem::Mesh> pv(
      "/tmp/test_paraview_rank0_warning_gate", smesh, /*order=*/1);
   ::mkdir("/tmp/test_paraview_rank0_warning_gate", 0755);

   pv.SetTotalRunTime(100.0);
   pv.GetSchedule().max_total_snapshots = 2;
   // Drop dt_interseismic to a small value so the cap can actually
   // exhaust within the 100 s simulation window.  At the default 1 yr
   // (~3.15e7 s) cadence, 100 s would only fit a single interseismic
   // write and remaining_budget would never reach 0.
   pv.GetSchedule().dt_interseismic = 10.0;
   pv.GetSchedule().Validate();

   // Drive enough interseismic ticks to exhaust the cap (K=2) and
   // then re-enter the cap-exhausted branch (which fires the warning).
   // V_max stays interseismic so we only exercise the cap-aware path
   // — coseismic / nucleation early-return `base` per the soft-cap
   // semantics (REVIEW.md 2026-05-16 round-2 revert).
   for (int i = 0; i < 200; ++i)
   {
      pv.ShouldWrite(i, mfem::real_t(i) * 0.5, 1e-12);
   }

   mfem::out.rdbuf(orig);

   const bool found = (buf.str().find("exhausted") != std::string::npos);
   int local = found ? 1 : 0;
   int global = 0;
   MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   bool pass = true;
   if (rank == 0)
   {
      std::cout << "rank=0  printed='" << (found ? "yes" : "no")
                << "'  global_printers=" << global << "\n";
   }

   if (global != 1)
   {
      if (rank == 0)
      {
         std::cout << "FAIL [R-104]: expected exactly one rank to print "
                      "the cap-exhausted warning, got " << global
                   << ".  The rank-0 gate in "
                      "`RecomputeIntervalForCap` is either missing or "
                      "ineffective.\n";
      }
      pass = false;
   }
   if (rank == 0 && !found)
   {
      std::cout << "FAIL [R-104]: rank 0 did not print.  The gate may "
                   "be wrongly routing the print to a non-zero rank.\n";
      pass = false;
   }

   if (pass && rank == 0)
   {
      std::cout << "PASS [R-104]: cap-exhausted warning fires exactly "
                   "once across " << nranks << " ranks, on rank 0.\n";
   }

   MPI_Finalize();
   return pass ? 0 : 1;
}

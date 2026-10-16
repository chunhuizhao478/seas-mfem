// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// checkpoint.hpp — fully dynamic BP5 cycle (plan rev 3, Section 3.3 i):
// the QD-regime checkpoint.
//
// One file per rank, <prefix>_checkpoint_r<rank>.txt, made of three blocks:
//   1. io/checkpoint.hpp WriteCheckpoint (SEAS_CHECKPOINT_V1): time, dt,
//      step, event counters, state, displacement, traction, slip rate;
//   2. io/petsc_ts_checkpoint.hpp WritePetscTSCheckpoint (PETSC_TS_V2):
//      time, next dt, TS step number, rejections;
//   3. the cycle block (FD_CYCLE_V1), appended here: the commit of the binary
//      that wrote it (rule M10), the IC set, the regime, the step, and
//      FdCycleOutput's state (cadence clock, event state, file sizes).
// Blocks 1 and 2 are written by the existing functions unchanged, so the
// readers of seas_bp5_full can read the first two blocks.
//
// The files are written under <prefix>.partial and renamed once every rank
// has written, so a job killed while writing keeps the previous checkpoint.
// A restart needs the same rank count (block 1 checks it) and, unless
// allowed explicitly, the same commit.
//
// No solver header (rule M4).

#ifndef MFEM_SEAS_CYCLE_CHECKPOINT_HPP
#define MFEM_SEAS_CYCLE_CHECKPOINT_HPP

#include "mfem.hpp"
#include "fd_cycle_output.hpp"
#include "../io/checkpoint.hpp"
#include "../io/petsc_ts_checkpoint.hpp"
#include "../common/mpi_context.hpp"

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>

namespace mfem
{
namespace seas
{
namespace cycle
{

/// Everything a QD checkpoint holds besides the vectors.
struct CheckpointInfo
{
   real_t t = 0.0;
   real_t dt_next = 0.0;
   int step = 0;
   int rejections = 0;
   std::string commit;
   std::string ic;
   std::string regime = "QD";
   FdCycleOutput::State out;
};

/// The vectors of a QD checkpoint (this rank's owned data).
struct CheckpointData
{
   Vector state, displacement, traction, slip_rate;
};

#ifdef MFEM_USE_PETSC

inline void WriteCycleCheckpoint(const std::string &prefix, const CheckpointInfo &info,
                                 const CheckpointData &d, const MPIContext &mpi)
{
   const std::string tmp = prefix + ".partial";
   Vector no_k0;
   WriteCheckpoint(tmp, info.t, info.dt_next, info.step, info.out.num_events,
                   info.out.in_event, d.state, d.displacement, d.traction,
                   d.slip_rate, false, no_k0, &mpi);
   WritePetscTSCheckpoint(tmp, info.t, info.dt_next, info.step, info.rejections,
                          0, -1e30, 0.0, 0, std::numeric_limits<int>::min(), -1e30, &mpi);
   {
      const std::string fn = CheckpointFilename(tmp, mpi.Rank());
      std::ofstream out(fn, std::ios::app);
      MFEM_VERIFY(out.good(), "WriteCycleCheckpoint: cannot append to " << fn);
      const FdCycleOutput::State &o = info.out;
      out << std::setprecision(17) << std::scientific;
      out << "FD_CYCLE_V1\n";
      out << "commit " << info.commit << "\n";
      out << "ic " << info.ic << "\n";
      out << "regime " << info.regime << "\n";
      out << "step " << info.step << "\n";
      out << "out_last_write_time " << o.last_write_time << "\n";
      out << "out_num_events " << o.num_events << "\n";
      out << "out_in_event " << (o.in_event ? 1 : 0) << "\n";
      out << "out_event_start " << o.event_start << "\n";
      out << "out_event_peak_V " << o.event_peak_V << "\n";
      out << "out_below_since " << o.below_since << "\n";
      out << "out_files " << o.files.size() << "\n";
      for (size_t i = 0; i < o.files.size(); i++)
      {
         out << o.files[i] << " " << o.sizes[i] << "\n";
      }
      out << "FD_CYCLE_END\n";
      out.close();
      MFEM_VERIFY(!out.fail(), "WriteCycleCheckpoint: write failed for " << fn);
   }
   mpi.Barrier();
   const std::string from = CheckpointFilename(tmp, mpi.Rank());
   const std::string to = CheckpointFilename(prefix, mpi.Rank());
   MFEM_VERIFY(std::rename(from.c_str(), to.c_str()) == 0,
               "WriteCycleCheckpoint: cannot rename " << from << " to " << to);
   mpi.Barrier();
}

/// Reads all three blocks; aborts on any mismatch.  `allow_other_commit`
/// accepts a checkpoint written by another commit.
inline void ReadCycleCheckpoint(const std::string &prefix, const std::string &commit,
                                bool allow_other_commit, CheckpointInfo &info,
                                CheckpointData &d, const MPIContext &mpi)
{
   int num_eq = 0;
   bool in_eq = false;
   bool fsal = false;
   Vector k0;
   real_t dt_v1 = 0.0;
   MFEM_VERIFY(ReadCheckpoint(prefix, info.t, dt_v1, info.step, num_eq, in_eq, d.state,
                              d.displacement, d.traction, d.slip_rate, fsal, k0, &mpi),
               "ReadCycleCheckpoint: cannot read " << prefix);
   real_t ts_t = 0.0, pv_w = 0.0, pv_v = 0.0, pv_vol = 0.0;
   int ts_step = 0, pv_s = 0, pv_r = 0, pv_c = 0;
   MFEM_VERIFY(ReadPetscTSCheckpoint(prefix, ts_t, info.dt_next, ts_step, info.rejections,
                                     pv_s, pv_w, pv_v, pv_r, pv_c, pv_vol, &mpi),
               "ReadCycleCheckpoint: " << prefix << " has no PETSC_TS_V2 block");
   MFEM_VERIFY(ts_t == info.t && ts_step == info.step,
               "ReadCycleCheckpoint: blocks 1 and 2 disagree (t " << info.t << " vs "
               << ts_t << ", step " << info.step << " vs " << ts_step << ")");

   const std::string fn = CheckpointFilename(prefix, mpi.Rank());
   std::ifstream in(fn);
   std::string tok;
   bool found = false;
   while (in >> tok) { if (tok == "FD_CYCLE_V1") { found = true; break; } }
   MFEM_VERIFY(found, "ReadCycleCheckpoint: " << fn << " has no FD_CYCLE_V1 block "
               "(not written by seas_fd_cycle_driver)");
   auto tag = [&](const char *want)
   {
      std::string t;
      in >> t;
      MFEM_VERIFY(t == want, "ReadCycleCheckpoint: " << fn << ": expected '" << want
                  << "', got '" << t << "'");
   };
   int step = -1, in_event = 0;
   size_t nfiles = 0;
   FdCycleOutput::State &o = info.out;
   tag("commit"); in >> info.commit;
   tag("ic"); in >> info.ic;
   tag("regime"); in >> info.regime;
   tag("step"); in >> step;
   tag("out_last_write_time"); in >> o.last_write_time;
   tag("out_num_events"); in >> o.num_events;
   tag("out_in_event"); in >> in_event;
   tag("out_event_start"); in >> o.event_start;
   tag("out_event_peak_V"); in >> o.event_peak_V;
   tag("out_below_since"); in >> o.below_since;
   tag("out_files"); in >> nfiles;
   o.in_event = (in_event != 0);
   o.files.resize(nfiles);
   o.sizes.resize(nfiles);
   for (size_t i = 0; i < nfiles; i++) { in >> o.files[i] >> o.sizes[i]; }
   tag("FD_CYCLE_END");
   MFEM_VERIFY(step == info.step, "ReadCycleCheckpoint: cycle block step " << step
               << " differs from block 1 step " << info.step);
   MFEM_VERIFY(info.regime == "QD", "ReadCycleCheckpoint: regime " << info.regime);
   MFEM_VERIFY(info.commit == commit || allow_other_commit,
               "ReadCycleCheckpoint: the checkpoint was written by commit " << info.commit
               << ", this binary is " << commit << " (rule M10).  Rebuild that commit, "
               "or pass --allow-other-commit once the change is approved.");
}

#endif // MFEM_USE_PETSC

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_CHECKPOINT_HPP

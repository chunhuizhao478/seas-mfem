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
//   3. the cycle block (FD_CYCLE_V4), appended here: the commit of the binary
//      that wrote it (rule M10), the IC set, the regime, the step, the
//      number of dynamic phases entered so far (the regime controller's
//      state in the QD regime), FdCycleOutput's state (cadence clock,
//      event state, file sizes), each snapshot file's count, last QD
//      mark and last time (snapshots.hpp), a digest of this rank's owned
//      fault DOF positions (the partition: an LTS-aware cut must come out
//      the same), and BenchmarkOutput's state (bench_output.hpp: profile
//      clock, file sizes, the open event's per-DOF arrays).  V1 (P1) had no
//      phase count, which reads as 0, and no catalog last-above time (see
//      ReadCycleCheckpoint); V1 and V2 have no snapshot record; V1 to V3
//      have no digest and no benchmark state.
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
#include "bench_output.hpp"
#include "../io/checkpoint.hpp"
#include "../io/petsc_ts_checkpoint.hpp"
#include "../common/mpi_context.hpp"

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

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
   int dyn_phases = 0;   ///< dynamic phases entered before this checkpoint
   FdCycleOutput::State out;
   std::vector<SnapshotFileState> snapshots;   ///< empty: none (or V1, V2)
   /// V4: this rank's owned fault DOF positions, FaultDigest (NaN: V1-V3).
   real_t fault_digest = std::numeric_limits<real_t>::quiet_NaN();
   bool bench = false;                         ///< V4: [cycle.benchmark] was on
   BenchmarkOutput::State bench_state;
};

/// A number that changes when this rank's owned fault DOFs (positions, 3 N)
/// or their order change: the checkpoint's check that a restart has the
/// same partition.
inline real_t FaultDigest(const Vector &X)
{
   real_t d = 0.0;
   for (int i = 0; i + 2 < X.Size(); i += 3)
   {
      d += (i / 3 + 1) * (X(i) + 2.0 * X(i + 1) + 3.0 * X(i + 2));
   }
   return d;
}

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
      out << "FD_CYCLE_V4\n";
      out << "commit " << info.commit << "\n";
      out << "ic " << info.ic << "\n";
      out << "regime " << info.regime << "\n";
      out << "step " << info.step << "\n";
      out << "dyn_phases " << info.dyn_phases << "\n";
      out << "out_last_write_time " << o.last_write_time << "\n";
      out << "out_num_events " << o.num_events << "\n";
      out << "out_in_event " << (o.in_event ? 1 : 0) << "\n";
      out << "out_event_start " << o.event_start << "\n";
      out << "out_event_peak_V " << o.event_peak_V << "\n";
      out << "out_below_since " << o.below_since << "\n";
      out << "out_last_above " << o.last_above << "\n";
      out << "out_files " << o.files.size() << "\n";
      for (size_t i = 0; i < o.files.size(); i++)
      {
         out << o.files[i] << " " << o.sizes[i] << "\n";
      }
      out << "snapshots " << info.snapshots.size() << "\n";
      for (const SnapshotFileState &sf : info.snapshots)
      {
         out << sf.kind << " " << sf.steps << " " << sf.qd_mark << " " << sf.last_t << "\n";
      }
      out << "fault_digest " << info.fault_digest << "\n";
      out << "bench " << (info.bench ? 1 : 0) << "\n";
      if (info.bench)
      {
         const BenchmarkOutput::State &b = info.bench_state;
         out << "bench_profile_last " << b.profile_last << "\n";
         out << "bench_files " << b.files.size() << "\n";
         for (size_t i = 0; i < b.files.size(); i++)
         {
            out << b.files[i] << " " << b.sizes[i] << "\n";
         }
         out << "bench_event " << (b.open ? 1 : 0) << " " << (b.have_end ? 1 : 0) << "\n";
         const char *names[5] = {"trupt", "slip0", "tau0", "slip_end", "tau_end"};
         const Vector *vs[5] = {&b.trupt, &b.slip0, &b.tau0, &b.slip_end, &b.tau_end};
         for (int k = 0; k < 5; k++)
         {
            out << names[k] << " " << vs[k]->Size();
            for (int i = 0; i < vs[k]->Size(); i++) { out << " " << (*vs[k])(i); }
            out << "\n";
         }
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
   int version = 0;
   while (in >> tok)
   {
      if (tok == "FD_CYCLE_V1") { version = 1; break; }
      if (tok == "FD_CYCLE_V2") { version = 2; break; }
      if (tok == "FD_CYCLE_V3") { version = 3; break; }
      if (tok == "FD_CYCLE_V4") { version = 4; break; }
   }
   MFEM_VERIFY(version > 0, "ReadCycleCheckpoint: " << fn << " has no FD_CYCLE_V1-V4 "
               "block (not written by seas_fd_cycle_driver)");
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
   info.dyn_phases = 0;
   if (version >= 2) { tag("dyn_phases"); in >> info.dyn_phases; }
   tag("out_last_write_time"); in >> o.last_write_time;
   tag("out_num_events"); in >> o.num_events;
   tag("out_in_event"); in >> in_event;
   tag("out_event_start"); in >> o.event_start;
   tag("out_event_peak_V"); in >> o.event_peak_V;
   tag("out_below_since"); in >> o.below_since;
   if (version >= 2) { tag("out_last_above"); in >> o.last_above; }
   else
   {
      // V1 has no last time above: inside an event and below, its first time
      // below (the old dwell, one step longer); above, the checkpoint's t.
      o.last_above = (in_event != 0 && o.below_since >= 0.0) ? o.below_since : info.t;
   }
   tag("out_files"); in >> nfiles;
   o.in_event = (in_event != 0);
   o.files.resize(nfiles);
   o.sizes.resize(nfiles);
   for (size_t i = 0; i < nfiles; i++) { in >> o.files[i] >> o.sizes[i]; }
   info.snapshots.clear();
   if (version >= 3)
   {
      size_t ns = 0;
      tag("snapshots"); in >> ns;
      info.snapshots.resize(ns);
      for (SnapshotFileState &sf : info.snapshots)
      {
         in >> sf.kind >> sf.steps >> sf.qd_mark >> sf.last_t;
      }
   }
   info.bench = false;
   info.fault_digest = std::numeric_limits<real_t>::quiet_NaN();
   if (version >= 4)
   {
      int bench = 0;
      // As a token: NaN (no digest: a checkpoint written outside the driver,
      // as the P3 test's) prints as "nan", which >> for a double cannot read.
      std::string dg;
      tag("fault_digest"); in >> dg;
      info.fault_digest = std::strtod(dg.c_str(), nullptr);
      tag("bench"); in >> bench;
      info.bench = (bench != 0);
      if (info.bench)
      {
         BenchmarkOutput::State &b = info.bench_state;
         size_t nb = 0;
         int open = 0, have_end = 0;
         tag("bench_profile_last"); in >> b.profile_last;
         tag("bench_files"); in >> nb;
         b.files.resize(nb);
         b.sizes.resize(nb);
         for (size_t i = 0; i < nb; i++) { in >> b.files[i] >> b.sizes[i]; }
         tag("bench_event"); in >> open >> have_end;
         b.open = (open != 0);
         b.have_end = (have_end != 0);
         const char *names[5] = {"trupt", "slip0", "tau0", "slip_end", "tau_end"};
         Vector *vs[5] = {&b.trupt, &b.slip0, &b.tau0, &b.slip_end, &b.tau_end};
         for (int k = 0; k < 5; k++)
         {
            int n = 0;
            tag(names[k]); in >> n;
            MFEM_VERIFY(n >= 0, "ReadCycleCheckpoint: " << fn << ": " << names[k] << " size " << n);
            vs[k]->SetSize(n);
            for (int i = 0; i < n; i++) { in >> (*vs[k])(i); }
         }
      }
   }
   MFEM_VERIFY(!in.fail(), "ReadCycleCheckpoint: " << fn << ": cannot read the cycle block");
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

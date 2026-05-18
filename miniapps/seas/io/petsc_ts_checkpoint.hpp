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

#ifndef MFEM_SEAS_PETSC_TS_CHECKPOINT_HPP
#define MFEM_SEAS_PETSC_TS_CHECKPOINT_HPP

#include "mfem.hpp"
#include "../common/mpi_context.hpp"
#include "checkpoint.hpp"   // CheckpointFilename — V1 file layout

#include <fstream>
#include <iomanip>
#include <limits>
#include <string>

#ifdef MFEM_USE_PETSC

namespace mfem
{
namespace seas
{

/// @brief Append PETSc TS internal state to an already-written V1
/// checkpoint.
///
/// Must be called AFTER `WriteCheckpoint` on the same prefix.  Opens
/// the per-rank file in append mode and writes the PETSC_TS_V2
/// trailing block (10 fields, plan §1):
///
///   PETSC_TS_V2
///   petsc_ts_time                    <real_t>
///   petsc_ts_dt_next                 <real_t>
///   petsc_ts_step                    <int>
///   petsc_ts_rejections              <int>   (cumulative across restarts; R-005)
///   paraview_snapshots               <int>
///   paraview_last_write_time         <real_t>  (R-304)
///   paraview_last_v_max              <real_t>  (R-304)
///   paraview_current_regime          <int>     (R-304)
///   paraview_last_committed_cycle    <int>     (R-004)
///   paraview_last_volume_write_time  <real_t>  (R-006)
///
/// `petsc_ts_dt_next` is the dt PETSc will use for the NEXT step, as
/// returned by `TSGetTimeStep` AFTER the just-completed step.
///
/// `petsc_ts_rejections` is the CUMULATIVE rejection count across the
/// entire restart chain — the WRITE-side caller is responsible for
/// adding `restart_rejections_carryover + TSGetStepRejections(ts)`
/// before passing the sum here (R-005).
///
/// The three R-304 schedule-state fields capture the ParaViewOutput
/// adaptive schedule at checkpoint time so the first ShouldWrite call
/// after restart does not fire unconditionally (the default
/// `last_write_time_ = -1e30` would otherwise make `time - last` look
/// like +infinity).
///
/// `paraview_last_committed_cycle` (R-004) is the dedup key used by
/// `ParaViewOutput::CommitSchedule`'s same-step guard.
///
/// `paraview_last_volume_write_time` (R-006) preserves the
/// independent volume-PV cadence across restart for callers that use
/// `--volume-pv-dt`.
inline void WritePetscTSCheckpoint(const std::string &prefix,
                                   real_t t, real_t dt_next,
                                   int step, int rejections,
                                   int paraview_snapshots,
                                   real_t paraview_last_write_time,
                                   real_t paraview_last_v_max,
                                   int paraview_current_regime,
                                   int paraview_last_committed_cycle,
                                   real_t paraview_last_volume_write_time,
                                   const MPIContext *mpi)
{
   const int rank = mpi ? mpi->Rank() : 0;
   const std::string filename = CheckpointFilename(prefix, rank);

   // R-006 (REVIEW.md round 4): probe for the pre-existing V1 file
   // BEFORE opening in append mode.  `std::ios::app` silently CREATES
   // the file if missing, which would leave a header-less file
   // containing only the PETSC_TS_V2 trailing block — `ReadCheckpoint`
   // would later fail to parse it with a confusing "expected
   // 'SEAS_CHECKPOINT_V1', got 'PETSC_TS_V2'" message that points at
   // the reader rather than at the misordered writer.  The explicit
   // probe surfaces the precondition violation at the right place.
   {
      std::ifstream probe(filename);
      MFEM_VERIFY(probe.good(),
                  "WritePetscTSCheckpoint: prerequisite V1 checkpoint "
                  "file " << filename << " does not exist; call "
                  "WriteCheckpoint on the same prefix first.");
      // R-003 (REVIEW.md round 5): probe.good() is true for empty
      // files too; without this size check, a truncated V1 file
      // (e.g., from a crashed prior write) would have V2 appended to
      // it and later confuse ReadCheckpoint with "expected
      // 'SEAS_CHECKPOINT_V1', got 'PETSC_TS_V2'".
      probe.seekg(0, std::ios::end);
      const auto file_size = probe.tellg();
      MFEM_VERIFY(file_size > 0,
                  "WritePetscTSCheckpoint: prerequisite V1 checkpoint "
                  "file " << filename << " exists but is empty "
                  "(prior WriteCheckpoint crashed mid-stream?).  "
                  "Cannot safely append the V2 trailing block.");
   }

   std::ofstream out(filename, std::ios::app);
   MFEM_VERIFY(out.good(),
               "WritePetscTSCheckpoint: cannot open " << filename
               << " for appending the PETSC_TS_V2 trailing block.");

   out << std::setprecision(17) << std::scientific;
   out << "PETSC_TS_V2\n";
   out << "petsc_ts_time " << t << "\n";
   out << "petsc_ts_dt_next " << dt_next << "\n";
   out << "petsc_ts_step " << step << "\n";
   out << "petsc_ts_rejections " << rejections << "\n";
   out << "paraview_snapshots " << paraview_snapshots << "\n";
   out << "paraview_last_write_time " << paraview_last_write_time << "\n";
   out << "paraview_last_v_max " << paraview_last_v_max << "\n";
   out << "paraview_current_regime " << paraview_current_regime << "\n";
   out << "paraview_last_committed_cycle "
       << paraview_last_committed_cycle << "\n";        // R-004
   out << "paraview_last_volume_write_time "
       << paraview_last_volume_write_time << "\n";      // R-006

   // R-012: mirror the WriteCheckpoint pattern — explicit close, then
   // barrier — so the "after this returns, all ranks have committed"
   // invariant holds for the V2 trailing block too.
   out.close();
   if (mpi) { mpi->Barrier(); }
}

/// @brief Read the PETSC_TS_V2 trailing block from a checkpoint file.
///
/// Must be called AFTER `ReadCheckpoint` (which reads the V1 prefix
/// and validates rank-count).  Returns `true` if the trailing
/// PETSC_TS_V2 block was found and parsed; returns `false` (out-params
/// untouched) if the file is V1-only (no trailing block).  Aborts
/// (`MFEM_VERIFY`) on a malformed V2 block.
///
/// The scan past the V1 body uses `in >> tag` to read whitespace-
/// separated tokens until either `PETSC_TS_V2` is matched or EOF is
/// reached.  Numeric tokens in the V1 body cannot match the
/// alpha-numeric tag, so the scan is robust.
///
/// Caller MUST treat a `false` return as a hard error when
/// `--petsc-ts` is active: the V1-only fallback ("restart MFEM state
/// but let PETSc TS start fresh from dt_init") produces a misleading
/// trajectory.  The driver enforces this with `MFEM_VERIFY` in the
/// V2 restart block (see plan §4).
///
/// @return true iff the PETSC_TS_V2 trailing block was parsed.
inline bool ReadPetscTSCheckpoint(const std::string &prefix,
                                  real_t &t, real_t &dt_next,
                                  int &step, int &rejections,
                                  int &paraview_snapshots,
                                  real_t &paraview_last_write_time,
                                  real_t &paraview_last_v_max,
                                  int &paraview_current_regime,
                                  int &paraview_last_committed_cycle,
                                  real_t &paraview_last_volume_write_time,
                                  const MPIContext *mpi)
{
   const int rank = mpi ? mpi->Rank() : 0;
   const std::string filename = CheckpointFilename(prefix, rank);

   std::ifstream in(filename);
   if (!in.good()) { return false; }

   // Scan tokens until PETSC_TS_V2 tag is found or EOF.
   {
      std::string tok;
      bool found = false;
      while (in >> tok)
      {
         if (tok == "PETSC_TS_V2") { found = true; break; }
      }
      if (!found) { return false; }
   }

   // Parse the 10 fields in §1 order.  read_tag aborts on mismatch.
   auto read_tag = [&](const std::string &expected)
   {
      std::string tag;
      in >> tag;
      MFEM_VERIFY(tag == expected,
                  "ReadPetscTSCheckpoint: malformed PETSC_TS_V2 block in "
                  << filename << "; expected tag '" << expected
                  << "', got '" << tag << "'");
   };

   read_tag("petsc_ts_time");                  in >> t;
   read_tag("petsc_ts_dt_next");               in >> dt_next;
   read_tag("petsc_ts_step");                  in >> step;
   read_tag("petsc_ts_rejections");            in >> rejections;
   read_tag("paraview_snapshots");             in >> paraview_snapshots;
   read_tag("paraview_last_write_time");       in >> paraview_last_write_time;
   read_tag("paraview_last_v_max");            in >> paraview_last_v_max;
   read_tag("paraview_current_regime");        in >> paraview_current_regime;
   read_tag("paraview_last_committed_cycle");  in >> paraview_last_committed_cycle;  // R-004
   read_tag("paraview_last_volume_write_time");
   in >> paraview_last_volume_write_time;                                            // R-006

   MFEM_VERIFY(!in.fail(),
               "ReadPetscTSCheckpoint: stream error parsing PETSC_TS_V2 "
               "fields in " << filename);

   in.close();
   if (mpi) { mpi->Barrier(); }
   return true;
}

} // namespace seas
} // namespace mfem

#endif // MFEM_USE_PETSC

#endif // MFEM_SEAS_PETSC_TS_CHECKPOINT_HPP

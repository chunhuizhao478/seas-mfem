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

#ifndef MFEM_SEAS_TPV104_CHECKPOINT_HPP
#define MFEM_SEAS_TPV104_CHECKPOINT_HPP

// TPV104 dynamic-rupture checkpoint format (Phase-4 port of the BP5 V1
// pattern).  Distinct file format from BP5's SEAS_CHECKPOINT_V1 — the
// magic tag is TPV104_CHECKPOINT_V1 — so cross-driver restarts fail
// loudly at the header check rather than corrupting state.
//
// MVP scope (this file):
//   - Wave-field state Q (bulk fluctuation, NUM_STATE * ndof_total)
//   - Per-fault-DOF DYNAMIC state from DOFData: psi, slip_rate,
//     V1/V2, slip1/slip2, tau1_nuc/tau2_nuc/sigma_n_nuc (9 fields)
//   - Scalars: time, dt, step
//   - Per-rank files like BP5 (one file per MPI rank)
//
// Out of scope (Phase-4b follow-up / V3):
//   - ParaView schedule state for BOTH pv_out AND pv_bulk_out
//     collections.  V2 carries the PRIMARY collection only; V3 will
//     extend to per-collection state.
//   - DOFData static fields (impedances, a, Dc, background prestress,
//     LSW params) — these are re-initialized by InitializeFaultDOFs_TPV104
//     and don't need to round-trip.
//   - DOFData corrected-traction fields (tau1_corr, tau2_corr,
//     sigma_n_corr) — these are recomputed by FaultFaceFlux::Evaluate
//     on the first post-restart step.

#include "mfem.hpp"
#include "../common/mpi_context.hpp"
#include "../dynamic/fault_face_flux.hpp"  // DOFData
#include "checkpoint.hpp"                  // CheckpointFilename helper

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{

/// @brief Per-rank file name for a TPV104 checkpoint.
///
/// Reuses BP5's CheckpointFilename helper — the on-disk naming scheme
/// is identical (`{prefix}_checkpoint_r{rank}.txt`), only the MAGIC
/// TAG inside the file differs.  This lets the BP5 vs TPV104 driver
/// distinguish "we own this file" from "this is a foreign-format
/// file" at the first read_tag check.
inline std::string Tpv104CheckpointFilename(const std::string &prefix,
                                            int rank)
{
   return CheckpointFilename(prefix, rank);
}

namespace internal
{

/// Single body for WriteTpv104Checkpoint (raw rank/size).  Both public
/// overloads forward here — R-007: deduplicate so R-002's fix lives in
/// exactly one place.
inline void WriteTpv104CheckpointImpl(const std::string &prefix,
                                      real_t t, real_t dt, int step,
                                      const Vector &Q,
                                      const std::vector<DOFData> &dof_data,
                                      int rank, int size,
                                      const std::string &driver_tag = "")
{
   // Phase 5b of spatial_dynamic_rupture_plan.md (rev-3): the optional
   // DRIVER_TAG_V1 trailer is appended at EOF ONLY when `driver_tag`
   // is non-empty.  When empty (the default), the function returns
   // after the existing V1 body — the output file is byte-identical
   // to what the pre-extension writer produced.  Validator rules: tag
   // must contain no whitespace AND be <= 31 chars.
   MFEM_VERIFY(driver_tag.empty()
               || (driver_tag.size() <= 31
                   && driver_tag.find_first_of(" \t\n\r") == std::string::npos),
               "WriteTpv104Checkpoint: driver_tag '" << driver_tag
               << "' must be <= 31 chars with no whitespace.");

   const std::string filename = Tpv104CheckpointFilename(prefix, rank);
   std::ofstream out(filename);
   MFEM_VERIFY(out.good(),
               "WriteTpv104Checkpoint: cannot open " << filename);

   out << std::setprecision(17) << std::scientific;

   out << "TPV104_CHECKPOINT_V1\n";
   out << "num_ranks " << size << "\n";
   out << "rank " << rank << "\n";
   out << "time " << t << "\n";
   out << "dt " << dt << "\n";
   out << "step " << step << "\n";

   out << "Q_size " << Q.Size() << "\n";
   for (int i = 0; i < Q.Size(); ++i)
   {
      out << Q(i) << "\n";
   }

   const int nd = static_cast<int>(dof_data.size());
   out << "dof_data_size " << nd << "\n";
   for (int i = 0; i < nd; ++i)
   {
      const DOFData &d = dof_data[i];
      out << d.psi         << "\n";
      out << d.slip_rate   << "\n";
      out << d.V1          << "\n";
      out << d.V2          << "\n";
      out << d.slip1       << "\n";
      out << d.slip2       << "\n";
      out << d.tau1_nuc    << "\n";
      out << d.tau2_nuc    << "\n";
      out << d.sigma_n_nuc << "\n";
   }

   // Phase 5b: DRIVER_TAG_V1 trailer.  Skipped when driver_tag is
   // empty so the file stays byte-identical to the pre-extension V1
   // body (TPV/BP5 byte-exact contract).
   if (!driver_tag.empty())
   {
      out << "DRIVER_TAG_V1\n" << driver_tag << "\n";
   }

   out.close();
}

/// Single body for ReadTpv104Checkpoint (raw rank/size).  Both public
/// overloads forward here.  `expected_Q_size` enforces R-002:
/// MUST be >= 0; the file's Q_size MUST equal it or the function
/// aborts via MFEM_VERIFY.  All callers pass the live wave-field
/// size (NUM_STATE * ndof_total) so wrong-mesh restart fails loudly.
/// Sentinel-defaulted "skip the check" mode was removed in R-104 to
/// close the silent-bypass hole.
inline bool ReadTpv104CheckpointImpl(const std::string &prefix,
                                     real_t &t, real_t &dt, int &step,
                                     Vector &Q, int expected_Q_size,
                                     std::vector<DOFData> &dof_data,
                                     int rank, int size,
                                     std::string *driver_tag = nullptr)
{
   const std::string filename = Tpv104CheckpointFilename(prefix, rank);
   std::ifstream in(filename);
   if (!in.good()) { return false; }

   auto read_tag = [&](const std::string &expected)
   {
      std::string tag;
      in >> tag;
      MFEM_VERIFY(tag == expected,
                  "TPV104 checkpoint parse error: expected '"
                  << expected << "', got '" << tag << "' in "
                  << filename);
   };

   read_tag("TPV104_CHECKPOINT_V1");

   int file_num_ranks = 0, file_rank = 0;
   read_tag("num_ranks");  in >> file_num_ranks;
   read_tag("rank");       in >> file_rank;
   MFEM_VERIFY(file_num_ranks == size,
               "TPV104 checkpoint num_ranks mismatch: file has "
               << file_num_ranks << " but running with " << size);
   MFEM_VERIFY(file_rank == rank,
               "TPV104 checkpoint rank mismatch: file has " << file_rank
               << " but this is rank " << rank);

   read_tag("time"); in >> t;
   read_tag("dt");   in >> dt;
   read_tag("step"); in >> step;

   int Q_size = 0;
   read_tag("Q_size"); in >> Q_size;
   MFEM_VERIFY(expected_Q_size >= 0,
               "ReadTpv104Checkpoint: caller passed expected_Q_size="
               << expected_Q_size << " (must be >= 0; this parameter "
               "is REQUIRED to gate wrong-mesh restart per R-002 / "
               "R-104; the sentinel-skip mode was removed).");
   MFEM_VERIFY(Q_size == expected_Q_size,
               "ReadTpv104Checkpoint: Q size mismatch: file has "
               << Q_size << " doubles but the current driver expects "
               << expected_Q_size << " (NUM_STATE * ndof_total). "
               "Different mesh, polynomial order, or partition?");
   Q.SetSize(Q_size);
   for (int i = 0; i < Q_size; ++i) { in >> Q(i); }

   int nd = 0;
   read_tag("dof_data_size"); in >> nd;
   // Resize WITHOUT clearing — preserve any caller-set static fields
   // that already live in dof_data[i].  std::vector::resize() with
   // grow uses default-construct, which would zero the static fields;
   // restart relies on InitializeFaultDOFs_TPV104 having been called
   // first so dof_data is the right size and the static fields are
   // valid.  Enforce size invariant here.
   MFEM_VERIFY(static_cast<int>(dof_data.size()) == nd,
               "TPV104 checkpoint dof_data_size mismatch: file has "
               << nd << " but dof_data.size()=" << dof_data.size()
               << " (caller must InitializeFaultDOFs_TPV104 BEFORE "
               "calling ReadTpv104Checkpoint, with the same fault "
               "DOF count the checkpoint was written from)");

   for (int i = 0; i < nd; ++i)
   {
      DOFData &d = dof_data[i];
      in >> d.psi;
      in >> d.slip_rate;
      in >> d.V1;
      in >> d.V2;
      in >> d.slip1;
      in >> d.slip2;
      in >> d.tau1_nuc;
      in >> d.tau2_nuc;
      in >> d.sigma_n_nuc;
   }

   MFEM_VERIFY(!in.fail(),
               "TPV104 checkpoint stream error in " << filename
               << " after reading " << nd << " DOF data blocks");

   // Phase 5b of spatial_dynamic_rupture_plan.md (rev-3): peek for
   // the optional DRIVER_TAG_V1 trailer.  Pre-extension files end
   // here; peek-extract returns "" in that case (BACK-COMPAT).
   if (driver_tag) { driver_tag->clear(); }
   {
      std::string maybe_tag;
      if (in >> maybe_tag)
      {
         MFEM_VERIFY(maybe_tag == "DRIVER_TAG_V1",
                     "TPV104 checkpoint trailer parse error: got '"
                     << maybe_tag << "' (expected EOF or "
                     "'DRIVER_TAG_V1') in " << filename);
         std::string tag;
         in >> tag;
         MFEM_VERIFY(!in.fail() && !tag.empty(),
                     "TPV104 checkpoint DRIVER_TAG_V1 trailer is "
                     "malformed in " << filename
                     << " (expected a non-empty tag on the line "
                     "after DRIVER_TAG_V1).");
         if (driver_tag) { *driver_tag = tag; }
      }
   }

   in.close();
   return true;
}

// ---------------------------------------------------------------------------
// (Checkpoint V2, LTS Phase 2)  Same text body as V1 but magic
// TPV104_CHECKPOINT_V2 and two extra header fields after `step` (the
// sync-interval counter): `lts_mode` and `layout_hash` (the LtsLayoutHash of the
// clustering, so a restart whose recomputed layout differs is REFUSED).
// dof_data is empty on the fault-free bulk path.
// ---------------------------------------------------------------------------
// `dof_canonical_perm` (LTS Phase 3, P-006): when non-empty, the per-fault-QP
// CANONICAL permutation — `perm[mem]` is the on-disk (canonical, pre-reorder)
// position of the in-memory fault QP `mem`.  `DOFData` is then serialized in
// CANONICAL order so the on-disk layout is independent of whether LTS reordered
// the in-memory fault stack (a restart at a different cluster count / lts mode
// reads the same bytes into its own order).  Empty ⇒ raw in-memory order
// (byte-identical to pre-P3), which is correct for the fault-free bulk path and
// for GTS.  Size (when non-empty) MUST equal dof_data.size().
inline void WriteTpv104CheckpointV2Impl(const std::string &prefix,
                                        real_t t, real_t dt, int sync_step,
                                        int lts_mode, std::uint64_t layout_hash,
                                        const Vector &Q,
                                        const std::vector<DOFData> &dof_data,
                                        int rank, int size,
                                        const std::string &driver_tag = "",
                                        const std::vector<int> *dof_canonical_perm
                                           = nullptr)
{
   MFEM_VERIFY(driver_tag.empty()
               || (driver_tag.size() <= 31
                   && driver_tag.find_first_of(" \t\n\r") == std::string::npos),
               "WriteTpv104CheckpointV2: bad driver_tag '" << driver_tag << "'.");
   const std::string filename = Tpv104CheckpointFilename(prefix, rank);
   std::ofstream out(filename);
   MFEM_VERIFY(out.good(), "WriteTpv104CheckpointV2: cannot open " << filename);
   out << std::setprecision(17) << std::scientific;
   out << "TPV104_CHECKPOINT_V2\n";
   out << "num_ranks " << size << "\n";
   out << "rank " << rank << "\n";
   out << "time " << t << "\n";
   out << "dt " << dt << "\n";
   out << "step " << sync_step << "\n";
   out << "lts_mode " << lts_mode << "\n";
   out << "layout_hash " << layout_hash << "\n";
   out << "Q_size " << Q.Size() << "\n";
   for (int i = 0; i < Q.Size(); ++i) { out << Q(i) << "\n"; }
   const int nd = static_cast<int>(dof_data.size());
   out << "dof_data_size " << nd << "\n";
   // LTS P-006: emit in canonical (on-disk) order.  `inv[canonical] = in-memory`
   // is the inverse of `dof_canonical_perm[in-memory] = canonical`.
   const bool reordered =
      dof_canonical_perm != nullptr && !dof_canonical_perm->empty();
   MFEM_VERIFY(!reordered || static_cast<int>(dof_canonical_perm->size()) == nd,
               "WriteTpv104CheckpointV2: dof_canonical_perm size "
               << dof_canonical_perm->size() << " != dof_data_size " << nd);
   std::vector<int> inv;
   if (reordered)
   {
      inv.assign(nd, -1);
      for (int mem = 0; mem < nd; ++mem)
      {
         const int c = (*dof_canonical_perm)[mem];
         MFEM_VERIFY(c >= 0 && c < nd && inv[c] == -1,
                     "WriteTpv104CheckpointV2: dof_canonical_perm is not a "
                     "permutation of [0, " << nd << ")");
         inv[c] = mem;
      }
   }
   for (int c = 0; c < nd; ++c)
   {
      const int mem = reordered ? inv[c] : c;
      const DOFData &d = dof_data[mem];
      out << d.psi << "\n" << d.slip_rate << "\n" << d.V1 << "\n" << d.V2 << "\n"
          << d.slip1 << "\n" << d.slip2 << "\n" << d.tau1_nuc << "\n"
          << d.tau2_nuc << "\n" << d.sigma_n_nuc << "\n";
   }
   if (!driver_tag.empty()) { out << "DRIVER_TAG_V1\n" << driver_tag << "\n"; }
   out.close();
}

/// Peek the checkpoint magic tag: 1 (V1/GTS), 2 (V2/LTS), or 0 (missing/unknown).
inline int PeekTpv104CheckpointVersion(const std::string &prefix, int rank)
{
   std::ifstream in(Tpv104CheckpointFilename(prefix, rank));
   if (!in.good()) { return 0; }
   std::string tag; in >> tag;
   if (tag == "TPV104_CHECKPOINT_V1") { return 1; }
   if (tag == "TPV104_CHECKPOINT_V2") { return 2; }
   return 0;
}

// `dof_canonical_perm` (LTS Phase 3, P-006): mirror of the write path — when
// non-empty, `perm[mem]` is the on-disk (canonical) position of in-memory fault
// QP `mem`; the on-disk records are read in canonical order and scattered back
// to in-memory order.  Empty ⇒ raw in-memory order (byte-identical to pre-P3).
inline bool ReadTpv104CheckpointV2Impl(const std::string &prefix,
                                       real_t &t, real_t &dt, int &sync_step,
                                       int &lts_mode, std::uint64_t &layout_hash,
                                       Vector &Q, int expected_Q_size,
                                       std::vector<DOFData> &dof_data,
                                       int rank, int size,
                                       std::string *driver_tag = nullptr,
                                       const std::vector<int> *dof_canonical_perm
                                          = nullptr)
{
   const std::string filename = Tpv104CheckpointFilename(prefix, rank);
   std::ifstream in(filename);
   if (!in.good()) { return false; }
   auto read_tag = [&](const std::string &e)
   {
      std::string tag; in >> tag;
      MFEM_VERIFY(tag == e, "TPV104 V2 checkpoint parse: expected '" << e
                  << "', got '" << tag << "' in " << filename);
   };
   read_tag("TPV104_CHECKPOINT_V2");
   int fnr = 0, fr = 0;
   read_tag("num_ranks"); in >> fnr;
   read_tag("rank");      in >> fr;
   MFEM_VERIFY(fnr == size, "V2 checkpoint num_ranks mismatch: " << fnr << " != " << size);
   MFEM_VERIFY(fr == rank, "V2 checkpoint rank mismatch: " << fr << " != " << rank);
   read_tag("time"); in >> t;
   read_tag("dt");   in >> dt;
   read_tag("step"); in >> sync_step;
   read_tag("lts_mode");    in >> lts_mode;
   read_tag("layout_hash"); in >> layout_hash;
   int Q_size = 0; read_tag("Q_size"); in >> Q_size;
   MFEM_VERIFY(expected_Q_size >= 0 && Q_size == expected_Q_size,
               "ReadTpv104CheckpointV2: Q size mismatch: file " << Q_size
               << " != expected " << expected_Q_size);
   Q.SetSize(Q_size);
   for (int i = 0; i < Q_size; ++i) { in >> Q(i); }
   int nd = 0; read_tag("dof_data_size"); in >> nd;
   MFEM_VERIFY(static_cast<int>(dof_data.size()) == nd,
               "V2 checkpoint dof_data_size mismatch: " << nd << " != " << dof_data.size());
   const bool reordered =
      dof_canonical_perm != nullptr && !dof_canonical_perm->empty();
   MFEM_VERIFY(!reordered || static_cast<int>(dof_canonical_perm->size()) == nd,
               "ReadTpv104CheckpointV2: dof_canonical_perm size "
               << dof_canonical_perm->size() << " != dof_data_size " << nd);
   // REVIEW R-004: mirror the WRITE guard — validate a genuine bijection
   // (range AND no duplicates), not just per-element range.  A duplicate would
   // otherwise pass the range check below and SILENTLY scatter wrong (drop one
   // canonical record, double-read another) with no abort.
   if (reordered)
   {
      std::vector<char> seen(nd, 0);
      for (int mem = 0; mem < nd; ++mem)
      {
         const int c = (*dof_canonical_perm)[mem];
         MFEM_VERIFY(c >= 0 && c < nd && !seen[c],
                     "ReadTpv104CheckpointV2: dof_canonical_perm is not a "
                     "permutation of [0, " << nd << ")");
         seen[c] = 1;
      }
   }
   // Read the canonical-order records, then scatter to in-memory order:
   // dof_data[mem].dyn = disk[perm[mem]].dyn.  Only the 9 dynamic fields are
   // read; the static fields the caller pre-seeded are left untouched.
   struct DynRec { real_t psi, slip_rate, V1, V2, slip1, slip2,
                          tau1_nuc, tau2_nuc, sigma_n_nuc; };
   std::vector<DynRec> disk(nd);
   for (int c = 0; c < nd; ++c)
   {
      DynRec &r = disk[c];
      in >> r.psi >> r.slip_rate >> r.V1 >> r.V2 >> r.slip1 >> r.slip2
         >> r.tau1_nuc >> r.tau2_nuc >> r.sigma_n_nuc;
   }
   for (int mem = 0; mem < nd; ++mem)
   {
      const int c = reordered ? (*dof_canonical_perm)[mem] : mem;
      MFEM_VERIFY(c >= 0 && c < nd,
                  "ReadTpv104CheckpointV2: canonical index " << c
                  << " out of range [0, " << nd << ")");
      DOFData &d = dof_data[mem];
      const DynRec &r = disk[c];
      d.psi = r.psi; d.slip_rate = r.slip_rate; d.V1 = r.V1; d.V2 = r.V2;
      d.slip1 = r.slip1; d.slip2 = r.slip2; d.tau1_nuc = r.tau1_nuc;
      d.tau2_nuc = r.tau2_nuc; d.sigma_n_nuc = r.sigma_n_nuc;
   }
   MFEM_VERIFY(!in.fail(), "V2 checkpoint stream error in " << filename);
   if (driver_tag) { driver_tag->clear(); }
   {
      std::string mt;
      if (in >> mt)
      {
         MFEM_VERIFY(mt == "DRIVER_TAG_V1",
                     "V2 checkpoint trailer parse: got '" << mt << "' in " << filename);
         std::string tag; in >> tag;
         // REVIEW CK-1: mirror the V1 reader's fail-loud guard so a truncated /
         // empty-tag trailer aborts rather than silently yielding driver_tag=""
         // (which would bypass the driver's provenance refusal).
         MFEM_VERIFY(!in.fail() && !tag.empty(),
                     "V2 checkpoint DRIVER_TAG_V1 trailer malformed in " << filename
                     << " (expected a non-empty tag after DRIVER_TAG_V1).");
         if (driver_tag) { *driver_tag = tag; }
      }
   }
   in.close();
   return true;
}

} // namespace internal

/// @brief Write a TPV104 V1 checkpoint (one file per MPI rank).
///
/// Format (plain text, 17-digit scientific for floats — round-trips
/// IEEE 754 doubles exactly):
///
///   TPV104_CHECKPOINT_V1
///   num_ranks   N
///   rank        R
///   time        T
///   dt          DT
///   step        S
///   Q_size      NQ
///   <NQ real_t values, one per line>
///   dof_data_size  ND
///   <ND blocks, each block is 9 reals on consecutive lines in the
///    canonical order: psi, slip_rate, V1, V2, slip1, slip2,
///                     tau1_nuc, tau2_nuc, sigma_n_nuc>
///
/// All ranks call collectively.  Each rank writes its own file at
/// CheckpointFilename(prefix, rank).  A barrier at the end guarantees
/// every rank has committed before the function returns.
inline void WriteTpv104Checkpoint(const std::string &prefix,
                                  real_t t, real_t dt, int step,
                                  const Vector &Q,
                                  const std::vector<DOFData> &dof_data,
                                  const MPIContext *mpi,
                                  const std::string &driver_tag = "")
{
   const int rank = mpi ? mpi->Rank() : 0;
   const int size = mpi ? mpi->Size() : 1;

   internal::WriteTpv104CheckpointImpl(prefix, t, dt, step, Q, dof_data,
                                       rank, size, driver_tag);

   if (mpi) { mpi->Barrier(); }

   if (!mpi || mpi->IsRoot())
   {
      mfem::out << "TPV104 checkpoint written: " << prefix
                << " (t=" << t << " s, step=" << step << ")\n";
   }
}

/// @brief Raw-MPI overload of WriteTpv104Checkpoint.
///
/// The TPV104 driver uses raw MPI directly (not the seas::MPIContext
/// wrapper) and the MPIContext class is the serial stub when the
/// driver's translation unit doesn't define SEAS_USE_MPI.  This
/// overload accepts (rank, size, comm) directly so the TPV104 driver
/// can call into the checkpoint machinery without pulling in
/// SEAS_USE_MPI globally.  BP5 callers continue to use the
/// MPIContext* signature above.
inline void WriteTpv104Checkpoint(const std::string &prefix,
                                  real_t t, real_t dt, int step,
                                  const Vector &Q,
                                  const std::vector<DOFData> &dof_data,
                                  int rank, int size
#ifdef MFEM_USE_MPI
                                  , MPI_Comm comm = MPI_COMM_NULL
#endif
                                  , const std::string &driver_tag = ""
                                  )
{
   internal::WriteTpv104CheckpointImpl(prefix, t, dt, step, Q, dof_data,
                                       rank, size, driver_tag);

#ifdef MFEM_USE_MPI
   if (comm != MPI_COMM_NULL) { MPI_Barrier(comm); }
#endif

   if (rank == 0)
   {
      mfem::out << "TPV104 checkpoint written: " << prefix
                << " (t=" << t << " s, step=" << step << ")\n";
   }
}

/// @brief Read a TPV104 V1 checkpoint (one file per MPI rank).
///
/// Counterpart to WriteTpv104Checkpoint.  Returns true on success;
/// false if the per-rank file does not exist (so the caller can decide
/// whether to fall back to a fresh run or abort).  ABORTS via
/// MFEM_VERIFY if the file exists but has a wrong magic tag, wrong
/// rank metadata, wrong Q size (when `expected_Q_size >= 0`), or
/// wrong field sequence.
///
/// `expected_Q_size` (R-002 / R-104): callers MUST pass the live
/// mesh's `NUM_STATE * ndof_total` (or, for unit tests, the size of
/// the synthetic Q they wrote).  Must be >= 0; the file's Q_size
/// must equal it or the function aborts via MFEM_VERIFY.  The
/// sentinel-skip mode that allowed -1 was removed in R-104.
///
/// `Q` is resized to the on-disk size (caller's pre-resize value is
/// discarded after the size check).  `dof_data` is resized to the
/// on-disk size; each element's DYNAMIC fields are overwritten while
/// the STATIC fields (impedances, a, Dc, prestress, LSW params) are
/// LEFT ALONE — the caller must have already initialised them via
/// InitializeFaultDOFs_TPV104 before calling this.  This invariant
/// mirrors BP5's pattern where ReadCheckpoint loads slip+psi but
/// leaves seas_op's static mesh/material setup intact.
inline bool ReadTpv104Checkpoint(const std::string &prefix,
                                 real_t &t, real_t &dt, int &step,
                                 Vector &Q, int expected_Q_size,
                                 std::vector<DOFData> &dof_data,
                                 const MPIContext *mpi,
                                 std::string *driver_tag = nullptr)
{
   const int rank = mpi ? mpi->Rank() : 0;
   const int size = mpi ? mpi->Size() : 1;

   const bool ok =
      internal::ReadTpv104CheckpointImpl(prefix, t, dt, step, Q,
                                         expected_Q_size, dof_data,
                                         rank, size, driver_tag);
   if (!ok) { return false; }

   if (mpi) { mpi->Barrier(); }

   if (!mpi || mpi->IsRoot())
   {
      mfem::out << "TPV104 checkpoint loaded: " << prefix
                << " (t=" << t << " s, step=" << step << ")\n";
   }
   return true;
}

/// @brief Raw-MPI overload of ReadTpv104Checkpoint.  See
/// WriteTpv104Checkpoint raw-MPI overload above for the rationale.
inline bool ReadTpv104Checkpoint(const std::string &prefix,
                                 real_t &t, real_t &dt, int &step,
                                 Vector &Q, int expected_Q_size,
                                 std::vector<DOFData> &dof_data,
                                 int rank, int size
#ifdef MFEM_USE_MPI
                                 , MPI_Comm comm = MPI_COMM_NULL
#endif
                                 , std::string *driver_tag = nullptr
                                 )
{
   const bool ok =
      internal::ReadTpv104CheckpointImpl(prefix, t, dt, step, Q,
                                         expected_Q_size, dof_data,
                                         rank, size, driver_tag);
   if (!ok) { return false; }

#ifdef MFEM_USE_MPI
   if (comm != MPI_COMM_NULL) { MPI_Barrier(comm); }
#endif

   if (rank == 0)
   {
      mfem::out << "TPV104 checkpoint loaded: " << prefix
                << " (t=" << t << " s, step=" << step << ")\n";
   }
   return true;
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_TPV104_CHECKPOINT_HPP

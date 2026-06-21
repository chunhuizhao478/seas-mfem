// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// seas_qd_checkpoint.hpp — checkpoint/restart for the quasi-dynamic
// spatial_seas driver (Phase 6 of
// document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md).
//
// QD-NATIVE (R-009): this header deliberately does NOT reuse
// io/tpv104_checkpoint.hpp.  That one #includes dynamic/fault_face_flux.hpp
// (DOFData) and is Q + std::vector<DOFData>-centric — reusing it would violate
// the plan's "No dynamic-code dependency" constraint and does not fit the QD
// state, which is a single flat Vector (no bulk wave field Q).  Only the
// per-rank ASCII framing + DRIVER_TAG_V1 trailer *logic* is mirrored here.
//
// Persists, per MPI rank, to "{prefix}_checkpoint_r{rank}.txt":
//   - a MAGIC line (cross-driver restarts fail on mismatch),
//   - (t, dt, step) and n_owned,
//   - the flat fault state vector ([s_dip, s_strike, psi] * n_owned),
//   - per-DOF static params (a, Dc, eta, V_init, sigma_n; tau_pre is 2*n_owned),
//   - a DRIVER_TAG_V1 trailer (default "spatial_seas").
//
// Full IEEE round-trip precision (max_digits10) so a restart reproduces the
// un-checkpointed trajectory TO ROUND-OFF (Phase-6 acceptance).

#ifndef MFEM_SEAS_QD_CHECKPOINT_HPP
#define MFEM_SEAS_QD_CHECKPOINT_HPP

#include "mfem.hpp"

#include <fstream>
#include <iomanip>
#include <limits>
#include <string>

namespace mfem
{
namespace seas
{

/// Magic line — bump the suffix on any incompatible format change.
inline const char *SpatialSeasCheckpointMagic() { return "SPATIAL_SEAS_CHECKPOINT_V1"; }

/// Per-rank checkpoint filename: "{prefix}_checkpoint_r{rank}.txt".
inline std::string SpatialSeasCheckpointFilename(const std::string &prefix, int rank)
{
   return prefix + "_checkpoint_r" + std::to_string(rank) + ".txt";
}

/// The full QD checkpoint payload for one rank.  `state` is the flat fault
/// state ([s_dip, s_strike, psi] * n_owned, size == fault_op.StateSize()).
/// The static-param vectors are size n_owned each, except `tau_pre` (2*n_owned).
struct SpatialSeasCheckpoint
{
   real_t t    = 0.0;
   real_t dt   = 0.0;
   int    step = 0;
   int    n_owned = 0;
   Vector state;     ///< 3 * n_owned
   Vector a;         ///< n_owned
   Vector Dc;        ///< n_owned
   Vector eta;       ///< n_owned
   Vector V_init;    ///< n_owned
   Vector sigma_n;   ///< n_owned
   Vector tau_pre;   ///< 2 * n_owned (dip, strike per DOF)
};

namespace detail
{
/// Write a labelled vector: "<label> <size>\n" then one value per line.
inline void WriteLabelledVector(std::ostream &out, const char *label,
                                const Vector &v)
{
   out << label << ' ' << v.Size() << '\n';
   for (int i = 0; i < v.Size(); ++i) { out << v(i) << '\n'; }
}

/// Read a labelled vector; verifies the label and (if >=0) the expected size.
inline bool ReadLabelledVector(std::istream &in, const char *label, Vector &v,
                               int expected_size = -1)
{
   std::string got; int n = -1;
   in >> got >> n;
   if (!in || got != label || n < 0) { return false; }
   if (expected_size >= 0 && n != expected_size) { return false; }
   v.SetSize(n);
   for (int i = 0; i < n; ++i) { in >> v(i); }
   return static_cast<bool>(in);
}
}  // namespace detail

/// Write the QD checkpoint for `rank`.  `driver_tag` (default "spatial_seas")
/// is written as a DRIVER_TAG_V1 trailer; ReadSpatialSeasCheckpoint refuses a
/// mismatched tag (prevents cross-driver restarts).
inline void WriteSpatialSeasCheckpoint(const std::string &prefix, int rank,
                                       const SpatialSeasCheckpoint &c,
                                       const std::string &driver_tag = "spatial_seas")
{
   MFEM_VERIFY(c.state.Size() == 3 * c.n_owned,
               "WriteSpatialSeasCheckpoint: state size " << c.state.Size()
               << " != 3*n_owned (" << 3 * c.n_owned << ")");
   MFEM_VERIFY(c.tau_pre.Size() == 2 * c.n_owned,
               "WriteSpatialSeasCheckpoint: tau_pre size " << c.tau_pre.Size()
               << " != 2*n_owned (" << 2 * c.n_owned << ")");
   MFEM_VERIFY(driver_tag.empty()
               || (driver_tag.size() <= 31
                   && driver_tag.find_first_of(" \t\n\r") == std::string::npos),
               "WriteSpatialSeasCheckpoint: driver_tag '" << driver_tag
               << "' must be <=31 chars and whitespace-free.");

   const std::string fname = SpatialSeasCheckpointFilename(prefix, rank);
   std::ofstream out(fname);
   MFEM_VERIFY(out.good(),
               "WriteSpatialSeasCheckpoint: cannot open " << fname << " for write.");
   out << std::setprecision(std::numeric_limits<real_t>::max_digits10);

   out << SpatialSeasCheckpointMagic() << '\n';
   out << "t " << c.t << '\n';
   out << "dt " << c.dt << '\n';
   out << "step " << c.step << '\n';
   out << "n_owned " << c.n_owned << '\n';
   detail::WriteLabelledVector(out, "state",   c.state);
   detail::WriteLabelledVector(out, "a",       c.a);
   detail::WriteLabelledVector(out, "Dc",      c.Dc);
   detail::WriteLabelledVector(out, "eta",     c.eta);
   detail::WriteLabelledVector(out, "V_init",  c.V_init);
   detail::WriteLabelledVector(out, "sigma_n", c.sigma_n);
   detail::WriteLabelledVector(out, "tau_pre", c.tau_pre);
   if (!driver_tag.empty()) { out << "DRIVER_TAG_V1\n" << driver_tag << '\n'; }
   MFEM_VERIFY(out.good(),
               "WriteSpatialSeasCheckpoint: write to " << fname << " failed.");
}

/// Read the QD checkpoint for `rank`.  Returns false (and leaves `c`
/// unspecified) on a missing file, a magic mismatch, or a malformed body.
/// If `expected_driver_tag` is non-empty, a mismatched trailer ⇒ false.
/// If `driver_tag_out` is non-null, the read tag is stored there.
inline bool ReadSpatialSeasCheckpoint(const std::string &prefix, int rank,
                                      SpatialSeasCheckpoint &c,
                                      const std::string &expected_driver_tag = "spatial_seas",
                                      std::string *driver_tag_out = nullptr)
{
   const std::string fname = SpatialSeasCheckpointFilename(prefix, rank);
   std::ifstream in(fname);
   if (!in.good()) { return false; }

   std::string magic;
   std::getline(in, magic);
   if (magic != SpatialSeasCheckpointMagic()) { return false; }

   std::string key;
   auto read_scalar = [&](const char *label, auto &dst) -> bool
   {
      std::string got;
      in >> got >> dst;
      return static_cast<bool>(in) && got == label;
   };
   if (!read_scalar("t", c.t))           { return false; }
   if (!read_scalar("dt", c.dt))         { return false; }
   if (!read_scalar("step", c.step))     { return false; }
   if (!read_scalar("n_owned", c.n_owned)) { return false; }
   if (c.n_owned < 0) { return false; }

   if (!detail::ReadLabelledVector(in, "state",   c.state,   3 * c.n_owned)) { return false; }
   if (!detail::ReadLabelledVector(in, "a",       c.a,       c.n_owned))     { return false; }
   if (!detail::ReadLabelledVector(in, "Dc",      c.Dc,      c.n_owned))     { return false; }
   if (!detail::ReadLabelledVector(in, "eta",     c.eta,     c.n_owned))     { return false; }
   if (!detail::ReadLabelledVector(in, "V_init",  c.V_init,  c.n_owned))     { return false; }
   if (!detail::ReadLabelledVector(in, "sigma_n", c.sigma_n, c.n_owned))     { return false; }
   if (!detail::ReadLabelledVector(in, "tau_pre", c.tau_pre, 2 * c.n_owned)) { return false; }

   // Optional DRIVER_TAG_V1 trailer.
   std::string tag;
   if (in >> key && key == "DRIVER_TAG_V1") { in >> tag; }
   if (driver_tag_out) { *driver_tag_out = tag; }
   if (!expected_driver_tag.empty() && tag != expected_driver_tag) { return false; }
   return true;
}

}  // namespace seas
}  // namespace mfem

#endif  // MFEM_SEAS_QD_CHECKPOINT_HPP

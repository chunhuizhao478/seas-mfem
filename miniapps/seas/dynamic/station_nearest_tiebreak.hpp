// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/station_nearest_tiebreak.hpp — deterministic, partition-invariant
// station -> fault-QP assignment, shared by the TPV102 / TPV104 / TPV205 /
// TPV31 station writers.
//
// WHY (debug_document/tpv104_debug_document/np4_attractor_root_cause_
// 2026-07-10.md): a station that sits exactly equidistant between two or
// more fault QPs (e.g. every x2 = 0 SCEC station on an x-symmetric mesh —
// on tpv104_symmirror_1000m the hypocenter station ties between the QPs at
// x = -166.7 m and x = +166.7 m, 400 m apart) used to be resolved by the
// rank-LOCAL scan order of FindNearestDOF_* plus the lowest-candidate-RANK
// pick in the MPI owner reduction.  Both are partition-dependent, so
// different rank counts reported different — individually correct —
// physical sample points.  That readout artifact masqueraded for months as
// a "np>=4 attractor" (a 4% hypocenter strike-slip "jump"); the simulated
// fields themselves are bit-exactly partition-invariant.
//
// RULE: among QPs whose station distance ties within
//     tie_tol(d) = max(1e-10, 1e-9 * max(d, 1))          [R4-007 formula]
// of the minimum, pick the lexicographically smallest PHYSICAL coordinate
// (x, then z, then y).  Rank id enters only as the very last tie among
// BIT-EQUAL coordinates (the two duplicated copies of one shared fault QP).
// A station with a unique nearest QP resolves to the identical QP as the
// historical strict-`<` scan, so non-tied traces are byte-unchanged.

#ifndef MFEM_SEAS_STATION_NEAREST_TIEBREAK_HPP
#define MFEM_SEAS_STATION_NEAREST_TIEBREAK_HPP

#include "mfem.hpp"

#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{

/// R4-007 tie tolerance on a sqrt-distance (see header comment).
inline real_t StationTieTol(real_t dist)
{
   return std::max<real_t>(static_cast<real_t>(1e-10),
                           static_cast<real_t>(1e-9)
                             * std::max<real_t>(dist,
                                                static_cast<real_t>(1.0)));
}

/// @brief Nearest fault QP to (along_strike, down_dip) with the
/// deterministic lexicographic tie-break.
///
/// Distance metric (identical to every historical FindNearestDOF_*):
/// dx = x - along_strike, dz = |z| - down_dip, dist2 = dx^2 + dz^2.
///
/// Pass 1 is the historical strict-`<` scan; pass 2 replaces the pick by
/// the lexicographic (x, z, y) minimum among QPs inside the tie window.
///
/// Window anchoring (R-201, REVIEW_station_tiebreak_2026-07-10.md):
/// - `anchor_dist < 0` (default; the serial path): window anchored at the
///   rank-LOCAL best distance, `d - d_best < StationTieTol(d_best)`.  For
///   a non-tied station the window is {best} and the historical pick is
///   returned unchanged.  Never returns -1 for `ndof > 0`.
/// - `anchor_dist >= 0` (the MPI re-scan): window `|d - anchor_dist| <
///   window_tol` (with `window_tol <= 0` defaulting to
///   `StationTieTol(anchor_dist)`).  Returns -1 when NO local QP is in
///   the window — anchoring every rank at the GLOBAL minimum makes the
///   candidate set ≡ the true global tie set, closing the osculating-band
///   partition dependence of local-best anchoring.
///
/// @return Index of the selected DOF; -1 if `ndof <= 0`, or (anchored
///         mode only) if the window is empty on this rank.
inline int FindNearestFaultDOFLex(const std::vector<Vector> &fault_coords,
                                  int ndof,
                                  real_t along_strike,
                                  real_t down_dip,
                                  real_t anchor_dist = static_cast<real_t>(-1.0),
                                  real_t window_tol  = static_cast<real_t>(-1.0))
{
   if (ndof <= 0) { return -1; }

   int best = 0;
   real_t best_dist2 = std::numeric_limits<real_t>::max();
   for (int i = 0; i < ndof; ++i)
   {
      const real_t dx = fault_coords[i](0) - along_strike;
      const real_t dz = std::abs(fault_coords[i](2)) - down_dip;
      const real_t dist2 = dx * dx + dz * dz;
      if (dist2 < best_dist2)
      {
         best_dist2 = dist2;
         best = i;
      }
   }

   const real_t d_best = std::sqrt(best_dist2);
   const real_t anchor = (anchor_dist >= static_cast<real_t>(0.0))
                            ? anchor_dist : d_best;
   const real_t tol = (window_tol > static_cast<real_t>(0.0))
                         ? window_tol : StationTieTol(anchor);
   int lex = -1;
   for (int i = 0; i < ndof; ++i)
   {
      const real_t dx = fault_coords[i](0) - along_strike;
      const real_t dz = std::abs(fault_coords[i](2)) - down_dip;
      const real_t d = std::sqrt(dx * dx + dz * dz);
      if (std::abs(d - anchor) >= tol) { continue; }
      if (lex < 0)
      {
         lex = i;
         continue;
      }
      const Vector &ci = fault_coords[i];
      const Vector &cl = fault_coords[lex];
      if (ci(0) < cl(0) ||
          (ci(0) == cl(0) &&
           (ci(2) < cl(2) ||
            (ci(2) == cl(2) && ci(1) < cl(1)))))
      {
         lex = i;
      }
   }
   // Local-best anchoring always contains `best` (|d_best - d_best| = 0 <
   // tol), so lex >= 0 there; -1 escapes only from an empty anchored window.
   return lex;
}

#ifdef MFEM_USE_MPI
/// @brief MPI station-owner resolution by lexicographic COORDINATE.
///
/// Among the ranks whose local best distance is a candidate (caller-side
/// gate, unchanged per writer), the owner is the rank holding the
/// lexicographically smallest (x, z, y) candidate QP; rank id breaks the
/// final tie only among bit-equal coordinates (duplicated shared-QP
/// copies).  COLLECTIVE on `comm`: every rank must call for every station,
/// candidates or not (non-candidates pass any coordinates; they are
/// masked out).
///
/// @return true iff the calling rank owns the station.
inline bool ResolveStationOwnerLex(MPI_Comm comm,
                                   bool is_candidate,
                                   real_t x, real_t z, real_t y,
                                   const char *writer,
                                   const std::string &station_name)
{
   int my_rank, nprocs;
   MPI_Comm_rank(comm, &my_rank);
   MPI_Comm_size(comm, &nprocs);
   const real_t RMAX = std::numeric_limits<real_t>::max();

   real_t cand_x = is_candidate ? x : RMAX;
   real_t win_x;
   MPI_Allreduce(&cand_x, &win_x, 1, MPITypeMap<real_t>::mpi_type,
                 MPI_MIN, comm);

   real_t cand_z = (is_candidate && x == win_x) ? z : RMAX;
   real_t win_z;
   MPI_Allreduce(&cand_z, &win_z, 1, MPITypeMap<real_t>::mpi_type,
                 MPI_MIN, comm);

   real_t cand_y = (is_candidate && x == win_x && z == win_z) ? y : RMAX;
   real_t win_y;
   MPI_Allreduce(&cand_y, &win_y, 1, MPITypeMap<real_t>::mpi_type,
                 MPI_MIN, comm);

   const bool coord_winner =
      is_candidate && x == win_x && z == win_z && y == win_y;
   const int candidate_rank = coord_winner ? my_rank : nprocs;
   int winning_rank;
   MPI_Allreduce(&candidate_rank, &winning_rank, 1, MPI_INT, MPI_MIN, comm);

   // R4-007 fail-loud contract: SOME rank must win, otherwise the station
   // is silently dropped and the trace file never opens.
   MFEM_VERIFY(winning_rank < nprocs,
               writer << "::Open: station " << station_name
               << " has no winning rank (no distance candidate submitted "
               "coordinates).  Most likely: no rank owns any fault QP "
               "(fault-less mesh / wrong fault attribute), or a NaN fault-QP "
               "coordinate (NaN fails the coordinate-equality rounds).");
   return my_rank == winning_rank;
}
#endif // MFEM_USE_MPI

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_STATION_NEAREST_TIEBREAK_HPP

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

#ifndef MFEM_SEAS_CHECKPOINT_HPP
#define MFEM_SEAS_CHECKPOINT_HPP

#include "mfem.hpp"
#include "../common/mpi_context.hpp"

#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

namespace mfem
{
namespace seas
{

/// @brief Generate checkpoint filename for a given rank.
///
/// Format: {prefix}_checkpoint_r{rank}.txt
inline std::string CheckpointFilename(const std::string &prefix, int rank)
{
   std::ostringstream oss;
   oss << prefix << "_checkpoint_r" << rank << ".txt";
   return oss.str();
}

/// @brief Write simulation checkpoint to disk.
///
/// Each MPI rank writes its own file containing local data.
/// All floating-point values use 17-digit precision for exact round-trip.
///
/// @param[in] prefix     Output file prefix (rank suffix appended automatically)
/// @param[in] t          Current simulation time [s]
/// @param[in] dt         Current adaptive time step [s]
/// @param[in] step       Number of accepted steps
/// @param[in] num_eq     Number of seismic events detected
/// @param[in] in_eq      Whether currently in a seismic event
/// @param[in] state      ODE state vector [slip_0, theta_0, ...]
/// @param[in] displacement  Domain displacement solution
/// @param[in] traction   Shear stress at fault DOFs
/// @param[in] slip_rate  Slip rate at fault DOFs
/// @param[in] fsal_init  Whether FSAL stage k_[0] is valid
/// @param[in] k0         FSAL stage vector from time stepper
/// @param[in] mpi        MPI context (may be nullptr for serial)
inline void WriteCheckpoint(const std::string &prefix,
                            real_t t, real_t dt,
                            int step, int num_eq, bool in_eq,
                            const Vector &state,
                            const Vector &displacement,
                            const Vector &traction,
                            const Vector &slip_rate,
                            bool fsal_init, const Vector &k0,
                            const MPIContext *mpi)
{
   int rank = mpi ? mpi->Rank() : 0;
   int size = mpi ? mpi->Size() : 1;

   std::string filename = CheckpointFilename(prefix, rank);
   std::ofstream out(filename);
   MFEM_VERIFY(out.good(), "Cannot open checkpoint file: " << filename);

   out << std::setprecision(17) << std::scientific;

   // Header
   out << "SEAS_CHECKPOINT_V1\n";
   out << "num_ranks " << size << "\n";
   out << "rank " << rank << "\n";
   out << "time " << t << "\n";
   out << "dt " << dt << "\n";
   out << "step " << step << "\n";
   out << "num_seismic_events " << num_eq << "\n";
   out << "in_seismic_event " << (in_eq ? 1 : 0) << "\n";

   // State vector
   out << "state_size " << state.Size() << "\n";
   for (int i = 0; i < state.Size(); i++)
   {
      out << state(i) << "\n";
   }

   // Displacement
   out << "displacement_size " << displacement.Size() << "\n";
   for (int i = 0; i < displacement.Size(); i++)
   {
      out << displacement(i) << "\n";
   }

   // Traction
   out << "traction_size " << traction.Size() << "\n";
   for (int i = 0; i < traction.Size(); i++)
   {
      out << traction(i) << "\n";
   }

   // Slip rate
   out << "slip_rate_size " << slip_rate.Size() << "\n";
   for (int i = 0; i < slip_rate.Size(); i++)
   {
      out << slip_rate(i) << "\n";
   }

   // FSAL state
   out << "fsal_initialized " << (fsal_init ? 1 : 0) << "\n";
   out << "k0_size " << k0.Size() << "\n";
   for (int i = 0; i < k0.Size(); i++)
   {
      out << k0(i) << "\n";
   }

   out.close();

   if (mpi) { mpi->Barrier(); }

   if (!mpi || mpi->IsRoot())
   {
      mfem::out << "Checkpoint written: " << prefix
                << " (t=" << t << " s, step=" << step << ")\n";
   }
}

/// @brief Read simulation checkpoint from disk.
///
/// Each MPI rank reads its own file.  Returns false if the file
/// does not exist or cannot be parsed.
///
/// @param[in]  prefix     Checkpoint file prefix
/// @param[out] t          Simulation time [s]
/// @param[out] dt         Adaptive time step [s]
/// @param[out] step       Step counter
/// @param[out] num_eq     Seismic event counter
/// @param[out] in_eq      In-seismic-event flag
/// @param[out] state      ODE state vector
/// @param[out] displacement  Domain displacement
/// @param[out] traction   Shear stress at fault
/// @param[out] slip_rate  Slip rate at fault
/// @param[out] fsal_init  FSAL initialized flag
/// @param[out] k0         FSAL stage vector
/// @param[in]  mpi        MPI context (may be nullptr for serial)
/// @return true if checkpoint was loaded successfully
inline bool ReadCheckpoint(const std::string &prefix,
                           real_t &t, real_t &dt,
                           int &step, int &num_eq, bool &in_eq,
                           Vector &state,
                           Vector &displacement,
                           Vector &traction,
                           Vector &slip_rate,
                           bool &fsal_init, Vector &k0,
                           const MPIContext *mpi)
{
   int rank = mpi ? mpi->Rank() : 0;
   int size = mpi ? mpi->Size() : 1;

   std::string filename = CheckpointFilename(prefix, rank);
   std::ifstream in(filename);
   if (!in.good()) { return false; }

   auto read_tag = [&](const std::string &expected_tag)
   {
      std::string tag;
      in >> tag;
      MFEM_VERIFY(tag == expected_tag,
                   "Checkpoint parse error: expected '" << expected_tag
                   << "', got '" << tag << "' in " << filename);
   };

   auto read_vector = [&](const std::string &size_tag, Vector &vec)
   {
      read_tag(size_tag);
      int n;
      in >> n;
      vec.SetSize(n);
      for (int i = 0; i < n; i++)
      {
         in >> vec(i);
      }
   };

   // Header
   read_tag("SEAS_CHECKPOINT_V1");

   int file_num_ranks, file_rank;
   read_tag("num_ranks"); in >> file_num_ranks;
   read_tag("rank");      in >> file_rank;

   MFEM_VERIFY(file_num_ranks == size,
               "Checkpoint num_ranks mismatch: file has " << file_num_ranks
               << " but running with " << size);
   MFEM_VERIFY(file_rank == rank,
               "Checkpoint rank mismatch: file has " << file_rank
               << " but this is rank " << rank);

   read_tag("time"); in >> t;
   read_tag("dt");   in >> dt;
   read_tag("step"); in >> step;

   read_tag("num_seismic_events"); in >> num_eq;
   int in_eq_int;
   read_tag("in_seismic_event"); in >> in_eq_int;
   in_eq = (in_eq_int != 0);

   // Vectors
   read_vector("state_size", state);
   read_vector("displacement_size", displacement);
   read_vector("traction_size", traction);
   read_vector("slip_rate_size", slip_rate);

   // FSAL
   int fsal_int;
   read_tag("fsal_initialized"); in >> fsal_int;
   fsal_init = (fsal_int != 0);
   read_vector("k0_size", k0);

   in.close();

   if (mpi) { mpi->Barrier(); }

   if (!mpi || mpi->IsRoot())
   {
      mfem::out << "Checkpoint loaded: " << prefix
                << " (t=" << t << " s, step=" << step << ")\n";
   }

   return true;
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CHECKPOINT_HPP

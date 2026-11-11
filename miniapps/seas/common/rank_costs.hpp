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

// rank_costs.hpp: where each rank spends a dynamic step (speed-up options,
// D0 design rev1, 2026-11-11).
//
// Why: "time per dynamic step, rank 0 / slowest rank" times dyn.Step() only,
// and a rank's time there is its own work plus its wait for the others at
// the step's reductions and halo exchanges.  To weight the mesh cut by the
// real cost of a band element and a fault element, the work and the wait of
// every rank must be told apart.
//
// What: wall-clock sums per rank, in named slots.  A "wait" slot is a
// blocking MPI call (the rank idles there until its partners arrive); a
// "busy" slot is a block of the rank's own work (its total, any wait inside
// it is subtracted by the report).
//
//   WSync    wait  the per-step reduction at the end of the fault advance
//   WPml     wait  the band's halo exchange (PmlDynamicTerms::Gather_)
//   WStiff   wait  MPI_Waitall of the stiffness halo (DenseBlockStiffness)
//   WFault   wait  the fault faces' neighbour exchange (FaultNbrExchange_)
//   WOut     wait  the "is a row due" reductions of the output classes
//   Pml      busy  PML Advance + AddForce (contains WPml)
//   Stiff    busy  the dense stiffness product (contains WStiff)
//   Rhs      busy  the Dirichlet wall loads of K and K_eta
//   Damp     busy  the damping blocks (absorbing faces, band C_d, quench)
//   Mass     busy  the inverse mass
//
// It only measures.  Off (the default) no clock is read; on or off no
// result, file or screen line of the run changes except the report itself.
// fd_cycle --rank-costs turns it on and prints the table.

#ifndef MFEM_SEAS_RANK_COSTS_HPP
#define MFEM_SEAS_RANK_COSTS_HPP

#include "mfem.hpp"

namespace mfem
{
namespace seas
{

struct RankCosts
{
   enum Slot { WSync = 0, WPml, WStiff, WFault, WOut, Pml, Stiff, Rhs, Damp, Mass, NumSlots };

   static bool &On()
   {
      static bool on = false;
      return on;
   }
   /// The sums [s] of this rank, by Slot.
   static double *Sum()
   {
      static double s[NumSlots] = {};
      return s;
   }
   static void Reset()
   {
      for (int i = 0; i < NumSlots; i++) { Sum()[i] = 0.0; }
   }

   /// Adds the wall time from the constructor to the destructor to `slot`.
   class Scope
   {
   public:
      explicit Scope(Slot slot) : slot_(slot), t0_(On() ? MPI_Wtime() : 0.0) { }
      ~Scope() { if (On()) { Sum()[slot_] += MPI_Wtime() - t0_; } }
      Scope(const Scope &) = delete;
      Scope &operator=(const Scope &) = delete;
   private:
      Slot slot_;
      double t0_;
   };
};

} // namespace seas
} // namespace mfem

#endif

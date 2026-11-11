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

// Timers for the QD loop's work outside the TS step (speed-up options note
// rev1, Q0; 2026-11-10).
//
// Why: in the profile of the sjf1 implicit QD run (PETSc -log_view, Mac,
// 494e2092) 116 s of the 303 s loop lay outside TSStep, in no named event.
// QDTimer puts each piece of that work in a PETSc log event of its own, so
// the same -log_view names it:
//
//   QDPostState    option I: the accepted U into State() and the gap
//   QDPostKernel   option I: the kernel evaluation at the accepted state
//   QDMonitor      the whole monitor call of one step; inside it
//   QDMonScan        the NaN scan and its reduction
//   QDMonRows        events, station and fault rows, the global row
//   QDMonSnaps       snapshots and the regime test
//   QDMonCkpt        checkpoints and the wall-clock test
//
// It only measures: no result, file or screen line changes.  PETSc logs an
// event only when logging is on (-log_view in the options file); otherwise
// Begin and End return at once.  Without PETSc the class does nothing.

#ifndef MFEM_SEAS_CYCLE_QD_TIMERS_HPP
#define MFEM_SEAS_CYCLE_QD_TIMERS_HPP

#include "mfem.hpp"

#ifdef MFEM_USE_PETSC
#include "petsc.h"
#include <map>
#include <string>
#endif

namespace mfem
{
namespace seas
{
namespace cycle
{

/// One PETSc log event from the constructor to the destructor (or Stop()).
/// Input: the event's name (one of the list above; registered on first use).
class QDTimer
{
public:
#ifdef MFEM_USE_PETSC
   explicit QDTimer(const char *name) : id_(Id_(name))
   {
      MFEM_VERIFY(PetscLogEventBegin(id_, 0, 0, 0, 0) == PETSC_SUCCESS,
                  "PetscLogEventBegin failed");
   }
   ~QDTimer() { Stop(); }
   void Stop()
   {
      if (!on_) { return; }
      on_ = false;
      MFEM_VERIFY(PetscLogEventEnd(id_, 0, 0, 0, 0) == PETSC_SUCCESS,
                  "PetscLogEventEnd failed");
   }

private:
   static PetscLogEvent Id_(const char *name)
   {
      static std::map<std::string, PetscLogEvent> ids;
      auto it = ids.find(name);
      if (it != ids.end()) { return it->second; }
      PetscLogEvent e = 0;
      MFEM_VERIFY(PetscLogEventRegister(name, TS_CLASSID, &e) == PETSC_SUCCESS,
                  "PetscLogEventRegister failed");
      ids[name] = e;
      return e;
   }
   PetscLogEvent id_;
   bool on_ = true;
#else
   explicit QDTimer(const char *) {}
   void Stop() {}
#endif
   QDTimer(const QDTimer &) = delete;
   QDTimer &operator=(const QDTimer &) = delete;
};

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif

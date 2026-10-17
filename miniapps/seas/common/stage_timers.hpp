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

// stage_timers.hpp — cumulative wall-clock timers for the parts of one
// quasi-dynamic stage (SEASQuasiDynamicOperator::Mult and the domain solve it
// calls) and of the accepted-step monitor.  Enabled by seas_bp5_full
// --stage-timers N (GPU plan G0: where a stage spends its time on the CPU).
//
// Timers are rank-local and never synchronize across ranks; the driver
// reduces them when it prints.  On a GPU the driver sets a device sync
// (SetDeviceSync), called at every start and stop, so asynchronous kernels
// are charged to the part that launched them.  An operator holding a null
// StageTimers pointer does no timing at all.
//
// GPU plan G4: in a CUDA build each timed part is also an NVTX range named
// "seas:<part>", so an Nsight Systems trace shows which part made each
// host-device copy (jobs/vista/g4_trace_check.py).  NVTX v3 is header-only
// and costs next to nothing when no profiler is attached.

#ifndef MFEM_SEAS_STAGE_TIMERS_HPP
#define MFEM_SEAS_STAGE_TIMERS_HPP

#include "mfem.hpp"

#include <array>
#include <chrono>

#if defined(MFEM_USE_CUDA) && defined(__has_include)
#if __has_include(<nvtx3/nvToolsExt.h>)
#include <nvtx3/nvToolsExt.h>
#define SEAS_HAVE_NVTX 1
#endif
#endif

namespace mfem
{
namespace seas
{

class StageTimers
{
public:
   enum Part
   {
      MatrixSetup = 0,   ///< K assembly + solver setup/factorization (once)
      SlipExpand,        ///< owned slip -> local fault slip (halo exchange)
      SlipRHS,           ///< RHS vector + fault-slip face terms
      DirichletRHS,      ///< Dirichlet loading face terms
      Solve,             ///< linear solve K u = b
      SolvePost,         ///< residual/NaN checks and copy to displacement
      Traction,          ///< fault traction (shear and normal) + restriction
      Friction,          ///< fault RHS: friction solve and state rate
      StageTotal,        ///< whole SEASQuasiDynamicOperator::Mult
      MonitorTotal,      ///< whole accepted-step monitor
      MonitorOutput,     ///< station, ParaView and global output
      MonitorCheckpoint, ///< checkpoint writes
      NumParts
   };

   static const char *Name(int p)
   {
      static const char *names[NumParts] =
      {
         "matrix_setup", "slip_expand", "slip_rhs", "dirichlet_rhs",
         "solve", "solve_post", "traction", "friction", "stage_total",
         "monitor_total", "monitor_output", "monitor_checkpoint"
      };
      return names[p];
   }

   /// NVTX range name of part @a p: "seas:" + Name(p).
   static const char *NvtxName(int p)
   {
      static const char *names[NumParts] =
      {
         "seas:matrix_setup", "seas:slip_expand", "seas:slip_rhs",
         "seas:dirichlet_rhs", "seas:solve", "seas:solve_post",
         "seas:traction", "seas:friction", "seas:stage_total",
         "seas:monitor_total", "seas:monitor_output", "seas:monitor_checkpoint"
      };
      return names[p];
   }

   /// Wait for the device at every Start and Stop (nullptr: never).
   void SetDeviceSync(void (*sync)()) { sync_ = sync; }

   void Start(Part p)
   {
      if (sync_) { sync_(); }
#ifdef SEAS_HAVE_NVTX
      nvtx_[p] = nvtxRangeStartA(NvtxName(p));
#endif
      start_[p] = Now();
   }

   void Stop(Part p)
   {
      if (sync_) { sync_(); }
      seconds_[p] += Now() - start_[p];
      calls_[p]++;
#ifdef SEAS_HAVE_NVTX
      nvtxRangeEnd(nvtx_[p]);
#endif
   }

   double Seconds(int p) const { return seconds_[p]; }
   long long Calls(int p) const { return calls_[p]; }

   /// Zero every part except MatrixSetup, which happens once at start-up.
   void ResetStageParts()
   {
      for (int p = 0; p < NumParts; p++)
      {
         if (p == MatrixSetup) { continue; }
         seconds_[p] = 0.0;
         calls_[p] = 0;
      }
   }

private:
   static double Now()
   {
      using clock = std::chrono::steady_clock;
      return std::chrono::duration<double>(
                clock::now().time_since_epoch()).count();
   }

   void (*sync_)() = nullptr;
   std::array<double, NumParts> start_{};
   std::array<double, NumParts> seconds_{};
   std::array<long long, NumParts> calls_{};
#ifdef SEAS_HAVE_NVTX
   std::array<nvtxRangeId_t, NumParts> nvtx_{};
#endif
};

/// RAII helper: times one part when @a timers is non-null.
class ScopedStageTimer
{
public:
   ScopedStageTimer(StageTimers *timers, StageTimers::Part part)
      : timers_(timers), part_(part)
   {
      if (timers_) { timers_->Start(part_); }
   }
   ~ScopedStageTimer() { if (timers_) { timers_->Stop(part_); } }

   ScopedStageTimer(const ScopedStageTimer &) = delete;
   ScopedStageTimer &operator=(const ScopedStageTimer &) = delete;

private:
   StageTimers *timers_;
   StageTimers::Part part_;
};

/// RAII helper: an NVTX range around a span that is not a timed part (GPU
/// plan G4); nothing outside a CUDA build.
class NvtxRange
{
public:
   explicit NvtxRange(const char *name)
   {
#ifdef SEAS_HAVE_NVTX
      id_ = nvtxRangeStartA(name);
#else
      (void) name;
#endif
   }
   ~NvtxRange()
   {
#ifdef SEAS_HAVE_NVTX
      nvtxRangeEnd(id_);
#endif
   }

   NvtxRange(const NvtxRange &) = delete;
   NvtxRange &operator=(const NvtxRange &) = delete;

#ifdef SEAS_HAVE_NVTX
private:
   nvtxRangeId_t id_;
#endif
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_STAGE_TIMERS_HPP

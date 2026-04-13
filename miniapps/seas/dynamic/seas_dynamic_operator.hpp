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

#ifndef MFEM_SEAS_DYNAMIC_OPERATOR_HPP
#define MFEM_SEAS_DYNAMIC_OPERATOR_HPP

#include "mfem.hpp"
#include "wave_operator.hpp"

namespace mfem
{
namespace seas
{

// Forward declaration for Phase 3 (fault-face Riemann solver).
class FaultFaceFlux;

/// @brief Coupling operator for dynamic rupture simulations.
///
/// Wraps WaveOperator* + FaultFaceFlux* (analogous to SEASQuasiDynamicOperator
/// which wraps DomainOperator* + RateStateFault*). Provides the coupled
/// Mult() that the time integrator calls.
///
/// Architecture (R-001 fix): Composition, not inheritance. WaveOperator does
/// not inherit from DomainOperator. SEASDynamicOperator holds non-owning
/// pointers to both WaveOperator and FaultFaceFlux.
///
/// For hybrid QD+dynamic (Phase 5), SEASHybridOperator switches between
/// SEASQuasiDynamicOperator* and SEASDynamicOperator* — trivial pointer swap.
class SEASDynamicOperator : public TimeDependentOperator
{
public:
   /// @brief Construct from a WaveOperator and optional FaultFaceFlux.
   ///
   /// @param[in] wave  Non-owning pointer to the wave operator.
   /// @param[in] fault  Non-owning pointer to the fault flux (nullptr if no fault).
   SEASDynamicOperator(WaveOperator *wave, FaultFaceFlux *fault = nullptr)
      : TimeDependentOperator(wave ? wave->Height() : 0),
        wave_(wave), fault_(fault)
   {
      MFEM_VERIFY(wave, "WaveOperator must not be null");
   }

   /// @brief Compute dQ/dt by delegating to the wave operator.
   ///
   /// If a FaultFaceFlux is set, the wave operator's face flux loop
   /// dispatches fault faces to it (Phase 3). For now (Phase 1),
   /// this simply calls wave_->Mult().
   void Mult(const Vector &Q, Vector &dQdt) const override
   {
      wave_->Mult(Q, dQdt);
   }

   /// Set the time for the operator and its sub-components.
   void SetTime(const real_t t_) override
   {
      TimeDependentOperator::SetTime(t_);
      wave_->SetTime(t_);
   }

   /// Access the wave operator.
   WaveOperator *GetWaveOperator() { return wave_; }
   const WaveOperator *GetWaveOperator() const { return wave_; }

   /// Access the fault flux (may be null in Phase 1).
   FaultFaceFlux *GetFaultFlux() { return fault_; }
   const FaultFaceFlux *GetFaultFlux() const { return fault_; }

   /// Set/replace the fault flux (for Phase 3 integration).
   void SetFaultFlux(FaultFaceFlux *fault) { fault_ = fault; }

private:
   WaveOperator *wave_;       ///< Non-owning
   FaultFaceFlux *fault_;     ///< Non-owning (null until Phase 3)
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_DYNAMIC_OPERATOR_HPP

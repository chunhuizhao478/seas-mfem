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
/// Wraps WaveOperator<MeshType>* + FaultFaceFlux*. Provides the coupled
/// Mult() that the time integrator calls.
///
/// Templated on MeshType to match WaveOperator<MeshType>.
template <typename MeshType = Mesh>
class SEASDynamicOperator : public TimeDependentOperator
{
public:
   SEASDynamicOperator(WaveOperator<MeshType> *wave,
                       FaultFaceFlux *fault = nullptr)
      : TimeDependentOperator(wave ? wave->Height() : 0),
        wave_(wave), fault_(fault)
   {
      MFEM_VERIFY(wave, "WaveOperator must not be null");
   }

   void Mult(const Vector &Q, Vector &dQdt) const override
   {
      wave_->Mult(Q, dQdt);
   }

   void SetTime(const real_t t_) override
   {
      TimeDependentOperator::SetTime(t_);
      wave_->SetTime(t_);
   }

   WaveOperator<MeshType> *GetWaveOperator() { return wave_; }
   const WaveOperator<MeshType> *GetWaveOperator() const { return wave_; }

   FaultFaceFlux *GetFaultFlux() { return fault_; }
   const FaultFaceFlux *GetFaultFlux() const { return fault_; }

   void SetFaultFlux(FaultFaceFlux *fault) { fault_ = fault; }

private:
   WaveOperator<MeshType> *wave_;
   FaultFaceFlux *fault_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_DYNAMIC_OPERATOR_HPP

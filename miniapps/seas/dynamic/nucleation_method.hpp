// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/nucleation_method.hpp — Phase 7 of
// PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.
//
// The method-oriented nucleation strategy interface (`INucleationMethod`,
// plan §5.2) and its four concrete methods, all benchmark-agnostic:
//
//   GaussianGradualOverstress           — SAFS + TPV option (per-sub-step)
//   CompactCircularGradualOverstress    — TPV102/104 SCEC bell (per-sub-step)
//   InstantaneousOverstressCircular     — TPV31 one-shot patch (ApplyOnce)
//   StaticOverstress                    — TPV205 (both hooks no-op)
//
// The driver builds one of these via `MakeNucleation` (nucleation_factory.hpp)
// and calls `ApplyIncrement` once per ADER sub-step (the `nuc_cb` body) and
// `ApplyOnce` once after init (for the one-shot kinds).  The concretes are
// thin wrappers over the resolved per-DOF params + the apply free functions in
// spatial_nucleation.{hpp,cpp}; each exposes its resolved params via `Params()`
// so the driver's diagnostics / ParaView fields can read them.

#ifndef MFEM_SEAS_NUCLEATION_METHOD_HPP
#define MFEM_SEAS_NUCLEATION_METHOD_HPP

#include "mfem.hpp"

#include "fault_face_flux.hpp"     // DOFData
#include "spatial_nucleation.hpp"  // resolved per-DOF param types + apply fns

#include <utility>
#include <vector>

namespace mfem
{
namespace seas
{

/// @brief Nucleation strategy interface (plan §5.2).  `ApplyIncrement` is the
/// per-sub-step `nuc_cb` body the driver passes into
/// `IFrictionIterator::Advance`; `ApplyOnce` seeds a one-shot patch after init
/// (default no-op for the per-sub-step kinds).  `IsPerSubStep()` lets the
/// driver decide whether the per-sub-step callback does any work.
class INucleationMethod
{
public:
   virtual void ApplyIncrement(std::vector<DOFData>& dof,
                               real_t t, real_t dt) = 0;
   /// Phase 14.3 — ABSOLUTE forcing at stage time `t` (set, not accumulate).
   /// The RK stepper calls this at the start of each RK stage with
   /// `t = t_step_start + c_i·dt`; the ADER path keeps using `ApplyIncrement`
   /// via `nuc_cb` (unchanged).  Setting the absolute SCEC-ramp value avoids
   /// the double-apply that the telescoped per-sub-step increment would incur
   /// across RK stages revisiting a sub-interval.  Gradual kinds set
   /// `SmoothStep(t)·F(r)·Δτ`; static / one-shot kinds set their
   /// time-independent value and ignore `t`.
   virtual void ApplyAbsolute(std::vector<DOFData>& dof, real_t t) = 0;
   virtual void ApplyOnce(std::vector<DOFData>& dof) { (void)dof; }
   virtual bool IsPerSubStep() const = 0;
   virtual ~INucleationMethod() = default;
};

/// @brief TPV205 / no-[nucleation] case: pre-stress (and any patches) already
/// live in `tau_pre_`; both hooks no-op.
class StaticOverstress : public INucleationMethod
{
public:
   void ApplyIncrement(std::vector<DOFData>& dof, real_t t, real_t dt) override
   { (void)dof; (void)t; (void)dt; }
   /// Pre-stress already lives in `tau_pre_`; nothing time-varying to set.
   void ApplyAbsolute(std::vector<DOFData>& dof, real_t t) override
   { (void)dof; (void)t; }
   bool IsPerSubStep() const override { return false; }
};

/// @brief SAFS / TPV Gaussian gradual overstress.  Wraps the existing
/// `ResolveGradualOverstress` / `ApplyGradualOverstressIncrement` (byte-
/// identical to the driver's pre-Phase-7 inline call).
class GaussianGradualOverstress : public INucleationMethod
{
public:
   GaussianGradualOverstress(spatial::GradualOverstressPerDOFParams params,
                             real_t T_nuc_s)
      : params_(std::move(params)), T_nuc_s_(T_nuc_s) {}

   void ApplyIncrement(std::vector<DOFData>& dof, real_t t, real_t dt) override
   {
      spatial::ApplyGradualOverstressIncrement(dof, params_, T_nuc_s_, t, dt);
   }
   void ApplyAbsolute(std::vector<DOFData>& dof, real_t t) override
   {
      spatial::ApplyGradualOverstressAbsolute(dof, params_, T_nuc_s_, t);
   }
   bool IsPerSubStep() const override { return true; }

   /// Resolved per-DOF amplitudes/radial (for driver diagnostics + ParaView).
   const spatial::GradualOverstressPerDOFParams& Params() const
   { return params_; }

private:
   spatial::GradualOverstressPerDOFParams params_;
   real_t                                 T_nuc_s_;
};

/// @brief TPV102/104 compact-circular gradual overstress (SCEC bell).  Wraps
/// `ResolveGradualOverstressCompactCircular` /
/// `ApplyGradualOverstressCompactCircularIncrement`.
class CompactCircularGradualOverstress : public INucleationMethod
{
public:
   CompactCircularGradualOverstress(
      spatial::CompactCircularPerDOFParams params, real_t T_nuc_s)
      : params_(std::move(params)), T_nuc_s_(T_nuc_s) {}

   void ApplyIncrement(std::vector<DOFData>& dof, real_t t, real_t dt) override
   {
      spatial::ApplyGradualOverstressCompactCircularIncrement(
         dof, params_, T_nuc_s_, t, dt);
   }
   void ApplyAbsolute(std::vector<DOFData>& dof, real_t t) override
   {
      spatial::ApplyGradualOverstressCompactCircularAbsolute(
         dof, params_, T_nuc_s_, t);
   }
   bool IsPerSubStep() const override { return true; }

   const spatial::CompactCircularPerDOFParams& Params() const
   { return params_; }

private:
   spatial::CompactCircularPerDOFParams params_;
   real_t                               T_nuc_s_;
};

/// @brief TPV31 instantaneous circular overstress (cosine-tapered, one-shot).
/// `ApplyOnce` seeds the strike component (`tau2_nuc`) once at t=0;
/// `ApplyIncrement` is a no-op; `IsPerSubStep() == false`.
class InstantaneousOverstressCircular : public INucleationMethod
{
public:
   explicit InstantaneousOverstressCircular(
      spatial::InstantaneousOverstressPerDOFParams params)
      : params_(std::move(params)) {}

   void ApplyIncrement(std::vector<DOFData>& dof, real_t t, real_t dt) override
   { (void)dof; (void)t; (void)dt; }   // one-shot: no per-sub-step work

   /// Phase 14.3 — instantaneous: the patch is full from t=0, so SET the
   /// absolute value (idempotent, ignores t).  Unlike `ApplyOnce` (which `+=`s
   /// once at init), this assigns `tau2_nuc = amplitude_strike(i)` so repeated
   /// per-stage RK calls do not re-seed.  No-op when disabled / no fault DOFs.
   void ApplyAbsolute(std::vector<DOFData>& dof, real_t t) override
   {
      (void)t;
      if (params_.amplitude_strike.Size() == 0) { return; }
      const int n = static_cast<int>(dof.size());
      MFEM_VERIFY(params_.amplitude_strike.Size() == n,
                  "InstantaneousOverstressCircular::ApplyAbsolute: "
                  "amplitude_strike size ("
                  << params_.amplitude_strike.Size()
                  << ") != dof.size() (" << n << ")");
      for (int i = 0; i < n; ++i)
      {
         dof[i].tau2_nuc = params_.amplitude_strike(i);
      }
   }

   /// Seed the patch ONCE.  Adds `amplitude_strike(i)` to `tau2_nuc` (pure
   /// strike-slip; tau1_nuc / sigma_n_nuc untouched).  No-op when disabled
   /// (zero-sized params) or this rank has no fault DOFs.  The driver must
   /// call this exactly once after initialization.
   void ApplyOnce(std::vector<DOFData>& dof) override
   {
      if (params_.amplitude_strike.Size() == 0) { return; }
      const int n = static_cast<int>(dof.size());
      MFEM_VERIFY(params_.amplitude_strike.Size() == n,
                  "InstantaneousOverstressCircular::ApplyOnce: "
                  "amplitude_strike size ("
                  << params_.amplitude_strike.Size()
                  << ") != dof.size() (" << n << ")");
      for (int i = 0; i < n; ++i)
      {
         dof[i].tau2_nuc += params_.amplitude_strike(i);
      }
   }
   bool IsPerSubStep() const override { return false; }

   const spatial::InstantaneousOverstressPerDOFParams& Params() const
   { return params_; }

private:
   spatial::InstantaneousOverstressPerDOFParams params_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_NUCLEATION_METHOD_HPP

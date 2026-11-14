// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// normal_stress_relaxation.hpp — opt-in sigma_n' relaxation of the friction
// law ([cycle.relaxation] t_relax_yr; Chunhui Zhao, 2026-11-02: "implement
// the relaxation, test and verify it, and apply to UCERF 2D plane strain,
// with 10MPa cap").  Design rev 0:
// /mnt/project-files/2d-inplane/relaxation/relaxation-design-rev0.md.
//
// The law at every owned fault DOF, in both regimes:
//
//    sigma_law = max(sigma0 + dsig - r, sigma_floor)
//    dr/dt     = (dsig - r) / t_r                 (QD regime; r = 0 at t = 0)
//
//    sigma0  prestress sigma_n' (compression positive)
//    dsig    elastic normal traction change of the QD solve (slip + load),
//            the kernel's normal_traction
//    r       the relaxed part of dsig [Pa]
//    t_r     relaxation time [s]
//
// Above the floor, d(sigma_law)/dt = d(dsig)/dt - (sigma_law - sigma0)/t_r:
// Ozawa, Ando & Dunham (2023, EPSL 617, 118276), eq. 6, for the normal
// stress only (their shear relaxation is not done: our tau0 / sigma_n' of
// 0.29-0.57 would pull tau down about as fast as the walls load it).  The
// bulk keeps the elastic stress: r changes only the strength (as the floor
// and the normal-stress hook H-2D1 do).
//
// QD time stepping (r is not a TS unknown).  Between accepted steps n and
// n + 1 every stage at time t sees
//
//    red(t) = r_n + (dsig_n - r_n) (1 - exp(-(t - t_n) / t_r))
//
// (exact for dsig frozen at its step-start value), and after the accepted
// step, at t_{n+1}:  r_{n+1} = red(t_{n+1}),  dsig_{n+1} = the accepted
// state's normal traction,  t_n = t_{n+1}.  The stage that ends step n and
// starts step n + 1 (FSAL) therefore sees the r of the new step's start,
// and a stage at t_n sees r_n exactly (exp(-0) = 1), so a checkpoint and a
// restart repeat the run bit for bit.  Error: of order (t_{n+1} - t_n) /
// t_r of (dsig - r) per step.
//
// Dynamic regime: r stays as it was at the switch (phases last seconds to
// an hour, t_r thousands of years).  FrozenNormalStressReduction hands it
// to the Nitsche operator's normal-stress hook (H-2D1): the law sees
// sigma_eff - r, then the floor.  The switch back keeps r and restarts the
// clock at the re-entry time (Restart), with dsig from the re-entry solve.
//
// [cycle.relaxation] opening = "relaxed" (Chunhui Zhao, 2026-11-13: "bring
// the relaxation review"; with [cycle.qd] integrator = "esdirk" only): the
// relaxed stress sigma0 + dsig - r is the fault's normal stress for the
// opening rule too, not only for the strength.  A point opens where the
// relaxed stress would go below zero, and at an open point the relaxed
// stress is 0 (the elastic one is r).  r keeps following dsig at open
// points (nothing is frozen), so an open point's gap closes again over t_r
// once the unclamping stops.  Why: on a non-planar trace the uniform plate
// loading changes the elastic normal stress of every locked section at a
// constant rate set by its strike (UCERF 2.5D, Vista 1142230: -0.7 MPa per
// kyr at strike 142 deg, +0.75 at 90-105 deg, slip direction 123 deg), so
// with the opening rule on the elastic stress every section turned one way
// opens after sigma0 / rate (52 MPa: about 73 kyr), whatever the strength
// law sees.  With the rule on the relaxed stress the change is bounded by
// rate * t_r.  The default "elastic" is the rule above, every bit unchanged.
//
// No solver header (rule M4).

#ifndef MFEM_SEAS_CYCLE_NORMAL_STRESS_RELAXATION_HPP
#define MFEM_SEAS_CYCLE_NORMAL_STRESS_RELAXATION_HPP

#include "mfem.hpp"
#include "../fault/fault_normal_stress_modifier.hpp"

#include <cmath>
#include <ostream>

namespace mfem
{
namespace seas
{
namespace cycle
{

class NormalStressRelaxation
{
public:
   /// n owned fault DOFs; t_relax [s] > 0.  r = dsig = 0 at t_n = 0.
   NormalStressRelaxation(int n, real_t t_relax) : t_r_(t_relax)
   {
      MFEM_VERIFY(n >= 0, "NormalStressRelaxation: " << n << " DOFs");
      MFEM_VERIFY(std::isfinite(t_relax) && t_relax > 0.0,
                  "NormalStressRelaxation: t_relax = " << t_relax << " s (must be > 0)");
      r_.SetSize(n);
      r_ = 0.0;
      d_.SetSize(n);
      d_ = 0.0;
      red_.SetSize(n);
      red_ = 0.0;
   }

   int Size() const { return r_.Size(); }
   real_t TRelax() const { return t_r_; }
   /// The time of r_n and dsig_n (the last accepted step, or the restart).
   real_t Time() const { return t_n_; }
   /// r_n [Pa], the value at Time().
   const Vector &R() const { return r_; }
   /// dsig_n [Pa], the elastic normal traction change at Time().
   const Vector &Delta() const { return d_; }
   /// What the QD kernel subtracts (RateStateFaultOperator::SetSigmaNReduction):
   /// red at the time of the last Fill.
   const Vector &Reduction() const { return red_; }

   /// red = r_n + (dsig_n - r_n) (1 - exp(-(t - t_n) / t_r)) (the file
   /// comment); t >= t_n.
   void Fill(real_t t)
   {
      const real_t E = -std::expm1(-(t - t_n_) / t_r_);
      const int n = r_.Size();
      const real_t *r = r_.HostRead();
      const real_t *d = d_.HostRead();
      real_t *red = red_.HostWrite();
      for (int i = 0; i < n; i++) { red[i] = r[i] + (d[i] - r[i]) * E; }
      if (frozen_)
      {
         // [cycle.opening] decision D5: r frozen at open points
         const real_t *f = frozen_->HostRead();
         for (int i = 0; i < n; i++) { if (f[i] != 0.0) { red[i] = r[i]; } }
      }
   }

   /// red of one DOF at time t >= t_n, the same expression as Fill's (and
   /// its frozen rule), without touching Reduction().  For the stage solves
   /// of [cycle.qd] integrator = "esdirk" (cycle/qd_implicit.hpp), which
   /// ask per point and stage time.
   real_t ReductionAt(real_t t, int i) const
   {
      const real_t r = r_(i);
      if (frozen_ && (*frozen_)(i) != 0.0) { return r; }
      const real_t E = -std::expm1(-(t - t_n_) / t_r_);
      return r + (d_(i) - r) * E;
   }

   /// [cycle.opening] (fault opening design rev 1.1, decision D5): r does not
   /// change at a DOF with mask(i) != 0 (an open point), so the law resumes
   /// with the r it had when the point opened.  The caller owns the mask.
   /// Null (the default): every bit unchanged.
   void SetFrozen(const Vector *mask)
   {
      MFEM_VERIFY(mask == nullptr || mask->Size() == r_.Size(),
                  "NormalStressRelaxation::SetFrozen: " << mask->Size() << " values for "
                  << r_.Size() << " DOFs");
      frozen_ = mask;
   }

   /// After an accepted QD step at t: r = red(t), dsig = `dsig` (the
   /// accepted state's elastic normal traction), t_n = t.
   void Accept(real_t t, const Vector &dsig)
   {
      Fill(t);
      r_ = red_;
      SetDelta(dsig);
      t_n_ = t;
   }

   /// Restart the clock at t with r unchanged (the switch back): dsig = r
   /// until SetDelta, so red = r at any t until then.
   void Restart(real_t t)
   {
      t_n_ = t;
      d_ = r_;
      red_ = r_;
   }

   /// dsig_n from a solve at t_n (the start, the switch, the re-entry).
   void SetDelta(const Vector &dsig)
   {
      MFEM_VERIFY(dsig.Size() == d_.Size(), "NormalStressRelaxation::SetDelta: "
                  << dsig.Size() << " values for " << d_.Size() << " DOFs");
      dsig.HostRead();
      d_ = dsig;
   }

   /// A checkpoint's state (restart): t_n, r_n, dsig_n; red = r_n.
   void SetState(real_t t_n, const Vector &r, const Vector &dsig)
   {
      MFEM_VERIFY(r.Size() == r_.Size() && dsig.Size() == d_.Size(),
                  "NormalStressRelaxation::SetState: " << r.Size() << " and " << dsig.Size()
                  << " values for " << r_.Size() << " DOFs");
      t_n_ = t_n;
      r_ = r;
      d_ = dsig;
      red_ = r_;
   }

private:
   real_t t_r_;
   real_t t_n_ = 0.0;
   Vector r_, d_, red_;
   const Vector *frozen_ = nullptr;   ///< SetFrozen
};

/// The dynamic regime's view: the hook H-2D1 takes r off sigma_eff at every
/// owned DOF; r does not change during a phase (Advance does nothing).  The
/// operator's owned fault DOFs are the QD table's, point for point
/// (dyn_phase.hpp), so r is indexed the same way.
class FrozenNormalStressReduction : public IFaultNormalStressModifier
{
public:
   explicit FrozenNormalStressReduction(const Vector &r) : r_(r) {}
   int NumPoints() const override { return r_.Size(); }
   real_t NormalStressReduction(int j) const override { return r_(j); }
   void Advance(real_t, const real_t *, const real_t *) override {}

private:
   const Vector &r_;
};

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_NORMAL_STRESS_RELAXATION_HPP

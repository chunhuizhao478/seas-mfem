// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// regime_controller.hpp — fully dynamic BP5 cycle (plan rev 3, Section
// 3.3 g): RegimeController, which decides when the driver switches regime.
// In: the global max slip rate and the dynamic phase's own clock tau.  Out:
// the next regime, or a stop.
//
//   Start(Vmax):      Vmax > V_on ? DYN : QD
//   GoDynamic(Vmax):  Vmax > V_on                     (after each QD step)
//   EnterDynamic(t):  the phase starts at time t
//   AfterDynStep(Vmax, t, tau):
//     below = (Vmax < V_off) ? t - t_above : 0,  t_above = time of the last
//             step with Vmax >= V_off, or of the phase start
//     Vmax < V_off, below >= T_dwell and tau >= T_min  ->  END_OF_PHASE (back to QD)
//     tau >= T_max                       ->  TOO_LONG (stop and report)
//     else                               ->  STAY
//
// "below" is the time max V has stayed below V_off, so the phase ends
// T_dwell after the last step with max V >= V_off.  It is the plan's
// `below = (Vmax < V_off) ? below + dt : 0` written as a difference of the
// driver's t: FdCycleOutput::UpdateEvents measures the catalog's dwell the
// same way from the same t, so with event_V = V_off and event_dwell =
// T_dwell the catalog closes the event on the step the phase ends.  T_max
// stops the run instead of forcing a switch (plan R9).  Defaults in
// SwitchConfig.
//
// The switch probe of plan P3 runs one controller per (V_off, T_dwell) pair
// on the same max V history; a checkpoint keeps Phases() (Restore).
//
// No solver header (rule M4).

#ifndef MFEM_SEAS_CYCLE_REGIME_CONTROLLER_HPP
#define MFEM_SEAS_CYCLE_REGIME_CONTROLLER_HPP

#include "mfem.hpp"
#include "cycle_config.hpp"

namespace mfem
{
namespace seas
{
namespace cycle
{

enum class Regime { QD, Dynamic };

enum class DynVerdict { Stay, EndOfPhase, TooLong };

inline const char *RegimeName(Regime r) { return r == Regime::QD ? "QD" : "dynamic"; }

class RegimeController
{
public:
   explicit RegimeController(const SwitchConfig &cfg) : cfg_(cfg) {}

   Regime Start(real_t Vmax) const { return GoDynamic(Vmax) ? Regime::Dynamic : Regime::QD; }
   bool GoDynamic(real_t Vmax) const { return Vmax > cfg_.V_on; }

   /// A dynamic phase starts at time t: the dwell counter restarts.
   void EnterDynamic(real_t t)
   {
      t_above_ = t;
      below_ = 0.0;
      phases_++;
   }

   /// After each dynamic step, at time t and phase time tau.
   DynVerdict AfterDynStep(real_t Vmax, real_t t, real_t tau)
   {
      const bool below = (Vmax < cfg_.V_off);
      if (!below) { t_above_ = t; }
      below_ = below ? t - t_above_ : 0.0;
      // `below` too: with T_dwell = 0 the counter alone (0 >= 0) would end
      // the phase at a step above V_off.
      if (below && below_ >= cfg_.T_dwell && tau >= cfg_.T_min) { return DynVerdict::EndOfPhase; }
      if (tau >= cfg_.T_max) { return DynVerdict::TooLong; }
      return DynVerdict::Stay;
   }

   /// Time max V has stayed below V_off in the current phase [s].
   real_t BelowTime() const { return below_; }
   /// Dynamic phases entered so far.
   int Phases() const { return phases_; }
   /// A restart from a QD checkpoint: the phases entered before it.
   void Restore(int phases)
   {
      MFEM_VERIFY(phases >= 0, "RegimeController::Restore: " << phases << " phases");
      phases_ = phases;
      below_ = 0.0;
   }
   /// True once the deck's max_dynamic_phases (> 0) have all been entered.
   bool LastPhase() const
   {
      return cfg_.max_dynamic_phases > 0 && phases_ >= cfg_.max_dynamic_phases;
   }
   const SwitchConfig &Config() const { return cfg_; }

private:
   SwitchConfig cfg_;
   real_t t_above_ = 0.0;   ///< last step with max V >= V_off, or the phase start
   real_t below_ = 0.0;
   int phases_ = 0;
};

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_REGIME_CONTROLLER_HPP

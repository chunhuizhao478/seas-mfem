// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// regime_controller.hpp — fully dynamic BP5 cycle (plan rev 3, Section
// 3.3 g): RegimeController, which decides when the driver switches regime.
// In: the global max slip rate and the dynamic phase's own clock tau.  Out:
// the next regime, or a stop.
//
//   Start(Vmax):      Vmax > V_on ? DYN : QD
//   GoDynamic(Vmax):  Vmax > V_on                     (after each QD step)
//   AfterDynStep(Vmax, tau, dt):
//     below = (Vmax < V_off) ? below + dt : 0
//     below >= T_dwell and tau >= T_min  ->  END_OF_PHASE (back to QD)
//     tau >= T_max                       ->  TOO_LONG (stop and report)
//     else                               ->  STAY
//
// "below" counts the time max V has stayed below V_off, so the phase ends
// T_dwell after the last step with max V >= V_off.  T_max stops the run
// instead of forcing a switch (plan R9).  Defaults in SwitchConfig.
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

   /// A dynamic phase starts: the dwell counter restarts.
   void EnterDynamic()
   {
      below_ = 0.0;
      phases_++;
   }

   /// After each dynamic step of length dt, at phase time tau.
   DynVerdict AfterDynStep(real_t Vmax, real_t tau, real_t dt)
   {
      below_ = (Vmax < cfg_.V_off) ? below_ + dt : 0.0;
      if (below_ >= cfg_.T_dwell && tau >= cfg_.T_min) { return DynVerdict::EndOfPhase; }
      if (tau >= cfg_.T_max) { return DynVerdict::TooLong; }
      return DynVerdict::Stay;
   }

   /// Time max V has stayed below V_off in the current phase [s].
   real_t BelowTime() const { return below_; }
   /// Dynamic phases entered so far.
   int Phases() const { return phases_; }
   /// True once the deck's max_dynamic_phases (> 0) have all been entered.
   bool LastPhase() const
   {
      return cfg_.max_dynamic_phases > 0 && phases_ >= cfg_.max_dynamic_phases;
   }
   const SwitchConfig &Config() const { return cfg_; }

private:
   SwitchConfig cfg_;
   real_t below_ = 0.0;
   int phases_ = 0;
};

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_REGIME_CONTROLLER_HPP

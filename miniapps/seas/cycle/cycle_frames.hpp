// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// cycle_frames.hpp — 2D in-plane design rev 3, phase A: the fully dynamic
// cycle's ONE fault orientation system.  The QD regime's operator builds the
// FaultFrameField with this configuration; the dynamic regime's Nitsche
// operator checks that its own fault setup asks for the same one and reads
// that field (fd_cycle_driver hands it over), so both regimes orient every
// fault face, frame and prescribed-slip face the same way.
//
// Scope (Chunhui Zhao, 2026-10-24).

#ifndef MFEM_SEAS_CYCLE_FRAMES_HPP
#define MFEM_SEAS_CYCLE_FRAMES_HPP

#include "mfem.hpp"
#include "cycle_config.hpp"
#include "../fault/fault_frame_field.hpp"

namespace mfem
{
namespace seas
{
namespace cycle
{

/// [cycle.dynamic] ref_normal, up, pin_min_margin (the deck keys the
/// dynamic regime read before phase A), with the reference pin and the
/// fault attribute.
inline FaultFrameConfig CycleFrameConfig(const CycleConfig &cfg)
{
   const DynamicConfig &dc = cfg.dyn;
   FaultFrameConfig f;
   FaultOrientationConfig &oc = f.orientation;
   for (int i = 0; i < 3; i++)
   {
      oc.ref_normal[i] = dc.ref_normal[i];
      oc.up[i] = dc.up[i];
      f.up[i] = dc.up[i];
   }
   oc.pin = FaultOrientationConfig::Pin::Reference;
   if (dc.pin_min_margin > 0.0) { oc.c_min = dc.pin_min_margin; }
   oc.fault_attrs = {cfg.fault_attr};
   return f;
}

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_FRAMES_HPP

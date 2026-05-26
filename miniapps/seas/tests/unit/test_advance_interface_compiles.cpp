// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_advance_interface_compiles.cpp — Phase 2 R-010 compile guard.
//
// Plan: PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24, §9.2.
//
// The retyped `AdvanceADERWithSubStep_Spatial(WaveOperator&, IFrictionIterator&,
// ...)` body invokes exactly these methods on the iterator: SetSubSteps,
// GetDeltaT, GetTimeWeights, Advance, SetDiagNumLocalFaultQPs (plus
// WaveOpLaw at the construction site).  This test mirrors that usage
// against the interface and a minimal stub IFrictionIterator whose only
// sub-step entry is Advance(...).  It must COMPILE — if the interface is
// missing any of those methods (e.g. the pre-fix call site used
// AdvanceWithSubStepStates, absent from the interface), the build fails
// here, at the interface, rather than deep in the driver.

#include "mfem.hpp"

#include "../../dynamic/friction_iterator.hpp"
#include "../../dynamic/fault_face_flux.hpp"   // DOFData

#include <functional>
#include <iostream>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace
{
// Minimal IFrictionIterator with no owned iterator — only the contract.
// Advance is the single sub-step entry point (mirrors the driver's :439
// call after the R-010 rename); the rest are trivial.
class StubFrictionIterator : public IFrictionIterator
{
public:
   void SetSubSteps(std::vector<real_t> deltaT,
                    std::vector<real_t> time_weights) override
   {
      deltaT_ = std::move(deltaT);
      weights_ = std::move(time_weights);
   }
   const std::vector<real_t> &GetDeltaT() const override { return deltaT_; }
   const std::vector<real_t> &GetTimeWeights() const override { return weights_; }

   void Advance(std::vector<DOFData> & /*dof_data*/,
                const std::vector<Vector> & /*fault_coords*/,
                const std::vector<std::vector<real_t>> & /*Qp*/,
                const std::vector<std::vector<real_t>> & /*Qm*/,
                real_t /*dt_macro*/, real_t /*t_macro_start*/,
                real_t * /*I_imp_plus_flat*/, real_t * /*I_imp_minus_flat*/,
                const std::function<void(real_t, real_t)> & /*nuc_callback*/)
      override
   { /* compile guard: no work */ }

   void SetDiagNumLocalFaultQPs(int /*n*/) override {}
   FaultFrictionLaw WaveOpLaw() const override { return FaultFrictionLaw::LSW; }

private:
   std::vector<real_t> deltaT_;
   std::vector<real_t> weights_;
};

// Mirror of the iterator calls in AdvanceADERWithSubStep_Spatial: if the
// interface drops any of these methods this function fails to compile.
void DriveLikeTheDriver(IFrictionIterator &iterator)
{
   iterator.SetDiagNumLocalFaultQPs(0);
   iterator.SetSubSteps({1.0}, {1.0});
   const std::vector<real_t> &dt = iterator.GetDeltaT();
   const std::vector<real_t> &w  = iterator.GetTimeWeights();
   (void) dt;
   (void) w;
   (void) iterator.WaveOpLaw();

   std::vector<DOFData> dof_data;
   std::vector<Vector> fault_coords;
   std::vector<std::vector<real_t>> Qp(1), Qm(1);
   auto noop = [](real_t, real_t) {};
   iterator.Advance(dof_data, fault_coords, Qp, Qm,
                    /*dt_macro=*/1.0, /*t_macro_start=*/0.0,
                    /*I_imp_plus_flat=*/nullptr, /*I_imp_minus_flat=*/nullptr,
                    noop);
}
}  // namespace

int main(int /*argc*/, char** /*argv*/)
{
   StubFrictionIterator stub;
   DriveLikeTheDriver(stub);
   std::cout << "test_advance_interface_compiles: IFrictionIterator covers "
                "every method AdvanceADERWithSubStep_Spatial invokes "
                "(SetSubSteps/GetDeltaT/GetTimeWeights/Advance/"
                "SetDiagNumLocalFaultQPs/WaveOpLaw). PASS\n";
   return 0;
}

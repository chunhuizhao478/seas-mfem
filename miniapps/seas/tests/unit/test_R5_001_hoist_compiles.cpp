// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 Phase 2a (§6.2.hoist) compile-time gate:
//   test_R5_001_dispatch_hoist_compiles
//
// The §6.2 dispatch restructure branches on `is_fault` at the outermost
// scope.  R5-001 warned that the pre-patch file declared `is_fault`,
// `dof_offset2`, `shape2`, and `I_nbr` INSIDE the `else (interior)`
// branch, so the restructured dispatch would fail to compile without
// the §6.2.hoist block.  This test is "build-only": the act of linking
// exercises the hoist block end-to-end in a TU that instantiates
// `WaveOperator<Mesh>` and calls AdvanceADER / UsePrecomputedFaceFluxes.
//
// If any of the four hoisted locals is missing from outer scope, the
// build fails with
//   "'is_fault' was not declared in this scope"
// and Phase 2a is blocked at the STOP gate.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../domain/boundary_config.hpp"

#include <iostream>

using namespace mfem;
using namespace mfem::seas;

namespace mfem { namespace seas { int g_seas_my_rank = 0; } }  // NOLINT

int main()
{
   std::cout << "\n=== TPV102 Phase 2a R5-001 hoist compile gate ===\n";

   // Minimal instantiation — we only need the binary to link.  The
   // restructured template is instantiated in wave_operator.o; this
   // test TU verifies the public API is callable with the new members.
   Mesh mesh = Mesh::MakeCartesian3D(2, 2, 2, Element::TETRAHEDRON,
                                     1000.0, 1000.0, 1000.0);
   BoundaryConfig bc;
   for (int i = 1; i <= 6; i++) { bc.absorbing_attrs.insert(i); }
   bc.fault_attr = 0;

   WaveOperator<Mesh> wave(mesh, /*order=*/1,
                           32.04e9, 32.04e9, 2670.0, bc);

   real_t Q_bg_zero[NUM_STATE] = {0};
   wave.SetAbsorbingBackground(Q_bg_zero);

   // Exercise the opt-in plumbing (compile + link only; no numerics).
   wave.UsePrecomputedFaceFluxes(false);
   wave.UsePrecomputedFaceFluxes(true);
   (void)wave.UsingPrecomputedFaceFluxes();
   (void)wave.GetPrecomputedFaceFluxes();
   (void)wave.GetFaultFaceSet();

   // One ADER step to exercise ComputeADERFaceFluxRHS's hoist block.
   Vector Q(wave.Height()); Q = 0.0;
   Vector Q1;
   wave.AdvanceADER(Q, /*dt=*/1e-5, /*order=*/2, Q1);

   std::cout << "  PASSED: Phase 2a dispatch hoist links and runs "
             "to completion\n";
   std::cout << "===================================================\n";
   std::cout << "Total:  1\n";
   std::cout << "Passed: 1\n";
   std::cout << "Failed: 0\n";
   std::cout << "===================================================\n";
   return 0;
}

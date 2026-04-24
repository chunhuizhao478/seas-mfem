// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// Phase 2D Step B sanity gate: verify BuildD4Mesh(add_fault=true)
// produces a usable fault at y = L/2.
//
// Checks:
//   1. Boundary triangles tagged attr=3 exist.
//   2. All attr=3 triangles have centroid y == L/2 (within tolerance).
//   3. WaveOperator constructed on this mesh populates fault_interior_faces_.
//   4. Fault face count matches the expected 8 (2x2 hex = 4 hex faces
//      split into 8 tet triangles).
//
// If any check fails, Phase 2D MUST stop and patch the fault tagging
// before proceeding.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/d4_tet_mesh.hpp"
#include "../../domain/boundary_config.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace mfem;
using namespace mfem::seas;

namespace mfem { namespace seas { int g_seas_my_rank = 0; } }  // NOLINT

int main()
{
   constexpr real_t kL      = 1000.0;
   constexpr real_t kRho    = 2670.0;
   constexpr real_t kLambda = 32.04e9;
   constexpr real_t kMu     = 32.04e9;

   std::cout << "=== Phase 2D Step B: BuildD4Mesh(true) fault sanity ===\n";

   Mesh mesh = BuildD4Mesh(/*add_fault=*/true, kL);

   int n_fault_bdr = 0, n_fault_misplaced = 0;
   for (int b = 0; b < mesh.GetNBE(); b++)
   {
      if (mesh.GetBdrAttribute(b) != 3) { continue; }
      n_fault_bdr++;
      Array<int> bv; mesh.GetBdrElementVertices(b, bv);
      real_t cy = 0.0;
      for (int v = 0; v < bv.Size(); v++) { cy += mesh.GetVertex(bv[v])[1]; }
      cy /= bv.Size();
      if (std::abs(cy - 0.5*kL) > 1e-6) { n_fault_misplaced++; }
   }
   std::cout << "  attr=3 boundary triangles      : " << n_fault_bdr
             << "   (expected 8)\n";
   std::cout << "  attr=3 triangles not at y=L/2  : "
             << n_fault_misplaced << "   (expected 0)\n";

   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr    = 3;
   WaveOperator<Mesh> wave(mesh, /*order=*/1, kLambda, kMu, kRho, bc);

   std::cout << "  GetFaultInteriorFaces().Size() : "
             << wave.GetFaultInteriorFaces().Size()
             << "   (expected 8)\n";
   std::cout << "  GetFaultSharedFaces().Size()   : "
             << wave.GetFaultSharedFaces().Size()
             << "   (expected 0 serial)\n";

   bool ok = (n_fault_bdr == 8)
          && (n_fault_misplaced == 0)
          && (wave.GetFaultInteriorFaces().Size() == 8);

   std::cout << "\n  overall: " << (ok ? "PASS" : "FAIL")
             << "   (Phase 2D " << (ok ? "may proceed" : "MUST STOP")
             << ")\n";

   return ok ? 0 : 1;
}

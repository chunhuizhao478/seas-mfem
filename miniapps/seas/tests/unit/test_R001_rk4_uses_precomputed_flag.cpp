// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.5.0 — REVIEW R-001 regression gate.
//
// Purpose: assert that `WaveOperator::Mult` (RK4 path) honors the
// `UsePrecomputedFaceFluxes(bool)` flag, symmetric with the ADER
// dispatch.  Pre-fix `ComputeFaceFluxRHS` ignored the flag, so the
// Arm 1 G_CONST_BFACE_LIFT probe was a vacuous no-op.
//
// Gate: on a Kuhn-split fault-less cartesian mesh, ||k_runtime
// - k_precomp||_inf > 0 for a nonzero input Q.  The Kuhn split's
// orbit-asymmetric vertex ordering makes the per-QP runtime
// BuildFrame(CalcOrtho) normals differ from the per-cell topology
// frame used by the precomputed-flux tables, so the two dispatches
// produce physically meaningful different rhs vectors.  Pre-fix the
// diff is 0 (flag ignored); post-fix the diff is strictly positive.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>

using namespace mfem;
using namespace mfem::seas;

namespace mfem { namespace seas { int g_seas_my_rank = 0; } }  // NOLINT

namespace
{

constexpr real_t kL      = 1000.0;
constexpr real_t kRho    = 2670.0;
constexpr real_t kLambda = 32.04e9;
constexpr real_t kMu     = 32.04e9;

real_t InfNormDiff(const Vector &a, const Vector &b)
{
   MFEM_VERIFY(a.Size() == b.Size(),
               "InfNormDiff: size mismatch");
   real_t m = 0.0;
   for (int i = 0; i < a.Size(); i++)
   {
      m = std::max(m, std::abs(a(i) - b(i)));
   }
   return m;
}

} // anonymous

int main()
{
   std::cout << "\n=== R-001 regression: Mult honors the precomputed "
             << "flux flag ===\n";

   // 2x2x2 Kuhn-split tetrahedral cartesian mesh, no fault.  All six
   // boundary attributes are absorbing; the interior faces are all
   // non-fault.  This exercises both the boundary and interior
   // non-fault dispatches added to ComputeFaceFluxRHS by R-001.
   const int n = 2;
   Mesh mesh = Mesh::MakeCartesian3D(n, n, n, Element::TETRAHEDRON,
                                     kL, kL, kL, /*sfc_ordering=*/false);

   BoundaryConfig bc;
   for (int i = 1; i <= 6; i++) { bc.absorbing_attrs.insert(i); }
   bc.fault_attr = 0;

   WaveOperator<Mesh> wave(mesh, /*order=*/1, kLambda, kMu, kRho, bc);

   real_t Q_bg_zero[NUM_STATE] = {0};
   wave.SetAbsorbingBackground(Q_bg_zero);

   const FiniteElementSpace &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();

   // Nonzero constant state in SXX (trivial to reason about; any flux
   // artefact manifests visibly).  A purely-zero Q would give k = 0
   // on both paths independent of the flag.
   Vector Q(NUM_STATE * ndof_total);
   Q = 0.0;
   for (int i = 0; i < ndof_total; i++)
   {
      Q(SXX * ndof_total + i) = 1.0e6;
   }

   // Runtime path (flag defaults to off).
   Vector k_runtime(NUM_STATE * ndof_total);
   wave.Mult(Q, k_runtime);

   // Precomputed path (flag on).
   wave.UsePrecomputedFaceFluxes(true);
   Vector k_precomp(NUM_STATE * ndof_total);
   wave.Mult(Q, k_precomp);

   const real_t path_diff = InfNormDiff(k_runtime, k_precomp);
   std::cout << "  || k_runtime - k_precomp ||_inf = "
             << std::scientific << std::setprecision(6) << path_diff
             << "\n";

   // Sanity: runtime rhs is nonzero — otherwise the test would pass
   // trivially for the wrong reason (if both paths produced all-zero
   // rhs, the diff would be zero and we would report "no diff" rather
   // than "flag ignored").
   real_t k_runtime_inf = 0.0;
   for (int i = 0; i < k_runtime.Size(); i++)
   {
      k_runtime_inf = std::max(k_runtime_inf, std::abs(k_runtime(i)));
   }
   std::cout << "  || k_runtime ||_inf             = "
             << std::scientific << std::setprecision(6) << k_runtime_inf
             << "\n";

   const bool runtime_nonzero = (k_runtime_inf > 0.0);
   const bool paths_differ    = (path_diff > 0.0);
   const bool verdict         = (runtime_nonzero && paths_differ);

   std::cout << "  runtime rhs nonzero : "
             << (runtime_nonzero ? "PASS" : "FAIL") << "\n";
   std::cout << "  runtime != precomp  : "
             << (paths_differ ? "PASS" : "FAIL")
             << "  (pre-fix would be 0 because flag was ignored in Mult)\n";
   std::cout << "  overall verdict     : "
             << (verdict ? "PASS" : "FAIL") << "\n";

   return verdict ? 0 : 1;
}

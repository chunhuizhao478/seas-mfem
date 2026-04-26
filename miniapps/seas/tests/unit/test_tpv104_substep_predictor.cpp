// Round-7 R-602: WaveOperator::ComputeADERSubStepStates round-trip test.
//
// Contract:
//   The CK Taylor expansion Q(τ) = Σ_{k=0}^{O-1} (τ^k / k!) · D(k) at any
//   set of nodes τ_o on [0, dt] must reproduce the macro-step time
//   integral I = ∫_0^{dt} Q(τ) dτ when integrated by an exact polynomial
//   quadrature on those same nodes.
//
//   For a polynomial of degree O − 1 in τ, an O-point Gauss-Legendre
//   quadrature on [0, dt] is exact.  Then:
//     Σ_o w_o · Q_per_node[o]  ==  I_macro / dt   (FP-bit equality up to
//                                                  the same recursion's
//                                                  rounding noise)
//   Equivalently: dt · Σ_o w_o · Q_per_node[o]  ==  I_macro.
//
//   We test this gate at ADER orders O ∈ {2, 3, 4} on a smooth plane-wave
//   IC (no fault, no boundary effects) — same wave-operator setup as
//   test_ader_linear_wave_equivalence.
//
// The test ALSO checks Q(0) ≈ Q (sub-step node τ=0 returns the input
// state to FP precision — the k=0 Taylor term is the only one with τ^0).

#include "mfem.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_LE(v, tol, msg) do { \
   num_tests++; \
   double vv = (v), tt = (tol); \
   if (vv <= tt) { num_passed++; \
      std::cout << "  PASSED: " << msg << "  (got " << std::scientific \
                << std::setprecision(3) << vv << ", tol " << tt << ")\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << "  (got " << std::scientific << std::setprecision(3) \
                << vv << ", expected <= " << tt << ")\n"; } \
} while (0)

namespace
{
constexpr real_t kRho    = 2670.0;
constexpr real_t kLambda = 32.04e9;
constexpr real_t kMu     = 32.04e9;

void FillPlaneWave(const FiniteElementSpace &fes, real_t k_wave,
                   real_t cp, real_t t, Vector &Q)
{
   const int ndof_total = fes.GetNDofs();
   Q.SetSize(NUM_STATE * ndof_total);
   Q = 0.0;
   const real_t lp = kLambda + 2.0 * kMu;
   real_t Q0[NUM_STATE] = {0};
   Q0[SXX] = -lp;
   Q0[SYY] = -kLambda;
   Q0[SZZ] = -kLambda;
   Q0[VX]  =  cp;
   for (int e = 0; e < fes.GetNE(); e++)
   {
      const FiniteElement *fe = fes.GetFE(e);
      ElementTransformation *Tr = fes.GetElementTransformation(e);
      const IntegrationRule &nodes = fe->GetNodes();
      Array<int> edofs; fes.GetElementDofs(e, edofs);
      for (int j = 0; j < fe->GetDof(); j++)
      {
         const IntegrationPoint &ip = nodes.IntPoint(j);
         Vector x(3);
         Tr->Transform(ip, x);
         const real_t phase = k_wave * (x(0) - cp * t);
         const real_t s = std::sin(phase);
         for (int c = 0; c < NUM_STATE; c++)
         {
            Q[c * ndof_total + edofs[j]] = Q0[c] * s;
         }
      }
   }
}

real_t L2Diff(const Vector &a, const Vector &b)
{
   real_t num = 0;
   for (int i = 0; i < a.Size(); i++)
   {
      real_t d = a(i) - b(i);
      num += d * d;
   }
   return std::sqrt(num);
}

real_t L2Norm(const Vector &a)
{
   real_t num = 0;
   for (int i = 0; i < a.Size(); i++) { num += a(i) * a(i); }
   return std::sqrt(num);
}

// O-point Gauss-Legendre quadrature on [0, dt].  Returns nodes (positions)
// and weights such that ∫_0^dt f(τ) dτ ≈ Σ w_o · f(node_o) for f a
// polynomial of degree ≤ 2O − 1.  We need f of degree O − 1, which is
// strictly within the exactness window.
void GaussLegendreOnInterval(real_t dt, int O,
                             std::vector<real_t> &nodes,
                             std::vector<real_t> &weights)
{
   // MFEM's standard segment quadrature already provides Gauss-Legendre.
   // Order argument 2*O - 1 picks O-point GL.
   const IntegrationRule &ir = IntRules.Get(Geometry::SEGMENT, 2 * O - 1);
   const int n = ir.GetNPoints();
   nodes.resize(n);
   weights.resize(n);
   for (int o = 0; o < n; o++)
   {
      const IntegrationPoint &ip = ir.IntPoint(o);
      // ip.x ∈ [0, 1] in MFEM segment convention; scale to [0, dt].
      nodes[o]   = dt * ip.x;
      weights[o] = dt * ip.weight;
   }
}

bool RunOrderTest(WaveOperator<Mesh> &wave, const FiniteElementSpace &fes,
                  real_t k_wave, real_t cp, real_t dt_test, int order)
{
   std::cout << "\n-- ADER-O" << order << "  dt = " << std::scientific
             << std::setprecision(3) << dt_test << " --\n";

   Vector Q; FillPlaneWave(fes, k_wave, cp, /*t=*/0.0, Q);

   // (1) Macro-step time integral via existing ComputeADERTimeIntegrated.
   Vector I_macro;
   wave.ComputeADERTimeIntegrated(Q, dt_test, order, I_macro);

   // (2) Per-sub-step states at Gauss-Legendre nodes.
   std::vector<real_t> tau_nodes, weights;
   GaussLegendreOnInterval(dt_test, order, tau_nodes, weights);
   std::vector<Vector> Q_per_node;
   wave.ComputeADERSubStepStates(Q, dt_test, order, tau_nodes, Q_per_node);

   if (static_cast<int>(Q_per_node.size()) != order)
   {
      std::cout << "  FAILED: Q_per_node.size() = " << Q_per_node.size()
                << ", expected " << order << "\n";
      num_tests++; num_failed++;
      return false;
   }

   // (3) Quadrature: Σ_o w_o · Q_per_node[o] should equal I_macro
   //     (since the Taylor poly Q(τ) is degree O−1 in τ, and the O-point
   //     GL rule integrates up to 2O−1, exactly).
   Vector I_quad(I_macro.Size());
   I_quad = 0.0;
   for (int o = 0; o < order; o++)
   {
      I_quad.Add(weights[o], Q_per_node[o]);
   }

   const real_t err  = L2Diff(I_quad, I_macro);
   const real_t norm = L2Norm(I_macro);
   const real_t rel  = (norm > 0.0) ? err / norm : err;
   std::cout << "  |I_quad − I_macro|_2 = " << err
             << "  (rel = " << rel << ")\n";

   // FP precision floor: predictor recursion has ~order·ULP rounding per
   // term; compounded with the per-DOF wave amplitude (~kLambda*1.0 ~3e10),
   // the absolute floor is ~order·1e-16·1e10 = 1e-6 in absolute terms or
   // ~1e-12 relative for our smooth IC.  Use 1e-10 relative tolerance —
   // robust to the MFEM kernel's rounding without being so loose as to
   // accept a real bug.
   TEST_LE(rel, 1.0e-10,
           ("Quadrature roundtrip rel err at O=" +
            std::to_string(order) + " ≤ 1e-10").c_str());

   // (4) Q(0) ≈ Q.  Insert τ = 0 as the first node; check identity.
   {
      std::vector<real_t> tau0_nodes(1, 0.0);
      std::vector<Vector> Q_at_zero;
      wave.ComputeADERSubStepStates(Q, dt_test, order, tau0_nodes, Q_at_zero);
      const real_t err0 = L2Diff(Q_at_zero[0], Q);
      std::cout << "  |Q(τ=0) − Q|_2 = " << err0 << "\n";
      TEST_LE(err0, 0.0,
              ("Q(τ=0) bit-equal to Q at O=" +
               std::to_string(order)).c_str());
   }
   return true;
}

} // anonymous

int main()
{
   std::cout << "\n=== Round-7 R-602: ComputeADERSubStepStates round-trip ===\n";

   const real_t L = 1000.0;
   const int n = 4;
   const int p = 4;
   Mesh mesh = Mesh::MakeCartesian3D(n, n, n, Element::HEXAHEDRON, L, L, L);

   BoundaryConfig bc;
   for (int i = 1; i <= 6; i++) { bc.absorbing_attrs.insert(i); }
   bc.fault_attr = 0;

   WaveOperator<Mesh> wave(mesh, p, kLambda, kMu, kRho, bc);
   {
      real_t Q_bg_zero[NUM_STATE] = {0};
      wave.SetAbsorbingBackground(Q_bg_zero);
   }
   const FiniteElementSpace &fes = wave.GetFESpace();

   const real_t lp = kLambda + 2.0 * kMu;
   const real_t cp = std::sqrt(lp / kRho);
   const real_t lambda_wave = 8.0 * L;
   const real_t k_wave = 2.0 * M_PI / lambda_wave;
   const real_t h = L / n;
   const real_t dt_cfl = h / cp;
   const real_t dt_test = 0.01 * dt_cfl;

   for (int order = 2; order <= 4; order++)
   {
      RunOrderTest(wave, fes, k_wave, cp, dt_test, order);
   }

   std::cout << "\n=== Summary: " << num_passed << "/" << num_tests
             << " passed, " << num_failed << " failed ===\n";
   return (num_failed == 0) ? 0 : 1;
}

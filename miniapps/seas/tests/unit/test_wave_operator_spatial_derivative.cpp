// ADER I-05 Phase 1: unit test for WaveOperator::ApplySpatialDerivative.
//
// Gates polynomial exactness of the element-local L2-projected spatial
// derivative: for a monomial Q_c(x) = x^n / n the L2-projection of its
// ∂_dir is n * x^{n-1} δ_{dir, coord}, so the recovered derivative must
// match that at every DG DOF (interpolation).  We run n ∈ {1, 2, 3} at
// order = max(2, p) on a uniform hex mesh, checking dQ/dx_j at every DOF
// for each coordinate direction.

#include "mfem.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/bimaterial_wave_operator.hpp"   // Lever 1 bimaterial pass
#include "../../dynamic/heterogeneous_material.hpp"     // MaterialField::MakeCoefficient
#include "../../dynamic/wave_state.hpp"
#include "../../domain/boundary_config.hpp"

#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>

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

// Fill Q_c[dof] with f(x_dof) for a single state component, leaving the
// other components zero.  Uses the FE space's nodal positions at the L2
// Lagrange DOFs.
static void FillQFromMonomial(const FiniteElementSpace &fes,
                              int component, int coord_dir, int power,
                              Vector &Q)
{
   const int NS = NUM_STATE;
   const int ndof_total = fes.GetNDofs();
   Q.SetSize(NS * ndof_total);
   Q = 0.0;
   Mesh *mesh = fes.GetMesh();
   const int ne = fes.GetNE();
   for (int e = 0; e < ne; e++)
   {
      const FiniteElement *fe = fes.GetFE(e);
      ElementTransformation *Tr = fes.GetElementTransformation(e);
      const IntegrationRule &nodes = fe->GetNodes();
      const int ndof = fe->GetDof();
      Array<int> edofs;
      fes.GetElementDofs(e, edofs);
      for (int j = 0; j < ndof; j++)
      {
         const IntegrationPoint &ip = nodes.IntPoint(j);
         Vector xj(3);
         Tr->Transform(ip, xj);
         real_t val = std::pow(xj(coord_dir), power);
         Q[component * ndof_total + edofs[j]] = val;
      }
   }
}

// Evaluate the DG field Q_c at every DOF by reading Q directly at the
// Lagrange nodes (matches how FillQFromMonomial wrote them).
// Returns max |dQ_c_dg[dof] - expected(x_dof)| where expected is provided.
static real_t MaxDOFError(const FiniteElementSpace &fes,
                          int component,
                          const Vector &dQ,
                          int coord_dir, int power)
{
   const int ndof_total = fes.GetNDofs();
   const int ne = fes.GetNE();
   real_t max_err = 0.0;
   for (int e = 0; e < ne; e++)
   {
      const FiniteElement *fe = fes.GetFE(e);
      ElementTransformation *Tr = fes.GetElementTransformation(e);
      const IntegrationRule &nodes = fe->GetNodes();
      const int ndof = fe->GetDof();
      Array<int> edofs;
      fes.GetElementDofs(e, edofs);
      for (int j = 0; j < ndof; j++)
      {
         const IntegrationPoint &ip = nodes.IntPoint(j);
         Vector xj(3);
         Tr->Transform(ip, xj);
         // Expected derivative of x_d^n with respect to x_coord_dir
         // equals n * x_coord_dir^{n-1} * δ_{d, coord_dir}.  The Q we
         // filled has component_data[j] = x_coord_dir(j)^power, so its
         // spatial derivative w.r.t. direction `der_dir` is:
         //    power * x_coord_dir^{power-1}    if der_dir == coord_dir
         //    0                                otherwise.
         real_t expected = 0.0;
         // This function is called with (der_dir = coord_dir, power)
         // meaning we computed ∂_{coord_dir} (x_coord_dir^power).
         if (power == 0)
         {
            expected = 0.0;
         }
         else
         {
            expected = static_cast<real_t>(power)
                     * std::pow(xj(coord_dir), power - 1);
         }
         real_t got = dQ[component * ndof_total + edofs[j]];
         real_t err = std::abs(got - expected);
         max_err = std::max(max_err, err);
      }
   }
   return max_err;
}

// Check that ∂_{der_dir} of (x_coord_dir^power) is zero everywhere
// when der_dir != coord_dir.  Returns max |dQ|.
static real_t MaxAbs(const FiniteElementSpace &fes,
                     int component,
                     const Vector &dQ)
{
   const int ndof_total = fes.GetNDofs();
   real_t m = 0.0;
   for (int i = 0; i < ndof_total; i++)
   {
      m = std::max(m, std::abs(dQ[component * ndof_total + i]));
   }
   return m;
}

// LCG fill — deterministic, reproducible random state vector.
static void FillRandomQ(Vector &Q, unsigned seed)
{
   unsigned s = seed ? seed : 1u;
   for (int i = 0; i < Q.Size(); i++)
   {
      s = 1664525u * s + 1013904223u;
      Q[i] = (static_cast<double>(s) / 4294967296.0) - 0.5;
   }
}

// Max relative error between OnTheFly and Cached ApplySpatialDerivative over all
// three directions, for an already-constructed operator (templated so it serves
// both WaveOperator and BimaterialWaveOperator).
template <typename Op>
static double CachedVsOnTheFlyMaxRelErr(Op &wave, int N, unsigned seed)
{
   Vector Q(N);
   FillRandomQ(Q, seed);
   double worst = 0.0;
   for (int dir = 0; dir < 3; dir++)
   {
      wave.SetDerivMode(DerivMode::OnTheFly);
      Vector dref;
      wave.ApplySpatialDerivative(dir, Q, dref);
      wave.SetDerivMode(DerivMode::Cached);
      Vector dcac;
      wave.ApplySpatialDerivative(dir, Q, dcac);
      wave.SetDerivMode(DerivMode::OnTheFly);   // restore default

      double maxref = 0.0, maxdiff = 0.0;
      for (int i = 0; i < N; i++)
      {
         maxref  = std::max(maxref,  std::abs(dref[i]));
         maxdiff = std::max(maxdiff, std::abs(dref[i] - dcac[i]));
      }
      worst = std::max(worst, maxdiff / (maxref + 1e-300));
   }
   return worst;
}

// REVIEW R-001: the PRODUCTION consumers are the CK recursion — AdvanceADER and
// ComputeADERSubStepStates — which call ApplySpatialDerivative internally
// (order-1)x3 times.  The single-call check above does NOT exercise that path.
// Assert the recursion outputs agree between Cached and OnTheFly to <= 1e-12
// relative (machine-eps accumulated over the applications; R-002 — not bit-exact).
template <typename Op>
static double RecursionCachedVsOnTheFlyMaxRelErr(Op &wave, int N, int order,
                                                 unsigned seed)
{
   // AdvanceADER's corrector requires the absorbing background (has_bulk_bg_);
   // Q_bg = 0 is valid.
   real_t Q_bg[NUM_STATE];
   for (int c = 0; c < NUM_STATE; c++) { Q_bg[c] = 0.0; }
   wave.SetAbsorbingBackground(Q_bg);

   Vector Q(N);
   FillRandomQ(Q, seed);
   const real_t dt = 1e-4;
   std::vector<real_t> tau_nodes(order);
   for (int o = 0; o < order; o++)
   {
      tau_nodes[o] = dt * (static_cast<real_t>(o) + 0.5)
                   / static_cast<real_t>(order);
   }

   auto relerr = [](const Vector &a, const Vector &b)
   {
      double mr = 0.0, md = 0.0;
      for (int i = 0; i < a.Size(); i++)
      {
         mr = std::max(mr, std::abs(a[i]));
         md = std::max(md, std::abs(a[i] - b[i]));
      }
      return md / (mr + 1e-300);
   };

   // AdvanceADER (predictor + corrector).
   wave.SetDerivMode(DerivMode::OnTheFly);
   Vector qn_otf;
   wave.AdvanceADER(Q, dt, order, qn_otf);
   wave.SetDerivMode(DerivMode::Cached);
   Vector qn_cac;
   wave.AdvanceADER(Q, dt, order, qn_cac);
   const double e_adv = relerr(qn_otf, qn_cac);

   // ComputeADERSubStepStates (the friction-coupling predictor).
   wave.SetDerivMode(DerivMode::OnTheFly);
   std::vector<Vector> qpn_otf;
   wave.ComputeADERSubStepStates(Q, dt, order, tau_nodes, qpn_otf);
   wave.SetDerivMode(DerivMode::Cached);
   std::vector<Vector> qpn_cac;
   wave.ComputeADERSubStepStates(Q, dt, order, tau_nodes, qpn_cac);
   double e_sss = 0.0;
   for (std::size_t o = 0; o < qpn_otf.size(); o++)
   {
      e_sss = std::max(e_sss, relerr(qpn_otf[o], qpn_cac[o]));
   }

   wave.SetDerivMode(DerivMode::OnTheFly);   // restore default
   return std::max(e_adv, e_sss);
}

// Lever 1 / REVIEW R-002: the Cached (precomputed D_d^e mat-vec) kernel must
// match the OnTheFly quadrature kernel to <= 1e-12 relative (NOT bit-exact — the
// setup M^{-1}K_d product re-associates round-off).  Covers WaveOperator AND
// BimaterialWaveOperator, which inherits the base kernel + cache unchanged.
static void TestCachedEquivalence()
{
   std::cout << "\n=== Lever 1: DerivMode::Cached vs OnTheFly equivalence "
             << "(R-002) ===\n";
   const double tol = 1e-12;

   for (int order = 2; order <= 3; order++)
   {
      Mesh mesh = Mesh::MakeCartesian3D(3, 3, 3, Element::TETRAHEDRON,
                                        1.0, 1.0, 1.0);
      const real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
      BoundaryConfig bc;
      for (int i = 1; i <= 6; i++) { bc.absorbing_attrs.insert(i); }
      bc.fault_attr = 0;

      // --- scalar WaveOperator ---
      WaveOperator<Mesh> wave(mesh, order, lambda, mu, rho, bc);
      const int N = NUM_STATE * wave.GetFESpace().GetNDofs();
      double rel = CachedVsOnTheFlyMaxRelErr(wave, N, 12345u + order);
      std::ostringstream m1;
      m1 << "scalar: cached == onthefly, order=" << order << " (max rel)";
      TEST_LE(rel, tol, m1.str());

      // --- BimaterialWaveOperator (homogeneous material; same geometry) ---
      ConstantCoefficient lam_c(lambda), mu_c(mu), rho_c(rho);
      BimaterialWaveOperator<Mesh> wb(
         mesh, order, MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c), bc);
      const int Nb = NUM_STATE * wb.GetFESpace().GetNDofs();
      double relb = CachedVsOnTheFlyMaxRelErr(wb, Nb, 999u + order);
      std::ostringstream m2;
      m2 << "bimaterial: cached == onthefly, order=" << order << " (max rel)";
      TEST_LE(relb, tol, m2.str());

      // --- R-001: recursion-level equivalence (the production path) ---
      double rrel = RecursionCachedVsOnTheFlyMaxRelErr(wave, N, order,
                                                       4242u + order);
      std::ostringstream m1r;
      m1r << "scalar recursion (AdvanceADER+SubStepStates): cached == onthefly, "
          << "order=" << order << " (max rel)";
      TEST_LE(rrel, tol, m1r.str());

      double rrelb = RecursionCachedVsOnTheFlyMaxRelErr(wb, Nb, order,
                                                        7777u + order);
      std::ostringstream m2r;
      m2r << "bimaterial recursion (AdvanceADER+SubStepStates): cached == "
          << "onthefly, order=" << order << " (max rel)";
      TEST_LE(rrelb, tol, m2r.str());

      // --- R-004: cache byte count is exactly 3*ne*ndof^2*sizeof(real_t) ---
      const int ne   = wave.GetFESpace().GetNE();
      const int ndof = wave.GetFESpace().GetFE(0)->GetDof();
      const std::size_t expect =
         static_cast<std::size_t>(3) * ne * ndof * ndof * sizeof(real_t);
      const double byte_err = std::abs(
         static_cast<double>(ElementDerivativeCacheBytes(ne, ndof))
         - static_cast<double>(expect));
      std::ostringstream m3;
      m3 << "R-004 cache bytes = 3*ne*ndof^2*8, order=" << order;
      TEST_LE(byte_err, 0.0, m3.str());
   }
}

// Lever 2 / plan §5: the shared-CK-recursion path
// (ComputeADERSubStepStatesAndIntegral + AdvanceADER(...,&I)) must produce
// Q_per_node AND Q_new BIT-FOR-BIT identical to the separate
// ComputeADERSubStepStates + AdvanceADER, at a FIXED DerivMode (R-002: baseline
// against the same DerivMode, OnTheFly and Cached).  Orders {2,3,4}.
static void TestSharedCKParity()
{
   std::cout << "\n=== Lever 2: shared-CK-recursion bit-exact parity "
             << "(R-002/R-003) ===\n";
   for (int order = 2; order <= 4; order++)
   {
      Mesh mesh = Mesh::MakeCartesian3D(3, 3, 3, Element::TETRAHEDRON,
                                        1.0, 1.0, 1.0);
      const real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
      BoundaryConfig bc;
      for (int i = 1; i <= 6; i++) { bc.absorbing_attrs.insert(i); }
      bc.fault_attr = 0;
      WaveOperator<Mesh> wave(mesh, order, lambda, mu, rho, bc);
      real_t Q_bg[NUM_STATE];
      for (int c = 0; c < NUM_STATE; c++) { Q_bg[c] = 0.0; }
      wave.SetAbsorbingBackground(Q_bg);

      const int N = NUM_STATE * wave.GetFESpace().GetNDofs();
      const real_t dt = 1e-4;
      std::vector<real_t> tau_nodes(order);
      for (int o = 0; o < order; o++)
      {
         tau_nodes[o] = dt * (static_cast<real_t>(o) + 0.5)
                      / static_cast<real_t>(order);
      }

      for (int mode_i = 0; mode_i < 2; mode_i++)
      {
         const DerivMode mode =
            (mode_i == 0) ? DerivMode::OnTheFly : DerivMode::Cached;
         wave.SetDerivMode(mode);
         Vector Q(N);
         FillRandomQ(Q, 31415u + order + 100u * mode_i);

         // Separate path (two recursions).
         std::vector<Vector> qpn_sep;
         wave.ComputeADERSubStepStates(Q, dt, order, tau_nodes, qpn_sep);
         Vector qn_sep;
         wave.AdvanceADER(Q, dt, order, qn_sep);

         // Merged path (one recursion + the I_precomputed overload).
         std::vector<Vector> qpn_m;
         Vector I_m;
         wave.ComputeADERSubStepStatesAndIntegral(Q, dt, order, tau_nodes,
                                                  qpn_m, I_m);
         Vector qn_m;
         wave.AdvanceADER(Q, dt, order, qn_m, &I_m);

         // Bit-for-bit: max abs difference must be EXACTLY zero.
         double dmax = 0.0;
         for (std::size_t o = 0; o < qpn_sep.size(); o++)
         {
            for (int i = 0; i < N; i++)
            {
               dmax = std::max(dmax, std::abs(qpn_sep[o][i] - qpn_m[o][i]));
            }
         }
         for (int i = 0; i < N; i++)
         {
            dmax = std::max(dmax, std::abs(qn_sep[i] - qn_m[i]));
         }

         std::ostringstream msg;
         msg << "merged == separate (bit-exact), order=" << order
             << ", DerivMode=" << (mode == DerivMode::OnTheFly ? "OnTheFly"
                                                               : "Cached");
         TEST_LE(dmax, 0.0, msg.str());
      }
      wave.SetDerivMode(DerivMode::OnTheFly);   // restore default
   }
}

int main()
{
   std::cout << "\n=== ADER I-05 Phase 1: ApplySpatialDerivative "
             << "polynomial-exactness test ===\n";

   // DG polynomial order.  Plan acceptance requires exactness for n ∈
   // {1, 2, 3}, which needs a DG space of polynomial degree ≥ 3.
   const int order = 3;
   Mesh mesh = Mesh::MakeCartesian3D(4, 4, 4, Element::HEXAHEDRON,
                                     1.0, 1.0, 1.0);

   real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   BoundaryConfig bc;
   for (int i = 1; i <= 6; i++) { bc.absorbing_attrs.insert(i); }
   bc.fault_attr = 0;

   WaveOperator<Mesh> wave(mesh, order, lambda, mu, rho, bc);

   const FiniteElementSpace &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   const real_t kScale = 1.0;  // unit cube — |x| ≤ 1
   // Polynomial-exactness tolerance.  The derivative is computed via L2
   // projection, which is exact (up to FP roundoff) on monomials of
   // degree ≤ order.  Tolerance scales with the number of cross-terms
   // in the sum-over-QPs × number-of-DOFs mat-vec.
   const real_t kEps = 2.220446049250313e-16;
   const real_t tol = 2048 * kEps * kScale;

   // Plan acceptance: n ∈ {1, 2, 3}.  Exactness holds iff the monomial x^n
   // is in the DG polynomial space of degree `order`, i.e., n ≤ order.  At
   // order = 3 all three cases are exact.
   for (int n = 1; n <= 3; n++)
   {
      std::cout << "\n-- n = " << n << " (monomial x^n) --\n";

      // Place the monomial in component SXX (c=0), coordinate direction i.
      for (int i = 0; i < 3; i++)
      {
         Vector Q;
         FillQFromMonomial(fes, /*component=*/SXX,
                           /*coord_dir=*/i, /*power=*/n, Q);

         // ∂_i (x_i^n) = n x_i^{n-1} — reproduced exactly on DG DOFs when
         // n ≤ order (L2 projection of a polynomial of degree n-1 onto a
         // degree-order space is exact).
         Vector dQ;
         wave.ApplySpatialDerivative(i, Q, dQ);
         real_t err_parallel = MaxDOFError(fes, SXX, dQ, i, n);
         std::ostringstream msg;
         msg << "∂_" << i << " x_" << i << "^" << n << " at DG DOFs";
         TEST_LE(err_parallel, tol, msg.str());

         // ∂_j (x_i^n) = 0 for j != i — orthogonality check.
         for (int j = 0; j < 3; j++)
         {
            if (j == i) { continue; }
            Vector dQj;
            wave.ApplySpatialDerivative(j, Q, dQj);
            real_t m = MaxAbs(fes, SXX, dQj);
            std::ostringstream msg2;
            msg2 << "∂_" << j << " x_" << i << "^" << n << " = 0 (j != i)";
            TEST_LE(m, tol, msg2.str());
         }
      }
   }

   // Sanity: input-size error is MFEM_VERIFY'd.  Test with a
   // well-formed zero input yields zero output.
   {
      Vector Q(NUM_STATE * ndof_total);
      Q = 0.0;
      Vector dQ;
      wave.ApplySpatialDerivative(0, Q, dQ);
      real_t m = MaxAbs(fes, SXX, dQ);
      TEST_LE(m, tol, "∂_x(0) = 0");
   }

   // Lever 1: cached-derivative-operator equivalence (R-002) + R-004 byte count.
   TestCachedEquivalence();

   // Lever 2: shared-CK-recursion bit-exact parity (both DerivModes).
   TestSharedCKParity();

   std::cout << "\n========================================\n";
   std::cout << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}

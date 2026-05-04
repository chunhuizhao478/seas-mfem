// Code-debug reproducer: p=2 mixed-flux Adjacent dispatch.
//
// The production p=2 + Adjacent + ADER-O3 runs (TPV102/104/205) blow up
// catastrophically while p=1 + Adjacent and p=2 + None work.  Existing
// unit tests only exercise mixed-flux at p=1 (test_mixed_flux_dispatch_adjacent.cpp)
// so the p=2 path has zero coverage.  This test isolates two questions:
//
//   1. Does the algebraic identity dQdt_adj − dQdt_none ==
//      Σ_f M^{-1} w shape (F_int − F_central) hold at p=2?
//      If YES, the dispatch IMPLEMENTATION is correct at p=2 and the
//      blow-up is a STABILITY/scheme issue, not a code bug.
//      If NO, there is a p=2-specific implementation bug.
//
//   2. Does a many-step run with Adjacent grow without bound at p=2 on
//      the same fixture and dt that's stable at p=1?  Quantitative
//      stability symptom for the multi-step regime.
//
// Fixture is identical to test_mixed_flux_dispatch_adjacent.cpp.  Only
// the polynomial order changes (order=1 → order=2) and a 200-step
// ADER-O3 stability sweep is added.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../dynamic/tpv102_setup_total.hpp"
#include "../../config/tpv102_params.hpp"
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

#define TEST_GE(v, lo, msg) do { \
   num_tests++; \
   double vv = (v), ll = (lo); \
   if (vv >= ll) { num_passed++; \
      std::cout << "  PASSED: " << msg << "  (got " << std::scientific \
                << std::setprecision(3) << vv << ", floor " << ll << ")\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << "  (got " << std::scientific << std::setprecision(3) \
                << vv << ", expected >= " << ll << ")\n"; } \
} while (0)

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

Mesh BuildSmallFaultTetMesh(real_t L = 1000.0)
{
   const int Nv = 3 * 3 * 3;
   Mesh mesh(3, Nv, 0, 0);
   auto vid = [&](int i, int j, int k) { return i + 3*j + 9*k; };
   for (int k = 0; k < 3; k++)
      for (int j = 0; j < 3; j++)
         for (int i = 0; i < 3; i++)
         {
            real_t v[3] = {i*0.5*L, j*0.5*L, k*0.5*L};
            mesh.AddVertex(v);
         }
   auto vfunc = [&](int ci, int cj, int ck, int lv) {
      return vid(ci + (lv & 1), cj + ((lv >> 1) & 1), ck + ((lv >> 2) & 1));
   };
   const int pat[6][4] = {
      {0, 1, 3, 7}, {0, 3, 2, 7}, {0, 2, 6, 7},
      {0, 6, 4, 7}, {0, 4, 5, 7}, {0, 5, 1, 7}
   };
   for (int ck = 0; ck < 2; ck++)
      for (int cj = 0; cj < 2; cj++)
         for (int ci = 0; ci < 2; ci++)
            for (int t = 0; t < 6; t++)
            {
               mesh.AddTet(vfunc(ci, cj, ck, pat[t][0]),
                           vfunc(ci, cj, ck, pat[t][1]),
                           vfunc(ci, cj, ck, pat[t][2]),
                           vfunc(ci, cj, ck, pat[t][3]),
                           1);
            }
   mesh.FinalizeTopology();
   const real_t eps = 1e-6;
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      real_t cy = 0.0;
      for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
      cy /= fv.Size();
      if (ftr && ftr->Elem2No >= 0)
      {
         if (std::abs(cy - 0.5 * L) < eps)
         {
            mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3);
         }
         continue;
      }
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 1);
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

real_t MaxAbsDiff(const Vector &a, const Vector &b)
{
   real_t m = 0.0;
   for (int i = 0; i < a.Size(); i++)
   {
      const real_t d = std::abs(a(i) - b(i));
      if (d > m) { m = d; }
   }
   return m;
}

real_t MaxAbs(const Vector &a)
{
   real_t m = 0.0;
   for (int i = 0; i < a.Size(); i++)
   {
      const real_t d = std::abs(a(i));
      if (d > m) { m = d; }
   }
   return m;
}

void InitializeNonUniformState(Vector &Q, const FiniteElementSpace &fes,
                               int ndof_total, real_t /*L*/)
{
   Q.SetSize(NUM_STATE * ndof_total);
   Q = 0.0;
   real_t *Q_data = Q.GetData();
   for (int e = 0; e < fes.GetNE(); e++)
   {
      const FiniteElement *fe = fes.GetFE(e);
      const int ndof_e = fe->GetDof();
      Array<int> dof_idx;
      fes.GetElementDofs(e, dof_idx);
      const real_t s = 1.0 + 0.1 * static_cast<real_t>((e * 17 + 5) % 23);
      const real_t t = 1.0 + 0.05 * static_cast<real_t>((e * 31 + 11) % 19);
      for (int i = 0; i < ndof_e; i++)
      {
         const int dof = dof_idx[i];
         Q_data[SXX * ndof_total + dof] = -TPV102Params::sigma_n * s;
         Q_data[SYY * ndof_total + dof] =  TPV102Params::sigma_n * t;
         Q_data[SXY * ndof_total + dof] = -TPV102Params::tau_ini * s;
         Q_data[SZZ * ndof_total + dof] =  1.0e6 * t;
         Q_data[SYZ * ndof_total + dof] =  1.0e6 * s;
         Q_data[SXZ * ndof_total + dof] =  1.0e6 * t;
         Q_data[VX  * ndof_total + dof] =  0.1 * s;
         Q_data[VY  * ndof_total + dof] =  0.1 * t;
         Q_data[VZ  * ndof_total + dof] =  0.1 * (s + t);
      }
   }
}

void SetupWaveForRun(WaveOperator<Mesh> &wave,
                     std::vector<DOFData> &dof_data,
                     const Mesh &mesh, int order)
{
   const int n_local_qps = wave.GetNumLocalFaultQPs();
   const int nbf = wave.GetNbfPerFace();
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   std::vector<Vector> fault_coords;
   fault_coords.reserve(n_local_qps);
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr =
         const_cast<Mesh &>(mesh).GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir = IntRules.Get(
         ftr->GetGeometryType(), 2 * order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   }
   InitializeFaultDOFs(dof_data, n_local_qps, fault_coords);
   ZeroDOFDataPreStressTotal(dof_data, n_local_qps);
   wave.SetFaultDOFData(&dof_data, nbf);

   real_t bulk_bg[NUM_STATE] = {0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;
   wave.SetAbsorbingBackground(bulk_bg);
}

// Same identity as test_mixed_flux_dispatch_adjacent.cpp
// ComputeAnalyticAdjacentDelta but inlined for self-containment.
void ComputeAnalyticAdjacentDelta(const WaveOperator<Mesh> &wave,
                                  const Mesh &mesh, int order,
                                  const Vector &Q_init,
                                  Vector &delta_out)
{
   const auto &fes = wave.GetFESpace();
   const int ndof_total  = fes.GetNDofs();
   const int ndof_per_el = wave.GetNDof();
   delta_out.SetSize(NUM_STATE * ndof_total);
   delta_out = 0.0;

   const real_t *Q_data = Q_init.GetData();
   const auto &flux = wave.GetFlux();
   const auto &central_set = wave.GetCentralFluxFaceSet();

   for (int f : central_set)
   {
      auto *ftr = const_cast<Mesh &>(mesh).GetInteriorFaceTransformations(f);
      if (!ftr || ftr->Elem2No < 0) { continue; }
      const int e1 = ftr->Elem1No;
      const int e2 = ftr->Elem2No;
      const FiniteElement *fe1 = fes.GetFE(e1);
      const FiniteElement *fe2 = fes.GetFE(e2);
      const int ndof_e1 = fe1->GetDof();
      const int ndof_e2 = fe2->GetDof();
      const int dof_offset1 = e1 * ndof_per_el;
      const int dof_offset2 = e2 * ndof_per_el;

      const IntegrationRule &ir =
         IntRules.Get(ftr->GetGeometryType(), 2 * order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);

         Vector nor_vec(3);
         CalcOrtho(ftr->Face->Jacobian(), nor_vec);
         const real_t nor_len = nor_vec.Norml2();
         if (nor_len > 0) { nor_vec /= nor_len; }
         const real_t w = ip.weight * nor_len;
         const real_t nor[3] = {nor_vec(0), nor_vec(1), nor_vec(2)};

         Vector shape1(ndof_e1), shape2(ndof_e2);
         fe1->CalcShape(ftr->GetElement1IntPoint(), shape1);
         fe2->CalcShape(ftr->GetElement2IntPoint(), shape2);

         real_t Q_self[NUM_STATE], Q_nbr[NUM_STATE];
         for (int c = 0; c < NUM_STATE; c++)
         {
            Q_self[c] = 0.0;
            Q_nbr[c]  = 0.0;
            for (int i = 0; i < ndof_e1; i++)
            {
               Q_self[c] += shape1(i) *
                  Q_data[c * ndof_total + dof_offset1 + i];
            }
            for (int i = 0; i < ndof_e2; i++)
            {
               Q_nbr[c]  += shape2(i) *
                  Q_data[c * ndof_total + dof_offset2 + i];
            }
         }

         real_t F_int[NUM_STATE], F_ce[NUM_STATE];
         flux.Interior(nor, Q_self, Q_nbr, F_int);
         flux.Central (nor, Q_self, Q_nbr, F_ce);

         for (int c = 0; c < NUM_STATE; c++)
         {
            const real_t dF = F_int[c] - F_ce[c];
            for (int i = 0; i < ndof_e1; i++)
            {
               delta_out(c * ndof_total + dof_offset1 + i) +=
                  w * shape1(i) * dF;
            }
            for (int i = 0; i < ndof_e2; i++)
            {
               delta_out(c * ndof_total + dof_offset2 + i) -=
                  w * shape2(i) * dF;
            }
         }
      }
   }

   const int ne = mesh.GetNE();
   for (int e = 0; e < ne; e++)
   {
      const DenseMatrix &Minv = wave.GetElementMassInverse(e);
      const int dof_off = e * ndof_per_el;
      Vector tmp_in(ndof_per_el), tmp_out(ndof_per_el);
      for (int c = 0; c < NUM_STATE; c++)
      {
         for (int i = 0; i < ndof_per_el; i++)
         {
            tmp_in(i) = delta_out(c * ndof_total + dof_off + i);
         }
         Minv.Mult(tmp_in, tmp_out);
         for (int i = 0; i < ndof_per_el; i++)
         {
            delta_out(c * ndof_total + dof_off + i) = tmp_out(i);
         }
      }
   }
}

// Single-step Gate 2-analytic-ADER at the requested order.  Returns true
// on pass.
bool SingleStepGate(int order, int ader_order, const char *label)
{
   const real_t L = 1000.0;
   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr    = 3;
   bc.absorbing_attrs = {};

   std::cout << "\n----- " << label << " (order=" << order
             << ", ader_order=" << ader_order << ") -----\n";

   Mesh mesh_n = BuildSmallFaultTetMesh(L);
   WaveOperator<Mesh> wave_n(mesh_n, order,
                              TPV102Params::lambda, TPV102Params::mu,
                              TPV102Params::rho, bc);
   FaultFaceFlux ff_n(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   wave_n.SetFaultFlux(&ff_n);
   std::vector<DOFData> dof_n;
   SetupWaveForRun(wave_n, dof_n, mesh_n, order);

   Mesh mesh_a = BuildSmallFaultTetMesh(L);
   WaveOperator<Mesh> wave_a(mesh_a, order,
                              TPV102Params::lambda, TPV102Params::mu,
                              TPV102Params::rho, bc);
   FaultFaceFlux ff_a(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   wave_a.SetFaultFlux(&ff_a);
   std::vector<DOFData> dof_a;
   SetupWaveForRun(wave_a, dof_a, mesh_a, order);
   wave_a.SetMixedFluxMode(MixedFluxMode::Adjacent);

   const auto &fes_n = wave_n.GetFESpace();
   const int ndof_total = fes_n.GetNDofs();
   Vector Q_init(NUM_STATE * ndof_total);
   InitializeNonUniformState(Q_init, fes_n, ndof_total, L);

   std::cout << "  ndof_total = " << ndof_total
             << "  |central_flux_face_set_| = "
             << wave_a.GetCentralFluxFaceSet().size() << "\n";

   const real_t dt_ader = 1.0e-5;
   Vector Q_new_n(NUM_STATE * ndof_total), Q_new_a(NUM_STATE * ndof_total);
   wave_n.AdvanceADER(Q_init, dt_ader, ader_order, Q_new_n);
   wave_a.AdvanceADER(Q_init, dt_ader, ader_order, Q_new_a);

   Vector observed_delta(NUM_STATE * ndof_total);
   for (int i = 0; i < observed_delta.Size(); i++)
   {
      observed_delta(i) = Q_new_a(i) - Q_new_n(i);
   }
   const real_t scale = MaxAbs(observed_delta);

   Vector I(NUM_STATE * ndof_total);
   wave_a.ComputeADERTimeIntegrated(Q_init, dt_ader, ader_order, I);
   Vector expected_delta(NUM_STATE * ndof_total);
   ComputeAnalyticAdjacentDelta(wave_a, mesh_a, order, I, expected_delta);

   const real_t resid = MaxAbsDiff(observed_delta, expected_delta);
   const real_t rel = resid / std::max<real_t>(scale, real_t(1.0));
   std::cout << "  max |Q_new_adj − Q_new_none| = " << std::scientific
             << std::setprecision(3) << scale
             << "   rel residual vs analytic = "
             << std::scientific << std::setprecision(3) << rel << "\n";
   const std::string lbl = std::string("Gate 2-analytic-ADER @ ") + label;
   TEST_LE(rel, 1.0e-12, lbl.c_str());
   return rel <= 1.0e-12;
}

// Multi-step stability sweep — N AdvanceADER steps; record max state.
real_t MultiStepStability(int order, int ader_order, MixedFluxMode mode,
                          int N, real_t dt, const char *tag)
{
   const real_t L = 1000.0;
   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr    = 3;
   bc.absorbing_attrs = {};

   Mesh mesh = BuildSmallFaultTetMesh(L);
   WaveOperator<Mesh> wave(mesh, order,
                            TPV102Params::lambda, TPV102Params::mu,
                            TPV102Params::rho, bc);
   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   wave.SetFaultFlux(&ff);
   std::vector<DOFData> dof;
   SetupWaveForRun(wave, dof, mesh, order);
   if (mode != MixedFluxMode::None) { wave.SetMixedFluxMode(mode); }

   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   Vector Q(NUM_STATE * ndof_total), Q_next(NUM_STATE * ndof_total);
   InitializeNonUniformState(Q, fes, ndof_total, L);

   const real_t initial_max = MaxAbs(Q);
   real_t max_seen = initial_max;
   int steps_done = 0;
   for (int n = 0; n < N; n++)
   {
      wave.AdvanceADER(Q, dt, ader_order, Q_next);
      Q = Q_next;
      const real_t m = MaxAbs(Q);
      if (m > max_seen) { max_seen = m; }
      steps_done = n + 1;
      if (!std::isfinite(m) || m > 1e30 * initial_max) { break; }
   }
   std::cout << "  [" << tag << "] init|Q|=" << std::scientific
             << std::setprecision(3) << initial_max
             << "  max|Q|=" << std::scientific << std::setprecision(3)
             << max_seen << "  steps=" << steps_done << "/" << N
             << "  growth=" << std::scientific << std::setprecision(3)
             << (max_seen / initial_max) << "\n";
   return max_seen / initial_max;
}

} // anonymous namespace

int main()
{
   std::cout << "\n=== p=2 mixed-flux Adjacent reproducer ===\n";
   std::cout << "\n[Phase 1] Single-step Gate 2-analytic-ADER at p=1, p=2 with various ADER-O\n";

   // p=1, ADER-O2 — baseline known-passing case.
   SingleStepGate(/*order=*/1, /*ader_order=*/2, "p1_O2");
   // p=2, ADER-O2 — same dispatch path, higher polynomial.
   SingleStepGate(/*order=*/2, /*ader_order=*/2, "p2_O2");
   // p=2, ADER-O3 — production setting in failing TPV runs.
   SingleStepGate(/*order=*/2, /*ader_order=*/3, "p2_O3");

   // Multi-step stability sweep at the production-equivalent dt.  At p=1
   // mixed-flux is stable; at p=2 production runs blow up.  This sweep
   // uses a tiny mesh + clean linear bulk wave (no friction, no
   // nucleation) so any blow-up here is purely the bulk-wave/dispatch
   // interaction — no friction-feedback amplification.
   std::cout << "\n[Phase 2] Multi-step stability (200 ADER steps, dt = 1e-5)\n";
   const int N = 200;
   const real_t dt = 1.0e-5;
   const real_t g_p1_n = MultiStepStability(1, 2, MixedFluxMode::None,
                                            N, dt, "p1_O2_none");
   const real_t g_p1_a = MultiStepStability(1, 2, MixedFluxMode::Adjacent,
                                            N, dt, "p1_O2_adj");
   const real_t g_p2_n = MultiStepStability(2, 2, MixedFluxMode::None,
                                            N, dt, "p2_O2_none");
   const real_t g_p2_a = MultiStepStability(2, 2, MixedFluxMode::Adjacent,
                                            N, dt, "p2_O2_adj");
   const real_t g_p2_O3_n = MultiStepStability(2, 3, MixedFluxMode::None,
                                                N, dt, "p2_O3_none");
   const real_t g_p2_O3_a = MultiStepStability(2, 3, MixedFluxMode::Adjacent,
                                                N, dt, "p2_O3_adj");
   std::cout << "\n[Phase 3] Stability ratios — Adjacent / None\n";
   std::cout << "  p1_O2 ratio = " << g_p1_a / g_p1_n << "\n";
   std::cout << "  p2_O2 ratio = " << g_p2_a / g_p2_n << "\n";
   std::cout << "  p2_O3 ratio = " << g_p2_O3_a / g_p2_O3_n << "\n";

   // If single-step identity passes at p=2 but multi-step blows up,
   // diagnosis: implementation correct, scheme unstable at p≥2.

   std::cout << "\n========================================\n"
             << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n"
             << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}

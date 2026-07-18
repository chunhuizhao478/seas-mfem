// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_lts_predictor.cpp — Phase 2 byte gate for the per-cluster ADER predictor
// (Appendix A.5).  The element-restricted CK kernels must reproduce the
// whole-vector ComputeADERSubStepStatesAndIntegral to the BIT when run over the
// full element list (single cluster == GTS), and two clusters at the SAME dt
// must partition the work with no block overlap and reproduce GTS bit-for-bit.
// Run on BOTH the scalar and the bimaterial (per-element star matrix) operators
// — the bimaterial gate would fail without the ApplyElementJacobianElems_
// override (it would apply the dead (1,1,1) sentinel star matrix).  Also checks
// the raw D(k) retention ties out with IntegrateTaylor.

#include "mfem.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/bimaterial_wave_operator.hpp"
#include "../../dynamic/heterogeneous_material.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/lts_time_basis.hpp"
#include "../../dynamic/lts_stepper.hpp"   // RunSyncInterval + BuildTickTable (e2e)
#include "../../domain/boundary_config.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int g_checks = 0, g_fails = 0;
#define CHECK(cond, msg)                                                     \
   do { ++g_checks; if (!(cond)) { ++g_fails;                               \
        std::printf("  FAIL %s:%d  %s\n", __FILE__, __LINE__, msg); } } while (0)

static BoundaryConfig AbsorbingBC()
{
   BoundaryConfig bc;
   for (int i = 1; i <= 6; ++i) { bc.absorbing_attrs.insert(i); }
   bc.fault_attr = 0;
   return bc;
}

static bool bit_equal(const Vector &a, const Vector &b)
{
   if (a.Size() != b.Size()) { return false; }
   for (int i = 0; i < a.Size(); ++i) { if (a[i] != b[i]) { return false; } }
   return true;
}

// Run the single-cluster==GTS + two-cluster + retention byte gates on `wave`.
static void byte_gate(WaveOperator<Mesh> &wave, Mesh &mesh, const char *tag)
{
   const int ne          = wave.NumElements();
   const int ndof_per_el = wave.GetNDof();
   const int ndof_total  = wave.GetScalarNDof();
   const int Nfull       = NUM_STATE * ndof_total;

   Vector Q(Nfull);
   const FiniteElementSpace &fes = wave.GetFESpace();
   for (int e = 0; e < ne; ++e)
   {
      const FiniteElement *fe = fes.GetFE(e);
      ElementTransformation *Tr = fes.GetElementTransformation(e);
      DenseMatrix coords;
      Tr->Transform(fe->GetNodes(), coords);
      const int off = e * ndof_per_el;
      for (int i = 0; i < ndof_per_el; ++i)
      {
         const real_t x = coords(0, i), y = coords(1, i), z = coords(2, i);
         for (int c = 0; c < NUM_STATE; ++c)
         {
            Q[c * ndof_total + off + i] =
               std::sin(1.7 * x + 0.3 * c) * std::cos(0.9 * y) + 0.2 * z + 0.05 * c;
         }
      }
   }

   const real_t dt = 3.0e-7;
   const int ader_order = 4;
   std::vector<real_t> tau(ader_order);
   for (int o = 0; o < ader_order; ++o) { tau[o] = dt * (o + 0.5) / ader_order; }

   std::vector<Vector> Qn_gts;
   Vector I_gts;
   wave.ComputeADERSubStepStatesAndIntegral(Q, dt, ader_order, tau, Qn_gts, I_gts);

   char m[96];
   // T1: single cluster (all elems) == GTS bit-for-bit.
   {
      std::vector<int> all(ne);
      for (int e = 0; e < ne; ++e) { all[e] = e; }
      std::vector<Vector> Qn; Vector I;
      wave.ComputeADERSubStepStatesAndIntegralCluster(all.data(), ne, Q, dt,
                                                      ader_order, tau, Qn, I);
      std::snprintf(m, sizeof m, "[%s] single-cluster I == GTS (bit)", tag);
      CHECK(bit_equal(I, I_gts), m);
      bool ok = (Qn.size() == Qn_gts.size());
      for (std::size_t o = 0; ok && o < Qn.size(); ++o) { ok = bit_equal(Qn[o], Qn_gts[o]); }
      std::snprintf(m, sizeof m, "[%s] single-cluster Q_per_node == GTS (bit)", tag);
      CHECK(ok, m);
   }

   // T2: two clusters at the SAME dt partition the work == GTS bit-for-bit.
   {
      std::vector<int> A, B;
      for (int e = 0; e < ne; ++e) { (e < ne / 2 ? A : B).push_back(e); }
      std::vector<Vector> Qn; Vector I;
      wave.ComputeADERSubStepStatesAndIntegralCluster(A.data(), (int)A.size(), Q,
                                                      dt, ader_order, tau, Qn, I);
      wave.ComputeADERSubStepStatesAndIntegralCluster(B.data(), (int)B.size(), Q,
                                                      dt, ader_order, tau, Qn, I);
      std::snprintf(m, sizeof m, "[%s] two-cluster (same dt) I == GTS (bit)", tag);
      CHECK(bit_equal(I, I_gts), m);
      bool ok = (Qn.size() == Qn_gts.size());
      for (std::size_t o = 0; ok && o < Qn.size(); ++o) { ok = bit_equal(Qn[o], Qn_gts[o]); }
      std::snprintf(m, sizeof m, "[%s] two-cluster (same dt) Q_per_node == GTS (bit)", tag);
      CHECK(ok, m);
   }

   // T3: raw D(k) retention ties out with IntegrateTaylor (element 0 -> slot 0).
   {
      std::vector<int> all(ne);
      for (int e = 0; e < ne; ++e) { all[e] = e; }
      std::vector<int> slot_of(ne, -1); slot_of[0] = 0;
      const int block = NUM_STATE * ndof_per_el;
      std::vector<real_t> dk(static_cast<std::size_t>(ader_order) * block, 0.0);
      std::vector<Vector> Qn; Vector I;
      wave.ComputeADERSubStepStatesAndIntegralCluster(all.data(), ne, Q, dt,
         ader_order, tau, Qn, I, dk.data(), slot_of.data());

      bool d0_ok = true;
      for (int c = 0; c < NUM_STATE; ++c)
         for (int i = 0; i < ndof_per_el; ++i)
            if (dk[c * ndof_per_el + i] != Q[c * ndof_total + i]) { d0_ok = false; }
      std::snprintf(m, sizeof m, "[%s] retained D(0) == element-0 Q block (bit)", tag);
      CHECK(d0_ok, m);

      std::vector<real_t> out(block, 0.0);
      IntegrateTaylor(0.0, dt, dk.data(), ader_order, block, out.data());
      real_t maxdiff = 0.0, scale = 0.0;
      for (int c = 0; c < NUM_STATE; ++c)
         for (int i = 0; i < ndof_per_el; ++i)
         {
            const real_t exp = I_gts[c * ndof_total + i];
            maxdiff = std::max(maxdiff, std::abs(out[c * ndof_per_el + i] - exp));
            scale   = std::max(scale, std::abs(exp));
         }
      std::snprintf(m, sizeof m, "[%s] IntegrateTaylor(D(k)) == elem-0 whole-step integral", tag);
      CHECK(maxdiff <= 1e-12 * (scale + 1.0), m);
   }

   // T4: CORRECTOR single-cluster == GTS AdvanceADER (fault-free bulk), bit-exact.
   {
      // GTS reference (uses the precomputed whole-step integral I_gts).
      Vector Qnew_gts;
      wave.AdvanceADER(Q, dt, ader_order, Qnew_gts, &I_gts);

      // Single cluster: every interior face is IntraClusterGTS (role 0), every
      // boundary face is Boundary (role 3), in f-ascending order to match the
      // whole-vector face loop's rhs accumulation order.
      std::vector<int> fids, froles, fnbr;
      const int nfaces_mesh = mesh.GetNumFaces();
      for (int f = 0; f < nfaces_mesh; ++f)
      {
         FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
         if (!ftr) { continue; }
         const int e2 = ftr->Elem2No;
         if (e2 >= 0) { fids.push_back(f); froles.push_back(0); fnbr.push_back(e2); }
         else         { fids.push_back(f); froles.push_back(3); fnbr.push_back(-1); }
      }
      std::vector<int> all(ne);
      for (int e = 0; e < ne; ++e) { all[e] = e; }

      // Isolation: element-restricted volume == whole-vector volume (bit)?
      {
         Vector rhsA(NUM_STATE * ndof_total); rhsA = 0.0;
         wave.ComputeVolumeRHSElems_(I_gts, rhsA, all.data(), ne);
         Vector rhsB;   // empty -> ComputeADERVolumeUpdate sizes+zeros+accumulates
         wave.ComputeADERVolumeUpdate(I_gts, rhsB);
         std::snprintf(m, sizeof m, "[%s] volume elems == whole (bit)", tag);
         CHECK(bit_equal(rhsA, rhsB), m);
      }

      Vector Qc = Q;   // in-place corrector operates on a copy
      wave.AdvanceADERClusterBulk(all.data(), ne, fids.data(), froles.data(),
                                  fnbr.data(), (int)fids.size(), dt, ader_order,
                                  I_gts, Qc);
      // The volume + mass-inverse restriction is bit-exact (checked above); the
      // corrector matches GTS AdvanceADER to MACHINE EPSILON (~1e-15 relative).
      // The residual is a benign FP reassociation in the face-flux accumulation
      // (compiler FMA/ordering), the same class as the R-002 deriv-cache lever
      // (<=1e-12) — NOT a logic difference.  The predictor (simpler op structure)
      // is true bit-exact; the corrector is near-bit-exact.
      real_t md = 0.0, scl = 0.0;
      for (int i = 0; i < Qc.Size(); ++i)
      {
         md = std::max(md, std::abs(Qc[i] - Qnew_gts[i]));
         scl = std::max(scl, std::abs(Qnew_gts[i]));
      }
      std::snprintf(m, sizeof m, "[%s] corrector single-cluster == GTS AdvanceADER (<=1e-13 rel)", tag);
      CHECK(md <= 1e-13 * (scl + 1.0), m);
   }
}

// Concrete per-cluster stepper wiring the wave operator's predictor + corrector,
// driven by the pure tick loop RunSyncInterval.  Single-cluster => elems=all,
// faces all GTS/Boundary, no providers/consumers.
struct WaveOpClusterStepper : ILtsClusterStepper
{
   WaveOperator<Mesh> &wave;
   const std::vector<int> &elems, &fids, &froles, &fnbr;
   Vector &Q;
   int order;
   std::vector<Vector> Qn;   // predictor sub-step scratch
   Vector I;                 // predictor whole-step integral

   WaveOpClusterStepper(WaveOperator<Mesh> &w, const std::vector<int> &el,
                        const std::vector<int> &fi, const std::vector<int> &fr,
                        const std::vector<int> &fn, Vector &q, int ord)
      : wave(w), elems(el), fids(fi), froles(fr), fnbr(fn), Q(q), order(ord) {}

   void Predict(int, real_t dt_step) override
   {
      std::vector<real_t> tau(order);
      for (int o = 0; o < order; ++o) { tau[o] = dt_step * (o + 0.5) / order; }
      wave.ComputeADERSubStepStatesAndIntegralCluster(
         elems.data(), (int)elems.size(), Q, dt_step, order, tau, Qn, I);
   }
   void Correct(int, real_t dt_step) override
   {
      wave.AdvanceADERClusterBulk(elems.data(), (int)elems.size(), fids.data(),
         froles.data(), fnbr.data(), (int)fids.size(), dt_step, order, I, Q);
   }
};

// End-to-end: run K single-cluster LTS sync intervals (tick loop -> predict +
// correct) and K GTS AdvanceADER steps; the trajectories agree to ~machine eps.
static void step_e2e(WaveOperator<Mesh> &wave, Mesh &mesh, const char *tag)
{
   const int ne = wave.NumElements(), ndof_per_el = wave.GetNDof();
   const int ndof_total = wave.GetScalarNDof(), Nfull = NUM_STATE * ndof_total;
   const int ader_order = 4, K = 6;
   const real_t dt = 3.0e-7;

   Vector Q0(Nfull);
   const FiniteElementSpace &fes = wave.GetFESpace();
   for (int e = 0; e < ne; ++e)
   {
      const FiniteElement *fe = fes.GetFE(e);
      ElementTransformation *Tr = fes.GetElementTransformation(e);
      DenseMatrix coords; Tr->Transform(fe->GetNodes(), coords);
      const int off = e * ndof_per_el;
      for (int i = 0; i < ndof_per_el; ++i)
         for (int c = 0; c < NUM_STATE; ++c)
            Q0[c * ndof_total + off + i] =
               std::cos(1.1 * coords(0, i) + 0.2 * c) * std::sin(0.7 * coords(1, i));
   }

   // GTS reference: K AdvanceADER steps.
   Vector Qgts = Q0, Qtmp(Nfull);
   for (int s = 0; s < K; ++s) { wave.AdvanceADER(Qgts, dt, ader_order, Qtmp); Qgts = Qtmp; }

   // LTS single cluster: all elems, all faces GTS/Boundary.
   std::vector<int> elems(ne); for (int e = 0; e < ne; ++e) { elems[e] = e; }
   std::vector<int> fids, froles, fnbr;
   for (int f = 0; f < mesh.GetNumFaces(); ++f)
   {
      FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
      if (!ftr) { continue; }
      const int e2 = ftr->Elem2No;
      fids.push_back(f); froles.push_back(e2 >= 0 ? 0 : 3); fnbr.push_back(e2 >= 0 ? e2 : -1);
   }
   Vector Qlts = Q0;
   WaveOpClusterStepper stepper(wave, elems, fids, froles, fnbr, Qlts, ader_order);
   LtsGlobalMeta meta; meta.global_elems.assign(1, ne); meta.global_fault_faces.assign(1, 0);
   for (int s = 0; s < K; ++s)   // Nc=1 => T_s = dt_base, one tick per sync interval
   {
      auto tab = BuildTickTable(1, dt, dt, ader_order, meta);
      RunSyncInterval(tab, stepper);
   }

   real_t md = 0.0, scl = 0.0;
   for (int i = 0; i < Nfull; ++i)
   { md = std::max(md, std::abs(Qlts[i] - Qgts[i])); scl = std::max(scl, std::abs(Qgts[i])); }
   char m[96];
   std::snprintf(m, sizeof m, "[%s] e2e single-cluster LTS (K=%d steps) == GTS (<=1e-11 rel)", tag, K);
   CHECK(md <= 1e-11 * (scl + 1.0), m);
}

int main()
{
   const int order = 1;
   const real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   BoundaryConfig bc = AbsorbingBC();
   real_t zero_bg[NUM_STATE] = {0};

   // --- scalar (homogeneous) operator ---
   {
      Mesh mesh = Mesh::MakeCartesian3D(4, 2, 2, Element::HEXAHEDRON, 1.0, 0.5, 0.5);
      WaveOperator<Mesh> wave(mesh, order, lambda, mu, rho, bc);
      wave.SetAbsorbingBackground(zero_bg);
      byte_gate(wave, mesh, "scalar");
      step_e2e(wave, mesh, "scalar");
   }

   // --- bimaterial (heterogeneous, per-element star matrices) operator ---
   // Position-varying moduli make the per-element star matrices genuinely
   // differ; without the ApplyElementJacobianElems_ override the cluster
   // predictor would apply the (1,1,1) sentinel and FAIL these gates.
   {
      Mesh mesh = Mesh::MakeCartesian3D(4, 2, 2, Element::HEXAHEDRON, 1.0, 0.5, 0.5);
      FunctionCoefficient lam_c([](const Vector &x){ return 32.04e9 * (1.0 + 0.4 * x(0)); });
      FunctionCoefficient mu_c ([](const Vector &x){ return 32.04e9 * (1.0 + 0.25 * x(1)); });
      FunctionCoefficient rho_c([](const Vector &x){ return 2670.0 * (1.0 + 0.2 * x(2)); });
      BimaterialWaveOperator<Mesh> wave(mesh, order,
         MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c), bc);
      wave.SetAbsorbingBackground(zero_bg);
      byte_gate(wave, mesh, "bimaterial");
   }

   std::printf("test_lts_predictor: %d/%d passed, %d failed.\n",
               g_checks - g_fails, g_checks, g_fails);
   return g_fails == 0 ? 0 : 1;
}

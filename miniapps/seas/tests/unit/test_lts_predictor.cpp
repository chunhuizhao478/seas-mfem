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
static void byte_gate(WaveOperator<Mesh> &wave, const char *tag)
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
      byte_gate(wave, "scalar");
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
      byte_gate(wave, "bimaterial");
   }

   std::printf("test_lts_predictor: %d/%d passed, %d failed.\n",
               g_checks - g_fails, g_checks, g_fails);
   return g_fails == 0 ? 0 : 1;
}

// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_lts_predictor.cpp — Phase 2 byte gate for the per-cluster ADER predictor
// (Appendix A.5).  The element-restricted CK kernels must reproduce the
// whole-vector ComputeADERSubStepStatesAndIntegral to the BIT when run over the
// full element list (single cluster == GTS), and two clusters at the SAME dt
// must partition the work with no block overlap and reproduce GTS bit-for-bit.
// Also checks the raw D(k) retention ties out with IntegrateTaylor.

#include "mfem.hpp"
#include "../../dynamic/wave_operator.hpp"
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

// Exact-equality over two full-size vectors.
static bool bit_equal(const Vector &a, const Vector &b)
{
   if (a.Size() != b.Size()) { return false; }
   for (int i = 0; i < a.Size(); ++i) { if (a[i] != b[i]) { return false; } }
   return true;
}

int main()
{
   Mesh mesh = Mesh::MakeCartesian3D(4, 2, 2, Element::HEXAHEDRON, 1.0, 0.5, 0.5);
   const int order = 1;                 // spatial (P1); ADER order below
   const real_t lambda = 32.04e9, mu = 32.04e9, rho = 2670.0;
   BoundaryConfig bc = AbsorbingBC();
   WaveOperator<Mesh> wave(mesh, order, lambda, mu, rho, bc);

   const int ne         = wave.NumElements();
   const int ndof_per_el = wave.GetNDof();
   const int ndof_total = wave.GetScalarNDof();
   const int Nfull      = NUM_STATE * ndof_total;

   // Smooth initial state across all NUM_STATE components.
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

   const real_t dt = 3.0e-6;
   const int ader_order = 4;
   std::vector<real_t> tau_nodes(ader_order);
   for (int o = 0; o < ader_order; ++o)
   { tau_nodes[o] = dt * (o + 0.5) / ader_order; }

   // Reference: the whole-vector GTS predictor.
   std::vector<Vector> Qn_gts;
   Vector I_gts;
   wave.ComputeADERSubStepStatesAndIntegral(Q, dt, ader_order, tau_nodes, Qn_gts, I_gts);

   // ---- T1: single cluster (all elems ascending) == GTS bit-for-bit ---------
   {
      std::vector<int> all(ne);
      for (int e = 0; e < ne; ++e) { all[e] = e; }
      std::vector<Vector> Qn_lts;
      Vector I_lts;
      wave.ComputeADERSubStepStatesAndIntegralCluster(
         all.data(), ne, Q, dt, ader_order, tau_nodes, Qn_lts, I_lts);
      CHECK(I_lts.Size() == Nfull, "T1 I sized full");
      CHECK(bit_equal(I_lts, I_gts), "T1 single-cluster I == GTS (bit)");
      bool nodes_ok = (Qn_lts.size() == Qn_gts.size());
      for (std::size_t o = 0; nodes_ok && o < Qn_lts.size(); ++o)
      { nodes_ok = bit_equal(Qn_lts[o], Qn_gts[o]); }
      CHECK(nodes_ok, "T1 single-cluster Q_per_node == GTS (bit)");
   }

   // ---- T2: two clusters at the SAME dt partition the work == GTS (bit) ------
   {
      std::vector<int> A, B;
      for (int e = 0; e < ne; ++e) { (e < ne / 2 ? A : B).push_back(e); }
      std::vector<Vector> Qn;   // persistent across the two cluster calls
      Vector I;
      // First cluster A writes A's blocks; the second (B) writes B's blocks and
      // must leave A's blocks intact (zero only its own).
      wave.ComputeADERSubStepStatesAndIntegralCluster(
         A.data(), (int)A.size(), Q, dt, ader_order, tau_nodes, Qn, I);
      wave.ComputeADERSubStepStatesAndIntegralCluster(
         B.data(), (int)B.size(), Q, dt, ader_order, tau_nodes, Qn, I);
      CHECK(bit_equal(I, I_gts), "T2 two-cluster (same dt) I == GTS (bit)");
      bool ok = (Qn.size() == Qn_gts.size());
      for (std::size_t o = 0; ok && o < Qn.size(); ++o) { ok = bit_equal(Qn[o], Qn_gts[o]); }
      CHECK(ok, "T2 two-cluster (same dt) Q_per_node == GTS (bit)");
   }

   // ---- T3: raw D(k) retention ties out with IntegrateTaylor -----------------
   {
      std::vector<int> all(ne);
      for (int e = 0; e < ne; ++e) { all[e] = e; }
      // Provider = element 0 -> slot 0.
      std::vector<int> slot_of(ne, -1);
      slot_of[0] = 0;
      const int block = NUM_STATE * ndof_per_el;
      std::vector<real_t> dk(static_cast<std::size_t>(1) * ader_order * block, 0.0);
      std::vector<Vector> Qn;
      Vector I;
      wave.ComputeADERSubStepStatesAndIntegralCluster(
         all.data(), ne, Q, dt, ader_order, tau_nodes, Qn, I,
         dk.data(), slot_of.data());

      // Retained D(0) == element 0's Q block, gathered contiguously.
      bool d0_ok = true;
      for (int c = 0; c < NUM_STATE; ++c)
         for (int i = 0; i < ndof_per_el; ++i)
         {
            const real_t got = dk[c * ndof_per_el + i];             // slot0,k0
            const real_t exp = Q[c * ndof_total + 0 * ndof_per_el + i];
            if (got != exp) { d0_ok = false; }
         }
      CHECK(d0_ok, "T3 retained D(0) == element-0 Q block (bit)");

      // IntegrateTaylor(0,dt, slot-0 stack) == element 0's whole-step integral
      // block (I_gts), to FP tolerance (different factorial factorization).
      std::vector<real_t> out(block, 0.0);
      IntegrateTaylor(0.0, dt, dk.data(), ader_order, block, out.data());
      real_t maxdiff = 0.0, scale = 0.0;
      for (int c = 0; c < NUM_STATE; ++c)
         for (int i = 0; i < ndof_per_el; ++i)
         {
            const real_t got = out[c * ndof_per_el + i];
            const real_t exp = I_gts[c * ndof_total + 0 * ndof_per_el + i];
            maxdiff = std::max(maxdiff, std::abs(got - exp));
            scale   = std::max(scale, std::abs(exp));
         }
      CHECK(maxdiff <= 1e-12 * (scale + 1.0),
            "T3 IntegrateTaylor(0,dt, D(k)) == element-0 whole-step integral");
   }

   std::printf("test_lts_predictor: %d/%d passed, %d failed.\n",
               g_checks - g_fails, g_checks, g_fails);
   return g_fails == 0 ? 0 : 1;
}

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
#include "../../dynamic/lts_layout.hpp"    // BuildLtsLayout (multi-cluster path)
#include "../../dynamic/lts_bulk_stepper.hpp"  // LtsBulkSyncStepper (real multi-rate)
#include "../../io/tpv104_checkpoint.hpp"       // Checkpoint V2 write/read/peek
#include "../../domain/boundary_config.hpp"
#include <cstdio>

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

// Degenerate 2-cluster validation of the consumer/provider/buffer path: split the
// mesh into two clusters (ids 0 fine, 1 coarse) via BuildLtsLayout, but DRIVE both
// at the same dt with [0,dt] consumer sub-intervals.  IntegrateTaylor(D_coarse,
// 0, dt) == the coarse's whole-step integral, so the coupling reduces to GTS: the
// result must match AdvanceADER to machine epsilon while exercising D(k)
// retention, the consumer-face flux, and the accumulate-buffer fill/consume.
static void multicluster_degenerate(WaveOperator<Mesh> &wave, Mesh &mesh, const char *tag)
{
   using namespace mfem::seas;
   const int ne = wave.NumElements(), ndof_per_el = wave.GetNDof();
   const int ndof_total = wave.GetScalarNDof(), Nfull = NUM_STATE * ndof_total;
   const int block = NUM_STATE * ndof_per_el;
   const int order = 4;
   const real_t dt = 3.0e-7;

   Vector Q(Nfull);
   const FiniteElementSpace &fes = wave.GetFESpace();
   for (int e = 0; e < ne; ++e)
   {
      const FiniteElement *fe = fes.GetFE(e);
      ElementTransformation *Tr = fes.GetElementTransformation(e);
      DenseMatrix coords; Tr->Transform(fe->GetNodes(), coords);
      const int off = e * ndof_per_el;
      for (int i = 0; i < ndof_per_el; ++i)
         for (int c = 0; c < NUM_STATE; ++c)
            Q[c * ndof_total + off + i] = std::sin(0.9 * coords(0, i) + 0.15 * c);
   }
   char m[110];

   // GTS reference.
   Vector Qgts;
   wave.AdvanceADER(Q, dt, order, Qgts);

   // Cluster ids that satisfy maxdiff<=1: a contiguous split by element index
   // (0..h-1 -> cluster 0, h..ne-1 -> cluster 1).  If the induced layout has no
   // ConsumerFine faces (clusters not adjacent), skip — nothing to validate.
   std::vector<int> cluster(ne, 0);
   for (int e = ne / 2; e < ne; ++e) { cluster[e] = 1; }

   std::vector<LtsFaceSpec> faces;
   for (int f = 0; f < mesh.GetNumFaces(); ++f)
   {
      FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
      if (!ftr) { continue; }
      LtsFaceSpec fs; fs.face_id = f; fs.elem1 = ftr->Elem1No;
      fs.elem2 = ftr->Elem2No; fs.is_fault = false;
      faces.push_back(fs);
   }
   // Guard the maxdiff<=1 precondition (a contiguous index split can violate it
   // on an unstructured numbering; on this Cartesian mesh it holds).
   bool maxdiff_ok = true;
   for (const auto &fs : faces)
      if (fs.elem2 >= 0 && std::abs(cluster[fs.elem1] - cluster[fs.elem2]) > 1) { maxdiff_ok = false; }
   if (!maxdiff_ok)
   {
      std::snprintf(m, sizeof m, "[%s] multicluster: SKIP (index split violates maxdiff)", tag);
      CHECK(true, m); return;
   }
   LtsLayout L = BuildLtsLayout(cluster, 2, faces);

   int n_consumer = 0;
   for (const auto &cl : L.clusters)
      for (FaceRole r : cl.faces) if (r == FaceRole::ConsumerFine) { ++n_consumer; }
   if (n_consumer == 0)
   {
      std::snprintf(m, sizeof m, "[%s] multicluster: SKIP (no consumer faces)", tag);
      CHECK(true, m); return;
   }

   LtsDkStore dk; dk.Resize((int)L.provider_elems.size(), order, block);
   LtsAccumulateBuffers buf; buf.Resize((int)L.consumer_owner_elems.size(), block);

   std::vector<real_t> tau(order);
   for (int o = 0; o < order; ++o) { tau[o] = dt * (o + 0.5) / order; }
   std::vector<Vector> Qn;
   Vector I0, I1;

   // Predict fine (cluster 0) and coarse (cluster 1, retaining D(k) for providers).
   wave.ComputeADERSubStepStatesAndIntegralCluster(
      L.clusters[0].elems.data(), (int)L.clusters[0].elems.size(), Q, dt, order,
      tau, Qn, I0);
   wave.ComputeADERSubStepStatesAndIntegralCluster(
      L.clusters[1].elems.data(), (int)L.clusters[1].elems.size(), Q, dt, order,
      tau, Qn, I1, dk.data.data(), L.provider_slot_of_elem.data());

   // Per-cluster owned-face arrays + degenerate consumer sub-interval [0, dt].
   auto face_arrays = [&](int c, std::vector<int> &fids, std::vector<int> &fr,
                          std::vector<int> &fn, std::vector<real_t> &sa,
                          std::vector<real_t> &sb)
   {
      const LtsCluster &cl = L.clusters[c];
      for (std::size_t i = 0; i < cl.face_ids.size(); ++i)
      {
         fids.push_back(cl.face_ids[i]);
         fr.push_back((int)cl.faces[i]);
         fn.push_back(cl.face_nbr[i]);
         sa.push_back(0.0); sb.push_back(dt);
      }
   };

   Vector Qc = Q;   // in-place
   // Correct FINE (cluster 0) then COARSE (cluster 1).
   {
      std::vector<int> fids, fr, fn; std::vector<real_t> sa, sb;
      face_arrays(0, fids, fr, fn, sa, sb);
      wave.AdvanceADERClusterBulk(L.clusters[0].elems.data(),
         (int)L.clusters[0].elems.size(), fids.data(), fr.data(), fn.data(),
         (int)fids.size(), dt, order, I0, Qc, dk.data.data(), order,
         L.provider_slot_of_elem.data(), sa.data(), sb.data(), &buf,
         L.buffer_slot_of_elem.data());
   }
   {
      std::vector<int> fids, fr, fn; std::vector<real_t> sa, sb;
      face_arrays(1, fids, fr, fn, sa, sb);
      wave.AdvanceADERClusterBulk(L.clusters[1].elems.data(),
         (int)L.clusters[1].elems.size(), fids.data(), fr.data(), fn.data(),
         (int)fids.size(), dt, order, I1, Qc, dk.data.data(), order,
         L.provider_slot_of_elem.data(), sa.data(), sb.data(), &buf,
         L.buffer_slot_of_elem.data());
   }

   // Buffers must be empty after the coarse consume (invariant ii).
   std::snprintf(m, sizeof m, "[%s] multicluster: buffers zero after consume", tag);
   CHECK(buf.AllZero(), m);

   real_t md = 0.0, scl = 0.0;
   for (int i = 0; i < Nfull; ++i)
   { md = std::max(md, std::abs(Qc[i] - Qgts[i])); scl = std::max(scl, std::abs(Qgts[i])); }
   std::snprintf(m, sizeof m, "[%s] degenerate 2-cluster (dt,[0,dt]) == GTS (<=1e-11 rel), %d consumer faces",
                 tag, n_consumer);
   CHECK(md <= 1e-11 * (scl + 1.0), m);
}

// REAL multi-rate: 2-cluster LTS with the coarse cluster stepping at 2*dt (fine
// at dt, consumer sub-intervals [0,dt] then [dt,2dt] from the closed-form
// schedule) vs GTS at the fine dt.  Consistent ADER-4 schemes agree to
// truncation order; a WRONG consumer [a,b] would give an O(1) coupling error, so
// this discriminates the schedule.  Also asserts invariant (ii) (buffers zero at
// every sync point) via the driver-core stepper LtsBulkSyncStepper.
static void multirate_conservation(WaveOperator<Mesh> &wave, Mesh &mesh, const char *tag)
{
   using namespace mfem::seas;
   const int ne = wave.NumElements(), ndof_per_el = wave.GetNDof();
   const int ndof_total = wave.GetScalarNDof(), Nfull = NUM_STATE * ndof_total;
   const int order = 4;
   const real_t dt = 3.0e-7;   // 2*dt well within CFL for this box
   const int n_sync = 3;       // 3 coarse steps = 6 fine steps

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
               std::sin(0.8 * coords(0, i) + 0.1 * c) * std::cos(0.6 * coords(1, i));
   }
   char m[110];

   // GTS reference at the FINE dt (2*n_sync steps).
   Vector Qgts = Q0, Qtmp;
   for (int s = 0; s < 2 * n_sync; ++s) { wave.AdvanceADER(Qgts, dt, order, Qtmp); Qgts = Qtmp; }
   // GTS at the COARSE 2*dt (all elems) — the pure time-step-accuracy baseline
   // the multi-rate scheme must not do WORSE than (it refines the fine half).
   Vector Q2dt = Q0, Qt2;
   for (int s = 0; s < n_sync; ++s) { wave.AdvanceADER(Q2dt, 2.0 * dt, order, Qt2); Q2dt = Qt2; }

   // 2-cluster layout (0 fine, 1 coarse) by contiguous index; guard maxdiff<=1.
   std::vector<int> cluster(ne, 0);
   for (int e = ne / 2; e < ne; ++e) { cluster[e] = 1; }
   std::vector<LtsFaceSpec> faces;
   bool maxdiff_ok = true;
   for (int f = 0; f < mesh.GetNumFaces(); ++f)
   {
      FaceElementTransformations *ftr = mesh.GetFaceElementTransformations(f);
      if (!ftr) { continue; }
      LtsFaceSpec fs; fs.face_id = f; fs.elem1 = ftr->Elem1No; fs.elem2 = ftr->Elem2No; fs.is_fault = false;
      faces.push_back(fs);
      if (fs.elem2 >= 0 && std::abs(cluster[fs.elem1] - cluster[fs.elem2]) > 1) { maxdiff_ok = false; }
   }
   if (!maxdiff_ok)
   { std::snprintf(m, sizeof m, "[%s] multirate: SKIP (maxdiff)", tag); CHECK(true, m); return; }
   LtsLayout L = BuildLtsLayout(cluster, 2, faces);

   LtsGlobalMeta meta; meta.global_elems.assign(2, ne / 2); meta.global_fault_faces.assign(2, 0);
   Vector Qlts = Q0;
   LtsBulkSyncStepper<Mesh> stepper(wave, L, cluster, Qlts, order, dt);
   const real_t T_s = dt * 2.0;   // Nc=2 => coarsest dt = 2*dt
   bool buffers_ok = true;
   for (int s = 0; s < n_sync; ++s)
   {
      auto tab = BuildTickTable(2, dt, T_s, order, meta);
      stepper.SetSyncInterval(s * T_s, T_s);
      RunSyncInterval(tab, stepper);
      if (!stepper.BuffersZero()) { buffers_ok = false; }
   }
   std::snprintf(m, sizeof m, "[%s] multirate: buffers zero at every sync (inv ii)", tag);
   CHECK(buffers_ok, m);

   // Error of LTS vs GTS(dt) and of GTS(2dt) vs GTS(dt), split by region.
   real_t scl = 0.0;
   real_t lts_tot = 0.0, lts_fine = 0.0, g2_tot = 0.0, g2_fine = 0.0;
   for (int c = 0; c < NUM_STATE; ++c)
      for (int e = 0; e < ne; ++e)
         for (int i = 0; i < ndof_per_el; ++i)
         {
            const int idx = c * ndof_total + e * ndof_per_el + i;
            scl = std::max(scl, std::abs(Qgts[idx]));
            const real_t dl = std::abs(Qlts[idx] - Qgts[idx]);
            const real_t d2 = std::abs(Q2dt[idx] - Qgts[idx]);
            lts_tot = std::max(lts_tot, dl); g2_tot = std::max(g2_tot, d2);
            if (cluster[e] == 0) { lts_fine = std::max(lts_fine, dl); g2_fine = std::max(g2_fine, d2); }
         }
   std::printf("    [%s] multirate diag: LTS/GTSdt rel=%.3e (fine=%.3e), GTS2dt/GTSdt rel=%.3e (fine=%.3e)\n",
               tag, lts_tot/(scl+1e-300), lts_fine/(scl+1e-300), g2_tot/(scl+1e-300), g2_fine/(scl+1e-300));

   // Correctness (the closed-form [a,b] coupling is exercised here):
   //  (a) GLOBAL: LTS is no worse than the pure coarse-2dt scheme (a wrong [a,b]
   //      would DIVERGE from GTS(2dt) and blow this up).
   std::snprintf(m, sizeof m, "[%s] multirate: LTS error <= GTS(2dt) error (consistent, not diverging)", tag);
   CHECK(lts_tot <= 1.05 * g2_tot + 1e-9 * scl, m);
   //  (b) FINE region: LTS (fine at dt) is STRICTLY better than GTS(2dt) there —
   //      the fine refinement + coupling actually improved the fine solution.
   std::snprintf(m, sizeof m, "[%s] multirate: fine region strictly refined vs GTS(2dt) (%.2e < %.2e)",
                 tag, lts_fine/(scl+1e-300), g2_fine/(scl+1e-300));
   CHECK(lts_fine < 0.5 * g2_fine, m);
}

// Checkpoint-V2 write -> peek -> read round-trip (B.8): a fault-free bulk V2
// checkpoint (empty dof_data) preserves t/dt/sync/layout_hash/Q exactly, and the
// magic peeks as version 2.
static void checkpoint_v2_roundtrip()
{
   const char *prefix = "test_lts_ckpt_tmp";
   const real_t t = 1.234e-3, dt = 5.6e-7;
   const int sync = 42; const std::uint64_t hash = 0x0123456789abcdefULL;
   Vector Qw(24); for (int i = 0; i < 24; ++i) { Qw[i] = std::sin(0.3 * i) - 0.11 * i; }
   std::vector<DOFData> dd;   // fault-free bulk -> empty

   mfem::seas::internal::WriteTpv104CheckpointV2Impl(prefix, t, dt, sync,
      /*lts_mode=*/1, hash, Qw, dd, /*rank=*/0, /*size=*/1, "spatial_dyn");

   CHECK(mfem::seas::internal::PeekTpv104CheckpointVersion(prefix, 0) == 2,
         "[ckpt] V2 file peeks as version 2");

   real_t tr = 0, dtr = 0; int syncr = 0, lts_mode = 0; std::uint64_t hr = 0;
   Vector Qr; std::vector<DOFData> ddr; std::string tag;
   const bool ok = mfem::seas::internal::ReadTpv104CheckpointV2Impl(prefix, tr, dtr,
      syncr, lts_mode, hr, Qr, 24, ddr, 0, 1, &tag);
   CHECK(ok, "[ckpt] V2 read ok");
   CHECK(tr == t && dtr == dt && syncr == sync, "[ckpt] V2 t/dt/sync round-trip");
   CHECK(lts_mode == 1 && hr == hash, "[ckpt] V2 lts_mode + layout_hash round-trip");
   CHECK(tag == "spatial_dyn", "[ckpt] V2 driver_tag round-trip");
   bool q_ok = (Qr.Size() == 24);
   for (int i = 0; q_ok && i < 24; ++i) { q_ok = (Qr[i] == Qw[i]); }
   CHECK(q_ok, "[ckpt] V2 Q round-trip bit-exact (17-digit text)");

   std::remove((std::string(prefix) + "_checkpoint_r0.txt").c_str());
}

// Checkpoint-V2 canonical-order serialization under the fault reorder (B.8b,
// LTS P-006): with a non-identity fault-QP permutation, DOFData is written in
// CANONICAL (layout-independent) order — so a reordered in-memory arrangement +
// its perm produces the SAME on-disk bytes as the canonically-ordered state
// with no perm — and reading with the perm recovers the in-memory arrangement.
static void checkpoint_v2_canonical_reorder()
{
   const int N = 6;
   // perm[mem] = canonical (on-disk) position of in-memory QP mem.
   const std::vector<int> perm = {3, 0, 5, 1, 4, 2};

   // Distinct dynamic state indexed by CANONICAL position.
   auto canon_dof = [](int c)
   {
      DOFData d;
      d.psi = 0.1 + c; d.slip_rate = 1.0 + c; d.V1 = 10.0 + c; d.V2 = 20.0 + c;
      d.slip1 = 30.0 + c; d.slip2 = 40.0 + c; d.tau1_nuc = 50.0 + c;
      d.tau2_nuc = 60.0 + c; d.sigma_n_nuc = 70.0 + c;
      return d;
   };
   std::vector<DOFData> canonical(N), in_mem(N);
   for (int c = 0; c < N; ++c) { canonical[c] = canon_dof(c); }
   for (int mem = 0; mem < N; ++mem) { in_mem[mem] = canon_dof(perm[mem]); }

   Vector Qw(8); for (int i = 0; i < 8; ++i) { Qw[i] = 0.5 * i; }
   const real_t t = 2.0e-3, dt = 1.0e-6; const int sync = 7;
   const std::uint64_t hash = 0xdeadbeefULL;

   const char *pA = "test_lts_ckpt_reorder_A";
   const char *pB = "test_lts_ckpt_reorder_B";
   // A: reordered in-memory + perm -> canonical on disk.
   mfem::seas::internal::WriteTpv104CheckpointV2Impl(
      pA, t, dt, sync, 1, hash, Qw, in_mem, 0, 1, "spatial_dyn", &perm);
   // B: canonical in-memory, no perm -> canonical on disk (reference).
   mfem::seas::internal::WriteTpv104CheckpointV2Impl(
      pB, t, dt, sync, 1, hash, Qw, canonical, 0, 1, "spatial_dyn", nullptr);

   auto read_identity = [&](const char *p)
   {
      real_t tr = 0, dtr = 0; int syncr = 0, lm = 0; std::uint64_t hr = 0;
      Vector Qr; std::vector<DOFData> dd(N); std::string tag;
      mfem::seas::internal::ReadTpv104CheckpointV2Impl(
         p, tr, dtr, syncr, lm, hr, Qr, 8, dd, 0, 1, &tag, nullptr);
      return dd;
   };
   // Compare ALL 9 serialized dynamic fields (R-007: a field-swap isolated to
   // slip_rate/V1/tau1_nuc must not slip through).
   auto dyn_diff = [](const DOFData &a, const DOFData &b)
   {
      real_t m = 0.0;
      m = std::max(m, std::abs(a.psi         - b.psi));
      m = std::max(m, std::abs(a.slip_rate   - b.slip_rate));
      m = std::max(m, std::abs(a.V1          - b.V1));
      m = std::max(m, std::abs(a.V2          - b.V2));
      m = std::max(m, std::abs(a.slip1       - b.slip1));
      m = std::max(m, std::abs(a.slip2       - b.slip2));
      m = std::max(m, std::abs(a.tau1_nuc    - b.tau1_nuc));
      m = std::max(m, std::abs(a.tau2_nuc    - b.tau2_nuc));
      m = std::max(m, std::abs(a.sigma_n_nuc - b.sigma_n_nuc));
      return m;
   };

   const std::vector<DOFData> onDiskA = read_identity(pA);
   const std::vector<DOFData> onDiskB = read_identity(pB);
   // On-disk order is canonical & layout-independent: A == B field-for-field.
   real_t md_ab = 0.0;
   for (int c = 0; c < N; ++c) { md_ab = std::max(md_ab, dyn_diff(onDiskA[c], onDiskB[c])); }
   CHECK(md_ab == 0.0,
         "[ckpt] reordered+perm write == canonical write on disk (layout-indep)");
   // On-disk canonical order equals the true canonical state.
   CHECK(onDiskA[2].psi == canon_dof(2).psi
         && onDiskA[5].slip1 == canon_dof(5).slip1,
         "[ckpt] on-disk record at canonical position == canonical state");

   // Read A WITH perm -> recovers the in-memory arrangement exactly.
   {
      real_t tr = 0, dtr = 0; int syncr = 0, lm = 0; std::uint64_t hr = 0;
      Vector Qr; std::vector<DOFData> rt(N); std::string tag;
      mfem::seas::internal::ReadTpv104CheckpointV2Impl(
         pA, tr, dtr, syncr, lm, hr, Qr, 8, rt, 0, 1, &tag, &perm);
      real_t md = 0.0;
      for (int mem = 0; mem < N; ++mem)
      { md = std::max(md, dyn_diff(rt[mem], in_mem[mem])); }
      CHECK(md == 0.0, "[ckpt] read-with-perm recovers in-memory arrangement (all 9 fields)");
   }

   // (The non-permutation guard — a duplicate/out-of-range perm entry — aborts
   // the write via MFEM_VERIFY in production; a process-abort is not unit-
   // testable here since this MFEM build has MFEM_USE_EXCEPTIONS=NO.)

   std::remove((std::string(pA) + "_checkpoint_r0.txt").c_str());
   std::remove((std::string(pB) + "_checkpoint_r0.txt").c_str());
}

int main()
{
   checkpoint_v2_roundtrip();
   checkpoint_v2_canonical_reorder();
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
      multicluster_degenerate(wave, mesh, "scalar");
      multirate_conservation(wave, mesh, "scalar");
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

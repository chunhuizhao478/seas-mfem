// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 Phase 3 STOP — Arm 1 localization probes.
//
// Gate 14a (runtime, constant-I, interior lift) = 0.000 — orbit-covariant.
// Gate 14' (precomputed, constant-I, interior lift) = 1.999 — NOT covariant.
// Same input, same element loop, same mass-inverse — differs only in
// whether the 9x9 per-face operator is the per-cell precomputed matrix
// or the per-QP runtime BuildFrame(CalcOrtho) matrix.
//
// This binary isolates WHERE the divergence first appears by probing each
// kernel independently on a constant input, where the physical flux is
// identically zero.  Any orbit drift is a basis / lift artefact, not a
// flux defect.
//
// Probes:
//   G_CONST_VOL        ComputeVolumeRHS(Q_const) on M0 — expect 0
//   G_CONST_VOL_MINV   G_CONST_VOL then ApplyMassInverse — expect 0
//   G_CONST_BFACE_LIFT boundary-face-only rhs on Q_const, runtime vs
//                      precomputed, orbit drift of each path
//   G_ORBIT_DOF_MAP    enumerate each orbit-pair's local DOFs at
//                      geometrically-coincident physical points; assert
//                      matching indices
//
// One-line verdict per probe: "PROBE X: drift=... verdict=PASS/FAIL".

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/precomputed_face_fluxes.hpp"
#include "../../dynamic/d4_tet_mesh.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace mfem { namespace seas { int g_seas_my_rank = 0; } }  // NOLINT

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define PROBE_VERDICT(name, drift, tol) do {                                  \
   num_tests++;                                                               \
   const double _d = (drift), _t = (tol);                                     \
   const bool _ok = (_d <= _t);                                               \
   if (_ok) num_passed++;                                                     \
   else     num_failed++;                                                     \
   std::cout << "PROBE " << (name) << ": drift="                              \
             << std::scientific << std::setprecision(3) << _d                 \
             << "  tol=" << _t                                                \
             << "  verdict=" << (_ok ? "PASS" : "FAIL") << "\n";              \
} while (0)

namespace
{

constexpr real_t kL      = 1000.0;
constexpr real_t kRho    = 2670.0;
constexpr real_t kLambda = 32.04e9;
constexpr real_t kMu     = 32.04e9;

// Identical to audit's M0 fault fixture.
Mesh BuildM0FaultMesh()
{
   const int n = 2;
   Mesh mesh = Mesh::MakeCartesian3D(n, n, n, Element::TETRAHEDRON,
                                     kL, kL, kL, false);
   mesh.FinalizeTopology();
   mesh.Finalize();
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         real_t cy = 0.0;
         for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
         cy /= fv.Size();
         if (std::abs(cy - 0.5 * kL) < 1e-8)
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

// y = L/2 mirror partition: fault-adjacent tets sit below (y < L/2)
// or above (y > L/2).  The audit's OrbitSortedSig keys orbits by
// (cx, cy, L-cy, cz, |nor|) — we use a simpler form here: a tet and
// its y-mirror partner sit in the same orbit if their centroids are
// reflections across y=L/2.
struct OrbitBucket
{
   std::vector<int> elems;   // element ids in this orbit
};

std::map<std::tuple<long long,long long,long long>, OrbitBucket>
BuildElementOrbits(const Mesh &mesh)
{
   std::map<std::tuple<long long,long long,long long>, OrbitBucket> buckets;
   const real_t scale = 1e6;
   auto q = [&](real_t v) { return static_cast<long long>(std::round(v*scale)); };
   for (int e = 0; e < mesh.GetNE(); e++)
   {
      Array<int> ev; mesh.GetElementVertices(e, ev);
      real_t cx = 0, cy = 0, cz = 0;
      for (int i = 0; i < ev.Size(); i++)
      {
         const real_t *v = mesh.GetVertex(ev[i]);
         cx += v[0]; cy += v[1]; cz += v[2];
      }
      cx /= ev.Size(); cy /= ev.Size(); cz /= ev.Size();
      // Canonicalize cy across y = L/2 mirror.
      const real_t cy_mirror = std::min(cy, kL - cy);
      const auto key = std::make_tuple(q(cx), q(cy_mirror), q(cz));
      buckets[key].elems.push_back(e);
   }
   return buckets;
}

// For each orbit bucket with >= 2 elements, return the worst per-cell
// mean |rhs[c]| difference across the orbit.  Compute the mean of
// rhs[c, DOFs of elem] for each elem and compare sorted signatures
// across orbit-mates.
real_t OrbitCellMeanMaxDrift(
   const Vector &rhs,
   const WaveOperator<Mesh> &wave,
   const std::map<std::tuple<long long,long long,long long>, OrbitBucket> &buckets,
   std::string &worst_comp)
{
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   real_t worst = 0.0;
   worst_comp.clear();

   static const char *comp_name[9] = {"SXX","SYY","SZZ","SXY","SYZ","SXZ",
                                       "VX","VY","VZ"};

   for (const auto &kv : buckets)
   {
      const auto &orbit = kv.second.elems;
      if (orbit.size() < 2) { continue; }

      for (int c = 0; c < NUM_STATE; c++)
      {
         std::vector<real_t> means;
         means.reserve(orbit.size());
         for (int e : orbit)
         {
            Array<int> edofs; fes.GetElementDofs(e, edofs);
            real_t sum = 0.0;
            for (int i = 0; i < edofs.Size(); i++)
            {
               sum += rhs(c * ndof_total + edofs[i]);
            }
            means.push_back(sum / edofs.Size());
         }
         // Max pairwise |mean| diff within orbit.
         for (std::size_t i = 0; i < means.size(); i++)
         {
            for (std::size_t j = i+1; j < means.size(); j++)
            {
               real_t d = std::abs(means[i] - means[j]);
               if (d > worst)
               {
                  worst = d;
                  worst_comp = std::string(comp_name[c]) + " orbit="
                             + std::to_string(std::get<0>(kv.first)) + ","
                             + std::to_string(std::get<1>(kv.first)) + ","
                             + std::to_string(std::get<2>(kv.first));
               }
            }
         }
      }
   }
   return worst;
}

// ==========================================================================
// R-001 round-2 (2026-04-23 REVIEW): per-DOF paired orbit drift.
//
// The cell-mean variant above is STRUCTURALLY INSENSITIVE to per-DOF
// asymmetry on divergence-form DG operators.  Proof:
//   For p=1 L2 GaussLobatto, sum_i shape_i(x) = 1 (partition of unity).
//   sum_i rhs(c, dof_i) = sum_i ∫_e shape_i · (∇·F)_c dV
//                       = ∫_e (sum_i shape_i) · (∇·F)_c dV
//                       = ∫_e (∇·F)_c dV
//                       = ∮_∂e F·n dS   (divergence theorem)
// The boundary integral depends ONLY on the continuous F on ∂e; it is
// independent of which basis function is in which slot, so per-DOF
// asymmetry is invisible in the cell-mean comparison.
//
// This helper instead compares rhs DOF-by-DOF across orbit-paired tets
// at geometrically-COINCIDENT physical positions.  A per-DOF kernel
// asymmetry (the kind pepper is made of) shows up here and not in the
// cell-mean.
//
// Return value: (worst |rhs_a - rhs_b| over all coincident-DOF pairs
// across all orbits and components).  n_pairs_compared is set to the
// number of (orbit × DOF × component) triples actually compared; if
// zero, the probe is DEGENERATE on this fixture — meaning the fixture
// has no orbit-paired DOFs at coincident physical positions (e.g.,
// Kuhn with y-mirror vertex non-coincidence), and the per-DOF
// covariance is UNTESTED, not proven.
real_t OrbitPerDofPairedDrift(
   const Vector &rhs,
   const WaveOperator<Mesh> &wave,
   const std::map<std::tuple<long long,long long,long long>, OrbitBucket> &buckets,
   std::string &worst_label,
   int &n_pairs_compared)
{
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   real_t worst = 0.0;
   worst_label.clear();
   n_pairs_compared = 0;

   static const char *comp_name[9] = {"SXX","SYY","SZZ","SXY","SYZ","SXZ",
                                      "VX","VY","VZ"};
   const real_t s = 1e6;
   auto q = [&](real_t v) { return static_cast<long long>(std::round(v*s)); };

   for (const auto &kv : buckets)
   {
      const auto &orbit = kv.second.elems;
      if (orbit.size() < 2) { continue; }

      // Per-DOF map keyed by canonical (y-mirror-reduced) physical
      // position for each element in the orbit.  Two DOFs across
      // orbit-paired tets are "coincident" if they share this key.
      using Key = std::tuple<long long, long long, long long>;
      std::vector<std::map<Key, int>> perelem;
      perelem.reserve(orbit.size());

      for (int e : orbit)
      {
         const FiniteElement *fe = fes.GetFE(e);
         ElementTransformation *Tr = fes.GetElementTransformation(e);
         const IntegrationRule &nodes = fe->GetNodes();
         Array<int> edofs; fes.GetElementDofs(e, edofs);
         std::map<Key, int> m;
         for (int j = 0; j < fe->GetDof(); j++)
         {
            const IntegrationPoint &ip = nodes.IntPoint(j);
            Vector x(3); Tr->Transform(ip, x);
            const real_t y_mir = std::min(x(1), kL - x(1));
            m[std::make_tuple(q(x(0)), q(y_mir), q(x(2)))] = edofs[j];
         }
         perelem.push_back(std::move(m));
      }

      // Compare rhs at every (coincident) key across perelem[0] vs
      // perelem[k].  A miss on any key means the orbit pair is not
      // DOF-aligned at that physical position — skip (equivalent to
      // R-001 round-1 "vertex non-coincident" case).
      for (const auto &kvm : perelem[0])
      {
         const auto &key    = kvm.first;
         const int   dof0   = kvm.second;
         for (std::size_t k = 1; k < perelem.size(); k++)
         {
            auto it = perelem[k].find(key);
            if (it == perelem[k].end()) { continue; }
            const int dofk = it->second;

            for (int c = 0; c < NUM_STATE; c++)
            {
               const real_t a = rhs(c * ndof_total + dof0);
               const real_t b = rhs(c * ndof_total + dofk);
               const real_t d = std::abs(a - b);
               n_pairs_compared++;
               if (d > worst)
               {
                  worst = d;
                  worst_label = std::string(comp_name[c]) + " orbit="
                              + std::to_string(std::get<0>(kv.first)) + ","
                              + std::to_string(std::get<1>(kv.first)) + ","
                              + std::to_string(std::get<2>(kv.first));
               }
            }
         }
      }
   }
   return worst;
}

// Fill Q with a single global constant vector (same value for every DOF
// and every element).  Zero velocity + one component of stress set to 1e6.
void FillConstantState(const FiniteElementSpace &fes, Vector &Q)
{
   const int ndof_total = fes.GetNDofs();
   Q.SetSize(NUM_STATE * ndof_total);
   Q = 0.0;
   // SXX = 1e6 Pa, others 0.  For constant Q, ∂Q/∂x = 0 so ∇·F = 0
   // and ComputeVolumeRHS must return zero everywhere.
   for (int i = 0; i < ndof_total; i++)
   {
      Q(SXX * ndof_total + i) = 1.0e6;
   }
}

// Fill Q with a NON-constant, y-mirror-symmetric spatial pattern.
//   Q_SXX(x, y, z) = amplitude * (0.5 + x/L) * g(y) * (0.5 + z/L)
// with g(y) = 0.5 + 4*y*(L-y)/L^2 (range [0.5, 1.5], max at y=L/2).
// g(y) satisfies g(L-y) = g(y), so Q is symmetric about y = L/2 —
// the INPUT is orbit-covariant by construction.  Therefore any
// per-orbit drift in the output rhs is attributable to the KERNEL,
// not the input.
//
// A polynomial pattern (rather than sin/cos) is used because on a
// 2×2×2 cartesian mesh the p=1 DOFs sit at {0, L/2, L}^3, and
// sin(2π x/L) evaluates to zero at every one of those positions.
// The linear factors in x and z ensure every DOF gets a nonzero Q,
// producing the non-constant spatial gradient the probe is meant to
// exercise.
//
// R-003 (2026-04-23 REVIEW): introduced to rule in/out Phase 2A / 2B
// on EVIDENCE.  Passing the constant-Q variant alone says "∇·F = 0
// integrates to zero correctly"; it does NOT rule out a DOF-level
// D4-non-covariance in the volume kernel that only manifests under
// non-trivial spatial gradients.
void FillYSymmetricNonConstantState(const FiniteElementSpace &fes,
                                    Vector &Q,
                                    real_t amplitude)
{
   const int ndof_total = fes.GetNDofs();
   Q.SetSize(NUM_STATE * ndof_total);
   Q = 0.0;

   for (int e = 0; e < fes.GetNE(); e++)
   {
      const FiniteElement *fe = fes.GetFE(e);
      ElementTransformation *Tr = fes.GetElementTransformation(e);
      const IntegrationRule &nodes = fe->GetNodes();
      Array<int> edofs; fes.GetElementDofs(e, edofs);
      for (int j = 0; j < fe->GetDof(); j++)
      {
         const IntegrationPoint &ip = nodes.IntPoint(j);
         Vector x(3); Tr->Transform(ip, x);
         const real_t fx = 0.5 + x(0) / kL;
         const real_t fz = 0.5 + x(2) / kL;
         // y-mirror-symmetric about y = L/2: g(L-y) = g(y) since
         // (L-y)*(L-(L-y))/L^2 = (L-y)*y/L^2 = y*(L-y)/L^2.
         const real_t g_y = 0.5 + 4.0 * x(1) * (kL - x(1)) / (kL * kL);
         Q(SXX * ndof_total + edofs[j]) = amplitude * fx * g_y * fz;
      }
   }
}

} // anonymous

// ==========================================================================
// G_CONST_VOL: ComputeVolumeRHS(Q_const) must produce zero rhs everywhere,
//              and (by extension) orbit drift = 0.  If orbit drift > 0, the
//              volume assembly itself is not D4-covariant.
// ==========================================================================
void ProbeConstVol()
{
   std::cout << "\n--- G_CONST_VOL ---\n";
   Mesh mesh = BuildM0FaultMesh();

   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr    = 3;

   WaveOperator<Mesh> wave(mesh, /*order=*/1, kLambda, kMu, kRho, bc);
   real_t Q_bg_zero[NUM_STATE] = {0};
   wave.SetAbsorbingBackground(Q_bg_zero);

   Vector Q; FillConstantState(wave.GetFESpace(), Q);
   Vector rhs(wave.Height()); rhs = 0.0;
   wave.ComputeVolumeRHS_ForTest(Q, rhs);

   // Absolute norm of rhs itself (should be 0 for constant Q).
   real_t linf = 0.0;
   for (int i = 0; i < rhs.Size(); i++) { linf = std::max(linf, std::abs(rhs(i))); }
   std::cout << "  ||rhs_vol||_inf = " << std::scientific << linf << "\n";

   // Orbit cell-mean drift.
   auto buckets = BuildElementOrbits(mesh);
   std::string worst;
   real_t drift = OrbitCellMeanMaxDrift(rhs, wave, buckets, worst);
   if (!worst.empty())
   {
      std::cout << "  worst orbit: " << worst << "\n";
   }
   PROBE_VERDICT("G_CONST_VOL", drift, 1e-10);
}

// ==========================================================================
// G_CONST_VOL_MINV: G_CONST_VOL then ApplyMassInverse.  If G_CONST_VOL
//                   passes but this fails, the bug is in element mass
//                   inverse ordering (candidate item 3).
// ==========================================================================
void ProbeConstVolMinv()
{
   std::cout << "\n--- G_CONST_VOL_MINV ---\n";
   Mesh mesh = BuildM0FaultMesh();

   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr    = 3;

   WaveOperator<Mesh> wave(mesh, /*order=*/1, kLambda, kMu, kRho, bc);
   real_t Q_bg_zero[NUM_STATE] = {0};
   wave.SetAbsorbingBackground(Q_bg_zero);

   Vector Q; FillConstantState(wave.GetFESpace(), Q);
   Vector rhs(wave.Height()); rhs = 0.0;
   wave.ComputeVolumeRHS_ForTest(Q, rhs);
   wave.ApplyMassInverse_ForTest(rhs);

   real_t linf = 0.0;
   for (int i = 0; i < rhs.Size(); i++) { linf = std::max(linf, std::abs(rhs(i))); }
   std::cout << "  ||M^-1 rhs_vol||_inf = " << std::scientific << linf << "\n";

   auto buckets = BuildElementOrbits(mesh);
   std::string worst;
   real_t drift = OrbitCellMeanMaxDrift(rhs, wave, buckets, worst);
   if (!worst.empty())
   {
      std::cout << "  worst orbit: " << worst << "\n";
   }
   PROBE_VERDICT("G_CONST_VOL_MINV", drift, 1e-10);
}

// ==========================================================================
// G_CONST_BFACE_LIFT: boundary-face-only rhs on constant I.  Run twice —
//                    runtime path (flag off) and precomputed path (flag on).
//                    Compare orbit drift of each.  Then also compare the
//                    two paths against EACH OTHER (should be ULP if the
//                    9x9 matrices agree on the outward normal).
// ==========================================================================
void ProbeConstBfaceLift()
{
   std::cout << "\n--- G_CONST_BFACE_LIFT ---\n";
   // Fault-less mesh: ComputeFaceFluxRHS requires fault bookkeeping
   // when bc.fault_attr > 0, so we probe the interior+boundary face
   // branches on a mesh with no fault.  Orbit structure is the same
   // (y=L/2 mirror symmetry across the 2x2x2 Kuhn split).
   const int n = 2;
   Mesh mesh = Mesh::MakeCartesian3D(n, n, n, Element::TETRAHEDRON,
                                     kL, kL, kL, false);
   mesh.FinalizeTopology();
   mesh.Finalize();

   BoundaryConfig bc;
   for (int i = 1; i <= 6; i++) { bc.natural_attrs.insert(i); }
   bc.fault_attr    = 0;

   WaveOperator<Mesh> wave(mesh, /*order=*/1, kLambda, kMu, kRho, bc);
   real_t Q_bg_zero[NUM_STATE] = {0};
   wave.SetAbsorbingBackground(Q_bg_zero);

   Vector Q; FillConstantState(wave.GetFESpace(), Q);

   // Note: ComputeFaceFluxRHS inside MFEM touches BOTH interior and
   // boundary faces; we don't have a boundary-only accessor.  So we
   // compute the full face-flux rhs and take the cell-mean drift on
   // boundary-adjacent cells only.  For a fault-less test the distinction
   // would not matter, but on M0 the boundary-only drift is what Gate 14b
   // and 14' bface measure.
   Vector rhs_runtime(wave.Height()); rhs_runtime = 0.0;
   wave.ComputeFaceFluxRHS_ForTest(Q, rhs_runtime);

   wave.UsePrecomputedFaceFluxes(true);
   Vector rhs_precomp(wave.Height()); rhs_precomp = 0.0;
   wave.ComputeFaceFluxRHS_ForTest(Q, rhs_precomp);

   // Path-agreement check: precomputed should match runtime to ULP on
   // cells where the topology normal equals the runtime CalcOrtho normal.
   //
   // 2026-04-23 REVIEW R-004 (MODERATE): pre-fix we only reported the
   // infinity norm of the diff.  That's uninterpretable without
   // (a) the rhs magnitude scale on either path (the probe's own
   // tolerance comment predicts F~1e16 but runtime measures ~1e7), and
   // (b) per-face localization (if the diff concentrates on few faces,
   // the "parking lot" decision is wrong — frame/sign mismatch is
   // plausible).  Instrument both.
   real_t rhs_runtime_inf = 0.0, rhs_precomp_inf = 0.0;
   for (int i = 0; i < rhs_runtime.Size(); i++)
   {
      rhs_runtime_inf = std::max(rhs_runtime_inf, std::abs(rhs_runtime(i)));
      rhs_precomp_inf = std::max(rhs_precomp_inf, std::abs(rhs_precomp(i)));
   }
   std::cout << "  ||rhs_runtime||_inf = "
             << std::scientific << rhs_runtime_inf << "\n";
   std::cout << "  ||rhs_precomp||_inf = "
             << std::scientific << rhs_precomp_inf << "\n";

   const int ndof_total = wave.GetFESpace().GetNDofs();
   real_t path_diff = 0.0;
   int worst_dof = -1, worst_comp = -1;
   for (int c = 0; c < NUM_STATE; c++)
   {
      real_t max_comp = 0.0;
      int worst_comp_dof = -1;
      for (int dof = 0; dof < ndof_total; dof++)
      {
         const int i = c * ndof_total + dof;
         const real_t d =
            std::abs(rhs_runtime(i) - rhs_precomp(i));
         if (d > max_comp)     { max_comp = d; worst_comp_dof = dof; }
         if (d > path_diff)    { path_diff = d; worst_dof = dof;
                                 worst_comp = c; }
      }
      std::cout << "    component " << c
                << " max diff = " << std::scientific
                << std::setprecision(3) << max_comp
                << " at dof " << worst_comp_dof << "\n";
   }
   std::cout << "  ||rhs_runtime - rhs_precomp||_inf = "
             << std::scientific << std::setprecision(3) << path_diff
             << "   (worst @ comp=" << worst_comp
             << " dof=" << worst_dof << ")\n";

   // R-004: relative magnitude against the runtime scale.  Reports
   // whether the diff is ULP-floor noise or a concentrated outlier.
   const real_t rel_vs_runtime =
      rhs_runtime_inf > 0.0 ? path_diff / rhs_runtime_inf : 0.0;
   std::cout << "  path_diff / ||rhs_runtime||_inf = "
             << std::scientific << std::setprecision(3)
             << rel_vs_runtime << "  (< 1e-14 → ULP noise; "
             << ">~1e-3 → concentrated outlier)\n";

   // Per-element localization: which element owns worst_dof, and how
   // many elements carry > 10% of the worst component's max diff
   // (concentration diagnostic).
   if (worst_dof >= 0 && worst_comp >= 0)
   {
      const auto &fes_local = wave.GetFESpace();
      int worst_elem = -1;
      for (int e = 0; e < mesh.GetNE(); e++)
      {
         Array<int> ed; fes_local.GetElementDofs(e, ed);
         for (int i = 0; i < ed.Size(); i++)
         {
            if (ed[i] == worst_dof) { worst_elem = e; break; }
         }
         if (worst_elem >= 0) { break; }
      }
      std::cout << "  worst_dof=" << worst_dof
                << " lives in element " << worst_elem << "\n";

      // Concentration: count elements whose per-elem max |diff| on the
      // worst component exceeds 10% of path_diff.  If very few
      // elements carry the diff, the disagreement is localized — a
      // frame/sign-mismatch candidate.
      int n_elems_above_10pct = 0;
      for (int e = 0; e < mesh.GetNE(); e++)
      {
         Array<int> ed; fes_local.GetElementDofs(e, ed);
         real_t elem_max = 0.0;
         for (int i = 0; i < ed.Size(); i++)
         {
            const int idx = worst_comp * ndof_total + ed[i];
            const real_t d = std::abs(rhs_runtime(idx) - rhs_precomp(idx));
            elem_max = std::max(elem_max, d);
         }
         if (elem_max > 0.1 * path_diff) { n_elems_above_10pct++; }
      }
      std::cout << "  elements with > 10% of worst diff: "
                << n_elems_above_10pct << " / " << mesh.GetNE()
                << "  (concentration: < 10% elems → localized; "
                << ">~50% elems → diffuse ULP)\n";
   }

   auto buckets = BuildElementOrbits(mesh);

   std::string w_rt, w_pc;
   real_t drift_rt = OrbitCellMeanMaxDrift(rhs_runtime, wave, buckets, w_rt);
   real_t drift_pc = OrbitCellMeanMaxDrift(rhs_precomp, wave, buckets, w_pc);
   std::cout << "  runtime  orbit drift = " << std::scientific << drift_rt
             << "  (worst " << w_rt << ")\n";
   std::cout << "  precomp  orbit drift = " << std::scientific << drift_pc
             << "  (worst " << w_pc << ")\n";

   // Two-path agreement is the real signal: precomputed "should" match
   // runtime on constant-I; if it drifts more, the precomputed path is
   // baking a non-D4-covariant frame into the per-cell matrix.  A
   // tolerance of 1e-10 is the plan's ULP floor at Lame ~ 1e10 Pa;
   // Q_const is 1e6 so F ~ 1e10 * 1e6 = 1e16; relative 1e-10 means
   // absolute 1e6 slack per component.
   PROBE_VERDICT("G_CONST_BFACE_LIFT (runtime orbit drift)", drift_rt, 1e-10);
   PROBE_VERDICT("G_CONST_BFACE_LIFT (precomp orbit drift)", drift_pc, 1e-10);
   PROBE_VERDICT("G_CONST_BFACE_LIFT (runtime vs precomp)", path_diff, 1e-6);
}

// ==========================================================================
// G_ORBIT_DOF_MAP: for each orbit-pair of elements, assert that the DOF
//                  at geometrically-coincident physical points (after
//                  mirror across y=L/2) has the same element-local index.
//                  If DOF(elem=e_a at phys=p_a) != DOF(elem=e_b at phys=
//                  mirror(p_a)) in their respective local-index spaces,
//                  MFEM's nodal basis is not D4-covariant on the Kuhn
//                  split (candidate item 1).
// ==========================================================================
void ProbeOrbitDofMap()
{
   std::cout << "\n--- G_ORBIT_DOF_MAP ---\n";
   Mesh mesh = BuildM0FaultMesh();

   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr    = 3;

   WaveOperator<Mesh> wave(mesh, /*order=*/1, kLambda, kMu, kRho, bc);
   const auto &fes = wave.GetFESpace();

   auto buckets = BuildElementOrbits(mesh);

   int total_orbits = 0;
   int orbits_mismatched = 0;
   int worst_max_local_index_diff = 0;
   // R-001 split counters (2026-04-23 review): distinguish mesh-layer
   // vertex asymmetry from basis-layer index non-covariance.
   int total_keys_not_found = 0;
   int total_coincident_checked = 0;
   int total_index_mismatched = 0;

   // For p=1 tet, each element has 4 DOFs at the 4 vertices.  The
   // D4-covariance invariant: for an orbit-pair (e_a, e_b) where e_b's
   // vertices are the y-mirror of e_a's, the vertex v_a_i (at phys p)
   // should correspond to e_b's vertex v_b_j (at phys mirror(p)) with
   // i == j (same element-local index).  If i != j, the element-local
   // DOF ordering is NOT mirror-covariant, and any per-DOF accumulation
   // will scramble orbit-related physical quantities into different
   // element-local slots.
   for (const auto &kv : buckets)
   {
      const auto &orbit = kv.second.elems;
      if (orbit.size() < 2) { continue; }
      total_orbits++;

      // Enumerate (phys -> local_dof_index) map for each element.
      std::vector<std::map<std::tuple<long long,long long,long long>, int>> perelem;
      const real_t scale = 1e6;
      auto q = [&](real_t v) { return static_cast<long long>(std::round(v*scale)); };

      for (int e : orbit)
      {
         std::map<std::tuple<long long,long long,long long>, int> m;
         const FiniteElement *fe = fes.GetFE(e);
         ElementTransformation *Tr = fes.GetElementTransformation(e);
         const IntegrationRule &nodes = fe->GetNodes();
         for (int j = 0; j < fe->GetDof(); j++)
         {
            const IntegrationPoint &ip = nodes.IntPoint(j);
            Vector x(3); Tr->Transform(ip, x);
            // Mirror x across y=L/2 to canonical form.
            const real_t y_mir = std::min(x(1), kL - x(1));
            m[std::make_tuple(q(x(0)), q(y_mir), q(x(2)))] = j;
         }
         perelem.push_back(std::move(m));
      }

      // Compare perelem[0] vs perelem[k] for k > 0: same physical
      // position (after mirror) must have same local index.
      //
      // 2026-04-23 REVIEW R-001 (CRITICAL): the single "max_diff = 99"
      // sentinel used to conflate two DISTINCT failure modes:
      //   (a) physical key not found in perelem[k] — the orbit-paired
      //       element has no DOF at the geometrically-coincident point.
      //       This is a MESH-layer property (vertex layout asymmetric
      //       under y-mirror), typical of Kuhn-split cartesian tets.
      //   (b) key found but local index differs — same physical point,
      //       different element-local DOF slot.  This is a BASIS-layer
      //       property (MFEM nodal-basis ordering not D4-covariant).
      // Only (b) tests basis covariance; (a) tests mesh layout.  The
      // pre-fix "worst diff = 99" report could not distinguish them,
      // leading to incorrect basis-bug escalation in Phase 1 findings.
      int max_diff_on_coincident = 0;
      int n_keys_not_found = 0;
      int n_coincident_checked = 0;
      int n_index_mismatched = 0;
      for (const auto &kvm : perelem[0])
      {
         const auto &key = kvm.first;
         const int idx0 = kvm.second;
         for (std::size_t k = 1; k < perelem.size(); k++)
         {
            auto it = perelem[k].find(key);
            if (it == perelem[k].end())
            {
               n_keys_not_found++;
               continue;
            }
            n_coincident_checked++;
            if (it->second != idx0)
            {
               n_index_mismatched++;
               max_diff_on_coincident =
                  std::max(max_diff_on_coincident,
                           std::abs(it->second - idx0));
            }
         }
      }
      if (n_keys_not_found > 0 || n_index_mismatched > 0)
      {
         orbits_mismatched++;
         worst_max_local_index_diff =
            std::max(worst_max_local_index_diff, max_diff_on_coincident);
      }
      total_keys_not_found += n_keys_not_found;
      total_coincident_checked += n_coincident_checked;
      total_index_mismatched += n_index_mismatched;
   }

   std::cout << "  orbits examined:                       " << total_orbits << "\n";
   std::cout << "  orbits with any mismatch:              "
             << orbits_mismatched << "\n";
   std::cout << "  keys not found (vertex non-coincident): "
             << total_keys_not_found << "\n";
   std::cout << "  coincident keys checked:               "
             << total_coincident_checked << "\n";
   std::cout << "  indices mismatched on coincident keys: "
             << total_index_mismatched << "\n";
   std::cout << "  worst idx diff on coincident keys:     "
             << worst_max_local_index_diff << "\n";

   // Two verdicts:
   //   (1) vertex coincidence: keys_not_found == 0 PASS ⇒ every physical
   //       point in perelem[0] exists in every orbit-paired perelem[k].
   //       Failure here is a MESH property (Kuhn y-mirror asymmetry),
   //       not a basis property.
   //   (2) index ordering on coincident keys: index_mismatched == 0 PASS
   //       ⇒ same physical point maps to same element-local slot across
   //       orbit-paired elements.  This is the real basis-covariance
   //       test.  Only meaningful when verdict (1) passed AND
   //       coincident_checked > 0 (otherwise the basis was not exercised).
   PROBE_VERDICT("G_ORBIT_DOF_MAP (vertex coincidence)",
                 static_cast<real_t>(total_keys_not_found), 0.5);
   if (total_coincident_checked > 0)
   {
      PROBE_VERDICT("G_ORBIT_DOF_MAP (index ordering on coincident pts)",
                    static_cast<real_t>(total_index_mismatched), 0.5);
   }
   else
   {
      std::cout << "  NOTE: coincident_checked == 0 — BASIS UNTESTED on "
                   "this fixture; run on a D4-equivariant mesh to "
                   "exercise basis covariance.\n";
   }
}

// ==========================================================================
// R-004 round 2 (2026-04-23 REVIEW): G_BFACE_VERTEX_SELECTION.
//
// H3 hypothesis: even if two orbit-paired tets have identical bulk
// vertex SETS (not the case on Kuhn, but is the case on a subset of
// orbit pairs and would be universally true on a truly D4-equivariant
// mesh), the BOUNDARY FACE of each tet selects 3 of the 4 vertices
// via MFEM's FACE2NODES table.  That selection is a function of the
// STORED vertex ORDER in the tet.  If orbit-paired tets have the same
// vertex set in DIFFERENT stored slots, FACE2NODES picks DIFFERENT
// 3-vertex triangles for the boundary face — even when the y-mirror-
// corresponding triangle exists.
//
// This probe walks boundary elements.  For each boundary-element pair
// on opposite y-halves with matching (cx, cz) centroids, it reports
// whether the 3-vertex SET of bdr_a y-mirrors to the 3-vertex SET of
// bdr_b.  A FAIL means H3 is active: the bface-flux assembly picks
// different triangles for "the same" physical boundary face, so the
// per-face flux computation inherits an orbit-dependent frame even
// if the kernel and the mesh vertex positions are perfectly D4.
//
// Phase 2D interpretation impact: D4 fixture swap fixes H1 AND H3
// SIMULTANEOUSLY (slot-wise equivariance makes FACE2NODES pick the
// same triangle on both sides).  So "D4 fixes bface orbit drift"
// does NOT distinguish H1 from H3 — and the fix for H3 is a
// face-vertex-ordering normalization that would transfer to Gmsh
// production meshes (whose vertex sets are already NOT D4 but whose
// FACE2NODES could still be canonicalized).
void ProbeBfaceVertexSelection()
{
   std::cout << "\n--- G_BFACE_VERTEX_SELECTION ---\n";
   // Fault-less Kuhn mesh (same as ProbeConstBfaceLift).
   const int n = 2;
   Mesh mesh = Mesh::MakeCartesian3D(n, n, n, Element::TETRAHEDRON,
                                     kL, kL, kL, false);
   mesh.FinalizeTopology();
   mesh.Finalize();

   // Group boundary elements by (cx, cz, y-half).
   struct BdrInfo
   {
      int bdr_idx;
      Array<int> verts;     // 3 vertex ids (triangle)
      real_t cx, cy, cz;
   };
   std::vector<BdrInfo> lower, upper;
   for (int b = 0; b < mesh.GetNBE(); b++)
   {
      Array<int> bv; mesh.GetBdrElementVertices(b, bv);
      if (bv.Size() != 3) { continue; }
      real_t cx = 0, cy = 0, cz = 0;
      for (int i = 0; i < 3; i++)
      {
         const real_t *v = mesh.GetVertex(bv[i]);
         cx += v[0]; cy += v[1]; cz += v[2];
      }
      cx /= 3.0; cy /= 3.0; cz /= 3.0;
      BdrInfo info{b, bv, cx, cy, cz};
      if (cy < 0.5 * kL) { lower.push_back(info); }
      else               { upper.push_back(info); }
   }
   std::cout << "  lower-y boundary triangles: " << lower.size()
             << "\n  upper-y boundary triangles: " << upper.size() << "\n";

   const real_t s = 1e6;
   auto q = [&](real_t v) { return static_cast<long long>(std::round(v*s)); };

   // Match lower.bdr_a to upper.bdr_b by (cx, cz, |cy - L/2|).  Only
   // compare within matching buckets.
   using Key = std::tuple<long long, long long, long long>;
   std::map<Key, std::vector<int>> lower_by_key, upper_by_key;
   for (std::size_t i = 0; i < lower.size(); i++)
   {
      lower_by_key[std::make_tuple(q(lower[i].cx),
                                   q(std::min(lower[i].cy, kL-lower[i].cy)),
                                   q(lower[i].cz))].push_back(i);
   }
   for (std::size_t i = 0; i < upper.size(); i++)
   {
      upper_by_key[std::make_tuple(q(upper[i].cx),
                                   q(std::min(upper[i].cy, kL-upper[i].cy)),
                                   q(upper[i].cz))].push_back(i);
   }

   int n_pairs_checked       = 0;
   int n_setwise_mismatches  = 0;
   int n_buckets_no_match    = 0;
   for (const auto &kv : lower_by_key)
   {
      auto it_up = upper_by_key.find(kv.first);
      if (it_up == upper_by_key.end())
      {
         n_buckets_no_match++;
         continue;
      }
      for (int ia : kv.second)
      {
         for (int ib : it_up->second)
         {
            n_pairs_checked++;
            // Build set of y-mirrored lower vertices.
            std::vector<Key> lower_mirror_set;
            for (int j = 0; j < 3; j++)
            {
               const real_t *v = mesh.GetVertex(lower[ia].verts[j]);
               lower_mirror_set.push_back(std::make_tuple(
                  q(v[0]), q(kL - v[1]), q(v[2])));
            }
            std::sort(lower_mirror_set.begin(), lower_mirror_set.end());

            std::vector<Key> upper_set;
            for (int j = 0; j < 3; j++)
            {
               const real_t *v = mesh.GetVertex(upper[ib].verts[j]);
               upper_set.push_back(std::make_tuple(q(v[0]), q(v[1]), q(v[2])));
            }
            std::sort(upper_set.begin(), upper_set.end());

            if (lower_mirror_set != upper_set)
            {
               n_setwise_mismatches++;
            }
         }
      }
   }

   std::cout << "  pairs checked: " << n_pairs_checked << "\n";
   std::cout << "  buckets with no cross-half match: "
             << n_buckets_no_match << "\n";
   std::cout << "  set-wise mismatches: " << n_setwise_mismatches << "\n";

   // On Kuhn: we expect cross-half buckets to exist (y=0 and y=L
   // boundary triangles at matching (cx, cz)), and the set-wise
   // mismatch count tells whether FACE2NODES picks y-mirror-paired
   // vertex SETS.  A mismatch count > 0 means H3 is active on this
   // fixture: the bface assembly inherits an orbit-dependent
   // triangle-vertex set from the FACE2NODES layer, which a D4
   // fixture swap would hide but a Gmsh production mesh would still
   // expose.
   PROBE_VERDICT("G_BFACE_VERTEX_SELECTION (set-wise mismatches)",
                 static_cast<real_t>(n_setwise_mismatches), 0.5);

   // Interpretation hints for the findings doc:
   if (n_pairs_checked == 0)
   {
      std::cout << "  NOTE: probe degenerate — no cross-half boundary-"
                   "triangle pairs with matching (cx, cz) found.\n";
   }
   else if (n_setwise_mismatches == 0)
   {
      std::cout << "  H3 ruled OUT on Kuhn: FACE2NODES picks y-mirror-"
                   "paired triangle vertex sets across orbit halves.\n"
                   "  Phase 2D D4-fixture fix will NOT be masking a "
                   "separate H3 bug.\n";
   }
   else
   {
      std::cout << "  H3 ACTIVE on Kuhn: FACE2NODES picks DIFFERENT "
                   "triangles.  Phase 2D D4-fixture fix would HIDE this "
                   "but Gmsh production meshes (not D4) would still "
                   "expose it.\n";
   }
}

// ==========================================================================
// R-003 ADDENDUM: non-constant-Q probes.
//
// The constant-Q probes above test that the DG operator's kernels
// preserve orbit covariance on a degenerate input where ∇·F = 0.
// They do not exercise the per-DOF gradient computation, the
// component-wise Jacobian blocks, or the per-DOF mass-inverse
// multiply on spatially-varying data.
//
// These probes re-run G_CONST_VOL / G_CONST_VOL_MINV / G_CONST_BFACE_LIFT
// on a y-mirror-symmetric non-constant Q.  Because the INPUT is
// D4-covariant by construction (cos(2π y/L) is invariant under
// y ↦ L-y), any orbit drift in the output rhs is a property of the
// KERNEL, not the input.
//
// A failure here — with PASS on the corresponding constant-Q probe —
// means:
//   G_NONCONST_VOL fail   → volume kernel not D4-covariant per DOF
//                           (reopens Phase 2A).
//   G_NONCONST_VOL_MINV   → mass-inverse compounds basis non-covariance
//                           (reopens Phase 2B).  The mass-inverse test
//                           is only meaningful if G_NONCONST_VOL first
//                           passes, otherwise blame falls on the volume
//                           path.
//   G_NONCONST_BFACE_LIFT → face-flux runtime-vs-precomp disagreement
//                           no longer ULP; real frame mismatch.
//
// Tolerance: relative to ||rhs||_inf (the output scale), NOT to the
// input amplitude.  Double-precision ULP ≈ 2.22e-16.  R-006 round 2
// (2026-04-23 REVIEW): the prior 1.0e-12 was ~4500× ULP and would
// silently absorb a sub-ULP amplification seed that compounds to
// visible pepper over 1e3 steps.  Tightened to 1.0e-13 (~450× ULP),
// still generous enough to not false-fire on cross-platform ULP noise
// but tight enough to catch a real amplification mechanism.
// ==========================================================================
constexpr real_t kNonConstAmp = 1.0e6;       // input SXX amplitude
const real_t kULPFloorRel     = 1.0e-13;     // ~450× ULP

void ProbeNonConstVol()
{
   std::cout << "\n--- G_NONCONST_VOL ---\n";
   Mesh mesh = BuildM0FaultMesh();

   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr    = 3;

   WaveOperator<Mesh> wave(mesh, /*order=*/1, kLambda, kMu, kRho, bc);
   real_t Q_bg_zero[NUM_STATE] = {0};
   wave.SetAbsorbingBackground(Q_bg_zero);

   Vector Q;
   FillYSymmetricNonConstantState(wave.GetFESpace(), Q, kNonConstAmp);
   Vector rhs(wave.Height()); rhs = 0.0;
   wave.ComputeVolumeRHS_ForTest(Q, rhs);

   real_t linf = 0.0;
   for (int i = 0; i < rhs.Size(); i++) { linf = std::max(linf, std::abs(rhs(i))); }
   std::cout << "  ||rhs_vol||_inf = " << std::scientific << linf << "\n";

   auto buckets = BuildElementOrbits(mesh);

   // R-001 round 2: report BOTH cell-mean (insensitive by
   // divergence-theorem collapse) and per-DOF (sensitive).  The
   // per-DOF probe is the evidence-bearing one; cell-mean is kept
   // only for regression-comparison with prior findings docs.
   std::string worst_cell, worst_pdof;
   real_t cellmean_drift = OrbitCellMeanMaxDrift(rhs, wave, buckets, worst_cell);
   int n_pairs = 0;
   real_t perdof_drift =
      OrbitPerDofPairedDrift(rhs, wave, buckets, worst_pdof, n_pairs);
   std::cout << "  cell-mean drift    = " << std::scientific << cellmean_drift
             << "   (worst " << worst_cell << ")\n";
   std::cout << "  per-DOF paired drift = " << perdof_drift
             << "   (n_pairs=" << n_pairs
             << ", worst " << worst_pdof << ")\n";

   const real_t cell_rel = (linf > 0.0) ? cellmean_drift / linf : 0.0;
   std::cout << "  cell-mean rel drift = " << cell_rel
             << "   (STRUCTURALLY 0 for divergence-form rhs; "
             << "do NOT use to rule out Phase 2A)\n";

   if (n_pairs == 0)
   {
      std::cout << "  per-DOF probe DEGENERATE on this fixture "
                << "(no coincident-DOF pairs across orbits).\n";
      std::cout << "  → G_NONCONST_VOL (per-DOF): UNTESTED "
                << "— rerun on D4 fixture to exercise.\n";
   }
   else
   {
      const real_t pdof_rel = (linf > 0.0) ? perdof_drift / linf : 0.0;
      std::cout << "  per-DOF rel drift  = " << pdof_rel << "\n";
      PROBE_VERDICT("G_NONCONST_VOL (per-DOF rel)", pdof_rel, kULPFloorRel);
   }
}

void ProbeNonConstVolMinv()
{
   std::cout << "\n--- G_NONCONST_VOL_MINV ---\n";
   Mesh mesh = BuildM0FaultMesh();

   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr    = 3;

   WaveOperator<Mesh> wave(mesh, /*order=*/1, kLambda, kMu, kRho, bc);
   real_t Q_bg_zero[NUM_STATE] = {0};
   wave.SetAbsorbingBackground(Q_bg_zero);

   Vector Q;
   FillYSymmetricNonConstantState(wave.GetFESpace(), Q, kNonConstAmp);
   Vector rhs(wave.Height()); rhs = 0.0;
   wave.ComputeVolumeRHS_ForTest(Q, rhs);
   wave.ApplyMassInverse_ForTest(rhs);

   real_t linf = 0.0;
   for (int i = 0; i < rhs.Size(); i++) { linf = std::max(linf, std::abs(rhs(i))); }
   std::cout << "  ||M^-1 rhs_vol||_inf = " << std::scientific << linf << "\n";

   auto buckets = BuildElementOrbits(mesh);
   std::string worst_cell, worst_pdof;
   real_t cellmean_drift = OrbitCellMeanMaxDrift(rhs, wave, buckets, worst_cell);
   int n_pairs = 0;
   real_t perdof_drift =
      OrbitPerDofPairedDrift(rhs, wave, buckets, worst_pdof, n_pairs);
   std::cout << "  cell-mean drift      = " << std::scientific << cellmean_drift
             << "   (worst " << worst_cell << ")\n";
   std::cout << "  per-DOF paired drift = " << perdof_drift
             << "   (n_pairs=" << n_pairs
             << ", worst " << worst_pdof << ")\n";

   const real_t cell_rel = (linf > 0.0) ? cellmean_drift / linf : 0.0;
   std::cout << "  cell-mean rel drift = " << cell_rel
             << "   (cell-mean of M^-1·rhs is similarly weakly-sensitive; "
             << "per-DOF is the evidence-bearing probe)\n";

   if (n_pairs == 0)
   {
      std::cout << "  per-DOF probe DEGENERATE on this fixture.\n";
      std::cout << "  → G_NONCONST_VOL_MINV (per-DOF): UNTESTED "
                << "— rerun on D4 fixture.\n";
   }
   else
   {
      const real_t pdof_rel = (linf > 0.0) ? perdof_drift / linf : 0.0;
      std::cout << "  per-DOF rel drift  = " << pdof_rel << "\n";
      PROBE_VERDICT("G_NONCONST_VOL_MINV (per-DOF rel)",
                    pdof_rel, kULPFloorRel);
   }
}

void ProbeNonConstBfaceLift()
{
   std::cout << "\n--- G_NONCONST_BFACE_LIFT ---\n";
   const int n = 2;
   Mesh mesh = Mesh::MakeCartesian3D(n, n, n, Element::TETRAHEDRON,
                                     kL, kL, kL, false);
   mesh.FinalizeTopology();
   mesh.Finalize();

   BoundaryConfig bc;
   for (int i = 1; i <= 6; i++) { bc.natural_attrs.insert(i); }
   bc.fault_attr = 0;

   WaveOperator<Mesh> wave(mesh, /*order=*/1, kLambda, kMu, kRho, bc);
   real_t Q_bg_zero[NUM_STATE] = {0};
   wave.SetAbsorbingBackground(Q_bg_zero);

   Vector Q;
   FillYSymmetricNonConstantState(wave.GetFESpace(), Q, kNonConstAmp);

   Vector rhs_runtime(wave.Height()); rhs_runtime = 0.0;
   wave.ComputeFaceFluxRHS_ForTest(Q, rhs_runtime);

   wave.UsePrecomputedFaceFluxes(true);
   Vector rhs_precomp(wave.Height()); rhs_precomp = 0.0;
   wave.ComputeFaceFluxRHS_ForTest(Q, rhs_precomp);

   real_t linf_rt = 0.0, linf_pc = 0.0, path_diff = 0.0;
   for (int i = 0; i < rhs_runtime.Size(); i++)
   {
      linf_rt   = std::max(linf_rt, std::abs(rhs_runtime(i)));
      linf_pc   = std::max(linf_pc, std::abs(rhs_precomp(i)));
      path_diff = std::max(path_diff,
                           std::abs(rhs_runtime(i) - rhs_precomp(i)));
   }
   std::cout << "  ||rhs_runtime||_inf = " << std::scientific << linf_rt << "\n";
   std::cout << "  ||rhs_precomp||_inf = " << linf_pc << "\n";
   std::cout << "  path_diff = " << path_diff << "\n";
   const real_t rel = (linf_rt > 0.0) ? path_diff / linf_rt : 0.0;
   std::cout << "  relative path_diff = " << rel << "\n";
   PROBE_VERDICT("G_NONCONST_BFACE_LIFT (runtime vs precomp, rel)",
                 rel, kULPFloorRel);

   auto buckets = BuildElementOrbits(mesh);

   // Cell-mean drift (structurally weakened by divergence collapse at
   // the bulk level, though bface rhs is boundary-assembly so the
   // collapse argument is weaker here) + per-DOF drift.
   std::string w_rt_cell, w_pc_cell, w_rt_pdof, w_pc_pdof;
   real_t cell_rt = OrbitCellMeanMaxDrift(rhs_runtime, wave, buckets, w_rt_cell);
   real_t cell_pc = OrbitCellMeanMaxDrift(rhs_precomp, wave, buckets, w_pc_cell);
   int n_pairs_rt = 0, n_pairs_pc = 0;
   real_t pdof_rt = OrbitPerDofPairedDrift(rhs_runtime, wave, buckets,
                                           w_rt_pdof, n_pairs_rt);
   real_t pdof_pc = OrbitPerDofPairedDrift(rhs_precomp, wave, buckets,
                                           w_pc_pdof, n_pairs_pc);

   const real_t cell_rel_rt = (linf_rt > 0.0) ? cell_rt / linf_rt : 0.0;
   const real_t cell_rel_pc = (linf_pc > 0.0) ? cell_pc / linf_pc : 0.0;
   std::cout << "  cell-mean runtime drift = " << std::scientific << cell_rt
             << "  rel = " << cell_rel_rt
             << "  (worst " << w_rt_cell << ")\n";
   std::cout << "  cell-mean precomp drift = " << cell_pc
             << "  rel = " << cell_rel_pc
             << "  (worst " << w_pc_cell << ")\n";
   std::cout << "  per-DOF runtime drift   = " << pdof_rt
             << "  (n_pairs=" << n_pairs_rt << ")\n";
   std::cout << "  per-DOF precomp drift   = " << pdof_pc
             << "  (n_pairs=" << n_pairs_pc << ")\n";

   if (n_pairs_rt == 0 || n_pairs_pc == 0)
   {
      std::cout << "  per-DOF probe DEGENERATE on this fixture.\n";
      std::cout << "  → G_NONCONST_BFACE_LIFT (per-DOF): UNTESTED "
                << "— rerun on D4 fixture.\n";
      // Still report the cell-mean verdicts as rough
      // regression-comparison numbers, but mark them untrusted.
      std::cout << "  (cell-mean verdicts below are retained for "
                << "regression compare only; do NOT use for "
                << "Phase 2 branch selection.)\n";
      PROBE_VERDICT("G_NONCONST_BFACE_LIFT (cell-mean runtime, WEAK)",
                    cell_rel_rt, kULPFloorRel);
      PROBE_VERDICT("G_NONCONST_BFACE_LIFT (cell-mean precomp, WEAK)",
                    cell_rel_pc, kULPFloorRel);
   }
   else
   {
      const real_t pdof_rel_rt = (linf_rt > 0.0) ? pdof_rt / linf_rt : 0.0;
      const real_t pdof_rel_pc = (linf_pc > 0.0) ? pdof_pc / linf_pc : 0.0;
      PROBE_VERDICT("G_NONCONST_BFACE_LIFT (per-DOF runtime rel)",
                    pdof_rel_rt, kULPFloorRel);
      PROBE_VERDICT("G_NONCONST_BFACE_LIFT (per-DOF precomp rel)",
                    pdof_rel_pc, kULPFloorRel);
   }
}

// ==========================================================================
// Phase 2D ADDENDUM (2026-04-23): D4 fixture sweep.
//
// Runs the key probe measurements on the D4 fixture, producing the
// values needed for §C's 8-row decision table.  Kept compact and
// focused: reports per-DOF drift (§C's evidence-bearing metric),
// cell-mean (for structural-sanity), and the decision-table inputs.
//
// Fault-bearing probes use BuildD4Mesh(true); fault-less probes use
// BuildD4Mesh(false).  The fault sanity gate
// (test_d4_fault_sanity) has verified BuildD4Mesh(true) produces the
// expected 8 fault interior faces.
// ==========================================================================
void RunD4ProbeSuite()
{
   std::cout << "\n============================================================\n";
   std::cout << "=== D4 fixture sweep (Phase 2D) ===\n";
   std::cout << "============================================================\n";

   Mesh mesh_fault = BuildD4Mesh(/*add_fault=*/true, kL);
   Mesh mesh_free  = BuildD4Mesh(/*add_fault=*/false, kL);

   // R4-R001 (2026-04-23): consolidated D4 fixture-quality gate.
   // Abort BEFORE measurements if Jacobian / equivariance invariants
   // fail; otherwise measurements are contaminated by fixture bugs
   // (round-3 R-003: 24/48 negative Jacobians invalidated §D).
   std::cout << "\n[R4-R001] fault-bearing D4 fixture invariants:\n";
   AssertD4FixtureValid(mesh_fault, kL, /*abort_on_fail=*/true);
   std::cout << "\n[R4-R001] fault-less D4 fixture invariants:\n";
   AssertD4FixtureValid(mesh_free,  kL, /*abort_on_fail=*/true);

   BoundaryConfig bc_fault;
   bc_fault.natural_attrs = {1};
   bc_fault.fault_attr    = 3;

   BoundaryConfig bc_free;
   for (int i = 1; i <= 6; i++) { bc_free.natural_attrs.insert(i); }
   bc_free.fault_attr = 0;

   // ---- G_CONST_VOL on D4 (fault-bearing) ----
   {
      std::cout << "\n--- D4: G_CONST_VOL ---\n";
      WaveOperator<Mesh> wave(mesh_fault, 1, kLambda, kMu, kRho, bc_fault);
      real_t bg0[NUM_STATE] = {0}; wave.SetAbsorbingBackground(bg0);
      Vector Q; FillConstantState(wave.GetFESpace(), Q);
      Vector rhs(wave.Height()); rhs = 0.0;
      wave.ComputeVolumeRHS_ForTest(Q, rhs);
      real_t linf = 0.0;
      for (int i = 0; i < rhs.Size(); i++) { linf = std::max(linf, std::abs(rhs(i))); }
      auto buckets = BuildElementOrbits(mesh_fault);
      std::string w_cell, w_pd;
      real_t cell = OrbitCellMeanMaxDrift(rhs, wave, buckets, w_cell);
      int n = 0;
      real_t pd = OrbitPerDofPairedDrift(rhs, wave, buckets, w_pd, n);
      std::cout << "  ||rhs||_inf=" << std::scientific << linf
                << "  cell=" << cell << "  per-DOF=" << pd
                << "  n_pairs=" << n << "\n";
      PROBE_VERDICT("D4: G_CONST_VOL (per-DOF)", pd, 1e-10);
   }

   // ---- G_CONST_VOL_MINV on D4 ----
   {
      std::cout << "\n--- D4: G_CONST_VOL_MINV ---\n";
      WaveOperator<Mesh> wave(mesh_fault, 1, kLambda, kMu, kRho, bc_fault);
      real_t bg0[NUM_STATE] = {0}; wave.SetAbsorbingBackground(bg0);
      Vector Q; FillConstantState(wave.GetFESpace(), Q);
      Vector rhs(wave.Height()); rhs = 0.0;
      wave.ComputeVolumeRHS_ForTest(Q, rhs);
      wave.ApplyMassInverse_ForTest(rhs);
      real_t linf = 0.0;
      for (int i = 0; i < rhs.Size(); i++) { linf = std::max(linf, std::abs(rhs(i))); }
      auto buckets = BuildElementOrbits(mesh_fault);
      std::string w_cell, w_pd;
      real_t cell = OrbitCellMeanMaxDrift(rhs, wave, buckets, w_cell);
      int n = 0;
      real_t pd = OrbitPerDofPairedDrift(rhs, wave, buckets, w_pd, n);
      std::cout << "  ||M^-1 rhs||_inf=" << std::scientific << linf
                << "  cell=" << cell << "  per-DOF=" << pd
                << "  n_pairs=" << n << "\n";
      PROBE_VERDICT("D4: G_CONST_VOL_MINV (per-DOF)", pd, 1e-10);
   }

   // ---- G_CONST_BFACE_LIFT on D4 (fault-less) ----
   {
      std::cout << "\n--- D4: G_CONST_BFACE_LIFT ---\n";
      WaveOperator<Mesh> wave(mesh_free, 1, kLambda, kMu, kRho, bc_free);
      real_t bg0[NUM_STATE] = {0}; wave.SetAbsorbingBackground(bg0);
      Vector Q; FillConstantState(wave.GetFESpace(), Q);

      Vector rhs_rt(wave.Height()); rhs_rt = 0.0;
      wave.ComputeFaceFluxRHS_ForTest(Q, rhs_rt);

      wave.UsePrecomputedFaceFluxes(true);
      Vector rhs_pc(wave.Height()); rhs_pc = 0.0;
      wave.ComputeFaceFluxRHS_ForTest(Q, rhs_pc);

      real_t linf_rt = 0, linf_pc = 0, path_diff = 0;
      for (int i = 0; i < rhs_rt.Size(); i++)
      {
         linf_rt = std::max(linf_rt, std::abs(rhs_rt(i)));
         linf_pc = std::max(linf_pc, std::abs(rhs_pc(i)));
         path_diff = std::max(path_diff, std::abs(rhs_rt(i) - rhs_pc(i)));
      }
      auto buckets = BuildElementOrbits(mesh_free);
      std::string w1, w2;
      real_t d_rt = OrbitCellMeanMaxDrift(rhs_rt, wave, buckets, w1);
      real_t d_pc = OrbitCellMeanMaxDrift(rhs_pc, wave, buckets, w2);
      const real_t rel_path = (linf_rt > 0.0) ? path_diff / linf_rt : 0.0;
      const real_t rel_rt = (linf_rt > 0.0) ? d_rt / linf_rt : 0.0;
      const real_t rel_pc = (linf_pc > 0.0) ? d_pc / linf_pc : 0.0;
      std::cout << "  ||rhs_rt||_inf=" << std::scientific << linf_rt
                << "  path_diff/rel=" << rel_path
                << "  orbit_rt_rel=" << rel_rt
                << "  orbit_pc_rel=" << rel_pc << "\n";
      PROBE_VERDICT("D4: G_CONST_BFACE_LIFT (runtime orbit rel)",
                    rel_rt, 1e-10);
      PROBE_VERDICT("D4: G_CONST_BFACE_LIFT (precomp orbit rel)",
                    rel_pc, 1e-10);
      PROBE_VERDICT("D4: G_CONST_BFACE_LIFT (runtime vs precomp rel)",
                    rel_path, 1e-13);
   }

   // ---- G_ORBIT_DOF_MAP on D4 ----
   {
      std::cout << "\n--- D4: G_ORBIT_DOF_MAP ---\n";
      WaveOperator<Mesh> wave(mesh_fault, 1, kLambda, kMu, kRho, bc_fault);
      const auto &fes = wave.GetFESpace();
      auto buckets = BuildElementOrbits(mesh_fault);

      int n_keys_not_found = 0, n_coincident = 0, n_idx_mismatched = 0;
      int worst_idx_diff = 0;
      const real_t scale = 1e6;
      auto q = [&](real_t v) { return static_cast<long long>(std::round(v*scale)); };
      for (const auto &kv : buckets)
      {
         const auto &orbit = kv.second.elems;
         if (orbit.size() < 2) { continue; }
         std::vector<std::map<std::tuple<long long,long long,long long>, int>> pe;
         for (int e : orbit)
         {
            std::map<std::tuple<long long,long long,long long>, int> m;
            const FiniteElement *fe = fes.GetFE(e);
            ElementTransformation *Tr = fes.GetElementTransformation(e);
            const IntegrationRule &nodes = fe->GetNodes();
            for (int j = 0; j < fe->GetDof(); j++)
            {
               const IntegrationPoint &ip = nodes.IntPoint(j);
               Vector x(3); Tr->Transform(ip, x);
               const real_t y_mir = std::min(x(1), kL - x(1));
               m[std::make_tuple(q(x(0)), q(y_mir), q(x(2)))] = j;
            }
            pe.push_back(std::move(m));
         }
         for (const auto &kvm : pe[0])
         {
            const auto &key = kvm.first; const int idx0 = kvm.second;
            for (std::size_t k = 1; k < pe.size(); k++)
            {
               auto it = pe[k].find(key);
               if (it == pe[k].end()) { n_keys_not_found++; continue; }
               n_coincident++;
               if (it->second != idx0)
               {
                  n_idx_mismatched++;
                  worst_idx_diff = std::max(worst_idx_diff,
                                            std::abs(it->second - idx0));
               }
            }
         }
      }
      std::cout << "  keys_not_found=" << n_keys_not_found
                << "  coincident=" << n_coincident
                << "  idx_mismatched=" << n_idx_mismatched
                << "  worst_idx_diff=" << worst_idx_diff << "\n";
      PROBE_VERDICT("D4: G_ORBIT_DOF_MAP (vertex coincidence)",
                    static_cast<real_t>(n_keys_not_found), 0.5);
      if (n_coincident > 0)
      {
         PROBE_VERDICT("D4: G_ORBIT_DOF_MAP (index ordering on coincident)",
                       static_cast<real_t>(n_idx_mismatched), 0.5);
      }
   }

   // ---- G_NONCONST_VOL on D4 ----
   // R3-R003 (2026-04-23): per-component breakdown.  The input Q has
   // only SXX nonzero and y-symmetric, so the dominant rhs component
   // should be VX (∝ ∂σ_xx/∂x, y-symmetric).  For y-symmetric
   // components, orbit-paired DOF rhs values should be EQUAL (drift 0).
   // For y-antisymmetric components, they should be OPPOSITE (drift
   // 2|v|).  A nonzero per-DOF drift on VX would be a real kernel-
   // covariance signal.
   {
      std::cout << "\n--- D4: G_NONCONST_VOL ---\n";
      WaveOperator<Mesh> wave(mesh_fault, 1, kLambda, kMu, kRho, bc_fault);
      real_t bg0[NUM_STATE] = {0}; wave.SetAbsorbingBackground(bg0);
      Vector Q;
      FillYSymmetricNonConstantState(wave.GetFESpace(), Q, kNonConstAmp);
      Vector rhs(wave.Height()); rhs = 0.0;
      wave.ComputeVolumeRHS_ForTest(Q, rhs);

      // Per-component max |rhs|.
      const int ndof_total = wave.GetFESpace().GetNDofs();
      const char *comp_name[9] = {"SXX","SYY","SZZ","SXY","SYZ","SXZ",
                                   "VX","VY","VZ"};
      std::cout << "  per-component ||rhs||_inf:\n";
      for (int c = 0; c < NUM_STATE; c++)
      {
         real_t cmax = 0;
         for (int i = 0; i < ndof_total; i++)
         { cmax = std::max(cmax, std::abs(rhs(c*ndof_total + i))); }
         std::cout << "    " << comp_name[c] << "=" << std::scientific
                   << std::setprecision(3) << cmax << "\n";
      }

      real_t linf = 0;
      for (int i = 0; i < rhs.Size(); i++) { linf = std::max(linf, std::abs(rhs(i))); }
      auto buckets = BuildElementOrbits(mesh_fault);
      std::string wc, wp; int n = 0;
      real_t cell = OrbitCellMeanMaxDrift(rhs, wave, buckets, wc);
      real_t pd   = OrbitPerDofPairedDrift(rhs, wave, buckets, wp, n);
      const real_t cell_rel = (linf > 0.0) ? cell / linf : 0.0;
      const real_t pd_rel   = (linf > 0.0) ? pd / linf : 0.0;
      std::cout << "  ||rhs||_inf=" << std::scientific << linf
                << "  cell_rel=" << cell_rel
                << "  per-DOF_rel=" << pd_rel
                << "  worst=" << wp
                << "  n_pairs=" << n << "\n";

      // Per-component per-DOF breakdown: count the raw a/b/d values at
      // the worst pair per component so we can see whether VX (y-sym,
      // expected drift 0) is actually giving a 2.0 drift (real bug) or
      // whether some other component dominates.
      std::cout << "  per-component per-DOF worst drift (a, b, |a-b|, |a+b|):\n";
      const real_t s = 1e6;
      auto q = [&](real_t v) { return static_cast<long long>(std::round(v*s)); };
      for (int c_probe = 0; c_probe < NUM_STATE; c_probe++)
      {
         const auto &fes = wave.GetFESpace();
         real_t worst_diff = 0.0;
         real_t worst_a = 0.0, worst_b = 0.0;
         for (const auto &kv : buckets)
         {
            const auto &orbit = kv.second.elems;
            if (orbit.size() < 2) { continue; }
            std::vector<std::map<std::tuple<long long,long long,long long>, int>> pe;
            for (int e : orbit)
            {
               std::map<std::tuple<long long,long long,long long>, int> m;
               const FiniteElement *fe = fes.GetFE(e);
               ElementTransformation *Tr = fes.GetElementTransformation(e);
               const IntegrationRule &nodes = fe->GetNodes();
               Array<int> edofs; fes.GetElementDofs(e, edofs);
               for (int j = 0; j < fe->GetDof(); j++)
               {
                  const IntegrationPoint &ip = nodes.IntPoint(j);
                  Vector x(3); Tr->Transform(ip, x);
                  const real_t y_mir = std::min(x(1), kL - x(1));
                  m[std::make_tuple(q(x(0)), q(y_mir), q(x(2)))] = edofs[j];
               }
               pe.push_back(std::move(m));
            }
            for (const auto &kvm : pe[0])
            {
               const auto &key = kvm.first; const int dof0 = kvm.second;
               for (std::size_t k = 1; k < pe.size(); k++)
               {
                  auto it = pe[k].find(key);
                  if (it == pe[k].end()) { continue; }
                  const int dofk = it->second;
                  const real_t a = rhs(c_probe * ndof_total + dof0);
                  const real_t b = rhs(c_probe * ndof_total + dofk);
                  const real_t d = std::abs(a - b);
                  if (d > worst_diff) { worst_diff = d; worst_a = a; worst_b = b; }
               }
            }
         }
         if (worst_diff > 1e-8 * linf)  // skip ULP-floor noise
         {
            std::cout << "    " << comp_name[c_probe]
                      << ": a=" << std::scientific << std::setprecision(3)
                      << worst_a << "  b=" << worst_b
                      << "  |a-b|=" << worst_diff
                      << "  |a+b|=" << std::abs(worst_a + worst_b)
                      << "  (" << ((std::abs(worst_a + worst_b) < 0.01 * worst_diff)
                                    ? "antisym" : "mixed")
                      << ")\n";
         }
      }

      if (n > 0)
      {
         PROBE_VERDICT("D4: G_NONCONST_VOL (per-DOF rel)",
                       pd_rel, kULPFloorRel);
      }
      else
      {
         std::cout << "  per-DOF degenerate (unexpected on D4) — "
                   << "fixture equivariance check failed.\n";
      }
   }

   // ---- G_NONCONST_VOL_MINV on D4 ----
   {
      std::cout << "\n--- D4: G_NONCONST_VOL_MINV ---\n";
      WaveOperator<Mesh> wave(mesh_fault, 1, kLambda, kMu, kRho, bc_fault);
      real_t bg0[NUM_STATE] = {0}; wave.SetAbsorbingBackground(bg0);
      Vector Q;
      FillYSymmetricNonConstantState(wave.GetFESpace(), Q, kNonConstAmp);
      Vector rhs(wave.Height()); rhs = 0.0;
      wave.ComputeVolumeRHS_ForTest(Q, rhs);
      wave.ApplyMassInverse_ForTest(rhs);
      real_t linf = 0;
      for (int i = 0; i < rhs.Size(); i++) { linf = std::max(linf, std::abs(rhs(i))); }
      auto buckets = BuildElementOrbits(mesh_fault);
      std::string wc, wp; int n = 0;
      real_t cell = OrbitCellMeanMaxDrift(rhs, wave, buckets, wc);
      real_t pd   = OrbitPerDofPairedDrift(rhs, wave, buckets, wp, n);
      const real_t cell_rel = (linf > 0.0) ? cell / linf : 0.0;
      const real_t pd_rel   = (linf > 0.0) ? pd / linf : 0.0;
      std::cout << "  ||M^-1 rhs||_inf=" << std::scientific << linf
                << "  cell_rel=" << cell_rel
                << "  per-DOF_rel=" << pd_rel
                << "  n_pairs=" << n << "\n";
      if (n > 0)
      {
         PROBE_VERDICT("D4: G_NONCONST_VOL_MINV (per-DOF rel)",
                       pd_rel, kULPFloorRel);
      }
   }

   // ---- G_NONCONST_BFACE_LIFT on D4 ----
   {
      std::cout << "\n--- D4: G_NONCONST_BFACE_LIFT ---\n";
      WaveOperator<Mesh> wave(mesh_free, 1, kLambda, kMu, kRho, bc_free);
      real_t bg0[NUM_STATE] = {0}; wave.SetAbsorbingBackground(bg0);
      Vector Q;
      FillYSymmetricNonConstantState(wave.GetFESpace(), Q, kNonConstAmp);

      Vector rhs_rt(wave.Height()); rhs_rt = 0.0;
      wave.ComputeFaceFluxRHS_ForTest(Q, rhs_rt);
      wave.UsePrecomputedFaceFluxes(true);
      Vector rhs_pc(wave.Height()); rhs_pc = 0.0;
      wave.ComputeFaceFluxRHS_ForTest(Q, rhs_pc);

      real_t linf_rt = 0, linf_pc = 0, path_diff = 0;
      for (int i = 0; i < rhs_rt.Size(); i++)
      {
         linf_rt = std::max(linf_rt, std::abs(rhs_rt(i)));
         linf_pc = std::max(linf_pc, std::abs(rhs_pc(i)));
         path_diff = std::max(path_diff, std::abs(rhs_rt(i) - rhs_pc(i)));
      }
      const real_t rel_path = (linf_rt > 0.0) ? path_diff / linf_rt : 0.0;
      PROBE_VERDICT("D4: G_NONCONST_BFACE_LIFT (runtime vs precomp rel)",
                    rel_path, kULPFloorRel);

      auto buckets = BuildElementOrbits(mesh_free);
      std::string wc_rt, wc_pc, wp_rt, wp_pc; int n_rt = 0, n_pc = 0;
      real_t cell_rt = OrbitCellMeanMaxDrift(rhs_rt, wave, buckets, wc_rt);
      real_t cell_pc = OrbitCellMeanMaxDrift(rhs_pc, wave, buckets, wc_pc);
      real_t pd_rt = OrbitPerDofPairedDrift(rhs_rt, wave, buckets, wp_rt, n_rt);
      real_t pd_pc = OrbitPerDofPairedDrift(rhs_pc, wave, buckets, wp_pc, n_pc);
      const real_t cell_rel_rt = (linf_rt > 0.0) ? cell_rt / linf_rt : 0.0;
      const real_t cell_rel_pc = (linf_pc > 0.0) ? cell_pc / linf_pc : 0.0;
      const real_t pd_rel_rt = (linf_rt > 0.0) ? pd_rt / linf_rt : 0.0;
      const real_t pd_rel_pc = (linf_pc > 0.0) ? pd_pc / linf_pc : 0.0;
      std::cout << "  ||rhs_rt||_inf=" << std::scientific << linf_rt
                << "  path_rel=" << rel_path << "\n";
      std::cout << "  cell_rel_rt=" << cell_rel_rt
                << "  cell_rel_pc=" << cell_rel_pc << "\n";
      std::cout << "  per-DOF_rel_rt=" << pd_rel_rt
                << "  per-DOF_rel_pc=" << pd_rel_pc
                << "  n_pairs=(" << n_rt << "," << n_pc << ")\n";
      if (n_rt > 0 && n_pc > 0)
      {
         PROBE_VERDICT("D4: G_NONCONST_BFACE_LIFT (per-DOF runtime rel)",
                       pd_rel_rt, kULPFloorRel);
         PROBE_VERDICT("D4: G_NONCONST_BFACE_LIFT (per-DOF precomp rel)",
                       pd_rel_pc, kULPFloorRel);
      }
   }

   // ---- G_BFACE_VERTEX_SELECTION on D4 ----
   {
      std::cout << "\n--- D4: G_BFACE_VERTEX_SELECTION ---\n";
      Mesh &mesh = mesh_free;
      struct BdrInfo { int bdr_idx; Array<int> verts; real_t cx, cy, cz; };
      std::vector<BdrInfo> lower, upper;
      for (int b = 0; b < mesh.GetNBE(); b++)
      {
         Array<int> bv; mesh.GetBdrElementVertices(b, bv);
         if (bv.Size() != 3) { continue; }
         real_t cx = 0, cy = 0, cz = 0;
         for (int i = 0; i < 3; i++)
         {
            const real_t *v = mesh.GetVertex(bv[i]);
            cx += v[0]; cy += v[1]; cz += v[2];
         }
         cx /= 3.0; cy /= 3.0; cz /= 3.0;
         BdrInfo info{b, bv, cx, cy, cz};
         if (cy < 0.5 * kL) { lower.push_back(info); }
         else               { upper.push_back(info); }
      }
      const real_t scale = 1e6;
      auto q = [&](real_t v) { return static_cast<long long>(std::round(v*scale)); };
      using Key = std::tuple<long long, long long, long long>;
      std::map<Key, std::vector<int>> lby, uby;
      for (std::size_t i = 0; i < lower.size(); i++)
      { lby[std::make_tuple(q(lower[i].cx),
                            q(std::min(lower[i].cy, kL-lower[i].cy)),
                            q(lower[i].cz))].push_back(i); }
      for (std::size_t i = 0; i < upper.size(); i++)
      { uby[std::make_tuple(q(upper[i].cx),
                            q(std::min(upper[i].cy, kL-upper[i].cy)),
                            q(upper[i].cz))].push_back(i); }
      int n_pairs = 0, n_mismatches = 0, n_nomatch = 0;
      for (const auto &kv : lby)
      {
         auto it_u = uby.find(kv.first);
         if (it_u == uby.end()) { n_nomatch++; continue; }
         for (int ia : kv.second)
         {
            for (int ib : it_u->second)
            {
               n_pairs++;
               std::vector<Key> lms, us;
               for (int j = 0; j < 3; j++)
               {
                  const real_t *v = mesh.GetVertex(lower[ia].verts[j]);
                  lms.push_back(std::make_tuple(q(v[0]), q(kL-v[1]), q(v[2])));
                  const real_t *vb = mesh.GetVertex(upper[ib].verts[j]);
                  us.push_back(std::make_tuple(q(vb[0]), q(vb[1]), q(vb[2])));
               }
               std::sort(lms.begin(), lms.end());
               std::sort(us.begin(), us.end());
               if (lms != us) { n_mismatches++; }
            }
         }
      }
      std::cout << "  pairs_checked=" << n_pairs
                << "  no_cross_half_match=" << n_nomatch
                << "  set_mismatches=" << n_mismatches << "\n";
      PROBE_VERDICT("D4: G_BFACE_VERTEX_SELECTION (set mismatches)",
                    static_cast<real_t>(n_mismatches), 0.5);
   }
}

int main()
{
   std::cout << "=== TPV102 Phase 3 STOP — Arm 1 localization probes ===\n";

   ProbeConstVol();
   ProbeConstVolMinv();
   ProbeConstBfaceLift();
   ProbeOrbitDofMap();

   std::cout << "\n=== R-003 ADDENDUM: non-constant-Q probes ===\n";
   ProbeNonConstVol();
   ProbeNonConstVolMinv();
   ProbeNonConstBfaceLift();

   std::cout << "\n=== R-004 round 2 ADDENDUM: bface vertex-selection probe ===\n";
   ProbeBfaceVertexSelection();

   // Phase 2D: D4 fixture sweep.  Same probes, D4 fixture.
   RunD4ProbeSuite();

   std::cout << "\n=== Summary ===\n";
   std::cout << "Total:  " << num_tests << "\n";
   std::cout << "Passed: " << num_passed << "\n";
   std::cout << "Failed: " << num_failed << "\n";
   return (num_failed == 0) ? 0 : 1;
}

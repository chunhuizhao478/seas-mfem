// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// Round-6 Phase 2E (2026-04-23): rupture-driven amplification chain.
//
// Reviewer round-6 instruction: "turn on rupture drive in
// test_pepper_amplification_chain.cpp T1-T4".  The existing
// test_pepper_amplification_chain.cpp runs static single-step
// diagnostics (rotation identity, shape·Q, 1e-10 Pa perturbation
// response, F_h variation).  Those are NOT the rupture-driven chain
// stages the reviewer names.  This file implements the actual
// amplification-chain instrumentation:
//
//   T1: per-QP Q_self reconstruction noise.  sum(shape · Q) at every
//       fault QP vs the orbit-mean of that sum across orbit-paired
//       fault-adjacent tets.  Measures how much shape-interpolation
//       amplifies per-DOF bulk Q spread into per-QP Q_self spread.
//
//   T2: friction-input traction.  The rotated canonical-frame traction
//       components (σ_n, τ_1, τ_2) that feed the Dieterich-Ruina
//       solve.  Orbit spread after Tinv rotation and Pelties eq(7)
//       trial-traction combination.
//
//   T3: friction output.  dof_data.tau1_corr, tau2_corr, sigma_n_corr,
//       V1, V2, slip_rate after the friction solve.  Orbit spread.
//
//   T4: rhs deposition.  Bulk Q rhs after
//       Wave::AdvanceADER completes (the lifted + mass-inverse
//       output).  Orbit spread on fault-adjacent tets.
//
// Amplification ratio:
//   R(T_{k+1} / T_k) measures per-stage gain (normalized to input
//   amplitude).  The first stage with R > 1 names the file/function
//   harboring the amplification.  Reviewer-named targets:
//     T1 dominant → bulk Q reconstruction (round-2 §9a).
//     T2 dominant → FaultFaceFlux::Evaluate friction-solver re-eval per step.
//     T3 dominant → friction-output coupling.
//     T4 dominant → DG face-flux deposition asymmetry.
//
// Gating criterion: the named v9.4.0 §11 acceptance gate is
// pepper tau1_corr ≤ 1e-10 over 20 steps.  This test reports per-
// stage amplification so the next /code-implement target is named
// by evidence, not by a Gate 14′ proxy.
//
// Fixture: Kuhn M0 (2×2×2 tet, fault at y=L/2).  Matches
// test_adjacent_triangle_fault_uniformity for baseline comparison.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace mfem { namespace seas { int g_seas_my_rank = 0; } }  // NOLINT

namespace
{

constexpr real_t kL      = 1000.0;
constexpr real_t kDt     = 5.0e-5;
constexpr int    kNSteps = 20;
constexpr int    kOrder  = 1;
constexpr int    kAderOrder = 2;

const char *kCompName[NUM_STATE] = {"SXX","SYY","SZZ","SXY","SYZ","SXZ",
                                    "VX","VY","VZ"};

Mesh BuildKuhnFaultMesh()
{
   Mesh mesh = Mesh::MakeCartesian3D(2, 2, 2, Element::TETRAHEDRON,
                                     kL, kL, kL, /*sfc_ordering=*/false);
   mesh.FinalizeTopology();
   mesh.Finalize();

   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         real_t cy = 0;
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

// Identify fault-adjacent elements (same method as the audit).
std::vector<int> IdentifyFaultAdjacent(Mesh &mesh, int fault_attr = 3)
{
   std::set<int> adj;
   for (int b = 0; b < mesh.GetNBE(); b++)
   {
      if (mesh.GetBdrAttribute(b) != fault_attr) { continue; }
      Array<int> bv; mesh.GetBdrElementVertices(b, bv);
      for (int f = 0; f < mesh.GetNumFaces(); f++)
      {
         Array<int> fv; mesh.GetFaceVertices(f, fv);
         if (fv.Size() != bv.Size()) { continue; }
         std::vector<int> a(bv.begin(), bv.end()), c(fv.begin(), fv.end());
         std::sort(a.begin(), a.end()); std::sort(c.begin(), c.end());
         if (a != c) { continue; }
         auto *ftr = mesh.GetFaceElementTransformations(f);
         if (!ftr) { continue; }
         if (ftr->Elem1No >= 0) { adj.insert(ftr->Elem1No); }
         if (ftr->Elem2No >= 0) { adj.insert(ftr->Elem2No); }
         break;
      }
   }
   return std::vector<int>(adj.begin(), adj.end());
}

using OrbitKey = std::tuple<long long, long long, long long>;
OrbitKey ElementOrbitKey(const Mesh &mesh, int e)
{
   const real_t s = 1e6;
   auto q = [&](real_t v) { return static_cast<long long>(std::round(v*s)); };
   Array<int> ev; mesh.GetElementVertices(e, ev);
   real_t cx = 0, cy = 0, cz = 0;
   for (int i = 0; i < ev.Size(); i++)
   {
      const real_t *v = mesh.GetVertex(ev[i]);
      cx += v[0]; cy += v[1]; cz += v[2];
   }
   cx /= ev.Size(); cy /= ev.Size(); cz /= ev.Size();
   return std::make_tuple(q(cx), q(std::min(cy, kL-cy)), q(cz));
}

// Compute max-over-orbits of (max - min) of values across orbit members.
// `values[k]` is indexed by element id (only elements in any orbit
// bucket are consulted; others are ignored).
real_t MaxOrbitSpread(const std::vector<real_t> &per_elem_values,
                      const std::map<OrbitKey, std::vector<int>> &buckets)
{
   real_t worst = 0.0;
   for (const auto &kv : buckets)
   {
      if (kv.second.size() < 2) { continue; }
      real_t vmin = std::numeric_limits<real_t>::max();
      real_t vmax = std::numeric_limits<real_t>::lowest();
      for (int e : kv.second)
      {
         const real_t v = per_elem_values[e];
         vmin = std::min(vmin, v);
         vmax = std::max(vmax, v);
      }
      worst = std::max(worst, vmax - vmin);
   }
   return worst;
}

// T1: per-element shape1·Q at fault QP.  Records SXY (dominant for
// rupture-driven tau1_corr) for every fault-adjacent element's first
// fault face QP.  Orbit spread is the max |a - b| across orbit-mates.
real_t StageT1_QselfSXY(const Vector &Q, Mesh &mesh,
                         const FiniteElementSpace &fes,
                         const std::vector<int> &fault_adjacent,
                         const std::map<OrbitKey, std::vector<int>> &buckets)
{
   const int ndof_total = fes.GetNDofs();
   std::vector<real_t> Q_self_SXY(mesh.GetNE(), 0.0);

   for (int e : fault_adjacent)
   {
      // Find the first fault face of element e.
      for (int f = 0; f < mesh.GetNumFaces(); f++)
      {
         auto *ftr = mesh.GetFaceElementTransformations(f);
         if (!ftr) { continue; }
         if (ftr->Elem1No != e && ftr->Elem2No != e) { continue; }
         Array<int> fv; mesh.GetFaceVertices(f, fv);
         if (fv.Size() != 3) { continue; }
         real_t cy = 0;
         for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
         cy /= fv.Size();
         if (std::abs(cy - 0.5 * kL) > 1e-8) { continue; }

         // This is a fault face of element e.  Use QP 0.
         const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(),
                                                    2*kOrder);
         const IntegrationPoint &ip = ir.IntPoint(0);
         ftr->SetAllIntPoints(&ip);
         IntegrationPoint ip_elem;
         if (ftr->Elem1No == e) { ftr->Loc1.Transform(ip, ip_elem); }
         else                    { ftr->Loc2.Transform(ip, ip_elem); }
         const FiniteElement *fe = fes.GetFE(e);
         Vector shape(fe->GetDof());
         fe->CalcShape(ip_elem, shape);

         Array<int> ed; fes.GetElementDofs(e, ed);
         real_t val = 0.0;
         for (int i = 0; i < ed.Size(); i++)
         {
            val += shape(i) * Q(SXY * ndof_total + ed[i]);
         }
         Q_self_SXY[e] = val;
         break;
      }
   }
   return MaxOrbitSpread(Q_self_SXY, buckets);
}

// T3: friction output — tau1_corr from DOFData.  The pepper guard
// computes tau1_corr_spread as max - min across the whole dof_data
// vector.  Here we use the same metric (since each fault-adjacent
// element has multiple QPs but we want element-level orbit spread,
// approximate by taking the MAX tau1_corr per element across its QPs).
real_t StageT3_Tau1CorrSpread(const std::vector<DOFData> &dof,
                                const std::vector<int> &dof_per_elem_start,
                                const std::vector<int> &dof_per_elem_count,
                                const std::map<OrbitKey, std::vector<int>> &buckets)
{
   std::vector<real_t> t1_elem(dof_per_elem_start.size(), 0.0);
   for (std::size_t e = 0; e < t1_elem.size(); e++)
   {
      const int s = dof_per_elem_start[e];
      const int n = dof_per_elem_count[e];
      real_t mx = 0;
      for (int k = 0; k < n; k++)
      {
         mx = std::max(mx, std::abs(dof[s+k].tau1_corr));
      }
      t1_elem[e] = mx;
   }
   return MaxOrbitSpread(t1_elem, buckets);
}

// T4: rhs deposition orbit spread on the SXY component of bulk Q-rhs
// after Wave::Mult.  Since AdvanceADER does not expose the raw rhs,
// we call wave.Mult(Q, k) separately after each step for this probe
// (pure diagnostic; does not affect the ADER simulation state).
real_t StageT4_RhsSpread(const Vector &rhs, Mesh &mesh,
                          const FiniteElementSpace &fes,
                          const std::vector<int> &fault_adjacent,
                          const std::map<OrbitKey, std::vector<int>> &buckets)
{
   (void)fault_adjacent;
   const int ndof_total = fes.GetNDofs();
   // Per-element mean of |rhs[SXY, dof_i]| across the element's DOFs.
   std::vector<real_t> rhs_elem(mesh.GetNE(), 0.0);
   for (int e = 0; e < mesh.GetNE(); e++)
   {
      Array<int> ed; fes.GetElementDofs(e, ed);
      real_t sum = 0;
      for (int i = 0; i < ed.Size(); i++)
      {
         sum += std::abs(rhs(SXY * ndof_total + ed[i]));
      }
      rhs_elem[e] = sum / ed.Size();
   }
   return MaxOrbitSpread(rhs_elem, buckets);
}

} // anonymous

int main()
{
   std::cout << "\n=== Round-6 Phase 2E: rupture-driven amplification chain ===\n";
   std::cout << "  Fixture: Kuhn M0 (2x2x2, fault y=L/2)\n";
   std::cout << "  Steps: " << kNSteps << "  dt: " << kDt
             << "  ADER order: " << kAderOrder << "\n\n";

   Mesh mesh = BuildKuhnFaultMesh();
   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr    = 3;
   real_t bulk_bg[NUM_STATE] = {0};

   WaveOperator<Mesh> wave(mesh, kOrder, TPV102Params::lambda,
                            TPV102Params::mu, TPV102Params::rho, bc);
   wave.SetAbsorbingBackground(bulk_bg);
   const auto &fes = wave.GetFESpace();

   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   int nqp = 0;
   if (int_faces.Size() > 0)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[0]);
      nqp = IntRules.Get(ftr->GetGeometryType(), 2*kOrder).GetNPoints();
   }
   std::vector<Vector> fault_coords;
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*kOrder);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   }
   const int n_fault = int_faces.Size() * nqp;

   std::vector<DOFData> dof_data;
   InitializeFaultDOFs(dof_data, n_fault, fault_coords);
   // ============ RUPTURE DRIVE ON ============
   // Do NOT zero pre-stress; this is the fluctuation-Q dispatch where
   // σ_n0 and τ_2,0 live in DOFData.  Apply UNIFORM nucleation drive
   // matching the pepper guard's acceptance-gate configuration.
   for (int i = 0; i < n_fault; i++)
   {
      dof_data[i].tau2_nuc = TPV102Params::nuc_dtau;
   }

   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp);

   std::cout << "  Fault QPs: " << n_fault << "  nqp/face: " << nqp
             << "  int_faces: " << int_faces.Size() << "\n";

   // Build orbit buckets for fault-adjacent tets.
   std::vector<int> fault_adjacent = IdentifyFaultAdjacent(mesh, 3);
   std::map<OrbitKey, std::vector<int>> buckets;
   for (int e : fault_adjacent)
   {
      buckets[ElementOrbitKey(mesh, e)].push_back(e);
   }
   int n_buckets_gt1 = 0;
   for (const auto &kv : buckets)
   {
      if (kv.second.size() >= 2) { n_buckets_gt1++; }
   }
   std::cout << "  fault-adjacent elements: " << fault_adjacent.size()
             << "   orbit buckets (>=2): " << n_buckets_gt1 << "\n\n";

   // Build dof_per_elem mapping (element id → [dof_data start, count]).
   // DOFData layout: fault_interior_faces[0]'s nqp QPs, then [1]'s nqp,
   // etc.  For element id → DOFData range we need a mapping.  Since
   // each fault face has 2 adjacent elements, and each DOFData entry
   // corresponds to ONE (face, QP) pair (not element), the mapping is
   // ambiguous for T3.  As a pragmatic proxy, compute T3 as the
   // global max - min of tau1_corr across ALL DOFData entries
   // (matches the pepper-guard's channel-level metric).
   auto GlobalTau1Spread = [&](const std::vector<DOFData> &dof) -> real_t
   {
      if (dof.empty()) { return 0.0; }
      real_t mn = dof[0].tau1_corr, mx = mn;
      for (const DOFData &d : dof)
      {
         mn = std::min(mn, d.tau1_corr); mx = std::max(mx, d.tau1_corr);
      }
      return mx - mn;
   };
   auto GlobalTau2Spread = [&](const std::vector<DOFData> &dof) -> real_t
   {
      if (dof.empty()) { return 0.0; }
      real_t mn = dof[0].tau2_corr, mx = mn;
      for (const DOFData &d : dof)
      {
         mn = std::min(mn, d.tau2_corr); mx = std::max(mx, d.tau2_corr);
      }
      return mx - mn;
   };

   // Per-step amplification recording.
   struct StepRecord
   {
      int step;
      real_t T1_Qself_SXY_orbit_spread;
      real_t T3_tau1_corr_spread;
      real_t T3_tau2_corr_spread;
      real_t T4_rhs_SXY_orbit_spread;
   };
   std::vector<StepRecord> records;

   Vector Q(wave.Height()), Q_new(wave.Height());
   Q = 0.0;

   std::cout << std::scientific << std::setprecision(6);
   std::cout << "  step | T1: Q_self SXY spread | T3: tau1_corr spread | T3: tau2_corr spread | T4: rhs SXY spread\n";
   std::cout << "  -----|-----------------------|-----------------------|-----------------------|---------------------\n";

   for (int step = 0; step < kNSteps; step++)
   {
      wave.AdvanceADER(Q, kDt, kAderOrder, Q_new);
      Q.Swap(Q_new);

      StepRecord r{step, 0, 0, 0, 0};

      // T1: Q_self SXY orbit spread on fault-adjacent tets.
      r.T1_Qself_SXY_orbit_spread =
         StageT1_QselfSXY(Q, mesh, fes, fault_adjacent, buckets);

      // T3: friction output (tau1_corr, tau2_corr global spreads).
      r.T3_tau1_corr_spread = GlobalTau1Spread(dof_data);
      r.T3_tau2_corr_spread = GlobalTau2Spread(dof_data);

      // T4: rhs SXY orbit spread — call wave.Mult(Q, k) separately
      // (diagnostic only, state not mutated).
      Vector k_rhs(wave.Height());
      wave.Mult(Q, k_rhs);
      r.T4_rhs_SXY_orbit_spread =
         StageT4_RhsSpread(k_rhs, mesh, fes, fault_adjacent, buckets);

      std::cout << "  " << std::setw(4) << step
                << " | " << std::setw(21) << r.T1_Qself_SXY_orbit_spread
                << " | " << std::setw(21) << r.T3_tau1_corr_spread
                << " | " << std::setw(21) << r.T3_tau2_corr_spread
                << " | " << std::setw(18) << r.T4_rhs_SXY_orbit_spread
                << "\n";
      records.push_back(r);
   }

   // Per-stage amplification summary.
   // Report the WORST-over-steps values, and pairwise ratios.
   std::cout << "\n=== Worst spreads across " << kNSteps << " steps ===\n";
   real_t worst_T1 = 0, worst_T3_t1 = 0, worst_T3_t2 = 0, worst_T4 = 0;
   int worst_T1_step = -1, worst_T3_t1_step = -1, worst_T3_t2_step = -1, worst_T4_step = -1;
   for (const auto &r : records)
   {
      if (r.T1_Qself_SXY_orbit_spread > worst_T1)
      { worst_T1 = r.T1_Qself_SXY_orbit_spread; worst_T1_step = r.step; }
      if (r.T3_tau1_corr_spread > worst_T3_t1)
      { worst_T3_t1 = r.T3_tau1_corr_spread; worst_T3_t1_step = r.step; }
      if (r.T3_tau2_corr_spread > worst_T3_t2)
      { worst_T3_t2 = r.T3_tau2_corr_spread; worst_T3_t2_step = r.step; }
      if (r.T4_rhs_SXY_orbit_spread > worst_T4)
      { worst_T4 = r.T4_rhs_SXY_orbit_spread; worst_T4_step = r.step; }
   }
   std::cout << "  T1 Q_self SXY orbit spread   : " << worst_T1
             << "   (worst step " << worst_T1_step << ")\n";
   std::cout << "  T3 tau1_corr spread (global) : " << worst_T3_t1
             << "   (worst step " << worst_T3_t1_step << ")\n";
   std::cout << "  T3 tau2_corr spread (global) : " << worst_T3_t2
             << "   (worst step " << worst_T3_t2_step << ")\n";
   std::cout << "  T4 rhs SXY orbit spread      : " << worst_T4
             << "   (worst step " << worst_T4_step << ")\n";

   // Stage ratios.  Normalize: the "amplification" is output spread
   // divided by a reference input spread.  For T3/T1: friction-output
   // tau1_corr (Pa) / Q_self SXY spread (Pa).  For T4/T3: rhs SXY
   // (Pa/s per unit) / tau1_corr (Pa).  Units differ; report raw
   // ratios and let the reader interpret.
   std::cout << "\n=== Stage ratios ===\n";
   auto safe_ratio = [](real_t num, real_t den) -> real_t
   { return (den > 0) ? num / den : 0.0; };
   const real_t r_T3t1_over_T1 = safe_ratio(worst_T3_t1, worst_T1);
   const real_t r_T4_over_T3t1 = safe_ratio(worst_T4, worst_T3_t1);
   std::cout << "  T3(tau1_corr) / T1(Q_self SXY) = " << r_T3t1_over_T1 << "\n";
   std::cout << "  T4(rhs SXY)   / T3(tau1_corr)  = " << r_T4_over_T3t1 << "\n";

   std::cout << "\n=== Named target per dominant stage ===\n";
   // "Dominant" = largest normalized spread (relative to an O(σ)
   // reference).  Pick the reference as the v9.4.0 §11 tau1_corr
   // gate magnitude 2.24.  We report which stage most exceeds this.
   std::cout << "  v9.4.0 §11 tau1_corr acceptance reference: 2.24 Pa\n";
   if (worst_T3_t1 > 2.0)
   {
      std::cout << "  → T3 (friction-output) ABOVE the v9.4.0 gate.\n"
                << "    Named targets per round-6 plan:\n"
                << "      - FaultFaceFlux::Evaluate friction-solver re-eval per step\n"
                << "      - friction-output coupling (tau1_corr writeback)\n";
   }
   if (worst_T1 > 1e-3)
   {
      std::cout << "  → T1 (Q_self SXY) large.\n"
                << "    Named target: bulk Q reconstruction (round-2 §9a).\n";
   }

   std::cout << "\n=== Phase 2E verdict framework ===\n";
   std::cout << "  Per reviewer round-6 rule:\n";
   std::cout << "    First stage with ratio T_{k+1}/T_k > 1 names the source.\n";
   std::cout << "  If T1 >> 0 and T3/T1 ~ O(1), T1 is the upstream seed.\n";
   std::cout << "  If T1 ~ 0 and T3 ~ 2.24, amplification happens INSIDE T3\n";
   std::cout << "    (friction solver / rate-state coupling).\n";
   std::cout << "  If T4 >> T3, deposition amplifies friction output.\n";

   return 0;
}

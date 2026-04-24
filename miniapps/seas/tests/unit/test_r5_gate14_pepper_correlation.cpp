// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// Round-5 review Option 1 (2026-04-23): Gate-14′-to-pepper correlation test.
//
// Purpose
// -------
// The §C decision-table row PASS/PASS/FAIL points at "Phase 2C-fault"
// (a kernel-layer fix to PrecomputedFaceFluxes::BuildInteriorMatrices)
// as the next implement target.  That recommendation assumes Gate 14′
// residual → tau1_corr pepper.  The round-5 review observed Gate 14′
// does NOT correlate with pepper across the three fixtures measured
// (Kuhn, broken D4, fixed D4): Gate 14′ iface varies 1.77-2.00 while
// pepper tau1_corr varies 1.61-4.80.  If Gate 14′ is not a pepper
// driver, Phase 2C-fault is unjustified and would be a sixth
// iteration of the v8→v9.4 "fix this gate, named gate stays unmoved"
// pattern.
//
// This test runs two 20-step pepper-guard simulations on Kuhn M0
// (the clean baseline, matches v9.4.0 §11 acceptance gate):
//   (A) Baseline: standard AdvanceADER.
//   (B) Projected: at each step, orbit-project Q on fault-adjacent
//       tets onto the orbit-uniform subspace (zero all per-DOF drift
//       across same-orbit fault-adjacent tets).
//
// The projection is an UPPER BOUND on what any Gate 14′-targeted fix
// (including Phase 2C-fault) could achieve — because it zeros ALL
// orbit drift on fault-adjacent tets, not just the Gate 14′
// residual.  If even the upper-bound intervention fails to close
// pepper, no Phase 2C-fault variant can.
//
// Decision rule
// -------------
//   (B) tau1_corr spread ≤ 0.1 × (A)  → Gate 14′ drives pepper;
//                                       Phase 2C-fault JUSTIFIED
//   (B) tau1_corr spread ≥ 0.5 × (A)  → Gate 14′ NOT a pepper driver;
//                                       Phase 2C-fault NOT justified
//   0.1 × (A) < (B) < 0.5 × (A)       → partial signal; inconclusive
//
// Fixture
// -------
// Kuhn M0 (2×2×2 Cartesian tet, fault at y=L/2).  This is
// intentional: (1) Kuhn uses MakeCartesian3D which is guaranteed
// positive-orientation (no R-003 fixture contamination), (2) Kuhn
// tau1_corr = 2.236 matches the v9.4.0 §11 published acceptance gate
// value (2.24) to 0.2%, so results are directly comparable.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
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

// Fault-adjacent tet identification (same logic as the audit).
std::vector<int> IdentifyFaultAdjacent(const Mesh &mesh,
                                        int fault_attr = 3)
{
   std::set<int> fault_elems_set;
   for (int b = 0; b < mesh.GetNBE(); b++)
   {
      if (mesh.GetBdrAttribute(b) != fault_attr) { continue; }
      int e = -1;
      // bdr element e1 for this triangle:
      Array<int> bv; mesh.GetBdrElementVertices(b, bv);
      // For interior fault boundary triangles, the parent face's
      // adjacent elements are what we want.  Find matching interior
      // face by vertex set.
      for (int f = 0; f < mesh.GetNumFaces(); f++)
      {
         Array<int> fv; const_cast<Mesh&>(mesh).GetFaceVertices(f, fv);
         if (fv.Size() != bv.Size()) { continue; }
         std::vector<int> a(bv.begin(), bv.end()), c(fv.begin(), fv.end());
         std::sort(a.begin(), a.end()); std::sort(c.begin(), c.end());
         if (a != c) { continue; }
         auto *ftr = const_cast<Mesh&>(mesh).GetFaceElementTransformations(f);
         if (!ftr) { continue; }
         if (ftr->Elem1No >= 0) { fault_elems_set.insert(ftr->Elem1No); }
         if (ftr->Elem2No >= 0) { fault_elems_set.insert(ftr->Elem2No); }
         break;
         (void)e;
      }
   }
   return std::vector<int>(fault_elems_set.begin(), fault_elems_set.end());
}

// Orbit bucket key: (quantized cx, quantized min(cy, L-cy), quantized cz).
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

// Project Q onto orbit-uniform subspace on the specified elements.
// For each orbit bucket (containing 2+ elements), for each component
// and per-slot local-DOF index, compute the mean across the bucket
// and write it back to every element in the bucket.  This zeros the
// within-orbit per-DOF drift — the signal Gate 14′ measures.
//
// For Kuhn M0 with slot-equivariant vertex ordering (which holds on
// MakeCartesian3D before any CheckElementOrientation swap), slot-wise
// pairing is the right average.  If orbit-paired tets have permuted
// slots, this projection mixes them; that's acceptable for the
// correlation experiment (it still zeros the MEASURED drift).
void OrbitProjectOnElements(Vector &Q,
                             const FiniteElementSpace &fes,
                             const std::vector<int> &target_elems,
                             const std::map<OrbitKey, std::vector<int>> &buckets)
{
   const int ndof_total = fes.GetNDofs();
   std::set<int> target_set(target_elems.begin(), target_elems.end());

   for (const auto &kv : buckets)
   {
      // Restrict bucket to target elements.
      std::vector<int> bucket;
      for (int e : kv.second)
      {
         if (target_set.count(e)) { bucket.push_back(e); }
      }
      if (bucket.size() < 2) { continue; }

      // For each component, for each local-slot j, compute mean
      // across bucket and write back.
      const FiniteElement *fe = fes.GetFE(bucket[0]);
      const int nloc = fe->GetDof();
      for (int c = 0; c < NUM_STATE; c++)
      {
         for (int j = 0; j < nloc; j++)
         {
            real_t sum = 0.0;
            for (int e : bucket)
            {
               Array<int> ed; fes.GetElementDofs(e, ed);
               sum += Q(c * ndof_total + ed[j]);
            }
            const real_t mean = sum / bucket.size();
            for (int e : bucket)
            {
               Array<int> ed; fes.GetElementDofs(e, ed);
               Q(c * ndof_total + ed[j]) = mean;
            }
         }
      }
   }
}

struct PepperStats
{
   real_t slip_rate_spread = 0;
   real_t tau1_corr_spread = 0;
   real_t tau2_corr_spread = 0;
   real_t sigma_n_corr_spread = 0;
};

PepperStats CollectPepper(const std::vector<DOFData> &dof)
{
   PepperStats s;
   if (dof.empty()) { return s; }
   real_t sr_min = dof[0].slip_rate, sr_max = sr_min;
   real_t t1_min = dof[0].tau1_corr, t1_max = t1_min;
   real_t t2_min = dof[0].tau2_corr, t2_max = t2_min;
   real_t sn_min = dof[0].sigma_n_corr, sn_max = sn_min;
   for (const DOFData &d : dof)
   {
      sr_min = std::min(sr_min, d.slip_rate); sr_max = std::max(sr_max, d.slip_rate);
      t1_min = std::min(t1_min, d.tau1_corr); t1_max = std::max(t1_max, d.tau1_corr);
      t2_min = std::min(t2_min, d.tau2_corr); t2_max = std::max(t2_max, d.tau2_corr);
      sn_min = std::min(sn_min, d.sigma_n_corr); sn_max = std::max(sn_max, d.sigma_n_corr);
   }
   s.slip_rate_spread    = sr_max - sr_min;
   s.tau1_corr_spread    = t1_max - t1_min;
   s.tau2_corr_spread    = t2_max - t2_min;
   s.sigma_n_corr_spread = sn_max - sn_min;
   return s;
}

// Run the pepper guard 20 ADER steps; optionally project Q onto
// orbit-uniform on fault-adjacent tets after each step.
// Returns the WORST spread over all steps per channel.
PepperStats RunPepperGuard(const char *label,
                            bool project_every_step,
                            const std::vector<int> &fault_adjacent,
                            const std::map<OrbitKey, std::vector<int>> &buckets)
{
   std::cout << "\n--- Running: " << label
             << "  (project=" << (project_every_step ? "yes" : "no") << ") ---\n";

   Mesh mesh = BuildKuhnFaultMesh();
   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr    = 3;
   real_t bulk_bg[NUM_STATE] = {0};

   WaveOperator<Mesh> wave(mesh, kOrder, TPV102Params::lambda,
                            TPV102Params::mu, TPV102Params::rho, bc);
   wave.SetAbsorbingBackground(bulk_bg);

   // Build fault DOFData (same pattern as the pepper guard).
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
   for (int i = 0; i < n_fault; i++)
   {
      dof_data[i].tau2_nuc = TPV102Params::nuc_dtau;
   }

   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp);

   Vector Q(wave.Height()), Q_new(wave.Height());
   Q = 0.0;

   PepperStats worst;
   for (int step = 0; step < kNSteps; step++)
   {
      wave.AdvanceADER(Q, kDt, kAderOrder, Q_new);
      Q.Swap(Q_new);

      if (project_every_step)
      {
         OrbitProjectOnElements(Q, wave.GetFESpace(),
                                 fault_adjacent, buckets);
      }

      PepperStats s = CollectPepper(dof_data);
      worst.slip_rate_spread = std::max(worst.slip_rate_spread, s.slip_rate_spread);
      worst.tau1_corr_spread = std::max(worst.tau1_corr_spread, s.tau1_corr_spread);
      worst.tau2_corr_spread = std::max(worst.tau2_corr_spread, s.tau2_corr_spread);
      worst.sigma_n_corr_spread = std::max(worst.sigma_n_corr_spread,
                                             s.sigma_n_corr_spread);
   }

   std::cout << "  worst spreads over " << kNSteps << " steps:\n";
   std::cout << "    slip_rate    : " << std::scientific
             << std::setprecision(6) << worst.slip_rate_spread << "\n";
   std::cout << "    tau1_corr    : " << worst.tau1_corr_spread << "\n";
   std::cout << "    tau2_corr    : " << worst.tau2_corr_spread << "\n";
   std::cout << "    sigma_n_corr : " << worst.sigma_n_corr_spread << "\n";
   return worst;
}

} // anonymous

int main()
{
   std::cout << "=== Round-5 Option 1: Gate-14′-to-pepper correlation test ===\n";
   std::cout << "  Fixture: Kuhn M0 (2x2x2 tet, fault at y=L/2)\n";
   std::cout << "  Baseline tau1_corr ≈ 2.24 (v9.4.0 §11 acceptance gate)\n";
   std::cout << "  Test: does orbit-projecting Q on fault-adjacent tets "
             << "close the pepper?\n";

   Mesh mesh_for_pairing = BuildKuhnFaultMesh();
   const std::vector<int> fault_adjacent =
      IdentifyFaultAdjacent(mesh_for_pairing, /*fault_attr=*/3);

   // Build orbit buckets for fault-adjacent elements only.
   std::map<OrbitKey, std::vector<int>> buckets;
   for (int e : fault_adjacent)
   {
      buckets[ElementOrbitKey(mesh_for_pairing, e)].push_back(e);
   }
   std::cout << "  fault-adjacent elements: " << fault_adjacent.size() << "\n";
   int n_orbit_buckets_gt1 = 0;
   for (const auto &kv : buckets)
   {
      if (kv.second.size() >= 2) { n_orbit_buckets_gt1++; }
   }
   std::cout << "  orbit buckets (size >= 2): " << n_orbit_buckets_gt1 << "\n";

   // (A) BASELINE: standard pepper guard.
   PepperStats baseline =
      RunPepperGuard("A. BASELINE (Kuhn, no projection)",
                      /*project_every_step=*/false,
                      fault_adjacent, buckets);

   // (B) PROJECTED: same guard, orbit-project Q on fault-adjacent
   // tets at every step.  This eliminates ALL per-DOF orbit drift on
   // those tets, an upper bound on what any Gate 14′-targeted fix
   // could achieve.
   PepperStats projected =
      RunPepperGuard("B. PROJECTED (Kuhn, orbit-uniform on fault-adjacent)",
                      /*project_every_step=*/true,
                      fault_adjacent, buckets);

   std::cout << "\n=== Correlation verdict (tau1_corr) ===\n";
   std::cout << "  baseline  tau1_corr spread: "
             << std::scientific << std::setprecision(6)
             << baseline.tau1_corr_spread << " Pa\n";
   std::cout << "  projected tau1_corr spread: "
             << projected.tau1_corr_spread << " Pa\n";

   const real_t ratio =
      (baseline.tau1_corr_spread > 0)
      ? projected.tau1_corr_spread / baseline.tau1_corr_spread
      : 0.0;
   std::cout << "  projected / baseline       : " << std::fixed
             << std::setprecision(4) << ratio << "\n\n";

   std::cout << "  decision rule:\n";
   std::cout << "    ratio ≤ 0.1  → Phase 2C-fault JUSTIFIED\n";
   std::cout << "    ratio ≥ 0.5  → Phase 2C-fault NOT justified\n";
   std::cout << "    0.1 < ratio < 0.5 → inconclusive\n\n";

   std::string verdict;
   int exit_code;
   if (ratio <= 0.1)
   {
      verdict = "Phase 2C-fault JUSTIFIED "
                "(orbit projection eliminates >= 90% of pepper)";
      exit_code = 0;
   }
   else if (ratio >= 0.5)
   {
      verdict = "Phase 2C-fault NOT justified "
                "(orbit projection leaves >= 50% of pepper untouched; "
                "the Gate 14′ residual is NOT a primary pepper driver)";
      exit_code = 1;
   }
   else
   {
      verdict = "INCONCLUSIVE "
                "(partial reduction; Gate 14′ may contribute but other "
                "sources dominate)";
      exit_code = 2;
   }
   std::cout << "  VERDICT: " << verdict << "\n";

   // Also report the other channels for completeness.
   std::cout << "\n=== Other channels (informational) ===\n";
   auto rel = [&](real_t a, real_t b)
   { return (a > 0) ? b / a : 0.0; };
   std::cout << "  slip_rate    baseline=" << std::scientific
             << std::setprecision(3) << baseline.slip_rate_spread
             << "  projected=" << projected.slip_rate_spread
             << "  ratio=" << std::fixed << std::setprecision(4)
             << rel(baseline.slip_rate_spread, projected.slip_rate_spread)
             << "\n";
   std::cout << "  tau2_corr    baseline=" << std::scientific
             << std::setprecision(3) << baseline.tau2_corr_spread
             << "  projected=" << projected.tau2_corr_spread
             << "  ratio=" << std::fixed << std::setprecision(4)
             << rel(baseline.tau2_corr_spread, projected.tau2_corr_spread)
             << "\n";
   std::cout << "  sigma_n_corr baseline=" << std::scientific
             << std::setprecision(3) << baseline.sigma_n_corr_spread
             << "  projected=" << projected.sigma_n_corr_spread
             << "  ratio=" << std::fixed << std::setprecision(4)
             << rel(baseline.sigma_n_corr_spread, projected.sigma_n_corr_spread)
             << "\n";

   return exit_code;
}

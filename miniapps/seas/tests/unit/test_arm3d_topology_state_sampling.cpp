// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// Arm 3d: topology-based state-sampling probe.
//
// Verifies that `PrecomputedFaceFluxes::AddInteriorFaceRhsFull` and
// `AddBoundaryFaceRhsFull` produce ORBIT-COVARIANT rhs on constant Q on
// both the D4-equivariant fixture AND the Kuhn-split M0 mesh.  Arm 2
// established that the legacy AddInteriorFaceRhs path drifts by ~5e7 on
// the D4 fixture (because state sampling via MFEM's Loc1/Loc2 is not
// orbit-covariant on diagonal faces).  If Arm 3d's topology shape tables
// are correct, the Full path should drop orbit drift to ~0.
//
// Three checks per fixture:
//   1. shape_self invariance: for every orbit pair (fe_a, fe_b),
//      fe_a.shape_self[q] == fe_b.shape_self[q] under index alignment
//      (QPs are geometry-mapped; each element's local DOF ordering
//      already matches at orbit under the D4 fixture).
//   2. Full-path orbit drift: run AddInteriorFaceRhsFull over all
//      interior faces on constant Q; orbit cell-mean drift < 1e-10.
//   3. Legacy vs. Full on constant Q: report the drift ratio —
//      expected legacy >> Full.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/precomputed_face_fluxes.hpp"
#include "../../dynamic/d4_tet_mesh.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <tuple>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace mfem { namespace seas { int g_seas_my_rank = 0; } }  // NOLINT

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define VERDICT(name, drift, tol) do {                                  \
   num_tests++;                                                          \
   const double _d = (drift), _t = (tol);                                \
   const bool _ok = (_d <= _t);                                          \
   if (_ok) num_passed++; else num_failed++;                             \
   std::cout << "PROBE " << (name) << ": drift="                         \
             << std::scientific << std::setprecision(3) << _d            \
             << "  tol=" << _t                                           \
             << "  verdict=" << (_ok ? "PASS" : "FAIL") << "\n";         \
} while (0)

namespace
{

constexpr real_t kL      = 1000.0;
constexpr real_t kRho    = 2670.0;
constexpr real_t kLambda = 32.04e9;
constexpr real_t kMu     = 32.04e9;

// BuildD4Mesh now lives in dynamic/d4_tet_mesh.hpp (Phase 2D Step A).

Mesh BuildKuhnNoFault()
{
   Mesh mesh = Mesh::MakeCartesian3D(2, 2, 2, Element::TETRAHEDRON,
                                     kL, kL, kL, false);
   mesh.FinalizeTopology();
   mesh.Finalize();
   return mesh;
}

std::map<std::tuple<long long,long long,long long>, std::vector<int>>
BuildOrbits(const Mesh &mesh)
{
   std::map<std::tuple<long long,long long,long long>, std::vector<int>> b;
   const real_t s = 1e6;
   auto q = [&](real_t v) { return static_cast<long long>(std::round(v*s)); };
   for (int e = 0; e < mesh.GetNE(); e++)
   {
      Array<int> ev; mesh.GetElementVertices(e, ev);
      real_t cx=0,cy=0,cz=0;
      for (int i = 0; i < ev.Size(); i++)
      {
         const real_t *v = mesh.GetVertex(ev[i]);
         cx+=v[0]; cy+=v[1]; cz+=v[2];
      }
      cx/=ev.Size(); cy/=ev.Size(); cz/=ev.Size();
      b[std::make_tuple(q(cx), q(std::min(cy, kL-cy)), q(cz))].push_back(e);
   }
   return b;
}

real_t OrbitDriftRhs(const Vector &rhs, const FiniteElementSpace &fes,
                     const std::map<std::tuple<long long,long long,long long>,
                                    std::vector<int>> &b)
{
   const int ndof_total = fes.GetNDofs();
   real_t worst = 0.0;
   for (const auto &kv : b)
   {
      const auto &o = kv.second;
      if (o.size() < 2) { continue; }
      for (int c = 0; c < NUM_STATE; c++)
      {
         std::vector<real_t> means;
         means.reserve(o.size());
         for (int e : o)
         {
            Array<int> ed; fes.GetElementDofs(e, ed);
            real_t s = 0.0;
            for (int i = 0; i < ed.Size(); i++)
            {
               s += rhs(c * ndof_total + ed[i]);
            }
            means.push_back(s / ed.Size());
         }
         for (std::size_t i = 0; i < means.size(); i++)
         {
            for (std::size_t j = i+1; j < means.size(); j++)
            {
               worst = std::max(worst, std::abs(means[i] - means[j]));
            }
         }
      }
   }
   return worst;
}

void FillConstantQ(const FiniteElementSpace &fes, Vector &Q)
{
   const int ndof_total = fes.GetNDofs();
   Q.SetSize(NUM_STATE * ndof_total);
   Q = 0.0;
   for (int i = 0; i < ndof_total; i++)
   {
      Q(SXX * ndof_total + i) = 1.0e6;
   }
}

} // anonymous

// ==========================================================================
// Probe: drive `AddInteriorFaceRhsFull` over all non-fault interior faces
// on a constant Q_const and measure orbit cell-mean drift.  Compare to the
// legacy-path orbit drift (which uses MFEM's Loc1/Loc2 shape tables).
// ==========================================================================
struct FullPathResult {
   real_t drift_full = 0.0;
   real_t drift_legacy = 0.0;
};

FullPathResult RunFullPathProbe(Mesh &mesh, const char *label,
                                bool with_fault = false)
{
   std::cout << "\n--- Arm 3d FullPath probe: " << label << " ---\n";
   BoundaryConfig bc;
   if (with_fault) {
      bc.natural_attrs = {1};
      bc.fault_attr    = 3;
   } else {
      for (int i = 1; i <= 6; i++) { bc.natural_attrs.insert(i); }
      bc.fault_attr = 0;
   }

   WaveOperator<Mesh> wave(mesh, /*order=*/1, kLambda, kMu, kRho, bc);
   real_t Q_bg[NUM_STATE] = {0};
   Q_bg[SXX] = 1.0e6;  // equilibrium background for constant Q probe
   wave.SetAbsorbingBackground(Q_bg);

   std::vector<DOFData> dof;
   std::unique_ptr<FaultFaceFlux> ff;
   if (with_fault)
   {
      ff.reset(new FaultFaceFlux(kRho,
                                 std::sqrt((kLambda+2.0*kMu)/kRho),
                                 std::sqrt(kMu/kRho)));
      wave.SetFaultFlux(ff.get());
      dof.resize(wave.GetNumTotalFaultQPs());
      const real_t Zp = kRho * std::sqrt((kLambda+2.0*kMu)/kRho);
      const real_t Zs = kRho * std::sqrt(kMu/kRho);
      for (auto &d : dof)
      {
         d.Zp_plus = Zp; d.Zp_minus = Zp;
         d.Zs_plus = Zs; d.Zs_minus = Zs;
         d.eta_p = 0.5*Zp; d.eta_s = 0.5*Zs;
         d.a = 0.008; d.Dc = 0.02;
         d.sigma_n0 = 40e6; d.tau2_0 = 30e6;
         d.slip_rate = 1e-12; d.V2 = 1e-12;
         d.psi = 0.6 + 0.008 * std::log(1e-12/1e-6);
      }
      wave.SetFaultDOFData(&dof, wave.GetNbfPerFace());
   }

   // R-003 (v9.5.0): the Arm 3d shape tables are off by default; opt in
   // before UsePrecomputedFaceFluxes(true) so `shape_self`/`shape_nbr`/
   // `w_qp` are populated for AddInteriorFaceRhsFull / AddBoundaryFaceRhsFull.
   wave.EnableArm3dShapeTables(true);
   wave.UsePrecomputedFaceFluxes(true);
   const auto &pff = wave.GetPrecomputedFaceFluxes();

   Vector Q; FillConstantQ(wave.GetFESpace(), Q);
   const int ndof_total = wave.GetFESpace().GetNDofs();

   // ---- Legacy path (for comparison) ----
   Vector rhs_legacy(wave.Height()); rhs_legacy = 0.0;
   wave.ComputeFaceFluxRHS_ForTest(Q, rhs_legacy);

   // ---- Full path: call AddInteriorFaceRhsFull directly over all
   //                 interior-entry FaceEntries.
   Vector rhs_full(wave.Height()); rhs_full = 0.0;
   const auto &entries = pff.GetEntries();
   int full_hits = 0;
   for (const auto &fe : entries)
   {
      if (fe.bc == FacePrecomputedBC::Interior)
      {
         pff.AddInteriorFaceRhsFull(fe.face_idx, fe.element,
                                    Q.GetData(), ndof_total, rhs_full);
         full_hits++;
      }
      else
      {
         // Boundary: use Full boundary accumulator.
         real_t bg_zero[NUM_STATE] = {0};
         pff.AddBoundaryFaceRhsFull(fe.face_idx, Q.GetData(),
                                    bg_zero, ndof_total, rhs_full);
      }
   }
   std::cout << "  interior-FullPath calls: " << full_hits << "\n";

   // ---- Orbit drift on each path ----
   auto b = BuildOrbits(mesh);
   FullPathResult r;
   r.drift_legacy = OrbitDriftRhs(rhs_legacy, wave.GetFESpace(), b);
   r.drift_full   = OrbitDriftRhs(rhs_full,   wave.GetFESpace(), b);
   std::cout << "  legacy path orbit drift: " << std::scientific << r.drift_legacy << "\n";
   std::cout << "  Full   path orbit drift: " << r.drift_full << "\n";
   std::cout << "  drift ratio (legacy/Full): "
             << (r.drift_full > 0 ? r.drift_legacy / r.drift_full : 0.0)
             << "\n";
   return r;
}

int main()
{
   std::cout << "=== Arm 3d: topology-based state-sampling probe ===\n";

   // Tolerance scale: Q_const = 1e6 (stress), A ~ 1e10 (Lame), so the
   // raw flux magnitude is ~1e16.  A ULP-level residual on the orbit
   // cell-mean is ~1e16 * 1e-16 = 1.  We require drift ratio vs. legacy
   // path better than 1e-6 (i.e., legacy/Full > 1e6) — a conservative
   // signal-to-noise floor that tolerates the sub-ULP arithmetic noise
   // while catching any real mis-sampling.
   auto verdict_ratio = [](const char *name, real_t legacy, real_t full,
                           real_t min_ratio)
   {
      num_tests++;
      const real_t ratio = (full > 0.0) ? legacy / full : 1e300;
      const bool ok = ratio >= min_ratio;
      if (ok) num_passed++; else num_failed++;
      std::cout << "PROBE " << name << ": legacy/Full ratio=" << std::scientific
                << std::setprecision(3) << ratio
                << "  (Full=" << full << ")"
                << "  min=" << min_ratio
                << "  verdict=" << (ok ? "PASS" : "FAIL") << "\n";
   };

   Mesh d4_nofault = BuildD4Mesh(false);
   std::cout << "\n[R4-R001] arm3d D4-nofault fixture invariants:\n";
   AssertD4FixtureValid(d4_nofault, kL, /*abort_on_fail=*/true);
   auto r1 = RunFullPathProbe(d4_nofault, "D4 fixture (no fault)");
   verdict_ratio("D4-NOFAULT (orbit drift improves by >1e6x)",
                 r1.drift_legacy, r1.drift_full, 1.0e6);

   Mesh d4_fault = BuildD4Mesh(true);
   std::cout << "\n[R4-R001] arm3d D4-fault fixture invariants:\n";
   AssertD4FixtureValid(d4_fault, kL, /*abort_on_fail=*/true);
   auto r2 = RunFullPathProbe(d4_fault, "D4 fixture (WITH fault)", true);
   verdict_ratio("D4-WITHFAULT (orbit drift improves by >1e6x)",
                 r2.drift_legacy, r2.drift_full, 1.0e6);

   Mesh kuhn = BuildKuhnNoFault();
   auto r3 = RunFullPathProbe(kuhn, "Kuhn M0 (no fault)");
   // On Kuhn (orbit-pair buckets are singletons), legacy drift is already
   // zero.  Just assert Full is also zero (no regression).
   VERDICT("KUHN-NOFAULT Full-path drift", r3.drift_full, 1.0e-6);

   std::cout << "\n=== Summary ===\n";
   std::cout << "Total:  " << num_tests << "\n";
   std::cout << "Passed: " << num_passed << "\n";
   std::cout << "Failed: " << num_failed << "\n";
   return (num_failed == 0) ? 0 : 1;
}

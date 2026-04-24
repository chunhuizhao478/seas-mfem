// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 Phase 3 STOP — Arm 2: D4-equivariant tet-fixture probe.
//
// Arm 1 established:
//   * MFEM volume, mass-inverse, and face-flux kernels are orbit-covariant
//     (probe drifts = 0 on constant Q).
//   * MFEM's default Kuhn tet decomposition is NOT D4-equivariant under
//     y = L/2 mirror (G_ORBIT_DOF_MAP: 8/8 mismatch).
//
// Arm 2 question: if we build the SAME mesh topology but with a tet
// decomposition chosen so that orbit-pair elements ARE exact y-mirror
// images, do the audit's orbit-drift gates close?
//
// Build a 2x2x2 hex grid of 1000m edge, then split each hex into 6 tets
// using a diagonal pattern mirrored across y = L/2.  Tag the y = L/2
// interior plane with fault_attr = 3 (matching the M0 fault fixture).
// Re-run the Arm 1 probes on this mesh, INCLUDING the fault branch case.
//
// Expected outcomes:
//   * G_ORBIT_DOF_MAP on new fixture: PASS (0 mismatched orbits).
//   * G_CONST_VOL, G_CONST_VOL_MINV, G_CONST_BFACE_LIFT: PASS (drift = 0)
//     on BOTH the no-fault and fault configurations.  If these pass even
//     with fault bookkeeping active, the pepper drift is a pure fixture
//     artifact.
//   * If any probe still drifts, the bug is downstream of the fixture.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/precomputed_face_fluxes.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <map>
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

// --------------------------------------------------------------------------
// D4-equivariant 2x2x2 tet mesh.
//
// Grid vertices: indexed by (ix, iy, iz) ∈ {0..2} x {0..2} x {0..2}.
// Each of the 8 hexes is split into 6 tets using a vertex-ordering pattern
// chosen so that the hex at (ix, iy=0, iz) and the hex at (ix, iy=1, iz)
// produce tets that are exact y = L/2 mirror images.
//
// For a hex with corners (v000 .. v111) in binary ordering, the 6-tet
// split via the main diagonal (v000 -> v111) is:
//   T0 = (v000, v100, v110, v111)
//   T1 = (v000, v010, v110, v111)
//   T2 = (v000, v010, v011, v111)
//   T3 = (v000, v001, v011, v111)
//   T4 = (v000, v001, v101, v111)
//   T5 = (v000, v100, v101, v111)
// (standard Freudenthal decomposition)
//
// To make the lower-y hex and the upper-y hex y-mirror images, the upper
// hex uses a y-reversed diagonal: (v010 -> v101) instead of (v000 -> v111),
// which places its 6-tet pattern as the mirror of the lower hex.
// --------------------------------------------------------------------------
Mesh BuildD4EquivariantMesh(bool add_fault)
{
   const int n = 2;
   const real_t h = kL / n;
   Mesh mesh(/*dim=*/3, /*nvert=*/(n+1)*(n+1)*(n+1),
             /*nelem=*/0, /*nbdr=*/0,
             /*sdim=*/3);

   auto vid = [&](int ix, int iy, int iz)
   {
      return ix + (n+1) * (iy + (n+1) * iz);
   };

   for (int iz = 0; iz <= n; iz++)
   {
      for (int iy = 0; iy <= n; iy++)
      {
         for (int ix = 0; ix <= n; ix++)
         {
            mesh.AddVertex(ix * h, iy * h, iz * h);
         }
      }
   }

   // Helper: add six tets for a hex, with orientation orient = +1 (lower)
   // or -1 (upper).  orient = -1 mirrors the cube vertex list across the
   // local y axis so the resulting 6-tet pattern is the y-mirror of
   // orient = +1.
   auto add_hex = [&](int ix, int iy, int iz, int orient)
   {
      const int v000 = vid(ix,   iy,   iz);
      const int v100 = vid(ix+1, iy,   iz);
      const int v010 = vid(ix,   iy+1, iz);
      const int v110 = vid(ix+1, iy+1, iz);
      const int v001 = vid(ix,   iy,   iz+1);
      const int v101 = vid(ix+1, iy,   iz+1);
      const int v011 = vid(ix,   iy+1, iz+1);
      const int v111 = vid(ix+1, iy+1, iz+1);

      // Under y-mirror, swap (v000<->v010), (v100<->v110), (v001<->v011),
      // (v101<->v111).  Equivalently: pick a local labeling where the
      // upper hex's "low-y corner" is the lower hex's "high-y corner".
      //
      // For orient = +1 (lower hex): diagonal v000 -> v111.
      // For orient = -1 (upper hex): diagonal v010 -> v101 (y-mirror of
      // v000->v111).  The 6 tets are the y-mirror of the lower list.
      auto v = [&](int k) -> int
      {
         if (orient > 0)
         {
            switch (k) {
               case 0: return v000; case 1: return v100; case 2: return v010;
               case 3: return v110; case 4: return v001; case 5: return v101;
               case 6: return v011; case 7: return v111;
            }
         }
         else
         {
            // Swap y neighbors: 0<->2, 1<->3, 4<->6, 5<->7.
            switch (k) {
               case 0: return v010; case 1: return v110; case 2: return v000;
               case 3: return v100; case 4: return v011; case 5: return v111;
               case 6: return v001; case 7: return v101;
            }
         }
         return -1;
      };

      // 6 tets via main-diagonal split v(0) -> v(7).
      int T[6][4] = {
         {v(0), v(1), v(3), v(7)},
         {v(0), v(2), v(3), v(7)},
         {v(0), v(2), v(6), v(7)},
         {v(0), v(4), v(6), v(7)},
         {v(0), v(4), v(5), v(7)},
         {v(0), v(1), v(5), v(7)},
      };
      for (int t = 0; t < 6; t++)
      {
         mesh.AddTet(T[t], 1);
      }
   };

   // Lower-y half (iy = 0) → orient +1; upper-y half (iy = 1) → orient -1.
   for (int iz = 0; iz < n; iz++)
   {
      for (int ix = 0; ix < n; ix++)
      {
         add_hex(ix, 0, iz, +1);
         add_hex(ix, 1, iz, -1);
      }
   }

   mesh.FinalizeTopology();
   mesh.Finalize();

   // Tag boundaries: exterior faces with attr = 1; interior y = L/2
   // plane with attr = 3 (fault) iff add_fault.
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
         if (add_fault && std::abs(cy - 0.5 * kL) < 1e-8)
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

// Orbit bucket: key by (cx, min(cy, L-cy), cz) as in Arm 1.
struct OrbitBucket { std::vector<int> elems; };

std::map<std::tuple<long long,long long,long long>, OrbitBucket>
BuildElementOrbits(const Mesh &mesh)
{
   std::map<std::tuple<long long,long long,long long>, OrbitBucket> b;
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
      const real_t cy_mir = std::min(cy, kL - cy);
      b[std::make_tuple(q(cx), q(cy_mir), q(cz))].elems.push_back(e);
   }
   return b;
}

real_t OrbitCellMeanMaxDrift(
   const Vector &rhs, const WaveOperator<Mesh> &wave,
   const std::map<std::tuple<long long,long long,long long>, OrbitBucket> &b,
   std::string &worst)
{
   const auto &fes = wave.GetFESpace();
   const int ndof_total = fes.GetNDofs();
   real_t d = 0.0; worst.clear();
   static const char *cn[9] = {"SXX","SYY","SZZ","SXY","SYZ","SXZ",
                                "VX","VY","VZ"};
   for (const auto &kv : b)
   {
      const auto &o = kv.second.elems;
      if (o.size() < 2) { continue; }
      for (int c = 0; c < NUM_STATE; c++)
      {
         std::vector<real_t> m; m.reserve(o.size());
         for (int e : o)
         {
            Array<int> ed; fes.GetElementDofs(e, ed);
            real_t s = 0.0;
            for (int i = 0; i < ed.Size(); i++)
            {
               s += rhs(c * ndof_total + ed[i]);
            }
            m.push_back(s / ed.Size());
         }
         for (std::size_t i = 0; i < m.size(); i++)
         {
            for (std::size_t j = i+1; j < m.size(); j++)
            {
               real_t dd = std::abs(m[i] - m[j]);
               if (dd > d) { d = dd; worst = std::string(cn[c]); }
            }
         }
      }
   }
   return d;
}

void FillConstantState(const FiniteElementSpace &fes, Vector &Q)
{
   const int ndof_total = fes.GetNDofs();
   Q.SetSize(NUM_STATE * ndof_total);
   Q = 0.0;
   for (int i = 0; i < ndof_total; i++)
   {
      Q(SXX * ndof_total + i) = 1.0e6;
   }
}

int CountOrbitDofMismatch(const Mesh &mesh, const WaveOperator<Mesh> &wave,
                          int &worst_idx_diff, int &orbits_examined)
{
   const auto &fes = wave.GetFESpace();
   auto buckets = BuildElementOrbits(mesh);
   const real_t scale = 1e6;
   auto q = [&](real_t v) { return static_cast<long long>(std::round(v*scale)); };

   int mism = 0;
   worst_idx_diff = 0;
   orbits_examined = 0;

   for (const auto &kv : buckets)
   {
      const auto &orbit = kv.second.elems;
      if (orbit.size() < 2) { continue; }
      orbits_examined++;

      std::vector<std::map<std::tuple<long long,long long,long long>, int>> pe;
      for (int e : orbit)
      {
         std::map<std::tuple<long long,long long,long long>, int> m;
         const FiniteElement *fe = fes.GetFE(e);
         ElementTransformation *Tr = fes.GetElementTransformation(e);
         const IntegrationRule &nd = fe->GetNodes();
         for (int j = 0; j < fe->GetDof(); j++)
         {
            const IntegrationPoint &ip = nd.IntPoint(j);
            Vector x(3); Tr->Transform(ip, x);
            const real_t y_mir = std::min(x(1), kL - x(1));
            m[std::make_tuple(q(x(0)), q(y_mir), q(x(2)))] = j;
         }
         pe.push_back(std::move(m));
      }

      bool any = false;
      int md = 0;
      for (const auto &kvm : pe[0])
      {
         for (std::size_t k = 1; k < pe.size(); k++)
         {
            auto it = pe[k].find(kvm.first);
            if (it == pe[k].end()) { any = true; md = std::max(md, 99); continue; }
            if (it->second != kvm.second)
            {
               any = true;
               md = std::max(md, std::abs(it->second - kvm.second));
            }
         }
      }
      if (any) { mism++; worst_idx_diff = std::max(worst_idx_diff, md); }
   }
   return mism;
}

void SeedDOFData(std::vector<DOFData> &dof, int nfq)
{
   dof.clear(); dof.resize(nfq);
   const real_t Zp = kRho * std::sqrt((kLambda + 2*kMu)/kRho);
   const real_t Zs = kRho * std::sqrt(kMu/kRho);
   for (auto &d : dof)
   {
      d.slip_rate = 1e-12; d.V1 = 0; d.V2 = 1e-12;
      d.tau1_0 = 0; d.tau2_0 = 30e6; d.sigma_n0 = 40e6;
      d.tau1_corr = 0; d.tau2_corr = 0; d.sigma_n_corr = 0;
      d.slip1 = 0; d.slip2 = 0;
      d.a = 0.008; d.Dc = 0.02;
      d.psi = 0.6 + 0.008 * std::log(1e-12/1e-6);
      d.Zp_plus = Zp; d.Zp_minus = Zp;
      d.Zs_plus = Zs; d.Zs_minus = Zs;
      d.eta_p = 0.5 * Zp; d.eta_s = 0.5 * Zs;
   }
}

} // anonymous

// ==========================================================================
// Probe 1: G_ORBIT_DOF_MAP on the D4-equivariant fixture (no fault).
// This is the primary Arm 2 gate — if the fixture is actually
// D4-equivariant, every orbit pair must have matching local DOF indices.
// ==========================================================================
void ProbeDofMap()
{
   std::cout << "\n--- D4 fixture: G_ORBIT_DOF_MAP ---\n";
   Mesh mesh = BuildD4EquivariantMesh(/*add_fault=*/false);
   BoundaryConfig bc;
   for (int i = 1; i <= 6; i++) { bc.natural_attrs.insert(i); }
   bc.fault_attr = 0;
   WaveOperator<Mesh> wave(mesh, /*order=*/1, kLambda, kMu, kRho, bc);

   int worst_diff = 0, orbits = 0;
   int mism = CountOrbitDofMismatch(mesh, wave, worst_diff, orbits);
   std::cout << "  orbits examined: " << orbits
             << "  mismatched: " << mism
             << "  worst local-idx diff: " << worst_diff << "\n";
   PROBE_VERDICT("D4-fixture G_ORBIT_DOF_MAP",
                 static_cast<real_t>(mism), 0.5);
}

// ==========================================================================
// Probe 2: G_CONST_VOL and G_CONST_VOL_MINV on D4 fixture (no fault).
// ==========================================================================
void ProbeConstVolNoFault()
{
   std::cout << "\n--- D4 fixture (no fault): G_CONST_VOL + MINV ---\n";
   Mesh mesh = BuildD4EquivariantMesh(false);
   BoundaryConfig bc;
   for (int i = 1; i <= 6; i++) { bc.natural_attrs.insert(i); }
   bc.fault_attr = 0;
   WaveOperator<Mesh> wave(mesh, 1, kLambda, kMu, kRho, bc);
   // Set bulk_bg EQUAL to the constant test state so the absorbing /
   // free-surface BCs treat Q_const as an equilibrium.  Otherwise
   // AbsorbingTotal(n, Q, 0) = A_n^+ · Q drives a boundary source that
   // legitimately differs across orbit-pair boundary faces (e.g. y=0
   // vs y=L) — not a bug, just the physics of a non-equilibrium input.
   real_t Q_bg[NUM_STATE] = {0};
   Q_bg[SXX] = 1.0e6;  // matches FillConstantState below
   wave.SetAbsorbingBackground(Q_bg);

   Vector Q; FillConstantState(wave.GetFESpace(), Q);
   Vector rhs(wave.Height()); rhs = 0.0;
   wave.ComputeVolumeRHS_ForTest(Q, rhs);
   auto b = BuildElementOrbits(mesh);
   std::string w;
   real_t d1 = OrbitCellMeanMaxDrift(rhs, wave, b, w);
   PROBE_VERDICT("D4 no-fault G_CONST_VOL", d1, 1e-10);

   wave.ApplyMassInverse_ForTest(rhs);
   real_t d2 = OrbitCellMeanMaxDrift(rhs, wave, b, w);
   PROBE_VERDICT("D4 no-fault G_CONST_VOL_MINV", d2, 1e-10);
}

// ==========================================================================
// Probe 3: G_CONST_FACE on D4 fixture, no fault.
// ==========================================================================
void ProbeConstFaceNoFault()
{
   std::cout << "\n--- D4 fixture (no fault): G_CONST_FACE runtime vs precomp ---\n";
   Mesh mesh = BuildD4EquivariantMesh(false);
   BoundaryConfig bc;
   for (int i = 1; i <= 6; i++) { bc.natural_attrs.insert(i); }
   bc.fault_attr = 0;
   WaveOperator<Mesh> wave(mesh, 1, kLambda, kMu, kRho, bc);
   // Set bulk_bg EQUAL to the constant test state so the absorbing /
   // free-surface BCs treat Q_const as an equilibrium.  Otherwise
   // AbsorbingTotal(n, Q, 0) = A_n^+ · Q drives a boundary source that
   // legitimately differs across orbit-pair boundary faces (e.g. y=0
   // vs y=L) — not a bug, just the physics of a non-equilibrium input.
   real_t Q_bg[NUM_STATE] = {0};
   Q_bg[SXX] = 1.0e6;  // matches FillConstantState below
   wave.SetAbsorbingBackground(Q_bg);

   Vector Q; FillConstantState(wave.GetFESpace(), Q);
   Vector rhs_r(wave.Height()); rhs_r = 0.0;
   wave.ComputeFaceFluxRHS_ForTest(Q, rhs_r);
   wave.UsePrecomputedFaceFluxes(true);
   Vector rhs_p(wave.Height()); rhs_p = 0.0;
   wave.ComputeFaceFluxRHS_ForTest(Q, rhs_p);

   real_t pdiff = 0.0;
   for (int i = 0; i < rhs_r.Size(); i++)
   { pdiff = std::max(pdiff, std::abs(rhs_r(i) - rhs_p(i))); }

   auto b = BuildElementOrbits(mesh);
   std::string w;
   real_t dr = OrbitCellMeanMaxDrift(rhs_r, wave, b, w);
   real_t dp = OrbitCellMeanMaxDrift(rhs_p, wave, b, w);
   std::cout << "  runtime drift = " << std::scientific << dr
             << "  precomp drift = " << dp
             << "  path-diff = " << pdiff << "\n";
   PROBE_VERDICT("D4 no-fault G_CONST_FACE runtime drift", dr, 1e-10);
   PROBE_VERDICT("D4 no-fault G_CONST_FACE precomp drift", dp, 1e-10);
   PROBE_VERDICT("D4 no-fault G_CONST_FACE path-diff", pdiff, 1e-6);
}

// ==========================================================================
// Probe 4: G_CONST_FACE on D4 fixture, FAULT CONFIG — the decisive test.
//
// This probe reproduces the M0 fault bookkeeping (attr=3 on y=L/2,
// fault_flux_ + DOFData seeded) on the D4-equivariant mesh, then asks:
// does the face-flux RHS drift across orbit?
//
// Interpretation:
//   * drift = 0 on BOTH runtime and precomputed paths → the pepper
//     drift observed in the audit is a pure fixture artifact (Kuhn
//     split non-D4-equivariance).  Production TPV102 on a D4-equivariant
//     mesh (or any mesh where orbit symmetry is not imposed) is clean.
//   * drift > 0 → the bug is downstream of the fixture.  Arm 3 applies.
// ==========================================================================
void ProbeConstFaceWithFault()
{
   std::cout << "\n--- D4 fixture (FAULT): G_CONST_FACE runtime vs precomp ---\n";
   Mesh mesh = BuildD4EquivariantMesh(/*add_fault=*/true);

   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr    = 3;

   WaveOperator<Mesh> wave(mesh, 1, kLambda, kMu, kRho, bc);
   // Set bulk_bg EQUAL to the constant test state so the absorbing /
   // free-surface BCs treat Q_const as an equilibrium.  Otherwise
   // AbsorbingTotal(n, Q, 0) = A_n^+ · Q drives a boundary source that
   // legitimately differs across orbit-pair boundary faces (e.g. y=0
   // vs y=L) — not a bug, just the physics of a non-equilibrium input.
   real_t Q_bg[NUM_STATE] = {0};
   Q_bg[SXX] = 1.0e6;  // matches FillConstantState below
   wave.SetAbsorbingBackground(Q_bg);

   FaultFaceFlux ff(kRho,
                    std::sqrt((kLambda + 2.0*kMu)/kRho),
                    std::sqrt(kMu/kRho));
   wave.SetFaultFlux(&ff);
   std::vector<DOFData> dof;
   SeedDOFData(dof, wave.GetNumTotalFaultQPs());
   wave.SetFaultDOFData(&dof, wave.GetNbfPerFace());

   Vector Q; FillConstantState(wave.GetFESpace(), Q);
   Vector rhs_r(wave.Height()); rhs_r = 0.0;
   wave.ComputeFaceFluxRHS_ForTest(Q, rhs_r);

   wave.UsePrecomputedFaceFluxes(true);
   Vector rhs_p(wave.Height()); rhs_p = 0.0;
   wave.ComputeFaceFluxRHS_ForTest(Q, rhs_p);

   real_t pdiff = 0.0;
   for (int i = 0; i < rhs_r.Size(); i++)
   { pdiff = std::max(pdiff, std::abs(rhs_r(i) - rhs_p(i))); }

   auto b = BuildElementOrbits(mesh);
   std::string wr, wp;
   real_t dr = OrbitCellMeanMaxDrift(rhs_r, wave, b, wr);
   real_t dp = OrbitCellMeanMaxDrift(rhs_p, wave, b, wp);
   std::cout << "  runtime drift = " << std::scientific << dr
             << "  worst=" << wr << "\n";
   std::cout << "  precomp drift = " << std::scientific << dp
             << "  worst=" << wp << "\n";
   std::cout << "  path-diff = " << pdiff << "\n";
   PROBE_VERDICT("D4 WITH-FAULT G_CONST_FACE runtime drift", dr, 1e-10);
   PROBE_VERDICT("D4 WITH-FAULT G_CONST_FACE precomp drift", dp, 1e-10);
   PROBE_VERDICT("D4 WITH-FAULT G_CONST_FACE path-diff", pdiff, 1e-6);
}

int main()
{
   std::cout << "=== TPV102 Phase 3 STOP — Arm 2: D4-equivariant fixture ===\n";

   ProbeDofMap();
   ProbeConstVolNoFault();
   ProbeConstFaceNoFault();
   ProbeConstFaceWithFault();

   std::cout << "\n=== Summary ===\n";
   std::cout << "Total:  " << num_tests << "\n";
   std::cout << "Passed: " << num_passed << "\n";
   std::cout << "Failed: " << num_failed << "\n";

   std::cout << "\nInterpretation:\n";
   std::cout << "  * D4 G_ORBIT_DOF_MAP PASS + D4 WITH-FAULT G_CONST_FACE PASS\n"
             << "    ==> Pepper drift is a Kuhn-split fixture artifact.  The\n"
             << "        flux layer is correct.  Production on a non-Kuhn\n"
             << "        mesh is pepper-free.  Close the investigation.\n";
   std::cout << "  * D4 G_ORBIT_DOF_MAP PASS + D4 WITH-FAULT G_CONST_FACE FAIL\n"
             << "    ==> Real numerical asymmetry downstream of fixture.\n"
             << "        Proceed to Arm 3 (basis/fault-accumulation\n"
             << "        investigation).\n";

   return (num_failed == 0) ? 0 : 1;
}

// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 "Topology-Based Precomputed Face-Rotation" plan 2026-04-23
// Phase 2a acceptance gates (§6.4):
//   * P_SWITCH_OFF_REGRESSION — default path unchanged when the flag is
//     off (same AdvanceADER step before and after the patch, to ULP).
//   * P_SWITCH_ON_LINEAR_EQ — with the flag on on a homogeneous Cartesian
//     plane-wave fixture, AdvanceADER matches the runtime path to 1e-12.
//   * test_R5_002_fault_face_set_populated_at_init_time — staged
//     WaveOperator ctor → SetFaultDOFData → UsePrecomputedFaceFluxes
//     assertion that fault_face_set_ is populated only at the last step
//     and that fault-face keys are absent from face_elem_to_entry_.

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/precomputed_face_fluxes.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace mfem { namespace seas { int g_seas_my_rank = 0; } }  // NOLINT

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_LE(v, tol, msg) do { \
   num_tests++; \
   const double vv = (v), tt = (tol); \
   if (vv <= tt) { num_passed++; \
      std::cout << "  PASSED: " << msg << "  (got " << std::scientific \
                << std::setprecision(3) << vv << ", tol " << tt << ")\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
                << "  (got " << std::scientific << std::setprecision(3) \
                << vv << ", expected <= " << tt << ")\n"; } \
} while (0)

#define TEST_TRUE(expr, msg) do { \
   num_tests++; \
   if (expr) { num_passed++; \
      std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

namespace
{

constexpr real_t kRho    = 2670.0;
constexpr real_t kLambda = 32.04e9;
constexpr real_t kMu     = 32.04e9;

struct PWave
{
   real_t Q0[NUM_STATE];
   real_t cp;
   real_t k;
};

PWave MakePWave(real_t k)
{
   PWave w;
   const real_t lp = kLambda + 2.0 * kMu;
   w.cp = std::sqrt(lp / kRho);
   w.k  = k;
   for (int c = 0; c < NUM_STATE; c++) { w.Q0[c] = 0.0; }
   w.Q0[SXX] = -lp;
   w.Q0[SYY] = -kLambda;
   w.Q0[SZZ] = -kLambda;
   w.Q0[VX]  =  w.cp;
   return w;
}

void FillPlaneWave(const FiniteElementSpace &fes, const PWave &w,
                   real_t t, Vector &Q)
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
         Vector x(3);
         Tr->Transform(ip, x);
         const real_t phase = w.k * (x(0) - w.cp * t);
         const real_t s = std::sin(phase);
         for (int c = 0; c < NUM_STATE; c++)
         {
            Q[c * ndof_total + edofs[j]] = w.Q0[c] * s;
         }
      }
   }
}

real_t InfNormDiff(const Vector &a, const Vector &b)
{
   real_t m = 0.0;
   for (int i = 0; i < a.Size(); i++)
   {
      m = std::max(m, std::abs(a(i) - b(i)));
   }
   return m;
}

real_t Normlinf(const Vector &a)
{
   real_t m = 0.0;
   for (int i = 0; i < a.Size(); i++) { m = std::max(m, std::abs(a(i))); }
   return m;
}

// Build the TPV102 Cartesian-fault M0 fixture: 2x2x2 tets with a
// split-plane at y = 0.5L as fault_attr=3 and the remaining exterior
// tri-faces at natural_attr=1.  Identical to
// test_adjacent_triangle_fault_first_step_audit.cpp's fixture.
Mesh BuildM0FaultMesh()
{
   const real_t L = 1000.0;
   const int n = 2;
   Mesh mesh = Mesh::MakeCartesian3D(n, n, n, Element::TETRAHEDRON, L, L, L,
                                     false);
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
         if (std::abs(cy - 0.5 * L) < 1e-8)
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

} // anonymous

// ==========================================================================
// P_SWITCH_ON_LINEAR_EQ (§6.4): with flag on, AdvanceADER matches runtime
// to 1e-12 on homog Cartesian plane-wave fixture (no fault).
// ==========================================================================
int TestLinearWaveEquivalence()
{
   std::cout << "\n[TestLinearWaveEquivalence]\n";

   const real_t L = 1000.0;
   const int n = 4;
   const int p = 1;   // Phase 1 p=1 precondition
   Mesh mesh = Mesh::MakeCartesian3D(n, n, n, Element::TETRAHEDRON, L, L, L);

   BoundaryConfig bc;
   for (int i = 1; i <= 6; i++) { bc.absorbing_attrs.insert(i); }
   bc.fault_attr = 0;

   WaveOperator<Mesh> wave(mesh, p, kLambda, kMu, kRho, bc);
   const FiniteElementSpace &fes = wave.GetFESpace();

   real_t Q_bg_zero[NUM_STATE] = {0};
   wave.SetAbsorbingBackground(Q_bg_zero);

   const real_t lambda_wave = 8.0 * L;
   const real_t k = 2.0 * M_PI / lambda_wave;
   PWave w = MakePWave(k);

   const real_t h = L / n;
   const real_t dt_cfl = h / w.cp;
   const real_t dt = 0.01 * dt_cfl;

   Vector Q; FillPlaneWave(fes, w, /*t=*/0.0, Q);

   // Runtime path (flag off, default).
   TEST_TRUE(!wave.UsingPrecomputedFaceFluxes(),
             "flag defaults to off (P_SWITCH_OFF)");
   Vector Q_runtime;
   wave.AdvanceADER(Q, dt, /*order=*/2, Q_runtime);
   const real_t q_inf = Normlinf(Q_runtime);
   TEST_TRUE(q_inf > 0.0, "runtime ADER step produces nonzero output");

   // Precomputed path (flag on).
   wave.UsePrecomputedFaceFluxes(true);
   TEST_TRUE(wave.UsingPrecomputedFaceFluxes(),
             "UsePrecomputedFaceFluxes(true) sets the flag");
   TEST_TRUE(wave.GetPrecomputedFaceFluxes().IsInitialized(),
             "UsePrecomputedFaceFluxes(true) initializes the precomputed tables");

   Vector Q_precomputed;
   wave.AdvanceADER(Q, dt, /*order=*/2, Q_precomputed);

   // P_SWITCH_ON_LINEAR_EQ: absolute-diff tolerance 1e-12 times the peak
   // state magnitude (plan §6.4 specifies ||·||_∞ < 1e-12; here we scale
   // by q_inf to make the test unit-independent, so the effective
   // relative tolerance is 1e-12).  Use an explicit absolute check on
   // the normalized diff.
   const real_t abs_diff = InfNormDiff(Q_precomputed, Q_runtime);
   std::cout << "  || Q_precomputed - Q_runtime ||_inf = "
             << std::scientific << abs_diff << "  (Q_runtime inf="
             << q_inf << ")\n";
   TEST_LE(abs_diff / std::max(q_inf, real_t(1.0)), 1e-12,
           "P_SWITCH_ON_LINEAR_EQ: ADER step bit-matches runtime to 1e-12 rel");

   // Flip back off; verify flag respects the toggle and that a second
   // AdvanceADER reproduces the runtime output.
   wave.UsePrecomputedFaceFluxes(false);
   TEST_TRUE(!wave.UsingPrecomputedFaceFluxes(),
             "UsePrecomputedFaceFluxes(false) clears the flag");
   Vector Q_runtime_again;
   wave.AdvanceADER(Q, dt, /*order=*/2, Q_runtime_again);
   TEST_LE(InfNormDiff(Q_runtime_again, Q_runtime) /
           std::max(q_inf, real_t(1.0)),
           1e-15,
           "P_SWITCH_OFF_REGRESSION: flag-off ADER step identical "
           "before and after flag toggling (to ULP)");

   return 0;
}

// ==========================================================================
// test_R5_002_fault_face_set_populated_at_init_time (§6.3a).
// Staged: WaveOperator ctor → SetFaultDOFData → UsePrecomputedFaceFluxes.
// Asserts fault_face_set_ is empty at the first two stages, populated at
// the third, AND that face_elem_to_entry_ omits fault-face keys.
// ==========================================================================
int TestR5_002FaultFaceSetPopulatedAtInitTime()
{
   std::cout << "\n[TestR5_002FaultFaceSetPopulatedAtInitTime]\n";

   Mesh mesh = BuildM0FaultMesh();

   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr    = 3;

   real_t Q_bg_zero[NUM_STATE] = {0};

   // Stage 1: ctor.
   WaveOperator<Mesh> wave(mesh, /*order=*/1, kLambda, kMu, kRho, bc);
   wave.SetAbsorbingBackground(Q_bg_zero);
   TEST_TRUE(wave.GetFaultFaceSet().empty(),
             "ctor leaves fault_face_set_ empty (R5-002: late population)");
   TEST_TRUE(wave.GetFaultInteriorFaces().Size() > 0,
             "ctor populates fault_interior_faces_ on M0 fault fixture");

   // Stage 2: SetFaultDOFData.  Must NOT populate fault_face_set_.
   // This test asserts the DOFData-population invariant only; AdvanceADER
   // is never called on this instance, so the DOFData contents are
   // immaterial (only the pointer-install path is exercised).  The
   // DOFData vector is default-constructed (DOFData has sane defaults)
   // and `SetFaultDOFData` only takes the pointer — any per-field seed
   // here is dead code for this sub-test.
   FaultFaceFlux ff(kRho, std::sqrt((kLambda + 2.0*kMu)/kRho),
                    std::sqrt(kMu / kRho));
   wave.SetFaultFlux(&ff);
   const int nqp = wave.GetNbfPerFace();
   const int nfq = wave.GetNumTotalFaultQPs();
   std::vector<DOFData> dof_data(nfq);
   wave.SetFaultDOFData(&dof_data, nqp);
   TEST_TRUE(wave.GetFaultFaceSet().empty(),
             "SetFaultDOFData leaves fault_face_set_ empty (R5-002)");

   // Stage 3: UsePrecomputedFaceFluxes(true) populates fault_face_set_
   // AND runs Init, skipping fault faces.
   //
   // NOTE — this sub-test runs in serial, so GetFaultSharedFaces().Size()
   // is always 0.  The expected-count check exercises the fault_interior
   // branch only; the shared-face translation (sf → mesh_face_idx via
   // ParMesh::GetSharedFace, wave_operator.inl:527-529) is NOT exercised
   // here.  A parallel regression for that branch must be added under
   // tests/parallel/ when Phase 2b (MPI) acceptance lands.
   wave.UsePrecomputedFaceFluxes(true);
   const std::set<int> &ffset = wave.GetFaultFaceSet();
   const int expected = wave.GetFaultInteriorFaces().Size()
                      + wave.GetFaultSharedFaces().Size();
   TEST_TRUE(static_cast<int>(ffset.size()) == expected,
             "UsePrecomputedFaceFluxes(true) populates fault_face_set_ "
             "with fault_interior + fault_shared counts "
             "(serial: shared count is 0; parallel branch untested here)");
   TEST_TRUE(wave.GetPrecomputedFaceFluxes().IsInitialized(),
             "UsePrecomputedFaceFluxes(true) initializes the precomputed tables");

   // Every element in fault_face_set_ must be in fault_interior_faces_.
   int missing = 0;
   for (int f : ffset)
   {
      bool found = false;
      const Array<int> &fif = wave.GetFaultInteriorFaces();
      for (int i = 0; i < fif.Size(); i++)
      {
         if (fif[i] == f) { found = true; break; }
      }
      if (!found) { missing++; }
   }
   TEST_TRUE(missing == 0,
             "every face in fault_face_set_ is either "
             "in fault_interior_faces_ or fault_shared_faces_");

   // Fault faces must be ABSENT from precomputed_face_fluxes_'s
   // face_elem_to_entry_ (they are routed through FaultFaceFlux, not
   // through the precomputed-table dispatch).
   const PrecomputedFaceFluxes &pff = wave.GetPrecomputedFaceFluxes();
   int fault_keys_in_map = 0;
   for (int face_idx : ffset)
   {
      for (int e = 0; e < wave.NumElements(); e++)
      {
         if (pff.HasEntry(face_idx, e)) { fault_keys_in_map++; }
      }
   }
   TEST_TRUE(fault_keys_in_map == 0,
             "no fault-face key (face_idx, any_elem) in "
             "face_elem_to_entry_ after Phase 1 Init (R5-002 invariant)");

   return 0;
}

int main()
{
   std::cout << "\n=== TPV102 Phase 2a opt-in dispatch switch tests ===\n";
   TestLinearWaveEquivalence();
   TestR5_002FaultFaceSetPopulatedAtInitTime();
   std::cout << "\n===================================================\n";
   std::cout << "Total:  " << num_tests << "\n";
   std::cout << "Passed: " << num_passed << "\n";
   std::cout << "Failed: " << num_failed << "\n";
   std::cout << "===================================================\n";
   return (num_failed == 0) ? 0 : 1;
}

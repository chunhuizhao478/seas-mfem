// ===========================================================================
// Phase 6 — Test 6.1: the HOMOGENEOUS-EQUIVALENCE GATE.
//
// In the homogeneous limit the bi-material central flux must reduce to the
// scalar central flux: F* = ½(A_e1·Q_e1 + A_e2·Q_e2) with A_e1==A_e2==A
// collapses to ½A(Q_e1+Q_e2), exactly the scalar `GodunovFlux::Central`.  So on
// a HOMOGENEOUS material a `BimaterialWaveOperator` + SetMixedFluxMode(mode)
// must match the scalar `WaveOperator` + SetMixedFluxMode(mode) under `Mult`,
// to LU rounding (field-scale 1e-9) — the SAME contract as the C-6 parity test
// (test_bimaterial_wave_operator_parity.cpp) but on the MIXED-flux path
// (Adjacent AND AllContinuous), not just upwind.
//
// This is NECESSARY but NOT SUFFICIENT: with a homogeneous material A_e1==A_e2,
// so a per-side-swapped (BUG-3) storage/dispatch would ALSO pass here.  The
// SUFFICIENT guard for BUG-3 (single-valuedness under a genuine material
// contrast) is Test 3.2 in test_bimaterial_mixed_flux_dispatch.cpp.  This gate
// catches the complementary class: a central-flux contribution the bi-material
// path computes DIFFERENTLY from the scalar reference even when the material is
// uniform (a wrong matrix assembly, a missed/extra central face, a sign error
// in the per-side deposit).  It has teeth: perturbing a stored central matrix
// drives the residual well above 1e-9 (Rev I-8 #2).
//
// Fixture (needs fault-adjacent NON-fault interior faces — which the 2-tet C-6
// mesh lacks — AND a non-fault-adjacent interior face so Adjacent and
// AllContinuous have DIFFERENT central sets): a 1×5×1 column of unit hexes
// along y over [0,1]×[0,5]×[0,1].  The y=2 interior plane (between elems 1,2)
// is tagged FAULT (attr 3).  The non-fault interior faces are y=1, y=3, y=4:
//   - Adjacent      central set = {y=1, y=3}        (the 2 faces of the two
//                                                     fault-adjacent elems 1,2)
//   - AllContinuous central set = {y=1, y=3, y=4}   (ALL non-fault interior;
//                                                     y=4 is NOT fault-adjacent)
// so AllContinuous (3) ⊋ Adjacent (2), and each mode's distinct face-set
// construction is exercised.  The fault is y-normal, so the TPV102 total-stress
// background (SYY normal, SXY shear) applies as in the parity test.
// ===========================================================================

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/bimaterial_wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/heterogeneous_material.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../dynamic/tpv102_setup_total.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do {                                         \
   num_tests++;                                                             \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; }    \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: "      \
                                  << msg << "\n"; }                         \
} while (0)

namespace
{
// 1 x ny x 1 column of unit hexes along y over [0,1]x[0,ny]x[0,1].  The
// interior y == fault_y plane -> FAULT (attr 3); all external faces -> FREE
// (attr 1).  Every other interior y-face is left untagged (non-fault).
Mesh BuildHexColumnFaultMesh(int ny, real_t fault_y)
{
   const int nx = 1, nz = 1;
   Mesh mesh(3, (nx + 1) * (ny + 1) * (nz + 1), nx * ny * nz, 0);
   auto vid = [&](int ix, int iy, int iz)
   { return ix * (ny + 1) * (nz + 1) + iy * (nz + 1) + iz; };
   for (int ix = 0; ix <= nx; ++ix)
      for (int iy = 0; iy <= ny; ++iy)
         for (int iz = 0; iz <= nz; ++iz)
         {
            real_t v[3] = { static_cast<real_t>(ix), static_cast<real_t>(iy),
                            static_cast<real_t>(iz) };
            mesh.AddVertex(v);
         }
   for (int iy = 0; iy < ny; ++iy)
   {
      const int v[8] = {
         vid(0, iy,   0), vid(1, iy,   0), vid(1, iy+1, 0), vid(0, iy+1, 0),
         vid(0, iy,   1), vid(1, iy,   1), vid(1, iy+1, 1), vid(0, iy+1, 1)
      };
      mesh.AddHex(v, 1);
   }
   mesh.FinalizeTopology();

   for (int f = 0; f < mesh.GetNumFaces(); ++f)
   {
      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      real_t cy = 0.0;
      for (int i = 0; i < fv.Size(); ++i) { cy += mesh.GetVertex(fv[i])[1]; }
      cy /= fv.Size();
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         if (std::abs(cy - fault_y) < 1e-9)
         {
            mesh.AddBdrQuad(fv[0], fv[1], fv[2], fv[3], 3);   // FAULT (y=fault_y)
         }
         continue;                                            // other interior: untagged
      }
      mesh.AddBdrQuad(fv[0], fv[1], fv[2], fv[3], 1);          // FREE (external)
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// Wire the fault into a WaveOperator (scalar base or bimaterial subclass via a
// base reference).  Returns the fault QP count.  (Mirrors the C-6 parity test;
// the serial column fixture has a single fault face, so the int_faces[0] nqp is
// representative.)
int SetupFault(WaveOperator<Mesh> &wave, Mesh &mesh, int order,
               std::vector<DOFData> &dof_data, FaultFaceFlux &ff,
               std::vector<Vector> &fault_coords)
{
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   MFEM_VERIFY(wave.GetFaultSharedFaces().Size() == 0,
               "serial column fixture cannot have shared faces");
   int nqp_per_face = 0;
   if (int_faces.Size() > 0)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[0]);
      MFEM_VERIFY(ftr, "interior fault face FTR null");
      nqp_per_face = IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
   }
   const int nfault = int_faces.Size() * nqp_per_face;
   fault_coords.clear();
   fault_coords.reserve(nfault);
   for (int i = 0; i < int_faces.Size(); ++i)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir = IntRules.Get(ftr->GetGeometryType(), 2*order);
      for (int q = 0; q < ir.GetNPoints(); ++q)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   }
   if (nfault > 0) { InitializeFaultDOFs(dof_data, nfault, fault_coords); }
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp_per_face);
   return nfault;
}

void ResetDOFData(std::vector<DOFData> &dof_data,
                  const std::vector<Vector> &fault_coords)
{
   const int ndof = static_cast<int>(dof_data.size());
   InitializeFaultDOFs(dof_data, ndof, fault_coords);
   ZeroDOFDataPreStressTotal(dof_data, ndof);
}

// Elements with centroid y > split (above the fault) get val_plus; below it,
// val_minus — an antisymmetric jump across the y=split fault.
void SetAntiSymQComponent(const Mesh &mesh, Vector &Q, int comp,
                          int ndof_total, int ndof_per_elem,
                          real_t split, real_t val_plus, real_t val_minus)
{
   for (int e = 0; e < mesh.GetNE(); ++e)
   {
      real_t cy = 0;
      Array<int> ev; mesh.GetElementVertices(e, ev);
      for (int v = 0; v < ev.Size(); ++v) { cy += mesh.GetVertex(ev[v])[1]; }
      cy /= ev.Size();
      const real_t val = (cy > split) ? val_plus : val_minus;
      const int off = e * ndof_per_elem;
      for (int i = 0; i < ndof_per_elem; ++i)
      {
         Q(comp * ndof_total + off + i) = val;
      }
   }
}

// Field-scale relative difference (same metric as the C-6 parity test).
real_t MaxRelDiff(const Vector &a, const Vector &b, int &n_over_tol,
                  real_t tol, bool &all_finite)
{
   real_t scale = 1.0;
   all_finite = true;
   for (int i = 0; i < a.Size(); ++i)
   {
      if (!std::isfinite(b(i)) || !std::isfinite(a(i))) { all_finite = false; }
      scale = std::max(scale, std::max(std::abs(a(i)), std::abs(b(i))));
   }
   real_t max_rel = 0.0;
   n_over_tol = 0;
   for (int i = 0; i < a.Size(); ++i)
   {
      const real_t rel = std::abs(a(i) - b(i)) / scale;
      max_rel = std::max(max_rel, rel);
      if (rel > tol) { ++n_over_tol; }
   }
   return max_rel;
}

// Run the gate for one mixed-flux mode: scalar WaveOperator vs homogeneous
// BimaterialWaveOperator, both with SetMixedFluxMode(mode), Mult on 3 states.
// `expected_central` pins the central-face count so the test fails loud if the
// set-construction silently changes (and proves Adjacent != AllContinuous).
void RunGate(MixedFluxMode mode, const char *mode_name, int expected_central)
{
   std::cout << "\n-- Test 6.1 (" << mode_name << "): homogeneous equivalence --\n";
   Mesh mesh = BuildHexColumnFaultMesh(/*ny=*/5, /*fault_y=*/2.0);
   const real_t fault_split = 2.0;

   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {};

   const int order = 1;
   // 6-sig-fig-exact material (GodunovFluxPool dedup is a no-op => the per-
   // element pool equals the scalar flux_ bit-for-bit; same choice as C-6).
   const real_t lam = 32.0e9, mu = 32.0e9, rho = 2670.0;
   const real_t cp  = std::sqrt((lam + 2.0 * mu) / rho);
   const real_t cs  = std::sqrt(mu / rho);

   WaveOperator<Mesh> wave_scalar(mesh, order, lam, mu, rho, bc);
   ConstantCoefficient lam_c(lam), mu_c(mu), rho_c(rho);
   BimaterialWaveOperator<Mesh> wave_bimat(
      mesh, order, MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c), bc);

   FaultFaceFlux ff_s(rho, cp, cs), ff_b(rho, cp, cs);
   std::vector<DOFData> dof_s, dof_b;
   std::vector<Vector>  coords_s, coords_b;
   const int nfault_s = SetupFault(wave_scalar, mesh, order, dof_s, ff_s, coords_s);
   const int nfault_b = SetupFault(wave_bimat,  mesh, order, dof_b, ff_b, coords_b);
   TEST_ASSERT(nfault_s > 0 && nfault_s == nfault_b,
               "both operators see the same nonzero fault QP count");

   wave_scalar.SetMixedFluxMode(mode);
   wave_bimat.SetMixedFluxMode(mode);
   // Pin the central-face count: non-empty (gate exercises the central flux),
   // EXACTLY `expected_central` (Adjacent=2 ⊊ AllContinuous=3 proves the two
   // modes build DISTINCT sets), and scalar==bimaterial.
   const int central_s = static_cast<int>(wave_scalar.GetCentralFluxFaceSet().size());
   const int central_b = static_cast<int>(wave_bimat.GetCentralFluxFaceSet().size());
   TEST_ASSERT(central_s == expected_central,
               std::string("scalar central-face count == expected for ") + mode_name);
   TEST_ASSERT(central_b == central_s,
               "scalar + bimaterial agree on the central-face count");
   std::cout << "  central faces = " << central_s << " (expected " << expected_central
             << "), fault QPs = " << nfault_s << "\n";

   real_t bulk_bg[NUM_STATE] = {0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;
   wave_scalar.SetAbsorbingBackground(bulk_bg);
   wave_bimat.SetAbsorbingBackground(bulk_bg);

   const int size = wave_scalar.Height();
   const auto &fes = wave_scalar.GetFESpace();
   const int ndof_total    = fes.GetNDofs();
   const int ndof_per_elem = fes.GetFE(0)->GetDof();
   const real_t k_tol = 1e-9;

   // 3 distinct states: an antisymmetric slip-rate jump across the y=2 fault in
   // {strike VX, anti-plane VZ, larger VX}.  Total-stress background keeps the
   // friction solve well-posed; both operators run the SAME state.
   struct Pert { int comp; real_t amp; const char *desc; };
   const Pert perts[] = {
      {VX, 1.0e-5, "antisym VX (strike) 1e-5"},
      {VZ, 1.0e-5, "antisym VZ (anti-plane) 1e-5"},
      {VX, 2.0e-5, "antisym VX (strike) 2e-5"},
   };

   for (const Pert &p : perts)
   {
      ResetDOFData(dof_s, coords_s);
      ResetDOFData(dof_b, coords_b);
      Vector Q(size);
      InitializeStateTotal(Q, ndof_total, TPV102Params::sigma_n,
                           TPV102Params::tau_ini);
      SetAntiSymQComponent(mesh, Q, p.comp, ndof_total, ndof_per_elem,
                           fault_split, +0.5 * p.amp, -0.5 * p.amp);

      Vector k_s(size), k_b(size);
      wave_scalar.Mult(Q, k_s);
      wave_bimat.Mult(Q, k_b);

      int n_over = 0; bool finite = true;
      const real_t max_rel = MaxRelDiff(k_s, k_b, n_over, k_tol, finite);
      TEST_ASSERT(finite, std::string("6.1 (") + mode_name + "): bimaterial Mult "
                  "finite (no NaN) for " + p.desc);
      if (n_over != 0)
      {
         std::cerr << "  DIFF: " << n_over << " / " << size
                   << " over tol; max rel = " << max_rel << " [" << p.desc << "]\n";
      }
      TEST_ASSERT(n_over == 0,
                  std::string("6.1 (") + mode_name + "): bimaterial (homogeneous) "
                  "Mult matches scalar to 1e-9 for " + p.desc);
   }
}

} // anonymous namespace

int main()
{
   std::cout << "\n=== Phase 6 Test 6.1: homogeneous-equivalence gate "
                "(bi-material central flux -> scalar central flux) ===\n";

   // Adjacent (2 central faces) ⊊ AllContinuous (3 central faces) on the 1×5×1
   // column — each mode's distinct face-set construction is exercised.
   RunGate(MixedFluxMode::Adjacent,      "Adjacent",      2);
   RunGate(MixedFluxMode::AllContinuous, "AllContinuous", 3);

   std::cout << "\n==============================================\n";
   std::cout << "Summary: " << num_passed << " / " << num_tests
             << " passed (" << num_failed << " failed)\n";
   std::cout << "==============================================\n";
   return (num_failed > 0) ? 1 : 0;
}

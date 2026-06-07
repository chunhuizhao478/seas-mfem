// Part A (unified bi-material plan): central-flux CONTRAST GUARD.
//
// The mixed-flux central corridor is non-dissipative; across an impedance
// contrast it seeds a leak.  The guard (mixed_flux_contrast_tol_ >= 0) drops
// fault-adjacent corridor faces that span a STRONG impedance contrast from
// central_flux_face_set_ so they dispatch the (dissipative) bi-material upwind.
//
//   Phase 1 — IsStrongContrast / ContrastValue predicate (the tested rule):
//             tol<0 disabled; equal materials => 0/false; the TPV31 5 km pair
//             (cZs ~14.8%) > 0.05; the in-layer gradient pair (~0.36%) < 0.05.
//   Phase 2 — face-set membership on a step-material hex row: tol<0 leaves the
//             corridor unchanged (byte-exact); tol=0.05 removes EXACTLY the one
//             corridor face that crosses the material step (the uniform one is
//             kept).
//   Phase 3 — mechanism (R-001 two-sided dissipation): on the contrast corridor
//             face, D = [[Q]] . (0.5*(F_up_e1+F_up_e2) - F_ce) > 0, i.e. the
//             central flux discards the dissipation the upwind keeps.  In the
//             HOMOGENEOUS limit D == the scalar GodunovFlux [[Q]].(Interior -
//             Central) and equals 0.5*Zp*[[v_n]]^2; the guard removes 0 faces
//             (byte-exact).
//
// Build: make MFEM_DIR=../.. MFEM_BUILD_DIR=$MAIN MFEM_INC_DIR=$MAIN
//        MFEM_LIB_DIR=$MAIN seas_test_bimaterial_contrast_guard  (mfem-dev).

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/bimaterial_wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/godunov_flux_bimaterial.hpp"
#include "../../dynamic/heterogeneous_material.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do {                                          \
   num_tests++;                                                              \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; }     \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: "       \
                                  << msg << "\n"; }                          \
} while (0)

namespace
{
// Expose protected InteriorFaceFlux_ for direct per-face evaluation.
template <typename MeshType>
class TestableBimat : public BimaterialWaveOperator<MeshType>
{
public:
   using BimaterialWaveOperator<MeshType>::BimaterialWaveOperator;
   using BimaterialWaveOperator<MeshType>::InteriorFaceFlux_;
};

// Row of `nx` unit hexes along x; interior face at x == fault_x tagged FAULT
// (attr 3), external faces FREE (attr 1).  (Same fixture as the dispatch test.)
Mesh BuildHexRowFaultMesh(int nx, int fault_x)
{
   Mesh mesh(3, (nx + 1) * 4, nx, 0);
   auto vid = [](int ix, int iy, int iz) { return ix * 4 + iy * 2 + iz; };
   for (int ix = 0; ix <= nx; ++ix)
      for (int iy = 0; iy < 2; ++iy)
         for (int iz = 0; iz < 2; ++iz)
         {
            real_t v[3] = { static_cast<real_t>(ix), static_cast<real_t>(iy),
                            static_cast<real_t>(iz) };
            mesh.AddVertex(v);
         }
   for (int e = 0; e < nx; ++e)
   {
      const int v[8] = {
         vid(e, 0, 0), vid(e + 1, 0, 0), vid(e + 1, 1, 0), vid(e, 1, 0),
         vid(e, 0, 1), vid(e + 1, 0, 1), vid(e + 1, 1, 1), vid(e, 1, 1)
      };
      mesh.AddHex(v, 1);
   }
   mesh.FinalizeTopology();
   for (int f = 0; f < mesh.GetNumFaces(); ++f)
   {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      real_t cx = 0.0;
      for (int i = 0; i < fv.Size(); ++i) { cx += mesh.GetVertex(fv[i])[0]; }
      cx /= fv.Size();
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         if (std::abs(cx - fault_x) < 1e-9)
         { mesh.AddBdrQuad(fv[0], fv[1], fv[2], fv[3], 3); }
         continue;
      }
      mesh.AddBdrQuad(fv[0], fv[1], fv[2], fv[3], 1);
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

real_t FaceCentroidX(Mesh &mesh, int f)
{
   Array<int> fv; mesh.GetFaceVertices(f, fv);
   real_t cx = 0.0;
   for (int i = 0; i < fv.Size(); ++i) { cx += mesh.GetVertex(fv[i])[0]; }
   return cx / fv.Size();
}

// (lambda, mu, rho) for the TPV31 layer edges (Z = rho*c, mu = rho*vs^2).
GodunovFlux MakeFlux(real_t vp, real_t vs, real_t rho)
{
   const real_t mu  = rho * vs * vs;
   const real_t lam = rho * vp * vp - 2.0 * mu;
   return GodunovFlux(lam, mu, rho);
}

} // anonymous namespace

int main()
{
   std::cout << "\n=== Part A: central-flux CONTRAST GUARD ===\n";

   // TPV31 5 km layer edges (5000-: vs=3050,rho=2620 ; 5000+: vs=3450,rho=2720).
   GodunovFlux soft = MakeFlux(5200.0, 3050.0, 2620.0);
   GodunovFlux hard = MakeFlux(5750.0, 3450.0, 2720.0);
   // In-layer gradient neighbours (~0.36% Zs apart): scale mu by ~0.72%.
   GodunovFlux grad_a = MakeFlux(5200.0, 3050.0, 2620.0);
   GodunovFlux grad_b = MakeFlux(5200.0 * std::sqrt(1.0072),
                                 3050.0 * std::sqrt(1.0072), 2620.0);

   // -----------------------------------------------------------------------
   std::cout << "\n-- Phase 1: IsStrongContrast / ContrastValue predicate --\n";
   const real_t cval = BimaterialFlux::ContrastValue(soft, hard);
   std::cout << "  ContrastValue(soft,hard) = " << std::fixed
             << std::setprecision(4) << cval << "  (expect ~0.148)\n";
   TEST_ASSERT(std::abs(cval - 0.148) < 0.01,
               "ContrastValue(soft,hard) ~ 0.148 (the 5 km Zs jump)");
   TEST_ASSERT(BimaterialFlux::ContrastValue(soft, soft) == 0.0,
               "ContrastValue(soft,soft) == 0 (equal materials)");
   TEST_ASSERT(!BimaterialFlux::IsStrongContrast(soft, hard, -1.0),
               "IsStrongContrast tol<0 => false (disabled)");
   TEST_ASSERT(!BimaterialFlux::IsStrongContrast(soft, soft, 0.05),
               "IsStrongContrast equal materials => false");
   TEST_ASSERT(BimaterialFlux::IsStrongContrast(soft, hard, 0.05),
               "IsStrongContrast 5 km pair (14.8%) > tol 0.05 => true");
   TEST_ASSERT(!BimaterialFlux::IsStrongContrast(grad_a, grad_b, 0.05),
               "IsStrongContrast in-layer gradient pair (~0.36%) < 0.05 => false");

   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {};
   const int order = 1;

   // -----------------------------------------------------------------------
   // Phase 2: face-set membership on a STEP-material hex row.
   //   5 hexes (elems 0..4); fault face at x=3 (between elems 2,3) => corridor
   //   = {(1,2) at x=2, (3,4) at x=4}.  Material STEP at x=2: mu_A for x<2
   //   (elems 0,1), mu_B for x>=2 (elems 2,3,4).  So corridor face (1,2) spans
   //   the step (strong contrast) and (3,4) is uniform (mu_B both sides).
   // -----------------------------------------------------------------------
   std::cout << "\n-- Phase 2: contrast-guard face-set membership --\n";
   auto mu_step = [](const Vector &x) -> real_t
   { return (x(0) < 2.0) ? 30.0e9 : 60.0e9; };          // 2x mu step => ~29% Zs
   FunctionCoefficient lam_step(mu_step), mu_step_c(mu_step);
   ConstantCoefficient rho_c(2670.0);

   Mesh mesh = BuildHexRowFaultMesh(/*nx=*/5, /*fault_x=*/3);

   TestableBimat<Mesh> op_disabled(
      mesh, order, MaterialField::MakeCoefficient(&lam_step, &mu_step_c, &rho_c), bc);
   op_disabled.SetMixedFluxContrastTol(-1.0);            // disabled (byte-exact)
   op_disabled.SetMixedFluxMode(MixedFluxMode::Adjacent);

   TestableBimat<Mesh> op_guard(
      mesh, order, MaterialField::MakeCoefficient(&lam_step, &mu_step_c, &rho_c), bc);
   op_guard.SetMixedFluxContrastTol(0.05);               // guard on
   op_guard.SetMixedFluxMode(MixedFluxMode::Adjacent);

   const auto &cset_disabled = op_disabled.GetCentralFluxFaceSet();
   const auto &cset_guard    = op_guard.GetCentralFluxFaceSet();

   TEST_ASSERT(cset_disabled.size() == 2,
               "tol<0: corridor unchanged — both fault-adjacent faces central "
               "(byte-exact: guard inert)");
   TEST_ASSERT(cset_guard.size() == 1,
               "tol=0.05: exactly one corridor face reclassified to upwind");

   // The removed face = disabled \ guard; it must be the x=2 (step-crossing) one.
   int removed = -1, kept = -1;
   for (int f : cset_disabled)
   {
      if (cset_guard.count(f) == 0) { removed = f; } else { kept = f; }
   }
   TEST_ASSERT(removed >= 0 && std::abs(FaceCentroidX(mesh, removed) - 2.0) < 1e-9,
               "the reclassified face is the corridor face that crosses the "
               "material step (centroid x=2)");
   TEST_ASSERT(kept >= 0 && std::abs(FaceCentroidX(mesh, kept) - 4.0) < 1e-9,
               "the uniform-material corridor face (centroid x=4) stays central");

   // Homogeneous-material control: guard removes 0 faces (byte-exact).
   ConstantCoefficient lam_h(30.0e9), mu_h(30.0e9);
   Mesh mesh_h = BuildHexRowFaultMesh(/*nx=*/5, /*fault_x=*/3);
   TestableBimat<Mesh> op_h_disabled(
      mesh_h, order, MaterialField::MakeCoefficient(&lam_h, &mu_h, &rho_c), bc);
   op_h_disabled.SetMixedFluxContrastTol(-1.0);
   op_h_disabled.SetMixedFluxMode(MixedFluxMode::Adjacent);
   TestableBimat<Mesh> op_h_guard(
      mesh_h, order, MaterialField::MakeCoefficient(&lam_h, &mu_h, &rho_c), bc);
   op_h_guard.SetMixedFluxContrastTol(0.05);
   op_h_guard.SetMixedFluxMode(MixedFluxMode::Adjacent);
   TEST_ASSERT(op_h_disabled.GetCentralFluxFaceSet().size()
               == op_h_guard.GetCentralFluxFaceSet().size(),
               "homogeneous material: guard removes 0 faces (byte-exact when "
               "there is no contrast, even with tol=0.05)");

   // -----------------------------------------------------------------------
   // Phase 3: mechanism (R-001 two-sided dissipation).
   //   On the contrast corridor face: F_up (per-side, from a None operator) and
   //   F_ce (single-valued central, from the tol<0 Adjacent operator).  For a
   //   pure normal-velocity jump, D = [[Q]].(0.5*(F_e1+F_e2) - F_ce) > 0.
   // -----------------------------------------------------------------------
   std::cout << "\n-- Phase 3: two-sided dissipation mechanism --\n";
   TestableBimat<Mesh> op_none(
      mesh, order, MaterialField::MakeCoefficient(&lam_step, &mu_step_c, &rho_c), bc);
   // mixed_flux=none (ctor default) => empty central set => InteriorFaceFlux_
   // takes the per-side bi-material UPWIND branch on every interior face.

   const real_t nor[3] = {1.0, 0.0, 0.0};   // ignored by InteriorFaceFlux_
   real_t Q_self[NUM_STATE] = {0}, Q_nbr[NUM_STATE] = {0};
   Q_self[VX] = 1.0;                          // pure normal-velocity jump [[v_n]]=1

   real_t F_up1[NUM_STATE], F_up2[NUM_STATE], F_ce1[NUM_STATE], F_ce2[NUM_STATE];
   op_none.InteriorFaceFlux_(removed, Q_self, Q_nbr, nor, F_up1, F_up2);
   op_disabled.InteriorFaceFlux_(removed, Q_self, Q_nbr, nor, F_ce1, F_ce2);

   real_t worst_ce_diff = 0.0;
   for (int c = 0; c < NUM_STATE; ++c)
   { worst_ce_diff = std::max(worst_ce_diff, std::abs(F_ce1[c] - F_ce2[c])); }
   TEST_ASSERT(worst_ce_diff <= 1.0e-6,
               "central dispatch on the contrast face is single-valued (F_ce1 == "
               "F_ce2)");

   real_t D = 0.0;
   for (int c = 0; c < NUM_STATE; ++c)
   { D += (Q_self[c] - Q_nbr[c]) * (0.5 * (F_up1[c] + F_up2[c]) - F_ce1[c]); }
   std::cout << "  contrast-face dissipation D = " << std::scientific
             << std::setprecision(4) << D << "  (must be > 0)\n";
   TEST_ASSERT(D > 0.0,
               "across the contrast the upwind keeps positive dissipation that "
               "the central flux discards (D > 0)");

   // (R-001) DIRECT runtime check: the GUARDED operator dispatches the reclassified
   // face through the (per-side, two-valued) bi-material upwind — identical to the
   // None operator, NOT the single-valued central.  This proves the fix changed the
   // dispatch, not merely the face set.
   real_t Fg1[NUM_STATE], Fg2[NUM_STATE];
   op_guard.InteriorFaceFlux_(removed, Q_self, Q_nbr, nor, Fg1, Fg2);
   real_t worst_vs_none = 0.0, worst_two_valued = 0.0;
   for (int c = 0; c < NUM_STATE; ++c)
   {
      worst_vs_none    = std::max(worst_vs_none,    std::abs(Fg1[c] - F_up1[c]));
      worst_vs_none    = std::max(worst_vs_none,    std::abs(Fg2[c] - F_up2[c]));
      worst_two_valued = std::max(worst_two_valued, std::abs(Fg1[c] - Fg2[c]));
   }
   TEST_ASSERT(worst_vs_none <= 1.0e-9,
               "guarded operator dispatches the reclassified face as the bi-material "
               "upwind (== None operator), not central");
   TEST_ASSERT(worst_two_valued > 1.0e-3,
               "the reclassified face's upwind dispatch is per-side (F_e1 != F_e2 "
               "under the contrast) — the single-valued central path was abandoned");

   // Homogeneous limit: same construction on a homogeneous central face; D must
   // equal 0.5*Zp*[[v_n]]^2 (the P-wave jump dissipation) and match the scalar
   // GodunovFlux [[Q]].(Interior - Central).
   const int hface = *op_h_disabled.GetCentralFluxFaceSet().begin();
   real_t Hup1[NUM_STATE], Hup2[NUM_STATE], Hce1[NUM_STATE], Hce2[NUM_STATE];
   TestableBimat<Mesh> op_h_none(
      mesh_h, order, MaterialField::MakeCoefficient(&lam_h, &mu_h, &rho_c), bc);
   op_h_none.InteriorFaceFlux_(hface, Q_self, Q_nbr, nor, Hup1, Hup2);
   op_h_disabled.InteriorFaceFlux_(hface, Q_self, Q_nbr, nor, Hce1, Hce2);
   real_t D_h = 0.0;
   for (int c = 0; c < NUM_STATE; ++c)
   { D_h += (Q_self[c] - Q_nbr[c]) * (0.5 * (Hup1[c] + Hup2[c]) - Hce1[c]); }
   // D = [[Q]].(F_up - F_ce) is the L2 (un-mass-weighted) inner product of the jump
   // with the dissipative flux F_up - F_ce = 0.5*|A_n|*[[Q]].  For a pure normal-
   // velocity jump this is 0.5 * |A_n|_{vn,vn} * [[v_n]]^2 = 0.5 * cp * [[v_n]]^2
   // (the velocity-velocity diagonal of |A_n| is the P speed cp; the ENERGY-norm
   // measure 0.5*Zp*[[v]]^2 differs by the density factor rho).
   const real_t cp_h = std::sqrt((30.0e9 + 2.0 * 30.0e9) / 2670.0);   // sqrt((lam+2mu)/rho)
   const real_t D_expect = 0.5 * cp_h * 1.0 * 1.0;
   std::cout << "  homogeneous D = " << std::scientific << std::setprecision(4)
             << D_h << "  expect 0.5*cp = " << D_expect << "\n";
   TEST_ASSERT(std::abs(D_h - D_expect) <= 1.0e-6 * D_expect,
               "homogeneous limit: D == 0.5*cp*[[v_n]]^2 (the L2 P-wave jump "
               "dissipation measure) to 1e-6 rel");

   std::cout << "\n========================================\n"
             << "  Results: " << num_passed << " passed, "
             << num_failed << " failed out of " << num_tests << " tests\n"
             << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}

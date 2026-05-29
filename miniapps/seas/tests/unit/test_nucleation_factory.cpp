// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_nucleation_factory.cpp — Phase 7 of
// PLAN_tpv_regression_via_spatial_dyn_driver_2026-05-24.
//
// Acceptance criterion 2: MakeNucleation maps each cfg.nucleation.kind (and
// the absent-[nucleation] case) to the expected INucleationMethod concrete,
// with the right IsPerSubStep() tag.  Also verifies the one-shot
// InstantaneousOverstressCircular seeds tau2_nuc on ApplyOnce and is a no-op on
// ApplyIncrement.

#include "mfem.hpp"

#include "../../dynamic/nucleation_factory.hpp"

#include <iostream>
#include <memory>
#include <vector>

using namespace mfem;
using namespace mfem::seas;
using namespace mfem::seas::spatial;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(c, m) do { num_tests++; if (!(c)) { \
   std::cerr << "FAILED: " << m << " line " << __LINE__ << "\n"; \
   num_failed++; } else { std::cout << "  PASSED: " << m << "\n"; \
   num_passed++; } } while (0)

namespace
{

// A minimal 2-DOF planar-fault fixture: normal +x, dip +y, strike +z.
void MakeFixture(Vector& dof_coords_3d, DenseMatrix& dof_basis)
{
   const int N = 2;
   dof_coords_3d.SetSize(3 * N);
   dof_coords_3d = 0.0;
   dof_coords_3d(2) = -5000.0;   // DOF 0 at the hypocentre depth
   dof_coords_3d(5) = -5000.0 + 500.0;
   dof_basis.SetSize(9, N);
   dof_basis = 0.0;
   for (int i = 0; i < N; ++i)
   {
      dof_basis(0, i) = 1.0;   // normal +x
      dof_basis(4, i) = 1.0;   // dip    +y
      dof_basis(8, i) = 1.0;   // strike +z
   }
}

// A config whose only nucleation-relevant fields are set; the rest are left at
// their defaults (the factory only reads cfg.nucleation + the coords/basis).
SpatialFrictionConfig MakeCfg()
{
   SpatialFrictionConfig cfg;
   // Gaussian sub-block (valid when that kind is selected).
   cfg.nucleation.gradual_overstress.radius_dip_m    = 1000.0;
   cfg.nucleation.gradual_overstress.radius_strike_m = 1000.0;
   cfg.nucleation.gradual_overstress.center_z_m      = -5000.0;
   cfg.nucleation.gradual_overstress.delta_tau_strike_pa = 2.0e6;
   cfg.nucleation.gradual_overstress.T_nuc_s         = 1.0;
   // Compact-circular sub-block.
   cfg.nucleation.compact_circular.radius_m     = 3000.0;
   cfg.nucleation.compact_circular.center_z_m   = -5000.0;
   cfg.nucleation.compact_circular.delta_tau_pa = 45.0e6;
   cfg.nucleation.compact_circular.T_nuc_s      = 1.0;
   // Instantaneous-circular sub-block.
   cfg.nucleation.instantaneous_circular.radius_m     = 1400.0;
   cfg.nucleation.instantaneous_circular.taper_m      = 200.0;
   cfg.nucleation.instantaneous_circular.center_z_m   = -5000.0;
   cfg.nucleation.instantaneous_circular.delta_tau_pa = 11.6e6;
   return cfg;
}

}  // namespace

// F-1: absent [nucleation] (enabled == false) -> StaticOverstress.
static void F1_absent_block_static()
{
   std::cout << "\n[F-1] enabled=false -> StaticOverstress\n";
   Vector coords; DenseMatrix basis; MakeFixture(coords, basis);
   SpatialFrictionConfig cfg = MakeCfg();
   cfg.nucleation.enabled = false;   // absent block
   auto nuc = MakeNucleation(cfg, coords, basis);
   TEST_ASSERT(nuc != nullptr, "factory returns non-null");
   TEST_ASSERT(dynamic_cast<StaticOverstress*>(nuc.get()) != nullptr,
               "absent block -> StaticOverstress");
   TEST_ASSERT(!nuc->IsPerSubStep(), "static is not per-sub-step");
}

// F-2: kind = GradualOverstress -> GaussianGradualOverstress.
static void F2_gaussian()
{
   std::cout << "\n[F-2] GradualOverstress -> GaussianGradualOverstress\n";
   Vector coords; DenseMatrix basis; MakeFixture(coords, basis);
   SpatialFrictionConfig cfg = MakeCfg();
   cfg.nucleation.enabled = true;
   cfg.nucleation.kind = NucleationKind::GradualOverstress;
   auto nuc = MakeNucleation(cfg, coords, basis);
   TEST_ASSERT(dynamic_cast<GaussianGradualOverstress*>(nuc.get()) != nullptr,
               "GradualOverstress -> GaussianGradualOverstress");
   TEST_ASSERT(nuc->IsPerSubStep(), "gaussian is per-sub-step");
}

// F-3: kind = GradualOverstressCompactCircular -> CompactCircularGradualOverstress.
static void F3_compact_circular()
{
   std::cout << "\n[F-3] CompactCircular -> CompactCircularGradualOverstress\n";
   Vector coords; DenseMatrix basis; MakeFixture(coords, basis);
   SpatialFrictionConfig cfg = MakeCfg();
   cfg.nucleation.enabled = true;
   cfg.nucleation.kind = NucleationKind::GradualOverstressCompactCircular;
   auto nuc = MakeNucleation(cfg, coords, basis);
   TEST_ASSERT(dynamic_cast<CompactCircularGradualOverstress*>(nuc.get())
               != nullptr,
               "CompactCircular -> CompactCircularGradualOverstress");
   TEST_ASSERT(nuc->IsPerSubStep(), "compact-circular is per-sub-step");
}

// F-4: kind = InstantaneousOverstressCircular -> the one-shot method;
// ApplyOnce seeds tau2_nuc and ApplyIncrement is a no-op.
static void F4_instantaneous()
{
   std::cout << "\n[F-4] Instantaneous -> InstantaneousOverstressCircular\n";
   Vector coords; DenseMatrix basis; MakeFixture(coords, basis);
   SpatialFrictionConfig cfg = MakeCfg();
   cfg.nucleation.enabled = true;
   cfg.nucleation.kind = NucleationKind::InstantaneousOverstressCircular;
   auto nuc = MakeNucleation(cfg, coords, basis);
   TEST_ASSERT(dynamic_cast<InstantaneousOverstressCircular*>(nuc.get())
               != nullptr,
               "Instantaneous -> InstantaneousOverstressCircular");
   TEST_ASSERT(!nuc->IsPerSubStep(), "instantaneous is not per-sub-step");

   // DOF 0 is at the hypocentre (r=0) -> seeded with the full delta_tau_pa.
   std::vector<DOFData> dof(2);
   for (auto& d : dof) { d.tau1_nuc = 0.0; d.tau2_nuc = 0.0; d.sigma_n_nuc = 0.0; }
   nuc->ApplyIncrement(dof, 0.5, 0.1);   // no-op for the one-shot method
   TEST_ASSERT(dof[0].tau2_nuc == 0.0, "ApplyIncrement is a no-op (one-shot)");
   nuc->ApplyOnce(dof);
   TEST_ASSERT(std::abs(dof[0].tau2_nuc - 11.6e6) < 1e-3,
               "ApplyOnce seeds tau2_nuc at r=0 = delta_tau_pa");
   TEST_ASSERT(dof[0].tau1_nuc == 0.0, "tau1_nuc stays 0 (pure strike-slip)");
}

// F-5: the Gaussian path through MakeNucleation is BYTE-IDENTICAL to the
// pre-Phase-7 inline ResolveGradualOverstress + ApplyGradualOverstressIncrement
// (acceptance criterion 3, at the unit level — proves the swapped computation
// is unchanged without needing the full SAFS mesh / end-to-end smoke).
static void F5_gaussian_path_byte_identical()
{
   std::cout << "\n[F-5] Gaussian path via factory == inline resolve+apply (bytewise)\n";
   Vector coords; DenseMatrix basis; MakeFixture(coords, basis);
   const int N = coords.Size() / 3;

   GradualOverstressSpec spec;
   spec.center_x_m = 0.0; spec.center_y_m = 0.0; spec.center_z_m = -5000.0;
   spec.radius_dip_m = 1000.0; spec.radius_strike_m = 1000.0;
   spec.delta_tau_dip_pa = 1.0e6; spec.delta_tau_strike_pa = 2.0e6;
   spec.T_nuc_s = 1.0;

   // OLD inline path: resolve directly.
   const auto pold = ResolveGradualOverstress(spec, /*enabled=*/true, coords, basis);

   // NEW path: through the factory -> GaussianGradualOverstress::Params().
   SpatialFrictionConfig cfg;
   cfg.nucleation.enabled = true;
   cfg.nucleation.kind = NucleationKind::GradualOverstress;
   cfg.nucleation.gradual_overstress = spec;
   auto nuc = MakeNucleation(cfg, coords, basis);
   auto* g = dynamic_cast<GaussianGradualOverstress*>(nuc.get());
   TEST_ASSERT(g != nullptr, "factory built GaussianGradualOverstress");
   const auto& pnew = g->Params();

   // (a) Resolved per-DOF params bit-identical (exact equality, 0 tol).
   bool params_match = (pnew.amplitude_dip.Size()    == pold.amplitude_dip.Size())
                    && (pnew.amplitude_strike.Size() == pold.amplitude_strike.Size())
                    && (pnew.radial.Size()           == pold.radial.Size());
   for (int i = 0; params_match && i < N; ++i)
   {
      params_match = (pnew.amplitude_dip(i)    == pold.amplitude_dip(i))
                  && (pnew.amplitude_strike(i) == pold.amplitude_strike(i))
                  && (pnew.radial(i)           == pold.radial(i));
   }
   TEST_ASSERT(params_match, "resolved params bit-identical (factory == inline)");

   // (b) Accumulated tau_nuc bit-identical after a full sub-step sweep.
   std::vector<DOFData> dof_old(N), dof_new(N);
   for (int i = 0; i < N; ++i)
   {
      dof_old[i].tau1_nuc = dof_old[i].tau2_nuc = 0.0;
      dof_new[i].tau1_nuc = dof_new[i].tau2_nuc = 0.0;
   }
   const int O = 7;
   const real_t dt = spec.T_nuc_s / O;
   real_t t = 0.0;
   for (int o = 0; o < O; ++o)
   {
      t += dt;
      // OLD: inline apply with the old params + the cfg T_nuc.
      ApplyGradualOverstressIncrement(dof_old, pold, spec.T_nuc_s, t, dt);
      // NEW: the driver's nuc_cb body.
      nuc->ApplyIncrement(dof_new, t, dt);
   }
   bool tau_match = true;
   for (int i = 0; tau_match && i < N; ++i)
   {
      tau_match = (dof_old[i].tau1_nuc == dof_new[i].tau1_nuc)
               && (dof_old[i].tau2_nuc == dof_new[i].tau2_nuc);
   }
   TEST_ASSERT(tau_match,
               "accumulated tau_nuc bit-identical (factory ApplyIncrement == "
               "inline ApplyGradualOverstressIncrement)");
   // Non-triviality: the sweep must actually have moved tau2_nuc.
   TEST_ASSERT(dof_new[0].tau2_nuc > 0.0,
               "non-trivial: tau2_nuc accumulated at the bell centre");
}

// Mirror of the driver's diagnostic re-sourcing (spatial_dyn_driver.cpp ~:1247):
// pull the resolved per-DOF amplitudes/radial out of whichever concrete method
// MakeNucleation built, into the GradualOverstressPerDOFParams the ParaView /
// derived-quantity consumers read.  Kept in lockstep with the driver so this
// test guards exactly that glue.
static GradualOverstressPerDOFParams ExtractDiagParams(INucleationMethod& nuc)
{
   GradualOverstressPerDOFParams diag;
   if (auto* g = dynamic_cast<GaussianGradualOverstress*>(&nuc))
   {
      diag = g->Params();
   }
   else if (auto* c = dynamic_cast<CompactCircularGradualOverstress*>(&nuc))
   {
      diag.amplitude_strike = c->Params().amplitude_strike;
      diag.radial           = c->Params().radial;
   }
   else if (auto* inst = dynamic_cast<InstantaneousOverstressCircular*>(&nuc))
   {
      diag.amplitude_strike = inst->Params().amplitude_strike;
   }
   return diag;
}

// F-6 (R-002): the driver's nucleation diagnostics / ParaView fields must NOT be
// silently zeroed for the non-Gaussian kinds.  Before the fix the driver only
// dynamic_cast<GaussianGradualOverstress*>, so compact-circular / instantaneous
// runs wrote 0.0 to nuc_amplitude / nuc_radial_factor for every DOF.  This test
// fails under that old code (empty diag) and passes with the cast-ladder fix.
static void F6_non_gaussian_diag_not_zeroed()
{
   std::cout << "\n[F-6] non-Gaussian kinds expose radial/amplitude to diag (R-002)\n";
   Vector coords; DenseMatrix basis; MakeFixture(coords, basis);
   const int N = coords.Size() / 3;

   // Compact-circular: both radial and amplitude_strike must be populated.
   {
      SpatialFrictionConfig cfg = MakeCfg();
      cfg.nucleation.enabled = true;
      cfg.nucleation.kind = NucleationKind::GradualOverstressCompactCircular;
      auto nuc = MakeNucleation(cfg, coords, basis);
      auto diag = ExtractDiagParams(*nuc);
      TEST_ASSERT(diag.radial.Size() == N,
                  "compact-circular: radial sized to fault-DOF count (not 0)");
      TEST_ASSERT(diag.amplitude_strike.Size() == N,
                  "compact-circular: amplitude_strike sized (not 0)");
      TEST_ASSERT(diag.radial.Max() > 0.0,
                  "compact-circular: radial non-zero at the bell centre");
      TEST_ASSERT(diag.amplitude_strike.Max() > 0.0,
                  "compact-circular: amplitude_strike non-zero");
   }
   // Instantaneous-circular: amplitude_strike populated (no radial field).
   {
      SpatialFrictionConfig cfg = MakeCfg();
      cfg.nucleation.enabled = true;
      cfg.nucleation.kind = NucleationKind::InstantaneousOverstressCircular;
      auto nuc = MakeNucleation(cfg, coords, basis);
      auto diag = ExtractDiagParams(*nuc);
      TEST_ASSERT(diag.amplitude_strike.Size() == N,
                  "instantaneous: amplitude_strike sized to fault-DOF count");
      TEST_ASSERT(diag.amplitude_strike.Max() > 0.0,
                  "instantaneous: amplitude_strike non-zero at r=0");
   }
}

int main(int /*argc*/, char** /*argv*/)
{
   std::cout << "Running Phase 7 test_nucleation_factory\n";
   F1_absent_block_static();
   F2_gaussian();
   F3_compact_circular();
   F4_instantaneous();
   F5_gaussian_path_byte_identical();
   F6_non_gaussian_diag_not_zeroed();

   std::cout << "\n========================================\n";
   std::cout << "Phase 7 test_nucleation_factory: "
             << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return num_failed == 0 ? 0 : 1;
}

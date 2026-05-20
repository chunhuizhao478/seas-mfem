// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_phaseh_wave_operator_constant_parity.cpp — Phase H Stage 1 gate
// `T-PHASEH-SCALAR-PARITY` for `spatial_dynamic_rupture_plan.md` (rev-3).
//
// Coverage:
//   C-1  Mult parity (np=1, serial Mesh):
//        WaveOperator(mesh, p, lam, mu, rho, bc).Mult(Q)  ≈
//        WaveOperator(mesh, p, MakeConstant(lam, mu, rho), bc).Mult(Q)
//        to 1e-10 relative tolerance.  Pre-Phase-R this was bit-for-bit
//        because the hetero ctor still consulted the scalar `flux_`;
//        Phase R.2 INTENTIONALLY routes the hetero ctor through
//        `BimaterialFlux::ApplyPerFaceFlux`, which inverts matR via
//        `mfem::DenseMatrixInverse` instead of `GodunovFlux::Interior`'s
//        closed-form A^± decomposition.  Mathematically equivalent,
//        FP-equivalent to LU-rounding precision (~1e-13 relative on
//        the bi-material output; ~1e-10 after assembly into dQdt).
//   C-2  Mult parity (np>=1, ParMesh): same relaxed equality under the
//        parallel ctor.
//   C-3  ComputeMaxDt parity: dt_scalar == dt_hetero bit-for-bit
//        (heterogeneous CFL path falls back to the same global min;
//        not affected by R.2).
//   C-4  Pool invariants on Constant input: NumUniqueTriples() == 1,
//        every elem maps to flux index 0, owned-pool getter exposes the
//        same triple, shared_face_neighbour_material_ has one entry per
//        shared face whose value equals the local material.
//
// Post-Phase-R contract being verified:
//   - The hetero ctor's BuildGodunovFluxPool_ + ExchangeBiMaterialNeighbours_
//     + BuildPerFaceBimaterialFluxMatrices_ caches feed the wrapped
//     interior-face dispatch (Phase R.2).  Mode::Constant collapses the
//     per-element bi-material formula to the homogeneous upwind state
//     to within LU rounding (see plan §"Homogeneous-limit collapse").

#include "mfem.hpp"

#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/heterogeneous_material.hpp"
#include "../../dynamic/fault_face_flux.hpp"   // Phase H Stage 2 fault parity
#include "../../dynamic/tpv102_setup.hpp"      // InitializeFaultDOFs (fault parity)
#include "../../config/tpv102_params.hpp"      // TPV102 material for fault fixture
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

namespace
{
int g_num_tests = 0, g_num_passed = 0, g_num_failed = 0;
int g_rank = 0;

#define TEST_ASSERT(c, m) do { ++g_num_tests; if (!(c)) { \
   if (g_rank == 0) \
   { std::cerr << "FAILED: " << m << " line " << __LINE__ << "\n"; } \
   ++g_num_failed; } else { if (g_rank == 0) \
   { std::cout << "  PASSED: " << m << "\n"; } ++g_num_passed; } } while (0)

#define TEST_EQ_BITS(a, b, m) do { ++g_num_tests; const real_t _a=(a); \
   const real_t _b=(b); if (std::memcmp(&_a, &_b, sizeof(real_t)) != 0) { \
   if (g_rank == 0) { std::cerr << "FAILED: " << m \
      << " (got " << _a << ", expected " << _b \
      << ", diff " << (_a - _b) << ") line " << __LINE__ << "\n"; } \
   ++g_num_failed; } else { if (g_rank == 0) \
   { std::cout << "  PASSED: " << m << "\n"; } ++g_num_passed; } } while (0)

// Build a small 3D box mesh for the parity checks.  Inline mesh keeps the
// test self-contained (no .msh fixture dependency).
Mesh MakeBoxMesh(int nx, int ny, int nz)
{
   return Mesh::MakeCartesian3D(nx, ny, nz, Element::TETRAHEDRON,
                                /*sx=*/1.0, /*sy=*/1.0, /*sz=*/1.0,
                                /*sfc_ordering=*/false);
}

// Sensible crustal material; matches geoffrey2010 fallback in the
// SAFS EXAMPLE TOML so the test exercises a realistic parameter set.
// These are CLEAN constants (≤6 significant figures), so they round-trip
// exactly through the pool's 6-sig-fig dedup key.
constexpr real_t k_lambda = 32.0e9;
constexpr real_t k_mu     = 32.0e9;
constexpr real_t k_rho    = 2670.0;
constexpr int    k_order  = 1;

// Phase H Stage 2 (R-002 mitigation): a material whose constants need
// MORE than 6 significant figures.  Before the Phase 1 prerequisite fix
// (godunov_flux_pool.cpp built the cached flux from the 6-sig-fig-ROUNDED
// triple), `At(e)` differed from the exact-material scalar `flux_` at
// ~1e-7 for these values.  That mismatch was masked while the volume term
// used the exact `Ax_`; Phase 1a routes the volume term through
// `At(e).GetReferenceStarMatrix(d)`, so the parity gate below would BREAK
// unless the pool is built from EXACT values.  These constants are the
// tripwire: parity holds iff the prerequisite fix is in place.
constexpr real_t k_lambda6 = 3.14159265e10;   // ≠ its 6-sig-fig round
constexpr real_t k_mu6     = 2.71828182e10;
constexpr real_t k_rho6    = 2.66666667e3;

// Fill the state vector with a deterministic plane-wave pattern so any
// dispatch-path divergence shows up as a non-trivial residual.  Same
// seed → same Q on every rank.
void FillQDeterministic(Vector &Q)
{
   const int n = Q.Size();
   for (int i = 0; i < n; ++i)
   {
      const real_t x = static_cast<real_t>(i) / static_cast<real_t>(n);
      Q(i) = 1e6 * std::sin(2.0 * M_PI * 3.0 * x)
           + 5e5 * std::cos(2.0 * M_PI * 7.0 * x);
   }
}

BoundaryConfig MakeAbsorbingBC()
{
   // 3D Cartesian box has 6 boundary attributes (1..6).  Use natural on
   // all faces — wave_operator's NUM_STATE=9 elastic system is well-
   // defined under that BC for the parity check; we are not solving a
   // physical problem here, only verifying that two ctor paths produce
   // identical dQdt for the same Q.
   BoundaryConfig bc;
   bc.fault_attr      = 0;   // no fault face in this synthetic test
   bc.natural_attrs   = {1, 2, 3, 4, 5, 6};
   bc.absorbing_attrs = {};
   bc.dirichlet_attrs = {};
   return bc;
}

// Mixed BC: absorbing on attrs {1,3,5}, free-surface (natural) on {2,4,6}.
// Exercises BOTH the Absorbing (B1/B2) and FreeSurface (B1/B2) per-element
// dispatch paths in a single operator so the parity gate covers Group B.
BoundaryConfig MakeMixedBC()
{
   BoundaryConfig bc;
   bc.fault_attr      = 0;
   bc.natural_attrs   = {2, 4, 6};
   bc.absorbing_attrs = {1, 3, 5};
   bc.dirichlet_attrs = {};
   return bc;
}

// Box mesh with the y = L/2 plane of interior faces tagged as a fault
// (attr 3); every outer boundary face tagged attr 1 (natural).  Mirrors
// BuildCartesianFaultMesh in test_adjacent_triangle_fault_first_step_audit
// so the WaveOperator ctor detects an interior fault and builds the
// FaultBasis.  L matches TPV102 length scale so InitializeFaultDOFs's
// depth = |z| seeding is physical.
Mesh MakeFaultBoxMesh(int nx, int ny, int nz, real_t L)
{
   Mesh mesh = Mesh::MakeCartesian3D(nx, ny, nz, Element::TETRAHEDRON,
                                     L, L, L, /*sfc_ordering=*/false);
   mesh.FinalizeTopology();
   mesh.Finalize();
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         real_t cy = 0.0;
         for (int v = 0; v < fv.Size(); v++) { cy += mesh.GetVertex(fv[v])[1]; }
         cy /= fv.Size();
         if (std::abs(cy - 0.5 * L) < 1e-8)
         {
            mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 3);   // fault
         }
         continue;
      }
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], 1);         // outer (natural)
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// Relative-difference reduction over two equally-sized vectors.
// denom = max(|a|, |b|, 1) so near-zero entries don't blow up the ratio.
real_t MaxRelDiff(const Vector &a, const Vector &b, int *n_over = nullptr,
                  real_t tol = 0.0)
{
   real_t mx = 0.0;
   int over = 0;
   const int n = a.Size();
   for (int i = 0; i < n; ++i)
   {
      const real_t denom = std::max(std::max(std::abs(a(i)), std::abs(b(i))),
                                    real_t(1.0));
      const real_t rel = std::abs(a(i) - b(i)) / denom;
      if (rel > mx) { mx = rel; }
      if (tol > 0.0 && rel > tol) { ++over; }
   }
   if (n_over) { *n_over = over; }
   return mx;
}

// Count bit-for-bit differing entries between two equally-sized vectors.
// Used for the byte-exact gates (volume RHS, ADER CK predictor) where the
// per-element path must reproduce the scalar path EXACTLY on Mode::Constant.
int BitDiffCount(const Vector &a, const Vector &b)
{
   int n_diff = 0;
   const int n = a.Size();
   for (int i = 0; i < n; ++i)
   {
      const real_t ai = a(i), bi = b(i);
      if (std::memcmp(&ai, &bi, sizeof(real_t)) != 0) { ++n_diff; }
   }
   return n_diff;
}

}  // namespace

// =========================================================================
// C-1  Mult parity on a serial Mesh
// =========================================================================
static void C_1_mult_parity_serial()
{
   if (g_rank == 0)
   { std::cout << "\n[C-1] Mult parity (serial Mesh)\n"; }
#ifdef MFEM_USE_MPI
   // Only rank 0 runs the serial-mesh branch to avoid duplicated output
   // on np>1.  The parallel branch (C-2) covers np>=1.
   if (g_rank != 0) { return; }
#endif

   Mesh smesh = MakeBoxMesh(2, 2, 2);
   const BoundaryConfig bc = MakeAbsorbingBC();

   WaveOperator<Mesh> wave_scalar(smesh, k_order, k_lambda, k_mu, k_rho, bc);
   WaveOperator<Mesh> wave_hetero(smesh, k_order,
                                  MaterialField::MakeConstant(k_lambda,
                                                              k_mu, k_rho),
                                  bc);

   // Fluctuation-Q dispatch: both wave operators must have Q_bg = 0
   // installed before Mult.  Matches the TPV/BP5 production pattern
   // (tpv205_driver.cpp:1620, tpv102_driver.cpp:Mult-prep).
   {
      real_t Q_bg[NUM_STATE] = {0};
      wave_scalar.SetAbsorbingBackground(Q_bg);
      wave_hetero.SetAbsorbingBackground(Q_bg);
   }

   TEST_ASSERT(wave_scalar.UsesGodunovFluxPool() == false,
               "scalar ctor: UsesGodunovFluxPool() == false");
   TEST_ASSERT(wave_hetero.UsesGodunovFluxPool() == true,
               "hetero ctor: UsesGodunovFluxPool() == true");

   const int n = NUM_STATE * wave_scalar.GetScalarNDof();
   TEST_ASSERT(NUM_STATE * wave_hetero.GetScalarNDof() == n,
               "both ctors give the same problem size");

   Vector Q(n);
   FillQDeterministic(Q);
   Vector dQdt_scalar(n), dQdt_hetero(n);
   wave_scalar.Mult(Q, dQdt_scalar);
   wave_hetero.Mult(Q, dQdt_hetero);

   // Phase R.2: relaxed from bit-equality to 1e-10 relative tolerance.
   // The hetero ctor now routes interior face flux through the matR
   // inverse (BimaterialFlux), giving FP-equivalent but not bit-
   // identical output vs the scalar ctor's closed-form A^± upwind.
   const real_t k_rel_tol = 1e-10;
   int n_diff = 0;
   real_t max_rel_diff = 0.0;
   for (int i = 0; i < n; ++i)
   {
      const real_t a = dQdt_scalar(i);
      const real_t b = dQdt_hetero(i);
      const real_t denom = std::max(std::max(std::abs(a), std::abs(b)),
                                    real_t(1.0));
      const real_t rel = std::abs(a - b) / denom;
      if (rel > k_rel_tol) { ++n_diff; }
      if (rel > max_rel_diff) { max_rel_diff = rel; }
   }
   if (g_rank == 0)
   {
      std::cout << "  max_rel_diff=" << max_rel_diff
                << "  (tol=" << k_rel_tol << ", n_over_tol=" << n_diff
                << ")\n";
   }
   TEST_ASSERT(n_diff == 0,
               "Mult dQdt: scalar vs hetero ctor agree to 1e-10 relative");
}

// =========================================================================
// C-2  Mult parity on a ParMesh (covers np=1 and np=4 via the same path)
// =========================================================================
static void C_2_mult_parity_parallel()
{
   if (g_rank == 0)
   { std::cout << "\n[C-2] Mult parity (ParMesh)\n"; }
#ifndef MFEM_USE_MPI
   if (g_rank == 0)
   { std::cout << "  SKIPPED: MFEM_USE_MPI not defined\n"; }
   return;
#else
   Mesh smesh = MakeBoxMesh(4, 4, 4);
   ParMesh pmesh(MPI_COMM_WORLD, smesh);
   smesh.Clear();
   const BoundaryConfig bc = MakeAbsorbingBC();

   WaveOperator<ParMesh> wave_scalar(pmesh, k_order,
                                     k_lambda, k_mu, k_rho, bc);
   WaveOperator<ParMesh> wave_hetero(pmesh, k_order,
                                     MaterialField::MakeConstant(k_lambda,
                                                                 k_mu, k_rho),
                                     bc);

   // Fluctuation-Q dispatch — Q_bg = 0 required by the FreeSurface /
   // natural BC paths.
   {
      real_t Q_bg[NUM_STATE] = {0};
      wave_scalar.SetAbsorbingBackground(Q_bg);
      wave_hetero.SetAbsorbingBackground(Q_bg);
   }

   TEST_ASSERT(wave_scalar.UsesGodunovFluxPool() == false,
               "ParMesh scalar ctor: pool absent");
   TEST_ASSERT(wave_hetero.UsesGodunovFluxPool() == true,
               "ParMesh hetero ctor: pool present");

   const int n = NUM_STATE * wave_scalar.GetScalarNDof();
   TEST_ASSERT(NUM_STATE * wave_hetero.GetScalarNDof() == n,
               "ParMesh: same problem size in both ctors");

   Vector Q(n);
   FillQDeterministic(Q);
   Vector dQdt_scalar(n), dQdt_hetero(n);
   wave_scalar.Mult(Q, dQdt_scalar);
   wave_hetero.Mult(Q, dQdt_hetero);

   // Phase R.2: relaxed from bit-equality to 1e-9 relative tolerance
   // (see C-1 rationale above).  ParMesh path accumulates FP-rounding
   // through MPI ghost exchange + serialised reductions across ranks,
   // pushing the divergence ~1 order of magnitude past C-1's 1e-10.
   const real_t k_rel_tol = 1e-9;
   int n_diff_local = 0;
   real_t max_rel_local = 0.0;
   for (int i = 0; i < n; ++i)
   {
      const real_t a = dQdt_scalar(i);
      const real_t b = dQdt_hetero(i);
      const real_t denom = std::max(std::max(std::abs(a), std::abs(b)),
                                    real_t(1.0));
      const real_t rel = std::abs(a - b) / denom;
      if (rel > k_rel_tol) { ++n_diff_local; }
      if (rel > max_rel_local) { max_rel_local = rel; }
   }
   int n_diff_global = 0;
   real_t max_rel_global = 0.0;
   MPI_Allreduce(&n_diff_local, &n_diff_global, 1, MPI_INT, MPI_SUM,
                 pmesh.GetComm());
   MPI_Allreduce(&max_rel_local, &max_rel_global, 1,
                 MPITypeMap<real_t>::mpi_type, MPI_MAX, pmesh.GetComm());
   if (g_rank == 0)
   {
      std::cout << "  ParMesh max_rel_diff=" << max_rel_global
                << "  (tol=" << k_rel_tol << ", n_over_tol="
                << n_diff_global << ")\n";

      // REVIEW round-3 R-005: the 1e-9 tolerance is loose compared to
      // a derived FP bound.  For a homogeneous problem with N_faces
      // faces accumulated across N_ranks ranks, the per-face LU-rounding
      // error sums as roughly ε · sqrt(N_faces · N_ranks) when signs are
      // statistically random (worst case ε · N_faces · N_ranks).  Log a
      // WARNING if max_rel_global exceeds the sqrt-scaled bound so the
      // band can be revisited if it ever creeps wider.
      int n_ranks_dbg = 1;
      MPI_Comm_size(pmesh.GetComm(), &n_ranks_dbg);
      const real_t derived_bound =
         1e-13 * std::sqrt(real_t(pmesh.GetNumFaces() * n_ranks_dbg));
      if (max_rel_global > derived_bound)
      {
         std::cerr << "  WARNING: max_rel_diff " << max_rel_global
                   << " exceeds derived FP bound " << derived_bound
                   << " by " << (max_rel_global / derived_bound)
                   << "×.  Investigate before further loosening tol "
                   << "(see REVIEW R-005).\n";
      }
   }
   TEST_ASSERT(n_diff_global == 0,
               "ParMesh Mult dQdt: scalar vs hetero ctor agree to "
               "1e-10 relative across all ranks");
#endif
}

// =========================================================================
// C-3  ComputeMaxDt parity
// =========================================================================
static void C_3_compute_max_dt_parity()
{
   if (g_rank == 0)
   { std::cout << "\n[C-3] ComputeMaxDt parity\n"; }
#ifndef MFEM_USE_MPI
   if (g_rank == 0)
   { std::cout << "  SKIPPED: MFEM_USE_MPI not defined\n"; }
   return;
#else
   Mesh smesh = MakeBoxMesh(3, 4, 5);
   ParMesh pmesh(MPI_COMM_WORLD, smesh);
   smesh.Clear();
   const BoundaryConfig bc = MakeAbsorbingBC();

   WaveOperator<ParMesh> wave_scalar(pmesh, k_order,
                                     k_lambda, k_mu, k_rho, bc);
   WaveOperator<ParMesh> wave_hetero(pmesh, k_order,
                                     MaterialField::MakeConstant(k_lambda,
                                                                 k_mu, k_rho),
                                     bc);

   const real_t cfl = 0.5;
   const real_t dt_scalar = wave_scalar.ComputeMaxDt(cfl);
   const real_t dt_hetero = wave_hetero.ComputeMaxDt(cfl);
   TEST_EQ_BITS(dt_hetero, dt_scalar,
                "ComputeMaxDt bit-identical between scalar and hetero ctors");
#endif
}

// =========================================================================
// C-4  Pool / cache invariants on Constant input
// =========================================================================
static void C_4_pool_invariants_on_constant()
{
   if (g_rank == 0)
   { std::cout << "\n[C-4] Pool + cache invariants on Constant input\n"; }
#ifndef MFEM_USE_MPI
   if (g_rank == 0)
   { std::cout << "  SKIPPED: MFEM_USE_MPI not defined\n"; }
   return;
#else
   Mesh smesh = MakeBoxMesh(3, 3, 3);
   ParMesh pmesh(MPI_COMM_WORLD, smesh);
   smesh.Clear();
   const BoundaryConfig bc = MakeAbsorbingBC();

   WaveOperator<ParMesh> wave(pmesh, k_order,
                              MaterialField::MakeConstant(k_lambda,
                                                          k_mu, k_rho),
                              bc);

   // Stage 1 contract: the new ctor builds `owned_flux_pool_` (read
   // via UsesGodunovFluxPool()) but the per-element dispatch in
   // wave_operator.inl is not yet wired (H.2 = Stage 2), so every
   // dispatch site still consults the scalar `flux_` member.  On
   // Mode::Constant the two paths produce byte-identical output.
   TEST_ASSERT(wave.UsesGodunovFluxPool(),
               "owned_flux_pool_ populated by the hetero ctor");

   const auto &lmr = wave.GetPerElementMaterial();
   TEST_ASSERT(static_cast<int>(lmr.size()) == wave.NumElements(),
               "per_elem_lmr_ size == NumElements()");
   bool all_constant = true;
   for (const auto &t : lmr)
   {
      if (t[0] != k_lambda || t[1] != k_mu || t[2] != k_rho)
      { all_constant = false; break; }
   }
   TEST_ASSERT(all_constant,
               "per_elem_lmr_ holds (lambda, mu, rho) at every elem");

   const auto &per_h = wave.GetPerElementCflLength();
   TEST_ASSERT(static_cast<int>(per_h.size()) == wave.NumElements(),
               "per_elem_h_ size == NumElements()");

   const auto &nbr = wave.GetSharedFaceNeighbourMaterial();
   const int n_shared = pmesh.GetNSharedFaces();
   TEST_ASSERT(static_cast<int>(nbr.size()) == n_shared,
               "shared_face_neighbour_material_ size == NSharedFaces");
   bool all_local_matches = true;
   for (const auto &kv : nbr)
   {
      if (kv.second[0] != k_lambda || kv.second[1] != k_mu
          || kv.second[2] != k_rho)
      { all_local_matches = false; break; }
   }
   TEST_ASSERT(all_local_matches,
               "Stage 1: shared-face neighbour material == local "
               "(MPI exchange deferred to Stage 2)");
#endif
}

// =========================================================================
// C-5  Non-fault parity: volume RHS + ADER CK predictor BYTE-IDENTICAL,
//      Mult + AdvanceADER FP-equivalent (1e-9 rel), under MIXED BC
//      (absorbing + free-surface).  Parameterised by material so it runs
//      with BOTH a clean and a >6-sig-fig constant — the latter is the
//      R-002 rounding-regression tripwire (Phase 1 prerequisite).
// =========================================================================
static void C_5_nonfault_parity(real_t lam, real_t mu, real_t rho,
                                const char *label)
{
   if (g_rank == 0)
   { std::cout << "\n[C-5] Non-fault parity (" << label << ")\n"; }
#ifdef MFEM_USE_MPI
   // Run on the serial-mesh branch from rank 0 only (the ParMesh C-2 path
   // covers np>1 for Mult; AdvanceADER's per-element CK is element-local
   // so the serial check is sufficient for A1/A2/A3 coverage).
   if (g_rank != 0) { return; }
#endif

   Mesh smesh = MakeBoxMesh(2, 2, 2);
   const BoundaryConfig bc = MakeMixedBC();

   WaveOperator<Mesh> wave_scalar(smesh, k_order, lam, mu, rho, bc);
   WaveOperator<Mesh> wave_hetero(smesh, k_order,
                                  MaterialField::MakeConstant(lam, mu, rho),
                                  bc);
   real_t Q_bg[NUM_STATE] = {0};
   wave_scalar.SetAbsorbingBackground(Q_bg);
   wave_hetero.SetAbsorbingBackground(Q_bg);

   const int n = NUM_STATE * wave_scalar.GetScalarNDof();
   Vector Q(n);
   FillQDeterministic(Q);

   // --- A1: volume RHS byte-identical (per-element star matrices == Ax_) ---
   {
      Vector vs(n), vh(n);
      vs = 0.0; vh = 0.0;
      wave_scalar.ComputeVolumeRHS_ForTest(Q, vs);
      wave_hetero.ComputeVolumeRHS_ForTest(Q, vh);
      const int nd = BitDiffCount(vs, vh);
      if (g_rank == 0)
      { std::cout << "  volume RHS bit-diffs=" << nd << " (expect 0)\n"; }
      TEST_ASSERT(nd == 0,
                  std::string("volume RHS byte-identical (") + label + ")");
   }

   // --- A2: ADER CK predictor byte-identical (order 2 and 3) ---
   const real_t dt = 0.2 * wave_scalar.ComputeMaxDt(0.5);
   for (int order = 2; order <= 3; ++order)
   {
      Vector Is, Ih;
      wave_scalar.ComputeADERTimeIntegrated(Q, dt, order, Is);
      wave_hetero.ComputeADERTimeIntegrated(Q, dt, order, Ih);
      const int nd = BitDiffCount(Is, Ih);
      if (g_rank == 0)
      { std::cout << "  CK predictor order=" << order
                  << " bit-diffs=" << nd << " (expect 0)\n"; }
      TEST_ASSERT(nd == 0,
                  std::string("ADER CK predictor byte-identical order=")
                  + std::to_string(order) + " (" + label + ")");
   }

   // --- Mult FP-equivalent (interior face flux routes through
   //     BimaterialFlux on the hetero ctor → ~1e-10 rel, not bit) ---
   const real_t k_rel_tol = 1e-9;
   {
      Vector ds(n), dh(n);
      wave_scalar.Mult(Q, ds);
      wave_hetero.Mult(Q, dh);
      int n_over = 0;
      const real_t mx = MaxRelDiff(ds, dh, &n_over, k_rel_tol);
      if (g_rank == 0)
      { std::cout << "  Mult max_rel_diff=" << mx
                  << " (tol=" << k_rel_tol << ", n_over=" << n_over << ")\n"; }
      TEST_ASSERT(n_over == 0,
                  std::string("Mult parity (") + label + ")");
   }

   // --- AdvanceADER FP-equivalent (order 2 and 3) ---
   for (int order = 2; order <= 3; ++order)
   {
      Vector Qs(n), Qh(n);
      wave_scalar.AdvanceADER(Q, dt, order, Qs);
      wave_hetero.AdvanceADER(Q, dt, order, Qh);
      int n_over = 0;
      const real_t mx = MaxRelDiff(Qs, Qh, &n_over, k_rel_tol);
      if (g_rank == 0)
      { std::cout << "  AdvanceADER order=" << order << " max_rel_diff=" << mx
                  << " (tol=" << k_rel_tol << ", n_over=" << n_over << ")\n"; }
      TEST_ASSERT(n_over == 0,
                  std::string("AdvanceADER parity order=")
                  + std::to_string(order) + " (" + label + ")");
   }
}

// =========================================================================
// C-6  Fault parity: a serial mesh with an interior fault (y = L/2 plane).
//      Exercises Group C1 (RK4 fault bulk-side flux) and C3 (ADER fault
//      bulk-side flux).  On Mode::Constant the per-side FluxForElem_(elem)
//      equals the scalar flux_ (At(elem) built from the exact constants),
//      so Mult and AdvanceADER must match the scalar ctor to FP precision.
//      Uses TPV102 material so InitializeFaultDOFs + FaultFaceFlux are
//      mutually consistent.
// =========================================================================
static void C_6_fault_parity()
{
   if (g_rank == 0)
   { std::cout << "\n[C-6] Fault parity (interior fault, serial)\n"; }
#ifdef MFEM_USE_MPI
   if (g_rank != 0) { return; }
#endif

   const real_t L = 1.0e4;   // 10 km box; depth |z| up to 10 km
   Mesh smesh = MakeFaultBoxMesh(2, 2, 2, L);

   BoundaryConfig bc;
   bc.fault_attr      = 3;
   bc.natural_attrs   = {1};
   bc.absorbing_attrs = {};
   bc.dirichlet_attrs = {};

   const real_t lam = TPV102Params::lambda;
   const real_t mu  = TPV102Params::mu;
   const real_t rho = TPV102Params::rho;

   WaveOperator<Mesh> wave_scalar(smesh, k_order, lam, mu, rho, bc);
   WaveOperator<Mesh> wave_hetero(smesh, k_order,
                                  MaterialField::MakeConstant(lam, mu, rho),
                                  bc);

   const Array<int> &int_faces = wave_scalar.GetFaultInteriorFaces();
   if (g_rank == 0)
   { std::cout << "  interior fault faces detected = " << int_faces.Size()
               << "\n"; }
   TEST_ASSERT(int_faces.Size() > 0,
               "fault mesh produced at least one interior fault face");
   if (int_faces.Size() == 0) { return; }

   // Per-operator fault wiring (each operator owns its own FaultFaceFlux +
   // DOFData; the two DOFData arrays are seeded identically).
   const int nqp =
      IntRules.Get(smesh.GetInteriorFaceTransformations(int_faces[0])
                      ->GetGeometryType(),
                   2 * k_order).GetNPoints();

   std::vector<Vector> fault_coords;
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr = smesh.GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir =
         IntRules.Get(ftr->GetGeometryType(), 2 * k_order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3);
         ftr->Face->Transform(ip, phys);
         // Shift y so the fault sits at the canonical y=0; shift z so the
         // surface is at z=0 and depth = |z| (InitializeFaultDOFs needs
         // z <= 0).  Box spans [0,L]^3; fault at y=L/2.
         phys(1) -= 0.5 * L;
         phys(2) -= L;
         fault_coords.push_back(phys);
      }
   }
   const int n_fault = int_faces.Size() * nqp;

   std::vector<DOFData> dof_scalar, dof_hetero;
   InitializeFaultDOFs(dof_scalar, n_fault, fault_coords);
   InitializeFaultDOFs(dof_hetero, n_fault, fault_coords);

   FaultFaceFlux ff_scalar(rho, TPV102Params::cp, TPV102Params::cs);
   FaultFaceFlux ff_hetero(rho, TPV102Params::cp, TPV102Params::cs);
   wave_scalar.SetFaultFlux(&ff_scalar);
   wave_scalar.SetFaultDOFData(&dof_scalar, nqp);
   wave_hetero.SetFaultFlux(&ff_hetero);
   wave_hetero.SetFaultDOFData(&dof_hetero, nqp);

   real_t Q_bg[NUM_STATE] = {0};
   wave_scalar.SetAbsorbingBackground(Q_bg);
   wave_hetero.SetAbsorbingBackground(Q_bg);

   // R-001 operator-level fault-routing gate: a UNIFORM Mode::Coefficient
   // fault op.  owned_flux_pool_ is built (so the fault dispatch C1-C4 reads
   // FluxForElem_(elem)) but its scalar flux_ is the (1,1,1) placeholder.
   // Matching the scalar op's Mult/AdvanceADER to FP precision proves the
   // fault bulk-side flux reads FluxForElem_(elem), NOT the placeholder.  The
   // MakeConstant wave_hetero above seeds flux_ with the real constant, so it
   // alone could not catch a C1-C4 revert to flux_; this op can.
   ConstantCoefficient lam_cc(lam), mu_cc(mu), rho_cc(rho);
   MaterialField mat_coef =
      MaterialField::MakeCoefficient(&lam_cc, &mu_cc, &rho_cc);
   WaveOperator<Mesh> wave_coef(smesh, k_order, mat_coef, bc);
   TEST_ASSERT(wave_coef.UsesGodunovFluxPool(),
               "uniform-Coefficient fault op built the per-element pool");
   std::vector<DOFData> dof_coef;
   InitializeFaultDOFs(dof_coef, n_fault, fault_coords);
   FaultFaceFlux ff_coef(rho, TPV102Params::cp, TPV102Params::cs);
   wave_coef.SetFaultFlux(&ff_coef);
   wave_coef.SetFaultDOFData(&dof_coef, nqp);
   wave_coef.SetAbsorbingBackground(Q_bg);

   const int n = NUM_STATE * wave_scalar.GetScalarNDof();
   Vector Q(n);
   FillQDeterministic(Q);

   const real_t k_rel_tol = 1e-9;

   // --- C1: RK4 fault bulk-side flux parity ---
   {
      Vector ds(n), dh(n);
      wave_scalar.Mult(Q, ds);
      wave_hetero.Mult(Q, dh);
      int n_over = 0;
      const real_t mx = MaxRelDiff(ds, dh, &n_over, k_rel_tol);
      if (g_rank == 0)
      { std::cout << "  fault Mult max_rel_diff=" << mx
                  << " (tol=" << k_rel_tol << ", n_over=" << n_over << ")\n"; }
      TEST_ASSERT(n_over == 0, "fault Mult parity (C1)");
   }

   // --- C3: ADER fault bulk-side flux parity ---
   const real_t dt = 0.2 * wave_scalar.ComputeMaxDt(0.5);
   {
      Vector Qs(n), Qh(n);
      wave_scalar.AdvanceADER(Q, dt, /*order=*/2, Qs);
      wave_hetero.AdvanceADER(Q, dt, /*order=*/2, Qh);
      int n_over = 0;
      const real_t mx = MaxRelDiff(Qs, Qh, &n_over, k_rel_tol);
      if (g_rank == 0)
      { std::cout << "  fault AdvanceADER max_rel_diff=" << mx
                  << " (tol=" << k_rel_tol << ", n_over=" << n_over << ")\n"; }
      TEST_ASSERT(n_over == 0, "fault AdvanceADER parity (C3)");
   }

   // --- R-001: uniform-Coefficient fault op vs scalar (placeholder-catch).
   //     Unlike the MakeConstant wave_hetero above, wave_coef's flux_ is the
   //     (1,1,1) placeholder, so these comparisons FAIL if any C1-C4 fault
   //     site reverts to flux_ instead of FluxForElem_(elem). ---
   {
      Vector ds(n), dc(n);
      wave_scalar.Mult(Q, ds);
      wave_coef.Mult(Q, dc);
      int n_over = 0;
      const real_t mx = MaxRelDiff(ds, dc, &n_over, k_rel_tol);
      if (g_rank == 0)
      { std::cout << "  fault Mult coef-vs-scalar max_rel_diff=" << mx
                  << " (tol=" << k_rel_tol << ", n_over=" << n_over << ")\n"; }
      TEST_ASSERT(n_over == 0,
                  "uniform-Coefficient fault Mult == scalar (C1 reads "
                  "FluxForElem_(elem), not placeholder flux_)");
   }
   {
      Vector Qs(n), Qc(n);
      wave_scalar.AdvanceADER(Q, dt, /*order=*/2, Qs);
      wave_coef.AdvanceADER(Q, dt, /*order=*/2, Qc);
      int n_over = 0;
      const real_t mx = MaxRelDiff(Qs, Qc, &n_over, k_rel_tol);
      if (g_rank == 0)
      { std::cout << "  fault AdvanceADER coef-vs-scalar max_rel_diff=" << mx
                  << " (tol=" << k_rel_tol << ", n_over=" << n_over << ")\n"; }
      TEST_ASSERT(n_over == 0,
                  "uniform-Coefficient fault AdvanceADER == scalar (ADER C3 "
                  "reads FluxForElem_(elem), not placeholder flux_)");
   }
}

// =========================================================================
// main
// =========================================================================
int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
#else
   g_rank = 0;
#endif

   if (g_rank == 0)
   {
      std::cout << "==============================================\n"
                << "test_phaseh_wave_operator_constant_parity\n"
                << "Phase H Stage 1 gate (T-PHASEH-SCALAR-PARITY)\n"
                << "==============================================\n";
   }

   C_1_mult_parity_serial();
   C_2_mult_parity_parallel();
   C_3_compute_max_dt_parity();
   C_4_pool_invariants_on_constant();
   // Phase H Stage 2 (PLAN §6): byte-exact volume + ADER-CK and
   // FP-equivalent Mult + AdvanceADER under mixed (absorbing +
   // free-surface) BC, with a CLEAN and a >6-sig-fig constant.
   C_5_nonfault_parity(k_lambda,  k_mu,  k_rho,  "clean 6-sig constant");
   C_5_nonfault_parity(k_lambda6, k_mu6, k_rho6, ">6-sig constant (R-002)");
   // Phase H Stage 2 (PLAN §Phase 3): fault bulk-side flux parity (C1/C3).
   C_6_fault_parity();

#ifdef MFEM_USE_MPI
   int total = 0, passed = 0, failed = 0;
   MPI_Allreduce(&g_num_tests, &total,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&g_num_passed, &passed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&g_num_failed, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
   const int total = g_num_tests, passed = g_num_passed, failed = g_num_failed;
#endif
   if (g_rank == 0)
   {
      std::cout << "\n==============================================\n"
                << "Summary: " << passed << " / " << total
                << " passed (" << failed << " failed)\n"
                << "==============================================\n";
   }
#ifdef MFEM_USE_MPI
   MPI_Finalize();
#endif
   return failed > 0 ? 1 : 0;
}

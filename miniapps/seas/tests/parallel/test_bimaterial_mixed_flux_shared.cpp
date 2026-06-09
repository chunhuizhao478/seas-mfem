// Phase 3 (PLAN_mixed_flux_hetero_riemann.md), P3-1 fix: MPI (np=2) coverage of
// the bi-material SHARED central-face path that the serial dispatch test
// (tests/unit/test_bimaterial_mixed_flux_dispatch.cpp) structurally cannot reach.
//
// The serial test exercises the fully-local 2-sided interior central arm of
// BuildPerFaceCentralFluxMatrices_ + InteriorFaceFlux_.  THIS test covers the
// cross-rank SHARED arm: the `sf` loop in BuildPerFaceCentralFluxMatrices_
// (dynamic/bimaterial_wave_operator.inl) + the central branch of
// SharedInteriorFaceFlux_ — the production case when a fault-adjacent interior
// face lands on a partition seam.
//
//   Test 3.4a (Dispatch_SharedCentralFace_Build):
//     On an np=2 ParMesh where a fault-adjacent non-fault interior face is split
//     across the partition seam (a SHARED central face), SetMixedFluxMode(Adjacent)
//     with a het Mode::Coefficient material builds the side-0 shared central
//     matrices and SharedInteriorFaceFlux_ dispatches the central F*
//     (== ApplyPerFaceFlux(per_face_central_flux_[mf], ...) bit-for-bit).  Since
//     cross-rank Phase 2 the matrices are built from the TRUE peer material, so
//     check (c) asserts c[0] != c[1] for a lateral mu(x) (see below).
//     (SetSeamContinuous(true) is still called but is now a NO-OP — deprecated.)
//
//   Test 3.4c (Dispatch_DepthProfile_NoFalseAbort):
//     A depth-only (mu(z)) Mode::Coefficient material on a 2D (x,z) grid whose
//     partition seam carries a fault-adjacent z-normal central SHARED face — and
//     whose material genuinely varies in depth across that seam — BUILDS with no
//     abort.  Cross-rank Phase 2 REMOVED the seam-continuity affirmation guard
//     entirely (the neighbour material is the true peer), so this builds for ANY
//     material, depth-varying or lateral.  (The fault itself is x-normal so the
//     FaultBasis strike/dip frame, up=(0,0,1), is well defined; the seam face is
//     z-normal.)
//
//   Test 3.4b (Dispatch_SharedCentralFace_RequiresSeamContinuous) — OBSOLETE.
//     The former guard (MFEM_VERIFY(mode==Constant || seam_continuous_) in the
//     shared arm of BuildPerFaceCentralFluxMatrices_) was REMOVED in cross-rank
//     Phase 2: the central build now reads the TRUE peer material, so no
//     seam-continuity affirmation is needed and there is no abort to test.
//     `seam_continuous` is deprecated (parsed for back-compat, never read).
//
// Determinism (ParMETIS is nondeterministic): we build a small serial mesh and
// pass an EXPLICIT partition array to ParMesh so a chosen fault-adjacent interior
// face is guaranteed to straddle the rank seam.
//
// Run: mpirun -np 2 ./seas_test_bimaterial_mixed_flux_shared

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

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

namespace
{
int g_num_tests = 0, g_num_passed = 0, g_num_failed = 0;
int g_rank = 0;

#define TEST_ASSERT(cond, msg) do {                                          \
   ++g_num_tests;                                                            \
   if (cond) { ++g_num_passed;                                               \
      if (g_rank == 0) { std::cout << "  PASSED: " << msg << "\n"; } }       \
   else { ++g_num_failed;                                                    \
      /* print failures on EVERY rank — an MPI dispatch defect may surface   \
         only on the non-zero rank that owns the seam face. */               \
      std::cerr << "  FAILED [rank " << g_rank << ", line " << __LINE__       \
                << "]: " << msg << "\n"; }                                   \
} while (0)

// Expose the protected SharedInteriorFaceFlux_ for direct per-face testing
// (mirrors the serial dispatch test's TestableBimat shim).
template <typename MeshType>
class TestableBimat : public BimaterialWaveOperator<MeshType>
{
public:
   using BimaterialWaveOperator<MeshType>::BimaterialWaveOperator;     // inherit ctor
   using BimaterialWaveOperator<MeshType>::InteriorFaceFlux_;          // protected -> public
   using BimaterialWaveOperator<MeshType>::SharedInteriorFaceFlux_;    // protected -> public
};

// Row of `n` unit hexes along x; the interior face at x == fault_x is tagged
// FAULT (attr 3) (x-normal), all external faces FREE (attr 1).  Used by 3.4a.
Mesh BuildHexRowFaultMesh(int n, int fault_x)
{
   Mesh mesh(3, (n + 1) * 4, n, 0);
   auto vid = [](int ix, int iy, int iz) { return ix * 4 + iy * 2 + iz; };
   for (int ix = 0; ix <= n; ++ix)
      for (int iy = 0; iy < 2; ++iy)
         for (int iz = 0; iz < 2; ++iz)
         {
            real_t v[3] = { static_cast<real_t>(ix), static_cast<real_t>(iy),
                            static_cast<real_t>(iz) };
            mesh.AddVertex(v);
         }
   for (int e = 0; e < n; ++e)
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
      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      real_t cx = 0.0;
      for (int i = 0; i < fv.Size(); ++i) { cx += mesh.GetVertex(fv[i])[0]; }
      cx /= fv.Size();
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         if (std::abs(cx - fault_x) < 1e-9)
         { mesh.AddBdrQuad(fv[0], fv[1], fv[2], fv[3], 3); }   // FAULT (x-normal)
         continue;
      }
      mesh.AddBdrQuad(fv[0], fv[1], fv[2], fv[3], 1);           // FREE (external)
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// nx-by-nz grid of unit hexes (ny=1).  The interior face at x == fault_x is the
// FAULT (attr 3) — x-normal, so the FaultBasis strike/dip frame (up=(0,0,1)) is
// well defined.  Interior z-normal faces are left untagged; partitioning along z
// turns a fault-adjacent z-normal face into a SHARED central face.  Used by 3.4c.
Mesh BuildHexGridXZFaultMesh(int nx, int nz, int fault_x)
{
   const int nvx = nx + 1, nvz = nz + 1;
   Mesh mesh(3, nvx * 2 * nvz, nx * nz, 0);
   auto vid = [&](int ix, int iy, int iz)
   { return (ix * 2 + iy) * nvz + iz; };
   for (int ix = 0; ix < nvx; ++ix)
      for (int iy = 0; iy < 2; ++iy)
         for (int iz = 0; iz < nvz; ++iz)
         {
            real_t v[3] = { static_cast<real_t>(ix), static_cast<real_t>(iy),
                            static_cast<real_t>(iz) };
            mesh.AddVertex(v);
         }
   for (int ix = 0; ix < nx; ++ix)
      for (int iz = 0; iz < nz; ++iz)
      {
         const int v[8] = {
            vid(ix,   0, iz),   vid(ix+1, 0, iz),   vid(ix+1, 1, iz),   vid(ix,   1, iz),
            vid(ix,   0, iz+1), vid(ix+1, 0, iz+1), vid(ix+1, 1, iz+1), vid(ix,   1, iz+1)
         };
         mesh.AddHex(v, 1);
      }
   mesh.FinalizeTopology();
   for (int f = 0; f < mesh.GetNumFaces(); ++f)
   {
      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      real_t cx = 0.0, cz = 0.0;
      for (int i = 0; i < fv.Size(); ++i)
      { cx += mesh.GetVertex(fv[i])[0]; cz += mesh.GetVertex(fv[i])[2]; }
      cx /= fv.Size(); cz /= fv.Size();
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         // x-normal interior face at x==fault_x => FAULT (constant cx across it).
         if (std::abs(cx - fault_x) < 1e-9)
         { mesh.AddBdrQuad(fv[0], fv[1], fv[2], fv[3], 3); }
         continue;   // other interior (incl. z-normal seam faces): untagged
      }
      mesh.AddBdrQuad(fv[0], fv[1], fv[2], fv[3], 1);   // FREE (external)
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

real_t WorstRel(const real_t *a, const real_t *b)
{
   real_t scale = 0.0;
   for (int i = 0; i < NUM_STATE; ++i)
   {
      scale = std::max(scale, std::max(std::abs(a[i]), std::abs(b[i])));
   }
   if (scale == 0.0) { return 0.0; }
   real_t m = 0.0;
   for (int i = 0; i < NUM_STATE; ++i)
   {
      m = std::max(m, std::abs(a[i] - b[i]) / scale);
   }
   return m;
}

struct StatePair { real_t self[NUM_STATE]; real_t nbr[NUM_STATE]; };

std::vector<StatePair> MakeStatePairs()
{
   std::vector<StatePair> ps(4);
   for (auto &p : ps)
   {
      for (int c = 0; c < NUM_STATE; ++c) { p.self[c] = 0.0; p.nbr[c] = 0.0; }
   }
   for (int c = 0; c < SXY; ++c) { ps[0].self[c] = 1.0e7; ps[0].nbr[c] = 0.5e7; }
   ps[1].self[SXY] = 2.0e6; ps[1].nbr[SXY] = -1.0e6;
   ps[1].self[SXZ] = 1.0e6; ps[1].nbr[SXZ] =  3.0e6;
   for (int c = VX; c < NUM_STATE; ++c) { ps[2].self[c] = 1.0; ps[2].nbr[c] = -0.5; }
   ps[3].self[SXX] = 7.0e6; ps[3].self[VX] = 0.8; ps[3].self[SYZ] = -2.0e6;
   ps[3].nbr[SXX]  = -3.0e6; ps[3].nbr[VY]  = 0.4; ps[3].nbr[SXZ]  =  1.5e6;
   return ps;
}

}  // namespace

// =========================================================================
// Drive the shared central-face path on an explicitly-partitioned ParMesh.
// `seam_axis` selects the partition direction (0=x, 2=z); `seam_pos` is the
// coordinate of the interior face placed ON the rank seam.  Building the
// operator with Adjacent exercises the shared central arm (SetSeamContinuous is a
// deprecated no-op since cross-rank Phase 2); when `value_check` is set the
// dispatched side-0 F* is checked against the central store, and the stored halves
// c[0] != c[1] reflect the TRUE peer material.
//
// Returns the number of SHARED central faces found across ranks.
// =========================================================================
#ifdef MFEM_USE_MPI
static int RunSharedCentralCase(const char *label, Mesh &&smesh_in,
                                int seam_axis, double seam_pos,
                                MaterialField mat, bool value_check)
{
   if (g_rank == 0) { std::cout << "\n-- " << label << " --\n"; }

   Mesh smesh(std::move(smesh_in));

   // EXPLICIT partition: elements with seam-axis centroid < seam_pos -> rank 0,
   // others -> rank 1.  The interior face at `seam_pos` then straddles the seam.
   const int ne = smesh.GetNE();
   std::vector<int> part(ne, 0);
   for (int e = 0; e < ne; ++e)
   {
      Array<int> ev;
      smesh.GetElementVertices(e, ev);
      real_t ca = 0.0;
      for (int v = 0; v < ev.Size(); ++v)
      { ca += smesh.GetVertex(ev[v])[seam_axis]; }
      ca /= ev.Size();
      part[e] = (ca < seam_pos) ? 0 : 1;
   }
   ParMesh pmesh(MPI_COMM_WORLD, smesh, part.data());
   smesh.Clear();

   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {};

   const int order = 1;
   TestableBimat<ParMesh> op(pmesh, order, mat, bc);
   op.SetSeamContinuous(true);                       // NO-OP since cross-rank Phase 2
                                                     // (deprecated; the central
                                                     // build reads the TRUE peer)
   op.SetMixedFluxMode(MixedFluxMode::Adjacent);     // builds central matrices

   // Lifecycle invariant (IMPL-8): map size == set size.
   const auto &cset = op.GetCentralFluxFaceSet();
   const auto &cmap = op.GetPerFaceCentralFlux();
   TEST_ASSERT(cmap.size() == cset.size(),
               std::string(label) + ": per_face_central_flux_.size() == "
               "central_flux_face_set_.size() (IMPL-8 lifecycle)");

   // Locate a SHARED central face: a shared-face index sf whose mesh face is in
   // the central set and is not a fault face.
   int sf_central = -1, mf_central = -1, local_elem = -1;
   const int n_shared = pmesh.GetNSharedFaces();
   const Array<int> &fault_shared = op.GetFaultSharedFaces();
   auto is_fault_sf = [&](int sf) {
      for (int i = 0; i < fault_shared.Size(); ++i)
      { if (fault_shared[i] == sf) { return true; } }
      return false;
   };
   for (int sf = 0; sf < n_shared; ++sf)
   {
      const int mf = pmesh.GetSharedFace(sf);
      if (cset.count(mf) == 0) { continue; }
      if (is_fault_sf(sf)) { continue; }
      auto *ftr = pmesh.GetSharedFaceTransformations(sf);
      if (!ftr) { continue; }
      sf_central = sf; mf_central = mf; local_elem = ftr->Elem1No;
      break;
   }

   int local_found = (sf_central >= 0) ? 1 : 0;
   int global_found = 0;
   MPI_Allreduce(&local_found, &global_found, 1, MPI_INT, MPI_SUM,
                 MPI_COMM_WORLD);
   TEST_ASSERT(global_found > 0,
               std::string(label) + ": at least one rank found a SHARED "
               "central (fault-adjacent, non-fault) face on the seam");

   if (value_check && sf_central >= 0)
   {
      const std::vector<StatePair> pairs = MakeStatePairs();
      const real_t dummy_nor[3] = {1.0, 0.0, 0.0};   // ignored by the dispatch

      // (a) The dispatched side-0 F* must route to the central store: equal to
      // ApplyPerFaceFlux(per_face_central_flux_[mf], ...) bit for bit.
      const auto &c = cmap.at(mf_central);
      real_t worst_store = 0.0;
      for (const auto &p : pairs)
      {
         real_t Fdisp[NUM_STATE], Fstore[NUM_STATE];
         op.SharedInteriorFaceFlux_(mf_central, p.self, p.nbr, dummy_nor, Fdisp);
         BimaterialFlux::ApplyPerFaceFlux(c[0], c[1], p.self, p.nbr, Fstore);
         worst_store = std::max(worst_store, WorstRel(Fdisp, Fstore));
      }
      TEST_ASSERT(worst_store == 0.0,
                  std::string(label) + ": shared central dispatch == "
                  "ApplyPerFaceFlux(per_face_central_flux_[mf], ...) bit-for-bit");

      // (b) The shared face must carry a neighbour-material entry — now the TRUE
      // peer (Cross-rank Phase 2: ExchangeBiMaterialNeighbours_ reads it via
      // MaterialAtNbr_, keyed by sf), not the former R-004 local-side stub.
      const auto &nbr_map = op.GetSharedFaceNeighbourMaterial();
      auto nbr_it = nbr_map.find(sf_central);
      TEST_ASSERT(nbr_it != nbr_map.end(),
                  std::string(label) + ": shared face has a neighbour-material "
                  "entry");

      // (c) Cross-rank Phase 2 content check.  This case's material is a LATERAL
      // mu(x) (== lam(x)) that varies ACROSS the seam, so the TRUE peer-centroid
      // material differs from the local-centroid material.  The former local-side
      // stub forced neighbour==local (=> c[0]==c[1]); with the real exchange the
      // two stored half-Jacobians c[0]=½A_local and c[1]=½A_peer must now DIFFER,
      // proving the central build consumed the genuine cross-rank peer material.
      if (nbr_it != nbr_map.end())
      {
         const auto &nbrmat = nbr_it->second;                 // true peer (l,m,r)
         const auto &locmat = op.GetPerElementMaterial()[local_elem];
         const real_t dmat = std::max({ std::abs(nbrmat[0] - locmat[0]),
                                        std::abs(nbrmat[1] - locmat[1]),
                                        std::abs(nbrmat[2] - locmat[2]) });
         TEST_ASSERT(dmat > 0.0,
                     std::string(label) + ": peer material != local (real "
                     "cross-rank exchange, not the local-side stub)");

         const DenseMatrix &c0 = c[0];   // ½A_local
         const DenseMatrix &c1 = c[1];   // ½A_peer (true peer material)
         real_t scale = 0.0, worst = 0.0;
         for (int i = 0; i < NUM_STATE; ++i)
            for (int j = 0; j < NUM_STATE; ++j)
            {
               scale = std::max(scale, std::max(std::abs(c0(i, j)),
                                                std::abs(c1(i, j))));
            }
         for (int i = 0; i < NUM_STATE; ++i)
            for (int j = 0; j < NUM_STATE; ++j)
            {
               worst = std::max(worst, std::abs(c0(i, j) - c1(i, j)));
            }
         const real_t rel = (scale > 0.0) ? worst / scale : worst;
         TEST_ASSERT(rel > 1.0e-9,
                     std::string(label) + ": stored shared central halves "
                     "c[0] != c[1] (built from the TRUE peer material via the "
                     "cross-rank exchange; rel " + std::to_string(rel) + ")");
      }
   }

   return global_found;
}
#endif

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
   MPI_Init(&argc, &argv);
   int nprocs = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &g_rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   if (nprocs != 2)
   {
      if (g_rank == 0)
      {
         std::cerr << "[P3-1] requires np=2, got " << nprocs
                   << " — TEST FAILED.  Run via `mpirun -np 2 "
                      "./seas_test_bimaterial_mixed_flux_shared`.\n";
      }
      MPI_Finalize();
      return 1;
   }

   if (g_rank == 0)
   {
      std::cout << "\n=== Phase 3 P3-1: BimaterialWaveOperator SHARED "
                   "central-face dispatch (np=2) ===\n";
   }

   // -----------------------------------------------------------------------
   // Test 3.4a: het mu(x) Mode::Coefficient (lateral).  5-hex row along x; fault
   // at x=3 (between elems 2,3); seam at x=4 (between elems 3,4) => the
   // fault-adjacent non-fault face (3,4) is a SHARED central face.  Cross-rank
   // Phase 2: the central halves are built from the TRUE peer (c[0] != c[1]).
   // -----------------------------------------------------------------------
   {
      auto lam_fn = [](const Vector &x) -> real_t
      { return 30.0e9 * (1.0 + 0.2 * x(0)); };
      auto mu_fn = [](const Vector &x) -> real_t
      { return 30.0e9 * (1.0 + 0.2 * x(0)); };
      FunctionCoefficient lam_c(lam_fn), mu_c(mu_fn);
      ConstantCoefficient rho_c(2670.0);
      RunSharedCentralCase("Test 3.4a (lateral mu(x), true-peer central)",
                           BuildHexRowFaultMesh(/*n=*/5, /*fault_x=*/3),
                           /*seam_axis=*/0, /*seam_pos=*/4.0,
                           MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c),
                           /*value_check=*/true);
   }

   // -----------------------------------------------------------------------
   // Test 3.4c: depth-only mu(z) Mode::Coefficient on a 2-by-5 (x,z) grid.  Fault
   // is the x-normal face at x=1 (well-defined frame); the seam splits along z at
   // z=3, so a fault-adjacent z-normal interior face is a SHARED central face whose
   // neighbouring elements have DISTINCT depth materials.  Must BUILD with no abort
   // (cross-rank Phase 2 removed the seam-continuity guard entirely).
   // -----------------------------------------------------------------------
   {
      auto lam_fn = [](const Vector &x) -> real_t
      { return 30.0e9 * (1.0 + 0.2 * x(2)); };
      auto mu_fn = [](const Vector &x) -> real_t
      { return 30.0e9 * (1.0 + 0.2 * x(2)); };
      FunctionCoefficient lam_c(lam_fn), mu_c(mu_fn);
      ConstantCoefficient rho_c(2670.0);
      RunSharedCentralCase("Test 3.4c (depth-only mu(z), no abort)",
                           BuildHexGridXZFaultMesh(/*nx=*/2, /*nz=*/5,
                                                   /*fault_x=*/1),
                           /*seam_axis=*/2, /*seam_pos=*/3.0,
                           MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c),
                           /*value_check=*/false);
   }

   // Test 3.4b is OBSOLETE: cross-rank Phase 2 removed the seam-continuity abort
   // entirely (the central build reads the TRUE peer material) — there is no
   // abort to test.  See the file header.

   int total = 0, passed = 0, failed = 0;
   MPI_Allreduce(&g_num_tests, &total,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&g_num_passed, &passed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&g_num_failed, &failed, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

   if (g_rank == 0)
   {
      std::cout << "\n==============================================\n"
                << "Summary (across ranks): " << passed << " / " << total
                << " passed (" << failed << " failed)\n"
                << "==============================================\n";
   }
   MPI_Finalize();
   return failed > 0 ? 1 : 0;
#else
   (void)argc; (void)argv;
   std::cout << "SKIPPED: MFEM_USE_MPI not defined\n";
   return 0;
#endif
}

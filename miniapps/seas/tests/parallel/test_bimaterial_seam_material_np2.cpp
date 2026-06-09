// Cross-rank material exchange — Phase 1 gate (PLAN_cross_rank_material_exchange
// _2026-06-06.md).  np=2 coverage of the ONE uniform material accessor
// (BimaterialWaveOperator::MaterialAtLocal_ / MaterialAtNbr_) reading the LOCAL
// element and the PEER (face-neighbour, across a partition seam) element in all
// three material representations (Constant, Coefficient, GridFunction).
//
//   accessor_local_peer_all_modes (the headline P1 gate):
//     On an np=2 ParMesh with an EXPLICIT partition that places a non-fault
//     interior face ON the rank seam, for EACH mode assert
//       MaterialAtLocal_(local_elem, ip) == the LOCAL element's material, and
//       MaterialAtNbr_(ftr, ip_peer)    == the PEER element's material,
//     per component to 1e-12, computed INDEPENDENTLY from the known mesh
//     geometry (not by re-reading the same accessor).  The Coefficient case is
//     checked at a smooth profile (exact centroid value) AND at a halfspace JUMP
//     (multiple interior IPs ⇒ each side reads its PURE material, no mix — the
//     R-301-class no-projection / no-pollution check).  The GridFunction case
//     uses a per-element step (constant within each element, jump across the
//     seam) so the local GetValue and the peer shape×FaceNbrData interpolation
//     must each return the owning element's value.
//
//   accessor_coeff_no_material_mpi:
//     Constant + Coefficient perform NO material ParGridFunction exchange
//     (GetMaterialGfExchangeCount() == 0); GridFunction performs exactly 3 (one
//     per source GF).  Confirms the analytic modes are globally evaluable.
//
// Determinism (ParMETIS is nondeterministic): a small serial mesh + an EXPLICIT
// partition array guarantees the chosen interior face straddles the rank seam.
//
// Run: mpirun -np 2 ./seas_test_bimaterial_seam_material_np2

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
#include <functional>
#include <iostream>
#include <memory>
#include <string>
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
      std::cerr << "  FAILED [rank " << g_rank << ", line " << __LINE__       \
                << "]: " << msg << "\n"; }                                   \
} while (0)

// Expose the protected accessors for direct per-face testing.
template <typename MeshType>
class TestableBimat : public BimaterialWaveOperator<MeshType>
{
public:
   using BimaterialWaveOperator<MeshType>::BimaterialWaveOperator;   // inherit ctor
   using BimaterialWaveOperator<MeshType>::MaterialAtLocal_;         // protected -> public
   using BimaterialWaveOperator<MeshType>::MaterialAtNbr_;           // protected -> public
};

// (P2-T4) A Coefficient that reads T.GetIntPoint() instead of the explicit `ip` —
// the class P2-005 hardens against.  Returns a positive, monotone-in-x value so
// the operator's rho>0 / lambda+2mu>0 asserts pass.  If a caller forgets to set
// the transform's int point before Eval, this reads the WRONG physical point.
struct IPProbeCoefficient : public mfem::Coefficient
{
   real_t base;
   explicit IPProbeCoefficient(real_t b) : base(b) {}
   real_t Eval(mfem::ElementTransformation &T,
               const mfem::IntegrationPoint & /*ip*/) override
   {
      mfem::Vector x(3);
      T.Transform(T.GetIntPoint(), x);   // deliberately uses GetIntPoint(), not ip
      return base * (1.0 + 0.1 * x(0));
   }
};

// Row of `n` unit hexes along x; the interior face at x == fault_x is tagged
// FAULT (attr 3, x-normal), all external faces FREE (attr 1).
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

real_t RelErr(real_t a, real_t b)
{
   const real_t s = std::max(std::abs(a), std::abs(b));
   return (s == 0.0) ? std::abs(a - b) : std::abs(a - b) / s;
}

// Worst entrywise relative difference between two equally-sized DenseMatrices.
real_t MatRel(const DenseMatrix &a, const DenseMatrix &b)
{
   real_t scale = 0.0, worst = 0.0;
   for (int i = 0; i < a.Height(); ++i)
      for (int j = 0; j < a.Width(); ++j)
      {
         scale = std::max(scale, std::max(std::abs(a(i, j)), std::abs(b(i, j))));
         worst = std::max(worst, std::abs(a(i, j) - b(i, j)));
      }
   return (scale > 0.0) ? worst / scale : worst;
}

}  // namespace

#ifdef MFEM_USE_MPI
// Locate the SHARED interior face whose physical centroid x ≈ seam_x.  Returns
// the shared-face index sf (>=0) on this rank, or -1 if not present.
static int FindSeamSharedFace(ParMesh &pmesh, double seam_x)
{
   const int n_shared = pmesh.GetNSharedFaces();
   for (int sf = 0; sf < n_shared; ++sf)
   {
      FaceElementTransformations *ftr = pmesh.GetSharedFaceTransformations(sf);
      if (!ftr) { continue; }
      Vector fc(3);
      ftr->Face->Transform(
         Geometries.GetCenter(ftr->Face->GetGeometryType()), fc);
      if (std::abs(fc(0) - seam_x) < 1e-9) { return sf; }
   }
   return -1;
}

// One mode's accessor check.  `expect(cx, lam, mu, rho)` returns the material at
// a physical centroid x INDEPENDENTLY of the accessor (from the known mesh).
// `piecewise_const` ⇒ the material is constant within each element (smooth modes
// false): then additionally probe off-centroid IPs and require the SAME pure
// per-side value (no mix / no pollution).
static void RunModeCase(
   const char *label, ParMesh &pmesh, int sf_seam, double seam_x,
   MaterialField mat, bool piecewise_const, std::size_t expect_gf_exch,
   const std::function<void(real_t, real_t &, real_t &, real_t &)> &expect)
{
   if (g_rank == 0) { std::cout << "\n-- " << label << " --\n"; }

   const int order = 1;
   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {};

   TestableBimat<ParMesh> op(pmesh, order, mat, bc);

   // accessor_coeff_no_material_mpi: analytic modes exchange no material GFs.
   TEST_ASSERT(op.GetMaterialGfExchangeCount() == expect_gf_exch,
               std::string(label) + ": material GF exchange count == "
               + std::to_string(expect_gf_exch));

   // The local rank may legitimately lack the seam face only if the partition
   // put both seam elements on one rank — here it never does, so require it.
   int local_ok = (sf_seam >= 0) ? 1 : 0;
   if (sf_seam < 0)
   {
      // Still participate in the reduction below.
      int g = 0; MPI_Allreduce(&local_ok, &g, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      TEST_ASSERT(g >= 1, std::string(label) + ": some rank owns the seam face");
      return;
   }

   // Re-fetch ftr AFTER construction (the ctor calls GetSharedFaceTransformations
   // internally, overwriting the mesh's shared scratch transform).  Do NOT call
   // any other GetSharedFaceTransformations until done with ftr.
   FaceElementTransformations *ftr = pmesh.GetSharedFaceTransformations(sf_seam);
   TEST_ASSERT(ftr != nullptr && ftr->Elem2 != nullptr,
               std::string(label) + ": seam face has a valid peer transform");
   if (!ftr || !ftr->Elem2)
   {
      int g = 0; MPI_Allreduce(&local_ok, &g, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      return;
   }

   const int local_elem = ftr->Elem1No;
   const IntegrationPoint &cref =
      Geometries.GetCenter(ftr->Elem1->GetGeometryType());

   // Independent local element transform: the overload that FILLS a caller-owned
   // object does NOT touch the mesh shared scratch that ftr->Elem1 aliases, so
   // ftr stays intact across every accessor call below (P1-001).
   IsoparametricTransformation Tloc;
   pmesh.GetElementTransformation(local_elem, &Tloc);
   Vector lc(3);
   Tloc.Transform(cref, lc);
   const real_t local_cx = lc(0);
   // P1-005: the mirror peer_cx AND reusing the Elem1 reference centroid `cref` as
   // the peer reference IP are valid ONLY because the mesh is a row of unit hexes
   // symmetric about the seam (peer centroid == 2*seam - local; both ref centroids
   // == (0.5,0.5,0.5)).  The expectations are still computed INDEPENDENTLY of the
   // accessor (the analytic profile / per-element step at the physical x).
   const real_t peer_cx  = 2.0 * seam_x - local_cx;

   real_t lam_lx, mu_lx, rho_lx, lam_px, mu_px, rho_px;
   expect(local_cx, lam_lx, mu_lx, rho_lx);
   expect(peer_cx,  lam_px, mu_px, rho_px);

   // PEER read (uses ftr->Elem2 / Elem2No).
   real_t lam_n, mu_n, rho_n;
   op.MaterialAtNbr_(ftr, cref, lam_n, mu_n, rho_n);
   const real_t worst_peer = std::max({ RelErr(lam_n, lam_px),
                                        RelErr(mu_n, mu_px),
                                        RelErr(rho_n, rho_px) });
   TEST_ASSERT(worst_peer <= 1e-12,
               std::string(label) + ": MaterialAtNbr_ == PEER material "
               "(worst rel " + std::to_string(worst_peer) + ")");

   // LOCAL read (caller-owned transform; ftr untouched).
   real_t lam_l, mu_l, rho_l;
   op.MaterialAtLocal_(local_elem, Tloc, cref, lam_l, mu_l, rho_l);
   const real_t worst_local = std::max({ RelErr(lam_l, lam_lx),
                                         RelErr(mu_l, mu_lx),
                                         RelErr(rho_l, rho_lx) });
   TEST_ASSERT(worst_local <= 1e-12,
               std::string(label) + ": MaterialAtLocal_ == LOCAL material "
               "(worst rel " + std::to_string(worst_local) + ")");

   // (Phase 2) flux_pool_gridfunction_constructs (+ all modes): the per-element
   // flux pool was built with no abort, one entry per local element.
   TEST_ASSERT(op.GetPerElementMaterial().size()
               == static_cast<size_t>(pmesh.GetNE()),
               std::string(label) + ": per-element flux pool built (size == NE)");

   // (Phase 2) neighbour_matches_accessor_all_modes: ExchangeBiMaterialNeighbours_
   // stored the TRUE peer material (via MaterialAtNbr_), == the independent PEER
   // expectation — NOT the former local-side stub.
   const auto &nbrmap = op.GetSharedFaceNeighbourMaterial();
   auto nit = nbrmap.find(sf_seam);
   TEST_ASSERT(nit != nbrmap.end(),
               std::string(label) + ": shared neighbour-material entry exists");
   if (nit != nbrmap.end())
   {
      const auto &nm = nit->second;
      const real_t wn = std::max({ RelErr(nm[0], lam_px),
                                   RelErr(nm[1], mu_px),
                                   RelErr(nm[2], rho_px) });
      TEST_ASSERT(wn <= 1e-12,
                  std::string(label) + ": shared_face_neighbour_material_ == PEER "
                  "(rel " + std::to_string(wn) + ")");

      // (Phase 2) bulk_*_byteexact: where the material is seam-continuous at the
      // seam (peer expectation == local expectation), the stored neighbour equals
      // the LOCAL element material — bit-identical to the former local-side stub.
      const bool seam_cont =
         (lam_lx == lam_px && mu_lx == mu_px && rho_lx == rho_px);
      if (seam_cont)
      {
         const auto &lm = op.GetPerElementMaterial()[local_elem];
         const real_t wl = std::max({ std::abs(nm[0] - lm[0]),
                                      std::abs(nm[1] - lm[1]),
                                      std::abs(nm[2] - lm[2]) });
         TEST_ASSERT(wl == 0.0,
                     std::string(label) + ": seam-continuous => neighbour == local "
                     "(byte-exact with the former stub; |dl|="
                     + std::to_string(wl) + ")");
      }
   }

   // P1-001 guard: MaterialAtLocal_ must NOT clobber a held ftr->Elem1, even when
   // called for a DIFFERENT local element.  Record ftr->Elem1's physical centroid,
   // call the accessor on a second local element (via its OWN transform), and
   // require ftr->Elem1 still maps to the same point.  Reverting the accessor to an
   // internal Mesh::GetElementTransformation would reconfigure the shared scratch
   // (== ftr->Elem1) to the other element and break this.
   if (pmesh.GetNE() >= 2)
   {
      Vector x_before(3);
      ftr->Elem1->Transform(cref, x_before);
      const int other = (local_elem != 0) ? 0 : 1;
      IsoparametricTransformation Tother;
      pmesh.GetElementTransformation(other, &Tother);
      real_t a, b, c;
      op.MaterialAtLocal_(other, Tother, cref, a, b, c);
      Vector x_after(3);
      ftr->Elem1->Transform(cref, x_after);
      const real_t move = std::sqrt(std::pow(x_after(0) - x_before(0), 2)
                                  + std::pow(x_after(1) - x_before(1), 2)
                                  + std::pow(x_after(2) - x_before(2), 2));
      TEST_ASSERT(move == 0.0,
                  std::string(label) + ": MaterialAtLocal_ leaves a held "
                  "ftr->Elem1 intact (P1-001; |move|="
                  + std::to_string(move) + ")");
   }

   // No-projection purity: for a per-element-constant (jump/step) material, EVERY
   // interior IP reads the owning element's PURE value (no mix — no node-on-jump
   // pollution by construction, since v4 never projects).  ftr stays valid because
   // Tloc is independent of the mesh scratch.
   if (piecewise_const)
   {
      const real_t xs[3] = { 0.10, 0.50, 0.90 };
      real_t worst_pure_local = 0.0, worst_pure_peer = 0.0;
      for (real_t xr : xs)
      {
         IntegrationPoint ip; ip.Set3(xr, 0.5, 0.5);
         real_t a, b, c;
         op.MaterialAtLocal_(local_elem, Tloc, ip, a, b, c);
         worst_pure_local = std::max(worst_pure_local,
            std::max({ RelErr(a, lam_lx), RelErr(b, mu_lx), RelErr(c, rho_lx) }));
         real_t d, e, f;
         op.MaterialAtNbr_(ftr, ip, d, e, f);
         worst_pure_peer = std::max(worst_pure_peer,
            std::max({ RelErr(d, lam_px), RelErr(e, mu_px), RelErr(f, rho_px) }));
      }
      TEST_ASSERT(worst_pure_local <= 1e-12,
                  std::string(label) + ": LOCAL pure per-side value at all "
                  "interior IPs (no mix; rel " + std::to_string(worst_pure_local) + ")");
      TEST_ASSERT(worst_pure_peer <= 1e-12,
                  std::string(label) + ": PEER pure per-side value at all "
                  "interior IPs (no mix; rel " + std::to_string(worst_pure_peer) + ")");
   }

   int g = 0; MPI_Allreduce(&local_ok, &g, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   TEST_ASSERT(g >= 1, std::string(label) + ": some rank owns the seam face");
}

// =========================================================================
// (Phase 2) shared_central_matches_serial: a fault-adjacent CONTRAST face built
// central in SERIAL (2-sided interior) vs at np=2 (on the seam) must produce
// IDENTICAL central half-Jacobians (c[0]=½A_self, c[1]=½A_peer).  Lateral mu(x)
// makes the peer differ from the local, so this would FAIL under the former
// local-side stub; with the real exchange the np=2 build (true peer via
// MaterialAtNbr_) reproduces the serial build (both materials local) to 1e-12.
// Compared at the MATRIX level (not the deposit) on the rank owning the self
// (elem3) side, so it is robust to face-normal orientation.
// =========================================================================
static void RunSharedCentralMatchesSerial()
{
   if (g_rank == 0) { std::cout << "\n-- shared_central_matches_serial --\n"; }

   auto fn = [](const Vector &x) -> real_t
   { return 30.0e9 * (1.0 + 0.2 * x(0)); };          // lateral lambda(x)==mu(x)
   FunctionCoefficient lam_c(fn), mu_c(fn);
   ConstantCoefficient rho_c(2670.0);

   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {};
   const int order = 1;
   const double seam_x = 4.0;   // 5-hex row, fault at x=3 => x=4 is fault-adjacent

   // --- serial reference (each rank builds its own; serial Mesh is cheap) ------
   Mesh sm = BuildHexRowFaultMesh(/*n=*/5, /*fault_x=*/3);
   BimaterialWaveOperator<Mesh> op_s(
      sm, order, MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c), bc);
   op_s.SetMixedFluxMode(MixedFluxMode::Adjacent);
   int fser = -1;
   for (int f = 0; f < sm.GetNumFaces(); ++f)
   {
      auto *ft = sm.GetFaceElementTransformations(f);
      if (!ft || ft->Elem2No < 0) { continue; }
      Array<int> fv; sm.GetFaceVertices(f, fv);
      real_t cx = 0.0;
      for (int i = 0; i < fv.Size(); ++i) { cx += sm.GetVertex(fv[i])[0]; }
      cx /= fv.Size();
      if (std::abs(cx - seam_x) < 1e-9) { fser = f; break; }
   }

   // --- parallel np=2 (seam at x=4) -------------------------------------------
   Mesh pm0 = BuildHexRowFaultMesh(5, 3);
   const int pne = pm0.GetNE();
   std::vector<int> part(pne, 0);
   for (int e = 0; e < pne; ++e)
   {
      Array<int> ev; pm0.GetElementVertices(e, ev);
      real_t cx = 0.0;
      for (int v = 0; v < ev.Size(); ++v) { cx += pm0.GetVertex(ev[v])[0]; }
      cx /= ev.Size();
      part[e] = (cx < seam_x) ? 0 : 1;
   }
   ParMesh ppm(MPI_COMM_WORLD, pm0, part.data());
   pm0.Clear();
   ppm.ExchangeFaceNbrData();
   BimaterialWaveOperator<ParMesh> op_p(
      ppm, order, MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c), bc);
   op_p.SetMixedFluxMode(MixedFluxMode::Adjacent);
   const int sf = FindSeamSharedFace(ppm, seam_x);

   int did = 0; real_t worst = 0.0;
   if (sf >= 0 && fser >= 0)
   {
      FaceElementTransformations *ftr = ppm.GetSharedFaceTransformations(sf);
      Vector cc(3);
      ftr->Elem1->Transform(
         Geometries.GetCenter(ftr->Elem1->GetGeometryType()), cc);
      if (cc(0) < seam_x)   // this rank owns the self (elem3, lower-x) side
      {
         const int mf = ppm.GetSharedFace(sf);
         const auto &cp = op_p.GetPerFaceCentralFlux();
         const auto &cs = op_s.GetPerFaceCentralFlux();
         auto itp = cp.find(mf);
         auto its = cs.find(fser);
         if (itp != cp.end() && its != cs.end())
         {
            did = 1;
            worst = std::max(MatRel(itp->second[0], its->second[0]),
                             MatRel(itp->second[1], its->second[1]));
         }
      }
   }
   int gdid = 0;
   MPI_Allreduce(&did, &gdid, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   TEST_ASSERT(gdid >= 1,
               "shared_central_matches_serial: a rank owns the self seam side "
               "with both central entries present");
   if (did)
   {
      TEST_ASSERT(worst <= 1e-12,
                  "shared_central_matches_serial: np=2 shared central matrices "
                  "== serial interior central matrices (worst rel "
                  + std::to_string(worst) + ")");
   }
}

// =========================================================================
// (Phase 2, P2-T1) shared_upwind_matches_serial: the BULK shared UPWIND Riemann
// matrix (per_face_bimaterial_flux_, side 0) on a fault-adjacent CONTRAST seam
// must equal the serial 2-sided upwind matrix to 1e-12.  This is the path
// TPV31/TPV102 production uses (interior_flux=matrix, mixed_flux=none) — the
// upwind matrices are built at construction, no mixed flux needed.  Lateral mu(x)
// gives a genuine contrast (peer != local) that would have FAILED under the
// former local-side stub; with the real exchange the np=2 shared upwind == the
// serial 2-sided upwind, proving the re-baseline is the CORRECT bi-material flux.
// `mat` is shared by both ops (its coefficients/constants are globally evaluable);
// the caller keeps any Coefficient objects alive for the call's duration.
// =========================================================================
static void RunSharedUpwindMatchesSerial(const char *label,
                                         const MaterialField &mat)
{
   if (g_rank == 0)
   { std::cout << "\n-- shared_upwind_matches_serial: " << label << " --\n"; }

   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {};
   const int order = 1;
   const double seam_x = 4.0;

   // serial reference: 5-hex row, fault x=3 => x=4 is a non-fault interior face.
   Mesh sm = BuildHexRowFaultMesh(/*n=*/5, /*fault_x=*/3);
   BimaterialWaveOperator<Mesh> op_s(sm, order, mat, bc);   // upwind built in ctor
   int fser = -1;
   for (int f = 0; f < sm.GetNumFaces(); ++f)
   {
      auto *ft = sm.GetFaceElementTransformations(f);
      if (!ft || ft->Elem2No < 0) { continue; }
      Array<int> fv; sm.GetFaceVertices(f, fv);
      real_t cx = 0.0;
      for (int i = 0; i < fv.Size(); ++i) { cx += sm.GetVertex(fv[i])[0]; }
      cx /= fv.Size();
      if (std::abs(cx - seam_x) < 1e-9) { fser = f; break; }
   }

   // parallel np=2, seam at x=4.
   Mesh pm0 = BuildHexRowFaultMesh(5, 3);
   const int pne = pm0.GetNE();
   std::vector<int> part(pne, 0);
   for (int e = 0; e < pne; ++e)
   {
      Array<int> ev; pm0.GetElementVertices(e, ev);
      real_t cx = 0.0;
      for (int v = 0; v < ev.Size(); ++v) { cx += pm0.GetVertex(ev[v])[0]; }
      cx /= ev.Size();
      part[e] = (cx < seam_x) ? 0 : 1;
   }
   ParMesh ppm(MPI_COMM_WORLD, pm0, part.data());
   pm0.Clear();
   ppm.ExchangeFaceNbrData();
   BimaterialWaveOperator<ParMesh> op_p(ppm, order, mat, bc);
   const int sf = FindSeamSharedFace(ppm, seam_x);

   int did = 0; real_t worst = 0.0;
   if (sf >= 0 && fser >= 0)
   {
      FaceElementTransformations *ftr = ppm.GetSharedFaceTransformations(sf);
      Vector cc(3);
      ftr->Elem1->Transform(
         Geometries.GetCenter(ftr->Elem1->GetGeometryType()), cc);
      if (cc(0) < seam_x)   // this rank owns the self (elem3) side
      {
         const int mf = ppm.GetSharedFace(sf);
         const auto &us = op_s.GetPerFaceBimaterialFlux();
         const auto &up = op_p.GetPerFaceBimaterialFlux();
         if (fser < static_cast<int>(us.size())
             && mf < static_cast<int>(up.size()))
         {
            did = 1;
            // side 0 = Elem1 (self): [0][0]=½A_local half, [0][1]=A_nbr half.
            worst = std::max(MatRel(up[mf][0][0], us[fser][0][0]),
                             MatRel(up[mf][0][1], us[fser][0][1]));
         }
      }
   }
   int gdid = 0;
   MPI_Allreduce(&did, &gdid, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   TEST_ASSERT(gdid >= 1,
               std::string("shared_upwind_matches_serial: ") + label
               + ": a rank owns the self seam side with the upwind matrix present");
   if (did)
   {
      TEST_ASSERT(worst <= 1e-12,
                  std::string("shared_upwind_matches_serial: ") + label
                  + ": np=2 shared upwind matrix == serial 2-sided upwind matrix "
                  "(worst rel " + std::to_string(worst) + ")");
   }
}

// =========================================================================
// (Phase 2, P2-T4) coeff_peer_read_independent_of_stale_intpoint: a Coefficient
// that reads T.GetIntPoint() (not the explicit ip) must, through MaterialAtNbr_,
// be evaluated at ip_peer — even if the peer transform carries a stale/wrong int
// point.  Guards the P2-005 SetIntPoint fix: WITHOUT it the read would use the
// poisoned int point and return the wrong material at the seam.
// =========================================================================
static void RunCoeffPeerReadIgnoresStaleIntPoint()
{
   if (g_rank == 0)
   {
      std::cout << "\n-- coeff_peer_read_independent_of_stale_intpoint (P2-T4) --\n";
   }
   IPProbeCoefficient lam_c(30.0e9), mu_c(24.0e9), rho_c(2670.0);
   MaterialField mat = MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c);
   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {};
   const int order = 1;
   const double seam_x = 2.0;

   Mesh sm = BuildHexRowFaultMesh(/*n=*/4, /*fault_x=*/1);
   const int ne = sm.GetNE();
   std::vector<int> part(ne, 0);
   for (int e = 0; e < ne; ++e)
   {
      Array<int> ev; sm.GetElementVertices(e, ev);
      real_t cx = 0.0;
      for (int v = 0; v < ev.Size(); ++v) { cx += sm.GetVertex(ev[v])[0]; }
      cx /= ev.Size();
      part[e] = (cx < seam_x) ? 0 : 1;
   }
   ParMesh ppm(MPI_COMM_WORLD, sm, part.data());
   sm.Clear();
   ppm.ExchangeFaceNbrData();
   TestableBimat<ParMesh> op(ppm, order, mat, bc);
   const int sf = FindSeamSharedFace(ppm, seam_x);

   int did = 0; real_t worst = 0.0;
   if (sf >= 0)
   {
      FaceElementTransformations *ftr = ppm.GetSharedFaceTransformations(sf);
      if (ftr && ftr->Elem2)
      {
         const IntegrationPoint &cref =
            Geometries.GetCenter(ftr->Elem2->GetGeometryType());
         Vector pc(3);
         ftr->Elem2->Transform(cref, pc);             // peer centroid physical
         const real_t exp_lam = 30.0e9 * (1.0 + 0.1 * pc(0));
         IntegrationPoint wrong; wrong.Set3(0.0, 0.0, 0.0);
         ftr->Elem2->SetIntPoint(&wrong);             // poison the peer int point
         real_t lam, mu, rho;
         op.MaterialAtNbr_(ftr, cref, lam, mu, rho);
         worst = RelErr(lam, exp_lam);
         did = 1;
      }
   }
   int gdid = 0;
   MPI_Allreduce(&did, &gdid, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   TEST_ASSERT(gdid >= 1,
               "coeff_peer_read: a rank owns the seam peer side");
   if (did)
   {
      TEST_ASSERT(worst <= 1e-12,
                  "coeff_peer_read_independent_of_stale_intpoint: MaterialAtNbr_ "
                  "evaluates at ip_peer, not the poisoned int point (rel "
                  + std::to_string(worst) + ")");
   }
}

// =========================================================================
// (Phase 3) shared_contrast_reclassify_and_impl8 / shared_weak_contrast_stays_central.
// A halfspace modulus jump AT the seam (x=4) makes the fault-adjacent SHARED
// corridor face a bi-material contrast (now correct via Phase 2's true peer).  With
// the contrast guard ON (tol>=0):
//   - STRONG jump  => the shared corridor face is reclassified OUT of
//     central_flux_face_set_ (dispatches dissipative upwind), and IMPL-8
//     (per_face_central_flux_.size() == central_flux_face_set_.size()) holds.
//   - WEAK jump    => it STAYS in the central set.
// A baseline operator with the guard OFF establishes that the face IS a central
// corridor face on this rank (only the rank whose fault-adjacent element abuts the
// seam owns it), so the check runs on the right rank.
// =========================================================================
static void RunSharedContrastReclassify(const char *label, real_t mod_hi,
                                        real_t tol, bool expect_reclassify)
{
   if (g_rank == 0) { std::cout << "\n-- " << label << " --\n"; }

   const double seam_x = 4.0;
   const real_t mod_lo = 30.0e9;
   auto mod_fn = [=](const Vector &x) -> real_t
   { return (x(0) < seam_x) ? mod_lo : mod_hi; };   // lam = mu = modulus(x)
   FunctionCoefficient lam_c(mod_fn), mu_c(mod_fn);
   ConstantCoefficient rho_c(2670.0);
   MaterialField mat = MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c);

   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {};
   const int order = 1;

   Mesh sm = BuildHexRowFaultMesh(/*n=*/5, /*fault_x=*/3);
   const int ne = sm.GetNE();
   std::vector<int> part(ne, 0);
   for (int e = 0; e < ne; ++e)
   {
      Array<int> ev; sm.GetElementVertices(e, ev);
      real_t cx = 0.0;
      for (int v = 0; v < ev.Size(); ++v) { cx += sm.GetVertex(ev[v])[0]; }
      cx /= ev.Size();
      part[e] = (cx < seam_x) ? 0 : 1;
   }
   ParMesh ppm(MPI_COMM_WORLD, sm, part.data());
   sm.Clear();
   ppm.ExchangeFaceNbrData();
   const int sf = FindSeamSharedFace(ppm, seam_x);
   const int mf = (sf >= 0) ? ppm.GetSharedFace(sf) : -1;

   // Baseline (guard OFF): does this rank own the shared corridor central face?
   BimaterialWaveOperator<ParMesh> op_off(ppm, order, mat, bc);
   op_off.SetMixedFluxMode(MixedFluxMode::Adjacent);
   const bool baseline_central =
      (mf >= 0) && (op_off.GetCentralFluxFaceSet().count(mf) > 0);

   // Guarded operator.
   BimaterialWaveOperator<ParMesh> op_on(ppm, order, mat, bc);
   op_on.SetMixedFluxContrastTol(tol);
   op_on.SetMixedFluxMode(MixedFluxMode::Adjacent);

   // IMPL-8 holds on every rank regardless of ownership.
   TEST_ASSERT(op_on.GetPerFaceCentralFlux().size()
               == op_on.GetCentralFluxFaceSet().size(),
               std::string(label) + ": IMPL-8 (per_face_central_flux_.size() == "
               "central_flux_face_set_.size())");

   const bool on_central =
      (mf >= 0) && (op_on.GetCentralFluxFaceSet().count(mf) > 0);

   int did = baseline_central ? 1 : 0;
   int gdid = 0;
   MPI_Allreduce(&did, &gdid, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   TEST_ASSERT(gdid >= 1,
               std::string(label) + ": a rank owns the shared corridor central face");
   if (did)
   {
      if (expect_reclassify)
      {
         TEST_ASSERT(!on_central,
                     std::string(label) + ": STRONG-contrast shared corridor face "
                     "reclassified OUT of the central set");
         // (P3-002) assert the per-rank shared reclassify COUNT, not just the
         // effect: this mesh has exactly ONE shared corridor face on the owning
         // rank, so the guard reclassifies exactly one.
         TEST_ASSERT(op_on.GetNReclassShared() == 1,
                     std::string(label) + ": GetNReclassShared() == 1 (per-rank "
                     "shared reclassify count; got "
                     + std::to_string(op_on.GetNReclassShared()) + ")");
      }
      else
      {
         TEST_ASSERT(on_central,
                     std::string(label) + ": WEAK-contrast shared corridor face "
                     "STAYS in the central set");
         TEST_ASSERT(op_on.GetNReclassShared() == 0,
                     std::string(label) + ": GetNReclassShared() == 0 (no shared "
                     "reclassify for a weak contrast; got "
                     + std::to_string(op_on.GetNReclassShared()) + ")");
      }
   }
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
         std::cerr << "[cross-rank P1] requires np=2, got " << nprocs
                   << " — TEST FAILED.  Run via `mpirun -np 2 "
                      "./seas_test_bimaterial_seam_material_np2`.\n";
      }
      MPI_Finalize();
      return 1;
   }

   if (g_rank == 0)
   {
      std::cout << "\n=== Cross-rank Phase 1: uniform material accessor "
                   "(MaterialAtLocal_ / MaterialAtNbr_) np=2 ===\n";
   }

   // 4-hex row [0,1],[1,2],[2,3],[3,4]; fault at x=1 (interior to rank 0).
   // EXPLICIT partition at seam_x=2: elems with centroid x<2 -> rank 0 {0,1},
   // else rank 1 {2,3}.  The non-fault interior face at x=2 (between elems 1 and
   // 2) straddles the seam.  Local seam elem: rank0 -> elem 1 (cx=1.5), rank1 ->
   // elem 2 (cx=2.5).  BOTH ranks own 2 elements, so the P1-001 no-clobber guard
   // (which needs a SECOND local element distinct from ftr->Elem1) runs on both.
   const double seam_x = 2.0;
   Mesh smesh = BuildHexRowFaultMesh(/*n=*/4, /*fault_x=*/1);
   const int ne = smesh.GetNE();
   std::vector<int> part(ne, 0);
   for (int e = 0; e < ne; ++e)
   {
      Array<int> ev; smesh.GetElementVertices(e, ev);
      real_t cx = 0.0;
      for (int v = 0; v < ev.Size(); ++v) { cx += smesh.GetVertex(ev[v])[0]; }
      cx /= ev.Size();
      part[e] = (cx < seam_x) ? 0 : 1;
   }
   ParMesh pmesh(MPI_COMM_WORLD, smesh, part.data());
   smesh.Clear();

   // GetSharedFaceTransformations needs the face-neighbour data populated; the
   // operator ctor does this too (idempotent), but FindSeamSharedFace runs first.
   pmesh.ExchangeFaceNbrData();
   const int sf_seam = FindSeamSharedFace(pmesh, seam_x);

   // --- Constant mode ---------------------------------------------------------
   {
      const real_t L = 32.04e9, M = 32.04e9, R = 2670.0;
      auto expect = [=](real_t /*cx*/, real_t &l, real_t &m, real_t &r)
      { l = L; m = M; r = R; };
      RunModeCase("Constant", pmesh, sf_seam, seam_x,
                  MaterialField::MakeConstant(L, M, R),
                  /*piecewise_const=*/false, /*expect_gf_exch=*/0, expect);
   }

   // --- Coefficient (smooth, lateral lambda(x)/mu(x)) -------------------------
   {
      auto lam_fn = [](const Vector &x) -> real_t
      { return 30.0e9 * (1.0 + 0.2 * x(0)); };
      auto mu_fn = [](const Vector &x) -> real_t
      { return 24.0e9 * (1.0 + 0.1 * x(0)); };
      FunctionCoefficient lam_c(lam_fn), mu_c(mu_fn);
      ConstantCoefficient rho_c(2670.0);
      auto expect = [](real_t cx, real_t &l, real_t &m, real_t &r)
      { l = 30.0e9 * (1.0 + 0.2 * cx); m = 24.0e9 * (1.0 + 0.1 * cx); r = 2670.0; };
      RunModeCase("Coefficient (smooth lateral)", pmesh, sf_seam, seam_x,
                  MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c),
                  /*piecewise_const=*/false, /*expect_gf_exch=*/0, expect);
   }

   // --- Coefficient (halfspace JUMP at the seam: pure per side) ---------------
   {
      const real_t LA = 30.0e9, MA = 24.0e9, LB = 50.0e9, MB = 40.0e9;
      auto lam_fn = [=](const Vector &x) -> real_t
      { return (x(0) < seam_x) ? LA : LB; };
      auto mu_fn = [=](const Vector &x) -> real_t
      { return (x(0) < seam_x) ? MA : MB; };
      FunctionCoefficient lam_c(lam_fn), mu_c(mu_fn);
      ConstantCoefficient rho_c(2670.0);
      auto expect = [=](real_t cx, real_t &l, real_t &m, real_t &r)
      { l = (cx < seam_x) ? LA : LB; m = (cx < seam_x) ? MA : MB; r = 2670.0; };
      RunModeCase("Coefficient (halfspace jump)", pmesh, sf_seam, seam_x,
                  MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c),
                  /*piecewise_const=*/true, /*expect_gf_exch=*/0, expect);
   }

   // --- Coefficient (halfspace, jump AWAY from the seam => seam-continuous) ----
   // Jump at x=3 (an element boundary): the seam (x=2) elements [1,2] and [2,3]
   // are BOTH on side A, so the peer material == the local material at the seam.
   // The stored neighbour must then be byte-identical to the former local stub
   // (the bulk_*_halfspace_byteexact gate).
   {
      const real_t LA = 30.0e9, MA = 24.0e9, LB = 50.0e9, MB = 40.0e9;
      const double jump_x = 3.0;
      auto lam_fn = [=](const Vector &x) -> real_t
      { return (x(0) < jump_x) ? LA : LB; };
      auto mu_fn = [=](const Vector &x) -> real_t
      { return (x(0) < jump_x) ? MA : MB; };
      FunctionCoefficient lam_c(lam_fn), mu_c(mu_fn);
      ConstantCoefficient rho_c(2670.0);
      auto expect = [=](real_t cx, real_t &l, real_t &m, real_t &r)
      { l = (cx < jump_x) ? LA : LB; m = (cx < jump_x) ? MA : MB; r = 2670.0; };
      RunModeCase("Coefficient (halfspace seam-continuous, byte-exact)",
                  pmesh, sf_seam, seam_x,
                  MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c),
                  /*piecewise_const=*/true, /*expect_gf_exch=*/0, expect);
   }

   // --- GridFunction (per-element step across the seam) ------------------------
   {
      const real_t LA = 30.0e9, MA = 24.0e9, LB = 50.0e9, MB = 40.0e9;
      const real_t RHO = 2670.0;
      auto lam_step = [=](real_t cx) { return (cx < seam_x) ? LA : LB; };
      auto mu_step  = [=](real_t cx) { return (cx < seam_x) ? MA : MB; };

      // Scalar L2 space at the state order (R-402 contract), on the SAME pmesh.
      L2_FECollection mfec(1, pmesh.Dimension(), BasisType::GaussLobatto);
      ParFiniteElementSpace mfes(&pmesh, &mfec);   // vdim = 1
      auto lam_gf = std::make_shared<ParGridFunction>(&mfes);
      auto mu_gf  = std::make_shared<ParGridFunction>(&mfes);
      auto rho_gf = std::make_shared<ParGridFunction>(&mfes);
      for (int e = 0; e < pmesh.GetNE(); ++e)
      {
         ElementTransformation *T = pmesh.GetElementTransformation(e);
         Vector cc(3);
         T->Transform(Geometries.GetCenter(T->GetGeometryType()), cc);
         Array<int> dofs; mfes.GetElementDofs(e, dofs);
         for (int d = 0; d < dofs.Size(); ++d)
         {
            (*lam_gf)(dofs[d]) = lam_step(cc(0));
            (*mu_gf)(dofs[d])  = mu_step(cc(0));
            (*rho_gf)(dofs[d]) = RHO;
         }
      }
      MaterialField mat = MaterialField::MakeGridFunction(rho_gf, lam_gf, mu_gf);
      auto expect = [=](real_t cx, real_t &l, real_t &m, real_t &r)
      { l = lam_step(cx); m = mu_step(cx); r = RHO; };
      RunModeCase("GridFunction (per-element step)", pmesh, sf_seam, seam_x,
                  mat, /*piecewise_const=*/true, /*expect_gf_exch=*/3, expect);
   }

   // Phase 2: serial-vs-parallel central matrices on a fault-adjacent contrast face.
   RunSharedCentralMatchesSerial();

   // Phase 2 (P2-T1): serial-vs-parallel BULK UPWIND matrices — the production
   // path (interior_flux=matrix, mixed_flux=none).  Coefficient contrast + a
   // Constant byte-exact control.
   {
      auto fn = [](const Vector &x) -> real_t
      { return 30.0e9 * (1.0 + 0.2 * x(0)); };
      FunctionCoefficient lam_c(fn), mu_c(fn);
      ConstantCoefficient rho_c(2670.0);
      MaterialField mat = MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c);
      RunSharedUpwindMatchesSerial("Coefficient mu(x) contrast", mat);
   }
   {
      MaterialField mat = MaterialField::MakeConstant(32.04e9, 32.04e9, 2670.0);
      RunSharedUpwindMatchesSerial("Constant (byte-exact)", mat);
   }

   // Phase 2 (P2-T4): peer Coefficient read must use ip_peer, not a stale int point.
   RunCoeffPeerReadIgnoresStaleIntPoint();

   // Phase 3: shared-face contrast guard (reclassify strong; keep weak) + IMPL-8.
   RunSharedContrastReclassify(
      "shared_contrast_reclassify_and_impl8 (strong)",
      /*mod_hi=*/120.0e9, /*tol=*/0.1, /*expect_reclassify=*/true);
   RunSharedContrastReclassify(
      "shared_weak_contrast_stays_central (weak)",
      /*mod_hi=*/33.0e9, /*tol=*/0.1, /*expect_reclassify=*/false);

   int total = 0, passed = 0, failed = 0;
   MPI_Allreduce(&g_num_tests,  &total,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
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

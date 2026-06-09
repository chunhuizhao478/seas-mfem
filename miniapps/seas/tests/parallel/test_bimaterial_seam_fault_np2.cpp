// Cross-rank material exchange — Phase 4 gates (PLAN_cross_rank_material_exchange
// _2026-06-06.md).  np=2 coverage of the per-side FAULT impedance assignment on a
// SHARED (cross-rank) fault face: BimaterialWaveOperator::
// AssignFaultSidePerMaterialImpedances reads the LOCAL side via MaterialAtLocal_
// and the PEER side across the seam via MaterialAtNbr_, BOTH at the eps-offset
// fault QP (in2 = -in1 ⇒ same fault QP, opposite normal ⇒ exact symmetry).
//
//   shared_fault_constant_byteexact:   Constant  ⇒ Zp_plus == Zp_minus == const.
//   shared_fault_symmetric_depthprofile: depth f(z) on an ASYMMETRIC (skewed-peer)
//                                        mesh ⇒ |Zp_plus - Zp_minus| <= 1e-10*Zp
//                                        (the eps-offset is along the fault normal,
//                                        z unchanged, even with a skewed peer — the
//                                        R-401/R-302 peer-Jacobian check).
//   shared_fault_contrast:             halfspace across the fault ⇒ Zp_plus !=
//                                        Zp_minus, and the {min,max} pair matches
//                                        the SERIAL interior-fault values to 1e-12.
//   shared_fault_gridfunction:         per-element step across the fault ⇒ correct
//                                        per-side contrast (NEW; legacy aborted).
//   shared_fault_index_mapping:        TWO shared fault faces with DISTINCT
//                                        contrasts ⇒ each DOF carries ITS face's
//                                        values (R-204 no index transposition).
//
// A separate file from test_bimaterial_seam_material_np2.cpp (the bulk-material
// gates) because the fault gates need the fault-DOF setup + shared-fault meshes;
// both are the plan's Phase-4 acceptance gates.
//
// Run: mpirun -np 2 ./seas_test_bimaterial_seam_fault_np2

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/bimaterial_wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/heterogeneous_material.hpp"
#include "../../dynamic/fault_face_flux.hpp"
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

real_t RelErr(real_t a, real_t b)
{
   const real_t s = std::max(std::abs(a), std::abs(b));
   return (s == 0.0) ? std::abs(a - b) : std::abs(a - b) / s;
}

// Two tets sharing the triangular fault face v0v1v2 (the y=0 plane).  `skew`
// displaces the +y apex (v3) so the +y (peer) tet has a DIFFERENT shape than the
// -y tet — an asymmetric mesh that exercises the peer eps-offset Jacobian (R-401).
Mesh BuildTwoTetSharedFault(bool skew)
{
   Mesh mesh(3, 5, 2, 0);
   // `skew` shifts the +y apex (v3) in x/z (NOT y) so the +y (peer) tet has a
   // different shape than the -y tet — the fault face v0v1v2 stays in the y=0 plane.
   real_t verts[5][3] = {
      {0.0, 0.0, 0.0},                                  // v0 (fault, y=0)
      {1.0, 0.0, 0.0},                                  // v1 (fault)
      {0.0, 0.0, 1.0},                                  // v2 (fault)
      {skew ? 0.4 : 0.0, 1.0, skew ? 0.3 : 0.0},        // v3 apex (+y)
      {0.0, -1.0, 0.0},                                 // v4 apex (-y)
   };
   for (int v = 0; v < 5; ++v) { mesh.AddVertex(verts[v]); }
   mesh.AddTet(0, 1, 2, 4, 1);   // elem 0: y<0
   mesh.AddTet(0, 1, 2, 3, 1);   // elem 1: y>0
   mesh.FinalizeTopology();
   for (int f = 0; f < mesh.GetNumFaces(); ++f)
   {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      if (fv.Size() != 3) { continue; }
      auto *ftr = mesh.GetFaceElementTransformations(f);
      const int attr = (ftr && ftr->Elem2No >= 0) ? 3 : 1;  // interior => fault
      mesh.AddBdrTriangle(fv[0], fv[1], fv[2], attr);
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// 4 unit hexes: x columns [0,1] and [1,2]; y halves [-1,0] and [0,1]; z [0,1].
// The y=0 plane carries TWO interior faces (one per x column), both tagged FAULT.
Mesh BuildFourHexTwoFault()
{
   const real_t xs[3] = {0.0, 1.0, 2.0};
   const real_t ys[3] = {-1.0, 0.0, 1.0};
   const real_t zs[2] = {0.0, 1.0};
   Mesh mesh(3, 3 * 3 * 2, 4, 0);
   auto vid = [](int ix, int iy, int iz) { return (ix * 3 + iy) * 2 + iz; };
   for (int ix = 0; ix < 3; ++ix)
      for (int iy = 0; iy < 3; ++iy)
         for (int iz = 0; iz < 2; ++iz)
         {
            real_t v[3] = { xs[ix], ys[iy], zs[iz] };
            mesh.AddVertex(v);
         }
   for (int ix = 0; ix < 2; ++ix)
      for (int iy = 0; iy < 2; ++iy)
      {
         const int v[8] = {
            vid(ix,   iy,   0), vid(ix+1, iy,   0), vid(ix+1, iy+1, 0), vid(ix, iy+1, 0),
            vid(ix,   iy,   1), vid(ix+1, iy,   1), vid(ix+1, iy+1, 1), vid(ix, iy+1, 1)
         };
         mesh.AddHex(v, 1);
      }
   mesh.FinalizeTopology();
   for (int f = 0; f < mesh.GetNumFaces(); ++f)
   {
      Array<int> fv; mesh.GetFaceVertices(f, fv);
      real_t cy = 0.0;
      for (int i = 0; i < fv.Size(); ++i) { cy += mesh.GetVertex(fv[i])[1]; }
      cy /= fv.Size();
      auto *ftr = mesh.GetFaceElementTransformations(f);
      if (ftr && ftr->Elem2No >= 0)
      {
         if (std::abs(cy) < 1e-9)   // y=0 interior face => FAULT
         { mesh.AddBdrQuad(fv[0], fv[1], fv[2], fv[3], 3); }
         continue;
      }
      mesh.AddBdrQuad(fv[0], fv[1], fv[2], fv[3], 1);   // FREE (external)
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

}  // namespace

#ifdef MFEM_USE_MPI
// Partition elements by centroid y-sign (y<0 -> rank 0, else rank 1) so every
// y=0 fault face straddles the seam.
static ParMesh PartitionByYSign(Mesh &&smesh_in)
{
   Mesh smesh(std::move(smesh_in));
   const int ne = smesh.GetNE();
   std::vector<int> part(ne, 0);
   for (int e = 0; e < ne; ++e)
   {
      Array<int> ev; smesh.GetElementVertices(e, ev);
      real_t cy = 0.0;
      for (int v = 0; v < ev.Size(); ++v) { cy += smesh.GetVertex(ev[v])[1]; }
      cy /= ev.Size();
      part[e] = (cy < 0.0) ? 0 : 1;
   }
   ParMesh pmesh(MPI_COMM_WORLD, smesh, part.data());
   smesh.Clear();
   return pmesh;
}

// Construct a parallel bimaterial operator, set up the fault DOFData, run
// AssignFaultSidePerMaterialImpedances, and return the dof_data (by value).
// out_nqp/out_nint/out_nshr report the fault layout.
static std::vector<DOFData> RunAssignFaultParallel(
   ParMesh &pmesh, const MaterialField &mat, int order,
   int &out_nqp, int &out_nint, int &out_nshr)
{
   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {};
   BimaterialWaveOperator<ParMesh> op(pmesh, order, mat, bc);
   out_nint = op.GetFaultInteriorFaces().Size();
   out_nshr = op.GetFaultSharedFaces().Size();
   int nqp = 0;
   if (out_nshr > 0)
   {
      auto *f = pmesh.GetSharedFaceTransformations(op.GetFaultSharedFaces()[0]);
      nqp = IntRules.Get(f->GetGeometryType(), 2 * order).GetNPoints();
   }
   else if (out_nint > 0)
   {
      auto *f = pmesh.GetInteriorFaceTransformations(op.GetFaultInteriorFaces()[0]);
      nqp = IntRules.Get(f->GetGeometryType(), 2 * order).GetNPoints();
   }
   int gnqp = nqp;
   MPI_Allreduce(&nqp, &gnqp, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   nqp = gnqp;
   out_nqp = nqp;
   const int ntot = (out_nint + out_nshr) * nqp;
   std::vector<DOFData> dof(ntot > 0 ? ntot : 0);
   if (ntot > 0) { op.SetFaultDOFData(&dof, nqp); }
   op.AssignFaultSidePerMaterialImpedances(dof);
   return dof;
}

// Serial variant (both elements local; fault is an interior face).
static std::vector<DOFData> RunAssignFaultSerial(
   Mesh &smesh, const MaterialField &mat, int order, int &out_nqp, int &out_nint)
{
   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {};
   BimaterialWaveOperator<Mesh> op(smesh, order, mat, bc);
   out_nint = op.GetFaultInteriorFaces().Size();
   int nqp = 0;
   if (out_nint > 0)
   {
      auto *f = smesh.GetInteriorFaceTransformations(op.GetFaultInteriorFaces()[0]);
      nqp = IntRules.Get(f->GetGeometryType(), 2 * order).GetNPoints();
   }
   out_nqp = nqp;
   const int ntot = out_nint * nqp;
   std::vector<DOFData> dof(ntot > 0 ? ntot : 0);
   if (ntot > 0) { op.SetFaultDOFData(&dof, nqp); }
   op.AssignFaultSidePerMaterialImpedances(dof);
   return dof;
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
         std::cerr << "[cross-rank P4] requires np=2, got " << nprocs
                   << " — TEST FAILED.\n";
      }
      MPI_Finalize();
      return 1;
   }
   if (g_rank == 0)
   {
      std::cout << "\n=== Cross-rank Phase 4: per-side FAULT material across a "
                   "seam (np=2) ===\n";
   }
   const int order = 1;

   // --- shared_fault_constant_byteexact --------------------------------------
   {
      if (g_rank == 0) { std::cout << "\n-- shared_fault_constant_byteexact --\n"; }
      const real_t L = 32.04e9, M = 32.04e9, R = 2670.0;
      const real_t Zp_const = std::sqrt((L + 2.0 * M) * R);
      ParMesh pmesh = PartitionByYSign(BuildTwoTetSharedFault(/*skew=*/false));
      int nqp = 0, nint = 0, nshr = 0;
      std::vector<DOFData> dof = RunAssignFaultParallel(
         pmesh, MaterialField::MakeConstant(L, M, R), order, nqp, nint, nshr);
      int local_shr = nshr, glob_shr = 0;
      MPI_Allreduce(&local_shr, &glob_shr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      TEST_ASSERT(glob_shr >= 1, "constant: a rank owns the shared fault face");
      if (nshr > 0)   // (P4-002) non-vacuous: only assert on a rank that owns one
      {
         real_t worst = 0.0;
         int n_checked = 0;
         for (int i = nint * nqp; i < static_cast<int>(dof.size()); ++i)
         {
            worst = std::max(worst, RelErr(dof[i].Zp_plus, Zp_const));
            worst = std::max(worst, RelErr(dof[i].Zp_minus, Zp_const));
            ++n_checked;
         }
         TEST_ASSERT(n_checked > 0 && worst <= 1e-12,
                     std::string("constant: shared fault Zp_plus == Zp_minus == "
                     "const (worst rel ") + std::to_string(worst) + ")");
      }
   }

   // --- shared_fault_symmetric_depthprofile (asymmetric/skewed peer) ---------
   {
      if (g_rank == 0)
      { std::cout << "\n-- shared_fault_symmetric_depthprofile (skewed peer) --\n"; }
      auto fz = [](const Vector &x) -> real_t
      { return 30.0e9 * (1.0 + 0.2 * x(2)); };   // depth profile in z
      FunctionCoefficient lam_c(fz), mu_c(fz);
      ConstantCoefficient rho_c(2670.0);
      ParMesh pmesh = PartitionByYSign(BuildTwoTetSharedFault(/*skew=*/true));
      int nqp = 0, nint = 0, nshr = 0;
      std::vector<DOFData> dof = RunAssignFaultParallel(
         pmesh, MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c),
         order, nqp, nint, nshr);
      int local_shr = nshr, glob_shr = 0;
      MPI_Allreduce(&local_shr, &glob_shr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      TEST_ASSERT(glob_shr >= 1, "depthprofile: a rank owns the shared fault face");
      if (nshr > 0)   // (P4-002) non-vacuous
      {
         real_t worst = 0.0;
         int n_checked = 0;
         for (int i = nint * nqp; i < static_cast<int>(dof.size()); ++i)
         {
            worst = std::max(worst, RelErr(dof[i].Zp_plus, dof[i].Zp_minus));
            ++n_checked;
         }
         TEST_ASSERT(n_checked > 0 && worst <= 1e-10,
                     std::string("depthprofile: shared fault Zp_plus == Zp_minus on "
                     "the asymmetric mesh (eps-offset along the normal; worst rel ")
                     + std::to_string(worst) + ")");
      }
   }

   // --- shared_fault_contrast (halfspace; matches serial interior) -----------
   {
      if (g_rank == 0) { std::cout << "\n-- shared_fault_contrast --\n"; }
      auto half = [](const Vector &x) -> real_t
      { return (x(1) >= 0.0) ? 120.0e9 : 30.0e9; };   // strong on +y, weak on -y
      FunctionCoefficient lam_c(half), mu_c(half);
      ConstantCoefficient rho_c(2670.0);
      // parallel (shared fault)
      ParMesh pmesh = PartitionByYSign(BuildTwoTetSharedFault(/*skew=*/false));
      int nqp = 0, nint = 0, nshr = 0;
      std::vector<DOFData> dof = RunAssignFaultParallel(
         pmesh, MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c),
         order, nqp, nint, nshr);
      // serial (interior fault) — same mesh, both tets local.
      FunctionCoefficient lam_s(half), mu_s(half);
      ConstantCoefficient rho_s(2670.0);
      Mesh smesh = BuildTwoTetSharedFault(/*skew=*/false);
      int snqp = 0, snint = 0;
      std::vector<DOFData> sdof = RunAssignFaultSerial(
         smesh, MaterialField::MakeCoefficient(&lam_s, &mu_s, &rho_s),
         order, snqp, snint);

      int local_shr = nshr, glob_shr = 0;
      MPI_Allreduce(&local_shr, &glob_shr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      TEST_ASSERT(glob_shr >= 1, "contrast: a rank owns the shared fault face");

      // halfspace is piecewise-constant per side => every QP has the same pair.
      const real_t strongZp = std::sqrt((120e9 + 2 * 120e9) * 2670.0);
      const real_t weakZp   = std::sqrt((30e9 + 2 * 30e9) * 2670.0);
      // serial pair (sorted) must equal the analytic weak/strong — a sanity that
      // the serial ORACLE itself is right before using it as the +/- reference.
      real_t s_lo = 0.0, s_hi = 0.0;
      if (!sdof.empty())
      {
         s_lo = std::min(sdof[0].Zp_plus, sdof[0].Zp_minus);
         s_hi = std::max(sdof[0].Zp_plus, sdof[0].Zp_minus);
      }
      TEST_ASSERT(!sdof.empty()
                  && RelErr(s_lo, weakZp) <= 1e-12
                  && RelErr(s_hi, strongZp) <= 1e-12,
                  "contrast: serial interior fault recovers weak/strong Zp");

      // (P4-001) PLUS/MINUS-SENSITIVE: compare the parallel shared fault's
      // Zp_plus/Zp_minus to the SERIAL interior oracle's UNSORTED pair (same
      // geometry => same FaultBasis +/- convention).  A plus/minus transposition
      // in the shared loop's e1_plus swap would make parallel Zp_plus match the
      // serial Zp_MINUS => this fails (the old sorted {min,max} compare could not
      // catch it).  Also assert a real contrast (Zp_plus != Zp_minus).
      const real_t s_zp_plus  = sdof.empty() ? 0.0 : sdof[0].Zp_plus;
      const real_t s_zp_minus = sdof.empty() ? 0.0 : sdof[0].Zp_minus;
      if (nshr > 0)
      {
         real_t worst_pm = 0.0;
         bool all_contrast = true;
         int n_checked = 0;
         for (int i = nint * nqp; i < static_cast<int>(dof.size()); ++i)
         {
            worst_pm = std::max(worst_pm,
               std::max(RelErr(dof[i].Zp_plus,  s_zp_plus),
                        RelErr(dof[i].Zp_minus, s_zp_minus)));
            if (RelErr(dof[i].Zp_plus, dof[i].Zp_minus) <= 1e-6)
            { all_contrast = false; }
            ++n_checked;
         }
         TEST_ASSERT(!sdof.empty() && n_checked > 0 && worst_pm <= 1e-12,
                     std::string("contrast: shared fault Zp_plus/Zp_minus == serial "
                     "interior UNSORTED (catches +/- transposition; worst rel ")
                     + std::to_string(worst_pm) + ")");
         TEST_ASSERT(all_contrast,
                     "contrast: shared fault Zp_plus != Zp_minus (real contrast)");
      }
   }

   // --- shared_fault_gridfunction (per-element step across the fault) ---------
   {
      if (g_rank == 0) { std::cout << "\n-- shared_fault_gridfunction --\n"; }
      const real_t MOD_WEAK = 30.0e9, MOD_STRONG = 120.0e9, RHO = 2670.0;
      const real_t weakZp   = std::sqrt(3.0 * MOD_WEAK * RHO);
      const real_t strongZp = std::sqrt(3.0 * MOD_STRONG * RHO);
      ParMesh pmesh = PartitionByYSign(BuildTwoTetSharedFault(/*skew=*/false));
      L2_FECollection mfec(order, pmesh.Dimension(), BasisType::GaussLobatto);
      ParFiniteElementSpace mfes(&pmesh, &mfec);
      auto lam_gf = std::make_shared<ParGridFunction>(&mfes);
      auto mu_gf  = std::make_shared<ParGridFunction>(&mfes);
      auto rho_gf = std::make_shared<ParGridFunction>(&mfes);
      for (int e = 0; e < pmesh.GetNE(); ++e)
      {
         ElementTransformation *T = pmesh.GetElementTransformation(e);
         Vector cc(3);
         T->Transform(Geometries.GetCenter(T->GetGeometryType()), cc);
         const real_t mod = (cc(1) >= 0.0) ? MOD_STRONG : MOD_WEAK;  // +y strong
         Array<int> dofs; mfes.GetElementDofs(e, dofs);
         for (int d = 0; d < dofs.Size(); ++d)
         {
            (*lam_gf)(dofs[d]) = mod;
            (*mu_gf)(dofs[d])  = mod;
            (*rho_gf)(dofs[d]) = RHO;
         }
      }
      MaterialField mat = MaterialField::MakeGridFunction(rho_gf, lam_gf, mu_gf);
      int nqp = 0, nint = 0, nshr = 0;
      std::vector<DOFData> dof =
         RunAssignFaultParallel(pmesh, mat, order, nqp, nint, nshr);
      int local_shr = nshr, glob_shr = 0;
      MPI_Allreduce(&local_shr, &glob_shr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      TEST_ASSERT(glob_shr >= 1, "gridfunction: a rank owns the shared fault face");
      if (nshr > 0)
      {
         real_t worst_match = 0.0;
         bool all_contrast = true;
         for (int i = nint * nqp; i < static_cast<int>(dof.size()); ++i)
         {
            const real_t lo = std::min(dof[i].Zp_plus, dof[i].Zp_minus);
            const real_t hi = std::max(dof[i].Zp_plus, dof[i].Zp_minus);
            worst_match = std::max(worst_match,
                                   std::max(RelErr(lo, weakZp), RelErr(hi, strongZp)));
            if (RelErr(dof[i].Zp_plus, dof[i].Zp_minus) <= 1e-6)
            { all_contrast = false; }
         }
         TEST_ASSERT(worst_match <= 1e-12 && all_contrast,
                     std::string("gridfunction: shared fault per-side contrast "
                     "(weak/strong) correct (worst rel ")
                     + std::to_string(worst_match) + ")");
      }
   }

   // --- shared_fault_index_mapping (2 faces, distinct contrasts) -------------
   {
      if (g_rank == 0) { std::cout << "\n-- shared_fault_index_mapping --\n"; }
      // Contrast differs by x-column: weak on -y always; strong on +y is bigger
      // for the x>1 column.  So fault face in column [0,1] has a SMALLER contrast
      // than the face in column [1,2]; each fault DOF must carry ITS face's pair.
      auto modf = [](const Vector &x) -> real_t
      {
         const real_t base = (x(1) >= 0.0) ? 1.0 : 0.0;     // +y strong, -y weak
         const real_t col  = (x(0) >= 1.0) ? 2.0 : 1.0;     // x-column scale
         return 30.0e9 * (1.0 + base * col);                // -y:30e9; +y:30*(1+col)
      };
      FunctionCoefficient lam_c(modf), mu_c(modf);
      ConstantCoefficient rho_c(2670.0);
      ParMesh pmesh = PartitionByYSign(BuildFourHexTwoFault());
      int nqp = 0, nint = 0, nshr = 0;
      std::vector<DOFData> dof = RunAssignFaultParallel(
         pmesh, MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c),
         order, nqp, nint, nshr);
      int local_shr = nshr, glob_shr = 0;
      MPI_Allreduce(&local_shr, &glob_shr, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      TEST_ASSERT(glob_shr >= 2, "index_mapping: >= 2 shared fault faces globally");

      // For each shared fault face si (DOFs at (nint+si)*nqp .. +nqp), the face's
      // x-centroid selects the expected per-side Zp pair; assert every DOF in the
      // block matches ITS face's pair (no cross-face leak).
      auto &op_pmesh = pmesh;
      const Array<int> *shr_faces = nullptr;   // re-derive via a throwaway op below
      // Re-fetch the shared fault face list + offsets from a fresh operator (the
      // RunAssignFaultParallel op is gone; rebuild to read GetFaultSharedFaces()).
      BoundaryConfig bc; bc.natural_attrs={1}; bc.fault_attr=3; bc.absorbing_attrs={};
      BimaterialWaveOperator<ParMesh> op2(op_pmesh, order,
         MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c), bc);
      const Array<int> &sf_list = op2.GetFaultSharedFaces();
      (void)shr_faces;
      real_t worst = 0.0; int n_blocks = 0;
      for (int si = 0; si < sf_list.Size(); ++si)
      {
         auto *ftr = op_pmesh.GetSharedFaceTransformations(sf_list[si]);
         if (!ftr) { continue; }
         Vector fc(3);
         ftr->Face->Transform(
            Geometries.GetCenter(ftr->Face->GetGeometryType()), fc);
         const real_t col = (fc(0) >= 1.0) ? 2.0 : 1.0;
         const real_t modStrong = 30.0e9 * (1.0 + col);
         const real_t expHi = std::sqrt(3.0 * modStrong * 2670.0);
         const real_t expLo = std::sqrt(3.0 * 30.0e9 * 2670.0);
         const int base = (nint + si) * nqp;
         for (int q = 0; q < nqp; ++q)
         {
            const int idx = base + q;
            if (idx >= static_cast<int>(dof.size())) { continue; }
            const real_t lo = std::min(dof[idx].Zp_plus, dof[idx].Zp_minus);
            const real_t hi = std::max(dof[idx].Zp_plus, dof[idx].Zp_minus);
            worst = std::max(worst, std::max(RelErr(lo, expLo), RelErr(hi, expHi)));
         }
         ++n_blocks;
      }
      if (nshr > 0)
      {
         TEST_ASSERT(n_blocks == nshr && worst <= 1e-12,
                     std::string("index_mapping: each shared fault DOF carries "
                     "ITS face's per-side Zp (worst rel ")
                     + std::to_string(worst) + ")");
      }
   }

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

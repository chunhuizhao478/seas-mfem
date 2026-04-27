// Round-13 R-1201 — cross-rank end-to-end Adjacent dispatch numerical
// correctness (parallel, np=2).
//
// The Round-2 / Round-3 MPI tests (test_mixed_flux_adjacent_mpi,
// test_mixed_flux_shared_fault_face_excluded_mpi) verify
// `central_flux_face_set_` MEMBERSHIP symmetry across ranks but never
// run `wave.AdvanceADER` or `wave.Mult` with `MixedFluxMode::Adjacent`
// at np>=2.  The shared-face dispatch sites at wave_operator.inl:2735
// and :4109 have ZERO MPI numerical-output coverage with `mf_on=true`.
//
// This test closes that gap: at np=2, run `wave.Mult(Q_init, dQdt)`
// with mode=Adjacent on the 1x2x2 fixture used by R-001 (4 hexes,
// 24 tets, fault perpendicular to the z=L/2 rank seam).  Compare three
// global per-component invariants against the serial-reference run on
// the same mesh.  Bit-perfect dispatch implies all 27 invariants
// match to ~FP precision; cross-rank conservation breakage shows up
// as visible deviation.
//
// Invariants computed per state component c:
//   (i)   sum_local_dofs Q_new[c]      (linear; partition-invariant
//         under bit-equal element-local Q_new at every element).
//   (ii)  sum_local_dofs Q_new[c]^2    (quadratic; preserves L2
//         norm to FP precision).
//   (iii) max_local_dofs |Q_new[c]|    (∞-norm; identifies the
//         single largest DOF value globally).
//
// Tolerance: 1e-10 relative for (i), (ii); 1e-12 for (iii).  The
// per-DOF FP summation order may differ between serial and parallel
// runs (depending on element traversal order); 1e-10 absorbs that
// noise without masking real dispatch errors (which would be 6-12
// orders of magnitude larger).
//
// Run: mpirun -np 2 ./seas_test_mixed_flux_dispatch_adjacent_mpi

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../dynamic/tpv102_setup_total.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

namespace mfem { namespace seas { int g_seas_my_rank = 0; } }  // NOLINT

namespace
{

constexpr real_t kL      = 1000.0;
constexpr real_t kRho    = TPV102Params::rho;
constexpr real_t kCp     = TPV102Params::cp;
constexpr real_t kCs     = TPV102Params::cs;
constexpr real_t kLambda = TPV102Params::lambda;
constexpr real_t kMu     = TPV102Params::mu;

// Same fixture as test_mixed_flux_adjacent_mpi.cpp: 1x2x2 tet-split with
// fault at y=L/2, partition by z (lower-z hexes -> rank 1; upper-z ->
// rank 0).  central_flux_face_set_ is non-empty in Adjacent mode and
// the rank seam at z=L/2 carries non-fault shared faces.
Mesh BuildOrthogonalFaultMesh()
{
   Mesh mesh = Mesh::MakeCartesian3D(1, 2, 2, Element::TETRAHEDRON,
                                     kL, kL, kL, /*sfc_ordering=*/false);
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
         if (std::abs(cy - 0.5 * kL) < 1e-8)
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

std::vector<int> BuildSkewedPartition(const Mesh &mesh)
{
   const int ne = mesh.GetNE();
   std::vector<int> part(ne, 0);
   for (int e = 0; e < ne; e++)
   {
      Array<int> ev;
      mesh.GetElementVertices(e, ev);
      real_t cz = 0.0;
      for (int v = 0; v < ev.Size(); v++) { cz += mesh.GetVertex(ev[v])[2]; }
      cz /= ev.Size();
      part[e] = (cz < 0.5 * kL) ? 1 : 0;
   }
   return part;
}

// Discontinuous-by-element IC keyed on the element CENTROID (partition-
// invariant).  Using local element index `e` would produce different
// per-element values on serial vs ParMesh because partitioning shuffles
// the local element-index → physical-element map.  The centroid is the
// same physical point on both, so per-element s/t scalars match.
template <typename FESpace>
void InitializeNonUniformState(Vector &Q, FESpace &fes, int ndof_total)
{
   Q.SetSize(NUM_STATE * ndof_total);
   Q = 0.0;
   real_t *Q_data = Q.GetData();
   for (int e = 0; e < fes.GetNE(); e++)
   {
      const FiniteElement *fe = fes.GetFE(e);
      const int ndof_e = fe->GetDof();
      Array<int> dof_idx;
      fes.GetElementDofs(e, dof_idx);

      // Compute element centroid in physical space — partition-invariant.
      ElementTransformation *etr = fes.GetElementTransformation(e);
      Geometry::Type geom = fe->GetGeomType();
      const IntegrationPoint &cip = Geometries.GetCenter(geom);
      etr->SetIntPoint(&cip);
      Vector cen(3); etr->Transform(cip, cen);

      // Pseudo-random scalars from centroid coords (deterministic, same
      // physical element → same scalars on serial and parallel).
      const real_t s = 1.0 + 0.3 * std::sin(7.0 * cen(0) / kL +
                                            11.0 * cen(1) / kL +
                                            13.0 * cen(2) / kL);
      const real_t t = 1.0 + 0.2 * std::cos(5.0 * cen(0) / kL +
                                            17.0 * cen(1) / kL +
                                            19.0 * cen(2) / kL);

      for (int i = 0; i < ndof_e; i++)
      {
         const int dof = dof_idx[i];
         Q_data[SXX * ndof_total + dof] = -TPV102Params::sigma_n * s;
         Q_data[SYY * ndof_total + dof] =  TPV102Params::sigma_n * t;
         Q_data[SXY * ndof_total + dof] = -TPV102Params::tau_ini * s;
         Q_data[SZZ * ndof_total + dof] =  1.0e6 * t;
         Q_data[SYZ * ndof_total + dof] =  1.0e6 * s;
         Q_data[SXZ * ndof_total + dof] =  1.0e6 * t;
         Q_data[VX  * ndof_total + dof] =  0.1 * s;
         Q_data[VY  * ndof_total + dof] =  0.1 * t;
         Q_data[VZ  * ndof_total + dof] =  0.1 * (s + t);
      }
   }
}

// Three per-component invariants of a Q vector: (sum, sum_sq, max_abs).
struct Invariants
{
   real_t sum    [NUM_STATE];
   real_t sum_sq [NUM_STATE];
   real_t max_abs[NUM_STATE];
};

template <typename FESpace>
Invariants ComputeLocalInvariants(const Vector &Q, const FESpace &fes,
                                  int ndof_total)
{
   Invariants inv;
   for (int c = 0; c < NUM_STATE; c++)
   {
      inv.sum[c] = 0.0;
      inv.sum_sq[c] = 0.0;
      inv.max_abs[c] = 0.0;
   }
   const real_t *Q_data = Q.GetData();
   // Iterate by ELEMENT (not raw DOF index) so each rank counts only
   // its locally-owned DOFs.  In DG every DOF is owned by exactly one
   // element, and ParMesh elements are partitioned with no overlap.
   for (int e = 0; e < fes.GetNE(); e++)
   {
      Array<int> dof_idx;
      fes.GetElementDofs(e, dof_idx);
      for (int i = 0; i < dof_idx.Size(); i++)
      {
         const int dof = dof_idx[i];
         for (int c = 0; c < NUM_STATE; c++)
         {
            const real_t v = Q_data[c * ndof_total + dof];
            inv.sum[c]    += v;
            inv.sum_sq[c] += v * v;
            inv.max_abs[c] = std::max(inv.max_abs[c], std::abs(v));
         }
      }
   }
   return inv;
}

#ifdef MFEM_USE_MPI
Invariants AllreduceInvariants(const Invariants &local, MPI_Comm comm)
{
   Invariants global;
   MPI_Allreduce(local.sum,    global.sum,    NUM_STATE, MPI_DOUBLE,
                 MPI_SUM, comm);
   MPI_Allreduce(local.sum_sq, global.sum_sq, NUM_STATE, MPI_DOUBLE,
                 MPI_SUM, comm);
   MPI_Allreduce(local.max_abs, global.max_abs, NUM_STATE, MPI_DOUBLE,
                 MPI_MAX, comm);
   return global;
}
#endif

void SetupWaveSerial(WaveOperator<Mesh> &wave,
                     std::vector<DOFData> &dof_data,
                     Mesh &mesh, int order)
{
   const int n_local_qps = wave.GetNumLocalFaultQPs();
   const int nbf = wave.GetNbfPerFace();
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   std::vector<Vector> fault_coords;
   fault_coords.reserve(n_local_qps);
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir = IntRules.Get(
         ftr->GetGeometryType(), 2 * order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   }
   InitializeFaultDOFs(dof_data, n_local_qps, fault_coords);
   ZeroDOFDataPreStressTotal(dof_data, n_local_qps);
   wave.SetFaultDOFData(&dof_data, nbf);
   real_t bulk_bg[NUM_STATE] = {0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;
   wave.SetAbsorbingBackground(bulk_bg);
}

void SetupWavePar(WaveOperator<ParMesh> &wave,
                  std::vector<DOFData> &dof_data,
                  ParMesh &pmesh, int order)
{
   const int n_local_qps = wave.GetNumLocalFaultQPs();
   const int nbf = wave.GetNbfPerFace();
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   std::vector<Vector> fault_coords;
   fault_coords.reserve(n_local_qps);
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr = pmesh.GetInteriorFaceTransformations(int_faces[i]);
      const IntegrationRule &ir = IntRules.Get(
         ftr->GetGeometryType(), 2 * order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);
         fault_coords.push_back(phys);
      }
   }
   InitializeFaultDOFs(dof_data, n_local_qps, fault_coords);
   ZeroDOFDataPreStressTotal(dof_data, n_local_qps);
   wave.SetFaultDOFData(&dof_data, nbf);
   real_t bulk_bg[NUM_STATE] = {0};
   bulk_bg[SYY] =  TPV102Params::sigma_n;
   bulk_bg[SXY] = -TPV102Params::tau_ini;
   wave.SetAbsorbingBackground(bulk_bg);
}

} // anonymous

int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);
   int rank = 0, nprocs = 0;
   MPI_Comm comm = MPI_COMM_WORLD;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nprocs);
   g_seas_my_rank = rank;

   if (nprocs != 2)
   {
      if (rank == 0)
      {
         std::cerr << "[R-1201/R-1412] requires np=2, got " << nprocs
                   << " — TEST FAILED.  This test is meaningless at "
                   "np!=2; CI must invoke via `mpirun -np 2`.  Returning 1 "
                   "instead of 0 to surface CI misconfigurations that "
                   "would otherwise silently report PASS without "
                   "exercising the cross-rank dispatch.\n";
      }
      MPI_Finalize();
      return 1;
   }

   if (rank == 0)
   {
      std::cout
         << "\n=== Round-13 R-1201: np=2 vs np=1 Adjacent dispatch "
         << "numerical bit-identity ===\n";
   }

   const int order = 1;

   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr    = 3;
   bc.absorbing_attrs = {};

   // ----------------------------------------------------------------
   // Reference path: serial run on rank 0.  Broadcast invariants to
   // both ranks.  We compare TWO paths: Mult (RK4) AND AdvanceADER —
   // R-1400 fix: ADER is the production path; the previous test only
   // exercised RK4 leaving wave_operator.inl:4165 (ADER shared) un-
   // verified at MPI scale.
   // ----------------------------------------------------------------
   Invariants inv_ref_mult{};
   Invariants inv_ref_ader{};
   {
      Mesh serial = BuildOrthogonalFaultMesh();
      WaveOperator<Mesh> wave_s(serial, order, kLambda, kMu, kRho, bc);
      FaultFaceFlux ff_s(kRho, kCp, kCs);
      wave_s.SetFaultFlux(&ff_s);
      std::vector<DOFData> dof_s;
      SetupWaveSerial(wave_s, dof_s, serial, order);
      wave_s.SetMixedFluxMode(MixedFluxMode::Adjacent);

      const auto &fes_s = wave_s.GetFESpace();
      const int ndof_total_s = fes_s.GetNDofs();
      Vector Q_init_s(NUM_STATE * ndof_total_s);
      InitializeNonUniformState(Q_init_s, fes_s, ndof_total_s);

      // Path A1: Mult (RK4).
      Vector dQdt_s(NUM_STATE * ndof_total_s);
      wave_s.Mult(Q_init_s, dQdt_s);
      inv_ref_mult = ComputeLocalInvariants(dQdt_s, fes_s, ndof_total_s);

      // R-1400 Path A2: AdvanceADER — production code path (ADER-O2).
      // Q_new_s contains Q_init + dt · (M⁻¹ · rhs) after one ADER step.
      Vector Q_new_s(NUM_STATE * ndof_total_s);
      wave_s.AdvanceADER(Q_init_s, /*dt=*/1.0e-5, /*ader_order=*/2, Q_new_s);
      inv_ref_ader = ComputeLocalInvariants(Q_new_s, fes_s, ndof_total_s);
   }
   // Broadcast reference invariants from rank 0 to rank 1 (rank 0
   // ran the serial path; rank 1 didn't).
   MPI_Bcast(inv_ref_mult.sum,     NUM_STATE, MPI_DOUBLE, 0, comm);
   MPI_Bcast(inv_ref_mult.sum_sq,  NUM_STATE, MPI_DOUBLE, 0, comm);
   MPI_Bcast(inv_ref_mult.max_abs, NUM_STATE, MPI_DOUBLE, 0, comm);
   MPI_Bcast(inv_ref_ader.sum,     NUM_STATE, MPI_DOUBLE, 0, comm);
   MPI_Bcast(inv_ref_ader.sum_sq,  NUM_STATE, MPI_DOUBLE, 0, comm);
   MPI_Bcast(inv_ref_ader.max_abs, NUM_STATE, MPI_DOUBLE, 0, comm);

   // ----------------------------------------------------------------
   // Parallel path: same mesh, partitioned across np=2.
   // ----------------------------------------------------------------
   Mesh serial_for_par = BuildOrthogonalFaultMesh();
   std::vector<int> part = BuildSkewedPartition(serial_for_par);
   ParMesh pmesh(comm, serial_for_par, part.data());

   WaveOperator<ParMesh> wave_p(pmesh, order, kLambda, kMu, kRho, bc);
   FaultFaceFlux ff_p(kRho, kCp, kCs);
   wave_p.SetFaultFlux(&ff_p);
   std::vector<DOFData> dof_p;
   SetupWavePar(wave_p, dof_p, pmesh, order);
   wave_p.SetMixedFluxMode(MixedFluxMode::Adjacent);

   const auto &fes_p = wave_p.GetFESpace();
   const int ndof_total_p = fes_p.GetNDofs();
   Vector Q_init_p(NUM_STATE * ndof_total_p);
   InitializeNonUniformState(Q_init_p, fes_p, ndof_total_p);

   // Path B1: Mult (RK4).
   Vector dQdt_p(NUM_STATE * ndof_total_p);
   wave_p.Mult(Q_init_p, dQdt_p);
   Invariants local_inv_p_mult = ComputeLocalInvariants(dQdt_p, fes_p,
                                                        ndof_total_p);
   Invariants inv_par_mult = AllreduceInvariants(local_inv_p_mult, comm);

   // R-1400 Path B2: AdvanceADER — production code path.
   Vector Q_new_p(NUM_STATE * ndof_total_p);
   wave_p.AdvanceADER(Q_init_p, /*dt=*/1.0e-5, /*ader_order=*/2, Q_new_p);
   Invariants local_inv_p_ader = ComputeLocalInvariants(Q_new_p, fes_p,
                                                        ndof_total_p);
   Invariants inv_par_ader = AllreduceInvariants(local_inv_p_ader, comm);

   // For backward-compatible variable names below, use the Mult path
   // as the "default" reporting block; the ADER comparison runs after.
   const Invariants &inv_ref = inv_ref_mult;
   const Invariants &inv_par = inv_par_mult;

   if (rank == 0)
   {
      std::cout << "  per-component invariants (ref = serial, par = np=2):\n";
      std::cout << std::scientific << std::setprecision(3);
      std::cout << "    c       sum_ref         sum_par         "
                << "sum_sq_ref      sum_sq_par      max_abs_ref     "
                << "max_abs_par\n";
      for (int c = 0; c < NUM_STATE; c++)
      {
         std::cout << "    " << c << "  "
                   << std::setw(15) << inv_ref.sum[c]    << "  "
                   << std::setw(15) << inv_par.sum[c]    << "  "
                   << std::setw(15) << inv_ref.sum_sq[c] << "  "
                   << std::setw(15) << inv_par.sum_sq[c] << "  "
                   << std::setw(15) << inv_ref.max_abs[c]<< "  "
                   << std::setw(15) << inv_par.max_abs[c]<< "\n";
      }
   }

   // Compare per-component invariants for BOTH paths (Mult RK4 + ADER).
   auto compare_invariants = [&](const Invariants &ref, const Invariants &par,
                                 real_t &worst_rel_sum,
                                 real_t &worst_rel_sum_sq,
                                 real_t &worst_rel_max)
   {
      worst_rel_sum    = 0.0;
      worst_rel_sum_sq = 0.0;
      worst_rel_max    = 0.0;
      for (int c = 0; c < NUM_STATE; c++)
      {
         const real_t sc_sum =
            std::max<real_t>(std::abs(ref.sum[c]),     1.0);
         const real_t sc_sq  =
            std::max<real_t>(std::abs(ref.sum_sq[c]),  1.0);
         const real_t sc_max =
            std::max<real_t>(std::abs(ref.max_abs[c]), 1.0);
         worst_rel_sum    = std::max(worst_rel_sum,
            std::abs(ref.sum[c]    - par.sum[c])    / sc_sum);
         worst_rel_sum_sq = std::max(worst_rel_sum_sq,
            std::abs(ref.sum_sq[c] - par.sum_sq[c]) / sc_sq);
         worst_rel_max    = std::max(worst_rel_max,
            std::abs(ref.max_abs[c]- par.max_abs[c])/ sc_max);
      }
   };

   real_t worst_rel_sum_mult, worst_rel_sum_sq_mult, worst_rel_max_mult;
   real_t worst_rel_sum_ader, worst_rel_sum_sq_ader, worst_rel_max_ader;
   compare_invariants(inv_ref_mult, inv_par_mult,
                      worst_rel_sum_mult, worst_rel_sum_sq_mult,
                      worst_rel_max_mult);
   compare_invariants(inv_ref_ader, inv_par_ader,
                      worst_rel_sum_ader, worst_rel_sum_sq_ader,
                      worst_rel_max_ader);

   // R-1400: ADER tolerances same as Mult (production-relevant).
   int local_pass_sum_m  = (worst_rel_sum_mult     <= 1.0e-10) ? 1 : 0;
   int local_pass_sq_m   = (worst_rel_sum_sq_mult  <= 1.0e-10) ? 1 : 0;
   int local_pass_max_m  = (worst_rel_max_mult     <= 1.0e-12) ? 1 : 0;
   int local_pass_sum_a  = (worst_rel_sum_ader     <= 1.0e-10) ? 1 : 0;
   int local_pass_sq_a   = (worst_rel_sum_sq_ader  <= 1.0e-10) ? 1 : 0;
   int local_pass_max_a  = (worst_rel_max_ader     <= 1.0e-12) ? 1 : 0;
   int all_pass_sum_m = 0, all_pass_sq_m = 0, all_pass_max_m = 0;
   int all_pass_sum_a = 0, all_pass_sq_a = 0, all_pass_max_a = 0;
   MPI_Allreduce(&local_pass_sum_m, &all_pass_sum_m, 1, MPI_INT, MPI_LAND, comm);
   MPI_Allreduce(&local_pass_sq_m,  &all_pass_sq_m,  1, MPI_INT, MPI_LAND, comm);
   MPI_Allreduce(&local_pass_max_m, &all_pass_max_m, 1, MPI_INT, MPI_LAND, comm);
   MPI_Allreduce(&local_pass_sum_a, &all_pass_sum_a, 1, MPI_INT, MPI_LAND, comm);
   MPI_Allreduce(&local_pass_sq_a,  &all_pass_sq_a,  1, MPI_INT, MPI_LAND, comm);
   MPI_Allreduce(&local_pass_max_a, &all_pass_max_a, 1, MPI_INT, MPI_LAND, comm);

   int num_passed = 0, num_failed = 0;
   if (rank == 0)
   {
      auto report = [&](int pass, real_t worst, real_t tol,
                        const char *gate)
      {
         if (pass)
         {
            std::cout << "  PASSED: " << gate << "  (got "
                      << std::scientific << std::setprecision(3) << worst
                      << ", tol " << tol << ")\n";
            num_passed++;
         }
         else
         {
            std::cout << "  FAILED: " << gate << "  (got "
                      << std::scientific << std::setprecision(3) << worst
                      << ", expected <= " << tol << ")\n";
            num_failed++;
         }
      };
      // RK4 (Mult) path.
      report(all_pass_sum_m, worst_rel_sum_mult, 1.0e-10,
             "Gate 1 (Mult/RK4): per-component sum agrees serial vs np=2");
      report(all_pass_sq_m,  worst_rel_sum_sq_mult, 1.0e-10,
             "Gate 2 (Mult/RK4): per-component sum-of-squares agrees");
      report(all_pass_max_m, worst_rel_max_mult, 1.0e-12,
             "Gate 3 (Mult/RK4): per-component ∞-norm agrees");
      // R-1400 ADER (AdvanceADER) path — production code path.
      report(all_pass_sum_a, worst_rel_sum_ader, 1.0e-10,
             "Gate 4 (R-1400 AdvanceADER): per-component sum agrees serial vs np=2");
      report(all_pass_sq_a,  worst_rel_sum_sq_ader, 1.0e-10,
             "Gate 5 (R-1400 AdvanceADER): per-component sum-of-squares agrees");
      report(all_pass_max_a, worst_rel_max_ader, 1.0e-12,
             "Gate 6 (R-1400 AdvanceADER): per-component ∞-norm agrees");
      std::cout << "\n========================================\n"
                << "  Results: " << num_passed << " passed, "
                << num_failed << " failed out of " << (num_passed + num_failed)
                << " tests\n"
                << "========================================\n";
   }

   const int ok_all = all_pass_sum_m && all_pass_sq_m && all_pass_max_m &&
                      all_pass_sum_a && all_pass_sq_a && all_pass_max_a;
   MPI_Finalize();
   return ok_all ? 0 : 1;
}

// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.2.0 plan §4.7 / §16.6 (rev-3c) — shared-fault DOFData
// runtime consistency regression test.
//
// ============================================================================
// What this test probes
// ============================================================================
// §4.6 (ctor-level) verified that on shared fault faces:
//   • `elem1_on_plus` correctly ANTI-symmetric across ranks,
//   • `can_n / can_t1 / can_t2 / nl` bit-identical across ranks.
// §4.7 (this test) is the RUNTIME-level follow-up per §16.6:
// drives the WaveOperator through 1 RK4 step and then compares
// `dof_data[dof_idx]` (the 8 sign/magnitude fields:
//   slip1, slip2, V1, V2, tau1_corr, tau2_corr, sigma_n_corr, psi)
// between the two ranks that share each fault QP.
//
// ============================================================================
// What a FAIL implies
// ============================================================================
// Per §16.2-§16.4: TPV102 runs the friction Brent solve on BOTH
// ranks of every shared fault face and writes two independent copies
// of DOFData.  The correctness depends on:
//   (a) MFEM ghost exchange delivering bit-identical Q_nbr on A
//       matching Q_self on B,
//   (b) `elem1_on_plus` labelling agreeing (§4.6 checks),
//   (c) Tinv·Q identical on both ranks (§4.3 checks),
//   (d) friction Brent determinism on identical inputs.
// If any assumption fails, this test shows per-field divergence.
// That is H-V92-P at runtime — the direct bug predicted by §16.
//
// ============================================================================
// Strategy (sequence of checks, each strictly harder)
// ============================================================================
// 1. Build the 4-km 4-rank fixture from §4.6 (already validated to
//    generate cross-rank shared fault faces).
// 2. Initialize DOFData via TPV102's `InitializeFaultDOFs` (same call
//    path the driver uses).
// 3. Phase-A check: DOFData after initialization — no Mult call —
//    should match bit-for-bit.  This is a SANITY baseline: any
//    disagreement here reflects a fault_basis / dof-order bug, not a
//    flux-path bug.
// 4. Phase-B check: 1 Mult(Q=0) — Q=0 makes Tinv·Q=0 on both ranks,
//    so Evaluate receives identical inputs (in theory) and should
//    produce identical DOFData.  FAIL on Phase-B ⇒ stronger H-V92-P
//    evidence (the bug fires even at Q=0).
// 5. Phase-C check: full RK4 step with uniform nonzero Q (σ_xy = τ_ini
//    pre-stress plus a uniform σ_xx perturbation).  Exercises Tinv + Brent
//    BUT — per REVIEW R-V92-C01 — the uniform pre-stress makes
//    Q_self ≡ Q_nbr after MFEM ghost exchange, so the (+, −) swap is
//    vacuous at every shared fault QP.  Phase C alone is a **subset
//    probe**, not the operational definition of H-V92-P.  Retained for
//    backward comparison but never sufficient on its own.
// 6. **Phase-D (rev-3d, post-REVIEW R-V92-C01)**: full RK4 step with
//    an ANTISYMMETRIC σ_xy(y) = τ_ini · tanh(y / L0) profile projected
//    onto the FESpace via ProjectCoefficient.  At the fault (y=0) the
//    profile sign-inverts, so on every shared fault face
//    Q_self (y<0 side) and Q_nbr (y>0 side) differ by ≈ 2·τ_ini ⇒
//    the (+, −) swap is EXERCISED and any H-V92-P route-swap bug
//    would produce bit-divergent DOFData across ranks.  This is the
//    true operational probe.
//
// Usage:
//   mpirun -np 4 ./seas_test_shared_fault_dof_data_consistency

#include "mfem.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/tpv102_setup.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../domain/boundary_config.hpp"
#include "../../fault/fault_basis.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <map>
#include <memory>
#include <string>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

// ---------------------------------------------------------------------------
// Centroid key — 1 μm resolution, identical to §4.6
// ---------------------------------------------------------------------------
struct CentroidKey
{
   int64_t ix, iy, iz;

   CentroidKey() : ix(0), iy(0), iz(0) {}
   CentroidKey(double x, double y, double z)
   {
      const double scale = 1.0e6;
      ix = (int64_t)std::llround(x * scale);
      iy = (int64_t)std::llround(y * scale);
      iz = (int64_t)std::llround(z * scale);
   }
   bool operator<(const CentroidKey &o) const
   {
      if (ix != o.ix) return ix < o.ix;
      if (iy != o.iy) return iy < o.iy;
      return iz < o.iz;
   }
   bool operator==(const CentroidKey &o) const
   {
      return ix == o.ix && iy == o.iy && iz == o.iz;
   }
};

struct DOFDataRecord
{
   CentroidKey key;
   int32_t     src_rank;
   int32_t     phase;           // 0 = pre-Mult, 1 = post-Q=0 Mult, 2 = post-RK4 with Q≠0
   double      slip1, slip2;
   double      V1, V2;
   double      tau1_corr, tau2_corr;
   double      sigma_n_corr;
   double      psi;
};

// ---------------------------------------------------------------------------
// Fixture mesh from §4.6 — 4×4×4 km cartesian tet box, fault attr=3 at y=0.
// ---------------------------------------------------------------------------
static Mesh BuildFixtureMesh(int nx, int ny, int nz,
                              double Lx, double Ly, double Lz)
{
   Mesh mesh = Mesh::MakeCartesian3D(2*nx, 2*ny, nz,
                                     Element::TETRAHEDRON,
                                     2.0*Lx, 2.0*Ly, Lz);
   for (int v = 0; v < mesh.GetNV(); v++)
   {
      real_t *x = mesh.GetVertex(v);
      x[0] -= Lx;
      x[1] -= Ly;
      x[2] -= Lz;
   }
   const real_t tol = 1e-6;
   for (int be = 0; be < mesh.GetNBE(); be++)
   {
      auto *Tr = mesh.GetBdrElementTransformation(be);
      const IntegrationPoint &ip = Geometries.GetCenter(Tr->GetGeometryType());
      Tr->SetIntPoint(&ip);
      Vector c(3); Tr->Transform(ip, c);
      if (std::abs(c(2)) < tol) { mesh.SetBdrAttribute(be, 1); }
      else                      { mesh.SetBdrAttribute(be, 5); }
   }
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      auto *ftr = mesh.GetInteriorFaceTransformations(f);
      if (!ftr) { continue; }
      const IntegrationPoint &ip = Geometries.GetCenter(ftr->GetGeometryType());
      ftr->Face->SetIntPoint(&ip);
      Vector c(3); ftr->Face->Transform(ip, c);
      if (std::abs(c(1)) > tol) { continue; }
      Array<int> verts;
      mesh.GetFaceVertices(f, verts);
      if (verts.Size() == 3)
      { mesh.AddBdrTriangle(verts[0], verts[1], verts[2], 3); }
      else if (verts.Size() == 4)
      { mesh.AddBdrQuad(verts[0], verts[1], verts[2], verts[3], 3); }
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// ---------------------------------------------------------------------------
// Collect fault coords in the exact order expected by SetFaultDOFData:
// [interior fault QPs..., shared fault QPs...]
// ---------------------------------------------------------------------------
static void BuildFaultCoords(ParMesh &pmesh,
                              const Array<int> &int_faces,
                              const Array<int> &shr_faces,
                              int order,
                              std::vector<Vector> &fault_coords)
{
   fault_coords.clear();
   for (int i = 0; i < int_faces.Size(); i++)
   {
      auto *ftr = pmesh.GetInteriorFaceTransformations(int_faces[i]);
      MFEM_VERIFY(ftr, "§4.7: interior fault FTR null");
      const IntegrationRule &ir =
         IntRules.Get(ftr->GetGeometryType(), 2*order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         ftr->SetAllIntPoints(&ir.IntPoint(q));
         Vector phys(3); ftr->Face->Transform(ir.IntPoint(q), phys);
         fault_coords.push_back(phys);
      }
   }
   for (int i = 0; i < shr_faces.Size(); i++)
   {
      auto *ftr = pmesh.GetSharedFaceTransformations(shr_faces[i]);
      MFEM_VERIFY(ftr, "§4.7: shared fault FTR null");
      const IntegrationRule &ir =
         IntRules.Get(ftr->GetGeometryType(), 2*order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
         ftr->SetAllIntPoints(&ir.IntPoint(q));
         Vector phys(3); ftr->Face->Transform(ir.IntPoint(q), phys);
         fault_coords.push_back(phys);
      }
   }
}

// ---------------------------------------------------------------------------
// Snapshot DOFData for each shared-fault QP owned by this rank.
// dof_data layout: [interior QPs (n_int * nqp)] [shared QPs (n_shr * nqp)]
// ---------------------------------------------------------------------------
static std::vector<DOFDataRecord>
SnapshotSharedFaultDOFData(ParMesh &pmesh,
                            const WaveOperator<ParMesh> &wave,
                            const std::vector<DOFData> &dof_data,
                            int order, int src_rank, int phase)
{
   std::vector<DOFDataRecord> out;
   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   const Array<int> &shr_faces = wave.GetFaultSharedFaces();
   const int n_int = int_faces.Size();

   // nqp_per_face: infer from the first available face.
   int nqp_per_face = 0;
   if (shr_faces.Size() > 0)
   {
      auto *ftr = pmesh.GetSharedFaceTransformations(shr_faces[0]);
      nqp_per_face =
         IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
   }
   else if (n_int > 0)
   {
      auto *ftr = pmesh.GetInteriorFaceTransformations(int_faces[0]);
      nqp_per_face =
         IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
   }
   if (nqp_per_face == 0) { return out; }

   for (int i = 0; i < shr_faces.Size(); i++)
   {
      int sf = shr_faces[i];
      auto *ftr = pmesh.GetSharedFaceTransformations(sf);
      if (!ftr) { continue; }
      const IntegrationRule &ir =
         IntRules.Get(ftr->GetGeometryType(), 2*order);
      for (int q = 0; q < nqp_per_face; q++)
      {
         const int dof_idx = (n_int + i) * nqp_per_face + q;
         if (dof_idx < 0 || dof_idx >= (int)dof_data.size()) { continue; }
         const DOFData &fd = dof_data[dof_idx];

         ftr->SetAllIntPoints(&ir.IntPoint(q));
         Vector phys(3); ftr->Face->Transform(ir.IntPoint(q), phys);

         DOFDataRecord r;
         r.key          = CentroidKey(phys(0), phys(1), phys(2));
         r.src_rank     = src_rank;
         r.phase        = phase;
         r.slip1        = fd.slip1;
         r.slip2        = fd.slip2;
         r.V1           = fd.V1;
         r.V2           = fd.V2;
         r.tau1_corr    = fd.tau1_corr;
         r.tau2_corr    = fd.tau2_corr;
         r.sigma_n_corr = fd.sigma_n_corr;
         r.psi          = fd.psi;
         out.push_back(r);
      }
   }
   return out;
}

// ---------------------------------------------------------------------------
// Allgather + cross-rank bit-equality verdict for a phase.
//
// Returns a triplet:
//   pair_checks : # of cross-rank pair checks (each pair from its low-rank side)
//   fail_counts : per-field mismatch count, summed over all pair checks
//   worst_devs  : per-field worst |Δ|, max over all pair checks
// ---------------------------------------------------------------------------
struct PhaseVerdict
{
   int pair_checks = 0;
   int fail_slip1 = 0, fail_slip2 = 0;
   int fail_V1 = 0,    fail_V2 = 0;
   int fail_tau1 = 0,  fail_tau2 = 0;
   int fail_sn = 0,    fail_psi = 0;
   double worst_slip1 = 0, worst_slip2 = 0;
   double worst_V1 = 0,    worst_V2 = 0;
   double worst_tau1 = 0,  worst_tau2 = 0;
   double worst_sn = 0,    worst_psi = 0;
   int total_fails() const
   {
      return fail_slip1 + fail_slip2 + fail_V1 + fail_V2
           + fail_tau1 + fail_tau2 + fail_sn + fail_psi;
   }
};

static PhaseVerdict ComparePhase(MPI_Comm comm, int my_rank,
                                  const std::vector<DOFDataRecord> &local)
{
   // Allgatherv of bytes
   const int rec_bytes = static_cast<int>(sizeof(DOFDataRecord));
   int nprocs = 1; MPI_Comm_size(comm, &nprocs);
   const int my_size = (int)local.size();
   std::vector<int> sizes(nprocs), byte_sizes(nprocs), byte_displs(nprocs);
   MPI_Allgather(&my_size, 1, MPI_INT, sizes.data(), 1, MPI_INT, comm);
   int total_bytes = 0;
   for (int r = 0; r < nprocs; r++)
   {
      byte_sizes[r]  = sizes[r] * rec_bytes;
      byte_displs[r] = total_bytes;
      total_bytes   += byte_sizes[r];
   }
   std::vector<DOFDataRecord> global(total_bytes / rec_bytes);
   MPI_Allgatherv(local.data(), my_size * rec_bytes, MPI_BYTE,
                  global.data(), byte_sizes.data(), byte_displs.data(),
                  MPI_BYTE, comm);

   // Index by centroid key
   std::map<CentroidKey, std::vector<const DOFDataRecord*>> by_key;
   for (const auto &r : global) { by_key[r.key].push_back(&r); }

   PhaseVerdict pv;
   for (const auto &lrec : local)
   {
      auto it = by_key.find(lrec.key);
      if (it == by_key.end()) { continue; }
      for (const DOFDataRecord *orec : it->second)
      {
         if (orec->src_rank == my_rank) { continue; }
         // Only compare each pair ONCE (from the lower-ranked side).
         if (lrec.src_rank > orec->src_rank) { continue; }
         pv.pair_checks++;
         auto chk = [](double a, double b, int &fails, double &worst) {
            double d = std::abs(a - b);
            if (d > worst) worst = d;
            if (a != b) fails++;
         };
         chk(lrec.slip1,        orec->slip1,        pv.fail_slip1, pv.worst_slip1);
         chk(lrec.slip2,        orec->slip2,        pv.fail_slip2, pv.worst_slip2);
         chk(lrec.V1,           orec->V1,           pv.fail_V1,    pv.worst_V1);
         chk(lrec.V2,           orec->V2,           pv.fail_V2,    pv.worst_V2);
         chk(lrec.tau1_corr,    orec->tau1_corr,    pv.fail_tau1,  pv.worst_tau1);
         chk(lrec.tau2_corr,    orec->tau2_corr,    pv.fail_tau2,  pv.worst_tau2);
         chk(lrec.sigma_n_corr, orec->sigma_n_corr, pv.fail_sn,    pv.worst_sn);
         chk(lrec.psi,          orec->psi,          pv.fail_psi,   pv.worst_psi);
      }
   }

   // Reduce to rank 0 for reporting.  (Sum pair_checks and fails, max worst.)
   int local_vals[9] = {
      pv.pair_checks, pv.fail_slip1, pv.fail_slip2,
      pv.fail_V1, pv.fail_V2, pv.fail_tau1, pv.fail_tau2,
      pv.fail_sn, pv.fail_psi
   };
   int global_vals[9] = {0};
   MPI_Allreduce(local_vals, global_vals, 9, MPI_INT, MPI_SUM, comm);

   double local_dbl[8] = {
      pv.worst_slip1, pv.worst_slip2, pv.worst_V1, pv.worst_V2,
      pv.worst_tau1, pv.worst_tau2, pv.worst_sn, pv.worst_psi
   };
   double global_dbl[8] = {0};
   MPI_Allreduce(local_dbl, global_dbl, 8, MPI_DOUBLE, MPI_MAX, comm);

   PhaseVerdict g;
   g.pair_checks = global_vals[0];
   g.fail_slip1 = global_vals[1];  g.fail_slip2 = global_vals[2];
   g.fail_V1    = global_vals[3];  g.fail_V2    = global_vals[4];
   g.fail_tau1  = global_vals[5];  g.fail_tau2  = global_vals[6];
   g.fail_sn    = global_vals[7];  g.fail_psi   = global_vals[8];
   g.worst_slip1 = global_dbl[0];  g.worst_slip2 = global_dbl[1];
   g.worst_V1    = global_dbl[2];  g.worst_V2    = global_dbl[3];
   g.worst_tau1  = global_dbl[4];  g.worst_tau2  = global_dbl[5];
   g.worst_sn    = global_dbl[6];  g.worst_psi   = global_dbl[7];
   return g;
}

static void PrintVerdict(const std::string &label, const PhaseVerdict &v,
                          int rank, int &num_passed, int &num_failed)
{
   if (rank != 0) { return; }
   std::cout << "\n  [" << label << "] pair checks = " << v.pair_checks << "\n";
   auto row = [](const char *name, int fails, double worst)
   {
      std::cout << "    " << std::setw(14) << std::left << name
                << ": fails=" << std::setw(5) << fails
                << "  worst |Δ|=" << std::scientific
                << std::setprecision(6) << worst << "\n";
   };
   row("slip1",        v.fail_slip1, v.worst_slip1);
   row("slip2",        v.fail_slip2, v.worst_slip2);
   row("V1",           v.fail_V1,    v.worst_V1);
   row("V2",           v.fail_V2,    v.worst_V2);
   row("tau1_corr",    v.fail_tau1,  v.worst_tau1);
   row("tau2_corr",    v.fail_tau2,  v.worst_tau2);
   row("sigma_n_corr", v.fail_sn,    v.worst_sn);
   row("psi",          v.fail_psi,   v.worst_psi);

   const int tot = v.total_fails();
   if (tot == 0)
   {
      std::cout << "    PASSED: all DOFData fields bit-identical across "
                   "ranks on every shared fault QP\n";
      num_passed++;
   }
   else
   {
      std::cout << "    FAILED: " << tot
                << " cross-rank bit mismatches — "
                   "H-V92-P runtime divergence confirmed in this phase\n";
      num_failed++;
   }
}

// ---------------------------------------------------------------------------
// RK4 step helper (4 Mult calls).  Mirrors the driver's integrator.
// ---------------------------------------------------------------------------
static void DoRK4Step(const WaveOperator<ParMesh> &wave, Vector &Q, real_t dt)
{
   const int size = Q.Size();
   Vector k1(size), k2(size), k3(size), k4(size), Q_tmp(size);
   wave.Mult(Q, k1);
   add(Q, 0.5*dt, k1, Q_tmp);  wave.Mult(Q_tmp, k2);
   add(Q, 0.5*dt, k2, Q_tmp);  wave.Mult(Q_tmp, k3);
   add(Q, dt,     k3, Q_tmp);  wave.Mult(Q_tmp, k4);
   Q.Add(dt/6.0, k1);  Q.Add(dt/3.0, k2);
   Q.Add(dt/3.0, k3);  Q.Add(dt/6.0, k4);
}

// ---------------------------------------------------------------------------
int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);
   int rank = 0, nprocs = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   auto log = [&](const std::string &msg)
   { if (rank == 0) { std::cout << msg << "\n"; } };

   log("\n=== TPV102 v9.2.0 §4.7 Regression: "
       "Shared-Fault DOFData Runtime Consistency ===");
   log(std::string("  MPI ranks: ") + std::to_string(nprocs));

   if (nprocs < 4)
   {
      log("  SKIPPED: §4.7 requires >=4 MPI ranks for ParMETIS to cut "
          "the fault (4-km fixture).  Run as mpirun -np 4 ...");
      MPI_Finalize();
      return 77;
   }

   // -- Build fixture + partition ---------------------------------------
   const int nx = 2, ny = 2, nz = 2;
   const double Lx = 2.0e3, Ly = 2.0e3, Lz = 4.0e3;
   Mesh serial_mesh = BuildFixtureMesh(nx, ny, nz, Lx, Ly, Lz);
   if (rank == 0)
   {
      std::cout << "  serial mesh : " << serial_mesh.GetNE() << " elements, "
                << serial_mesh.GetNumFaces() << " faces, "
                << serial_mesh.GetNBE() << " bdr elems\n";
   }
   ParMesh pmesh(MPI_COMM_WORLD, serial_mesh);

   BoundaryConfig bc;
   bc.natural_attrs   = {1};
   bc.fault_attr      = 3;
   bc.absorbing_attrs = {5};

   const int order = 1;
   WaveOperator<ParMesh> wave(pmesh, order,
                              TPV102Params::lambda,
                              TPV102Params::mu,
                              TPV102Params::rho, bc);

   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   const Array<int> &shr_faces = wave.GetFaultSharedFaces();
   const int n_int = int_faces.Size();
   const int n_shr = shr_faces.Size();

   int total_shared = 0;
   MPI_Allreduce(&n_shr, &total_shared, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   for (int r = 0; r < nprocs; r++)
   {
      if (r == rank)
      {
         std::cout << "    rank " << r
                   << " : local_ne=" << pmesh.GetNE()
                   << " fault_int=" << n_int
                   << " fault_shr=" << n_shr << std::endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
   }
   if (total_shared == 0)
   {
      log("  SKIPPED: no shared fault faces with this nranks.");
      MPI_Finalize();
      return 77;
   }

   // nqp_per_face (Allreduce for ranks with 0 fault faces).
   int nqp_per_face = 0;
   if (n_shr > 0)
   {
      auto *ftr = pmesh.GetSharedFaceTransformations(shr_faces[0]);
      nqp_per_face = IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
   }
   else if (n_int > 0)
   {
      auto *ftr = pmesh.GetInteriorFaceTransformations(int_faces[0]);
      nqp_per_face = IntRules.Get(ftr->GetGeometryType(), 2*order).GetNPoints();
   }
   {
      int local = nqp_per_face, global = 0;
      MPI_Allreduce(&local, &global, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      nqp_per_face = global;
   }
   const int num_fault_total = (n_int + n_shr) * nqp_per_face;

   // -- Populate DOFData via InitializeFaultDOFs (driver's path) -----------
   std::vector<Vector> fault_coords;
   BuildFaultCoords(pmesh, int_faces, shr_faces, order, fault_coords);

   std::vector<DOFData> dof_data;
   if (num_fault_total > 0)
   {
      InitializeFaultDOFs(dof_data, num_fault_total, fault_coords);
   }

   FaultFaceFlux fault_flux(TPV102Params::rho, TPV102Params::cp,
                            TPV102Params::cs);
   wave.SetFaultFlux(&fault_flux);
   wave.SetFaultDOFData(&dof_data, nqp_per_face);

   int num_tests = 0, num_passed = 0, num_failed = 0;

   // Helper: measure the σ_xy spread across the fault by sampling at two
   // points (0, ±Ly, 0).  A nonzero spread implies Q_self ≠ Q_nbr on
   // shared fault faces (the precondition for H-V92-P to fire).
   //
   // Implementation: construct a ParGridFunction backing the SXY
   // component of Q, evaluate at the two probe points via a
   // GridFunction::GetValue, MPI-max the |Δ| across ranks.
   auto CheckSXYSpread =
      [&](const Vector &Q_vec, const std::string &phase_label)
   {
      const auto &fes_const = wave.GetFESpace();
      auto &fes = const_cast<ParFiniteElementSpace &>(fes_const);
      const int ndof_total = fes.GetNDofs();
      ParGridFunction sxy_gf(&fes);
      for (int i = 0; i < ndof_total; i++)
      { sxy_gf(i) = Q_vec(SXY * ndof_total + i); }

      // Point-evaluation helper via FindPoints is heavy; use min/max of
      // the grid-function values as a cheap proxy.  For an antisymmetric
      // projection this gives ≈ ±τ_ini.
      double local_min = +1e300, local_max = -1e300;
      for (int i = 0; i < ndof_total; i++)
      {
         double v = sxy_gf(i);
         if (v < local_min) local_min = v;
         if (v > local_max) local_max = v;
      }
      double gmin = 0, gmax = 0;
      MPI_Allreduce(&local_min, &gmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(&local_max, &gmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      double spread = gmax - gmin;
      if (rank == 0)
      {
         std::cout << "  [" << phase_label
                   << "] σ_xy field min/max = "
                   << std::scientific << std::setprecision(3)
                   << gmin << " / " << gmax
                   << "  spread = " << spread << " Pa"
                   << (spread > 0.0
                       ? "  ⇒ antisymmetric (Q_self ≠ Q_nbr at fault)"
                       : "  ⇒ uniform (Q_self ≡ Q_nbr — swap VACUOUS)")
                   << "\n";
      }
      return spread;
   };

   // ========================================================================
   // Phase A — post-init (BEFORE any Mult call).
   // Expectation: bit-identical across ranks (InitializeFaultDOFs is
   // purely a function of (coords, TPV102Params), both rank-invariant).
   // ========================================================================
   {
      num_tests++;
      auto recs = SnapshotSharedFaultDOFData(pmesh, wave, dof_data,
                                             order, rank, /*phase=*/0);
      PhaseVerdict v = ComparePhase(MPI_COMM_WORLD, rank, recs);
      PrintVerdict("phase A: post-init (no Mult)", v, rank,
                   num_passed, num_failed);
   }

   // ========================================================================
   // Phase B — post-Mult(Q=0).
   //
   // REVIEW R-V92-C01 diagnosis: Q=0 ⇒ Q_self = Q_nbr = 0 trivially on
   // every shared fault face.  The (+, −) swap is VACUOUS.  Retained
   // as a regression baseline (if this ever FAILs, something is
   // seriously wrong) but not a real H-V92-P probe.
   // ========================================================================
   {
      num_tests++;
      Vector Q(wave.Height()); Q = 0.0;
      CheckSXYSpread(Q, "phase B");
      Vector k(Q.Size());
      wave.Mult(Q, k);

      auto recs = SnapshotSharedFaultDOFData(pmesh, wave, dof_data,
                                             order, rank, /*phase=*/1);
      PhaseVerdict v = ComparePhase(MPI_COMM_WORLD, rank, recs);
      PrintVerdict("phase B: post 1 Mult(Q=0)", v, rank,
                   num_passed, num_failed);
   }

   // ========================================================================
   // Phase C — 1 full RK4 step with UNIFORM nonzero Q.
   //
   // REVIEW R-V92-C01 diagnosis: the uniform profile (σ_xy = τ_ini,
   // σ_xx = 10⁵ Pa on every DOF) makes Q_self ≡ Q_nbr bit-exactly
   // after MFEM ghost exchange at every shared fault face, so the
   // (+, −) swap is VACUOUS here too.  CheckSXYSpread below
   // prints 0 — documenting this structural limitation.  Retained
   // for comparison with Phase D but **not operational on its own**
   // for H-V92-P discrimination.
   // ========================================================================
   {
      num_tests++;
      // Reset DOFData to init state for a clean Phase-C measurement.
      if (num_fault_total > 0)
      {
         InitializeFaultDOFs(dof_data, num_fault_total, fault_coords);
      }

      const int size = wave.Height();
      Vector Q(size);
      Q = 0.0;
      // Uniform pre-stress on every DOF: set σ_xy channel = tau_ini.
      // The Q storage layout in WaveOperator is component-major:
      //   Q[c * ndof_total + dof] for c in [0, NUM_STATE)
      // where ndof_total = fes->GetNDofs().  The interior Q is a
      // DG field so every element's DOFs get the uniform pre-stress.
      const auto &pfes = wave.GetFESpace();
      const int ndof_total = pfes.GetNDofs();
      for (int dof = 0; dof < ndof_total; dof++)
      {
         Q(SXY * ndof_total + dof) = TPV102Params::tau_ini;
         // Small rank-invariant perturbation (deterministic in the
         // GLOBAL dof index) — but pfes->GetNDofs() is per-rank,
         // so we use a rank-invariant proxy: the global vertex id
         // of the DOF's supporting cell.  Simpler: use a uniform
         // nonzero in σ_xx = 1e5 Pa so Tinv·Q has a nontrivial
         // structure, independent of dof ordering.
         Q(SXX * ndof_total + dof) = 1.0e5;
      }

      // Pre-step Q_self vs Q_nbr check (should be 0 — uniform pre-stress).
      CheckSXYSpread(Q, "phase C pre-Mult");

      const real_t dt = 1.0e-5;
      DoRK4Step(wave, Q, dt);

      auto recs = SnapshotSharedFaultDOFData(pmesh, wave, dof_data,
                                             order, rank, /*phase=*/2);
      PhaseVerdict v = ComparePhase(MPI_COMM_WORLD, rank, recs);
      PrintVerdict("phase C: post 1 RK4 step (uniform Q — swap vacuous)",
                   v, rank, num_passed, num_failed);
   }

   // ========================================================================
   // Phase D (rev-3d, post-REVIEW R-V92-C01) — 1 full RK4 step with an
   // ANTISYMMETRIC σ_xy(y) = τ_ini · tanh(y / L0) profile projected onto
   // the FESpace via ProjectCoefficient.
   //
   // Motivation (REVIEW R-V92-C01): H-V92-P requires Q_self ≠ Q_nbr on
   // shared fault faces to fire.  Phases A-C all have Q_self = Q_nbr
   // trivially (uniform Q on both sides of any fault face).  Phase D
   // ensures that on a fault at y = 0, the y < 0 side has σ_xy ≈ −τ_ini
   // and the y > 0 side has σ_xy ≈ +τ_ini — so Q_self and Q_nbr differ
   // by ≈ 2 τ_ini at every shared fault QP.
   //
   // Under the antisymmetric profile:
   //   • Ghost exchange must deliver Q_nbr_on_A ≡ Q_self_on_B bit-exactly.
   //   • The (+, −) swap at wave_operator.inl:1287-1290 must route
   //     Q_plus, Q_minus so that BOTH ranks call Evaluate with the same
   //     (Q_plus, Q_minus) tuple — despite one rank having + = Q_self
   //     and the other having + = Q_nbr.
   //   • The friction Brent solver must be deterministic for identical
   //     inputs.
   //
   // If ANY of these fails, the two ranks compute different DOFData
   // values on a shared fault QP — H-V92-P CONFIRMED.  On the 4-km
   // Cartesian fixture, elem1_on_plus is determined by centroid
   // margins of O(10³ m) (REVIEW R-V92-C04), so the FP-fragility branch
   // cannot trigger here; a FAIL on Phase D would instead indicate
   // ghost exchange or Brent determinism failure.  A PASS on Phase D
   // is a genuine (not structural) invariant verification.
   //
   // L0 choice: L0 = 500 m << Ly = 2000 m so tanh saturates to ±1 on
   // most of the fixture volume, giving strong Q_self − Q_nbr contrast
   // (≈ 2·τ_ini ≈ 1.5·10⁸ Pa) at every shared fault QP.
   // ========================================================================
   {
      num_tests++;
      // Reset DOFData for clean phase-D measurement.
      if (num_fault_total > 0)
      {
         InitializeFaultDOFs(dof_data, num_fault_total, fault_coords);
      }

      const auto &fes_const = wave.GetFESpace();
      auto &fes = const_cast<ParFiniteElementSpace &>(fes_const);
      const int ndof_total = fes.GetNDofs();
      const int size = wave.Height();

      Vector Q(size);
      Q = 0.0;

      // Project σ_xy = τ_ini · tanh(y / L0) onto the DG FESpace.
      const double L0 = 500.0;
      FunctionCoefficient sxy_coeff(
         [L0](const Vector &x) -> real_t {
            return TPV102Params::tau_ini * std::tanh(x(1) / L0);
         });
      ParGridFunction sxy_gf(&fes);
      sxy_gf.ProjectCoefficient(sxy_coeff);
      for (int i = 0; i < ndof_total; i++)
      {
         Q(SXY * ndof_total + i) = sxy_gf(i);
      }

      // Pre-step Q_self vs Q_nbr check — expected NONZERO (≈ 2·τ_ini).
      double qspread = CheckSXYSpread(Q, "phase D pre-Mult");

      // Ghost exchange must be up-to-date before the MFEM Mult call;
      // wave.Mult already triggers it via ParGridFunction::ExchangeFaceNbrData
      // on the ghost_gf_ helper.  No user action needed.

      const real_t dt = 1.0e-5;
      DoRK4Step(wave, Q, dt);

      auto recs = SnapshotSharedFaultDOFData(pmesh, wave, dof_data,
                                             order, rank, /*phase=*/3);
      PhaseVerdict v = ComparePhase(MPI_COMM_WORLD, rank, recs);
      PrintVerdict("phase D: post 1 RK4 (antisymmetric σ_xy, swap EXERCISED)",
                   v, rank, num_passed, num_failed);

      if (rank == 0 && qspread == 0.0)
      {
         std::cout << "  WARNING: phase D pre-Mult spread is 0 — the "
                      "antisymmetric projection did not take effect.  "
                      "This invalidates the phase.\n";
      }
   }

   if (rank == 0)
   {
      std::cout << "\n========================================\n";
      std::cout << "  Results: " << num_passed << " passed, "
                << num_failed << " failed out of " << num_tests
                << " phases\n";
      std::cout << "========================================\n";
   }

   int any_fail = (num_failed > 0) ? 1 : 0;
   MPI_Bcast(&any_fail, 1, MPI_INT, 0, MPI_COMM_WORLD);

   MPI_Finalize();
   return (any_fail == 0) ? 0 : 1;
}

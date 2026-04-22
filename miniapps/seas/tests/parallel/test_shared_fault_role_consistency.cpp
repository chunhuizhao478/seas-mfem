// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// TPV102 v9.2.0 plan §4.6 / §15.5.1 (rev-3b) — shared-fault role
// consistency regression test.
//
// ============================================================================
// What this test probes
// ============================================================================
// H-V92-P (rev-3b): the two ranks that share a fault QP must produce
// a bit-identical (+, −) routing for the Pelties-9 flux.  Contract:
//   (i)   elem1_on_plus MUST be ANTI-SYMMETRIC across ranks:
//         exactly ONE rank has true, the OTHER has false.  This is
//         the designed behaviour — each rank's local "Elem1" is a
//         different physical element (the local one), and one is on
//         the canonical + side, the other on the − side.  If BOTH
//         ranks have the SAME bit, the (+, −) routing at
//         wave_operator.inl:1287-1290 disagrees on the same face →
//         different Q_plus/Q_minus on the two ranks → different
//         σ_n_trial → H-V92-P CONFIRMED.
//   (ii)  can_n[3]  must match BIT-for-BIT across ranks (canonical
//         frame is ref_normal-aligned and rank-invariant per R-802).
//   (iii) can_t1[3] must match BIT-for-BIT across ranks.
//   (iv)  can_t2[3] must match BIT-for-BIT across ranks.
//   (v)   nl (|n_raw|) must match BIT-for-BIT across ranks.
//
// ============================================================================
// How it works
// ============================================================================
// 1. Build a tiny inline tet fixture mesh (4×4×4 km box, 4×4×2 tets
//    in each partitioned half, fault at y=0 with attribute 3 on every
//    interior y=0 face).  The fixture is small (few hundred elements)
//    so ParMETIS has freedom to partition across the fault.
// 2. Partition to 2 MPI ranks via ParMesh.
// 3. Instantiate WaveOperator<ParMesh> (same ctor path that fires in
//    production).  This populates fault_basis_ (per-QP) for both
//    interior and shared fault faces.
// 4. For each shared fault face, each rank writes a record:
//      key   = face centroid (x, y, z) rounded to 1e-6 m
//      value = (elem1_on_plus | sign_flipped, can_n, can_t1, can_t2)
// 5. MPI_Alltoallv swap of records; each rank matches its own
//    records against the partner's by centroid-key and asserts
//    bit-for-bit equality on EVERY field.
// 6. FAIL on any mismatch (reports the disagreeing fields + per-rank
//    values).  PASS iff every shared face matches across ranks.
//
// ============================================================================
// Expected outcome
// ============================================================================
// If v9.0.0 R-802 (canonical normal via FaultBasis) holds for EVERY
// sign-determining quantity: PASS.  The test then functions as a
// regression gate.
//
// If H-V92-P is the root cause of R-V92: FAIL, with the failing field
// named (elem1_on_plus bit, or ULP-level can_n drift, etc.).
//
// ============================================================================
// Usage
// ============================================================================
//   mpirun -np 2 ./seas_test_shared_fault_role_consistency
//   mpirun -np 4 ./seas_test_shared_fault_role_consistency  (more
//                  partition patterns; not required for the gate)

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
#include <iostream>
#include <map>
#include <memory>
#include <vector>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;
using namespace mfem::seas;

// ---------------------------------------------------------------------------
// Round centroid key to 1e-6 m to paper over ULP-level partition-side
// differences in the Jacobian's computed centroid, while keeping the key
// unique across the (~hundred) shared faces in the fixture.
// ---------------------------------------------------------------------------
struct CentroidKey
{
   int64_t ix, iy, iz;

   CentroidKey() : ix(0), iy(0), iz(0) {}
   CentroidKey(double x, double y, double z)
   {
      const double scale = 1.0e6;  // 1 μm resolution
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

// Per-QP record exchanged across ranks.
struct QPRecord
{
   CentroidKey key;
   int32_t     src_rank;          // originating rank (for self-vs-partner skip)
   int32_t     elem1_on_plus;     // 0 or 1
   int32_t     sign_flipped;      // 0 or 1
   double      n[3];
   double      t1[3];
   double      t2[3];
   double      nl;                // |n_raw| (|Jacobian|)
};

// ---------------------------------------------------------------------------
// Fixture mesh: 4×4×4 km box, Nx × Ny × Nz cells split into tets, fault
// at y=0 gets attribute 3 on every interior face.
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
      x[2] -= Lz;           // z in [-Lz, 0]
   }

   // Boundary attributes: 1 = top (z=0), 5 = sides+bottom, 3 = fault (y=0)
   const real_t tol = 1e-6;
   for (int be = 0; be < mesh.GetNBE(); be++)
   {
      auto *Tr = mesh.GetBdrElementTransformation(be);
      const IntegrationPoint &ip = Geometries.GetCenter(Tr->GetGeometryType());
      Tr->SetIntPoint(&ip);
      Vector c(3); Tr->Transform(ip, c);
      if (std::abs(c(2)) < tol)
      { mesh.SetBdrAttribute(be, 1); }
      else
      { mesh.SetBdrAttribute(be, 5); }
   }

   // Add fault boundary elements on every interior y=0 face.
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
      {
         mesh.AddBdrTriangle(verts[0], verts[1], verts[2], 3);
      }
      else if (verts.Size() == 4)
      {
         mesh.AddBdrQuad(verts[0], verts[1], verts[2], verts[3], 3);
      }
   }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return mesh;
}

// ---------------------------------------------------------------------------
// Extract (elem1_on_plus, can_n, can_t1, can_t2) for each shared-fault QP
// on this rank.  Uses the same logic as wave_operator.inl:
//   elem1_on_plus = !sign_flipped   (canonical frame aligned with
//                                    FaultBasis ref_normal)
//   can_{n,t1,t2} = ±(stored qp_data) depending on sign_flipped
// ---------------------------------------------------------------------------
static std::vector<QPRecord>
CollectSharedFaultQPRecords(ParMesh &pmesh,
                             const WaveOperator<ParMesh> &wave,
                             int order, int src_rank)
{
   std::vector<QPRecord> out;
   const FaultBasis *fb = wave.GetFaultBasis();
   if (fb == nullptr) { return out; }

   const Array<int> &int_faces = wave.GetFaultInteriorFaces();
   const Array<int> &shr_faces = wave.GetFaultSharedFaces();
   const int n_int = int_faces.Size();

   for (int i = 0; i < shr_faces.Size(); i++)
   {
      int sf = shr_faces[i];
      auto *ftr = pmesh.GetSharedFaceTransformations(sf);
      if (!ftr) { continue; }

      const int fb_idx = n_int + i;
      MFEM_VERIFY(fb_idx < fb->NumFaces(),
                  "§4.6: shared-face basis index out of range");
      const FaultBasisData &bd = fb->GetBasis(fb_idx);

      const IntegrationRule &ir =
         IntRules.Get(ftr->GetGeometryType(), 2*order);
      const int nq = ir.GetNPoints();
      MFEM_VERIFY((int)bd.qp_data.size() == nq,
                  "§4.6: qp_data size mismatch on shared face " << sf);

      for (int q = 0; q < nq; q++)
      {
         const IntegrationPoint &ip = ir.IntPoint(q);
         ftr->SetAllIntPoints(&ip);
         Vector phys(3); ftr->Face->Transform(ip, phys);

         const FaultBasisQPData &qpd = bd.qp_data[q];

         QPRecord rec;
         rec.key           = CentroidKey(phys(0), phys(1), phys(2));
         rec.src_rank      = src_rank;
         rec.sign_flipped  = qpd.sign_flipped ? 1 : 0;
         rec.elem1_on_plus = qpd.sign_flipped ? 0 : 1;
         rec.nl            = qpd.nl;
         for (int d = 0; d < 3; d++)
         {
            // Canonical (ref-aligned) frame = stored frame "undone" by
            // the sign_flipped flip that FaultBasis::ComputeOrientedFrame
            // applied in Step 5.  On the two ranks sharing the face,
            // the canonical frame is rank-INVARIANT (by construction,
            // v9.0.0 R-802).  The stored basis may differ by a global
            // sign between ranks; the canonical must not.
            rec.n[d]  = qpd.sign_flipped ? -qpd.normal[d]
                                         :  qpd.normal[d];
            rec.t1[d] = qpd.sign_flipped ? -qpd.tangent1[d]
                                         :  qpd.tangent1[d];
            rec.t2[d] = qpd.sign_flipped ? -qpd.tangent2[d]
                                         :  qpd.tangent2[d];
         }
         out.push_back(rec);
      }
   }
   return out;
}

// ---------------------------------------------------------------------------
// MPI exchange: pack my QP records (with key), send to the "partner" rank.
//
// For simplicity in a 2-rank fixture: pack+send to ALL other ranks via
// Allgatherv, then each rank matches its own records against the union
// and verifies pairwise equality.  O(N²) in face count but N is small
// (100s of shared faces in the fixture).
// ---------------------------------------------------------------------------
static void AllgatherRecords(MPI_Comm comm,
                              const std::vector<QPRecord> &local,
                              std::vector<QPRecord> &global)
{
   int nprocs = 1;
   MPI_Comm_size(comm, &nprocs);
   const int my_size = static_cast<int>(local.size());

   std::vector<int> sizes(nprocs), displs(nprocs);
   MPI_Allgather(&my_size, 1, MPI_INT,
                 sizes.data(), 1, MPI_INT, comm);

   const int rec_bytes = static_cast<int>(sizeof(QPRecord));
   std::vector<int> byte_sizes(nprocs), byte_displs(nprocs);
   int total_bytes = 0;
   for (int r = 0; r < nprocs; r++)
   {
      byte_sizes[r] = sizes[r] * rec_bytes;
      byte_displs[r] = total_bytes;
      total_bytes   += byte_sizes[r];
   }

   global.resize(total_bytes / rec_bytes);
   MPI_Allgatherv(local.data(), my_size * rec_bytes, MPI_BYTE,
                  global.data(), byte_sizes.data(), byte_displs.data(),
                  MPI_BYTE, comm);
}

// ---------------------------------------------------------------------------
// Main test body
// ---------------------------------------------------------------------------
int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);

   int rank = 0, nprocs = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   int num_tests = 0, num_passed = 0, num_failed = 0;

   auto log = [&](const std::string &msg)
   {
      if (rank == 0) { std::cout << msg << "\n"; }
   };

   log("\n=== TPV102 v9.2.0 §4.6 Regression: "
       "Shared-Fault Role Consistency ===");
   log(std::string("  MPI ranks: ") + std::to_string(nprocs));

   if (nprocs < 2)
   {
      log("  SKIPPED: this test requires >=2 MPI ranks "
          "(run as mpirun -np 2 ...).");
      MPI_Finalize();
      return 77;  // skipped-sentinel (automake convention)
   }

   // -------- Build fixture mesh + partition to ParMesh --------------------
   // 4×4×4 km box.  Tet resolution: (nx, ny, nz) = (2, 2, 2) gives
   // 4×4×2 tets along each axis (since MakeCartesian3D takes 2*nx, 2*ny),
   // total = 2·4·4·2 tets per fundamental cell × 6 = hundreds of tets.
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

   // -------- Build WaveOperator (exercises the real ctor) -----------------
   BoundaryConfig bc;
   bc.natural_attrs = {1};
   bc.fault_attr    = 3;
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

   // Print per-rank layout (serialised for readability).
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

   int total_shared = 0;
   MPI_Allreduce(&n_shr, &total_shared, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if (total_shared == 0)
   {
      log("  SKIPPED: no shared fault faces with this nranks.  "
          "Try increasing nprocs or enlarging the fixture.");
      MPI_Finalize();
      return 77;
   }
   if (rank == 0)
   {
      std::cout << "  total shared fault face entries (summed over ranks): "
                << total_shared << "\n";
   }

   // -------- Collect per-QP records + Allgatherv across ranks ------------
   std::vector<QPRecord> local =
      CollectSharedFaultQPRecords(pmesh, wave, order, rank);
   if (rank == 0)
   {
      std::cout << "  my shared-fault QP records: " << local.size() << "\n";
   }

   std::vector<QPRecord> global;
   AllgatherRecords(MPI_COMM_WORLD, local, global);

   // -------- Bit-for-bit comparison -------------------------------------
   // For every record I OWN (key in my local set), find the partner(s) in
   // the global union (any rank with the same centroid key); assert bit-
   // equality on (elem1_on_plus, n, t1, t2, nl).  sign_flipped bit is
   // NOT required to match — FaultBasis may produce opposite stored
   // signs on the two ranks (§v61 analysis), but the canonical frame
   // (encoded here as n/t1/t2 after undoing sign_flipped) MUST match.
   std::map<CentroidKey, std::vector<const QPRecord*>> by_key;
   for (const auto &r : global) { by_key[r.key].push_back(&r); }

   int bad_elem1_on_plus = 0;  // pairs that violate ANTI-SYMMETRY
   int mismatches_n      = 0;
   int mismatches_t1     = 0;
   int mismatches_t2     = 0;
   int mismatches_nl     = 0;
   int pair_checks       = 0;

   // Per-field worst deviation for the diagnostic print.
   double worst_n_dev = 0.0, worst_t1_dev = 0.0, worst_t2_dev = 0.0;
   double worst_nl_dev = 0.0;

   for (const auto &lrec : local)
   {
      auto it = by_key.find(lrec.key);
      if (it == by_key.end()) { continue; }  // uniquely mine
      for (const QPRecord *orec : it->second)
      {
         // Skip self-compare: orec from my own rank is a copy of lrec
         // (Allgatherv packs my records into `global` too).  Compare
         // only across ranks.
         if (orec->src_rank == rank) { continue; }
         pair_checks++;

         // Contract: on a shared fault face, the two ranks' local
         // "Elem1" instances sit on OPPOSITE physical sides, so
         // elem1_on_plus MUST invert.  Equal bits across ranks mean
         // both ranks would route Q_plus_local = Q_self_can — the
         // H-V92-P bug signature.  Count "bad" as equal-bit pairs.
         if (lrec.elem1_on_plus == orec->elem1_on_plus)
         { bad_elem1_on_plus++; }

         for (int d = 0; d < 3; d++)
         {
            const double dn  = std::abs(lrec.n[d]  - orec->n[d]);
            const double dt1 = std::abs(lrec.t1[d] - orec->t1[d]);
            const double dt2 = std::abs(lrec.t2[d] - orec->t2[d]);
            if (dn  > worst_n_dev)  { worst_n_dev  = dn; }
            if (dt1 > worst_t1_dev) { worst_t1_dev = dt1; }
            if (dt2 > worst_t2_dev) { worst_t2_dev = dt2; }
            if (lrec.n [d] != orec->n [d]) { mismatches_n++;  }
            if (lrec.t1[d] != orec->t1[d]) { mismatches_t1++; }
            if (lrec.t2[d] != orec->t2[d]) { mismatches_t2++; }
         }

         const double dnl = std::abs(lrec.nl - orec->nl);
         if (dnl > worst_nl_dev) { worst_nl_dev = dnl; }
         if (lrec.nl != orec->nl) { mismatches_nl++; }
      }
   }

   // Reduce per-rank mismatch counts to rank 0 for the verdict.
   int total_bad_elem1 = 0, total_mis_n = 0, total_mis_t1 = 0;
   int total_mis_t2    = 0, total_mis_nl = 0, total_pair_checks = 0;
   double g_worst_n = 0, g_worst_t1 = 0, g_worst_t2 = 0, g_worst_nl = 0;

   MPI_Allreduce(&bad_elem1_on_plus,        &total_bad_elem1, 1,
                 MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&mismatches_n,             &total_mis_n,     1,
                 MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&mismatches_t1,            &total_mis_t1,    1,
                 MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&mismatches_t2,            &total_mis_t2,    1,
                 MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&mismatches_nl,            &total_mis_nl,    1,
                 MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&pair_checks,              &total_pair_checks, 1,
                 MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&worst_n_dev,  &g_worst_n,  1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   MPI_Allreduce(&worst_t1_dev, &g_worst_t1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   MPI_Allreduce(&worst_t2_dev, &g_worst_t2, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   MPI_Allreduce(&worst_nl_dev, &g_worst_nl, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

   // -------- Report --------------------------------------------------------
   num_tests += 5;
   if (rank == 0)
   {
      std::cout << "\n  Shared-fault role audit (summed over ranks):\n";
      std::cout << "    pair-wise comparisons  : " << total_pair_checks << "\n";
      std::cout << "    elem1_on_plus EQUAL    : " << total_bad_elem1
                << " (expect 0 — should be ANTI-symmetric; non-zero ⇒ H-V92-P)\n";
      std::cout << "    can_n mismatches       : " << total_mis_n
                << " entries; worst |Δ| = " << g_worst_n << "\n";
      std::cout << "    can_t1 mismatches      : " << total_mis_t1
                << " entries; worst |Δ| = " << g_worst_t1 << "\n";
      std::cout << "    can_t2 mismatches      : " << total_mis_t2
                << " entries; worst |Δ| = " << g_worst_t2 << "\n";
      std::cout << "    nl (|n_raw|) mismatches: " << total_mis_nl
                << " entries; worst |Δ| = " << g_worst_nl << "\n";

      auto verdict = [&](const char *name, int bad)
      {
         std::cout << "    " << (bad == 0 ? "PASSED: " : "FAILED: ")
                   << name << "\n";
         return (bad == 0);
      };
      if (verdict("elem1_on_plus is ANTI-symmetric across ranks (0 equal-bit pairs)",
                  total_bad_elem1)) { num_passed++; } else { num_failed++; }
      if (verdict("can_n  bit-identical across ranks",
                  total_mis_n))    { num_passed++; } else { num_failed++; }
      if (verdict("can_t1 bit-identical across ranks",
                  total_mis_t1))   { num_passed++; } else { num_failed++; }
      if (verdict("can_t2 bit-identical across ranks",
                  total_mis_t2))   { num_passed++; } else { num_failed++; }
      if (verdict("|n_raw| bit-identical across ranks",
                  total_mis_nl))   { num_passed++; } else { num_failed++; }

      std::cout << "\n========================================\n";
      std::cout << "  Results: " << num_passed << " passed, "
                << num_failed << " failed out of " << num_tests << " tests\n";
      std::cout << "========================================\n";
   }

   int any_fail = (total_bad_elem1 + total_mis_n + total_mis_t1
                   + total_mis_t2 + total_mis_nl > 0) ? 1 : 0;
   MPI_Bcast(&any_fail, 1, MPI_INT, 0, MPI_COMM_WORLD);

   MPI_Finalize();
   return (any_fail == 0) ? 0 : 1;
}

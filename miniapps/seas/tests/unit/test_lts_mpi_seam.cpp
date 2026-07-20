// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_lts_mpi_seam.cpp — LTS Phase 4a Stage 1 gate: the fault-free bulk
// multi-rate corrector's rank-SEAM flux (ComputeADERClusterSeamFaceFluxRHS)
// must make an np>1 run reproduce the np=1 run to ~1e-10 for a SINGLE-cluster
// (Nc=1) mesh — the P-022 "LTS(1 cluster) vs GTS given the LTS partition" gate.
//
// Why Nc=1: every rank-seam bulk face is then a SAME-CLUSTER (diff-0) seam, the
// only kind Stage 1 handles (a cross-cluster diff-1 seam FAILS LOUD, Stage 2).
// With one cluster the tick loop runs one Correct per sync, driving the seam
// exchange each sync — so this exercises the new ghost-`I` exchange + the
// one-sided shared flux + the matched-collective count without needing a
// multi-rate schedule.
//
// The gate is partition-invariant: each run computes global reductions over the
// bulk state Q after K sync intervals (total energy Σ‖Q_e‖², a CENTROID-weighted
// energy Σ w(centroid_e)·‖Q_e‖² — sensitive to WHICH element holds which value,
// and partition-invariant because the physical centroid is — and the max
// per-element L2 norm, order-independent).  Run at np=1 first: it writes the
// reference file.  Run at np>1: it reads the reference and asserts agreement.
//
// Usage (the Makefile target runs both):
//   mpirun -np 1 seas_test_lts_mpi_seam        # writes the reference
//   mpirun -np 2 seas_test_lts_mpi_seam        # compares against it
// Reference path: $SEAM_REF_FILE (default ./seam_lts_ref.txt).

#include "mfem.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/lts_layout.hpp"
#include "../../dynamic/lts_stepper.hpp"
#include "../../dynamic/lts_bulk_stepper.hpp"
#include "../../domain/boundary_config.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int g_checks = 0, g_fails = 0;
#define CHECK(cond, msg)                                                      \
   do { ++g_checks; if (!(cond)) { ++g_fails;                                \
        std::printf("  FAIL %s:%d  %s\n", __FILE__, __LINE__, msg); } } while (0)

static BoundaryConfig AbsorbingBC()
{
   BoundaryConfig bc;
   for (int i = 1; i <= 6; ++i) { bc.absorbing_attrs.insert(i); }
   bc.fault_attr = 0;   // fault-free (bulk stepper)
   return bc;
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   const int nprocs = Mpi::WorldSize();
   const int rank   = Mpi::WorldRank();

   // --- serial mesh (identical on every rank) -------------------------------
   const int order = 2;
   const real_t lam = 32.04e9, mu = 32.04e9, rho = 2670.0;   // ~ SCEC bulk
   Mesh smesh = Mesh::MakeCartesian3D(4, 2, 2, Element::TETRAHEDRON, 4000.0,
                                      2000.0, 2000.0);
   const int serial_ne = smesh.GetNE();

   // --- ParMesh (default METIS partition at np>1; global elem numbering kept) -
   ParMesh pmesh(MPI_COMM_WORLD, smesh);
   const int ne = pmesh.GetNE();

   // Nc=1: every element in cluster 0 (all rank-seams are diff-0).
   std::vector<int> cluster_id(static_cast<std::size_t>(ne), 0);

   WaveOperator<ParMesh> wave(pmesh, order, lam, mu, rho, AbsorbingBC(),
                              &cluster_id);
   const int ndof_per_el = wave.GetNDof();
   const int ndof_total  = wave.GetScalarNDof();
   const int Nfull       = NUM_STATE * ndof_total;

   // --- rank-local layout (Nc=1): every interior face IntraClusterGTS, every
   //     boundary / rank-seam face Boundary (the seam corrector handles seams) --
   std::vector<LtsFaceSpec> faces;
   faces.reserve(static_cast<std::size_t>(pmesh.GetNumFaces()));
   for (int f = 0; f < pmesh.GetNumFaces(); ++f)
   {
      int e1 = -1, e2 = -1;
      pmesh.GetFaceElements(f, &e1, &e2);
      if (e1 < 0) { continue; }
      LtsFaceSpec fs; fs.face_id = f; fs.elem1 = e1; fs.elem2 = e2;
      fs.is_fault = false;
      faces.push_back(fs);
   }
   LtsLayout layout = BuildLtsLayout(cluster_id, 1, faces);

   // --- initial bulk state: a smooth field of the PHYSICAL coordinates, so it
   //     is identical regardless of how the mesh is partitioned. --------------
   Vector Q(Nfull);
   const FiniteElementSpace &fes = wave.GetFESpace();
   for (int e = 0; e < ne; ++e)
   {
      const FiniteElement *fe = fes.GetFE(e);
      ElementTransformation *Tr = fes.GetElementTransformation(e);
      DenseMatrix coords; Tr->Transform(fe->GetNodes(), coords);
      const int off = e * ndof_per_el;
      for (int i = 0; i < ndof_per_el; ++i)
      {
         const real_t x = coords(0, i), y = coords(1, i), z = coords(2, i);
         for (int c = 0; c < NUM_STATE; ++c)
         {
            Q[c * ndof_total + off + i] =
               std::cos(2.0e-3 * x + 0.2 * c) * std::sin(1.3e-3 * y)
               * std::cos(0.9e-3 * z);
         }
      }
   }

   // Initial global energy — for the "field actually evolved" negative control
   // (a no-op stepper would pass the np==np gate trivially otherwise).
   double init_energy = 0.0;
   for (int i = 0; i < Nfull; ++i) { init_energy += Q[i] * Q[i]; }
   MPI_Allreduce(MPI_IN_PLACE, &init_energy, 1, MPI_DOUBLE, MPI_SUM,
                 MPI_COMM_WORLD);

   // --- drive K sync intervals of the fault-free bulk multi-rate stepper -----
   const int ader_order = 4, K = 8;
   const real_t dt_base = 2.0e-4;   // Nc=1 => one tick per sync, T_s = dt_base
   LtsBulkSyncStepper<ParMesh> stepper(wave, layout, cluster_id, Q, ader_order,
                                       dt_base);
   LtsGlobalMeta meta = ReduceGlobalMeta(
      layout,
      (nprocs > 1)
         ? std::function<void(std::vector<long long>&)>(
              [](std::vector<long long> &buf)
              {
                 MPI_Allreduce(MPI_IN_PLACE, buf.data(),
                               static_cast<int>(buf.size()),
                               MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
              })
         : nullptr);

   for (int s = 0; s < K; ++s)
   {
      auto tab = BuildTickTable(1, dt_base, dt_base, ader_order, meta);
      stepper.SetSyncInterval(static_cast<real_t>(s) * dt_base, dt_base);
      wave.SetTime(static_cast<real_t>(s) * dt_base);
      wave.ResetGhostExchangeCount();
      long long expected = 0;
      for (const auto &tk : tab) { expected += tk.n_collectives; }
      RunSyncInterval(tab, stepper);
      CHECK(nprocs == 1 || wave.GhostExchangeCount() == expected,
            "per-sync ghost-exchange count == tick-table n_collectives");
      CHECK(stepper.BuffersZero(), "accumulate buffers zero at sync point");
      CHECK(Q.CheckFinite() == 0, "Q finite after sync");
   }

   // --- partition-invariant reductions over the final Q ----------------------
   // NOTE: weight the "which-element-holds-what" reduction by the element's
   // PHYSICAL CENTROID, not ParMesh::GetGlobalElementNum — the latter is a
   // rank-blocked parallel numbering that depends on the partition, so it is NOT
   // partition-invariant (it would flag a false mismatch even when the field is
   // bit-identical).  The centroid is physical and partition-invariant.
   double energy = 0.0, cen_energy = 0.0, max_norm = 0.0;
   for (int e = 0; e < ne; ++e)
   {
      const int off = e * ndof_per_el;
      double en = 0.0;
      for (int c = 0; c < NUM_STATE; ++c)
         for (int i = 0; i < ndof_per_el; ++i)
         {
            const double v = Q[c * ndof_total + off + i];
            en += v * v;
         }
      Array<int> vtx; pmesh.GetElementVertices(e, vtx);
      double cx = 0.0, cy = 0.0, cz = 0.0;
      for (int vi = 0; vi < vtx.Size(); ++vi)
      {
         const real_t *X = pmesh.GetVertex(vtx[vi]);
         cx += X[0]; cy += X[1]; cz += X[2];
      }
      const double inv = (vtx.Size() > 0) ? 1.0 / vtx.Size() : 0.0;
      const double w = 1.0e-3 * (cx * inv) + 2.0e-3 * (cy * inv)
                     + 3.0e-3 * (cz * inv);
      energy     += en;
      cen_energy += w * en;
      max_norm    = std::max(max_norm, std::sqrt(en));
   }
   double g_energy = energy, g_cen_energy = cen_energy, g_max = max_norm;
   MPI_Allreduce(MPI_IN_PLACE, &g_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(MPI_IN_PLACE, &g_cen_energy, 1, MPI_DOUBLE, MPI_SUM,
                 MPI_COMM_WORLD);
   MPI_Allreduce(MPI_IN_PLACE, &g_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

   // Negative control: the field must actually have evolved, else a no-op stepper
   // would pass the np==np gate trivially (REVIEW p4a LOW).
   CHECK(std::abs(g_energy - init_energy) > 1e-6 * init_energy,
         "bulk field evolved over K syncs (gate is not trivially green)");

   const char *ref_env = std::getenv("SEAM_REF_FILE");
   const std::string ref_path = ref_env ? ref_env : "seam_lts_ref.txt";

   if (rank == 0)
   {
      std::printf("[seam] np=%d serial_ne=%d after K=%d syncs: "
                  "energy=%.15e cen_energy=%.15e max_elem_norm=%.15e\n",
                  nprocs, serial_ne, K, g_energy, g_cen_energy, g_max);
   }

   if (nprocs == 1)
   {
      if (rank == 0)
      {
         // Stamp serial_ne + K so an np>1 run can reject a stale/mismatched
         // reference (REVIEW p4a LOW) rather than silently trusting it.
         std::ofstream of(ref_path);
         of.precision(17);
         of << serial_ne << " " << K << " "
            << g_energy << " " << g_cen_energy << " " << g_max << "\n";
         std::printf("[seam] wrote np=1 reference to %s\n", ref_path.c_str());
      }
   }
   else
   {
      // np>1: compare against the np=1 reference (must match mesh + K).
      long long r_ne = -1, r_K = -1;
      double r_energy = 0.0, r_cen = 0.0, r_max = 0.0;
      int hr = 0;   // 1 = readable + same (serial_ne, K)
      if (rank == 0)
      {
         std::ifstream inf(ref_path);
         if (inf && (inf >> r_ne >> r_K >> r_energy >> r_cen >> r_max))
         {
            hr = (r_ne == serial_ne && r_K == K) ? 1 : 0;
            if (!hr)
            {
               std::printf("[seam] reference %s is stale: serial_ne=%lld K=%lld "
                           "(expected %d, %d)\n", ref_path.c_str(), r_ne, r_K,
                           serial_ne, K);
            }
         }
      }
      MPI_Bcast(&hr, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&r_energy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&r_cen, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&r_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      CHECK(hr == 1, "np=1 reference present + matches (serial_ne, K) — run np=1 first");
      if (hr == 1)
      {
         auto rel = [](double a, double b)
         { return std::abs(a - b) / (std::abs(b) + 1e-30); };
         const double tol = 1e-10;
         CHECK(rel(g_energy, r_energy) <= tol, "np>1 energy == np=1 energy (1e-10)");
         CHECK(rel(g_cen_energy, r_cen) <= tol,
               "np>1 centroid-weighted energy == np=1 (1e-10)");
         CHECK(rel(g_max, r_max) <= tol, "np>1 max elem norm == np=1 (1e-10)");
         if (rank == 0)
         {
            std::printf("[seam] np=%d vs np=1: d(energy)=%.3e d(cen)=%.3e "
                        "d(max)=%.3e\n", nprocs,
                        rel(g_energy, r_energy), rel(g_cen_energy, r_cen),
                        rel(g_max, r_max));
         }
      }
   }

   int local_fails = g_fails, total_fails = 0, total_checks = 0;
   MPI_Allreduce(&local_fails, &total_fails, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&g_checks, &total_checks, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if (rank == 0)
   {
      std::printf("test_lts_mpi_seam (np=%d): %d/%d checks passed, %d failed.\n",
                  nprocs, total_checks - total_fails, total_checks, total_fails);
   }
   const int rc = (total_fails == 0) ? 0 : 1;
   return rc;
}

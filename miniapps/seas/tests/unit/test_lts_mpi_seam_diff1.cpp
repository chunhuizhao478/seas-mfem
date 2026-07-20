// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_lts_mpi_seam_diff1.cpp — LTS Phase 4a Stage 2 gate: the CROSS-CLUSTER
// (diff-1) rank-seam corrector must make an np>1 run reproduce the np=1 run.
//
// Design of the differential (strong): a fault-free mesh with a 2-cluster split
// at x=mid (cluster 0 fine for x<mid, cluster 1 coarse for x>=mid; maxdiff 1).
//   - np=1: the whole mesh is on rank 0, so the mid boundary is a LOCAL diff-1
//     face handled by the ALREADY-VALIDATED Phase-2 buffer path
//     (ConsumerFine/ProviderCoarseSkip + accumulate buffers).
//   - np=2: an explicit partition puts x<mid on rank 0 and x>=mid on rank 1, so
//     the mid boundary becomes a RANK-SEAM diff-1 face handled by the NEW Stage-2
//     path (fine-side reads the exchanged coarse forecast D(k); coarse-side fills
//     its local seam in-tray from the ghost fine I).
// So np=2==np=1 proves the new rank-seam diff-1 == the validated local diff-1.
//
// Run np=1 first (writes the reference), then np>1 (compares).  Partition-
// invariant CENTROID-weighted reductions; field-evolved control; per-sync ghost-
// exchange counter == tick-table n_collectives; SeamCoarseBuffersZero at sync.

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
   bc.fault_attr = 0;
   return bc;
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   const int nprocs = Mpi::WorldSize();
   const int rank   = Mpi::WorldRank();

   const int order = 2;
   const real_t lam = 32.04e9, mu = 32.04e9, rho = 2670.0;
   // NClust clusters banded along x (adjacent bands differ by 1 => maxdiff 1).
   // Requires np == 1 (reference; bands are LOCAL diff-1 faces) or np == NClust
   // (each rank owns one band => every band boundary is a RANK-SEAM diff-1 face,
   // and a MIDDLE band is simultaneously a provider and a consumer — the B.11
   // 3-cluster-chain scenario at NClust=3).
   const int NClust = std::getenv("SEAM_DIFF1_NC")
                      ? std::atoi(std::getenv("SEAM_DIFF1_NC")) : 2;
   MFEM_VERIFY(NClust >= 2 && (nprocs == 1 || nprocs == NClust),
               "test_lts_mpi_seam_diff1: run with np==1 or np==NClust "
               "(SEAM_DIFF1_NC).");
   const real_t Lx = 2000.0 * NClust, Ly = 2000.0, Lz = 2000.0;
   Mesh smesh = Mesh::MakeCartesian3D(4 * NClust, 2, 2, Element::TETRAHEDRON,
                                      Lx, Ly, Lz);
   const int serial_ne = smesh.GetNE();

   // Serial-element centroid x (for the forced partition at np>1).
   auto centroid_x = [](Mesh &m, int e) -> double
   {
      Array<int> v; m.GetElementVertices(e, v);
      double cx = 0.0;
      for (int i = 0; i < v.Size(); ++i) { cx += m.GetVertex(v[i])[0]; }
      return (v.Size() > 0) ? cx / v.Size() : 0.0;
   };
   const real_t band = Lx / NClust;
   auto band_of = [&](double cx) -> int
   { int b = static_cast<int>(cx / band); return std::min(std::max(b, 0), NClust - 1); };

   std::vector<int> part(static_cast<std::size_t>(serial_ne), 0);
   if (nprocs > 1)
   {
      for (int e = 0; e < serial_ne; ++e)
      { part[e] = band_of(centroid_x(smesh, e)); }   // one band per rank
   }

   ParMesh pmesh(MPI_COMM_WORLD, smesh, nprocs > 1 ? part.data() : nullptr);
   const int ne = pmesh.GetNE();

   // Cluster ids from the LOCAL centroid (partition-invariant): band index in
   // [0, NClust).  Adjacent bands differ by 1.
   std::vector<int> cluster_id(static_cast<std::size_t>(ne), 0);
   for (int e = 0; e < ne; ++e)
   { cluster_id[e] = band_of(centroid_x(pmesh, e)); }

   WaveOperator<ParMesh> wave(pmesh, order, lam, mu, rho, AbsorbingBC(),
                              &cluster_id);
   const int ndof_per_el = wave.GetNDof();
   const int ndof_total  = wave.GetScalarNDof();
   const int Nfull       = NUM_STATE * ndof_total;

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
   LtsLayout layout = BuildLtsLayout(cluster_id, NClust, faces);

   // Initial bulk state (physical-coordinate field, partition-invariant).
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
            Q[c * ndof_total + off + i] =
               std::cos(2.0e-3 * x + 0.2 * c) * std::sin(1.3e-3 * y)
               * std::cos(0.9e-3 * z);
      }
   }
   double init_energy = 0.0;
   for (int i = 0; i < Nfull; ++i) { init_energy += Q[i] * Q[i]; }
   MPI_Allreduce(MPI_IN_PLACE, &init_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   const int ader_order = 4, K = 8;
   const real_t dt_base = 1.0e-4;   // T_s = 2^(NClust-1)*dt_base, 2^(NClust-1) ticks/sync
   LtsBulkSyncStepper<ParMesh> stepper(wave, layout, cluster_id, Q, ader_order,
                                       dt_base);
   stepper.SetExchangeProviderDk(nprocs > 1);   // Stage 2 diff-1 D(k) exchange

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
   meta.exchange_bulk_provider_dk = (nprocs > 1);

   real_t t_lts = 0.0;
   for (int s = 0; s < K; ++s)
   {
      const real_t T_s = dt_base * static_cast<real_t>(1LL << (NClust - 1));
      auto tab = BuildTickTable(NClust, dt_base, T_s, ader_order, meta);
      stepper.SetSyncInterval(t_lts, T_s);
      wave.SetTime(t_lts);
      wave.ResetGhostExchangeCount();
      long long sched = 0;
      for (const auto &tk : tab) { sched += tk.n_collectives; }
      const long long expected = wave.HasSharedFaces() ? sched : 0;
      RunSyncInterval(tab, stepper);
      CHECK(nprocs == 1 || wave.GhostExchangeCount() == expected,
            "per-sync ghost-exchange count == tick-table n_collectives");
      CHECK(stepper.BuffersZero(), "local accumulate buffers zero at sync");
      CHECK(wave.SeamCoarseBuffersZero(), "coarse-side seam in-trays zero at sync");
      CHECK(Q.CheckFinite() == 0, "Q finite after sync");
      t_lts += T_s;
   }

   double energy = 0.0, cen_energy = 0.0, max_norm = 0.0;
   for (int e = 0; e < ne; ++e)
   {
      const int off = e * ndof_per_el;
      double en = 0.0;
      for (int c = 0; c < NUM_STATE; ++c)
         for (int i = 0; i < ndof_per_el; ++i)
         { const double v = Q[c * ndof_total + off + i]; en += v * v; }
      Array<int> vtx; pmesh.GetElementVertices(e, vtx);
      double cx = 0.0, cy = 0.0, cz = 0.0;
      for (int vi = 0; vi < vtx.Size(); ++vi)
      { const real_t *X = pmesh.GetVertex(vtx[vi]); cx += X[0]; cy += X[1]; cz += X[2]; }
      const double inv = (vtx.Size() > 0) ? 1.0 / vtx.Size() : 0.0;
      const double w = 1.0e-3 * (cx * inv) + 2.0e-3 * (cy * inv) + 3.0e-3 * (cz * inv);
      energy += en; cen_energy += w * en; max_norm = std::max(max_norm, std::sqrt(en));
   }
   double g_energy = energy, g_cen = cen_energy, g_max = max_norm;
   MPI_Allreduce(MPI_IN_PLACE, &g_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(MPI_IN_PLACE, &g_cen, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(MPI_IN_PLACE, &g_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

   CHECK(std::abs(g_energy - init_energy) > 1e-6 * init_energy,
         "bulk field evolved over K syncs (gate is not trivially green)");

   const char *ref_env = std::getenv("SEAM_DIFF1_REF_FILE");
   const std::string ref_path = ref_env ? ref_env : "seam_lts_diff1_ref.txt";
   if (rank == 0)
   {
      std::printf("[seam-diff1] np=%d serial_ne=%d K=%d: energy=%.15e cen=%.15e "
                  "max=%.15e\n", nprocs, serial_ne, K, g_energy, g_cen, g_max);
   }

   if (nprocs == 1)
   {
      if (rank == 0)
      {
         std::ofstream of(ref_path);
         of.precision(17);
         of << serial_ne << " " << K << " "
            << g_energy << " " << g_cen << " " << g_max << "\n";
         std::printf("[seam-diff1] wrote np=1 reference to %s\n", ref_path.c_str());
      }
   }
   else
   {
      long long r_ne = -1, r_K = -1;
      double r_energy = 0.0, r_cen = 0.0, r_max = 0.0;
      int hr = 0;
      if (rank == 0)
      {
         std::ifstream inf(ref_path);
         if (inf && (inf >> r_ne >> r_K >> r_energy >> r_cen >> r_max))
         { hr = (r_ne == serial_ne && r_K == K) ? 1 : 0; }
      }
      MPI_Bcast(&hr, 1, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&r_energy, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&r_cen, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&r_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      CHECK(hr == 1, "np=1 reference present + matches (serial_ne, K)");
      if (hr == 1)
      {
         auto rel = [](double a, double b)
         { return std::abs(a - b) / (std::abs(b) + 1e-30); };
         const double tol = 1e-10;
         CHECK(rel(g_energy, r_energy) <= tol, "np>1 energy == np=1 (1e-10)");
         CHECK(rel(g_cen, r_cen) <= tol, "np>1 centroid energy == np=1 (1e-10)");
         CHECK(rel(g_max, r_max) <= tol, "np>1 max elem norm == np=1 (1e-10)");
         if (rank == 0)
            std::printf("[seam-diff1] np=%d vs np=1: d(energy)=%.3e d(cen)=%.3e "
                        "d(max)=%.3e\n", nprocs, rel(g_energy, r_energy),
                        rel(g_cen, r_cen), rel(g_max, r_max));
      }
   }

   int local_fails = g_fails, total_fails = 0, total_checks = 0;
   MPI_Allreduce(&local_fails, &total_fails, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(&g_checks, &total_checks, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if (rank == 0)
      std::printf("test_lts_mpi_seam_diff1 (np=%d): %d/%d checks passed, %d failed.\n",
                  nprocs, total_checks - total_fails, total_checks, total_fails);
   return (total_fails == 0) ? 0 : 1;
}

// LTS Track-A A1: BITWISE ON-vs-OFF gate on the FAULT-INTERLEAVE stepper.
//
// Companion to test_lts_a1_merge_np2.cpp, which gates the BULK stepper.  The
// production TPV104 deck runs LtsFaultSyncStepper, and before this test that class
// had NO coverage anywhere in the repo -- so A1's wiring there was enabled but
// unproven on exactly the path production uses.  This closes that gap.
//
// WHY ZERO FAULT FACES IS THE RIGHT FIXTURE (and not a shortcut): A1 does not touch
// the fault half at all.  Under D-2 fault faces are rank-interior, so the fault half
// fires ZERO seam collectives; everything A1 changed in this stepper is its
// BeforeCorrects/AfterCorrects override and the BULK rank-seam merge it drives.  A
// fault-free fault-stepper therefore exercises 100% of A1's fault-path delta, while
// needing no friction physics, no DOFData population and no nucleation.  The
// friction iterator below is a mock whose Advance is never reached (no fault QPs).
//
// Run: mpirun -np 2 ./seas_test_lts_a1_merge_fault_np2
//      SEAM_DIFF1_NC=3 mpirun -np 3 ./seas_test_lts_a1_merge_fault_np2

#include "mfem.hpp"
#include "../../dynamic/wave_operator.hpp"
#include "../../dynamic/lts_layout.hpp"
#include "../../dynamic/lts_stepper.hpp"
#include "../../dynamic/lts_fault_stepper.hpp"
#include "../../dynamic/friction_iterator.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/lts_bulk_stepper.hpp"
#include "../../domain/boundary_config.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
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


// Deterministic no-op friction iterator.  The fixture has NO fault faces, so
// Advance is never invoked; SetSubSteps/GetDeltaT/GetTimeWeights are still called
// by the stepper's ctor + Predict, so they must behave.
class MockFrictionIterator : public IFrictionIterator
{
public:
   /// The stepper's ctor reads GetDeltaT()/GetTimeWeights() to size the ADER
   /// sub-step quadrature and VERIFIES O >= 1, so seed them like a real iterator:
   /// `order` sub-steps of equal weight.  Values are irrelevant here (no fault
   /// QPs consume them) but the SIZE is load-bearing.
   explicit MockFrictionIterator(int order)
      : dT_(static_cast<std::size_t>(order), real_t(1) / real_t(order)),
        w_(static_cast<std::size_t>(order), real_t(1) / real_t(order))
   { MFEM_VERIFY(order >= 1, "MockFrictionIterator: order must be >= 1."); }

   void SetSubSteps(std::vector<real_t> deltaT,
                    std::vector<real_t> time_weights) override
   { dT_ = std::move(deltaT); w_ = std::move(time_weights); }
   const std::vector<real_t> &GetDeltaT()      const override { return dT_; }
   const std::vector<real_t> &GetTimeWeights() const override { return w_; }
   void Advance(std::vector<DOFData> &, const std::vector<Vector> &,
                const std::vector<std::vector<real_t>> &,
                const std::vector<std::vector<real_t>> &,
                real_t, real_t, real_t *, real_t *,
                const std::function<void(real_t, real_t)> &) override
   {
      // Unreachable with zero fault QPs; fail loudly if the fixture ever changes.
      MFEM_ABORT("MockFrictionIterator::Advance called -- the A1 fault fixture is "
                 "supposed to have NO fault faces.");
   }
   void SetDiagNumLocalFaultQPs(int n) override
   {
      // The fixture is fault-free; the stepper still reports its (zero) count.
      MFEM_VERIFY(n == 0, "A1 fault fixture expects ZERO local fault QPs, got " << n);
      diag_n_ = n;
   }
   FaultFrictionLaw WaveOpLaw() const override { return FaultFrictionLaw::LSW; }
private:
   std::vector<real_t> dT_, w_;
   int diag_n_ = 0;
};

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
               "test_lts_a1_merge_np2: run with np==1 or np==NClust "
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

   const Vector Q0(Q);   // pristine initial state for both A/B runs
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


   // ---- run the identical sync sequence on one path -----------------------
   // Returns the final state; records the total ghost-exchange count.
   auto run_sequence = [&](bool merge, Vector &Q_out, long long &n_exch,
                           long long &n_sched) -> void
   {
      Q_out = Q0;                       // same initial state, byte for byte
      // Fault-interleave stepper with an EMPTY fault set (see the header note).
      MockFrictionIterator iter(ader_order);
      std::vector<DOFData> dof_data;                 // no fault QPs
      std::vector<Vector>  fault_coords;             // no fault QPs
      LtsFaultSyncStepper<ParMesh> st(wave, layout, cluster_id, Q_out, ader_order,
                                      dt_base, iter, dof_data, fault_coords,
                                      /*nuc=*/nullptr);
      st.SetExchangeProviderDk(nprocs > 1);
      st.SetMergeTickExchanges(merge);
      LtsGlobalMeta m = meta;
      m.merge_tick_seam_exchanges = merge;   // table must predict what we fire
      n_exch = 0; n_sched = 0;
      real_t t = 0.0;
      for (int s = 0; s < K; ++s)
      {
         const real_t T_full = dt_base * static_cast<real_t>(1LL << (NClust - 1));
         const real_t T_s = (s == K - 1) ? (T_full - dt_base) : T_full;
         auto tab = BuildTickTable(NClust, dt_base, T_s, ader_order, m);
         st.SetSyncInterval(t, T_s);
         wave.SetTime(t);
         wave.ResetGhostExchangeCount();
         long long sched = 0;
         for (const auto &tk : tab) { sched += tk.n_collectives; }
         RunSyncInterval(tab, st);
         const long long expected = wave.HasSharedFaces() ? sched : 0;
         CHECK(nprocs == 1 || wave.GhostExchangeCount() == expected,
               "per-sync ghost-exchange count == tick-table n_collectives");
         CHECK(st.BuffersZero(), "local accumulate buffers zero at sync");
         CHECK(wave.SeamCoarseBuffersZero(), "coarse-side seam in-trays zero at sync");
         CHECK(Q_out.CheckFinite() == 0, "Q finite after sync");
         n_exch  += wave.GhostExchangeCount();
         n_sched += sched;
         t += T_s;
      }
   };

   Vector Q_off, Q_on;
   long long e_off = 0, e_on = 0, s_off = 0, s_on = 0;
   run_sequence(false, Q_off, e_off, s_off);
   run_sequence(true,  Q_on,  e_on,  s_on);

   // ---- THE GATE: byte-for-byte identical state ---------------------------
   CHECK(Q_off.Size() == Q_on.Size(), "ON/OFF state vectors same size");
   const bool bitwise =
      (Q_off.Size() == Q_on.Size())
      && std::memcmp(Q_off.GetData(), Q_on.GetData(),
                     static_cast<std::size_t>(Q_off.Size()) * sizeof(real_t)) == 0;
   CHECK(bitwise, "A1 merged path is BITWISE IDENTICAL to the per-correct path");
   if (!bitwise && Q_off.Size() == Q_on.Size())
   {
      double worst = 0.0; int at = -1;
      for (int i = 0; i < Q_off.Size(); ++i)
      {
         const double d = std::abs(static_cast<double>(Q_off[i] - Q_on[i]));
         if (d > worst) { worst = d; at = i; }
      }
      std::printf("[rank %d] FIRST DIVERGENCE: max |off-on| = %.17g at index %d\n",
                  rank, worst, at);
   }
   // Identical on EVERY rank, not just this one.
   int local_ok = bitwise ? 1 : 0, all_ok = 0;
   MPI_Allreduce(&local_ok, &all_ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
   CHECK(all_ok == 1, "bitwise identity holds on ALL ranks");

   // ---- and the merge actually merged -------------------------------------
   if (nprocs > 1)
   {
      CHECK(e_on < e_off, "merged path fires FEWER collectives than per-correct");
      CHECK(e_on == s_on && e_off == s_off,
            "both paths fire exactly what their tick table predicted");
      // Per tick the merged path pays 1 (I) + 1 (forecast) regardless of how many
      // clusters correct, so the saving grows with cluster count.
      if (rank == 0)
      {
         std::printf("A1 merge (FAULT stepper): collectives %lld -> %lld over %d syncs "
                     "(NClust=%d, np=%d) = %.2fx fewer\n",
                     e_off, e_on, K, NClust, nprocs,
                     e_on > 0 ? static_cast<double>(e_off) / static_cast<double>(e_on)
                              : 0.0);
      }
   }

   int total_checks = g_checks, total_fails = g_fails;
   MPI_Allreduce(MPI_IN_PLACE, &total_checks, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   MPI_Allreduce(MPI_IN_PLACE, &total_fails,  1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
   if (rank == 0)
   {
      std::printf("test_lts_a1_merge_fault_np2 (np=%d): %d/%d passed, %d failed.\n",
                  nprocs, total_checks - total_fails, total_checks, total_fails);
   }
   return total_fails == 0 ? 0 : 1;
}

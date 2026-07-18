// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_lts_layout.cpp — Appendix B.2 of PLAN_clustered_lts_ader_2026-07-18.md.
//
// Covers the Phase-1 run-side layout (dynamic/lts_layout.{hpp,cpp}):
//   * face-role assignment on a 3-cluster chain,
//   * provider / consumer coarse-side element sets,
//   * dense + complete slot maps,
//   * LOCAL and reduced GLOBAL per-cluster metadata counts,
//   * fault-face bookkeeping (diff-0),
// plus a sanity check of the tick-table generator (dynamic/lts_stepper.hpp,
// Appendix A.3) that Phase 1 also delivers: predict/correct predicates,
// FINE->COARSE correct order, per-cluster correct cadence, dt_step, and the
// per-tick collective count.

#include "../../dynamic/lts_layout.hpp"
#include "../../dynamic/lts_stepper.hpp"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace mfem::seas;

// ---- tiny check framework (same style as test_lts_clustering) -------------
static int g_checks = 0;
static int g_fails  = 0;

#define CHECK(cond)                                                        \
   do {                                                                    \
      ++g_checks;                                                          \
      if (!(cond))                                                         \
      {                                                                    \
         ++g_fails;                                                        \
         std::printf("  FAIL %s:%d  CHECK(%s)\n", __FILE__, __LINE__, #cond); \
      }                                                                    \
   } while (0)

#define CHECK_EQ(a, b)                                                     \
   do {                                                                    \
      ++g_checks;                                                          \
      const long long _va = static_cast<long long>(a);                    \
      const long long _vb = static_cast<long long>(b);                    \
      if (_va != _vb)                                                      \
      {                                                                    \
         ++g_fails;                                                        \
         std::printf("  FAIL %s:%d  %s (%lld) != %s (%lld)\n",            \
                     __FILE__, __LINE__, #a, _va, #b, _vb);               \
      }                                                                    \
   } while (0)

// Count occurrences of a role in a cluster's face list.
static int count_role(const LtsCluster& c, FaceRole r)
{
   int n = 0;
   for (FaceRole f : c.faces) { if (f == r) { ++n; } }
   return n;
}

// Find the role assigned to a given face id within a cluster (or -1 sentinel).
static int role_of_face(const LtsCluster& c, int face_id)
{
   for (std::size_t i = 0; i < c.face_ids.size(); ++i)
   {
      if (c.face_ids[i] == face_id) { return static_cast<int>(c.faces[i]); }
   }
   return -1;
}

// ---------------------------------------------------------------------------
// T1: 3-cluster chain e0(c0)-e1(c1)-e2(c2), faces:
//   f10 = boundary at e0, f0 = (e0,e1), f1 = (e1,e2), f20 = boundary at e2.
// ---------------------------------------------------------------------------
static void test_chain_roles()
{
   std::vector<int> cluster = {0, 1, 2};
   std::vector<LtsFaceSpec> faces = {
      {10, 0, -1, false},   // boundary at e0
      { 0, 0,  1, false},   // interior (e0 fine, e1 coarse)
      { 1, 1,  2, false},   // interior (e1 fine, e2 coarse)
      {20, 2, -1, false},   // boundary at e2
   };

   LtsLayout L = BuildLtsLayout(cluster, 3, faces);

   CHECK_EQ(L.num_clusters, 3);
   CHECK_EQ(L.clusters.size(), 3u);

   // Element lists.
   CHECK_EQ(L.clusters[0].elems.size(), 1u);
   CHECK_EQ(L.clusters[0].elems[0], 0);
   CHECK_EQ(L.clusters[1].elems[0], 1);
   CHECK_EQ(L.clusters[2].elems[0], 2);

   // Face roles from each side's POV.
   // Cluster 0 owns: f0 as ConsumerFine (fine side toward coarser e1) + f10 Boundary.
   CHECK_EQ(role_of_face(L.clusters[0], 0), (int)FaceRole::ConsumerFine);
   CHECK_EQ(role_of_face(L.clusters[0], 10), (int)FaceRole::Boundary);
   // Cluster 1 owns: f0 ProviderCoarseSkip (coarse toward finer e0) + f1 ConsumerFine.
   CHECK_EQ(role_of_face(L.clusters[1], 0), (int)FaceRole::ProviderCoarseSkip);
   CHECK_EQ(role_of_face(L.clusters[1], 1), (int)FaceRole::ConsumerFine);
   // Cluster 2 owns: f1 ProviderCoarseSkip + f20 Boundary.
   CHECK_EQ(role_of_face(L.clusters[2], 1), (int)FaceRole::ProviderCoarseSkip);
   CHECK_EQ(role_of_face(L.clusters[2], 20), (int)FaceRole::Boundary);

   // No GTS or Fault faces in this all-distinct-cluster chain.
   CHECK_EQ(count_role(L.clusters[0], FaceRole::IntraClusterGTS), 0);
   CHECK_EQ(count_role(L.clusters[1], FaceRole::IntraClusterGTS), 0);

   // Neighbour bookkeeping: cluster 1's ConsumerFine face f1 points at e2.
   for (std::size_t i = 0; i < L.clusters[1].face_ids.size(); ++i)
   {
      if (L.clusters[1].face_ids[i] == 1)
      { CHECK_EQ(L.clusters[1].face_nbr[i], 2); }
      if (L.clusters[1].face_ids[i] == 0)
      { CHECK_EQ(L.clusters[1].face_nbr[i], 0); }
   }

   // Coarse-side sets: e1 and e2 are each the coarse side of a diff-1 face.
   CHECK_EQ(L.provider_elems.size(), 2u);
   CHECK_EQ(L.provider_elems[0], 1);
   CHECK_EQ(L.provider_elems[1], 2);
   CHECK_EQ(L.consumer_owner_elems.size(), 2u);
   CHECK_EQ(L.consumer_owner_elems[0], 1);
   CHECK_EQ(L.consumer_owner_elems[1], 2);

   // Dense slot maps: complete over the coarse set, -1 elsewhere.
   CHECK_EQ(L.provider_slot_of_elem.size(), 3u);
   CHECK_EQ(L.provider_slot_of_elem[0], -1);
   CHECK_EQ(L.provider_slot_of_elem[1], 0);
   CHECK_EQ(L.provider_slot_of_elem[2], 1);
   CHECK_EQ(L.buffer_slot_of_elem[0], -1);
   CHECK_EQ(L.buffer_slot_of_elem[1], 0);
   CHECK_EQ(L.buffer_slot_of_elem[2], 1);
   // Slot round-trip: provider_elems[slot] == elem.
   for (std::size_t s = 0; s < L.provider_elems.size(); ++s)
   {
      const int e = L.provider_elems[s];
      CHECK_EQ(L.provider_slot_of_elem[e], (int)s);
   }

   // Metadata.
   CHECK_EQ(L.local_elems_per_cluster[0], 1);
   CHECK_EQ(L.local_elems_per_cluster[1], 1);
   CHECK_EQ(L.local_elems_per_cluster[2], 1);
   CHECK_EQ(L.local_fault_faces_per_cluster[0], 0);
}

// ---------------------------------------------------------------------------
// T2: GTS faces + a fault face (diff 0).  Two same-cluster elements joined by a
// fault face; a plain interior face between equal clusters is IntraClusterGTS.
// ---------------------------------------------------------------------------
static void test_gts_and_fault()
{
   // e0,e1 in cluster 0 (GTS neighbours); e2,e3 in cluster 1; fault between e0,e1.
   std::vector<int> cluster = {0, 0, 1, 1};
   std::vector<LtsFaceSpec> faces = {
      {100, 0, 1, true },   // fault face, both cluster 0 (diff 0)
      {  0, 1, 2, false},   // interior diff-1 (e1 c0 fine, e2 c1 coarse)
      {  1, 2, 3, false},   // interior GTS (both cluster 1)
   };

   LtsLayout L = BuildLtsLayout(cluster, 2, faces);

   // Fault face is recorded once, in cluster 0's fault_faces (not a bulk face).
   CHECK_EQ(L.clusters[0].fault_faces.size(), 1u);
   CHECK_EQ(L.clusters[0].fault_faces[0], 100);
   CHECK_EQ(L.local_fault_faces_per_cluster[0], 1);
   CHECK_EQ(L.local_fault_faces_per_cluster[1], 0);
   // The fault face id must NOT appear in any bulk face list.
   CHECK_EQ(role_of_face(L.clusters[0], 100), -1);

   // GTS face f1 lives once in cluster 1 as IntraClusterGTS.
   CHECK_EQ(role_of_face(L.clusters[1], 1), (int)FaceRole::IntraClusterGTS);

   // Diff-1 face f0: ConsumerFine in cluster 0 (e1 fine), ProviderCoarseSkip in
   // cluster 1 (e2 coarse).  e2 is the sole coarse-side element.
   CHECK_EQ(role_of_face(L.clusters[0], 0), (int)FaceRole::ConsumerFine);
   CHECK_EQ(role_of_face(L.clusters[1], 0), (int)FaceRole::ProviderCoarseSkip);
   CHECK_EQ(L.provider_elems.size(), 1u);
   CHECK_EQ(L.provider_elems[0], 2);
   CHECK_EQ(L.buffer_slot_of_elem[2], 0);
   CHECK_EQ(L.buffer_slot_of_elem[0], -1);

   CHECK_EQ(L.local_elems_per_cluster[0], 2);
   CHECK_EQ(L.local_elems_per_cluster[1], 2);
}

// ---------------------------------------------------------------------------
// T3: global-metadata reduction.  nullptr reducer == serial identity; a mock
// "2-rank" reducer that doubles the buffer must double the counts.
// ---------------------------------------------------------------------------
static void test_global_meta()
{
   std::vector<int> cluster = {0, 1, 2};
   std::vector<LtsFaceSpec> faces = { {0, 0, 1, false}, {1, 1, 2, false} };

   LtsLayout L = BuildLtsLayout(cluster, 3, faces);

   // Serial: global == local.
   LtsGlobalMeta serial = ReduceGlobalMeta(L, nullptr);
   CHECK_EQ(serial.global_elems.size(), 3u);
   CHECK_EQ(serial.global_elems[0], 1);
   CHECK_EQ(serial.global_elems[1], 1);
   CHECK_EQ(serial.global_elems[2], 1);
   CHECK_EQ(serial.global_fault_faces[0], 0);

   // Mock 2-rank: doubling reducer.
   auto doubler = [](std::vector<long long>& b)
   { for (long long& v : b) { v *= 2; } };
   LtsGlobalMeta two = ReduceGlobalMeta(L, doubler);
   CHECK_EQ(two.global_elems[0], 2);
   CHECK_EQ(two.global_elems[1], 2);
   CHECK_EQ(two.global_elems[2], 2);
}

// ---------------------------------------------------------------------------
// T4: tick-table generator (Appendix A.3) on a 3-cluster, T=4*dt_base interval.
// ---------------------------------------------------------------------------
static void test_tick_table()
{
   LtsGlobalMeta meta;
   meta.global_elems.assign(3, 1);
   meta.global_fault_faces.assign(3, 0);   // no fault faces -> no predictor exchanges
   meta.num_state = 9;
   meta.drop_fault_predictor_exchange = true;

   std::vector<LtsTick> tab = BuildTickTable(/*num_clusters=*/3, /*dt_base=*/1.0,
                                             /*T_actual=*/4.0, /*ader_order=*/4,
                                             meta);
   CHECK_EQ(tab.size(), 4u);

   auto has = [](const std::vector<int>& v, int x)
   { for (int e : v) { if (e == x) { return true; } } return false; };

   // Predict sets: c predicts when tick % 2^c == 0.
   CHECK(has(tab[0].predict_clusters, 0) && has(tab[0].predict_clusters, 1)
         && has(tab[0].predict_clusters, 2));
   CHECK_EQ(tab[1].predict_clusters.size(), 1u);   // only c0
   CHECK(has(tab[1].predict_clusters, 0));
   CHECK_EQ(tab[2].predict_clusters.size(), 2u);   // c0,c1
   CHECK(has(tab[2].predict_clusters, 0) && has(tab[2].predict_clusters, 1));
   CHECK_EQ(tab[3].predict_clusters.size(), 1u);   // only c0

   // Correct sets: (tick+1)%2^c==0 OR tick+1==nticks.
   CHECK_EQ(tab[0].correct_clusters.size(), 1u);   // c0
   CHECK(has(tab[0].correct_clusters, 0));
   CHECK_EQ(tab[1].correct_clusters.size(), 2u);   // c0,c1
   CHECK_EQ(tab[2].correct_clusters.size(), 1u);   // c0
   CHECK_EQ(tab[3].correct_clusters.size(), 3u);   // c0,c1,c2 (final sync)

   // FINE->COARSE order (ascending cluster id) in every correct list.
   for (const LtsTick& tk : tab)
   {
      for (std::size_t i = 1; i < tk.correct_clusters.size(); ++i)
      { CHECK(tk.correct_clusters[i - 1] < tk.correct_clusters[i]); }
   }

   // Per-cluster correct cadence over the interval == 2^(maxC - c).
   int corr[3] = {0, 0, 0};
   for (const LtsTick& tk : tab)
   { for (int c : tk.correct_clusters) { ++corr[c]; } }
   CHECK_EQ(corr[0], 4);   // 2^(2-0)
   CHECK_EQ(corr[1], 2);   // 2^(2-1)
   CHECK_EQ(corr[2], 1);   // 2^(2-2)

   // dt_step at tick 0: full untruncated steps.
   CHECK(std::abs(tab[0].dt_step[0] - 1.0) < 1e-14);
   CHECK(std::abs(tab[0].dt_step[1] - 2.0) < 1e-14);
   CHECK(std::abs(tab[0].dt_step[2] - 4.0) < 1e-14);

   // Collective count = num_state * |correct| (fault predictor dropped).
   CHECK_EQ(tab[0].n_collectives, 9);    // 1 correcting cluster
   CHECK_EQ(tab[1].n_collectives, 18);   // 2
   CHECK_EQ(tab[2].n_collectives, 9);    // 1
   CHECK_EQ(tab[3].n_collectives, 27);   // 3
}

// ---------------------------------------------------------------------------
// T5: ragged final interval (T = 3.5*dt_base) — the coarse (c1) step truncates.
// ---------------------------------------------------------------------------
static void test_tick_table_ragged()
{
   LtsGlobalMeta meta;
   meta.global_elems.assign(2, 1);
   meta.global_fault_faces.assign(2, 0);

   // 2 clusters, dt_base=1, T=3.5 -> ticks = ceil(3.5) = 4.
   std::vector<LtsTick> tab = BuildTickTable(2, 1.0, 3.5, 4, meta);
   CHECK_EQ(tab.size(), 4u);

   // Every cluster must correct on the final tick (t+1 == nticks clause).
   CHECK(tab[3].correct_clusters.size() == 2u);

   // c1's step starting at tick 2 spans [2,4] nominally but truncates at 3.5:
   // dt_step[1] at tick 2/3 == 3.5 - 2 = 1.5.
   CHECK(std::abs(tab[2].dt_step[1] - 1.5) < 1e-14);
   CHECK(std::abs(tab[3].dt_step[1] - 1.5) < 1e-14);
   // c0's last step [3,4] also truncates to 0.5.
   CHECK(std::abs(tab[3].dt_step[0] - 0.5) < 1e-14);
}

int main()
{
   test_chain_roles();
   test_gts_and_fault();
   test_global_meta();
   test_tick_table();
   test_tick_table_ragged();

   std::printf("test_lts_layout: %d/%d passed, %d failed.\n",
               g_checks - g_fails, g_checks, g_fails);
   return g_fails == 0 ? 0 : 1;
}

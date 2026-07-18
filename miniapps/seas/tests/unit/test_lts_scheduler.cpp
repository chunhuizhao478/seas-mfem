// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_lts_scheduler.cpp — Appendix B.4 of PLAN_clustered_lts_ader_2026-07-18.md.
//
// Golden tick tables for {2,3} clusters x ticks_per_sync {4,8}; the per-cluster
// correct cadence 2^(maxC-c); an instrumented mock stepper driven by
// RunSyncInterval that verifies FINE->COARSE fill-before-consume, buffers zero
// at sync, and no consume of an unfilled buffer; and a table-walk proving no
// corrector reads a stale (wrong-epoch) provider predict.

#include "../../dynamic/lts_stepper.hpp"

#include <cstdio>
#include <vector>

using namespace mfem::seas;
using mfem::real_t;

static int g_checks = 0, g_fails = 0;
#define CHECK(cond, msg)                                                     \
   do { ++g_checks; if (!(cond)) { ++g_fails;                               \
        std::printf("  FAIL %s:%d  %s\n", __FILE__, __LINE__, msg); } } while (0)

static LtsGlobalMeta meta_no_fault(int nc)
{
   LtsGlobalMeta m;
   m.global_elems.assign(nc, 1);
   m.global_fault_faces.assign(nc, 0);   // no fault predictor exchanges
   m.num_state = 9;
   m.drop_fault_predictor_exchange = true;
   return m;
}

static bool has(const std::vector<int>& v, int x)
{ for (int e : v) { if (e == x) { return true; } } return false; }

// ---------------------------------------------------------------------------
// Mock chain stepper: a linear cluster chain where cluster c is the FINE side of
// the c/(c+1) face and the COARSE side of the (c-1)/c face.  On Correct(c) it
// fills coarse (c+1)'s buffer and consumes its own (c's) buffer.  A consume that
// finds an empty buffer means a corrector ran before its filler (FINE->COARSE
// violated).
// ---------------------------------------------------------------------------
struct MockChainStepper : ILtsClusterStepper
{
   int num_clusters;
   LtsAccumulateBuffers acc;
   int stale_consume = 0;
   std::vector<int> consume_fill;

   explicit MockChainStepper(int nc) : num_clusters(nc)
   { acc.Resize(nc, 1); }

   void Predict(int, real_t) override {}
   void Correct(int c, real_t) override
   {
      const real_t contrib = 1.0;
      if (c + 1 < num_clusters) { acc.AddInto(c + 1, &contrib); }  // fill coarse neighbour
      if (c >= 1)                                                   // consume own buffer
      {
         real_t out = 0.0;
         const int f = acc.ConsumeZero(c, &out);
         consume_fill.push_back(f);
         if (f == 0) { ++stale_consume; }
      }
   }
};

int main()
{
   // ---- T1: golden predict/correct sets, 3 clusters, ticks_per_sync = 4 ------
   {
      auto tab = BuildTickTable(3, 1.0, 4.0, 4, meta_no_fault(3));
      CHECK(tab.size() == 4u, "T1 4 ticks");
      // predicts
      CHECK(has(tab[0].predict_clusters, 0) && has(tab[0].predict_clusters, 1)
            && has(tab[0].predict_clusters, 2), "T1 t0 predict {0,1,2}");
      CHECK(tab[1].predict_clusters.size() == 1u && has(tab[1].predict_clusters, 0),
            "T1 t1 predict {0}");
      CHECK(tab[2].predict_clusters.size() == 2u && has(tab[2].predict_clusters, 1),
            "T1 t2 predict {0,1}");
      CHECK(tab[3].predict_clusters.size() == 1u && has(tab[3].predict_clusters, 0),
            "T1 t3 predict {0}");
      // corrects (FINE->COARSE)
      CHECK(tab[0].correct_clusters == std::vector<int>({0}),       "T1 t0 correct {0}");
      CHECK(tab[1].correct_clusters == std::vector<int>({0, 1}),    "T1 t1 correct {0,1}");
      CHECK(tab[2].correct_clusters == std::vector<int>({0}),       "T1 t2 correct {0}");
      CHECK(tab[3].correct_clusters == std::vector<int>({0, 1, 2}), "T1 t3 correct {0,1,2}");
   }

   // ---- T2: golden sets, 2 clusters, ticks_per_sync = 4 ---------------------
   {
      auto tab = BuildTickTable(2, 1.0, 4.0, 4, meta_no_fault(2));
      CHECK(tab.size() == 4u, "T2 4 ticks");
      CHECK(tab[0].correct_clusters == std::vector<int>({0}),    "T2 t0 correct {0}");
      CHECK(tab[1].correct_clusters == std::vector<int>({0, 1}), "T2 t1 correct {0,1}");
      CHECK(tab[2].correct_clusters == std::vector<int>({0}),    "T2 t2 correct {0}");
      CHECK(tab[3].correct_clusters == std::vector<int>({0, 1}), "T2 t3 correct {0,1}");
   }

   // ---- T3: per-cluster correct cadence == ticks_per_sync / 2^c -------------
   // (For a natural single sync interval ticks_per_sync == 2^maxC this reduces
   // to the plan's 2^(maxC-c); for a longer power-of-2 table it is ticks/2^c.)
   auto cadence_ok = [&](int nc, double T) -> bool
   {
      auto tab = BuildTickTable(nc, 1.0, T, 4, meta_no_fault(nc));
      const int ticks = static_cast<int>(tab.size());
      std::vector<int> corr(nc, 0);
      for (const auto& tk : tab) { for (int c : tk.correct_clusters) { ++corr[c]; } }
      for (int c = 0; c < nc; ++c) { if (corr[c] != ticks / (1 << c)) { return false; } }
      return true;
   };
   CHECK(cadence_ok(2, 4.0), "T3 cadence 2c x T4");   // ticks 4: {4,2}
   CHECK(cadence_ok(3, 4.0), "T3 cadence 3c x T4");   // ticks 4: {4,2,1} == 2^(maxC-c)
   CHECK(cadence_ok(3, 8.0), "T3 cadence 3c x T8");   // ticks 8: {8,4,2}
   CHECK(cadence_ok(2, 8.0), "T3 cadence 2c x T8");   // ticks 8: {8,4}

   // ---- T4: mock stepper via RunSyncInterval — fill/consume + zero-at-sync ---
   for (int nc : {2, 3}) for (double T : {4.0, 8.0})
   {
      auto tab = BuildTickTable(nc, 1.0, T, 4, meta_no_fault(nc));
      MockChainStepper s(nc);
      RunSyncInterval(tab, s);
      char m[80];
      std::snprintf(m, sizeof m, "T4 nc=%d T=%g: no stale consume", nc, T);
      CHECK(s.stale_consume == 0, m);
      std::snprintf(m, sizeof m, "T4 nc=%d T=%g: buffers zero at sync", nc, T);
      CHECK(s.acc.AllZero(), m);
      // Power-of-two interval => every coarse consume closes exactly 2 fine sub-steps.
      bool all_two = true;
      for (int f : s.consume_fill) { if (f != 2) { all_two = false; } }
      std::snprintf(m, sizeof m, "T4 nc=%d T=%g: every consume fill == 2", nc, T);
      CHECK(all_two, m);
   }

   // ---- T5: no corrector reads a stale provider predict ----------------------
   // Walk the table tracking each cluster's last predict tick; on every correct
   // of cluster c, its own last predict must be the CURRENT step-open tick
   // floor(t/2^c)*2^c (the expansion point a consumer integrates from).
   {
      const int nc = 3;
      auto tab = BuildTickTable(nc, 1.0, 8.0, 4, meta_no_fault(nc));
      std::vector<long long> last_predict(nc, -1);
      bool stale = false;
      for (long long t = 0; t < static_cast<long long>(tab.size()); ++t)
      {
         for (int c : tab[t].predict_clusters) { last_predict[c] = t; }
         for (int c : tab[t].correct_clusters)
         {
            const long long period = 1LL << c;
            const long long open = (t / period) * period;
            if (last_predict[c] != open) { stale = true; }
         }
      }
      CHECK(!stale, "T5 every corrector's predict is its current step-open (no stale)");
   }

   std::printf("test_lts_scheduler: %d/%d passed, %d failed.\n",
               g_checks - g_fails, g_checks, g_fails);
   return g_fails == 0 ? 0 : 1;
}

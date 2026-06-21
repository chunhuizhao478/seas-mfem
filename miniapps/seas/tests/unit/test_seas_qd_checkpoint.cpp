// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_seas_qd_checkpoint.cpp — Phase 6 of
// document/spatial_seas_dev/PLAN_spatial_seas_quasidynamic_driver_2026-05-31.md
//
// Serialization round-trip for io/seas_qd_checkpoint.hpp.  Pure (no mesh, no
// MPI, no SAFS coupling), so it runs locally — it is the local half of the
// Phase-6 checkpoint acceptance ("restart reproduces the trajectory TO
// ROUND-OFF"): a write→read must return bit-identical (t, dt, step, n_owned)
// and bit-identical Vectors, using full IEEE precision.  The full restart-
// reproduces-trajectory run is Frontera (SAFS coupling needs the BP5 mesh).
//
// Also checks the refusal paths: missing file, magic mismatch, driver_tag
// mismatch (no cross-driver restart).

#include "mfem.hpp"
#include "../../io/seas_qd_checkpoint.hpp"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

using namespace mfem;
using namespace mfem::seas;

namespace
{
int g_fail = 0;
int g_checks = 0;

void Check(bool ok, const std::string &msg)
{
   g_checks++;
   if (!ok) { g_fail++; std::cerr << "FAILED: " << msg << "\n"; }
}

// Fill v with deliberately precision-demanding values (need 17 sig digits).
void FillTricky(Vector &v, int n, real_t base)
{
   v.SetSize(n);
   for (int i = 0; i < n; ++i)
   {
      v(i) = base + static_cast<real_t>(i + 1) / 3.0 - 1.0 / 7.0
             + 1.234567890123456e-9 * (i + 1);
   }
}

bool VecEqual(const Vector &a, const Vector &b)
{
   if (a.Size() != b.Size()) { return false; }
   for (int i = 0; i < a.Size(); ++i) { if (a(i) != b(i)) { return false; } }
   return true;
}
}  // namespace

int main()
{
   const std::string prefix = "/tmp/seas_qd_ckpt_test";
   const int rank = 0;

   // ---- build a populated checkpoint ----
   const int N = 11;
   SpatialSeasCheckpoint c;
   c.t       = 3.15576e9 + 1.0 / 3.0;     // ~100 yr + a non-representable fraction
   c.dt      = 0.0123456789012345;
   c.step    = 4242;
   c.n_owned = N;
   FillTricky(c.state,   3 * N, 1.0e-9);   // [s_dip, s_strike, psi]*N
   FillTricky(c.a,       N, 0.010);
   FillTricky(c.Dc,      N, 0.004);
   FillTricky(c.eta,     N, 4.6e6);
   FillTricky(c.V_init,  N, 1.0e-9);
   FillTricky(c.sigma_n, N, 50.0e6);
   FillTricky(c.tau_pre, 2 * N, 3.0e7);

   // ---- write + read back ----
   WriteSpatialSeasCheckpoint(prefix, rank, c);
   SpatialSeasCheckpoint r;
   std::string tag;
   const bool ok = ReadSpatialSeasCheckpoint(prefix, rank, r, "spatial_seas", &tag);
   Check(ok, "Read returns true on a valid checkpoint");
   Check(tag == "spatial_seas", "driver_tag trailer round-trips");

   // ---- bit-exact scalars ----
   Check(r.t == c.t, "t bit-exact (full precision)");
   Check(r.dt == c.dt, "dt bit-exact");
   Check(r.step == c.step, "step exact");
   Check(r.n_owned == c.n_owned, "n_owned exact");

   // ---- bit-exact vectors ----
   Check(VecEqual(r.state,   c.state),   "state bit-exact");
   Check(VecEqual(r.a,       c.a),       "a bit-exact");
   Check(VecEqual(r.Dc,      c.Dc),      "Dc bit-exact");
   Check(VecEqual(r.eta,     c.eta),     "eta bit-exact");
   Check(VecEqual(r.V_init,  c.V_init),  "V_init bit-exact");
   Check(VecEqual(r.sigma_n, c.sigma_n), "sigma_n bit-exact");
   Check(VecEqual(r.tau_pre, c.tau_pre), "tau_pre bit-exact");

   // ---- refusal: missing file ----
   SpatialSeasCheckpoint dummy;
   Check(!ReadSpatialSeasCheckpoint("/tmp/seas_qd_ckpt_nonexistent_xyz", rank, dummy),
         "Read returns false on a missing file");

   // ---- refusal: driver_tag mismatch (no cross-driver restart) ----
   Check(!ReadSpatialSeasCheckpoint(prefix, rank, dummy, "some_other_driver"),
         "Read refuses a mismatched driver_tag");

   // ---- refusal: magic mismatch ----
   {
      const std::string bad = "/tmp/seas_qd_ckpt_badmagic";
      std::ofstream out(SpatialSeasCheckpointFilename(bad, rank));
      out << "NOT_THE_MAGIC\nt 0\n";
      out.close();
      Check(!ReadSpatialSeasCheckpoint(bad, rank, dummy),
            "Read refuses a magic mismatch");
      std::remove(SpatialSeasCheckpointFilename(bad, rank).c_str());
   }

   std::remove(SpatialSeasCheckpointFilename(prefix, rank).c_str());

   std::cout << "=== test_seas_qd_checkpoint ===\n"
             << "  checks: " << g_checks << ", failures: " << g_fail << "\n"
             << (g_fail == 0 ? "PASS\n" : "FAIL\n");
   return (g_fail == 0) ? 0 : 1;
}

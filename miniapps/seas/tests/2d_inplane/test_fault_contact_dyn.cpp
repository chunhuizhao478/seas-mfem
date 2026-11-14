// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/2d_inplane/test_fault_contact_dyn.cpp — fault opening version B,
// phase B2 (plan 2026-11-05, Section 4): the dynamic contact projection in
// NitscheWaveOperator::CloseDOF, on the in-plane slip-weakening fault of
// test_sigma_floor.cpp (p = 1 and 2).
//
//   D1  the sign s of plan eq. 4.3, pinned through the ASSEMBLED operator
//       (the plan leaves s to this test): with u = 0 frozen, a gap
//       g = delta at every fault DOF (d = delta x oriented normal) moves
//       the elastic sigma_eff = -(T . n) by slope x delta per DOF.  The
//       projection (4.4) as implemented assumes the plan's orientation,
//       sigma^{n+1} = sigma^X + gamma V_g with V_g = V . n, i.e. slope
//       = +2 gamma / dt = +sigma_F.  A NEGATIVE measured slope means the
//       projection is ill-posed as coded: a major plan deviation — stop
//       and report, do not trust D2.
//   D2  CloseDOF branches (only meaningful when D1 measures s = +1):
//       a  tension at gap 0 opens: sigma_eff returns 0, V . n = -sigma^X
//          / gamma, and the shear closes to X + gamma V_t = 0;
//       b  compression at gap 0 stays closed, bit for bit the slip rate
//          and sigma_eff of an operator without any opening mode;
//       c  compression with a gap g* closes it: V . n = -2 g*/dt exactly,
//          sigma_eff returns the projected sigma^X - 2 gamma g*/dt;
//       d  a large gap under mild compression stays open: sigma_eff = 0,
//          V . n = -sigma^X/gamma (the same open formula);
//       e  complementarity sweep: g^{n+1} >= 0, sigma^{n+1} >= 0 and
//          g^{n+1} sigma^{n+1} = 0 across tension/compression x gaps.
//   D3  version A untouched: SetShearFreeOpening still gives V = -X/gamma
//       with a GLUED normal (V . n = 0) and sigma_eff = the trial value.
//   D4  [cycle.relaxation] opening = "relaxed" (SetRelaxedContact with the
//       normal-stress hook holding a fixed r per DOF): CloseDOF at a trial
//       sigma^X is CloseDOF of the plain contact operator at sigma^X - r
//       (the same V; the returned bulk sigma_eff is that one's plus r; the
//       same open flag), over tension/compression x gaps x r of both
//       signs; and with r = 0 it is the plain contact operator bit for bit.
//
// Usage (from miniapps/seas): mpirun -np N ./seas_test_fault_contact_dyn

#include "mfem.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../solver/nitsche_cfl.hpp"
#include "../../friction/slip_weakening_closure.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../cycle/normal_stress_relaxation.hpp"

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   bool ok_ = (cond); \
   int bad_ = ok_ ? 0 : 1; \
   MPI_Allreduce(MPI_IN_PLACE, &bad_, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); \
   if (bad_ == 0) { if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (!ok_) std::cout << "  FAILED [" << __LINE__ << ", rank " \
                                            << g_rank << "]: " << msg << "\n"; } \
} while(0)

static const real_t kL = 30.0e3;
static const real_t kMuS = 0.677, kMuD = 0.525, kDc = 0.40;
static const real_t kTau0 = 70.0e6, kSigma0 = 120.0e6;
static const real_t kLf = 15.0e3;
static const int kNf = 24;

/// [-L, L]^2 in kNf x kNf cells split into 4 triangles through the centre;
/// walls attr 5, the grid edges on y = 0 with |x| < kLf attr 3 (the mesh of
/// test_sigma_floor.cpp / test_normal_stress_hook.cpp).
static Mesh *FaultBox()
{
   const int nx = kNf, ny = kNf;
   const real_t h = 2 * kL / kNf, x0 = -kL, y0 = -kL;
   const int ngv = (nx + 1) * (ny + 1);
   Mesh *m = new Mesh(2, ngv + nx * ny, 4 * nx * ny, 2 * (nx + ny), 2);
   auto gv = [&](int i, int j) { return i + (nx + 1) * j; };
   auto cv = [&](int i, int j) { return ngv + i + nx * j; };
   for (int j = 0; j <= ny; j++)
      for (int i = 0; i <= nx; i++) { m->AddVertex(x0 + i * h, y0 + j * h); }
   for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++) { m->AddVertex(x0 + (i + 0.5) * h, y0 + (j + 0.5) * h); }
   for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++)
      {
         const int c = cv(i, j);
         m->AddTriangle(gv(i, j), gv(i + 1, j), c, 1);
         m->AddTriangle(gv(i + 1, j), gv(i + 1, j + 1), c, 1);
         m->AddTriangle(gv(i + 1, j + 1), gv(i, j + 1), c, 1);
         m->AddTriangle(gv(i, j + 1), gv(i, j), c, 1);
      }
   for (int i = 0; i < nx; i++)
   {
      m->AddBdrSegment(gv(i, 0), gv(i + 1, 0), 5);
      m->AddBdrSegment(gv(i + 1, ny), gv(i, ny), 5);
   }
   for (int j = 0; j < ny; j++)
   {
      m->AddBdrSegment(gv(nx, j), gv(nx, j + 1), 5);
      m->AddBdrSegment(gv(0, j + 1), gv(0, j), 5);
   }
   const int jf = ny / 2;   // y = 0
   for (int i = 0; i < nx; i++)
   {
      if (std::abs(-kL + (i + 0.5) * (2 * kL / kNf)) < kLf)
      {
         m->AddBdrSegment(gv(i, jf), gv(i + 1, jf), 3);
      }
   }
   m->FinalizeTopology();
   m->Finalize(false, true);
   m->SetAttributes();
   return m;
}

static NitscheFaultSetup FaultSetup(const IFaultClosure &cl)
{
   NitscheFaultSetup fs;
   fs.closure = &cl;
   const real_t ref[3] = {0.0, -1.0, 0.0}, up[3] = {0.0, 0.0, 1.0};
   for (int d = 0; d < 3; d++)
   {
      fs.orientation.ref_normal[d] = ref[d];
      fs.orientation.up[d] = up[d];
      fs.up[d] = up[d];
   }
   fs.orientation.fault_attrs = {3};
   fs.prestress = [](const real_t[3], real_t s[9]) {
      for (int i = 0; i < 9; i++) { s[i] = 0.0; }
      s[1] = s[3] = kTau0;
      s[4] = -kSigma0;
   };
   fs.parameters = [](const real_t[3], FaultPointState &st) {
      st.mu_s = kMuS;
      st.mu_d = kMuD;
      st.Dc = kDc;
      st.cohesion = 0.0;
   };
   return fs;
}

/// The background traction of DOF j moved along the normal to sigma_n_eff = s.
static void MovedTraction(const NitscheWaveOperator &op, const Vector &S, int j,
                          real_t s, real_t TS[3])
{
   const real_t *nv = &op.FaultFrames()(9 * j);
   real_t sig = 0.0;
   for (int c = 0; c < 3; c++) { sig -= S(3 * j + c) * nv[c]; }
   for (int c = 0; c < 3; c++) { TS[c] = S(3 * j + c) + (sig - s) * nv[c]; }
}

static real_t Dot3(const real_t *a, const real_t *b)
{ return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();

   for (int p : {1, 2})
   {
      if (g_rank == 0) { std::cout << "\n=== B2: contact projection, p = " << p << " ===\n"; }
      std::unique_ptr<Mesh> m(FaultBox());
      ParMesh pm(MPI_COMM_WORLD, *m);
      m.reset();
      BoundaryConfig bc;
      bc.fault_attr = 3;
      bc.absorbing_attrs = {5};
      MaterialField mat = MaterialField::MakeConstant(
         TPV102Params::lambda, TPV102Params::mu, TPV102Params::rho);
      SlipWeakeningClosure cl;
      NitscheFaultSetup f_ct = FaultSetup(cl), f_pl = FaultSetup(cl),
                        f_va = FaultSetup(cl);
      NitscheWaveOperator ct(pm, p, mat, bc, DomainConfig{}, &f_ct);
      NitscheWaveOperator pl(pm, p, mat, bc, DomainConfig{}, &f_pl);
      NitscheWaveOperator va(pm, p, mat, bc, DomainConfig{}, &f_va);
      ct.SetContactOpening(true);
      va.SetShearFreeOpening(true);
      const real_t dt = ComputeMaxDt(ct.Sample(), pm, p,
                                     0.5 * NitscheStableConstant(2, p));
      for (NitscheWaveOperator *op : {&ct, &pl, &va}) { op->SetTimeStep(dt); }
      const int no = ct.NumOwnedFaultDOFs();

      // ---- D1: the one-step slope of sigma_eff in the gap (pins s) ----
      // u = 0 frozen; d = delta x n per DOF; T from the assembled operator.
      const real_t delta = 1.0e-4;
      Vector U(ct.Height()), d0(3 * no), dg(3 * no), T0, T1;
      U = 0.0;
      d0 = 0.0;
      for (int j = 0; j < no; j++)
      {
         const real_t *nv = &ct.FaultFrames()(9 * j);
         for (int c = 0; c < 3; c++) { dg(3 * j + c) = delta * nv[c]; }
      }
      ct.TrialTraction(U, d0, T0);
      ct.TrialTraction(U, dg, T1);
      real_t rmin = 1e300, rmax = -1e300;
      for (int j = 0; j < no; j++)
      {
         const real_t *nv = &ct.FaultFrames()(9 * j);
         const real_t s0 = -Dot3(&T0(3 * j), nv), s1 = -Dot3(&T1(3 * j), nv);
         const real_t slope = (s1 - s0) / delta;
         const real_t r = slope / ct.FaultSigmaF()(j);   // +1 expected
         rmin = std::min(rmin, r);
         rmax = std::max(rmax, r);
      }
      MPI_Allreduce(MPI_IN_PLACE, &rmin, 1, MPITypeMap<real_t>::mpi_type,
                    MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &rmax, 1, MPITypeMap<real_t>::mpi_type,
                    MPI_MAX, MPI_COMM_WORLD);
      if (g_rank == 0)
      {
         std::cout << "  D1: slope / sigma_F in [" << rmin << ", " << rmax
                   << "] (the plan's s = +1 means +1)\n";
      }
      TEST_ASSERT(rmin > 0.0, "D1 p = " << p << ": the gap slope of sigma_eff "
                  "through the assembled operator is POSITIVE (plan eq. 4.3, "
                  "s = +1; measured slope/sigma_F in [" << rmin << ", " << rmax
                  << "])");
      TEST_ASSERT(std::abs(rmin - 1.0) < 0.5 && std::abs(rmax - 1.0) < 0.5,
                  "D1 p = " << p << ": the slope is sigma_F = 2 gamma / dt "
                  "to 50% (plan eq. 4.3 cross-check of gamma)");
      const bool s_plus = rmin > 0.0;

      // ---- D2: CloseDOF branches ----
      if (!s_plus)
      {
         if (g_rank == 0)
         {
            std::cout << "  D2 SKIPPED: s = -1 measured; the projection (4.4) "
                      "does not transfer - STOP AND REPORT (plan Section 4)\n";
         }
      }
      else
      {
         Vector S;
         ct.BackgroundTraction(0.0, S);
         long bit_bad = 0;
         real_t e_open = 0.0, e_close = 0.0, e_shear = 0.0;
         int n_open = 0, n_closed = 0;
         for (int j = 0; j < no; j++)
         {
            const real_t *nv = &ct.FaultFrames()(9 * j);
            const real_t *dv = &ct.FaultFrames()(9 * j + 3);
            const real_t *sv = &ct.FaultFrames()(9 * j + 6);
            const real_t gamma = ct.ClosureGamma(j);
            real_t TS[3], Va[3], Vb[3], X[2], Xb[2], sa, sb;
            bool open = false;

            // a: tension at gap 0 opens
            const real_t st_a = -3.0e6;
            MovedTraction(ct, S, j, st_a, TS);
            ct.CloseDOF(j, TS, 0.0, Va, sa, X, nullptr, 0.0, 1.0, 0.0, &open, 0.0);
            if (!open || sa != 0.0) { bit_bad += 1000000; }
            n_open += open ? 1 : 0;
            e_open = std::max(e_open, std::abs(Dot3(Va, nv) - (-st_a / gamma)) /
                              (std::abs(st_a) / gamma));
            e_shear = std::max(e_shear, std::abs(X[0] + gamma * Dot3(Va, dv)) /
                               std::max(real_t(1.0), std::abs(X[0])));
            e_shear = std::max(e_shear, std::abs(X[1] + gamma * Dot3(Va, sv)) /
                               std::max(real_t(1.0), std::abs(X[1])));

            // b: compression at gap 0 = the plain operator, bit for bit
            const real_t st_b = 20.0e6;
            MovedTraction(ct, S, j, st_b, TS);
            ct.CloseDOF(j, TS, 0.0, Va, sa, X, nullptr, 0.0, 1.0, 0.0, &open, 0.0);
            pl.CloseDOF(j, TS, 0.0, Vb, sb, Xb);
            if (open) { bit_bad += 1000; }
            n_closed += open ? 0 : 1;
            for (int c = 0; c < 3; c++) { if (Va[c] != Vb[c]) { bit_bad++; } }
            if (sa != sb || X[0] != Xb[0] || X[1] != Xb[1]) { bit_bad++; }

            // c: compression with a gap g* closes it this step
            const real_t gs = 1.0e-3;
            ct.CloseDOF(j, TS, 0.0, Va, sa, X, nullptr, 0.0, 1.0, 0.0, &open, gs);
            if (open) { bit_bad += 1000; }
            e_close = std::max(e_close, std::abs(Dot3(Va, nv) - (-2.0 * gs / dt)) /
                               (2.0 * gs / dt));
            e_close = std::max(e_close, std::abs(sa - (st_b - 2.0 * gamma * gs / dt)) /
                               st_b);

            // d: a large gap under mild compression stays open
            const real_t st_d = 1.0e5, gd = 1.0;
            MovedTraction(ct, S, j, st_d, TS);
            ct.CloseDOF(j, TS, 0.0, Va, sa, X, nullptr, 0.0, 1.0, 0.0, &open, gd);
            if (!open || sa != 0.0) { bit_bad += 1000000; }
            e_open = std::max(e_open, std::abs(Dot3(Va, nv) - (-st_d / gamma)) /
                              (st_d / gamma));

            // e: complementarity sweep
            for (real_t st : {-5.0e6, -1.0e5, 0.0, 1.0e5, 5.0e6, 2.0e7})
            {
               for (real_t g0 : {0.0, 1.0e-6, 1.0e-3})
               {
                  MovedTraction(ct, S, j, st, TS);
                  ct.CloseDOF(j, TS, 0.0, Va, sa, X, nullptr, 0.0, 1.0, 0.0,
                              &open, g0);
                  const real_t gn = g0 + 0.5 * dt * Dot3(Va, nv);
                  if (gn < -1e-12 || sa < -1e-6 ||
                      std::abs(gn * sa) > 1e-6) { bit_bad += 1000000; }
               }
            }
         }
         MPI_Allreduce(MPI_IN_PLACE, &bit_bad, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
         real_t em[3] = {e_open, e_close, e_shear};
         MPI_Allreduce(MPI_IN_PLACE, em, 3, MPITypeMap<real_t>::mpi_type,
                       MPI_MAX, MPI_COMM_WORLD);
         long nn[2] = {n_open, n_closed};
         MPI_Allreduce(MPI_IN_PLACE, nn, 2, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
         TEST_ASSERT(nn[0] > 0 && nn[1] > 0 && (bit_bad / 1000000) == 0,
                     "D2 p = " << p << ": branch flags and sigma_eff ends (open: "
                     "sigma 0; closed under compression; complementarity "
                     "g sigma = 0 over the sweep); " << nn[0] << " open, "
                     << nn[1] << " closed DOF cases");
         TEST_ASSERT((bit_bad % 1000000) == 0,
                     "D2b p = " << p << ": compression at gap 0 = the plain "
                     "operator bit for bit (V, sigma_eff, X)");
         TEST_ASSERT(em[0] <= 1e-12, "D2a/d p = " << p << ": open V . n = "
                     "-sigma^X/gamma to 1e-12 (rel " << em[0] << ")");
         TEST_ASSERT(em[1] <= 1e-12, "D2c p = " << p << ": closing V . n = "
                     "-2 g*/dt and the projected sigma to 1e-12 (rel " << em[1] << ")");
         TEST_ASSERT(em[2] <= 1e-9, "D2a p = " << p << ": open shear closes to "
                     "X + gamma V_t = 0 (rel " << em[2] << ")");
      }

      // ---- D3: version A untouched ----
      {
         Vector S;
         va.BackgroundTraction(0.0, S);
         real_t e_glue = 0.0, e_sig = 0.0, e_vt = 0.0;
         for (int j = 0; j < no; j++)
         {
            const real_t *nv = &va.FaultFrames()(9 * j);
            const real_t *dv = &va.FaultFrames()(9 * j + 3);
            const real_t gamma = va.ClosureGamma(j);
            real_t TS[3], V[3], X[2], sig;
            bool open = false;
            MovedTraction(va, S, j, -3.0e6, TS);
            va.CloseDOF(j, TS, 0.0, V, sig, X, nullptr, 0.0, 1.0, 0.0, &open);
            e_glue = std::max(e_glue, std::abs(Dot3(V, nv)));
            e_sig = std::max(e_sig, std::abs(sig - (-3.0e6)));
            e_vt = std::max(e_vt, std::abs(Dot3(V, dv) - (-X[0] / gamma)));
            if (!open) { e_glue += 1.0; }
         }
         real_t em[3] = {e_glue, e_sig, e_vt};
         MPI_Allreduce(MPI_IN_PLACE, em, 3, MPITypeMap<real_t>::mpi_type,
                       MPI_MAX, MPI_COMM_WORLD);
         TEST_ASSERT(em[0] <= 1e-9 && em[1] <= 1e-6 && em[2] <= 1e-9,
                     "D3 p = " << p << ": version A's open rule unchanged "
                     "(glued normal " << em[0] << ", sigma_eff = trial to "
                     << em[1] << " Pa abs, V_t = -X/gamma " << em[2] << ")");
      }

      // ---- D4: the contact condition on the relaxed stress ----
      if (s_plus)
      {
         NitscheFaultSetup f_rc = FaultSetup(cl), f_r0 = FaultSetup(cl);
         NitscheWaveOperator rc(pm, p, mat, bc, DomainConfig{}, &f_rc);
         NitscheWaveOperator r0(pm, p, mat, bc, DomainConfig{}, &f_r0);
         Vector rv(no), zv(no);
         zv = 0.0;
         for (int j = 0; j < no; j++) { rv(j) = (j % 2 == 0) ? -8.0e6 : 4.0e6; }
         cycle::FrozenNormalStressReduction mod(rv), mod0(zv);
         for (NitscheWaveOperator *op : {&rc, &r0})
         {
            op->SetContactOpening(true);
            op->SetTimeStep(dt);
         }
         rc.SetNormalStressModifier(&mod);
         r0.SetNormalStressModifier(&mod0);
         rc.SetRelaxedContact(true);
         r0.SetRelaxedContact(true);
         Vector S;
         ct.BackgroundTraction(0.0, S);
         long bad_flag = 0, bad_bits = 0, n_open = 0, n_closed = 0;
         real_t eV = 0.0, eS = 0.0;
         for (int j = 0; j < no; j++)
         {
            const real_t gamma = ct.ClosureGamma(j);
            for (real_t st : {-5.0e6, -1.0e5, 0.0, 1.0e5, 5.0e6, 2.0e7})
            {
               for (real_t g0 : {0.0, 1.0e-6, 1.0e-3})
               {
                  real_t TS[3], TR[3], Va[3], Vb[3], Xa[2], Xb[2], sa, sb;
                  bool oa = false, ob = false;
                  // relaxed trial st: the bulk's trial is st + r
                  MovedTraction(ct, S, j, st, TS);
                  MovedTraction(ct, S, j, st + rv(j), TR);
                  ct.CloseDOF(j, TS, 0.0, Va, sa, Xa, nullptr, 0.0, 1.0, 0.0, &oa, g0);
                  rc.CloseDOF(j, TR, 0.0, Vb, sb, Xb, nullptr, 0.0, 1.0, 0.0, &ob, g0);
                  if (oa != ob) { bad_flag++; }
                  (oa ? n_open : n_closed)++;
                  const real_t vs = std::max(real_t(1e-12), std::sqrt(Dot3(Va, Va)));
                  for (int c = 0; c < 3; c++)
                  {
                     eV = std::max(eV, std::abs(Va[c] - Vb[c]) / vs);
                  }
                  eS = std::max(eS, std::abs(sb - (sa + rv(j))) /
                                std::max(real_t(1.0e6), gamma * vs));
                  // r = 0: the plain contact operator, bit for bit
                  real_t Vc[3], Xc[2], sc;
                  bool oc = false;
                  r0.CloseDOF(j, TS, 0.0, Vc, sc, Xc, nullptr, 0.0, 1.0, 0.0, &oc, g0);
                  if (oc != oa || sc != sa || Vc[0] != Va[0] || Vc[1] != Va[1] ||
                      Vc[2] != Va[2] || Xc[0] != Xa[0] || Xc[1] != Xa[1]) { bad_bits++; }
               }
            }
         }
         long nn[4] = {bad_flag, bad_bits, n_open, n_closed};
         MPI_Allreduce(MPI_IN_PLACE, nn, 4, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
         real_t em[2] = {eV, eS};
         MPI_Allreduce(MPI_IN_PLACE, em, 2, MPITypeMap<real_t>::mpi_type,
                       MPI_MAX, MPI_COMM_WORLD);
         TEST_ASSERT(nn[0] == 0 && nn[2] > 0 && nn[3] > 0 && em[0] <= 1e-9 && em[1] <= 1e-9,
                     "D4 p = " << p << ": relaxed contact at sigma^X = plain contact at "
                     "sigma^X - r (flags; V rel " << em[0] << "; bulk sigma_eff = that "
                     "one's + r, rel " << em[1] << "); " << nn[2] << " open, " << nn[3]
                     << " closed cases");
         TEST_ASSERT(nn[1] == 0, "D4 p = " << p << ": with r = 0 the plain contact "
                     "operator bit for bit (V, sigma_eff, X, flag)");
      }
   }

   int failed = num_failed;
   MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   if (g_rank == 0)
   {
      std::cout << "\n" << (failed > 0 ? "FAILED" : "PASSED") << " ("
                << num_tests << " checks, " << failed << " failed)\n";
   }
   return failed > 0 ? 1 : 0;
}

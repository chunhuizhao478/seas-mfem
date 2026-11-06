// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/2d_inplane/test_fault_contact_qd.cpp — fault opening version B,
// phases B3 and B4 (plan 2026-11-05, Sections 5 and 6): the QD contact solve
// (OpenSetQD in contact mode), the active-set rule with hysteresis, the
// tapered floor, and the B4 gap carriers (the bridge, checkpoint V7), on
// the miniature BP3 layout of test_fault_contact.cpp
// (ShearBox, dip 60, p = 1 and 2, stage operators with slip_comps = 3 and
// the dense Green's function, the production lambdas of
// QDPhase::SetupOpening_):
//
//   Q1  the solve (plan eq. 5.3): with a tensile patch in the normal
//       pre-stress, one stage drives tau_pre + T and sigma_pre + T_n at the
//       open DOFs to round-off, with g > 0 (the measured orientation: an
//       opening raises sigma, so relieving tension needs separation).  A
//       second stage at the solved state changes nothing (idempotent).
//   Q2  the settle loop (plan 5.4): from all closed under the same patch,
//       stage + rule until nothing flips; afterwards g >= 0 everywhere,
//       sigma >= -sigma_tol at closed DOFs, g sigma = 0 to round-off, the
//       strongly tensile DOFs are open, and the set is stable.
//   Q3  hysteresis (plan 5.5, D-B7): a patch inside the band (|sigma| <
//       sigma_tol) flips nothing; beyond the band it opens; a strong
//       reversal closes it again (the solved g goes negative), the gap is
//       projected to 0 and a repeat settle flips nothing.
//   Q4  off = untouched: with the open set empty, the stage hook leaves the
//       traction bit for bit and no solve runs.
//   Q5  SetGap: a handed-over gap reads back, enters the next stage (the
//       carried-gap columns), and a SetFlags close zeroes it.
//   Q6  QDUnitColumn: linear (col(e_i + e_j) = col(e_i) + col(e_j) to
//       solver round-off) and consistent with the dense G on a tangential
//       unit slip.
//   Q7  the tapered floor (plan eq. 6.2): RSLawSigmaN with taper off equals
//       RSFlooredSigmaN bit for bit (value and status); with taper on the
//       arm is min(floor, max(sigma, 0)) — identical to the hard floor at
//       sigma >= floor, equal to sigma in the band 0 < sigma < floor (for a
//       relaxed point), 0 at sigma <= 0, and continuous at both edges.
//       The whole-law form (r < 0, the unclamping relaxation of version
//       A's T6 flips): = eq. 6.2 bit for bit (value and status) for r >= 0
//       and for sigma >= floor; (sigma - r) w below the floor with r < 0,
//       never below sigma; continuous; 0 at and below opening.
//   Q8  a valid contact deck (D0 + mode = "contact", a floor, no
//       check_every) loads: Contact(), sigma_tol_pa and floor_taper
//       defaults.
//   Q9  the bridge with the gap (B4, cycle/state_bridge.hpp): on synthetic
//       frames (eps = +-1, rotated triads) ToNitsche with a gap puts g on
//       d's normal component (d . eps n_k = g), the tangential pairs read
//       back through ToQD unchanged, GapFromGlobal returns the gap, and V
//       keeps no normal part; without a gap d is ToGlobal's alone, bit for
//       bit, and its normal part is round-off (an ulp, not exact zero: an
//       FMA-contracting compiler rounds the dot differently).
//   Q10 checkpoint FD_CYCLE_V7 (B4, cycle/checkpoint.hpp): the open flags
//       and the gap round-trip bit for bit; a shear-free info still writes
//       V6.
//
// Usage (from miniapps/seas; PETSc/MUMPS needed):
//   mpirun -np N ./seas_test_fault_contact_qd

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../fault/fault_frame_field.hpp"
#include "../../constitutive/linear_elastic.hpp"
#include "../../cycle/fault_opening.hpp"
#include "../../cycle/cycle_config.hpp"
#include "../../cycle/state_bridge.hpp"
#include "../../cycle/checkpoint.hpp"
#include "../../friction/device_rate_state.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const real_t kLambda = 32.04e9, kMu = 32.04e9;
static const real_t kLs = 20e3, kD = 20e3, kWf = 12e3;
static const int kNs = 8, kNy = 10;

static real_t RMax(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, c); return v; }
static long RSuml(long v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_LONG, MPI_SUM, c); return v; }

// The BP3 layout in miniature (test_2d_fault_layout.cpp): a grid (s, y)
// mapped by x = s - y cot(dip); the grid line s = 0 is the fault, attr 3
// for depth -y < Wf, attr 5 below; top attr 1, other sides attr 5.
static Mesh *ShearBox(real_t dip_deg)
{
   const real_t c = (dip_deg == 90.0) ? 0.0 : 1.0 / std::tan(dip_deg * M_PI / 180.0);
   const int nvx = kNs + 1, nvy = kNy + 1;
   Mesh *m = new Mesh(2, nvx * nvy, 2 * kNs * kNy, 0, 2);
   for (int j = 0; j < nvy; j++)
      for (int i = 0; i < nvx; i++)
      {
         const real_t s = -kLs + 2.0 * kLs * i / kNs, y = -kD + kD * j / kNy;
         m->AddVertex(s - y * c, y);
      }
   auto v = [&](int i, int j) { return i + nvx * j; };
   for (int j = 0; j < kNy; j++)
      for (int i = 0; i < kNs; i++)
      {
         m->AddTriangle(v(i, j), v(i + 1, j), v(i + 1, j + 1), 1);
         m->AddTriangle(v(i, j), v(i + 1, j + 1), v(i, j + 1), 1);
      }
   for (int i = 0; i < kNs; i++)
   {
      m->AddBdrSegment(v(i, 0), v(i + 1, 0), 5);
      m->AddBdrSegment(v(i, kNy), v(i + 1, kNy), 1);
   }
   for (int j = 0; j < kNy; j++)
   {
      m->AddBdrSegment(v(0, j), v(0, j + 1), 5);
      m->AddBdrSegment(v(kNs, j), v(kNs, j + 1), 5);
      const real_t ymid = -kD + kD * (j + 0.5) / kNy;
      m->AddBdrSegment(v(kNs / 2, j), v(kNs / 2, j + 1), (-ymid < kWf) ? 3 : 5);
   }
   m->FinalizeTopology();
   m->Finalize(false, true);
   m->SetAttributes();
   return m;
}

static std::shared_ptr<const FaultFrameConfig> FrameConfig(const real_t ref[3],
                                                           const real_t up[3])
{
   auto c = std::make_shared<FaultFrameConfig>();
   for (int d = 0; d < 3; d++)
   {
      c->orientation.ref_normal[d] = ref[d];
      c->orientation.up[d] = up[d];
      c->up[d] = up[d];
   }
   c->orientation.fault_attrs = {3};
   return c;
}

// Linear in t, so g(0) = 0 and the Green build accepts it.
static void PlateLoading2D(const Vector &x, real_t t, Vector &u)
{
   const real_t n[2] = {std::sqrt(3.0) / 2.0, 0.5}, sd[2] = {-0.5, std::sqrt(3.0) / 2.0};
   const real_t Vp = 1e-9, d = x(0) * n[0] + x(1) * n[1];
   const real_t Vh = (d > 1.0) ? 0.5 * Vp * t : (d < -1.0) ? -0.5 * Vp * t : Vp * t;
   u.SetSize(2);
   u(0) = Vh * sd[0];
   u(1) = Vh * sd[1];
}

using Op = ElasticityDomainOperator<ParMesh>;

static std::unique_ptr<Op> MakeOp(ParMesh &pm, int order)
{
   static const real_t ref[3] = {1, 0, 0}, up[3] = {0, 1, 0};
   LinearElastic *model = new LinearElastic(kLambda, kMu);   // leaked: test lifetime
   DomainConfig dc;
   FaultFrameSource src;
   src.config = FrameConfig(ref, up);
   dc.fault_frames = src;
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = PlateLoading2D;
   return std::make_unique<Op>(pm, order, *model, 0.0, kWf, 2.0 * kLs, bc,
                               DGMethod::IP, SolverType::MUMPS, dc);
}

// One production stage at t = 0 (no loading): the Green stage, then the
// open-set hook, as SEASQuasiDynamicOperator::Mult runs them.
static void StageAt0(Op &op, cycle::OpenSetQD &os, Vector &slip, Vector &T,
                     Vector &Tn)
{
   op.GreenStage(0.0, slip, T, &Tn);
   os.Stage(slip, T, &Tn);
}

// Stage + rule until nothing flips (QDPhase::SettleOpenSet_'s loop, without
// the TS).  Returns the number of passes with flips (-1: did not settle).
static int Settle(Op &op, cycle::OpenSetQD &os, const Vector &sigma_pre,
                  Vector &slip, Vector &T, Vector &Tn)
{
   const int nown = sigma_pre.Size();
   Vector sig(nown);
   for (int pass = 0; pass < 25; pass++)
   {
      StageAt0(op, os, slip, T, Tn);
      for (int j = 0; j < nown; j++) { sig(j) = sigma_pre(j) + Tn(j); }
      if (os.Update(sig) == 0) { return pass; }
   }
   return -1;
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   MPI_Comm comm = MPI_COMM_WORLD;
   g_rank = Mpi::WorldRank();

   // ---- Q7: the tapered floor law (host, no mesh) ----
   {
      const real_t fl = 10e6;
      real_t worst = 0.0;
      bool status_ok = true;
      for (real_t eff : {-5e6, -1e3, 0.0, 1e3, 5e6, 9.999e6, 10e6, 10.001e6, 15e6})
      {
         int s1 = 0, s2 = 0;
         const real_t a = RSFlooredSigmaN(eff, fl, s1);
         const real_t b = RSLawSigmaN(eff, eff, fl, false, s2);   // r = 0
         worst = std::max(worst, std::abs(a - b));
         status_ok = status_ok && (s1 == s2);
      }
      TEST_ASSERT(worst == 0.0 && status_ok, "Q7: taper off = the hard floor, bit for "
                  "bit, value and status (worst " << worst << ")");
      int st = 0;
      const real_t r = 4e6;   // a relaxed point: eff = sigma - r
      const real_t above = RSLawSigmaN(15e6 - r, 15e6, fl, true, st);
      const real_t band = RSLawSigmaN(5e6 - r, 5e6, fl, true, st);
      const real_t at0 = RSLawSigmaN(0.0 - r, 0.0, fl, true, st);
      const real_t below = RSLawSigmaN(-1e3 - r, -1e3, fl, true, st);
      TEST_ASSERT(above == 15e6 - r && band == 5e6 && at0 == 0.0 && below == 0.0,
                  "Q7: taper on — sigma - r above the floor, sigma in the band, 0 at "
                  "and below opening (" << above << ", " << band << ", " << at0
                  << ", " << below << ")");
      const real_t eps = 1.0;   // 1 Pa against 10 MPa
      const real_t e1 = std::abs(RSLawSigmaN(fl + eps - r, fl + eps, fl, true, st) -
                                 RSLawSigmaN(fl - eps - r, fl - eps, fl, true, st));
      const real_t e2 = std::abs(RSLawSigmaN(eps - r, eps, fl, true, st) -
                                 RSLawSigmaN(-eps - r, -eps, fl, true, st));
      TEST_ASSERT(e1 <= 2.0 * eps && e2 <= 2.0 * eps, "Q7: the tapered strength is "
                  "continuous at sigma = floor and at sigma = 0 (" << e1 << ", "
                  << e2 << " Pa per " << 2.0 * eps << ")");
      int sn = 0;
      const real_t nan_through = RSLawSigmaN(std::nan(""), std::nan(""), fl, true, sn);
      TEST_ASSERT(std::isnan(nan_through), "Q7: a NaN stays NaN for the input checks");
      // the whole-law form against eq. 6.2 as B3 wrote it (the reference below
      // is that helper's taper branch, line for line): equal bit for bit
      // (value and status) for r >= 0 at every sigma, and for any r at
      // sigma >= floor
      auto eq62 = [&](real_t eff, real_t sig, int &stt) -> real_t
      {
         const real_t arm = std::fmin(fl, std::fmax(sig, (real_t) 0.0));
         if (eff < arm) { stt |= kRSSigmaNFloored; return arm; }
         return eff;
      };
      long n62 = 0, n62_diff = 0;
      for (real_t r0 : {0.0, 1e3, 4.0e6, 12.0e6, 40.0e6, -1e3, -4.0e6, -30.0e6})
      {
         for (real_t sig : {-5e6, -1e3, 0.0, 1e3, 5e6, 9.999e6, 10e6, 10.001e6, 15e6, 60e6})
         {
            if (r0 < 0.0 && sig < fl) { continue; }   // the whole-law rows, below
            int s1 = 0, s2 = 0;
            const real_t a = eq62(sig - r0, sig, s1);
            const real_t b = RSLawSigmaN(sig - r0, sig, fl, true, s2);
            n62++;
            if (std::memcmp(&a, &b, sizeof(real_t)) != 0 || s1 != s2) { n62_diff++; }
         }
      }
      TEST_ASSERT(n62_diff == 0, "Q7 whole law: = eq. 6.2 bit for bit (value and status) "
                  "for r >= 0 at every sigma and for r < 0 at sigma >= floor (" << n62
                  << " points, " << n62_diff << " differ)");
      // r < 0 (unclamping: the relaxation holds sigma - r up), below the floor:
      // max((sigma - r) w, sigma), w = sigma / floor; 0 at and below opening
      const real_t rn = -30.0e6;
      int swn = 0;
      const real_t wb = RSLawSigmaN(5e6 - rn, 5e6, fl, true, swn);   // 35 MPa x 0.5
      int s62 = 0;
      const real_t wb62 = eq62(5e6 - rn, 5e6, s62);                  // 35 MPa
      const real_t wl = RSLawSigmaN(1e5 - rn, 1e5, fl, true, swn);   // 30.1 MPa x 0.01
      const real_t wlo = RSLawSigmaN(5e6 + 2e6, 5e6, fl, true, swn); // r = -2: 7 x 0.5 < 5
      const real_t w0 = RSLawSigmaN(0.0 - rn, 0.0, fl, true, swn);
      const real_t wneg = RSLawSigmaN(-1e3 - rn, -1e3, fl, true, swn);
      TEST_ASSERT(wb == 17.5e6 && wb62 == 35.0e6 && std::abs(wl - 0.301e6) <= 1e-9 * 0.301e6 &&
                  wlo == 5e6 && w0 == 0.0 && wneg == 0.0 && (swn & kRSSigmaNFloored) != 0,
                  "Q7 whole law: r < 0 below the floor: (sigma - r) w, never below sigma ("
                  << wb / 1e6 << " MPa at sigma = 5 MPa, r = -30 MPa, where eq. 6.2 keeps "
                  << wb62 / 1e6 << " MPa), 0 at and below opening");
      const real_t f1 = std::abs(RSLawSigmaN(fl + eps - rn, fl + eps, fl, true, st) -
                                 RSLawSigmaN(fl - eps - rn, fl - eps, fl, true, st));
      const real_t f2 = std::abs(RSLawSigmaN(eps - rn, eps, fl, true, st) -
                                 RSLawSigmaN(-eps - rn, -eps, fl, true, st));
      const real_t lip = 2.0 + std::abs(rn) / fl;   // the slope of (sigma - r) sigma / floor
      TEST_ASSERT(f1 <= 2.0 * eps * lip && f2 <= 2.0 * eps * lip, "Q7 whole law: continuous "
                  "at sigma = floor and at sigma = 0 with r = -30 MPa (" << f1 << ", " << f2
                  << " Pa per 2 Pa; slope <= 2 + |r| / floor = 5)");
      int sn_st2 = 0;
      TEST_ASSERT(std::isnan(RSLawSigmaN(std::nan(""), 5e6, fl, true, sn_st2)),
                  "Q7 whole law: a NaN sigma - r stays NaN in the band");
   }

   // ---- Q9: the bridge with the gap (host, synthetic frames) ----
   {
      const int N = 5;
      Vector packed(cycle::FrameMap::kPerDOF * N);
      for (int j = 0; j < N; j++)
      {
         const real_t th = 0.37 * j - 0.6, c = std::cos(th), s = std::sin(th);
         real_t *f = packed.GetData() + cycle::FrameMap::kPerDOF * j;
         f[0] = (j % 2 == 0) ? 1.0 : -1.0;                   // eps_f
         f[1] = c; f[2] = s; f[3] = 0.0;                     // n_k
         f[4] = 0.0; f[5] = 0.0; f[6] = 1.0;                 // t1_k
         f[7] = s; f[8] = -c; f[9] = 0.0;                    // t2_k = n x t1
      }
      cycle::FrameMap f(packed);
      Vector sk(2 * N), vk(2 * N), psi_qd(N), gap(N);
      for (int j = 0; j < N; j++)
      {
         sk(2 * j) = 0.5 * j - 1.0;
         sk(2 * j + 1) = 0.25 * j + 0.1;
         vk(2 * j) = 1e-9 * (j + 1);
         vk(2 * j + 1) = -2e-9 * j;
         psi_qd(j) = 0.7;
         gap(j) = (j % 2 == 1) ? 0.01 * (j + 1) : 0.0;      // 0 at "closed" DOFs
      }
      Vector d, V, psi, sk2, psi2, gap2;
      cycle::bridge::ToNitsche(f, sk, vk, psi_qd, d, V, psi, &gap);
      real_t e_norm = 0.0, e_vnorm = 0.0;
      for (int j = 0; j < N; j++)
      {
         const real_t dn = f.Eps(j) * fault_frame::Dot(&d(3 * j), f.Nk(j));
         const real_t vn = f.Eps(j) * fault_frame::Dot(&V(3 * j), f.Nk(j));
         e_norm = std::max(e_norm, std::abs(dn - gap(j)));
         e_vnorm = std::max(e_vnorm, std::abs(vn));
      }
      TEST_ASSERT(e_norm <= 1e-14, "Q9: d . n = g after ToNitsche with the gap ("
                  << e_norm << " m)");
      // the two tangential products round independently, so the normal dot
      // is zero only to an ulp of |V| ~ 1e-9 m/s
      TEST_ASSERT(e_vnorm <= 1e-23, "Q9: V keeps no normal part (" << e_vnorm
                  << " m/s against 1e-23)");
      cycle::bridge::ToQD(f, d, psi, sk2, psi2, &gap2);
      real_t e_slip = 0.0, e_gap = 0.0;
      for (int j = 0; j < N; j++)
      {
         e_slip = std::max(e_slip, std::abs(sk2(2 * j) - sk(2 * j)));
         e_slip = std::max(e_slip, std::abs(sk2(2 * j + 1) - sk(2 * j + 1)));
         e_gap = std::max(e_gap, std::abs(gap2(j) - gap(j)));
      }
      TEST_ASSERT(e_slip <= 1e-14, "Q9: the tangential pairs cross the bridge "
                  "unchanged (" << e_slip << " m)");
      TEST_ASSERT(e_gap <= 1e-14, "Q9: GapFromGlobal reads the gap back ("
                  << e_gap << " m)");
      Vector d0, V0, psi0, d_ref;
      cycle::bridge::ToNitsche(f, sk, vk, psi_qd, d0, V0, psi0);   // no gap
      f.ToGlobal(sk, d_ref);
      bool same = d0.Size() == d_ref.Size();
      for (int i = 0; same && i < d0.Size(); i++) { same = d0(i) == d_ref(i); }
      TEST_ASSERT(same, "Q9: without a gap d is the tangential map alone, bit "
                  "for bit");
      real_t e_zero = 0.0;
      for (int j = 0; j < N; j++)
      {
         e_zero = std::max(e_zero, std::abs(fault_frame::Dot(&d0(3 * j), f.Nk(j))));
      }
      // d = a t1 + b t2 rounds its products on their own, and a compiler that
      // contracts the dot into an FMA (Apple clang on arm64 does by default)
      // rounds one of them less, so d . n is zero only to an ulp of |d| ~ 1 m
      TEST_ASSERT(e_zero <= 1e-15, "Q9: without a gap d keeps no normal part ("
                  << e_zero << " m against 1e-15)");
   }

   // ---- Q10: checkpoint FD_CYCLE_V7 round trip ----
   {
      MPIContext mpi(MPI_COMM_WORLD);
      const int n = 7 + g_rank;
      cycle::CheckpointInfo w;
      w.t = 1.25e9;
      w.dt_next = 3.5;
      w.step = 42;
      w.commit = "contact-test";
      w.ic = "scec_fd";
      w.fault_digest = 1.0;
      w.opening = true;
      w.open.resize(n);
      w.contact = true;
      w.gap.SetSize(n);
      for (int i = 0; i < n; i++)
      {
         w.open[i] = (i % 3 == 1) ? 1 : 0;
         w.gap(i) = w.open[i] ? 1.0e-3 * (i + 1) + 1.0e-13 : 0.0;
      }
      cycle::CheckpointData d;
      d.state.SetSize(3 * n);
      d.state = 0.5;
      d.displacement.SetSize(4);
      d.displacement = 0.0;
      d.traction.SetSize(2 * n);
      d.traction = 1.0;
      d.slip_rate.SetSize(2 * n);
      d.slip_rate = 1e-9;
      const std::string prefix = "test_fault_contact_qd_ckpt";
      cycle::WriteCycleCheckpoint(prefix, w, d, mpi);
      cycle::CheckpointInfo r;
      cycle::CheckpointData e;
      cycle::ReadCycleCheckpoint(prefix, "contact-test", false, r, e, mpi);
      bool ok = r.opening && r.contact && r.open == w.open && r.gap.Size() == n;
      for (int i = 0; ok && i < n; i++) { ok = (r.gap(i) == w.gap(i)); }
      int all = ok ? 1 : 0;
      MPI_Allreduce(MPI_IN_PLACE, &all, 1, MPI_INT, MPI_MIN, comm);
      TEST_ASSERT(all == 1, "Q10: V7 flags and gap read back bit for bit");
      std::remove(CheckpointFilename(prefix, mpi.Rank()).c_str());

      w.contact = false;        // shear-free stays V6
      w.gap.SetSize(0);
      const std::string p6 = "test_fault_contact_qd_ckpt6";
      cycle::WriteCycleCheckpoint(p6, w, d, mpi);
      cycle::CheckpointInfo r6;
      cycle::CheckpointData e6;
      cycle::ReadCycleCheckpoint(p6, "contact-test", false, r6, e6, mpi);
      int v6 = (r6.opening && !r6.contact && r6.gap.Size() == 0 && r6.open == w.open)
               ? 1 : 0;
      MPI_Allreduce(MPI_IN_PLACE, &v6, 1, MPI_INT, MPI_MIN, comm);
      TEST_ASSERT(v6 == 1, "Q10: without contact the block stays V6 (no gap)");
      std::remove(CheckpointFilename(p6, mpi.Rank()).c_str());
   }

   const real_t sigma_tol = 1000.0;   // the D-B7 default
   for (int p : {1, 2})
   {
      if (g_rank == 0) { std::cout << "\n=== B3: QD contact, p = " << p << " ===\n"; }
      std::unique_ptr<Mesh> ser(ShearBox(60.0));
      ParMesh pm(comm, *ser);
      auto op = MakeOp(pm, p);
      op->BuildQDStageOperators(true, 3);
      op->BuildQDGreen(QDGreenConfig{});
      Op *dom = op.get();
      const int nown = op->GetNumOwnedFaultDOFs();
      const int slot = op->FrameField()->InPlaneSlot();
      long long first = 0, mine = nown;
      MPI_Exscan(&mine, &first, 1, MPI_LONG_LONG, MPI_SUM, comm);
      if (g_rank == 0) { first = 0; }
      long long nf = first + nown;
      MPI_Allreduce(MPI_IN_PLACE, &nf, 1, MPI_LONG_LONG, MPI_MAX, comm);

      // The pre-stress the test owns (held by reference in OpenSetQD):
      // modest shear, strong compression, a tensile patch of 3 global DOFs
      // in the middle.
      Vector tau_pre(2 * nown), sigma_pre(nown);
      tau_pre = 0.0;
      const long long patch0 = nf / 2 - 1, patch1 = nf / 2 + 1;
      auto in_patch = [&](int j) {
         const long long gj = first + j;
         return gj >= patch0 && gj <= patch1;
      };
      for (int j = 0; j < nown; j++)
      {
         tau_pre(2 * j + slot) = 1e6 * std::sin(0.61 * (first + j) + 0.3);
         sigma_pre(j) = in_patch(j) ? -2e6 : 30e6;
      }

      cycle::OpenSetQD os(
         comm, nown, slot,
         [dom](const Vector &e, Vector &T, Vector &Tn) { dom->GreenStage(0.0, e, T, &Tn); },
         tau_pre,
         [dom](const Vector &e, Vector &T, Vector &Tn) { dom->QDUnitColumn(e, T, &Tn); },
         sigma_pre, sigma_tol);

      // ---- Q4: empty set = untouched ----
      {
         Vector slip(2 * nown), T0, Tn0, T1, Tn1;
         slip = 0.0;
         op->GreenStage(0.0, slip, T0, &Tn0);
         T1 = T0;
         Tn1 = Tn0;
         os.Stage(slip, T1, &Tn1);
         real_t d = 0.0;
         for (int i = 0; i < T0.Size(); i++) { d = std::max(d, std::abs(T0(i) - T1(i))); }
         for (int i = 0; i < Tn0.Size(); i++) { d = std::max(d, std::abs(Tn0(i) - Tn1(i))); }
         d = RMax(d, comm);
         TEST_ASSERT(d == 0.0 && os.NumOpenGlobal() == 0, "Q4 p = " << p << ": an empty "
                     "open set leaves the stage bit for bit (" << d << ")");
      }

      // ---- Q2 (and Q1 inside it): settle from all closed ----
      Vector slip(2 * nown), T, Tn;
      slip = 0.0;
      const int passes = Settle(*op, os, sigma_pre, slip, T, Tn);
      TEST_ASSERT(passes >= 1, "Q2 p = " << p << ": the set settles ("
                  << passes << " passes with flips)");
      {
         // one more stage at the settled set: residuals and idempotence
         Vector slip1(slip), T1, Tn1;
         StageAt0(*op, os, slip1, T1, Tn1);
         real_t res_t = 0.0, res_n = 0.0, scale = 0.0, gneg = 0.0;
         real_t gmin = std::numeric_limits<real_t>::infinity();   // over OPEN DOFs
         real_t comp = 0.0, sig_closed_min = 0.0, dslip = 0.0;
         long n_open_patch = 0;
         const std::vector<real_t> &gap = os.Gap();
         for (int j = 0; j < nown; j++)
         {
            scale = std::max(scale, std::abs(sigma_pre(j)));
            const real_t sig = sigma_pre(j) + Tn1(j);
            if (os.IsOpen(j))
            {
               res_t = std::max(res_t, std::abs(tau_pre(2 * j + slot) + T1(2 * j + slot)));
               res_n = std::max(res_n, std::abs(sig));
               gmin = std::min(gmin, gap[j]);
               if (in_patch(j)) { n_open_patch++; }
            }
            else
            {
               sig_closed_min = std::min(sig_closed_min, sig);
               gneg = std::max(gneg, std::abs(gap[j]));
            }
            comp = std::max(comp, std::abs(gap[j] * sig));
            dslip = std::max(dslip, std::abs(slip1(2 * j + slot) - slip(2 * j + slot)));
         }
         res_t = RMax(res_t, comm);
         res_n = RMax(res_n, comm);
         scale = RMax(scale, comm);
         gmin = -RMax(-gmin, comm);
         gneg = RMax(gneg, comm);
         comp = RMax(comp, comm);
         sig_closed_min = -RMax(-sig_closed_min, comm);
         dslip = RMax(dslip, comm);
         n_open_patch = RSuml(n_open_patch, comm);
         if (g_rank == 0)
         {
            std::cout << "  Q2 p = " << p << ": " << os.NumOpenGlobal() << " open ("
                      << n_open_patch << " of 3 patch DOFs), residual tau " << res_t
                      << " Pa, sigma " << res_n << " Pa, min g "
                      << gmin << " m, min closed sigma " << sig_closed_min / 1e6
                      << " MPa, max |g sigma| " << comp << "\n";
         }
         TEST_ASSERT(res_t <= 1e-8 * scale && res_n <= 1e-8 * scale,
                     "Q1 p = " << p << ": tau and sigma vanish at the open DOFs ("
                     << res_t << ", " << res_n << " Pa per " << scale << ")");
         TEST_ASSERT(n_open_patch == 3 && gmin > 0.0,
                     "Q1 p = " << p << ": the tensile patch is open with g > 0 (min g "
                     << gmin << " m)");
         TEST_ASSERT(gneg == 0.0 && sig_closed_min >= -sigma_tol - 1e-8 * scale &&
                     comp <= 1e-8 * scale * 1.0 /* g ~ mm, sigma ~ Pa residual */,
                     "Q2 p = " << p << ": g = 0 at closed DOFs, closed sigma >= "
                     "-sigma_tol, g sigma = 0 (" << gneg << ", "
                     << sig_closed_min << ", " << comp << ")");
         TEST_ASSERT(dslip <= 1e-12 * (1.0 + slip.Normlinf()),
                     "Q1 p = " << p << ": a second stage at the solved state is "
                     "idempotent (max slip change " << dslip << " m)");
      }

      // ---- Q3: hysteresis band, then a strong reversal ----
      {
         // inside the band at one closed DOF: no flip
         Vector sig(nown);
         Vector slip3(slip), T3, Tn3;
         StageAt0(*op, os, slip3, T3, Tn3);
         long flips = 0;
         {
            for (int j = 0; j < nown; j++)
            {
               sig(j) = sigma_pre(j) + Tn3(j);
               if (!os.IsOpen(j) && first + j == patch0 - 2)
               {
                  sig(j) = -0.5 * sigma_tol;   // tensile, inside the band
               }
            }
            flips = os.Update(sig);
         }
         TEST_ASSERT(flips == 0, "Q3 p = " << p << ": a grazing point inside the "
                     "band (sigma = -sigma_tol/2) does not flip");
         // beyond the band: it opens (drive it with the pre-stress so the
         // settled state is consistent)
         const long n_before = os.NumOpenGlobal();
         for (int j = 0; j < nown; j++)
         {
            if (first + j == patch0 - 2) { sigma_pre(j) = -5.0 * sigma_tol; }
         }
         const int p2 = Settle(*op, os, sigma_pre, slip, T, Tn);
         TEST_ASSERT(p2 >= 1 && os.NumOpenGlobal() == n_before + 1,
                     "Q3 p = " << p << ": beyond the band (sigma = -5 sigma_tol) the "
                     "point opens (" << n_before << " -> " << os.NumOpenGlobal() << ")");
         // a strong reversal: the whole patch recompresses, g wants < 0,
         // every open DOF closes and its gap is projected to 0
         for (int j = 0; j < nown; j++) { sigma_pre(j) = 30e6; }
         const int p3 = Settle(*op, os, sigma_pre, slip, T, Tn);
         real_t gleft = 0.0;
         for (int j = 0; j < nown; j++) { gleft = std::max(gleft, std::abs(os.Gap()[j])); }
         gleft = RMax(gleft, comm);
         TEST_ASSERT(p3 >= 1 && os.NumOpenGlobal() == 0 && gleft == 0.0,
                     "Q3 p = " << p << ": a strong reversal closes every DOF and "
                     "projects g to 0 (" << os.NumOpenGlobal() << " open, max |g| "
                     << gleft << ")");
         const int p4 = Settle(*op, os, sigma_pre, slip, T, Tn);
         TEST_ASSERT(p4 == 0, "Q3 p = " << p << ": the closed state is stable "
                     "(no flips on a repeat settle)");
      }

      // ---- Q5: SetGap round trip and the carried gap ----
      {
         std::vector<char> flags(nown, 0);
         std::vector<real_t> g(nown, 0.0);
         for (int j = 0; j < nown; j++)
         {
            if (first + j == patch0) { flags[j] = 1; g[j] = 1e-3; }
         }
         os.SetFlags(flags);
         os.SetGap(g);
         real_t gback = 0.0;
         for (int j = 0; j < nown; j++)
         {
            if (first + j == patch0) { gback = os.Gap()[j]; }
         }
         gback = RMax(gback, comm);
         TEST_ASSERT(gback == 1e-3, "Q5 p = " << p << ": SetGap reads back (g = "
                     << gback << " m)");
         // the handed-over gap enters the stage: with slip = 0 and the
         // tension relieved by the carried g, the solve still zeroes the
         // open DOF's traction
         Vector slip5(2 * nown), T5, Tn5;
         slip5 = 0.0;
         StageAt0(*op, os, slip5, T5, Tn5);
         real_t res = 0.0, scale = 0.0;
         for (int j = 0; j < nown; j++)
         {
            scale = std::max(scale, std::abs(sigma_pre(j)));
            if (!os.IsOpen(j)) { continue; }
            res = std::max(res, std::abs(tau_pre(2 * j + slot) + T5(2 * j + slot)));
            res = std::max(res, std::abs(sigma_pre(j) + Tn5(j)));
         }
         res = RMax(res, comm);
         scale = RMax(scale, comm);
         TEST_ASSERT(res <= 1e-8 * scale, "Q5 p = " << p << ": the stage solves on "
                     "top of the carried gap (residual " << res << " Pa)");
         // a SetFlags close zeroes the gap
         std::fill(flags.begin(), flags.end(), 0);
         os.SetFlags(flags);
         real_t gz = 0.0;
         for (int j = 0; j < nown; j++) { gz = std::max(gz, std::abs(os.Gap()[j])); }
         gz = RMax(gz, comm);
         TEST_ASSERT(gz == 0.0, "Q5 p = " << p << ": a SetFlags close zeroes the "
                     "gap (max |g| " << gz << ")");
      }

      // ---- Q6: QDUnitColumn linearity and dense-G consistency ----
      {
         Vector ei(3 * nown), ej(3 * nown), eij(3 * nown);
         ei = 0.0;
         ej = 0.0;
         const long long gi = nf / 2, gj = nf / 2 - 1;
         if (gi >= first && gi < first + nown)
         {
            ei(3 * static_cast<int>(gi - first) + 2) = 1.0;
         }
         if (gj >= first && gj < first + nown)
         {
            ej(3 * static_cast<int>(gj - first) + slot) = 1.0;
         }
         eij = ei;
         eij += ej;
         Vector Ti, Tni, Tj, Tnj, Tij, Tnij;
         op->QDUnitColumn(ei, Ti, &Tni);
         op->QDUnitColumn(ej, Tj, &Tnj);
         op->QDUnitColumn(eij, Tij, &Tnij);
         real_t dlin = 0.0, scale = 0.0;
         for (int i = 0; i < 2 * nown; i++)
         {
            dlin = std::max(dlin, std::abs(Tij(i) - Ti(i) - Tj(i)));
            scale = std::max(scale, std::abs(Tij(i)));
         }
         for (int i = 0; i < nown; i++)
         {
            dlin = std::max(dlin, std::abs(Tnij(i) - Tni(i) - Tnj(i)));
            scale = std::max(scale, std::abs(Tnij(i)));
         }
         dlin = RMax(dlin, comm);
         scale = RMax(scale, comm);
         TEST_ASSERT(dlin <= 1e-9 * scale, "Q6 p = " << p << ": QDUnitColumn is "
                     "linear (" << dlin << " per " << scale << " Pa)");
         // the tangential unit column equals the dense G's stage product
         Vector e2(2 * nown);
         e2 = 0.0;
         if (gj >= first && gj < first + nown)
         {
            e2(2 * static_cast<int>(gj - first) + slot) = 1.0;
         }
         Vector Tg, Tng;
         op->GreenStage(0.0, e2, Tg, &Tng);
         real_t dg = 0.0;
         for (int i = 0; i < 2 * nown; i++) { dg = std::max(dg, std::abs(Tj(i) - Tg(i))); }
         for (int i = 0; i < nown; i++) { dg = std::max(dg, std::abs(Tnj(i) - Tng(i))); }
         dg = RMax(dg, comm);
         TEST_ASSERT(dg <= 1e-6 * scale, "Q6 p = " << p << ": a tangential unit "
                     "column equals the dense G's (" << dg << " per " << scale
                     << " Pa)");
      }
   }

   // ---- Q8: a valid contact deck loads with the defaults ----
   {
      const char *base = "simple_cases/configs/D0_straight.toml";
      const std::string fn = "test_fault_contact_qd_deck.toml";
      if (g_rank == 0)
      {
         std::ifstream in(base);
         MFEM_VERIFY(in.good(), "cannot read " << base);
         std::stringstream ss;
         ss << in.rdbuf();
         std::string deck = ss.str();
         const std::string g = "[cycle.qd.green]\ncheck_every = 10000\n";
         const size_t at = deck.find(g);
         MFEM_VERIFY(at != std::string::npos, "D0 has no check_every block");
         deck.erase(at, g.size());
         const std::string cy = "[cycle]\n";
         const size_t cat = deck.find(cy);
         MFEM_VERIFY(cat != std::string::npos, "D0 has no [cycle] section");
         deck.insert(cat + cy.size(), "sigma_n_floor_pa = 10.0e6\n");
         deck += "\n[cycle.qd.green]\nblock = 8\n"
                 "\n[cycle.opening]\nmode = \"contact\"\n";
         std::ofstream out(fn);
         out << deck;
      }
      MPI_Barrier(comm);
      cycle::CycleConfig c = cycle::CycleConfig::Load(fn);
      MPI_Barrier(comm);
      if (g_rank == 0) { std::remove(fn.c_str()); }
      TEST_ASSERT(c.opening.Contact() && c.opening.On() &&
                  c.opening.sigma_tol == 1000.0 && c.opening.floor_taper &&
                  c.opening.max_passes == 20,
                  "Q8: a contact deck loads (mode contact, sigma_tol_pa 1000, "
                  "floor_taper on, max_passes 20)");
   }

   int failed = num_failed;
   MPI_Allreduce(MPI_IN_PLACE, &failed, 1, MPI_INT, MPI_MAX, comm);
   if (g_rank == 0)
   {
      std::cout << "\n" << (failed ? "FAILED" : "PASSED") << " (" << num_tests
                << " checks, " << failed << " failed)\n";
   }
   return failed ? 1 : 0;
}

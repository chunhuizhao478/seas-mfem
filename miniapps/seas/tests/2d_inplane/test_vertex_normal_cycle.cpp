// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/2d_inplane/test_vertex_normal_cycle.cpp — vertex-normal plan rev1,
// phase P3: the fully dynamic cycle's two regimes with [fault_geometry]
// vertex_normals = "averaged" on the simple cases' map-view decks
// (seas_fd_cycle_driver's QDPhase and DynPhase, built here without the
// driver; no time stepping).
//
//   C1 frames, D5f deck (simple_cases/configs/D5f_arc015L10_floor.toml: a
//      1.5 deg arc over 10 km, 500 m, p = 1, absorbing walls) and its
//      averaged form (D5f_arc015L10_avg_floor.toml, read with the option
//      on): both regimes build with the option on; DynPhase's max |eps_f n_k - n| (n the
//      frame field's FrameAt) and its ghost mismatch <= 1e-12; every owned
//      DOF's frame (NitscheWaveOperator::FaultFrames) is orthonormal to
//      1e-13 and eps_f times its FrameMap triad to 1e-12; at every smooth
//      vertex the DOFs of the two faces carry one normal, the vertex
//      table's n_v, to 1e-12 (p = 1: every DOF sits on a vertex).  Face
//      mode: there the two normals differ by 2 sin(turn / 2) to 1e-12.
//   C2 the vertex table (cycle::FaultVertexRows), D5f deck: one row per
//      fault vertex, 2 tips, no kink; every row equals a serial
//      recomputation from the mesh on the root (face normals oriented by
//      ref_normal, n_v their normalized sum, the turn angle; 1e-12), so the
//      table does not depend on the partition; the turn angles add up to
//      the angle between the two end faces (1e-9 deg), 1.5 deg.
//   C3 the QD balance of a tapered slip is the dynamic operator's static
//      state for that slip with the option on (the slip-held residual of
//      test_2d5_width.cpp U5): max |M^-1 (F - K u_s + G d_s)| / a_K <=
//      max(1e-8, 10 x face mode's).  Absorbing walls off, as U5.
//   C4 straight fault, D0 deck (simple_cases/configs/D0_straight.toml): on
//      vs off give the same FrameMap, FaultFrames and kernel basis bit for
//      bit, and the same QD balance u_s and slip-held residual: bit for bit
//      on 1 rank; on more, u_s to relative 1e-12 and the residual (itself
//      round-off, about 1e-10 a_K) to 1e-12 a_K, since there two builds
//      with the option off already differ in the last bits of u_s (about
//      1e-15 relative) from one run to the next; the vertex table: 2 tips,
//      every other vertex smooth with turn 0.
//
//   make -f tests/2d_inplane/2d_inplane.mk test-vertex-normal-cycle
//
// Scope (Chunhui Zhao, 2026-11-01).

#include "mfem.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../cycle/scaled_mesh.hpp"
#include "../../cycle/cycle_config.hpp"
#include "../../cycle/cycle_frames.hpp"
#include "../../cycle/fault_fields.hpp"
#include "../../cycle/prestress.hpp"
#include "../../cycle/qd_phase.hpp"
#include "../../cycle/dyn_phase.hpp"
#include "../../cycle/state_bridge.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const char *kArcDeck = "simple_cases/configs/D5f_arc015L10_floor.toml";
static const char *kArcAvgDeck = "simple_cases/configs/D5f_arc015L10_avg_floor.toml";
static const char *kStraightDeck = "simple_cases/configs/D0_straight.toml";

static real_t RMax(real_t v, MPI_Comm c)
{
   MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, c);
   return v;
}
static real_t MaxAbs(const Vector &x, MPI_Comm c)
{
   real_t m = 0.0;
   for (int i = 0; i < x.Size(); i++) { m = std::max(m, std::abs(x(i))); }
   return RMax(m, c);
}
/// Same size and the same bits on every rank.
static bool SameBits(const Vector &a, const Vector &b, MPI_Comm c)
{
   int bad = (a.Size() != b.Size()) ? 1 : 0;
   if (!bad && a.Size() > 0)
   {
      bad = std::memcmp(a.GetData(), b.GetData(), sizeof(real_t) * a.Size()) ? 1 : 0;
   }
   MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_MAX, c);
   return bad == 0;
}
/// The root's verdict on every rank.
static bool RootSays(bool ok, MPI_Comm c)
{
   int v = ok ? 1 : 0;
   MPI_Bcast(&v, 1, MPI_INT, 0, c);
   return v == 1;
}

static std::vector<std::vector<real_t>> GatherRows(const std::vector<std::vector<real_t>> &rows,
                                                   int width, MPI_Comm c)
{
   std::vector<double> loc;
   for (const auto &r : rows) { loc.insert(loc.end(), r.begin(), r.end()); }
   int nl = static_cast<int>(loc.size()), np = 0;
   MPI_Comm_size(c, &np);
   std::vector<int> cnt(np), dsp(np, 0);
   MPI_Allgather(&nl, 1, MPI_INT, cnt.data(), 1, MPI_INT, c);
   for (int r = 1; r < np; r++) { dsp[r] = dsp[r - 1] + cnt[r - 1]; }
   std::vector<double> all(dsp[np - 1] + cnt[np - 1]);
   MPI_Allgatherv(loc.data(), nl, MPI_DOUBLE, all.data(), cnt.data(), dsp.data(), MPI_DOUBLE, c);
   std::vector<std::vector<real_t>> out;
   for (size_t o = 0; o + width <= all.size(); o += width)
   {
      out.emplace_back(all.begin() + o, all.begin() + o + width);
   }
   return out;
}

using PosKey = std::pair<long long, long long>;
/// (x, y) rounded to 1 mm.
static PosKey Pos(real_t x, real_t y) { return {std::llround(x * 1e3), std::llround(y * 1e3)}; }

static cycle::CycleConfig Deck(const char *path, bool averaged, bool absorbing = true)
{
   cycle::CycleConfig c = cycle::CycleConfig::Load(path);
   c.dyn.averaged_vertex_normals = averaged;
   if (!absorbing) { c.absorbing.faces = "off"; }
   return c;
}

static std::unique_ptr<Mesh> SerialMesh(const cycle::CycleConfig &c)
{
   std::unique_ptr<Mesh> sm = cycle::LoadScaledMesh(c.shared.mesh.path, c.mesh_scale);
   c.CheckDimension(sm->Dimension());
   MFEM_VERIFY(sm->Dimension() == 2, c.shared.mesh.path << " is not a 2-D mesh");
   return sm;
}

// The QD regime of a deck (test_2d5_width.cpp U5).
struct Setup
{
   std::unique_ptr<cycle::QDPhase> qd;
   std::unique_ptr<cycle::FaultParams> par;
   std::unique_ptr<cycle::Prestress> pre;
};

static void Build(ParMesh &pmesh, MPIContext &mpi, const cycle::CycleConfig &cfg, Setup &S)
{
   S.qd = std::make_unique<cycle::QDPhase>(pmesh, mpi, cfg);
   const spatial::SpatialFrictionConfig &s = cfg.shared;
   const MaterialField material = MaterialField::MakeConstant(
      s.material_fallback.lambda, s.material_fallback.mu, s.material_fallback.rho);
   const spatial::PorePressureSpec no_pp;
   S.par = std::make_unique<cycle::FaultParams>(cycle::BuildFaultParams(
      S.qd->Table(), *s.rate_state, cfg.fields, material, pmesh, no_pp));
   cycle::PrestressConfig pc = cfg.prestress;
   pc.a_max = cycle::GlobalMaxA(*S.par, mpi.GetComm());
   S.pre = std::make_unique<cycle::Prestress>(cycle::BuildPrestress(*S.par, pc));
   S.qd->Setup(*S.par, *S.pre);
}

static std::unique_ptr<cycle::DynPhase> MakeDyn(ParMesh &pmesh, const cycle::CycleConfig &cfg,
                                                Setup &S)
{
   return std::make_unique<cycle::DynPhase>(pmesh, cfg, S.qd->Table(), *S.par, *S.pre, 1.0,
                                            nullptr, S.qd->Domain().FrameField());
}

struct Balance
{
   Vector u, sk, Vqd, psi, sn;
   real_t t = 0.0, Vmax = 0.0;
};

/// The QD balance at t_s of a strike slip S0 cos^2(pi x / (2 Lf)) (kernel
/// pairs (dip, strike); the map-view decks' slip slot) with psi0.
static Balance Solve(cycle::QDPhase &qd, const cycle::Prestress &pre, real_t ts, MPIContext &mpi)
{
   const Vector &X = qd.Table().coords;
   const int N = qd.Table().Size();
   real_t Lf = 0.0;
   for (int i = 0; i < N; i++) { Lf = std::max(Lf, std::abs(X(3 * i))); }
   Lf = mpi.GlobalMax(Lf);
   Vector s(2 * N), psi(N);
   for (int i = 0; i < N; i++)
   {
      const real_t x = X(3 * i);
      const real_t c = (std::abs(x) < Lf) ? std::cos(0.5 * M_PI * x / Lf) : 0.0;
      s(2 * i) = 0.0;
      s(2 * i + 1) = 0.02 * c * c;
      psi(i) = pre.psi0(i);
   }
   Balance b;
   qd.SetFaultState(s, psi);
   b.Vmax = qd.Solve(ts);
   qd.Slip(b.sk);
   qd.Psi(b.psi);
   qd.SigmaN(b.sn);
   b.Vqd = qd.SlipRate();
   b.u = qd.Displacement();
   b.t = ts;
   return b;
}

/// max |M^-1 (F(t_s) - K u_s + G d_s)| / a_K, a_K = max |M^-1 K u_s|
/// (test_2d5_width.cpp U5): the dynamic operator's residual at the QD
/// balance with the slip held.
static real_t SlipHeldResidual(cycle::DynPhase &dyn, const Balance &b, MPI_Comm c)
{
   Vector d, V, psi, k, zero(b.u.Size()), Ku, MKu;
   cycle::bridge::ToNitsche(dyn.Frames(), b.sk, b.Vqd, b.psi, d, V, psi);
   dyn.Enter(b.u, d, V, psi, b.sn, b.t);
   zero = 0.0;
   dyn.Op().AccelerationWithSlip(b.u, zero, d, k);
   dyn.Op().StiffnessMult(b.u, Ku);
   dyn.Op().Mass().Apply(Ku, MKu);
   return MaxAbs(k, c) / MaxAbs(MKu, c);
}

/// Per owned DOF: x, y, n_x, n_y, the frame's orthonormality error, max
/// |frame - eps_f FrameMap triad|.
static std::vector<std::vector<real_t>> FrameRows(cycle::DynPhase &dyn)
{
   const NitscheWaveOperator &W = dyn.Op();
   const Vector &F = W.FaultFrames(), &X = W.FaultPositions();
   const cycle::FrameMap &M = dyn.Frames();
   std::vector<std::vector<real_t>> rows;
   for (int j = 0; j < W.NumOwnedFaultDOFs(); j++)
   {
      const real_t *v[3] = {&F(9 * j), &F(9 * j + 3), &F(9 * j + 6)};
      const real_t *w[3] = {M.Nk(j), M.T1k(j), M.T2k(j)};
      real_t orth = 0.0, tri = 0.0;
      for (int a = 0; a < 3; a++)
      {
         for (int b = 0; b < 3; b++)
         {
            const real_t dot = v[a][0] * v[b][0] + v[a][1] * v[b][1] + v[a][2] * v[b][2];
            orth = std::max(orth, std::abs(dot - (a == b ? 1.0 : 0.0)));
         }
         for (int c = 0; c < 3; c++) { tri = std::max(tri, std::abs(M.Eps(j) * w[a][c] - v[a][c])); }
      }
      rows.push_back({X(3 * j), X(3 * j + 1), v[0][0], v[0][1], orth, tri});
   }
   return rows;
}

/// The vertex table (FaultVertexRows) on the root, by position.
static std::map<PosKey, std::vector<real_t>> VertexTable(const FaultFrameField &f, MPI_Comm c)
{
   Vector rows;
   cycle::FaultVertexRows(f, c, rows);
   const int nc = static_cast<int>(cycle::FaultVertexColumns().size());
   std::map<PosKey, std::vector<real_t>> t;
   for (int i = 0; i + nc <= rows.Size(); i += nc)
   {
      t[Pos(rows(i + 1), rows(i + 2))] = std::vector<real_t>(rows.GetData() + i,
                                                             rows.GetData() + i + nc);
   }
   return t;
}

// ------------------------------------------------------------------ C1
static void TestFrames(MPIContext &mpi)
{
   MPI_Comm comm = mpi.GetComm();
   if (g_rank == 0) { std::cout << "\n=== C1 frames: " << kArcDeck << " ===\n"; }
   const cycle::CycleConfig on = cycle::CycleConfig::Load(kArcAvgDeck), off = Deck(kArcDeck, false);
   TEST_ASSERT(on.dyn.averaged_vertex_normals && on.dyn.vertex_normal_max_turn_deg == 5.0
               && on.shared.mesh.path == off.shared.mesh.path, "C1 " << kArcAvgDeck
               << ": vertex_normals = \"averaged\", 5 deg, D5f's mesh");
   std::unique_ptr<Mesh> sm = SerialMesh(on);
   ParMesh pm(comm, *sm);
   sm.reset();
   Setup Son, Soff;
   Build(pm, mpi, on, Son);
   Build(pm, mpi, off, Soff);
   auto dyn_on = MakeDyn(pm, on, Son);
   auto dyn_off = MakeDyn(pm, off, Soff);
   TEST_ASSERT(dyn_on->NormalMismatch() <= 1e-12 && dyn_on->GhostFrameMismatch() <= 1e-12,
               "C1 averaged: DynPhase max |eps_f n_k - n| = " << std::scientific
               << std::setprecision(2) << dyn_on->NormalMismatch() << ", ghost mismatch "
               << dyn_on->GhostFrameMismatch() << " (<= 1e-12)");
   const auto table = VertexTable(*Son.qd->Domain().FrameField(), comm);
   const auto ron = GatherRows(FrameRows(*dyn_on), 6, comm);
   const auto roff = GatherRows(FrameRows(*dyn_off), 6, comm);
   bool ok_frame = true, ok_on = true, ok_off = true;
   real_t orth = 0.0, tri = 0.0, spread_on = 0.0, nv_err = 0.0, off_err = 0.0, off_max = 0.0;
   int smooth_pairs = 0, groups = 0;
   if (g_rank == 0)
   {
      for (const auto &r : ron)
      {
         orth = std::max(orth, r[4]);
         tri = std::max(tri, r[5]);
      }
      ok_frame = (orth <= 1e-13 && tri <= 1e-12);
      std::map<PosKey, std::vector<std::vector<real_t>>> gon, goff;
      for (const auto &r : ron) { gon[Pos(r[0], r[1])].push_back(r); }
      for (const auto &r : roff) { goff[Pos(r[0], r[1])].push_back(r); }
      groups = static_cast<int>(gon.size());
      ok_on = (gon.size() == table.size() && goff.size() == table.size());
      for (const auto &g : gon)
      {
         auto it = table.find(g.first);
         if (it == table.end()) { ok_on = false; continue; }
         const std::vector<real_t> &v = it->second;   // key x y kind turn nv_x nv_y
         const int kind = static_cast<int>(v[3]);
         if (kind != 1) { continue; }
         if (g.second.size() != 2 || goff[g.first].size() != 2) { ok_on = false; continue; }
         smooth_pairs++;
         const auto &a = g.second[0], &b = g.second[1];
         spread_on = std::max(spread_on, std::hypot(a[2] - b[2], a[3] - b[3]));
         for (const auto *r : {&a, &b})
         {
            nv_err = std::max(nv_err, std::hypot((*r)[2] - v[5], (*r)[3] - v[6]));
         }
         const auto &fa = goff[g.first][0], &fb = goff[g.first][1];
         const real_t chord = 2.0 * std::sin(0.5 * v[4] * M_PI / 180.0);
         const real_t got = std::hypot(fa[2] - fb[2], fa[3] - fb[3]);
         off_err = std::max(off_err, std::abs(got - chord));
         off_max = std::max(off_max, got);
      }
      ok_on = ok_on && smooth_pairs > 0 && spread_on <= 1e-12 && nv_err <= 1e-12;
      ok_off = (off_err <= 1e-12 && off_max > 1e-4);
      std::cout << std::scientific << std::setprecision(3) << "  " << groups << " DOF positions, "
                << smooth_pairs << " smooth vertices; averaged: frame orthonormality " << orth
                << ", |frame - eps_f triad| " << tri << ", spread at a vertex " << spread_on
                << ", |n - n_v| " << nv_err << "; face: largest spread " << off_max
                << ", |spread - 2 sin(turn/2)| " << off_err << "\n";
   }
   TEST_ASSERT(RootSays(ok_frame, comm), "C1 averaged: every owned DOF's frame orthonormal "
               "(1e-13) and eps_f times its FrameMap triad (1e-12)");
   TEST_ASSERT(RootSays(ok_on, comm), "C1 averaged: at every smooth vertex the DOFs of both "
               "faces carry one normal, the vertex table's n_v (1e-12)");
   TEST_ASSERT(RootSays(ok_off, comm), "C1 face mode: the two normals at a vertex differ by "
               "2 sin(turn / 2) of the table's turn angle (1e-12)");
}

// ------------------------------------------------------------------ C2
static void TestVertexTable(MPIContext &mpi)
{
   MPI_Comm comm = mpi.GetComm();
   if (g_rank == 0) { std::cout << "\n=== C2 the vertex table: " << kArcDeck << " ===\n"; }
   const cycle::CycleConfig on = Deck(kArcDeck, true);
   std::unique_ptr<Mesh> sm = SerialMesh(on);
   ParMesh pm(comm, *sm);
   cycle::QDPhase qd(pm, mpi, on);
   const auto table = VertexTable(*qd.Domain().FrameField(), comm);
   bool ok_rows = true, ok_sum = true;
   if (g_rank == 0)
   {
      // serial recomputation: the fault segments of the mesh, oriented by ref_normal
      const real_t *ref = on.dyn.ref_normal.data();
      std::map<int, std::vector<std::array<real_t, 2>>> at;   // vertex -> face normals
      std::vector<std::array<int, 2>> faces;
      Array<int> v;
      for (int be = 0; be < sm->GetNBE(); be++)
      {
         if (sm->GetBdrAttribute(be) != on.fault_attr) { continue; }
         sm->GetBdrElementVertices(be, v);
         const real_t *a = sm->GetVertex(v[0]), *b = sm->GetVertex(v[1]);
         real_t n[2] = {b[1] - a[1], -(b[0] - a[0])};
         const real_t l = std::hypot(n[0], n[1]);
         n[0] /= l;
         n[1] /= l;
         if (n[0] * ref[0] + n[1] * ref[1] < 0.0) { n[0] = -n[0]; n[1] = -n[1]; }
         at[v[0]].push_back({n[0], n[1]});
         at[v[1]].push_back({n[0], n[1]});
         faces.push_back({v[0], v[1]});
      }
      int tips = 0, kinks = 0, smooth = 0;
      real_t err_turn = 0.0, err_nv = 0.0, sum_turn = 0.0;
      std::vector<std::array<real_t, 2>> tip_n;
      for (const auto &kv : at)
      {
         const real_t *x = sm->GetVertex(kv.first);
         auto it = table.find(Pos(x[0], x[1]));
         if (it == table.end() || kv.second.size() > 2) { ok_rows = false; continue; }
         const std::vector<real_t> &r = it->second;
         if (kv.second.size() == 1)
         {
            tips++;
            tip_n.push_back(kv.second[0]);
            ok_rows = ok_rows && r[3] == 0.0 && r[4] == 0.0 && r[5] == 0.0 && r[6] == 0.0;
            continue;
         }
         const auto &A = kv.second[0], &B = kv.second[1];
         const real_t turn = std::atan2(std::abs(A[0] * B[1] - A[1] * B[0]),
                                        A[0] * B[0] + A[1] * B[1]) * 180.0 / M_PI;
         if (turn > on.dyn.vertex_normal_max_turn_deg) { kinks++; ok_rows = ok_rows && r[3] == 2.0; continue; }
         smooth++;
         real_t nv[2] = {A[0] + B[0], A[1] + B[1]};
         const real_t l = std::hypot(nv[0], nv[1]);
         nv[0] /= l;
         nv[1] /= l;
         ok_rows = ok_rows && r[3] == 1.0;
         err_turn = std::max(err_turn, std::abs(r[4] - turn));
         err_nv = std::max(err_nv, std::hypot(r[5] - nv[0], r[6] - nv[1]));
         sum_turn += r[4];
      }
      ok_rows = ok_rows && table.size() == at.size() && tips == 2 && kinks == 0
                && err_turn <= 1e-12 && err_nv <= 1e-12;
      real_t end_angle = -1.0;
      if (tip_n.size() == 2)
      {
         const auto &A = tip_n[0], &B = tip_n[1];
         end_angle = std::atan2(std::abs(A[0] * B[1] - A[1] * B[0]), A[0] * B[0] + A[1] * B[1])
                     * 180.0 / M_PI;
      }
      ok_sum = std::abs(sum_turn - end_angle) <= 1e-9 && std::abs(end_angle - 1.5) <= 1e-6;
      std::cout << std::setprecision(12) << "  " << table.size() << " rows (" << faces.size()
                << " fault faces): " << tips << " tips, " << smooth << " smooth, " << kinks
                << " kinks; |turn - serial| " << std::scientific << std::setprecision(2)
                << err_turn << " deg, |n_v - serial| " << err_nv << "; sum of turns "
                << std::fixed << std::setprecision(12) << sum_turn << " deg, end faces "
                << end_angle << " deg\n";
   }
   TEST_ASSERT(RootSays(ok_rows, comm), "C2 one row per fault vertex, 2 tips, no kink; kind, "
               "turn and n_v equal the serial recomputation (1e-12)");
   TEST_ASSERT(RootSays(ok_sum, comm), "C2 the turn angles add up to the angle between the end "
               "faces (1e-9 deg), the arc's 1.5 deg");
}

// ------------------------------------------------------------------ C3
static void TestBalance(MPIContext &mpi)
{
   MPI_Comm comm = mpi.GetComm();
   if (g_rank == 0) { std::cout << "\n=== C3 the QD balance, slip held: " << kArcDeck << " ===\n"; }
   const real_t ts = 3.0e7;
   real_t r[2] = {0.0, 0.0};
   for (int a = 0; a < 2; a++)
   {
      const cycle::CycleConfig cfg = Deck(kArcDeck, a == 1, false);
      std::unique_ptr<Mesh> sm = SerialMesh(cfg);
      ParMesh pm(comm, *sm);
      sm.reset();
      Setup S;
      Build(pm, mpi, cfg, S);
      const Balance b = Solve(*S.qd, *S.pre, ts, mpi);
      auto dyn = MakeDyn(pm, cfg, S);
      r[a] = SlipHeldResidual(*dyn, b, comm);
   }
   if (g_rank == 0)
   {
      std::cout << std::scientific << std::setprecision(3) << "  max |M^-1 (F - K u_s + G d_s)| = "
                << r[0] << " a_K (face), " << r[1] << " a_K (averaged)\n";
   }
   TEST_ASSERT(r[1] <= std::max(1e-8, 10.0 * r[0]), "C3 averaged: the QD balance is the "
               "dynamic operator's static state for its slip (residual " << std::scientific
               << std::setprecision(2) << r[1] << " a_K <= max(1e-8, 10 x face " << r[0] << ")");
}

// ------------------------------------------------------------------ C4
static void TestStraight(MPIContext &mpi)
{
   MPI_Comm comm = mpi.GetComm();
   if (g_rank == 0) { std::cout << "\n=== C4 straight fault, on vs off: " << kStraightDeck << " ===\n"; }
   const real_t ts = 3.0e7;
   Vector frames[2], ff[2], basis[2], u[2];
   real_t r[2] = {0.0, 0.0};
   bool ok_table = true;
   for (int a = 0; a < 2; a++)
   {
      const cycle::CycleConfig cfg = Deck(kStraightDeck, a == 1, false);
      std::unique_ptr<Mesh> sm = SerialMesh(cfg);
      ParMesh pm(comm, *sm);
      sm.reset();
      Setup S;
      Build(pm, mpi, cfg, S);
      const Balance b = Solve(*S.qd, *S.pre, ts, mpi);
      auto dyn = MakeDyn(pm, cfg, S);
      frames[a] = dyn->Frames().Packed();
      ff[a] = dyn->Op().FaultFrames();
      DenseMatrix B;
      dyn->Op().Elasticity().GetFaultDOFBasis(B);
      basis[a].SetSize(B.Height() * B.Width());
      for (int i = 0; i < basis[a].Size(); i++) { basis[a](i) = B.Data()[i]; }
      u[a] = b.u;
      r[a] = SlipHeldResidual(*dyn, b, comm);
      if (a == 1)
      {
         const auto table = VertexTable(*S.qd->Domain().FrameField(), comm);
         if (g_rank == 0)
         {
            int tips = 0, smooth = 0;
            for (const auto &kv : table)
            {
               const std::vector<real_t> &v = kv.second;
               if (v[3] == 0.0) { tips++; }
               else if (v[3] == 1.0 && v[4] == 0.0) { smooth++; }
            }
            ok_table = (tips == 2 && smooth + tips == static_cast<int>(table.size()));
            std::cout << "  vertex table: " << table.size() << " rows, " << tips << " tips, "
                      << smooth << " smooth with turn 0\n";
         }
      }
   }
   const bool same_map = SameBits(frames[0], frames[1], comm),
              same_ff = SameBits(ff[0], ff[1], comm), same_b = SameBits(basis[0], basis[1], comm);
   TEST_ASSERT(same_map && same_ff && same_b, "C4 straight: FrameMap (" << same_map
               << "), FaultFrames (" << same_ff << ") and the kernel basis (" << same_b
               << "), on == off bit for bit");
   if (mpi.Size() == 1)
   {
      TEST_ASSERT(SameBits(u[0], u[1], comm) && r[0] == r[1], "C4 straight: the QD balance u_s "
                  "and the slip-held residual (" << std::scientific << std::setprecision(2) << r[0]
                  << " a_K), on == off bit for bit (1 rank)");
   }
   else
   {
      Vector du(u[0]);
      du -= u[1];
      const real_t eu = MaxAbs(du, comm) / MaxAbs(u[0], comm);
      const real_t er = std::abs(r[1] - r[0]);
      TEST_ASSERT(eu <= 1e-12 && er <= 1e-12, "C4 straight: on == off, the QD balance u_s to "
                  "relative 1e-12 (" << std::scientific << std::setprecision(2) << eu << ") and "
                  "the slip-held residual to 1e-12 a_K (" << r[0] << " vs " << r[1] << "; "
                  << mpi.Size() << " ranks)");
   }
   TEST_ASSERT(RootSays(ok_table, comm), "C4 straight: the vertex table has 2 tips and every "
               "other vertex smooth with turn 0");
}

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);
   g_rank = mpi.Rank();
#ifndef MFEM_USE_PETSC
   if (mpi.IsRoot()) { std::cout << "test_vertex_normal_cycle needs PETSc (the QD regime)\n"; }
   return 1;
#else
   {
      const cycle::CycleConfig c = cycle::CycleConfig::Load(kArcDeck);
      MFEMInitializePetsc(&argc, &argv, c.qd.petsc_options.c_str(), NULL);
   }
   if (mpi.IsRoot())
   {
      std::cout << "test_vertex_normal_cycle (vertex-normal plan rev1, P3) on " << mpi.Size()
                << " rank(s)\n";
   }
   TestFrames(mpi);
   TestVertexTable(mpi);
   TestBalance(mpi);
   TestStraight(mpi);
   if (mpi.IsRoot())
   {
      std::cout << "\ntest_vertex_normal_cycle: " << (num_tests - num_failed) << " / " << num_tests
                << " passed" << (num_failed ? " -- FAILED" : "") << "\n";
   }
   MFEMFinalizePetsc();
   return num_failed ? 1 : 0;
#endif
}

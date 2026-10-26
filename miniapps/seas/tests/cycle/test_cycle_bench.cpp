// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_cycle_bench.cpp — fully dynamic BP5 cycle: the SCEC benchmark writer
// (cycle/bench_output.hpp) and the event updates it follows
// (FdCycleOutput::UpdateEvents), on made-up fields with known answers.  No
// solve; 1 or 2 ranks (with 2, each rank owns one of the two fault faces).
//   events   UpdateEvents returns Started, Below, Above, Ended with the
//            right t_end, and Ended with t_end = t when event_dwell = 0
//   source   max |V| over the DOFs in A and the moment rate over A
//            (face 1 lies outside A and moves faster)
//   catalog  an event row: R = the area weights of the DOFs that reached
//            rupture_V, tau_avg at the start and at the first time below,
//            the slip in between; the rupture file (1e9 = never)
//   restart  a writer built from GetState() in the middle of an event (the
//            files trimmed) ends it with the same row
//   profiles the first row and the cadence (dt_slow, dt_fast while
//            max V > fast_V)
//   off-fault a linear field on a small hex mesh, read at a station and
//            projected on rotated axes, both regimes' row layout
//   locate   LocateOnFaces and StationOwners: a point where several faces
//            meet reads the same face on every layout of the faces (order,
//            vertex order, ranks)
//
//   mpirun -np N ./seas_test_cycle_bench WORK     (N = 1 or 2; WORK replaced)

#include "mfem.hpp"
#include "../../cycle/bench_output.hpp"

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;
using namespace mfem::seas::cycle;

namespace
{

int rank_ = 0, nranks_ = 1, failures_ = 0;

void Check(const std::string &name, bool ok_local, const std::string &detail = "")
{
   int ok = ok_local ? 1 : 0, all = 0;
   MPI_Allreduce(&ok, &all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
   if (!all) { failures_++; }
   if (rank_ == 0)
   {
      std::cout << "  " << (all ? "ok  " : "FAIL") << " " << std::left << std::setw(40) << name
                << std::right << " " << detail << "\n";
   }
}

/// The data rows of a RowFile.
std::vector<std::vector<double>> Rows(const std::string &path)
{
   std::vector<std::vector<double>> out;
   std::ifstream f(path);
   std::string line;
   while (std::getline(f, line))
   {
      if (line.empty() || line[0] == '#') { continue; }
      std::istringstream ls(line);
      std::vector<double> r;
      double v;
      while (ls >> v) { r.push_back(v); }
      out.push_back(r);
   }
   return out;
}

bool Near(double a, double b, double tol = 1e-12)
{
   return std::abs(a - b) <= tol * std::max(1.0, std::max(std::abs(a), std::abs(b)));
}

// Two triangular faces of p = 1 (3 DOFs each): face 0 at x2 in [0, 1] km
// (inside A), face 1 at x2 in [50, 51] km (outside).  The 3-point edge
// midpoint rule: e_q(k, q) = 1/2 on the two ends of edge q, weight
// |face| / 3, so a_k = |face| / 3.
const double kArea = 0.5e6;   // m^2 per face
const double kMu = 3.0e10;

BenchFaultGeometry Geometry(const std::vector<int> &faces)
{
   BenchFaultGeometry g;
   g.nbf = 3;
   g.e_q.SetSize(3, 3);
   g.e_q = 0.0;
   const int ends[3][2] = {{0, 1}, {1, 2}, {2, 0}};
   for (int q = 0; q < 3; q++) { g.e_q(ends[q][0], q) = g.e_q(ends[q][1], q) = 0.5; }
   const int N = 3 * static_cast<int>(faces.size());
   g.X.SetSize(3 * N);
   g.x2.SetSize(N);
   g.x3.SetSize(N);
   g.wJ.SetSize(3 * static_cast<int>(faces.size()));
   g.q2.SetSize(g.wJ.Size());
   g.q3.SetSize(g.wJ.Size());
   for (size_t i = 0; i < faces.size(); i++)
   {
      const double x0 = faces[i] == 0 ? 0.0 : 50.0e3;
      const double p2[3] = {x0, x0 + 1000.0, x0}, p3[3] = {0.0, 0.0, 1000.0};
      for (int k = 0; k < 3; k++)
      {
         const int d = 3 * static_cast<int>(i) + k;
         g.x2(d) = p2[k];
         g.x3(d) = p3[k];
         g.X(3 * d) = p2[k];
         g.X(3 * d + 1) = 0.0;
         g.X(3 * d + 2) = -p3[k];
      }
      for (int q = 0; q < 3; q++)
      {
         const int j = 3 * static_cast<int>(i) + q;
         g.wJ(j) = kArea / 3.0;
         g.q2(j) = 0.5 * (p2[ends[q][0]] + p2[ends[q][1]]);
         g.q3(j) = 0.5 * (p3[ends[q][0]] + p3[ends[q][1]]);
      }
   }
   return g;
}

/// Fields on this rank's DOFs: V_strike per face (face 0: v0, face 1: v1),
/// slip_strike s, shear stress tau [MPa] (strike).
cycle::StationFields Fields(const std::vector<int> &faces, double v0, double v1, double s,
                            double tau)
{
   cycle::StationFields f;
   const int N = 3 * static_cast<int>(faces.size());
   for (Vector *v : {&f.slip_dip, &f.slip_strike, &f.V_dip, &f.V_strike, &f.tau_a_dip,
                     &f.tau_a_strike, &f.tau_b_dip, &f.tau_b_strike, &f.theta, &f.sigma_n})
   {
      v->SetSize(N);
      *v = 0.0;
   }
   for (size_t i = 0; i < faces.size(); i++)
   {
      for (int k = 0; k < 3; k++)
      {
         const int d = 3 * static_cast<int>(i) + k;
         f.V_strike(d) = faces[i] == 0 ? v0 : v1;
         f.slip_strike(d) = s;
         f.tau_b_strike(d) = -tau * 1e6;
      }
   }
   return f;
}

void EventCases(const std::string &work)
{
   cycle::OutputConfig oc;
   oc.dir = work;
   oc.event_V = 1e-3;
   {
      oc.prefix = "ev";
      oc.event_dwell = 2.5;
      cycle::FdCycleOutput out(oc, MPI_COMM_WORLD, {}, "test");
      using C = cycle::EventChange;
      const double V[] = {1e-4, 1e-2, 1e-2, 5e-4, 2e-3, 5e-4, 5e-4, 5e-4, 1e-4};
      const C want[] = {C::None, C::Started, C::None, C::Below, C::Above, C::Below, C::None,
                        C::Ended, C::None};
      bool ok = true;
      std::ostringstream d;
      for (int k = 0; k < 9; k++)
      {
         const cycle::EventUpdate u = out.UpdateEvents(k, V[k]);
         ok = ok && u.change == want[k];
         if (u.change == C::Ended) { ok = ok && u.t_end == 5.0 && u.t_start == 1.0 && u.event == 1; }
         d << static_cast<int>(u.change);
      }
      Check("UpdateEvents sequence", ok, d.str() + " (0 none 1 started 2 below 3 above 4 ended)");
   }
   {
      oc.prefix = "ev0";
      oc.event_dwell = 0.0;
      cycle::FdCycleOutput out(oc, MPI_COMM_WORLD, {}, "test");
      out.UpdateEvents(0.0, 1e-2);
      const cycle::EventUpdate u = out.UpdateEvents(0.5, 1e-4);
      Check("UpdateEvents dwell 0: Ended at t", u.change == cycle::EventChange::Ended &&
            u.t_end == 0.5);
   }
}

void WriterCases(const std::string &work)
{
   std::vector<int> faces;
   if (nranks_ == 1) { faces = {0, 1}; }
   else { faces = {rank_}; }
   const BenchFaultGeometry g0 = Geometry(faces);

   // A tiny hex mesh for the off-fault station: [0, 1]^3, 2 x 2 x 2.
   Mesh sm = Mesh::MakeCartesian3D(2, 2, 2, Element::HEXAHEDRON, 1.0, 1.0, 1.0);
   ParMesh pm(MPI_COMM_WORLD, sm);
   H1_FECollection fec(1, 3);
   ParFiniteElementSpace fes(&pm, &fec, 3);
   ParGridFunction u(&fes), v(&fes);
   VectorFunctionCoefficient cu(3, [](const Vector &x, Vector &y)
   { y(0) = 1 + 2 * x(0); y(1) = 3 * x(1) - x(2); y(2) = 4 * x(2) + x(0); });
   VectorFunctionCoefficient cv(3, [](const Vector &x, Vector &y)
   { y(0) = -x(1); y(1) = 0.5; y(2) = 7.0; });
   u.ProjectCoefficient(cu);
   v.ProjectCoefficient(cv);

   cycle::CycleConfig cfg;
   cfg.shared.material_fallback.mu = kMu;
   cfg.output.dir = work;
   cfg.output.prefix = "b";
   cfg.output.fast_V = 1e-3;
   cfg.output.event_V = 1e-3;
   cfg.output.event_dwell = 1.0;
   cycle::BenchmarkConfig &bc = cfg.bench;
   bc.enabled = true;
   bc.region_x2 = {{-34e3, 34e3}};
   bc.region_x3 = {{0.0, 20e3}};
   bc.rupture_V = 1e-3;
   bc.profile_dt_slow = 10.0;
   bc.profile_dt_fast = 1.0;
   bc.volume_axes = {{{{0.0, 1.0, 0.0}}, {{1.0, 0.0, 0.0}}, {{0.0, 0.0, -1.0}}}};
   cycle::StationSpec st;
   st.name = "blkst_test";
   st.x = {{0.3, 0.4, 0.5}};
   bc.volume_stations = {st};
   cycle::ProfileSpec ps;
   ps.name = "line";
   ps.from = {{0.0, 0.0}};
   ps.to = {{50e3, 0.0}};
   ps.points = 2;
   ps.coordinate = 0;
   bc.profiles = {ps};
   // Profile point 0 at face 0's first DOF, point 1 at face 1's.
   BenchFaultGeometry g = g0;
   g.profile_locs.resize(2);
   for (size_t i = 0; i < faces.size(); i++)
   {
      g.profile_locs[faces[i]].distance = 0.0;
      g.profile_locs[faces[i]].nearest = 3 * static_cast<int>(i);
   }

   auto fill_of = [&](double v0, double v1, double s, double tau)
   {
      return [&faces, v0, v1, s, tau](cycle::StationFields &f) { f = Fields(faces, v0, v1, s, tau); };
   };
   cycle::VolumeFieldsFn vol = [&](const Vector *&a, const Vector *&b) { a = &u; b = &v; };
   // |V| per DOF: face 0's DOF k at sp[k], face 1 at 0.
   auto speed_of = [&](double s0, double s1, double s2)
   {
      return [&faces, s0, s1, s2](Vector &out)
      {
         out.SetSize(3 * static_cast<int>(faces.size()));
         out = 0.0;
         for (size_t i = 0; i < faces.size(); i++)
         {
            if (faces[i] != 0) { continue; }
            out(3 * i) = s0;
            out(3 * i + 1) = s1;
            out(3 * i + 2) = s2;
         }
      };
   };

   cycle::BenchmarkOutput::State mid;
   {
      cycle::BenchmarkOutput b(cfg, MPI_COMM_WORLD, pm, fes, g, "test");
      Check("areas", Near(b.FaultArea(), 2 * kArea) && Near(b.AreaA(), kArea),
            "fault " + std::to_string(b.FaultArea()) + " m^2, A " + std::to_string(b.AreaA()));
      // Source and off-fault rows; profiles at t = 0 (first row), 5 (not due), 10.
      b.Rows(0.0, fill_of(2.0, 10.0, 0.0, 5.0), vol);
      b.Profiles(0.0, 1e-9, fill_of(2.0, 10.0, 0.0, 5.0));
      b.Profiles(5.0, 1e-9, fill_of(2.0, 10.0, 0.5, 5.0));
      b.Profiles(10.0, 1e-9, fill_of(2.0, 10.0, 1.0, 5.0));
      b.Profiles(10.5, 1e-2, fill_of(2.0, 10.0, 1.5, 5.0));   // fast: 0.5 s < 1 s, not due
      b.Profiles(11.0, 1e-2, fill_of(2.0, 10.0, 2.0, 5.0));   // due
      // An event: started at t = 20 (slip 1, tau 5 MPa); DOF 0 ruptures at
      // 20, DOF 1 at 21; first below at 22 (slip 3, tau 3 MPa); ends at 24.
      cycle::EventUpdate e;
      e.event = 1;
      e.t_start = 20.0;
      e.change = cycle::EventChange::Started;
      b.OnEvents(e, 20.0, fill_of(1e-2, 0.0, 1.0, 5.0), speed_of(1e-2, 0.0, 0.0));
      e.change = cycle::EventChange::None;
      b.Track(20.5, [&] { Vector s; speed_of(5e-4, 2e-3, 0.0)(s); return s; }());
      b.OnEvents(e, 21.0, fill_of(1e-2, 0.0, 2.0, 4.0), speed_of(5e-4, 2e-3, 0.0));
      mid = b.GetState();
      b.Flush();
   }
   // The rest from a writer built on the mid-event state (files trimmed).
   {
      cycle::BenchmarkOutput b(cfg, MPI_COMM_WORLD, pm, fes, g, "test", &mid);
      cycle::EventUpdate e;
      e.event = 1;
      e.t_start = 20.0;
      e.change = cycle::EventChange::Below;
      b.OnEvents(e, 22.0, fill_of(5e-4, 0.0, 3.0, 3.0), speed_of(5e-4, 0.0, 0.0));
      e.change = cycle::EventChange::Ended;
      e.t_end = 22.0;
      b.OnEvents(e, 24.0, fill_of(1e-4, 0.0, 3.5, 2.0), speed_of(1e-4, 0.0, 0.0));
      b.Flush();
   }
   MPI_Barrier(MPI_COMM_WORLD);

   if (rank_ == 0)
   {
      const auto src = Rows(work + "/b_source.txt");
      const bool ok = src.size() == 1 && src[0].size() == 3 && src[0][0] == 0.0 &&
                      Near(src[0][1], std::log10(2.0)) && Near(src[0][2], kMu * 2.0 * kArea);
      std::ostringstream d;
      if (src.size() == 1) { d << "log10 Vmax_A " << src[0][1] << ", moment rate " << src[0][2]; }
      Check("source row (face 1 outside A)", ok, d.str());
   }
   else { Check("source row (face 1 outside A)", true); }
   {
      // u.e1 = u_y, u.e2 = u_x, u.e3 = -u_z (whichever rank wrote the file).
      const std::string p = work + "/b_blkst_test.txt";
      bool ok = true;
      std::ostringstream d;
      if (rank_ == 0)
      {
         const auto r = Rows(p);
         const double x = 0.3, y = 0.4, z = 0.5;
         const double want[7] = {0.0, 3 * y - z, 1 + 2 * x, -(4 * z + x), 0.5, -y, -7.0};
         ok = r.size() == 1 && r[0].size() == 7;
         for (int c = 0; ok && c < 7; c++) { ok = Near(r[0][c], want[c], 1e-12); }
         d << (ok ? "u, v on (e1, e2, e3) exact" : "row differs");
      }
      Check("off-fault row (lowest rank writes)", ok, d.str());
   }
   if (rank_ == 0)
   {
      const auto pr = Rows(work + "/b_profile_line_slip_2.txt");
      // first row, t = 0, 10, 11
      bool ok = pr.size() == 4 && pr[0].size() == 4 && pr[0][0] == 0 && pr[0][1] == 0 &&
                pr[0][2] == 0.0 && pr[0][3] == 50e3;
      ok = ok && pr[1][0] == 0.0 && pr[2][0] == 10.0 && pr[3][0] == 11.0;
      ok = ok && Near(pr[2][2], 1.0) && Near(pr[2][3], 1.0) && Near(pr[3][1], -2.0);
      Check("profile rows and cadence", ok, std::to_string(pr.size()) + " rows");
      const auto ev = Rows(work + "/b_events.txt");
      const double a = kArea / 3.0, R = 2 * a;
      ok = ev.size() == 1 && ev[0].size() == 7;
      std::ostringstream d;
      if (ok)
      {
         d << "R " << ev[0][3] << ", tau " << ev[0][4] << " -> " << ev[0][5] << " MPa, slip "
           << ev[0][6] << " m";
         ok = ev[0][0] == 1 && ev[0][1] == 20.0 && ev[0][2] == 22.0 && Near(ev[0][3], R) &&
              Near(ev[0][4], 5.0) && Near(ev[0][5], 3.0) && Near(ev[0][6], 2.0);
      }
      Check("catalog row across the restart", ok, d.str());
      const auto rt = Rows(work + "/b_rupture_event1.txt");
      ok = rt.size() == 6;
      int ruptured = 0;
      for (const auto &r : rt)
      {
         ok = ok && r.size() == 6;
         if (!ok) { break; }
         if (r[5] != 1e9) { ruptured++; }
         if (r[3] == 0.0 && r[4] == 0.0) { ok = ok && r[5] == 20.0; }
         if (r[3] == 1000.0 && r[4] == 0.0) { ok = ok && r[5] == 20.5; }
      }
      Check("rupture times (1e9 = never)", ok && ruptured == 2,
            std::to_string(ruptured) + " of " + std::to_string(rt.size()) + " DOFs ruptured");
   }
   else
   {
      Check("profile rows and cadence", true);
      Check("catalog row across the restart", true);
      Check("rupture times (1e9 = never)", true);
   }
}

// Four p = 1 faces around the vertex O = (0, 0) km: T0 = (O, A, B),
// T1 = (O, B, C), T2 = (O, C, D), T3 = (O, D, A), with A = (2, 0),
// B = (0, 1), C = (-1, 0), D = (0, -1); vertex centres (2/3, 1/3),
// (-1/3, 1/3), (-1/3, -1/3), (2/3, -1/3).  The face each point must read on
// every layout (the nearest centre, then the smaller x2, then x3):
//   vertex   O, on all four faces: T1 and T2 tie on the centre's distance
//            and x2, T2 has the smaller x3                        -> T2
//   edge     (0, 0.25) on T0 and T1: T1's centre is nearer         -> T1
//   inside   (0.5, 0.2), on T0 only                                -> T0
//   outside  (-1.5, 0), on no face: the nearest DOF is C, on T1 and T2,
//            which tie as at O                                     -> T2
// Layouts: on 1 rank the faces in order, then reversed with each face's
// vertices rotated; on 2 ranks {T0, T1} | {T2, T3}, then {T2, T0} | {T3, T1}
// rotated, so T2 changes rank.  The fields read: 100 k on face k, and x2;
// the owner's key (StationOwners' gkey, the face centre of the run log) is
// the wanted face's centre.  The first-face rule (Probe2DInterpolator)
// reads T0 at O on the first 1-rank layout.
void LocateCases()
{
   const double km = 1000.0;
   const double V[5][2] = {{0, 0}, {2, 0}, {0, 1}, {-1, 0}, {0, -1}};   // O A B C D
   const int T[4][3] = {{0, 1, 2}, {0, 2, 3}, {0, 3, 4}, {0, 4, 1}};
   const std::vector<std::array<real_t, 2>> pts = {{{0.0, 0.0}}, {{0.0, 0.25 * km}},
                                                  {{0.5 * km, 0.2 * km}}, {{-1.5 * km, 0.0}}};
   const char *names[4] = {"vertex", "edge", "inside", "outside"};
   const int want[4] = {2, 1, 0, 2};
   const double want_x2[4] = {0.0, 0.0, 0.5 * km, -1.0 * km};
   std::vector<std::vector<std::vector<int>>> layouts;   // [layout][rank]: faces
   if (nranks_ == 1) { layouts = {{{0, 1, 2, 3}}, {{3, 2, 1, 0}}}; }
   else { layouts = {{{0, 1}, {2, 3}}, {{2, 0}, {3, 1}}}; }
   for (size_t l = 0; l < layouts.size(); l++)
   {
      const std::vector<int> &mine = layouts[l][rank_];
      const int N = 3 * static_cast<int>(mine.size());
      Vector x2(N), x3(N), id(N);
      for (size_t i = 0; i < mine.size(); i++)
      {
         for (int k = 0; k < 3; k++)
         {
            const int v = T[mine[i]][l == 1 ? (k + 1) % 3 : k];
            const int d = 3 * static_cast<int>(i) + k;
            x2(d) = V[v][0] * km;
            x3(d) = V[v][1] * km;
            id(d) = 100.0 * mine[i];
         }
      }
      const std::vector<StationLocation> locs = LocateOnFaces(x2, x3, 3, pts);
      std::vector<real_t> gmin;
      std::vector<std::array<real_t, 3>> gkey;
      const std::vector<int> owner = StationOwners(locs, MPI_COMM_WORLD, gmin, &gkey);
      for (size_t s = 0; s < pts.size(); s++)
      {
         double c2 = 0.0, c3 = 0.0;   // the wanted face's centre
         for (int k = 0; k < 3; k++)
         {
            c2 += V[T[want[s]][k]][0] * km / 3.0;
            c3 += V[T[want[s]][k]][1] * km / 3.0;
         }
         double got[4] = {0.0, 0.0, 0.0, 0.0}, sum[4];
         if (owner[s] == rank_)
         {
            got[0] = FdCycleOutput::Eval(id, locs[s]);
            got[1] = FdCycleOutput::Eval(x2, locs[s]);
            got[2] = 1.0;
            got[3] = locs[s].exact ? 1.0 : 0.0;
         }
         MPI_Allreduce(got, sum, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         const bool ok = sum[2] == 1.0 && Near(sum[0], 100.0 * want[s])
                         && Near(sum[1] / km, want_x2[s] / km) && sum[3] == (s < 3 ? 1.0 : 0.0)
                         && Near(gkey[s][1] / km, c2 / km) && Near(gkey[s][2] / km, c3 / km);
         std::ostringstream d;
         d << "layout " << l << ": face " << sum[0] / 100.0 << " (want " << want[s] << "), x2 "
           << sum[1] << " m, " << (sum[3] == 1.0 ? "on a face" : "nearest DOF") << ", "
           << sum[2] << " owner(s), centre (" << gkey[s][1] << ", " << gkey[s][2] << ") m";
         Check(std::string("location ") + names[s] + " " + std::to_string(l), ok, d.str());
      }
   }
}

} // namespace

int main(int argc, char **argv)
{
   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
   MPI_Comm_size(MPI_COMM_WORLD, &nranks_);
   if (argc < 2 || nranks_ > 2)
   {
      if (rank_ == 0) { std::cerr << "usage: mpirun -np 1|2 seas_test_cycle_bench WORK\n"; }
      MPI_Finalize();
      return 2;
   }
   const std::string work = argv[1];
   if (rank_ == 0)
   {
      std::filesystem::remove_all(work);
      std::filesystem::create_directories(work);
   }
   MPI_Barrier(MPI_COMM_WORLD);
   if (rank_ == 0) { std::cout << "test_cycle_bench on " << nranks_ << " rank(s)\n"; }
   EventCases(work);
   WriterCases(work);
   LocateCases();
   if (rank_ == 0)
   {
      std::cout << "test_cycle_bench: " << (failures_ ? "FAILED" : "PASSED") << " ("
                << failures_ << " failures)\n";
   }
   MPI_Finalize();
   return failures_ ? 1 : 0;
}

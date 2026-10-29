// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_snapshots_volume.cpp — the volume snapshots of [cycle.snapshots.volume]
// (cycle/snapshots.hpp; Chunhui Zhao, 2026-10-28: "we should add the volume
// vtkdhf option, which saves the velocity and displacement data in the 2D
// domain, this is particularly important for later UCERF 2D fault regional
// runs").  A unit test of cycle::Snapshots with only the volume on, on a DG
// GaussLobatto space as the QD regime's (vdim = dim, byNODES), and of the
// deck loader:
//
//   grid_2d_p1, grid_2d_p2, grid_3d_p1, grid_3d_p2
//                A 4 km x 3 km triangle box (4 x 3 squares, 2 triangles
//                each) and a 2 x 2 x 2 tetrahedral box, p = 1 and 2.  Seven
//                snapshots through QD and Dyn: QD at t = 0, 5, 12 s (dt_qd
//                10 s: 0 and 12 are written), the switch at 25 s, a dynamic
//                phase with tau = 0.05, 0.1, 0.15, 0.25 (dt_dyn 0.1 s: 0.1
//                and 0.25 written) and its end at tau = 0.3, then QD at 26
//                and 31 s (31 written).  Read back with HDF5:
//                  layout   VTKHDF 2.0 UnstructuredGrid; (p+1)(p+2)/2
//                           points (2-D) or (p+1)(p+2)(p+3)/6 (3-D) and p^d
//                           cells per element, all VTK_TRIANGLE (5) or
//                           VTK_TETRA (10), each cell's points inside its
//                           element's block, no degenerate cell; the step
//                           offsets; /fd_cycle per snapshot
//                  values   at every point of every snapshot, displacement
//                           and velocity = the projected polynomial fields
//                           (degree p: exact in the DG space) at the
//                           point's coordinates (z = 0 and zero z
//                           components in 2-D), to 1e-12 of their largest
//                           magnitude
//                  schedule the seven times, regimes and reasons
//                  restart  a second writer from the state after the 12 s
//                           snapshot, the same calls after it: the file
//                           equals the straight one dataset by dataset,
//                           bit for bit
//   box_2d, box_3d
//                box_min .. box_max: the elements whose centroid is in the
//                box, counted independently from the mesh vertices; every
//                element block of the file has its centroid in the box.
//
// Modes for the make target's refusals: --empty-box (a box that holds no
// element centroid), --quads (a quadrilateral mesh), --load-deck d.toml
// (CycleConfig::Load; prints [cycle.snapshots.volume]).
//
// Usage (from miniapps/seas):
//   mpirun -np N ./seas_test_snapshots_volume [--empty-box | --quads |
//                                              --load-deck d.toml]
// Exit code 0 = PASSED.

#include "mfem.hpp"
#include "../../cycle/cycle_config.hpp"
#include "../../cycle/snapshots.hpp"

#include <hdf5.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace
{

int g_failures = 0;
int g_rank = 0;

void Check(bool ok, const std::string &name, const std::string &detail = "")
{
   if (!ok) { g_failures++; }
   if (g_rank == 0)
   {
      std::cout << "  " << (ok ? "PASSED" : "FAILED") << ": " << name
                << (detail.empty() ? "" : " (" + detail + ")") << "\n";
   }
}

std::string Tmp(const std::string &name)
{
   const char *t = std::getenv("TMPDIR");
   return std::string(t ? t : "/tmp") + "/" + name;
}

/// The displacement at x and time scale s: degree p in x (exact in the DG
/// space of order p); the velocity is the same with other coefficients.
void Field(const Vector &x, int p, double s, double shift, Vector &u)
{
   const int dim = x.Size();
   for (int d = 0; d < dim; d++)
   {
      double v = 0.3 + shift + 0.1 * d;
      for (int k = 0; k < dim; k++) { v += (1.0 + d + 2 * k + shift) * 1e-4 * x(k); }
      if (p >= 2) { v += (0.5 + d + shift) * 1e-8 * x(0) * x(dim - 1); }
      u(d) = s * v;
   }
}

// ---------------------------------------------------------------- HDF5 reads
std::vector<double> ReadD(hid_t f, const std::string &name)
{
   const hid_t d = H5Dopen2(f, name.c_str(), H5P_DEFAULT);
   MFEM_VERIFY(d >= 0, "cannot open " << name);
   const hid_t sp = H5Dget_space(d);
   std::vector<double> v(static_cast<size_t>(H5Sget_simple_extent_npoints(sp)));
   MFEM_VERIFY(H5Dread(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()) >= 0,
               "cannot read " << name);
   H5Sclose(sp);
   H5Dclose(d);
   return v;
}

std::vector<long long> ReadI(hid_t f, const std::string &name)
{
   const hid_t d = H5Dopen2(f, name.c_str(), H5P_DEFAULT);
   MFEM_VERIFY(d >= 0, "cannot open " << name);
   const hid_t sp = H5Dget_space(d);
   std::vector<long long> v(static_cast<size_t>(H5Sget_simple_extent_npoints(sp)));
   MFEM_VERIFY(H5Dread(d, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()) >= 0,
               "cannot read " << name);
   H5Sclose(sp);
   H5Dclose(d);
   return v;
}

long long AttrNSteps(hid_t f)
{
   const hid_t g = H5Gopen2(f, "/VTKHDF/Steps", H5P_DEFAULT);
   const hid_t a = H5Aopen(g, "NSteps", H5P_DEFAULT);
   long long n = -1;
   H5Aread(a, H5T_NATIVE_LLONG, &n);
   H5Aclose(a);
   H5Gclose(g);
   return n;
}

/// Every dataset the restart check compares.
const char *kDatasets[] = {"/VTKHDF/Points", "/VTKHDF/Connectivity", "/VTKHDF/Offsets",
                           "/VTKHDF/PointData/displacement", "/VTKHDF/PointData/velocity",
                           "/VTKHDF/Steps/Values", "/VTKHDF/Steps/PointDataOffsets/displacement",
                           "/VTKHDF/Steps/PointDataOffsets/velocity", "/VTKHDF/Steps/PartOffsets",
                           "/VTKHDF/Steps/NumberOfParts", "/VTKHDF/Steps/PointOffsets",
                           "/VTKHDF/Steps/CellOffsets", "/VTKHDF/Steps/ConnectivityIdOffsets",
                           "/fd_cycle/regime", "/fd_cycle/phase", "/fd_cycle/reason",
                           "/fd_cycle/step"};

/// The raw bytes of a dataset (as doubles; every dataset here is 8 bytes wide).
std::vector<double> RawD(hid_t f, const std::string &name)
{
   const hid_t d = H5Dopen2(f, name.c_str(), H5P_DEFAULT);
   MFEM_VERIFY(d >= 0, "cannot open " << name);
   const hid_t ty = H5Dget_type(d);
   const hid_t sp = H5Dget_space(d);
   std::vector<double> v(static_cast<size_t>(H5Sget_simple_extent_npoints(sp)));
   const hid_t mt = H5Tget_native_type(ty, H5T_DIR_ASCEND);
   MFEM_VERIFY(H5Tget_size(mt) == 8, name << " is not 8 bytes wide");
   MFEM_VERIFY(H5Dread(d, mt, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()) >= 0, "read " << name);
   H5Tclose(mt);
   H5Sclose(sp);
   H5Tclose(ty);
   H5Dclose(d);
   return v;
}

// ----------------------------------------------------------- the writer runs
struct Grid
{
   std::unique_ptr<ParMesh> pmesh;
   std::unique_ptr<L2_FECollection> fec;
   std::unique_ptr<ParFiniteElementSpace> fes;
};

Grid MakeGrid(MPI_Comm comm, int dim, int p, bool quads = false)
{
   Mesh m = (dim == 2)
            ? Mesh::MakeCartesian2D(4, 3, quads ? Element::QUADRILATERAL : Element::TRIANGLE,
                                    true, 4000.0, 3000.0)
            : Mesh::MakeCartesian3D(2, 2, 2, Element::TETRAHEDRON, 2000.0, 2000.0, 2000.0);
   Grid g;
   g.pmesh = std::make_unique<ParMesh>(comm, m);
   g.fec = std::make_unique<L2_FECollection>(p, dim, BasisType::GaussLobatto);
   g.fes = std::make_unique<ParFiniteElementSpace>(g.pmesh.get(), g.fec.get(), dim,
                                                   Ordering::byNODES);
   return g;
}

cycle::CycleConfig VolumeConfig(const std::string &dir, const std::string &prefix, int p)
{
   cycle::CycleConfig cfg;
   cfg.shared.mesh.order = p;
   cfg.output.dir = dir;
   cfg.output.prefix = prefix;
   cfg.snapshots.volume.enabled = true;
   cfg.snapshots.volume.dt_qd = 10.0;
   cfg.snapshots.volume.dt_dyn = 0.1;
   return cfg;
}

/// One call of the sequence: QD (regime 0) or Dyn (1) at time t (tau in the
/// phase), with the fields projected at t.
struct Call
{
   int regime;
   double t, tau;
   cycle::SnapReason reason;   ///< QD: Cadence or ToDynamic; Dyn: end = EndOfDynamic
   bool enter = false;         ///< EnterDynamic before it
};

const std::vector<Call> kCalls = {
   {0, 0.0, 0.0, cycle::SnapReason::Cadence},
   {0, 5.0, 0.0, cycle::SnapReason::Cadence},
   {0, 12.0, 0.0, cycle::SnapReason::Cadence},
   {0, 25.0, 0.0, cycle::SnapReason::ToDynamic},
   {1, 25.05, 0.05, cycle::SnapReason::Cadence, true},
   {1, 25.1, 0.1, cycle::SnapReason::Cadence},
   {1, 25.15, 0.15, cycle::SnapReason::Cadence},
   {1, 25.25, 0.25, cycle::SnapReason::Cadence},
   {1, 25.3, 0.3, cycle::SnapReason::EndOfDynamic},
   {0, 26.0, 0.0, cycle::SnapReason::Cadence},
   {0, 31.0, 0.0, cycle::SnapReason::Cadence}};
const std::vector<double> kTimes = {0.0, 12.0, 25.0, 25.1, 25.25, 25.3, 31.0};
const std::vector<long long> kRegimes = {0, 0, 0, 1, 1, 1, 0};
const std::vector<long long> kReasons = {0, 0, 2, 0, 0, 3, 0};
/// The fields at time t are the polynomials times this.
double TimeScale(double t) { return 1.0 + 0.01 * t; }

/// Runs calls [from, to) on `snaps`; `state_after_12` (if given) gets the
/// writer's state right after the call at t = 12 s.
void RunCalls(cycle::Snapshots &snaps, Grid &g, int p, size_t from, size_t to,
              std::vector<cycle::SnapshotFileState> *state_after_12)
{
   const int dim = g.pmesh->Dimension();
   ParGridFunction u(g.fes.get()), v(g.fes.get());
   auto fill = [](cycle::StationFields &) {};
   for (size_t i = from; i < to; i++)
   {
      const Call &c = kCalls[i];
      const double s = TimeScale(c.t);
      VectorFunctionCoefficient uc(dim, [p, s](const Vector &x, Vector &y) { Field(x, p, s, 0.0, y); });
      VectorFunctionCoefficient vc(dim, [p, s](const Vector &x, Vector &y) { Field(x, p, s, 1.0, y); });
      u.ProjectCoefficient(uc);
      v.ProjectCoefficient(vc);
      if (c.regime == 0)
      {
         const Vector *uu = &u, *vv = &v;
         cycle::VolumeFieldsFn vol = [uu, vv](const Vector *&a, const Vector *&b) { a = uu; b = vv; };
         snaps.QD(c.t, static_cast<long long>(i), 0, c.reason, fill, vol);
      }
      else
      {
         if (c.enter) { snaps.EnterDynamic(); }
         snaps.Dyn(c.t, c.tau, static_cast<long long>(i), 1, fill, u, v,
                   c.reason == cycle::SnapReason::EndOfDynamic);
      }
      if (state_after_12 && c.t == 12.0) { *state_after_12 = snaps.GetState(); }
   }
}

void TestGrid(MPI_Comm comm, int dim, int p)
{
   const std::string name = "grid_" + std::to_string(dim) + "d_p" + std::to_string(p);
   if (g_rank == 0) { std::cout << "\n=== " << name << " ===\n"; }
   Grid g = MakeGrid(comm, dim, p);
   const std::string dir = Tmp("seas_test_snapshots_volume");
   if (g_rank == 0) { std::filesystem::create_directories(dir); }
   MPI_Barrier(comm);
   const std::string prefix = name;
   const cycle::CycleConfig cfg = VolumeConfig(dir, prefix, p);
   cycle::FaultDOFTable table;
   std::vector<cycle::SnapshotFileState> st12;
   long long npts = 0, ncells = 0;
   {
      cycle::Snapshots snaps(cfg, comm, *g.pmesh, *g.fes, table, 0, "test", nullptr);
      npts = snaps.VolumePoints();
      ncells = snaps.VolumeCells();
      RunCalls(snaps, g, p, 0, kCalls.size(), &st12);
      Check(snaps.VolumeSteps() == 7, name + ": 7 volume snapshots",
            std::to_string(snaps.VolumeSteps()));
      Check(snaps.FaultSteps() == 0 && snaps.SurfaceSteps() == 0 && snaps.GetState().size() == 1,
            name + ": no fault or surface file, one state record (volume)");
   }
   const std::string path = dir + "/" + prefix + "_volume_snapshots.vtkhdf";
   const long long ne = g.pmesh->GetGlobalNE();
   const int ppe = (dim == 2) ? (p + 1) * (p + 2) / 2 : (p + 1) * (p + 2) * (p + 3) / 6;
   const int cpe = (dim == 2) ? p * p : p * p * p;
   const int per = dim + 1;
   if (g_rank == 0)
   {
      const hid_t f = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      Check(f >= 0, name + ": the file opens");
      const std::vector<double> X = ReadD(f, "/VTKHDF/Points");
      const std::vector<long long> conn = ReadI(f, "/VTKHDF/Connectivity");
      const std::vector<long long> off = ReadI(f, "/VTKHDF/Offsets");
      const std::vector<long long> types = ReadI(f, "/VTKHDF/Types");
      const long long np = ReadI(f, "/VTKHDF/NumberOfPoints")[0];
      const long long nc = ReadI(f, "/VTKHDF/NumberOfCells")[0];
      Check(np == ne * ppe && npts == np, name + ": " + std::to_string(ppe) +
            " points per element", std::to_string(np) + " points, " + std::to_string(ne) +
            " elements");
      Check(nc == ne * cpe && ncells == nc, name + ": " + std::to_string(cpe) +
            " cells per element", std::to_string(nc) + " cells");
      const long long want_type = (dim == 2) ? 5 : 10;
      bool ty_ok = static_cast<long long>(types.size()) == nc &&
                   std::all_of(types.begin(), types.end(),
                               [want_type](long long t) { return t == want_type; });
      Check(ty_ok, name + ": every cell " + std::string(dim == 2 ? "VTK_TRIANGLE" : "VTK_TETRA"));
      bool off_ok = static_cast<long long>(off.size()) == nc + 1;
      for (long long c = 0; off_ok && c <= nc; c++) { off_ok = (off[c] == per * c); }
      Check(off_ok && static_cast<long long>(conn.size()) == per * nc,
            name + ": offsets " + std::to_string(per) + " per cell");
      // Each cell's points in one element block (blocks of ppe points, cpe
      // cells each, in order); no degenerate cell.
      bool block_ok = true, measure_ok = true;
      double zmax = 0.0;
      for (long long c = 0; c < nc; c++)
      {
         const long long b = c / cpe;
         double P[4][3];
         for (int k = 0; k < per; k++)
         {
            const long long j = conn[per * c + k];
            block_ok = block_ok && j >= b * ppe && j < (b + 1) * ppe;
            for (int d = 0; d < 3; d++) { P[k][d] = X[3 * j + d]; }
         }
         double meas = 0.0;
         if (dim == 2)
         {
            meas = 0.5 * std::abs((P[1][0] - P[0][0]) * (P[2][1] - P[0][1]) -
                                  (P[2][0] - P[0][0]) * (P[1][1] - P[0][1]));
         }
         else
         {
            double a[3], bb[3], cc[3];
            for (int d = 0; d < 3; d++)
            {
               a[d] = P[1][d] - P[0][d]; bb[d] = P[2][d] - P[0][d]; cc[d] = P[3][d] - P[0][d];
            }
            meas = std::abs(a[0] * (bb[1] * cc[2] - bb[2] * cc[1]) -
                            a[1] * (bb[0] * cc[2] - bb[2] * cc[0]) +
                            a[2] * (bb[0] * cc[1] - bb[1] * cc[0])) / 6.0;
         }
         measure_ok = measure_ok && meas > 1e-6;
      }
      if (dim == 2) { for (long long j = 0; j < np; j++) { zmax = std::max(zmax, std::abs(X[3 * j + 2])); } }
      Check(block_ok, name + ": each cell's points in its element's block");
      Check(measure_ok, name + ": no degenerate cell");
      if (dim == 2) { Check(zmax == 0.0, name + ": points at z = 0"); }
      // Steps and /fd_cycle.
      const long long ns = AttrNSteps(f);
      const std::vector<double> t = ReadD(f, "/VTKHDF/Steps/Values");
      const std::vector<long long> reg = ReadI(f, "/fd_cycle/regime");
      const std::vector<long long> rsn = ReadI(f, "/fd_cycle/reason");
      Check(ns == 7 && t == kTimes, name + ": the 7 snapshot times 0, 12, 25, 25.1, 25.25, 25.3, 31 s");
      Check(reg == kRegimes && rsn == kReasons,
            name + ": regimes QD QD QD dyn dyn dyn QD, reasons cadence cadence to-dynamic "
            "cadence cadence end-of-dynamic cadence");
      bool so_ok = true;
      for (const char *n : {"displacement", "velocity"})
      {
         const std::vector<long long> o = ReadI(f, std::string("/VTKHDF/Steps/PointDataOffsets/") + n);
         so_ok = so_ok && static_cast<long long>(o.size()) == ns;
         for (long long k = 0; so_ok && k < ns; k++) { so_ok = (o[k] == k * np); }
      }
      for (const char *n : {"PartOffsets", "PointOffsets", "CellOffsets", "ConnectivityIdOffsets"})
      {
         const std::vector<long long> o = ReadI(f, std::string("/VTKHDF/Steps/") + n);
         so_ok = so_ok && static_cast<long long>(o.size()) == ns &&
                 std::all_of(o.begin(), o.end(), [](long long x) { return x == 0; });
      }
      const std::vector<long long> parts = ReadI(f, "/VTKHDF/Steps/NumberOfParts");
      so_ok = so_ok && std::all_of(parts.begin(), parts.end(), [](long long x) { return x == 1; });
      Check(so_ok, name + ": step offsets (PointDataOffsets k * npts, one part)");
      // Values.
      const std::vector<double> U = ReadD(f, "/VTKHDF/PointData/displacement");
      const std::vector<double> V = ReadD(f, "/VTKHDF/PointData/velocity");
      Check(static_cast<long long>(U.size()) == 3 * ns * np && U.size() == V.size(),
            name + ": PointData 3 x npts per snapshot");
      double err = 0.0, big = 0.0, zc = 0.0;
      Vector x(dim), uu(dim), vv(dim);
      for (long long k = 0; k < ns && static_cast<long long>(U.size()) == 3 * ns * np; k++)
      {
         const double s = TimeScale(t[k]);
         for (long long j = 0; j < np; j++)
         {
            for (int d = 0; d < dim; d++) { x(d) = X[3 * j + d]; }
            Field(x, p, s, 0.0, uu);
            Field(x, p, s, 1.0, vv);
            const size_t b = static_cast<size_t>(3 * (k * np + j));
            for (int d = 0; d < dim; d++)
            {
               err = std::max(err, std::abs(U[b + d] - uu(d)));
               err = std::max(err, std::abs(V[b + d] - vv(d)));
               big = std::max(big, std::max(std::abs(uu(d)), std::abs(vv(d))));
            }
            if (dim == 2) { zc = std::max(zc, std::max(std::abs(U[b + 2]), std::abs(V[b + 2]))); }
         }
      }
      std::ostringstream o;
      o << "max error " << err << " of " << big;
      Check(err <= 1e-12 * big, name + ": displacement and velocity = the fields at every "
            "point of every snapshot", o.str());
      if (dim == 2) { Check(zc == 0.0, name + ": zero z components"); }
      H5Fclose(f);
      std::filesystem::copy_file(path, path + ".straight",
                                 std::filesystem::copy_options::overwrite_existing);
   }
   MPI_Barrier(comm);
   // Restart from the state after the 12 s snapshot: the same calls after it.
   {
      Check(st12.size() == 1 && st12[0].kind == "volume" && st12[0].steps == 2,
            name + ": the state after 12 s: volume, 2 snapshots");
      cycle::Snapshots snaps(cfg, comm, *g.pmesh, *g.fes, table, 0, "test", &st12);
      Check(snaps.VolumeSteps() == 2, name + ": the restarted writer starts at 2 snapshots");
      RunCalls(snaps, g, p, 3, kCalls.size(), nullptr);
      Check(snaps.VolumeSteps() == 7, name + ": the restarted writer ends at 7 snapshots");
   }
   if (g_rank == 0)
   {
      const hid_t a = H5Fopen((path + ".straight").c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      const hid_t b = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      bool same = (a >= 0 && b >= 0 && AttrNSteps(a) == AttrNSteps(b));
      std::string first_diff;
      for (const char *n : kDatasets)
      {
         const std::vector<double> x = RawD(a, n), y = RawD(b, n);
         const bool eq = x.size() == y.size() &&
                         std::equal(x.begin(), x.end(), y.begin(), [](double p1, double p2)
                         { return std::memcmp(&p1, &p2, sizeof(double)) == 0; });
         if (!eq && first_diff.empty()) { first_diff = n; }
         same = same && eq;
      }
      H5Fclose(a);
      H5Fclose(b);
      Check(same, name + ": restarted file = straight file, every dataset bit for bit",
            first_diff.empty() ? "" : "first difference " + first_diff);
      std::filesystem::remove(path + ".straight");
      std::filesystem::remove(path);
   }
   MPI_Barrier(comm);
}

void TestBox(MPI_Comm comm, int dim)
{
   const std::string name = "box_" + std::to_string(dim) + "d";
   if (g_rank == 0) { std::cout << "\n=== " << name << " ===\n"; }
   const int p = 1;
   Grid g = MakeGrid(comm, dim, p);
   const std::string dir = Tmp("seas_test_snapshots_volume");
   if (g_rank == 0) { std::filesystem::create_directories(dir); }
   MPI_Barrier(comm);
   cycle::CycleConfig cfg = VolumeConfig(dir, name, p);
   cfg.snapshots.volume_box = true;
   // 2-D: x in [500, 2600], y in [0, 1500]; 3-D: the half x <= 1000 of the box.
   cfg.snapshots.volume_box_min = (dim == 2) ? std::array<real_t, 3>{{500.0, 0.0, 0.0}}
                                             : std::array<real_t, 3>{{0.0, 0.0, 0.0}};
   cfg.snapshots.volume_box_max = (dim == 2) ? std::array<real_t, 3>{{2600.0, 1500.0, 0.0}}
                                             : std::array<real_t, 3>{{1000.0, 2000.0, 2000.0}};
   const auto &lo = cfg.snapshots.volume_box_min, &hi = cfg.snapshots.volume_box_max;
   // The expected elements, from the vertices (centroid = vertex average).
   long long want = 0;
   for (int e = 0; e < g.pmesh->GetNE(); e++)
   {
      Array<int> vv;
      g.pmesh->GetElementVertices(e, vv);
      double c[3] = {0.0, 0.0, 0.0};
      for (int k = 0; k < vv.Size(); k++)
      {
         const real_t *q = g.pmesh->GetVertex(vv[k]);
         for (int d = 0; d < dim; d++) { c[d] += q[d] / vv.Size(); }
      }
      bool in = true;
      for (int d = 0; d < dim; d++) { in = in && c[d] >= lo[d] && c[d] <= hi[d]; }
      want += in ? 1 : 0;
   }
   MPI_Allreduce(MPI_IN_PLACE, &want, 1, MPI_LONG_LONG, MPI_SUM, comm);
   cycle::FaultDOFTable table;
   long long npts = 0;
   {
      cycle::Snapshots snaps(cfg, comm, *g.pmesh, *g.fes, table, 0, "test", nullptr);
      npts = snaps.VolumePoints();
      RunCalls(snaps, g, p, 0, 3, nullptr);
   }
   const int ppe = dim + 1;
   Check(want > 0 && want < g.pmesh->GetGlobalNE() && npts == want * ppe,
         name + ": the elements whose centroid is in the box",
         std::to_string(npts / ppe) + " of " + std::to_string(g.pmesh->GetGlobalNE()) +
         " elements, expected " + std::to_string(want));
   if (g_rank == 0)
   {
      const std::string path = dir + "/" + name + "_volume_snapshots.vtkhdf";
      const hid_t f = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      const std::vector<double> X = ReadD(f, "/VTKHDF/Points");
      bool in_ok = true;
      for (long long b = 0; b < npts / ppe; b++)
      {
         double c[3] = {0.0, 0.0, 0.0};
         for (int k = 0; k < ppe; k++)
         {
            for (int d = 0; d < 3; d++) { c[d] += X[3 * (b * ppe + k) + d] / ppe; }
         }
         for (int d = 0; d < dim; d++)
         {
            in_ok = in_ok && c[d] >= lo[d] - 1e-9 && c[d] <= hi[d] + 1e-9;
         }
      }
      H5Fclose(f);
      Check(in_ok, name + ": every element block's centroid in the box");
      std::filesystem::remove(path);
   }
   MPI_Barrier(comm);
}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   g_rank = Mpi::WorldRank();
   const char *deck = "";
   bool empty_box = false, quads = false;
   OptionsParser args(argc, argv);
   args.AddOption(&deck, "-ld", "--load-deck", "Only load this deck and print [cycle.snapshots.volume].");
   args.AddOption(&empty_box, "-eb", "--empty-box", "-no-eb", "--no-empty-box",
                  "Only build a volume writer whose box holds no element centroid.");
   args.AddOption(&quads, "-q", "--quads", "-no-q", "--no-quads",
                  "Only build a volume writer on a quadrilateral mesh.");
   args.Parse();
   if (!args.Good())
   {
      if (g_rank == 0) { args.PrintUsage(std::cout); }
      return 1;
   }
   if (std::string(deck).size() > 0)
   {
      const cycle::CycleConfig cfg = cycle::CycleConfig::Load(deck);
      const cycle::SnapshotConfig &s = cfg.snapshots;
      if (g_rank == 0)
      {
         std::cout << "[cycle.snapshots.volume] enabled " << s.volume.enabled << ", dt_qd "
                   << s.volume.dt_qd << ", dt_dyn " << s.volume.dt_dyn << ", box "
                   << s.volume_box << " [" << s.volume_box_min[0] << ", " << s.volume_box_min[1]
                   << ", " << s.volume_box_min[2] << "] .. [" << s.volume_box_max[0] << ", "
                   << s.volume_box_max[1] << ", " << s.volume_box_max[2] << "]\n";
      }
      return 0;
   }
   if (empty_box || quads)
   {
      Grid g = MakeGrid(MPI_COMM_WORLD, 2, 1, quads);
      cycle::CycleConfig cfg = VolumeConfig(Tmp("seas_test_snapshots_volume"), "refuse", 1);
      if (empty_box)
      {
         cfg.snapshots.volume_box = true;
         cfg.snapshots.volume_box_min = {{-500.0, -500.0, 0.0}};
         cfg.snapshots.volume_box_max = {{-100.0, -100.0, 0.0}};
      }
      cycle::FaultDOFTable table;
      cycle::Snapshots snaps(cfg, MPI_COMM_WORLD, *g.pmesh, *g.fes, table, 0, "test", nullptr);
      return 0;   // not reached: the writer refuses
   }
   if (g_rank == 0)
   {
      std::cout << "test_snapshots_volume: ranks " << Mpi::WorldSize() << "\n";
   }
   for (int dim : {2, 3})
   {
      for (int p : {1, 2}) { TestGrid(MPI_COMM_WORLD, dim, p); }
   }
   TestBox(MPI_COMM_WORLD, 2);
   TestBox(MPI_COMM_WORLD, 3);
   if (g_rank == 0)
   {
      std::cout << "\ntest_snapshots_volume: " << (g_failures ? "FAILED" : "PASSED") << " ("
                << g_failures << " failing checks)\n";
   }
   return g_failures ? 1 : 0;
}

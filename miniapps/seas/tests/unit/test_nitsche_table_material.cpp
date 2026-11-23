// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_nitsche_table_material.cpp — the 3-D table material of the Nitsche
// driver (dynamic/nitsche_table_material.hpp; TPV34 plan rev 2, step S2).
//
//   T1  a table of functions linear in (x, y, z): the element-centroid sample
//       (NitscheMaterialSample) equals rho, rho Vs^2, rho (Vp^2 - 2 Vs^2) of
//       the functions at the centroid, with the table in the mesh frame.
//   T2  the same with a rotated, shifted frame ([material.frame]).
//   T3  a table that depends on depth only, piecewise linear with kinks on
//       grid planes (a layered profile): exact at every centroid.
//   T4  mu_at / eval_at_xyz at points; ContainsMeshBox true for the mesh,
//       false for a larger box.
//   T5  a point outside the table: aborts with "abort" (in a child process
//       without MPI), takes the edge value with "edge_value".
//   T6  crs "mesh-local" and "EPSG:32611" are read; another label aborts.
//
// Run on 1 and 2 ranks (make test-nitsche-table-material); with 2 ranks the
// table is held once per node (shared-memory mode).

#include "mfem.hpp"

#include "../../domain/nitsche_material_sample.hpp"
#include "../../dynamic/nitsche_table_material.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include <hdf5.h>
#include <unistd.h>

using namespace mfem;
using namespace mfem::seas;

namespace
{

int num_passed = 0, num_failed = 0;

void Check(bool ok, const std::string &what)
{
   if (ok) { num_passed++; }
   else { num_failed++; }
   if (Mpi::Root()) { std::cout << (ok ? "  PASSED: " : "FAILED: ") << what << "\n"; }
}

void StrAttr(hid_t loc, const char *name, const std::string &v)
{
   hid_t sp = H5Screate(H5S_SCALAR), ty = H5Tcopy(H5T_C_S1);
   H5Tset_size(ty, v.size() + 1);
   H5Tset_strpad(ty, H5T_STR_NULLTERM);
   hid_t a = H5Acreate2(loc, name, ty, sp, H5P_DEFAULT, H5P_DEFAULT);
   H5Awrite(a, ty, v.c_str());
   H5Aclose(a); H5Tclose(ty); H5Sclose(sp);
}
void DblAttr(hid_t loc, const char *name, double v)
{
   hid_t sp = H5Screate(H5S_SCALAR);
   hid_t a = H5Acreate2(loc, name, H5T_NATIVE_DOUBLE, sp, H5P_DEFAULT, H5P_DEFAULT);
   H5Awrite(a, H5T_NATIVE_DOUBLE, &v);
   H5Aclose(a); H5Sclose(sp);
}

using Fn = std::function<double(double, double, double)>;

/// A data_projection_v1 table with Vp, Vs, density sampled on the axes.
void WriteTable(const std::string &path, const std::string &crs, const std::vector<double> &x,
                const std::vector<double> &y, const std::vector<double> &z, const Fn &vp,
                const Fn &vs, const Fn &rho)
{
   hid_t f = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   StrAttr(f, "schema_version", "data_projection_v1");
   StrAttr(f, "crs", crs);
   StrAttr(f, "units", "m");
   StrAttr(f, "z_positive", "elevation");
   StrAttr(f, "created_at", "test");
   hid_t g = H5Gcreate2(f, "/grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   const std::vector<double> *ax[3] = {&x, &y, &z};
   const char *an[3] = {"x", "y", "z"};
   for (int i = 0; i < 3; i++)
   {
      hsize_t n = ax[i]->size();
      hid_t sp = H5Screate_simple(1, &n, nullptr);
      hid_t d = H5Dcreate2(g, an[i], H5T_NATIVE_DOUBLE, sp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ax[i]->data());
      H5Dclose(d); H5Sclose(sp);
   }
   H5Gclose(g);
   g = H5Gcreate2(f, "/fields", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   const Fn *fn[3] = {&vp, &vs, &rho};
   const char *fname[3] = {"Vp", "Vs", "density"};
   const char *unit[3] = {"m/s", "m/s", "kg/m3"};
   hsize_t shape[3] = {x.size(), y.size(), z.size()};
   for (int k = 0; k < 3; k++)
   {
      std::vector<double> v(x.size() * y.size() * z.size());
      double lo = 1e300, hi = -1e300;
      size_t n = 0;
      for (double xi : x) { for (double yi : y) { for (double zi : z)
      {
         v[n] = (*fn[k])(xi, yi, zi);
         lo = std::min(lo, v[n]); hi = std::max(hi, v[n]);
         n++;
      } } }
      hid_t sp = H5Screate_simple(3, shape, nullptr);
      hid_t d = H5Dcreate2(g, fname[k], H5T_NATIVE_DOUBLE, sp, H5P_DEFAULT, H5P_DEFAULT,
                           H5P_DEFAULT);
      H5Dwrite(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
      StrAttr(d, "units", unit[k]);
      DblAttr(d, "min_value", lo);
      DblAttr(d, "max_value", hi);
      H5Dclose(d); H5Sclose(sp);
   }
   H5Gclose(g);
   H5Fclose(f);
}

std::vector<double> Axis(double a, double b, int n)
{
   std::vector<double> v(n + 1);
   for (int i = 0; i <= n; i++) { v[i] = a + (b - a) * i / n; }
   return v;
}

/// Largest relative error of the centroid sample against the functions
/// evaluated at the table point of each centroid (all ranks).
real_t SampleError(ParMesh &pmesh, const NitscheTableMaterial &tm, const Fn &vp, const Fn &vs,
                   const Fn &rho)
{
   NitscheMaterialSample smp(pmesh, tm.field);
   real_t err = 0.0;
   Vector c(3);
   for (int e = 0; e < pmesh.GetNE(); e++)
   {
      ElementTransformation *T = pmesh.GetElementTransformation(e);
      const IntegrationPoint &ip = Geometries.GetCenter(pmesh.GetElementBaseGeometry(e));
      T->Transform(ip, c);
      real_t q[3];
      tm.frame.Map(c.GetData(), q);
      const real_t r = rho(q[0], q[1], q[2]), a = vp(q[0], q[1], q[2]), b = vs(q[0], q[1], q[2]);
      const real_t ref[3] = {r * (a * a - 2.0 * b * b), r * b * b, r};
      const real_t got[3] = {smp.Lambda().Eval(*T, ip), smp.Mu().Eval(*T, ip),
                             smp.Rho().Eval(*T, ip)};
      for (int k = 0; k < 3; k++) { err = std::max(err, std::abs(got[k] - ref[k]) / std::abs(ref[k])); }
   }
   MPI_Allreduce(MPI_IN_PLACE, &err, 1, MFEM_MPI_REAL_T, MPI_MAX, MPI_COMM_WORLD);
   return err;
}

/// A case that may abort, run by a child process that never starts MPI (so
/// MFEM_ABORT ends in abort(), not MPI_Abort): 0 = a point outside the table
/// with "abort", 1 = only read the table.
int ChildCase(int which, const std::string &table)
{
   auto tm = MakeNitscheTableMaterial(table, TableFrame(), false);
   if (which == 0)
   {
      const real_t xout[3] = {9e3, 0.0, -1e3};
      real_t l, m, r;
      tm->eval_at_xyz(xout, l, m, r);
   }
   return 0;
}

/// True when the child case exits with a failure.
bool ChildAborts(const std::string &self, int which, const std::string &table)
{
   const std::string cmd = self + " --child " + std::to_string(which) + " " + table
                           + " > /dev/null 2>&1";
   return std::system(cmd.c_str()) != 0;
}

} // namespace

int main(int argc, char *argv[])
{
   if (argc == 4 && std::string(argv[1]) == "--child")
   {
      return ChildCase(std::atoi(argv[2]), argv[3]);
   }
   Mpi::Init(argc, argv);
   const int rank = Mpi::WorldRank();
   const std::string dir = "/tmp/seas_test_nitsche_table_material_" + std::to_string(getuid());
   const std::string t_lin = dir + "_lin.h5", t_rot = dir + "_rot.h5", t_lay = dir + "_lay.h5",
                     t_bad = dir + "_bad.h5", t_utm = dir + "_utm.h5";

   // mesh: box [-4, 4] x [-3, 3] x [-6, 0] km, tetrahedra
   Mesh smesh = Mesh::MakeCartesian3D(6, 5, 6, Element::TETRAHEDRON, 8e3, 6e3, 6e3);
   smesh.Transform([](const Vector &x, Vector &y) {
      y.SetSize(3);
      y(0) = x(0) - 4e3; y(1) = x(1) - 3e3; y(2) = x(2) - 6e3;
   });
   ParMesh pmesh(MPI_COMM_WORLD, smesh);

   const Fn vp = [](double x, double y, double z) { return 5000.0 + 0.02 * x - 0.03 * y - 0.25 * z; };
   const Fn vs = [](double x, double y, double z) { return 2800.0 - 0.01 * x + 0.02 * y - 0.15 * z; };
   const Fn rho = [](double x, double y, double z) { return 2500.0 + 0.004 * x + 0.003 * y - 0.05 * z; };

   // rotated frame: table point q_xy = o + ax X + ay Y (a rotation by 35 degrees)
   TableFrame rot;
   const double th = 35.0 * M_PI / 180.0;
   rot.origin[0] = 648446.0; rot.origin[1] = 3625237.0;
   rot.x_axis[0] = std::cos(th); rot.x_axis[1] = std::sin(th);
   rot.y_axis[0] = -std::sin(th); rot.y_axis[1] = std::cos(th);
   // functions of the table point, centred on the origin so the values stay physical
   const Fn vp_r = [&](double x, double y, double z) { return vp(x - rot.origin[0], y - rot.origin[1], z); };
   const Fn vs_r = [&](double x, double y, double z) { return vs(x - rot.origin[0], y - rot.origin[1], z); };
   const Fn rho_r = [&](double x, double y, double z) { return rho(x - rot.origin[0], y - rot.origin[1], z); };

   // layered: kinks at depths 1000 and 2500 m, both grid planes
   auto prof = [](double z, double v0, double g1, double g2, double g3) {
      const double d = -z;
      if (d <= 1000.0) { return v0 + g1 * d; }
      if (d <= 2500.0) { return v0 + g1 * 1000.0 + g2 * (d - 1000.0); }
      return v0 + g1 * 1000.0 + g2 * 1500.0 + g3 * (d - 2500.0);
   };
   const Fn vp_l = [&](double, double, double z) { return prof(z, 3000.0, 1.0, 0.2, 0.05); };
   const Fn vs_l = [&](double, double, double z) { return prof(z, 1400.0, 0.6, 0.1, 0.03); };
   const Fn rho_l = [&](double, double, double z) { return prof(z, 2200.0, 0.2, 0.05, 0.01); };

   if (rank == 0)
   {
      // non-uniform axes on purpose (two spacing zones)
      std::vector<double> x = Axis(-4e3, 0.0, 4), y = Axis(-3e3, 3e3, 5), z = Axis(-6e3, 0.0, 9);
      for (double v : Axis(500.0, 4e3, 7)) { x.push_back(v); }
      WriteTable(t_lin, "mesh-local", x, y, z, vp, vs, rho);
      WriteTable(t_utm, "EPSG:32611", x, y, z, vp, vs, rho);
      WriteTable(t_bad, "EPSG:4326", x, y, z, vp, vs, rho);
      WriteTable(t_rot, "EPSG:32611", Axis(rot.origin[0] - 6e3, rot.origin[0] + 6e3, 7),
                 Axis(rot.origin[1] - 6e3, rot.origin[1] + 6e3, 5), z, vp_r, vs_r, rho_r);
      // depth axis with planes at every 500 m (holds the kinks)
      WriteTable(t_lay, "mesh-local", Axis(-4e3, 4e3, 2), Axis(-3e3, 3e3, 2), Axis(-6e3, 0.0, 12),
                 vp_l, vs_l, rho_l);
   }
   MPI_Barrier(MPI_COMM_WORLD);

   MPI_Comm node_comm;
   MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, rank, MPI_INFO_NULL, &node_comm);
   const real_t tol = 1e-12;
   {
      auto tm = MakeNitscheTableMaterial(t_lin, TableFrame(), false, node_comm);
      const real_t e1 = SampleError(pmesh, *tm, vp, vs, rho);
      Check(e1 < tol, "T1 linear table, mesh frame: centroid sample exact (rel err "
            + std::to_string(e1) + ")");
      Check(tm->vp->Crs() == "mesh-local" && tm->frame.IsIdentity(), "T1 crs mesh-local, identity frame");

      // T4
      const real_t xp[3] = {1234.0, -2100.0, -3456.0};
      real_t l, m, r;
      tm->eval_at_xyz(xp, l, m, r);
      const real_t rr = rho(xp[0], xp[1], xp[2]), bb = vs(xp[0], xp[1], xp[2]),
                   aa = vp(xp[0], xp[1], xp[2]);
      Check(std::abs(r - rr) < tol * rr && std::abs(m - rr * bb * bb) < tol * m
            && std::abs(l - rr * (aa * aa - 2 * bb * bb)) < tol * l, "T4 eval_at_xyz at a point");
      Vector lo, hi;
      pmesh.GetBoundingBox(lo, hi);
      Check(tm->ContainsMeshBox(lo.GetData(), hi.GetData()), "T4 ContainsMeshBox: the mesh box");
      hi(0) += 1.0;
      Check(!tm->ContainsMeshBox(lo.GetData(), hi.GetData()), "T4 ContainsMeshBox: a larger box");
   }
   {
      auto tm = MakeNitscheTableMaterial(t_rot, rot, false, node_comm);
      const real_t e2 = SampleError(pmesh, *tm, vp_r, vs_r, rho_r);
      Check(e2 < 1e-10, "T2 linear table, rotated and shifted frame (rel err " + std::to_string(e2) + ")");
      Vector lo, hi;
      pmesh.GetBoundingBox(lo, hi);
      Check(tm->ContainsMeshBox(lo.GetData(), hi.GetData()), "T2 ContainsMeshBox through the frame");
   }
   {
      auto tm = MakeNitscheTableMaterial(t_lay, TableFrame(), false, node_comm);
      // the mesh has planes every 1000 m; the kinks (1000, 2500 m) are table
      // planes, and every centroid lies inside one 500 m table cell in z
      const real_t e3 = SampleError(pmesh, *tm, vp_l, vs_l, rho_l);
      Check(e3 < tol, "T3 layered depth profile with kinks on table planes (rel err "
            + std::to_string(e3) + ")");
   }
   {
      auto tm = MakeNitscheTableMaterial(t_utm, TableFrame(), true, node_comm);
      const real_t xin[3] = {4e3, 0.0, -1e3}, xout[3] = {9e3, 0.0, -1e3};
      real_t l0, m0, r0, l1, m1, r1;
      tm->eval_at_xyz(xin, l0, m0, r0);
      tm->eval_at_xyz(xout, l1, m1, r1);
      Check(l0 == l1 && m0 == m1 && r0 == r1, "T5 edge_value: outside point = edge point");
      Check(tm->vp->Crs() == "EPSG:32611", "T6 crs EPSG:32611 is read");
   }
   MPI_Comm_free(&node_comm);
   if (rank == 0 && Mpi::WorldSize() == 1)
   {
      Check(ChildAborts(argv[0], 0, t_lin), "T5 abort: a point outside the table aborts");
      Check(!ChildAborts(argv[0], 1, t_lin), "T5 control: reading the table does not abort");
      Check(ChildAborts(argv[0], 1, t_bad), "T6 another crs label aborts");
   }
   if (rank == 0)
   {
      for (const std::string &f : {t_lin, t_rot, t_lay, t_bad, t_utm}) { std::remove(f.c_str()); }
   }
   int tot[2] = {num_passed, num_failed};
   if (rank == 0)
   {
      std::cout << "test_nitsche_table_material (" << Mpi::WorldSize() << " rank(s)): " << tot[0]
                << " passed, " << tot[1] << " failed\n";
   }
   int fail = num_failed;
   MPI_Allreduce(MPI_IN_PLACE, &fail, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   return fail ? 1 : 0;
}

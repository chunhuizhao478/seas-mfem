// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_grid_field_2d.cpp — unit tests of io/grid_field_2d.hpp
// (heterogeneity plan rev4, phase H1, gates U1-U3).
//
//   make seas_test_grid_field_2d && ./seas_test_grid_field_2d [real_file.nc]
//
// U1  a linear field f = a + b x + c y is reproduced to 1e-12 (relative) at
//     random points; a constant field is reproduced bit for bit; grid nodes
//     return their stored values bit for bit; a non-uniform axis works.
// U2  outside the grid: Clamp returns the edge value; Contains() is false.
//     (Abort is MFEM_ABORT; it is checked in a forked child.)
// U3  the values read equal the values written, bit for bit (HDF5 file with
//     the fixed-length string attribute netCDF4-python writes); a wrong
//     schema is refused (forked child).
// With an argument: the real UCERF file is opened and its bbox and a node
// printed (the build_hetero_grid.py output).

#include "mfem.hpp"
#include "../../io/grid_field_2d.hpp"

#include <hdf5.h>
#include <sys/wait.h>
#include <unistd.h>

#include <cmath>
#include <cstdio>
#include <iostream>
#include <functional>
#include <random>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_failed = 0;
#define CHECK(cond, msg) do { num_tests++; if (!(cond)) { num_failed++; \
   std::cerr << "FAILED: " << msg << " (line " << __LINE__ << ")\n"; } \
   else { std::cout << "  PASSED: " << msg << "\n"; } } while (0)

static void WriteFile(const std::string &path, const std::vector<double> &x, const std::vector<double> &y,
                      const std::vector<std::pair<std::string, std::vector<double>>> &fields,
                      const std::string &schema)
{
   hid_t f = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   auto w1 = [&](const char *name, const std::vector<double> &v)
   {
      hsize_t n = v.size(); hid_t sp = H5Screate_simple(1, &n, nullptr);
      hid_t d = H5Dcreate2(f, name, H5T_IEEE_F64LE, sp, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()); H5Dclose(d); H5Sclose(sp);
   };
   w1("x", x); w1("y", y);
   for (const auto &fv : fields)
   {
      hsize_t dims[2] = {y.size(), x.size()}; hid_t sp = H5Screate_simple(2, dims, nullptr);
      hid_t pl = H5Pcreate(H5P_DATASET_CREATE); hsize_t ch[2] = {std::min<hsize_t>(4, dims[0]), std::min<hsize_t>(4, dims[1])};
      H5Pset_chunk(pl, 2, ch); H5Pset_deflate(pl, 4);                       // compressed, as netCDF4 zlib=True
      hid_t d = H5Dcreate2(f, fv.first.c_str(), H5T_IEEE_F64LE, sp, H5P_DEFAULT, pl, H5P_DEFAULT);
      H5Dwrite(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, fv.second.data());
      H5Dclose(d); H5Pclose(pl); H5Sclose(sp);
   }
   hid_t st = H5Tcopy(H5T_C_S1); H5Tset_size(st, schema.size()); H5Tset_strpad(st, H5T_STR_NULLTERM);
   hid_t asp = H5Screate(H5S_SCALAR);
   hid_t a = H5Acreate2(f, "schema", st, asp, H5P_DEFAULT, H5P_DEFAULT);
   H5Awrite(a, st, schema.c_str()); H5Aclose(a); H5Sclose(asp); H5Tclose(st);
   H5Fclose(f);
}

static bool AbortsInChild(const std::function<void()> &fn)
{
   std::cout.flush(); std::cerr.flush();
   pid_t pid = fork();
   if (pid == 0) { freopen("/dev/null", "w", stderr); freopen("/dev/null", "w", stdout); fn(); _exit(0); }
   int status = 0; waitpid(pid, &status, 0);
   return !(WIFEXITED(status) && WEXITSTATUS(status) == 0);
}

int main(int argc, char **argv)
{
   const std::string path = "/tmp/seas_test_grid_field_2d.h5";
   std::vector<double> x = {-3000.0, -1000.0, 0.0, 500.0, 2500.0, 4000.0};   // non-uniform
   std::vector<double> y = {-2000.0, 0.0, 1000.0, 3000.0};
   const double A = 3.2e10, Bx = 1.7e5, Cy = -2.3e5, K = 32038120320.0;
   std::vector<double> lin(x.size() * y.size()), cst(x.size() * y.size(), K), rnd(x.size() * y.size());
   std::mt19937_64 g(7); std::uniform_real_distribution<double> u(-1.0, 1.0);
   for (size_t j = 0; j < y.size(); j++)
      for (size_t i = 0; i < x.size(); i++)
      {
         lin[j * x.size() + i] = A + Bx * x[i] + Cy * y[j];
         rnd[j * x.size() + i] = 1e7 * u(g);
      }
   WriteFile(path, x, y, {{"lin", lin}, {"cst", cst}, {"rnd", rnd}}, "seas_grid2d_v1");

   std::cout << "U3: read back\n";
   GridField2D G(path, {"lin", "cst", "rnd"});
   CHECK(G.NumX() == 6 && G.NumY() == 4, "axis sizes");
   bool same = true;
   for (int j = 0; j < G.NumY(); j++)
      for (int i = 0; i < G.NumX(); i++)
      {
         same &= G.Node(0, i, j) == lin[j * x.size() + i];
         same &= G.Node(2, i, j) == rnd[j * x.size() + i];
      }
   CHECK(same, "stored values bit for bit (deflate-compressed datasets)");
   CHECK(G.Index("cst") == 1, "field index by name");

   std::cout << "U1: interpolation\n";
   bool nodes_exact = true;
   for (int j = 0; j < G.NumY(); j++)
      for (int i = 0; i < G.NumX(); i++)
      {
         nodes_exact &= G.Evaluate(2, x[i], y[j]) == rnd[j * x.size() + i];
      }
   CHECK(nodes_exact, "grid nodes return the stored value bit for bit");
   std::uniform_real_distribution<double> ux(x.front(), x.back()), uy(y.front(), y.back());
   double err = 0.0; bool cst_exact = true;
   for (int k = 0; k < 10000; k++)
   {
      const double px = ux(g), py = uy(g);
      const double f = G.Evaluate(0, px, py), fe = A + Bx * px + Cy * py;
      err = std::max(err, std::abs(f - fe) / std::abs(fe));
      cst_exact &= G.Evaluate(1, px, py) == K;
      double all[3]; G.EvaluateAll(px, py, all);
      cst_exact &= all[1] == K && all[0] == f;
   }
   std::cout << "    max relative error of the linear field: " << err << "\n";
   CHECK(err < 1e-12, "linear field to 1e-12 relative");
   CHECK(cst_exact, "constant field bit for bit at 10,000 random points; EvaluateAll == Evaluate");
   CHECK(G.Evaluate(0, x.back(), y.back()) == lin.back(), "the far corner (last cell) exact");

   std::cout << "U2: outside the grid\n";
   GridField2D Gc(path, {"lin"}, OOBPolicy::Clamp);
   CHECK(!Gc.Contains(x.back() + 10.0, 0.0), "Contains() false outside");
   CHECK(Gc.Evaluate(0, x.back() + 1e4, y.back() + 1e4) == lin.back(), "Clamp returns the corner value");
   CHECK(AbortsInChild([&] { GridField2D Ga(path, {"lin"}); Ga.Evaluate(0, x.front() - 1.0, 0.0); }),
         "Abort policy aborts outside the grid");
   const std::string bad = "/tmp/seas_test_grid_field_2d_bad.h5";
   WriteFile(bad, x, y, {{"lin", lin}}, "something_else");
   CHECK(AbortsInChild([&] { GridField2D Gb(bad, {"lin"}); }), "wrong schema is refused");
   CHECK(AbortsInChild([&] { GridField2D Gm(path, {"missing"}); }), "missing field is refused");

   if (argc > 1)
   {
      std::cout << "real file: " << argv[1] << "\n";
      GridField2D R(argv[1], {"vp", "vs", "rho", "s_xx", "s_yy", "s_xy", "s_zz"});
      const auto b = R.BBox();
      std::cout << "    " << R.NumX() << " x " << R.NumY() << ", x [" << b[0] << ", " << b[1] << "], y ["
                << b[2] << ", " << b[3] << "]\n";
      std::cout.precision(17);
      std::cout << "    vs(0,0) = " << R.Node(1, 0, 0) << ", s_zz(0,0) = " << R.Node(6, 0, 0) << "\n";
      CHECK(R.NumX() > 1 && R.NumY() > 1, "real file opens with its 7 fields");
   }

   std::remove(path.c_str()); std::remove(bad.c_str());
   std::cout << (num_failed ? "FAIL" : "PASS") << ": " << num_tests - num_failed << " / " << num_tests << "\n";
   return num_failed ? 1 : 0;
}

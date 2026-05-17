// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Unit tests for StressField3D (Phase 6 §1 of
// PLAN_onfaultstress.md).
//
// Builds a synthetic schema-v1 sidecar with the six canonical
// stress components, constructs a StressField3D, and asserts:
//   - all six fields load
//   - Evaluate() returns a 3x3 symmetric tensor
//   - sign convention: pass-through (no flip)
//   - SetInterpMode propagates to all six readers
//   - mismatched-grid sidecar aborts
//   - missing-field sidecar aborts
//   - bbox accessor and ContainsBBox
//
// The synthetic sidecar helpers below are written from scratch to
// avoid coupling to the test_data_field_3d.cpp fixtures.  Each
// component field can be assigned an independent linear-in-(x,y,z)
// function so the symmetry of Evaluate() can be verified through
// non-trivial off-diagonal entries.

#include "mfem.hpp"

#include "../../io/stress_field_3d.hpp"
#include "../../io/stress_field_coefficient.hpp"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#include <hdf5.h>
#include <sys/wait.h>
#include <unistd.h>

using namespace mfem;
using namespace mfem::seas;

// ----------------------------------------------------------------------
// Tiny test framework (matches test_data_field_3d.cpp)
// ----------------------------------------------------------------------

static int num_tests = 0;
static int num_passed = 0;
static int num_failed = 0;

#define TEST_ASSERT(condition, message) \
   do { \
      num_tests++; \
      if (!(condition)) { \
         std::cerr << "FAILED: " << message << " (line " << __LINE__ \
                   << ")\n"; \
         num_failed++; \
      } else { \
         std::cout << "  PASSED: " << message << "\n"; \
         num_passed++; \
      } \
   } while (0)

#define TEST_NEAR(value, expected, tol, message) \
   do { \
      num_tests++; \
      const real_t _v = (value); \
      const real_t _e = (expected); \
      const real_t _t = (tol); \
      if (std::abs(_v - _e) > _t) { \
         std::cerr << "FAILED: " << message << " (got " << _v \
                   << ", expected " << _e << ", tol " << _t \
                   << ", line " << __LINE__ << ")\n"; \
         num_failed++; \
      } else { \
         std::cout << "  PASSED: " << message << "\n"; \
         num_passed++; \
      } \
   } while (0)


// ----------------------------------------------------------------------
// HDF5 fixture helpers (mirrored from test_data_field_3d.cpp)
// ----------------------------------------------------------------------

namespace
{

void write_string_attr(hid_t loc, const char* name, const std::string& v)
{
   hid_t aspace = H5Screate(H5S_SCALAR);
   hid_t atype = H5Tcopy(H5T_C_S1);
   H5Tset_size(atype, v.size() + 1);
   H5Tset_strpad(atype, H5T_STR_NULLTERM);
   hid_t aid = H5Acreate2(loc, name, atype, aspace, H5P_DEFAULT, H5P_DEFAULT);
   H5Awrite(aid, atype, v.c_str());
   H5Aclose(aid);
   H5Tclose(atype);
   H5Sclose(aspace);
}

void write_double_attr(hid_t loc, const char* name, double v)
{
   hid_t aspace = H5Screate(H5S_SCALAR);
   hid_t aid = H5Acreate2(loc, name, H5T_NATIVE_DOUBLE, aspace,
                          H5P_DEFAULT, H5P_DEFAULT);
   H5Awrite(aid, H5T_NATIVE_DOUBLE, &v);
   H5Aclose(aid);
   H5Sclose(aspace);
}

void write_axis(hid_t group, const char* name, const std::vector<double>& v)
{
   hsize_t dim = v.size();
   hid_t sid = H5Screate_simple(1, &dim, nullptr);
   hid_t did = H5Dcreate2(group, name, H5T_NATIVE_DOUBLE, sid,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
   H5Dclose(did);
   H5Sclose(sid);
}

void write_field(hid_t flds,
                 const char* name,
                 const std::vector<double>& field_flat,
                 const std::vector<double>& x,
                 const std::vector<double>& y,
                 const std::vector<double>& z,
                 double min_value,
                 double max_value,
                 const std::string& units = "Pa")
{
   hsize_t shape[3] = { x.size(), y.size(), z.size() };
   hid_t sid = H5Screate_simple(3, shape, nullptr);
   hid_t did = H5Dcreate2(flds, name, H5T_NATIVE_DOUBLE, sid,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            field_flat.data());
   write_string_attr(did, "units", units);
   write_double_attr(did, "min_value", min_value);
   write_double_attr(did, "max_value", max_value);
   H5Dclose(did);
   H5Sclose(sid);
}

bool expect_abort(const std::function<void(void)>& f)
{
   pid_t pid = fork();
   if (pid == 0)
   {
      ::close(1);
      ::close(2);
      try { f(); }
      catch (...) { std::_Exit(1); }
      std::_Exit(0);
   }
   int status = 0;
   ::waitpid(pid, &status, 0);
   if (WIFSIGNALED(status))             { return true; }
   if (WIFEXITED(status) && WEXITSTATUS(status) != 0) { return true; }
   return false;
}

std::string make_tmp_path(const char* stem)
{
   const char* tmpdir = std::getenv("TMPDIR");
   if (!tmpdir || !*tmpdir) { tmpdir = "/tmp"; }
   std::string path = std::string(tmpdir) + "/seas_test_stress_field_3d_"
                      + std::string(stem) + "_"
                      + std::to_string(::getpid()) + ".h5";
   ::unlink(path.c_str());
   return path;
}

/// Fill a (Nx, Ny, Nz) row-major field with f(x, y, z) = a*x + b*y + c*z + d.
std::vector<double> linear_field(const std::vector<double>& x,
                                 const std::vector<double>& y,
                                 const std::vector<double>& z,
                                 double a, double b, double c, double d)
{
   std::vector<double> out(x.size() * y.size() * z.size(), 0.0);
   for (size_t i = 0; i < x.size(); ++i)
   {
      for (size_t j = 0; j < y.size(); ++j)
      {
         for (size_t k = 0; k < z.size(); ++k)
         {
            const size_t flat = (i * y.size() + j) * z.size() + k;
            out[flat] = a * x[i] + b * y[j] + c * z[k] + d;
         }
      }
   }
   return out;
}

/// Coefficients for the six per-component fields written to the
/// synthetic sidecar.  Each component is a different linear-in-(x,y,z)
/// function, so any silent sign flip or component swap is caught.
struct ComponentCoeffs
{
   double a{}, b{}, c{}, d{};
   double eval(double x, double y, double z) const
   {
      return a * x + b * y + c * z + d;
   }
};

struct StressSidecarSpec
{
   std::vector<double> x, y, z;
   ComponentCoeffs xx, yy, zz, xy, yz, xz;
   // Per-component bounds; default wide.
   double v_min = -1.0e12;
   double v_max =  1.0e12;
   // Overrides for intentionally-broken fixtures.
   bool   omit_sigma_xy = false;
   bool   shrink_sigma_yz_grid = false;     // axis-size mismatch
};

/// Build a synthetic schema-v1 sidecar carrying the six canonical
/// sigma_* components.
void write_stress_sidecar(const std::string& path,
                          const StressSidecarSpec& s)
{
   hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                          H5P_DEFAULT);
   write_string_attr(file, "schema_version", "data_projection_v1");
   write_string_attr(file, "crs",            "EPSG:32611");
   write_string_attr(file, "units",          "m");
   write_string_attr(file, "z_positive",     "elevation");
   write_string_attr(file, "created_at",     "1970-01-01T00:00:00Z");

   hid_t grid = H5Gcreate2(file, "/grid", H5P_DEFAULT, H5P_DEFAULT,
                           H5P_DEFAULT);
   write_axis(grid, "x", s.x);
   write_axis(grid, "y", s.y);
   write_axis(grid, "z", s.z);
   H5Gclose(grid);

   hid_t flds = H5Gcreate2(file, "/fields", H5P_DEFAULT, H5P_DEFAULT,
                           H5P_DEFAULT);

   const std::vector<double> data_xx =
      linear_field(s.x, s.y, s.z, s.xx.a, s.xx.b, s.xx.c, s.xx.d);
   const std::vector<double> data_yy =
      linear_field(s.x, s.y, s.z, s.yy.a, s.yy.b, s.yy.c, s.yy.d);
   const std::vector<double> data_zz =
      linear_field(s.x, s.y, s.z, s.zz.a, s.zz.b, s.zz.c, s.zz.d);
   const std::vector<double> data_xy =
      linear_field(s.x, s.y, s.z, s.xy.a, s.xy.b, s.xy.c, s.xy.d);

   write_field(flds, "sigma_xx", data_xx, s.x, s.y, s.z, s.v_min, s.v_max);
   write_field(flds, "sigma_yy", data_yy, s.x, s.y, s.z, s.v_min, s.v_max);
   write_field(flds, "sigma_zz", data_zz, s.x, s.y, s.z, s.v_min, s.v_max);
   if (!s.omit_sigma_xy)
   {
      write_field(flds, "sigma_xy", data_xy, s.x, s.y, s.z,
                  s.v_min, s.v_max);
   }
   // sigma_yz with optionally-shrunken grid for the consistency test.
   if (s.shrink_sigma_yz_grid)
   {
      std::vector<double> z_short(s.z.begin(), s.z.end() - 1);
      const std::vector<double> data_yz_short =
         linear_field(s.x, s.y, z_short, s.yz.a, s.yz.b, s.yz.c, s.yz.d);
      write_field(flds, "sigma_yz", data_yz_short, s.x, s.y, z_short,
                  s.v_min, s.v_max);
   }
   else
   {
      const std::vector<double> data_yz =
         linear_field(s.x, s.y, s.z, s.yz.a, s.yz.b, s.yz.c, s.yz.d);
      write_field(flds, "sigma_yz", data_yz, s.x, s.y, s.z,
                  s.v_min, s.v_max);
   }
   const std::vector<double> data_xz =
      linear_field(s.x, s.y, s.z, s.xz.a, s.xz.b, s.xz.c, s.xz.d);
   write_field(flds, "sigma_xz", data_xz, s.x, s.y, s.z,
               s.v_min, s.v_max);

   H5Gclose(flds);
   H5Fclose(file);
}

StressSidecarSpec make_default_spec()
{
   StressSidecarSpec s;
   s.x = {0.0, 1.0, 2.0};
   s.y = {0.0, 1.0, 2.0};
   s.z = {-2.0, -1.0, 0.0};
   // Six distinct linear functions so component swaps are caught.
   s.xx = {1.0, 0.0, 0.0, 100.0};   // varies with x
   s.yy = {0.0, 2.0, 0.0, 200.0};   // varies with y
   s.zz = {0.0, 0.0, 3.0, 300.0};   // varies with z
   s.xy = {0.0, 0.0, 0.0,  10.0};   // constant
   s.yz = {0.0, 0.0, 0.0,  20.0};   // constant
   s.xz = {0.0, 0.0, 0.0,  30.0};   // constant
   return s;
}

} // anonymous namespace


// ======================================================================
// Tests
// ======================================================================

static void T_6_1_construct_and_evaluate()
{
   std::cout << "\n[T-6-1] Construct StressField3D and evaluate at a grid corner\n";

   const auto spec = make_default_spec();
   const std::string path = make_tmp_path("ctor");
   write_stress_sidecar(path, spec);

   StressField3D field(path);

   const DenseMatrix S = field.Evaluate(1.0, 1.0, -1.0);

   // sigma_xx = 1*1 + 0*1 + 0*-1 + 100 = 101
   // sigma_yy = 0*1 + 2*1 + 0*-1 + 200 = 202
   // sigma_zz = 0*1 + 0*1 + 3*-1 + 300 = 297
   // sigma_xy = 10, sigma_yz = 20, sigma_xz = 30
   TEST_NEAR(S(0, 0), 101.0, 1e-9, "sigma_xx at (1, 1, -1)");
   TEST_NEAR(S(1, 1), 202.0, 1e-9, "sigma_yy at (1, 1, -1)");
   TEST_NEAR(S(2, 2), 297.0, 1e-9, "sigma_zz at (1, 1, -1)");
   TEST_NEAR(S(0, 1),  10.0, 1e-9, "sigma_xy at (1, 1, -1)");
   TEST_NEAR(S(1, 2),  20.0, 1e-9, "sigma_yz at (1, 1, -1)");
   TEST_NEAR(S(0, 2),  30.0, 1e-9, "sigma_xz at (1, 1, -1)");

   // Symmetry assertions.
   TEST_NEAR(S(1, 0), S(0, 1), 1e-12, "sigma_yx == sigma_xy (symmetry)");
   TEST_NEAR(S(2, 1), S(1, 2), 1e-12, "sigma_zy == sigma_yz (symmetry)");
   TEST_NEAR(S(2, 0), S(0, 2), 1e-12, "sigma_zx == sigma_xz (symmetry)");

   ::unlink(path.c_str());
}


static void T_6_2_field_accessor_and_names()
{
   std::cout << "\n[T-6-2] Field(component_index) and ComponentName lookups\n";

   const auto spec = make_default_spec();
   const std::string path = make_tmp_path("acc");
   write_stress_sidecar(path, spec);

   StressField3D field(path);

   TEST_ASSERT(field.NumComponents() == 6, "NumComponents == 6");
   TEST_ASSERT(field.Field(0).FieldName() == "sigma_xx", "Field(0) is sigma_xx");
   TEST_ASSERT(field.Field(1).FieldName() == "sigma_yy", "Field(1) is sigma_yy");
   TEST_ASSERT(field.Field(2).FieldName() == "sigma_zz", "Field(2) is sigma_zz");
   TEST_ASSERT(field.Field(3).FieldName() == "sigma_xy", "Field(3) is sigma_xy");
   TEST_ASSERT(field.Field(4).FieldName() == "sigma_yz", "Field(4) is sigma_yz");
   TEST_ASSERT(field.Field(5).FieldName() == "sigma_xz", "Field(5) is sigma_xz");

   TEST_ASSERT(StressField3D::ComponentName(0) == "sigma_xx",
               "ComponentName(0)");
   TEST_ASSERT(StressField3D::ComponentName(5) == "sigma_xz",
               "ComponentName(5)");

   ::unlink(path.c_str());
}


static void T_6_3_bbox_and_contains()
{
   std::cout << "\n[T-6-3] BBox accessor and ContainsBBox\n";

   const auto spec = make_default_spec();
   const std::string path = make_tmp_path("bbox");
   write_stress_sidecar(path, spec);

   StressField3D field(path);
   const auto& bb = field.BBox();
   TEST_NEAR(bb[0], 0.0, 1e-12, "bbox xmin");
   TEST_NEAR(bb[1], 2.0, 1e-12, "bbox xmax");
   TEST_NEAR(bb[2], 0.0, 1e-12, "bbox ymin");
   TEST_NEAR(bb[3], 2.0, 1e-12, "bbox ymax");
   TEST_NEAR(bb[4], -2.0, 1e-12, "bbox zmin");
   TEST_NEAR(bb[5], 0.0, 1e-12, "bbox zmax");

   TEST_ASSERT(field.ContainsBBox(0.1, 1.9, 0.1, 1.9, -1.9, -0.1),
               "ContainsBBox: interior bbox contained");
   TEST_ASSERT(!field.ContainsBBox(-0.1, 1.0, 0.0, 1.0, -1.0, 0.0),
               "ContainsBBox: bbox extending outside x rejected");

   ::unlink(path.c_str());
}


static void T_6_4_set_interp_mode_propagates()
{
   std::cout << "\n[T-6-4] SetInterpMode propagates to all six readers\n";

   const auto spec = make_default_spec();
   const std::string path = make_tmp_path("interp");
   write_stress_sidecar(path, spec);

   StressField3D field(path);
   TEST_ASSERT(field.GetInterpMode() == InterpMode::Trilinear,
               "default InterpMode is Trilinear");

   field.SetInterpMode(InterpMode::CatmullRom);
   TEST_ASSERT(field.GetInterpMode() == InterpMode::CatmullRom,
               "GetInterpMode returns CatmullRom after Set");
   for (int c = 0; c < 6; ++c)
   {
      TEST_ASSERT(field.Field(c).GetInterpMode() == InterpMode::CatmullRom,
                  std::string("component ")
                  + StressField3D::ComponentName(c)
                  + " has CatmullRom");
   }

   ::unlink(path.c_str());
}


static void T_6_5_pass_through_no_sign_flip()
{
   std::cout << "\n[T-6-5] Reader is pass-through: negative on-disk values "
             << "stay negative\n";

   StressSidecarSpec spec = make_default_spec();
   // Override the diagonal components to be NEGATIVE (which would be
   // H&Z compression-negative on disk if the writer hadn't done its
   // job).  The reader must NOT silently flip the sign.  This is the
   // R-501/R-502 contract guard.
   spec.xx = {0.0, 0.0, 0.0, -50.0};
   spec.yy = {0.0, 0.0, 0.0, -30.0};
   spec.zz = {0.0, 0.0, 0.0, -20.0};
   const std::string path = make_tmp_path("pass");
   write_stress_sidecar(path, spec);

   StressField3D field(path);
   const DenseMatrix S = field.Evaluate(1.0, 1.0, -1.0);
   TEST_NEAR(S(0, 0), -50.0, 1e-9, "sigma_xx pass-through (no flip)");
   TEST_NEAR(S(1, 1), -30.0, 1e-9, "sigma_yy pass-through (no flip)");
   TEST_NEAR(S(2, 2), -20.0, 1e-9, "sigma_zz pass-through (no flip)");

   ::unlink(path.c_str());
}


static void T_6_6_missing_component_aborts()
{
   std::cout << "\n[T-6-6] Missing sigma_xy in sidecar aborts ctor\n";

   StressSidecarSpec spec = make_default_spec();
   spec.omit_sigma_xy = true;
   const std::string path = make_tmp_path("missing");
   write_stress_sidecar(path, spec);

   const bool aborted = expect_abort([&]() { StressField3D f(path); });
   TEST_ASSERT(aborted,
               "StressField3D ctor aborts on missing sigma_xy field");

   ::unlink(path.c_str());
}


static void T_6_7_mismatched_grid_aborts()
{
   std::cout << "\n[T-6-7] Mismatched per-component grid aborts ctor\n";

   StressSidecarSpec spec = make_default_spec();
   spec.shrink_sigma_yz_grid = true;
   const std::string path = make_tmp_path("mismatch");
   write_stress_sidecar(path, spec);

   const bool aborted = expect_abort([&]() { StressField3D f(path); });
   TEST_ASSERT(aborted,
               "StressField3D ctor aborts on size-mismatched sigma_yz");

   ::unlink(path.c_str());
}


static void T_6_8_evaluate_returns_dense_matrix_3x3()
{
   std::cout << "\n[T-6-8] Evaluate returns a 3x3 DenseMatrix\n";

   const auto spec = make_default_spec();
   const std::string path = make_tmp_path("shape");
   write_stress_sidecar(path, spec);

   StressField3D field(path);
   const DenseMatrix S = field.Evaluate(0.5, 0.5, -0.5);
   TEST_ASSERT(S.Height() == 3, "DenseMatrix height == 3");
   TEST_ASSERT(S.Width() == 3, "DenseMatrix width == 3");
   // The off-diagonals come from the constant sigma_xy/yz/xz at this
   // point: 10 / 20 / 30.  Mirror equality is the symmetric-tensor
   // contract.
   TEST_NEAR(S(1, 0), 10.0, 1e-9, "sigma_yx mirrors sigma_xy");
   TEST_NEAR(S(2, 1), 20.0, 1e-9, "sigma_zy mirrors sigma_yz");
   TEST_NEAR(S(2, 0), 30.0, 1e-9, "sigma_zx mirrors sigma_xz");

   ::unlink(path.c_str());
}


static void T_6_9_oob_evaluate_aborts()
{
   std::cout << "\n[T-6-9] Out-of-bbox Evaluate aborts\n";

   const auto spec = make_default_spec();
   const std::string path = make_tmp_path("oob");
   write_stress_sidecar(path, spec);

   const bool aborted = expect_abort([&]() {
      StressField3D field(path);
      // x = 100 is outside x bbox [0, 2].
      (void) field.Evaluate(100.0, 1.0, -1.0);
   });
   TEST_ASSERT(aborted, "Evaluate aborts on out-of-bbox x");

   ::unlink(path.c_str());
}


static void T_6_10_field_index_out_of_range_aborts()
{
   std::cout << "\n[T-6-10] Field(idx) with idx outside [0, 6) aborts\n";

   const auto spec = make_default_spec();
   const std::string path = make_tmp_path("idx");
   write_stress_sidecar(path, spec);

   const bool aborted = expect_abort([&]() {
      StressField3D field(path);
      (void) field.Field(6);
   });
   TEST_ASSERT(aborted, "Field(6) aborts");

   ::unlink(path.c_str());
}


static void T_6_11_stress_field_coefficient_eval()
{
   // R-901: StressFieldCoefficient::Eval emits the six components
   // in schema-v1 canonical order (xx, yy, zz, xy, yz, xz) — NOT
   // Voigt — and applies scale/offset uniformly.  A silent reorder
   // bug in the Eval body would compile and pass every existing
   // T-6-* test, then poison the static-equilibration pipeline.
   std::cout << "\n[T-6-11] StressFieldCoefficient::Eval order + "
             << "scale/offset (R-901)\n";

   const auto spec = make_default_spec();   // bbox [0,2] x [0,2] x [-2,0]
   const std::string path = make_tmp_path("scfeval");
   write_stress_sidecar(path, spec);

   StressField3D field(path);
   StressFieldCoefficient coef(field);

   // A unit-hex serial mesh, located strictly inside the synthetic
   // sidecar bbox so the trilinear interpolation never aborts.  The
   // 8 vertices span [0.5, 1.5]^2 x [-1.5, -0.5].
   mfem::Mesh mesh(/*dim=*/3, /*NVert=*/8, /*NElem=*/1);
   const double verts[8][3] = {
      {0.5, 0.5, -1.5}, {1.5, 0.5, -1.5},
      {1.5, 1.5, -1.5}, {0.5, 1.5, -1.5},
      {0.5, 0.5, -0.5}, {1.5, 0.5, -0.5},
      {1.5, 1.5, -0.5}, {0.5, 1.5, -0.5},
   };
   for (int v = 0; v < 8; ++v)
   {
      mesh.AddVertex(verts[v]);
   }
   const int hex_idx[8] = {0, 1, 2, 3, 4, 5, 6, 7};
   mesh.AddHex(hex_idx);
   mesh.FinalizeHexMesh(/*generate_edges=*/1);

   mfem::ElementTransformation& T = *mesh.GetElementTransformation(0);
   mfem::IntegrationPoint ip;
   ip.Set3(0.5, 0.5, 0.5);   // element-reference center

   mfem::Vector x(3);
   T.Transform(ip, x);
   mfem::Vector v(6);
   coef.Eval(v, T, ip);

   const real_t xc = x[0], yc = x[1], zc = x[2];
   const real_t xx_ref = field.Field(0).Evaluate(xc, yc, zc);
   const real_t yy_ref = field.Field(1).Evaluate(xc, yc, zc);
   const real_t zz_ref = field.Field(2).Evaluate(xc, yc, zc);
   const real_t xy_ref = field.Field(3).Evaluate(xc, yc, zc);
   const real_t yz_ref = field.Field(4).Evaluate(xc, yc, zc);
   const real_t xz_ref = field.Field(5).Evaluate(xc, yc, zc);

   TEST_ASSERT(coef.GetVDim() == 6,
               "StressFieldCoefficient::GetVDim() == 6");
   TEST_NEAR(v(0), xx_ref, 1e-9,
             "v[0] == sigma_xx (schema-v1 canonical)");
   TEST_NEAR(v(1), yy_ref, 1e-9,
             "v[1] == sigma_yy");
   TEST_NEAR(v(2), zz_ref, 1e-9,
             "v[2] == sigma_zz");
   TEST_NEAR(v(3), xy_ref, 1e-9,
             "v[3] == sigma_xy (NOT Voigt order!)");
   TEST_NEAR(v(4), yz_ref, 1e-9,
             "v[4] == sigma_yz");
   TEST_NEAR(v(5), xz_ref, 1e-9,
             "v[5] == sigma_xz");

   // Scale / offset apply uniformly.
   StressFieldCoefficient coef2(field, /*scale=*/2.0, /*offset=*/100.0);
   mfem::Vector v2(6);
   coef2.Eval(v2, T, ip);
   TEST_NEAR(v2(0), 2.0 * xx_ref + 100.0, 1e-9,
             "scale * sigma_xx + offset");
   TEST_NEAR(v2(3), 2.0 * xy_ref + 100.0, 1e-9,
             "scale * sigma_xy + offset");
   TEST_NEAR(v2(5), 2.0 * xz_ref + 100.0, 1e-9,
             "scale * sigma_xz + offset");

   ::unlink(path.c_str());
}


// ======================================================================
// Main
// ======================================================================

int main(int argc, char** argv)
{
   (void) argc;
   (void) argv;
   std::cout << "Running StressField3D unit tests "
             << "(Phase 6 §1 of PLAN_onfaultstress.md)\n";

   T_6_1_construct_and_evaluate();
   T_6_2_field_accessor_and_names();
   T_6_3_bbox_and_contains();
   T_6_4_set_interp_mode_propagates();
   T_6_5_pass_through_no_sign_flip();
   T_6_6_missing_component_aborts();
   T_6_7_mismatched_grid_aborts();
   T_6_8_evaluate_returns_dense_matrix_3x3();
   T_6_9_oob_evaluate_aborts();
   T_6_10_field_index_out_of_range_aborts();
   T_6_11_stress_field_coefficient_eval();

   std::cout << "\n=========================================\n";
   std::cout << "Total: " << num_tests << "\n";
   std::cout << "Passed: " << num_passed << "\n";
   std::cout << "Failed: " << num_failed << "\n";
   std::cout << "=========================================\n";
   return (num_failed == 0) ? 0 : 1;
}

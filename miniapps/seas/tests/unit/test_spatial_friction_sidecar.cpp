// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_spatial_friction_sidecar.cpp — Phase 3 of
// PLAN_thermal_case2_mixedflux_port_2026-07-08.md.
//
// Covers `[friction.rate_state.sidecar]`: per-DOF `a` and `V_w` (and optionally
// `b` / `Dc`) sampled from a 3-D `data_projection_v1` HDF5 sidecar, the channel
// that carries the SCEC Community Thermal Model zoning into the resolver.
//
// Three plan-mandated tests:
//   S-1 round-trip : synthetic 3x3x3 sidecar with a known linear a(x,y,z);
//                    resolve at the 8 corners + the centre; trilinear to 1e-12.
//   S-2 clamp      : out-of-hull query with far_field_clamp=true -> edge value;
//                    with false -> abort.
//   S-3 precedence : sidecar + a `box` rule setting V_w -> the rule wins inside
//                    the box, the sidecar wins outside.
//
// Plus the guards the plan's design implies:
//   S-4 legacy     : a null sidecar leaves the resolver byte-identical.
//   S-5 optional   : b / Dc evaluators override only when supplied.
//   S-6 ctor guard : a half-built RateStateSidecarFields aborts.
//
// The resolver consumes `RateStateSidecarFields::FieldFn` (std::function), not
// DataField3D directly (deviation D-A, see spatial_friction.hpp).  These tests
// therefore build the callables exactly the way the driver's
// MakeRateStateSidecarFields does: a shared_ptr<DataField3D> captured by value.

#include "mfem.hpp"

#include "../../spatial/code/spatial_friction.hpp"
#include "../../dynamic/heterogeneous_material.hpp"
#include "../../io/data_field_3d.hpp"

#include <hdf5.h>

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

using namespace mfem;
using namespace mfem::seas;
using namespace mfem::seas::spatial;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(c, m) do { num_tests++; if (!(c)) { \
   std::cerr << "FAILED: " << m << " line " << __LINE__ << "\n"; \
   num_failed++; } else { std::cout << "  PASSED: " << m << "\n"; \
   num_passed++; } } while (0)

#define TEST_NEAR(v, e, t, m) do { num_tests++; const real_t _v=(v); \
   const real_t _e=(e); const real_t _t=(t); \
   if (std::abs(_v - _e) > _t) { \
      std::cerr << "FAILED: " << m << " (got " << _v << ", expected " << _e \
                << ", tol " << _t << ") line " << __LINE__ << "\n"; \
      num_failed++; } else { std::cout << "  PASSED: " << m << "\n"; \
      num_passed++; } } while (0)

namespace
{

// ----------------------------------------------------------------------
// schema-v1 HDF5 fixture writer (same layout as tests/unit/test_data_field_3d.cpp)
// ----------------------------------------------------------------------

void write_string_attr(hid_t loc, const char* name, const std::string& v)
{
   hid_t aspace = H5Screate(H5S_SCALAR);
   hid_t atype  = H5Tcopy(H5T_C_S1);
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

struct FieldSpec
{
   std::string name;
   std::string units;
   double      min_value;
   double      max_value;
   std::vector<double> data;   ///< (Nx, Ny, Nz) row-major
};

/// Write a multi-field schema-v1 sidecar (the friction sidecar carries at
/// least rs_a + rs_srW, so the single-field helper in test_data_field_3d.cpp
/// does not suffice).
void write_sidecar(const std::string& path,
                   const std::vector<double>& x,
                   const std::vector<double>& y,
                   const std::vector<double>& z,
                   const std::vector<FieldSpec>& fields)
{
   hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   write_string_attr(file, "schema_version", "data_projection_v1");
   write_string_attr(file, "crs",            "EPSG:32611");
   write_string_attr(file, "units",          "m");
   write_string_attr(file, "z_positive",     "elevation");
   write_string_attr(file, "created_at",     "1970-01-01T00:00:00Z");

   hid_t grid = H5Gcreate2(file, "/grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   write_axis(grid, "x", x);
   write_axis(grid, "y", y);
   write_axis(grid, "z", z);
   H5Gclose(grid);

   hid_t flds = H5Gcreate2(file, "/fields", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   hsize_t shape[3] = { x.size(), y.size(), z.size() };
   for (const auto& f : fields)
   {
      hid_t sid = H5Screate_simple(3, shape, nullptr);
      hid_t did = H5Dcreate2(flds, f.name.c_str(), H5T_NATIVE_DOUBLE, sid,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
               f.data.data());
      write_string_attr(did, "units", f.units);
      write_double_attr(did, "min_value", f.min_value);
      write_double_attr(did, "max_value", f.max_value);
      H5Dclose(did);
      H5Sclose(sid);
   }
   H5Gclose(flds);
   H5Fclose(file);
}

std::string tmp_path(const char* stem)
{
   const char* tmpdir = std::getenv("TMPDIR");
   if (!tmpdir || !*tmpdir) { tmpdir = "/tmp"; }
   std::string p = std::string(tmpdir) + "/seas_test_spatial_friction_sidecar_"
                   + stem + "_" + std::to_string(::getpid()) + ".h5";
   ::unlink(p.c_str());
   return p;
}

// The fixture grid.  Deliberately NOT centred on the origin so a sign error in
// the axis handling shows up.
const std::vector<double> kX = { 0.0, 1.0, 2.0 };
const std::vector<double> kY = { -1.0, 0.0, 1.0 };
const std::vector<double> kZ = { -2.0, -1.0, 0.0 };

/// f(x,y,z) = c0 + cx*x + cy*y + cz*z, sampled on the fixture grid.
/// Trilinear interpolation reproduces a linear function EXACTLY, so any
/// deviation is an indexing / axis-order bug, not interpolation error.
std::vector<double> linear_field(double c0, double cx, double cy, double cz)
{
   std::vector<double> d(kX.size() * kY.size() * kZ.size(), 0.0);
   for (size_t i = 0; i < kX.size(); ++i)
      for (size_t j = 0; j < kY.size(); ++j)
         for (size_t k = 0; k < kZ.size(); ++k)
         {
            const size_t flat = (i * kY.size() + j) * kZ.size() + k;
            d[flat] = c0 + cx * kX[i] + cy * kY[j] + cz * kZ[k];
         }
   return d;
}

real_t linear_eval(double c0, double cx, double cy, double cz,
                   real_t x, real_t y, real_t z)
{
   return c0 + cx * x + cy * y + cz * z;
}

// a(x,y,z) coefficients: a in [0.010, 0.030] over the fixture grid.
constexpr double kA0 = 0.020, kAx = 0.002, kAy = 0.001, kAz = 0.0005;
// V_w(x,y,z) coefficients: V_w in [0.05, ~5] over the fixture grid.
constexpr double kV0 = 2.0,  kVx = 0.5,   kVy = 0.25,  kVz = 0.125;

std::string make_friction_sidecar()
{
   const std::string path = tmp_path("fields");
   std::vector<FieldSpec> f;
   f.push_back({ "rs_a",   "1",   -1.0e3, 1.0e3, linear_field(kA0, kAx, kAy, kAz) });
   f.push_back({ "rs_srW", "m/s", -1.0e3, 1.0e3, linear_field(kV0, kVx, kVy, kVz) });
   f.push_back({ "rs_b",   "1",   -1.0e3, 1.0e3, linear_field(0.030, 0.0, 0.0, 0.0) });
   f.push_back({ "rs_sl0", "m",   -1.0e3, 1.0e3, linear_field(0.250, 0.0, 0.0, 0.0) });
   write_sidecar(path, kX, kY, kZ, f);
   return path;
}

/// Build the resolver's callables the same way the driver's
/// MakeRateStateSidecarFields does.
RateStateSidecarFields::FieldFn bind_field(const std::string& path,
                                           const std::string& field,
                                           OOBPolicy oob)
{
   auto reader = std::make_shared<DataField3D>(path, field, oob);
   return [reader](real_t x, real_t y, real_t z)
          { return reader->Evaluate(x, y, z); };
}

/// Mirrors the driver exactly, INCLUDING the R-001 bbox / clamp-flag capture:
/// a test helper that skipped those would not exercise the code path the
/// diagnostic depends on.
std::shared_ptr<RateStateSidecarFields>
make_fields(const std::string& path, OOBPolicy oob,
            bool with_b = false, bool with_Dc = false)
{
   auto f = std::make_shared<RateStateSidecarFields>();

   auto a_reader = std::make_shared<DataField3D>(path, "rs_a", oob);
   f->bbox               = a_reader->BBox();
   f->clamps_out_of_hull = (oob == OOBPolicy::Clamp);
   f->a = [a_reader](real_t x, real_t y, real_t z)
          { return a_reader->Evaluate(x, y, z); };

   f->V_w = bind_field(path, "rs_srW", oob);
   if (with_b)  { f->b  = bind_field(path, "rs_b",   oob); }
   if (with_Dc) { f->Dc = bind_field(path, "rs_sl0", oob); }
   return f;
}

class TinyMeshHolder
{
public:
   TinyMeshHolder() : mesh_(mfem::Mesh::MakeCartesian3D(
                              1, 1, 1, mfem::Element::HEXAHEDRON,
                              1.0, 1.0, 1.0)) { }
   mfem::Mesh& mesh() { return mesh_; }
private:
   mfem::Mesh mesh_;
};

/// Run `body` in a forked child.  True iff the child aborted.
bool RunInChild(const std::function<void()>& body)
{
   ::fflush(stdout);
   ::fflush(stderr);
   const pid_t pid = ::fork();
   if (pid < 0) { std::cerr << "fork() failed\n"; return false; }
   if (pid == 0)
   {
      ::close(1);
      ::close(2);
      body();
      std::_Exit(0);
   }
   int status = 0;
   while (::waitpid(pid, &status, 0) < 0 && errno == EINTR) { /* retry */ }
   if (WIFEXITED(status))   { return WEXITSTATUS(status) != 0; }
   if (WIFSIGNALED(status)) { return true; }
   return true;
}

/// Minimal SRW rate-state block.  a_default / V_w_default are deliberately
/// NOT the sidecar values, so a test that reads them instead of the sidecar
/// fails loudly.
RateStateBlock make_cfg()
{
   RateStateBlock cfg;
   cfg.a_default       = 0.0111;   // sentinel: never expected once sidecar wins
   cfg.b_default       = 0.019;
   cfg.Dc_default      = 0.10;
   cfg.V_init_default  = 1.0e-12;
   cfg.f_0_default     = 0.6;
   cfg.V_0_default     = 1.0e-6;
   cfg.sigma_n_default = 50.0e6;
   cfg.eta_auto        = false;
   cfg.eta_default     = 5.0e6;
   cfg.state_evolution = StateEvolutionKind::SlipLawStrongRateWeakening;
   cfg.f_w_default     = 0.0;      // the SAFS v3_4_1 value (Phase 2)
   cfg.V_w_default     = 0.7777;   // sentinel
   return cfg;
}

/// Resolve at the supplied coordinates.  `pts` is a flat (3N) coordinate list.
RateStatePerDOFParams resolve_at(const std::vector<real_t>& pts,
                                 const RateStateBlock& cfg,
                                 std::shared_ptr<const RateStateSidecarFields> f)
{
   const int N = static_cast<int>(pts.size() / 3);
   Vector dofs(3 * N);
   for (int i = 0; i < 3 * N; ++i) { dofs(i) = pts[i]; }
   Array<int> attr(N), elem(N);
   for (int i = 0; i < N; ++i) { attr[i] = 101; elem[i] = 0; }
   Vector sn_total(N);  sn_total = 50.0e6;

   auto mat = MaterialField::MakeConstant(32.0e9, 32.0e9, 2670.0);
   TinyMeshHolder mh;
   PorePressureSpec pp;   // all zero
   SpatialFrictionResolver R(std::move(f));
   mfem::Mesh& srl = mh.mesh();
   return R.ResolveRateState(cfg, dofs, elem, attr, mat, srl, pp, sn_total);
}

}  // namespace

// ==========================================================================
// S-1  Round-trip: trilinear reproduction of a linear field, to 1e-12
// ==========================================================================
static void S_1_roundtrip_trilinear()
{
   std::cout << "\n[S-1] sidecar round-trip: 8 corners + centre, trilinear exact\n";
   const std::string path = make_friction_sidecar();
   auto fields = make_fields(path, OOBPolicy::Abort);

   // 8 grid corners + the geometric centre (a genuine interpolation point).
   std::vector<real_t> pts;
   for (double x : { kX.front(), kX.back() })
      for (double y : { kY.front(), kY.back() })
         for (double z : { kZ.front(), kZ.back() })
         { pts.push_back(x); pts.push_back(y); pts.push_back(z); }
   pts.push_back(1.0); pts.push_back(0.0); pts.push_back(-1.0);   // centre
   pts.push_back(0.5); pts.push_back(-0.5); pts.push_back(-1.5);  // off-node

   const auto cfg = make_cfg();
   const auto p   = resolve_at(pts, cfg, fields);

   const int N = static_cast<int>(pts.size() / 3);
   int a_ok = 0, v_ok = 0;
   for (int i = 0; i < N; ++i)
   {
      const real_t x = pts[3*i], y = pts[3*i+1], z = pts[3*i+2];
      if (std::abs(p.a(i)   - linear_eval(kA0, kAx, kAy, kAz, x, y, z)) <= 1e-12)
      { ++a_ok; }
      if (std::abs(p.V_w(i) - linear_eval(kV0, kVx, kVy, kVz, x, y, z)) <= 1e-12)
      { ++v_ok; }
   }
   TEST_ASSERT(a_ok == N,  "a(x,y,z) trilinear-exact at all 10 probes");
   TEST_ASSERT(v_ok == N,  "V_w(x,y,z) trilinear-exact at all 10 probes");

   // The scalar defaults must NOT leak through when the sidecar is present.
   TEST_ASSERT(std::abs(p.a(0)   - cfg.a_default)   > 1e-6,
               "a is NOT the scalar a_default (sidecar wins)");
   TEST_ASSERT(std::abs(p.V_w(0) - cfg.V_w_default) > 1e-6,
               "V_w is NOT the scalar V_w_default (sidecar wins)");
   // b / Dc were not supplied -> the resolved defaults survive.
   TEST_NEAR(p.b(0),  cfg.b_default,  0.0, "b keeps b_default (no b_field)");
   TEST_NEAR(p.Dc(0), cfg.Dc_default, 0.0, "Dc keeps Dc_default (no Dc_field)");

   ::unlink(path.c_str());
}

// ==========================================================================
// S-2  Clamp: out-of-hull query holds the edge value; Abort policy aborts
// ==========================================================================
static void S_2_far_field_clamp()
{
   std::cout << "\n[S-2] out-of-hull: clamp holds the edge value, abort aborts\n";
   const std::string path = make_friction_sidecar();

   // A point well outside the hull on every axis (hull is x[0,2] y[-1,1] z[-2,0]).
   const std::vector<real_t> out = { 10.0, 10.0, 10.0 };
   // Clamped coordinates: the nearest hull corner.
   const real_t cx = kX.back(), cy = kY.back(), cz = kZ.back();

   {
      auto fields = make_fields(path, OOBPolicy::Clamp);
      const auto p = resolve_at(out, make_cfg(), fields);
      TEST_NEAR(p.a(0),   linear_eval(kA0, kAx, kAy, kAz, cx, cy, cz), 1e-12,
                "clamp: a = value at the nearest hull corner");
      TEST_NEAR(p.V_w(0), linear_eval(kV0, kVx, kVy, kVz, cx, cy, cz), 1e-12,
                "clamp: V_w = value at the nearest hull corner");
   }
   {
      const bool aborted = RunInChild([&]()
      {
         auto fields = make_fields(path, OOBPolicy::Abort);
         (void)resolve_at(out, make_cfg(), fields);
      });
      TEST_ASSERT(aborted, "abort policy: out-of-hull query aborts");
   }
   // In-hull queries must be unaffected by the policy (bit-identical).
   {
      const std::vector<real_t> in = { 0.5, -0.5, -1.5 };
      const auto pc = resolve_at(in, make_cfg(), make_fields(path, OOBPolicy::Clamp));
      const auto pa = resolve_at(in, make_cfg(), make_fields(path, OOBPolicy::Abort));
      TEST_ASSERT(pc.a(0) == pa.a(0) && pc.V_w(0) == pa.V_w(0),
                  "in-hull query bit-identical under clamp vs abort");
   }
   ::unlink(path.c_str());
}

// ==========================================================================
// S-3  Precedence: an explicit `box` rule beats the sidecar inside the box
// ==========================================================================
static void S_3_rule_beats_sidecar_inside_box()
{
   std::cout << "\n[S-3] precedence: box rule wins inside, sidecar wins outside\n";
   const std::string path = make_friction_sidecar();
   auto fields = make_fields(path, OOBPolicy::Abort);

   auto cfg = make_cfg();
   SpatialRule r;
   r.kind    = SpatialRule::Kind::Box;
   r.x_min_m = 0.9;  r.x_max_m = 1.1;
   r.y_min_m = -0.1; r.y_max_m = 0.1;
   r.z_min_m = -1.1; r.z_max_m = -0.9;
   r.V_w     = 0.05;      // the SAFS seismogenic V_w
   cfg.spatial.push_back(r);

   // p0 inside the box, p1 outside it (same z-plane).
   const std::vector<real_t> pts = { 1.0, 0.0, -1.0,
                                     2.0, 1.0,  0.0 };
   const auto p = resolve_at(pts, cfg, fields);

   TEST_NEAR(p.V_w(0), 0.05, 0.0, "V_w inside box = rule value (rule wins)");
   TEST_NEAR(p.V_w(1), linear_eval(kV0, kVx, kVy, kVz, 2.0, 1.0, 0.0), 1e-12,
             "V_w outside box = sidecar value");
   // `a` was not overridden by the rule -> the sidecar governs at BOTH points.
   TEST_NEAR(p.a(0), linear_eval(kA0, kAx, kAy, kAz, 1.0, 0.0, -1.0), 1e-12,
             "a inside box still from sidecar (rule sets only V_w)");
   TEST_NEAR(p.a(1), linear_eval(kA0, kAx, kAy, kAz, 2.0, 1.0, 0.0), 1e-12,
             "a outside box from sidecar");

   ::unlink(path.c_str());
}

// ==========================================================================
// S-4  Legacy: a null sidecar leaves the resolver on the scalar path
// ==========================================================================
static void S_4_null_sidecar_is_legacy()
{
   std::cout << "\n[S-4] null sidecar -> scalar defaults (legacy path unchanged)\n";
   const std::vector<real_t> pts = { 1.0, 0.0, -1.0 };
   const auto cfg = make_cfg();

   const auto p_null    = resolve_at(pts, cfg, nullptr);
   // The default-constructed resolver must agree bit-for-bit with the
   // explicitly-null one (the shared_ptr overload is the only new entry point).
   Vector dofs(3); dofs(0) = 1.0; dofs(1) = 0.0; dofs(2) = -1.0;
   Array<int> attr(1), elem(1); attr[0] = 101; elem[0] = 0;
   Vector sn_total(1); sn_total = 50.0e6;
   auto mat = MaterialField::MakeConstant(32.0e9, 32.0e9, 2670.0);
   TinyMeshHolder mh;
   PorePressureSpec pp;
   SpatialFrictionResolver R_default;
   const auto p_default = R_default.ResolveRateState(cfg, dofs, elem, attr, mat,
                                                     mh.mesh(), pp, sn_total);

   TEST_NEAR(p_null.a(0),   cfg.a_default,   0.0, "null sidecar -> a_default");
   TEST_NEAR(p_null.V_w(0), cfg.V_w_default, 0.0, "null sidecar -> V_w_default");
   TEST_ASSERT(p_null.a(0)   == p_default.a(0) &&
               p_null.V_w(0) == p_default.V_w(0) &&
               p_null.b(0)   == p_default.b(0)   &&
               p_null.Dc(0)  == p_default.Dc(0),
               "explicit-null ctor == default ctor, bit-for-bit");
}

// ==========================================================================
// S-5  Optional b / Dc evaluators override only when supplied
// ==========================================================================
static void S_5_optional_b_and_Dc()
{
   std::cout << "\n[S-5] optional b / Dc fields override the scalar defaults\n";
   const std::string path = make_friction_sidecar();
   const std::vector<real_t> pts = { 1.0, 0.0, -1.0 };
   const auto cfg = make_cfg();

   auto fields = make_fields(path, OOBPolicy::Abort, /*with_b=*/true,
                             /*with_Dc=*/true);
   const auto p = resolve_at(pts, cfg, fields);
   TEST_NEAR(p.b(0),  0.030, 1e-12, "b  from the rs_b sidecar field");
   TEST_NEAR(p.Dc(0), 0.250, 1e-12, "Dc from the rs_sl0 sidecar field");

   // Supplying only b leaves Dc on the default.
   auto only_b = make_fields(path, OOBPolicy::Abort, true, false);
   const auto q = resolve_at(pts, cfg, only_b);
   TEST_NEAR(q.b(0),  0.030,          1e-12, "b overridden when b_field given");
   TEST_NEAR(q.Dc(0), cfg.Dc_default, 0.0,   "Dc keeps default when Dc_field absent");

   ::unlink(path.c_str());
}

// ==========================================================================
// S-6  Constructor guard: a half-built sidecar aborts
// ==========================================================================
static void S_6_ctor_rejects_incomplete_sidecar()
{
   std::cout << "\n[S-6] resolver ctor rejects a sidecar missing a or V_w\n";
   const std::string path = make_friction_sidecar();

   const bool no_a = RunInChild([&]()
   {
      auto f = std::make_shared<RateStateSidecarFields>();
      f->V_w = bind_field(path, "rs_srW", OOBPolicy::Abort);   // a left empty
      SpatialFrictionResolver R(f);
      (void)R;
   });
   TEST_ASSERT(no_a, "sidecar without an `a` evaluator aborts in the ctor");

   const bool no_vw = RunInChild([&]()
   {
      auto f = std::make_shared<RateStateSidecarFields>();
      f->a = bind_field(path, "rs_a", OOBPolicy::Abort);       // V_w left empty
      SpatialFrictionResolver R(f);
      (void)R;
   });
   TEST_ASSERT(no_vw, "sidecar without a `V_w` evaluator aborts in the ctor");

   ::unlink(path.c_str());
}

// ==========================================================================
// S-7  Resolver rejects sidecar + depth_profile (the parser guard's twin)
// ==========================================================================
static void S_7_resolver_rejects_sidecar_plus_depth_profile()
{
   std::cout << "\n[S-7] resolver rejects sidecar + depth_profile together\n";
   const std::string path = make_friction_sidecar();

   const bool aborted = RunInChild([&]()
   {
      auto cfg = make_cfg();
      // Hand-build the ambiguous pair, bypassing the parser entirely.
      cfg.depth_profile.enabled = true;
      cfg.depth_profile.profile.a_of_depth.x   = { 0.0, 20000.0 };
      cfg.depth_profile.profile.a_of_depth.y   = { 0.010, 0.020 };
      cfg.depth_profile.profile.amb_of_depth.x = { 0.0, 20000.0 };
      cfg.depth_profile.profile.amb_of_depth.y = { -0.004, 0.010 };
      const std::vector<real_t> pts = { 1.0, 0.0, -1.0 };
      (void)resolve_at(pts, cfg, make_fields(path, OOBPolicy::Abort));
   });
   TEST_ASSERT(aborted, "ResolveRateState aborts on sidecar + depth_profile");

   // Control: the SAME depth_profile with NO sidecar must resolve fine, so the
   // abort above is attributable to the pair, not to the hand-built profile.
   auto cfg = make_cfg();
   cfg.depth_profile.enabled = true;
   cfg.depth_profile.profile.a_of_depth.x   = { 0.0, 20000.0 };
   cfg.depth_profile.profile.a_of_depth.y   = { 0.010, 0.020 };
   cfg.depth_profile.profile.amb_of_depth.x = { 0.0, 20000.0 };
   cfg.depth_profile.profile.amb_of_depth.y = { -0.004, 0.010 };
   const std::vector<real_t> pts = { 1.0, 0.0, -1.0 };
   const auto p = resolve_at(pts, cfg, nullptr);
   TEST_ASSERT(p.a(0) > 0.0 && p.b(0) > 0.0,
               "depth_profile alone (no sidecar) resolves");

   ::unlink(path.c_str());
}

// ==========================================================================
// S-8  R-001: bbox is captured, and InHull() distinguishes a real sample from
//             a clamped edge value.
//
// The regression: under OOBPolicy::Clamp, an out-of-hull read returns the
// nearest-edge value with no way to tell it apart from a genuine sample.  For
// the production CTM CASE2 field the edge values ARE the hypocentre anchors
// (a = 0.015, V_w = 0.05), so `--print-derived` printed a perfect "pass" for a
// point 300 km outside the hull.  InHull() is what makes the gate falsifiable.
// ==========================================================================
static void S_8_bbox_and_in_hull_gate_clamped_reads()
{
   std::cout << "\n[S-8] bbox captured; InHull() separates real samples from clamped edges\n";
   const std::string path = make_friction_sidecar();
   auto f = make_fields(path, OOBPolicy::Clamp);

   // The fixture hull is x[0,2] y[-1,1] z[-2,0].
   TEST_NEAR(f->bbox[0], kX.front(), 0.0, "bbox xmin captured from DataField3D");
   TEST_NEAR(f->bbox[1], kX.back(),  0.0, "bbox xmax captured");
   TEST_NEAR(f->bbox[2], kY.front(), 0.0, "bbox ymin captured");
   TEST_NEAR(f->bbox[3], kY.back(),  0.0, "bbox ymax captured");
   TEST_NEAR(f->bbox[4], kZ.front(), 0.0, "bbox zmin captured");
   TEST_NEAR(f->bbox[5], kZ.back(),  0.0, "bbox zmax captured");

   TEST_ASSERT(f->InHull(1.0, 0.0, -1.0),   "hull centre is in-hull");
   TEST_ASSERT(f->InHull(kX.front(), kY.front(), kZ.front()),
               "min corner is in-hull (inclusive)");
   TEST_ASSERT(f->InHull(kX.back(), kY.back(), kZ.back()),
               "max corner is in-hull (inclusive)");
   TEST_ASSERT(!f->InHull(10.0, 10.0, 10.0), "far point is out of hull");
   TEST_ASSERT(!f->InHull(kX.back() + 1e-6, 0.0, -1.0),
               "just past x_max is out of hull");

   // THE BUG: a clamped out-of-hull read is byte-identical to the edge sample.
   // Without InHull() there is no way for a diagnostic to notice.
   const real_t clamped   = f->a(10.0, 10.0, 10.0);
   const real_t at_corner = f->a(kX.back(), kY.back(), kZ.back());
   TEST_ASSERT(clamped == at_corner,
               "clamped out-of-hull read == edge value (this is WHY InHull must gate it)");

   TEST_ASSERT(f->clamps_out_of_hull, "clamp flag recorded (Clamp policy)");
   auto g = make_fields(path, OOBPolicy::Abort);
   TEST_ASSERT(!g->clamps_out_of_hull, "clamp flag recorded (Abort policy)");

   // A default-constructed struct must NOT claim points are in-hull: a caller
   // that forgets to set bbox gets "outside", never a spurious pass.
   RateStateSidecarFields bare;
   TEST_ASSERT(!bare.InHull(1.0, 0.0, -1.0),
               "default (all-zero) bbox rejects a real point (conservative)");

   ::unlink(path.c_str());
}

int main(int argc, char* argv[])
{
   // Serial-only, like tests/unit/test_spatial_friction_resolver.cpp and
   // test_data_field_3d.cpp: the resolver's serial-mesh overload needs no MPI,
   // and fork()-based abort tests are cleaner without an initialized runtime.
   (void)argc; (void)argv;

   std::cout << "Running test_spatial_friction_sidecar "
                "(PLAN_thermal_case2 Phase 3)\n";

   S_1_roundtrip_trilinear();
   S_2_far_field_clamp();
   S_3_rule_beats_sidecar_inside_box();
   S_4_null_sidecar_is_legacy();
   S_5_optional_b_and_Dc();
   S_6_ctor_rejects_incomplete_sidecar();
   S_7_resolver_rejects_sidecar_plus_depth_profile();
   S_8_bbox_and_in_hull_gate_clamped_reads();

   std::cout << "\n========================================\n";
   std::cout << "test_spatial_friction_sidecar: " << num_passed << " / "
             << num_tests << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}

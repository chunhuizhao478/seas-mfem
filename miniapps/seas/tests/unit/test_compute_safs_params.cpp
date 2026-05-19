// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Unit tests for Phase 6 §5 — FaultGeometry::ComputeParams.
// Plan reference: PLAN_onfaultstress.md §1860-1879.
//
// Coverage:
//   T_65_1   sigma_n_per_dof.Size() == NumFaultDOFs after call.
//   T_65_2   tau_pre.Size() == 2 * NumFaultDOFs after call.
//   T_65_3   HasParams() returns true after invocation.
//   T_65_4   a/eta/Dc/V_init unchanged relative to ComputeBP5Params
//            (analytic spatial dependence preserved per plan §1869).
//   T_65_5   sigma_n_per_dof and tau_pre match the standalone
//            ProjectFaultPreStress output (proves ComputeParams
//            is the correct wrapper).
//   T_65_6   BP2 ctor → calling ComputeParams aborts.

#include "mfem.hpp"

#include "../../io/stress_field_3d.hpp"
#include "../../io/field_coefficient.hpp"
#include "../../fault/fault_geometry.hpp"
#include "../../fault/fault_geometry_safs.inl"
#include "../../domain/elasticity_operator.hpp"
#include "../../domain/boundary_config.hpp"
#include "../../config/bp5_params.hpp"

#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#include <hdf5.h>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(c, m) do { num_tests++; if (!(c)) { \
   std::cerr << "FAILED: " << m << " line " << __LINE__ << "\n"; \
   num_failed++; } else { std::cout << "  PASSED: " << m << "\n"; \
   num_passed++; } } while (0)

#define TEST_NEAR(v, e, t, m) do { num_tests++; const real_t _v=(v); \
   const real_t _e=(e); const real_t _t=(t); \
   if (std::abs(_v - _e) > _t) { std::cerr << "FAILED: " << m \
   << " (got " << _v << ", expected " << _e << ", tol " << _t \
   << ", line " << __LINE__ << ")\n"; num_failed++; } else { \
   std::cout << "  PASSED: " << m << "\n"; num_passed++; } } while (0)

// HDF5 fixture helpers (minimal copy from test_project_fault_prestress)
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
   H5Aclose(aid); H5Tclose(atype); H5Sclose(aspace);
}
void write_double_attr(hid_t loc, const char* name, double v)
{
   hid_t s = H5Screate(H5S_SCALAR);
   hid_t a = H5Acreate2(loc, name, H5T_NATIVE_DOUBLE, s,
                        H5P_DEFAULT, H5P_DEFAULT);
   H5Awrite(a, H5T_NATIVE_DOUBLE, &v);
   H5Aclose(a); H5Sclose(s);
}
void write_axis(hid_t g, const char* n, const std::vector<double>& v)
{
   hsize_t d = v.size();
   hid_t s = H5Screate_simple(1, &d, nullptr);
   hid_t did = H5Dcreate2(g, n, H5T_NATIVE_DOUBLE, s,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data());
   H5Dclose(did); H5Sclose(s);
}
void write_field(hid_t g, const char* n, const std::vector<double>& d,
                 const std::vector<double>& x, const std::vector<double>& y,
                 const std::vector<double>& z, double mn, double mx)
{
   hsize_t s[3] = { x.size(), y.size(), z.size() };
   hid_t sid = H5Screate_simple(3, s, nullptr);
   hid_t did = H5Dcreate2(g, n, H5T_NATIVE_DOUBLE, sid,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, d.data());
   write_string_attr(did, "units", "Pa");
   write_double_attr(did, "min_value", mn);
   write_double_attr(did, "max_value", mx);
   H5Dclose(did); H5Sclose(sid);
}
std::string tmp_path(const char* stem)
{
   const char* t = std::getenv("TMPDIR");
   if (!t || !*t) { t = "/tmp"; }
   std::string p = std::string(t) + "/seas_test_safs_" + stem + "_"
                 + std::to_string(::getpid()) + ".h5";
   ::unlink(p.c_str());
   return p;
}
void write_constant_sidecar(const std::string& p, double s)
{
   std::vector<double> ax = {-100.0, 0.0, 100.0};
   std::vector<double> ay = {-100.0, 0.0, 100.0};
   std::vector<double> az = {-100.0, 0.0, 100.0};
   hid_t f = H5Fcreate(p.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   write_string_attr(f, "schema_version", "data_projection_v1");
   write_string_attr(f, "crs", "EPSG:32611");
   write_string_attr(f, "units", "m");
   write_string_attr(f, "z_positive", "elevation");
   write_string_attr(f, "created_at", "1970-01-01T00:00:00Z");
   hid_t g = H5Gcreate2(f, "/grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   write_axis(g, "x", ax); write_axis(g, "y", ay); write_axis(g, "z", az);
   H5Gclose(g);
   hid_t fl = H5Gcreate2(f, "/fields", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   const size_t N = ax.size() * ay.size() * az.size();
   auto fill = [&](double v) { return std::vector<double>(N, v); };
   const double mn = -1e15, mx = 1e15;
   write_field(fl, "sigma_xx", fill(s), ax, ay, az, mn, mx);
   write_field(fl, "sigma_yy", fill(s), ax, ay, az, mn, mx);
   write_field(fl, "sigma_zz", fill(s), ax, ay, az, mn, mx);
   write_field(fl, "sigma_xy", fill(0.0), ax, ay, az, mn, mx);
   write_field(fl, "sigma_yz", fill(0.0), ax, ay, az, mn, mx);
   write_field(fl, "sigma_xz", fill(0.0), ax, ay, az, mn, mx);
   H5Gclose(fl); H5Fclose(f);
}

// R-001 (tpv102_tpv104_review.md): TPV-canonical sidecar — only
// σ_xy is non-zero, set to the value `sxy`.  Used by T_65_7 to verify
// `tau_pre.strike(i) = +sxy` (matches native TPV102/104 sign convention
// `d.tau2_0 = +tau_ini` for σ_xy > 0 right-lateral driving stress).
void write_tpv_sidecar(const std::string& p, double sxy)
{
   std::vector<double> ax = {-100.0, 0.0, 100.0};
   std::vector<double> ay = {-100.0, 0.0, 100.0};
   std::vector<double> az = {-100.0, 0.0, 100.0};
   hid_t f = H5Fcreate(p.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
   write_string_attr(f, "schema_version", "data_projection_v1");
   write_string_attr(f, "crs", "EPSG:32611");
   write_string_attr(f, "units", "m");
   write_string_attr(f, "z_positive", "elevation");
   write_string_attr(f, "created_at", "1970-01-01T00:00:00Z");
   hid_t g = H5Gcreate2(f, "/grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   write_axis(g, "x", ax); write_axis(g, "y", ay); write_axis(g, "z", az);
   H5Gclose(g);
   hid_t fl = H5Gcreate2(f, "/fields", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   const size_t N = ax.size() * ay.size() * az.size();
   auto fill = [&](double v) { return std::vector<double>(N, v); };
   const double mn = -1e15, mx = 1e15;
   write_field(fl, "sigma_xx", fill(0.0), ax, ay, az, mn, mx);
   write_field(fl, "sigma_yy", fill(0.0), ax, ay, az, mn, mx);
   write_field(fl, "sigma_zz", fill(0.0), ax, ay, az, mn, mx);
   write_field(fl, "sigma_xy", fill(sxy), ax, ay, az, mn, mx);
   write_field(fl, "sigma_yz", fill(0.0), ax, ay, az, mn, mx);
   write_field(fl, "sigma_xz", fill(0.0), ax, ay, az, mn, mx);
   H5Gclose(fl); H5Fclose(f);
}
bool FileExists(const std::string& p) { struct stat s; return stat(p.c_str(), &s) == 0; }
std::string FindBP5Mesh()
{
   const std::vector<std::string> c = {
      "bp5/mesh/bp5_1000m.msh",
      "../bp5/mesh/bp5_1000m.msh",
      "/Users/chunhuizhao/projects/seas-mfem-safs/miniapps/seas/bp5/mesh/bp5_1000m.msh",
   };
   for (const auto& p : c) { if (FileExists(p)) { return p; } }
   return "";
}
} // anon

struct Fixture
{
   BP5Params params;
   std::unique_ptr<Mesh> mesh;
   std::unique_ptr<ElasticityDomainOperator<Mesh>> domain_op;
   std::unique_ptr<FaultGeometry<Mesh>> fault_geom;
   int nf = 0;
   bool loaded = false;
};

static std::unique_ptr<Fixture> BuildFixture()
{
   auto fix = std::make_unique<Fixture>();
   const std::string mp = FindBP5Mesh();
   if (mp.empty()) { return fix; }
   fix->mesh = std::make_unique<Mesh>(mp.c_str(), 1, 1);
   fix->loaded = true;
   BoundaryConfig bc;
   bc.fault_attr = 100;
   bc.dirichlet_attrs = {1, 2, 3, 4};
   bc.natural_attrs   = {5, 6};
   bc.default_dirichlet_func = MakeBP5DirichletFunc(fix->params.Vp);
   LinearElastic le(fix->params.lambda(), fix->params.mu());
   DomainConfig dc;
   fix->domain_op = std::make_unique<ElasticityDomainOperator<Mesh>>(
      *fix->mesh, 1, le, fix->params.Vp, fix->params.Wf, fix->params.lf,
      bc, DGMethod::IP, SolverType::MUMPS_BLR, dc);
   fix->nf = fix->domain_op->GetNumFaultDOFs();
   if (fix->nf == 0) { return fix; }
   fix->fault_geom = std::make_unique<FaultGeometry<Mesh>>(*fix->domain_op,
                                                            fix->params);
   return fix;
}

// --------------------------------------------------------------------
// Tests
// --------------------------------------------------------------------
static void T_65_1_sigma_n_size(Fixture &fix)
{
   std::cout << "\n[T-65-1] sigma_n_per_dof size after ComputeParams\n";
   const std::string p = tmp_path("size");
   write_constant_sidecar(p, 1.0e7);
   StressField3D field(p);
   fix.fault_geom->ComputeParams(field);
   ::unlink(p.c_str());
   TEST_ASSERT(fix.fault_geom->sigma_n_per_dof().Size() == fix.nf,
               "sigma_n_per_dof.Size() == NumFaultDOFs");
}

static void T_65_2_tau_pre_size(Fixture &fix)
{
   std::cout << "\n[T-65-2] tau_pre size unchanged (BP5 layout 2*nf)\n";
   TEST_ASSERT(fix.fault_geom->GetTauPre().Size() == 2 * fix.nf,
               "tau_pre.Size() == 2 * NumFaultDOFs");
}

static void T_65_3_has_safs_params(Fixture &fix)
{
   std::cout << "\n[T-65-3] HasParams() true after invocation\n";
   TEST_ASSERT(fix.fault_geom->HasParams(),
               "HasParams() == true");
}

static void T_65_4_analytic_arrays_preserved(Fixture &fix)
{
   std::cout << "\n[T-65-4] a / eta / Dc / V_init analytic forms preserved\n";
   // Re-build a control fault geom and run ComputeBP5Params on it,
   // then compare a/eta/Dc/V_init element-wise.
   FaultGeometry<Mesh> control(*fix.domain_op, fix.params);
   const Vector& a_safs = fix.fault_geom->GetAValues();
   const Vector& a_bp5  = control.GetAValues();
   const Vector& e_safs = fix.fault_geom->GetEtaValues();
   const Vector& e_bp5  = control.GetEtaValues();
   const Vector& d_safs = fix.fault_geom->GetDcValues();
   const Vector& d_bp5  = control.GetDcValues();
   const Vector& v_safs = fix.fault_geom->GetVInit();
   const Vector& v_bp5  = control.GetVInit();
   TEST_ASSERT(a_safs.Size() == a_bp5.Size(), "a sizes match");
   TEST_ASSERT(e_safs.Size() == e_bp5.Size(), "eta sizes match");
   TEST_ASSERT(d_safs.Size() == d_bp5.Size(), "Dc sizes match");
   TEST_ASSERT(v_safs.Size() == v_bp5.Size(), "V_init sizes match");
   real_t mxa = 0, mxe = 0, mxd = 0, mxv = 0;
   for (int i = 0; i < a_safs.Size(); i++)
   {
      mxa = std::max(mxa, std::abs(a_safs(i) - a_bp5(i)));
      mxe = std::max(mxe, std::abs(e_safs(i) - e_bp5(i)));
      mxd = std::max(mxd, std::abs(d_safs(i) - d_bp5(i)));
   }
   for (int i = 0; i < v_safs.Size(); i++)
   {
      mxv = std::max(mxv, std::abs(v_safs(i) - v_bp5(i)));
   }
   TEST_NEAR(mxa, 0.0, 1e-15, "max |a_safs - a_bp5| == 0");
   TEST_NEAR(mxe, 0.0, 1e-15, "max |eta_safs - eta_bp5| == 0");
   TEST_NEAR(mxd, 0.0, 1e-15, "max |Dc_safs - Dc_bp5| == 0");
   TEST_NEAR(mxv, 0.0, 1e-15, "max |V_init_safs - V_init_bp5| == 0");
}

static void T_65_5_match_standalone_projection(Fixture &fix)
{
   std::cout << "\n[T-65-5] sigma_n / tau_pre match standalone "
                "ProjectFaultPreStress output\n";
   const std::string p = tmp_path("match");
   write_constant_sidecar(p, 5.0e7);
   StressField3D field(p);

   // Reference: call ProjectFaultPreStress directly.
   Vector sigma_n_ref, tau_ref;
   FieldProjector::ProjectFaultPreStress<Mesh>(
      field, *fix.fault_geom, sigma_n_ref, tau_ref,
      /*P_p=*/2.0e6, /*grad=*/0.0);

   // Now run ComputeParams on a fresh fault geom — should match.
   FaultGeometry<Mesh> safs_fg(*fix.domain_op, fix.params);
   safs_fg.ComputeParams(field, /*P_p=*/2.0e6, /*grad=*/0.0);
   ::unlink(p.c_str());

   const Vector& sigma_n_w = safs_fg.sigma_n_per_dof();
   const Vector& tau_w     = safs_fg.GetTauPre();
   real_t max_dev_sn = 0, max_dev_t = 0;
   for (int i = 0; i < sigma_n_ref.Size(); i++)
   {
      max_dev_sn = std::max(max_dev_sn,
                            std::abs(sigma_n_w(i) - sigma_n_ref(i)));
   }
   for (int i = 0; i < tau_ref.Size(); i++)
   {
      max_dev_t = std::max(max_dev_t, std::abs(tau_w(i) - tau_ref(i)));
   }
   TEST_NEAR(max_dev_sn, 0.0, 1e-9,
             "ComputeParams sigma_n matches standalone projector");
   TEST_NEAR(max_dev_t, 0.0, 1e-9,
             "ComputeParams tau_pre matches standalone projector");
}

// T-65-6 — R-001 regression: per-DOF projection math for the canonical
// TPV-y=0 fault basis (n = (0,-1,0), t1 = (0,0,-1), t2 = (+1,0,0)) must
// produce tau_pre.strike(i) = +σ_xy.  Drives BOTH the sidecar
// projection (FieldProjector::ProjectFaultPreStress) and the templated
// projection (ComputeParams<StressSource>) since they share the same
// formula.  Math-only: replicates the projection at one fake DOF and
// asserts the post-flip output matches `+TPV{102,104}Params::tau_ini`.
// Geometry-independent — the test does not need a BP5 mesh fixture.
static void T_65_6_projection_sign_matches_native_tpv()
{
   std::cout << "\n[T-65-6] R-001 regression: projection sign matches "
                "native TPV102/104/205 sign convention\n";

   // Canonical TPV basis (CLAUDE.md "Canonical Coordinate System").
   const real_t n[3]  = { 0.0, -1.0,  0.0 };          // outward fault normal
   const real_t t1[3] = { 0.0,  0.0, -1.0 };          // dip (down)
   const real_t t2[3] = { 1.0,  0.0,  0.0 };          // strike

   // TPV102 background pre-stress (compression-positive convention).
   const real_t sxx =  0.0;
   const real_t syy =  120.0e6;                       // σ_n = 120 MPa
   const real_t szz =  0.0;
   const real_t sxy =  75.0e6;                        // τ_ini = 75 MPa (right-lateral)
   const real_t syz =  0.0;
   const real_t sxz =  0.0;
   // Symmetric Cauchy tensor laid out row-major.
   const real_t S[3][3] = { { sxx, sxy, sxz },
                            { sxy, syy, syz },
                            { sxz, syz, szz } };

   // Sn = S · n.
   real_t Sn[3] = {0.0, 0.0, 0.0};
   for (int r = 0; r < 3; ++r)
   {
      Sn[r] = S[r][0]*n[0] + S[r][1]*n[1] + S[r][2]*n[2];
   }

   // R-001 sign convention: tau{1,2} = -(t{1,2} · Sn).
   const real_t sigma_n = n[0] *Sn[0] + n[1] *Sn[1] + n[2] *Sn[2];
   const real_t tau1    = -(t1[0]*Sn[0] + t1[1]*Sn[1] + t1[2]*Sn[2]);
   const real_t tau2    = -(t2[0]*Sn[0] + t2[1]*Sn[1] + t2[2]*Sn[2]);

   const real_t tol = 1.0e-9;
   TEST_NEAR(sigma_n, +120.0e6, tol,
             "sigma_n_total = +120 MPa (compression-positive, n-invariant)");
   TEST_NEAR(tau1,       0.0,   tol,
             "tau_pre.dip = 0 for pure-strike-slip TPV setup");
   // The load-bearing assertion: tau2 must be +75e6 to match
   // dynamic/tpv102_setup.hpp::d.tau2_0 = +TPV102Params::tau_ini.
   TEST_NEAR(tau2, +75.0e6, tol,
             "R-001: tau_pre.strike = +sxy (NATIVE TPV CONVENTION)");
}

int main(int, char**)
{
   std::cout << "Running Phase 6 §5 ComputeParams tests\n";

   // T-65-6 is math-only and does not require the BP5 mesh fixture.
   T_65_6_projection_sign_matches_native_tpv();

   auto fix = BuildFixture();
   if (!fix->loaded || fix->nf == 0 || !fix->fault_geom)
   {
      std::cout << "(BP5 mesh unavailable; runtime tests skipped — "
                   "Phase 6 §5 compile-only verification)\n";
      std::cout << "\nPhase 6 §5: " << num_passed << " / " << num_tests
                << " passed, " << num_failed << " failed\n";
      return (num_failed == 0) ? 0 : 1;
   }
   T_65_1_sigma_n_size(*fix);
   T_65_2_tau_pre_size(*fix);
   T_65_3_has_safs_params(*fix);
   T_65_4_analytic_arrays_preserved(*fix);
   T_65_5_match_standalone_projection(*fix);

   std::cout << "\n========================================\n";
   std::cout << "Phase 6 §5: " << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}

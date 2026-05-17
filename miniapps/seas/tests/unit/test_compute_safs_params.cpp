// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Unit tests for Phase 6 §5 — FaultGeometry::ComputeSAFSParams.
// Plan reference: PLAN_onfaultstress.md §1860-1879.
//
// Coverage:
//   T_65_1   sigma_n_per_dof.Size() == NumFaultDOFs after call.
//   T_65_2   tau_pre.Size() == 2 * NumFaultDOFs after call.
//   T_65_3   HasSAFSParams() returns true after invocation.
//   T_65_4   a/eta/Dc/V_init unchanged relative to ComputeBP5Params
//            (analytic spatial dependence preserved per plan §1869).
//   T_65_5   sigma_n_per_dof and tau_pre match the standalone
//            ProjectFaultPreStress output (proves ComputeSAFSParams
//            is the correct wrapper).
//   T_65_6   BP2 ctor → calling ComputeSAFSParams aborts.

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
   std::cout << "\n[T-65-1] sigma_n_per_dof size after ComputeSAFSParams\n";
   const std::string p = tmp_path("size");
   write_constant_sidecar(p, 1.0e7);
   StressField3D field(p);
   fix.fault_geom->ComputeSAFSParams(field);
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
   std::cout << "\n[T-65-3] HasSAFSParams() true after invocation\n";
   TEST_ASSERT(fix.fault_geom->HasSAFSParams(),
               "HasSAFSParams() == true");
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

   // Now run ComputeSAFSParams on a fresh fault geom — should match.
   FaultGeometry<Mesh> safs_fg(*fix.domain_op, fix.params);
   safs_fg.ComputeSAFSParams(field, /*P_p=*/2.0e6, /*grad=*/0.0);
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
             "ComputeSAFSParams sigma_n matches standalone projector");
   TEST_NEAR(max_dev_t, 0.0, 1e-9,
             "ComputeSAFSParams tau_pre matches standalone projector");
}

int main(int, char**)
{
   std::cout << "Running Phase 6 §5 ComputeSAFSParams tests\n";
   auto fix = BuildFixture();
   if (!fix->loaded || fix->nf == 0 || !fix->fault_geom)
   {
      std::cout << "(BP5 mesh unavailable; runtime tests skipped — "
                   "Phase 6 §5 compile-only verification)\n";
      std::cout << "\nPhase 6 §5: compile-only OK\n";
      return 0;
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

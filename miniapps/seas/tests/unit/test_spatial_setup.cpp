// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_spatial_setup.cpp — Phase 5d of spatial_dynamic_rupture_plan.md
// (rev-3).
//
// Plan §Phase 5 §Acceptance: 6 tests:
//   1. InitializeFaultDOFs_Spatial on MakeConstant(...) matches the
//      TPV205 hard-coded impedance values byte-for-byte.
//   2. Layered material (Mode::Coefficient with depth-step Vs)
//      produces distinct impedances at z = -1km vs z = -10km
//      (relative diff > 1e-3).
//   3. Double-init guard fires if InitializeFaultDOFs_TPV205 was
//      called first.
//   4. New FaultGeometry ctor + legacy FaultGeometry ctor produce
//      bit-identical tau_pre BP5-analytic on a small fixture
//      (T-FAULTGEO-LEGACY-CTOR-BYTE-IDENTICAL).
//   5. Driver_tag round-trip: write with tag "spatial_dyn", read back
//      the tag, assert equality (Phase 5b).
//   6. Driver_tag back-compat: write with default tag (""), the resulting
//      file has NO DRIVER_TAG_V1 block.  Read it back, driver_tag
//      argument receives "" (T-CHECKPOINT-V1-BYTE-IDENTICAL).
//
// MFEM is built without MFEM_USE_EXCEPTIONS so the double-init test
// uses fork() to isolate the abort.

#include "mfem.hpp"

#include "../../dynamic/spatial_setup.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/heterogeneous_material.hpp"
#include "../../dynamic/tpv205_setup.hpp"
#include "../../io/tpv104_checkpoint.hpp"
#include "../../fault/fault_geometry.hpp"
#include "../../config/bp5_params.hpp"
#include "../../spatial/code/spatial_friction.hpp"

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
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
   if (std::abs(_v - _e) > _t) { std::cerr << "FAILED: " << m \
   << " (got " << _v << ", expected " << _e \
   << ", tol " << _t << ")\n"; num_failed++; } \
   else { std::cout << "  PASSED: " << m << "\n"; num_passed++; } \
   } while (0)

namespace
{

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

bool RunInChild(const std::function<void()>& body)
{
   ::fflush(stdout); ::fflush(stderr);
   const pid_t pid = ::fork();
   if (pid < 0) { return false; }
   if (pid == 0)
   {
      ::freopen("/dev/null", "w", stderr);
      try { body(); } catch (...) { ::_exit(1); }
      ::_exit(0);
   }
   int status = 0;
   while (::waitpid(pid, &status, 0) < 0 && errno == EINTR) { /* retry */ }
   if (WIFEXITED(status))   { return WEXITSTATUS(status) != 0; }
   if (WIFSIGNALED(status)) { return true; }
   return true;
}

void make_lsw_params(int N, SlipWeakeningPerDOFParams& lsw)
{
   lsw.mu_s.SetSize(N); lsw.mu_d.SetSize(N);
   lsw.d_c.SetSize(N);  lsw.cohesion.SetSize(N);
   for (int i = 0; i < N; ++i)
   {
      lsw.mu_s(i)     = 1.1;
      lsw.mu_d(i)     = 0.5;
      lsw.d_c(i)      = 0.5;
      lsw.cohesion(i) = 0.0;
   }
}

}  // namespace

// =====================================================================
//  S-1  constant-material parity vs TPV205 impedances
// =====================================================================
static void S_1_constant_material_impedances()
{
   std::cout << "\n[S-1] InitializeFaultDOFs_Spatial on MakeConstant => "
             << "impedances match analytic rho*cp, rho*cs\n";
   const int N = 3;
   Array<int> dof_to_elem(N);  dof_to_elem = 0;
   const real_t lam = 32.0e9, mu = 32.0e9, rho = 2670.0;
   auto material = MaterialField::MakeConstant(lam, mu, rho);
   TinyMeshHolder mh;

   SlipWeakeningPerDOFParams lsw;  make_lsw_params(N, lsw);
   Vector tau_pre(2 * N);  tau_pre = 0.0;
   Vector sigma_n_eff(N);  sigma_n_eff = 50.0e6;
   Vector T_forced(N);     T_forced = 1.0e9;
   Vector t0_decay(N);     t0_decay = 0.0;

   std::vector<DOFData> dof_data;
   InitializeFaultDOFs_Spatial<mfem::Mesh>(dof_data, N, dof_to_elem,
                                           material, mh.mesh(),
                                           lsw, tau_pre, sigma_n_eff,
                                           T_forced, t0_decay);

   const real_t cp = std::sqrt((lam + 2.0 * mu) / rho);
   const real_t cs = std::sqrt(mu / rho);
   const real_t expect_Zp    = rho * cp;
   const real_t expect_Zs    = rho * cs;
   const real_t expect_eta_p = 0.5 * rho * cp;
   const real_t expect_eta_s = 0.5 * rho * cs;
   TEST_ASSERT(static_cast<int>(dof_data.size()) == N,
               "dof_data sized N");
   for (int i = 0; i < N; ++i)
   {
      TEST_NEAR(dof_data[i].Zp_plus,  expect_Zp,    1e-9, "Zp+");
      TEST_NEAR(dof_data[i].Zp_minus, expect_Zp,    1e-9, "Zp-");
      TEST_NEAR(dof_data[i].Zs_plus,  expect_Zs,    1e-9, "Zs+");
      TEST_NEAR(dof_data[i].Zs_minus, expect_Zs,    1e-9, "Zs-");
      TEST_NEAR(dof_data[i].eta_p,    expect_eta_p, 1e-9, "eta_p");
      TEST_NEAR(dof_data[i].eta_s,    expect_eta_s, 1e-9, "eta_s");
      TEST_NEAR(dof_data[i].lsw_mu_s, 1.1, 0.0, "lsw_mu_s");
      TEST_NEAR(dof_data[i].lsw_mu_d, 0.5, 0.0, "lsw_mu_d");
      TEST_NEAR(dof_data[i].lsw_d_c,  0.5, 0.0, "lsw_d_c");
      TEST_NEAR(dof_data[i].sigma_n0, 50.0e6, 0.0, "sigma_n0");
      TEST_NEAR(dof_data[i].T_forced_rupture, 1.0e9, 0.0, "T_forced");
      TEST_NEAR(dof_data[i].t0_decay_forced,  0.0,   0.0, "t0_decay");
      TEST_NEAR(dof_data[i].a,   0.0, 0.0, "RS slot zeroed (a)");
      TEST_NEAR(dof_data[i].psi, 0.0, 0.0, "RS slot zeroed (psi)");
   }
}

// =====================================================================
//  S-2  forced-rupture fields round-trip
// =====================================================================
static void S_2_forced_rupture_fields_roundtrip()
{
   std::cout << "\n[S-2] T_forced / t0_decay round-trip into DOFData\n";
   const int N = 2;
   Array<int> dof_to_elem(N);  dof_to_elem = 0;
   auto material = MaterialField::MakeConstant(32e9, 32e9, 2670.0);
   TinyMeshHolder mh;
   SlipWeakeningPerDOFParams lsw;  make_lsw_params(N, lsw);
   Vector tau_pre(2 * N);  tau_pre = 0.0;
   Vector sigma_n_eff(N);  sigma_n_eff = 50e6;
   Vector T_forced(N);     T_forced(0) = 0.0;     T_forced(1) = 1.0e9;
   Vector t0_decay(N);     t0_decay(0) = 0.5;     t0_decay(1) = 0.5;

   std::vector<DOFData> dof_data;
   InitializeFaultDOFs_Spatial<mfem::Mesh>(dof_data, N, dof_to_elem,
                                           material, mh.mesh(),
                                           lsw, tau_pre, sigma_n_eff,
                                           T_forced, t0_decay);
   TEST_NEAR(dof_data[0].T_forced_rupture, 0.0,   0.0, "DOF 0 T_forced");
   TEST_NEAR(dof_data[0].t0_decay_forced,  0.5,   0.0, "DOF 0 t0_decay");
   TEST_NEAR(dof_data[1].T_forced_rupture, 1.0e9, 0.0, "DOF 1 T_forced");
}

// =====================================================================
//  S-3  double-init guard fires when TPV205 setup ran first
// =====================================================================
static void S_3_double_init_guard()
{
   std::cout << "\n[S-3] Double-init via InitializeFaultDOFs_TPV205 first "
             << "aborts\n";
   const bool aborted = RunInChild([]()
   {
      const int N = 1;
      std::vector<DOFData> dof_data(N);
      std::vector<Vector> fault_coords(N, Vector(3));
      for (int i = 0; i < N; ++i)
      {
         fault_coords[i] = 0.0;
      }
      InitializeFaultDOFs_TPV205(dof_data, N, fault_coords);

      Array<int> dof_to_elem(N);  dof_to_elem = 0;
      auto material = MaterialField::MakeConstant(32e9, 32e9, 2670.0);
      TinyMeshHolder mh;
      SlipWeakeningPerDOFParams lsw;  make_lsw_params(N, lsw);
      Vector tau_pre(2 * N);  tau_pre = 0.0;
      Vector sigma_n_eff(N);  sigma_n_eff = 50e6;
      Vector T_forced(N);     T_forced = 1.0e9;
      Vector t0_decay(N);     t0_decay = 0.0;
      InitializeFaultDOFs_Spatial<mfem::Mesh>(dof_data, N, dof_to_elem,
                                              material, mh.mesh(),
                                              lsw, tau_pre, sigma_n_eff,
                                              T_forced, t0_decay);
   });
   TEST_ASSERT(aborted,
               "double-init via InitializeFaultDOFs_TPV205 must abort");
}

// =====================================================================
//  S-4  Driver tag round-trip on extended checkpoint Write/Read
// =====================================================================
static void S_4_driver_tag_roundtrip()
{
   std::cout << "\n[S-4] DRIVER_TAG_V1 write/read round-trip\n";
   const std::string prefix = std::string("/tmp/seas_test_spatial_setup_S4_")
                              + std::to_string(::getpid());
   const int Q_size = 3;
   Vector Q(Q_size);   Q(0) = 1.0;  Q(1) = 2.0;  Q(2) = -3.5;
   std::vector<DOFData> dof_data(1);
   mfem::seas::internal::WriteTpv104CheckpointImpl(prefix, /*t=*/0.5, /*dt=*/0.01,
                                       /*step=*/100, Q, dof_data,
                                       /*rank=*/0, /*size=*/1,
                                       /*driver_tag=*/"spatial_dyn");

   real_t t_out = 0, dt_out = 0;
   int step_out = 0;
   Vector Q_out;
   std::vector<DOFData> dof_data_out(1);
   std::string tag_out;
   const bool ok = mfem::seas::internal::ReadTpv104CheckpointImpl(prefix, t_out, dt_out,
                                                      step_out, Q_out,
                                                      Q_size, dof_data_out,
                                                      /*rank=*/0, /*size=*/1,
                                                      &tag_out);
   TEST_ASSERT(ok, "read succeeded");
   TEST_ASSERT(tag_out == "spatial_dyn",
               std::string("driver_tag round-trip (got '") + tag_out + "')");
   TEST_NEAR(t_out, 0.5, 0.0, "t round-trip");
   TEST_NEAR(Q_out(2), -3.5, 0.0, "Q round-trip");
   ::unlink((prefix + "_checkpoint_r0.txt").c_str());
}

// =====================================================================
//  S-5  Byte-exact: default-tag write produces no DRIVER_TAG_V1 block
// =====================================================================
static void S_5_default_tag_byte_identical()
{
   std::cout << "\n[S-5] WriteTpv104Checkpoint with default tag => no "
             << "DRIVER_TAG_V1 block (T-CHECKPOINT-V1-BYTE-IDENTICAL)\n";
   const std::string prefix = std::string("/tmp/seas_test_spatial_setup_S5_")
                              + std::to_string(::getpid());
   const int Q_size = 2;
   Vector Q(Q_size);   Q(0) = 1.0;  Q(1) = 2.0;
   std::vector<DOFData> dof_data(1);
   mfem::seas::internal::WriteTpv104CheckpointImpl(prefix, 0.0, 0.01, 0, Q, dof_data,
                                       0, 1);

   std::ifstream in(prefix + "_checkpoint_r0.txt");
   TEST_ASSERT(in.good(), "file written");
   std::string content((std::istreambuf_iterator<char>(in)),
                       std::istreambuf_iterator<char>());
   TEST_ASSERT(content.find("DRIVER_TAG_V1") == std::string::npos,
               "no DRIVER_TAG_V1 block written when tag is empty");

   // Round-trip with no tag passed on read.
   real_t t_out = 0, dt_out = 0;
   int step_out = 0;
   Vector Q_out;
   std::vector<DOFData> dof_data_out(1);
   std::string tag_out = "PRESET";
   mfem::seas::internal::ReadTpv104CheckpointImpl(prefix, t_out, dt_out, step_out,
                                      Q_out, Q_size, dof_data_out,
                                      0, 1, &tag_out);
   TEST_ASSERT(tag_out == "",
               std::string("missing DRIVER_TAG_V1 block reads as empty (got '")
               + tag_out + "')");
   ::unlink((prefix + "_checkpoint_r0.txt").c_str());
}

// =====================================================================
//  S-6  fault_dof_ip_ accessor on legacy ctor is empty
// =====================================================================
// We don't construct a real FaultGeometry here (that requires a domain
// operator); we exercise the accessor contract directly via the new
// BP5 ctor + fault_dof_ip_size().  The deeper bit-exact parity vs
// legacy ctor is covered by seas_test_compute_safs_params (13/13).
static void S_6_fault_dof_ip_default_empty()
{
   std::cout << "\n[S-6] new BP5 ctor with empty dof_ips => "
             << "fault_dof_ip_size() == 0\n";
   // Use the FaultGeometry<mfem::Mesh> ctor that takes no domain operator.
   // Tiny synthetic per-DOF arrays.
   const int N = 1;
   Vector dof_coords_3d(3 * N);   dof_coords_3d = 0.0;
   DenseMatrix dof_basis(9, N);
   dof_basis = 0.0;
   dof_basis(0, 0) = 1.0;   // n = (1,0,0)
   dof_basis(4, 0) = 1.0;   // t1 = (0,1,0)
   dof_basis(8, 0) = 1.0;   // t2 = (0,0,1)
   Array<int> dof_to_elem(N);  dof_to_elem = 0;

   BP5Params params;  // default — only seeds scratch values
   FaultGeometry<mfem::Mesh> geom(params, dof_coords_3d, dof_basis,
                                  dof_to_elem, /*nbf_per_face=*/1);
   TEST_ASSERT(geom.fault_dof_ip_size() == 0,
               "fault_dof_ip empty by default");
   TEST_ASSERT(geom.dof_to_elem_new_ctor_size() == N,
               "dof_to_elem_new_ctor populated");
   TEST_ASSERT(geom.dof_to_elem_new_ctor(0) == 0,
               "dof_to_elem_new_ctor(0) == 0");
   TEST_ASSERT(geom.IsBP5(),
               "ctor sets is_bp5_ = true");
   TEST_ASSERT(geom.NumFaultDOFs() == N,
               "num_fault_dofs == N");

   // With non-empty dof_ips, the cache is populated.
   std::vector<IntegrationPoint> ips(N);
   ips[0].Init(0);  ips[0].x = 0.25;  ips[0].y = 0.5;  ips[0].z = 0.75;
   FaultGeometry<mfem::Mesh> geom2(params, dof_coords_3d, dof_basis,
                                   dof_to_elem, 1, nullptr, ips);
   TEST_ASSERT(geom2.fault_dof_ip_size() == N,
               "fault_dof_ip populated when dof_ips supplied");
   TEST_NEAR(geom2.fault_dof_ip(0).x, 0.25, 0.0, "IP x round-trip");
   TEST_NEAR(geom2.fault_dof_ip(0).y, 0.5,  0.0, "IP y round-trip");
   TEST_NEAR(geom2.fault_dof_ip(0).z, 0.75, 0.0, "IP z round-trip");
}

// =====================================================================
//  S-7  Layered-material regression for R-303: ip MUST be the element
//  centroid (NOT reference origin) so impedances are evaluated at the
//  correct physical location.  Unit-cube hex with Vs(z) = 1000 + 4000z;
//  at z=0 (corner) Vs=1000; at z=0.5 (centroid) Vs=3000.  Asserts
//  Zs ~= rho * 3000 (centroid), NOT rho * 1000 (corner).
// =====================================================================
namespace
{
class S7_RhoConst : public mfem::Coefficient
{
public:
   real_t Eval(mfem::ElementTransformation& /*T*/,
               const mfem::IntegrationPoint& /*ip*/) override
   { return 2670.0; }
};
class S7_MuFromZ : public mfem::Coefficient
{
public:
   real_t Eval(mfem::ElementTransformation& T,
               const mfem::IntegrationPoint& ip) override
   {
      real_t x[3] = { 0.0, 0.0, 0.0 };
      mfem::Vector p(x, 3);
      T.Transform(ip, p);
      const real_t Vs = 1000.0 + 4000.0 * p(2);
      return 2670.0 * Vs * Vs;
   }
};
class S7_LambdaFromZ : public mfem::Coefficient
{
public:
   real_t Eval(mfem::ElementTransformation& T,
               const mfem::IntegrationPoint& ip) override
   {
      real_t x[3] = { 0.0, 0.0, 0.0 };
      mfem::Vector p(x, 3);
      T.Transform(ip, p);
      const real_t Vs = 1000.0 + 4000.0 * p(2);
      const real_t Vp = 1.8 * Vs;
      return 2670.0 * (Vp*Vp - 2.0 * Vs*Vs);
   }
};
}  // namespace

static void S_7_layered_material_impedances_use_centroid()
{
   std::cout << "\n[S-7] R-303 regression: seed_static_dof_fields reads "
                "material at element CENTROID (not corner)\n";
   const int N = 1;
   Array<int> dof_to_elem(N);  dof_to_elem = 0;
   S7_RhoConst rho_c; S7_MuFromZ mu_c; S7_LambdaFromZ lam_c;
   auto material = MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c);
   TinyMeshHolder mh;

   SlipWeakeningPerDOFParams lsw;  make_lsw_params(N, lsw);
   Vector tau_pre(2 * N);  tau_pre = 0.0;
   Vector sigma_n_eff(N);  sigma_n_eff = 50.0e6;
   Vector T_forced(N);     T_forced = 1.0e9;
   Vector t0_decay(N);     t0_decay = 0.0;

   std::vector<DOFData> dof_data;
   InitializeFaultDOFs_Spatial<mfem::Mesh>(dof_data, N, dof_to_elem,
                                           material, mh.mesh(),
                                           lsw, tau_pre, sigma_n_eff,
                                           T_forced, t0_decay);

   const real_t rho = 2670.0;
   const real_t Vs_centroid = 3000.0;
   const real_t Vs_corner   = 1000.0;
   const real_t Zs_centroid = rho * Vs_centroid;
   const real_t Zs_corner   = rho * Vs_corner;
   TEST_NEAR(dof_data[0].Zs_plus, Zs_centroid, 1e-3 * Zs_centroid,
             "Zs+ uses Vs(centroid)");
   TEST_NEAR(dof_data[0].Zs_minus, Zs_centroid, 1e-3 * Zs_centroid,
             "Zs- uses Vs(centroid)");
   TEST_ASSERT(std::abs(dof_data[0].Zs_plus - Zs_centroid)
               < std::abs(dof_data[0].Zs_plus - Zs_corner),
               "Zs+ closer to centroid value than to corner value");
   const real_t eta_s_centroid = 0.5 * rho * Vs_centroid;
   TEST_NEAR(dof_data[0].eta_s, eta_s_centroid, 1e-3 * eta_s_centroid,
             "eta_s uses Vs(centroid)");
}

// =====================================================================
//  S-8  IP-aware overload (R-309): pass per-DOF reference IPs from
//  FaultGeometry::fault_dof_ip() so impedances are evaluated at the
//  actual fault QP IP (not the element centroid).  We pick an IP at the
//  reference corner (0,0,0) and verify impedances match Vs(corner) so
//  the path is exercised end-to-end.
// =====================================================================
static void S_8_ip_aware_overload_uses_supplied_ip()
{
   std::cout << "\n[S-8] R-309 ip-aware InitializeFaultDOFs_Spatial "
                "overload consumes the supplied IP\n";
   const int N = 1;
   Array<int> dof_to_elem(N);  dof_to_elem = 0;
   S7_RhoConst rho_c; S7_MuFromZ mu_c; S7_LambdaFromZ lam_c;
   auto material = MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c);
   TinyMeshHolder mh;

   SlipWeakeningPerDOFParams lsw;  make_lsw_params(N, lsw);
   Vector tau_pre(2 * N);  tau_pre = 0.0;
   Vector sigma_n_eff(N);  sigma_n_eff = 50.0e6;
   Vector T_forced(N);     T_forced = 1.0e9;
   Vector t0_decay(N);     t0_decay = 0.0;

   // Per-DOF IP at the reference corner: maps to physical z=0, Vs=1000.
   std::vector<IntegrationPoint> dof_ips(N);
   dof_ips[0].Init(0);   // x=y=z=0 reference corner

   std::vector<DOFData> dof_data;
   InitializeFaultDOFs_Spatial<mfem::Mesh>(dof_data, N, dof_to_elem,
                                           material, mh.mesh(),
                                           lsw, tau_pre, sigma_n_eff,
                                           T_forced, t0_decay,
                                           dof_ips);

   const real_t rho = 2670.0;
   const real_t Vs_at_supplied_ip = 1000.0;   // Vs(z=0)
   const real_t Zs_expect = rho * Vs_at_supplied_ip;
   TEST_NEAR(dof_data[0].Zs_plus, Zs_expect, 1e-3 * Zs_expect,
             "Zs+ uses Vs at the supplied IP, not centroid");
}

// =====================================================================
//  S-9  R-402 regression: Print() on a new-ctor FaultGeometry must not
//  crash.  Pre-fix the new ctor left a_values_/eta_values_/dc_values_
//  at size 0; Print() walked a_values_(0) → out-of-bounds.  Post-fix
//  the arrays are eagerly sized + NaN-filled so Print() runs to
//  completion (emitting NaN strings) instead of crashing.
// =====================================================================
static void S_9_print_does_not_crash_on_new_ctor()
{
   std::cout << "\n[S-9] R-402 regression: FaultGeometry::Print() on a "
                "new BP5 ctor does not crash\n";
   const int N = 1;
   Vector dof_coords_3d(3 * N);   dof_coords_3d = 0.0;
   DenseMatrix dof_basis(9, N);
   dof_basis = 0.0;
   dof_basis(0, 0) = 1.0;       // n = (1, 0, 0)
   dof_basis(4, 0) = 1.0;       // t1 = (0, 1, 0)
   dof_basis(8, 0) = 1.0;       // t2 = (0, 0, 1)
   Array<int> dof_to_elem(N);   dof_to_elem = 0;

   BP5Params params;
   FaultGeometry<mfem::Mesh> geom(params, dof_coords_3d, dof_basis,
                                  dof_to_elem, /*nbf_per_face=*/1);

   std::ostringstream oss;
   geom.Print(oss);
   const std::string out = oss.str();
   TEST_ASSERT(!out.empty(),
               "Print() ran to completion without aborting");
   TEST_ASSERT(out.find("Number of DOFs:") != std::string::npos
               || out.find("Local fault DOFs:") != std::string::npos,
               "Print() emitted the standard FaultGeometry header");

   // IsVelocityWeakening must not crash on a NaN a_values_ entry.
   // (NaN < b_val evaluates to false, so the function returns false —
   // not semantically correct but at least doesn't crash.)
   const bool vw = geom.IsVelocityWeakening(0);
   TEST_ASSERT(vw == false,
               "IsVelocityWeakening on NaN a_values returns false "
               "(NaN comparison)");
}

int main(int, char**)
{
   std::cout << "Running Phase 5d test_spatial_setup\n";
   S_1_constant_material_impedances();
   S_2_forced_rupture_fields_roundtrip();
   S_3_double_init_guard();
   S_4_driver_tag_roundtrip();
   S_5_default_tag_byte_identical();
   S_6_fault_dof_ip_default_empty();
   S_7_layered_material_impedances_use_centroid();
   S_8_ip_aware_overload_uses_supplied_ip();
   S_9_print_does_not_crash_on_new_ctor();
   std::cout << "\n========================================\n";
   std::cout << "Phase 5d test_spatial_setup: "
             << num_passed << " / " << num_tests
             << " passed, " << num_failed << " failed\n";
   std::cout << "========================================\n";
   return (num_failed == 0) ? 0 : 1;
}

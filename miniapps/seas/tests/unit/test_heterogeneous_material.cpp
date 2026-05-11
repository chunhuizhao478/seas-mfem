// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// Unit tests for MaterialField.
// Test catalog T-5-1 .. T-5-3 — original Constant + GridFunction modes
//                               (data_projection_feature_plan_v2.md).
// Test catalog T-5-4 .. T-5-7 — new Coefficient mode and the three
//                               sidecar-backed Coefficient subclasses
//                               (heterogeneous_material_plan.md, Phase 1).
//
// All fixtures are synthetic and serial.  GridFunction mode requires a
// `ParFiniteElementSpace`, which in turn needs `mfem::Mpi` initialised
// — we run via mpirun -np 1.

#include "mfem.hpp"

#include "../../dynamic/heterogeneous_material.hpp"
#include "../../io/material_coefficients.hpp"
#include "../../io/data_field_3d.hpp"

#include <hdf5.h>

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <iostream>
#include <memory>
#include <string>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0;
static int num_passed = 0;
static int num_failed = 0;

#define TEST_ASSERT(cond, msg) \
   do { num_tests++; \
        if (!(cond)) { \
           std::cerr << "FAILED: " << msg << " (line " << __LINE__ \
                     << ")\n"; num_failed++; } \
        else { std::cout << "  PASSED: " << msg << "\n"; num_passed++; } \
   } while (0)

#define TEST_NEAR(value, expected, tol, msg) \
   do { num_tests++; \
        const real_t _v = (value); const real_t _e = (expected); \
        const real_t _t = (tol); \
        if (std::abs(_v - _e) > _t) { \
           std::cerr << "FAILED: " << msg << " (got " << _v \
                     << ", exp " << _e << ", tol " << _t \
                     << ", line " << __LINE__ << ")\n"; num_failed++; } \
        else { std::cout << "  PASSED: " << msg << "\n"; num_passed++; } \
   } while (0)


// ----------------------------------------------------------------------
// T-5-1 — Constant mode
// ----------------------------------------------------------------------

static void T_5_1_constant_mode()
{
   std::cout << "\n[T-5-1] MaterialField::At in Constant mode\n";
   auto m = MaterialField::MakeConstant(/*lambda=*/10.0,
                                        /*mu=*/20.0,
                                        /*rho=*/30.0);
   TEST_ASSERT(m.mode == MaterialField::Mode::Constant,
               "T-5-1 mode == Constant");

   real_t la, mu, rho;
   m.At(/*elem=*/0, /*dof=*/0, la, mu, rho);
   TEST_NEAR(la,  10.0, 1.0e-15, "T-5-1 lambda");
   TEST_NEAR(mu,  20.0, 1.0e-15, "T-5-1 mu");
   TEST_NEAR(rho, 30.0, 1.0e-15, "T-5-1 rho");

   // Constant mode is element/dof independent — verify a different
   // (elem, dof) returns the same triple.
   m.At(/*elem=*/12345, /*dof=*/3, la, mu, rho);
   TEST_NEAR(la,  10.0, 1.0e-15, "T-5-1 lambda invariant under (elem, dof)");
   TEST_NEAR(mu,  20.0, 1.0e-15, "T-5-1 mu invariant under (elem, dof)");
   TEST_NEAR(rho, 30.0, 1.0e-15, "T-5-1 rho invariant under (elem, dof)");

   // MaxCpInElement: c_p = sqrt((lambda + 2 mu) / rho)
   //                     = sqrt((10 + 40) / 30) = sqrt(50 / 30)
   const real_t cp = m.MaxCpInElement(0);
   TEST_NEAR(cp, std::sqrt(50.0 / 30.0), 1.0e-12,
             "T-5-1 MaxCpInElement matches closed form");
}


// ----------------------------------------------------------------------
// T-5-2 — GridFunction mode with constant GFs
// ----------------------------------------------------------------------

static void T_5_2_grid_function_constant()
{
   std::cout << "\n[T-5-2] MaterialField::At in GridFunction mode (const GFs)\n";

   auto base = std::make_shared<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON,
                                  1.0, 1.0, 1.0));
   auto pmesh = std::make_shared<mfem::ParMesh>(MPI_COMM_WORLD, *base);
   mfem::H1_FECollection fec(/*p=*/1, /*dim=*/3);
   auto fes = std::make_shared<mfem::ParFiniteElementSpace>(pmesh.get(),
                                                            &fec);
   auto rho_gf    = std::make_shared<mfem::ParGridFunction>(fes.get());
   auto lambda_gf = std::make_shared<mfem::ParGridFunction>(fes.get());
   auto mu_gf     = std::make_shared<mfem::ParGridFunction>(fes.get());
   *rho_gf    = 30.0;
   *lambda_gf = 10.0;
   *mu_gf     = 20.0;

   auto m = MaterialField::MakeGridFunction(rho_gf, lambda_gf, mu_gf);
   TEST_ASSERT(m.mode == MaterialField::Mode::GridFunction,
               "T-5-2 mode == GridFunction");

   real_t la, mu, rho;
   m.At(/*elem=*/0, /*dof=*/0, la, mu, rho);
   TEST_NEAR(la,  10.0, 1.0e-12, "T-5-2 lambda matches constant GF");
   TEST_NEAR(mu,  20.0, 1.0e-12, "T-5-2 mu matches constant GF");
   TEST_NEAR(rho, 30.0, 1.0e-12, "T-5-2 rho matches constant GF");

   // Different (elem, dof) — same value because GFs are constant.
   const int last_elem = pmesh->GetNE() - 1;
   m.At(last_elem, /*dof=*/0, la, mu, rho);
   TEST_NEAR(la,  10.0, 1.0e-12, "T-5-2 lambda matches at last element");
   TEST_NEAR(mu,  20.0, 1.0e-12, "T-5-2 mu matches at last element");
   TEST_NEAR(rho, 30.0, 1.0e-12, "T-5-2 rho matches at last element");

   const real_t cp = m.MaxCpInElement(0);
   TEST_NEAR(cp, std::sqrt(50.0 / 30.0), 1.0e-12,
             "T-5-2 MaxCpInElement(constant GF) matches closed form");

   // Cross-check: GridFunction mode with constant fields equals
   // Constant mode with the same scalars, to within float precision.
   auto m_const = MaterialField::MakeConstant(10.0, 20.0, 30.0);
   real_t la_c, mu_c, rho_c;
   m_const.At(0, 0, la_c, mu_c, rho_c);
   TEST_NEAR(la,  la_c,  1.0e-12, "T-5-2 GF-mode == Constant-mode (lambda)");
   TEST_NEAR(mu,  mu_c,  1.0e-12, "T-5-2 GF-mode == Constant-mode (mu)");
   TEST_NEAR(rho, rho_c, 1.0e-12, "T-5-2 GF-mode == Constant-mode (rho)");
}


// ----------------------------------------------------------------------
// T-5-3 — MaxCpInElement with non-constant GFs
// ----------------------------------------------------------------------

static void T_5_3_max_cp_in_element()
{
   std::cout << "\n[T-5-3] MaterialField::MaxCpInElement (heterogeneous)\n";
   auto base = std::make_shared<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON,
                                  1.0, 1.0, 1.0));
   auto pmesh = std::make_shared<mfem::ParMesh>(MPI_COMM_WORLD, *base);
   mfem::H1_FECollection fec(1, 3);
   auto fes = std::make_shared<mfem::ParFiniteElementSpace>(pmesh.get(),
                                                            &fec);
   auto rho_gf    = std::make_shared<mfem::ParGridFunction>(fes.get());
   auto lambda_gf = std::make_shared<mfem::ParGridFunction>(fes.get());
   auto mu_gf     = std::make_shared<mfem::ParGridFunction>(fes.get());

   // Linear gradient in x: lambda(x) = 10 + 90*x, mu = 20, rho = 30.
   //   c_p(x) = sqrt(((10 + 90 x) + 40) / 30) = sqrt((50 + 90 x) / 30).
   class LambdaCoeff : public mfem::Coefficient {
   public:
      real_t Eval(mfem::ElementTransformation& T,
                  const mfem::IntegrationPoint& ip) override
      {
         mfem::Vector x(3);
         T.Transform(ip, x);
         return 10.0 + 90.0 * x[0];
      }
   } lambda_coef;
   *rho_gf = 30.0;
   *mu_gf  = 20.0;
   lambda_gf->ProjectCoefficient(lambda_coef);

   auto m = MaterialField::MakeGridFunction(rho_gf, lambda_gf, mu_gf);

   // Find an element somewhere in the middle (not at x = 0).  Since
   // the test mesh is structured 2x2x2 over [0, 1]^3, every element
   // has DOFs spanning 0.5 m in x — so c_p_max varies across DOFs
   // within each element.  For element 0 (x ∈ [0, 0.5]):
   //   c_p_max = sqrt((50 + 90 * 0.5) / 30) = sqrt(95/30).
   const real_t cp0 = m.MaxCpInElement(0);
   const real_t expected_cp0 = std::sqrt(95.0 / 30.0);
   TEST_NEAR(cp0, expected_cp0, 1.0e-10,
             "T-5-3 MaxCpInElement(elem=0) matches max over DOFs");

   // For an element in the right half (x ∈ [0.5, 1.0]), c_p_max
   // = sqrt((50 + 90 * 1.0) / 30) = sqrt(140/30).  Find such an element
   // by inspecting per-element x ranges.
   real_t cp_max_global = 0.0;
   for (int e = 0; e < pmesh->GetNE(); ++e)
   {
      const real_t cp = m.MaxCpInElement(e);
      if (cp > cp_max_global) { cp_max_global = cp; }
   }
   const real_t expected_global = std::sqrt(140.0 / 30.0);
   TEST_NEAR(cp_max_global, expected_global, 1.0e-10,
             "T-5-3 MaxCpInElement global max matches sqrt(140/30)");
}


// ----------------------------------------------------------------------
// Synthetic HDF5 sidecar helpers for T-5-6 and T-5-7.
//
// Schema-v1 layout — same as `tests/unit/test_data_field_3d.cpp` but
// inlined here so the test file is self-contained.  A single sidecar
// can hold multiple fields (Vp, Vs, density) — we write all three in
// one call to `write_three_field_sidecar` so a single HDF5 file
// satisfies the three `DataField3D` constructors.
// ----------------------------------------------------------------------

static void hdf5_write_string_attr(hid_t loc, const char* name,
                                    const std::string& value)
{
   hid_t atype = H5Tcopy(H5T_C_S1);
   H5Tset_size(atype, value.size() + 1);
   H5Tset_strpad(atype, H5T_STR_NULLTERM);
   hid_t aspace = H5Screate(H5S_SCALAR);
   hid_t aid = H5Acreate2(loc, name, atype, aspace,
                          H5P_DEFAULT, H5P_DEFAULT);
   H5Awrite(aid, atype, value.c_str());
   H5Aclose(aid);
   H5Sclose(aspace);
   H5Tclose(atype);
}

static void hdf5_write_double_attr(hid_t loc, const char* name,
                                    double v)
{
   hid_t aspace = H5Screate(H5S_SCALAR);
   hid_t aid = H5Acreate2(loc, name, H5T_NATIVE_DOUBLE, aspace,
                          H5P_DEFAULT, H5P_DEFAULT);
   H5Awrite(aid, H5T_NATIVE_DOUBLE, &v);
   H5Aclose(aid);
   H5Sclose(aspace);
}

static void hdf5_write_axis(hid_t loc, const char* name,
                             const std::vector<double>& v)
{
   hsize_t n = v.size();
   hid_t sid = H5Screate_simple(1, &n, nullptr);
   hid_t did = H5Dcreate2(loc, name, H5T_NATIVE_DOUBLE, sid,
                          H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
   H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
            v.data());
   H5Dclose(did);
   H5Sclose(sid);
}

struct FieldSpec
{
   std::string name;
   double      value;     // constant value for the whole grid
   double      min_value;
   double      max_value;
   std::string units;
};

static void write_three_field_sidecar(const std::string& path,
                                       const std::vector<double>& x,
                                       const std::vector<double>& y,
                                       const std::vector<double>& z,
                                       const std::vector<FieldSpec>& fields)
{
   hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                          H5P_DEFAULT);
   hdf5_write_string_attr(file, "schema_version", "data_projection_v1");
   hdf5_write_string_attr(file, "crs",            "EPSG:32611");
   hdf5_write_string_attr(file, "units",          "m");
   hdf5_write_string_attr(file, "z_positive",     "elevation");
   hdf5_write_string_attr(file, "created_at",     "1970-01-01T00:00:00Z");

   hid_t grid = H5Gcreate2(file, "/grid", H5P_DEFAULT, H5P_DEFAULT,
                            H5P_DEFAULT);
   hdf5_write_axis(grid, "x", x);
   hdf5_write_axis(grid, "y", y);
   hdf5_write_axis(grid, "z", z);
   H5Gclose(grid);

   hid_t flds = H5Gcreate2(file, "/fields", H5P_DEFAULT, H5P_DEFAULT,
                            H5P_DEFAULT);
   hsize_t shape[3] = { x.size(), y.size(), z.size() };
   const size_t N = x.size() * y.size() * z.size();
   for (const FieldSpec& fs : fields)
   {
      std::vector<double> data(N, fs.value);
      hid_t sid = H5Screate_simple(3, shape, nullptr);
      hid_t did = H5Dcreate2(flds, fs.name.c_str(),
                             H5T_NATIVE_DOUBLE, sid,
                             H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
               data.data());
      hdf5_write_string_attr(did, "units",     fs.units);
      hdf5_write_double_attr(did, "min_value", fs.min_value);
      hdf5_write_double_attr(did, "max_value", fs.max_value);
      H5Dclose(did);
      H5Sclose(sid);
   }
   H5Gclose(flds);
   H5Fclose(file);
}

static std::string make_tmp_sidecar_path(const char* stem)
{
   const char* tmpdir = std::getenv("TMPDIR");
   if (!tmpdir || !*tmpdir) { tmpdir = "/tmp"; }
   std::string path = std::string(tmpdir)
                      + "/seas_test_heterogeneous_material_"
                      + std::string(stem) + "_"
                      + std::to_string(::getpid()) + ".h5";
   ::unlink(path.c_str());
   return path;
}

// Build a single-hex unit-cube serial Mesh covering [0, 1]^3 with one
// element.  Provides a well-defined `ElementTransformation` for
// element 0 with its centroid at physical (0.5, 0.5, 0.5).
static std::unique_ptr<mfem::Mesh> make_single_hex_mesh()
{
   return std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::HEXAHEDRON,
                                   1.0, 1.0, 1.0));
}

// Helper: fork-and-check that a callable triggers a MFEM_ABORT.
// Adapted from tests/unit/test_data_field_3d.cpp.  Returns true iff
// the child terminated abnormally (signal or non-zero exit).
static bool expect_abort(const std::function<void(void)>& f)
{
   pid_t pid = ::fork();
   if (pid == 0)
   {
      // Child: silence stderr/stdout, abort is expected.
      ::close(1);
      ::close(2);
      try { f(); }
      catch (...) { std::_Exit(1); }
      std::_Exit(0);
   }
   int status = 0;
   ::waitpid(pid, &status, 0);
   if (WIFSIGNALED(status))                              { return true; }
   if (WIFEXITED(status) && WEXITSTATUS(status) != 0)    { return true; }
   return false;
}


// ----------------------------------------------------------------------
// T-5-4 — Mode::Coefficient round-trip (using ConstantCoefficient)
// ----------------------------------------------------------------------

static void T_5_4_coefficient_mode_round_trip()
{
   std::cout << "\n[T-5-4] MaterialField::EvalAt in Coefficient mode\n";

   mfem::ConstantCoefficient lambda_c(1.0);
   mfem::ConstantCoefficient mu_c    (2.0);
   mfem::ConstantCoefficient rho_c   (3.0);

   auto m = MaterialField::MakeCoefficient(&lambda_c, &mu_c, &rho_c);
   TEST_ASSERT(m.mode == MaterialField::Mode::Coefficient,
               "T-5-4 mode == Coefficient");
   TEST_ASSERT(m.lambda_coef == &lambda_c, "T-5-4 lambda_coef stored");
   TEST_ASSERT(m.mu_coef     == &mu_c,     "T-5-4 mu_coef stored");
   TEST_ASSERT(m.rho_coef    == &rho_c,    "T-5-4 rho_coef stored");

   // Need an ElementTransformation + IntegrationPoint to evaluate.
   auto mesh = make_single_hex_mesh();
   mfem::ElementTransformation* T = mesh->GetElementTransformation(0);
   mfem::IntegrationPoint ip;
   ip.Set3(0.5, 0.5, 0.5);
   T->SetIntPoint(&ip);

   real_t la, mu, rho;
   m.EvalAt(0, *T, ip, la, mu, rho);
   TEST_NEAR(la,  1.0, 1.0e-12, "T-5-4 EvalAt → lambda");
   TEST_NEAR(mu,  2.0, 1.0e-12, "T-5-4 EvalAt → mu");
   TEST_NEAR(rho, 3.0, 1.0e-12, "T-5-4 EvalAt → rho");
}


// ----------------------------------------------------------------------
// T-5-5 — abort paths in Mode::Coefficient
//
// Both sub-tests fork-and-check.  They MUST run BEFORE mfem::Mpi::Init,
// otherwise the child's MFEM_ABORT call lands inside `mfem_error`'s
// MPI-initialised branch and invokes `MPI_Abort(MPI_COMM_WORLD, 1)`,
// which propagates the abort to the parent rank (mpirun -np 1) and
// hangs the test.  Before Mpi::Init, MFEM_ABORT dispatches to
// std::abort() → child exits non-zero → parent observes it cleanly.
// ----------------------------------------------------------------------

static void T_5_5_abort_paths_in_coefficient_mode()
{
   std::cout << "\n[T-5-5] MaterialField abort paths in Coefficient mode\n";

   // (a) At(elem, dof, ...) in Mode::Coefficient.
   const bool aborted_on_at = expect_abort([](){
      mfem::ConstantCoefficient lambda_c(1.0);
      mfem::ConstantCoefficient mu_c    (2.0);
      mfem::ConstantCoefficient rho_c   (3.0);
      auto m = MaterialField::MakeCoefficient(&lambda_c, &mu_c, &rho_c);
      real_t la, mu, rho;
      m.At(/*elem=*/0, /*dof=*/0, la, mu, rho);
   });
   TEST_ASSERT(aborted_on_at,
               "T-5-5(a) MaterialField::At in Mode::Coefficient aborts");

   // (b) MaxCpInElement(elem, nullptr) in Mode::Coefficient.
   const bool aborted_on_null_T = expect_abort([](){
      mfem::ConstantCoefficient lambda_c(1.0);
      mfem::ConstantCoefficient mu_c    (2.0);
      mfem::ConstantCoefficient rho_c   (3.0);
      auto m = MaterialField::MakeCoefficient(&lambda_c, &mu_c, &rho_c);
      m.MaxCpInElement(/*elem=*/0, /*T=*/nullptr);
   });
   TEST_ASSERT(aborted_on_null_T,
               "T-5-5(b) MaxCpInElement(elem, nullptr) aborts in "
               "Coefficient mode");

   // (c) R-001 round-4: EvalAt called on a GridFunction-mode
   // MaterialField MUST abort (plan Phase 1 Detailed Req. 2).  We do
   // NOT need a valid ParGridFunction here — the EvalAt mode-branch
   // hits MFEM_ABORT BEFORE any GF dereference.  Construct the struct
   // by hand (bypassing the MakeGridFunction factory's non-null guard)
   // and pass a default-constructed IsoparametricTransformation: the
   // T/ip arguments are never inspected.
   const bool aborted_on_eval_gf = expect_abort([](){
      MaterialField m;
      m.mode = MaterialField::Mode::GridFunction;
      mfem::IsoparametricTransformation T_local;
      mfem::IntegrationPoint ip;
      ip.Set3(0.5, 0.5, 0.5);
      real_t la, mu, rho;
      m.EvalAt(/*elem=*/0, T_local, ip, la, mu, rho);
   });
   TEST_ASSERT(aborted_on_eval_gf,
               "T-5-5(c) EvalAt in Mode::GridFunction aborts");

   // (d) R-004 round-4: MakeCoefficient with any null pointer must
   // abort (plan Phase 1 Detailed Req. 6).  Pin the MFEM_VERIFY so a
   // future refactor cannot silently downgrade it to MFEM_ASSERT.
   const bool aborted_on_null_lambda = expect_abort([](){
      mfem::ConstantCoefficient mu_c (2.0);
      mfem::ConstantCoefficient rho_c(3.0);
      MaterialField::MakeCoefficient(/*lambda=*/nullptr, &mu_c, &rho_c);
   });
   TEST_ASSERT(aborted_on_null_lambda,
               "T-5-5(d) MakeCoefficient(nullptr, mu, rho) aborts");

   // (e) R-005 round-4: MakeGridFunction with any null shared_ptr must
   // abort (plan Phase 1 Detailed Req. 7, R-006 round-3).  Empty
   // shared_ptrs are bool-false; the MFEM_VERIFY fires WITHOUT ever
   // constructing a ParGridFunction, so no MPI is required.
   const bool aborted_on_null_gf = expect_abort([](){
      std::shared_ptr<mfem::ParGridFunction> rho_gf;     // empty
      std::shared_ptr<mfem::ParGridFunction> lambda_gf;  // empty
      std::shared_ptr<mfem::ParGridFunction> mu_gf;      // empty
      MaterialField::MakeGridFunction(rho_gf, lambda_gf, mu_gf);
   });
   TEST_ASSERT(aborted_on_null_gf,
               "T-5-5(e) MakeGridFunction(null shared_ptrs) aborts");
}


// ----------------------------------------------------------------------
// T-5-6 — MaxCpInElement in Mode::Coefficient
// ----------------------------------------------------------------------

static void T_5_6_max_cp_in_element_coefficient()
{
   std::cout << "\n[T-5-6] MaterialField::MaxCpInElement in Coefficient mode\n";

   // Build a sidecar with constant fields:
   //   Vp = 5000 m/s, Vs = 3000 m/s, rho = 2700 kg/m^3.
   //   μ = ρ Vs² = 2700 · 9e6 = 2.43e10
   //   λ = ρ Vp² − 2μ = 2700 · 25e6 − 2·2.43e10 = 6.75e10 − 4.86e10
   //                  = 1.89e10
   //   c_p = sqrt((λ + 2μ) / ρ) = sqrt(6.75e10 / 2700) = sqrt(2.5e7)
   //                            = 5000 m/s (matches Vp by definition).
   const std::string sidecar = make_tmp_sidecar_path("T-5-6");
   const std::vector<double> x = { 0.0, 0.5, 1.0 };
   const std::vector<double> y = { 0.0, 0.5, 1.0 };
   const std::vector<double> z = { 0.0, 0.5, 1.0 };
   const std::vector<FieldSpec> fields = {
      { "Vp",      5000.0,  100.0,  9000.0, "m/s"    },
      { "Vs",      3000.0,  100.0,  5000.0, "m/s"    },
      { "density", 2700.0, 1000.0,  3500.0, "kg/m^3" },
   };
   write_three_field_sidecar(sidecar, x, y, z, fields);

   DataField3D vp_field  (sidecar, "Vp");
   DataField3D vs_field  (sidecar, "Vs");
   DataField3D rho_field (sidecar, "density");

   LambdaFromSidecar lambda_coef(vp_field, vs_field, rho_field);
   MuFromSidecar     mu_coef    (vs_field, rho_field);
   RhoFromSidecar    rho_coef   (rho_field);

   auto m = MaterialField::MakeCoefficient(&lambda_coef, &mu_coef,
                                            &rho_coef);

   auto mesh = make_single_hex_mesh();
   mfem::ElementTransformation* T = mesh->GetElementTransformation(0);

   const real_t cp = m.MaxCpInElement(/*elem=*/0, T);
   const real_t expected_cp = 5000.0;
   TEST_NEAR(cp, expected_cp, 1.0e-6,
             "T-5-6 MaxCpInElement matches sqrt((lambda+2mu)/rho) = Vp");

   // The abort-on-nullptr-T case is exercised in T-5-5(b).

   ::unlink(sidecar.c_str());
}


// ----------------------------------------------------------------------
// T-5-6b — R-002 round-4: MaxCpInElement in Mode::Coefficient must NOT
// leave T's IntegrationPoint pointer pointing into the internal
// IntRules rule.  Pins the save/restore contract introduced in
// heterogeneous_material.cpp R-002 fix.
// ----------------------------------------------------------------------

static void T_5_6b_max_cp_preserves_T_ip()
{
   std::cout << "\n[T-5-6b] MaxCpInElement preserves caller's T IntPoint\n";

   const std::string sidecar = make_tmp_sidecar_path("T-5-6b");
   const std::vector<double> x = { 0.0, 0.5, 1.0 };
   const std::vector<double> y = { 0.0, 0.5, 1.0 };
   const std::vector<double> z = { 0.0, 0.5, 1.0 };
   const std::vector<FieldSpec> fields = {
      { "Vp",      5000.0,  100.0,  9000.0, "m/s"    },
      { "Vs",      3000.0,  100.0,  5000.0, "m/s"    },
      { "density", 2700.0, 1000.0,  3500.0, "kg/m^3" },
   };
   write_three_field_sidecar(sidecar, x, y, z, fields);

   DataField3D vp_field  (sidecar, "Vp");
   DataField3D vs_field  (sidecar, "Vs");
   DataField3D rho_field (sidecar, "density");

   LambdaFromSidecar lambda_coef(vp_field, vs_field, rho_field);
   MuFromSidecar     mu_coef    (vs_field, rho_field);
   RhoFromSidecar    rho_coef   (rho_field);
   auto m = MaterialField::MakeCoefficient(&lambda_coef, &mu_coef,
                                            &rho_coef);

   auto mesh = make_single_hex_mesh();
   mfem::ElementTransformation* T = mesh->GetElementTransformation(0);

   // Caller's deliberate IntegrationPoint, distinct from any qpoint of
   // the order-4 hex rule MaxCpInElement walks internally.
   mfem::IntegrationPoint ip_user;
   ip_user.Set3(0.25, 0.25, 0.25);
   T->SetIntPoint(&ip_user);

   const mfem::IntegrationPoint* const before = &T->GetIntPoint();
   (void) m.MaxCpInElement(/*elem=*/0, T);
   const mfem::IntegrationPoint* const after  = &T->GetIntPoint();

   TEST_ASSERT(before == after,
               "T-5-6b MaxCpInElement preserves T's IntegrationPoint");

   ::unlink(sidecar.c_str());
}


// ----------------------------------------------------------------------
// T-5-7 — LambdaFromSidecar / MuFromSidecar / RhoFromSidecar numerical
// ----------------------------------------------------------------------

static void T_5_7_sidecar_coefficients_numerical()
{
   std::cout << "\n[T-5-7] LambdaFromSidecar / MuFromSidecar / "
                "RhoFromSidecar Eval\n";

   const std::string sidecar = make_tmp_sidecar_path("T-5-7");
   const std::vector<double> x = { 0.0, 0.5, 1.0 };
   const std::vector<double> y = { 0.0, 0.5, 1.0 };
   const std::vector<double> z = { 0.0, 0.5, 1.0 };
   const std::vector<FieldSpec> fields = {
      { "Vp",      5000.0,  100.0,  9000.0, "m/s"    },
      { "Vs",      3000.0,  100.0,  5000.0, "m/s"    },
      { "density", 2700.0, 1000.0,  3500.0, "kg/m^3" },
   };
   write_three_field_sidecar(sidecar, x, y, z, fields);

   DataField3D vp_field  (sidecar, "Vp");
   DataField3D vs_field  (sidecar, "Vs");
   DataField3D rho_field (sidecar, "density");

   LambdaFromSidecar lambda_coef(vp_field, vs_field, rho_field);
   MuFromSidecar     mu_coef    (vs_field, rho_field);
   RhoFromSidecar    rho_coef   (rho_field);

   auto mesh = make_single_hex_mesh();
   mfem::ElementTransformation* T = mesh->GetElementTransformation(0);
   mfem::IntegrationPoint ip;
   ip.Set3(0.5, 0.5, 0.5);    // → physical (0.5, 0.5, 0.5)
   T->SetIntPoint(&ip);

   const real_t rho = rho_coef   .Eval(*T, ip);
   const real_t mu  = mu_coef    .Eval(*T, ip);
   const real_t la  = lambda_coef.Eval(*T, ip);

   // Per the plan: ρ = 2700, μ = ρ Vs² = 2700·9e6 = 2.43e10,
   //               λ = ρ Vp² − 2μ = 2700·(25e6 − 18e6) = 1.89e10.
   TEST_NEAR(rho, 2700.0,    1.0e-6,
             "T-5-7 RhoFromSidecar::Eval returns rho");
   TEST_NEAR(mu,  2.43e10,   1.0e-2,
             "T-5-7 MuFromSidecar::Eval = rho * Vs^2");
   TEST_NEAR(la,  1.89e10,   1.0e-2,
             "T-5-7 LambdaFromSidecar::Eval = rho * Vp^2 - 2 * mu");

   ::unlink(sidecar.c_str());
}


// ----------------------------------------------------------------------
// Main
// ----------------------------------------------------------------------

int main(int argc, char* argv[])
{
   std::cout << "===============================================\n";
   std::cout << "MaterialField unit tests (T-5-1..T-5-7)\n";
   std::cout << "===============================================\n";

   // Run fork+abort tests FIRST, before mfem::Mpi::Init.  See the
   // T-5-5 docstring for why: MPI_Abort propagates to the parent
   // rank if we abort while MPI is initialised.
   T_5_4_coefficient_mode_round_trip();
   T_5_5_abort_paths_in_coefficient_mode();
   T_5_6_max_cp_in_element_coefficient();   // serial Mesh + DataField3D; no MPI needed
   T_5_6b_max_cp_preserves_T_ip();          // R-002 round-4 regression
   T_5_7_sidecar_coefficients_numerical();  // serial Mesh + DataField3D; no MPI needed

   // GridFunction-mode tests need ParMesh, so init MPI now.
   mfem::Mpi::Init(argc, argv);
   T_5_1_constant_mode();
   T_5_2_grid_function_constant();
   T_5_3_max_cp_in_element();

   std::cout << "\n===============================================\n";
   std::cout << "Total tests: " << num_tests << "\n";
   std::cout << "Passed:      " << num_passed << "\n";
   std::cout << "Failed:      " << num_failed << "\n";
   if (num_failed) { std::cout << "\nSOME TESTS FAILED!\n"; return 1; }
   std::cout << "\nALL TESTS PASSED!\n";
   return 0;
}

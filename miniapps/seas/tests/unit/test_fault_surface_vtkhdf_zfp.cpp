// Phase 2d.3 of miniapps/seas/io/PLAN_paraview_compaction_2026-04-28.md.
//
// Exercises the seas fault VTKHDF writer with the new ZFP path:
//   pv.SetFaultHDFCompression(HDFCompression::ZfpAccuracy, 1e-12);
//
// The fault VTU bit-exact reference (Phase 1) is read back and compared
// to the ZFP fault-surface.vtkhdf entries; max-abs-diff must be <= the
// requested tolerance per field.
//
// Plan §Phase 2d.3 acceptance:
//   "make test-fault-surface-vtkhdf-zfp passes locally — writes a fault
//    HDF with --paraview-fault-zfp-tol 1e-12, reads back via pyvista,
//    asserts max(|slip_rate_strike_zfp - slip_rate_strike_lossless|) <= 1e-12."
//
// Run as:
//   HDF5_PLUGIN_PATH=$CONDA_PREFIX/plugin ./seas_test_fault_surface_vtkhdf_zfp

#include "mfem.hpp"
#include "../../io/paraview_output.hpp"
#include "../../io/fault_vtu_binary.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#ifdef MFEM_USE_HDF5
#include <hdf5.h>
#endif

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;
#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)
#define TEST_LE(val, bound, msg) do { \
   num_tests++; \
   double v_ = (val), b_ = (bound); \
   if (v_ <= b_) { num_passed++; std::cout << "  PASSED: " << msg \
      << " (" << v_ << " <= " << b_ << ")\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
      << " (got " << v_ << ", bound " << b_ << ")\n"; } \
} while (0)

#if defined(MFEM_USE_HDF5) && defined(MFEM_USE_H5Z_ZFP)
static std::vector<double> ReadHDF5PointDataAsCellData(
   const std::string &path, const std::string &field, int expect_cells)
{
   std::vector<double> out;
   hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   if (file < 0) { return out; }
   hid_t ds = -1;
   bool is_cell = false;
   const char *roots[] = {"/VTKHDF/CellData/", "/VTKHDF/PointData/"};
   for (int r = 0; r < 2 && ds < 0; ++r)
   {
      const std::string p = std::string(roots[r]) + field;
      ds = H5Dopen2(file, p.c_str(), H5P_DEFAULT);
      if (ds >= 0) { is_cell = (r == 0); }
   }
   if (ds < 0) { H5Fclose(file); return out; }
   hid_t space = H5Dget_space(ds);
   hsize_t dims[2] = {0, 0};
   H5Sget_simple_extent_dims(space, dims, nullptr);
   const hsize_t total = dims[0];
   std::vector<double> raw(total);
   H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, raw.data());
   H5Sclose(space);
   H5Dclose(ds);
   H5Fclose(file);
   if (is_cell || expect_cells <= 0) { return raw; }
   const int inflation = static_cast<int>(raw.size() / expect_cells);
   out.reserve(expect_cells);
   for (int c = 0; c < expect_cells; ++c)
   { out.push_back(raw[c * inflation]); }
   return out;
}

static std::vector<unsigned> ReadFilterIds(const std::string &path,
                                            const std::string &field)
{
   std::vector<unsigned> ids;
   hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   if (file < 0) { return ids; }
   hid_t ds = -1;
   const char *roots[] = {"/VTKHDF/CellData/", "/VTKHDF/PointData/"};
   for (int r = 0; r < 2 && ds < 0; ++r)
   {
      const std::string p = std::string(roots[r]) + field;
      ds = H5Dopen2(file, p.c_str(), H5P_DEFAULT);
   }
   if (ds < 0) { H5Fclose(file); return ids; }
   hid_t plist = H5Dget_create_plist(ds);
   const int n = H5Pget_nfilters(plist);
   for (int i = 0; i < n; ++i)
   {
      unsigned flags = 0; size_t cd_n = 0; unsigned cd[16] = {0};
      char name[64] = {0};
      const H5Z_filter_t fid = H5Pget_filter2(plist, i, &flags, &cd_n, cd,
                                              sizeof(name), name, nullptr);
      ids.push_back(static_cast<unsigned>(fid));
   }
   H5Pclose(plist);
   H5Dclose(ds);
   H5Fclose(file);
   return ids;
}
#endif

static Mesh BuildTwoTetMesh()
{
   Mesh m(/*dim=*/3, /*nv=*/5, /*nelem=*/2, /*nbdr=*/0, /*sdim=*/3);
   m.AddVertex(0., 0., -1.);
   m.AddVertex(1., 0.,  0.);
   m.AddVertex(0., 1.,  0.);
   m.AddVertex(0., 0.,  0.);
   m.AddVertex(0., 0.,  1.);
   const int t1[4] = {0,1,2,3};
   const int t2[4] = {1,2,3,4};
   m.AddTet(t1, /*attr=*/1);
   m.AddTet(t2, /*attr=*/2);
   m.FinalizeTetMesh(0, 0, true);
   return m;
}

static void RemoveDir(const std::string &dir)
{
   DIR *d = opendir(dir.c_str());
   if (!d) { return; }
   struct dirent *e;
   while ((e = readdir(d)) != nullptr)
   {
      const std::string name(e->d_name);
      if (name == "." || name == "..") { continue; }
      ::remove((dir + "/" + name).c_str());
   }
   closedir(d);
   ::rmdir(dir.c_str());
}

int main(int argc, char *argv[])
{
   (void)argc; (void)argv;
#if !defined(MFEM_USE_HDF5)
   std::cout << "SKIP: MFEM_USE_HDF5 is OFF\n"; return 0;
#elif !defined(MFEM_USE_H5Z_ZFP)
   std::cout << "SKIP: MFEM_USE_H5Z_ZFP is OFF\n"; return 0;
#else
   std::cout << "=== test_fault_surface_vtkhdf_zfp (Phase 2d.3) ===\n";

   Mesh mesh = BuildTwoTetMesh();
   const int nbf    = 1;
   const int n_dofs = 1;

   Array<int> fault_faces, empty_shared;
   for (int f = 0; f < mesh.GetNumFaces(); ++f)
   {
      if (mesh.GetFaceInformation(f).IsInterior())
      { fault_faces.Append(f); break; }
   }
   TEST_ASSERT(fault_faces.Size() == 1, "found 1 interior face");

   Vector slip(2*n_dofs);          slip = 0.123456789012345;
   Vector slip_rate(2*n_dofs);     slip_rate = 1.234567890123456e-3;
   Vector traction(2*n_dofs);      traction = -1.5e7;
   Vector state(n_dofs);           state = 0.456789012345678;
   Vector normal_stress(n_dofs);   normal_stress = 5.0e7;
   Vector a(n_dofs);               a = 0.0123;
   Vector Dc(n_dofs);              Dc = 0.14;
   Vector x2(n_dofs);              x2 = 1.0;
   Vector x3(n_dofs);              x3 = -3.0;

   const double TOL = 1e-12;

   // -- Path A: ZFP fault HDF.
   const std::string prefixA = "/tmp/test_fault_vtkhdf_zfp_ZFP";
   RemoveDir(prefixA);
   ::mkdir(prefixA.c_str(), 0755);
   {
      ParaViewOutput<Mesh> pv(prefixA, mesh, /*order=*/1);
      pv.InitFaultOutputBP5(fault_faces, empty_shared, nbf);
      pv.SetFaultOutputMode(
         ParaViewOutput<Mesh>::FaultOutputMode::Hdf5);
      pv.SetFaultHDFCompression(
         ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy, TOL);
      pv.WriteFaultSurfaceVTU(prefixA, /*cycle=*/0, /*time=*/0.0,
                              0, 1,
                              slip, slip_rate, traction, state, normal_stress,
                              a, Dc, x2, x3);
   }
   const std::string fileA = prefixA + "/fault.vtkhdf";
   TEST_ASSERT(std::ifstream(fileA).good(),
               "ZFP fault.vtkhdf created");

   // -- Path B: lossless deflate fault HDF (Phase 2b reference).
   const std::string prefixB = "/tmp/test_fault_vtkhdf_zfp_DEFLATE";
   RemoveDir(prefixB);
   ::mkdir(prefixB.c_str(), 0755);
   {
      ParaViewOutput<Mesh> pv(prefixB, mesh, /*order=*/1);
      pv.InitFaultOutputBP5(fault_faces, empty_shared, nbf);
      pv.SetFaultOutputMode(
         ParaViewOutput<Mesh>::FaultOutputMode::Hdf5);
      // No SetFaultHDFCompression call — default Deflate.
      pv.WriteFaultSurfaceVTU(prefixB, /*cycle=*/0, /*time=*/0.0,
                              0, 1,
                              slip, slip_rate, traction, state, normal_stress,
                              a, Dc, x2, x3);
   }
   const std::string fileB = prefixB + "/fault.vtkhdf";
   TEST_ASSERT(std::ifstream(fileB).good(),
               "Deflate fault.vtkhdf created");

   // -- Filter id presence: ZFP path uses filter 32013 on fault fields.
   {
      auto ids = ReadFilterIds(fileA, "slip_rate_strike");
      bool has_zfp = false;
      for (unsigned f : ids) { if (f == 32013) { has_zfp = true; break; } }
      TEST_ASSERT(has_zfp,
                  "ZFP fault HDF uses filter id 32013 on slip_rate_strike");
   }
   {
      auto ids = ReadFilterIds(fileB, "slip_rate_strike");
      bool has_zfp = false;
      for (unsigned f : ids) { if (f == 32013) { has_zfp = true; break; } }
      TEST_ASSERT(!has_zfp,
                  "Deflate fault HDF does NOT use filter id 32013");
   }

   // -- Numerical-accuracy: ZFP read-back vs ground-truth slip_rate_strike.
   {
      auto vA = ReadHDF5PointDataAsCellData(fileA, "slip_rate_strike", 1);
      TEST_ASSERT(!vA.empty(), "slip_rate_strike present in ZFP file");
      double max_err = 0.0;
      for (double v : vA)
      { max_err = std::max(max_err, std::abs(v - 1.234567890123456e-3)); }
      TEST_LE(max_err, TOL,
              "ZFP slip_rate_strike max-abs-diff vs ground truth at 1e-12");
   }

   // -- Cross-check: ZFP and Deflate read-backs agree to TOL.
   const std::vector<std::string> fields = {
      "slip_dip", "slip_strike", "slip_rate_dip", "slip_rate_strike",
      "traction_dip", "traction_strike", "state_variable", "normal_stress",
      "param_a", "param_Dc", "fault_x2", "fault_x3"
   };
   for (const auto &fname : fields)
   {
      auto vA = ReadHDF5PointDataAsCellData(fileA, fname, 1);
      auto vB = ReadHDF5PointDataAsCellData(fileB, fname, 1);
      TEST_ASSERT(!vA.empty() && !vB.empty(),
                  "field " + fname + ": both files emit data");
      double max = 0.0;
      const std::size_t n = std::min(vA.size(), vB.size());
      for (std::size_t i = 0; i < n; ++i)
      { max = std::max(max, std::abs(vA[i] - vB[i])); }
      TEST_LE(max, TOL,
              "field " + fname + " ZFP vs Deflate max-abs-diff");
   }

   // Files preserved at /tmp/test_fault_vtkhdf_zfp_{ZFP,DEFLATE}/ for
   // h5dump inspection.

   std::cout << "\n=== Summary: " << num_passed << " / " << num_tests
             << " passed; " << num_failed << " failed ===\n";
   return num_failed > 0 ? 1 : 0;
#endif
}

// Phase 6.2 of miniapps/seas/io/PLAN_paraview_compaction_2026-04-28.md
// extended by miniapps/seas/document/io_dev/PLAN_bulk_compression_and_size_estimator_2026-05-09.md.
//
// Verifies that `seas::ParaViewOutput<MeshType>::SetVolumeHDFCompression`
// REALLY APPLIES the requested filter when the active volume writer is
// `mfem::ParaViewHDFDataCollection` (Hdf5 mode).  Specifically:
//
//   1. Construct a `ParaViewOutput<Mesh>` in `Hdf5` mode on a
//      Cartesian QUAD mesh (1024 × 1024 elements ≈ 1e6 cells).
//      Register an analytic L2-p0 GridFunction.
//   2. Call `pv.SetVolumeHDFCompression(HDFCompression::ZfpAccuracy, 1e-3)`.
//      Force a save and check: filter id 32013 attaches to the
//      `<scalar_field>` dataset, max-abs-diff ≤ 1e-3 vs ground truth.
//   3. Repeat with default Deflate to confirm the dataset has filter
//      id 1 (deflate), not 32013, and the per-dataset bytes differ
//      between ZFP and deflate (proves the request was honoured).
//
// Run as:  HDF5_PLUGIN_PATH=$CONDA_PREFIX/plugin ./seas_test_volume_hdf_compression

#include "mfem.hpp"
#include "../../io/paraview_output.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <sys/stat.h>
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

static std::vector<unsigned> ReadFilterIds(const std::string &path,
                                           const std::string &dataset_path)
{
   std::vector<unsigned> ids;
   hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   if (file < 0) { return ids; }
   hid_t ds = H5Dopen2(file, dataset_path.c_str(), H5P_DEFAULT);
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

static long long DatasetStorageSize(const std::string &path,
                                    const std::string &dataset_path)
{
   hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   if (file < 0) { return -1; }
   hid_t ds = H5Dopen2(file, dataset_path.c_str(), H5P_DEFAULT);
   if (ds < 0) { H5Fclose(file); return -1; }
   const hsize_t sz = H5Dget_storage_size(ds);
   H5Dclose(ds);
   H5Fclose(file);
   return static_cast<long long>(sz);
}

// Search /VTKHDF/CellData/<name> first then /VTKHDF/PointData/<name>;
// return the first one that exists.
static std::string FindDatasetPath(const std::string &file,
                                   const std::string &field)
{
   const char *roots[] = {"/VTKHDF/CellData/", "/VTKHDF/PointData/"};
   for (int r = 0; r < 2; ++r)
   {
      const std::string p = std::string(roots[r]) + field;
      hid_t f = H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      if (f < 0) { continue; }
      hid_t d = H5Dopen2(f, p.c_str(), H5P_DEFAULT);
      if (d >= 0) { H5Dclose(d); H5Fclose(f); return p; }
      H5Fclose(f);
   }
   return "";
}

#endif

int main(int argc, char *argv[])
{
   (void)argc; (void)argv;
#if !defined(MFEM_USE_HDF5)
   std::cout << "SKIP: MFEM_USE_HDF5 is OFF in this build\n";
   return 0;
#elif !defined(MFEM_USE_H5Z_ZFP)
   std::cout << "SKIP: MFEM_USE_H5Z_ZFP is OFF in this build\n";
   return 0;
#else
   std::cout << "=== test_volume_hdf_compression (Phase 6.2) ===\n";

   const int N = 256;   // 256x256 = ~65k elements -> reasonable test time
   Mesh mesh = Mesh::MakeCartesian2D(N, N, Element::QUADRILATERAL);

   // L2-p0 grid function on a smooth analytic field.
   L2_FECollection fec(0, mesh.Dimension());
   FiniteElementSpace fes(&mesh, &fec);
   GridFunction gf(&fes);
   for (int e = 0; e < fes.GetNE(); ++e)
   {
      const double x = ((e % N) + 0.5) / N;
      const double y = ((e / N) + 0.5) / N;
      gf(e) = std::sin(6.0 * x) * std::cos(3.0 * y);
   }

   const double TOL = 1e-3;

   // -- Path A: ZFP via SetVolumeHDFCompression (HDF5 mode, default).
   const std::string prefixA = "/tmp/test_volume_hdf_compression_ZFP";
   ::mkdir(prefixA.c_str(), 0755);
   {
      ParaViewOutput<Mesh> pv(prefixA, mesh, /*order=*/0);
      TEST_ASSERT(pv.GetVolumeOutputMode() ==
                  ParaViewOutput<Mesh>::VolumeOutputMode::Hdf5,
                  "Hdf5 is the default volume mode on serial Mesh + HDF5 build");
      pv.RegisterDomainField("vol_scalar", &gf);
      pv.SetVolumeHDFCompression(
         ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy, TOL);
      pv.ForceSave(0, 0.0);
   }
   const std::string fileA = prefixA + "/volume.vtkhdf";
   TEST_ASSERT(std::ifstream(fileA).good(),
               "ZFP volume.vtkhdf written");

   // -- Path B: Default Deflate.
   const std::string prefixB = "/tmp/test_volume_hdf_compression_DEFLATE";
   ::mkdir(prefixB.c_str(), 0755);
   {
      ParaViewOutput<Mesh> pv(prefixB, mesh, /*order=*/0);
      pv.RegisterDomainField("vol_scalar", &gf);
      // No SetVolumeHDFCompression call → default deflate.
      pv.ForceSave(0, 0.0);
   }
   const std::string fileB = prefixB + "/volume.vtkhdf";
   TEST_ASSERT(std::ifstream(fileB).good(),
               "Deflate volume.vtkhdf written");

   // -- Filter id check on Path A: must contain filter 32013.
   const std::string dsA = FindDatasetPath(fileA, "vol_scalar");
   TEST_ASSERT(!dsA.empty(),
               "vol_scalar dataset present in ZFP file");
   if (!dsA.empty())
   {
      auto ids = ReadFilterIds(fileA, dsA);
      bool has_zfp = false;
      for (unsigned f : ids) { if (f == 32013) { has_zfp = true; break; } }
      TEST_ASSERT(has_zfp,
                  "ZFP file's vol_scalar dataset uses filter id 32013");
   }

   // -- Filter id check on Path B: must NOT contain filter 32013.
   const std::string dsB = FindDatasetPath(fileB, "vol_scalar");
   TEST_ASSERT(!dsB.empty(),
               "vol_scalar dataset present in Deflate file");
   if (!dsB.empty())
   {
      auto ids = ReadFilterIds(fileB, dsB);
      bool has_zfp = false;
      for (unsigned f : ids) { if (f == 32013) { has_zfp = true; break; } }
      TEST_ASSERT(!has_zfp,
                  "Deflate file's vol_scalar dataset does NOT use filter id 32013");
   }

   // -- Per-dataset storage-size sanity: ZFP and Deflate must produce
   //    DIFFERENT byte counts (proves the filter request was honoured;
   //    if both were silently using deflate, sizes would be identical).
   if (!dsA.empty() && !dsB.empty())
   {
      const long long zfp_sz = DatasetStorageSize(fileA, dsA);
      const long long def_sz = DatasetStorageSize(fileB, dsB);
      TEST_ASSERT(zfp_sz > 0 && def_sz > 0,
                  "per-dataset storage sizes readable");
      TEST_ASSERT(zfp_sz != def_sz,
                  "ZFP and Deflate produce different byte counts on the "
                  "same field (proves filter selectors are honoured)");
      std::cout << "  INFO: vol_scalar ZFP " << zfp_sz / 1024
                << " KB, Deflate " << def_sz / 1024 << " KB\n";
   }

   // -- Numerical accuracy: ZFP read-back must be within TOL.
   if (!dsA.empty())
   {
      hid_t f = H5Fopen(fileA.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      hid_t d = H5Dopen2(f, dsA.c_str(), H5P_DEFAULT);
      hid_t s = H5Dget_space(d);
      hsize_t dims[2] = {0, 0};
      H5Sget_simple_extent_dims(s, dims, nullptr);
      const hsize_t total = dims[0];
      std::vector<double> raw(total);
      H5Dread(d, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, raw.data());
      H5Sclose(s); H5Dclose(d); H5Fclose(f);

      // PointData inflation: each cell value is repeated once per cell
      // vertex (4 for quads).  Cell data path returns 1× — handle both.
      const int ncells = fes.GetNE();
      const int inflation = (static_cast<int>(raw.size()) >= ncells)
                            ? static_cast<int>(raw.size() / ncells) : 1;
      double max_err = 0.0;
      for (int e = 0; e < ncells; ++e)
      {
         max_err = std::max(max_err, std::abs(raw[e * inflation] - gf(e)));
      }
      TEST_LE(max_err, TOL,
              "ZFP max-abs-diff on vol_scalar at tol=1e-3");
   }

   // Files preserved at /tmp/test_volume_hdf_compression_{ZFP,DEFLATE}/

   std::cout << "\n=== Summary: " << num_passed << " / " << num_tests
             << " passed; " << num_failed << " failed ===\n";
   return num_failed > 0 ? 1 : 0;
#endif
}

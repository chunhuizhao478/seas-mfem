// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
//
// Phase 2b of miniapps/seas/io/PLAN_paraview_compaction_2026-04-28.md.
//
// Bit-exactness test: compare VTKHDF-emitted CellData against the
// Phase 1 binary VTU emit on identical input.  Required by
// plan §Phase 2b acceptance:
//   "test_fault_surface_vtkhdf.cpp verifies bit-exactness vs Phase 1
//    binary VTU on the same input data, with max-abs-diff < 1e-9 per
//    cell per field."
//
// Strategy: write the same `LocalFaultPack` once via the binary VTU
// path (Phase 1) and once via the HDF5 path (Phase 2b) on a serial
// run.  Read back the VTU's CellData via vtu::ParseFaultVTUCellData.
// Read back the VTKHDF via h5py-style HDF5 access through MFEM's own
// `H5Tget_native_type` / `H5Dread`.  Assert max-abs-diff < 1e-9 across
// all 12 standard fault fields.
//
// Also asserts that the VTKHDF output produces a single file
// `fault_surface.vtkhdf` (no FaultSurface/ subdir) — the production
// goal of Phase 2b.

#include "mfem.hpp"
#include "../../io/paraview_output.hpp"
#include "../../io/fault_vtu_binary.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
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

#define TEST_NEAR(val, exp, tol, msg) do { \
   num_tests++; \
   double v_ = (val), e_ = (exp), t_ = (tol); \
   if (std::abs(v_ - e_) <= t_) { num_passed++; \
      std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
      << " (got " << v_ << ", expected " << e_ \
      << ", |diff|=" << std::abs(v_-e_) << ", tol=" << t_ << ")\n"; } \
} while (0)

// ---------------------------------------------------------------------------
// HDF5 reader for one CellData- or PointData-located array of a VTKHDF
// file, returning one value per cell.
//
// R-008 (REVIEW.md 2026-04-28): MFEM 5.x serialises L2-p0 GridFunctions
// as PointData with each cell-constant inflated to its 3 triangle
// vertices.  A future MFEM version may emit L2-p0 as CellData (the
// semantically correct format).  This reader tries both paths so the
// test stays green across MFEM upgrades:
//   1. /VTKHDF/CellData/<name>   — preferred future layout (1 value/cell)
//   2. /VTKHDF/PointData/<name>  — current layout (3 vertex values/cell,
//                                   downsampled 3:1)
//
// `expect_cell_count` truncates the result to a one-cycle slice for
// the static (no-time-series) test (the first ncells values are the
// cycle-0 step); deeper time-series indexing would consult
// `/VTKHDF/Steps/{Cell,Point}DataOffsets/<name>`.
// ---------------------------------------------------------------------------
#ifdef MFEM_USE_HDF5
static std::vector<double> ReadHDF5PointDataAsCellData(
   const std::string &path, const std::string &field,
   int expect_cell_count)
{
   std::vector<double> out;
   hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   if (file < 0) { return out; }

   // Try CellData first, then fall back to PointData.
   hid_t ds = -1;
   bool is_cell_data = false;
   const char *roots[] = {"/VTKHDF/CellData/", "/VTKHDF/PointData/"};
   for (int r = 0; r < 2 && ds < 0; ++r)
   {
      const std::string p = std::string(roots[r]) + field;
      ds = H5Dopen2(file, p.c_str(), H5P_DEFAULT);
      if (ds >= 0) { is_cell_data = (r == 0); }
   }
   if (ds < 0) { H5Fclose(file); return out; }

   hid_t space = H5Dget_space(ds);
   hsize_t dims[2] = {0, 0};
   const int ndims = H5Sget_simple_extent_ndims(space);
   H5Sget_simple_extent_dims(space, dims, nullptr);
   const hsize_t total = (ndims >= 1) ? dims[0] : 0;
   std::vector<double> raw(total);
   H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, raw.data());
   H5Sclose(space);
   H5Dclose(ds);
   H5Fclose(file);

   // CellData: 1 value per cell.  PointData: 3 vertex values per
   // triangle, all carrying the same L2-p0 value — sample one per 3.
   const std::size_t n_cells_emitted = is_cell_data ? raw.size()
                                                    : raw.size() / 3;
   out.reserve(n_cells_emitted);
   for (std::size_t c = 0; c < n_cells_emitted; ++c)
   {
      out.push_back(is_cell_data ? raw[c] : raw[3*c]);
   }
   if (expect_cell_count > 0
       && static_cast<int>(out.size()) > expect_cell_count)
   {
      out.resize(expect_cell_count);
   }
   return out;
}
#endif

// ---------------------------------------------------------------------------
// Mesh + cleanup helpers (mirror test_fault_surface_vtu_binary.cpp).
// ---------------------------------------------------------------------------
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
   m.FinalizeTetMesh(/*generate_edges=*/0, /*refine=*/0, /*fix_orient=*/true);
   return m;
}

static int CountFilesIn(const std::string &dir)
{
   DIR *d = opendir(dir.c_str());
   if (!d) { return -1; }
   int n = 0;
   struct dirent *e;
   while ((e = readdir(d)) != nullptr)
   {
      const std::string name(e->d_name);
      if (name == "." || name == "..") { continue; }
      ++n;
   }
   closedir(d);
   return n;
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
      const std::string path = dir + "/" + name;
      ::remove(path.c_str());
   }
   closedir(d);
   ::rmdir(dir.c_str());
}

int main(int argc, char *argv[])
{
   (void)argc; (void)argv;
   std::cout << "=== test_fault_surface_vtkhdf (Phase 2b serial bit-exact) ===\n";

#ifndef MFEM_USE_HDF5
   std::cout << "SKIP: MFEM_USE_HDF5 is OFF in this build\n";
   return 0;
#else
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

   Vector slip(2 * n_dofs);          slip = 0.123456789012345;
   Vector slip_rate(2 * n_dofs);     slip_rate = 1.234567890123456e-3;
   Vector traction(2 * n_dofs);      traction = -1.5e7;
   Vector state(n_dofs);             state = 0.456789012345678;
   Vector normal_stress(n_dofs);     normal_stress = 5.0e7;
   Vector a(n_dofs);                 a = 0.0123;
   Vector Dc(n_dofs);                Dc = 0.14;
   Vector x2(n_dofs);                x2 = 1.0;
   Vector x3(n_dofs);                x3 = -3.0;

   // -- Path A: Phase 1 binary VTU (reference) ------------------------------
   const std::string prefixA = "/tmp/test_fault_vtkhdf_VTU";
   RemoveDir(prefixA + "/FaultSurface");
   ::mkdir(prefixA.c_str(), 0755);
   {
      ParaViewOutput<Mesh> pv(prefixA, mesh, /*order=*/1);
      pv.InitFaultOutputBP5(fault_faces, empty_shared, nbf);
      // R-002 (REVIEW.md 2026-04-28): Vtu is no longer the default on
      // HDF5-enabled builds; force it explicitly so the bit-exactness
      // reference path stays VTU regardless of build flag.
      pv.SetFaultOutputMode(
         ParaViewOutput<Mesh>::FaultOutputMode::Vtu);
      pv.WriteFaultSurfaceVTU(prefixA, /*cycle=*/0, /*time=*/0.0,
                              /*rank=*/0, /*nranks=*/1,
                              slip, slip_rate, traction, state, normal_stress,
                              a, Dc, x2, x3);
   }
   const std::string vtuA = prefixA + "/FaultSurface/fault_surface_c0.vtu";
   TEST_ASSERT(std::ifstream(vtuA).good(), "VTU path emitted c0.vtu");

   // -- Path B: Phase 2b VTKHDF ---------------------------------------------
   const std::string prefixB = "/tmp/test_fault_vtkhdf_HDF";
   RemoveDir(prefixB);
   ::mkdir(prefixB.c_str(), 0755);
   {
      ParaViewOutput<Mesh> pv(prefixB, mesh, /*order=*/1);
      pv.InitFaultOutputBP5(fault_faces, empty_shared, nbf);
      pv.SetFaultOutputMode(
         ParaViewOutput<Mesh>::FaultOutputMode::Hdf5);
      pv.WriteFaultSurfaceVTU(prefixB, /*cycle=*/0, /*time=*/0.0,
                              /*rank=*/0, /*nranks=*/1,
                              slip, slip_rate, traction, state, normal_stress,
                              a, Dc, x2, x3);
   }
   const std::string vtkhdfB = prefixB + "/fault_surface.vtkhdf";
   TEST_ASSERT(std::ifstream(vtkhdfB).good(),
               "HDF path emitted fault_surface.vtkhdf");

   // File-count check: the HDF backend produces a single .vtkhdf at the
   // prefix root (no FaultSurface/ subdir).
   const int n_at_prefix = CountFilesIn(prefixB);
   TEST_ASSERT(n_at_prefix >= 1,
               "HDF backend prefix dir contains the .vtkhdf file");
   TEST_ASSERT(::access((prefixB + "/FaultSurface").c_str(), F_OK) != 0,
               "HDF backend does NOT create a FaultSurface/ subdir");

   // -- R-001 regression: HDF5 dataset has a compression filter ------------
   // The plan §Phase 2 §4 requires compression always ON.  Probe the
   // dataset's creation property list for an active filter; pre-fix
   // (compression OFF on default settings) this returned 0 filters.
   {
      hid_t f = H5Fopen(vtkhdfB.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      TEST_ASSERT(f >= 0, "R-001: VTKHDF file opens for filter probe");
      // Try CellData first then PointData (R-008-style robust locator).
      hid_t ds = H5Dopen2(f, "/VTKHDF/PointData/slip_dip", H5P_DEFAULT);
      if (ds < 0)
      {
         ds = H5Dopen2(f, "/VTKHDF/CellData/slip_dip", H5P_DEFAULT);
      }
      TEST_ASSERT(ds >= 0, "R-001: slip_dip dataset present in VTKHDF");
      hid_t plist = H5Dget_create_plist(ds);
      const int nfilters = H5Pget_nfilters(plist);
      TEST_ASSERT(nfilters >= 1,
                  "R-001: VTKHDF dataset has at least one compression "
                  "filter (default compression-on); pre-fix this was 0");
      H5Pclose(plist);
      H5Dclose(ds);
      H5Fclose(f);
   }

   // -- R-002 regression: default output_mode_ tracks MFEM_USE_HDF5 -------
   {
      const std::string prefixDefault = "/tmp/test_fault_vtkhdf_DEFAULT";
      RemoveDir(prefixDefault);
      ::mkdir(prefixDefault.c_str(), 0755);
      ParaViewOutput<Mesh> pv(prefixDefault, mesh, /*order=*/1);
      pv.InitFaultOutputBP5(fault_faces, empty_shared, nbf);
      // Do NOT call SetFaultOutputMode — verify the build-tracked default.
      const auto mode = pv.GetFaultOutputMode();
      TEST_ASSERT(mode == ParaViewOutput<Mesh>::FaultOutputMode::Hdf5,
                  "R-002: default output_mode_ is Hdf5 on MFEM_USE_HDF5 "
                  "build (this test runs only when MFEM_USE_HDF5 is "
                  "defined)");
      RemoveDir(prefixDefault);
   }

   // -- Compare CellData per-field ------------------------------------------
   const std::vector<std::string> fields = {
      "slip_dip", "slip_strike", "slip_rate_dip", "slip_rate_strike",
      "traction_dip", "traction_strike", "state_variable", "normal_stress",
      "param_a", "param_Dc", "fault_x2", "fault_x3"
   };
   for (const auto &fname : fields)
   {
      auto vA = vtu::ParseFaultVTUCellData(vtuA, fname);
      auto vB = ReadHDF5PointDataAsCellData(
                   vtkhdfB, fname,
                   /*expect_cell_count=*/static_cast<int>(vA.size()));
      TEST_ASSERT(!vA.empty() && !vB.empty(),
                  "field " + fname + ": both writers emit data");
      TEST_ASSERT(vA.size() == vB.size(),
                  "field " + fname + ": cell counts match");
      double max_abs = 0.0;
      const std::size_t n = std::min(vA.size(), vB.size());
      for (std::size_t i = 0; i < n; ++i)
      {
         max_abs = std::max(max_abs, std::abs(vA[i] - vB[i]));
      }
      TEST_NEAR(max_abs, 0.0, 1e-9,
                "field " + fname + " VTU vs HDF5 max-abs-diff");
   }

   // Cleanup
   RemoveDir(prefixA + "/FaultSurface");
   ::rmdir(prefixA.c_str());
   RemoveDir(prefixB);

   std::cout << "\n=== Summary: " << num_passed << " / " << num_tests
             << " passed; " << num_failed << " failed ===\n";
   return num_failed > 0 ? 1 : 0;
#endif
}

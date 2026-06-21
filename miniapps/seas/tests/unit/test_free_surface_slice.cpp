// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_free_surface_slice.cpp — Phase 4 of
// document/free_surface_slice_dev/PLAN_free_surface_slice.md.
//
// Class-level coverage for seas::FreeSurfaceOutput (io/free_surface_output.hpp).
// Run with `mpirun -np 1` AND `mpirun -np 2 seas_test_free_surface_slice`
// (the np=2 case exercises the empty-local-submesh ranks and the collective
// CreateFromBoundary / ParTransferMap::Transfer / DataCollection::Save paths).
//
// Builds a 2x2x2 unit-cube hex mesh, tags the z==z_max boundary faces with a
// distinct attribute, and asserts:
//   1. the slice is codim-1 (Dimension()==2) with GlobalNE()==4 (the four top
//      quads, reduced across ranks);
//   2. after UpdateVelocity with a parent field v(x)=(x0,x1,x2), the sliced
//      velocity equals the analytic trace to < 1e-12 at every sub node
//      (GaussLobatto nodal coincidence; linear field is interpolated exactly);
//   3. the fixed-dt schedule (ShouldWrite) honours the 0.99 tolerance;
//   4. Save writes ONE consolidated fs_c<cycle>.vtu + fs.pvd (like the fault
//      output) — content-checked (8 triangles, 3-component velocity, PVD
//      reference) with NO per-rank proc*.vtu;
//   5. the standalone fsvtu:: writer emits correct point/cell counts + values.

#include "mfem.hpp"
#include "../../io/free_surface_output.hpp"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

using namespace mfem;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; \
      std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while (0)

#ifdef MFEM_USE_MPI

namespace
{

// Build a 2x2x2 unit-cube hex mesh and tag the z==z_max boundary faces with
// attribute `kTop` (all other boundary faces get attribute 1).
Mesh MakeTaggedBox(int kTop)
{
   Mesh mesh = Mesh::MakeCartesian3D(2, 2, 2, Element::HEXAHEDRON);
   Array<int> v;
   for (int be = 0; be < mesh.GetNBE(); ++be)
   {
      mesh.GetBdrElementVertices(be, v);
      real_t zc = 0.0;
      for (int k = 0; k < v.Size(); ++k) { zc += mesh.GetVertex(v[k])[2]; }
      zc /= static_cast<real_t>(v.Size());
      mesh.SetBdrAttribute(be, (std::abs(zc - 1.0) < 1e-9) ? kTop : 1);
   }
   mesh.SetAttributes();
   return mesh;
}

// Read an entire file into a string (binary-safe; we only grep the ASCII XML
// header of the VTU, which precedes any appended binary block).
std::string Slurp(const std::string &path)
{
   std::ifstream f(path, std::ios::binary);
   std::stringstream ss; ss << f.rdbuf();
   return ss.str();
}

// Parse an integer XML attribute value, e.g. AttrInt(s, "NumberOfCells").
int AttrInt(const std::string &s, const std::string &key)
{
   const std::string m = key + "=\"";
   const auto p = s.find(m);
   if (p == std::string::npos) { return -1; }
   const auto q = s.find('"', p + m.size());
   if (q == std::string::npos) { return -1; }
   return std::stoi(s.substr(p + m.size(), q - (p + m.size())));
}

}  // namespace

int main(int argc, char *argv[])
{
   MPI_Init(&argc, &argv);
   int rank = 0, nranks = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nranks);

   if (rank == 0)
   {
      std::cout << "Running Phase 4 test_free_surface_slice (np=" << nranks
                << ")\n";
   }

   const int    kTop  = 7;
   const int    order = 2;
   const real_t dt    = 0.05;
   const std::string tmpdir = "/tmp/test_free_surface_slice";

   // Fresh output dir (rank 0 creates; all ranks wait).
   if (rank == 0)
   {
      std::error_code ec;
      std::filesystem::remove_all(tmpdir, ec);
      std::filesystem::create_directories(tmpdir, ec);
   }
   MPI_Barrier(MPI_COMM_WORLD);

   // --- 0. Direct writer (serial pack -> VTU -> value check) --------------
   // Two triangles with constant per-triangle velocity: verify point/cell
   // counts, the 3-component PointData velocity (with values), the CellData
   // mpi_rank, and that BINARY uses the fault-identical raw-appended UInt64.
   if (rank == 0)
   {
      seas::fsvtu::LocalSurfacePack pk;
      pk.vertices  = { {0,0,0},{1,0,0},{0,1,0}, {1,1,0},{2,1,0},{1,2,0} };
      pk.velocity  = { {1,2,3},{1,2,3},{1,2,3}, {4,5,6},{4,5,6},{4,5,6} };
      pk.rank_cell = { 0.0, 0.0 };
      seas::fsvtu::GatheredSurfacePack gg =
         seas::fsvtu::GatherSurfacePackToRoot(pk, 0, 1);
      TEST_ASSERT(gg.ntri() == 2 && gg.vertices.size() == 6,
                  "serial gather identity: 2 triangles / 6 points");

      const std::string ap = tmpdir + "/direct_ascii.vtu";
      seas::fsvtu::WriteSurfacePackVTU(ap, gg, VTKFormat::ASCII);
      const std::string a = Slurp(ap);
      TEST_ASSERT(AttrInt(a, "NumberOfPoints") == 6, "ASCII VTU NumberOfPoints==6");
      TEST_ASSERT(AttrInt(a, "NumberOfCells")  == 2, "ASCII VTU NumberOfCells==2");
      TEST_ASSERT(a.find("Name=\"velocity\"") != std::string::npos &&
                  a.find("NumberOfComponents=\"3\"") != std::string::npos,
                  "velocity is a 3-component PointData array");
      TEST_ASSERT(a.find("4 5 6") != std::string::npos,
                  "tri-1 velocity value (4 5 6) written to PointData");
      TEST_ASSERT(a.find("Name=\"mpi_rank\"") != std::string::npos,
                  "mpi_rank CellData present");

      const std::string bp = tmpdir + "/direct_bin.vtu";
      seas::fsvtu::WriteSurfacePackVTU(bp, gg, VTKFormat::BINARY);
      const std::string b = Slurp(bp);
      TEST_ASSERT(b.find("AppendedData encoding=\"raw\"") != std::string::npos &&
                  b.find("header_type=\"UInt64\"") != std::string::npos,
                  "BINARY VTU uses raw-appended UInt64 (fault-identical)");
      TEST_ASSERT(AttrInt(b, "NumberOfCells") == 2, "BINARY VTU NumberOfCells==2");
   }

   Mesh mesh = MakeTaggedBox(kTop);
   ParMesh pmesh(MPI_COMM_WORLD, mesh);

   Array<int> top_attrs(1);
   top_attrs[0] = kTop;

   seas::FreeSurfaceOutput fs(tmpdir, pmesh, top_attrs, order, rank, dt,
                              "fs", seas::FreeSurfaceOutput::Mode::Vtu);

   // --- 1. Geometry: codim-1 slice with 4 top quads (global) -------------
   TEST_ASSERT(fs.Dimension() == 2, "slice is codim-1 (Dimension()==2)");
   TEST_ASSERT(fs.GlobalNE() == 4,
               "GlobalNE()==4 (four top quads of a 2x2x2 box, reduced)");

   // --- 2. Trace equivalence: sub velocity == analytic v(x)=(x0,x1,x2) ----
   // Parent velocity GF on the SAME space the class builds (L2 order-2 GLL,
   // vdim=3, byNODES) so UpdateVelocity's memcpy is meaningful.  A linear
   // field is interpolated exactly on a GaussLobatto nodal basis.
   L2_FECollection      pfec(order, 3, BasisType::GaussLobatto);
   ParFiniteElementSpace pfes(&pmesh, &pfec, 3, Ordering::byNODES);
   ParGridFunction       pgf(&pfes);
   VectorFunctionCoefficient vc(3, [](const Vector &x, Vector &val)
   {
      val = x;   // v(x) = (x0, x1, x2)
   });
   pgf.ProjectCoefficient(vc);

   fs.UpdateVelocity(pgf.GetData(), pfes.GetNDofs());

   // Independent reference: project the SAME analytic field onto a fresh
   // velocity space built on the class's own submesh.  Same FEC kind / vdim /
   // ordering => identical dof layout, so a direct vector comparison is valid.
   L2_FECollection      rfec(order, fs.SubMesh().Dimension(),
                             BasisType::GaussLobatto);
   ParFiniteElementSpace rfes(&fs.SubMesh(), &rfec, 3, Ordering::byNODES);
   ParGridFunction       rgf(&rfes);
   rgf.ProjectCoefficient(vc);

   const ParGridFunction &sub = fs.SubVelocity();
   real_t local_err = 0.0;
   TEST_ASSERT(sub.Size() == rgf.Size(),
               "sub velocity dof count matches reference space");
   const int n = std::min(sub.Size(), rgf.Size());
   for (int i = 0; i < n; ++i)
   {
      local_err = std::max(local_err, std::abs(sub(i) - rgf(i)));
   }
   real_t global_err = 0.0;
   MPI_Allreduce(&local_err, &global_err, 1, MPI_DOUBLE, MPI_MAX,
                 MPI_COMM_WORLD);
   TEST_ASSERT(global_err < 1e-12,
               "sliced velocity == analytic trace to < 1e-12 "
               "(max err over all ranks)");

   // --- 3. Schedule: 0.99-tolerance fixed-dt cadence ---------------------
   TEST_ASSERT(fs.ShouldWrite(0.0),
               "ShouldWrite(0.0) true initially (last_write_time_=-1e30)");
   fs.Save(0, 0.0);
   TEST_ASSERT(!fs.ShouldWrite(0.01),
               "ShouldWrite(0.01) false after Save(0,0.0) (0.01 < 0.99*0.05)");
   TEST_ASSERT(fs.ShouldWrite(0.05),
               "ShouldWrite(0.05) true after Save(0,0.0) (0.05 >= 0.99*0.05)");

   // --- 4. Files written: a non-empty .pvd + a non-empty .vtu/.pvtu ------
   MPI_Barrier(MPI_COMM_WORLD);   // ensure Save's writes are flushed
   if (rank == 0)
   {
      bool found_pvd = false, found_vtu = false;
      std::error_code ec;
      for (const auto &p :
           std::filesystem::recursive_directory_iterator(tmpdir, ec))
      {
         if (!p.is_regular_file()) { continue; }
         const std::string ext = p.path().extension().string();
         if (ext == ".pvd"  && p.file_size() > 0) { found_pvd = true; }
         if ((ext == ".vtu" || ext == ".pvtu") && p.file_size() > 0)
         {
            found_vtu = true;
         }
      }
      TEST_ASSERT(found_pvd, "Save wrote a non-empty .pvd");
      TEST_ASSERT(found_vtu, "Save wrote a non-empty .vtu/.pvtu");

      // 4b. Single CONSOLIDATED VTU (flat fs_c0.vtu), PVD points at it, and
      //     NO per-rank proc*.vtu — the whole point of this writer.
      const std::string vtu_c0 = tmpdir + "/fs_c0.vtu";
      const std::string pvd_c0 = tmpdir + "/fs.pvd";
      TEST_ASSERT(std::filesystem::exists(vtu_c0),
                  "single consolidated fs_c0.vtu written (flat, like the fault)");
      bool any_proc = false;
      std::error_code ec2;
      for (const auto &p :
           std::filesystem::recursive_directory_iterator(tmpdir, ec2))
      {
         if (p.is_regular_file() &&
             p.path().filename().string().rfind("proc", 0) == 0)
         { any_proc = true; }
      }
      TEST_ASSERT(!any_proc, "NO per-rank proc*.vtu emitted (single-file writer)");
      if (std::filesystem::exists(vtu_c0))
      {
         const std::string vc = Slurp(vtu_c0);
         // 4 top quads x 2 triangles = 8, gathered across all ranks.
         TEST_ASSERT(AttrInt(vc, "NumberOfCells") == 8,
                     "fs_c0.vtu has 8 triangles (4 top quads x 2 tris, gathered)");
         TEST_ASSERT(vc.find("Name=\"velocity\"") != std::string::npos &&
                     vc.find("NumberOfComponents=\"3\"") != std::string::npos,
                     "fs_c0.vtu carries a 3-component velocity vector");
      }
      const std::string pc = Slurp(pvd_c0);
      TEST_ASSERT(pc.find("file=\"fs_c0.vtu\"") != std::string::npos,
                  "fs.pvd references fs_c0.vtu");
   }

   // Cleanup (rank 0 only, after all ranks finished reading/writing).
   MPI_Barrier(MPI_COMM_WORLD);
   if (rank == 0)
   {
      std::error_code ec;
      std::filesystem::remove_all(tmpdir, ec);
   }

   // Reduce pass/fail across ranks so every rank returns the same code.
   int total_failed = 0;
   MPI_Allreduce(&num_failed, &total_failed, 1, MPI_INT, MPI_SUM,
                 MPI_COMM_WORLD);
   if (rank == 0)
   {
      std::cout << "\n========================================\n";
      std::cout << "Phase 4 test_free_surface_slice: "
                << (total_failed == 0 ? "ALL PASSED" : "FAILURES")
                << " (" << total_failed << " failed across " << nranks
                << " ranks)\n";
      std::cout << "========================================\n";
   }

   MPI_Finalize();
   return (total_failed == 0) ? 0 : 1;
}

#else  // !MFEM_USE_MPI

int main(int, char **)
{
   std::cout << "test_free_surface_slice: MFEM_USE_MPI not defined — skipping.\n";
   return 0;
}

#endif // MFEM_USE_MPI

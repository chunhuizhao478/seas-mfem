// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// project_stress_to_mesh — driver that loads a gmsh .msh, loads a
// schema-v1 HDF5 stress sidecar (six components: sigma_xx, sigma_yy,
// sigma_zz, sigma_xy, sigma_yz, sigma_xz, all in Pa, compression
// POSITIVE), projects each component onto the mesh's H1(p) finite-
// element space via FieldProjector::ProjectStress, and writes a
// ParaView .pvd / .vtu pair so the user can verify the projection
// visually before launching a long SAFS simulation.
//
// Plan reference: PLAN_onfaultstress.md Phase 7 §1986-2052.
//
// Usage:
//   mpirun -np 1 seas_project_stress_to_mesh \
//       --mesh    safs_fault_box_nwcut_2000m.msh \
//       --sidecar stress_safs.h5 \
//       --out     projected_stress_2000m \
//       [--order P] \
//       [--interp {trilinear,catmull-rom}] \
//       [--ascii | --binary]
//
// Output: <out>/<out>.pvd referencing per-rank .vtu files with six
// PointData scalar fields {sigma_xx, sigma_yy, sigma_zz, sigma_xy,
// sigma_yz, sigma_xz} in Pa.

#include "mfem.hpp"

#include "../io/data_field_3d.hpp"
#include "../io/field_coefficient.hpp"
#include "../io/stress_field_3d.hpp"

#include <cstring>
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
   // Match project_velocity_to_mesh.cpp MPI init pattern.
   MPI_Init(&argc, &argv);
   int rank = 0, nprocs = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   const char* mesh_path    = "";
   const char* sidecar_path = "";
   const char* out_stem     = "projected_stress";
   int order = 2;
   bool ascii_output = false;
   const char* interp_str = "trilinear";

   mfem::OptionsParser args(argc, argv);
   args.AddOption(&mesh_path,    "-m", "--mesh",
                  "Path to the gmsh .msh file (required).");
   args.AddOption(&sidecar_path, "-s", "--sidecar",
                  "Path to a schema-v1 HDF5 stress sidecar (required).");
   args.AddOption(&out_stem,     "-o", "--out",
                  "ParaView output stem (default: projected_stress); "
                  "the resulting .pvd lives at <stem>/<stem>.pvd.");
   args.AddOption(&order, "-p", "--order",
                  "H1 finite-element order (default 2).");
   args.AddOption(&ascii_output, "-a", "--ascii", "-b", "--binary",
                  "Emit ASCII (.vtu) instead of binary (default binary).  "
                  "ASCII is needed for meshio round-trip.");
   args.AddOption(&interp_str, "-i", "--interp",
                  "Within-sidecar interpolation: 'trilinear' (default) "
                  "or 'catmull-rom' (64-voxel cardinal cubic).");
   args.Parse();
   if (!args.Good())
   {
      if (rank == 0) { args.PrintUsage(std::cout); }
      MPI_Finalize();
      return 1;
   }
   if (std::strlen(mesh_path) == 0 || std::strlen(sidecar_path) == 0)
   {
      if (rank == 0)
      {
         std::cerr << "ERROR: --mesh and --sidecar are both required.\n";
         args.PrintUsage(std::cerr);
      }
      MPI_Finalize();
      return 1;
   }

   const std::string interp_name(interp_str);
   mfem::seas::InterpMode interp_mode = mfem::seas::InterpMode::Trilinear;
   if (interp_name == "trilinear")
   {
      interp_mode = mfem::seas::InterpMode::Trilinear;
   }
   else if (interp_name == "catmull-rom" || interp_name == "catmull_rom")
   {
      interp_mode = mfem::seas::InterpMode::CatmullRom;
   }
   else
   {
      if (rank == 0)
      {
         std::cerr << "ERROR: --interp must be 'trilinear' or "
                   << "'catmull-rom'; got '" << interp_name << "'.\n";
      }
      MPI_Finalize();
      return 1;
   }

   if (rank == 0)
   {
      std::cout << "loading mesh   : " << mesh_path << "\n"
                << "loading sidecar: " << sidecar_path << "\n"
                << "FE order       : " << order << "\n"
                << "interp mode    : " << interp_name << "\n"
                << "output stem    : " << out_stem << "\n";
   }

   mfem::Mesh smesh(mesh_path, /*generate_edges=*/1);
   mfem::ParMesh pmesh(MPI_COMM_WORLD, smesh);
   smesh.Clear();

   mfem::H1_FECollection fec(order, pmesh.Dimension());
   mfem::ParFiniteElementSpace fes(&pmesh, &fec);

   const long long global_ne = pmesh.GetGlobalNE();
   if (rank == 0)
   {
      std::cout << "global ne = " << global_ne << "\n";
   }

   // R-105 pre-flight: open the sidecar once and check mesh-bbox
   // containment BEFORE touching the filesystem with
   // ParaViewDataCollection (which mkdirs <out>/ on construction).
   // Without this check, a containment failure inside ProjectStress
   // aborts AFTER the output directory has been created, leaving an
   // empty leftover on disk.
   {
      mfem::real_t mxmin = +std::numeric_limits<mfem::real_t>::infinity();
      mfem::real_t mymin = +std::numeric_limits<mfem::real_t>::infinity();
      mfem::real_t mzmin = +std::numeric_limits<mfem::real_t>::infinity();
      mfem::real_t mxmax = -std::numeric_limits<mfem::real_t>::infinity();
      mfem::real_t mymax = -std::numeric_limits<mfem::real_t>::infinity();
      mfem::real_t mzmax = -std::numeric_limits<mfem::real_t>::infinity();
      for (int v = 0; v < pmesh.GetNV(); v++)
      {
         const mfem::real_t *p = pmesh.GetVertex(v);
         mxmin = std::min(mxmin, p[0]); mxmax = std::max(mxmax, p[0]);
         mymin = std::min(mymin, p[1]); mymax = std::max(mymax, p[1]);
         mzmin = std::min(mzmin, p[2]); mzmax = std::max(mzmax, p[2]);
      }
#ifdef MFEM_USE_MPI
      mfem::real_t lo[3] = {mxmin, mymin, mzmin};
      mfem::real_t hi[3] = {mxmax, mymax, mzmax};
      MPI_Allreduce(MPI_IN_PLACE, lo, 3, mfem::MPITypeMap<mfem::real_t>::mpi_type,
                    MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, hi, 3, mfem::MPITypeMap<mfem::real_t>::mpi_type,
                    MPI_MAX, MPI_COMM_WORLD);
      mxmin = lo[0]; mymin = lo[1]; mzmin = lo[2];
      mxmax = hi[0]; mymax = hi[1]; mzmax = hi[2];
#endif
      mfem::seas::StressField3D probe(sidecar_path);
      probe.SetInterpMode(interp_mode);
      if (!probe.ContainsBBox(mxmin, mxmax, mymin, mymax, mzmin, mzmax))
      {
         if (rank == 0)
         {
            const auto &bb = probe.BBox();
            std::cerr << "ERROR: mesh bbox [" << mxmin << ", " << mxmax
                      << "] x [" << mymin << ", " << mymax
                      << "] x [" << mzmin << ", " << mzmax
                      << "] is NOT contained in sidecar bbox ["
                      << bb[0] << ", " << bb[1] << "] x ["
                      << bb[2] << ", " << bb[3] << "] x ["
                      << bb[4] << ", " << bb[5] << "].\n";
         }
         MPI_Finalize();
         return 1;
      }
   }

   // Six-component bulk Cauchy stress projection.  Pre-flight bbox
   // check above already ensured containment; ProjectStress's
   // per-component checks are redundant but kept for defence in depth.
   auto sf = mfem::seas::FieldProjector::ProjectStress(
      sidecar_path, fes, interp_mode);

   // Split out_stem into directory + basename so the
   // ParaViewDataCollection writes <dir>/<base>/<base>.pvd instead of
   // doubling the path under SetPrefixPath.
   std::string out_str(out_stem);
   std::string out_dir = ".";
   std::string out_base = out_str;
   {
      const std::size_t slash = out_str.find_last_of('/');
      if (slash != std::string::npos)
      {
         out_dir = out_str.substr(0, slash);
         out_base = out_str.substr(slash + 1);
      }
   }
   mfem::ParaViewDataCollection pv(out_base, &pmesh);
   pv.SetPrefixPath(out_dir);
   pv.SetLevelsOfDetail(order);
   pv.SetHighOrderOutput(order > 1);
   pv.SetDataFormat(ascii_output ? mfem::VTKFormat::ASCII
                                 : mfem::VTKFormat::BINARY);
   // Six PointData fields in schema-v1 canonical order
   // (xx, yy, zz, xy, yz, xz).  Same order as the sidecar /fields/
   // group, so a downstream consumer can use either source.
   pv.RegisterField("sigma_xx", sf.sigma_xx.get());
   pv.RegisterField("sigma_yy", sf.sigma_yy.get());
   pv.RegisterField("sigma_zz", sf.sigma_zz.get());
   pv.RegisterField("sigma_xy", sf.sigma_xy.get());
   pv.RegisterField("sigma_yz", sf.sigma_yz.get());
   pv.RegisterField("sigma_xz", sf.sigma_xz.get());
   pv.SetCycle(0);
   pv.SetTime(0.0);
   pv.Save();

   if (rank == 0)
   {
      std::cout << "wrote ParaView collection: " << out_dir << "/"
                << out_base << "/" << out_base << ".pvd\n"
                << "  sigma_xx: [" << sf.min_sigma_xx << ", "
                << sf.max_sigma_xx << "] Pa\n"
                << "  sigma_yy: [" << sf.min_sigma_yy << ", "
                << sf.max_sigma_yy << "] Pa\n"
                << "  sigma_zz: [" << sf.min_sigma_zz << ", "
                << sf.max_sigma_zz << "] Pa\n"
                << "  sigma_xy: [" << sf.min_sigma_xy << ", "
                << sf.max_sigma_xy << "] Pa\n"
                << "  sigma_yz: [" << sf.min_sigma_yz << ", "
                << sf.max_sigma_yz << "] Pa\n"
                << "  sigma_xz: [" << sf.min_sigma_xz << ", "
                << sf.max_sigma_xz << "] Pa\n";
   }
   MPI_Finalize();
   return 0;
}

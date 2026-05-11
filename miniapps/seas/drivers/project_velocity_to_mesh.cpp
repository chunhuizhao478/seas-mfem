// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// project_velocity_to_mesh — driver that loads a gmsh .msh, loads a
// schema-v1 HDF5 velocity sidecar, projects (Vp, Vs, density, lambda,
// mu) onto the mesh's H1(p) finite-element space via
// FieldProjector::ProjectVelocity, and writes a ParaView .pvd / .vtu
// pair so the user can verify the projection visually before committing
// to a long simulation run.
//
// Usage:
//   mpirun -np 1 seas_project_velocity_to_mesh \
//       --mesh    safs_fault_box_nwcut_2000m.msh \
//       --sidecar velocity_safs.h5 \
//       --out     projected_velocity_2000m \
//       [--order P]
//
// Output: <out>.pvd referencing per-rank .vtu files with five PointData
// scalar fields {Vp, Vs, density, lambda, mu}.

#include "mfem.hpp"

#include "../io/data_field_3d.hpp"
#include "../io/field_coefficient.hpp"

#include <cstring>
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
   // Match the existing seas-driver MPI init pattern (raw MPI_Init);
   // mfem::Mpi::Init + Hypre::Init segfaults on this MFEM build path
   // when constructing the first ParFiniteElementSpace.
   MPI_Init(&argc, &argv);
   int rank = 0, nprocs = 1;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

   const char* mesh_path    = "";
   const char* sidecar_path = "";
   const char* out_stem     = "projected_velocity";
   int order = 2;
   bool ascii_output = false;
   const char* interp_str = "trilinear";

   mfem::OptionsParser args(argc, argv);
   args.AddOption(&mesh_path,    "-m", "--mesh",
                  "Path to the gmsh .msh file (required).");
   args.AddOption(&sidecar_path, "-s", "--sidecar",
                  "Path to a schema-v1 HDF5 velocity sidecar (required).");
   args.AddOption(&out_stem,     "-o", "--out",
                  "ParaView output stem (default: projected_velocity); "
                  "the resulting .pvd lives at <stem>/<stem>.pvd.");
   args.AddOption(&order, "-p", "--order",
                  "H1 finite-element order (default 2; G-2 in "
                  "fault_zone_projection_plan_v3.  Pass --order 1 to "
                  "reproduce pre-G-2 H1-P1 output).");
   args.AddOption(&ascii_output, "-a", "--ascii", "-b", "--binary",
                  "Emit ASCII (.vtu) instead of binary (default binary).  "
                  "ASCII is needed when the downstream consumer is meshio, "
                  "which does not parse MFEM's binary VTU 2.2 dialect.");
   args.AddOption(&interp_str, "-i", "--interp",
                  "Within-sidecar interpolation: 'trilinear' (default, "
                  "8-voxel O(h^2)) or 'catmull-rom' (64-voxel "
                  "tensor-product cardinal cubic, C^1, O(h^4) on "
                  "smooth fields, no overshoot at C^1 kinks; auto-"
                  "fallback to trilinear within 1 voxel of the "
                  "sidecar bbox).");
   args.Parse();
   if (!args.Good())
   {
      if (rank == 0) { args.PrintUsage(std::cout); }
      return 1;
   }
   if (std::strlen(mesh_path) == 0 || std::strlen(sidecar_path) == 0)
   {
      if (rank == 0)
      {
         std::cerr << "ERROR: --mesh and --sidecar are both required.\n";
         args.PrintUsage(std::cerr);
      }
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

   // GetGlobalNE() is collective (does an internal MPI_Allreduce),
   // so it MUST be called on every rank — wrapping it in
   // `if (rank == 0)` causes rank 0 to call MPI_Allreduce(1 long_long)
   // alone, while other ranks proceed to ProjectVelocity and call
   // MPI_Allreduce(6 doubles).  OpenMPI matches them together and
   // returns MPI_ERR_TRUNCATE -> abort.
   const long long global_ne = pmesh.GetGlobalNE();
   if (rank == 0)
   {
      std::cout << "global ne = " << global_ne << "\n";
   }

   auto vf = mfem::seas::FieldProjector::ProjectVelocity(
      sidecar_path, fes,
      /*lambda_min_pa=*/0.0,
      /*lambda_max_pa=*/1.0e12,
      /*mu_min_pa=*/    1.0e7,
      /*mu_max_pa=*/    1.0e11,
      interp_mode);

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
   pv.RegisterField("Vp",      vf.vp.get());
   pv.RegisterField("Vs",      vf.vs.get());
   pv.RegisterField("density", vf.rho.get());
   pv.RegisterField("lambda",  vf.lambda.get());
   pv.RegisterField("mu",      vf.mu.get());
   pv.SetCycle(0);
   pv.SetTime(0.0);
   pv.Save();

   if (rank == 0)
   {
      std::cout << "wrote ParaView collection: " << out_dir << "/"
                << out_base << "/" << out_base << ".pvd\n"
                << "  Vp:      [" << vf.min_vp     << ", "
                << vf.max_vp     << "] m/s\n"
                << "  Vs:      [" << vf.min_vs     << ", "
                << vf.max_vs     << "] m/s\n"
                << "  density: [" << vf.min_rho    << ", "
                << vf.max_rho    << "] kg/m^3\n"
                << "  mu:      [" << vf.min_mu     << ", "
                << vf.max_mu     << "] Pa\n"
                << "  lambda:  [" << vf.min_lambda << ", "
                << vf.max_lambda << "] Pa\n";
   }
   MPI_Finalize();
   return 0;
}

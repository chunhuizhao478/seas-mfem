// pff_partition: split a gmsh mesh once into per-rank MFEM mesh files, so that
// pff_driver ranks read only their own part instead of every rank loading the
// whole serial mesh (14.6 M tets of the SPE cube: > 13 GB per rank).
//
//   pff_partition -m mesh.msh -np 288 -o parts_np288
//
// Writes  <out>/mesh.<rank, 6 digits>   (MeshPart::Print, read by ParMesh(comm, istream))
//         <out>/physical_names.msh      ($PhysicalNames block of the input)
//         <out>/info.txt                (np, source, sizes)
// The partition is Mesh::GeneratePartitioning(np, 1) (METIS k-way), the same
// call ParMesh(comm, serial) makes, so a run on the parts matches a run that
// loads the serial mesh on the same number of ranks.
// Serial program: run it once on a node with enough memory.

#include "mfem.hpp"
#include "io/gmsh_physical_names.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

using namespace mfem;

int main(int argc, char *argv[])
{
   const char *mesh_file = "";
   const char *out_dir = "";
   int np = 0;
   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh", "gmsh mesh file to partition");
   args.AddOption(&np, "-np", "--num-parts", "number of parts (= MPI ranks of the run)");
   args.AddOption(&out_dir, "-o", "--output", "output directory (created; must not hold parts)");
   args.Parse();
   if (!args.Good() || np < 1 || std::string(mesh_file).empty() || std::string(out_dir).empty())
   {
      args.PrintUsage(std::cout);
      return 1;
   }
   namespace fs = std::filesystem;
   const fs::path dir(out_dir);
   fs::create_directories(dir);
   MFEM_VERIFY(!fs::exists(dir / "info.txt"), "pff_partition: " << dir << " already holds a partition");

   // Physical names: copy the $PhysicalNames block verbatim.
   {
      std::ifstream in(mesh_file);
      MFEM_VERIFY(in, "pff_partition: cannot open " << mesh_file);
      std::ofstream out(dir / "physical_names.msh");
      std::string line;
      bool on = false;
      while (std::getline(in, line))
      {
         if (line.rfind("$PhysicalNames", 0) == 0) { on = true; }
         if (on) { out << line << "\n"; }
         if (line.rfind("$EndPhysicalNames", 0) == 0) { break; }
      }
      MFEM_VERIFY(on, "pff_partition: no $PhysicalNames in " << mesh_file);
   }

   // Same load as pff_driver: Mesh(file, generate_edges = 1, refine = 1).
   Mesh mesh(mesh_file, 1, 1);
   std::cout << "pff_partition: " << mesh_file << ": dim " << mesh.Dimension() << ", "
             << mesh.GetNE() << " elements, " << mesh.GetNV() << " vertices, "
             << mesh.GetNBE() << " boundary elements -> " << np << " parts" << std::endl;
   int *part = mesh.GeneratePartitioning(np, 1);
   {
      Array<int> count(np);
      count = 0;
      for (int e = 0; e < mesh.GetNE(); e++) { count[part[e]]++; }
      int mn = count[0], mx = count[0];
      for (int p = 0; p < np; p++) { mn = std::min(mn, count[p]); mx = std::max(mx, count[p]); }
      MFEM_VERIFY(mn > 0, "pff_partition: METIS left an empty part");
      std::cout << "  elements per part: min " << mn << ", max " << mx << std::endl;
   }
   MeshPartitioner partitioner(mesh, np, part);
   MeshPart mesh_part;
   for (int p = 0; p < np; p++)
   {
      partitioner.ExtractPart(p, mesh_part);
      std::ofstream omesh(dir / MakeParFilename("mesh.", p));
      omesh.precision(16);
      mesh_part.Print(omesh);
      MFEM_VERIFY(omesh.good(), "pff_partition: write failed for part " << p);
   }
   delete [] part;

   std::ofstream info(dir / "info.txt");
   info << "np " << np << "\n"
        << "source " << fs::absolute(mesh_file).string() << "\n"
        << "dim " << mesh.Dimension() << "\n"
        << "elements " << mesh.GetNE() << "\n"
        << "vertices " << mesh.GetNV() << "\n"
        << "boundary_elements " << mesh.GetNBE() << "\n"
        << "max_bdr_attribute " << mesh.bdr_attributes.Max() << "\n"
        << "partitioning GeneratePartitioning(np, 1)\n";
   std::cout << "  wrote " << np << " parts to " << dir << std::endl;
   return 0;
}

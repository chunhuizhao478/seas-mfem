// Standalone tool to load a Gmsh BP2 mesh and export VTK for ParaView.
//
// Usage:
//   gmsh -2 mesh/bp2.geo -o bp2.msh
//   ./seas_generate_mesh --mesh bp2.msh [--scale 1000] [--output bp2_mesh.vtk]
//
// Outputs a VTK file that can be opened in ParaView to verify:
// - Mesh refinement near the fault (x=0)
// - Smooth grading toward far-field
// - Boundary attributes

#include "mfem.hpp"
#include "../domain/bp2_mesh.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>

using namespace mfem;
using namespace mfem::seas;

int main(int argc, char *argv[])
{
   std::string mesh_file;
   real_t scale = 1000.0;  // km → m
   std::string output = "bp2_mesh.vtk";

   for (int i = 1; i < argc; i++)
   {
      std::string arg(argv[i]);
      if (arg == "--mesh" && i+1 < argc) { mesh_file = argv[++i]; }
      if (arg == "--scale" && i+1 < argc) { scale = std::atof(argv[++i]); }
      if (arg == "--output" && i+1 < argc) { output = argv[++i]; }
   }

   if (mesh_file.empty())
   {
      std::cerr << "Usage: " << argv[0]
                << " --mesh <file.msh> [--scale 1000] [--output bp2_mesh.vtk]\n";
      std::cerr << "\nGenerate mesh first: gmsh -2 mesh/bp2.geo -o bp2.msh\n";
      return 1;
   }

   std::cout << "Loading BP2 mesh: " << mesh_file << "\n";
   std::cout << "  Scale: " << scale << "\n";

   auto mesh = BP2MeshGenerator::LoadGmshMesh(mesh_file, scale);

   std::cout << "\nMesh statistics:\n";
   std::cout << "  Elements: " << mesh->GetNE() << "\n";
   std::cout << "  Vertices: " << mesh->GetNV() << "\n";
   std::cout << "  Element type: "
             << (mesh->GetElementType(0) == Element::TRIANGLE ? "Triangle" :
                 mesh->GetElementType(0) == Element::QUADRILATERAL ? "Quad" :
                 "Unknown") << "\n";

   // Compute min/max element sizes
   real_t h_min = 1e30, h_max = 0.0;
   for (int e = 0; e < mesh->GetNE(); e++)
   {
      real_t h = mesh->GetElementSize(e);
      h_min = std::min(h_min, h);
      h_max = std::max(h_max, h);
   }
   std::cout << "  Min element size: " << h_min << " m\n";
   std::cout << "  Max element size: " << h_max << " m\n";
   std::cout << "  Size ratio: " << h_max / h_min << "\n";

   // Count elements touching the fault (x=0)
   int n_fault_elems = 0;
   for (int e = 0; e < mesh->GetNE(); e++)
   {
      Array<int> verts;
      mesh->GetElementVertices(e, verts);
      bool touches_fault = false;
      for (int v = 0; v < verts.Size(); v++)
      {
         const real_t *coords = mesh->GetVertex(verts[v]);
         if (std::abs(coords[0]) < 1.0)  // within 1m of x=0
         {
            touches_fault = true;
            break;
         }
      }
      if (touches_fault) { n_fault_elems++; }
   }
   std::cout << "  Elements touching fault (x~0): " << n_fault_elems << "\n";

   // Print boundary attribute summary
   std::cout << "\nBoundary attributes:\n";
   int nbe = mesh->GetNBE();
   int counts[5] = {0, 0, 0, 0, 0};
   for (int be = 0; be < nbe; be++)
   {
      int attr = mesh->GetBdrAttribute(be);
      if (attr >= 1 && attr <= 4) { counts[attr]++; }
      else { counts[0]++; }
   }
   std::cout << "  1 (FARFIELD_LEFT):  " << counts[1] << "\n";
   std::cout << "  2 (FARFIELD_RIGHT): " << counts[2] << "\n";
   std::cout << "  3 (FREE_SURFACE):   " << counts[3] << "\n";
   std::cout << "  4 (BOTTOM):         " << counts[4] << "\n";
   if (counts[0] > 0)
   {
      std::cout << "  Unknown:            " << counts[0] << "\n";
   }

   // Save VTK
   BP2MeshGenerator::SaveVTK(*mesh, output);
   std::cout << "\nSaved VTK to: " << output << "\n";

   return 0;
}

// remesh_fault_cgal.cpp — Phase 3 C of
// meshing/docs/PLAN_fault_triangle_quality_2026-05-26.md.
//
// Isotropic surface remesh of a cut SAFS fault that DOES clear the short
// z=0 clip segments mmgs cannot: every open-boundary edge (z=0 trace,
// NW-cut, deep perimeter) is marked constrained, then isotropic_remeshing
// runs with protect_constraints=false + collapse_constraints=true +
// relax_constraints=true.  That collapses the short trace segments and
// slides the survivors ALONG the (planar) trace polyline, so the trace
// stays a 1-D feature in the z=0 plane while the needles disappear.
//
// I/O is OFF (ASCII double precision) to avoid the float32 STL coord drift
// that would break z=0 trace conformality (EXPLORE_cgal_corefine.md).  The
// Python wrapper (remesh_fault_stl.py --engine cgal) converts STL<->OFF and
// snaps the z=0 trace to exactly 0 after.
//
// usage: remesh_fault_cgal in.off out.off target_edge [iters]
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/boost/graph/iterator.h>

#include <cstdlib>
#include <iostream>
#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3>                       Mesh;
typedef boost::graph_traits<Mesh>::edge_descriptor          edge_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor      halfedge_descriptor;
namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char** argv) {
  if (argc < 4) {
    std::cerr << "usage: remesh_fault_cgal in.off out.off target_edge [iters]\n";
    return 1;
  }
  const std::string in = argv[1], out = argv[2];
  const double target = std::atof(argv[3]);
  const unsigned iters = (argc > 4) ? static_cast<unsigned>(std::atoi(argv[4])) : 5u;
  if (!(target > 0.0)) {
    std::cerr << "ERROR: target_edge must be > 0 (got " << target << ")\n";
    return 1;
  }

  Mesh mesh;
  if (!PMP::IO::read_polygon_mesh(in, mesh) || is_empty(mesh)) {
    std::cerr << "ERROR: cannot read polygon mesh: " << in << "\n";
    return 2;
  }
  std::cout << "read  : " << num_vertices(mesh) << " verts, "
            << num_faces(mesh) << " faces\n";

  // Constrain every open-boundary edge (the fault is an open sheet: its
  // borders are the z=0 trace + NW-cut + deep perimeter).  collapse_constraints
  // then removes the short clip segments; relax_constraints redistributes the
  // survivors along the polyline.  protect_constraints MUST be false (true
  // would lock the slivers — the documented trap).
  typedef Mesh::Property_map<edge_descriptor, bool> ECM;
  ECM ecm = mesh.add_property_map<edge_descriptor, bool>("e:constrained", false).first;
  std::size_t n_border = 0;
  for (edge_descriptor e : edges(mesh)) {
    const halfedge_descriptor h = halfedge(e, mesh);
    if (is_border(h, mesh) || is_border(opposite(h, mesh), mesh)) {
      put(ecm, e, true);
      ++n_border;
    }
  }
  std::cout << "border (constrained) edges: " << n_border << "\n";

  PMP::isotropic_remeshing(
      faces(mesh), target, mesh,
      PMP::parameters::number_of_iterations(iters)
                     .edge_is_constrained_map(ecm)
                     .protect_constraints(false)
                     .collapse_constraints(true)
                     .relax_constraints(true));

  std::cout << "remesh: " << num_vertices(mesh) << " verts, "
            << num_faces(mesh) << " faces (iters=" << iters
            << ", target=" << target << ")\n";

  if (!CGAL::IO::write_polygon_mesh(out, mesh,
                                    CGAL::parameters::stream_precision(17))) {
    std::cerr << "ERROR: cannot write: " << out << "\n";
    return 3;
  }
  std::cout << "wrote : " << out << "\n";
  return 0;
}

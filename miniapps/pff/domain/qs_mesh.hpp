// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// qs_mesh.hpp — meshes of the quasi-static benchmarks (phase R).
//
// Reproduces the MOOSE mesh generators of the RACCOON tutorials:
//
//  - GeneratedMeshGenerator (dim 2, QUAD4): MFEM's Cartesian mesh.  The
//    boundary names map to MFEM's attributes bottom 1, right 2, top 3, left 4.
//  - BoundingBoxNodeSetGenerator + construct_side_list_from_node_list: the
//    boundary sides of a named boundary whose vertices all lie in the box get
//    a new attribute (5, 6, ...), named after the node set.  It is applied on
//    the coarse mesh, so refined sides inherit it.
//  - Adaptivity with initial_steps = n and a BoxMarker / OrientedBoxMarker
//    (inside = REFINE): n passes, each refining the elements whose vertex
//    average lies in the box (BoxMarker::computeElementMarker uses
//    Elem::vertex_average and an inclusive BoundingBox test).  MFEM
//    nonconforming refinement with nc_limit = 1 gives libMesh's 2:1 face level
//    limit (MeshRefinement::face_level_mismatch_limit = 1).  For
//    mode1_brittle_fracture this reproduces the gold mesh exactly: 888
//    elements and 994 nodes at the same coordinates.
//  - StitchedMeshGenerator of two halves joined on part of y = 0 (mode II):
//    MakeSlitMesh.

#ifndef MFEM_PFF_QS_MESH_HPP
#define MFEM_PFF_QS_MESH_HPP

#include "mfem.hpp"
#include "../config/qs_config.hpp"

#include <array>
#include <cmath>
#include <vector>
#include <map>
#include <stdexcept>
#include <string>

namespace mfem
{
namespace pff
{

using BoundaryNames = std::map<std::string, int>;

inline bool InAxisBox(const RefineBox &b, const real_t *x, int dim)
{
   for (int i = 0; i < dim; i++)
   {
      if (x[i] < b.lo[i] || x[i] > b.hi[i]) { return false; }
   }
   return true;
}

/// MOOSE OrientedBoxMarker: the box of the given length, width and height
/// centred at `center`, with axes length_direction, width_direction and their
/// cross product (both directions are normalised).
inline bool InOrientedBox(const RefineBox &b, const real_t *x, int dim)
{
   real_t L[3], W[3], H[3], r[3] = {0.0, 0.0, 0.0};
   real_t nl = 0.0, nw = 0.0;
   for (int i = 0; i < 3; i++)
   {
      L[i] = b.length_direction[i];
      W[i] = b.width_direction[i];
      nl += L[i] * L[i];
      nw += W[i] * W[i];
   }
   nl = std::sqrt(nl);
   nw = std::sqrt(nw);
   for (int i = 0; i < 3; i++) { L[i] /= nl; W[i] /= nw; }
   H[0] = L[1] * W[2] - L[2] * W[1];
   H[1] = L[2] * W[0] - L[0] * W[2];
   H[2] = L[0] * W[1] - L[1] * W[0];
   for (int i = 0; i < dim; i++) { r[i] = x[i] - b.center[i]; }
   real_t pl = 0.0, pw = 0.0, ph = 0.0;
   for (int i = 0; i < 3; i++) { pl += r[i] * L[i]; pw += r[i] * W[i]; ph += r[i] * H[i]; }
   return std::abs(pl) <= 0.5 * b.length && std::abs(pw) <= 0.5 * b.width &&
          std::abs(ph) <= 0.5 * b.height;
}

inline void ElementVertexAverage(const Mesh &mesh, int e, real_t *c)
{
   Array<int> v;
   mesh.GetElementVertices(e, v);
   const int dim = mesh.SpaceDimension();
   for (int i = 0; i < dim; i++) { c[i] = 0.0; }
   for (int k = 0; k < v.Size(); k++)
   {
      const real_t *x = mesh.GetVertex(v[k]);
      for (int i = 0; i < dim; i++) { c[i] += x[i]; }
   }
   for (int i = 0; i < dim; i++) { c[i] /= v.Size(); }
}

/// `steps` passes of nonconforming refinement of the elements in any box.
inline void RefineInBoxes(Mesh &mesh, const std::vector<RefineBox> &boxes, int steps)
{
   if (steps <= 0 || boxes.empty()) { return; }
   mesh.EnsureNCMesh();
   const int dim = mesh.SpaceDimension();
   for (int s = 0; s < steps; s++)
   {
      Array<int> marked;
      for (int e = 0; e < mesh.GetNE(); e++)
      {
         real_t c[3] = {0.0, 0.0, 0.0};
         ElementVertexAverage(mesh, e, c);
         for (const RefineBox &b : boxes)
         {
            if (b.oriented ? InOrientedBox(b, c, dim) : InAxisBox(b, c, dim))
            {
               marked.Append(e);
               break;
            }
         }
      }
      mesh.GeneralRefinement(marked, 1, 1);
   }
}

/// Give the sides of `ns.boundary` with all vertices in [lo, hi] attribute
/// `attr`.  Returns the number of sides changed.
inline int ApplyNodeSet(Mesh &mesh, const NodeSetBox &ns, const BoundaryNames &names, int attr)
{
   const auto it = names.find(ns.boundary);
   if (it == names.end())
   {
      throw std::invalid_argument("node set '" + ns.name + "': unknown boundary '" +
                                  ns.boundary + "'");
   }
   const int dim = mesh.SpaceDimension();
   int count = 0;
   for (int be = 0; be < mesh.GetNBE(); be++)
   {
      if (mesh.GetBdrAttribute(be) != it->second) { continue; }
      Array<int> v;
      mesh.GetBdrElementVertices(be, v);
      bool inside = true;
      for (int k = 0; k < v.Size() && inside; k++)
      {
         const real_t *x = mesh.GetVertex(v[k]);
         for (int i = 0; i < dim; i++)
         {
            if (x[i] < ns.lo[i] || x[i] > ns.hi[i]) { inside = false; }
         }
      }
      if (inside) { mesh.SetBdrAttribute(be, attr); count++; }
   }
   mesh.SetAttributes();
   return count;
}

/// The mode II slit mesh (two GeneratedMeshGenerator halves joined by
/// StitchedMeshGenerator on the BoundingBoxNodeSets x >= slit_x of y = 0):
/// [xmin,xmax] x [0,ymax] (top_half) and [xmin,xmax] x [ymin,0]
/// (bottom_half), nx x ny QUAD4 each. The y = 0 vertices with x >= slit_x
/// are shared; those with x < slit_x are duplicated, so the halves are
/// disconnected along the slit. Boundaries: <half>_{bottom,right,top,left}
/// (boundary_name_prefix), the y = 0 sides only on the slit.
inline Mesh MakeSlitMesh(const QsCaseConfig &c, BoundaryNames &names)
{
   const int nx = c.nx, ny = c.ny;
   const real_t hx = (c.xmax - c.xmin) / nx;
   MFEM_VERIFY(c.ymin < 0.0 && c.ymax > 0.0, "slit mesh: need ymin < 0 < ymax");
   // Vertex ids: top half T(i, j), j = 0..ny at y = j ymax/ny; bottom half
   // B(i, j), j = 0..ny at y = ymin + j (-ymin)/ny, B(i, ny) on y = 0.
   const int nvx = nx + 1;
   std::vector<int> T(nvx * (ny + 1)), B(nvx * (ny + 1));
   int nv = 0;
   std::vector<std::array<real_t, 2>> X;
   for (int j = 0; j <= ny; j++)
      for (int i = 0; i <= nx; i++)
      {
         T[i + nvx * j] = nv++;
         X.push_back({c.xmin + i * hx, j * c.ymax / ny});
      }
   const real_t tol = 1e-12 * (c.xmax - c.xmin);
   for (int j = 0; j <= ny; j++)
      for (int i = 0; i <= nx; i++)
      {
         const real_t x = c.xmin + i * hx;
         if (j == ny && x >= c.slit_x - tol) { B[i + nvx * j] = T[i]; continue; }
         B[i + nvx * j] = nv++;
         X.push_back({x, c.ymin + j * (-c.ymin) / ny});
      }
   names = {{"top_half_bottom", 1}, {"top_half_right", 2}, {"top_half_top", 3},
            {"top_half_left", 4}, {"bottom_half_bottom", 5}, {"bottom_half_right", 6},
            {"bottom_half_top", 7}, {"bottom_half_left", 8}};
   int nbe = 0;
   int nslit = 0;
   for (int i = 0; i < nx; i++) { if (c.xmin + (i + 0.5) * hx < c.slit_x) { nslit++; } }
   // outer boundary 2 (2 nx + 2 ny) sides, the slit 2 nslit sides
   Mesh mesh(2, nv, 2 * nx * ny, 2 * (nx + 2 * ny) + 2 * nslit, 2);
   for (const auto &x : X) { mesh.AddVertex(x[0], x[1]); }
   auto add_half = [&](const std::vector<int> &V, int off, bool top)
   {
      for (int j = 0; j < ny; j++)
         for (int i = 0; i < nx; i++)
         {
            mesh.AddQuad(V[i + nvx * j], V[i + 1 + nvx * j], V[i + 1 + nvx * (j + 1)],
                         V[i + nvx * (j + 1)], 1);
         }
      for (int i = 0; i < nx; i++)
      {
         // bottom side (y = 0 for the top half: slit part only)
         const real_t xm = c.xmin + (i + 0.5) * hx;
         if (!top || xm < c.slit_x) { mesh.AddBdrSegment(V[i + 1], V[i], off + 1); nbe++; }
         // top side (y = 0 for the bottom half: slit part only)
         if (top || xm < c.slit_x)
         {
            mesh.AddBdrSegment(V[i + nvx * ny], V[i + 1 + nvx * ny], off + 3);
            nbe++;
         }
      }
      for (int j = 0; j < ny; j++)
      {
         mesh.AddBdrSegment(V[nx + nvx * j], V[nx + nvx * (j + 1)], off + 2);
         mesh.AddBdrSegment(V[nvx * (j + 1)], V[nvx * j], off + 4);
         nbe += 2;
      }
   };
   add_half(T, 0, true);
   add_half(B, 4, false);
   mesh.FinalizeQuadMesh(1, 1, true);
   return mesh;
}

/// The serial mesh of a case, refined, and the boundary name map.
inline Mesh BuildQsMesh(const QsCaseConfig &c, BoundaryNames &names)
{
   if (c.mesh_type == "slit")
   {
      Mesh mesh = MakeSlitMesh(c, names);
      int next = names.size() + 1;
      for (const NodeSetBox &ns : c.nodesets)
      {
         if (names.count(ns.name))
         {
            throw std::invalid_argument("node set name '" + ns.name + "' is already a boundary");
         }
         if (ApplyNodeSet(mesh, ns, names, next) == 0)
         {
            throw std::invalid_argument("node set '" + ns.name + "' selects no boundary side");
         }
         names[ns.name] = next++;
      }
      RefineInBoxes(mesh, c.refine, c.refine_steps);
      mesh.Finalize();
      return mesh;
   }
   Mesh mesh = Mesh::MakeCartesian2D(c.nx, c.ny, Element::QUADRILATERAL, true,
                                     c.xmax - c.xmin, c.ymax - c.ymin, false);
   for (int i = 0; i < mesh.GetNV(); i++)
   {
      real_t *x = mesh.GetVertex(i);
      x[0] += c.xmin;
      x[1] += c.ymin;
   }
   names = {{"bottom", 1}, {"right", 2}, {"top", 3}, {"left", 4}};
   int next = 5;
   for (const NodeSetBox &ns : c.nodesets)
   {
      if (names.count(ns.name))
      {
         throw std::invalid_argument("node set name '" + ns.name + "' is already a boundary");
      }
      if (ApplyNodeSet(mesh, ns, names, next) == 0)
      {
         throw std::invalid_argument("node set '" + ns.name + "' selects no boundary side");
      }
      names[ns.name] = next++;
   }
   RefineInBoxes(mesh, c.refine, c.refine_steps);
   mesh.Finalize();
   return mesh;
}

/// Attribute numbers of a list of boundary names (or numbers as strings).
inline Array<int> BoundaryMarker(const BoundaryNames &names, const std::vector<std::string> &list,
                                 int max_attr)
{
   Array<int> marker(max_attr);
   marker = 0;
   for (const std::string &s : list)
   {
      int a = -1;
      const auto it = names.find(s);
      if (it != names.end()) { a = it->second; }
      else
      {
         try { a = std::stoi(s); }
         catch (...) { throw std::invalid_argument("unknown boundary '" + s + "'"); }
      }
      if (a < 1 || a > max_attr) { throw std::invalid_argument("boundary '" + s + "' out of range"); }
      marker[a - 1] = 1;
   }
   return marker;
}

} // namespace pff
} // namespace mfem

#endif

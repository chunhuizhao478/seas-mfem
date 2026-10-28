// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_absorbing_all.cpp — [cycle.absorbing] faces = "all" (Chunhui Zhao,
// 2026-10-27: the UCERF 2-D fault traces sit in a box that QD loads through
// its walls, and in the dynamic regime every wall of the box absorbs; there is
// no free surface, unlike BP3).  A unit test of RetagOuterFaces
// (cycle/absorbing_faces.hpp) with keep_free_surface = false and of
// AbsorbingConfig; the dynamic regime's frozen wall load and first
// acceleration with "all" are checked by test_2d_absorbing.cpp.
//
// Each case re-tags copies of one ParMesh (ParMesh(pm, true): the same local
// elements and boundary elements in the same order), so the copies compare
// boundary element by boundary element:
//
//   config       AbsorbingConfig: "all" is on, has walls and keeps no free
//                surface; "outer", "bottom" and "off" as before.
//   map_view     The TPV205-2D 100 m mesh (a box in the (x, y) plane, the
//                interior fault 103, the outer boundary 105; its frame's up
//                is e_z, which "outer" refuses in 2-D).  With 105 as the QD
//                loading walls (Dirichlet), "all" makes every exterior face
//                a wall; with 105 traction free, every exterior face absorbs
//                with no load.  Both with up = e_z and with up = 0 (not
//                read): the same attributes.  The fault faces, the element
//                attributes and the number of boundary elements unchanged.
//   bp3          The BP3 d60 1 km mesh (the surface and the bottom traction
//                free, 1; the sides and the interior creep segment
//                Dirichlet, 5; up = [0, 1]): "outer" moves the bottom and the
//                sides and keeps the surface; "all" moves the same faces to
//                the same attributes and the surface faces too (all at
//                y = 0, n . up = 1), and nothing else.  The creep segment
//                (interior) keeps 5 in both.
//   box_3d       A 2 x 2 x 2 tetrahedral box (MFEM's boundary attributes 1-6,
//                the top z = 1 is 6): the top traction free, the rest walls,
//                up = e_z.  "outer" keeps the top (8 triangles) and makes 40
//                walls; "all" absorbs the top too.
//
// Deck mode: --load-deck d.toml loads a deck (CycleConfig::Load) and prints
// its [cycle.absorbing]; the make target uses it for the loader's refusals.
//
// Usage (from miniapps/seas):
//   mpirun -np N ./seas_test_absorbing_all [--load-deck d.toml]
// Exit code 0 = PASSED.

#include "mfem.hpp"
#include "../../cycle/cycle_config.hpp"
#include "../../cycle/absorbing_faces.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace
{

int g_failures = 0;
int g_rank = 0;

void Check(bool ok, const std::string &name, const std::string &detail = "")
{
   if (!ok) { g_failures++; }
   if (g_rank == 0)
   {
      std::cout << "  " << (ok ? "PASSED" : "FAILED") << ": " << name
                << (detail.empty() ? "" : " (" + detail + ")") << "\n";
   }
}

long long Sum(long long x, MPI_Comm comm)
{
   MPI_Allreduce(MPI_IN_PLACE, &x, 1, MPI_LONG_LONG, MPI_SUM, comm);
   return x;
}

bool All(bool ok, MPI_Comm comm)
{
   int v = ok ? 1 : 0;
   MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_INT, MPI_MIN, comm);
   return v == 1;
}

bool Exterior(const ParMesh &m, int be)
{
   return const_cast<ParMesh &>(m).GetFaceInformation(m.GetBdrElementFaceIndex(be)).IsBoundary();
}

std::vector<int> BdrAttrs(const ParMesh &m)
{
   std::vector<int> a(m.GetNBE());
   for (int be = 0; be < m.GetNBE(); be++) { a[be] = m.GetBdrAttribute(be); }
   return a;
}

std::vector<int> ElemAttrs(const ParMesh &m)
{
   std::vector<int> a(m.GetNE());
   for (int e = 0; e < m.GetNE(); e++) { a[e] = m.GetAttribute(e); }
   return a;
}

/// The vertices' largest |x_c - v| over boundary element be (0: the face lies on x_c = v).
real_t OffPlane(const ParMesh &m, int be, int c, real_t v)
{
   Array<int> vs;
   m.GetBdrElementVertices(be, vs);
   real_t d = 0.0;
   for (int k = 0; k < vs.Size(); k++) { d = std::max(d, std::abs(m.GetVertex(vs[k])[c] - v)); }
   return d;
}

std::string Counts(const cycle::OuterFaceTags &t)
{
   return std::to_string(t.natural) + " absorbing, " + std::to_string(t.dirichlet) + " walls, "
          + std::to_string(t.free_surface) + " free surface";
}

void TestConfig()
{
   if (g_rank == 0) { std::cout << "\nconfig\n"; }
   cycle::AbsorbingConfig a;
   const bool off = !a.On() && !a.Walls() && a.FreeSurface();
   a.faces = "bottom";
   const bool bottom = a.On() && !a.Walls() && a.FreeSurface();
   a.faces = "outer";
   const bool outer = a.On() && a.Walls() && a.FreeSurface();
   a.faces = "all";
   const bool all = a.On() && a.Walls() && !a.FreeSurface();
   Check(off && bottom && outer, "config: \"off\", \"bottom\", \"outer\" as before");
   Check(all, "config: \"all\" is on, has walls, keeps no free surface");
}

/// The TPV205-2D box: every exterior face of the outer boundary moves, as a
/// wall (Dirichlet in QD) or as an absorbing traction-free face.
void TestMapView(MPI_Comm comm, const std::string &path)
{
   if (g_rank == 0) { std::cout << "\nmap_view: " << path << "\n"; }
   Mesh smesh(path.c_str(), 1, 1);
   MFEM_VERIFY(smesh.Dimension() == 2, path << " is not a 2-D mesh");
   ParMesh pm(comm, smesh);
   const int fault = 103, outer = 105;
   long long ext = 0, fault_faces = 0;
   for (int be = 0; be < pm.GetNBE(); be++)
   {
      if (pm.GetBdrAttribute(be) == outer && Exterior(pm, be)) { ext++; }
      if (pm.GetBdrAttribute(be) == fault) { fault_faces++; }
   }
   ext = Sum(ext, comm);
   fault_faces = Sum(fault_faces, comm);
   const std::vector<int> a0 = BdrAttrs(pm), e0 = ElemAttrs(pm);
   const std::array<real_t, 3> ez = {0.0, 0.0, 1.0}, none = {0.0, 0.0, 0.0};

   ParMesh w(pm, true), w0(pm, true), n(pm, true);
   const cycle::OuterFaceTags tw = cycle::RetagOuterFaces(w, {}, {outer}, true, ez, 0.5, false);
   const cycle::OuterFaceTags tw0 = cycle::RetagOuterFaces(w0, {}, {outer}, true, none, 0.5, false);
   const cycle::OuterFaceTags tn = cycle::RetagOuterFaces(n, {outer}, {}, true, ez, 0.5, false);
   if (g_rank == 0)
   {
      std::cout << "  " << ext << " exterior faces (attribute " << outer << "), " << fault_faces
                << " fault faces (" << fault << ")\n  walls:         " << Counts(tw)
                << "\n  traction free: " << Counts(tn) << "\n";
   }
   Check(ext > 0 && fault_faces > 0, "map_view: the box has exterior and fault faces");
   Check(tw.dirichlet == ext && tw.natural == 0 && tw.free_surface == 0,
         "map_view: loading walls, every exterior face is a wall", Counts(tw));
   Check(tn.natural == ext && tn.dirichlet == 0 && tn.free_surface == 0,
         "map_view: traction-free box, every exterior face absorbs", Counts(tn));
   bool wall_ok = true, nat_ok = true, same = true, kept = true;
   for (int be = 0; be < pm.GetNBE(); be++)
   {
      const bool moved = a0[be] == outer && Exterior(pm, be);
      wall_ok = wall_ok && w.GetBdrAttribute(be) == (moved ? tw.dirichlet_to : a0[be]);
      nat_ok = nat_ok && n.GetBdrAttribute(be) == (moved ? tn.natural_to : a0[be]);
      same = same && w0.GetBdrAttribute(be) == w.GetBdrAttribute(be);
      kept = kept && (a0[be] != fault || (w.GetBdrAttribute(be) == fault
                                          && n.GetBdrAttribute(be) == fault));
   }
   Check(All(wall_ok && nat_ok, comm), "map_view: each face's new attribute, face by face");
   Check(All(same, comm) && tw0.dirichlet == tw.dirichlet,
         "map_view: up is not read (up = e_z and up = 0 give the same attributes)");
   Check(All(kept, comm), "map_view: the fault faces keep attribute 103");
   const bool elems = ElemAttrs(w) == e0 && ElemAttrs(n) == e0;
   const bool nbe = w.GetNBE() == pm.GetNBE() && n.GetNBE() == pm.GetNBE();
   Check(All(elems && nbe, comm), "map_view: element attributes and boundary element count kept");
   Check(All(w.bdr_attributes.Find(tw.dirichlet_to) >= 0 && w.bdr_attributes.Find(outer) < 0
             && n.bdr_attributes.Find(tn.natural_to) >= 0, comm),
         "map_view: the copies' attribute lists hold the new attributes, not 105");
}

/// BP3: "all" moves what "outer" moves, to the same attributes, and the surface.
void TestBP3(MPI_Comm comm, const std::string &path)
{
   if (g_rank == 0) { std::cout << "\nbp3: " << path << "\n"; }
   Mesh smesh(path.c_str(), 1, 1);
   MFEM_VERIFY(smesh.Dimension() == 2, path << " is not a 2-D mesh");
   ParMesh pm(comm, smesh);
   const std::array<real_t, 3> up = {0.0, 1.0, 0.0};
   ParMesh o(pm, true), a(pm, true);
   const cycle::OuterFaceTags to = cycle::RetagOuterFaces(o, {1}, {5}, true, up, 0.5, true);
   const cycle::OuterFaceTags ta = cycle::RetagOuterFaces(a, {1}, {5}, true, up, 0.5, false);
   if (g_rank == 0)
   {
      std::cout << "  \"outer\": " << Counts(to) << "\n  \"all\":   " << Counts(ta) << "\n";
   }
   Check(to.natural > 0 && to.dirichlet > 0 && to.free_surface > 0,
         "bp3: \"outer\" finds bottom, side and surface faces", Counts(to));
   Check(ta.natural == to.natural + to.free_surface && ta.dirichlet == to.dirichlet
         && ta.free_surface == 0, "bp3: \"all\" = \"outer\" plus the surface", Counts(ta));
   Check(ta.natural_to == to.natural_to && ta.dirichlet_to == to.dirichlet_to,
         "bp3: the same new attributes");
   bool surface = true, creep = true;
   long long surf = 0;
   for (int be = 0; be < pm.GetNBE(); be++)
   {
      const int p = pm.GetBdrAttribute(be), ao = o.GetBdrAttribute(be), aa = a.GetBdrAttribute(be);
      if (ao == aa) { continue; }
      // the only faces "all" treats otherwise: "outer"'s free surface, at y = 0
      const bool s = p == 1 && ao == 1 && aa == ta.natural_to && Exterior(pm, be)
                     && OffPlane(pm, be, 1, 0.0) == 0.0;
      surface = surface && s;
      surf += s ? 1 : 0;
   }
   for (int be = 0; be < pm.GetNBE(); be++)
   {
      if (pm.GetBdrAttribute(be) == 5 && !Exterior(pm, be))
      {
         creep = creep && o.GetBdrAttribute(be) == 5 && a.GetBdrAttribute(be) == 5;
      }
   }
   surf = Sum(surf, comm);
   Check(All(surface, comm) && surf == to.free_surface,
         "bp3: faces that differ are exactly the surface faces (y = 0)",
         std::to_string(surf) + " faces");
   Check(All(creep, comm), "bp3: the creep segment (interior) keeps attribute 5");
}

/// 3-D: a buried box; "outer" keeps the top, "all" absorbs it.
void Test3D(MPI_Comm comm)
{
   if (g_rank == 0) { std::cout << "\nbox_3d\n"; }
   Mesh smesh = Mesh::MakeCartesian3D(2, 2, 2, Element::TETRAHEDRON, 1.0, 1.0, 1.0);
   ParMesh pm(comm, smesh);
   const std::array<real_t, 3> up = {0.0, 0.0, 1.0};
   ParMesh o(pm, true), a(pm, true);
   const std::set<int> walls = {1, 2, 3, 4, 5};
   const cycle::OuterFaceTags to = cycle::RetagOuterFaces(o, {6}, walls, true, up, 0.5, true);
   const cycle::OuterFaceTags ta = cycle::RetagOuterFaces(a, {6}, walls, true, up, 0.5, false);
   long long top = 0;
   bool on_top = true;
   for (int be = 0; be < pm.GetNBE(); be++)
   {
      if (pm.GetBdrAttribute(be) != 6) { continue; }
      top++;
      on_top = on_top && OffPlane(pm, be, 2, 1.0) == 0.0;
   }
   top = Sum(top, comm);
   if (g_rank == 0)
   {
      std::cout << "  " << top << " top faces\n  \"outer\": " << Counts(to) << "\n  \"all\":   "
                << Counts(ta) << "\n";
   }
   Check(All(on_top, comm) && top == 8, "box_3d: attribute 6 is the top z = 1, 8 triangles");
   Check(to.free_surface == top && to.natural == 0 && to.dirichlet == 40,
         "box_3d: \"outer\" keeps the top, 40 walls", Counts(to));
   Check(ta.natural == top && ta.free_surface == 0 && ta.dirichlet == 40,
         "box_3d: \"all\" absorbs the top too", Counts(ta));
}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   g_rank = Mpi::WorldRank();
   const char *deck = "";
   const char *map_mesh = "tpv205/mesh/tpv205_2d_100m.msh";
   const char *bp3_mesh = "bp3/mesh/bp3_d60_1km.msh";
   OptionsParser args(argc, argv);
   args.AddOption(&deck, "-ld", "--load-deck", "Only load this deck and print [cycle.absorbing].");
   args.AddOption(&map_mesh, "-mm", "--map-mesh", "A 2-D box with fault 103 and outer boundary 105.");
   args.AddOption(&bp3_mesh, "-bm", "--bp3-mesh", "A BP3 mesh (1: traction free, 5: Dirichlet).");
   args.Parse();
   if (!args.Good())
   {
      if (g_rank == 0) { args.PrintUsage(std::cout); }
      return 1;
   }
   if (std::string(deck).size() > 0)
   {
      const cycle::CycleConfig cfg = cycle::CycleConfig::Load(deck);
      if (g_rank == 0)
      {
         std::cout << "[cycle.absorbing] faces = \"" << cfg.absorbing.faces << "\": on "
                   << cfg.absorbing.On() << ", walls " << cfg.absorbing.Walls()
                   << ", free surface " << cfg.absorbing.FreeSurface() << "\n";
      }
      return 0;
   }
   if (g_rank == 0)
   {
      std::cout << "test_absorbing_all: ranks " << Mpi::WorldSize() << "\n";
   }
   TestConfig();
   TestMapView(MPI_COMM_WORLD, map_mesh);
   TestBP3(MPI_COMM_WORLD, bp3_mesh);
   Test3D(MPI_COMM_WORLD);
   if (g_rank == 0)
   {
      std::cout << "\ntest_absorbing_all: " << (g_failures ? "FAILED" : "PASSED") << " ("
                << g_failures << " failing checks)\n";
   }
   return g_failures ? 1 : 0;
}

// Nitsche dynamic rupture, the parallel mesh from parts sent by rank 0
// (dynamic/mesh_parts.hpp; Chunhui Zhao 2026-10-22, partition-memory report
// rev1 option D, speed-next S1).  On each rank the ParMesh that
// DistributeMeshParts builds is compared, entity by entity, with
// ParMesh(comm, serial mesh, partition) for the same serial mesh and
// partition:
//
//   MP-1  local vertices (coordinates bit for bit, in order), elements and
//         boundary elements (attributes and vertex lists, in order), the
//         attribute and boundary attribute lists;
//   MP-2  edges and faces (vertex lists, face elements and infos; a shared
//         face's second element is left out, its numbering follows the
//         shared-face order);
//   MP-3  shared vertices, edges and faces as sets, with the ranks that
//         share them; the shared-face lists are allowed to be permutations;
//   MP-4  face-neighbour elements (as a set: attribute and vertex
//         coordinates) after ExchangeFaceNbrData;
//   MP-5  the global true DOF number of every local DOF of an H1 (p = 2,
//         three components) and an L2 (p = 2) space;
//   MP-6  the element values sent with the parts are those of the local
//         elements of ParMesh(comm, serial mesh, partition);
//   MP-7  an interior-penalty DG operator applied to the same vector agrees
//         to round-off (the shared-face order changes the order of the sum);
//   MP-8  the part files (WriteMeshParts on rank 0, ReadMeshPart on every
//         rank) give the ParMesh and the element values DistributeMeshParts
//         gives, bit for bit (ParPrint), and mesh_parts.txt reads back.
//
// Two partitions: METIS (Mesh::GeneratePartitioning, what the driver uses
// without LTS) and stripes that cut every fault face (the fault is an
// interior boundary, attribute 7, at x = 1 of the built-in mesh).  Run on
// 1-4 ranks (make test-mesh-parts, which then removes the part files in
// test_mesh_parts_out/).  --mesh m.msh tests a Gmsh mesh instead
// (METIS partition only; --fault-attr a counts its fault faces, default 7).

#include "mfem.hpp"
#include "../../dynamic/mesh_parts.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <sys/stat.h>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;
static int g_rank = 0, g_size = 1;
static int g_fault_attr = 7;

// every rank checks its own mesh; the test passes when it passes on all ranks
#define TEST_ASSERT(cond, msg) do { \
   int ok_ = (cond) ? 1 : 0, all_ = 0; \
   MPI_Allreduce(&ok_, &all_, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD); \
   num_tests++; \
   if (all_) { num_passed++; if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

/// [0,2]x[0,1]x[0,1] in n x n x n cubes of 6 tets, with the interior faces
/// at x = 1 as boundary elements of attribute 7 (a fault), finalized as
/// Mesh(path, 1, 1) does (edges, tets marked for refinement, orientation).
static Mesh FaultBox(int n)
{
   Mesh box = Mesh::MakeCartesian3D(2 * n, n, n, Element::TETRAHEDRON, 2.0, 1.0, 1.0);
   Mesh m(3, box.GetNV(), box.GetNE(), box.GetNBE() + 2 * n * n);
   for (int v = 0; v < box.GetNV(); v++) { m.AddVertex(box.GetVertex(v)); }
   Array<int> vs;
   for (int e = 0; e < box.GetNE(); e++)
   {
      box.GetElementVertices(e, vs);
      m.AddTet(vs.GetData(), box.GetAttribute(e));
   }
   for (int b = 0; b < box.GetNBE(); b++)
   {
      box.GetBdrElementVertices(b, vs);
      m.AddBdrTriangle(vs.GetData(), box.GetBdrAttribute(b));
   }
   int nfault = 0;
   for (int f = 0; f < box.GetNFaces(); f++)
   {
      int e1, e2;
      box.GetFaceElements(f, &e1, &e2);
      if (e2 < 0) { continue; }
      box.GetFaceVertices(f, vs);
      bool on = true;
      for (int j = 0; j < vs.Size(); j++) { on = on && box.GetVertex(vs[j])[0] == 1.0; }
      if (on) { m.AddBdrTriangle(vs.GetData(), 7); nfault++; }
   }
   MFEM_VERIFY(nfault == 2 * n * n, "fault faces: " << nfault);
   m.FinalizeTetMesh(1, 1, true);
   return m;
}

static std::array<real_t, 3> Centroid(const Mesh &m, int e)
{
   Array<int> v;
   m.GetElementVertices(e, v);
   std::array<real_t, 3> c{0.0, 0.0, 0.0};
   for (int k = 0; k < v.Size(); k++)
      for (int d = 0; d < 3; d++) { c[d] += m.GetVertex(v[k])[d] / v.Size(); }
   return c;
}

static bool SameVertexLists(const Array<int> &a, const Array<int> &b)
{
   if (a.Size() != b.Size()) { return false; }
   for (int i = 0; i < a.Size(); i++) { if (a[i] != b[i]) { return false; } }
   return true;
}

/// The members of the group of shared vertex / edge / face entity, as world ranks.
static std::vector<int> GroupRanks(const ParMesh &m, int group)
{
   const GroupTopology &gt = m.gtopo;
   std::vector<int> r;
   for (int j = 0; j < gt.GetGroupSize(group); j++)
   {
      r.push_back(gt.GetNeighborRank(gt.GetGroup(group)[j]));
   }
   std::sort(r.begin(), r.end());
   return r;
}

/// Shared entities as a set of (sorted local vertices, group ranks).
static std::set<std::pair<std::vector<int>, std::vector<int>>>
SharedSet(const ParMesh &m, int kind)
{
   std::set<std::pair<std::vector<int>, std::vector<int>>> s;
   for (int g = 1; g < m.GetNGroups(); g++)
   {
      const int n = kind == 0 ? m.GroupNVertices(g) : kind == 1 ? m.GroupNEdges(g)
                    : m.GroupNTriangles(g);
      for (int i = 0; i < n; i++)
      {
         std::vector<int> v;
         if (kind == 0) { v.push_back(m.GroupVertex(g, i)); }
         else if (kind == 1)
         {
            int le, o;
            m.GroupEdge(g, i, le, o);
            Array<int> ev;
            m.GetEdgeVertices(le, ev);
            v.assign(ev.begin(), ev.end());
         }
         else
         {
            int lf, o;
            m.GroupTriangle(g, i, lf, o);
            Array<int> fv;
            m.GetFaceVertices(lf, fv);
            v.assign(fv.begin(), fv.end());
         }
         std::sort(v.begin(), v.end());
         s.insert({v, GroupRanks(m, g)});
      }
   }
   return s;
}

static const std::string kOut = "test_mesh_parts_out";

static std::string ParText(const ParMesh &m)
{
   std::ostringstream os;
   os.precision(17);
   m.ParPrint(os);
   return os.str();
}

static void Compare(Mesh &smesh, const std::vector<int> &part, const char *name,
                    const std::string &tag)
{
   if (g_rank == 0)
   {
      std::cout << "--- " << name << ": " << smesh.GetNE() << " tets, " << smesh.GetNBE()
                << " boundary elements, " << g_size << " rank(s)\n";
   }
   std::vector<int> values(smesh.GetNE());
   for (int e = 0; e < smesh.GetNE(); e++) { values[e] = (37 * e + 11) % 5; }

   ParMesh a(MPI_COMM_WORLD, smesh, part.data());
   std::vector<int> va;   // the values of a's local elements (serial order)
   for (int e = 0; e < smesh.GetNE(); e++) { if (part[e] == g_rank) { va.push_back(values[e]); } }

   std::vector<int> vb;
   std::unique_ptr<ParMesh> bp = DistributeMeshParts(
      MPI_COMM_WORLD, g_rank == 0 ? &smesh : nullptr, g_rank == 0 ? part.data() : nullptr,
      true, g_rank == 0 ? values.data() : nullptr, &vb);
   ParMesh &b = *bp;

   // MP-1
   bool vert = a.GetNV() == b.GetNV();
   for (int v = 0; vert && v < a.GetNV(); v++)
      for (int d = 0; d < 3; d++) { vert = vert && a.GetVertex(v)[d] == b.GetVertex(v)[d]; }
   TEST_ASSERT(vert, "MP-1 local vertices: same count and coordinates, bit for bit");
   bool el = a.GetNE() == b.GetNE();
   Array<int> x, y;
   for (int e = 0; el && e < a.GetNE(); e++)
   {
      a.GetElementVertices(e, x);
      b.GetElementVertices(e, y);
      el = a.GetAttribute(e) == b.GetAttribute(e) && SameVertexLists(x, y);
   }
   TEST_ASSERT(el, "MP-1 elements: same attributes and vertex lists, in order");
   bool be = a.GetNBE() == b.GetNBE();
   for (int i = 0; be && i < a.GetNBE(); i++)
   {
      a.GetBdrElementVertices(i, x);
      b.GetBdrElementVertices(i, y);
      be = a.GetBdrAttribute(i) == b.GetBdrAttribute(i) && SameVertexLists(x, y);
   }
   TEST_ASSERT(be, "MP-1 boundary elements: same attributes and vertex lists, in order");
   int nfault = 0;
   for (int i = 0; i < a.GetNBE(); i++) { nfault += a.GetBdrAttribute(i) == g_fault_attr; }
   TEST_ASSERT(SameVertexLists(a.attributes, b.attributes)
               && SameVertexLists(a.bdr_attributes, b.bdr_attributes),
               "MP-1 attribute and boundary attribute lists");

   // MP-2
   bool fc = a.GetNEdges() == b.GetNEdges() && a.GetNFaces() == b.GetNFaces();
   for (int i = 0; fc && i < a.GetNEdges(); i++)
   {
      a.GetEdgeVertices(i, x);
      b.GetEdgeVertices(i, y);
      fc = SameVertexLists(x, y);
   }
   std::set<int> sfa;
   for (int i = 0; i < a.GetNSharedFaces(); i++) { sfa.insert(a.GetSharedFace(i)); }
   for (int f = 0; fc && f < a.GetNFaces(); f++)
   {
      a.GetFaceVertices(f, x);
      b.GetFaceVertices(f, y);
      int a1, a2, b1, b2, ai1, ai2, bi1, bi2;
      a.GetFaceElements(f, &a1, &a2);
      b.GetFaceElements(f, &b1, &b2);
      a.GetFaceInfos(f, &ai1, &ai2);
      b.GetFaceInfos(f, &bi1, &bi2);
      const bool shared = sfa.count(f) > 0;
      fc = SameVertexLists(x, y) && a1 == b1 && ai1 == bi1 && (shared || (a2 == b2 && ai2 == bi2));
   }
   TEST_ASSERT(fc, "MP-2 edges and faces: same vertex lists, face elements and infos");

   // MP-3
   std::set<int> sfb;
   int moved = 0;
   for (int i = 0; i < b.GetNSharedFaces(); i++)
   {
      sfb.insert(b.GetSharedFace(i));
      moved += i < a.GetNSharedFaces() && a.GetSharedFace(i) != b.GetSharedFace(i);
   }
   TEST_ASSERT(a.GetNSharedFaces() == b.GetNSharedFaces() && sfa == sfb,
               "MP-3 shared faces: the same set of local faces");
   TEST_ASSERT(SharedSet(a, 0) == SharedSet(b, 0) && SharedSet(a, 1) == SharedSet(b, 1)
               && SharedSet(a, 2) == SharedSet(b, 2),
               "MP-3 shared vertices, edges and faces with the same ranks");
   long long stat[3] = {a.GetNSharedFaces(), moved, nfault}, tot[3];
   MPI_Reduce(stat, tot, 3, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
   if (g_rank == 0)
   {
      std::cout << "  (shared faces " << tot[0] << ", of which " << tot[1]
                << " at another list position; fault boundary elements " << tot[2] << ")\n";
   }

   // MP-4
   a.ExchangeFaceNbrData();
   b.ExchangeFaceNbrData();
   auto nbr_set = [](ParMesh &m) {
      std::multiset<std::vector<real_t>> s;
      for (int j = 0; j < m.GetNFaceNeighborElements(); j++)
      {
         const Element *e = m.face_nbr_elements[j];
         std::vector<real_t> k{static_cast<real_t>(e->GetAttribute())};
         for (int i = 0; i < e->GetNVertices(); i++)
            for (int d = 0; d < 3; d++) { k.push_back(m.face_nbr_vertices[e->GetVertices()[i]](d)); }
         s.insert(k);
      }
      return s;
   };
   TEST_ASSERT(a.GetNFaceNeighbors() == b.GetNFaceNeighbors()
               && a.GetNFaceNeighborElements() == b.GetNFaceNeighborElements()
               && nbr_set(a) == nbr_set(b),
               "MP-4 face-neighbour elements: same set (attribute, vertex coordinates)");

   // MP-5
   H1_FECollection h1(2, 3);
   L2_FECollection l2(2, 3);
   ParFiniteElementSpace ha(&a, &h1, 3), hb(&b, &h1, 3), la(&a, &l2), lb(&b, &l2);
   bool dof = ha.GlobalTrueVSize() == hb.GlobalTrueVSize() && ha.GetVSize() == hb.GetVSize()
              && la.GlobalTrueVSize() == lb.GlobalTrueVSize() && la.GetVSize() == lb.GetVSize();
   for (int i = 0; dof && i < ha.GetVSize(); i++)
   {
      dof = ha.GetGlobalTDofNumber(i) == hb.GetGlobalTDofNumber(i);
   }
   for (int i = 0; dof && i < la.GetVSize(); i++)
   {
      dof = la.GetGlobalTDofNumber(i) == lb.GetGlobalTDofNumber(i);
   }
   TEST_ASSERT(dof, "MP-5 H1 (p = 2, 3 components) and L2 (p = 2) global true DOF numbers");

   // MP-6
   TEST_ASSERT(va == vb, "MP-6 the element values sent with each part");

   // MP-7
   L2_FECollection dg(1, 3);
   ParFiniteElementSpace da(&a, &dg), db(&b, &dg);
   ConstantCoefficient one(1.0);
   auto apply = [&](ParFiniteElementSpace &fes) {
      ParBilinearForm k(&fes);
      k.AddDomainIntegrator(new DiffusionIntegrator(one));
      k.AddInteriorFaceIntegrator(new DGDiffusionIntegrator(one, -1.0, 10.0));
      k.Assemble();
      k.Finalize();
      std::unique_ptr<HypreParMatrix> K(k.ParallelAssemble());
      ParGridFunction u(&fes);
      FunctionCoefficient f([](const Vector &p) { return std::sin(3.0 * p(0)) + p(1) * p(2); });
      u.ProjectCoefficient(f);
      std::unique_ptr<HypreParVector> U(u.ParallelProject());
      HypreParVector R(*U);
      K->Mult(*U, R);
      return Vector(R);
   };
   Vector ra = apply(da), rb = apply(db);
   real_t dmax = 0.0, rmax = 0.0;
   for (int i = 0; i < ra.Size(); i++)
   {
      dmax = std::max(dmax, std::abs(ra(i) - rb(i)));
      rmax = std::max(rmax, std::abs(ra(i)));
   }
   real_t g[2] = {dmax, rmax}, gm[2];
   MPI_Allreduce(g, gm, 2, MPITypeMap<real_t>::mpi_type, MPI_MAX, MPI_COMM_WORLD);
   if (g_rank == 0)
   {
      std::cout << "  (DG operator: max |difference| " << gm[0] << " against max |entry| "
                << gm[1] << ")\n";
   }
   TEST_ASSERT(ra.Size() == rb.Size() && gm[0] <= 1e-12 * gm[1],
               "MP-7 interior-penalty DG operator: same result to round-off");

   // MP-8
   const std::string dir = kOut + "/" + tag;
   if (g_rank == 0)
   {
      ::mkdir(kOut.c_str(), 0755);
      ::mkdir(dir.c_str(), 0755);
      WriteMeshParts(dir, smesh, part.data(), g_size, true, values.data());
      WriteMeshPartsInfo(dir, {{"format", "test 1"}, {"name", name}});
   }
   MPI_Barrier(MPI_COMM_WORLD);
   std::vector<int> vc;
   std::unique_ptr<ParMesh> cp = ReadMeshPart(MPI_COMM_WORLD, dir, true, &vc);
   std::unique_ptr<ParMesh> bq = DistributeMeshParts(
      MPI_COMM_WORLD, g_rank == 0 ? &smesh : nullptr, g_rank == 0 ? part.data() : nullptr,
      true, g_rank == 0 ? values.data() : nullptr, &vb);
   TEST_ASSERT(ParText(*cp) == ParText(*bq) && vc == vb,
               "MP-8 part files: the same ParMesh (ParPrint) and element values as the sent parts");
   const std::map<std::string, std::string> kv = ReadMeshPartsInfo(MPI_COMM_WORLD, dir);
   TEST_ASSERT(kv.size() == 2 && kv.at("format") == "test 1" && kv.at("name") == name,
               "MP-8 mesh_parts.txt: the entries read back");
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   g_size = Mpi::WorldSize();
   std::string mesh_file;
   for (int i = 1; i + 1 < argc; i++)
   {
      if (std::string(argv[i]) == "--mesh") { mesh_file = argv[i + 1]; }
      if (std::string(argv[i]) == "--fault-attr") { g_fault_attr = std::atoi(argv[i + 1]); }
   }
   if (g_rank == 0) { std::cout << "test_mesh_parts on " << g_size << " rank(s)\n"; }

   Mesh smesh = mesh_file.empty() ? FaultBox(4) : Mesh(mesh_file.c_str(), 1, 1);
   const int ne = smesh.GetNE();
   {
      int *p = smesh.GeneratePartitioning(g_size, 1);
      std::vector<int> part(p, p + ne);
      delete [] p;
      Compare(smesh, part, "METIS partition", "metis");
   }
   if (mesh_file.empty())
   {
      // stripes in y, the two sides of the fault on different ranks
      std::vector<int> part(ne);
      for (int e = 0; e < ne; e++)
      {
         const std::array<real_t, 3> c = Centroid(smesh, e);
         part[e] = (2 * static_cast<int>(c[1] * g_size) + (c[0] < 1.0 ? 0 : 1)) % g_size;
      }
      Compare(smesh, part, "stripes across the fault", "stripes");
   }

   if (g_rank == 0)
   {
      std::cout << "\n" << num_passed << " of " << num_tests << " passed, " << num_failed
                << " failed\n";
   }
   return num_failed ? 1 : 0;
}

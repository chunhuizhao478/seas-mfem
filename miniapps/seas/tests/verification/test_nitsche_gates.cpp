// Nitsche dynamic rupture, phase P0 (plan gates G0.2 and G0.5).
//
// Read-only probes of the existing quasi-static ElasticityDomainOperator on
// the TPV102 inline box (60 x 60 x 30 km, fault at Y=0, free surface at Z=0):
//   G0.2-1: every face tagged 3 is classified as a fault face
//   G0.2-2: fault DOF count / area match the tagged fault area
//   G0.2-3: GetFaultDOFBasis frames are orthonormal; frames on the Z=0
//           edge equal the mid-plane frame up to the per-face sign
//   G0.2-4: ComputeTraction on a manufactured displacement matches the
//           analytic traction at surface-adjacent fault faces
//   G0.5:   count fault faces whose stored normal is anti-aligned with
//           ref_normal (0,-1,0); on each, check t1 x t2 = -n (left-handed)
//   Deck:   tpv102_spatial.toml with a constant_tensor [stress] block
//           resolves to the TPV102Params material (cs, cp, rho).
//
// The checks run twice: on the hex box with DGMethod::BR2 (one DOF per
// face), and on the same box split into 24 tets per hex with DGMethod::IP at
// p = 1, 2, 3 (the IP fault path is TRIANGLE-only; nbf = (p+1)(p+2)/2), plus
// a Y -> -Y mirror-symmetry check of the tet mesh.
//
// No production code is touched.  FacetBC/face_bc_ are private, so G0.2-1
// is checked through the public GetFaultInteriorFaces() (the operator itself
// MFEM_VERIFYs that every face in that list is FacetBC::Fault).

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../spatial/code/spatial_friction.hpp"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <filesystem>
#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <vector>

using namespace mfem;
using namespace mfem::seas;
namespace sp = mfem::seas::spatial;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

#define TEST_NEAR(val, exp, tol, msg) do { \
   num_tests++; \
   double v_ = (val), e_ = (exp), t_ = (tol); \
   if (std::abs(v_ - e_) < t_) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg \
      << " (got " << v_ << ", expected " << e_ << ", diff " << std::abs(v_-e_) << ")\n"; } \
} while(0)

// ---------------------------------------------------------------------------
// Helper: Create a small inline hex mesh for testing
// (60km x 60km x 30km, with fault at Y=0, free surface at Z=0)
// Copied verbatim from tests/verification/test_tpv102_local.cpp.
// ---------------------------------------------------------------------------
static Mesh *CreateTestMesh(int nx, int ny, int nz,
                            real_t Lx, real_t Ly, real_t Lz)
{
   // Create [-Lx, Lx] x [-Ly, Ly] x [-Lz, 0] hex mesh.
   // Note: Y=0 is an interior plane in a Cartesian mesh, so NO fault boundary
   // faces are possible. The inline mesh tests wave propagation and BCs only.
   // Fault coupling requires a split mesh (from Gmsh .geo files).
   Mesh *mesh = new Mesh(Mesh::MakeCartesian3D(
      nx, ny, nz, Element::HEXAHEDRON, 2*Lx, 2*Ly, Lz));

   // Shift vertices directly (not nodes, so GetVertex() returns correct coords)
   for (int v = 0; v < mesh->GetNV(); v++)
   {
      real_t *coord = mesh->GetVertex(v);
      coord[0] -= Lx;   // X: [-Lx, Lx]
      coord[1] -= Ly;   // Y: [-Ly, Ly]
      coord[2] -= Lz;   // Z: [-Lz, 0]
   }

   // Set boundary attributes:
   //   Z = 0 (top) -> attr 1 (free surface)
   //   Everything else -> attr 5 (absorbing)
   // No fault faces (Y=0 is interior in Cartesian mesh)
   for (int b = 0; b < mesh->GetNBE(); b++)
   {
      Vector center(3);
      center = 0.0;
      Array<int> verts;
      mesh->GetBdrElementVertices(b, verts);
      for (int v = 0; v < verts.Size(); v++)
      {
         const real_t *coord = mesh->GetVertex(verts[v]);
         for (int d = 0; d < 3; d++) { center(d) += coord[d]; }
      }
      center /= verts.Size();

      real_t tol = Lz / nz * 0.1;

      if (std::abs(center(2)) < tol)
      {
         // Z = 0: free surface
         mesh->SetBdrAttribute(b, 1);
      }
      else
      {
         // All other boundaries: absorbing
         mesh->SetBdrAttribute(b, 5);
      }
   }

   return mesh;
}

// ---------------------------------------------------------------------------
// P0 addition (approved deviation "A"): the verbatim builder has no fault
// faces, so tag every interior face on the Y=0 plane with attr 3 via an
// interior boundary quad (same pattern as test_bimaterial_fault_perside_
// material.cpp:66).  Returns the number of faces tagged.
// ---------------------------------------------------------------------------
static int TagFaultPlaneY0(Mesh &mesh, real_t tol)
{
   std::vector<std::array<int, 4>> quads;
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      if (mesh.GetFaceInformation(f).IsBoundary()) { continue; }
      Array<int> fv;
      mesh.GetFaceVertices(f, fv);
      bool on_plane = true;
      for (int i = 0; i < fv.Size(); i++)
      {
         if (std::abs(mesh.GetVertex(fv[i])[1]) > tol) { on_plane = false; }
      }
      if (on_plane) { quads.push_back({fv[0], fv[1], fv[2], fv[3]}); }
   }
   for (const auto &q : quads) { mesh.AddBdrQuad(q[0], q[1], q[2], q[3], 3); }
   mesh.FinalizeTopology();
   mesh.Finalize();
   mesh.SetAttributes();
   return static_cast<int>(quads.size());
}

// Manufactured displacement (continuous across Y=0, zero slip), quadratic so
// the fault traction is linear in (x, z) and its face average equals the value
// at the face centroid:
//   u_x = A y z,  u_y = B x z,  u_z = C z^2 + D x y
static constexpr real_t kA = 1.0e-9, kB = 2.0e-9, kC = 0.5e-9, kD = -1.0e-9;

static void ManufacturedU(const Vector &x, Vector &u)
{
   u.SetSize(3);
   u(0) = kA * x(1) * x(2);
   u(1) = kB * x(0) * x(2);
   u(2) = kC * x(2) * x(2) + kD * x(0) * x(1);
}

static void ManufacturedStress(const real_t x[3], real_t lam, real_t mu,
                               real_t S[3][3])
{
   real_t G[3][3] = {
      {0.0,          kA * x[2],    kA * x[1]},
      {kB * x[2],    0.0,          kB * x[0]},
      {kD * x[1],    kD * x[0],    2.0 * kC * x[2]}};
   const real_t tr = G[0][0] + G[1][1] + G[2][2];
   for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
      {
         S[i][j] = mu * (G[i][j] + G[j][i]) + (i == j ? lam * tr : 0.0);
      }
}

static void Cross(const real_t a[3], const real_t b[3], real_t c[3])
{
   c[0] = a[1] * b[2] - a[2] * b[1];
   c[1] = a[2] * b[0] - a[0] * b[2];
   c[2] = a[0] * b[1] - a[1] * b[0];
}

static real_t Dot(const real_t a[3], const real_t b[3])
{
   return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static void FaceCentroid(Mesh &mesh, int face, real_t c[3])
{
   Array<int> fv;
   mesh.GetFaceVertices(face, fv);
   c[0] = c[1] = c[2] = 0.0;
   for (int i = 0; i < fv.Size(); i++)
      for (int d = 0; d < 3; d++) { c[d] += mesh.GetVertex(fv[i])[d]; }
   for (int d = 0; d < 3; d++) { c[d] /= fv.Size(); }
}

// Linear manufactured displacement (constant gradient) for p = 1, where the
// quadratic field above is not in the DG space.
static constexpr real_t kGL[3][3] = {{0.0,     1.0e-5, 0.0},
                                     {2.0e-5,  0.0,    0.7e-5},
                                     {0.0,    -1.0e-5, 0.5e-5}};

static void ManufacturedULinear(const Vector &x, Vector &u)
{
   u.SetSize(3);
   for (int i = 0; i < 3; i++)
   {
      u(i) = kGL[i][0] * x(0) + kGL[i][1] * x(1) + kGL[i][2] * x(2);
   }
}

static void ManufacturedStressLinear(real_t lam, real_t mu, real_t S[3][3])
{
   const real_t tr = kGL[0][0] + kGL[1][1] + kGL[2][2];
   for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
      {
         S[i][j] = mu * (kGL[i][j] + kGL[j][i]) + (i == j ? lam * tr : 0.0);
      }
}

// ---------------------------------------------------------------------------
// The same box split into 24 tets per hex (each hex face -> 4 triangles
// through the face centre), so the Y=0 plane carries triangular fault faces
// and the IP fault path (FaceQuadrature is TRIANGLE-only) applies.  Vertex
// shift and attrs 1/5 follow CreateTestMesh; fault attr 3 added on Y=0.
// ---------------------------------------------------------------------------
static Mesh *CreateTestMesh24Tet(int nx, int ny, int nz,
                                 real_t Lx, real_t Ly, real_t Lz)
{
   Mesh *mesh = new Mesh(Mesh::MakeCartesian3DWith24TetsPerHex(
      nx, ny, nz, 2*Lx, 2*Ly, Lz));
   for (int v = 0; v < mesh->GetNV(); v++)
   {
      real_t *coord = mesh->GetVertex(v);
      coord[0] -= Lx;
      coord[1] -= Ly;
      coord[2] -= Lz;
   }
   for (int b = 0; b < mesh->GetNBE(); b++)
   {
      Vector center(3);
      center = 0.0;
      Array<int> verts;
      mesh->GetBdrElementVertices(b, verts);
      for (int v = 0; v < verts.Size(); v++)
      {
         const real_t *coord = mesh->GetVertex(verts[v]);
         for (int d = 0; d < 3; d++) { center(d) += coord[d]; }
      }
      center /= verts.Size();
      real_t tol = Lz / nz * 0.1;
      mesh->SetBdrAttribute(b, (std::abs(center(2)) < tol) ? 1 : 5);
   }

   const real_t ytol = 1e-6 * Ly;
   std::vector<std::array<int, 3>> tris;
   for (int f = 0; f < mesh->GetNumFaces(); f++)
   {
      if (mesh->GetFaceInformation(f).IsBoundary()) { continue; }
      Array<int> fv;
      mesh->GetFaceVertices(f, fv);
      bool on_plane = true;
      for (int i = 0; i < fv.Size(); i++)
      {
         if (std::abs(mesh->GetVertex(fv[i])[1]) > ytol) { on_plane = false; }
      }
      if (on_plane) { tris.push_back({fv[0], fv[1], fv[2]}); }
   }
   for (const auto &t : tris) { mesh->AddBdrTriangle(t[0], t[1], t[2], 3); }
   mesh->FinalizeTopology();
   mesh->Finalize();
   mesh->SetAttributes();
   return mesh;
}

// Mirror symmetry Y -> -Y: every element centroid (x,y,z) has a partner at
// (x,-y,z).  Returns the number of centroids without a partner.
static int CountUnmirroredCentroids(Mesh &mesh, real_t h)
{
   auto key = [h](real_t x, real_t y, real_t z) {
      const real_t q = 1e-6 * h;
      return std::array<long long, 3>{std::llround(x / q), std::llround(y / q),
                                      std::llround(z / q)};
   };
   std::set<std::array<long long, 3>> cents;
   std::vector<std::array<real_t, 3>> list;
   for (int e = 0; e < mesh.GetNE(); e++)
   {
      Array<int> ev;
      mesh.GetElementVertices(e, ev);
      real_t c[3] = {0.0, 0.0, 0.0};
      for (int i = 0; i < ev.Size(); i++)
         for (int d = 0; d < 3; d++) { c[d] += mesh.GetVertex(ev[i])[d]; }
      for (int d = 0; d < 3; d++) { c[d] /= ev.Size(); }
      cents.insert(key(c[0], c[1], c[2]));
      list.push_back({c[0], c[1], c[2]});
   }
   int missing = 0;
   for (const auto &c : list)
   {
      if (!cents.count(key(c[0], -c[1], c[2]))) { missing++; }
   }
   return missing;
}

// ---------------------------------------------------------------------------
// IP fault path on the 24-tet box at order p: checks G0.2-1..4 and G0.5.
// ---------------------------------------------------------------------------
static void RunTetIP(int p)
{
   std::cout << "\n##### 24-tet box, DGMethod::IP, p = " << p << " #####\n";
   const int nx = 12, ny = 12, nz = 6;
   const real_t Lx = 30.0e3, Ly = 30.0e3, Lz = 30.0e3;
   const real_t dx = 2.0 * Lx / nx, dz = Lz / nz;

   Mesh *mesh = CreateTestMesh24Tet(nx, ny, nz, Lx, Ly, Lz);
   std::set<int> tagged_faces;
   for (int be = 0; be < mesh->GetNBE(); be++)
   {
      if (mesh->GetBdrAttribute(be) == 3)
      {
         tagged_faces.insert(mesh->GetBdrElementFaceIndex(be));
      }
   }
   std::cout << "  mesh: " << mesh->GetNE() << " tets, " << mesh->GetNBE()
             << " bdr elements, " << tagged_faces.size()
             << " faces tagged attr 3\n";

   const int unmirrored = CountUnmirroredCentroids(*mesh, dx);
   std::cout << "  Y -> -Y mirror: " << unmirrored << " of " << mesh->GetNE()
             << " element centroids without a partner\n";
   TEST_ASSERT(unmirrored == 0, "mesh is mirror-symmetric under Y -> -Y");

   LinearElastic model(TPV102Params::lambda, TPV102Params::mu);
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.dirichlet_attrs = {5};
   ElasticityDomainOperator<Mesh> op(*mesh, p, model,
                                     0.0, Lz, 2.0 * Lx, bc,
                                     DGMethod::IP, SolverType::CG_AMG);

   const Array<int> &ffaces = op.GetFaultInteriorFaces();
   const int nfault = ffaces.Size();
   const int nbf = op.GetNbfPerFace();
   const int ndof = op.GetNumFaultDOFs();

   // ---------------------------------------------------------------- G0.2-1
   std::cout << "\n=== [IP p=" << p << "] G0.2-1: tagged faces classify as fault ===\n";
   {
      std::set<int> fault_set(ffaces.begin(), ffaces.end());
      std::cout << "  tagged = " << tagged_faces.size()
                << ", GetFaultInteriorFaces = " << nfault
                << ", GetFaultSharedFaces = " << op.GetFaultSharedFaces().Size()
                << "\n";
      TEST_ASSERT(static_cast<int>(tagged_faces.size()) == 4 * nx * nz,
                  "tagged face count == 4*nx*nz = " << 4 * nx * nz);
      TEST_ASSERT(fault_set == tagged_faces,
                  "fault face set == tagged attr-3 face set");
      TEST_ASSERT(nfault == static_cast<int>(fault_set.size()),
                  "no duplicate fault faces");
      TEST_ASSERT(op.GetFaultSharedFaces().Size() == 0,
                  "serial mesh: no shared fault faces");
   }

   // ---------------------------------------------------------------- G0.2-2
   std::cout << "\n=== [IP p=" << p << "] G0.2-2: fault DOF count vs tagged area ===\n";
   {
      real_t area = 0.0;
      for (int fi = 0; fi < nfault; fi++)
      {
         auto *T = mesh->GetFaceTransformation(ffaces[fi]);
         const IntegrationRule &ir = IntRules.Get(T->GetGeometryType(), 2);
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            T->SetIntPoint(&ir.IntPoint(q));
            area += ir.IntPoint(q).weight * T->Weight();
         }
      }
      const real_t area_exact = 2.0 * Lx * Lz;
      std::cout << std::setprecision(15)
                << "  nbf_per_face = " << nbf << " (expected "
                << (p + 1) * (p + 2) / 2 << "), num_fault_dofs = " << ndof
                << ", area = " << area << " m^2 (exact " << area_exact << ")\n";
      TEST_ASSERT(nbf == (p + 1) * (p + 2) / 2, "nbf_per_face == (p+1)(p+2)/2");
      TEST_ASSERT(ndof == nfault * nbf, "num_fault_dofs == faces * nbf");
      TEST_NEAR(area, area_exact, 1e-9 * area_exact, "fault area == 2*Lx*Lz");
      TEST_ASSERT(nfault == static_cast<int>(std::lround(4.0 * area / (dx * dz))),
                  "fault faces == 4 * area/(dx*dz)");
   }

   DenseMatrix B;
   op.GetFaultDOFBasis(B);
   Vector xyz;
   op.GetFaultDOFCoords3D(xyz);

   // ---------------------------------------------------------------- G0.2-3
   std::cout << "\n=== [IP p=" << p << "] G0.2-3: per-DOF fault frame ===\n";
   {
      real_t max_orth = 0.0;
      for (int k = 0; k < ndof; k++)
      {
         real_t n[3], t1[3], t2[3];
         for (int d = 0; d < 3; d++)
         {
            n[d] = B(d, k); t1[d] = B(3 + d, k); t2[d] = B(6 + d, k);
         }
         max_orth = std::max({max_orth,
                              std::abs(Dot(n, n) - 1.0),
                              std::abs(Dot(t1, t1) - 1.0),
                              std::abs(Dot(t2, t2) - 1.0),
                              std::abs(Dot(n, t1)), std::abs(Dot(n, t2)),
                              std::abs(Dot(t1, t2))});
      }
      std::cout << "  max orthonormality defect = " << max_orth << "\n";
      TEST_ASSERT(max_orth < 1e-12, "every fault DOF frame is orthonormal");

      // Reference: the fault DOF nodal point closest to (0, 0, -Lz/2).
      int rk = -1;
      real_t best = 1e300;
      for (int k = 0; k < ndof; k++)
      {
         const real_t x = xyz(3 * k), z = xyz(3 * k + 2);
         const real_t d2 = x * x + (z + 0.5 * Lz) * (z + 0.5 * Lz);
         if (d2 < best) { best = d2; rk = k; }
      }
      real_t rn[3], rt1[3], rt2[3];
      for (int d = 0; d < 3; d++)
      {
         rn[d] = B(d, rk); rt1[d] = B(3 + d, rk); rt2[d] = B(6 + d, rk);
      }
      std::cout << "  mid-plane ref DOF at (" << xyz(3 * rk) << ", "
                << xyz(3 * rk + 1) << ", " << xyz(3 * rk + 2) << ")\n"
                << "    n  = (" << rn[0] << ", " << rn[1] << ", " << rn[2] << ")\n"
                << "    t1 = (" << rt1[0] << ", " << rt1[1] << ", " << rt1[2] << ")\n"
                << "    t2 = (" << rt2[0] << ", " << rt2[1] << ", " << rt2[2] << ")\n";

      // DOFs whose nodal point lies on the Z=0 edge itself.
      int n_edge = 0, n_plus = 0, n_minus = 0;
      std::set<int> edge_faces;
      real_t max_dev = 0.0;
      for (int k = 0; k < ndof; k++)
      {
         if (std::abs(xyz(3 * k + 2)) > 1e-9 * Lz) { continue; }
         real_t n[3], t1[3], t2[3];
         for (int d = 0; d < 3; d++)
         {
            n[d] = B(d, k); t1[d] = B(3 + d, k); t2[d] = B(6 + d, k);
         }
         const real_t s = (Dot(n, rn) >= 0.0) ? 1.0 : -1.0;
         (s > 0 ? n_plus : n_minus)++;
         n_edge++;
         edge_faces.insert(k / nbf);
         for (int d = 0; d < 3; d++)
         {
            max_dev = std::max({max_dev, std::abs(n[d] - s * rn[d]),
                                std::abs(t1[d] - s * rt1[d]),
                                std::abs(t2[d] - s * rt2[d])});
         }
      }
      std::cout << "  DOFs with nodal z = 0: " << n_edge << " on "
                << edge_faces.size() << " faces (sign +1: " << n_plus
                << ", sign -1: " << n_minus << "), max |frame - s*ref| = "
                << max_dev << "\n";
      // Each hex face on Y=0 touching Z=0 splits into 4 triangles through
      // its centre: the top one has an edge on Z=0 (p+1 nodal points), the
      // two side ones touch Z=0 at one corner vertex each (1 nodal point).
      TEST_ASSERT(static_cast<int>(edge_faces.size()) == 3 * nx,
                  "faces with a nodal point on Z=0 == 3*nx");
      TEST_ASSERT(n_edge == nx * (p + 3), "Z=0 nodal DOFs == nx*(p+3)");
      TEST_ASSERT(max_dev < 1e-12,
                  "Z=0-edge frames == mid-plane frame up to per-face sign");
   }

   // ---------------------------------------------------------------- G0.2-4
   std::cout << "\n=== [IP p=" << p << "] G0.2-4: ComputeTraction vs analytic at"
             << " fault DOF nodes (" << (p >= 2 ? "quadratic" : "linear")
             << " manufactured u) ===\n";
   {
      GridFunction u(&op.GetFESpace());
      VectorFunctionCoefficient ucoef(3, p >= 2 ? ManufacturedU
                                                : ManufacturedULinear);
      u.ProjectCoefficient(ucoef);

      Vector slip(2 * ndof);
      slip = 0.0;
      Vector trac, sn;
      op.ComputeTraction(u, slip, trac, &sn);

      const FaultBasis *fb = op.GetFaultBasis();
      const real_t lam = TPV102Params::lambda, mu = TPV102Params::mu;
      real_t max_rel_all = 0.0, max_rel_surf = 0.0, scale = 0.0;
      int n_surf = 0, sample = -1;
      for (int k = 0; k < ndof; k++)
      {
         const real_t x[3] = {xyz(3 * k), xyz(3 * k + 1), xyz(3 * k + 2)};
         real_t S[3][3];
         if (p >= 2) { ManufacturedStress(x, lam, mu, S); }
         else { ManufacturedStressLinear(lam, mu, S); }
         // Projection frame = per-face basis (flat faces: same at every QP).
         const FaultBasisData &b = fb->GetBasis(k / nbf);
         real_t T[3];
         for (int i = 0; i < 3; i++)
         {
            T[i] = S[i][0] * b.normal[0] + S[i][1] * b.normal[1]
                   + S[i][2] * b.normal[2];
         }
         const real_t ex_dip = Dot(T, b.tangent1);
         const real_t ex_str = Dot(T, b.tangent2);
         const real_t ex_sn = -Dot(T, b.normal);
         const real_t sc = std::max({std::abs(T[0]), std::abs(T[1]),
                                     std::abs(T[2]), 1.0});
         scale = std::max(scale, sc);
         const real_t err = std::max({std::abs(trac(2 * k) - ex_dip),
                                      std::abs(trac(2 * k + 1) - ex_str),
                                      std::abs(sn(k) - ex_sn)});
         max_rel_all = std::max(max_rel_all, err);
         if (x[2] > -dz)
         {
            n_surf++;
            max_rel_surf = std::max(max_rel_surf, err);
            if (sample < 0 && std::abs(x[2]) < 1e-9 * Lz && x[0] > 0.0)
            {
               sample = k;
            }
         }
         if (k == ndof - 1 && sample >= 0)
         {
            const int s = sample;
            const real_t xs[3] = {xyz(3 * s), xyz(3 * s + 1), xyz(3 * s + 2)};
            real_t Ss[3][3];
            if (p >= 2) { ManufacturedStress(xs, lam, mu, Ss); }
            else { ManufacturedStressLinear(lam, mu, Ss); }
            const FaultBasisData &bs = fb->GetBasis(s / nbf);
            real_t Ts[3];
            for (int i = 0; i < 3; i++)
            {
               Ts[i] = Ss[i][0] * bs.normal[0] + Ss[i][1] * bs.normal[1]
                       + Ss[i][2] * bs.normal[2];
            }
            std::cout << std::setprecision(15)
                      << "  sample DOF on Z=0 at (" << xs[0] << ", " << xs[1]
                      << ", " << xs[2] << ")\n"
                      << "    computed: tau_dip = " << trac(2 * s)
                      << ", tau_strike = " << trac(2 * s + 1)
                      << ", sigma_n = " << sn(s) << "\n"
                      << "    analytic: tau_dip = " << Dot(Ts, bs.tangent1)
                      << ", tau_strike = " << Dot(Ts, bs.tangent2)
                      << ", sigma_n = " << -Dot(Ts, bs.normal) << "\n";
         }
      }
      // Field-scale-relative metric: max absolute error over the global
      // traction scale (a per-DOF denominator blows up where T ~ 0).
      std::cout << "  max abs err (surface) = " << max_rel_surf
                << " Pa, max abs err (all) = " << max_rel_all << " Pa\n";
      max_rel_surf /= scale;
      max_rel_all /= scale;
      std::cout << "  DOFs with z > -dz = " << n_surf
                << ", max rel err (surface) = " << max_rel_surf
                << ", max rel err (all " << ndof << ") = " << max_rel_all
                << ", traction scale = " << scale << " Pa\n";
      TEST_ASSERT(sample >= 0, "found a sample DOF on the Z=0 edge");
      TEST_ASSERT(max_rel_surf < 1e-9,
                  "traction matches analytic at surface-adjacent DOFs");
      TEST_ASSERT(max_rel_all < 1e-9, "traction matches analytic at all DOFs");
   }

   // ------------------------------------------------------------------ G0.5
   std::cout << "\n=== [IP p=" << p << "] G0.5: stored-normal orientation ===\n";
   {
      const FaultBasis *fb = op.GetFaultBasis();
      const real_t ref[3] = {0.0, -1.0, 0.0};
      int n_anti = 0, n_anti_left = 0, n_aligned = 0, n_aligned_right = 0;
      for (int fi = 0; fi < fb->NumFaces(); fi++)
      {
         const FaultBasisData &b = fb->GetBasis(fi);
         real_t c[3];
         Cross(b.tangent1, b.tangent2, c);
         bool left = true, right = true;
         for (int d = 0; d < 3; d++)
         {
            if (std::abs(c[d] + b.normal[d]) > 1e-12) { left = false; }
            if (std::abs(c[d] - b.normal[d]) > 1e-12) { right = false; }
         }
         if (Dot(b.normal, ref) < 0.0) { n_anti++; if (left) { n_anti_left++; } }
         else { n_aligned++; if (right) { n_aligned_right++; } }
      }
      std::cout << "  fault faces = " << fb->NumFaces()
                << ", anti-aligned with (0,-1,0) = " << n_anti
                << " (left-handed: " << n_anti_left << ")"
                << ", aligned = " << n_aligned
                << " (right-handed: " << n_aligned_right << ")\n";
      TEST_ASSERT(n_anti_left == n_anti,
                  "every anti-aligned face has a left-handed triad");
      TEST_ASSERT(n_aligned_right == n_aligned,
                  "every aligned face has a right-handed triad");
   }

   delete mesh;
}

// ---------------------------------------------------------------------------
// Deck check: tpv102_spatial.toml with [stress] replaced by constant_tensor.
// ---------------------------------------------------------------------------
static void TestDeckMaterial()
{
   std::cout << "\n=== Deck: tpv102_spatial.toml + constant_tensor stress ===\n";
   const std::string src = "tpv102/configs/tpv102_spatial.toml";
   std::ifstream in(src);
   TEST_ASSERT(in.good(), "tracked deck readable: " << src);
   if (!in.good()) { return; }

   // Replace the [stress] block (header through the line before the next
   // [section]) with the constant_tensor block.  The tracked deck is not
   // modified; the edited copy goes to the system temp directory.
   std::ostringstream out;
   std::string line;
   bool in_stress = false, replaced = false;
   while (std::getline(in, line))
   {
      if (line.rfind("[stress]", 0) == 0)
      {
         in_stress = true;
         replaced = true;
         out << "[stress]\n"
             << "kind        = \"constant_tensor\"\n"
             << "sigma_xx_pa = 0.0\n"
             << "sigma_yy_pa = -120.0e6\n"
             << "sigma_zz_pa = 0.0\n"
             << "sigma_xy_pa = 75.0e6\n"
             << "sigma_yz_pa = 0.0\n"
             << "sigma_xz_pa = 0.0\n\n";
         continue;
      }
      if (in_stress && !line.empty() && line[0] == '[') { in_stress = false; }
      if (in_stress) { continue; }
      out << line << "\n";
   }
   TEST_ASSERT(replaced, "found [stress] block in deck");

   const std::filesystem::path tmp =
      std::filesystem::temp_directory_path()
      / "nitsche_gates_tpv102_constant_tensor.toml";
   {
      std::ofstream o(tmp);
      o << out.str();
   }
   std::cout << "  edited copy: " << tmp.string() << "\n";

   sp::SpatialFrictionConfig cfg = sp::LoadSpatialFrictionConfig(tmp.string());

   const real_t lam = cfg.material_fallback.lambda;
   const real_t mu = cfg.material_fallback.mu;
   const real_t rho = cfg.material_fallback.rho;
   const real_t cs = std::sqrt(mu / rho);
   const real_t cp = std::sqrt((lam + 2.0 * mu) / rho);
   std::cout << std::setprecision(12)
             << "  stress kind = "
             << (cfg.stress.kind == sp::StressSourceKind::ConstantTensor
                 ? "constant_tensor" : "OTHER")
             << ", sigma_xy = " << cfg.stress.sigma_xy_pa
             << ", sigma_yy = " << cfg.stress.sigma_yy_pa << "\n"
             << "  material kind = "
             << (cfg.material.kind == sp::MaterialKind::Constant ? "Constant" : "OTHER")
             << ", velocity.use_sidecar = " << cfg.velocity.use_sidecar << "\n"
             << "  resolved: lambda = " << lam << ", mu = " << mu
             << ", rho = " << rho << "\n"
             << "  resolved: cs = " << cs << ", cp = " << cp << "\n"
             << "  TPV102Params: cs = " << TPV102Params::cs
             << ", cp = " << TPV102Params::cp
             << ", rho = " << TPV102Params::rho << "\n";

   TEST_ASSERT(cfg.stress.kind == sp::StressSourceKind::ConstantTensor,
               "stress kind parsed as constant_tensor");
   TEST_ASSERT(cfg.material.kind == sp::MaterialKind::Constant,
               "material kind is Constant ([material_constant_fallback])");
   TEST_NEAR(rho, TPV102Params::rho, 1e-9, "rho == TPV102Params::rho");
   TEST_NEAR(cs, TPV102Params::cs, 1e-6, "cs == TPV102Params::cs");
   TEST_NEAR(cp, TPV102Params::cp, 1e-6, "cp == TPV102Params::cp");
}

int main()
{
   std::cout << "=== Nitsche P0 gates (G0.2, G0.5) ===\n";

   const int nx = 12, ny = 12, nz = 6;
   const real_t Lx = 30.0e3, Ly = 30.0e3, Lz = 30.0e3;
   const real_t dx = 2.0 * Lx / nx, dz = Lz / nz;
   const int order = 2;

   Mesh *mesh = CreateTestMesh(nx, ny, nz, Lx, Ly, Lz);
   const int n_tagged = TagFaultPlaneY0(*mesh, 1e-6 * Ly);
   std::cout << "  mesh: " << mesh->GetNE() << " hexes, " << mesh->GetNBE()
             << " bdr elements, " << n_tagged << " faces tagged attr 3\n";

   std::set<int> tagged_faces;
   for (int be = 0; be < mesh->GetNBE(); be++)
   {
      if (mesh->GetBdrAttribute(be) == 3)
      {
         tagged_faces.insert(mesh->GetBdrElementFaceIndex(be));
      }
   }

   LinearElastic model(TPV102Params::lambda, TPV102Params::mu);
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.dirichlet_attrs = {5};

   ElasticityDomainOperator<Mesh> op(*mesh, order, model,
                                     0.0, Lz, 2.0 * Lx, bc,
                                     DGMethod::BR2, SolverType::CG_AMG);

   const Array<int> &ffaces = op.GetFaultInteriorFaces();
   const int nfault = ffaces.Size();
   const int nbf = op.GetNbfPerFace();
   const int ndof = op.GetNumFaultDOFs();

   // ---------------------------------------------------------------- G0.2-1
   std::cout << "\n=== G0.2-1: tagged faces classify as fault ===\n";
   {
      std::set<int> fault_set(ffaces.begin(), ffaces.end());
      std::cout << "  tagged = " << tagged_faces.size()
                << ", GetFaultInteriorFaces = " << nfault
                << ", GetFaultSharedFaces = " << op.GetFaultSharedFaces().Size()
                << "\n";
      TEST_ASSERT(static_cast<int>(tagged_faces.size()) == nx * nz,
                  "tagged face count == nx*nz = " << nx * nz);
      TEST_ASSERT(fault_set == tagged_faces,
                  "fault face set == tagged attr-3 face set");
      TEST_ASSERT(nfault == static_cast<int>(fault_set.size()),
                  "no duplicate fault faces");
      TEST_ASSERT(op.GetFaultSharedFaces().Size() == 0,
                  "serial mesh: no shared fault faces");
   }

   // ---------------------------------------------------------------- G0.2-2
   std::cout << "\n=== G0.2-2: fault DOF count vs tagged area ===\n";
   real_t area = 0.0;
   {
      for (int fi = 0; fi < nfault; fi++)
      {
         auto *T = mesh->GetFaceTransformation(ffaces[fi]);
         const IntegrationRule &ir = IntRules.Get(T->GetGeometryType(), 2);
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            T->SetIntPoint(&ir.IntPoint(q));
            area += ir.IntPoint(q).weight * T->Weight();
         }
      }
      const real_t area_exact = 2.0 * Lx * Lz;
      std::cout << std::setprecision(15)
                << "  nbf_per_face = " << nbf << ", num_fault_dofs = " << ndof
                << ", area = " << area << " m^2 (exact " << area_exact << ")\n";
      TEST_NEAR(area, area_exact, 1e-9 * area_exact, "fault area == 2*Lx*Lz");
      TEST_ASSERT(ndof == nbf * static_cast<int>(std::lround(area / (dx * dz))),
                  "num_fault_dofs == nbf * area/(dx*dz)");
   }

   // ---------------------------------------------------------------- G0.2-3
   std::cout << "\n=== G0.2-3: per-DOF fault frame ===\n";
   {
      DenseMatrix B;
      op.GetFaultDOFBasis(B);
      Vector xyz;
      op.GetFaultDOFCoords3D(xyz);

      real_t max_orth = 0.0;
      for (int k = 0; k < ndof; k++)
      {
         real_t n[3], t1[3], t2[3];
         for (int d = 0; d < 3; d++)
         {
            n[d] = B(d, k); t1[d] = B(3 + d, k); t2[d] = B(6 + d, k);
         }
         max_orth = std::max({max_orth,
                              std::abs(Dot(n, n) - 1.0),
                              std::abs(Dot(t1, t1) - 1.0),
                              std::abs(Dot(t2, t2) - 1.0),
                              std::abs(Dot(n, t1)), std::abs(Dot(n, t2)),
                              std::abs(Dot(t1, t2))});
      }
      std::cout << "  max orthonormality defect = " << max_orth << "\n";
      TEST_ASSERT(max_orth < 1e-12, "every fault DOF frame is orthonormal");

      // Reference: the fault face whose centroid is closest to (0, 0, -Lz/2).
      int ref_fi = -1;
      real_t best = 1e300;
      for (int fi = 0; fi < nfault; fi++)
      {
         real_t c[3];
         FaceCentroid(*mesh, ffaces[fi], c);
         const real_t d2 = c[0] * c[0] + (c[2] + 0.5 * Lz) * (c[2] + 0.5 * Lz);
         if (d2 < best) { best = d2; ref_fi = fi; }
      }
      real_t rc[3];
      FaceCentroid(*mesh, ffaces[ref_fi], rc);
      const int rk = ref_fi * nbf;
      real_t rn[3], rt1[3], rt2[3];
      for (int d = 0; d < 3; d++)
      {
         rn[d] = B(d, rk); rt1[d] = B(3 + d, rk); rt2[d] = B(6 + d, rk);
      }
      std::cout << "  mid-plane ref face centroid = (" << rc[0] << ", " << rc[1]
                << ", " << rc[2] << ")\n"
                << "    n  = (" << rn[0] << ", " << rn[1] << ", " << rn[2] << ")\n"
                << "    t1 = (" << rt1[0] << ", " << rt1[1] << ", " << rt1[2] << ")\n"
                << "    t2 = (" << rt2[0] << ", " << rt2[1] << ", " << rt2[2] << ")\n";

      // DOFs on faces touching the Z=0 edge.
      int n_edge = 0, n_plus = 0, n_minus = 0;
      real_t max_dev = 0.0;
      for (int fi = 0; fi < nfault; fi++)
      {
         Array<int> fv;
         mesh->GetFaceVertices(ffaces[fi], fv);
         real_t zmax = -1e300;
         for (int i = 0; i < fv.Size(); i++)
         {
            zmax = std::max(zmax, mesh->GetVertex(fv[i])[2]);
         }
         if (std::abs(zmax) > 1e-6 * Lz) { continue; }
         for (int kk = 0; kk < nbf; kk++)
         {
            const int k = fi * nbf + kk;
            real_t n[3], t1[3], t2[3];
            for (int d = 0; d < 3; d++)
            {
               n[d] = B(d, k); t1[d] = B(3 + d, k); t2[d] = B(6 + d, k);
            }
            const real_t s = (Dot(n, rn) >= 0.0) ? 1.0 : -1.0;
            (s > 0 ? n_plus : n_minus)++;
            n_edge++;
            for (int d = 0; d < 3; d++)
            {
               max_dev = std::max({max_dev, std::abs(n[d] - s * rn[d]),
                                   std::abs(t1[d] - s * rt1[d]),
                                   std::abs(t2[d] - s * rt2[d])});
            }
         }
      }
      std::cout << "  Z=0-edge fault DOFs = " << n_edge
                << " (sign +1: " << n_plus << ", sign -1: " << n_minus
                << "), max |frame - s*ref| = " << max_dev << "\n";
      TEST_ASSERT(n_edge == nx * nbf, "Z=0-edge DOF count == nx*nbf");
      TEST_ASSERT(max_dev < 1e-12,
                  "Z=0-edge frames == mid-plane frame up to per-face sign");
   }

   // ---------------------------------------------------------------- G0.2-4
   std::cout << "\n=== G0.2-4: ComputeTraction vs analytic (manufactured u) ===\n";
   {
      GridFunction u(&op.GetFESpace());
      VectorFunctionCoefficient ucoef(3, ManufacturedU);
      u.ProjectCoefficient(ucoef);

      Vector slip(2 * ndof);
      slip = 0.0;
      Vector trac, sn;
      op.ComputeTraction(u, slip, trac, &sn);

      const FaultBasis *fb = op.GetFaultBasis();
      const real_t lam = TPV102Params::lambda, mu = TPV102Params::mu;
      real_t max_rel_all = 0.0, max_rel_surf = 0.0, scale = 0.0;
      int n_surf = 0;
      for (int fi = 0; fi < nfault; fi++)
      {
         real_t c[3];
         FaceCentroid(*mesh, ffaces[fi], c);
         real_t S[3][3];
         ManufacturedStress(c, lam, mu, S);
         const FaultBasisData &b = fb->GetBasis(fi);
         real_t T[3];
         for (int i = 0; i < 3; i++)
         {
            T[i] = S[i][0] * b.normal[0] + S[i][1] * b.normal[1]
                   + S[i][2] * b.normal[2];
         }
         const real_t ex_dip = Dot(T, b.tangent1);
         const real_t ex_str = Dot(T, b.tangent2);
         const real_t ex_sn = -Dot(T, b.normal);
         const real_t sc = std::max({std::abs(T[0]), std::abs(T[1]),
                                     std::abs(T[2]), 1.0});
         scale = std::max(scale, sc);
         // BR2: nbf == 1, so DOF index == face index.
         const real_t err = std::max({std::abs(trac(2 * fi) - ex_dip),
                                      std::abs(trac(2 * fi + 1) - ex_str),
                                      std::abs(sn(fi) - ex_sn)});
         max_rel_all = std::max(max_rel_all, err);
         if (c[2] > -dz)
         {
            n_surf++;
            max_rel_surf = std::max(max_rel_surf, err);
         }
      }

      // Print one surface-adjacent sample (the face at x centroid ~ +dx/2).
      for (int fi = 0; fi < nfault; fi++)
      {
         real_t c[3];
         FaceCentroid(*mesh, ffaces[fi], c);
         if (!(c[2] > -dz && std::abs(c[0] - 0.5 * dx) < 1e-6 * Lx)) { continue; }
         real_t S[3][3];
         ManufacturedStress(c, lam, mu, S);
         const FaultBasisData &b = fb->GetBasis(fi);
         real_t T[3];
         for (int i = 0; i < 3; i++)
         {
            T[i] = S[i][0] * b.normal[0] + S[i][1] * b.normal[1]
                   + S[i][2] * b.normal[2];
         }
         std::cout << std::setprecision(15)
                   << "  surface-adjacent face centroid (" << c[0] << ", "
                   << c[1] << ", " << c[2] << ")\n"
                   << "    computed: tau_dip = " << trac(2 * fi)
                   << ", tau_strike = " << trac(2 * fi + 1)
                   << ", sigma_n = " << sn(fi) << "\n"
                   << "    analytic: tau_dip = " << Dot(T, b.tangent1)
                   << ", tau_strike = " << Dot(T, b.tangent2)
                   << ", sigma_n = " << -Dot(T, b.normal) << "\n";
         break;
      }
      // Field-scale-relative metric: max absolute error over the global
      // traction scale (a per-DOF denominator blows up where T ~ 0).
      std::cout << "  max abs err (surface) = " << max_rel_surf
                << " Pa, max abs err (all) = " << max_rel_all << " Pa\n";
      max_rel_surf /= scale;
      max_rel_all /= scale;
      std::cout << "  surface-adjacent faces = " << n_surf
                << ", max rel err (surface) = " << max_rel_surf
                << ", max rel err (all " << nfault << ") = " << max_rel_all
                << ", traction scale = " << scale << " Pa\n";
      TEST_ASSERT(n_surf == nx, "surface-adjacent fault faces == nx");
      TEST_ASSERT(max_rel_surf < 1e-9,
                  "traction matches analytic at surface-adjacent faces");
      TEST_ASSERT(max_rel_all < 1e-9, "traction matches analytic on all faces");
   }

   // ------------------------------------------------------------------ G0.5
   std::cout << "\n=== G0.5: stored-normal orientation and handedness ===\n";
   {
      const FaultBasis *fb = op.GetFaultBasis();
      const real_t ref[3] = {0.0, -1.0, 0.0};
      int n_anti = 0, n_anti_left = 0, n_aligned = 0, n_aligned_right = 0;
      for (int fi = 0; fi < fb->NumFaces(); fi++)
      {
         const FaultBasisData &b = fb->GetBasis(fi);
         real_t c[3];
         Cross(b.tangent1, b.tangent2, c);
         if (Dot(b.normal, ref) < 0.0)
         {
            n_anti++;
            bool left = true;
            for (int d = 0; d < 3; d++)
            {
               if (std::abs(c[d] + b.normal[d]) > 1e-12) { left = false; }
            }
            if (left) { n_anti_left++; }
         }
         else
         {
            n_aligned++;
            bool right = true;
            for (int d = 0; d < 3; d++)
            {
               if (std::abs(c[d] - b.normal[d]) > 1e-12) { right = false; }
            }
            if (right) { n_aligned_right++; }
         }
      }
      std::cout << "  fault faces = " << fb->NumFaces()
                << ", anti-aligned with (0,-1,0) = " << n_anti
                << " (left-handed t1 x t2 = -n: " << n_anti_left << ")"
                << ", aligned = " << n_aligned
                << " (right-handed t1 x t2 = +n: " << n_aligned_right << ")\n";
      TEST_ASSERT(n_anti_left == n_anti,
                  "every anti-aligned face has a left-handed triad");
   }

   delete mesh;

   for (int p = 1; p <= 3; p++) { RunTetIP(p); }

   TestDeckMaterial();

   std::cout << "\n=== Summary: " << num_passed << "/" << num_tests
             << " passed, " << num_failed << " failed ===\n";
   return (num_failed == 0) ? 0 : 1;
}

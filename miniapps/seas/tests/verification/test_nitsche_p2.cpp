// Nitsche dynamic rupture, phase P2 (plan C-02, test T-02): NitscheMass.
//
//   T-02(a): contrast.  Two tets sharing a face, rho_1 = 2670, rho_2 = 1000,
//            p = 1, 2, 3, random r:  ||M * Apply(r) - r|| / ||r|| < 1e-12,
//            with M assembled independently by MFEM (BilinearForm +
//            VectorMassIntegrator over a PWConstCoefficient).
//   T-02(b): bit-identity.  rho_1 = rho_2 = rho: Apply(r) equals the
//            separately computed (1/rho) * M_e^{-1} r_c, 0 differing entries.
//   Layout:  on the P0/P1 24-tet box (IP fault), NitscheMass's vector size and
//            per-element dofs match ElasticityDomainOperator's vector space,
//            and T-02(a) holds there with a per-element rho (contrast <= 10x)
//            sampled from a MaterialField through NitscheMaterialSample.
//
// Run at np = 1 and np = 2.

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../domain/nitsche_mass.hpp"
#include "../../config/tpv102_params.hpp"

#include <array>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static real_t GlobalNorm(const Vector &v, MPI_Comm comm)
{
   real_t loc = v * v, glob = 0.0;
   MPI_Allreduce(&loc, &glob, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);
   return std::sqrt(glob);
}

static int GlobalSum(int v, MPI_Comm comm)
{
   int g = 0;
   MPI_Allreduce(&v, &g, 1, MPI_INT, MPI_SUM, comm);
   return g;
}

// Two tets sharing the face (1,2,3); attributes 1 and 2.
static Mesh TwoTetMesh()
{
   Mesh mesh(3, 5, 2, 0);
   const real_t v[5][3] = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
                           {1, 1, 1}};
   for (int i = 0; i < 5; i++) { mesh.AddVertex(v[i]); }
   const int t1[4] = {0, 1, 2, 3}, t2[4] = {1, 2, 3, 4};
   mesh.AddTet(t1, 1);
   mesh.AddTet(t2, 2);
   mesh.FinalizeTopology();
   mesh.Finalize(false, true);
   mesh.SetAttributes();
   return mesh;
}

// ||M a - r|| / ||r|| with M from MFEM's BilinearForm on the vector DG space.
static real_t IndependentResidual(ParMesh &pmesh, int p, Coefficient &rho,
                                  const Vector &r, const Vector &a)
{
   DG_FECollection fec(p, 3, BasisType::GaussLobatto);
   ParFiniteElementSpace vfes(&pmesh, &fec, 3, Ordering::byNODES);
   ParBilinearForm m(&vfes);
   m.AddDomainIntegrator(new VectorMassIntegrator(rho));
   m.Assemble();
   m.Finalize();
   MFEM_VERIFY(m.SpMat().Height() == r.Size(), "layout mismatch");
   Vector Ma(r.Size());
   m.SpMat().Mult(a, Ma);
   Ma -= r;
   return GlobalNorm(Ma, pmesh.GetComm()) / GlobalNorm(r, pmesh.GetComm());
}

static void TestTwoTet(int p)
{
   if (g_rank == 0) { std::cout << "\n=== Two-tet mesh, p = " << p << " ===\n"; }
   Mesh smesh = TwoTetMesh();
   ParMesh pmesh(MPI_COMM_WORLD, smesh);
   const int ne = pmesh.GetNE();

   // ---- T-02(a): contrast rho_1 = 2670, rho_2 = 1000.
   {
      Vector rho_attr(2);
      rho_attr(0) = 2670.0;
      rho_attr(1) = 1000.0;
      PWConstCoefficient rho_pw(rho_attr);
      Vector rho_e(ne);
      for (int e = 0; e < ne; e++)
      {
         rho_e(e) = rho_attr(pmesh.GetAttribute(e) - 1);
      }
      NitscheMass mass(pmesh, p, rho_e);
      Vector r(mass.VSize());
      r.Randomize(11 + p + 100 * g_rank);
      r -= 0.5;
      Vector a;
      mass.Apply(r, a);
      const real_t rel = IndependentResidual(pmesh, p, rho_pw, r, a);
      const int gsize = GlobalSum(mass.VSize(), pmesh.GetComm());
      if (g_rank == 0)
      {
         std::cout << std::setprecision(6) << "  T-02(a): global vector size "
                   << gsize
                   << ", ||M Apply(r) - r|| / ||r|| = " << rel << "\n";
      }
      TEST_ASSERT(rel < 1e-12, "T-02(a) two-tet contrast (2670 / 1000)");
   }

   // ---- T-02(b): equal densities, bitwise vs separate (1/rho) M_e^{-1} r_c.
   {
      const real_t rho = 2670.0;
      Vector rho_e(ne);
      rho_e = rho;
      NitscheMass mass(pmesh, p, rho_e);
      Vector r(mass.VSize());
      r.Randomize(21 + p + 100 * g_rank);
      r -= 0.5;
      Vector a;
      mass.Apply(r, a);

      // Separate computation: same scalar space, same rule, same inverse.
      DG_FECollection fec(p, 3, BasisType::GaussLobatto);
      ParFiniteElementSpace fes(&pmesh, &fec);
      const int n = fes.GetNDofs();
      const real_t inv_rho = 1.0 / rho;
      Vector ref(r.Size());
      Array<int> dofs;
      for (int e = 0; e < ne; e++)
      {
         const FiniteElement *fe = fes.GetFE(e);
         ElementTransformation *T = pmesh.GetElementTransformation(e);
         const int nd = fe->GetDof();
         DenseMatrix M(nd), Minv(nd);
         M = 0.0;
         Vector shape(nd);
         const IntegrationRule &ir = IntRules.Get(fe->GetGeomType(),
                                                  2 * fe->GetOrder());
         for (int j = 0; j < ir.GetNPoints(); j++)
         {
            const IntegrationPoint &ip = ir.IntPoint(j);
            T->SetIntPoint(&ip);
            fe->CalcShape(ip, shape);
            const real_t w = ip.weight * T->Weight();
            for (int k = 0; k < nd; k++)
               for (int l = 0; l < nd; l++) { M(k, l) += w * shape(k) * shape(l); }
         }
         DenseMatrixInverse Mi(M);
         Mi.GetInverseMatrix(Minv);
         fes.GetElementDofs(e, dofs);
         Vector x(nd), y(nd);
         for (int c = 0; c < 3; c++)
         {
            for (int k = 0; k < nd; k++) { x(k) = r(c * n + dofs[k]); }
            Minv.Mult(x, y);
            for (int k = 0; k < nd; k++) { ref(c * n + dofs[k]) = y(k) * inv_rho; }
         }
      }
      int ndiff = 0;
      for (int i = 0; i < a.Size(); i++)
      {
         if (std::memcmp(&a(i), &ref(i), sizeof(real_t)) != 0) { ndiff++; }
      }
      ndiff = GlobalSum(ndiff, pmesh.GetComm());
      const int gsize = GlobalSum(a.Size(), pmesh.GetComm());
      if (g_rank == 0)
      {
         std::cout << "  T-02(b): differing entries = " << ndiff << " of "
                   << gsize << "\n";
      }
      TEST_ASSERT(ndiff == 0, "T-02(b) equal density bit-identical");
   }
}

// ---------------------------------------------------------------------------
// 24-tet TPV102 box (as test_nitsche_p1.cpp), 6x6x3 hexes.
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

// Pseudo-random density in [1000, 10000] (contrast up to 10x).
static real_t RhoField(const Vector &x)
{
   const real_t s = std::sin(1.3e-3 * x(0) + 2.1e-3 * x(1) + 0.7e-3 * x(2)
                             + 1.7e-7 * x(0) * x(2));
   const real_t u = s * s;              // [0, 1]
   return 1000.0 + 9000.0 * u;
}

static void TestBox(ParMesh &pmesh, int p)
{
   if (g_rank == 0) { std::cout << "\n=== 24-tet box, IP, p = " << p << " ===\n"; }
   const real_t L = 30.0e3;
   LinearElastic model(TPV102Params::lambda, TPV102Params::mu);
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = [](const Vector &, real_t, Vector &u)
   { u.SetSize(3); u = 0.0; };
   ElasticityDomainOperator<ParMesh> op(pmesh, p, model, 0.0, L, 2.0 * L, bc,
                                        DGMethod::IP, SolverType::MUMPS);

   ConstantCoefficient lam_c(TPV102Params::lambda), mu_c(TPV102Params::mu);
   FunctionCoefficient rho_c(RhoField);
   MaterialField mat = MaterialField::MakeCoefficient(&lam_c, &mu_c, &rho_c);
   NitscheMaterialSample sample(pmesh, mat);
   Vector rho_e;
   RhoPerElement(sample, rho_e);
   real_t rmin = rho_e.Min(), rmax = rho_e.Max();
   MPI_Allreduce(MPI_IN_PLACE, &rmin, 1, MPITypeMap<real_t>::mpi_type, MPI_MIN,
                 pmesh.GetComm());
   MPI_Allreduce(MPI_IN_PLACE, &rmax, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX,
                 pmesh.GetComm());

   NitscheMass mass(pmesh, p, rho_e);

   // ---- Layout vs the operator's vector space.
   const auto &vfes = op.GetFESpace();
   const int n = mass.NumScalarDofs();
   bool size_ok = (vfes.GetVSize() == mass.VSize())
                  && (vfes.GetNDofs() == n)
                  && (vfes.GetOrdering() == Ordering::byNODES);
   int bad = 0;
   Array<int> vdofs;
   for (int e = 0; e < pmesh.GetNE(); e++)
   {
      vfes.GetElementVDofs(e, vdofs);
      const Array<int> &sd = mass.ElementScalarDofs(e);
      if (vdofs.Size() != 3 * sd.Size()) { bad++; continue; }
      for (int c = 0; c < 3; c++)
         for (int k = 0; k < sd.Size(); k++)
         {
            if (vdofs[c * sd.Size() + k] != c * n + sd[k]) { bad++; }
         }
   }
   bad = GlobalSum(bad, pmesh.GetComm());
   size_ok = GlobalSum(size_ok ? 0 : 1, pmesh.GetComm()) == 0;
   const int gop = GlobalSum(vfes.GetVSize(), pmesh.GetComm());
   const int gms = GlobalSum(mass.VSize(), pmesh.GetComm());
   if (g_rank == 0)
   {
      std::cout << "  operator VSize (global) = "
                << gop
                << ", NitscheMass VSize (global) = "
                << gms
                << ", per-element dof mismatches = " << bad << "\n";
   }
   TEST_ASSERT(size_ok, "vector size / scalar ndofs / byNODES match operator");
   TEST_ASSERT(bad == 0, "per-element dof lists match operator (c*ndofs + i)");

   // ---- T-02(a) on the box with the sampled rho.
   GridFunctionCoefficient rho_gfc(&sample.RhoGF());
   Vector r(mass.VSize());
   r.Randomize(31 + p + 100 * g_rank);
   r -= 0.5;
   Vector a;
   mass.Apply(r, a);
   const real_t rel = IndependentResidual(pmesh, p, rho_gfc, r, a);
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  sampled rho in [" << rmin << ", "
                << rmax << "] (contrast " << rmax / rmin << "x)"
                << ", ||M Apply(r) - r|| / ||r|| = " << rel << "\n";
   }
   TEST_ASSERT(rmax / rmin > 5.0, "sampled density contrast > 5x");
   TEST_ASSERT(rel < 1e-12, "T-02(a) on the box with sampled rho");
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   if (g_rank == 0)
   {
      std::cout << "=== Nitsche P2 (NitscheMass, T-02) on "
                << Mpi::WorldSize() << " rank(s) ===\n";
   }

   for (int p = 1; p <= 3; p++) { TestTwoTet(p); }

   Mesh *mesh = CreateTestMesh24Tet(6, 6, 3, 30.0e3, 30.0e3, 30.0e3);
   ParMesh pmesh(MPI_COMM_WORLD, *mesh);
   delete mesh;
   for (int p = 1; p <= 3; p++) { TestBox(pmesh, p); }

   if (g_rank == 0)
   {
      std::cout << "\n=== Summary: " << num_passed << "/" << num_tests
                << " passed, " << num_failed << " failed ===\n";
   }
   return (num_failed == 0) ? 0 : 1;
}

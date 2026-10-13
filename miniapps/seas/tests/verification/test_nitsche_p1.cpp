// Nitsche dynamic rupture, phase P1 (C-01 ApplyStiffness, C-09 spatially
// varying Lame coefficients, per-element material sample).
//
//   T-01:    Solve then ApplyStiffness reproduces the slip RHS:
//            ||K x - b|| / ||b|| < 1e-10 with b = G d, t = 0 (zero Dirichlet).
//   T-19(a): the Coefficient constructor fed by NitscheMaterialSample of a
//            constant MaterialField is BITWISE equal to the LinearElastic
//            constructor on K u, G d and ComputeTraction.
//   Penalty: per-side ComputePenalty matches the hand formula for
//            lam2 != lam1, and is bitwise equal to the former single-moduli
//            formula when both sides are equal.
//
// Mesh: the P0 TPV102 box (60 x 60 x 30 km, fault on Y=0, free surface at
// Z=0) split into 24 tets per hex, DGMethod::IP, p = 1, 2, 3.  To keep the
// MUMPS solves and p=3 assemblies laptop-sized the box uses 6x6x3 hexes
// (10 km cells, 5184 tets) instead of P0's 12x12x6.  Run with mpirun -np 1.

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
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

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

// 24-tet TPV102 box: same construction as test_nitsche_gates.cpp
// (CreateTestMesh24Tet): shift to [-Lx,Lx]x[-Ly,Ly]x[-Lz,0], attr 1 at Z=0,
// attr 5 elsewhere, attr 3 triangles on the interior Y=0 plane.
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

static BoundaryConfig MakeBC()
{
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.dirichlet_attrs = {5};
   // Zero far-field displacement (t = 0 loading).
   bc.default_dirichlet_func = [](const Vector &, real_t, Vector &u)
   {
      u.SetSize(3);
      u = 0.0;
   };
   return bc;
}

static real_t GlobalNorm(const Vector &v, MPI_Comm comm)
{
   real_t loc = v * v, glob = 0.0;
   MPI_Allreduce(&loc, &glob, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);
   return std::sqrt(glob);
}

static bool BitwiseEqual(const Vector &a, const Vector &b, int &n_diff,
                         real_t &max_abs)
{
   n_diff = 0;
   max_abs = 0.0;
   if (a.Size() != b.Size()) { n_diff = -1; return false; }
   for (int i = 0; i < a.Size(); i++)
   {
      if (std::memcmp(&a(i), &b(i), sizeof(real_t)) != 0)
      {
         n_diff++;
         max_abs = std::max(max_abs, std::abs(a(i) - b(i)));
      }
   }
   return n_diff == 0;
}

static const int kNx = 6, kNy = 6, kNz = 3;
static const real_t kL = 30.0e3;

// ---------------------------------------------------------------------- T-01
static void TestT01(ParMesh &pmesh, int p)
{
   std::cout << "\n=== T-01 (p=" << p << "): ||K x - G d|| / ||G d|| ===\n";
   LinearElastic model(TPV102Params::lambda, TPV102Params::mu);
   ElasticityDomainOperator<ParMesh> op(pmesh, p, model, 0.0, kL, 2.0 * kL,
                                        MakeBC(), DGMethod::IP,
                                        SolverType::MUMPS);
   const int nslip = 2 * op.GetNumFaultDOFs();
   Vector d(nslip);
   d.Randomize(1234 + p);
   d -= 0.5;   // nonzero, mixed-sign slip [m]

   Vector b;
   op.AssembleSlipOnlyRHS(b, d);

   ParGridFunction x(&op.GetFESpace());
   op.Solve(0.0, d, x);

   Vector Kx;
   op.ApplyStiffness(x, Kx);
   Vector r(Kx);
   r -= b;
   const real_t nb = GlobalNorm(b, pmesh.GetComm());
   const real_t rel = GlobalNorm(r, pmesh.GetComm()) / nb;
   std::cout << std::setprecision(6) << "  true dofs = " << x.Size()
             << ", fault dofs = " << op.GetNumFaultDOFs()
             << ", ||b|| = " << nb << ", ||Kx - b||/||b|| = " << rel << "\n";
   TEST_ASSERT(nb > 0.0, "slip RHS is nonzero");
   TEST_ASSERT(rel < 1e-10, "||K x - b|| / ||b|| < 1e-10");
}

// ------------------------------------------------------------------- T-19(a)
static void TestT19a(ParMesh &pmesh, int p)
{
   std::cout << "\n=== T-19(a) (p=" << p << "): constant material, old vs new "
             << "constructor, bitwise ===\n";
   const real_t lam = 32043759360.0, mu = 32038120320.0, rho = 2670.0;

   LinearElastic model(lam, mu);
   ElasticityDomainOperator<ParMesh> op_old(pmesh, p, model, 0.0, kL, 2.0 * kL,
                                            MakeBC(), DGMethod::IP,
                                            SolverType::MUMPS);

   MaterialField mat = MaterialField::MakeConstant(lam, mu, rho);
   NitscheMaterialSample sample(pmesh, mat);
   ElasticityDomainOperator<ParMesh> op_new(pmesh, p, sample.Lambda(),
                                            sample.Mu(), lam, mu,
                                            0.0, kL, 2.0 * kL, MakeBC(),
                                            DGMethod::IP, SolverType::MUMPS);

   // Random displacement and slip shared by both operators.
   ParGridFunction u(&op_old.GetFESpace());
   u.Randomize(77 + p);
   const int nslip = 2 * op_old.GetNumFaultDOFs();
   TEST_ASSERT(nslip == 2 * op_new.GetNumFaultDOFs(), "same fault DOF count");
   Vector d(nslip);
   d.Randomize(99 + p);
   d -= 0.5;

   int nd; real_t mx;

   Vector Ku_old, Ku_new;
   op_old.ApplyStiffness(u, Ku_old);
   op_new.ApplyStiffness(u, Ku_new);
   const bool eqK = BitwiseEqual(Ku_old, Ku_new, nd, mx);
   std::cout << "  K u:  size " << Ku_old.Size() << ", differing entries "
             << nd << ", max |diff| " << mx << "\n";
   TEST_ASSERT(eqK, "K u bitwise equal");

   Vector G_old, G_new;
   op_old.AssembleSlipOnlyRHS(G_old, d);
   op_new.AssembleSlipOnlyRHS(G_new, d);
   const bool eqG = BitwiseEqual(G_old, G_new, nd, mx);
   std::cout << "  G d:  size " << G_old.Size() << ", differing entries "
             << nd << ", max |diff| " << mx << "\n";
   TEST_ASSERT(eqG, "G d bitwise equal");

   Vector T_old, T_new, sn_old, sn_new;
   op_old.ComputeTraction(u, d, T_old, &sn_old);
   op_new.ComputeTraction(u, d, T_new, &sn_new);
   const bool eqT = BitwiseEqual(T_old, T_new, nd, mx);
   std::cout << "  traction (dip, strike): size " << T_old.Size()
             << ", differing entries " << nd << ", max |diff| " << mx << "\n";
   const bool eqN = BitwiseEqual(sn_old, sn_new, nd, mx);
   std::cout << "  normal traction: size " << sn_old.Size()
             << ", differing entries " << nd << ", max |diff| " << mx << "\n";
   TEST_ASSERT(eqT && eqN, "ComputeTraction bitwise equal");

   // The sample itself reproduces the constants exactly.
   const real_t lmin = sample.LambdaGF().Min(), lmax = sample.LambdaGF().Max();
   const real_t mmin = sample.MuGF().Min(), mmax = sample.MuGF().Max();
   const real_t rmin = sample.RhoGF().Min(), rmax = sample.RhoGF().Max();
   TEST_ASSERT(lmin == lam && lmax == lam && mmin == mu && mmax == mu
               && rmin == rho && rmax == rho,
               "per-element sample == (lambda, mu, rho) on every element");
}

// ------------------------------------------------------------------ Penalty
static void TestPenalty(ParMesh &pmesh, int p)
{
   std::cout << "\n=== Penalty per side (p=" << p << ") ===\n";
   const int dim = 3;
   L2_FECollection fec(p, dim, BasisType::GaussLobatto);
   ParFiniteElementSpace fes(&pmesh, &fec);
   const FiniteElement &fe = *fes.GetFE(0);

   ConstantCoefficient lc(1.0), mc(1.0);
   const real_t pf = 1.0;
   DGElasticityIPCombinedIntegrator integ(lc, mc, dim, -1.0, pf);

   const real_t detJ1 = 3.7e10, detJ2 = 5.1e10, nl = 2.4e7;
   const real_t lam1 = 32043759360.0, mu1 = 32038120320.0;
   const real_t lam2 = 2.0e10, mu2 = 1.5e10;

   // Hand formula, per side.
   const real_t cN = p * (p + dim - 1.0) / dim;
   auto side = [&](real_t detJ, real_t lam, real_t mu)
   {
      const real_t c0 = 2.0 * mu, c1 = dim * lam + 2.0 * mu;
      return (dim + 1) * cN * (real_t(dim) * nl / detJ) * (c1 * c1 / c0);
   };
   const real_t hand = pf * (side(detJ1, lam1, mu1) + side(detJ2, lam2, mu2))
                       / 4.0;
   const real_t got = integ.GetPenalty(fe, fe, detJ1, detJ2,
                                       lam1, mu1, lam2, mu2, nl);
   const real_t rel = std::abs(got - hand) / hand;
   std::cout << std::setprecision(17) << "  lam2 != lam1: got " << got
             << ", hand " << hand << ", rel diff " << rel << "\n";
   TEST_ASSERT(rel < 1e-14, "per-side penalty == hand formula");

   // Equal sides: bitwise equal to the former single-moduli formula,
   // reproduced here with its original operation order.
   const real_t c0 = 2.0 * mu1, c1 = dim * lam1 + 2.0 * mu1;
   const real_t ratio = c1 * c1 / c0;
   const real_t p0 = (dim + 1) * cN * (real_t(dim) * nl / detJ1) * ratio;
   const real_t p1 = (dim + 1) * cN * (real_t(dim) * nl / detJ2) * ratio;
   const real_t old_val = pf * (p0 + p1) / 4.0;
   const real_t eq9 = integ.GetPenalty(fe, fe, detJ1, detJ2,
                                       lam1, mu1, lam1, mu1, nl);
   const real_t eq7 = integ.GetPenalty(fe, fe, detJ1, detJ2, lam1, mu1, nl);
   std::cout << "  equal sides: per-side " << eq9 << ", legacy accessor "
             << eq7 << ", former formula " << old_val << "\n";
   TEST_ASSERT(std::memcmp(&eq9, &old_val, sizeof(real_t)) == 0,
               "equal sides: per-side penalty bitwise == former formula");
   TEST_ASSERT(std::memcmp(&eq7, &old_val, sizeof(real_t)) == 0,
               "legacy GetPenalty bitwise == former formula");
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();

   std::cout << "=== Nitsche P1 (C-01, C-09, material sample) ===\n";
   Mesh *mesh = CreateTestMesh24Tet(kNx, kNy, kNz, kL, kL, kL);
   std::cout << "  serial mesh: " << mesh->GetNE() << " tets, "
             << mesh->GetNBE() << " bdr elements\n";
   ParMesh pmesh(MPI_COMM_WORLD, *mesh);
   delete mesh;

   for (int p = 1; p <= 3; p++) { TestPenalty(pmesh, p); }
   for (int p = 1; p <= 3; p++) { TestT19a(pmesh, p); }
   for (int p = 1; p <= 3; p++) { TestT01(pmesh, p); }

   std::cout << "\n=== Summary: " << num_passed << "/" << num_tests
             << " passed, " << num_failed << " failed ===\n";
   return (num_failed == 0) ? 0 : 1;
}

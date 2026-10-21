// Poroelastic plan rev 3, phase P1 (steps P1.1 and P1.2): the Biot material
// sample, the closed-form wave speeds and the time-step bound.
//
//   U1-1  DeriveBiotConstants (derivation (4.4), (8.4)) for tpv5-poro, the
//         TPV205-equivalent set and the V-F3 contrast set, against
//         poro_test_refs.hpp (NumPy and SciPy): 1e-13 relative (class R).
//   U1-2  Identities, each side computed its own way, 1e-14 relative:
//         K_u = lambda_u + 2 mu / 3 against K_d + alpha^2 M; B = alpha M / K_u
//         against (K_u - K_d) / (alpha K_u); lambda_u - lambda against
//         alpha^2 M; Gassmann's K_u; tau against (m - rho_f^2 / rho) / b;
//         SeisSol's M against 1 / (phi / K_f + (alpha - phi) / K_s).
//   U1-3  Limits: K_s = 1e30 gives alpha = 1 and M = K_f / phi; the
//         TPV205-equivalent and V-F3 lambda_u equal rho (c_p^2 - 2 c_s^2) of
//         TPV205 (1e-15).
//   U1-4  PoroMaterialSample on F-C (tpv5-poro for x < 0, TPV205-equivalent
//         for x > 0): every element and every face neighbour holds its side's
//         constants bit for bit; the undrained sample's (lambda, mu, rho) are
//         the sample's (lambda_u, mu, rho) bit for bit; the deck path
//         (PoroConfig) on F-B gives the same constants as DeriveBiotConstants.
//   U1-5  ComputeBiotSpeeds (derivation (4.6)): fast and slow P, S at high
//         and low frequency, undrained and drained P against
//         poro_test_refs.hpp, 1e-13 relative.
//   U1-6  PoroMaxDt (derivation (9.6)) on F-B, F-C and, with --mesh1000, the
//         gmsh 1000 m TPV205 mesh, p = 1-3, against (9.6) evaluated here on
//         the serial mesh: bit for bit on every rank count.
//
// Usage: mpirun -np 1|2|4 ./seas_test_poro_p1 [--mesh1000 <tpv2053d_1000m.msh>]

#include "mfem.hpp"
#include "../poro_material.hpp"
#include "../poro_dispersion.hpp"
#include "../poro_cfl.hpp"
#include "poro_fixtures.hpp"
#include "poro_test_refs.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;
using namespace mfem::seas::poro;

/// One material set of poro_test_refs.hpp: its ten inputs and expected values.
struct Ref
{
   std::string name;
   BiotFluidInputs f;
   real_t lam, mu;
   real_t Kd, alpha, M, lam_u, Ku, B, rho, m, b, tau;
   real_t c_fast, c_slow, c_s_hf, c_s_lf, c_p_undrained, c_p_drained;
};

#define PORO_REF(NS) Ref{#NS, BiotFluidInputs{NS::phi, NS::kappa, NS::eta, NS::T, \
   NS::rho_f, NS::rho_s, NS::Kf, NS::Ks}, NS::lam, NS::mu, NS::Kd, NS::alpha, NS::M, \
   NS::lam_u, NS::Ku, NS::B, NS::rho, NS::m, NS::b, NS::tau, NS::c_fast, NS::c_slow, \
   NS::c_s_hf, NS::c_s_lf, NS::c_p_undrained, NS::c_p_drained}

static std::vector<Ref> Refs()
{
   using namespace poro_refs;
   return {PORO_REF(tpv5_poro), PORO_REF(tpv205_equivalent), PORO_REF(vf3_contrast)};
}

static BiotConstants Derive(const Ref &r)
{
   return DeriveBiotConstants(r.lam, r.mu, BulkDensity(r.f), r.f);
}

static real_t Rel(real_t a, real_t b) { return std::abs(a - b) / std::abs(b); }

static bool SameBits(real_t a, real_t b) { return std::memcmp(&a, &b, sizeof(real_t)) == 0; }

struct Named { const char *name; real_t got, expected; };

/// Prints the relative error of each entry and returns the worst.
static real_t Worst(const std::string &label, const std::vector<Named> &v, std::string &which)
{
   real_t worst = 0.0;
   which = v.front().name;
   std::ostringstream line;
   line.precision(2);
   for (const Named &n : v)
   {
      const real_t e = Rel(n.got, n.expected);
      line << " " << n.name << " " << std::scientific << e;
      if (!(e <= worst)) { worst = e; which = n.name; }
   }
   if (g_rank == 0) { std::cout << "  " << label << ": rel err" << line.str() << "\n"; }
   return worst;
}

// ------------------------------------------------------------------- U1-1
static void TestU11()
{
   for (const Ref &r : Refs())
   {
      TEST_ASSERT(CheckBiotInputs(r.lam, r.mu, BulkDensity(r.f), r.f).empty(),
                  "U1-1 " << r.name << ": the inputs make a valid Biot medium");
      const BiotConstants c = Derive(r);
      std::string which;
      const real_t worst = Worst("U1-1 " + r.name,
         {{"K_d", c.Kd, r.Kd}, {"alpha", c.alpha, r.alpha}, {"M", c.M, r.M},
          {"lambda_u", c.lambda_u, r.lam_u}, {"K_u", c.Ku, r.Ku}, {"B", c.B, r.B},
          {"rho", c.rho, r.rho}, {"m", c.m, r.m}, {"b", c.b, r.b}, {"tau", c.tau, r.tau}},
         which);
      TEST_ASSERT(worst <= 1e-13, "U1-1 " << r.name << ": K_d, alpha, M, lambda_u, K_u, B, "
                  "rho, m, b, tau match the reference to 1e-13 (worst " << which << " "
                  << worst << ")");
   }
}

// ------------------------------------------------------------------- U1-2
static void TestU12()
{
   for (const Ref &r : Refs())
   {
      const BiotConstants c = Derive(r);
      const BiotFluidInputs &f = r.f;
      const real_t phi = f.porosity, Kf = f.bulk_fluid, Ks = f.bulk_solid;
      const std::string n = "U1-2 " + r.name + ": ";
      real_t e;
      e = Rel(c.lambda_u + 2.0 / 3.0 * c.mu, c.Kd + c.alpha * c.alpha * c.M);
      TEST_ASSERT(e <= 1e-14, n << "K_u = lambda_u + 2 mu / 3 = K_d + alpha^2 M (" << e << ")");
      e = Rel(c.alpha * c.M / c.Ku, (c.Ku - c.Kd) / (c.alpha * c.Ku));
      TEST_ASSERT(e <= 1e-14, n << "B = alpha M / K_u = (K_u - K_d) / (alpha K_u) (" << e << ")");
      e = Rel(c.lambda_u - c.lambda, c.alpha * c.alpha * c.M);
      TEST_ASSERT(e <= 1e-14, n << "lambda_u - lambda = alpha^2 M (" << e << ")");
      const real_t Kd = c.Kd;
      const real_t Ku_gassmann = Kd + (1.0 - Kd / Ks) * (1.0 - Kd / Ks)
                                 / (phi / Kf + (1.0 - phi) / Ks - Kd / (Ks * Ks));
      e = Rel(c.Ku, Ku_gassmann);
      TEST_ASSERT(e <= 1e-14, n << "K_u equals Gassmann's K_u (" << e << ")");
      e = Rel(c.tau, (c.m - c.rho_f * c.rho_f / c.rho) / c.b);
      TEST_ASSERT(e <= 1e-14, n << "tau = (m - rho_f^2 / rho) / b (" << e << ")");
      e = Rel(c.M, 1.0 / (phi / Kf + (c.alpha - phi) / Ks));
      TEST_ASSERT(e <= 1e-14, n << "M = 1 / (phi / K_f + (alpha - phi) / K_s) (" << e << ")");
   }
}

// ------------------------------------------------------------------- U1-3
static void TestU13()
{
   const Ref r = Refs()[0];
   BiotFluidInputs f = r.f;
   f.bulk_solid = 1e30;
   const BiotConstants c = DeriveBiotConstants(r.lam, r.mu, BulkDensity(f), f);
   TEST_ASSERT(std::abs(c.alpha - 1.0) <= 1e-15,
               "U1-3 K_s = 1e30: alpha = 1 (|alpha - 1| " << std::abs(c.alpha - 1.0) << ")");
   const real_t e = Rel(c.M, f.bulk_fluid / f.porosity);
   TEST_ASSERT(e <= 1e-15, "U1-3 K_s = 1e30: M = K_f / phi (" << e << ")");
   // TPV205: rho 2670 kg/m^3, c_s 3464 m/s, c_p 6000 m/s (plan section 3.2)
   const real_t lam_tpv205 = 2670.0 * (6000.0 * 6000.0 - 2.0 * 3464.0 * 3464.0);
   const std::vector<Ref> R = Refs();
   for (int k = 1; k <= 2; k++)
   {
      const Ref &s = R[k];
      const real_t ek = Rel(Derive(s).lambda_u, lam_tpv205);
      TEST_ASSERT(ek <= 1e-15, "U1-3 " << s.name << ": lambda_u = rho (c_p^2 - 2 c_s^2) of "
                  "TPV205 (" << ek << ")");
   }
}

// ------------------------------------------------------------------- U1-4
/// Number of values that differ in bits from DeriveBiotConstants of the side.
static long CheckElement(const PoroMaterialSample &s, int e, const BiotConstants &c)
{
   Array<int> dofs;
   s.AlphaGF().ParFESpace()->GetElementDofs(e, dofs);
   const int k = dofs[0];
   Array<int> ddofs;
   s.DrainedSample().RhoGF().ParFESpace()->GetElementDofs(e, ddofs);
   const int dk = ddofs[0];
   const BiotConstants &h = s.Constants(e);
   const real_t got[] = {s.DrainedSample().LambdaGF()(dk), s.DrainedSample().MuGF()(dk),
                         s.DrainedSample().RhoGF()(dk), s.RhoFluidGF()(k), s.FluidMassGF()(k),
                         s.DragGF()(k), s.AlphaGF()(k), s.BiotModulusGF()(k), s.TauGF()(k),
                         s.LambdaUGF()(k), h.lambda, h.mu, h.rho, h.rho_f, h.Kd, h.alpha, h.M,
                         h.lambda_u, h.Ku, h.B, h.m, h.b, h.tau};
   const real_t exp[] = {c.lambda, c.mu, c.rho, c.rho_f, c.m, c.b, c.alpha, c.M, c.tau,
                         c.lambda_u, c.lambda, c.mu, c.rho, c.rho_f, c.Kd, c.alpha, c.M,
                         c.lambda_u, c.Ku, c.B, c.m, c.b, c.tau};
   long n = 0;
   for (size_t i = 0; i < sizeof(got) / sizeof(got[0]); i++)
   {
      if (!SameBits(got[i], exp[i])) { n++; }
   }
   return n;
}

static void TestU14()
{
   const std::vector<Ref> R = Refs();
   const Ref &A = R[0], &B = R[1];
   const BiotConstants cA = Derive(A), cB = Derive(B);

   std::unique_ptr<Mesh> m(TPV102Box(4, 4, 2));
   MarkSidesFC(*m);
   ParMesh pm(MPI_COMM_WORLD, *m);
   Vector lv(2), mv(2), rv(2);
   lv(0) = A.lam; lv(1) = B.lam;
   mv(0) = A.mu;  mv(1) = B.mu;
   rv(0) = BulkDensity(A.f); rv(1) = BulkDensity(B.f);
   PWConstCoefficient lamc(lv), muc(mv), rhoc(rv);
   MaterialField mf;
   mf.mode = MaterialField::Mode::Coefficient;
   mf.lambda_coef = &lamc;
   mf.mu_coef = &muc;
   mf.rho_coef = &rhoc;
   PoroMaterialSample s(pm, mf, [&](int, const Vector &x) { return x(0) < 0.0 ? A.f : B.f; });

   long nA = 0, nB = 0, bad = 0;
   Vector x(3);
   for (int e = 0; e < pm.GetNE(); e++)
   {
      pm.GetElementTransformation(e)->Transform(
         Geometries.GetCenter(pm.GetElementBaseGeometry(e)), x);
      const bool a = x(0) < 0.0;
      (a ? nA : nB)++;
      bad += CheckElement(s, e, a ? cA : cB);
   }
   nA = LSum(nA, pm.GetComm());
   nB = LSum(nB, pm.GetComm());
   bad = LSum(bad, pm.GetComm());
   TEST_ASSERT(nA == 384 && nB == 384, "U1-4 F-C has 384 elements on each side ("
               << nA << ", " << nB << ")");
   TEST_ASSERT(bad == 0, "U1-4 every element holds its side's constants bit for bit, in the "
               "grid functions and in Constants(e) (" << bad << " values differ)");

   // Face neighbours (other ranks' elements across shared faces).
   long nnbr = 0, nnbrA = 0, badnbr = 0;
   const ParGridFunction *gfs[] = {&s.AlphaGF(), &s.BiotModulusGF(), &s.LambdaUGF(),
                                   &s.RhoFluidGF(), &s.FluidMassGF(), &s.DragGF(),
                                   &s.TauGF(), &s.UndrainedSample().LambdaGF(),
                                   &s.UndrainedSample().MuGF(), &s.UndrainedSample().RhoGF()};
   for (int j = 0; j < pm.GetNFaceNeighborElements(); j++)
   {
      pm.GetFaceNbrElementTransformation(j)->Transform(
         Geometries.GetCenter(Geometry::TETRAHEDRON), x);
      const BiotConstants &c = x(0) < 0.0 ? cA : cB;
      if (x(0) < 0.0) { nnbrA++; }
      const real_t exp[] = {c.alpha, c.M, c.lambda_u, c.rho_f, c.m, c.b, c.tau, c.lambda_u,
                            c.mu, c.rho};
      for (int g = 0; g < 10; g++)
      {
         Array<int> vd;
         gfs[g]->ParFESpace()->GetFaceNbrElementVDofs(j, vd);
         if (!SameBits(gfs[g]->FaceNbrData()(vd[0]), exp[g])) { badnbr++; }
      }
      nnbr++;
   }
   nnbr = LSum(nnbr, pm.GetComm());
   nnbrA = LSum(nnbrA, pm.GetComm());
   badnbr = LSum(badnbr, pm.GetComm());
   TEST_ASSERT(badnbr == 0, "U1-4 face-neighbour data of the sample and the undrained sample "
               "hold the neighbour's constants bit for bit (" << nnbr << " neighbours, "
               << nnbrA << " on side A, " << badnbr << " values differ)");

   const NitscheMaterialSample &u = s.UndrainedSample();
   const long du = BitDiff(u.LambdaGF(), s.LambdaUGF())
                   + BitDiff(u.MuGF(), s.DrainedSample().MuGF())
                   + BitDiff(u.RhoGF(), s.DrainedSample().RhoGF());
   TEST_ASSERT(LSum(du, pm.GetComm()) == 0, "U1-4 the undrained sample is (lambda_u, mu, rho) "
               "of the sample bit for bit");

   // The deck path: uniform tpv5-poro from a PoroConfig on F-B.
   std::unique_ptr<Mesh> mb(TPV102Box(4, 4, 2));
   ParMesh pb(MPI_COMM_WORLD, *mb);
   PoroConfig cfg;
   cfg.porosity = A.f.porosity;
   cfg.permeability_m2 = A.f.permeability;
   cfg.viscosity_pa_s = A.f.viscosity;
   cfg.tortuosity = A.f.tortuosity;
   cfg.rho_fluid = A.f.rho_fluid;
   cfg.rho_solid = A.f.rho_solid;
   cfg.bulk_fluid_pa = A.f.bulk_fluid;
   cfg.bulk_solid_pa = A.f.bulk_solid;
   PoroMaterialSample sb(pb, MaterialField::MakeConstant(A.lam, A.mu, cfg.BulkDensity()), cfg);
   long badb = 0;
   for (int e = 0; e < pb.GetNE(); e++) { badb += CheckElement(sb, e, cA); }
   TEST_ASSERT(LSum(badb, pb.GetComm()) == 0, "U1-4 the deck path (PoroConfig, uniform) on F-B "
               "gives DeriveBiotConstants bit for bit");
}

// ------------------------------------------------------------------- U1-5
static void TestU15()
{
   for (const Ref &r : Refs())
   {
      const BiotSpeeds s = ComputeBiotSpeeds(Derive(r));
      std::string which;
      const real_t worst = Worst("U1-5 " + r.name,
         {{"fast P", s.c_fast, r.c_fast}, {"slow P", s.c_slow, r.c_slow},
          {"S hf", s.c_s_hf, r.c_s_hf}, {"S lf", s.c_s_lf, r.c_s_lf},
          {"undrained P", s.c_p_undrained, r.c_p_undrained},
          {"drained P", s.c_p_drained, r.c_p_drained}}, which);
      TEST_ASSERT(worst <= 1e-13, "U1-5 " << r.name << ": the six speeds match the reference "
                  "to 1e-13 (worst " << which << " " << worst << ")");
   }
}

// ------------------------------------------------------------------- U1-6
/// PoroMaxDt on pm against (9.6) evaluated on the serial mesh, p = 1-3.
static void CheckMaxDt(const std::string &label, Mesh &serial, ParMesh &pm,
                       const MaterialField &mf, const PoroMaterialSample::FluidField &fluid,
                       const std::function<BiotConstants(const Vector &)> &constants_at)
{
   const real_t cfl = 0.5;
   PoroMaterialSample s(pm, mf, fluid);
   Vector x(3);
   for (int p = 1; p <= 3; p++)
   {
      const real_t dt = PoroMaxDt(s, pm, p, cfl);
      const real_t p2 = real_t(p) * real_t(p);
      real_t dt_ref = std::numeric_limits<real_t>::infinity();
      for (int e = 0; e < serial.GetNE(); e++)
      {
         serial.GetElementTransformation(e)->Transform(
            Geometries.GetCenter(serial.GetElementBaseGeometry(e)), x);
         const real_t c_e = ComputeBiotSpeeds(constants_at(x)).c_fast;
         const real_t h_e = NitscheElementLength(serial, e);
         dt_ref = std::min(dt_ref, cfl * h_e / (p2 * c_e));
      }
      char hex[64];
      std::snprintf(hex, sizeof(hex), "%a", dt);
      if (g_rank == 0)
      {
         std::cout.precision(17);
         std::cout << "  U1-6 " << label << " p = " << p << ": dt = " << dt << " (" << hex
                   << ")\n";
      }
      TEST_ASSERT(SameBits(dt, dt_ref), "U1-6 " << label << " p = " << p
                  << ": PoroMaxDt equals eq. (9.6) on the serial mesh bit for bit");
   }
}

static void TestU16(const std::string &mesh1000)
{
   const std::vector<Ref> R = Refs();
   const Ref &A = R[0], &B = R[1];
   const BiotConstants cA = Derive(A), cB = Derive(B);

   std::unique_ptr<Mesh> fb(TPV102Box(4, 4, 2));
   {
      ParMesh pm(MPI_COMM_WORLD, *fb);
      CheckMaxDt("F-B tpv5-poro", *fb, pm,
                 MaterialField::MakeConstant(A.lam, A.mu, BulkDensity(A.f)),
                 PoroMaterialSample::UniformField(A.f),
                 [&](const Vector &) { return cA; });
   }
   std::unique_ptr<Mesh> fc(TPV102Box(4, 4, 2));
   MarkSidesFC(*fc);
   {
      ParMesh pm(MPI_COMM_WORLD, *fc);
      Vector lv(2), mv(2), rv(2);
      lv(0) = A.lam; lv(1) = B.lam;
      mv(0) = A.mu;  mv(1) = B.mu;
      rv(0) = BulkDensity(A.f); rv(1) = BulkDensity(B.f);
      PWConstCoefficient lamc(lv), muc(mv), rhoc(rv);
      MaterialField mf;
      mf.mode = MaterialField::Mode::Coefficient;
      mf.lambda_coef = &lamc;
      mf.mu_coef = &muc;
      mf.rho_coef = &rhoc;
      CheckMaxDt("F-C", *fc, pm, mf,
                 [&](int, const Vector &x) { return x(0) < 0.0 ? A.f : B.f; },
                 [&](const Vector &x) { return x(0) < 0.0 ? cA : cB; });
   }
   if (mesh1000.empty())
   {
      if (g_rank == 0)
      {
         std::cout << "  U1-6 gmsh 1000 m TPV205 mesh: not run (no --mesh1000)\n";
      }
      return;
   }
   Mesh m1(mesh1000.c_str(), 1, 1);
   if (g_rank == 0)
   {
      std::cout << "  U1-6 " << mesh1000 << ": " << m1.GetNE() << " tets\n";
   }
   ParMesh pm(MPI_COMM_WORLD, m1);
   CheckMaxDt("gmsh 1000 m TPV205-equivalent", m1, pm,
              MaterialField::MakeConstant(B.lam, B.mu, BulkDensity(B.f)),
              PoroMaterialSample::UniformField(B.f),
              [&](const Vector &) { return cB; });
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   std::string mesh1000;
   for (int i = 1; i < argc; i++)
   {
      const std::string a = argv[i];
      if (a == "--mesh1000" && i + 1 < argc) { mesh1000 = argv[++i]; }
   }
   if (g_rank == 0)
   {
      std::cout << "=== Poroelastic P1 (plan rev 3) on " << Mpi::WorldSize() << " rank(s) ===\n";
   }
   TestU11();
   TestU12();
   TestU13();
   TestU14();
   TestU15();
   TestU16(mesh1000);
   int nf = num_failed;
   MPI_Allreduce(MPI_IN_PLACE, &nf, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
   if (g_rank == 0)
   {
      std::cout << "=== Summary: " << num_passed << "/" << num_tests << " passed, "
                << num_failed << " failed ===\n";
   }
   return nf == 0 ? 0 : 1;
}

// Nitsche dynamic rupture, phase P5: friction closure interface (N-03) and
// the rate-and-state adapter (N-04).
//
//   T-06(a): RateStateClosure::SolveSlipRate == SolveSlipRateVectorPsi with
//            gamma in the eta slot, bit for bit (both components, |V|).
//   T-06(b): StateRate == AgingLawPsi::Rate / SlipLawPsi::Rate, bit for bit.
//   T-06(c): SeedEquilibrium agrees with ComputeInitialPsi(a) to 1e-14 rel.
//   T-07:    V is anti-parallel to X for non-axis-aligned X.
//   T-11:    frozen dT_h/d delta == sigma_F; quasi-static energy form
//            delta^T M_f [T(delta) - T(0)] > 0 (Schur complement identity).
//
// Scope (Chunhui Zhao, 2026-10-14).

#include "mfem.hpp"
#include "../../friction/rate_state_closure.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../solver/nitsche_cfl.hpp"
#include "../../config/tpv102_params.hpp"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
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

static const real_t kL = 30.0e3;

// 24-tet TPV102 box with the fault on Y=0 (attr 3), free surface Z=0
// (attr 1), Dirichlet elsewhere (attr 5).
static Mesh *FaultBox(int nx, int ny, int nz)
{
   Mesh *mesh = new Mesh(Mesh::MakeCartesian3DWith24TetsPerHex(
      nx, ny, nz, 2*kL, 2*kL, kL));
   for (int v = 0; v < mesh->GetNV(); v++)
   {
      real_t *c = mesh->GetVertex(v);
      c[0] -= kL; c[1] -= kL; c[2] -= kL;
   }
   const real_t tol = kL / nz * 0.1;
   for (int b = 0; b < mesh->GetNBE(); b++)
   {
      Array<int> verts;
      mesh->GetBdrElementVertices(b, verts);
      real_t zc = 0.0;
      for (int i = 0; i < verts.Size(); i++) { zc += mesh->GetVertex(verts[i])[2]; }
      zc /= verts.Size();
      mesh->SetBdrAttribute(b, (std::abs(zc) < tol) ? 1 : 5);
   }
   const real_t ytol = 1e-6 * kL;
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

// Nitsche coupling gamma = 1/2 sigma_F dt at one fault face of the 24-tet box
// (p = 2): sigma_F from the IP integrator's GetPenalty, dt from the P3 CFL
// step (C_use = 0.11 at p = 2).
static real_t NitscheGamma(ParMesh &pm, real_t &sigmaF, real_t &dt)
{
   const int p = 2;
   const real_t lam = TPV102Params::lambda, mu = TPV102Params::mu;
   MaterialField mat = MaterialField::MakeConstant(lam, mu, TPV102Params::rho);
   NitscheMaterialSample sample(pm, mat);
   dt = ComputeMaxDt(sample, pm, p, 0.11);
   L2_FECollection fec(p, 3, BasisType::GaussLobatto);
   ParFiniteElementSpace fes(&pm, &fec);
   ConstantCoefficient lc(lam), mc(mu);
   DGElasticityIPCombinedIntegrator integ(lc, mc, 3, -1.0, 1.0);
   sigmaF = 0.0;
   for (int be = 0; be < pm.GetNBE(); be++)
   {
      if (pm.GetBdrAttribute(be) != 3) { continue; }
      const int f = pm.GetBdrElementFaceIndex(be);
      FaceElementTransformations *Tr = pm.GetInteriorFaceTransformations(f);
      if (!Tr) { continue; }
      const IntegrationPoint &ip = Geometries.GetCenter(Tr->GetGeometryType());
      Tr->SetAllIntPoints(&ip);
      Vector nor(3);
      CalcOrtho(Tr->Jacobian(), nor);
      const real_t nl = nor.Norml2();
      sigmaF = integ.GetPenalty(*fes.GetFE(Tr->Elem1No), *fes.GetFE(Tr->Elem2No),
                                Tr->Elem1->Weight(), Tr->Elem2->Weight(),
                                lam, mu, nl);
      break;
   }
   MPI_Allreduce(MPI_IN_PLACE, &sigmaF, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX,
                 pm.GetComm());
   return 0.5 * sigmaF * dt;
}

struct FrictionSet
{
   const char *name;
   DieterichRuinaFriction::Constants c;
   std::vector<real_t> a, psi, sigma;
};

// ---------------------------------------------------------- T-06(a,b), T-07
static void TestAdapter(const std::vector<real_t> &gammas)
{
   if (g_rank == 0) { std::cout << "\n=== T-06(a,b) adapter bit-identity, T-07 anti-parallel ===\n"; }
   const real_t lam = TPV102Params::lambda, mu = TPV102Params::mu;
   (void) lam;
   std::vector<FrictionSet> sets = {
      {"TPV102", {TPV102Params::V0, TPV102Params::f0, TPV102Params::b, TPV102Params::Dc},
       {0.008, 0.01, 0.016}, {0.3, 0.6, 0.77, 1.0}, {50.0e6, 120.0e6}},
      {"BP5-like", {1.0e-6, 0.6, 0.03, 0.13},
       {0.004, 0.01, 0.04}, {0.3, 0.6, 0.9}, {25.0e6}}};
   const real_t taus[] = {1.0e3, 1.0e6, 3.0e7, 7.5e7, 8.1e7, 1.0e8, 1.5e8};
   const real_t angs[] = {0.17, 0.523599, 1.97222, 4.38078};   // rad, non-axis
   long n_cases = 0, n_bad_a = 0, n_bad_b = 0, n_bad_7 = 0, n_inf = 0, n_inf_g = 0;
   real_t Vmin = 1e300, Vmax = 0.0, worst_cross = 0.0;
   for (const auto &S : sets)
   {
      DieterichRuinaFriction fr(S.c);
      RateStateClosure ag(fr, RateStateClosure::StateLaw::Aging);
      RateStateClosure sl(fr, RateStateClosure::StateLaw::Slip);
      AgingLawPsi ref_ag(S.c.b, S.c.V0, S.c.f0);
      SlipLawPsi ref_sl(S.c.b, S.c.V0, S.c.f0);
      for (real_t a : S.a)
         for (real_t psi : S.psi)
            for (real_t sig : S.sigma)
               for (real_t g : gammas)
                  for (real_t ta : taus)
                     for (real_t th : angs)
                     {
                        const real_t X[2] = {ta * std::cos(th), ta * std::sin(th)};
                        FaultPointState s; s.psi = psi; s.a = a; s.Dc = S.c.Dc;
                        real_t V[2], Vr[2];
                        const real_t Vabs = ag.SolveSlipRate(X, g, sig, s, V);
                        fr.SolveSlipRateVectorPsi(X, psi, sig, g, a, Vr);
                        const real_t Vabs_r = std::sqrt(Vr[0] * Vr[0] + Vr[1] * Vr[1]);
                        n_cases++;
                        if (!(V[0] == Vr[0] && V[1] == Vr[1] && Vabs == Vabs_r)) { n_bad_a++; }
                        // T-06(b)
                        if (!(ag.StateRate(Vabs, s) == ref_ag.Rate(Vabs, psi, S.c.Dc))
                            || !(sl.StateRate(Vabs, s) == ref_sl.Rate(Vabs, psi, S.c.Dc)))
                        {
                           n_bad_b++;
                        }
                        // T-07
                        const real_t dot = V[0] * X[0] + V[1] * X[1];
                        const real_t cross = std::abs(V[0] * X[1] - V[1] * X[0]);
                        const real_t rel = (Vabs > 0) ? cross / (Vabs * ta) : 0.0;
                        worst_cross = std::max(worst_cross, rel);
                        // gamma = 0 with |X| far above the strength overflows
                        // |V| to inf in the existing solver (the adapter
                        // reproduces it bitwise); T-07 is judged on finite |V|.
                        if (!std::isfinite(Vabs))
                        {
                           n_inf++;
                           if (g != 0.0) { n_inf_g++; }
                           continue;
                        }
                        if (Vabs > 0 && !(dot < 0.0 && rel < 1e-15)) { n_bad_7++; }
                        if (Vabs > 0) { Vmin = std::min(Vmin, Vabs); }
                        Vmax = std::max(Vmax, Vabs);
                     }
   }
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  cases = " << n_cases << ", |V| range ["
                << Vmin << ", " << Vmax << "] m/s\n"
                << "  T-06(a) mismatches = " << n_bad_a << ", T-06(b) mismatches = "
                << n_bad_b << ", T-07 failures = " << n_bad_7
                << ", worst |V x X|/(|V||X|) = " << worst_cross << "\n"
                << "  non-finite |V| (existing solver overflow): " << n_inf
                << " cases, of which with gamma > 0: " << n_inf_g << "\n";
   }
   TEST_ASSERT(n_bad_a == 0, "T-06(a) SolveSlipRate == SolveSlipRateVectorPsi(gamma as eta), bitwise");
   TEST_ASSERT(n_bad_b == 0, "T-06(b) StateRate == AgingLawPsi/SlipLawPsi::Rate, bitwise");
   TEST_ASSERT(n_bad_7 == 0, "T-07 V anti-parallel to X (V.X < 0, |V x X| < 1e-15 |V||X|), finite |V|");
   TEST_ASSERT(n_inf_g == 0, "non-finite |V| only at gamma = 0 (bounded by |X|/gamma otherwise)");
}

// ------------------------------------------------------------------ T-06(c)
static void TestSeed()
{
   if (g_rank == 0) { std::cout << "\n=== T-06(c) SeedEquilibrium vs ComputeInitialPsi ===\n"; }
   DieterichRuinaFriction fr({TPV102Params::V0, TPV102Params::f0, TPV102Params::b,
                              TPV102Params::Dc});
   RateStateClosure cl(fr, RateStateClosure::StateLaw::Aging);
   real_t worst = 0.0;
   int n = 0;
   for (real_t a = 0.008; a <= 0.0161; a += 0.0005)
   {
      FaultPointState s; s.a = a; s.Dc = TPV102Params::Dc;
      cl.SeedEquilibrium(TPV102Params::tau_ini, TPV102Params::V_ini,
                         TPV102Params::sigma_n, s);
      const real_t ref = ComputeInitialPsi(a);
      worst = std::max(worst, std::abs(s.psi - ref) / std::abs(ref));
      n++;
   }
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  a in [0.008, 0.016], " << n
                << " values: worst relative difference = " << worst << "\n";
   }
   TEST_ASSERT(worst < 1e-14, "T-06(c) SeedEquilibrium == ComputeInitialPsi to 1e-14 rel");
}



// delta^T M_f T over the fault (fault mass matrix from the fault face basis:
// M_face(k,l) = |J_face| sum_q w_q e_k(q) e_l(q)).  Shared fault faces are
// seen by both ranks and weighted 1/2.  Global sum.
static real_t FaultWork(ElasticityDomainOperator<ParMesh> &op, ParMesh &pm,
                        const Vector &d, const Vector &T)
{
   const FaceQuadrature *fq = op.GetFaceQuadrature();
   const int nbf = fq->NumBasisFunctions(), nq = fq->NumQuadPoints();
   const DenseMatrix &E = fq->BasisAtQuadPoints();
   const IntegrationRule &ir = fq->GetQuadRule();
   DenseMatrix Mref(nbf);
   Mref = 0.0;
   for (int q = 0; q < nq; q++)
      for (int k = 0; k < nbf; k++)
         for (int l = 0; l < nbf; l++) { Mref(k, l) += ir.IntPoint(q).weight * E(k, q) * E(l, q); }
   auto face_nl = [](FaceElementTransformations *Tr) {
      const IntegrationPoint &ip = Geometries.GetCenter(Tr->GetGeometryType());
      Tr->SetAllIntPoints(&ip);
      Vector nor(3);
      CalcOrtho(Tr->Jacobian(), nor);
      return nor.Norml2();
   };
   const Array<int> &fint = op.GetFaultInteriorFaces();
   const Array<int> &fsh = op.GetFaultSharedFaces();
   real_t w = 0.0;
   auto add_face = [&](int fi, real_t nl, real_t weight) {
      for (int c = 0; c < 2; c++)
         for (int k = 0; k < nbf; k++)
            for (int l = 0; l < nbf; l++)
            {
               w += weight * nl * Mref(k, l) * d(2 * (fi * nbf + k) + c)
                    * T(2 * (fi * nbf + l) + c);
            }
   };
   for (int i = 0; i < fint.Size(); i++)
   {
      add_face(i, face_nl(pm.GetInteriorFaceTransformations(fint[i])), 1.0);
   }
   for (int i = 0; i < fsh.Size(); i++)
   {
      add_face(fint.Size() + i, face_nl(pm.GetSharedFaceTransformations(fsh[i])), 0.5);
   }
   MPI_Allreduce(MPI_IN_PLACE, &w, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, pm.GetComm());
   return w;
}

// ---------------------------------------------------------------------- T-11
// Sign test on the 24-tet box with its fault (IP, p = 2, zero wall data):
// (1) frozen u_h = 0, uniform strike slip delta: dT_h/d delta per fault DOF
//     (ComputeTractionDiagnostics), compared with sigma_F from GetPenalty;
// (2) along the quasi-static solution K u = G delta (Solve), dT/d delta.
static void TestSign(ParMesh &pm)
{
   if (g_rank == 0) { std::cout << "\n=== T-11 sign test (24-tet 4x4x2 box, fault y=0, IP p=2) ===\n"; }
   MPI_Comm comm = pm.GetComm();
   const int p = 2;
   const real_t lam = TPV102Params::lambda, mu = TPV102Params::mu;
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = [](const Vector &, real_t, Vector &u) { u.SetSize(3); u = 0.0; };
   LinearElastic model(lam, mu);
   ElasticityDomainOperator<ParMesh> op(pm, p, model, 0.0, kL, 2.0 * kL, bc,
                                        DGMethod::IP, SolverType::MUMPS);
   const int nd = op.GetNumFaultDOFs();
   const real_t delta = 1.0e-3;            // [m], strike component
   Vector s0(2 * nd), s1(2 * nd);
   s0 = 0.0; s1 = 0.0;
   for (int k = 0; k < nd; k++) { s1(2 * k + 1) = delta; }

   auto traction = [&](const ParGridFunction &u, const Vector &slip, Vector &T) {
      Vector Ts, Tc, R;
      op.ComputeTractionDiagnostics(u, slip, T, Ts, Tc, R);
   };
   ParGridFunction u0(&op.GetFESpace());
   u0 = 0.0;
   Vector T0, T1;
   traction(u0, s0, T0);
   traction(u0, s1, T1);
   // (1) frozen
   real_t fmin = 1e300, fmax = -1e300, fsum = 0, dipmax = 0;
   for (int k = 0; k < nd; k++)
   {
      const real_t d = (T1(2 * k + 1) - T0(2 * k + 1)) / delta;
      fmin = std::min(fmin, d); fmax = std::max(fmax, d); fsum += d;
      dipmax = std::max(dipmax, std::abs(T1(2 * k) - T0(2 * k)) / delta);
   }
   real_t red[4] = {-fmin, fmax, fsum, dipmax};
   MPI_Allreduce(MPI_IN_PLACE, red, 2, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
   MPI_Allreduce(MPI_IN_PLACE, red + 2, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);
   MPI_Allreduce(MPI_IN_PLACE, red + 3, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
   fmin = -red[0]; fmax = red[1]; fsum = red[2]; dipmax = red[3];
   int ndg = nd;
   MPI_Allreduce(MPI_IN_PLACE, &ndg, 1, MPI_INT, MPI_SUM, comm);

   // sigma_F over the fault faces (GetPenalty, face centre).
   L2_FECollection fec(p, 3, BasisType::GaussLobatto);
   ParFiniteElementSpace sfes(&pm, &fec);
   ConstantCoefficient lc(lam), mc(mu);
   DGElasticityIPCombinedIntegrator integ(lc, mc, 3, -1.0, 1.0);
   real_t smin = 1e300, smax = -1e300;
   for (int f : op.GetFaultInteriorFaces())
   {
      FaceElementTransformations *Tr = pm.GetInteriorFaceTransformations(f);
      const IntegrationPoint &ip = Geometries.GetCenter(Tr->GetGeometryType());
      Tr->SetAllIntPoints(&ip);
      Vector nor(3);
      CalcOrtho(Tr->Jacobian(), nor);
      const real_t sF = integ.GetPenalty(*sfes.GetFE(Tr->Elem1No), *sfes.GetFE(Tr->Elem2No),
                                         Tr->Elem1->Weight(), Tr->Elem2->Weight(),
                                         lam, mu, nor.Norml2());
      smin = std::min(smin, sF); smax = std::max(smax, sF);
   }
   real_t sr[2] = {-smin, smax};
   MPI_Allreduce(MPI_IN_PLACE, sr, 2, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
   smin = -sr[0]; smax = sr[1];

   // (1) assertions: frozen dT_h/d delta equals sigma_F to round-off.
   const real_t sF = 0.5 * (smin + smax);
   const real_t frel = std::max(std::abs(fmin - sF), std::abs(fmax - sF)) / sF;
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  fault DOFs = " << ndg << ", delta = "
                << delta << " m (strike, all DOFs)\n"
                << "  frozen u_h: dT_strike/d delta in [" << fmin << ", " << fmax
                << "] Pa/m; |dT_dip/d delta| max " << dipmax << "\n"
                << "  sigma_F (GetPenalty, fault faces) in [" << smin << ", " << smax
                << "] Pa/m; max |dT/d delta - sigma_F| / sigma_F = " << frel << "\n";
   }
   TEST_ASSERT(smax - smin <= 1e-12 * sF && frel < 1e-12,
               "T-11 frozen dT_h/d delta == sigma_F (GetPenalty) to round-off");

   // (2) quasi-static energy form, for uniform and random slip fields:
   //   W = delta^T M_f [T(delta) - T(0)] > 0, and the identity
   //   W = delta^T M_f T_frozen(delta) - u^T K u   (K u = G delta),
   //   i.e. the Schur complement delta^T (P - G^T K^{-1} G) delta.
   ParGridFunction uz(&op.GetFESpace());
   op.Solve(0.0, s0, uz);
   Vector Tq0;
   traction(uz, s0, Tq0);
   Vector xyz;
   op.GetFaultDOFCoords3D(xyz);
   // Slip fields: uniform strike, a cosine-tapered strike patch (zero well
   // inside the clamped walls: |x| < 18 km, -24 < z < -3 km), and 5 random
   // fields (both components, fixed seeds).
   std::vector<std::string> names;
   std::vector<Vector> fields;
   names.push_back("uniform strike 1 mm"); fields.push_back(s1);
   {
      Vector d(2 * nd);
      d = 0.0;
      for (int k = 0; k < nd; k++)
      {
         const real_t x = xyz(3 * k), z = xyz(3 * k + 2);
         const real_t ax = 18.0e3, zc = -13.5e3, az = 10.5e3;
         if (std::abs(x) < ax && std::abs(z - zc) < az)
         {
            d(2 * k + 1) = delta * 0.25 * (1 + std::cos(M_PI * x / ax))
                           * (1 + std::cos(M_PI * (z - zc) / az));
         }
      }
      names.push_back("tapered strike patch"); fields.push_back(d);
   }
   for (int r = 1; r <= 5; r++)
   {
      Vector d(2 * nd);
      d.Randomize(1000 * r + g_rank);
      d -= 0.5;
      d *= 2.0 * delta;
      names.push_back("random " + std::to_string(r)); fields.push_back(d);
   }
   bool all_pos = true;
   std::vector<Vector> dTs;
   for (size_t c = 0; c < fields.size(); c++)
   {
      const Vector &d = fields[c];
      ParGridFunction ud(&op.GetFESpace());
      op.Solve(0.0, d, ud);
      Vector Tq, Tf, dT;
      traction(ud, d, Tq);
      traction(u0, d, Tf);
      dT = Tq;
      dT -= Tq0;
      dTs.push_back(dT);
      const real_t W = FaultWork(op, pm, d, dT);
      const real_t Wf = FaultWork(op, pm, d, Tf);
      Vector Ku;
      op.ApplyStiffness(ud, Ku);
      real_t uKu = ud * Ku;
      MPI_Allreduce(MPI_IN_PLACE, &uKu, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);
      const real_t idrel = std::abs(W - (Wf - uKu)) / std::abs(Wf);
      if (!(W > 0.0)) { all_pos = false; }
      if (g_rank == 0)
      {
         std::cout << std::setprecision(6) << "  " << names[c] << ": W = delta.M_f.dT = " << W
                   << "  [diag: delta.M_f.T_frozen = " << Wf << ", u.K.u = " << uKu
                   << ", |W - (that - u.K.u)| / that = " << idrel << "]\n";
      }
      if (c == 0)
      {
         // Pointwise signs on the static solution (diagnostic only).
         int nneg = 0, nneg_wall = 0;
         real_t qmin = 1e300, qmax = -1e300;
         for (int k = 0; k < nd; k++)
         {
            const real_t v = dT(2 * k + 1) / delta;
            qmin = std::min(qmin, v); qmax = std::max(qmax, v);
            if (v < 0)
            {
               nneg++;
               const real_t x = xyz(3 * k), z = xyz(3 * k + 2);
               if (std::abs(std::abs(x) - kL) < 7.6e3 || z < -kL + 7.6e3) { nneg_wall++; }
            }
         }
         int ci[2] = {nneg, nneg_wall};
         MPI_Allreduce(MPI_IN_PLACE, ci, 2, MPI_INT, MPI_SUM, comm);
         real_t mm[2] = {-qmin, qmax};
         MPI_Allreduce(MPI_IN_PLACE, mm, 2, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm);
         if (g_rank == 0)
         {
            std::cout << "    diagnostic (uniform slip, pointwise): dT_strike/d delta in ["
                      << -mm[0] << ", " << mm[1] << "] Pa/m; " << ci[0]
                      << " negative DOFs, " << ci[1]
                      << " of them within one hex of a clamped wall or the bottom\n";
         }
      }
   }
   TEST_ASSERT(all_pos, "T-11 energy form delta^T M_f [T(delta) - T(0)] > 0 (uniform, tapered patch, 5 random)");

   // Symmetry diagnostic: d1^T M_f dT(d2) vs d2^T M_f dT(d1) over random pairs.
   real_t worst_sym = 0.0;
   for (size_t i = 2; i < fields.size(); i++)
      for (size_t j = i + 1; j < fields.size(); j++)
      {
         const real_t a12 = FaultWork(op, pm, fields[i], dTs[j]);
         const real_t a21 = FaultWork(op, pm, fields[j], dTs[i]);
         const real_t aii = FaultWork(op, pm, fields[i], dTs[i]);
         const real_t ajj = FaultWork(op, pm, fields[j], dTs[j]);
         worst_sym = std::max(worst_sym, std::abs(a12 - a21) / std::sqrt(aii * ajj));
      }
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  symmetry (diagnostic, 10 random pairs): max |d1.M_f.dT(d2) - d2.M_f.dT(d1)| / sqrt(W1 W2) = "
                << worst_sym << "\n";
   }
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   if (g_rank == 0)
   {
      std::cout << "=== Nitsche P5 (friction closure) on " << Mpi::WorldSize()
                << " rank(s) ===\n";
   }
   std::unique_ptr<Mesh> m(FaultBox(4, 4, 2));
   ParMesh pm(MPI_COMM_WORLD, *m);
   real_t sigmaF, dt;
   const real_t gN = NitscheGamma(pm, sigmaF, dt);
   const real_t cs = std::sqrt(TPV102Params::mu / TPV102Params::rho);
   const real_t etaBP5 = TPV102Params::mu / (2.0 * cs);
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  gamma values: 0, BP5 eta = mu/(2 cs) = "
                << etaBP5 << ", Nitsche 1/2 sigma_F dt = " << gN << " (sigma_F = "
                << sigmaF << " Pa/m, dt = " << dt << " s, p = 2, 4x4x2 box), 1e10\n";
   }
   TestAdapter({0.0, etaBP5, gN, 1.0e10});
   TestSeed();
   TestSign(pm);
   if (g_rank == 0)
   {
      std::cout << "\n=== Summary: " << num_passed << "/" << num_tests
                << " passed, " << num_failed << " failed ===\n";
   }
   return (num_failed == 0) ? 0 : 1;
}

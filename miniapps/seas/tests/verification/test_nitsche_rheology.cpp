// Nitsche dynamic rupture with off-fault rheology: operator tests O1-O7 of
// the rheology plan rev 1 (section 10.2; phase R3).
//
// Mesh: the mirror-symmetric 60 x 60 x 30 km box of the P6 fixture (24 tets
// per hex, fault on the plane Y = 0 (attr 3), free surface Z = 0 (attr 1),
// absorbing walls (attr 5)), as test_nitsche_lsw_inline.cpp; O2 makes the
// walls X = +-L and the bottom Dirichlet (attr 6).  Material as TPV102;
// linear slip weakening with the TPV205 values and the smooth nucleation of
// test_nitsche_lsw_inline.cpp; prestress and sigma_0' sigma_XY = 70 MPa,
// sigma_YY = -120 MPa.
//
//   O1  Drucker-Prager with a cohesion that never yields (1e12 Pa): 200
//       steps identical to the run without the coupler, bit for bit (U,
//       Udot, d, V, slip path), p = 1 and 2.
//   O2  patch test (plan Proposition A.1): a uniform dsig with six distinct
//       components gives P dsig = -int_{free, absorbing} (dsig n) . w and
//       Q dsig = dsig nhat at every fault DOF, 1e-12 relative, p = 1 and 2,
//       with fault, interior, shared and Dirichlet faces.
//   O3  S exact on displacements of degree p, p = 1, 2, 3: 1e-12 relative.
//   O4  the linear part: for a continuous degree-p displacement v and
//       dsig = C_0 : eps(v) at the nodes, P dsig = -K v (the total-stress
//       weak form of a fully relaxed field gives no force), 1e-12 relative,
//       p = 1 and 2 (no Dirichlet faces: there K also holds the symmetry and
//       penalty terms on the boundary data); with the same dsig, Q dsig =
//       (C_0 : eps(v)) nhat at every owned fault DOF (at p = 2 a linear field,
//       shared fault faces included on 2 ranks).
//   O5  rank invariance with yielding (cohesion 91 MPa, bulk friction 0.6,
//       Tv 0.03 s: the initial sqrt(J2) = 98.49 MPa is 0.12 MPa below the yield
//       stress), 2.5 s: fault data matched by position and face centroid
//       against a one-rank run on MPI_COMM_SELF, 1e-12 relative; same yield
//       and active counts at the end.  On 2 and 4 ranks the partition cuts
//       the fault (every fault face shared) and the yielding zone.
//   O6  a block in simple shear: every face Dirichlet with u = (gdot t y,
//       0, 0) and the matching initial velocity, sigma_0' = -50 MPa I: the
//       stress path at every node matches the pointwise model driven by the
//       exact strain history, 1e-10 relative, Tv = 0 and 0.03 s.
//   O7  refusals: --lts rate2, energy_balance, output_dt_s without parallel
//       HDF5, [stress] / [[stress_region]] / P_p next to [initial_stress]:
//       the message names the key.
//
// Usage: mpirun -np 1|2|4 ./seas_test_nitsche_rheology
//
// Scope (Chunhui Zhao, 2026-10-20).

#include "mfem.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../solver/nitsche_cfl.hpp"
#include "../../friction/slip_weakening_closure.hpp"
#include "../../config/tpv102_params.hpp"
#include "../../rheology/rheology_coupler.hpp"
#include "../../rheology/rheology_config.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <map>
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
static const real_t kCuse[5] = {0.0, 0.057, 0.11, 0.14, 0.166};   // 0.5 C_eig (P3)
static const real_t kMuS = 0.677, kMuD = 0.525, kDc = 0.40;
static const real_t kTau0 = 70.0e6, kSigma0 = 120.0e6;
static const real_t kNucR = 9.0e3, kNucDtau = 15.0e6;

static int ISum(int v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_INT, MPI_SUM, c); return v; }
static real_t RMax(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, c); return v; }

// ---------------------------------------------------------------- fixture

/// The P6 box; dirichlet: walls X = +-L and the bottom get attr 6.
static Mesh *Box(int nx, int ny, int nz, bool dirichlet)
{
   Mesh *mesh = new Mesh(Mesh::MakeCartesian3DWith24TetsPerHex(
      nx, ny, nz, 2 * kL, 2 * kL, kL));
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
      real_t xc[3] = {0.0, 0.0, 0.0};
      for (int i = 0; i < verts.Size(); i++)
      {
         for (int d = 0; d < 3; d++) { xc[d] += mesh->GetVertex(verts[i])[d] / verts.Size(); }
      }
      int attr = (std::abs(xc[2]) < tol) ? 1 : 5;
      if (dirichlet && (std::abs(std::abs(xc[0]) - kL) < tol || std::abs(xc[2] + kL) < tol))
      {
         attr = 6;
      }
      mesh->SetBdrAttribute(b, attr);
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

static BoundaryConfig BC(bool dirichlet)
{
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.absorbing_attrs = {5};
   if (dirichlet)
   {
      bc.dirichlet_attrs = {6};
      bc.default_dirichlet_func = [](const Vector &, real_t, Vector &u) { u.SetSize(3); u = 0.0; };
   }
   return bc;
}

static MaterialField Material()
{
   return MaterialField::MakeConstant(TPV102Params::lambda, TPV102Params::mu,
                                      TPV102Params::rho);
}

static void Prestress(const real_t *, real_t s[9])
{
   for (int i = 0; i < 9; i++) { s[i] = 0.0; }
   s[1] = s[3] = kTau0;
   s[4] = -kSigma0;
}

static real_t NucF(real_t r)
{
   if (r >= kNucR) { return 0.0; }
   const real_t r2 = r * r, R2 = kNucR * kNucR;
   return std::exp(r2 / (r2 - R2));
}

static NitscheFaultSetup LSWSetup(const IFaultClosure &cl)
{
   NitscheFaultSetup fs;
   fs.closure = &cl;
   fs.prestress = Prestress;
   fs.nucleation_spatial = [](const real_t x[3], real_t s[9]) {
      for (int i = 0; i < 9; i++) { s[i] = 0.0; }
      const real_t v = kNucDtau * NucF(std::hypot(x[0] - TPV102Params::hypo_along_strike,
                                                 x[2] + TPV102Params::hypo_down_dip));
      s[1] = s[3] = v;
   };
   fs.nucleation_temporal = [](real_t t) { return NucleationTemporal(t); };
   fs.parameters = [](const real_t[3], FaultPointState &st) {
      st.mu_s = kMuS;
      st.mu_d = kMuD;
      st.Dc = kDc;
      st.cohesion = 0.0;
   };
   fs.sigma_min = 0.0;
   return fs;
}

static void InitialState(NitscheWaveOperator &op)
{
   const int no = op.NumOwnedFaultDOFs();
   Vector d(3 * no), V(3 * no), psi(no);
   d = 0.0;
   V = 0.0;
   for (int j = 0; j < no; j++)
   {
      FaultPointState st = op.FaultParams(j);
      op.Closure().SeedEquilibrium(kTau0, 0.0, kSigma0, st);
      psi(j) = st.psi;
   }
   op.SetFaultState(d, V, psi);
}

static std::unique_ptr<RheologyStack> DPStack(real_t c, real_t nu, real_t tv)
{
   RheologyConfig cfg;
   cfg.drucker_prager = true;
   cfg.dp.cohesion = c;
   cfg.dp.bulk_friction = nu;
   cfg.dp.relaxation_time = tv;
   return BuildRheologyStack(cfg);
}

static std::unique_ptr<RheologyCoupler> MakeCoupler(NitscheWaveOperator &op,
                                                    std::unique_ptr<RheologyStack> st,
                                                    const RheologyField::StressFn &s0)
{
   return std::make_unique<RheologyCoupler>(op.Elasticity().GetFESpace(), op.Sample(),
                                            op.Boundary().dirichlet_attrs, std::move(st), s0);
}

// ------------------------------------------------------------- one LSW run

struct LSWRun
{
   Vector U, Ud, d, V, psi, sig, X;
   long long yielding = 0, active = 0;
   int shared_fault = 0, active_shared = 0;   ///< shared fault faces; shared face sides on active elements
   real_t eta_max = 0.0;
   long steps = 0;
   bool finite = true;
   std::vector<std::array<long long, 6>> keys;   ///< position, face centroid (mm) per owned DOF
};

/// Face centroid per owned fault DOF (3 no): with the position it keys a DG
/// fault DOF across rank counts (as test_nitsche_lsw_inline.cpp).
static void OwnedFaceCentroids(NitscheWaveOperator &op, Vector &c_own)
{
   ElasticityDomainOperator<ParMesh> &E = op.Elasticity();
   const int nd = E.GetNumFaultDOFs(), nbf = E.GetNbfPerFace();
   Vector c_loc(3 * nd);
   for (int k = 0; k < nd; k++)
   {
      long long key[3];
      real_t xv[3][3];
      op.Orientation().FaceVertexData(k / nbf, key, xv);
      for (int d = 0; d < 3; d++) { c_loc(3 * k + d) = (xv[0][d] + xv[1][d] + xv[2][d]) / 3.0; }
   }
   E.RestrictToOwnedFault(c_loc, c_own, 3);
}

/// mode 0: no coupler; 1: never-yielding Drucker-Prager; 2: yielding.
/// nsteps > 0 steps, else to tfinal.
static LSWRun RunLSW(MPI_Comm comm, int p, long nsteps, int mode, real_t tfinal = 0.0)
{
   LSWRun R;
   std::unique_ptr<Mesh> m(Box(4, 4, 2, false));
   // on 2 and 4 ranks the cut runs through the fault (y < 0 | y > 0, and
   // x < 0 | x > 0 on 4 ranks), so every fault face is a shared face
   int nranks;
   MPI_Comm_size(comm, &nranks);
   std::vector<int> part;
   if (nranks == 2 || nranks == 4)
   {
      part.resize(m->GetNE());
      Vector c(3);
      for (int e = 0; e < m->GetNE(); e++)
      {
         m->GetElementCenter(e, c);
         part[e] = (c(1) > 0.0 ? 1 : 0) + ((nranks == 4 && c(0) > 0.0) ? 2 : 0);
      }
   }
   ParMesh pm(comm, *m, part.empty() ? nullptr : part.data());
   m.reset();
   SlipWeakeningClosure cl;
   NitscheFaultSetup fs = LSWSetup(cl);
   NitscheWaveOperator op(pm, p, Material(), BC(false), {}, &fs);
   R.shared_fault = op.NumSharedFaultFacesGlobal();
   const real_t dt = ComputeMaxDt(op.Sample(), pm, p, kCuse[p]);
   op.SetTimeStep(dt);
   InitialState(op);
   std::unique_ptr<RheologyCoupler> rc;
   if (mode == 1) { rc = MakeCoupler(op, DPStack(1e12, 0.6, 0.03), Prestress); }
   if (mode == 2) { rc = MakeCoupler(op, DPStack(91e6, 0.6, 0.03), Prestress); }
   if (rc) { op.SetRheology(rc.get()); }
   const int n = op.Height();
   Vector U(n), Ud(n);
   U = 0.0;
   Ud = 0.0;
   real_t t = 0.0;
   op.SetTime(t);
   CentralDifferenceSolver ode;
   ode.Init(op);
   if (nsteps <= 0) { nsteps = std::lround(std::ceil(tfinal / dt - 1e-9)); }
   R.steps = nsteps;
   for (long s = 0; s < nsteps; s++)
   {
      real_t h = dt;
      ode.Step(U, Ud, t, h);
      if (!std::isfinite(op.LastStep().max_V)) { R.finite = false; break; }
   }
   R.finite = ISum(R.finite ? 0 : 1, comm) == 0;
   R.U = U; R.Ud = Ud;
   R.d = op.Slip(); R.V = op.SlipRate(); R.psi = op.State();
   R.sig = op.LastSigmaEff(); R.X = op.LastShear();
   if (rc)
   {
      R.yielding = rc->YieldingNodes();
      R.active = rc->ActiveElements();
      R.eta_max = rc->MaxDiagnostic("eta");
      R.active_shared = ISum(rc->Operators().NumActiveFaces(rc->Field(), 1), comm);
   }
   const Vector &Xp = op.FaultPositions();
   Vector cf;
   OwnedFaceCentroids(op, cf);
   for (int j = 0; j < op.NumOwnedFaultDOFs(); j++)
   {
      R.keys.push_back({std::llround(Xp(3 * j) * 1e3), std::llround(Xp(3 * j + 1) * 1e3),
                        std::llround(Xp(3 * j + 2) * 1e3), std::llround(cf(3 * j) * 1e3),
                        std::llround(cf(3 * j + 1) * 1e3), std::llround(cf(3 * j + 2) * 1e3)});
   }
   return R;
}

static bool BitEqual(const Vector &a, const Vector &b)
{
   if (a.Size() != b.Size()) { return false; }
   for (int i = 0; i < a.Size(); i++)
   {
      if (std::memcmp(&a(i), &b(i), sizeof(real_t)) != 0) { return false; }
   }
   return true;
}

// --------------------------------------------------------------------- O1
static void TestO1()
{
   if (g_rank == 0) { std::cout << "\n=== O1 never-yielding Drucker-Prager is the elastic run, bit for bit ===\n"; }
   for (int p = 1; p <= 2; p++)
   {
      LSWRun a = RunLSW(MPI_COMM_WORLD, p, 200, 0);
      LSWRun b = RunLSW(MPI_COMM_WORLD, p, 200, 1);
      int same = (BitEqual(a.U, b.U) && BitEqual(a.Ud, b.Ud) && BitEqual(a.d, b.d)
                  && BitEqual(a.V, b.V) && BitEqual(a.psi, b.psi)) ? 0 : 1;
      same = ISum(same, MPI_COMM_WORLD);
      const real_t vmax = RMax(a.V.Normlinf(), MPI_COMM_WORLD);
      if (g_rank == 0)
      {
         std::cout << "  p = " << p << ": 200 steps, max|V| " << vmax << " m/s, yielding nodes "
                   << b.yielding << ", active elements " << b.active << "\n";
      }
      TEST_ASSERT(a.finite && b.finite && vmax > 0.0, "O1 p = " << p << " runs slip");
      TEST_ASSERT(same == 0 && b.active == 0, "O1 p = " << p
                  << " U, Udot, d, V, slip path identical bit for bit");
   }
}

// --------------------------------------------------------------------- O2
static void TestO2()
{
   if (g_rank == 0) { std::cout << "\n=== O2 patch test: uniform dsig ===\n"; }
   const real_t D[6] = {3.0e6, -2.0e6, 1.5e6, 0.7e6, -1.1e6, 0.4e6};   // xx yy zz xy yz xz
   for (int p = 1; p <= 2; p++)
   {
      std::unique_ptr<Mesh> m(Box(2, 2, 2, true));
      ParMesh pm(MPI_COMM_WORLD, *m);
      m.reset();
      SlipWeakeningClosure cl;
      NitscheFaultSetup fs = LSWSetup(cl);
      NitscheWaveOperator op(pm, p, Material(), BC(true), {}, &fs);
      auto rc = MakeCoupler(op, DPStack(1e12, 0.6, 0.0), Prestress);
      op.SetRheology(rc.get());
      RheologyField &F = rc->Field();
      for (int e = 0; e < F.NumElements(); e++)
      {
         for (int i = 0; i < F.NodesPerElement(); i++)
         {
            for (int c = 0; c < 6; c++) { F.Dsig(e, i)[c] = D[c]; }
         }
      }
      rc->Refresh();
      ParFiniteElementSpace &ufes = op.Elasticity().GetFESpace();
      Vector Fp(ufes.GetVSize()), Fr(ufes.GetVSize());
      Fp = 0.0;
      rc->AddLoad(Fp);
      // reference: -int_{free, absorbing} (dsig n) . w
      Fr = 0.0;
      real_t S[9];
      sym6::ToRowMajor(D, S);
      Array<int> vd;
      Vector nor(3), sh;
      int nN = 0;
      for (int be = 0; be < pm.GetNBE(); be++)
      {
         const int a = pm.GetBdrAttribute(be);
         if (a != 1 && a != 5) { continue; }
         FaceElementTransformations *Tr = pm.GetBdrFaceTransformations(be);
         if (!Tr) { continue; }
         nN++;
         const FiniteElement &fe = *ufes.GetFE(Tr->Elem1No);
         const int nd = fe.GetDof();
         sh.SetSize(nd);
         ufes.GetElementVDofs(Tr->Elem1No, vd);
         const IntegrationRule &ir = IntRules.Get(Tr->GetGeometryType(), 2 * p + 1);
         for (int q = 0; q < ir.GetNPoints(); q++)
         {
            const IntegrationPoint &ip = ir.IntPoint(q);
            Tr->SetAllIntPoints(&ip);
            CalcOrtho(Tr->Jacobian(), nor);
            fe.CalcShape(Tr->GetElement1IntPoint(), sh);
            for (int c = 0; c < 3; c++)
            {
               const real_t t = ip.weight * (S[3 * c] * nor(0) + S[3 * c + 1] * nor(1)
                                             + S[3 * c + 2] * nor(2));
               for (int j = 0; j < nd; j++) { Fr(vd[c * nd + j]) -= sh(j) * t; }
            }
         }
      }
      real_t scale = 0.0, err = 0.0;
      for (int i = 0; i < Fr.Size(); i++)
      {
         scale = std::max(scale, std::abs(Fr(i)));
         err = std::max(err, std::abs(Fp(i) - Fr(i)));
      }
      scale = RMax(scale, MPI_COMM_WORLD);
      err = RMax(err, MPI_COMM_WORLD);
      const int nN_all = ISum(nN, MPI_COMM_WORLD);
      const int nD = ISum(rc->Operators().NumFaces(2), MPI_COMM_WORLD);
      const int nS = ISum(rc->Operators().NumFaces(1), MPI_COMM_WORLD);
      // Q: dsig nhat at each owned fault DOF
      Vector T;
      op.RheologyTraction(T);
      real_t qerr = 0.0, qscale = 0.0;
      for (int j = 0; j < op.NumOwnedFaultDOFs(); j++)
      {
         const real_t *n = &op.FaultFrames()(9 * j);
         for (int c = 0; c < 3; c++)
         {
            const real_t ex = S[3 * c] * n[0] + S[3 * c + 1] * n[1] + S[3 * c + 2] * n[2];
            qscale = std::max(qscale, std::abs(ex));
            qerr = std::max(qerr, std::abs(T(3 * j + c) - ex));
         }
      }
      qscale = RMax(qscale, MPI_COMM_WORLD);
      qerr = RMax(qerr, MPI_COMM_WORLD);
      if (g_rank == 0)
      {
         std::cout << std::setprecision(4) << "  p = " << p << ": " << nN_all
                   << " free/absorbing faces, " << nD << " Dirichlet faces, " << nS
                   << " shared face sides; |P dsig - ref| / max|ref| = " << err / scale
                   << ", |Q dsig - dsig nhat| / |dsig nhat| = " << qerr / qscale << "\n";
      }
      TEST_ASSERT(nD > 0 && nN_all > 0, "O2 p = " << p << " has Dirichlet and free/absorbing faces");
      TEST_ASSERT(err <= 1e-12 * scale, "O2 p = " << p << " P dsig = -int_N (dsig n) . w");
      TEST_ASSERT(qerr <= 1e-12 * qscale, "O2 p = " << p << " Q dsig = dsig nhat");
   }
}

// --------------------------------------------------------------- O3 and O4
/// A displacement of degree p with all coefficients nonzero (scaled by 1/L
/// per power so the strain is O(1e-4)).
static void PolyU(int p, const real_t x[3], real_t u[3], real_t H[9])
{
   const real_t a = 1e-4, X = x[0] / kL, Y = x[1] / kL, Z = x[2] / kL;
   // u_c = a L (c1 X + c2 Y + c3 Z + [p>=2] quad + [p>=3] cubic)
   const real_t k[3][3] = {{0.3, -0.7, 0.2}, {0.5, 0.1, -0.4}, {-0.2, 0.6, 0.9}};
   for (int c = 0; c < 3; c++)
   {
      u[c] = a * kL * (k[c][0] * X + k[c][1] * Y + k[c][2] * Z);
      H[3 * c] = a * k[c][0];
      H[3 * c + 1] = a * k[c][1];
      H[3 * c + 2] = a * k[c][2];
      if (p >= 2)
      {
         const real_t q = 0.3 + 0.1 * c;
         u[c] += a * kL * q * (X * Y + 0.5 * Z * Z - 0.25 * X * X);
         H[3 * c] += a * q * (Y - 0.5 * X);
         H[3 * c + 1] += a * q * X;
         H[3 * c + 2] += a * q * Z;
      }
      if (p >= 3)
      {
         const real_t r = 0.2 - 0.05 * c;
         u[c] += a * kL * r * (X * Y * Z + Y * Y * Y);
         H[3 * c] += a * r * (Y * Z);
         H[3 * c + 1] += a * r * (X * Z + 3.0 * Y * Y);
         H[3 * c + 2] += a * r * (X * Y);
      }
   }
}

static void TestO3O4()
{
   if (g_rank == 0) { std::cout << "\n=== O3 strain sampling, O4 the linear part ===\n"; }
   for (int p = 1; p <= 3; p++)
   {
      std::unique_ptr<Mesh> m(Box(2, 2, 2, false));
      ParMesh pm(MPI_COMM_WORLD, *m);
      m.reset();
      SlipWeakeningClosure cl;
      NitscheFaultSetup fs = LSWSetup(cl);
      NitscheWaveOperator op(pm, p, Material(), BC(false), {}, &fs);
      auto rc = MakeCoupler(op, DPStack(1e12, 0.6, 0.0), Prestress);
      ParFiniteElementSpace &ufes = op.Elasticity().GetFESpace();
      ParGridFunction v(&ufes);
      VectorFunctionCoefficient vc(3, [p](const Vector &x, Vector &u) {
         real_t uu[3], H[9];
         PolyU(p, x.GetData(), uu, H);
         u.SetSize(3);
         for (int c = 0; c < 3; c++) { u(c) = uu[c]; }
      });
      v.ProjectCoefficient(vc);
      RheologyField &F = rc->Field();
      const int nr = F.NodesPerElement();
      std::vector<real_t> eps(6 * nr);
      real_t err = 0.0, scale = 0.0;
      for (int e = 0; e < F.NumElements(); e++)
      {
         rc->Operators().Strain(e, v, eps.data());
         for (int i = 0; i < nr; i++)
         {
            real_t uu[3], H[9], ex[6];
            PolyU(p, F.Context(e, i).x, uu, H);
            sym6::FromRowMajor(H, ex);
            for (int c = 0; c < 6; c++)
            {
               scale = std::max(scale, std::abs(ex[c]));
               err = std::max(err, std::abs(eps[6 * i + c] - ex[c]));
            }
            // O4: dsig = C_0 : eps(v) at the node
            sym6::ApplyIsotropic(F.Context(e, i).lambda0, F.Context(e, i).mu0, ex, F.Dsig(e, i));
         }
      }
      err = RMax(err, MPI_COMM_WORLD);
      scale = RMax(scale, MPI_COMM_WORLD);
      if (g_rank == 0)
      {
         std::cout << std::setprecision(4) << "  p = " << p << ": " << nr
                   << " nodes per element, max |S v - eps(v)| / max|eps| = " << err / scale << "\n";
      }
      TEST_ASSERT(err <= 1e-12 * scale, "O3 p = " << p << " S exact on degree-p displacements");
      if (p > 2) { continue; }
      rc->Refresh();
      Vector Fp(ufes.GetVSize()), Kv(ufes.GetVSize());
      Fp = 0.0;
      rc->AddLoad(Fp);
      op.Elasticity().ApplyStiffness(v, Kv);
      real_t e4 = 0.0, s4 = 0.0;
      for (int i = 0; i < Kv.Size(); i++)
      {
         s4 = std::max(s4, std::abs(Kv(i)));
         e4 = std::max(e4, std::abs(Fp(i) + Kv(i)));
      }
      e4 = RMax(e4, MPI_COMM_WORLD);
      s4 = RMax(s4, MPI_COMM_WORLD);
      if (g_rank == 0)
      {
         std::cout << "  p = " << p << ": max |P dsig + K v| / max|K v| = " << e4 / s4 << "\n";
      }
      TEST_ASSERT(e4 <= 1e-12 * s4, "O4 p = " << p << " P (C_0 : eps(v)) = -K v");
      // Q with the same dsig (degree p - 1, the same polynomial on both
      // sides of every face): Q dsig = (C_0 : eps(v)) nhat at each owned
      // fault DOF, shared fault faces included
      op.SetRheology(rc.get());
      Vector T;
      op.RheologyTraction(T);
      const Vector &X = op.FaultPositions();
      real_t qe = 0.0, qs = 0.0;
      for (int j = 0; j < op.NumOwnedFaultDOFs(); j++)
      {
         real_t uu[3], H[9], ex[6], s6[6], S9[9];
         PolyU(p, &X(3 * j), uu, H);
         sym6::FromRowMajor(H, ex);
         sym6::ApplyIsotropic(TPV102Params::lambda, TPV102Params::mu, ex, s6);
         sym6::ToRowMajor(s6, S9);
         const real_t *n = &op.FaultFrames()(9 * j);
         for (int c = 0; c < 3; c++)
         {
            const real_t t = S9[3 * c] * n[0] + S9[3 * c + 1] * n[1] + S9[3 * c + 2] * n[2];
            qs = std::max(qs, std::abs(t));
            qe = std::max(qe, std::abs(T(3 * j + c) - t));
         }
      }
      qe = RMax(qe, MPI_COMM_WORLD);
      qs = RMax(qs, MPI_COMM_WORLD);
      const int nsf = ISum(rc->Operators().NumFaces(1), MPI_COMM_WORLD);
      if (g_rank == 0)
      {
         std::cout << "  p = " << p << ": max |Q dsig - (C_0 : eps(v)) nhat| / max = " << qe / qs
                   << " (" << nsf << " shared face sides)\n";
      }
      TEST_ASSERT(qe <= 1e-12 * qs, "O4 p = " << p << " Q (C_0 : eps(v)) = (C_0 : eps(v)) nhat");
   }
}

// --------------------------------------------------------------------- O5
static void TestO5()
{
   if (g_rank == 0) { std::cout << "\n=== O5 rank invariance with yielding ===\n"; }
   int nranks;
   MPI_Comm_size(MPI_COMM_WORLD, &nranks);
   for (int p = 1; p <= 2; p++)
   {
      LSWRun a = RunLSW(MPI_COMM_WORLD, p, 0, 2, 2.5);
      const long nsteps = a.steps;
      // keyed rows: position, face centroid (mm) -> d, V, psi, sigma
      std::vector<double> buf;
      for (size_t j = 0; j < a.keys.size(); j++)
      {
         for (int c = 0; c < 6; c++) { buf.push_back(static_cast<double>(a.keys[j][c])); }
         for (int c = 0; c < 3; c++) { buf.push_back(a.d(3 * j + c)); }
         for (int c = 0; c < 3; c++) { buf.push_back(a.V(3 * j + c)); }
         buf.push_back(a.psi(j));
         buf.push_back(a.sig(j));
      }
      int nloc = static_cast<int>(buf.size());
      std::vector<int> counts(nranks), displs(nranks, 0);
      MPI_Gather(&nloc, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
      std::vector<double> all;
      if (g_rank == 0)
      {
         for (int r = 1; r < nranks; r++) { displs[r] = displs[r - 1] + counts[r - 1]; }
         all.resize(displs[nranks - 1] + counts[nranks - 1]);
      }
      MPI_Gatherv(buf.data(), nloc, MPI_DOUBLE, all.data(), counts.data(), displs.data(),
                  MPI_DOUBLE, 0, MPI_COMM_WORLD);
      int ok = 1;
      if (g_rank == 0)
      {
         LSWRun b = RunLSW(MPI_COMM_SELF, p, nsteps, 2);
         // a DG fault DOF repeats per face: the face centroid tells them apart
         std::map<std::array<long long, 6>, std::vector<std::array<double, 8>>> A, B;
         for (size_t k = 0; k + 14 <= all.size(); k += 14)
         {
            std::array<long long, 6> key;
            for (int c = 0; c < 6; c++) { key[c] = std::llround(all[k + c]); }
            std::array<double, 8> v;
            for (int c = 0; c < 8; c++) { v[c] = all[k + 6 + c]; }
            A[key].push_back(v);
         }
         for (size_t j = 0; j < b.keys.size(); j++)
         {
            std::array<double, 8> v;
            for (int c = 0; c < 3; c++) { v[c] = b.d(3 * j + c); v[3 + c] = b.V(3 * j + c); }
            v[6] = b.psi(j);
            v[7] = b.sig(j);
            B[b.keys[j]].push_back(v);
         }
         real_t sc[8] = {0, 0, 0, 0, 0, 0, 0, 0};
         for (auto &kv : B)
         {
            for (auto &v : kv.second) { for (int c = 0; c < 8; c++) { sc[c] = std::max(sc[c], std::abs(v[c])); } }
         }
         const real_t sd = std::max({sc[0], sc[1], sc[2]}), sv = std::max({sc[3], sc[4], sc[5]});
         const real_t scale[8] = {sd, sd, sd, sv, sv, sv, sc[6], sc[7]};
         real_t worst = 0.0;
         bool keys_ok = A.size() == B.size();
         for (auto &kv : B)
         {
            auto it = A.find(kv.first);
            if (it == A.end() || it->second.size() != kv.second.size()) { keys_ok = false; continue; }
            auto x = kv.second, y = it->second;
            std::sort(x.begin(), x.end());
            std::sort(y.begin(), y.end());
            for (size_t r = 0; r < x.size(); r++)
            {
               for (int c = 0; c < 8; c++)
               {
                  if (scale[c] > 0.0) { worst = std::max(worst, std::abs(x[r][c] - y[r][c]) / scale[c]); }
               }
            }
         }
         std::cout << std::setprecision(4) << "  p = " << p << ", " << nsteps << " steps to 2.5 s on "
                   << nranks << " rank(s): yielding nodes " << a.yielding << " (1 rank " << b.yielding
                   << "), active elements " << a.active << " (" << b.active << "), max eta "
                   << a.eta_max << " (" << b.eta_max << "), max|V| " << sv
                   << " m/s; " << a.shared_fault << " shared fault faces, " << a.active_shared
                   << " shared face sides on active elements; worst relative difference " << worst << "\n";
         ok = keys_ok && a.active > 0 && a.active == b.active && a.yielding == b.yielding
              && worst <= 1e-12 && std::abs(a.eta_max - b.eta_max) <= 1e-12 * b.eta_max;
         if (nranks > 1 && (a.shared_fault == 0 || a.active_shared == 0))
         {
            std::cout << "  the run must cut the fault and the yielding zone\n";
            ok = 0;
         }
         if (!keys_ok) { std::cout << "  fault DOF keys differ between the runs\n"; }
      }
      MPI_Bcast(&ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
      TEST_ASSERT(a.finite, "O5 p = " << p << " finite");
      TEST_ASSERT(ok, "O5 p = " << p << " yielding run equals the one-rank run (1e-12)");
   }
}

// --------------------------------------------------------------------- O6
static void TestO6()
{
   if (g_rank == 0) { std::cout << "\n=== O6 simple shear block ===\n"; }
   const real_t gdot = 2.0e-4, p0 = 50.0e6;   // 1/s, Pa
   for (real_t tv : {0.0, 0.03})
   {
      const int p = 2;
      // 6 km box, every face Dirichlet
      Mesh *mm = new Mesh(Mesh::MakeCartesian3DWith24TetsPerHex(2, 2, 2, 6.0e3, 6.0e3, 6.0e3));
      for (int v = 0; v < mm->GetNV(); v++)
      {
         real_t *c = mm->GetVertex(v);
         c[0] -= 3.0e3; c[1] -= 3.0e3; c[2] -= 6.0e3;
      }
      for (int b = 0; b < mm->GetNBE(); b++) { mm->SetBdrAttribute(b, 6); }
      mm->SetAttributes();
      std::unique_ptr<Mesh> m(mm);
      ParMesh pm(MPI_COMM_WORLD, *m);
      m.reset();
      BoundaryConfig bc;
      bc.dirichlet_attrs = {6};
      bc.default_dirichlet_func = [gdot](const Vector &x, real_t t, Vector &u)
      {
         u.SetSize(3);
         u = 0.0;
         u(0) = gdot * t * x(1);
      };
      NitscheWaveOperator op(pm, p, Material(), bc, {}, nullptr);
      const real_t dt = ComputeMaxDt(op.Sample(), pm, p, kCuse[p]);
      op.SetTimeStep(dt);
      auto s0 = [p0](const real_t *, real_t s[9]) {
         for (int i = 0; i < 9; i++) { s[i] = 0.0; }
         s[0] = s[4] = s[8] = -p0;
      };
      auto rc = MakeCoupler(op, DPStack(1.0e6, 0.5, tv), s0);
      op.SetRheology(rc.get());
      ParFiniteElementSpace &ufes = op.Elasticity().GetFESpace();
      ParGridFunction U(&ufes), Ud(&ufes);
      U = 0.0;
      VectorFunctionCoefficient vc(3, [gdot](const Vector &x, Vector &v) {
         v.SetSize(3);
         v = 0.0;
         v(0) = gdot * x(1);
      });
      Ud.ProjectCoefficient(vc);
      real_t t = 0.0;
      op.SetTime(t);
      CentralDifferenceSolver ode;
      ode.Init(op);
      // the pointwise model, same context as node (0, 0) (all nodes alike)
      const RheologyField &F = rc->Field();
      const RheologyStack &st = rc->Stack();
      std::vector<real_t> q(st.NumStateVars());
      st.InitState(F.Context(0, 0), q.data());
      const long nsteps = std::lround(std::ceil(8.0 / dt));
      real_t worst = 0.0, scale = p0;
      long long yielded = 0;
      for (long s = 0; s < nsteps; s++)
      {
         real_t h = dt;
         ode.Step(U, Ud, t, h);
         RheologyPointContext ctx = F.Context(0, 0);
         ctx.t = t;
         const real_t eps[6] = {0.0, 0.0, 0.0, 0.5 * gdot * t, 0.0, 0.0};
         real_t ds[6];
         st.Update(ctx, eps, dt, q.data(), ds);
         yielded = std::max(yielded, rc->YieldingNodes());
         for (int e = 0; e < F.NumElements(); e++)
         {
            for (int i = 0; i < F.NodesPerElement(); i++)
            {
               for (int c = 0; c < 6; c++)
               {
                  worst = std::max(worst, std::abs(F.Dsig(e, i)[c] - ds[c]));
               }
            }
         }
      }
      worst = RMax(worst, MPI_COMM_WORLD);
      const real_t sxy_el = TPV102Params::mu * gdot * t;
      if (g_rank == 0)
      {
         std::cout << std::setprecision(4) << "  Tv = " << tv << " s: " << nsteps
                   << " steps to t = " << t << " s, elastic shear stress would be " << sxy_el
                   << " Pa, final dsig_xy " << F.Dsig(0, 0)[3] << " Pa, most yielding nodes in a step "
                   << yielded << "; max |dsig - pointwise| / p0 = " << worst / scale << "\n";
      }
      TEST_ASSERT(yielded > 0, "O6 Tv = " << tv << " the block yields");
      TEST_ASSERT(worst <= 1e-10 * scale, "O6 Tv = " << tv << " stress path = pointwise model");
   }
}

// --------------------------------------------------------------------- O7
static void TestO7()
{
   if (g_rank == 0) { std::cout << "\n=== O7 refusals ===\n"; }
   RheologyConfig on;
   on.drucker_prager = true;
   on.dp.cohesion = 1.36e6;
   on.dp.bulk_friction = 0.1934;
   on.dp.relaxation_time = 0.03;
   RheologyConfig off;
   auto has = [](const std::string &e, const char *k) { return e.find(k) != std::string::npos; };
   TEST_ASSERT(CheckRheologyRun(on, false, false, false).empty(), "O7 plain run accepted");
   TEST_ASSERT(CheckRheologyRun(off, true, true, false).empty(), "O7 nothing enabled: no refusal");
   TEST_ASSERT(has(CheckRheologyRun(on, true, false, true), "--lts rate2"), "O7 --lts rate2 refused");
   TEST_ASSERT(has(CheckRheologyRun(on, false, true, true), "energy_balance"),
               "O7 energy_balance refused");
   RheologyConfig out = on;
   out.output_dt = 0.5;
   TEST_ASSERT(has(CheckRheologyRun(out, false, false, false), "output_dt_s"),
               "O7 output_dt_s without parallel HDF5 refused");
   TEST_ASSERT(CheckRheologyRun(out, false, false, true).empty(),
               "O7 output_dt_s with parallel HDF5 accepted");
   TEST_ASSERT(CheckInitialStressRun(false, false, true, false).empty(),
               "O7 no [initial_stress]: no refusal");
   TEST_ASSERT(has(CheckInitialStressRun(true, false, false, true), "[stress]"),
               "O7 [stress] next to [initial_stress] refused");
   TEST_ASSERT(has(CheckInitialStressRun(true, true, true, true), "[[stress_region]]"),
               "O7 [[stress_region]] next to [initial_stress] refused");
   TEST_ASSERT(has(CheckInitialStressRun(true, true, false, false), "P_p"),
               "O7 P_p next to [initial_stress] refused");
   TEST_ASSERT(CheckInitialStressRun(true, true, false, true).empty(),
               "O7 [initial_stress] alone accepted");
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   if (g_rank == 0)
   {
      std::cout << "=== Nitsche off-fault rheology: operator tests O1-O7 on "
                << Mpi::WorldSize() << " rank(s) ===\n";
   }
   TestO7();
   TestO3O4();
   TestO2();
   TestO6();
   TestO1();
   TestO5();
   if (g_rank == 0)
   {
      std::cout << "\n" << num_passed << " passed, " << num_failed << " failed, "
                << num_tests << " total\n";
   }
   return num_failed == 0 ? 0 : 1;
}

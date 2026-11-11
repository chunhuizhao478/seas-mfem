// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// tests/2d_inplane/test_qd_implicit_green.cpp — option I, phase I2
// (cycle/qd_implicit.hpp, PlanarProblem): the implicit fault system on the
// real planar Green's function of the miniature BP3 layout of
// test_fault_contact_qd.cpp (ShearBox, dip 60, p = 1 and 2, stage operators
// with slip_comps = 3, the dense G packed planar), plan "Option I:
// implementation and test plan", phase I2 (T0 on the real G).
//
//   G1  tractions: QDImplicitFault's tau and sigma at a slip field and
//       t = 100 yr equal the production stage (GreenStage) plus the
//       pre-stress to 1e-12; the shear diagonal of G is positive (slip
//       relieves its own shear) on every point.
//   G2  gap columns: PlanarProblem's gap_column(l) equals the solve path's
//       QDUnitColumn of a unit opening at l bit for bit, twice the same, and
//       the opening raises its own sigma (gamma > 0).
//   G3  the Newton matrix matches central differences of F and the
//       preconditioner inverts it, with an open point and a closed point
//       with g > 0, no floor / taper / hard floor, two shifts.
//   G4  consistent start with a tensile patch: g >= 0, sigma >= 0,
//       g sigma = 0, and gap columns for exactly the open points.
//   G5  one event (velocity-weakening top, strengthening below, a faster
//       patch to start it): ARK3 and ARK4 (atol 1e-7) against RK 5dp
//       (atol 1e-11) on the production stage + the kernel's friction
//       solve; event time within 1%, slip after the event within 1%.
//   G6  two shallow points pulled into tension and back (test-only extra
//       load): sigma >= 0, g >= 0, g sigma = 0 at every accepted state, psi
//       frozen while open, closed with g = 0 at the end; the same run on
//       all ranks and on one rank (the gathered G and columns) bit for bit.
//   G7  cost sign: one Jacobian (assembly + dense LU) and one solve at
//       N = 2,638 (the UCERF Type-A count), 0 and 50 open points, on one
//       rank.
//   G8  a non-finite Newton iterate (NaN gap at a point without a gap
//       column, NaN slip): Residual, LocalResidual, Evaluate and Jacobian
//       return non-finite values on every rank instead of stopping the run,
//       so the solve fails and TS rejects the step (UCERF T4b, Vista
//       1133642, stopped in Tractions); a finite state stays finite.
//
// Usage (from miniapps/seas; PETSc and MUMPS needed):
//   mpirun -np N ./seas_test_qd_implicit_green

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../fault/fault_frame_field.hpp"
#include "../../constitutive/linear_elastic.hpp"
#include "../../cycle/qd_implicit.hpp"
#include "../../friction/device_rate_state.hpp"

#include <cmath>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;
using mfem::seas::cycle::QDImplicitFault;
using mfem::seas::cycle::QDImplicitSettings;
using mfem::seas::cycle::PlanarFaultInputs;

namespace
{

int rank = 0, nranks = 1, fails = 0;
std::string only;      ///< -only <letters>: run these tests only, e.g. -only 56 (development)
bool trace = false;    ///< -trace: one line per accepted step of G6 (development)
bool newton_local = QDImplicitSettings().newton_local;   ///< -qd_newton_local (cost study)
std::string dense_lu = "serial";   ///< -qd_dense_lu shared: the whole suite on the shared LU

bool Wanted(char c) { return only.empty() || only.find(c) != std::string::npos; }

void Check(bool ok, const std::string &what)
{
   int all = ok ? 1 : 0;
   MPI_Allreduce(MPI_IN_PLACE, &all, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
   if (rank == 0) { std::cout << (all ? "  PASSED: " : "  FAILED: ") << what << std::endl; }
   if (!all) { fails++; }
}

void Say(const std::string &s) { if (rank == 0) { std::cout << s << std::endl; } }

std::string Fmt(double x, int prec = 6)
{
   std::ostringstream o;
   o << std::setprecision(prec) << x;
   return o.str();
}

double MaxAll(double v, MPI_Comm c = MPI_COMM_WORLD)
{
   MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPI_DOUBLE, MPI_MAX, c);
   return v;
}

std::vector<real_t> GatherAll(MPI_Comm comm, const real_t *mine, int nmine)
{
   int p = 1;
   MPI_Comm_size(comm, &p);
   std::vector<int> c(p), d(p);
   MPI_Allgather(&nmine, 1, MPI_INT, c.data(), 1, MPI_INT, comm);
   int tot = 0;
   for (int r = 0; r < p; r++) { d[r] = tot; tot += c[r]; }
   std::vector<real_t> all(tot);
   MPI_Allgatherv(mine, nmine, MPITypeMap<real_t>::mpi_type, all.data(), c.data(), d.data(),
                  MPITypeMap<real_t>::mpi_type, comm);
   return all;
}

bool AllFinite(const std::vector<real_t> &v)
{
   for (real_t x : v) { if (!std::isfinite(x)) { return false; } }
   return true;
}

// ---- the miniature BP3 layout (test_fault_contact_qd.cpp) ---------------------

const real_t kLambda = 32.04e9, kMu = 32.04e9, kRho = 2670.0;
const real_t kLs = 20e3, kD = 20e3, kWf = 12e3;
const int kNs = 8, kNy = 10;

Mesh *ShearBox(real_t dip_deg)
{
   const real_t c = (dip_deg == 90.0) ? 0.0 : 1.0 / std::tan(dip_deg * M_PI / 180.0);
   const int nvx = kNs + 1, nvy = kNy + 1;
   Mesh *m = new Mesh(2, nvx * nvy, 2 * kNs * kNy, 0, 2);
   for (int j = 0; j < nvy; j++)
      for (int i = 0; i < nvx; i++)
      {
         const real_t s = -kLs + 2.0 * kLs * i / kNs, y = -kD + kD * j / kNy;
         m->AddVertex(s - y * c, y);
      }
   auto v = [&](int i, int j) { return i + nvx * j; };
   for (int j = 0; j < kNy; j++)
      for (int i = 0; i < kNs; i++)
      {
         m->AddTriangle(v(i, j), v(i + 1, j), v(i + 1, j + 1), 1);
         m->AddTriangle(v(i, j), v(i + 1, j + 1), v(i, j + 1), 1);
      }
   for (int i = 0; i < kNs; i++)
   {
      m->AddBdrSegment(v(i, 0), v(i + 1, 0), 5);
      m->AddBdrSegment(v(i, kNy), v(i + 1, kNy), 1);
   }
   for (int j = 0; j < kNy; j++)
   {
      m->AddBdrSegment(v(0, j), v(0, j + 1), 5);
      m->AddBdrSegment(v(kNs, j), v(kNs, j + 1), 5);
      const real_t ymid = -kD + kD * (j + 0.5) / kNy;
      m->AddBdrSegment(v(kNs / 2, j), v(kNs / 2, j + 1), (-ymid < kWf) ? 3 : 5);
   }
   m->FinalizeTopology();
   m->Finalize(false, true);
   m->SetAttributes();
   return m;
}

std::shared_ptr<const FaultFrameConfig> FrameConfig(const real_t ref[3], const real_t up[3])
{
   auto c = std::make_shared<FaultFrameConfig>();
   for (int d = 0; d < 3; d++)
   {
      c->orientation.ref_normal[d] = ref[d];
      c->orientation.up[d] = up[d];
      c->up[d] = up[d];
   }
   c->orientation.fault_attrs = {3};
   return c;
}

// Linear in t, so g(0) = 0 and the Green build accepts it.
void PlateLoading2D(const Vector &x, real_t t, Vector &u)
{
   const real_t n[2] = {std::sqrt(3.0) / 2.0, 0.5}, sd[2] = {-0.5, std::sqrt(3.0) / 2.0};
   const real_t Vp = 1e-9, d = x(0) * n[0] + x(1) * n[1];
   const real_t Vh = (d > 1.0) ? 0.5 * Vp * t : (d < -1.0) ? -0.5 * Vp * t : Vp * t;
   u.SetSize(2);
   u(0) = Vh * sd[0];
   u(1) = Vh * sd[1];
}

using Op = ElasticityDomainOperator<ParMesh>;

/// The operator with its stage operators (slip_comps 3) and the planar G.
struct Box
{
   std::unique_ptr<ParMesh> pm;
   std::unique_ptr<LinearElastic> model;
   std::unique_ptr<Op> op;
   int n = 0, slot = 0, N = 0;
   long long first = 0;
   Vector depth;   ///< owned, m
};

std::unique_ptr<Box> MakeBox(MPI_Comm comm, int order)
{
   static const real_t ref[3] = {1, 0, 0}, up[3] = {0, 1, 0};
   auto B = std::make_unique<Box>();
   std::unique_ptr<Mesh> ser(ShearBox(60.0));
   B->pm = std::make_unique<ParMesh>(comm, *ser);
   B->model = std::make_unique<LinearElastic>(kLambda, kMu);
   DomainConfig dc;
   FaultFrameSource src;
   src.config = FrameConfig(ref, up);
   dc.fault_frames = src;
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.natural_attrs = {1};
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = PlateLoading2D;
   B->op = std::make_unique<Op>(*B->pm, order, *B->model, 0.0, kWf, 2.0 * kLs, bc,
                                DGMethod::IP, SolverType::MUMPS, dc);
   B->op->BuildQDStageOperators(true, 3);
   B->op->BuildQDGreen(QDGreenConfig{});
   B->n = B->op->GetNumOwnedFaultDOFs();
   B->slot = B->op->FrameField()->InPlaneSlot();
   const QDGreenOperator *G = B->op->GetQDGreen();
   MFEM_VERIFY(G && G->Planar(), "the planar Green's function");
   B->N = static_cast<int>(G->NumFaultDofs());
   B->first = G->FirstOwned();
   Vector dl;
   B->op->GetFaultDepths(dl);
   B->op->RestrictToOwnedFault(dl, B->depth, 1);
   MFEM_VERIFY(B->depth.Size() == B->n, "owned depths");
   return B;
}

cycle::UnitResponse UnitG(Op *op)
{
   return [op](const Vector &e, Vector &T, Vector &Tn) { op->QDUnitColumn(e, T, &Tn); };
}

// ---- physics: BP3-like rate and state on the miniature fault ------------------

struct Phys
{
   double a_vw = 0.01, a_vs = 0.025, b = 0.015, Dc = 0.01, V0 = 1e-6, f0 = 0.6;
   double sigma = 50e6, Vpl = 1e-9;
   double eta() const { return kMu / (2.0 * std::sqrt(kMu / kRho)); }
   /// a(depth): weakening to 6 km, strengthening below 9 km, linear between
   double a(double z) const
   {
      if (z <= 6e3) { return a_vw; }
      if (z >= 9e3) { return a_vs; }
      return a_vw + (a_vs - a_vw) * (z - 6e3) / 3e3;
   }
   double psi_ss(double V) const { return f0 - b * std::log(V / V0); }
};
const Phys ph;

/// The kernel's signed slip rate: tau + sigma f(V, psi) + eta V = 0, V
/// anti-parallel to tau.
double SlipRate(double tau, double psi, double sigma, double a)
{
   int st = 0;
   const double Va = RSSolveSlipRatePsi(std::abs(tau), psi, std::max(sigma, 0.0), ph.eta(), a,
                                        ph.V0, st);
   return (tau > 0.0) ? -Va : ((tau < 0.0) ? Va : 0.0);
}

/// Inputs with tau_pre = s (sigma0 f(Vi, psi0) + eta Vi) (so |V(0)| = Vi
/// at delta = 0, t = 0) and the given normal pre-stress.
PlanarFaultInputs Inputs(const Box &B, double s,
                         const std::function<double(long long j, double z)> &Vi,
                         const std::function<double(long long j, double z)> &sigma0)
{
   PlanarFaultInputs in;
   const int n = B.n;
   in.tau_pre.SetSize(2 * n);
   in.tau_pre = 0.0;
   in.sigma0.SetSize(n);
   in.a.SetSize(n);
   in.eta.SetSize(n);
   in.Dc.SetSize(n);
   const double psi0 = ph.psi_ss(ph.Vpl);
   for (int o = 0; o < n; o++)
   {
      const long long j = B.first + o;
      const double z = B.depth(o), a = ph.a(z), sg = sigma0(j, z), V = Vi(j, z);
      in.sigma0(o) = sg;
      in.a(o) = a;
      in.eta(o) = ph.eta();
      in.Dc(o) = ph.Dc;
      const double fr = RSFrictionCoefficientPsi(V, psi0, a, ph.V0);
      in.tau_pre(2 * o + B.slot) = s * (std::max(sg, 0.0) * fr + ph.eta() * V);
   }
   in.b = ph.b;
   in.f0 = ph.f0;
   in.V0 = ph.V0;
   return in;
}

/// The loading's sign: -1 when T_D drives the shear negative (then a
/// negative tau_pre slips forward, V > 0, and the loading adds to it).
double LoadSign(const QDImplicitFault::Problem &P)
{
   double s = 0.0;
   for (real_t v : P.TDt) { s += v; }
   return (s < 0.0) ? -1.0 : 1.0;
}

// ---- the TS driver (as test_qd_esdirk_points.cpp) ------------------------------

struct RunSpec
{
   double t0 = 0.0, h0 = 1.0, tmax = 1.0;
   PetscInt max_steps = 1000000;
   std::function<void(TS, QDImplicitFault &, double, const std::vector<real_t> &)> post;
};

struct RunOut
{
   TSConvergedReason reason = TS_CONVERGED_ITERATING;
   PetscErrorCode err = PETSC_SUCCESS;
   PetscInt steps = 0, rejections = 0, snes_fails = 0, snes_its = 0;
   double t_end = 0.0, seconds = 0.0;
   std::vector<real_t> final_all;
   QDImplicitFault::Counts counts;
};

struct Ctx
{
   QDImplicitFault *f = nullptr;
   const RunSpec *spec = nullptr;
};

PetscErrorCode PostStep(TS ts)
{
   PetscFunctionBeginUser;
   void *c = nullptr;
   PetscCall(TSGetApplicationContext(ts, &c));
   Ctx *ctx = static_cast<Ctx *>(c);
   PetscReal t;
   Vec U;
   PetscCall(TSGetTime(ts, &t));
   PetscCall(TSGetSolution(ts, &U));
   const PetscScalar *u;
   PetscInt n;
   PetscCall(VecGetLocalSize(U, &n));
   PetscCall(VecGetArrayRead(U, &u));
   const std::vector<real_t> all = GatherAll(ctx->f->Comm(), reinterpret_cast<const real_t *>(u),
                                             static_cast<int>(n));
   PetscCall(VecRestoreArrayRead(U, &u));
   if (ctx->spec->post) { ctx->spec->post(ts, *ctx->f, static_cast<double>(t), all); }
   PetscFunctionReturn(PETSC_SUCCESS);
}

#define TSCHK(call) MFEM_VERIFY((call) == PETSC_SUCCESS, #call)

RunOut Run(QDImplicitFault &f, const std::vector<real_t> &U0own, const RunSpec &spec)
{
   RunOut r;
   MPI_Comm comm = f.Comm();
   const int nl = QDImplicitFault::kPer * f.NumOwned();
   MFEM_VERIFY(U0own.size() == size_t(nl), "Run: U0 has " << U0own.size() << " entries");
   Vec U;
   TS ts;
   TSCHK(VecCreateMPI(comm, nl, PETSC_DETERMINE, &U));
   PetscScalar *u;
   TSCHK(VecGetArray(U, &u));
   for (int i = 0; i < nl; i++) { u[i] = U0own[i]; }
   TSCHK(VecRestoreArray(U, &u));
   TSCHK(TSCreate(comm, &ts));
   TSCHK(f.SetupTS(ts));
   TSCHK(f.SetTolerances(ts));
   TSCHK(TSSetTime(ts, spec.t0));
   TSCHK(TSSetTimeStep(ts, spec.h0));
   TSCHK(TSSetMaxTime(ts, spec.tmax));
   TSCHK(TSSetMaxSteps(ts, spec.max_steps));
   TSCHK(TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP));
   TSCHK(TSSetMaxStepRejections(ts, -1));
   TSAdapt adapt;
   TSCHK(TSGetAdapt(ts, &adapt));
   TSCHK(TSAdaptSetType(adapt, TSADAPTBASIC));
   TSCHK(TSSetFromOptions(ts));
   Ctx ctx;
   ctx.f = &f;
   ctx.spec = &spec;
   TSCHK(TSSetApplicationContext(ts, &ctx));
   TSCHK(TSSetPostStep(ts, PostStep));
   TSCHK(TSSetSolution(ts, U));
   TSCHK(TSSetUp(ts));
   f.ResetCounts();
   const double t_start = MPI_Wtime();
   PetscPushErrorHandler(PetscReturnErrorHandler, nullptr);
   r.err = TSSolve(ts, U);
   PetscPopErrorHandler();
   r.seconds = MPI_Wtime() - t_start;
   PetscReal t;
   TSGetConvergedReason(ts, &r.reason);
   TSGetStepNumber(ts, &r.steps);
   TSGetStepRejections(ts, &r.rejections);
   TSGetSNESFailures(ts, &r.snes_fails);
   TSGetSNESIterations(ts, &r.snes_its);
   TSGetTime(ts, &t);
   r.t_end = static_cast<double>(t);
   const PetscScalar *uu;
   VecGetArrayRead(U, &uu);
   r.final_all = GatherAll(comm, reinterpret_cast<const real_t *>(uu), nl);
   VecRestoreArrayRead(U, &uu);
   r.counts = f.GetCounts();
   TSDestroy(&ts);
   VecDestroy(&U);
   return r;
}

std::string Summary(const RunOut &r)
{
   std::ostringstream o;
   o << "reason " << TSConvergedReasons[r.reason] << ", t " << std::setprecision(8) << r.t_end
     << ", steps " << r.steps << ", rejections " << r.rejections << ", SNES failures "
     << r.snes_fails << ", Newton its " << r.snes_its << ", factors " << r.counts.factors
     << ", residuals " << r.counts.residuals << ", clipped its " << r.counts.clipped_its
     << ", branch changes " << r.counts.branch_changes << ", gap columns "
     << r.counts.columns << ", most open " << r.counts.max_open << ", plain updates "
     << r.counts.plain_applies << ", "
     << std::setprecision(3) << r.seconds << " s";
   return o.str();
}

std::vector<real_t> Owned(const QDImplicitFault &f, const std::vector<real_t> &all)
{
   const int k = QDImplicitFault::kPer;
   return std::vector<real_t>(all.begin() + k * f.FirstOwned(),
                              all.begin() + k * (f.FirstOwned() + f.NumOwned()));
}

/// Largest |V| of the state (global, collective).
double VMax(QDImplicitFault &f, double t, const std::vector<real_t> &own)
{
   std::vector<real_t> V, tau, sig;
   std::vector<char> open;
   f.Evaluate(t, own.data(), V, tau, sig, open);
   double v = 0.0;
   for (real_t x : V) { v = std::max(v, std::abs(static_cast<double>(x))); }
   return MaxAll(v, f.Comm());
}

// ---- G1, G2 ------------------------------------------------------------------------

void TestG1G2(Box &B, int p)
{
   MPI_Comm comm = MPI_COMM_WORLD;
   const QDGreenOperator &G = *B.op->GetQDGreen();
   const int n = B.n, N = B.N;
   PlanarFaultInputs in = Inputs(B, -1.0, [](long long, double) { return ph.Vpl; },
                                 [](long long j, double) { return 40e6 + 1e5 * j; });
   QDImplicitFault::Problem P = cycle::PlanarProblem(comm, G, in, UnitG(B.op.get()));
   Say("  N = " + std::to_string(N) + " fault points, in-plane slot " + std::to_string(B.slot));

   // G1: tractions against the production stage
   {
      QDImplicitFault f(comm, P, QDImplicitSettings());
      const double t = 100.0 * 365.25 * 86400.0;
      Vector slip(2 * n);
      slip = 0.0;
      std::vector<real_t> U(3 * n, 0.0);
      for (int o = 0; o < n; o++)
      {
         const double j = static_cast<double>(B.first + o);
         slip(2 * o + B.slot) = std::sin(0.37 * j + 0.1) + 0.5 * std::cos(0.11 * j);
         U[3 * o] = slip(2 * o + B.slot);
         U[3 * o + 1] = 0.7;
      }
      std::vector<real_t> V, tau, sig;
      std::vector<char> open;
      f.Evaluate(t, U.data(), V, tau, sig, open);
      Vector T, Tn;
      B.op->GreenStage(t, slip, T, &Tn);
      double dt = 0.0, ds = 0.0, st = 0.0, ss = 0.0;
      for (int o = 0; o < n; o++)
      {
         const double rt = in.tau_pre(2 * o + B.slot) + T(2 * o + B.slot);
         const double rs = in.sigma0(o) + Tn(o);
         dt = std::max(dt, std::abs(tau[o] - rt));
         ds = std::max(ds, std::abs(sig[o] - rs));
         st = std::max(st, std::abs(rt));
         ss = std::max(ss, std::abs(rs));
      }
      dt = MaxAll(dt);
      ds = MaxAll(ds);
      st = MaxAll(st);
      ss = MaxAll(ss);
      double dmin = std::numeric_limits<double>::infinity(), dmax = 0.0, tdmin = dmin,
             tdmax = -dmin, gsd = 0.0;
      for (int j = 0; j < N; j++)
      {
         const double d = P.Gt[size_t(j) * N + j];
         dmin = std::min(dmin, d);
         dmax = std::max(dmax, d);
         tdmin = std::min(tdmin, static_cast<double>(P.TDt[j]));
         tdmax = std::max(tdmax, static_cast<double>(P.TDt[j]));
         gsd = std::max(gsd, std::abs(static_cast<double>(P.Gs[size_t(j) * N + j])));
      }
      Say("  G1 p = " + std::to_string(p) + ": |tau - stage| " + Fmt(dt, 3) + " of " +
          Fmt(st, 3) + " Pa, |sigma - stage| " + Fmt(ds, 3) + " of " + Fmt(ss, 3) +
          " Pa; shear diagonal " + Fmt(dmin, 4) + " .. " + Fmt(dmax, 4) + " Pa/m, |normal "
          "diagonal| <= " + Fmt(gsd, 3) + " Pa/m, T_D shear " + Fmt(tdmin, 4) + " .. " +
          Fmt(tdmax, 4) + " Pa/s");
      Check(dt <= 1e-12 * st && ds <= 1e-12 * ss, "G1 p = " + std::to_string(p) +
            ": tau and sigma = the production stage + pre-stress to 1e-12");
      Check(dmin > 0.0, "G1 p = " + std::to_string(p) + ": the shear diagonal of G is "
            "positive at every point (slip relieves its own shear)");
   }

   // G2: gap columns against the solve path
   {
      bool same = true, twice = true, gam_pos = true;
      std::string line;
      for (int l : {0, N / 2, N - 1})
      {
         std::vector<real_t> t1, s1, t2, s2;
         P.gap_column(l, t1, s1);
         P.gap_column(l, t2, s2);
         Vector e(3 * n), T, Tn;
         e = 0.0;
         if (B.first <= l && l < B.first + n) { e(3 * (l - static_cast<int>(B.first)) + 2) = 1.0; }
         B.op->QDUnitColumn(e, T, &Tn);
         std::vector<real_t> mt(n), ms(n);
         for (int o = 0; o < n; o++) { mt[o] = T(2 * o + B.slot); ms[o] = Tn(o); }
         const std::vector<real_t> rt = GatherAll(comm, mt.data(), n);
         const std::vector<real_t> rs = GatherAll(comm, ms.data(), n);
         same = same && (rt == t1) && (rs == s1);
         twice = twice && (t1 == t2) && (s1 == s2);
         gam_pos = gam_pos && s1[l] > 0.0;
         line += " " + std::to_string(l) + ": gamma " + Fmt(s1[l], 4) + " Pa/m, shear " +
                 Fmt(t1[l], 4) + " Pa/m;";
      }
      Say("  G2 p = " + std::to_string(p) + ":" + line);
      Check(same, "G2 p = " + std::to_string(p) + ": gap_column = QDUnitColumn bit for bit");
      Check(twice, "G2 p = " + std::to_string(p) + ": a second fetch gives the same bits");
      Check(gam_pos, "G2 p = " + std::to_string(p) + ": an opening raises its own sigma "
            "(gamma > 0)");
   }
}

// ---- G3, G4 ------------------------------------------------------------------------

void TestG3G4(Box &B, int p)
{
   MPI_Comm comm = MPI_COMM_WORLD;
   const QDGreenOperator &G = *B.op->GetQDGreen();
   const int n = B.n, N = B.N;
   const int ja = N / 2, jb = N / 2 + 2;   // ja open, jb closed with g > 0
   const std::string ps = "p = " + std::to_string(p);

   // G3: Newton matrix and preconditioner
   struct Mode { double floor; bool taper; const char *name; };
   for (const Mode m : {Mode{0.0, false, "no floor"}, Mode{10e6, true, "floor 10 MPa, taper"},
                        Mode{10e6, false, "hard floor 10 MPa"}})
   {
      for (const double shift : {2.3e-3, 10.0})
      for (const bool local : {false, true})
      {
         // local: Newton's own residual G (each point's stage solve) and
         // its matrix; else F and dF/dU + shift dF/dUdot
         PlanarFaultInputs in = Inputs(B, -1.0,
                                       [](long long j, double) { return 1e-8 * (1.0 + 0.3 * j); },
                                       [=](long long j, double) { return j == ja ? -3e5 : 20e6; });
         in.sigma_floor = m.floor;
         in.floor_taper = m.taper;
         QDImplicitSettings set;
         set.newton_local = local;
         set.dense_lu = dense_lu;
         QDImplicitFault f(comm, cycle::PlanarProblem(comm, G, in, UnitG(B.op.get())), set);
         f.EnsureColumn(ja);
         f.EnsureColumn(jb);
         const double gam_a = f.Gamma(ja);
         std::vector<real_t> Uall(3 * N), Udall(3 * N), xall(3 * N), ball(3 * N);
         for (int j = 0; j < N; j++)
         {
            Uall[3 * j] = 1e-4 * std::sin(0.5 * j + 0.2);
            Uall[3 * j + 1] = ph.psi_ss(ph.Vpl) + 0.002 * std::sin(0.9 * j);
            Uall[3 * j + 2] = (j == ja) ? 2.0 * 3e5 / gam_a : ((j == jb) ? 1e-6 : 0.0);
            Udall[3 * j] = 1e-6 * std::sin(0.4 * j + 0.1);
            Udall[3 * j + 1] = 1e-4 * std::cos(0.7 * j);
            Udall[3 * j + 2] = 1e-5 * std::sin(1.1 * j + 0.3);
            xall[3 * j] = 1e-4 * std::sin(1.3 * j + 0.7);
            xall[3 * j + 1] = 1e-3 * std::cos(0.9 * j + 0.2);
            xall[3 * j + 2] = (j == ja || j == jb) ? 1e-4 * (1.0 + 0.1 * j) : 0.0;
            ball[3 * j] = 1e-6 * std::cos(0.3 * j + 1.0);
            ball[3 * j + 1] = 1e-5 * std::sin(0.8 * j);
            ball[3 * j + 2] = (j == ja || j == jb) ? 3e3 * (1.0 + 0.2 * j) : 0.0;
         }
         const std::vector<real_t> U = Owned(f, Uall), Ud = Owned(f, Udall), x = Owned(f, xall),
                                   b = Owned(f, ball);
         std::vector<real_t> V, tau, sig;
         std::vector<char> open;
         f.Evaluate(0.0, U.data(), V, tau, sig, open);
         std::vector<real_t> openr(open.begin(), open.end());
         const std::vector<real_t> open_all = GatherAll(comm, openr.data(), n);
         double vmin = std::numeric_limits<double>::infinity();
         for (real_t v : V) { vmin = std::min(vmin, std::abs(static_cast<double>(v))); }
         vmin = -MaxAll(-vmin);
         const bool branches = open_all[ja] == 1.0 && open_all[jb] == 0.0;
         f.Jacobian(0.0, U.data(), Ud.data(), shift);
         std::vector<real_t> Jx(3 * n), Fp(3 * n), Fm(3 * n), Up(3 * n), Um(3 * n),
             Udp(3 * n), Udm(3 * n);
         if (local) { f.ApplyLocalJacobian(x.data(), Jx.data()); }
         else { f.ApplyJacobian(x.data(), Jx.data()); }
         const double eps = 1e-4;
         for (int i = 0; i < 3 * n; i++)
         {
            Up[i] = U[i] + eps * x[i];
            Um[i] = U[i] - eps * x[i];
            Udp[i] = Ud[i] + shift * eps * x[i];
            Udm[i] = Ud[i] - shift * eps * x[i];
         }
         if (local)
         {
            f.LocalResidual(0.0, Up.data(), Udp.data(), shift, Fp.data());
            f.LocalResidual(0.0, Um.data(), Udm.data(), shift, Fm.data());
         }
         else
         {
            f.Residual(0.0, Up.data(), Udp.data(), Fp.data());
            f.Residual(0.0, Um.data(), Udm.data(), Fm.data());
         }
         double num = 0.0, den = 0.0;
         for (int i = 0; i < 3 * n; i++)
         {
            num = std::max(num, std::abs((Fp[i] - Fm[i]) / (2.0 * eps) - Jx[i]));
            den = std::max(den, std::abs(static_cast<double>(Jx[i])));
         }
         const double rel = MaxAll(num) / MaxAll(den);
         std::vector<real_t> y(3 * n), Jy(3 * n);
         if (local)
         {
            f.ApplyLocalInverse(b.data(), y.data());
            f.ApplyLocalJacobian(y.data(), Jy.data());
         }
         else
         {
            f.ApplyInverse(b.data(), y.data());
            f.ApplyJacobian(y.data(), Jy.data());
         }
         double pn = 0.0, pd = 0.0;
         for (int i = 0; i < 3 * n; i++)
         {
            pn = std::max(pn, std::abs(static_cast<double>(Jy[i] - b[i])));
            pd = std::max(pd, std::abs(static_cast<double>(b[i])));
         }
         const double prel = MaxAll(pn) / MaxAll(pd);
         const std::string w = ps + ", " + m.name + ", shift " + Fmt(shift, 3) +
                               (local ? ", G (point stage solves)" : ", F");
         Say("    G3 " + w + ": min |V| " + Fmt(vmin, 3) + " m/s, |J x - FD| / |J x| = " +
             Fmt(rel, 3) + ", |J P b - b| / |b| = " + Fmt(prel, 3));
         Check(branches && vmin > 0.0, "G3 " + w + ": point " + std::to_string(ja) +
               " open, " + std::to_string(jb) + " closed with g > 0, every V nonzero");
         Check(rel < 1e-6, "G3 " + w + ": Newton matrix = central differences to 1e-6");
         Check(prel < 1e-10, "G3 " + w + ": the preconditioner inverts it to 1e-10");
      }
   }

   // G4: consistent start with a tensile patch of 3 points
   {
      const long long p0 = N / 2 - 1, p1 = N / 2 + 1;
      PlanarFaultInputs in = Inputs(B, -1.0, [](long long, double) { return ph.Vpl; },
                                    [=](long long j, double)
      { return (j >= p0 && j <= p1) ? -2e6 : 30e6; });
      QDImplicitFault f(comm, cycle::PlanarProblem(comm, G, in, UnitG(B.op.get())),
                        QDImplicitSettings());
      std::vector<real_t> U(3 * n, 0.0);
      for (int o = 0; o < n; o++) { U[3 * o + 1] = ph.psi_ss(ph.Vpl); }
      const long opened = f.ConsistentStart(0.0, U.data());
      std::vector<real_t> V, tau, sig;
      std::vector<char> open;
      f.Evaluate(0.0, U.data(), V, tau, sig, open);
      std::vector<real_t> g(n);
      for (int o = 0; o < n; o++) { g[o] = U[3 * o + 2]; }
      const std::vector<real_t> g_all = GatherAll(comm, g.data(), n);
      const std::vector<real_t> s_all = GatherAll(comm, sig.data(), n);
      double min_sig = 1e30, min_g = 1e30, max_gs = 0.0;
      long n_patch = 0;
      bool cols = true;
      for (int j = 0; j < N; j++)
      {
         min_sig = std::min(min_sig, static_cast<double>(s_all[j]));
         min_g = std::min(min_g, static_cast<double>(g_all[j]));
         max_gs = std::max(max_gs, std::abs(static_cast<double>(g_all[j] * s_all[j])));
         if (j >= p0 && j <= p1 && g_all[j] > 0.0) { n_patch++; }
         cols = cols && (f.HasColumn(j) == (g_all[j] > 0.0));
      }
      Say("  G4 " + ps + ": " + std::to_string(opened) + " open (" + std::to_string(n_patch) +
          " of the 3 patch points), min g " + Fmt(min_g, 3) + " m, min sigma " +
          Fmt(min_sig, 3) + " Pa, max |g sigma| " + Fmt(max_gs, 3) + " N/m");
      Check(opened >= 1 && n_patch >= 1 && min_g >= 0.0 && min_sig >= -1e-6 * 30e6 &&
            max_gs <= 1e-6 * 30e6 * 1e-2,
            "G4 " + ps + ": g >= 0, sigma >= 0 and g sigma = 0 at every point, the patch open");
      Check(cols, "G4 " + ps + ": gap columns for exactly the open points");
   }
}

// ---- G8: a non-finite iterate fails the solve, it does not stop the run ----------

void TestG8(Box &B, int p)
{
   MPI_Comm comm = MPI_COMM_WORLD;
   const QDGreenOperator &G = *B.op->GetQDGreen();
   const int n = B.n, N = B.N;
   const int jb = N / 2 + 2;   // no gap column is fetched for it
   const std::string ps = "p = " + std::to_string(p);
   PlanarFaultInputs in = Inputs(B, -1.0,
                                 [](long long j, double) { return 1e-8 * (1.0 + 0.3 * j); },
                                 [](long long, double) { return 20e6; });
   for (const bool local : {false, true})
   {
      QDImplicitSettings set;
      set.newton_local = local;
      set.dense_lu = dense_lu;
      QDImplicitFault f(comm, cycle::PlanarProblem(comm, G, in, UnitG(B.op.get())), set);
      const std::string w = ps + (local ? ", newton_local" : ", F");
      std::vector<real_t> Uall(3 * N), Udall(3 * N, 0.0);
      for (int j = 0; j < N; j++)
      {
         Uall[3 * j] = 1e-4 * std::sin(0.5 * j + 0.2);
         Uall[3 * j + 1] = ph.psi_ss(ph.Vpl) + 0.002 * std::sin(0.9 * j);
         Uall[3 * j + 2] = 0.0;
         Udall[3 * j] = 1e-9;
      }
      struct Case { int comp; const char *name; bool bad; };
      for (const Case c : {Case{-1, "finite state", false}, Case{2, "NaN gap", true},
                           Case{0, "NaN slip", true}})
      {
         std::vector<real_t> Ua = Uall;
         if (c.comp >= 0) { Ua[3 * jb + c.comp] = std::numeric_limits<real_t>::quiet_NaN(); }
         const std::vector<real_t> U = Owned(f, Ua), Ud = Owned(f, Udall);
         std::vector<real_t> F(3 * n), Gl(3 * n), V, tau, sig;
         std::vector<char> open;
         f.Residual(0.0, U.data(), Ud.data(), F.data());
         f.LocalResidual(0.0, U.data(), Ud.data(), 10.0, Gl.data());
         f.Evaluate(0.0, U.data(), V, tau, sig, open);
         f.Jacobian(0.0, U.data(), Ud.data(), 10.0);
         // every rank: all slip and psi rows finite (a finite state) or
         // all non-finite (a bad iterate)
         int all_bad = 1, all_fin = 1;
         for (int o = 0; o < n; o++)
         {
            for (const real_t v : {F[3 * o], F[3 * o + 1], Gl[3 * o], Gl[3 * o + 1], V[o],
                                   tau[o], sig[o]})
            {
               all_bad = all_bad && !std::isfinite(v);
               all_fin = all_fin && std::isfinite(v);
            }
         }
         MPI_Allreduce(MPI_IN_PLACE, &all_bad, 1, MPI_INT, MPI_MIN, comm);
         MPI_Allreduce(MPI_IN_PLACE, &all_fin, 1, MPI_INT, MPI_MIN, comm);
         Say("  G8 " + w + ", " + c.name + ": returned; rows " +
             (all_fin ? "all finite" : (all_bad ? "all non-finite" : "mixed")));
         Check(c.bad ? all_bad == 1 : all_fin == 1,
               std::string("G8 ") + w + ", " + c.name + ": " +
               (c.bad ? "non-finite residual on every rank, no stop" : "finite residual"));
      }
   }
}

// ---- G5: one event against RK 5dp --------------------------------------------------

struct RefCtx
{
   Op *op = nullptr;
   const PlanarFaultInputs *in = nullptr;
   int n = 0, slot = 0;
   double vth = 1e-3, t_event = -1.0, vmax = 0.0;
   long steps = 0;
};

/// (delta, psi) per owned point -> slip rate per owned point (the
/// production stage and the kernel's solve).
void RefRates(RefCtx &c, double t, const PetscScalar *u, std::vector<double> &V,
              std::vector<double> &dpsi)
{
   const int n = c.n;
   Vector slip(2 * n), T, Tn;
   slip = 0.0;
   for (int o = 0; o < n; o++) { slip(2 * o + c.slot) = u[2 * o]; }
   c.op->GreenStage(t, slip, T, &Tn);
   V.resize(n);
   dpsi.resize(n);
   const RSDeviceParams prm{ph.V0, ph.b, ph.V0, ph.f0};
   for (int o = 0; o < n; o++)
   {
      const double tau = c.in->tau_pre(2 * o + c.slot) + T(2 * o + c.slot);
      const double sig = c.in->sigma0(o) + Tn(o);
      V[o] = SlipRate(tau, u[2 * o + 1], sig, c.in->a(o));
      dpsi[o] = RSAgingRatePsi(std::abs(V[o]), u[2 * o + 1], c.in->Dc(o), prm);
   }
}

PetscErrorCode RefRHS(TS, PetscReal t, Vec U, Vec F, void *cv)
{
   PetscFunctionBeginUser;
   RefCtx *c = static_cast<RefCtx *>(cv);
   const PetscScalar *u;
   PetscScalar *f;
   PetscCall(VecGetArrayRead(U, &u));
   PetscCall(VecGetArray(F, &f));
   std::vector<double> V, dpsi;
   RefRates(*c, t, u, V, dpsi);
   for (int o = 0; o < c->n; o++) { f[2 * o] = V[o]; f[2 * o + 1] = dpsi[o]; }
   PetscCall(VecRestoreArray(F, &f));
   PetscCall(VecRestoreArrayRead(U, &u));
   PetscFunctionReturn(PETSC_SUCCESS);
}

PetscErrorCode RefPost(TS ts)
{
   PetscFunctionBeginUser;
   void *cv;
   PetscCall(TSGetApplicationContext(ts, &cv));
   RefCtx *c = static_cast<RefCtx *>(cv);
   PetscReal t;
   Vec U;
   PetscCall(TSGetTime(ts, &t));
   PetscCall(TSGetSolution(ts, &U));
   const PetscScalar *u;
   PetscCall(VecGetArrayRead(U, &u));
   std::vector<double> V, dpsi;
   RefRates(*c, t, u, V, dpsi);
   PetscCall(VecRestoreArrayRead(U, &u));
   double vm = 0.0;
   for (double v : V) { vm = std::max(vm, std::abs(v)); }
   vm = MaxAll(vm);
   c->vmax = std::max(c->vmax, vm);
   if (c->t_event < 0.0 && vm >= c->vth) { c->t_event = t; }
   c->steps++;
   PetscFunctionReturn(PETSC_SUCCESS);
}

/// RK 5dp from delta = 0, psi = psi0 to tmax (atol 1e-11); the gathered
/// final (delta, psi).
std::vector<real_t> Reference(RefCtx &c, double tmax)
{
   const int n = c.n;
   Vec U;
   TS ts;
   TSCHK(VecCreateMPI(MPI_COMM_WORLD, 2 * n, PETSC_DETERMINE, &U));
   PetscScalar *u;
   TSCHK(VecGetArray(U, &u));
   for (int o = 0; o < n; o++) { u[2 * o] = 0.0; u[2 * o + 1] = ph.psi_ss(ph.Vpl); }
   TSCHK(VecRestoreArray(U, &u));
   TSCHK(TSCreate(MPI_COMM_WORLD, &ts));
   TSCHK(TSSetType(ts, TSRK));
   TSCHK(TSRKSetType(ts, TSRK5DP));
   TSCHK(TSSetRHSFunction(ts, nullptr, RefRHS, &c));
   TSCHK(TSSetTimeStep(ts, 1e3));
   TSCHK(TSSetMaxTime(ts, tmax));
   TSCHK(TSSetMaxSteps(ts, 10000000));
   TSCHK(TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP));
   TSCHK(TSSetTolerances(ts, 1e-11, nullptr, 0.0, nullptr));
   TSCHK(TSSetMaxStepRejections(ts, -1));
   TSCHK(TSSetFromOptions(ts));
   TSCHK(TSSetApplicationContext(ts, &c));
   TSCHK(TSSetPostStep(ts, RefPost));
   TSCHK(TSSetSolution(ts, U));
   TSCHK(TSSolve(ts, U));
   const PetscScalar *uu;
   TSCHK(VecGetArrayRead(U, &uu));
   const std::vector<real_t> all = GatherAll(MPI_COMM_WORLD, reinterpret_cast<const real_t *>(uu),
                                             2 * n);
   TSCHK(VecRestoreArrayRead(U, &uu));
   TSDestroy(&ts);
   VecDestroy(&U);
   return all;
}

void TestG5(Box &B, int p)
{
   MPI_Comm comm = MPI_COMM_WORLD;
   const QDGreenOperator &G = *B.op->GetQDGreen();
   const int n = B.n, N = B.N;
   const std::string ps = "p = " + std::to_string(p);
   // the loading's sign from G's T_D (a throw-away problem)
   double s = 0.0;
   {
      PlanarFaultInputs in0 = Inputs(B, -1.0, [](long long, double) { return ph.Vpl; },
                                     [](long long, double) { return ph.sigma; });
      s = LoadSign(cycle::PlanarProblem(comm, G, in0, UnitG(B.op.get())));
   }
   // a faster patch at 3-5 km depth starts the event
   const PlanarFaultInputs in = Inputs(B, s, [](long long, double z)
   { return (z >= 3e3 && z <= 5e3) ? 1e3 * ph.Vpl : ph.Vpl; },
   [](long long, double) { return ph.sigma; });
   const double vth = 1e-3;
   // reference, pass 1: the event time (stop well past it)
   RefCtx c1;
   c1.op = B.op.get();
   c1.in = &in;
   c1.n = n;
   c1.slot = B.slot;
   c1.vth = vth;
   const double horizon = 5.0 * 365.25 * 86400.0;
   (void) Reference(c1, horizon);
   Say("  G5 " + ps + ": loading sign " + Fmt(s, 2) + "; RK 5dp (atol 1e-11): |V| >= 1e-3 m/s "
       "first at t = " + Fmt(c1.t_event, 10) + " s, max |V| " + Fmt(c1.vmax, 4) + " m/s, " +
       std::to_string(c1.steps) + " steps to " + Fmt(horizon, 4) + " s");
   Check(c1.t_event > 0.0, "G5 " + ps + ": the reference has an event within 5 yr");
   if (c1.t_event <= 0.0) { return; }
   // pass 2: to the end time, after the event
   const double t_end = 1.5 * c1.t_event;
   RefCtx c2 = c1;
   c2.t_event = -1.0;
   c2.vmax = 0.0;
   c2.steps = 0;
   const std::vector<real_t> ref = Reference(c2, t_end);
   double slip_ref = 0.0;
   for (int j = 0; j < N; j++) { slip_ref = std::max(slip_ref, std::abs(static_cast<double>(ref[2 * j]))); }
   Say("    RK 5dp to t = " + Fmt(t_end, 8) + " s: " + std::to_string(c2.steps) +
       " steps, max |slip| " + Fmt(slip_ref, 6) + " m");
   for (const char *type : {"3", "4"})
   {
      QDImplicitSettings set;
      set.esdirk_type = type;
      set.newton_local = newton_local;
      set.dense_lu = dense_lu;
      QDImplicitFault f(comm, cycle::PlanarProblem(comm, G, in, UnitG(B.op.get())), set);
      std::vector<real_t> U0(3 * n, 0.0);
      for (int o = 0; o < n; o++) { U0[3 * o + 1] = ph.psi_ss(ph.Vpl); }
      double t_event = -1.0, vmax = 0.0;
      RunSpec spec;
      spec.h0 = 1e3;
      spec.tmax = t_end;
      spec.post = [&](TS, QDImplicitFault &ff, double t, const std::vector<real_t> &all)
      {
         const double vm = VMax(ff, t, Owned(ff, all));
         vmax = std::max(vmax, vm);
         if (t_event < 0.0 && vm >= vth) { t_event = t; }
      };
      const RunOut r = Run(f, U0, spec);
      double dslip = 0.0, dpsi = 0.0;
      for (int j = 0; j < N && r.final_all.size() == size_t(3 * N); j++)
      {
         dslip = std::max(dslip, std::abs(static_cast<double>(r.final_all[3 * j] - ref[2 * j])));
         dpsi = std::max(dpsi, std::abs(static_cast<double>(r.final_all[3 * j + 1] -
                                                              ref[2 * j + 1])));
      }
      const double rel_t = std::abs(t_event - c1.t_event) / c1.t_event;
      const std::string w = "G5 " + ps + " ARK" + type;
      Say("    ARK" + std::string(type) + " (atol 1e-7): event at t = " + Fmt(t_event, 10) +
          " s (relative difference " + Fmt(rel_t, 3) + "), max |V| " + Fmt(vmax, 4) +
          " m/s; at the end max |slip difference| " + Fmt(dslip, 3) + " m of " +
          Fmt(slip_ref, 4) + " m, max |psi difference| " + Fmt(dpsi, 3) + "; " + Summary(r));
      Check(r.reason == TS_CONVERGED_TIME && AllFinite(r.final_all) && t_event > 0.0,
            w + ": reaches the end time through the event, finite");
      Check(rel_t < 0.01, w + ": event time within 1% of RK 5dp");
      Check(dslip < 0.01 * slip_ref, w + ": slip after the event within 1% of RK 5dp");
   }
}

// ---- G6: opening on the real G, all ranks vs one rank ------------------------------

struct OpenRun
{
   RunOut r;
   bool comp_ok = true;
   double min_sig = 0.0, min_g = 0.0, max_gs = 0.0, max_dpsi_open = 0.0;
   long open_steps = 0;
};

OpenRun OpeningRun(MPI_Comm comm, QDImplicitFault::Problem P, int N,
                   const std::vector<int> &pulled, double T, double sigma)
{
   P.extra_load = [=](real_t t, int j, real_t &dt, real_t &ds)
   {
      dt = 0.0;
      bool hit = false;
      for (int q : pulled) { hit = hit || (q == j); }
      ds = hit ? -1.2 * sigma * std::max(0.0, std::sin(M_PI * (3.0 * t / T - 1.0))) : 0.0;
   };
   OpenRun o;
   QDImplicitSettings set;
   set.newton_local = newton_local;
   set.dense_lu = dense_lu;
   QDImplicitFault f(comm, P, set);
   std::vector<real_t> U0(3 * f.NumOwned(), 0.0);
   for (int k = 0; k < f.NumOwned(); k++) { U0[3 * k + 1] = ph.psi_ss(ph.Vpl); }
   std::vector<real_t> prev_psi(N, 0.0);
   std::vector<char> prev_open(N, 0);
   QDImplicitFault::Counts c_prev;
   double t_prev = 0.0;
   RunSpec spec;
   spec.h0 = 1e3;
   spec.tmax = T;
   spec.max_steps = 20000;
   spec.post = [&](TS, QDImplicitFault &ff, double t, const std::vector<real_t> &all)
   {
      std::vector<real_t> V, tau, sig;
      std::vector<char> open;
      const std::vector<real_t> own = Owned(ff, all);
      ff.Evaluate(t, own.data(), V, tau, sig, open);
      const std::vector<real_t> s_all = GatherAll(ff.Comm(), sig.data(), ff.NumOwned());
      std::vector<real_t> openr(open.begin(), open.end());
      const std::vector<real_t> o_all = GatherAll(ff.Comm(), openr.data(), ff.NumOwned());
      for (int j = 0; j < N; j++)
      {
         const double g = all[3 * j + 2], sg = s_all[j];
         o.min_sig = std::min(o.min_sig, sg);
         o.min_g = std::min(o.min_g, g);
         o.max_gs = std::max(o.max_gs, std::abs(g * sg));
         if (o_all[j] == 1.0)
         {
            o.open_steps++;
            if (prev_open[j])
            {
               o.max_dpsi_open = std::max(o.max_dpsi_open,
                                          std::abs(static_cast<double>(all[3 * j + 1] -
                                                                       prev_psi[j])));
            }
         }
         prev_open[j] = (o_all[j] == 1.0);
         prev_psi[j] = all[3 * j + 1];
      }
      if (trace && ff.Comm() == MPI_COMM_WORLD)
      {
         const QDImplicitFault::Counts &c = ff.GetCounts();
         int nopen = 0;
         double smin = 1e300, vmax = 0.0;
         for (int j = 0; j < N; j++)
         {
            nopen += (o_all[j] == 1.0);
            smin = std::min(smin, static_cast<double>(s_all[j]));
         }
         for (real_t v : V) { vmax = std::max(vmax, std::abs(static_cast<double>(v))); }
         vmax = MaxAll(vmax, ff.Comm());
         std::ostringstream w;
         w << "    trace t " << std::setprecision(12) << t << " h " << std::setprecision(4)
           << (t - t_prev) << " open " << nopen << " min sigma " << smin << " Vmax " << vmax
           << " | factors +" << (c.factors - c_prev.factors) << " Newton its +"
           << (c.newton_its - c_prev.newton_its) << " fails +"
           << (c.newton_fails - c_prev.newton_fails) << " branch changes +"
           << (c.branch_changes - c_prev.branch_changes) << " clipped +"
           << (c.clipped_its - c_prev.clipped_its);
         Say(w.str());
         c_prev = c;
         t_prev = t;
      }
   };
   o.r = Run(f, U0, spec);
   o.comp_ok = o.min_sig >= -1e-6 * sigma && o.min_g >= 0.0 && o.max_gs <= 1e-6 * sigma * 1e-2;
   return o;
}

void TestG6(Box &B, int p)
{
   MPI_Comm comm = MPI_COMM_WORLD;
   const QDGreenOperator &G = *B.op->GetQDGreen();
   const int N = B.N;
   const std::string ps = "p = " + std::to_string(p);
   double s = 0.0;
   {
      PlanarFaultInputs in0 = Inputs(B, -1.0, [](long long, double) { return ph.Vpl; },
                                     [](long long, double) { return ph.sigma; });
      s = LoadSign(cycle::PlanarProblem(comm, G, in0, UnitG(B.op.get())));
   }
   const PlanarFaultInputs in = Inputs(B, s, [](long long, double) { return ph.Vpl; },
                                       [](long long, double) { return ph.sigma; });
   QDImplicitFault::Problem P = cycle::PlanarProblem(comm, G, in, UnitG(B.op.get()));
   // the two shallowest points
   std::vector<double> dloc(B.n);
   for (int o = 0; o < B.n; o++) { dloc[o] = B.depth(o); }
   const std::vector<real_t> d_all = GatherAll(comm, dloc.data(), B.n);
   std::vector<int> order(N);
   for (int j = 0; j < N; j++) { order[j] = j; }
   std::stable_sort(order.begin(), order.end(), [&](int x, int y) { return d_all[x] < d_all[y]; });
   const std::vector<int> pulled = {order[0], order[1]};
   const double T = 3e8;
   const OpenRun ow = OpeningRun(comm, P, N, pulled, T, ph.sigma);
   double gend = 0.0;
   for (int j = 0; j < N && ow.r.final_all.size() == size_t(3 * N); j++)
   {
      gend = std::max(gend, std::abs(static_cast<double>(ow.r.final_all[3 * j + 2])));
   }
   Say("  G6 " + ps + ": points " + std::to_string(pulled[0]) + ", " +
       std::to_string(pulled[1]) + " pulled to -0.2 sigma; " + Summary(ow.r) +
       "; point-steps open " + std::to_string(ow.open_steps) + ", min sigma " +
       Fmt(ow.min_sig, 3) + " Pa, min g " + Fmt(ow.min_g, 3) + " m, max |g sigma| " +
       Fmt(ow.max_gs, 3) + " N/m, max |psi change| over a step open at both ends " +
       Fmt(ow.max_dpsi_open, 3) + ", max |g| at the end " + Fmt(gend, 3) + " m");
   Check(ow.r.reason == TS_CONVERGED_TIME && AllFinite(ow.r.final_all) && ow.open_steps > 0,
         "G6 " + ps + ": reaches T with the pulled points opened on the way, finite");
   Check(ow.comp_ok, "G6 " + ps + ": sigma >= 0, g >= 0, g sigma = 0 at every accepted state");
   Check(ow.max_dpsi_open <= 1e-12, "G6 " + ps + ": psi frozen while open (to 1e-12)");
   Check(gend == 0.0, "G6 " + ps + ": closed with g = 0 at the end");

   // one rank: the gathered G and every gap column the run could need
   std::vector<std::vector<real_t>> ct(N), cs(N);
   for (int l = 0; l < N; l++) { P.gap_column(l, ct[l], cs[l]); }
   int steps = 0, nfinal = 0;
   double tend = 0.0;
   std::vector<real_t> fin;
   // gather the owned inputs (collective), then run on rank 0 alone
   std::vector<real_t> tp_own(B.n), s0_own(B.n), a_own(B.n), e_own(B.n), dc_own(B.n);
   for (int o = 0; o < B.n; o++)
   {
      tp_own[o] = in.tau_pre(2 * o + B.slot);
      s0_own[o] = in.sigma0(o);
      a_own[o] = in.a(o);
      e_own[o] = in.eta(o);
      dc_own[o] = in.Dc(o);
   }
   const std::vector<real_t> tp_all = GatherAll(comm, tp_own.data(), B.n);
   const std::vector<real_t> s0_all = GatherAll(comm, s0_own.data(), B.n);
   const std::vector<real_t> a_all = GatherAll(comm, a_own.data(), B.n);
   const std::vector<real_t> e_all = GatherAll(comm, e_own.data(), B.n);
   const std::vector<real_t> dc_all = GatherAll(comm, dc_own.data(), B.n);
   if (rank == 0)
   {
      QDImplicitFault::Problem Q = P;
      Q.lo = 0;
      Q.n = N;
      Q.tau_pre = tp_all;
      Q.sigma0 = s0_all;
      Q.a = a_all;
      Q.eta = e_all;
      Q.Dc = dc_all;
      Q.gap_column = [ct, cs](int l, std::vector<real_t> &t, std::vector<real_t> &s)
      {
         t = ct[l];
         s = cs[l];
      };
      const OpenRun os = OpeningRun(MPI_COMM_SELF, Q, N, pulled, T, ph.sigma);
      std::cout << "    1 rank: " << Summary(os.r) << std::endl;
      steps = static_cast<int>(os.r.steps);
      fin = os.r.final_all;
      nfinal = static_cast<int>(fin.size());
      tend = os.r.t_end;
   }
   MPI_Bcast(&steps, 1, MPI_INT, 0, comm);
   MPI_Bcast(&nfinal, 1, MPI_INT, 0, comm);
   MPI_Bcast(&tend, 1, MPI_DOUBLE, 0, comm);
   fin.resize(nfinal);
   MPI_Bcast(fin.data(), nfinal, MPITypeMap<real_t>::mpi_type, 0, comm);
   double maxdiff = 0.0;
   for (int i = 0; i < nfinal && i < static_cast<int>(ow.r.final_all.size()); i++)
   {
      maxdiff = std::max(maxdiff, std::abs(static_cast<double>(fin[i] - ow.r.final_all[i])));
   }
   Say("    max |state difference| " + Fmt(maxdiff, 3) + ", steps " +
       std::to_string(ow.r.steps) + " vs " + std::to_string(steps));
   Check(ow.r.steps == steps && ow.r.final_all == fin && ow.r.t_end == tend,
         "G6 " + ps + ": same steps and final state bit for bit on " +
         std::to_string(nranks) + " rank(s) and on 1 rank");
}

// ---- G7: Jacobian and LU time at N = 2,638 -----------------------------------------

void TestG7()
{
   Say("G7: one Jacobian (assembly + dense LU) and one solve at N = 2,638, rank 0 alone");
   if (rank == 0)
   {
      const int N = 2638;
      for (const int m : {0, 50})
      {
         QDImplicitFault::Problem P;
         P.N = N;
         P.Gt.assign(size_t(N) * N, 0.0);
         P.Gs.assign(size_t(N) * N, 0.0);
         const double k = 1.4e7;
         for (int j = 0; j < N; j++)
         {
            for (int l = 0; l < N; l++)
            {
               const double d = std::abs(j - l);
               P.Gt[size_t(j) * N + l] = (j == l) ? k : -0.25 * k / (d * d);
               P.Gs[size_t(j) * N + l] = (j == l) ? 0.0 : 0.05 * k * ((l > j) ? 1.0 : -1.0) / (d * d);
            }
         }
         P.TDt.assign(N, -1e-2);
         P.TDs.assign(N, 0.0);
         P.lo = 0;
         P.n = N;
         P.tau_pre.resize(N);
         P.sigma0.resize(N);
         for (int j = 0; j < N; j++)
         {
            const bool tens = m > 0 && j % (N / m) == 0 && j / (N / m) < m;
            P.sigma0[j] = tens ? -2e5 : 50e6;
            const double sg = std::max(static_cast<double>(P.sigma0[j]), 0.0);
            P.tau_pre[j] = -(sg * RSFrictionCoefficientPsi(1e-9, ph.psi_ss(1e-9), 0.01, ph.V0) +
                             ph.eta() * 1e-9);
         }
         P.a.assign(N, 0.01);
         P.eta.assign(N, ph.eta());
         P.Dc.assign(N, 0.01);
         P.b = ph.b;
         P.f0 = ph.f0;
         P.V0 = ph.V0;
         const double gam = 2.0 * k;
         P.gap_column = [=](int l, std::vector<real_t> &t, std::vector<real_t> &s)
         {
            t.assign(N, 0.0);
            s.assign(N, 0.0);
            for (int j = 0; j < N; j++)
            {
               const double d = std::abs(j - l);
               t[j] = (j == l) ? 0.1 * k : 0.02 * k / (d * d);
               s[j] = (j == l) ? gam : -0.2 * gam / (d * d);
            }
         };
         QDImplicitFault f(MPI_COMM_SELF, P, QDImplicitSettings());
         std::vector<real_t> U(3 * N, 0.0), Ud(3 * N, 0.0);
         for (int j = 0; j < N; j++) { U[3 * j + 1] = ph.psi_ss(1e-9); }
         const long opened = f.ConsistentStart(0.0, U.data());
         const double t0 = MPI_Wtime();
         f.Jacobian(0.0, U.data(), Ud.data(), 1e-3);
         const double t1 = MPI_Wtime();
         std::vector<real_t> b(3 * N), y(3 * N);
         for (int i = 0; i < 3 * N; i++)
         {
            b[i] = (i % 3 == 2 && !f.HasColumn(i / 3)) ? 0.0 : 1e-6 * std::sin(0.1 * i);
         }
         f.ApplyInverse(b.data(), y.data());
         const double t2 = MPI_Wtime();
         const int M = N + static_cast<int>(opened);
         const double mib = (2.0 * double(N) * N + double(M) * M) * sizeof(real_t) / 1048576.0;
         std::cout << "    " << opened << " open (matrix " << M << " x " << M << "): Jacobian "
                   << std::setprecision(3) << (t1 - t0) << " s (assembly + LU), solve "
                   << (t2 - t1) * 1e3 << " ms, factors " << f.GetCounts().factors
                   << "; G and LU " << std::setprecision(4) << mib << " MiB per rank"
                   << std::endl;
      }
   }
   MPI_Barrier(MPI_COMM_WORLD);
}

} // namespace

int main(int argc, char **argv)
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   MFEMInitializePetsc(&argc, &argv, nullptr, nullptr);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &nranks);
   Say("test_qd_implicit_green on " + std::to_string(nranks) + " rank(s), PETSc " +
       std::to_string(PETSC_VERSION_MAJOR) + "." + std::to_string(PETSC_VERSION_MINOR) + "." +
       std::to_string(PETSC_VERSION_SUBMINOR));
   MFEM_VERIFY(PetscOptionsSetValue(nullptr, "-ts_adapt_wnormtype", "infinity")
               == PETSC_SUCCESS, "options");
   char buf[64] = "";
   PetscBool has = PETSC_FALSE, tr = PETSC_FALSE;
   MFEM_VERIFY(PetscOptionsGetString(nullptr, nullptr, "-only", buf, sizeof(buf), &has)
               == PETSC_SUCCESS, "-only");
   if (has) { only = buf; }
   MFEM_VERIFY(PetscOptionsGetBool(nullptr, nullptr, "-trace", &tr, nullptr) == PETSC_SUCCESS,
               "-trace");
   trace = (tr == PETSC_TRUE);
   PetscBool nl = newton_local ? PETSC_TRUE : PETSC_FALSE, has_nl = PETSC_FALSE;
   MFEM_VERIFY(PetscOptionsGetBool(nullptr, nullptr, "-qd_newton_local", &nl, &has_nl)
               == PETSC_SUCCESS, "-qd_newton_local");
   newton_local = (nl == PETSC_TRUE);
   char dl[32] = "serial";
   PetscBool has_dl = PETSC_FALSE;
   MFEM_VERIFY(PetscOptionsGetString(nullptr, nullptr, "-qd_dense_lu", dl, sizeof(dl), &has_dl)
               == PETSC_SUCCESS, "-qd_dense_lu");
   dense_lu = dl;
   if (has_dl) { Say("dense_lu " + dense_lu + " in every QDImplicitFault of this run (-qd_dense_lu)"); }
   if (has_nl) { Say(std::string("newton_local ") + (newton_local ? "on" : "off") +
                     " in G5, G6 (-qd_newton_local)"); }
   for (int p : {1, 2})
   {
      Say("\n=== the miniature BP3 layout, p = " + std::to_string(p) + " ===");
      std::unique_ptr<Box> B = MakeBox(MPI_COMM_WORLD, p);
      if (Wanted('1') || Wanted('2')) { TestG1G2(*B, p); }
      if (Wanted('3') || Wanted('4')) { TestG3G4(*B, p); }
      if (Wanted('5')) { TestG5(*B, p); }
      if (Wanted('6')) { TestG6(*B, p); }
      if (Wanted('8')) { TestG8(*B, p); }
   }
   if (Wanted('7')) { TestG7(); }
   Say(fails ? "test_qd_implicit_green: " + std::to_string(fails) + " FAILED" :
       "test_qd_implicit_green: all PASSED");
   MFEMFinalizePetsc();
   Mpi::Finalize();
   return fails ? 1 : 0;
}

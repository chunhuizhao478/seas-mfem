// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_2d5_width.cpp — the 2.5D finite-width term of a map-view cycle
// ([cycle.finite_width], cycle/finite_width.hpp; Chunhui Zhao, 2026-10-30:
// "yes send this paper and equation to 2D thread and start implementation
// and tests"), the Level 1 checks of the derivation note's verification plan
// (/mnt/project-files/2d-inplane/nonplanar-literature/inplane-2p5d-derivation-rev8.pdf,
// Sections 6-7):
//
//   rho u'' = div sigma(u) - beta (u - u_pl),   beta = 4 c^2 mu / W^2 = mu / l^2
//
// W = 15 km, c = (1 - nu) pi / 2 ("wa2020", nu = 1/4), so l = W / (2 c) =
// 6.37 km, on lambda = mu = 32.04 GPa, rho = 2670 kg/m^3 (the map-view
// deck's material).  Meshes: grids of squares, each split into 4 triangles
// through its centre (test_2d_nitsche.cpp's CrossGrid), the fault on a grid
// line (attribute 3), map-view frames (ref_normal (0, -1), up (0, 0, 1): the
// in-plane slip in the strike slot).
//
//   U1  manufactured static solution: u* = A (cos kx cos ky, sin kx sin ky),
//       k = 2 / l on [-l, l]^2, div u* = 0, so div sigma(u*) = -2 mu k^2 u*
//       and u_pl = (1 + 2 k^2 l^2) u* makes u* exact; Dirichlet u* on the
//       walls.  The QD operator (ElasticityDomainOperator, MUMPS) with
//       SetFiniteWidth(1 / l^2, U_pl1): L2 order on the finest pair above
//       p + 0.5 (p = 1, 2), without and with a zero-slip fault across the box.
//       Linearity in time and the time offset of the dynamic regime: the load
//       (t_off + t) M_beta U_pl1 with t_off = 0.5, t = 1.5 and the walls at
//       2 u* gives the t = 1 solution times 2 (1e-12).
//   U2  long fault: the fault spans the box [-2 l, 2 l] x [-6 l, 6 l] on
//       y = 0, uniform strike slip D, u_pl = 0, the walls carry the exact
//       field u = sgn (D / 2) s e^{-|y| / l} (sgn = +1 on the plus side):
//       the profile (L2 order above p + 0.5) and the plus-side traction
//       tau_s = mu D / (2 l) = c mu D / W (the note's eq. 6; relative RMS
//       error falling with h), the normal traction 0, the dip slot exactly 0.
//       Control: the same walls without the term miss the profile by more
//       than 10 times.
//   U3  stress drop -> slip: a fault of half-length a = 10 l in a clamped
//       box, the strike-strike traction matrix G from unit slips (MUMPS
//       solves), G delta = dtau 1: the slip in |x| < l is W dtau / (c mu) =
//       2 l dtau / mu to 1%; the plane-strain value (no term) printed.
//       One rank only.
//   U4  the dynamic operator (NitscheWaveOperator::SetFiniteWidth):
//       (a) (K + M_beta) x - K x = (beta / rho) M x for a random x (1e-12);
//       (b) a uniform displacement on a traction-free box oscillates at the
//           cutoff, U^n = U^0 cos(n theta), cos theta = 1 - beta dt^2 / (2
//           rho), over one period (1e-9);
//       (c) dispersion: plane P and S waves with omega^2 = c^2 k^2 + beta /
//           rho (wavelength 120 km, l_k = 1 / k) and, below the cutoff,
//           evanescent P and S waves d e^{-kappa n.x} cos(omega t) at omega =
//           sqrt(beta / rho) / 2, kappa^2 = (beta / rho - omega^2) / c^2,
//           exact Dirichlet data: L2 order above p + 0.5 (p = 2); stepped
//           without the term they miss by more than 10 times the finest
//           error;
//       (d) rule C3: Kelvin-Voigt (whole domain) on K alone, M^-1 C V = 0
//           for a uniform V (1e-10 of (beta / rho) tau_kv |V|);
//       (e) rule C6: the largest eigenvalue of M^-1 (K + M_beta) is that of
//           M^-1 K plus beta / rho (1e-6), and with beta / rho = 2 omega_max^2
//           a step 0.95 x 2 / omega_max blows up while the step cut by
//           1 / sqrt(1 + beta dt^2 / (4 rho)) (DynPhase) stays bounded.
//   U5  the regime switch on the map-view deck (default config/fd_cycle/
//       mapview_rs_tensor_fd.toml, 500 m, p = 1) with [cycle.finite_width]
//       W 15 km, wa2020, for each substrate ("loading_plane",
//       "locked_reference"): a QD balance at t_s = 3e7 s with a tapered strike
//       slip 0.02 m cos^2(pi x / 2 L_f) and psi0, handed to the dynamic regime
//       (faces "off": the walls clamped in both regimes):
//         C1  (K_W - K_W0) U_pl1 = M_beta U_pl1 (the f_hat of SetFiniteWidth,
//             1e-12); K_W u_s = K_QD u_s off the fault elements (1e-12);
//             with the slip held, M^-1 (F - K u_s + G d_s) = 0 at every DOF
//             (the QD balance is where the damped dynamic operator comes to
//             rest; 1e-8 a_K, or 10 x the residual without the term);
//             control: W0 leaves more than 1e3 x that;
//         C2  the dynamic load F_W(tau) = F_QD(t_s + tau) at tau = 0 and after
//             n steps (1e-13);
//         A0  the first acceleration off the fault elements below 1e-8 a_K
//             (a_K = max |M^-1 K u_s|); control: the dynamic regime without
//             the term (W0) is off there by about (beta / rho)(u_s - u_pl)
//             (above half of it);
//         C6  dt_cfl = dt_cfl(W0) / sqrt(1 + beta dt_cfl(W0)^2 / (4 rho));
//         hand-back: after n steps, ToQD + Reenter at t_s + tau: the hand-back
//             table's largest |tau_qd - tau_dyn| and |sigma_qd - sigma_dyn|
//             at most twice the same pair's without the term (+ 1e-9 of max
//             |tau|); control: W0 handed back to the QD regime with the term
//             jumps by more;
//         all  the deck's faces = "all" (the frozen wall load): A0 off the
//             fault elements below 1e-8 a_K.
//   U7  normal stiffness (the note's Step 7b): a uniform normal jump D_n,
//       solved on the plus half (the fault slip has no opening slot): the
//       strip [-2 l_n, 2 l_n] x [0, 6 l_n], every wall at u_y = (D_n / 2)
//       e^{-y / l_n}, l_n = sqrt(M / beta), M = lambda + 2 mu: the profile
//       (L2 order above p + 0.5) and the face stress sigma_yy(0+) = -(D_n /
//       2) sqrt(M beta) (error falling with h, below 1% at h = l_n / 8 for
//       p = 2); the ratio to the 3D mode I stiffness printed (isotropic beta,
//       Chunhui Zhao, 2026-10-30: "we will start with beta = beta (normal
//       direction)").
//   U8  the "loading_plane" substrate (the note rev8, Section 5, items 2
//       and 3) on U2's strip, U_pl1 from cycle::MakePlateSubstrate and
//       cycle::InterpolateSubstrate as QDPhase builds it, Vp = 1e-9 m/s, t =
//       1 yr: a locked fault (walls at u_pl (1 - e^{-|y| / l})) is loaded at
//       tau_s = -c mu Vp t / W (U2's sign; error falling with h, below 2% at
//       h = l / 8, p = 1, 2); steady creep at Vp (slip Vp t, walls at u_pl)
//       keeps u = u_pl and tau_s = 0 (1e-10, 1e-8); control: without the
//       substrate load the walls do not load the fault (error above 0.25,
//       the walls at x = +-2 l reach in by l).
//   R   refusals (--refuse N, each must exit non-zero with its message):
//       1 beta <= 0, 2 U_pl1 of the wrong size, 3 a 3-D mesh, 4
//       NitscheWaveOperator::SetFiniteWidth after the dense stiffness read K,
//       5 DynPhase with the term but no substrate.  The deck refusals run
//       through seas_fd_cycle_driver in the make target.
//
// Usage (from miniapps/seas; PETSc needed for U5 and refusal 5):
//   mpirun -np 1 ./seas_test_2d5_width [--only U1,U2,...,U8] [--deck ...]
//          [--dyn-steps 20] [--refuse N]
//   make -f tests/2d_inplane/2d_inplane.mk test-2d5-width
// Exit code 0 = PASSED.
//
// Scope (Chunhui Zhao, 2026-10-30).

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../solver/nitsche_cfl.hpp"
#include "../../cycle/scaled_mesh.hpp"
#include "../../cycle/cycle_config.hpp"
#include "../../cycle/fault_fields.hpp"
#include "../../cycle/prestress.hpp"
#include "../../cycle/qd_phase.hpp"
#include "../../cycle/dyn_phase.hpp"
#include "../../cycle/state_bridge.hpp"
#include "../../cycle/finite_width.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const real_t kLam = 32038120320.0, kMu = 32038120320.0, kRho = 2670.0;
static const real_t kW = 15.0e3;
static real_t Nu() { return kLam / (2.0 * (kLam + kMu)); }
static real_t CWa() { return (1.0 - Nu()) * M_PI / 2.0; }
static real_t Ell() { return kW / (2.0 * CWa()); }          ///< l = W / (2 c)
static real_t BetaOverMu() { return 1.0 / (Ell() * Ell()); }
static real_t Cs() { return std::sqrt(kMu / kRho); }
static real_t Cp() { return std::sqrt((kLam + 2.0 * kMu) / kRho); }
static real_t CUse(int p) { return 0.5 * NitscheStableConstant(2, p); }

static real_t RSum(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, c); return v; }
static real_t RMax(real_t v, MPI_Comm c) { MPI_Allreduce(MPI_IN_PLACE, &v, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, c); return v; }
static real_t Dot(const Vector &a, const Vector &b, MPI_Comm c) { return RSum(a * b, c); }
static real_t Norm(const Vector &a, MPI_Comm c) { return std::sqrt(Dot(a, a, c)); }
static real_t MaxAbs(const Vector &a, MPI_Comm c) { return RMax(a.Normlinf(), c); }
static real_t MaxDiff(const Vector &a, const Vector &b, MPI_Comm c)
{
   MFEM_VERIFY(a.Size() == b.Size(), "MaxDiff: sizes " << a.Size() << " vs " << b.Size());
   real_t d = 0.0;
   for (int i = 0; i < a.Size(); i++) { d = std::max(d, std::abs(a(i) - b(i))); }
   return RMax(d, c);
}

using ExactFn = std::function<void(const Vector &, real_t, Vector &)>;
using Op = ElasticityDomainOperator<ParMesh>;

// ------------------------------------------------------------------ meshes
// test_2d_nitsche.cpp's CrossGrid: [x0, x0 + nx hx] x [y0, y0 + ny hy], each
// cell split into 4 triangles through its centre; outer edges get
// bdr_attr(midpoint), interior grid edges with fault(midpoint) attribute 3.
using AttrFn = std::function<int(real_t, real_t)>;
using FaultFn = std::function<bool(real_t, real_t)>;

static Mesh *CrossGrid(int nx, int ny, real_t x0, real_t y0, real_t hx, real_t hy,
                       const AttrFn &bdr_attr, const FaultFn &fault = nullptr)
{
   const int ngv = (nx + 1) * (ny + 1);
   Mesh *m = new Mesh(2, ngv + nx * ny, 4 * nx * ny, 2 * (nx + ny), 2);
   auto gv = [&](int i, int j) { return i + (nx + 1) * j; };
   auto cv = [&](int i, int j) { return ngv + i + nx * j; };
   for (int j = 0; j <= ny; j++)
      for (int i = 0; i <= nx; i++) { m->AddVertex(x0 + i * hx, y0 + j * hy); }
   for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++) { m->AddVertex(x0 + (i + 0.5) * hx, y0 + (j + 0.5) * hy); }
   for (int j = 0; j < ny; j++)
      for (int i = 0; i < nx; i++)
      {
         const int c = cv(i, j);
         m->AddTriangle(gv(i, j), gv(i + 1, j), c, 1);
         m->AddTriangle(gv(i + 1, j), gv(i + 1, j + 1), c, 1);
         m->AddTriangle(gv(i + 1, j + 1), gv(i, j + 1), c, 1);
         m->AddTriangle(gv(i, j + 1), gv(i, j), c, 1);
      }
   for (int i = 0; i < nx; i++)
   {
      const real_t xm = x0 + (i + 0.5) * hx;
      m->AddBdrSegment(gv(i, 0), gv(i + 1, 0), bdr_attr(xm, y0));
      m->AddBdrSegment(gv(i + 1, ny), gv(i, ny), bdr_attr(xm, y0 + ny * hy));
   }
   for (int j = 0; j < ny; j++)
   {
      const real_t ym = y0 + (j + 0.5) * hy;
      m->AddBdrSegment(gv(nx, j), gv(nx, j + 1), bdr_attr(x0 + nx * hx, ym));
      m->AddBdrSegment(gv(0, j + 1), gv(0, j), bdr_attr(x0, ym));
   }
   if (fault)
   {
      for (int j = 1; j < ny; j++)
         for (int i = 0; i < nx; i++)
         {
            if (fault(x0 + (i + 0.5) * hx, y0 + j * hy)) { m->AddBdrSegment(gv(i, j), gv(i + 1, j), 3); }
         }
   }
   m->FinalizeTopology();
   m->Finalize(false, true);
   m->SetAttributes();
   return m;
}

static AttrFn Walls(int a) { return [a](real_t, real_t) { return a; }; }

/// Map-view frames on the fault (attribute 3): ref_normal (0, -1), up e_z.
static DomainConfig MapViewFrames()
{
   auto fc = std::make_shared<FaultFrameConfig>();
   const real_t ref[3] = {0.0, -1.0, 0.0}, up[3] = {0.0, 0.0, 1.0};
   for (int d = 0; d < 3; d++)
   {
      fc->orientation.ref_normal[d] = ref[d];
      fc->orientation.up[d] = up[d];
      fc->up[d] = up[d];
   }
   fc->orientation.fault_attrs = {3};
   DomainConfig dc;
   dc.fault_frames.config = fc;
   return dc;
}

/// The plus side's unit normal p = -n and the slip direction s of every
/// fault frame; false unless all frames agree (a straight fault).
static bool UniformFrame(Op &op, real_t p[2], real_t s[2])
{
   auto F = op.FrameField();
   if (!F) { return false; }
   bool ok = true, set = false;
   for (int f = 0; f < F->NumFaultFaces(); f++)
   {
      const fault_frame::Frame &fr = F->Frame(f);
      const real_t pp[2] = {-fr.n[0], -fr.n[1]}, ss[2] = {fr.s[0], fr.s[1]};
      if (!set) { p[0] = pp[0]; p[1] = pp[1]; s[0] = ss[0]; s[1] = ss[1]; set = true; }
      for (int i = 0; i < 2; i++)
      {
         if (std::abs(pp[i] - p[i]) > 1e-12 || std::abs(ss[i] - s[i]) > 1e-12) { ok = false; }
      }
      if (std::abs(fr.n[2]) > 1e-12 || std::abs(fr.s[2]) > 1e-12) { ok = false; }
   }
   // the same on every rank
   real_t loc[4] = {set ? p[0] : -2.0, set ? p[1] : -2.0, set ? s[0] : -2.0, set ? s[1] : -2.0};
   real_t mx[4], mn[4];
   MPI_Allreduce(loc, mx, 4, MPITypeMap<real_t>::mpi_type, MPI_MAX, op.GetMesh().GetComm());
   for (int i = 0; i < 4; i++) { if (!set) { loc[i] = 2.0; } }
   MPI_Allreduce(loc, mn, 4, MPITypeMap<real_t>::mpi_type, MPI_MIN, op.GetMesh().GetComm());
   for (int i = 0; i < 4; i++) { if (mx[i] != mn[i]) { ok = false; } }
   p[0] = mx[0]; p[1] = mx[1]; s[0] = mx[2]; s[1] = mx[3];
   int okl = ok ? 1 : 0;
   MPI_Allreduce(MPI_IN_PLACE, &okl, 1, MPI_INT, MPI_MIN, op.GetMesh().GetComm());
   return okl == 1;
}

static real_t RelL2(ParGridFunction &u, VectorCoefficient &ex)
{
   ParGridFunction zero(u.ParFESpace());
   zero = 0.0;
   return u.ComputeL2Error(ex) / zero.ComputeL2Error(ex);
}

// ====================================================================== U1
struct Manufactured
{
   real_t A = 0.1, k = 2.0 / Ell();
   void U(const Vector &x, Vector &u) const
   {
      u.SetSize(2);
      u(0) = A * std::cos(k * x(0)) * std::cos(k * x(1));
      u(1) = A * std::sin(k * x(0)) * std::sin(k * x(1));
   }
   real_t PlFactor() const { return 1.0 + 2.0 * k * k / BetaOverMu(); }   ///< u_pl = this u*
};

/// The static solve of U1 at time t (walls (t_off + t) u*, load (t_off + t)
/// M_beta U_pl1); returns the relative L2 error and the field.
static real_t RunManufactured(MPI_Comm comm, int p, int n, bool fault, real_t t_off, real_t t,
                              Vector *U = nullptr)
{
   const real_t L = Ell(), h = 2.0 * L / n;
   std::unique_ptr<Mesh> m(CrossGrid(n, n, -L, -L, h, h, Walls(5),
                                     fault ? FaultFn([](real_t, real_t y) { return std::abs(y) < 1e-6; })
                                           : FaultFn()));
   ParMesh pm(comm, *m);
   m.reset();
   const Manufactured ms;
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = [ms, t_off](const Vector &x, real_t tt, Vector &u)
   {
      ms.U(x, u);
      u *= (t_off + tt);
   };
   LinearElastic model(kLam, kMu);
   Op op(pm, p, model, 0.0, 1e4, 1e4, bc, DGMethod::IP, SolverType::MUMPS,
         fault ? MapViewFrames() : DomainConfig());
   ParFiniteElementSpace &fes = op.GetFESpace();
   VectorFunctionCoefficient upl(2, [ms](const Vector &x, Vector &u) { ms.U(x, u); u *= ms.PlFactor(); });
   ParGridFunction gpl(&fes);
   gpl.ProjectCoefficient(upl);
   op.SetFiniteWidth(BetaOverMu(), gpl, t_off != 0.0 ? std::make_shared<const real_t>(t_off) : nullptr);
   Vector slip(2 * op.GetNumFaultDOFs());
   slip = 0.0;
   ParGridFunction u(&fes);
   op.Solve(t, slip, u);
   if (U) { *U = u; }
   VectorFunctionCoefficient ex(2, [ms, t_off, t](const Vector &x, Vector &v) { ms.U(x, v); v *= (t_off + t); });
   return RelL2(u, ex);
}

static void TestManufactured(MPI_Comm comm)
{
   if (g_rank == 0)
   {
      std::cout << "\n=== U1 manufactured static solution (l = " << Ell() << " m, k = 2 / l, u_pl = "
                << Manufactured().PlFactor() << " u*) ===\n";
   }
   for (const bool fault : {false, true})
      for (int p = 1; p <= 2; p++)
      {
         const int n0 = (p == 1) ? 8 : 4;
         real_t e[3];
         for (int r = 0; r < 3; r++) { e[r] = RunManufactured(comm, p, n0 << r, fault, 0.0, 1.0); }
         const real_t o1 = std::log2(e[0] / e[1]), o2 = std::log2(e[1] / e[2]);
         if (g_rank == 0)
         {
            std::cout << std::scientific << std::setprecision(3) << "  p = " << p
                      << (fault ? ", zero-slip fault across" : ", no fault") << ": rel L2 error "
                      << e[0] << ", " << e[1] << ", " << e[2] << " on " << n0 << ", " << 2 * n0
                      << ", " << 4 * n0 << " cells per side; order " << std::fixed
                      << std::setprecision(2) << o1 << ", " << o2 << "\n";
         }
         TEST_ASSERT(o2 > p + 0.5, "U1 p = " << p << (fault ? " (fault)" : "")
                     << ": L2 order on the finest pair " << o2 << " > p + 0.5");
      }
   // linearity in t and the time offset
   Vector U1, U2, U3;
   RunManufactured(comm, 2, 8, true, 0.0, 1.0, &U1);
   RunManufactured(comm, 2, 8, true, 0.5, 1.5, &U2);
   RunManufactured(comm, 2, 8, true, 0.0, 2.0, &U3);
   Vector U1x2(U1);
   U1x2 *= 2.0;
   const real_t s = MaxAbs(U1x2, comm);
   const real_t d2 = MaxDiff(U2, U1x2, comm) / s, d3 = MaxDiff(U3, U1x2, comm) / s;
   TEST_ASSERT(d2 <= 1e-12 && d3 <= 1e-12, "U1 load (t_off + t) f_hat: t_off 0.5, t 1.5 and t 2 give "
               "2 x the t = 1 field (" << std::scientific << std::setprecision(2) << d2 << ", "
               << d3 << " of max |u|)");
}

// ====================================================================== U2
struct LongFault
{
   real_t D = 1.0;
   real_t p[2] = {0.0, 1.0}, s[2] = {1.0, 0.0};   ///< plus-side normal, slip direction
   void U(const Vector &x, Vector &u) const
   {
      const real_t d = x(0) * p[0] + x(1) * p[1];   // signed distance into the plus side
      const real_t side = (d > 0.0) ? 1.0 : (d < 0.0) ? -1.0 : 0.0;
      const real_t a = side * 0.5 * D * std::exp(-std::abs(d) / Ell());
      u.SetSize(2);
      u(0) = a * s[0];
      u(1) = a * s[1];
   }
};

struct U2Result { real_t e_u = 0, e_tau = 0, e_sn = 0, dip = 0; long dofs = 0; bool frame_ok = false; };

static U2Result RunLongFault(MPI_Comm comm, int p, int nper, bool with_beta)
{
   const real_t l = Ell(), h = l / nper, Lx = 2.0 * l, Ly = 6.0 * l;
   const int nx = (int) std::lround(2.0 * Lx / h), ny = (int) std::lround(2.0 * Ly / h);
   std::unique_ptr<Mesh> m(CrossGrid(nx, ny, -Lx, -Ly, h, h, Walls(5),
                                     [](real_t, real_t y) { return std::abs(y) < 1e-6; }));
   ParMesh pm(comm, *m);
   m.reset();
   auto lf = std::make_shared<LongFault>();
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = [lf](const Vector &x, real_t, Vector &u) { lf->U(x, u); };
   LinearElastic model(kLam, kMu);
   Op op(pm, p, model, 0.0, 1e4, 1e4, bc, DGMethod::IP, SolverType::MUMPS, MapViewFrames());
   U2Result r;
   r.frame_ok = UniformFrame(op, lf->p, lf->s);
   ParFiniteElementSpace &fes = op.GetFESpace();
   if (with_beta)
   {
      Vector zero(fes.GetVSize());
      zero = 0.0;
      op.SetFiniteWidth(BetaOverMu(), zero, nullptr);
   }
   const int nd = op.GetNumFaultDOFs();
   Vector slip(2 * nd);
   for (int k = 0; k < nd; k++) { slip(2 * k) = 0.0; slip(2 * k + 1) = lf->D; }
   ParGridFunction u(&fes);
   op.Solve(0.0, slip, u);
   Vector trac, sn, xyz;
   op.ComputeTraction(u, slip, trac, &sn);
   op.GetFaultDOFCoords3D(xyz);
   VectorFunctionCoefficient ex(2, [lf](const Vector &x, Vector &v) { lf->U(x, v); });
   r.e_u = RelL2(u, ex);
   const real_t tref = kMu * lf->D / (2.0 * l);
   real_t s2[4] = {0.0, 0.0, 0.0, 0.0};   // tau err^2, sn^2, count, dip max
   for (int k = 0; k < nd; k++)
   {
      s2[3] = std::max(s2[3], std::abs(trac(2 * k)));
      if (std::abs(xyz(3 * k)) > Lx - h) { continue; }   // the faces at the walls
      s2[0] += (trac(2 * k + 1) - tref) * (trac(2 * k + 1) - tref);
      s2[1] += sn(k) * sn(k);
      s2[2] += 1.0;
   }
   MPI_Allreduce(MPI_IN_PLACE, s2, 3, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);
   r.e_tau = std::sqrt(s2[0] / s2[2]) / tref;
   r.e_sn = std::sqrt(s2[1] / s2[2]) / tref;
   r.dip = RMax(s2[3], comm) / tref;
   r.dofs = (long) RSum((real_t) fes.GlobalTrueVSize(), MPI_COMM_SELF);
   return r;
}

static void TestLongFault(MPI_Comm comm)
{
   if (g_rank == 0)
   {
      std::cout << "\n=== U2 long fault: profile e^{-|y| / l}, tau_s = mu D / (2 l) = c mu D / W"
                << " (l = " << Ell() << " m) ===\n";
   }
   for (int p = 1; p <= 2; p++)
   {
      const int npers[3] = {2, 4, 8};
      U2Result r[3];
      for (int i = 0; i < 3; i++) { r[i] = RunLongFault(comm, p, npers[i], true); }
      for (int i = 0; i < 3; i++)
      {
         if (g_rank == 0)
         {
            std::cout << std::scientific << std::setprecision(3) << "  p = " << p << ", h = l / "
                      << npers[i] << " (" << r[i].dofs << " DOFs): rel L2 error of u " << r[i].e_u
                      << ", RMS |tau_s - mu D / 2 l| / (mu D / 2 l) " << r[i].e_tau
                      << ", RMS |sigma_n| / (mu D / 2 l) " << r[i].e_sn << ", max |dip slot| "
                      << r[i].dip << "\n";
         }
      }
      const real_t o = std::log2(r[1].e_u / r[2].e_u);
      TEST_ASSERT(r[0].frame_ok && r[1].frame_ok && r[2].frame_ok,
                  "U2 p = " << p << ": one straight frame (n along y, s along x)");
      TEST_ASSERT(o > p + 0.5, "U2 p = " << p << ": profile L2 order on the finest pair " << o
                  << " > p + 0.5");
      TEST_ASSERT(r[2].e_tau < r[1].e_tau && r[1].e_tau < r[0].e_tau && r[2].e_tau < 0.02,
                  "U2 p = " << p << ": tau_s = c mu D / W, the error falls with h and is "
                  << r[2].e_tau << " (< 2%) at h = l / 8");
      TEST_ASSERT(r[2].e_sn < 0.02, "U2 p = " << p << ": normal traction " << r[2].e_sn
                  << " of mu D / 2 l (< 2%) at h = l / 8");
      TEST_ASSERT(r[0].dip == 0.0 && r[1].dip == 0.0 && r[2].dip == 0.0,
                  "U2 p = " << p << ": the dip slot is exactly 0");
   }
   const U2Result c = RunLongFault(comm, 2, 4, false), b = RunLongFault(comm, 2, 4, true);
   if (g_rank == 0)
   {
      std::cout << std::scientific << std::setprecision(3) << "  control without the term, p = 2,"
                << " h = l / 4: rel L2 error " << c.e_u << ", tau_s error " << c.e_tau << "\n";
   }
   TEST_ASSERT(c.e_u > 10.0 * b.e_u, "U2 control: without the term the exact walls miss the "
               "profile by " << c.e_u / b.e_u << " x the error with it (> 10)");
}

// ====================================================================== U8
// The "loading_plane" substrate on U2's strip (the note rev8, Section 5,
// items 2 and 3): u_pl = sgn(y) (Vp t / 2) e_x from cycle::MakePlateSubstrate
// and cycle::InterpolateSubstrate, as QDPhase builds it.
//   locked  zero slip, walls at u_pl (1 - e^{-|y| / l}): the locked fault is
//           loaded at tau_s = -c mu Vp t / W (U2's sign: slip D gives +c mu D
//           / W, and the locked fault carries the stress of slip -Vp t);
//   creep   slip Vp t, walls at u_pl: u = u_pl, tau_s = 0 (steady creep at
//           the plate rate is stress-free);
//   control the locked case with the substrate load left out (U_pl1 = 0):
//           the walls y = +-6 l reach the fault only through e^{-6}, tau_s
//           misses (the walls x = +-2 l reach in by l).
struct U8Result { real_t e_tau = 0, tau_max = 0, e_u = 0; bool frame_ok = false; };

static U8Result RunPlateSubstrate(MPI_Comm comm, int p, int nper, int mode, real_t T, real_t Vp)
{
   const real_t l = Ell(), h = l / nper, Lx = 2.0 * l, Ly = 6.0 * l;
   const int nx = (int) std::lround(2.0 * Lx / h), ny = (int) std::lround(2.0 * Ly / h);
   std::unique_ptr<Mesh> m(CrossGrid(nx, ny, -Lx, -Ly, h, h, Walls(5),
                                     [](real_t, real_t y) { return std::abs(y) < 1e-6; }));
   ParMesh pm(comm, *m);
   m.reset();
   const bool creep = (mode == 1);
   auto exact = [l, Vp, creep](const Vector &x, real_t t, Vector &u)
   {
      u.SetSize(2);
      const real_t side = (x(1) > 0.0) ? 1.0 : (x(1) < 0.0) ? -1.0 : 0.0;
      u(0) = side * 0.5 * Vp * t * (creep ? 1.0 : 1.0 - std::exp(-std::abs(x(1)) / l));
      u(1) = 0.0;
   };
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = exact;
   LinearElastic model(kLam, kMu);
   Op op(pm, p, model, 0.0, 1e4, 1e4, bc, DGMethod::IP, SolverType::MUMPS, MapViewFrames());
   U8Result r;
   real_t pn[2] = {0.0, 1.0}, sd[2] = {1.0, 0.0};
   r.frame_ok = UniformFrame(op, pn, sd);
   ParFiniteElementSpace &fes = op.GetFESpace();
   Vector U_pl1(fes.GetVSize());
   U_pl1 = 0.0;
   if (mode != 2)
   {
      cycle::LoadingConfig L;   // not "plate": the plane y = 0, s = e_x (the "bp5" layout)
      cycle::InterpolateSubstrate(fes, cycle::MakePlateSubstrate(L, Vp), U_pl1);
   }
   op.SetFiniteWidth(BetaOverMu(), U_pl1, nullptr);
   const int nd = op.GetNumFaultDOFs();
   Vector slip(2 * nd);
   for (int k = 0; k < nd; k++) { slip(2 * k) = 0.0; slip(2 * k + 1) = creep ? Vp * T : 0.0; }
   ParGridFunction u(&fes);
   op.Solve(T, slip, u);
   Vector trac, sn, xyz;
   op.ComputeTraction(u, slip, trac, &sn);
   op.GetFaultDOFCoords3D(xyz);
   VectorFunctionCoefficient ex(2, [exact, T](const Vector &x, Vector &v) { exact(x, T, v); });
   r.e_u = RelL2(u, ex);
   const real_t tref = -CWa() * kMu * Vp * T / kW;   // = -mu Vp T / (2 l)
   real_t s2[3] = {0.0, 0.0, 0.0};   // tau err^2, count, max |tau|
   for (int k = 0; k < nd; k++)
   {
      if (std::abs(xyz(3 * k)) > Lx - h) { continue; }   // the faces at the walls
      const real_t want = creep ? 0.0 : tref;
      s2[0] += (trac(2 * k + 1) - want) * (trac(2 * k + 1) - want);
      s2[1] += 1.0;
      s2[2] = std::max(s2[2], std::abs(trac(2 * k + 1)));
   }
   MPI_Allreduce(MPI_IN_PLACE, s2, 2, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);
   r.e_tau = std::sqrt(s2[0] / s2[1]) / std::abs(tref);
   r.tau_max = RMax(s2[2], comm) / std::abs(tref);
   return r;
}

static void TestPlateSubstrate(MPI_Comm comm)
{
   const real_t Vp = 1.0e-9, T = 3.15576e7;   // the map-view deck's plate rate, one year
   const real_t rate = CWa() * kMu * Vp * T / kW;
   if (g_rank == 0)
   {
      std::cout << "\n=== U8 loading_plane substrate: a locked fault is loaded at c mu Vp / W = "
                << std::setprecision(6) << rate / 1e3 << " kPa/yr, creep at Vp is stress-free ===\n";
   }
   for (int p = 1; p <= 2; p++)
   {
      const int npers[3] = {2, 4, 8};
      U8Result r[3];
      for (int i = 0; i < 3; i++) { r[i] = RunPlateSubstrate(comm, p, npers[i], 0, T, Vp); }
      for (int i = 0; i < 3; i++)
      {
         if (g_rank == 0)
         {
            std::cout << std::scientific << std::setprecision(3) << "  locked, p = " << p
                      << ", h = l / " << npers[i] << ": rel L2 error of u " << r[i].e_u
                      << ", RMS |tau_s + c mu Vp T / W| / (c mu Vp T / W) " << r[i].e_tau << "\n";
         }
      }
      TEST_ASSERT(r[0].frame_ok && r[2].frame_ok, "U8 p = " << p
                  << ": one straight frame (n along y, s along x)");
      TEST_ASSERT(r[2].e_tau < r[1].e_tau && r[1].e_tau < r[0].e_tau && r[2].e_tau < 0.02,
                  "U8 p = " << p << ": a locked fault on the plane is loaded at tau_s = -c mu Vp t"
                  " / W, the error falls with h and is " << r[2].e_tau << " (< 2%) at h = l / 8");
      const U8Result c = RunPlateSubstrate(comm, p, 4, 1, T, Vp);
      if (g_rank == 0)
      {
         std::cout << std::scientific << std::setprecision(3) << "  creep at Vp, p = " << p
                   << ", h = l / 4: rel L2 error of u = u_pl " << c.e_u << ", max |tau_s| / (c mu Vp T"
                   << " / W) " << c.tau_max << "\n";
      }
      TEST_ASSERT(c.e_u < 1e-10 && c.tau_max < 1e-8, "U8 p = " << p << ": steady creep at Vp keeps"
                  " u = u_pl (" << c.e_u << ") and tau_s = 0 (" << c.tau_max << " of c mu Vp T / W)");
   }
   const U8Result k = RunPlateSubstrate(comm, 2, 4, 2, T, Vp), b = RunPlateSubstrate(comm, 2, 4, 0, T, Vp);
   if (g_rank == 0)
   {
      std::cout << std::scientific << std::setprecision(3) << "  control without the substrate"
                << " load, p = 2, h = l / 4: tau_s error " << k.e_tau << " (with it " << b.e_tau << ")\n";
   }
   TEST_ASSERT(k.e_tau > 0.25 && k.e_tau > 10.0 * b.e_tau, "U8 control: without the substrate load"
               " the walls do not load the fault (tau_s error " << k.e_tau << ", > 0.25 and > 10 x "
               << b.e_tau << ")");
}

// ====================================================================== U3
static void TestStressDrop(MPI_Comm comm)
{
   if (g_rank == 0) { std::cout << "\n=== U3 stress drop -> slip, a = 10 l ===\n"; }
   int np = 1;
   MPI_Comm_size(comm, &np);
   if (np != 1)
   {
      if (g_rank == 0) { std::cout << "  (one rank only; skipped on " << np << ")\n"; }
      return;
   }
   const real_t l = Ell(), a = 10.0 * l, dtau = 1.0e6;
   const int p = 2;
   const real_t h = l / 2.0, Lx = a + 4.0 * l, Ly = 4.0 * l;
   const int nx = (int) std::lround(2.0 * Lx / h), ny = (int) std::lround(2.0 * Ly / h);
   const real_t D_long = 2.0 * l * dtau / kMu;   // = W dtau / (c mu)
   real_t mid[2] = {0.0, 0.0};
   for (int b = 0; b < 2; b++)
   {
      std::unique_ptr<Mesh> m(CrossGrid(nx, ny, -Lx, -Ly, h, h, Walls(5),
                                        [a](real_t x, real_t y) { return std::abs(y) < 1e-6 && std::abs(x) < a; }));
      ParMesh pm(comm, *m);
      m.reset();
      BoundaryConfig bc;
      bc.fault_attr = 3;
      bc.dirichlet_attrs = {5};
      bc.default_dirichlet_func = [](const Vector &, real_t, Vector &u) { u.SetSize(2); u = 0.0; };
      LinearElastic model(kLam, kMu);
      Op op(pm, p, model, 0.0, 1e4, 1e4, bc, DGMethod::IP, SolverType::MUMPS, MapViewFrames());
      ParFiniteElementSpace &fes = op.GetFESpace();
      if (b == 0)
      {
         Vector zero(fes.GetVSize());
         zero = 0.0;
         op.SetFiniteWidth(BetaOverMu(), zero, nullptr);
      }
      const int nd = op.GetNumFaultDOFs();
      DenseMatrix G(nd, nd);
      Vector slip(2 * nd), trac, xyz;
      ParGridFunction u(&fes);
      for (int j = 0; j < nd; j++)
      {
         slip = 0.0;
         slip(2 * j + 1) = 1.0;
         op.Solve(0.0, slip, u);
         op.ComputeTraction(u, slip, trac, nullptr);
         for (int i = 0; i < nd; i++) { G(i, j) = trac(2 * i + 1); }
      }
      op.GetFaultDOFCoords3D(xyz);
      Vector rhs(nd), delta(nd);
      rhs = dtau;
      DenseMatrixInverse Ginv(G);
      Ginv.Mult(rhs, delta);
      real_t sum = 0.0, cnt = 0.0, dmin = 1e300, dmax = -1e300;
      for (int i = 0; i < nd; i++)
      {
         if (std::abs(xyz(3 * i)) < l) { sum += delta(i); cnt += 1.0; dmin = std::min(dmin, delta(i)); dmax = std::max(dmax, delta(i)); }
      }
      mid[b] = sum / cnt;
      if (g_rank == 0)
      {
         std::cout << std::scientific << std::setprecision(4) << "  " << (b == 0 ? "with the term" : "plane strain (no term)")
                   << ": " << nd << " fault DOFs, slip in |x| < l from " << dmin << " to " << dmax
                   << " m, mean " << mid[b] << " m; W dtau / (c mu) = " << D_long << " m (ratio "
                   << std::fixed << std::setprecision(4) << mid[b] / D_long << ")\n";
      }
      if (b == 0)
      {
         TEST_ASSERT(std::abs(mid[0] / D_long - 1.0) < 0.01 && std::abs(dmax / dmin - 1.0) < 0.01,
                     "U3 the slip away from the ends is W dtau / (c mu) to 1% (ratio "
                     << mid[0] / D_long << ", flat to " << dmax / dmin - 1.0 << ")");
      }
   }
   if (g_rank == 0)
   {
      std::cout << std::fixed << std::setprecision(2) << "  (the term reduces the middle slip by "
                << mid[1] / mid[0] << " x against plane strain in this clamped box)\n";
   }
}

// ====================================================================== U4
static MaterialField Material() { return MaterialField::MakeConstant(kLam, kMu, kRho); }

static std::unique_ptr<NitscheWaveOperator> FreeBoxOp(ParMesh &pm, int p)
{
   BoundaryConfig bc;
   bc.fault_attr = 3;   // absent
   bc.natural_attrs = {1};
   return std::make_unique<NitscheWaveOperator>(pm, p, Material(), bc);
}

static Mesh *Box(int n, real_t L, int attr)
{
   return CrossGrid(n, n, -L, -L, 2 * L / n, 2 * L / n, Walls(attr));
}

struct MassMatrix
{
   std::unique_ptr<ParBilinearForm> form;
   void Build(NitscheWaveOperator &op)
   {
      form = std::make_unique<ParBilinearForm>(&op.Elasticity().GetFESpace());
      form->AddDomainIntegrator(new VectorMassIntegrator(
         const_cast<NitscheMaterialSample &>(op.Sample()).Rho()));
      form->Assemble();
      form->Finalize();
   }
   void Mult(const Vector &x, Vector &y) const { y.SetSize(x.Size()); form->SpMat().Mult(x, y); }
};

static real_t FwFactor(real_t dt, real_t beta_over_rho)
{
   return 1.0 / std::sqrt(1.0 + 0.25 * beta_over_rho * dt * dt);
}

/// lambda_max of M^-1 K (K through StiffnessMult, M^-1 the operator's
/// NitscheMass) by power iteration with the Rayleigh quotient v.Kv / v.Mv.
static real_t LambdaMax(NitscheWaveOperator &op, MPI_Comm comm)
{
   MassMatrix M;
   M.Build(op);
   const int n = op.Height();
   Vector v(n), Kv(n), Mv(n), w(n);
   v.Randomize(911 + g_rank);
   v -= 0.5;
   v /= Norm(v, comm);
   real_t lam = 0.0, lam_old = 0.0;
   for (int it = 0; it < 50000; it++)
   {
      op.StiffnessMult(v, Kv);
      M.Mult(v, Mv);
      lam = Dot(v, Kv, comm) / Dot(v, Mv, comm);
      if (it > 20 && std::abs(lam - lam_old) <= 1e-14 * std::abs(lam)) { break; }
      lam_old = lam;
      op.Mass().Apply(Kv, w);   // M^-1 K v
      v = w;
      v /= Norm(w, comm);
   }
   return lam;
}

static void TestDynamic(MPI_Comm comm)
{
   if (g_rank == 0) { std::cout << "\n=== U4 the dynamic operator ===\n"; }
   const real_t bor = BetaOverMu() * kMu / kRho;   // beta / rho
   const int p = 2;
   // (a) the matrix
   {
      std::unique_ptr<Mesh> m(Box(8, 30.0e3, 1));
      ParMesh pm(comm, *m);
      m.reset();
      auto op0 = FreeBoxOp(pm, p), op1 = FreeBoxOp(pm, p);
      Vector zero(op1->Height());
      zero = 0.0;
      op1->SetFiniteWidth(BetaOverMu(), zero, nullptr);
      MassMatrix M;
      M.Build(*op0);
      Vector x(op0->Height()), K0x, K1x, Mx;
      x.Randomize(17 + g_rank);
      op0->StiffnessMult(x, K0x);
      op1->StiffnessMult(x, K1x);
      M.Mult(x, Mx);
      Vector d(K1x);
      d -= K0x;
      Mx *= bor;
      const real_t e = MaxDiff(d, Mx, comm) / MaxAbs(Mx, comm);
      TEST_ASSERT(e <= 1e-12, "U4a (K + M_beta) x - K x = (beta / rho) M x (" << std::scientific
                  << std::setprecision(2) << e << " of max)");
   }
   // (b) the cutoff: a uniform displacement oscillates at sqrt(beta / rho)
   {
      std::unique_ptr<Mesh> m(Box(8, 30.0e3, 1));
      ParMesh pm(comm, *m);
      m.reset();
      auto op = FreeBoxOp(pm, p);
      Vector zero(op->Height());
      zero = 0.0;
      op->SetFiniteWidth(BetaOverMu(), zero, nullptr);
      real_t dt = CUse(p) * ComputeMaxDt(op->Sample(), pm, p, 1.0);
      dt *= FwFactor(dt, bor);
      const real_t th = std::acos(1.0 - 0.5 * bor * dt * dt);
      const int n = (int) std::ceil(2.0 * M_PI / th);
      ParGridFunction g(&op->Elasticity().GetFESpace());
      Vector c(2);
      c(0) = 1.0;
      c(1) = 0.5;
      VectorConstantCoefficient cc(c);
      g.ProjectCoefficient(cc);
      Vector U(g), V(U.Size()), U0(g);
      V = 0.0;
      op->SetTime(0.0);
      op->SetTimeStep(dt);
      CentralDifferenceSolver ode;
      ode.Init(*op);
      real_t t = 0.0, e = 0.0;
      for (int s = 1; s <= n; s++)
      {
         real_t dts = dt;
         ode.Step(U, V, t, dts);
         Vector ex(U0);
         ex *= std::cos(s * th);
         e = std::max(e, MaxDiff(U, ex, comm));
      }
      if (g_rank == 0)
      {
         std::cout << std::scientific << std::setprecision(4) << "  cutoff sqrt(beta / rho) = "
                   << std::sqrt(bor) << " rad/s (f_c " << std::sqrt(bor) / (2 * M_PI) << " Hz), dt "
                   << dt << " s, " << n << " steps\n";
      }
      TEST_ASSERT(e <= 1e-9 * 1.0, "U4b a uniform displacement oscillates as cos(n theta), "
                  "cos theta = 1 - beta dt^2 / (2 rho), over one period (" << std::scientific
                  << std::setprecision(2) << e << " m of 1 m)");
   }
   // (c) dispersion above the cutoff, evanescence below it
   {
      real_t nrm[2] = {1.0, 0.7};
      const real_t nn = std::hypot(nrm[0], nrm[1]);
      nrm[0] /= nn;
      nrm[1] /= nn;
      const real_t kmag = 2.0 * M_PI / 120.0e3, bom_k = kmag * kmag;   // l_k = 1 / k
      const real_t bor_k = bom_k * kMu / kRho;
      std::unique_ptr<ParMesh> pm[2];
      for (int r = 0; r < 2; r++)
      {
         std::unique_ptr<Mesh> sm(Box(8 << r, 30.0e3, 5));
         pm[r] = std::make_unique<ParMesh>(comm, *sm);
      }
      for (int kind = 0; kind < 4; kind++)
      {
         const bool pw = (kind % 2 == 0), evan = (kind >= 2);
         const real_t c = pw ? Cp() : Cs();
         real_t d[2];
         if (pw) { d[0] = nrm[0]; d[1] = nrm[1]; }
         else { d[0] = -nrm[1]; d[1] = nrm[0]; }
         // propagating: u = d sin(k s - omega t), omega^2 = c^2 k^2 + beta / rho;
         // evanescent: u = d e^{-kappa s} cos(omega t), omega = sqrt(beta / rho) / 2,
         // kappa^2 = (beta / rho - omega^2) / c^2 (s = n . x)
         const real_t om = evan ? 0.5 * std::sqrt(bor_k) : std::sqrt(c * c * kmag * kmag + bor_k);
         const real_t kap = evan ? std::sqrt((bor_k - om * om) / (c * c)) : 0.0;
         ExactFn ue = [=](const Vector &x, real_t t, Vector &v)
         {
            const real_t s = nrm[0] * x(0) + nrm[1] * x(1);
            const real_t a = evan ? std::exp(-kap * s) * std::cos(om * t) : std::sin(kmag * s - om * t);
            v.SetSize(2);
            v(0) = d[0] * a;
            v(1) = d[1] * a;
         };
         ExactFn ve = [=](const Vector &x, real_t t, Vector &v)
         {
            const real_t s = nrm[0] * x(0) + nrm[1] * x(1);
            const real_t a = evan ? -om * std::exp(-kap * s) * std::sin(om * t)
                                  : -om * std::cos(kmag * s - om * t);
            v.SetSize(2);
            v(0) = d[0] * a;
            v(1) = d[1] * a;
         };
         // the exact data, stepped with (with_term) or without the term
         auto run = [&](ParMesh &m, real_t dt, int nsteps, bool with_term)
         {
            BoundaryConfig bc;
            bc.fault_attr = 3;
            bc.dirichlet_attrs = {5};
            bc.default_dirichlet_func = ue;
            NitscheWaveOperator op(m, p, Material(), bc);
            Vector zero(op.Height());
            zero = 0.0;
            if (with_term) { op.SetFiniteWidth(bom_k, zero, nullptr); }
            auto &fes = op.Elasticity().GetFESpace();
            ParGridFunction gu(&fes), gv(&fes);
            VectorFunctionCoefficient cu(2, ue), cv(2, ve);
            cu.SetTime(0.0);
            cv.SetTime(0.0);
            gu.ProjectCoefficient(cu);
            gv.ProjectCoefficient(cv);
            Vector U(gu), V(gv);
            op.SetTime(0.0);
            op.SetTimeStep(dt);
            real_t t = 0.0;
            CentralDifferenceSolver ode;
            ode.Init(op);
            for (int s = 0; s < nsteps; ++s)
            {
               real_t dts = dt;
               ode.Step(U, V, t, dts);
            }
            gu = U;
            cu.SetTime(t);
            return RelL2(gu, cu);
         };
         const real_t T = evan ? M_PI / om : 2.0;   // evanescent: half a period
         real_t dt;
         {
            BoundaryConfig bc;
            bc.fault_attr = 3;
            bc.dirichlet_attrs = {5};
            bc.default_dirichlet_func = ue;
            NitscheWaveOperator op(*pm[1], p, Material(), bc);
            dt = CUse(p) * ComputeMaxDt(op.Sample(), *pm[1], p, 1.0);
         }
         dt *= FwFactor(dt, bor_k);
         int nsteps = (int) std::ceil(T / dt);
         dt = T / nsteps;
         real_t e_dt = run(*pm[1], dt, nsteps, true), change = 1.0;
         for (int hh = 0; hh < 4; ++hh)
         {
            const real_t e_half = run(*pm[1], 0.5 * dt, 2 * nsteps, true);
            change = std::abs(e_dt - e_half) / e_half;
            if (change < 0.01) { break; }
            dt *= 0.5;
            nsteps *= 2;
            e_dt = e_half;
         }
         const real_t e8 = run(*pm[0], dt, nsteps, true), e16 = e_dt;
         const real_t e_nt = run(*pm[1], dt, nsteps, false);
         const real_t o = std::log2(e8 / e16);
         const std::string name = std::string(evan ? "evanescent " : "plane ") + (pw ? "P" : "S");
         if (g_rank == 0)
         {
            std::cout << std::scientific << std::setprecision(3) << "  " << name << ": omega = "
                      << om << " rad/s (cutoff " << std::sqrt(bor_k) << ", c k = " << c * kmag
                      << ")" << (evan ? ", decay 1 / kappa = " : "");
            if (evan) { std::cout << 1.0 / kap << " m"; }
            std::cout << "; rel L2 error at T = " << T << " s on 8 x 8 / 16 x 16: " << e8 << ", "
                      << e16 << " (order " << std::fixed << std::setprecision(2) << o
                      << "); stepped without the term: " << std::scientific << e_nt << "\n";
         }
         TEST_ASSERT(change < 0.01 && o > p + 0.5, "U4c " << name << (evan ? " (omega = f_c / 2)"
                     : " omega^2 = c^2 k^2 + beta / rho") << ": L2 order " << o << " > p + 0.5");
         TEST_ASSERT(e_nt > 10.0 * e16, "U4c " << name << " control: stepped without the term it "
                     "misses by " << e_nt / e16 << " x (> 10)");
      }
   }
   // (d) rule C3: Kelvin-Voigt damps K alone
   {
      std::unique_ptr<Mesh> m(Box(8, 30.0e3, 1));
      ParMesh pm(comm, *m);
      m.reset();
      auto op = FreeBoxOp(pm, p);
      const real_t eta = 0.1;
      op->SetKelvinVoigt(eta, NitscheWaveOperator::KVScale::Node, -1.0);
      Vector zero(op->Height());
      zero = 0.0;
      op->SetFiniteWidth(BetaOverMu(), zero, nullptr);
      ParGridFunction g(&op->Elasticity().GetFESpace());
      Vector c(2);
      c(0) = 1.0;
      c(1) = -0.3;
      VectorConstantCoefficient cc(c);
      g.ProjectCoefficient(cc);
      Vector V(g), A;
      op->SetTime(0.0);
      op->Mult(zero, V, A);
      const real_t h = 2 * 30.0e3 / 8 / 2.0;   // the smallest triangle's legs ~ h / 2
      const real_t tau_kv = eta * NitscheWaveOperator::kKVGamma * (h / p) / Cs();
      const real_t scale = bor * tau_kv * 1.0;
      const real_t e = MaxAbs(A, comm) / scale;
      TEST_ASSERT(e <= 1e-10, "U4d rule C3: Kelvin-Voigt on K alone, M^-1 C V = 0 for a uniform V ("
                  << std::scientific << std::setprecision(2) << e
                  << " of (beta / rho) tau_kv |V|)");
   }
   // (e) rule C6
   {
      std::unique_ptr<Mesh> m(Box(4, 30.0e3, 1));
      ParMesh pm(comm, *m);
      m.reset();
      auto op0 = FreeBoxOp(pm, p);
      const real_t lam0 = LambdaMax(*op0, comm);
      const real_t bom_big = 2.0 * lam0 * kRho / kMu, bor_big = 2.0 * lam0;
      auto op1 = FreeBoxOp(pm, p);
      Vector zero(op1->Height());
      zero = 0.0;
      op1->SetFiniteWidth(bom_big, zero, nullptr);
      const real_t lam1 = LambdaMax(*op1, comm);
      const real_t e = std::abs((lam1 - lam0) / bor_big - 1.0);
      if (g_rank == 0)
      {
         std::cout << std::scientific << std::setprecision(6) << "  omega_max^2 = " << lam0
                   << " (no term), " << lam1 << " with beta / rho = " << bor_big << "\n";
      }
      TEST_ASSERT(e <= 1e-6, "U4e rule C6: lambda_max(M^-1 (K + M_beta)) = lambda_max(M^-1 K) + "
                  "beta / rho (" << std::scientific << std::setprecision(2) << e << ")");
      const real_t dt0 = 0.95 * 2.0 / std::sqrt(lam0), dtc = dt0 * FwFactor(dt0, bor_big);
      auto grows = [&](real_t dt)
      {
         auto op = FreeBoxOp(pm, p);
         op->SetFiniteWidth(bom_big, zero, nullptr);
         Vector U(op->Height()), V(op->Height());
         U.Randomize(5 + g_rank);
         U -= 0.5;
         V = 0.0;
         const real_t n0 = Norm(U, comm);
         op->SetTime(0.0);
         op->SetTimeStep(dt);
         CentralDifferenceSolver ode;
         ode.Init(*op);
         real_t t = 0.0, nmax = n0;
         for (int s = 0; s < 2000; s++)
         {
            real_t dts = dt;
            ode.Step(U, V, t, dts);
            nmax = std::max(nmax, Norm(U, comm));
            if (!std::isfinite(nmax) || nmax > 1e12 * n0) { break; }
         }
         return nmax / n0;
      };
      const real_t g0 = grows(dt0), gc = grows(dtc);
      if (g_rank == 0)
      {
         std::cout << std::scientific << std::setprecision(3) << "  2000 steps at 0.95 x 2 / omega_max = "
                   << dt0 << " s: max |U| / |U0| = " << g0 << "; at the cut step " << dtc
                   << " s: " << gc << "\n";
      }
      TEST_ASSERT(!(g0 < 1e6) && gc < 10.0, "U4e rule C6: the uncut step blows up (" << g0
                  << "), the cut step stays bounded (" << gc << ")");
   }
}

// ====================================================================== U7
// The normal response of the isotropic term (the note's Step 7b, eqs.
// 7.5-7.7): a uniform normal jump D_n on a long straight fault.  The fault
// slip carries no opening, so the plus half is solved on its own: the strip
// [-2 l_n, 2 l_n] x [0, 6 l_n], the face y = 0 displaced by D_n / 2 normal
// to itself, every wall at the exact u_y = (D_n / 2) e^{-y / l_n}, l_n =
// sqrt(M / beta), M = lambda + 2 mu.  The profile (L2 order above p + 0.5)
// and the face stress sigma_yy(0+) = -(D_n / 2) sqrt(M beta) from the
// solution's gradient on the faces |x| < l_n (relative RMS error falling
// with h, below 1% at h = l_n / 8 for p = 2); printed against the 3D mode I
// value 4 mu / (pi (1 - nu) W) per unit jump (eq. 7.9).
static void TestNormalStiffness(MPI_Comm comm)
{
   const real_t M = kLam + 2.0 * kMu, beta = BetaOverMu() * kMu, Dn = 1.0;
   const real_t ln = std::sqrt(M / beta);
   const real_t sref = -0.5 * Dn * std::sqrt(M * beta);
   if (g_rank == 0)
   {
      std::cout << "\n=== U7 normal stiffness: l_n = sqrt(M / beta) = " << ln << " m, sigma_yy(0+) = "
                << "-(D_n / 2) sqrt(M beta) = " << sref << " Pa for D_n = 1 m ===\n";
   }
   for (int p = 1; p <= 2; p++)
   {
      const int npers[3] = {2, 4, 8};
      real_t eu[3], es[3];
      for (int r = 0; r < 3; r++)
      {
         const real_t h = ln / npers[r], Lx = 2.0 * ln, Ly = 6.0 * ln;
         const int nx = (int) std::lround(2.0 * Lx / h), ny = (int) std::lround(Ly / h);
         std::unique_ptr<Mesh> m(CrossGrid(nx, ny, -Lx, 0.0, h, h, Walls(5)));
         ParMesh pm(comm, *m);
         m.reset();
         auto exact = [ln, Dn](const Vector &x, Vector &u)
         {
            u.SetSize(2);
            u(0) = 0.0;
            u(1) = 0.5 * Dn * std::exp(-x(1) / ln);
         };
         BoundaryConfig bc;
         bc.fault_attr = 3;
         bc.dirichlet_attrs = {5};
         bc.default_dirichlet_func = [exact](const Vector &x, real_t, Vector &u) { exact(x, u); };
         LinearElastic model(kLam, kMu);
         Op op(pm, p, model, 0.0, 1e4, 1e4, bc, DGMethod::IP, SolverType::MUMPS);
         ParFiniteElementSpace &fes = op.GetFESpace();
         Vector zero(fes.GetVSize());
         zero = 0.0;
         op.SetFiniteWidth(BetaOverMu(), zero, nullptr);
         Vector slip(2 * op.GetNumFaultDOFs());
         slip = 0.0;
         ParGridFunction u(&fes);
         op.Solve(0.0, slip, u);
         VectorFunctionCoefficient ex(2, exact);
         eu[r] = RelL2(u, ex);
         // sigma_yy on the faces y = 0, |x| < l_n, at their quadrature points
         real_t acc[2] = {0.0, 0.0};
         DenseMatrix grad;
         for (int be = 0; be < pm.GetNBE(); be++)
         {
            FaceElementTransformations *Tr = pm.GetBdrFaceTransformations(be);
            if (!Tr) { continue; }
            Tr->SetAllIntPoints(&Geometries.GetCenter(Tr->GetGeometryType()));
            Vector xc;
            Tr->Transform(Geometries.GetCenter(Tr->GetGeometryType()), xc);
            if (std::abs(xc(1)) > 1e-6 * ln || std::abs(xc(0)) > ln) { continue; }
            const IntegrationRule &ir = IntRules.Get(Tr->GetGeometryType(), 2 * p + 2);
            for (int q = 0; q < ir.GetNPoints(); q++)
            {
               Tr->SetAllIntPoints(&ir.IntPoint(q));
               u.GetVectorGradient(*Tr->Elem1, grad);
               const real_t syy = kLam * (grad(0, 0) + grad(1, 1)) + 2.0 * kMu * grad(1, 1);
               acc[0] += (syy - sref) * (syy - sref);
               acc[1] += 1.0;
            }
         }
         MPI_Allreduce(MPI_IN_PLACE, acc, 2, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm);
         es[r] = std::sqrt(acc[0] / acc[1]) / std::abs(sref);
         if (g_rank == 0)
         {
            std::cout << std::scientific << std::setprecision(3) << "  p = " << p << ", h = l_n / "
                      << npers[r] << ": rel L2 error of u " << eu[r] << ", RMS |sigma_yy(0+) - ref| / |ref| "
                      << es[r] << "\n";
         }
      }
      const real_t o = std::log2(eu[1] / eu[2]);
      TEST_ASSERT(o > p + 0.5, "U7 p = " << p << ": normal-jump profile e^{-y / l_n}, L2 order " << o
                  << " > p + 0.5");
      TEST_ASSERT(es[2] < es[1] && es[1] < es[0] && (p == 1 || es[2] < 0.01),
                  "U7 p = " << p << ": face stress (D_n / 2) sqrt(M beta), error falls with h ("
                  << es[2] << " at h = l_n / 8" << (p == 2 ? ", < 1%" : "") << ")");
   }
   const real_t k25 = 0.5 * std::sqrt(M * beta), k3d = 4.0 * kMu / (M_PI * (1.0 - Nu()) * kW);
   if (g_rank == 0)
   {
      std::cout << std::fixed << std::setprecision(3) << "  normal stiffness per unit jump: 2.5D "
                << k25 * kW / kMu << " mu / W, 3D mode I " << k3d * kW / kMu << " mu / W: ratio "
                << k25 / k3d << " (the isotropic beta's known model error, Step 7e)\n";
   }
}

// ====================================================================== U5
namespace u5
{

/// A smooth tapered strike slip on the deck's straight fault, psi0.
void FaultState(const cycle::QDPhase &qd, const cycle::Prestress &pre, real_t S0, real_t Lf,
                Vector &s, Vector &psi)
{
   const Vector &X = qd.Table().coords;
   const int N = qd.Table().Size();
   s.SetSize(2 * N);
   psi.SetSize(N);
   for (int i = 0; i < N; i++)
   {
      const real_t x = X(3 * i);
      const real_t c = (std::abs(x) < Lf) ? std::cos(0.5 * M_PI * x / Lf) : 0.0;
      s(2 * i) = 0.0;                     // kernel pairs: (dip, strike)
      s(2 * i + 1) = S0 * c * c;
      psi(i) = pre.psi0(i);
   }
}

struct Balance
{
   Vector u, sk, Vqd, psi, sn;
   real_t t = 0.0, Vmax = 0.0;
};

/// 1 on the local DOFs of the elements with a vertex on the fault line
/// (y = 0, x_min <= x <= x_max of the fault DOFs), 0 elsewhere.
std::vector<char> FaultElementMask(NitscheWaveOperator &W, real_t xmin, real_t xmax)
{
   ParMesh &m = W.GetMesh();
   const ParFiniteElementSpace &fes = W.Elasticity().GetFESpace();
   std::vector<char> on(fes.GetVSize(), 0);
   Array<int> v, vd;
   const real_t tol = 1e-6;
   for (int e = 0; e < m.GetNE(); e++)
   {
      m.GetElementVertices(e, v);
      bool f = false;
      for (int k = 0; k < v.Size(); k++)
      {
         const real_t *X = m.GetVertex(v[k]);
         if (std::abs(X[1]) <= tol && X[0] >= xmin - tol && X[0] <= xmax + tol) { f = true; }
      }
      if (!f) { continue; }
      fes.GetElementVDofs(e, vd);
      for (int i = 0; i < vd.Size(); i++) { on[vd[i] >= 0 ? vd[i] : -1 - vd[i]] = 1; }
   }
   return on;
}

real_t MaxAbsOff(const Vector &x, const std::vector<char> &mask, MPI_Comm c)
{
   real_t m = 0.0;
   for (int i = 0; i < x.Size(); i++) { if (!mask[i]) { m = std::max(m, std::abs(x(i))); } }
   return RMax(m, c);
}

real_t MaxDiffOff(const Vector &x, const Vector &y, const std::vector<char> &mask, MPI_Comm c)
{
   real_t m = 0.0;
   for (int i = 0; i < x.Size(); i++) { if (!mask[i]) { m = std::max(m, std::abs(x(i) - y(i))); } }
   return RMax(m, c);
}

/// The hand-back table's largest |tau_qd - tau_dyn| and |sigma_qd - sigma_dyn|
/// (the driver's WriteHandbackTable) and max |tau_qd|.
void HandbackJump(const cycle::StationFields &fd, const cycle::StationFields &fq, MPI_Comm c,
                  real_t &dtau, real_t &dsig, real_t &tmax)
{
   auto tau = [](const cycle::StationFields &f, int j, bool strike)
   {
      return strike ? -(f.tau_a_strike(j) + f.tau_b_strike(j) + f.eta * f.V_strike(j))
                    : -(f.tau_a_dip(j) + f.tau_b_dip(j) + f.eta * f.V_dip(j));
   };
   dtau = dsig = tmax = 0.0;
   for (int j = 0; j < fd.V_strike.Size(); j++)
   {
      dtau = std::max(dtau, std::hypot(tau(fq, j, true) - tau(fd, j, true),
                                       tau(fq, j, false) - tau(fd, j, false)));
      dsig = std::max(dsig, std::abs(fq.sigma_n(j) - fd.sigma_n(j)));
      tmax = std::max(tmax, std::hypot(tau(fq, j, true), tau(fq, j, false)));
   }
   dtau = RMax(dtau, c);
   dsig = RMax(dsig, c);
   tmax = RMax(tmax, c);
}

struct Setup
{
   std::unique_ptr<cycle::QDPhase> qd;
   std::unique_ptr<cycle::FaultParams> par;
   std::unique_ptr<cycle::Prestress> pre;
};

void Build(ParMesh &pmesh, MPIContext &mpi, const cycle::CycleConfig &cfg, Setup &S)
{
   S.qd = std::make_unique<cycle::QDPhase>(pmesh, mpi, cfg);
   const spatial::SpatialFrictionConfig &s = cfg.shared;
   const MaterialField material = MaterialField::MakeConstant(
      s.material_fallback.lambda, s.material_fallback.mu, s.material_fallback.rho);
   const spatial::PorePressureSpec no_pp;
   S.par = std::make_unique<cycle::FaultParams>(cycle::BuildFaultParams(
      S.qd->Table(), *s.rate_state, cfg.fields, material, pmesh, no_pp));
   cycle::PrestressConfig pc = cfg.prestress;
   pc.a_max = cycle::GlobalMaxA(*S.par, mpi.GetComm());
   S.pre = std::make_unique<cycle::Prestress>(cycle::BuildPrestress(*S.par, pc));
   S.qd->Setup(*S.par, *S.pre);
}

Balance Solve(cycle::QDPhase &qd, const cycle::Prestress &pre, real_t ts, real_t Lf)
{
   Balance b;
   Vector s0, psi0;
   FaultState(qd, pre, 0.02, Lf, s0, psi0);
   qd.SetFaultState(s0, psi0);
   b.Vmax = qd.Solve(ts);
   qd.Slip(b.sk);
   qd.Psi(b.psi);
   qd.SigmaN(b.sn);
   b.Vqd = qd.SlipRate();
   b.u = qd.Displacement();
   b.t = ts;
   return b;
}

void Enter(cycle::DynPhase &dyn, const Balance &b, Vector *d_out = nullptr)
{
   Vector d, V, psi;
   cycle::bridge::ToNitsche(dyn.Frames(), b.sk, b.Vqd, b.psi, d, V, psi);
   dyn.Enter(b.u, d, V, psi, b.sn, b.t);
   if (d_out) { *d_out = d; }
}

/// max |M^-1 (F(t_s) - K u_s + G d_s)| over every DOF, fault elements
/// included: the dynamic operator's residual at the QD balance with the
/// slip held (no friction traction), relative to a_K = max |M^-1 K u_s|.
real_t SlipHeldResidual(cycle::DynPhase &dyn, const Balance &b, MPI_Comm c, real_t &aK)
{
   Vector d, k, zero(b.u.Size()), Ku, MKu;
   zero = 0.0;
   Enter(dyn, b, &d);
   dyn.Op().AccelerationWithSlip(b.u, zero, d, k);
   dyn.Op().StiffnessMult(b.u, Ku);
   dyn.Op().Mass().Apply(Ku, MKu);
   aK = MaxAbs(MKu, c);
   return MaxAbs(k, c) / aK;
}

/// n steps, then the hand-back to `qd` at t_s + tau (ToQD, Reenter): the
/// jump of the hand-back table.
void Handback(cycle::DynPhase &dyn, cycle::QDPhase &qd, const cycle::CycleConfig &cfg, int n,
              MPI_Comm c, real_t &dtau, real_t &dsig, real_t &tmax)
{
   for (int k = 0; k < n; k++) { dyn.Step(); }
   cycle::StationFields fd, fq;
   dyn.FillStationFields(fd);
   Vector sk, psi;
   cycle::bridge::ToQD(dyn.Frames(), dyn.Slip(), dyn.Psi(), sk, psi);
   qd.Reenter(sk, psi, dyn.Time(), 0, cfg.dt_init);
   qd.FillStationFields(fq);
   HandbackJump(fd, fq, c, dtau, dsig, tmax);
}

} // namespace u5

static void TestRegimeSwitch(MPIContext &mpi, const char *deck, int n_dyn)
{
#ifndef MFEM_USE_PETSC
   (void) mpi; (void) deck; (void) n_dyn;
   TEST_ASSERT(false, "U5 needs PETSc");
#else
   MPI_Comm comm = mpi.GetComm();
   cycle::CycleConfig base = cycle::CycleConfig::Load(deck);
   MFEM_VERIFY(base.switching && !base.lts.On() && !base.finite_width.On(),
               deck << ": U5 needs a switching, global-stepping deck without [cycle.finite_width]");
   std::unique_ptr<Mesh> smesh = cycle::LoadScaledMesh(base.shared.mesh.path, base.mesh_scale);
   base.CheckDimension(smesh->Dimension());
   MFEM_VERIFY(smesh->Dimension() == 2, base.shared.mesh.path << " is not a 2-D mesh");
   ParMesh pmesh(comm, *smesh);
   smesh.reset();
   if (g_rank == 0)
   {
      std::cout << "\n=== U5 regime switch: " << deck << ", " << base.shared.mesh.path << ", p = "
                << base.shared.mesh.order << ", ranks " << mpi.Size() << " ===\n";
   }
   const real_t ts = 3.0e7;
   const auto &mf = base.shared.material_fallback;
   const real_t nu = mf.lambda / (2.0 * (mf.lambda + mf.mu));

   // the reference without the term (the same state and steps)
   cycle::CycleConfig off = base;
   off.absorbing.faces = "off";
   u5::Setup S0;
   u5::Build(pmesh, mpi, off, S0);
   real_t Lf = 0.0, xmin = 1e300, xmax = -1e300;
   {
      const Vector &X = S0.qd->Table().coords;
      for (int i = 0; i < S0.qd->Table().Size(); i++)
      {
         xmin = std::min(xmin, X(3 * i));
         xmax = std::max(xmax, X(3 * i));
      }
      xmin = mpi.GlobalMin(xmin);
      xmax = mpi.GlobalMax(xmax);
      Lf = std::max(std::abs(xmin), std::abs(xmax));
   }
   const u5::Balance b0 = u5::Solve(*S0.qd, *S0.pre, ts, Lf);
   real_t dtau0 = 0.0, dsig0 = 0.0, tmax0 = 0.0, dt_cfl0 = 0.0, res0 = 0.0;
   {
      cycle::DynPhase dyn(pmesh, off, S0.qd->Table(), *S0.par, *S0.pre, 1.0, nullptr,
                          S0.qd->Domain().FrameField());
      dt_cfl0 = dyn.DtCFL();
      real_t aK0 = 0.0;
      res0 = u5::SlipHeldResidual(dyn, b0, comm, aK0);
      u5::Enter(dyn, b0);
      u5::Handback(dyn, *S0.qd, off, n_dyn, comm, dtau0, dsig0, tmax0);
   }
   if (g_rank == 0)
   {
      std::cout << std::scientific << std::setprecision(3) << "  without the term: QD balance at t_s = "
                << ts << " s, max V " << b0.Vmax << " m/s; hand-back after " << n_dyn
                << " steps: max |tau_qd - tau_dyn| " << dtau0 << " Pa, |sigma_qd - sigma_dyn| "
                << dsig0 << " Pa (max |tau| " << tmax0 << " Pa); slip-held residual "
                << res0 << " a_K\n";
   }

   for (const std::string sub : {"loading_plane", "locked_reference"})
   {
      cycle::CycleConfig cfg = off;
      cfg.finite_width.W = 15.0e3;
      cfg.finite_width.calibration = "wa2020";
      cfg.finite_width.c = (1.0 - nu) * M_PI / 2.0;
      cfg.finite_width.substrate = sub;
      const real_t bor = cfg.finite_width.BetaOverMu() * mf.mu / mf.rho;
      if (g_rank == 0)
      {
         std::cout << "\n  substrate \"" << sub << "\": " << cycle::FiniteWidthReport(cfg) << "\n";
      }
      u5::Setup S;
      u5::Build(pmesh, mpi, cfg, S);
      const Vector *Upl = S.qd->SubstrateShape();
      TEST_ASSERT(Upl != nullptr && Upl->Size() == S.qd->Domain().GetFESpace().GetVSize(),
                  "U5 " << sub << ": QDPhase::SubstrateShape on the displacement space");
      const u5::Balance b = u5::Solve(*S.qd, *S.pre, ts, Lf);
      cycle::DynPhase dyn(pmesh, cfg, S.qd->Table(), *S.par, *S.pre, 1.0, nullptr,
                          S.qd->Domain().FrameField(), Upl);
      cycle::DynPhase dyn0(pmesh, off, S.qd->Table(), *S.par, *S.pre, 1.0, nullptr,
                           S.qd->Domain().FrameField());
      NitscheWaveOperator &W = dyn.Op(), &W0 = dyn0.Op();
      const std::vector<char> fm = u5::FaultElementMask(W, xmin, xmax);

      // C1: one M_beta
      Vector KU, K0U, KQ, Ku, KQu;
      W.StiffnessMult(*Upl, KU);
      W0.StiffnessMult(*Upl, K0U);
      KU -= K0U;
      const Vector &fhat = W.Elasticity().SubstrateLoad();
      const real_t c1a = MaxDiff(KU, fhat, comm) / MaxAbs(fhat, comm);
      const real_t c1q = MaxDiff(S.qd->Domain().SubstrateLoad(), fhat, comm) / MaxAbs(fhat, comm);
      TEST_ASSERT(c1a <= 1e-12 && c1q <= 1e-12, "U5 " << sub << " C1: (K_W - K_W0) U_pl1 = "
                  "M_beta U_pl1 = the QD regime's f_hat (" << std::scientific << std::setprecision(2)
                  << c1a << ", " << c1q << ")");
      W.StiffnessMult(b.u, Ku);
      S.qd->Domain().ApplyStiffness(b.u, KQu);
      const real_t c1k = u5::MaxDiffOff(Ku, KQu, fm, comm) / MaxAbs(KQu, comm);
      TEST_ASSERT(c1k <= 1e-12, "U5 " << sub << " C1: K_W u_s = K_QD u_s off the fault elements ("
                  << std::scientific << std::setprecision(2) << c1k << ")");

      // C1 with the slip held: the QD balance is the dynamic operator's static
      // state for that slip, fault elements included (the damped dynamic run
      // with the slip held comes to rest there)
      {
         real_t aK1 = 0.0, aK0 = 0.0;
         const real_t r1 = u5::SlipHeldResidual(dyn, b, comm, aK1);
         const real_t r0 = u5::SlipHeldResidual(dyn0, b, comm, aK0);
         if (g_rank == 0)
         {
            std::cout << std::scientific << std::setprecision(3) << "  slip held at the QD balance:"
                      << " max |M^-1 (F - K u_s + G d_s)| = " << r1 << " a_K (W), " << r0
                      << " a_K (W0, no term)\n";
         }
         TEST_ASSERT(r1 <= std::max(1e-8, 10.0 * res0), "U5 " << sub << " C1: the QD balance is "
                     "the dynamic operator's static state for its slip (residual " << r1
                     << " a_K, every DOF)");
         TEST_ASSERT(r0 > 1e3 * r1, "U5 " << sub << " control: W0 (no term) leaves " << r0
                     << " a_K (> 1e3 x)");
      }

      // C6
      const real_t want = dt_cfl0 * FwFactor(dt_cfl0, bor);
      const real_t c6 = std::abs(dyn.DtCFL() / want - 1.0);
      TEST_ASSERT(c6 <= 1e-14 && std::abs(dyn.FiniteWidthDtFactor() - FwFactor(dt_cfl0, bor)) <= 1e-15,
                  "U5 " << sub << " C6: dt_cfl = " << dyn.DtCFL() << " s = dt_cfl(W0) / sqrt(1 + "
                  "beta dt^2 / (4 rho)) (1 - factor = " << std::scientific << std::setprecision(3)
                  << 1.0 - dyn.FiniteWidthDtFactor() << ")");

      // A0 and C2 at tau = 0
      u5::Enter(dyn, b);
      u5::Enter(dyn0, b);
      Vector zero(b.u.Size()), A0, A00, MKu, Fw, Fq;
      zero = 0.0;
      W.Mult(b.u, zero, A0);
      W0.Mult(b.u, zero, A00);
      W.Mass().Apply(Ku, MKu);
      const real_t aK = MaxAbs(MKu, comm);
      Vector dU(b.u);
      dU.Add(-ts, *Upl);
      dU *= bor;   // (beta / rho)(u_s - u_pl(t_s))
      const real_t aB = u5::MaxAbsOff(dU, fm, comm);
      const real_t a0 = u5::MaxAbsOff(A0, fm, comm), a00 = u5::MaxAbsOff(A00, fm, comm);
      if (g_rank == 0)
      {
         std::cout << std::scientific << std::setprecision(3) << "  QD balance max V " << b.Vmax
                   << " m/s; a_K = " << aK << " m/s^2, max off-fault (beta / rho)|u_s - u_pl| = "
                   << aB << " m/s^2; A0 off the fault elements: " << a0 << " (W), " << a00
                   << " (W0, no term)\n";
      }
      TEST_ASSERT(a0 <= 1e-8 * aK, "U5 " << sub << " A0 off the fault elements " << a0 / aK
                  << " a_K (<= 1e-8)");
      TEST_ASSERT(aB > 1e-6 * aK && a00 >= 0.5 * aB, "U5 " << sub << " control: W0 (no term) is "
                  "off by " << a00 / aB << " x (beta / rho)|u_s - u_pl| there (>= 0.5)");
      W.Elasticity().AssembleDirichletOnlyRHS(Fw, W.GetTime());
      S.qd->Domain().AssembleDirichletOnlyRHS(Fq, ts);
      const real_t c2a = MaxDiff(Fw, Fq, comm) / MaxAbs(Fq, comm);

      // hand-back (and C2 at the end of the steps)
      real_t dtau = 0.0, dsig = 0.0, tmax = 0.0;
      u5::Handback(dyn, *S.qd, cfg, n_dyn, comm, dtau, dsig, tmax);
      W.Elasticity().AssembleDirichletOnlyRHS(Fw, W.GetTime());
      S.qd->Domain().AssembleDirichletOnlyRHS(Fq, dyn.Time());
      const real_t c2b = MaxDiff(Fw, Fq, comm) / MaxAbs(Fq, comm);
      TEST_ASSERT(c2a <= 1e-13 && c2b <= 1e-13, "U5 " << sub << " C2: F_W(tau) = F_QD(t_s + tau) "
                  "at tau = 0 and " << dyn.Tau() << " s (" << std::scientific << std::setprecision(2)
                  << c2a << ", " << c2b << ")");
      if (g_rank == 0)
      {
         std::cout << std::scientific << std::setprecision(3) << "  hand-back after " << n_dyn
                   << " steps (tau = " << dyn.Tau() << " s): max |tau_qd - tau_dyn| " << dtau
                   << " Pa, |sigma_qd - sigma_dyn| " << dsig << " Pa (without the term " << dtau0
                   << ", " << dsig0 << ")\n";
      }
      TEST_ASSERT(dtau <= 2.0 * dtau0 + 1e-9 * tmax && dsig <= 2.0 * dsig0 + 1e-9 * tmax,
                  "U5 " << sub << " hand-back: the term adds no jump (tau " << dtau << " vs "
                  << dtau0 << " Pa, sigma_n " << dsig << " vs " << dsig0 << " Pa)");
      {
         // control: W0 handed back to the QD regime with the term
         real_t ct = 0.0, cs = 0.0, cm = 0.0;
         u5::Enter(dyn0, b);
         u5::Handback(dyn0, *S.qd, cfg, n_dyn, comm, ct, cs, cm);
         if (g_rank == 0)
         {
            std::cout << std::scientific << std::setprecision(3) << "  control (W0 -> QD with the "
                      << "term): max |tau_qd - tau_dyn| " << ct << " Pa, |sigma_qd - sigma_dyn| "
                      << cs << " Pa\n";
         }
         TEST_ASSERT(std::max(ct, cs) > 10.0 * std::max(dtau, dsig), "U5 " << sub
                     << " control: without the term in the dynamic regime the hand-back jumps "
                     << std::max(ct, cs) / std::max(dtau, dsig) << " x more (> 10)");
      }

      // the deck's faces = "all" (frozen wall load)
      {
         cycle::CycleConfig cfg_all = cfg;
         cfg_all.absorbing.faces = "all";
         cycle::DynPhase dyna(pmesh, cfg_all, S.qd->Table(), *S.par, *S.pre, 1.0, nullptr,
                              S.qd->Domain().FrameField(), Upl);
         u5::Enter(dyna, b);
         NitscheWaveOperator &Wa = dyna.Op();
         const std::vector<char> fma = u5::FaultElementMask(Wa, xmin, xmax);
         Vector Aa;
         Wa.Mult(b.u, zero, Aa);
         const real_t aa = u5::MaxAbsOff(Aa, fma, comm);
         TEST_ASSERT(aa <= 1e-8 * aK, "U5 " << sub << " faces \"all\": A0 off the fault elements "
                     << aa / aK << " a_K (<= 1e-8), the walls' frozen load included");
      }
   }
#endif
}

// ====================================================================== R
static int Refuse(int which, MPIContext &mpi, const char *deck)
{
   MPI_Comm comm = mpi.GetComm();
   LinearElastic model(kLam, kMu);
   BoundaryConfig bc;
   bc.fault_attr = 3;
   bc.dirichlet_attrs = {5};
   bc.default_dirichlet_func = [](const Vector &, real_t, Vector &u) { u.SetSize(2); u = 0.0; };
   if (which == 1 || which == 2)
   {
      std::unique_ptr<Mesh> m(Box(2, 1.0e3, 5));
      ParMesh pm(comm, *m);
      Op op(pm, 1, model, 0.0, 1e4, 1e4, bc, DGMethod::IP, SolverType::MUMPS);
      Vector z(op.GetFESpace().GetVSize() + (which == 2 ? 1 : 0));
      z = 0.0;
      op.SetFiniteWidth(which == 1 ? 0.0 : BetaOverMu(), z, nullptr);
   }
   else if (which == 3)
   {
      Mesh m = Mesh::MakeCartesian3D(1, 1, 1, Element::TETRAHEDRON, 1e3, 1e3, 1e3);
      ParMesh pm(comm, m);
      BoundaryConfig b3;
      b3.fault_attr = 99;
      b3.dirichlet_attrs = {1, 2, 3, 4, 5, 6};
      b3.default_dirichlet_func = [](const Vector &, real_t, Vector &u) { u.SetSize(3); u = 0.0; };
      Op op(pm, 1, model, 0.0, 1e4, 1e4, b3, DGMethod::IP, SolverType::MUMPS);
      Vector z(op.GetFESpace().GetVSize());
      z = 0.0;
      op.SetFiniteWidth(BetaOverMu(), z, nullptr);
   }
   else if (which == 4)
   {
      std::unique_ptr<Mesh> m(Box(2, 1.0e3, 1));
      ParMesh pm(comm, *m);
      auto op = FreeBoxOp(pm, 1);
      Vector x(op->Height()), Kx;
      x = 1.0;
      op->StiffnessMult(x, Kx);   // the dense blocks read K
      Vector z(op->Height());
      z = 0.0;
      op->SetFiniteWidth(BetaOverMu(), z, nullptr);
   }
   else if (which == 5)
   {
#ifdef MFEM_USE_PETSC
      cycle::CycleConfig cfg = cycle::CycleConfig::Load(deck);
      cfg.finite_width.W = 15.0e3;
      cfg.finite_width.c = CWa();
      cfg.finite_width.substrate = "loading_plane";
      std::unique_ptr<Mesh> smesh = cycle::LoadScaledMesh(cfg.shared.mesh.path, cfg.mesh_scale);
      ParMesh pmesh(comm, *smesh);
      smesh.reset();
      u5::Setup S;
      u5::Build(pmesh, mpi, cfg, S);
      cycle::DynPhase dyn(pmesh, cfg, S.qd->Table(), *S.par, *S.pre, 1.0, nullptr,
                          S.qd->Domain().FrameField(), nullptr);
#endif
   }
   if (g_rank == 0) { std::cout << "refusal " << which << ": NOT refused\n"; }
   return 0;
}

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);
   g_rank = mpi.Rank();
   MPI_Comm comm = mpi.GetComm();
   const char *deck = "config/fd_cycle/mapview_rs_tensor_fd.toml";
   const char *only = "";
   int n_dyn = 20, refuse = 0;
   OptionsParser args(argc, argv);
   args.AddOption(&deck, "-d", "--deck", "Map-view deck for U5 (without [cycle.finite_width]).");
   args.AddOption(&only, "-only", "--only", "Comma-separated parts (U1,...,U5, U7, U8); empty = all.");
   args.AddOption(&n_dyn, "-nd", "--dyn-steps", "U5: dynamic steps before the hand-back.");
   args.AddOption(&refuse, "-r", "--refuse", "Run refusal N (1-5); must abort.");
   args.Parse();
   if (!args.Good())
   {
      if (mpi.IsRoot()) { args.PrintUsage(std::cout); }
      return 1;
   }
   std::set<std::string> parts;
   {
      std::stringstream ss(only);
      std::string t;
      while (std::getline(ss, t, ',')) { if (!t.empty()) { parts.insert(t); } }
   }
   auto want = [&](const char *p) { return parts.empty() || parts.count(p) > 0; };
   const bool petsc = (refuse == 5) || (refuse == 0 && want("U5"));
#ifdef MFEM_USE_PETSC
   if (petsc)
   {
      const cycle::CycleConfig c = cycle::CycleConfig::Load(deck);
      MFEMInitializePetsc(&argc, &argv, c.qd.petsc_options.c_str(), NULL);
   }
#endif
   int rc = 0;
   if (refuse > 0) { rc = Refuse(refuse, mpi, deck); }
   else
   {
      if (mpi.IsRoot())
      {
         std::cout << "test_2d5_width on " << mpi.Size() << " rank(s): W = " << kW << " m, c = "
                   << CWa() << " (wa2020, nu = " << Nu() << "), l = " << Ell() << " m, commit "
                   << SEAS_CYCLE_COMMIT << "\n";
      }
      if (want("U1")) { TestManufactured(comm); }
      if (want("U2")) { TestLongFault(comm); }
      if (want("U3")) { TestStressDrop(comm); }
      if (want("U4")) { TestDynamic(comm); }
      if (want("U5")) { TestRegimeSwitch(mpi, deck, n_dyn); }
      if (want("U7")) { TestNormalStiffness(comm); }
      if (want("U8")) { TestPlateSubstrate(comm); }
      if (mpi.IsRoot())
      {
         std::cout << "\ntest_2d5_width: " << (num_tests - num_failed) << " / " << num_tests
                   << " passed" << (num_failed ? " -- FAILED" : "") << "\n";
      }
      rc = num_failed ? 1 : 0;
   }
#ifdef MFEM_USE_PETSC
   if (petsc) { MFEMFinalizePetsc(); }
#endif
   return rc;
}

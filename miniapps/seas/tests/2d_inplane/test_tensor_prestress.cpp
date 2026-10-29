// tests/2d_inplane/test_tensor_prestress.cpp — the cycle's tensor prestress,
// [stress] kind = "constant_tensor" (Chunhui Zhao, 2026-10-29: "we should
// keep the function the same as in dynamic phase, to make it consistent and
// work for 2D and 3D"; "test on the rotated planar fault line (modified from
// TPV205-2D) in a certain angle").
//
// The deck path of seas_fd_cycle_driver on meshes built here: CycleConfig::
// Load, QDPhase, BuildFaultParams, QDPhase::ProjectTensor (the dynamic
// operator's projection, domain/fault_stress_projection.hpp),
// BuildTensorPrestress, QDPhase::Setup, and DynPhase (hook H2).  The
// unrotated layouts: the fault in the plane y = 0 with ref_normal (0, -1, 0),
// so n0 = (0, -1, 0); every case is that layout turned rigidly by R (mesh
// nodes, deck vectors, and the tensor S = R S0 R^T).
//   2D: TPV205-2D's map view, the square |x|, |y| <= 20 km in 16 x 16 cells
//       (two triangles each), the fault y = 0, |x| < 10 km with two free
//       tips, the sides attr 5; up = [0, 0, 1] (the strike slot).  S0:
//       sigma_xx -30 MPa, sigma_yy -25 MPa, sigma_xy +19.490040931338385 MPa
//       (the steady state at V_p of the velocity-weakening zone below), so
//       sigma_n = 25 MPa and the shear 19.49 MPa, right lateral.  a = 0.004
//       for |x| <= 5 km, rising to 0.04 at 7 km.  R: 0, 30 and 135 degrees
//       about z.
//   3D: tests/2d_inplane/frame_box.hpp's box (kind 0), its band (attr 6)
//       turned into prescribed-slip faces (attr 5); up = R e_z.  S0:
//       sigma_xx -30, sigma_yy -25, sigma_zz -28, sigma_xy +10, sigma_yz -6,
//       sigma_xz +2 MPa, so on y = 0 T = S0 n0 = (-10, 25, 6) MPa, sigma_n =
//       25 MPa and the shear (-10, 0, 6) MPa in both slots.  R: I and 40
//       degrees about (1, 2, 3).
//
//   T1  every case at p = 1 and 2 (ic = "scec_fd", switching on): IcLabel
//       has "+tensor"; at every owned fault DOF sigma0 = -n.S.n and |tau0| =
//       |S n - (n.S n) n| (1e-12 rel.) with n = R n0; V_init = -|V| tau0 /
//       |tau0| with |V| = 1e-9 m/s (1e-12); psi0 = a ln[(2 V0 / V)
//       sinh(|tau0| / (sigma0 a))] written out here (1e-12 abs.); in 2D the
//       out-of-plane (dip) slot of tau0 and V_init is exactly 0; at R = I in
//       2D (sigma_xy > 0) tau0's strike slot is < 0 and V_init's > 0 (right
//       lateral, as config/fd_cycle/mapview_rs_fd.toml)
//   T2  after QDPhase::Setup: the state's psi is psi0 bit for bit; the
//       kernel's slip rate is along V_init and solves |tau0| - eta |V| =
//       sigma0 a asinh(|V| / (2 V0) exp(psi0 / a)) with the DOF's own sigma0
//       (1e-12 rel.), so the kernel runs on the per-DOF sigma_n and tau0
//       (SetSAFSMode); QDPhase::SigmaN = sigma0 (1e-12 rel.)
//   T3  QD against dynamic: DynPhase's background traction (hook H2:
//       FrameMap::TractionToGlobal(tau0, sigma0)) is S n with the dynamic
//       operator's oriented normal n, which is R n0 (1e-12 rel.); and it
//       equals the background traction of a NitscheWaveOperator built with
//       fs.prestress = S (the dynamic phase's own road, as the Nitsche driver
//       for TPV205-2D) to 1e-14 rel.
//   T4  rotation: each rotated case gives the R = I case's values at R x
//       (kernel tau0 and V_init, sigma0, psi0) to 1e-12 rel.
//   T5  ic = "gold" (switching off, eta_V_init = true) at 30 degrees in 2D,
//       p = 1 and 2: the kernel's own gold inversion on the per-DOF values
//       gives BuildTensorPrestress's psi0 (1e-12) and V = V_init (1e-9 rel.)
//   On 2 ranks the two sides of the fault on different ranks (every fault
//   face a rank-boundary face).
//   T6  2D at R = I against [stress] kind = "fault_local_prestress" (the
//       steady state, BuildPrestress, sigma_n 25 MPa): at the DOFs with a =
//       0.004 the tensor's tau0 and psi0 are the steady state's (1e-12 rel.)
//
// Modes: (default) T1-T6;  --refuse N (N = 1..7) must abort with:
//   1 sigma_yy = +25 MPa (the fault in tension)  "the fault must be in compression"
//   2 sigma_xy = 0 in 2D (no shear)              "the tensor drives no slip there"
//   3 [cycle.prestress] a = "amax"               "[cycle.prestress].a = \"amax\" with [stress] kind = \"constant_tensor\""
//   4 [friction.rate_state].sigma_n_default      "[friction.rate_state].sigma_n_default with [stress] kind = \"constant_tensor\""
//   5 a [[stress_region]] table                  "[[stress_region]] with [stress] kind = \"constant_tensor\""
//   6 sigma_xz = 1 MPa in 2D (anti-plane)        "the in-plane model carries no anti-plane stress"
//   7 eta_V_init = true, sigma_xy = 1e-3 Pa      "no state balances it at V_init"
//
//   make -f tests/2d_inplane/2d_inplane.mk test-tensor-prestress
//
// Scope (Chunhui Zhao, 2026-10-29).

#include "mfem.hpp"
#include "frame_box.hpp"
#include "../../cycle/cycle_config.hpp"
#include "../../cycle/cycle_frames.hpp"
#include "../../cycle/fault_fields.hpp"
#include "../../cycle/prestress.hpp"
#include "../../cycle/qd_phase.hpp"
#include "../../cycle/dyn_phase.hpp"
#include "../../cycle/state_bridge.hpp"

#include <unistd.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;
using frame_box::Rot;

static int num_tests = 0, num_failed = 0;
static int g_rank = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { if (g_rank == 0) std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; if (g_rank == 0) std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

// The friction (BP5's scalars, as config/fd_cycle/mapview_rs_fd.toml).
static const real_t kRho = 2670.0, kMu = 32038120320.0;
static const real_t kB = 0.03, kV0 = 1.0e-6, kF0 = 0.6, kVp = 1.0e-9;
static const real_t kAvw = 0.004, kAvs = 0.04;
static const real_t kL2 = 20.0e3;              // 2D box half-width [m]
static const real_t kTauSS = 19490040.931338385;   // sigma_n a asinh(...) at V_p, a = 0.004

static real_t Eta() { return 0.5 * std::sqrt(kMu * kRho); }   // "auto": mu / (2 cs)

// ------------------------------------------------------------------ tensors
using Tensor = std::array<real_t, 9>;   // row major, tension positive

static Tensor Sym(real_t xx, real_t yy, real_t zz, real_t xy, real_t yz, real_t xz)
{
   return {{xx, xy, xz, xy, yy, yz, xz, yz, zz}};
}

/// R S R^T.
static Tensor Rotate(const Rot &r, const Tensor &S)
{
   Tensor out{};
   for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
      {
         real_t v = 0.0;
         for (int k = 0; k < 3; k++)
            for (int l = 0; l < 3; l++) { v += r.R[i][k] * S[3 * k + l] * r.R[j][l]; }
         out[3 * i + j] = v;
      }
   return out;
}

static void MatVec(const Tensor &S, const real_t n[3], real_t T[3])
{
   for (int i = 0; i < 3; i++) { T[i] = S[3 * i] * n[0] + S[3 * i + 1] * n[1] + S[3 * i + 2] * n[2]; }
}

static real_t Dot(const real_t a[3], const real_t b[3]) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }

static Tensor Base2D() { return Sym(-30.0e6, -25.0e6, -25.0e6, kTauSS, 0.0, 0.0); }
static Tensor Base3D() { return Sym(-30.0e6, -25.0e6, -28.0e6, 10.0e6, -6.0e6, 2.0e6); }

// ------------------------------------------------------------------- meshes
/// 2D: the square |x|, |y| <= kL2 in n x n cells (diagonal (i, j)-(i+1, j+1)),
/// the fault y = 0, |x| < kL2 / 2 (attr 3), free tips, the sides attr 5; the
/// vertices turned by R (about z) after tagging.
static Mesh *MapLine(const Rot &r)
{
   const int n = 16;
   Mesh *m = new Mesh(2, (n + 1) * (n + 1), 2 * n * n, 0, 2);
   for (int j = 0; j <= n; j++)
      for (int i = 0; i <= n; i++) { m->AddVertex(-kL2 + 2.0 * kL2 * i / n, -kL2 + 2.0 * kL2 * j / n); }
   auto v = [&](int i, int j) { return i + (n + 1) * j; };
   for (int j = 0; j < n; j++)
      for (int i = 0; i < n; i++)
      {
         m->AddTriangle(v(i, j), v(i + 1, j), v(i + 1, j + 1), 1);
         m->AddTriangle(v(i, j), v(i + 1, j + 1), v(i, j + 1), 1);
      }
   for (int k = 0; k < n; k++)
   {
      m->AddBdrSegment(v(k, 0), v(k + 1, 0), 5);
      m->AddBdrSegment(v(k, n), v(k + 1, n), 5);
      m->AddBdrSegment(v(0, k), v(0, k + 1), 5);
      m->AddBdrSegment(v(n, k), v(n, k + 1), 5);
      const real_t xm = -kL2 + 2.0 * kL2 * (k + 0.5) / n;
      if (std::abs(xm) < 0.5 * kL2) { m->AddBdrSegment(v(k, n / 2), v(k + 1, n / 2), 3); }
   }
   m->FinalizeTopology();
   m->Finalize(false, true);
   m->SetAttributes();
   for (int i = 0; i < m->GetNV(); i++)
   {
      real_t *c = m->GetVertex(i);
      const real_t x[3] = {c[0], c[1], 0.0};
      real_t y[3];
      r.Apply(x, y);
      c[0] = y[0];
      c[1] = y[1];
   }
   return m;
}

/// 3D: FrameBox(0, R) with its band (attr 6) as prescribed-slip faces (attr 5).
static Mesh *Box3D(const Rot &r)
{
   Mesh *m = frame_box::FrameBox(0, r);
   for (int b = 0; b < m->GetNBE(); b++)
   {
      if (m->GetBdrAttribute(b) == 6) { m->SetBdrAttribute(b, 5); }
   }
   m->SetAttributes();
   return m;
}

// -------------------------------------------------------------------- decks
struct Spec
{
   int dim = 2;
   int p = 1;
   Rot R;
   Tensor S{};                  ///< the tensor in the global (rotated) frame
   std::string ic = "scec_fd";
   bool eta_V_init = false;
   bool steady = false;         ///< T6: [stress] fault_local_prestress instead
   std::string prestress_a;     ///< "" = local (default)
   std::string rs_extra;        ///< lines appended to [friction.rate_state]
   std::string tail;            ///< text appended to the deck
};

static std::string Num(real_t v)
{
   std::ostringstream o;
   o << std::setprecision(17) << v;
   std::string s = o.str();
   if (s.find_first_of(".eE") == std::string::npos) { s += ".0"; }
   return s;
}

/// R v as a deck vector of dim entries.
static std::string DeckVec(const Rot &r, real_t x, real_t y, real_t z, int dim)
{
   const real_t a[3] = {x, y, z};
   real_t b[3];
   r.Apply(a, b);
   std::string s = "[" + Num(b[0]) + ", " + Num(b[1]);
   if (dim == 3) { s += ", " + Num(b[2]); }
   return s + "]";
}

/// The rectangle-taper field: 2D one axis, 3D two (along strike x and down
/// dip z of the unrotated plane).
static std::string Field(const Spec &sp, const std::string &name, real_t half, real_t taper,
                         real_t inside, real_t outside)
{
   std::ostringstream o;
   o << "[cycle.fields." << name << "]\n"
     << "center  = " << DeckVec(sp.R, 0.0, 0.0, 0.0, sp.dim) << "\n"
     << "axis1   = " << DeckVec(sp.R, 1.0, 0.0, 0.0, sp.dim) << "\n"
     << "half1   = " << Num(half) << "\n";
   if (sp.dim == 3)
   {
      o << "axis2   = " << DeckVec(sp.R, 0.0, 0.0, 1.0, 3) << "\n"
        << "half2   = " << Num(half) << "\n";
   }
   if (taper > 0.0) { o << "taper   = " << Num(taper) << "\n"; }
   o << "inside  = " << Num(inside) << "\n"
     << "outside = " << Num(outside) << "\n\n";
   return o.str();
}

static std::string DeckText(const Spec &sp)
{
   const int d = sp.dim;
   const bool gold = (sp.ic == "gold");
   std::ostringstream o;
   o << "[meta]\nschema_version = 1\nlaw = \"rate_state\"\n"
     << "description = \"test_tensor_prestress\"\n\n"
     << "[material_constant_fallback]\nlambda = " << Num(kMu) << "\nmu = " << Num(kMu)
     << "\nrho = " << Num(kRho) << "\n\n"
     << "[pore_pressure]\nP_p_pa = 0.0\nP_p_grad_pa_per_m = 0.0\nmin_sigma_n_pa = 0.0\n\n"
     << "[mesh]\npath = \"built_in_the_test.msh\"\norder = " << sp.p << "\n\n"
     << "[velocity]\nuse_sidecar = false\n\n";
   if (sp.steady)
   {
      o << "[stress]\nkind = \"fault_local_prestress\"\ntau_strike_pa = 0.0\n"
        << "tau_dip_pa = 0.0\nsigma_n_pa = 25.0e6\n\n";
   }
   else
   {
      const Tensor &S = sp.S;
      o << "[stress]\nkind = \"constant_tensor\"\n"
        << "sigma_xx_pa = " << Num(S[0]) << "\nsigma_yy_pa = " << Num(S[4])
        << "\nsigma_zz_pa = " << Num(S[8]) << "\nsigma_xy_pa = " << Num(S[1])
        << "\nsigma_yz_pa = " << Num(S[5]) << "\nsigma_xz_pa = " << Num(S[2]) << "\n\n";
   }
   o << "[numerics]\ncfl = 0.5\nkv_eta = 0.1\nkv_scale = \"node\"\nkv_length = \"inscribed\"\n\n"
     << "[time]\ntfinal = 1.0e9\ndt_initial = 0.01\n\n"
     << "[output]\noutput_dir = \".\"\n\n"
     << "[boundary]\nfault_attr = 3\nnatural_attrs = [1]\n\n"
     << "[fault_geometry]\nref_normal = " << DeckVec(sp.R, 0.0, -1.0, 0.0, d) << "\n"
     << "up = " << (d == 2 ? std::string("[0.0, 0.0, 1.0]") : DeckVec(sp.R, 0.0, 0.0, 1.0, 3)) << "\n"
     << "kind = \"vertical_strike_slip\"\norientation = \"propagated\"\npin_by = \"reference\"\n"
     << "tangent = \"pointwise\"\njunction_policy = \"refuse\"\npin_min_margin = 0.1\n\n"
     << "[friction.rate_state]\na_default = 0.004\nb_default = 0.03\nDc_default = 0.14\n"
     << "V_0_default = 1.0e-6\nf_0_default = 0.6\nV_init_default = 1.0e-9\n"
     << (sp.steady ? "sigma_n_default = 25.0e6\n" : "")
     << "state_evolution = \"aging_law\"\neta = \"auto\"\n" << sp.rs_extra << "\n"
     << "[cycle]\nic = \"" << sp.ic << "\"\nswitching = " << (gold ? "false" : "true") << "\n"
     << "qd_radiation_damping = true\nmesh_scale = 1.0\nplate_rate = 1.0e-9\n"
     << "penalty_factor = 1.0\n\n"
     << "[cycle.qd]\nsolver = \"mumps\"\n"
     << "petsc_options = \"tests/verification/petsc_ts_rk45_tandem.cfg\"\n"
     << "ts_atol = 1.0e-7\nts_rtol = 1.0e-50\n\n"
     << "[cycle.prestress]\nV_background = 1.0e-9\neta_V_init = "
     << (sp.eta_V_init ? "true" : "false") << "\n"
     << (sp.prestress_a.empty() ? "" : "a = \"" + sp.prestress_a + "\"\n") << "\n"
     << "[cycle.loading]\nkind = \"plate\"\nplane_point = " << DeckVec(sp.R, 0.0, 0.0, 0.0, d) << "\n"
     << "plane_normal = " << DeckVec(sp.R, 0.0, 1.0, 0.0, d) << "\n"
     << "slip_dir = " << DeckVec(sp.R, 1.0, 0.0, 0.0, d) << "\n"
     << "band = " << (d == 2 ? "0.0" : "1.0") << "\n\n"
     << "[cycle.fields]\nV_init_dip = 0.0\n\n"
     << Field(sp, "a", 5.0e3, 2.0e3, kAvw, kAvs)
     << Field(sp, "Dc", 0.0, 0.0, 0.14, 0.14)
     << Field(sp, "V_init", 0.0, 0.0, kVp, kVp)
     << "[cycle.output]\nprefix = \"t\"\nplane_origin = " << DeckVec(sp.R, 0.0, 0.0, 0.0, d) << "\n"
     << "plane_axis1 = " << DeckVec(sp.R, 1.0, 0.0, 0.0, d) << "\n"
     << "plane_axis2 = " << DeckVec(sp.R, 0.0, d == 2 ? 1.0 : 0.0, d == 2 ? 0.0 : 1.0, d) << "\n"
     << sp.tail;
   return o.str();
}

static std::string WriteDeck(const std::string &tag, const std::string &text)
{
   const char *tmp = std::getenv("TMPDIR");
   const std::string path = std::string(tmp ? tmp : "/tmp") + "/seas_test_tensor_prestress_" + tag
                            + "_" + std::to_string(getpid()) + "_" + std::to_string(g_rank) + ".toml";
   std::ofstream f(path);
   f << text;
   MFEM_VERIFY(f, "cannot write " << path);
   return path;
}

static cycle::CycleConfig LoadSpec(const std::string &tag, const Spec &sp)
{
   const std::string path = WriteDeck(tag, DeckText(sp));
   cycle::CycleConfig cfg = cycle::CycleConfig::Load(path);
   std::remove(path.c_str());
   return cfg;
}

// ----------------------------------------------------------------- one case
/// Per owned DOF, gathered on rank 0: x (3), tau0 (2), sigma0, psi0, V_init
/// (2), a.
using Row = std::array<real_t, 10>;

static std::vector<Row> Gather(const std::vector<Row> &mine, MPI_Comm comm)
{
   int size = 1;
   MPI_Comm_size(comm, &size);
   const int n = 10 * static_cast<int>(mine.size());
   std::vector<int> counts(size), displ(size, 0);
   MPI_Gather(&n, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, comm);
   std::vector<real_t> flat;
   for (const Row &r : mine) { flat.insert(flat.end(), r.begin(), r.end()); }
   int total = 0;
   if (g_rank == 0)
   {
      for (int i = 0; i < size; i++) { displ[i] = total; total += counts[i]; }
   }
   std::vector<real_t> all(std::max(total, 1));
   MPI_Gatherv(flat.data(), n, MPI_DOUBLE, all.data(), counts.data(), displ.data(), MPI_DOUBLE,
               0, comm);
   std::vector<Row> out;
   for (int k = 0; k + 10 <= total; k += 10)
   {
      Row r;
      std::copy(all.begin() + k, all.begin() + k + 10, r.begin());
      out.push_back(r);
   }
   return out;
}

static std::string CaseName(const Spec &sp, const std::string &rot)
{
   std::ostringstream o;
   o << sp.dim << "D " << rot << " p = " << sp.p << (sp.ic == "gold" ? " gold" : "")
     << (sp.steady ? " steady state" : "") << ": ";
   return o.str();
}

/// T1-T3 (switching) or T5 (gold) on one spec; returns the rows (rank 0).
static std::vector<Row> RunCase(const Spec &sp, const std::string &rot, MPIContext &mpi)
{
   const std::string tag = CaseName(sp, rot);
   cycle::CycleConfig cfg = LoadSpec("case", sp);
   std::unique_ptr<Mesh> smesh(sp.dim == 2 ? MapLine(sp.R) : Box3D(sp.R));
   cfg.CheckDimension(smesh->Dimension());
   // On 2 ranks the two sides of the fault plane on different ranks, so that
   // every fault face is a rank-boundary face (the projection's shared faces).
   std::vector<int> part;
   if (mpi.Size() == 2)
   {
      const real_t m0[3] = {0.0, 1.0, 0.0};
      real_t m[3];
      sp.R.Apply(m0, m);
      part.resize(smesh->GetNE());
      for (int e = 0; e < smesh->GetNE(); e++)
      {
         Vector c;
         smesh->GetElementCenter(e, c);
         real_t side = 0.0;
         for (int k = 0; k < c.Size(); k++) { side += c(k) * m[k]; }
         part[e] = (side < 0.0) ? 1 : 0;
      }
   }
   ParMesh pmesh(mpi.GetComm(), *smesh, part.empty() ? nullptr : part.data());
   smesh.reset();
   const bool gold = (sp.ic == "gold");
   if (!sp.steady)
   {
      const std::string want = sp.ic + "+tensor" + (sp.eta_V_init && !gold ? "+eta_V_init" : "");
      TEST_ASSERT(cfg.TensorPrestress() && cfg.IcLabel() == want && std::isnan(cfg.sigma_n),
                  tag << "deck loads: IcLabel \"" << cfg.IcLabel() << "\", no scalar sigma_n");
   }

   cycle::QDPhase qd(pmesh, mpi, cfg);
   const spatial::SpatialFrictionConfig &s = cfg.shared;
   const MaterialField material = MaterialField::MakeConstant(
      s.material_fallback.lambda, s.material_fallback.mu, s.material_fallback.rho);
   const spatial::PorePressureSpec no_pp;
   cycle::FaultParams par = cycle::BuildFaultParams(
      qd.Table(), *s.rate_state, cfg.fields, material, pmesh, no_pp);
   cycle::Prestress pre;
   if (sp.steady) { pre = cycle::BuildPrestress(par, cfg.prestress); }
   else
   {
      Vector K;
      qd.ProjectTensor(cfg.stress_tensor, K);
      pre = cycle::BuildTensorPrestress(par, K, cfg.prestress, qd.Table());
   }
   const int N = qd.Table().Size();
   const int nglob = mpi.GlobalSumInt(N);

   // ---- T1: against the analytic projection ----
   const real_t n0[3] = {0.0, -1.0, 0.0};
   real_t n[3], T[3], tv[3];
   sp.R.Apply(n0, n);
   MatVec(sp.S, n, T);
   const real_t Tn = Dot(T, n), sn_a = -Tn;
   for (int c = 0; c < 3; c++) { tv[c] = T[c] - Tn * n[c]; }
   const real_t tau_a = std::sqrt(Dot(tv, tv));
   real_t e_sn = 0.0, e_tau = 0.0, e_psi = 0.0, e_v = 0.0, e_dir = 0.0, out_slot = 0.0;
   int right_lateral = 0;
   std::vector<Row> rows;
   for (int i = 0; i < N; i++)
   {
      const real_t t1 = pre.tau0(2 * i), t2 = pre.tau0(2 * i + 1);
      const real_t tau = std::sqrt(t1 * t1 + t2 * t2), sn = pre.sigma0(i), a = par.rs.a(i);
      const real_t v1 = par.V_init(2 * i), v2 = par.V_init(2 * i + 1);
      const real_t vabs = std::sqrt(v1 * v1 + v2 * v2);
      const real_t eta = sp.eta_V_init ? par.rs.eta(i) : 0.0;
      const real_t psi = a * std::log(2.0 * kV0 / kVp * std::sinh((tau - eta * kVp) / (sn * a)));
      e_sn = std::max(e_sn, std::abs(sn - sn_a) / sn_a);
      e_tau = std::max(e_tau, std::abs(tau - tau_a) / tau_a);
      e_psi = std::max(e_psi, std::abs(pre.psi0(i) - psi));
      e_v = std::max(e_v, std::abs(vabs - kVp) / kVp);
      // V_init anti-parallel to tau0: V + |V| tau0 / |tau0| = 0
      e_dir = std::max(e_dir, std::hypot(v1 + vabs * t1 / tau, v2 + vabs * t2 / tau) / vabs);
      if (sp.dim == 2) { out_slot = std::max({out_slot, std::abs(t1), std::abs(v1)}); }
      if (t2 < 0.0 && v2 > 0.0) { right_lateral++; }
      const real_t *x = qd.Table().coords.GetData() + 3 * i;
      rows.push_back({{x[0], x[1], x[2], t1, t2, sn, pre.psi0(i), v1, v2, a}});
   }
   e_sn = mpi.GlobalMax(e_sn);
   e_tau = mpi.GlobalMax(e_tau);
   e_psi = mpi.GlobalMax(e_psi);
   e_v = mpi.GlobalMax(e_v);
   e_dir = mpi.GlobalMax(e_dir);
   out_slot = mpi.GlobalMax(out_slot);
   right_lateral = mpi.GlobalSumInt(right_lateral);
   if (!sp.steady)
   {
      TEST_ASSERT(e_sn <= 1e-12 && e_tau <= 1e-12, tag << nglob << " fault DOFs: sigma0 = -n.S.n = "
                  << sn_a / 1e6 << " MPa (rel. " << e_sn << "), |tau0| = |S n - (n.S n) n| = "
                  << tau_a / 1e6 << " MPa (rel. " << e_tau << ")");
      TEST_ASSERT(e_v <= 1e-12 && e_dir <= 1e-12, tag << "V_init = -|V| tau0 / |tau0|, |V| = 1e-9 "
                  "(rel. " << e_v << ", direction " << e_dir << ")");
      TEST_ASSERT(e_psi <= 1e-12, tag << "psi0 = a ln[(2 V0 / V) sinh((|tau0| - eta V) / "
                  "(sigma0 a))] (abs. " << e_psi << ")");
      if (sp.dim == 2)
      {
         TEST_ASSERT(out_slot == 0.0, tag << "the out-of-plane (dip) slot of tau0 and V_init "
                     "is exactly 0");
      }
      if (sp.dim == 2 && rot == "0 deg")
      {
         TEST_ASSERT(right_lateral == nglob, tag << "sigma_xy > 0: tau0 strike < 0 and V_init "
                     "strike > 0 at " << right_lateral << " of " << nglob << " DOFs (right lateral)");
      }
   }

   // ---- T2 / T5: the kernel after Setup ----
   qd.Setup(par, pre);
   Vector psi_k, sn_k;
   qd.Psi(psi_k);
   qd.SigmaN(sn_k);
   const Vector &V = qd.SlipRate();
   if (!sp.steady)
   {
      real_t e_state = 0.0, e_bal = 0.0, e_vd = 0.0, e_snk = 0.0, e_vg = 0.0;
      for (int i = 0; i < N; i++)
      {
         const real_t a = par.rs.a(i), sn = pre.sigma0(i);
         const real_t vi1 = par.V_init(2 * i), vi2 = par.V_init(2 * i + 1);
         const real_t vk = std::hypot(V(2 * i), V(2 * i + 1));
         const real_t tau = std::hypot(pre.tau0(2 * i), pre.tau0(2 * i + 1));
         e_state = std::max(e_state, std::abs(psi_k(i) - pre.psi0(i)));
         e_snk = std::max(e_snk, std::abs(sn_k(i) - sn) / sn);
         if (gold)
         {
            e_vg = std::max(e_vg, std::hypot(V(2 * i) - vi1, V(2 * i + 1) - vi2) / kVp);
         }
         else
         {
            const real_t F = sn * a * std::asinh(vk / (2.0 * kV0) * std::exp(psi_k(i) / a));
            e_bal = std::max(e_bal, std::abs(F + par.rs.eta(i) * vk - tau) / tau);
            // along V_init: V |V_init| - V_init |V| = 0
            e_vd = std::max(e_vd, std::hypot(V(2 * i) * kVp - vi1 * vk,
                                             V(2 * i + 1) * kVp - vi2 * vk) / (kVp * vk));
         }
      }
      e_state = mpi.GlobalMax(e_state);
      e_bal = mpi.GlobalMax(e_bal);
      e_vd = mpi.GlobalMax(e_vd);
      e_snk = mpi.GlobalMax(e_snk);
      e_vg = mpi.GlobalMax(e_vg);
      if (gold)
      {
         TEST_ASSERT(e_state <= 1e-12, tag << "the kernel's gold inversion on the per-DOF "
                     "values = BuildTensorPrestress psi0 (abs. " << e_state << ")");
         TEST_ASSERT(e_vg <= 1e-9, tag << "kernel slip rate at t = 0 = V_init (rel. " << e_vg << ")");
      }
      else
      {
         TEST_ASSERT(e_state == 0.0, tag << "state psi = psi0 bit for bit (" << e_state << ")");
         TEST_ASSERT(e_bal <= 1e-12 && e_vd <= 1e-12, tag << "kernel slip rate at t = 0: |tau0| - "
                     "eta |V| = sigma0 a asinh(|V| / (2 V0) exp(psi0 / a)) with the DOF's sigma0 "
                     "(rel. " << e_bal << "), along V_init (" << e_vd << ")");
      }
      TEST_ASSERT(e_snk <= 1e-12, tag << "QDPhase::SigmaN = sigma0 per DOF (rel. " << e_snk << ")");
   }
   if (gold || sp.steady) { return Gather(rows, mpi.GetComm()); }

   // ---- T3: QD against dynamic ----
   auto frames = qd.Domain().FrameField();
   cycle::DynPhase dyn(pmesh, cfg, qd.Table(), par, pre, 1.0, nullptr, frames);
   const Vector &S0 = dyn.Op().PrestressTraction();
   const Vector &F = dyn.Op().FaultFrames();
   real_t e_n = 0.0, e_T = 0.0;
   const real_t Tm = std::sqrt(Dot(T, T));
   for (int j = 0; j < N; j++)
   {
      const real_t *nk = &F(9 * j);
      real_t Tj[3];
      MatVec(sp.S, nk, Tj);
      e_n = std::max(e_n, std::sqrt((nk[0] - n[0]) * (nk[0] - n[0]) + (nk[1] - n[1]) * (nk[1] - n[1])
                                    + (nk[2] - n[2]) * (nk[2] - n[2])));
      e_T = std::max(e_T, std::sqrt((S0(3 * j) - Tj[0]) * (S0(3 * j) - Tj[0])
                                    + (S0(3 * j + 1) - Tj[1]) * (S0(3 * j + 1) - Tj[1])
                                    + (S0(3 * j + 2) - Tj[2]) * (S0(3 * j + 2) - Tj[2])) / Tm);
   }
   e_n = mpi.GlobalMax(e_n);
   e_T = mpi.GlobalMax(e_T);
   TEST_ASSERT(e_n <= 1e-12 && e_T <= 1e-12, tag << "H2 background traction = S n, n the dynamic "
               "operator's oriented normal = R n0 (" << e_n << "), rel. " << e_T);
   if (mpi.Size() == 2)
   {
      TEST_ASSERT(dyn.NumSharedFaultFaces() > 0, tag << dyn.NumSharedFaultFaces()
                  << " fault faces on the rank boundary");
   }

   // the dynamic phase's own road: the operator built with fs.prestress = S
   const spatial::RateStateBlock &rsb = *s.rate_state;
   DieterichRuinaFriction fr(DieterichRuinaFriction::Constants{
      rsb.V_0_default, rsb.f_0_default, rsb.b_default, rsb.Dc_default});
   RateStateClosure closure(fr, RateStateClosure::StateLaw::Aging);
   NitscheFaultSetup fs;
   fs.closure = &closure;
   {
      const FaultFrameConfig fc = cycle::CycleFrameConfig(cfg);
      fs.orientation = fc.orientation;
      for (int i = 0; i < 3; i++) { fs.up[i] = fc.up[i]; }
   }
   fs.sigma_min = cfg.dyn.sigma_min;
   const Tensor S = cfg.stress_tensor;
   fs.prestress = [S](const real_t[3], real_t sig[9]) { for (int i = 0; i < 9; i++) { sig[i] = S[i]; } };
   std::map<std::array<real_t, 3>, int> where;
   for (int i = 0; i < N; i++)
   {
      const real_t *x = qd.Table().coords.GetData() + 3 * i;
      where[{x[0], x[1], x[2]}] = i;
   }
   const Vector a_dof = par.rs.a, Dc_dof = par.rs.Dc;
   fs.parameters = [where, a_dof, Dc_dof](const real_t x[3], FaultPointState &st)
   {
      const auto it = where.find({x[0], x[1], x[2]});
      MFEM_VERIFY(it != where.end(), "test_tensor_prestress: fault DOF not in the table");
      st.a = a_dof(it->second);
      st.Dc = Dc_dof(it->second);
   };
   BoundaryConfig bc = dyn.Boundaries();
   DomainConfig dcfg;
   dcfg.penalty_factor = cfg.penalty_factor;
   dcfg.fault_frames.field = frames;
   NitscheWaveOperator W(pmesh, sp.p, material, bc, dcfg, &fs);
   MFEM_VERIFY(W.NumOwnedFaultDOFs() == N, "test_tensor_prestress: owned fault DOFs differ");
   const Vector &Sw = W.PrestressTraction();
   real_t e_w = 0.0, e_pos = 0.0;
   for (int k = 0; k < 3 * N; k++)
   {
      e_w = std::max(e_w, std::abs(Sw(k) - S0(k)) / Tm);
      e_pos = std::max(e_pos, std::abs(W.FaultPositions()(k) - dyn.FaultPositions()(k)));
   }
   e_w = mpi.GlobalMax(e_w);
   e_pos = mpi.GlobalMax(e_pos);
   TEST_ASSERT(e_pos == 0.0 && e_w <= 1e-14, tag << "QD tau0, sigma0 through H2 = the Nitsche "
               "operator's own projection of S (fs.prestress), rel. " << e_w);
   return Gather(rows, mpi.GetComm());
}

// ------------------------------------------------------------------ T4, T6
/// B's rows matched to A's by position (R^T x_B against x_A, within 1e-6 m);
/// the largest differences of (tau0 pair, sigma0, psi0, V_init pair)
/// relative to A's largest |value| of each.  Rank 0.
static void Compare(const std::vector<Row> &A, const std::vector<Row> &B, const Rot &r,
                    real_t &e_tau, real_t &e_sn, real_t &e_psi, real_t &e_v, int &matched,
                    bool vw_only = false)
{
   e_tau = e_sn = e_psi = e_v = 0.0;
   matched = 0;
   real_t s_tau = 0.0, s_sn = 0.0, s_psi = 0.0;
   for (const Row &a : A)
   {
      s_tau = std::max(s_tau, std::hypot(a[3], a[4]));
      s_sn = std::max(s_sn, std::abs(a[5]));
      s_psi = std::max(s_psi, std::abs(a[6]));
   }
   for (const Row &a : A)
   {
      if (vw_only && a[9] != kAvw) { continue; }
      int best = -1;
      real_t dbest = 1e300;
      for (size_t k = 0; k < B.size(); k++)
      {
         const real_t xb[3] = {B[k][0], B[k][1], B[k][2]};
         real_t xr[3] = {0, 0, 0};   // R^T x_B
         for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) { xr[i] += r.R[j][i] * xb[j]; }
         const real_t d = std::max({std::abs(xr[0] - a[0]), std::abs(xr[1] - a[1]),
                                    std::abs(xr[2] - a[2])});
         if (d <= 1e-6 && d < dbest) { dbest = d; best = static_cast<int>(k); }
      }
      if (best < 0) { e_tau = e_sn = e_psi = e_v = 1e300; return; }
      const Row &b = B[best];
      matched++;
      e_tau = std::max(e_tau, std::hypot(b[3] - a[3], b[4] - a[4]) / s_tau);
      e_sn = std::max(e_sn, std::abs(b[5] - a[5]) / s_sn);
      e_psi = std::max(e_psi, std::abs(b[6] - a[6]) / s_psi);
      e_v = std::max(e_v, std::hypot(b[7] - a[7], b[8] - a[8]) / kVp);
   }
}

static void Report(const std::string &what, const std::vector<Row> &A, const std::vector<Row> &B,
                   const Rot &r, bool vw_only, MPIContext &mpi)
{
   real_t e[4] = {0, 0, 0, 0};
   int m = 0;
   if (g_rank == 0) { Compare(A, B, r, e[0], e[1], e[2], e[3], m, vw_only); }
   MPI_Bcast(e, 4, MPI_DOUBLE, 0, mpi.GetComm());
   MPI_Bcast(&m, 1, MPI_INT, 0, mpi.GetComm());
   int sizes[2] = {static_cast<int>(A.size()), static_cast<int>(B.size())};
   MPI_Bcast(sizes, 2, MPI_INT, 0, mpi.GetComm());
   // every DOF of A matched (T6: every a = 0.004 DOF), B no larger
   const bool all = sizes[0] == sizes[1] && (vw_only ? m > 0 : m == sizes[0]);
   TEST_ASSERT(all && e[0] <= 1e-12 && e[1] <= 1e-12 && e[2] <= 1e-12 && e[3] <= 1e-12,
               what << ": " << m << " DOFs matched by position: tau0 " << e[0] << ", sigma0 "
               << e[1] << ", psi0 " << e[2] << ", V_init " << e[3] << " (rel.)");
}

// --------------------------------------------------------------- refusals
static int Refuse(int n, MPIContext &mpi)
{
   Spec sp;
   sp.S = Base2D();
   switch (n)
   {
      case 1: sp.S[4] = 25.0e6; break;
      case 2: sp.S[1] = sp.S[3] = 0.0; break;
      case 3: sp.prestress_a = "amax"; break;
      case 4: sp.rs_extra = "sigma_n_default = 25.0e6\n"; break;
      case 5:
         sp.tail = "\n[[stress_region]]\nx_min_m = -1500.0\nx_max_m = 1500.0\n"
                   "sigma_xx_pa = 0.0\nsigma_yy_pa = -25.0e6\nsigma_zz_pa = 0.0\n"
                   "sigma_xy_pa = 20.0e6\nsigma_yz_pa = 0.0\nsigma_xz_pa = 0.0\n";
         break;
      case 6: sp.S[2] = sp.S[6] = 1.0e6; break;
      case 7: sp.eta_V_init = true; sp.S[1] = sp.S[3] = 1.0e-3; break;
      default: MFEM_ABORT("--refuse N: N = 1..7");
   }
   if (g_rank == 0) { std::cout << "--refuse " << n << ": loading the deck\n"; }
   cycle::CycleConfig cfg = LoadSpec("refuse" + std::to_string(n), sp);   // 3, 4, 5 abort here
   std::unique_ptr<Mesh> smesh(MapLine(sp.R));
   ParMesh pmesh(mpi.GetComm(), *smesh);
   smesh.reset();
   cycle::QDPhase qd(pmesh, mpi, cfg);
   const spatial::SpatialFrictionConfig &s = cfg.shared;
   const MaterialField material = MaterialField::MakeConstant(
      s.material_fallback.lambda, s.material_fallback.mu, s.material_fallback.rho);
   const spatial::PorePressureSpec no_pp;
   cycle::FaultParams par = cycle::BuildFaultParams(
      qd.Table(), *s.rate_state, cfg.fields, material, pmesh, no_pp);
   Vector K;
   if (g_rank == 0) { std::cout << "--refuse " << n << ": projecting the tensor\n"; }
   qd.ProjectTensor(cfg.stress_tensor, K);                                 // 6 aborts here
   if (g_rank == 0) { std::cout << "--refuse " << n << ": BuildTensorPrestress\n"; }
   cycle::BuildTensorPrestress(par, K, cfg.prestress, qd.Table());        // 1, 2, 7 abort here
   if (g_rank == 0) { std::cout << "--refuse " << n << ": accepted\n"; }
   return 0;
}

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);
   g_rank = mpi.Rank();
   int refuse = 0;
   for (int i = 1; i < argc; i++)
   {
      if (std::strcmp(argv[i], "--refuse") == 0 && i + 1 < argc) { refuse = std::atoi(argv[++i]); }
   }
   MFEMInitializePetsc(&argc, &argv, nullptr, nullptr);
   if (refuse)
   {
      const int rc = Refuse(refuse, mpi);
      MFEMFinalizePetsc();
      return rc;
   }
   if (g_rank == 0)
   {
      std::cout << "test_tensor_prestress ([stress] kind = \"constant_tensor\" in the cycle) on "
                << mpi.Size() << " rank(s)\n";
   }
   for (int p : {1, 2})
   {
      // 2D map view, the rotated line
      if (g_rank == 0) { std::cout << "\n=== 2D map view, p = " << p << " ===\n"; }
      std::vector<Row> ref;
      for (real_t deg : {0.0, 30.0, 135.0})
      {
         Spec sp;
         sp.p = p;
         sp.R = Rot::Axis(0, 0, 1, deg);
         sp.S = Rotate(sp.R, Base2D());
         std::ostringstream rot;
         rot << deg << " deg";
         const std::vector<Row> rows = RunCase(sp, rot.str(), mpi);
         if (deg == 0.0) { ref = rows; }
         else { Report("T4 2D p = " + std::to_string(p) + " " + rot.str() + " against 0 deg",
                       ref, rows, sp.R, false, mpi); }
      }
      {
         Spec sp;   // T5
         sp.p = p;
         sp.R = Rot::Axis(0, 0, 1, 30.0);
         sp.S = Rotate(sp.R, Base2D());
         sp.ic = "gold";
         sp.eta_V_init = true;
         RunCase(sp, "30 deg", mpi);
      }
      {
         Spec sp;   // T6
         sp.p = p;
         sp.steady = true;
         sp.S = Base2D();
         const std::vector<Row> st = RunCase(sp, "0 deg", mpi);
         Report("T6 2D p = " + std::to_string(p) + " tensor against the steady state, a = 0.004 "
                "DOFs", ref, st, Rot(), true, mpi);
      }
      // 3D, the rotated plane
      if (g_rank == 0) { std::cout << "\n=== 3D, p = " << p << " ===\n"; }
      std::vector<Row> ref3;
      for (int k = 0; k < 2; k++)
      {
         Spec sp;
         sp.dim = 3;
         sp.p = p;
         sp.R = (k == 0) ? Rot() : Rot::Axis(1, 2, 3, 40.0);
         sp.S = Rotate(sp.R, Base3D());
         const std::string rot = (k == 0) ? "R = I" : "R = 40 deg about (1, 2, 3)";
         const std::vector<Row> rows = RunCase(sp, rot, mpi);
         if (k == 0) { ref3 = rows; }
         else { Report("T4 3D p = " + std::to_string(p) + " " + rot + " against R = I",
                       ref3, rows, sp.R, false, mpi); }
      }
   }
   MFEMFinalizePetsc();
   if (g_rank == 0)
   {
      std::cout << "\n=== Summary: " << num_tests - num_failed << "/" << num_tests
                << " passed ===\n";
   }
   return num_failed ? 1 : 0;
}

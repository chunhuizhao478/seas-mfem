// Nitsche profiling plan item 7 (clustered local time stepping), phases P1
// and P2 of lts-kernels-design-rev1 (§3.7): the operator split by element
// owner and the serial multirate leapfrog of solver/nitsche_lts_stepper.hpp,
// pure elastic (P1), then with Kelvin-Voigt damping and absorbing faces (P2).
//
// Base suite: graded tet boxes (6 tets per hex, cell widths growing
// geometrically along x, y and z), every outer face Dirichlet with zero data,
// TPV102 material, clusters from BuildLtsClustering (rules, lambda = 1, no
// cap) on the Nitsche length h_e = 6V/A.  dt_CFL = 2/omega_max (power
// iteration, omega_max^2 = max eig M^{-1}K) is the global-stepping limit;
// cluster 0 steps at dt0 = r dt_CFL and cluster c at 2^c dt0.
//   L-1      split: sum_c K^(c) equals the assembled K (max entry, 1e-13).
//   L-2 (b)  every G_e positive semidefinite: Cholesky of G_e + tau I with
//            tau = 1e-12 lambda_max(G_e) on every element; Jacobi
//            eigenvalues of the 20 worst-shaped elements (smallest
//            h_e / longest edge).
//   L-3 (a)  one cluster equals CentralDifferenceSolver to round-off
//            (1000 steps, smooth data).
//   L-4 (c)  growth scan r = 0.30, 0.31, ..., 1.00: random data, 2000 sync
//            intervals; energy at sync points, max over the last 10% of the
//            intervals over the max over the first 10%.  Gate: no growth
//            above 5% for r <= 0.60.  The design asked for r <= 1; the
//            multirate leapfrog has resonance bands from about 0.65, and
//            0.60 was accepted (Chunhui Zhao, 2026-10-19).  The scan still
//            runs to 1.00 and prints the onset.  (The design's spectral
//            radius: --phi and nitsche_lts_eig.py.)
//   L-5 (d)  energy at sync points, E = v.Mv/2 + u.Ku/2, over 10^4 intervals
//            (--long: 10^5), r = 0.5, 0.6 and 0.9, smooth and random data.
//            Gate: bounded for r <= 0.60 (L-4); r = 0.9 is printed only.
//   L-6 (e)  temporal order: dt0, dt0/2, dt0/4 against one cluster at
//            dt0/64; error ratios in [3.5, 4.5].
// P2 suite: the same boxes with the faces x = -L and x = +L absorbing (attr
// 6; the fine corner of the grading touches x = -L) and Kelvin-Voigt damping
// kv_scale = "node" in every element (eta_bar = 0.1 unless --kv).
//   L-7      split: sum_c K_eta^(c) = K_eta and sum_c C^(c) = C (products
//            with random vectors, 1e-13); (b) damping power: v.K_eta^(c) v
//            >= 0 and v.C^(c) v >= 0 for every cluster (random v), every
//            dashpot block C_e PSD (G_e of K_eta = beta_e G_e of K: L-2).
//   L-8 (a)  one cluster equals CentralDifferenceSolver on the damped
//            NitscheWaveOperator to round-off (1000 steps).
//   L-9 (a)  damped step limit: power iteration on the sync-interval map in
//            the energy norm, r = 0.30, 0.32, ..., 1.00, growth rate over the
//            last quarter; the first r with rate > 1e-6 per interval, one
//            cluster (global stepping) and clusters, for eta_bar = 0.05,
//            0.1, 0.2.  Gate: the clusters' limit is at least the global
//            limit, or 0.60 (the accepted multirate limit, L-4), whichever
//            is lower.  (Dense spectra: --phi --kv --absorb.)
//   L-10 (c) absorbing faces, no damping: energy history E(t)/E0 over 30 s
//            against global stepping at the same dt0; the absorbed fraction
//            within 0.01 of E0; second order in time.
//
// Modes:
//   --phi r file [--gts] [--order p]   the sync-interval map of the 3x3x3 box
//        (p = 1; the 2x2x2 box for p = 2)
//        (dense, state (u, dt0 v), column after column, float64, preceded by
//        its size as int64); --gts: one cluster at dt0 (the baseline).
//        tests/verification/nitsche_lts_eig.py computes its eigenvalues.
//   --mesh m.msh --fault-attr a [--order p] [--cap n] [--long]  L-1 to L-5 on
//        a mesh file: natural outer faces, the fault faces ordinary interior
//        faces (no closure), but both sides of a fault face in one cluster.
//   --psd-only with --mesh: L-2 only (large meshes).
//   --energy-only [--nint N] with --mesh: L-5 only, over N sync intervals.
//   --energy-r a,b,...: the r of L-5 (default 0.5, 0.6, 0.9).
//   --seam x [--order p]: L-4 on the box with two clusters, cluster 1 = the
//        elements with dt_e >= x min dt_e (x >= 2), to test where the
//        multirate limit comes from (the seam resonance, see L-4).
//   --kv eta_bar, --absorb: Kelvin-Voigt damping and the absorbing faces of
//        the P2 suite, for --phi; with --mesh, --kv runs L-9 at --limit-r
//        a,b,... (default 0.38, 0.45, 0.5, 0.55, 0.6, 0.65) over --nint
//        intervals, and --absorb-attr a makes boundary attribute a
//        absorbing.
//   --p1-only, --p2-only: one of the two suites.
//
// Scope (Chunhui Zhao, 2026-10-13): tetrahedral meshes only.

#include "mfem.hpp"
#include "../../domain/nitsche_wave_operator.hpp"
#include "../../solver/nitsche_cfl.hpp"
#include "../../solver/nitsche_lts_setup.hpp"
#include "../../solver/nitsche_lts_stepper.hpp"
#include "../../config/tpv102_params.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

static int num_tests = 0, num_passed = 0, num_failed = 0;

#define TEST_ASSERT(cond, msg) do { \
   num_tests++; \
   if (cond) { num_passed++; std::cout << "  PASSED: " << msg << "\n"; } \
   else { num_failed++; std::cout << "  FAILED [" << __LINE__ << "]: " << msg << "\n"; } \
} while(0)

static const real_t kL = 30.0e3;   // box [-L, L] x [-L, L] x [-L, 0]
// --seam x: two clusters, cluster 1 = the elements with dt_e >= x min dt_e
static real_t g_seam = 0.0;
static const int kAbsorbAttr = 6;   // P2 boxes: faces x = -L and x = +L

static double Now()
{
   return std::chrono::duration<double>(
      std::chrono::steady_clock::now().time_since_epoch()).count();
}

static std::string Arg(int argc, char *argv[], const char *flag, const std::string &def)
{
   for (int i = 1; i + 1 < argc; i++) { if (std::string(argv[i]) == flag) { return argv[i + 1]; } }
   return def;
}
static bool Has(int argc, char *argv[], const char *flag)
{
   for (int i = 1; i < argc; i++) { if (std::string(argv[i]) == flag) { return true; } }
   return false;
}

// Graded box: n^3 hexes of 6 tets; along each axis the cell widths grow by
// R^(1/n) per cell (last/first = R^((n-1)/n)).  Outer faces attr 5.
static Mesh *GradedBox(int n, real_t R)
{
   Mesh *mesh = new Mesh(Mesh::MakeCartesian3D(n, n, n, Element::TETRAHEDRON,
                                               2 * kL, 2 * kL, kL));
   const real_t ext[3] = {2 * kL, 2 * kL, kL};
   for (int v = 0; v < mesh->GetNV(); v++)
   {
      real_t *c = mesh->GetVertex(v);
      for (int d = 0; d < 3; d++)
      {
         const real_t s = c[d] / ext[d];
         c[d] = ext[d] * (std::pow(R, s) - 1.0) / (R - 1.0);
      }
      c[0] -= kL; c[1] -= kL; c[2] -= kL;
   }
   for (int b = 0; b < mesh->GetNBE(); b++) { mesh->SetBdrAttribute(b, 5); }
   mesh->SetAttributes();
   return mesh;
}

// P2: the box faces x = -L and x = +L absorbing (attr kAbsorbAttr).
static Mesh *MarkAbsorbing(Mesh *mesh)
{
   Array<int> v;
   for (int b = 0; b < mesh->GetNBE(); b++)
   {
      mesh->GetBdrElementVertices(b, v);
      bool lo = true, hi = true;
      for (int i = 0; i < v.Size(); i++)
      {
         const real_t x = mesh->GetVertex(v[i])[0];
         lo = lo && std::abs(x + kL) < 1e-6 * kL;
         hi = hi && std::abs(x - kL) < 1e-6 * kL;
      }
      if (lo || hi) { mesh->SetBdrAttribute(b, kAbsorbAttr); }
   }
   mesh->SetAttributes();
   return mesh;
}

static MaterialField TPV102Material()
{
   return MaterialField::MakeConstant(TPV102Params::lambda, TPV102Params::mu,
                                      TPV102Params::rho);
}

// Assembled density-weighted mass, for energies (as test_nitsche_p3.cpp).
struct MassMatrix
{
   std::unique_ptr<ParBilinearForm> form;
   void Build(NitscheWaveOperator &op)
   {
      auto &fes = op.Elasticity().GetFESpace();
      form = std::make_unique<ParBilinearForm>(&fes);
      form->AddDomainIntegrator(new VectorMassIntegrator(
         const_cast<NitscheMaterialSample &>(op.Sample()).Rho()));
      form->Assemble();
      form->Finalize();
   }
   void Mult(const Vector &x, Vector &y) const
   {
      y.SetSize(x.Size());
      form->SpMat().Mult(x, y);
   }
};

// Smooth data vanishing on the box boundary (test_nitsche_p3.cpp).
static void SmoothU(const Vector &x, Vector &u)
{
   const real_t s = std::sin(M_PI * (x(0) + kL) / (2 * kL))
                    * std::sin(M_PI * (x(1) + kL) / (2 * kL))
                    * std::sin(M_PI * (x(2) + kL) / kL);
   u.SetSize(3);
   u(0) = 1.0 * s; u(1) = 0.5 * s; u(2) = -0.3 * s;
}
static void SmoothV(const Vector &x, Vector &v)
{
   const real_t s = std::sin(2 * M_PI * (x(0) + kL) / (2 * kL))
                    * std::sin(M_PI * (x(1) + kL) / (2 * kL))
                    * std::sin(M_PI * (x(2) + kL) / kL);
   v.SetSize(3);
   v(0) = -0.2 * s; v(1) = 0.1 * s; v(2) = 0.05 * s;
}

// One problem: mesh, operator (K, M^{-1}), clusters, split, dt_CFL.
struct Problem
{
   std::unique_ptr<Mesh> smesh;
   std::unique_ptr<ParMesh> pmesh;
   std::unique_ptr<NitscheWaveOperator> op;
   NitscheLtsMesh lm;
   LtsClustering cl;
   std::unique_ptr<NitscheLtsSplit> split;   // per cluster
   std::unique_ptr<NitscheLtsSplit> whole;   // one cluster: sum_e G_e
   // P2: K_eta split the same way (on the beta-weighted moduli)
   std::unique_ptr<ParGridFunction> kv_lam, kv_mu;
   std::unique_ptr<GridFunctionCoefficient> kv_lam_c, kv_mu_c;
   std::unique_ptr<NitscheLtsSplit> kvsplit, kvwhole;
   Vector rho_e;
   MassMatrix M;
   real_t dt_cfl = 0.0;
   real_t eta = 0.0;      // Kelvin-Voigt eta_bar (0: none)
   bool absorb = false;   // absorbing faces
   int p = 1;
   bool dirichlet = true;
   std::string name;

   int N() const { return op->Height(); }

   // E = v.Mv/2 + u.Ku/2
   real_t Energy(const Vector &u, const Vector &v) const
   {
      Vector Ku, Mv;
      op->Elasticity().ApplyStiffness(u, Ku);
      M.Mult(v, Mv);
      return 0.5 * (v * Mv) + 0.5 * (u * Ku);
   }

   std::vector<const Operator *> Clusters() const
   {
      std::vector<const Operator *> k;
      for (int c = 0; c < split->NumClusters(); c++) { k.push_back(&split->K(c)); }
      return k;
   }
   std::vector<const Operator *> One() const { return {&whole->K(0)}; }

   // The stepper at dt0, clusters or one cluster, with this problem's
   // damping and absorbing faces; B_k = M + 1/2 dt_{c_e} C on the elements
   // with c_e <= k (the stepper header) are built here.
   struct Lts
   {
      std::vector<std::unique_ptr<NitscheMass>> B;
      std::unique_ptr<NitscheMultirateLeapfrog> lf;
      NitscheMultirateLeapfrog *operator->() { return lf.get(); }
   };
   Lts MakeLts(real_t dt0, bool one) const
   {
      Lts out;
      const NitscheLtsSplit &s = one ? *whole : *split;
      std::vector<const Operator *> kv, c;
      if (eta > 0.0)
      {
         const NitscheLtsSplit &ks = one ? *kvwhole : *kvsplit;
         for (int k = 0; k < ks.NumClusters(); k++) { kv.push_back(&ks.K(k)); }
      }
      std::vector<const NitscheMass *> b;
      if (absorb)
      {
         for (int k = 0; k < s.NumClusters(); k++) { c.push_back(&s.Absorbing(k)); }
         Vector hdt(pmesh->GetNE());
         for (int e = 0; e < hdt.Size(); e++)
         {
            hdt(e) = 0.5 * (one ? dt0 : std::ldexp(dt0, cl.cluster[e]));
         }
         for (int k = 0; k < s.NumClusters(); k++)
         {
            std::map<int, DenseMatrix> blk;
            for (const auto &e : op->DampingBlocks())
            {
               if (one || cl.cluster[e.first] <= k) { blk.emplace(e.first, e.second); }
            }
            out.B.push_back(std::make_unique<NitscheMass>(*pmesh, p, rho_e, &blk, 0.0, &hdt));
            b.push_back(out.B.back().get());
         }
      }
      out.lf = std::make_unique<NitscheMultirateLeapfrog>(one ? One() : Clusters(), op->Mass(),
                                                          dt0, kv, c, b);
      return out;
   }

   void SmoothData(Vector &U, Vector &V) const
   {
      ParGridFunction gu(&op->Elasticity().GetFESpace());
      ParGridFunction gv(&op->Elasticity().GetFESpace());
      VectorFunctionCoefficient cu(3, SmoothU), cv(3, SmoothV);
      gu.ProjectCoefficient(cu);
      gv.ProjectCoefficient(cv);
      U = gu;
      V = gv;
   }
};

static std::unique_ptr<Problem> MakeProblem(Mesh *mesh, int p, bool dirichlet, int fault_attr,
                                            int cap, const std::string &name,
                                            real_t eta = 0.0, int absorb_attr = -1)
{
   auto P = std::make_unique<Problem>();
   P->smesh.reset(mesh);
   P->p = p;
   P->dirichlet = dirichlet;
   P->name = name;
   P->eta = eta;
   P->absorb = (absorb_attr > 0);
   const real_t cp = std::sqrt((TPV102Params::lambda + 2.0 * TPV102Params::mu)
                               / TPV102Params::rho);
   P->lm = BuildNitscheLtsMesh(*P->smesh, p, 1.0, cp, fault_attr);
   LtsClusteringOptions opt;
   opt.wiggle_scan = false;
   opt.nc_cap = cap;
   opt.merge_loss_tol = std::numeric_limits<double>::infinity();
   P->cl = BuildLtsClustering(P->lm.dt, P->smesh->ElementToElementTable(), P->lm.fault_pairs, opt);
   if (g_seam > 0.0)
   {
      MFEM_VERIFY(g_seam >= 2.0, "--seam x needs x >= 2 (cluster 1 steps at 2 dt0)");
      const real_t dmin = *std::min_element(P->lm.dt.begin(), P->lm.dt.end());
      for (std::size_t e = 0; e < P->lm.dt.size(); e++)
      {
         P->cl.cluster[e] = (P->lm.dt[e] >= g_seam * dmin) ? 1 : 0;
      }
      P->cl.num_clusters = 2;
      P->cl.lambda = 1.0;
      P->cl.dt_base = dmin;
      // the fastest cluster-0 element on a seam, in units of min dt_e
      real_t s0 = std::numeric_limits<real_t>::infinity();
      for (const auto &f : P->lm.faces)
      {
         const int a = P->cl.cluster[f.first], b = P->cl.cluster[f.second];
         if (a != b) { s0 = std::min(s0, P->lm.dt[a == 0 ? f.first : f.second] / dmin); }
      }
      std::cout << "  --seam " << g_seam << ": fastest cluster-0 seam element dt_e = " << s0
                << " min dt_e; seam resonance predicted at r = pi/4 x " << s0 << " = "
                << M_PI / 4 * s0 << "\n";
   }
   P->pmesh = std::make_unique<ParMesh>(MPI_COMM_WORLD, *P->smesh);

   BoundaryConfig bc;
   bc.fault_attr = 9999;   // absent: the fault faces are ordinary interior faces
   Array<int> dmark;
   if (dirichlet)
   {
      bc.dirichlet_attrs = {5};
      bc.default_dirichlet_func = [](const Vector &, real_t, Vector &u)
      { u.SetSize(3); u = 0.0; };
      dmark.SetSize(P->pmesh->bdr_attributes.Max());
      dmark = 0;
      dmark[4] = 1;
   }
   if (P->absorb) { bc.absorbing_attrs = {absorb_attr}; }
   P->op = std::make_unique<NitscheWaveOperator>(*P->pmesh, p, TPV102Material(), bc);
   auto &sample = const_cast<NitscheMaterialSample &>(P->op->Sample());
   auto &fes = P->op->Elasticity().GetFESpace();
   const std::vector<int> one(P->smesh->GetNE(), 0);
   P->split = std::make_unique<NitscheLtsSplit>(fes, sample.Lambda(), sample.Mu(), dmark, 1.0);
   P->split->AssembleClusters(P->cl.cluster, P->cl.num_clusters);
   P->whole = std::make_unique<NitscheLtsSplit>(fes, sample.Lambda(), sample.Mu(), dmark, 1.0);
   P->whole->AssembleClusters(one, 1);
   RhoPerElement(sample, P->rho_e);
   if (P->absorb)
   {
      MFEM_VERIFY(P->op->HasDamping(), "MakeProblem: no absorbing face with attribute "
                  << absorb_attr);
      P->split->AssembleAbsorbingClusters(P->op->DampingBlocks(), P->cl.cluster,
                                          P->cl.num_clusters);
      P->whole->AssembleAbsorbingClusters(P->op->DampingBlocks(), one, 1);
   }
   if (eta > 0.0)
   {
      // the operator's own K_eta (whole domain, kv_scale = "node", 6V/A),
      // and the same beta_e for the split
      P->op->SetKelvinVoigt(eta, NitscheWaveOperator::KVScale::Node, -1.0);
      ParFiniteElementSpace *mfes = sample.RhoGF().ParFESpace();
      P->kv_lam = std::make_unique<ParGridFunction>(mfes);
      P->kv_mu = std::make_unique<ParGridFunction>(mfes);
      Array<int> md;
      for (int e = 0; e < P->pmesh->GetNE(); e++)
      {
         mfes->GetElementDofs(e, md);
         const real_t cs = std::sqrt(sample.MuGF()(md[0]) / sample.RhoGF()(md[0]));
         const real_t h = P->op->KVElementLength(e, NitscheWaveOperator::KVLength::Inscribed);
         const real_t beta = eta * NitscheWaveOperator::kKVGamma * (h / p) / cs;
         (*P->kv_lam)(md[0]) = beta * sample.LambdaGF()(md[0]);
         (*P->kv_mu)(md[0]) = beta * sample.MuGF()(md[0]);
      }
      P->kv_lam_c = std::make_unique<GridFunctionCoefficient>(P->kv_lam.get());
      P->kv_mu_c = std::make_unique<GridFunctionCoefficient>(P->kv_mu.get());
      P->kvsplit = std::make_unique<NitscheLtsSplit>(fes, *P->kv_lam_c, *P->kv_mu_c, dmark, 1.0);
      P->kvsplit->AssembleClusters(P->cl.cluster, P->cl.num_clusters);
      P->kvwhole = std::make_unique<NitscheLtsSplit>(fes, *P->kv_lam_c, *P->kv_mu_c, dmark, 1.0);
      P->kvwhole->AssembleClusters(one, 1);
   }
   P->M.Build(*P->op);

   // dt_CFL = 2 / omega_max, omega_max^2 by power iteration (Rayleigh
   // quotient v.Kv / v.Mv), relative change below 1e-9
   const int n = P->N();
   Vector v(n), Kv(n), Mv(n), w(n);
   v.Randomize(4242);
   v -= 0.5;
   v /= v.Norml2();
   real_t lam = 0.0, lam_old = 0.0;
   int it = 0;
   for (it = 1; it <= 50000; it++)
   {
      P->op->Elasticity().ApplyStiffness(v, Kv);
      P->M.Mult(v, Mv);
      lam = (v * Kv) / (v * Mv);
      if (it > 1 && std::abs(lam - lam_old) < 1e-9 * lam) { break; }
      lam_old = lam;
      P->op->Mass().Apply(Kv, w);
      v = w;
      v /= v.Norml2();
   }
   P->dt_cfl = 2.0 / std::sqrt(lam);

   const auto cells = P->cl.cells_per_cluster();
   std::cout << std::setprecision(6) << "  [" << name << ", p=" << p;
   if (eta > 0.0) { std::cout << ", Kelvin-Voigt eta_bar " << eta; }
   if (P->absorb) { std::cout << ", " << P->op->DampingBlocks().size() << " absorbing elements"; }
   std::cout << "] " << P->smesh->GetNE()
             << " tets, " << n << " DOFs, " << P->cl.num_clusters << " clusters (";
   for (std::size_t c = 0; c < cells.size(); c++) { std::cout << (c ? " " : "") << cells[c]; }
   std::cout << "), ideal " << P->cl.lambda * HarmonicUpdateSpeedup(P->cl)
             << "x; dt_CFL = " << P->dt_cfl << " s (" << it << " power iterations), h-based min dt_e"
             << " x C(p) = " << P->cl.dt_base * (p == 1 ? 0.113 : p == 2 ? 0.212 : 0.282) << " s\n";
   return P;
}

// ----------------------------------------------------------------- helpers
static bool CholeskyPD(DenseMatrix A)
{
   const int n = A.Height();
   for (int j = 0; j < n; j++)
   {
      real_t s = A(j, j);
      for (int k = 0; k < j; k++) { s -= A(j, k) * A(j, k); }
      if (!(s > 0.0)) { return false; }
      const real_t d = std::sqrt(s);
      A(j, j) = d;
      for (int i = j + 1; i < n; i++)
      {
         real_t t = A(i, j);
         for (int k = 0; k < j; k++) { t -= A(i, k) * A(j, k); }
         A(i, j) = t / d;
      }
   }
   return true;
}

// Eigenvalues of a symmetric matrix, cyclic Jacobi.
static void JacobiEig(DenseMatrix A, Vector &ev)
{
   const int n = A.Height();
   for (int sweep = 0; sweep < 100; sweep++)
   {
      real_t off = 0.0, tot = 0.0;
      for (int i = 0; i < n; i++)
      {
         for (int j = 0; j < n; j++)
         {
            tot += A(i, j) * A(i, j);
            if (i != j) { off += A(i, j) * A(i, j); }
         }
      }
      if (off <= 1e-32 * tot) { break; }
      for (int p = 0; p < n - 1; p++)
      {
         for (int q = p + 1; q < n; q++)
         {
            const real_t apq = A(p, q);
            if (std::abs(apq) <= 1e-300) { continue; }
            const real_t theta = (A(q, q) - A(p, p)) / (2.0 * apq);
            const real_t t = (theta >= 0 ? 1.0 : -1.0)
                             / (std::abs(theta) + std::sqrt(theta * theta + 1.0));
            const real_t c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
            for (int k = 0; k < n; k++)
            {
               const real_t akp = A(k, p), akq = A(k, q);
               A(k, p) = c * akp - s * akq;
               A(k, q) = s * akp + c * akq;
            }
            for (int k = 0; k < n; k++)
            {
               const real_t apk = A(p, k), aqk = A(q, k);
               A(p, k) = c * apk - s * aqk;
               A(q, k) = s * apk + c * aqk;
            }
         }
      }
   }
   ev.SetSize(n);
   for (int i = 0; i < n; i++) { ev(i) = A(i, i); }
}

static real_t PowerMax(const DenseMatrix &A, int iters)
{
   const int n = A.Height();
   Vector x(n), y(n);
   for (int i = 0; i < n; i++) { x(i) = 1.0 + 0.01 * (i % 7); }
   real_t lam = 0.0;
   for (int k = 0; k < iters; k++)
   {
      A.Mult(x, y);
      lam = y.Norml2() / x.Norml2();
      x = y;
      x /= x.Norml2();
   }
   return lam;
}

// ---------------------------------------------------------------- L-1 split
static void TestSplit(Problem &P)
{
   std::cout << "\n=== L-1 split, " << P.name << ", p = " << P.p << " ===\n";
   SparseMatrix Kd;
   P.op->Elasticity().ParStiffness().GetDiag(Kd);
   SparseMatrix *sum = new SparseMatrix(P.split->K(0));
   for (int c = 1; c < P.split->NumClusters(); c++)
   {
      SparseMatrix *s2 = Add(1.0, *sum, 1.0, P.split->K(c));
      delete sum;
      sum = s2;
   }
   SparseMatrix *d = Add(1.0, *sum, -1.0, Kd);
   SparseMatrix *dw = Add(1.0, P.whole->K(0), -1.0, Kd);
   const real_t kmax = Kd.MaxNorm(), e = d->MaxNorm() / kmax, ew = dw->MaxNorm() / kmax;
   long long nnz_split = 0;
   for (int c = 0; c < P.split->NumClusters(); c++) { nnz_split += P.split->K(c).NumNonZeroElems(); }
   std::cout << std::setprecision(3) << "  max |sum_c K^(c) - K| / max |K| = " << e
             << "; one cluster " << ew << "; stored entries: clusters " << nnz_split << ", K "
             << Kd.NumNonZeroElems() << " (" << 100.0 * (double(nnz_split) / Kd.NumNonZeroElems() - 1.0)
             << "% more)\n";
   TEST_ASSERT(e <= 1e-13 && ew <= 1e-13, "sum_c K^(c) == K to 1e-13 (max entry)");
   delete sum;
   delete d;
   delete dw;
}

// ------------------------------------------------------------------ L-2 PSD
static void TestPSD(Problem &P)
{
   std::cout << "\n=== L-2 (b) every G_e positive semidefinite, " << P.name << ", p = "
             << P.p << " ===\n";
   const double t0 = Now();
   Mesh &m = *P.smesh;
   const int ne = m.GetNE();
   int fails = 0, patch_max = 0;
   std::vector<int> failed;
   Array<int> vd;
   DenseMatrix G;
   for (int e = 0; e < ne; e++)
   {
      P.split->ElementPatch(e, vd, G);
      patch_max = std::max(patch_max, vd.Size());
      const real_t lmax = PowerMax(G, 60);
      DenseMatrix Gs(G);
      for (int i = 0; i < Gs.Height(); i++) { Gs(i, i) += 1e-12 * lmax; }
      if (!CholeskyPD(Gs)) { fails++; if (failed.size() < 20) { failed.push_back(e); } }
   }
   // the 20 worst-shaped elements: smallest h_e / longest edge
   std::vector<real_t> q(ne);
   for (int e = 0; e < ne; e++)
   {
      Array<int> v;
      m.GetElementVertices(e, v);
      real_t lmax = 0.0;
      for (int i = 0; i < 4; i++)
      {
         for (int j = i + 1; j < 4; j++)
         {
            const real_t *a = m.GetVertex(v[i]), *b = m.GetVertex(v[j]);
            lmax = std::max(lmax, std::sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1])
                                            + (a[2]-b[2])*(a[2]-b[2])));
         }
      }
      q[e] = P.lm.h[e] / lmax;
   }
   std::vector<int> idx(ne);
   std::iota(idx.begin(), idx.end(), 0);
   const int ns = std::min(ne, 20);
   std::partial_sort(idx.begin(), idx.begin() + ns, idx.end(),
                     [&](int a, int b) { return q[a] < q[b]; });
   std::vector<int> sample(idx.begin(), idx.begin() + ns);
   sample.insert(sample.end(), failed.begin(), failed.end());
   real_t worst = std::numeric_limits<real_t>::infinity();
   int worst_e = -1;
   Vector ev;
   for (int e : sample)
   {
      P.split->ElementPatch(e, vd, G);
      JacobiEig(G, ev);
      const real_t r = ev.Min() / ev.Max();
      if (r < worst) { worst = r; worst_e = e; }
   }
   std::cout << std::setprecision(4) << "  " << ne << " elements (patches up to " << patch_max
             << " DOFs): Cholesky of G_e + 1e-12 lambda_max I failed on "
             << fails << "; worst h_e/longest edge " << q[idx[0]] << "; over those " << ns
             << " and the failures, min lambda_min/lambda_max = " << worst << " (element "
             << worst_e << ", h_e/edge " << q[worst_e] << "); " << Now() - t0 << " s\n";
   TEST_ASSERT(fails == 0 && worst >= -1e-12,
               "every G_e PSD (lambda_min >= -1e-12 lambda_max)");
}

// --------------------------------------------------------------- L-3 Nc = 1
static void TestOneCluster(Problem &P)
{
   const bool p2 = P.eta > 0.0 || P.absorb;
   std::cout << "\n=== " << (p2 ? "L-8" : "L-3") << " (a) one cluster == CentralDifferenceSolver"
             << (p2 ? " (damped operator), " : ", ") << P.name << ", p = " << P.p << " ===\n";
   const real_t dt = 0.5 * P.dt_cfl;
   const int n = P.N(), nsteps = 1000;
   Vector U0(n), V0(n);
   P.SmoothData(U0, V0);
   P.op->SetTimeStep(dt);
   Vector Ua(U0), Va(V0);
   real_t ta = 0.0;
   RunCentralDifference(*P.op, Ua, Va, ta, dt, nsteps);
   auto lf = P.MakeLts(dt, true);
   Vector Ub(U0), Vb(V0);
   real_t tb = 0.0;
   lf->Init(Ub, Vb);
   for (int s = 0; s < nsteps; s++) { lf->Step(Ub, Vb, tb); }
   Vector dU(Ua), dV(Va);
   dU -= Ub;
   dV -= Vb;
   const real_t eu = dU.Normlinf() / Ua.Normlinf(), evv = dV.Normlinf() / Va.Normlinf();
   std::cout << std::setprecision(3) << "  dt = 0.5 dt_CFL, " << nsteps
             << " steps: max |u_lts - u_cd| / max |u_cd| = " << eu << ", velocity " << evv
             << "; t " << ta << " vs " << tb << "\n";
   TEST_ASSERT(eu <= 1e-10 && evv <= 1e-10, "one cluster == CentralDifferenceSolver to 1e-10");
}

// Energy at sync points over nint intervals from (U, V); returns
// max(E over the last 10%) / max(E over the first 10%), and the max and
// min of E/E0 (non-finite: infinity).
struct EnergyRun { real_t growth = 0, emax = 0, emin = 0; };
static EnergyRun RunEnergy(Problem &P, real_t r, int nint, Vector U, Vector V, bool one = false)
{
   EnergyRun out;
   auto lf = P.MakeLts(r * P.dt_cfl, one);
   real_t t = 0.0;
   lf->Init(U, V);
   const real_t E0 = P.Energy(U, V);
   const int w = std::max(1, nint / 10);
   real_t first = 0.0, last = 0.0;
   out.emax = 1.0;
   out.emin = 1.0;
   for (int k = 1; k <= nint; k++)
   {
      lf->Step(U, V, t);
      if (k <= w || k > nint - w || k % 64 == 0)
      {
         const real_t E = P.Energy(U, V) / E0;
         if (!std::isfinite(E) || E > 1e30)
         {
            out.growth = out.emax = std::numeric_limits<real_t>::infinity();
            return out;
         }
         if (k <= w) { first = std::max(first, E); }
         if (k > nint - w) { last = std::max(last, E); }
         out.emax = std::max(out.emax, E);
         out.emin = std::min(out.emin, E);
      }
   }
   out.growth = last / first;
   return out;
}

// omega_c^2 = max eig M^{-1} K^(c) by power iteration (Rayleigh quotient
// v.K^(c)v / v.Mv); 2/omega_c is cluster c's own leapfrog limit.
static real_t MaxEigMinvA(Problem &P, const SparseMatrix &K, int seed);
static real_t ClusterOmega2(Problem &P, int c) { return MaxEigMinvA(P, P.split->K(c), 4243 + c); }
static real_t MaxEigMinvA(Problem &P, const SparseMatrix &K, int seed)
{
   const int n = P.N();
   Vector v(n), Kv(n), Mv(n), w(n);
   v.Randomize(seed);
   v -= 0.5;
   v /= v.Norml2();
   real_t lam = 0.0, lam_old = 0.0;
   for (int it = 1; it <= 50000; it++)
   {
      K.Mult(v, Kv);
      P.M.Mult(v, Mv);
      lam = (v * Kv) / (v * Mv);
      if (it > 1 && std::abs(lam - lam_old) < 1e-9 * lam) { break; }
      lam_old = lam;
      P.op->Mass().Apply(Kv, w);
      v = w;
      v /= v.Norml2();
   }
   return lam;
}

// --------------------------------------------------------- L-4 growth scan
// The accepted multirate step limit as a fraction of dt_CFL (L-4, L-9).
static constexpr real_t kAcceptedR = 0.60;

static void TestGrowthScan(Problem &P, int nint)
{
   std::cout << "\n=== L-4 (c) growth scan, " << P.name << ", p = " << P.p << ", " << nint
             << " sync intervals per r ===\n";
   const double t0 = Now();
   {
      // each cluster's own limit against its step 2^c dt0 at r = 1
      real_t rpred = std::numeric_limits<real_t>::infinity();
      std::cout << std::setprecision(4) << "  cluster limits 2/omega_c / (2^c dt_CFL):";
      for (int c = 0; c < P.split->NumClusters(); c++)
      {
         const real_t q = 2.0 / std::sqrt(ClusterOmega2(P, c)) / (std::ldexp(1.0, c) * P.dt_cfl);
         std::cout << " c" << c << " " << q;
         rpred = std::min(rpred, q);
      }
      std::cout << "; smallest " << rpred << "\n";
   }
   const int n = P.N();
   Vector U(n), V(n);
   U.Randomize(777);
   U -= 0.5;
   V = 0.0;
   real_t r_first_bad = 0.0, worst = 0.0, worst_r = 0.0;
   std::cout << std::setprecision(4) << "  r: growth";
   for (int k = 30; k <= 100; k++)
   {
      const real_t r = 0.01 * k;
      const EnergyRun er = RunEnergy(P, r, nint, U, V);
      if (k % 5 == 0 || er.growth > 1.05) { std::cout << "  " << r << ": " << er.growth; }
      if (er.growth > worst) { worst = er.growth; worst_r = r; }
      if (er.growth > 1.05 && r_first_bad == 0.0) { r_first_bad = r; }
   }
   std::cout << "\n  largest growth " << worst << " at r = " << worst_r << "; first r with growth > 1.05: "
             << (r_first_bad > 0 ? std::to_string(r_first_bad) : std::string("none")) << "; "
             << Now() - t0 << " s\n";
   TEST_ASSERT(r_first_bad == 0.0 || r_first_bad > kAcceptedR + 1e-12,
               "no energy growth > 5% for r = 0.30 .. 0.60 (the accepted limit)");
}

// ------------------------------------------------------------- L-5 energy
static std::vector<real_t> g_energy_r = {0.5, 0.6, 0.9};   // --energy-r a,b,...

static void TestEnergy(Problem &P, int nint)
{
   std::cout << "\n=== L-5 (d) energy at sync points, " << P.name << ", p = " << P.p << ", "
             << nint << " sync intervals ===\n";
   const double t0 = Now();
   const int n = P.N();
   bool ok = true;
   for (real_t r : g_energy_r)
   {
      // above the accepted limit the runs are printed, not gated
      const bool gated = r <= kAcceptedR + 1e-12;
      const char *tag = gated ? "" : " (above the accepted limit, not gated)";
      Vector U(n), V(n);
      if (P.dirichlet)
      {
         P.SmoothData(U, V);
         const EnergyRun s = RunEnergy(P, r, nint, U, V);
         std::cout << std::setprecision(6) << "  r = " << r << ", smooth data: E/E0 in [" << s.emin
                   << ", " << s.emax << "], growth " << s.growth << tag << "\n";
         ok = ok && (!gated || (s.emax <= 1.01 && s.emin >= 0.99));
      }
      U.Randomize(778);
      U -= 0.5;
      V = 0.0;
      const EnergyRun rr = RunEnergy(P, r, nint, U, V);
      std::cout << std::setprecision(6) << "  r = " << r << ", random data: E/E0 in [" << rr.emin
                << ", " << rr.emax << "], growth " << rr.growth << tag << "\n";
      ok = ok && (!gated || rr.growth <= 1.05);
   }
   std::cout << "  " << Now() - t0 << " s\n";
   TEST_ASSERT(ok, "energy bounded for r <= 0.60 (smooth: within 1%; random: last/first window max <= 1.05)");
}

// --------------------------------------------------------------- L-6 order
static void TestOrder(Problem &P)
{
   std::cout << "\n=== " << (P.absorb ? "L-10 (c) temporal order with absorbing faces"
                                       : "L-6 (e) temporal order")
             << ", " << P.name << ", p = " << P.p << " ===\n";
   const int n = P.N(), L = P.cl.num_clusters - 1;
   const real_t dt0 = 0.5 * P.dt_cfl;
   const int nint0 = 10;
   Vector U0(n), V0(n);
   P.SmoothData(U0, V0);
   // reference: one cluster at dt0 / 64 over the same time
   Vector Ur(U0), Vr(V0);
   {
      auto lf = P.MakeLts(dt0 / 64, true);
      real_t t = 0.0;
      lf->Init(Ur, Vr);
      for (int s = 0; s < nint0 * (1 << L) * 64; s++) { lf->Step(Ur, Vr, t); }
   }
   real_t err[3];
   for (int k = 0; k < 3; k++)
   {
      auto lf = P.MakeLts(dt0 / (1 << k), false);
      Vector U(U0), V(V0);
      real_t t = 0.0;
      lf->Init(U, V);
      for (int s = 0; s < nint0 * (1 << k); s++) { lf->Step(U, V, t); }
      U -= Ur;
      err[k] = U.Norml2() / Ur.Norml2();
   }
   const real_t r1 = err[0] / err[1], r2 = err[1] / err[2];
   std::cout << std::setprecision(4) << "  T = " << nint0 * (1 << L) * dt0 << " s; errors "
             << err[0] << ", " << err[1] << ", " << err[2] << "; ratios " << r1 << ", " << r2 << "\n";
   TEST_ASSERT(r1 >= 3.5 && r1 <= 4.5 && r2 >= 3.5 && r2 <= 4.5,
               "second order: error ratios in [3.5, 4.5]");
}


// ------------------------------------------------------- L-7 P2 split, (b)
static void TestDampingSplit(Problem &P)
{
   std::cout << "\n=== L-7 split of K_eta and C, (b) damping power, " << P.name << ", p = "
             << P.p << " ===\n";
   const int n = P.N(), nc = P.split->NumClusters();
   Vector x(n), y(n), z(n), w(n);
   real_t ekv = 0.0, ec = 0.0;
   real_t pkv = std::numeric_limits<real_t>::infinity(), pc = pkv;
   for (int trial = 0; trial < 5; trial++)
   {
      x.Randomize(900 + trial);
      x -= 0.5;
      if (P.eta > 0.0)
      {
         P.op->ApplyKVStiffness(x, y);   // the operator's K_eta x
         z = 0.0;
         for (int c = 0; c < nc; c++)
         {
            P.kvsplit->K(c).Mult(x, w);
            z += w;
            pkv = std::min(pkv, (x * w) / (x * y));
         }
         z -= y;
         ekv = std::max(ekv, z.Normlinf() / y.Normlinf());
      }
      if (P.absorb)
      {
         P.op->ApplyDamping(x, y);   // the operator's C x
         z = 0.0;
         for (int c = 0; c < nc; c++)
         {
            P.split->Absorbing(c).Mult(x, w);
            z += w;
            if (x * y > 0.0) { pc = std::min(pc, (x * w) / (x * y)); }
         }
         z -= y;
         ec = std::max(ec, z.Normlinf() / y.Normlinf());
      }
   }
   int cfail = 0;
   for (const auto &b : P.op->DampingBlocks())
   {
      DenseMatrix Cs(b.second);
      const real_t lmax = PowerMax(Cs, 60);
      for (int i = 0; i < Cs.Height(); i++) { Cs(i, i) += 1e-12 * lmax; }
      if (!CholeskyPD(Cs)) { cfail++; }
   }
   std::cout << std::setprecision(3) << "  max |sum_c K_eta^(c) x - K_eta x| / max |K_eta x| = "
             << ekv << ", same for C: " << ec << " (5 random x)\n"
             << "  smallest x.K_eta^(c) x / x.K_eta x = " << pkv << ", x.C^(c) x / x.C x = "
             << pc << "; dashpot blocks not PSD: " << cfail << " of "
             << P.op->DampingBlocks().size() << "\n";
   TEST_ASSERT(ekv <= 1e-13 && ec <= 1e-13, "sum_c K_eta^(c) == K_eta and sum_c C^(c) == C");
   TEST_ASSERT(pkv >= 0.0 && pc >= 0.0 && cfail == 0,
               "(b) damping power <= 0 in every cluster (v.K_eta^(c) v, v.C^(c) v >= 0, C_e PSD)");
}

// ----------------------------------------------- L-9 (a) damped step limit
// Growth rate per sync interval of the energy norm sqrt(E), from random
// data over nint intervals: mean of log(sqrt(E_k / E_{k-1})) over the last
// quarter (no renormalisation: the stepper keeps a_c).  Infinity once E is
// not finite or above 1e200.
static real_t GrowthRate(Problem &P, real_t r, int nint, bool one)
{
   const int n = P.N();
   Vector U(n), V(n);
   U.Randomize(779);
   U -= 0.5;
   V = 0.0;
   auto lf = P.MakeLts(r * P.dt_cfl, one);
   real_t t = 0.0;
   lf->Init(U, V);
   real_t E = P.Energy(U, V), acc = 0.0;
   int cnt = 0;
   for (int k = 1; k <= nint; k++)
   {
      lf->Step(U, V, t);
      const real_t E1 = P.Energy(U, V);
      if (!std::isfinite(E1) || E1 > 1e200) { return std::numeric_limits<real_t>::infinity(); }
      if (k > nint - nint / 4) { acc += 0.5 * std::log(E1 / E); cnt++; }
      E = E1;
   }
   return acc / cnt;
}

// Kelvin-Voigt stable fraction of one operator pair (KV derivation rev 2.1
// §6): r = 1/2 max eig(M^{-1} K_eta) / omega_max, fraction sqrt(1 + r^2) - r
// of the undamped leapfrog limit 2/omega_max.
static real_t KVFraction(Problem &P, const SparseMatrix &K, const SparseMatrix &Kv, int seed)
{
   const real_t w = std::sqrt(MaxEigMinvA(P, K, seed));
   const real_t rr = 0.5 * MaxEigMinvA(P, Kv, seed + 100) / w;
   return std::sqrt(1.0 + rr * rr) - rr;
}

static std::vector<real_t> g_limit_r;   // --limit-r a,b,... (mesh mode)

// For one problem with damping: (a1) every cluster's own damped fraction
// against the whole operator's; (a2) the first r with growth rate > 1e-6
// per interval, global stepping (one cluster) and clusters.
struct DampedLimit { real_t frac_whole = 0, frac_min = 0, frac_max = 0, r_gts = 0, r_lts = 0; };
static DampedLimit MeasureDampedLimit(Problem &P, int nint, const std::vector<real_t> &rs)
{
   DampedLimit out;
   out.frac_whole = KVFraction(P, P.whole->K(0), P.kvwhole->K(0), 5000);
   out.frac_min = std::numeric_limits<real_t>::infinity();
   out.frac_max = 0.0;
   std::cout << std::setprecision(4) << "  eta_bar " << P.eta << ": damped fraction whole "
             << out.frac_whole << "; per cluster";
   for (int c = 0; c < P.split->NumClusters(); c++)
   {
      const real_t f = KVFraction(P, P.split->K(c), P.kvsplit->K(c), 5100 + 10 * c);
      out.frac_min = std::min(out.frac_min, f);
      out.frac_max = std::max(out.frac_max, f);
      std::cout << " c" << c << " " << f;
   }
   std::cout << "\n  r: rate one cluster / clusters";
   for (real_t r : rs)
   {
      const real_t g = (out.r_gts == 0.0) ? GrowthRate(P, r, nint, true) : 0.0;
      const real_t l = (out.r_lts == 0.0) ? GrowthRate(P, r, nint, false) : 0.0;
      std::cout << std::setprecision(3) << "  " << r << ": " << g << " / " << l;
      if (g > 1e-6 && out.r_gts == 0.0) { out.r_gts = r; }
      if (l > 1e-6 && out.r_lts == 0.0) { out.r_lts = r; }
      if (out.r_gts > 0.0 && out.r_lts > 0.0) { break; }
   }
   std::cout << "\n  first unstable r: one cluster "
             << (out.r_gts > 0 ? std::to_string(out.r_gts) : std::string("none")) << ", clusters "
             << (out.r_lts > 0 ? std::to_string(out.r_lts) : std::string("none")) << "\n";
   return out;
}

static void TestDampedLimit(int p, int nint)
{
   std::cout << "\n=== L-9 (a) damped step limit, Kelvin-Voigt + absorbing, p = " << p
             << ", growth rate over " << nint << " sync intervals ===\n";
   const double t0 = Now();
   std::vector<real_t> rs;
   for (int k = 15; k <= 50; k++) { rs.push_back(0.02 * k); }
   const real_t doc[3] = {0.78, 0.61, 0.41}, etas[3] = {0.05, 0.1, 0.2};
   bool ok = true;
   for (int i = 0; i < 3; i++)
   {
      Mesh *m = MarkAbsorbing(p == 1 ? GradedBox(3, 8.0) : GradedBox(2, 8.0));
      auto P = MakeProblem(m, p, true, -1, 0, p == 1 ? "3x3x3 graded box" : "2x2x2 graded box",
                           etas[i], kAbsorbAttr);
      const DampedLimit d = MeasureDampedLimit(*P, nint, rs);
      std::cout << std::setprecision(3) << "  KV rev 2.1 §6 measured fraction (p = 1 meshes): "
                << doc[i] << "\n";
      // (a1) one fraction for every cluster; (a2) clusters no worse than
      // global stepping or the accepted multirate limit 0.60, whichever
      // is lower (minus the 0.02 grid step)
      const real_t rg = d.r_gts > 0 ? d.r_gts : 1.0, rl = d.r_lts > 0 ? d.r_lts : 1.02;
      ok = ok && d.frac_max - d.frac_min <= 0.05 && rl >= std::min(rg, kAcceptedR) - 0.02 + 1e-12;
   }
   std::cout << "  " << Now() - t0 << " s\n";
   TEST_ASSERT(ok, "(a) every cluster's damped fraction within 0.05 of the others; clusters "
               "stable to min(global limit, 0.60)");
}

// --------------------------------------------- L-10 (c) absorbing faces
static void TestAbsorbing(Problem &P)
{
   std::cout << "\n=== L-10 (c) absorbing faces, energy history against global stepping, "
             << P.name << ", p = " << P.p << " ===\n";
   const int n = P.N(), L = P.cl.num_clusters - 1;
   const real_t dt0 = 0.5 * P.dt_cfl, T = 30.0;
   const int nint = int(std::ceil(T / std::ldexp(dt0, L)));
   Vector U0(n), V0(n);
   P.SmoothData(U0, V0);
   const real_t E0 = P.Energy(U0, V0);
   Vector Ul(U0), Vl(V0), Ug(U0), Vg(V0);
   auto lts = P.MakeLts(dt0, false);
   auto gts = P.MakeLts(dt0, true);
   real_t tl = 0.0, tg = 0.0, dmax = 0.0;
   lts->Init(Ul, Vl);
   gts->Init(Ug, Vg);
   std::cout << std::setprecision(4) << "  E/E0 (clusters, one cluster):";
   for (int k = 1; k <= nint; k++)
   {
      lts->Step(Ul, Vl, tl);
      for (int s = 0; s < (1 << L); s++) { gts->Step(Ug, Vg, tg); }
      const real_t el = P.Energy(Ul, Vl) / E0, eg = P.Energy(Ug, Vg) / E0;
      dmax = std::max(dmax, std::abs(el - eg));
      if (k % std::max(1, nint / 6) == 0 || k == nint)
      {
         std::cout << "  t " << tl << ": " << el << ", " << eg;
      }
   }
   const real_t al = 1.0 - P.Energy(Ul, Vl) / E0, ag = 1.0 - P.Energy(Ug, Vg) / E0;
   std::cout << "\n  absorbed by t = " << tl << " s: clusters " << al << ", one cluster " << ag
             << "; max |E_lts - E_gts| / E0 = " << dmax << "\n";
   TEST_ASSERT(std::abs(al - ag) <= 0.01, "absorbed energy within 0.01 E0 of global stepping");
}

// ------------------------------------------------------------ --phi mode
static int WritePhi(int p, real_t r, const std::string &path, bool gts, real_t eta, bool absorb)
{
   Mesh *m = (p == 1) ? GradedBox(3, 8.0) : GradedBox(2, 8.0);
   if (absorb) { MarkAbsorbing(m); }
   auto P = MakeProblem(m, p, true, -1, 0, p == 1 ? "3x3x3 graded box" : "2x2x2 graded box", eta,
                        absorb ? kAbsorbAttr : -1);
   const int n = P->N();
   const real_t dt0 = r * P->dt_cfl;
   auto lf = P->MakeLts(dt0, gts);
   std::ofstream out(path, std::ios::binary);
   const std::int64_t n2 = 2 * n;
   out.write(reinterpret_cast<const char *>(&n2), sizeof(n2));
   Vector U(n), V(n);
   std::vector<double> col(2 * n);
   for (int j = 0; j < 2 * n; j++)
   {
      U = 0.0;
      V = 0.0;
      if (j < n) { U(j) = 1.0; } else { V(j - n) = 1.0 / dt0; }
      real_t t = 0.0;
      lf->StepKicked(U, V, t);   // state (u, v+): closes the recursion with damping
      for (int i = 0; i < n; i++) { col[i] = U(i); col[n + i] = dt0 * V(i); }
      out.write(reinterpret_cast<const char *>(col.data()), sizeof(double) * col.size());
   }
   std::cout << "  wrote the " << 2 * n << " x " << 2 * n << " map at r = " << r
             << (gts ? " (one cluster)" : "") << " to " << path << "\n";
   return 0;
}

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   MFEM_VERIFY(Mpi::WorldSize() == 1, "seas_test_nitsche_lts: one rank only (P1 is serial)");
   const int p = std::atoi(Arg(argc, argv, "--order", "1").c_str());
   std::cout << "=== Nitsche LTS P1/P2 (operator split, serial multirate leapfrog) ===\n";
   const real_t eta = std::atof(Arg(argc, argv, "--kv", "0").c_str());
   const bool absorb = Has(argc, argv, "--absorb");
   if (Has(argc, argv, "--phi"))
   {
      std::string r, path;
      for (int i = 1; i + 2 < argc; i++)
      {
         if (std::string(argv[i]) == "--phi") { r = argv[i + 1]; path = argv[i + 2]; }
      }
      return WritePhi(p, std::atof(r.c_str()), path, Has(argc, argv, "--gts"), eta, absorb);
   }
   if (Has(argc, argv, "--seam"))
   {
      g_seam = std::atof(Arg(argc, argv, "--seam", "0").c_str());
      auto P = (p == 1) ? MakeProblem(GradedBox(3, 8.0), p, true, -1, 0, "3x3x3 graded box")
                        : MakeProblem(GradedBox(2, 8.0), p, true, -1, 0, "2x2x2 graded box");
      TestGrowthScan(*P, 2000);
      std::cout << "\n=== Summary: " << num_passed << "/" << num_tests << " passed, " << num_failed
                << " failed ===\n";
      return (num_failed == 0) ? 0 : 1;
   }
   const int nlong = Has(argc, argv, "--long") ? 100000 : 10000;
   if (Has(argc, argv, "--energy-r"))
   {
      g_energy_r.clear();
      std::stringstream ss(Arg(argc, argv, "--energy-r", ""));
      for (std::string tok; std::getline(ss, tok, ',');) { g_energy_r.push_back(std::atof(tok.c_str())); }
   }
   const std::string mesh_path = Arg(argc, argv, "--mesh", "");
   if (!mesh_path.empty())
   {
      const int fa = std::atoi(Arg(argc, argv, "--fault-attr", "-1").c_str());
      const int cap = std::atoi(Arg(argc, argv, "--cap", "4").c_str());
      const int aa = std::atoi(Arg(argc, argv, "--absorb-attr", "-1").c_str());
      Mesh *m = new Mesh(mesh_path.c_str(), 1, 1);
      auto P = MakeProblem(m, p, false, fa, cap, mesh_path, eta, aa);
      if (eta > 0.0)
      {
         g_limit_r = {0.38, 0.45, 0.5, 0.55, 0.6, 0.65};
         if (Has(argc, argv, "--limit-r"))
         {
            g_limit_r.clear();
            std::stringstream ss(Arg(argc, argv, "--limit-r", ""));
            for (std::string tok; std::getline(ss, tok, ',');) { g_limit_r.push_back(std::atof(tok.c_str())); }
         }
         std::cout << "\n=== L-9 (a) damped step limit, " << mesh_path << ", p = " << p << " ===\n";
         TestDampingSplit(*P);
         MeasureDampedLimit(*P, std::atoi(Arg(argc, argv, "--nint", "1000").c_str()), g_limit_r);
         std::cout << "\n=== Summary: " << num_passed << "/" << num_tests << " passed, " << num_failed
                   << " failed ===\n";
         return (num_failed == 0) ? 0 : 1;
      }
      if (Has(argc, argv, "--energy-only"))
      {
         TestEnergy(*P, std::atoi(Arg(argc, argv, "--nint", "1000").c_str()));
         std::cout << "\n=== Summary: " << num_passed << "/" << num_tests << " passed, " << num_failed
                   << " failed ===\n";
         return (num_failed == 0) ? 0 : 1;
      }
      TestPSD(*P);
      if (!Has(argc, argv, "--psd-only"))
      {
         TestSplit(*P);
         TestGrowthScan(*P, std::max(200, nlong / 50));
         TestEnergy(*P, nlong / 10);
      }
   }
   else
   {
      for (int q = 1; q <= 2 && !Has(argc, argv, "--p2-only"); q++)
      {
         auto P = (q == 1) ? MakeProblem(GradedBox(3, 8.0), q, true, -1, 0, "3x3x3 graded box")
                           : MakeProblem(GradedBox(2, 8.0), q, true, -1, 0, "2x2x2 graded box");
         TestSplit(*P);
         TestPSD(*P);
         TestOneCluster(*P);
         TestGrowthScan(*P, 2000);
         TestEnergy(*P, nlong);
         TestOrder(*P);
      }
      for (int q = 1; q <= 2 && !Has(argc, argv, "--p1-only"); q++)
      {
         const std::string nm = (q == 1) ? "3x3x3 graded box" : "2x2x2 graded box";
         auto box = [q]() { return MarkAbsorbing(q == 1 ? GradedBox(3, 8.0) : GradedBox(2, 8.0)); };
         {
            auto P = MakeProblem(box(), q, true, -1, 0, nm, eta > 0.0 ? eta : 0.1, kAbsorbAttr);
            TestDampingSplit(*P);
            TestOneCluster(*P);
         }
         TestDampedLimit(q, 3000);
         {
            auto P = MakeProblem(box(), q, true, -1, 0, nm, 0.0, kAbsorbAttr);
            TestAbsorbing(*P);
            TestOrder(*P);
         }
      }
   }
   std::cout << "\n=== Summary: " << num_passed << "/" << num_tests << " passed, " << num_failed
             << " failed ===\n";
   return (num_failed == 0) ? 0 : 1;
}

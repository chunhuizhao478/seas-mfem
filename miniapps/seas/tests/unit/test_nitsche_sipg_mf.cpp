// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// test_nitsche_sipg_mf.cpp — GPU dynamic-rupture plan D2 (gpu-dyn-nitsche
// plan rev2 §4.3): the matrix-free SIPG product of
// domain/nitsche_sipg_operator.hpp against the assembled products of
// ElasticityDomainOperator (the operator the Nitsche driver steps with), and
// a microbenchmark of its two passes.
//
// Checks, on random u and v, ||y - y_asm||_inf <= tol ||y_asm||_inf
// (plan rev2 §4.3: tol = 1e-13), with the norms over all ranks:
//   [K]    y = K u
//   [E]    y = K_eta v (u = 0)
//   [KE]   y = K u + K_eta v
// K and K_eta are built as NitscheWaveOperator builds them: the IP
// ElasticityDomainOperator on an order-0 L2 material, and the Kelvin-Voigt
// moduli beta_e lambda_e, beta_e mu_e with beta_e = eta (1/3) (h_e / p) / c_s,
// h_e = 6 V / A (SetKelvinVoigt, KVScale::Node).  --hetero h scales lambda
// and mu (and beta) per element by smooth factors in [1 - h, 1 + h] of the
// centroid, so the two sides of a face differ (no rank-dependent data).
// --dirichlet attr makes that boundary attribute a Dirichlet boundary of
// both operators (its faces enter K); the TPV decks have none.
//
// Microbenchmark (--reps R > 0): the median of R calls of each pass (volume,
// faces, both) for K and for K + K_eta, as time per call, time per tet, and
// GFLOP/s of the operator's model flop count (NitscheSipgOperator::Flops);
// on a GPU every call is synchronised.  --dense also times the dense blocks
// (plan item 8a, DenseBlockStiffness) on the same matrices (CPU only).
// --no-check skips the assembled matrices (timing on large meshes).
// --layout element|team|both (default both) picks the kernel layouts the
// checks and the benchmark run (NitscheSipgOperator::SetLayout); --batch
// 2,4,8 the team layout's elements per thread block (default 4).
//
// The rest of the D2 step on the device (with the checks): random r and
// slip, the device result against the host one at --mtol (default 1e-14,
// plan rev2 §4.3 for the mass inverse):
//   [M]    NitscheMass::Apply, a = M^{-1} r
//   [MC]   the same with damped blocks (rho_e M_e (x) I_3 + dt/2 C_e)^{-1}
//          on every fifth element, C_e = c_e A^T A with a random A
//   [Gd]   DeviceSlipLoad::AddMult against AddSlipLoadFast (rhs += G s)
//   [T]    ComputeTractionFast on the values FaultGather brought back from
//          the device against the same call on the host vector (exact)
// and the fault of plan D3 (--device-fault 1) on the device:
//   [DM]   DeviceFaultMaps::ToKernelLocal (owned three-vectors to local
//          kernel slip, with the shared-face ghost copies) against
//          ExpandOwnedToLocalFault and fault_frame::ToKernel on the host,
//          random triads and eps_f (--mtol)
//   [DT]   DeviceFaultTraction::Mult (with the neighbour exchange) against
//          ComputeTractionFast on the host (--tol)
// On the CPU and the debug device [DM] and [DT] repeat the host arithmetic
// and their errors are 0.
//
// Usage (from miniapps/seas):
//   mpirun -np N ./seas_test_nitsche_sipg_mf --mesh tpv205/mesh/tpv2053d_1000m.msh
//       --order 1|2 [--device cpu|cuda|debug] [--device-id n] [--eta 0.1]
//       [--hetero 0.2] [--dirichlet 105] [--reps 20] [--dense] [--no-check]
//       [--fault-attr 103] [--natural 101] [--absorbing 105] [--tol 1e-13]
//       [--layout element|team|both] [--batch 4[,2,8]] [--mtol 1e-14]
// Exit status 0 when every check passes.

#include "mfem.hpp"
#include "../../common/seas_device.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../domain/dense_block_stiffness.hpp"
#include "../../domain/nitsche_sipg_operator.hpp"
#include "../../domain/nitsche_mass.hpp"
#include "../../domain/nitsche_device_fault.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

namespace
{

int g_rank = 0, g_size = 1;
int num_checks = 0, num_failed = 0;

// TPV205 / TPV104 rock (rho 2670 kg/m^3, c_s 3464 m/s, c_p 6000 m/s), as
// bench_nitsche_kernels.cpp
const real_t kRho = 2670.0, kMu = 32038120320.0, kLam = 32043759360.0;

void Check(bool ok, const std::string &msg)
{
   num_checks++;
   if (!ok) { num_failed++; }
   if (g_rank == 0) { std::cout << (ok ? "  PASSED: " : "  FAILED: ") << msg << std::endl; }
}

std::string ArgStr(int argc, char *argv[], const char *flag, const std::string &def)
{
   for (int i = 1; i + 1 < argc; i++) { if (std::string(argv[i]) == flag) { return argv[i + 1]; } }
   return def;
}
bool ArgFlag(int argc, char *argv[], const char *flag)
{
   for (int i = 1; i < argc; i++) { if (std::string(argv[i]) == flag) { return true; } }
   return false;
}

std::string Sci(double x)
{
   std::ostringstream s;
   s << std::scientific << std::setprecision(2) << x;
   return s.str();
}

/// max |a - b| / max |b| over all ranks (host reads).
real_t RelErr(const Vector &a, const Vector &b)
{
   const real_t *pa = a.HostRead(), *pb = b.HostRead();
   real_t m[2] = {0.0, 0.0};
   for (int i = 0; i < a.Size(); i++)
   {
      m[0] = std::max(m[0], std::abs(pa[i] - pb[i]));
      m[1] = std::max(m[1], std::abs(pb[i]));
   }
   MPI_Allreduce(MPI_IN_PLACE, m, 2, MPITypeMap<real_t>::mpi_type, MPI_MAX, MPI_COMM_WORLD);
   return m[1] > 0.0 ? m[0] / m[1] : m[0];
}

/// Inscribed length 6 V / A of element e (SetKelvinVoigt, KVLength::Inscribed).
real_t InscribedLength(ParMesh &pm, int e)
{
   Array<int> v;
   pm.GetElementVertices(e, v);
   const real_t *x[4];
   for (int i = 0; i < 4; i++) { x[i] = pm.GetVertex(v[i]); }
   auto sub = [](const real_t *a, const real_t *b, real_t *c) { for (int d = 0; d < 3; d++) { c[d] = a[d] - b[d]; } };
   auto cross = [](const real_t *a, const real_t *b, real_t *c)
   { c[0] = a[1] * b[2] - a[2] * b[1]; c[1] = a[2] * b[0] - a[0] * b[2]; c[2] = a[0] * b[1] - a[1] * b[0]; };
   real_t e1[3], e2[3], e3[3], c[3];
   sub(x[1], x[0], e1); sub(x[2], x[0], e2); sub(x[3], x[0], e3);
   cross(e1, e2, c);
   const real_t V = std::abs(c[0] * e3[0] + c[1] * e3[1] + c[2] * e3[2]) / 6.0;
   real_t A = 0.0;
   const int F[4][3] = {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}};
   for (auto &f : F)
   {
      real_t a[3], b[3], n[3];
      sub(x[f[1]], x[f[0]], a); sub(x[f[2]], x[f[0]], b); cross(a, b, n);
      A += 0.5 * std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
   }
   return 6.0 * V / A;
}

/// Median and minimum of reps timed calls of f (each synchronised).
template <typename F>
void Time(int reps, F &&f, double &med, double &mn)
{
   f();   // warm-up
   MFEM_DEVICE_SYNC;
   std::vector<double> t(reps);
   for (int r = 0; r < reps; r++)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      const double t0 = MPI_Wtime();
      f();
      MFEM_DEVICE_SYNC;
      double dt = MPI_Wtime() - t0;
      MPI_Allreduce(MPI_IN_PLACE, &dt, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      t[r] = dt;
   }
   std::sort(t.begin(), t.end());
   med = t[reps / 2];
   mn = t[0];
}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   g_size = Mpi::WorldSize();

   const std::string device_str = ArgStr(argc, argv, "--device", "cpu");
   const int device_id = std::atoi(ArgStr(argc, argv, "--device-id", "-1").c_str());
   Device device(device_str, SelectDeviceId(device_str, device_id, MPI_COMM_WORLD));
#if defined(HYPRE_USING_GPU) && (MFEM_HYPRE_VERSION >= 23100)
   // hypre on the host, as the Nitsche driver keeps it (plan D1): the
   // assembled reference products run there
   HYPRE_SetMemoryLocation(HYPRE_MEMORY_HOST);
   HYPRE_SetExecutionPolicy(HYPRE_EXEC_HOST);
#endif
   MFEM_VERIFY(!HypreUsingGPU(), "hypre must run on the host (hypre >= 2.31)");

   const std::string mesh_path = ArgStr(argc, argv, "--mesh", "tpv205/mesh/tpv2053d_1000m.msh");
   const int p = std::atoi(ArgStr(argc, argv, "--order", "1").c_str());
   const int fault_attr = std::atoi(ArgStr(argc, argv, "--fault-attr", "103").c_str());
   const int natural_attr = std::atoi(ArgStr(argc, argv, "--natural", "101").c_str());
   const int absorbing_attr = std::atoi(ArgStr(argc, argv, "--absorbing", "105").c_str());
   const int dirichlet_attr = std::atoi(ArgStr(argc, argv, "--dirichlet", "0").c_str());
   const real_t eta = std::atof(ArgStr(argc, argv, "--eta", "0.1").c_str());
   const real_t hetero = std::atof(ArgStr(argc, argv, "--hetero", "0").c_str());
   const real_t tol = std::atof(ArgStr(argc, argv, "--tol", "1e-13").c_str());
   const real_t mtol = std::atof(ArgStr(argc, argv, "--mtol", "1e-14").c_str());
   const int reps = std::atoi(ArgStr(argc, argv, "--reps", "0").c_str());
   const bool check = !ArgFlag(argc, argv, "--no-check");
   const bool dense = ArgFlag(argc, argv, "--dense");
   const std::string layout_str = ArgStr(argc, argv, "--layout", "both");
   MFEM_VERIFY(layout_str == "element" || layout_str == "team" || layout_str == "both",
               "--layout " << layout_str << ": not element, team or both");
   // the layouts to check and time: label, layout, team batch
   struct LayoutRun { std::string label; NitscheSipgOperator::Layout layout; int batch; };
   std::vector<LayoutRun> layouts;
   if (layout_str != "team") { layouts.push_back({"element", NitscheSipgOperator::kElement, 4}); }
   if (layout_str != "element")
   {
      std::stringstream bs(ArgStr(argc, argv, "--batch", "4"));
      std::string b;
      while (std::getline(bs, b, ','))
      {
         const int bz = std::atoi(b.c_str());
         layouts.push_back({"team" + std::to_string(bz), NitscheSipgOperator::kTeam, bz});
      }
   }

   std::unique_ptr<ParMesh> pmp;
   {
      Mesh smesh(mesh_path.c_str(), 1, 1);
      pmp = std::make_unique<ParMesh>(MPI_COMM_WORLD, smesh);
   }
   ParMesh &pm = *pmp;
   long long ne_glob = pm.GetNE();
   MPI_Allreduce(MPI_IN_PLACE, &ne_glob, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
   if (g_rank == 0)
   {
      std::cout << "=== test_nitsche_sipg_mf: " << mesh_path << ", " << ne_glob << " tets, p = " << p
                << ", " << g_size << " rank(s), device " << device_str << ", eta " << eta
                << ", hetero " << hetero << ", Dirichlet attr " << dirichlet_attr << " ===" << std::endl;
      device.Print();
   }

   // ------------------------------------------------------------------ material
   // order-0 L2, as NitscheMaterialSample; the Kelvin-Voigt moduli on the
   // same space, as SetKelvinVoigt
   L2_FECollection p0(0, 3);
   ParFiniteElementSpace p0fes(&pm, &p0);
   ParGridFunction lam(&p0fes), mu(&p0fes), kv_lam(&p0fes), kv_mu(&p0fes);
   const real_t cs = std::sqrt(kMu / kRho);
   for (int e = 0; e < pm.GetNE(); e++)
   {
      Vector xc(3);
      pm.GetElementTransformation(e)->Transform(Geometries.GetCenter(Geometry::TETRAHEDRON), xc);
      // smooth in space, so the same element gets the same factors on any
      // partition
      const real_t fl = 1.0 + hetero * std::sin(7.3e-4 * xc(0) + 1.1e-3 * xc(1) - 5.9e-4 * xc(2));
      const real_t fm = 1.0 + hetero * std::cos(5.3e-4 * xc(0) - 8.7e-4 * xc(1) + 9.1e-4 * xc(2));
      const real_t fb = 1.0 + hetero * std::sin(1.7e-3 * xc(0) + 2.3e-4 * xc(2));
      Array<int> d;
      p0fes.GetElementDofs(e, d);
      lam(d[0]) = fl * kLam;
      mu(d[0]) = fm * kMu;
      const real_t beta = fb * eta * (1.0 / 3.0) * (InscribedLength(pm, e) / p) / cs;
      kv_lam(d[0]) = beta * lam(d[0]);
      kv_mu(d[0]) = beta * mu(d[0]);
   }
   for (ParGridFunction *g : {&lam, &mu, &kv_lam, &kv_mu}) { g->ExchangeFaceNbrData(); g->HostRead(); }
   GridFunctionCoefficient lam_c(&lam), mu_c(&mu), kv_lam_c(&kv_lam), kv_mu_c(&kv_mu);

   BoundaryConfig bc;
   bc.fault_attr = fault_attr;
   if (natural_attr != dirichlet_attr) { bc.natural_attrs = {natural_attr}; }
   if (absorbing_attr != dirichlet_attr) { bc.absorbing_attrs = {absorbing_attr}; }
   if (dirichlet_attr > 0) { bc.dirichlet_attrs = {dirichlet_attr}; }
   DomainConfig cfg;
   Array<int> dmark(pm.bdr_attributes.Size() ? pm.bdr_attributes.Max() : 0);
   dmark = 0;
   if (dirichlet_attr > 0) { dmark[dirichlet_attr - 1] = 1; }

   // ------------------------------------------------------------------ operators
   std::unique_ptr<ElasticityDomainOperator<ParMesh>> elast, kv;
   elast = std::make_unique<ElasticityDomainOperator<ParMesh>>(
              pm, p, lam_c, mu_c, kLam, kMu, 0.0, 0.0, 0.0, bc, DGMethod::IP,
              SolverType::OperatorOnly, cfg);
   if (check || dense)
   {
      kv = std::make_unique<ElasticityDomainOperator<ParMesh>>(
              pm, p, kv_lam_c, kv_mu_c, kLam, kMu, 0.0, 0.0, 0.0, bc, DGMethod::IP,
              SolverType::OperatorOnly, cfg);
   }
   auto *pfes = dynamic_cast<ParFiniteElementSpace *>(&elast->GetFESpace());
   MFEM_VERIFY(pfes, "no ParFiniteElementSpace");
   double t0 = MPI_Wtime();
   NitscheSipgOperator mf(*pfes, lam, mu, &kv_lam, &kv_mu, dmark, cfg.penalty_factor);
   const double t_setup = MPI_Wtime() - t0;
   long long counts[5] = {mf.NumInteriorSides(), mf.NumSharedSides(), mf.NumDirichletSides(),
                          mf.Bytes(), pm.GetNE()};
   MPI_Allreduce(MPI_IN_PLACE, counts, 5, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
   if (g_rank == 0)
   {
      std::cout << "matrix-free setup " << std::setprecision(4) << t_setup << " s (rank max not "
                << "taken); element sides: " << counts[0] << " interior, " << counts[1] << " shared, "
                << counts[2] << " Dirichlet; data " << double(counts[3]) / counts[4] << " bytes per tet"
                << std::endl;
   }

   const int n = pfes->GetVSize();
   Vector uh(n), vh(n), zero(n);
   {
      std::mt19937_64 rng(1234 + g_rank);
      std::uniform_real_distribution<real_t> U(-1.0, 1.0);
      for (int i = 0; i < n; i++) { uh(i) = 1e-5 * U(rng); vh(i) = 1e-5 * U(rng); }
   }
   zero = 0.0;
   // device copies for the matrix-free operator; the reference stays host
   Vector ud(n), vd(n), zd(n), y(n);
   for (Vector *w : {&ud, &vd, &zd, &y}) { w->UseDevice(true); }
   ud = uh;
   vd = vh;
   zd = 0.0;

   // ------------------------------------------------------------------ checks
   if (check)
   {
      if (g_rank == 0) { std::cout << "\nchecks against the assembled products (tol " << Sci(tol) << ")" << std::endl; }
      Vector Ku(n), Ev(n), KE(n);
      t0 = MPI_Wtime();
      elast->ApplyStiffness(uh, Ku);
      kv->ApplyStiffness(vh, Ev);
      const double t_asm = MPI_Wtime() - t0;
      add(Ku, Ev, KE);
      for (const LayoutRun &L : layouts)
      {
         mf.SetLayout(L.layout, L.batch);
         const std::string lab = " (" + L.label + ")";
         mf.Mult(ud, y);
         real_t err = RelErr(y, Ku);
         Check(err <= tol, "[K]  K u" + lab + ": rel err " + Sci(err));
         if (L.layout == NitscheSipgOperator::kTeam && !mf.TeamBatchFits(L.batch, true))
         {
            if (g_rank == 0)
            {
               std::cout << "  (skipped [E], [KE]" << lab << ": the block does not fit in shared "
                            "memory with Kelvin-Voigt)" << std::endl;
            }
            continue;
         }
         mf.Mult(zd, vd, y);
         err = RelErr(y, Ev);
         Check(err <= tol, "[E]  K_eta v" + lab + ": rel err " + Sci(err));
         mf.Mult(ud, vd, y);
         err = RelErr(y, KE);
         Check(err <= tol, "[KE] K u + K_eta v" + lab + ": rel err " + Sci(err));
      }
      if (g_rank == 0)
      {
         std::cout << "  (assembly of K and K_eta and the first products: " << std::setprecision(4)
                   << t_asm << " s)" << std::endl;
      }
   }

   // ------------------------------------------------------------ mass, fault
   if (check)
   {
      if (g_rank == 0)
      {
         std::cout << "\nmass inverse and fault transfers, device against host (tol " << Sci(mtol)
                   << ")" << std::endl;
      }
      std::mt19937_64 rng(4321 + g_rank);
      std::uniform_real_distribution<real_t> U(-1.0, 1.0);
      const int ne = pm.GetNE();
      Vector rho(ne);
      for (int e = 0; e < ne; e++)
      {
         Vector xc(3);
         pm.GetElementTransformation(e)->Transform(Geometries.GetCenter(Geometry::TETRAHEDRON), xc);
         rho(e) = kRho * (1.0 + hetero * std::sin(3.1e-4 * xc(0) - 6.7e-4 * xc(2)));
      }
      std::map<int, DenseMatrix> blocks;
      const int nd = ne > 0 ? pfes->GetFE(0)->GetDof() : 0;
      const real_t half_dt = 1e-3;
      for (int e = 0; e < ne; e += 5)
      {
         DenseMatrix A(3 * nd), &C = blocks[e];
         for (int i = 0; i < 3 * nd; i++)
            for (int j = 0; j < 3 * nd; j++) { A(i, j) = U(rng); }
         C.SetSize(3 * nd);
         MultAtB(A, A, C);
         C *= 100.0 * rho(e) * pm.GetElementVolume(e) / (3 * nd);
      }
      NitscheMass m0(pm, p, rho), m1(pm, p, rho, &blocks, half_dt);
      MFEM_VERIFY(m0.VSize() == n, "NitscheMass size " << m0.VSize() << " != " << n);
      Vector rh(n), ah(n), rd(n), ad(n);
      for (int i = 0; i < n; i++) { rh(i) = U(rng); }
      for (Vector *w : {&rd, &ad}) { w->UseDevice(true); }
      rd = rh;
      for (int k = 0; k < 2; k++)
      {
         const NitscheMass &m = k ? m1 : m0;
         m.Apply(rh, ah);   // host vectors: the host loop
         m.Apply(rd, ad);   // device vectors: the kernel (the host loop on the cpu device)
         const real_t err = RelErr(ad, ah);
         Check(err <= mtol, std::string(k ? "[MC] damped M^-1 r, " + std::to_string(m1.NumDampedElements())
                                        + " damped elements on rank 0"
                                        : "[M]  M^-1 r") + ": rel err " + Sci(err));
      }

      int nfd = elast->GetNumFaultDOFs();
      MPI_Allreduce(MPI_IN_PLACE, &nfd, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      if (nfd > 0)
      {
         elast->BuildFaultFaceOps();
         const int ns = 2 * elast->GetNumFaultDOFs();
         Vector sl(ns), Gh(rh), Gd(n);
         for (int i = 0; i < ns; i++) { sl(i) = 1e-3 * U(rng); }
         Gd.UseDevice(true);
         Gd = rh;
         elast->AddSlipLoadFast(sl, Gh);
         DeviceSlipLoad dsl(*elast);
         dsl.AddMult(sl, Gd);
         real_t err = RelErr(Gd, Gh);
         Check(err <= mtol, "[Gd] slip load r + G s: rel err " + Sci(err));

         std::vector<int> gv;
         AppendFaultFaceVDofs(*elast, gv);
         FaultGather fg(gv);
         Vector xf(n), T1, N1, T2, N2;
         xf = 0.0;
         fg.Gather(ud, xf);
         elast->ComputeTractionFast(uh, sl, T1, N1);
         elast->ComputeTractionFast(xf, sl, T2, N2);
         err = std::max(RelErr(T2, T1), RelErr(N2, N1));
         long long ng = fg.Size();
         MPI_Allreduce(MPI_IN_PLACE, &ng, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
         Check(err == 0.0, "[T]  fault traction from the " + std::to_string(ng)
               + " gathered values (all ranks): rel err " + Sci(err));

         // plan D3: the owned -> local kernel slip map and the traction on
         // the device
         const int ndl = elast->GetNumFaultDOFs(), nol = elast->GetNumOwnedFaultDOFs();
         DenseMatrix Bk(9, ndl);
         std::vector<real_t> eps(ndl);
         for (int k = 0; k < ndl; k++)
         {
            for (int r = 0; r < 9; r++) { Bk(r, k) = U(rng); }
            eps[k] = U(rng) < 0.0 ? -1.0 : 1.0;
         }
         Vector dow(3 * nol), dloc, skh(2 * ndl), skd;
         for (int i = 0; i < 3 * nol; i++) { dow(i) = 1e-2 * U(rng); }
         elast->ExpandOwnedToLocalFault(dow, dloc, 3);
         for (int k = 0; k < ndl; k++)
         {
            real_t t1[3], t2[3];
            for (int c = 0; c < 3; c++)
            {
               t1[c] = Bk(3 + c, k);
               t2[c] = Bk(6 + c, k);
            }
            fault_frame::ToKernel(&dloc(3 * k), eps[k], t1, t2, skh(2 * k), skh(2 * k + 1));
         }
         DeviceFaultMaps maps(*elast, Bk, eps, MPI_COMM_WORLD);
         maps.ToKernelLocal(dow, skd);
         err = RelErr(skd, skh);
         Check(err <= mtol, "[DM] owned -> local kernel slip on the device: rel err " + Sci(err));

         DeviceFaultTraction dtr(*elast, MPI_COMM_WORLD);
         Vector sld(ns), T3, N3;
         sld.UseDevice(true);
         sld = sl;
         dtr.Mult(ud, sld, T3, N3);
         err = std::max(RelErr(T3, T1), RelErr(N3, N1));
         Check(err <= tol, "[DT] fault traction on the device (" + std::to_string(dtr.Bytes())
               + " bytes of face blocks on rank 0): rel err " + Sci(err));
      }
   }

   // ------------------------------------------------------------------ bench
   if (reps > 0)
   {
      if (g_rank == 0)
      {
         std::cout << "\nmicrobenchmark: median of " << reps << " calls (min in brackets), "
                   << g_size << " rank(s), " << ne_glob << " tets" << std::endl
                   << std::left << std::setw(36) << "  layout: product, pass" << std::right << std::setw(12)
                   << "ms/call" << std::setw(12) << "(min)" << std::setw(12) << "ns/tet" << std::setw(12)
                   << "GFLOP/s" << std::setw(14) << "kflop/tet" << std::endl;
      }
      const char *pname[4] = {"", "volume", "faces", "both"};
      for (const LayoutRun &L : layouts)
         for (int kvon = 0; kvon < 2; kvon++)
            for (int pass : {NitscheSipgOperator::kVolume, NitscheSipgOperator::kFaces,
                             NitscheSipgOperator::kAll})
         {
            mf.SetLayout(L.layout, L.batch);
            if (L.layout == NitscheSipgOperator::kTeam && !mf.TeamBatchFits(L.batch, kvon == 1))
            {
               continue;
            }
            double med, mn;
            if (kvon) { Time(reps, [&]() { mf.Mult(ud, vd, y, pass); }, med, mn); }
            else { Time(reps, [&]() { mf.Mult(ud, y, pass); }, med, mn); }
            double fl = mf.Flops(kvon == 1, pass);
            MPI_Allreduce(MPI_IN_PLACE, &fl, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            if (g_rank == 0)
            {
               std::ostringstream lab;
               lab << "  " << L.label << ": " << (kvon ? "K u + K_eta v, " : "K u, ")
                   << pname[pass];
               std::cout << std::left << std::setw(36) << lab.str() << std::right << std::fixed
                         << std::setprecision(3) << std::setw(12) << med * 1e3 << std::setw(12) << mn * 1e3
                         << std::setprecision(1) << std::setw(12) << med / ne_glob * g_size * 1e9
                         << std::setw(12) << fl / med / 1e9 << std::setw(14) << fl / ne_glob / 1e3
                         << std::defaultfloat << std::endl;
            }
         }
      if (g_rank == 0)
      {
         std::cout << "  (ns/tet is per rank: time x ranks / tets; GFLOP/s is the total over the "
                      "ranks, of the model count)" << std::endl;
      }
      if (dense && device_str == "cpu")
      {
         DenseBlockStiffness db;
         int ok = db.Init(*pfes, elast->ParStiffness(), &kv->ParStiffness()) ? 1 : 0;
         MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
         if (ok)
         {
            Vector yd(n);
            db.Mult(uh, vh, yd);
            mf.SetLayout(layouts[0].layout, layouts[0].batch);
            mf.Mult(ud, vd, y);
            const real_t err = RelErr(y, yd);
            Check(err <= tol, "dense blocks and matrix-free K u + K_eta v agree: rel err " + Sci(err));
            double med, mn;
            Time(reps, [&]() { db.Mult(uh, vh, yd); }, med, mn);
            if (g_rank == 0)
            {
               std::cout << "  dense blocks (item 8a), K u + K_eta v: " << std::setprecision(4)
                         << med * 1e3 << " ms (" << mn * 1e3 << "), " << med / ne_glob * g_size * 1e9
                         << " ns per tet-rank" << std::endl;
            }
         }
         else { Check(false, "dense blocks could not be built"); }
      }
   }

   if (g_rank == 0)
   {
      std::cout << "\n" << num_checks - num_failed << " PASSED, " << num_failed << " FAILED" << std::endl;
   }
   return num_failed == 0 ? 0 : 1;
}

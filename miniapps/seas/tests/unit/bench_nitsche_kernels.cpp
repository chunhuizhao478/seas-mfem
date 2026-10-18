// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// bench_nitsche_kernels.cpp — profiling plan (nitsche-profiling, Chunhui Zhao
// 2026-10-17): timings and equivalence checks for the per-step kernels of the
// Nitsche dynamic-rupture step (NitscheWaveOperator::ImplicitSolve) and for
// the speedup prototypes the profiling report proposes.  No solver code is
// changed; every prototype below lives in this file and is checked against the
// production call it would replace.
//
//   [S] STREAM-like triad: the single-core memory bandwidth of this machine.
//   [K] the production stiffness products, ElasticityDomainOperator::
//       ApplyStiffness for K (elastic) and K_eta (Kelvin-Voigt, moduli
//       beta_e lambda, beta_e mu with the driver's beta_e = eta_bar (1/3)
//       (h_e/p)/c_s, h_e = 6V/A): time per call, nnz, bytes moved, GB/s.
//   [F] prototype: fused K u + K_eta v in one CSR sweep (the two matrices
//       share one sparsity pattern), checked against the two products.
//   [V] the volume part alone (MFEM ElasticityIntegrator): full assembly
//       (CSR), MFEM partial assembly (AssemblyLevel::PARTIAL, non-tensor
//       kernels, valid on tets) and dense element blocks; checked against
//       each other.  Runs after [T2] (partial assembly gives the mesh nodes,
//       which deletes its face-neighbour data).
//   [T] fault traction: ComputeTractionDiagnostics as the Nitsche path calls
//       it (shear + normal + stress/correction/jump diagnostics),
//       ComputeTraction (shear + normal only), and a prototype precomputed
//       per-face dense operator T_f = A_f [u_1; u_2; slip_f] built once by
//       probing ComputeTraction; checked against ComputeTraction.  Serial
//       only (the probe needs every fault face local).
//   [F2] the production fused sweep (domain/fused_stiffness_pair.hpp, plan
//       item 4) on the matrices AssembleStiffness keeps, any rank count.
//   [D] the production dense-block products (domain/dense_block_stiffness.hpp,
//       plan item 8a): K u alone and K u + K_eta v, against the hypre
//       products and the item 4 sweep, any rank count; block bytes per tet
//       and the largest asymmetry of the assembled matrices.
//   [T2] the production per-face fault operators (ComputeTractionFast,
//       AddSlipLoadFast; plan item 3) against ComputeTractionDiagnostics and
//       AssembleSlipOnlyRHS, any rank count.
//   [P] partition balance (serial run, --partitions n1,n2,...): METIS
//       partitions of the mesh, tets and fault faces per part, and the
//       predicted step-time imbalance max/mean from the measured costs per
//       tet and per fault face.
//   [H] fault output back end (--hdf-dir): per-frame cost of
//       ParaViewHDFDataCollection::Save for the fault writer's layout, deflate
//       level 1 (the decks), 0 and no filter; then the two production writers
//       on the same triangles, 1 rank, deflate level 1 and no filter: the SAFS
//       writer vtkhdf::WriteFaultPackHdf and vtkhdf::FaultHDFStream (plan item 1).
//
// Usage (from miniapps/seas):
//   ./seas_bench_nitsche_kernels --mesh tpv205/mesh/tpv2053d_1000m.msh --order 1
//       [--fault-attr 103] [--natural 101] [--absorbing 105] [--reps 20] [--partitions 576,2304]
//       [--eta 0.1] [--skip-volume]
//   ./seas_bench_nitsche_kernels --mesh big.msh --partition-only --partitions 576,2304
//       --c-tet s --c-face s     ([P] alone, costs from a smaller-mesh run)
//   ./seas_bench_nitsche_kernels --hdf-dir DIR [--hdf-cells 100000] [--hdf-fields 12]
//       [--hdf-frames 30]        ([H] alone: per-frame cost of the fault .vtkhdf writer)
// Exit status 0 when every equivalence check passes.

#include "mfem.hpp"
#include "../../domain/elasticity_operator.hpp"
#include "../../domain/dense_block_stiffness.hpp"
#include "../../domain/fused_stiffness_pair.hpp"
#include "../../dynamic/lts_partition.hpp"
#include "../../integrator/dg_elasticity_ip_combined_integrator.hpp"
#ifdef MFEM_USE_HDF5
#include "../../io/fault_vtkhdf_writer.hpp"
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <fstream>
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

void Check(bool ok, const std::string &msg)
{
   num_checks++;
   if (!ok) { num_failed++; }
   if (g_rank == 0) { std::cout << (ok ? "  PASSED: " : "  FAILED: ") << msg << "\n"; }
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

/// Resident set size of this process (GB), from /proc/self/status.
double RssGB()
{
   std::ifstream f("/proc/self/status");
   std::string k;
   while (f >> k)
   {
      if (k == "VmRSS:") { double kb; f >> kb; return kb / 1048576.0; }
   }
   return -1.0;
}

double Now()
{
   return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}

/// Median of `reps` timed calls of f (after one warm-up call), max over ranks.
template <typename F>
double TimeMedian(int reps, F &&f)
{
   f();
   std::vector<double> t(reps);
   for (int r = 0; r < reps; r++)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      const double t0 = Now();
      f();
      t[r] = Now() - t0;
   }
   std::sort(t.begin(), t.end());
   double med = t[reps / 2], glob = 0.0;
   MPI_Allreduce(&med, &glob, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   return glob;
}

/// A and B timed in alternation (A, B, A, B, ...) so that both see the same
/// machine state; min and median of each, max over ranks.  Used where the
/// difference is small against the run-to-run noise of a shared machine.
template <typename FA, typename FB>
void TimeAlternating(int reps, FA &&fa, FB &&fb, double &a_min, double &a_med, double &b_min,
                     double &b_med)
{
   fa();
   fb();
   std::vector<double> ta(reps), tb(reps);
   for (int r = 0; r < reps; r++)
   {
      MPI_Barrier(MPI_COMM_WORLD);
      double t0 = Now();
      fa();
      ta[r] = Now() - t0;
      MPI_Barrier(MPI_COMM_WORLD);
      t0 = Now();
      fb();
      tb[r] = Now() - t0;
   }
   std::sort(ta.begin(), ta.end());
   std::sort(tb.begin(), tb.end());
   double loc[4] = {ta[0], ta[reps / 2], tb[0], tb[reps / 2]}, glob[4];
   MPI_Allreduce(loc, glob, 4, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   a_min = glob[0];
   a_med = glob[1];
   b_min = glob[2];
   b_med = glob[3];
}

/// x in scientific notation (std::to_string prints 0.000000 for round-off).
std::string Sci(double x)
{
   std::ostringstream os;
   os << std::scientific << std::setprecision(2) << x;
   return os.str();
}

real_t RelDiff(const Vector &a, const Vector &b)
{
   real_t num = 0.0, den = 0.0;
   for (int i = 0; i < a.Size(); i++)
   {
      num = std::max<real_t>(num, std::abs(a(i) - b(i)));
      den = std::max<real_t>(den, std::abs(b(i)));
   }
   real_t g[2] = {num, den};
   MPI_Allreduce(MPI_IN_PLACE, g, 2, MPITypeMap<real_t>::mpi_type, MPI_MAX, MPI_COMM_WORLD);
   return g[1] > 0 ? g[0] / g[1] : g[0];
}

// TPV205/TPV104 material (rho 2670, cs 3464, cp 6000), as in the decks.
const real_t kRho = 2670.0, kMu = 32038120320.0, kLam = 32043759360.0;

/// beta_e of KVScale::Node, KVLength::Inscribed: eta (1/3) (h_e / p) / c_s,
/// h_e = 6 V / A (nitsche_wave_operator.hpp, SetKelvinVoigt).
void KVBeta(ParMesh &pm, int p, real_t eta, Vector &beta)
{
   const real_t cs = std::sqrt(kMu / kRho);
   beta.SetSize(pm.GetNE());
   for (int e = 0; e < pm.GetNE(); e++)
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
      beta(e) = eta / 3.0 * (6.0 * V / A / p) / cs;
   }
}

/// The production K assembly (elasticity_operator_assembly.inl,
/// AssembleStiffness, IP branch, no Dirichlet faces), rebuilt here so that
/// the CSR arrays are reachable.
HypreParMatrix *AssembleK(ParFiniteElementSpace &fes, Coefficient &lam, Coefficient &mu,
                          real_t penalty_factor)
{
   ParBilinearForm a(&fes);
   a.KeepNbrBlock(true);
   a.AddDomainIntegrator(new ElasticityIntegrator(lam, mu));
   a.AddInteriorFaceIntegrator(new DGElasticityIPCombinedIntegrator(lam, mu, 3, -1.0,
                                                                    penalty_factor));
   a.Assemble(0);
   a.Finalize(0);
   return a.ParallelAssemble();
}

/// Prototype [F]: y = K u + E v for two hypre matrices with one sparsity
/// pattern (checked), diag part only (serial).
struct FusedCSR
{
   const HYPRE_Int *I = nullptr, *J = nullptr;
   const HYPRE_Real *A = nullptr, *B = nullptr;
   int n = 0;
   bool Init(HypreParMatrix &K, HypreParMatrix &E)
   {
      hypre_CSRMatrix *k = hypre_ParCSRMatrixDiag((hypre_ParCSRMatrix *)K);
      hypre_CSRMatrix *e = hypre_ParCSRMatrixDiag((hypre_ParCSRMatrix *)E);
      n = hypre_CSRMatrixNumRows(k);
      if (n != hypre_CSRMatrixNumRows(e)) { return false; }
      const int nnz = hypre_CSRMatrixNumNonzeros(k);
      if (nnz != hypre_CSRMatrixNumNonzeros(e)) { return false; }
      I = hypre_CSRMatrixI(k);
      J = hypre_CSRMatrixJ(k);
      for (int i = 0; i <= n; i++) { if (I[i] != hypre_CSRMatrixI(e)[i]) { return false; } }
      for (int i = 0; i < nnz; i++) { if (J[i] != hypre_CSRMatrixJ(e)[i]) { return false; } }
      A = hypre_CSRMatrixData(k);
      B = hypre_CSRMatrixData(e);
      return true;
   }
   void Mult(const Vector &u, const Vector &v, Vector &y) const
   {
      const real_t *pu = u.GetData(), *pv = v.GetData();
      real_t *py = y.GetData();
      for (int i = 0; i < n; i++)
      {
         real_t s = 0.0;
         for (HYPRE_Int k = I[i]; k < I[i + 1]; k++) { s += A[k] * pu[J[k]] + B[k] * pv[J[k]]; }
         py[i] = s;
      }
   }
};

/// Prototype [T]: per-face dense traction operator.  Column layout of A_f:
/// [u on the "+" element (3 nd, byNODES component blocks), u on the "-"
/// element (3 nd), slip_f (2 nbf: dip, strike per kernel DOF, interleaved)].
/// Row layout: [shear dip, shear strike, normal] per kernel DOF.
struct FaceTractionOp
{
   int nd = 0, nbf = 0, nf = 0;
   std::vector<Array<int>> vd_plus, vd_minus;   // element vdofs per face
   std::vector<DenseMatrix> A;
   void Apply(const Vector &u, const Vector &slip, Vector &T, Vector &Tn) const
   {
      const int nc = 6 * nd + 2 * nbf;
      Vector x(nc), y(3 * nbf);
      T.SetSize(2 * nf * nbf);
      Tn.SetSize(nf * nbf);
      for (int f = 0; f < nf; f++)
      {
         for (int i = 0; i < 3 * nd; i++)
         {
            x(i) = u(vd_plus[f][i]);
            x(3 * nd + i) = u(vd_minus[f][i]);
         }
         for (int i = 0; i < 2 * nbf; i++) { x(6 * nd + i) = slip(2 * f * nbf + i); }
         A[f].Mult(x, y);
         for (int kk = 0; kk < nbf; kk++)
         {
            T(2 * (f * nbf + kk)) = y(3 * kk);
            T(2 * (f * nbf + kk) + 1) = y(3 * kk + 1);
            Tn(f * nbf + kk) = y(3 * kk + 2);
         }
      }
   }
};

/// [P]: METIS partitions of the serial mesh; tets and fault-face sides per
/// part; predicted step time max/mean for per-step costs c_tet, c_face (s).
void PartitionBalance(Mesh &smesh, int fault_attr, const std::vector<int> &part_list,
                      double c_tet, double c_face)
{
   std::cout << std::setprecision(4) << "\n[P] partition balance (METIS); c_tet = "
             << c_tet * 1e9 << " ns, c_face = " << c_face * 1e9 << " ns per step (ratio "
             << c_face / c_tet << ")\n";
   std::vector<int> nff(smesh.GetNE(), 0);
   long long nfault = 0;
   for (int be = 0; be < smesh.GetNBE(); be++)
   {
      if (smesh.GetBdrAttribute(be) != fault_attr) { continue; }
      int e1, e2;
      smesh.GetFaceElements(smesh.GetBdrElementFaceIndex(be), &e1, &e2);
      if (e1 >= 0) { nff[e1]++; }
      if (e2 >= 0) { nff[e2]++; }
      nfault++;
   }
   std::cout << "  mesh: " << smesh.GetNE() << " tets, " << nfault << " fault faces\n";
   // the driver's weighted cut (plan item 6): METIS weight 1 + w per
   // fault-face side, w = c_side / c_tet
   const double w = 0.5 * c_face / c_tet;
   const int ne = smesh.GetNE();
   const Table &e2e = smesh.ElementToElementTable();
   const std::vector<int> xadj(e2e.GetI(), e2e.GetI() + ne + 1);
   const std::vector<int> adj(e2e.GetJ(), e2e.GetJ() + e2e.Size_of_connections());
   const std::vector<int> cluster(ne, 0);
   std::vector<double> cost(ne);
   for (int e = 0; e < ne; e++) { cost[e] = 100.0 * (1.0 + w * nff[e]); }
   for (int np : part_list)
   {
      for (int weighted = 0; weighted < 2; weighted++)
      {
         std::vector<int> part;
         if (weighted)
         {
            LtsPartitionOptions popt;
            popt.imbalance = 1.03;   // METIS's default for k-way, as MFEM's cut
            BuildLtsAwarePartition(ne, xadj, adj, cluster, 1, np, &cost, {}, popt, part);
         }
         else
         {
            int *pa = smesh.GeneratePartitioning(np, 1);
            part.assign(pa, pa + ne);
            delete [] pa;
         }
         std::vector<long long> ne_p(np, 0), nf_p(np, 0);
         for (int e = 0; e < ne; e++) { ne_p[part[e]]++; nf_p[part[e]] += nff[e]; }
         double tmax = 0.0, tsum = 0.0, emax = 0.0, tmax_bulk = 0.0;
         long long fmax = 0, nfr = 0;
         for (int r = 0; r < np; r++)
         {
            // each side of a fault face costs half a face (the two sides of a
            // face split over two ranks are both computed)
            const double tr = ne_p[r] * c_tet + 0.5 * nf_p[r] * c_face;
            tmax = std::max(tmax, tr);
            tmax_bulk = std::max(tmax_bulk, ne_p[r] * c_tet);
            tsum += tr;
            emax = std::max<double>(emax, ne_p[r]);
            fmax = std::max(fmax, nf_p[r]);
            nfr += nf_p[r] > 0 ? 1 : 0;
         }
         std::cout << "  " << np << " parts, " << (weighted ? "weighted w = " : "MFEM default cut")
                   << (weighted ? std::to_string(w) : std::string()) << ": tets/part mean "
                   << double(ne) / np << " max " << emax << "; parts with fault faces " << nfr << " ("
                   << 100.0 * nfr / np << "%), fault-face sides per part max " << fmax
                   << "; predicted step time max/mean " << tmax / (tsum / np) << " (bulk only "
                   << tmax_bulk / (tsum / np) << ")\n";
      }
   }
}

/// [H]: per-frame cost of the fault writer's back end, ParaViewHDFDataCollection
/// (one rank-0 serial .vtkhdf, L2 p = 0 cell fields, as vtkhdf::FaultHDFState),
/// on a triangle mesh with about `cells` cells and `nfields` fields, for
/// deflate levels 1 (the decks), 0 and none.
void HdfFrameCost(const std::string &dir, int cells, int nfields, int frames)
{
   const int nx = std::max(1, static_cast<int>(std::sqrt(cells / 2.0)));
   Mesh tri = Mesh::MakeCartesian2D(nx, nx, Element::TRIANGLE, true, 30e3, 15e3);
   L2_FECollection fec(0, 2);
   FiniteElementSpace fes(&tri, &fec);
   std::vector<std::unique_ptr<GridFunction>> gfs;
   for (int k = 0; k < nfields; k++) { gfs.emplace_back(new GridFunction(&fes)); }
   std::cout << std::setprecision(4) << "\n[H] fault writer back end: " << tri.GetNE() << " triangles, "
             << nfields << " cell fields, " << frames << " frames\n";
   const struct { const char *name; double level; } modes[] = {{"deflate 1 (decks)", 1.0},
      {"deflate 0", 0.0}, {"no filter", -1.0}};
   for (const auto &m : modes)
   {
      const std::string sub = dir + "/hdf_" + std::to_string(static_cast<int>(m.level + 1));
      ParaViewHDFDataCollection dc("fault", &tri);
      dc.SetPrefixPath(sub);
      dc.SetDataFormat(VTKFormat::BINARY);
      dc.SetHDFCompression(ParaViewHDFDataCollection::HDFCompression::Deflate, m.level);
      dc.SetHighOrderOutput(false);
      for (int k = 0; k < nfields; k++) { dc.RegisterField("f" + std::to_string(k), gfs[k].get()); }
      std::vector<double> t(frames);
      for (int f = 0; f < frames; f++)
      {
         for (int k = 0; k < nfields; k++)
            for (int e = 0; e < tri.GetNE(); e++) { (*gfs[k])(e) = std::sin(0.001 * e + 0.1 * f + k); }
         dc.SetCycle(f);
         dc.SetTime(0.02 * f);
         const double t0 = Now();
         dc.Save();
         t[f] = Now() - t0;
      }
      double tail = 0.0;
      for (int f = 1; f < frames; f++) { tail += t[f]; }
      std::cout << "  " << std::left << std::setw(18) << m.name << std::right << " first frame " << t[0] * 1e3
                << " ms, later frames " << tail / (frames - 1) * 1e3 << " ms each\n";
   }
#ifdef MFEM_USE_HDF5
   // The production writers on the same triangles, as the driver calls them
   // (the last two fields constant, as param_a, param_Dc).
   vtu::LocalFaultPack pack;
   for (int e = 0; e < tri.GetNE(); e++)
   {
      Array<int> v;
      tri.GetElementVertices(e, v);
      const int b = static_cast<int>(pack.vertices.size());
      for (int i = 0; i < 3; i++)
      {
         const real_t *x = tri.GetVertex(v[i]);
         pack.vertices.push_back({x[0], 0.0, -x[1]});
      }
      pack.triangles.push_back({b, b + 1, b + 2});
   }
   std::vector<char> constant(nfields, 0);
   for (int k = 0; k < nfields; k++)
   {
      pack.field_names.push_back("f" + std::to_string(10 + k));
      pack.field_arrays.emplace_back(tri.GetNE(), 0.0);
      if (k >= nfields - 2) { constant[k] = 1; }
   }
   const auto fill = [&](int f)
   {
      for (int k = 0; k < nfields; k++)
         for (int e = 0; e < tri.GetNE(); e++)
         {
            pack.field_arrays[k][e] = constant[k] ? std::cos(0.001 * e + k)
                                      : std::sin(0.001 * e + 0.1 * f + k);
         }
   };
   for (int run = 0; run < 4; run++)
   {
      // SAFS writer / FaultHDFStream, each with deflate 1 and with no filter
      const int w = run / 2;
      const double level = (run % 2 == 0) ? 1.0 : -1.0;
      const std::string sub = dir + (w == 0 ? "/hdf_safs" : "/hdf_stream") + (level < 0 ? "_none" : "_d1");
      vtkhdf::FaultHDFState state;
      state.compression_alg = ParaViewHDFDataCollection::HDFCompression::Deflate;
      state.compression_param = level;
      vtkhdf::FaultHDFStream stream;
      stream.SetCompression(ParaViewHDFDataCollection::HDFCompression::Deflate, level);
      std::vector<double> t(frames);
      for (int f = 0; f < frames; f++)
      {
         fill(f);
         const double t0 = Now();
         if (w == 0) { vtkhdf::WriteFaultPackHdf(state, sub, f, 0.02 * f, pack, 0, 1, MPI_COMM_SELF); }
         else { stream.Write(sub, 0.02 * f, pack, constant, 0, 1, MPI_COMM_SELF); }
         t[f] = Now() - t0;
      }
      double tail = 0.0;
      for (int f = 1; f < frames; f++) { tail += t[f]; }
      std::ifstream fs(sub + "/fault.vtkhdf", std::ios::binary | std::ios::ate);
      const double mb = fs ? static_cast<double>(fs.tellg()) / 1e6 : -1.0;
      std::cout << "  " << std::left << std::setw(18) << (w == 0 ? "SAFS writer" : "FaultHDFStream")
                << std::right << " first frame " << t[0] * 1e3 << " ms, later frames "
                << tail / (frames - 1) * 1e3 << " ms each (" << (level < 0 ? "no filter" : "deflate 1")
                << ", production call; file " << mb << " MB)\n";
   }
#endif
}

}  // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   g_size = Mpi::WorldSize();

   const std::string mesh_path = ArgStr(argc, argv, "--mesh", "tpv205/mesh/tpv2053d_1000m.msh");
   const int p = std::atoi(ArgStr(argc, argv, "--order", "1").c_str());
   const int fault_attr = std::atoi(ArgStr(argc, argv, "--fault-attr", "103").c_str());
   const int natural_attr = std::atoi(ArgStr(argc, argv, "--natural", "101").c_str());
   const int absorbing_attr = std::atoi(ArgStr(argc, argv, "--absorbing", "105").c_str());
   const int reps = std::atoi(ArgStr(argc, argv, "--reps", "20").c_str());
   const real_t eta = std::atof(ArgStr(argc, argv, "--eta", "0.1").c_str());
   const std::string parts = ArgStr(argc, argv, "--partitions", "");
   const bool skip_volume = ArgFlag(argc, argv, "--skip-volume");

   Mesh smesh(mesh_path.c_str(), 1, 1);
   std::vector<int> part_list;
   {
      std::stringstream ss(parts);
      std::string tok;
      while (std::getline(ss, tok, ',')) { if (!tok.empty()) { part_list.push_back(std::atoi(tok.c_str())); } }
   }
   // [P] needs the serial mesh; done at the end with the measured costs.
   if (!ArgStr(argc, argv, "--hdf-dir", "").empty())
   {
      HdfFrameCost(ArgStr(argc, argv, "--hdf-dir", ""), std::atoi(ArgStr(argc, argv, "--hdf-cells", "100000").c_str()),
                   std::atoi(ArgStr(argc, argv, "--hdf-fields", "12").c_str()),
                   std::atoi(ArgStr(argc, argv, "--hdf-frames", "30").c_str()));
      return 0;
   }
   if (ArgFlag(argc, argv, "--partition-only"))
   {
      // [P] alone, with per-step costs measured on a smaller mesh
      PartitionBalance(smesh, fault_attr, part_list, std::atof(ArgStr(argc, argv, "--c-tet", "0").c_str()),
                       std::atof(ArgStr(argc, argv, "--c-face", "0").c_str()));
      return 0;
   }
   ParMesh pm(MPI_COMM_WORLD, smesh);
   long long ne_glob = pm.GetNE();
   MPI_Allreduce(MPI_IN_PLACE, &ne_glob, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
   if (g_rank == 0)
   {
      std::cout << "=== bench_nitsche_kernels: " << mesh_path << ", " << ne_glob << " tets, p = " << p
                << ", " << g_size << " rank(s), reps " << reps << " ===\n";
   }

   // ------------------------------------------------------------------ [S]
   double bw_triad = 0.0;
   {
      const int n = 25'000'000;   // 3 x 200 MB
      std::vector<double> a(n, 1.0), b(n, 2.0), c(n, 0.5);
      const double t = TimeMedian(5, [&]() { for (int i = 0; i < n; i++) { a[i] = b[i] + 3.0 * c[i]; } });
      bw_triad = 3.0 * 8.0 * n / t / 1e9;
      if (g_rank == 0)
      {
         std::cout << "\n[S] triad a = b + 3c, 3 x 200 MB per rank: " << std::setprecision(4) << bw_triad
                   << " GB/s per rank (" << g_size << " ranks running at once; checksum " << a[n / 2] << ")\n";
      }
   }

   // ------------------------------------------------------------------ operators
   BoundaryConfig bc;
   bc.fault_attr = fault_attr;
   bc.natural_attrs = {natural_attr};
   bc.absorbing_attrs = {absorbing_attr};
   DomainConfig cfg;
   ConstantCoefficient lam_c(kLam), mu_c(kMu);
   L2_FECollection p0(0, 3);
   ParFiniteElementSpace p0fes(&pm, &p0);
   Vector beta;
   KVBeta(pm, p, eta, beta);
   ParGridFunction kv_lam(&p0fes), kv_mu(&p0fes);
   for (int e = 0; e < pm.GetNE(); e++) { kv_lam(e) = beta(e) * kLam; kv_mu(e) = beta(e) * kMu; }
   kv_lam.ExchangeFaceNbrData();
   kv_mu.ExchangeFaceNbrData();
   GridFunctionCoefficient kv_lam_c(&kv_lam), kv_mu_c(&kv_mu);

   const double rss0 = RssGB();
   double t0 = Now();
   ElasticityDomainOperator<ParMesh> elast(pm, p, lam_c, mu_c, kLam, kMu, 0.0, 0.0, 0.0, bc,
                                           DGMethod::IP, SolverType::CG_AMG, cfg);
   ElasticityDomainOperator<ParMesh> kv(pm, p, kv_lam_c, kv_mu_c, beta.Sum() / beta.Size() * kLam,
                                        beta.Sum() / beta.Size() * kMu, 0.0, 0.0, 0.0, bc,
                                        DGMethod::IP, SolverType::CG_AMG, cfg);
   const double t_ctor = Now() - t0;
   const int n = elast.GetFESpace().GetVSize();
   Vector u(n), v(n), Ku(n), Ev(n), y(n);
   {
      std::mt19937_64 rng(1234 + g_rank);
      std::uniform_real_distribution<real_t> U(-1.0, 1.0);
      for (int i = 0; i < n; i++) { u(i) = 1e-5 * U(rng); v(i) = 1e-5 * U(rng); }   // small: no TRACTION BLOWUP prints
   }
   t0 = Now();
   elast.ApplyStiffness(u, Ku);   // lazy assembly of K (+ the unused AMG setup)
   const double t_asmK = Now() - t0;
   t0 = Now();
   kv.ApplyStiffness(v, Ev);
   const double t_asmE = Now() - t0;
   const double rss1 = RssGB();

   // ------------------------------------------------------------------ [K]
   auto *pfes = dynamic_cast<ParFiniteElementSpace *>(&elast.GetFESpace());
   std::unique_ptr<HypreParMatrix> K(AssembleK(*pfes, lam_c, mu_c, cfg.penalty_factor));
   std::unique_ptr<HypreParMatrix> E(AssembleK(*pfes, kv_lam_c, kv_mu_c, cfg.penalty_factor));
   {
      Vector z(n);
      K->Mult(u, z);
      Check(RelDiff(z, Ku) < 1e-13, "rebuilt K equals ElasticityDomainOperator::ApplyStiffness "
            "(rel diff " + std::to_string(RelDiff(z, Ku)) + ")");
      E->Mult(v, z);
      Check(RelDiff(z, Ev) < 1e-13, "rebuilt K_eta equals the Kelvin-Voigt operator's ApplyStiffness "
            "(rel diff " + std::to_string(RelDiff(z, Ev)) + ")");
   }
   const double tK = TimeMedian(reps, [&]() { elast.ApplyStiffness(u, Ku); });
   const double tE = TimeMedian(reps, [&]() { kv.ApplyStiffness(v, Ev); });
   long long nnz = K->NNZ(), rows = K->GetGlobalNumRows();
   const double bytes = 12.0 * nnz + 4.0 * rows + 3.0 * 8.0 * rows;   // values+cols, row ptr, x read (>=1x), y write
   const double bytes_mat = 12.0 * nnz + 4.0 * rows;
   if (g_rank == 0)
   {
      std::cout << std::setprecision(4)
                << "\n[K] stiffness products (production ApplyStiffness)\n"
                << "  rows " << rows << ", nnz " << nnz << " (" << double(nnz) / rows << " per row, "
                << double(nnz) / ne_glob << " per tet); one matrix " << bytes_mat / 1e9 << " GB, "
                << bytes_mat / ne_glob / 1e3 << " kB per tet\n"
                << "  operator construction (both) " << t_ctor << " s; first ApplyStiffness (assembly + "
                   "unused AMG) K " << t_asmK << " s, K_eta " << t_asmE << " s\n"
                << "  rank-0 RSS: " << rss0 << " GB after the mesh, " << rss1 << " GB after both operators "
                   "are assembled (" << (rss1 - rss0) / bytes_mat * 1e9 * g_size << " x the bytes of one CSR matrix)\n"
                << "  K u     " << tK * 1e3 << " ms  = " << tK / ne_glob * g_size * 1e9 << " ns per tet-core, "
                << bytes / tK / 1e9 / g_size << " GB/s per rank (" << 100 * bytes / tK / 1e9 / g_size / bw_triad
                << "% of triad)\n"
                << "  K_eta v " << tE * 1e3 << " ms\n";
   }

   // ------------------------------------------------------------------ [F]
   if (g_size == 1)
   {
      FusedCSR fused;
      const bool same = fused.Init(*K, *E);
      Check(same, "K and K_eta share one CSR pattern (row pointers and columns identical)");
      if (same)
      {
         Vector z1(n), z2(n);
         K->Mult(u, z1);
         E->Mult(v, z2);
         z1 += z2;
         fused.Mult(u, v, y);
         Check(RelDiff(y, z1) < 1e-13, "fused K u + K_eta v equals the two products (rel diff "
               + std::to_string(RelDiff(y, z1)) + ")");
         const double tF = TimeMedian(reps, [&]() { fused.Mult(u, v, y); });
         const double tKE = TimeMedian(reps, [&]() { K->Mult(u, z1); E->Mult(v, z2); });
         std::cout << std::setprecision(4) << "\n[F] fused sweep: K u + K_eta v\n"
                   << "  two hypre products " << tKE * 1e3 << " ms; fused " << tF * 1e3 << " ms; speedup "
                   << tKE / tF << "x\n";
      }
   }

   // ------------------------------------------------------------------ [F2]
   {
      // the production fused sweep (FusedStiffnessPair, plan item 4) on the
      // matrices AssembleStiffness keeps (finalized with skip_zeros = 1, so
      // the two patterns may differ: union pattern), any rank count
      FusedStiffnessPair pair;
      int ok = pair.Init(elast.ParStiffness(), kv.ParStiffness()) ? 1 : 0;
      MPI_Allreduce(MPI_IN_PLACE, &ok, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      Check(ok == 1, "FusedStiffnessPair::Init accepts the production K and K_eta");
      if (ok)
      {
         long long c[2] = {pair.ExtraEntries(), pair.EntriesK()};
         MPI_Allreduce(MPI_IN_PLACE, c, 2, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
         Vector z(n);
         elast.ApplyStiffness(u, Ku);
         kv.ApplyStiffness(v, Ev);
         z = Ku;
         z += Ev;
         pair.Mult(u, v, y);
         const double rd = RelDiff(y, z);
         Check(rd < 1e-13, "FusedStiffnessPair K u + K_eta v equals the two ApplyStiffness calls "
               "(rel diff " + std::to_string(rd) + ")");
         double t2min, t2med, tPmin, tPmed;
         TimeAlternating(3 * reps, [&]() { elast.ApplyStiffness(u, Ku); kv.ApplyStiffness(v, Ev); },
                         [&]() { pair.Mult(u, v, y); }, t2min, t2med, tPmin, tPmed);
         if (g_rank == 0)
         {
            std::cout << std::setprecision(4) << "\n[F2] FusedStiffnessPair on the production K, K_eta\n"
                      << "  union pattern adds " << c[0] << " of " << c[1] << " entries (" << 100.0 * c[0] / c[1]
                      << "%)\n  timed in alternation, " << 3 * reps << " pairs: two ApplyStiffness min "
                      << t2min * 1e3 << " ms, median " << t2med * 1e3 << " ms; fused min " << tPmin * 1e3
                      << " ms, median " << tPmed * 1e3 << " ms; speedup " << t2min / tPmin << "x (min), "
                      << t2med / tPmed << "x (median)\n";
         }
      }
   }

   // ------------------------------------------------------------------ [D]
   {
      // the production dense blocks (plan item 8a) on the same matrices:
      // K alone (the undamped decks) and K with K_eta in one sweep
      DenseBlockStiffness d1, d2;
      int ok[2] = {d1.Init(*pfes, elast.ParStiffness(), nullptr) ? 1 : 0,
                   d2.Init(*pfes, elast.ParStiffness(), &kv.ParStiffness()) ? 1 : 0};
      MPI_Allreduce(MPI_IN_PLACE, ok, 2, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
      Check(ok[0] == 1, "DenseBlockStiffness::Init accepts the production K (every entry in a block)");
      Check(ok[1] == 1, "DenseBlockStiffness::Init accepts the production K and K_eta");
      if (ok[0] && ok[1])
      {
         Vector z(n), z2(n);
         elast.ApplyStiffness(u, Ku);
         kv.ApplyStiffness(v, Ev);
         d1.Mult(u, z);
         const double r1 = RelDiff(z, Ku);
         Check(r1 < 1e-13, "dense blocks K u equals ApplyStiffness (rel diff " + Sci(r1) + ")");
         z = Ku;
         z += Ev;
         d2.Mult(u, v, y);
         const double r2 = RelDiff(y, z);
         Check(r2 < 1e-13, "dense blocks K u + K_eta v equals the two ApplyStiffness calls (rel diff "
               + Sci(r2) + ")");
         // with u = v the fused sweep must equal K u + K_eta u
         kv.ApplyStiffness(u, Ev);
         z = Ku;
         z += Ev;
         d2.Mult(u, u, z2);
         const double r3 = RelDiff(z2, z);
         Check(r3 < 1e-13, "dense blocks K u + K_eta u (same vector twice) equals the two products (rel diff "
               + Sci(r3) + ")");
         real_t as[4] = {d2.Asymmetry(0), d2.Asymmetry(1), d2.MaxEntry(0), d2.MaxEntry(1)};
         MPI_Allreduce(MPI_IN_PLACE, as, 4, MPITypeMap<real_t>::mpi_type, MPI_MAX, MPI_COMM_WORLD);
         long long c[5] = {d2.Bytes(), d1.Bytes(), d2.NumFaces(), d2.NumSharedBlocks(), d2.NumElements()};
         MPI_Allreduce(MPI_IN_PLACE, c, 5, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
         FusedStiffnessPair pair;
         int okp = pair.Init(elast.ParStiffness(), kv.ParStiffness()) ? 1 : 0;
         MPI_Allreduce(MPI_IN_PLACE, &okp, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
         double tHmin, tHmed, tDmin, tDmed, tPmin = 0, tPmed = 0, tQmin = 0, tQmed = 0, tKmin, tKmed, t1min, t1med;
         TimeAlternating(3 * reps, [&]() { elast.ApplyStiffness(u, Ku); kv.ApplyStiffness(v, Ev); },
                         [&]() { d2.Mult(u, v, y); }, tHmin, tHmed, tDmin, tDmed);
         if (okp)
         {
            TimeAlternating(3 * reps, [&]() { pair.Mult(u, v, z); }, [&]() { d2.Mult(u, v, y); },
                            tPmin, tPmed, tQmin, tQmed);
         }
         TimeAlternating(3 * reps, [&]() { elast.ApplyStiffness(u, Ku); }, [&]() { d1.Mult(u, z); },
                         tKmin, tKmed, t1min, t1med);
         if (g_rank == 0)
         {
            const double csr2 = 2.0 * (12.0 * nnz + 4.0 * rows);   // two CSR matrices (values + cols + row ptr)
            std::cout << std::setprecision(4) << "\n[D] dense blocks (plan item 8a), block size "
                      << d2.BlockSize() << "\n  " << c[4] << " elements, " << c[2] << " interior faces, "
                      << c[3] << " shared-face blocks; largest asymmetry |A_ij - A_ji| / max|A|: K "
                      << as[0] / as[2] << ", K_eta " << as[1] / as[3] << "\n  bytes per tet: K and K_eta "
                      << double(c[0]) / ne_glob / 1e3 << " kB, K alone " << double(c[1]) / ne_glob / 1e3
                      << " kB; the two CSR matrices " << csr2 / ne_glob / 1e3 << " kB ("
                      << csr2 / c[0] << "x more)\n"
                      << "  K u + K_eta v, timed in alternation, " << 3 * reps << " pairs: two ApplyStiffness min "
                      << tHmin * 1e3 << " ms, median " << tHmed * 1e3 << " ms; dense min " << tDmin * 1e3
                      << " ms, median " << tDmed * 1e3 << " ms; speedup " << tHmin / tDmin << "x (min), "
                      << tHmed / tDmed << "x (median)\n";
            if (okp)
            {
               std::cout << "  against the item 4 sweep: fused CSR min " << tPmin * 1e3 << " ms, median "
                         << tPmed * 1e3 << " ms; dense min " << tQmin * 1e3 << " ms, median " << tQmed * 1e3
                         << " ms; speedup " << tPmin / tQmin << "x (min), " << tPmed / tQmed << "x (median)\n";
            }
            std::cout << "  K u alone: ApplyStiffness min " << tKmin * 1e3 << " ms, median " << tKmed * 1e3
                      << " ms; dense min " << t1min * 1e3 << " ms, median " << t1med * 1e3 << " ms; speedup "
                      << tKmin / t1min << "x (min), " << tKmed / t1med << "x (median)\n";
         }
      }
   }

   // ------------------------------------------------------------------ [T]
   const int nfd = elast.GetNumFaultDOFs(), nbf = elast.GetNbfPerFace();
   const int nfi = elast.GetFaultInteriorFaces().Size(), nfs = elast.GetFaultSharedFaces().Size();
   long long nface_glob = nfi + nfs;
   MPI_Allreduce(MPI_IN_PLACE, &nface_glob, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
   ParGridFunction ugf(pfes);
   ugf = u;
   Vector slip(2 * nfd), T, Ts, Tc, jr, Tn, T2, Tn2;
   {
      std::mt19937_64 rng(77 + g_rank);
      std::uniform_real_distribution<real_t> U(-1.0, 1.0);
      for (int i = 0; i < slip.Size(); i++) { slip(i) = 1e-5 * U(rng); }
   }
   const double tTd = TimeMedian(reps, [&]() {
      elast.ComputeTractionDiagnostics(ugf, slip, T, Ts, Tc, jr, &Tn);
   });
   const double tTp = TimeMedian(reps, [&]() { elast.ComputeTraction(ugf, slip, T2, &Tn2); });
   Check(RelDiff(T2, T) < 1e-13 && RelDiff(Tn2, Tn) < 1e-13,
         "ComputeTraction (no diagnostics) gives the Nitsche path's shear and normal traction");
   const double tGd = TimeMedian(reps, [&]() { elast.AssembleSlipOnlyRHS(y, slip); });
   // the per-call fixed part of ComputeTractionImpl: a fresh ParGridFunction
   // over the whole space, a full copy of u (and, in parallel, the full
   // face-neighbour exchange, not timed here)
   const double tCopy = TimeMedian(reps, [&]() { ParGridFunction tmp(pfes); tmp = ugf; });
   double tTpre = -1.0;
   if (g_size == 1 && nfs == 0)
   {
      // Build the per-face operator by probing ComputeTraction.  Sides: the
      // two elements of a planar-fault face lie on opposite sides of the
      // fault plane (sign of (centroid - face centroid) . face normal), so
      // one probe per (side, component, element DOF) never mixes faces.
      FaceTractionOp op;
      op.nd = pfes->GetFE(0)->GetDof();
      op.nbf = nbf;
      op.nf = nfi;
      op.vd_plus.resize(nfi);
      op.vd_minus.resize(nfi);
      op.A.resize(nfi);
      const int nd = op.nd, nc = 6 * nd + 2 * nbf;
      std::vector<int> elem_side(pm.GetNE(), 0);
      std::vector<std::array<int, 2>> fe(nfi);
      bool sides_ok = true;
      for (int fi = 0; fi < nfi; fi++)
      {
         const int face = elast.GetFaultInteriorFaces()[fi];
         FaceElementTransformations *Tr = pm.GetInteriorFaceTransformations(face);
         Tr->SetAllIntPoints(&Geometries.GetCenter(Tr->GetGeometryType()));
         Vector nor(3), xf(3), c1(3), c2(3);
         CalcOrtho(Tr->Jacobian(), nor);
         Tr->Transform(Geometries.GetCenter(Tr->GetGeometryType()), xf);
         Tr->Elem1->Transform(Geometries.GetCenter(Tr->Elem1->GetGeometryType()), c1);
         Tr->Elem2->Transform(Geometries.GetCenter(Tr->Elem2->GetGeometryType()), c2);
         c1 -= xf;
         c2 -= xf;
         const int s1 = (c1 * nor > 0) ? 1 : -1, s2 = (c2 * nor > 0) ? 1 : -1;
         // fix the side by the reference normal sign (a global plane)
         const int e_plus = (s1 > 0) ? Tr->Elem1No : Tr->Elem2No;
         const int e_minus = (s1 > 0) ? Tr->Elem2No : Tr->Elem1No;
         if (s1 == s2) { sides_ok = false; }
         // global plane check: compare against the first face's normal
         static Vector n0;
         if (fi == 0) { n0 = nor; }
         const bool flip = (nor * n0) < 0;
         fe[fi] = flip ? std::array<int, 2>{e_minus, e_plus} : std::array<int, 2>{e_plus, e_minus};
         for (int s = 0; s < 2; s++)
         {
            int &es = elem_side[fe[fi][s]];
            const int want = s == 0 ? 1 : -1;
            if (es != 0 && es != want) { sides_ok = false; }
            es = want;
         }
         pfes->GetElementVDofs(fe[fi][0], op.vd_plus[fi]);
         pfes->GetElementVDofs(fe[fi][1], op.vd_minus[fi]);
         op.A[fi].SetSize(3 * nbf, nc);
      }
      Check(sides_ok, "every fault-adjacent element lies on one side of the fault plane");
      t0 = Now();
      const real_t ps = 1e-6;   // probe amplitude (keeps the traction blow-up check quiet)
      ParGridFunction pr(pfes);
      Vector zero_slip(2 * nfd), Tp, Tnp;
      zero_slip = 0.0;
      for (int s = 0; s < 2; s++)
         for (int i = 0; i < 3 * nd; i++)
         {
            pr = 0.0;
            for (int fi = 0; fi < nfi; fi++)
            {
               pr((s == 0 ? op.vd_plus : op.vd_minus)[fi][i]) = ps;
            }
            elast.ComputeTraction(pr, zero_slip, Tp, &Tnp);
            for (int fi = 0; fi < nfi; fi++)
               for (int kk = 0; kk < nbf; kk++)
               {
                  op.A[fi](3 * kk, s * 3 * nd + i) = Tp(2 * (fi * nbf + kk)) / ps;
                  op.A[fi](3 * kk + 1, s * 3 * nd + i) = Tp(2 * (fi * nbf + kk) + 1) / ps;
                  op.A[fi](3 * kk + 2, s * 3 * nd + i) = Tnp(fi * nbf + kk) / ps;
               }
         }
      pr = 0.0;
      for (int j = 0; j < 2 * nbf; j++)
      {
         Vector sl(2 * nfd);
         sl = 0.0;
         for (int fi = 0; fi < nfi; fi++) { sl(2 * fi * nbf + j) = ps; }
         elast.ComputeTraction(pr, sl, Tp, &Tnp);
         for (int fi = 0; fi < nfi; fi++)
            for (int kk = 0; kk < nbf; kk++)
            {
               op.A[fi](3 * kk, 6 * nd + j) = Tp(2 * (fi * nbf + kk)) / ps;
               op.A[fi](3 * kk + 1, 6 * nd + j) = Tp(2 * (fi * nbf + kk) + 1) / ps;
               op.A[fi](3 * kk + 2, 6 * nd + j) = Tnp(fi * nbf + kk) / ps;
            }
      }
      const double t_probe = Now() - t0;
      Vector T3, Tn3;
      op.Apply(u, slip, T3, Tn3);
      const real_t e1 = RelDiff(T3, T2), e2 = RelDiff(Tn3, Tn2);
      std::ostringstream m;
      m << "precomputed per-face traction operator equals ComputeTraction (shear rel diff " << e1
        << ", normal " << e2 << ")";
      Check(e1 < 1e-11 && e2 < 1e-11, m.str());
      tTpre = TimeMedian(reps, [&]() { op.Apply(u, slip, T3, Tn3); });
      std::cout << std::setprecision(4) << "  per-face operator: " << nfi << " dense " << 3 * nbf << " x " << nc
                << " blocks (" << 8.0 * nfi * 3 * nbf * nc / 1e6 << " MB), built by " << 2 * 3 * nd + 2 * nbf
                << " probes in " << t_probe << " s\n";
   }
   if (g_rank == 0)
   {
      std::cout << std::setprecision(4) << "\n[T] fault traction (" << nface_glob << " fault faces, " << nbf
                << " DOFs per face)\n"
                << "  ComputeTractionDiagnostics (as the Nitsche step calls it) " << tTd * 1e3 << " ms = "
                << tTd / nface_glob * g_size * 1e6 << " us per face\n"
                << "  ComputeTraction (shear + normal only)                    " << tTp * 1e3 << " ms ("
                << tTd / tTp << "x)\n"
                << "  of which the whole-space ParGridFunction alloc + copy    " << tCopy * 1e3 << " ms; per-face "
                   "marginal cost of the Nitsche call " << (tTd - tCopy) / nface_glob * g_size * 1e6 << " us\n";
      if (tTpre > 0)
      {
         std::cout << "  precomputed per-face operator                           " << tTpre * 1e3 << " ms ("
                   << tTd / tTpre << "x vs the Nitsche call)\n";
      }
      std::cout << "  AssembleSlipOnlyRHS (G d)                               " << tGd * 1e3 << " ms\n";
   }

   // ------------------------------------------------------------------ [T2]
   {
      // the production per-face operators (plan item 3) against the calls
      // they replace in the Nitsche step, any rank count (shared fault faces
      // exercise the face-neighbour exchange)
      MPI_Barrier(MPI_COMM_WORLD);
      double t0 = Now();
      elast.BuildFaultFaceOps();
      double t_build = Now() - t0;
      MPI_Allreduce(MPI_IN_PLACE, &t_build, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      Vector Tf, Tnf, g1(n), g2(n);
      elast.ComputeTractionDiagnostics(ugf, slip, T, Ts, Tc, jr, &Tn);
      elast.ComputeTractionFast(u, slip, Tf, Tnf);
      const real_t e1 = RelDiff(Tf, T), e2 = RelDiff(Tnf, Tn);
      std::ostringstream m1;
      m1 << "ComputeTractionFast equals ComputeTractionDiagnostics (shear rel diff " << e1
         << ", normal " << e2 << ")";
      Check(e1 < 1e-11 && e2 < 1e-11, m1.str());
      elast.AssembleSlipOnlyRHS(g1, slip);
      g2 = 0.0;
      elast.AddSlipLoadFast(slip, g2);
      const real_t e3 = RelDiff(g2, g1);
      std::ostringstream m2;
      m2 << "AddSlipLoadFast equals AssembleSlipOnlyRHS (rel diff " << e3 << ")";
      Check(e3 < 1e-11, m2.str());
      const double tTf = TimeMedian(reps, [&]() { elast.ComputeTractionFast(u, slip, Tf, Tnf); });
      const double tGf = TimeMedian(reps, [&]() { elast.AddSlipLoadFast(slip, g2); });
      if (g_rank == 0)
      {
         std::cout << std::setprecision(4) << "\n[T2] production per-face fault operators (plan item 3)\n"
                   << "  build " << t_build << " s; ComputeTractionFast " << tTf * 1e3 << " ms ("
                   << tTd / tTf << "x vs ComputeTractionDiagnostics); AddSlipLoadFast " << tGf * 1e3
                   << " ms (" << tGd / tGf << "x vs AssembleSlipOnlyRHS)\n";
      }
   }

   // ------------------------------------------------------------------ [V]
   // after [T] and [T2]: MFEM partial assembly calls Mesh::EnsureNodes, and
   // ParMesh::SetCurvature then deletes the face-neighbour data that the
   // traction calls need on ranks with shared fault faces
   if (!skip_volume)
   {
      ParBilinearForm af(pfes), ap(pfes);
      af.AddDomainIntegrator(new ElasticityIntegrator(lam_c, mu_c));
      af.Assemble(0);
      af.Finalize(0);
      std::unique_ptr<HypreParMatrix> KV(af.ParallelAssemble());
      ap.SetAssemblyLevel(AssemblyLevel::PARTIAL);
      ap.AddDomainIntegrator(new ElasticityIntegrator(lam_c, mu_c));
      ap.Assemble();
      // dense element blocks, element vdofs (byNODES within the element)
      const int ne = pm.GetNE();
      std::vector<DenseMatrix> blk(ne);
      std::vector<Array<int>> evd(ne);
      {
         ElasticityIntegrator ei(lam_c, mu_c);
         for (int e = 0; e < ne; e++)
         {
            pfes->GetElementVDofs(e, evd[e]);
            ei.AssembleElementMatrix(*pfes->GetFE(e), *pm.GetElementTransformation(e), blk[e]);
         }
      }
      auto block_mult = [&](const Vector &x, Vector &z) {
         z = 0.0;
         Vector xe, ze;
         for (int e = 0; e < ne; e++)
         {
            x.GetSubVector(evd[e], xe);
            ze.SetSize(xe.Size());
            blk[e].Mult(xe, ze);
            z.AddElementVector(evd[e], ze);
         }
      };
      Vector z1(n), z2(n), z3(n);
      KV->Mult(u, z1);
      ap.Mult(u, z2);
      block_mult(u, z3);
      Check(RelDiff(z2, z1) < 1e-12, "volume part: MFEM partial assembly equals full assembly (rel diff "
            + std::to_string(RelDiff(z2, z1)) + ")");
      Check(RelDiff(z3, z1) < 1e-12, "volume part: dense element blocks equal full assembly (rel diff "
            + std::to_string(RelDiff(z3, z1)) + ")");
      const double tVf = TimeMedian(reps, [&]() { KV->Mult(u, z1); });
      const double tVp = TimeMedian(reps, [&]() { ap.Mult(u, z2); });
      const double tVb = TimeMedian(reps, [&]() { block_mult(u, z3); });
      if (g_rank == 0)
      {
         std::cout << std::setprecision(4) << "\n[V] volume term only (ElasticityIntegrator)\n"
                   << "  nnz " << KV->NNZ() << " (" << 100.0 * KV->NNZ() / nnz << "% of K)\n"
                   << "  full assembly (CSR) " << tVf * 1e3 << " ms; partial assembly " << tVp * 1e3
                   << " ms (" << tVf / tVp << "x); dense element blocks " << tVb * 1e3 << " ms ("
                   << tVf / tVb << "x)\n"
                   << "  so the face part of K u is about " << (tK - tVf) * 1e3 << " ms of the "
                   << tK * 1e3 << " ms\n";
      }
   }

   // ------------------------------------------------------------------ [P]
   if (!part_list.empty() && g_size == 1)
   {
      // per-step costs measured above: per tet (K u + K_eta v), per fault
      // face (two traction calls + two slip loads)
      PartitionBalance(smesh, fault_attr, part_list, (tK + tE) / ne_glob,
                       (2.0 * (tTd - tCopy) + 2.0 * tGd) / std::max<long long>(nface_glob, 1));
   }

   if (g_rank == 0)
   {
      std::cout << "\n=== " << num_checks - num_failed << " of " << num_checks << " checks passed ===\n";
   }
   return num_failed == 0 ? 0 : 1;
}

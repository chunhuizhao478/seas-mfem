// bench_ader_kernel_variants.cpp
//
// Microbenchmark for the ADER kernel-efficiency plan's Phase-2 bet:
// can MFEM-native batching (contiguous operator packs + tiled/fused CK
// recursion) approach SeisSol per-update cost, versus the current
// per-element-DenseMatrix / whole-vector-pass structure?
//
// Simulates the ADER-O4 p3 predictor over NE synthetic elements:
//   ndof = 20 (p3 tet), NUM_STATE = 9, 3 directions, 3 CK levels,
//   sparse 9x9 flux-Jacobian pass (~10 nnz per direction), and
//   nodal (Taylor) + time-integral accumulations.
//
// Variants:
//   A) CURRENT-STRUCTURE REPLICA: std::vector<std::array<DenseMatrix,3>>
//      operators, component-major state Q[c*Nt + e*20 + i], per CK level
//      3 whole-vector matvec passes with the TRANSPOSED inner loop
//      s += D(i,j)*Qc[j] (wave_operator.inl:1244-49 Cached branch),
//      then a strided whole-vector Jacobian AXPY pass, then whole-vector
//      accumulation AXPYs.
//   B) TILED-FUSED: operators repacked contiguous (DenseTensor, one per
//      direction, row-major-packed slice per element), tiles of E=64
//      elements, per element ALL 3 levels back-to-back, unit-stride inner
//      loops, Jacobian + both factorial accumulations fused per element,
//      state gathered to a (20x9) panel per element at tile start,
//      scattered once at end.
//   C) BatchedLinAlg NATIVE: same (column-major) DenseTensor operators,
//      per level per direction one BatchedLinAlg::Get(NATIVE).Mult over
//      the whole (20,9,NE) batch; Jacobian pass separate.
//   D) per-element LAPACK dgemm via mfem::Mult(DenseMatrix...) (20x20)x(20x9),
//      rest fused like B.
//
// All variants compute the same math from the same deterministic
// pseudo-random data; checksums agree to round-off (A and B bit-identical
// by construction: same summation order).
//
// Single-threaded (MFEM_USE_OPENMP=NO); run with OPENBLAS_NUM_THREADS=1
// VECLIB_MAXIMUM_THREADS=1 so variant D's dgemm is 1-core too.

#include "mfem.hpp"
#include "linalg/batched/batched.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <vector>

using namespace mfem;

// ----------------------------------------------------------------------------
// Problem shape
// ----------------------------------------------------------------------------
static const int NE    = 20000; // elements
static const int NDOF  = 20;    // p3 tet L2 dofs
static const int NS    = 9;     // NUM_STATE (velocity-stress)
static const int NDIR  = 3;     // spatial directions
static const int NLEV  = 3;     // CK levels for ADER-O4
static const int NNZJ  = 10;    // nnz per 9x9 flux Jacobian
static const int TILE  = 64;    // tile size for variant B/D
static const int NT    = NE * NDOF;      // dofs per component
static const int PANEL = NDOF * NS;      // 180 doubles per element panel

static const double DT = 1.0e-3;

// FLOP count per element per predictor (stated in report):
//   matvecs : 3 dirs x 3 levels x 9 comps x (20x20 MAC x 2) = 64,800
//   Jacobian: 3 levels x 3 dirs x 10 nnz x 20 dofs x 2      =  3,600
//   accums  : 2 accumulators x 4 levels x 180 x 2           =  2,880
static const double FLOPS_PER_ELEM = 64800.0 + 3600.0 + 2880.0; // 71,280

// ----------------------------------------------------------------------------
// Deterministic pseudo-random data (splitmix64 hash -> [-0.5, 0.5])
// ----------------------------------------------------------------------------
static inline double hval(uint64_t x)
{
   x += 0x9e3779b97f4a7c15ULL;
   x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
   x ^= x >> 27; x *= 0x94d049bb133111ebULL;
   x ^= x >> 31;
   return (double)(x >> 11) * (1.0 / 9007199254740992.0) - 0.5;
}

// element derivative operator D_d^e(i,j), scaled to keep CK recursion bounded
static inline double opval(int e, int d, int i, int j)
{
   uint64_t key = ((((uint64_t)e * NDIR + d) * NDOF + i) * NDOF + j);
   return 0.05 * hval(key ^ 0x1111111111111111ULL);
}

// initial state Q(e, c, i)
static inline double qval(int e, int c, int i)
{
   uint64_t key = (((uint64_t)e * NS + c) * NDOF + i);
   return hval(key ^ 0x2222222222222222ULL);
}

// Sparse 9x9 flux Jacobians: fixed pattern, 10 nnz per direction.
struct JacNnz { int co, ci; double a; };
static std::array<std::array<JacNnz, NNZJ>, NDIR> BuildJacobians()
{
   std::array<std::array<JacNnz, NNZJ>, NDIR> J;
   for (int d = 0; d < NDIR; d++)
   {
      int k = 0;
      for (int r = 0; r < 6; r++)          // stress rows <- velocity cols
      {
         J[d][k].co = r; J[d][k].ci = 6 + ((r + d) % 3); k++;
      }
      for (int r = 0; r < 3; r++)          // velocity rows <- stress cols
      {
         J[d][k].co = 6 + r; J[d][k].ci = (r + 2 * d) % 6; k++;
      }
      J[d][k].co = 6 + (d % 3); J[d][k].ci = (d + 3) % 6; k++; // 10th
      for (int m = 0; m < NNZJ; m++)
      {
         uint64_t key = (uint64_t)d * 64 + m;
         J[d][m].a = 0.3 * hval(key ^ 0x3333333333333333ULL);
      }
   }
   return J;
}

// Taylor / time-integral weights (levels 0..NLEV)
static void BuildWeights(double *w_nod, double *w_int)
{
   double fact = 1.0, dtk = 1.0;
   for (int k = 0; k <= NLEV; k++)
   {
      if (k > 0) { fact *= k; dtk *= DT; }
      w_nod[k] = dtk / fact;                    // dt^k / k!
      w_int[k] = dtk * DT / (fact * (k + 1));   // dt^(k+1) / (k+1)!
   }
}

// canonical-order checksum: iterate (e, c, i) identically for every layout
template <typename IDX>
static double Checksum(const double *I, const double *N, IDX idx)
{
   double s = 0.0;
   for (int e = 0; e < NE; e++)
      for (int c = 0; c < NS; c++)
         for (int i = 0; i < NDOF; i++)
         {
            const long t = idx(e, c, i);
            s += I[t] + N[t];
         }
   return s;
}

// median-of-5 timing with 2 warmup reps
template <typename F>
static double TimeMedian(const char *name, F &&run)
{
   run(); run(); // warmup (steady state: caches, allocator, BLAS init)
   std::array<double, 5> t;
   for (int r = 0; r < 5; r++)
   {
      auto t0 = std::chrono::steady_clock::now();
      run();
      auto t1 = std::chrono::steady_clock::now();
      t[r] = std::chrono::duration<double>(t1 - t0).count();
   }
   std::sort(t.begin(), t.end());
   std::printf("  %s reps (s): %.4f %.4f %.4f %.4f %.4f  -> median %.4f\n",
               name, t[0], t[1], t[2], t[3], t[4], t[2]);
   return t[2];
}

static void Report(const char *name, double sec, double secA, double chk)
{
   const double us_per_elem = sec / NE * 1e6;
   const double gflops = FLOPS_PER_ELEM * NE / sec / 1e9;
   std::printf("%-14s %10.3f us/elem  %8.2f GFLOP/s/core  speedup vs A: %6.2fx"
               "  checksum % .15e\n",
               name, us_per_elem, gflops, secA / sec, chk);
}

// ============================================================================
int main()
{
   std::printf("bench_ader_kernel_variants: NE=%d ndof=%d NUM_STATE=%d "
               "dirs=%d CK levels=%d nnzJ=%d tile=%d\n",
               NE, NDOF, NS, NDIR, NLEV, NNZJ, TILE);
   std::printf("FLOPs/elem used for GFLOP/s: %.0f "
               "(matvec 64800 + Jacobian 3600 + accum 2880)\n\n",
               FLOPS_PER_ELEM);

   const auto J = BuildJacobians();
   double w_nod[NLEV + 1], w_int[NLEV + 1];
   BuildWeights(w_nod, w_int);

   // ------------------------------------------------------------------
   // Variant A data: per-element heap DenseMatrix ops, component-major Q
   // ------------------------------------------------------------------
   std::vector<std::array<DenseMatrix, NDIR>> opsA(NE);
   for (int e = 0; e < NE; e++)
      for (int d = 0; d < NDIR; d++)
      {
         opsA[e][d].SetSize(NDOF, NDOF); // separate heap allocation each
         for (int i = 0; i < NDOF; i++)
            for (int j = 0; j < NDOF; j++)
            {
               opsA[e][d](i, j) = opval(e, d, i, j);
            }
      }

   Vector Q0cm(NS * NT); // component-major Q[c*NT + e*NDOF + i]
   for (int c = 0; c < NS; c++)
      for (int e = 0; e < NE; e++)
         for (int i = 0; i < NDOF; i++)
         {
            Q0cm[c * NT + e * NDOF + i] = qval(e, c, i);
         }

   Vector Qk(NS * NT), Qn(NS * NT), Iv(NS * NT), Nv(NS * NT);
   std::array<Vector, NDIR> dQ;
   for (int d = 0; d < NDIR; d++) { dQ[d].SetSize(NS * NT); }

   // ------------------------------------------------------------------
   // Variant A: current-structure replica
   // ------------------------------------------------------------------
   auto runA = [&]()
   {
      const double *q0 = Q0cm.GetData();
      double *qk = Qk.GetData(), *qn = Qn.GetData();
      double *iv = Iv.GetData(), *nv = Nv.GetData();
      Qk = Q0cm; // predictor starts from state
      // level-0 accumulation (whole-vector)
      for (long t = 0; t < (long)NS * NT; t++)
      {
         iv[t] = w_int[0] * q0[t];
         nv[t] = w_nod[0] * q0[t];
      }
      for (int k = 1; k <= NLEV; k++)
      {
         // 3 whole-vector spatial-derivative passes
         // (replicates wave_operator.inl Cached branch: column-major D,
         //  transposed inner loop s += D(i,j)*Qc[j])
         for (int d = 0; d < NDIR; d++)
         {
            double *dq = dQ[d].GetData();
            for (int e = 0; e < NE; e++)
            {
               const DenseMatrix &D = opsA[e][d];
               const int off = e * NDOF;
               for (int c = 0; c < NS; c++)
               {
                  const double *Qc = qk + c * NT + off;
                  double *dQc = dq + c * NT + off;
                  for (int i = 0; i < NDOF; i++)
                  {
                     double s = 0.0;
                     for (int j = 0; j < NDOF; j++) { s += D(i, j) * Qc[j]; }
                     dQc[i] = s;
                  }
               }
            }
         }
         // sparse-Jacobian strided whole-vector AXPY pass
         for (long t = 0; t < (long)NS * NT; t++) { qn[t] = 0.0; }
         for (int d = 0; d < NDIR; d++)
         {
            const double *dq = dQ[d].GetData();
            for (int m = 0; m < NNZJ; m++)
            {
               const double a = -J[d][m].a;
               double *y = qn + (long)J[d][m].co * NT;
               const double *x = dq + (long)J[d][m].ci * NT;
               for (long t = 0; t < NT; t++) { y[t] += a * x[t]; }
            }
         }
         // whole-vector accumulation AXPYs (nodal + integral)
         const double wi = w_int[k], wn = w_nod[k];
         for (long t = 0; t < (long)NS * NT; t++) { iv[t] += wi * qn[t]; }
         for (long t = 0; t < (long)NS * NT; t++) { nv[t] += wn * qn[t]; }
         std::swap(qk, qn); // vector-handle swap, as the driver does
      }
   };

   const double secA = TimeMedian("A", runA);
   const double chkA = Checksum(Iv.GetData(), Nv.GetData(),
                                [](int e, int c, int i)
                                { return (long)c * NT + (long)e * NDOF + i; });

   // ------------------------------------------------------------------
   // Variant B data: contiguous DenseTensor per direction, row-major pack
   // Drow[d](j,i,e) = D_d^e(i,j) -> unit-stride inner loop over j
   // ------------------------------------------------------------------
   std::array<DenseTensor, NDIR> Drow;
   for (int d = 0; d < NDIR; d++)
   {
      Drow[d].SetSize(NDOF, NDOF, NE);
      for (int e = 0; e < NE; e++)
         for (int i = 0; i < NDOF; i++)
            for (int j = 0; j < NDOF; j++)
            {
               Drow[d](j, i, e) = opval(e, d, i, j);
            }
   }

   // ------------------------------------------------------------------
   // Variant B: tiled-fused (the plan's design)
   // ------------------------------------------------------------------
   std::vector<double> qpan_buf(TILE * PANEL), ipan_buf(TILE * PANEL),
       npan_buf(TILE * PANEL); // tile buffers, allocated once
   auto runB = [&]()
   {
      const double *q0 = Q0cm.GetData();
      double *iv = Iv.GetData(), *nv = Nv.GetData();
      std::vector<double> &qpan = qpan_buf, &ipan = ipan_buf, &npan = npan_buf;
      double qk[PANEL], qn[PANEL], dq[NDIR][PANEL];
      for (int e0 = 0; e0 < NE; e0 += TILE)
      {
         const int nel = std::min(TILE, NE - e0);
         // gather: contiguous (20x9) panel per element, once per tile
         for (int t = 0; t < nel; t++)
         {
            const int e = e0 + t;
            double *qp = qpan.data() + t * PANEL;
            for (int c = 0; c < NS; c++)
               for (int i = 0; i < NDOF; i++)
               {
                  qp[c * NDOF + i] = q0[c * NT + e * NDOF + i];
               }
         }
         // per element: ALL 3 CK levels back-to-back, operators read once
         for (int t = 0; t < nel; t++)
         {
            const int e = e0 + t;
            const double *qp = qpan.data() + t * PANEL;
            double *ip = ipan.data() + t * PANEL;
            double *np = npan.data() + t * PANEL;
            for (int p = 0; p < PANEL; p++)
            {
               qk[p] = qp[p];
               ip[p] = w_int[0] * qp[p];
               np[p] = w_nod[0] * qp[p];
            }
            const double *De[NDIR];
            for (int d = 0; d < NDIR; d++)
            {
               De[d] = Drow[d].Data() + (long)e * NDOF * NDOF;
            }
            for (int k = 1; k <= NLEV; k++)
            {
               for (int d = 0; d < NDIR; d++)
               {
                  const double *Dt = De[d];
                  double *dqd = dq[d];
                  for (int c = 0; c < NS; c++)
                  {
                     const double *qc = qk + c * NDOF;
                     double *dc = dqd + c * NDOF;
                     for (int i = 0; i < NDOF; i++)
                     {
                        const double *Di = Dt + i * NDOF; // unit stride in j
                        double s = 0.0;
                        for (int j = 0; j < NDOF; j++) { s += Di[j] * qc[j]; }
                        dc[i] = s;
                     }
                  }
               }
               // fused Jacobian + both factorial accumulations
               for (int p = 0; p < PANEL; p++) { qn[p] = 0.0; }
               for (int d = 0; d < NDIR; d++)
               {
                  for (int m = 0; m < NNZJ; m++)
                  {
                     const double a = -J[d][m].a;
                     double *y = qn + J[d][m].co * NDOF;
                     const double *x = dq[d] + J[d][m].ci * NDOF;
                     for (int i = 0; i < NDOF; i++) { y[i] += a * x[i]; }
                  }
               }
               const double wi = w_int[k], wn = w_nod[k];
               for (int p = 0; p < PANEL; p++)
               {
                  ip[p] += wi * qn[p];
                  np[p] += wn * qn[p];
                  qk[p] = qn[p];
               }
            }
         }
         // scatter once at tile end
         for (int t = 0; t < nel; t++)
         {
            const int e = e0 + t;
            const double *ip = ipan.data() + t * PANEL;
            const double *np = npan.data() + t * PANEL;
            for (int c = 0; c < NS; c++)
               for (int i = 0; i < NDOF; i++)
               {
                  iv[c * NT + e * NDOF + i] = ip[c * NDOF + i];
                  nv[c * NT + e * NDOF + i] = np[c * NDOF + i];
               }
         }
      }
   };

   const double secB = TimeMedian("B", runB);
   const double chkB = Checksum(Iv.GetData(), Nv.GetData(),
                                [](int e, int c, int i)
                                { return (long)c * NT + (long)e * NDOF + i; });

   // ------------------------------------------------------------------
   // Variant C data: column-major DenseTensor ops + panel-major state
   // x[i + c*20 + e*180]  == BatchedLinAlg shape (20, 9, NE)
   // ------------------------------------------------------------------
   std::array<DenseTensor, NDIR> Dcol;
   for (int d = 0; d < NDIR; d++)
   {
      Dcol[d].SetSize(NDOF, NDOF, NE);
      for (int e = 0; e < NE; e++)
         for (int i = 0; i < NDOF; i++)
            for (int j = 0; j < NDOF; j++)
            {
               Dcol[d](i, j, e) = opval(e, d, i, j);
            }
   }
   Vector Q0p((long)NE * PANEL); // panel-major (layout conversion outside timer)
   for (int e = 0; e < NE; e++)
      for (int c = 0; c < NS; c++)
         for (int i = 0; i < NDOF; i++)
         {
            Q0p[e * PANEL + c * NDOF + i] = qval(e, c, i);
         }
   Vector Qkp((long)NE * PANEL), Qnp((long)NE * PANEL),
          Ivp((long)NE * PANEL), Nvp((long)NE * PANEL);
   std::array<Vector, NDIR> dQp;
   for (int d = 0; d < NDIR; d++) { dQp[d].SetSize((long)NE * PANEL); }

   MFEM_VERIFY(BatchedLinAlg::IsAvailable(BatchedLinAlg::NATIVE),
               "NATIVE batched backend unavailable");
   const BatchedLinAlgBase &batched = BatchedLinAlg::Get(BatchedLinAlg::NATIVE);

   auto runC = [&]()
   {
      const double *q0 = Q0p.GetData();
      double *iv = Ivp.GetData(), *nv = Nvp.GetData();
      Qkp = Q0p;
      Vector *qkv = &Qkp, *qnv = &Qnp;
      for (long t = 0; t < (long)NE * PANEL; t++)
      {
         iv[t] = w_int[0] * q0[t];
         nv[t] = w_nod[0] * q0[t];
      }
      for (int k = 1; k <= NLEV; k++)
      {
         // whole-batch (20,9,NE) matvec per direction per level
         for (int d = 0; d < NDIR; d++)
         {
            batched.Mult(Dcol[d], *qkv, dQp[d]);
         }
         // Jacobian pass, separate (panel layout)
         double *qn = qnv->GetData();
         for (long t = 0; t < (long)NE * PANEL; t++) { qn[t] = 0.0; }
         for (int d = 0; d < NDIR; d++)
         {
            const double *dq = dQp[d].GetData();
            for (int m = 0; m < NNZJ; m++)
            {
               const double a = -J[d][m].a;
               const int co = J[d][m].co * NDOF, ci = J[d][m].ci * NDOF;
               for (int e = 0; e < NE; e++)
               {
                  double *y = qn + (long)e * PANEL + co;
                  const double *x = dq + (long)e * PANEL + ci;
                  for (int i = 0; i < NDOF; i++) { y[i] += a * x[i]; }
               }
            }
         }
         const double wi = w_int[k], wn = w_nod[k];
         for (long t = 0; t < (long)NE * PANEL; t++) { iv[t] += wi * qn[t]; }
         for (long t = 0; t < (long)NE * PANEL; t++) { nv[t] += wn * qn[t]; }
         std::swap(qkv, qnv); // handle swap only
      }
   };

   const double secC = TimeMedian("C", runC);
   const double chkC = Checksum(Ivp.GetData(), Nvp.GetData(),
                                [](int e, int c, int i)
                                { return (long)e * PANEL + c * NDOF + i; });

   // ------------------------------------------------------------------
   // Variant D: per-element LAPACK dgemm (20x20)x(20x9), rest fused like B
   // ------------------------------------------------------------------
   auto runD = [&]()
   {
      const double *q0 = Q0cm.GetData();
      double *iv = Iv.GetData(), *nv = Nv.GetData();
      double qkbuf[PANEL], qn[PANEL], dqbuf[NDIR][PANEL];
      DenseMatrix Dw, Qw(qkbuf, NDOF, NS);
      std::array<DenseMatrix, NDIR> dQw;
      for (int d = 0; d < NDIR; d++)
      {
         dQw[d].UseExternalData(dqbuf[d], NDOF, NS);
      }
      for (int e0 = 0; e0 < NE; e0 += TILE)
      {
         const int nel = std::min(TILE, NE - e0);
         for (int t = 0; t < nel; t++)
         {
            const int e = e0 + t;
            double ip[PANEL], np[PANEL];
            for (int c = 0; c < NS; c++)
               for (int i = 0; i < NDOF; i++)
               {
                  const double q = q0[c * NT + e * NDOF + i];
                  qkbuf[c * NDOF + i] = q;
                  ip[c * NDOF + i] = w_int[0] * q;
                  np[c * NDOF + i] = w_nod[0] * q;
               }
            for (int k = 1; k <= NLEV; k++)
            {
               for (int d = 0; d < NDIR; d++)
               {
                  Dw.UseExternalData(Dcol[d].Data() + (long)e * NDOF * NDOF,
                                     NDOF, NDOF);
                  mfem::Mult(Dw, Qw, dQw[d]); // dgemm_ (MFEM_USE_LAPACK)
                  Dw.ClearExternalData();
               }
               for (int p = 0; p < PANEL; p++) { qn[p] = 0.0; }
               for (int d = 0; d < NDIR; d++)
                  for (int m = 0; m < NNZJ; m++)
                  {
                     const double a = -J[d][m].a;
                     double *y = qn + J[d][m].co * NDOF;
                     const double *x = dqbuf[d] + J[d][m].ci * NDOF;
                     for (int i = 0; i < NDOF; i++) { y[i] += a * x[i]; }
                  }
               const double wi = w_int[k], wn = w_nod[k];
               for (int p = 0; p < PANEL; p++)
               {
                  ip[p] += wi * qn[p];
                  np[p] += wn * qn[p];
                  qkbuf[p] = qn[p];
               }
            }
            for (int c = 0; c < NS; c++)
               for (int i = 0; i < NDOF; i++)
               {
                  iv[c * NT + e * NDOF + i] = ip[c * NDOF + i];
                  nv[c * NT + e * NDOF + i] = np[c * NDOF + i];
               }
         }
      }
   };

   const double secD = TimeMedian("D", runD);
   const double chkD = Checksum(Iv.GetData(), Nv.GetData(),
                                [](int e, int c, int i)
                                { return (long)c * NT + (long)e * NDOF + i; });

   // ------------------------------------------------------------------
   std::printf("\n");
   Report("A current",   secA, secA, chkA);
   Report("B tiled",     secB, secA, chkB);
   Report("C batched",   secC, secA, chkC);
   Report("D dgemm",     secD, secA, chkD);
   std::printf("\nB>=3x A: %s   C beats A: %s\n",
               (secA / secB >= 3.0) ? "YES" : "NO",
               (secC < secA) ? "YES" : "NO");
   return 0;
}

// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// nitsche_lts_setup.hpp — Nitsche profiling plan item 7 (clustered local
// time stepping), phase P0 of the design lts-kernels-design-rev1 (§3.7): the
// cluster report.
//
// Everything here works on the serial mesh before ParMesh, because cluster
// ids must be a pure function of it (lts_clustering.hpp:12-15).
//
//   BuildNitscheLtsMesh: per-element CFL step
//
//        dt_e = cfl_eff h_e / (p^2 c_p),   cfl_eff = cfl C(p),
//
//     with h_e = 6 V / A of NitscheElementLength and the constant c_p, the
//     expression of ComputeMaxDt, so min_e dt_e is the run's dt.  Plus the
//     interior face pairs (fault faces included) and the fault face pairs.
//
//   PrintNitscheLtsReport: the rate-2 clustering of BuildLtsClustering
//     (SeisSol's rules: face neighbours differ by at most one cluster, the two
//     tets of a fault face share one) in the variants the design compares:
//     raw binning, the rules, the wiggle scan, hard caps, the library's
//     default options, the fault-adjacent tets forced into cluster 0.  Then
//     fault faces per cluster, the seam-face share, the block-count model of
//     design §7.6 (lts_seam_sums.py), and, for given rank counts, the
//     per-cluster imbalance of MFEM's default cut and of the multi-constraint
//     cut of BuildLtsAwarePartition, with a tick model of the parallel gain.
//
// "Ideal" is element updates per simulated second, global stepping at
// min_e dt_e over clustered stepping: lambda N / sum_c N_c 2^-c, which is
// HarmonicUpdateSpeedup times the wiggle factor lambda (cluster 0 steps at
// lambda min_e dt_e).  Information only: nothing here changes a run.

#ifndef MFEM_SEAS_NITSCHE_LTS_SETUP_HPP
#define MFEM_SEAS_NITSCHE_LTS_SETUP_HPP

#include "mfem.hpp"
#include "nitsche_cfl.hpp"
#include "../dynamic/lts_clustering.hpp"
#include "../dynamic/lts_partition.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mfem
{
namespace seas
{

#ifdef MFEM_USE_MPI

struct NitscheLtsMesh
{
   std::vector<double> h;                         ///< 6 V / A per element (m)
   std::vector<double> dt;                        ///< CFL step per element (s)
   std::vector<std::pair<int, int>> faces;        ///< interior faces (e1, e2)
   std::vector<std::pair<int, int>> fault_pairs;  ///< fault faces (e1, e2)
};

/// Per-element CFL step, face pairs and fault pairs of the serial mesh.
inline NitscheLtsMesh BuildNitscheLtsMesh(Mesh &mesh, int order, real_t cfl_eff,
                                          real_t cp, int fault_attr)
{
   NitscheLtsMesh m;
   const int ne = mesh.GetNE();
   const real_t p2 = real_t(order) * real_t(order);
   m.h.resize(ne);
   m.dt.resize(ne);
   for (int e = 0; e < ne; e++)
   {
      m.h[e] = NitscheElementLength(mesh, e);
      m.dt[e] = cfl_eff * m.h[e] / (p2 * cp);
   }
   for (int f = 0; f < mesh.GetNumFaces(); f++)
   {
      int e1, e2;
      mesh.GetFaceElements(f, &e1, &e2);
      if (e1 >= 0 && e2 >= 0) { m.faces.emplace_back(e1, e2); }
   }
   for (int be = 0; be < mesh.GetNBE(); be++)
   {
      if (mesh.GetBdrAttribute(be) != fault_attr) { continue; }
      int f, o, e1, e2;
      mesh.GetBdrElementFace(be, &f, &o);
      mesh.GetFaceElements(f, &e1, &e2);
      MFEM_VERIFY(e1 >= 0 && e2 >= 0, "BuildNitscheLtsMesh: fault boundary element "
                  << be << " lies on face " << f << " with one element; the "
                  "Nitsche fault must be an interior face");
      m.fault_pairs.emplace_back(e1, e2);
   }
   return m;
}

namespace lts_report_detail
{

/// lambda N / sum_c N_c 2^-c (see the file comment).
inline double Ideal(const LtsClustering &cl)
{
   return cl.lambda * HarmonicUpdateSpeedup(cl);
}

inline std::string Hist(const LtsClustering &cl)
{
   std::string s;
   for (long long n : cl.cells_per_cluster()) { s += " " + std::to_string(n); }
   return s;
}

/// numpy's default (linear) percentile of an unsorted copy.
inline double Percentile(std::vector<double> v, double q)
{
   std::sort(v.begin(), v.end());
   const double pos = q * double(v.size() - 1);
   const std::size_t lo = std::size_t(pos);
   const std::size_t hi = std::min(lo + 1, v.size() - 1);
   return v[lo] + (pos - double(lo)) * (v[hi] - v[lo]);
}

/// Per-rank tets of each cluster for one cut, with its imbalance figures and
/// the tick model: every tick waits for the slowest rank, and a tick costs
/// rank r the tets of the clusters that end a step there.  Over one sync
/// interval of 2^L ticks (L = Nc - 1), tick m ends clusters 0..k with k the
/// trailing zero bits of m + 1 (at most L), so the interval costs
/// sum_{k<L} 2^(L-k-1) M_k + M_L, M_k = max_r sum_{c<=k} N_{c,r}.  Global
/// stepping on MFEM's cut costs 2^L lambda max_r N_r over the same time.
inline void PrintCut(std::ostream &os, const char *name, const std::vector<int> &part,
                     int nparts, const LtsClustering &cl,
                     const std::vector<std::pair<int, int>> &fault_pairs, double gts_max)
{
   const int nc = cl.num_clusters, L = nc - 1;
   std::vector<double> N(std::size_t(nc) * nparts, 0.0);   // N[c * nparts + r]
   for (std::size_t e = 0; e < cl.cluster.size(); e++)
   {
      N[std::size_t(cl.cluster[e]) * nparts + part[e]] += 1.0;
   }
   const LtsPartitionImbalance imb =
      ComputeLtsPartitionImbalance(cl.cluster, nc, part, nparts, nullptr, fault_pairs);
   std::vector<double> cum(nparts, 0.0);
   double interval = 0.0;
   for (int k = 0; k <= L; k++)
   {
      double mk = 0.0;
      for (int r = 0; r < nparts; r++)
      {
         cum[r] += N[std::size_t(k) * nparts + r];
         mk = std::max(mk, cum[r]);
      }
      interval += (k < L ? std::ldexp(1.0, L - k - 1) : 1.0) * mk;
   }
   double tmax = 0.0;
   for (double x : cum) { tmax = std::max(tmax, x); }
   const double mean = double(cl.cluster.size()) / nparts;
   os << "    " << name << ": tets per rank max/mean " << tmax / mean
      << "; per cluster (tets per rank, max/mean):";
   const auto cells = cl.cells_per_cluster();
   for (int c = 0; c < nc; c++)
   {
      os << " " << double(cells[c]) / nparts << ", " << imb.per_cluster[c] << ";";
   }
   os << " fault tets max/mean " << imb.fault << "; tick model "
      << std::ldexp(1.0, L) * cl.lambda * gts_max / interval << "x\n";
}

} // namespace lts_report_detail

/// The P0 report (see the file comment).  cap: the cluster cap whose
/// clustering the partition section uses (<= 0: no cap); nparts_list: rank
/// counts for the partition section (empty: no partition section).
inline void PrintNitscheLtsReport(std::ostream &os, Mesh &mesh, const NitscheLtsMesh &m,
                                  int cap, const std::vector<int> &nparts_list)
{
   using namespace lts_report_detail;
   const int ne = int(m.dt.size());
   const Table &e2e = mesh.ElementToElementTable();
   Table no_adj;   // no neighbours: BuildLtsClustering then only bins
   no_adj.MakeI(ne);
   no_adj.MakeJ();
   no_adj.ShiftUpI();

   auto build = [&](const std::vector<double> &dt, const Table &adj,
                    const std::vector<std::pair<int, int>> &fp, bool scan,
                    int nc_cap, double tol) -> LtsClustering
   {
      LtsClusteringOptions opt;
      opt.wiggle_scan = scan;
      opt.lambda_fixed = 1.0;
      opt.nc_cap = nc_cap;
      opt.merge_loss_tol = tol;
      return BuildLtsClustering(dt, adj, fp, opt);
   };
   // a hard cap: merge the top cluster down whatever it costs
   const double kAnyCost = std::numeric_limits<double>::infinity();
   const std::vector<std::pair<int, int>> none;

   const double h_min = *std::min_element(m.h.begin(), m.h.end());
   const double dt_min = *std::min_element(m.dt.begin(), m.dt.end());
   std::vector<char> on_fault(ne, 0);
   for (const auto &pr : m.fault_pairs) { on_fault[pr.first] = on_fault[pr.second] = 1; }
   std::vector<double> h_fault;
   for (int e = 0; e < ne; e++) { if (on_fault[e]) { h_fault.push_back(m.h[e]); } }

   os << std::setprecision(4)
      << "=== LTS report (plan item 7, P0; information only) ===\n"
      << "  mesh: " << ne << " tets, " << m.faces.size() << " interior faces, "
      << m.fault_pairs.size() << " fault faces\n"
      << "  h = 6V/A: min " << h_min << " m, median " << Percentile(m.h, 0.5) << " m";
   if (!h_fault.empty())
   {
      os << "; fault-adjacent tets " << h_fault.size() << " (" << 100.0 * h_fault.size() / ne
         << "%): h min " << *std::min_element(h_fault.begin(), h_fault.end()) << " m (= "
         << *std::min_element(h_fault.begin(), h_fault.end()) / h_min << " h_min), 1% "
         << Percentile(h_fault, 0.01) << " m, median " << Percentile(h_fault, 0.5) << " m";
   }
   os << "\n  min dt_e " << dt_min << " s (the run's dt)\n";

   const LtsClustering raw = build(m.dt, no_adj, none, false, 0, 0.0);
   os << "  raw rate-2 binning: " << raw.num_clusters << " clusters, ideal " << Ideal(raw)
      << "x; per cluster" << Hist(raw) << "\n";
   const LtsClustering rules = build(m.dt, e2e, m.fault_pairs, false, 0, 0.0);
   os << "  rules (neighbours differ by <= 1, fault faces in one cluster): "
      << rules.num_clusters << " clusters, ideal " << Ideal(rules) << "x; per cluster"
      << Hist(rules) << "\n";
   const LtsClustering wiggle = build(m.dt, e2e, m.fault_pairs, true, 0, 0.0);
   os << "  rules + wiggle scan: lambda " << wiggle.lambda << ", " << wiggle.num_clusters
      << " clusters, ideal " << Ideal(wiggle) << "x; per cluster" << Hist(wiggle) << "\n";
   os << "  rules, hard cap:";
   for (int c : {2, 3, 4, 6})
   {
      os << " Nc <= " << c << " " << Ideal(build(m.dt, e2e, m.fault_pairs, false, c, kAnyCost))
         << "x;";
   }
   os << "\n";
   const LtsClustering lib = BuildLtsClustering(m.dt, e2e, m.fault_pairs);
   os << "  library defaults (scan, nc_cap 6 within 5% of the cost): lambda " << lib.lambda
      << ", " << lib.num_clusters << " clusters, ideal " << Ideal(lib) << "x\n";
   {
      std::vector<double> dt0 = m.dt;
      for (int e = 0; e < ne; e++) { if (on_fault[e]) { dt0[e] = dt_min; } }
      const LtsClustering f0 = build(dt0, e2e, m.fault_pairs, false, 0, 0.0);
      os << "  fault-adjacent tets forced into cluster 0: ideal " << Ideal(f0)
         << "x; per cluster" << Hist(f0) << "\n";
   }
   {
      std::vector<long long> nf(rules.num_clusters, 0);
      for (const auto &pr : m.fault_pairs) { nf[rules.cluster[pr.first]]++; }
      os << "  fault faces per cluster (rules):";
      for (long long n : nf) { os << " " << n; }
      os << "\n";
   }

   // Block-count model of design §7.6 on the rules clustering: element
   // blocks at their cluster's rate, an intra-cluster face's two off-diagonal
   // blocks at its rate, each side of a seam face its half-form blocks at its
   // own rate; global stepping is one element block and two face blocks per
   // element and face every dt_0.
   {
      double se = 0.0, si = 0.0, ss = 0.0;
      long long nseam = 0;
      for (int c : rules.cluster) { se += std::ldexp(1.0, -c); }
      for (const auto &pr : m.faces)
      {
         const int c1 = rules.cluster[pr.first], c2 = rules.cluster[pr.second];
         if (c1 == c2) { si += std::ldexp(1.0, -c1); }
         else { ss += std::ldexp(1.0, -c1) + std::ldexp(1.0, -c2); nseam++; }
      }
      const double nface = double(m.faces.size());
      os << "  seam faces (two clusters): " << 100.0 * nseam / nface
         << "% of interior faces; block-count model:";
      for (int p : {1, 2})
      {
         const int np = (p + 1) * (p + 2) * (p + 3) / 6, nfn = (p + 1) * (p + 2) / 2;
         const double n = 3.0 * np, d = n * n;
         // full: dense 3np x 3np blocks as the CSR stores them; trace: the
         // face-node pattern (off-diagonal (3np)^2 - (3(np-nf))^2, seam side
         // 2 (3nf)(3np) + 3 nf^2)
         const double o_full = n * n, s_full = 3.0 * n * n;
         const double o_tr = n * n - std::pow(3.0 * (np - nfn), 2);
         const double s_tr = 2.0 * (3.0 * nfn) * n + 3.0 * nfn * nfn;
         auto speed = [&](double o, double s)
         { return (ne * d + nface * 2.0 * o) / (se * d + si * 2.0 * o + ss * s); };
         os << " p=" << p << " trace " << speed(o_tr, s_tr) << "x, full "
            << speed(o_full, s_full) << "x;";
      }
      os << "\n";
   }

   if (nparts_list.empty()) { return; }
   const LtsClustering pc = (cap > 0) ? build(m.dt, e2e, m.fault_pairs, false, cap, kAnyCost)
                                      : rules;
   os << "  partition, rules clustering" << (cap > 0 ? " with hard cap " + std::to_string(cap)
                                                    : std::string(", no cap"))
      << " (" << pc.num_clusters << " clusters, ideal " << Ideal(pc) << "x):\n";
   std::vector<int> xadj(e2e.GetI(), e2e.GetI() + ne + 1);
   std::vector<int> adj(e2e.GetJ(), e2e.GetJ() + e2e.Size_of_connections());
   for (int nparts : nparts_list)
   {
      int *mp = mesh.GeneratePartitioning(nparts);   // ParMesh's default cut
      std::vector<int> part(mp, mp + ne);
      delete [] mp;
      std::vector<double> nr(nparts, 0.0);
      for (int r : part) { nr[r] += 1.0; }
      const double gts_max = *std::max_element(nr.begin(), nr.end());
      os << "  " << nparts << " ranks (" << double(ne) / nparts << " tets per rank):\n";
      PrintCut(os, "MFEM's default cut", part, nparts, pc, m.fault_pairs, gts_max);
      try
      {
         LtsPartitionOptions popt;
         BuildLtsAwarePartition(ne, xadj, adj, pc.cluster, pc.num_clusters, nparts,
                                nullptr, {}, popt, part);
         PrintCut(os, "multi-constraint cut", part, nparts, pc, m.fault_pairs, gts_max);
      }
      catch (const std::exception &ex)
      {
         os << "    multi-constraint cut: failed: " << ex.what() << "\n";
      }
   }
}

#endif // MFEM_USE_MPI

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_NITSCHE_LTS_SETUP_HPP

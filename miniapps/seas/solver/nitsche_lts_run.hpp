// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// nitsche_lts_run.hpp — Nitsche profiling plan item 7 (clustered local time
// stepping), phases P3 and P4 of the design lts-kernels-design-rev1 (§3.4,
// §3.5, §3.7): the multirate leapfrog of nitsche_lts_stepper.hpp driving a
// NitscheWaveOperator with its fault, on any number of ranks.  The driver
// uses it with --lts rate2 (or [numerics] lts = "rate2"); off by default, and
// a run with it off does not reach this file.
//
// BuildNitscheLtsAssignment (serial mesh, one rank): per-element steps dt_e
// of BuildNitscheLtsMesh (capped by dt_max), the rate-2 clustering of
// BuildLtsClustering with the rules (face neighbours differ by at most one
// cluster, both tets of a fault face in one cluster), no wiggle (cluster 0
// steps at min_e dt_e, the global step) and a hard cluster cap.  Then the
// whole fault goes into one cluster c_F, the fastest one any fault tet is in
// (design §3.4 step 1): every fault tet takes dt_e = min(dt_e, dt0 2^c_F) and
// the clustering is repeated until no rule moves a fault tet below c_F.
//
// NitscheLtsBlocks: K^(c) = sum_{e in c} G_e and K_eta^(c) as dense blocks
// (the layout of DenseBlockStiffness, plan item 8a), built from the
// integrators of ElasticityDomainOperator::AssembleStiffness with the element-
// owner split of NitscheLtsSplit: the volume and Dirichlet terms of e, and
// the side-e half of each face of e (DGElasticityIPCombinedIntegrator::
// SetSide).  The rows K^(c) touches are the patch of c: its elements and
// their face neighbours.  Per cluster: one upper-triangular diagonal block
// per patch element, one block per interior face with a side in c (rows of
// the lower-numbered element; the other rows are its transpose, the side
// forms being exactly symmetric), one block per shared face with a side in c
// (this rank's rows against the neighbour element's DOFs).  Fault faces are
// ordinary interior faces of K, as in the assembled K; both of their tets
// are in c_F, so their slip loads belong to c_F alone.
//
// Halo (design §3.5 "Parallel"): a cluster reads the face-neighbour elements
// of its shared faces.  Level of a neighbour element = the lowest cluster of
// the two sides over its shared faces with this rank; at a tick where
// clusters 0..k end, one merged exchange of u (and v) moves the neighbour
// elements of level <= k, posted before the fault closure and the interior
// blocks and completed before the shared blocks.  Both ranks of a pair sort
// the elements they exchange by (level, position in MFEM's face-neighbour
// lists), so a level-k message is a prefix of the full one.  The positions
// match as in ElasticityDomainOperator::FaultNbrExchange_.
//
// NitscheLtsRun: one tick of design §3.4 step 3 at a time, dt_c = 2^c dt0:
//
//    for c starting a step at t:    v[P_c] += dt_c/2 a_c         (start kicks)
//    u += dt0 v,  t += dt0                                       (drift)
//    E = {0 .. top}: the clusters ending a step at t
//    post the level-top halo of (u, v)
//    if c_F in E: AdvanceFault_(u, t, v)    (boxes 3-11 with dt = dt_F)
//    for c in E: r_c = -(K^(c) u + K_eta^(c) v) [+ G d + G_eta V if c = c_F]
//                     + F_D(t) on c's elements - C_e v on c's absorbing elements
//                a_c = B_top^{-1} r_c on P_c
//    for c in E:                    v[P_c] += dt_c/2 a_c         (end kicks)
//
// B_top = M + 1/2 dt_(c_e) C_e on the absorbing elements with c_e <= top and
// M elsewhere (nitsche_lts_stepper.hpp, P2).  The fault closes after the
// drift and before the forces, as in NitscheWaveOperator::ImplicitSolve, with
// u at t and v before the end kicks (the fault cluster's mid-step velocity
// on its own rows) in the roles of U^{n+1} and Udot*; GetTimeStep() must be
// dt_F.  Init computes every a_c with M^{-1}, as NewmarkSolver's first Mult.
// With one cluster this is central difference in kick-drift-kick form, equal
// to CentralDifferenceSolver on the operator up to round-off (P3 test (a)).
// same_step (test only, P3 test (b)): the same split and clusters, but every
// cluster steps at dt0 and ends at every tick; equal to global stepping up to
// round-off.

#ifndef MFEM_SEAS_NITSCHE_LTS_RUN_HPP
#define MFEM_SEAS_NITSCHE_LTS_RUN_HPP

#include "mfem.hpp"
#include "nitsche_lts_setup.hpp"
#include "../domain/nitsche_mass.hpp"
#include "../domain/nitsche_wave_operator.hpp"
#include "../dynamic/lts_clustering.hpp"
#include "../integrator/dg_elasticity_ip_combined_integrator.hpp"

#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

namespace mfem
{
namespace seas
{

#ifdef MFEM_USE_MPI

struct NitscheLtsAssignment
{
   std::vector<int> cluster;        ///< per serial element
   int num_clusters = 1;
   int fault_cluster = 0;           ///< c_F
   double dt0 = 0.0;                ///< cluster 0's step (s)
   double ideal = 1.0;              ///< element updates, global over clustered
   int passes = 0;                  ///< clusterings until the fault sat in one cluster
   long long fault_tets = 0;
   std::vector<long long> cells;    ///< per cluster
   std::vector<std::pair<int, int>> fault_pairs;
};

/// See the file comment.  cap: hard cluster cap (<= 0: none).
inline NitscheLtsAssignment BuildNitscheLtsAssignment(Mesh &mesh, int order, real_t cfl_eff,
                                                      real_t cp, int fault_attr, real_t dt_max,
                                                      int cap)
{
   NitscheLtsMesh m = BuildNitscheLtsMesh(mesh, order, cfl_eff, cp, fault_attr);
   const int ne = static_cast<int>(m.dt.size());
   std::vector<double> dt = m.dt;
   for (double &d : dt) { d = std::min(d, static_cast<double>(dt_max)); }
   std::vector<char> on_fault(ne, 0);
   for (const auto &pr : m.fault_pairs) { on_fault[pr.first] = on_fault[pr.second] = 1; }
   LtsClusteringOptions opt;
   opt.wiggle_scan = false;
   opt.lambda_fixed = 1.0;
   opt.nc_cap = cap;
   opt.merge_loss_tol = std::numeric_limits<double>::infinity();   // a hard cap
   const Table &e2e = mesh.ElementToElementTable();
   NitscheLtsAssignment a;
   LtsClustering cl;
   for (a.passes = 1;; a.passes++)
   {
      cl = BuildLtsClustering(dt, e2e, m.fault_pairs, opt);
      int cmin = INT_MAX, cmax = -1;
      for (int e = 0; e < ne; e++)
      {
         if (!on_fault[e]) { continue; }
         cmin = std::min(cmin, cl.cluster[e]);
         cmax = std::max(cmax, cl.cluster[e]);
      }
      a.fault_cluster = (cmax < 0) ? 0 : cmin;
      if (cmax < 0 || cmin == cmax) { break; }
      MFEM_VERIFY(a.passes < 64, "BuildNitscheLtsAssignment: the fault tets do not "
                  "settle in one cluster after " << a.passes << " clusterings");
      const double dtf = std::ldexp(cl.dt_base, cmin);
      for (int e = 0; e < ne; e++) { if (on_fault[e]) { dt[e] = std::min(dt[e], dtf); } }
   }
   a.cluster = cl.cluster;
   a.num_clusters = cl.num_clusters;
   a.dt0 = cl.dt_base;
   a.ideal = cl.lambda * HarmonicUpdateSpeedup(cl);
   a.cells = cl.cells_per_cluster();
   for (char f : on_fault) { a.fault_tets += f; }
   a.fault_pairs = m.fault_pairs;
   // every element steps no faster than its own CFL step (BuildLtsClustering
   // asserts it against the lowered dt; this is the original one)
   for (int e = 0; e < ne; e++)
   {
      MFEM_VERIFY(std::ldexp(a.dt0, a.cluster[e]) <= m.dt[e],
                  "BuildNitscheLtsAssignment: element " << e << " in cluster "
                  << a.cluster[e] << " steps at " << std::ldexp(a.dt0, a.cluster[e])
                  << " s > its CFL step " << m.dt[e] << " s");
   }
   return a;
}

/// K^(c) and K_eta^(c) as dense blocks, with the per-level halo (file comment).
class NitscheLtsBlocks
{
public:
   /// cl: cluster of each local element; cl_nbr: of each face-neighbour
   /// element (ParMesh::ExchangeFaceNbrData order); klam, kmu: the Kelvin-
   /// Voigt moduli (both null: no K_eta).  Collective.
   void Build(ParFiniteElementSpace &fes, Coefficient &lam, Coefficient &mu, Coefficient *klam,
              Coefficient *kmu, const Array<int> &dir_marker, real_t penalty,
              const std::vector<int> &cl, const std::vector<int> &cl_nbr, int nc)
   {
      MFEM_VERIFY((klam == nullptr) == (kmu == nullptr), "NitscheLtsBlocks: give both "
                  "Kelvin-Voigt moduli or neither");
      ParMesh &pm = *fes.GetParMesh();
      pm.ExchangeFaceNbrData();
      fes.ExchangeFaceNbrData();
      nm_ = klam ? 2 : 1;
      nc_ = nc;
      ne_ = pm.GetNE();
      n_ = fes.GetVSize();
      comm_ = pm.GetComm();
      MFEM_VERIFY(static_cast<int>(cl.size()) == ne_, "NitscheLtsBlocks: " << cl.size()
                  << " cluster ids for " << ne_ << " elements");
      MFEM_VERIFY(static_cast<int>(cl_nbr.size()) == pm.GetNFaceNeighborElements(),
                  "NitscheLtsBlocks: " << cl_nbr.size() << " neighbour cluster ids for "
                  << pm.GetNFaceNeighborElements() << " face-neighbour elements");
      // element DOFs: every local DOF in exactly one element, one block size
      Array<int> vd;
      nb_ = 0;
      if (ne_ > 0) { fes.GetElementVDofs(0, vd); nb_ = vd.Size(); }
      int nb_all = nb_;
      {
         int mx = nb_;
         MPI_Allreduce(MPI_IN_PLACE, &mx, 1, MPI_INT, MPI_MAX, comm_);
         nb_all = mx;
      }
      if (ne_ == 0) { nb_ = nb_all; }
      MFEM_VERIFY(nb_ == nb_all && nb_ <= kMaxNB, "NitscheLtsBlocks: block size " << nb_
                  << " (largest over the ranks " << nb_all << ", at most " << kMaxNB << ")");
      tri_ = nb_ * (nb_ + 1) / 2;
      idx_.assign(static_cast<std::size_t>(ne_) * nb_, 0);
      std::vector<char> seen(n_, 0);
      for (int e = 0; e < ne_; e++)
      {
         fes.GetElementVDofs(e, vd);
         MFEM_VERIFY(vd.Size() == nb_, "NitscheLtsBlocks: element " << e << " has "
                     << vd.Size() << " DOFs, expected " << nb_ << " (one block size)");
         for (int a = 0; a < nb_; a++)
         {
            MFEM_VERIFY(vd[a] >= 0 && vd[a] < n_ && !seen[vd[a]], "NitscheLtsBlocks: "
                        "element DOFs are not a partition of the DG space");
            seen[vd[a]] = 1;
            idx_[static_cast<std::size_t>(e) * nb_ + a] = vd[a];
         }
      }
      for (int i = 0; i < n_; i++)
      {
         MFEM_VERIFY(seen[i], "NitscheLtsBlocks: DOF " << i << " lies in no element");
      }

      // which clusters touch the rows of each element: its own, and those of
      // its face neighbours (their side halves reach its rows)
      const int nsf = pm.GetNSharedFaces();
      std::vector<unsigned> touch(ne_, 0u);
      MFEM_VERIFY(nc_ >= 1 && nc_ <= 30, "NitscheLtsBlocks: " << nc_ << " clusters");
      for (int e = 0; e < ne_; e++)
      {
         MFEM_VERIFY(cl[e] >= 0 && cl[e] < nc_, "NitscheLtsBlocks: element " << e
                     << " in cluster " << cl[e]);
         touch[e] |= 1u << cl[e];
      }
      std::vector<std::array<int, 3>> ifaces;   // interior faces: (e1, e2, face)
      for (int f = 0; f < pm.GetNumFaces(); f++)
      {
         int e1, e2;
         pm.GetFaceElements(f, &e1, &e2);
         if (e1 < 0 || e2 < 0) { continue; }
         ifaces.push_back({e1, e2, f});
         touch[e1] |= 1u << cl[e2];
         touch[e2] |= 1u << cl[e1];
      }
      std::vector<std::array<int, 3>> sfaces;   // shared faces: (e1, nbr j, shared index)
      for (int sf = 0; sf < nsf; sf++)
      {
         int e1, e2;
         pm.GetFaceElements(pm.GetSharedFace(sf), &e1, &e2);
         MFEM_VERIFY(e1 >= 0 && e2 < 0, "NitscheLtsBlocks: shared face " << sf
                     << " has elements " << e1 << ", " << e2);
         const int j = -1 - e2;
         MFEM_VERIFY(cl_nbr[j] >= 0 && cl_nbr[j] < nc_, "NitscheLtsBlocks: face-neighbour "
                     "element " << j << " in cluster " << cl_nbr[j]);
         sfaces.push_back({e1, j, sf});
         touch[e1] |= 1u << cl_nbr[j];
      }
      C_.assign(nc_, Cluster());
      std::vector<std::vector<int>> pos(nc_, std::vector<int>(ne_, -1));
      for (int e = 0; e < ne_; e++)
      {
         for (int c = 0; c < nc_; c++)
         {
            if (touch[e] & (1u << c))
            {
               pos[c][e] = static_cast<int>(C_[c].pe.size());
               C_[c].pe.push_back(e);
            }
         }
      }

      // halo: face-neighbour groups, levels, the sorted exchange lists
      BuildHalo_(pm, cl, cl_nbr, sfaces);

      // face slots per cluster: interior faces grouped by their lower
      // element (the order of the sweep), shared faces by their element
      std::vector<std::array<int, 4>> islot(ifaces.size(), {-1, -1, -1, -1});
      {
         std::vector<std::vector<std::array<int, 3>>> ent(nc_);   // (low pos, face #, high)
         for (int q = 0; q < static_cast<int>(ifaces.size()); q++)
         {
            const int e1 = ifaces[q][0], e2 = ifaces[q][1];
            const int lo = std::min(e1, e2), hi = std::max(e1, e2);
            for (int c : {cl[e1], cl[e2]})
            {
               if (!ent[c].empty() && ent[c].back()[1] == q) { continue; }   // c1 == c2
               ent[c].push_back({pos[c][lo], q, hi});
            }
         }
         for (int c = 0; c < nc_; c++)
         {
            std::sort(ent[c].begin(), ent[c].end());
            Cluster &K = C_[c];
            K.fown.assign(K.pe.size() + 1, 0);
            K.fhi.resize(ent[c].size());
            for (std::size_t s = 0; s < ent[c].size(); s++)
            {
               K.fown[ent[c][s][0] + 1]++;
               K.fhi[s] = ent[c][s][2];
               std::array<int, 4> &sl = islot[ent[c][s][1]];
               if (sl[0] < 0) { sl[0] = c; sl[1] = static_cast<int>(s); }
               else { sl[2] = c; sl[3] = static_cast<int>(s); }
            }
            for (std::size_t i = 0; i < K.pe.size(); i++) { K.fown[i + 1] += K.fown[i]; }
            K.B.assign(ent[c].size() * nb_ * nb_ * nm_, 0.0);
            K.D.assign(K.pe.size() * tri_ * nm_, 0.0);
         }
      }
      std::vector<std::array<int, 4>> sslot(sfaces.size(), {-1, -1, -1, -1});
      {
         std::vector<std::vector<std::array<int, 3>>> ent(nc_);   // (elem, face #, halo slot)
         for (int q = 0; q < static_cast<int>(sfaces.size()); q++)
         {
            const int e1 = sfaces[q][0], j = sfaces[q][1];
            for (int c : {cl[e1], cl_nbr[j]})
            {
               if (!ent[c].empty() && ent[c].back()[1] == q) { continue; }
               ent[c].push_back({e1, q, nbr_slot_[j]});
            }
         }
         for (int c = 0; c < nc_; c++)
         {
            std::sort(ent[c].begin(), ent[c].end());
            Cluster &K = C_[c];
            K.se.resize(ent[c].size());
            K.sh.resize(ent[c].size());
            for (std::size_t s = 0; s < ent[c].size(); s++)
            {
               K.se[s] = ent[c][s][0];
               K.sh[s] = ent[c][s][2];
               std::array<int, 4> &sl = sslot[ent[c][s][1]];
               if (sl[0] < 0) { sl[0] = c; sl[1] = static_cast<int>(s); }
               else { sl[2] = c; sl[3] = static_cast<int>(s); }
            }
            K.S.assign(ent[c].size() * nb_ * nb_ * nm_, 0.0);
         }
      }

      // the blocks: the integrators of AssembleStiffness, split by owner
      ElasticityIntegrator vol(lam, mu);
      DGElasticityIPCombinedIntegrator face(lam, mu, 3, -1.0, penalty);   // SIPG: epsilon -1
      std::unique_ptr<ElasticityIntegrator> kvol;
      std::unique_ptr<DGElasticityIPCombinedIntegrator> kface;
      if (nm_ == 2)
      {
         kvol = std::make_unique<ElasticityIntegrator>(*klam, *kmu);
         kface = std::make_unique<DGElasticityIPCombinedIntegrator>(*klam, *kmu, 3, -1.0, penalty);
      }
      DenseMatrix elmat;
      for (int e = 0; e < ne_; e++)
      {
         for (int m = 0; m < nm_; m++)
         {
            ElasticityIntegrator &vi = (m == 0) ? vol : *kvol;
            vi.AssembleElementMatrix(*fes.GetFE(e), *fes.GetElementTransformation(e), elmat);
            AddDiag_(cl[e], pos[cl[e]][e], elmat, 0, m);
         }
      }
      for (int be = 0; be < pm.GetNBE(); be++)
      {
         const int attr = pm.GetBdrAttribute(be);
         if (attr > dir_marker.Size() || !dir_marker[attr - 1]) { continue; }
         FaceElementTransformations *tr = pm.GetBdrFaceTransformations(be);
         if (!tr) { continue; }   // an interior face, skipped as BilinearForm does
         const int e = tr->Elem1No;
         for (int m = 0; m < nm_; m++)
         {
            DGElasticityIPCombinedIntegrator &fi = (m == 0) ? face : *kface;
            fi.SetSide(0);
            fi.AssembleFaceMatrix(*fes.GetFE(e), *fes.GetFE(e), *tr, elmat);
            AddDiag_(cl[e], pos[cl[e]][e], elmat, 0, m);
            ndir_++;
         }
      }
      for (int q = 0; q < static_cast<int>(ifaces.size()); q++)
      {
         const int e1 = ifaces[q][0], e2 = ifaces[q][1];
         FaceElementTransformations *tr = pm.GetInteriorFaceTransformations(ifaces[q][2]);
         MFEM_VERIFY(tr && tr->Elem1No == e1 && tr->Elem2No == e2, "NitscheLtsBlocks: "
                     "interior face " << ifaces[q][2] << " transformation");
         for (int s = 1; s <= 2; s++)
         {
            const int cs = cl[s == 1 ? e1 : e2];
            const int slot = (islot[q][0] == cs) ? islot[q][1] : islot[q][3];
            MFEM_VERIFY((islot[q][0] == cs || islot[q][2] == cs) && slot >= 0,
                        "NitscheLtsBlocks: no slot for interior face " << q);
            for (int m = 0; m < nm_; m++)
            {
               DGElasticityIPCombinedIntegrator &fi = (m == 0) ? face : *kface;
               fi.SetSide(s);
               fi.AssembleFaceMatrix(*fes.GetFE(e1), *fes.GetFE(e2), *tr, elmat);
               AddDiag_(cs, pos[cs][e1], elmat, 0, m);
               AddDiag_(cs, pos[cs][e2], elmat, nb_, m);
               if (e1 < e2) { AddFace_(C_[cs].B, slot, elmat, 0, nb_, m); }
               else { AddFace_(C_[cs].B, slot, elmat, nb_, 0, m); }
            }
         }
      }
      for (int q = 0; q < static_cast<int>(sfaces.size()); q++)
      {
         const int e1 = sfaces[q][0], j = sfaces[q][1];
         FaceElementTransformations *tr = pm.GetSharedFaceTransformations(sfaces[q][2]);
         MFEM_VERIFY(tr && tr->Elem1No == e1 && tr->Elem2No == ne_ + j, "NitscheLtsBlocks: "
                     "shared face " << sfaces[q][2] << " transformation");
         for (int s = 1; s <= 2; s++)
         {
            const int cs = (s == 1) ? cl[e1] : cl_nbr[j];
            const int slot = (sslot[q][0] == cs) ? sslot[q][1] : sslot[q][3];
            MFEM_VERIFY((sslot[q][0] == cs || sslot[q][2] == cs) && slot >= 0,
                        "NitscheLtsBlocks: no slot for shared face " << q);
            for (int m = 0; m < nm_; m++)
            {
               DGElasticityIPCombinedIntegrator &fi = (m == 0) ? face : *kface;
               fi.SetSide(s);
               fi.AssembleFaceMatrix(*fes.GetFE(e1), *fes.GetFaceNbrFE(j), *tr, elmat);
               AddDiag_(cs, pos[cs][e1], elmat, 0, m);
               AddFace_(C_[cs].S, slot, elmat, 0, nb_, m);
            }
         }
      }
      face.SetSide(0);
      if (kface) { kface->SetSide(0); }
      // the stored upper triangles stand for exactly symmetric blocks
      real_t g[2] = {asym_, amax_};
      MPI_Allreduce(MPI_IN_PLACE, g, 2, MPITypeMap<real_t>::mpi_type, MPI_MAX, comm_);
      asym_ = g[0];
      amax_ = g[1];
      MFEM_VERIFY(asym_ <= 1e-12 * amax_, "NitscheLtsBlocks: element or face matrix not "
                  "symmetric: largest |A_ij - A_ji| " << asym_ << ", largest entry " << amax_);
   }

   int NumClusters() const { return nc_; }
   int BlockSize() const { return nb_; }
   int NumMatrices() const { return nm_; }
   /// Local elements whose rows K^(c) touches (sorted).
   const std::vector<int> &Patch(int c) const { return C_[c].pe; }
   /// Vector DOFs of local element e (nb, byNODES as GetElementVDofs).
   const int *ElementDofs(int e) const { return idx_.data() + static_cast<std::size_t>(e) * nb_; }
   long long Bytes() const
   {
      long long b = 0;
      for (const Cluster &K : C_)
      {
         b += static_cast<long long>(sizeof(real_t)) * (K.D.size() + K.B.size() + K.S.size())
              + static_cast<long long>(sizeof(int)) * (K.pe.size() + K.fown.size() + K.fhi.size()
                                                      + K.se.size() + K.sh.size());
      }
      return b + static_cast<long long>(sizeof(int)) * idx_.size();
   }
   /// Halo rounds posted, and the values received in them, since Build.
   long long ExchangeRounds() const { return rounds_; }
   long long ExchangeValues() const { return recv_values_; }
   /// Seconds this rank spent in WaitHalo (the exposed part of the exchange).
   double WaitSeconds() const { return wait_sec_; }
   /// Face-neighbour elements this rank receives at level <= k (all ranks' max).
   int NumHaloElements() const { return nh_; }

   /// Post the exchange of the halo clusters 0..top read: u (and v for K_eta).
   void PostHalo(int top, const Vector &u, const Vector *v)
   {
      MFEM_VERIFY(!posted_, "NitscheLtsBlocks::PostHalo: the previous exchange is open");
      MFEM_VERIFY(top >= 0 && top < nc_, "NitscheLtsBlocks::PostHalo: level " << top);
      MFEM_VERIFY(nm_ == 1 || v, "NitscheLtsBlocks::PostHalo: K_eta needs v");
      const real_t *x[2] = {u.GetData(), v ? v->GetData() : nullptr};
      const int nbm = nb_ * nm_;
      const int tag = 7331;
      int q = 0;
      reqs_.resize(2 * groups_.size());
      for (const Group &g : groups_)
      {
         const int nr = g.nrecv[top], ns = g.nsend[top];
         if (nr > 0)
         {
            MPI_Irecv(recv_.data() + static_cast<std::size_t>(g.rbase) * nbm, nr * nbm,
                      MPITypeMap<real_t>::mpi_type, g.rank, tag, comm_, &reqs_[q++]);
            recv_values_ += static_cast<long long>(nr) * nbm;
         }
         if (ns > 0)
         {
            real_t *b = send_.data() + static_cast<std::size_t>(g.sbase) * nbm;
            for (int i = 0; i < ns; i++)
            {
               const int *ie = ElementDofs(g.send_elems[i]);
               for (int a = 0; a < nb_; a++)
               {
                  for (int m = 0; m < nm_; m++) { *b++ = x[m][ie[a]]; }
               }
            }
            MPI_Isend(send_.data() + static_cast<std::size_t>(g.sbase) * nbm, ns * nbm,
                      MPITypeMap<real_t>::mpi_type, g.rank, tag, comm_, &reqs_[q++]);
         }
      }
      nreq_ = q;
      posted_ = true;
      rounds_++;
   }
   void WaitHalo()
   {
      MFEM_VERIFY(posted_, "NitscheLtsBlocks::WaitHalo: nothing posted");
      const double w0 = MPI_Wtime();
      MPI_Waitall(nreq_, reqs_.data(), MPI_STATUSES_IGNORE);
      wait_sec_ += MPI_Wtime() - w0;
      posted_ = false;
   }

   /// y = K^(c) u (+ K_eta^(c) v) on the patch rows of c, from the element
   /// and interior-face blocks; the shared blocks follow in AddShared after
   /// WaitHalo.  Other rows of y are not touched.
   void MultInterior(int c, const Vector &u, const Vector *v, Vector &y) const
   {
      const real_t *x[2] = {u.GetData(), v ? v->GetData() : nullptr};
      if (nm_ == 2) { Dispatch_<2>(c, x, y.GetData(), false); }
      else { Dispatch_<1>(c, x, y.GetData(), false); }
   }
   void AddShared(int c, Vector &y) const
   {
      MFEM_VERIFY(!posted_, "NitscheLtsBlocks::AddShared: call WaitHalo first");
      if (nm_ == 2) { Dispatch_<2>(c, nullptr, y.GetData(), true); }
      else { Dispatch_<1>(c, nullptr, y.GetData(), true); }
   }

private:
   static constexpr int kMaxNB = 168;   // 3 x 56 nodes: p = 5 tets

   struct Cluster
   {
      std::vector<int> pe;            ///< patch elements, sorted
      std::vector<real_t> D;          ///< their diagonal blocks (upper, packed)
      std::vector<int> fown, fhi;     ///< faces per patch position (lower element); upper element
      std::vector<real_t> B;          ///< interior face blocks, rows of the lower element
      std::vector<int> se, sh;        ///< shared blocks: local element, halo element slot
      std::vector<real_t> S;          ///< shared blocks, rows of se
   };
   struct Group
   {
      int rank = -1;
      int rbase = 0, sbase = 0;       ///< first slot in recv_ / send_ (elements)
      std::vector<int> nrecv, nsend;  ///< elements with level <= k, per k
      std::vector<int> send_elems;    ///< local elements sent, in level order
   };

   void BuildHalo_(ParMesh &pm, const std::vector<int> &cl, const std::vector<int> &cl_nbr,
                   const std::vector<std::array<int, 3>> &sfaces)
   {
      const int ng = pm.GetNFaceNeighbors();
      const int nnbr = pm.GetNFaceNeighborElements();
      std::vector<int> grp(nnbr, -1);
      for (int g = 0; g < ng; g++)
      {
         for (int j = pm.face_nbr_elements_offset[g]; j < pm.face_nbr_elements_offset[g + 1]; j++)
         {
            grp[j] = g;
         }
      }
      std::vector<int> lev_r(nnbr, INT_MAX);
      std::vector<std::unordered_map<int, int>> lev_s(ng);   // local element -> level
      for (const auto &f : sfaces)
      {
         const int e1 = f[0], j = f[1];
         const int lv = std::min(cl[e1], cl_nbr[j]);
         lev_r[j] = std::min(lev_r[j], lv);
         auto it = lev_s[grp[j]].find(e1);
         if (it == lev_s[grp[j]].end()) { lev_s[grp[j]][e1] = lv; }
         else { it->second = std::min(it->second, lv); }
      }
      groups_.assign(ng, Group());
      nbr_slot_.assign(nnbr, -1);
      int rb = 0, sb = 0;
      for (int g = 0; g < ng; g++)
      {
         Group &G = groups_[g];
         G.rank = pm.GetFaceNbrRank(g);
         const int off = pm.face_nbr_elements_offset[g];
         const int nr = pm.face_nbr_elements_offset[g + 1] - off;
         std::vector<std::pair<int, int>> r(nr);   // (level, position)
         for (int p = 0; p < nr; p++)
         {
            const int lv = lev_r[off + p];
            MFEM_VERIFY(lv != INT_MAX, "NitscheLtsBlocks: face-neighbour element " << off + p
                        << " (rank " << G.rank << ") shares no face with this rank");
            r[p] = {lv, p};
         }
         std::sort(r.begin(), r.end());
         G.rbase = rb;
         G.nrecv.assign(nc_, 0);
         for (int i = 0; i < nr; i++)
         {
            nbr_slot_[off + r[i].second] = rb + i;
            for (int k = r[i].first; k < nc_; k++) { G.nrecv[k]++; }
         }
         rb += nr;
         const int *row = pm.send_face_nbr_elements.GetRow(g);
         const int ns = pm.send_face_nbr_elements.RowSize(g);
         std::vector<std::pair<int, int>> s(ns);
         for (int p = 0; p < ns; p++)
         {
            auto it = lev_s[g].find(row[p]);
            MFEM_VERIFY(it != lev_s[g].end(), "NitscheLtsBlocks: element " << row[p]
                        << " sent to rank " << G.rank << " shares no face with it");
            s[p] = {it->second, p};
         }
         std::sort(s.begin(), s.end());
         G.sbase = sb;
         G.nsend.assign(nc_, 0);
         G.send_elems.resize(ns);
         for (int i = 0; i < ns; i++)
         {
            G.send_elems[i] = row[s[i].second];
            for (int k = s[i].first; k < nc_; k++) { G.nsend[k]++; }
         }
         sb += ns;
      }
      nh_ = rb;
      recv_.assign(static_cast<std::size_t>(rb) * nb_ * nm_, 0.0);
      send_.assign(static_cast<std::size_t>(sb) * nb_ * nm_, 0.0);
   }

   /// Upper triangle of the (off, off) nb x nb block of M into the diagonal
   /// block at patch position i of cluster c, matrix m; the mirror entries
   /// are compared (asym_).
   void AddDiag_(int c, int i, const DenseMatrix &M, int off, int m)
   {
      MFEM_VERIFY(i >= 0, "NitscheLtsBlocks: element outside the patch of cluster " << c);
      real_t *d = C_[c].D.data() + static_cast<std::size_t>(i) * tri_ * nm_;
      for (int a = 0; a < nb_; a++)
      {
         const int len = nb_ - a;
         real_t *dr = d + static_cast<std::size_t>(nm_) * (a * nb_ - a * (a - 1) / 2) + m * len;
         for (int b = a; b < nb_; b++)
         {
            const real_t x = M(off + a, off + b);
            dr[b - a] += x;
            asym_ = std::max(asym_, std::abs(x - M(off + b, off + a)));
            amax_ = std::max(amax_, std::abs(x));
         }
      }
   }
   /// The (r0, c0) nb x nb block of M into face block `slot` of B, matrix m.
   void AddFace_(std::vector<real_t> &B, int slot, const DenseMatrix &M, int r0, int c0, int m)
   {
      for (int a = 0; a < nb_; a++)
      {
         real_t *br = B.data() + ((static_cast<std::size_t>(slot) * nb_ + a) * nm_ + m) * nb_;
         for (int b = 0; b < nb_; b++)
         {
            br[b] += M(r0 + a, c0 + b);
            amax_ = std::max(amax_, std::abs(M(r0 + a, c0 + b)));
         }
      }
   }

   // two-wide vectors, as DenseBlockStiffness
   typedef real_t V2 __attribute__((vector_size(2 * sizeof(real_t))));
   static V2 Load2_(const real_t *p) { V2 v; __builtin_memcpy(&v, p, sizeof(V2)); return v; }
   static void Store2_(real_t *p, V2 v) { __builtin_memcpy(p, &v, sizeof(V2)); }

   template <int NBC>
   static real_t RowDotAxpy_(int nb, const real_t *r, const real_t *x, real_t xa, real_t *y)
   {
      const int n = NBC > 0 ? NBC : nb;
      const V2 xav = {xa, xa};
      V2 s0 = {0.0, 0.0}, s1 = {0.0, 0.0};
      int b = 0;
      for (; b + 4 <= n; b += 4)
      {
         const V2 r0 = Load2_(r + b), r1 = Load2_(r + b + 2);
         s0 += r0 * Load2_(x + b);
         s1 += r1 * Load2_(x + b + 2);
         Store2_(y + b, Load2_(y + b) + r0 * xav);
         Store2_(y + b + 2, Load2_(y + b + 2) + r1 * xav);
      }
      for (; b + 2 <= n; b += 2)
      {
         const V2 r0 = Load2_(r + b);
         s0 += r0 * Load2_(x + b);
         Store2_(y + b, Load2_(y + b) + r0 * xav);
      }
      real_t s = 0.0;
      for (; b < n; b++)
      {
         s += r[b] * x[b];
         y[b] += r[b] * xa;
      }
      const V2 t = s0 + s1;
      return s + (t[0] + t[1]);
   }

   template <int NBC>
   static real_t RowDot_(int nb, const real_t *r, const real_t *x)
   {
      const int n = NBC > 0 ? NBC : nb;
      V2 s0 = {0.0, 0.0}, s1 = {0.0, 0.0};
      int b = 0;
      for (; b + 4 <= n; b += 4)
      {
         s0 += Load2_(r + b) * Load2_(x + b);
         s1 += Load2_(r + b + 2) * Load2_(x + b + 2);
      }
      for (; b + 2 <= n; b += 2) { s0 += Load2_(r + b) * Load2_(x + b); }
      real_t s = 0.0;
      for (; b < n; b++) { s += r[b] * x[b]; }
      const V2 t = s0 + s1;
      return s + (t[0] + t[1]);
   }

   template <int NM>
   void Dispatch_(int c, const real_t *const *x, real_t *y, bool shared) const
   {
      if (nb_ == 12) { shared ? Shared_<NM, 12>(c, y) : Sweep_<NM, 12>(c, x, y); }
      else if (nb_ == 30) { shared ? Shared_<NM, 30>(c, y) : Sweep_<NM, 30>(c, x, y); }
      else { shared ? Shared_<NM, 0>(c, y) : Sweep_<NM, 0>(c, x, y); }
   }

   /// The element-order sweep of DenseBlockStiffness::Apply_ over the patch
   /// of c: when element e is reached, y_e holds the transposed face terms of
   /// the lower patch elements; e adds its diagonal block and its own faces.
   template <int NM, int NBC>
   void Sweep_(int c, const real_t *const *x, real_t *y) const
   {
      const Cluster &K = C_[c];
      const int nb = NBC > 0 ? NBC : nb_;
      constexpr int L = NBC > 0 ? NBC : kMaxNB;
      const int np = static_cast<int>(K.pe.size());
      for (int i = 0; i < np; i++)
      {
         const int *ie = ElementDofs(K.pe[i]);
         for (int a = 0; a < nb; a++) { y[ie[a]] = 0.0; }
      }
      real_t xe[NM][L], xh[NM][L], acc[L], yh[L];
      const real_t *Dp = K.D.data(), *Bp = K.B.data();
      const std::size_t nbb = static_cast<std::size_t>(nb) * nb * NM;
      for (int i = 0; i < np; i++)
      {
         const int *ie = ElementDofs(K.pe[i]);
         for (int a = 0; a < nb; a++)
         {
            for (int m = 0; m < NM; m++) { xe[m][a] = x[m][ie[a]]; }
            acc[a] = y[ie[a]];
         }
         const real_t *d = Dp + static_cast<std::size_t>(i) * tri_ * NM;
         for (int a = 0; a < nb; a++)
         {
            const int len = nb - a;
            V2 sv = {0.0, 0.0};
            real_t s0 = 0.0;
            for (int m = 0; m < NM; m++)
            {
               const real_t *dr = d + m * len;
               const real_t xa = xe[m][a], *xm = xe[m] + a;
               const V2 xav = {xa, xa};
               real_t *am = acc + a;
               s0 += dr[0] * xa;   // the diagonal entry
               int b = 1;
               for (; b + 1 < len; b += 2)
               {
                  const V2 dv = Load2_(dr + b);
                  Store2_(am + b, Load2_(am + b) + dv * xav);
                  sv += dv * Load2_(xm + b);
               }
               if (b < len)
               {
                  am[b] += dr[b] * xa;
                  s0 += dr[b] * xm[b];
               }
            }
            acc[a] += s0 + (sv[0] + sv[1]);
            d += NM * len;
         }
         for (int f = K.fown[i]; f < K.fown[i + 1]; f++)
         {
            const int *ih = ElementDofs(K.fhi[f]);
            for (int b = 0; b < nb; b++)
            {
               for (int m = 0; m < NM; m++) { xh[m][b] = x[m][ih[b]]; }
               yh[b] = 0.0;
            }
            const real_t *br = Bp + static_cast<std::size_t>(f) * nbb;
            for (int a = 0; a < nb; a++)
            {
               real_t sa = 0.0;
               for (int m = 0; m < NM; m++, br += nb) { sa += RowDotAxpy_<NBC>(nb, br, xh[m], xe[m][a], yh); }
               acc[a] += sa;
            }
            for (int b = 0; b < nb; b++) { y[ih[b]] += yh[b]; }
         }
         for (int a = 0; a < nb; a++) { y[ie[a]] = acc[a]; }
      }
   }

   template <int NM, int NBC>
   void Shared_(int c, real_t *y) const
   {
      const Cluster &K = C_[c];
      const int nb = NBC > 0 ? NBC : nb_;
      constexpr int L = NBC > 0 ? NBC : kMaxNB;
      real_t xh[NM][L];
      const real_t *h = recv_.data(), *Sp = K.S.data();
      const std::size_t nbb = static_cast<std::size_t>(nb) * nb * NM;
      for (std::size_t s = 0; s < K.se.size(); s++)
      {
         const int *ie = ElementDofs(K.se[s]);
         const real_t *hs = h + static_cast<std::size_t>(K.sh[s]) * nb * NM;
         for (int k = 0; k < nb; k++)
         {
            for (int m = 0; m < NM; m++) { xh[m][k] = hs[NM * k + m]; }
         }
         const real_t *br = Sp + s * nbb;
         for (int a = 0; a < nb; a++)
         {
            real_t sa = 0.0;
            for (int m = 0; m < NM; m++, br += nb) { sa += RowDot_<NBC>(nb, br, xh[m]); }
            y[ie[a]] += sa;
         }
      }
   }

   int nm_ = 1, nc_ = 1, ne_ = 0, n_ = 0, nb_ = 0, tri_ = 0, nh_ = 0;
   long long ndir_ = 0;
   real_t asym_ = 0.0, amax_ = 0.0;
   std::vector<int> idx_;
   std::vector<Cluster> C_;
   std::vector<Group> groups_;
   std::vector<int> nbr_slot_;          ///< face-neighbour element -> halo slot
   MPI_Comm comm_ = MPI_COMM_NULL;
   std::vector<real_t> recv_, send_;
   std::vector<MPI_Request> reqs_;
   int nreq_ = 0;
   bool posted_ = false;
   long long rounds_ = 0, recv_values_ = 0;
   double wait_sec_ = 0.0;
};

/// The multirate leapfrog on a NitscheWaveOperator (file comment).
class NitscheLtsRun
{
public:
   /// cl, cl_nbr: clusters of the local and face-neighbour elements; c_fault:
   /// the fault cluster (every fault tet in it; ignored without a fault);
   /// penalty: the operator's penalty_factor.  op.GetTimeStep() must be
   /// dt0 2^c_fault (dt0 with same_step).  Collective.
   NitscheLtsRun(NitscheWaveOperator &op, const std::vector<int> &cl,
                 const std::vector<int> &cl_nbr, int nc, int c_fault, real_t dt0,
                 real_t penalty, bool same_step = false)
      : op_(op), cl_(cl), nc_(nc), cf_(c_fault), dt0_(dt0), same_(same_step)
   {
      ParMesh &pm = op.GetMesh();
      comm_ = pm.GetComm();
      L_ = nc_ - 1;
      MFEM_VERIFY(nc_ >= 1 && nc_ <= 30, "NitscheLtsRun: " << nc_ << " clusters");
      MFEM_VERIFY(!op.HasFault() || (cf_ >= 0 && cf_ < nc_), "NitscheLtsRun: fault cluster "
                  << cf_ << " of " << nc_);
      const real_t dtF = Dt_(op.HasFault() ? cf_ : 0);
      MFEM_VERIFY(!op.HasFault() || std::abs(op.GetTimeStep() - dtF) <= 1e-14 * dtF,
                  "NitscheLtsRun: the operator's step " << op.GetTimeStep()
                  << " s is not the fault cluster's " << dtF << " s; call SetTimeStep(dt_F)");
      MFEM_VERIFY(!op.KelvinVoigtScaleDt(), "NitscheLtsRun: kv_scale = \"dt\" has no single "
                  "dt under local time stepping (design §3.4 step 5); use kv_scale = \"node\"");
      MFEM_VERIFY(!op.KelvinVoigtHasDirichletVelocity(), "NitscheLtsRun: Kelvin-Voigt "
                  "Dirichlet velocity data is not supported");
      // the fault step record is reduced once per sync interval (FaultStep)
      if (op.HasFault()) { op.LtsDeferRecordReduce(true); }
      auto &sample = const_cast<NitscheMaterialSample &>(op.Sample());
      ParFiniteElementSpace &fes = op.Elasticity().GetFESpace();
      Array<int> dir;
      {
         const int nbdr = pm.bdr_attributes.Size() ? pm.bdr_attributes.Max() : 0;
         dir.SetSize(nbdr);
         dir = 0;
         for (int a : op.Boundary().dirichlet_attrs) { dir[a - 1] = 1; }
         has_dirichlet_ = !op.Boundary().dirichlet_attrs.empty();
      }
      blocks_.Build(fes, sample.Lambda(), sample.Mu(), op.KelvinVoigtLambda(),
                    op.KelvinVoigtMu(), dir, penalty, cl, cl_nbr, nc_);
      kv_ = blocks_.NumMatrices() == 2;
      nb_ = blocks_.BlockSize();
      const int ne = pm.GetNE();

      // elements per cluster; Dirichlet elements per cluster
      own_.assign(nc_, {});
      for (int e = 0; e < ne; e++) { own_[cl_[e]].push_back(e); }
      dir_elems_.assign(nc_, {});
      if (has_dirichlet_)
      {
         std::vector<char> isd(ne, 0);
         for (int be = 0; be < pm.GetNBE(); be++)
         {
            const int a = pm.GetBdrAttribute(be);
            if (!dir[a - 1]) { continue; }
            int f, o, e1, e2;
            pm.GetBdrElementFace(be, &f, &o);
            pm.GetFaceElements(f, &e1, &e2);
            if (e1 >= 0 && e2 < 0 && !pm.FaceIsInterior(f)) { isd[e1] = 1; }
         }
         for (int e = 0; e < ne; e++) { if (isd[e]) { dir_elems_[cl_[e]].push_back(e); } }
      }

      // the fault: every fault tet in c_F (design §3.4 step 4)
      if (op.HasFault())
      {
         const auto &el = op.Elasticity();
         int bad = 0;
         for (int f : el.GetFaultInteriorFaces())
         {
            int e1, e2;
            pm.GetFaceElements(f, &e1, &e2);
            if (cl_[e1] != cf_ || cl_[e2] != cf_) { bad++; }
         }
         for (int sf : el.GetFaultSharedFaces())
         {
            int e1, e2;
            pm.GetFaceElements(pm.GetSharedFace(sf), &e1, &e2);
            if (cl_[e1] != cf_ || cl_nbr[-1 - e2] != cf_) { bad++; }
         }
         MPI_Allreduce(MPI_IN_PLACE, &bad, 1, MPI_INT, MPI_SUM, comm_);
         MFEM_VERIFY(bad == 0, "NitscheLtsRun: " << bad << " fault faces have a tet outside "
                     "the fault cluster " << cf_);
      }

      // mass: M^{-1} (the operator's) and, on the absorbing elements,
      // (rho M_e + 1/2 dt_(c_e) C_e)^{-1}
      const NitscheMass &M = op.Mass();
      Vector rho;
      RhoPerElement(sample, rho);
      nd_ = nb_ / 3;
      minv_.resize(static_cast<std::size_t>(ne) * nd_ * nd_);
      inv_rho_.resize(ne);
      for (int e = 0; e < ne; e++)
      {
         const DenseMatrix &Mi = M.ElementMassInverse(e);
         MFEM_VERIFY(Mi.Height() == nd_, "NitscheLtsRun: element " << e << " mass block "
                     << Mi.Height() << " != " << nd_);
         for (int k = 0; k < nd_; k++)
            for (int l = 0; l < nd_; l++) { minv_[(static_cast<std::size_t>(e) * nd_ + k) * nd_ + l] = Mi(k, l); }
         inv_rho_[e] = 1.0 / rho(e);
      }
      absorb_.assign(ne, -1);
      int local_abs = op.HasDamping() ? 1 : 0, any_abs = 0;
      MPI_Allreduce(&local_abs, &any_abs, 1, MPI_INT, MPI_LOR, comm_);
      if (any_abs)
      {
         Vector hdt(ne);
         for (int e = 0; e < ne; e++) { hdt(e) = 0.5 * Dt_(cl_[e]); }
         NitscheMass B(pm, op.GetOrder(), rho, &op.DampingBlocks(), 0.0, &hdt);   // collective
         int q = 0;
         for (const auto &kv : op.DampingBlocks())
         {
            const int e = kv.first;
            const DenseMatrix *Bi = B.ElementBlockInverse(e);
            MFEM_VERIFY(Bi && Bi->Height() == nb_, "NitscheLtsRun: no damped block for element " << e);
            absorb_[e] = q++;
            for (int a = 0; a < nb_; a++)
               for (int b = 0; b < nb_; b++) { binv_.push_back((*Bi)(a, b)); }
            for (int a = 0; a < nb_; a++)
               for (int b = 0; b < nb_; b++) { cblk_.push_back(kv.second(a, b)); }
            abs_elems_[cl_[e]].push_back(e);
         }
      }

      const int n = op.Height();
      a_.assign(nc_, Vector(n));
      for (Vector &a : a_) { a = 0.0; }
      if (has_dirichlet_) { FD_.SetSize(n); }
   }

   int NumClusters() const { return nc_; }
   int FaultCluster() const { return cf_; }
   real_t Dt0() const { return dt0_; }
   /// Step of cluster c.
   real_t ClusterDt(int c) const { return Dt_(c); }
   /// Ticks per sync interval (2^L; 1 with same_step).
   int TicksPerSync() const { return same_ ? 1 : (1 << L_); }
   real_t SyncInterval() const { return TicksPerSync() * dt0_; }
   /// Ticks between fault closures (2^c_F; 1 with same_step).
   int TicksPerFaultStep() const { return same_ ? 1 : (1 << cf_); }
   const NitscheLtsBlocks &Blocks() const { return blocks_; }
   long long Ticks() const { return ticks_; }
   /// Patch elements kicked since construction (start and end kicks), the
   /// work count of the element updates.
   long long PatchUpdates() const { return patch_updates_; }

   /// a_c at (u, v, t) for every cluster with M^{-1} (NewmarkSolver's first
   /// Mult); t starts the first sync interval.
   void Init(const Vector &u, const Vector &v, real_t t)
   {
      t0_ = t;
      m_ = 0;
      blocks_.PostHalo(L_, u, kv_ ? &v : nullptr);
      for (int c = 0; c <= L_; c++) { blocks_.MultInterior(c, u, kv_ ? &v : nullptr, a_[c]); }
      blocks_.WaitHalo();
      for (int c = 0; c <= L_; c++) { FinishForce_(c, -1, v, t); }
   }

   /// Ticks up to and including the next fault closure (every tick without a
   /// fault).  True when the step ended at a sync point (v synchronised).
   bool FaultStep(Vector &u, Vector &v, real_t &t)
   {
      const int need = op_.HasFault() ? cf_ : 0;
      int top;
      do { top = Tick_(u, v, t); } while (top < need);
      if (m_ == 0 && op_.HasFault()) { op_.LtsReduceRecord(); }   // once per sync interval
      return m_ == 0;
   }

private:
   real_t Dt_(int c) const { return same_ ? dt0_ : std::ldexp(dt0_, c); }

   /// One tick; returns the top cluster that ended a step.
   int Tick_(Vector &u, Vector &v, real_t &t)
   {
      for (int c = 0; c <= L_; c++)
      {
         if (same_ || m_ % (1 << c) == 0) { Kick_(c, v); }
      }
      u.Add(dt0_, v);
      t = t0_ + (m_ + 1) * dt0_;
      int top = 0;
      if (same_) { top = L_; }
      else { while (top < L_ && (m_ + 1) % (1 << (top + 1)) == 0) { top++; } }
      blocks_.PostHalo(top, u, kv_ ? &v : nullptr);
      if (op_.HasFault() && cf_ <= top) { op_.LtsAdvanceFault(u, t, kv_ ? &v : nullptr); }
      for (int c = 0; c <= top; c++) { blocks_.MultInterior(c, u, kv_ ? &v : nullptr, a_[c]); }
      blocks_.WaitHalo();
      if (has_dirichlet_) { op_.Elasticity().AssembleDirichletOnlyRHS(FD_, t); }
      for (int c = 0; c <= top; c++) { FinishForce_(c, top, v, t); }
      for (int c = 0; c <= top; c++) { Kick_(c, v); }
      ticks_++;
      if (++m_ == TicksPerSync()) { m_ = 0; t0_ = t; }
      return top;
   }

   /// a_c holds K^(c) u + K_eta^(c) v from the interior blocks: add the
   /// shared blocks, negate, add the loads, apply B_top^{-1} (M^{-1} with
   /// top < 0).  Dirichlet data F_D(t) is in FD_ (Tick_).
   void FinishForce_(int c, int top, const Vector &v, real_t t)
   {
      Vector &a = a_[c];
      blocks_.AddShared(c, a);
      const std::vector<int> &pe = blocks_.Patch(c);
      for (int e : pe)
      {
         const int *ie = blocks_.ElementDofs(e);
         for (int k = 0; k < nb_; k++) { a(ie[k]) = -a(ie[k]); }
      }
      if (op_.HasFault() && c == cf_) { op_.LtsAddFaultLoads(a); }
      if (has_dirichlet_)
      {
         if (top < 0) { op_.Elasticity().AssembleDirichletOnlyRHS(FD_, t); }
         for (int e : dir_elems_[c])
         {
            const int *ie = blocks_.ElementDofs(e);
            for (int k = 0; k < nb_; k++) { a(ie[k]) += FD_(ie[k]); }
         }
      }
      auto it = abs_elems_.find(c);
      if (it != abs_elems_.end())
      {
         for (int e : it->second)
         {
            const int *ie = blocks_.ElementDofs(e);
            const real_t *C = cblk_.data() + static_cast<std::size_t>(absorb_[e]) * nb_ * nb_;
            for (int i = 0; i < nb_; i++)
            {
               real_t s = 0.0;
               for (int j = 0; j < nb_; j++) { s += C[i * nb_ + j] * v(ie[j]); }
               a(ie[i]) -= s;
            }
         }
      }
      ApplyMass_(pe, top, a);
   }

   void ApplyMass_(const std::vector<int> &pe, int top, Vector &a) const
   {
      real_t x[kMaxNB], y[kMaxNB];
      for (int e : pe)
      {
         const int *ie = blocks_.ElementDofs(e);
         for (int k = 0; k < nb_; k++) { x[k] = a(ie[k]); }
         if (top >= 0 && absorb_[e] >= 0 && cl_[e] <= top)
         {
            const real_t *B = binv_.data() + static_cast<std::size_t>(absorb_[e]) * nb_ * nb_;
            for (int i = 0; i < nb_; i++)
            {
               real_t s = 0.0;
               for (int j = 0; j < nb_; j++) { s += B[i * nb_ + j] * x[j]; }
               y[i] = s;
            }
            for (int k = 0; k < nb_; k++) { a(ie[k]) = y[k]; }
            continue;
         }
         const real_t *Mi = minv_.data() + static_cast<std::size_t>(e) * nd_ * nd_;
         for (int comp = 0; comp < 3; comp++)
         {
            const real_t *xc = x + comp * nd_;
            for (int i = 0; i < nd_; i++)
            {
               real_t s = 0.0;
               for (int j = 0; j < nd_; j++) { s += Mi[i * nd_ + j] * xc[j]; }
               a(ie[comp * nd_ + i]) = s * inv_rho_[e];
            }
         }
      }
   }

   void Kick_(int c, Vector &v)
   {
      const real_t h = 0.5 * Dt_(c);
      const Vector &a = a_[c];
      const std::vector<int> &pe = blocks_.Patch(c);
      for (int e : pe)
      {
         const int *ie = blocks_.ElementDofs(e);
         for (int k = 0; k < nb_; k++) { v(ie[k]) += h * a(ie[k]); }
      }
      patch_updates_ += static_cast<long long>(pe.size());
   }

   static constexpr int kMaxNB = 168;

   NitscheWaveOperator &op_;
   std::vector<int> cl_;
   int nc_ = 1, L_ = 0, cf_ = 0, nb_ = 0, nd_ = 0;
   real_t dt0_ = 0.0;
   bool same_ = false, kv_ = false, has_dirichlet_ = false;
   MPI_Comm comm_ = MPI_COMM_NULL;
   NitscheLtsBlocks blocks_;
   std::vector<std::vector<int>> own_, dir_elems_;
   std::map<int, std::vector<int>> abs_elems_;   ///< cluster -> its absorbing elements
   std::vector<int> absorb_;                    ///< element -> absorbing block, or -1
   std::vector<real_t> minv_, inv_rho_, binv_, cblk_;
   std::vector<Vector> a_;
   Vector FD_;
   real_t t0_ = 0.0;
   int m_ = 0;
   long long ticks_ = 0, patch_updates_ = 0;
};

#endif // MFEM_USE_MPI

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_NITSCHE_LTS_RUN_HPP

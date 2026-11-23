// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// fault_side_stations.hpp — per-side fault stations for the Nitsche
// dynamic-rupture driver ([output.stations] per_side = true; TPV6/TPV7/TPV31
// plan rev 1, item A4; Chunhui Zhao, 2026-10-19).
//
// A station on the fault has a displacement and a velocity on each side.  At
// the station's snapped fault node P every local fault face whose closure
// holds P contributes, for each of its elements on this rank, u and du/dt of
// that element evaluated at P (the shape functions at P are kept from the
// setup, so a step costs one small dot product per copy).  The side of each
// element is the fault orientation's (FaultOrientation::PlusElement):
//
//   far side  = the plus side (nhat points from it into the near side),
//   near side = the minus side,
//
// so u(far) - u(near) = delta, the slip of <prefix>_station_<name>.dat
// (<prefix>_station_conventions.txt).  Each side's value is the face-area-
// weighted average over its copies, summed on rank 0 in face-centroid order,
// so it does not depend on the partition.  Components in each copy's face
// frame (the operator's N-09 frame from the oriented normal and up):
//
//   h = u.shat,   v = u.dhat,   n = -u.nhat,
//
// -nhat pointing from the near side into the far side.  For a vertical fault
// with up = +z this is the split-node convention of the SCEC TPV6/TPV7
// description: h-slip = u_h(far) - u_h(near) > 0 for right-lateral slip,
// v positive down-dip, n positive from near to far (normal slip > 0 opening).
//
// Files <prefix>_station_<name>_farside.dat and _nearside.dat, columns
//   t h_disp h_vel h_stress v_disp v_vel v_stress n_disp n_vel n_stress
// (m, m/s, Pa) in the order of the DRDG3D reference files.  The stress
// columns hold the station's traction from the slip station file, the same
// in both files: h_stress = tau2, v_stress = tau1 (SCEC shear tractions) and
// n_stress = -sigma_n (tension positive, as in the SCEC per-side files).
// The header names each side's material (the operator's element centroid
// samples), the check that the weak material is on the side the deck means.

#ifndef MFEM_SEAS_FAULT_SIDE_STATIONS_HPP
#define MFEM_SEAS_FAULT_SIDE_STATIONS_HPP

#include "mfem.hpp"
#include "fault_station_file.hpp"
#include "../domain/nitsche_wave_operator.hpp"
#include "../fault/fault_frame.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{

class FaultSideStationWriter
{
public:
   static constexpr real_t kTie = 1e-6;   // m, as StationCopyAverage
   enum Side { kFar = 0, kNear = 1 };

   /// Collective.  `snap` from SnapFaultStations (the node P of every station
   /// on every rank); `up` the deck's [fault_geometry] up.
   FaultSideStationWriter(NitscheWaveOperator &op, const std::vector<FaultStation> &st,
                          const std::vector<FaultStationSnap> &snap, const real_t up[3],
                          const std::string &dir, const std::string &prefix)
      : st_(st), snap_(snap)
   {
      ns_ = static_cast<int>(st.size());
      MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
      ElasticityDomainOperator<ParMesh> &E = op.Elasticity();
      ParMesh &pm = op.GetMesh();
      const ParFiniteElementSpace &fes = E.GetFESpace();
      const FaultOrientation &orient = op.Orientation();
      const NitscheMaterialSample &sample = op.Sample();
      const int n_int = E.GetFaultInteriorFaces().Size();
      const int nloc = n_int + E.GetFaultSharedFaces().Size();
      MFEM_VERIFY(nloc == orient.NumLocalFaces(), "FaultSideStationWriter: " << nloc
                  << " local fault faces, the orientation has " << orient.NumLocalFaces());
      std::vector<real_t> info;   // per copy: s side cx cy cz area mu rho cp n d s (18)
      for (int fi = 0; fi < nloc; fi++)
      {
         FaceElementTransformations *Tr =
            (fi < n_int) ? pm.GetInteriorFaceTransformations(E.GetFaultInteriorFaces()[fi])
                         : pm.GetSharedFaceTransformations(E.GetFaultSharedFaces()[fi - n_int]);
         MFEM_VERIFY(Tr, "FaultSideStationWriter: no transformation for local fault face " << fi);
         // face corners in the face map's own order, then centroid and area
         // from the canonical (sorted) vertices, as StationCopyAverage
         real_t A[3][3];
         for (int a = 0; a < 3; a++)
         {
            IntegrationPoint ip;
            ip.Set2(a == 1 ? 1.0 : 0.0, a == 2 ? 1.0 : 0.0);
            Vector x(3);
            Tr->Transform(ip, x);
            for (int c = 0; c < 3; c++) { A[a][c] = x(c); }
         }
         long long key[3];
         real_t xv[3][3];
         orient.FaceVertexData(fi, key, xv);
         std::array<std::array<real_t, 3>, 3> v;
         for (int i = 0; i < 3; i++) { for (int c = 0; c < 3; c++) { v[i][c] = xv[i][c]; } }
         std::sort(v.begin(), v.end());
         real_t cen[3], e1[3], e2[3];
         for (int c = 0; c < 3; c++)
         {
            cen[c] = (v[0][c] + v[1][c] + v[2][c]) / 3.0;
            e1[c] = v[1][c] - v[0][c];
            e2[c] = v[2][c] - v[0][c];
         }
         const real_t cx = e1[1] * e2[2] - e1[2] * e2[1], cy = e1[2] * e2[0] - e1[0] * e2[2],
                      cz = e1[0] * e2[1] - e1[1] * e2[0];
         const real_t area = 0.5 * std::sqrt(cx * cx + cy * cy + cz * cz);
         fault_frame::FaceVertices fv;
         fv.nv = 3;
         orient.FaceVertexData(fi, fv.key, fv.x);
         fault_frame::Frame fr;
         MFEM_VERIFY(fault_frame::ComputeFrame(orient.Normal(fi), up,
                                               fault_frame::DegeneracyThreshold(1.0), fv, fr),
                     "FaultSideStationWriter: degenerate frame on local fault face " << fi);
         for (int s = 0; s < ns_; s++)
         {
            real_t xi, eta;
            if (!PointInTriangle_(A, snap[s].node, xi, eta)) { continue; }
            IntegrationPoint fip;
            fip.Set2(xi, eta);
            Tr->SetAllIntPoints(&fip);
            const int plus = orient.PlusElement(fi);
            const int ne = (fi < n_int) ? 2 : 1;
            for (int k = 0; k < ne; k++)
            {
               const int e = (k == 0) ? Tr->Elem1No : Tr->Elem2No;
               const IntegrationPoint &eip =
                  (k == 0) ? Tr->GetElement1IntPoint() : Tr->GetElement2IntPoint();
               Copy c;
               c.s = s;
               c.side = (e == plus) ? kFar : kNear;
               for (int d = 0; d < 3; d++)
               {
                  c.key[d] = cen[d];
                  c.n[d] = fr.n[d];
                  c.d[d] = fr.d[d];
                  c.t[d] = fr.s[d];
               }
               c.area = area;
               const FiniteElement *fe = fes.GetFE(e);
               c.shape.SetSize(fe->GetDof());
               fe->CalcShape(eip, c.shape);
               fes.GetElementVDofs(e, c.vdofs);
               copies_.push_back(c);
               const real_t lam = sample.LambdaGF()(e), mu = sample.MuGF()(e),
                            rho = sample.RhoGF()(e);
               const real_t row[18] = {real_t(s), real_t(c.side), cen[0], cen[1], cen[2], area,
                                       mu, rho, std::sqrt((lam + 2.0 * mu) / rho),
                                       fr.n[0], fr.n[1], fr.n[2], fr.d[0], fr.d[1], fr.d[2],
                                       fr.s[0], fr.s[1], fr.s[2]};
               info.insert(info.end(), row, row + 18);
            }
         }
      }
      std::vector<real_t> all = GatherToRoot_(info);
      std::vector<int> cnt(2 * ns_, 0);
      if (rank_ == 0)
      {
         const int w = 18, nrow = static_cast<int>(all.size()) / w;
         std::vector<int> ord = SortedRows_(all, w, nrow);
         // per station and side: copies, material range, the first copy's frame
         std::vector<std::array<real_t, 6>> mat(2 * ns_, {1e300, -1e300, 1e300, -1e300, 1e300, -1e300});
         std::vector<std::array<real_t, 9>> frame(2 * ns_);
         for (int i : ord)
         {
            const real_t *r = &all[w * i];
            const int k = 2 * static_cast<int>(r[0]) + static_cast<int>(r[1]);
            if (cnt[k] == 0) { for (int c = 0; c < 9; c++) { frame[k][c] = r[9 + c]; } }
            cnt[k]++;
            for (int q = 0; q < 3; q++)
            {
               mat[k][2 * q] = std::min(mat[k][2 * q], r[6 + q]);
               mat[k][2 * q + 1] = std::max(mat[k][2 * q + 1], r[6 + q]);
            }
         }
         files_.resize(2 * ns_);
         for (int s = 0; s < ns_; s++)
         {
            MFEM_VERIFY(cnt[2 * s] > 0 && cnt[2 * s] == cnt[2 * s + 1],
                        "per-side station " << st[s].name << ": " << cnt[2 * s] << " far-side and "
                        << cnt[2 * s + 1] << " near-side copies at the node (" << snap[s].node[0]
                        << ", " << snap[s].node[1] << ", " << snap[s].node[2] << "); every fault "
                        "face through the node must give one element on each side");
            for (int side = 0; side < 2; side++)
            {
               const int k = 2 * s + side;
               const std::string name = dir + "/" + prefix + "_station_" + st[s].name
                                        + (side == kFar ? "_farside.dat" : "_nearside.dat");
               files_[k].open(name);
               MFEM_VERIFY(files_[k].is_open(), "per-side station file " << name
                           << " cannot be written");
               const auto &f = frame[k];
               files_[k] << std::setprecision(12)
                         << "# fault station " << st[s].name << ", "
                         << (side == kFar ? "far" : "near") << " side: requested (" << st[s].x[0]
                         << ", " << st[s].x[1] << ", " << st[s].x[2] << ") m, at the fault node ("
                         << snap[s].node[0] << ", " << snap[s].node[1] << ", " << snap[s].node[2]
                         << ") m\n# " << cnt[k] << " copies (elements of the fault faces through "
                         "the node), face-area-weighted average\n# side material (element "
                         "samples, min max): mu " << mat[k][0] << " " << mat[k][1] << " Pa, rho "
                         << mat[k][2] << " " << mat[k][3] << " kg/m^3, cp " << mat[k][4] << " "
                         << mat[k][5] << " m/s\n"
                         << "# far side = plus side of the fault orientation (nhat points from it "
                         "into the near side); u(far) - u(near) = slip of " << prefix << "_station_"
                         << st[s].name << ".dat\n"
                         << "# h = shat = (" << f[6] << ", " << f[7] << ", " << f[8] << "), v = dhat"
                         " = (" << f[3] << ", " << f[4] << ", " << f[5] << "), n = -nhat = ("
                         << -f[0] << ", " << -f[1] << ", " << -f[2] << ") (near side to far side)\n"
                         << "# Columns: t h_disp h_vel h_stress v_disp v_vel v_stress n_disp n_vel "
                         "n_stress\n# units m, m/s, Pa; h_stress = tau2, v_stress = tau1 (SCEC "
                         "shear tractions), n_stress = -sigma_n (tension positive)\n";
            }
         }
      }
      MPI_Bcast(cnt.data(), 2 * ns_, MPI_INT, 0, MPI_COMM_WORLD);
      ncopies_ = cnt;
   }

   /// Copies of station s on side `side` (all ranks).
   int NumCopies(int s, int side) const { return ncopies_[2 * s + side]; }

   /// Side values of station s from the last WriteStep (rank 0): h v n of
   /// displacement then velocity.
   const std::array<real_t, 6> &LastValues(int s, int side) const { return last_[2 * s + side]; }

   /// Collective.  U, Ud: the displacement and velocity (the operator's
   /// space); dd: the station data the slip station files write.
   void WriteStep(real_t t, const Vector &U, const Vector &Ud, const std::vector<DOFData> &dd)
   {
      constexpr int w = 12;   // s side cx cy cz area h v n (u) h v n (du/dt)
      std::vector<real_t> loc;
      loc.reserve(w * (copies_.size() + ns_));
      for (const Copy &c : copies_)
      {
         const int nd = c.shape.Size();
         real_t u[3] = {0.0, 0.0, 0.0}, ud[3] = {0.0, 0.0, 0.0};
         for (int comp = 0; comp < 3; comp++)
         {
            for (int k = 0; k < nd; k++)
            {
               const int vd = c.vdofs[comp * nd + k];
               u[comp] += c.shape(k) * U(vd);
               ud[comp] += c.shape(k) * Ud(vd);
            }
         }
         auto dot = [](const real_t *a, const real_t *b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; };
         const real_t row[w] = {real_t(c.s), real_t(c.side), c.key[0], c.key[1], c.key[2], c.area,
                                dot(u, c.t), dot(u, c.d), -dot(u, c.n),
                                dot(ud, c.t), dot(ud, c.d), -dot(ud, c.n)};
         loc.insert(loc.end(), row, row + w);
      }
      for (int s = 0; s < ns_; s++)
      {
         if (!snap_[s].owner) { continue; }
         const DOFData &q = dd[snap_[s].dof];
         const real_t row[w] = {real_t(s), 2.0, 0.0, 0.0, 0.0, 0.0,
                                q.tau2_corr, q.tau1_corr, -q.sigma_n_corr, 0.0, 0.0, 0.0};
         loc.insert(loc.end(), row, row + w);
      }
      std::vector<real_t> all = GatherToRoot_(loc);
      if (rank_ != 0) { return; }
      const int nrow = static_cast<int>(all.size()) / w;
      std::vector<int> ord = SortedRows_(all, w, nrow);
      std::vector<std::array<real_t, 6>> sum(2 * ns_, {0, 0, 0, 0, 0, 0});
      std::vector<real_t> wsum(2 * ns_, 0.0);
      std::vector<std::array<real_t, 3>> trac(ns_, {0, 0, 0});
      for (int i : ord)
      {
         const real_t *r = &all[w * i];
         const int s = static_cast<int>(r[0]), side = static_cast<int>(r[1]);
         if (side == 2) { for (int c = 0; c < 3; c++) { trac[s][c] = r[6 + c]; } continue; }
         const int k = 2 * s + side;
         for (int c = 0; c < 6; c++) { sum[k][c] += r[5] * r[6 + c]; }
         wsum[k] += r[5];
      }
      last_.assign(2 * ns_, {0, 0, 0, 0, 0, 0});
      for (int k = 0; k < 2 * ns_; k++)
      {
         for (int c = 0; c < 6; c++) { last_[k][c] = sum[k][c] / wsum[k]; }
         const auto &a = last_[k];
         const auto &tr = trac[k / 2];
         files_[k] << std::scientific << std::setprecision(10) << t << " "
                   << a[0] << " " << a[3] << " " << tr[0] << " "
                   << a[1] << " " << a[4] << " " << tr[1] << " "
                   << a[2] << " " << a[5] << " " << tr[2] << "\n";
         // flushed every write, so a job killed at its wall-time limit keeps its rows
         files_[k].flush();
      }
   }

   void Close()
   {
      for (auto &f : files_) { if (f.is_open()) { f.close(); } }
   }

private:
   struct Copy
   {
      int s = 0, side = 0;
      real_t key[3], area = 0.0, n[3], d[3], t[3];
      Vector shape;
      Array<int> vdofs;
   };

   /// P on the flat triangle A[0..2] (within kTie): its coordinates (xi, eta)
   /// in the face map, P = A0 + xi (A1 - A0) + eta (A2 - A0).
   static bool PointInTriangle_(const real_t A[3][3], const real_t P[3], real_t &xi, real_t &eta)
   {
      real_t a[3], b[3], p[3];
      for (int c = 0; c < 3; c++)
      {
         a[c] = A[1][c] - A[0][c];
         b[c] = A[2][c] - A[0][c];
         p[c] = P[c] - A[0][c];
      }
      const real_t aa = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
      const real_t bb = b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
      const real_t ab = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
      const real_t ap = a[0] * p[0] + a[1] * p[1] + a[2] * p[2];
      const real_t bp = b[0] * p[0] + b[1] * p[1] + b[2] * p[2];
      const real_t det = aa * bb - ab * ab;
      if (det <= 0.0) { return false; }
      xi = (bb * ap - ab * bp) / det;
      eta = (aa * bp - ab * ap) / det;
      real_t r2 = 0.0;
      for (int c = 0; c < 3; c++)
      {
         const real_t q = p[c] - xi * a[c] - eta * b[c];
         r2 += q * q;
      }
      // barycentric tolerance: kTie over the longer edge
      const real_t tb = kTie / std::sqrt(std::max(aa, bb));
      if (std::sqrt(r2) > kTie || xi < -tb || eta < -tb || xi + eta > 1.0 + tb) { return false; }
      xi = std::min(std::max(xi, real_t(0.0)), real_t(1.0));
      eta = std::min(std::max(eta, real_t(0.0)), real_t(1.0) - xi);
      return true;
   }

   std::vector<real_t> GatherToRoot_(const std::vector<real_t> &loc) const
   {
      int nr = 0, n = static_cast<int>(loc.size());
      MPI_Comm_size(MPI_COMM_WORLD, &nr);
      std::vector<int> cnt(nr), dsp(nr, 0);
      MPI_Gather(&n, 1, MPI_INT, cnt.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
      std::vector<real_t> all;
      if (rank_ == 0)
      {
         for (int r = 1; r < nr; r++) { dsp[r] = dsp[r - 1] + cnt[r - 1]; }
         all.resize(dsp[nr - 1] + cnt[nr - 1]);
      }
      MPI_Gatherv(loc.data(), n, MPITypeMap<real_t>::mpi_type, all.data(), cnt.data(),
                  dsp.data(), MPITypeMap<real_t>::mpi_type, 0, MPI_COMM_WORLD);
      return all;
   }

   /// Row order by (station, side, face centroid): partition independent.
   static std::vector<int> SortedRows_(const std::vector<real_t> &all, int w, int nrow)
   {
      std::vector<int> ord(nrow);
      for (int i = 0; i < nrow; i++) { ord[i] = i; }
      std::sort(ord.begin(), ord.end(), [&](int a, int b) {
         return std::lexicographical_compare(&all[w * a], &all[w * a] + 5,
                                             &all[w * b], &all[w * b] + 5);
      });
      return ord;
   }

   std::vector<FaultStation> st_;
   std::vector<FaultStationSnap> snap_;
   int ns_ = 0, rank_ = 0;
   std::vector<Copy> copies_;
   std::vector<int> ncopies_;
   std::vector<std::ofstream> files_;
   std::vector<std::array<real_t, 6>> last_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_FAULT_SIDE_STATIONS_HPP

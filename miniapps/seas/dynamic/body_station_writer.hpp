// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// body_station_writer.hpp — off-fault ("body") stations for the Nitsche
// dynamic-rupture driver ([output.body_stations]; TPV34 plan rev 2, step S3;
// Chunhui Zhao, 2026-11-12).
//
// A body station is a point of the volume (a surface point or a borehole
// point).  The writer records the displacement and the velocity there:
//
//   setup   every tetrahedron holding the point (one inside an element; two
//           or more on a shared face, edge or vertex, because the DG space
//           has one copy of the field per element) gives one "copy": the
//           element's shape functions at the point, kept for the run;
//   step    each copy's u and du/dt at the point (one small dot product),
//           kept in memory on the rank that owns the element;
//   flush   every `flush_every` rows the ranks send their rows to rank 0,
//           which writes, per station, the mean over its copies, summed in
//           element-centroid order, so the files do not depend on the
//           partition (to round-off).
//
// The flush is the only collective (one MPI_Gatherv), so a step costs no
// MPI.  A run that is killed loses the rows since the last flush (at most
// flush_every - 1); Close() flushes and must be called before MPI ends.
//
// Components along three deck axes (unit vectors of the mesh frame):
//
//   h = u . h_axis,   v = u . v_axis,   n = u . n_axis.
//
// Files <prefix>_offfault_<name>.dat, columns (the SCEC off-fault layout)
//   t h_disp h_vel v_disp v_vel n_disp n_vel       (s, m, m/s)
//
// Limits (the writer stops instead of giving a wrong value): straight-sided
// tetrahedra only; a station outside the mesh stops the run, naming it.

#ifndef MFEM_SEAS_BODY_STATION_WRITER_HPP
#define MFEM_SEAS_BODY_STATION_WRITER_HPP

#include "mfem.hpp"
#include "fault_station_file.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{

class BodyStationWriter
{
public:
   /// Barycentric tolerance: a point counts as inside a tetrahedron when its
   /// four barycentric coordinates are >= -kTol.
   static constexpr real_t kTol = 1e-9;

   /// Vp, Vs, rho of local element e, for the file header (optional).
   using MaterialFn = std::function<void(int e, real_t &vp, real_t &vs, real_t &rho)>;

   /// Barycentric coordinates of P in the tetrahedron with corners V[0..3];
   /// false when the tetrahedron is degenerate.
   static bool Barycentric(const real_t V[4][3], const real_t P[3], real_t lam[4])
   {
      real_t a[3], b[3], c[3], p[3];
      for (int d = 0; d < 3; d++)
      {
         a[d] = V[1][d] - V[0][d];
         b[d] = V[2][d] - V[0][d];
         c[d] = V[3][d] - V[0][d];
         p[d] = P[d] - V[0][d];
      }
      auto det3 = [](const real_t *x, const real_t *y, const real_t *z)
      {
         return x[0] * (y[1] * z[2] - y[2] * z[1]) - x[1] * (y[0] * z[2] - y[2] * z[0])
                + x[2] * (y[0] * z[1] - y[1] * z[0]);
      };
      const real_t det = det3(a, b, c);
      if (det == 0.0) { return false; }
      lam[1] = det3(p, b, c) / det;
      lam[2] = det3(a, p, c) / det;
      lam[3] = det3(a, b, p) / det;
      lam[0] = 1.0 - lam[1] - lam[2] - lam[3];
      return true;
   }

   /// Collective.  `fes`: the displacement space (vector, on tetrahedra);
   /// `axes`: h, v, n unit vectors (rows).
   BodyStationWriter(const ParFiniteElementSpace &fes, const std::vector<FaultStation> &st,
                     const real_t axes[3][3], const std::string &dir, const std::string &prefix,
                     int flush_every, const MaterialFn &material = MaterialFn())
      : ns_(static_cast<int>(st.size())), every_(flush_every)
   {
      MFEM_VERIFY(flush_every >= 1, "body stations: flush every " << flush_every
                  << " rows; it must be >= 1");
      ParMesh &pm = *fes.GetParMesh();
      comm_ = pm.GetComm();
      MPI_Comm_rank(comm_, &rank_);
      MPI_Comm_size(comm_, &nranks_);
      MFEM_VERIFY(pm.Dimension() == 3 && fes.GetVDim() == 3,
                  "body stations: a 3-D mesh and a 3-component space are needed");
      MFEM_VERIFY(pm.GetNodes() == nullptr || pm.GetNodalFESpace()->GetMaxElementOrder() == 1,
                  "body stations: the mesh has curved elements; only straight-sided tetrahedra "
                  "are supported");
      for (int a = 0; a < 3; a++) { for (int d = 0; d < 3; d++) { ax_[a][d] = axes[a][d]; } }

      // local copies
      constexpr int w = 8;   // s cx cy cz vp vs rho (local copy index)
      std::vector<real_t> info;
      for (int e = 0; e < pm.GetNE(); e++)
      {
         MFEM_VERIFY(pm.GetElementGeometry(e) == Geometry::TETRAHEDRON,
                     "body stations: element " << e << " is not a tetrahedron");
         Array<int> vert;
         pm.GetElementVertices(e, vert);
         real_t V[4][3], lo[3], hi[3];
         for (int i = 0; i < 4; i++)
         {
            const real_t *x = pm.GetVertex(vert[i]);
            for (int d = 0; d < 3; d++)
            {
               V[i][d] = x[d];
               lo[d] = (i == 0) ? x[d] : std::min(lo[d], x[d]);
               hi[d] = (i == 0) ? x[d] : std::max(hi[d], x[d]);
            }
         }
         const real_t pad = 1e-6 * std::max({hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]});
         for (int s = 0; s < ns_; s++)
         {
            const real_t *P = st[s].x;
            if (P[0] < lo[0] - pad || P[0] > hi[0] + pad || P[1] < lo[1] - pad
                || P[1] > hi[1] + pad || P[2] < lo[2] - pad || P[2] > hi[2] + pad) { continue; }
            real_t lam[4];
            if (!Barycentric(V, P, lam)) { continue; }
            if (std::min({lam[0], lam[1], lam[2], lam[3]}) < -kTol) { continue; }
            // the reference point: MFEM's tetrahedron has vertex 0 at the
            // origin and vertices 1, 2, 3 on the axes
            IntegrationPoint ip;
            ip.Set3(std::max(lam[1], real_t(0.0)), std::max(lam[2], real_t(0.0)),
                    std::max(lam[3], real_t(0.0)));
            Copy c;
            c.s = s;
            const FiniteElement *fe = fes.GetFE(e);
            c.shape.SetSize(fe->GetDof());
            fe->CalcShape(ip, c.shape);
            fes.GetElementVDofs(e, c.vdofs);
            real_t vp = 0.0, vs = 0.0, rho = 0.0;
            if (material) { material(e, vp, vs, rho); }
            const real_t row[w] = {real_t(s),
                                   0.25 * (V[0][0] + V[1][0] + V[2][0] + V[3][0]),
                                   0.25 * (V[0][1] + V[1][1] + V[2][1] + V[3][1]),
                                   0.25 * (V[0][2] + V[1][2] + V[2][2] + V[3][2]),
                                   vp, vs, rho, real_t(copies_.size())};
            info.insert(info.end(), row, row + w);
            copies_.push_back(c);
         }
      }

      // rank 0: every copy of every station, in (station, centroid) order
      std::vector<int> cnt;
      const std::vector<real_t> all = GatherToRoot_(info, cnt);
      std::vector<int> ncop(ns_, 0);
      if (rank_ == 0)
      {
         const int nrow = static_cast<int>(all.size()) / w;
         rank_copies_.assign(nranks_, 0);
         std::vector<int> row_rank(nrow);
         for (int r = 0, i = 0; r < nranks_; r++)
         {
            rank_copies_[r] = cnt[r] / w;
            for (int k = 0; k < rank_copies_[r]; k++) { row_rank[i++] = r; }
         }
         std::vector<int> ord(nrow);
         for (int i = 0; i < nrow; i++) { ord[i] = i; }
         std::sort(ord.begin(), ord.end(), [&](int a, int b) {
            return std::lexicographical_compare(&all[w * a], &all[w * a] + 4,
                                                &all[w * b], &all[w * b] + 4);
         });
         std::vector<std::array<real_t, 6>> mat(ns_, {1e300, -1e300, 1e300, -1e300, 1e300, -1e300});
         order_.resize(nrow);
         for (int j = 0; j < nrow; j++)
         {
            const real_t *r = &all[w * ord[j]];
            const int s = static_cast<int>(r[0]);
            order_[j] = {s, row_rank[ord[j]], static_cast<int>(r[7])};
            ncop[s]++;
            for (int q = 0; q < 3; q++)
            {
               mat[s][2 * q] = std::min(mat[s][2 * q], r[4 + q]);
               mat[s][2 * q + 1] = std::max(mat[s][2 * q + 1], r[4 + q]);
            }
         }
         files_.resize(ns_);
         for (int s = 0; s < ns_; s++)
         {
            MFEM_VERIFY(ncop[s] > 0, "body station " << st[s].name << " (" << st[s].x[0] << ", "
                        << st[s].x[1] << ", " << st[s].x[2] << ") m is in no element of the mesh");
            const std::string name = dir + "/" + prefix + "_offfault_" + st[s].name + ".dat";
            files_[s].open(name);
            MFEM_VERIFY(files_[s].is_open(), "body station file " << name << " cannot be written");
            files_[s] << std::setprecision(12) << "# body station " << st[s].name << " at ("
                      << st[s].x[0] << ", " << st[s].x[1] << ", " << st[s].x[2] << ") m\n# "
                      << ncop[s] << " element(s) hold the point; the values are their mean\n";
            if (material)
            {
               files_[s] << "# material of the element(s) (min max): Vp " << mat[s][0] << " "
                         << mat[s][1] << " m/s, Vs " << mat[s][2] << " " << mat[s][3]
                         << " m/s, rho " << mat[s][4] << " " << mat[s][5] << " kg/m^3\n";
            }
            files_[s] << "# h = u . (" << ax_[0][0] << ", " << ax_[0][1] << ", " << ax_[0][2]
                      << "), v = u . (" << ax_[1][0] << ", " << ax_[1][1] << ", " << ax_[1][2]
                      << "), n = u . (" << ax_[2][0] << ", " << ax_[2][1] << ", " << ax_[2][2]
                      << ")\n# Columns: t h_disp h_vel v_disp v_vel n_disp n_vel\n"
                      "# units s, m, m/s\n";
            files_[s].flush();
         }
      }
      MPI_Bcast(ncop.data(), ns_, MPI_INT, 0, comm_);
      ncopies_ = ncop;
   }

   int NumStations() const { return ns_; }
   /// Elements holding station s (all ranks).
   int NumCopies(int s) const { return ncopies_[s]; }
   /// Local copies on this rank.
   int LocalCopies() const { return static_cast<int>(copies_.size()); }

   /// Not collective, except every `flush_every`-th call (Flush).  U, Ud: the
   /// displacement and the velocity in the space given to the constructor.
   void WriteStep(real_t t, const Vector &U, const Vector &Ud)
   {
      if (rank_ == 0) { times_.push_back(t); }
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
         for (int a = 0; a < 3; a++)
         {
            buf_.push_back(u[0] * ax_[a][0] + u[1] * ax_[a][1] + u[2] * ax_[a][2]);
            buf_.push_back(ud[0] * ax_[a][0] + ud[1] * ax_[a][1] + ud[2] * ax_[a][2]);
         }
      }
      if (++rows_ % every_ == 0) { Flush(); }
   }

   /// Collective: write the rows held in memory.
   void Flush()
   {
      if (closed_) { return; }
      const int nrows = rows_ - rows_written_;
      std::vector<int> cnt;
      const std::vector<real_t> all = GatherToRoot_(buf_, cnt);
      buf_.clear();
      rows_written_ = rows_;
      if (rank_ != 0 || nrows == 0) { return; }
      std::vector<int> off(nranks_, 0);
      for (int r = 1; r < nranks_; r++) { off[r] = off[r - 1] + cnt[r - 1]; }
      for (int r = 0; r < nranks_; r++)
      {
         MFEM_VERIFY(cnt[r] == nrows * rank_copies_[r] * 6, "body stations: rank " << r
                     << " sent " << cnt[r] << " values for " << nrows << " row(s) of "
                     << rank_copies_[r] << " copies");
      }
      for (int i = 0; i < nrows; i++)
      {
         std::vector<std::array<real_t, 6>> sum(ns_, {0, 0, 0, 0, 0, 0});
         for (const Slot &q : order_)
         {
            const real_t *v = &all[off[q.rank] + (i * rank_copies_[q.rank] + q.local) * 6];
            for (int c = 0; c < 6; c++) { sum[q.s][c] += v[c]; }
         }
         const real_t t = times_[i];
         for (int s = 0; s < ns_; s++)
         {
            files_[s] << std::scientific << std::setprecision(10) << t;
            for (int c = 0; c < 6; c++) { files_[s] << " " << sum[s][c] / ncopies_[s]; }
            files_[s] << "\n";
         }
      }
      times_.clear();
      for (auto &f : files_) { f.flush(); }
   }

   /// Collective.
   void Close()
   {
      if (closed_) { return; }
      Flush();
      closed_ = true;
      for (auto &f : files_) { if (f.is_open()) { f.close(); } }
   }

private:
   struct Copy
   {
      int s = 0;
      Vector shape;
      Array<int> vdofs;
   };
   struct Slot { int s, rank, local; };

   std::vector<real_t> GatherToRoot_(const std::vector<real_t> &loc, std::vector<int> &cnt) const
   {
      int n = static_cast<int>(loc.size());
      cnt.assign(nranks_, 0);
      std::vector<int> dsp(nranks_, 0);
      MPI_Gather(&n, 1, MPI_INT, cnt.data(), 1, MPI_INT, 0, comm_);
      std::vector<real_t> all;
      if (rank_ == 0)
      {
         for (int r = 1; r < nranks_; r++) { dsp[r] = dsp[r - 1] + cnt[r - 1]; }
         all.resize(dsp[nranks_ - 1] + cnt[nranks_ - 1]);
      }
      MPI_Gatherv(loc.data(), n, MPITypeMap<real_t>::mpi_type, all.data(), cnt.data(),
                  dsp.data(), MPITypeMap<real_t>::mpi_type, 0, comm_);
      return all;
   }

   int ns_ = 0, every_ = 1, rank_ = 0, nranks_ = 1;
   MPI_Comm comm_ = MPI_COMM_NULL;
   real_t ax_[3][3];
   std::vector<Copy> copies_;
   std::vector<int> ncopies_;
   std::vector<real_t> buf_;          // this rank: rows x local copies x 6
   int rows_ = 0, rows_written_ = 0;
   bool closed_ = false;
   // rank 0
   std::vector<int> rank_copies_;     // copies per rank
   std::vector<Slot> order_;          // every copy, (station, centroid) order
   std::vector<real_t> times_;
   std::vector<std::ofstream> files_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_BODY_STATION_WRITER_HPP

// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/station_copy_average.hpp — the station copy average of the Nitsche
// dynamic-rupture driver (Chunhui Zhao, 2026-10-14; moved here unchanged from
// drivers/nitsche_dyn_driver.cpp, 2026-10-20).  The driver supplies the face
// centroid and area of every owned fault DOF (StationCopyGeometry there).
//
// Members (speed-next O2, Chunhui Zhao 2026-10-20): with sub = true (the
// driver's --station-comm sub, its default) each step involves only the
// ranks that hold a copy or an overwritten DOF of some station, in their own
// communicator; its lowest world rank (the root) sums and writes the
// envelope.  The order of the sum is computed once, at construction, as the
// sort of the static rows (station, face centroid) that the world version
// (sub = false, the code before O2) repeats every step; the rows arrive in
// world-rank order either way, so the sort, the sums and the files are the
// same bit for bit.  The envelope files are StationRowFiles.
//
// Extra values (2D in-plane design rev 3, phase D: T and p of thermal
// pressurization): with nx > 0 each owned fault DOF carries nx more values,
// averaged over the copies with the same weights and written back to the
// overwritten DOFs like the others (not in the envelope); nx = 0 is the
// code before, bit for bit.  On a 2-D mesh a fault face is an edge: the
// driver supplies its midpoint and length as the centroid and area.

#ifndef MFEM_SEAS_STATION_COPY_AVERAGE_HPP
#define MFEM_SEAS_STATION_COPY_AVERAGE_HPP

#include "mfem.hpp"
#include "fault_face_flux.hpp"      // DOFData
#include "fault_station_file.hpp"   // FaultStation, FaultStationDistance

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{

/// Station copy average (header, "Station copies").  Every owned DOF within
/// kTie of the globally nearest station distance (3-D, the snapping metric)
/// gets the face-area-weighted average over all DOFs at the station point P
/// (P = the lexicographically smallest such position, copies within kTie of
/// P); rank 0 (with sub, the members' root) sums in face-centroid order, so
/// the sum does not depend on the partition, and writes the envelope
/// (min/max over copies).  Collective (with sub, over the members: header
/// comment "Members").
class StationCopyAverage
{
public:
   static constexpr real_t kTie = 1e-6;   // m
   /// ca: per owned fault DOF, the centroid (3) and the area of its face,
   /// with the face vertices in a canonical order; X: the owned fault DOF
   /// positions (3 per DOF); the envelope rows reach the files every
   /// `flush_every` rows; nx: extra values per owned DOF (header).
   StationCopyAverage(const Vector &ca, const Vector &X, const std::vector<FaultStation> &st,
                      const std::string &dir, const std::string &prefix, bool sub = true,
                      int flush_every = 1, int nx = 0)
      : ns_(static_cast<int>(st.size())), nx_(nx), nt_(nf + nx), ca_(ca), sub_(sub)
   {
      MFEM_VERIFY(nx >= 0, "station copy average: " << nx << " extra values per DOF");
      MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
      const int no = X.Size() / 3;
      MFEM_VERIFY(ca_.Size() == 4 * no, "station copy average: " << ca_.Size() / 4
                  << " face records for " << no << " owned fault DOFs");
      auto dist = [&](int s, int j) { return FaultStationDistance(st[s], &X(3 * j)); };
      std::vector<real_t> dmin(ns_, std::numeric_limits<real_t>::max());
      for (int s = 0; s < ns_; s++)
      {
         for (int j = 0; j < no; j++) { dmin[s] = std::min(dmin[s], dist(s, j)); }
      }
      MPI_Allreduce(MPI_IN_PLACE, dmin.data(), ns_, MPITypeMap<real_t>::mpi_type, MPI_MIN,
                    MPI_COMM_WORLD);
      cand_.assign(ns_, {});
      copy_.assign(ns_, {});
      std::vector<int> ncop;    // copies per station, all ranks
      std::vector<real_t> pt;   // station points P, all ranks
      for (int s = 0; s < ns_; s++)
      {
         // P: the lexicographically smallest candidate position over all ranks
         real_t P[3] = {std::numeric_limits<real_t>::max(), 0.0, 0.0};
         for (int j = 0; j < no; j++)
         {
            if (dist(s, j) > dmin[s] + kTie) { continue; }
            cand_[s].push_back(j);
            const real_t *x = &X(3 * j);
            if (std::lexicographical_compare(x, x + 3, P, P + 3)) { for (int c = 0; c < 3; c++) { P[c] = x[c]; } }
         }
         int nr = 0;
         MPI_Comm_size(MPI_COMM_WORLD, &nr);
         std::vector<real_t> all(3 * nr);
         MPI_Allgather(P, 3, MPITypeMap<real_t>::mpi_type, all.data(), 3,
                       MPITypeMap<real_t>::mpi_type, MPI_COMM_WORLD);
         for (int r = 0; r < nr; r++)
         {
            if (std::lexicographical_compare(&all[3 * r], &all[3 * r] + 3, P, P + 3))
            {
               for (int c = 0; c < 3; c++) { P[c] = all[3 * r + c]; }
            }
         }
         for (int j = 0; j < no; j++)
         {
            const real_t *x = &X(3 * j);
            if (std::hypot(x[0] - P[0], x[1] - P[1], x[2] - P[2]) <= kTie) { copy_[s].push_back(j); }
         }
         int nc = static_cast<int>(copy_[s].size());
         MPI_Allreduce(MPI_IN_PLACE, &nc, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
         MFEM_VERIFY(nc > 0, "station copy average: no DOF at the point of station " << st[s].name);
         ncop.push_back(nc);
         pt.insert(pt.end(), P, P + 3);
      }
      // the members and their communicator; the root is its rank 0
      bool member = false;
      for (int s = 0; s < ns_; s++) { member = member || !cand_[s].empty() || !copy_[s].empty(); }
      if (sub_)
      {
         MPI_Comm_split(MPI_COMM_WORLD, member ? 0 : MPI_UNDEFINED, rank_, &comm_);
         int r = -1;
         if (comm_ != MPI_COMM_NULL) { MPI_Comm_rank(comm_, &r); }
         root_ = (r == 0);
      }
      else { root_ = (rank_ == 0); }
      if (root_)
      {
         for (int s = 0; s < ns_; s++)
         {
            env_.emplace_back();
            env_.back().Open(dir + "/" + prefix + "_station_" + st[s].name + "_copies.dat",
                             flush_every);
            env_.back().Text() << "# station " << st[s].name << ": " << ncop[s] << " DG copies at ("
                               << std::setprecision(12) << pt[3 * s] << ", " << pt[3 * s + 1]
                               << ", " << pt[3 * s + 2]
                               << "); the station file holds their face-area-weighted average\n"
                               << "# t ncopies sigma_n_min sigma_n_max tau1_min tau1_max tau2_min tau2_max"
                               << " V1_min V1_max V2_min V2_max  (station columns' units and signs)\n";
         }
      }
      if (sub_ && comm_ != MPI_COMM_NULL) { BuildOrder_(); }
   }

   ~StationCopyAverage()
   {
      int fin = 0;
      MPI_Finalized(&fin);
      if (!fin && comm_ != MPI_COMM_NULL) { MPI_Comm_free(&comm_); }
   }
   StationCopyAverage(const StationCopyAverage &) = delete;
   StationCopyAverage &operator=(const StationCopyAverage &) = delete;

   /// Replace the station DOFs' data by the copy average; write the envelope.
   /// `extra`: nx values per owned DOF (header, "Extra values"), averaged and
   /// replaced likewise; required when nx > 0, ignored when nx = 0.
   void Apply(real_t t, std::vector<DOFData> &dd, std::vector<real_t> *extra = nullptr)
   {
      MFEM_VERIFY(nx_ == 0 || (extra && extra->size() >= static_cast<size_t>(nx_) * ca_.Size() / 4),
                  "station copy average: " << nx_ << " extra values per DOF, none given");
      if (sub_) { ApplySub_(t, dd, extra); return; }
      std::vector<real_t> loc;
      for (int s = 0; s < ns_; s++)
      {
         for (int j : copy_[s])
         {
            const DOFData &q = dd[j];
            const real_t row[5 + nf] = {real_t(s), ca_(4 * j), ca_(4 * j + 1), ca_(4 * j + 2),
                                        ca_(4 * j + 3), q.slip1, q.slip2, q.V1, q.V2,
                                        q.tau1_corr, q.tau2_corr, q.sigma_n_corr, q.psi};
            loc.insert(loc.end(), row, row + 5 + nf);
            for (int c = 0; c < nx_; c++) { loc.push_back((*extra)[nx_ * j + c]); }
         }
      }
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
      MPI_Gatherv(loc.data(), n, MPITypeMap<real_t>::mpi_type, all.data(), cnt.data(), dsp.data(),
                  MPITypeMap<real_t>::mpi_type, 0, MPI_COMM_WORLD);
      std::vector<real_t> avg(ns_ * nt_, 0.0);
      if (rank_ == 0)
      {
         const int w = 5 + nt_, nrow = static_cast<int>(all.size()) / w;
         std::vector<int> ord(nrow);
         for (int i = 0; i < nrow; i++) { ord[i] = i; }
         std::sort(ord.begin(), ord.end(), [&](int a, int b) {
            return std::lexicographical_compare(&all[w * a], &all[w * a] + 4, &all[w * b], &all[w * b] + 4);
         });
         for (int s = 0; s < ns_; s++)
         {
            real_t wsum = 0.0, lo[nf], hi[nf];
            for (int f = 0; f < nf; f++) { lo[f] = std::numeric_limits<real_t>::max(); hi[f] = -lo[f]; }
            int nc = 0;
            for (int i : ord)
            {
               const real_t *r = &all[w * i];
               if (static_cast<int>(r[0]) != s) { continue; }
               nc++;
               wsum += r[4];
               for (int f = 0; f < nf; f++)
               {
                  avg[s * nt_ + f] += r[4] * r[5 + f];
                  lo[f] = std::min(lo[f], r[5 + f]);
                  hi[f] = std::max(hi[f], r[5 + f]);
               }
               for (int f = nf; f < nt_; f++) { avg[s * nt_ + f] += r[4] * r[5 + f]; }
            }
            for (int f = 0; f < nt_; f++) { avg[s * nt_ + f] /= wsum; }
            WriteEnvelope_(s, t, nc, lo, hi);
         }
      }
      MPI_Bcast(avg.data(), ns_ * nt_, MPITypeMap<real_t>::mpi_type, 0, MPI_COMM_WORLD);
      Overwrite_(avg, dd, extra);
   }

   void SetFlushEvery(int n) { for (auto &f : env_) { f.SetFlushEvery(n); } }
   void Flush() { for (auto &f : env_) { f.Flush(); } }

   /// The owned DOFs Apply reads or overwrites, sorted.
   std::vector<int> Dofs() const
   {
      std::vector<int> d;
      for (int s = 0; s < ns_; s++)
      {
         d.insert(d.end(), cand_[s].begin(), cand_[s].end());
         d.insert(d.end(), copy_[s].begin(), copy_[s].end());
      }
      std::sort(d.begin(), d.end());
      d.erase(std::unique(d.begin(), d.end()), d.end());
      return d;
   }
   /// Ranks in the per-step communicator (sub; the world size otherwise),
   /// and whether this rank writes the envelope.
   int Members() const
   {
      int n = 0;
      if (!sub_) { MPI_Comm_size(MPI_COMM_WORLD, &n); }
      else if (comm_ != MPI_COMM_NULL) { MPI_Comm_size(comm_, &n); }
      return n;
   }
   bool Root() const { return root_; }

private:
   static constexpr int nf = 8;   // slip1 slip2 V1 V2 tau1 tau2 sigma_n psi

   /// Members: gather the static rows (station, centroid, area) once and, on
   /// the root, the order Apply's world version sorts them into.
   void BuildOrder_()
   {
      std::vector<real_t> loc;
      for (int s = 0; s < ns_; s++)
      {
         for (int j : copy_[s])
         {
            const real_t row[5] = {real_t(s), ca_(4 * j), ca_(4 * j + 1), ca_(4 * j + 2),
                                   ca_(4 * j + 3)};
            loc.insert(loc.end(), row, row + 5);
         }
      }
      int nr = 0, n = static_cast<int>(loc.size()) / 5;
      MPI_Comm_size(comm_, &nr);
      if (root_) { cnt_.assign(nr, 0); dsp_.assign(nr, 0); }
      MPI_Gather(&n, 1, MPI_INT, root_ ? cnt_.data() : nullptr, 1, MPI_INT, 0, comm_);
      std::vector<int> c5, d5;
      std::vector<real_t> all;
      if (root_)
      {
         c5.resize(nr);
         d5.resize(nr);
         int tot = 0;   // rows before rank r
         for (int r = 0; r < nr; r++)
         {
            c5[r] = 5 * cnt_[r];
            d5[r] = 5 * tot;
            dsp_[r] = nt_ * tot;
            tot += cnt_[r];
            cnt_[r] *= nt_;
         }
         all.resize(5 * static_cast<size_t>(tot));
      }
      MPI_Gatherv(loc.data(), 5 * n, MPITypeMap<real_t>::mpi_type, all.data(), c5.data(), d5.data(),
                  MPITypeMap<real_t>::mpi_type, 0, comm_);
      nloc_ = n;
      if (!root_) { return; }
      const int w = 5, nrow = static_cast<int>(all.size()) / w;
      std::vector<int> ord(nrow);
      for (int i = 0; i < nrow; i++) { ord[i] = i; }
      std::sort(ord.begin(), ord.end(), [&](int a, int b) {
         return std::lexicographical_compare(&all[w * a], &all[w * a] + 4, &all[w * b], &all[w * b] + 4);
      });
      rows_.assign(ns_, {});
      area_.resize(nrow);
      for (int i = 0; i < nrow; i++) { area_[i] = all[w * i + 4]; }
      for (int s = 0; s < ns_; s++)
      {
         for (int i : ord) { if (static_cast<int>(all[w * i]) == s) { rows_[s].push_back(i); } }
      }
      vals_.resize(static_cast<size_t>(nt_) * nrow);
   }

   /// Apply with the members only (header comment, "Members").
   void ApplySub_(real_t t, std::vector<DOFData> &dd, std::vector<real_t> *extra)
   {
      if (comm_ == MPI_COMM_NULL) { return; }
      std::vector<real_t> loc;
      loc.reserve(static_cast<size_t>(nt_) * nloc_);
      for (int s = 0; s < ns_; s++)
      {
         for (int j : copy_[s])
         {
            const DOFData &q = dd[j];
            const real_t row[nf] = {q.slip1, q.slip2, q.V1, q.V2,
                                    q.tau1_corr, q.tau2_corr, q.sigma_n_corr, q.psi};
            loc.insert(loc.end(), row, row + nf);
            for (int c = 0; c < nx_; c++) { loc.push_back((*extra)[nx_ * j + c]); }
         }
      }
      MPI_Gatherv(loc.data(), static_cast<int>(loc.size()), MPITypeMap<real_t>::mpi_type,
                  vals_.data(), cnt_.data(), dsp_.data(), MPITypeMap<real_t>::mpi_type, 0, comm_);
      std::vector<real_t> avg(ns_ * nt_, 0.0);
      if (root_)
      {
         for (int s = 0; s < ns_; s++)
         {
            real_t wsum = 0.0, lo[nf], hi[nf];
            for (int f = 0; f < nf; f++) { lo[f] = std::numeric_limits<real_t>::max(); hi[f] = -lo[f]; }
            int nc = 0;
            for (int i : rows_[s])
            {
               const real_t ar = area_[i], *v = &vals_[nt_ * static_cast<size_t>(i)];
               nc++;
               wsum += ar;
               for (int f = 0; f < nf; f++)
               {
                  avg[s * nt_ + f] += ar * v[f];
                  lo[f] = std::min(lo[f], v[f]);
                  hi[f] = std::max(hi[f], v[f]);
               }
               for (int f = nf; f < nt_; f++) { avg[s * nt_ + f] += ar * v[f]; }
            }
            for (int f = 0; f < nt_; f++) { avg[s * nt_ + f] /= wsum; }
            WriteEnvelope_(s, t, nc, lo, hi);
         }
      }
      MPI_Bcast(avg.data(), ns_ * nt_, MPITypeMap<real_t>::mpi_type, 0, comm_);
      Overwrite_(avg, dd, extra);
   }

   void WriteEnvelope_(int s, real_t t, int nc, const real_t *lo, const real_t *hi)
   {
      env_[s].Text() << std::scientific << std::setprecision(10) << t << " " << nc << " "
                     << lo[6] << " " << hi[6] << " " << lo[4] << " " << hi[4] << " " << lo[5] << " "
                     << hi[5] << " " << lo[2] << " " << hi[2] << " " << lo[3] << " " << hi[3] << "\n";
      env_[s].EndRow();
   }

   void Overwrite_(const std::vector<real_t> &avg, std::vector<DOFData> &dd,
                   std::vector<real_t> *extra) const
   {
      for (int s = 0; s < ns_; s++)
      {
         const real_t *a = &avg[s * nt_];
         for (int j : cand_[s])
         {
            DOFData &q = dd[j];
            q.slip1 = a[0]; q.slip2 = a[1]; q.V1 = a[2]; q.V2 = a[3];
            q.tau1_corr = a[4]; q.tau2_corr = a[5]; q.sigma_n_corr = a[6]; q.psi = a[7];
            q.slip_rate = std::hypot(a[2], a[3]);
            for (int c = 0; c < nx_; c++) { (*extra)[nx_ * j + c] = a[nf + c]; }
         }
      }
   }

   int ns_, nx_, nt_, rank_ = 0;   // stations, extra values per DOF, nf + nx
   Vector ca_;                                   // owned: centroid (3) + area
   bool sub_ = true, root_ = false;
   MPI_Comm comm_ = MPI_COMM_NULL;               // sub: the members
   std::vector<std::vector<int>> cand_, copy_;   // owned DOFs overwritten / averaged
   std::vector<StationRowFile> env_;             // the root
   // sub: copies on this rank; the root's gather counts and offsets (values),
   // rows per station in summation order, areas and the gathered values
   int nloc_ = 0;
   std::vector<int> cnt_, dsp_;
   std::vector<std::vector<int>> rows_;
   std::vector<real_t> area_, vals_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_STATION_COPY_AVERAGE_HPP

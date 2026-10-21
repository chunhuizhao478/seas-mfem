// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// dynamic/station_copy_average.hpp — the station copy average of the Nitsche
// dynamic-rupture driver (Chunhui Zhao, 2026-10-14; moved here unchanged from
// drivers/nitsche_dyn_driver.cpp, 2026-10-20).  The driver supplies the face
// centroid and area of every owned fault DOF (StationCopyGeometry there).

#ifndef MFEM_SEAS_STATION_COPY_AVERAGE_HPP
#define MFEM_SEAS_STATION_COPY_AVERAGE_HPP

#include "mfem.hpp"
#include "fault_face_flux.hpp"      // DOFData
#include "fault_station_file.hpp"   // FaultStation, FaultStationDistance

#include <algorithm>
#include <cmath>
#include <fstream>
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
/// P); rank 0 sums in face-centroid order, so the sum does not depend on the
/// partition, and writes the envelope (min/max over copies).  Collective.
class StationCopyAverage
{
public:
   static constexpr real_t kTie = 1e-6;   // m
   /// ca: per owned fault DOF, the centroid (3) and the area of its face,
   /// with the face vertices in a canonical order; X: the owned fault DOF
   /// positions (3 per DOF).
   StationCopyAverage(const Vector &ca, const Vector &X, const std::vector<FaultStation> &st,
                      const std::string &dir, const std::string &prefix)
      : ns_(static_cast<int>(st.size())), ca_(ca)
   {
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
         if (rank_ == 0)
         {
            env_.emplace_back(dir + "/" + prefix + "_station_" + st[s].name + "_copies.dat");
            env_.back() << "# station " << st[s].name << ": " << nc << " DG copies at ("
                        << std::setprecision(12) << P[0] << ", " << P[1] << ", " << P[2]
                        << "); the station file holds their face-area-weighted average\n"
                        << "# t ncopies sigma_n_min sigma_n_max tau1_min tau1_max tau2_min tau2_max"
                        << " V1_min V1_max V2_min V2_max  (station columns' units and signs)\n";
         }
      }
   }

   /// Replace the station DOFs' data by the copy average; write the envelope.
   void Apply(real_t t, std::vector<DOFData> &dd)
   {
      constexpr int nf = 8;   // slip1 slip2 V1 V2 tau1 tau2 sigma_n psi
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
      std::vector<real_t> avg(ns_ * nf, 0.0);
      if (rank_ == 0)
      {
         const int w = 5 + nf, nrow = static_cast<int>(all.size()) / w;
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
                  avg[s * nf + f] += r[4] * r[5 + f];
                  lo[f] = std::min(lo[f], r[5 + f]);
                  hi[f] = std::max(hi[f], r[5 + f]);
               }
            }
            for (int f = 0; f < nf; f++) { avg[s * nf + f] /= wsum; }
            env_[s] << std::scientific << std::setprecision(10) << t << " " << nc << " "
                    << lo[6] << " " << hi[6] << " " << lo[4] << " " << hi[4] << " " << lo[5] << " "
                    << hi[5] << " " << lo[2] << " " << hi[2] << " " << lo[3] << " " << hi[3] << "\n";
            env_[s].flush();
         }
      }
      MPI_Bcast(avg.data(), ns_ * nf, MPITypeMap<real_t>::mpi_type, 0, MPI_COMM_WORLD);
      for (int s = 0; s < ns_; s++)
      {
         const real_t *a = &avg[s * nf];
         for (int j : cand_[s])
         {
            DOFData &q = dd[j];
            q.slip1 = a[0]; q.slip2 = a[1]; q.V1 = a[2]; q.V2 = a[3];
            q.tau1_corr = a[4]; q.tau2_corr = a[5]; q.sigma_n_corr = a[6]; q.psi = a[7];
            q.slip_rate = std::hypot(a[2], a[3]);
         }
      }
   }

private:
   int ns_, rank_ = 0;
   Vector ca_;                                   // owned: centroid (3) + area
   std::vector<std::vector<int>> cand_, copy_;   // owned DOFs overwritten / averaged
   std::vector<std::ofstream> env_;              // rank 0
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_STATION_COPY_AVERAGE_HPP

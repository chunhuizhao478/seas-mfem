// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// cycle/pml_partition.hpp — a mesh cut that shares the PML band's work.
//
// Why (Chunhui Zhao, 2026-11-10: "we should look at code, not threshold"):
// in the dynamic regime an element inside a [cycle.pml] band costs several
// times an element outside (the memory terms of domain/pml_memory.hpp).
// MFEM's default cut gives every rank the same number of elements, so the
// ranks that hold the band are slow and the others wait for them at every
// step (San Jacinto alone, 8 ranks, commit 494e2092: the band update was
// about 50 % of the step on a band rank and the fault rank waited about
// 30 % of the time).
//
// What: [cycle.pml] partition_weight = w > 0 cuts the serial mesh with METIS
// vertex weights 1 on an element outside the bands and 1 + w on an element
// whose centroid is inside one (the weighted cut the LTS path uses,
// dynamic/lts_partition.hpp, with one cluster and no fault lock: fault faces
// may be cut, as with the default cut).  w = 0 (default) returns an empty
// vector and the caller keeps MFEM's default cut.
//
// In:  the serial mesh (the same on every rank), the configuration, the
//      communicator.  Out: one rank per element (empty = default cut) and a
//      one-line report.  Collective (rank 0 cuts, the others receive).
//
// A restart needs the cut of the run it continues (the checkpoints and the
// snapshot files are per rank): keep the key as it was.

#ifndef MFEM_SEAS_CYCLE_PML_PARTITION_HPP
#define MFEM_SEAS_CYCLE_PML_PARTITION_HPP

#include "mfem.hpp"
#include "cycle_config.hpp"
#include "../domain/pml_band.hpp"
#include "../dynamic/lts_partition.hpp"
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{
namespace cycle
{

struct PmlPartition
{
   std::vector<int> part;   ///< rank of each serial element; empty = default cut
   std::string report;      ///< for the driver's start-up lines (rank 0)
};

inline PmlPartition PlanPmlPartition(Mesh &smesh, const CycleConfig &cfg, MPI_Comm comm)
{
   PmlPartition R;
   const PmlConfig &P = cfg.pml;
   int rank = 0, nranks = 1;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nranks);
   if (!P.On() || P.partition_weight <= 0.0 || nranks == 1) { return R; }
   const int ne = smesh.GetNE(), dim = smesh.Dimension();
   R.part.assign(ne, 0);
   if (rank == 0)
   {
      // the band of the serial mesh: only the walls and the thickness matter
      // for InBand (the profile numbers are the deck's, as QDPhase sets them)
      PmlBandParams prm;
      prm.thickness = P.thickness;
      prm.power = P.power;
      prm.R0 = P.R0;
      prm.f_alpha = P.Sponge() ? 1.0 : P.f_alpha;
      prm.c_ref = P.c_ref;
      const PmlBand band(smesh, P.wall_attrs, prm);
      // integer METIS weights: 100 per element, 100 (1 + w) inside a band
      std::vector<double> cost(ne, 100.0);
      std::vector<int> cluster(ne, 0);
      long inside = 0;
      Vector c(dim);
      for (int e = 0; e < ne; e++)
      {
         smesh.GetElementCenter(e, c);
         real_t x[3] = {0.0, 0.0, 0.0};
         for (int i = 0; i < dim; i++) { x[i] = c(i); }
         if (band.InBand(x))
         {
            cost[e] = 100.0 * (1.0 + P.partition_weight);
            inside++;
         }
      }
      const Table &e2e = smesh.ElementToElementTable();
      std::vector<int> xadj(e2e.GetI(), e2e.GetI() + ne + 1);
      std::vector<int> adj(e2e.GetJ(), e2e.GetJ() + e2e.Size_of_connections());
      const std::vector<std::pair<int, int>> no_pairs;
      LtsPartitionOptions popt;
      try
      {
         BuildLtsAwarePartition(ne, xadj, adj, cluster, 1, nranks, &cost, no_pairs, popt,
                                R.part);
      }
      catch (const std::exception &ex)
      {
         MFEM_ABORT("[cycle.pml] partition_weight: the weighted METIS cut failed on "
                    << nranks << " ranks: " << ex.what());
      }
      // what the cut achieved: elements and weighted work per rank
      std::vector<long> cnt(nranks, 0), bnd(nranks, 0);
      std::vector<double> work(nranks, 0.0);
      double total = 0.0;
      for (int e = 0; e < ne; e++)
      {
         const int r = R.part[e];
         MFEM_VERIFY(r >= 0 && r < nranks, "[cycle.pml] partition_weight: rank " << r
                     << " of element " << e << " out of range");
         cnt[r]++;
         if (cost[e] > 100.0) { bnd[r]++; }
         work[r] += cost[e];
         total += cost[e];
      }
      long cmin = ne, cmax = 0, bmax = 0;
      double wmax = 0.0;
      for (int r = 0; r < nranks; r++)
      {
         MFEM_VERIFY(cnt[r] > 0, "[cycle.pml] partition_weight: rank " << r << " got no "
                     "element (fewer ranks or a smaller weight)");
         cmin = std::min(cmin, cnt[r]);
         cmax = std::max(cmax, cnt[r]);
         bmax = std::max(bmax, bnd[r]);
         wmax = std::max(wmax, work[r]);
      }
      std::ostringstream s;
      s << std::setprecision(4) << "[cycle.pml] partition_weight " << P.partition_weight
        << ": weighted METIS cut, " << inside << " of " << ne << " elements inside a band; "
        << "elements per rank " << cmin << " to " << cmax << ", most band elements on a rank "
        << bmax << ", weighted work max/mean " << wmax / (total / nranks);
      R.report = s.str();
   }
   MPI_Bcast(R.part.data(), ne, MPI_INT, 0, comm);
   return R;
}

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif

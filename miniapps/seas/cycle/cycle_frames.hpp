// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// cycle_frames.hpp — 2D in-plane design rev 3, phase A: the fully dynamic
// cycle's ONE fault orientation system.  The QD regime's operator builds the
// FaultFrameField with this configuration; the dynamic regime's Nitsche
// operator checks that its own fault setup asks for the same one and reads
// that field (fd_cycle_driver hands it over), so both regimes orient every
// fault face, frame and prescribed-slip face the same way.  FaultVertexRows
// (vertex-normal plan rev1, P3) tabulates the field's fault vertices for the
// output's _vertices.txt.
//
// Scope (Chunhui Zhao, 2026-10-24).

#ifndef MFEM_SEAS_CYCLE_FRAMES_HPP
#define MFEM_SEAS_CYCLE_FRAMES_HPP

#include "mfem.hpp"
#include "cycle_config.hpp"
#include "../fault/fault_frame_field.hpp"

#include <cmath>
#include <map>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{
namespace cycle
{

/// [cycle.dynamic] ref_normal, up, pin_min_margin (the deck keys the
/// dynamic regime read before phase A), with the reference pin and the
/// fault attribute.
inline FaultFrameConfig CycleFrameConfig(const CycleConfig &cfg)
{
   const DynamicConfig &dc = cfg.dyn;
   FaultFrameConfig f;
   FaultOrientationConfig &oc = f.orientation;
   for (int i = 0; i < 3; i++)
   {
      oc.ref_normal[i] = dc.ref_normal[i];
      oc.up[i] = dc.up[i];
      f.up[i] = dc.up[i];
   }
   oc.pin = FaultOrientationConfig::Pin::Reference;
   if (dc.pin_min_margin > 0.0) { oc.c_min = dc.pin_min_margin; }
   oc.fault_attrs = {cfg.fault_attr};
   f.averaged_vertex_normals = dc.averaged_vertex_normals;
   f.vertex_normal_max_turn_deg = dc.vertex_normal_max_turn_deg;
   return f;
}

/// Columns of FaultVertexRows (the output's _vertices.txt).
inline std::vector<std::string> FaultVertexColumns()
{
   return {"key", "x(m)", "y(m)", "kind", "turn(deg)", "nv_x", "nv_y"};
}

/// Vertex-normal plan rev1, P3 (Chunhui Zhao, 2026-11-01): one row per
/// vertex of the 2-D fault, for checking the averaged vertex normals: the
/// census vertex key, the position, the kind (0 tip, 1 smooth, 2 kink; a
/// straight line is smooth), the turn angle between the two faces [deg]
/// (0 at a tip) and the vertex normal n_v, with (0, 0) where each face keeps
/// its own normal (a tip or a kink).  Collective: on the root the rows of
/// every rank, one per key in increasing key order (a vertex shared by
/// faces or ranks is written once), so the table does not depend on the
/// partition; empty on the other ranks.  The field must be 2-D with the
/// option on.
inline void FaultVertexRows(const FaultFrameField &field, MPI_Comm comm, Vector &rows)
{
   MFEM_VERIFY(field.AveragedVertexNormals(), "FaultVertexRows: vertex_normals = \"face\"");
   const FaultOrientation &census = field.FaultCensus();
   MFEM_VERIFY(census.FaceVertices() == 2, "FaultVertexRows: 2-D faults only");
   const int ncol = static_cast<int>(FaultVertexColumns().size());
   std::vector<real_t> loc;
   for (int fi = 0; fi < field.NumFaultFaces(); fi++)
   {
      long long key[3];
      real_t x[3][3];
      census.FaceVertexData(fi, key, x);
      for (int a = 0; a < 2; a++)
      {
         const int kind = field.VertexKind(fi, a);
         const real_t *xv = field.VertexX(fi, a);
         const real_t *nv = field.VertexNormal(fi, a);
         const bool own = (kind == 1);
         loc.insert(loc.end(), {static_cast<real_t>(key[a]), xv[0], xv[1],
                                static_cast<real_t>(kind), field.VertexTurnDeg(fi, a),
                                own ? nv[0] : 0.0, own ? nv[1] : 0.0});
      }
   }
   int rank = 0, nr = 1;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nr);
   int nloc = static_cast<int>(loc.size());
   std::vector<int> cnt(nr), dsp(nr, 0);
   MPI_Gather(&nloc, 1, MPI_INT, cnt.data(), 1, MPI_INT, 0, comm);
   std::vector<real_t> all;
   if (rank == 0)
   {
      for (int r = 1; r < nr; r++) { dsp[r] = dsp[r - 1] + cnt[r - 1]; }
      all.resize(static_cast<size_t>(dsp[nr - 1] + cnt[nr - 1]));
   }
   MPI_Gatherv(loc.data(), nloc, MPITypeMap<real_t>::mpi_type, all.data(), cnt.data(),
               dsp.data(), MPITypeMap<real_t>::mpi_type, 0, comm);
   rows.SetSize(0);
   if (rank != 0) { return; }
   std::map<long long, std::vector<real_t>> by_key;
   for (size_t i = 0; i + ncol <= all.size(); i += ncol)
   {
      const long long k = std::llround(all[i]);
      std::vector<real_t> r(all.begin() + i, all.begin() + i + ncol);
      auto it = by_key.find(k);
      if (it == by_key.end()) { by_key.emplace(k, r); continue; }
      // a vertex seen from two faces or ranks: the field gave it one value
      MFEM_VERIFY(it->second == r, "FaultVertexRows: vertex " << k << " has two "
                  "different rows");
   }
   rows.SetSize(static_cast<int>(by_key.size()) * ncol);
   int o = 0;
   for (const auto &kv : by_key)
   {
      for (real_t v : kv.second) { rows(o++) = v; }
   }
}

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_FRAMES_HPP

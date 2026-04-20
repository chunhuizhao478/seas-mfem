// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC.
// SPDX-License-Identifier: BSD-3-Clause
//
// Shared-fault face-vertex key — canonical integer key for matching a
// triangular fault face across MPI ranks.
//
// Usage today:
//   Only TPV102 (dynamic-rupture) includes this header (via
//   `wave_operator.hpp`).  The R-101 shared-fault verifier pairs
//   cross-rank QP entries by the integer key defined here, replacing the
//   earlier physical-centroid-plus-tolerance pairing that proved fragile
//   at 50-rank production scale (v8.0.0, job 7666233).
//
// Usage plan (BP5 merge):
//   BP5's `ElasticityOperator` currently owns an equivalent, bit-identical
//   `FaceVertexKey` nested struct and `MakeFaceKey` member function in
//   `miniapps/seas/domain/elasticity_operator.hpp` (circa lines 483-514).
//   The long-term plan is for BP5 to ALSO include this header and drop
//   its nested copies — the shared file becomes the single source of
//   truth across quasi-dynamic (BP5) and dynamic (TPV102) paths, eliminating
//   the duplication we are carrying now.
//
//   That BP5-side change is deliberately deferred so the current R-101
//   fix does not touch BP5 source (BP5 is mature and in debug-stable
//   state; pairing the TPV102 fix with a BP5 refactor would mix scopes).
//
// **Invariant until unification:** any change to `FaceVertexKey` or
// `MakeFaceKey` here MUST be mirrored in `elasticity_operator.hpp`
// (and vice versa).  Keeping them bit-identical ensures the merge is a
// simple delete-the-nested-copy operation.

#ifndef MFEM_SEAS_DYNAMIC_SHARED_FAULT_KEY_HPP
#define MFEM_SEAS_DYNAMIC_SHARED_FAULT_KEY_HPP

#include "mfem.hpp"
#include <algorithm>

namespace mfem
{
namespace seas
{
namespace dynamic
{

// Canonical vertex-id triple identifying a triangular fault face across
// MPI ranks.  Keep this struct bit-identical to
// `ElasticityOperator::FaceVertexKey`.
struct FaceVertexKey
{
   HYPRE_BigInt v[3] = {-1, -1, -1};

   bool operator<(const FaceVertexKey &o) const
   {
      if (v[0] != o.v[0]) { return v[0] < o.v[0]; }
      if (v[1] != o.v[1]) { return v[1] < o.v[1]; }
      return v[2] < o.v[2];
   }

   bool operator==(const FaceVertexKey &o) const
   {
      return v[0] == o.v[0] && v[1] == o.v[1] && v[2] == o.v[2];
   }
};

// Build a canonical face key from a local face index.  `gvert` is the
// rank-local table produced by `ParMesh::GetGlobalVertexIndices`.
//
// Keep the body bit-identical to `ElasticityOperator::MakeFaceKey`
// (`miniapps/seas/domain/elasticity_operator.hpp`).
inline FaceVertexKey MakeFaceKey(int local_face,
                                 const Array<HYPRE_BigInt> &gvert,
                                 const Mesh &mesh)
{
   Array<int> verts;
   mesh.GetFaceVertices(local_face, verts);
   FaceVertexKey key;
   for (int j = 0; j < 3 && j < verts.Size(); j++)
   {
      key.v[j] = gvert[verts[j]];
   }
   if (key.v[0] > key.v[1]) { std::swap(key.v[0], key.v[1]); }
   if (key.v[1] > key.v[2]) { std::swap(key.v[1], key.v[2]); }
   if (key.v[0] > key.v[1]) { std::swap(key.v[0], key.v[1]); }
   return key;
}

} // namespace dynamic
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_DYNAMIC_SHARED_FAULT_KEY_HPP

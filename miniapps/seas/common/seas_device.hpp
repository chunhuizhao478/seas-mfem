// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_SEAS_DEVICE_HPP
#define MFEM_SEAS_DEVICE_HPP

#include "mfem.hpp"
#include <string>

namespace mfem
{
namespace seas
{

/// The GPU a rank uses for `--device STR --device-id N` (GPU plans G1, D1):
/// N when N >= 0, otherwise the node-local rank modulo the number of visible
/// GPUs, so several ranks may share a GPU.  The GPUs are counted only for a
/// cuda or hip device string: a GPU build may run `--device cpu` on a node
/// without a GPU, where the count aborts.  Collective on `comm` when it
/// counts.  Make the mfem::Device from it right after MPI starts, before any
/// mesh, space or vector (and before PETSc):
///
///   Device device(device_str, SelectDeviceId(device_str, id_override, comm));
inline int SelectDeviceId(const std::string &device_str, int device_id_override,
                          MPI_Comm comm)
{
   if (device_id_override >= 0) { return device_id_override; }
   int device_id = 0;
#if defined(MFEM_USE_CUDA) || defined(MFEM_USE_HIP)
   // One GPU per rank on a node: node-local rank modulo the visible GPUs.
   const bool gpu_backend = device_str.find("cuda") != std::string::npos ||
                            device_str.find("hip") != std::string::npos;
   const int num_gpu = gpu_backend ? Device::GetDeviceCount() : 0;
   if (num_gpu > 1)
   {
      MPI_Comm node_comm;
      MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node_comm);
      int node_rank = 0;
      MPI_Comm_rank(node_comm, &node_rank);
      MPI_Comm_free(&node_comm);
      device_id = node_rank % num_gpu;
   }
#else
   (void)device_str;
   (void)comm;
#endif
   return device_id;
}

/// Plan D1: host code that is handed state that may live on the device
/// (U, Udot, the acceleration) works on non-owning host views of it, so
/// nothing inside reaches a stale copy.  A read view syncs v to the host
/// (v.HostRead()) and stays valid until v is next written on the device.
/// With the cpu device both are plain aliases: the same bits.
///
/// A view is for host code only.  An MFEM kernel that reads it on the device
/// (Vector::Read(), e.g. in SparseMatrix::Mult) registers v's pointer a
/// second time, and the view's destructor then frees v's device copy; the
/// debug device turns such a read into an abort.
inline Vector HostReadView(const Vector &v)
{
   return Vector(const_cast<real_t *>(v.HostRead()), v.Size());
}

/// A write view: v.HostWrite(), so the device copy is stale until v's next
/// device Read, which copies the host values over.  The view must not be
/// resized; size v first.
inline Vector HostWriteView(Vector &v)
{
   return Vector(v.HostWrite(), v.Size());
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_DEVICE_HPP

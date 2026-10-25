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

// qd_green_operator.hpp — QD speed item 3: the discrete Green's function of
// a quasi-dynamic stage (opt-in, --qd-green).  With the stage operators of
// qd_stage_operators.hpp, one stage maps slip δ to traction linearly:
//
//   [T; T_n] = G·δ + t·T_D,   G   = [C_u; N_u]·K⁻¹·B + [C_δ; N_δ]
//                             T_D = [C_u; N_u]·K⁻¹·D·ĝ,  g(t) = t·ĝ
//
// G is dense, (2 or 3)·N_f rows × 2·N_f columns.  ElasticityDomainOperator
// builds it column by column with blocked solves (BuildQDGreen in
// elasticity_operator_qd.inl); this file holds what it is built into and
// the stage product.  Valid only while K is constant (linear elastic bulk,
// fixed mesh) and the loading has the one time shape g(t) = t·ĝ.
//
// Row layout: global fault index g in [0, N_f), rows ncomp·g + c with
// c = 0 dip, 1 strike, 2 normal (with elastic sigma_n).  Rows are split
// evenly over ALL ranks (not only the ranks that own fault DOFs), so the
// product does not inherit the fault ranks' load imbalance.  Columns are
// the global slip index 2·g + c.  A stage:
//   1. all-gather the owned slip into δ_all (2·N_f values),
//   2. y = G_rows·δ_all + t·T_D,rows on this rank's rows,
//   3. send each value to the rank that owns its fault DOF.
// The owned layouts are those of the stage operators: shear 2·o + c,
// normal o, for owned fault DOF o (global index fault_first + o).
//
// The global index g follows the partition (rank by rank, owned order), so
// it changes with the number of ranks.  The cache file therefore also
// stores a position key per fault DOF (its node, then its face's centroid,
// 6 coordinates), and a read matches every owned DOF to a file row by its
// key: one file serves any number of ranks.
//
// GPU plan G3c: with a device backend, UseDevice() moves G and T_D to the
// device once (after the build or the read; on a CUDA device the host copy
// of G is freed), and Apply() then reads the owned slip and writes the
// tractions on the device: cuBLAS dgemv on a CUDA device, an mfem::forall
// row loop on the debug device.  One rank needs no MPI and no host copy;
// several ranks pass the all-gather and the return to the owners device
// pointers with CUDA-aware MPI (--qd-green-gpu-mpi), or stage them through
// the host.  Double precision only.

#ifndef MFEM_SEAS_QD_GREEN_OPERATOR_HPP
#define MFEM_SEAS_QD_GREEN_OPERATOR_HPP

#include "mfem.hpp"

#include <string>

namespace mfem
{
namespace seas
{

/// Settings of the Green's-function stage (--qd-green*).
struct QDGreenConfig
{
   int block = 64;           ///< unit slips per blocked solve while building
   std::string cache_path;   ///< read G here if present, else build + write
   double tol = 1e-8;        ///< validation: max rel. difference to a solve
   bool gpu_mpi = false;     ///< G3c: device pointers to MPI (CUDA-aware)
};

} // namespace seas
} // namespace mfem

#ifdef MFEM_USE_MPI

#include "../common/stage_timers.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>

namespace mfem
{
namespace seas
{

class QDGreenOperator
{
public:
   /// @a nf_glob fault DOFs in all; this rank owns @a nown of them, global
   /// indices [fault_first, fault_first + nown).  Collective.
   QDGreenOperator(MPI_Comm comm, HYPRE_BigInt nf_glob,
                   HYPRE_BigInt fault_first, int nown, bool with_normal)
      : comm_(comm), nf_(nf_glob), ncols_(2 * nf_glob),
        ncomp_(with_normal ? 3 : 2), nown_(nown)
   {
      MPI_Comm_size(comm_, &nranks_);
      MPI_Comm_rank(comm_, &rank_);
      // Owned ranges of every rank
      long long mine[2] = {static_cast<long long>(fault_first), nown};
      std::vector<long long> all(2 * nranks_);
      MPI_Allgather(mine, 2, MPI_LONG_LONG, all.data(), 2, MPI_LONG_LONG,
                    comm_);
      own_lo_.resize(nranks_);
      own_n_.resize(nranks_);
      for (int r = 0; r < nranks_; r++)
      {
         own_lo_[r] = all[2 * r];
         own_n_[r] = all[2 * r + 1];
      }
      // Balanced row ranges over all ranks
      row_lo_.resize(nranks_ + 1);
      for (int r = 0; r <= nranks_; r++)
      {
         row_lo_[r] = (static_cast<long long>(nf_) * r) / nranks_;
      }
      nb_ = static_cast<int>(row_lo_[rank_ + 1] - row_lo_[rank_]);
      g_size_ = static_cast<size_t>(nb_) * ncomp_ * ncols_;
      G_.assign(g_size_, 0.0);
      TD_.assign(static_cast<size_t>(nb_) * ncomp_, 0.0);
      have_cols_.assign(static_cast<size_t>(ncols_), 0);

      // Slip all-gather: rank r's owned slip is global [2·lo, 2·(lo + n))
      slip_counts_.resize(nranks_);
      slip_displs_.resize(nranks_);
      for (int r = 0; r < nranks_; r++)
      {
         slip_counts_[r] = static_cast<int>(2 * own_n_[r]);
         slip_displs_[r] = static_cast<int>(2 * own_lo_[r]);
      }
      // Fault DOFs this rank owns that rank r holds rows of (to_rows_),
      // and fault DOFs rank r owns that this rank holds rows of
      // (from_owners_); both contiguous in g.
      to_rows_lo_.resize(nranks_);
      to_rows_n_.resize(nranks_);
      from_own_lo_.resize(nranks_);
      from_own_n_.resize(nranks_);
      for (int r = 0; r < nranks_; r++)
      {
         Overlap(own_lo_[rank_], own_n_[rank_], row_lo_[r],
                 row_lo_[r + 1] - row_lo_[r], to_rows_lo_[r], to_rows_n_[r]);
         Overlap(own_lo_[r], own_n_[r], row_lo_[rank_], nb_,
                 from_own_lo_[r], from_own_n_[r]);
      }
      delta_all_.SetSize(static_cast<int>(ncols_));
   }

   ~QDGreenOperator()
   {
#ifdef MFEM_USE_CUDA
      if (d_G_alloc_) { CuMemFree(d_G_alloc_); }
#endif
   }

   QDGreenOperator(const QDGreenOperator &) = delete;
   QDGreenOperator &operator=(const QDGreenOperator &) = delete;

   int NComp() const { return ncomp_; }
   HYPRE_BigInt NumFaultDofs() const { return nf_; }
   HYPRE_BigInt NumColumns() const { return ncols_; }
   /// Rows (fault DOFs × ncomp) this rank holds.
   int NumLocalRows() const { return nb_ * ncomp_; }
   /// Bytes of G on this rank.
   double LocalBytes() const { return double(g_size_) * sizeof(real_t); }

   /// GPU plan G3c: move G and T_D to the device.  From then on Apply()
   /// reads the owned slip and writes the tractions on the device, and G
   /// can no longer be set, read or written (call WriteFile() first).  On
   /// a CUDA device G goes into one device allocation (it can pass 2^31
   /// values, beyond an mfem::Vector) and its host copy is freed; on the
   /// debug device it stays where it is.  @a gpu_mpi: give MPI device
   /// pointers (CUDA-aware MPI) when there is more than one rank, else
   /// stage the fault-sized messages through the host.  Collective.
   void UseDevice(bool gpu_mpi)
   {
      MFEM_VERIFY(!on_device_, "QD Green: UseDevice() called twice");
      MFEM_VERIFY(Device::Allows(Backend::DEVICE_MASK),
                  "QD Green: UseDevice() needs a device backend (--device)");
      MFEM_VERIFY(!Device::Allows(Backend::HIP_MASK), "QD Green: the "
                  "device product is implemented for CUDA (and the debug "
                  "device), not HIP");
#ifdef MFEM_USE_SINGLE
      MFEM_ABORT("QD Green: the device product is double precision only");
#endif
      const int nrows = nb_ * ncomp_;
      TD_dev_.UseDevice(true);
      TD_dev_.SetSize(nrows);
      {
         real_t *td = TD_dev_.HostWrite();
         for (int i = 0; i < nrows; i++) { td[i] = TD_[i]; }
      }
      TD_dev_.Read();
      y_dev_.UseDevice(true);
      y_dev_.SetSize(nrows);
      recv_dev_.UseDevice(true);
      recv_dev_.SetSize(nown_ * ncomp_);
      delta_all_.UseDevice(true);

      // The return to the owners sends y as it is: this rank's rows, in g
      // order, are the owners' pieces in rank order; and receives the owned
      // rows in g order.  Check that, and keep the counts.
      a2a_scount_.assign(nranks_, 0);
      a2a_sdispl_.assign(nranks_, 0);
      a2a_rcount_.assign(nranks_, 0);
      a2a_rdispl_.assign(nranks_, 0);
      int stot = 0, rtot = 0;
      for (int r = 0; r < nranks_; r++)
      {
         a2a_sdispl_[r] = stot;
         a2a_scount_[r] = static_cast<int>(from_own_n_[r] * ncomp_);
         MFEM_VERIFY(from_own_n_[r] == 0 ||
                     (from_own_lo_[r] - row_lo_[rank_]) * ncomp_ == stot,
                     "QD Green: rows of rank " << rank_ << " are not in the "
                     "owners' rank order");
         stot += a2a_scount_[r];
         a2a_rdispl_[r] = rtot;
         a2a_rcount_[r] = static_cast<int>(to_rows_n_[r] * ncomp_);
         MFEM_VERIFY(to_rows_n_[r] == 0 ||
                     (to_rows_lo_[r] - own_lo_[rank_]) * ncomp_ == rtot,
                     "QD Green: owned rows of rank " << rank_ << " do not "
                     "arrive in g order");
         rtot += a2a_rcount_[r];
      }
      MFEM_VERIFY(stot == nrows && rtot == nown_ * ncomp_, "QD Green: the "
                  "return to the owners sends " << stot << " of " << nrows
                  << " rows and receives " << rtot << " of "
                  << nown_ * ncomp_);

#ifdef MFEM_USE_CUDA
      if (Device::Allows(Backend::CUDA_MASK))
      {
         const size_t bytes = g_size_ * sizeof(real_t);
         if (bytes > 0)
         {
            CuMemAlloc(&d_G_alloc_, bytes);
            CuMemcpyHtoD(d_G_alloc_, G_.data(), bytes);
         }
         std::vector<real_t>().swap(G_);
         dG_ = static_cast<const real_t *>(d_G_alloc_);
         product_ = "cuBLAS dgemv";
      }
      else
#endif
      {
         dG_ = G_.data();
         product_ = "mfem::forall row loop";
      }
      gpu_mpi_ = gpu_mpi;
      on_device_ = true;
   }

   /// Whether UseDevice() was called.
   bool OnDevice() const { return on_device_; }
   /// How Apply() forms the product: "host loop", "cuBLAS dgemv" or
   /// "mfem::forall row loop".
   const char *ProductName() const { return product_; }

   /// Position keys of the owned fault DOFs, 6 per DOF (node x, y, z, then
   /// the centroid x, y, z of its fault face), for WriteFile and ReadFile.
   void SetKeys(const Vector &keys_owned)
   {
      MFEM_VERIFY(keys_owned.Size() == kKey * nown_, "QD Green: "
                  << keys_owned.Size() << " key values, expected "
                  << kKey * nown_);
      const real_t *k = keys_owned.HostRead();
      keys_.assign(k, k + kKey * nown_);
   }

   /// Columns [j0, j0 + m) of G from the owned tractions of m unit slips:
   /// @a Ts[i] shear (2·nown), @a Tn[i] normal (nown, with elastic
   /// sigma_n).  Collective; every rank passes the same j0 and m.
   void SetColumns(HYPRE_BigInt j0, int m, const std::vector<Vector> &Ts,
                   const std::vector<Vector> &Tn)
   {
      MFEM_VERIFY(!on_device_, "QD Green: SetColumns() after UseDevice()");
      MFEM_VERIFY(j0 >= 0 && j0 + m <= ncols_,
                  "QD Green: columns [" << j0 << ", " << j0 + m
                  << ") outside [0, " << ncols_ << ")");
      std::vector<real_t> recv;
      ToRows(m, Ts, Tn, recv);
      size_t k = 0;
      for (int r = 0; r < nranks_; r++)
      {
         for (long long g = from_own_lo_[r];
              g < from_own_lo_[r] + from_own_n_[r]; g++)
         {
            const size_t row0 = static_cast<size_t>(g - row_lo_[rank_]) *
                                ncomp_;
            for (int i = 0; i < m; i++)
            {
               for (int c = 0; c < ncomp_; c++)
               {
                  G_[(row0 + c) * ncols_ + j0 + i] = recv[k++];
               }
            }
         }
      }
      for (int i = 0; i < m; i++) { have_cols_[j0 + i] = 1; }
   }

   /// T_D from the owned tractions of the loading ĝ with zero slip.
   /// Collective.
   void SetLoading(const Vector &TDs, const Vector &TDn)
   {
      MFEM_VERIFY(!on_device_, "QD Green: SetLoading() after UseDevice()");
      std::vector<Vector> Ts(1), Tn(1);
      Ts[0] = TDs;
      Tn[0] = TDn;
      std::vector<real_t> recv;
      ToRows(1, Ts, Tn, recv);
      size_t k = 0;
      for (int r = 0; r < nranks_; r++)
      {
         for (long long g = from_own_lo_[r];
              g < from_own_lo_[r] + from_own_n_[r]; g++)
         {
            for (int c = 0; c < ncomp_; c++)
            {
               TD_[static_cast<size_t>(g - row_lo_[rank_]) * ncomp_ + c] =
                  recv[k++];
            }
         }
      }
      have_loading_ = true;
   }

   /// True when every column and T_D are set (built or read).
   bool Complete() const
   {
      int ok = have_loading_ ? 1 : 0;
      for (char h : have_cols_) { if (!h) { ok = 0; break; } }
      int all_ok = 0;
      MPI_Allreduce(&ok, &all_ok, 1, MPI_INT, MPI_MIN, comm_);
      return all_ok == 1;
   }

   /// [T; T_n] = G·δ + t·T_D on owned fault DOFs.  @a slip_owned: 2·nown;
   /// @a T_owned: 2·nown; @a Tn_owned: nown (with elastic sigma_n) or
   /// nullptr.  Times the all-gather as SlipExpand, the product as Solve
   /// and the return to the owners as Traction.  Collective.
   void Apply(real_t t, const Vector &slip_owned, Vector &T_owned,
              Vector *Tn_owned, StageTimers *timers) const
   {
      MFEM_VERIFY(slip_owned.Size() == 2 * nown_,
                  "QD Green: slip size " << slip_owned.Size()
                  << ", expected " << 2 * nown_);
      MFEM_VERIFY(Tn_owned == nullptr || ncomp_ == 3,
                  "QD Green: normal traction requested, but G was built "
                  "without the normal rows");
      if (on_device_)
      {
         ApplyDevice(t, slip_owned, T_owned, Tn_owned, timers);
         return;
      }
      {
         ScopedStageTimer timer(timers, StageTimers::SlipExpand);
         const real_t *s = slip_owned.HostRead();
         MPI_Allgatherv(s, 2 * nown_, MPITypeMap<real_t>::mpi_type,
                        delta_all_.HostWrite(), slip_counts_.data(),
                        slip_displs_.data(), MPITypeMap<real_t>::mpi_type,
                        comm_);
      }
      {
         ScopedStageTimer timer(timers, StageTimers::Solve);
         y_.resize(static_cast<size_t>(nb_) * ncomp_);
         const real_t *d = delta_all_.HostRead();
         const size_t nrows = y_.size();
         const size_t n = static_cast<size_t>(ncols_);
         for (size_t row = 0; row < nrows; row++)
         {
            const real_t *a = G_.data() + row * n;
            // Four partial sums: independent chains, so the loop streams
            // G at memory speed.
            real_t s0 = 0.0, s1 = 0.0, s2 = 0.0, s3 = 0.0;
            size_t j = 0;
            for (; j + 4 <= n; j += 4)
            {
               s0 += a[j] * d[j];
               s1 += a[j + 1] * d[j + 1];
               s2 += a[j + 2] * d[j + 2];
               s3 += a[j + 3] * d[j + 3];
            }
            for (; j < n; j++) { s0 += a[j] * d[j]; }
            y_[row] = (s0 + s1) + (s2 + s3) + t * TD_[row];
         }
      }
      ScopedStageTimer timer(timers, StageTimers::Traction);
      FromRows(T_owned, Tn_owned);
   }

   /// Write G and T_D to @a path: a header, then the rows of every fault
   /// DOF in global order, then T_D.  The layout does not depend on the
   /// number of ranks.  Collective.
   ///
   /// No two ranks write at the same time where their blocks meet.  With
   /// Open MPI 5's MPI-IO (OMPIO), independent writes by neighbouring ranks
   /// into the same file page lost one rank's part of that page: Vista job
   /// 1080407 (/scratch on NFS) got two runs of zero keys (447 of 27,936
   /// fault DOFs), and a Mac retest got 1,792 zero bytes at the start of
   /// rank 5's G rows, up to the next 4 KiB page boundary.  So even ranks
   /// write their rows, then odd ranks, then rank 0 writes the header, the
   /// keys and T_D (gathered), with sync, barrier, sync between the phases.
   /// Afterwards every rank reads back the first and last 64 KiB of its
   /// rows, and rank 0 the keys and T_D; on any difference the file is
   /// deleted and the run stops.
   void WriteFile(const std::string &path, const std::string &info) const
   {
      MFEM_VERIFY(!on_device_, "QD Green: WriteFile() after UseDevice() (the "
                  "host copy of G may be gone)");
      MFEM_VERIFY(keys_.size() == static_cast<size_t>(kKey) * nown_,
                  "QD Green: SetKeys() was not called before WriteFile()");
      for (int o = 0; o < nown_; o++)
      {
         const real_t *k = &keys_[static_cast<size_t>(kKey) * o];
         bool finite = true, apart = false;
         for (int c = 0; c < kKey; c++)
         {
            finite = finite && std::isfinite(k[c]);
         }
         for (int d = 0; d < 3; d++) { apart = apart || k[d] != k[3 + d]; }
         MFEM_VERIFY(finite && apart, "QD Green: the position key of owned "
                     "fault DOF " << own_lo_[rank_] + o << " on rank "
                     << rank_ << " is " << (finite ? "its face centroid"
                                            : "not finite")
                     << " (node " << k[0] << ", " << k[1] << ", " << k[2]
                     << "); not writing '" << path << "'");
      }
      // Keys in owned order and T_D in row order, on rank 0
      std::vector<real_t> keys_all, td_all;
      Gather(keys_, own_lo_, own_n_, kKey, keys_all);
      std::vector<long long> nrow(nranks_);
      for (int r = 0; r < nranks_; r++)
      {
         nrow[r] = row_lo_[r + 1] - row_lo_[r];
      }
      Gather(TD_, row_lo_, nrow, ncomp_, td_all);

      MPI_File fh;
      int err = MPI_File_open(comm_, path.c_str(),
                              MPI_MODE_CREATE | MPI_MODE_WRONLY,
                              MPI_INFO_NULL, &fh);
      MFEM_VERIFY(err == MPI_SUCCESS, "QD Green: cannot create cache file '"
                  << path << "' (MPI_File_open error " << err << ")");
      MPI_File_set_size(fh, 0);
      auto settle = [&]()
      {
         CheckIO(MPI_File_sync(fh), path, "sync");
         MPI_Barrier(comm_);
         CheckIO(MPI_File_sync(fh), path, "sync");
      };
      for (int phase = 0; phase < 2; phase++)
      {
         if (rank_ % 2 == phase)
         {
            WriteChunked(fh, RowOffset(row_lo_[rank_]), G_.data(), G_.size(),
                         path);
         }
         settle();
      }
      if (rank_ == 0)
      {
         Header h = MakeHeader();
         std::strncpy(h.info, info.c_str(), sizeof(h.info) - 1);
         CheckIO(MPI_File_write_at(fh, 0, &h, sizeof(Header), MPI_BYTE,
                                   MPI_STATUS_IGNORE), path, "header");
         WriteChunked(fh, KeyOffset(0), keys_all.data(), keys_all.size(),
                      path);
         WriteChunked(fh, TDOffset(0), td_all.data(), td_all.size(), path);
      }
      settle();
      CheckIO(MPI_File_close(&fh), path, "close");

      // Read back the ends of every rank's rows, and the keys and T_D
      CheckIO(MPI_File_open(comm_, path.c_str(), MPI_MODE_RDONLY,
                            MPI_INFO_NULL, &fh), path, "reopen");
      int bad = 0;
      const size_t edge = std::min(G_.size(), size_t(8192));
      std::vector<real_t> back(edge);
      if (edge > 0)
      {
         ReadChunked(fh, RowOffset(row_lo_[rank_]), back.data(), edge, path);
         bad |= std::memcmp(back.data(), G_.data(),
                            edge * sizeof(real_t)) != 0;
         const size_t last = G_.size() - edge;
         ReadChunked(fh, RowOffset(row_lo_[rank_]) +
                     MPI_Offset(last * sizeof(real_t)), back.data(), edge,
                     path);
         bad |= std::memcmp(back.data(), G_.data() + last,
                            edge * sizeof(real_t)) != 0;
      }
      if (rank_ == 0)
      {
         back.resize(keys_all.size());
         ReadChunked(fh, KeyOffset(0), back.data(), back.size(), path);
         bad |= 2 * (std::memcmp(back.data(), keys_all.data(),
                                 back.size() * sizeof(real_t)) != 0);
         back.resize(td_all.size());
         ReadChunked(fh, TDOffset(0), back.data(), back.size(), path);
         bad |= 4 * (std::memcmp(back.data(), td_all.data(),
                                 back.size() * sizeof(real_t)) != 0);
      }
      MPI_File_close(&fh);
      int bad_all = 0, first_bad = bad ? rank_ : nranks_;
      MPI_Allreduce(&bad, &bad_all, 1, MPI_INT, MPI_BOR, comm_);
      MPI_Allreduce(MPI_IN_PLACE, &first_bad, 1, MPI_INT, MPI_MIN, comm_);
      if (bad_all != 0 && rank_ == 0)
      {
         MPI_File_delete(path.c_str(), MPI_INFO_NULL);
      }
      MFEM_VERIFY(bad_all == 0, "QD Green: read back from '" << path
                  << "', " << ((bad_all & 1) ? "G rows " : "")
                  << ((bad_all & 2) ? "keys " : "")
                  << ((bad_all & 4) ? "T_D " : "")
                  << "differ from what was written (first on rank "
                  << first_bad << "); the file was deleted");
   }

   /// Read G and T_D from @a path, written by WriteFile for the same
   /// number of fault DOFs and components (any number of ranks).  Returns
   /// false when the file does not exist; aborts when it does not fit.
   /// Collective.
   bool ReadFile(const std::string &path, std::string &info)
   {
      MFEM_VERIFY(!on_device_, "QD Green: ReadFile() after UseDevice()");
      MPI_File fh;
      int err = MPI_File_open(comm_, path.c_str(), MPI_MODE_RDONLY,
                              MPI_INFO_NULL, &fh);
      if (err != MPI_SUCCESS) { return false; }
      Header h;
      if (rank_ == 0)
      {
         CheckIO(MPI_File_read_at(fh, 0, &h, sizeof(Header), MPI_BYTE,
                                  MPI_STATUS_IGNORE), path, "header");
      }
      MPI_Bcast(&h, sizeof(Header), MPI_BYTE, 0, comm_);
      const Header want = MakeHeader();
      MFEM_VERIFY(std::memcmp(h.magic, want.magic, sizeof(h.magic)) == 0,
                  "QD Green: '" << path << "' is not a QD Green cache file of "
                  "this version (SEAS-QD-GREEN-2)");
      MFEM_VERIFY(h.nf == want.nf && h.ncomp == want.ncomp &&
                  h.ncols == want.ncols,
                  "QD Green: cache file '" << path << "' holds N_f = "
                  << h.nf << ", " << h.ncomp << " components, " << h.ncols
                  << " columns; this model has N_f = " << want.nf << ", "
                  << want.ncomp << ", " << want.ncols
                  << ".  Remove the file or pass another --qd-green-cache.");
      MPI_Offset size = 0;
      MPI_File_get_size(fh, &size);
      MFEM_VERIFY(size == TDOffset(nf_),
                  "QD Green: cache file '" << path << "' has " << size
                  << " bytes, expected " << TDOffset(nf_)
                  << " (incomplete write?)");
      h.info[sizeof(h.info) - 1] = '\0';
      info = h.info;
      // File row of every fault DOF of this partition, by key
      std::vector<long long> perm;
      MatchKeys(fh, path, perm);
      bool identity = true;
      for (long long g = 0; g < nf_; g++)
      {
         if (perm[g] != g) { identity = false; break; }
      }
      if (identity)
      {
         ReadChunked(fh, RowOffset(row_lo_[rank_]), G_.data(), G_.size(),
                     path);
         ReadChunked(fh, TDOffset(row_lo_[rank_]), TD_.data(), TD_.size(),
                     path);
      }
      else
      {
         // Row block of fault DOF g is file block perm[g]; column 2·h + c
         // is file column 2·perm[h] + c.
         const size_t n = static_cast<size_t>(ncols_);
         std::vector<real_t> buf(static_cast<size_t>(ncomp_) * n);
         for (int b = 0; b < nb_; b++)
         {
            const long long gf = perm[row_lo_[rank_] + b];
            ReadChunked(fh, RowOffset(gf), buf.data(), buf.size(), path);
            for (int c = 0; c < ncomp_; c++)
            {
               real_t *row = G_.data() + (static_cast<size_t>(b) * ncomp_ + c)
                             * n;
               const real_t *src = buf.data() + static_cast<size_t>(c) * n;
               for (long long hh = 0; hh < nf_; hh++)
               {
                  const size_t jf = static_cast<size_t>(2 * perm[hh]);
                  row[2 * hh] = src[jf];
                  row[2 * hh + 1] = src[jf + 1];
               }
            }
            ReadChunked(fh, TDOffset(gf),
                        TD_.data() + static_cast<size_t>(b) * ncomp_,
                        static_cast<size_t>(ncomp_), path);
         }
      }
      MPI_File_close(&fh);
      std::fill(have_cols_.begin(), have_cols_.end(), 1);
      have_loading_ = true;
      return true;
   }

private:
   struct Header
   {
      char magic[16];
      long long nf, ncomp, ncols;
      long long reserved[5];
      char info[432];
   };
   static_assert(sizeof(Header) == 512, "QD Green cache header is 512 B");

   Header MakeHeader() const
   {
      Header h;
      std::memset(&h, 0, sizeof(Header));
      std::memcpy(h.magic, "SEAS-QD-GREEN-2", 15);
      h.nf = nf_;
      h.ncomp = ncomp_;
      h.ncols = ncols_;
      return h;
   }

   // File: header, keys (kKey per fault DOF), G rows, T_D; all in the
   // writer's global order.
   static constexpr int kKey = 6;
   MPI_Offset KeyOffset(long long g) const
   {
      return static_cast<MPI_Offset>(sizeof(Header)) +
             static_cast<MPI_Offset>(g) * kKey * sizeof(real_t);
   }
   MPI_Offset RowOffset(long long g) const
   {
      return KeyOffset(nf_) +
             static_cast<MPI_Offset>(g) * ncomp_ * ncols_ * sizeof(real_t);
   }
   MPI_Offset TDOffset(long long g) const
   {
      return RowOffset(nf_) +
             static_cast<MPI_Offset>(g) * ncomp_ * sizeof(real_t);
   }

   static void CheckIO(int err, const std::string &path, const char *what)
   {
      MFEM_VERIFY(err == MPI_SUCCESS, "QD Green: I/O error " << err
                  << " on " << what << " of '" << path << "'");
   }

   /// @a all on rank 0 = every rank's @a mine (@a n[r] entries of @a w
   /// values each), placed at entry @a lo[r].  Collective.
   void Gather(const std::vector<real_t> &mine,
               const std::vector<long long> &lo,
               const std::vector<long long> &n, int w,
               std::vector<real_t> &all) const
   {
      std::vector<int> cnt(nranks_), dsp(nranks_);
      long long total = 0;
      for (int r = 0; r < nranks_; r++)
      {
         cnt[r] = static_cast<int>(w * n[r]);
         dsp[r] = static_cast<int>(w * lo[r]);
         total += n[r];
      }
      MFEM_VERIFY(static_cast<long long>(mine.size()) == w * n[rank_],
                  "QD Green: gather size mismatch on rank " << rank_);
      all.assign(rank_ == 0 ? static_cast<size_t>(w * total) : 0, 0.0);
      MPI_Gatherv(mine.data(), cnt[rank_], MPITypeMap<real_t>::mpi_type,
                  all.data(), cnt.data(), dsp.data(),
                  MPITypeMap<real_t>::mpi_type, 0, comm_);
   }

   // Independent writes in pieces of at most 2^27 values (1 GiB), so no
   // count overflows an int.
   static constexpr size_t kChunk = size_t(1) << 27;
   static void WriteChunked(MPI_File fh, MPI_Offset off, const real_t *p,
                            size_t n, const std::string &path)
   {
      for (size_t k = 0; k < n; k += kChunk)
      {
         const int c = static_cast<int>(std::min(kChunk, n - k));
         CheckIO(MPI_File_write_at(fh, off + MPI_Offset(k * sizeof(real_t)),
                                   p + k, c, MPITypeMap<real_t>::mpi_type,
                                   MPI_STATUS_IGNORE), path, "G rows");
      }
   }
   static void ReadChunked(MPI_File fh, MPI_Offset off, real_t *p, size_t n,
                           const std::string &path)
   {
      for (size_t k = 0; k < n; k += kChunk)
      {
         const int c = static_cast<int>(std::min(kChunk, n - k));
         CheckIO(MPI_File_read_at(fh, off + MPI_Offset(k * sizeof(real_t)),
                                  p + k, c, MPITypeMap<real_t>::mpi_type,
                                  MPI_STATUS_IGNORE), path, "G rows");
      }
   }

   /// @a perm[g] = file fault DOF of this partition's fault DOF g, found
   /// by matching keys to within 1e-9 of the largest coordinate.  Aborts
   /// unless every fault DOF matches exactly one file DOF.  Collective.
   void MatchKeys(MPI_File fh, const std::string &path,
                  std::vector<long long> &perm) const
   {
      MFEM_VERIFY(keys_.size() == static_cast<size_t>(kKey) * nown_,
                  "QD Green: SetKeys() was not called before ReadFile()");
      std::vector<real_t> fk(static_cast<size_t>(kKey) * nf_);
      if (rank_ == 0)
      {
         ReadChunked(fh, KeyOffset(0), fk.data(), fk.size(), path);
      }
      for (size_t k = 0; k < fk.size(); k += kChunk)
      {
         MPI_Bcast(fk.data() + k, static_cast<int>(std::min(kChunk,
                                                            fk.size() - k)),
                   MPITypeMap<real_t>::mpi_type, 0, comm_);
      }
      real_t ext = 1.0;
      for (real_t v : fk) { ext = std::max(ext, std::abs(v)); }
      const real_t tol = 1e-9 * ext;
      // Sort the file DOFs along a generic direction, then search a window
      static const real_t dir[kKey] = {0.5381, 0.3112, 0.7834,
                                       0.1297, 0.2269, 0.4513};
      real_t dsum = 0.0;
      for (int c = 0; c < kKey; c++) { dsum += dir[c]; }
      auto proj = [&](const real_t *k)
      {
         real_t p = 0.0;
         for (int c = 0; c < kKey; c++) { p += dir[c] * k[c]; }
         return p;
      };
      std::vector<real_t> pf(nf_);
      std::vector<long long> order(nf_);
      for (long long g = 0; g < nf_; g++)
      {
         pf[g] = proj(&fk[static_cast<size_t>(kKey) * g]);
         order[g] = g;
      }
      std::sort(order.begin(), order.end(),
                [&](long long a, long long b) { return pf[a] < pf[b]; });
      std::vector<real_t> ps(nf_);
      for (long long i = 0; i < nf_; i++) { ps[i] = pf[order[i]]; }
      std::vector<long long> mine(nown_, -1);
      long long bad = 0;
      for (int o = 0; o < nown_; o++)
      {
         const real_t *ko = &keys_[static_cast<size_t>(kKey) * o];
         const real_t po = proj(ko);
         int nmatch = 0;
         for (auto it = std::lower_bound(ps.begin(), ps.end(),
                                         po - 2.0 * tol * dsum);
              it != ps.end() && *it <= po + 2.0 * tol * dsum; ++it)
         {
            const long long gf = order[it - ps.begin()];
            const real_t *kf = &fk[static_cast<size_t>(kKey) * gf];
            real_t dmax = 0.0;
            for (int c = 0; c < kKey; c++)
            {
               dmax = std::max(dmax, std::abs(ko[c] - kf[c]));
            }
            if (dmax <= tol) { mine[o] = gf; nmatch++; }
         }
         if (nmatch != 1) { bad++; }
      }
      long long bad_all = 0;
      MPI_Allreduce(&bad, &bad_all, 1, MPI_LONG_LONG, MPI_SUM, comm_);
      MFEM_VERIFY(bad_all == 0, "QD Green: " << bad_all << " of " << nf_
                  << " fault DOFs of this mesh have no unique match (by "
                  "position, tolerance " << tol << ") in cache file '"
                  << path << "'.  Remove the file or pass another "
                  "--qd-green-cache.");
      perm.assign(nf_, -1);
      std::vector<int> cnt(nranks_), dsp(nranks_);
      for (int r = 0; r < nranks_; r++)
      {
         cnt[r] = static_cast<int>(own_n_[r]);
         dsp[r] = static_cast<int>(own_lo_[r]);
      }
      MPI_Allgatherv(mine.data(), nown_, MPI_LONG_LONG, perm.data(),
                     cnt.data(), dsp.data(), MPI_LONG_LONG, comm_);
      std::vector<char> seen(nf_, 0);
      for (long long g = 0; g < nf_; g++)
      {
         MFEM_VERIFY(perm[g] >= 0 && perm[g] < nf_ && !seen[perm[g]],
                     "QD Green: cache file '" << path << "' maps two fault "
                     "DOFs of this mesh to one file row.  Remove the file "
                     "or pass another --qd-green-cache.");
         seen[perm[g]] = 1;
      }
   }

   static void Overlap(long long a, long long na, long long b, long long nb,
                       long long &lo, long long &n)
   {
      lo = std::max(a, b);
      n = std::max(0LL, std::min(a + na, b + nb) - lo);
   }

   /// Owned tractions of m columns -> this rank's rows.  @a recv holds, for
   /// each source rank in order, for g in from_own, for i < m, ncomp values.
   void ToRows(int m, const std::vector<Vector> &Ts,
               const std::vector<Vector> &Tn, std::vector<real_t> &recv) const
   {
      MFEM_VERIFY(static_cast<int>(Ts.size()) >= m &&
                  (ncomp_ == 2 || static_cast<int>(Tn.size()) >= m),
                  "QD Green: " << m << " columns, got " << Ts.size()
                  << " shear and " << Tn.size() << " normal vectors");
      std::vector<int> scount(nranks_), sdispl(nranks_), rcount(nranks_),
          rdispl(nranks_);
      std::vector<real_t> send;
      send.reserve(static_cast<size_t>(nown_) * m * ncomp_);
      for (int i = 0; i < m; i++) { Ts[i].HostRead(); }
      if (ncomp_ == 3) { for (int i = 0; i < m; i++) { Tn[i].HostRead(); } }
      for (int r = 0; r < nranks_; r++)
      {
         sdispl[r] = static_cast<int>(send.size());
         for (long long g = to_rows_lo_[r]; g < to_rows_lo_[r] + to_rows_n_[r];
              g++)
         {
            const int o = static_cast<int>(g - own_lo_[rank_]);
            for (int i = 0; i < m; i++)
            {
               send.push_back(Ts[i](2 * o));
               send.push_back(Ts[i](2 * o + 1));
               if (ncomp_ == 3) { send.push_back(Tn[i](o)); }
            }
         }
         scount[r] = static_cast<int>(send.size()) - sdispl[r];
      }
      int total = 0;
      for (int r = 0; r < nranks_; r++)
      {
         rdispl[r] = total;
         rcount[r] = static_cast<int>(from_own_n_[r] * m * ncomp_);
         total += rcount[r];
      }
      recv.resize(total);
      MPI_Alltoallv(send.data(), scount.data(), sdispl.data(),
                    MPITypeMap<real_t>::mpi_type, recv.data(), rcount.data(),
                    rdispl.data(), MPITypeMap<real_t>::mpi_type, comm_);
   }

public:
   // Public only because nvcc does not allow an extended __host__ __device__
   // lambda inside a private member function; call Apply().

   /// Apply() after UseDevice(): the same three steps on device vectors.
   void ApplyDevice(real_t t, const Vector &slip_owned, Vector &T_owned,
                    Vector *Tn_owned, StageTimers *timers) const
   {
      const int nrows = nb_ * ncomp_;
      const real_t *d = nullptr;
      {
         ScopedStageTimer timer(timers, StageTimers::SlipExpand);
         if (nranks_ == 1)
         {
            d = slip_owned.Read(true);
         }
         else if (gpu_mpi_)
         {
            // Kernels that wrote the slip run asynchronously; MPI reads the
            // buffer on its own streams.
            const real_t *s = slip_owned.Read(true);
            MFEM_STREAM_SYNC;
            MPI_Allgatherv(s, 2 * nown_,
                           MPITypeMap<real_t>::mpi_type,
                           delta_all_.Write(true), slip_counts_.data(),
                           slip_displs_.data(), MPITypeMap<real_t>::mpi_type,
                           comm_);
            d = delta_all_.Read(true);
         }
         else
         {
            MPI_Allgatherv(slip_owned.HostRead(), 2 * nown_,
                           MPITypeMap<real_t>::mpi_type,
                           delta_all_.HostWrite(), slip_counts_.data(),
                           slip_displs_.data(), MPITypeMap<real_t>::mpi_type,
                           comm_);
            d = delta_all_.Read(true);
         }
      }
      {
         ScopedStageTimer timer(timers, StageTimers::Solve);
         real_t *y = y_dev_.Write(true);
         const real_t *td = TD_dev_.Read(true);
         mfem::forall(nrows, [=] MFEM_HOST_DEVICE (int i)
         {
            y[i] = t * td[i];
         });
         RowProductAdd(d, y);
      }
      ScopedStageTimer timer(timers, StageTimers::Traction);
      const real_t *src = nullptr;
      if (nranks_ == 1)
      {
         src = y_dev_.Read(true);
      }
      else if (gpu_mpi_)
      {
         const real_t *yd = y_dev_.Read(true);
         MFEM_STREAM_SYNC;
         MPI_Alltoallv(yd, a2a_scount_.data(),
                       a2a_sdispl_.data(), MPITypeMap<real_t>::mpi_type,
                       recv_dev_.Write(true), a2a_rcount_.data(),
                       a2a_rdispl_.data(), MPITypeMap<real_t>::mpi_type,
                       comm_);
         src = recv_dev_.Read(true);
      }
      else
      {
         MPI_Alltoallv(y_dev_.HostRead(), a2a_scount_.data(),
                       a2a_sdispl_.data(), MPITypeMap<real_t>::mpi_type,
                       recv_dev_.HostWrite(), a2a_rcount_.data(),
                       a2a_rdispl_.data(), MPITypeMap<real_t>::mpi_type,
                       comm_);
         src = recv_dev_.Read(true);
      }
      T_owned.SetSize(2 * nown_);
      real_t *T = T_owned.Write(true);
      real_t *Tn = nullptr;
      if (Tn_owned)
      {
         Tn_owned->SetSize(nown_);
         Tn = Tn_owned->Write(true);
      }
      const int nc = ncomp_;
      mfem::forall(nown_, [=] MFEM_HOST_DEVICE (int o)
      {
         T[2 * o] = src[nc * o];
         T[2 * o + 1] = src[nc * o + 1];
         if (Tn) { Tn[o] = src[nc * o + 2]; }
      });
   }

   /// y += G_rows·d on the device (d: 2·N_f values).
   void RowProductAdd(const real_t *d, real_t *y) const
   {
      const int nrows = nb_ * ncomp_;
      if (nrows == 0) { return; }
#ifdef MFEM_USE_CUDA
      if (Device::Allows(Backend::CUDA_MASK))
      {
#ifndef MFEM_USE_SINGLE
         // G is row-major, so as a column-major matrix it is G^T with
         // leading dimension N_cols: y = op(A)·d with op = transpose.
         const double one = 1.0;
         const cublasStatus_t st = cublasDgemv_64(
            GPUBlas::Handle(), CUBLAS_OP_T, static_cast<int64_t>(ncols_),
            static_cast<int64_t>(nrows), &one, dG_,
            static_cast<int64_t>(ncols_), d, 1, &one, y, 1);
         MFEM_VERIFY(st == CUBLAS_STATUS_SUCCESS, "QD Green: cublasDgemv_64 "
                     "failed with status " << static_cast<int>(st));
#endif
         return;
      }
#endif
      const real_t *G = dG_;
      const long long n = static_cast<long long>(ncols_);
      mfem::forall(nrows, [=] MFEM_HOST_DEVICE (int row)
      {
         const real_t *a = G + static_cast<size_t>(row) * n;
         real_t s = 0.0;
         for (long long j = 0; j < n; j++) { s += a[j] * d[j]; }
         y[row] += s;
      });
   }

private:
   /// y_ (this rank's rows) -> owned tractions.
   void FromRows(Vector &T_owned, Vector *Tn_owned) const
   {
      std::vector<int> scount(nranks_), sdispl(nranks_), rcount(nranks_),
          rdispl(nranks_);
      int stot = 0, rtot = 0;
      for (int r = 0; r < nranks_; r++)
      {
         sdispl[r] = stot;
         scount[r] = static_cast<int>(from_own_n_[r] * ncomp_);
         stot += scount[r];
         rdispl[r] = rtot;
         rcount[r] = static_cast<int>(to_rows_n_[r] * ncomp_);
         rtot += rcount[r];
      }
      // Rows go out in g order, which is source order for every target.
      send_.resize(stot);
      {
         size_t k = 0;
         for (int r = 0; r < nranks_; r++)
         {
            const size_t row0 = static_cast<size_t>(from_own_lo_[r] -
                                                    row_lo_[rank_]) * ncomp_;
            for (long long q = 0; q < from_own_n_[r] * ncomp_; q++)
            {
               send_[k++] = y_[row0 + q];
            }
         }
      }
      recv_.resize(rtot);
      MPI_Alltoallv(send_.data(), scount.data(), sdispl.data(),
                    MPITypeMap<real_t>::mpi_type, recv_.data(), rcount.data(),
                    rdispl.data(), MPITypeMap<real_t>::mpi_type, comm_);
      T_owned.SetSize(2 * nown_);
      real_t *T = T_owned.HostWrite();
      real_t *Tn = nullptr;
      if (Tn_owned)
      {
         Tn_owned->SetSize(nown_);
         Tn = Tn_owned->HostWrite();
      }
      size_t k = 0;
      for (int r = 0; r < nranks_; r++)
      {
         for (long long g = to_rows_lo_[r]; g < to_rows_lo_[r] + to_rows_n_[r];
              g++)
         {
            const int o = static_cast<int>(g - own_lo_[rank_]);
            T[2 * o] = recv_[k];
            T[2 * o + 1] = recv_[k + 1];
            if (Tn) { Tn[o] = recv_[k + 2]; }
            k += ncomp_;
         }
      }
   }

   MPI_Comm comm_;
   int nranks_ = 1, rank_ = 0;
   HYPRE_BigInt nf_, ncols_;
   int ncomp_;
   int nown_;
   std::vector<long long> own_lo_, own_n_, row_lo_;
   int nb_ = 0;
   std::vector<long long> to_rows_lo_, to_rows_n_, from_own_lo_, from_own_n_;
   std::vector<int> slip_counts_, slip_displs_;
   std::vector<real_t> G_;      ///< nb_·ncomp_ rows × ncols_, row-major
   std::vector<real_t> TD_;     ///< nb_·ncomp_
   std::vector<char> have_cols_;
   bool have_loading_ = false;
   std::vector<real_t> keys_;   ///< 6 per owned fault DOF (SetKeys)
   mutable Vector delta_all_;
   mutable std::vector<real_t> y_, send_, recv_;
   size_t g_size_ = 0;          ///< entries of G on this rank

   // GPU plan G3c (UseDevice)
   bool on_device_ = false;
   bool gpu_mpi_ = false;
   const char *product_ = "host loop";
   const real_t *dG_ = nullptr;   ///< G for the device product
   void *d_G_alloc_ = nullptr;    ///< CUDA allocation behind dG_, or null
   Vector TD_dev_;
   mutable Vector y_dev_, recv_dev_;
   std::vector<int> a2a_scount_, a2a_sdispl_, a2a_rcount_, a2a_rdispl_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_USE_MPI

#endif // MFEM_SEAS_QD_GREEN_OPERATOR_HPP

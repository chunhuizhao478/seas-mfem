// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// nitsche_device_fault.hpp — GPU dynamic-rupture plan D2 (plan rev2 §4.3,
// "Fault on the host, for this phase only"): the host fault fed from device
// vectors.
//
// In D2 the volume work of a step runs on the device and the fault (boxes
// 3-11) stays on the host.  Two small transfers per step connect them, both
// built on the per-fault-face operators of ElasticityDomainOperator
// (BuildFaultFaceOps, profiling plan item 3):
//
//   FaultGather     the values of U^{n+1} (and Udot*) that the fault face
//                   operators read: every fault face's element vdofs and the
//                   vdofs sent to fault-face neighbours (FaultNbrExchange_).
//                   They are gathered on the device into one buffer, copied
//                   to the host and scattered into host vectors of full
//                   size, of which the fault reads only these entries.
//   DeviceSlipLoad  rhs += G s on the device for the local kernel slip s
//                   (2 nbf values per fault face), which the host copies to
//                   the device, with the stored blocks G_f of
//                   AddSlipLoadFast.  The rows are summed per target vdof,
//                   one thread per vdof, in face order, so no entry is
//                   written twice at once and each entry gets the same sums
//                   in the same order as AddSlipLoadFast (the debug device
//                   gives its bits).
//
// Neither needs an assembled matrix.
//
// D3 (plan rev2 §4.4) moves the fault itself to the device and removes
// both copies.  Its pieces, all built from the same host structures:
//
//   DeviceFaultTraction  ComputeTractionFast: the face blocks A_f and their
//                        vdofs, and the fault-face neighbour exchange
//                        (FaultNbrExchange_) packed and unpacked on the
//                        device.
//   DeviceFaultMaps      the owned <-> local fault DOF maps:
//                        RestrictToOwnedFault's index, and ToKernelLocal_ of
//                        NitscheWaveOperator (ExpandOwnedToLocalFault with
//                        the FaultScatter ghost copies, then
//                        fault_frame::ToKernel).
//   DeviceSlipLoad       AddMultDevice: G s for a slip already on the device.
//
// Their messages leave from device memory with GPU-aware MPI
// (Device::GetGPUAwareMPI, as NitscheSipgOperator's halo), else through host
// buffers.  Every sum keeps the host's order, so the debug device gives the
// host's bits.

#ifndef MFEM_SEAS_NITSCHE_DEVICE_FAULT_HPP
#define MFEM_SEAS_NITSCHE_DEVICE_FAULT_HPP

#include "mfem.hpp"
#include "../common/mpi_tags.hpp"
#include "../fault/fault_frame.hpp"

#include <algorithm>
#include <vector>

namespace mfem
{
namespace seas
{

/// Appends the local vdofs that ComputeTractionFast of `op` reads (the side
/// vdofs of every fault face and the vdofs it sends to fault-face
/// neighbours).  `op` is an ElasticityDomainOperator with its face operators
/// built.
template <class Op>
inline void AppendFaultFaceVDofs(const Op &op, std::vector<int> &vd)
{
   const auto &f = op.GetFaultFaceOps();
   MFEM_VERIFY(f.built, "AppendFaultFaceVDofs: call BuildFaultFaceOps first");
   for (int v : f.vd) { if (v >= 0) { vd.push_back(v); } }   // side 2 of a shared face: -1
   for (int v : f.send_vd)
   {
      MFEM_VERIFY(v >= 0, "AppendFaultFaceVDofs: signed vdof (not DG?)");
      vd.push_back(v);
   }
}

class FaultGather
{
public:
   /// `vdofs`: the local vdofs to gather (any order, repeats allowed).
   explicit FaultGather(std::vector<int> vdofs)
   {
      std::sort(vdofs.begin(), vdofs.end());
      vdofs.erase(std::unique(vdofs.begin(), vdofs.end()), vdofs.end());
      host_ = vdofs;
      const int n = static_cast<int>(host_.size());
      idx_.SetSize(n);
      int *h = idx_.HostWrite();
      for (int i = 0; i < n; i++) { h[i] = host_[i]; }
      buf_.UseDevice(true);
   }

   int Size() const { return static_cast<int>(host_.size()); }

   /// xh(i) = x(i) for the gathered i (and vh(i) = v(i) when v is given);
   /// x and v may live on the device, xh and vh are host vectors of x's
   /// size whose other entries are left as they are.
   void Gather(const Vector &x, Vector &xh, const Vector *v = nullptr,
               Vector *vh = nullptr) const
   {
      const int n = Size(), nv = v ? 2 : 1;
      MFEM_VERIFY(xh.Size() == x.Size() && (!v || (vh && vh->Size() == v->Size())),
                  "FaultGather: host vectors must have the sizes of x and v");
      buf_.SetSize(nv * n);
      if (n > 0)
      {
         const int *id = idx_.Read();
         const real_t *dx = x.Read(), *dv = v ? v->Read() : nullptr;
         real_t *b = buf_.Write();
         mfem::forall(n, [=] MFEM_HOST_DEVICE (int i)
         {
            b[i] = dx[id[i]];
            if (dv) { b[n + i] = dv[id[i]]; }
         });
      }
      const real_t *hb = buf_.HostRead();
      for (int i = 0; i < n; i++) { xh(host_[i]) = hb[i]; }
      if (v)
      {
         for (int i = 0; i < n; i++) { (*vh)(host_[i]) = hb[n + i]; }
      }
   }

private:
   std::vector<int> host_;
   Array<int> idx_;
   mutable Vector buf_;
};

class DeviceSlipLoad
{
public:
   /// From the face operators of `op` (an ElasticityDomainOperator with
   /// BuildFaultFaceOps done); the blocks G_f are copied.
   template <class Op>
   explicit DeviceSlipLoad(const Op &op)
   {
      const auto &f = op.GetFaultFaceOps();
      MFEM_VERIFY(f.built, "DeviceSlipLoad: call BuildFaultFaceOps first");
      const int nbf = op.GetNbfPerFace();
      m_ = 2 * nbf;
      ns_ = 2 * f.nf * nbf;
      // one row per (face, side vdof) in AddSlipLoadFast's order, then
      // grouped by target vdof keeping that order
      struct Row { int t, g, s; };
      std::vector<Row> rows;
      for (int fc = 0; fc < f.nf; fc++)
      {
         const int ng = (fc < f.nfi) ? f.n1[fc] + f.n2[fc] : f.n1[fc];
         for (int r = 0; r < ng; r++)
         {
            rows.push_back({f.vd[f.off_vd[fc] + r], f.off_G[fc] + r * m_, 2 * fc * nbf});
         }
      }
      std::stable_sort(rows.begin(), rows.end(),
                       [](const Row &a, const Row &b) { return a.t < b.t; });
      std::vector<int> tgt, ptr(1, 0);
      for (std::size_t i = 0; i < rows.size(); i++)
      {
         if (i == 0 || rows[i].t != rows[i - 1].t) { tgt.push_back(rows[i].t); ptr.push_back(0); }
         ptr.back() = static_cast<int>(i) + 1;
      }
      const int nt = static_cast<int>(tgt.size()), nr = static_cast<int>(rows.size());
      tgt_.SetSize(nt);
      ptr_.SetSize(nt + 1);
      goff_.SetSize(nr);
      soff_.SetSize(nr);
      int *ht = tgt_.HostWrite(), *hp = ptr_.HostWrite(), *hg = goff_.HostWrite(),
           *hs = soff_.HostWrite();
      for (int i = 0; i < nt; i++) { ht[i] = tgt[i]; }
      for (int i = 0; i <= nt; i++) { hp[i] = ptr[i]; }
      for (int i = 0; i < nr; i++) { hg[i] = rows[i].g; hs[i] = rows[i].s; }
      G_.SetSize(static_cast<int>(f.G.size()));
      real_t *hG = G_.HostWrite();
      for (std::size_t i = 0; i < f.G.size(); i++) { hG[i] = f.G[i]; }
      s_.UseDevice(true);
   }

   /// rhs += G s, s the local kernel slip (host, 2 nbf per local fault face,
   /// AddSlipLoadFast's slip_bc); rhs may live on the device.
   void AddMult(const Vector &s, Vector &rhs) const
   {
      MFEM_VERIFY(s.Size() == ns_, "DeviceSlipLoad: slip size " << s.Size() << " != " << ns_);
      if (tgt_.Size() == 0) { return; }
      s_.SetSize(ns_);
      const real_t *hs = s.HostRead();
      real_t *hd = s_.HostWrite();
      for (int i = 0; i < ns_; i++) { hd[i] = hs[i]; }
      AddMultDevice(s_, rhs);
   }

   /// rhs += G s for a slip s on the device (plan D3: the device fault
   /// step's ToKernelLocal), as AddMult.
   void AddMultDevice(const Vector &s, Vector &rhs) const
   {
      MFEM_VERIFY(s.Size() == ns_, "DeviceSlipLoad: slip size " << s.Size() << " != " << ns_);
      const int nt = tgt_.Size(), m = m_;
      if (nt == 0) { return; }
      const int *t = tgt_.Read(), *p = ptr_.Read(), *go = goff_.Read(), *so = soff_.Read();
      const real_t *G = G_.Read(), *ds = s.Read();
      real_t *r = rhs.ReadWrite();
      mfem::forall(nt, [=] MFEM_HOST_DEVICE (int i)
      {
         real_t acc = r[t[i]];
         for (int j = p[i]; j < p[i + 1]; j++)
         {
            const real_t *g = G + go[j], *sj = ds + so[j];
            real_t v = 0.0;
            for (int c = 0; c < m; c++) { v += g[c] * sj[c]; }
            acc += v;
         }
         r[t[i]] = acc;
      });
   }

   /// Bytes of the device copy (G and the index arrays).
   long long Bytes() const
   {
      return 8LL * G_.Size() + 4LL * (tgt_.Size() + ptr_.Size() + goff_.Size() + soff_.Size());
   }

private:
   int m_ = 0, ns_ = 0;
   Array<int> tgt_, ptr_, goff_, soff_;
   Vector G_;
   mutable Vector s_;
};

/// The messages of DeviceFaultTraction and DeviceFaultMaps: one receive per
/// rank of recv_ranks into recv (values recv_off[k] .. recv_off[k + 1]) and
/// one send per rank of send_ranks from send (send_off[k] ..), then a wait
/// for all.  With GPU-aware MPI both buffers stay in device memory, else the
/// messages go through their host copies; recv is valid where the next
/// Read asks for it either way.
inline void DeviceFaultMessages(MPI_Comm comm, int tag, const std::vector<int> &recv_ranks,
                                const std::vector<int> &recv_off, Vector &recv,
                                const std::vector<int> &send_ranks,
                                const std::vector<int> &send_off, const Vector &send,
                                std::vector<MPI_Request> &reqs)
{
   const int nr = static_cast<int>(recv_ranks.size()), ns = static_cast<int>(send_ranks.size());
   if (nr + ns == 0) { return; }
   const bool gpu_mpi = Device::IsEnabled() && Device::GetGPUAwareMPI();
   const real_t *sp = gpu_mpi ? send.Read() : send.HostRead();
   real_t *rp = gpu_mpi ? recv.Write() : recv.HostWrite();
   if (gpu_mpi) { MFEM_STREAM_SYNC; }
   reqs.resize(nr + ns);
   for (int k = 0; k < nr; k++)
   {
      MPI_Irecv(rp + recv_off[k], recv_off[k + 1] - recv_off[k], MPITypeMap<real_t>::mpi_type,
                recv_ranks[k], tag, comm, &reqs[k]);
   }
   for (int k = 0; k < ns; k++)
   {
      MPI_Isend(sp + send_off[k], send_off[k + 1] - send_off[k], MPITypeMap<real_t>::mpi_type,
                send_ranks[k], tag, comm, &reqs[nr + k]);
   }
   MPI_Waitall(nr + ns, reqs.data(), MPI_STATUSES_IGNORE);
}

/// ComputeTractionFast of an ElasticityDomainOperator on the device (plan
/// D3).  One thread per output value (local fault DOF k, row r = dip,
/// strike, normal) sums its row of A_f against [u_1, u_2, slip] in the
/// host's order; side 2 of a shared face reads the neighbour exchange,
/// packed and unpacked on the device in FaultNbrExchange_'s layout.
class DeviceFaultTraction
{
public:
   template <class Op>
   DeviceFaultTraction(const Op &op, MPI_Comm comm) : comm_(comm)
   {
      const auto &f = op.GetFaultFaceOps();
      MFEM_VERIFY(f.built, "DeviceFaultTraction: call BuildFaultFaceOps first");
      nbf_ = op.GetNbfPerFace();
      nf_ = f.nf;
      nfi_ = f.nfi;
      auto copy = [](const std::vector<int> &h, Array<int> &d)
      {
         d.SetSize(static_cast<int>(h.size()));
         int *p = d.HostWrite();
         for (std::size_t i = 0; i < h.size(); i++) { p[i] = h[i]; }
      };
      copy(f.n1, n1_);
      copy(f.n2, n2_);
      copy(f.off_vd, off_vd_);
      copy(f.off_A, off_A_);
      copy(f.vd, vd_);
      copy(f.nbr_off, nbr_off_);
      copy(f.send_vd, send_vd_);
      A_.SetSize(static_cast<int>(f.A.size()));
      real_t *hA = A_.HostWrite();
      for (std::size_t i = 0; i < f.A.size(); i++) { hA[i] = f.A[i]; }
      // FaultNbrExchange_: send_off counts elements, recv_off values
      const int vdim = f.send_elems.empty() ? 0
                       : static_cast<int>(f.send_vd.size() / f.send_elems.size());
      recv_ranks_ = f.recv_ranks;
      recv_off_ = f.recv_off;
      send_ranks_ = f.send_ranks;
      send_off_.resize(f.send_off.size());
      for (std::size_t k = 0; k < f.send_off.size(); k++) { send_off_[k] = f.send_off[k] * vdim; }
      send_buf_.UseDevice(true);
      recv_buf_.UseDevice(true);
      send_buf_.SetSize(send_vd_.Size());
      recv_buf_.SetSize(static_cast<int>(f.recv_buf.size()));
   }

   /// Tk (2 per local fault DOF: dip, strike) and sn (1 per local DOF,
   /// positive in compression) from the local displacement u and the local
   /// kernel slip s (2 nbf per local face), as ComputeTractionFast(u, s, Tk,
   /// sn).  Collective (the neighbour exchange).
   void Mult(const Vector &u, const Vector &s, Vector &Tk, Vector &sn) const
   {
      const int nd = nf_ * nbf_;
      MFEM_VERIFY(s.Size() == 2 * nd, "DeviceFaultTraction: slip size " << s.Size()
                  << " != " << 2 * nd);
      Tk.SetSize(2 * nd);
      sn.SetSize(nd);
      const real_t *du = u.Read(true);
      Exchange_(du);
      if (nd == 0) { return; }
      const int nbf = nbf_, nfi = nfi_;
      const int *n1 = n1_.Read(true), *n2 = n2_.Read(true), *ov = off_vd_.Read(true),
                 *oa = off_A_.Read(true), *vd = vd_.Read(true), *no = nbr_off_.Read(true);
      const real_t *A = A_.Read(true), *ds = s.Read(true), *rb = recv_buf_.Read(true);
      real_t *dT = Tk.Write(true), *dn = sn.Write(true);
      mfem::forall(3 * nd, [=] MFEM_HOST_DEVICE (int i)
      {
         const int r = i % 3, k = i / 3, f = k / nbf, kk = k % nbf;
         const int m1 = n1[f], m2 = n2[f], nc = m1 + m2 + 2 * nbf;
         const int *v = vd + ov[f];
         const real_t *a = A + oa[f] + static_cast<std::size_t>(r * nbf + kk) * nc;
         real_t t = 0.0;
         for (int c = 0; c < m1; c++) { t += a[c] * du[v[c]]; }
         if (f < nfi)
         {
            for (int c = 0; c < m2; c++) { t += a[m1 + c] * du[v[m1 + c]]; }
         }
         else
         {
            const real_t *nb = rb + no[f - nfi];
            for (int c = 0; c < m2; c++) { t += a[m1 + c] * nb[c]; }
         }
         const real_t *sf = ds + 2 * f * nbf;
         for (int c = 0; c < 2 * nbf; c++) { t += a[m1 + m2 + c] * sf[c]; }
         if (r < 2) { dT[2 * k + r] = t; }
         else { dn[k] = t; }
      });
   }

   /// Bytes of the device copy (A and the index arrays).
   long long Bytes() const
   {
      return 8LL * (A_.Size() + send_buf_.Size() + recv_buf_.Size())
             + 4LL * (n1_.Size() + n2_.Size() + off_vd_.Size() + off_A_.Size() + vd_.Size()
                      + nbr_off_.Size() + send_vd_.Size());
   }

   // The function holding the kernel is public: nvcc refuses an extended
   // __host__ __device__ lambda inside a private or protected member.
   void Exchange_(const real_t *du) const
   {
      const int ns = send_vd_.Size();
      if (recv_ranks_.empty() && send_ranks_.empty()) { return; }
      if (ns > 0)
      {
         const int *sv = send_vd_.Read(true);
         real_t *sb = send_buf_.Write(true);
         mfem::forall(ns, [=] MFEM_HOST_DEVICE (int i) { sb[i] = du[sv[i]]; });
      }
      DeviceFaultMessages(comm_, MPITag::kFaultNbrDevice, recv_ranks_, recv_off_, recv_buf_,
                          send_ranks_, send_off_, send_buf_, reqs_);
   }

private:
   MPI_Comm comm_;
   int nbf_ = 1, nf_ = 0, nfi_ = 0;
   Array<int> n1_, n2_, off_vd_, off_A_, vd_, nbr_off_, send_vd_;
   Vector A_;
   std::vector<int> recv_ranks_, recv_off_, send_ranks_, send_off_;
   mutable Vector send_buf_, recv_buf_;
   mutable std::vector<MPI_Request> reqs_;
};

/// The owned <-> local fault DOF maps of an ElasticityDomainOperator on the
/// device (plan D3).  OwnedToLocal() is RestrictToOwnedFault's index;
/// ToKernelLocal is NitscheWaveOperator::ToKernelLocal_: each local DOF
/// takes its owned value or, on a shared fault face, the ghost copy that
/// ExpandOwnedToLocalFault receives (FaultScatter, canonical -> local
/// permutation), or zero, and then fault_frame::ToKernel with eps_f and the
/// kernel triad of that DOF.
class DeviceFaultMaps
{
public:
   /// Bk: the kernel triads (9 x local DOFs: n, t1, t2), eps: eps_f per
   /// local DOF, as NitscheWaveOperator holds them.
   template <class Op>
   DeviceFaultMaps(const Op &op, const DenseMatrix &Bk, const std::vector<real_t> &eps,
                   MPI_Comm comm) : comm_(comm)
   {
      no_ = op.GetNumOwnedFaultDOFs();
      nd_ = op.GetNumFaultDOFs();
      const int nbf = op.GetNbfPerFace();
      MFEM_VERIFY(Bk.Height() == 9 && Bk.Width() == nd_ && static_cast<int>(eps.size()) == nd_,
                  "DeviceFaultMaps: triads " << Bk.Height() << " x " << Bk.Width()
                  << " and " << eps.size() << " eps_f for " << nd_ << " local fault DOFs");
      const Array<int> &own = op.OwnedFaultDofToLocalDof();
      MFEM_VERIFY(own.Size() == no_, "DeviceFaultMaps: " << own.Size()
                  << " owned -> local entries for " << no_ << " owned fault DOFs");
      own_loc_.SetSize(no_);
      int *ho = own_loc_.HostWrite();
      // src >= 0: owned DOF; -1: zero; <= -2: entry -2 - src of the receive buffer
      std::vector<int> src(nd_, -1);
      for (int j = 0; j < no_; j++) { ho[j] = own[j]; src[own[j]] = j; }
      const Array<int> &perm = op.CanonicalToLocalPerm();
      int rpos = 0, spos = 0;
      std::vector<int> send_own;
      recv_off_.assign(1, 0);
      send_off_.assign(1, 0);
      for (const auto &b : op.SharedFaultCommBlocks())
      {
         if (!b.recv_local_faces.empty())
         {
            for (std::size_t j = 0; j < b.recv_local_faces.size(); j++)
            {
               const int lf = b.recv_local_faces[j];
               for (int ckk = 0; ckk < nbf; ckk++)
               {
                  src[lf * nbf + perm[lf * nbf + ckk]] = -2 - (rpos + static_cast<int>(j) * nbf + ckk);
               }
            }
            rpos += static_cast<int>(b.recv_local_faces.size()) * nbf;
            recv_ranks_.push_back(b.neighbor_rank);
            recv_off_.push_back(3 * rpos);
         }
         if (!b.send_owned_faces.empty())
         {
            for (int of : b.send_owned_faces)
            {
               for (int kk = 0; kk < nbf; kk++) { send_own.push_back(of * nbf + kk); }
            }
            spos = static_cast<int>(send_own.size());
            send_ranks_.push_back(b.neighbor_rank);
            send_off_.push_back(3 * spos);
         }
      }
      src_.SetSize(nd_);
      int *hs = src_.HostWrite();
      for (int k = 0; k < nd_; k++) { hs[k] = src[k]; }
      send_own_.SetSize(static_cast<int>(send_own.size()));
      int *hso = send_own_.HostWrite();
      for (std::size_t i = 0; i < send_own.size(); i++) { hso[i] = send_own[i]; }
      bk_.SetSize(9 * nd_);
      eps_.SetSize(nd_);
      real_t *hb = bk_.HostWrite(), *he = eps_.HostWrite();
      for (int i = 0; i < 9 * nd_; i++) { hb[i] = Bk.GetData()[i]; }
      for (int k = 0; k < nd_; k++) { he[k] = eps[k]; }
      send_buf_.UseDevice(true);
      recv_buf_.UseDevice(true);
      send_buf_.SetSize(3 * spos);
      recv_buf_.SetSize(3 * rpos);
   }

   int NumOwned() const { return no_; }
   int NumLocal() const { return nd_; }
   /// Owned fault DOF -> local fault DOF (RestrictToOwnedFault).
   const Array<int> &OwnedToLocal() const { return own_loc_; }
   /// Kernel triads (n_k, t1_k, t2_k; 9 per local DOF) and eps_f per local DOF.
   const Vector &Triads() const { return bk_; }
   const Vector &EpsF() const { return eps_; }

   /// slip_k (2 per local DOF: dip, strike) of the owned three-vectors d
   /// (3 per owned DOF), as ToKernelLocal_(d, slip_k).  Collective (the
   /// ghost copies).
   void ToKernelLocal(const Vector &d, Vector &slip_k) const
   {
      MFEM_VERIFY(d.Size() == 3 * no_, "DeviceFaultMaps: size " << d.Size() << " != 3 x "
                  << no_ << " owned fault DOFs");
      slip_k.SetSize(2 * nd_);
      const real_t *dd = d.Read(true);
      Scatter_(dd);
      if (nd_ == 0) { return; }
      const int *src = src_.Read(true);
      const real_t *rb = recv_buf_.Read(true), *bk = bk_.Read(true), *ep = eps_.Read(true);
      real_t *sk = slip_k.Write(true);
      mfem::forall(nd_, [=] MFEM_HOST_DEVICE (int k)
      {
         const int q = src[k];
         const real_t zero[3] = {0.0, 0.0, 0.0};
         const real_t *dv = (q >= 0) ? dd + 3 * q : (q == -1) ? zero : rb + 3 * (-2 - q);
         fault_frame::ToKernel(dv, ep[k], bk + 9 * k + 3, bk + 9 * k + 6, sk[2 * k],
                               sk[2 * k + 1]);
      });
   }

   // The function holding the kernel is public (nvcc extended lambdas).
   void Scatter_(const real_t *dd) const
   {
      if (recv_ranks_.empty() && send_ranks_.empty()) { return; }
      const int ns = send_own_.Size();
      if (ns > 0)
      {
         const int *so = send_own_.Read(true);
         real_t *sb = send_buf_.Write(true);
         mfem::forall(ns, [=] MFEM_HOST_DEVICE (int i)
         {
            for (int c = 0; c < 3; c++) { sb[3 * i + c] = dd[3 * so[i] + c]; }
         });
      }
      DeviceFaultMessages(comm_, MPITag::kFaultScatterDevice, recv_ranks_, recv_off_, recv_buf_,
                          send_ranks_, send_off_, send_buf_, reqs_);
   }

private:
   MPI_Comm comm_;
   int no_ = 0, nd_ = 0;
   Array<int> own_loc_, src_, send_own_;
   Vector bk_, eps_;
   std::vector<int> recv_ranks_, recv_off_, send_ranks_, send_off_;
   mutable Vector send_buf_, recv_buf_;
   mutable std::vector<MPI_Request> reqs_;
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_NITSCHE_DEVICE_FAULT_HPP

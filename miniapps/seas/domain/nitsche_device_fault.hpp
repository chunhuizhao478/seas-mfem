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
// Neither needs an assembled matrix.  D3 moves the fault itself to the
// device and removes both copies.

#ifndef MFEM_SEAS_NITSCHE_DEVICE_FAULT_HPP
#define MFEM_SEAS_NITSCHE_DEVICE_FAULT_HPP

#include "mfem.hpp"

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
      const int nt = tgt_.Size(), m = m_;
      if (nt == 0) { return; }
      s_.SetSize(ns_);
      const real_t *hs = s.HostRead();
      real_t *hd = s_.HostWrite();
      for (int i = 0; i < ns_; i++) { hd[i] = hs[i]; }
      const int *t = tgt_.Read(), *p = ptr_.Read(), *go = goff_.Read(), *so = soff_.Read();
      const real_t *G = G_.Read(), *ds = s_.Read();
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

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_NITSCHE_DEVICE_FAULT_HPP

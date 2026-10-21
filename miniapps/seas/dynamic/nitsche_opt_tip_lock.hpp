// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// nitsche_opt_tip_lock.hpp — option O5 of the complex-geometries plan
// (rev 2, Section 5.6; Chunhui Zhao 2026-10-20): locked fault tips, for
// TPV22/23 ("border nodes do not slip", except along the free-surface trace).
//
// Deck (Nitsche driver only; absent, or enabled = false, = every bit as before):
//
//   [fault_tip_lock]
//   enabled       = true
//   exclude_attrs = [101]   # optional; boundary attributes whose edges stay
//                           # free (default: [boundary] natural_attrs)
//
// A tip edge is an edge of exactly one fault face (boundary element with the
// fault attribute) that is not an edge of a boundary element with an
// excluded attribute.  The edges come from the SERIAL mesh before the
// partition, so the set does not depend on the rank count, and they are
// needed before the operator exists: NitscheWaveOperator reads
// NitscheFaultSetup::parameters in its constructor.  WrapTipLock wraps that
// function and locks every fault DOF whose position lies on a tip edge
// (within 1e-6 of the edge length), so the DG copies of every face that
// touches a tip, even at one vertex, are locked; the vertex where a side tip
// meets the trace lies on the side tip and is locked too.
//
// Locked = mu_s = kLockedMuS and cohesion = kLockedCohesion: the strength
// exceeds any traction, so V = 0 and the slip stays 0 (the slip-weakening
// closure has no barrier sentinel; a large mu_s never yields).  The face
// still transmits traction like a welded interface; sigma_n and the step
// size are untouched.  Linear slip weakening only: the driver refuses the
// option with rate and state.  The driver also refuses an attribute listed in
// exclude_attrs that shares no edge with the fault (a misspelt attribute would
// lock the trace it was meant to free) and a lock that locks no DOF.

#ifndef MFEM_SEAS_NITSCHE_OPT_TIP_LOCK_HPP
#define MFEM_SEAS_NITSCHE_OPT_TIP_LOCK_HPP

#include "mfem.hpp"
#include "nitsche_opt_common.hpp"
#include "../friction/fault_closure.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace mfem
{
namespace seas
{
namespace nitsche_opt
{

struct TipLockConfig
{
   static constexpr const char *kTable = "fault_tip_lock";
   static constexpr real_t kLockedMuS = 1.0e4;
   static constexpr real_t kLockedCohesion = 1.0e12;   // Pa
   bool enabled = false;
   bool has_exclude = false;          ///< exclude_attrs given (else the natural attrs)
   std::vector<int> exclude_attrs;
};

inline TipLockConfig ParseTipLock(const toml::value &t)
{
   using T = TipLockConfig;
   OptRefuseUnknownKeys(t, T::kTable, {"enabled", "exclude_attrs"});
   T c;
   c.enabled = OptBool(t, T::kTable, "enabled");
   if (t.contains("exclude_attrs"))
   {
      c.has_exclude = true;
      c.exclude_attrs = OptIntArray(t, T::kTable, "exclude_attrs");
      for (int a : c.exclude_attrs)
      {
         MFEM_VERIFY(a > 0, "[" << T::kTable << "] exclude_attrs: attribute " << a << " <= 0");
      }
   }
   return c;
}

/// The tip edges of a fault surface as segments, and a point query.
class TipEdges
{
public:
   /// Tip edges of the fault faces with attribute `fault_attr` of the serial
   /// mesh `mesh`; edges of boundary elements with an attribute in
   /// `exclude` are not tips.  Aborts if a fault face appears twice.
   TipEdges(const Mesh &mesh, int fault_attr, const std::set<int> &exclude)
   {
      std::map<std::pair<int, int>, int> count;
      std::set<std::pair<int, int>> excl;
      std::set<std::pair<std::pair<int, int>, int>> excl_by_attr;   // (edge, attribute)
      std::set<std::array<int, 3>> faces;
      Array<int> v;
      auto edge = [](int a, int b) { return std::make_pair(std::min(a, b), std::max(a, b)); };
      for (int be = 0; be < mesh.GetNBE(); be++)
      {
         const int attr = mesh.GetBdrAttribute(be);
         const bool is_fault = (attr == fault_attr), is_excl = exclude.count(attr) > 0;
         if (!is_fault && !is_excl) { continue; }
         mesh.GetBdrElementVertices(be, v);
         MFEM_VERIFY(v.Size() == 3, "[fault_tip_lock]: boundary element " << be
                     << " is not a triangle (tetrahedral meshes only)");
         if (is_fault)
         {
            std::array<int, 3> f = {v[0], v[1], v[2]};
            std::sort(f.begin(), f.end());
            MFEM_VERIFY(faces.insert(f).second, "[fault_tip_lock]: fault face with vertices "
                        << f[0] << " " << f[1] << " " << f[2] << " appears twice among the "
                        "boundary elements; its border edges cannot be found");
            nfaces_++;
         }
         for (int k = 0; k < 3; k++)
         {
            const auto e = edge(v[k], v[(k + 1) % 3]);
            if (is_fault) { count[e]++; }
            if (is_excl) { excl.insert(e); excl_by_attr.insert({e, attr}); }
         }
      }
      for (const auto &ea : excl_by_attr)
      {
         if (count.count(ea.first)) { shared_[ea.second]++; }
      }
      for (const auto &kv : count)
      {
         if (kv.second != 1 || excl.count(kv.first)) { continue; }
         const real_t *a = mesh.GetVertex(kv.first.first), *b = mesh.GetVertex(kv.first.second);
         seg_.push_back({a[0], a[1], a[2], b[0], b[1], b[2]});
      }
      BuildGrid_();
   }

   int NumFaultFaces() const { return nfaces_; }
   int NumEdges() const { return static_cast<int>(seg_.size()); }
   /// Fault-face edges that are also edges of boundary elements with the
   /// excluded attribute `attr` (0 when they share none, e.g. a misspelt
   /// attribute: then the trace it was meant to free is a tip).
   int SharedEdges(int attr) const
   {
      const auto it = shared_.find(attr);
      return it == shared_.end() ? 0 : it->second;
   }
   const std::vector<std::array<real_t, 6>> &Segments() const { return seg_; }

   /// Distance from x to segment i.
   real_t Distance(int i, const real_t x[3]) const
   {
      const auto &s = seg_[i];
      real_t ab[3], ax[3], ab2 = 0.0, t = 0.0;
      for (int c = 0; c < 3; c++)
      {
         ab[c] = s[3 + c] - s[c];
         ax[c] = x[c] - s[c];
         ab2 += ab[c] * ab[c];
         t += ab[c] * ax[c];
      }
      t = (ab2 > 0.0) ? std::min(std::max(t / ab2, real_t(0.0)), real_t(1.0)) : 0.0;
      real_t d2 = 0.0;
      for (int c = 0; c < 3; c++) { const real_t r = ax[c] - t * ab[c]; d2 += r * r; }
      return std::sqrt(d2);
   }

   /// True when x lies on a tip edge, within kRelTol of that edge's length.
   bool OnTip(const real_t x[3]) const
   {
      if (seg_.empty()) { return false; }
      const auto it = grid_.find(Key_(x));
      if (it == grid_.end()) { return false; }
      for (int i : it->second)
      {
         if (Distance(i, x) <= kRelTol * len_[i]) { return true; }
      }
      return false;
   }

   static constexpr real_t kRelTol = 1e-6;

private:
   std::vector<std::array<real_t, 6>> seg_;
   std::vector<real_t> len_;
   std::map<int, int> shared_;   // excluded attribute -> edges shared with the fault
   int nfaces_ = 0;
   real_t cell_ = 1.0;
   std::unordered_map<std::int64_t, std::vector<int>> grid_;

   std::int64_t Cell_(real_t x) const { return static_cast<std::int64_t>(std::floor(x / cell_)); }
   static std::int64_t Pack_(std::int64_t i, std::int64_t j, std::int64_t k)
   {
      const std::int64_t o = std::int64_t(1) << 20;
      return ((i + o) << 42) | ((j + o) << 21) | (k + o);
   }
   std::int64_t Key_(const real_t x[3]) const { return Pack_(Cell_(x[0]), Cell_(x[1]), Cell_(x[2])); }

   /// Uniform grid of cell size = the longest tip edge; each segment goes into
   /// every cell its tolerance-padded box touches.
   void BuildGrid_()
   {
      len_.resize(seg_.size());
      real_t lmax = 0.0;
      for (size_t i = 0; i < seg_.size(); i++)
      {
         const auto &s = seg_[i];
         len_[i] = std::sqrt((s[3] - s[0]) * (s[3] - s[0]) + (s[4] - s[1]) * (s[4] - s[1])
                             + (s[5] - s[2]) * (s[5] - s[2]));
         lmax = std::max(lmax, len_[i]);
      }
      if (seg_.empty()) { return; }
      cell_ = lmax > 0.0 ? lmax : 1.0;
      for (size_t i = 0; i < seg_.size(); i++)
      {
         const auto &s = seg_[i];
         const real_t pad = 2.0 * kRelTol * len_[i];
         std::int64_t lo[3], hi[3];
         for (int c = 0; c < 3; c++)
         {
            lo[c] = Cell_(std::min(s[c], s[3 + c]) - pad);
            hi[c] = Cell_(std::max(s[c], s[3 + c]) + pad);
         }
         for (std::int64_t a = lo[0]; a <= hi[0]; a++)
            for (std::int64_t b = lo[1]; b <= hi[1]; b++)
               for (std::int64_t c = lo[2]; c <= hi[2]; c++)
               {
                  grid_[Pack_(a, b, c)].push_back(static_cast<int>(i));
               }
      }
   }
};

/// `params` with every DOF on a tip edge locked (header).
inline std::function<void(const real_t x[3], FaultPointState &st)>
WrapTipLock(std::function<void(const real_t x[3], FaultPointState &st)> params,
            std::shared_ptr<const TipEdges> tips)
{
   return [params, tips](const real_t x[3], FaultPointState &st) {
      params(x, st);
      if (tips->OnTip(x))
      {
         st.mu_s = TipLockConfig::kLockedMuS;
         st.cohesion = TipLockConfig::kLockedCohesion;
      }
   };
}

} // namespace nitsche_opt
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_NITSCHE_OPT_TIP_LOCK_HPP

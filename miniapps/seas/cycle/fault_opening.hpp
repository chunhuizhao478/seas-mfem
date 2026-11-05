// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// fault_opening.hpp — [cycle.opening] mode = "shear_free": version A of the
// fault opening design (rev 1.1, 2026-11-04,
// /mnt/project-files/fully-dynamic-pml/fault-opening-design/; Chunhui Zhao,
// 2026-11-04: "start version A", "build and test here"), the QD regime's
// part.  The dynamic regime's open rule is in NitscheWaveOperator::CloseDOF.
//
// The rule (decisions D2-D5 of the design, as recommended):
//   * a fault DOF is OPEN while its true elastic normal stress
//     sigma = sigma0 + T_n (compression positive; no relaxation r, no floor)
//     is below 0, and CLOSED otherwise;
//   * an open DOF carries no shear traction, tau = 0, in both regimes; the
//     faces stay glued in the normal direction (sigma stays tensile);
//   * psi and the relaxation r are frozen while open; the floor acts on
//     closed DOFs only; max V, the switch thresholds and the exit-7 guard
//     are over the closed DOFs (the kernel gives open DOFs V = 0).
//
// QD open-set solve (design Section 3.2, eq. 3.5), at every stage, on the
// Green's-function path T = G delta + t T_D (qd_green_operator.hpp):
//
//    O        the open unknowns (fault DOF, in-plane kernel slot): 2-D only
//    s_O      the state's slip there: the integrator holds it (rate 0), and
//             QDPhase writes the solved slip back after every accepted step
//    tau_O    tau_pre + T at O from the stage, computed with s_O
//    G_OO     the |O| x |O| block of G (shear rows of O, slip columns of O)
//
//    G_OO Delta = -tau_O,   delta_O = s_O + Delta
//    T += G_{:,O} Delta,    T_n += G^n_{:,O} Delta
//
// so after the stage tau = 0 at O to round-off and the closed DOFs see the
// Schur complement S = G_CC - G_CO G_OO^-1 G_OC of eq. 3.6 without forming
// it.  The columns G_{:,O} (owned rows of T and T_n) are unit products of
// the Green's function (t = 0: G e exactly), fetched once per DOF the first
// time it opens and kept.  G_OO is gathered on every rank and factored
// (LU, partial pivoting) on every rank alike, so every rank solves the same
// Delta bit for bit.  Cost per stage: one gather of |O| values, one |O| x
// |O| solve, |O| x 3 N_own multiply-adds.
//
// Set changes (Update) are made after accepted steps, at the switch back
// and at the start (QDPhase::SettleOpenSet_): a newly open DOF drops its
// shear at once (a static slip step), a closing DOF keeps the slip it has.
//
// No solver header (rule M4).

#ifndef MFEM_SEAS_CYCLE_FAULT_OPENING_HPP
#define MFEM_SEAS_CYCLE_FAULT_OPENING_HPP

#include "mfem.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <vector>

namespace mfem
{
namespace seas
{
namespace cycle
{

class OpenSetQD
{
public:
   /// e (2 N_own, kernel slip pairs on owned DOFs) -> T (2 N_own), T_n (N_own):
   /// the elastic traction of the slip e with no loading (G e).  Collective.
   using UnitProduct = std::function<void(const Vector &e, Vector &T, Vector &Tn)>;

   /// @a nown owned fault DOFs (the Green's function's owned order); @a slot
   /// the in-plane kernel slot (0 dip, 1 strike); @a tau_pre the shear
   /// pre-stress the kernel adds (2 N_own, held by reference).  Collective.
   OpenSetQD(MPI_Comm comm, int nown, int slot, UnitProduct unit, const Vector &tau_pre)
      : comm_(comm), nown_(nown), slot_(slot), unit_(std::move(unit)), tau_pre_(tau_pre)
   {
      MFEM_VERIFY(slot == 0 || slot == 1, "OpenSetQD: slot " << slot);
      MFEM_VERIFY(tau_pre.Size() == 2 * nown, "OpenSetQD: " << tau_pre.Size()
                  << " pre-stress values for " << nown << " DOFs");
      MPI_Comm_rank(comm_, &rank_);
      MPI_Comm_size(comm_, &nranks_);
      // The global fault index follows the partition, rank by rank in owned
      // order (qd_green_operator.hpp): the first is the sum of the lower
      // ranks' counts.
      long long mine = nown, first = 0;
      MPI_Exscan(&mine, &first, 1, MPI_LONG_LONG, MPI_SUM, comm_);
      first_ = (rank_ == 0) ? 0 : first;
      open_.assign(nown_, 0);
      mask_.SetSize(nown_);
      mask_ = 0.0;
      counts_.assign(nranks_, 0);
      displs_.assign(nranks_, 0);
   }

   int NumOwned() const { return nown_; }
   int Slot() const { return slot_; }
   /// Per owned DOF: 1 open, 0 closed.
   const std::vector<char> &Flags() const { return open_; }
   bool IsOpen(int j) const { return open_[j] != 0; }
   /// The kernel's and the relaxation's mask (1.0 open, 0.0 closed).
   const Vector &Mask() const { return mask_; }
   /// Open DOFs over all ranks.
   int NumOpenGlobal() const { return static_cast<int>(list_.size()); }
   /// Unit products and set changes so far (this rank's view, the same on all).
   long UnitProducts() const { return unit_products_; }
   long SetChanges() const { return set_changes_; }
   long Flips() const { return flips_; }

   /// Replace the flags (the switch back, a restart).  Collective.
   void SetFlags(const std::vector<char> &f)
   {
      MFEM_VERIFY(static_cast<int>(f.size()) == nown_, "OpenSetQD::SetFlags: " << f.size()
                  << " flags for " << nown_ << " DOFs");
      long n = 0;
      for (int j = 0; j < nown_; j++)
      {
         const char v = f[j] ? 1 : 0;
         if (v != open_[j]) { n++; }
         open_[j] = v;
      }
      MPI_Allreduce(MPI_IN_PLACE, &n, 1, MPI_LONG, MPI_SUM, comm_);
      if (n > 0) { set_changes_++; flips_ += n; }
      Rebuild_();
   }

   /// The rule at one state: closed with sigma < 0 opens, open with sigma >=
   /// 0 closes (sigma: true elastic normal stress per owned DOF, compression
   /// positive).  Returns the global number of DOFs that changed; with any,
   /// the set is rebuilt.  Collective.
   long Update(const Vector &sigma)
   {
      MFEM_VERIFY(sigma.Size() == nown_, "OpenSetQD::Update: " << sigma.Size()
                  << " values for " << nown_ << " DOFs");
      long n = 0;
      last_flipped_.clear();
      for (int j = 0; j < nown_; j++)
      {
         MFEM_VERIFY(std::isfinite(sigma(j)), "OpenSetQD::Update: sigma_n' at owned DOF "
                     << j << " on rank " << rank_ << " is " << sigma(j));
         const char want = (sigma(j) < 0.0) ? 1 : 0;
         if (want != open_[j])
         {
            open_[j] = want;
            last_flipped_.push_back(j);
            n++;
         }
      }
      MPI_Allreduce(MPI_IN_PLACE, &n, 1, MPI_LONG, MPI_SUM, comm_);
      if (n > 0)
      {
         set_changes_++;
         flips_ += n;
         Rebuild_();
      }
      return n;
   }

   /// The stage hook: @a slip (2 N_own), @a T (2 N_own), @a Tn (N_own) of a
   /// Green's-function stage, computed with the state's slip at the open
   /// unknowns.  Afterwards tau_pre + T = 0 there (round-off) and slip, T,
   /// T_n are those of the solved open slip (file comment).  Collective.
   void Stage(Vector &slip, Vector &T, Vector *Tn)
   {
      const int m = static_cast<int>(list_.size());
      if (m == 0) { return; }
      MFEM_VERIFY(Tn != nullptr, "OpenSetQD::Stage: the stage has no normal traction "
                  "(elastic sigma_n is required)");
      const int mo = static_cast<int>(own_.size());
      const real_t *tp = tau_pre_.HostRead();
      real_t *s = slip.HostReadWrite();
      real_t *t = T.HostReadWrite();
      real_t *tn = Tn->HostReadWrite();
      r_own_.resize(std::max(mo, 1));
      for (int k = 0; k < mo; k++)
      {
         const int i = 2 * own_[k] + slot_;
         r_own_[k] = tp[i] + t[i];
      }
      r_.SetSize(m);
      MPI_Allgatherv(r_own_.data(), mo, MPITypeMap<real_t>::mpi_type, r_.HostWrite(),
                     counts_.data(), displs_.data(), MPITypeMap<real_t>::mpi_type, comm_);
      delta_.SetSize(m);
      inv_->Mult(r_, delta_);
      delta_.Neg();
      for (int k = 0; k < mo; k++) { s[2 * own_[k] + slot_] += delta_(own_pos_[k]); }
      const int n2 = 2 * nown_;
      for (int q = 0; q < m; q++)
      {
         const real_t dq = delta_(q);
         if (dq == 0.0) { continue; }
         const real_t *c = cols_.at(list_[q]).data();
         for (int i = 0; i < n2; i++) { t[i] += dq * c[i]; }
         for (int i = 0; i < nown_; i++) { tn[i] += dq * c[n2 + i]; }
      }
      stages_++;
   }

   /// The last stage's |Delta| max (the static step size of the open set).
   real_t LastStepMax() const { return delta_.Size() ? delta_.Normlinf() : 0.0; }
   /// The owned DOFs the last Update flipped (this rank's; their new state
   /// is IsOpen).  For the settle log.
   const std::vector<int> &LastFlipped() const { return last_flipped_; }
   /// Global index of owned DOF @a j (the open list's numbering).
   long long GlobalIndex(int j) const { return first_ + j; }

private:
   /// The global open list, its columns and the factored G_OO.  Collective.
   void Rebuild_()
   {
      delta_.SetSize(0);   // LastStepMax: no stale solve after a change
      real_t *mk = mask_.HostWrite();   // the kernel reads it with Read(use_dev)
      for (int j = 0; j < nown_; j++) { mk[j] = open_[j] ? 1.0 : 0.0; }
      // own_: the owned open DOFs in owned order (= global order)
      own_.clear();
      for (int j = 0; j < nown_; j++) { if (open_[j]) { own_.push_back(j); } }
      const int mo = static_cast<int>(own_.size());
      MPI_Allgather(&mo, 1, MPI_INT, counts_.data(), 1, MPI_INT, comm_);
      int m = 0;
      for (int r = 0; r < nranks_; r++)
      {
         displs_[r] = m;
         m += counts_[r];
      }
      std::vector<long long> mine(std::max(mo, 1));
      for (int k = 0; k < mo; k++) { mine[k] = first_ + own_[k]; }
      list_.assign(m, 0);
      MPI_Allgatherv(mine.data(), mo, MPI_LONG_LONG, list_.data(), counts_.data(),
                     displs_.data(), MPI_LONG_LONG, comm_);
      own_pos_.resize(mo);
      for (int k = 0; k < mo; k++) { own_pos_[k] = displs_[rank_] + k; }
      inv_.reset();
      if (m == 0) { return; }

      // Columns not fetched yet: one unit product each (the same list on
      // every rank, so the collective products match).
      Vector e(2 * nown_), T, Tn;
      for (int q = 0; q < m; q++)
      {
         const long long g = list_[q];
         if (cols_.count(g)) { continue; }
         e = 0.0;
         if (g >= first_ && g < first_ + nown_)
         {
            e(2 * static_cast<int>(g - first_) + slot_) = 1.0;
         }
         unit_(e, T, Tn);
         unit_products_++;
         MFEM_VERIFY(T.Size() == 2 * nown_ && Tn.Size() == nown_, "OpenSetQD: the unit "
                     "product gave " << T.Size() << " shear and " << Tn.Size()
                     << " normal values for " << nown_ << " DOFs");
         std::vector<real_t> &c = cols_[g];
         c.resize(3 * static_cast<size_t>(nown_));
         T.HostRead();
         Tn.HostRead();
         for (int i = 0; i < 2 * nown_; i++) { c[i] = T(i); }
         for (int i = 0; i < nown_; i++) { c[2 * nown_ + i] = Tn(i); }
      }

      // G_OO: this rank's rows (its open DOFs' in-plane shear), gathered
      std::vector<real_t> rows(std::max(static_cast<size_t>(mo) * m, size_t(1)));
      for (int k = 0; k < mo; k++)
      {
         const int i = 2 * own_[k] + slot_;
         for (int q = 0; q < m; q++) { rows[static_cast<size_t>(k) * m + q] = cols_.at(list_[q])[i]; }
      }
      std::vector<int> rc(nranks_), rd(nranks_);
      for (int r = 0; r < nranks_; r++)
      {
         rc[r] = counts_[r] * m;
         rd[r] = displs_[r] * m;
      }
      std::vector<real_t> all(static_cast<size_t>(m) * m);
      MPI_Allgatherv(rows.data(), mo * m, MPITypeMap<real_t>::mpi_type, all.data(),
                     rc.data(), rd.data(), MPITypeMap<real_t>::mpi_type, comm_);
      A_.SetSize(m);
      for (int p = 0; p < m; p++)
      {
         for (int q = 0; q < m; q++) { A_(p, q) = all[static_cast<size_t>(p) * m + q]; }
      }
      // Each diagonal entry is the self-stiffness of an open DOF; zero means
      // a column with no response (a DOF outside the Green's function's slip
      // space, e.g. the out-of-plane slot).
      for (int p = 0; p < m; p++)
      {
         MFEM_VERIFY(std::isfinite(A_(p, p)) && A_(p, p) != 0.0, "OpenSetQD: G_OO("
                     << p << ", " << p << ") = " << A_(p, p) << " for global fault DOF "
                     << list_[p]);
      }
      inv_ = std::make_unique<DenseMatrixInverse>(A_);
   }

   MPI_Comm comm_;
   int rank_ = 0, nranks_ = 1;
   int nown_ = 0, slot_ = 0;
   long long first_ = 0;
   UnitProduct unit_;
   const Vector &tau_pre_;
   std::vector<char> open_;
   Vector mask_;
   std::vector<int> own_;          ///< owned open DOFs (owned index)
   std::vector<int> own_pos_;      ///< their positions in list_
   std::vector<long long> list_;   ///< global open DOFs, ascending
   std::vector<int> counts_, displs_;
   std::map<long long, std::vector<real_t>> cols_;   ///< [T (2 N_own) | T_n (N_own)]
   DenseMatrix A_;
   std::unique_ptr<DenseMatrixInverse> inv_;
   std::vector<real_t> r_own_;
   std::vector<int> last_flipped_;   ///< owned DOFs flipped by the last Update
   Vector r_, delta_;
   long unit_products_ = 0, set_changes_ = 0, flips_ = 0, stages_ = 0;
};

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_FAULT_OPENING_HPP

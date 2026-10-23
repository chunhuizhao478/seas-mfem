// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// bench_output.hpp — fully dynamic BP5 cycle: the outputs of the SCEC BP5
// description (section 4) besides the on-fault station files, on with
// [cycle.benchmark] (cycle_config.hpp).  Files, all under <dir>/<prefix>:
//   _source.txt              t, log10 max |V| over the domain A, moment rate
//                            over A [N m/s]; a row with every station row
//   _<name>.txt              one per [[cycle.benchmark.volume_station]]: t,
//                            u.e1, u.e2, u.e3 [m], v.e1, v.e2, v.e3 [m/s]
//                            (e_i: volume_axes); a row with every station row
//   _profile_<name>_<q>.txt  one per [[cycle.benchmark.profile]] and q in
//                            slip_2, slip_3, stress_2, stress_3: a first row
//                            "0 0 s_1 .. s_n" (s: the coordinate the profile
//                            lists, x2 or x3 [m]), then rows t, log10 max V
//                            (whole fault), q at each point; a row every
//                            profile_dt_slow, profile_dt_fast while max V >
//                            fast_V ([cycle.output])
//   _events.txt              a row per event: event, t_start, t_end, rupture
//                            area R [m^2], tau_avg at t_start and at t_end
//                            [MPa], s_avg [m]
//   _rupture_event<k>.txt    at the end of event k, a row per fault DOF: x,
//                            y, z, x2, x3, t_rupture (1e9: never)
// Plane coordinates (x2, x3) are those of [cycle.output]: (x - plane_origin)
// . plane_axis1 and . plane_axis2.  The SCEC upload files (headers, units,
// the time steps of global.dat) are written from these by
// cycle/scec_bp5_export.py.
//
// Definitions (BP5 section 4; components and shear stress are the station
// files', term by term, so the QD regime matches seas_bp5_full's sign):
//   |V| = sqrt(V_strike^2 + V_dip^2) at a fault DOF; max |V| over A: over
//   the DOFs in A.
//   moment rate = int_A mu |V| dA, with the fault face quadrature of the
//   traction (2p + 1 rule): sum over faces and points x_q in A of
//   w_q |J_q| mu |sum_k phi_k(x_q) V_k|.
//   An event is FdCycleOutput's (it starts when max V > event_V and ends
//   event_dwell s below it; t_end is the first time below).  Its rupture
//   area Omega_R holds the DOFs whose |V| reached rupture_V at a step of the
//   event; with a_k = int phi_k dA (zero at the vertices for p = 2: the
//   integral of the interpolant),
//      R          = sum_{k in Omega_R} a_k
//      tau_avg(t) = | sum_k a_k tau_k(t) | / R     (tau = (tau_strike, tau_dip))
//      s_avg      = | sum_k a_k (s_k(t_end) - s_k(t_start)) | / R
//   (0 when R = 0).  t_rupture of a DOF: the first step of the event with
//   |V| >= rupture_V.  The driver calls Track at every step of both regimes
//   (every fault step under local time stepping), so the rupture times have
//   that resolution.
// Off-fault values: the displacement and velocity fields of the regime at
// each station point (ParMesh::FindPoints; the lowest rank among those that
// find it writes the file).  In the QD regime the velocity is the exact
// quasi-static du/dt (QDPhase::Velocity), in the dynamic regime U'.
//
// Restart: GetState holds the profile clock, the file sizes and the open
// event's per-DOF arrays (checkpoint block FD_CYCLE_V4); a writer built with
// a State truncates the files to those sizes and continues.  Same rank count
// as the run that wrote it.
//
// No solver header (rule M4): fields arrive through the regime adapters.

#ifndef MFEM_SEAS_CYCLE_BENCH_OUTPUT_HPP
#define MFEM_SEAS_CYCLE_BENCH_OUTPUT_HPP

#include "mfem.hpp"
#include "cycle_config.hpp"
#include "fd_cycle_output.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <functional>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{
namespace cycle
{

/// The fault data the benchmark output needs, on this rank's owned fault
/// DOFs (QDPhase::BenchGeometry fills it).  The owned DOFs are grouped by
/// face, nbf per face.
struct BenchFaultGeometry
{
   Vector X;            ///< positions, 3 N
   Vector x2, x3;       ///< plane coordinates, N
   int nbf = 1;         ///< DOFs per fault face
   DenseMatrix e_q;     ///< face basis at the quadrature points, nbf x nq
   Vector wJ;           ///< w_q |J_q| per owned face and point (faces x nq)
   Vector q2, q3;       ///< plane coordinates of those points
   /// Every profile point in deck order (profile by profile).
   std::vector<StationLocation> profile_locs;
};

class BenchmarkOutput
{
public:
   using FieldsFn = FdCycleOutput::FieldsFn;
   using VolumeFn = VolumeFieldsFn;
   /// |V| on the owned fault DOFs.
   using SpeedFn = std::function<void(Vector &)>;

   static constexpr real_t kNever = 1.0e9;
   static constexpr int kVolumeColumns = 7;
   static constexpr int kEventColumns = 7;
   static const char *const kProfileQuantities[4];

   struct State
   {
      real_t profile_last = -1e30;
      std::vector<std::string> files;   ///< this rank's files, by name
      std::vector<long long> sizes;
      // The open event on this rank's owned DOFs (empty outside an event).
      bool open = false;
      bool have_end = false;
      Vector trupt, slip0, tau0, slip_end, tau_end;   ///< N, 2N, 2N, 2N, 2N
   };

   /// `fes`: the displacement space of both regimes; `geom`: the owned fault
   /// DOFs.  With `restart`, the files continue from the checkpointed sizes.
   /// Collective.
   BenchmarkOutput(const CycleConfig &cfg, MPI_Comm comm, ParMesh &pmesh,
                   const ParFiniteElementSpace &fes, const BenchFaultGeometry &geom,
                   const std::string &commit, const State *restart = nullptr)
      : cfg_(cfg.bench), out_(cfg.output), comm_(comm), fes_(fes), geom_(geom),
        mu_(cfg.shared.material_fallback.mu)
   {
      MPI_Comm_rank(comm_, &rank_);
      MFEM_VERIFY(cfg_.enabled, "BenchmarkOutput: [cycle.benchmark] is off");
      rupture_V_ = cfg_.rupture_V > 0.0 ? cfg_.rupture_V : out_.event_V;
      const int N = geom_.x2.Size();
      MFEM_VERIFY(geom_.x3.Size() == N && geom_.X.Size() == 3 * N && N % geom_.nbf == 0,
                  "BenchmarkOutput: fault geometry of " << N << " DOFs is inconsistent");
      const int nq = geom_.e_q.Width();
      const int nf = N / geom_.nbf;
      MFEM_VERIFY(geom_.e_q.Height() == geom_.nbf && geom_.wJ.Size() == nf * nq &&
                  geom_.q2.Size() == nf * nq && geom_.q3.Size() == nf * nq,
                  "BenchmarkOutput: face quadrature of " << nf << " faces is inconsistent");

      // Area weight a_k of each DOF, the DOFs in A and the points in A.
      area_.SetSize(N);
      area_ = 0.0;
      inA_.assign(N, 0);
      wA_.SetSize(nf * nq);
      for (int f = 0; f < nf; f++)
      {
         for (int q = 0; q < nq; q++)
         {
            const real_t w = geom_.wJ(f * nq + q);
            for (int k = 0; k < geom_.nbf; k++) { area_(f * geom_.nbf + k) += w * geom_.e_q(k, q); }
            wA_(f * nq + q) = InA_(geom_.q2(f * nq + q), geom_.q3(f * nq + q)) ? w : 0.0;
         }
      }
      for (int k = 0; k < N; k++) { inA_[k] = InA_(geom_.x2(k), geom_.x3(k)) ? 1 : 0; }
      {
         real_t loc[3] = {0.0, 0.0, 0.0}, glob[3];
         for (int k = 0; k < N; k++) { loc[0] += area_(k); loc[1] += inA_[k]; }
         for (int i = 0; i < wA_.Size(); i++) { loc[2] += wA_(i); }
         MPI_Allreduce(loc, glob, 3, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm_);
         fault_area_ = glob[0];
         area_A_ = glob[2];
         MFEM_VERIFY(glob[1] > 0.0 && glob[2] > 0.0, "BenchmarkOutput: no fault DOF "
                     "or quadrature point lies in the domain A (region_x2, region_x3)");
      }

      if (restart) { state_ = *restart; }
      const State *rs = restart;
      Open_(rs, commit);
      if (restart && state_.open)
      {
         MFEM_VERIFY(state_.trupt.Size() == N && state_.slip0.Size() == 2 * N &&
                     state_.tau0.Size() == 2 * N, "BenchmarkOutput restart: the open "
                     "event's arrays do not fit this rank's " << N << " fault DOFs");
      }
   }

   real_t FaultArea() const { return fault_area_; }
   real_t AreaA() const { return area_A_; }
   int NumVolumeStations() const { return static_cast<int>(cfg_.volume_stations.size()); }
   /// Rank that writes volume station s.
   int VolumeOwner(int s) const { return vol_owner_[s]; }

   /// The rows that go with a station row at t: the source row and the
   /// off-fault rows.  Collective: `fill` and `volume` run on every rank.
   void Rows(real_t t, const FieldsFn &fill, const VolumeFn &volume)
   {
      StationFields f;
      fill(f);
      const int N = geom_.x2.Size();
      MFEM_VERIFY(f.V_dip.Size() == N, "BenchmarkOutput: the regime filled "
                  << f.V_dip.Size() << " values for " << N << " fault DOFs");
      real_t vmax = 0.0, mdot = 0.0;
      for (int k = 0; k < N; k++)
      {
         if (inA_[k]) { vmax = std::max(vmax, Speed_(f, k)); }
      }
      const int nq = geom_.e_q.Width(), nbf = geom_.nbf;
      for (int fc = 0; fc < N / nbf; fc++)
      {
         for (int q = 0; q < nq; q++)
         {
            const real_t w = wA_(fc * nq + q);
            if (w == 0.0) { continue; }
            real_t vs = 0.0, vd = 0.0;
            for (int k = 0; k < nbf; k++)
            {
               vs += geom_.e_q(k, q) * f.V_strike(fc * nbf + k);
               vd += geom_.e_q(k, q) * f.V_dip(fc * nbf + k);
            }
            mdot += w * mu_ * std::sqrt(vs * vs + vd * vd);
         }
      }
      real_t gmax = 0.0, gsum = 0.0;
      MPI_Reduce(&vmax, &gmax, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, 0, comm_);
      MPI_Reduce(&mdot, &gsum, 1, MPITypeMap<real_t>::mpi_type, MPI_SUM, 0, comm_);
      if (source_)
      {
         source_->Row({t, gmax > 0.0 ? std::log10(gmax) : -300.0, gsum});
         source_->Flush();
      }

      if (!cfg_.volume_stations.empty())
      {
         const Vector *u = nullptr, *v = nullptr;
         volume(u, v);
         MFEM_VERIFY(u && v && u->Size() == fes_.GetVSize() && v->Size() == fes_.GetVSize(),
                     "BenchmarkOutput: the regime's displacement or velocity is not a "
                     "vector of the displacement space (" << fes_.GetVSize() << ")");
         ParGridFunction ug(const_cast<ParFiniteElementSpace *>(&fes_),
                            const_cast<real_t *>(u->GetData()));
         ParGridFunction vg(const_cast<ParFiniteElementSpace *>(&fes_),
                            const_cast<real_t *>(v->GetData()));
         Vector uv(3), vv(3);
         for (size_t i = 0; i < vol_files_.size(); i++)
         {
            const int s = vol_mine_[i];
            ug.GetVectorValue(vol_elem_[s], vol_ip_[s], uv);
            vg.GetVectorValue(vol_elem_[s], vol_ip_[s], vv);
            std::vector<real_t> row(kVolumeColumns);
            row[0] = t;
            for (int a = 0; a < 3; a++)
            {
               const std::array<real_t, 3> &e = cfg_.volume_axes[a];
               row[1 + a] = uv(0) * e[0] + uv(1) * e[1] + uv(2) * e[2];
               row[4 + a] = vv(0) * e[0] + vv(1) * e[1] + vv(2) * e[2];
            }
            vol_files_[i]->Row(row);
            vol_files_[i]->Flush();
         }
      }
   }

   /// The profile rows on their cadence (`force`: now).  Collective; `fill`
   /// runs on every rank when a row is due.
   bool Profiles(real_t t, real_t Vmax, const FieldsFn &fill, bool force = false)
   {
      if (cfg_.profiles.empty()) { return false; }
      const real_t dt = (Vmax > out_.fast_V) ? cfg_.profile_dt_fast : cfg_.profile_dt_slow;
      int due = (force || t - state_.profile_last >= dt * FdCycleOutput::kOutputTimeTolerance)
                ? 1 : 0;
      MPI_Allreduce(MPI_IN_PLACE, &due, 1, MPI_INT, MPI_MAX, comm_);
      if (!due || t == state_.profile_last) { return false; }
      StationFields f;
      fill(f);
      const int np = static_cast<int>(geom_.profile_locs.size());
      std::vector<real_t> mine(4 * static_cast<size_t>(np), 0.0), all(mine.size(), 0.0);
      for (int i = 0; i < np; i++)
      {
         if (prof_owner_[i] != rank_) { continue; }
         const StationLocation &L = geom_.profile_locs[i];
         real_t *r = &mine[4 * static_cast<size_t>(i)];
         r[0] = FdCycleOutput::Eval(f.slip_strike, L);
         r[1] = FdCycleOutput::Eval(f.slip_dip, L);
         r[2] = -(FdCycleOutput::Eval(f.tau_a_strike, L) + FdCycleOutput::Eval(f.tau_b_strike, L)
                  + f.eta * FdCycleOutput::Eval(f.V_strike, L)) / 1e6;
         r[3] = -(FdCycleOutput::Eval(f.tau_a_dip, L) + FdCycleOutput::Eval(f.tau_b_dip, L)
                  + f.eta * FdCycleOutput::Eval(f.V_dip, L)) / 1e6;
      }
      // Each point has one owner, so the sum is that owner's values.
      MPI_Reduce(mine.data(), all.data(), static_cast<int>(mine.size()),
                 MPITypeMap<real_t>::mpi_type, MPI_SUM, 0, comm_);
      if (rank_ == 0)
      {
         const real_t lv = Vmax > 0.0 ? std::log10(Vmax) : -300.0;
         int first = 0;
         for (size_t p = 0; p < cfg_.profiles.size(); p++)
         {
            const int n = cfg_.profiles[p].points;
            for (int q = 0; q < 4; q++)
            {
               std::vector<real_t> row(2 + n);
               row[0] = t;
               row[1] = lv;
               for (int i = 0; i < n; i++) { row[2 + i] = all[4 * static_cast<size_t>(first + i) + q]; }
               prof_files_[4 * p + q]->Row(row);
               prof_files_[4 * p + q]->Flush();
            }
            first += n;
         }
      }
      state_.profile_last = t;
      return true;
   }

   /// Rupture times: |V| of the owned DOFs at t, during an open event.
   void Track(real_t t, const Vector &speed)
   {
      if (!state_.open) { return; }
      for (int k = 0; k < speed.Size(); k++)
      {
         if (state_.trupt(k) < 0.0 && speed(k) >= rupture_V_) { state_.trupt(k) = t; }
      }
   }

   /// After FdCycleOutput::UpdateEvents at t, with what it returned.
   /// Collective when an event ends (the row and the rupture file).
   void OnEvents(const EventUpdate &ev, real_t t, const FieldsFn &fill, const SpeedFn &speed)
   {
      const int N = geom_.x2.Size();
      const EventChange change = ev.change;
      if (change == EventChange::Started)
      {
         StationFields f;
         fill(f);
         state_.open = true;
         state_.have_end = false;
         state_.trupt.SetSize(N);
         state_.trupt = -1.0;
         SlipTau_(f, state_.slip0, state_.tau0);
      }
      if (!state_.open) { return; }
      {
         Vector s;
         speed(s);
         Track(t, s);
      }
      if (change == EventChange::Below)
      {
         StationFields f;
         fill(f);
         SlipTau_(f, state_.slip_end, state_.tau_end);
         state_.have_end = true;
      }
      else if (change == EventChange::Above) { state_.have_end = false; }
      else if (change == EventChange::Ended)
      {
         if (!state_.have_end && ev.t_end == t)
         {
            // Below and Ended on one step (event_dwell = 0): the end is now.
            StationFields f;
            fill(f);
            SlipTau_(f, state_.slip_end, state_.tau_end);
            state_.have_end = true;
         }
         EndEvent_(ev);
      }
   }

   void Flush()
   {
      if (source_) { source_->Flush(); }
      if (events_) { events_->Flush(); }
      for (auto &f : vol_files_) { f->Flush(); }
      for (auto &f : prof_files_) { f->Flush(); }
   }

   State GetState()
   {
      State s = state_;
      s.files.clear();
      s.sizes.clear();
      auto add = [&s](RowFile &f)
      {
         s.files.push_back(std::filesystem::path(f.Path()).filename().string());
         s.sizes.push_back(f.Size());
      };
      if (source_) { add(*source_); }
      if (events_) { add(*events_); }
      for (auto &f : vol_files_) { add(*f); }
      for (auto &f : prof_files_) { add(*f); }
      if (!s.open)
      {
         for (Vector *v : {&s.trupt, &s.slip0, &s.tau0, &s.slip_end, &s.tau_end}) { v->SetSize(0); }
      }
      return s;
   }

private:
   bool InA_(real_t x2, real_t x3) const
   {
      // A boundary point counts as inside: tolerance relative to the region.
      const real_t tol = 1e-9 * (std::abs(cfg_.region_x2[0]) + std::abs(cfg_.region_x2[1]) +
                                 std::abs(cfg_.region_x3[0]) + std::abs(cfg_.region_x3[1]) + 1.0);
      return x2 >= cfg_.region_x2[0] - tol && x2 <= cfg_.region_x2[1] + tol &&
             x3 >= cfg_.region_x3[0] - tol && x3 <= cfg_.region_x3[1] + tol;
   }

   static real_t Speed_(const StationFields &f, int k)
   {
      return std::sqrt(f.V_strike(k) * f.V_strike(k) + f.V_dip(k) * f.V_dip(k));
   }

   /// Slip (strike, dip) [m] and shear stress (strike, dip) [MPa] per DOF.
   static void SlipTau_(const StationFields &f, Vector &slip, Vector &tau)
   {
      const int N = f.slip_dip.Size();
      slip.SetSize(2 * N);
      tau.SetSize(2 * N);
      for (int k = 0; k < N; k++)
      {
         slip(2 * k) = f.slip_strike(k);
         slip(2 * k + 1) = f.slip_dip(k);
         tau(2 * k) = -(f.tau_a_strike(k) + f.tau_b_strike(k) + f.eta * f.V_strike(k)) / 1e6;
         tau(2 * k + 1) = -(f.tau_a_dip(k) + f.tau_b_dip(k) + f.eta * f.V_dip(k)) / 1e6;
      }
   }

   void EndEvent_(const EventUpdate &ev)
   {
      const int event = ev.event;
      MFEM_VERIFY(state_.have_end, "BenchmarkOutput: event " << event << " ended (t_end = "
                  << ev.t_end << ") without a recorded first time below");
      const int N = geom_.x2.Size();
      // R, sum a tau(start), sum a tau(end), sum a ds: 7 sums.
      real_t loc[7] = {0, 0, 0, 0, 0, 0, 0}, glob[7];
      for (int k = 0; k < N; k++)
      {
         if (state_.trupt(k) < 0.0) { continue; }
         const real_t a = area_(k);
         loc[0] += a;
         for (int c = 0; c < 2; c++)
         {
            loc[1 + c] += a * state_.tau0(2 * k + c);
            loc[3 + c] += a * state_.tau_end(2 * k + c);
            loc[5 + c] += a * (state_.slip_end(2 * k + c) - state_.slip0(2 * k + c));
         }
      }
      MPI_Allreduce(loc, glob, 7, MPITypeMap<real_t>::mpi_type, MPI_SUM, comm_);
      const real_t R = glob[0];
      auto avg = [R](real_t a, real_t b) { return R > 0.0 ? std::hypot(a, b) / R : 0.0; };
      if (events_)
      {
         events_->Row({static_cast<real_t>(event), ev.t_start, ev.t_end, R,
                       avg(glob[1], glob[2]), avg(glob[3], glob[4]), avg(glob[5], glob[6])});
         events_->Flush();
      }
      // The rupture times of this event.
      Vector tab(6 * N);
      for (int k = 0; k < N; k++)
      {
         real_t *r = &tab(6 * k);
         r[0] = geom_.X(3 * k);
         r[1] = geom_.X(3 * k + 1);
         r[2] = geom_.X(3 * k + 2);
         r[3] = geom_.x2(k);
         r[4] = geom_.x3(k);
         r[5] = state_.trupt(k) >= 0.0 ? state_.trupt(k) : kNever;
      }
      std::ostringstream d;
      d << "fd_cycle rupture times of event " << event << " (first |V| >= " << rupture_V_
        << " m/s in the event; " << kNever << " = never), commit " << commit_;
      WriteFaultTable(Path_("_rupture_event" + std::to_string(event) + ".txt"), comm_, d.str(),
                      {"x(m)", "y(m)", "z(m)", "x2(m)", "x3(m)", "t_rupture(s)"}, tab);
      state_.open = false;
      state_.have_end = false;
   }

   std::string Path_(const std::string &suffix) const
   {
      return out_.dir + "/" + out_.prefix + suffix;
   }

   std::unique_ptr<RowFile> NewOrReopen_(const State *rs, const std::string &suffix,
                                         const std::string &desc,
                                         const std::vector<std::string> &cols)
   {
      const std::string path = Path_(suffix);
      if (!rs) { return std::make_unique<RowFile>(path, desc, cols); }
      const std::string name = std::filesystem::path(path).filename().string();
      for (size_t i = 0; i < rs->files.size(); i++)
      {
         if (rs->files[i] == name)
         {
            return std::make_unique<RowFile>(path, static_cast<int>(cols.size()), rs->sizes[i]);
         }
      }
      MFEM_ABORT("BenchmarkOutput restart: the checkpoint of rank " << rank_ << " does not "
                 "list " << name << "; the rank count or [cycle.benchmark] differs from the "
                 "run that wrote it");
      return nullptr;
   }

   void Open_(const State *rs, const std::string &commit)
   {
      commit_ = commit;
      std::ostringstream reg;
      reg << "x2 in [" << cfg_.region_x2[0] << ", " << cfg_.region_x2[1] << "] m, x3 in ["
          << cfg_.region_x3[0] << ", " << cfg_.region_x3[1] << "] m";
      if (rank_ == 0)
      {
         source_ = NewOrReopen_(rs, "_source.txt",
                                "fd_cycle source time series over the domain A (" + reg.str() +
                                "; area " + Num_(area_A_) + " m^2, mu " + Num_(mu_) +
                                " Pa), commit " + commit,
                                {"time(s)", "log10(Vmax_A)(m/s)", "moment_rate_A(N*m/s)"});
         events_ = NewOrReopen_(rs, "_events.txt",
                                "fd_cycle events (start: max V > " + Num_(out_.event_V) +
                                " m/s; end: " + Num_(out_.event_dwell) + " s below it; rupture "
                                "area: DOFs with |V| >= " + Num_(rupture_V_) + " m/s), commit " +
                                commit,
                                {"event", "t_start(s)", "t_end(s)", "rupture_area(m^2)",
                                 "tau_avg_start(MPa)", "tau_avg_end(MPa)", "slip_avg(m)"});
      }

      // Off-fault stations: the owner of each point.
      const int nv = static_cast<int>(cfg_.volume_stations.size());
      vol_owner_.assign(nv, -1);
      if (nv > 0)
      {
         ParMesh &pm = *fes_.GetParMesh();
         DenseMatrix pts(3, nv);
         for (int s = 0; s < nv; s++)
         {
            for (int d = 0; d < 3; d++) { pts(d, s) = cfg_.volume_stations[s].x[d]; }
         }
         Array<int> elem;
         Array<IntegrationPoint> ips;
         pm.FindPoints(pts, elem, ips, false);
         vol_elem_.assign(nv, -1);
         vol_ip_.resize(nv);
         std::vector<int> cand(nv);
         for (int s = 0; s < nv; s++)
         {
            cand[s] = (elem[s] >= 0) ? rank_ : std::numeric_limits<int>::max();
            if (elem[s] >= 0) { vol_elem_[s] = elem[s]; vol_ip_[s] = ips[s]; }
         }
         MPI_Allreduce(cand.data(), vol_owner_.data(), nv, MPI_INT, MPI_MIN, comm_);
         for (int s = 0; s < nv; s++)
         {
            const StationSpec &st = cfg_.volume_stations[s];
            MFEM_VERIFY(vol_owner_[s] != std::numeric_limits<int>::max(),
                        "BenchmarkOutput: volume station " << st.name << " at (" << st.x[0]
                        << ", " << st.x[1] << ", " << st.x[2] << ") m lies in no element");
            if (vol_owner_[s] != rank_) { continue; }
            vol_mine_.push_back(s);
            std::ostringstream d;
            d << "fd_cycle off-fault time series at " << st.name << " (x=" << st.x[0] << " m, y="
              << st.x[1] << " m, z=" << st.x[2] << " m); components along e1 = ("
              << Vec_(cfg_.volume_axes[0]) << "), e2 = (" << Vec_(cfg_.volume_axes[1])
              << "), e3 = (" << Vec_(cfg_.volume_axes[2]) << "); QD velocity: quasi-static "
              "du/dt, commit " << commit;
            vol_files_.push_back(NewOrReopen_(rs, "_" + st.name + ".txt", d.str(),
                                              {"time(s)", "disp_1(m)", "disp_2(m)", "disp_3(m)",
                                               "vel_1(m/s)", "vel_2(m/s)", "vel_3(m/s)"}));
         }
      }

      // Profile points: the stations' rule (StationOwners).
      const int np = static_cast<int>(geom_.profile_locs.size());
      int want = 0;
      for (const ProfileSpec &p : cfg_.profiles) { want += p.points; }
      MFEM_VERIFY(np == want, "BenchmarkOutput: " << np << " located profile points for "
                  << want << " in the deck");
      prof_owner_.assign(np, -1);
      if (np > 0)
      {
         std::vector<real_t> gmin;
         prof_owner_ = StationOwners(geom_.profile_locs, comm_, gmin);
         for (int i = 0; i < np; i++)
         {
            MFEM_VERIFY(gmin[i] < std::numeric_limits<real_t>::max(),
                        "BenchmarkOutput: profile point " << i << " found no fault DOF");
         }
         if (rank_ == 0)
         {
            real_t worst = 0.0;
            for (int i = 0; i < np; i++) { worst = std::max(worst, gmin[i]); }
            std::cout << "  profile points: " << np << ", largest match distance " << worst
                      << " m\n";
         }
      }
      if (rank_ == 0)
      {
         for (const ProfileSpec &p : cfg_.profiles)
         {
            std::ostringstream pos;
            for (int i = 0; i < p.points; i++)
            {
               const real_t s = static_cast<real_t>(i) / (p.points - 1);
               const real_t a = p.from[p.coordinate] + s * (p.to[p.coordinate] - p.from[p.coordinate]);
               pos << (i ? " " : "") << a;
            }
            for (int q = 0; q < 4; q++)
            {
               std::ostringstream d;
               d << "fd_cycle profile " << p.name << ", " << kProfileQuantities[q]
                 << (q < 2 ? " [m]" : " [MPa]") << " at " << p.points << " points from (x2, x3) = ("
                 << p.from[0] << ", " << p.from[1] << ") m to (" << p.to[0] << ", " << p.to[1]
                 << ") m; first row: 0 0 and each point's " << (p.coordinate == 0 ? "x2" : "x3")
                 << " [m], commit " << commit;
               std::vector<std::string> cols = {"time(s)", "log10(Vmax)(m/s)"};
               for (int i = 0; i < p.points; i++)
               {
                  cols.push_back(std::string(kProfileQuantities[q]) + "_" + std::to_string(i));
               }
               prof_files_.push_back(NewOrReopen_(
                  rs, "_profile_" + p.name + "_" + kProfileQuantities[q] + ".txt", d.str(), cols));
               if (!rs)
               {
                  std::vector<real_t> first = {0.0, 0.0};
                  for (int i = 0; i < p.points; i++)
                  {
                     const real_t s = static_cast<real_t>(i) / (p.points - 1);
                     first.push_back(p.from[p.coordinate] +
                                     s * (p.to[p.coordinate] - p.from[p.coordinate]));
                  }
                  prof_files_.back()->Row(first);
                  prof_files_.back()->Flush();
               }
            }
         }
      }
   }

   static std::string Num_(real_t v)
   {
      std::ostringstream o;
      o << v;
      return o.str();
   }
   static std::string Vec_(const std::array<real_t, 3> &e)
   {
      std::ostringstream o;
      o << e[0] << ", " << e[1] << ", " << e[2];
      return o.str();
   }

   BenchmarkConfig cfg_;
   OutputConfig out_;
   MPI_Comm comm_;
   const ParFiniteElementSpace &fes_;
   BenchFaultGeometry geom_;
   real_t mu_ = 0.0;
   real_t rupture_V_ = 0.0;
   int rank_ = 0;
   std::string commit_;
   Vector area_, wA_;
   std::vector<char> inA_;
   real_t fault_area_ = 0.0, area_A_ = 0.0;
   std::unique_ptr<RowFile> source_, events_;
   std::vector<int> vol_owner_, vol_elem_, vol_mine_;
   std::vector<IntegrationPoint> vol_ip_;
   std::vector<std::unique_ptr<RowFile>> vol_files_;
   std::vector<int> prof_owner_;
   std::vector<std::unique_ptr<RowFile>> prof_files_;
   State state_;
};

inline const char *const BenchmarkOutput::kProfileQuantities[4] = {"slip_2", "slip_3",
                                                                   "stress_2", "stress_3"};

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_BENCH_OUTPUT_HPP

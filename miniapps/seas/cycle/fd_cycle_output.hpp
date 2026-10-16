// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// fd_cycle_output.hpp — fully dynamic BP5 cycle (plan rev 3, Section 3.3 h):
// FdCycleOutput, one continuous set of output files across both regimes.
//
// Files, all under <dir>/<prefix>:
//   _<station>.txt   one per station, SCEC on-fault columns
//   _global.txt      time and log10 of the global max slip rate, every step
//   _catalog.txt     one row per event (root)
// The station and global files have the format of seas_bp5_full
// (io/bp5_parallel_output.hpp, io/probe_output.hpp): two '#' header lines,
// then rows of std::scientific, uppercase, 12 digits, space separated.  The
// station rows are computed term by term as ParallelBP5BenchmarkOutput::Write
// does, so the QD regime reproduces seas_bp5_full's rows:
//   slip_s = E(slip_strike)            slip_d = E(slip_dip)
//   V_s    = |E(V_strike)|             V_d    = |E(V_dip)|
//   tau_s  = -(E(tau_a_s) + E(tau_b_s) + eta E(V_strike)) / 1e6   [MPa]
//   tau_d  = -(E(tau_a_d) + E(tau_b_d) + eta E(V_dip)) / 1e6
//   row    = t, slip_s, slip_d, log10 V_s, log10 V_d, tau_s, tau_d, log10 theta
// E is the station evaluation on owned fault DOFs given by a regime adapter
// (StationLocation): a face-basis interpolation when the station lies on an
// owned fault face, else the nearest owned DOF.  Each station is written by
// the rank with the smallest match distance (ties: lowest rank).
//
// Cadence: a row when t - t_last >= 0.99 dt_out, dt_out = dt_fast while the
// global max V > fast_V, else dt_slow (deck values; the decision is reduced
// over ranks).  ForceWrite writes regardless.
//
// Restart: State() flushes and records each file's size; a writer built
// with a State truncates its files to those sizes and appends, so a run
// restarted from a checkpoint continues the files where the checkpoint
// left them.  The checkpoint must come from the same rank count, which
// keeps the station owners.
//
// No solver header (rule M4): fields and station locations arrive as plain
// vectors from the regime adapters.

#ifndef MFEM_SEAS_CYCLE_FD_CYCLE_OUTPUT_HPP
#define MFEM_SEAS_CYCLE_FD_CYCLE_OUTPUT_HPP

#include "mfem.hpp"
#include "cycle_config.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
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

/// Where one station sits among this rank's owned fault DOFs.
struct StationLocation
{
   real_t distance = std::numeric_limits<real_t>::max();  ///< max: no owned DOF
   bool exact = false;   ///< on an owned face: weights over that face's DOFs
   int start = -1;       ///< first DOF of that face
   Vector weights;       ///< face basis at the station
   int nearest = -1;     ///< nearest owned DOF otherwise
};

/// The station fields on owned fault DOFs, one value per DOF.
struct StationFields
{
   Vector slip_dip, slip_strike;
   Vector V_dip, V_strike;
   Vector tau_a_dip, tau_a_strike;   ///< e.g. the prestress
   Vector tau_b_dip, tau_b_strike;   ///< e.g. the elastic traction
   Vector theta;                     ///< state variable theta [s]
   real_t eta = 0.0;                 ///< multiplies E(V) in tau
};

/// One output file, appended row by row in the ProbeOutput format.
class RowFile
{
public:
   /// Fresh file: writes the two header lines.
   RowFile(const std::string &path, const std::string &description,
           const std::vector<std::string> &columns)
      : path_(path), ncol_(static_cast<int>(columns.size()))
   {
      file_.open(path, std::ios::out | std::ios::trunc);
      MFEM_VERIFY(file_.is_open(), "FdCycleOutput: cannot open " << path);
      file_ << "# " << description << "\n" << "# Columns:";
      for (size_t i = 0; i < columns.size(); i++)
      {
         file_ << " " << columns[i];
         if (i + 1 < columns.size()) { file_ << ","; }
      }
      file_ << "\n";
      file_ << std::scientific << std::uppercase << std::setprecision(12);
   }

   /// Existing file from a checkpoint: truncated to `size` bytes, appended.
   RowFile(const std::string &path, int ncol, long long size)
      : path_(path), ncol_(ncol)
   {
      MFEM_VERIFY(std::filesystem::exists(path),
                  "FdCycleOutput restart: " << path << " does not exist");
      const long long have = static_cast<long long>(std::filesystem::file_size(path));
      MFEM_VERIFY(have >= size, "FdCycleOutput restart: " << path << " has " << have
                  << " bytes, the checkpoint recorded " << size);
      std::filesystem::resize_file(path, static_cast<std::uintmax_t>(size));
      file_.open(path, std::ios::out | std::ios::app);
      MFEM_VERIFY(file_.is_open(), "FdCycleOutput: cannot reopen " << path);
      file_ << std::scientific << std::uppercase << std::setprecision(12);
   }

   void Row(const std::vector<real_t> &v)
   {
      MFEM_ASSERT(static_cast<int>(v.size()) == ncol_, "RowFile: column count");
      for (int i = 0; i < ncol_; i++)
      {
         if (i > 0) { file_ << " "; }
         file_ << v[i];
      }
      file_ << "\n";
   }

   void Flush() { file_.flush(); }
   long long Size()
   {
      file_.flush();
      return static_cast<long long>(std::filesystem::file_size(path_));
   }
   const std::string &Path() const { return path_; }

private:
   std::string path_;
   int ncol_;
   std::ofstream file_;
};

class FdCycleOutput
{
public:
   static constexpr real_t kOutputTimeTolerance = 0.99;

   /// What a checkpoint keeps to continue the files.
   struct State
   {
      real_t last_write_time = -1e30;
      int num_events = 0;
      bool in_event = false;
      real_t event_start = 0.0;
      real_t event_peak_V = 0.0;
      real_t below_since = -1.0;            ///< < 0: not below event_V in an event
      std::vector<std::string> files;       ///< this rank's files, by name
      std::vector<long long> sizes;         ///< their sizes in bytes
   };

   using FieldsFn = std::function<void(StationFields &)>;

   /// `locs[s]` is station s located on this rank's owned DOFs.  With
   /// `restart`, files continue from the checkpointed sizes.
   FdCycleOutput(const OutputConfig &cfg, MPI_Comm comm,
                 const std::vector<StationLocation> &locs,
                 const std::string &commit, const State *restart = nullptr)
      : cfg_(cfg), comm_(comm), locs_(locs)
   {
      MPI_Comm_rank(comm_, &rank_);
      const int ns = static_cast<int>(cfg_.stations.size());
      MFEM_VERIFY(static_cast<int>(locs_.size()) == ns,
                  "FdCycleOutput: " << locs_.size() << " station locations for "
                  << ns << " stations");

      // Ownership: smallest distance, then lowest rank (as seas_bp5_full).
      std::vector<real_t> dist(ns), gmin(ns);
      for (int s = 0; s < ns; s++) { dist[s] = locs_[s].distance; }
      MPI_Allreduce(dist.data(), gmin.data(), ns, MPITypeMap<real_t>::mpi_type,
                    MPI_MIN, comm_);
      std::vector<int> cand(ns), owner(ns);
      for (int s = 0; s < ns; s++)
      {
         cand[s] = (dist[s] == gmin[s]) ? rank_ : std::numeric_limits<int>::max();
      }
      MPI_Allreduce(cand.data(), owner.data(), ns, MPI_INT, MPI_MIN, comm_);

      if (restart) { state_ = *restart; }
      const std::vector<std::string> columns = {
         "time(s)", "slip_strike(m)", "slip_dip(m)",
         "log10(V_strike)(m/s)", "log10(V_dip)(m/s)",
         "tau_strike(MPa)", "tau_dip(MPa)", "log10(state)(s)"};
      for (int s = 0; s < ns; s++)
      {
         MFEM_VERIFY(gmin[s] < std::numeric_limits<real_t>::max(),
                     "FdCycleOutput: station " << cfg_.stations[s].name
                     << " found no fault DOF on any rank");
         if (owner[s] != rank_) { continue; }
         owned_.push_back(s);
         const std::string path = Path("_" + cfg_.stations[s].name + ".txt");
         if (restart)
         {
            stations_.push_back(Reopen(path, static_cast<int>(columns.size())));
         }
         else
         {
            const StationSpec &st = cfg_.stations[s];
            std::ostringstream d;
            d << "fd_cycle time series at " << st.name << " (x=" << st.x[0]
              << " m, y=" << st.x[1] << " m, z=" << st.x[2] << " m), commit " << commit;
            stations_.push_back(std::make_unique<RowFile>(path, d.str(), columns));
         }
      }
      if (rank_ == 0)
      {
         const std::string gpath = Path("_global.txt");
         const std::string cpath = Path("_catalog.txt");
         if (restart)
         {
            global_ = Reopen(gpath, 2);
            catalog_ = Reopen(cpath, 4);
         }
         else
         {
            global_ = std::make_unique<RowFile>(
               gpath, "fd_cycle global output, commit " + commit,
               std::vector<std::string>{"time(s)", "log10(Vmax)(m/s)"});
            catalog_ = std::make_unique<RowFile>(
               cpath, "fd_cycle event catalog (start: max V > " + Num(cfg_.event_V) +
               " m/s; end: " + Num(cfg_.event_dwell) + " s below it), commit " + commit,
               std::vector<std::string>{"event", "t_start(s)", "t_end(s)",
                                        "peak_Vmax(m/s)"});
         }
         std::cout << "  Station owners:\n";
         for (int s = 0; s < ns; s++)
         {
            std::cout << "    " << cfg_.stations[s].name << " -> rank " << owner[s]
                      << " (dist " << gmin[s] << " m)\n";
         }
      }
   }

   /// Station rows on the cadence; collective.  `fill` runs only on ranks
   /// that own a station, and only when a row is due.
   bool Write(real_t t, real_t Vmax, const FieldsFn &fill)
   {
      const real_t dt_out = (Vmax > cfg_.fast_V) ? cfg_.dt_fast : cfg_.dt_slow;
      int due = (t - state_.last_write_time >= dt_out * kOutputTimeTolerance) ? 1 : 0;
      MPI_Allreduce(MPI_IN_PLACE, &due, 1, MPI_INT, MPI_MAX, comm_);
      if (!due) { return false; }
      if (!owned_.empty())
      {
         StationFields f;
         fill(f);
         auto safe_log10 = [](real_t v) -> real_t
         {
            if (!std::isfinite(v)) { return -999.0; }
            return v > 0.0 ? std::log10(v) : -300.0;
         };
         for (size_t k = 0; k < owned_.size(); k++)
         {
            const StationLocation &L = locs_[owned_[k]];
            const real_t sd = Eval(f.slip_dip, L);
            const real_t ss = Eval(f.slip_strike, L);
            const real_t vd = std::abs(Eval(f.V_dip, L));
            const real_t vs = std::abs(Eval(f.V_strike, L));
            const real_t tau_d = -(Eval(f.tau_a_dip, L) + Eval(f.tau_b_dip, L)
                                   + f.eta * Eval(f.V_dip, L)) / 1e6;
            const real_t tau_s = -(Eval(f.tau_a_strike, L) + Eval(f.tau_b_strike, L)
                                   + f.eta * Eval(f.V_strike, L)) / 1e6;
            const real_t th = Eval(f.theta, L);
            stations_[k]->Row({t, ss, sd, safe_log10(vs), safe_log10(vd),
                               tau_s, tau_d, safe_log10(th)});
         }
      }
      state_.last_write_time = t;
      Flush();
      return true;
   }

   /// A station row now, whatever the cadence; collective.
   void ForceWrite(real_t t, real_t Vmax, const FieldsFn &fill)
   {
      state_.last_write_time = -1e30;
      Write(t, Vmax, fill);
   }

   /// One global row (root writes; others return).
   void WriteGlobal(real_t t, real_t Vmax)
   {
      if (!global_) { return; }
      global_->Row({t, Vmax > 0.0 ? std::log10(Vmax) : -300.0});
   }

   /// Event catalog on the global max V; the same on every rank.
   void UpdateEvents(real_t t, real_t Vmax)
   {
      State &s = state_;
      if (!s.in_event)
      {
         if (Vmax > cfg_.event_V)
         {
            s.in_event = true;
            s.num_events++;
            s.event_start = t;
            s.event_peak_V = Vmax;
            s.below_since = -1.0;
         }
         return;
      }
      s.event_peak_V = std::max(s.event_peak_V, Vmax);
      if (Vmax > cfg_.event_V) { s.below_since = -1.0; return; }
      if (s.below_since < 0.0) { s.below_since = t; }
      if (t - s.below_since >= cfg_.event_dwell)
      {
         if (catalog_)
         {
            catalog_->Row({static_cast<real_t>(s.num_events), s.event_start,
                           s.below_since, s.event_peak_V});
            catalog_->Flush();
         }
         s.in_event = false;
         s.below_since = -1.0;
      }
   }

   void Flush()
   {
      for (auto &f : stations_) { f->Flush(); }
      if (global_) { global_->Flush(); }
      if (catalog_) { catalog_->Flush(); }
   }

   /// The schedule, the event state and this rank's file sizes, flushed.
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
      for (auto &f : stations_) { add(*f); }
      if (global_) { add(*global_); }
      if (catalog_) { add(*catalog_); }
      return s;
   }

   int NumEvents() const { return state_.num_events; }
   bool InEvent() const { return state_.in_event; }

   /// Station value as Probe2DInterpolator::EvaluateScalar computes it.
   static real_t Eval(const Vector &field, const StationLocation &L)
   {
      if (L.exact)
      {
         real_t value = 0.0;
         for (int k = 0; k < L.weights.Size(); k++)
         {
            value += L.weights(k) * field(L.start + k);
         }
         return value;
      }
      return (L.nearest >= 0) ? field(L.nearest) : 0.0;
   }

private:
   std::string Path(const std::string &suffix) const
   {
      return cfg_.dir + "/" + cfg_.prefix + suffix;
   }

   static std::string Num(real_t v)
   {
      std::ostringstream o;
      o << v;
      return o.str();
   }

   std::unique_ptr<RowFile> Reopen(const std::string &path, int ncol)
   {
      const std::string name = std::filesystem::path(path).filename().string();
      for (size_t i = 0; i < state_.files.size(); i++)
      {
         if (state_.files[i] == name)
         {
            return std::make_unique<RowFile>(path, ncol, state_.sizes[i]);
         }
      }
      MFEM_ABORT("FdCycleOutput restart: the checkpoint of rank " << rank_
                 << " does not list " << name << "; the rank count or the "
                 "station list differs from the run that wrote it");
      return nullptr;
   }

   OutputConfig cfg_;
   MPI_Comm comm_;
   int rank_ = 0;
   std::vector<StationLocation> locs_;
   std::vector<int> owned_;
   std::vector<std::unique_ptr<RowFile>> stations_;
   std::unique_ptr<RowFile> global_;
   std::unique_ptr<RowFile> catalog_;
   State state_;
};

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_FD_CYCLE_OUTPUT_HPP

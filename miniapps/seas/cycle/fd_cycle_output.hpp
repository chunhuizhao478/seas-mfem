// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// fd_cycle_output.hpp — fully dynamic BP5 cycle (plan rev 3, Section 3.3 h):
// FdCycleOutput, one continuous set of output files across both regimes.
//
// Files, all under <dir>/<prefix>:
//   _<station>.txt   one per station, SCEC on-fault columns
//   _global.txt      time and log10 of the global max slip rate: every QD
//                    step; in the dynamic regime with each station row
//                    (the driver calls WriteGlobal when Write wrote)
//   _catalog.txt     one row per event (root)
//   _fault.txt       WriteFaultTable: one row per fault DOF (P2: rupture
//                    time, final slip, peak V and psi at the end of the
//                    first dynamic phase)
//   _switch.txt      OpenSwitchLog / SwitchLog (root, switching on): one
//                    line per regime change or stop
// Sample returns the station rows Write would write, without writing (the
// switch probe of plan P3).  The SCEC outputs besides the station files
// ([cycle.benchmark]) are BenchmarkOutput's (bench_output.hpp); it follows
// this writer's events through the EventUpdate UpdateEvents returns.
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
// (StationLocation, from LocateOnFaces): a face-basis interpolation when the
// station lies on an owned fault face, else the nearest owned DOF.  Each
// station is written by one rank (StationOwners).  A station where several
// fault faces meet (an edge or a vertex; BP5's depth-0 stations sit on
// vertices of 3 faces on the 1000 m mesh) takes the face LocateOnFaces
// picks from the faces' geometry, so its rows do not depend on the
// partition; there seas_bp5_full and Tandem take the first face found and
// the lowest rank, so their rows can differ by the jump of the traction
// between the faces (fully dynamic BP5 cycle: check 1092146, station
// strk-16dp+00, 576 against 144 ranks).
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

#include <algorithm>
#include <array>
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
   /// The face used (the exact match's, else the nearest DOF's): the
   /// distance from the station to the centre of the face's vertices, then
   /// that centre's x2 and x3.  It comes from the face's geometry alone, so
   /// it is the same on any partition; after the distance it orders the
   /// candidates (LocateOnFaces, StationOwners).
   std::array<real_t, 3> key = {{std::numeric_limits<real_t>::max(),
                                 std::numeric_limits<real_t>::max(),
                                 std::numeric_limits<real_t>::max()}};
};

/// Locate points (x2, x3) of the fault plane among this rank's owned fault
/// DOFs, given by their plane coordinates x2, x3 in blocks of nbf per face,
/// in the face basis's order (GaussLobatto triangle: DOFs 0, 1, 2 at the
/// vertices).  A point on an owned face (barycentric coordinates >= -1e-10,
/// as Probe2DInterpolator) gets that face's basis at the point and distance
/// 0; any other point gets the nearest owned DOF.  Where several faces
/// qualify (a point on an edge or a vertex the faces share, or DOFs of
/// several faces at the nearest position), the face whose vertex centre is
/// nearest the point wins, then the one with the smaller centre x2, then
/// x3 (StationLocation::key).  The centres are sums of the vertices in
/// sorted order, so the choice depends on neither the partition nor the
/// order of the faces or of their vertices.  Probe2DInterpolator takes the
/// first face it finds instead.  Without face blocks (nbf < 3, or N not a
/// multiple of nbf) every point gets the nearest DOF, keyed (0, x2, x3) by
/// that DOF's position.
inline std::vector<StationLocation> LocateOnFaces(
   const Vector &x2, const Vector &x3, int nbf,
   const std::vector<std::array<real_t, 2>> &pts)
{
   MFEM_VERIFY(x2.Size() == x3.Size(), "LocateOnFaces: " << x2.Size() << " x2 and "
               << x3.Size() << " x3 values");
   const int N = x2.Size();
   std::vector<StationLocation> out(pts.size());
   if (N == 0) { return out; }
   int order = -1;
   if (nbf >= 3 && N % nbf == 0)
   {
      const int p = static_cast<int>(std::llround((std::sqrt(1.0 + 8.0 * nbf) - 3.0) / 2.0));
      if ((p + 1) * (p + 2) / 2 == nbf) { order = p; }
   }
   const int nf = (order >= 1) ? N / nbf : 0;
   std::vector<std::array<real_t, 2>> centre(nf);
   for (int f = 0; f < nf; f++)
   {
      std::array<std::array<real_t, 2>, 3> v;
      for (int k = 0; k < 3; k++) { v[k] = {{x2(f * nbf + k), x3(f * nbf + k)}}; }
      std::sort(v.begin(), v.end());
      centre[f] = {{(v[0][0] + v[1][0] + v[2][0]) / 3.0, (v[0][1] + v[1][1] + v[2][1]) / 3.0}};
   }
   auto key_of = [&](int f, const std::array<real_t, 2> &p)
   {
      const real_t d2 = p[0] - centre[f][0], d3 = p[1] - centre[f][1];
      return std::array<real_t, 3> {{std::sqrt(d2 * d2 + d3 * d3), centre[f][0],
                                     centre[f][1]}};
   };
   std::unique_ptr<H1_TriangleElement> fe;
   Vector shape;
   if (nf > 0)
   {
      fe.reset(new H1_TriangleElement(order, BasisType::GaussLobatto));
      shape.SetSize(nbf);
   }
   const real_t inf = std::numeric_limits<real_t>::max();
   for (size_t s = 0; s < pts.size(); s++)
   {
      const std::array<real_t, 2> &p = pts[s];
      StationLocation &L = out[s];
      // The nearest DOF: the smallest distance, then its face's key.
      std::array<real_t, 4> best = {{inf, inf, inf, inf}};
      for (int i = 0; i < N; i++)
      {
         const real_t d2 = x2(i) - p[0], d3 = x3(i) - p[1];
         const std::array<real_t, 3> k = (nf > 0) ? key_of(i / nbf, p)
                                         : std::array<real_t, 3> {{0.0, x2(i), x3(i)}};
         const std::array<real_t, 4> c = {{std::sqrt(d2 * d2 + d3 * d3), k[0], k[1], k[2]}};
         if (c < best) { best = c; L.nearest = i; }
      }
      L.distance = best[0];
      L.key = {{best[1], best[2], best[3]}};
      // The faces that contain the point: the smallest key.
      for (int f = 0; f < nf; f++)
      {
         const int b = f * nbf;
         const real_t ax = x2(b + 1) - x2(b), az = x3(b + 1) - x3(b);
         const real_t bx = x2(b + 2) - x2(b), bz = x3(b + 2) - x3(b);
         const real_t px = p[0] - x2(b), pz = p[1] - x3(b);
         const real_t det = ax * bz - az * bx;
         if (std::abs(det) < 1e-14) { continue; }
         const real_t r = (px * bz - pz * bx) / det, q = (ax * pz - az * px) / det;
         const real_t tol = 1e-10;
         if (1.0 - r - q < -tol || r < -tol || q < -tol) { continue; }
         const std::array<real_t, 3> k = key_of(f, p);
         if (L.exact && !(k < L.key)) { continue; }
         IntegrationPoint ip;
         ip.x = r;
         ip.y = q;
         ip.weight = 0.0;
         fe->CalcShape(ip, shape);
         L.exact = true;
         L.start = b;
         L.weights = shape;
         L.distance = 0.0;
         L.key = k;
      }
   }
   return out;
}

/// The rank that writes each located point: the smallest distance, then the
/// smallest key (LocateOnFaces), then the lowest rank.  A face belongs to
/// one rank, so on any partition the owner holds the same face.  gmin: the
/// smallest distance per point (max: no rank found the point); gkey, if
/// given: the owner's key.  Collective.
inline std::vector<int> StationOwners(const std::vector<StationLocation> &locs, MPI_Comm comm,
                                      std::vector<real_t> &gmin,
                                      std::vector<std::array<real_t, 3>> *gkey = nullptr)
{
   int rank = 0;
   MPI_Comm_rank(comm, &rank);
   const int n = static_cast<int>(locs.size());
   const real_t inf = std::numeric_limits<real_t>::max();
   std::vector<char> cand(n, 1);
   std::vector<real_t> mine(n), all(n);
   for (int c = 0; c < 4; c++)
   {
      for (int i = 0; i < n; i++)
      {
         mine[i] = cand[i] ? (c == 0 ? locs[i].distance : locs[i].key[c - 1]) : inf;
      }
      MPI_Allreduce(mine.data(), all.data(), n, MPITypeMap<real_t>::mpi_type, MPI_MIN, comm);
      if (c == 0) { gmin = all; }
      if (c > 0 && gkey)
      {
         gkey->resize(n);
         for (int i = 0; i < n; i++) { (*gkey)[i][c - 1] = all[i]; }
      }
      for (int i = 0; i < n; i++) { cand[i] = cand[i] && mine[i] == all[i]; }
   }
   std::vector<int> r(n), owner(n);
   for (int i = 0; i < n; i++) { r[i] = cand[i] ? rank : std::numeric_limits<int>::max(); }
   MPI_Allreduce(r.data(), owner.data(), n, MPI_INT, MPI_MIN, comm);
   return owner;
}

/// The station fields on owned fault DOFs, one value per DOF.
struct StationFields
{
   Vector slip_dip, slip_strike;
   Vector V_dip, V_strike;
   Vector tau_a_dip, tau_a_strike;   ///< e.g. the prestress
   Vector tau_b_dip, tau_b_strike;   ///< e.g. the elastic traction
   Vector theta;                     ///< state variable theta [s]
   Vector sigma_n;                   ///< sigma_n_eff, compression positive [Pa]
   real_t eta = 0.0;                 ///< multiplies E(V) in tau
};

/// The regime's displacement u and velocity v, vectors of the QD regime's
/// displacement space (both regimes use it).  Collective: the QD velocity
/// is a solve (QDPhase::Velocity).
using VolumeFieldsFn = std::function<void(const Vector *&u, const Vector *&v)>;

/// What one UpdateEvents call changed.
enum class EventChange
{
   None,      ///< nothing (inside or outside an event)
   Started,   ///< max V passed event_V: a new event
   Below,     ///< in an event, the first time below event_V since the last above
   Above,     ///< in an event, above event_V again after Below
   Ended      ///< event_dwell s below: the event ended (its catalog row is written)
};

/// UpdateEvents' result: the change, the event number (1, 2, ...; the open
/// or just ended event) and its start and end times (t_end: the first time
/// below, valid with Ended).
struct EventUpdate
{
   EventChange change = EventChange::None;
   int event = 0;
   real_t t_start = 0.0;
   real_t t_end = 0.0;
};

/// What a checkpoint keeps of one snapshot file (cycle/snapshots.hpp).
struct SnapshotFileState
{
   std::string kind;        ///< "fault" or "surface"
   long long steps = 0;     ///< snapshots in the file
   long long qd_mark = -1;  ///< last QD mark written
   real_t last_t = -1e30;   ///< time of the last snapshot
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

   /// A line of text (the switch log's rows).
   void Text(const std::string &line) { file_ << line << "\n"; }

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
   /// Columns of a station row: t, slip_s, slip_d, log10 V_s, log10 V_d,
   /// tau_s, tau_d (MPa), log10 theta.
   static constexpr int kStationColumns = 8;

   /// What a checkpoint keeps to continue the files.
   struct State
   {
      real_t last_write_time = -1e30;
      int num_events = 0;
      bool in_event = false;
      real_t event_start = 0.0;
      real_t event_peak_V = 0.0;
      real_t below_since = -1.0;            ///< < 0: not below event_V in an event
      real_t last_above = -1.0;             ///< in an event: last max V >= event_V
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

      // Ownership: smallest distance, then the face key, then lowest rank.
      std::vector<real_t> gmin;
      std::vector<std::array<real_t, 3>> gkey;
      const std::vector<int> owner = StationOwners(locs_, comm_, gmin, &gkey);

      if (restart) { state_ = *restart; }
      restarted_ = (restart != nullptr);
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
                      << " (dist " << gmin[s] << " m, face centre (" << gkey[s][1] << ", "
                      << gkey[s][2] << ") m)\n";
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
         for (size_t k = 0; k < owned_.size(); k++)
         {
            stations_[k]->Row(StationRow_(t, f, locs_[owned_[k]]));
         }
      }
      state_.last_write_time = t;
      Flush();
      return true;
   }

   /// The rows Write would write now, for every station in deck order
   /// (kStationColumns values each), on the root; zeros on the other ranks.
   /// Writes nothing and leaves the cadence alone.  Collective.
   std::vector<real_t> Sample(real_t t, const FieldsFn &fill)
   {
      const int ns = static_cast<int>(cfg_.stations.size());
      std::vector<real_t> mine(static_cast<size_t>(kStationColumns * ns), 0.0);
      if (!owned_.empty())
      {
         StationFields f;
         fill(f);
         for (int s : owned_)
         {
            const std::vector<real_t> r = StationRow_(t, f, locs_[s]);
            std::copy(r.begin(), r.end(), mine.begin() + kStationColumns * s);
         }
      }
      std::vector<real_t> all(mine.size(), 0.0);
      // Each station has one owner, so the sum is that owner's row.
      MPI_Reduce(mine.data(), all.data(), static_cast<int>(mine.size()),
                 MPITypeMap<real_t>::mpi_type, MPI_SUM, 0, comm_);
      return all;
   }

   /// The switch log _switch.txt on the root: a fresh file, or with a
   /// restart State the checkpointed one, continued.  Collective (only the
   /// root opens a file).
   void OpenSwitchLog(const std::string &commit)
   {
      if (rank_ != 0) { return; }
      const std::string path = Path("_switch.txt");
      const std::string name = std::filesystem::path(path).filename().string();
      if (restarted_)
      {
         bool listed = false;
         for (const std::string &f : state_.files) { listed = listed || (f == name); }
         MFEM_VERIFY(listed, "FdCycleOutput restart: the checkpoint does not list " << name
                     << "; it was written with [cycle].switching = false");
         switch_ = Reopen(path, 5);
      }
      else
      {
         switch_ = std::make_unique<RowFile>(
            path, "fd_cycle switch log, commit " + commit,
            std::vector<std::string>{"time(s)", "Vmax(m/s)", "from", "to", "reason"});
      }
   }

   /// One switch-log line (root; no-op elsewhere or without OpenSwitchLog).
   void SwitchLog(real_t t, real_t Vmax, const std::string &from, const std::string &to,
                  const std::string &reason)
   {
      if (!switch_) { return; }
      std::ostringstream o;
      o << std::scientific << std::uppercase << std::setprecision(12) << t << " " << Vmax
        << " " << from << " " << to << " " << reason;
      switch_->Text(o.str());
      switch_->Flush();
   }

   /// A station row now, whatever the cadence; collective.  True (a row).
   bool ForceWrite(real_t t, real_t Vmax, const FieldsFn &fill)
   {
      state_.last_write_time = -1e30;
      return Write(t, Vmax, fill);
   }

   /// One global row (root writes; others return).
   void WriteGlobal(real_t t, real_t Vmax)
   {
      if (!global_) { return; }
      global_->Row({t, Vmax > 0.0 ? std::log10(Vmax) : -300.0});
   }

   /// Event catalog on the global max V; the same on every rank.  The row's
   /// t_end is the first time below event_V.  The dwell counts from the last
   /// time at or above it, as RegimeController's does from the same t, so a
   /// dynamic phase that ends on the same threshold and dwell closes the
   /// event on the same step.  Returns what changed (EventUpdate).
   EventUpdate UpdateEvents(real_t t, real_t Vmax)
   {
      State &s = state_;
      EventUpdate u;
      if (!s.in_event)
      {
         if (Vmax > cfg_.event_V)
         {
            s.in_event = true;
            s.num_events++;
            s.event_start = t;
            s.event_peak_V = Vmax;
            s.below_since = -1.0;
            s.last_above = t;
            u.change = EventChange::Started;
         }
         u.event = s.num_events;
         u.t_start = s.event_start;
         return u;
      }
      u.event = s.num_events;
      u.t_start = s.event_start;
      s.event_peak_V = std::max(s.event_peak_V, Vmax);
      if (Vmax >= cfg_.event_V)
      {
         if (s.below_since >= 0.0) { u.change = EventChange::Above; }
         s.below_since = -1.0;
         s.last_above = t;
         return u;
      }
      if (s.below_since < 0.0)
      {
         s.below_since = t;
         u.change = EventChange::Below;
      }
      if (t - s.last_above >= cfg_.event_dwell)
      {
         if (catalog_)
         {
            catalog_->Row({static_cast<real_t>(s.num_events), s.event_start,
                           s.below_since, s.event_peak_V});
            catalog_->Flush();
         }
         // Below and Ended on one step (event_dwell 0): the end fields are
         // this step's, which BenchmarkOutput takes with Ended too.
         if (u.change == EventChange::Below) { u.change = EventChange::Ended; u.t_end = t; }
         else { u.change = EventChange::Ended; u.t_end = s.below_since; }
         s.in_event = false;
         s.below_since = -1.0;
      }
      return u;
   }

   void Flush()
   {
      for (auto &f : stations_) { f->Flush(); }
      if (global_) { global_->Flush(); }
      if (catalog_) { catalog_->Flush(); }
      if (switch_) { switch_->Flush(); }
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
      if (switch_) { add(*switch_); }
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
   /// One station row, term by term as ParallelBP5BenchmarkOutput::Write.
   static std::vector<real_t> StationRow_(real_t t, const StationFields &f,
                                          const StationLocation &L)
   {
      auto safe_log10 = [](real_t v) -> real_t
      {
         if (!std::isfinite(v)) { return -999.0; }
         return v > 0.0 ? std::log10(v) : -300.0;
      };
      const real_t sd = Eval(f.slip_dip, L);
      const real_t ss = Eval(f.slip_strike, L);
      const real_t vd = std::abs(Eval(f.V_dip, L));
      const real_t vs = std::abs(Eval(f.V_strike, L));
      const real_t tau_d = -(Eval(f.tau_a_dip, L) + Eval(f.tau_b_dip, L)
                             + f.eta * Eval(f.V_dip, L)) / 1e6;
      const real_t tau_s = -(Eval(f.tau_a_strike, L) + Eval(f.tau_b_strike, L)
                             + f.eta * Eval(f.V_strike, L)) / 1e6;
      const real_t th = Eval(f.theta, L);
      return {t, ss, sd, safe_log10(vs), safe_log10(vd), tau_s, tau_d, safe_log10(th)};
   }

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
   std::unique_ptr<RowFile> switch_;
   bool restarted_ = false;
   State state_;
};

/// One file for the whole fault: `ncol` values per owned fault DOF on each
/// rank (`data`, DOF-major), gathered to the root and written in rank order
/// in the RowFile format.  Collective.
inline void WriteFaultTable(const std::string &path, MPI_Comm comm,
                            const std::string &description,
                            const std::vector<std::string> &columns, const Vector &data)
{
   const int ncol = static_cast<int>(columns.size());
   MFEM_VERIFY(ncol > 0 && data.Size() % ncol == 0, "WriteFaultTable: " << data.Size()
               << " values is not a multiple of " << ncol << " columns");
   int rank = 0, nr = 1;
   MPI_Comm_rank(comm, &rank);
   MPI_Comm_size(comm, &nr);
   int nloc = data.Size();
   std::vector<int> cnt(nr), dsp(nr, 0);
   MPI_Gather(&nloc, 1, MPI_INT, cnt.data(), 1, MPI_INT, 0, comm);
   std::vector<real_t> all;
   if (rank == 0)
   {
      for (int r = 1; r < nr; r++) { dsp[r] = dsp[r - 1] + cnt[r - 1]; }
      all.resize(static_cast<size_t>(dsp[nr - 1] + cnt[nr - 1]));
   }
   MPI_Gatherv(data.GetData(), nloc, MPITypeMap<real_t>::mpi_type, all.data(), cnt.data(),
               dsp.data(), MPITypeMap<real_t>::mpi_type, 0, comm);
   if (rank != 0) { return; }
   RowFile f(path, description, columns);
   std::vector<real_t> row(static_cast<size_t>(ncol));
   for (size_t i = 0; i + ncol <= all.size(); i += ncol)
   {
      for (int c = 0; c < ncol; c++) { row[c] = all[i + c]; }
      f.Row(row);
   }
   f.Flush();
}

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_FD_CYCLE_OUTPUT_HPP

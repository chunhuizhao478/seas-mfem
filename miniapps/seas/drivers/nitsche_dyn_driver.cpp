// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// nitsche_dyn_driver.cpp — Nitsche dynamic rupture, plan item N-07 (build
// sequence rev 5.2, phase P7): the driver.
//
// Reads the SAFS spatial deck through spatial::LoadSpatialFrictionConfig (the
// SAFS driver's own parser) plus the Nitsche-only keys the SAFS parser does
// not know ([numerics] penalty_factor, the fault-coordinate specification's
// [fault_geometry] keys, [diagnostics] energy_balance), builds the mesh, the
// material, the boundary configuration and the fault setup, constructs the
// NitscheWaveOperator (N-01) and advances it with MFEM's
// CentralDifferenceSolver (N-02).  Nothing benchmark-specific lives here: the
// law comes from a factory keyed on [meta] law and [friction.rate_state]
// state_evolution, the friction parameters from the SAFS resolver, the
// prestress and the nucleation from the deck as global tensors (spec I1),
// the stations from a station file the deck names.
//
// Outputs (in --output-dir):
//   <prefix>_station_<name>.dat   FaultStationWriter (the column layout of the
//                                 existing runs) through a DOFData adapter;
//   <prefix>_station_snap.txt     per station: requested point, snapped fault
//                                 node, distance;
//   <prefix>_station_conventions.txt  frame, pin and sign map (I3);
//   <prefix>_station_<name>_copies.dat  per station, the min/max over the DG
//                                 copies at the station point (--station-copies
//                                 average, the default);
//   <prefix>_fault_r<rank>.csv    fault-field dump every --dump-every steps
//                                 (with Prakash-Clifton one more column, sigma*);
//   <prefix>_rupture_time_r<rank>.csv  first |V| > 1 mm/s per fault DOF;
//   stdout log every --log-every steps.
//
// Station copies (Chunhui Zhao, 2026-10-14): a station sits on a fault point
// that carries one DOF per adjacent fault face (DG copies).  With
// --station-copies average (default) the station files hold the face-area-
// weighted average over those copies (summed on rank 0 in face-centroid
// order); --station-copies first keeps the writer's single nearest copy.
//
// Station output cost (speed-next O2, Chunhui Zhao 2026-10-20; the files are
// the same bit for bit): --station-comm sub (default) averages the copies
// over the ranks that hold them, in their own communicator, with the order
// of the sum computed once (dynamic/station_copy_average.hpp); world is the
// code before O2, a collective of every rank each station step.  The
// station and envelope rows reach the files every --station-flush-every n
// rows (default 100; 1 = every row, as before) and at the end; SIGTERM (the
// wall-time limit), SIGINT or SIGUSR1 makes every rank flush its files at
// the end of its step and write every row from there on, until the job is
// killed (a second SIGTERM or SIGINT ends it at once).  An abort inside the
// run (MFEM_VERIFY) loses up to n - 1 rows per file; --station-flush-every
// 1 keeps them.  Only the DOFs the stations read get their station data
// each station step.
//
// Fault frame filter (speed-next O1, Chunhui Zhao 2026-10-20):
// --paraview-fault-deflate-level n (0..9) sets the deflate level of
// fault.vtkhdf in place of the deck's [output] paraview_fault_deflate_level;
// none writes its datasets with no filter at all (no compression work on
// rank 0, a larger file).  Deflate is lossless, so the values are the same
// either way.  Not with [output] paraview_fault_zfp_tol.  Absent = the deck.
//
// Fault frame compression (fault-writer note rev1, W1, Chunhui Zhao
// 2026-10-24): with --fault-frame-compress spread (default) and deflate on,
// rank 0 cuts each changing field of a fault frame into the dataset's chunks
// and sends them to up to nranks - 1 other ranks, which shuffle and deflate
// them as HDF5 would; rank 0 writes the returned bytes as they are
// (vtkhdf::FaultHDFStream::SetSpread).  The file is the one rank 0 alone
// writes (rank0), byte for byte in its chunks; only the compression work
// moves off rank 0.  The fast fault writer only (--fast-fault-output 1);
// with no filter, ZFP or one rank, rank 0 writes as before.
//
// Rank balance (fault-writer note rev1, L1, Chunhui Zhao 2026-10-24): with
// --lts rate2 the run ends with one row per rank in
// <output-dir>/<prefix>_rank_balance.csv: the time in ode.Step after step 1,
// the halo wait, the non-wait time (their difference), the step loop's
// time outside ode.Step, the fault time, and the work counts (patch element
// kicks, tets, fault-face sides, halo values; per cluster tets, patch
// elements, interior and shared face blocks) with the node name.  The log
// gives max/mean of each, their correlation with the non-wait time, and a
// least-squares fit non-wait = a kicks + b sides + c (8 ranks or more).
// Slow-parts study (Chunhui Zhao 2026-10-25): the rows also hold the
// seconds of each piece of the ticks after step 1 (NitscheLtsRun phase
// timers, ph_* columns; ph_untimed is the step time they leave out) and the
// absorbing elements per cluster; the log gives each piece's mean, max and
// the slowest rank's value.
//
// Partition weights (L2 note rev1, Chunhui Zhao 2026-10-25): with --lts
// rate2 on more than one part, --lts-partition-weight tets (default) gives
// each tet weight 1 in its own cluster's METIS constraint, as before; blocks
// gives constraint c of a tet the dense blocks of cluster c's sweep it
// carries (its diagonal block when it is in the patch of c, one unit per
// face with a side in c; BuildLtsBlockWeights), so each cluster's block work
// is balanced per rank instead of its tets.  --fault-weight w adds 100 w per
// fault-face side on the tet's own cluster with either weight.  The
// partition line prints the block work per cluster (max/mean) and the tick
// model of the cut made, whichever weights made it.  blocks changes only the
// cut, so the stations move at round-off; the mesh part files record the
// weight.
//
// Setup phases (speed-next S1, Chunhui Zhao 2026-10-20): before step 1 the
// log prints the wall time of each setup phase, max and mean over the ranks
// (SetupTimer; one reduction, no barrier), and (2026-10-22) each rank's
// VmData at the end of each phase, max and mean over the ranks: the private
// writable memory that a node with strict overcommit (Vista) charges against
// its commit limit, and does not give back when the heap shrinks.  Where
// /proc/self/status is absent (macOS) that line is left out.
//
// Mesh setup (Chunhui Zhao 2026-10-22: partition-memory report rev1, option
// D; speed-next S1): with --mesh-setup rank0 (default) only rank 0 reads the
// serial mesh; it computes the partition (the same calls as before: the LTS
// multi-constraint cut, the --fault-weight cut, or MFEM's default METIS cut
// that ParMesh made on every rank), cuts the mesh into one part per rank and
// sends each rank its part with the LTS clusters of its elements
// (dynamic/mesh_parts.hpp).  The other ranks never hold the serial mesh, so
// a rank needs memory for its own part only.  The parallel mesh has the same
// elements, vertices, faces and DOF numbers as before; its shared faces are
// listed in another order, which can change sums over them in the last bit.
// --mesh-setup every-rank is the code before: every rank reads the serial
// mesh and builds its ParMesh from it.  On one rank both are the same.
//
// Mesh part files (speed-next S1 step 2, Chunhui Zhao 2026-10-22): run on
// one rank, --mesh-parts-out d --mesh-parts-n N reads the deck and the
// mesh, computes the partition into N parts as a run on N ranks would (LTS
// clusters included), writes the N parts, their LTS clusters and
// d/mesh_parts.txt (the mesh's size and FNV-1a hash, the options the
// partition depends on, the LTS header) and stops.  A run on N ranks with
// --mesh-parts d then reads its own part; no rank reads the serial mesh.
// The run refuses parts cut for another rank count, from another mesh, or
// with other values of the options the partition depends on (the fault
// weight and attribute; with LTS also the order, cfl, the LTS cap and fault
// buffer, dt_max and the material).  The parts are the ones --mesh-setup
// rank0 sends, so the results are the same bit for bit.
//
// Station columns (SCEC map, plan "The SCEC conventions, mapped"):
//   slip1 = delta.dhat (v-slip), slip2 = delta.shat (h-slip),
//   V1 = V.dhat, V2 = V.shat,
//   tau1 = -T.dhat, tau2 = -T.shat (SCEC shear tractions, Pa),
//   sigma_n = -T.nhat (compression positive, Pa),
// with T the total traction (Nitsche traction at the end of the step plus
// prestress plus nucleation) on the face with normal nhat.
//
// Optional [friction.prakash_clifton] table (Chunhui Zhao, 2026-10-14): v_star
// (m/s) and length (m), both required and > 0; its presence switches on the
// operator's Prakash-Clifton regularization of the normal stress the law
// sees, seeded with the initial sigma_n per DOF.  Absent = every bit as
// before.  The station files are unchanged (sigma_n there stays sigma_eff);
// sigma* goes to the fault CSV dump and to fault.vtkhdf (normal_stress_pc).
//
// Scope (Chunhui Zhao, 2026-10-14): tetrahedral meshes; constant material
// (the sidecar / depth-profile / halfspace selections are P8).
//
// Laws (benchmark plan rev 1.3, Chunhui Zhao, 2026-10-16): rate and state
// with the aging law (RateStateClosure); rate and state with the slip law and
// strong rate weakening (RateStateSRWClosure; the deck must type
// [friction.rate_state] f_w_default, because the parser fills 0.2 when the
// key is absent); linear slip weakening (SlipWeakeningClosure,
// [friction.slip_weakening], starting at rest: V = 0, slip path 0).
//
// Prestress regions (plan 4.5): optional top-level [[stress_region]] tables,
// read by this driver.  Each is a box (x_min_m, x_max_m, y_min_m, y_max_m,
// z_min_m, z_max_m; an absent bound is open, at least one is required) and a
// full global tensor (sigma_xx_pa, sigma_yy_pa, sigma_zz_pa, sigma_xy_pa,
// sigma_yz_pa, sigma_xz_pa, all six, tension positive as in [stress]) that
// replaces the background inside the box; the last match wins.  The
// operator projects the prestress face by face at quadrature points inside
// each face, so a node on a box edge carries each adjacent face's value in
// that face's DG copy.
//
// Stations (Chunhui Zhao, 2026-10-16): [output.stations] file = "..." (CSV
// name,x,y,z in the global frame, written by scripts/make_fault_stations.py)
// and snap_tolerance_m.  Each station snaps to the nearest fault node in
// 3-D; the run aborts, naming the station, when that node is farther than
// the tolerance.  No table, no station files.  The ninth column is
// log10 theta (aging law, with V_0, f_0, b from [friction.rate_state]), psi
// (strong rate weakening) or the slip path (slip weakening).
//
// --lts-report 1 (plan item 7, local time stepping, phase P0): rank 0 prints
// the rate-2 cluster report of solver/nitsche_lts_setup.hpp for the serial
// mesh and the run stops before the partition.  --lts-report-cap n picks the
// clustering of its partition section (hard cap n; 0 = no cap; default 4)
// and --lts-report-ranks a,b,... the rank counts there (default: this run's
// rank count when it is above 1).
//
// Local time stepping (plan item 7, P3 and P4; solver/nitsche_lts_run.hpp):
// [numerics] lts = "rate2" (or --lts rate2) steps cluster c at 2^c dt0 with
// the multirate leapfrog, dt0 = the global step; "off" (default) is
// CentralDifferenceSolver, every bit as before.  The clusters come from the
// serial mesh (rate-2 rules, no wiggle, hard cap [numerics] lts_nc_cap or
// --lts-cap, default 4 here), with the whole fault and --lts-fault-buffer
// layers of face neighbours around it (default 3) in one cluster c_F; on
// more than one rank the partition is BuildLtsAwarePartition (every cluster
// balanced per rank, fault faces on one rank, --fault-weight as cost).  The
// fault closes every 2^c_F ticks, and the stations, dumps and fault frames
// follow those fault steps.  The step record (max|V|, sigma ranges) is
// reduced over the ranks once per sync interval, as the extremes over its
// fault steps, so the log lines come at sync points.  With energy_balance
// the energy is evaluated at sync points and without the modified-energy
// term (the multirate scheme has no known discrete energy), so rho measures
// the time discretisation error.  Each rank posts one halo exchange per tick
// (checked at the end).  LTS needs cfl <= 0.6 (the accepted multirate
// limit, Chunhui Zhao 2026-10-19) and kv_scale = "node".  --lts-same-step 1
// (test only): the same clusters, every one at dt0.

#include "mfem.hpp"
#include "../domain/nitsche_wave_operator.hpp"
#include "../domain/nitsche_energy.hpp"
#include "../solver/nitsche_cfl.hpp"
#include "../friction/rate_state_closure.hpp"
#include "../friction/rate_state_srw_closure.hpp"
#include "../friction/slip_weakening_closure.hpp"
#include "../fault/fault_frame.hpp"
#include "../spatial/code/spatial_friction.hpp"
#include "../dynamic/spatial_nucleation.hpp"
#include "../dynamic/fault_station_file.hpp"
#include "../dynamic/station_copy_average.hpp"
#include "../dynamic/lts_partition.hpp"   // METIS with element weights (plan item 6)
#include "../dynamic/mesh_parts.hpp"      // --mesh-setup rank0 (header, "Mesh setup")
#include "../solver/nitsche_lts_setup.hpp"   // --lts-report (plan item 7, P0)
#include "../solver/nitsche_lts_run.hpp"     // --lts rate2 (plan item 7, P3/P4)
#ifdef MFEM_USE_HDF5
#include "../io/fault_vtkhdf_writer.hpp"
#include "../io/paraview_output.hpp"   // ParaViewOutput::ProbeH5ZZfpPluginOrAbort
#endif

#include <toml.hpp>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <utility>
#include <vector>

using namespace mfem;
using namespace mfem::seas;
namespace ff = mfem::seas::fault_frame;

namespace
{

/// P3's measured stable constants C(p), p = 1..4 (dt = cfl C(p) h/(p^2 c_p)).
const real_t kCeig[5] = {0.0, 0.113, 0.212, 0.282, 0.332};

int g_rank = 0;

std::string ArgStr(int argc, char *argv[], const char *flag, const std::string &def)
{
   for (int i = 1; i + 1 < argc; i++) { if (std::string(argv[i]) == flag) { return argv[i + 1]; } }
   return def;
}
/// --fault-weight default (plan item 6): extra METIS weight per fault-face
/// side, in units of one tet.
constexpr real_t kDefaultFaultWeight = 0.0;

real_t ArgReal(int argc, char *argv[], const char *flag, real_t def)
{
   const std::string v = ArgStr(argc, argv, flag, "");
   return v.empty() ? def : std::atof(v.c_str());
}
int ArgInt(int argc, char *argv[], const char *flag, int def)
{
   const std::string v = ArgStr(argc, argv, flag, "");
   return v.empty() ? def : std::atoi(v.c_str());
}

const toml::value *DeckTable(const toml::value &root, const char *name)
{
   return root.contains(name) ? &root.at(name) : nullptr;
}
real_t TomlReal(const toml::value *t, const char *key, real_t def)
{
   if (!t || !t->contains(key)) { return def; }
   const auto &v = t->at(key);
   if (v.is_floating()) { return static_cast<real_t>(v.as_floating()); }
   if (v.is_integer()) { return static_cast<real_t>(v.as_integer()); }
   MFEM_ABORT("deck key '" << key << "' must be a number");
   return def;
}
std::string TomlStr(const toml::value *t, const char *key, const std::string &def)
{
   if (!t || !t->contains(key)) { return def; }
   return t->at(key).as_string();
}
bool TomlBool(const toml::value *t, const char *key, bool def)
{
   if (!t || !t->contains(key)) { return def; }
   return t->at(key).as_boolean();
}

/// The station copy average's face data (StationCopyAverage): per owned
/// fault DOF, the centroid and the area of its face, vertices in a canonical
/// order.
Vector StationCopyGeometry(NitscheWaveOperator &op)
{
   ElasticityDomainOperator<ParMesh> &E = op.Elasticity();
   const int nd = E.GetNumFaultDOFs(), nbf = E.GetNbfPerFace();
   // face centroid and area per local DOF, vertices in a canonical order
   Vector cl(4 * nd);
   for (int k = 0; k < nd; k++)
   {
      long long key[3];
      real_t xv[3][3];
      op.Orientation().FaceVertexData(k / nbf, key, xv);
      std::array<std::array<real_t, 3>, 3> v;
      for (int i = 0; i < 3; i++) { for (int c = 0; c < 3; c++) { v[i][c] = xv[i][c]; } }
      std::sort(v.begin(), v.end());
      real_t e1[3], e2[3];
      for (int c = 0; c < 3; c++)
      {
         cl(4 * k + c) = (v[0][c] + v[1][c] + v[2][c]) / 3.0;
         e1[c] = v[1][c] - v[0][c];
         e2[c] = v[2][c] - v[0][c];
      }
      const real_t cx = e1[1] * e2[2] - e1[2] * e2[1], cy = e1[2] * e2[0] - e1[0] * e2[2],
                   cz = e1[0] * e2[1] - e1[1] * e2[0];
      cl(4 * k + 3) = 0.5 * std::sqrt(cx * cx + cy * cy + cz * cz);
   }
   Vector ca;
   E.RestrictToOwnedFault(cl, ca, 4);
   return ca;
}

/// Friction-law factory keyed on [meta] law and [friction.rate_state]
/// state_evolution (plan N-07, Generality; benchmark plan rev 1.3, 4.5).
enum class LawKind { Aging, StrongRateWeakening, SlipWeakening };

struct ClosureBundle
{
   std::unique_ptr<DieterichRuinaFriction> fr;
   std::unique_ptr<IFaultClosure> closure;
   std::string name;
   LawKind kind = LawKind::Aging;
};

/// `t_rs` is the raw [friction.rate_state] table: the strong-rate-weakening
/// branch requires f_w_default typed in the deck, because the parsed struct
/// cannot tell an absent key from the 0.2 the parser fills in (plan 3.2).
ClosureBundle MakeClosure(const spatial::SpatialFrictionConfig &cfg, const toml::value *t_rs)
{
   ClosureBundle B;
   if (cfg.law == spatial::FrictionLawKind::RateState)
   {
      MFEM_VERIFY(cfg.rate_state.has_value(), "[meta] law = rate_state needs "
                  "[friction.rate_state]");
      const auto &rs = *cfg.rate_state;
      B.fr = std::make_unique<DieterichRuinaFriction>(
         DieterichRuinaFriction::Constants{rs.V_0_default, rs.f_0_default,
                                        rs.b_default, rs.Dc_default});
      if (rs.state_evolution == spatial::StateEvolutionKind::SlipLawStrongRateWeakening)
      {
         MFEM_VERIFY(t_rs && t_rs->is_table() && t_rs->contains("f_w_default"),
                     "Nitsche driver: state_evolution = slip-law strong rate weakening "
                     "needs [friction.rate_state] f_w_default in the deck; when the key is "
                     "absent the parser fills 0.2 without a word, so the driver refuses");
         B.closure = std::make_unique<RateStateSRWClosure>(*B.fr, rs.f_w_default);
         std::ostringstream nm;
         nm << "rate_state (slip law, strong rate weakening, f_w = " << rs.f_w_default << ")";
         B.name = nm.str();
         B.kind = LawKind::StrongRateWeakening;
         return B;
      }
      B.closure = std::make_unique<RateStateClosure>(*B.fr, RateStateClosure::StateLaw::Aging);
      B.name = "rate_state (aging law)";
      return B;
   }
   MFEM_VERIFY(cfg.slip_weakening.has_value(), "[meta] law = slip_weakening needs "
               "[friction.slip_weakening]");
   B.closure = std::make_unique<SlipWeakeningClosure>();
   B.name = "slip_weakening (linear)";
   B.kind = LawKind::SlipWeakening;
   return B;
}

/// Unit vector helpers.
void Normalize(real_t v[3])
{
   const real_t l = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
   MFEM_VERIFY(l > 0.0, "zero-length vector in [fault_geometry]");
   for (int i = 0; i < 3; i++) { v[i] /= l; }
}

void MakeDirs(const std::string &d)
{
   if (g_rank == 0) { ::mkdir(d.c_str(), 0755); }
   MPI_Barrier(MPI_COMM_WORLD);
}

/// This process's VmData in MiB (header, "Setup phases"); -1 without /proc.
double VmDataMiB()
{
   std::ifstream f("/proc/self/status");
   for (std::string line; std::getline(f, line);)
   {
      if (line.compare(0, 7, "VmData:") == 0) { return std::atof(line.c_str() + 7) / 1024.0; }
   }
   return -1.0;
}

/// Wall time of each setup phase on this rank (header, "Setup phases").
class SetupTimer
{
public:
   explicit SetupTimer(std::chrono::steady_clock::time_point start) : t0_(start) {}
   /// The phase that ends now.
   void Mark(const char *name)
   {
      const auto t1 = std::chrono::steady_clock::now();
      names_.push_back(name);
      sec_.push_back(std::chrono::duration<double>(t1 - t0_).count());
      vm_.push_back(VmDataMiB());
      t0_ = t1;
   }
   /// Collective: rank 0 prints each phase's max and mean over the ranks.
   void Print() const
   {
      const int n = static_cast<int>(sec_.size());
      std::vector<double> mx(sec_), sm(sec_), vmx(vm_), vsm(vm_), vmn(vm_);
      MPI_Allreduce(MPI_IN_PLACE, mx.data(), n, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, sm.data(), n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, vmx.data(), n, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, vsm.data(), n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, vmn.data(), n, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      if (g_rank != 0) { return; }
      const int nr = Mpi::WorldSize();
      double tmx = 0.0, tsm = 0.0;
      const std::ios::fmtflags fl = std::cout.flags();
      const std::streamsize pr = std::cout.precision();
      std::cout << std::fixed << std::setprecision(2) << "  setup (wall s, max/mean over ranks):";
      for (int i = 0; i < n; i++)
      {
         std::cout << (i ? ", " : " ") << names_[i] << " " << mx[i] << "/" << sm[i] / nr;
         tmx += mx[i];
         tsm += sm[i];
      }
      std::cout << "; sum " << tmx << "/" << tsm / nr << "\n";
      if (n > 0 && *std::min_element(vmn.begin(), vmn.end()) >= 0.0)
      {
         std::cout << std::setprecision(0) << "  setup memory (VmData MiB at the end of each phase,"
                   " max/mean over ranks):";
         for (int i = 0; i < n; i++)
         {
            std::cout << (i ? ", " : " ") << names_[i] << " " << vmx[i] << "/" << vsm[i] / nr;
         }
         std::cout << "\n";
      }
      std::cout.flags(fl);
      std::cout.precision(pr);
   }

private:
   std::chrono::steady_clock::time_point t0_;
   std::vector<std::string> names_;
   std::vector<double> sec_;
   std::vector<double> vm_;   // VmData (MiB) at the end of each phase
};

#ifdef MFEM_USE_HDF5
/// Fault ParaView output, [output] paraview_fault = "hdf5": ONE
/// <output_dir>/fault.vtkhdf per run through the SAFS writer's generic
/// entry point vtkhdf::WriteFaultPackHdf (rank-0 gather, one appended step
/// per frame).  Each owned fault face is split at its edge midpoints into 4
/// sub-triangles; every field is evaluated through the face's own nodal
/// basis at the 3 vertices and 3 edge midpoints (at p = 2 these are exactly
/// the 6 nodes, at p = 1 the midpoints are linear interpolants) and each
/// sub-triangle carries the mean of its 3 corner values as cell data.
/// Every sub-triangle has its own 3 vertices (the writer's layout), so the
/// output stays discontinuous like the DG field.  Owned faces only: a
/// shared fault face is written once.  Sub-triangles are wound so that
/// their normal agrees with the fault frame nhat.  Scope: Chunhui Zhao,
/// 2026-10-14 (P7).
class FaultVTKHDFOutput
{
public:
   /// fast: plan item 1 (vtkhdf::FaultHDFStream: same file, less work per
   /// frame; param_a and param_Dc written once); else the SAFS writer.
   /// spread (fast only): the frames' deflate on the other ranks
   /// (FaultHDFStream::SetSpread).
   FaultVTKHDFOutput(NitscheWaveOperator &op, ParMesh &pmesh, int order,
                     const std::string &dir,
                     ParaViewHDFDataCollection::HDFCompression alg, double param,
                     bool fast = true, bool spread = false)
      : op_(op), E_(op.Elasticity()), dir_(dir), fe_(order, BasisType::GaussLobatto),
        fast_(fast)
   {
      state_.compression_alg = alg;
      state_.compression_param = param;
      stream_.SetCompression(alg, param);
      stream_.SetSpread(spread);
      nbf_ = E_.GetNbfPerFace();
      // the face basis must be the one the fault DOFs are nodal values of
      const IntegrationRule &nir = E_.GetFaceQuadrature()->GetNodalRule();
      MFEM_VERIFY(nir.GetNPoints() == nbf_ && fe_.GetDof() == nbf_,
                  "fault ParaView output: face basis has " << fe_.GetDof()
                  << " functions, fault DOFs per face " << nbf_);
      for (int k = 0; k < nbf_; k++)
      {
         MFEM_VERIFY(std::abs(nir.IntPoint(k).x - fe_.GetNodes().IntPoint(k).x) < 1e-12 &&
                     std::abs(nir.IntPoint(k).y - fe_.GetNodes().IntPoint(k).y) < 1e-12,
                     "fault ParaView output: face node " << k << " differs from "
                     "the fault DOF nodal rule");
      }
      // 3 vertices, then the midpoints of edges 01, 12, 20
      const real_t ref[6][2] = {{0, 0}, {1, 0}, {0, 1}, {0.5, 0}, {0.5, 0.5}, {0, 0.5}};
      phi_.SetSize(6, nbf_);
      Vector sh(nbf_);
      for (int c = 0; c < 6; c++)
      {
         IntegrationPoint ip;
         ip.Set2(ref[c][0], ref[c][1]);
         fe_.CalcShape(ip, sh);
         for (int k = 0; k < nbf_; k++) { phi_(c, k) = sh(k); }
      }
      // geometry and a check of the face map against the DOF positions
      Vector Xl, Fl;
      E_.ExpandOwnedToLocalFault(op.FaultPositions(), Xl, 3);
      E_.ExpandOwnedToLocalFault(op.FaultFrames(), Fl, 9);
      const Array<int> &omap = E_.GetOwnedFaultFaceMap();
      const Array<int> &fint = E_.GetFaultInteriorFaces();
      const Array<int> &fsh = E_.GetFaultSharedFaces();
      const int nint = fint.Size();
      local_faces_.resize(E_.GetNumOwnedFaultFaces());
      real_t pos_err = 0.0;
      int nflip = 0;
      for (int o = 0; o < E_.GetNumOwnedFaultFaces(); o++)
      {
         const int fi = omap[o];
         local_faces_[o] = fi;
         const int mf = (fi < nint) ? fint[fi] : pmesh.GetSharedFace(fsh[fi - nint]);
         ElementTransformation *T = pmesh.GetFaceTransformation(mf);
         real_t P[6][3];
         for (int c = 0; c < 6; c++)
         {
            IntegrationPoint ip;
            ip.Set2(ref[c][0], ref[c][1]);
            T->SetIntPoint(&ip);
            Vector x(3);
            T->Transform(ip, x);
            for (int i = 0; i < 3; i++) { P[c][i] = x(i); }
         }
         for (int k = 0; k < nbf_; k++)
         {
            Vector x(3);
            T->Transform(nir.IntPoint(k), x);
            for (int i = 0; i < 3; i++)
            {
               pos_err = std::max(pos_err, std::abs(x(i) - Xl(3 * (fi * nbf_ + k) + i)));
            }
         }
         real_t e1[3], e2[3];
         for (int i = 0; i < 3; i++) { e1[i] = P[1][i] - P[0][i]; e2[i] = P[2][i] - P[0][i]; }
         const real_t nx[3] = {e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2],
                               e1[0] * e2[1] - e1[1] * e2[0]};
         const real_t *nh = &Fl(9 * (fi * nbf_));
         const bool flip = (nx[0] * nh[0] + nx[1] * nh[1] + nx[2] * nh[2]) < 0.0;
         nflip += flip ? 1 : 0;
         for (int s = 0; s < 4; s++)
         {
            int c0 = kSub[s][0], c1 = kSub[s][1], c2 = kSub[s][2];
            if (flip) { std::swap(c1, c2); }
            const int b = static_cast<int>(pack_.vertices.size());
            for (int c : {c0, c1, c2}) { pack_.vertices.push_back({P[c][0], P[c][1], P[c][2]}); }
            pack_.triangles.push_back({b, b + 1, b + 2});
            corners_.push_back({c0, c1, c2});
         }
      }
      MPI_Allreduce(MPI_IN_PLACE, &pos_err, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX,
                    MPI_COMM_WORLD);
      long long cnt[2] = {E_.GetNumOwnedFaultFaces(), nflip};
      MPI_Allreduce(MPI_IN_PLACE, cnt, 2, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
      MFEM_VERIFY(pos_err < 1e-6, "fault ParaView output: the face map does not reproduce "
                  "the fault DOF positions (max error " << pos_err << " m)");
      if (g_rank == 0)
      {
         std::cout << "  fault ParaView output: " << dir_ << "/fault.vtkhdf, " << cnt[0]
                   << " fault faces x 4 sub-triangles (" << cnt[1] << " faces re-wound to"
                   << " nhat); face map vs DOF positions max error " << pos_err << " m\n";
      }
      pack_.field_names = {"slip_dip", "slip_strike", "slip_rate_dip", "slip_rate_strike",
                           "slip_rate", "traction_dip", "traction_strike", "state_variable",
                           "normal_stress", "param_a", "param_Dc"};
      if (op_.HasPrakashClifton()) { pack_.field_names.push_back("normal_stress_pc"); nc_ = 11; }
      pack_.field_arrays.assign(pack_.field_names.size(), std::vector<double>());
      constant_.assign(pack_.field_names.size(), 0);
      constant_[9] = constant_[10] = 1;   // param_a, param_Dc
   }

   /// Collective.  tau = -(X* + gamma V) and sigma_n_eff as the stations;
   /// at t = 0 (no step yet) the projected background traction.
   void Write(int cycle, real_t t, real_t dt)
   {
      const auto t0 = std::chrono::steady_clock::now();
      const int no = op_.FaultDc().Size();
      const Vector &d = op_.Slip(), &V = op_.SlipRate(), &psi = op_.State();
      const Vector &Fr = op_.FaultFrames();
      Vector Xs = op_.LastShear(), sg = op_.LastSigmaEff();
      const bool have_x = Xs.Size() == 2 * no;
      if (!have_x)
      {
         Vector S;
         op_.BackgroundTraction(op_.GetTime(), S);
         Xs.SetSize(2 * no);
         sg.SetSize(no);
         for (int j = 0; j < no; j++)
         {
            fault_frame::ToLocal(&S(3 * j), 0.0, &Fr(9 * j), &Fr(9 * j + 3), &Fr(9 * j + 6),
                                 Xs(2 * j), Xs(2 * j + 1), sg(j));
         }
      }
      const int nc = nc_;   // sd ss vd vs td ts psi sn a Dc [sigma*]
      Vector own(nc * no), loc;
      for (int j = 0; j < no; j++)
      {
         const real_t *dd = &Fr(9 * j + 3), *ss = &Fr(9 * j + 6);
         real_t *q = &own(nc * j);
         q[0] = fault_frame::Dot(&d(3 * j), dd);
         q[1] = fault_frame::Dot(&d(3 * j), ss);
         q[2] = fault_frame::Dot(&V(3 * j), dd);
         q[3] = fault_frame::Dot(&V(3 * j), ss);
         const real_t g = have_x ? op_.ClosureGamma(j) : 0.0;   // + sigma_F^eta with Kelvin-Voigt
         q[4] = -(Xs(2 * j) + g * q[2]);
         q[5] = -(Xs(2 * j + 1) + g * q[3]);
         q[6] = psi(j);
         q[7] = sg(j);
         q[8] = op_.FaultA()(j);
         q[9] = op_.FaultDc()(j);
         if (nc == 11) { q[10] = op_.PrakashCliftonSigma()(j); }
      }
      E_.ExpandOwnedToLocalFault(own, loc, nc);
      for (auto &a : pack_.field_arrays) { a.clear(); }
      for (std::size_t o = 0; o < local_faces_.size(); o++)
      {
         const int fi = local_faces_[o];
         real_t cv[6][11];
         for (int c = 0; c < 6; c++)
         {
            for (int m = 0; m < nc; m++)
            {
               real_t v = 0.0;
               for (int k = 0; k < nbf_; k++) { v += phi_(c, k) * loc(nc * (fi * nbf_ + k) + m); }
               cv[c][m] = v;
            }
         }
         for (int s = 0; s < 4; s++)
         {
            const auto &cr = corners_[4 * o + s];
            auto mean = [&](int m) { return (cv[cr[0]][m] + cv[cr[1]][m] + cv[cr[2]][m]) / 3.0; };
            auto speed = [&](int c) { return std::hypot(cv[c][2], cv[c][3]); };
            pack_.field_arrays[0].push_back(mean(0));
            pack_.field_arrays[1].push_back(mean(1));
            pack_.field_arrays[2].push_back(mean(2));
            pack_.field_arrays[3].push_back(mean(3));
            pack_.field_arrays[4].push_back((speed(cr[0]) + speed(cr[1]) + speed(cr[2])) / 3.0);
            for (int m = 4; m < nc; m++) { pack_.field_arrays[m + 1].push_back(mean(m)); }
         }
      }
      const auto t1 = std::chrono::steady_clock::now();
      if (fast_) { stream_.Write(dir_, t, pack_, constant_, g_rank, Mpi::WorldSize(), MPI_COMM_WORLD); }
      else
      {
         vtkhdf::WriteFaultPackHdf(state_, dir_, cycle, t, pack_, g_rank, Mpi::WorldSize(),
                                   MPI_COMM_WORLD);
      }
      frames_++;
      const auto t2 = std::chrono::steady_clock::now();
      pack_seconds_ += std::chrono::duration<double>(t1 - t0).count();
      seconds_ += std::chrono::duration<double>(t2 - t0).count();
   }

   int Frames() const { return frames_; }
   double Seconds() const { return seconds_; }
   double PackSeconds() const { return pack_seconds_; }
   bool Fast() const { return fast_; }
   /// fast writer only: gather and HDF5 write seconds on this rank
   double GatherSeconds() const { return stream_.GatherSeconds(); }
   double WriteSeconds() const { return stream_.WriteSeconds(); }
   /// fast writer, spread compression: on, helper ranks, chunk values, rank
   /// 0's scatter-to-last-chunk seconds, this rank's compress2 seconds
   bool SpreadActive() const { return stream_.SpreadActive(); }
   int SpreadHelpers() const { return stream_.SpreadHelpers(); }
   long long ChunkValues() const { return stream_.ChunkValues(); }
   double SpreadSeconds() const { return stream_.SpreadSeconds(); }
   double CompressSeconds() const { return stream_.CompressSeconds(); }

private:
   /// sub-triangles of the 6 corners (0,1,2 vertices; 3,4,5 midpoints of
   /// edges 01, 12, 20), each wound like the face (0,1,2)
   static constexpr int kSub[4][3] = {{0, 3, 5}, {3, 1, 4}, {5, 4, 2}, {3, 4, 5}};

   NitscheWaveOperator &op_;
   ElasticityDomainOperator<ParMesh> &E_;
   std::string dir_;
   H1_TriangleElement fe_;
   int nbf_ = 0;
   int nc_ = 10;   ///< packed fields per DOF (11 with Prakash-Clifton)
   DenseMatrix phi_;
   std::vector<int> local_faces_;
   std::vector<std::array<int, 3>> corners_;
   vtu::LocalFaultPack pack_;
   vtkhdf::FaultHDFState state_;
   vtkhdf::FaultHDFStream stream_;
   std::vector<char> constant_;
   bool fast_ = true;
   int frames_ = 0;
   double seconds_ = 0.0, pack_seconds_ = 0.0;
};
#endif // MFEM_USE_HDF5

/// Per-rank balance report (header, "Rank balance"): collective.  Every rank
/// sends its stepping times and work counts to rank 0, which writes one row
/// per rank to @a csv and logs how the work counts follow the time.
void RankBalanceReport(const NitscheLtsRun &run, const std::vector<int> &cl, int ne,
                       double step_sec, double step1_sec, double loop_sec, double fault_sec,
                       double fault_sides, const std::string &csv)
{
   const NitscheLtsBlocks &B = run.Blocks();
   const int nc = run.NumClusters();
   const int nranks = Mpi::WorldSize();
   // per rank: 11 scalars, then per cluster tets, patch elements, interior
   // face blocks, shared face blocks; then the phase seconds (11 scalars,
   // per cluster interior and finish) and the absorbing elements per cluster
   enum { kStep, kStep1, kWait, kLoop, kFault, kKicks, kTets, kSides, kHalo, kNbr, kRounds, kN };
   enum { pKick, pDrift, pPost, pFault, pWait, pDir, pShared, pFLoads, pAbsorb, pMass, pReduce, pN };
   const int kPh = kN + 4 * nc, kAbs = kPh + pN + 2 * nc;
   const int nv = kAbs + nc;
   std::vector<double> v(nv, 0.0);
   v[kStep] = step_sec;
   v[kStep1] = step1_sec;
   v[kWait] = B.WaitSeconds();
   v[kLoop] = loop_sec;
   v[kFault] = fault_sec;
   v[kKicks] = static_cast<double>(run.PatchUpdates());
   v[kTets] = ne;
   v[kSides] = fault_sides;
   v[kHalo] = static_cast<double>(B.ExchangeValues());
   v[kNbr] = B.NumNeighbourRanks();
   v[kRounds] = static_cast<double>(B.ExchangeRounds());
   for (int e = 0; e < ne; e++) { v[kN + cl[e]] += 1.0; }
   for (int c = 0; c < nc; c++)
   {
      const std::array<long long, 3> k = B.ClusterCounts(c);
      for (int j = 0; j < 3; j++) { v[kN + (j + 1) * nc + c] = static_cast<double>(k[j]); }
   }
   const NitscheLtsRun::PhaseSeconds &ph = run.Phases();
   const double phv[pN] = {ph.kick, ph.drift, ph.post, ph.fault, ph.wait, ph.dirichlet,
                           ph.shared, ph.fault_loads, ph.absorb, ph.mass, ph.reduce};
   for (int j = 0; j < pN; j++) { v[kPh + j] = phv[j]; }
   for (int c = 0; c < nc; c++)
   {
      v[kPh + pN + c] = ph.interior[c];
      v[kPh + pN + nc + c] = ph.finish[c];
      v[kAbs + c] = run.AbsorbingElements(c);
   }
   constexpr int kName = 64;
   char name[MPI_MAX_PROCESSOR_NAME + kName] = {0};
   int len = 0;
   MPI_Get_processor_name(name, &len);
   name[kName - 1] = '\0';
   std::vector<double> all(g_rank == 0 ? static_cast<std::size_t>(nv) * nranks : 0);
   std::vector<char> names(g_rank == 0 ? static_cast<std::size_t>(kName) * nranks : 0);
   MPI_Gather(v.data(), nv, MPI_DOUBLE, all.data(), nv, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   MPI_Gather(name, kName, MPI_CHAR, names.data(), kName, MPI_CHAR, 0, MPI_COMM_WORLD);
   if (g_rank != 0) { return; }

   const auto at = [&](int r, int j) { return all[static_cast<std::size_t>(r) * nv + j]; };
   std::ofstream f(csv);
   f << std::setprecision(10) << "# per-rank balance (fault-writer note rev1, L1): step_sec = "
     "ode.Step after step 1, wait_sec = WaitHalo (step 1 included), nonwait_sec = step_sec - "
     "wait_sec, outside_sec = loop_sec - step_sec - step1_sec, kicks = patch elements kicked, "
     "fault_sides = 2 interior + shared fault faces; ph_* = seconds of each piece of the ticks "
     "after step 1 (NitscheLtsRun phase timers; ph_fin_c holds ph_shared, ph_fault_loads, "
     "ph_absorb and ph_mass, which are summed over the clusters), ph_untimed = step_sec - the sum "
     "of ph_kick, ph_drift, ph_post, ph_fault, ph_int_c, ph_wait, ph_dirichlet, ph_fin_c and "
     "ph_reduce; abs_c = absorbing elements\n";
   f << "rank,node,step_sec,step1_sec,wait_sec,nonwait_sec,loop_sec,outside_sec,fault_sec,"
        "kicks,tets,fault_sides,halo_values,halo_rounds,neighbour_ranks";
   const char *kind[4] = {"tets", "patch", "iface", "sface"};
   for (int j = 0; j < 4; j++)
   {
      for (int c = 0; c < nc; c++) { f << "," << kind[j] << "_c" << c; }
   }
   const char *phname[pN] = {"kick", "drift", "post", "fault", "wait", "dirichlet", "shared",
                             "fault_loads", "absorb", "mass", "reduce"};
   for (int j = 0; j < pN; j++) { f << ",ph_" << phname[j]; }
   f << ",ph_untimed";
   for (int c = 0; c < nc; c++) { f << ",ph_int_c" << c; }
   for (int c = 0; c < nc; c++) { f << ",ph_fin_c" << c; }
   for (int c = 0; c < nc; c++) { f << ",abs_c" << c; }
   f << "\n";
   // the pieces that add up to the step time (ph_fin_c holds shared .. mass)
   const auto timed = [&](int r)
   {
      double sum = at(r, kPh + pKick) + at(r, kPh + pDrift) + at(r, kPh + pPost)
                   + at(r, kPh + pFault) + at(r, kPh + pWait) + at(r, kPh + pDir)
                   + at(r, kPh + pReduce);
      for (int j = kPh + pN; j < kAbs; j++) { sum += at(r, j); }
      return sum;
   };
   std::vector<double> work(nranks), kicks(nranks), tets(nranks), sides(nranks);
   for (int r = 0; r < nranks; r++)
   {
      work[r] = at(r, kStep) - at(r, kWait);
      kicks[r] = at(r, kKicks);
      tets[r] = at(r, kTets);
      sides[r] = at(r, kSides);
      f << r << "," << &names[static_cast<std::size_t>(r) * kName] << "," << at(r, kStep) << ","
        << at(r, kStep1) << "," << at(r, kWait) << "," << work[r] << "," << at(r, kLoop) << ","
        << at(r, kLoop) - at(r, kStep) - at(r, kStep1) << "," << at(r, kFault) << ","
        << static_cast<long long>(kicks[r]) << "," << static_cast<long long>(tets[r]) << ","
        << static_cast<long long>(sides[r]) << "," << static_cast<long long>(at(r, kHalo)) << ","
        << static_cast<long long>(at(r, kRounds)) << "," << static_cast<long long>(at(r, kNbr));
      for (int j = kN; j < kPh; j++) { f << "," << static_cast<long long>(at(r, j)); }
      for (int j = kPh; j < kPh + pN; j++) { f << "," << at(r, j); }
      f << "," << at(r, kStep) - timed(r);
      for (int j = kPh + pN; j < kAbs; j++) { f << "," << at(r, j); }
      for (int j = kAbs; j < nv; j++) { f << "," << static_cast<long long>(at(r, j)); }
      f << "\n";
   }
   // max / mean, the correlation of each count with the non-wait time, and
   // the least-squares fit work = a kicks + b sides + c
   const auto mean = [&](const std::vector<double> &x)
   {
      double s = 0.0;
      for (double y : x) { s += y; }
      return s / x.size();
   };
   const auto maxmean = [&](const std::vector<double> &x)
   {
      return *std::max_element(x.begin(), x.end()) / std::max(mean(x), 1e-300);
   };
   const auto corr = [&](const std::vector<double> &x, const std::vector<double> &y)
   {
      const double mx = mean(x), my = mean(y);
      double sxy = 0.0, sxx = 0.0, syy = 0.0;
      for (std::size_t i = 0; i < x.size(); i++)
      {
         sxy += (x[i] - mx) * (y[i] - my);
         sxx += (x[i] - mx) * (x[i] - mx);
         syy += (y[i] - my) * (y[i] - my);
      }
      return (sxx > 0.0 && syy > 0.0) ? sxy / std::sqrt(sxx * syy) : 0.0;
   };
   std::cout << std::setprecision(4) << "  rank balance (" << csv << "): max/mean non-wait "
             << maxmean(work) << ", tets " << maxmean(tets) << ", patch kicks " << maxmean(kicks)
             << ", fault sides " << maxmean(sides) << "; correlation of non-wait with tets "
             << corr(work, tets) << ", patch kicks " << corr(work, kicks) << ", fault sides "
             << corr(work, sides) << "\n";
   if (nranks >= 8)
   {
      // normal equations of [kicks sides 1] x = work, in units scaled to ~1
      const double sk = std::max(mean(kicks), 1e-300), ss = std::max(mean(sides), 1e-300);
      DenseMatrix A(3, 3);
      Vector rhs(3), x(3);
      A = 0.0;
      rhs = 0.0;
      for (int r = 0; r < nranks; r++)
      {
         const double row[3] = {kicks[r] / sk, sides[r] / ss, 1.0};
         for (int i = 0; i < 3; i++)
         {
            rhs(i) += row[i] * work[r];
            for (int j = 0; j < 3; j++) { A(i, j) += row[i] * row[j]; }
         }
      }
      if (std::abs(A.Det()) > 1e-12 * std::pow(A.MaxMaxNorm(), 3))
      {
         DenseMatrixInverse(A).Mult(rhs, x);
         const double mw = mean(work);
         double res = 0.0, tot = 0.0;
         for (int r = 0; r < nranks; r++)
         {
            const double fit = x(0) * kicks[r] / sk + x(1) * sides[r] / ss + x(2);
            res += (work[r] - fit) * (work[r] - fit);
            tot += (work[r] - mw) * (work[r] - mw);
         }
         std::cout << "  rank balance fit: non-wait = a kicks + b sides + c with a = "
                   << 1e6 * x(0) / sk << " us per patch kick, b = " << 1e3 * x(1) / ss
                   << " ms per fault side, c = " << x(2) << " s; R^2 "
                   << (tot > 0.0 ? 1.0 - res / tot : 0.0) << "\n";
      }
      else { std::cout << "  rank balance fit: singular (the counts do not vary apart)\n"; }
   }
   std::vector<int> order(nranks);
   for (int r = 0; r < nranks; r++) { order[r] = r; }
   std::sort(order.begin(), order.end(), [&](int a, int b) { return work[a] > work[b]; });
   {
      // each piece: mean over the ranks / max / the slowest rank's (by non-wait)
      const int rs = order[0];
      const auto piece = [&](const std::string &name, auto val)
      {
         double sm = 0.0, mx = 0.0;
         for (int r = 0; r < nranks; r++) { sm += val(r); mx = std::max(mx, val(r)); }
         std::cout << " " << name << " " << sm / nranks << "/" << mx << "/" << val(rs);
      };
      std::cout << std::setprecision(4) << "  rank balance phases (s after step 1, mean/max/slowest "
                "rank r" << rs << "):";
      for (int j = 0; j < pN; j++)
      {
         piece(phname[j], [&](int r) { return at(r, kPh + j); });
      }
      for (int c = 0; c < nc; c++)
      {
         piece("int_c" + std::to_string(c), [&](int r) { return at(r, kPh + pN + c); });
      }
      for (int c = 0; c < nc; c++)
      {
         piece("fin_c" + std::to_string(c), [&](int r) { return at(r, kPh + pN + nc + c); });
      }
      piece("untimed", [&](int r) { return at(r, kStep) - timed(r); });
      std::cout << "\n";
   }
   std::cout << "  slowest ranks (non-wait s, tets, patch kicks, fault sides, node):";
   for (int i = 0; i < std::min(nranks, 3); i++)
   {
      const int r = order[i];
      std::cout << " r" << r << " " << work[r] << " " << static_cast<long long>(tets[r]) << " "
                << static_cast<long long>(kicks[r]) << " " << static_cast<long long>(sides[r]) << " "
                << &names[static_cast<std::size_t>(r) * kName] << (i + 1 < std::min(nranks, 3) ? ";" : "");
   }
   std::cout << "\n";
}

} // namespace

int main(int argc, char *argv[])
{
   Mpi::Init(argc, argv);
   Hypre::Init();
   g_rank = Mpi::WorldRank();
   const int nranks = Mpi::WorldSize();
   auto wall0 = std::chrono::steady_clock::now();
   SetupTimer setup_timer(wall0);

   const std::string cfg_path = ArgStr(argc, argv, "--config", "");
   if (cfg_path.empty())
   {
      if (g_rank == 0)
      {
         std::cerr << "usage: seas_nitsche_dyn_driver --config deck.toml [--mesh m.msh]"
                      " [--order p] [--tfinal s] [--output-dir d] [--prefix name]"
                      " [--station-every n] [--dump-every n] [--log-every n]"
                      " [--station-copies average|first] [--station-comm sub|world]"
                      " [--station-flush-every n] [--fast-fault 0|1]"
                      " [--fused-spmv 0|1] [--dense-blocks 0|1] [--fault-weight w]"
                      " [--fast-fault-output 0|1] [--paraview-fault-deflate-level n|none]"
                      " [--fault-frame-compress spread|rank0]"
                      " [--lts-report 0|1]"
                      " [--lts-report-cap n] [--lts-report-ranks a,b,...]"
                      " [--lts off|rate2] [--lts-cap n] [--lts-fault-buffer n]"
                      " [--lts-partition-weight tets|blocks]"
                      " [--lts-same-step 0|1] [--mesh-setup rank0|every-rank]"
                      " [--mesh-parts-out d --mesh-parts-n N] [--mesh-parts d]\n";
      }
      return 1;
   }

   // ------------------------------------------------------------ deck
   spatial::SpatialFrictionConfig cfg = spatial::LoadSpatialFrictionConfig(cfg_path);
   const toml::value root = toml::parse(cfg_path);
   const toml::value *t_num = DeckTable(root, "numerics");
   const toml::value *t_fg = DeckTable(root, "fault_geometry");
   const toml::value *t_diag = DeckTable(root, "diagnostics");
   const toml::value *t_fric = DeckTable(root, "friction");
   const toml::value *t_pc = (t_fric && t_fric->contains("prakash_clifton"))
                             ? &t_fric->at("prakash_clifton") : nullptr;
   const toml::value *t_rs = (t_fric && t_fric->contains("rate_state"))
                             ? &t_fric->at("rate_state") : nullptr;
   const toml::value *t_out = DeckTable(root, "output");
   const toml::value *t_st = (t_out && t_out->contains("stations"))
                             ? &t_out->at("stations") : nullptr;

   const std::string mesh_path = ArgStr(argc, argv, "--mesh", cfg.mesh.path);
   const int p = ArgInt(argc, argv, "--order", cfg.mesh.order);
   const real_t tfinal = ArgReal(argc, argv, "--tfinal", cfg.time.tfinal);
   const std::string out_dir = ArgStr(argc, argv, "--output-dir", cfg.output.output_dir);
   const std::string prefix = ArgStr(argc, argv, "--prefix",
                                     cfg.problem.tag.empty() ? "nitsche" : cfg.problem.tag);
   const int station_every = ArgInt(argc, argv, "--station-every", 1);
   const int dump_every = ArgInt(argc, argv, "--dump-every", 0);
   const int log_every = ArgInt(argc, argv, "--log-every", 100);
   const std::string station_copies = ArgStr(argc, argv, "--station-copies", "average");
   // speed-next O2 (header, "Station output cost")
   const std::string station_comm = ArgStr(argc, argv, "--station-comm", "sub");
   const int station_flush_every = ArgInt(argc, argv, "--station-flush-every", 100);
   // Profiling plan (Chunhui Zhao, 2026-10-17) items 3, 4, 6: precomputed fault
   // face operators, fused K u + K_eta v, fault-weighted partition.  1 = on.
   const int fast_fault = ArgInt(argc, argv, "--fast-fault", 1);
   const int fused_spmv = ArgInt(argc, argv, "--fused-spmv", 1);
   // item 8a (2026-10-18): K u (+ K_eta v) from dense element and face
   // blocks; takes precedence over --fused-spmv, which applies when it is 0
   const int dense_blocks = ArgInt(argc, argv, "--dense-blocks", 1);
   const int fast_fault_output = ArgInt(argc, argv, "--fast-fault-output", 1);
   // speed-next O1 (header, "Fault frame filter"): "" = the deck
   const std::string fault_deflate = ArgStr(argc, argv, "--paraview-fault-deflate-level", "");
   MFEM_VERIFY(fault_deflate.empty() || fault_deflate == "none"
               || (fault_deflate.size() == 1 && fault_deflate[0] >= '0' && fault_deflate[0] <= '9'),
               "--paraview-fault-deflate-level " << fault_deflate << ": expected 0..9 or none");
   // fault-writer note rev1, W1 (header, "Fault frame compression")
   const std::string fault_compress = ArgStr(argc, argv, "--fault-frame-compress", "spread");
   MFEM_VERIFY(fault_compress == "spread" || fault_compress == "rank0",
               "--fault-frame-compress " << fault_compress << ": expected spread or rank0");
   const real_t fault_weight = ArgReal(argc, argv, "--fault-weight", kDefaultFaultWeight);
   const int lts_report = ArgInt(argc, argv, "--lts-report", 0);
   const int lts_report_cap = ArgInt(argc, argv, "--lts-report-cap", 4);
   std::vector<int> lts_report_ranks;
   {
      std::stringstream ss(ArgStr(argc, argv, "--lts-report-ranks", ""));
      for (std::string tok; std::getline(ss, tok, ',');)
      {
         if (!tok.empty()) { lts_report_ranks.push_back(std::atoi(tok.c_str())); }
      }
      if (lts_report_ranks.empty() && nranks > 1) { lts_report_ranks.push_back(nranks); }
      for (int r : lts_report_ranks)
      {
         MFEM_VERIFY(r >= 2, "--lts-report-ranks: rank counts must be >= 2; got " << r);
      }
   }
   // plan item 7: local time stepping (header)
   const std::string lts_mode = ArgStr(argc, argv, "--lts", cfg.numerics.lts);
   MFEM_VERIFY(lts_mode == "off" || lts_mode == "rate2", "--lts / [numerics] lts = '"
               << lts_mode << "': expected \"off\" or \"rate2\"");
   const bool lts_on = (lts_mode == "rate2");
   const int lts_cap = ArgInt(argc, argv, "--lts-cap",
                              (t_num && t_num->contains("lts_nc_cap")) ? cfg.numerics.lts_nc_cap : 4);
   const bool lts_same = ArgInt(argc, argv, "--lts-same-step", 0) != 0;
   const int lts_buffer = ArgInt(argc, argv, "--lts-fault-buffer", 3);
   MFEM_VERIFY(lts_buffer >= 0, "--lts-fault-buffer " << lts_buffer << ": expected >= 0");
   // Partition weights (header)
   const std::string lts_weight = ArgStr(argc, argv, "--lts-partition-weight", "tets");
   MFEM_VERIFY(lts_weight == "tets" || lts_weight == "blocks", "--lts-partition-weight "
               << lts_weight << ": expected tets or blocks");
   MFEM_VERIFY(lts_on || lts_weight == "tets", "--lts-partition-weight blocks needs --lts rate2");
   // Mesh setup (header)
   const std::string mesh_setup = ArgStr(argc, argv, "--mesh-setup", "rank0");
   MFEM_VERIFY(mesh_setup == "rank0" || mesh_setup == "every-rank", "--mesh-setup "
               << mesh_setup << ": expected rank0 or every-rank");
   // Mesh part files (header)
   const std::string parts_out = ArgStr(argc, argv, "--mesh-parts-out", "");
   const int parts_n = ArgInt(argc, argv, "--mesh-parts-n", 0);
   const std::string parts_in = ArgStr(argc, argv, "--mesh-parts", "");
   const bool parts_split = !parts_out.empty();
   const bool parts_load = !parts_in.empty();
   MFEM_VERIFY(!parts_split || (nranks == 1 && parts_n >= 1), "--mesh-parts-out writes the "
               "parts on one rank: run it on 1 rank with --mesh-parts-n N (N >= 1)");
   MFEM_VERIFY(parts_n == 0 || parts_split, "--mesh-parts-n needs --mesh-parts-out");
   MFEM_VERIFY(!(parts_split && parts_load), "--mesh-parts-out and --mesh-parts exclude each other");
   MFEM_VERIFY(!lts_report || (!parts_split && !parts_load), "--lts-report reads the serial "
               "mesh; not with --mesh-parts-out or --mesh-parts");
   // the number of parts of the partition
   const int nparts = parts_split ? parts_n : nranks;
   const bool mesh_rank0 = mesh_setup == "rank0" && nranks > 1 && !parts_load;
   // the ranks that hold the serial mesh
   const bool smesh_here = !parts_load && (!mesh_rank0 || g_rank == 0);
   MFEM_VERIFY(!lts_same || lts_on, "--lts-same-step needs --lts rate2");
   MFEM_VERIFY(station_copies == "average" || station_copies == "first",
               "--station-copies " << station_copies << ": expected average or first");
   MFEM_VERIFY(station_comm == "sub" || station_comm == "world",
               "--station-comm " << station_comm << ": expected sub or world");
   MFEM_VERIFY(station_flush_every >= 1, "--station-flush-every " << station_flush_every
               << ": expected >= 1");
   const real_t cfl = cfg.numerics.cfl;
   const real_t penalty_factor = TomlReal(t_num, "penalty_factor", 1.0);
   const bool energy_on = TomlBool(t_diag, "energy_balance", false);
   const bool pc_on = t_pc != nullptr;
   real_t pc_vstar = 0.0, pc_len = 0.0;
   if (pc_on)
   {
      MFEM_VERIFY(t_pc->is_table(), "[friction.prakash_clifton] must be a table");
      MFEM_VERIFY(t_pc->contains("v_star") && t_pc->contains("length"),
                  "[friction.prakash_clifton] needs both v_star (m/s) and length (m); "
                  "there are no defaults");
      pc_vstar = TomlReal(t_pc, "v_star", 0.0);
      pc_len = TomlReal(t_pc, "length", 0.0);
      MFEM_VERIFY(pc_vstar > 0.0 && pc_len > 0.0, "[friction.prakash_clifton] v_star = "
                  << pc_vstar << ", length = " << pc_len << "; both must be > 0");
   }
   MFEM_VERIFY(p >= 1 && p <= 4, "order p must be in 1..4 (P3 constants); got " << p);
   MFEM_VERIFY(cfl > 0.0 && cfl <= 1.0, "[numerics] cfl (fraction of C(p)) must be in (0,1]");
   MFEM_VERIFY(!lts_on || cfl <= 0.6, "[numerics] cfl = " << cfl << ": local time stepping "
               "needs cfl <= 0.6 (the accepted multirate step limit; undamped resonance bands "
               "start near 0.65)");
   MakeDirs(out_dir);
   // Mesh part files (header): the options the partition depends on, as
   // mesh_parts.txt records them; with LTS off only the first three matter
   auto exact = [](real_t x) {
      std::ostringstream o;
      o << std::setprecision(17) << x;
      return o.str();
   };
   const std::vector<std::pair<std::string, std::string>> parts_opts = {
      {"lts", lts_mode}, {"fault_weight", exact(fault_weight)},
      {"fault_attr", std::to_string(cfg.boundary.fault_attr)},
      {"order", std::to_string(p)}, {"cfl", exact(cfl)}, {"lts_cap", std::to_string(lts_cap)},
      {"lts_fault_buffer", std::to_string(lts_buffer)}, {"dt_max", exact(cfg.time.dt_max)},
      {"lambda", exact(cfg.material_fallback.lambda)}, {"mu", exact(cfg.material_fallback.mu)},
      {"rho", exact(cfg.material_fallback.rho)}};

   // [output] ParaView keys, SAFS semantics (spatial_dyn_driver): the master
   // gate paraview_enabled switches the ParaView collections; the fault
   // surface is written as ONE fault.vtkhdf per run.  Every key this driver
   // cannot honour aborts instead of being ignored.
   const auto &ocfg = cfg.output;
   const bool pv_fault = ocfg.paraview_enabled && ocfg.paraview_fault != "off";
   MFEM_VERIFY(ocfg.paraview_fault == "off" || ocfg.paraview_fault == "hdf5"
               || ocfg.paraview_fault == "vtu",
               "[output] paraview_fault = '" << ocfg.paraview_fault
               << "': expected \"off\", \"hdf5\" or \"vtu\"");
   MFEM_VERIFY(!(pv_fault && ocfg.paraview_fault == "vtu"),
               "[output] paraview_fault = \"vtu\" writes one file per snapshot; the "
               "Nitsche driver writes the fault as ONE file per run: use \"hdf5\" "
               "(VTKHDF) or \"off\"");
   MFEM_VERIFY(ocfg.paraview_free_surface == "off",
               "[output] paraview_free_surface = '" << ocfg.paraview_free_surface
               << "': the Nitsche driver writes no free-surface slice; set it to \"off\" "
               "(the parser's default is \"vtu\")");
   if (ocfg.paraview_enabled)
   {
      MFEM_VERIFY(ocfg.paraview_volume == "off" && ocfg.paraview_bulk == "off",
                  "[output] paraview_volume / paraview_bulk: the Nitsche driver writes "
                  "only the fault surface; set both to \"off\"");
      MFEM_VERIFY(!ocfg.paraview_fault_legacy_ascii,
                  "[output] paraview_fault_legacy_ascii: not supported by the Nitsche "
                  "driver (single VTKHDF file only)");
      MFEM_VERIFY(ocfg.paraview_coseismic_dt <= 0.0 && ocfg.paraview_nucleation_dt <= 0.0
                  && ocfg.paraview_interseismic_dt <= 0.0,
                  "[output] paraview_{coseismic,nucleation,interseismic}_dt: the Nitsche "
                  "driver writes at the fixed cadence paraview_fault_dt (or "
                  "paraview_every_steps); remove the regime-adaptive keys");
   }
#ifndef MFEM_USE_HDF5
   MFEM_VERIFY(!pv_fault, "[output] paraview_fault = \"hdf5\" needs MFEM built with "
               "MFEM_USE_HDF5=YES");
#endif

   if (g_rank == 0)
   {
      std::cout << "=== seas_nitsche_dyn_driver (plan N-07) on " << nranks << " rank(s) ===\n"
                << "  deck " << cfg_path << "\n  " << cfg.description << "\n";
      // velocity-stress numerics keys: listed once, ignored
      const char *vs_keys[] = {"ader_order", "time_integrator", "mixed_flux",
                               "interior_flux", "fault_iterator", "use_pml",
                               "cfl_safety", "cfl_dg_safety", "lts_wiggle",
                               "lts_merge_loss_tol", "lts_sync_dt", "lts_fault_maxdiff"};
      std::string ign;
      for (const char *k : vs_keys) { if (t_num && t_num->contains(k)) { ign += std::string(" ") + k; } }
      if (!ign.empty()) { std::cout << "  ignored [numerics] keys (velocity-stress path):" << ign << "\n"; }
   }

   setup_timer.Mark("deck");

   // ------------------------------------------------------------ material
   MFEM_VERIFY(cfg.material.kind == spatial::MaterialKind::Constant && !cfg.velocity.use_sidecar,
               "Nitsche driver (P7): only the constant material "
               "([velocity] use_sidecar = false, [material] kind constant) is "
               "wired; the sidecar / depth-profile / halfspace selections are P8.");
   MaterialField material = MaterialField::MakeConstant(
      cfg.material_fallback.lambda, cfg.material_fallback.mu, cfg.material_fallback.rho);

   // ------------------------------------------------------------ mesh
   std::unique_ptr<Mesh> smesh;
   long long ne_global = 0;
   std::map<std::string, std::string> parts_info;
   auto info = [&](const std::string &k) {
      const auto it = parts_info.find(k);
      MFEM_VERIFY(it != parts_info.end(), "--mesh-parts: " << parts_in
                  << "/mesh_parts.txt has no " << k);
      return it->second;
   };
   if (parts_load)
   {
      parts_info = ReadMeshPartsInfo(MPI_COMM_WORLD, parts_in);
      MFEM_VERIFY(info("format") == "seas-mesh-parts 1", "--mesh-parts: " << parts_in
                  << "/mesh_parts.txt: format '" << info("format") << "'");
      MFEM_VERIFY(std::stoi(info("parts")) == nranks, "--mesh-parts: the parts in " << parts_in
                  << " are for " << info("parts") << " ranks; this run has " << nranks);
      for (size_t i = 0; i < parts_opts.size(); i++)
      {
         if (i >= 3 && !lts_on) { break; }
         MFEM_VERIFY(info(parts_opts[i].first) == parts_opts[i].second, "--mesh-parts: the "
                     "parts in " << parts_in << " were cut with " << parts_opts[i].first << " = "
                     << info(parts_opts[i].first) << "; this run has " << parts_opts[i].second);
      }
      if (lts_on)
      {
         // parts cut before L2 carry no weight: theirs is the tet weight
         const auto it = parts_info.find("lts_partition_weight");
         const std::string w = it == parts_info.end() ? "tets" : it->second;
         MFEM_VERIFY(w == lts_weight, "--mesh-parts: the parts in " << parts_in << " were cut "
                     "with lts_partition_weight = " << w << "; this run has " << lts_weight);
      }
      if (g_rank == 0)
      {
         const std::string h = FileFnv1a64(mesh_path);
         MFEM_VERIFY(h == info("mesh_fnv1a64"), "--mesh-parts: the parts in " << parts_in
                     << " were cut from " << info("mesh") << " (FNV-1a " << info("mesh_fnv1a64")
                     << "); this run's mesh " << mesh_path << " has FNV-1a " << h);
      }
      ne_global = std::stoll(info("tets"));
   }
   else if (smesh_here)
   {
      smesh = std::make_unique<Mesh>(mesh_path.c_str(), 1, 1);
      MFEM_VERIFY(smesh->Dimension() == 3, "mesh must be 3-D");
      for (int e = 0; e < smesh->GetNE(); e++)
      {
         MFEM_VERIFY(smesh->GetElementBaseGeometry(e) == Geometry::TETRAHEDRON,
                     "Nitsche driver: tetrahedral meshes only; element " << e << " is "
                     << Geometry::Name[smesh->GetElementBaseGeometry(e)]);
      }
      ne_global = smesh->GetNE();
   }
   if (mesh_rank0) { MPI_Bcast(&ne_global, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD); }
   setup_timer.Mark("mesh read");
   if (lts_report)
   {
      MFEM_VERIFY(cfg.boundary.fault_attr > 0, "[boundary] fault_attr is required");
      if (g_rank == 0)
      {
         const auto &mf = cfg.material_fallback;
         const real_t cp = std::sqrt((mf.lambda + 2.0 * mf.mu) / mf.rho);
         const NitscheLtsMesh lm = BuildNitscheLtsMesh(*smesh, p, cfl * kCeig[p], cp,
                                                       cfg.boundary.fault_attr);
         std::cout << "  LTS report: p = " << p << ", cfl " << cfl << " x C(p) " << kCeig[p]
                   << ", c_p " << cp << " m/s\n";
         PrintNitscheLtsReport(std::cout, *smesh, lm, lts_report_cap, lts_report_ranks);
         std::cout << std::flush;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      return 0;
   }
   // Plan item 7: the clusters of the serial mesh, on rank 0 (header)
   NitscheLtsAssignment lts_asg;
   if (lts_on)
   {
      MFEM_VERIFY(cfg.boundary.fault_attr > 0, "[boundary] fault_attr is required");
      const int ne = static_cast<int>(ne_global);
      int hdr[4] = {0, 0, 0, 0};
      double dt0 = 0.0;
      if (parts_load)
      {
         hdr[0] = std::stoi(info("lts_num_clusters"));
         hdr[1] = std::stoi(info("lts_fault_cluster"));
         hdr[2] = std::stoi(info("lts_passes"));
         hdr[3] = std::stoi(info("lts_fault_pairs"));
         dt0 = std::stod(info("lts_dt0"));
         if (g_rank == 0)
         {
            std::cout << std::setprecision(6) << "  LTS (plan item 7): rate2" << (lts_same
                      ? " with --lts-same-step (every cluster at dt0; test)" : "") << ", "
                      << hdr[0] << " clusters (hard cap " << lts_cap << "), dt0 = " << dt0
                      << " s; fault in cluster " << hdr[1] << " (from the mesh parts in "
                      << parts_in << ")\n";
         }
      }
      else if (g_rank == 0)
      {
         const auto &mf = cfg.material_fallback;
         const real_t cp = std::sqrt((mf.lambda + 2.0 * mf.mu) / mf.rho);
         lts_asg = BuildNitscheLtsAssignment(*smesh, p, cfl * kCeig[p], cp, cfg.boundary.fault_attr,
                                             cfg.time.dt_max, lts_cap, lts_buffer);
         hdr[0] = lts_asg.num_clusters;
         hdr[1] = lts_asg.fault_cluster;
         hdr[2] = lts_asg.passes;
         hdr[3] = static_cast<int>(lts_asg.fault_pairs.size());
         dt0 = lts_asg.dt0;
         std::cout << std::setprecision(6) << "  LTS (plan item 7): rate2" << (lts_same ? " with "
                   "--lts-same-step (every cluster at dt0; test)" : "") << ", " << lts_asg.num_clusters
                   << " clusters (hard cap " << lts_cap << "), dt0 = " << dt0 << " s; tets per cluster";
         for (long long c : lts_asg.cells) { std::cout << " " << c; }
         std::cout << "; ideal " << lts_asg.ideal << "x; fault (" << lts_asg.fault_tets
                   << " tets, " << lts_asg.buffer_tets << " buffer tets in " << lts_buffer
                   << " layer(s)) in cluster " << lts_asg.fault_cluster << " after "
                   << lts_asg.passes << " clustering(s)\n";
      }
      MPI_Bcast(hdr, 4, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Bcast(&dt0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      if (!mesh_rank0 && !parts_load)   // rank0, parts: the clusters of the rank's part
      {
         lts_asg.cluster.resize(ne);
         MPI_Bcast(lts_asg.cluster.data(), ne, MPI_INT, 0, MPI_COMM_WORLD);
      }
      lts_asg.num_clusters = hdr[0];
      lts_asg.fault_cluster = hdr[1];
      lts_asg.passes = hdr[2];
      lts_asg.dt0 = dt0;
      MFEM_VERIFY(hdr[3] > 0, "--lts rate2: the mesh has no fault faces (attr "
                  << cfg.boundary.fault_attr << ")");
   }
   setup_timer.Mark("LTS clusters");
   // Plan item 6: METIS weight per element 1 + w (number of its fault faces),
   // so ranks that own fault faces get fewer tets (w = 0: MFEM's default cut).
   // With LTS (plan item 7) the same weight goes to the multi-constraint cut.
   std::vector<int> part;
   if (lts_on && nparts > 1 && !parts_load)
   {
      const int ne = static_cast<int>(ne_global);
      if (smesh_here) { part.resize(ne); }
      if (g_rank == 0)
      {
         const bool blocks = lts_weight == "blocks";
         std::vector<double> cost;
         if (fault_weight > 0.0)
         {
            // integer METIS weights: 100 per tet + 100 w per fault-face side;
            // blocks: w per fault-face side, which BuildLtsBlockWeights scales
            // by 100 and adds on the tet's own cluster to 100 per block
            const double base = blocks ? 0.0 : 100.0, side = blocks ? 1.0 : 100.0;
            cost.assign(ne, base);
            for (const auto &pr : lts_asg.fault_pairs)
            {
               cost[pr.first] += side * fault_weight;
               cost[pr.second] += side * fault_weight;
            }
         }
         const Table &e2e = smesh->ElementToElementTable();
         std::vector<int> xadj(e2e.GetI(), e2e.GetI() + ne + 1);
         std::vector<int> adj(e2e.GetJ(), e2e.GetJ() + e2e.Size_of_connections());
         LtsPartitionOptions popt;
         popt.block_weights = blocks;
         std::vector<int> bw;   // block work per element and cluster (L2)
         try
         {
            BuildLtsAwarePartition(ne, xadj, adj, lts_asg.cluster, lts_asg.num_clusters,
                                   nparts, cost.empty() ? nullptr : &cost,
                                   lts_asg.fault_pairs, popt, part);
            bw = BuildLtsBlockWeights(xadj, adj, lts_asg.cluster, lts_asg.num_clusters, nullptr);
         }
         catch (const std::exception &ex)
         {
            MFEM_ABORT("--lts rate2: the partition failed on " << nparts
                       << " ranks: " << ex.what());
         }
         // the tet balance, weighted as the cut with tets (the cost), not with blocks
         const LtsPartitionImbalance imb = ComputeLtsPartitionImbalance(
            lts_asg.cluster, lts_asg.num_clusters, part, nparts,
            cost.empty() || blocks ? nullptr : &cost, lts_asg.fault_pairs);
         const LtsBlockBalance bb =
            ComputeLtsBlockBalance(bw, lts_asg.num_clusters, part, nparts);
         const LtsBlockBalance tb = ComputeLtsBlockBalance(
            BuildLtsVertexWeights(lts_asg.cluster, lts_asg.num_clusters, nullptr, false),
            lts_asg.num_clusters, part, nparts);
         std::cout << std::setprecision(4) << "  partition: multi-constraint (BuildLtsAwarePartition, "
                   << (blocks ? "block" : "tet") << " weights, fault weight w = " << fault_weight
                   << "); max/mean overall " << imb.overall << ", per cluster";
         for (double x : imb.per_cluster) { std::cout << " " << x; }
         std::cout << ", fault tets " << imb.fault << "; block work per cluster";
         for (double x : bb.per_cluster) { std::cout << " " << x; }
         std::cout << "; tick model " << bb.tick_model << " on block work, " << tb.tick_model
                   << " on tets\n";
      }
      if (!mesh_rank0) { MPI_Bcast(part.data(), ne, MPI_INT, 0, MPI_COMM_WORLD); }
   }
   else if (fault_weight > 0.0 && nparts > 1 && smesh_here)
   {
      // every rank with the serial mesh computes the same cut
      MFEM_VERIFY(cfg.boundary.fault_attr > 0, "[boundary] fault_attr is required");
      const int ne = smesh->GetNE();
      std::vector<double> cost(ne, 1.0);
      long long nface = 0;
      for (int be = 0; be < smesh->GetNBE(); be++)
      {
         if (smesh->GetBdrAttribute(be) != cfg.boundary.fault_attr) { continue; }
         int f, o, e1, e2;
         smesh->GetBdrElementFace(be, &f, &o);
         smesh->GetFaceElements(f, &e1, &e2);
         if (e1 >= 0) { cost[e1] += fault_weight; }
         if (e2 >= 0) { cost[e2] += fault_weight; }
         nface++;
      }
      // integer METIS weights: 100 per tet + 100 w per fault-face side
      for (double &c : cost) { c *= 100.0; }
      const Table &e2e = smesh->ElementToElementTable();
      std::vector<int> xadj(e2e.GetI(), e2e.GetI() + ne + 1);
      std::vector<int> adj(e2e.GetJ(), e2e.GetJ() + e2e.Size_of_connections());
      std::vector<int> cluster(ne, 0);
      LtsPartitionOptions popt;
      popt.imbalance = 1.03;   // METIS's default for k-way, as MFEM's cut
      BuildLtsAwarePartition(ne, xadj, adj, cluster, 1, nparts, &cost, {}, popt, part);
      const LtsPartitionImbalance imb =
         ComputeLtsPartitionImbalance(cluster, 1, part, nparts, &cost, {});
      if (g_rank == 0)
      {
         std::cout << "  partition: METIS with fault weight w = " << fault_weight << " ("
                   << nface << " fault faces); weighted max/mean " << imb.overall << "\n";
      }
   }
   else if ((mesh_rank0 || parts_split) && g_rank == 0)
   {
      // the cut ParMesh(comm, serial mesh) makes on every rank
      int *q = smesh->GeneratePartitioning(nparts, 1);
      part.assign(q, q + smesh->GetNE());
      delete [] q;
   }
   setup_timer.Mark("partition");
   if (parts_split)
   {
      const auto t0 = std::chrono::steady_clock::now();
      MakeDirs(parts_out);
      if (part.empty()) { part.assign(smesh->GetNE(), 0); }   // N = 1
      WriteMeshParts(parts_out, *smesh, part.data(), nparts, lts_on, lts_asg.cluster.data());
      std::ifstream mf(mesh_path, std::ios::binary | std::ios::ate);
      std::vector<std::pair<std::string, std::string>> kv = {
         {"format", "seas-mesh-parts 1"}, {"parts", std::to_string(nparts)}, {"mesh", mesh_path},
         {"mesh_bytes", std::to_string(static_cast<long long>(mf.tellg()))},
         {"mesh_fnv1a64", FileFnv1a64(mesh_path)}, {"tets", std::to_string(ne_global)},
         {"partition", lts_on ? "BuildLtsAwarePartition (LTS)" : fault_weight > 0.0
                       ? "BuildLtsAwarePartition (fault weight)" : "Mesh::GeneratePartitioning (METIS)"}};
      kv.insert(kv.end(), parts_opts.begin(), parts_opts.end());
      if (lts_on)
      {
         kv.push_back({"lts_num_clusters", std::to_string(lts_asg.num_clusters)});
         kv.push_back({"lts_fault_cluster", std::to_string(lts_asg.fault_cluster)});
         kv.push_back({"lts_passes", std::to_string(lts_asg.passes)});
         kv.push_back({"lts_fault_pairs", std::to_string(lts_asg.fault_pairs.size())});
         kv.push_back({"lts_dt0", exact(lts_asg.dt0)});
         kv.push_back({"lts_partition_weight", lts_weight});
      }
      WriteMeshPartsInfo(parts_out, kv);
      std::cout << "  mesh parts: " << nparts << " parts of " << mesh_path << " written to "
                << parts_out << " in " << std::chrono::duration<double>(
                   std::chrono::steady_clock::now() - t0).count() << " s\n" << std::flush;
      return 0;
   }
   // LTS: the clusters of the local elements (the serial elements of this
   // rank, in serial order, as ParMesh builds them; with every-rank checked
   // here by centroid, with rank0 by rank 0 vertex by vertex before it sends
   // the part) and of the face-neighbour elements
   std::vector<int> lts_cl, lts_cl_nbr;
   std::unique_ptr<ParMesh> pmesh_owner;
   if (parts_load)
   {
      pmesh_owner = ReadMeshPart(MPI_COMM_WORLD, parts_in, lts_on, &lts_cl);
   }
   else if (mesh_rank0)
   {
      pmesh_owner = DistributeMeshParts(MPI_COMM_WORLD, smesh.get(), part.data(), lts_on,
                                        lts_asg.cluster.data(), &lts_cl);
      smesh.reset();
      part.clear();
      part.shrink_to_fit();
   }
   else
   {
      pmesh_owner = std::make_unique<ParMesh>(MPI_COMM_WORLD, *smesh,
                                              part.empty() ? nullptr : part.data());
   }
   ParMesh &pmesh = *pmesh_owner;
   if (g_rank == 0)
   {
      if (parts_load)
      {
         std::cout << "  mesh setup: every rank read its part from " << parts_in << " (cut from "
                   << info("mesh") << ", " << info("partition") << ")\n";
      }
      else
      {
         std::cout << "  mesh setup " << (mesh_rank0 ? "rank0: rank 0 read the mesh and sent "
                                          "each rank its part\n"
                                          : nranks == 1 ? "(one rank): ParMesh from the serial "
                                          "mesh\n" : "every-rank: every rank read the mesh\n");
      }
   }
   setup_timer.Mark("ParMesh");
   if (lts_on && !mesh_rank0 && !parts_load)
   {
      auto centroid = [](const Mesh &m, int e, real_t c[3]) {
         Array<int> v;
         m.GetElementVertices(e, v);
         c[0] = c[1] = c[2] = 0.0;
         for (int k = 0; k < v.Size(); k++)
            for (int d = 0; d < 3; d++) { c[d] += m.GetVertex(v[k])[d] / v.Size(); }
      };
      for (int i = 0; i < smesh->GetNE(); i++)
      {
         if (part.empty() || part[i] == g_rank)
         {
            const int k = static_cast<int>(lts_cl.size());
            MFEM_VERIFY(k < pmesh.GetNE(), "--lts: more serial elements on rank " << g_rank
                        << " than its " << pmesh.GetNE() << " local elements");
            real_t a[3], b[3];
            centroid(*smesh, i, a);
            centroid(pmesh, k, b);
            const real_t dd = std::hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
            MFEM_VERIFY(dd <= 1e-9 * (1.0 + std::hypot(a[0], a[1], a[2])), "--lts: local element "
                        << k << " of rank " << g_rank << " is not serial element " << i
                        << " (centroids " << dd << " m apart)");
            lts_cl.push_back(lts_asg.cluster[i]);
         }
      }
      MFEM_VERIFY(static_cast<int>(lts_cl.size()) == pmesh.GetNE(), "--lts: rank " << g_rank
                  << " has " << pmesh.GetNE() << " local elements but " << lts_cl.size()
                  << " serial elements");
   }
   if (lts_on)
   {
      pmesh.ExchangeFaceNbrData();
      DG_FECollection fec0(0, 3);
      ParFiniteElementSpace fes0(&pmesh, &fec0);
      ParGridFunction g0(&fes0);
      Array<int> d0;
      for (int e = 0; e < pmesh.GetNE(); e++)
      {
         fes0.GetElementDofs(e, d0);
         g0(d0[0]) = lts_cl[e];
      }
      g0.ExchangeFaceNbrData();
      lts_cl_nbr.resize(pmesh.GetNFaceNeighborElements());
      for (int j = 0; j < pmesh.GetNFaceNeighborElements(); j++)
      {
         fes0.GetFaceNbrElementVDofs(j, d0);
         lts_cl_nbr[j] = static_cast<int>(std::lround(g0.FaceNbrData()(d0[0])));
      }
      lts_asg.cluster.clear();
      lts_asg.cluster.shrink_to_fit();
   }
   smesh.reset();
   setup_timer.Mark("LTS cluster map");

   // ------------------------------------------------------------ boundary
   MFEM_VERIFY(cfg.boundary.fault_attr > 0, "[boundary] fault_attr is required");
   BoundaryConfig bc;
   bc.fault_attr = cfg.boundary.fault_attr;
   for (int a : cfg.boundary.natural_attrs) { bc.natural_attrs.insert(a); }
   for (int a : cfg.boundary.absorbing_attrs) { bc.absorbing_attrs.insert(a); }

   // ------------------------------------------------------------ fault geometry
   NitscheFaultSetup fs;
   {
      const std::string orient = TomlStr(t_fg, "orientation", "propagated");
      const std::string pin_by = TomlStr(t_fg, "pin_by", "reference");
      const std::string tangent = TomlStr(t_fg, "tangent", "pointwise");
      const std::string junction = TomlStr(t_fg, "junction_policy", "refuse");
      MFEM_VERIFY(orient == "propagated", "[fault_geometry] orientation = '" << orient
                  << "': only \"propagated\" is implemented (N-08)");
      MFEM_VERIFY(tangent == "pointwise", "[fault_geometry] tangent = '" << tangent
                  << "': only \"pointwise\" is implemented (N-09; F8 transport is not)");
      MFEM_VERIFY(junction == "refuse", "[fault_geometry] junction_policy = '" << junction
                  << "': only \"refuse\" is implemented");
      auto &oc = fs.orientation;
      for (int i = 0; i < 3; i++)
      {
         oc.ref_normal[i] = cfg.fault_geometry.ref_normal[i];
         oc.up[i] = cfg.fault_geometry.up[i];
         fs.up[i] = cfg.fault_geometry.up[i];
      }
      if (pin_by == "reference") { oc.pin = FaultOrientationConfig::Pin::Reference; }
      else
      {
         MFEM_ABORT("[fault_geometry] pin_by = '" << pin_by << "': only \"reference\" "
                    "is wired in the driver (Pin::Point needs a pin_point key)");
      }
      oc.c_min = TomlReal(t_fg, "pin_min_margin", oc.c_min);
      oc.fault_attrs = {bc.fault_attr};
   }

   // ------------------------------------------------------------ closure
   ClosureBundle law = MakeClosure(cfg, t_rs);
   fs.closure = law.closure.get();
   fs.sigma_min = cfg.stress.pore_pressure.min_sigma_n_pa;
   MFEM_VERIFY(cfg.stress.pore_pressure.P_p_pa == 0.0
               && cfg.stress.pore_pressure.P_p_grad_pa_per_m == 0.0,
               "Nitsche driver (P7): pore pressure is not wired (P_p must be 0)");

   // ------------------------------------------------------------ prestress (I1)
   MFEM_VERIFY(cfg.stress.kind == spatial::StressSourceKind::ConstantTensor,
               "Nitsche driver: [stress] kind must be \"constant_tensor\" (the global "
               "Cauchy tensor, tension positive; spec I1).  fault_local_prestress is "
               "the slot pattern I2 forbids.");
   const real_t s0[9] = {cfg.stress.sigma_xx_pa, cfg.stress.sigma_xy_pa, cfg.stress.sigma_xz_pa,
                         cfg.stress.sigma_xy_pa, cfg.stress.sigma_yy_pa, cfg.stress.sigma_yz_pa,
                         cfg.stress.sigma_xz_pa, cfg.stress.sigma_yz_pa, cfg.stress.sigma_zz_pa};
   // [[stress_region]] (header, "Prestress regions"): a box and its own global
   // tensor, replacing the background inside the box; the last match wins
   struct StressRegion { real_t lo[3], hi[3], s[9]; };
   std::vector<StressRegion> regions;
   if (root.contains("stress_region"))
   {
      const toml::value &arr = root.at("stress_region");
      MFEM_VERIFY(arr.is_array(), "[[stress_region]] must be an array of tables");
      const char *bk[3][2] = {{"x_min_m", "x_max_m"}, {"y_min_m", "y_max_m"},
                              {"z_min_m", "z_max_m"}};
      const char *sk[6] = {"sigma_xx_pa", "sigma_yy_pa", "sigma_zz_pa",
                           "sigma_xy_pa", "sigma_yz_pa", "sigma_xz_pa"};
      for (size_t r = 0; r < arr.size(); r++)
      {
         const toml::value &e = arr.at(r);
         MFEM_VERIFY(e.is_table(), "[[stress_region]] entry " << r + 1 << " must be a table");
         for (const auto &kv : e.as_table())
         {
            bool known = false;
            for (int c = 0; c < 3; c++) { known = known || kv.first == bk[c][0] || kv.first == bk[c][1]; }
            for (int k = 0; k < 6; k++) { known = known || kv.first == sk[k]; }
            MFEM_VERIFY(known, "[[stress_region]] entry " << r + 1 << ": unknown key '" << kv.first
                        << "' (bounds x/y/z_min_m, x/y/z_max_m; tensor sigma_xx_pa ... sigma_xz_pa)");
         }
         StressRegion R;
         int nb = 0;
         for (int c = 0; c < 3; c++)
         {
            R.lo[c] = -std::numeric_limits<real_t>::infinity();
            R.hi[c] = std::numeric_limits<real_t>::infinity();
            if (e.contains(bk[c][0])) { R.lo[c] = TomlReal(&e, bk[c][0], 0.0); nb++; }
            if (e.contains(bk[c][1])) { R.hi[c] = TomlReal(&e, bk[c][1], 0.0); nb++; }
            MFEM_VERIFY(R.lo[c] < R.hi[c], "[[stress_region]] entry " << r + 1 << ": "
                        << bk[c][0] << " must be < " << bk[c][1]);
         }
         MFEM_VERIFY(nb > 0, "[[stress_region]] entry " << r + 1 << " has no bound; the "
                     "background tensor belongs in [stress]");
         real_t v[6];
         for (int k = 0; k < 6; k++)
         {
            MFEM_VERIFY(e.contains(sk[k]), "[[stress_region]] entry " << r + 1 << " needs "
                        << sk[k] << " (all six components, tension positive; no defaults)");
            v[k] = TomlReal(&e, sk[k], 0.0);
         }
         const real_t full[9] = {v[0], v[3], v[5], v[3], v[1], v[4], v[5], v[4], v[2]};
         for (int i = 0; i < 9; i++) { R.s[i] = full[i]; }
         regions.push_back(R);
      }
      if (g_rank == 0)
      {
         std::cout << "  prestress regions (last match wins over the background):\n";
         for (size_t r = 0; r < regions.size(); r++)
         {
            const StressRegion &R = regions[r];
            std::cout << "    " << r + 1 << ": x [" << R.lo[0] << ", " << R.hi[0] << "] y ["
                      << R.lo[1] << ", " << R.hi[1] << "] z [" << R.lo[2] << ", " << R.hi[2]
                      << "] m: sigma_xx yy zz xy yz xz = " << R.s[0] << " " << R.s[4] << " "
                      << R.s[8] << " " << R.s[1] << " " << R.s[5] << " " << R.s[2] << " Pa\n";
         }
      }
   }
   fs.prestress = [s0, regions](const real_t x[3], real_t s[9]) {
      const real_t *src = s0;
      for (const StressRegion &R : regions)
      {
         if (x[0] >= R.lo[0] && x[0] <= R.hi[0] && x[1] >= R.lo[1] && x[1] <= R.hi[1]
             && x[2] >= R.lo[2] && x[2] <= R.hi[2]) { src = R.s; }
      }
      for (int i = 0; i < 9; i++) { s[i] = src[i]; }
   };

   // pin frame (F1 with the reference normal): nhat, shat, dhat
   real_t pn[3], pu[3], ps[3], pd[3];
   for (int i = 0; i < 3; i++) { pn[i] = cfg.fault_geometry.ref_normal[i]; pu[i] = cfg.fault_geometry.up[i]; }
   Normalize(pn);
   ff::Cross(pu, pn, ps);
   Normalize(ps);
   ff::Cross(ps, pn, pd);

   // ------------------------------------------------------------ nucleation
   if (cfg.nucleation.enabled)
   {
      MFEM_VERIFY(cfg.nucleation.kind == spatial::NucleationKind::GradualOverstressCompactCircular,
                  "Nitsche driver (P7): only [nucleation] kind = "
                  "\"gradual_overstress_compact_circular\" is wired");
      const auto spec = cfg.nucleation.compact_circular;
      MFEM_VERIFY(spec.radius_m > 0.0 && spec.T_nuc_s > 0.0, "compact-circular nucleation: "
                  "radius_m and T_nuc_s must be > 0");
      // strike overstress delta_tau F(r) G(t) as the global tensor
      //   dsigma = -delta_tau F G (shat (x) nhat + nhat (x) shat),
      // so that -(dsigma nhat).shat = +delta_tau F G.  r is the in-plane
      // distance from the centre measured with (dhat, shat), as the SAFS
      // resolver ResolveGradualOverstressCompactCircular measures it.
      std::array<real_t, 3> c = {spec.center_x_m, spec.center_y_m, spec.center_z_m};
      std::array<real_t, 3> sv = {ps[0], ps[1], ps[2]}, nv = {pn[0], pn[1], pn[2]},
                            dv = {pd[0], pd[1], pd[2]};
      fs.nucleation_spatial = [spec, c, sv, nv, dv](const real_t x[3], real_t s[9]) {
         real_t e[3] = {x[0] - c[0], x[1] - c[1], x[2] - c[2]};
         const real_t ed = e[0] * dv[0] + e[1] * dv[1] + e[2] * dv[2];
         const real_t es = e[0] * sv[0] + e[1] * sv[1] + e[2] * sv[2];
         const real_t amp = -spec.delta_tau_pa
                            * spatial::CompactBellFactor(std::sqrt(ed * ed + es * es), spec.radius_m);
         for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) { s[3 * i + j] = amp * (sv[i] * nv[j] + nv[i] * sv[j]); }
      };
      const real_t T = spec.T_nuc_s;
      fs.nucleation_temporal = [T](real_t t) { return spatial::SmoothStep(t, T); };
   }

   // ------------------------------------------------------------ friction parameters
   // split by law (plan 4.5): the slip-weakening branch never touches
   // cfg.rate_state, which an LSW deck does not have
   const bool is_rs = law.kind != LawKind::SlipWeakening;
   spatial::RateStateBlock rs;
   real_t b0 = 0.0, f00 = 0.0, V00 = 0.0;
   if (is_rs)
   {
      rs = *cfg.rate_state;
      rs.eta_auto = false;      // eta is quasi-dynamic radiation damping; the
      rs.eta_default = 0.0;     // Nitsche closure uses gamma = 1/2 sigma_F dt
      b0 = rs.b_default;
      f00 = rs.f_0_default;
      V00 = rs.V_0_default;
   }
   spatial::SpatialFrictionResolver resolver;
   if (is_rs)
   {
      fs.parameters = [&](const real_t x[3], FaultPointState &st) {
         Vector xc(3);
         for (int i = 0; i < 3; i++) { xc(i) = x[i]; }
         Array<int> el(1), at(1);
         el[0] = 0;
         at[0] = bc.fault_attr;
         Vector sn0;
         spatial::RateStatePerDOFParams r = resolver.ResolveRateState(
            rs, xc, el, at, material, pmesh, cfg.stress.pore_pressure, sn0);
         MFEM_VERIFY(r.b(0) == b0 && r.f_0(0) == f00 && r.V_0(0) == V00,
                     "Nitsche driver: per-DOF b / f_0 / V_0 differ from the defaults at ("
                     << x[0] << ", " << x[1] << ", " << x[2] << "); the closure takes "
                     "them as constants (spatially varying b is not wired)");
         st.a = r.a(0);
         st.Dc = r.Dc(0);
         if (law.kind == LawKind::StrongRateWeakening) { st.V_w = r.V_w(0); }
      };
   }
   else
   {
      fs.parameters = [&](const real_t x[3], FaultPointState &st) {
         Vector xc(3);
         for (int i = 0; i < 3; i++) { xc(i) = x[i]; }
         Array<int> at(1);
         at[0] = bc.fault_attr;
         spatial::SlipWeakeningPerDOFParams r =
            resolver.ResolveSlipWeakening(*cfg.slip_weakening, xc, at);
         st.mu_s = r.mu_s(0);
         st.mu_d = r.mu_d(0);
         st.Dc = r.d_c(0);
         st.cohesion = r.cohesion(0);
      };
   }

   setup_timer.Mark("fault setup");

   // ------------------------------------------------------------ operator
   DomainConfig dcfg;
   dcfg.penalty_factor = penalty_factor;
   NitscheWaveOperator op(pmesh, p, material, bc, dcfg, &fs);
   MFEM_VERIFY(op.Orientation().GetReport().junction_edges == 0,
               "[fault_geometry] junction_policy = \"refuse\": the fault has "
               << op.Orientation().GetReport().junction_edges << " junction edges");
   const real_t dt_cfl = ComputeMaxDt(op.Sample(), pmesh, p, cfl * kCeig[p]);
   const real_t dt_gts = std::min(dt_cfl, cfg.time.dt_max);
   // the step of the fault (and of the outputs): the global step, or with LTS
   // the fault cluster's step dt0 2^c_F
   const real_t dt = !lts_on ? dt_gts
                     : (lts_same ? lts_asg.dt0 : std::ldexp(lts_asg.dt0, lts_asg.fault_cluster));
   if (lts_on)
   {
      MFEM_VERIFY(lts_asg.dt0 <= dt_gts * (1.0 + 1e-12), "--lts: dt0 = " << lts_asg.dt0
                  << " s exceeds the global step " << dt_gts << " s");
      if (g_rank == 0)
      {
         std::cout << std::setprecision(12) << "  LTS: dt0 = " << lts_asg.dt0 << " s, global dt = "
                   << dt_gts << " s (relative difference "
                   << (dt_gts - lts_asg.dt0) / dt_gts << "); fault step dt_F = " << dt << " s\n";
      }
   }
   op.SetTimeStep(dt);

   // Kelvin-Voigt damping (Day et al. 2005 eq. A8): [numerics] kv_eta =
   // eta_bar (default 0 = off), kv_scale = "node" (default; beta_e =
   // eta_bar (1/3) (h_e/p) / c_s,e) | "dt" (beta = eta_bar dt), kv_layer_width
   // (m, optional; absent = every element), kv_length = "inscribed" (default; h_e
   // = 6 V / A, the CFL length) | "edge" (h_e = mean edge length).
   {
      MFEM_VERIFY(!(t_num && t_num->contains("kv_eta_over_dt")),
                  "[numerics] kv_eta_over_dt is gone; use kv_eta with kv_scale = \"dt\"");
      const real_t eta_bar = TomlReal(t_num, "kv_eta", 0.0);
      MFEM_VERIFY(eta_bar >= 0.0, "[numerics] kv_eta must be >= 0");
      const std::string scale_s = TomlStr(t_num, "kv_scale", "node");
      MFEM_VERIFY(scale_s == "node" || scale_s == "dt",
                  "[numerics] kv_scale = \"" << scale_s << "\"; use \"node\" or \"dt\"");
      const std::string len_s = TomlStr(t_num, "kv_length", "inscribed");
      MFEM_VERIFY(len_s == "edge" || len_s == "inscribed",
                  "[numerics] kv_length = \"" << len_s << "\"; use \"edge\" or \"inscribed\"");
      if (eta_bar > 0.0)
      {
         const bool has_w = t_num && t_num->contains("kv_layer_width");
         const real_t width = has_w ? TomlReal(t_num, "kv_layer_width", 0.0) : -1.0;
         MFEM_VERIFY(!has_w || width > 0.0, "[numerics] kv_layer_width must be > 0 m");
         op.SetKelvinVoigt(eta_bar, scale_s == "node" ? NitscheWaveOperator::KVScale::Node
                                                      : NitscheWaveOperator::KVScale::Dt,
                           width, len_s == "edge" ? NitscheWaveOperator::KVLength::MeanEdge
                                                  : NitscheWaveOperator::KVLength::Inscribed);
         if (g_rank == 0)
         {
            std::cout << "  Kelvin-Voigt: eta_bar = " << eta_bar << ", scale " << scale_s
                      << (scale_s == "node" ? " (h_e " + len_s + ")" : std::string())
                      << ", beta " << op.KelvinVoigtBetaMin() << " .. "
                      << op.KelvinVoigtBetaMax() << " s (dt = " << dt << " s), "
                      << (has_w ? "width " + std::to_string(width) + " m, " : std::string("whole domain, "))
                      << op.KelvinVoigtLayerElements() << " of " << ne_global << " elements\n";
         }
      }
   }

   setup_timer.Mark("operator");

   // initial fault state: d = 0, psi from the continuous friction condition
   // at the projected prestress, V = V_init along -X_tan (C-07)
   const int no = op.NumOwnedFaultDOFs();
   {
      Vector d(3 * no), V(3 * no), psi(no), S0, sig0(no);
      d = 0.0;
      op.BackgroundTraction(cfg.time.t_initial, S0);
      Vector xc(3);
      Array<int> el(1), at(1);
      el[0] = 0; at[0] = bc.fault_attr;
      for (int j = 0; j < no; j++)
      {
         const real_t *n = &op.FaultFrames()(9 * j), *dd = &op.FaultFrames()(9 * j + 3),
                       *ss = &op.FaultFrames()(9 * j + 6);
         real_t X[2], sig;
         ff::ToLocal(&S0(3 * j), 0.0, n, dd, ss, X[0], X[1], sig);
         const real_t xm = std::hypot(X[0], X[1]);
         for (int c = 0; c < 3; c++) { xc(c) = op.FaultPositions()(3 * j + c); }
         real_t Vi = 0.0;   // slip weakening starts at rest
         if (is_rs)
         {
            Vector sn0;
            spatial::RateStatePerDOFParams r = resolver.ResolveRateState(
               rs, xc, el, at, material, pmesh, cfg.stress.pore_pressure, sn0);
            Vi = r.V_init(0);
         }
         FaultPointState st = op.FaultParams(j);
         law.closure->SeedEquilibrium(xm, Vi, sig, st);
         psi(j) = st.psi;
         real_t Vl[2] = {0.0, 0.0};
         if (xm > 0.0 && Vi > 0.0) { Vl[0] = -Vi * X[0] / xm; Vl[1] = -Vi * X[1] / xm; }
         ff::Embed(Vl[0], Vl[1], dd, ss, &V(3 * j));
         sig0(j) = sig;
      }
      op.SetFaultState(d, V, psi);
      if (pc_on)
      {
         // sigma*^0 = the initial sigma_n per DOF
         op.SetPrakashClifton(pc_vstar, pc_len);
         op.SetPrakashCliftonState(sig0);
      }
      // initial fault state per owned DOF (generic; e.g. to check psi^0
      // against a benchmark's closed form outside the driver)
      std::ofstream ini(out_dir + "/" + prefix + "_initial_r" + std::to_string(g_rank) + ".csv");
      ini << std::setprecision(17) << "# x y z a Dc psi0 tau0_strike(-T.shat) sigma_n0 V0_strike";
      if (law.kind == LawKind::StrongRateWeakening) { ini << " V_w"; }
      if (law.kind == LawKind::SlipWeakening) { ini << " mu_s mu_d cohesion"; }
      ini << "\n";
      if (law.kind == LawKind::SlipWeakening)
      {
         ini << "# slip weakening: a is unused (0) and psi0 is the slip path\n";
      }
      for (int j = 0; j < no; j++)
      {
         const real_t *n = &op.FaultFrames()(9 * j), *dd = &op.FaultFrames()(9 * j + 3),
                       *ss = &op.FaultFrames()(9 * j + 6);
         real_t X[2], sig;
         ff::ToLocal(&S0(3 * j), 0.0, n, dd, ss, X[0], X[1], sig);
         const real_t *x = &op.FaultPositions()(3 * j);
         ini << x[0] << " " << x[1] << " " << x[2] << " " << op.FaultA()(j) << " "
             << op.FaultDc()(j) << " " << psi(j) << " " << -X[1] << " " << sig << " "
             << ff::Dot(&V(3 * j), ss);
         const FaultPointState &P = op.FaultParams(j);
         if (law.kind == LawKind::StrongRateWeakening) { ini << " " << P.V_w; }
         if (law.kind == LawKind::SlipWeakening)
         {
            ini << " " << P.mu_s << " " << P.mu_d << " " << P.cohesion;
         }
         ini << "\n";
      }
      // parameter ranges over the whole fault (D-02 / D-06 read the files)
      if (law.kind != LawKind::Aging)
      {
         real_t mm[12];
         for (int k = 0; k < 6; k++) { mm[2 * k] = 1e300; mm[2 * k + 1] = 1e300; }
         long long nbar = 0;
         for (int j = 0; j < no; j++)
         {
            const FaultPointState &P = op.FaultParams(j);
            const real_t v[6] = {P.a, P.Dc, P.mu_s, P.mu_d, P.cohesion, P.V_w};
            for (int k = 0; k < 6; k++)
            {
               mm[2 * k] = std::min(mm[2 * k], v[k]);
               mm[2 * k + 1] = std::min(mm[2 * k + 1], -v[k]);
            }
            if (P.mu_s >= 1.0e6) { nbar++; }
         }
         MPI_Allreduce(MPI_IN_PLACE, mm, 12, MPITypeMap<real_t>::mpi_type, MPI_MIN, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &nbar, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
         if (g_rank == 0)
         {
            const char *nm[6] = {"a", "Dc", "mu_s", "mu_d", "cohesion", "V_w"};
            std::cout << std::setprecision(10) << "  friction parameters over the owned fault DOFs:";
            for (int k = 0; k < 6; k++)
            {
               if (law.kind == LawKind::SlipWeakening && (k == 0 || k == 5)) { continue; }
               if (law.kind == LawKind::StrongRateWeakening && k >= 2 && k <= 4) { continue; }
               std::cout << " " << nm[k] << " [" << mm[2 * k] << ", " << -mm[2 * k + 1] << "]";
            }
            if (law.kind == LawKind::SlipWeakening)
            {
               std::cout << "; DOFs with a barrier mu_s (>= 1e6): " << nbar;
            }
            std::cout << "\n";
         }
      }
   }

   // ------------------------------------------------------------ nucleation check
   // The projected SCEC along-strike traction -(S.shat) at the fault DOF
   // nearest the deck's hypocentre against the analytic tau_0 + delta_tau
   // F(r) G(t) at that DOF's position (plan P7 decision 9).
   if (cfg.nucleation.enabled)
   {
      const auto &spec = cfg.nucleation.compact_circular;
      const real_t h[3] = {cfg.hypocenter.x_m, cfg.hypocenter.y_m, cfg.hypocenter.z_m};
      struct { double d; int r; } loc = {1e300, g_rank}, glob;
      int jh = -1;
      for (int j = 0; j < no; j++)
      {
         const real_t *x = &op.FaultPositions()(3 * j);
         const real_t dd = std::hypot(x[0] - h[0], x[1] - h[1], x[2] - h[2]);
         if (dd < loc.d) { loc.d = dd; jh = j; }
      }
      MPI_Allreduce(&loc, &glob, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
      if (glob.r == g_rank && jh >= 0)
      {
         const real_t *x = &op.FaultPositions()(3 * jh);
         const real_t *ss = &op.FaultFrames()(9 * jh + 6), *dd = &op.FaultFrames()(9 * jh + 3);
         real_t e[3] = {x[0] - spec.center_x_m, x[1] - spec.center_y_m, x[2] - spec.center_z_m};
         const real_t r = std::hypot(ff::Dot(e, dd), ff::Dot(e, ss));
         const real_t F = spatial::CompactBellFactor(r, spec.radius_m);
         Vector S0;
         op.BackgroundTraction(0.0, S0);
         const real_t tau0 = -ff::Dot(&S0(3 * jh), ss);
         std::cout << std::setprecision(8) << "  nucleation check at the fault DOF nearest the hypocentre ("
                   << x[0] << ", " << x[1] << ", " << x[2] << "), " << glob.d << " m away, r = "
                   << r << " m, F(r) = " << F << ":\n     t [s]   projected -T.shat [MPa]   analytic "
                   "tau_0 + dtau F G [MPa]   difference [Pa]\n";
         for (real_t tc : {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5})
         {
            Vector S;
            op.BackgroundTraction(tc, S);
            const real_t proj = -ff::Dot(&S(3 * jh), ss);
            const real_t ana = tau0 + spec.delta_tau_pa * F * spatial::SmoothStep(tc, spec.T_nuc_s);
            std::cout << "     " << std::setw(5) << tc << "   " << std::setw(20) << proj / 1e6
                      << "   " << std::setw(26) << ana / 1e6 << "   " << proj - ana << "\n";
         }
         std::cout << std::flush;
      }
      MPI_Barrier(MPI_COMM_WORLD);
   }

   setup_timer.Mark("initial state");

   // ------------------------------------------------------------ outputs
   // stations from the station file the deck names (header, "Stations")
   std::unique_ptr<FaultStationWriter> stations;
   std::unique_ptr<StationCopyAverage> station_avg;
   std::vector<int> station_dofs;   // the owned DOFs the stations read (sorted)
   std::string state_label = "log10_theta";
   if (t_st)
   {
      MFEM_VERIFY(t_st->is_table() && t_st->contains("file") && t_st->contains("snap_tolerance_m"),
                  "[output.stations] needs file (the station CSV) and snap_tolerance_m (m); "
                  "there are no defaults");
      const std::string st_path = TomlStr(t_st, "file", "");
      const real_t snap_tol = TomlReal(t_st, "snap_tolerance_m", 0.0);
      MFEM_VERIFY(snap_tol > 0.0, "[output.stations] snap_tolerance_m = " << snap_tol
                  << "; it must be > 0");
      const std::vector<FaultStation> st_list = ReadFaultStationFile(st_path);
      const std::vector<FaultStationSnap> snap =
         SnapFaultStations(st_list, op.FaultPositions(), MPI_COMM_WORLD);
      if (g_rank == 0)
      {
         std::ofstream sl(out_dir + "/" + prefix + "_station_snap.txt");
         sl << std::setprecision(12) << "# station file " << st_path << "; snap_tolerance_m "
            << snap_tol << "\n# name x y z (requested, m)  X Y Z (nearest fault node, m)  distance (m)\n";
         std::cout << "  stations from " << st_path << " (" << st_list.size()
                   << "), nearest fault node in 3-D, snap_tolerance_m " << snap_tol << ":\n";
         for (size_t s = 0; s < st_list.size(); s++)
         {
            const FaultStation &q = st_list[s];
            sl << q.name << " " << q.x[0] << " " << q.x[1] << " " << q.x[2] << "  "
               << snap[s].node[0] << " " << snap[s].node[1] << " " << snap[s].node[2] << "  "
               << snap[s].distance << "\n";
            std::cout << "    " << q.name << ": requested (" << q.x[0] << ", " << q.x[1] << ", "
                      << q.x[2] << "), snapped (" << snap[s].node[0] << ", " << snap[s].node[1]
                      << ", " << snap[s].node[2] << "), distance " << snap[s].distance << " m\n";
         }
         std::cout << std::flush;
      }
      for (size_t s = 0; s < st_list.size(); s++)
      {
         MFEM_VERIFY(snap[s].distance <= snap_tol, "station " << st_list[s].name << " ("
                     << st_list[s].x[0] << ", " << st_list[s].x[1] << ", " << st_list[s].x[2]
                     << "): the nearest fault node is " << snap[s].distance << " m away, more "
                     "than [output.stations] snap_tolerance_m = " << snap_tol);
      }
      FaultStationWriter::StateFn state_fn;
      if (law.kind == LawKind::Aging)
      {
         // theta = (Dc / V_0) exp((psi - f_0) / b), constants from the deck
         state_fn = [V00, f00, b0](const DOFData &d) {
            const real_t theta = (d.Dc / V00) * std::exp((d.psi - f00) / b0);
            return std::log10(std::max(theta, real_t(1e-300)));
         };
      }
      else
      {
         state_label = law.kind == LawKind::SlipWeakening ? "slip_path" : "psi";
         state_fn = [](const DOFData &d) { return d.psi; };
      }
      stations = std::make_unique<FaultStationWriter>();
      stations->Open(out_dir, prefix, st_list, snap, state_label, state_fn, station_flush_every);
      if (station_copies == "average")
      {
         station_avg = std::make_unique<StationCopyAverage>(
            StationCopyGeometry(op), op.FaultPositions(), st_list, out_dir, prefix,
            station_comm == "sub", station_flush_every);
      }
      station_dofs = stations->Dofs();
      if (station_avg)
      {
         const std::vector<int> d = station_avg->Dofs();
         station_dofs.insert(station_dofs.end(), d.begin(), d.end());
      }
      std::sort(station_dofs.begin(), station_dofs.end());
      station_dofs.erase(std::unique(station_dofs.begin(), station_dofs.end()), station_dofs.end());
      if (station_flush_every > 1) { InstallStationSignalHandlers(); }
      int nmem = station_avg ? station_avg->Members() : 0, nown = static_cast<int>(stations->Dofs().size());
      int root = (station_avg && station_avg->Root()) ? g_rank : -1;
      MPI_Allreduce(MPI_IN_PLACE, &nown, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &root, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &nmem, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if (g_rank == 0)
      {
         std::cout << "  station output: rows flushed every " << station_flush_every << " row(s)"
                   << (station_flush_every > 1 ? " and on SIGTERM/SIGINT/SIGUSR1" : "") << "; ";
         if (station_avg)
         {
            std::cout << "copy average over " << (station_comm == "sub" ? "the " : "all ") << nmem
                      << " rank(s) " << (station_comm == "sub" ? "holding copies" : "(world)")
                      << ", envelope files on rank " << root;
         }
         else { std::cout << "nearest copy only"; }
         std::cout << "; " << nown << " station file(s)\n";
      }
      if (g_rank == 0)
      {
         std::ofstream h(out_dir + "/" + prefix + "_station_conventions.txt");
         h << "# Nitsche dynamic-rupture station conventions (spec I3)\n"
           << "# frame: nhat = pin reference normal (" << pn[0] << ", " << pn[1] << ", " << pn[2]
           << "), up = (" << pu[0] << ", " << pu[1] << ", " << pu[2] << ")\n"
           << "#        shat = up x nhat / |up x nhat| = (" << ps[0] << ", " << ps[1] << ", " << ps[2]
           << "), dhat = shat x nhat = (" << pd[0] << ", " << pd[1] << ", " << pd[2] << ")\n"
           << "# pin: orientation propagated, pinned by the reference normal; plus side = side nhat points into\n"
           << "# columns: t slip1 slip2 V1 V2 tau1 tau2 sigma_n " << state_label << "\n"
           << "#   slip1 = delta.dhat (SCEC v-slip)      slip2 = delta.shat (SCEC h-slip)\n"
           << "#   V1    = V.dhat     (SCEC v-slip-rate) V2    = V.shat     (SCEC h-slip-rate)\n"
           << "#   tau1  = -T.dhat    (SCEC v-shear-stress, Pa)  tau2 = -T.shat (SCEC h-shear-stress, Pa)\n"
           << "#   sigma_n = -T.nhat (compression positive, Pa); T = total traction on the face with normal nhat\n"
           << "#   delta = u(+) - u(-), V = d delta/dt (SCEC slip definition)\n"
           << "#   " << state_label << (law.kind == LawKind::Aging
                                        ? ": log10 of theta = (Dc/V_0) exp((psi - f_0)/b), deck constants\n"
                                        : law.kind == LawKind::SlipWeakening
                                        ? ": accumulated slip path (m)\n" : ": state variable psi\n");
      }
   }
   else if (g_rank == 0)
   {
      std::cout << "  no [output.stations] table: no station files\n";
   }
   std::vector<DOFData> dofdata(no);
   auto fill_dofdata = [&]() {
      const Vector &d = op.Slip(), &V = op.SlipRate(), &psi = op.State();
      const bool have_x = op.LastShear().Size() == 2 * no;
      Vector Xb, sgb;
      if (!have_x)
      {
         // before the first step: no bulk traction yet, the total traction
         // is the projected background (prestress + nucleation at t)
         Vector S;
         op.BackgroundTraction(op.GetTime(), S);
         Xb.SetSize(2 * no);
         sgb.SetSize(no);
         for (int j = 0; j < no; j++)
         {
            const real_t *n = &op.FaultFrames()(9 * j), *dd = &op.FaultFrames()(9 * j + 3),
                          *ss = &op.FaultFrames()(9 * j + 6);
            ff::ToLocal(&S(3 * j), 0.0, n, dd, ss, Xb(2 * j), Xb(2 * j + 1), sgb(j));
         }
      }
      const Vector &Xs = have_x ? op.LastShear() : Xb, &sg = have_x ? op.LastSigmaEff() : sgb;
      for (int j : station_dofs)
      {
         const real_t *dd = &op.FaultFrames()(9 * j + 3), *ss = &op.FaultFrames()(9 * j + 6);
         DOFData &q = dofdata[j];
         q.slip1 = ff::Dot(&d(3 * j), dd);
         q.slip2 = ff::Dot(&d(3 * j), ss);
         q.V1 = ff::Dot(&V(3 * j), dd);
         q.V2 = ff::Dot(&V(3 * j), ss);
         const real_t g = have_x ? op.ClosureGamma(j) : 0.0;   // + sigma_F^eta with Kelvin-Voigt
         // endpoint traction tangential components in the frame: X* + gamma V
         const real_t X1 = Xs(2 * j) + g * q.V1;
         const real_t X2 = Xs(2 * j + 1) + g * q.V2;
         q.tau1_corr = -X1;
         q.tau2_corr = -X2;
         q.sigma_n_corr = sg(j);
         q.psi = psi(j);
         q.Dc = op.FaultDc()(j);
         q.slip_rate = std::hypot(q.V1, q.V2);
      }
   };
   std::unique_ptr<std::ofstream> dump;
   Vector cf;
   if (dump_every > 0)
   {
      dump = std::make_unique<std::ofstream>(out_dir + "/" + prefix + "_fault_r"
                                             + std::to_string(g_rank) + ".csv");
      *dump << std::setprecision(12);
      // face centroid per owned DOF (keys DG copies uniquely)
      ElasticityDomainOperator<ParMesh> &E = op.Elasticity();
      const int nd = E.GetNumFaultDOFs(), nbf = E.GetNbfPerFace();
      Vector cl(3 * nd);
      for (int k = 0; k < nd; k++)
      {
         long long key[3];
         real_t xv[3][3];
         op.Orientation().FaceVertexData(k / nbf, key, xv);
         for (int c = 0; c < 3; c++) { cl(3 * k + c) = (xv[0][c] + xv[1][c] + xv[2][c]) / 3.0; }
      }
      E.RestrictToOwnedFault(cl, cf, 3);
   }
   std::vector<real_t> trupt(no, -1.0);

   // fault ParaView output (one fault.vtkhdf per run)
#ifdef MFEM_USE_HDF5
   std::unique_ptr<FaultVTKHDFOutput> pv;
#endif
   int pv_every = 0;
   bool pv_capped = false;
   if (pv_fault)
   {
#ifdef MFEM_USE_HDF5
      pv_every = ocfg.paraview_every_steps > 0
                 ? ocfg.paraview_every_steps
                 : std::max(1, static_cast<int>(std::lround(ocfg.paraview_fault_dt / dt)));
      auto alg = ParaViewHDFDataCollection::HDFCompression::Deflate;
      double param = 6.0;   // the writer's default (lossless deflate level 6)
      MFEM_VERIFY(fault_deflate.empty() || ocfg.paraview_fault_zfp_tol <= 0.0,
                  "--paraview-fault-deflate-level " << fault_deflate << " with [output] "
                  "paraview_fault_zfp_tol = " << ocfg.paraview_fault_zfp_tol << ": choose one filter");
      if (ocfg.paraview_fault_zfp_tol > 0.0)
      {
         ParaViewOutput<ParMesh>::ProbeH5ZZfpPluginOrAbort();
         alg = ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy;
         param = ocfg.paraview_fault_zfp_tol;
      }
      else if (fault_deflate == "none") { param = -1.0; }   // no filter (both writers)
      else if (!fault_deflate.empty()) { param = fault_deflate[0] - '0'; }
      else if (ocfg.paraview_fault_deflate_level >= 0)
      {
         param = ocfg.paraview_fault_deflate_level;
      }
      pv = std::make_unique<FaultVTKHDFOutput>(op, pmesh, p, out_dir, alg, param,
                                               fast_fault_output != 0,
                                               fault_compress == "spread");
      if (g_rank == 0)
      {
         std::cout << "  fault ParaView output every " << pv_every << " steps (" << pv_every * dt
                   << " s), ";
         if (param < 0.0) { std::cout << "no filter"; }
         else
         {
            std::cout << (alg == ParaViewHDFDataCollection::HDFCompression::ZfpAccuracy
                          ? "ZFP accuracy " : "deflate level ") << param;
         }
         std::cout << ", at most " << ocfg.max_snapshots << " snapshots";
         if (fast_fault_output != 0 && param >= 0.0
             && alg == ParaViewHDFDataCollection::HDFCompression::Deflate)
         {
            std::cout << ", --fault-frame-compress " << fault_compress;
         }
         std::cout << "\n";
      }
#endif
   }

   std::unique_ptr<NitscheEnergyMonitor> mon;

   long nsteps = std::lround(std::ceil((tfinal - cfg.time.t_initial) / dt - 1e-9));
   // LTS: a whole number of sync intervals, nsteps fault steps
   int lts_ticks_sync = 1, lts_ticks_fault = 1;
   long lts_nsync = 0;
   if (lts_on)
   {
      lts_ticks_sync = lts_same ? 1 : (1 << (lts_asg.num_clusters - 1));
      lts_ticks_fault = lts_same ? 1 : (1 << lts_asg.fault_cluster);
      const real_t tsync = lts_ticks_sync * lts_asg.dt0;
      lts_nsync = std::lround(std::ceil((tfinal - cfg.time.t_initial) / tsync - 1e-9));
      nsteps = lts_nsync * (lts_ticks_sync / lts_ticks_fault);
   }
   if (g_rank == 0)
   {
      long long nloc = op.Height();
      std::cout << "  mesh " << mesh_path << ": " << ne_global << " tets; p = " << p
                << "; law " << law.name << "; fault faces " << op.Orientation().GetReport().faces
                << ", shared fault faces " << op.NumSharedFaultFacesGlobal() << "\n"
                << "  dt = " << dt << " s (CFL walk: cfl " << cfl << " x C(p) " << kCeig[p]
                << " -> " << dt_cfl << "; dt_max " << cfg.time.dt_max << "), " << nsteps
                << " steps to t = " << tfinal << " s; penalty_factor " << penalty_factor
                << "; energy balance " << (energy_on ? "on" : "off") << "\n";
      if (lts_on)
      {
         std::cout << "  LTS: " << lts_nsync << " sync intervals of " << lts_ticks_sync
                   << " ticks (" << lts_ticks_sync * lts_asg.dt0 << " s), the fault every "
                   << lts_ticks_fault << " tick(s): " << nsteps << " fault steps to t = "
                   << cfg.time.t_initial + lts_nsync * lts_ticks_sync * lts_asg.dt0 << " s\n";
      }
      if (pc_on)
      {
         std::cout << "  Prakash-Clifton regularization of the law's normal stress: v_star = "
                   << pc_vstar << " m/s, length = " << pc_len << " m (L/v_star = "
                   << pc_len / pc_vstar << " s); sigma*^0 = initial sigma_n\n";
      }
      (void)nloc;
   }

   // ------------------------------------------------------------ loop
   const int n = op.Height();
   Vector U(n), Ud(n);
   U = 0.0;
   Ud = 0.0;
   real_t t = cfg.time.t_initial;
   op.SetTime(t);
   // energy records: er is the latest, er_prev the one before it; W_fric is
   // the trapezoid of P_fric over the evaluation instants (every step, or
   // every sync point with LTS), the scale for sum rho in the closing line
   NitscheEnergyRecord er, er_prev;
   real_t W_fric = 0.0;
   if (energy_on)
   {
      mon = std::make_unique<NitscheEnergyMonitor>(op);
      Vector A0(n);
      if (lts_on) { A0 = 0.0; }   // LTS: the plain energy (header)
      else { op.Mult(U, Ud, A0); }
      er_prev = mon->Evaluate(U, Ud, A0, t);
      mon->Start(er_prev);
   }
   // --kv-station-sigma 1: write, per station, the Kelvin-Voigt fault penalty sigma_F^eta (Pa s/m)
   // averaged over the DG copies exactly as the station traction is (columns tau1 = tau2 =
   // sigma_F^eta, sigma_n = sigma_F dt/2), then stop.  Offline correction of runs whose station
   // traction missed the damping part (before ClosureGamma, 2026-10-15):
   //    tau_true = tau_written - sigma_F^eta V   (per component).
   if (ArgInt(argc, argv, "--kv-station-sigma", 0) != 0)
   {
      MFEM_VERIFY(stations && op.HasKelvinVoigt(),
                  "--kv-station-sigma needs the station writer and [numerics] kv_eta > 0");
      fill_dofdata();
      for (int j = 0; j < no; j++)
      {
         DOFData &q = dofdata[j];
         q.tau1_corr = q.tau2_corr = op.FaultSigmaFKV()(j);
         q.sigma_n_corr = 0.5 * op.FaultSigmaF()(j) * dt;
      }
      if (station_avg) { station_avg->Apply(t, dofdata); }
      stations->WriteStep(t, dofdata);
      if (g_rank == 0)
      {
         std::cout << "  --kv-station-sigma: station files hold sigma_F^eta (tau columns) and "
                      "sigma_F dt/2 (sigma_n column); stopping\n";
      }
      return 0;
   }
   setup_timer.Mark("outputs");
   op.SetFastFaultKernels(fast_fault != 0);
   op.SetFusedStiffness(fused_spmv != 0);
   op.SetDenseBlocks(dense_blocks != 0);
   CentralDifferenceSolver ode;
   std::unique_ptr<NitscheLtsRun> lts_run;
   if (lts_on)
   {
      const auto tl0 = std::chrono::steady_clock::now();
      lts_run = std::make_unique<NitscheLtsRun>(op, lts_cl, lts_cl_nbr, lts_asg.num_clusters,
                                                lts_asg.fault_cluster, lts_asg.dt0,
                                                penalty_factor, lts_same);
      lts_run->Init(U, Ud, t);
      const double tl = std::chrono::duration<double>(std::chrono::steady_clock::now() - tl0).count();
      long long b = lts_run->Blocks().Bytes();
      double tlm = tl;
      MPI_Allreduce(MPI_IN_PLACE, &b, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &tlm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      if (g_rank == 0)
      {
         std::cout << std::setprecision(6) << "  LTS: cluster blocks " << b / 1e6
                   << " MB over all ranks, built in " << tlm << " s (with the first forces)\n";
      }
   }
   else { ode.Init(op); }
   setup_timer.Mark(lts_on ? "LTS blocks" : "solver init");
   setup_timer.Print();
   Vector lts_zero_acc;
   if (lts_on && mon) { lts_zero_acc.SetSize(n); lts_zero_acc = 0.0; }
   if (stations)
   {
      fill_dofdata();
      if (station_avg) { station_avg->Apply(t, dofdata); }
      stations->WriteStep(t, dofdata);
   }
   real_t max_ratio_run = 0.0, max_V_run = 0.0, min_sig_run = 1e300, max_sig_run = -1e300;
   real_t min_pc_run = 1e300, max_pc_run = -1e300;
   // this rank's time inside ode.Step; step 1 is kept apart because it also
   // assembles K, K_eta, the fused pattern and the fault face operators
   double step_sec = 0.0, step1_sec = 0.0, fault_sec1 = 0.0;
   long lts_next_log = log_every;   // LTS: log at the first sync point from here
   bool station_signal_seen = false;
   const auto loop0 = std::chrono::steady_clock::now();
   for (long s = 1; s <= nsteps; s++)
   {
      real_t h = dt;
      const auto ts0 = std::chrono::steady_clock::now();
      bool sync = true;   // v synchronised (always without LTS)
      if (lts_run)
      {
         sync = lts_run->FaultStep(U, Ud, t);
         op.SetTime(t);
      }
      else { ode.Step(U, Ud, t, h); }
      const double ds = std::chrono::duration<double>(std::chrono::steady_clock::now() - ts0).count();
      if (s == 1)
      {
         step1_sec = ds;
         fault_sec1 = op.FaultSeconds();
         if (lts_run) { lts_run->ResetPhaseSeconds(); }
      }
      else { step_sec += ds; }
      // with LTS the record is reduced over the ranks once per sync interval
      // (extremes over its fault steps); between sync points it is this
      // rank's, so the run extremes and the log lines use sync points only
      const NitscheStepRecord &rec = op.LastStep();
      if (sync)
      {
         max_ratio_run = std::max(max_ratio_run, rec.max_V_dt_over_Dc);
         max_V_run = std::max(max_V_run, rec.max_V);
         min_sig_run = std::min(min_sig_run, rec.min_sigma_eff);
         max_sig_run = std::max(max_sig_run, rec.max_sigma_eff);
         if (pc_on)
         {
            min_pc_run = std::min(min_pc_run, rec.min_sigma_pc);
            max_pc_run = std::max(max_pc_run, rec.max_sigma_pc);
         }
      }
      MFEM_VERIFY(std::isfinite(rec.max_V), "non-finite slip rate at step " << s);
      const Vector &V = op.SlipRate();
      for (int j = 0; j < no; j++)
      {
         if (trupt[j] < 0.0 && std::sqrt(V(3 * j) * V(3 * j) + V(3 * j + 1) * V(3 * j + 1)
                                         + V(3 * j + 2) * V(3 * j + 2)) > 1e-3) { trupt[j] = t; }
      }
      real_t rho = 0.0;
      if (mon && sync)
      {
         er = mon->Evaluate(U, Ud, lts_run ? lts_zero_acc : op.LastAcceleration(), t);
         rho = mon->Step(er);
         W_fric += 0.5 * (er.t - er_prev.t) * (er.P_fric + er_prev.P_fric);
         er_prev = er;
      }
      if (stations && s % station_every == 0)
      {
         fill_dofdata();
         if (station_avg) { station_avg->Apply(t, dofdata); }
         stations->WriteStep(t, dofdata);
      }
      if (stations && !station_signal_seen && StationSignal() != 0)
      {
         // header, "Station output cost": this rank's files, from here every row
         station_signal_seen = true;
         stations->SetFlushEvery(1);
         stations->Flush();
         if (station_avg) { station_avg->SetFlushEvery(1); station_avg->Flush(); }
         if (g_rank == 0)
         {
            std::cout << "  signal " << StationSignal() << " at step " << s << " (t = " << t
                      << " s): station files flushed; every row is flushed from here\n" << std::flush;
         }
      }
      if (dump && s % dump_every == 0)
      {
         const Vector &d = op.Slip(), &psi = op.State();
         const Vector &Xs = op.LastShear(), &sg = op.LastSigmaEff();
         const Vector &X = op.FaultPositions();
         for (int j = 0; j < no; j++)
         {
            *dump << s << "," << t << "," << X(3 * j) << "," << X(3 * j + 1) << "," << X(3 * j + 2)
                  << "," << cf(3 * j) << "," << cf(3 * j + 1) << "," << cf(3 * j + 2)
                  << "," << d(3 * j) << "," << d(3 * j + 1) << "," << d(3 * j + 2)
                  << "," << V(3 * j) << "," << V(3 * j + 1) << "," << V(3 * j + 2)
                  << "," << psi(j) << "," << Xs(2 * j) << "," << Xs(2 * j + 1) << "," << sg(j);
            if (pc_on) { *dump << "," << op.PrakashCliftonSigma()(j); }
            *dump << "\n";
         }
      }
#ifdef MFEM_USE_HDF5
      if (pv && s % pv_every == 0)
      {
         // Frames() is the same on every rank, so the collective Write is uniform
         if (pv->Frames() < ocfg.max_snapshots) { pv->Write(static_cast<int>(s), t, dt); }
         else if (!pv_capped)
         {
            pv_capped = true;
            if (g_rank == 0)
            {
               std::cout << "  fault ParaView output: max_snapshots = " << ocfg.max_snapshots
                         << " reached at t = " << t << " s; no further frames\n";
            }
         }
      }
#endif
      const bool log_now = lts_run ? (sync && (s >= lts_next_log || s == nsteps))
                                   : (s % log_every == 0 || s == nsteps);
      if (lts_run && log_now) { lts_next_log = (s / log_every + 1) * log_every; }
      if (g_rank == 0 && log_now)
      {
         const real_t wall = std::chrono::duration<real_t>(std::chrono::steady_clock::now() - wall0).count();
         std::cout << std::setprecision(6) << "  step " << s << " t = " << t << "  max|V| = "
                   << rec.max_V << "  max|V|dt/Dc = " << rec.max_V_dt_over_Dc
                   << "  sigma_n_eff in [" << rec.min_sigma_eff << ", " << rec.max_sigma_eff << "]";
         if (pc_on)
         {
            std::cout << "  sigma_pc in [" << rec.min_sigma_pc << ", " << rec.max_sigma_pc << "]";
         }
         if (mon)
         {
            std::cout << "  E_h = " << er.E_h << "  sum rho = " << mon->SumRho()
                      << "  W_fric = " << W_fric;
         }
         std::cout << "  wall " << wall << " s\n" << std::flush;
         (void)rho;
      }
   }
   // this rank's whole step loop (stations, frames, logs included)
   const double loop_sec = std::chrono::duration<double>(std::chrono::steady_clock::now()
                                                         - loop0).count();
   if (stations) { stations->Close(); }
   if (mon && g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  energy balance: sum rho = " << mon->SumRho()
                << " J, sum |rho| = " << mon->SumAbsRho() << " J, friction work W_fric = "
                << W_fric << " J (trapezoid over " << (lts_run ? "sync points" : "steps")
                << "), sum rho / |W_fric| = "
                << (W_fric != 0.0 ? mon->SumRho() / std::abs(W_fric) : 0.0)
                << ", sum |rho| / |W_fric| = "
                << (W_fric != 0.0 ? mon->SumAbsRho() / std::abs(W_fric) : 0.0) << "\n";
   }

   // rupture times
   {
      std::ofstream rt(out_dir + "/" + prefix + "_rupture_time_r" + std::to_string(g_rank) + ".csv");
      rt << std::setprecision(12) << "# x y z t_rupture (first |V| > 1e-3 m/s; -1 = never)\n";
      const Vector &X = op.FaultPositions();
      for (int j = 0; j < no; j++)
      {
         rt << X(3 * j) << " " << X(3 * j + 1) << " " << X(3 * j + 2) << " " << trupt[j] << "\n";
      }
   }
#ifdef MFEM_USE_HDF5
   // the slowest helper's deflate time (collective: pv exists on every rank)
   double pv_compress_max = 0.0;
   if (pv)
   {
      pv_compress_max = pv->CompressSeconds();
      MPI_Allreduce(MPI_IN_PLACE, &pv_compress_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
   }
   if (pv && g_rank == 0)
   {
      std::cout << std::setprecision(6) << "  fault ParaView output: " << pv->Frames()
                << " frames in " << out_dir << "/fault.vtkhdf; on rank 0 "
                << pv->Seconds() << " s ("
                << (pv->Frames() ? pv->Seconds() / pv->Frames() : 0.0) << " s/frame: pack "
                << (pv->Frames() ? pv->PackSeconds() / pv->Frames() : 0.0)
                << ", gather + HDF write "
                << (pv->Frames() ? (pv->Seconds() - pv->PackSeconds()) / pv->Frames() : 0.0);
      if (pv->Fast() && pv->Frames() && pv->SpreadActive())
      {
         std::cout << " = gather " << pv->GatherSeconds() / pv->Frames() << " + spread "
                   << pv->SpreadSeconds() / pv->Frames() << " + write "
                   << pv->WriteSeconds() / pv->Frames() << "; deflate on " << pv->SpreadHelpers()
                   << " ranks in chunks of " << pv->ChunkValues() << " values, slowest "
                   << pv_compress_max / pv->Frames() << " s/frame";
      }
      else if (pv->Fast() && pv->Frames())
      {
         std::cout << " = gather " << pv->GatherSeconds() / pv->Frames() << " + write "
                   << pv->WriteSeconds() / pv->Frames();
      }
      std::cout << "; " << (pv->Fast() ? "fast" : "SAFS") << " writer)\n";
   }
#endif
   {
      // per-rank time inside ode.Step after step 1: the spread shows load imbalance
      const long ns = nsteps - 1;
      double mn = step_sec, mx = step_sec, sm = step_sec, s1 = step1_sec;
      MPI_Allreduce(MPI_IN_PLACE, &mn, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &mx, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &sm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(MPI_IN_PLACE, &s1, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      if (g_rank == 0)
      {
         std::cout << std::setprecision(6) << "  time in ode.Step per rank after step 1: min " << mn
                   << " s, mean " << sm / Mpi::WorldSize() << " s, max " << mx << " s ("
                   << (ns > 0 ? 1e3 * mx / ns : 0.0) << " ms/step at the max); step 1 " << s1
                   << " s (with the operator assembly); fast fault "
                   << (op.FastFaultKernelsActive() ? "on" : "off") << ", fused K/K_eta "
                   << (op.FusedStiffnessActive() ? "on" : "off");
         if (op.FusedStiffnessActive())
         {
            std::cout << " (union pattern adds " << op.FusedExtraEntries() << " of "
                      << op.FusedEntriesK() << " entries)";
         }
         std::cout << ", dense blocks " << (op.DenseBlocksActive() ? "on" : "off");
         if (op.DenseBlocksActive())
         {
            std::cout << " (" << op.DenseBlockBytes() / 1e6 << " MB over all ranks)";
         }
         std::cout << "\n";
      }
      if (lts_run)
      {
         // element updates: patch elements kicked, against the global step's
         // two kicks of every element per dt0 tick
         long long w[3] = {lts_run->PatchUpdates(), pmesh.GetNE(), lts_run->Blocks().ExchangeRounds()};
         long long hv = lts_run->Blocks().ExchangeValues();
         // P4 (b): one merged exchange per tick plus Init's, on every rank
         int bad_rounds = (w[2] != lts_run->Ticks() + 1) ? 1 : 0;
         MPI_Allreduce(MPI_IN_PLACE, &bad_rounds, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
         MFEM_VERIFY(bad_rounds == 0, "LTS: " << bad_rounds << " rank(s) posted a halo count "
                     "other than ticks + 1 (" << lts_run->Ticks() + 1 << "); this rank "
                     << w[2]);
         MPI_Allreduce(MPI_IN_PLACE, w, 3, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &hv, 1, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
         // exposed exchange: seconds in WaitHalo, and its share of this
         // rank's stepping time (step 1 included, as WaitHalo's count is)
         const double wsec = lts_run->Blocks().WaitSeconds();
         double wr[3] = {wsec, wsec, wsec / std::max(step_sec + step1_sec, 1e-300)};
         MPI_Allreduce(MPI_IN_PLACE, &wr[0], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &wr[1], 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         MPI_Allreduce(MPI_IN_PLACE, &wr[2], 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
         if (g_rank == 0)
         {
            const long long ticks = lts_run->Ticks();
            std::cout << std::setprecision(6) << "  LTS: " << ticks << " ticks; patch element kicks "
                      << w[0] << " = " << double(w[0]) / (2.0 * double(w[1]) * double(ticks))
                      << " of global stepping's per tick; halo rounds " << w[2] / Mpi::WorldSize()
                      << " per rank (= ticks + 1 on every rank), " << hv
                      << " values received over all ranks; halo wait max " << wr[0] << " s, mean "
                      << wr[1] / Mpi::WorldSize() << " s, largest share of a rank's stepping "
                      << wr[2] << "\n";
         }
      }
      // cost model for the partition weight (plan item 6), after step 1: per
      // step, c_tet per tet of the non-fault work and c_side per fault-face
      // side of the fault work (an interior face has two sides on this rank,
      // a shared face one); w = c_side / c_tet.  On more than one rank the
      // fault time also holds waits for neighbours, so 1-rank runs give the
      // cleanest w.
      const double fsec = op.FaultSeconds() - fault_sec1;
      const auto &el = op.Elasticity();
      double loc[3] = {step_sec - fsec, fsec,
                       2.0 * el.FaultFaceOpsInterior() + el.FaultFaceOpsShared()};
      double glob[3];
      MPI_Allreduce(loc, glob, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      double fmx = fsec;
      MPI_Allreduce(MPI_IN_PLACE, &fmx, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      if (g_rank == 0 && ns > 0)
      {
         const double c_tet = glob[0] / ns / ne_global;
         const double c_side = glob[2] > 0 ? glob[1] / ns / glob[2] : 0.0;
         std::cout << std::setprecision(4) << "  fault work: " << 1e3 * fmx / ns
                   << " ms/step on the busiest rank; cost per step " << 1e6 * c_tet << " us per tet, "
                   << 1e6 * c_side << " us per fault-face side (w = " << c_side / c_tet << " tets)\n";
      }
      if (lts_run)
      {
         RankBalanceReport(*lts_run, lts_cl, pmesh.GetNE(), step_sec, step1_sec, loop_sec, fsec,
                           2.0 * el.FaultFaceOpsInterior() + el.FaultFaceOpsShared(),
                           out_dir + "/" + prefix + "_rank_balance.csv");
      }
   }
   const real_t wall = std::chrono::duration<real_t>(std::chrono::steady_clock::now() - wall0).count();
   if (g_rank == 0)
   {
      std::cout << std::setprecision(6) << "=== done: " << nsteps << " steps, t = " << t
                << " s; max|V| " << max_V_run << " m/s; max|V|dt/Dc " << max_ratio_run
                << "; sigma_n_eff in [" << min_sig_run << ", " << max_sig_run << "] Pa";
      if (pc_on) { std::cout << "; sigma_pc in [" << min_pc_run << ", " << max_pc_run << "] Pa"; }
      std::cout << "; wall " << wall << " s ===\n";
   }
   return 0;
}

// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// cycle_config.hpp — fully dynamic BP5 cycle (plan rev 3, Section 3.3 a):
// CycleConfig::Load(deck).
//
// One TOML deck holds both parts of the run (rule M6):
//   - the shared sections, read by the existing loader
//     spatial::LoadSpatialFrictionConfig unchanged: [meta], [mesh],
//     [material_constant_fallback], [pore_pressure], [velocity], [stress],
//     [numerics], [time], [output], [boundary], [friction.rate_state];
//   - the new [cycle] section, read here.
//
// What the cycle takes from the shared sections:
//   [mesh] path, order;  [material_constant_fallback] lambda, mu, rho;
//   [pore_pressure];  [friction.rate_state] (the resolver's block);
//   [boundary] fault_attr, natural_attrs;  [time] tfinal, dt_initial;
//   [output] output_dir;  [stress] sigma_n_pa.
// [stress] must be kind = "fault_local_prestress" with tau_strike_pa =
// tau_dip_pa = 0 and no patches: the shear prestress is the steady state that
// BuildPrestress computes, and a non-zero deck value would be ignored.  Its
// sigma_n_pa must equal [friction.rate_state].sigma_n_default.
//
// [cycle] keys (defaults in brackets):
//   ic                   "gold" | "scec_fd"                         [required]
//   switching            false: the QD regime only; true: the regime
//                        controller, which needs ic = "scec_fd".  The run
//                        starts dynamic at t = 0 if the QD balance's max V
//                        is above V_on (else in QD), switches back to QD at
//                        the end of each dynamic phase and to dynamic when
//                        the QD max V passes V_on again (P4), to t_final or
//                        [cycle.switch].max_dynamic_phases          [false]
//   qd_radiation_damping true; false is refused with "gold"         [true]
//   mesh_scale           vertex scale before the partition          [1.0]
//   plate_rate           Vp of the far-field and creep loading     [required]
//   penalty_factor       IP penalty scale of BOTH regimes' K (one value,
//                        so the two stiffness matrices agree, plan 6.1) [1.0]
//   [cycle.qd]      solver ("mumps" | "mumps_blr"), blr_tol,
//                   petsc_options (PETSc TS options file; it must pick an RK
//                   tableau with the FSAL property, see qd_phase.hpp),
//                   ts_atol, ts_rtol,
//                   operators ["loops"]: "matrix" runs each stage through
//                   the domain's prebuilt stage operators (seas_bp5_full
//                   --qd-operators matrix, QD speed item 1)
//   [cycle.qd.green] the discrete Green's function stage (seas_bp5_full
//                   --qd-green, items 3 and 4); present means on unless
//                   enabled = false, and it needs operators = "matrix".
//                   Linear elastic bulk and loading g(t) = t g only, so it
//                   is never on by default.  compress [0 = dense G; eps in
//                   (0, 1): the H-matrix], cache, hcache (files read when
//                   they fit, else written; the driver's --qd-green-cache
//                   and --qd-green-hcache override them), block [64],
//                   tol [1e-8], leaf [64], eta [2], compress_mem_gib [64],
//                   check_every [0 = off]: seas_bp5_full's --qd-green-* values
//   [cycle.switch]  V_on [1e-2], V_off [1e-3], T_dwell [10], T_min [0],
//                   T_max [1200], max_dynamic_phases [0 = no limit; N > 0:
//                   the run stops at the end of the N-th dynamic phase,
//                   before switching back]
//   [cycle.probe]   V_off, T_dwell: arrays of positive numbers (both or
//                   neither; default none; a V_off above V_on is allowed,
//                   since the pairs are only measured).
//                   The switch probe of plan P3: the first dynamic phase
//                   runs on until every (V_off, T_dwell) pair would have
//                   switched back (or T_max), and at each pair's switch the
//                   driver solves the QD balance and writes the jump in tau
//                   and V at every station to _switch_probe.txt.  The run
//                   itself switches back with [cycle.switch]'s pair.
//   [cycle.prestress] V_background (V_bg of psi0)                   [required]
//                   eta_V_init: tau0 += eta |V_init| (BuildPrestress
//                   include_eta: eq. 20, and eq. 23 with dtau = eta V_i, the
//                   QD convention).  "gold" always adds it (false is
//                   refused); with "scec_fd" it is false by default (eq. 22,
//                   eq. 23 with dtau = 0) and true makes the dynamic regime
//                   start the BP5 patch at V_i = 0.03 m/s from u = u' = 0 (the
//                   balance tau0 - eta V = friction at t = 0 then holds at
//                   V_init, as in SBEM and BICyclE) instead of the 0.015 m/s
//                   of tau_i - eta V.  One tau0 for both regimes, the whole
//                   run; checkpoints record the IC set as
//                   "scec_fd+eta_V_init" (IcLabel)          [ic == "gold"]
//   [cycle.fields]  V_init_dip, and the RectTaperFields a, Dc, V_init:
//                   center, axis1, axis2 (3-vectors, global XYZ), half1,
//                   half2, taper, eps, inside, outside               [required]
//   [cycle.output]  prefix, fast_V, dt_fast, dt_slow, event_V,
//                   event_dwell, plane_origin, plane_axis1, plane_axis2,
//                   and [[cycle.output.station]] name, x (global XYZ)
//   [cycle.checkpoint] every_qd_steps (0 = no checkpoint; else also one at
//                   the end of a run, as seas_bp5_full) [0]
//   [cycle.snapshots.fault], [cycle.snapshots.surface]  whole-field
//                   snapshots, one VTKHDF file per kind (cycle/snapshots.hpp).
//                   A table that is present is on unless enabled = false.
//                   dt_qd [3.15576e8 s = 10 yr]: a snapshot at the first QD
//                   step at or after each multiple of dt_qd since t = 0;
//                   dt_dyn [0.1 s]: at the first dynamic step at or after
//                   each multiple of dt_dyn since the phase started.  Plus
//                   one at the start, at every switch and at the end.
//                   surface only: attrs [[1]] (boundary attributes) and the
//                   plane n . x = plane_offset (plane_normal [[0, 0, 1]],
//                   plane_offset [0], plane_tol [1e-3] m): a boundary face
//                   is taken when its attribute is listed and all its
//                   vertices lie on the plane (attribute 1 also holds the
//                   bottom of the BP5 mesh).
//   [cycle.snapshots.volume]  only enabled = false (not implemented)
//   [cycle.benchmark]  the SCEC BP5 outputs beyond the on-fault stations
//                   (BP5 description section 4; cycle/bench_output.hpp),
//                   present means on unless enabled = false.  Plane
//                   coordinates (x2, x3) are those of [cycle.output]:
//                   (x - plane_origin) . plane_axis1 and . plane_axis2.
//                   region_x2, region_x3: [min, max] of the domain A of the
//                   source time series (max |V| and moment rate over A)
//                                                                [required]
//                   rupture_V: the rupture time is the first time |V| >=
//                   rupture_V in an event, and the rupture area of the
//                   catalog holds the DOFs that reached it        [event_V]
//                   profile_dt_slow, profile_dt_fast: the profile row
//                   cadence, dt_fast while max V > fast_V   [5 yr, 1 s]
//                   volume_axes: 3 unit vectors e1, e2, e3 (global XYZ);
//                   an off-fault row holds u . e_i and v . e_i
//                                          [required with volume stations]
//                   [[cycle.benchmark.volume_station]] name, x (global XYZ)
//                   [[cycle.benchmark.profile]] name, from, to (plane
//                   coordinates [x2, x3]), points (>= 2), coordinate ("x2"
//                   or "x3": the one the first row lists)
//   [cycle.lts]     local time stepping of the dynamic regime (the Nitsche
//                   driver's --lts rate2, solver/nitsche_lts_run.hpp):
//                   mode ["off" | "rate2"], cluster_cap [4], fault_buffer
//                   [3] (layers of tets around the fault in its cluster),
//                   fault_weight [0] (METIS weight per fault-face side).
//                   rate2 needs [numerics] cfl <= 0.6 and kv_scale = "node",
//                   and partitions the mesh with BuildLtsAwarePartition
//                   (both regimes then use that partition)
//   [cycle.loading] the Dirichlet loading of both regimes (2D in-plane
//                   design rev 3, phase A; cycle/loading.hpp).  kind
//                   ["bp5": MakeBP5DirichletFunc, the fault plane y = 0,
//                   slip along x, band 1 km; no other key] or "plate":
//                   plane_point p [required], plane_normal m [the
//                   [fault_geometry] ref_normal], slip_dir s [required],
//                   band b [m, required]; unit m and s with s . m = 0, and m
//                   not in the fault plane (|m . ref_normal| >= 1e-6).
//                     d = (x - p) . m;  u = +(Vp t / 2) s  for d > b,
//                                       u = -(Vp t / 2) s  for d < -b,
//                                       u = sgn(-m . ref_normal) Vp t s inside
//                   (inside the band: the jump across the prescribed-slip
//                   faces, which ComputeSkeletonDirichletSign signs by the
//                   oriented normal; the factor keeps that jump equal to the
//                   side m points to moving by +(Vp t / 2) s)
//
// The dynamic regime reads the keys the Nitsche driver reads, from the same
// shared sections (drivers/nitsche_dyn_driver.cpp:616-648, 708-736,
// 849-884), with its defaults:
//   [numerics]  cfl [0.5] (fraction of C(p); the shared loader reads it),
//               kv_eta [0 = Kelvin-Voigt off], kv_scale ["node" | "dt"],
//               kv_length ["inscribed" | "edge"], kv_layer_width [m; absent
//               = the whole domain].  No other key: penalty_factor lives in
//               [cycle], the velocity-stress keys mean nothing here.
//   [friction.prakash_clifton]  v_star, length (both required if present)
//   [fault_geometry]  ref_normal, up (the shared loader), orientation
//               ["propagated"], pin_by ["reference"], tangent ["pointwise"],
//               junction_policy ["refuse"] (the only implemented values),
//               pin_min_margin [the orientation's default]
//   [pore_pressure] min_sigma_n_pa: the lower bound of sigma_n_eff [0]
//
// No solver header (rule M4).

#ifndef MFEM_SEAS_CYCLE_CYCLE_CONFIG_HPP
#define MFEM_SEAS_CYCLE_CYCLE_CONFIG_HPP

#include "mfem.hpp"
#include "fault_fields.hpp"
#include "prestress.hpp"
#include "../spatial/code/spatial_friction.hpp"

#ifdef SEAS_USE_TOML
#include <toml.hpp>
#endif

#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{
namespace cycle
{

/// One on-fault station, global XYZ [m].
struct StationSpec
{
   std::string name;
   std::array<real_t, 3> x = {{0.0, 0.0, 0.0}};
};

/// [cycle.output]
struct OutputConfig
{
   std::string dir = ".";      ///< from the shared [output].output_dir
   std::string prefix;         ///< file prefix inside dir
   /// Row cadence of the station files: dt_fast while max V > fast_V, else
   /// dt_slow; a row is written once t - t_last >= 0.99 dt (io/
   /// bp5_parallel_output.hpp).  The SCEC spec says 0.1 s above 1e-3 m/s.
   real_t fast_V = 1.0e-3;
   real_t dt_fast = 0.1;
   real_t dt_slow = 0.1 * 365.25 * 24.0 * 3600.0;
   /// Event catalog: an event starts when max V > event_V and ends after
   /// event_dwell seconds below it (the SCEC definition).
   real_t event_V = 1.0e-3;
   real_t event_dwell = 10.0;
   /// The fault plane frame used to locate stations on the fault: 2-D
   /// coordinates (x - o).e1, (x - o).e2 of stations and fault DOFs.
   std::array<real_t, 3> plane_origin = {{0.0, 0.0, 0.0}};
   std::array<real_t, 3> plane_axis1 = {{1.0, 0.0, 0.0}};
   std::array<real_t, 3> plane_axis2 = {{0.0, 0.0, -1.0}};
   std::vector<StationSpec> stations;
};

/// One [cycle.snapshots.<kind>] table.
struct SnapshotKindConfig
{
   bool enabled = false;
   real_t dt_qd = 10.0 * 365.25 * 24.0 * 3600.0;   ///< [s], marks from t = 0
   real_t dt_dyn = 0.1;                            ///< [s], marks from the phase start
};

/// [cycle.snapshots]: the fault and the free surface (volume: not implemented).
struct SnapshotConfig
{
   SnapshotKindConfig fault;
   SnapshotKindConfig surface;
   std::vector<int> surface_attrs = {1};
   std::array<real_t, 3> surface_normal = {{0.0, 0.0, 1.0}};   ///< unit
   real_t surface_offset = 0.0;   ///< the plane n . x = offset [m]
   real_t surface_tol = 1.0e-3;   ///< [m]
   bool Any() const { return fault.enabled || surface.enabled; }
};

/// [cycle.qd.green]: the Green's-function stage (seas_bp5_full --qd-green-*).
struct QDGreenDeck
{
   bool enabled = false;
   real_t compress = 0.0;        ///< 0: dense G; (0, 1): H-matrix tolerance
   std::string cache;            ///< dense G file (read if it fits, else written)
   std::string hcache;           ///< H-matrix file (likewise)
   int block = 64;
   real_t tol = 1e-8;
   int leaf = 64;
   real_t eta = 2.0;
   real_t compress_mem_gib = 64.0;
   int check_every = 0;
};

/// [cycle.qd]
struct QDConfig
{
   std::string solver = "mumps";   ///< "mumps" or "mumps_blr"
   /// MUMPS-BLR tolerance.  1e-12 is the operator's default, and the value
   /// seas_bp5_full's BLR runs use: its --blr-tol default of 1e-10 is passed
   /// to SetBLRTol only when it differs from 1e-10.
   real_t blr_tol = 1e-12;
   std::string petsc_options = "tests/verification/petsc_ts_rk45_tandem.cfg";
   real_t ts_atol = 1e-7;
   real_t ts_rtol = 1e-50;
   std::string operators = "loops";   ///< "loops" or "matrix"
   QDGreenDeck green;
};

/// One [[cycle.benchmark.profile]]: points evenly spaced from `from` to `to`
/// (plane coordinates x2, x3) [m].
struct ProfileSpec
{
   std::string name;
   std::array<real_t, 2> from = {{0.0, 0.0}};
   std::array<real_t, 2> to = {{0.0, 0.0}};
   int points = 2;
   int coordinate = 0;   ///< 0: the first row lists x2; 1: x3
};

/// [cycle.benchmark]: the SCEC BP5 outputs of section 4 besides the on-fault
/// stations (see the file comment).
struct BenchmarkConfig
{
   bool enabled = false;
   std::array<real_t, 2> region_x2 = {{0.0, 0.0}};
   std::array<real_t, 2> region_x3 = {{0.0, 0.0}};
   real_t rupture_V = -1.0;   ///< < 0: [cycle.output].event_V
   real_t profile_dt_slow = 5.0 * 365.25 * 24.0 * 3600.0;
   real_t profile_dt_fast = 1.0;
   std::array<std::array<real_t, 3>, 3> volume_axes = {};
   std::vector<StationSpec> volume_stations;
   std::vector<ProfileSpec> profiles;
};

/// [cycle.lts]: local time stepping of the dynamic regime.
struct LtsConfig
{
   std::string mode = "off";   ///< "off" or "rate2"
   int cluster_cap = 4;
   int fault_buffer = 3;
   real_t fault_weight = 0.0;
   bool On() const { return mode == "rate2"; }
};

/// [cycle.switch]: the regime controller (plan 3.3 g).  Defaults: Thrase's
/// BP1-FD thresholds, the spec's 10 s event end, T_max about six times the
/// longest Tandem p=1 phase (plan 6.5).
struct SwitchConfig
{
   real_t V_on = 1.0e-2;    ///< dynamic when the global max V exceeds it [m/s]
   real_t V_off = 1.0e-3;   ///< back to QD once max V stays below it ...
   real_t T_dwell = 10.0;   ///< ... for this long [s]
   real_t T_min = 0.0;      ///< no switch back before this phase time [s]
   real_t T_max = 1200.0;   ///< a longer dynamic phase stops the run [s]
   /// Stop at the end of this many dynamic phases; 0 = no limit.
   int max_dynamic_phases = 0;
};

/// [cycle.probe]: the switch probe of plan P3 (empty: off).  Every pair
/// (V_off[i], T_dwell[j]) is measured.
struct ProbeConfig
{
   std::vector<real_t> V_off;
   std::vector<real_t> T_dwell;
   bool Enabled() const { return !V_off.empty(); }
};

/// [cycle.loading] (see the file comment and cycle/loading.hpp).
struct LoadingConfig
{
   std::string kind = "bp5";                                ///< "bp5" | "plate"
   std::array<real_t, 3> plane_point = {{0.0, 0.0, 0.0}};   ///< p [m]
   std::array<real_t, 3> plane_normal = {{0.0, 1.0, 0.0}};  ///< m, unit
   std::array<real_t, 3> slip_dir = {{1.0, 0.0, 0.0}};      ///< s, unit, s . m = 0
   real_t band = 1000.0;                                    ///< b [m]
   real_t band_sign = 1.0;   ///< sgn(-m . ref_normal), set by Load
   bool Plate() const { return kind == "plate"; }
};

/// The dynamic regime's numerics (see the file comment for the keys).
struct DynamicConfig
{
   real_t cfl = 0.5;
   real_t kv_eta = 0.0;
   std::string kv_scale = "node";
   std::string kv_length = "inscribed";
   real_t kv_layer_width = -1.0;   ///< < 0: the whole domain
   bool pc = false;                ///< Prakash-Clifton on
   real_t pc_v_star = 0.0;         ///< [m/s]
   real_t pc_length = 0.0;         ///< [m]
   std::array<real_t, 3> ref_normal = {{0.0, -1.0, 0.0}};
   std::array<real_t, 3> up = {{0.0, 0.0, 1.0}};
   real_t pin_min_margin = -1.0;   ///< < 0: the orientation's default
   real_t sigma_min = 0.0;         ///< [Pa]
};

struct CycleConfig
{
   std::string deck_path;
   spatial::SpatialFrictionConfig shared;   ///< the shared sections, unchanged

   std::string ic = "gold";
   bool switching = false;
   bool qd_radiation_damping = true;
   real_t mesh_scale = 1.0;
   real_t plate_rate = 0.0;   ///< Vp [m/s] of the Dirichlet loading
   real_t penalty_factor = 1.0;   ///< IP penalty scale of both regimes

   QDConfig qd;
   SwitchConfig switch_cfg;
   ProbeConfig probe;
   DynamicConfig dyn;
   PrestressConfig prestress;     ///< include_eta: [cycle.prestress].eta_V_init
   FaultFieldsConfig fields;
   OutputConfig output;
   SnapshotConfig snapshots;
   BenchmarkConfig bench;
   LtsConfig lts;
   LoadingConfig loading;         ///< [cycle.loading] (default: the BP5 function)
   int checkpoint_every = 0;

   // Shortcuts into `shared`, checked by Load.
   real_t t_final = 0.0;
   real_t dt_init = 0.0;
   real_t sigma_n = 0.0;
   int fault_attr = 3;

   bool IsGold() const { return ic == "gold"; }
   /// The IC set as checkpoints record it: "scec_fd+eta_V_init" when
   /// [cycle.prestress].eta_V_init adds eta |V_init| to the scec_fd tau0, so
   /// a restart refuses a deck with the other prestress.
   std::string IcLabel() const
   {
      return (prestress.include_eta && !IsGold()) ? ic + "+eta_V_init" : ic;
   }

   static CycleConfig Load(const std::string &deck);
};

#ifdef SEAS_USE_TOML
namespace detail
{

inline real_t Real(const toml::value &t, const std::string &sec,
                   const std::string &key, real_t def, bool required = false)
{
   if (!t.contains(key))
   {
      MFEM_VERIFY(!required, "[" << sec << "]." << key << " is required");
      return def;
   }
   const toml::value &v = t.at(key);
   if (v.is_floating()) { return static_cast<real_t>(v.as_floating()); }
   if (v.is_integer()) { return static_cast<real_t>(v.as_integer()); }
   MFEM_ABORT("[" << sec << "]." << key << " must be a number");
   return def;
}

inline std::string Str(const toml::value &t, const std::string &sec,
                       const std::string &key, const std::string &def,
                       bool required = false)
{
   if (!t.contains(key))
   {
      MFEM_VERIFY(!required, "[" << sec << "]." << key << " is required");
      return def;
   }
   MFEM_VERIFY(t.at(key).is_string(), "[" << sec << "]." << key << " must be a string");
   return t.at(key).as_string();
}

inline bool Bool(const toml::value &t, const std::string &sec,
                 const std::string &key, bool def)
{
   if (!t.contains(key)) { return def; }
   MFEM_VERIFY(t.at(key).is_boolean(), "[" << sec << "]." << key << " must be true or false");
   return t.at(key).as_boolean();
}

inline int Int(const toml::value &t, const std::string &sec,
               const std::string &key, int def)
{
   if (!t.contains(key)) { return def; }
   MFEM_VERIFY(t.at(key).is_integer(), "[" << sec << "]." << key << " must be an integer");
   return static_cast<int>(t.at(key).as_integer());
}

inline std::array<real_t, 3> Vec3(const toml::value &t, const std::string &sec,
                                  const std::string &key,
                                  const std::array<real_t, 3> &def, bool required = false)
{
   if (!t.contains(key))
   {
      MFEM_VERIFY(!required, "[" << sec << "]." << key << " is required");
      return def;
   }
   const toml::value &v = t.at(key);
   MFEM_VERIFY(v.is_array() && v.as_array().size() == 3,
               "[" << sec << "]." << key << " must be an array of 3 numbers");
   std::array<real_t, 3> out;
   for (int d = 0; d < 3; d++)
   {
      const toml::value &c = v.as_array()[d];
      MFEM_VERIFY(c.is_floating() || c.is_integer(),
                  "[" << sec << "]." << key << "[" << d << "] must be a number");
      out[d] = c.is_floating() ? static_cast<real_t>(c.as_floating())
                               : static_cast<real_t>(c.as_integer());
   }
   return out;
}

/// A non-empty array of positive numbers.
inline std::vector<real_t> PositiveReals(const toml::value &t, const std::string &sec,
                                         const std::string &key)
{
   MFEM_VERIFY(t.contains(key), "[" << sec << "]." << key << " is required");
   const toml::value &v = t.at(key);
   MFEM_VERIFY(v.is_array() && !v.as_array().empty(),
               "[" << sec << "]." << key << " must be a non-empty array of numbers");
   std::vector<real_t> out;
   for (const toml::value &c : v.as_array())
   {
      MFEM_VERIFY(c.is_floating() || c.is_integer(),
                  "[" << sec << "]." << key << " must hold numbers only");
      out.push_back(c.is_floating() ? static_cast<real_t>(c.as_floating())
                                    : static_cast<real_t>(c.as_integer()));
      MFEM_VERIFY(out.back() > 0.0, "[" << sec << "]." << key << " values must be > 0");
   }
   return out;
}

/// A non-empty array of integers.
inline std::vector<int> Ints(const toml::value &t, const std::string &sec,
                             const std::string &key, const std::vector<int> &def)
{
   if (!t.contains(key)) { return def; }
   const toml::value &v = t.at(key);
   MFEM_VERIFY(v.is_array() && !v.as_array().empty(),
               "[" << sec << "]." << key << " must be a non-empty array of integers");
   std::vector<int> out;
   for (const toml::value &c : v.as_array())
   {
      MFEM_VERIFY(c.is_integer(), "[" << sec << "]." << key << " must hold integers only");
      out.push_back(static_cast<int>(c.as_integer()));
   }
   return out;
}

/// Refuses keys a table does not know, so a misspelt key cannot fall back to
/// its default unnoticed.
inline void OnlyKeys(const toml::value &t, const std::string &sec,
                     const std::vector<std::string> &known)
{
   MFEM_VERIFY(t.is_table(), "[" << sec << "] must be a table");
   for (const auto &kv : t.as_table())
   {
      bool ok = false;
      for (const auto &k : known) { if (kv.first == k) { ok = true; break; } }
      MFEM_VERIFY(ok, "[" << sec << "] has an unknown key '" << kv.first << "'");
   }
}

/// One [cycle.snapshots.<kind>] table: present means on unless enabled = false.
inline SnapshotKindConfig SnapshotKind(const toml::value &t, const std::string &sec,
                                       const std::vector<std::string> &extra_keys)
{
   std::vector<std::string> known = {"enabled", "dt_qd", "dt_dyn"};
   known.insert(known.end(), extra_keys.begin(), extra_keys.end());
   OnlyKeys(t, sec, known);
   SnapshotKindConfig k;
   k.enabled = Bool(t, sec, "enabled", true);
   k.dt_qd = Real(t, sec, "dt_qd", k.dt_qd);
   k.dt_dyn = Real(t, sec, "dt_dyn", k.dt_dyn);
   MFEM_VERIFY(k.dt_qd > 0.0 && k.dt_dyn > 0.0,
               "[" << sec << "] dt_qd and dt_dyn must be > 0");
   return k;
}

inline RectTaperSpec Field(const toml::value &fields, const std::string &name)
{
   const std::string sec = "cycle.fields." + name;
   RectTaperSpec f;
   if (!fields.contains(name)) { return f; }
   const toml::value &t = fields.at(name);
   OnlyKeys(t, sec, {"center", "axis1", "axis2", "half1", "half2", "taper", "eps",
                     "inside", "outside"});
   f.center = Vec3(t, sec, "center", f.center, true);
   f.axis1 = Vec3(t, sec, "axis1", f.axis1, true);
   f.axis2 = Vec3(t, sec, "axis2", f.axis2, true);
   f.half1 = Real(t, sec, "half1", 0.0, true);
   f.half2 = Real(t, sec, "half2", 0.0, true);
   f.taper = Real(t, sec, "taper", 0.0);
   f.eps = Real(t, sec, "eps", 0.0);
   f.inside = Real(t, sec, "inside", 0.0, true);
   f.outside = Real(t, sec, "outside", 0.0, true);
   f.set = true;
   f.Validate(name);
   return f;
}

} // namespace detail
#endif // SEAS_USE_TOML

inline CycleConfig CycleConfig::Load(const std::string &deck)
{
#ifndef SEAS_USE_TOML
   MFEM_ABORT("CycleConfig::Load: built without toml11 (extern/toml11); "
              "rebuild with the submodule initialised");
   return CycleConfig();
#else
   using namespace detail;
   CycleConfig c;
   c.deck_path = deck;
   c.shared = spatial::LoadSpatialFrictionConfig(deck);
   const spatial::SpatialFrictionConfig &s = c.shared;

   // ---- shared sections the cycle relies on ----
   MFEM_VERIFY(s.law == spatial::FrictionLawKind::RateState && s.rate_state.has_value(),
               "fd_cycle: [meta].law must be \"rate_state\" with a [friction.rate_state] block");
   MFEM_VERIFY(!s.velocity.use_sidecar && s.material.kind == spatial::MaterialKind::Constant,
               "fd_cycle: P1 supports the constant material only "
               "([velocity].use_sidecar = false, [material] kind constant)");
   MFEM_VERIFY(s.stress.kind == spatial::StressSourceKind::FaultLocalPrestress,
               "fd_cycle: [stress].kind must be \"fault_local_prestress\" (the cycle "
               "takes sigma_n from it; the shear prestress is BuildPrestress's steady state)");
   MFEM_VERIFY(s.stress.tau_strike_pa == 0.0 && s.stress.tau_dip_pa == 0.0 &&
               s.stress.fault_local_patches.empty(),
               "fd_cycle: [stress].tau_strike_pa and tau_dip_pa must be 0 and no "
               "[[stress.patch]] given: the shear prestress is the steady state at "
               "V_init (BuildPrestress), so a deck value would be ignored");
   c.sigma_n = s.stress.sigma_n_pa;
   MFEM_VERIFY(c.sigma_n == s.rate_state->sigma_n_default,
               "fd_cycle: [stress].sigma_n_pa (" << c.sigma_n << ") must equal "
               "[friction.rate_state].sigma_n_default (" << s.rate_state->sigma_n_default << ")");
   MFEM_VERIFY(s.stress.pore_pressure.P_p_pa == 0.0 &&
               s.stress.pore_pressure.P_p_grad_pa_per_m == 0.0,
               "fd_cycle: [pore_pressure] must be zero: [friction.rate_state]."
               "sigma_n_default is the effective normal stress the QD kernel uses");
   for (const auto &r : s.rate_state->spatial)
   {
      MFEM_VERIFY(std::isnan(r.V_init),
                  "fd_cycle: a [[friction.rate_state.spatial]] rule sets V_init; the "
                  "initial slip rate comes from [cycle.fields.V_init] only");
   }
   c.fault_attr = s.boundary.fault_attr;
   MFEM_VERIFY(c.fault_attr == 3 && s.boundary.natural_attrs == std::vector<int>{1} &&
               s.boundary.absorbing_attrs.empty(),
               "fd_cycle: P1 uses the QD operator's far-field boundary layout "
               "(ElasticityDomainOperator legacy ctor, BCMode::FarField): fault 3, "
               "natural [1], Dirichlet 5, no absorbing.  Set [boundary] fault_attr = 3, "
               "natural_attrs = [1], absorbing_attrs = []");
   c.t_final = s.time.tfinal;
   c.dt_init = s.time.dt_initial;
   MFEM_VERIFY(c.t_final > 0.0, "fd_cycle: [time].tfinal must be > 0");
   MFEM_VERIFY(c.dt_init > 0.0, "fd_cycle: [time].dt_initial must be an explicit "
               "value > 0 (the first QD step; seas_bp5_full --tandem-time-stepping uses 0.01 s)");
   c.output.dir = s.output.output_dir;

   // ---- [cycle] ----
   const toml::value root = toml::parse(deck);
   MFEM_VERIFY(root.contains("cycle"), "fd_cycle: the deck has no [cycle] section");
   const toml::value &cy = root.at("cycle");
   OnlyKeys(cy, "cycle", {"ic", "switching", "qd_radiation_damping", "mesh_scale",
                          "plate_rate", "penalty_factor", "qd", "switch", "probe",
                          "prestress", "fields", "output", "checkpoint", "snapshots",
                          "benchmark", "lts", "loading"});
   c.ic = Str(cy, "cycle", "ic", "", true);
   MFEM_VERIFY(c.ic == "gold" || c.ic == "scec_fd",
               "[cycle].ic must be \"gold\" or \"scec_fd\"; got \"" << c.ic << "\"");
   c.switching = Bool(cy, "cycle", "switching", false);
   MFEM_VERIFY(!c.switching || !c.IsGold(),
               "[cycle].switching = true is refused with ic = \"gold\": that set starts "
               "at V = V_i = 0.01 m/s = V_on, so the start regime would depend on "
               "round-off (plan 4.2); use ic = \"scec_fd\"");
   c.qd_radiation_damping = Bool(cy, "cycle", "qd_radiation_damping", true);
   MFEM_VERIFY(c.qd_radiation_damping || !c.IsGold(),
               "[cycle].qd_radiation_damping = false is refused with ic = \"gold\": "
               "that set's tau0 contains eta V (BP5 eq. 20)");
   MFEM_VERIFY(c.qd_radiation_damping, "[cycle].qd_radiation_damping = false is not "
               "implemented: the QD kernel always carries eta V");
   c.mesh_scale = Real(cy, "cycle", "mesh_scale", 1.0);
   MFEM_VERIFY(c.mesh_scale > 0.0, "[cycle].mesh_scale must be > 0");
   c.plate_rate = Real(cy, "cycle", "plate_rate", 0.0, true);
   MFEM_VERIFY(c.plate_rate > 0.0, "[cycle].plate_rate must be > 0");
   c.penalty_factor = Real(cy, "cycle", "penalty_factor", c.penalty_factor);
   MFEM_VERIFY(c.penalty_factor > 0.0, "[cycle].penalty_factor must be > 0");

   if (cy.contains("qd"))
   {
      const toml::value &q = cy.at("qd");
      MFEM_VERIFY(!q.contains("penalty_factor"), "[cycle.qd].penalty_factor is now "
                  "[cycle].penalty_factor: one value for both regimes' stiffness");
      OnlyKeys(q, "cycle.qd", {"solver", "blr_tol", "petsc_options", "ts_atol", "ts_rtol",
                               "operators", "green"});
      c.qd.solver = Str(q, "cycle.qd", "solver", c.qd.solver);
      c.qd.blr_tol = Real(q, "cycle.qd", "blr_tol", c.qd.blr_tol);
      c.qd.petsc_options = Str(q, "cycle.qd", "petsc_options", c.qd.petsc_options);
      c.qd.ts_atol = Real(q, "cycle.qd", "ts_atol", c.qd.ts_atol);
      c.qd.ts_rtol = Real(q, "cycle.qd", "ts_rtol", c.qd.ts_rtol);
      c.qd.operators = Str(q, "cycle.qd", "operators", c.qd.operators);
      if (q.contains("green"))
      {
         const std::string sec = "cycle.qd.green";
         const toml::value &g = q.at("green");
         OnlyKeys(g, sec, {"enabled", "compress", "cache", "hcache", "block", "tol",
                           "leaf", "eta", "compress_mem_gib", "check_every"});
         QDGreenDeck &gr = c.qd.green;
         gr.enabled = Bool(g, sec, "enabled", true);
         gr.compress = Real(g, sec, "compress", gr.compress);
         gr.cache = Str(g, sec, "cache", gr.cache);
         gr.hcache = Str(g, sec, "hcache", gr.hcache);
         gr.block = Int(g, sec, "block", gr.block);
         gr.tol = Real(g, sec, "tol", gr.tol);
         gr.leaf = Int(g, sec, "leaf", gr.leaf);
         gr.eta = Real(g, sec, "eta", gr.eta);
         gr.compress_mem_gib = Real(g, sec, "compress_mem_gib", gr.compress_mem_gib);
         gr.check_every = Int(g, sec, "check_every", gr.check_every);
         // seas_bp5_full's checks of the same values
         MFEM_VERIFY(gr.compress == 0.0 || (gr.compress > 0.0 && gr.compress < 1.0),
                     "[" << sec << "].compress must be 0 (dense G) or in (0, 1)");
         MFEM_VERIFY(gr.hcache.empty() || gr.compress > 0.0,
                     "[" << sec << "].hcache needs compress > 0");
         MFEM_VERIFY(gr.block >= 1 && gr.leaf >= 1 && gr.tol > 0.0 && gr.eta > 0.0 &&
                     gr.compress_mem_gib > 0.0 && gr.check_every >= 0,
                     "[" << sec << "] block, leaf, tol, eta and compress_mem_gib must be "
                     "> 0, check_every >= 0");
      }
   }
   MFEM_VERIFY(c.qd.solver == "mumps" || c.qd.solver == "mumps_blr",
               "[cycle.qd].solver must be \"mumps\" or \"mumps_blr\"");
   MFEM_VERIFY(c.qd.operators == "loops" || c.qd.operators == "matrix",
               "[cycle.qd].operators must be \"loops\" or \"matrix\"");
   MFEM_VERIFY(!c.qd.green.enabled || c.qd.operators == "matrix",
               "[cycle.qd.green] needs [cycle.qd].operators = \"matrix\" (the Green's "
               "function is built from the stage operators)");

   if (cy.contains("switch"))
   {
      const std::string sec = "cycle.switch";
      const toml::value &w = cy.at("switch");
      OnlyKeys(w, sec, {"V_on", "V_off", "T_dwell", "T_min", "T_max", "max_dynamic_phases"});
      SwitchConfig &sw = c.switch_cfg;
      sw.V_on = Real(w, sec, "V_on", sw.V_on);
      sw.V_off = Real(w, sec, "V_off", sw.V_off);
      sw.T_dwell = Real(w, sec, "T_dwell", sw.T_dwell);
      sw.T_min = Real(w, sec, "T_min", sw.T_min);
      sw.T_max = Real(w, sec, "T_max", sw.T_max);
      sw.max_dynamic_phases = Int(w, sec, "max_dynamic_phases", sw.max_dynamic_phases);
   }
   {
      const SwitchConfig &sw = c.switch_cfg;
      MFEM_VERIFY(sw.V_off > 0.0 && sw.V_on > sw.V_off,
                  "[cycle.switch] needs 0 < V_off < V_on; got V_off = " << sw.V_off
                  << ", V_on = " << sw.V_on);
      MFEM_VERIFY(sw.T_dwell >= 0.0 && sw.T_min >= 0.0 && sw.T_max > sw.T_min,
                  "[cycle.switch] needs T_dwell >= 0, T_min >= 0, T_max > T_min");
      MFEM_VERIFY(sw.max_dynamic_phases >= 0,
                  "[cycle.switch].max_dynamic_phases must be >= 0 (0 = no limit)");
   }

   if (cy.contains("probe"))
   {
      const std::string sec = "cycle.probe";
      const toml::value &pr = cy.at("probe");
      OnlyKeys(pr, sec, {"V_off", "T_dwell"});
      c.probe.V_off = PositiveReals(pr, sec, "V_off");
      c.probe.T_dwell = PositiveReals(pr, sec, "T_dwell");
      MFEM_VERIFY(c.switching, "[cycle.probe] needs [cycle].switching = true");
   }

   MFEM_VERIFY(cy.contains("prestress"), "[cycle.prestress] is required");
   {
      const toml::value &p = cy.at("prestress");
      OnlyKeys(p, "cycle.prestress", {"V_background", "eta_V_init"});
      c.prestress.V_background = Real(p, "cycle.prestress", "V_background", 0.0, true);
      MFEM_VERIFY(c.prestress.V_background > 0.0, "[cycle.prestress].V_background must be > 0");
      c.prestress.include_eta = Bool(p, "cycle.prestress", "eta_V_init", c.IsGold());
      MFEM_VERIFY(c.prestress.include_eta || !c.IsGold(),
                  "[cycle.prestress].eta_V_init = false with ic = \"gold\": the gold IC set "
                  "(Tandem convention) adds eta |V_init| to tau0");
   }

   MFEM_VERIFY(cy.contains("fields"), "[cycle.fields] is required");
   {
      const toml::value &f = cy.at("fields");
      OnlyKeys(f, "cycle.fields", {"V_init_dip", "a", "Dc", "V_init"});
      c.fields.V_init_dip = Real(f, "cycle.fields", "V_init_dip", 0.0, true);
      MFEM_VERIFY(c.fields.V_init_dip >= 0.0, "[cycle.fields].V_init_dip must be >= 0");
      c.fields.a = Field(f, "a");
      c.fields.Dc = Field(f, "Dc");
      c.fields.V_init = Field(f, "V_init");
      MFEM_VERIFY(c.fields.a.set && c.fields.Dc.set && c.fields.V_init.set,
                  "[cycle.fields.a], [cycle.fields.Dc] and [cycle.fields.V_init] are required");
   }

   MFEM_VERIFY(cy.contains("output"), "[cycle.output] is required");
   {
      const std::string sec = "cycle.output";
      const toml::value &o = cy.at("output");
      OnlyKeys(o, sec, {"prefix", "fast_V", "dt_fast", "dt_slow", "event_V", "event_dwell",
                        "plane_origin", "plane_axis1", "plane_axis2", "station"});
      OutputConfig &oc = c.output;
      oc.prefix = Str(o, sec, "prefix", "", true);
      oc.fast_V = Real(o, sec, "fast_V", oc.fast_V);
      oc.dt_fast = Real(o, sec, "dt_fast", oc.dt_fast);
      oc.dt_slow = Real(o, sec, "dt_slow", oc.dt_slow);
      oc.event_V = Real(o, sec, "event_V", oc.event_V);
      oc.event_dwell = Real(o, sec, "event_dwell", oc.event_dwell);
      oc.plane_origin = Vec3(o, sec, "plane_origin", oc.plane_origin);
      oc.plane_axis1 = Vec3(o, sec, "plane_axis1", oc.plane_axis1);
      oc.plane_axis2 = Vec3(o, sec, "plane_axis2", oc.plane_axis2);
      MFEM_VERIFY(oc.dt_fast > 0.0 && oc.dt_slow > 0.0 && oc.fast_V > 0.0 &&
                  oc.event_V > 0.0 && oc.event_dwell >= 0.0,
                  "[cycle.output] fast_V, dt_fast, dt_slow, event_V must be > 0, "
                  "event_dwell >= 0");
      auto dot = [](const std::array<real_t, 3> &u, const std::array<real_t, 3> &v)
      { return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]; };
      MFEM_VERIFY(std::abs(dot(oc.plane_axis1, oc.plane_axis1) - 1.0) < 1e-12 &&
                  std::abs(dot(oc.plane_axis2, oc.plane_axis2) - 1.0) < 1e-12 &&
                  std::abs(dot(oc.plane_axis1, oc.plane_axis2)) < 1e-12,
                  "[cycle.output] plane_axis1 and plane_axis2 must be orthonormal");
      if (o.contains("station"))
      {
         MFEM_VERIFY(o.at("station").is_array(),
                     "[[cycle.output.station]] must be an array of tables");
         for (const toml::value &st : o.at("station").as_array())
         {
            OnlyKeys(st, "cycle.output.station", {"name", "x"});
            StationSpec sp;
            sp.name = Str(st, "cycle.output.station", "name", "", true);
            sp.x = Vec3(st, "cycle.output.station", "x", sp.x, true);
            oc.stations.push_back(sp);
         }
      }
   }

   if (cy.contains("checkpoint"))
   {
      const toml::value &k = cy.at("checkpoint");
      OnlyKeys(k, "cycle.checkpoint", {"every_qd_steps"});
      c.checkpoint_every = Int(k, "cycle.checkpoint", "every_qd_steps", 0);
      MFEM_VERIFY(c.checkpoint_every >= 0, "[cycle.checkpoint].every_qd_steps must be >= 0");
   }

   if (cy.contains("snapshots"))
   {
      const toml::value &sn = cy.at("snapshots");
      OnlyKeys(sn, "cycle.snapshots", {"fault", "surface", "volume"});
      SnapshotConfig &sc = c.snapshots;
      if (sn.contains("fault"))
      {
         sc.fault = SnapshotKind(sn.at("fault"), "cycle.snapshots.fault", {});
      }
      if (sn.contains("surface"))
      {
         const std::string sec = "cycle.snapshots.surface";
         const toml::value &s2 = sn.at("surface");
         sc.surface = SnapshotKind(s2, sec, {"attrs", "plane_normal", "plane_offset",
                                             "plane_tol"});
         sc.surface_attrs = Ints(s2, sec, "attrs", sc.surface_attrs);
         sc.surface_normal = Vec3(s2, sec, "plane_normal", sc.surface_normal);
         sc.surface_offset = Real(s2, sec, "plane_offset", sc.surface_offset);
         sc.surface_tol = Real(s2, sec, "plane_tol", sc.surface_tol);
         const std::array<real_t, 3> &n = sc.surface_normal;
         MFEM_VERIFY(std::abs(n[0] * n[0] + n[1] * n[1] + n[2] * n[2] - 1.0) < 1e-12,
                     "[" << sec << "].plane_normal must be a unit vector");
         MFEM_VERIFY(sc.surface_tol > 0.0, "[" << sec << "].plane_tol must be > 0");
      }
      if (sn.contains("volume"))
      {
         const toml::value &v = sn.at("volume");
         OnlyKeys(v, "cycle.snapshots.volume", {"enabled"});
         MFEM_VERIFY(v.contains("enabled") && !Bool(v, "cycle.snapshots.volume", "enabled", true),
                     "[cycle.snapshots.volume] is not implemented: only enabled = false "
                     "is accepted");
      }
   }

   if (cy.contains("benchmark"))
   {
      const std::string sec = "cycle.benchmark";
      const toml::value &b = cy.at("benchmark");
      OnlyKeys(b, sec, {"enabled", "region_x2", "region_x3", "rupture_V", "profile_dt_slow",
                        "profile_dt_fast", "volume_axes", "volume_station", "profile"});
      BenchmarkConfig &bc = c.bench;
      bc.enabled = Bool(b, sec, "enabled", true);
      auto range = [&](const std::string &key)
      {
         MFEM_VERIFY(b.contains(key), "[" << sec << "]." << key << " is required");
         const toml::value &v = b.at(key);
         MFEM_VERIFY(v.is_array() && v.as_array().size() == 2,
                     "[" << sec << "]." << key << " must be [min, max]");
         std::array<real_t, 2> r;
         for (int i = 0; i < 2; i++)
         {
            const toml::value &e = v.as_array()[i];
            MFEM_VERIFY(e.is_floating() || e.is_integer(),
                        "[" << sec << "]." << key << " must hold numbers");
            r[i] = e.is_floating() ? static_cast<real_t>(e.as_floating())
                                   : static_cast<real_t>(e.as_integer());
         }
         MFEM_VERIFY(r[0] < r[1], "[" << sec << "]." << key << " needs min < max");
         return r;
      };
      bc.region_x2 = range("region_x2");
      bc.region_x3 = range("region_x3");
      bc.rupture_V = Real(b, sec, "rupture_V", c.output.event_V);
      bc.profile_dt_slow = Real(b, sec, "profile_dt_slow", bc.profile_dt_slow);
      bc.profile_dt_fast = Real(b, sec, "profile_dt_fast", bc.profile_dt_fast);
      MFEM_VERIFY(bc.rupture_V > 0.0 && bc.profile_dt_slow > 0.0 && bc.profile_dt_fast > 0.0,
                  "[" << sec << "] rupture_V, profile_dt_slow, profile_dt_fast must be > 0");
      if (b.contains("volume_station"))
      {
         MFEM_VERIFY(b.at("volume_station").is_array(),
                     "[[" << sec << ".volume_station]] must be an array of tables");
         for (const toml::value &st : b.at("volume_station").as_array())
         {
            const std::string ssec = sec + ".volume_station";
            OnlyKeys(st, ssec, {"name", "x"});
            StationSpec sp;
            sp.name = Str(st, ssec, "name", "", true);
            sp.x = Vec3(st, ssec, "x", sp.x, true);
            bc.volume_stations.push_back(sp);
         }
      }
      if (!bc.volume_stations.empty())
      {
         MFEM_VERIFY(b.contains("volume_axes") && b.at("volume_axes").is_array() &&
                     b.at("volume_axes").as_array().size() == 3,
                     "[" << sec << "].volume_axes (3 unit vectors) is required with "
                     "volume stations");
         const auto &ax = b.at("volume_axes").as_array();
         for (int i = 0; i < 3; i++)
         {
            MFEM_VERIFY(ax[i].is_array() && ax[i].as_array().size() == 3,
                        "[" << sec << "].volume_axes[" << i << "] must be 3 numbers");
            for (int d = 0; d < 3; d++)
            {
               const toml::value &e = ax[i].as_array()[d];
               MFEM_VERIFY(e.is_floating() || e.is_integer(),
                           "[" << sec << "].volume_axes must hold numbers");
               bc.volume_axes[i][d] = e.is_floating() ? static_cast<real_t>(e.as_floating())
                                                      : static_cast<real_t>(e.as_integer());
            }
         }
         for (int i = 0; i < 3; i++)
         {
            for (int j = 0; j < 3; j++)
            {
               const real_t d = bc.volume_axes[i][0] * bc.volume_axes[j][0] +
                                bc.volume_axes[i][1] * bc.volume_axes[j][1] +
                                bc.volume_axes[i][2] * bc.volume_axes[j][2];
               MFEM_VERIFY(std::abs(d - (i == j ? 1.0 : 0.0)) < 1e-12,
                           "[" << sec << "].volume_axes must be orthonormal");
            }
         }
      }
      else
      {
         MFEM_VERIFY(!b.contains("volume_axes"), "[" << sec << "].volume_axes without "
                     "any [[" << sec << ".volume_station]]");
      }
      if (b.contains("profile"))
      {
         MFEM_VERIFY(b.at("profile").is_array(),
                     "[[" << sec << ".profile]] must be an array of tables");
         for (const toml::value &pt : b.at("profile").as_array())
         {
            const std::string psec = sec + ".profile";
            OnlyKeys(pt, psec, {"name", "from", "to", "points", "coordinate"});
            ProfileSpec ps;
            ps.name = Str(pt, psec, "name", "", true);
            for (int k = 0; k < 2; k++)
            {
               const std::string key = k == 0 ? "from" : "to";
               MFEM_VERIFY(pt.contains(key) && pt.at(key).is_array() &&
                           pt.at(key).as_array().size() == 2,
                           "[" << psec << "]." << key << " must be [x2, x3]");
               std::array<real_t, 2> &dst = k == 0 ? ps.from : ps.to;
               for (int i = 0; i < 2; i++)
               {
                  const toml::value &e = pt.at(key).as_array()[i];
                  MFEM_VERIFY(e.is_floating() || e.is_integer(),
                              "[" << psec << "]." << key << " must hold numbers");
                  dst[i] = e.is_floating() ? static_cast<real_t>(e.as_floating())
                                           : static_cast<real_t>(e.as_integer());
               }
            }
            ps.points = Int(pt, psec, "points", 0);
            MFEM_VERIFY(ps.points >= 2, "[" << psec << "].points must be >= 2");
            const std::string coord = Str(pt, psec, "coordinate", "", true);
            MFEM_VERIFY(coord == "x2" || coord == "x3",
                        "[" << psec << "].coordinate must be \"x2\" or \"x3\"");
            ps.coordinate = coord == "x2" ? 0 : 1;
            for (const ProfileSpec &o : bc.profiles)
            {
               MFEM_VERIFY(o.name != ps.name, "[" << psec << "]: name '" << ps.name
                           << "' given twice");
            }
            bc.profiles.push_back(ps);
         }
      }
   }

   if (cy.contains("lts"))
   {
      const std::string sec = "cycle.lts";
      const toml::value &l = cy.at("lts");
      OnlyKeys(l, sec, {"mode", "cluster_cap", "fault_buffer", "fault_weight"});
      LtsConfig &lc = c.lts;
      lc.mode = Str(l, sec, "mode", lc.mode);
      lc.cluster_cap = Int(l, sec, "cluster_cap", lc.cluster_cap);
      lc.fault_buffer = Int(l, sec, "fault_buffer", lc.fault_buffer);
      lc.fault_weight = Real(l, sec, "fault_weight", lc.fault_weight);
      MFEM_VERIFY(lc.mode == "off" || lc.mode == "rate2",
                  "[" << sec << "].mode must be \"off\" or \"rate2\"");
      MFEM_VERIFY(lc.cluster_cap >= 1 && lc.fault_buffer >= 0 && lc.fault_weight >= 0.0,
                  "[" << sec << "] needs cluster_cap >= 1, fault_buffer >= 0, "
                  "fault_weight >= 0");
   }

   // ---- the dynamic regime's numerics (shared sections, Nitsche keys) ----
   {
      DynamicConfig &d = c.dyn;
      d.cfl = s.numerics.cfl;
      MFEM_VERIFY(d.cfl > 0.0 && d.cfl <= 1.0,
                  "[numerics].cfl (fraction of C(p)) must be in (0, 1]; got " << d.cfl);
      MFEM_VERIFY(s.mesh.order >= 1 && s.mesh.order <= 4,
                  "[mesh].order must be in 1..4 (the measured C(p) of the explicit step)");
      if (root.contains("numerics"))
      {
         const std::string sec = "numerics";
         const toml::value &n = root.at("numerics");
         MFEM_VERIFY(!n.contains("penalty_factor"), "[numerics].penalty_factor: the cycle "
                     "takes the penalty from [cycle].penalty_factor, one value for both regimes");
         OnlyKeys(n, sec, {"cfl", "kv_eta", "kv_scale", "kv_length", "kv_layer_width"});
         d.kv_eta = Real(n, sec, "kv_eta", d.kv_eta);
         d.kv_scale = Str(n, sec, "kv_scale", d.kv_scale);
         d.kv_length = Str(n, sec, "kv_length", d.kv_length);
         if (n.contains("kv_layer_width"))
         {
            d.kv_layer_width = Real(n, sec, "kv_layer_width", 0.0);
            MFEM_VERIFY(d.kv_layer_width > 0.0, "[numerics].kv_layer_width must be > 0 m");
         }
      }
      MFEM_VERIFY(d.kv_eta >= 0.0, "[numerics].kv_eta must be >= 0");
      MFEM_VERIFY(d.kv_scale == "node" || d.kv_scale == "dt",
                  "[numerics].kv_scale = \"" << d.kv_scale << "\"; use \"node\" or \"dt\"");
      MFEM_VERIFY(d.kv_length == "inscribed" || d.kv_length == "edge",
                  "[numerics].kv_length = \"" << d.kv_length
                  << "\"; use \"inscribed\" or \"edge\"");

      if (root.contains("friction") && root.at("friction").contains("prakash_clifton"))
      {
         const std::string sec = "friction.prakash_clifton";
         const toml::value &p = root.at("friction").at("prakash_clifton");
         OnlyKeys(p, sec, {"v_star", "length"});
         d.pc = true;
         d.pc_v_star = Real(p, sec, "v_star", 0.0, true);
         d.pc_length = Real(p, sec, "length", 0.0, true);
         MFEM_VERIFY(d.pc_v_star > 0.0 && d.pc_length > 0.0, "[" << sec << "] v_star = "
                     << d.pc_v_star << ", length = " << d.pc_length << "; both must be > 0");
      }

      d.ref_normal = s.fault_geometry.ref_normal;   // normalized by the shared loader
      d.up = s.fault_geometry.up;
      if (root.contains("fault_geometry"))
      {
         const std::string sec = "fault_geometry";
         const toml::value &g = root.at("fault_geometry");
         OnlyKeys(g, sec, {"ref_normal", "up", "kind", "orientation", "pin_by", "tangent",
                           "junction_policy", "pin_min_margin"});
         const std::string orient = Str(g, sec, "orientation", "propagated");
         const std::string pin_by = Str(g, sec, "pin_by", "reference");
         const std::string tangent = Str(g, sec, "tangent", "pointwise");
         const std::string junction = Str(g, sec, "junction_policy", "refuse");
         MFEM_VERIFY(orient == "propagated" && pin_by == "reference" &&
                     tangent == "pointwise" && junction == "refuse",
                     "[fault_geometry] orientation, pin_by, tangent, junction_policy: only "
                     "\"propagated\", \"reference\", \"pointwise\", \"refuse\" are implemented");
         if (g.contains("pin_min_margin"))
         {
            d.pin_min_margin = Real(g, sec, "pin_min_margin", 0.0);
            MFEM_VERIFY(d.pin_min_margin > 0.0, "[fault_geometry].pin_min_margin must be > 0");
         }
      }
      d.sigma_min = s.stress.pore_pressure.min_sigma_n_pa;
   }
   // ---- [cycle.loading] (after [fault_geometry]: m defaults to ref_normal) ----
   if (cy.contains("loading"))
   {
      const std::string sec = "cycle.loading";
      const toml::value &l = cy.at("loading");
      OnlyKeys(l, sec, {"kind", "plane_point", "plane_normal", "slip_dir", "band"});
      LoadingConfig &L = c.loading;
      L.kind = Str(l, sec, "kind", L.kind);
      MFEM_VERIFY(L.kind == "bp5" || L.kind == "plate",
                  "[" << sec << "].kind must be \"bp5\" or \"plate\"; got \"" << L.kind << "\"");
      if (L.kind == "bp5")
      {
         MFEM_VERIFY(!l.contains("plane_point") && !l.contains("plane_normal") &&
                     !l.contains("slip_dir") && !l.contains("band"),
                     "[" << sec << "] kind = \"bp5\" takes no other key (it is the fault "
                     "plane y = 0, slip along x, band 1 km); use kind = \"plate\"");
      }
      else
      {
         auto dot = [](const std::array<real_t, 3> &u, const std::array<real_t, 3> &v)
         { return u[0] * v[0] + u[1] * v[1] + u[2] * v[2]; };
         L.plane_point = Vec3(l, sec, "plane_point", L.plane_point, true);
         L.plane_normal = Vec3(l, sec, "plane_normal", c.dyn.ref_normal);
         L.slip_dir = Vec3(l, sec, "slip_dir", L.slip_dir, true);
         L.band = Real(l, sec, "band", 0.0, true);
         MFEM_VERIFY(std::abs(dot(L.plane_normal, L.plane_normal) - 1.0) < 1e-12 &&
                     std::abs(dot(L.slip_dir, L.slip_dir) - 1.0) < 1e-12,
                     "[" << sec << "] plane_normal and slip_dir must be unit vectors");
         MFEM_VERIFY(std::abs(dot(L.plane_normal, L.slip_dir)) < 1e-12,
                     "[" << sec << "] slip_dir must lie in the plane (slip_dir . "
                     "plane_normal = " << dot(L.plane_normal, L.slip_dir) << ")");
         MFEM_VERIFY(L.band >= 0.0, "[" << sec << "].band must be >= 0");
         const real_t mr = dot(L.plane_normal, c.dyn.ref_normal);
         MFEM_VERIFY(std::abs(mr) >= 1e-6, "[" << sec << "] plane_normal . "
                     "[fault_geometry].ref_normal = " << mr << ": the plane normal lies in "
                     "the fault plane, so the side of the prescribed jump is undefined");
         L.band_sign = (mr < 0.0) ? 1.0 : -1.0;
      }
   }
   if (!c.loading.Plate())
   {
      // The BP5 function's band value +Vp t x is the jump of the side y > 0
      // over the side y < 0, which ComputeSkeletonDirichletSign gives only
      // with the oriented normal along (0, -1, 0).
      const std::array<real_t, 3> &r = c.dyn.ref_normal;
      MFEM_VERIFY(std::abs(r[0]) < 1e-12 && std::abs(r[1] + 1.0) < 1e-12 &&
                  std::abs(r[2]) < 1e-12,
                  "[cycle.loading] kind = \"bp5\" (the default) is the BP5 function of "
                  "the plane y = 0 and needs [fault_geometry].ref_normal = (0, -1, 0); got ("
                  << r[0] << ", " << r[1] << ", " << r[2] << "): use kind = \"plate\"");
   }

   if (c.lts.On())
   {
      // The Nitsche driver's --lts rate2 conditions (nitsche_lts_run.hpp)
      MFEM_VERIFY(c.dyn.cfl <= 0.6, "[cycle.lts] mode = \"rate2\" needs [numerics].cfl "
                  "<= 0.6 (the rate-2 clusters' stability margin); got " << c.dyn.cfl);
      MFEM_VERIFY(c.dyn.kv_scale == "node", "[cycle.lts] mode = \"rate2\" needs "
                  "[numerics].kv_scale = \"node\" (a dt-scaled KV would differ per cluster)");
   }
   return c;
#endif
}

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_CYCLE_CONFIG_HPP

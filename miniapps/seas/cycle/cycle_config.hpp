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
//                        controller, which needs ic = "scec_fd".  In this
//                        build (P3) the run starts dynamic at t = 0, switches
//                        back to QD at the end of the phase and continues to
//                        t_final; a QD -> dynamic switch at t > 0 (P4) stops
//                        the run instead                            [false]
//   qd_radiation_damping true; false is refused with "gold"         [true]
//   mesh_scale           vertex scale before the partition          [1.0]
//   plate_rate           Vp of the far-field and creep loading     [required]
//   penalty_factor       IP penalty scale of BOTH regimes' K (one value,
//                        so the two stiffness matrices agree, plan 6.1) [1.0]
//   [cycle.qd]      solver ("mumps" | "mumps_blr"), blr_tol,
//                   petsc_options (PETSc TS options file), ts_atol, ts_rtol
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
//   [cycle.fields]  V_init_dip, and the RectTaperFields a, Dc, V_init:
//                   center, axis1, axis2 (3-vectors, global XYZ), half1,
//                   half2, taper, eps, inside, outside               [required]
//   [cycle.output]  prefix, fast_V, dt_fast, dt_slow, event_V,
//                   event_dwell, plane_origin, plane_axis1, plane_axis2,
//                   and [[cycle.output.station]] name, x (global XYZ)
//   [cycle.checkpoint] every_qd_steps (0 = no checkpoint; else also one at
//                   the end of a run, as seas_bp5_full) [0]
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
   PrestressConfig prestress;     ///< include_eta follows ic
   FaultFieldsConfig fields;
   OutputConfig output;
   int checkpoint_every = 0;

   // Shortcuts into `shared`, checked by Load.
   real_t t_final = 0.0;
   real_t dt_init = 0.0;
   real_t sigma_n = 0.0;
   int fault_attr = 3;

   bool IsGold() const { return ic == "gold"; }

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
                          "prestress", "fields", "output", "checkpoint"});
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
      OnlyKeys(q, "cycle.qd", {"solver", "blr_tol", "petsc_options", "ts_atol", "ts_rtol"});
      c.qd.solver = Str(q, "cycle.qd", "solver", c.qd.solver);
      c.qd.blr_tol = Real(q, "cycle.qd", "blr_tol", c.qd.blr_tol);
      c.qd.petsc_options = Str(q, "cycle.qd", "petsc_options", c.qd.petsc_options);
      c.qd.ts_atol = Real(q, "cycle.qd", "ts_atol", c.qd.ts_atol);
      c.qd.ts_rtol = Real(q, "cycle.qd", "ts_rtol", c.qd.ts_rtol);
   }
   MFEM_VERIFY(c.qd.solver == "mumps" || c.qd.solver == "mumps_blr",
               "[cycle.qd].solver must be \"mumps\" or \"mumps_blr\"");

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
      OnlyKeys(p, "cycle.prestress", {"V_background"});
      c.prestress.V_background = Real(p, "cycle.prestress", "V_background", 0.0, true);
      MFEM_VERIFY(c.prestress.V_background > 0.0, "[cycle.prestress].V_background must be > 0");
      c.prestress.include_eta = c.IsGold();
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
   return c;
#endif
}

} // namespace cycle
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CYCLE_CYCLE_CONFIG_HPP

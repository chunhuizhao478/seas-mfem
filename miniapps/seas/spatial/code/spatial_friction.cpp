// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// spatial_friction.cpp — Phase 1 implementation of
// spatial_dynamic_rupture_plan.md rev-3.
//
// All TOML I/O is gated by SEAS_USE_TOML.  Without toml11, the parser
// entry points abort at runtime; the resolver and time-helper paths
// compile unconditionally.

#include "spatial_friction.hpp"

// Phase 1 (QD): the full definition of `mfem::seas::SolverType` (forward-
// declared in spatial_friction.hpp) is needed here to define
// ParseQDSolverType.  Included ONLY in this .cpp — never the header — so the
// heavy elasticity-operator include does not propagate to the many TUs that
// include spatial_friction.hpp (incl. spatial_dyn_driver).  Read-only use of
// the SolverType enum; no ODR-use of any elasticity-operator symbol, so the
// minimal link set (spatial_friction.o + MFEM_LIBS) is preserved.
#include "../../domain/elasticity_operator.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <regex>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#ifdef SEAS_USE_TOML
#include <toml.hpp>
#endif

namespace mfem
{
namespace seas
{
namespace spatial
{

// =====================================================================
//  SpatialTimeParseSeconds
// =====================================================================

real_t SpatialTimeParseSeconds(const std::string& s_in)
{
   // Strip surrounding whitespace.
   auto lstrip = s_in.find_first_not_of(" \t\n\r");
   auto rstrip = s_in.find_last_not_of(" \t\n\r");
   std::string s = (lstrip == std::string::npos)
                   ? std::string()
                   : s_in.substr(lstrip, rstrip - lstrip + 1);

   // Lower-case copy for sentinel checks.
   std::string lower(s);
   std::transform(lower.begin(), lower.end(), lower.begin(),
                  [](unsigned char c) { return std::tolower(c); });

   if (lower == "auto") { return -1.0; }

   MFEM_VERIFY(!s.empty(),
               "SpatialTimeParseSeconds: empty time string");

   // Optional trailing "s" suffix; rest must parse as a floating-point
   // number with std::strtod.
   std::string numpart = s;
   if (!numpart.empty() && (numpart.back() == 's' || numpart.back() == 'S'))
   {
      numpart.pop_back();
   }
   // Trim trailing whitespace between the number and the 's'.
   while (!numpart.empty()
          && (numpart.back() == ' ' || numpart.back() == '\t'))
   {
      numpart.pop_back();
   }
   MFEM_VERIFY(!numpart.empty(),
               "SpatialTimeParseSeconds: empty numeric body in '" << s_in << "'");

   const char* cstr = numpart.c_str();
   char* end = nullptr;
   const double val = std::strtod(cstr, &end);
   MFEM_VERIFY(end != cstr && *end == '\0',
               "SpatialTimeParseSeconds: cannot parse '" << s_in
               << "' as a real with optional 's' suffix");
   return static_cast<real_t>(val);
}

// =====================================================================
//  SpatialRule::matches
// =====================================================================

bool SpatialRule::matches(real_t x, real_t y, real_t z, int attr) const
{
   const bool x_ok = (x >= x_min_m) && (x <= x_max_m);
   const bool y_ok = (y >= y_min_m) && (y <= y_max_m);
   const bool z_ok = (z >= z_min_m) && (z <= z_max_m);

   switch (kind)
   {
   case Kind::Depth:
      return z_ok;
   case Kind::Box:
      return x_ok && y_ok && z_ok;
   case Kind::RegionAttribute:
      return (region_attr >= 0) && (attr == region_attr);
   case Kind::Barrier:
      // Same as Depth, with any non-infinite x/y bounds also applied.
      return z_ok && x_ok && y_ok;
   case Kind::BoxcarTaper:
      // Phase 6 req 5: the taper geometry (not the coord bounds) governs the
      // region; the rule applies wherever the boxcar factor is non-zero.
      return BoxcarTaperFactor(x, y, z) > 0.0;
   }
   return false;
}

real_t SpatialRule::BoxcarTaperFactor(real_t x, real_t y, real_t z) const
{
   if (kind != Kind::BoxcarTaper) { return 1.0; }
   // Per-axis SCEC boxcar; a non-finite half makes that axis untapered.
   auto axis = [](real_t coord, real_t center, real_t half, real_t trans)
   {
      if (!std::isfinite(half)) { return static_cast<real_t>(1.0); }
      return SCECBoxcar(coord - center, half, trans);
   };
   return axis(x, boxcar_center_x_m, boxcar_half_x_m, boxcar_trans_x_m)
        * axis(y, boxcar_center_y_m, boxcar_half_y_m, boxcar_trans_y_m)
        * axis(z, boxcar_center_z_m, boxcar_half_z_m, boxcar_trans_z_m);
}

// =====================================================================
//  Depth-profile interpolant + two-CSV loader (Phase 11b)
//  (No toml11 dependency — the resolver uses these unconditionally.)
// =====================================================================

real_t PiecewiseLinear1D::operator()(real_t xq) const
{
   const std::size_t n = x.size();
   MFEM_ASSERT(n >= 2 && y.size() == n,
               "PiecewiseLinear1D: call Validate() before evaluating");
   // Flat (constant) clamp outside the sampled range.
   if (xq <= x.front()) { return y.front(); }
   if (xq >= x.back())  { return y.back(); }
   // Binary search for the bracketing interval [x[k-1], x[k]] (1 <= k <= n-1).
   const auto it = std::upper_bound(x.begin(), x.end(), xq);
   const std::size_t k = static_cast<std::size_t>(it - x.begin());
   const real_t x0 = x[k - 1], x1 = x[k];
   const real_t y0 = y[k - 1], y1 = y[k];
   const real_t t = (xq - x0) / (x1 - x0);
   return y0 + t * (y1 - y0);
}

void PiecewiseLinear1D::Validate() const
{
   MFEM_VERIFY(x.size() >= 2,
               "PiecewiseLinear1D: need >= 2 samples (got " << x.size() << ")");
   MFEM_VERIFY(x.size() == y.size(),
               "PiecewiseLinear1D: x/y size mismatch (" << x.size()
               << " vs " << y.size() << ")");
   for (std::size_t k = 0; k < x.size(); ++k)
   {
      MFEM_VERIFY(std::isfinite(x[k]) && std::isfinite(y[k]),
                  "PiecewiseLinear1D: non-finite sample at index " << k);
      if (k > 0)
      {
         MFEM_VERIFY(x[k] > x[k - 1],
                     "PiecewiseLinear1D: depths must be strictly increasing; x["
                     << k << "]=" << x[k] << " <= x[" << (k - 1) << "]="
                     << x[k - 1]);
      }
   }
}

namespace
{

// Parse one depth-profile CSV into a PiecewiseLinear1D.  Each data row is
// `value, depth_km` (value FIRST, depth SECOND), comma- OR whitespace-
// separated; `#` comments and blank lines are skipped.  depth_m = depth_csv *
// depth_to_m.  Rows are sorted ascending by depth; duplicate depths abort.
// When require_positive_y (the `a` curve), every value must be > 0; the `a-b`
// curve may be negative (VW), so only finiteness is checked there.
PiecewiseLinear1D load_one_depth_csv(const std::string& path,
                                     real_t depth_to_m,
                                     bool require_positive_y,
                                     const char* which)
{
   std::ifstream ifs(path);
   MFEM_VERIFY(ifs.good(),
               "[friction.rate_state.depth_profile] cannot open " << which
               << " file '" << path << "'");

   std::vector<std::pair<real_t, real_t>> rows;   // (depth_m, value)
   std::string line;
   int lineno = 0;
   while (std::getline(ifs, line))
   {
      ++lineno;
      const auto hash = line.find('#');                 // strip inline comment
      if (hash != std::string::npos) { line.erase(hash); }
      for (char& c : line) { if (c == ',') { c = ' '; } }  // comma -> space
      std::istringstream ss(line);
      real_t value, depth_km;
      if (!(ss >> value >> depth_km)) { continue; }     // blank / comment-only
      // R-003 (Phase 11 review): read the leftover as a string so a trailing
      // NON-numeric token (e.g. "0.1 50 junk") is rejected too, not just a
      // numeric 3rd field — enforces "exactly 2 fields" as documented.
      std::string extra;
      MFEM_VERIFY(!(ss >> extra),
                  "[friction.rate_state.depth_profile] " << which << " '" << path
                  << "' line " << lineno
                  << ": expected exactly 2 fields 'value, depth_km' (unexpected "
                  "trailing token '" << extra << "')");
      rows.emplace_back(depth_km * depth_to_m, value);
   }
   MFEM_VERIFY(rows.size() >= 2,
               "[friction.rate_state.depth_profile] " << which << " '" << path
               << "' needs >= 2 data rows (got " << rows.size() << ")");

   std::sort(rows.begin(), rows.end(),
             [](const std::pair<real_t, real_t>& l,
                const std::pair<real_t, real_t>& r) { return l.first < r.first; });

   PiecewiseLinear1D pl;
   pl.x.reserve(rows.size());
   pl.y.reserve(rows.size());
   for (std::size_t k = 0; k < rows.size(); ++k)
   {
      if (k > 0)
      {
         MFEM_VERIFY(rows[k].first != rows[k - 1].first,
                     "[friction.rate_state.depth_profile] " << which << " '"
                     << path << "' has a duplicate depth " << rows[k].first
                     << " m (ambiguous interpolation)");
      }
      if (require_positive_y)
      {
         MFEM_VERIFY(rows[k].second > 0.0,
                     "[friction.rate_state.depth_profile] " << which << " '"
                     << path << "': value " << rows[k].second << " at depth "
                     << rows[k].first << " m must be > 0");
      }
      pl.x.push_back(rows[k].first);
      pl.y.push_back(rows[k].second);
   }
   pl.Validate();
   return pl;
}

}  // namespace

FrictionDepthProfile1D LoadFrictionDepthProfileCSVs(
   const FrictionDepthProfileSpec& spec)
{
   FrictionDepthProfile1D prof;
   prof.a_of_depth = load_one_depth_csv(spec.param_a_csv, spec.depth_to_m,
                                        /*require_positive_y=*/true,
                                        "param_a_csv (a)");
   prof.amb_of_depth = load_one_depth_csv(spec.param_a_minus_b_csv,
                                          spec.depth_to_m,
                                          /*require_positive_y=*/false,
                                          "param_a_minus_b_csv (a-b)");
   return prof;
}

// =====================================================================
//  TOML parsing helpers
// =====================================================================

#ifdef SEAS_USE_TOML

namespace
{

// Set of keys we've already shown the d_o deprecation notice for, so
// each tag fires exactly once per process.
struct DoNoticeOnce
{
   bool default_warned = false;
   bool rule_warned    = false;
};
static DoNoticeOnce g_do_notice_once;

real_t toml_real(const toml::value& tbl, const std::string& key,
                 real_t default_val)
{
   if (!tbl.contains(key)) { return default_val; }
   const auto& v = tbl.at(key);
   if (v.is_floating()) { return static_cast<real_t>(v.as_floating()); }
   if (v.is_integer())
   {
      return static_cast<real_t>(v.as_integer());
   }
   MFEM_ABORT("TOML key '" << key << "' must be a number (float or int)");
   return default_val;
}

// Like toml_real but returns/keeps `double` regardless of MFEM's real_t.
// Used for the QD [time] knobs whose struct fields are `double` ON PURPOSE
// (rk45_rtol=1e-50 underflows single precision); routing them through
// toml_real would narrow the `double` default through a `real_t` parameter
// and, on a single-precision build, underflow it to 0 before it reached the
// field — tripping the rk45_rtol>0 guard for every config (R-101).  Accepts
// both float and int literals, mirroring toml_real.
double toml_double(const toml::value& tbl, const std::string& key,
                   double default_val)
{
   if (!tbl.contains(key)) { return default_val; }
   const auto& v = tbl.at(key);
   if (v.is_floating()) { return v.as_floating(); }
   if (v.is_integer())  { return static_cast<double>(v.as_integer()); }
   MFEM_ABORT("TOML key '" << key << "' must be a number (float or int)");
   return default_val;
}

int toml_int(const toml::value& tbl, const std::string& key, int default_val)
{
   if (!tbl.contains(key)) { return default_val; }
   const auto& v = tbl.at(key);
   if (v.is_integer())
   {
      return static_cast<int>(v.as_integer());
   }
   MFEM_ABORT("TOML key '" << key << "' must be an integer");
   return default_val;
}

bool toml_bool(const toml::value& tbl, const std::string& key, bool default_val)
{
   if (!tbl.contains(key)) { return default_val; }
   const auto& v = tbl.at(key);
   if (v.is_boolean()) { return v.as_boolean(); }
   MFEM_ABORT("TOML key '" << key << "' must be a boolean (true/false)");
   return default_val;
}

std::string toml_str(const toml::value& tbl, const std::string& key,
                     const std::string& default_val)
{
   if (!tbl.contains(key)) { return default_val; }
   return tbl.at(key).as_string();
}

/// Read a TOML value that may be either a string ("12s", "auto") or a
/// number (12.0).  Returns the parsed seconds.
real_t toml_time_seconds(const toml::value& tbl, const std::string& key,
                         real_t default_val)
{
   if (!tbl.contains(key)) { return default_val; }
   const auto& v = tbl.at(key);
   if (v.is_string())
   {
      return SpatialTimeParseSeconds(v.as_string());
   }
   if (v.is_floating())
   {
      return static_cast<real_t>(v.as_floating());
   }
   if (v.is_integer())
   {
      return static_cast<real_t>(v.as_integer());
   }
   MFEM_ABORT("TOML key '" << key << "' must be a string (\"12s\"/\"auto\") "
              "or a number");
   return default_val;
}

/// Phase 6 req 1: read a 3-element numeric TOML array (int or float
/// elements) into a std::array<real_t,3>.  Leaves `out` untouched (keeps
/// its struct default) when the key is absent.
void toml_vec3(const toml::value& tbl, const std::string& key,
               std::array<real_t, 3>& out)
{
   if (!tbl.contains(key)) { return; }
   const auto& arr = tbl.at(key).as_array();
   MFEM_VERIFY(arr.size() == 3,
               "TOML key '" << key << "' must be a 3-element array; got "
               << arr.size() << " elements");
   for (int i = 0; i < 3; ++i)
   {
      const auto& e = arr[i];
      if      (e.is_floating()) { out[i] = static_cast<real_t>(e.as_floating()); }
      else if (e.is_integer())  { out[i] = static_cast<real_t>(e.as_integer()); }
      else { MFEM_ABORT("TOML key '" << key << "'[" << i
                        << "] must be a number (float or int)"); }
   }
}

/// Phase 6 req 1: read an integer TOML array into a std::vector<int>
/// (appends).  No-op when the key is absent.
void toml_int_array(const toml::value& tbl, const std::string& key,
                    std::vector<int>& out)
{
   if (!tbl.contains(key)) { return; }
   for (const auto& e : tbl.at(key).as_array())
   {
      MFEM_VERIFY(e.is_integer(),
                  "TOML key '" << key << "' entries must all be integers");
      out.push_back(static_cast<int>(e.as_integer()));
   }
}

/// Parse `[time].dt_initial`.  Returns the -1 sentinel ONLY when the
/// value was the literal string "auto" (or the key is absent).  Any
/// numeric or string-numeric value must be > 0; otherwise abort.  This
/// prevents the validator collision where a literal `-1.0` or `"-1.0s"`
/// would silently be treated as the "auto" sentinel (round-2 R-202).
real_t parse_dt_initial(const toml::value& t)
{
   if (!t.contains("dt_initial")) { return -1.0; }   // missing => auto
   const auto& v = t.at("dt_initial");
   if (v.is_string())
   {
      const std::string s = v.as_string();
      if (s == "auto") { return -1.0; }
      const real_t parsed = SpatialTimeParseSeconds(s);
      MFEM_VERIFY(parsed > 0.0,
                  "[time].dt_initial = \"" << s << "\" must be > 0 or "
                  "the literal string \"auto\"");
      return parsed;
   }
   real_t val = 0.0;
   if      (v.is_floating()) { val = static_cast<real_t>(v.as_floating()); }
   else if (v.is_integer())  { val = static_cast<real_t>(v.as_integer()); }
   else
   {
      MFEM_ABORT("[time].dt_initial must be a number or a string "
                 "(\"<seconds>s\" or \"auto\")");
   }
   MFEM_VERIFY(val > 0.0,
               "[time].dt_initial = " << val << " must be > 0 or the "
               "literal string \"auto\"");
   return val;
}

VelocityModel parse_velocity_model(const std::string& s)
{
   if (s == "cvmh")                 { return VelocityModel::CVMH; }
   if (s == "cvm_s4.26.m01")        { return VelocityModel::CVMS_4_26_M01; }
   if (s == "multiscale_statewise") { return VelocityModel::MultiscaleStatewise; }
   MFEM_ABORT("velocity.model must be one of "
              "{cvmh, cvm_s4.26.m01, multiscale_statewise}; got '"
              << s << "'");
   return VelocityModel::CVMH;
}

StressSourceKind parse_stress_kind(const std::string& s)
{
   if (s == "constant_tensor")      { return StressSourceKind::ConstantTensor; }
   if (s == "sidecar_hdf5")         { return StressSourceKind::SidecarHDF5; }
   if (s == "fault_local_prestress"){ return StressSourceKind::FaultLocalPrestress; }
   if (s == "depth_proportional")   { return StressSourceKind::DepthProportionalToShearModulus; }
   MFEM_ABORT("stress.kind must be one of {constant_tensor, sidecar_hdf5, "
              "fault_local_prestress, depth_proportional}; got '" << s << "'");
   return StressSourceKind::ConstantTensor;
}

void parse_spatial_rule(const toml::value& rule_tbl, SpatialRule& out,
                        bool is_lsw)
{
   // Kind (required).
   MFEM_VERIFY(rule_tbl.contains("kind"),
               "[[spatial]] rule missing required key 'kind'");
   const std::string kind_s = rule_tbl.at("kind").as_string();
   if      (kind_s == "depth")             { out.kind = SpatialRule::Kind::Depth; }
   else if (kind_s == "box")               { out.kind = SpatialRule::Kind::Box; }
   else if (kind_s == "region_attribute")  { out.kind = SpatialRule::Kind::RegionAttribute; }
   else if (kind_s == "barrier")           { out.kind = SpatialRule::Kind::Barrier; }
   else if (kind_s == "boxcar_taper")      { out.kind = SpatialRule::Kind::BoxcarTaper; }
   else
   {
      MFEM_ABORT("Unknown spatial.kind '" << kind_s
                 << "'.  Valid kinds (rev-3 + Phase 6): depth, box, "
                 << "region_attribute, barrier, boxcar_taper.  (rev-1 "
                 << "'nucleation_box' removed by D-4; use the [nucleation] "
                 << "block instead.)");
   }

   // Coordinate bounds (optional; sentinels are ±inf).
   const real_t neg_inf = -std::numeric_limits<real_t>::infinity();
   const real_t pos_inf =  std::numeric_limits<real_t>::infinity();
   out.x_min_m = toml_real(rule_tbl, "x_min_m", neg_inf);
   out.x_max_m = toml_real(rule_tbl, "x_max_m", pos_inf);
   out.y_min_m = toml_real(rule_tbl, "y_min_m", neg_inf);
   out.y_max_m = toml_real(rule_tbl, "y_max_m", pos_inf);
   out.z_min_m = toml_real(rule_tbl, "z_min_m", neg_inf);
   out.z_max_m = toml_real(rule_tbl, "z_max_m", pos_inf);
   out.region_attr = toml_int(rule_tbl, "region_attr", -1);

   const real_t nan = std::numeric_limits<real_t>::quiet_NaN();

   // Phase 6 req 5: BoxcarTaper geometry (parsed for the boxcar_taper kind).
   // Half defaults to +inf (untapered axis); center/trans default 0.
   if (out.kind == SpatialRule::Kind::BoxcarTaper)
   {
      const real_t pinf = std::numeric_limits<real_t>::infinity();
      out.boxcar_center_x_m = toml_real(rule_tbl, "boxcar_center_x_m", 0.0);
      out.boxcar_center_y_m = toml_real(rule_tbl, "boxcar_center_y_m", 0.0);
      out.boxcar_center_z_m = toml_real(rule_tbl, "boxcar_center_z_m", 0.0);
      out.boxcar_half_x_m   = toml_real(rule_tbl, "boxcar_half_x_m",   pinf);
      out.boxcar_half_y_m   = toml_real(rule_tbl, "boxcar_half_y_m",   pinf);
      out.boxcar_half_z_m   = toml_real(rule_tbl, "boxcar_half_z_m",   pinf);
      out.boxcar_trans_x_m  = toml_real(rule_tbl, "boxcar_trans_x_m",  0.0);
      out.boxcar_trans_y_m  = toml_real(rule_tbl, "boxcar_trans_y_m",  0.0);
      out.boxcar_trans_z_m  = toml_real(rule_tbl, "boxcar_trans_z_m",  0.0);
      // Half-widths and transitions must be non-negative; a non-finite half
      // is the sentinel for "this axis untapered".
      for (const auto& kv : {
           std::pair<const char*, real_t>{"boxcar_half_x_m",  out.boxcar_half_x_m},
           std::pair<const char*, real_t>{"boxcar_half_y_m",  out.boxcar_half_y_m},
           std::pair<const char*, real_t>{"boxcar_half_z_m",  out.boxcar_half_z_m},
           std::pair<const char*, real_t>{"boxcar_trans_x_m", out.boxcar_trans_x_m},
           std::pair<const char*, real_t>{"boxcar_trans_y_m", out.boxcar_trans_y_m},
           std::pair<const char*, real_t>{"boxcar_trans_z_m", out.boxcar_trans_z_m}})
      {
         MFEM_VERIFY(kv.second >= 0.0,
                     "[[spatial]] boxcar_taper rule: " << kv.first
                     << " must be >= 0; got " << kv.second);
      }
      // A boxcar_taper rule with every axis untapered (all halves +inf) is a
      // no-op taper (factor 1 everywhere) — almost certainly a config error
      // (the author meant kind=\"depth\"/\"box\").
      MFEM_VERIFY(std::isfinite(out.boxcar_half_x_m)
                  || std::isfinite(out.boxcar_half_y_m)
                  || std::isfinite(out.boxcar_half_z_m),
                  "[[spatial]] boxcar_taper rule has no finite boxcar_half_* on "
                  "any axis (all untapered) — the taper would be 1 everywhere; "
                  "set at least one boxcar_half_{x,y,z}_m or use kind=\"box\".");
   }

   if (is_lsw)
   {
      out.mu_s     = toml_real(rule_tbl, "mu_s",     nan);
      out.mu_d     = toml_real(rule_tbl, "mu_d",     nan);
      // D-3 alias: accept `d_o` in spatial rules as a synonym for `d_c`.
      if (rule_tbl.contains("d_c"))
      {
         out.d_c = toml_real(rule_tbl, "d_c", nan);
         MFEM_VERIFY(!rule_tbl.contains("d_o"),
                     "[[spatial]] rule sets both 'd_c' and 'd_o'; supply only one");
      }
      else if (rule_tbl.contains("d_o"))
      {
         out.d_c = toml_real(rule_tbl, "d_o", nan);
         if (!g_do_notice_once.rule_warned)
         {
            mfem::out << "[spatial_friction] notice: accepting deprecated "
                         "alias 'd_o' for 'd_c' in [[spatial]] rule (D-3); "
                         "this notice fires once per process.\n";
            g_do_notice_once.rule_warned = true;
         }
      }
      else
      {
         out.d_c = nan;
      }
      out.cohesion = toml_real(rule_tbl, "cohesion", nan);
      // Phase 6 req 5: cohesion-taper endpoints for a BoxcarTaper LSW rule.
      out.cohesion_inner = toml_real(rule_tbl, "cohesion_inner", nan);
      out.cohesion_outer = toml_real(rule_tbl, "cohesion_outer", nan);
      // Phase 10 (TPV31): depth-linear cohesion taper (additive, NaN-disabled).
      out.cohesion_grad_pa_per_m =
         toml_real(rule_tbl, "cohesion_grad_pa_per_m", nan);
      out.cohesion_ref_depth_m =
         toml_real(rule_tbl, "cohesion_ref_depth_m", nan);
      out.cohesion_floor_pa = toml_real(rule_tbl, "cohesion_floor_pa", 0.0);
      {
         const std::string axis_s =
            toml_str(rule_tbl, "cohesion_taper_axis", "z");
         MFEM_VERIFY(axis_s.size() == 1 &&
                     (axis_s[0] == 'x' || axis_s[0] == 'y' || axis_s[0] == 'z'),
                     "[[friction.slip_weakening.spatial]] cohesion_taper_axis "
                     "must be 'x', 'y', or 'z'; got '" << axis_s << "'");
         out.cohesion_taper_axis = axis_s[0];
      }
      // A finite grad requires a finite ref_depth (the taper origin).
      MFEM_VERIFY(std::isnan(out.cohesion_grad_pa_per_m)
                  || !std::isnan(out.cohesion_ref_depth_m),
                  "[[friction.slip_weakening.spatial]] cohesion_grad_pa_per_m "
                  "is set but cohesion_ref_depth_m is missing (the depth-linear "
                  "cohesion taper needs both).");

      // R-114: barrier-sentinel guard.  Reject any user-supplied
      // mu_s > 1e5 in a spatial-rule override.  Barrier kind is the
      // sanctioned route; the resolver assigns mu_s = 1e6 internally.
      if (!std::isnan(out.mu_s) && out.mu_s > 1.0e5)
      {
         MFEM_ABORT("[[spatial]] rule sets mu_s = " << out.mu_s
                    << " > 1.0e5 — to mark a barrier region, use "
                    << "kind = \"barrier\" instead of a large mu_s "
                    << "(R-114).");
      }
   }
   else
   {
      out.a       = toml_real(rule_tbl, "a",       nan);
      out.b       = toml_real(rule_tbl, "b",       nan);
      out.Dc      = toml_real(rule_tbl, "Dc",      nan);
      out.V_init  = toml_real(rule_tbl, "V_init",  nan);
      out.f_0     = toml_real(rule_tbl, "f_0",     nan);
      out.V_0     = toml_real(rule_tbl, "V_0",     nan);
      out.sigma_n = toml_real(rule_tbl, "sigma_n", nan);
      out.eta     = toml_real(rule_tbl, "eta",     nan);
      out.V_w     = toml_real(rule_tbl, "V_w",     nan);  // Phase 6 req 4 (SRW per-QP V_w)
   }
}

void parse_slip_weakening(const toml::value& fw_tbl, SlipWeakeningBlock& out)
{
   out.mu_s_default     = toml_real(fw_tbl, "mu_s_default",     1.1);
   out.mu_d_default     = toml_real(fw_tbl, "mu_d_default",     0.5);
   // D-3 alias for d_c_default.
   if (fw_tbl.contains("d_c_default"))
   {
      out.d_c_default = toml_real(fw_tbl, "d_c_default", 0.5);
      MFEM_VERIFY(!fw_tbl.contains("d_o_default"),
                  "[friction.slip_weakening] sets both 'd_c_default' and "
                  "'d_o_default'; supply only one");
   }
   else if (fw_tbl.contains("d_o_default"))
   {
      out.d_c_default = toml_real(fw_tbl, "d_o_default", 0.5);
      if (!g_do_notice_once.default_warned)
      {
         mfem::out << "[spatial_friction] notice: accepting deprecated "
                      "alias 'd_o_default' for 'd_c_default' (D-3); this "
                      "notice fires once per process.\n";
         g_do_notice_once.default_warned = true;
      }
   }
   else
   {
      out.d_c_default = 0.5;
   }
   out.cohesion_default = toml_real(fw_tbl, "cohesion_default", 0.0);

   if (fw_tbl.contains("spatial"))
   {
      const auto& arr = fw_tbl.at("spatial").as_array();
      out.spatial.reserve(arr.size());
      for (const auto& r : arr)
      {
         SpatialRule rule;
         parse_spatial_rule(r, rule, /*is_lsw=*/true);
         out.spatial.push_back(rule);
      }
   }

   // D-3 validator on defaults (relaxed: no upper bound on mu_s).
   MFEM_VERIFY(out.mu_d_default > 0.0,
               "[friction.slip_weakening] mu_d_default ("
               << out.mu_d_default << ") must be > 0");
   MFEM_VERIFY(out.mu_d_default < out.mu_s_default,
               "[friction.slip_weakening] mu_d_default (" << out.mu_d_default
               << ") must be < mu_s_default (" << out.mu_s_default << ")");
   MFEM_VERIFY(out.d_c_default > 0.0,
               "[friction.slip_weakening] d_c_default ("
               << out.d_c_default << ") must be > 0");
   MFEM_VERIFY(out.cohesion_default >= 0.0,
               "[friction.slip_weakening] cohesion_default ("
               << out.cohesion_default << ") must be >= 0");
   // R-114 on defaults: a user-supplied mu_s_default > 1e5 is rejected.
   MFEM_VERIFY(out.mu_s_default <= 1.0e5,
               "[friction.slip_weakening] mu_s_default ("
               << out.mu_s_default << ") > 1.0e5 — to mark every DOF as a "
               "barrier, use a [[spatial]] rule with kind = \"barrier\" "
               "(R-114).");
}

void parse_rate_state(const toml::value& rs_tbl, RateStateBlock& out)
{
   out.f_0_default     = toml_real(rs_tbl, "f_0_default",     0.6);
   out.V_0_default     = toml_real(rs_tbl, "V_0_default",     1.0e-6);
   out.a_default       = toml_real(rs_tbl, "a_default",       0.010);
   out.b_default       = toml_real(rs_tbl, "b_default",       0.015);
   out.Dc_default      = toml_real(rs_tbl, "Dc_default",      0.004);
   out.V_init_default  = toml_real(rs_tbl, "V_init_default",  1.0e-9);
   out.sigma_n_default = toml_real(rs_tbl, "sigma_n_default", 50.0e6);
   // QD Phase 7 (R-007): BP5-native per-DOF rate-state fill (bypasses resolver).
   out.bp5_analytic    = toml_bool(rs_tbl, "bp5_analytic",    false);

   // Phase 6 req 4: state-evolution selector + SRW scalars.  Default
   // "aging_law" keeps existing RS configs byte-identical.
   // REVIEW R-006: f_w_default = 0.2 (SCEC TPV104 spec, post-2026-04-24 audit;
   // the old 0.1 was half the spec value and contributed to the 1.7-2.7x
   // V_strike over-shoot — see config/tpv104_params.hpp:58).  Fed straight to
   // SlipLawSRWPsi as muW by the factory, so an SRW config that omits the key
   // must default to the corrected value.
   out.f_w_default     = toml_real(rs_tbl, "f_w_default",     0.2);
   out.V_w_default     = toml_real(rs_tbl, "V_w_default",     0.1);
   {
      const std::string se = toml_str(rs_tbl, "state_evolution", "aging_law");
      if (se == "aging_law")
      {
         out.state_evolution = StateEvolutionKind::AgingLaw;
      }
      else if (se == "slip_law_strong_rate_weakening" || se == "slip_law_srw")
      {
         out.state_evolution = StateEvolutionKind::SlipLawStrongRateWeakening;
         MFEM_VERIFY(out.V_w_default > 0.0,
                     "[friction.rate_state] state_evolution=slip_law_strong_"
                     "rate_weakening requires V_w_default > 0; got "
                     << out.V_w_default);
         // R-004: f_w_default (the weakening friction muW) feeds
         // SlipLawSRWPsi directly; a non-physical value silently produces a
         // bad steady-state friction.  Guard the SRW path (aging ignores it).
         MFEM_VERIFY(out.f_w_default > 0.0 && out.f_w_default < 1.0,
                     "[friction.rate_state] state_evolution=slip_law_strong_"
                     "rate_weakening requires f_w_default in (0,1); got "
                     << out.f_w_default);
      }
      else
      {
         MFEM_ABORT("[friction.rate_state].state_evolution must be "
                    "\"aging_law\" or \"slip_law_strong_rate_weakening\"; got '"
                    << se << "'");
      }
   }

   if (rs_tbl.contains("eta"))
   {
      const auto& v = rs_tbl.at("eta");
      if (v.is_string())
      {
         const std::string s = v.as_string();
         if (s == "auto")
         {
            out.eta_auto = true;
            out.eta_default = 0.0;
         }
         else
         {
            MFEM_ABORT("[friction.rate_state] eta must be 'auto' or a "
                       "positive real; got string '" << s << "'");
         }
      }
      else if (v.is_floating() || v.is_integer())
      {
         out.eta_auto = false;
         out.eta_default = static_cast<real_t>(
                              v.is_floating() ? v.as_floating()
                                              : v.as_integer());
         MFEM_VERIFY(out.eta_default > 0.0,
                     "[friction.rate_state] eta (" << out.eta_default
                     << ") must be > 0 when supplied as a number");
      }
      else
      {
         MFEM_ABORT("[friction.rate_state] eta must be 'auto' or a positive real");
      }
   }
   else
   {
      out.eta_auto    = true;
      out.eta_default = 0.0;
   }

   if (rs_tbl.contains("spatial"))
   {
      const auto& arr = rs_tbl.at("spatial").as_array();
      out.spatial.reserve(arr.size());
      for (const auto& r : arr)
      {
         SpatialRule rule;
         parse_spatial_rule(r, rule, /*is_lsw=*/false);
         out.spatial.push_back(rule);
      }
   }

   // Phase 11b: optional depth profile for a(z) / b(z) from two CSV files.
   // When present, the resolver seeds per-DOF a/b from the profile and the
   // scalar a_default/b_default become an unused fallback (the a<b default
   // check below is gated off).  The CSVs are READ here at parse time so a
   // missing/bad file aborts at config load, not mid-run.
   if (rs_tbl.contains("depth_profile"))
   {
      const auto& dp = rs_tbl.at("depth_profile");
      out.depth_profile.enabled = true;
      out.depth_profile.param_a_csv =
         toml_str(dp, "param_a_csv", std::string());
      out.depth_profile.param_a_minus_b_csv =
         toml_str(dp, "param_a_minus_b_csv", std::string());
      MFEM_VERIFY(!out.depth_profile.param_a_csv.empty(),
                  "[friction.rate_state.depth_profile] 'param_a_csv' is required");
      MFEM_VERIFY(!out.depth_profile.param_a_minus_b_csv.empty(),
                  "[friction.rate_state.depth_profile] 'param_a_minus_b_csv' "
                  "is required");
      const std::string units = toml_str(dp, "depth_units", std::string("km"));
      if      (units == "km") { out.depth_profile.depth_to_m = 1000.0; }
      else if (units == "m")  { out.depth_profile.depth_to_m = 1.0; }
      else
      {
         MFEM_ABORT("[friction.rate_state.depth_profile] depth_units must be "
                    "\"km\" or \"m\"; got '" << units << "'");
      }
      out.depth_profile.profile =
         LoadFrictionDepthProfileCSVs(out.depth_profile);
   }

   // Defaults validator.
   MFEM_VERIFY(out.a_default > 0.0,
               "[friction.rate_state] a_default must be > 0; got "
               << out.a_default);
   MFEM_VERIFY(out.b_default > 0.0,
               "[friction.rate_state] b_default must be > 0; got "
               << out.b_default);
   // Phase 11b: a_default < b_default governs only the scalar fallback; when a
   // depth profile is enabled the scalars are unused for a/b and the profile
   // legitimately has a > b (a-b > 0, VS) at depth, so skip this check.
   if (!out.depth_profile.enabled)
   {
      MFEM_VERIFY(out.a_default < out.b_default,
                  "[friction.rate_state] a_default (" << out.a_default
                  << ") must be < b_default (" << out.b_default << ")");
   }
   MFEM_VERIFY(out.Dc_default > 0.0,
               "[friction.rate_state] Dc_default must be > 0");
   MFEM_VERIFY(out.V_0_default > 0.0,
               "[friction.rate_state] V_0_default must be > 0");
   MFEM_VERIFY(out.sigma_n_default > 0.0,
               "[friction.rate_state] sigma_n_default must be > 0");
   MFEM_VERIFY(out.f_0_default > 0.0 && out.f_0_default < 1.0,
               "[friction.rate_state] f_0_default must be in (0, 1)");
}

SpatialFrictionConfig parse_root(const toml::value& root)
{
   SpatialFrictionConfig cfg;

   MFEM_VERIFY(root.contains("meta"),
               "Top-level [meta] block missing");
   const auto& meta = root.at("meta");
   cfg.schema_version = toml_int(meta, "schema_version", 0);
   MFEM_VERIFY(cfg.schema_version == 1,
               "[meta].schema_version must be 1; got "
               << cfg.schema_version);
   {
      const std::string law_s = toml_str(meta, "law", std::string());
      if      (law_s == "slip_weakening") { cfg.law = FrictionLawKind::SlipWeakening; }
      else if (law_s == "rate_state")     { cfg.law = FrictionLawKind::RateState; }
      else
      {
         MFEM_ABORT("[meta].law must be 'slip_weakening' or 'rate_state'; "
                    "got '" << law_s << "'");
      }
   }
   cfg.description = toml_str(meta, "description", std::string());

   // Top-level mandatory blocks per spatial_friction_config_schema.md
   // §"Top-level layout".  Without these guards, an omitted block would
   // silently default to the struct defaults (e.g. zero pre-stress for
   // a missing [stress] block, tfinal=12s for a missing [time] block).
   MFEM_VERIFY(root.contains("material_constant_fallback"),
               "Top-level [material_constant_fallback] block missing");
   MFEM_VERIFY(root.contains("pore_pressure"),
               "Top-level [pore_pressure] block missing");
   MFEM_VERIFY(root.contains("mesh"),
               "Top-level [mesh] block missing");
   MFEM_VERIFY(root.contains("velocity"),
               "Top-level [velocity] block missing");
   MFEM_VERIFY(root.contains("stress"),
               "Top-level [stress] block missing");
   MFEM_VERIFY(root.contains("numerics"),
               "Top-level [numerics] block missing");
   MFEM_VERIFY(root.contains("time"),
               "Top-level [time] block missing");
   MFEM_VERIFY(root.contains("output"),
               "Top-level [output] block missing");

   if (root.contains("material_constant_fallback"))
   {
      const auto& m = root.at("material_constant_fallback");
      cfg.material_fallback.lambda = toml_real(m, "lambda", 32.0e9);
      cfg.material_fallback.mu     = toml_real(m, "mu",     32.0e9);
      cfg.material_fallback.rho    = toml_real(m, "rho",    2670.0);
   }
   // Material-fallback validator (D-3 doesn't relax these).
   MFEM_VERIFY(cfg.material_fallback.mu > 0.0,
               "[material_constant_fallback].mu must be > 0");
   MFEM_VERIFY(cfg.material_fallback.rho > 0.0,
               "[material_constant_fallback].rho must be > 0");
   MFEM_VERIFY(cfg.material_fallback.lambda + 2.0 * cfg.material_fallback.mu > 0.0,
               "[material_constant_fallback] (lambda + 2*mu) must be > 0");

   if (root.contains("pore_pressure"))
   {
      const auto& pp = root.at("pore_pressure");
      cfg.stress.pore_pressure.P_p_pa = toml_real(pp, "P_p_pa", 0.0);
      cfg.stress.pore_pressure.P_p_grad_pa_per_m = toml_real(pp, "P_p_grad_pa_per_m", 0.0);
      cfg.stress.pore_pressure.min_sigma_n_pa = toml_real(pp, "min_sigma_n_pa", 0.0);
   }
   MFEM_VERIFY(cfg.stress.pore_pressure.min_sigma_n_pa >= 0.0,
               "[pore_pressure].min_sigma_n_pa must be >= 0");
   if (cfg.stress.pore_pressure.P_p_grad_pa_per_m < 0.0)
   {
      mfem::out << "[spatial_friction] warning: [pore_pressure].P_p_grad_pa_per_m = "
                << cfg.stress.pore_pressure.P_p_grad_pa_per_m
                << " is negative (pore pressure decreases with depth) — "
                << "unusual but permitted.\n";
   }

   if (root.contains("mesh"))
   {
      const auto& m = root.at("mesh");
      cfg.mesh.path  = toml_str(m, "path", std::string());
      cfg.mesh.order = toml_int(m, "order", 1);
      cfg.mesh.lc_far_m = toml_real(m, "lc_far_m", -1.0);   // Phase 12.2 (PML thickness derivation)
   }
   MFEM_VERIFY(!cfg.mesh.path.empty(), "[mesh].path must be non-empty");
   MFEM_VERIFY(cfg.mesh.order >= 1,    "[mesh].order must be >= 1");

   if (root.contains("velocity"))
   {
      const auto& v = root.at("velocity");
      cfg.velocity.model = parse_velocity_model(
                              toml_str(v, "model", "cvmh"));
      cfg.velocity.dataset_root = toml_str(v, "dataset_root", std::string());
      cfg.velocity.override_path = toml_str(v, "override_path", std::string());
      cfg.velocity.use_sidecar  = toml_bool(v, "use_sidecar", true);
   }
   // When use_sidecar = false the driver consumes [material_constant_fallback]
   // and never reads model/dataset_root/override_path; relax the non-empty
   // check in that case.
   if (cfg.velocity.use_sidecar)
   {
      MFEM_VERIFY(!cfg.velocity.dataset_root.empty()
                  || !cfg.velocity.override_path.empty(),
                  "[velocity] must set either 'dataset_root' "
                  "(for model-based resolution) or 'override_path' when "
                  "'use_sidecar' is true (the default)");
   }

   if (root.contains("stress"))
   {
      const auto& s = root.at("stress");
      MFEM_VERIFY(s.contains("kind"),
                  "[stress].kind is required; must be "
                  "\"constant_tensor\", \"sidecar_hdf5\", "
                  "\"fault_local_prestress\", or \"depth_proportional\"");
      cfg.stress.kind = parse_stress_kind(toml_str(s, "kind", "constant_tensor"));
      // QD-only (R-002): BP5-analytic per-DOF prestress flag (no enum change).
      cfg.stress.bp5_analytic = toml_bool(s, "bp5_analytic", false);

      const bool has_sxx = s.contains("sigma_xx_pa");
      const bool has_syy = s.contains("sigma_yy_pa");
      const bool has_szz = s.contains("sigma_zz_pa");
      const bool has_sxy = s.contains("sigma_xy_pa");
      const bool has_syz = s.contains("sigma_yz_pa");
      const bool has_sxz = s.contains("sigma_xz_pa");
      const bool has_path = s.contains("sidecar_path")
                            && !toml_str(s, "sidecar_path", std::string()).empty();
      // R-004: fault-local-prestress keys; rejected by the other two kinds.
      const bool has_flp = s.contains("tau_strike_pa")
                           || s.contains("tau_dip_pa")
                           || s.contains("sigma_n_pa");

      if (cfg.stress.kind == StressSourceKind::ConstantTensor)
      {
         // bp5_analytic builds its own per-DOF Cauchy tensor from bp5_params
         // (tau0_vec + effective sigma_n) in the driver; the six sigma_*_pa
         // keys are provably unused on that path, so do not require them.
         MFEM_VERIFY(cfg.stress.bp5_analytic
                     || (has_sxx && has_syy && has_szz
                         && has_sxy && has_syz && has_sxz),
                     "[stress] kind=\"constant_tensor\" requires all six "
                     "sigma_*_pa keys to be present (unless bp5_analytic=true)");
         MFEM_VERIFY(!has_path,
                     "[stress] kind=\"constant_tensor\" must NOT set "
                     "sidecar_path (it is for kind=\"sidecar_hdf5\" only)");
         MFEM_VERIFY(!has_flp,
                     "[stress] kind=\"constant_tensor\" must NOT set "
                     "tau_strike_pa/tau_dip_pa/sigma_n_pa (those are for "
                     "kind=\"fault_local_prestress\" only)");
         cfg.stress.sigma_xx_pa = toml_real(s, "sigma_xx_pa", 0.0);
         cfg.stress.sigma_yy_pa = toml_real(s, "sigma_yy_pa", 0.0);
         cfg.stress.sigma_zz_pa = toml_real(s, "sigma_zz_pa", 0.0);
         cfg.stress.sigma_xy_pa = toml_real(s, "sigma_xy_pa", 0.0);
         cfg.stress.sigma_yz_pa = toml_real(s, "sigma_yz_pa", 0.0);
         cfg.stress.sigma_xz_pa = toml_real(s, "sigma_xz_pa", 0.0);
      }
      else if (cfg.stress.kind == StressSourceKind::FaultLocalPrestress)
      {
         // D3.2 (Phase 6): fault-local background pre-stress, right-lateral /
         // compression POSITIVE.  Seeded DIRECTLY (no Cauchy projection) via
         // FaultGeometry::ComputeParamsFaultLocal: tau2_0=tau_strike,
         // tau1_0=tau_dip, sigma_n0=sigma_n_pa - P_p.
         MFEM_VERIFY(!has_sxx && !has_syy && !has_szz
                     && !has_sxy && !has_syz && !has_sxz,
                     "[stress] kind=\"fault_local_prestress\" must NOT set any "
                     "Cauchy sigma_*_pa key (those are for "
                     "kind=\"constant_tensor\")");
         MFEM_VERIFY(!has_path,
                     "[stress] kind=\"fault_local_prestress\" must NOT set "
                     "sidecar_path (it is for kind=\"sidecar_hdf5\" only)");
         MFEM_VERIFY(s.contains("sigma_n_pa"),
                     "[stress] kind=\"fault_local_prestress\" requires "
                     "sigma_n_pa (compression POSITIVE)");
         // R-003: tau_strike is the shear driver; a missing key would silently
         // seed a zero-shear (locked) fault.
         MFEM_VERIFY(s.contains("tau_strike_pa"),
                     "[stress] kind=\"fault_local_prestress\" requires "
                     "tau_strike_pa (the shear driver; right-lateral POSITIVE)");
         cfg.stress.tau_strike_pa = toml_real(s, "tau_strike_pa", 0.0);
         cfg.stress.tau_dip_pa    = toml_real(s, "tau_dip_pa", 0.0);
         cfg.stress.sigma_n_pa    = toml_real(s, "sigma_n_pa", 0.0);
         MFEM_VERIFY(cfg.stress.sigma_n_pa > 0.0,
                     "[stress] kind=\"fault_local_prestress\" sigma_n_pa must "
                     "be > 0 (compression positive); got "
                     << cfg.stress.sigma_n_pa);
         // R-002 (req 2): optional rectangular tau_strike patches
         // ([[stress.patch]], last-match-wins).  center_* default NaN /
         // half_* default +inf => that axis is unconstrained.  tau_strike_pa
         // is required per patch (a patch with no override is meaningless).
         if (s.contains("patch"))
         {
            const auto& arr = s.at("patch").as_array();
            cfg.stress.fault_local_patches.reserve(arr.size());
            for (const auto& pn : arr)
            {
               FaultLocalPatch p;
               p.center_x_m = toml_real(pn, "center_x_m",
                                        std::numeric_limits<real_t>::quiet_NaN());
               p.center_y_m = toml_real(pn, "center_y_m",
                                        std::numeric_limits<real_t>::quiet_NaN());
               p.center_z_m = toml_real(pn, "center_z_m",
                                        std::numeric_limits<real_t>::quiet_NaN());
               p.half_x_m   = toml_real(pn, "half_x_m",
                                        std::numeric_limits<real_t>::infinity());
               p.half_y_m   = toml_real(pn, "half_y_m",
                                        std::numeric_limits<real_t>::infinity());
               p.half_z_m   = toml_real(pn, "half_z_m",
                                        std::numeric_limits<real_t>::infinity());
               p.tau_strike_pa = toml_real(pn, "tau_strike_pa",
                                           std::numeric_limits<real_t>::quiet_NaN());
               MFEM_VERIFY(!std::isnan(p.tau_strike_pa),
                           "[[stress.patch]] requires tau_strike_pa "
                           "(right-lateral POSITIVE)");
               cfg.stress.fault_local_patches.push_back(p);
            }
         }
      }
      else if (cfg.stress.kind ==
               StressSourceKind::DepthProportionalToShearModulus)
      {
         // Phase 10 (TPV31): depth-proportional Cauchy tensor, scaled per-point
         // by mu(point)/mu_ref.  Components are entered in MPa in the
         // [stress.depth_proportional] sub-table and converted to Pa here.
         MFEM_VERIFY(!has_sxx && !has_syy && !has_szz
                     && !has_sxy && !has_syz && !has_sxz,
                     "[stress] kind=\"depth_proportional\" must NOT set the "
                     "Cauchy sigma_*_pa keys (use [stress.depth_proportional] "
                     "with sigma_*_per_mu in MPa instead)");
         MFEM_VERIFY(!has_path && !has_flp,
                     "[stress] kind=\"depth_proportional\" must NOT set "
                     "sidecar_path or tau_*_pa/sigma_n_pa keys");
         MFEM_VERIFY(s.contains("depth_proportional"),
                     "[stress] kind=\"depth_proportional\" requires a "
                     "[stress.depth_proportional] sub-table with the six "
                     "sigma_*_per_mu components (MPa) + mu_ref_pa");
         const auto& dp = s.at("depth_proportional");
         const real_t mpa = 1.0e6;   // TOML values are in MPa
         auto& d = cfg.stress.depth_proportional;
         d.sigma_xx_per_mu = mpa * toml_real(dp, "sigma_xx_per_mu", 0.0);
         d.sigma_yy_per_mu = mpa * toml_real(dp, "sigma_yy_per_mu", 0.0);
         d.sigma_zz_per_mu = mpa * toml_real(dp, "sigma_zz_per_mu", 0.0);
         d.sigma_xy_per_mu = mpa * toml_real(dp, "sigma_xy_per_mu", 0.0);
         d.sigma_yz_per_mu = mpa * toml_real(dp, "sigma_yz_per_mu", 0.0);
         d.sigma_xz_per_mu = mpa * toml_real(dp, "sigma_xz_per_mu", 0.0);
         d.mu_ref_pa       = toml_real(dp, "mu_ref_pa", 32.03812032e9);
         MFEM_VERIFY(d.mu_ref_pa > 0.0,
                     "[stress.depth_proportional].mu_ref_pa must be > 0; got "
                     << d.mu_ref_pa);
      }
      else
      {
         MFEM_VERIFY(!has_sxx && !has_syy && !has_szz
                     && !has_sxy && !has_syz && !has_sxz,
                     "[stress] kind=\"sidecar_hdf5\" must NOT set any of "
                     "sigma_xx_pa..sigma_xz_pa (those are for "
                     "kind=\"constant_tensor\" only)");
         MFEM_VERIFY(has_path,
                     "[stress] kind=\"sidecar_hdf5\" requires a non-empty "
                     "sidecar_path");
         MFEM_VERIFY(!has_flp,
                     "[stress] kind=\"sidecar_hdf5\" must NOT set "
                     "tau_strike_pa/tau_dip_pa/sigma_n_pa (those are for "
                     "kind=\"fault_local_prestress\" only)");
         cfg.stress.sidecar_path = toml_str(s, "sidecar_path", std::string());
      }
   }

   if (root.contains("numerics"))
   {
      const auto& n = root.at("numerics");
      cfg.numerics.ader_order = toml_int(n, "ader_order", 2);
      cfg.numerics.mixed_flux = toml_str(n, "mixed_flux", "none");
      cfg.numerics.cfl        = toml_real(n, "cfl", 0.5);
      cfg.numerics.use_pml    = toml_bool(n, "use_pml", false);

      // Phase 6 req 3: cfl_safety / fault_iterator / interior_flux selectors.
      // REVIEW R-002: cfl_safety defaults to "dg" (NOT "raw") so a config that
      // omits the key keeps the driver's long-standing always-DG-factored
      // behavior (every SAFS config omits it).  The driver now honors the
      // selector (CflSafetyFactor), so a default of "raw" would have silently
      // un-factored every existing SAFS run (~9x larger dt).  "raw" is the
      // explicit experimental escape hatch.  fault_iterator/interior_flux keep
      // their one-shot/scalar defaults (current behavior).
      const std::string cs = toml_str(n, "cfl_safety", "dg");
      if      (cs == "raw") { cfg.numerics.cfl_safety = CflSafety::Raw; }
      else if (cs == "dg")  { cfg.numerics.cfl_safety = CflSafety::Dg; }
      else { MFEM_ABORT("[numerics].cfl_safety must be \"raw\" or \"dg\"; got '"
                        << cs << "'"); }

      // REVIEW R-001: default "substep" (NOT "one-shot") — the spatial driver
      // always sub-steps and MFEM_VERIFYs against "one-shot", so a config that
      // omits the key (all 8 SAFS configs) must default to the supported mode.
      const std::string fi = toml_str(n, "fault_iterator", "substep");
      if      (fi == "one-shot") { cfg.numerics.fault_iterator = FaultIteratorKind::OneShot; }
      else if (fi == "substep")  { cfg.numerics.fault_iterator = FaultIteratorKind::Substep; }
      else { MFEM_ABORT("[numerics].fault_iterator must be \"one-shot\" or "
                        "\"substep\"; got '" << fi << "'"); }

      const std::string ifx = toml_str(n, "interior_flux", "scalar");
      if      (ifx == "scalar") { cfg.numerics.interior_flux = InteriorFlux::Scalar; }
      else if (ifx == "matrix") { cfg.numerics.interior_flux = InteriorFlux::Matrix; }
      else { MFEM_ABORT("[numerics].interior_flux must be \"scalar\" or "
                        "\"matrix\"; got '" << ifx << "'"); }

      // Phase 14: time-integrator selector (default "ader" — byte-exact for
      // every existing config, which omits the key).  rk4 / rk45 select the
      // explicit coupled Runge–Kutta stepper (the driver enforces the
      // rate_state + scalar-interior-flux preconditions at setup).
      const std::string ti = toml_str(n, "time_integrator", "ader");
      if      (ti == "ader") { cfg.numerics.time_integrator = TimeIntegratorKind::ADER; }
      else if (ti == "rk4")  { cfg.numerics.time_integrator = TimeIntegratorKind::RK4; }
      else if (ti == "rk45") { cfg.numerics.time_integrator = TimeIntegratorKind::RK45; }
      else { MFEM_ABORT("[numerics].time_integrator must be one of "
                        "{ader, rk4, rk45}; got '" << ti << "'"); }

      // Phase 12.2: PML knobs (only consulted when use_pml; defaults derive
      // the thickness from pml_cells*lc_far and damp the bottom but not the
      // free surface).  Every existing config omits these — the defaults
      // reproduce a sensible SAFS half-space PML.
      cfg.numerics.pml_thickness_m = toml_real(n, "pml_thickness_m", -1.0);
      cfg.numerics.pml_target_R    = toml_real(n, "pml_target_R",    1.0e-3);
      cfg.numerics.pml_cells       = toml_int (n, "pml_cells",       4);
      cfg.numerics.pml_damp_bottom = toml_bool(n, "pml_damp_bottom", true);
      cfg.numerics.pml_damp_top    = toml_bool(n, "pml_damp_top",    false);
   }
   MFEM_VERIFY(cfg.numerics.ader_order >= 1,
               "[numerics].ader_order must be >= 1");
   MFEM_VERIFY(cfg.numerics.mixed_flux == "none"
               || cfg.numerics.mixed_flux == "adjacent"
               || cfg.numerics.mixed_flux == "all_continuous",
               "[numerics].mixed_flux must be one of "
               "{none, adjacent, all_continuous}; got '"
               << cfg.numerics.mixed_flux << "'");
   MFEM_VERIFY(cfg.numerics.cfl > 0.0 && cfg.numerics.cfl < 1.0,
               "[numerics].cfl must be in (0, 1)");
   // Phase 6 req 3 mutual-exclusion (R-1203 sibling): mixed-flux is a
   // scalar-path-only optimization, so interior_flux="matrix" forbids it.
   MFEM_VERIFY(cfg.numerics.interior_flux == InteriorFlux::Scalar
               || cfg.numerics.mixed_flux == "none",
               "[numerics] interior_flux=\"matrix\" is incompatible with "
               "mixed_flux=\"" << cfg.numerics.mixed_flux << "\" (mixed-flux "
               "is valid only on the scalar interior-flux path); set "
               "mixed_flux=\"none\" when using matrix.");
   // NOTE: the companion guard "interior_flux=matrix requires material.kind !=
   // Constant" is enforced after the [material] block is parsed below (it
   // needs cfg.material.kind, which is read further down in parse_root).

   // Phase 12.2: PML knob validators.  These hold unconditionally (the
   // defaults pass); the driver additionally aborts at construction if PML
   // is enabled but the thickness cannot be derived (lc_far_m unset and no
   // explicit thickness) or the shell would reach the fault.
   MFEM_VERIFY(cfg.numerics.pml_target_R > 0.0 && cfg.numerics.pml_target_R < 1.0,
               "[numerics].pml_target_R must be in (0,1); got "
               << cfg.numerics.pml_target_R);
   MFEM_VERIFY(cfg.numerics.pml_cells >= 1,
               "[numerics].pml_cells must be >= 1; got " << cfg.numerics.pml_cells);
   MFEM_VERIFY(cfg.numerics.pml_thickness_m < 0.0
               || cfg.numerics.pml_thickness_m > 0.0,
               "[numerics].pml_thickness_m must be > 0 (or < 0 to derive it "
               "from pml_cells*lc_far_m); got " << cfg.numerics.pml_thickness_m);

   if (root.contains("time"))
   {
      const auto& t = root.at("time");
      cfg.time.tfinal     = toml_time_seconds(t, "tfinal",     12.0);
      cfg.time.t_initial  = toml_real        (t, "t_initial",  0.0);
      cfg.time.dt_initial = parse_dt_initial (t);
      cfg.time.dt_max     = toml_time_seconds(t, "dt_max",     0.1);

      // Phase 1 (QD) knobs — read only by the spatial_seas (quasi-dynamic)
      // driver; absent keys keep the struct defaults.  Plain numbers (no
      // unit-string parsing) per the plan: dt_init [s], dt_max_years [yr],
      // plate_rate_vp [m/s], rk45_{atol,rtol}; use_petsc_ts bool.  Read via
      // toml_double (NOT toml_real) so the double defaults survive on a
      // single-precision MFEM build (R-101 — rk45_rtol=1e-50 < float min).
      cfg.time.rk45_atol     = toml_double(t, "rk45_atol",     cfg.time.rk45_atol);
      cfg.time.rk45_rtol     = toml_double(t, "rk45_rtol",     cfg.time.rk45_rtol);
      cfg.time.dt_init       = toml_double(t, "dt_init",       cfg.time.dt_init);
      cfg.time.dt_max_years  = toml_double(t, "dt_max_years",  cfg.time.dt_max_years);
      cfg.time.plate_rate_vp = toml_double(t, "plate_rate_vp", cfg.time.plate_rate_vp);
      cfg.time.use_petsc_ts  = toml_bool  (t, "use_petsc_ts",  cfg.time.use_petsc_ts);
   }
   MFEM_VERIFY(cfg.time.tfinal > 0.0,
               "[time].tfinal must be > 0; got " << cfg.time.tfinal);
   MFEM_VERIFY(cfg.time.dt_max > 0.0,
               "[time].dt_max must be > 0; got " << cfg.time.dt_max);
   MFEM_VERIFY(cfg.time.t_initial >= 0.0,
               "[time].t_initial must be >= 0; got " << cfg.time.t_initial);
   // (dt_initial is fully validated inside parse_dt_initial — it is
   // either > 0 or exactly the -1 "auto" sentinel set by the literal
   // string "auto".  Literal numeric -1.0 is rejected per R-202.)
   // QD knob validation (Phase 1).  dt_init / plate_rate_vp use a <0
   // "derive at runtime" sentinel, so only their POSITIVE range matters;
   // the tolerances and dt ceiling must be strictly positive.
   MFEM_VERIFY(cfg.time.rk45_atol > 0.0,
               "[time].rk45_atol must be > 0; got " << cfg.time.rk45_atol);
   MFEM_VERIFY(cfg.time.rk45_rtol > 0.0,
               "[time].rk45_rtol must be > 0; got " << cfg.time.rk45_rtol);
   MFEM_VERIFY(cfg.time.dt_max_years > 0.0,
               "[time].dt_max_years must be > 0; got " << cfg.time.dt_max_years);

   // Phase 1 (QD): optional [solver] block (linear solver + AMG/KSP knobs).
   // Read only by the spatial_seas (quasi-dynamic) driver; a missing block
   // keeps the SolverSpec defaults (type="cg_amg").  Additive — the dynamic
   // driver does not read [solver].
   if (root.contains("solver"))
   {
      const auto& s = root.at("solver");
      cfg.solver.type                   = toml_str (s, "type",                   cfg.solver.type);
      cfg.solver.ksp_rtol               = toml_real(s, "ksp_rtol",               cfg.solver.ksp_rtol);
      cfg.solver.ksp_atol               = toml_real(s, "ksp_atol",               cfg.solver.ksp_atol);
      cfg.solver.ksp_maxit              = toml_int (s, "ksp_maxit",              cfg.solver.ksp_maxit);
      cfg.solver.amg_elasticity_options = toml_bool(s, "amg_elasticity_options", cfg.solver.amg_elasticity_options);
      cfg.solver.amg_relax_type         = toml_int (s, "amg_relax_type",         cfg.solver.amg_relax_type);
      cfg.solver.amg_aggressive_levels  = toml_int (s, "amg_aggressive_levels",  cfg.solver.amg_aggressive_levels);
      cfg.solver.amg_print_level        = toml_int (s, "amg_print_level",        cfg.solver.amg_print_level);
      cfg.solver.blr_tol                = toml_real(s, "blr_tol",                cfg.solver.blr_tol);
      cfg.solver.residual_check         = toml_bool(s, "residual_check",         cfg.solver.residual_check);
   }
   // Validate solver.type at parse (Edge Case: unknown ⇒ hard error, never
   // silently defaulted).  ParseQDSolverType aborts on an unrecognized value;
   // this runs even when [solver] is absent so the default "cg_amg" is
   // confirmed valid.
   //
   // R-102: this is INTENTIONAL shared-parser validation.  It executes for
   // every spatial config, including the dynamic-rupture driver's, but is
   // dyn-safe: no spatial-schema config carries a [solver] table (the BP5
   // [solver] tables are the seas_driver schema, parsed elsewhere, and use
   // the key `solver_type`, not `type`), so cfg.solver.type is always the
   // default "cg_amg" there and never aborts.  Kept in the parser (not the QD
   // driver) because the plan's Edge Case mandates the hard error AT PARSE.
   (void)ParseQDSolverType(cfg.solver.type);
   MFEM_VERIFY(cfg.solver.ksp_rtol >= 0.0,
               "[solver].ksp_rtol must be >= 0; got " << cfg.solver.ksp_rtol);
   MFEM_VERIFY(cfg.solver.ksp_atol >= 0.0,
               "[solver].ksp_atol must be >= 0; got " << cfg.solver.ksp_atol);
   MFEM_VERIFY(cfg.solver.ksp_maxit >= 1,
               "[solver].ksp_maxit must be >= 1; got " << cfg.solver.ksp_maxit);
   MFEM_VERIFY(cfg.solver.blr_tol >= 0.0,
               "[solver].blr_tol must be >= 0; got " << cfg.solver.blr_tol);

   if (root.contains("output"))
   {
      const auto& o = root.at("output");
      cfg.output.output_dir              = toml_str (o, "output_dir",              std::string());
      cfg.output.restart_prefix          = toml_str (o, "restart_prefix",          "cp");
      // Parity Phase 2: default-flip to "off" for volume + bulk;
      // fault stays "hdf5" for backwards compatibility with existing
      // SAFS workflows that depend on fault output.
      cfg.output.paraview_volume         = toml_str (o, "paraview_volume",         "off");
      cfg.output.paraview_bulk           = toml_str (o, "paraview_bulk",           "off");
      cfg.output.paraview_fault          = toml_str (o, "paraview_fault",          "hdf5");
      cfg.output.paraview_volume_dt      = toml_time_seconds(o, "paraview_volume_dt", 0.05);
      cfg.output.paraview_bulk_dt        = toml_time_seconds(o, "paraview_bulk_dt",   0.05);
      cfg.output.paraview_fault_dt       = toml_time_seconds(o, "paraview_fault_dt",  0.001);
      cfg.output.paraview_volume_zfp_tol = toml_real(o, "paraview_volume_zfp_tol", 1e-3);
      cfg.output.paraview_bulk_zfp_tol   = toml_real(o, "paraview_bulk_zfp_tol",   1e-3);
      cfg.output.paraview_fault_zfp_tol  = toml_real(o, "paraview_fault_zfp_tol",  1e-12);
      cfg.output.max_snapshots           = toml_int (o, "max_snapshots",           5000);
      cfg.output.checkpoint_every_steps  = toml_int (o, "checkpoint_every_steps",  10000);
      cfg.output.bp5_stations            = toml_bool(o, "bp5_stations",            false);

      // Parity Phase 2: 9 new fields (master gate + every-step + legacy
      // ASCII + per-collection deflate levels + regime-adaptive cadences).
      cfg.output.paraview_enabled              = toml_bool(o, "paraview_enabled",              false);
      cfg.output.paraview_every_steps          = toml_int (o, "paraview_every_steps",          0);
      cfg.output.paraview_fault_legacy_ascii   = toml_bool(o, "paraview_fault_legacy_ascii",   false);
      cfg.output.paraview_volume_deflate_level = toml_int (o, "paraview_volume_deflate_level", -1);
      cfg.output.paraview_bulk_deflate_level   = toml_int (o, "paraview_bulk_deflate_level",   -1);
      cfg.output.paraview_fault_deflate_level  = toml_int (o, "paraview_fault_deflate_level",  -1);
      cfg.output.paraview_coseismic_dt         = toml_time_seconds(o, "paraview_coseismic_dt",    -1.0);
      cfg.output.paraview_nucleation_dt        = toml_time_seconds(o, "paraview_nucleation_dt",   -1.0);
      cfg.output.paraview_interseismic_dt      = toml_time_seconds(o, "paraview_interseismic_dt", -1.0);
   }
   MFEM_VERIFY(!cfg.output.output_dir.empty(),
               "[output].output_dir must be non-empty");
   MFEM_VERIFY(!cfg.output.restart_prefix.empty(),
               "[output].restart_prefix must be non-empty");
   for (const auto& mode_name : {
        std::pair<std::string, std::string>{"paraview_volume", cfg.output.paraview_volume},
        std::pair<std::string, std::string>{"paraview_bulk",   cfg.output.paraview_bulk},
        std::pair<std::string, std::string>{"paraview_fault",  cfg.output.paraview_fault}})
   {
      MFEM_VERIFY(mode_name.second == "hdf5"
                  || mode_name.second == "vtu"
                  || mode_name.second == "off",
                  "[output]." << mode_name.first << " must be one of "
                  "{hdf5, vtu, off}; got '" << mode_name.second << "'");
   }
   MFEM_VERIFY(cfg.output.paraview_volume_zfp_tol >= 0.0
               && cfg.output.paraview_bulk_zfp_tol >= 0.0
               && cfg.output.paraview_fault_zfp_tol >= 0.0,
               "[output] paraview_*_zfp_tol must all be >= 0");
   MFEM_VERIFY(cfg.output.paraview_volume_dt > 0.0,
               "[output].paraview_volume_dt must be > 0; got "
               << cfg.output.paraview_volume_dt);
   MFEM_VERIFY(cfg.output.paraview_bulk_dt > 0.0,
               "[output].paraview_bulk_dt must be > 0; got "
               << cfg.output.paraview_bulk_dt);
   MFEM_VERIFY(cfg.output.paraview_fault_dt > 0.0,
               "[output].paraview_fault_dt must be > 0; got "
               << cfg.output.paraview_fault_dt);
   MFEM_VERIFY(cfg.output.max_snapshots >= 1,
               "[output].max_snapshots must be >= 1");
   MFEM_VERIFY(cfg.output.checkpoint_every_steps >= 1,
               "[output].checkpoint_every_steps must be >= 1");
   // Parity Phase 2 validators.
   MFEM_VERIFY(cfg.output.paraview_every_steps >= 0,
               "[output].paraview_every_steps must be >= 0; got "
               << cfg.output.paraview_every_steps);
   for (const auto& kv : {
        std::pair<std::string, int>{"paraview_volume_deflate_level",
                                    cfg.output.paraview_volume_deflate_level},
        std::pair<std::string, int>{"paraview_bulk_deflate_level",
                                    cfg.output.paraview_bulk_deflate_level},
        std::pair<std::string, int>{"paraview_fault_deflate_level",
                                    cfg.output.paraview_fault_deflate_level}})
   {
      MFEM_VERIFY(kv.second == -1 || (kv.second >= 0 && kv.second <= 9),
                  "[output]." << kv.first
                  << " must be -1 (none) or in [0..9]; got " << kv.second);
   }
   for (const auto& kv : {
        std::pair<std::string, real_t>{"paraview_coseismic_dt",
                                       cfg.output.paraview_coseismic_dt},
        std::pair<std::string, real_t>{"paraview_nucleation_dt",
                                       cfg.output.paraview_nucleation_dt},
        std::pair<std::string, real_t>{"paraview_interseismic_dt",
                                       cfg.output.paraview_interseismic_dt}})
   {
      MFEM_VERIFY(kv.second < 0.0 || kv.second > 0.0,
                  "[output]." << kv.first
                  << " must be unset (negative) or > 0; got " << kv.second);
   }

   if (root.contains("nucleation"))
   {
      const auto& nuc = root.at("nucleation");
      cfg.nucleation.enabled = true;

      const std::string kind_s = toml_str(nuc, "kind", "");
      if (kind_s == "gradual_overstress")
      {
         cfg.nucleation.kind = NucleationKind::GradualOverstress;
         MFEM_VERIFY(nuc.contains("gradual_overstress"),
                     "[nucleation] kind=\"gradual_overstress\" requires a "
                     "[nucleation.gradual_overstress] sub-block");
         const auto& g = nuc.at("gradual_overstress");
         auto& gs = cfg.nucleation.gradual_overstress;
         gs.center_x_m          = toml_real(g, "center_x_m",          0.0);
         gs.center_y_m          = toml_real(g, "center_y_m",          0.0);
         gs.center_z_m          = toml_real(g, "center_z_m",          0.0);
         gs.radius_dip_m        = toml_real(g, "radius_dip_m",        0.0);
         gs.radius_strike_m     = toml_real(g, "radius_strike_m",     0.0);
         gs.delta_tau_dip_pa    = toml_real(g, "delta_tau_dip_pa",    0.0);
         gs.delta_tau_strike_pa = toml_real(g, "delta_tau_strike_pa", 0.0);
         gs.T_nuc_s             = toml_time_seconds(g, "T_nuc_s",     0.0);
         MFEM_VERIFY(gs.radius_dip_m    > 0.0,
                     "[nucleation.gradual_overstress].radius_dip_m must be > 0; "
                     "got " << gs.radius_dip_m);
         MFEM_VERIFY(gs.radius_strike_m > 0.0,
                     "[nucleation.gradual_overstress].radius_strike_m must be > 0; "
                     "got " << gs.radius_strike_m);
         MFEM_VERIFY(gs.T_nuc_s         > 0.0,
                     "[nucleation.gradual_overstress].T_nuc_s must be > 0; "
                     "got " << gs.T_nuc_s);
      }
      else if (kind_s == "gradual_overstress_compact_circular")
      {
         // Phase 6 req 6 (TPV102/104).  Config-only; the resolver is Phase 7.
         cfg.nucleation.kind = NucleationKind::GradualOverstressCompactCircular;
         MFEM_VERIFY(nuc.contains("gradual_overstress_compact_circular"),
                     "[nucleation] kind=\"gradual_overstress_compact_circular\" "
                     "requires a [nucleation.gradual_overstress_compact_circular]"
                     " sub-block");
         const auto& g = nuc.at("gradual_overstress_compact_circular");
         auto& cc = cfg.nucleation.compact_circular;
         cc.center_x_m   = toml_real(g, "center_x_m",   0.0);
         cc.center_y_m   = toml_real(g, "center_y_m",   0.0);
         cc.center_z_m   = toml_real(g, "center_z_m",   0.0);
         cc.radius_m     = toml_real(g, "radius_m",     0.0);
         cc.delta_tau_pa = toml_real(g, "delta_tau_pa", 0.0);
         cc.T_nuc_s      = toml_time_seconds(g, "T_nuc_s", 0.0);
         MFEM_VERIFY(cc.radius_m > 0.0,
                     "[nucleation.gradual_overstress_compact_circular].radius_m "
                     "must be > 0; got " << cc.radius_m);
         MFEM_VERIFY(cc.T_nuc_s  > 0.0,
                     "[nucleation.gradual_overstress_compact_circular].T_nuc_s "
                     "must be > 0; got " << cc.T_nuc_s);
      }
      else if (kind_s == "instantaneous_overstress_circular")
      {
         // Phase 6 req 6 (TPV31).  Config-only; the applicator is Phase 7.
         cfg.nucleation.kind = NucleationKind::InstantaneousOverstressCircular;
         MFEM_VERIFY(nuc.contains("instantaneous_overstress_circular"),
                     "[nucleation] kind=\"instantaneous_overstress_circular\" "
                     "requires a [nucleation.instantaneous_overstress_circular]"
                     " sub-block");
         const auto& g = nuc.at("instantaneous_overstress_circular");
         auto& ic = cfg.nucleation.instantaneous_circular;
         ic.center_x_m   = toml_real(g, "center_x_m",   0.0);
         ic.center_y_m   = toml_real(g, "center_y_m",   0.0);
         ic.center_z_m   = toml_real(g, "center_z_m",   0.0);
         ic.radius_m     = toml_real(g, "radius_m",     0.0);
         ic.taper_m      = toml_real(g, "taper_m",      0.0);
         ic.delta_tau_pa = toml_real(g, "delta_tau_pa", 0.0);
         // Phase 10 (TPV31 spec p. 7): optional per-DOF mu(depth)/mu_ref
         // amplitude scaling.  Default 0.0 ⇒ disabled (uniform delta_tau_pa);
         // TPV31 sets it to the spec reference modulus mu_0.
         ic.mu_ref_pa    = toml_real(g, "mu_ref_pa",    0.0);
         MFEM_VERIFY(ic.radius_m > 0.0,
                     "[nucleation.instantaneous_overstress_circular].radius_m "
                     "must be > 0; got " << ic.radius_m);
         MFEM_VERIFY(ic.taper_m >= 0.0,
                     "[nucleation.instantaneous_overstress_circular].taper_m "
                     "must be >= 0; got " << ic.taper_m);
         MFEM_VERIFY(ic.mu_ref_pa >= 0.0,
                     "[nucleation.instantaneous_overstress_circular].mu_ref_pa "
                     "must be >= 0 (0 disables mu-scaling); got "
                     << ic.mu_ref_pa);
      }
      else
      {
         MFEM_ABORT("[nucleation].kind must be one of {gradual_overstress, "
                    "gradual_overstress_compact_circular, "
                    "instantaneous_overstress_circular}; got '" << kind_s << "'");
      }
   }
   // else: enabled stays false; driver runs without nucleation perturbation.

   // Friction-law-specific block.  Both-or-neither rule (Validator §3).
   const bool has_lsw = root.contains("friction")
                        && root.at("friction").contains("slip_weakening");
   const bool has_rs  = root.contains("friction")
                        && root.at("friction").contains("rate_state");

   if (cfg.law == FrictionLawKind::SlipWeakening)
   {
      MFEM_VERIFY(has_lsw,
                  "[meta].law=\"slip_weakening\" requires "
                  "[friction.slip_weakening] block");
      MFEM_VERIFY(!has_rs,
                  "[meta].law=\"slip_weakening\" forbids "
                  "[friction.rate_state] block");
      SlipWeakeningBlock blk;
      parse_slip_weakening(root.at("friction").at("slip_weakening"), blk);
      cfg.slip_weakening = blk;
   }
   else
   {
      MFEM_VERIFY(has_rs,
                  "[meta].law=\"rate_state\" requires "
                  "[friction.rate_state] block");
      MFEM_VERIFY(!has_lsw,
                  "[meta].law=\"rate_state\" forbids "
                  "[friction.slip_weakening] block");
      RateStateBlock blk;
      parse_rate_state(root.at("friction").at("rate_state"), blk);
      cfg.rate_state = blk;
   }

   // Optional top-level `[friction].sigma_n_strength_floor_pa` (the σ_n
   // strength floor; sliver-blowup plan 2026-05-26).  Applies to BOTH
   // laws.  `[friction]` is guaranteed present here (one of the two law
   // sub-blocks above was validated as required).  Default -1.0 ⇒
   // disabled (each law keeps its exact current strength expression).
   // When the key IS present it must be a finite value >= 0 — a negative
   // floor is the disabled sentinel, not a usable value, so an explicit
   // negative is a configuration error.
   {
      const auto& fr = root.at("friction");
      // R-001: the floor is a TOP-LEVEL `[friction]` key, but the per-law
      // sub-block parsers do NOT reject unknown keys — so a floor mistakenly
      // nested under `[friction.slip_weakening]` / `[friction.rate_state]`
      // (the natural place, since every other friction parameter lives there)
      // would be SILENTLY ignored, disabling the safety feature with no error.
      // Catch that mis-nesting loudly.
      for (const char* sub : {"slip_weakening", "rate_state"})
      {
         if (fr.contains(sub)
             && fr.at(sub).contains("sigma_n_strength_floor_pa"))
         {
            MFEM_ABORT("[friction." << sub
                       << "].sigma_n_strength_floor_pa is mis-placed: the "
                       "strength floor is a TOP-LEVEL [friction] key (it "
                       "applies to both laws), not a per-law key.  Move it "
                       "directly under [friction].");
         }
      }
      cfg.sigma_n_strength_floor_pa =
         toml_real(fr, "sigma_n_strength_floor_pa", -1.0);
      if (fr.contains("sigma_n_strength_floor_pa"))
      {
         MFEM_VERIFY(std::isfinite(cfg.sigma_n_strength_floor_pa)
                     && cfg.sigma_n_strength_floor_pa >= 0.0,
                     "[friction].sigma_n_strength_floor_pa must be a finite "
                     "value >= 0 (Pa); a negative value is the disabled "
                     "sentinel and must be expressed by OMITTING the key, "
                     "not by setting it negative.  Got "
                     << cfg.sigma_n_strength_floor_pa);
         // R-003: no upper bound is enforced (the floor's usable range is
         // problem-dependent), but a floor far above any physical
         // cohesion-like value almost certainly indicates an exponent typo
         // (e.g. 10.0e9 for 10.0e6).  Because the strength uses
         // max(sigma_n, floor), a floor exceeding the ambient normal stress
         // silently LOCKS the fault (no rupture) — warn so the typo is caught.
         if (cfg.sigma_n_strength_floor_pa > 1.0e9)
         {
            mfem::out << "[spatial_friction] WARNING: "
                      << "[friction].sigma_n_strength_floor_pa = "
                      << cfg.sigma_n_strength_floor_pa << " Pa (> 1 GPa) is far "
                      << "above any physical cohesion-like floor; a floor "
                      << "exceeding the ambient normal stress LOCKS the fault "
                      << "(no rupture).  Did you mean "
                      << (cfg.sigma_n_strength_floor_pa / 1.0e3) << " Pa?\n";
         }
      }
   }

   // REVIEW R-005: σ_n double-source consistency.  For a
   // fault_local_prestress RS config, [stress].sigma_n_pa (seeded into
   // geom.sigma_n_per_dof by ComputeParamsFaultLocal, MINUS P_p) and
   // [friction.rate_state].sigma_n_default are two independent inputs for
   // the same effective normal stress.  The resolver uses the seeded value
   // (sigma_n_pa − P_p); sigma_n_default is the fallback when no seed is
   // passed.  Require them to agree so a future edit to one cannot silently
   // desync the RS resolver/initial-ψ fallback from the seeded σ_n.  Gated
   // on fault_local_prestress so SAFS (constant_tensor / sidecar) configs —
   // whose σ_n comes from a projection, not sigma_n_pa — are exempt.
   if (cfg.stress.kind == StressSourceKind::FaultLocalPrestress
       && cfg.rate_state.has_value())
   {
      const real_t sigma_n_eff =
         cfg.stress.sigma_n_pa - cfg.stress.pore_pressure.P_p_pa;
      MFEM_VERIFY(std::abs(cfg.rate_state->sigma_n_default - sigma_n_eff)
                  <= 1e-6 * cfg.stress.sigma_n_pa,
                  "[friction.rate_state].sigma_n_default ("
                  << cfg.rate_state->sigma_n_default << ") must equal "
                  "[stress].sigma_n_pa - [pore_pressure].P_p_pa ("
                  << sigma_n_eff << ") for a fault_local_prestress "
                  "rate-state config (σ_n double-source consistency, R-005).");
   }

   // =====================================================================
   //  Phase 6 req 1: optional TPV config blocks ([problem], [boundary],
   //  [fault_geometry], [hypocenter], [material]).  Each is OPTIONAL —
   //  an absent block keeps the struct defaults, so all existing SAFS
   //  TOMLs (which set none of these) parse unchanged.  req-7 guards
   //  validate each block only when it is present.
   // =====================================================================

   if (root.contains("problem"))
   {
      // Informational tag only (no code branches on it).
      cfg.problem.tag = toml_str(root.at("problem"), "tag", std::string());
   }

   if (root.contains("boundary"))
   {
      const auto& b = root.at("boundary");
      cfg.boundary.fault_attr = toml_int(b, "fault_attr", -1);
      toml_int_array(b, "natural_attrs",   cfg.boundary.natural_attrs);
      toml_int_array(b, "absorbing_attrs", cfg.boundary.absorbing_attrs);
      // Phase 2 (QD): far-field Dirichlet plate-loading walls (optional;
      // read only by the quasi-dynamic driver).  Additive — absent key ⇒
      // empty, so dynamic-rupture configs are unaffected.
      toml_int_array(b, "dirichlet_attrs", cfg.boundary.dirichlet_attrs);
      // QD far-field plate-loading function (quasi-dynamic driver only).
      cfg.boundary.plate_loading = toml_str(b, "plate_loading", "bp5");
      MFEM_VERIFY(cfg.boundary.plate_loading == "bp5" ||
                  cfg.boundary.plate_loading == "saf_recenter_y",
                  "[boundary].plate_loading must be \"bp5\" or "
                  "\"saf_recenter_y\"; got '"
                  << cfg.boundary.plate_loading << "'");

      // req-7 guards.  fault_attr must be a positive mesh attribute (0 is
      // the "interior" sentinel in MFEM; a missing/negative fault_attr would
      // silently disable fault detection).
      MFEM_VERIFY(cfg.boundary.fault_attr > 0,
                  "[boundary].fault_attr must be a positive mesh attribute; got "
                  << cfg.boundary.fault_attr);
      auto require_positive = [](const std::vector<int>& v, const char* name)
      {
         for (int a : v)
         {
            MFEM_VERIFY(a > 0, "[boundary]." << name << " entries must be "
                        "positive mesh attributes; got " << a);
         }
      };
      require_positive(cfg.boundary.natural_attrs,   "natural_attrs");
      require_positive(cfg.boundary.absorbing_attrs, "absorbing_attrs");
      require_positive(cfg.boundary.dirichlet_attrs, "dirichlet_attrs");

      // Disjointness: an attribute assigned to two different roles (e.g. both
      // natural and absorbing, or fault and natural) is an unresolvable BC
      // conflict — abort rather than silently pick one.
      auto has = [](const std::vector<int>& v, int a)
      { return std::find(v.begin(), v.end(), a) != v.end(); };
      MFEM_VERIFY(!has(cfg.boundary.natural_attrs,   cfg.boundary.fault_attr)
                  && !has(cfg.boundary.absorbing_attrs, cfg.boundary.fault_attr)
                  && !has(cfg.boundary.dirichlet_attrs, cfg.boundary.fault_attr),
                  "[boundary].fault_attr (" << cfg.boundary.fault_attr
                  << ") must not also appear in natural_attrs/absorbing_attrs/"
                  "dirichlet_attrs");
      for (int a : cfg.boundary.natural_attrs)
      {
         MFEM_VERIFY(!has(cfg.boundary.absorbing_attrs, a),
                     "[boundary] natural_attrs and absorbing_attrs must be "
                     "disjoint; attribute " << a << " appears in both");
         MFEM_VERIFY(!has(cfg.boundary.dirichlet_attrs, a),
                     "[boundary] natural_attrs and dirichlet_attrs must be "
                     "disjoint; attribute " << a << " appears in both");
      }
      for (int a : cfg.boundary.absorbing_attrs)
      {
         MFEM_VERIFY(!has(cfg.boundary.dirichlet_attrs, a),
                     "[boundary] absorbing_attrs and dirichlet_attrs must be "
                     "disjoint; attribute " << a << " appears in both");
      }
   }

   if (root.contains("fault_geometry"))
   {
      const auto& f = root.at("fault_geometry");
      toml_vec3(f, "ref_normal", cfg.fault_geometry.ref_normal);
      toml_vec3(f, "up",         cfg.fault_geometry.up);
      cfg.fault_geometry.kind = toml_str(f, "kind", std::string());

      // req-7 guards: ref_normal and up define the fault-local frame; the
      // strike axis is t2 = normalize(up × ref_normal), so up must not be
      // parallel to ref_normal (else the cross product is zero / frame
      // undefined).  R-003: rather than demand an EXACT unit-norm (which would
      // reject legitimately hand-entered direction vectors — e.g. a 45°
      // dipping-fault normal [0.577,0.577,0.577] is off-unit by ~6e-4), accept
      // any non-degenerate vector and NORMALIZE it in place.  Reject only a
      // near-zero vector, which carries no direction.
      auto norm3 = [](const std::array<real_t, 3>& v)
      { return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]); };
      const real_t nn = norm3(cfg.fault_geometry.ref_normal);
      const real_t un = norm3(cfg.fault_geometry.up);
      MFEM_VERIFY(nn > 1e-8,
                  "[fault_geometry].ref_normal must be a non-zero direction "
                  "vector; |ref_normal| = " << nn);
      MFEM_VERIFY(un > 1e-8,
                  "[fault_geometry].up must be a non-zero direction vector; "
                  "|up| = " << un);
      for (int k = 0; k < 3; ++k)
      {
         cfg.fault_geometry.ref_normal[k] /= nn;
         cfg.fault_geometry.up[k]         /= un;
      }
      const auto& n = cfg.fault_geometry.ref_normal;   // unit now
      const auto& u = cfg.fault_geometry.up;           // unit now
      const real_t dot = n[0]*u[0] + n[1]*u[1] + n[2]*u[2];
      MFEM_VERIFY(std::abs(dot) <= 1.0 - 1e-9,
                  "[fault_geometry].up must not be parallel to ref_normal "
                  "(|ref_normal·up| = " << std::abs(dot) << " ~ 1); the strike "
                  "axis t2 = normalize(up x ref_normal) would be undefined");
   }

   if (root.contains("hypocenter"))
   {
      const auto& h = root.at("hypocenter");
      cfg.hypocenter.x_m = toml_real(h, "x_m", 0.0);
      cfg.hypocenter.y_m = toml_real(h, "y_m", 0.0);
      cfg.hypocenter.z_m = toml_real(h, "z_m", 0.0);
      cfg.hypocenter.nucleation_radius_m = toml_real(h, "nucleation_radius_m", 0.0);
      cfg.hypocenter.nucleation_taper_m  = toml_real(h, "nucleation_taper_m",  0.0);

      MFEM_VERIFY(cfg.hypocenter.nucleation_radius_m >= 0.0,
                  "[hypocenter].nucleation_radius_m must be >= 0; got "
                  << cfg.hypocenter.nucleation_radius_m);
      MFEM_VERIFY(cfg.hypocenter.nucleation_taper_m >= 0.0,
                  "[hypocenter].nucleation_taper_m must be >= 0; got "
                  << cfg.hypocenter.nucleation_taper_m);
      // R-008 guard: with an up-pointing vertical axis (up[2] > 0, z increases
      // upward) a hypocenter must lie at or below the free surface (z <= 0).
      // A positive z would place the nucleation patch in the air.
      if (cfg.fault_geometry.up[2] > 0.0)
      {
         MFEM_VERIFY(cfg.hypocenter.z_m <= 0.0,
                     "[hypocenter].z_m must be <= 0 when [fault_geometry].up[2] "
                     "> 0 (z increases upward, free surface at z = 0); got z_m = "
                     << cfg.hypocenter.z_m);
      }
   }

   if (root.contains("material"))
   {
      const auto& m = root.at("material");
      const std::string mk = toml_str(m, "kind", "constant");
      if      (mk == "constant")         { cfg.material.kind = MaterialKind::Constant; }
      else if (mk == "depth_profile_1d") { cfg.material.kind = MaterialKind::DepthProfile1D; }
      else if (mk == "sidecar_hdf5")     { cfg.material.kind = MaterialKind::SidecarHDF5; }
      else
      {
         MFEM_ABORT("[material].kind must be one of {constant, "
                    "depth_profile_1d, sidecar_hdf5}; got '" << mk << "'");
      }
      cfg.material.profile_csv  = toml_str(m, "profile_csv",  std::string());
      cfg.material.sidecar_path = toml_str(m, "sidecar_path", std::string());

      if (cfg.material.kind == MaterialKind::DepthProfile1D)
      {
         // Phase 10 (TPV31): inline depth-profile layers + depth axis.
         {
            const std::string axis_s = toml_str(m, "depth_axis", "z");
            MFEM_VERIFY(axis_s.size() == 1 &&
                        (axis_s[0] == 'x' || axis_s[0] == 'y' || axis_s[0] == 'z'),
                        "[material].depth_axis must be 'x', 'y', or 'z'; got '"
                        << axis_s << "'");
            cfg.material.depth_axis = axis_s[0];
         }
         MFEM_VERIFY(root.contains("material_profile"),
                     "[material].kind=\"depth_profile_1d\" requires a "
                     "[[material_profile.layer]] array of layer tables");
         const auto& mp = root.at("material_profile");
         MFEM_VERIFY(mp.contains("layer"),
                     "[material_profile] must contain a `layer` array "
                     "(TOML syntax `[[material_profile.layer]]`)");
         const auto& layers_arr = mp.at("layer").as_array();
         MFEM_VERIFY(!layers_arr.empty(),
                     "[material_profile.layer] must have at least one entry");
         cfg.material.profile_layers.clear();
         cfg.material.profile_layers.reserve(layers_arr.size());
         for (std::size_t i = 0; i < layers_arr.size(); ++i)
         {
            const auto& L = layers_arr[i];
            DepthProfileLayer dpl;
            dpl.depth_top_m = toml_real(L, "depth_top_m", 0.0);
            dpl.depth_bot_m = toml_real(L, "depth_bot_m", 0.0);
            dpl.vp_ms       = toml_real(L, "vp_ms",       0.0);
            dpl.vs_ms       = toml_real(L, "vs_ms",       0.0);
            dpl.rho_kgm3    = toml_real(L, "rho_kgm3",    0.0);
            dpl.interp      = toml_str (L, "interp",      "constant");
            MFEM_VERIFY(dpl.vp_ms > 0.0,
                        "[[material_profile.layer]] " << i
                        << ".vp_ms must be > 0; got " << dpl.vp_ms);
            MFEM_VERIFY(dpl.vs_ms > 0.0,
                        "[[material_profile.layer]] " << i
                        << ".vs_ms must be > 0; got " << dpl.vs_ms);
            MFEM_VERIFY(dpl.rho_kgm3 > 0.0,
                        "[[material_profile.layer]] " << i
                        << ".rho_kgm3 must be > 0; got " << dpl.rho_kgm3);
            MFEM_VERIFY(dpl.interp == "constant" || dpl.interp == "linear",
                        "[[material_profile.layer]] " << i
                        << ".interp must be 'constant' or 'linear'; got '"
                        << dpl.interp << "'");
            cfg.material.profile_layers.push_back(dpl);
         }
      }
      if (cfg.material.kind == MaterialKind::SidecarHDF5)
      {
         MFEM_VERIFY(!cfg.material.sidecar_path.empty(),
                     "[material] kind=\"sidecar_hdf5\" requires a non-empty "
                     "sidecar_path");
      }
   }

   // Phase 6 req 3 (deferred guard, completed now that [material] is parsed):
   // the matrix interior-flux path needs per-element material contrast to be
   // meaningful — a spatially-Constant material has no contrast, so matrix is
   // a wasteful no-op there and almost certainly a config mistake.
   MFEM_VERIFY(cfg.numerics.interior_flux == InteriorFlux::Scalar
               || cfg.material.kind != MaterialKind::Constant,
               "[numerics].interior_flux=\"matrix\" requires a non-Constant "
               "[material].kind (depth_profile_1d or sidecar_hdf5); a Constant "
               "material has no element-to-element contrast for the matrix flux "
               "to resolve.");

   return cfg;
}

}  // namespace

#endif  // SEAS_USE_TOML

// =====================================================================
//  LoadSpatialFrictionConfig / ParseSpatialFrictionConfigString
// =====================================================================

SpatialFrictionConfig LoadSpatialFrictionConfig(const std::string& toml_path)
{
#ifdef SEAS_USE_TOML
   try
   {
      const toml::value root = toml::parse(toml_path);
      return parse_root(root);
   }
   catch (const std::exception& e)
   {
      MFEM_ABORT("LoadSpatialFrictionConfig('" << toml_path
                 << "'): toml parse error: " << e.what());
      return {};
   }
#else
   (void)toml_path;
   MFEM_ABORT("LoadSpatialFrictionConfig: SEAS_USE_TOML is not defined.  "
              "Rebuild with extern/toml11/toml.hpp present (the Makefile "
              "auto-enables SEAS_USE_TOML when the file exists).");
   return {};
#endif
}

SpatialFrictionConfig ParseSpatialFrictionConfigString(const std::string& toml_text)
{
#ifdef SEAS_USE_TOML
   try
   {
      std::istringstream iss(toml_text);
      const toml::value root = toml::parse(iss, "<inline-test>");
      return parse_root(root);
   }
   catch (const std::exception& e)
   {
      MFEM_ABORT("ParseSpatialFrictionConfigString: toml parse error: "
                 << e.what());
      return {};
   }
#else
   (void)toml_text;
   MFEM_ABORT("ParseSpatialFrictionConfigString: SEAS_USE_TOML not defined.");
   return {};
#endif
}

// =====================================================================
//  ParseQDSolverType  (Phase 1, QD)
// =====================================================================
// Maps a [solver].type string to mfem::seas::SolverType.  No TOML needed,
// so it is defined unconditionally (outside the SEAS_USE_TOML guard).  An
// unrecognized value is a hard error — the QD config schema must never
// silently fall back to a solver the user did not request.
mfem::seas::SolverType ParseQDSolverType(const std::string& type)
{
   if (type == "cg_amg")    { return SolverType::CG_AMG; }
   if (type == "gmres_amg") { return SolverType::GMRES_AMG; }
   if (type == "gmres_ilu") { return SolverType::GMRES_BlockILU; }
   if (type == "mumps")     { return SolverType::MUMPS; }
   if (type == "mumps_blr") { return SolverType::MUMPS_BLR; }
   if (type == "superlu")   { return SolverType::SUPERLU; }
   if (type == "strumpack") { return SolverType::STRUMPACK; }
   MFEM_ABORT("ParseQDSolverType: unknown [solver].type '" << type
              << "'.  Accepted: cg_amg | gmres_amg | gmres_ilu | mumps | "
                 "mumps_blr | superlu | strumpack.");
   return SolverType::CG_AMG;  // unreachable (MFEM_ABORT does not return)
}

// =====================================================================
//  SpatialFrictionResolver
// =====================================================================

SlipWeakeningPerDOFParams SpatialFrictionResolver::ResolveSlipWeakening(
   const SlipWeakeningBlock& cfg,
   const Vector&             dof_coords_3d,
   const Array<int>&         dof_to_attr) const
{
   const int N = dof_coords_3d.Size() / 3;
   MFEM_VERIFY(dof_coords_3d.Size() == 3 * N,
               "ResolveSlipWeakening: dof_coords_3d.Size() must be 3 * N; "
               "got " << dof_coords_3d.Size());
   MFEM_VERIFY(dof_to_attr.Size() == N,
               "ResolveSlipWeakening: dof_to_attr.Size() ("
               << dof_to_attr.Size() << ") must equal num_dofs (" << N << ")");

   // Validate every default is not NaN (R-014).
   MFEM_VERIFY(!std::isnan(cfg.mu_s_default)
               && !std::isnan(cfg.mu_d_default)
               && !std::isnan(cfg.d_c_default)
               && !std::isnan(cfg.cohesion_default),
               "ResolveSlipWeakening: every LSW default must be a number "
               "(no NaN sentinels at the block level)");

   // R-001 (Phase 6 req-5 deferral, justified): the boxcar_taper rule kind and
   // its SCECBoxcar/BoxcarTaperFactor helpers are CONFIG-LEVEL ONLY this phase.
   // req 5 specifies the kind + the helper functions but NOT how the resolver
   // consumes the taper — the blend semantics and which parameters taper are
   // unspecified, and req 6's sibling new kinds are likewise "config-only;
   // applicator wired in a later phase".  The generic per-DOF override loop
   // below would otherwise apply a matching boxcar_taper rule as a HARD region
   // across the entire boxcar+transition footprint (matches()>0), silently
   // dropping cohesion_inner/cohesion_outer and producing no taper — a
   // wrong-physics trap.  Reject it explicitly until the consumption phase.
   for (const auto& r : cfg.spatial)
   {
      MFEM_VERIFY(r.kind != SpatialRule::Kind::BoxcarTaper,
                  "ResolveSlipWeakening: a 'boxcar_taper' spatial rule is not "
                  "yet consumed by the resolver (Phase 6 req 5 is config-only — "
                  "the kind + SCECBoxcar/BoxcarTaperFactor exist and parse, but "
                  "the per-DOF cohesion/parameter taper blend is wired in a "
                  "later phase).  Remove the boxcar_taper rule or use "
                  "kind=\"box\"/\"depth\" for now.");
   }

   SlipWeakeningPerDOFParams p;
   p.mu_s.SetSize(N);
   p.mu_d.SetSize(N);
   p.d_c.SetSize(N);
   p.cohesion.SetSize(N);

   for (int i = 0; i < N; ++i)
   {
      const real_t x = dof_coords_3d(3*i + 0);
      const real_t y = dof_coords_3d(3*i + 1);
      const real_t z = dof_coords_3d(3*i + 2);

      real_t mu_s_i = cfg.mu_s_default;
      real_t mu_d_i = cfg.mu_d_default;
      real_t d_c_i  = cfg.d_c_default;
      real_t coh_i  = cfg.cohesion_default;

      for (const auto& r : cfg.spatial)
      {
         if (!r.matches(x, y, z, dof_to_attr[i])) { continue; }
         if (r.kind == SpatialRule::Kind::Barrier)
         {
            // R-114 sentinel: locks the DOF.  Users never type 1.0e6.
            mu_s_i = 1.0e6;
            // Leave mu_d, d_c, cohesion at their pre-barrier values:
            // the barrier short-circuit in
            // LSWFrictionCoefficient_ForcedRupture only reads mu_s.
            continue;
         }
         if (!std::isnan(r.mu_s))     { mu_s_i = r.mu_s; }
         if (!std::isnan(r.mu_d))     { mu_d_i = r.mu_d; }
         if (!std::isnan(r.d_c))      { d_c_i  = r.d_c; }
         if (!std::isnan(r.cohesion)) { coh_i  = r.cohesion; }
         // Phase 10 (TPV31): depth-linear cohesion taper.  Takes precedence
         // over any constant `cohesion` set above; clamped from below by
         // `cohesion_floor_pa` (typically 0).  NaN grad ⇒ disabled (other
         // LSW configs unchanged).  Ported verbatim from hrs-ref.
         if (!std::isnan(r.cohesion_grad_pa_per_m))
         {
            MFEM_VERIFY(!std::isnan(r.cohesion_ref_depth_m),
                        "ResolveSlipWeakening: cohesion_grad_pa_per_m set "
                        "but cohesion_ref_depth_m missing in spatial rule");
            // For axis 'z' (canonical SEAS), interpret the axis value as
            // POSITIVE depth (depth = max(0, -z)) so the taper formula
            // `cohesion = floor + grad * (ref_depth - depth)` can be
            // written with POSITIVE `ref_depth`/`grad` regardless of mesh-z
            // sign — matches MakeDepthProfile1DMaterial.  Axes 'x'/'y' keep
            // the literal "axis_val is used directly" semantic.
            real_t axis_val = (r.cohesion_taper_axis == 'x') ? x :
                              (r.cohesion_taper_axis == 'z') ? z : y;
            if (r.cohesion_taper_axis == 'z')
            {
               const real_t d = -axis_val;
               axis_val = (d > 0.0) ? d : 0.0;
            }
            const real_t raw = r.cohesion_floor_pa
               + r.cohesion_grad_pa_per_m
                 * (r.cohesion_ref_depth_m - axis_val);
            coh_i = (raw > r.cohesion_floor_pa) ? raw : r.cohesion_floor_pa;
         }
      }

      // Per-DOF validator.  Barrier DOFs have mu_s = 1e6 ≫ mu_d, so
      // the inequality is automatically satisfied.
      MFEM_VERIFY(mu_d_i > 0.0,
                  "ResolveSlipWeakening: mu_d <= 0 at DOF " << i
                  << " (got " << mu_d_i << ")");
      MFEM_VERIFY(mu_d_i < mu_s_i,
                  "ResolveSlipWeakening: mu_d (" << mu_d_i << ") >= mu_s ("
                  << mu_s_i << ") at DOF " << i);
      MFEM_VERIFY(d_c_i > 0.0,
                  "ResolveSlipWeakening: d_c <= 0 at DOF " << i
                  << " (got " << d_c_i << ")");
      MFEM_VERIFY(coh_i >= 0.0,
                  "ResolveSlipWeakening: cohesion < 0 at DOF " << i
                  << " (got " << coh_i << ")");

      p.mu_s(i)     = mu_s_i;
      p.mu_d(i)     = mu_d_i;
      p.d_c(i)      = d_c_i;
      p.cohesion(i) = coh_i;
   }

   return p;
}

namespace
{

template <typename MeshT>
RateStatePerDOFParams resolve_rs_impl(
   const RateStateBlock&     cfg,
   const Vector&             dof_coords_3d,
   const Array<int>&         dof_to_elem,
   const Array<int>&         dof_to_attr,
   const MaterialField&      material,
   MeshT&                    mesh,
   const PorePressureSpec&   pp,
   const Vector&             sigma_n_total_per_dof)
{
   const int N = dof_coords_3d.Size() / 3;
   MFEM_VERIFY(dof_coords_3d.Size() == 3 * N,
               "ResolveRateState: dof_coords_3d.Size() must be 3 * N");
   MFEM_VERIFY(dof_to_elem.Size() == N,
               "ResolveRateState: dof_to_elem.Size() != N");
   MFEM_VERIFY(dof_to_attr.Size() == N,
               "ResolveRateState: dof_to_attr.Size() != N");
   MFEM_VERIFY(sigma_n_total_per_dof.Size() == N || sigma_n_total_per_dof.Size() == 0,
               "ResolveRateState: sigma_n_total_per_dof.Size() must be 0 or N");

   // R-001 (Phase 6 req-5 deferral, justified): boxcar_taper is config-only
   // this phase (see the same guard in ResolveSlipWeakening).  The per-DOF loop
   // below would silently apply a matching boxcar_taper rule's a/b/Dc/... as a
   // HARD region over the whole boxcar+transition footprint; reject it until
   // the consumption phase wires the taper.
   for (const auto& r : cfg.spatial)
   {
      MFEM_VERIFY(r.kind != SpatialRule::Kind::BoxcarTaper,
                  "ResolveRateState: a 'boxcar_taper' spatial rule is not yet "
                  "consumed by the resolver (Phase 6 req 5 is config-only — the "
                  "kind + SCECBoxcar/BoxcarTaperFactor exist and parse, but the "
                  "per-DOF parameter taper blend is wired in a later phase).  "
                  "Remove the boxcar_taper rule or use kind=\"box\"/\"depth\" "
                  "for now.");
   }

   // R-001 (Phase 11 review, CRITICAL): per-DOF b is wired through the AGING
   // iterator (Tpv102SubStepIterator / RateStateAgingPolicy read d.b) and the
   // equilibrium seed (SeedEquilibriumPsi_RS uses rs.b(ii)) ONLY.  The slip-law
   // strong-rate-weakening (SRW) policy evolves psi with the SCALAR L.GetB()
   // (= blk.b_default; tpv104 left untouched, friction/state_policies.hpp).  If
   // a non-scalar b reaches the SRW path, the seed uses per-DOF b while the
   // dynamics use scalar b -> the fault is seeded OUT of the SRW iterator's
   // equilibrium at t=0 (silent disequilibrium / spurious transient).  Pre-
   // Phase-11 the seed used the scalar b too, so this is a Phase-11-introduced
   // regression that only bites once b is non-scalar.  Reject the combination
   // loudly until SRW is given a per-DOF b (out of Phase-11 scope, per the plan
   // §11a step 2 / §Constraints).  Aging is unaffected (it honours d.b).
   if (cfg.state_evolution == StateEvolutionKind::SlipLawStrongRateWeakening)
   {
      MFEM_VERIFY(!cfg.depth_profile.enabled,
                  "ResolveRateState: [friction.rate_state.depth_profile] "
                  "(per-DOF b) is NOT supported with state_evolution="
                  "slip_law_strong_rate_weakening; the SRW iterator evolves psi "
                  "with the scalar b_default (Phase 11 wired per-DOF b through "
                  "the aging iterator only), so a per-DOF b would seed psi out "
                  "of the SRW iterator's equilibrium at t=0.  Use the aging law, "
                  "or remove the depth profile and use a scalar b.");
      for (const auto& r : cfg.spatial)
      {
         MFEM_VERIFY(std::isnan(r.b),
                     "ResolveRateState: a per-DOF 'b' spatial override is NOT "
                     "supported with state_evolution="
                     "slip_law_strong_rate_weakening (the SRW iterator uses the "
                     "scalar b_default; per-DOF b is aging-law only).  Remove "
                     "the 'b' override from the spatial rule, or use the aging "
                     "law.");
      }
   }

   RateStatePerDOFParams p;
   p.a.SetSize(N);   p.b.SetSize(N);   p.Dc.SetSize(N);
   p.V_init.SetSize(N);
   p.f_0.SetSize(N); p.V_0.SetSize(N); p.eta.SetSize(N);
   p.sigma_n_eff.SetSize(N);
   p.V_w.SetSize(N);   // Phase 6 req 4 (SRW per-DOF weakening velocity)

   for (int i = 0; i < N; ++i)
   {
      const real_t x = dof_coords_3d(3*i + 0);
      const real_t y = dof_coords_3d(3*i + 1);
      const real_t z = dof_coords_3d(3*i + 2);

      // Vertical depth (metres) for the profile + pore pressure (z<0 below).
      const real_t depth_i = std::max(static_cast<real_t>(0.0), -z);

      // Seed defaults.  Phase 11b: when a depth profile is enabled, a/b come
      // from the profile (a = a(depth), b = a - (a-b)); otherwise the scalars.
      real_t a_i;
      real_t b_i;
      if (cfg.depth_profile.enabled)
      {
         a_i = cfg.depth_profile.profile.a(depth_i);
         b_i = cfg.depth_profile.profile.b(depth_i);
      }
      else
      {
         a_i = cfg.a_default;
         b_i = cfg.b_default;
      }
      real_t Dc_i     = cfg.Dc_default;
      real_t V_init_i = cfg.V_init_default;
      real_t f0_i     = cfg.f_0_default;
      real_t V0_i     = cfg.V_0_default;
      real_t eta_i    = cfg.eta_default;
      real_t sn_i     = cfg.sigma_n_default;
      real_t V_w_i    = cfg.V_w_default;   // Phase 6 req 4 (SRW)

      bool   eta_rule_set = false;

      for (const auto& r : cfg.spatial)
      {
         if (!r.matches(x, y, z, dof_to_attr[i])) { continue; }
         // RS resolver ignores Barrier (LSW-only concept).
         if (r.kind == SpatialRule::Kind::Barrier) { continue; }
         // Phase 11a relaxes R-006: per-DOF b is now supported.  DOFData.b
         // routes it through the Tpv102SubStepIterator ψ-update + the
         // equilibrium-ψ seed, so the original "knob that does nothing" reason
         // no longer holds for b.  f_0 / V_0 stay scalar — they are built once
         // from blk.{f_0,V_0}_default (the RateStateAgingFrictionIterator
         // ctor) and V_0 is pinned to FrictionSolver::V0 by the R-009 guard, so
         // a per-rule f_0 / V_0 would be silently dropped: reject it loudly.
         // Per-DOF a / b / Dc / V_init / sigma_n / eta are allowed.
         MFEM_VERIFY(std::isnan(r.f_0) && std::isnan(r.V_0),
                     "ResolveRateState: per-DOF f_0/V_0 spatial overrides are "
                     "not supported (scalar aging-law globals; V_0 pinned by "
                     "R-009); set them only in the [friction.rate_state] "
                     "defaults. Offending rule at DOF " << i);
         if (!std::isnan(r.a))       { a_i      = r.a; }
         if (!std::isnan(r.b))       { b_i      = r.b; }  // Phase 11a: per-DOF b
         if (!std::isnan(r.Dc))      { Dc_i     = r.Dc; }
         if (!std::isnan(r.V_init))  { V_init_i = r.V_init; }
         if (!std::isnan(r.sigma_n)) { sn_i     = r.sigma_n; }
         if (!std::isnan(r.eta))     { eta_i = r.eta;  eta_rule_set = true; }
         if (!std::isnan(r.V_w))     { V_w_i    = r.V_w; }  // Phase 6 req 4 (SRW)
      }

      // Effective normal stress (depth convention: z<0 below surface).
      real_t sigma_n_total;
      if (sigma_n_total_per_dof.Size() == N)
      {
         sigma_n_total = sigma_n_total_per_dof(i);
      }
      else
      {
         sigma_n_total = sn_i;
      }
      const real_t P_p   = pp.P_p_pa + pp.P_p_grad_pa_per_m * depth_i;
      real_t sigma_n_eff = sigma_n_total - P_p;
      if (pp.min_sigma_n_pa > 0.0 && sigma_n_eff < pp.min_sigma_n_pa)
      {
         sigma_n_eff = pp.min_sigma_n_pa;
      }

      // eta: auto path computes 0.5 * sqrt(mu * rho) per DOF.
      if (cfg.eta_auto)
      {
         MFEM_VERIFY(!eta_rule_set,
                     "ResolveRateState: eta='auto' and a spatial rule "
                     "also set 'eta' at DOF " << i
                     << " — one source of truth required");
         const int e = dof_to_elem[i];
         ElementTransformation* T = mesh.GetElementTransformation(e);
         // True reference centroid for the bulk element.  For
         // MaterialField::Mode::Coefficient, ip.Init(0) lands on the
         // reference origin (a CORNER), not the centroid, which biases
         // mu and rho by the heterogeneity scale of one element.  The
         // centroid is the best per-element stop-gap until Phase 5
         // wires FaultGeometry::fault_dof_ip(i) (R-111).
         const Geometry::Type gtype = mesh.GetElementBaseGeometry(e);
         const IntegrationPoint& ip = Geometries.GetCenter(gtype);
         real_t lam_e, mu_e, rho_e;
         material.EvalAt(e, *T, ip, lam_e, mu_e, rho_e);
         MFEM_VERIFY(mu_e > 0.0 && rho_e > 0.0,
                     "ResolveRateState: eta=auto needs mu>0 and rho>0 at "
                     "DOF " << i << " (got mu=" << mu_e << ", rho="
                     << rho_e << ")");
         eta_i = 0.5 * std::sqrt(mu_e * rho_e);
      }

      // Per-DOF validator.
      MFEM_VERIFY(a_i > 0.0, "ResolveRateState: a <= 0 at DOF " << i);
      MFEM_VERIFY(b_i > 0.0,
                  "ResolveRateState: b <= 0 at DOF " << i << " (depth "
                  << depth_i << " m)"
                  << (cfg.depth_profile.enabled
                      ? "; check param_a_minus_b.csv: a - (a-b) went non-positive"
                      : ""));
      // R-011: a > b (velocity-strengthening) is allowed — it is how aging-law
      // ruptures arrest at the fault edges (TPV102 uses a high-a border).  Only
      // require a, b finite and > 0; do NOT require a < b.
      MFEM_VERIFY(std::isfinite(a_i) && std::isfinite(b_i),
                  "ResolveRateState: a, b must be finite at DOF " << i);
      MFEM_VERIFY(Dc_i > 0.0, "ResolveRateState: Dc <= 0 at DOF " << i);
      MFEM_VERIFY(V0_i > 0.0, "ResolveRateState: V_0 <= 0 at DOF " << i);
      MFEM_VERIFY(sigma_n_eff > 0.0,
                  "ResolveRateState: sigma_n_eff <= 0 at DOF " << i);
      MFEM_VERIFY(f0_i > 0.0 && f0_i < 1.0,
                  "ResolveRateState: f_0 not in (0,1) at DOF " << i);
      // Phase 6 req 4: V_w is consumed only by the slip-law-SRW state update;
      // require it > 0 only on that path (the aging path fills but ignores it).
      if (cfg.state_evolution == StateEvolutionKind::SlipLawStrongRateWeakening)
      {
         MFEM_VERIFY(V_w_i > 0.0,
                     "ResolveRateState: V_w <= 0 at DOF " << i
                     << " (required for state_evolution="
                     "slip_law_strong_rate_weakening)");
      }

      p.a(i) = a_i;   p.b(i) = b_i;   p.Dc(i) = Dc_i;
      p.V_init(i) = V_init_i;
      p.f_0(i) = f0_i; p.V_0(i) = V0_i;
      p.eta(i) = eta_i;
      p.sigma_n_eff(i) = sigma_n_eff;
      p.V_w(i) = V_w_i;
   }

   return p;
}

}  // namespace

RateStatePerDOFParams SpatialFrictionResolver::ResolveRateState(
   const RateStateBlock&     cfg,
   const Vector&             dof_coords_3d,
   const Array<int>&         dof_to_elem,
   const Array<int>&         dof_to_attr,
   const MaterialField&      material,
   ParMesh&                  pmesh,
   const PorePressureSpec&   pp,
   const Vector&             sigma_n_total_per_dof) const
{
   return resolve_rs_impl<ParMesh>(cfg, dof_coords_3d, dof_to_elem,
                                   dof_to_attr, material, pmesh, pp,
                                   sigma_n_total_per_dof);
}

RateStatePerDOFParams SpatialFrictionResolver::ResolveRateState(
   const RateStateBlock&     cfg,
   const Vector&             dof_coords_3d,
   const Array<int>&         dof_to_elem,
   const Array<int>&         dof_to_attr,
   const MaterialField&      material,
   mfem::Mesh&               mesh,
   const PorePressureSpec&   pp,
   const Vector&             sigma_n_total_per_dof) const
{
   return resolve_rs_impl<mfem::Mesh>(cfg, dof_coords_3d, dof_to_elem,
                                      dof_to_attr, material, mesh, pp,
                                      sigma_n_total_per_dof);
}

}  // namespace spatial
}  // namespace seas
}  // namespace mfem

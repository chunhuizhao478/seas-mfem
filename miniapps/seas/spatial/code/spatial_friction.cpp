// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// spatial_friction.cpp — Phase 1 implementation of
// spatial_dynamic_rupture_plan.md rev-3.
//
// All TOML I/O is gated by SEAS_USE_TOML.  Without toml11, the parser
// entry points abort at runtime; the resolver and time-helper paths
// compile unconditionally.

#include "spatial_friction.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <fstream>
#include <regex>
#include <sstream>
#include <string>

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

// SCEC C∞ Boxcar B(s; W, w) — Eq. (5) of TPV101/102/104.  `signed_offset`
// is `(coord − center)`; the function takes `|signed_offset|`.
real_t SCECBoxcar(real_t signed_offset, real_t W, real_t w)
{
   MFEM_VERIFY(W >= 0.0,
               "SCECBoxcar: half_width must be >= 0; got " << W);
   MFEM_VERIFY(w >= 0.0,
               "SCECBoxcar: transition must be >= 0; got " << w);
   const real_t s = std::abs(signed_offset);
   if (s <= W)         { return 1.0; }
   if (s >= W + w)     { return 0.0; }
   // s in (W, W+w): smooth tanh transition.  The two divisor terms
   // (s − W − w) and (s − W) are both finite and non-zero throughout
   // the open interval.
   return 0.5 * (1.0 + std::tanh(w / (s - W - w) + w / (s - W)));
}

real_t SpatialRule::BoxcarTaperFactor(real_t x, real_t y, real_t z) const
{
   real_t B = 1.0;
   if (!std::isnan(boxcar_half_x_m))
   {
      const real_t w = std::isnan(boxcar_trans_x_m)
                       ? static_cast<real_t>(0.0) : boxcar_trans_x_m;
      B *= SCECBoxcar(x - boxcar_center_x_m, boxcar_half_x_m, w);
   }
   if (!std::isnan(boxcar_half_y_m))
   {
      const real_t w = std::isnan(boxcar_trans_y_m)
                       ? static_cast<real_t>(0.0) : boxcar_trans_y_m;
      B *= SCECBoxcar(y - boxcar_center_y_m, boxcar_half_y_m, w);
   }
   if (!std::isnan(boxcar_half_z_m))
   {
      const real_t w = std::isnan(boxcar_trans_z_m)
                       ? static_cast<real_t>(0.0) : boxcar_trans_z_m;
      B *= SCECBoxcar(z - boxcar_center_z_m, boxcar_half_z_m, w);
   }
   return B;
}

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
      // BoxcarTaper always matches; the per-DOF blend weight is
      // computed by BoxcarTaperFactor in the resolver.  Returning true
      // unconditionally avoids re-evaluating the boxcar in the
      // resolver's match-test path.
      return true;
   }
   return false;
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
   if (s == "constant_tensor")    { return StressSourceKind::ConstantTensor; }
   if (s == "sidecar_hdf5")       { return StressSourceKind::SidecarHDF5; }
   if (s == "depth_proportional") { return StressSourceKind::DepthProportionalToShearModulus; }
   if (s == "constant_tensor_with_patches")
   { return StressSourceKind::ConstantTensorWithPatches; }
   MFEM_ABORT("stress.kind must be one of {constant_tensor, sidecar_hdf5, "
              "depth_proportional, constant_tensor_with_patches}; got '"
              << s << "'");
   return StressSourceKind::ConstantTensor;
}

MaterialKind parse_material_kind(const std::string& s)
{
   if (s == "constant")          { return MaterialKind::Constant; }
   if (s == "depth_profile_1d")  { return MaterialKind::DepthProfile1D; }
   if (s == "sidecar_hdf5")      { return MaterialKind::SidecarHDF5; }
   MFEM_ABORT("material.kind must be one of "
              "{constant, depth_profile_1d, sidecar_hdf5}; got '"
              << s << "'");
   return MaterialKind::Constant;
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
                 << "'.  Valid kinds (rev-3): depth, box, region_attribute, "
                 << "barrier, boxcar_taper.  (rev-1 'nucleation_box' removed "
                 << "by D-4; use the [nucleation] block instead.)");
   }

   // BoxcarTaper geometry fields (validated only when this rule's kind
   // is BoxcarTaper).  Each axis is OPTIONAL — an unset `boxcar_half_*`
   // means the axis is unconstrained (B = 1 along that axis).
   if (out.kind == SpatialRule::Kind::BoxcarTaper)
   {
      const real_t nan = std::numeric_limits<real_t>::quiet_NaN();
      out.boxcar_center_x_m = toml_real(rule_tbl, "boxcar_center_x_m", 0.0);
      out.boxcar_center_y_m = toml_real(rule_tbl, "boxcar_center_y_m", 0.0);
      out.boxcar_center_z_m = toml_real(rule_tbl, "boxcar_center_z_m", 0.0);
      out.boxcar_half_x_m   = toml_real(rule_tbl, "boxcar_half_x_m",   nan);
      out.boxcar_half_y_m   = toml_real(rule_tbl, "boxcar_half_y_m",   nan);
      out.boxcar_half_z_m   = toml_real(rule_tbl, "boxcar_half_z_m",   nan);
      out.boxcar_trans_x_m  = toml_real(rule_tbl, "boxcar_trans_x_m",  nan);
      out.boxcar_trans_y_m  = toml_real(rule_tbl, "boxcar_trans_y_m",  nan);
      out.boxcar_trans_z_m  = toml_real(rule_tbl, "boxcar_trans_z_m",  nan);

      const bool has_x = !std::isnan(out.boxcar_half_x_m);
      const bool has_y = !std::isnan(out.boxcar_half_y_m);
      const bool has_z = !std::isnan(out.boxcar_half_z_m);
      MFEM_VERIFY(has_x || has_y || has_z,
                  "[[spatial]] kind='boxcar_taper' requires at least one "
                  "of boxcar_half_{x,y,z}_m to be set");
      if (has_x)
      {
         MFEM_VERIFY(out.boxcar_half_x_m >= 0.0,
                     "boxcar_half_x_m must be >= 0");
         MFEM_VERIFY(!std::isnan(out.boxcar_trans_x_m)
                     && out.boxcar_trans_x_m >= 0.0,
                     "boxcar_trans_x_m required and >= 0 when "
                     "boxcar_half_x_m is set");
      }
      if (has_y)
      {
         MFEM_VERIFY(out.boxcar_half_y_m >= 0.0,
                     "boxcar_half_y_m must be >= 0");
         MFEM_VERIFY(!std::isnan(out.boxcar_trans_y_m)
                     && out.boxcar_trans_y_m >= 0.0,
                     "boxcar_trans_y_m required and >= 0 when "
                     "boxcar_half_y_m is set");
      }
      if (has_z)
      {
         MFEM_VERIFY(out.boxcar_half_z_m >= 0.0,
                     "boxcar_half_z_m must be >= 0");
         MFEM_VERIFY(!std::isnan(out.boxcar_trans_z_m)
                     && out.boxcar_trans_z_m >= 0.0,
                     "boxcar_trans_z_m required and >= 0 when "
                     "boxcar_half_z_m is set");
      }
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

      // Depth-linear cohesion taper (TPV31-style).  All four fields
      // are optional; when `cohesion_grad_pa_per_m` is set the resolver
      // computes per-DOF cohesion via the linear ramp and the constant
      // `cohesion` field is ignored.
      out.cohesion_grad_pa_per_m = toml_real(rule_tbl,
                                             "cohesion_grad_pa_per_m", nan);
      out.cohesion_ref_depth_m   = toml_real(rule_tbl,
                                             "cohesion_ref_depth_m", nan);
      out.cohesion_floor_pa      = toml_real(rule_tbl,
                                             "cohesion_floor_pa", 0.0);
      const std::string axis_s = toml_str(rule_tbl, "cohesion_taper_axis",
                                          std::string("y"));
      MFEM_VERIFY(axis_s == "x" || axis_s == "y" || axis_s == "z",
                  "[[spatial]] cohesion_taper_axis must be 'x', 'y', or 'z'; "
                  "got '" << axis_s << "'");
      out.cohesion_taper_axis = axis_s[0];
      // Either both grad + ref_depth are set, or neither.
      const bool has_grad  = !std::isnan(out.cohesion_grad_pa_per_m);
      const bool has_ref   = !std::isnan(out.cohesion_ref_depth_m);
      MFEM_VERIFY(has_grad == has_ref,
                  "[[spatial]] cohesion_grad_pa_per_m and "
                  "cohesion_ref_depth_m must both be set or both omitted");
      MFEM_VERIFY(out.cohesion_floor_pa >= 0.0,
                  "[[spatial]] cohesion_floor_pa must be >= 0; got "
                  << out.cohesion_floor_pa);

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
      out.V_w     = toml_real(rule_tbl, "V_w",     nan);
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
   // State-evolution kind (default = AgingLaw, matches BP5 / TPV102).
   if (rs_tbl.contains("state_evolution"))
   {
      const std::string se = rs_tbl.at("state_evolution").as_string();
      if      (se == "aging_law")
      { out.state_evolution = StateEvolutionKind::AgingLaw; }
      else if (se == "slip_law_srw")
      { out.state_evolution = StateEvolutionKind::SlipLawStrongRateWeakening; }
      else
      {
         MFEM_ABORT("[friction.rate_state].state_evolution must be "
                    "\"aging_law\" or \"slip_law_srw\"; got '" << se << "'");
      }
   }

   out.f_0_default     = toml_real(rs_tbl, "f_0_default",     0.6);
   out.V_0_default     = toml_real(rs_tbl, "V_0_default",     1.0e-6);
   out.a_default       = toml_real(rs_tbl, "a_default",       0.010);
   out.b_default       = toml_real(rs_tbl, "b_default",       0.015);
   out.Dc_default      = toml_real(rs_tbl, "Dc_default",      0.004);
   out.V_init_default  = toml_real(rs_tbl, "V_init_default",  1.0e-9);
   out.sigma_n_default = toml_real(rs_tbl, "sigma_n_default", 50.0e6);
   out.f_w_default     = toml_real(rs_tbl, "f_w_default",     0.2);
   out.V_w_default     = toml_real(rs_tbl, "V_w_default",     1.0);

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

   // Defaults validator.
   MFEM_VERIFY(out.a_default > 0.0,
               "[friction.rate_state] a_default must be > 0; got "
               << out.a_default);
   MFEM_VERIFY(out.b_default > 0.0,
               "[friction.rate_state] b_default must be > 0; got "
               << out.b_default);
   // Note: `a < b` (velocity-weakening) is intentionally NOT required
   // at the default level — TPV102/104 use `a_vs > b > a_vw`, with the
   // outside-VW (VS) region stable and the inside-VW (VW) core
   // unstable.  The parser only requires positivity here; per-DOF
   // physics regime is determined by the spatial rules + defaults.
   MFEM_VERIFY(out.Dc_default > 0.0,
               "[friction.rate_state] Dc_default must be > 0");
   MFEM_VERIFY(out.V_0_default > 0.0,
               "[friction.rate_state] V_0_default must be > 0");
   MFEM_VERIFY(out.sigma_n_default > 0.0,
               "[friction.rate_state] sigma_n_default must be > 0");
   MFEM_VERIFY(out.f_0_default > 0.0 && out.f_0_default < 1.0,
               "[friction.rate_state] f_0_default must be in (0, 1)");
   if (out.state_evolution == StateEvolutionKind::SlipLawStrongRateWeakening)
   {
      MFEM_VERIFY(out.f_w_default > 0.0 && out.f_w_default < 1.0,
                  "[friction.rate_state] f_w_default ("
                  << out.f_w_default << ") must be in (0, 1) when "
                  "state_evolution=\"slip_law_srw\"");
      MFEM_VERIFY(out.V_w_default > 0.0,
                  "[friction.rate_state] V_w_default ("
                  << out.V_w_default << ") must be > 0 when "
                  "state_evolution=\"slip_law_srw\"");
   }
}

// ---------------------------------------------------------------------------
// Phase R.3 — helpers for the new top-level sections.
// ---------------------------------------------------------------------------

std::vector<int> toml_int_array(const toml::value& tbl, const std::string& key,
                                const std::vector<int>& default_val)
{
   if (!tbl.contains(key)) { return default_val; }
   const auto& arr = tbl.at(key).as_array();
   std::vector<int> out;
   out.reserve(arr.size());
   for (const auto& v : arr)
   {
      MFEM_VERIFY(v.is_integer(),
                  "TOML array key '" << key << "' must contain integers");
      out.push_back(static_cast<int>(v.as_integer()));
   }
   return out;
}

std::array<real_t, 3> toml_real3(const toml::value& tbl,
                                 const std::string& key,
                                 const std::array<real_t, 3>& default_val)
{
   if (!tbl.contains(key)) { return default_val; }
   const auto& arr = tbl.at(key).as_array();
   MFEM_VERIFY(arr.size() == 3,
               "TOML key '" << key << "' must be a length-3 array; got "
               << arr.size());
   std::array<real_t, 3> out{};
   for (int i = 0; i < 3; ++i)
   {
      MFEM_VERIFY(arr[i].is_floating() || arr[i].is_integer(),
                  "TOML key '" << key << "' element " << i
                  << " must be numeric");
      out[i] = arr[i].is_floating()
               ? static_cast<real_t>(arr[i].as_floating())
               : static_cast<real_t>(arr[i].as_integer());
   }
   return out;
}

void parse_problem(const toml::value& tbl, ProblemSpec& out)
{
   out.tag = toml_str(tbl, "tag", out.tag);
}

void parse_boundary(const toml::value& tbl, BoundarySpec& out)
{
   out.fault_attr      = toml_int      (tbl, "fault_attr",      out.fault_attr);
   out.natural_attrs   = toml_int_array(tbl, "natural_attrs",   out.natural_attrs);
   out.absorbing_attrs = toml_int_array(tbl, "absorbing_attrs", out.absorbing_attrs);

   MFEM_VERIFY(out.fault_attr > 0,
               "[boundary].fault_attr must be > 0; got " << out.fault_attr);
   // Mutex check: no attribute appears in more than one of
   // {fault, natural, absorbing}.
   auto in_list = [](int a, const std::vector<int>& v) {
      for (int x : v) { if (x == a) { return true; } }
      return false;
   };
   MFEM_VERIFY(!in_list(out.fault_attr, out.natural_attrs)
               && !in_list(out.fault_attr, out.absorbing_attrs),
               "[boundary].fault_attr (" << out.fault_attr
               << ") appears in natural_attrs or absorbing_attrs");
   for (int a : out.natural_attrs)
   {
      MFEM_VERIFY(!in_list(a, out.absorbing_attrs),
                  "[boundary] attr " << a
                  << " appears in both natural_attrs and absorbing_attrs");
   }
}

void parse_fault_geometry(const toml::value& tbl, FaultGeometrySpec& out)
{
   out.ref_normal = toml_real3(tbl, "ref_normal", out.ref_normal);
   out.up         = toml_real3(tbl, "up",         out.up);
   out.kind       = toml_str  (tbl, "kind",       out.kind);

   // Unit-norm check on ref_normal.
   const real_t n2 = out.ref_normal[0]*out.ref_normal[0]
                   + out.ref_normal[1]*out.ref_normal[1]
                   + out.ref_normal[2]*out.ref_normal[2];
   MFEM_VERIFY(std::abs(n2 - 1.0) < 1e-10,
               "[fault_geometry].ref_normal must be a unit vector; got "
               "|n|^2=" << n2);

   // up must NOT be parallel to ref_normal (otherwise FaultBasis::Compute
   // would have an ambiguous tangent frame).
   const real_t dot = out.ref_normal[0]*out.up[0]
                    + out.ref_normal[1]*out.up[1]
                    + out.ref_normal[2]*out.up[2];
   const real_t up2 = out.up[0]*out.up[0] + out.up[1]*out.up[1]
                    + out.up[2]*out.up[2];
   MFEM_VERIFY(up2 > 0.0,
               "[fault_geometry].up must be non-zero");
   const real_t cos2 = dot * dot / up2;   // cos^2(angle), n is unit
   MFEM_VERIFY(cos2 < 1.0 - 1e-10,
               "[fault_geometry].up is parallel (or anti-parallel) to "
               "ref_normal — FaultBasis cannot form a tangent frame");

   // [fault_geometry].kind is INFORMATIONAL only — no code branches on
   // this string.  Any non-empty value is accepted so users can tag
   // configs with the benchmark name (e.g., "tpv102_rs", "tpv104_srw")
   // without modifying the parser.
   MFEM_VERIFY(!out.kind.empty(),
               "[fault_geometry].kind must be a non-empty string");
}

void parse_hypocenter(const toml::value& tbl, HypocenterSpec& out)
{
   out.x                   = toml_real(tbl, "x",                   out.x);
   out.y                   = toml_real(tbl, "y",                   out.y);
   out.z                   = toml_real(tbl, "z",                   out.z);
   out.nucleation_radius_m = toml_real(tbl, "nucleation_radius_m", out.nucleation_radius_m);
   out.nucleation_taper_m  = toml_real(tbl, "nucleation_taper_m",  out.nucleation_taper_m);

   MFEM_VERIFY(out.nucleation_radius_m > 0.0,
               "[hypocenter].nucleation_radius_m must be > 0; got "
               << out.nucleation_radius_m);
   MFEM_VERIFY(out.nucleation_taper_m >= 0.0,
               "[hypocenter].nucleation_taper_m must be >= 0; got "
               << out.nucleation_taper_m);
}

void parse_material(const toml::value& tbl, MaterialSpec& out)
{
   const std::string kind_s = toml_str(tbl, "kind", "constant");
   out.kind = parse_material_kind(kind_s);

   if (out.kind == MaterialKind::DepthProfile1D)
   {
      MFEM_VERIFY(tbl.contains("depth_axis") == false
                  || toml_str(tbl, "depth_axis", "y").size() == 1,
                  "[material].depth_axis must be a single character "
                  "from {'x','y','z'}");
      const std::string axis_s = toml_str(tbl, "depth_axis", "y");
      MFEM_VERIFY(axis_s.size() == 1 &&
                  (axis_s[0] == 'x' || axis_s[0] == 'y' || axis_s[0] == 'z'),
                  "[material].depth_axis must be 'x', 'y', or 'z'; got '"
                  << axis_s << "'");
      out.depth_axis = axis_s[0];
   }
}

void parse_material_profile(const toml::value& root, MaterialSpec& out)
{
   if (out.kind != MaterialKind::DepthProfile1D) { return; }

   MFEM_VERIFY(root.contains("material_profile"),
               "[material].kind=\"depth_profile_1d\" requires a "
               "[[material_profile.layer]] array of layer tables");
   const auto& mp = root.at("material_profile");
   MFEM_VERIFY(mp.contains("layer"),
               "[material_profile] must contain a `layer` array "
               "(use TOML syntax `[[material_profile.layer]]`)");
   const auto& layers_arr = mp.at("layer").as_array();
   MFEM_VERIFY(!layers_arr.empty(),
               "[material_profile.layer] must have at least one entry");

   out.profile_layers.clear();
   out.profile_layers.reserve(layers_arr.size());
   for (std::size_t i = 0; i < layers_arr.size(); ++i)
   {
      const auto& L = layers_arr[i];
      DepthProfileLayer dpl;
      dpl.depth_top_m  = toml_real(L, "depth_top_m", 0.0);
      dpl.depth_bot_m  = toml_real(L, "depth_bot_m", 0.0);
      dpl.vp_ms        = toml_real(L, "vp_ms",       0.0);
      dpl.vs_ms        = toml_real(L, "vs_ms",       0.0);
      dpl.rho_kgm3     = toml_real(L, "rho_kgm3",    0.0);
      dpl.interp       = toml_str (L, "interp",      "constant");
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
      out.profile_layers.push_back(dpl);
   }
}

void parse_depth_proportional_stress(const toml::value& tbl,
                                     DepthProportionalStressSpec& out)
{
   // Component values are stored in MPa in the TOML (per plan §R.5
   // step 4: "Initial stress block... carrying the 6 scaled components
   // in MPa + mu_ref_pa = 32.03812032e9").  Multiply by 1e6 internally.
   const real_t mpa = 1.0e6;
   out.sigma_xx_per_mu = mpa * toml_real(tbl, "sigma_xx_per_mu", 0.0);
   out.sigma_yy_per_mu = mpa * toml_real(tbl, "sigma_yy_per_mu", 0.0);
   out.sigma_zz_per_mu = mpa * toml_real(tbl, "sigma_zz_per_mu", 0.0);
   out.sigma_xy_per_mu = mpa * toml_real(tbl, "sigma_xy_per_mu", 0.0);
   out.sigma_yz_per_mu = mpa * toml_real(tbl, "sigma_yz_per_mu", 0.0);
   out.sigma_xz_per_mu = mpa * toml_real(tbl, "sigma_xz_per_mu", 0.0);
   out.mu_ref_pa       = toml_real(tbl, "mu_ref_pa", 32.03812032e9);
   MFEM_VERIFY(out.mu_ref_pa > 0.0,
               "[stress.depth_proportional].mu_ref_pa must be > 0; got "
               << out.mu_ref_pa);
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
   // Phase R.3: [velocity] is now OPTIONAL — only required when
   // [material].kind = "sidecar_hdf5" (or when material is absent
   // and the legacy sidecar-velocity workflow is in effect).  See
   // material_kind dispatch below.
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
   }
   // Phase R.3: [velocity] requirement deferred until after the material
   // section is parsed — only required when material.kind == sidecar_hdf5.

   if (root.contains("stress"))
   {
      const auto& s = root.at("stress");
      MFEM_VERIFY(s.contains("kind"),
                  "[stress].kind is required; must be one of "
                  "\"constant_tensor\", \"sidecar_hdf5\", "
                  "\"depth_proportional\"");
      cfg.stress.kind = parse_stress_kind(toml_str(s, "kind", "constant_tensor"));

      const bool has_sxx = s.contains("sigma_xx_pa");
      const bool has_syy = s.contains("sigma_yy_pa");
      const bool has_szz = s.contains("sigma_zz_pa");
      const bool has_sxy = s.contains("sigma_xy_pa");
      const bool has_syz = s.contains("sigma_yz_pa");
      const bool has_sxz = s.contains("sigma_xz_pa");
      const bool has_path = s.contains("sidecar_path")
                            && !toml_str(s, "sidecar_path", std::string()).empty();
      const bool has_dp_block = s.contains("depth_proportional");

      const bool has_patches = s.contains("patch");
      if (cfg.stress.kind == StressSourceKind::ConstantTensor)
      {
         MFEM_VERIFY(has_sxx && has_syy && has_szz
                     && has_sxy && has_syz && has_sxz,
                     "[stress] kind=\"constant_tensor\" requires all six "
                     "sigma_*_pa keys to be present");
         MFEM_VERIFY(!has_path,
                     "[stress] kind=\"constant_tensor\" must NOT set "
                     "sidecar_path (it is for kind=\"sidecar_hdf5\" only)");
         MFEM_VERIFY(!has_dp_block,
                     "[stress] kind=\"constant_tensor\" must NOT set the "
                     "[stress.depth_proportional] sub-block (it is for "
                     "kind=\"depth_proportional\" only)");
         MFEM_VERIFY(!has_patches,
                     "[stress] kind=\"constant_tensor\" must NOT set any "
                     "[[stress.patch]] entries (use "
                     "kind=\"constant_tensor_with_patches\" instead)");
         cfg.stress.sigma_xx_pa = toml_real(s, "sigma_xx_pa", 0.0);
         cfg.stress.sigma_yy_pa = toml_real(s, "sigma_yy_pa", 0.0);
         cfg.stress.sigma_zz_pa = toml_real(s, "sigma_zz_pa", 0.0);
         cfg.stress.sigma_xy_pa = toml_real(s, "sigma_xy_pa", 0.0);
         cfg.stress.sigma_yz_pa = toml_real(s, "sigma_yz_pa", 0.0);
         cfg.stress.sigma_xz_pa = toml_real(s, "sigma_xz_pa", 0.0);
      }
      else if (cfg.stress.kind == StressSourceKind::ConstantTensorWithPatches)
      {
         // Same six-key requirement on the BACKGROUND tensor as
         // constant_tensor; patches override individual components on
         // top of it.
         MFEM_VERIFY(has_sxx && has_syy && has_szz
                     && has_sxy && has_syz && has_sxz,
                     "[stress] kind=\"constant_tensor_with_patches\" "
                     "requires all six background sigma_*_pa keys to "
                     "be present");
         MFEM_VERIFY(!has_path,
                     "[stress] kind=\"constant_tensor_with_patches\" "
                     "must NOT set sidecar_path");
         MFEM_VERIFY(!has_dp_block,
                     "[stress] kind=\"constant_tensor_with_patches\" "
                     "must NOT set the [stress.depth_proportional] "
                     "sub-block");
         MFEM_VERIFY(has_patches,
                     "[stress] kind=\"constant_tensor_with_patches\" "
                     "requires at least one [[stress.patch]] entry "
                     "(otherwise use kind=\"constant_tensor\")");

         cfg.stress.sigma_xx_pa = toml_real(s, "sigma_xx_pa", 0.0);
         cfg.stress.sigma_yy_pa = toml_real(s, "sigma_yy_pa", 0.0);
         cfg.stress.sigma_zz_pa = toml_real(s, "sigma_zz_pa", 0.0);
         cfg.stress.sigma_xy_pa = toml_real(s, "sigma_xy_pa", 0.0);
         cfg.stress.sigma_yz_pa = toml_real(s, "sigma_yz_pa", 0.0);
         cfg.stress.sigma_xz_pa = toml_real(s, "sigma_xz_pa", 0.0);

         // Parse the [[stress.patch]] array.  toml11 represents an array
         // of inline tables as a value with is_array()==true; each entry
         // is itself a table.
         const auto& parr = s.at("patch");
         MFEM_VERIFY(parr.is_array(),
                     "[stress.patch] must be a TOML array of tables");
         const auto& vec = parr.as_array();
         MFEM_VERIFY(!vec.empty(),
                     "[stress.patch] must contain at least one patch");
         const real_t nan = std::numeric_limits<real_t>::quiet_NaN();
         const real_t pos_inf =  std::numeric_limits<real_t>::infinity();
         cfg.stress.patches.clear();
         cfg.stress.patches.reserve(vec.size());
         for (std::size_t pi = 0; pi < vec.size(); ++pi)
         {
            const auto& pt = vec[pi];
            MFEM_VERIFY(pt.is_table(),
                        "[[stress.patch]] entry " << pi << " must be a table");
            StressPatch sp;
            sp.center_x_m  = toml_real(pt, "center_x_m", nan);
            sp.center_y_m  = toml_real(pt, "center_y_m", nan);
            sp.center_z_m  = toml_real(pt, "center_z_m", nan);
            // Missing half-* ⇒ +inf (no constraint along that axis).
            sp.half_x_m    = toml_real(pt, "half_x_m", pos_inf);
            sp.half_y_m    = toml_real(pt, "half_y_m", pos_inf);
            sp.half_z_m    = toml_real(pt, "half_z_m", pos_inf);
            sp.sigma_xx_pa = toml_real(pt, "sigma_xx_pa", nan);
            sp.sigma_yy_pa = toml_real(pt, "sigma_yy_pa", nan);
            sp.sigma_zz_pa = toml_real(pt, "sigma_zz_pa", nan);
            sp.sigma_xy_pa = toml_real(pt, "sigma_xy_pa", nan);
            sp.sigma_yz_pa = toml_real(pt, "sigma_yz_pa", nan);
            sp.sigma_xz_pa = toml_real(pt, "sigma_xz_pa", nan);

            // Required: center_* on every axis that has a finite half-*.
            // If half_x_m is +inf the patch spans all x and center_x_m is
            // unused (NaN is fine); same for y / z.
            const bool has_xc = !std::isnan(sp.center_x_m);
            const bool has_yc = !std::isnan(sp.center_y_m);
            const bool has_zc = !std::isnan(sp.center_z_m);
            const bool x_constrained = std::isfinite(sp.half_x_m);
            const bool y_constrained = std::isfinite(sp.half_y_m);
            const bool z_constrained = std::isfinite(sp.half_z_m);
            MFEM_VERIFY(!x_constrained || has_xc,
                        "[[stress.patch]] entry " << pi
                        << " sets half_x_m but not center_x_m");
            MFEM_VERIFY(!y_constrained || has_yc,
                        "[[stress.patch]] entry " << pi
                        << " sets half_y_m but not center_y_m");
            MFEM_VERIFY(!z_constrained || has_zc,
                        "[[stress.patch]] entry " << pi
                        << " sets half_z_m but not center_z_m");
            MFEM_VERIFY(x_constrained || y_constrained || z_constrained,
                        "[[stress.patch]] entry " << pi
                        << " has all three half_*_m = +inf — the patch "
                        "would cover the whole domain; supply at least "
                        "one finite half_*_m");
            MFEM_VERIFY(!x_constrained || sp.half_x_m >= 0.0,
                        "[[stress.patch]] entry " << pi
                        << " half_x_m must be >= 0; got " << sp.half_x_m);
            MFEM_VERIFY(!y_constrained || sp.half_y_m >= 0.0,
                        "[[stress.patch]] entry " << pi
                        << " half_y_m must be >= 0; got " << sp.half_y_m);
            MFEM_VERIFY(!z_constrained || sp.half_z_m >= 0.0,
                        "[[stress.patch]] entry " << pi
                        << " half_z_m must be >= 0; got " << sp.half_z_m);
            // Require at least one stress override per patch.
            const bool any_override =
               !std::isnan(sp.sigma_xx_pa) || !std::isnan(sp.sigma_yy_pa) ||
               !std::isnan(sp.sigma_zz_pa) || !std::isnan(sp.sigma_xy_pa) ||
               !std::isnan(sp.sigma_yz_pa) || !std::isnan(sp.sigma_xz_pa);
            MFEM_VERIFY(any_override,
                        "[[stress.patch]] entry " << pi
                        << " has no sigma_*_pa overrides — the patch "
                        "would be a no-op");
            cfg.stress.patches.push_back(sp);
         }
      }
      else if (cfg.stress.kind == StressSourceKind::DepthProportionalToShearModulus)
      {
         MFEM_VERIFY(!has_sxx && !has_syy && !has_szz
                     && !has_sxy && !has_syz && !has_sxz,
                     "[stress] kind=\"depth_proportional\" must NOT set "
                     "sigma_*_pa keys directly; put the 6 scaled "
                     "components inside the [stress.depth_proportional] "
                     "sub-block (in MPa)");
         MFEM_VERIFY(!has_path,
                     "[stress] kind=\"depth_proportional\" must NOT set "
                     "sidecar_path");
         MFEM_VERIFY(!has_patches,
                     "[stress] kind=\"depth_proportional\" must NOT set "
                     "[[stress.patch]] entries");
         MFEM_VERIFY(has_dp_block,
                     "[stress] kind=\"depth_proportional\" requires a "
                     "[stress.depth_proportional] sub-block with "
                     "sigma_xx_per_mu, ..., sigma_xz_per_mu (MPa) and "
                     "mu_ref_pa");
         parse_depth_proportional_stress(s.at("depth_proportional"),
                                         cfg.stress.depth_proportional);
      }
      else
      {
         MFEM_VERIFY(!has_sxx && !has_syy && !has_szz
                     && !has_sxy && !has_syz && !has_sxz,
                     "[stress] kind=\"sidecar_hdf5\" must NOT set any of "
                     "sigma_xx_pa..sigma_xz_pa (those are for "
                     "kind=\"constant_tensor\" only)");
         MFEM_VERIFY(!has_dp_block,
                     "[stress] kind=\"sidecar_hdf5\" must NOT set the "
                     "[stress.depth_proportional] sub-block");
         MFEM_VERIFY(!has_patches,
                     "[stress] kind=\"sidecar_hdf5\" must NOT set "
                     "[[stress.patch]] entries");
         MFEM_VERIFY(has_path,
                     "[stress] kind=\"sidecar_hdf5\" requires a non-empty "
                     "sidecar_path");
         cfg.stress.sidecar_path = toml_str(s, "sidecar_path", std::string());
      }
   }

   if (root.contains("numerics"))
   {
      const auto& n = root.at("numerics");
      cfg.numerics.ader_order     = toml_int (n, "ader_order",  2);
      cfg.numerics.mixed_flux     = toml_str (n, "mixed_flux",  "none");
      cfg.numerics.cfl            = toml_real(n, "cfl",         0.5);
      cfg.numerics.use_pml        = toml_bool(n, "use_pml",     false);
      // REVIEW R-004 / R-005 / R-006: TPV-parity opt-in knobs.
      cfg.numerics.cfl_safety     = toml_str (n, "cfl_safety",     "raw");
      cfg.numerics.fault_iterator = toml_str (n, "fault_iterator", "one-shot");
      cfg.numerics.interior_flux  = toml_str (n, "interior_flux",  "bimaterial");
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
   MFEM_VERIFY(cfg.numerics.cfl_safety == "raw"
               || cfg.numerics.cfl_safety == "dg",
               "[numerics].cfl_safety must be one of {raw, dg}; got '"
               << cfg.numerics.cfl_safety << "' (REVIEW R-004)");
   MFEM_VERIFY(cfg.numerics.fault_iterator == "one-shot"
               || cfg.numerics.fault_iterator == "substep",
               "[numerics].fault_iterator must be one of {one-shot, substep}; "
               "got '" << cfg.numerics.fault_iterator << "' (REVIEW R-005)");
   MFEM_VERIFY(cfg.numerics.interior_flux == "bimaterial"
               || cfg.numerics.interior_flux == "scalar",
               "[numerics].interior_flux must be one of {bimaterial, scalar}; "
               "got '" << cfg.numerics.interior_flux << "' (REVIEW R-006)");

   if (root.contains("time"))
   {
      const auto& t = root.at("time");
      cfg.time.tfinal     = toml_time_seconds(t, "tfinal",     12.0);
      // Drive-by: t_initial now accepts both numeric and string ("0s")
      // forms, matching the rest of the [time] block.  Previously a
      // string here aborted with "t_initial must be a number" even
      // though tpv205.toml writes it as "0s" for consistency.
      cfg.time.t_initial  = toml_time_seconds(t, "t_initial",  0.0);
      cfg.time.dt_initial = parse_dt_initial (t);
      cfg.time.dt_max     = toml_time_seconds(t, "dt_max",     0.1);
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
   for (const auto& mode_name : { std::pair<std::string, std::string>
                                  {"paraview_volume", cfg.output.paraview_volume},
                                  {"paraview_bulk",   cfg.output.paraview_bulk},
                                  {"paraview_fault",  cfg.output.paraview_fault} })
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
      MFEM_VERIFY(kind_s == "gradual_overstress"
                  || kind_s == "square_overstress"
                  || kind_s == "instantaneous_overstress_circular"
                  || kind_s == "gradual_overstress_compact_circular",
                  "[nucleation].kind must be \"gradual_overstress\", "
                  "\"square_overstress\", "
                  "\"instantaneous_overstress_circular\", or "
                  "\"gradual_overstress_compact_circular\"; got '"
                  << kind_s << "'");

      // Shared mutual-exclusion: only the sub-block matching `kind` may
      // appear.  We collect the four optional sub-block names and check
      // them once per branch below.
      auto reject_extra_subblocks =
         [&nuc](const std::string& allowed)
      {
         static const char* kAll[] = {
            "gradual_overstress", "square_overstress",
            "instantaneous_overstress_circular",
            "gradual_overstress_compact_circular"
         };
         for (const char* k : kAll)
         {
            if (allowed == k) { continue; }
            MFEM_VERIFY(!nuc.contains(k),
                        "[nucleation] kind=\"" << allowed
                        << "\" must NOT set the [nucleation." << k
                        << "] sub-block");
         }
      };

      if (kind_s == "gradual_overstress")
      {
         cfg.nucleation.kind = NucleationKind::GradualOverstress;
         MFEM_VERIFY(nuc.contains("gradual_overstress"),
                     "[nucleation] kind=\"gradual_overstress\" requires a "
                     "[nucleation.gradual_overstress] sub-block");
         reject_extra_subblocks("gradual_overstress");
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
      else if (kind_s == "square_overstress")
      {
         cfg.nucleation.kind = NucleationKind::SquareOverstress;
         MFEM_VERIFY(nuc.contains("square_overstress"),
                     "[nucleation] kind=\"square_overstress\" requires a "
                     "[nucleation.square_overstress] sub-block");
         reject_extra_subblocks("square_overstress");
         const auto& sq = nuc.at("square_overstress");
         auto& ss = cfg.nucleation.square_overstress;
         ss.T_nuc_s = toml_time_seconds(sq, "T_nuc_s", 0.0);
         MFEM_VERIFY(ss.T_nuc_s > 0.0,
                     "[nucleation.square_overstress].T_nuc_s must be > 0; "
                     "got " << ss.T_nuc_s);

         MFEM_VERIFY(sq.contains("patch"),
                     "[nucleation.square_overstress] requires at least one "
                     "[[nucleation.square_overstress.patch]] entry");
         const auto& parr = sq.at("patch");
         MFEM_VERIFY(parr.is_array(),
                     "[nucleation.square_overstress.patch] must be an array");
         const auto& vec = parr.as_array();
         MFEM_VERIFY(!vec.empty(),
                     "[nucleation.square_overstress.patch] must be non-empty");
         const real_t pos_inf = std::numeric_limits<real_t>::infinity();
         ss.patches.clear();
         ss.patches.reserve(vec.size());
         for (std::size_t pi = 0; pi < vec.size(); ++pi)
         {
            const auto& pt = vec[pi];
            MFEM_VERIFY(pt.is_table(),
                        "[[nucleation.square_overstress.patch]] entry "
                        << pi << " must be a table");
            SquareOverstressPatch sp;
            sp.center_x_m          = toml_real(pt, "center_x_m", 0.0);
            sp.center_y_m          = toml_real(pt, "center_y_m", 0.0);
            sp.center_z_m          = toml_real(pt, "center_z_m", 0.0);
            sp.half_x_m            = toml_real(pt, "half_x_m",   pos_inf);
            sp.half_y_m            = toml_real(pt, "half_y_m",   pos_inf);
            sp.half_z_m            = toml_real(pt, "half_z_m",   pos_inf);
            sp.delta_tau_dip_pa    = toml_real(pt, "delta_tau_dip_pa",    0.0);
            sp.delta_tau_strike_pa = toml_real(pt, "delta_tau_strike_pa", 0.0);
            const bool x_constrained = std::isfinite(sp.half_x_m);
            const bool y_constrained = std::isfinite(sp.half_y_m);
            const bool z_constrained = std::isfinite(sp.half_z_m);
            MFEM_VERIFY(x_constrained || y_constrained || z_constrained,
                        "[[nucleation.square_overstress.patch]] entry "
                        << pi << " has all three half_*_m = +inf");
            MFEM_VERIFY(!x_constrained || sp.half_x_m >= 0.0,
                        "[[nucleation.square_overstress.patch]] entry "
                        << pi << " half_x_m must be >= 0");
            MFEM_VERIFY(!y_constrained || sp.half_y_m >= 0.0,
                        "[[nucleation.square_overstress.patch]] entry "
                        << pi << " half_y_m must be >= 0");
            MFEM_VERIFY(!z_constrained || sp.half_z_m >= 0.0,
                        "[[nucleation.square_overstress.patch]] entry "
                        << pi << " half_z_m must be >= 0");
            ss.patches.push_back(sp);
         }
      }
      else if (kind_s == "instantaneous_overstress_circular")
      {
         cfg.nucleation.kind = NucleationKind::InstantaneousOverstressCircular;
         MFEM_VERIFY(nuc.contains("instantaneous_overstress_circular"),
                     "[nucleation] kind=\"instantaneous_overstress_circular\" "
                     "requires a [nucleation.instantaneous_overstress_circular] "
                     "sub-block");
         reject_extra_subblocks("instantaneous_overstress_circular");
         const auto& ic = nuc.at("instantaneous_overstress_circular");
         auto& is = cfg.nucleation.instantaneous_overstress_circular;
         is.center_x_m         = toml_real(ic, "center_x_m", 0.0);
         is.center_y_m         = toml_real(ic, "center_y_m", 0.0);
         is.center_z_m         = toml_real(ic, "center_z_m", 0.0);
         is.radius_inner_m     = toml_real(ic, "radius_inner_m", 0.0);
         is.radius_outer_m     = toml_real(ic, "radius_outer_m", 0.0);
         is.delta_tau_peak_pa  = toml_real(ic, "delta_tau_peak_pa", 0.0);
         is.mu_ref_pa          = toml_real(ic, "mu_ref_pa", 32.03812032e9);
         is.dip_fraction       = toml_real(ic, "dip_fraction",    0.0);
         is.strike_fraction    = toml_real(ic, "strike_fraction", 1.0);
         MFEM_VERIFY(is.radius_inner_m > 0.0,
                     "[nucleation.instantaneous_overstress_circular]"
                     ".radius_inner_m must be > 0; got "
                     << is.radius_inner_m);
         MFEM_VERIFY(is.radius_outer_m >= is.radius_inner_m,
                     "[nucleation.instantaneous_overstress_circular]"
                     ".radius_outer_m (" << is.radius_outer_m
                     << ") must be >= radius_inner_m ("
                     << is.radius_inner_m << ")");
         MFEM_VERIFY(is.mu_ref_pa > 0.0,
                     "[nucleation.instantaneous_overstress_circular]"
                     ".mu_ref_pa must be > 0; got " << is.mu_ref_pa);
      }
      else  // gradual_overstress_compact_circular
      {
         cfg.nucleation.kind =
            NucleationKind::GradualOverstressCompactCircular;
         MFEM_VERIFY(nuc.contains("gradual_overstress_compact_circular"),
                     "[nucleation] kind=\"gradual_overstress_compact_circular\" "
                     "requires a [nucleation.gradual_overstress_compact_circular] "
                     "sub-block");
         reject_extra_subblocks("gradual_overstress_compact_circular");
         const auto& g = nuc.at("gradual_overstress_compact_circular");
         auto& gs = cfg.nucleation.gradual_overstress_compact_circular;
         gs.center_x_m          = toml_real(g, "center_x_m",          0.0);
         gs.center_y_m          = toml_real(g, "center_y_m",          0.0);
         gs.center_z_m          = toml_real(g, "center_z_m",          0.0);
         gs.radius_m            = toml_real(g, "radius_m",            0.0);
         gs.delta_tau_dip_pa    = toml_real(g, "delta_tau_dip_pa",    0.0);
         gs.delta_tau_strike_pa = toml_real(g, "delta_tau_strike_pa", 0.0);
         gs.T_nuc_s             = toml_time_seconds(g, "T_nuc_s",     0.0);

         MFEM_VERIFY(gs.radius_m > 0.0,
                     "[nucleation.gradual_overstress_compact_circular]"
                     ".radius_m must be > 0; got " << gs.radius_m);
         MFEM_VERIFY(gs.T_nuc_s > 0.0,
                     "[nucleation.gradual_overstress_compact_circular]"
                     ".T_nuc_s must be > 0; got " << gs.T_nuc_s);
      }
   }
   // else: enabled stays false; driver runs without nucleation perturbation.

   // -----------------------------------------------------------------
   // Phase R.3 new top-level sections.  Every one is OPTIONAL; missing
   // sections fall back to the struct defaults, which match the SAFS
   // hardcoded values (backwards-compatibility per plan §R.3 step 11).
   // -----------------------------------------------------------------
   if (root.contains("problem"))
   {
      parse_problem(root.at("problem"), cfg.problem);
   }
   if (root.contains("boundary"))
   {
      parse_boundary(root.at("boundary"), cfg.boundary);
   }
   if (root.contains("fault_geometry"))
   {
      parse_fault_geometry(root.at("fault_geometry"), cfg.fault_geometry);
   }
   if (root.contains("hypocenter"))
   {
      parse_hypocenter(root.at("hypocenter"), cfg.hypocenter);
      // REVIEW R-008: catch the SPEC-depth-vs-mesh-z foot-gun.  Mesh z
      // is NEGATIVE below the free surface (CLAUDE.md "Canonical
      // Coordinate System"); hypocenter z should therefore be <= 0
      // when up has positive z (the canonical case).  A positive
      // hypocenter z with positive up_z would place the hypocenter
      // above the free surface — almost certainly a SPEC-depth value
      // copied without the sign flip.  The check is conditional on
      // up[2] > 0 so TPV31's rotated frame (up = (0, -1, 0)) is not
      // affected; TPV31 has up_z = 0 and hypocenter z = 0.
      if (cfg.fault_geometry.up[2] > 0.0)
      {
         MFEM_VERIFY(cfg.hypocenter.z <= 0.0,
                     "[hypocenter].z = " << cfg.hypocenter.z
                     << " is positive, but [fault_geometry].up_z > 0 "
                     "(canonical convention has z < 0 below the surface). "
                     "A positive hypocenter z places the hypocenter ABOVE "
                     "the free surface, almost certainly a SPEC-depth value "
                     "copied without the mesh-z sign flip.  Set z = -<depth_m> "
                     "in the TOML (e.g., z = -7500.0 for spec depth 7.5 km).  "
                     "See REVIEW R-008 and `CLAUDE.md` \"Canonical Coordinate "
                     "System\".");
      }
   }
   if (root.contains("material"))
   {
      parse_material(root.at("material"), cfg.material);
      parse_material_profile(root, cfg.material);
   }

   // Velocity requirement: now scoped to material.kind == sidecar_hdf5
   // (the legacy SAFS dataset_root path).  For Constant + DepthProfile1D,
   // the wave-operator material comes from the [material_constant_fallback]
   // or [material_profile.layer] blocks respectively, and the [velocity]
   // block is unused.
   if (cfg.material.kind == MaterialKind::SidecarHDF5)
   {
      MFEM_VERIFY(root.contains("velocity"),
                  "[material].kind=\"sidecar_hdf5\" requires a "
                  "[velocity] block (legacy CVM-H / CVM-S dataset path)");
      MFEM_VERIFY(!cfg.velocity.dataset_root.empty()
                  || !cfg.velocity.override_path.empty(),
                  "[velocity] must set either 'dataset_root' "
                  "(for model-based resolution) or 'override_path'");
   }

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
         if (r.kind == SpatialRule::Kind::BoxcarTaper)
         {
            // SCEC C∞ boxcar blend (symmetric with the RS path).
            const real_t B = r.BoxcarTaperFactor(x, y, z);
            if (B > 0.0)
            {
               if (!std::isnan(r.mu_s))     { mu_s_i += (r.mu_s     - mu_s_i) * B; }
               if (!std::isnan(r.mu_d))     { mu_d_i += (r.mu_d     - mu_d_i) * B; }
               if (!std::isnan(r.d_c))      { d_c_i  += (r.d_c      - d_c_i)  * B; }
               if (!std::isnan(r.cohesion)) { coh_i  += (r.cohesion - coh_i)  * B; }
            }
            continue;
         }
         if (!std::isnan(r.mu_s))     { mu_s_i = r.mu_s; }
         if (!std::isnan(r.mu_d))     { mu_d_i = r.mu_d; }
         if (!std::isnan(r.d_c))      { d_c_i  = r.d_c; }
         if (!std::isnan(r.cohesion)) { coh_i  = r.cohesion; }
         // Depth-linear cohesion taper (TPV31).  Takes precedence over
         // any constant `cohesion` set above; clamped from below by
         // `cohesion_floor_pa` (typically 0).
         if (!std::isnan(r.cohesion_grad_pa_per_m))
         {
            MFEM_VERIFY(!std::isnan(r.cohesion_ref_depth_m),
                        "ResolveSlipWeakening: cohesion_grad_pa_per_m set "
                        "but cohesion_ref_depth_m missing in spatial rule");
            // For axis 'z' (canonical SEAS), interpret the axis value as
            // POSITIVE depth (depth = max(0, -z)) so the taper formula
            // `cohesion = floor + grad * (ref_depth - depth)` can be
            // written with POSITIVE `ref_depth` and `grad` regardless of
            // mesh-z sign — matches the MakeDepthProfile1DMaterial
            // convention.  Axes 'x' and 'y' keep the legacy "axis_val is
            // used literally" semantic for backwards compatibility with
            // pre-rotation TPV31 (which encoded depth on +y).
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

   RateStatePerDOFParams p;
   p.a.SetSize(N);   p.b.SetSize(N);   p.Dc.SetSize(N);
   p.V_init.SetSize(N);
   p.f_0.SetSize(N); p.V_0.SetSize(N); p.eta.SetSize(N);
   p.sigma_n_eff.SetSize(N);
   const bool fvw =
      (cfg.state_evolution == StateEvolutionKind::SlipLawStrongRateWeakening);
   if (fvw) { p.V_w.SetSize(N); }

   for (int i = 0; i < N; ++i)
   {
      const real_t x = dof_coords_3d(3*i + 0);
      const real_t y = dof_coords_3d(3*i + 1);
      const real_t z = dof_coords_3d(3*i + 2);

      // Seed defaults.
      real_t a_i      = cfg.a_default;
      real_t b_i      = cfg.b_default;
      real_t Dc_i     = cfg.Dc_default;
      real_t V_init_i = cfg.V_init_default;
      real_t f0_i     = cfg.f_0_default;
      real_t V0_i     = cfg.V_0_default;
      real_t eta_i    = cfg.eta_default;
      real_t sn_i     = cfg.sigma_n_default;
      real_t Vw_i     = cfg.V_w_default;

      bool   eta_rule_set = false;

      for (const auto& r : cfg.spatial)
      {
         if (!r.matches(x, y, z, dof_to_attr[i])) { continue; }
         // RS resolver ignores Barrier (LSW-only concept).
         if (r.kind == SpatialRule::Kind::Barrier) { continue; }
         if (r.kind == SpatialRule::Kind::BoxcarTaper)
         {
            // SCEC C∞ boxcar blend: blend the current baseline toward
            // the override using `B = B(x)·B(y)·B(z)`.  `B = 1` ⇒
            // value = override (inside), `B = 0` ⇒ value unchanged
            // (outside).  No-op for parameters that this rule does
            // NOT set (override stays NaN).
            const real_t B = r.BoxcarTaperFactor(x, y, z);
            if (B > 0.0)
            {
               if (!std::isnan(r.a))       { a_i      += (r.a      - a_i)      * B; }
               if (!std::isnan(r.b))       { b_i      += (r.b      - b_i)      * B; }
               if (!std::isnan(r.Dc))      { Dc_i     += (r.Dc     - Dc_i)     * B; }
               if (!std::isnan(r.V_init))  { V_init_i += (r.V_init - V_init_i) * B; }
               if (!std::isnan(r.f_0))     { f0_i     += (r.f_0    - f0_i)     * B; }
               if (!std::isnan(r.V_0))     { V0_i     += (r.V_0    - V0_i)     * B; }
               if (!std::isnan(r.sigma_n)) { sn_i     += (r.sigma_n- sn_i)     * B; }
               if (!std::isnan(r.eta))
               { eta_i += (r.eta - eta_i) * B; eta_rule_set = true; }
               if (!std::isnan(r.V_w))     { Vw_i     += (r.V_w    - Vw_i)     * B; }
            }
            continue;
         }
         // Step-style rules (Depth / Box / RegionAttribute): direct
         // override assignment, last-match wins.
         if (!std::isnan(r.a))       { a_i      = r.a; }
         if (!std::isnan(r.b))       { b_i      = r.b; }
         if (!std::isnan(r.Dc))      { Dc_i     = r.Dc; }
         if (!std::isnan(r.V_init))  { V_init_i = r.V_init; }
         if (!std::isnan(r.f_0))     { f0_i     = r.f_0; }
         if (!std::isnan(r.V_0))     { V0_i     = r.V_0; }
         if (!std::isnan(r.sigma_n)) { sn_i     = r.sigma_n; }
         if (!std::isnan(r.eta))     { eta_i = r.eta;  eta_rule_set = true; }
         if (!std::isnan(r.V_w))     { Vw_i     = r.V_w; }
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
      const real_t depth = std::max(static_cast<real_t>(0.0), -z);
      const real_t P_p   = pp.P_p_pa + pp.P_p_grad_pa_per_m * depth;
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
      MFEM_VERIFY(b_i > 0.0, "ResolveRateState: b <= 0 at DOF " << i);
      // `a < b` is the velocity-weakening criterion; SCEC TPV102/104
      // intentionally use `a > b` outside the VW core (stable VS
      // regions).  Do not enforce `a < b` per DOF — the regime
      // (VW vs VS) is a physical input, not a validation rule.
      MFEM_VERIFY(Dc_i > 0.0, "ResolveRateState: Dc <= 0 at DOF " << i);
      MFEM_VERIFY(V0_i > 0.0, "ResolveRateState: V_0 <= 0 at DOF " << i);
      MFEM_VERIFY(sigma_n_eff > 0.0,
                  "ResolveRateState: sigma_n_eff <= 0 at DOF " << i);
      MFEM_VERIFY(f0_i > 0.0 && f0_i < 1.0,
                  "ResolveRateState: f_0 not in (0,1) at DOF " << i);

      if (fvw)
      {
         MFEM_VERIFY(Vw_i > 0.0,
                     "ResolveRateState: V_w <= 0 at DOF " << i
                     << " (state_evolution=slip_law_srw)");
         p.V_w(i) = Vw_i;
      }

      p.a(i) = a_i;   p.b(i) = b_i;   p.Dc(i) = Dc_i;
      p.V_init(i) = V_init_i;
      p.f_0(i) = f0_i; p.V_0(i) = V0_i;
      p.eta(i) = eta_i;
      p.sigma_n_eff(i) = sigma_n_eff;
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

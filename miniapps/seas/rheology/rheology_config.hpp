// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// rheology_config.hpp — the Nitsche deck tables [rheology] and
// [initial_stress] (rheology plan rev 1, section 7).  Nitsche-only tables:
// the SAFS parser and its byte-exact contract never see them.
//
//    [rheology]                       # absent = linear elastic, as before
//    representation = "nodal_pm1"     # the only option (plan 4.3)
//    [rheology.drucker_prager]
//    enabled           = true
//    cohesion_pa       = 1.36e6
//    bulk_friction     = 0.1934       # nu = tan(phi)
//    relaxation_time_s = 0.03         # Tv; 0 = rate-independent
//    [rheology.viscoelastic]          # enabled = true refused until R7
//    enabled = false
//    [rheology.damage_breakage]       # enabled = true refused until R8
//    enabled = false
//
//    [initial_stress]
//    kind = "lithostatic_ratio"       # or "constant_tensor" (sigma_xx_pa ...)
//    rho_rock_kg_m3, rho_fluid_kg_m3, gravity_m_s2, surface_elevation_m,
//    b_xx, b_yy, b_zz, b_xy, b_yz, b_xz, taper_top_depth_m, taper_bottom_depth_m
//
// Every model table needs `enabled`; an enabled model needs every key, with
// no defaults, as the [[stress_region]] keys.  Unknown keys, missing or
// mistyped values and refused combinations return a message naming the key;
// the driver stops with it.  The parsers return "" on success.
//
// Scope (Chunhui Zhao, 2026-10-20).

#ifndef MFEM_SEAS_RHEOLOGY_CONFIG_HPP
#define MFEM_SEAS_RHEOLOGY_CONFIG_HPP

#include "mfem.hpp"
#include "rheology_model.hpp"
#include "rheology_stack.hpp"
#include "drucker_prager.hpp"
#include "viscoelastic.hpp"
#include "damage_breakage.hpp"
#include "initial_stress.hpp"

#include <toml.hpp>

#include <cmath>
#include <initializer_list>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace mfem
{
namespace seas
{

struct RheologyConfig
{
   bool present = false;                     ///< the deck has [rheology]
   std::string representation = "nodal_pm1";
   bool drucker_prager = false;
   DruckerPragerParams dp;

   /// Any model enabled (an empty stack is linear elasticity).
   bool Any() const { return drucker_prager; }
};

namespace rheology_toml
{

/// "" when every key of table `t` is in `allowed`.
inline std::string CheckKeys(const toml::value &t, const std::string &name,
                             std::initializer_list<const char *> allowed)
{
   for (const auto &kv : t.as_table())
   {
      bool known = false;
      for (const char *a : allowed) { known = known || kv.first == a; }
      if (!known)
      {
         std::ostringstream e;
         e << name << " unknown key '" << kv.first << "' (allowed:";
         for (const char *a : allowed) { e << " " << a; }
         e << ")";
         return e.str();
      }
   }
   return "";
}

/// Required finite number.
inline std::string GetReal(const toml::value &t, const std::string &name, const char *key,
                           real_t &out)
{
   if (!t.contains(key)) { return name + " needs " + key + " (no default)"; }
   const toml::value &v = t.at(key);
   if (v.is_floating()) { out = static_cast<real_t>(v.as_floating()); }
   else if (v.is_integer()) { out = static_cast<real_t>(v.as_integer()); }
   else { return name + " " + key + " must be a number"; }
   if (!std::isfinite(out)) { return name + " " + key + " must be finite"; }
   return "";
}

/// Required boolean `enabled`.
inline std::string GetEnabled(const toml::value &t, const std::string &name, bool &out)
{
   if (!t.contains("enabled")) { return name + " needs enabled = true or false"; }
   const toml::value &v = t.at("enabled");
   if (!v.is_boolean()) { return name + " enabled must be true or false"; }
   out = v.as_boolean();
   return "";
}

/// A table that is not implemented yet: only `enabled`, and it must be false.
inline std::string CheckPlanned(const toml::value &t, const std::string &name,
                                const char *what, const char *phase)
{
   if (!t.is_table()) { return name + " must be a table"; }
   std::string e = CheckKeys(t, name, {"enabled"});
   if (!e.empty())
   {
      return e + "; the " + what + " keys are defined in rheology plan phase " + phase;
   }
   bool on = false;
   if (!(e = GetEnabled(t, name, on)).empty()) { return e; }
   if (on)
   {
      return name + " enabled = true: " + what + " is not implemented yet "
             "(rheology plan phase " + phase + ")";
   }
   return "";
}

} // namespace rheology_toml

/// Reads [rheology] from the deck root into `cfg`.
inline std::string ParseRheologyConfig(const toml::value &root, RheologyConfig &cfg)
{
   using namespace rheology_toml;
   cfg = RheologyConfig();
   if (!root.contains("rheology")) { return ""; }
   const toml::value &t = root.at("rheology");
   const std::string name = "[rheology]";
   if (!t.is_table()) { return name + " must be a table"; }
   cfg.present = true;
   std::string e = CheckKeys(t, name, {"representation", "drucker_prager", "viscoelastic",
                                       "damage_breakage"
                                      });
   if (!e.empty()) { return e; }
   if (t.contains("representation"))
   {
      const toml::value &v = t.at("representation");
      if (!v.is_string()) { return name + " representation must be a string"; }
      cfg.representation = v.as_string();
      if (cfg.representation != "nodal_pm1")
      {
         return name + " representation = \"" + cfg.representation
                + "\": only \"nodal_pm1\" is implemented";
      }
   }
   if (t.contains("viscoelastic"))
   {
      e = CheckPlanned(t.at("viscoelastic"), "[rheology.viscoelastic]",
                       "viscoelastic attenuation", kViscoelasticPhase);
      if (!e.empty()) { return e; }
   }
   if (t.contains("damage_breakage"))
   {
      e = CheckPlanned(t.at("damage_breakage"), "[rheology.damage_breakage]",
                       "continuum damage-breakage", kDamageBreakagePhase);
      if (!e.empty()) { return e; }
   }
   if (t.contains("drucker_prager"))
   {
      const toml::value &d = t.at("drucker_prager");
      const std::string dn = "[rheology.drucker_prager]";
      if (!d.is_table()) { return dn + " must be a table"; }
      e = CheckKeys(d, dn, {"enabled", "cohesion_pa", "bulk_friction", "relaxation_time_s"});
      if (!e.empty()) { return e; }
      if (!(e = GetEnabled(d, dn, cfg.drucker_prager)).empty()) { return e; }
      if (cfg.drucker_prager)
      {
         DruckerPragerParams &p = cfg.dp;
         if (!(e = GetReal(d, dn, "cohesion_pa", p.cohesion)).empty()) { return e; }
         if (!(e = GetReal(d, dn, "bulk_friction", p.bulk_friction)).empty()) { return e; }
         if (!(e = GetReal(d, dn, "relaxation_time_s", p.relaxation_time)).empty()) { return e; }
         if (p.cohesion < 0.0) { return dn + " cohesion_pa must be >= 0"; }
         if (p.bulk_friction < 0.0) { return dn + " bulk_friction must be >= 0"; }
         if (p.relaxation_time < 0.0) { return dn + " relaxation_time_s must be >= 0"; }
      }
   }
   std::vector<RheologyKind> kinds;
   if (cfg.drucker_prager) { kinds.push_back(RheologyKind::DruckerPrager); }
   return RheologyStack::CheckCombination(kinds);
}

/// The stack of the enabled models, in stack order (empty when none).
inline std::unique_ptr<RheologyStack> BuildRheologyStack(const RheologyConfig &cfg)
{
   std::unique_ptr<RheologyStack> s(new RheologyStack);
   if (cfg.drucker_prager)
   {
      s->Add(RheologyKind::DruckerPrager,
             std::unique_ptr<IRheologyModel>(new DruckerPrager(cfg.dp)));
   }
   return s;
}

/// Reads [initial_stress] from the deck root.  `present` is false (and the
/// spec untouched) when the table is absent; `up` is [fault_geometry] up.
inline std::string ParseInitialStressConfig(const toml::value &root, const real_t up[3],
                                            bool &present, InitialStressSpec &spec)
{
   using namespace rheology_toml;
   present = root.contains("initial_stress");
   if (!present) { return ""; }
   const toml::value &t = root.at("initial_stress");
   const std::string name = "[initial_stress]";
   if (!t.is_table()) { return name + " must be a table"; }
   if (!t.contains("kind") || !t.at("kind").is_string())
   {
      return name + " needs kind = \"lithostatic_ratio\" or \"constant_tensor\"";
   }
   const std::string kind = t.at("kind").as_string();
   InitialStressSpec s;
   for (int i = 0; i < 3; i++) { s.up[i] = up[i]; }
   std::string e;
   if (kind == "constant_tensor")
   {
      s.kind = InitialStressSpec::Kind::ConstantTensor;
      e = CheckKeys(t, name, {"kind", "sigma_xx_pa", "sigma_yy_pa", "sigma_zz_pa",
                              "sigma_xy_pa", "sigma_yz_pa", "sigma_xz_pa"
                             });
      if (!e.empty()) { return e; }
      const char *k[6] = {"sigma_xx_pa", "sigma_yy_pa", "sigma_zz_pa",
                          "sigma_xy_pa", "sigma_yz_pa", "sigma_xz_pa"
                         };
      real_t v[6];
      for (int i = 0; i < 6; i++) { if (!(e = GetReal(t, name, k[i], v[i])).empty()) { return e; } }
      sym6::ToRowMajor(v, s.sigma);
   }
   else if (kind == "lithostatic_ratio")
   {
      s.kind = InitialStressSpec::Kind::LithostaticRatio;
      e = CheckKeys(t, name, {"kind", "rho_rock_kg_m3", "rho_fluid_kg_m3", "gravity_m_s2",
                              "surface_elevation_m", "b_xx", "b_yy", "b_zz", "b_xy", "b_yz",
                              "b_xz", "taper_top_depth_m", "taper_bottom_depth_m"
                             });
      if (!e.empty()) { return e; }
      struct { const char *key; real_t *dst; } scalars[] =
      {
         {"rho_rock_kg_m3", &s.rho_rock}, {"rho_fluid_kg_m3", &s.rho_fluid},
         {"gravity_m_s2", &s.gravity}, {"surface_elevation_m", &s.surface_elevation},
         {"taper_top_depth_m", &s.taper_top_depth},
         {"taper_bottom_depth_m", &s.taper_bottom_depth}
      };
      for (auto &sc : scalars)
      {
         if (!(e = GetReal(t, name, sc.key, *sc.dst)).empty()) { return e; }
      }
      const char *k[6] = {"b_xx", "b_yy", "b_zz", "b_xy", "b_yz", "b_xz"};
      real_t b[6];
      for (int i = 0; i < 6; i++) { if (!(e = GetReal(t, name, k[i], b[i])).empty()) { return e; } }
      sym6::ToRowMajor(b, s.B);
   }
   else
   {
      return name + " kind = \"" + kind + "\": use \"lithostatic_ratio\" or \"constant_tensor\"";
   }
   if (!(e = InitialStress::Validate(s)).empty()) { return e; }
   spec = s;
   return "";
}

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_RHEOLOGY_CONFIG_HPP

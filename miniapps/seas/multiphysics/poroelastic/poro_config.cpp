// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// poro_config.cpp — parser of the poroelastic option (poro_config.hpp).

#include "poro_config.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <sstream>

namespace mfem
{
namespace seas
{
namespace poro
{

namespace
{

const char *kMain = "[multiphysics.poroelastic]";
const char *kInit = "[multiphysics.poroelastic.initial_pressure]";
const char *kSrc = "[multiphysics.poroelastic.source]";

std::string Msg(const char *table, const std::string &key, const std::string &what)
{
   return std::string(table) + " " + key + ": " + what;
}

/// First key of t (sorted, so the message does not depend on hash order)
/// that is not in `known`, or "".
std::string UnknownKey(const toml::value &t, const std::set<std::string> &known)
{
   std::vector<std::string> keys;
   for (const auto &kv : t.as_table()) { keys.push_back(kv.first); }
   std::sort(keys.begin(), keys.end());
   for (const auto &k : keys) { if (!known.count(k)) { return k; } }
   return "";
}

bool GetNumber(const toml::value &t, const char *table, const char *key, real_t &out,
               std::string &err)
{
   if (!t.contains(key)) { err = Msg(table, key, "missing (every key is required)"); return false; }
   const auto &v = t.at(key);
   if (v.is_floating()) { out = static_cast<real_t>(v.as_floating()); }
   else if (v.is_integer()) { out = static_cast<real_t>(v.as_integer()); }
   else { err = Msg(table, key, "must be a number"); return false; }
   if (!std::isfinite(out)) { err = Msg(table, key, "must be finite"); return false; }
   return true;
}

bool GetString(const toml::value &t, const char *table, const char *key,
               const std::vector<std::string> &allowed, std::string &out, std::string &err)
{
   if (!t.contains(key)) { err = Msg(table, key, "missing (every key is required)"); return false; }
   const auto &v = t.at(key);
   if (!v.is_string()) { err = Msg(table, key, "must be a string"); return false; }
   out = v.as_string();
   if (std::find(allowed.begin(), allowed.end(), out) == allowed.end())
   {
      std::ostringstream os;
      os << "\"" << out << "\" is not one of";
      for (const auto &a : allowed) { os << " \"" << a << "\""; }
      err = Msg(table, key, os.str());
      return false;
   }
   return true;
}

bool GetPoint(const toml::value &t, const char *table, const char *key, real_t out[3],
              std::string &err)
{
   if (!t.contains(key)) { err = Msg(table, key, "missing (every key is required)"); return false; }
   const auto &v = t.at(key);
   if (!v.is_array() || v.as_array().size() != 3)
   {
      err = Msg(table, key, "must be an array of three numbers [x, y, z] in m");
      return false;
   }
   for (int c = 0; c < 3; c++)
   {
      const auto &e = v.as_array()[c];
      if (e.is_floating()) { out[c] = static_cast<real_t>(e.as_floating()); }
      else if (e.is_integer()) { out[c] = static_cast<real_t>(e.as_integer()); }
      else { err = Msg(table, key, "must be an array of three numbers [x, y, z] in m"); return false; }
      if (!std::isfinite(out[c])) { err = Msg(table, key, "must be finite"); return false; }
   }
   return true;
}

bool Positive(const char *table, const char *key, real_t v, std::string &err)
{
   if (v > 0.0) { return true; }
   std::ostringstream os;
   os << "must be > 0; got " << v;
   err = Msg(table, key, os.str());
   return false;
}

std::string ParseInitialPressure(const toml::value &t, PoroInitialPressure &ip)
{
   std::string err;
   if (!t.is_table()) { return std::string(kInit) + " must be a table"; }
   const std::string unk = UnknownKey(t, {"origin_m", "magnitude_pa", "width_per_m2"});
   if (!unk.empty()) { return Msg(kInit, unk, "unknown key"); }
   if (!GetPoint(t, kInit, "origin_m", ip.origin_m, err)) { return err; }
   if (!GetNumber(t, kInit, "magnitude_pa", ip.magnitude_pa, err)) { return err; }
   if (!GetNumber(t, kInit, "width_per_m2", ip.width_per_m2, err)) { return err; }
   if (!Positive(kInit, "width_per_m2", ip.width_per_m2, err)) { return err; }
   ip.on = true;
   return "";
}

std::string ParseSource(const toml::value &t, PoroSource &s)
{
   std::string err;
   if (!t.is_table()) { return std::string(kSrc) + " must be a table"; }
   const std::string unk = UnknownKey(t, {"shape", "origin_m", "width_per_m2", "rate",
                                          "time_function", "rise_time_s"});
   if (!unk.empty()) { return Msg(kSrc, unk, "unknown key"); }
   if (!GetString(t, kSrc, "shape", {"gaussian", "well"}, s.shape, err)) { return err; }
   if (!GetPoint(t, kSrc, "origin_m", s.origin_m, err)) { return err; }
   if (s.shape == "gaussian")
   {
      if (!GetNumber(t, kSrc, "width_per_m2", s.width_per_m2, err)) { return err; }
      if (!Positive(kSrc, "width_per_m2", s.width_per_m2, err)) { return err; }
   }
   else if (t.contains("width_per_m2"))
   {
      return Msg(kSrc, "width_per_m2", "applies to shape = \"gaussian\" only");
   }
   if (!GetNumber(t, kSrc, "rate", s.rate, err)) { return err; }
   if (!GetString(t, kSrc, "time_function", {"step", "ramp", "smooth_step"}, s.time_function,
                  err)) { return err; }
   if (s.time_function == "step")
   {
      if (t.contains("rise_time_s"))
      {
         return Msg(kSrc, "rise_time_s", "applies to time_function = \"ramp\" or "
                    "\"smooth_step\" only");
      }
   }
   else
   {
      if (!GetNumber(t, kSrc, "rise_time_s", s.rise_time_s, err)) { return err; }
      if (!Positive(kSrc, "rise_time_s", s.rise_time_s, err)) { return err; }
   }
   s.on = true;
   return "";
}

std::string ParseMain(const toml::value &t, PoroConfig &c)
{
   std::string err;
   const std::string unk = UnknownKey(t, {"enabled", "porosity", "permeability_m2",
                                          "viscosity_pa_s", "tortuosity", "rho_fluid",
                                          "rho_solid", "bulk_fluid_pa", "bulk_solid_pa",
                                          "fault_hydraulics", "background_p_pa",
                                          "fluid_penalty", "drag_gamma",
                                          "velocity_coupling", "initial_pressure",
                                          "source"});
   if (!unk.empty()) { return Msg(kMain, unk, "unknown key"); }
   if (!GetNumber(t, kMain, "porosity", c.porosity, err)) { return err; }
   if (!(c.porosity > 0.0 && c.porosity < 1.0))
   {
      std::ostringstream os;
      os << "must be in (0, 1); got " << c.porosity;
      return Msg(kMain, "porosity", os.str());
   }
   if (!GetNumber(t, kMain, "permeability_m2", c.permeability_m2, err)) { return err; }
   if (!Positive(kMain, "permeability_m2", c.permeability_m2, err)) { return err; }
   if (!GetNumber(t, kMain, "viscosity_pa_s", c.viscosity_pa_s, err)) { return err; }
   if (!Positive(kMain, "viscosity_pa_s", c.viscosity_pa_s, err)) { return err; }
   if (!GetNumber(t, kMain, "tortuosity", c.tortuosity, err)) { return err; }
   if (!(c.tortuosity >= 1.0))
   {
      std::ostringstream os;
      os << "must be >= 1; got " << c.tortuosity;
      return Msg(kMain, "tortuosity", os.str());
   }
   if (!GetNumber(t, kMain, "rho_fluid", c.rho_fluid, err)) { return err; }
   if (!Positive(kMain, "rho_fluid", c.rho_fluid, err)) { return err; }
   if (!GetNumber(t, kMain, "rho_solid", c.rho_solid, err)) { return err; }
   if (!Positive(kMain, "rho_solid", c.rho_solid, err)) { return err; }
   if (!GetNumber(t, kMain, "bulk_fluid_pa", c.bulk_fluid_pa, err)) { return err; }
   if (!Positive(kMain, "bulk_fluid_pa", c.bulk_fluid_pa, err)) { return err; }
   if (!GetNumber(t, kMain, "bulk_solid_pa", c.bulk_solid_pa, err)) { return err; }
   if (!Positive(kMain, "bulk_solid_pa", c.bulk_solid_pa, err)) { return err; }
   if (!GetString(t, kMain, "fault_hydraulics", {"permeable", "sealed"}, c.fault_hydraulics,
                  err)) { return err; }
   if (c.fault_hydraulics == "sealed")
   {
      return Msg(kMain, "fault_hydraulics", "\"sealed\" arrives in step P6 of the plan; "
                 "only \"permeable\" is implemented");
   }
   if (!GetNumber(t, kMain, "background_p_pa", c.background_p_pa, err)) { return err; }
   if (!GetNumber(t, kMain, "fluid_penalty", c.fluid_penalty, err)) { return err; }
   if (!Positive(kMain, "fluid_penalty", c.fluid_penalty, err)) { return err; }
   if (!GetString(t, kMain, "drag_gamma", {"fitted", "half"}, c.drag_gamma, err)) { return err; }
   if (!GetString(t, kMain, "velocity_coupling", {"momentum", "diagonal"},
                  c.velocity_coupling, err)) { return err; }
   if (t.contains("initial_pressure"))
   {
      err = ParseInitialPressure(t.at("initial_pressure"), c.initial_pressure);
      if (!err.empty()) { return err; }
   }
   if (t.contains("source"))
   {
      err = ParseSource(t.at("source"), c.source);
      if (!err.empty()) { return err; }
   }
   return "";
}

} // namespace

PoroParseResult ParsePoroConfig(const toml::value &root, int cli_flag)
{
   PoroParseResult r;
   if (cli_flag != -1 && cli_flag != 0 && cli_flag != 1)
   {
      std::ostringstream os;
      os << "--poroelastic " << cli_flag << ": expected 0 or 1";
      r.error = os.str();
      return r;
   }
   const toml::value *t = nullptr;
   if (root.contains("multiphysics"))
   {
      const toml::value &mp = root.at("multiphysics");
      if (!mp.is_table()) { r.error = "[multiphysics] must be a table"; return r; }
      const std::string unk = UnknownKey(mp, {"poroelastic"});
      if (!unk.empty()) { r.error = "[multiphysics] " + unk + ": unknown key"; return r; }
      if (mp.contains("poroelastic"))
      {
         t = &mp.at("poroelastic");
         if (!t->is_table())
         {
            r.error = std::string(kMain) + " must be a table (switch it on with "
                      "enabled = true inside the table)";
            return r;
         }
      }
   }
   bool deck_on = false;
   if (t)
   {
      if (!t->contains("enabled"))
      {
         r.error = Msg(kMain, "enabled", "missing (every key is required)");
         return r;
      }
      if (!t->at("enabled").is_boolean())
      {
         r.error = Msg(kMain, "enabled", "must be true or false");
         return r;
      }
      deck_on = t->at("enabled").as_boolean();
   }
   const bool on = (cli_flag == -1) ? deck_on : (cli_flag == 1);
   if (on && !t)
   {
      r.error = "--poroelastic 1 needs the [multiphysics.poroelastic] table";
      return r;
   }
   if (!on)
   {
      if (t)
      {
         r.cfg.warning = std::string(kMain) + " is present but the option is off ("
                         + (cli_flag == 0 ? "--poroelastic 0" : "enabled = false")
                         + "); the table is ignored";
      }
      return r;
   }
   r.cfg.enabled = true;
   r.error = ParseMain(*t, r.cfg);
   if (!r.error.empty()) { r.cfg.enabled = false; }
   return r;
}

std::string CheckBulkDensity(const PoroConfig &cfg, real_t rho_fallback)
{
   const real_t rho = cfg.BulkDensity();
   if (std::abs(rho_fallback - rho) <= 1e-12 * rho) { return ""; }
   std::ostringstream os;
   os.precision(17);
   os << "[material_constant_fallback] rho = " << rho_fallback << " differs from the "
      "poroelastic bulk density (1 - porosity) rho_solid + porosity rho_fluid = " << rho
      << " (relative difference " << (rho_fallback - rho) / rho << " > 1e-12); with the "
      "option on that table holds the drained frame, so rho must be the bulk density";
   return os.str();
}

} // namespace poro
} // namespace seas
} // namespace mfem

// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// pff_config_parser.hpp — TOML -> CaseConfig (see pff_config.hpp for the tables).
//
// toml11 v3 (extern/toml11, the submodule commit seas uses).  Reals accept TOML
// integers too (toml11's find_or<double> would otherwise fall back to the
// default on `E = 50`).  Unknown tables and keys throw.

#ifndef MFEM_PFF_CONFIG_PARSER_HPP
#define MFEM_PFF_CONFIG_PARSER_HPP

#include "pff_config.hpp"
#include "toml_reader.hpp"

#include <set>
#include <sstream>
#include <stdexcept>

namespace mfem
{
namespace pff
{

class CaseConfigParser : private TomlReader
{
public:
   static CaseConfig ParseFile(const std::string &path)
   {
      return Parse(toml::parse(path));
   }

   static CaseConfig ParseString(const std::string &s, const std::string &name = "string")
   {
      std::istringstream iss(s);
      return Parse(toml::parse(iss, name));
   }

private:
   static CaseConfig Parse(const toml::value &data)
   {
      CaseConfig c;
      Check(data, "root", {"case", "mesh", "boundaries", "material", "phase_field", "pulse",
                           "confinement", "damping", "time", "coupling", "poro", "output",
                           "metrics"
                          });

      {
         const auto &t = data.at("case");
         Check(t, "case", {"name", "model", "source"});
         c.name = Str(t, "name", "");
         const std::string m = Str(t, "model", "solid");
         if (m == "solid") { c.model = CaseModel::Solid; }
         else if (m == "hydromech") { c.model = CaseModel::Hydromech; }
         else { throw std::invalid_argument("[case] model must be solid or hydromech"); }
         c.source = Str(t, "source", "");
      }
      {
         const auto &t = data.at("mesh");
         Check(t, "mesh", {"file", "dim", "uniform_refinements", "corner"});
         c.mesh_file = Str(t, "file", "");
         c.dim = Int(t, "dim", 3);
         c.uniform_refinements = Int(t, "uniform_refinements", 0);
         c.corner = Vec3(t, "corner", c.corner);
      }
      if (data.contains("boundaries"))
      {
         const auto &t = data.at("boundaries");
         Check(t, "boundaries", {"pressure", "absorbing", "fixed_z", "confinement", "drained"});
         c.pressure_bdr = Names(t, "pressure");
         c.absorbing_bdr = Names(t, "absorbing");
         c.fixed_z_bdr = Names(t, "fixed_z");
         c.confinement_bdr = Names(t, "confinement");
         c.drained_bdr = Names(t, "drained");
      }
      {
         const auto &t = data.at("material");
         Check(t, "material", {"E", "nu", "density"});
         c.E = Real(t, "E", 0.0);
         c.nu = Real(t, "nu", 0.0);
         c.density = Real(t, "density", 0.0);
      }
      {
         const auto &t = data.at("phase_field");
         Check(t, "phase_field", {"model", "Gc", "ft", "l", "eta", "p", "decomposition",
                                  "driving", "initial_damage"
                                 });
         c.pf.crack.model = ParseCrackModel(Str(t, "model", "AT1"));
         c.pf.l = Real(t, "l", 0.0);
         c.pf.deg.eta = Real(t, "eta", 1e-6);
         c.pf.deg.p = Real(t, "p", 2.0);
         c.ft = Real(t, "ft", 0.0);
         if (t.contains("Gc") && t.contains("ft"))
         {
            throw std::invalid_argument("[phase_field] give Gc or ft, not both");
         }
         c.pf.Gc = t.contains("Gc") ? Real(t, "Gc", 0.0)
                   : PhaseFieldModel::AT1ToughnessFromStrength(c.E, c.ft, c.pf.l);
         c.decomposition = ParseDecomposition(Str(t, "decomposition", "spectral"));
         c.pf.driving = ParseDrivingEnergy(Str(t, "driving", "active"));
         if (t.contains("initial_damage"))
         {
            for (const auto &b : toml::get<toml::array>(t.at("initial_damage")))
            {
               Check(b, "phase_field.initial_damage", {"lo", "hi", "value"});
               DamageBox box;
               box.lo = Vec3(b, "lo", box.lo);
               box.hi = Vec3(b, "hi", box.hi);
               box.value = Real(b, "value", 0.0);
               c.initial_damage.push_back(box);
            }
         }
      }
      if (data.contains("pulse"))
      {
         const auto &t = data.at("pulse");
         Check(t, "pulse", {"shape_alpha", "shape_beta", "rise_time", "pulse_duration", "EM",
                            "gap", "efficiency", "kappa", "m_exp", "base_factor", "center",
                            "num_pulses", "peak_pressure", "min_pressure"
                           });
         PulseLoadParams &p = c.pulse;
         p.shape_alpha = Real(t, "shape_alpha", 0.0);
         p.shape_beta = Real(t, "shape_beta", 0.0);
         p.rise_time = Real(t, "rise_time", 0.0);
         p.pulse_duration = Real(t, "pulse_duration", 0.0);
         p.EM = Real(t, "EM", 0.0);
         p.gap = Real(t, "gap", 0.0);
         p.efficiency = Real(t, "efficiency", 1.0);
         p.kappa = Real(t, "kappa", 0.35);
         p.m_exp = Real(t, "m_exp", 1.0);
         p.base_factor = Real(t, "base_factor", 9000.0);
         const auto ctr = Vec3(t, "center", {{0.0, 0.0, 0.0}});
         for (int i = 0; i < 3; i++) { p.center[i] = ctr[i]; }
         p.num_pulses = Int(t, "num_pulses", 1);
         p.peak_pressure = Real(t, "peak_pressure", 0.0);
         if (t.contains("min_pressure"))
         {
            p.use_min_pressure = true;
            p.min_pressure = Real(t, "min_pressure", 0.0);
         }
      }
      if (data.contains("confinement"))
      {
         const auto &t = data.at("confinement");
         Check(t, "confinement", {"pressure"});
         c.confinement_pressure = Real(t, "pressure", 0.0);
      }
      if (data.contains("damping"))
      {
         const auto &t = data.at("damping");
         Check(t, "damping", {"zeta_pulse", "zeta_quench", "t_free", "quench", "alpha_mass"});
         c.zeta_pulse = Real(t, "zeta_pulse", 0.0);
         c.zeta_quench = Real(t, "zeta_quench", 0.0);
         c.t_free = Real(t, "t_free", 0.0);
         c.alpha_mass = Real(t, "alpha_mass", 0.0);
         const std::string q = Str(t, "quench", "stiffness");
         if (q == "stiffness") { c.quench = QuenchMode::Stiffness; }
         else if (q == "mass") { c.quench = QuenchMode::Mass; }
         else if (q == "implicit") { c.quench = QuenchMode::Implicit; }
         else { throw std::invalid_argument("[damping] quench must be stiffness, mass or implicit"); }
      }
      {
         const auto &t = data.at("time");
         Check(t, "time", {"end_time", "cfl", "dt_max", "dt_quiet", "mass"});
         c.end_time = Real(t, "end_time", c.pulse.num_pulses * c.pulse.pulse_duration);
         c.cfl = Real(t, "cfl", 0.4);
         c.dt_max = Real(t, "dt_max", 0.0);
         c.dt_quiet = Real(t, "dt_quiet", c.dt_max);
         const std::string ms = Str(t, "mass", "consistent");
         if (ms == "consistent") { c.consistent_mass = true; }
         else if (ms == "lumped") { c.consistent_mass = false; }
         else { throw std::invalid_argument("[time] mass must be consistent or lumped"); }
      }
      if (data.contains("coupling"))
      {
         const auto &t = data.at("coupling");
         Check(t, "coupling", {"scheme", "max_passes", "tolerance"});
         const std::string s = Str(t, "scheme", "single_pass");
         if (s == "single_pass") { c.coupling = CouplingScheme::SinglePass; }
         else if (s == "fixed_point") { c.coupling = CouplingScheme::FixedPoint; }
         else { throw std::invalid_argument("[coupling] scheme must be single_pass or fixed_point"); }
         c.max_passes = Int(t, "max_passes", 1);
         c.coupling_tol = Real(t, "tolerance", 1e-6);
      }
      if (data.contains("poro"))
      {
         const auto &t = data.at("poro");
         Check(t, "poro", {"K_s", "K_f", "phi0", "phi_lb", "phi_ub", "g_min", "viscosity",
                           "fluid_density", "initial_pore_pressure", "porosity_model",
                           "permeability_model", "k0", "perm_exponent", "coeff_b", "wc",
                           "d_threshold", "fc", "w_res", "anisotropic", "normal_source",
                           "regularize_normal", "normal_reg_eps"
                          });
         PoroDamageParams &P = c.poro;
         P.K = c.K();
         P.K_s = Real(t, "K_s", 0.0);
         P.K_f = Real(t, "K_f", 0.0);
         P.phi0 = Real(t, "phi0", 0.0);
         P.phi_lb = Real(t, "phi_lb", P.phi0);
         P.phi_ub = Real(t, "phi_ub", 1.0);
         P.g_min = Real(t, "g_min", c.pf.deg.eta);
         P.viscosity = Real(t, "viscosity", 1e-3);
         P.porosity_model = ParsePorosityModel(Str(t, "porosity_model", "damage"));
         P.perm_model = ParsePermeabilityModel(Str(t, "permeability_model", "none"));
         P.k0 = Real(t, "k0", 0.0);
         P.perm_exponent = Real(t, "perm_exponent", 1.0);
         P.coeff_b = Real(t, "coeff_b", 0.0);
         // darcy_poiseuille default w_c = Gc / ft, ft from the AT1 relation (MOOSE inputs)
         P.wc = Real(t, "wc", c.pf.Gc / PhaseFieldModel::AT1Strength(c.E, c.pf.Gc, c.pf.l));
         P.d_threshold = Real(t, "d_threshold", 0.5);
         P.fc = Real(t, "fc", 1.0);
         P.w_res = Real(t, "w_res", 0.0);
         P.anisotropic = Bool(t, "anisotropic", true);
         P.normal_source = ParseCrackNormalSource(Str(t, "normal_source", "damage_gradient"));
         P.regularize_normal = Bool(t, "regularize_normal", false);
         P.normal_reg_eps = Real(t, "normal_reg_eps", 1e-8);
         c.fluid_density = Real(t, "fluid_density", 1000.0);
         c.initial_pore_pressure = Real(t, "initial_pore_pressure", 0.0);
      }
      else if (c.model == CaseModel::Hydromech)
      {
         throw std::invalid_argument("model = hydromech needs a [poro] table");
      }
      if (data.contains("output"))
      {
         const auto &t = data.at("output");
         Check(t, "output", {"directory", "interval", "fields"});
         c.output_dir = Str(t, "directory", c.output_dir);
         c.output_interval = Int(t, "interval", c.output_interval);
         if (t.contains("fields"))
         {
            c.output_fields.clear();
            for (const auto &v : toml::get<toml::array>(t.at("fields")))
            {
               c.output_fields.push_back(toml::get<std::string>(v));
            }
         }
      }
      if (data.contains("metrics"))
      {
         const auto &t = data.at("metrics");
         Check(t, "metrics", {"d_threshold", "r_core", "r_hole"});
         c.d_threshold = Real(t, "d_threshold", 0.5);
         c.r_core = Real(t, "r_core", 0.0);
         c.r_hole = Real(t, "r_hole", 0.0);
      }
      return c;
   }
};

} // namespace pff
} // namespace mfem

#endif

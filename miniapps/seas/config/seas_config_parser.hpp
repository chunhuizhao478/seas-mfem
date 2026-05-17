// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// This file is part of the MFEM library. For more information and source code
// availability visit https://mfem.org.
//
// MFEM is free software; you can redistribute it and/or modify it under the
// terms of the BSD-3 license. We welcome feedback and contributions, see file
// CONTRIBUTING.md for details.

#ifndef MFEM_SEAS_CONFIG_PARSER_HPP
#define MFEM_SEAS_CONFIG_PARSER_HPP

#include "seas_config.hpp"
#include "bp5_params.hpp"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#ifdef SEAS_USE_TOML
#include <toml.hpp>
#endif

namespace mfem
{
namespace seas
{

/// @brief Parse a TOML configuration file into SEASConfig.
///
/// If `benchmark = "bp5"` is set, fills defaults from BP5Params first,
/// then any explicit fields override them.
/// Unknown sections produce a warning, not an error.
///
/// Requires SEAS_USE_TOML to be defined. Without it, ParseSEASConfig
/// aborts with an error message.
class SEASConfigParser
{
public:
   static SEASConfig ParseFile(const std::string &filepath)
   {
#ifdef SEAS_USE_TOML
      auto data = toml::parse(filepath);
      return ParseToml(data);
#else
      MFEM_ABORT("SEAS TOML config requires SEAS_USE_TOML. "
                 "Rebuild with -DSEAS_USE_TOML and toml11 in extern/.");
      return {};
#endif
   }

   static SEASConfig ParseString(const std::string &toml_str)
   {
#ifdef SEAS_USE_TOML
      std::istringstream iss(toml_str);
      auto data = toml::parse(iss);
      return ParseToml(data);
#else
      MFEM_ABORT("SEAS TOML config requires SEAS_USE_TOML.");
      return {};
#endif
   }

   /// Apply CLI overrides in "key=value" format (e.g., "material.density=3000")
   static real_t ParseDouble(const std::string &key, const std::string &val)
   {
      try { return std::stod(val); }
      catch (const std::exception &)
      {
         MFEM_ABORT("Invalid numeric value for " << key << ": '" << val << "'");
         return 0;
      }
   }

   static int ParseInt(const std::string &key, const std::string &val)
   {
      try { return std::stoi(val); }
      catch (const std::exception &)
      {
         MFEM_ABORT("Invalid integer value for " << key << ": '" << val << "'");
         return 0;
      }
   }

   static void ApplyCLIOverrides(SEASConfig &config,
                                 const std::vector<std::string> &overrides)
   {
      for (const auto &ov : overrides)
      {
         auto eq = ov.find('=');
         if (eq == std::string::npos)
         {
            std::cerr << "WARNING: CLI override missing '=': " << ov << "\n";
            continue;
         }
         std::string key = ov.substr(0, eq);
         std::string val = ov.substr(eq + 1);

         if (key == "material.density") { config.material.density = ParseDouble(key, val); }
         else if (key == "material.cs") { config.material.cs = ParseDouble(key, val); }
         else if (key == "material.nu") { config.material.nu = ParseDouble(key, val); }
         else if (key == "friction.V0") { config.friction.V0 = ParseDouble(key, val); }
         else if (key == "friction.f0") { config.friction.f0 = ParseDouble(key, val); }
         else if (key == "friction.b") { config.friction.b = ParseDouble(key, val); }
         else if (key == "friction.L0") { config.friction.L0 = ParseDouble(key, val); }
         else if (key == "friction.sigma_n") { config.friction.sigma_n = ParseDouble(key, val); }
         else if (key == "loading.Vp") { config.loading.Vp = ParseDouble(key, val); }
         else if (key == "loading.V_nuc") { config.loading.V_nuc = ParseDouble(key, val); }
         else if (key == "time.t_final") { config.time.t_final = ParseDouble(key, val); }
         else if (key == "time.max_steps") { config.time.max_steps = ParseInt(key, val); }
         else if (key == "time.use_petsc_ts") { config.time.use_petsc_ts = (val == "true" || val == "1"); }
         else if (key == "time.petsc_ts_options") { config.time.petsc_ts_options = val; }
         else if (key == "mesh.file") { config.mesh.file = val; }
         else if (key == "mesh.scale") { config.mesh.scale = ParseDouble(key, val); }
         else if (key == "mesh.order") { config.mesh.order = ParseInt(key, val); }
         else if (key == "solver.dg_method") { config.solver.dg_method = val; }
         else if (key == "output.output_dir") { config.output.output_dir = val; }
         else if (key == "output.output_prefix") { config.output.output_prefix = val; }
         else
         {
            std::cerr << "WARNING: Unknown CLI override key: " << key << "\n";
         }
      }
   }

   /// Validate a parsed config. Aborts on invalid values.
   static void Validate(const SEASConfig &config)
   {
      MFEM_VERIFY(!config.mesh.file.empty(),
                  "Config validation: mesh.file is required");
      MFEM_VERIFY(config.material.density > 0,
                  "Config validation: material.density must be positive");
      MFEM_VERIFY(config.material.cs > 0,
                  "Config validation: material.cs must be positive");
      MFEM_VERIFY(config.material.nu > 0 && config.material.nu < 0.5,
                  "Config validation: material.nu must be in (0, 0.5)");
      MFEM_VERIFY(config.friction.V0 > 0,
                  "Config validation: friction.V0 must be positive");
      MFEM_VERIFY(config.friction.sigma_n > 0,
                  "Config validation: friction.sigma_n must be positive");
      MFEM_VERIFY(config.loading.Vp > 0,
                  "Config validation: loading.Vp must be positive");
      MFEM_VERIFY(config.time.t_final > 0,
                  "Config validation: time.t_final must be positive");
      MFEM_VERIFY(!config.boundary.dirichlet_attrs.empty(),
                  "Config validation: boundary.dirichlet must not be empty "
                  "(specify attrs or use benchmark preset)");
      MFEM_VERIFY(config.boundary.fault_attr > 0,
                  "Config validation: boundary.fault must be set "
                  "(specify attr or use benchmark preset)");
      MFEM_VERIFY(config.time.atol > 0,
                  "Config validation: time.atol must be positive");
      MFEM_VERIFY(config.solver.dg_method == "IP" ||
                  config.solver.dg_method == "BR2",
                  "Config validation: solver.dg_method must be \"IP\" or "
                  "\"BR2\", got \"" << config.solver.dg_method << "\"");
      MFEM_VERIFY(config.solver.solver_type == "cg" ||
                  config.solver.solver_type == "mumps" ||
                  config.solver.solver_type == "mumps-blr" ||
                  config.solver.solver_type == "superlu" ||
                  config.solver.solver_type == "strumpack" ||
                  config.solver.solver_type == "gmres",
                  "Config validation: solver.solver_type must be one of "
                  "{cg, mumps, mumps-blr, superlu, strumpack, gmres}, got \""
                  << config.solver.solver_type << "\"");

      if (config.simulation.mode != "qd")
      {
         MFEM_ABORT("Simulation mode '" << config.simulation.mode
                    << "' requires the dynamic solver module "
                       "(not yet implemented). Use mode = \"qd\".");
      }
   }

private:
#ifdef SEAS_USE_TOML
   static SEASConfig ParseToml(const toml::value &data)
   {
      SEASConfig config;

      // Check for benchmark preset
      config.benchmark = toml::find_or<std::string>(data, "benchmark", "");
      if (config.benchmark == "bp5")
      {
         ApplyBP5Defaults(config);
      }

      // Parse sections — explicit values override preset defaults
      if (data.contains("mesh"))
      {
         const auto &m = data.at("mesh");
         config.mesh.file = toml::find_or<std::string>(m, "file", config.mesh.file);
         config.mesh.scale = toml::find_or<double>(m, "scale", config.mesh.scale);
         config.mesh.order = static_cast<int>(toml::find_or<toml::integer>(m, "order", config.mesh.order));
      }

      if (data.contains("material"))
      {
         const auto &m = data.at("material");
         config.material.density = toml::find_or<double>(m, "density", config.material.density);
         config.material.cs = toml::find_or<double>(m, "cs", config.material.cs);
         config.material.nu = toml::find_or<double>(m, "nu", config.material.nu);
      }

      if (data.contains("friction"))
      {
         const auto &f = data.at("friction");
         config.friction.V0 = toml::find_or<double>(f, "V0", config.friction.V0);
         config.friction.f0 = toml::find_or<double>(f, "f0", config.friction.f0);
         config.friction.b = toml::find_or<double>(f, "b", config.friction.b);
         config.friction.L0 = toml::find_or<double>(f, "L0", config.friction.L0);
         config.friction.L_nuc = toml::find_or<double>(f, "L_nuc", config.friction.L_nuc);
         config.friction.a0 = toml::find_or<double>(f, "a0", config.friction.a0);
         config.friction.amax = toml::find_or<double>(f, "amax", config.friction.amax);
         config.friction.sigma_n = toml::find_or<double>(f, "sigma_n", config.friction.sigma_n);
      }

      if (data.contains("loading"))
      {
         const auto &l = data.at("loading");
         config.loading.Vp = toml::find_or<double>(l, "Vp", config.loading.Vp);
         config.loading.V_init = toml::find_or<double>(l, "V_init", config.loading.V_init);
         config.loading.V_nuc = toml::find_or<double>(l, "V_nuc", config.loading.V_nuc);
         config.loading.delta_tau_factor = toml::find_or<double>(l, "delta_tau_factor", config.loading.delta_tau_factor);
         config.loading.nucleation_eps = toml::find_or<double>(l, "nucleation_eps", config.loading.nucleation_eps);
         config.loading.smooth_nucleation = toml::find_or<bool>(l, "smooth_nucleation", config.loading.smooth_nucleation);
      }

      if (data.contains("fault_geometry"))
      {
         const auto &g = data.at("fault_geometry");
         config.fault_geom.Wf = toml::find_or<double>(g, "Wf", config.fault_geom.Wf);
         config.fault_geom.lf = toml::find_or<double>(g, "lf", config.fault_geom.lf);
         config.fault_geom.hs = toml::find_or<double>(g, "hs", config.fault_geom.hs);
         config.fault_geom.ht = toml::find_or<double>(g, "ht", config.fault_geom.ht);
         config.fault_geom.H = toml::find_or<double>(g, "H", config.fault_geom.H);
         config.fault_geom.l_vw = toml::find_or<double>(g, "l_vw", config.fault_geom.l_vw);
         config.fault_geom.w_nuc = toml::find_or<double>(g, "w_nuc", config.fault_geom.w_nuc);
      }

      if (data.contains("boundary"))
      {
         const auto &bd = data.at("boundary");
         if (bd.contains("dirichlet"))
         {
            config.boundary.dirichlet_attrs.clear();
            for (const auto &v : toml::find<std::vector<toml::integer>>(bd, "dirichlet"))
            {
               config.boundary.dirichlet_attrs.insert(v);
            }
         }
         if (bd.contains("natural"))
         {
            config.boundary.natural_attrs.clear();
            for (const auto &v : toml::find<std::vector<toml::integer>>(bd, "natural"))
            {
               config.boundary.natural_attrs.insert(v);
            }
         }
         config.boundary.fault_attr = static_cast<int>(toml::find_or<toml::integer>(bd, "fault", config.boundary.fault_attr));
      }

      if (data.contains("solver"))
      {
         const auto &s = data.at("solver");
         config.solver.dg_method = toml::find_or<std::string>(s, "dg_method", config.solver.dg_method);
         config.solver.solver_type = toml::find_or<std::string>(s, "solver_type", config.solver.solver_type);
         config.solver.penalty_factor = toml::find_or<double>(s, "penalty_factor", config.solver.penalty_factor);
         config.solver.blr_tol = toml::find_or<double>(s, "blr_tol", config.solver.blr_tol);
         config.solver.check_residual = toml::find_or<bool>(s, "check_residual", config.solver.check_residual);
         config.solver.face_basis_type = static_cast<int>(toml::find_or<toml::integer>(s, "face_basis_type", config.solver.face_basis_type));
      }

      if (data.contains("time"))
      {
         const auto &t = data.at("time");
         config.time.t_final = toml::find_or<double>(t, "t_final", config.time.t_final);
         config.time.atol = toml::find_or<double>(t, "atol", config.time.atol);
         config.time.rtol = toml::find_or<double>(t, "rtol", config.time.rtol);
         config.time.max_steps = static_cast<int>(toml::find_or<toml::integer>(t, "max_steps", config.time.max_steps));
         config.time.checkpoint_interval = static_cast<int>(toml::find_or<toml::integer>(t, "checkpoint_interval", config.time.checkpoint_interval));
         config.time.tandem_time_stepping = toml::find_or<bool>(t, "tandem_time_stepping", config.time.tandem_time_stepping);
         config.time.use_petsc_ts = toml::find_or<bool>(t, "use_petsc_ts", config.time.use_petsc_ts);
         config.time.petsc_ts_options = toml::find_or<std::string>(t, "petsc_ts_options", config.time.petsc_ts_options);
      }

      if (data.contains("output"))
      {
         const auto &o = data.at("output");
         config.output.output_dir = toml::find_or<std::string>(o, "output_dir", config.output.output_dir);
         config.output.output_prefix = toml::find_or<std::string>(o, "output_prefix", config.output.output_prefix);
         config.output.write_every_step = toml::find_or<bool>(o, "write_every_step", config.output.write_every_step);
         config.output.ref_dir = toml::find_or<std::string>(o, "ref_dir", config.output.ref_dir);
         config.output.regression_tolerance = toml::find_or<double>(o, "regression_tolerance", config.output.regression_tolerance);
      }

      if (data.contains("simulation"))
      {
         const auto &s = data.at("simulation");
         config.simulation.mode = toml::find_or<std::string>(s, "mode", config.simulation.mode);
      }

      // [stress] — SAFS sidecar configuration (Phase 6 §7 of
      // PLAN_onfaultstress.md).  Strict opt-in.  Missing section
      // leaves use_sidecar = false, which preserves every existing
      // BP5 / BP2 / TPV* path bit-exact.
      if (data.contains("stress"))
      {
         const auto &st = data.at("stress");
         config.stress.use_sidecar = toml::find_or<bool>(
            st, "use_sidecar", config.stress.use_sidecar);
         config.stress.sidecar_path = toml::find_or<std::string>(
            st, "sidecar_path", config.stress.sidecar_path);
         config.stress.P_p_pa = toml::find_or<double>(
            st, "P_p_pa", config.stress.P_p_pa);
         config.stress.P_p_grad_pa_per_m = toml::find_or<double>(
            st, "P_p_grad_pa_per_m", config.stress.P_p_grad_pa_per_m);
         config.stress.min_sigma_n_pa = toml::find_or<double>(
            st, "min_sigma_n_pa", config.stress.min_sigma_n_pa);

         if (config.stress.use_sidecar &&
             config.stress.sidecar_path.empty())
         {
            MFEM_ABORT("seas_config_parser: stress.use_sidecar = true "
                       "requires stress.sidecar_path to be set");
         }
      }

      return config;
   }

   static void ApplyBP5Defaults(SEASConfig &config)
   {
      BP5Params bp5;
      config.material.density = bp5.rho;
      config.material.cs = bp5.cs;
      config.material.nu = bp5.nu;
      config.friction.V0 = bp5.V0;
      config.friction.f0 = bp5.f0;
      config.friction.b = bp5.b;
      config.friction.L0 = bp5.L0;
      config.friction.L_nuc = bp5.L_nuc;
      config.friction.a0 = bp5.a0;
      config.friction.amax = bp5.amax;
      config.friction.sigma_n = bp5.sigma_n;
      config.loading.Vp = bp5.Vp;
      config.loading.V_init = bp5.V_init;
      config.loading.V_nuc = bp5.V_nuc;
      config.loading.delta_tau_factor = bp5.delta_tau_factor;
      config.loading.nucleation_eps = bp5.nucleation_eps;
      config.fault_geom.Wf = bp5.Wf;
      config.fault_geom.lf = bp5.lf;
      config.fault_geom.hs = bp5.hs;
      config.fault_geom.ht = bp5.ht;
      config.fault_geom.H = bp5.H;
      config.fault_geom.l_vw = bp5.l_vw;
      config.fault_geom.w_nuc = bp5.w_nuc;
      config.boundary.dirichlet_attrs = {5};
      config.boundary.natural_attrs = {1};
      config.boundary.fault_attr = 3;
      config.time.t_final = bp5.t_final;
   }
#endif // SEAS_USE_TOML
};

} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_CONFIG_PARSER_HPP

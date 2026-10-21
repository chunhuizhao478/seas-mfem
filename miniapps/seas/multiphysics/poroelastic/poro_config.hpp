// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// poro_config.hpp — the poroelastic option (poroelastic plan rev 3, §4.2,
// step P0.1; equations in poro-derivation-rev1.2).
//
//   [multiphysics.poroelastic]          absent = off; every key required,
//   enabled = true | false              unknown keys abort
//   porosity, permeability_m2, viscosity_pa_s, tortuosity, rho_fluid,
//   rho_solid, bulk_fluid_pa, bulk_solid_pa, fault_hydraulics,
//   background_p_pa, fluid_penalty, drag_gamma, velocity_coupling
//
//   [multiphysics.poroelastic.initial_pressure]   optional (SeisSol
//   origin_m, magnitude_pa, width_per_m2           'PressureInjection')
//
//   [multiphysics.poroelastic.source]   optional; injection rate Q(x, t) (F1)
//   shape = "gaussian" | "well", origin_m, width_per_m2 (gaussian only),
//   rate, time_function = "step" | "ramp" | "smooth_step",
//   rise_time_s (ramp and smooth_step only)
//
// The plan's layout ([multiphysics] poroelastic = true next to a
// [multiphysics.poroelastic] table) is not valid TOML: one key cannot be both
// a boolean and a table.  The switch is therefore the key `enabled` inside
// the table, as [rheology.<model>] enabled in the rheology plan rev 1.
//
// The command line --poroelastic 0|1 overrides the deck key; 1 without the
// table is an error.  A table present while the option is off is ignored with
// a warning.  With the option on, [material_constant_fallback] lambda, mu are
// the drained frame moduli and its rho must equal the bulk density
// (1 - phi) rho_s + phi rho_f to 1e-12 relative (CheckBulkDensity).
//
// ParsePoroConfig never aborts: it returns the first error as text, so the
// unit test (U0-4) can check every case in one process; the driver aborts
// with that text.

#ifndef MFEM_SEAS_PORO_CONFIG_HPP
#define MFEM_SEAS_PORO_CONFIG_HPP

#include "mfem.hpp"

#include <toml.hpp>

#include <string>

namespace mfem
{
namespace seas
{
namespace poro
{

/// SeisSol's starting pressure p_init = magnitude exp(-width |x - origin|^2),
/// with zero strain and zero total stress (derivation §4.2, §7.4).
struct PoroInitialPressure
{
   bool on = false;
   real_t origin_m[3] = {0.0, 0.0, 0.0};
   real_t magnitude_pa = 0.0;
   real_t width_per_m2 = 0.0;
};

/// The injection rate Q(x, t) = rate g(x) f(t) (derivation (4.3), F1).
struct PoroSource
{
   bool on = false;
   std::string shape;            ///< "gaussian" or "well"
   real_t origin_m[3] = {0.0, 0.0, 0.0};
   real_t width_per_m2 = 0.0;    ///< gaussian only
   real_t rate = 0.0;            ///< peak Q [1/s] (gaussian) or Q_w [m^3/s] (well)
   std::string time_function;    ///< "step", "ramp" or "smooth_step"
   real_t rise_time_s = 0.0;     ///< ramp and smooth_step only
};

struct PoroConfig
{
   bool enabled = false;
   // [multiphysics.poroelastic], SI units
   real_t porosity = 0.0;
   real_t permeability_m2 = 0.0;
   real_t viscosity_pa_s = 0.0;
   real_t tortuosity = 0.0;
   real_t rho_fluid = 0.0;
   real_t rho_solid = 0.0;
   real_t bulk_fluid_pa = 0.0;
   real_t bulk_solid_pa = 0.0;
   std::string fault_hydraulics;    ///< "permeable" ("sealed" arrives in P6)
   real_t background_p_pa = 0.0;    ///< p_0, uniform; enters sigma_n' only
   real_t fluid_penalty = 0.0;      ///< kappa_w of derivation (6.5)
   std::string drag_gamma;          ///< "fitted" or "half" (test U2-15 only)
   std::string velocity_coupling;   ///< "momentum" (Gamma_mom, F3) or "diagonal"
   PoroInitialPressure initial_pressure;
   PoroSource source;
   /// Set when the table is present but the option is off: the driver logs it.
   std::string warning;

   /// Bulk density (1 - phi) rho_s + phi rho_f.
   real_t BulkDensity() const
   {
      return (1.0 - porosity) * rho_solid + porosity * rho_fluid;
   }
};

/// Result of ParsePoroConfig: error is empty on success.
struct PoroParseResult
{
   PoroConfig cfg;
   std::string error;
};

/// Parse [multiphysics] and [multiphysics.poroelastic] from the deck root.
/// cli_flag: -1 = --poroelastic not given, 0 or 1 = its value.
PoroParseResult ParsePoroConfig(const toml::value &root, int cli_flag);

/// Empty when rho_fallback equals cfg.BulkDensity() to 1e-12 relative,
/// otherwise the error text.
std::string CheckBulkDensity(const PoroConfig &cfg, real_t rho_fallback);

} // namespace poro
} // namespace seas
} // namespace mfem

#endif // MFEM_SEAS_PORO_CONFIG_HPP

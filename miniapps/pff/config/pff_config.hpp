// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// pff_config.hpp — case configuration of the dynamic phase-field fracture
// miniapp, read from TOML (toml11, as miniapps/seas/config/seas_config_parser.hpp).
//
// One TOML file describes one MOOSE case of farms branch cdms; the case files
// live in <model>/configs/ and cite the MOOSE input they reproduce.  Tables:
//
//   [case]         name, model = "solid" | "hydromech", source (MOOSE input path)
//   [mesh]         file, dim, uniform_refinements, corner point
//   [boundaries]   attribute lists: pressure, absorbing, fixed_x/y/z, confinement,
//                  drained (hydromech: prescribed pore pressure)
//   [material]     E, nu, density
//   [phase_field]  model, Gc or ft, l, eta, p, decomposition, driving,
//                  initial_damage boxes
//   [pulse]        ElkPulseLoadExperiment parameters
//   [confinement]  static outer pressure; initial_state = "none" (applied from
//                  rest at t = 0) | "static" (start from the static equilibrium
//                  under it, MOOSE static_solve.i + SolutionIC)
//   [damping]      zeta_pulse, zeta_quench, t_free, quench ("stiffness" | "mass" | "implicit"),
//                  alpha_mass
//   [time]         end_time or num_pulses, cfl, dt_max, dt_quiet, mass ("consistent"
//                  default, MOOSE ADInertialForce; or "lumped")
//   [coupling]     scheme, max_passes, tolerance
//   [poro]         Biot/porosity/permeability parameters (hydromech only),
//                  drained_pressure = "zero" | "pulse" (value on [boundaries]
//                  drained), preconditioner = "amg" | "jacobi"
//   [output]       directory, interval, fields
//   [metrics]      d_threshold, r_core, r_hole
//
// Keys are checked: an unknown key or table is an error, so a misspelt
// parameter cannot silently fall back to its default.

#ifndef MFEM_PFF_CONFIG_HPP
#define MFEM_PFF_CONFIG_HPP

#include "mfem.hpp"
#include "../common/pff_types.hpp"
#include "../phasefield/phase_field_functions.hpp"
#include "../multiphysics/poroelastic/poro_damage.hpp"
#include "../loading/pulse_load.hpp"

#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace mfem
{
namespace pff
{

enum class CaseModel { Solid, Hydromech };
enum class QuenchMode { Stiffness, Mass, Implicit };
enum class CouplingScheme { SinglePass, FixedPoint };

/// Axis-aligned box with a lower bound on the damage (MOOSE initial_damage_aux
/// on a SubdomainBoundingBoxGenerator block).
struct DamageBox
{
   std::array<real_t, 3> lo{{0.0, 0.0, 0.0}};
   std::array<real_t, 3> hi{{0.0, 0.0, 0.0}};
   real_t value = 0.0;
};

struct CaseConfig
{
   // [case]
   std::string name;
   CaseModel model = CaseModel::Solid;
   std::string source;

   // [mesh]
   std::string mesh_file;
   int dim = 3;
   int uniform_refinements = 0;
   std::array<real_t, 3> corner{{0.0, 0.0, 0.0}};

   // [boundaries]: Gmsh physical-group names (as MOOSE uses them) or numeric
   // attributes written as strings; resolved against the mesh file's
   // $PhysicalNames by ResolveAttributes (io/gmsh_physical_names.hpp).
   std::vector<std::string> pressure_bdr, absorbing_bdr, fixed_z_bdr, confinement_bdr,
       drained_bdr, fixed_x_bdr, fixed_y_bdr;

   // [material]
   real_t E = 0.0, nu = 0.0, density = 0.0;

   // [phase_field]
   PhaseFieldModel pf;
   real_t ft = 0.0;                    ///< used when Gc is derived from the strength
   Decomposition decomposition = Decomposition::Spectral;
   std::vector<DamageBox> initial_damage;

   // [pulse]
   PulseLoadParams pulse;

   // [confinement]
   real_t confinement_pressure = 0.0;
   bool static_initial_state = false;   ///< [confinement] initial_state = "static"

   // [damping]
   real_t zeta_pulse = 0.0, zeta_quench = 0.0, t_free = 0.0, alpha_mass = 0.0;
   QuenchMode quench = QuenchMode::Stiffness;

   // [time]
   real_t end_time = 0.0, cfl = 0.4, dt_max = 0.0, dt_quiet = 0.0;
   bool consistent_mass = true;   ///< [time] mass (Andy 2026-11-22: consistent for production)

   // [coupling]
   CouplingScheme coupling = CouplingScheme::SinglePass;
   int max_passes = 1;
   real_t coupling_tol = 1e-6;

   // [poro]
   PoroDamageParams poro;
   real_t fluid_density = 1000.0;
   real_t initial_pore_pressure = 0.0;
   bool drained_pulse = false;          ///< [poro] drained_pressure = "pulse" (FunctionDirichletBC pp = pulse)
   bool poro_amg = true;                ///< [poro] preconditioner

   // [output]
   std::string output_dir = "pff_out";
   int output_interval = 100;
   std::vector<std::string> output_fields{"d", "vel_mag"};

   // [metrics]
   real_t d_threshold = 0.5, r_core = 0.0, r_hole = 0.0;

   // derived
   real_t K() const { return E / 3.0 / (1.0 - 2.0 * nu); }
   real_t G() const { return E / 2.0 / (1.0 + nu); }
   real_t Cp() const { return std::sqrt((K() + 4.0 / 3.0 * G()) / density); }
   real_t Cs() const { return std::sqrt(G() / density); }
};

} // namespace pff
} // namespace mfem

#endif

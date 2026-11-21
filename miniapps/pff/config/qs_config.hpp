// Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
//
// qs_config.hpp — configuration of a quasi-static staggered phase-field case
// (phase R: the RACCOON tutorials).  One TOML per RACCOON input in
// raccoon/<case>/configs/; see qs_config_parser.hpp for the tables.

#ifndef MFEM_PFF_QS_CONFIG_HPP
#define MFEM_PFF_QS_CONFIG_HPP

#include "mfem.hpp"
#include "../common/pff_types.hpp"
#include "../phasefield/phase_field_functions.hpp"

#include <array>
#include <string>
#include <vector>

namespace mfem
{
namespace pff
{

/// Initial refinement of the elements whose vertex average lies in a box
/// (MOOSE BoxMarker / OrientedBoxMarker with `initial_steps`).
struct RefineBox
{
   bool oriented = false;
   std::array<real_t, 3> lo{{0.0, 0.0, 0.0}}, hi{{0.0, 0.0, 0.0}};    // axis-aligned
   std::array<real_t, 3> center{{0.0, 0.0, 0.0}};                       // oriented
   real_t length = 0.0, width = 0.0, height = 0.0;
   std::array<real_t, 3> length_direction{{1.0, 0.0, 0.0}};
   std::array<real_t, 3> width_direction{{0.0, 1.0, 0.0}};
};

/// A node set made from the boundary sides of `boundary` whose vertices all
/// lie in [lo, hi] (BoundingBoxNodeSetGenerator with
/// construct_side_list_from_node_list).  It becomes a boundary attribute.
struct NodeSetBox
{
   std::string name, boundary;
   std::array<real_t, 3> lo{{0.0, 0.0, 0.0}}, hi{{0.0, 0.0, 0.0}};
};

/// u_component = value + rate * t on the named boundaries (DirichletBC,
/// FunctionDirichletBC with function = 'rate*t').
struct DirichletSpec
{
   std::vector<std::string> boundaries;
   int component = 0;
   real_t value = 0.0;
   real_t rate = 0.0;
};

struct QsCaseConfig
{
   std::string name, source;

   // [mesh]: "cartesian" (GeneratedMeshGenerator) or "slit" (two halves
   // [xmin,xmax] x [ymin,0] and [xmin,xmax] x [0,ymax] stitched only where
   // x >= slit_x, as the mode II StitchedMeshGenerator input).
   std::string mesh_type = "cartesian";
   int nx = 1, ny = 1;
   real_t xmin = 0.0, xmax = 1.0, ymin = 0.0, ymax = 1.0;
   real_t slit_x = 0.0;
   int refine_steps = 0;
   std::vector<RefineBox> refine;
   std::vector<NodeSetBox> nodesets;

   // [material]
   real_t E = 0.0, nu = 0.0;
   real_t K() const { return E / 3.0 / (1.0 - 2.0 * nu); }
   real_t G() const { return E / 2.0 / (1.0 + nu); }

   // [phase_field]
   PhaseFieldModel pf;
   Decomposition decomposition = Decomposition::None;
   /// RACCOON transfers psie_active as a CONSTANT MONOMIAL aux variable: the
   /// damage source sees its element average.  "quadrature" uses the
   /// quadrature-point values instead.
   bool driving_element_average = true;

   std::vector<DirichletSpec> dirichlet;

   // [reaction]: NodalSum of the save_in residual of one component.
   std::string reaction_name = "F";
   std::vector<std::string> reaction_boundaries;
   int reaction_component = 0;

   // [time]
   real_t dt = 0.0, end_time = 0.0;
   int num_steps = 0;   ///< 0: run to end_time

   // [solver]
   int fixed_point_max_its = 20;
   real_t fixed_point_rel_tol = 1e-8, fixed_point_abs_tol = 1e-10;
   bool accept_on_max_fixed_point_iteration = true;
   bool automatic_scaling = false;   ///< MOOSE Executioner automatic_scaling (default false)
   real_t nl_rel_tol = 1e-8, nl_abs_tol = 1e-10;   ///< MOOSE Executioner values
   real_t linear_rel_tol = 1e-12;
   int quadrature_order = 3;  ///< 2x2 Gauss on quads, as MOOSE's SECOND for P1

   // [output]
   std::string output_dir = "out";
   std::vector<int> nodal_output_steps;   ///< steps whose nodal fields are written
};

} // namespace pff
} // namespace mfem

#endif

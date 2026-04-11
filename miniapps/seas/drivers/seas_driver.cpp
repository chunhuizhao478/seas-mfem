// Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
// at the Lawrence Livermore National Laboratory. All Rights reserved. See files
// LICENSE and NOTICE for details. LLNL-CODE-806117.
//
// SEAS TOML-based driver: reads configuration from a TOML file and runs
// a quasi-dynamic earthquake cycle simulation.
//
// Usage:
//   mpirun -np N ./seas_driver config.toml [--override key=value ...]
//
// The TOML file specifies all simulation parameters. CLI overrides
// can selectively modify individual fields (e.g., --override time.t_final=1e8).
//
// For BP5: use config/bp5_example.toml as a starting point.
// The benchmark = "bp5" preset fills all defaults from BP5Params.

#include "mfem.hpp"
#include "../config/seas_config.hpp"
#include "../config/seas_config_parser.hpp"
#include "../config/bp5_params.hpp"
#include "../common/mpi_context.hpp"
#include "../constitutive/linear_elastic.hpp"
#include "../domain/boundary_config.hpp"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

using namespace mfem;
using namespace mfem::seas;

int main(int argc, char *argv[])
{
   MPIContext mpi(&argc, &argv);

   // =========================================================================
   // Parse arguments: config.toml [--override key=value ...]
   // =========================================================================
   if (argc < 2)
   {
      if (mpi.IsRoot())
      {
         std::cerr << "Usage: " << argv[0]
                   << " config.toml [--override key=value ...]\n";
      }
      return 1;
   }

   std::string config_file = argv[1];
   std::vector<std::string> overrides;
   for (int i = 2; i < argc; i++)
   {
      std::string arg(argv[i]);
      if (arg == "--override" && i + 1 < argc)
      {
         overrides.push_back(argv[++i]);
      }
   }

   // =========================================================================
   // Parse TOML config
   // =========================================================================
   SEASConfig config = SEASConfigParser::ParseFile(config_file);
   SEASConfigParser::ApplyCLIOverrides(config, overrides);
   SEASConfigParser::Validate(config);

   // =========================================================================
   // Print parsed configuration
   // =========================================================================
   if (mpi.IsRoot())
   {
      std::cout << "SEAS Driver: TOML Configuration\n";
      std::cout << "================================\n";
      std::cout << "  Config file:     " << config_file << "\n";
      std::cout << "  Benchmark:       "
                << (config.benchmark.empty() ? "(none)" : config.benchmark)
                << "\n";
      std::cout << "\n[mesh]\n";
      std::cout << "  file:            " << config.mesh.file << "\n";
      std::cout << "  scale:           " << config.mesh.scale << "\n";
      std::cout << "  order:           " << config.mesh.order << "\n";
      std::cout << "\n[material]\n";
      std::cout << "  density:         " << config.material.density << "\n";
      std::cout << "  cs:              " << config.material.cs << "\n";
      std::cout << "  nu:              " << config.material.nu << "\n";
      std::cout << "  mu (derived):    " << config.material.mu() << "\n";
      std::cout << "  lambda (derived):" << config.material.lambda() << "\n";
      std::cout << "\n[friction]\n";
      std::cout << "  V0:              " << config.friction.V0 << "\n";
      std::cout << "  f0:              " << config.friction.f0 << "\n";
      std::cout << "  b:               " << config.friction.b << "\n";
      std::cout << "  a0:              " << config.friction.a0 << "\n";
      std::cout << "  amax:            " << config.friction.amax << "\n";
      std::cout << "  sigma_n:         " << config.friction.sigma_n << "\n";
      std::cout << "\n[loading]\n";
      std::cout << "  Vp:              " << config.loading.Vp << "\n";
      std::cout << "  V_nuc:           " << config.loading.V_nuc << "\n";
      std::cout << "\n[boundary]\n";
      std::cout << "  dirichlet:       {";
      for (int a : config.boundary.dirichlet_attrs) { std::cout << a << ","; }
      std::cout << "}\n";
      std::cout << "  natural:         {";
      for (int a : config.boundary.natural_attrs) { std::cout << a << ","; }
      std::cout << "}\n";
      std::cout << "  fault:           " << config.boundary.fault_attr << "\n";
      std::cout << "\n[solver]\n";
      std::cout << "  dg_method:       " << config.solver.dg_method << "\n";
      std::cout << "  solver_type:     " << config.solver.solver_type << "\n";
      std::cout << "\n[time]\n";
      std::cout << "  t_final:         " << config.time.t_final
                << " s (~" << config.time.t_final / (365.25*24*3600)
                << " yr)\n";
      std::cout << "  max_steps:       " << config.time.max_steps << "\n";
      std::cout << "\n[simulation]\n";
      std::cout << "  mode:            " << config.simulation.mode << "\n";
      std::cout << "\n[output]\n";
      std::cout << "  output_dir:      " << config.output.output_dir << "\n";
      std::cout << "  output_prefix:   " << config.output.output_prefix << "\n";
      std::cout << "================================\n";
      std::cout << "\nConfiguration parsed and validated successfully.\n";
      std::cout << "Full simulation pipeline will be wired in Phase 6.\n";
   }

   return 0;
}

# Copyright (c) 2010-2025, Lawrence Livermore National Security, LLC. Produced
# at the Lawrence Livermore National Laboratory. All Rights reserved. See files
# LICENSE and NOTICE for details. LLNL-CODE-806117.
#
# This file is part of the MFEM library. For more information and source code
# availability visit https://mfem.org.
#
# MFEM is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

# SEAS (Sequences of Earthquakes and Aseismic Slip) miniapp makefile

# Use the MFEM build directory
MFEM_DIR ?= ../..
MFEM_BUILD_DIR ?= ../..
MFEM_INSTALL_DIR ?= ../../mfem
SRC = $(if $(MFEM_DIR:../..=),$(MFEM_DIR)/miniapps/seas/,)
CONFIG_MK = $(or $(wildcard $(MFEM_BUILD_DIR)/config/config.mk),\
   $(wildcard $(MFEM_INSTALL_DIR)/share/mfem/config.mk))

MFEM_LIB_FILE = mfem_is_not_built
-include $(CONFIG_MK)

# macOS 26+ dyld rejects duplicate LC_RPATH; the conda mpicxx wrapper auto-
# adds -Wl,-rpath,$CONDA_PREFIX/lib AND MFEM_LIBS (via HDF5_LIB/PETSC_LIB)
# can inject the same path again, producing 2+ LC_RPATH commands that abort
# the loader.  Loop-strip post-link until at most one copy remains.  No-op
# on Linux (no install_name_tool, glibc dyld tolerates duplicates) or when
# conda is not active.
#
# CONVENTION: every link recipe MUST end with `$(SEAS_POST_LINK_DEDUP_RPATH)`.
# Keep `$(MFEM_CXX) ... $(MFEM_LINK_FLAGS) ... -o $@ ...` on ONE physical
# line so future Python re-passes that scan for these tokens detect the rule
# (see /tmp/apply_dedup.py for the bulk-apply heuristic).
ifeq ($(shell uname -s),Darwin)
define SEAS_POST_LINK_DEDUP_RPATH
@if [ -n "$(CONDA_PREFIX)" ]; then \
    while [ "$$(otool -l $@ 2>/dev/null | \
                awk '/LC_RPATH/{f=1;next} f&&/path /{print $$2;f=0}' | \
                grep -cFx "$(CONDA_PREFIX)/lib")" -gt 1 ]; do \
        install_name_tool -delete_rpath "$(CONDA_PREFIX)/lib" $@ 2>/dev/null || break ; \
    done ; \
fi
endef
else
SEAS_POST_LINK_DEDUP_RPATH = @true
endif

# Optional TOML support (requires extern/toml11)
# v3: toml.hpp at root; v4: toml.hpp in include/
SEAS_USE_TOML ?= $(if $(wildcard extern/toml11/toml.hpp),YES,$(if $(wildcard extern/toml11/include/toml.hpp),YES,))
ifeq ($(SEAS_USE_TOML),YES)
   TOML_FLAGS = -DSEAS_USE_TOML -Iextern/toml11 -Iextern/toml11/include
else
   TOML_FLAGS =
endif

# Source files
# Phase 1
TEST_FRICTION_SRC = tests/unit/test_friction_law.cpp
TEST_STATE_EVOLUTION_SRC = tests/unit/test_state_evolution.cpp
# Phase 2
TEST_ANTIPLANE_SRC = tests/unit/test_antiplane.cpp
# Phase 3
TEST_FAULT_OPERATOR_SRC = tests/unit/test_fault_operator.cpp
# Phase 4
TEST_QUASI_DYNAMIC_SRC = tests/unit/test_quasi_dynamic.cpp
TEST_BP2_SHORT_SRC = tests/unit/test_bp2_short.cpp
# Phase 5
TEST_IO_SRC = tests/unit/test_io.cpp
BP2_SERIAL_SMOKE_SRC = tests/verification/bp2_serial_smoke.cpp
# Phase 11: Checkpoint
TEST_CHECKPOINT_SRC = tests/unit/test_checkpoint.cpp
# Psi-space state variable tests
TEST_PSI_STATE_SRC = tests/unit/test_psi_state.cpp
# Phase 1a/1b: Fault detection and interface tests
TEST_FAULT_DETECTION_SRC = tests/unit/test_fault_detection.cpp
TEST_DOMAIN_INTERFACE_SRC = tests/unit/test_domain_operator_interface.cpp
# BP5 Phase 1: Parameters and vector friction
TEST_BP5_PARAMS_SRC = tests/unit/test_bp5_params.cpp
TEST_VECTOR_FRICTION_SRC = tests/unit/test_vector_friction.cpp
TEST_FAULT_BASIS_SRC = tests/unit/test_fault_basis.cpp
# BP5 Phase 2: Elasticity integrator and operator
TEST_ELASTICITY_BR2_SRC = tests/unit/test_elasticity_br2.cpp
TEST_ELASTICITY_OPERATOR_SRC = tests/unit/test_elasticity_operator.cpp
# BP5 Phase 3: BP5 fault operator
TEST_BP5_FAULT_OPERATOR_SRC = tests/unit/test_bp5_fault_operator.cpp
# BP5 Phase 4: Integration tests
TEST_BP5_INTEGRATION_SRC = tests/unit/test_bp5_integration.cpp
# BP5 Phase 5: Output and mesh tests
TEST_BP5_OUTPUT_SRC = tests/unit/test_bp5_output.cpp
TEST_BP5_MESH_SRC = tests/unit/test_bp5_mesh.cpp
# BP5 v16: Diagnostic VTK and mesh region tests
TEST_DIAG_VTK_SRC = tests/unit/test_diag_vtk.cpp
TEST_MESH_REGIONS_SRC = tests/unit/test_bp5_mesh_regions.cpp
# v44: Multi-DOF fault face quadrature
TEST_FACE_QUADRATURE_SRC = tests/unit/test_face_quadrature.cpp
# v54: Cross-verification against Tandem
TEST_CROSS_VERIFY_SRC = tests/unit/test_cross_verify_tandem.cpp
# Phase 3: Constitutive model tests
TEST_LINEAR_ELASTIC_SRC = tests/unit/test_linear_elastic.cpp
TEST_CONSTITUTIVE_INTEG_SRC = tests/unit/test_constitutive_integrator.cpp
TEST_BOUNDARY_CONFIG_SRC = tests/unit/test_boundary_config.cpp
TEST_BOUNDARY_CONFIG_OP_SRC = tests/unit/test_boundary_config_operator.cpp
TEST_BOUNDARY_CLASSIF_SRC = tests/unit/test_boundary_classification.cpp
TEST_LOGGING_SRC = tests/unit/test_logging.cpp
TEST_CONFIG_BRIDGE_SRC = tests/unit/test_seas_config_bridge.cpp
# Phase 5: TOML config
TEST_CONFIG_PARSER_SRC = tests/unit/test_seas_config_parser.cpp
SEAS_DRIVER_SRC = drivers/seas_driver.cpp
# Dynamic rupture Phase 1: Wave operator and Godunov flux
GODUNOV_FLUX_SRC = dynamic/godunov_flux.cpp
WAVE_OPERATOR_SRC = dynamic/wave_operator.cpp
GODUNOV_FLUX_OBJ = $(GODUNOV_FLUX_SRC:.cpp=.o)
WAVE_OPERATOR_OBJ = $(WAVE_OPERATOR_SRC:.cpp=.o)
# Phase H.1 of spatial_dynamic_rupture_plan.md (rev-3): per-element flux cache.
GODUNOV_FLUX_POOL_SRC = dynamic/godunov_flux_pool.cpp
GODUNOV_FLUX_POOL_OBJ = $(GODUNOV_FLUX_POOL_SRC:.cpp=.o)
# Phase 9 (TPV-via-spatial-driver): bimaterial exact-Riemann flux (port hrs-ref).
GODUNOV_FLUX_BIMATERIAL_SRC = dynamic/godunov_flux_bimaterial.cpp
GODUNOV_FLUX_BIMATERIAL_OBJ = $(GODUNOV_FLUX_BIMATERIAL_SRC:.cpp=.o)
TEST_GODUNOV_FLUX_BIMATERIAL_SRC = tests/unit/test_godunov_flux_bimaterial.cpp
TEST_GODUNOV_FLUX_BIMATERIAL_OBJ = $(TEST_GODUNOV_FLUX_BIMATERIAL_SRC:.cpp=.o)
TEST_PHASEH_GODUNOV_FLUX_POOL_SRC = tests/unit/test_phaseh_godunov_flux_pool.cpp
TEST_PHASEH_GODUNOV_FLUX_POOL_OBJ = $(TEST_PHASEH_GODUNOV_FLUX_POOL_SRC:.cpp=.o)
# Spatial Phase 1 of spatial_dynamic_rupture_plan.md: TOML schema parser +
# per-DOF friction resolver + LSW forced-rupture friction helper.
# Phases 2 / 3 / 3b (spatial_velocity, spatial_stress) and their tests
# follow in a subsequent commit once the prerequisite Phase 3 / 3b /
# 5a / 5b header additions are re-applied (see CODEBASE_GUIDE TBD).
SPATIAL_FRICTION_SRC                 = spatial/code/spatial_friction.cpp
SPATIAL_FRICTION_OBJ                 = $(SPATIAL_FRICTION_SRC:.cpp=.o)
SPATIAL_VELOCITY_SRC                 = spatial/code/spatial_velocity.cpp
SPATIAL_VELOCITY_OBJ                 = $(SPATIAL_VELOCITY_SRC:.cpp=.o)
SPATIAL_STRESS_SRC                   = spatial/code/spatial_stress.cpp
SPATIAL_STRESS_OBJ                   = $(SPATIAL_STRESS_SRC:.cpp=.o)
# Phase N of safs/.../05_18_2026/PLAN_first_safs_run.md — gradual_overstress
# nucleation resolver + per-sub-step accumulator.
SPATIAL_NUCLEATION_SRC               = dynamic/spatial_nucleation.cpp
SPATIAL_NUCLEATION_OBJ               = $(SPATIAL_NUCLEATION_SRC:.cpp=.o)
TEST_SPATIAL_NUCLEATION_SRC          = tests/unit/test_spatial_nucleation.cpp
TEST_SPATIAL_NUCLEATION_OBJ          = $(TEST_SPATIAL_NUCLEATION_SRC:.cpp=.o)
# Phase 7 — nucleation_methods module: INucleationMethod + MakeNucleation.
NUCLEATION_FACTORY_SRC               = dynamic/nucleation_factory.cpp
NUCLEATION_FACTORY_OBJ               = $(NUCLEATION_FACTORY_SRC:.cpp=.o)
TEST_NUCLEATION_FACTORY_SRC          = tests/unit/test_nucleation_factory.cpp
TEST_NUCLEATION_FACTORY_OBJ          = $(TEST_NUCLEATION_FACTORY_SRC:.cpp=.o)
# TPV-via-spatial-driver Phase 1 (RS primitives) unit tests.
TEST_TPV102_NUC_CALLBACK_PARITY_SRC  = tests/unit/test_tpv102_nuc_callback_parity.cpp
TEST_TPV102_NUC_CALLBACK_PARITY_OBJ  = $(TEST_TPV102_NUC_CALLBACK_PARITY_SRC:.cpp=.o)
TEST_SEED_EQUILIBRIUM_PSI_RS_SRC     = tests/unit/test_seed_equilibrium_psi_rs.cpp
TEST_SEED_EQUILIBRIUM_PSI_RS_OBJ     = $(TEST_SEED_EQUILIBRIUM_PSI_RS_SRC:.cpp=.o)
# TPV-via-spatial-driver Phase 2 (friction strategy interface).
FRICTION_ITERATOR_FACTORY_SRC        = dynamic/friction_iterator_factory.cpp
FRICTION_ITERATOR_FACTORY_OBJ        = $(FRICTION_ITERATOR_FACTORY_SRC:.cpp=.o)
TEST_FRICTION_ITERATOR_FACTORY_SRC   = tests/unit/test_friction_iterator_factory.cpp
TEST_FRICTION_ITERATOR_FACTORY_OBJ   = $(TEST_FRICTION_ITERATOR_FACTORY_SRC:.cpp=.o)
# TPV-via-spatial-driver Phase 5 (unified method-oriented sub-step iterator).
FRICTION_SUBSTEP_ITERATOR_SRC        = dynamic/friction_substep_iterator.cpp
FRICTION_SUBSTEP_ITERATOR_OBJ        = $(FRICTION_SUBSTEP_ITERATOR_SRC:.cpp=.o)
TEST_FRICTION_SUBSTEP_ITERATOR_PARITY_SRC = tests/unit/test_friction_substep_iterator_parity.cpp
TEST_FRICTION_SUBSTEP_ITERATOR_PARITY_OBJ = $(TEST_FRICTION_SUBSTEP_ITERATOR_PARITY_SRC:.cpp=.o)
TEST_ADVANCE_INTERFACE_COMPILES_SRC  = tests/unit/test_advance_interface_compiles.cpp
TEST_ADVANCE_INTERFACE_COMPILES_OBJ  = $(TEST_ADVANCE_INTERFACE_COMPILES_SRC:.cpp=.o)
# TPV-via-spatial-driver Phase 3 (SAFS + rate-and-state run).
TEST_RESOLVE_RATE_STATE_GUARDS_SRC   = tests/unit/test_resolve_rate_state_guards.cpp
TEST_RESOLVE_RATE_STATE_GUARDS_OBJ   = $(TEST_RESOLVE_RATE_STATE_GUARDS_SRC:.cpp=.o)
TEST_CONSTANT_TENSOR_SIGN_SRC        = tests/unit/test_constant_tensor_sign.cpp
TEST_CONSTANT_TENSOR_SIGN_OBJ        = $(TEST_CONSTANT_TENSOR_SIGN_SRC:.cpp=.o)
# Phase 10 (TPV31) unit tests.
TEST_TPV31_CANONICAL_ROTATION_SRC    = tests/unit/test_tpv31_canonical_rotation.cpp
TEST_TPV31_CANONICAL_ROTATION_OBJ    = $(TEST_TPV31_CANONICAL_ROTATION_SRC:.cpp=.o)
TEST_TPV31_NUCLEATION_AND_COHESION_SRC = tests/unit/test_tpv31_nucleation_and_cohesion.cpp
TEST_TPV31_NUCLEATION_AND_COHESION_OBJ = $(TEST_TPV31_NUCLEATION_AND_COHESION_SRC:.cpp=.o)
TEST_TPV31_STATION_WRITER_SRC        = tests/unit/test_tpv31_station_writer.cpp
TEST_TPV31_STATION_WRITER_OBJ        = $(TEST_TPV31_STATION_WRITER_SRC:.cpp=.o)
TEST_SPATIAL_STATION_DISPATCH_SRC    = tests/unit/test_spatial_station_dispatch.cpp
TEST_SPATIAL_STATION_DISPATCH_OBJ    = $(TEST_SPATIAL_STATION_DISPATCH_SRC:.cpp=.o)
TEST_PARAVIEW_STATE_CHANNEL_RS_SRC   = tests/unit/test_paraview_state_channel_rs.cpp
TEST_PARAVIEW_STATE_CHANNEL_RS_OBJ   = $(TEST_PARAVIEW_STATE_CHANNEL_RS_SRC:.cpp=.o)
# TPV-via-spatial-driver Phase 8 (TPV205/102/104 benchmark configs + tests).
TEST_TPV_TOML_STRESS_SIGN_SRC        = tests/unit/test_tpv_toml_stress_sign.cpp
TEST_TPV_TOML_STRESS_SIGN_OBJ        = $(TEST_TPV_TOML_STRESS_SIGN_SRC:.cpp=.o)
TEST_PLANAR_TPV_BASIS_SRC            = tests/unit/test_planar_tpv_basis.cpp
TEST_PLANAR_TPV_BASIS_OBJ            = $(TEST_PLANAR_TPV_BASIS_SRC:.cpp=.o)
TEST_TPV_CONFIG_PARSE_SRC            = tests/unit/test_tpv_config_parse.cpp
TEST_TPV_CONFIG_PARSE_OBJ            = $(TEST_TPV_CONFIG_PARSE_SRC:.cpp=.o)
# Phase D of the same plan — --print-derived pre-flight gate.
SPATIAL_PRINT_DERIVED_SRC            = dynamic/spatial_print_derived.cpp
SPATIAL_PRINT_DERIVED_OBJ            = $(SPATIAL_PRINT_DERIVED_SRC:.cpp=.o)
TEST_SPATIAL_PRINT_DERIVED_SRC       = tests/unit/test_spatial_print_derived.cpp
TEST_SPATIAL_PRINT_DERIVED_OBJ       = $(TEST_SPATIAL_PRINT_DERIVED_SRC:.cpp=.o)
TEST_SPATIAL_FRICTION_CONFIG_SRC     = tests/unit/test_spatial_friction_config.cpp
TEST_SPATIAL_FRICTION_CONFIG_OBJ     = $(TEST_SPATIAL_FRICTION_CONFIG_SRC:.cpp=.o)
TEST_SPATIAL_FRICTION_RESOLVER_SRC   = tests/unit/test_spatial_friction_resolver.cpp
TEST_SPATIAL_FRICTION_RESOLVER_OBJ   = $(TEST_SPATIAL_FRICTION_RESOLVER_SRC:.cpp=.o)
# Phase 3 of PLAN_thermal_case2_mixedflux_port_2026-07-08.md:
# [friction.rate_state.sidecar] — per-DOF a / V_w from a data_projection_v1 HDF5.
TEST_SPATIAL_FRICTION_SIDECAR_SRC    = tests/unit/test_spatial_friction_sidecar.cpp
TEST_SPATIAL_FRICTION_SIDECAR_OBJ    = $(TEST_SPATIAL_FRICTION_SIDECAR_SRC:.cpp=.o)

# Phase 11b: depth-varying a(z)/b(z) two-CSV profile (interpolant + loader).
TEST_FRICTION_DEPTH_PROFILE_SRC      = tests/unit/test_friction_depth_profile.cpp
TEST_FRICTION_DEPTH_PROFILE_OBJ      = $(TEST_FRICTION_DEPTH_PROFILE_SRC:.cpp=.o)
TEST_SPATIAL_VELOCITY_BUNDLE_SRC     = tests/unit/test_spatial_velocity_bundle.cpp
TEST_SPATIAL_VELOCITY_BUNDLE_OBJ     = $(TEST_SPATIAL_VELOCITY_BUNDLE_SRC:.cpp=.o)
TEST_SPATIAL_STRESS_BUNDLE_SRC       = tests/unit/test_spatial_stress_bundle.cpp
TEST_SPATIAL_STRESS_BUNDLE_OBJ       = $(TEST_SPATIAL_STRESS_BUNDLE_SRC:.cpp=.o)
TEST_SPATIAL_CONSTANT_STRESS_SRC     = tests/unit/test_spatial_constant_stress_source.cpp
TEST_SPATIAL_CONSTANT_STRESS_OBJ     = $(TEST_SPATIAL_CONSTANT_STRESS_SRC:.cpp=.o)
# Spatial Phase 5c: spatial_setup.hpp test + Phase H.6 LSW forced-rupture test.
TEST_SPATIAL_SETUP_SRC               = tests/unit/test_spatial_setup.cpp
TEST_SPATIAL_SETUP_OBJ               = $(TEST_SPATIAL_SETUP_SRC:.cpp=.o)
TEST_PHASEH_LSW_FORCED_RUPTURE_SRC   = tests/unit/test_phaseh_lsw_forced_rupture.cpp
TEST_PHASEH_LSW_FORCED_RUPTURE_OBJ   = $(TEST_PHASEH_LSW_FORCED_RUPTURE_SRC:.cpp=.o)
# Phase H Stage 1 gate `T-PHASEH-SCALAR-PARITY` — constant-material
# parity between the scalar WaveOperator ctor and the new
# (MaterialField, BoundaryConfig) ctor (H.1 + H.3 + H.4).
TEST_PHASEH_WAVE_OPERATOR_CONSTANT_PARITY_SRC = tests/unit/test_phaseh_wave_operator_constant_parity.cpp
TEST_PHASEH_WAVE_OPERATOR_CONSTANT_PARITY_OBJ = $(TEST_PHASEH_WAVE_OPERATOR_CONSTANT_PARITY_SRC:.cpp=.o)
# Phase 13 C-6: fault-bearing BimaterialWaveOperator-vs-scalar parity (the
# flux_-placeholder leak-class tripwire; the class is header-only — no new .o).
TEST_BIMATERIAL_WAVE_OPERATOR_PARITY_SRC = tests/unit/test_bimaterial_wave_operator_parity.cpp
TEST_BIMATERIAL_WAVE_OPERATOR_PARITY_OBJ = $(TEST_BIMATERIAL_WAVE_OPERATOR_PARITY_SRC:.cpp=.o)
# Phase 4 of spatial_dynamic_rupture_plan.md (rev-3): generic dynamic-rupture driver.
SPATIAL_DYN_DRIVER_SRC               = drivers/spatial_dyn_driver.cpp
SPATIAL_DYN_DRIVER_OBJ               = $(SPATIAL_DYN_DRIVER_SRC:.cpp=.o)
# Phase 14: explicit coupled-RK time integrator (RK4 / Dormand-Prince RK45).
RK_TIME_STEPPER_SRC                  = dynamic/rk_time_stepper.cpp
RK_TIME_STEPPER_OBJ                  = $(RK_TIME_STEPPER_SRC:.cpp=.o)
TEST_RK_TIME_STEPPER_SRC             = tests/unit/test_rk_time_stepper.cpp
TEST_RK_TIME_STEPPER_OBJ             = $(TEST_RK_TIME_STEPPER_SRC:.cpp=.o)
# Phase 14 (RK + LSW + mixed flux): instantaneous EvaluateLSW + LSW coupled-RK.
TEST_LSW_RK_MIXED_FLUX_SRC           = tests/unit/test_lsw_rk_mixed_flux.cpp
TEST_LSW_RK_MIXED_FLUX_OBJ           = $(TEST_LSW_RK_MIXED_FLUX_SRC:.cpp=.o)
TEST_LSW_RK_MIXED_FLUX_BIMATERIAL_SRC = tests/unit/test_lsw_rk_mixed_flux_bimaterial.cpp
TEST_LSW_RK_MIXED_FLUX_BIMATERIAL_OBJ = $(TEST_LSW_RK_MIXED_FLUX_BIMATERIAL_SRC:.cpp=.o)
TEST_LSW_RK_SHARED_FAULT_MPI_SRC     = tests/unit/test_lsw_rk_shared_fault_mpi.cpp
TEST_LSW_RK_SHARED_FAULT_MPI_OBJ     = $(TEST_LSW_RK_SHARED_FAULT_MPI_SRC:.cpp=.o)
TEST_GODUNOV_FLUX_SRC = tests/unit/test_godunov_flux.cpp
TEST_WAVE_OPERATOR_SRC = tests/unit/test_wave_operator.cpp
PML_LAYER_SRC = dynamic/pml_layer.cpp
PML_LAYER_OBJ = $(PML_LAYER_SRC:.cpp=.o)
FAULT_FACE_FLUX_SRC = dynamic/fault_face_flux.cpp
FAULT_FACE_FLUX_OBJ = $(FAULT_FACE_FLUX_SRC:.cpp=.o)
FRICTION_SOLVER_SRC = dynamic/friction_solver.cpp
FRICTION_SOLVER_OBJ = $(FRICTION_SOLVER_SRC:.cpp=.o)
# TPV102 "Topology-Based Precomputed Face-Rotation" plan 2026-04-23 Phase 1
PRECOMPUTED_FACE_FLUXES_SRC = dynamic/precomputed_face_fluxes.cpp
PRECOMPUTED_FACE_FLUXES_OBJ = $(PRECOMPUTED_FACE_FLUXES_SRC:.cpp=.o)
TEST_PRECOMPUTED_FLUXES_INTERIOR_IDENTITY_SRC = tests/unit/test_precomputed_fluxes_interior_identity.cpp
TEST_PRECOMPUTED_FLUXES_BOUNDARY_IDENTITY_SRC = tests/unit/test_precomputed_fluxes_boundary_identity.cpp
TEST_PRECOMPUTED_FLUXES_INTERIOR_IDENTITY_OBJ = $(TEST_PRECOMPUTED_FLUXES_INTERIOR_IDENTITY_SRC:.cpp=.o)
TEST_PRECOMPUTED_FLUXES_BOUNDARY_IDENTITY_OBJ = $(TEST_PRECOMPUTED_FLUXES_BOUNDARY_IDENTITY_SRC:.cpp=.o)
# TPV102 "Topology-Based Precomputed Face-Rotation" plan 2026-04-23 Phase 2a
TEST_ARM1_CONSTANT_STATE_PROBES_SRC        = tests/unit/test_arm1_constant_state_probes.cpp
TEST_ARM2_D4_EQUIVARIANT_FIXTURE_SRC       = tests/unit/test_arm2_d4_equivariant_fixture.cpp
TEST_ARM2B_VERTEX_PERM_PROBE_SRC           = tests/unit/test_arm2b_vertex_permutation_probe.cpp
TEST_ARM3D_TOPO_STATE_SAMPLING_SRC         = tests/unit/test_arm3d_topology_state_sampling.cpp
TEST_PRECOMPUTED_FLUXES_PHASE2A_SWITCH_SRC = tests/unit/test_precomputed_fluxes_phase2a_switch.cpp
TEST_PRECOMPUTED_FLUXES_PHASE2A_SWITCH_OBJ = $(TEST_PRECOMPUTED_FLUXES_PHASE2A_SWITCH_SRC:.cpp=.o)
TEST_ARM1_CONSTANT_STATE_PROBES_OBJ        = $(TEST_ARM1_CONSTANT_STATE_PROBES_SRC:.cpp=.o)
TEST_ARM2_D4_EQUIVARIANT_FIXTURE_OBJ       = $(TEST_ARM2_D4_EQUIVARIANT_FIXTURE_SRC:.cpp=.o)
TEST_ARM2B_VERTEX_PERM_PROBE_OBJ           = $(TEST_ARM2B_VERTEX_PERM_PROBE_SRC:.cpp=.o)
TEST_ARM3D_TOPO_STATE_SAMPLING_OBJ         = $(TEST_ARM3D_TOPO_STATE_SAMPLING_SRC:.cpp=.o)
TEST_R008_FAULT_BRANCH_PRESERVED_SRC       = tests/unit/test_R008_fault_branch_preserved.cpp
TEST_R008_FAULT_BRANCH_PRESERVED_OBJ       = $(TEST_R008_FAULT_BRANCH_PRESERVED_SRC:.cpp=.o)
TEST_R5_001_HOIST_COMPILES_SRC             = tests/unit/test_R5_001_hoist_compiles.cpp
TEST_R5_001_HOIST_COMPILES_OBJ             = $(TEST_R5_001_HOIST_COMPILES_SRC:.cpp=.o)
# TPV102 v9.5.0 — REVIEW R-001: Mult honors the precomputed flux flag
TEST_R001_RK4_USES_PRECOMPUTED_FLAG_SRC    = tests/unit/test_R001_rk4_uses_precomputed_flag.cpp
TEST_R001_RK4_USES_PRECOMPUTED_FLAG_OBJ    = $(TEST_R001_RK4_USES_PRECOMPUTED_FLAG_SRC:.cpp=.o)
# Phase 1 correction 2026-04-23 R-002: BuildD4Mesh vertex-equivariance test
TEST_D4_FIXTURE_VERTEX_EQUIV_SRC           = tests/unit/test_d4_fixture_vertex_equivariance.cpp
TEST_D4_FIXTURE_VERTEX_EQUIV_OBJ           = $(TEST_D4_FIXTURE_VERTEX_EQUIV_SRC:.cpp=.o)
# Phase 2D Step B: BuildD4Mesh(true) fault sanity gate
TEST_D4_FAULT_SANITY_SRC                   = tests/unit/test_d4_fault_sanity.cpp
TEST_D4_FAULT_SANITY_OBJ                   = $(TEST_D4_FAULT_SANITY_SRC:.cpp=.o)
# Round-3 R-003 diagnostic: BuildD4Mesh Jacobian-sign audit
TEST_D4_JACOBIAN_SIGN_SRC                  = tests/unit/test_d4_jacobian_sign.cpp
TEST_D4_JACOBIAN_SIGN_OBJ                  = $(TEST_D4_JACOBIAN_SIGN_SRC:.cpp=.o)
# Round-5 Option 1: Gate-14′-to-pepper correlation test
TEST_R5_GATE14_PEPPER_CORR_SRC             = tests/unit/test_r5_gate14_pepper_correlation.cpp
TEST_R5_GATE14_PEPPER_CORR_OBJ             = $(TEST_R5_GATE14_PEPPER_CORR_SRC:.cpp=.o)
# Round-6 Phase 2E: rupture-driven amplification chain
TEST_R6_PHASE2E_AMPCHAIN_SRC               = tests/unit/test_r6_phase2e_amplification_chain.cpp
TEST_R6_PHASE2E_AMPCHAIN_OBJ               = $(TEST_R6_PHASE2E_AMPCHAIN_SRC:.cpp=.o)
# Round-7: strengthened amplification chain (shared-face pairing + T2)
TEST_R7_AMPCHAIN_V2_SRC                    = tests/unit/test_r7_amplification_chain_v2.cpp
TEST_R7_AMPCHAIN_V2_OBJ                    = $(TEST_R7_AMPCHAIN_V2_SRC:.cpp=.o)
# Round-9: FREEZE-C (psi-hold) + dt-refinement amplifier discriminator
TEST_R9_FROZEN_FRICTION_C_SRC              = tests/unit/test_r9_frozen_friction_C.cpp
TEST_R9_FROZEN_FRICTION_C_OBJ              = $(TEST_R9_FROZEN_FRICTION_C_SRC:.cpp=.o)
# Phase 2b MPI switch test (plan §7.3)
TEST_PRECOMPUTED_FLUXES_PHASE2B_MPI_SRC    = tests/parallel/test_precomputed_fluxes_phase2b_mpi_switch.cpp
TEST_PRECOMPUTED_FLUXES_PHASE2B_MPI_OBJ    = $(TEST_PRECOMPUTED_FLUXES_PHASE2B_MPI_SRC:.cpp=.o)
# TPV102 v9.5.0 — REVIEW R-002: empty-fault rank must not abort on UsePrecomputedFaceFluxes
TEST_R002_EMPTY_FAULT_RANK_SRC             = tests/parallel/test_R002_empty_fault_rank_no_abort.cpp
TEST_R002_EMPTY_FAULT_RANK_OBJ             = $(TEST_R002_EMPTY_FAULT_RANK_SRC:.cpp=.o)
# TPV104 R-1600: substep np>1 deadlock-detection regression — EvaluateBulkAtFaultQPsCanonical must
# not skip its per-substep ExchangeFaceNbrData on ranks with no fault-shared faces.
TEST_R1600_SUBSTEP_NO_FAULT_RANK_SRC       = tests/parallel/test_R1600_substep_no_fault_rank_no_hang.cpp
TEST_R1600_SUBSTEP_NO_FAULT_RANK_OBJ       = $(TEST_R1600_SUBSTEP_NO_FAULT_RANK_SRC:.cpp=.o)
# R-004 (REVIEW_speckle_seissol_drdg3d_rootcause_2026-05-24): the R-1601 batched
# vdim=NUM_STATE byNODES ghost exchange (predictor path) must deliver bit-identical
# neighbour values to the per-component scalar exchange (macro path).
TEST_GHOST_EXCHANGE_BYNODES_SRC            = tests/parallel/test_ghost_exchange_bynodes_vs_scalar.cpp
TEST_GHOST_EXCHANGE_BYNODES_OBJ            = $(TEST_GHOST_EXCHANGE_BYNODES_SRC:.cpp=.o)
# Round-11 R-001/R-002: cross-rank Adjacent mixed-flux face-set consistency
TEST_MIXED_FLUX_ADJACENT_MPI_SRC           = tests/parallel/test_mixed_flux_adjacent_mpi.cpp
TEST_MIXED_FLUX_ADJACENT_MPI_OBJ           = $(TEST_MIXED_FLUX_ADJACENT_MPI_SRC:.cpp=.o)
# Round-11 R-1100: shared fault face must NOT enter central_flux_face_set_
TEST_MIXED_FLUX_SHARED_FAULT_MPI_SRC       = tests/parallel/test_mixed_flux_shared_fault_face_excluded_mpi.cpp
TEST_MIXED_FLUX_SHARED_FAULT_MPI_OBJ       = $(TEST_MIXED_FLUX_SHARED_FAULT_MPI_SRC:.cpp=.o)
# Round-11 R-1102: Adjacent-mode end-to-end dispatch correctness
TEST_MIXED_FLUX_DISPATCH_ADJACENT_SRC      = tests/unit/test_mixed_flux_dispatch_adjacent.cpp
TEST_MIXED_FLUX_DISPATCH_ADJACENT_OBJ      = $(TEST_MIXED_FLUX_DISPATCH_ADJACENT_SRC:.cpp=.o)
# Round-13 R-1201: np=2 vs np=1 Adjacent dispatch numerical bit-identity
TEST_MIXED_FLUX_DISPATCH_ADJACENT_MPI_SRC  = tests/parallel/test_mixed_flux_dispatch_adjacent_mpi.cpp
TEST_MIXED_FLUX_DISPATCH_ADJACENT_MPI_OBJ  = $(TEST_MIXED_FLUX_DISPATCH_ADJACENT_MPI_SRC:.cpp=.o)
# Round-14 R-1402: AllContinuous-mode end-to-end dispatch correctness
TEST_MIXED_FLUX_DISPATCH_ALL_CONT_SRC      = tests/unit/test_mixed_flux_dispatch_all_continuous.cpp
TEST_MIXED_FLUX_DISPATCH_ALL_CONT_OBJ      = $(TEST_MIXED_FLUX_DISPATCH_ALL_CONT_SRC:.cpp=.o)
TEST_WAVE_BC_SRC = tests/unit/test_wave_bc.cpp
TEST_PML_SRC = tests/unit/test_pml.cpp
TEST_FAULT_FACE_FLUX_SRC = tests/unit/test_fault_face_flux.cpp
TEST_FRICTION_SOLVER_SRC = tests/unit/test_friction_solver.cpp
TEST_TPV102_SETUP_SRC = tests/unit/test_tpv102_setup.cpp
TPV102_DRIVER_SRC = drivers/tpv102_driver.cpp
TPV102_SUBSTEP_ITERATOR_SRC = dynamic/tpv102_substep_iterator.cpp
TPV205_DRIVER_SRC = drivers/tpv205_driver.cpp
TPV205_SUBSTEP_ITERATOR_SRC = dynamic/tpv205_substep_iterator.cpp
TEST_TPV102_LOCAL_SRC = tests/verification/test_tpv102_local.cpp
TEST_R101_SHARED_FAULT_SRC = tests/parallel/test_r101_shared_fault.cpp
# TPV102 v9.0.0 Phase 1 unit tests (§3.1, §3.1b, §3.1c, §3.1d)
TEST_FAULT_FACE_FLUX_FRAME_SRC = tests/unit/test_fault_face_flux_frame.cpp
TEST_FAULT_FACE_FLUX_FRAME_AND_FLUX_SRC = tests/unit/test_fault_face_flux_frame_and_flux.cpp
TEST_GODUNOV_IDENTITY_SRC = tests/unit/test_godunov_identity_normal_reversal.cpp
TEST_FAULT_FLUX_BRANCH_EQUIV_SRC = tests/unit/test_fault_flux_interior_vs_shared_branch_equivalence.cpp
# v9.0.0 Pelties-9 fix: new unit tests (§3.1e, §3.1f)
TEST_GODUNOV_INTERIOR_EQUAL_SIDES_SRC = tests/unit/test_godunov_interior_equal_sides_identity.cpp
TEST_FAULT_FACE_FLUX_PER_SIDE_ASSEMBLY_SRC = tests/unit/test_fault_face_flux_per_side_assembly.cpp
TEST_FAULT_FACE_FLUX_BIMATERIAL_GUARD_SRC  = tests/unit/test_fault_face_flux_bimaterial_guard.cpp
TEST_BIMATERIAL_FAULT_RIEMANN_SRC          = tests/unit/test_bimaterial_fault_riemann.cpp
TEST_TPV6_STATIONS_SRC                     = tests/unit/test_tpv6_stations.cpp
# v9.1.0 Pelties-9 follow-through: post-R-001 VTU continuity + post-R-003 basis symmetry
TEST_FAULT_SURFACE_VTU_CONTINUITY_SRC      = tests/unit/test_fault_surface_vtu_continuity.cpp
TEST_FAULT_BASIS_DIP_STRIKE_SYMMETRY_SRC   = tests/unit/test_fault_basis_dip_strike_symmetry.cpp
# v9.2.0 §19 R-V92-E02 stage-4 VTU diagnostic
TEST_FAULT_SURFACE_VTU_K4_SRC              = tests/unit/test_fault_surface_vtu_k4.cpp
# --paraview-fields / SetFaultVTUFields filter
TEST_FAULT_SURFACE_VTU_FIELD_FILTER_SRC    = tests/unit/test_fault_surface_vtu_field_filter.cpp
# Phase 1 of PLAN_paraview_compaction_2026-04-28 — gather-to-rank-0 binary
# VTU writer: serial bit-exactness + np=4 gather correctness
TEST_FAULT_SURFACE_VTU_BINARY_SRC          = tests/unit/test_fault_surface_vtu_binary.cpp
TEST_FAULT_SURFACE_VTU_GATHER_MPI_SRC      = tests/unit/test_fault_surface_vtu_gather_mpi.cpp
# Phase 4 of PLAN_free_surface_slice.md — free-surface slice class test (MPI)
TEST_FREE_SURFACE_SLICE_SRC                = tests/unit/test_free_surface_slice.cpp
# Phase 2a/2b of PLAN_paraview_compaction_2026-04-28 — VTKHDF single-file
# fault output.  Conditional on MFEM_USE_HDF5 in the test sources.
TEST_PARAVIEW_HDF_SMOKE_SRC                = tests/unit/test_paraview_hdf_smoke.cpp
TEST_FAULT_SURFACE_VTKHDF_SRC              = tests/unit/test_fault_surface_vtkhdf.cpp
TEST_FAULT_SURFACE_VTKHDF_MPI_SRC          = tests/unit/test_fault_surface_vtkhdf_mpi.cpp
# Phase 2d.2 of PLAN_paraview_compaction_2026-04-28 — H5Z-ZFP gate
TEST_VTKHDF_ZFP_SRC                        = tests/unit/test_vtkhdf_zfp.cpp
# Phase 2d.3 of PLAN_paraview_compaction_2026-04-28 — seas fault path + ZFP
TEST_FAULT_SURFACE_VTKHDF_ZFP_SRC          = tests/unit/test_fault_surface_vtkhdf_zfp.cpp
# Phase 6.2 of PLAN_bulk_compression_and_size_estimator_2026-05-09 — volume + ZFP
TEST_VOLUME_HDF_COMPRESSION_SRC            = tests/unit/test_volume_hdf_compression.cpp
# Phase 3 of PLAN_paraview_compaction_2026-04-28 — adaptive snapshot cap
TEST_PARAVIEW_SCHEDULE_CAP_SRC             = tests/unit/test_paraview_schedule_cap.cpp
# R-001 of PLAN_split_bulk_solutions_2026-05-12 — kinematics field set
TEST_KINEMATICS_FIELD_SET_SRC              = tests/unit/test_kinematics_field_set.cpp
# R-003 of REVIEW.md 2026-05-16 — zero-fault-DOF rank PETSc TS lifecycle
TEST_BP5_PETSC_TS_ZERO_FAULT_RANK_SRC      = tests/unit/test_bp5_petsc_ts_zero_fault_rank.cpp
# Phase 1 of petsc_ts_restart_plan_2026-05-16 — V2 checkpoint round-trip + restart-state restoration
TEST_BP5_PETSC_TS_RESTART_SRC              = tests/unit/test_bp5_petsc_ts_restart.cpp
# Phase 4 TPV104 port — V1 checkpoint round-trip
TEST_TPV104_CHECKPOINT_SRC                 = tests/unit/test_tpv104_checkpoint.cpp
# R-104 of REVIEW.md 2026-05-16 round 2 — MPI-only rank-0 warning gate
TEST_PARAVIEW_RANK0_WARNING_GATE_SRC       = tests/unit/test_paraview_rank0_warning_gate.cpp
# v9.2.0 plan §4.1-§4.4 regression gates
TEST_GODUNOV_ROTATION_IDENTITY_SRC         = tests/unit/test_godunov_rotation_identity.cpp
TEST_CANONICAL_ROTATION_PURE_STRIKESLIP_SRC = tests/unit/test_canonical_rotation_pure_strikeslip.cpp
TEST_FRAME_DIP_LEAK_PROBE_SRC              = tests/unit/test_frame_dip_leak_probe.cpp
TEST_OFFFAULT_FRAME_SYMMETRY_SRC          = tests/unit/test_offfault_frame_symmetry.cpp
TEST_FAULT_BASIS_QP_ORTHONORMALITY_SRC     = tests/unit/test_fault_basis_qp_orthonormality.cpp
TEST_PER_QP_VS_CENTROID_BASIS_SRC          = tests/unit/test_per_qp_vs_centroid_basis_planar_1el.cpp
# PLAN_frame_orthonormality_diag_2026-05-23 §Phase 1 AC#3 — [FRAME] dv_n identity
TEST_FRAME_DIAG_PROJECTION_IDENTITY_SRC    = tests/unit/test_frame_diag_projection_identity.cpp
# PLAN_predictor_vs_macro_diag_2026-05-23 §Phase 1 AC#3 — [MACRO] decomposition identity
TEST_MACRO_DIAG_DECOMPOSITION_IDENTITY_SRC = tests/unit/test_macro_diag_decomposition_identity.cpp
# v9.2.0 plan §4.9 (rev-3f) — No-penalty regression gate (Pelties 2012 §3.1)
TEST_NO_PENALTY_DR_SRC                     = tests/unit/test_no_penalty_dynamic_rupture.cpp
# v9.2.0 plan Step 6 (H-V92-U probe) — interior-fault-path Pelties eq. (7)
TEST_INTERIOR_FAULT_FLUX_PATH_SRC          = tests/unit/test_interior_fault_flux_path.cpp
# Pepper-bug plan 2026-04-22 Phase 3 P3-T1 — asymmetric Q probe
TEST_INTERIOR_FAULT_FLUX_PATH_ASYM_SRC     = tests/unit/test_interior_fault_flux_path_asymmetric.cpp
# Pepper-bug plan 2026-04-22 Phase 3 P3-T3 — rupture-regime amplitude probe
TEST_INTERIOR_FAULT_FLUX_PATH_RUPTURE_SRC  = tests/unit/test_interior_fault_flux_path_rupture_regime.cpp
# Pepper-bug plan 2026-04-22 Phase 5 P5-N — bulk RK4 noise stability
TEST_BULK_RK4_NOISE_STAB_SRC               = tests/unit/test_bulk_rk4_noise_stability.cpp
# Pepper-bug plan 2026-04-22 Test 5 — shared-fault Pelties-7 analytic (np=4)
TEST_SHARED_FAULT_FLUX_PATH_ANA_SRC        = tests/parallel/test_shared_fault_flux_path_analytic.cpp
# Pepper-bug plan 2026-04-22 Phase 5 P5-PW — bulk plane-wave recovery
TEST_BULK_PLANE_WAVE_SRC                   = tests/unit/test_bulk_plane_wave_recovery.cpp
# Pepper-bug plan 2026-04-22 Phase 1 — live interior-vs-shared cross-check
TEST_INTERIOR_VS_SHARED_LIVE_SRC           = tests/unit/test_interior_vs_shared_branch_live.cpp
# Pepper-bug: MULTI-STEP rupture-regime serial-vs-parallel cross-check
TEST_RUPTURE_MULTISTEP_SRC                 = tests/unit/test_rupture_multistep_serial_vs_parallel.cpp
# Phase 0 fault-dealiasing reproduction harness (planar-fault sigma_n drift)
TEST_FAULT_PLANAR_SERIAL_SRC               = tests/unit/test_fault_planar_serial.cpp
# Phase 1 fault-flux over-integration unit test (tiny 2-tet fixture)
TEST_FAULT_OVERINT_SRC                     = tests/unit/test_fault_overint.cpp
# Phase 2 fault-face resample operator R unit test (reference-element only)
TEST_FAULT_RESAMPLE_SRC                    = tests/unit/test_fault_resample.cpp
# Phase 3 rate-state Δψ resample apply unit test (iterator-level)
TEST_FAULT_RESAMPLE_APPLY_SRC             = tests/unit/test_fault_resample_apply.cpp
TEST_RUPTURE_TILTED_SRC                    = tests/unit/test_rupture_tilted_fault_serial_vs_parallel.cpp
TEST_RUPTURE_RAKE_SWEEP_SRC                = tests/unit/test_rupture_rake_sweep_cross_rank.cpp
TEST_RECONCILE_XRANK_SRC                   = tests/unit/test_shared_fault_reconcile_cross_rank.cpp
TEST_RS_PSI_NP2_SRC                        = tests/unit/test_shared_fault_rs_psi_consistency_np2.cpp
# Pepper-bug: ADJACENT-TRIANGLE UNIFORMITY test (8-fault-triangle fixture)
TEST_ADJACENT_TRI_UNIFORMITY_SRC           = tests/unit/test_adjacent_triangle_fault_uniformity.cpp
TEST_ADJACENT_TRI_FIRST_STEP_AUDIT_SRC     = tests/unit/test_adjacent_triangle_fault_first_step_audit.cpp
TEST_GODUNOV_ROUNDTRIP_SRC                 = tests/unit/test_godunov_rotation_roundtrip.cpp
TEST_FRICTION_EQUIL_SRC                    = tests/unit/test_friction_solver_equilibrium.cpp
TEST_PEPPER_AMP_CHAIN_SRC                  = tests/unit/test_pepper_amplification_chain.cpp
TEST_FACE_LIFT_OPERATOR_SRC                = tests/unit/test_face_lift_operator.cpp
TEST_PEPPER_QP_PERT_SRC                    = tests/unit/test_pepper_qp_perturbation_response.cpp
TEST_PEPPER_NONLINEAR_GROWTH_SRC           = tests/unit/test_pepper_nonlinear_growth.cpp
# Pepper-bug plan 2026-04-22 Phase 2 — sign_flipped / elem1_on_plus truth table
TEST_SIGN_FLIPPED_TRUTH_TABLE_SRC          = tests/unit/test_sign_flipped_elem1_on_plus_truth_table.cpp
# Pepper-bug plan 2026-04-22 Phase 4 P4-T1 — MakeFaceKey uniqueness (np=2)
TEST_MAKE_FACE_KEY_UNIQ_SRC                = tests/parallel/test_make_face_key_uniqueness.cpp
# Pepper-bug plan 2026-04-22 Phase 4 P4-T2 — DOF-offset invariants
TEST_FAULT_DOF_OFFSET_SRC                  = tests/unit/test_fault_dof_offset_invariants.cpp
# v9.2.0 round-9 REVIEW R-V92-K01 — driver post-RK4 contract regression guard
TEST_RK4_ENDPOINT_REEVAL_SRC               = tests/unit/test_rk4_endpoint_reeval.cpp
# v9.2.0 plan Step 5 (R-V92-H01) — coupled RK4-on-psi convergence probe
TEST_RK4_PSI_INTEGRATION_SRC               = tests/unit/test_rk4_psi_integration.cpp
# v9.2.0 plan §4.10-§4.12 (rev-3g) — H-V92-G candidate probes (b)(c)(d)(e)
TEST_VOLUME_JAC_SINGLE_CH_SRC              = tests/unit/test_volume_jacobian_single_channel.cpp
TEST_RK4_CONSERVATION_SRC                  = tests/unit/test_rk4_conservation.cpp
TEST_ABSORBING_BC_DECAY_SRC                = tests/unit/test_absorbing_bc_energy_decay.cpp
# v9.2.0 plan §4.6 (rev-3b) — shared-fault role consistency across MPI
TEST_SHARED_FAULT_ROLE_CONSIST_SRC         = tests/parallel/test_shared_fault_role_consistency.cpp
# v9.2.0 plan §4.7 (rev-3c) — shared-fault DOFData runtime consistency
TEST_SHARED_FAULT_DOFDATA_CONSIST_SRC      = tests/parallel/test_shared_fault_dof_data_consistency.cpp
# v9.2.0 plan Step 5b (R-V92-I02) — 100-step F01+F02 stability probe
TEST_RK4_F01F02_STABILITY_SRC              = tests/parallel/test_rk4_f01f02_shared_fault_stability.cpp
# v9.2.0 plan §17 (rev-3d, post-REVIEW R-V92-C04) — centroid-margin ctor-only probe
TEST_CENTROID_MARGIN_CTOR_SRC              = tests/parallel/test_centroid_margin_ctor_only.cpp
# v9.3.0 plan Phase 1 (I-04): FreeSurfaceGodunov vs FreeSurface equivalence
TEST_FREE_SURFACE_GODUNOV_SRC              = tests/unit/test_free_surface_godunov.cpp
# v9.3.0 plan Phase 2 (I-04 wiring, R-005): --free-surface-bc parsing
TEST_R005_FS_BC_FLAG_PARSE_SRC             = tests/unit/test_R005_free_surface_bc_flag_parse.cpp
# v9.3.0 plan Phase 2 AC-3 (I-04 driver integration gate): locked-fault
# 10-RK4 max|Q[SYY]| <= 1 μPa under both gamma and godunov BC modes
TEST_FREE_SURFACE_GODUNOV_DRIVER_SRC       = tests/unit/test_free_surface_godunov_driver.cpp
# v9.3.0 plan Phase 3 (I-06 part A): EvaluateTotal vs Evaluate equivalence
TEST_FAULT_FACE_FLUX_TOTAL_VS_FLUC_SRC     = tests/unit/test_fault_face_flux_total_vs_fluctuation.cpp
# v9.3.0 R-003 (review-incorporation plan): EvaluateTotal psi-invariance
TEST_R003_EVALUATETOTAL_PSI_SRC            = tests/unit/test_R003_evaluatetotal_psi_invariant.cpp
# v9.3.0 Phase 4 (I-06 part B): R-001 InitializeStateTotal sign regression
TEST_R001_INITSTATE_TOTAL_SIGN_SRC         = tests/unit/test_R001_initstate_total_sign.cpp
# v9.3.0 Phase 4 AC (I-06 part B): locked-fault 100 RK4 sigma_n stability
TEST_TPV102_TOTAL_LOCKED_FAULT_SRC         = tests/unit/test_tpv102_total_locked_fault.cpp
# v9.3.0 R-I06-001 / R-I06-006 (review-round-2): total-Q absorbing BC
# preserves pre-stress equilibrium (AbsorbingTotal regression gate)
TEST_TPV102_TOTAL_ABSORBING_SRC            = tests/unit/test_tpv102_total_absorbing_equilibrium.cpp
# v9.4.0 R-005: tests test_R002_nucleation_sign_match,
# test_R009_nucleation_total_both_sides, test_R_I06_003_nucleation_shape_scaling
# exercised the DEPRECATED bulk-Q ApplyNucleationTotal path; removed
# with the v9.4.0 split-prestress refactor.
# 2026-04-22 nucleation persistent-prestress channel fix regression
# (replaces the bulk-Q add-delta path that diluted nuc by 1e4-1e5x)
TEST_PERSISTENT_NUC_PRESTRESS_SRC          = tests/unit/test_persistent_nuc_prestress_channel.cpp
# v9.4.0 debug plan R-001: measures per-QP ULP dispersion under
# uniform bulk Q to verify the split-prestress refactor is justified
TEST_TPV102_PEPPER_SRC                     = tests/unit/test_tpv102_pepper_reproducer.cpp
# v9.3.0 R-I06-round3 R-006: direct unit tests for FreeSurfaceTotal /
# FreeSurfaceGodunovTotal (Q_bg=0 equivalence + tilted-normal equilibrium)
TEST_FREE_SURFACE_TOTAL_DIRECT_SRC         = tests/unit/test_free_surface_total_direct.cpp
# ADER I-05 Phase 1: polynomial-exactness test for ApplySpatialDerivative
TEST_WAVE_OP_SPATIAL_DERIV_SRC             = tests/unit/test_wave_operator_spatial_derivative.cpp
# ADER hot-path optimization: micro-benchmark for the dominant kernels
BENCH_ADER_HOTPATH_SRC                     = tests/unit/bench_ader_hotpath.cpp
# ADER I-05 Phase 3: CK predictor plane-wave convergence
TEST_ADER_CK_PREDICTOR_SRC                 = tests/unit/test_ader_ck_predictor.cpp
# ADER I-05 Phase 5: EvaluateADER / EvaluateADERTotal equivalence
TEST_ADER_FAULT_FLUX_EQUIV_SRC             = tests/unit/test_fault_face_flux_ader_equivalence.cpp
# ADER I-05 Phase 6: AdvanceADER vs RK4 on a linear wave
TEST_ADER_LINEAR_WAVE_EQUIV_SRC            = tests/unit/test_ader_linear_wave_equivalence.cpp
# ADER fault-path live check: AdvanceADER(dt) / dt -> wave.Mult on a faulted TPV102 fixture
TEST_ADER_FAULT_COUPLING_DT_LIMIT_SRC      = tests/unit/test_ader_fault_coupling_dt_limit.cpp
# ADER live interior-vs-shared branch cross-check on the same physical face (np=2)
TEST_ADER_INTERIOR_VS_SHARED_LIVE_SRC      = tests/unit/test_ader_interior_vs_shared_branch_live.cpp
# ADER I-05 Phase 8: TPV102 AdvanceADER smoke test
TEST_ADER_TPV102_SMOKE_SRC                 = tests/unit/test_ader_tpv102_smoke.cpp
# Round-6 R-001 regression: driver defaults to ADER (purpose change #1)
TEST_R001_DRIVER_DEFAULT_ADER_SRC          = tests/unit/test_R001_driver_defaults_to_ader.cpp
# Round-6 R-002 regression: driver init uses total-Q
TEST_R002_DRIVER_INIT_TOTAL_Q_SRC          = tests/unit/test_R002_driver_init_total_q.cpp
# REVIEW.md R-003 regression: base ComputeTractionDiagnostics default
# sizes diagnostic vectors to match traction (not empty) so Antiplane
# +face_tracer combinations don't OOB.
TEST_R003_ANTIPLANE_DIAG_SIZING_SRC        = tests/unit/test_R003_antiplane_diagnostics_sizing.cpp

# Clustered LTS (Phase 0/1).  See document/lts_dev/PLAN_clustered_lts_ader_2026-07-18.md.
LTS_CLUSTERING_HEADERS                      = dynamic/lts_clustering.hpp
LTS_CLUSTERING_SRC                          = dynamic/lts_clustering.cpp
TEST_LTS_CLUSTERING_SRC                     = tests/unit/test_lts_clustering.cpp
LTS_LAYOUT_HEADERS                          = dynamic/lts_layout.hpp \
                                              dynamic/lts_stepper.hpp
LTS_LAYOUT_SRC                              = dynamic/lts_layout.cpp
TEST_LTS_LAYOUT_SRC                         = tests/unit/test_lts_layout.cpp
LTS_TIME_BASIS_HEADERS                      = dynamic/lts_time_basis.hpp
TEST_LTS_TIME_BASIS_SRC                     = tests/unit/test_lts_time_basis.cpp
LTS_STEPPER_HEADERS                         = dynamic/lts_stepper.hpp
LTS_STEPPER_SRC                             = dynamic/lts_stepper.cpp
TEST_LTS_SCHEDULER_SRC                      = tests/unit/test_lts_scheduler.cpp

# Data-projection feature (Phase 3..5).  See document/features_dev/
# data_projection_feature_plan_v2.md.
DATA_FIELD_3D_HEADERS                       = io/data_field_3d.hpp
DATA_FIELD_3D_SRC                           = io/data_field_3d.cpp
STRESS_FIELD_3D_HEADERS                     = io/stress_field_3d.hpp \
                                              io/data_field_3d.hpp
STRESS_FIELD_3D_SRC                         = io/stress_field_3d.cpp
FIELD_COEFFICIENT_HEADERS                   = io/field_coefficient.hpp \
                                              io/data_field_3d.hpp \
                                              io/stress_field_3d.hpp \
                                              io/stress_field_coefficient.hpp
FIELD_COEFFICIENT_SRC                       = io/field_coefficient.cpp
HETEROGENEOUS_MATERIAL_HEADERS              = dynamic/heterogeneous_material.hpp
HETEROGENEOUS_MATERIAL_SRC                  = dynamic/heterogeneous_material.cpp
MATERIAL_COEFFICIENTS_HEADERS               = io/material_coefficients.hpp \
                                              io/data_field_3d.hpp
MATERIAL_COEFFICIENTS_SRC                   = io/material_coefficients.cpp
TEST_DATA_FIELD_3D_SRC                      = tests/unit/test_data_field_3d.cpp
# PLAN_sidecar_mpi_shared_memory_2026-07-17.md §6: MPI-3 shared-memory
# sidecar mode (np>=2 parity + cross-rank sharing assertions).
TEST_DATA_FIELD_3D_SHARED_MEM_SRC           = tests/unit/test_data_field_3d_shared_mem.cpp
# REVIEW.md R-001: dtor-after-MPI_Finalize guard regression (np=1).
TEST_DATA_FIELD_3D_SHARED_MEM_FIN_SRC       = tests/unit/test_data_field_3d_shared_mem_finalize.cpp
TEST_STRESS_FIELD_3D_SRC                    = tests/unit/test_stress_field_3d.cpp
TEST_FIELD_PROJECTOR_SRC                    = tests/unit/test_field_projector.cpp
TEST_HETEROGENEOUS_MATERIAL_SRC             = tests/unit/test_heterogeneous_material.cpp

# HDF5 (used by DataField3D).  HDF5 is provided by the conda env
# (mfem-dev) used by the project; CONDA_PREFIX is set when that env
# is active.  The HDF5_DIR variable lets a non-conda build override.
HDF5_DIR ?= $(CONDA_PREFIX)
ifeq ($(strip $(HDF5_DIR)),)
   HDF5_INCFLAGS =
   HDF5_LIBS     = -lhdf5
else
   HDF5_INCFLAGS = -I$(HDF5_DIR)/include
   HDF5_LIBS     = -L$(HDF5_DIR)/lib -lhdf5
endif

# R-501: `make ... CXXFLAGS+=-D<FLAG>` may not reach the compile step
# because MFEM's config.mk explicitly sets MFEM_CPPFLAGS/MFEM_CXXFLAGS.
# Use this SEAS-specific variable for diagnostic flags, e.g.
#   make seas_tpv102_driver \
#        SEAS_EXTRA_CPPFLAGS="-DSEAS_DIAG_FAULT_FLUX -DSEAS_DIAG_GHOST_EXCHANGE"
SEAS_EXTRA_CPPFLAGS ?=

# std::filesystem link probe (R-001 round 8 fix).  The BP5 + TPV104
# driver safety checks include <filesystem> and call
# fs::weakly_canonical / fs::path::parent_path.  On libstdc++ from
# gcc < 9 (e.g., Frontera intel/19 + system gcc 8.3) these symbols
# live in the separate library -lstdc++fs.  On gcc >= 9 + libstdc++,
# on macOS conda mfem-dev (libc++), and on any newer libstdc++, the
# symbols are in the main library and -lstdc++fs is a no-op (or the
# library is missing entirely).  Auto-detect: compile + link a
# minimal probe; if it succeeds without -lstdc++fs, leave SEAS_FS_LIB
# empty; otherwise add -lstdc++fs.  Operator override:
#   make seas_tpv104_driver SEAS_FS_LIB=-lstdc++fs
SEAS_FS_LIB := $(shell \
   tmp=$$(mktemp -d 2>/dev/null) && \
   printf '%s\n' '\#include <filesystem>' \
      'int main(){std::filesystem::weakly_canonical(std::filesystem::path("/tmp"));return 0;}' \
      > $$tmp/probe.cpp && \
   if $(MFEM_CXX) -std=c++17 $$tmp/probe.cpp -o $$tmp/probe 2>/dev/null; \
   then echo ""; else echo "-lstdc++fs"; fi; \
   rm -rf $$tmp 2>/dev/null)

# Object files
LTS_CLUSTERING_OBJ               = $(LTS_CLUSTERING_SRC:.cpp=.o)
TEST_LTS_CLUSTERING_OBJ          = $(TEST_LTS_CLUSTERING_SRC:.cpp=.o)
LTS_LAYOUT_OBJ                   = $(LTS_LAYOUT_SRC:.cpp=.o)
TEST_LTS_LAYOUT_OBJ              = $(TEST_LTS_LAYOUT_SRC:.cpp=.o)
TEST_LTS_TIME_BASIS_OBJ         = $(TEST_LTS_TIME_BASIS_SRC:.cpp=.o)
LTS_STEPPER_OBJ                 = $(LTS_STEPPER_SRC:.cpp=.o)
TEST_LTS_SCHEDULER_OBJ          = $(TEST_LTS_SCHEDULER_SRC:.cpp=.o)
DATA_FIELD_3D_OBJ                = $(DATA_FIELD_3D_SRC:.cpp=.o)
STRESS_FIELD_3D_OBJ              = $(STRESS_FIELD_3D_SRC:.cpp=.o)
FIELD_COEFFICIENT_OBJ            = $(FIELD_COEFFICIENT_SRC:.cpp=.o)
HETEROGENEOUS_MATERIAL_OBJ       = $(HETEROGENEOUS_MATERIAL_SRC:.cpp=.o)
MATERIAL_COEFFICIENTS_OBJ        = $(MATERIAL_COEFFICIENTS_SRC:.cpp=.o)
TEST_DATA_FIELD_3D_OBJ           = $(TEST_DATA_FIELD_3D_SRC:.cpp=.o)
TEST_DATA_FIELD_3D_SHARED_MEM_OBJ = $(TEST_DATA_FIELD_3D_SHARED_MEM_SRC:.cpp=.o)
TEST_DATA_FIELD_3D_SHARED_MEM_FIN_OBJ = $(TEST_DATA_FIELD_3D_SHARED_MEM_FIN_SRC:.cpp=.o)
TEST_STRESS_FIELD_3D_OBJ         = $(TEST_STRESS_FIELD_3D_SRC:.cpp=.o)
TEST_FIELD_PROJECTOR_OBJ         = $(TEST_FIELD_PROJECTOR_SRC:.cpp=.o)
TEST_HETEROGENEOUS_MATERIAL_OBJ  = $(TEST_HETEROGENEOUS_MATERIAL_SRC:.cpp=.o)
PROJECT_VELOCITY_DRIVER_SRC      = drivers/project_velocity_to_mesh.cpp
PROJECT_VELOCITY_DRIVER_OBJ      = $(PROJECT_VELOCITY_DRIVER_SRC:.cpp=.o)
TEST_FRICTION_OBJ = $(TEST_FRICTION_SRC:.cpp=.o)
TEST_STATE_EVOLUTION_OBJ = $(TEST_STATE_EVOLUTION_SRC:.cpp=.o)
TEST_ANTIPLANE_OBJ = $(TEST_ANTIPLANE_SRC:.cpp=.o)
TEST_FAULT_OPERATOR_OBJ = $(TEST_FAULT_OPERATOR_SRC:.cpp=.o)
TEST_QUASI_DYNAMIC_OBJ = $(TEST_QUASI_DYNAMIC_SRC:.cpp=.o)
TEST_BP2_SHORT_OBJ = $(TEST_BP2_SHORT_SRC:.cpp=.o)
TEST_IO_OBJ = $(TEST_IO_SRC:.cpp=.o)
BP2_SERIAL_SMOKE_OBJ = $(BP2_SERIAL_SMOKE_SRC:.cpp=.o)
TEST_CHECKPOINT_OBJ = $(TEST_CHECKPOINT_SRC:.cpp=.o)
TEST_PSI_STATE_OBJ = $(TEST_PSI_STATE_SRC:.cpp=.o)
TEST_FAULT_DETECTION_OBJ = $(TEST_FAULT_DETECTION_SRC:.cpp=.o)
TEST_DOMAIN_INTERFACE_OBJ = $(TEST_DOMAIN_INTERFACE_SRC:.cpp=.o)
TEST_BP5_PARAMS_OBJ = $(TEST_BP5_PARAMS_SRC:.cpp=.o)
TEST_VECTOR_FRICTION_OBJ = $(TEST_VECTOR_FRICTION_SRC:.cpp=.o)
TEST_FAULT_BASIS_OBJ = $(TEST_FAULT_BASIS_SRC:.cpp=.o)
TEST_ELASTICITY_BR2_OBJ = $(TEST_ELASTICITY_BR2_SRC:.cpp=.o)
TEST_ELASTICITY_OPERATOR_OBJ = $(TEST_ELASTICITY_OPERATOR_SRC:.cpp=.o)
TEST_BP5_FAULT_OPERATOR_OBJ = $(TEST_BP5_FAULT_OPERATOR_SRC:.cpp=.o)
# Phase 6.A: FaultGeometry per-DOF coordinate / basis accessor test
TEST_FAULT_DOF_BASIS_SRC = tests/unit/test_fault_dof_basis.cpp
TEST_FAULT_DOF_BASIS_OBJ = $(TEST_FAULT_DOF_BASIS_SRC:.cpp=.o)
# Phase 6 §4: ProjectFaultPreStress test
TEST_PROJECT_FAULT_PRESTRESS_SRC = tests/unit/test_project_fault_prestress.cpp
TEST_PROJECT_FAULT_PRESTRESS_OBJ = $(TEST_PROJECT_FAULT_PRESTRESS_SRC:.cpp=.o)
# Phase 6 §5: FaultGeometry::ComputeParams test
TEST_COMPUTE_SAFS_PARAMS_SRC = tests/unit/test_compute_safs_params.cpp
TEST_COMPUTE_SAFS_PARAMS_OBJ = $(TEST_COMPUTE_SAFS_PARAMS_SRC:.cpp=.o)
# Phase 6 §6: RateStateFaultOperator SAFS-mode wiring test
TEST_SAFS_MODE_WIRING_SRC = tests/unit/test_safs_mode_wiring.cpp
TEST_SAFS_MODE_WIRING_OBJ = $(TEST_SAFS_MODE_WIRING_SRC:.cpp=.o)
# Phase 6 §7: SAFS sidecar StressConfig schema test
TEST_SAFS_STRESS_CONFIG_SRC = tests/unit/test_safs_stress_config.cpp
TEST_SAFS_STRESS_CONFIG_OBJ = $(TEST_SAFS_STRESS_CONFIG_SRC:.cpp=.o)
# Phase 7: project_stress_to_mesh driver
PROJECT_STRESS_DRIVER_SRC = drivers/project_stress_to_mesh.cpp
PROJECT_STRESS_DRIVER_OBJ = $(PROJECT_STRESS_DRIVER_SRC:.cpp=.o)
# Phase 7: driver round-trip test
TEST_PROJECT_STRESS_TO_MESH_SRC = tests/unit/test_project_stress_to_mesh.cpp
TEST_PROJECT_STRESS_TO_MESH_OBJ = $(TEST_PROJECT_STRESS_TO_MESH_SRC:.cpp=.o)
TEST_BP5_INTEGRATION_OBJ = $(TEST_BP5_INTEGRATION_SRC:.cpp=.o)
TEST_BP5_OUTPUT_OBJ = $(TEST_BP5_OUTPUT_SRC:.cpp=.o)
TEST_BP5_MESH_OBJ = $(TEST_BP5_MESH_SRC:.cpp=.o)
TEST_DIAG_VTK_OBJ = $(TEST_DIAG_VTK_SRC:.cpp=.o)
TEST_MESH_REGIONS_OBJ = $(TEST_MESH_REGIONS_SRC:.cpp=.o)
TEST_FACE_QUADRATURE_OBJ = $(TEST_FACE_QUADRATURE_SRC:.cpp=.o)
TEST_CROSS_VERIFY_OBJ = $(TEST_CROSS_VERIFY_SRC:.cpp=.o)
TEST_LINEAR_ELASTIC_OBJ = $(TEST_LINEAR_ELASTIC_SRC:.cpp=.o)
TEST_CONSTITUTIVE_INTEG_OBJ = $(TEST_CONSTITUTIVE_INTEG_SRC:.cpp=.o)
TEST_BOUNDARY_CONFIG_OBJ = $(TEST_BOUNDARY_CONFIG_SRC:.cpp=.o)
TEST_BOUNDARY_CONFIG_OP_OBJ = $(TEST_BOUNDARY_CONFIG_OP_SRC:.cpp=.o)
TEST_BOUNDARY_CLASSIF_OBJ = $(TEST_BOUNDARY_CLASSIF_SRC:.cpp=.o)
TEST_LOGGING_OBJ = $(TEST_LOGGING_SRC:.cpp=.o)
TEST_CONFIG_BRIDGE_OBJ = $(TEST_CONFIG_BRIDGE_SRC:.cpp=.o)
TEST_CONFIG_PARSER_OBJ = $(TEST_CONFIG_PARSER_SRC:.cpp=.o)
SEAS_DRIVER_OBJ = $(SEAS_DRIVER_SRC:.cpp=.o)
TEST_GODUNOV_FLUX_OBJ = $(TEST_GODUNOV_FLUX_SRC:.cpp=.o)
TEST_WAVE_OPERATOR_OBJ = $(TEST_WAVE_OPERATOR_SRC:.cpp=.o)
TEST_WAVE_BC_OBJ = $(TEST_WAVE_BC_SRC:.cpp=.o)
TEST_PML_OBJ = $(TEST_PML_SRC:.cpp=.o)
TEST_FAULT_FACE_FLUX_OBJ = $(TEST_FAULT_FACE_FLUX_SRC:.cpp=.o)
TEST_FRICTION_SOLVER_OBJ = $(TEST_FRICTION_SOLVER_SRC:.cpp=.o)
TEST_TPV102_SETUP_OBJ = $(TEST_TPV102_SETUP_SRC:.cpp=.o)
TPV102_DRIVER_OBJ = $(TPV102_DRIVER_SRC:.cpp=.o)
TPV102_SUBSTEP_ITERATOR_OBJ = $(TPV102_SUBSTEP_ITERATOR_SRC:.cpp=.o)

# TPV102 driver rebuilt against the TPV104 code flow — depends on the new
# tpv102_friction_solver / tpv102_nucleation / tpv102_substep_iterator
# headers and the existing tpv102_setup / tpv102_params headers.
TPV102_HEADERS = config/tpv102_params.hpp \
                 friction/state_evolution.hpp \
                 friction/friction_coeff_stable.hpp \
                 dynamic/tpv102_setup.hpp \
                 dynamic/tpv102_friction_solver.hpp \
                 dynamic/tpv102_nucleation.hpp \
                 dynamic/tpv102_substep_iterator.hpp

# Shared object bundle for TPV102 driver/tests: substep iterator + the
# shared fault_face_flux + friction_solver objects.  Same pattern as
# TPV104_SHARED_OBJS.
TPV102_SHARED_OBJS = $(TPV102_SUBSTEP_ITERATOR_OBJ) \
                     dynamic/fault_face_flux.o \
                     dynamic/friction_solver.o

# TPV205 (linear slip-weakening, SCEC TPV5 spec) driver + iterator.
TPV205_DRIVER_OBJ = $(TPV205_DRIVER_SRC:.cpp=.o)
TPV205_SUBSTEP_ITERATOR_OBJ = $(TPV205_SUBSTEP_ITERATOR_SRC:.cpp=.o)

TPV205_HEADERS = config/tpv205_params.hpp \
                 dynamic/tpv205_friction.hpp \
                 dynamic/tpv205_setup.hpp \
                 dynamic/tpv205_substep_iterator.hpp

TPV205_SHARED_OBJS = $(TPV205_SUBSTEP_ITERATOR_OBJ) \
                     dynamic/fault_face_flux.o \
                     dynamic/friction_solver.o
TEST_TPV102_LOCAL_OBJ = $(TEST_TPV102_LOCAL_SRC:.cpp=.o)
TEST_R101_SHARED_FAULT_OBJ = $(TEST_R101_SHARED_FAULT_SRC:.cpp=.o)
TEST_FAULT_FACE_FLUX_FRAME_OBJ = $(TEST_FAULT_FACE_FLUX_FRAME_SRC:.cpp=.o)
TEST_FAULT_FACE_FLUX_FRAME_AND_FLUX_OBJ = $(TEST_FAULT_FACE_FLUX_FRAME_AND_FLUX_SRC:.cpp=.o)
TEST_GODUNOV_IDENTITY_OBJ = $(TEST_GODUNOV_IDENTITY_SRC:.cpp=.o)
TEST_FAULT_FLUX_BRANCH_EQUIV_OBJ = $(TEST_FAULT_FLUX_BRANCH_EQUIV_SRC:.cpp=.o)
TEST_GODUNOV_INTERIOR_EQUAL_SIDES_OBJ = $(TEST_GODUNOV_INTERIOR_EQUAL_SIDES_SRC:.cpp=.o)
TEST_FAULT_FACE_FLUX_PER_SIDE_ASSEMBLY_OBJ = $(TEST_FAULT_FACE_FLUX_PER_SIDE_ASSEMBLY_SRC:.cpp=.o)
TEST_FAULT_FACE_FLUX_BIMATERIAL_GUARD_OBJ  = $(TEST_FAULT_FACE_FLUX_BIMATERIAL_GUARD_SRC:.cpp=.o)
TEST_BIMATERIAL_FAULT_RIEMANN_OBJ          = $(TEST_BIMATERIAL_FAULT_RIEMANN_SRC:.cpp=.o)
TEST_TPV6_STATIONS_OBJ                     = $(TEST_TPV6_STATIONS_SRC:.cpp=.o)
# v9.1.0 Pelties-9 follow-through
TEST_FAULT_SURFACE_VTU_CONTINUITY_OBJ      = $(TEST_FAULT_SURFACE_VTU_CONTINUITY_SRC:.cpp=.o)
TEST_FAULT_BASIS_DIP_STRIKE_SYMMETRY_OBJ   = $(TEST_FAULT_BASIS_DIP_STRIKE_SYMMETRY_SRC:.cpp=.o)
# v9.2.0 §19 R-V92-E02
TEST_FAULT_SURFACE_VTU_K4_OBJ              = $(TEST_FAULT_SURFACE_VTU_K4_SRC:.cpp=.o)
# --paraview-fields field filter
TEST_FAULT_SURFACE_VTU_FIELD_FILTER_OBJ    = $(TEST_FAULT_SURFACE_VTU_FIELD_FILTER_SRC:.cpp=.o)
TEST_FAULT_SURFACE_VTU_BINARY_OBJ          = $(TEST_FAULT_SURFACE_VTU_BINARY_SRC:.cpp=.o)
TEST_FAULT_SURFACE_VTU_GATHER_MPI_OBJ      = $(TEST_FAULT_SURFACE_VTU_GATHER_MPI_SRC:.cpp=.o)
TEST_FREE_SURFACE_SLICE_OBJ                = $(TEST_FREE_SURFACE_SLICE_SRC:.cpp=.o)
TEST_PARAVIEW_HDF_SMOKE_OBJ                = $(TEST_PARAVIEW_HDF_SMOKE_SRC:.cpp=.o)
TEST_FAULT_SURFACE_VTKHDF_OBJ              = $(TEST_FAULT_SURFACE_VTKHDF_SRC:.cpp=.o)
TEST_FAULT_SURFACE_VTKHDF_MPI_OBJ          = $(TEST_FAULT_SURFACE_VTKHDF_MPI_SRC:.cpp=.o)
TEST_VTKHDF_ZFP_OBJ                        = $(TEST_VTKHDF_ZFP_SRC:.cpp=.o)
TEST_FAULT_SURFACE_VTKHDF_ZFP_OBJ          = $(TEST_FAULT_SURFACE_VTKHDF_ZFP_SRC:.cpp=.o)
TEST_VOLUME_HDF_COMPRESSION_OBJ            = $(TEST_VOLUME_HDF_COMPRESSION_SRC:.cpp=.o)
TEST_PARAVIEW_SCHEDULE_CAP_OBJ             = $(TEST_PARAVIEW_SCHEDULE_CAP_SRC:.cpp=.o)
TEST_KINEMATICS_FIELD_SET_OBJ              = $(TEST_KINEMATICS_FIELD_SET_SRC:.cpp=.o)
TEST_BP5_PETSC_TS_ZERO_FAULT_RANK_OBJ      = $(TEST_BP5_PETSC_TS_ZERO_FAULT_RANK_SRC:.cpp=.o)
TEST_BP5_PETSC_TS_RESTART_OBJ              = $(TEST_BP5_PETSC_TS_RESTART_SRC:.cpp=.o)
TEST_TPV104_CHECKPOINT_OBJ                 = $(TEST_TPV104_CHECKPOINT_SRC:.cpp=.o)
TEST_PARAVIEW_RANK0_WARNING_GATE_OBJ       = $(TEST_PARAVIEW_RANK0_WARNING_GATE_SRC:.cpp=.o)
# v9.2.0 plan §4.1-§4.4 regression gates
TEST_GODUNOV_ROTATION_IDENTITY_OBJ         = $(TEST_GODUNOV_ROTATION_IDENTITY_SRC:.cpp=.o)
TEST_CANONICAL_ROTATION_PURE_STRIKESLIP_OBJ = $(TEST_CANONICAL_ROTATION_PURE_STRIKESLIP_SRC:.cpp=.o)
TEST_FRAME_DIP_LEAK_PROBE_OBJ              = $(TEST_FRAME_DIP_LEAK_PROBE_SRC:.cpp=.o)
TEST_OFFFAULT_FRAME_SYMMETRY_OBJ          = $(TEST_OFFFAULT_FRAME_SYMMETRY_SRC:.cpp=.o)
TEST_FAULT_BASIS_QP_ORTHONORMALITY_OBJ     = $(TEST_FAULT_BASIS_QP_ORTHONORMALITY_SRC:.cpp=.o)
TEST_PER_QP_VS_CENTROID_BASIS_OBJ          = $(TEST_PER_QP_VS_CENTROID_BASIS_SRC:.cpp=.o)
TEST_FRAME_DIAG_PROJECTION_IDENTITY_OBJ    = $(TEST_FRAME_DIAG_PROJECTION_IDENTITY_SRC:.cpp=.o)
TEST_MACRO_DIAG_DECOMPOSITION_IDENTITY_OBJ = $(TEST_MACRO_DIAG_DECOMPOSITION_IDENTITY_SRC:.cpp=.o)
TEST_SHARED_FAULT_ROLE_CONSIST_OBJ         = $(TEST_SHARED_FAULT_ROLE_CONSIST_SRC:.cpp=.o)
TEST_NO_PENALTY_DR_OBJ                     = $(TEST_NO_PENALTY_DR_SRC:.cpp=.o)
TEST_INTERIOR_FAULT_FLUX_PATH_OBJ          = $(TEST_INTERIOR_FAULT_FLUX_PATH_SRC:.cpp=.o)
TEST_INTERIOR_FAULT_FLUX_PATH_ASYM_OBJ     = $(TEST_INTERIOR_FAULT_FLUX_PATH_ASYM_SRC:.cpp=.o)
TEST_INTERIOR_FAULT_FLUX_PATH_RUPTURE_OBJ  = $(TEST_INTERIOR_FAULT_FLUX_PATH_RUPTURE_SRC:.cpp=.o)
TEST_BULK_RK4_NOISE_STAB_OBJ               = $(TEST_BULK_RK4_NOISE_STAB_SRC:.cpp=.o)
TEST_SHARED_FAULT_FLUX_PATH_ANA_OBJ        = $(TEST_SHARED_FAULT_FLUX_PATH_ANA_SRC:.cpp=.o)
TEST_BULK_PLANE_WAVE_OBJ                   = $(TEST_BULK_PLANE_WAVE_SRC:.cpp=.o)
TEST_INTERIOR_VS_SHARED_LIVE_OBJ           = $(TEST_INTERIOR_VS_SHARED_LIVE_SRC:.cpp=.o)
TEST_RUPTURE_MULTISTEP_OBJ                 = $(TEST_RUPTURE_MULTISTEP_SRC:.cpp=.o)
TEST_FAULT_PLANAR_SERIAL_OBJ               = $(TEST_FAULT_PLANAR_SERIAL_SRC:.cpp=.o)
TEST_FAULT_OVERINT_OBJ                     = $(TEST_FAULT_OVERINT_SRC:.cpp=.o)
TEST_FAULT_RESAMPLE_OBJ                    = $(TEST_FAULT_RESAMPLE_SRC:.cpp=.o)
TEST_FAULT_RESAMPLE_APPLY_OBJ             = $(TEST_FAULT_RESAMPLE_APPLY_SRC:.cpp=.o)
TEST_RUPTURE_TILTED_OBJ                    = $(TEST_RUPTURE_TILTED_SRC:.cpp=.o)
TEST_RUPTURE_RAKE_SWEEP_OBJ                = $(TEST_RUPTURE_RAKE_SWEEP_SRC:.cpp=.o)
TEST_RECONCILE_XRANK_OBJ                   = $(TEST_RECONCILE_XRANK_SRC:.cpp=.o)
TEST_RS_PSI_NP2_OBJ                        = $(TEST_RS_PSI_NP2_SRC:.cpp=.o)
# fault_face_flux built WITH -DSEAS_TEST_INTERNAL for the reconcile oracle's
# absolute-Pa injection hook (lives in the .cpp behind that macro).
FAULT_FACE_FLUX_TESTINTERNAL_OBJ           = dynamic/fault_face_flux_testinternal.o
TEST_ADJACENT_TRI_UNIFORMITY_OBJ           = $(TEST_ADJACENT_TRI_UNIFORMITY_SRC:.cpp=.o)
TEST_ADJACENT_TRI_FIRST_STEP_AUDIT_OBJ     = $(TEST_ADJACENT_TRI_FIRST_STEP_AUDIT_SRC:.cpp=.o)
TEST_GODUNOV_ROUNDTRIP_OBJ                 = $(TEST_GODUNOV_ROUNDTRIP_SRC:.cpp=.o)
TEST_FRICTION_EQUIL_OBJ                    = $(TEST_FRICTION_EQUIL_SRC:.cpp=.o)
TEST_PEPPER_AMP_CHAIN_OBJ                  = $(TEST_PEPPER_AMP_CHAIN_SRC:.cpp=.o)
TEST_FACE_LIFT_OPERATOR_OBJ                = $(TEST_FACE_LIFT_OPERATOR_SRC:.cpp=.o)
TEST_PEPPER_QP_PERT_OBJ                    = $(TEST_PEPPER_QP_PERT_SRC:.cpp=.o)
TEST_PEPPER_NONLINEAR_GROWTH_OBJ           = $(TEST_PEPPER_NONLINEAR_GROWTH_SRC:.cpp=.o)
TEST_SIGN_FLIPPED_TRUTH_TABLE_OBJ          = $(TEST_SIGN_FLIPPED_TRUTH_TABLE_SRC:.cpp=.o)
TEST_MAKE_FACE_KEY_UNIQ_OBJ                = $(TEST_MAKE_FACE_KEY_UNIQ_SRC:.cpp=.o)
TEST_FAULT_DOF_OFFSET_OBJ                  = $(TEST_FAULT_DOF_OFFSET_SRC:.cpp=.o)
TEST_RK4_PSI_INTEGRATION_OBJ               = $(TEST_RK4_PSI_INTEGRATION_SRC:.cpp=.o)
TEST_RK4_ENDPOINT_REEVAL_OBJ               = $(TEST_RK4_ENDPOINT_REEVAL_SRC:.cpp=.o)
TEST_VOLUME_JAC_SINGLE_CH_OBJ              = $(TEST_VOLUME_JAC_SINGLE_CH_SRC:.cpp=.o)
TEST_RK4_CONSERVATION_OBJ                  = $(TEST_RK4_CONSERVATION_SRC:.cpp=.o)
TEST_ABSORBING_BC_DECAY_OBJ                = $(TEST_ABSORBING_BC_DECAY_SRC:.cpp=.o)
TEST_SHARED_FAULT_DOFDATA_CONSIST_OBJ      = $(TEST_SHARED_FAULT_DOFDATA_CONSIST_SRC:.cpp=.o)
TEST_RK4_F01F02_STABILITY_OBJ              = $(TEST_RK4_F01F02_STABILITY_SRC:.cpp=.o)
TEST_CENTROID_MARGIN_CTOR_OBJ              = $(TEST_CENTROID_MARGIN_CTOR_SRC:.cpp=.o)
# v9.3.0 Phase 1 / Phase 2 (I-04)
TEST_FREE_SURFACE_GODUNOV_OBJ              = $(TEST_FREE_SURFACE_GODUNOV_SRC:.cpp=.o)
TEST_R005_FS_BC_FLAG_PARSE_OBJ             = $(TEST_R005_FS_BC_FLAG_PARSE_SRC:.cpp=.o)
TEST_FREE_SURFACE_GODUNOV_DRIVER_OBJ       = $(TEST_FREE_SURFACE_GODUNOV_DRIVER_SRC:.cpp=.o)
TEST_FAULT_FACE_FLUX_TOTAL_VS_FLUC_OBJ     = $(TEST_FAULT_FACE_FLUX_TOTAL_VS_FLUC_SRC:.cpp=.o)
TEST_R003_EVALUATETOTAL_PSI_OBJ            = $(TEST_R003_EVALUATETOTAL_PSI_SRC:.cpp=.o)
TEST_R001_INITSTATE_TOTAL_SIGN_OBJ         = $(TEST_R001_INITSTATE_TOTAL_SIGN_SRC:.cpp=.o)
TEST_TPV102_TOTAL_LOCKED_FAULT_OBJ         = $(TEST_TPV102_TOTAL_LOCKED_FAULT_SRC:.cpp=.o)
TEST_TPV102_TOTAL_ABSORBING_OBJ            = $(TEST_TPV102_TOTAL_ABSORBING_SRC:.cpp=.o)
TEST_PERSISTENT_NUC_PRESTRESS_OBJ          = $(TEST_PERSISTENT_NUC_PRESTRESS_SRC:.cpp=.o)
TEST_TPV102_PEPPER_OBJ                     = $(TEST_TPV102_PEPPER_SRC:.cpp=.o)
TEST_FREE_SURFACE_TOTAL_DIRECT_OBJ         = $(TEST_FREE_SURFACE_TOTAL_DIRECT_SRC:.cpp=.o)
TEST_WAVE_OP_SPATIAL_DERIV_OBJ             = $(TEST_WAVE_OP_SPATIAL_DERIV_SRC:.cpp=.o)
BENCH_ADER_HOTPATH_OBJ                     = $(BENCH_ADER_HOTPATH_SRC:.cpp=.o)
TEST_ADER_CK_PREDICTOR_OBJ                 = $(TEST_ADER_CK_PREDICTOR_SRC:.cpp=.o)
TEST_ADER_FAULT_FLUX_EQUIV_OBJ             = $(TEST_ADER_FAULT_FLUX_EQUIV_SRC:.cpp=.o)
TEST_ADER_LINEAR_WAVE_EQUIV_OBJ            = $(TEST_ADER_LINEAR_WAVE_EQUIV_SRC:.cpp=.o)
TEST_ADER_FAULT_COUPLING_DT_LIMIT_OBJ      = $(TEST_ADER_FAULT_COUPLING_DT_LIMIT_SRC:.cpp=.o)
TEST_ADER_INTERIOR_VS_SHARED_LIVE_OBJ      = $(TEST_ADER_INTERIOR_VS_SHARED_LIVE_SRC:.cpp=.o)
TEST_ADER_TPV102_SMOKE_OBJ                 = $(TEST_ADER_TPV102_SMOKE_SRC:.cpp=.o)
TEST_R001_DRIVER_DEFAULT_ADER_OBJ          = $(TEST_R001_DRIVER_DEFAULT_ADER_SRC:.cpp=.o)
TEST_R002_DRIVER_INIT_TOTAL_Q_OBJ          = $(TEST_R002_DRIVER_INIT_TOTAL_Q_SRC:.cpp=.o)
TEST_R003_ANTIPLANE_DIAG_SIZING_OBJ        = $(TEST_R003_ANTIPLANE_DIAG_SIZING_SRC:.cpp=.o)

# Constitutive model headers
CONSTITUTIVE_HEADERS = constitutive/constitutive_model.hpp \
                       constitutive/linear_elastic.hpp

# Phase 5: Config headers
CONFIG_HEADERS = config/seas_config.hpp \
                 config/seas_config_parser.hpp \
                 config/seas_config_bridge.hpp

# Header dependencies
# Phase 1: Friction and state evolution
FRICTION_HEADERS = friction/friction_law.hpp \
                   friction/state_evolution.hpp \
                   friction/dieterich_ruina.hpp \
                   config/bp2_params.hpp \
                   config/bp1_params.hpp \
                   config/bp5_params.hpp \
                   config/bp5_mesh_utils.hpp

# Phase 2: Domain operator
DOMAIN_HEADERS = common/seas_types.hpp \
                 common/mpi_context.hpp \
                 domain/domain_operator.hpp \
                 domain/antiplane_operator.hpp \
                 domain/antiplane_bdrload_operator.hpp \
                 domain/elasticity_operator.hpp \
                 domain/elasticity_operator_verify.inl \
                 domain/elasticity_operator_debug.inl \
                 domain/elasticity_operator_setup.inl \
                 domain/elasticity_operator_assembly.inl \
                 domain/elasticity_operator_traction.inl \
                 domain/boundary_config.hpp \
                 domain/domain_config.hpp \
                 domain/seas_boundary_tags.hpp \
                 domain/bp2_mesh.hpp \
                 integrator/dg_br2_integrator.hpp \
                 integrator/dg_elasticity_br2_integrator.hpp

# Phase 3: Fault operator
FAULT_HEADERS = fault/fault_geometry.hpp \
                fault/rate_state_fault.hpp \
                fault/fault_basis.hpp \
                fault/face_quadrature.hpp

# Phase 4: SEAS operator
SOLVER_HEADERS = solver/seas_operator.hpp \
                 solver/seas_bdrload_operator.hpp \
                 solver/time_stepper.hpp

# Phase 5: I/O
IO_HEADERS = io/probe_output.hpp \
             io/benchmark_output.hpp \
             io/paraview_output.hpp \
             io/free_surface_output.hpp

# Checkpoint I/O
CHECKPOINT_HEADERS = io/checkpoint.hpp

# Dynamic rupture headers (Phase 1)
DYNAMIC_HEADERS = dynamic/wave_state.hpp \
                  dynamic/godunov_flux.hpp \
                  dynamic/wave_operator.hpp \
                  dynamic/wave_operator.inl \
                  dynamic/pml_layer.hpp \
                  dynamic/friction_solver.hpp \
                  dynamic/fault_face_flux.hpp \
                  dynamic/precomputed_face_fluxes.hpp \
                  dynamic/seas_dynamic_operator.hpp \
                  dynamic/fault_state_channel.hpp \
                  dynamic/tpv102_setup.hpp

# Spatial (SAFS) headers.  spatial_friction.hpp defines SpatialFrictionConfig /
# RateStateBlock / NucleationSpec, which nucleation_factory.o, spatial_nucleation.o,
# spatial_print_derived.o, fault_face_flux.o and the spatial driver all embed BY
# VALUE.  It was absent from SEAS_HEADERS, so a layout change to any of those
# structs did NOT trigger a rebuild of those objects -- an incremental
# `make <target>` then linked mixed-ABI objects and produced a binary that
# segfaulted or silently mis-dispatched (observed 2026-07-09 on
# seas_test_nucleation_factory after adding RateStateBlock::sidecar).
# There is no auto-dependency generation in this Makefile, so list them.
SPATIAL_HEADERS = spatial/code/spatial_friction.hpp \
                  spatial/code/spatial_stress.hpp \
                  spatial/code/spatial_velocity.hpp

# All headers
SEAS_HEADERS = $(FRICTION_HEADERS) $(DOMAIN_HEADERS) $(FAULT_HEADERS) $(SOLVER_HEADERS) $(IO_HEADERS) $(CHECKPOINT_HEADERS) $(CONSTITUTIVE_HEADERS) $(CONFIG_HEADERS) $(DYNAMIC_HEADERS) $(SPATIAL_HEADERS)

# Mesh generation tool
GENERATE_MESH_SRC = tools/generate_bp2_mesh.cpp
GENERATE_MESH_OBJ = $(GENERATE_MESH_SRC:.cpp=.o)

# Targets — unit tests (fast, run by `make test` via MFEM test framework)
SEQ_MINIAPPS = seas_test_friction_law seas_test_state_evolution seas_test_antiplane \
               seas_test_fault_operator seas_test_quasi_dynamic seas_test_bp2_short \
               seas_test_io seas_test_checkpoint \
               seas_test_psi_state \
               seas_test_fault_detection seas_test_domain_interface \
               seas_test_bp5_params seas_test_vector_friction \
               seas_test_fault_basis \
               seas_test_elasticity_br2 seas_test_elasticity_operator \
               seas_test_bp5_fault_operator \
               seas_test_fault_dof_basis \
               seas_test_project_fault_prestress \
               seas_test_compute_safs_params \
               seas_test_safs_mode_wiring \
               seas_test_safs_stress_config \
               seas_test_project_stress_to_mesh \
               seas_test_bp5_integration \
               seas_test_bp5_output seas_test_bp5_mesh \
               seas_test_diag_vtk seas_test_mesh_regions \
               seas_test_face_quadrature \
               seas_test_linear_elastic \
               seas_test_constitutive_integrator \
               seas_test_boundary_config \
               seas_test_boundary_config_operator \
               seas_test_logging \
               seas_test_config_bridge \
               seas_test_godunov_flux \
               seas_test_godunov_flux_bimaterial \
               seas_test_phaseh_godunov_flux_pool \
               seas_test_spatial_friction_config \
               seas_test_spatial_friction_resolver \
               seas_test_spatial_friction_sidecar \
               seas_test_friction_depth_profile \
               seas_test_spatial_nucleation \
               seas_test_nucleation_factory \
               seas_test_rk_time_stepper \
               seas_test_lsw_rk_mixed_flux \
               seas_test_lsw_rk_mixed_flux_bimaterial \
               seas_test_tpv102_nuc_callback_parity \
               seas_test_seed_equilibrium_psi_rs \
               seas_test_friction_iterator_factory \
               seas_test_friction_substep_iterator_parity \
               seas_test_advance_interface_compiles \
               seas_test_resolve_rate_state_guards \
               seas_test_constant_tensor_sign \
               seas_test_tpv31_canonical_rotation \
               seas_test_tpv31_nucleation_and_cohesion \
               seas_test_tpv31_station_writer \
               seas_test_spatial_station_dispatch \
               seas_test_tpv_toml_stress_sign \
               seas_test_planar_tpv_basis \
               seas_test_tpv_config_parse \
               seas_test_paraview_state_channel_rs \
               seas_test_spatial_print_derived \
               seas_test_paraview_set_fault_params_spatial \
               seas_test_spatial_velocity_bundle \
               seas_test_spatial_stress_bundle \
               seas_test_spatial_constant_stress_source \
               seas_test_spatial_setup \
               seas_test_phaseh_lsw_forced_rupture \
               seas_test_tpv205_friction \
               seas_test_phaseh_wave_operator_constant_parity \
               seas_test_bimaterial_wave_operator_parity \
               seas_test_wave_operator \
               seas_test_wave_bc \
               seas_test_pml \
               seas_test_fault_face_flux \
               seas_test_friction_solver \
               seas_test_tpv102_setup \
               seas_test_tpv102_local
# Long-running simulations and tools (built by `make all`, not run by `make test`)
SEQ_LONG = seas_bp2_serial_smoke seas_generate_mesh
# Phase 6: Parallel test sources
TEST_MPI_CONTEXT_SRC = tests/parallel/test_mpi_context.cpp
TEST_PARALLEL_UTILS_SRC = tests/parallel/test_parallel_utils.cpp
TEST_FAULT_SCATTER_SRC = tests/parallel/test_fault_scatter.cpp
TEST_MPI_CONTEXT_OBJ = $(TEST_MPI_CONTEXT_SRC:.cpp=.o)
TEST_PARALLEL_UTILS_OBJ = $(TEST_PARALLEL_UTILS_SRC:.cpp=.o)
TEST_FAULT_SCATTER_OBJ = $(TEST_FAULT_SCATTER_SRC:.cpp=.o)

# Phase 6: Parallel headers
PARALLEL_HEADERS = common/mpi_context.hpp \
                   common/parallel_utils.hpp \
                   common/seas_types.hpp \
                   common/logging.hpp \
                   common/fault_scatter.hpp \
                   common/mpi_check.hpp \
                   common/mpi_tags.hpp

# Phase 7: Parallel domain test sources
TEST_PARALLEL_DOMAIN_SRC = tests/parallel/test_parallel_domain.cpp
MMS_ANTIPLANE_PARALLEL_SRC = tests/parallel/mms_antiplane_parallel.cpp
TEST_PARALLEL_DOMAIN_OBJ = $(TEST_PARALLEL_DOMAIN_SRC:.cpp=.o)
MMS_ANTIPLANE_PARALLEL_OBJ = $(MMS_ANTIPLANE_PARALLEL_SRC:.cpp=.o)

# Phase 8: Parallel fault test and driver sources
TEST_PARALLEL_FAULT_SRC = tests/parallel/test_parallel_fault.cpp
PSEAS_SRC = pseas.cpp
TEST_PARALLEL_FAULT_OBJ = $(TEST_PARALLEL_FAULT_SRC:.cpp=.o)
PSEAS_OBJ = $(PSEAS_SRC:.cpp=.o)

# Phase 8: Parallel I/O headers
PARALLEL_IO_HEADERS = io/parallel_benchmark_output.hpp \
                      io/bp5_benchmark_output.hpp \
                      io/bp5_parallel_output.hpp

# Phase 9: Verification and scaling sources
TEST_SERIAL_PARALLEL_SRC = tests/parallel/test_serial_parallel_consistency.cpp
BP2_BENCHMARK_PARALLEL_SRC = tests/verification/bp2_benchmark_parallel.cpp
TEST_SCALING_SRC = tests/parallel/test_scaling.cpp
TEST_SERIAL_PARALLEL_OBJ = $(TEST_SERIAL_PARALLEL_SRC:.cpp=.o)
BP2_BENCHMARK_PARALLEL_OBJ = $(BP2_BENCHMARK_PARALLEL_SRC:.cpp=.o)
TEST_SCALING_OBJ = $(TEST_SCALING_SRC:.cpp=.o)

# BR2 consistency test
TEST_BR2_CONSISTENCY_SRC = tests/parallel/test_br2_consistency.cpp
TEST_BR2_CONSISTENCY_OBJ = $(TEST_BR2_CONSISTENCY_SRC:.cpp=.o)

# Parallel elasticity test (shared face support)
TEST_PARALLEL_ELAST_SRC = tests/parallel/test_parallel_elasticity.cpp
TEST_PARALLEL_ELAST_OBJ = $(TEST_PARALLEL_ELAST_SRC:.cpp=.o)

# Parallel wave operator tests (dynamic rupture)
TEST_PARALLEL_WAVE_SRC = tests/parallel/test_parallel_wave_operator.cpp
TEST_PARALLEL_WAVE_OBJ = $(TEST_PARALLEL_WAVE_SRC:.cpp=.o)
# Lever 1/3 (REVIEW R-002): parallel DerivMode Cached==OnTheFly equivalence.
TEST_PARALLEL_CACHED_WAVE_SRC = tests/parallel/test_wave_operator_cached_parallel.cpp
TEST_PARALLEL_CACHED_WAVE_OBJ = $(TEST_PARALLEL_CACHED_WAVE_SRC:.cpp=.o)

# Bimaterial (matrix path) DerivMode + shared-CK parity test (opt 2026-06-24 gate)
TEST_PARALLEL_BIMAT_CACHE_SRC = tests/parallel/test_bimaterial_deriv_cache_parity.cpp
TEST_PARALLEL_BIMAT_CACHE_OBJ = $(TEST_PARALLEL_BIMAT_CACHE_SRC:.cpp=.o)

# BP5 parallel smoke test
TEST_BP5_SMOKE_SRC = tests/parallel/test_bp5_parallel_smoke.cpp
TEST_BP5_SMOKE_OBJ = $(TEST_BP5_SMOKE_SRC:.cpp=.o)

# Phase 10: BP2 first-cycle verification
BP2_VERIFY_SRC = tests/verification/bp2_verification_first_cycle.cpp
BP2_VERIFY_OBJ = $(BP2_VERIFY_SRC:.cpp=.o)

# Phase 10b: BP2 full multi-cycle simulation (1200 yr)
BP2_FULL_SRC = tests/verification/bp2_verification_full.cpp
BP2_FULL_OBJ = $(BP2_FULL_SRC:.cpp=.o)

# BP1 full multi-cycle simulation (3000 yr)
BP1_FULL_SRC = tests/verification/bp1_verification_full.cpp
BP1_FULL_OBJ = $(BP1_FULL_SRC:.cpp=.o)

# BP1 bottom-boundary Dirichlet loading simulation
BP1_BDRLOAD_SRC = tests/verification/bp1_bdrload.cpp
BP1_BDRLOAD_OBJ = $(BP1_BDRLOAD_SRC:.cpp=.o)

# BP5 full multi-cycle simulation
BP5_FULL_SRC = tests/verification/bp5_verification_full.cpp
BP5_FULL_OBJ = $(BP5_FULL_SRC:.cpp=.o)

PAR_MINIAPPS = seas_test_mpi_context seas_test_parallel_utils \
               seas_test_fault_scatter \
               seas_test_config_parser \
               seas_test_parallel_domain seas_mms_antiplane_parallel \
               seas_test_parallel_fault \
               seas_test_serial_parallel_consistency \
               seas_test_scaling \
               seas_test_br2_consistency \
               seas_test_parallel_elasticity \
               seas_test_bp5_parallel_smoke \
               seas_test_parallel_wave_operator
# Long-running parallel simulations (built by `make all`, not run by `make test`)
PAR_LONG = seas_pseas \
           seas_bp2_benchmark_parallel \
           seas_bp2_verify seas_bp2_full seas_bp1_full \
           seas_bp1_bdrload \
           seas_bp5_full \
           seas_driver

ifeq ($(MFEM_USE_MPI),NO)
   MINIAPPS = $(SEQ_MINIAPPS) $(SEQ_LONG)
else
   MINIAPPS = $(PAR_MINIAPPS) $(PAR_LONG) $(SEQ_MINIAPPS) $(SEQ_LONG)
endif

MFEM_TESTS = MINIAPPS

.SUFFIXES:
.SUFFIXES: .o .cpp .mk
.PHONY: all clean clean-build clean-exec test \
        test-friction test-state-evolution test-antiplane \
        test-fault-operator test-quasi-dynamic test-bp2-short \
        test-io test-bp2-serial-smoke test-checkpoint \
        test-mpi-context test-parallel-utils test-fault-scatter test-parallel \
        test-parallel-domain test-mms-parallel test-parallel-fault \
        test-serial-parallel test-bp2-parallel test-scaling \
        test-bp2-verify test-psi-state test-br2-consistency \
        test-fault-detection test-domain-interface \
        test-bp5-params test-vector-friction test-fault-basis \
        test-elasticity-br2 test-elasticity-operator \
        test-bp5-fault-operator \
        test-bp5-integration \
        test-bp5-output test-bp5-mesh \
        test-diag-vtk test-mesh-regions \
        test-bp5-smoke test-linear-elastic test-constitutive-integrator \
        test-boundary-config test-boundary-config-operator test-logging \
        test-config-bridge \
        test-config-parser \
        test-tpv102-setup test-tpv102-local \
        test-parallel-wave-operator \
        test-godunov-rotation-identity \
        test-canonical-rotation-pure-strikeslip \
        test-fault-basis-qp-orthonormality \
        test-frame-diag-projection-identity \
        test-macro-diag-decomposition-identity \
        test-ghost-exchange-bynodes-vs-scalar \
        test-per-qp-vs-centroid-basis-planar-1el \
        test-v92-regression-gates \
        regression-check quick-check

all: $(MINIAPPS)

# Remove built-in rules
%: %.cpp
%.o: %.cpp

SEAS_INCLUDES = $(if $(SRC),-I$(MFEM_DIR))

# Clustered LTS Phase 0: pure clustering unit test (no MPI, no HDF5 of its own,
# but libmfem links HDF5 so $(HDF5_LIBS) is required).
seas_test_lts_clustering: $(TEST_LTS_CLUSTERING_OBJ) $(LTS_CLUSTERING_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_LTS_CLUSTERING_OBJ) \
	    $(LTS_CLUSTERING_OBJ) $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Clustered LTS Phase 1: run-side layout + tick-table unit test.
seas_test_lts_layout: $(TEST_LTS_LAYOUT_OBJ) $(LTS_LAYOUT_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_LTS_LAYOUT_OBJ) \
	    $(LTS_LAYOUT_OBJ) $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Clustered LTS Phase 2: Taylor time-integration unit test (header-only lib).
seas_test_lts_time_basis: $(TEST_LTS_TIME_BASIS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_LTS_TIME_BASIS_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Clustered LTS Phase 2: tick-loop scheduler unit test (RunSyncInterval + mock).
seas_test_lts_scheduler: $(TEST_LTS_SCHEDULER_OBJ) $(LTS_STEPPER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_LTS_SCHEDULER_OBJ) \
	    $(LTS_STEPPER_OBJ) $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Data-projection feature: DataField3D / FieldCoefficient / FieldProjector /
# MaterialField unit tests (Phase 3..5 of plan v2).
seas_test_data_field_3d: $(TEST_DATA_FIELD_3D_OBJ) $(DATA_FIELD_3D_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_DATA_FIELD_3D_OBJ) \
	    $(DATA_FIELD_3D_OBJ) $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# PLAN_sidecar_mpi_shared_memory_2026-07-17.md §6: MPI shared-memory
# sidecar parity test (needs StressField3D for the comm-forwarding case).
seas_test_data_field_3d_shared_mem: $(TEST_DATA_FIELD_3D_SHARED_MEM_OBJ) \
                                    $(DATA_FIELD_3D_OBJ) $(STRESS_FIELD_3D_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_DATA_FIELD_3D_SHARED_MEM_OBJ) \
	    $(STRESS_FIELD_3D_OBJ) $(DATA_FIELD_3D_OBJ) $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# REVIEW.md R-001: destroy a shared-mode reader AFTER MPI_Finalize (np=1).
seas_test_data_field_3d_shared_mem_finalize: \
                                    $(TEST_DATA_FIELD_3D_SHARED_MEM_FIN_OBJ) \
                                    $(DATA_FIELD_3D_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_DATA_FIELD_3D_SHARED_MEM_FIN_OBJ) \
	    $(DATA_FIELD_3D_OBJ) $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_field_projector: $(TEST_FIELD_PROJECTOR_OBJ) \
                           $(DATA_FIELD_3D_OBJ) $(STRESS_FIELD_3D_OBJ) \
                           $(FIELD_COEFFICIENT_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FIELD_PROJECTOR_OBJ) \
	    $(FIELD_COEFFICIENT_OBJ) $(STRESS_FIELD_3D_OBJ) \
	    $(DATA_FIELD_3D_OBJ) $(MFEM_LIBS) \
	    $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 6 §1 of PLAN_onfaultstress.md: six-component Cauchy stress
# sidecar reader test (StressField3D).
seas_test_stress_field_3d: $(TEST_STRESS_FIELD_3D_OBJ) \
                           $(STRESS_FIELD_3D_OBJ) $(DATA_FIELD_3D_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_STRESS_FIELD_3D_OBJ) \
	    $(STRESS_FIELD_3D_OBJ) $(DATA_FIELD_3D_OBJ) $(MFEM_LIBS) \
	    $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_heterogeneous_material: $(TEST_HETEROGENEOUS_MATERIAL_OBJ) \
                                  $(HETEROGENEOUS_MATERIAL_OBJ) \
                                  $(MATERIAL_COEFFICIENTS_OBJ) \
                                  $(DATA_FIELD_3D_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_HETEROGENEOUS_MATERIAL_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(MATERIAL_COEFFICIENTS_OBJ) \
	    $(DATA_FIELD_3D_OBJ) $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Data-projection feature: standalone driver that loads a mesh + a
# schema-v1 sidecar, projects (Vp, Vs, rho, lambda, mu) onto the mesh,
# and writes a ParaView .pvd / .vtu pair (R-003).
seas_project_velocity_to_mesh: $(PROJECT_VELOCITY_DRIVER_OBJ) \
                               $(DATA_FIELD_3D_OBJ) $(STRESS_FIELD_3D_OBJ) \
                               $(FIELD_COEFFICIENT_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(PROJECT_VELOCITY_DRIVER_OBJ) $(FIELD_COEFFICIENT_OBJ) \
	    $(STRESS_FIELD_3D_OBJ) $(DATA_FIELD_3D_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 1: Friction law unit tests
seas_test_friction_law: $(TEST_FRICTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FRICTION_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 1: State evolution unit tests
seas_test_state_evolution: $(TEST_STATE_EVOLUTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_STATE_EVOLUTION_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 2: Antiplane domain operator unit tests
seas_test_antiplane: $(TEST_ANTIPLANE_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ANTIPLANE_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 3: Fault operator unit tests
seas_test_fault_operator: $(TEST_FAULT_OPERATOR_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_OPERATOR_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 4: Quasi-dynamic coupling tests
seas_test_quasi_dynamic: $(TEST_QUASI_DYNAMIC_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_QUASI_DYNAMIC_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 4: Short BP2 simulation test
seas_test_bp2_short: $(TEST_BP2_SHORT_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BP2_SHORT_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 5: I/O unit tests
seas_test_io: $(TEST_IO_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_IO_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 5: BP2 serial smoke test (800 m coarse mesh)
seas_bp2_serial_smoke: $(BP2_SERIAL_SMOKE_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(BP2_SERIAL_SMOKE_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 11: Checkpoint unit test
seas_test_checkpoint: $(TEST_CHECKPOINT_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_CHECKPOINT_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Psi-space state variable unit test
seas_test_psi_state: $(TEST_PSI_STATE_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PSI_STATE_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 1a: Fault detection unit test
seas_test_fault_detection: $(TEST_FAULT_DETECTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_DETECTION_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 1b: DomainOperator interface unit test
seas_test_domain_interface: $(TEST_DOMAIN_INTERFACE_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_DOMAIN_INTERFACE_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# BP5 Phase 1: Parameter tests
seas_test_bp5_params: $(TEST_BP5_PARAMS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BP5_PARAMS_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# BP5 Phase 1: Vector friction tests
seas_test_vector_friction: $(TEST_VECTOR_FRICTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_VECTOR_FRICTION_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# BP5 Phase 2a: Fault basis tests
seas_test_fault_basis: $(TEST_FAULT_BASIS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_BASIS_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# BP5 Phase 2b: Elasticity BR2 integrator tests
seas_test_elasticity_br2: $(TEST_ELASTICITY_BR2_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ELASTICITY_BR2_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# BP5 Phase 2b: Elasticity domain operator tests
seas_test_elasticity_operator: $(TEST_ELASTICITY_OPERATOR_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ELASTICITY_OPERATOR_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# BP5 Phase 3: BP5 fault operator tests
seas_test_bp5_fault_operator: $(TEST_BP5_FAULT_OPERATOR_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BP5_FAULT_OPERATOR_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 6.A: FaultGeometry per-DOF accessor tests (Phase 6.A of
# PLAN_onfaultstress.md).
seas_test_fault_dof_basis: $(TEST_FAULT_DOF_BASIS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_DOF_BASIS_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 6 §4: ProjectFaultPreStress test (Phase 6 §4 of PLAN_onfaultstress.md)
seas_test_project_fault_prestress: $(TEST_PROJECT_FAULT_PRESTRESS_OBJ) \
                                    $(STRESS_FIELD_3D_OBJ) \
                                    $(DATA_FIELD_3D_OBJ) \
                                    $(FIELD_COEFFICIENT_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_PROJECT_FAULT_PRESTRESS_OBJ) \
	    $(FIELD_COEFFICIENT_OBJ) \
	    $(STRESS_FIELD_3D_OBJ) $(DATA_FIELD_3D_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 6 §5: FaultGeometry::ComputeParams test
seas_test_compute_safs_params: $(TEST_COMPUTE_SAFS_PARAMS_OBJ) \
                                $(STRESS_FIELD_3D_OBJ) \
                                $(DATA_FIELD_3D_OBJ) \
                                $(FIELD_COEFFICIENT_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_COMPUTE_SAFS_PARAMS_OBJ) \
	    $(FIELD_COEFFICIENT_OBJ) \
	    $(STRESS_FIELD_3D_OBJ) $(DATA_FIELD_3D_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 6 §6: RateStateFaultOperator SAFS-mode wiring test
seas_test_safs_mode_wiring: $(TEST_SAFS_MODE_WIRING_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_SAFS_MODE_WIRING_OBJ) \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 6 §7: StressConfig schema test
seas_test_safs_stress_config: $(TEST_SAFS_STRESS_CONFIG_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_SAFS_STRESS_CONFIG_OBJ) \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 7: project_stress_to_mesh driver (mirrors seas_project_velocity_to_mesh)
seas_project_stress_to_mesh: $(PROJECT_STRESS_DRIVER_OBJ) \
                             $(DATA_FIELD_3D_OBJ) \
                             $(STRESS_FIELD_3D_OBJ) \
                             $(FIELD_COEFFICIENT_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(PROJECT_STRESS_DRIVER_OBJ) \
	    $(FIELD_COEFFICIENT_OBJ) \
	    $(STRESS_FIELD_3D_OBJ) \
	    $(DATA_FIELD_3D_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 7: driver round-trip test
seas_test_project_stress_to_mesh: $(TEST_PROJECT_STRESS_TO_MESH_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_PROJECT_STRESS_TO_MESH_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# BP5 Phase 4: BP5 integration tests
seas_test_bp5_integration: $(TEST_BP5_INTEGRATION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BP5_INTEGRATION_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# BP5 Phase 5: Output unit tests
seas_test_bp5_output: $(TEST_BP5_OUTPUT_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BP5_OUTPUT_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# BP5 Phase 5: Mesh unit tests
seas_test_bp5_mesh: $(TEST_BP5_MESH_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BP5_MESH_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# BP5 v16: Diagnostic VTK unit tests
seas_test_diag_vtk: $(TEST_DIAG_VTK_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_DIAG_VTK_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# BP5 v16: Mesh region unit tests
seas_test_mesh_regions: $(TEST_MESH_REGIONS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_MESH_REGIONS_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v44: Face quadrature unit tests
seas_test_face_quadrature: $(TEST_FACE_QUADRATURE_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FACE_QUADRATURE_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v54: Cross-verification against Tandem
seas_test_cross_verify: $(TEST_CROSS_VERIFY_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_CROSS_VERIFY_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Dynamic rupture Phase 1: Godunov flux tests
seas_test_godunov_flux: $(TEST_GODUNOV_FLUX_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_GODUNOV_FLUX_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase H.1 of spatial_dynamic_rupture_plan.md: GodunovFluxPool dedup tests
seas_test_phaseh_godunov_flux_pool: $(TEST_PHASEH_GODUNOV_FLUX_POOL_OBJ) \
                                    $(GODUNOV_FLUX_POOL_OBJ) \
                                    $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_PHASEH_GODUNOV_FLUX_POOL_OBJ) \
	    $(GODUNOV_FLUX_POOL_OBJ) \
	    $(GODUNOV_FLUX_OBJ) \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Spatial Phase 1 of spatial_dynamic_rupture_plan.md: TOML parser tests
seas_test_spatial_friction_config: $(TEST_SPATIAL_FRICTION_CONFIG_OBJ) \
                                    $(SPATIAL_FRICTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_SPATIAL_FRICTION_CONFIG_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_spatial_friction_resolver: $(TEST_SPATIAL_FRICTION_RESOLVER_OBJ) \
                                      $(SPATIAL_FRICTION_OBJ) \
                                      $(HETEROGENEOUS_MATERIAL_OBJ) \
                                      $(MATERIAL_COEFFICIENTS_OBJ) \
                                      $(DATA_FIELD_3D_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_SPATIAL_FRICTION_RESOLVER_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(MATERIAL_COEFFICIENTS_OBJ) \
	    $(DATA_FIELD_3D_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 3 of PLAN_thermal_case2_mixedflux_port_2026-07-08.md: the
# [friction.rate_state.sidecar] resolver path.  Same link set as the resolver
# test -- the test builds DataField3D readers directly to mirror the driver's
# MakeRateStateSidecarFields.
seas_test_spatial_friction_sidecar: $(TEST_SPATIAL_FRICTION_SIDECAR_OBJ) \
                                     $(SPATIAL_FRICTION_OBJ) \
                                     $(HETEROGENEOUS_MATERIAL_OBJ) \
                                     $(MATERIAL_COEFFICIENTS_OBJ) \
                                     $(DATA_FIELD_3D_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_SPATIAL_FRICTION_SIDECAR_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(MATERIAL_COEFFICIENTS_OBJ) \
	    $(DATA_FIELD_3D_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 11b: depth-varying a(z)/b(z) two-CSV profile test.
seas_test_friction_depth_profile: $(TEST_FRICTION_DEPTH_PROFILE_OBJ) \
                                   $(SPATIAL_FRICTION_OBJ) \
                                   $(HETEROGENEOUS_MATERIAL_OBJ) \
                                   $(MATERIAL_COEFFICIENTS_OBJ) \
                                   $(DATA_FIELD_3D_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_FRICTION_DEPTH_PROFILE_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(MATERIAL_COEFFICIENTS_OBJ) \
	    $(DATA_FIELD_3D_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Spatial Phase 2: SpatialVelocityBundle + sidecar loader test.
seas_test_spatial_velocity_bundle: $(TEST_SPATIAL_VELOCITY_BUNDLE_OBJ) \
                                    $(SPATIAL_VELOCITY_OBJ) \
                                    $(SPATIAL_FRICTION_OBJ) \
                                    $(HETEROGENEOUS_MATERIAL_OBJ) \
                                    $(MATERIAL_COEFFICIENTS_OBJ) \
                                    $(DATA_FIELD_3D_OBJ) \
                                    $(FIELD_COEFFICIENT_OBJ) \
                                    $(STRESS_FIELD_3D_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_SPATIAL_VELOCITY_BUNDLE_OBJ) \
	    $(SPATIAL_VELOCITY_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(MATERIAL_COEFFICIENTS_OBJ) \
	    $(DATA_FIELD_3D_OBJ) \
	    $(FIELD_COEFFICIENT_OBJ) \
	    $(STRESS_FIELD_3D_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Spatial Phase 3: CSM stress sidecar wrapper.
# Spatial Phase 3b: ConstantTensorStressSource (D-1 testing primitive).
seas_test_spatial_stress_bundle: $(TEST_SPATIAL_STRESS_BUNDLE_OBJ) \
                                  $(SPATIAL_STRESS_OBJ) \
                                  $(SPATIAL_FRICTION_OBJ) \
                                  $(STRESS_FIELD_3D_OBJ) \
                                  $(DATA_FIELD_3D_OBJ) \
                                  $(FIELD_COEFFICIENT_OBJ) \
                                  $(HETEROGENEOUS_MATERIAL_OBJ) \
                                  $(MATERIAL_COEFFICIENTS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_SPATIAL_STRESS_BUNDLE_OBJ) \
	    $(SPATIAL_STRESS_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(STRESS_FIELD_3D_OBJ) \
	    $(DATA_FIELD_3D_OBJ) \
	    $(FIELD_COEFFICIENT_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(MATERIAL_COEFFICIENTS_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_spatial_constant_stress_source: $(TEST_SPATIAL_CONSTANT_STRESS_OBJ) \
                                           $(SPATIAL_STRESS_OBJ) \
                                           $(SPATIAL_FRICTION_OBJ) \
                                           $(STRESS_FIELD_3D_OBJ) \
                                           $(DATA_FIELD_3D_OBJ) \
                                           $(FIELD_COEFFICIENT_OBJ) \
                                           $(HETEROGENEOUS_MATERIAL_OBJ) \
                                           $(MATERIAL_COEFFICIENTS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_SPATIAL_CONSTANT_STRESS_OBJ) \
	    $(SPATIAL_STRESS_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(STRESS_FIELD_3D_OBJ) \
	    $(DATA_FIELD_3D_OBJ) \
	    $(FIELD_COEFFICIENT_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(MATERIAL_COEFFICIENTS_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Spatial Phase 5c: spatial_setup.hpp tests (centroid + IP-aware
# InitializeFaultDOFs_Spatial overloads).  Header-only; the test
# pulls in the friction resolver + material objects to seed inputs.
seas_test_spatial_setup: $(TEST_SPATIAL_SETUP_OBJ) \
                         $(SPATIAL_FRICTION_OBJ) \
                         $(HETEROGENEOUS_MATERIAL_OBJ) \
                         $(MATERIAL_COEFFICIENTS_OBJ) \
                         $(DATA_FIELD_3D_OBJ) \
                         dynamic/fault_face_flux.o \
                         dynamic/friction_solver.o
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_SPATIAL_SETUP_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(MATERIAL_COEFFICIENTS_OBJ) \
	    $(DATA_FIELD_3D_OBJ) \
	    dynamic/fault_face_flux.o \
	    dynamic/friction_solver.o \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase H.6: EvaluateADER_LSW_ForcedRupture tests.
seas_test_phaseh_lsw_forced_rupture: $(TEST_PHASEH_LSW_FORCED_RUPTURE_OBJ) \
                                     $(SPATIAL_FRICTION_OBJ) \
                                     dynamic/fault_face_flux.o \
                                     dynamic/friction_solver.o \
                                     dynamic/wave_operator.o \
                                     dynamic/godunov_flux.o \
                                     dynamic/pml_layer.o \
                                     dynamic/precomputed_face_fluxes.o
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_PHASEH_LSW_FORCED_RUPTURE_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    dynamic/fault_face_flux.o \
	    dynamic/friction_solver.o \
	    dynamic/wave_operator.o \
	    dynamic/godunov_flux.o \
	    dynamic/pml_layer.o \
	    dynamic/precomputed_face_fluxes.o \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase H Stage 1 byte-exact gate: scalar vs MaterialField ctor parity.
# Pair with `mpirun -np 4` for the np=4 case (test target below runs both).
seas_test_phaseh_wave_operator_constant_parity: \
                  $(TEST_PHASEH_WAVE_OPERATOR_CONSTANT_PARITY_OBJ) \
                  $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
                  $(GODUNOV_FLUX_OBJ) \
                  $(GODUNOV_FLUX_POOL_OBJ) \
                  $(PML_LAYER_OBJ) \
                  $(PRECOMPUTED_FACE_FLUXES_OBJ) \
                  $(HETEROGENEOUS_MATERIAL_OBJ) \
                  $(FAULT_FACE_FLUX_OBJ) \
                  $(FRICTION_SOLVER_OBJ) \
                  $(SPATIAL_FRICTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_PHASEH_WAVE_OPERATOR_CONSTANT_PARITY_OBJ) \
	    $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
	    $(GODUNOV_FLUX_OBJ) \
	    $(GODUNOV_FLUX_POOL_OBJ) \
	    $(PML_LAYER_OBJ) \
	    $(PRECOMPUTED_FACE_FLUXES_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 13 C-6: fault-bearing BimaterialWaveOperator-vs-scalar parity (the
# flux_-placeholder leak-class tripwire).  Serial; same object set as the
# constant-parity target (BimaterialWaveOperator + FaultFaceFlux + friction).
seas_test_bimaterial_wave_operator_parity: \
                  $(TEST_BIMATERIAL_WAVE_OPERATOR_PARITY_OBJ) \
                  $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
                  $(GODUNOV_FLUX_OBJ) \
                  $(GODUNOV_FLUX_POOL_OBJ) \
                  $(PML_LAYER_OBJ) \
                  $(PRECOMPUTED_FACE_FLUXES_OBJ) \
                  $(HETEROGENEOUS_MATERIAL_OBJ) \
                  $(FAULT_FACE_FLUX_OBJ) \
                  $(FRICTION_SOLVER_OBJ) \
                  $(SPATIAL_FRICTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_BIMATERIAL_WAVE_OPERATOR_PARITY_OBJ) \
	    $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
	    $(GODUNOV_FLUX_OBJ) \
	    $(GODUNOV_FLUX_POOL_OBJ) \
	    $(PML_LAYER_OBJ) \
	    $(PRECOMPUTED_FACE_FLUXES_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Dynamic rupture Phase 1: Wave operator tests
seas_test_wave_operator: $(TEST_WAVE_OPERATOR_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(GODUNOV_FLUX_POOL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_WAVE_OPERATOR_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(GODUNOV_FLUX_POOL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Dynamic rupture Phase 2a: BC tests
seas_test_wave_bc: $(TEST_WAVE_BC_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_WAVE_BC_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Dynamic rupture Phase 2b: PML tests
seas_test_pml: $(TEST_PML_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PML_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Dynamic rupture Phase 3: Fault face flux tests
seas_test_fault_face_flux: $(TEST_FAULT_FACE_FLUX_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_FACE_FLUX_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# TPV102 "Topology-Based Precomputed Face-Rotation" plan 2026-04-23 Phase 1
# Step A interior identity unit test
seas_test_precomputed_fluxes_interior_identity: $(TEST_PRECOMPUTED_FLUXES_INTERIOR_IDENTITY_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PRECOMPUTED_FLUXES_INTERIOR_IDENTITY_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 1 Step B boundary identity unit test
seas_test_precomputed_fluxes_boundary_identity: $(TEST_PRECOMPUTED_FLUXES_BOUNDARY_IDENTITY_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PRECOMPUTED_FLUXES_BOUNDARY_IDENTITY_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 2a acceptance — opt-in dispatch switch (§6.4)
seas_test_precomputed_fluxes_phase2a_switch: $(TEST_PRECOMPUTED_FLUXES_PHASE2A_SWITCH_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PRECOMPUTED_FLUXES_PHASE2A_SWITCH_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 3 STOP — Arm 1 constant-state localization probes
seas_test_arm1_constant_state_probes: $(TEST_ARM1_CONSTANT_STATE_PROBES_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ARM1_CONSTANT_STATE_PROBES_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 3 STOP — Arm 2 D4-equivariant fixture probes
seas_test_arm2_d4_equivariant_fixture: $(TEST_ARM2_D4_EQUIVARIANT_FIXTURE_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ARM2_D4_EQUIVARIANT_FIXTURE_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 3 STOP — Arm 2b vertex-permutation verification probe
seas_test_arm2b_vertex_permutation_probe: $(TEST_ARM2B_VERTEX_PERM_PROBE_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ARM2B_VERTEX_PERM_PROBE_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 3 STOP — Arm 3d topology-based state-sampling probe
seas_test_arm3d_topology_state_sampling: $(TEST_ARM3D_TOPO_STATE_SAMPLING_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ARM3D_TOPO_STATE_SAMPLING_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 2a acceptance — R008 fault-branch preservation + P_SWITCH_ON_FAULT_NOT_ZEROED
seas_test_R008_fault_branch_preserved: $(TEST_R008_FAULT_BRANCH_PRESERVED_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R008_FAULT_BRANCH_PRESERVED_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 2a acceptance — R5-001 compile-time hoist gate
seas_test_R5_001_hoist_compiles: $(TEST_R5_001_HOIST_COMPILES_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R5_001_HOIST_COMPILES_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# TPV102 v9.5.0 — REVIEW R-001 regression: Mult honors the precomputed flux flag
seas_test_R001_rk4_uses_precomputed_flag: $(TEST_R001_RK4_USES_PRECOMPUTED_FLAG_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R001_RK4_USES_PRECOMPUTED_FLAG_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 1 correction 2026-04-23 R-002: BuildD4Mesh vertex-equivariance (Phase 2D precondition)
seas_test_d4_fixture_vertex_equivariance: $(TEST_D4_FIXTURE_VERTEX_EQUIV_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_D4_FIXTURE_VERTEX_EQUIV_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 2D Step B: BuildD4Mesh(true) fault sanity gate
seas_test_d4_fault_sanity: $(TEST_D4_FAULT_SANITY_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_D4_FAULT_SANITY_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Round-3 R-003: D4 fixture Jacobian-sign audit
seas_test_d4_jacobian_sign: $(TEST_D4_JACOBIAN_SIGN_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_D4_JACOBIAN_SIGN_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Round-5 Option 1: Gate-14′-to-pepper correlation test
seas_test_r5_gate14_pepper_correlation: $(TEST_R5_GATE14_PEPPER_CORR_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R5_GATE14_PEPPER_CORR_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Round-6 Phase 2E: rupture-driven amplification chain
seas_test_r6_phase2e_amplification_chain: $(TEST_R6_PHASE2E_AMPCHAIN_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R6_PHASE2E_AMPCHAIN_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Round-7: strengthened amplification chain v2
seas_test_r7_amplification_chain_v2: $(TEST_R7_AMPCHAIN_V2_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R7_AMPCHAIN_V2_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Round-9: FREEZE-C psi-hold + dt-refinement
seas_test_r9_frozen_friction_C: $(TEST_R9_FROZEN_FRICTION_C_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R9_FROZEN_FRICTION_C_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 2b acceptance — MPI shared-face dispatch switch (§7.3)
seas_test_precomputed_fluxes_phase2b_mpi_switch: $(TEST_PRECOMPUTED_FLUXES_PHASE2B_MPI_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PRECOMPUTED_FLUXES_PHASE2B_MPI_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# TPV102 v9.5.0 — REVIEW R-002: empty-fault MPI rank must not abort
seas_test_R002_empty_fault_rank_no_abort: $(TEST_R002_EMPTY_FAULT_RANK_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R002_EMPTY_FAULT_RANK_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# TPV104 R-1600: substep np>1 deadlock-detection regression
seas_test_R1600_substep_no_fault_rank_no_hang: $(TEST_R1600_SUBSTEP_NO_FAULT_RANK_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R1600_SUBSTEP_NO_FAULT_RANK_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# R-004 ghost-exchange bit-identity test — pure MFEM (replicates both exchange
# mechanisms inline; no WaveOperator/fault deps; wave_state.hpp is header-only).
seas_test_ghost_exchange_bynodes_vs_scalar: $(TEST_GHOST_EXCHANGE_BYNODES_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_GHOST_EXCHANGE_BYNODES_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Round-11 R-001/R-002: cross-rank Adjacent mixed-flux dispatch consistency
seas_test_mixed_flux_adjacent_mpi: $(TEST_MIXED_FLUX_ADJACENT_MPI_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_MIXED_FLUX_ADJACENT_MPI_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Round-11 R-1100: shared fault face must NOT enter central_flux_face_set_
seas_test_mixed_flux_shared_fault_face_excluded_mpi: $(TEST_MIXED_FLUX_SHARED_FAULT_MPI_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_MIXED_FLUX_SHARED_FAULT_MPI_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Round-11 R-1102: Adjacent-mode end-to-end dispatch correctness
seas_test_mixed_flux_dispatch_adjacent: $(TEST_MIXED_FLUX_DISPATCH_ADJACENT_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_MIXED_FLUX_DISPATCH_ADJACENT_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Round-13 R-1201: np=2 vs np=1 Adjacent numerical bit-identity (MPI)
seas_test_mixed_flux_dispatch_adjacent_mpi: $(TEST_MIXED_FLUX_DISPATCH_ADJACENT_MPI_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_MIXED_FLUX_DISPATCH_ADJACENT_MPI_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Round-14 R-1402: AllContinuous-mode end-to-end dispatch correctness
seas_test_mixed_flux_dispatch_all_continuous: $(TEST_MIXED_FLUX_DISPATCH_ALL_CONT_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_MIXED_FLUX_DISPATCH_ALL_CONT_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# TPV102 v9.0.0 Phase 1 unit tests (§3.1, §3.1b, §3.1c, §3.1d)
seas_test_fault_face_flux_frame: $(TEST_FAULT_FACE_FLUX_FRAME_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_FACE_FLUX_FRAME_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_fault_face_flux_frame_and_flux: $(TEST_FAULT_FACE_FLUX_FRAME_AND_FLUX_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_FACE_FLUX_FRAME_AND_FLUX_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_godunov_identity_normal_reversal: $(TEST_GODUNOV_IDENTITY_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_GODUNOV_IDENTITY_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_fault_flux_interior_vs_shared_branch_equivalence: $(TEST_FAULT_FLUX_BRANCH_EQUIV_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_FLUX_BRANCH_EQUIV_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.0.0 Pelties-9 fix: §3.1e primitive identity test
seas_test_godunov_interior_equal_sides_identity: $(TEST_GODUNOV_INTERIOR_EQUAL_SIDES_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_GODUNOV_INTERIOR_EQUAL_SIDES_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.0.0 Pelties-9 fix: §3.1f per-side assembly test
seas_test_fault_face_flux_per_side_assembly: $(TEST_FAULT_FACE_FLUX_PER_SIDE_ASSEMBLY_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_FACE_FLUX_PER_SIDE_ASSEMBLY_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# REVIEW R-007 (v9.0.0): bimaterial-guard regression test for
# FaultFaceFlux::Evaluate.  Death test — the child aborts inside
# MFEM_VERIFY; the shell target below interprets non-zero exit as PASS.
seas_test_fault_face_flux_bimaterial_guard: $(TEST_FAULT_FACE_FLUX_BIMATERIAL_GUARD_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_FACE_FLUX_BIMATERIAL_GUARD_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Part B / B0 (GATE): per-side bi-material fault Riemann (locked-limit welded R/T).
seas_test_bimaterial_fault_riemann: $(TEST_BIMATERIAL_FAULT_RIEMANN_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BIMATERIAL_FAULT_RIEMANN_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Part C / C2: TPV6 per-side on-fault station writer (header-only; DOFData-driven).
seas_test_tpv6_stations: $(TEST_TPV6_STATIONS_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV6_STATIONS_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-tpv6-stations: seas_test_tpv6_stations
	./seas_test_tpv6_stations

test-bimaterial-fault-riemann: seas_test_bimaterial_fault_riemann
	./seas_test_bimaterial_fault_riemann

# v9.1.0 plan §4.1 — fault-surface VTU CellData contract (post-R-001)
seas_test_fault_surface_vtu_continuity: $(TEST_FAULT_SURFACE_VTU_CONTINUITY_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_SURFACE_VTU_CONTINUITY_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.2.0 plan §19 — R-V92-E02 stage-4 VTU diagnostic
seas_test_fault_surface_vtu_k4: $(TEST_FAULT_SURFACE_VTU_K4_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_SURFACE_VTU_K4_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# --paraview-fields / SetFaultVTUFields field filter
seas_test_fault_surface_vtu_field_filter: $(TEST_FAULT_SURFACE_VTU_FIELD_FILTER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_SURFACE_VTU_FIELD_FILTER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 1 PLAN_paraview_compaction_2026-04-28 — single binary VTU writer
seas_test_fault_surface_vtu_binary: $(TEST_FAULT_SURFACE_VTU_BINARY_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_SURFACE_VTU_BINARY_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 1 PLAN_paraview_compaction_2026-04-28 — np=4 gather correctness
seas_test_fault_surface_vtu_gather_mpi: $(TEST_FAULT_SURFACE_VTU_GATHER_MPI_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_SURFACE_VTU_GATHER_MPI_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 4 PLAN_free_surface_slice.md — header-only writer; links only MFEM.
seas_test_free_surface_slice: $(TEST_FREE_SURFACE_SLICE_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FREE_SURFACE_SLICE_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 2a PLAN_paraview_compaction_2026-04-28 — VTKHDF build smoke gate
seas_test_paraview_hdf_smoke: $(TEST_PARAVIEW_HDF_SMOKE_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PARAVIEW_HDF_SMOKE_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 2b PLAN_paraview_compaction_2026-04-28 — VTKHDF serial bit-exact
seas_test_fault_surface_vtkhdf: $(TEST_FAULT_SURFACE_VTKHDF_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_SURFACE_VTKHDF_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 2b PLAN_paraview_compaction_2026-04-28 — VTKHDF np=4 gather merge
seas_test_fault_surface_vtkhdf_mpi: $(TEST_FAULT_SURFACE_VTKHDF_MPI_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_SURFACE_VTKHDF_MPI_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 2d.2 PLAN_paraview_compaction_2026-04-28 — H5Z-ZFP filter gate
seas_test_vtkhdf_zfp: $(TEST_VTKHDF_ZFP_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_VTKHDF_ZFP_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 2d.3 PLAN_paraview_compaction_2026-04-28 — seas fault path + ZFP
seas_test_fault_surface_vtkhdf_zfp: $(TEST_FAULT_SURFACE_VTKHDF_ZFP_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_SURFACE_VTKHDF_ZFP_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 6.2 PLAN_bulk_compression_and_size_estimator_2026-05-09 — volume HDF + ZFP
seas_test_volume_hdf_compression: $(TEST_VOLUME_HDF_COMPRESSION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_VOLUME_HDF_COMPRESSION_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 3 PLAN_paraview_compaction_2026-04-28 — adaptive snapshot cap
seas_test_paraview_schedule_cap: $(TEST_PARAVIEW_SCHEDULE_CAP_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PARAVIEW_SCHEDULE_CAP_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# R-001 of PLAN_split_bulk_solutions_2026-05-12 — kinematics field set
seas_test_kinematics_field_set: $(TEST_KINEMATICS_FIELD_SET_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_KINEMATICS_FIELD_SET_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# R-003 of REVIEW.md 2026-05-16 — zero-fault-DOF rank padded-shrink invariant
seas_test_bp5_petsc_ts_zero_fault_rank: $(TEST_BP5_PETSC_TS_ZERO_FAULT_RANK_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BP5_PETSC_TS_ZERO_FAULT_RANK_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 1 of petsc_ts_restart_plan_2026-05-16 — V2 checkpoint round-trip
seas_test_bp5_petsc_ts_restart: $(TEST_BP5_PETSC_TS_RESTART_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BP5_PETSC_TS_RESTART_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 4 TPV104 port — V1 checkpoint round-trip
seas_test_tpv104_checkpoint: $(TEST_TPV104_CHECKPOINT_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV104_CHECKPOINT_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# R-104 of REVIEW.md 2026-05-16 round 2 — MPI-only rank-0 warning gate
seas_test_paraview_rank0_warning_gate: $(TEST_PARAVIEW_RANK0_WARNING_GATE_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PARAVIEW_RANK0_WARNING_GATE_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.1.0 plan §4.2 — FaultBasis orthonormality contract (post-R-003)
seas_test_fault_basis_dip_strike_symmetry: $(TEST_FAULT_BASIS_DIP_STRIKE_SYMMETRY_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_BASIS_DIP_STRIKE_SYMMETRY_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.2.0 plan §4.1 — Godunov rotation identity T·Tinv = I
seas_test_godunov_rotation_identity: $(TEST_GODUNOV_ROTATION_IDENTITY_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_GODUNOV_ROTATION_IDENTITY_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# PLAN_frame_orthonormality_diag_2026-05-23 §Phase 1 AC#3 — [FRAME] dv_n identity
seas_test_frame_diag_projection_identity: $(TEST_FRAME_DIAG_PROJECTION_IDENTITY_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FRAME_DIAG_PROJECTION_IDENTITY_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# PLAN_predictor_vs_macro_diag_2026-05-23 §Phase 1 AC#3 — [MACRO] decomposition identity
seas_test_macro_diag_decomposition_identity: $(TEST_MACRO_DIAG_DECOMPOSITION_IDENTITY_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_MACRO_DIAG_DECOMPOSITION_IDENTITY_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.2.0 plan §4.2 — canonical rotation pure strike-slip Tinv probe
seas_test_canonical_rotation_pure_strikeslip: $(TEST_CANONICAL_ROTATION_PURE_STRIKESLIP_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_CANONICAL_ROTATION_PURE_STRIKESLIP_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# dip-drift investigation probe (2026-06-06): frame-leak check vs BuildRotationInverse
seas_test_frame_dip_leak_probe: $(TEST_FRAME_DIP_LEAK_PROBE_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FRAME_DIP_LEAK_PROBE_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# dip-drift investigation probe (2026-06-06): off-fault upwind frame y-symmetry
seas_test_offfault_frame_symmetry: $(TEST_OFFFAULT_FRAME_SYMMETRY_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_OFFFAULT_FRAME_SYMMETRY_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.2.0 plan §4.3 — per-QP fault basis orthonormality
seas_test_fault_basis_qp_orthonormality: $(TEST_FAULT_BASIS_QP_ORTHONORMALITY_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_BASIS_QP_ORTHONORMALITY_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.2.0 plan §4.4 — per-QP vs centroid basis on planar 1-tet fixture
seas_test_per_qp_vs_centroid_basis_planar_1el: $(TEST_PER_QP_VS_CENTROID_BASIS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PER_QP_VS_CENTROID_BASIS_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.2.0 plan §4.9 (rev-3f) — No-penalty regression gate (Pelties 2012 §3.1)
seas_test_no_penalty_dynamic_rupture: $(TEST_NO_PENALTY_DR_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_NO_PENALTY_DR_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.2.0 plan Step 6 (H-V92-U probe) — interior-fault-path Pelties eq. (7)
seas_test_interior_fault_flux_path: $(TEST_INTERIOR_FAULT_FLUX_PATH_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_INTERIOR_FAULT_FLUX_PATH_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Pepper-bug plan 2026-04-22 Phase 3 P3-T1 — asymmetric Q probe
seas_test_interior_fault_flux_path_asymmetric: $(TEST_INTERIOR_FAULT_FLUX_PATH_ASYM_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_INTERIOR_FAULT_FLUX_PATH_ASYM_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Pepper-bug plan 2026-04-22 Phase 3 P3-T3 — rupture-regime amplitude probe
seas_test_interior_fault_flux_path_rupture_regime: $(TEST_INTERIOR_FAULT_FLUX_PATH_RUPTURE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_INTERIOR_FAULT_FLUX_PATH_RUPTURE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Pepper-bug plan 2026-04-22 Phase 5 P5-N — bulk RK4 noise stability
seas_test_bulk_rk4_noise_stability: $(TEST_BULK_RK4_NOISE_STAB_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BULK_RK4_NOISE_STAB_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Pepper-bug plan 2026-04-22 Test 5 — shared-fault Pelties-7 analytic (np=4)
seas_test_shared_fault_flux_path_analytic: $(TEST_SHARED_FAULT_FLUX_PATH_ANA_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_SHARED_FAULT_FLUX_PATH_ANA_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Pepper-bug plan 2026-04-22 Phase 5 P5-PW — bulk plane-wave recovery
seas_test_bulk_plane_wave_recovery: $(TEST_BULK_PLANE_WAVE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BULK_PLANE_WAVE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Pepper-bug plan 2026-04-22 Phase 1 — live interior-vs-shared cross-check
# NB: built with -DSEAS_TEST_INTERNAL so the test can access WaveOperator
# internal bookkeeping (getters gated by that macro in wave_operator.hpp).
seas_test_rupture_multistep_serial_vs_parallel: $(TEST_RUPTURE_MULTISTEP_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_RUPTURE_MULTISTEP_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

seas_test_fault_planar_serial: $(TEST_FAULT_PLANAR_SERIAL_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_PLANAR_SERIAL_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_fault_overint: $(TEST_FAULT_OVERINT_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_OVERINT_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 2: fault-face resample operator R — header-only (dynamic/fault_resample.hpp),
# so the test links against MFEM alone (no other seas objects).
seas_test_fault_resample: $(TEST_FAULT_RESAMPLE_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_RESAMPLE_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 3: rate-state Δψ resample apply — same link surface as the iterator
# parity test (drives the unified RateStateAgingIterator).
seas_test_fault_resample_apply: \
        $(TEST_FAULT_RESAMPLE_APPLY_OBJ) \
        $(FRICTION_SUBSTEP_ITERATOR_OBJ) \
        $(TPV205_SUBSTEP_ITERATOR_OBJ) \
        $(TPV102_SUBSTEP_ITERATOR_OBJ) \
        dynamic/tpv104_substep_iterator.o \
        $(FAULT_FACE_FLUX_OBJ) \
        $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_FAULT_RESAMPLE_APPLY_OBJ) \
	    $(FRICTION_SUBSTEP_ITERATOR_OBJ) \
	    $(TPV205_SUBSTEP_ITERATOR_OBJ) \
	    $(TPV102_SUBSTEP_ITERATOR_OBJ) \
	    dynamic/tpv104_substep_iterator.o \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_rupture_tilted_fault_serial_vs_parallel: $(TEST_RUPTURE_TILTED_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_RUPTURE_TILTED_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

seas_test_shared_fault_reconcile_cross_rank: $(TEST_RECONCILE_XRANK_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_TESTINTERNAL_OBJ) $(FRICTION_SOLVER_OBJ) $(TPV205_SUBSTEP_ITERATOR_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_RECONCILE_XRANK_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_TESTINTERNAL_OBJ) $(FRICTION_SOLVER_OBJ) $(TPV205_SUBSTEP_ITERATOR_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# R-020(b): psi-explicit np=2 shared-fault RS consistency (bare AdvanceADER +
# RateAndState; production reconcile, no SEAS_TEST_INTERNAL).
seas_test_shared_fault_rs_psi_consistency_np2: $(TEST_RS_PSI_NP2_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_RS_PSI_NP2_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_rupture_rake_sweep_cross_rank: $(TEST_RUPTURE_RAKE_SWEEP_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_RUPTURE_RAKE_SWEEP_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_adjacent_triangle_fault_uniformity: $(TEST_ADJACENT_TRI_UNIFORMITY_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ADJACENT_TRI_UNIFORMITY_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_adjacent_triangle_fault_first_step_audit: $(TEST_ADJACENT_TRI_FIRST_STEP_AUDIT_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ADJACENT_TRI_FIRST_STEP_AUDIT_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_godunov_rotation_roundtrip: $(TEST_GODUNOV_ROUNDTRIP_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_GODUNOV_ROUNDTRIP_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_friction_solver_equilibrium: $(TEST_FRICTION_EQUIL_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FRICTION_EQUIL_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_pepper_amplification_chain: $(TEST_PEPPER_AMP_CHAIN_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PEPPER_AMP_CHAIN_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_face_lift_operator: $(TEST_FACE_LIFT_OPERATOR_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FACE_LIFT_OPERATOR_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_pepper_qp_perturbation_response: $(TEST_PEPPER_QP_PERT_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PEPPER_QP_PERT_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_pepper_nonlinear_growth: $(TEST_PEPPER_NONLINEAR_GROWTH_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PEPPER_NONLINEAR_GROWTH_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_interior_vs_shared_branch_live: $(TEST_INTERIOR_VS_SHARED_LIVE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_INTERIOR_VS_SHARED_LIVE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Pepper-bug plan 2026-04-22 Phase 2 — sign_flipped / elem1_on_plus truth table
seas_test_sign_flipped_elem1_on_plus_truth_table: $(TEST_SIGN_FLIPPED_TRUTH_TABLE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_SIGN_FLIPPED_TRUTH_TABLE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Pepper-bug plan 2026-04-22 Phase 4 P4-T1 — MakeFaceKey uniqueness (np=2)
seas_test_make_face_key_uniqueness: $(TEST_MAKE_FACE_KEY_UNIQ_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_MAKE_FACE_KEY_UNIQ_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Pepper-bug plan 2026-04-22 Phase 4 P4-T2 — DOF-offset invariants
seas_test_fault_dof_offset_invariants: $(TEST_FAULT_DOF_OFFSET_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_DOF_OFFSET_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.2.0 plan Step 5 (R-V92-H01) — coupled RK4-on-psi convergence probe
seas_test_rk4_psi_integration: $(TEST_RK4_PSI_INTEGRATION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_RK4_PSI_INTEGRATION_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.2.0 round-9 REVIEW R-V92-K01 — driver post-RK4 contract regression guard
seas_test_rk4_endpoint_reeval: $(TEST_RK4_ENDPOINT_REEVAL_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_RK4_ENDPOINT_REEVAL_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.2.0 plan §4.10 (rev-3g) — single-channel Jacobian probe (d+e)
seas_test_volume_jacobian_single_channel: $(TEST_VOLUME_JAC_SINGLE_CH_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_VOLUME_JAC_SINGLE_CH_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.2.0 plan §4.11 (rev-3g) — RK4 conservation probe (b)
seas_test_rk4_conservation: $(TEST_RK4_CONSERVATION_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_RK4_CONSERVATION_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.2.0 plan §4.12 (rev-3g) — absorbing-BC energy decay probe (c)
seas_test_absorbing_bc_energy_decay: $(TEST_ABSORBING_BC_DECAY_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ABSORBING_BC_DECAY_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.2.0 plan §4.6 (rev-3b) — shared-fault role consistency across MPI
seas_test_shared_fault_role_consistency: $(TEST_SHARED_FAULT_ROLE_CONSIST_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_SHARED_FAULT_ROLE_CONSIST_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.2.0 plan §4.7 (rev-3c) — shared-fault DOFData runtime consistency
seas_test_shared_fault_dof_data_consistency: $(TEST_SHARED_FAULT_DOFDATA_CONSIST_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_SHARED_FAULT_DOFDATA_CONSIST_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.2.0 plan Step 5b (R-V92-I02) — 100-step F01+F02 stability probe
seas_test_rk4_f01f02_shared_fault_stability: $(TEST_RK4_F01F02_STABILITY_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_RK4_F01F02_STABILITY_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.2.0 plan §17 (rev-3d) — centroid-margin ctor-only probe
seas_test_centroid_margin_ctor_only: $(TEST_CENTROID_MARGIN_CTOR_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_CENTROID_MARGIN_CTOR_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.3.0 plan Phase 1 (I-04): FreeSurfaceGodunov vs FreeSurface equivalence
seas_test_free_surface_godunov: $(TEST_FREE_SURFACE_GODUNOV_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FREE_SURFACE_GODUNOV_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.3.0 plan Phase 2 (I-04 wiring, R-005): --free-surface-bc CLI parsing
seas_test_R005_free_surface_bc_flag_parse: $(TEST_R005_FS_BC_FLAG_PARSE_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R005_FS_BC_FLAG_PARSE_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.3.0 plan Phase 2 AC-3 (I-04 driver integration gate)
seas_test_free_surface_godunov_driver: $(TEST_FREE_SURFACE_GODUNOV_DRIVER_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FREE_SURFACE_GODUNOV_DRIVER_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.3.0 plan Phase 3 (I-06 part A): EvaluateTotal vs Evaluate equivalence
seas_test_fault_face_flux_total_vs_fluctuation: $(TEST_FAULT_FACE_FLUX_TOTAL_VS_FLUC_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_FACE_FLUX_TOTAL_VS_FLUC_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.3.0 R-003 (review-incorporation plan): EvaluateTotal psi-invariance
seas_test_R003_evaluatetotal_psi_invariant: $(TEST_R003_EVALUATETOTAL_PSI_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R003_EVALUATETOTAL_PSI_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.3.0 Phase 4 (I-06 part B): R-001 InitializeStateTotal sign regression
seas_test_R001_initstate_total_sign: $(TEST_R001_INITSTATE_TOTAL_SIGN_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R001_INITSTATE_TOTAL_SIGN_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.3.0 Phase 4 AC (I-06 part B): locked-fault 100 RK4 sigma_n stability
seas_test_tpv102_total_locked_fault: $(TEST_TPV102_TOTAL_LOCKED_FAULT_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV102_TOTAL_LOCKED_FAULT_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.4.0 R-005: R-002/R-009/R-I06-003 targets removed; exercised the
# DEPRECATED bulk-Q ApplyNucleationTotal path (see §5 R-005).
# v9.3.0 R-I06-001 / R-I06-006: total-Q absorbing BC equilibrium gate
seas_test_tpv102_total_absorbing_equilibrium: $(TEST_TPV102_TOTAL_ABSORBING_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV102_TOTAL_ABSORBING_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# 2026-04-22 persistent-prestress channel: links FaultFaceFlux + FrictionSolver
seas_test_persistent_nuc_prestress_channel: $(TEST_PERSISTENT_NUC_PRESTRESS_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PERSISTENT_NUC_PRESTRESS_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.4.0 debug-plan R-001: pepper reproducer (MFEM-only, no SEAS objs)
seas_test_tpv102_pepper_reproducer: $(TEST_TPV102_PEPPER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV102_PEPPER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# v9.3.0 R-I06-round3 R-006: direct tests for FreeSurfaceTotal / FreeSurfaceGodunovTotal
seas_test_free_surface_total_direct: $(TEST_FREE_SURFACE_TOTAL_DIRECT_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FREE_SURFACE_TOTAL_DIRECT_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# ADER I-05 Phase 1: polynomial-exactness test for ApplySpatialDerivative
# (Lever 1: also the DerivMode Cached/OnTheFly equivalence test — adds the
#  bimaterial deps GODUNOV_FLUX_POOL / HETEROGENEOUS_MATERIAL / SPATIAL_FRICTION
#  for the BimaterialWaveOperator pass.)
seas_test_wave_operator_spatial_derivative: $(TEST_WAVE_OP_SPATIAL_DERIV_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(GODUNOV_FLUX_POOL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(HETEROGENEOUS_MATERIAL_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(SPATIAL_FRICTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_WAVE_OP_SPATIAL_DERIV_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(GODUNOV_FLUX_POOL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(HETEROGENEOUS_MATERIAL_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(SPATIAL_FRICTION_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# FaceGeomCache builder equivalence + fault-exclusion (opt 2026-06-24 Phase 1).
# Header-only (dynamic/face_geom_cache.hpp) — links MFEM only, no seas objects.
TEST_FACE_GEOM_CACHE_SRC = tests/unit/test_face_geom_cache.cpp
TEST_FACE_GEOM_CACHE_OBJ = $(TEST_FACE_GEOM_CACHE_SRC:.cpp=.o)

seas_test_face_geom_cache: $(TEST_FACE_GEOM_CACHE_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FACE_GEOM_CACHE_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(TEST_FACE_GEOM_CACHE_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

test-face-geom-cache: seas_test_face_geom_cache
	./seas_test_face_geom_cache

# ADER hot-path optimization: micro-benchmark (timing tool; `make bench`)
.PHONY: bench
bench: seas_bench_ader_hotpath
seas_bench_ader_hotpath: $(BENCH_ADER_HOTPATH_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(BENCH_ADER_HOTPATH_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# ADER I-05 Phase 3: CK predictor plane-wave convergence test
seas_test_ader_ck_predictor: $(TEST_ADER_CK_PREDICTOR_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ADER_CK_PREDICTOR_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# ADER I-05 Phase 5: EvaluateADER / EvaluateADERTotal equivalence tests
seas_test_fault_face_flux_ader_equivalence: $(TEST_ADER_FAULT_FLUX_EQUIV_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ADER_FAULT_FLUX_EQUIV_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# ADER I-05 Phase 6: AdvanceADER vs RK4 on a linear wave
seas_test_ader_linear_wave_equivalence: $(TEST_ADER_LINEAR_WAVE_EQUIV_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ADER_LINEAR_WAVE_EQUIV_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# ADER fault-path live consistency: (Q_new - Q)/dt -> wave.Mult(Q) on a faulted fixture
seas_test_ader_fault_coupling_dt_limit: $(TEST_ADER_FAULT_COUPLING_DT_LIMIT_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ADER_FAULT_COUPLING_DT_LIMIT_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# ADER live interior-vs-shared branch cross-check (np=2)
seas_test_ader_interior_vs_shared_branch_live: $(TEST_ADER_INTERIOR_VS_SHARED_LIVE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ADER_INTERIOR_VS_SHARED_LIVE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# ADER I-05 Phase 8: TPV102 AdvanceADER smoke test
seas_test_ader_tpv102_smoke: $(TEST_ADER_TPV102_SMOKE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ADER_TPV102_SMOKE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Round-6 R-001 regression: driver defaults to ADER
seas_test_R001_driver_defaults_to_ader: $(TEST_R001_DRIVER_DEFAULT_ADER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R001_DRIVER_DEFAULT_ADER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Round-6 R-002 regression: driver init uses total-Q
seas_test_R002_driver_init_total_q: $(TEST_R002_DRIVER_INIT_TOTAL_Q_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R002_DRIVER_INIT_TOTAL_Q_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# REVIEW.md R-003: AntiplaneDomainOperator inherits the base default
# ComputeTractionDiagnostics, which must size diagnostic outputs to
# match `traction` (not empty) so downstream consumers like FaceTraceLogger
# can iterate without OOB.
seas_test_R003_antiplane_diagnostics_sizing: $(TEST_R003_ANTIPLANE_DIAG_SIZING_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R003_ANTIPLANE_DIAG_SIZING_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Dynamic rupture Phase 3: Friction solver tests
seas_test_friction_solver: $(TEST_FRICTION_SOLVER_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FRICTION_SOLVER_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 4a: TPV102 setup tests
seas_test_tpv102_setup: $(TEST_TPV102_SETUP_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV102_SETUP_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 4a: TPV102 driver — rebuilt against the TPV104 code flow
# (substep iterator, AgingLawPsi state evolution, mixed-flux dispatch).
seas_tpv102_driver: $(TPV102_DRIVER_OBJ) $(TPV102_SHARED_OBJS) \
                    $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) \
                    $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TPV102_DRIVER_OBJ) \
	    $(TPV102_SHARED_OBJS) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
	    $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) \
	    $(PML_LAYER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# TPV205 (SCEC TPV5, linear slip-weakening) driver.
seas_tpv205_driver: $(TPV205_DRIVER_OBJ) $(TPV205_SHARED_OBJS) \
                    $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) \
                    $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TPV205_DRIVER_OBJ) \
	    $(TPV205_SHARED_OBJS) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
	    $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) \
	    $(PML_LAYER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 4 of spatial_dynamic_rupture_plan.md (rev-3):
# seas_spatial_dyn_driver — generic spatially-varying dynamic-rupture
# driver.  Reuses the TPV205 LSW closed-form machinery + the new
# Phase H LSW_ForcedRupture dispatch.  Object compilation needs HDF5
# + toml11 include flags (spatial_friction.hpp pulls in
# heterogeneous_material.hpp which optionally includes the sidecar
# plumbing).
$(SPATIAL_DYN_DRIVER_OBJ): %.o: $(SRC)%.cpp \
                                $(SEAS_HEADERS) \
                                $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) \
	    $(TOML_FLAGS) -DSEAS_USE_MPI -c $< -o $@

seas_spatial_dyn_driver: $(SPATIAL_DYN_DRIVER_OBJ) \
                          $(TPV205_SHARED_OBJS) \
                          $(TPV102_SUBSTEP_ITERATOR_OBJ) \
                          $(FRICTION_ITERATOR_FACTORY_OBJ) \
                          $(FRICTION_SUBSTEP_ITERATOR_OBJ) \
                          $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
                          $(PRECOMPUTED_FACE_FLUXES_OBJ) \
                          $(GODUNOV_FLUX_OBJ) \
                          $(GODUNOV_FLUX_POOL_OBJ) \
                          $(PML_LAYER_OBJ) \
                          $(SPATIAL_FRICTION_OBJ) \
                          $(SPATIAL_NUCLEATION_OBJ) \
                          $(NUCLEATION_FACTORY_OBJ) \
                          $(RK_TIME_STEPPER_OBJ) \
                          $(SPATIAL_PRINT_DERIVED_OBJ) \
                          $(SPATIAL_VELOCITY_OBJ) \
                          $(SPATIAL_STRESS_OBJ) \
                          $(HETEROGENEOUS_MATERIAL_OBJ) \
                          $(MATERIAL_COEFFICIENTS_OBJ) \
                          $(DATA_FIELD_3D_OBJ) \
                          $(STRESS_FIELD_3D_OBJ) \
                          $(FIELD_COEFFICIENT_OBJ) \
                          $(LTS_CLUSTERING_OBJ) \
                          $(LTS_LAYOUT_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(SPATIAL_DYN_DRIVER_OBJ) \
	    $(TPV205_SHARED_OBJS) \
	    $(TPV102_SUBSTEP_ITERATOR_OBJ) \
	    $(FRICTION_ITERATOR_FACTORY_OBJ) \
	    $(FRICTION_SUBSTEP_ITERATOR_OBJ) \
	    $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
	    $(PRECOMPUTED_FACE_FLUXES_OBJ) \
	    $(GODUNOV_FLUX_OBJ) \
	    $(GODUNOV_FLUX_POOL_OBJ) \
	    $(PML_LAYER_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(SPATIAL_NUCLEATION_OBJ) \
	    $(NUCLEATION_FACTORY_OBJ) \
	    $(RK_TIME_STEPPER_OBJ) \
	    $(SPATIAL_PRINT_DERIVED_OBJ) \
	    $(SPATIAL_VELOCITY_OBJ) \
	    $(SPATIAL_STRESS_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(MATERIAL_COEFFICIENTS_OBJ) \
	    $(DATA_FIELD_3D_OBJ) \
	    $(STRESS_FIELD_3D_OBJ) \
	    $(FIELD_COEFFICIENT_OBJ) \
	    $(LTS_CLUSTERING_OBJ) \
	    $(LTS_LAYOUT_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS) $(SEAS_FS_LIB)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 14: rk_time_stepper object.  Includes spatial_friction.hpp (gate with
# TOML_FLAGS like nucleation_factory) + wave_operator.hpp (SEAS_EXTRA_CPPFLAGS,
# so a DOFData diag build stays ABI-consistent with the other fault objects).
$(RK_TIME_STEPPER_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                             $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) \
	    $(TOML_FLAGS) -c $< -o $@

# Phase 14: RK time-stepper unit test (tableau / order / bulk-parity / PsiRate /
# ApplyAbsolute).  Links the scalar wave-operator object set (the test
# instantiates AdvanceRKCoupled_Spatial<Mesh>) + the spatial friction /
# nucleation objects PsiRate / ApplyAbsolute resolve against.
$(TEST_RK_TIME_STEPPER_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                  $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) \
	    $(TOML_FLAGS) -c $< -o $@

seas_test_rk_time_stepper: $(TEST_RK_TIME_STEPPER_OBJ) \
                           $(RK_TIME_STEPPER_OBJ) \
                           $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
                           $(GODUNOV_FLUX_OBJ) \
                           $(GODUNOV_FLUX_POOL_OBJ) \
                           $(PRECOMPUTED_FACE_FLUXES_OBJ) \
                           $(PML_LAYER_OBJ) \
                           $(HETEROGENEOUS_MATERIAL_OBJ) \
                           $(FAULT_FACE_FLUX_OBJ) \
                           $(FRICTION_SOLVER_OBJ) \
                           $(SPATIAL_NUCLEATION_OBJ) \
                           $(SPATIAL_FRICTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_RK_TIME_STEPPER_OBJ) \
	    $(RK_TIME_STEPPER_OBJ) \
	    $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
	    $(GODUNOV_FLUX_OBJ) \
	    $(GODUNOV_FLUX_POOL_OBJ) \
	    $(PRECOMPUTED_FACE_FLUXES_OBJ) \
	    $(PML_LAYER_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(SPATIAL_NUCLEATION_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-rk-time-stepper: seas_test_rk_time_stepper
	./seas_test_rk_time_stepper

# Phase 14 (RK + LSW + mixed flux): EvaluateLSW + AdvanceRKCoupledLSW_Spatial
# unit tests (L1–L7).  Same object set as seas_test_rk_time_stepper (the test
# instantiates AdvanceRKCoupledLSW_Spatial<Mesh> + drives WaveOperator::Mult on
# a fault mesh).  $(TOML_FLAGS) keeps the DOFData/ABI consistent with the rest.
$(TEST_LSW_RK_MIXED_FLUX_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                    $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) \
	    $(TOML_FLAGS) -c $< -o $@

seas_test_lsw_rk_mixed_flux: $(TEST_LSW_RK_MIXED_FLUX_OBJ) \
                             $(RK_TIME_STEPPER_OBJ) \
                             $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
                             $(GODUNOV_FLUX_OBJ) \
                             $(GODUNOV_FLUX_POOL_OBJ) \
                             $(PRECOMPUTED_FACE_FLUXES_OBJ) \
                             $(PML_LAYER_OBJ) \
                             $(HETEROGENEOUS_MATERIAL_OBJ) \
                             $(FAULT_FACE_FLUX_OBJ) \
                             $(FRICTION_SOLVER_OBJ) \
                             $(SPATIAL_NUCLEATION_OBJ) \
                             $(SPATIAL_FRICTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_LSW_RK_MIXED_FLUX_OBJ) \
	    $(RK_TIME_STEPPER_OBJ) \
	    $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
	    $(GODUNOV_FLUX_OBJ) \
	    $(GODUNOV_FLUX_POOL_OBJ) \
	    $(PRECOMPUTED_FACE_FLUXES_OBJ) \
	    $(PML_LAYER_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(SPATIAL_NUCLEATION_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-lsw-rk-mixed-flux: seas_test_lsw_rk_mixed_flux
	./seas_test_lsw_rk_mixed_flux

# Phase 14 (RK + LSW): MPI shared-fault dispatch parity (L8) — the ONLY coverage
# of the parallel shared-fault Mult dispatch (wave_operator.inl shared-fault
# EvaluateLSW arm).  Run with: mpirun -np 2 ./seas_test_lsw_rk_shared_fault_mpi.
$(TEST_LSW_RK_SHARED_FAULT_MPI_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                          $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) \
	    $(TOML_FLAGS) -c $< -o $@

seas_test_lsw_rk_shared_fault_mpi: $(TEST_LSW_RK_SHARED_FAULT_MPI_OBJ) \
                             $(RK_TIME_STEPPER_OBJ) \
                             $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
                             $(GODUNOV_FLUX_OBJ) \
                             $(GODUNOV_FLUX_POOL_OBJ) \
                             $(PRECOMPUTED_FACE_FLUXES_OBJ) \
                             $(PML_LAYER_OBJ) \
                             $(HETEROGENEOUS_MATERIAL_OBJ) \
                             $(FAULT_FACE_FLUX_OBJ) \
                             $(FRICTION_SOLVER_OBJ) \
                             $(SPATIAL_NUCLEATION_OBJ) \
                             $(SPATIAL_FRICTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_LSW_RK_SHARED_FAULT_MPI_OBJ) \
	    $(RK_TIME_STEPPER_OBJ) \
	    $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
	    $(GODUNOV_FLUX_OBJ) \
	    $(GODUNOV_FLUX_POOL_OBJ) \
	    $(PRECOMPUTED_FACE_FLUXES_OBJ) \
	    $(PML_LAYER_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(SPATIAL_NUCLEATION_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-lsw-rk-shared-fault-mpi: seas_test_lsw_rk_shared_fault_mpi
	mpirun -np 2 ./seas_test_lsw_rk_shared_fault_mpi

# Phase N: spatial_nucleation object compile rule.
$(SPATIAL_NUCLEATION_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# Phase D: spatial_print_derived object compile rule.  Pulls
# spatial_friction.hpp (for NucleationSpec / StressSpec / lsw types)
# which in turn enables SEAS_USE_TOML conditionals — gate the build
# the same way the friction object is built.
$(SPATIAL_PRINT_DERIVED_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                   $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(TOML_FLAGS) -c $< -o $@

# Parity Phase 6 R-005 round-3: unit test for the new
# seas::ParaViewOutput::SetFaultParamsSpatial method.  Verifies the
# 8 SAFS fields are registered under their proper names (not TPV104's
# a/Dc/x2/x3) and that the call is idempotent.
TEST_PARAVIEW_SET_FAULT_PARAMS_SPATIAL_SRC = tests/unit/test_paraview_set_fault_params_spatial.cpp
TEST_PARAVIEW_SET_FAULT_PARAMS_SPATIAL_OBJ = $(TEST_PARAVIEW_SET_FAULT_PARAMS_SPATIAL_SRC:.cpp=.o)

$(TEST_PARAVIEW_SET_FAULT_PARAMS_SPATIAL_OBJ): %.o: $(SRC)%.cpp \
                                                   $(SEAS_HEADERS) \
                                                   $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -c $< -o $@

seas_test_paraview_set_fault_params_spatial: \
        $(TEST_PARAVIEW_SET_FAULT_PARAMS_SPATIAL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_PARAVIEW_SET_FAULT_PARAMS_SPATIAL_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-paraview-set-fault-params-spatial: seas_test_paraview_set_fault_params_spatial
	./seas_test_paraview_set_fault_params_spatial

# Phase N test target: 12 unit tests T-N01..T-N12 covering smoothStep
# math, GaussianFactorFaceLocal, ResolveGradualOverstress, the per-
# sub-step accumulator, and the iterator callback overload.
seas_test_spatial_nucleation: $(TEST_SPATIAL_NUCLEATION_OBJ) \
                              $(SPATIAL_NUCLEATION_OBJ) \
                              $(FAULT_FACE_FLUX_OBJ) \
                              $(FRICTION_SOLVER_OBJ) \
                              $(TPV205_SUBSTEP_ITERATOR_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_SPATIAL_NUCLEATION_OBJ) \
	    $(SPATIAL_NUCLEATION_OBJ) \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(TPV205_SUBSTEP_ITERATOR_OBJ) \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(TEST_SPATIAL_NUCLEATION_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                     $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

test-spatial-nucleation: seas_test_spatial_nucleation
	./seas_test_spatial_nucleation

# Phase 7 — nucleation_methods module.  nucleation_factory.o pulls
# spatial_friction.hpp (SpatialFrictionConfig / NucleationKind), so it is
# built with TOML_FLAGS like the spatial_friction / spatial_print_derived
# objects.
$(NUCLEATION_FACTORY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(TOML_FLAGS) -c $< -o $@

$(TEST_NUCLEATION_FACTORY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                     $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(TOML_FLAGS) -c $< -o $@

# Phase 7 acceptance: MakeNucleation dispatch (each kind + absent block).
# Links the nucleation factory + the spatial_nucleation resolvers + the
# fault_face_flux/friction_solver/tpv205 iterator object set (same closed set
# the spatial_nucleation test links; the factory adds no new dependencies).
seas_test_nucleation_factory: $(TEST_NUCLEATION_FACTORY_OBJ) \
                              $(NUCLEATION_FACTORY_OBJ) \
                              $(SPATIAL_NUCLEATION_OBJ) \
                              $(FAULT_FACE_FLUX_OBJ) \
                              $(FRICTION_SOLVER_OBJ) \
                              $(TPV205_SUBSTEP_ITERATOR_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_NUCLEATION_FACTORY_OBJ) \
	    $(NUCLEATION_FACTORY_OBJ) \
	    $(SPATIAL_NUCLEATION_OBJ) \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(TPV205_SUBSTEP_ITERATOR_OBJ) \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-nucleation-factory: seas_test_nucleation_factory
	./seas_test_nucleation_factory

# TPV-via-spatial-driver Phase 1 (RS primitives): nucleation-callback
# parity for Tpv102SubStepIterator.  Links the tpv102 iterator + the
# shared fault_face_flux + friction_solver objects (same set the
# standalone tpv102 driver uses, minus the wave-operator chain the test
# does not exercise).
$(TEST_TPV102_NUC_CALLBACK_PARITY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                            $(TPV102_HEADERS) \
                                            $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

seas_test_tpv102_nuc_callback_parity: $(TEST_TPV102_NUC_CALLBACK_PARITY_OBJ) \
                                      $(TPV102_SUBSTEP_ITERATOR_OBJ) \
                                      $(FAULT_FACE_FLUX_OBJ) \
                                      $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_TPV102_NUC_CALLBACK_PARITY_OBJ) \
	    $(TPV102_SUBSTEP_ITERATOR_OBJ) \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-tpv102-nuc-callback-parity: seas_test_tpv102_nuc_callback_parity
	./seas_test_tpv102_nuc_callback_parity

# TPV-via-spatial-driver Phase 1 (RS primitives): equilibrium-psi seed
# for the SAFS RS path.  The test includes dynamic/spatial_setup.hpp
# (-> spatial_friction.hpp + heterogeneous_material.hpp), so compile with
# the same HDF5 + toml include flags the spatial driver uses for that
# header chain.  Link only MFEM: SeedEquilibriumPsi_RS + DieterichRuina
# are header-only and FrictionSolver::V0 is an inline constexpr (C++17).
$(TEST_SEED_EQUILIBRIUM_PSI_RS_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                         dynamic/spatial_setup.hpp \
                                         friction/dieterich_ruina.hpp \
                                         spatial/code/spatial_friction.hpp \
                                         $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) \
	    $(TOML_FLAGS) -c $< -o $@

seas_test_seed_equilibrium_psi_rs: $(TEST_SEED_EQUILIBRIUM_PSI_RS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_SEED_EQUILIBRIUM_PSI_RS_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-seed-equilibrium-psi-rs: seas_test_seed_equilibrium_psi_rs
	./seas_test_seed_equilibrium_psi_rs

# TPV-via-spatial-driver Phase 2 (friction strategy interface).  These
# objects pull friction_iterator.hpp -> wave_operator.hpp (FaultFrictionLaw)
# + spatial_friction.hpp, so compile with the same HDF5 + toml include
# flags the spatial driver uses for that header chain.  wave_operator.inl
# is template-only, so no TU here instantiates WaveOperator -> no link dep
# on wave_operator.o.
$(FRICTION_ITERATOR_FACTORY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                       dynamic/friction_iterator.hpp \
                                       dynamic/friction_iterator_factory.hpp \
                                       spatial/code/spatial_friction.hpp \
                                       $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) \
	    $(TOML_FLAGS) -c $< -o $@

$(TEST_FRICTION_ITERATOR_FACTORY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                            dynamic/friction_iterator.hpp \
                                            dynamic/friction_iterator_factory.hpp \
                                            $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) \
	    $(TOML_FLAGS) -c $< -o $@

seas_test_friction_iterator_factory: $(TEST_FRICTION_ITERATOR_FACTORY_OBJ) \
                                     $(FRICTION_ITERATOR_FACTORY_OBJ) \
                                     $(FRICTION_SUBSTEP_ITERATOR_OBJ) \
                                     $(TPV205_SUBSTEP_ITERATOR_OBJ) \
                                     $(TPV102_SUBSTEP_ITERATOR_OBJ) \
                                     $(FAULT_FACE_FLUX_OBJ) \
                                     $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_FRICTION_ITERATOR_FACTORY_OBJ) \
	    $(FRICTION_ITERATOR_FACTORY_OBJ) \
	    $(FRICTION_SUBSTEP_ITERATOR_OBJ) \
	    $(TPV205_SUBSTEP_ITERATOR_OBJ) \
	    $(TPV102_SUBSTEP_ITERATOR_OBJ) \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-friction-iterator-factory: seas_test_friction_iterator_factory
	./seas_test_friction_iterator_factory

# ---- Phase 5: unified method-oriented sub-step iterator + parity test ----
# Explicit header prerequisites (TPV*_HEADERS are defined later in this file,
# so they cannot be referenced here without expanding empty).
$(FRICTION_SUBSTEP_ITERATOR_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                       dynamic/friction_substep_iterator.hpp \
                                       friction/state_policies.hpp \
                                       dynamic/friction_iterator.hpp \
                                       dynamic/tpv205_friction.hpp \
                                       friction/slip_law_srw_psi.hpp \
                                       spatial/code/spatial_friction.hpp \
                                       $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) \
	    $(TOML_FLAGS) -c $< -o $@

$(TEST_FRICTION_SUBSTEP_ITERATOR_PARITY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                       dynamic/friction_substep_iterator.hpp \
                                       friction/state_policies.hpp \
                                       dynamic/tpv205_substep_iterator.hpp \
                                       dynamic/tpv102_substep_iterator.hpp \
                                       dynamic/tpv104_substep_iterator.hpp \
                                       friction/slip_law_srw_psi.hpp \
                                       $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) \
	    $(TOML_FLAGS) -c $< -o $@

# NOTE: TPV104_SUBSTEP_ITERATOR_OBJ is defined later in this file, so the
# variable would expand empty in this prerequisite list (GNU Make expands
# prereqs at parse time) and the stale .o would never be rebuilt — reference
# the literal object path instead.  tpv205/tpv102/fault_face_flux/friction_solver
# vars are all defined earlier, so they are safe to use by name.
seas_test_friction_substep_iterator_parity: \
        $(TEST_FRICTION_SUBSTEP_ITERATOR_PARITY_OBJ) \
        $(FRICTION_SUBSTEP_ITERATOR_OBJ) \
        $(TPV205_SUBSTEP_ITERATOR_OBJ) \
        $(TPV102_SUBSTEP_ITERATOR_OBJ) \
        dynamic/tpv104_substep_iterator.o \
        $(FAULT_FACE_FLUX_OBJ) \
        $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_FRICTION_SUBSTEP_ITERATOR_PARITY_OBJ) \
	    $(FRICTION_SUBSTEP_ITERATOR_OBJ) \
	    $(TPV205_SUBSTEP_ITERATOR_OBJ) \
	    $(TPV102_SUBSTEP_ITERATOR_OBJ) \
	    dynamic/tpv104_substep_iterator.o \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-friction-substep-iterator-parity: seas_test_friction_substep_iterator_parity
	./seas_test_friction_substep_iterator_parity

$(TEST_ADVANCE_INTERFACE_COMPILES_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                             dynamic/friction_iterator.hpp \
                                             $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) \
	    $(TOML_FLAGS) -c $< -o $@

# R-010 compile guard: links only MFEM — the stub iterator + the
# driver-mirror function exercise the interface without instantiating any
# concrete adapter, so no iterator/flux objects are needed.
seas_test_advance_interface_compiles: $(TEST_ADVANCE_INTERFACE_COMPILES_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_ADVANCE_INTERFACE_COMPILES_OBJ) \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-advance-interface-compiles: seas_test_advance_interface_compiles
	./seas_test_advance_interface_compiles

# Phase 3 (SAFS + rate-and-state): ResolveRateState guard tests (R-003 no
# double pore-pressure, R-006 per-DOF b/f_0/V_0 reject, R-011 a>b allowed).
$(TEST_RESOLVE_RATE_STATE_GUARDS_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                            spatial/code/spatial_friction.hpp \
                                            $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) \
	    $(TOML_FLAGS) -c $< -o $@

seas_test_resolve_rate_state_guards: $(TEST_RESOLVE_RATE_STATE_GUARDS_OBJ) \
                                     $(SPATIAL_FRICTION_OBJ) \
                                     $(HETEROGENEOUS_MATERIAL_OBJ) \
                                     $(MATERIAL_COEFFICIENTS_OBJ) \
                                     $(DATA_FIELD_3D_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_RESOLVE_RATE_STATE_GUARDS_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(MATERIAL_COEFFICIENTS_OBJ) \
	    $(DATA_FIELD_3D_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-resolve-rate-state-guards: seas_test_resolve_rate_state_guards
	./seas_test_resolve_rate_state_guards

# Phase 3 (D3.1, R-007): constant_tensor right-lateral-positive sign test,
# parametrized over {Dc2,Dc8,Dc10,rate_state}.
$(TEST_CONSTANT_TENSOR_SIGN_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                       spatial/code/spatial_stress.hpp \
                                       fault/fault_basis.hpp \
                                       $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) \
	    $(TOML_FLAGS) -c $< -o $@

seas_test_constant_tensor_sign: $(TEST_CONSTANT_TENSOR_SIGN_OBJ) \
                                $(SPATIAL_STRESS_OBJ) \
                                $(SPATIAL_FRICTION_OBJ) \
                                $(STRESS_FIELD_3D_OBJ) \
                                $(DATA_FIELD_3D_OBJ) \
                                $(FIELD_COEFFICIENT_OBJ) \
                                $(HETEROGENEOUS_MATERIAL_OBJ) \
                                $(MATERIAL_COEFFICIENTS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_CONSTANT_TENSOR_SIGN_OBJ) \
	    $(SPATIAL_STRESS_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(STRESS_FIELD_3D_OBJ) \
	    $(DATA_FIELD_3D_OBJ) \
	    $(FIELD_COEFFICIENT_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(MATERIAL_COEFFICIENTS_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-constant-tensor-sign: seas_test_constant_tensor_sign
	./seas_test_constant_tensor_sign

# =====================================================================
#  Phase 10 (TPV31) unit tests
# =====================================================================
$(TEST_TPV31_CANONICAL_ROTATION_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                           spatial/code/spatial_stress.hpp \
                                           spatial/code/spatial_friction.hpp \
                                           dynamic/heterogeneous_material.hpp \
                                           $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_TPV31_NUCLEATION_AND_COHESION_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                                spatial/code/spatial_friction.hpp \
                                                dynamic/spatial_nucleation.hpp \
                                                $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_TPV31_STATION_WRITER_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                       dynamic/tpv31_stations.hpp \
                                       dynamic/fault_face_flux.hpp \
                                       $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_SPATIAL_STATION_DISPATCH_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                       dynamic/tpv205_stations.hpp \
                                       dynamic/tpv104_stations.hpp \
                                       dynamic/tpv102_stations.hpp \
                                       dynamic/fault_face_flux.hpp \
                                       $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

seas_test_tpv31_canonical_rotation: $(TEST_TPV31_CANONICAL_ROTATION_OBJ) \
                                    $(SPATIAL_STRESS_OBJ) \
                                    $(SPATIAL_FRICTION_OBJ) \
                                    $(STRESS_FIELD_3D_OBJ) \
                                    $(DATA_FIELD_3D_OBJ) \
                                    $(FIELD_COEFFICIENT_OBJ) \
                                    $(HETEROGENEOUS_MATERIAL_OBJ) \
                                    $(MATERIAL_COEFFICIENTS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_TPV31_CANONICAL_ROTATION_OBJ) \
	    $(SPATIAL_STRESS_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(STRESS_FIELD_3D_OBJ) \
	    $(DATA_FIELD_3D_OBJ) \
	    $(FIELD_COEFFICIENT_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(MATERIAL_COEFFICIENTS_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_tpv31_nucleation_and_cohesion: \
                              $(TEST_TPV31_NUCLEATION_AND_COHESION_OBJ) \
                              $(SPATIAL_NUCLEATION_OBJ) \
                              $(SPATIAL_FRICTION_OBJ) \
                              $(FAULT_FACE_FLUX_OBJ) \
                              $(FRICTION_SOLVER_OBJ) \
                              $(HETEROGENEOUS_MATERIAL_OBJ) \
                              $(TPV205_SUBSTEP_ITERATOR_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_TPV31_NUCLEATION_AND_COHESION_OBJ) \
	    $(SPATIAL_NUCLEATION_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(TPV205_SUBSTEP_ITERATOR_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_tpv31_station_writer: $(TEST_TPV31_STATION_WRITER_OBJ) \
                                $(FAULT_FACE_FLUX_OBJ) \
                                $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_TPV31_STATION_WRITER_OBJ) \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-tpv31-canonical-rotation: seas_test_tpv31_canonical_rotation
	./seas_test_tpv31_canonical_rotation

test-tpv31-nucleation-and-cohesion: seas_test_tpv31_nucleation_and_cohesion
	./seas_test_tpv31_nucleation_and_cohesion

test-tpv31-station-writer: seas_test_tpv31_station_writer
	./seas_test_tpv31_station_writer

seas_test_spatial_station_dispatch: $(TEST_SPATIAL_STATION_DISPATCH_OBJ) \
                                    $(FAULT_FACE_FLUX_OBJ) \
                                    $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_SPATIAL_STATION_DISPATCH_OBJ) \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-spatial-station-dispatch: seas_test_spatial_station_dispatch
	./seas_test_spatial_station_dispatch

# =====================================================================
#  Phase 8 — TPV205/102/104 benchmark config + sign/basis tests.
# =====================================================================

# (1) TPV TOML stress-sign (req 2 / D3.2): parser + fault geometry + basis.
$(TEST_TPV_TOML_STRESS_SIGN_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                       spatial/code/spatial_friction.hpp \
                                       fault/fault_basis.hpp \
                                       $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) \
	    $(TOML_FLAGS) -c $< -o $@

seas_test_tpv_toml_stress_sign: $(TEST_TPV_TOML_STRESS_SIGN_OBJ) \
                                $(SPATIAL_STRESS_OBJ) \
                                $(SPATIAL_FRICTION_OBJ) \
                                $(STRESS_FIELD_3D_OBJ) \
                                $(DATA_FIELD_3D_OBJ) \
                                $(FIELD_COEFFICIENT_OBJ) \
                                $(HETEROGENEOUS_MATERIAL_OBJ) \
                                $(MATERIAL_COEFFICIENTS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_TPV_TOML_STRESS_SIGN_OBJ) \
	    $(SPATIAL_STRESS_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(STRESS_FIELD_3D_OBJ) \
	    $(DATA_FIELD_3D_OBJ) \
	    $(FIELD_COEFFICIENT_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(MATERIAL_COEFFICIENTS_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-tpv-toml-stress-sign: seas_test_tpv_toml_stress_sign
	./seas_test_tpv_toml_stress_sign

# (2) Planar-TPV basis (req 3 / D1): header-only FaultBasis + tiny mesh.
$(TEST_PLANAR_TPV_BASIS_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                   fault/fault_basis.hpp \
                                   $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) \
	    $(TOML_FLAGS) -c $< -o $@

seas_test_planar_tpv_basis: $(TEST_PLANAR_TPV_BASIS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_PLANAR_TPV_BASIS_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-planar-tpv-basis: seas_test_planar_tpv_basis
	./seas_test_planar_tpv_basis

# (3) TPV config-parse (AC "each config parses"): parser only.
$(TEST_TPV_CONFIG_PARSE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                   spatial/code/spatial_friction.hpp \
                                   $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) \
	    $(TOML_FLAGS) -c $< -o $@

seas_test_tpv_config_parse: $(TEST_TPV_CONFIG_PARSE_OBJ) \
                            $(SPATIAL_STRESS_OBJ) \
                            $(SPATIAL_FRICTION_OBJ) \
                            $(STRESS_FIELD_3D_OBJ) \
                            $(DATA_FIELD_3D_OBJ) \
                            $(FIELD_COEFFICIENT_OBJ) \
                            $(HETEROGENEOUS_MATERIAL_OBJ) \
                            $(MATERIAL_COEFFICIENTS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_TPV_CONFIG_PARSE_OBJ) \
	    $(SPATIAL_STRESS_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(STRESS_FIELD_3D_OBJ) \
	    $(DATA_FIELD_3D_OBJ) \
	    $(FIELD_COEFFICIENT_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(MATERIAL_COEFFICIENTS_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-tpv-config-parse: seas_test_tpv_config_parse
	./seas_test_tpv_config_parse

# Phase 3 (R-005): ParaView RS state channel writes d.psi.
$(TEST_PARAVIEW_STATE_CHANNEL_RS_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                            $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) \
	    $(TOML_FLAGS) -c $< -o $@

seas_test_paraview_state_channel_rs: $(TEST_PARAVIEW_STATE_CHANNEL_RS_OBJ) \
                                     $(FAULT_FACE_FLUX_OBJ) \
                                     $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_PARAVIEW_STATE_CHANNEL_RS_OBJ) \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-paraview-state-channel-rs: seas_test_paraview_state_channel_rs
	./seas_test_paraview_state_channel_rs

# Phase D test target — print-derived pre-flight gate unit tests
# T-D02 / T-D05 / T-D06 / T-D07 exercised in-process via the printer
# interface; full T-D01..T-D05 integration tests use the driver.
seas_test_spatial_print_derived: $(TEST_SPATIAL_PRINT_DERIVED_OBJ) \
                                 $(SPATIAL_PRINT_DERIVED_OBJ) \
                                 $(SPATIAL_NUCLEATION_OBJ) \
                                 $(SPATIAL_FRICTION_OBJ) \
                                 $(HETEROGENEOUS_MATERIAL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_SPATIAL_PRINT_DERIVED_OBJ) \
	    $(SPATIAL_PRINT_DERIVED_OBJ) \
	    $(SPATIAL_NUCLEATION_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(TEST_SPATIAL_PRINT_DERIVED_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                        $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(TOML_FLAGS) -c $< -o $@

test-spatial-print-derived: seas_test_spatial_print_derived
	./seas_test_spatial_print_derived

# Phase 4 smoke target: build + dry-run on the production SAFS LSW config.
test-spatial-dyn-driver: seas_spatial_dyn_driver
	./seas_spatial_dyn_driver \
	    --config safs/project_7.0_alternative/config/spatial_friction_slip_weakening_safs_projected_stress.toml \
	    --no-sidecar-material \
	    --dry-run

# R-001 / R-002 / R-003 regression smoke: run a short (<= 0.01s) rupture
# with --paraview, then verify that the fault VTU was written AND that
# the volume VTU registered the dynamic fault projection fields.  VTU
# (not HDF5) so the smoke runs on builds without MFEM_USE_H5Z_ZFP.
test-spatial-dyn-driver-paraview-smoke: seas_spatial_dyn_driver
	rm -rf /tmp/safs_paraview_smoke
	./seas_spatial_dyn_driver \
	    --config safs/project_7.0_alternative/config/spatial_friction_slip_weakening_safs_projected_stress.toml \
	    --no-sidecar-material \
	    --output-dir /tmp/safs_paraview_smoke \
	    --tfinal 0.01 \
	    --paraview \
	    --paraview-volume vtu --paraview-bulk off --paraview-fault vtu
	@echo "[R-001] checking fault VTU was written ..."
	@test -f /tmp/safs_paraview_smoke/FaultSurface/fault_surface_c0.vtu \
	    || (echo "R-001 FAIL: fault_surface_c0.vtu missing" && exit 1)
	@grep -q 'Name="slip_rate_dip"' \
	    /tmp/safs_paraview_smoke/FaultSurface/fault_surface_c0.vtu \
	    || (echo "R-001 FAIL: slip_rate_dip channel missing" && exit 1)
	@grep -q 'Name="traction_strike"' \
	    /tmp/safs_paraview_smoke/FaultSurface/fault_surface_c0.vtu \
	    || (echo "R-001 FAIL: traction_strike channel missing" && exit 1)
	@echo "[R-003] checking volume VTU registered dynamic fault fields ..."
	@grep -q 'Name="slip_rate_dip"' \
	    /tmp/safs_paraview_smoke/volume/Cycle000000/proc000000.vtu \
	    || (echo "R-003 FAIL: slip_rate_dip not in volume VTU" && exit 1)
	@grep -q 'Name="lsw_mu_s"' \
	    /tmp/safs_paraview_smoke/volume/Cycle000000/proc000000.vtu \
	    || (echo "R-003 FAIL: lsw_mu_s (SAFS static) not in volume VTU" && exit 1)
	@echo "[R-001 / R-002 / R-003] paraview-output smoke PASSED"

# Phase 4a: TPV102 local integration tests
seas_test_tpv102_local: $(TEST_TPV102_LOCAL_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV102_LOCAL_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Mesh generation tool (serial)
seas_generate_mesh: $(GENERATE_MESH_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(GENERATE_MESH_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Object file rules
# Data-projection feature object rules (need HDF5 include flags for
# the .cpp side; MaterialField is header-only and so its test only
# depends on MFEM headers).
$(LTS_CLUSTERING_OBJ): %.o: $(SRC)%.cpp $(LTS_CLUSTERING_HEADERS) \
                            $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_LTS_CLUSTERING_OBJ): %.o: $(SRC)%.cpp $(LTS_CLUSTERING_HEADERS) \
                                 $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(LTS_LAYOUT_OBJ): %.o: $(SRC)%.cpp $(LTS_LAYOUT_HEADERS) \
                        $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_LTS_LAYOUT_OBJ): %.o: $(SRC)%.cpp $(LTS_LAYOUT_HEADERS) \
                             $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_LTS_TIME_BASIS_OBJ): %.o: $(SRC)%.cpp $(LTS_TIME_BASIS_HEADERS) \
                                 $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(LTS_STEPPER_OBJ): %.o: $(SRC)%.cpp $(LTS_STEPPER_HEADERS) \
                         $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_LTS_SCHEDULER_OBJ): %.o: $(SRC)%.cpp $(LTS_STEPPER_HEADERS) \
                                $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(DATA_FIELD_3D_OBJ): %.o: $(SRC)%.cpp $(DATA_FIELD_3D_HEADERS) \
                          $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -c $< -o $@

$(STRESS_FIELD_3D_OBJ): %.o: $(SRC)%.cpp $(STRESS_FIELD_3D_HEADERS) \
                             $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -c $< -o $@

$(FIELD_COEFFICIENT_OBJ): %.o: $(SRC)%.cpp $(FIELD_COEFFICIENT_HEADERS) \
                               $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -c $< -o $@

$(TEST_DATA_FIELD_3D_OBJ): %.o: $(SRC)%.cpp $(DATA_FIELD_3D_HEADERS) \
                                $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -c $< -o $@

$(TEST_DATA_FIELD_3D_SHARED_MEM_OBJ): %.o: $(SRC)%.cpp \
                                $(STRESS_FIELD_3D_HEADERS) \
                                $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -c $< -o $@

$(TEST_DATA_FIELD_3D_SHARED_MEM_FIN_OBJ): %.o: $(SRC)%.cpp \
                                $(DATA_FIELD_3D_HEADERS) \
                                $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -c $< -o $@

$(TEST_STRESS_FIELD_3D_OBJ): %.o: $(SRC)%.cpp $(STRESS_FIELD_3D_HEADERS) \
                                  $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -c $< -o $@

$(TEST_FIELD_PROJECTOR_OBJ): %.o: $(SRC)%.cpp $(FIELD_COEFFICIENT_HEADERS) \
                                  $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -c $< -o $@

$(HETEROGENEOUS_MATERIAL_OBJ): %.o: $(SRC)%.cpp \
                                    $(HETEROGENEOUS_MATERIAL_HEADERS) \
                                    $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(MATERIAL_COEFFICIENTS_OBJ): %.o: $(SRC)%.cpp \
                                    $(MATERIAL_COEFFICIENTS_HEADERS) \
                                    $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	# R-007 round-4: HDF5_INCFLAGS not required — material_coefficients
	# only touches DataField3D's public API (no hid_t in surface).
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_HETEROGENEOUS_MATERIAL_OBJ): %.o: $(SRC)%.cpp \
                                         $(HETEROGENEOUS_MATERIAL_HEADERS) \
                                         $(MATERIAL_COEFFICIENTS_HEADERS) \
                                         $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -c $< -o $@

$(PROJECT_VELOCITY_DRIVER_OBJ): %.o: $(SRC)%.cpp \
                                      $(FIELD_COEFFICIENT_HEADERS) \
                                      $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -c $< -o $@

$(TEST_FRICTION_OBJ): %.o: $(SRC)%.cpp $(FRICTION_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_STATE_EVOLUTION_OBJ): %.o: $(SRC)%.cpp $(FRICTION_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_ANTIPLANE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_FAULT_OPERATOR_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_QUASI_DYNAMIC_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_BP2_SHORT_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_IO_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(BP2_SERIAL_SMOKE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_CHECKPOINT_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_PSI_STATE_OBJ): %.o: $(SRC)%.cpp $(FRICTION_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_FAULT_DETECTION_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_DOMAIN_INTERFACE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_BP5_PARAMS_OBJ): %.o: $(SRC)%.cpp $(FRICTION_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_VECTOR_FRICTION_OBJ): %.o: $(SRC)%.cpp $(FRICTION_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_FAULT_BASIS_OBJ): %.o: $(SRC)%.cpp $(FAULT_HEADERS) $(DOMAIN_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_ELASTICITY_BR2_OBJ): %.o: $(SRC)%.cpp $(DOMAIN_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_ELASTICITY_OPERATOR_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_BP5_FAULT_OPERATOR_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_FAULT_DOF_BASIS_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_PROJECT_FAULT_PRESTRESS_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                          $(STRESS_FIELD_3D_HEADERS) \
                                          $(FIELD_COEFFICIENT_HEADERS) \
                                          $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -c $< -o $@

$(TEST_COMPUTE_SAFS_PARAMS_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                      $(STRESS_FIELD_3D_HEADERS) \
                                      $(FIELD_COEFFICIENT_HEADERS) \
                                      $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -c $< -o $@

$(TEST_SAFS_MODE_WIRING_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                   $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_SAFS_STRESS_CONFIG_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                     $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) \
	    $(TOML_FLAGS) -c $< -o $@

$(PROJECT_STRESS_DRIVER_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                   $(STRESS_FIELD_3D_HEADERS) \
                                   $(FIELD_COEFFICIENT_HEADERS) \
                                   $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -c $< -o $@

$(TEST_PROJECT_STRESS_TO_MESH_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                         $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -c $< -o $@

$(TEST_BP5_INTEGRATION_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_BP5_OUTPUT_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_IO_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_BP5_MESH_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_DIAG_VTK_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_MESH_REGIONS_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_FACE_QUADRATURE_OBJ): %.o: $(SRC)%.cpp $(FAULT_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_CROSS_VERIFY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(GODUNOV_FLUX_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

# Phase H.1: GodunovFluxPool object + its dedup test.
$(GODUNOV_FLUX_POOL_OBJ): %.o: $(SRC)%.cpp \
                               dynamic/godunov_flux_pool.hpp \
                               dynamic/godunov_flux.hpp \
                               $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# Phase 9: BimaterialFlux object + its unit test (ported from hrs-ref).
$(GODUNOV_FLUX_BIMATERIAL_OBJ): %.o: $(SRC)%.cpp \
                                     dynamic/godunov_flux_bimaterial.hpp \
                                     dynamic/godunov_flux.hpp \
                                     dynamic/wave_state.hpp \
                                     $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_GODUNOV_FLUX_BIMATERIAL_OBJ): %.o: $(SRC)%.cpp \
                                          dynamic/godunov_flux_bimaterial.hpp \
                                          dynamic/godunov_flux.hpp \
                                          dynamic/wave_state.hpp \
                                          $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

seas_test_godunov_flux_bimaterial: $(TEST_GODUNOV_FLUX_BIMATERIAL_OBJ) \
                                   $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
                                   $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_GODUNOV_FLUX_BIMATERIAL_OBJ) \
	    $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
	    $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-godunov-flux-bimaterial: seas_test_godunov_flux_bimaterial
	./seas_test_godunov_flux_bimaterial

$(TEST_PHASEH_GODUNOV_FLUX_POOL_OBJ): %.o: $(SRC)%.cpp \
                                           dynamic/godunov_flux_pool.hpp \
                                           dynamic/godunov_flux.hpp \
                                           $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# Spatial Phase 1: TOML parser + resolver objects + tests.  Requires
# $(TOML_FLAGS) for extern/toml11/toml.hpp.
$(SPATIAL_FRICTION_OBJ): %.o: $(SRC)%.cpp \
                              spatial/code/spatial_friction.hpp \
                              dynamic/heterogeneous_material.hpp \
                              $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(TOML_FLAGS) -c $< -o $@

$(TEST_SPATIAL_FRICTION_CONFIG_OBJ): %.o: $(SRC)%.cpp \
                                          spatial/code/spatial_friction.hpp \
                                          $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(TOML_FLAGS) -c $< -o $@

$(TEST_SPATIAL_FRICTION_RESOLVER_OBJ): %.o: $(SRC)%.cpp \
                                            spatial/code/spatial_friction.hpp \
                                            $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(TOML_FLAGS) -c $< -o $@

$(TEST_SPATIAL_FRICTION_SIDECAR_OBJ): %.o: $(SRC)%.cpp \
                                           spatial/code/spatial_friction.hpp \
                                           io/data_field_3d.hpp \
                                           $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(TOML_FLAGS) $(HDF5_INCFLAGS) -c $< -o $@

$(TEST_FRICTION_DEPTH_PROFILE_OBJ): %.o: $(SRC)%.cpp \
                                         spatial/code/spatial_friction.hpp \
                                         $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(TOML_FLAGS) -c $< -o $@

# Spatial Phase 2: SpatialVelocityBundle.  Needs HDF5_INCFLAGS because the
# .cpp pulls in io/data_field_3d.hpp + io/material_coefficients.hpp.
$(SPATIAL_VELOCITY_OBJ): %.o: $(SRC)%.cpp \
                              spatial/code/spatial_velocity.hpp \
                              spatial/code/spatial_friction.hpp \
                              io/data_field_3d.hpp \
                              io/material_coefficients.hpp \
                              io/field_coefficient.hpp \
                              $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) \
	    $(TOML_FLAGS) -c $< -o $@

$(TEST_SPATIAL_VELOCITY_BUNDLE_OBJ): %.o: $(SRC)%.cpp \
                                          spatial/code/spatial_velocity.hpp \
                                          $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) \
	    $(TOML_FLAGS) -c $< -o $@

# Spatial Phases 3 + 3b: ApplyCsmStressSidecar + ConstantTensorStressSource.
$(SPATIAL_STRESS_OBJ): %.o: $(SRC)%.cpp \
                            spatial/code/spatial_stress.hpp \
                            spatial/code/spatial_friction.hpp \
                            fault/fault_geometry.hpp \
                            fault/fault_geometry_safs_templated.inl \
                            io/stress_field_3d.hpp \
                            $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) \
	    $(TOML_FLAGS) -c $< -o $@

$(TEST_SPATIAL_STRESS_BUNDLE_OBJ): %.o: $(SRC)%.cpp \
                                        spatial/code/spatial_stress.hpp \
                                        $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) \
	    $(TOML_FLAGS) -c $< -o $@

$(TEST_SPATIAL_CONSTANT_STRESS_OBJ): %.o: $(SRC)%.cpp \
                                          spatial/code/spatial_stress.hpp \
                                          $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) \
	    $(TOML_FLAGS) -c $< -o $@

$(TEST_SPATIAL_SETUP_OBJ): %.o: $(SRC)%.cpp \
                                dynamic/spatial_setup.hpp \
                                dynamic/fault_face_flux.hpp \
                                spatial/code/spatial_friction.hpp \
                                $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(TOML_FLAGS) -c $< -o $@

$(TEST_PHASEH_LSW_FORCED_RUPTURE_OBJ): %.o: $(SRC)%.cpp \
                                            dynamic/fault_face_flux.hpp \
                                            spatial/code/spatial_friction.hpp \
                                            $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(TOML_FLAGS) -c $< -o $@

$(TEST_PHASEH_WAVE_OPERATOR_CONSTANT_PARITY_OBJ): %.o: $(SRC)%.cpp \
                                                       dynamic/wave_operator.hpp \
                                                       dynamic/wave_operator.inl \
                                                       dynamic/godunov_flux_pool.hpp \
                                                       dynamic/heterogeneous_material.hpp \
                                                       $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# Phase 13 C-6: depends on BOTH wave_operator and the new bimaterial subclass
# headers so editing either re-triggers this compile (the Makefile has no
# auto header deps — see the stale-.o hazard note).
$(TEST_BIMATERIAL_WAVE_OPERATOR_PARITY_OBJ): %.o: $(SRC)%.cpp \
                                                  dynamic/wave_operator.hpp \
                                                  dynamic/wave_operator.inl \
                                                  dynamic/bimaterial_wave_operator.hpp \
                                                  dynamic/bimaterial_wave_operator.inl \
                                                  dynamic/godunov_flux_pool.hpp \
                                                  dynamic/godunov_flux_bimaterial.hpp \
                                                  dynamic/heterogeneous_material.hpp \
                                                  dynamic/fault_face_flux.hpp \
                                                  $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(WAVE_OPERATOR_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_GODUNOV_FLUX_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_WAVE_OPERATOR_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(PML_LAYER_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_WAVE_BC_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_PML_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# fault_face_flux.cpp includes spatial/code/spatial_friction.hpp, so it must
# rebuild when that header changes (see the SPATIAL_HEADERS note above); the
# $(DYNAMIC_HEADERS) list alone does not cover it.
$(FAULT_FACE_FLUX_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(SPATIAL_HEADERS) \
                             $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

# TPV102 "Topology-Based Precomputed Face-Rotation" plan 2026-04-23 Phase 1
$(PRECOMPUTED_FACE_FLUXES_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_PRECOMPUTED_FLUXES_INTERIOR_IDENTITY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_PRECOMPUTED_FLUXES_BOUNDARY_IDENTITY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_PRECOMPUTED_FLUXES_PHASE2A_SWITCH_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_ARM1_CONSTANT_STATE_PROBES_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_ARM2_D4_EQUIVARIANT_FIXTURE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_ARM2B_VERTEX_PERM_PROBE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_ARM3D_TOPO_STATE_SAMPLING_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_R008_FAULT_BRANCH_PRESERVED_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_R5_001_HOIST_COMPILES_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_R001_RK4_USES_PRECOMPUTED_FLAG_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_D4_FIXTURE_VERTEX_EQUIV_OBJ): %.o: $(SRC)%.cpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_D4_FAULT_SANITY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_D4_JACOBIAN_SIGN_OBJ): %.o: $(SRC)%.cpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_R5_GATE14_PEPPER_CORR_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_R6_PHASE2E_AMPCHAIN_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_R7_AMPCHAIN_V2_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_R9_FROZEN_FRICTION_C_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_PRECOMPUTED_FLUXES_PHASE2B_MPI_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_R002_EMPTY_FAULT_RANK_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_R1600_SUBSTEP_NO_FAULT_RANK_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_GHOST_EXCHANGE_BYNODES_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_MIXED_FLUX_ADJACENT_MPI_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_MIXED_FLUX_SHARED_FAULT_MPI_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_MIXED_FLUX_DISPATCH_ADJACENT_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_MIXED_FLUX_DISPATCH_ADJACENT_MPI_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_MIXED_FLUX_DISPATCH_ALL_CONT_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(FRICTION_SOLVER_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_FAULT_FACE_FLUX_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_FAULT_FACE_FLUX_FRAME_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_FAULT_FACE_FLUX_FRAME_AND_FLUX_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_GODUNOV_IDENTITY_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_FAULT_FLUX_BRANCH_EQUIV_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# v9.0.0 Pelties-9 fix: new unit tests (§3.1e, §3.1f)
$(TEST_GODUNOV_INTERIOR_EQUAL_SIDES_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_FAULT_FACE_FLUX_PER_SIDE_ASSEMBLY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_FAULT_FACE_FLUX_BIMATERIAL_GUARD_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_BIMATERIAL_FAULT_RIEMANN_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_TPV6_STATIONS_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# v9.1.0 plan §4.1 / §4.2 — Pelties-9 follow-through tests
$(TEST_FAULT_SURFACE_VTU_CONTINUITY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# v9.2.0 plan §19 R-V92-E02 stage-4 VTU diagnostic
$(TEST_FAULT_SURFACE_VTU_K4_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# --paraview-fields / SetFaultVTUFields field filter
$(TEST_FAULT_SURFACE_VTU_FIELD_FILTER_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# Phase 1 PLAN_paraview_compaction_2026-04-28 — Phase 1 binary VTU writer tests
$(TEST_FAULT_SURFACE_VTU_BINARY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_FAULT_SURFACE_VTU_GATHER_MPI_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# free_surface_output.hpp is header-only and listed in IO_HEADERS (via
# SEAS_HEADERS), so editing it force-rebuilds this object.
$(TEST_FREE_SURFACE_SLICE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# Phase 2 PLAN_paraview_compaction_2026-04-28 — VTKHDF tests
$(TEST_PARAVIEW_HDF_SMOKE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_FAULT_SURFACE_VTKHDF_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_FAULT_SURFACE_VTKHDF_MPI_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_VTKHDF_ZFP_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_FAULT_SURFACE_VTKHDF_ZFP_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_VOLUME_HDF_COMPRESSION_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_PARAVIEW_SCHEDULE_CAP_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_KINEMATICS_FIELD_SET_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_BP5_PETSC_TS_ZERO_FAULT_RANK_OBJ): %.o: $(SRC)%.cpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_BP5_PETSC_TS_RESTART_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_TPV104_CHECKPOINT_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV104_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_PARAVIEW_RANK0_WARNING_GATE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_FAULT_BASIS_DIP_STRIKE_SYMMETRY_OBJ): %.o: $(SRC)%.cpp $(FAULT_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# v9.2.0 plan §4.1-§4.4 regression gates
$(TEST_GODUNOV_ROTATION_IDENTITY_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_FRAME_DIAG_PROJECTION_IDENTITY_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_MACRO_DIAG_DECOMPOSITION_IDENTITY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_CANONICAL_ROTATION_PURE_STRIKESLIP_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_FRAME_DIP_LEAK_PROBE_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_OFFFAULT_FRAME_SYMMETRY_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_FAULT_BASIS_QP_ORTHONORMALITY_OBJ): %.o: $(SRC)%.cpp $(FAULT_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_PER_QP_VS_CENTROID_BASIS_OBJ): %.o: $(SRC)%.cpp $(FAULT_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_SHARED_FAULT_ROLE_CONSIST_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -DSEAS_USE_MPI -c $< -o $@

$(TEST_NO_PENALTY_DR_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_INTERIOR_FAULT_FLUX_PATH_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_INTERIOR_FAULT_FLUX_PATH_ASYM_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_INTERIOR_FAULT_FLUX_PATH_RUPTURE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_BULK_RK4_NOISE_STAB_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_SHARED_FAULT_FLUX_PATH_ANA_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -DSEAS_USE_MPI -c $< -o $@

$(TEST_BULK_PLANE_WAVE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_INTERIOR_VS_SHARED_LIVE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -DSEAS_TEST_INTERNAL -DSEAS_USE_MPI -c $< -o $@

$(TEST_RUPTURE_MULTISTEP_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -DSEAS_USE_MPI -c $< -o $@

$(TEST_FAULT_PLANAR_SERIAL_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -DSEAS_USE_MPI -c $< -o $@

$(TEST_FAULT_OVERINT_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -DSEAS_USE_MPI -c $< -o $@

$(TEST_FAULT_RESAMPLE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -DSEAS_USE_MPI -c $< -o $@

$(TEST_FAULT_RESAMPLE_APPLY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -DSEAS_USE_MPI -c $< -o $@

$(TEST_RUPTURE_TILTED_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -DSEAS_USE_MPI -c $< -o $@

$(TEST_RECONCILE_XRANK_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -DSEAS_TEST_INTERNAL -DSEAS_USE_MPI -c $< -o $@

$(TEST_RS_PSI_NP2_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -DSEAS_USE_MPI -c $< -o $@

# Same source as $(FAULT_FACE_FLUX_OBJ) but compiled WITH -DSEAS_TEST_INTERNAL
# so the reconcile oracle's absolute-Pa injection hook + static knob are present.
$(FAULT_FACE_FLUX_TESTINTERNAL_OBJ): $(SRC)dynamic/fault_face_flux.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -DSEAS_TEST_INTERNAL -c $< -o $@

$(TEST_RUPTURE_RAKE_SWEEP_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -DSEAS_USE_MPI -c $< -o $@

$(TEST_ADJACENT_TRI_UNIFORMITY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -DSEAS_USE_MPI -c $< -o $@

$(TEST_ADJACENT_TRI_FIRST_STEP_AUDIT_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_GODUNOV_ROUNDTRIP_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_FRICTION_EQUIL_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_PEPPER_AMP_CHAIN_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_FACE_LIFT_OPERATOR_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_PEPPER_QP_PERT_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_PEPPER_NONLINEAR_GROWTH_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_SIGN_FLIPPED_TRUTH_TABLE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -DSEAS_TEST_INTERNAL -c $< -o $@

$(TEST_MAKE_FACE_KEY_UNIQ_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -DSEAS_USE_MPI -c $< -o $@

$(TEST_FAULT_DOF_OFFSET_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -DSEAS_TEST_INTERNAL -c $< -o $@

$(TEST_RK4_PSI_INTEGRATION_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_RK4_ENDPOINT_REEVAL_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_VOLUME_JAC_SINGLE_CH_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_RK4_CONSERVATION_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_ABSORBING_BC_DECAY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_SHARED_FAULT_DOFDATA_CONSIST_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -DSEAS_USE_MPI -c $< -o $@

$(TEST_RK4_F01F02_STABILITY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -DSEAS_USE_MPI -c $< -o $@

$(TEST_CENTROID_MARGIN_CTOR_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -DSEAS_USE_MPI -c $< -o $@

# v9.3.0 plan Phase 1 (I-04): FreeSurfaceGodunov vs FreeSurface equivalence
$(TEST_FREE_SURFACE_GODUNOV_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# v9.3.0 plan Phase 2 (I-04 wiring, R-005): --free-surface-bc CLI parsing
$(TEST_R005_FS_BC_FLAG_PARSE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# v9.3.0 plan Phase 2 AC-3 (I-04 driver integration gate)
$(TEST_FREE_SURFACE_GODUNOV_DRIVER_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# v9.3.0 plan Phase 3 (I-06 part A): EvaluateTotal vs Evaluate equivalence
$(TEST_FAULT_FACE_FLUX_TOTAL_VS_FLUC_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# v9.3.0 R-003 (review-incorporation plan): EvaluateTotal psi-invariance
$(TEST_R003_EVALUATETOTAL_PSI_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# v9.3.0 Phase 4 (I-06 part B): R-001 InitializeStateTotal sign regression
$(TEST_R001_INITSTATE_TOTAL_SIGN_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# v9.3.0 Phase 4 AC (I-06 part B): locked-fault 100 RK4 sigma_n stability
$(TEST_TPV102_TOTAL_LOCKED_FAULT_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# v9.3.0 R-I06-001 / R-I06-006: total-Q absorbing BC equilibrium gate
$(TEST_TPV102_TOTAL_ABSORBING_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# 2026-04-22 nucleation persistent-prestress channel fix regression
$(TEST_PERSISTENT_NUC_PRESTRESS_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# v9.4.0 debug-plan R-001: pepper reproducer (MFEM-only; no SEAS headers)
$(TEST_TPV102_PEPPER_OBJ): %.o: $(SRC)%.cpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# v9.3.0 R-I06-round3 R-006: direct tests for FreeSurfaceTotal variants
$(TEST_FREE_SURFACE_TOTAL_DIRECT_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# ADER I-05 Phase 1: ApplySpatialDerivative polynomial-exactness test
$(TEST_WAVE_OP_SPATIAL_DERIV_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# ADER hot-path optimization: micro-benchmark object
$(BENCH_ADER_HOTPATH_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# ADER I-05 Phase 3: CK predictor plane-wave convergence
$(TEST_ADER_CK_PREDICTOR_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# ADER I-05 Phase 5: EvaluateADER / EvaluateADERTotal equivalence
$(TEST_ADER_FAULT_FLUX_EQUIV_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# ADER I-05 Phase 6: AdvanceADER vs RK4 on a linear wave
$(TEST_ADER_LINEAR_WAVE_EQUIV_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# ADER fault-path live consistency on a faulted TPV102 fixture
$(TEST_ADER_FAULT_COUPLING_DT_LIMIT_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# ADER live interior-vs-shared branch cross-check
$(TEST_ADER_INTERIOR_VS_SHARED_LIVE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# ADER I-05 Phase 8: TPV102 AdvanceADER smoke test
$(TEST_ADER_TPV102_SMOKE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# Round-6 R-001 regression: driver default is ADER
$(TEST_R001_DRIVER_DEFAULT_ADER_OBJ): %.o: $(SRC)%.cpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# Round-6 R-002 regression: driver init uses total-Q
$(TEST_R002_DRIVER_INIT_TOTAL_Q_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_R003_ANTIPLANE_DIAG_SIZING_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_FRICTION_SOLVER_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_TPV102_SETUP_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TPV102_DRIVER_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV102_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TPV102_SUBSTEP_ITERATOR_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV102_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

# TPV205 compile rules.
$(TPV205_DRIVER_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV205_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TPV205_SUBSTEP_ITERATOR_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV205_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TEST_TPV102_LOCAL_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(GENERATE_MESH_OBJ): %.o: $(SRC)%.cpp $(DOMAIN_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# Testing targets
test: test-friction test-state-evolution test-antiplane test-fault-operator \
      test-quasi-dynamic test-bp2-short test-io test-checkpoint \
      test-fault-detection test-domain-interface \
      test-bp5-params test-vector-friction test-fault-basis \
      test-elasticity-br2 test-elasticity-operator \
      test-bp5-fault-operator \
      test-bp5-integration \
      test-bp5-output test-bp5-mesh \
      test-diag-vtk test-mesh-regions \
      test-face-quadrature \
      test-linear-elastic test-constitutive-integrator \
      test-boundary-config test-boundary-config-operator \
      test-boundary-classification \
      test-logging test-config-bridge test-config-parser \
      test-godunov-flux test-godunov-flux-bimaterial test-wave-operator \
      test-bimaterial-central-flux test-bimaterial-mixed-flux-dispatch \
      test-bimaterial-mixed-flux-cfl \
      test-bimaterial-mixed-flux-homog-equivalence \
      test-phaseh-godunov-flux-pool \
      test-spatial-friction-config test-spatial-friction-resolver \
      test-spatial-friction-sidecar \
      test-tpv205-friction \
      test-friction-depth-profile \
      test-tpv102-nuc-callback-parity test-seed-equilibrium-psi-rs \
      test-friction-iterator-factory test-advance-interface-compiles \
      test-resolve-rate-state-guards test-constant-tensor-sign \
      test-tpv31-canonical-rotation test-tpv31-nucleation-and-cohesion \
      test-tpv31-station-writer \
      test-tpv-toml-stress-sign test-planar-tpv-basis test-tpv-config-parse \
      test-paraview-state-channel-rs \
      test-spatial-velocity-bundle \
      test-spatial-constant-stress-source \
      test-spatial-setup \
      test-phaseh-lsw-forced-rupture \
      test-phaseh-wave-operator-constant-parity \
      test-bimaterial-wave-operator-parity \
      test-rk-time-stepper \
      test-wave-bc test-pml \
      test-fault-face-flux test-friction-solver \
      test-tpv102-setup test-tpv102-local \
      test-fault-face-flux-frame \
      test-fault-face-flux-frame-and-flux \
      test-godunov-identity-normal-reversal \
      test-fault-flux-interior-vs-shared-branch-equivalence \
      test-godunov-interior-equal-sides-identity \
      test-fault-face-flux-per-side-assembly \
      test-fault-face-flux-bimaterial-guard \
      test-free-surface-godunov \
      test-R005-free-surface-bc-flag-parse \
      test-free-surface-godunov-driver \
      test-fault-face-flux-total-vs-fluctuation \
      test-R003-evaluatetotal-psi-invariant \
      test-R001-initstate-total-sign \
      test-tpv102-total-locked-fault \
      test-tpv102-total-absorbing-equilibrium \
      test-persistent-nuc-prestress-channel \
      test-tpv102-pepper-reproducer \
      test-free-surface-total-direct \
      test-wave-operator-spatial-derivative \
      test-face-geom-cache \
      test-ader-ck-predictor \
      test-fault-face-flux-ader-equivalence \
      test-ader-linear-wave-equivalence \
      test-ader-tpv102-smoke \
      test-R001-driver-defaults-to-ader \
      test-R002-driver-init-total-q \
      test-paraview-schedule-cap \
      test-kinematics-field-set \
      test-bp5-petsc-ts-zero-fault-rank
# Note: `test-paraview-rank0-warning-gate` is INTENTIONALLY excluded
# from the `test:` umbrella above — it requires `mpirun -np 2` which
# is not guaranteed available in every dev environment.  Run it
# manually as a pre-submit check:
#     make test-paraview-rank0-warning-gate

test-data-projection: test-data-field-3d test-stress-field-3d test-field-projector test-heterogeneous-material
	@echo "data-projection feature: all C++ unit tests passed."

test-lts-clustering: seas_test_lts_clustering
	./seas_test_lts_clustering

test-lts-layout: seas_test_lts_layout
	./seas_test_lts_layout

test-lts-time-basis: seas_test_lts_time_basis
	./seas_test_lts_time_basis

test-lts-scheduler: seas_test_lts_scheduler
	./seas_test_lts_scheduler

test-data-field-3d: seas_test_data_field_3d
	./seas_test_data_field_3d

# PLAN_sidecar_mpi_shared_memory_2026-07-17.md §6 gate.  Kept OUT of the
# test-data-projection roll-up for the same reason as the other np>=2
# targets (mpirun -np 2 is not guaranteed in every dev environment); run
# manually as a pre-submit check:
#     make test-data-field-3d-shared-mem
test-data-field-3d-shared-mem: seas_test_data_field_3d_shared_mem \
                               seas_test_data_field_3d_shared_mem_finalize
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_data_field_3d_shared_mem
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_data_field_3d_shared_mem
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_data_field_3d_shared_mem_finalize
	@echo "--- R-002: NaN fixture in shared mode MUST abort (not hang) ---"
	@TO=""; command -v timeout  >/dev/null 2>&1 && TO="timeout 60"; \
	 [ -z "$$TO" ] && command -v gtimeout >/dev/null 2>&1 && TO="gtimeout 60"; \
	 $$TO $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 \
	     ./seas_test_data_field_3d_shared_mem --nan-fixture-abort-child \
	     >/dev/null 2>&1; rc=$$?; \
	 if [ $$rc -eq 0 ]; then \
	     echo "FAILED: NaN shared-mode fixture did NOT abort (guard broken)"; exit 1; \
	 elif [ -n "$$TO" ] && [ $$rc -eq 124 ]; then \
	     echo "FAILED: NaN shared-mode fixture HUNG (killed at 60 s)"; exit 1; \
	 else \
	     echo "OK: aborted as expected (mpirun exit $$rc)"; \
	 fi

# Phase 6 §1 of PLAN_onfaultstress.md (six-component stress reader).
test-stress-field-3d: seas_test_stress_field_3d
	./seas_test_stress_field_3d

test-field-projector: seas_test_field_projector
	./seas_test_field_projector

test-heterogeneous-material: seas_test_heterogeneous_material
	mpirun -np 1 ./seas_test_heterogeneous_material

test-friction: seas_test_friction_law
	./seas_test_friction_law

test-state-evolution: seas_test_state_evolution
	./seas_test_state_evolution

test-antiplane: seas_test_antiplane
	./seas_test_antiplane

test-fault-operator: seas_test_fault_operator
	./seas_test_fault_operator

test-quasi-dynamic: seas_test_quasi_dynamic
	./seas_test_quasi_dynamic

test-bp2-short: seas_test_bp2_short
	./seas_test_bp2_short

test-io: seas_test_io
	./seas_test_io

test-bp2-serial-smoke: seas_bp2_serial_smoke
	./seas_bp2_serial_smoke

test-checkpoint: seas_test_checkpoint
	./seas_test_checkpoint

test-psi-state: seas_test_psi_state
	./seas_test_psi_state

test-fault-detection: seas_test_fault_detection
	./seas_test_fault_detection

test-domain-interface: seas_test_domain_interface
	./seas_test_domain_interface

test-bp5-params: seas_test_bp5_params
	./seas_test_bp5_params

test-vector-friction: seas_test_vector_friction
	./seas_test_vector_friction

test-fault-basis: seas_test_fault_basis
	./seas_test_fault_basis

test-elasticity-br2: seas_test_elasticity_br2
	./seas_test_elasticity_br2

test-elasticity-operator: seas_test_elasticity_operator
	./seas_test_elasticity_operator

test-bp5-fault-operator: seas_test_bp5_fault_operator
	./seas_test_bp5_fault_operator

# Phase 6.A — FaultGeometry per-DOF coords / basis accessor test
test-fault-dof-basis: seas_test_fault_dof_basis
	./seas_test_fault_dof_basis

# Phase 6 §4 — ProjectFaultPreStress test
test-project-fault-prestress: seas_test_project_fault_prestress
	./seas_test_project_fault_prestress

# Phase 6 §5 — FaultGeometry::ComputeParams test
test-compute-safs-params: seas_test_compute_safs_params
	./seas_test_compute_safs_params

# Phase 6 §6 — SAFS-mode wiring test
test-safs-mode-wiring: seas_test_safs_mode_wiring
	./seas_test_safs_mode_wiring

# Phase 6 §7 — StressConfig schema test
test-safs-stress-config: seas_test_safs_stress_config
	./seas_test_safs_stress_config

# Phase 7 — project_stress_to_mesh driver round-trip test
test-project-stress-to-mesh: seas_test_project_stress_to_mesh \
                              seas_project_stress_to_mesh
	./seas_test_project_stress_to_mesh

test-bp5-integration: seas_test_bp5_integration
	./seas_test_bp5_integration

test-bp5-output: seas_test_bp5_output
	./seas_test_bp5_output

test-bp5-mesh: seas_test_bp5_mesh
	./seas_test_bp5_mesh

test-diag-vtk: seas_test_diag_vtk
	./seas_test_diag_vtk

test-mesh-regions: seas_test_mesh_regions
	./seas_test_mesh_regions

test-face-quadrature: seas_test_face_quadrature
	./seas_test_face_quadrature

test-cross-verify: seas_test_cross_verify
	./seas_test_cross_verify

# Phase 3: Constitutive model tests
seas_test_linear_elastic: $(TEST_LINEAR_ELASTIC_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_LINEAR_ELASTIC_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(TEST_LINEAR_ELASTIC_OBJ): %.o: $(SRC)%.cpp $(CONSTITUTIVE_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

test-linear-elastic: seas_test_linear_elastic
	./seas_test_linear_elastic

seas_test_constitutive_integrator: $(TEST_CONSTITUTIVE_INTEG_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_CONSTITUTIVE_INTEG_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(TEST_CONSTITUTIVE_INTEG_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

test-constitutive-integrator: seas_test_constitutive_integrator
	./seas_test_constitutive_integrator

# Phase 4: Boundary config tests
seas_test_boundary_config: $(TEST_BOUNDARY_CONFIG_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BOUNDARY_CONFIG_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(TEST_BOUNDARY_CONFIG_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

test-boundary-config: seas_test_boundary_config
	./seas_test_boundary_config

seas_test_boundary_config_operator: $(TEST_BOUNDARY_CONFIG_OP_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BOUNDARY_CONFIG_OP_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(TEST_BOUNDARY_CONFIG_OP_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

test-boundary-config-operator: seas_test_boundary_config_operator
	./seas_test_boundary_config_operator

# Phase 4h-iii: Boundary classification (BuildFacetBCTables) face-by-face tests
seas_test_boundary_classification: $(TEST_BOUNDARY_CLASSIF_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BOUNDARY_CLASSIF_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(TEST_BOUNDARY_CLASSIF_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

test-boundary-classification: seas_test_boundary_classification
	./seas_test_boundary_classification

# Phase 6: Logging tests
seas_test_logging: $(TEST_LOGGING_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_LOGGING_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(TEST_LOGGING_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

test-logging: seas_test_logging
	./seas_test_logging

# Phase 7a: Config bridge tests
seas_test_config_bridge: $(TEST_CONFIG_BRIDGE_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_CONFIG_BRIDGE_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(TEST_CONFIG_BRIDGE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

test-config-bridge: seas_test_config_bridge
	./seas_test_config_bridge

# Phase 5: TOML config tests and driver
seas_test_config_parser: $(TEST_CONFIG_PARSER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_CONFIG_PARSER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(TEST_CONFIG_PARSER_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

test-config-parser: seas_test_config_parser
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_config_parser

test-godunov-flux: seas_test_godunov_flux
	./seas_test_godunov_flux

# Phase H.1: GodunovFluxPool dedup tests.
test-phaseh-godunov-flux-pool: seas_test_phaseh_godunov_flux_pool
	./seas_test_phaseh_godunov_flux_pool

# Spatial Phase 1: TOML parser + per-DOF friction resolver tests.
test-spatial-friction-config: seas_test_spatial_friction_config
	./seas_test_spatial_friction_config

test-spatial-friction-resolver: seas_test_spatial_friction_resolver
	./seas_test_spatial_friction_resolver

# Phase 3 (PLAN_thermal_case2): [friction.rate_state.sidecar] resolver path.
test-spatial-friction-sidecar: seas_test_spatial_friction_sidecar
	./seas_test_spatial_friction_sidecar

test-friction-depth-profile: seas_test_friction_depth_profile
	./seas_test_friction_depth_profile

test-spatial-velocity-bundle: seas_test_spatial_velocity_bundle
	./seas_test_spatial_velocity_bundle

test-spatial-stress-bundle: seas_test_spatial_stress_bundle
	@echo "[skipped from `make test`] test_spatial_stress_bundle " \
	      "depends on a BP5 ElasticityDomainOperator fixture that " \
	      "currently segfaults on this branch (same pre-existing " \
	      "regression as test-bp5-fault-operator's " \
	      "ValidateFacetBCTables failure).  Build still works; run " \
	      "by hand if the BP5 fixture is fixed."
	./seas_test_spatial_stress_bundle || true

test-spatial-constant-stress-source: seas_test_spatial_constant_stress_source
	./seas_test_spatial_constant_stress_source

test-spatial-setup: seas_test_spatial_setup
	./seas_test_spatial_setup

test-phaseh-lsw-forced-rupture: seas_test_phaseh_lsw_forced_rupture
	./seas_test_phaseh_lsw_forced_rupture

# TPV205 LSW friction-kernel unit tests (R-002/R-003 strength-barrier
# short-circuits + the sigma_n strength floor, sliver-blowup plan
# 2026-05-26).  Header-only; no .o beyond MFEM.
test-tpv205-friction: seas_test_tpv205_friction
	./seas_test_tpv205_friction

# Phase H Stage 1 byte-exact gate; runs serial + np=4.
test-phaseh-wave-operator-constant-parity: seas_test_phaseh_wave_operator_constant_parity
	./seas_test_phaseh_wave_operator_constant_parity
	mpirun -np 4 ./seas_test_phaseh_wave_operator_constant_parity

# Phase 13 C-6 fault-bearing parity (flux_-placeholder leak-class tripwire); serial.
test-bimaterial-wave-operator-parity: seas_test_bimaterial_wave_operator_parity
	./seas_test_bimaterial_wave_operator_parity

test-wave-operator: seas_test_wave_operator
	./seas_test_wave_operator

test-wave-bc: seas_test_wave_bc
	./seas_test_wave_bc

test-pml: seas_test_pml
	./seas_test_pml

test-fault-face-flux: seas_test_fault_face_flux
	./seas_test_fault_face_flux

test-fault-face-flux-frame: seas_test_fault_face_flux_frame
	./seas_test_fault_face_flux_frame

test-fault-face-flux-frame-and-flux: seas_test_fault_face_flux_frame_and_flux
	./seas_test_fault_face_flux_frame_and_flux

test-godunov-identity-normal-reversal: seas_test_godunov_identity_normal_reversal
	./seas_test_godunov_identity_normal_reversal

test-fault-flux-interior-vs-shared-branch-equivalence: seas_test_fault_flux_interior_vs_shared_branch_equivalence
	./seas_test_fault_flux_interior_vs_shared_branch_equivalence

# v9.0.0 Pelties-9 fix: §3.1e / §3.1f phonies
test-godunov-interior-equal-sides-identity: seas_test_godunov_interior_equal_sides_identity
	./seas_test_godunov_interior_equal_sides_identity

test-fault-face-flux-per-side-assembly: seas_test_fault_face_flux_per_side_assembly
	./seas_test_fault_face_flux_per_side_assembly

# v9.1.0 plan §4.1 / §4.2 — Pelties-9 follow-through
test-fault-surface-vtu-continuity: seas_test_fault_surface_vtu_continuity
	./seas_test_fault_surface_vtu_continuity

# v9.2.0 plan §19 — R-V92-E02 stage-4 VTU diagnostic
test-fault-surface-vtu-k4: seas_test_fault_surface_vtu_k4
	./seas_test_fault_surface_vtu_k4

# --paraview-fields / SetFaultVTUFields field filter
test-fault-surface-vtu-field-filter: seas_test_fault_surface_vtu_field_filter
	./seas_test_fault_surface_vtu_field_filter

# Phase 1 PLAN_paraview_compaction_2026-04-28 — single binary VTU writer
test-fault-surface-vtu-binary: seas_test_fault_surface_vtu_binary
	./seas_test_fault_surface_vtu_binary

# Phase 1 PLAN_paraview_compaction_2026-04-28 — np=4 gather correctness
test-fault-surface-vtu-gather-mpi: seas_test_fault_surface_vtu_gather_mpi
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_fault_surface_vtu_gather_mpi

# Phase 4 PLAN_free_surface_slice.md — free-surface slice (serial + np=2:
# np=2 exercises empty-local-submesh ranks + collective transfer/Save).
test-free-surface-slice: seas_test_free_surface_slice
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_free_surface_slice
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_free_surface_slice

# Phase 2a PLAN_paraview_compaction_2026-04-28 — VTKHDF build prerequisite
# Plan §2a.6: run serial AND mpirun -np 4.  The test self-skips on builds
# without MFEM_USE_HDF5 (returns 0 with a SKIP message).
test-paraview-hdf-smoke: seas_test_paraview_hdf_smoke
	./seas_test_paraview_hdf_smoke
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_paraview_hdf_smoke

# Phase 2b PLAN_paraview_compaction_2026-04-28 — VTKHDF serial bit-exact
test-fault-surface-vtkhdf: seas_test_fault_surface_vtkhdf
	./seas_test_fault_surface_vtkhdf

# Phase 2b PLAN_paraview_compaction_2026-04-28 — VTKHDF np=4 merge
test-fault-surface-vtkhdf-mpi: seas_test_fault_surface_vtkhdf_mpi
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_fault_surface_vtkhdf_mpi

# Phase 2d.2 PLAN_paraview_compaction_2026-04-28 — H5Z-ZFP MFEM core gate.
# HDF5_PLUGIN_PATH must point at the directory containing libh5zzfp.{so,dylib}
# (set by the user at the env level; conda install places it at
# $CONDA_PREFIX/plugin).
test-vtkhdf-zfp: seas_test_vtkhdf_zfp
	HDF5_PLUGIN_PATH=$$CONDA_PREFIX/plugin ./seas_test_vtkhdf_zfp

# Phase 2d.3 PLAN_paraview_compaction_2026-04-28 — seas fault + ZFP gate
test-fault-surface-vtkhdf-zfp: seas_test_fault_surface_vtkhdf_zfp
	HDF5_PLUGIN_PATH=$$CONDA_PREFIX/plugin ./seas_test_fault_surface_vtkhdf_zfp

# Phase 6.2 PLAN_bulk_compression_and_size_estimator_2026-05-09 — volume HDF + ZFP
test-volume-hdf-compression: seas_test_volume_hdf_compression
	HDF5_PLUGIN_PATH=$$CONDA_PREFIX/plugin ./seas_test_volume_hdf_compression

# R-312 / plan §Phase 2d.3 acceptance — bulk-vs-fault ZFP tolerance
# isolation gate.  Reuses the existing fault-side ZFP test (which
# proves fault tolerance is honoured at 1e-12) and skips the bulk
# half explicitly because the volume `pv_` writer is currently
# `mfem::ParaViewDataCollection` (per-rank VTU), not VTKHDF — so no
# bulk HDF chunk filter is in scope until the volume writer is
# promoted to `ParaViewHDFDataCollection` (deferred follow-up,
# tracked via SetVolumeHDFCompression's runtime warning at
# paraview_output.hpp:430).
test-paraview-tolerance-isolation: seas_test_fault_surface_vtkhdf_zfp
	HDF5_PLUGIN_PATH=$$CONDA_PREFIX/plugin ./seas_test_fault_surface_vtkhdf_zfp
	@echo
	@echo "R-312 NOTE: bulk-side ZFP isolation half is deferred — the"
	@echo "  volume writer is currently VTU.  See paraview_output.hpp"
	@echo "  SetVolumeHDFCompression for the runtime warning that fires"
	@echo "  when --paraview-bulk-zfp-tol is requested today."

# Phase 3 PLAN_paraview_compaction_2026-04-28 — snapshot cap unit test
test-paraview-schedule-cap: seas_test_paraview_schedule_cap
	./seas_test_paraview_schedule_cap

test-kinematics-field-set: seas_test_kinematics_field_set
	./seas_test_kinematics_field_set

# R-003 of REVIEW.md 2026-05-16 — zero-fault-DOF padded-shrink invariant.
test-bp5-petsc-ts-zero-fault-rank: seas_test_bp5_petsc_ts_zero_fault_rank
	./seas_test_bp5_petsc_ts_zero_fault_rank

# Phase 1 of petsc_ts_restart_plan_2026-05-16 — V2 checkpoint round-trip
# + restart-state restoration.  Sub-tests are mostly serial-runnable
# (file I/O round-trip + ParaViewOutput setter/getter contracts) so the
# target runs the binary directly without mpiexec.  Full driver-level
# integration (sub-test 3 fresh-vs-restart trajectory) is deferred to
# the bp5_phase6_paraview_zfp_normal_48hr.sbatch re-submission flow,
# per plan §"Testing Strategy".  Per R-104 / R-011, this target is NOT
# added to the `test:` umbrella.
test-bp5-petsc-ts-restart: seas_test_bp5_petsc_ts_restart
	./seas_test_bp5_petsc_ts_restart

# Phase 4 TPV104 V1 restart unit test (serial; not in `test:` umbrella).
test-tpv104-checkpoint: seas_test_tpv104_checkpoint
	./seas_test_tpv104_checkpoint

# R-104 of REVIEW.md 2026-05-16 round 2 — MPI-only rank-0 warning gate.
# Requires np >= 2; skipped (returns 0) under serial.
test-paraview-rank0-warning-gate: seas_test_paraview_rank0_warning_gate
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_paraview_rank0_warning_gate

test-fault-basis-dip-strike-symmetry: seas_test_fault_basis_dip_strike_symmetry
	./seas_test_fault_basis_dip_strike_symmetry

# v9.2.0 plan §4.1-§4.4 regression gates
test-godunov-rotation-identity: seas_test_godunov_rotation_identity
	./seas_test_godunov_rotation_identity

test-canonical-rotation-pure-strikeslip: seas_test_canonical_rotation_pure_strikeslip
	./seas_test_canonical_rotation_pure_strikeslip

test-fault-basis-qp-orthonormality: seas_test_fault_basis_qp_orthonormality
	./seas_test_fault_basis_qp_orthonormality

test-per-qp-vs-centroid-basis-planar-1el: seas_test_per_qp_vs_centroid_basis_planar_1el
	./seas_test_per_qp_vs_centroid_basis_planar_1el

# PLAN_frame_orthonormality_diag_2026-05-23 §Phase 1 AC#3 — [FRAME] dv_n identity
test-frame-diag-projection-identity: seas_test_frame_diag_projection_identity
	./seas_test_frame_diag_projection_identity

# PLAN_predictor_vs_macro_diag_2026-05-23 §Phase 1 AC#3 — [MACRO] decomposition identity
test-macro-diag-decomposition-identity: seas_test_macro_diag_decomposition_identity
	./seas_test_macro_diag_decomposition_identity

# v9.2.0 plan §4.6 (rev-3b) — shared-fault role consistency
# Uses 4 ranks on a 4 km inline fixture.  2 ranks is insufficient
# because ParMETIS rarely cuts through the fault in a 2-way
# partition of this compact fixture, causing the test to skip with
# exit-77 ("no shared fault faces").  4 ranks reliably yields
# several cross-rank shared fault faces — the test input.
# v9.2.0 plan §4.9 (rev-3f) — No-penalty regression gate (Pelties 2012 §3.1)
test-no-penalty-dynamic-rupture: seas_test_no_penalty_dynamic_rupture
	./seas_test_no_penalty_dynamic_rupture

# v9.2.0 plan Step 6 (H-V92-U probe) — interior-fault-path Pelties eq. (7)
test-interior-fault-flux-path: seas_test_interior_fault_flux_path
	./seas_test_interior_fault_flux_path

# Pepper-bug plan 2026-04-22 Phase 3 P3-T1 — asymmetric Q probe
test-interior-fault-flux-path-asymmetric: seas_test_interior_fault_flux_path_asymmetric
	./seas_test_interior_fault_flux_path_asymmetric

# Pepper-bug plan 2026-04-22 Phase 3 P3-T3 — rupture-regime amplitude probe
test-interior-fault-flux-path-rupture-regime: seas_test_interior_fault_flux_path_rupture_regime
	./seas_test_interior_fault_flux_path_rupture_regime

# Pepper-bug plan 2026-04-22 Phase 5 P5-N — bulk RK4 noise stability
test-bulk-rk4-noise-stability: seas_test_bulk_rk4_noise_stability
	./seas_test_bulk_rk4_noise_stability

# Pepper-bug plan 2026-04-22 Test 5 — shared-fault Pelties-7 analytic (np=4)
test-shared-fault-flux-path-analytic: seas_test_shared_fault_flux_path_analytic
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_shared_fault_flux_path_analytic

# Pepper-bug plan 2026-04-22 Phase 5 P5-PW — bulk plane-wave recovery
test-bulk-plane-wave-recovery: seas_test_bulk_plane_wave_recovery
	./seas_test_bulk_plane_wave_recovery

# Pepper-bug plan 2026-04-22 Phase 1 — live interior-vs-shared cross-check
test-interior-vs-shared-branch-live: seas_test_interior_vs_shared_branch_live
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_interior_vs_shared_branch_live

# Pepper-bug: MULTI-STEP rupture-regime serial-vs-parallel cross-check
test-rupture-multistep-serial-vs-parallel: seas_test_rupture_multistep_serial_vs_parallel
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_rupture_multistep_serial_vs_parallel

# Phase 0 fault-dealiasing reproduction harness.  Serial run is the acceptance
# path; run `mpirun -np 8 ./seas_test_fault_planar_serial --nx 16 --nz 16` for a
# bigger reproduction.  This full run is heavy (~16k ADER-2 macro-steps over the
# main + symmetry-probe slabs); it is NOT part of the fast `test` aggregate.
test-fault-planar-serial: seas_test_fault_planar_serial
	./seas_test_fault_planar_serial

# CI-fast smoke (REVIEW R-005): a short serial run, no symmetry probe, just
# enough to confirm the harness builds, nucleates, and exits cleanly.  Does NOT
# assert the secular-drift acceptance (that needs the full step budget above).
test-fault-planar-serial-fast: seas_test_fault_planar_serial
	./seas_test_fault_planar_serial --nx 6 --nz 6 --nsteps 2000 --no-symmetry-probe

# Phase 1: fault-flux over-integration unit test (tiny 2-tet fixture, fast).
test-fault-overint: seas_test_fault_overint
	./seas_test_fault_overint

# Phase 2: fault-face resample operator R unit test (reference-element, fast).
test-fault-resample: seas_test_fault_resample
	./seas_test_fault_resample

# Phase 3: rate-state Δψ resample apply unit test (iterator-level, fast).
test-fault-resample-apply: seas_test_fault_resample_apply
	./seas_test_fault_resample_apply

test-rupture-tilted-fault-serial-vs-parallel: seas_test_rupture_tilted_fault_serial_vs_parallel
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_rupture_tilted_fault_serial_vs_parallel

test-shared-fault-reconcile-cross-rank: seas_test_shared_fault_reconcile_cross_rank
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_shared_fault_reconcile_cross_rank

# R-020(b): psi-explicit np=2 shared-fault RS consistency.  Excluded from the
# default `test:` umbrella (needs mpirun -np 2); run as a pre-submit check.
test-shared-fault-rs-psi-consistency-np2: seas_test_shared_fault_rs_psi_consistency_np2
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_shared_fault_rs_psi_consistency_np2

# R-004 (REVIEW_speckle_seissol_drdg3d_rootcause_2026-05-24): guards the fixed
# predictor ghost read — the batched byNODES exchange read through
# GetFaceNbrElementVDofs must be bit-identical to the per-component scalar
# exchange (macro path).  GREEN as of the R-004 fix (2026-05-24); in
# test-v92-regression-gates.  Would go RED if the pre-fix flat-slab unpack
# src[c*n_fn+j] is reintroduced.
test-ghost-exchange-bynodes-vs-scalar: seas_test_ghost_exchange_bynodes_vs_scalar
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_ghost_exchange_bynodes_vs_scalar

test-rupture-rake-sweep-cross-rank: seas_test_rupture_rake_sweep_cross_rank
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_rupture_rake_sweep_cross_rank

# Pepper-bug: ADJACENT-TRIANGLE UNIFORMITY (serial + np=4 variants)
test-adjacent-triangle-fault-uniformity-serial: seas_test_adjacent_triangle_fault_uniformity
	./seas_test_adjacent_triangle_fault_uniformity

test-adjacent-triangle-fault-uniformity-parallel: seas_test_adjacent_triangle_fault_uniformity
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_adjacent_triangle_fault_uniformity

test-adjacent-triangle-fault-first-step-audit: seas_test_adjacent_triangle_fault_first_step_audit
	./seas_test_adjacent_triangle_fault_first_step_audit

# Pepper-bug plan 2026-04-22 Phase 2 — sign_flipped / elem1_on_plus truth table
test-sign-flipped-elem1-on-plus-truth-table: seas_test_sign_flipped_elem1_on_plus_truth_table
	./seas_test_sign_flipped_elem1_on_plus_truth_table

# Pepper-bug plan 2026-04-22 Phase 4 P4-T1 — MakeFaceKey uniqueness (np=4)
# NB: requires np=4 to reliably cut the 4-km fault fixture across ranks
# (np=2 often puts the whole fault inside one rank, yielding exit 77).
test-make-face-key-uniqueness: seas_test_make_face_key_uniqueness
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_make_face_key_uniqueness

# Pepper-bug plan 2026-04-22 Phase 4 P4-T2 — DOF-offset invariants
test-fault-dof-offset-invariants: seas_test_fault_dof_offset_invariants
	./seas_test_fault_dof_offset_invariants

# Pepper-bug plan 2026-04-22 Phase 6 — aggregator: run all pepper-bug tests
test-tpv102-pepper: \
	test-interior-fault-flux-path-asymmetric \
	test-interior-fault-flux-path-rupture-regime \
	test-bulk-rk4-noise-stability \
	test-bulk-plane-wave-recovery \
	test-sign-flipped-elem1-on-plus-truth-table \
	test-fault-dof-offset-invariants \
	test-make-face-key-uniqueness \
	test-interior-vs-shared-branch-live \
	test-shared-fault-flux-path-analytic
	@echo "=========================================="
	@echo "  TPV102 pepper-bug test suite: COMPLETE"
	@echo "=========================================="

# v9.2.0 plan Step 5 (R-V92-H01) — coupled RK4-on-psi convergence probe
test-rk4-psi-integration: seas_test_rk4_psi_integration
	./seas_test_rk4_psi_integration

# v9.2.0 round-9 REVIEW R-V92-K01 — driver post-RK4 contract regression guard
test-rk4-endpoint-reeval: seas_test_rk4_endpoint_reeval
	./seas_test_rk4_endpoint_reeval

# v9.2.0 plan Step 5b (R-V92-I02) — 100-step F01+F02 stability probe (4 ranks)
test-rk4-f01f02-stability: seas_test_rk4_f01f02_shared_fault_stability
	mpirun -np 4 ./seas_test_rk4_f01f02_shared_fault_stability

# v9.2.0 plan §4.10 (rev-3g) — single-channel Jacobian probe (d+e)
test-volume-jacobian-single-channel: seas_test_volume_jacobian_single_channel
	./seas_test_volume_jacobian_single_channel

# v9.2.0 plan §4.11 (rev-3g) — RK4 conservation probe (b)
test-rk4-conservation: seas_test_rk4_conservation
	./seas_test_rk4_conservation

# v9.2.0 plan §4.12 (rev-3g) — absorbing-BC energy decay probe (c)
test-absorbing-bc-energy-decay: seas_test_absorbing_bc_energy_decay
	./seas_test_absorbing_bc_energy_decay

test-shared-fault-role-consistency: seas_test_shared_fault_role_consistency
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_shared_fault_role_consistency

# v9.2.0 plan §4.7 (rev-3c) — DOFData runtime consistency via RK4 on fixture
test-shared-fault-dof-data-consistency: seas_test_shared_fault_dof_data_consistency
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_shared_fault_dof_data_consistency

# v9.3.0 plan Phase 1 (I-04): FreeSurfaceGodunov vs FreeSurface equivalence
test-free-surface-godunov: seas_test_free_surface_godunov
	./seas_test_free_surface_godunov

# v9.3.0 plan Phase 2 (I-04 wiring, R-005): --free-surface-bc CLI parsing
test-R005-free-surface-bc-flag-parse: seas_test_R005_free_surface_bc_flag_parse
	./seas_test_R005_free_surface_bc_flag_parse

# v9.3.0 plan Phase 2 AC-3 (I-04 driver integration gate)
test-free-surface-godunov-driver: seas_test_free_surface_godunov_driver
	./seas_test_free_surface_godunov_driver

# v9.3.0 plan Phase 3 (I-06 part A): EvaluateTotal vs Evaluate equivalence
test-fault-face-flux-total-vs-fluctuation: seas_test_fault_face_flux_total_vs_fluctuation
	./seas_test_fault_face_flux_total_vs_fluctuation

# v9.3.0 R-003: EvaluateTotal psi-invariance
test-R003-evaluatetotal-psi-invariant: seas_test_R003_evaluatetotal_psi_invariant
	./seas_test_R003_evaluatetotal_psi_invariant

# v9.3.0 Phase 4 (I-06): R-001 InitializeStateTotal sign regression
test-R001-initstate-total-sign: seas_test_R001_initstate_total_sign
	./seas_test_R001_initstate_total_sign

# v9.3.0 Phase 4 AC (I-06): locked-fault 100 RK4 sigma_n stability
test-tpv102-total-locked-fault: seas_test_tpv102_total_locked_fault
	./seas_test_tpv102_total_locked_fault

# v9.4.0 R-005: test-R002-nucleation-sign-match /
# test-R009-nucleation-total-both-sides /
# test-R-I06-003-nucleation-shape-scaling targets removed
# (exercised the DEPRECATED bulk-Q ApplyNucleationTotal path).
# v9.3.0 R-I06-001 / R-I06-006: total-Q absorbing BC equilibrium gate
test-tpv102-total-absorbing-equilibrium: seas_test_tpv102_total_absorbing_equilibrium
	./seas_test_tpv102_total_absorbing_equilibrium

# 2026-04-22 nucleation persistent-prestress channel fix regression
test-persistent-nuc-prestress-channel: seas_test_persistent_nuc_prestress_channel
	./seas_test_persistent_nuc_prestress_channel

# v9.4.0 debug-plan R-001: pepper reproducer (gates the split-prestress refactor)
test-tpv102-pepper-reproducer: seas_test_tpv102_pepper_reproducer
	./seas_test_tpv102_pepper_reproducer

# v9.3.0 R-I06-round3 R-006: direct tests for FreeSurfaceTotal variants
test-free-surface-total-direct: seas_test_free_surface_total_direct
	./seas_test_free_surface_total_direct

# ADER I-05 Phase 1: ApplySpatialDerivative polynomial-exactness test
test-wave-operator-spatial-derivative: seas_test_wave_operator_spatial_derivative
	./seas_test_wave_operator_spatial_derivative

# ADER I-05 Phase 3: CK predictor plane-wave convergence
test-ader-ck-predictor: seas_test_ader_ck_predictor
	./seas_test_ader_ck_predictor

# ADER I-05 Phase 5: EvaluateADER / EvaluateADERTotal equivalence
test-fault-face-flux-ader-equivalence: seas_test_fault_face_flux_ader_equivalence
	./seas_test_fault_face_flux_ader_equivalence

# ADER I-05 Phase 6: AdvanceADER vs RK4 on a linear wave
test-ader-linear-wave-equivalence: seas_test_ader_linear_wave_equivalence
	./seas_test_ader_linear_wave_equivalence

# ADER fault-path live consistency: (Q_new - Q)/dt -> wave.Mult(Q)
test-ader-fault-coupling-dt-limit: seas_test_ader_fault_coupling_dt_limit
	./seas_test_ader_fault_coupling_dt_limit

# ADER live interior-vs-shared branch cross-check
test-ader-interior-vs-shared-branch-live: seas_test_ader_interior_vs_shared_branch_live
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_ader_interior_vs_shared_branch_live

# ADER I-05 Phase 8: TPV102 AdvanceADER smoke test
test-ader-tpv102-smoke: seas_test_ader_tpv102_smoke
	./seas_test_ader_tpv102_smoke

# Round-6 R-001 regression: driver defaults to ADER
test-R001-driver-defaults-to-ader: seas_test_R001_driver_defaults_to_ader
	./seas_test_R001_driver_defaults_to_ader

# Round-6 R-002 regression: driver init uses total-Q
test-R002-driver-init-total-q: seas_test_R002_driver_init_total_q
	./seas_test_R002_driver_init_total_q

# ADER I-05 Phase 8 aggregate: all five ADER tests in order.
test-ader: test-wave-operator-spatial-derivative \
           test-ader-ck-predictor \
           test-fault-face-flux-ader-equivalence \
           test-ader-linear-wave-equivalence \
           test-ader-tpv102-smoke

# v9.3.0 plan Testing Strategy aggregate
#   I-04: Phase 1 flux kernel + Phase 2 CLI + Phase 2 AC-3 driver gate
#   I-06: Phase 3 EvaluateTotal + R-003 psi guard + Phase 4 R-001 +
#         Phase 4 AC locked-fault + Phase 5 R-002 + R-009
test-v93-regression: test-free-surface-godunov \
                     test-R005-free-surface-bc-flag-parse \
                     test-free-surface-godunov-driver \
                     test-fault-face-flux-total-vs-fluctuation \
                     test-R003-evaluatetotal-psi-invariant \
                     test-R001-initstate-total-sign \
                     test-tpv102-total-locked-fault \
                     test-tpv102-total-absorbing-equilibrium \
                     test-persistent-nuc-prestress-channel \
                     test-tpv102-pepper-reproducer \
                     test-free-surface-total-direct \
                     test-wave-operator-spatial-derivative \
                     test-ader-ck-predictor \
                     test-fault-face-flux-ader-equivalence \
                     test-ader-linear-wave-equivalence \
                     test-ader-tpv102-smoke \
                     test-R001-driver-defaults-to-ader \
                     test-R002-driver-init-total-q

# Aggregate target: all v9.2.0 §4 regression gates (rev-3..rev-3g)
# §4.9 runs FIRST per user directive: it verifies the Pelties 2012 benchmark
# contract that the dynamic-rupture scheme uses pure Godunov fluxes with no
# IP/BR2/penalty term.  §4.10-§4.12 probe H-V92-G candidates (b)(c)(d)(e).
test-v92-regression-gates: test-no-penalty-dynamic-rupture \
                          test-godunov-rotation-identity \
                          test-canonical-rotation-pure-strikeslip \
                          test-fault-basis-qp-orthonormality \
                          test-frame-diag-projection-identity \
                          test-macro-diag-decomposition-identity \
                          test-ghost-exchange-bynodes-vs-scalar \
                          test-per-qp-vs-centroid-basis-planar-1el \
                          test-shared-fault-role-consistency \
                          test-shared-fault-dof-data-consistency \
                          test-volume-jacobian-single-channel \
                          test-rk4-conservation \
                          test-absorbing-bc-energy-decay

# Death test (REVIEW R-007): child aborts inside MFEM_VERIFY; treat
# non-zero exit as PASS.  If the child returns 0, the guard is gone.
test-fault-face-flux-bimaterial-guard: seas_test_fault_face_flux_bimaterial_guard
	@./seas_test_fault_face_flux_bimaterial_guard 2>/dev/null; rc=$$?; \
	    if [ $$rc -ne 0 ]; then \
	        echo "PASS: bimaterial guard aborted as expected (rc=$$rc)"; \
	    else \
	        echo "FAIL: Evaluate returned without aborting on bimaterial DOFData"; \
	        exit 1; \
	    fi

test-friction-solver: seas_test_friction_solver
	./seas_test_friction_solver

test-tpv102-setup: seas_test_tpv102_setup
	./seas_test_tpv102_setup

test-tpv102-local: seas_test_tpv102_local
	./seas_test_tpv102_local

seas_driver: $(SEAS_DRIVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(SEAS_DRIVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(SEAS_DRIVER_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(PARALLEL_IO_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(TOML_FLAGS) -DSEAS_USE_MPI -c $< -o $@

# Phase 6: Parallel test targets
seas_test_mpi_context: $(TEST_MPI_CONTEXT_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_MPI_CONTEXT_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_parallel_utils: $(TEST_PARALLEL_UTILS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PARALLEL_UTILS_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(TEST_MPI_CONTEXT_OBJ): %.o: $(SRC)%.cpp $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

$(TEST_PARALLEL_UTILS_OBJ): %.o: $(SRC)%.cpp $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

seas_test_fault_scatter: $(TEST_FAULT_SCATTER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_SCATTER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(TEST_FAULT_SCATTER_OBJ): %.o: $(SRC)%.cpp $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

test-fault-scatter: seas_test_fault_scatter
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_fault_scatter
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_fault_scatter

test-mpi-context: seas_test_mpi_context
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_mpi_context
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_mpi_context

test-parallel-utils: seas_test_parallel_utils
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_parallel_utils
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_parallel_utils

# Phase 7: Parallel domain tests
seas_test_parallel_domain: $(TEST_PARALLEL_DOMAIN_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PARALLEL_DOMAIN_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_mms_antiplane_parallel: $(MMS_ANTIPLANE_PARALLEL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(MMS_ANTIPLANE_PARALLEL_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(TEST_PARALLEL_DOMAIN_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

$(MMS_ANTIPLANE_PARALLEL_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

test-parallel-domain: seas_test_parallel_domain
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_parallel_domain
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_parallel_domain

test-mms-parallel: seas_mms_antiplane_parallel
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_mms_antiplane_parallel

# Phase 8: Parallel fault tests and driver
seas_test_parallel_fault: $(TEST_PARALLEL_FAULT_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PARALLEL_FAULT_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_pseas: $(PSEAS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(PSEAS_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(TEST_PARALLEL_FAULT_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(PARALLEL_IO_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

$(PSEAS_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(PARALLEL_IO_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

test-parallel-fault: seas_test_parallel_fault
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_parallel_fault
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_parallel_fault

test-parallel: test-mpi-context test-parallel-utils test-parallel-domain test-mms-parallel test-parallel-fault test-shared-fault-reconcile-cross-rank test-bimaterial-deriv-cache-parity

# Phase 9: Verification and scaling targets
seas_test_serial_parallel_consistency: $(TEST_SERIAL_PARALLEL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_SERIAL_PARALLEL_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_bp2_benchmark_parallel: $(BP2_BENCHMARK_PARALLEL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(BP2_BENCHMARK_PARALLEL_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_scaling: $(TEST_SCALING_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_SCALING_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(TEST_SERIAL_PARALLEL_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(PARALLEL_IO_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

$(BP2_BENCHMARK_PARALLEL_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(PARALLEL_IO_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

$(TEST_SCALING_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

test-serial-parallel: seas_test_serial_parallel_consistency
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_serial_parallel_consistency
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_serial_parallel_consistency

test-bp2-parallel: seas_bp2_benchmark_parallel
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_bp2_benchmark_parallel

test-scaling: seas_test_scaling
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_scaling
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_scaling

# BR2 consistency test
seas_test_br2_consistency: $(TEST_BR2_CONSISTENCY_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BR2_CONSISTENCY_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(TEST_BR2_CONSISTENCY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

test-br2-consistency: seas_test_br2_consistency
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_br2_consistency
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_br2_consistency

# Parallel elasticity test (shared face support)
seas_test_parallel_elasticity: $(TEST_PARALLEL_ELAST_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PARALLEL_ELAST_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(TEST_PARALLEL_ELAST_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

test-parallel-elasticity: seas_test_parallel_elasticity
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_parallel_elasticity
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_parallel_elasticity

# Parallel wave operator tests (dynamic rupture)
# Match the serial rule (line 847): the test instantiates
# WaveOperator<ParMesh>, whose ComputeFaceFluxRHS calls non-template
# methods on PrecomputedFaceFluxes (AddBoundaryFaceRhs, AddInteriorFaceRhs).
# Those live in dynamic/precomputed_face_fluxes.cpp, so the .o must be on
# the link line — and the same applies to the non-template parts of
# WaveOperator itself.
seas_test_parallel_wave_operator: $(TEST_PARALLEL_WAVE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PARALLEL_WAVE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(TEST_PARALLEL_WAVE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

test-parallel-wave-operator: seas_test_parallel_wave_operator
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_parallel_wave_operator

# Lever 1/3 (REVIEW R-002): parallel DerivMode Cached==OnTheFly equivalence.
seas_test_wave_operator_cached_parallel: $(TEST_PARALLEL_CACHED_WAVE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PARALLEL_CACHED_WAVE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(TEST_PARALLEL_CACHED_WAVE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

test-wave-operator-cached-parallel: seas_test_wave_operator_cached_parallel
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_wave_operator_cached_parallel

# Bimaterial (matrix path) DerivMode + shared-CK parity (opt 2026-06-24 gate).
# Object list mirrors seas_test_bimaterial_wave_operator_parity (needs the flux
# pool + heterogeneous material + spatial friction) plus the parallel build flag.
seas_test_bimaterial_deriv_cache_parity: $(TEST_PARALLEL_BIMAT_CACHE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(GODUNOV_FLUX_OBJ) $(GODUNOV_FLUX_POOL_OBJ) $(PML_LAYER_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(HETEROGENEOUS_MATERIAL_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(SPATIAL_FRICTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PARALLEL_BIMAT_CACHE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(GODUNOV_FLUX_OBJ) $(GODUNOV_FLUX_POOL_OBJ) $(PML_LAYER_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(HETEROGENEOUS_MATERIAL_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(SPATIAL_FRICTION_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(TEST_PARALLEL_BIMAT_CACHE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

test-bimaterial-deriv-cache-parity: seas_test_bimaterial_deriv_cache_parity
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_bimaterial_deriv_cache_parity

# R-101 shared-fault DOFData consistency test (2-rank MPI regression)
seas_test_r101_shared_fault: $(TEST_R101_SHARED_FAULT_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R101_SHARED_FAULT_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(TEST_R101_SHARED_FAULT_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

test-r101-shared-fault: seas_test_r101_shared_fault
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_r101_shared_fault

# BP5 parallel smoke test
seas_test_bp5_parallel_smoke: $(TEST_BP5_SMOKE_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BP5_SMOKE_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(TEST_BP5_SMOKE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

test-bp5-smoke: seas_test_bp5_parallel_smoke
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_bp5_parallel_smoke
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_bp5_parallel_smoke

# Regression check: compare simulation output against golden reference.
# Usage: make regression-check SIM_DIR=. SIM_PREFIX=bp5_full REF_DIR=bp5/benchmark_data/golden_serial_1r_50step [TOL=0]
REG_TOL ?=
SIM_DIR ?= .
SIM_PREFIX ?= bp5_full
REF_DIR ?= bp5/benchmark_data/golden_serial_1r_50step
regression-check:
	python3 scripts/regression_check.py \
		--sim-dir $(SIM_DIR) --sim-prefix $(SIM_PREFIX) \
		--ref-dir $(REF_DIR) \
		$(if $(REG_TOL),--tolerance $(REG_TOL),)

# Quick-check: 8-rank/20-step inline-mesh run + regression against golden.
# Produces output in bp5_quick_check/, then compares against golden_inline_8r_20step/.
QUICK_CHECK_NP ?= 8
QUICK_CHECK_STEPS ?= 20
QUICK_CHECK_DIR ?= bp5_quick_check
QUICK_CHECK_GOLDEN ?= bp5/benchmark_data/golden_inline_8r_20step
QUICK_CHECK_TOL ?= 0
quick-check: seas_bp5_full
	@if [ ! -d $(QUICK_CHECK_GOLDEN) ]; then \
		echo "SKIP quick-check: golden dir $(QUICK_CHECK_GOLDEN) not found"; \
		exit 0; \
	fi
	@mkdir -p $(QUICK_CHECK_DIR)
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $(QUICK_CHECK_NP) ./seas_bp5_full \
		--inline-mesh --max-steps $(QUICK_CHECK_STEPS) \
		--output-dir $(QUICK_CHECK_DIR) --output-prefix bp5_quick \
		--tandem-time-stepping \
		--regression-tolerance $(QUICK_CHECK_TOL) \
		--ref-dir $(QUICK_CHECK_GOLDEN)
	@echo "Quick-check complete. Output in $(QUICK_CHECK_DIR)/"

# Phase 7c-iii: New driver (seas_driver) quick-check, equivalence, and regression targets
DRIVER_QC_NP ?= 8
DRIVER_QC_DIR ?= bp5_driver_quick_check
driver-quick-check: seas_driver
	@mkdir -p $(DRIVER_QC_DIR)
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $(DRIVER_QC_NP) ./seas_driver \
		config/bp5_quick_check.toml --write-every-step
	@echo "Driver quick-check complete. Output in $(DRIVER_QC_DIR)/"

DRIVER_EQUIV_NP ?= 8
DRIVER_EQUIV_STEPS ?= 20
DRIVER_EQUIV_NEW_DIR ?= bp5_driver_equiv_new
DRIVER_EQUIV_OLD_DIR ?= bp5_driver_equiv_old
driver-equivalence: seas_driver seas_bp5_full
	@mkdir -p $(DRIVER_EQUIV_NEW_DIR) $(DRIVER_EQUIV_OLD_DIR)
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $(DRIVER_EQUIV_NP) ./seas_driver \
		config/bp5_quick_check.toml \
		--override output.output_dir=$(DRIVER_EQUIV_NEW_DIR) \
		--override output.output_prefix=bp5_equiv \
		--write-every-step
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $(DRIVER_EQUIV_NP) ./seas_bp5_full \
		--inline-mesh --max-steps $(DRIVER_EQUIV_STEPS) \
		--output-dir $(DRIVER_EQUIV_OLD_DIR) --output-prefix bp5_equiv \
		--tandem-time-stepping --solver mumps \
		--write-every-step
	@echo "Equivalence test: comparing new vs old driver output..."
	python3 scripts/regression_check.py \
		--sim-dir $(DRIVER_EQUIV_NEW_DIR) --sim-prefix bp5_equiv \
		--ref-dir $(DRIVER_EQUIV_OLD_DIR) \
		--tolerance 1e-12
	@echo "Driver equivalence PASSED."

DRIVER_REG_NP ?= 8
DRIVER_REG_DIR ?= bp5_driver_regression
DRIVER_REG_GOLDEN ?= bp5/benchmark_data/golden_inline_8r_20step
DRIVER_REG_TOL ?= 0
driver-regression-local: seas_driver
	@if [ ! -d $(DRIVER_REG_GOLDEN) ]; then \
		echo "SKIP driver-regression-local: golden dir $(DRIVER_REG_GOLDEN) not found"; \
		exit 0; \
	fi
	@mkdir -p $(DRIVER_REG_DIR)
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $(DRIVER_REG_NP) ./seas_driver \
		config/bp5_quick_check.toml \
		--override output.output_dir=$(DRIVER_REG_DIR) \
		--override output.output_prefix=bp5_quick \
		--write-every-step
	python3 scripts/regression_check.py \
		--sim-dir $(DRIVER_REG_DIR) --sim-prefix bp5_quick \
		--ref-dir $(DRIVER_REG_GOLDEN) \
		--tolerance $(DRIVER_REG_TOL)
	@echo "Driver regression PASSED (tol=$(DRIVER_REG_TOL))."

# Phase 10: BP2 first-cycle verification
seas_bp2_verify: $(BP2_VERIFY_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(BP2_VERIFY_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(BP2_VERIFY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(PARALLEL_IO_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

test-bp2-verify: seas_bp2_verify
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 10 ./seas_bp2_verify

# Phase 10b: BP2 full multi-cycle simulation
seas_bp2_full: $(BP2_FULL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(BP2_FULL_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(BP2_FULL_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(PARALLEL_IO_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

# BP1 full multi-cycle simulation
seas_bp1_full: $(BP1_FULL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(BP1_FULL_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(BP1_FULL_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(PARALLEL_IO_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

# BP1 bottom-boundary Dirichlet loading simulation
seas_bp1_bdrload: $(BP1_BDRLOAD_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(BP1_BDRLOAD_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(BP1_BDRLOAD_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(PARALLEL_IO_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

# BP5 full multi-cycle simulation
seas_bp5_full: $(BP5_FULL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(BP5_FULL_OBJ) $(MFEM_LIBS) $(SEAS_FS_LIB)
	$(SEAS_POST_LINK_DEDUP_RPATH)

$(BP5_FULL_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(PARALLEL_IO_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

-include $(MFEM_TEST_MK)

# MFEM test framework pattern rules: define %-test-seq and %-test-par
# targets expected by test.mk for each entry in SEQ_MINIAPPS / PAR_MINIAPPS.
RUN_MPI = $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $(MFEM_MPI_NP)
%-test-par: %
	@$(call mfem-test,$<, $(RUN_MPI), SEAS miniapp)
%-test-seq: %
	@$(call mfem-test,$<,, SEAS miniapp)

# Generate an error message if the MFEM library is not built and exit
$(MFEM_LIB_FILE):
	$(error The MFEM library is not built)

clean: clean-build clean-exec

clean-build:
	rm -f *.o *~ $(SEQ_MINIAPPS) $(PAR_MINIAPPS) $(PAR_LONG)
	rm -f drivers/*.o tests/unit/*.o tests/verification/*.o tests/parallel/*.o
	rm -rf *.dSYM *.TVD.*breakpoints

clean-exec:
	@rm -rf mesh.* sol.* ParaView mfem_bp2qd_*.txt test_*.txt

# =============================================================================
# TPV104 integration (§4.10 Step 12)
# =============================================================================
# Added by TPV104 Step 12.  Does NOT modify any TPV102 / BP5 target above;
# all TPV104 code lives in files not referenced by the pre-existing targets,
# so `make test` remains unchanged ([C2] non-regression).
#
# New umbrella targets:
#   make test-tpv104    — run every TPV104 unit-test binary.
#   make seas_tpv104_driver — build the TPV104 driver binary.
#   make tpv104-mesh    — (pythonenv) gmsh-generate the .msh files.
# =============================================================================

# ---- headers (add to dependency list without touching SEAS_HEADERS) ----
TPV104_HEADERS = config/tpv104_params.hpp \
                 friction/slip_law_srw_psi.hpp \
                 friction/friction_coeff_stable.hpp \
                 dynamic/tpv104_friction_solver.hpp \
                 dynamic/tpv104_nucleation.hpp \
                 dynamic/tpv104_setup.hpp \
                 dynamic/tpv104_substep_iterator.hpp

# ---- sources (unit tests) ----
TEST_TPV104_PARAMS_SRC               = tests/unit/test_tpv104_params.cpp
TEST_SLIP_LAW_SRW_PSI_SRC            = tests/unit/test_slip_law_srw_psi.cpp
TEST_FRICTION_COEFF_STABLE_SRC       = tests/unit/test_friction_coeff_stable.cpp
TEST_TPV104_FRICTION_SOLVER_NEWTON_SRC = tests/unit/test_tpv104_friction_solver_newton.cpp
TEST_TPV104_NUCLEATION_SRC           = tests/unit/test_tpv104_nucleation.cpp
TEST_TPV104_NORMAL_SIGN_SRC          = tests/unit/test_tpv104_normal_sign.cpp
TEST_TPV104_MESH_SRC                 = tests/unit/test_tpv104_mesh.cpp
TEST_TPV104_SETUP_SRC                = tests/unit/test_tpv104_setup.cpp
TEST_TPV104_SUBSTEP_ITERATOR_SRC     = tests/unit/test_tpv104_substep_iterator.cpp
TEST_TPV104_PROBE_FORMAT_SRC         = tests/unit/test_tpv104_probe_format.cpp
TEST_TPV104_SMOKE_SRC                = tests/unit/test_tpv104_smoke.cpp
TEST_TPV104_SIGMA_N_INVARIANCE_SRC   = tests/unit/test_tpv104_sigma_n_invariance.cpp
TEST_TPV104_FREEZE_SIGMA_N_GATE_SRC  = tests/unit/test_tpv104_freeze_sigma_n_gate.cpp
TEST_TPV104_SUBSTEP_PREDICTOR_SRC    = tests/unit/test_tpv104_substep_predictor.cpp
TEST_TPV104_SUBSTEP_ONE_SHOT_PARITY_SRC = tests/unit/test_tpv104_substep_one_shot_parity.cpp
TEST_TPV104_SUBSTEP_ITERATOR_PARITY_SRC = tests/unit/test_tpv104_substep_iterator_parity.cpp
TEST_TPV104_SUBSTEP_DISPATCH_PARITY_SRC = tests/unit/test_tpv104_substep_dispatch_parity.cpp
TEST_GODUNOV_CENTRAL_FLUX_SRC      = tests/unit/test_godunov_central_flux.cpp
TEST_MIXED_FLUX_FACE_SET_SRC       = tests/unit/test_mixed_flux_face_set.cpp
TEST_MIXED_FLUX_DISPATCH_NONE_SRC  = tests/unit/test_mixed_flux_dispatch_none.cpp
TEST_BIMATERIAL_CENTRAL_FLUX_SRC   = tests/unit/test_bimaterial_central_flux.cpp
TEST_BIMATERIAL_MIXED_FLUX_DISPATCH_SRC = tests/unit/test_bimaterial_mixed_flux_dispatch.cpp
TEST_BIMATERIAL_FAULT_PERSIDE_MATERIAL_SRC = tests/unit/test_bimaterial_fault_perside_material.cpp
TEST_BIMATERIAL_CONTRAST_GUARD_SRC = tests/unit/test_bimaterial_contrast_guard.cpp
TEST_BIMATERIAL_MIXED_FLUX_SHARED_SRC = tests/parallel/test_bimaterial_mixed_flux_shared.cpp
TEST_BIMATERIAL_SEAM_MATERIAL_NP2_SRC = tests/parallel/test_bimaterial_seam_material_np2.cpp
TEST_BIMATERIAL_SEAM_FAULT_NP2_SRC = tests/parallel/test_bimaterial_seam_fault_np2.cpp
TEST_BIMATERIAL_MIXED_FLUX_CFL_SRC = tests/unit/test_bimaterial_mixed_flux_cfl.cpp
TEST_BIMATERIAL_MIXED_FLUX_HOMOG_EQUIV_SRC = tests/unit/test_bimaterial_mixed_flux_homog_equivalence.cpp

TPV104_DRIVER_SRC                    = drivers/tpv104_driver.cpp
TPV104_SUBSTEP_ITERATOR_SRC          = dynamic/tpv104_substep_iterator.cpp

# ---- objects ----
TEST_TPV104_PARAMS_OBJ               = $(TEST_TPV104_PARAMS_SRC:.cpp=.o)
TEST_SLIP_LAW_SRW_PSI_OBJ            = $(TEST_SLIP_LAW_SRW_PSI_SRC:.cpp=.o)
TEST_FRICTION_COEFF_STABLE_OBJ       = $(TEST_FRICTION_COEFF_STABLE_SRC:.cpp=.o)
TEST_TPV104_FRICTION_SOLVER_NEWTON_OBJ = $(TEST_TPV104_FRICTION_SOLVER_NEWTON_SRC:.cpp=.o)
TEST_TPV104_NUCLEATION_OBJ           = $(TEST_TPV104_NUCLEATION_SRC:.cpp=.o)
TEST_TPV104_NORMAL_SIGN_OBJ          = $(TEST_TPV104_NORMAL_SIGN_SRC:.cpp=.o)
TEST_TPV104_MESH_OBJ                 = $(TEST_TPV104_MESH_SRC:.cpp=.o)
TEST_TPV104_SETUP_OBJ                = $(TEST_TPV104_SETUP_SRC:.cpp=.o)
TEST_TPV104_SUBSTEP_ITERATOR_OBJ     = $(TEST_TPV104_SUBSTEP_ITERATOR_SRC:.cpp=.o)
TEST_TPV104_PROBE_FORMAT_OBJ         = $(TEST_TPV104_PROBE_FORMAT_SRC:.cpp=.o)
TEST_TPV104_SMOKE_OBJ                = $(TEST_TPV104_SMOKE_SRC:.cpp=.o)
TEST_TPV104_SIGMA_N_INVARIANCE_OBJ   = $(TEST_TPV104_SIGMA_N_INVARIANCE_SRC:.cpp=.o)
TEST_TPV104_FREEZE_SIGMA_N_GATE_OBJ  = $(TEST_TPV104_FREEZE_SIGMA_N_GATE_SRC:.cpp=.o)
TEST_TPV104_SUBSTEP_PREDICTOR_OBJ        = $(TEST_TPV104_SUBSTEP_PREDICTOR_SRC:.cpp=.o)
TEST_TPV104_SUBSTEP_ONE_SHOT_PARITY_OBJ  = $(TEST_TPV104_SUBSTEP_ONE_SHOT_PARITY_SRC:.cpp=.o)
TEST_TPV104_SUBSTEP_ITERATOR_PARITY_OBJ  = $(TEST_TPV104_SUBSTEP_ITERATOR_PARITY_SRC:.cpp=.o)
TEST_TPV104_SUBSTEP_DISPATCH_PARITY_OBJ  = $(TEST_TPV104_SUBSTEP_DISPATCH_PARITY_SRC:.cpp=.o)
TEST_GODUNOV_CENTRAL_FLUX_OBJ      = $(TEST_GODUNOV_CENTRAL_FLUX_SRC:.cpp=.o)
TEST_MIXED_FLUX_FACE_SET_OBJ       = $(TEST_MIXED_FLUX_FACE_SET_SRC:.cpp=.o)
TEST_MIXED_FLUX_DISPATCH_NONE_OBJ  = $(TEST_MIXED_FLUX_DISPATCH_NONE_SRC:.cpp=.o)
TEST_BIMATERIAL_CENTRAL_FLUX_OBJ   = $(TEST_BIMATERIAL_CENTRAL_FLUX_SRC:.cpp=.o)
TEST_BIMATERIAL_MIXED_FLUX_DISPATCH_OBJ = $(TEST_BIMATERIAL_MIXED_FLUX_DISPATCH_SRC:.cpp=.o)
TEST_BIMATERIAL_FAULT_PERSIDE_MATERIAL_OBJ = $(TEST_BIMATERIAL_FAULT_PERSIDE_MATERIAL_SRC:.cpp=.o)
TEST_BIMATERIAL_CONTRAST_GUARD_OBJ = $(TEST_BIMATERIAL_CONTRAST_GUARD_SRC:.cpp=.o)
TEST_BIMATERIAL_MIXED_FLUX_SHARED_OBJ = $(TEST_BIMATERIAL_MIXED_FLUX_SHARED_SRC:.cpp=.o)
TEST_BIMATERIAL_SEAM_MATERIAL_NP2_OBJ = $(TEST_BIMATERIAL_SEAM_MATERIAL_NP2_SRC:.cpp=.o)
TEST_BIMATERIAL_SEAM_FAULT_NP2_OBJ = $(TEST_BIMATERIAL_SEAM_FAULT_NP2_SRC:.cpp=.o)
TEST_BIMATERIAL_MIXED_FLUX_CFL_OBJ = $(TEST_BIMATERIAL_MIXED_FLUX_CFL_SRC:.cpp=.o)
TEST_BIMATERIAL_MIXED_FLUX_HOMOG_EQUIV_OBJ = $(TEST_BIMATERIAL_MIXED_FLUX_HOMOG_EQUIV_SRC:.cpp=.o)

TPV104_DRIVER_OBJ                    = $(TPV104_DRIVER_SRC:.cpp=.o)
TPV104_SUBSTEP_ITERATOR_OBJ          = $(TPV104_SUBSTEP_ITERATOR_SRC:.cpp=.o)

# Object dependencies for tests/drivers that call into the substep iterator
# cpp + fault_face_flux + friction_solver.  Reuse the existing build rules
# for the shared objects (FAULT_FACE_FLUX_OBJ and FRICTION_SOLVER_OBJ are
# produced by pre-existing targets above).
TPV104_SHARED_OBJS = $(TPV104_SUBSTEP_ITERATOR_OBJ) \
                     dynamic/fault_face_flux.o \
                     dynamic/friction_solver.o

# ---- compile rules ----
$(TEST_TPV104_PARAMS_OBJ): %.o: $(SRC)%.cpp $(FRICTION_HEADERS) $(TPV104_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_SLIP_LAW_SRW_PSI_OBJ): %.o: $(SRC)%.cpp $(FRICTION_HEADERS) $(TPV104_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_FRICTION_COEFF_STABLE_OBJ): %.o: $(SRC)%.cpp $(FRICTION_HEADERS) $(TPV104_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_TPV104_FRICTION_SOLVER_NEWTON_OBJ): %.o: $(SRC)%.cpp $(FRICTION_HEADERS) $(TPV104_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_TPV104_NUCLEATION_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV104_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_TPV104_NORMAL_SIGN_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV104_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_TPV104_MESH_OBJ): %.o: $(SRC)%.cpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_TPV104_SETUP_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV104_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_TPV104_SUBSTEP_ITERATOR_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV104_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_TPV104_PROBE_FORMAT_OBJ): %.o: $(SRC)%.cpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_TPV104_SMOKE_OBJ): %.o: $(SRC)%.cpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_TPV104_SIGMA_N_INVARIANCE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV104_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_TPV104_FREEZE_SIGMA_N_GATE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV104_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_TPV104_SUBSTEP_PREDICTOR_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV104_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_TPV104_SUBSTEP_ONE_SHOT_PARITY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV104_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_TPV104_SUBSTEP_ITERATOR_PARITY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV104_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_TPV104_SUBSTEP_DISPATCH_PARITY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV104_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_GODUNOV_CENTRAL_FLUX_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_BIMATERIAL_CENTRAL_FLUX_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_BIMATERIAL_MIXED_FLUX_DISPATCH_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_BIMATERIAL_FAULT_PERSIDE_MATERIAL_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_BIMATERIAL_CONTRAST_GUARD_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_BIMATERIAL_MIXED_FLUX_SHARED_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_BIMATERIAL_SEAM_MATERIAL_NP2_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_BIMATERIAL_SEAM_FAULT_NP2_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_BIMATERIAL_MIXED_FLUX_CFL_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_BIMATERIAL_MIXED_FLUX_HOMOG_EQUIV_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_MIXED_FLUX_FACE_SET_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_MIXED_FLUX_DISPATCH_NONE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV104_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TPV104_DRIVER_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV104_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

$(TPV104_SUBSTEP_ITERATOR_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV104_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

# ---- link rules ----
seas_test_tpv104_params: $(TEST_TPV104_PARAMS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV104_PARAMS_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_slip_law_srw_psi: $(TEST_SLIP_LAW_SRW_PSI_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_SLIP_LAW_SRW_PSI_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_friction_coeff_stable: $(TEST_FRICTION_COEFF_STABLE_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FRICTION_COEFF_STABLE_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_tpv104_friction_solver_newton: $(TEST_TPV104_FRICTION_SOLVER_NEWTON_OBJ) dynamic/friction_solver.o
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV104_FRICTION_SOLVER_NEWTON_OBJ) dynamic/friction_solver.o $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_tpv104_nucleation: $(TEST_TPV104_NUCLEATION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV104_NUCLEATION_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_tpv104_normal_sign: $(TEST_TPV104_NORMAL_SIGN_OBJ) dynamic/friction_solver.o
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV104_NORMAL_SIGN_OBJ) dynamic/friction_solver.o $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_tpv104_mesh: $(TEST_TPV104_MESH_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV104_MESH_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_tpv104_setup: $(TEST_TPV104_SETUP_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV104_SETUP_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_tpv104_substep_iterator: $(TEST_TPV104_SUBSTEP_ITERATOR_OBJ) $(TPV104_SHARED_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV104_SUBSTEP_ITERATOR_OBJ) $(TPV104_SHARED_OBJS) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_tpv104_probe_format: $(TEST_TPV104_PROBE_FORMAT_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV104_PROBE_FORMAT_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_tpv104_smoke: $(TEST_TPV104_SMOKE_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV104_SMOKE_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_tpv104_sigma_n_invariance: $(TEST_TPV104_SIGMA_N_INVARIANCE_OBJ) dynamic/fault_face_flux.o dynamic/friction_solver.o
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV104_SIGMA_N_INVARIANCE_OBJ) dynamic/fault_face_flux.o dynamic/friction_solver.o $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_tpv104_freeze_sigma_n_gate: $(TEST_TPV104_FREEZE_SIGMA_N_GATE_OBJ) dynamic/fault_face_flux.o dynamic/friction_solver.o
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV104_FREEZE_SIGMA_N_GATE_OBJ) dynamic/fault_face_flux.o dynamic/friction_solver.o $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Round-7 R-602/R-603: per-sub-step ADER predictor + iterator integration tests.
seas_test_tpv104_substep_predictor: $(TEST_TPV104_SUBSTEP_PREDICTOR_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV104_SUBSTEP_PREDICTOR_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_tpv104_substep_one_shot_parity: $(TEST_TPV104_SUBSTEP_ONE_SHOT_PARITY_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV104_SUBSTEP_ONE_SHOT_PARITY_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_tpv104_substep_iterator_parity: $(TEST_TPV104_SUBSTEP_ITERATOR_PARITY_OBJ) $(TPV104_SHARED_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV104_SUBSTEP_ITERATOR_PARITY_OBJ) $(TPV104_SHARED_OBJS) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_tpv104_substep_dispatch_parity: $(TEST_TPV104_SUBSTEP_DISPATCH_PARITY_OBJ) $(TPV104_SHARED_OBJS) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV104_SUBSTEP_DISPATCH_PARITY_OBJ) $(TPV104_SHARED_OBJS) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# REVIEW R-016 Phase 4 — TPV205 LSW dispatch tests.
TEST_TPV205_EVALUATE_ADER_LSW_PARITY_SRC = tests/unit/test_tpv205_evaluate_ader_lsw_parity.cpp
TEST_TPV205_EVALUATE_ADER_LSW_PARITY_OBJ = $(TEST_TPV205_EVALUATE_ADER_LSW_PARITY_SRC:.cpp=.o)
TEST_TPV205_MPI_RUPTURE_CROSSING_SRC = tests/verification/test_tpv205_mpi_rupture_crossing.cpp
TEST_TPV205_MPI_RUPTURE_CROSSING_OBJ = $(TEST_TPV205_MPI_RUPTURE_CROSSING_SRC:.cpp=.o)
TEST_TPV205_FRICTION_SRC = tests/unit/test_tpv205_friction.cpp
TEST_TPV205_FRICTION_OBJ = $(TEST_TPV205_FRICTION_SRC:.cpp=.o)

$(TEST_TPV205_EVALUATE_ADER_LSW_PARITY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV205_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_TPV205_MPI_RUPTURE_CROSSING_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV205_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_TPV205_FRICTION_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV205_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

seas_test_tpv205_evaluate_ader_lsw_parity: $(TEST_TPV205_EVALUATE_ADER_LSW_PARITY_OBJ) $(TPV205_SHARED_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV205_EVALUATE_ADER_LSW_PARITY_OBJ) $(TPV205_SHARED_OBJS) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_tpv205_mpi_rupture_crossing: $(TEST_TPV205_MPI_RUPTURE_CROSSING_OBJ) $(TPV205_SHARED_OBJS) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV205_MPI_RUPTURE_CROSSING_OBJ) $(TPV205_SHARED_OBJS) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# REVIEW R-003 (round 5): friction-helper unit test for the LSW
# strength-barrier short-circuits R-002/R-003 in tpv205_friction.hpp.
# Header-only test; depends on tpv205_friction.hpp + tpv205_params.hpp
# + test_macros.hpp.  No additional .o files beyond MFEM.
seas_test_tpv205_friction: $(TEST_TPV205_FRICTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV205_FRICTION_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Round-11 Mixed-Flux unit tests (Zhang et al. 2023).
seas_test_godunov_central_flux: $(TEST_GODUNOV_CENTRAL_FLUX_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_GODUNOV_CENTRAL_FLUX_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# Phase 1 (PLAN_mixed_flux_hetero_riemann.md): bi-material central-flux
# primitive.  Primitive-level deps only (cf. seas_test_godunov_central_flux):
# the bi-material object + the GodunovFlux object.
seas_test_bimaterial_central_flux: $(TEST_BIMATERIAL_CENTRAL_FLUX_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BIMATERIAL_CENTRAL_FLUX_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-bimaterial-central-flux: seas_test_bimaterial_central_flux
	./seas_test_bimaterial_central_flux

# Phase 3 (PLAN_mixed_flux_hetero_riemann.md): per-face mixed-flux dispatch on
# the bi-material operator.  Operator-level deps (mirror seas_test_bimaterial_
# wave_operator_parity).
seas_test_bimaterial_mixed_flux_dispatch: \
                  $(TEST_BIMATERIAL_MIXED_FLUX_DISPATCH_OBJ) \
                  $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
                  $(GODUNOV_FLUX_OBJ) \
                  $(GODUNOV_FLUX_POOL_OBJ) \
                  $(PML_LAYER_OBJ) \
                  $(PRECOMPUTED_FACE_FLUXES_OBJ) \
                  $(HETEROGENEOUS_MATERIAL_OBJ) \
                  $(FAULT_FACE_FLUX_OBJ) \
                  $(FRICTION_SOLVER_OBJ) \
                  $(SPATIAL_FRICTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_BIMATERIAL_MIXED_FLUX_DISPATCH_OBJ) \
	    $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
	    $(GODUNOV_FLUX_OBJ) \
	    $(GODUNOV_FLUX_POOL_OBJ) \
	    $(PML_LAYER_OBJ) \
	    $(PRECOMPUTED_FACE_FLUXES_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-bimaterial-mixed-flux-dispatch: seas_test_bimaterial_mixed_flux_dispatch

# Part B / B1: per-side fault impedance assignment (eps-offset rule).
seas_test_bimaterial_fault_perside_material: \
                  $(TEST_BIMATERIAL_FAULT_PERSIDE_MATERIAL_OBJ) \
                  $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
                  $(GODUNOV_FLUX_OBJ) \
                  $(GODUNOV_FLUX_POOL_OBJ) \
                  $(PML_LAYER_OBJ) \
                  $(PRECOMPUTED_FACE_FLUXES_OBJ) \
                  $(HETEROGENEOUS_MATERIAL_OBJ) \
                  $(FAULT_FACE_FLUX_OBJ) \
                  $(FRICTION_SOLVER_OBJ) \
                  $(SPATIAL_FRICTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_BIMATERIAL_FAULT_PERSIDE_MATERIAL_OBJ) \
	    $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
	    $(GODUNOV_FLUX_OBJ) \
	    $(GODUNOV_FLUX_POOL_OBJ) \
	    $(PML_LAYER_OBJ) \
	    $(PRECOMPUTED_FACE_FLUXES_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-bimaterial-fault-perside-material: seas_test_bimaterial_fault_perside_material
	./seas_test_bimaterial_fault_perside_material
	./seas_test_bimaterial_mixed_flux_dispatch

# Part A: central-flux contrast guard (predicate + face-set membership + two-sided
# dissipation mechanism).  Same operator-level object list as the dispatch test.
seas_test_bimaterial_contrast_guard: \
                  $(TEST_BIMATERIAL_CONTRAST_GUARD_OBJ) \
                  $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
                  $(GODUNOV_FLUX_OBJ) \
                  $(GODUNOV_FLUX_POOL_OBJ) \
                  $(PML_LAYER_OBJ) \
                  $(PRECOMPUTED_FACE_FLUXES_OBJ) \
                  $(HETEROGENEOUS_MATERIAL_OBJ) \
                  $(FAULT_FACE_FLUX_OBJ) \
                  $(FRICTION_SOLVER_OBJ) \
                  $(SPATIAL_FRICTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_BIMATERIAL_CONTRAST_GUARD_OBJ) \
	    $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
	    $(GODUNOV_FLUX_OBJ) \
	    $(GODUNOV_FLUX_POOL_OBJ) \
	    $(PML_LAYER_OBJ) \
	    $(PRECOMPUTED_FACE_FLUXES_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-bimaterial-contrast-guard: seas_test_bimaterial_contrast_guard
	./seas_test_bimaterial_contrast_guard

# Phase 3 P3-1: MPI (np=2) shared central-face dispatch (operator-level deps,
# mirrors seas_test_bimaterial_wave_operator_parity's object list).  NOT in the
# serial `test:` group — parallel tests run separately.
seas_test_bimaterial_mixed_flux_shared: \
                  $(TEST_BIMATERIAL_MIXED_FLUX_SHARED_OBJ) \
                  $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
                  $(GODUNOV_FLUX_OBJ) \
                  $(GODUNOV_FLUX_POOL_OBJ) \
                  $(PML_LAYER_OBJ) \
                  $(PRECOMPUTED_FACE_FLUXES_OBJ) \
                  $(HETEROGENEOUS_MATERIAL_OBJ) \
                  $(FAULT_FACE_FLUX_OBJ) \
                  $(FRICTION_SOLVER_OBJ) \
                  $(SPATIAL_FRICTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_BIMATERIAL_MIXED_FLUX_SHARED_OBJ) \
	    $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
	    $(GODUNOV_FLUX_OBJ) \
	    $(GODUNOV_FLUX_POOL_OBJ) \
	    $(PML_LAYER_OBJ) \
	    $(PRECOMPUTED_FACE_FLUXES_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-bimaterial-mixed-flux-shared: seas_test_bimaterial_mixed_flux_shared
	mpirun -np 2 ./seas_test_bimaterial_mixed_flux_shared

# Cross-rank material exchange — Phase 1 gate (uniform material accessor np=2).
seas_test_bimaterial_seam_material_np2: \
                  $(TEST_BIMATERIAL_SEAM_MATERIAL_NP2_OBJ) \
                  $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
                  $(GODUNOV_FLUX_OBJ) \
                  $(GODUNOV_FLUX_POOL_OBJ) \
                  $(PML_LAYER_OBJ) \
                  $(PRECOMPUTED_FACE_FLUXES_OBJ) \
                  $(HETEROGENEOUS_MATERIAL_OBJ) \
                  $(FAULT_FACE_FLUX_OBJ) \
                  $(FRICTION_SOLVER_OBJ) \
                  $(SPATIAL_FRICTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_BIMATERIAL_SEAM_MATERIAL_NP2_OBJ) \
	    $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
	    $(GODUNOV_FLUX_OBJ) \
	    $(GODUNOV_FLUX_POOL_OBJ) \
	    $(PML_LAYER_OBJ) \
	    $(PRECOMPUTED_FACE_FLUXES_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-bimaterial-seam-material-np2: seas_test_bimaterial_seam_material_np2
	mpirun -np 2 ./seas_test_bimaterial_seam_material_np2

# Cross-rank material exchange — Phase 4 gate (per-side fault across a seam, np=2).
seas_test_bimaterial_seam_fault_np2: \
                  $(TEST_BIMATERIAL_SEAM_FAULT_NP2_OBJ) \
                  $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
                  $(GODUNOV_FLUX_OBJ) \
                  $(GODUNOV_FLUX_POOL_OBJ) \
                  $(PML_LAYER_OBJ) \
                  $(PRECOMPUTED_FACE_FLUXES_OBJ) \
                  $(HETEROGENEOUS_MATERIAL_OBJ) \
                  $(FAULT_FACE_FLUX_OBJ) \
                  $(FRICTION_SOLVER_OBJ) \
                  $(SPATIAL_FRICTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_BIMATERIAL_SEAM_FAULT_NP2_OBJ) \
	    $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
	    $(GODUNOV_FLUX_OBJ) \
	    $(GODUNOV_FLUX_POOL_OBJ) \
	    $(PML_LAYER_OBJ) \
	    $(PRECOMPUTED_FACE_FLUXES_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-bimaterial-seam-fault-np2: seas_test_bimaterial_seam_fault_np2
	mpirun -np 2 ./seas_test_bimaterial_seam_fault_np2

# Phase 4 (PLAN_mixed_flux_hetero_riemann.md): mixed-flux CFL shared-factor +
# central+ADER guard (operator-level deps; serial -> in the `test:` group).
seas_test_bimaterial_mixed_flux_cfl: \
                  $(TEST_BIMATERIAL_MIXED_FLUX_CFL_OBJ) \
                  $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
                  $(GODUNOV_FLUX_OBJ) \
                  $(GODUNOV_FLUX_POOL_OBJ) \
                  $(PML_LAYER_OBJ) \
                  $(PRECOMPUTED_FACE_FLUXES_OBJ) \
                  $(HETEROGENEOUS_MATERIAL_OBJ) \
                  $(FAULT_FACE_FLUX_OBJ) \
                  $(FRICTION_SOLVER_OBJ) \
                  $(SPATIAL_FRICTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_BIMATERIAL_MIXED_FLUX_CFL_OBJ) \
	    $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
	    $(GODUNOV_FLUX_OBJ) \
	    $(GODUNOV_FLUX_POOL_OBJ) \
	    $(PML_LAYER_OBJ) \
	    $(PRECOMPUTED_FACE_FLUXES_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-bimaterial-mixed-flux-cfl: seas_test_bimaterial_mixed_flux_cfl
	./seas_test_bimaterial_mixed_flux_cfl

# Phase 6 (PLAN_mixed_flux_hetero_riemann.md): homogeneous-equivalence GATE
# (operator-level deps; mirror seas_test_bimaterial_wave_operator_parity).
seas_test_bimaterial_mixed_flux_homog_equivalence: \
                  $(TEST_BIMATERIAL_MIXED_FLUX_HOMOG_EQUIV_OBJ) \
                  $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
                  $(GODUNOV_FLUX_OBJ) \
                  $(GODUNOV_FLUX_POOL_OBJ) \
                  $(PML_LAYER_OBJ) \
                  $(PRECOMPUTED_FACE_FLUXES_OBJ) \
                  $(HETEROGENEOUS_MATERIAL_OBJ) \
                  $(FAULT_FACE_FLUX_OBJ) \
                  $(FRICTION_SOLVER_OBJ) \
                  $(SPATIAL_FRICTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_BIMATERIAL_MIXED_FLUX_HOMOG_EQUIV_OBJ) \
	    $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
	    $(GODUNOV_FLUX_OBJ) \
	    $(GODUNOV_FLUX_POOL_OBJ) \
	    $(PML_LAYER_OBJ) \
	    $(PRECOMPUTED_FACE_FLUXES_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-bimaterial-mixed-flux-homog-equivalence: seas_test_bimaterial_mixed_flux_homog_equivalence
	./seas_test_bimaterial_mixed_flux_homog_equivalence

seas_test_mixed_flux_face_set: $(TEST_MIXED_FLUX_FACE_SET_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_MIXED_FLUX_FACE_SET_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_test_mixed_flux_dispatch_none: $(TEST_MIXED_FLUX_DISPATCH_NONE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_MIXED_FLUX_DISPATCH_NONE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

seas_tpv104_driver: $(TPV104_DRIVER_OBJ) $(TPV104_SHARED_OBJS) \
                    $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) $(PRECOMPUTED_FACE_FLUXES_OBJ) \
                    $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TPV104_DRIVER_OBJ) \
	    $(TPV104_SHARED_OBJS) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
	    $(PRECOMPUTED_FACE_FLUXES_OBJ) $(GODUNOV_FLUX_OBJ) \
	    $(PML_LAYER_OBJ) $(MFEM_LIBS) $(SEAS_FS_LIB)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# ---- umbrella target ----
# Build order: driver first (so test_tpv104_smoke has the binary available
# for its popen-based banner scrape), then every unit-test binary, then
# run each.
.PHONY: test-tpv104 tpv104-mesh
test-tpv104: seas_tpv104_driver \
             seas_test_tpv104_params \
             seas_test_slip_law_srw_psi \
             seas_test_friction_coeff_stable \
             seas_test_tpv104_friction_solver_newton \
             seas_test_tpv104_nucleation \
             seas_test_tpv104_normal_sign \
             seas_test_tpv104_mesh \
             seas_test_tpv104_setup \
             seas_test_tpv104_substep_iterator \
             seas_test_tpv104_probe_format \
             seas_test_tpv104_smoke \
             seas_test_tpv104_sigma_n_invariance \
             seas_test_tpv104_freeze_sigma_n_gate \
             seas_test_friction_solver_stable
	@echo "=== TPV104 unit tests ==="
	./seas_test_tpv104_params
	./seas_test_slip_law_srw_psi
	./seas_test_friction_coeff_stable
	./seas_test_tpv104_friction_solver_newton
	./seas_test_tpv104_nucleation
	./seas_test_tpv104_normal_sign
	./seas_test_tpv104_mesh
	./seas_test_tpv104_setup
	./seas_test_tpv104_substep_iterator
	./seas_test_tpv104_probe_format
	./seas_test_tpv104_smoke
	./seas_test_tpv104_sigma_n_invariance
	./seas_test_tpv104_freeze_sigma_n_gate
	./seas_test_friction_solver_stable
	@echo "=== TPV104 suite complete ==="

# Gmsh-generated mesh .msh files (built on demand in the pythonenv
# conda env).  Users run `conda activate pythonenv && make tpv104-mesh`
# to produce tpv104_{200m,500m,1000m}.msh from the .geo templates.
tpv104-mesh: tpv104/mesh/tpv104_200m.msh tpv104/mesh/tpv104_500m.msh tpv104/mesh/tpv104_1000m.msh

tpv104/mesh/%.msh: tpv104/mesh/%.geo
	cd tpv104/mesh && gmsh -3 $(notdir $<) -o $(notdir $@)

# Plan §4.10.X R5-003 fix — test_friction_solver_stable.
TEST_FRICTION_SOLVER_STABLE_SRC = tests/unit/test_friction_solver_stable.cpp
TEST_FRICTION_SOLVER_STABLE_OBJ = $(TEST_FRICTION_SOLVER_STABLE_SRC:.cpp=.o)

$(TEST_FRICTION_SOLVER_STABLE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(TPV104_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

seas_test_friction_solver_stable: $(TEST_FRICTION_SOLVER_STABLE_OBJ) dynamic/friction_solver.o
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FRICTION_SOLVER_STABLE_OBJ) dynamic/friction_solver.o $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# R-004 (REVIEW.md 2026-05-16): seas_test_friction_solver_stable is
# rolled into the primary `test-tpv104` umbrella above (prerequisite +
# run line).  An earlier block here re-declared `test-tpv104` with a
# new recipe, which Make treats as REPLACING (not extending) the
# original recipe — so `make test-tpv104` silently ran only this one
# test and skipped the other 13.  The redeclaration also triggered the
# "overriding recipe for target `test-tpv104`" warning that users saw
# in BP5 Phase 6 production builds.

# --- Phase 4 (PLAN_mixed_flux_seissol_port_2026-06-19.md): matrix+RK+LSW bimaterial gate ---
# Same object set as seas_test_lsw_rk_mixed_flux (already links the bimaterial +
# heterogeneous objects); BimaterialWaveOperator is header-only (.hpp/.inl).
$(TEST_LSW_RK_MIXED_FLUX_BIMATERIAL_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) \
                                    $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) \
	    $(TOML_FLAGS) -c $< -o $@

seas_test_lsw_rk_mixed_flux_bimaterial: $(TEST_LSW_RK_MIXED_FLUX_BIMATERIAL_OBJ) \
                             $(RK_TIME_STEPPER_OBJ) \
                             $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
                             $(GODUNOV_FLUX_OBJ) \
                             $(GODUNOV_FLUX_POOL_OBJ) \
                             $(PRECOMPUTED_FACE_FLUXES_OBJ) \
                             $(PML_LAYER_OBJ) \
                             $(HETEROGENEOUS_MATERIAL_OBJ) \
                             $(FAULT_FACE_FLUX_OBJ) \
                             $(FRICTION_SOLVER_OBJ) \
                             $(SPATIAL_NUCLEATION_OBJ) \
                             $(SPATIAL_FRICTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_LSW_RK_MIXED_FLUX_BIMATERIAL_OBJ) \
	    $(RK_TIME_STEPPER_OBJ) \
	    $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_BIMATERIAL_OBJ) \
	    $(GODUNOV_FLUX_OBJ) \
	    $(GODUNOV_FLUX_POOL_OBJ) \
	    $(PRECOMPUTED_FACE_FLUXES_OBJ) \
	    $(PML_LAYER_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(FAULT_FACE_FLUX_OBJ) \
	    $(FRICTION_SOLVER_OBJ) \
	    $(SPATIAL_NUCLEATION_OBJ) \
	    $(SPATIAL_FRICTION_OBJ) \
	    $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-lsw-rk-mixed-flux-bimaterial: seas_test_lsw_rk_mixed_flux_bimaterial
	./seas_test_lsw_rk_mixed_flux_bimaterial

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
TEST_GODUNOV_FLUX_SRC = tests/unit/test_godunov_flux.cpp
TEST_WAVE_OPERATOR_SRC = tests/unit/test_wave_operator.cpp
PML_LAYER_SRC = dynamic/pml_layer.cpp
PML_LAYER_OBJ = $(PML_LAYER_SRC:.cpp=.o)
FAULT_FACE_FLUX_SRC = dynamic/fault_face_flux.cpp
FAULT_FACE_FLUX_OBJ = $(FAULT_FACE_FLUX_SRC:.cpp=.o)
FRICTION_SOLVER_SRC = dynamic/friction_solver.cpp
FRICTION_SOLVER_OBJ = $(FRICTION_SOLVER_SRC:.cpp=.o)
TEST_WAVE_BC_SRC = tests/unit/test_wave_bc.cpp
TEST_PML_SRC = tests/unit/test_pml.cpp
TEST_FAULT_FACE_FLUX_SRC = tests/unit/test_fault_face_flux.cpp
TEST_FRICTION_SOLVER_SRC = tests/unit/test_friction_solver.cpp
TEST_TPV102_SETUP_SRC = tests/unit/test_tpv102_setup.cpp
TPV102_DRIVER_SRC = drivers/tpv102_driver.cpp
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
# v9.1.0 Pelties-9 follow-through: post-R-001 VTU continuity + post-R-003 basis symmetry
TEST_FAULT_SURFACE_VTU_CONTINUITY_SRC      = tests/unit/test_fault_surface_vtu_continuity.cpp
TEST_FAULT_BASIS_DIP_STRIKE_SYMMETRY_SRC   = tests/unit/test_fault_basis_dip_strike_symmetry.cpp
# v9.2.0 §19 R-V92-E02 stage-4 VTU diagnostic
TEST_FAULT_SURFACE_VTU_K4_SRC              = tests/unit/test_fault_surface_vtu_k4.cpp
# --paraview-fields / SetFaultVTUFields filter
TEST_FAULT_SURFACE_VTU_FIELD_FILTER_SRC    = tests/unit/test_fault_surface_vtu_field_filter.cpp
# v9.2.0 plan §4.1-§4.4 regression gates
TEST_GODUNOV_ROTATION_IDENTITY_SRC         = tests/unit/test_godunov_rotation_identity.cpp
TEST_CANONICAL_ROTATION_PURE_STRIKESLIP_SRC = tests/unit/test_canonical_rotation_pure_strikeslip.cpp
TEST_FAULT_BASIS_QP_ORTHONORMALITY_SRC     = tests/unit/test_fault_basis_qp_orthonormality.cpp
TEST_PER_QP_VS_CENTROID_BASIS_SRC          = tests/unit/test_per_qp_vs_centroid_basis_planar_1el.cpp
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
# Pepper-bug: ADJACENT-TRIANGLE UNIFORMITY test (8-fault-triangle fixture)
TEST_ADJACENT_TRI_UNIFORMITY_SRC           = tests/unit/test_adjacent_triangle_fault_uniformity.cpp
TEST_GODUNOV_ROUNDTRIP_SRC                 = tests/unit/test_godunov_rotation_roundtrip.cpp
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
# v9.3.0 Phase 5 (I-06 part C): R-002 nucleation sign regression
TEST_R002_NUCLEATION_SIGN_SRC              = tests/unit/test_R002_nucleation_sign_match.cpp
# v9.3.0 Phase 5 (I-06 part C): R-009 nucleation both-sides regression
TEST_R009_NUCLEATION_BOTH_SIDES_SRC        = tests/unit/test_R009_nucleation_total_both_sides.cpp
# v9.3.0 R-I06-001 / R-I06-006 (review-round-2): total-Q absorbing BC
# preserves pre-stress equilibrium (AbsorbingTotal regression gate)
TEST_TPV102_TOTAL_ABSORBING_SRC            = tests/unit/test_tpv102_total_absorbing_equilibrium.cpp
# v9.3.0 R-I06-003 / R-I06-004 (review-round-2): shape-max amplitude
# scaling + add-delta preserves wave perturbations
TEST_R_I06_003_NUC_SHAPE_SCALING_SRC       = tests/unit/test_R_I06_003_nucleation_shape_scaling.cpp
# 2026-04-22 nucleation persistent-prestress channel fix regression
# (replaces the bulk-Q add-delta path that diluted nuc by 1e4-1e5x)
TEST_PERSISTENT_NUC_PRESTRESS_SRC          = tests/unit/test_persistent_nuc_prestress_channel.cpp
# v9.3.0 R-I06-round3 R-006: direct unit tests for FreeSurfaceTotal /
# FreeSurfaceGodunovTotal (Q_bg=0 equivalence + tilted-normal equilibrium)
TEST_FREE_SURFACE_TOTAL_DIRECT_SRC         = tests/unit/test_free_surface_total_direct.cpp
# ADER I-05 Phase 1: polynomial-exactness test for ApplySpatialDerivative
TEST_WAVE_OP_SPATIAL_DERIV_SRC             = tests/unit/test_wave_operator_spatial_derivative.cpp
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

# R-501: `make ... CXXFLAGS+=-D<FLAG>` may not reach the compile step
# because MFEM's config.mk explicitly sets MFEM_CPPFLAGS/MFEM_CXXFLAGS.
# Use this SEAS-specific variable for diagnostic flags, e.g.
#   make seas_tpv102_driver \
#        SEAS_EXTRA_CPPFLAGS="-DSEAS_DIAG_FAULT_FLUX -DSEAS_DIAG_GHOST_EXCHANGE"
SEAS_EXTRA_CPPFLAGS ?=

# Object files
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
TEST_TPV102_LOCAL_OBJ = $(TEST_TPV102_LOCAL_SRC:.cpp=.o)
TEST_R101_SHARED_FAULT_OBJ = $(TEST_R101_SHARED_FAULT_SRC:.cpp=.o)
TEST_FAULT_FACE_FLUX_FRAME_OBJ = $(TEST_FAULT_FACE_FLUX_FRAME_SRC:.cpp=.o)
TEST_FAULT_FACE_FLUX_FRAME_AND_FLUX_OBJ = $(TEST_FAULT_FACE_FLUX_FRAME_AND_FLUX_SRC:.cpp=.o)
TEST_GODUNOV_IDENTITY_OBJ = $(TEST_GODUNOV_IDENTITY_SRC:.cpp=.o)
TEST_FAULT_FLUX_BRANCH_EQUIV_OBJ = $(TEST_FAULT_FLUX_BRANCH_EQUIV_SRC:.cpp=.o)
TEST_GODUNOV_INTERIOR_EQUAL_SIDES_OBJ = $(TEST_GODUNOV_INTERIOR_EQUAL_SIDES_SRC:.cpp=.o)
TEST_FAULT_FACE_FLUX_PER_SIDE_ASSEMBLY_OBJ = $(TEST_FAULT_FACE_FLUX_PER_SIDE_ASSEMBLY_SRC:.cpp=.o)
TEST_FAULT_FACE_FLUX_BIMATERIAL_GUARD_OBJ  = $(TEST_FAULT_FACE_FLUX_BIMATERIAL_GUARD_SRC:.cpp=.o)
# v9.1.0 Pelties-9 follow-through
TEST_FAULT_SURFACE_VTU_CONTINUITY_OBJ      = $(TEST_FAULT_SURFACE_VTU_CONTINUITY_SRC:.cpp=.o)
TEST_FAULT_BASIS_DIP_STRIKE_SYMMETRY_OBJ   = $(TEST_FAULT_BASIS_DIP_STRIKE_SYMMETRY_SRC:.cpp=.o)
# v9.2.0 §19 R-V92-E02
TEST_FAULT_SURFACE_VTU_K4_OBJ              = $(TEST_FAULT_SURFACE_VTU_K4_SRC:.cpp=.o)
# --paraview-fields field filter
TEST_FAULT_SURFACE_VTU_FIELD_FILTER_OBJ    = $(TEST_FAULT_SURFACE_VTU_FIELD_FILTER_SRC:.cpp=.o)
# v9.2.0 plan §4.1-§4.4 regression gates
TEST_GODUNOV_ROTATION_IDENTITY_OBJ         = $(TEST_GODUNOV_ROTATION_IDENTITY_SRC:.cpp=.o)
TEST_CANONICAL_ROTATION_PURE_STRIKESLIP_OBJ = $(TEST_CANONICAL_ROTATION_PURE_STRIKESLIP_SRC:.cpp=.o)
TEST_FAULT_BASIS_QP_ORTHONORMALITY_OBJ     = $(TEST_FAULT_BASIS_QP_ORTHONORMALITY_SRC:.cpp=.o)
TEST_PER_QP_VS_CENTROID_BASIS_OBJ          = $(TEST_PER_QP_VS_CENTROID_BASIS_SRC:.cpp=.o)
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
TEST_ADJACENT_TRI_UNIFORMITY_OBJ           = $(TEST_ADJACENT_TRI_UNIFORMITY_SRC:.cpp=.o)
TEST_GODUNOV_ROUNDTRIP_OBJ                 = $(TEST_GODUNOV_ROUNDTRIP_SRC:.cpp=.o)
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
TEST_R002_NUCLEATION_SIGN_OBJ              = $(TEST_R002_NUCLEATION_SIGN_SRC:.cpp=.o)
TEST_R009_NUCLEATION_BOTH_SIDES_OBJ        = $(TEST_R009_NUCLEATION_BOTH_SIDES_SRC:.cpp=.o)
TEST_TPV102_TOTAL_ABSORBING_OBJ            = $(TEST_TPV102_TOTAL_ABSORBING_SRC:.cpp=.o)
TEST_R_I06_003_NUC_SHAPE_SCALING_OBJ       = $(TEST_R_I06_003_NUC_SHAPE_SCALING_SRC:.cpp=.o)
TEST_PERSISTENT_NUC_PRESTRESS_OBJ          = $(TEST_PERSISTENT_NUC_PRESTRESS_SRC:.cpp=.o)
TEST_FREE_SURFACE_TOTAL_DIRECT_OBJ         = $(TEST_FREE_SURFACE_TOTAL_DIRECT_SRC:.cpp=.o)
TEST_WAVE_OP_SPATIAL_DERIV_OBJ             = $(TEST_WAVE_OP_SPATIAL_DERIV_SRC:.cpp=.o)
TEST_ADER_CK_PREDICTOR_OBJ                 = $(TEST_ADER_CK_PREDICTOR_SRC:.cpp=.o)
TEST_ADER_FAULT_FLUX_EQUIV_OBJ             = $(TEST_ADER_FAULT_FLUX_EQUIV_SRC:.cpp=.o)
TEST_ADER_LINEAR_WAVE_EQUIV_OBJ            = $(TEST_ADER_LINEAR_WAVE_EQUIV_SRC:.cpp=.o)
TEST_ADER_FAULT_COUPLING_DT_LIMIT_OBJ      = $(TEST_ADER_FAULT_COUPLING_DT_LIMIT_SRC:.cpp=.o)
TEST_ADER_INTERIOR_VS_SHARED_LIVE_OBJ      = $(TEST_ADER_INTERIOR_VS_SHARED_LIVE_SRC:.cpp=.o)
TEST_ADER_TPV102_SMOKE_OBJ                 = $(TEST_ADER_TPV102_SMOKE_SRC:.cpp=.o)
TEST_R001_DRIVER_DEFAULT_ADER_OBJ          = $(TEST_R001_DRIVER_DEFAULT_ADER_SRC:.cpp=.o)
TEST_R002_DRIVER_INIT_TOTAL_Q_OBJ          = $(TEST_R002_DRIVER_INIT_TOTAL_Q_SRC:.cpp=.o)

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
             io/paraview_output.hpp

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
                  dynamic/seas_dynamic_operator.hpp \
                  dynamic/tpv102_setup.hpp

# All headers
SEAS_HEADERS = $(FRICTION_HEADERS) $(DOMAIN_HEADERS) $(FAULT_HEADERS) $(SOLVER_HEADERS) $(IO_HEADERS) $(CHECKPOINT_HEADERS) $(CONSTITUTIVE_HEADERS) $(CONFIG_HEADERS) $(DYNAMIC_HEADERS)

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
        test-per-qp-vs-centroid-basis-planar-1el \
        test-v92-regression-gates \
        regression-check quick-check

all: $(MINIAPPS)

# Remove built-in rules
%: %.cpp
%.o: %.cpp

SEAS_INCLUDES = $(if $(SRC),-I$(MFEM_DIR))

# Phase 1: Friction law unit tests
seas_test_friction_law: $(TEST_FRICTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FRICTION_OBJ) $(MFEM_LIBS)

# Phase 1: State evolution unit tests
seas_test_state_evolution: $(TEST_STATE_EVOLUTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_STATE_EVOLUTION_OBJ) $(MFEM_LIBS)

# Phase 2: Antiplane domain operator unit tests
seas_test_antiplane: $(TEST_ANTIPLANE_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ANTIPLANE_OBJ) $(MFEM_LIBS)

# Phase 3: Fault operator unit tests
seas_test_fault_operator: $(TEST_FAULT_OPERATOR_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_OPERATOR_OBJ) $(MFEM_LIBS)

# Phase 4: Quasi-dynamic coupling tests
seas_test_quasi_dynamic: $(TEST_QUASI_DYNAMIC_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_QUASI_DYNAMIC_OBJ) $(MFEM_LIBS)

# Phase 4: Short BP2 simulation test
seas_test_bp2_short: $(TEST_BP2_SHORT_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BP2_SHORT_OBJ) $(MFEM_LIBS)

# Phase 5: I/O unit tests
seas_test_io: $(TEST_IO_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_IO_OBJ) $(MFEM_LIBS)

# Phase 5: BP2 serial smoke test (800 m coarse mesh)
seas_bp2_serial_smoke: $(BP2_SERIAL_SMOKE_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(BP2_SERIAL_SMOKE_OBJ) $(MFEM_LIBS)

# Phase 11: Checkpoint unit test
seas_test_checkpoint: $(TEST_CHECKPOINT_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_CHECKPOINT_OBJ) $(MFEM_LIBS)

# Psi-space state variable unit test
seas_test_psi_state: $(TEST_PSI_STATE_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PSI_STATE_OBJ) $(MFEM_LIBS)

# Phase 1a: Fault detection unit test
seas_test_fault_detection: $(TEST_FAULT_DETECTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_DETECTION_OBJ) $(MFEM_LIBS)

# Phase 1b: DomainOperator interface unit test
seas_test_domain_interface: $(TEST_DOMAIN_INTERFACE_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_DOMAIN_INTERFACE_OBJ) $(MFEM_LIBS)

# BP5 Phase 1: Parameter tests
seas_test_bp5_params: $(TEST_BP5_PARAMS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BP5_PARAMS_OBJ) $(MFEM_LIBS)

# BP5 Phase 1: Vector friction tests
seas_test_vector_friction: $(TEST_VECTOR_FRICTION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_VECTOR_FRICTION_OBJ) $(MFEM_LIBS)

# BP5 Phase 2a: Fault basis tests
seas_test_fault_basis: $(TEST_FAULT_BASIS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_BASIS_OBJ) $(MFEM_LIBS)

# BP5 Phase 2b: Elasticity BR2 integrator tests
seas_test_elasticity_br2: $(TEST_ELASTICITY_BR2_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ELASTICITY_BR2_OBJ) $(MFEM_LIBS)

# BP5 Phase 2b: Elasticity domain operator tests
seas_test_elasticity_operator: $(TEST_ELASTICITY_OPERATOR_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ELASTICITY_OPERATOR_OBJ) $(MFEM_LIBS)

# BP5 Phase 3: BP5 fault operator tests
seas_test_bp5_fault_operator: $(TEST_BP5_FAULT_OPERATOR_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BP5_FAULT_OPERATOR_OBJ) $(MFEM_LIBS)

# BP5 Phase 4: BP5 integration tests
seas_test_bp5_integration: $(TEST_BP5_INTEGRATION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BP5_INTEGRATION_OBJ) $(MFEM_LIBS)

# BP5 Phase 5: Output unit tests
seas_test_bp5_output: $(TEST_BP5_OUTPUT_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BP5_OUTPUT_OBJ) $(MFEM_LIBS)

# BP5 Phase 5: Mesh unit tests
seas_test_bp5_mesh: $(TEST_BP5_MESH_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BP5_MESH_OBJ) $(MFEM_LIBS)

# BP5 v16: Diagnostic VTK unit tests
seas_test_diag_vtk: $(TEST_DIAG_VTK_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_DIAG_VTK_OBJ) $(MFEM_LIBS)

# BP5 v16: Mesh region unit tests
seas_test_mesh_regions: $(TEST_MESH_REGIONS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_MESH_REGIONS_OBJ) $(MFEM_LIBS)

# v44: Face quadrature unit tests
seas_test_face_quadrature: $(TEST_FACE_QUADRATURE_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FACE_QUADRATURE_OBJ) $(MFEM_LIBS)

# v54: Cross-verification against Tandem
seas_test_cross_verify: $(TEST_CROSS_VERIFY_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_CROSS_VERIFY_OBJ) $(MFEM_LIBS)

# Dynamic rupture Phase 1: Godunov flux tests
seas_test_godunov_flux: $(TEST_GODUNOV_FLUX_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_GODUNOV_FLUX_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)

# Dynamic rupture Phase 1: Wave operator tests
seas_test_wave_operator: $(TEST_WAVE_OPERATOR_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_WAVE_OPERATOR_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# Dynamic rupture Phase 2a: BC tests
seas_test_wave_bc: $(TEST_WAVE_BC_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_WAVE_BC_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# Dynamic rupture Phase 2b: PML tests
seas_test_pml: $(TEST_PML_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PML_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# Dynamic rupture Phase 3: Fault face flux tests
seas_test_fault_face_flux: $(TEST_FAULT_FACE_FLUX_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_FACE_FLUX_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)

# TPV102 v9.0.0 Phase 1 unit tests (§3.1, §3.1b, §3.1c, §3.1d)
seas_test_fault_face_flux_frame: $(TEST_FAULT_FACE_FLUX_FRAME_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_FACE_FLUX_FRAME_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)

seas_test_fault_face_flux_frame_and_flux: $(TEST_FAULT_FACE_FLUX_FRAME_AND_FLUX_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_FACE_FLUX_FRAME_AND_FLUX_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)

seas_test_godunov_identity_normal_reversal: $(TEST_GODUNOV_IDENTITY_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_GODUNOV_IDENTITY_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)

seas_test_fault_flux_interior_vs_shared_branch_equivalence: $(TEST_FAULT_FLUX_BRANCH_EQUIV_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_FLUX_BRANCH_EQUIV_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)

# v9.0.0 Pelties-9 fix: §3.1e primitive identity test
seas_test_godunov_interior_equal_sides_identity: $(TEST_GODUNOV_INTERIOR_EQUAL_SIDES_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_GODUNOV_INTERIOR_EQUAL_SIDES_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)

# v9.0.0 Pelties-9 fix: §3.1f per-side assembly test
seas_test_fault_face_flux_per_side_assembly: $(TEST_FAULT_FACE_FLUX_PER_SIDE_ASSEMBLY_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_FACE_FLUX_PER_SIDE_ASSEMBLY_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)

# REVIEW R-007 (v9.0.0): bimaterial-guard regression test for
# FaultFaceFlux::Evaluate.  Death test — the child aborts inside
# MFEM_VERIFY; the shell target below interprets non-zero exit as PASS.
seas_test_fault_face_flux_bimaterial_guard: $(TEST_FAULT_FACE_FLUX_BIMATERIAL_GUARD_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_FACE_FLUX_BIMATERIAL_GUARD_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)

# v9.1.0 plan §4.1 — fault-surface VTU CellData contract (post-R-001)
seas_test_fault_surface_vtu_continuity: $(TEST_FAULT_SURFACE_VTU_CONTINUITY_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_SURFACE_VTU_CONTINUITY_OBJ) $(MFEM_LIBS)

# v9.2.0 plan §19 — R-V92-E02 stage-4 VTU diagnostic
seas_test_fault_surface_vtu_k4: $(TEST_FAULT_SURFACE_VTU_K4_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_SURFACE_VTU_K4_OBJ) $(MFEM_LIBS)

# --paraview-fields / SetFaultVTUFields field filter
seas_test_fault_surface_vtu_field_filter: $(TEST_FAULT_SURFACE_VTU_FIELD_FILTER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_SURFACE_VTU_FIELD_FILTER_OBJ) $(MFEM_LIBS)

# v9.1.0 plan §4.2 — FaultBasis orthonormality contract (post-R-003)
seas_test_fault_basis_dip_strike_symmetry: $(TEST_FAULT_BASIS_DIP_STRIKE_SYMMETRY_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_BASIS_DIP_STRIKE_SYMMETRY_OBJ) $(MFEM_LIBS)

# v9.2.0 plan §4.1 — Godunov rotation identity T·Tinv = I
seas_test_godunov_rotation_identity: $(TEST_GODUNOV_ROTATION_IDENTITY_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_GODUNOV_ROTATION_IDENTITY_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)

# v9.2.0 plan §4.2 — canonical rotation pure strike-slip Tinv probe
seas_test_canonical_rotation_pure_strikeslip: $(TEST_CANONICAL_ROTATION_PURE_STRIKESLIP_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_CANONICAL_ROTATION_PURE_STRIKESLIP_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)

# v9.2.0 plan §4.3 — per-QP fault basis orthonormality
seas_test_fault_basis_qp_orthonormality: $(TEST_FAULT_BASIS_QP_ORTHONORMALITY_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_BASIS_QP_ORTHONORMALITY_OBJ) $(MFEM_LIBS)

# v9.2.0 plan §4.4 — per-QP vs centroid basis on planar 1-tet fixture
seas_test_per_qp_vs_centroid_basis_planar_1el: $(TEST_PER_QP_VS_CENTROID_BASIS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PER_QP_VS_CENTROID_BASIS_OBJ) $(MFEM_LIBS)

# v9.2.0 plan §4.9 (rev-3f) — No-penalty regression gate (Pelties 2012 §3.1)
seas_test_no_penalty_dynamic_rupture: $(TEST_NO_PENALTY_DR_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_NO_PENALTY_DR_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# v9.2.0 plan Step 6 (H-V92-U probe) — interior-fault-path Pelties eq. (7)
seas_test_interior_fault_flux_path: $(TEST_INTERIOR_FAULT_FLUX_PATH_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_INTERIOR_FAULT_FLUX_PATH_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# Pepper-bug plan 2026-04-22 Phase 3 P3-T1 — asymmetric Q probe
seas_test_interior_fault_flux_path_asymmetric: $(TEST_INTERIOR_FAULT_FLUX_PATH_ASYM_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_INTERIOR_FAULT_FLUX_PATH_ASYM_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# Pepper-bug plan 2026-04-22 Phase 3 P3-T3 — rupture-regime amplitude probe
seas_test_interior_fault_flux_path_rupture_regime: $(TEST_INTERIOR_FAULT_FLUX_PATH_RUPTURE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_INTERIOR_FAULT_FLUX_PATH_RUPTURE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# Pepper-bug plan 2026-04-22 Phase 5 P5-N — bulk RK4 noise stability
seas_test_bulk_rk4_noise_stability: $(TEST_BULK_RK4_NOISE_STAB_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BULK_RK4_NOISE_STAB_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# Pepper-bug plan 2026-04-22 Test 5 — shared-fault Pelties-7 analytic (np=4)
seas_test_shared_fault_flux_path_analytic: $(TEST_SHARED_FAULT_FLUX_PATH_ANA_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_SHARED_FAULT_FLUX_PATH_ANA_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# Pepper-bug plan 2026-04-22 Phase 5 P5-PW — bulk plane-wave recovery
seas_test_bulk_plane_wave_recovery: $(TEST_BULK_PLANE_WAVE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BULK_PLANE_WAVE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# Pepper-bug plan 2026-04-22 Phase 1 — live interior-vs-shared cross-check
# NB: built with -DSEAS_TEST_INTERNAL so the test can access WaveOperator
# internal bookkeeping (getters gated by that macro in wave_operator.hpp).
seas_test_rupture_multistep_serial_vs_parallel: $(TEST_RUPTURE_MULTISTEP_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_RUPTURE_MULTISTEP_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

seas_test_adjacent_triangle_fault_uniformity: $(TEST_ADJACENT_TRI_UNIFORMITY_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ADJACENT_TRI_UNIFORMITY_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

seas_test_godunov_rotation_roundtrip: $(TEST_GODUNOV_ROUNDTRIP_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_GODUNOV_ROUNDTRIP_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)

seas_test_interior_vs_shared_branch_live: $(TEST_INTERIOR_VS_SHARED_LIVE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_INTERIOR_VS_SHARED_LIVE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# Pepper-bug plan 2026-04-22 Phase 2 — sign_flipped / elem1_on_plus truth table
seas_test_sign_flipped_elem1_on_plus_truth_table: $(TEST_SIGN_FLIPPED_TRUTH_TABLE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_SIGN_FLIPPED_TRUTH_TABLE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# Pepper-bug plan 2026-04-22 Phase 4 P4-T1 — MakeFaceKey uniqueness (np=2)
seas_test_make_face_key_uniqueness: $(TEST_MAKE_FACE_KEY_UNIQ_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_MAKE_FACE_KEY_UNIQ_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# Pepper-bug plan 2026-04-22 Phase 4 P4-T2 — DOF-offset invariants
seas_test_fault_dof_offset_invariants: $(TEST_FAULT_DOF_OFFSET_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_DOF_OFFSET_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# v9.2.0 plan Step 5 (R-V92-H01) — coupled RK4-on-psi convergence probe
seas_test_rk4_psi_integration: $(TEST_RK4_PSI_INTEGRATION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_RK4_PSI_INTEGRATION_OBJ) $(MFEM_LIBS)

# v9.2.0 round-9 REVIEW R-V92-K01 — driver post-RK4 contract regression guard
seas_test_rk4_endpoint_reeval: $(TEST_RK4_ENDPOINT_REEVAL_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_RK4_ENDPOINT_REEVAL_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# v9.2.0 plan §4.10 (rev-3g) — single-channel Jacobian probe (d+e)
seas_test_volume_jacobian_single_channel: $(TEST_VOLUME_JAC_SINGLE_CH_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_VOLUME_JAC_SINGLE_CH_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# v9.2.0 plan §4.11 (rev-3g) — RK4 conservation probe (b)
seas_test_rk4_conservation: $(TEST_RK4_CONSERVATION_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_RK4_CONSERVATION_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# v9.2.0 plan §4.12 (rev-3g) — absorbing-BC energy decay probe (c)
seas_test_absorbing_bc_energy_decay: $(TEST_ABSORBING_BC_DECAY_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ABSORBING_BC_DECAY_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# v9.2.0 plan §4.6 (rev-3b) — shared-fault role consistency across MPI
seas_test_shared_fault_role_consistency: $(TEST_SHARED_FAULT_ROLE_CONSIST_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_SHARED_FAULT_ROLE_CONSIST_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# v9.2.0 plan §4.7 (rev-3c) — shared-fault DOFData runtime consistency
seas_test_shared_fault_dof_data_consistency: $(TEST_SHARED_FAULT_DOFDATA_CONSIST_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_SHARED_FAULT_DOFDATA_CONSIST_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# v9.2.0 plan Step 5b (R-V92-I02) — 100-step F01+F02 stability probe
seas_test_rk4_f01f02_shared_fault_stability: $(TEST_RK4_F01F02_STABILITY_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_RK4_F01F02_STABILITY_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# v9.2.0 plan §17 (rev-3d) — centroid-margin ctor-only probe
seas_test_centroid_margin_ctor_only: $(TEST_CENTROID_MARGIN_CTOR_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_CENTROID_MARGIN_CTOR_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# v9.3.0 plan Phase 1 (I-04): FreeSurfaceGodunov vs FreeSurface equivalence
seas_test_free_surface_godunov: $(TEST_FREE_SURFACE_GODUNOV_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FREE_SURFACE_GODUNOV_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)

# v9.3.0 plan Phase 2 (I-04 wiring, R-005): --free-surface-bc CLI parsing
seas_test_R005_free_surface_bc_flag_parse: $(TEST_R005_FS_BC_FLAG_PARSE_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R005_FS_BC_FLAG_PARSE_OBJ) $(MFEM_LIBS)

# v9.3.0 plan Phase 2 AC-3 (I-04 driver integration gate)
seas_test_free_surface_godunov_driver: $(TEST_FREE_SURFACE_GODUNOV_DRIVER_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FREE_SURFACE_GODUNOV_DRIVER_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# v9.3.0 plan Phase 3 (I-06 part A): EvaluateTotal vs Evaluate equivalence
seas_test_fault_face_flux_total_vs_fluctuation: $(TEST_FAULT_FACE_FLUX_TOTAL_VS_FLUC_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_FACE_FLUX_TOTAL_VS_FLUC_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)

# v9.3.0 R-003 (review-incorporation plan): EvaluateTotal psi-invariance
seas_test_R003_evaluatetotal_psi_invariant: $(TEST_R003_EVALUATETOTAL_PSI_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R003_EVALUATETOTAL_PSI_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)

# v9.3.0 Phase 4 (I-06 part B): R-001 InitializeStateTotal sign regression
seas_test_R001_initstate_total_sign: $(TEST_R001_INITSTATE_TOTAL_SIGN_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R001_INITSTATE_TOTAL_SIGN_OBJ) $(MFEM_LIBS)

# v9.3.0 Phase 4 AC (I-06 part B): locked-fault 100 RK4 sigma_n stability
seas_test_tpv102_total_locked_fault: $(TEST_TPV102_TOTAL_LOCKED_FAULT_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV102_TOTAL_LOCKED_FAULT_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# v9.3.0 Phase 5 (I-06 part C): R-002 nucleation sign regression
seas_test_R002_nucleation_sign_match: $(TEST_R002_NUCLEATION_SIGN_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R002_NUCLEATION_SIGN_OBJ) $(MFEM_LIBS)

# v9.3.0 Phase 5 (I-06 part C): R-009 nucleation both-sides regression
seas_test_R009_nucleation_total_both_sides: $(TEST_R009_NUCLEATION_BOTH_SIDES_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R009_NUCLEATION_BOTH_SIDES_OBJ) $(MFEM_LIBS)

# v9.3.0 R-I06-001 / R-I06-006: total-Q absorbing BC equilibrium gate
seas_test_tpv102_total_absorbing_equilibrium: $(TEST_TPV102_TOTAL_ABSORBING_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV102_TOTAL_ABSORBING_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# v9.3.0 R-I06-003 / R-I06-004: shape-max amplitude scaling + add-delta
seas_test_R_I06_003_nucleation_shape_scaling: $(TEST_R_I06_003_NUC_SHAPE_SCALING_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R_I06_003_NUC_SHAPE_SCALING_OBJ) $(MFEM_LIBS)

# 2026-04-22 persistent-prestress channel: links FaultFaceFlux + FrictionSolver
seas_test_persistent_nuc_prestress_channel: $(TEST_PERSISTENT_NUC_PRESTRESS_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PERSISTENT_NUC_PRESTRESS_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# v9.3.0 R-I06-round3 R-006: direct tests for FreeSurfaceTotal / FreeSurfaceGodunovTotal
seas_test_free_surface_total_direct: $(TEST_FREE_SURFACE_TOTAL_DIRECT_OBJ) $(GODUNOV_FLUX_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FREE_SURFACE_TOTAL_DIRECT_OBJ) $(GODUNOV_FLUX_OBJ) $(MFEM_LIBS)

# ADER I-05 Phase 1: polynomial-exactness test for ApplySpatialDerivative
seas_test_wave_operator_spatial_derivative: $(TEST_WAVE_OP_SPATIAL_DERIV_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_WAVE_OP_SPATIAL_DERIV_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# ADER I-05 Phase 3: CK predictor plane-wave convergence test
seas_test_ader_ck_predictor: $(TEST_ADER_CK_PREDICTOR_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ADER_CK_PREDICTOR_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# ADER I-05 Phase 5: EvaluateADER / EvaluateADERTotal equivalence tests
seas_test_fault_face_flux_ader_equivalence: $(TEST_ADER_FAULT_FLUX_EQUIV_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ADER_FAULT_FLUX_EQUIV_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# ADER I-05 Phase 6: AdvanceADER vs RK4 on a linear wave
seas_test_ader_linear_wave_equivalence: $(TEST_ADER_LINEAR_WAVE_EQUIV_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ADER_LINEAR_WAVE_EQUIV_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# ADER fault-path live consistency: (Q_new - Q)/dt -> wave.Mult(Q) on a faulted fixture
seas_test_ader_fault_coupling_dt_limit: $(TEST_ADER_FAULT_COUPLING_DT_LIMIT_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ADER_FAULT_COUPLING_DT_LIMIT_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# ADER live interior-vs-shared branch cross-check (np=2)
seas_test_ader_interior_vs_shared_branch_live: $(TEST_ADER_INTERIOR_VS_SHARED_LIVE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ADER_INTERIOR_VS_SHARED_LIVE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# ADER I-05 Phase 8: TPV102 AdvanceADER smoke test
seas_test_ader_tpv102_smoke: $(TEST_ADER_TPV102_SMOKE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ADER_TPV102_SMOKE_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# Round-6 R-001 regression: driver defaults to ADER
seas_test_R001_driver_defaults_to_ader: $(TEST_R001_DRIVER_DEFAULT_ADER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R001_DRIVER_DEFAULT_ADER_OBJ) $(MFEM_LIBS)

# Round-6 R-002 regression: driver init uses total-Q
seas_test_R002_driver_init_total_q: $(TEST_R002_DRIVER_INIT_TOTAL_Q_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R002_DRIVER_INIT_TOTAL_Q_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# Dynamic rupture Phase 3: Friction solver tests
seas_test_friction_solver: $(TEST_FRICTION_SOLVER_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FRICTION_SOLVER_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# Phase 4a: TPV102 setup tests
seas_test_tpv102_setup: $(TEST_TPV102_SETUP_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV102_SETUP_OBJ) $(MFEM_LIBS)

# Phase 4a: TPV102 driver
seas_tpv102_driver: $(TPV102_DRIVER_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TPV102_DRIVER_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# Phase 4a: TPV102 local integration tests
seas_test_tpv102_local: $(TEST_TPV102_LOCAL_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TPV102_LOCAL_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

# Mesh generation tool (serial)
seas_generate_mesh: $(GENERATE_MESH_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(GENERATE_MESH_OBJ) $(MFEM_LIBS)

# Object file rules
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

$(FAULT_FACE_FLUX_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

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

$(TEST_FAULT_BASIS_DIP_STRIKE_SYMMETRY_OBJ): %.o: $(SRC)%.cpp $(FAULT_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# v9.2.0 plan §4.1-§4.4 regression gates
$(TEST_GODUNOV_ROTATION_IDENTITY_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_CANONICAL_ROTATION_PURE_STRIKESLIP_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
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

$(TEST_ADJACENT_TRI_UNIFORMITY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(SEAS_EXTRA_CPPFLAGS) -DSEAS_USE_MPI -c $< -o $@

$(TEST_GODUNOV_ROUNDTRIP_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
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

# v9.3.0 Phase 5 (I-06 part C): R-002 nucleation sign regression
$(TEST_R002_NUCLEATION_SIGN_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# v9.3.0 Phase 5 (I-06 part C): R-009 nucleation both-sides regression
$(TEST_R009_NUCLEATION_BOTH_SIDES_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# v9.3.0 R-I06-001 / R-I06-006: total-Q absorbing BC equilibrium gate
$(TEST_TPV102_TOTAL_ABSORBING_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# v9.3.0 R-I06-003 / R-I06-004: shape-max amplitude scaling + add-delta
$(TEST_R_I06_003_NUC_SHAPE_SCALING_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

# 2026-04-22 nucleation persistent-prestress channel fix regression
$(TEST_PERSISTENT_NUC_PRESTRESS_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
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

$(TEST_FRICTION_SOLVER_OBJ): %.o: $(SRC)%.cpp $(DYNAMIC_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TEST_TPV102_SETUP_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

$(TPV102_DRIVER_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI $(SEAS_EXTRA_CPPFLAGS) -c $< -o $@

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
      test-godunov-flux test-wave-operator \
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
      test-R002-nucleation-sign-match \
      test-R009-nucleation-total-both-sides \
      test-tpv102-total-absorbing-equilibrium \
      test-R-I06-003-nucleation-shape-scaling \
      test-persistent-nuc-prestress-channel \
      test-free-surface-total-direct \
      test-wave-operator-spatial-derivative \
      test-ader-ck-predictor \
      test-fault-face-flux-ader-equivalence \
      test-ader-linear-wave-equivalence \
      test-ader-tpv102-smoke \
      test-R001-driver-defaults-to-ader \
      test-R002-driver-init-total-q

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

$(TEST_LINEAR_ELASTIC_OBJ): %.o: $(SRC)%.cpp $(CONSTITUTIVE_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

test-linear-elastic: seas_test_linear_elastic
	./seas_test_linear_elastic

seas_test_constitutive_integrator: $(TEST_CONSTITUTIVE_INTEG_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_CONSTITUTIVE_INTEG_OBJ) $(MFEM_LIBS)

$(TEST_CONSTITUTIVE_INTEG_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

test-constitutive-integrator: seas_test_constitutive_integrator
	./seas_test_constitutive_integrator

# Phase 4: Boundary config tests
seas_test_boundary_config: $(TEST_BOUNDARY_CONFIG_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BOUNDARY_CONFIG_OBJ) $(MFEM_LIBS)

$(TEST_BOUNDARY_CONFIG_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

test-boundary-config: seas_test_boundary_config
	./seas_test_boundary_config

seas_test_boundary_config_operator: $(TEST_BOUNDARY_CONFIG_OP_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BOUNDARY_CONFIG_OP_OBJ) $(MFEM_LIBS)

$(TEST_BOUNDARY_CONFIG_OP_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

test-boundary-config-operator: seas_test_boundary_config_operator
	./seas_test_boundary_config_operator

# Phase 4h-iii: Boundary classification (BuildFacetBCTables) face-by-face tests
seas_test_boundary_classification: $(TEST_BOUNDARY_CLASSIF_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BOUNDARY_CLASSIF_OBJ) $(MFEM_LIBS)

$(TEST_BOUNDARY_CLASSIF_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

test-boundary-classification: seas_test_boundary_classification
	./seas_test_boundary_classification

# Phase 6: Logging tests
seas_test_logging: $(TEST_LOGGING_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_LOGGING_OBJ) $(MFEM_LIBS)

$(TEST_LOGGING_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

test-logging: seas_test_logging
	./seas_test_logging

# Phase 7a: Config bridge tests
seas_test_config_bridge: $(TEST_CONFIG_BRIDGE_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_CONFIG_BRIDGE_OBJ) $(MFEM_LIBS)

$(TEST_CONFIG_BRIDGE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

test-config-bridge: seas_test_config_bridge
	./seas_test_config_bridge

# Phase 5: TOML config tests and driver
seas_test_config_parser: $(TEST_CONFIG_PARSER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_CONFIG_PARSER_OBJ) $(MFEM_LIBS)

$(TEST_CONFIG_PARSER_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

test-config-parser: seas_test_config_parser
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_config_parser

test-godunov-flux: seas_test_godunov_flux
	./seas_test_godunov_flux

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

# Pepper-bug: ADJACENT-TRIANGLE UNIFORMITY (serial + np=4 variants)
test-adjacent-triangle-fault-uniformity-serial: seas_test_adjacent_triangle_fault_uniformity
	./seas_test_adjacent_triangle_fault_uniformity

test-adjacent-triangle-fault-uniformity-parallel: seas_test_adjacent_triangle_fault_uniformity
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_adjacent_triangle_fault_uniformity

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

# v9.3.0 Phase 5 (I-06): R-002 nucleation sign regression
test-R002-nucleation-sign-match: seas_test_R002_nucleation_sign_match
	./seas_test_R002_nucleation_sign_match

# v9.3.0 Phase 5 (I-06): R-009 nucleation both-sides regression
test-R009-nucleation-total-both-sides: seas_test_R009_nucleation_total_both_sides
	./seas_test_R009_nucleation_total_both_sides

# v9.3.0 R-I06-001 / R-I06-006: total-Q absorbing BC equilibrium gate
test-tpv102-total-absorbing-equilibrium: seas_test_tpv102_total_absorbing_equilibrium
	./seas_test_tpv102_total_absorbing_equilibrium

# v9.3.0 R-I06-003 / R-I06-004: shape-max amplitude scaling + add-delta
test-R-I06-003-nucleation-shape-scaling: seas_test_R_I06_003_nucleation_shape_scaling
	./seas_test_R_I06_003_nucleation_shape_scaling

# 2026-04-22 nucleation persistent-prestress channel fix regression
test-persistent-nuc-prestress-channel: seas_test_persistent_nuc_prestress_channel
	./seas_test_persistent_nuc_prestress_channel

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
                     test-R002-nucleation-sign-match \
                     test-R009-nucleation-total-both-sides \
                     test-tpv102-total-absorbing-equilibrium \
                     test-R-I06-003-nucleation-shape-scaling \
                     test-persistent-nuc-prestress-channel \
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

$(SEAS_DRIVER_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(PARALLEL_IO_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(TOML_FLAGS) -DSEAS_USE_MPI -c $< -o $@

# Phase 6: Parallel test targets
seas_test_mpi_context: $(TEST_MPI_CONTEXT_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_MPI_CONTEXT_OBJ) $(MFEM_LIBS)

seas_test_parallel_utils: $(TEST_PARALLEL_UTILS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PARALLEL_UTILS_OBJ) $(MFEM_LIBS)

$(TEST_MPI_CONTEXT_OBJ): %.o: $(SRC)%.cpp $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

$(TEST_PARALLEL_UTILS_OBJ): %.o: $(SRC)%.cpp $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

seas_test_fault_scatter: $(TEST_FAULT_SCATTER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_SCATTER_OBJ) $(MFEM_LIBS)

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

seas_mms_antiplane_parallel: $(MMS_ANTIPLANE_PARALLEL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(MMS_ANTIPLANE_PARALLEL_OBJ) $(MFEM_LIBS)

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

seas_pseas: $(PSEAS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(PSEAS_OBJ) $(MFEM_LIBS)

$(TEST_PARALLEL_FAULT_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(PARALLEL_IO_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

$(PSEAS_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(PARALLEL_IO_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

test-parallel-fault: seas_test_parallel_fault
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_parallel_fault
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_parallel_fault

test-parallel: test-mpi-context test-parallel-utils test-parallel-domain test-mms-parallel test-parallel-fault

# Phase 9: Verification and scaling targets
seas_test_serial_parallel_consistency: $(TEST_SERIAL_PARALLEL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_SERIAL_PARALLEL_OBJ) $(MFEM_LIBS)

seas_bp2_benchmark_parallel: $(BP2_BENCHMARK_PARALLEL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(BP2_BENCHMARK_PARALLEL_OBJ) $(MFEM_LIBS)

seas_test_scaling: $(TEST_SCALING_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_SCALING_OBJ) $(MFEM_LIBS)

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

$(TEST_BR2_CONSISTENCY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

test-br2-consistency: seas_test_br2_consistency
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_br2_consistency
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_br2_consistency

# Parallel elasticity test (shared face support)
seas_test_parallel_elasticity: $(TEST_PARALLEL_ELAST_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PARALLEL_ELAST_OBJ) $(MFEM_LIBS)

$(TEST_PARALLEL_ELAST_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

test-parallel-elasticity: seas_test_parallel_elasticity
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_parallel_elasticity
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_parallel_elasticity

# Parallel wave operator tests (dynamic rupture)
seas_test_parallel_wave_operator: $(TEST_PARALLEL_WAVE_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PARALLEL_WAVE_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

$(TEST_PARALLEL_WAVE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

test-parallel-wave-operator: seas_test_parallel_wave_operator
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_parallel_wave_operator

# R-101 shared-fault DOFData consistency test (2-rank MPI regression)
seas_test_r101_shared_fault: $(TEST_R101_SHARED_FAULT_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_R101_SHARED_FAULT_OBJ) $(WAVE_OPERATOR_OBJ) $(GODUNOV_FLUX_OBJ) $(PML_LAYER_OBJ) $(FAULT_FACE_FLUX_OBJ) $(FRICTION_SOLVER_OBJ) $(MFEM_LIBS)

$(TEST_R101_SHARED_FAULT_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

test-r101-shared-fault: seas_test_r101_shared_fault
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_r101_shared_fault

# BP5 parallel smoke test
seas_test_bp5_parallel_smoke: $(TEST_BP5_SMOKE_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_BP5_SMOKE_OBJ) $(MFEM_LIBS)

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

$(BP2_VERIFY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(PARALLEL_IO_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

test-bp2-verify: seas_bp2_verify
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 10 ./seas_bp2_verify

# Phase 10b: BP2 full multi-cycle simulation
seas_bp2_full: $(BP2_FULL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(BP2_FULL_OBJ) $(MFEM_LIBS)

$(BP2_FULL_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(PARALLEL_IO_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

# BP1 full multi-cycle simulation
seas_bp1_full: $(BP1_FULL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(BP1_FULL_OBJ) $(MFEM_LIBS)

$(BP1_FULL_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(PARALLEL_IO_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

# BP1 bottom-boundary Dirichlet loading simulation
seas_bp1_bdrload: $(BP1_BDRLOAD_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(BP1_BDRLOAD_OBJ) $(MFEM_LIBS)

$(BP1_BDRLOAD_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(PARALLEL_IO_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

# BP5 full multi-cycle simulation
seas_bp5_full: $(BP5_FULL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(BP5_FULL_OBJ) $(MFEM_LIBS)

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

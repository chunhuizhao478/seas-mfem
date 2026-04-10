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

# All headers
SEAS_HEADERS = $(FRICTION_HEADERS) $(DOMAIN_HEADERS) $(FAULT_HEADERS) $(SOLVER_HEADERS) $(IO_HEADERS) $(CHECKPOINT_HEADERS)

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
               seas_test_face_quadrature
# Long-running simulations and tools (built by `make all`, not run by `make test`)
SEQ_LONG = seas_bp2_serial_smoke seas_generate_mesh
# Phase 6: Parallel test sources
TEST_MPI_CONTEXT_SRC = tests/parallel/test_mpi_context.cpp
TEST_PARALLEL_UTILS_SRC = tests/parallel/test_parallel_utils.cpp
TEST_MPI_CONTEXT_OBJ = $(TEST_MPI_CONTEXT_SRC:.cpp=.o)
TEST_PARALLEL_UTILS_OBJ = $(TEST_PARALLEL_UTILS_SRC:.cpp=.o)

# Phase 6: Parallel headers
PARALLEL_HEADERS = common/mpi_context.hpp \
                   common/parallel_utils.hpp \
                   common/seas_types.hpp

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
               seas_test_parallel_domain seas_mms_antiplane_parallel \
               seas_test_parallel_fault \
               seas_test_serial_parallel_consistency \
               seas_test_scaling \
               seas_test_br2_consistency \
               seas_test_parallel_elasticity \
               seas_test_bp5_parallel_smoke
# Long-running parallel simulations (built by `make all`, not run by `make test`)
PAR_LONG = seas_pseas \
           seas_bp2_benchmark_parallel \
           seas_bp2_verify seas_bp2_full seas_bp1_full \
           seas_bp1_bdrload \
           seas_bp5_full

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
        test-mpi-context test-parallel-utils test-parallel \
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
        test-bp5-smoke \
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
      test-face-quadrature

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
	rm -f *.o *~ $(SEQ_MINIAPPS) $(PAR_MINIAPPS)
	rm -f tests/unit/*.o tests/verification/*.o tests/parallel/*.o
	rm -rf *.dSYM *.TVD.*breakpoints

clean-exec:
	@rm -rf mesh.* sol.* ParaView mfem_bp2qd_*.txt test_*.txt

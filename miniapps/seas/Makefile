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

# Object files
TEST_FRICTION_OBJ = $(TEST_FRICTION_SRC:.cpp=.o)
TEST_STATE_EVOLUTION_OBJ = $(TEST_STATE_EVOLUTION_SRC:.cpp=.o)
TEST_ANTIPLANE_OBJ = $(TEST_ANTIPLANE_SRC:.cpp=.o)
TEST_FAULT_OPERATOR_OBJ = $(TEST_FAULT_OPERATOR_SRC:.cpp=.o)
TEST_QUASI_DYNAMIC_OBJ = $(TEST_QUASI_DYNAMIC_SRC:.cpp=.o)
TEST_BP2_SHORT_OBJ = $(TEST_BP2_SHORT_SRC:.cpp=.o)
TEST_IO_OBJ = $(TEST_IO_SRC:.cpp=.o)
BP2_SERIAL_SMOKE_OBJ = $(BP2_SERIAL_SMOKE_SRC:.cpp=.o)

# Header dependencies
# Phase 1: Friction and state evolution
FRICTION_HEADERS = friction/friction_law.hpp \
                   friction/state_evolution.hpp \
                   friction/dieterich_ruina.hpp \
                   config/bp2_params.hpp

# Phase 2: Domain operator
DOMAIN_HEADERS = common/seas_types.hpp \
                 common/mpi_context.hpp \
                 domain/domain_operator.hpp \
                 domain/antiplane_operator.hpp \
                 domain/bp2_mesh.hpp

# Phase 3: Fault operator
FAULT_HEADERS = fault/fault_geometry.hpp \
                fault/rate_state_fault.hpp

# Phase 4: SEAS operator
SOLVER_HEADERS = solver/seas_operator.hpp \
                 solver/time_stepper.hpp

# Phase 5: I/O
IO_HEADERS = io/probe_output.hpp \
             io/benchmark_output.hpp \
             io/paraview_output.hpp

# All headers
SEAS_HEADERS = $(FRICTION_HEADERS) $(DOMAIN_HEADERS) $(FAULT_HEADERS) $(SOLVER_HEADERS) $(IO_HEADERS)

# Targets
SEQ_MINIAPPS = seas_test_friction_law seas_test_state_evolution seas_test_antiplane \
               seas_test_fault_operator seas_test_quasi_dynamic seas_test_bp2_short \
               seas_test_io seas_bp2_serial_smoke
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

PAR_MINIAPPS = seas_test_mpi_context seas_test_parallel_utils \
               seas_test_parallel_domain seas_mms_antiplane_parallel

ifeq ($(MFEM_USE_MPI),NO)
   MINIAPPS = $(SEQ_MINIAPPS)
else
   MINIAPPS = $(PAR_MINIAPPS) $(SEQ_MINIAPPS)
endif

MFEM_TESTS = MINIAPPS

.SUFFIXES:
.SUFFIXES: .o .cpp .mk
.PHONY: all clean clean-build clean-exec test \
        test-friction test-state-evolution test-antiplane \
        test-fault-operator test-quasi-dynamic test-bp2-short \
        test-io test-bp2-serial-smoke \
        test-mpi-context test-parallel-utils test-parallel \
        test-parallel-domain test-mms-parallel

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

# Testing targets
test: test-friction test-state-evolution test-antiplane test-fault-operator \
      test-quasi-dynamic test-bp2-short test-io

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

test-parallel: test-mpi-context test-parallel-utils test-parallel-domain test-mms-parallel

-include $(MFEM_TEST_MK)

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

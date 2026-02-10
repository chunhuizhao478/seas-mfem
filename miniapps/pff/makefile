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

# Phase Field Fracture (PFF) miniapp makefile

# Use the MFEM build directory
MFEM_DIR ?= ../..
MFEM_BUILD_DIR ?= ../..
MFEM_INSTALL_DIR ?= ../../mfem
SRC = $(if $(MFEM_DIR:../..=),$(MFEM_DIR)/miniapps/pff/,)
CONFIG_MK = $(or $(wildcard $(MFEM_BUILD_DIR)/config/config.mk),\
   $(wildcard $(MFEM_INSTALL_DIR)/share/mfem/config.mk))

MFEM_LIB_FILE = mfem_is_not_built
-include $(CONFIG_MK)

# Source files
UNIT_TESTS_PHASE1_SRC = tests/unit_tests_phase1.cpp
UNIT_TESTS_PHASE2_SRC = tests/unit_tests_phase2.cpp
UNIT_TESTS_PHASE3_SRC = tests/unit_tests_phase3.cpp
UNIT_TESTS_PHASE4_SRC = tests/unit_tests_phase4.cpp

# Object files
UNIT_TESTS_PHASE1_OBJ = $(UNIT_TESTS_PHASE1_SRC:.cpp=.o)
UNIT_TESTS_PHASE2_OBJ = $(UNIT_TESTS_PHASE2_SRC:.cpp=.o)
UNIT_TESTS_PHASE3_OBJ = $(UNIT_TESTS_PHASE3_SRC:.cpp=.o)
UNIT_TESTS_PHASE4_OBJ = $(UNIT_TESTS_PHASE4_SRC:.cpp=.o)

# Header dependencies
PFF_HEADERS = materials/degradation_function.hpp \
              materials/spectral_decomposition.hpp \
              materials/pff_material.hpp \
              integrators/degraded_elasticity_integrator.hpp \
              integrators/damage_diffusion_integrator.hpp \
              integrators/damage_source_integrator.hpp \
              operators/elasticity_operator.hpp \
              operators/damage_operator.hpp \
              pff_solver.hpp

# Targets
SEQ_MINIAPPS = pff_unit_tests_phase1 pff_unit_tests_phase2 pff_unit_tests_phase3 pff_unit_tests_phase4 pff
PAR_MINIAPPS =

ifeq ($(MFEM_USE_MPI),NO)
   MINIAPPS = $(SEQ_MINIAPPS)
else
   MINIAPPS = $(PAR_MINIAPPS) $(SEQ_MINIAPPS)
endif

MFEM_TESTS = MINIAPPS

.SUFFIXES:
.SUFFIXES: .o .cpp .mk
.PHONY: all clean clean-build clean-exec test test-phase1 test-phase2 test-phase3 test-phase4

all: $(MINIAPPS)

# Remove built-in rules
%: %.cpp
%.o: %.cpp

PFF_INCLUDES = $(if $(SRC),-I$(MFEM_DIR))

# Phase 1 unit tests
pff_unit_tests_phase1: $(UNIT_TESTS_PHASE1_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(UNIT_TESTS_PHASE1_OBJ) $(MFEM_LIBS)

# Phase 2 unit tests
pff_unit_tests_phase2: $(UNIT_TESTS_PHASE2_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(UNIT_TESTS_PHASE2_OBJ) $(MFEM_LIBS)

# Phase 3 unit tests
pff_unit_tests_phase3: $(UNIT_TESTS_PHASE3_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(UNIT_TESTS_PHASE3_OBJ) $(MFEM_LIBS)

# Phase 4 unit tests
pff_unit_tests_phase4: $(UNIT_TESTS_PHASE4_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(UNIT_TESTS_PHASE4_OBJ) $(MFEM_LIBS)

# Phase 5 main driver
pff: pff.o
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ pff.o $(MFEM_LIBS)

# Object file rules
$(UNIT_TESTS_PHASE1_OBJ): %.o: $(SRC)%.cpp $(PFF_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(PFF_INCLUDES) -c $< -o $@

$(UNIT_TESTS_PHASE2_OBJ): %.o: $(SRC)%.cpp $(PFF_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(PFF_INCLUDES) -c $< -o $@

$(UNIT_TESTS_PHASE3_OBJ): %.o: $(SRC)%.cpp $(PFF_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(PFF_INCLUDES) -c $< -o $@

$(UNIT_TESTS_PHASE4_OBJ): %.o: $(SRC)%.cpp $(PFF_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(PFF_INCLUDES) -c $< -o $@

pff.o: $(SRC)pff.cpp $(PFF_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(PFF_INCLUDES) -c $< -o $@

# Testing targets
test: test-phase1 test-phase2 test-phase3 test-phase4

test-phase1: pff_unit_tests_phase1
	./pff_unit_tests_phase1

test-phase2: pff_unit_tests_phase2
	./pff_unit_tests_phase2

test-phase3: pff_unit_tests_phase3
	./pff_unit_tests_phase3

test-phase4: pff_unit_tests_phase4
	./pff_unit_tests_phase4

-include $(MFEM_TEST_MK)

# Generate an error message if the MFEM library is not built and exit
$(MFEM_LIB_FILE):
	$(error The MFEM library is not built)

clean: clean-build clean-exec

clean-build:
	rm -f *.o *~ $(SEQ_MINIAPPS) $(PAR_MINIAPPS)
	rm -f tests/*.o
	rm -rf *.dSYM *.TVD.*breakpoints

clean-exec:
	@rm -rf mesh.* sol.* ParaView

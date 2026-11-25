# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# Dynamic phase-field fracture (pff) miniapp makefile, laid out like
# miniapps/seas/Makefile: it reads the MFEM build's config.mk.
#
#   make test                       build and run the unit tests
#   make MFEM_BUILD_DIR=/path test  against an out-of-tree MFEM build (CMake or make)
#
# Targets are added phase by phase (see PLAN.md).

MFEM_DIR ?= ../..
MFEM_BUILD_DIR ?= ../..
MFEM_INSTALL_DIR ?= ../../mfem
CONFIG_MK = $(or $(wildcard $(MFEM_BUILD_DIR)/config/config.mk),\
   $(wildcard $(MFEM_INSTALL_DIR)/share/mfem/config.mk))

MFEM_LIB_FILE = mfem_is_not_built
-include $(CONFIG_MK)

PFF_INC = -I. -Iextern/toml11
CXXFLAGS_PFF = $(MFEM_FLAGS) $(PFF_INC)
LIBS_PFF = $(MFEM_LIBS)

# -----------------------------------------------------------------------------
# Unit tests (serial, no mesh): tests/unit/test_*.cpp, one executable each
# -----------------------------------------------------------------------------
UNIT_TESTS = test_phase_field_functions test_strain_split test_poro_damage \
             test_pulse_load test_config test_explicit_operator test_damage_step
UNIT_BINS = $(addprefix tests/unit/,$(UNIT_TESTS))

HEADERS = $(wildcard common/*.hpp config/*.hpp constitutive/*.hpp phasefield/*.hpp \
                     loading/*.hpp io/*.hpp multiphysics/poroelastic/*.hpp \
                     domain/*.hpp integrator/*.hpp kernels/*.hpp solver/*.hpp tests/*.hpp)

# Drivers: drivers/*.cpp, one executable each in the miniapp directory
DRIVERS = qs_driver pff_driver pff_partition

.PHONY: all test clean

all: $(UNIT_BINS) $(DRIVERS)

qs_driver: drivers/qs_driver.cpp $(HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	$(MFEM_CXX) $(CXXFLAGS_PFF) $< -o $@ $(LIBS_PFF)

pff_driver: drivers/pff_driver.cpp $(HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	$(MFEM_CXX) $(CXXFLAGS_PFF) $< -o $@ $(LIBS_PFF)

pff_partition: drivers/pff_partition.cpp $(HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	$(MFEM_CXX) $(CXXFLAGS_PFF) $< -o $@ $(LIBS_PFF)

# Verification programs (tests/verification/*.cpp)
VERIFY_PROGS = tests/verification/dashpot_reflection tests/verification/cavity_hoop

tests/verification/%: tests/verification/%.cpp $(HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	$(MFEM_CXX) $(CXXFLAGS_PFF) $< -o $@ $(LIBS_PFF)

tests/unit/%: tests/unit/%.cpp $(HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	$(MFEM_CXX) $(CXXFLAGS_PFF) $< -o $@ $(LIBS_PFF)

test: $(UNIT_BINS)
	@fail=0; for t in $(UNIT_BINS); do \
	   echo "=== $$t"; ./$$t > $$t.log 2>&1 && tail -1 $$t.log || { cat $$t.log; fail=1; }; \
	done; exit $$fail

clean:
	rm -f $(UNIT_BINS) $(DRIVERS) $(VERIFY_PROGS) tests/unit/*.log

# Make sure the MFEM library is built
$(MFEM_LIB_FILE):
	@echo "MFEM is not built: build it first, or set MFEM_BUILD_DIR"; exit 1

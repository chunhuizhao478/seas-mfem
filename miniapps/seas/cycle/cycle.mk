# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# cycle/cycle.mk — build file of the fully dynamic BP5 cycle branch
# (plan rev 3, rule M2: zero edits to the shared Makefile).
#
# Run from miniapps/seas:
#   make -f cycle/cycle.mk seas_fd_cycle_driver   # the driver
#   make -f cycle/cycle.mk test-cycle-api-pin     # rule M5, compile-time pins
#   make -f cycle/cycle.mk test-cycle-params-parity  # P1, deck vs BP5Params
#   make -f cycle/cycle.mk cycle-owned-check      # rule M1, owned paths only
#   make -f cycle/cycle.mk cycle-all              # all builds and checks (default)
# Without the conda HDF5 (CONDA_PREFIX), point HDF5_INCFLAGS and HDF5_LIBS at
# the system one, e.g. on Ubuntu with Open MPI:
#   make -f cycle/cycle.mk HDF5_INCFLAGS=-I/usr/include/hdf5/openmpi \
#        HDF5_LIBS="-L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -lhdf5"
# Every shared target (seas_bp5_full, regression-check, ...) works through
# this file too, because it includes the shared Makefile unchanged.
#
# Rule M10: the driver prints the commit it was built from.  CYCLE_COMMIT is
# `git rev-parse HEAD`, with "-dirty" when a tracked file differs from it or
# an owned path holds an untracked file (a new source not committed yet).
# cycle/.commit.cache (ignored by the root .gitignore pattern *.cache*) holds
# the last value and is rewritten only when it changes, so every object that
# carries the hash is rebuilt after a commit or checkout.

include Makefile

.DEFAULT_GOAL := cycle-all

CYCLE_OWNED_PATHS = cycle drivers/fd_cycle_driver.cpp tests/cycle config/fd_cycle \
                    jobs/fd_cycle document/fd_cycle_dev
CYCLE_GIT_HEAD  := $(shell git rev-parse HEAD 2>/dev/null)
CYCLE_GIT_DIRTY := $(shell { git status --porcelain --untracked-files=no; \
                             git status --porcelain --untracked-files=all -- $(CYCLE_OWNED_PATHS); \
                           } 2>/dev/null | head -n 1)
CYCLE_COMMIT    := $(if $(CYCLE_GIT_HEAD),$(CYCLE_GIT_HEAD)$(if $(CYCLE_GIT_DIRTY),-dirty),unknown)
CYCLE_STAMP     := cycle/.commit.cache
$(shell [ "$$(cat $(CYCLE_STAMP) 2>/dev/null)" = "$(CYCLE_COMMIT)" ] || \
        printf '%s\n' "$(CYCLE_COMMIT)" > $(CYCLE_STAMP))
CYCLE_COMMIT_FLAG = -DSEAS_CYCLE_COMMIT='"$(CYCLE_COMMIT)"'

CYCLE_HEADERS = $(wildcard cycle/*.hpp)
CYCLE_CXXFLAGS = $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) $(TOML_FLAGS) \
                 -DSEAS_USE_MPI $(CYCLE_COMMIT_FLAG)
# Shared headers the cycle includes (QD stack, deck loader, material field).
CYCLE_SHARED_HEADERS = $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(PARALLEL_IO_HEADERS) \
                       spatial/code/spatial_friction.hpp dynamic/heterogeneous_material.hpp
# Shared objects: the deck loader and resolver (spatial_friction) and what it
# links (material field, coefficients, HDF5 data fields).
CYCLE_LINK_OBJS = $(SPATIAL_FRICTION_OBJ) $(HETEROGENEOUS_MATERIAL_OBJ) \
                  $(MATERIAL_COEFFICIENTS_OBJ) $(DATA_FIELD_3D_OBJ)
CYCLE_LIBS = $(MFEM_LIBS) $(HDF5_LIBS) $(SEAS_FS_LIB)

# The driver (plan item j).
FD_CYCLE_DRIVER_SRC = drivers/fd_cycle_driver.cpp
FD_CYCLE_DRIVER_OBJ = $(FD_CYCLE_DRIVER_SRC:.cpp=.o)
$(FD_CYCLE_DRIVER_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                        $(CYCLE_SHARED_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_fd_cycle_driver: $(FD_CYCLE_DRIVER_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(FD_CYCLE_DRIVER_OBJ) $(CYCLE_LINK_OBJS) \
	    $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# P1 test: deck parameters against BP5Params, kernel psi0 against
# BuildPrestress (no domain solve; runs on the mesh of the decks).
CYCLE_PARITY_SRC = tests/cycle/test_cycle_params_parity.cpp
CYCLE_PARITY_OBJ = $(CYCLE_PARITY_SRC:.cpp=.o)
$(CYCLE_PARITY_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                     $(CYCLE_SHARED_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_cycle_params_parity: $(CYCLE_PARITY_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(CYCLE_PARITY_OBJ) $(CYCLE_LINK_OBJS) \
	    $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

CYCLE_TEST_NP ?= 1
test-cycle-params-parity: seas_test_cycle_params_parity
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $(CYCLE_TEST_NP) ./seas_test_cycle_params_parity

# Rule M5: the API pin test.  Compile-time only; it includes both solver
# stacks in one translation unit, as the driver will.
CYCLE_API_PIN_SRC = tests/cycle/test_cycle_api_pin.cpp
CYCLE_API_PIN_OBJ = $(CYCLE_API_PIN_SRC:.cpp=.o)
$(CYCLE_API_PIN_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(SEAS_HEADERS) \
                      $(PARALLEL_HEADERS) $(PARALLEL_IO_HEADERS) \
                      domain/nitsche_wave_operator.hpp solver/nitsche_cfl.hpp \
                      friction/rate_state_closure.hpp spatial/code/spatial_friction.hpp \
                      $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_cycle_api_pin: $(CYCLE_API_PIN_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(CYCLE_API_PIN_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-cycle-api-pin: seas_test_cycle_api_pin
	./seas_test_cycle_api_pin

# Rule M1: the branch changes only its own paths.  CYCLE_BASE is the shared
# branch this one merges from (M7); the check diffs HEAD against their merge
# base, so a merge of that branch adds nothing to the list.
CYCLE_BASE ?= origin/nitsche-dynamic-rupture
CYCLE_OWNED_RE = ^miniapps/seas/(cycle/|drivers/fd_cycle_driver\.cpp$$|tests/cycle/|config/fd_cycle/|jobs/fd_cycle/|document/fd_cycle_dev/)
cycle-owned-check:
	@base=$$(git merge-base HEAD $(CYCLE_BASE)) || { echo "ERROR: no merge base with $(CYCLE_BASE)"; exit 1; }; \
	echo "Branch point: $$base ($(CYCLE_BASE))"; \
	git diff --stat $$base HEAD; \
	bad=$$(git diff --name-only $$base HEAD | grep -Ev '$(CYCLE_OWNED_RE)'); \
	if [ -n "$$bad" ]; then echo "FAILED: files outside the owned paths:"; echo "$$bad"; exit 1; fi; \
	echo "PASSED: every changed file is under an owned path (rule M1)"

cycle-all: seas_fd_cycle_driver seas_test_cycle_params_parity test-cycle-api-pin \
           cycle-owned-check

cycle-clean:
	rm -f $(FD_CYCLE_DRIVER_OBJ) $(CYCLE_API_PIN_OBJ) $(CYCLE_PARITY_OBJ) \
	      seas_fd_cycle_driver seas_test_cycle_api_pin seas_test_cycle_params_parity \
	      $(CYCLE_STAMP)

.PHONY: cycle-all cycle-clean cycle-owned-check test-cycle-api-pin \
        test-cycle-params-parity

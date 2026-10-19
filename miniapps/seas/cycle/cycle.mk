# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# cycle/cycle.mk — build file of the fully dynamic BP5 cycle branch
# (plan rev 3, rule M2: zero edits to the shared Makefile).
#
# Run from miniapps/seas:
#   make -f cycle/cycle.mk seas_fd_cycle_driver   # the driver
#   make -f cycle/cycle.mk test-cycle-api-pin     # rule M5, compile-time pins
#   make -f cycle/cycle.mk test-cycle-params-parity  # P1, deck vs BP5Params
#   make -f cycle/cycle.mk test-cycle-p2          # P2, regime hand-over tests
#   make -f cycle/cycle.mk test-cycle-p3          # P3, the switch back to QD
#   make -f cycle/cycle.mk test-cycle-p4          # P4, QD -> dynamic at t > 0
#   make -f cycle/cycle.mk test-cycle-snapshots   # the snapshot files, 7 runs
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
# Shared headers the cycle includes (QD stack, deck loader, material field,
# and the Nitsche dynamic stack that dyn_phase.hpp includes).
CYCLE_DYN_HEADERS = domain/nitsche_wave_operator.hpp domain/nitsche_material_sample.hpp \
                    domain/nitsche_mass.hpp domain/nitsche_energy.hpp solver/nitsche_cfl.hpp \
                    friction/rate_state_closure.hpp friction/fault_closure.hpp
CYCLE_SHARED_HEADERS = $(SEAS_HEADERS) $(PARALLEL_HEADERS) $(PARALLEL_IO_HEADERS) \
                       spatial/code/spatial_friction.hpp dynamic/heterogeneous_material.hpp \
                       $(CYCLE_DYN_HEADERS)
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
$(CYCLE_API_PIN_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(CYCLE_SHARED_HEADERS) \
                      $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_cycle_api_pin: $(CYCLE_API_PIN_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(CYCLE_API_PIN_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-cycle-api-pin: seas_test_cycle_api_pin
	./seas_test_cycle_api_pin

# P2 tests: the frame map and its ghost copies on 1, 2, 4 and 8 ranks, on
# the METIS partition and on one that puts the two sides of the fault on
# different ranks (every fault face shared; even rank counts); then, on the
# first rank count above 1, same K, prestress sign, load clock and the static
# hand-over.  BP5 mesh at 1000 m; the energy balance is a separate, longer
# run (seas_test_cycle_p2 --energy 2).  Hosts with fewer cores need
# CYCLE_P2_MPIFLAGS=--oversubscribe for 8 ranks.
CYCLE_P2_SRC = tests/cycle/test_cycle_p2.cpp
CYCLE_P2_OBJ = $(CYCLE_P2_SRC:.cpp=.o)
$(CYCLE_P2_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                 $(CYCLE_SHARED_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_cycle_p2: $(CYCLE_P2_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(CYCLE_P2_OBJ) $(CYCLE_LINK_OBJS) \
	    $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

CYCLE_P2_NPS ?= 1 2 4 8
CYCLE_P2_MPIFLAGS ?=
test-cycle-p2: seas_test_cycle_p2
	@for n in $(CYCLE_P2_NPS); do \
	   $(MFEM_MPIEXEC) $(CYCLE_P2_MPIFLAGS) $(MFEM_MPIEXEC_NP) $$n ./seas_test_cycle_p2 --frames-only || exit 1; \
	   if [ $$((n % 2)) -eq 0 ]; then \
	      $(MFEM_MPIEXEC) $(CYCLE_P2_MPIFLAGS) $(MFEM_MPIEXEC_NP) $$n ./seas_test_cycle_p2 --frames-only --split-y || exit 1; \
	   fi; \
	done
	$(MFEM_MPIEXEC) $(CYCLE_P2_MPIFLAGS) $(MFEM_MPIEXEC_NP) $(firstword $(filter-out 1,$(CYCLE_P2_NPS)) 1) ./seas_test_cycle_p2

# P3 tests: the switch back to QD (bridge round trip; Reenter against a
# fresh QD run and against a restart from the switch checkpoint).  Needs
# PETSc; three QD factorizations one after the other.
CYCLE_P3_SRC = tests/cycle/test_cycle_p3.cpp
CYCLE_P3_OBJ = $(CYCLE_P3_SRC:.cpp=.o)
$(CYCLE_P3_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                 $(CYCLE_SHARED_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_cycle_p3: $(CYCLE_P3_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(CYCLE_P3_OBJ) $(CYCLE_LINK_OBJS) \
	    $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# 1 rank: bit for bit; 2 ranks: to MUMPS's parallel round-off (see the
# test's file comment).
CYCLE_P3_NPS ?= 1 2
CYCLE_P3_DECK ?= config/fd_cycle/bp5_fd_first_event.toml
CYCLE_P3_WORK ?= $(or $(TMPDIR),/tmp)/seas_test_cycle_p3
test-cycle-p3: seas_test_cycle_p3
	@for n in $(CYCLE_P3_NPS); do \
	   $(MFEM_MPIEXEC) $(CYCLE_P2_MPIFLAGS) $(MFEM_MPIEXEC_NP) $$n ./seas_test_cycle_p3 \
	      --deck $(CYCLE_P3_DECK) --work-dir $(CYCLE_P3_WORK)/np$$n || exit 1; \
	done

# P4 tests: the QD -> dynamic switch at t > 0 (the hand-over closure; one
# DynPhase re-entered).  Needs PETSc; one QD factorization.
CYCLE_P4_SRC = tests/cycle/test_cycle_p4.cpp
CYCLE_P4_OBJ = $(CYCLE_P4_SRC:.cpp=.o)
$(CYCLE_P4_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                 $(CYCLE_SHARED_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_cycle_p4: $(CYCLE_P4_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(CYCLE_P4_OBJ) $(CYCLE_LINK_OBJS) \
	    $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

CYCLE_P4_NPS ?= 1 2
CYCLE_P4_DECK ?= config/fd_cycle/bp5_fd_first_event.toml
test-cycle-p4: seas_test_cycle_p4
	@for n in $(CYCLE_P4_NPS); do \
	   $(MFEM_MPIEXEC) $(CYCLE_P2_MPIFLAGS) $(MFEM_MPIEXEC_NP) $$n ./seas_test_cycle_p4 \
	      --deck $(CYCLE_P4_DECK) || exit 1; \
	done

# The fault and surface snapshot files (tests/cycle/test_cycle_snapshots.sh):
# seven short runs on CYCLE_SNAP_NP ranks in CYCLE_SNAP_WORK (replaced).
CYCLE_SNAP_NP ?= 2
CYCLE_SNAP_WORK ?= ../../../runs_test_cycle_snapshots
test-cycle-snapshots: seas_fd_cycle_driver
	MPIRUN="$(MFEM_MPIEXEC) $(CYCLE_P2_MPIFLAGS)" \
	   tests/cycle/test_cycle_snapshots.sh $(CYCLE_SNAP_WORK) $(CYCLE_SNAP_NP)

# Rule M1: the branch changes only its own paths, plus the shared files that
# hold its hooks (document/fd_cycle_dev/HOOKS.md), and those by additions
# only.  CYCLE_BASE is the shared branch this one merges from (M7); the check
# diffs HEAD against their merge base, so a merge of that branch adds nothing
# to the list.
CYCLE_BASE ?= origin/nitsche-dynamic-rupture
CYCLE_OWNED_RE = ^miniapps/seas/(cycle/|drivers/fd_cycle_driver\.cpp$$|tests/cycle/|config/fd_cycle/|jobs/fd_cycle/|document/fd_cycle_dev/)
CYCLE_HOOK_FILES = miniapps/seas/fault/fault_geometry.hpp \
                   miniapps/seas/domain/nitsche_wave_operator.hpp
CYCLE_HOOK_GREP = grep -xF $(foreach f,$(CYCLE_HOOK_FILES),-e $(f))
# Fixes taken whole from other branches (HOOKS.md, "Shared changes merged from
# other tracks"), cherry-picked with -x: the branch may change the files such
# a commit touches only by exactly that commit's change (same patch id).
CYCLE_IMPORTS = ef70517101ab3824e0955d7c7e3f9013e8eacf8c
cycle-owned-check:
	@base=$$(git merge-base HEAD $(CYCLE_BASE)) || { echo "ERROR: no merge base with $(CYCLE_BASE)"; exit 1; }; \
	echo "Branch point: $$base ($(CYCLE_BASE))"; \
	git diff --stat $$base HEAD; \
	imp=$$(for c in $(CYCLE_IMPORTS); do git diff-tree --no-commit-id --name-only -r $$c; done); \
	bad=$$(git diff --name-only $$base HEAD | grep -Ev '$(CYCLE_OWNED_RE)' | $(CYCLE_HOOK_GREP) -v | \
	       IMP="$$imp" awk 'BEGIN { n = split(ENVIRON["IMP"], f, "\n"); for (i = 1; i <= n; i++) skip[f[i]] = 1 } !($$0 in skip)') \
	   || { echo "ERROR: cycle-owned-check could not list the changed files"; exit 1; }; \
	if [ -n "$$bad" ]; then echo "FAILED: files outside the owned paths, hook files and imports:"; echo "$$bad"; exit 1; fi; \
	del=$$(git diff --numstat $$base HEAD | awk '$$2 != 0 {print $$3}' | $(CYCLE_HOOK_GREP)); \
	if [ -n "$$del" ]; then echo "FAILED: hook files must change by additions only:"; echo "$$del"; exit 1; fi; \
	for c in $(CYCLE_IMPORTS); do \
	   f=$$(git diff-tree --no-commit-id --name-only -r $$c | sed 's/^/:(top)/'); \
	   a=$$(git show $$c -- $$f | git patch-id --stable | cut -d' ' -f1); \
	   b=$$(git diff $$base HEAD -- $$f | git patch-id --stable | cut -d' ' -f1); \
	   if [ -z "$$a" ] || [ "$$a" != "$$b" ]; then echo "FAILED: the branch changes the files of import $$c differently from $$c"; exit 1; fi; \
	   echo "Import $$c: the same change to its $$(echo "$$f" | wc -l) files"; \
	done; \
	echo "PASSED: every changed file is under an owned path, a hook file changed by additions only, or an import (rule M1)"

cycle-all: seas_fd_cycle_driver seas_test_cycle_params_parity seas_test_cycle_p2 \
           seas_test_cycle_p3 seas_test_cycle_p4 test-cycle-api-pin cycle-owned-check

cycle-clean:
	rm -f $(FD_CYCLE_DRIVER_OBJ) $(CYCLE_API_PIN_OBJ) $(CYCLE_PARITY_OBJ) $(CYCLE_P2_OBJ) \
	      $(CYCLE_P3_OBJ) $(CYCLE_P4_OBJ) seas_fd_cycle_driver seas_test_cycle_api_pin \
	      seas_test_cycle_params_parity seas_test_cycle_p2 seas_test_cycle_p3 \
	      seas_test_cycle_p4 $(CYCLE_STAMP)

.PHONY: cycle-all cycle-clean cycle-owned-check test-cycle-api-pin \
        test-cycle-params-parity test-cycle-p2 test-cycle-p3 test-cycle-p4 \
        test-cycle-snapshots

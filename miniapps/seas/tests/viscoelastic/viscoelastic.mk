# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# tests/viscoelastic/viscoelastic.mk — build file of the viscoelastic
# attenuation branch (attenuation plan rev6, section 5.11).  It includes
# tests/2d_inplane/2d_inplane.mk (which includes cycle/cycle.mk and the
# shared Makefile unchanged), so every earlier target works through it too.
#
# Run from miniapps/seas:
#   make -f tests/viscoelastic/viscoelastic.mk test-rheology          # V0: ported pointwise/config tests
#   make -f tests/viscoelastic/viscoelastic.mk test-nitsche-rheology  # V0: ported operator tests O1-O7
#   make -f tests/viscoelastic/viscoelastic.mk test-viscoelastic      # every V phase test
include tests/2d_inplane/2d_inplane.mk

RHEOLOGY_HEADERS = $(wildcard rheology/*.hpp)
VISCO_LOG_DIR = $(or $(TMPDIR),/tmp)
VISCO_CXXFLAGS = $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) $(TOML_FLAGS) -DSEAS_USE_MPI

# V0: the pointwise and deck tests of the rheology plan rev 1 (d8327aad),
# without the forced-rupture test (friction/forced_rupture.hpp is TPV26/27
# work that is not ported).
TEST_RHEOLOGY_SRC = tests/unit/test_rheology_pointwise.cpp \
                    tests/unit/test_rheology_initial_stress.cpp \
                    tests/unit/test_rheology_config.cpp
TEST_RHEOLOGY_OBJ = $(TEST_RHEOLOGY_SRC:.cpp=.o)
TEST_RHEOLOGY_BIN = seas_test_rheology_pointwise seas_test_rheology_initial_stress \
                    seas_test_rheology_config
$(TEST_RHEOLOGY_OBJ): %.o: $(SRC)%.cpp $(RHEOLOGY_HEADERS) tests/unit/test_macros.hpp \
                           $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(TOML_FLAGS) -c $< -o $@

$(TEST_RHEOLOGY_BIN): seas_test_%: tests/unit/test_%.o
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $< $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-rheology: $(TEST_RHEOLOGY_BIN)
	./seas_test_rheology_pointwise
	./seas_test_rheology_initial_stress
	./seas_test_rheology_config

# V0: the operator tests O1-O7 of the rheology plan rev 1 on 1, 2 and 4 ranks.
TEST_NITSCHE_RHEOLOGY_SRC = tests/verification/test_nitsche_rheology.cpp
TEST_NITSCHE_RHEOLOGY_OBJ = $(TEST_NITSCHE_RHEOLOGY_SRC:.cpp=.o)
$(TEST_NITSCHE_RHEOLOGY_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(RHEOLOGY_HEADERS) \
                             domain/nitsche_wave_operator.hpp domain/fault_stress_projection.hpp \
                             friction/fault_closure.hpp solver/nitsche_cfl.hpp \
                             $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(VISCO_CXXFLAGS) -c $< -o $@

seas_test_nitsche_rheology: $(TEST_NITSCHE_RHEOLOGY_OBJ) $(HETEROGENEOUS_MATERIAL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_NITSCHE_RHEOLOGY_OBJ) $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-nitsche-rheology: seas_test_nitsche_rheology
	for np in 1 2 4; do $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$np ./seas_test_nitsche_rheology || exit 1; done

VISCO_TESTS = test-rheology test-nitsche-rheology

# V1: FitGMB and ViscoelasticModel (no mesh, no MPI); the refusal must abort.
TEST_VISCO_POINT_SRC = tests/unit/test_viscoelastic_point.cpp
TEST_VISCO_POINT_OBJ = $(TEST_VISCO_POINT_SRC:.cpp=.o)
$(TEST_VISCO_POINT_OBJ): %.o: $(SRC)%.cpp $(RHEOLOGY_HEADERS) tests/unit/test_macros.hpp \
                             $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

seas_test_viscoelastic_point: $(TEST_VISCO_POINT_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $< $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-viscoelastic-point: seas_test_viscoelastic_point
	./seas_test_viscoelastic_point
	@if ./seas_test_viscoelastic_point --refuse 1 > $(VISCO_LOG_DIR)/visco_refuse.log 2>&1; then \
	   echo "  FAILED: refusal 1 did not abort"; exit 1; fi; \
	 grep -qF 'moduli = "reference" needs the relaxed moduli' $(VISCO_LOG_DIR)/visco_refuse.log \
	   && echo "  PASSED: refusal 1 (moduli reference)" \
	   || { tail -5 $(VISCO_LOG_DIR)/visco_refuse.log; echo "  FAILED: refusal 1 message"; exit 1; }

VISCO_TESTS += test-viscoelastic-point
test-viscoelastic: $(VISCO_TESTS)

# V2: the viscoelastic model in the coupler and the Nitsche operator, 3-D
# (T-V2 plane waves, T-V3 energy, T-V4 rest round trip, T-V5 two phases);
# T-V1 and O2 at p = 3 are in test_nitsche_rheology above.
TEST_VISCO_COUPLER_SRC = tests/viscoelastic/test_visco_coupler.cpp
TEST_VISCO_COUPLER_OBJ = $(TEST_VISCO_COUPLER_SRC:.cpp=.o)
$(TEST_VISCO_COUPLER_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(RHEOLOGY_HEADERS) \
                             domain/nitsche_wave_operator.hpp domain/fault_stress_projection.hpp \
                             solver/nitsche_cfl.hpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(VISCO_CXXFLAGS) -c $< -o $@

seas_test_visco_coupler: $(TEST_VISCO_COUPLER_OBJ) $(HETEROGENEOUS_MATERIAL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_VISCO_COUPLER_OBJ) $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-visco-coupler: seas_test_visco_coupler
	for np in 1 2; do $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$np ./seas_test_visco_coupler || exit 1; done

VISCO_TESTS += test-visco-coupler

# V3: plane strain on triangles (O2-2D, out-of-plane memory, T-V6).
TEST_VISCO_PS_SRC = tests/viscoelastic/test_visco_plane_strain.cpp
TEST_VISCO_PS_OBJ = $(TEST_VISCO_PS_SRC:.cpp=.o)
$(TEST_VISCO_PS_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(RHEOLOGY_HEADERS) \
                         domain/nitsche_wave_operator.hpp domain/fault_stress_projection.hpp \
                         solver/nitsche_cfl.hpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(VISCO_CXXFLAGS) -c $< -o $@

seas_test_visco_plane_strain: $(TEST_VISCO_PS_OBJ) $(HETEROGENEOUS_MATERIAL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_VISCO_PS_OBJ) $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-visco-plane-strain: seas_test_visco_plane_strain
	for np in 1 2; do $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$np ./seas_test_visco_plane_strain || exit 1; done

VISCO_TESTS += test-visco-plane-strain

# V4: the 2.5D width memory (T-W1 uniform mode, T-W2 2.5D waves, T-W3
# energy, T-W5 plate motion; T-W4 the two refusals).
TEST_VISCO_WIDTH_SRC = tests/viscoelastic/test_visco_width.cpp
TEST_VISCO_WIDTH_OBJ = $(TEST_VISCO_WIDTH_SRC:.cpp=.o)
$(TEST_VISCO_WIDTH_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(RHEOLOGY_HEADERS) \
                            domain/nitsche_wave_operator.hpp domain/fault_stress_projection.hpp \
                            solver/nitsche_cfl.hpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(VISCO_CXXFLAGS) -c $< -o $@

seas_test_visco_width: $(TEST_VISCO_WIDTH_OBJ) $(HETEROGENEOUS_MATERIAL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ \
	    $(TEST_VISCO_WIDTH_OBJ) $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-visco-width: seas_test_visco_width
	for np in 1 2; do $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$np ./seas_test_visco_width || exit 1; done
	@if $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_visco_width --refuse 1 > $(VISCO_LOG_DIR)/visco_width_refuse.log 2>&1; then \
	   echo "  FAILED: T-W4 refusal 1 did not abort"; exit 1; fi; \
	 grep -qF 'must be > 0 (without the 2.5D term there is no width memory' $(VISCO_LOG_DIR)/visco_width_refuse.log \
	   && echo "  PASSED: T-W4 refusal 1 (b = 0)" \
	   || { tail -5 $(VISCO_LOG_DIR)/visco_width_refuse.log; echo "  FAILED: T-W4 refusal 1 message"; exit 1; }
	@if $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_visco_width --refuse 2 > $(VISCO_LOG_DIR)/visco_width_refuse.log 2>&1; then \
	   echo "  FAILED: T-W4 refusal 2 did not abort"; exit 1; fi; \
	 grep -qF 'the operator has no 2.5D term' $(VISCO_LOG_DIR)/visco_width_refuse.log \
	   && echo "  PASSED: T-W4 refusal 2 (no 2.5D term)" \
	   || { tail -5 $(VISCO_LOG_DIR)/visco_width_refuse.log; echo "  FAILED: T-W4 refusal 2 message"; exit 1; }

VISCO_TESTS += test-visco-width

# V5: the cycle wiring ([rheology.viscoelastic] in the fd_cycle deck, the
# s_Y step cut, the hand-back report and handback_tol_pa stop, no memory in
# the PML band, the refusals).  One 800-step dynamic phase per deck, then
# the switch back (exit 7, the guard, as the control deck).
VISCO_CYCLE_TMP = $(or $(TMPDIR),/tmp)/seas_test_visco_cycle
VISCO_CYCLE_STEPS = 800
VISCO_CYCLE_ARGS = --end-phase-after $(VISCO_CYCLE_STEPS) --max-dyn-steps $(VISCO_CYCLE_STEPS) --max-steps 3
test-visco-cycle: seas_fd_cycle_driver
	@rm -rf $(VISCO_CYCLE_TMP); mkdir -p $(VISCO_CYCLE_TMP); \
	 python3 tests/viscoelastic/visco_cycle_check.py decks $(LTS_QUENCH_DECK) $(VISCO_CYCLE_TMP) || exit 1; \
	 t=$(VISCO_CYCLE_TMP); \
	 for v in off ve venow vetol vepml r_fw r_lts r_sw r_probe r_dp r_odt; do \
	    $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_fd_cycle_driver -d $$t/$$v.toml -o $$t/$$v \
	       $(VISCO_CYCLE_ARGS) > $$t/$$v.log 2>&1; echo $$? > $$t/$$v.rc; \
	 done; \
	 $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_fd_cycle_driver -d $$t/ve.toml -o $$t/ve_np2 \
	    $(VISCO_CYCLE_ARGS) > $$t/ve_np2.log 2>&1; echo $$? > $$t/ve_np2.rc; \
	 python3 tests/viscoelastic/visco_cycle_check.py check $$t

VISCO_TESTS += test-visco-cycle

# V6: attenuation under [cycle.lts] rate2 (per-cluster memory, owner-split
# force, ghost memory), on test-lts-2d5's 2.5D deck: global stepping,
# same_step, rate2 and the control without attenuation (the same steps as
# test-lts-2d5), on 1 and 2 ranks; check_ghosts in the rate2 runs; and a
# rate2 start-up run with the PML band.
VISCO_LTS_TMP = $(or $(TMPDIR),/tmp)/seas_test_visco_lts
VISCO_LTS_STEPS = 6400
VISCO_LTS_BAND_STEPS = 100
test-visco-lts: seas_fd_cycle_driver
	@mkdir -p $(VISCO_LTS_TMP); \
	 python3 tests/viscoelastic/visco_lts_check.py decks $(LTS_QUENCH_DECK) $(VISCO_LTS_TMP) || exit 1; \
	 for np in 1 2; do t=$(VISCO_LTS_TMP)/np$$np; mkdir -p $$t; \
	    for v in off glob same lts band; do n=$(VISCO_LTS_STEPS); \
	       [ $$v = lts ] && n=$$(($(VISCO_LTS_STEPS) / 8)); \
	       [ $$v = band ] && n=$(VISCO_LTS_BAND_STEPS); \
	       $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$np ./seas_fd_cycle_driver -d $(VISCO_LTS_TMP)/$$v.toml \
	          -o $$t/$$v --end-phase-after $$n --max-dyn-steps $$n --max-steps 3 > $$t/$$v.log 2>&1; rc=$$?; \
	       echo $$rc > $$t/$$v.rc; \
	       if [ $$rc -ne 0 ] && [ $$rc -ne 7 ]; then tail -15 $$t/$$v.log; \
	          echo "  FAILED: $$np rank(s), $$v: exit $$rc"; exit 1; fi; \
	    done; \
	    echo "  $$np rank(s):"; \
	    python3 tests/viscoelastic/visco_lts_check.py check $$t $(LTS_QUENCH_PREFIX) || exit 1; \
	 done

VISCO_TESTS += test-visco-lts

# tests/pml/pml.mk - build file of the fully-dynamic-pml branch's tests
# (PML design plan rev2, Section 7).  It includes the 2D in-plane test file
# (which includes cycle/cycle.mk and the shared Makefile unchanged), so every
# target there works through it too.
#
# Run from miniapps/seas:
#   make -f tests/pml/pml.mk test-pml-band     # P1: PmlBand, C0, band integrators, patch test
#   make -f tests/pml/pml.mk test-pml-qd       # P1: the QD band on Q0 (T1, D3), 1 and 2 ranks
#   make -f tests/pml/pml.mk test-pml-dyn      # P2: the dynamic PML (S0, S1, T2-T4), 1 and 2 ranks
#   make -f tests/pml/pml.mk test-pml-switch   # P2/P3: DynPhase with the band on D0 (also KV, sponge)
#   make -f tests/pml/pml.mk test-pml-decks    # P3: the hand-back decks' switch set-up (Mac)

include tests/2d_inplane/2d_inplane.mk

PML_HEADERS = domain/pml_band.hpp integrator/pml_elasticity_integrators.hpp \
              integrator/dg_elasticity_ip_combined_integrator.hpp domain/pml_memory.hpp

TEST_PML_BAND_SRC = tests/pml/test_pml_band.cpp
TEST_PML_BAND_OBJ = $(TEST_PML_BAND_SRC:.cpp=.o)
$(TEST_PML_BAND_OBJ): %.o: $(SRC)%.cpp $(PML_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

seas_test_pml_band: $(TEST_PML_BAND_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PML_BAND_OBJ) $(MFEM_LIBS)

test-pml-band: seas_test_pml_band
	./seas_test_pml_band

TEST_PML_QD_SRC = tests/pml/test_pml_qd.cpp
TEST_PML_QD_OBJ = $(TEST_PML_QD_SRC:.cpp=.o)
$(TEST_PML_QD_OBJ): %.o: $(SRC)%.cpp $(PML_HEADERS) $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                     $(CYCLE_SHARED_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_pml_qd: $(TEST_PML_QD_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PML_QD_OBJ) $(CYCLE_LINK_OBJS) $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-pml-qd: seas_test_pml_qd
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_pml_qd
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_pml_qd

TEST_PML_DYN_SRC = tests/pml/test_pml_dyn.cpp
TEST_PML_DYN_OBJ = $(TEST_PML_DYN_SRC:.cpp=.o)
$(TEST_PML_DYN_OBJ): %.o: $(SRC)%.cpp $(PML_HEADERS) $(SEAS_HEADERS) $(INPLANE_HEADERS) \
                     $(NITSCHE_2D_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -DSEAS_USE_MPI -c $< -o $@

seas_test_pml_dyn: $(TEST_PML_DYN_OBJ) $(HETEROGENEOUS_MATERIAL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PML_DYN_OBJ) $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# S0, S1, T2, T2KV on 1 and 2 ranks; T3/T4 (the pulse runs) on 1 rank
test-pml-dyn: seas_test_pml_dyn
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_pml_dyn
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_pml_dyn --only S0,S1,T2,T2KV

TEST_PML_SWITCH_SRC = tests/pml/test_pml_switch.cpp
TEST_PML_SWITCH_OBJ = $(TEST_PML_SWITCH_SRC:.cpp=.o)
$(TEST_PML_SWITCH_OBJ): %.o: $(SRC)%.cpp $(PML_HEADERS) $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                        $(CYCLE_SHARED_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_pml_switch: $(TEST_PML_SWITCH_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PML_SWITCH_OBJ) $(CYCLE_LINK_OBJS) $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-pml-switch: seas_test_pml_switch
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_pml_switch
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_pml_switch
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_pml_switch --kv-eta 0.1
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_pml_switch --kind sponge --kv-eta 0.1

# P3: the hand-back decks (simple_cases/configs/D0_straight_hb_*): every band
# deck read (cloud), and its switch set-up and 20 steps (the Mac: a start-up check)
PML_P3_DECKS = D0_straight_hb_pml0p005 D0_straight_hb_pml0p01 D0_straight_hb_pml0p02 \
               D0_straight_hb_sponge
test-pml-decks-parse: seas_test_pml_switch
	@for d in $(PML_P3_DECKS); do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_pml_switch --deck simple_cases/configs/$$d.toml \
	      --as-is --parse-only || exit 1; done
test-pml-decks: seas_test_pml_switch
	@for d in $(PML_P3_DECKS); do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_pml_switch --deck simple_cases/configs/$$d.toml \
	      --as-is --dyn-steps 20 || exit 1; done

.PHONY: test-pml-band test-pml-qd test-pml-dyn test-pml-switch test-pml-decks-parse test-pml-decks

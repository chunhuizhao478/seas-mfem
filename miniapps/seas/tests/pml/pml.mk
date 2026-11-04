# tests/pml/pml.mk - build file of the fully-dynamic-pml branch's tests
# (PML design plan rev2, Section 7).  It includes the 2D in-plane test file
# (which includes cycle/cycle.mk and the shared Makefile unchanged), so every
# target there works through it too.
#
# Run from miniapps/seas:
#   make -f tests/pml/pml.mk test-pml-band     # P1: PmlBand, C0, band integrators, patch test
#   make -f tests/pml/pml.mk test-pml-qd       # P1: the QD band on Q0 (T1, D3), 1 and 2 ranks

include tests/2d_inplane/2d_inplane.mk

PML_HEADERS = domain/pml_band.hpp integrator/pml_elasticity_integrators.hpp \
              integrator/dg_elasticity_ip_combined_integrator.hpp

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

.PHONY: test-pml-band test-pml-qd

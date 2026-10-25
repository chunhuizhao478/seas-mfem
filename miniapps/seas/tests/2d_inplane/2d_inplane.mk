# Copyright (c) 2010-2026, Lawrence Livermore National Security, LLC.
#
# tests/2d_inplane/2d_inplane.mk — build file of the 2D in-plane branch's
# tests (design rev 3).  It includes cycle/cycle.mk (which includes the shared
# Makefile unchanged), so every cycle and shared target works through it too.
#
# Run from miniapps/seas:
#   make -f tests/2d_inplane/2d_inplane.mk test-fault-frame-field   # phase A: A1 + A6 (rotation)
#   make -f tests/2d_inplane/2d_inplane.mk test-plate-loading       # phase A: A4 ([cycle.loading])
#   make -f tests/2d_inplane/2d_inplane.mk test-okada-dislocation   # phase A: the Okada check
#   make -f tests/2d_inplane/2d_inplane.mk 2d-inplane-all            # every phase A test above
# The 3D parity set (every phase must leave it bit for bit unchanged) is the
# script tests/2d_inplane/parity_3d.sh (see its header).

include cycle/cycle.mk

INPLANE_HEADERS = fault/fault_frame_field.hpp fault/fault_orientation.hpp fault/fault_frame.hpp \
                  fault/fault_basis.hpp domain/domain_config.hpp domain/elasticity_operator.hpp \
                  domain/elasticity_operator_setup.inl domain/elasticity_operator_traction.inl \
                  domain/elasticity_operator_assembly.inl tests/2d_inplane/frame_box.hpp

TEST_FAULT_FRAME_FIELD_SRC = tests/2d_inplane/test_fault_frame_field.cpp
TEST_FAULT_FRAME_FIELD_OBJ = $(TEST_FAULT_FRAME_FIELD_SRC:.cpp=.o)
$(TEST_FAULT_FRAME_FIELD_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(INPLANE_HEADERS) \
                                    $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -DSEAS_USE_MPI \
	    -c $< -o $@

seas_test_fault_frame_field: $(TEST_FAULT_FRAME_FIELD_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_FRAME_FIELD_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# F1-F5 on 1, 2 and 4 ranks; the seam rule: the consistent L-shaped surface
# passes, the inconsistent one must exit non-zero with the seam-rule message.
test-fault-frame-field: seas_test_fault_frame_field
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_fault_frame_field
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_fault_frame_field
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_fault_frame_field
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_fault_frame_field --seam good
	@out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_fault_frame_field --seam bad 2>&1); rc=$$?; \
	 if [ $$rc -ne 0 ] && echo "$$out" | grep -q "seam rule"; then \
	    echo "  PASSED: seam bad: refused (exit $$rc, seam-rule message)"; \
	 else echo "$$out" | tail -20; echo "  FAILED: seam bad: exit $$rc, seam-rule message not found"; exit 1; fi

# A4: the plate loading (needs the deck loader: the cycle's flags and objects).
TEST_PLATE_LOADING_SRC = tests/2d_inplane/test_plate_loading.cpp
TEST_PLATE_LOADING_OBJ = $(TEST_PLATE_LOADING_SRC:.cpp=.o)
$(TEST_PLATE_LOADING_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                                $(CYCLE_SHARED_HEADERS) $(INPLANE_HEADERS) \
                                $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_plate_loading: $(TEST_PLATE_LOADING_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PLATE_LOADING_OBJ) $(CYCLE_LINK_OBJS) \
	    $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# L1-L6 on 1 and 2 ranks; then each refusal must exit non-zero with its message.
PLATE_REFUSALS = 1:takes no other key|2:must be unit vectors|3:slip_dir must lie in the plane|4:lies in the fault plane|5:needs [fault_geometry].ref_normal|6:must be "bp5" or "plate"
test-plate-loading: seas_test_plate_loading
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_plate_loading
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_plate_loading
	@list='$(PLATE_REFUSALS)'; set -f; IFS='|'; for r in $$list; do n=$${r%%:*}; msg=$${r#*:}; \
	   out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_plate_loading --refuse $$n 2>&1); rc=$$?; \
	   if [ $$rc -ne 0 ] && printf '%s' "$$out" | grep -qF "$$msg"; then \
	      echo "  PASSED: refusal $$n: exit $$rc, \"$$msg\""; \
	   else printf '%s\n' "$$out" | tail -15; echo "  FAILED: refusal $$n: exit $$rc, \"$$msg\" not found"; exit 1; fi; \
	done

# The static dislocation check against Okada (1992): the full set on 2 ranks
# (largest solve: n = 8, p = 2, 0.29 M DOFs, about 7 GB with MUMPS BLR), the
# n = 4 cases on 1 rank.
TEST_OKADA_SRC = tests/2d_inplane/test_okada_dislocation.cpp
TEST_OKADA_OBJ = $(TEST_OKADA_SRC:.cpp=.o)
$(TEST_OKADA_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(INPLANE_HEADERS) \
                        tests/2d_inplane/okada92.hpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -DSEAS_USE_MPI \
	    -c $< -o $@

seas_test_okada_dislocation: $(TEST_OKADA_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_OKADA_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-okada-dislocation: seas_test_okada_dislocation
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_okada_dislocation
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_okada_dislocation --quick

2d-inplane-all: test-fault-frame-field test-plate-loading test-okada-dislocation

.PHONY: test-fault-frame-field test-plate-loading test-okada-dislocation 2d-inplane-all

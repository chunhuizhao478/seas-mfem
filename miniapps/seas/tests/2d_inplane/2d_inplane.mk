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
#   make -f tests/2d_inplane/2d_inplane.mk test-2d-fault-layout     # phase B: the 2D switch
#   make -f tests/2d_inplane/2d_inplane.mk test-2d-deck             # phase B: 2D deck vectors
#   make -f tests/2d_inplane/2d_inplane.mk 2d-inplane-all            # every test above
# The 3D parity set (every phase must leave it bit for bit unchanged) is the
# script tests/2d_inplane/parity_3d.sh (see its header).

include cycle/cycle.mk

INPLANE_HEADERS = fault/fault_frame_field.hpp fault/fault_orientation.hpp fault/fault_frame.hpp \
                  fault/fault_basis.hpp domain/domain_config.hpp domain/elasticity_operator.hpp \
                  domain/elasticity_operator_setup.inl domain/elasticity_operator_traction.inl \
                  domain/elasticity_operator_assembly.inl domain/elasticity_operator_fault_ops.inl \
                  domain/elasticity_operator_penalty.inl domain/elasticity_operator_qd.inl \
                  fault/face_quadrature.hpp cycle/scaled_mesh.hpp tests/2d_inplane/frame_box.hpp

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

# Phase B: the 2D switch on 2D triangle meshes built in the test (loader,
# censuses, slot rule, fault DOF layout).
TEST_2D_FAULT_LAYOUT_SRC = tests/2d_inplane/test_2d_fault_layout.cpp
TEST_2D_FAULT_LAYOUT_OBJ = $(TEST_2D_FAULT_LAYOUT_SRC:.cpp=.o)
$(TEST_2D_FAULT_LAYOUT_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(INPLANE_HEADERS) \
                                  $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -DSEAS_USE_MPI \
	    -c $< -o $@

seas_test_2d_fault_layout: $(TEST_2D_FAULT_LAYOUT_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_2D_FAULT_LAYOUT_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# L1, C1-C3, S1, D1 on 1, 2 and 4 ranks; then each refusal must exit
# non-zero with its message.
LAYOUT_2D_REFUSALS = 1:triangle meshes only|2:needs the fault frame field|3:must lie in the mesh plane ([fault_geometry] up|4:FaultFrameField (2D): ref_normal|5:not available on a 2D mesh yet|6:up vector is nearly collinear|7:BR2 is 3D only
test-2d-fault-layout: seas_test_2d_fault_layout
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_2d_fault_layout
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_2d_fault_layout
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_2d_fault_layout
	@list='$(LAYOUT_2D_REFUSALS)'; set -f; IFS='|'; for r in $$list; do n=$${r%%:*}; msg=$${r#*:}; \
	   out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_2d_fault_layout --refuse $$n 2>&1); rc=$$?; \
	   if [ $$rc -ne 0 ] && printf '%s' "$$out" | grep -qF "$$msg"; then \
	      echo "  PASSED: refusal $$n: exit $$rc, \"$$msg\""; \
	   else printf '%s\n' "$$out" | tail -15; echo "  FAILED: refusal $$n: exit $$rc, \"$$msg\" not found"; exit 1; fi; \
	done; rm -f test_2d_fault_layout_quad.mesh

# Phase B (part 2): the deck's vectors in 2D (needs the deck loader: the
# cycle's flags and objects, as test-plate-loading).
TEST_2D_DECK_SRC = tests/2d_inplane/test_2d_deck.cpp
TEST_2D_DECK_OBJ = $(TEST_2D_DECK_SRC:.cpp=.o)
$(TEST_2D_DECK_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                          $(CYCLE_SHARED_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_2d_deck: $(TEST_2D_DECK_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_2D_DECK_OBJ) $(CYCLE_LINK_OBJS) \
	    $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# K1-K4 on 1 rank; each refusal must exit non-zero with its message; then
# both drivers, given a 3-D mesh and a 2-entry [fault_geometry].ref_normal,
# must refuse right after reading the mesh (the Nitsche deck's fault
# ParaView output off, so that an MFEM without HDF5 gets that far).
DECK_2D_REFUSALS = 1:ref_normal has 2 entries, but the mesh is 3-D|2:up has 2 entries, but the mesh is 3-D|3:x has 2 entries, but the mesh is 3-D|4:plane_origin has 3 entries, but the mesh is 2-D|5:up is required on a 2-D mesh|6:must be [0, 0, +-1]|7:ref_normal has 3 entries, but the mesh is 2-D|8:ref_normal is required on a 2-D mesh|9:must be a 2- or 3-element array|10:must be an array of 3 numbers (2 on a 2-D mesh)
DECK_2D_TMP = $(or $(TMPDIR),/tmp)/seas_test_2d_deck_drivers
DECK_2D_DRIVER_MSG = ref_normal has 2 entries, but the mesh is 3-D
test-2d-deck: seas_test_2d_deck seas_fd_cycle_driver seas_nitsche_dyn_driver
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_2d_deck
	@list='$(DECK_2D_REFUSALS)'; set -f; IFS='|'; for r in $$list; do n=$${r%%:*}; msg=$${r#*:}; \
	   out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_2d_deck --refuse $$n 2>&1); rc=$$?; \
	   if [ $$rc -ne 0 ] && printf '%s' "$$out" | grep -qF "$$msg"; then \
	      echo "  PASSED: refusal $$n: exit $$rc, \"$$msg\""; \
	   else printf '%s\n' "$$out" | tail -15; echo "  FAILED: refusal $$n: exit $$rc, \"$$msg\" not found"; exit 1; fi; \
	done
	@rm -rf $(DECK_2D_TMP); mkdir -p $(DECK_2D_TMP); \
	 sed 's/^ref_normal      = \[0.0, -1.0, 0.0\]$$/ref_normal      = [0.0, -1.0]/' \
	    config/fd_cycle/bp5_fd_first_event.toml > $(DECK_2D_TMP)/fd.toml; \
	 sed -e 's/^ref_normal        = \[0.0, -1.0, 0.0\]$$/ref_normal        = [0.0, -1.0]/' \
	     -e 's/^paraview_fault         = "hdf5"$$/paraview_fault         = "off"/' \
	    tpv205/configs/tpv205_nitsche.toml > $(DECK_2D_TMP)/nitsche.toml; \
	 for f in fd nitsche; do grep -q '^ref_normal *= \[0.0, -1.0\]$$' $(DECK_2D_TMP)/$$f.toml \
	    || { echo "  FAILED: the $$f deck's ref_normal line moved"; exit 1; }; done; \
	 out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_fd_cycle_driver -d $(DECK_2D_TMP)/fd.toml \
	    -o $(DECK_2D_TMP)/fd_out 2>&1); rc=$$?; \
	 if [ $$rc -ne 0 ] && printf '%s' "$$out" | grep -qF "$(DECK_2D_DRIVER_MSG)" \
	    && printf '%s' "$$out" | grep -q "^  mesh: "; then \
	    echo "  PASSED: seas_fd_cycle_driver: exit $$rc, \"$(DECK_2D_DRIVER_MSG)\""; \
	 else printf '%s\n' "$$out" | tail -15; echo "  FAILED: seas_fd_cycle_driver: exit $$rc"; exit 1; fi; \
	 out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_nitsche_dyn_driver --config $(DECK_2D_TMP)/nitsche.toml \
	    --mesh bp5/mesh/reference/bp5_tandem_exact.msh --output-dir $(DECK_2D_TMP)/nitsche_out 2>&1); rc=$$?; \
	 if [ $$rc -ne 0 ] && printf '%s' "$$out" | grep -qF "$(DECK_2D_DRIVER_MSG)"; then \
	    echo "  PASSED: seas_nitsche_dyn_driver: exit $$rc, \"$(DECK_2D_DRIVER_MSG)\""; \
	 else printf '%s\n' "$$out" | tail -15; echo "  FAILED: seas_nitsche_dyn_driver: exit $$rc"; exit 1; fi; \
	 rm -rf $(DECK_2D_TMP)

2d-inplane-all: test-fault-frame-field test-plate-loading test-okada-dislocation test-2d-fault-layout \
                test-2d-deck

.PHONY: test-fault-frame-field test-plate-loading test-okada-dislocation test-2d-fault-layout \
        test-2d-deck 2d-inplane-all

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
#   make -f tests/2d_inplane/2d_inplane.mk test-2d-nitsche          # phase C: 2D Nitsche operator
#   make -f tests/2d_inplane/2d_inplane.mk test-2d-driver           # phase C: the driver, TPV205-2D
#   make -f tests/2d_inplane/2d_inplane.mk test-thermal-pressurization  # phase D: the TP module
#   make -f tests/2d_inplane/2d_inplane.mk test-normal-stress-hook  # phase D: the hook, TP on a rupture
#   make -f tests/2d_inplane/2d_inplane.mk test-2d-driver-tp        # phase D: the driver, TPV105-2D
#   make -f tests/2d_inplane/2d_inplane.mk test-2d-dislocation      # phase E: plane-strain dislocation
#   make -f tests/2d_inplane/2d_inplane.mk test-prestress-bp3       # phases E, F: BP3 prestress (QD, FD), decks, outputs
#   make -f tests/2d_inplane/2d_inplane.mk test-2d-absorbing        # BP3's absorbing walls ([cycle.absorbing], 2D)
#   make -f tests/2d_inplane/2d_inplane.mk test-absorbing-all       # faces = "all" (every outer face absorbs)
#   make -f tests/2d_inplane/2d_inplane.mk test-snapshots-volume    # [cycle.snapshots.volume] (2-D and 3-D)
#   make -f tests/2d_inplane/2d_inplane.mk test-tensor-prestress    # [stress] constant_tensor in the cycle
#   make -f tests/2d_inplane/2d_inplane.mk test-handback-table      # [cycle.output].handback_table
#   make -f tests/2d_inplane/2d_inplane.mk 2d-inplane-all            # every test above
# The 3D parity set (every phase must leave it bit for bit unchanged) is the
# script tests/2d_inplane/parity_3d.sh (see its header).  The tensor
# prestress's rotation run (the map-view deck and mesh turned by 30 degrees,
# the QD regime and the first event) is tests/2d_inplane/tensor_rotation_run.sh.

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

# L1-L8 on 1 and 2 ranks; then each refusal must exit non-zero with its message.
PLATE_REFUSALS = 1:takes no other key|2:must be unit vectors|3:slip_dir must lie in the plane|4:lies in the fault plane|5:needs [fault_geometry].ref_normal|6:must be "bp5" or "plate"|7:have a vertex beyond the band
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

# Phase E: the static dislocation check in 2D, plane strain, against Okada
# (1992) in the middle plane of a 2000 km long patch: the full set on 2
# ranks, the n = 4 cases on 1 rank.
TEST_DISL2D_SRC = tests/2d_inplane/test_2d_dislocation.cpp
TEST_DISL2D_OBJ = $(TEST_DISL2D_SRC:.cpp=.o)
$(TEST_DISL2D_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(INPLANE_HEADERS) \
                         tests/2d_inplane/okada92.hpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -DSEAS_USE_MPI \
	    -c $< -o $@

seas_test_2d_dislocation: $(TEST_DISL2D_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_DISL2D_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-2d-dislocation: seas_test_2d_dislocation
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_2d_dislocation
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_2d_dislocation --quick

# Phases E and F: BP3's prestress ([cycle.prestress].a = "amax", with and
# without eta |V|), the initial state of the BP3-QD and BP3-FD decks, the
# QD decks' outputs at t = 0, and the station
# location on segment faces (the deck loader, the QD phase and the output
# writers: the cycle's flags and objects, as test-2d-deck).
TEST_PRESTRESS_BP3_SRC = tests/2d_inplane/test_prestress_bp3.cpp
TEST_PRESTRESS_BP3_OBJ = $(TEST_PRESTRESS_BP3_SRC:.cpp=.o)
$(TEST_PRESTRESS_BP3_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                                $(CYCLE_SHARED_HEADERS) $(INPLANE_HEADERS) \
                                $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_prestress_bp3: $(TEST_PRESTRESS_BP3_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PRESTRESS_BP3_OBJ) $(CYCLE_LINK_OBJS) \
	    $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

# B1-B7 on 1 and 2 ranks; then each refusal must exit non-zero with its message.
PRESTRESS_BP3_REFUSALS = 1:is refused with ic = "gold"|2:must be "local" or "amax"|3:needs a_max > 0|4:eta_V_init = false with ic = "gold"|5:takes no axis2 or half2 on a 2-D mesh|6:needs axis2 and half2 on a 3-D mesh|7:half2 is required|8:from its node (DOF order != basis order)
test-prestress-bp3: seas_test_prestress_bp3
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_prestress_bp3
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_prestress_bp3
	@list='$(PRESTRESS_BP3_REFUSALS)'; set -f; IFS='|'; for r in $$list; do n=$${r%%:*}; msg=$${r#*:}; \
	   out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_prestress_bp3 --refuse $$n 2>&1); rc=$$?; \
	   if [ $$rc -ne 0 ] && printf '%s' "$$out" | grep -qF "$$msg"; then \
	      echo "  PASSED: refusal $$n: exit $$rc, \"$$msg\""; \
	   else printf '%s\n' "$$out" | tail -15; echo "  FAILED: refusal $$n: exit $$rc, \"$$msg\" not found"; exit 1; fi; \
	done

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
LAYOUT_2D_REFUSALS = 1:triangle meshes only|2:needs the fault frame field|3:must lie in the mesh plane ([fault_geometry] up|4:FaultFrameField (2D): ref_normal|5:need a 3-D mesh (they key fault DOFs by position)|6:up vector is nearly collinear|7:BR2 is 3D only
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

# Phase C (part 1): the Nitsche dynamic operator on 2D triangle meshes built
# in the test (plane waves, the CFL constant on triangles, the step, the
# dashpot, an in-plane rupture).  N1-N5 on 1 rank, then the parts that must
# not depend on the partition (N3, N5) on 2 and 4 ranks.
NITSCHE_2D_HEADERS = domain/nitsche_wave_operator.hpp domain/nitsche_mass.hpp \
                     domain/nitsche_material_sample.hpp solver/nitsche_cfl.hpp \
                     integrator/dg_absorbing_bc_integrator.hpp friction/slip_weakening_closure.hpp
TEST_2D_NITSCHE_SRC = tests/2d_inplane/test_2d_nitsche.cpp
TEST_2D_NITSCHE_OBJ = $(TEST_2D_NITSCHE_SRC:.cpp=.o)
$(TEST_2D_NITSCHE_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(INPLANE_HEADERS) $(NITSCHE_2D_HEADERS) \
                             $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -DSEAS_USE_MPI \
	    -c $< -o $@

seas_test_2d_nitsche: $(TEST_2D_NITSCHE_OBJ) $(HETEROGENEOUS_MATERIAL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_2D_NITSCHE_OBJ) $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-2d-nitsche: seas_test_2d_nitsche
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_2d_nitsche
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_2d_nitsche --np
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_2d_nitsche --np

# Phase C (part 2): the Nitsche driver on the TPV205-2D deck and its
# committed 100 m mesh, 0.1 s: 1 rank and 4 ranks must run (the triangle
# count and the 2-D stable constant in the log) and give the same station
# files to DRIVER_2D_RTOL; 4 ranks with local time stepping must run; then
# each 2-D refusal (a station off the plane, a z bound in [[stress_region]],
# the fault ParaView output, a nucleation centre off the plane) must exit
# non-zero with its message.
DRIVER_2D_DECK = tpv205/configs/tpv205_2d_nitsche.toml
DRIVER_2D_TMP = $(or $(TMPDIR),/tmp)/seas_test_2d_driver
DRIVER_2D_RTOL = 1e-9
DRIVER_2D_RUN = ./seas_nitsche_dyn_driver --config $(DRIVER_2D_DECK) --tfinal 0.1 --log-every 1000
DRIVER_2D_REFUSALS = station:2-D stations lie in the (x, y) plane|region:z_min_m / z_max_m on a 2-D mesh|pvfault:on a 2-D mesh: fault.vtkhdf holds triangles|nucleation:both must be 0 (the (x, y) plane)
test-2d-driver: seas_nitsche_dyn_driver
	@rm -rf $(DRIVER_2D_TMP); mkdir -p $(DRIVER_2D_TMP)
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 $(DRIVER_2D_RUN) --output-dir $(DRIVER_2D_TMP)/np1 > $(DRIVER_2D_TMP)/np1.log 2>&1 \
	   || { tail -20 $(DRIVER_2D_TMP)/np1.log; echo "  FAILED: TPV205-2D on 1 rank"; exit 1; }
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 $(DRIVER_2D_RUN) --output-dir $(DRIVER_2D_TMP)/np4 > $(DRIVER_2D_TMP)/np4.log 2>&1 \
	   || { tail -20 $(DRIVER_2D_TMP)/np4.log; echo "  FAILED: TPV205-2D on 4 ranks"; exit 1; }
	@for f in np1 np4; do grep -q ": 40526 triangles; p = 1;" $(DRIVER_2D_TMP)/$$f.log \
	    && grep -q "x C(p) 0.215 ->" $(DRIVER_2D_TMP)/$$f.log \
	    || { echo "  FAILED: $$f log lacks the triangle count or the 2-D constant"; exit 1; }; done; \
	 echo "  PASSED: TPV205-2D on 1 and 4 ranks (40526 triangles, C(1) = 0.215)"
	python3 tests/2d_inplane/compare_station_runs.py $(DRIVER_2D_TMP)/np1 $(DRIVER_2D_TMP)/np4 --rtol $(DRIVER_2D_RTOL)
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 $(DRIVER_2D_RUN) --lts rate2 --output-dir $(DRIVER_2D_TMP)/lts > $(DRIVER_2D_TMP)/lts.log 2>&1 \
	   || { tail -20 $(DRIVER_2D_TMP)/lts.log; echo "  FAILED: TPV205-2D with --lts rate2"; exit 1; }
	@grep -q "^  LTS: .* ticks; patch element kicks" $(DRIVER_2D_TMP)/lts.log \
	    || { echo "  FAILED: the LTS run has no tick report"; exit 1; }; \
	 echo "  PASSED: TPV205-2D with --lts rate2 on 4 ranks"
	@d=$(DRIVER_2D_TMP); printf 'name,x,y,z\nup_dip,0,0,-7500\n' > $$d/st.csv; \
	 sed -e "s|^file             = .*|file             = \"$$d/st.csv\"|" $(DRIVER_2D_DECK) > $$d/station.toml; \
	 awk '{print} /^x_max_m     =  1500.0$$/ && !n {print "z_min_m     = -9000.0"; n = 1}' $(DRIVER_2D_DECK) > $$d/region.toml; \
	 sed -e 's/^paraview_enabled       = false$$/paraview_enabled       = true/' \
	     -e 's/^paraview_fault         = "off".*/paraview_fault         = "hdf5"/' $(DRIVER_2D_DECK) > $$d/pvfault.toml; \
	 { cat $(DRIVER_2D_DECK); printf '\n[nucleation]\nkind = "gradual_overstress_compact_circular"\n'; \
	   printf '[nucleation.gradual_overstress_compact_circular]\ncenter_x_m = 0.0\ncenter_y_m = 0.0\n'; \
	   printf 'center_z_m = -7500.0\nradius_m = 1500.0\ndelta_tau_pa = 1.0e6\nT_nuc_s = "1.0s"\n'; } > $$d/nucleation.toml; \
	 for f in station region pvfault nucleation; do cmp -s $(DRIVER_2D_DECK) $$d/$$f.toml \
	    && { echo "  FAILED: the $$f edit did not change the deck"; exit 1; }; done; \
	 list='$(DRIVER_2D_REFUSALS)'; set -f; IFS='|'; for r in $$list; do n=$${r%%:*}; msg=$${r#*:}; \
	   out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_nitsche_dyn_driver --config $$d/$$n.toml \
	      --tfinal 0.01 --output-dir $$d/out_$$n 2>&1); rc=$$?; \
	   if [ $$rc -ne 0 ] && printf '%s' "$$out" | grep -qF "$$msg"; then \
	      echo "  PASSED: refusal $$n: exit $$rc, \"$$msg\""; \
	   else printf '%s\n' "$$out" | tail -15; echo "  FAILED: refusal $$n: exit $$rc, \"$$msg\" not found"; exit 1; fi; \
	 done; rm -rf $$d

# Phase D (part 1): the thermal pressurization module alone (P1-P8 against
# the analytic solutions of the 1-D problem), then each refusal must exit
# non-zero with its message.
TEST_TP_SRC = tests/2d_inplane/test_thermal_pressurization.cpp
TEST_TP_OBJ = $(TEST_TP_SRC:.cpp=.o)
$(TEST_TP_OBJ): %.o: $(SRC)%.cpp fault/thermal_pressurization.hpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -c $< -o $@

seas_test_thermal_pressurization: $(TEST_TP_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TP_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

TP_REFUSALS = 1:ThermalPressurization::Advance: dt = 0|2:the shear heating must be finite and >= 0|3:the shear heating must be finite and >= 0|4:need dz0 > 0, growth >= 1 and z_max > dz0|5:need alpha_th, alpha_hy, rho_c, h > 0|6:has no parameters (SetParameters)|7:ThermalPressurization::SetState: size
test-thermal-pressurization: seas_test_thermal_pressurization
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_thermal_pressurization
	@list='$(TP_REFUSALS)'; set -f; IFS='|'; for r in $$list; do n=$${r%%:*}; msg=$${r#*:}; \
	   out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_thermal_pressurization --refuse $$n 2>&1); rc=$$?; \
	   if [ $$rc -ne 0 ] && printf '%s' "$$out" | grep -qF "$$msg"; then \
	      echo "  PASSED: refusal $$n: exit $$rc, \"$$msg\""; \
	   else printf '%s\n' "$$out" | tail -15; echo "  FAILED: refusal $$n: exit $$rc, \"$$msg\" not found"; exit 1; fi; \
	done

# Phase D (part 2): the Nitsche operator's normal-stress hook and thermal
# pressurization behind it, on the in-plane rupture of N5 (H1-H4), then each
# refusal must exit non-zero with its message.
TP_HOOK_HEADERS = fault/fault_normal_stress_modifier.hpp fault/thermal_pressurization.hpp \
                  fault/thermal_pressurization_modifier.hpp
TEST_NSH_SRC = tests/2d_inplane/test_normal_stress_hook.cpp
TEST_NSH_OBJ = $(TEST_NSH_SRC:.cpp=.o)
$(TEST_NSH_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(INPLANE_HEADERS) $(NITSCHE_2D_HEADERS) \
                      $(TP_HOOK_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -DSEAS_USE_MPI \
	    -c $< -o $@

seas_test_normal_stress_hook: $(TEST_NSH_OBJ) $(HETEROGENEOUS_MATERIAL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_NSH_OBJ) $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

NSH_REFUSALS = 1:SetNormalStressModifier: the modifier covers|2:the normal stress the law sees left its range
test-normal-stress-hook: seas_test_normal_stress_hook
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_normal_stress_hook
	@list='$(NSH_REFUSALS)'; set -f; IFS='|'; for r in $$list; do n=$${r%%:*}; msg=$${r#*:}; \
	   out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_normal_stress_hook --refuse $$n 2>&1); rc=$$?; \
	   if [ $$rc -ne 0 ] && printf '%s' "$$out" | grep -qF "$$msg"; then \
	      echo "  PASSED: refusal $$n: exit $$rc, \"$$msg\""; \
	   else printf '%s\n' "$$out" | tail -15; echo "  FAILED: refusal $$n: exit $$rc, \"$$msg\" not found"; exit 1; fi; \
	done

# Phase D (part 2): the Nitsche driver on the TPV105-2D deck and its
# committed 100 m mesh, 1 s (the rupture is under way, T and p rise): 1 rank
# and 4 ranks give the same station files, T and p included, to
# DRIVER_2D_RTOL; Lambda = 0 gives the deck without
# [friction.thermal_pressurization] bit for bit in every column they share
# (4 ranks both); then each deck refusal must exit non-zero with its message.
DRIVER_TP_DECK = tpv105/configs/tpv105_2d_nitsche.toml
DRIVER_TP_TMP = $(or $(TMPDIR),/tmp)/seas_test_2d_driver_tp
DRIVER_TP_RUN = ./seas_nitsche_dyn_driver --tfinal 1.0 --log-every 1000
DRIVER_TP_REFUSALS = unknown:[friction.thermal_pressurization]: unknown key|missing:[friction.thermal_pressurization] needs h; there are no defaults|scheme:scheme = "crank-nicolson"; use "backward_euler" or "crank_nicolson"
test-2d-driver-tp: seas_nitsche_dyn_driver
	@rm -rf $(DRIVER_TP_TMP); mkdir -p $(DRIVER_TP_TMP)
	@d=$(DRIVER_TP_TMP); \
	 sed -e 's/^Lambda   = 0.1e6 .*/Lambda   = 0.0/' $(DRIVER_TP_DECK) > $$d/lam0.toml; \
	 awk '/^\[friction.thermal_pressurization\]/ {skip = 1; next} /^\[/ {skip = 0} !skip' $(DRIVER_TP_DECK) > $$d/notp.toml; \
	 sed -e 's/^Lambda   = 0.1e6/lambda   = 0.1e6/' $(DRIVER_TP_DECK) > $$d/unknown.toml; \
	 sed -e '/^h        = 0.01 /d' $(DRIVER_TP_DECK) > $$d/missing.toml; \
	 sed -e 's/^scheme   = "backward_euler"/scheme   = "crank-nicolson"/' $(DRIVER_TP_DECK) > $$d/scheme.toml; \
	 for f in lam0 notp unknown missing scheme; do cmp -s $(DRIVER_TP_DECK) $$d/$$f.toml \
	    && { echo "  FAILED: the $$f edit did not change the deck"; exit 1; }; done; \
	 grep -q '^\[friction.thermal_pressurization\]' $$d/notp.toml && { echo "  FAILED: notp still has the table"; exit 1; }; \
	 grep -q '^\[nucleation\]' $$d/notp.toml || { echo "  FAILED: notp lost [nucleation]"; exit 1; }; true
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 $(DRIVER_TP_RUN) --config $(DRIVER_TP_DECK) --output-dir $(DRIVER_TP_TMP)/np1 \
	   > $(DRIVER_TP_TMP)/np1.log 2>&1 || { tail -20 $(DRIVER_TP_TMP)/np1.log; echo "  FAILED: TPV105-2D on 1 rank"; exit 1; }
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 $(DRIVER_TP_RUN) --config $(DRIVER_TP_DECK) --output-dir $(DRIVER_TP_TMP)/np4 \
	   > $(DRIVER_TP_TMP)/np4.log 2>&1 || { tail -20 $(DRIVER_TP_TMP)/np4.log; echo "  FAILED: TPV105-2D on 4 ranks"; exit 1; }
	@for f in np1 np4; do grep -q "77 unknowns per fault DOF" $(DRIVER_TP_TMP)/$$f.log \
	    && grep -q "# Columns: time slip1 slip2 V1 V2 tau1 tau2 sigma_n psi T p" $(DRIVER_TP_TMP)/$$f/tpv105_2d_station_x2_-9_x3_7.5.dat \
	    || { echo "  FAILED: $$f lacks the thermal pressurization log line or the T p columns"; exit 1; }; done; \
	 echo "  PASSED: TPV105-2D on 1 and 4 ranks (77 unknowns per fault DOF, columns T and p)"
	python3 tests/2d_inplane/compare_station_runs.py $(DRIVER_TP_TMP)/np1 $(DRIVER_TP_TMP)/np4 --rtol $(DRIVER_2D_RTOL)
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 $(DRIVER_TP_RUN) --config $(DRIVER_TP_TMP)/lam0.toml --output-dir $(DRIVER_TP_TMP)/lam0 \
	   > $(DRIVER_TP_TMP)/lam0.log 2>&1 || { tail -20 $(DRIVER_TP_TMP)/lam0.log; echo "  FAILED: Lambda = 0"; exit 1; }
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 $(DRIVER_TP_RUN) --config $(DRIVER_TP_TMP)/notp.toml --output-dir $(DRIVER_TP_TMP)/notp \
	   > $(DRIVER_TP_TMP)/notp.log 2>&1 || { tail -20 $(DRIVER_TP_TMP)/notp.log; echo "  FAILED: no thermal pressurization"; exit 1; }
	python3 tests/2d_inplane/compare_station_runs.py $(DRIVER_TP_TMP)/lam0 $(DRIVER_TP_TMP)/notp --common --rtol 0
	@d=$(DRIVER_TP_TMP); list='$(DRIVER_TP_REFUSALS)'; set -f; IFS='|'; for r in $$list; do n=$${r%%:*}; msg=$${r#*:}; \
	   out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_nitsche_dyn_driver --config $$d/$$n.toml \
	      --tfinal 0.01 --output-dir $$d/out_$$n 2>&1); rc=$$?; \
	   if [ $$rc -ne 0 ] && printf '%s' "$$out" | grep -qF "$$msg"; then \
	      echo "  PASSED: refusal $$n: exit $$rc, \"$$msg\""; \
	   else printf '%s\n' "$$out" | tail -15; echo "  FAILED: refusal $$n: exit $$rc, \"$$msg\" not found"; exit 1; fi; \
	 done; rm -rf $$d

# BP3's absorbing walls on a 2-D mesh ([cycle.absorbing] = "outer" and
# "bottom", cycle/absorbing_faces.hpp, hook H3; the 3-D test is
# tests/cycle/test_cycle_absorbing.cpp): the BP3-FD dip 60 thrust deck on its
# 1 km mesh at p = 1, on 1 and 2 ranks.
TEST_2D_ABS_SRC = tests/2d_inplane/test_2d_absorbing.cpp
TEST_2D_ABS_OBJ = $(TEST_2D_ABS_SRC:.cpp=.o)
$(TEST_2D_ABS_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                         $(CYCLE_SHARED_HEADERS) $(INPLANE_HEADERS) \
                         $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_2d_absorbing: $(TEST_2D_ABS_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_2D_ABS_OBJ) $(CYCLE_LINK_OBJS) \
	    $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-2d-absorbing: seas_test_2d_absorbing
	@for n in 1 2; do \
	   $(MFEM_MPIEXEC) $(CYCLE_P2_MPIFLAGS) $(MFEM_MPIEXEC_NP) $$n ./seas_test_2d_absorbing || exit 1; \
	done

# [cycle.absorbing] faces = "all": every outer face absorbs, no free surface
# (the UCERF 2-D box, loaded through its walls in QD; Chunhui Zhao,
# 2026-10-27).  RetagOuterFaces on the TPV205-2D box, the BP3 1 km mesh and a
# 3-D box, on 1 and 2 ranks; then the deck loader on BP3-FD d60 thrust copies:
# "all" loads (on, walls, no free surface), and "all" with
# free_surface_min_up and an unknown faces value are refused.
TEST_ABS_ALL_SRC = tests/2d_inplane/test_absorbing_all.cpp
TEST_ABS_ALL_OBJ = $(TEST_ABS_ALL_SRC:.cpp=.o)
$(TEST_ABS_ALL_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                         $(CYCLE_SHARED_HEADERS) $(INPLANE_HEADERS) \
                         $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_absorbing_all: $(TEST_ABS_ALL_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_ABS_ALL_OBJ) $(CYCLE_LINK_OBJS) \
	    $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

ABS_ALL_DECK = config/fd_cycle/bp3_fd_d60_thrust.toml
ABS_ALL_TMP = $(or $(TMPDIR),/tmp)/seas_test_absorbing_all
ABS_ALL_REFUSALS = minup:[cycle.absorbing] faces = "all" keeps no free surface; remove free_surface_min_up|unknown:[cycle.absorbing].faces must be "off", "bottom", "outer" or "all"; got "every"
test-absorbing-all: seas_test_absorbing_all
	@for n in 1 2; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_absorbing_all || exit 1; \
	done
	@rm -rf $(ABS_ALL_TMP); mkdir -p $(ABS_ALL_TMP); d=$(ABS_ALL_TMP); \
	 sed -e 's/^faces = "outer"/faces = "all"/' $(ABS_ALL_DECK) > $$d/all.toml; \
	 awk '{print} /^faces = "all"/ {print "free_surface_min_up = 0.5"}' $$d/all.toml > $$d/minup.toml; \
	 sed -e 's/^faces = "outer"/faces = "every"/' $(ABS_ALL_DECK) > $$d/unknown.toml; \
	 for f in all unknown; do cmp -s $(ABS_ALL_DECK) $$d/$$f.toml \
	    && { echo "  FAILED: the $$f edit did not change the deck"; exit 1; }; done; \
	 cmp -s $$d/all.toml $$d/minup.toml && { echo "  FAILED: the minup edit did not change the deck"; exit 1; }; \
	 out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_absorbing_all --load-deck $$d/all.toml 2>&1); rc=$$?; \
	 if [ $$rc -eq 0 ] && printf '%s' "$$out" | grep -qF 'faces = "all": on 1, walls 1, free surface 0'; then \
	    echo "  PASSED: deck faces = \"all\" loads: on, walls, no free surface"; \
	 else printf '%s\n' "$$out" | tail -15; echo "  FAILED: deck faces = \"all\""; exit 1; fi; \
	 list='$(ABS_ALL_REFUSALS)'; set -f; IFS='|'; for r in $$list; do n=$${r%%:*}; msg=$${r#*:}; \
	   out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_absorbing_all --load-deck $$d/$$n.toml 2>&1); rc=$$?; \
	   if [ $$rc -ne 0 ] && printf '%s' "$$out" | grep -qF "$$msg"; then \
	      echo "  PASSED: refusal $$n: exit $$rc, \"$$msg\""; \
	   else printf '%s\n' "$$out" | tail -15; echo "  FAILED: refusal $$n: exit $$rc, \"$$msg\" not found"; exit 1; fi; \
	 done; rm -rf $$d

# The volume snapshots ([cycle.snapshots.volume], cycle/snapshots.hpp; Chunhui
# Zhao, 2026-10-28: "we should add the volume vtkdhf option").  The writer on a
# triangle and a tetrahedron box at p = 1 and 2 (layout, values, schedule,
# restart, the box) on 1 and 2 ranks; the writer's refusals (a box with no
# element centroid, a quadrilateral mesh); then the deck loader on a BP3-FD
# d60 thrust copy with the volume on and a box, and its refusals.
TEST_SNAP_VOL_SRC = tests/2d_inplane/test_snapshots_volume.cpp
TEST_SNAP_VOL_OBJ = $(TEST_SNAP_VOL_SRC:.cpp=.o)
$(TEST_SNAP_VOL_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                          $(CYCLE_SHARED_HEADERS) $(INPLANE_HEADERS) \
                          $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_snapshots_volume: $(TEST_SNAP_VOL_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_SNAP_VOL_OBJ) $(CYCLE_LINK_OBJS) \
	    $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

SNAP_VOL_DECK = config/fd_cycle/bp3_fd_d60_thrust.toml
SNAP_VOL_TMP = $(or $(TMPDIR),/tmp)/seas_test_snapshots_volume_decks
SNAP_VOL_REFUSALS = empty-box:no element centroid lies in [cycle.snapshots.volume] box_min .. box_max|quads:the volume handles triangles only on a 2-D mesh (element 0 is Square)
SNAP_VOL_DECK_REFUSALS = onlymin:[cycle.snapshots.volume] box_min and box_max go together|order:[cycle.snapshots.volume] box_min must not exceed box_max (entry 1)
test-snapshots-volume: seas_test_snapshots_volume
	@for n in 1 2; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_snapshots_volume || exit 1; \
	done
	@list='$(SNAP_VOL_REFUSALS)'; set -f; IFS='|'; for r in $$list; do n=$${r%%:*}; msg=$${r#*:}; \
	   out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_snapshots_volume --$$n 2>&1); rc=$$?; \
	   if [ $$rc -ne 0 ] && printf '%s' "$$out" | grep -qF "$$msg"; then \
	      echo "  PASSED: refusal $$n: exit $$rc, \"$$msg\""; \
	   else printf '%s\n' "$$out" | tail -15; echo "  FAILED: refusal $$n: exit $$rc, \"$$msg\" not found"; exit 1; fi; \
	 done
	@rm -rf $(SNAP_VOL_TMP); mkdir -p $(SNAP_VOL_TMP); d=$(SNAP_VOL_TMP); \
	 [ "$$(tail -2 $(SNAP_VOL_DECK) | tr '\n' '|')" = '[cycle.snapshots.volume]|enabled = false|' ] \
	    || { echo "  FAILED: $(SNAP_VOL_DECK) does not end with the volume table off"; exit 1; }; \
	 head -n -1 $(SNAP_VOL_DECK) > $$d/base.toml; \
	 { cat $$d/base.toml; printf 'dt_qd   = 3.1536e8\ndt_dyn  = 1.0\nbox_min = [-20000.0, -30000.0]\nbox_max = [20000.0, 0.0]\n'; } > $$d/on.toml; \
	 { cat $$d/base.toml; printf 'box_min = [-20000.0, -30000.0]\n'; } > $$d/onlymin.toml; \
	 { cat $$d/base.toml; printf 'box_min = [-20000.0, 10.0]\nbox_max = [20000.0, 0.0]\n'; } > $$d/order.toml; \
	 out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_snapshots_volume --load-deck $$d/on.toml 2>&1); rc=$$?; \
	 want='[cycle.snapshots.volume] enabled 1, dt_qd 3.1536e+08, dt_dyn 1, box 1 [-20000, -30000, 0] .. [20000, 0, 0]'; \
	 if [ $$rc -eq 0 ] && printf '%s' "$$out" | grep -qF "$$want"; then \
	    echo "  PASSED: deck with the volume on and a box loads: $$want"; \
	 else printf '%s\n' "$$out" | tail -15; echo "  FAILED: deck with the volume on"; exit 1; fi; \
	 out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_snapshots_volume --load-deck $(SNAP_VOL_DECK) 2>&1); rc=$$?; \
	 if [ $$rc -eq 0 ] && printf '%s' "$$out" | grep -qF '[cycle.snapshots.volume] enabled 0,'; then \
	    echo "  PASSED: $(SNAP_VOL_DECK) (enabled = false) loads with the volume off"; \
	 else printf '%s\n' "$$out" | tail -15; echo "  FAILED: $(SNAP_VOL_DECK)"; exit 1; fi; \
	 list='$(SNAP_VOL_DECK_REFUSALS)'; set -f; IFS='|'; for r in $$list; do n=$${r%%:*}; msg=$${r#*:}; \
	   out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_snapshots_volume --load-deck $$d/$$n.toml 2>&1); rc=$$?; \
	   if [ $$rc -ne 0 ] && printf '%s' "$$out" | grep -qF "$$msg"; then \
	      echo "  PASSED: deck refusal $$n: exit $$rc, \"$$msg\""; \
	   else printf '%s\n' "$$out" | tail -15; echo "  FAILED: deck refusal $$n: exit $$rc, \"$$msg\" not found"; exit 1; fi; \
	 done; rm -rf $$d

# The cycle's tensor prestress ([stress] kind = "constant_tensor",
# cycle/prestress.hpp BuildTensorPrestress, the dynamic operator's projection
# domain/fault_stress_projection.hpp; Chunhui Zhao, 2026-10-29): T1-T6 on the
# rotated map-view line (2-D) and the rotated plane (3-D) at p = 1 and 2, on
# 1 and 2 ranks; then each refusal must exit non-zero with its message.
TEST_TENSOR_SRC = tests/2d_inplane/test_tensor_prestress.cpp
TEST_TENSOR_OBJ = $(TEST_TENSOR_SRC:.cpp=.o)
$(TEST_TENSOR_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                         $(CYCLE_SHARED_HEADERS) $(INPLANE_HEADERS) \
                         $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_tensor_prestress: $(TEST_TENSOR_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TENSOR_OBJ) $(CYCLE_LINK_OBJS) \
	    $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

TENSOR_REFUSALS = 1:the fault must be in compression|2:the tensor drives no slip there|3:[cycle.prestress].a = "amax" with [stress] kind = "constant_tensor"|4:[friction.rate_state].sigma_n_default with [stress] kind = "constant_tensor"|5:[[stress_region]] with [stress] kind = "constant_tensor"|6:the in-plane model carries no anti-plane stress|7:no state balances it at V_init
test-tensor-prestress: seas_test_tensor_prestress
	@for n in 1 2; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_tensor_prestress || exit 1; \
	done
	@list='$(TENSOR_REFUSALS)'; set -f; IFS='|'; for r in $$list; do n=$${r%%:*}; msg=$${r#*:}; \
	   out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_tensor_prestress --refuse $$n 2>&1); rc=$$?; \
	   if [ $$rc -ne 0 ] && printf '%s' "$$out" | grep -qF "$$msg"; then \
	      echo "  PASSED: refusal $$n: exit $$rc, \"$$msg\""; \
	   else printf '%s\n' "$$out" | tail -15; echo "  FAILED: refusal $$n: exit $$rc, \"$$msg\" not found"; exit 1; fi; \
	done

# [cycle.output].handback_table (the UCERF 2-D run stopped at the switch back,
# Chunhui Zhao 2026-10-30): the map-view tensor deck with its event started at
# t = 0 and the phase ended after HANDBACK_STEPS steps, on 1 rank with and
# without the table and on 2 ranks with it.  tests/2d_inplane/
# handback_table_check.py writes the decks and checks the files (see its
# header); then the deck refusal.  After 200 steps the QD balance is above
# V_on, so each run stops at the guard (exit 7) with the table written.
HANDBACK_TMP = $(or $(TMPDIR),/tmp)/seas_test_handback_table
HANDBACK_STEPS = 200
HANDBACK_DECK = config/fd_cycle/mapview_rs_tensor_fd.toml
HANDBACK_MSG = [cycle.output].handback_table must be true or false
test-handback-table: seas_fd_cycle_driver
	@rm -rf $(HANDBACK_TMP); mkdir -p $(HANDBACK_TMP); \
	 python3 tests/2d_inplane/handback_table_check.py decks $(HANDBACK_DECK) $(HANDBACK_TMP) || exit 1; \
	 for c in on:on:1 on2:on:2 off:off:1; do d=$${c%%:*}; r=$${c#*:}; deck=$${r%%:*}; np=$${r#*:}; \
	    $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$np ./seas_fd_cycle_driver -d $(HANDBACK_TMP)/$$deck.toml \
	       -o $(HANDBACK_TMP)/$$d --end-phase-after $(HANDBACK_STEPS) > $(HANDBACK_TMP)/$$d.log 2>&1; rc=$$?; \
	    if [ $$rc -eq 0 ] || [ $$rc -eq 7 ]; then \
	       echo "  PASSED: seas_fd_cycle_driver $$deck.toml on $$np rank(s): exit $$rc"; \
	    else tail -15 $(HANDBACK_TMP)/$$d.log; \
	       echo "  FAILED: seas_fd_cycle_driver $$deck.toml on $$np rank(s): exit $$rc"; exit 1; fi; \
	 done; \
	 python3 tests/2d_inplane/handback_table_check.py check $(HANDBACK_TMP) mapview_rs_tensor_fd \
	    || exit 1; \
	 out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_fd_cycle_driver -d $(HANDBACK_TMP)/bad.toml \
	    -o $(HANDBACK_TMP)/bad 2>&1); rc=$$?; \
	 if [ $$rc -ne 0 ] && printf '%s' "$$out" | grep -qF "$(HANDBACK_MSG)"; then \
	    echo "  PASSED: refusal: exit $$rc, \"$(HANDBACK_MSG)\""; \
	 else printf '%s\n' "$$out" | tail -15; echo "  FAILED: refusal: exit $$rc"; exit 1; fi; \
	 rm -rf $(HANDBACK_TMP)

2d-inplane-all: test-fault-frame-field test-plate-loading test-okada-dislocation test-2d-fault-layout \
                test-2d-deck test-2d-nitsche test-2d-driver test-thermal-pressurization \
                test-normal-stress-hook test-2d-driver-tp test-2d-dislocation test-prestress-bp3 \
                test-2d-absorbing test-absorbing-all test-snapshots-volume test-tensor-prestress \
                test-handback-table

.PHONY: test-fault-frame-field test-plate-loading test-okada-dislocation test-2d-fault-layout \
        test-2d-deck test-2d-nitsche test-2d-driver test-thermal-pressurization \
        test-normal-stress-hook test-2d-driver-tp test-2d-dislocation test-prestress-bp3 \
        test-2d-absorbing test-absorbing-all test-snapshots-volume test-tensor-prestress \
        test-handback-table 2d-inplane-all

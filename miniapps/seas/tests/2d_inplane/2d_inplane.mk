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
#   make -f tests/2d_inplane/2d_inplane.mk test-tip-zone            # [cycle.fields.a] kind = "tip_zone"
#   make -f tests/2d_inplane/2d_inplane.mk test-sigma-floor         # [cycle].sigma_n_floor_pa
#   make -f tests/2d_inplane/2d_inplane.mk test-relaxation          # [cycle.relaxation] (sigma_n' relaxation)
#   make -f tests/2d_inplane/2d_inplane.mk test-2d5-width           # [cycle.finite_width] (2.5D)
#   make -f tests/2d_inplane/2d_inplane.mk test-2d5-dyn-driver      # [finite_width] of the Nitsche driver
#   make -f tests/2d_inplane/2d_inplane.mk test-vertex-normals      # vertex-normal plan P1 (T2, T3)
#   make -f tests/2d_inplane/2d_inplane.mk test-vertex-normal-kernels # vertex-normal plan P2 (T1, T4)
#   make -f tests/2d_inplane/2d_inplane.mk test-vertex-normal-cycle # vertex-normal plan P3 (the cycle)
#   make -f tests/2d_inplane/2d_inplane.mk test-fault-opening       # [cycle.opening] version A (shear-free)
#   make -f tests/2d_inplane/2d_inplane.mk test-fault-contact       # fault opening B, B1 (opening columns)
#   make -f tests/2d_inplane/2d_inplane.mk test-fault-contact-dyn   # fault opening B, B2 (dynamic projection)
#   make -f tests/2d_inplane/2d_inplane.mk test-fault-contact-qd    # fault opening B, B3 (QD contact + taper)
#   make -f tests/2d_inplane/2d_inplane.mk test-exact-retry         # [cycle.qd] exact_retry (PETSc RK rollback)
#   make -f tests/2d_inplane/2d_inplane.mk test-qd-esdirk-points    # option I, I1 (implicit QD fault system)
#   make -f tests/2d_inplane/2d_inplane.mk test-qd-implicit-green   # option I, I2 (on the planar Green's function)
#   make -f tests/2d_inplane/2d_inplane.mk test-qd-esdirk-cycle     # option I, I3 (in the cycle driver: restart)
#   make -f tests/2d_inplane/2d_inplane.mk test-qd-esdirk-relax     # option I with [cycle.relaxation] opening = "relaxed"
#   make -f tests/2d_inplane/2d_inplane.mk test-par-dense-lu        # [cycle.qd] dense_lu = "shared" (ScaLAPACK LU)
#   make -f tests/2d_inplane/2d_inplane.mk test-pml-partition       # [cycle.pml] partition_weight (band-weighted cut)
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

# Vertex-normal plan rev1, phase P1 (2026-11-01): the averaged vertex normals
# of FaultFrameField, T2 and T3 on 1, 2 and 4 ranks; the VN-DIGEST line (every
# vertex normal in hex) must be the same at every rank count; the three
# refusals must exit non-zero with their messages.
TEST_VERTEX_NORMALS_SRC = tests/2d_inplane/test_vertex_normals.cpp
TEST_VERTEX_NORMALS_OBJ = $(TEST_VERTEX_NORMALS_SRC:.cpp=.o)
$(TEST_VERTEX_NORMALS_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(INPLANE_HEADERS) \
                                 $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -DSEAS_USE_MPI \
	    -c $< -o $@

seas_test_vertex_normals: $(TEST_VERTEX_NORMALS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_VERTEX_NORMALS_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

VN_TMP = $(or $(TMPDIR),/tmp)/seas_test_vertex_normals
test-vertex-normals: seas_test_vertex_normals
	@rm -rf $(VN_TMP) && mkdir -p $(VN_TMP)
	@set -e; for n in 1 2 4; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_vertex_normals > $(VN_TMP)/np$$n.log 2>&1 \
	      || { cat $(VN_TMP)/np$$n.log; echo "  FAILED: $$n rank(s)"; exit 1; }; \
	   grep -E "PASSED|FAILED|tests passed" $(VN_TMP)/np$$n.log; \
	   grep "^VN-DIGEST" $(VN_TMP)/np$$n.log > $(VN_TMP)/digest_np$$n.txt; done; \
	 d=$(VN_TMP)/digest_np; \
	 if cmp -s $${d}1.txt $${d}2.txt && cmp -s $${d}1.txt $${d}4.txt \
	    && [ "$$(sort -u $${d}1.txt | wc -l)" -eq 1 ] && [ "$$(wc -l < $${d}1.txt)" -eq 2 ]; \
	 then echo "  PASSED: T2g the vertex normals are the same bit for bit on 1, 2 and 4 ranks, METIS's partition and one cut along the fault"; \
	 else echo "  FAILED: T2g the VN-DIGEST lines differ between rank counts"; exit 1; fi
	@for c in "1 junction" "2 2-D meshes only" "3 must be > 0"; do \
	   n=$${c%% *}; msg=$${c#* }; \
	   out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_vertex_normals --refuse $$n 2>&1); rc=$$?; \
	   if [ $$rc -ne 0 ] && echo "$$out" | grep -q "$$msg"; then \
	      echo "  PASSED: refusal $$n: exit $$rc, message \"$$msg\""; \
	   else echo "$$out" | tail -20; echo "  FAILED: refusal $$n: exit $$rc, \"$$msg\" not found"; exit 1; fi; done

# Vertex-normal plan rev1, phase P2 (2026-11-01; decision "Per-QP frame from
# eq. 2"): the kernels with averaged vertex normals (frames, T1 regression, T4
# adjoint and the per-face operators, a static sawtooth check) on 1, 2 and 4
# ranks.
TEST_VN_KERNELS_SRC = tests/2d_inplane/test_vertex_normal_kernels.cpp
TEST_VN_KERNELS_OBJ = $(TEST_VN_KERNELS_SRC:.cpp=.o)
$(TEST_VN_KERNELS_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(INPLANE_HEADERS) \
                             $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -DSEAS_USE_MPI \
	    -c $< -o $@

seas_test_vertex_normal_kernels: $(TEST_VN_KERNELS_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_VN_KERNELS_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-vertex-normal-kernels: seas_test_vertex_normal_kernels
	@for n in 1 2 4; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_vertex_normal_kernels \
	      || { echo "  FAILED: $$n rank(s)"; exit 1; }; done

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

# K1-K5 on 1 rank; each refusal must exit non-zero with its message; then
# both drivers, given a 3-D mesh and a 2-entry [fault_geometry].ref_normal,
# must refuse right after reading the mesh (the Nitsche deck's fault
# ParaView output off, so that an MFEM without HDF5 gets that far).
DECK_2D_REFUSALS = 1:ref_normal has 2 entries, but the mesh is 3-D|2:up has 2 entries, but the mesh is 3-D|3:x has 2 entries, but the mesh is 3-D|4:plane_origin has 3 entries, but the mesh is 2-D|5:up is required on a 2-D mesh|6:must be [0, 0, +-1]|7:ref_normal has 3 entries, but the mesh is 2-D|8:ref_normal is required on a 2-D mesh|9:must be a 2- or 3-element array|10:must be an array of 3 numbers (2 on a 2-D mesh)|11:is implemented for 2-D meshes only|12:"face" or "averaged"|13:vertex_normal_max_turn_deg must be > 0
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

# The velocity-strengthening tip zone ([cycle.fields.a] kind = "tip_zone",
# cycle/fault_chains.hpp; Chunhui Zhao, 2026-10-29): G1, D1-D5 on 1 and 2
# ranks; then each refusal must exit non-zero with its message.  The check
# on a mesh file (the UCERF 2-D faults) is tests/2d_inplane/tip_zone_check.py.
TEST_TIP_ZONE_SRC = tests/2d_inplane/test_tip_zone.cpp
TEST_TIP_ZONE_OBJ = $(TEST_TIP_ZONE_SRC:.cpp=.o)
$(TEST_TIP_ZONE_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                           $(CYCLE_SHARED_HEADERS) $(INPLANE_HEADERS) \
                           $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_tip_zone: $(TEST_TIP_ZONE_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TIP_ZONE_OBJ) $(CYCLE_LINK_OBJS) \
	    $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

TIP_ZONE_REFUSALS = 1:three or more fault faces|2:form closed chains|3:has no free tip|4:"tip_zone" is for [cycle.fields.a] only|5:needs a 2-D mesh (fault lines with tips)|6:is shorter than the taper|7:taper must be in [0, width]|8:has an unknown key|9:lies on no fault segment|10:max_fraction must be > 0
test-tip-zone: seas_test_tip_zone
	@for n in 1 2; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_tip_zone || exit 1; \
	done
	@list='$(TIP_ZONE_REFUSALS)'; set -f; IFS='|'; for r in $$list; do n=$${r%%:*}; msg=$${r#*:}; \
	   out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_tip_zone --refuse $$n 2>&1); rc=$$?; \
	   if [ $$rc -ne 0 ] && printf '%s' "$$out" | grep -qF "$$msg"; then \
	      echo "  PASSED: refusal $$n: exit $$rc, \"$$msg\""; \
	   else printf '%s\n' "$$out" | tail -15; echo "  FAILED: refusal $$n: exit $$rc, \"$$msg\" not found"; exit 1; fi; \
	done

# The opt-in sigma_n' floor ([cycle].sigma_n_floor_pa; Chunhui Zhao,
# 2026-10-30): K1 (the QD kernels) and D1 (the dynamic operator's CloseDOF,
# p = 1, 2), then each refusal must exit non-zero with its message; then the
# driver on the BP3-FD d60 deck (1 km mesh): each deck refusal, and a floor
# of 1 MPa accepted and reported (3 QD steps).
TEST_SIGMA_FLOOR_SRC = tests/2d_inplane/test_sigma_floor.cpp
TEST_SIGMA_FLOOR_OBJ = $(TEST_SIGMA_FLOOR_SRC:.cpp=.o)
$(TEST_SIGMA_FLOOR_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(INPLANE_HEADERS) $(NITSCHE_2D_HEADERS) \
                              friction/device_rate_state.hpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -DSEAS_USE_MPI \
	    -c $< -o $@

seas_test_sigma_floor: $(TEST_SIGMA_FLOOR_OBJ) $(HETEROGENEOUS_MATERIAL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_SIGMA_FLOOR_OBJ) $(HETEROGENEOUS_MATERIAL_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

SIGMA_FLOOR_REFUSALS = 1:NitscheFaultSetup: sigma_floor = -1|2:the floor replaces the range check|3:sigma_n_eff is not finite|4:sigma_n_eff left its range
SIGMA_FLOOR_DECK = config/fd_cycle/bp3_fd_d60_thrust.toml
SIGMA_FLOOR_TMP = $(or $(TMPDIR),/tmp)/seas_test_sigma_floor
SIGMA_FLOOR_DECKS = zero:0.0:0.0:must be > 0 Pa|min:1.0e6:1.0e5:the floor replaces the dynamic regime|high:60.0e6:0.0:is above the prestress sigma_n
test-sigma-floor: seas_test_sigma_floor seas_fd_cycle_driver
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_sigma_floor
	@list='$(SIGMA_FLOOR_REFUSALS)'; set -f; IFS='|'; for r in $$list; do n=$${r%%:*}; msg=$${r#*:}; \
	   out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_sigma_floor --refuse $$n 2>&1); rc=$$?; \
	   if [ $$rc -ne 0 ] && printf '%s' "$$out" | grep -qF "$$msg"; then \
	      echo "  PASSED: refusal $$n: exit $$rc, \"$$msg\""; \
	   else printf '%s\n' "$$out" | tail -15; echo "  FAILED: refusal $$n: exit $$rc, \"$$msg\" not found"; exit 1; fi; \
	done
	@rm -rf $(SIGMA_FLOOR_TMP); mkdir -p $(SIGMA_FLOOR_TMP); \
	 deck() { sed -e 's|^path  = "bp3/mesh/bp3_d60_25m.msh"$$|path  = "bp3/mesh/bp3_d60_1km.msh"|' \
	              -e "s|^min_sigma_n_pa    = 0.0$$|min_sigma_n_pa    = $$3|" $(SIGMA_FLOOR_DECK) \
	          | awk -v f="sigma_n_floor_pa     = $$2" '{ print } /^penalty_factor  / { print f }' > $$1; \
	          grep -q '^path  = "bp3/mesh/bp3_d60_1km.msh"$$' $$1 && grep -q "^sigma_n_floor_pa     = $$2$$" $$1 \
	          && grep -q "^min_sigma_n_pa    = $$3$$" $$1 || { echo "  FAILED: the deck lines moved"; exit 1; }; }; \
	 n=0; list='$(SIGMA_FLOOR_DECKS)'; set -f; IFS='|'; for r in $$list; do IFS=':'; set -- $$r; IFS='|'; \
	    deck $(SIGMA_FLOOR_TMP)/$$1.toml $$2 $$3; n=$$((n + 1)); \
	    out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_fd_cycle_driver -d $(SIGMA_FLOOR_TMP)/$$1.toml \
	       -o $(SIGMA_FLOOR_TMP)/$$1_out 2>&1); rc=$$?; \
	    if [ $$rc -ne 0 ] && printf '%s' "$$out" | grep -qF "$$4"; then \
	       echo "  PASSED: seas_fd_cycle_driver, floor $$2, min_sigma_n_pa $$3: exit $$rc, \"$$4\""; \
	    else printf '%s\n' "$$out" | tail -15; echo "  FAILED: seas_fd_cycle_driver $$1: exit $$rc"; exit 1; fi; \
	 done; IFS=' '; \
	 [ $$n -eq 3 ] || { echo "  FAILED: $$n of the 3 deck refusals ran"; exit 1; }; \
	 deck $(SIGMA_FLOOR_TMP)/ok.toml 1.0e6 0.0; \
	 out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_fd_cycle_driver -d $(SIGMA_FLOOR_TMP)/ok.toml \
	    -o $(SIGMA_FLOOR_TMP)/ok_out --max-steps 3 2>&1); rc=$$?; \
	 if [ $$rc -eq 0 ] && printf '%s' "$$out" | grep -qF "sigma_n' floor: 1 MPa in the friction law of both regimes" \
	    && printf '%s' "$$out" | grep -qF "sigma_n' floor in the QD regime since the last switch: never acted"; then \
	    echo "  PASSED: seas_fd_cycle_driver, floor 1 MPa: accepted, reported, never acted in 3 QD steps"; \
	 else printf '%s\n' "$$out" | tail -15; echo "  FAILED: seas_fd_cycle_driver ok: exit $$rc"; exit 1; fi; \
	 rm -rf $(SIGMA_FLOOR_TMP)

# The opt-in sigma_n' relaxation ([cycle.relaxation] t_relax_yr,
# cycle/normal_stress_relaxation.hpp; Chunhui Zhao, 2026-11-02: "implement
# the relaxation, test and verify it"): K1 (the QD kernels with the
# reduction) and R1-R4 (the update rule, the stage contract, the dynamic
# view) on 1 rank, then each refusal must exit non-zero with its message;
# then the driver on the BP3-FD d60 deck (1 km mesh): the deck refusals, the
# checkpoint versions, and a stop + restart through a switch back repeating
# the run byte for byte (tests/2d_inplane/relaxation_decks.sh).
TEST_RELAXATION_SRC = tests/2d_inplane/test_relaxation.cpp
TEST_RELAXATION_OBJ = $(TEST_RELAXATION_SRC:.cpp=.o)
$(TEST_RELAXATION_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) friction/device_rate_state.hpp \
                             cycle/normal_stress_relaxation.hpp fault/fault_normal_stress_modifier.hpp \
                             $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) -DSEAS_USE_MPI -c $< -o $@

seas_test_relaxation: $(TEST_RELAXATION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_RELAXATION_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

RELAXATION_REFUSALS = 1:NormalStressRelaxation: t_relax = 0.0000000000000000e+00 s (must be > 0)|2:SetDelta: 3 values for 2 DOFs
test-relaxation: seas_test_relaxation seas_fd_cycle_driver
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_relaxation
	@list='$(RELAXATION_REFUSALS)'; set -f; IFS='|'; for r in $$list; do n=$${r%%:*}; msg=$${r#*:}; \
	   out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_relaxation --refuse $$n 2>&1); rc=$$?; \
	   if [ $$rc -ne 0 ] && printf '%s' "$$out" | grep -qF "$$msg"; then \
	      echo "  PASSED: refusal $$n: exit $$rc, \"$$msg\""; \
	   else printf '%s\n' "$$out" | tail -15; echo "  FAILED: refusal $$n: exit $$rc, \"$$msg\" not found"; exit 1; fi; \
	done
	MPIEXEC="$(MFEM_MPIEXEC)" MPIEXEC_NP="$(MFEM_MPIEXEC_NP)" tests/2d_inplane/relaxation_decks.sh

# The opt-in 2.5D finite width ([cycle.finite_width], cycle/finite_width.hpp;
# Chunhui Zhao, 2026-10-30): U1-U5 and U7 of the derivation note's Level 1 on 1 rank
# and U1, U2, U4, U5, U7 on 2 ranks; each API refusal must exit non-zero with its
# message; then the deck refusals and two accepted tables through the driver
# (tests/2d_inplane/finite_width_decks.sh, the map-view deck).
TEST_2D5_SRC = tests/2d_inplane/test_2d5_width.cpp
TEST_2D5_OBJ = $(TEST_2D5_SRC:.cpp=.o)
$(TEST_2D5_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                      $(CYCLE_SHARED_HEADERS) $(INPLANE_HEADERS) $(NITSCHE_2D_HEADERS) \
                      cycle/finite_width.hpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_2d5_width: $(TEST_2D5_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_2D5_OBJ) $(CYCLE_LINK_OBJS) \
	    $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

FW_REFUSALS = 1:SetFiniteWidth: beta / mu =|2:SetFiniteWidth: U_pl1 has|3:the 2.5D term is for 2-D map-view meshes|4:SetFiniteWidth after the dense or fused stiffness read K|5:substrate shape (QDPhase::SubstrateShape) with
test-2d5-width: seas_test_2d5_width seas_fd_cycle_driver
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_2d5_width
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_2d5_width --only U1,U2,U4,U5,U7,U8
	@list='$(FW_REFUSALS)'; set -f; IFS='|'; for r in $$list; do n=$${r%%:*}; msg=$${r#*:}; \
	   out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_2d5_width --refuse $$n 2>&1); rc=$$?; \
	   if [ $$rc -ne 0 ] && printf '%s' "$$out" | grep -qF "$$msg"; then \
	      echo "  PASSED: refusal $$n: exit $$rc, \"$$msg\""; \
	   else printf '%s\n' "$$out" | tail -15; echo "  FAILED: refusal $$n: exit $$rc, \"$$msg\" not found"; exit 1; fi; \
	done
	MPIEXEC="$(MFEM_MPIEXEC)" MPIEXEC_NP="$(MFEM_MPIEXEC_NP)" tests/2d_inplane/finite_width_decks.sh

# The Nitsche driver's [finite_width] table (the 2.5D note's Level 2 test,
# elongated2p5d/; Chunhui Zhao, 2026-10-30): the deck refusals, three
# accepted tables (c and beta in the log), and the term changing the
# stations (tests/2d_inplane/nitsche_finite_width_decks.sh, the TPV205-2D
# deck; the 3-D refusal runs when MESHES holds tpv2053d_1000m.msh).
test-2d5-dyn-driver: seas_nitsche_dyn_driver
	MPIEXEC="$(MFEM_MPIEXEC)" MPIEXEC_NP="$(MFEM_MPIEXEC_NP)" tests/2d_inplane/nitsche_finite_width_decks.sh

# Vertex-normal plan rev1, P3 (Chunhui Zhao, 2026-11-01): the cycle's two
# regimes with vertex_normals = "averaged" on the simple cases' D5f and D0
# decks (QDPhase and DynPhase, no stepping), C1-C4 on 1, 2 and 4 ranks.
TEST_VN_CYCLE_SRC = tests/2d_inplane/test_vertex_normal_cycle.cpp
TEST_VN_CYCLE_OBJ = $(TEST_VN_CYCLE_SRC:.cpp=.o)
$(TEST_VN_CYCLE_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                          $(CYCLE_SHARED_HEADERS) $(INPLANE_HEADERS) $(NITSCHE_2D_HEADERS) \
                          $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_vertex_normal_cycle: $(TEST_VN_CYCLE_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_VN_CYCLE_OBJ) $(CYCLE_LINK_OBJS) \
	    $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-vertex-normal-cycle: seas_test_vertex_normal_cycle
	@for n in 1 2 4; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_vertex_normal_cycle < /dev/null || exit 1; \
	done

# [cycle.opening] mode = "shear_free", version A of the fault opening design
# rev 1.1 (Chunhui Zhao, 2026-11-04: "start version A", "build and test
# here"; cycle/fault_opening.hpp): T1 the static mode II crack (D0), T2 open,
# re-close, a QD run and its restart, T4 the switch both ways (D9f), C1 the
# checkpoint V6 round trip, on 1 and 2 ranks; then the refusals 1-7 and a
# read of the two opening decks (T5 D9f, T6 UCERF; no solve).
TEST_FAULT_OPENING_SRC = tests/2d_inplane/test_fault_opening.cpp
TEST_FAULT_OPENING_OBJ = $(TEST_FAULT_OPENING_SRC:.cpp=.o)
$(TEST_FAULT_OPENING_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                              $(CYCLE_SHARED_HEADERS) $(INPLANE_HEADERS) $(NITSCHE_2D_HEADERS) \
                              $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_fault_opening: $(TEST_FAULT_OPENING_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_OPENING_OBJ) $(CYCLE_LINK_OBJS) \
	    $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

FAULT_OPENING_DECKS = simple_cases/configs/D9f_arc08L2p5_h250_avg_floor10_open.toml \
                      config/fd_cycle/ucerf2d_saf_sjf_500m_R70_avg_tip_relax_open_fd.toml
FAULT_OPENING_REFUSALS = 1:check_every: the solve-path check does not carry|2:use "shear_free"|3:max_passes = 0 must be >= 1|4:[cycle.opening] needs [cycle.qd.green]|5:[cycle.opening] is built and tested on 2-D meshes only|6:[cycle.opening] with [friction.prakash_clifton]|7:[cycle.probe] with [cycle.opening]|8:H-matrix path carries no opening columns|9:belong to mode = "contact"|10:floor_taper needs [cycle].sigma_n_floor_pa > 0|11:must be >= 0 (0 = the pure law, D-B7)
test-fault-opening: seas_test_fault_opening
	@for n in 1 2; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_fault_opening < /dev/null || exit 1; \
	done
	@list='$(FAULT_OPENING_REFUSALS)'; set -f; IFS='|'; fail=0; for r in $$list; do n=$${r%%:*}; msg=$${r#*:}; \
	   out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_fault_opening --refuse $$n 2>&1); rc=$$?; \
	   if [ $$rc -eq 0 ] || ! printf '%s' "$$out" | grep -qF -- "$$msg"; then \
	      echo "refusal $$n: FAIL (rc $$rc)"; printf '%s\n' "$$out" | tail -5; fail=1; \
	   else echo "refusal $$n: PASS ($$msg)"; fi; done; \
	   ./seas_test_fault_opening --clean; \
	   for d in $(FAULT_OPENING_DECKS); do \
	      out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_fault_opening --load-deck $$d 2>&1); rc=$$?; \
	      if [ $$rc -eq 0 ] && printf '%s' "$$out" | grep -qF '[cycle.opening] mode shear_free, max_passes 20; [cycle.qd.green] 1, switching 1'; then \
	         echo "deck $$d: PASS (loads, opening on)"; \
	      else echo "deck $$d: FAIL (rc $$rc)"; printf '%s\n' "$$out" | tail -5; fail=1; fi; done; \
	   exit $$fail

# Fault opening version B, B1 (plan 2026-11-05, Section 9): the opening as
# the third kernel slip component of the stage operators, operator-level
# checks on the miniature BP3 layout (g = 0 equivalence is exact; the
# opening column's sign, reciprocity and a uniform opening), on 1 and 2
# ranks.  No cycle objects: links like the layout test.
TEST_FAULT_CONTACT_SRC = tests/2d_inplane/test_fault_contact.cpp
TEST_FAULT_CONTACT_OBJ = $(TEST_FAULT_CONTACT_SRC:.cpp=.o)
$(TEST_FAULT_CONTACT_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(INPLANE_HEADERS) \
                              $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -DSEAS_USE_MPI \
	    -c $< -o $@

seas_test_fault_contact: $(TEST_FAULT_CONTACT_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_CONTACT_OBJ) \
	    $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-fault-contact: seas_test_fault_contact
	@for n in 1 2; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_fault_contact < /dev/null || exit 1; \
	done

# Fault opening version B, phase B2 (plan Section 4): the dynamic contact
# projection in CloseDOF, the kernel widening (slip_comps = 3) and the sign
# of eq. 4.3 pinned through the assembled operator.
TEST_FAULT_CONTACT_DYN_SRC = tests/2d_inplane/test_fault_contact_dyn.cpp
TEST_FAULT_CONTACT_DYN_OBJ = $(TEST_FAULT_CONTACT_DYN_SRC:.cpp=.o)
$(TEST_FAULT_CONTACT_DYN_OBJ): %.o: $(SRC)%.cpp $(SEAS_HEADERS) $(INPLANE_HEADERS) \
                              $(NITSCHE_2D_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(MFEM_FLAGS) $(SEAS_INCLUDES) $(HDF5_INCFLAGS) -DSEAS_USE_MPI \
	    -c $< -o $@

seas_test_fault_contact_dyn: $(TEST_FAULT_CONTACT_DYN_OBJ) $(HETEROGENEOUS_MATERIAL_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_CONTACT_DYN_OBJ) \
	    $(HETEROGENEOUS_MATERIAL_OBJ) $(MFEM_LIBS) $(HDF5_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-fault-contact-dyn: seas_test_fault_contact_dyn
	@for n in 1 2; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_fault_contact_dyn < /dev/null || exit 1; \
	done

# Fault opening version B, phase B3 (plan Sections 5-6): the QD contact
# solve (OpenSetQD 2x2 blocks), the active-set rule with hysteresis, the
# tapered floor and the contact deck.
TEST_FAULT_CONTACT_QD_SRC = tests/2d_inplane/test_fault_contact_qd.cpp
TEST_FAULT_CONTACT_QD_OBJ = $(TEST_FAULT_CONTACT_QD_SRC:.cpp=.o)
$(TEST_FAULT_CONTACT_QD_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                              $(CYCLE_SHARED_HEADERS) $(INPLANE_HEADERS) $(NITSCHE_2D_HEADERS) \
                              $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_fault_contact_qd: $(TEST_FAULT_CONTACT_QD_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_FAULT_CONTACT_QD_OBJ) $(CYCLE_LINK_OBJS) \
	    $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-fault-contact-qd: seas_test_fault_contact_qd
	@for n in 1 2; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_fault_contact_qd < /dev/null || exit 1; \
	done

# [cycle.qd] exact_retry (cycle/exact_retry.hpp; the UCERF quench run's QD
# crash, 2026-11-07): a PETSc RK 5dp toy, PETSc's rollback vs the saved copy.
TEST_EXACT_RETRY_SRC = tests/2d_inplane/test_exact_retry.cpp
TEST_EXACT_RETRY_OBJ = $(TEST_EXACT_RETRY_SRC:.cpp=.o)
$(TEST_EXACT_RETRY_OBJ): %.o: $(SRC)%.cpp cycle/exact_retry.hpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_exact_retry: $(TEST_EXACT_RETRY_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_EXACT_RETRY_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-exact-retry: seas_test_exact_retry
	@for n in 1 2; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_exact_retry < /dev/null || exit 1; \
	done

# Option I of the QD step at sigma_n' -> 0, phase I1 (cycle/qd_implicit.hpp;
# plan "Option I: implementation and test plan", 2026-11-09): the implicit
# fault system on small dense G's, PETSc arkimex fully implicit (T0, T1).
TEST_QD_ESDIRK_SRC = tests/2d_inplane/test_qd_esdirk_points.cpp
TEST_QD_ESDIRK_OBJ = $(TEST_QD_ESDIRK_SRC:.cpp=.o)
$(TEST_QD_ESDIRK_OBJ): %.o: $(SRC)%.cpp cycle/qd_implicit.hpp friction/device_rate_state.hpp \
                       domain/qd_green_operator.hpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_qd_esdirk_points: $(TEST_QD_ESDIRK_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_QD_ESDIRK_OBJ) $(MFEM_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-qd-esdirk-points: seas_test_qd_esdirk_points
	@for n in 1 2 4; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_qd_esdirk_points < /dev/null || exit 1; \
	done

# Option I, phase I2 (cycle/qd_implicit.hpp, PlanarProblem): the implicit
# fault system on the planar Green's function of the miniature BP3 layout
# (tractions, gap columns, Newton matrix, consistent start, one event vs
# RK 5dp, opening, all ranks vs one rank, LU time at N = 2,638).
TEST_QD_IMPLICIT_GREEN_SRC = tests/2d_inplane/test_qd_implicit_green.cpp
TEST_QD_IMPLICIT_GREEN_OBJ = $(TEST_QD_IMPLICIT_GREEN_SRC:.cpp=.o)
$(TEST_QD_IMPLICIT_GREEN_OBJ): %.o: $(SRC)%.cpp cycle/qd_implicit.hpp $(CYCLE_HEADERS) \
                               $(CYCLE_STAMP) $(CYCLE_SHARED_HEADERS) $(INPLANE_HEADERS) \
                               $(NITSCHE_2D_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_qd_implicit_green: $(TEST_QD_IMPLICIT_GREEN_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_QD_IMPLICIT_GREEN_OBJ) $(CYCLE_LINK_OBJS) \
	    $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-qd-implicit-green: seas_test_qd_implicit_green
	@for n in 1 2 4; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_qd_implicit_green < /dev/null || exit 1; \
	done
	@echo "the same suite with [cycle.qd] stage_guess = extrapolate"
	@for n in 1 2; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_qd_implicit_green \
	      -qd_stage_guess extrapolate < /dev/null || exit 1; \
	done

# cycle/par_dense_lu.hpp (ParDenseLU): the dense LU shared by all ranks
# (ScaLAPACK) against the serial LU of the implicit QD fault solve, and the
# time of both (-n N -reps R; UCERF: -n 3290).
TEST_PAR_DENSE_LU_SRC = tests/2d_inplane/test_par_dense_lu.cpp
TEST_PAR_DENSE_LU_OBJ = $(TEST_PAR_DENSE_LU_SRC:.cpp=.o)
$(TEST_PAR_DENSE_LU_OBJ): %.o: $(SRC)%.cpp cycle/par_dense_lu.hpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_par_dense_lu: $(TEST_PAR_DENSE_LU_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PAR_DENSE_LU_OBJ) $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-par-dense-lu: seas_test_par_dense_lu
	@for n in 1 2 3 4; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_par_dense_lu < /dev/null || exit 1; \
	done

# cycle/pml_partition.hpp ([cycle.pml] partition_weight): the mesh cut that
# shares the PML band's work between the ranks.
TEST_PML_PARTITION_SRC = tests/2d_inplane/test_pml_partition.cpp
TEST_PML_PARTITION_OBJ = $(TEST_PML_PARTITION_SRC:.cpp=.o)
$(TEST_PML_PARTITION_OBJ): %.o: $(SRC)%.cpp cycle/pml_partition.hpp cycle/cycle_config.hpp \
                            domain/pml_band.hpp dynamic/lts_partition.hpp $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_pml_partition: $(TEST_PML_PARTITION_OBJ) $(LTS_PARTITION_OBJ)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PML_PARTITION_OBJ) $(LTS_PARTITION_OBJ) $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-pml-partition: seas_test_pml_partition
	@for n in 1 2 4 8; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_pml_partition < /dev/null || exit 1; \
	done

# [cycle.pml] dt_rule = "element" (cycle/dyn_phase.hpp, PmlDynamicTerms::
# ElementKappa): the dynamic step from each element's own penalty factor
# instead of the whole mesh's step / sqrt(kappa_max).  The straight-fault band
# deck of the simple cases with an event started at t = 0 runs
# PML_DT_RULE_TIME s with both rules on 1 rank (the default rule first; the
# element rule takes sqrt(kappa_max) times fewer steps), and the fault
# snapshots must agree as two step sizes do (tests/2d_inplane/
# pml_dt_rule_check.py).  PML_DT_RULE_STEPS: the default rule's steps for
# that time (dt 0.00136101 s); the element rule's are found from the two
# steps in the first log.  Each run stops at the guard (exit 7) or at
# --max-dyn-steps.  The runs overwrite an earlier run's files.
PML_DT_RULE_TMP = $(or $(TMPDIR),/tmp)/seas_test_pml_dt_rule
PML_DT_RULE_DECK = simple_cases/configs/D0_straight_hb_pml0p01.toml
PML_DT_RULE_PREFIX = D0_straight_hb_pml0p01
PML_DT_RULE_STEPS = 2205
PML_DT_RULE_ESTEPS = 800
test-pml-dt-rule: seas_fd_cycle_driver
	@mkdir -p $(PML_DT_RULE_TMP); t=$(PML_DT_RULE_TMP); \
	 python3 tests/2d_inplane/pml_dt_rule_check.py decks $(PML_DT_RULE_DECK) $$t || exit 1; \
	 for v in max element; do n=$(PML_DT_RULE_STEPS); [ $$v = element ] && n=$(PML_DT_RULE_ESTEPS); \
	    $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_fd_cycle_driver -d $$t/$$v.toml -o $$t/$$v \
	       --end-phase-after $$n --max-dyn-steps $$n > $$t/$$v.log 2>&1; rc=$$?; \
	    if [ $$rc -ne 0 ] && [ $$rc -ne 7 ]; then tail -15 $$t/$$v.log; \
	       echo "  FAILED: dt_rule $$v: exit $$rc"; exit 1; fi; \
	 done; \
	 python3 tests/2d_inplane/pml_dt_rule_check.py check $$t $(PML_DT_RULE_PREFIX)

# [cycle.lts] static_split (NitscheLtsRun::SetStaticSplit; LTS with PML design
# rev1, Section 4.1, piece P-a, the rest test R1): on 1 and 2 ranks, the
# straight-fault deck of the simple cases without its band is loaded by
# LTS_SPLIT_QD QD steps (global stepping and rate2: the mesh cut differs), then
# restarted with V_on below the fault's max V, so the dynamic regime starts
# from the loaded state with the fault at rest: global stepping
# (LTS_SPLIT_STEPS steps), rate2 and rate2 with static_split = true
# (LTS_SPLIT_STEPS / 8 sync intervals: 4 clusters).  The split must keep the
# fault stresses of global stepping to 1e-7 MPa; rate2 alone must not
# (tests/2d_inplane/lts_static_split_check.py).  Each restart stops at the
# guard (exit 7) or at --max-dyn-steps.  The runs overwrite an earlier run's
# files.
LTS_SPLIT_TMP = $(or $(TMPDIR),/tmp)/seas_test_lts_static_split
LTS_SPLIT_DECK = simple_cases/configs/D0_straight_hb_pml0p01.toml
LTS_SPLIT_PREFIX = D0_straight_hb_pml0p01
LTS_SPLIT_QD = 300
LTS_SPLIT_STEPS = 1600
test-lts-static-split: seas_fd_cycle_driver
	@mkdir -p $(LTS_SPLIT_TMP); \
	 python3 tests/2d_inplane/lts_static_split_check.py decks $(LTS_SPLIT_DECK) $(LTS_SPLIT_TMP) || exit 1; \
	 for np in 1 2; do t=$(LTS_SPLIT_TMP)/np$$np; mkdir -p $$t; \
	    for v in glob lts; do \
	       $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$np ./seas_fd_cycle_driver -d $(LTS_SPLIT_TMP)/prep_$$v.toml \
	          -o $$t/$$v --max-steps $(LTS_SPLIT_QD) --checkpoint-every $(LTS_SPLIT_QD) > $$t/prep_$$v.log 2>&1 \
	          || { tail -15 $$t/prep_$$v.log; echo "  FAILED: $$np rank(s), QD loading ($$v)"; exit 1; }; \
	    done; \
	    mkdir -p $$t/split; cp $$t/lts/* $$t/split/; \
	    for v in glob lts split; do n=$$(($(LTS_SPLIT_STEPS) / 8)); [ $$v = glob ] && n=$(LTS_SPLIT_STEPS); \
	       $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$np ./seas_fd_cycle_driver -d $(LTS_SPLIT_TMP)/$$v.toml \
	          -o $$t/$$v --restart $$t/$$v/$(LTS_SPLIT_PREFIX) --end-phase-after $$n --max-dyn-steps $$n \
	          > $$t/$$v.log 2>&1; rc=$$?; \
	       if [ $$rc -ne 0 ] && [ $$rc -ne 7 ]; then tail -15 $$t/$$v.log; \
	          echo "  FAILED: $$np rank(s), restart $$v: exit $$rc"; exit 1; fi; \
	    done; \
	    grep -q "LTS rate2.*static split" $$t/split.log \
	       || { echo "  FAILED: $$np rank(s): static_split did not reach the driver"; exit 1; }; \
	    echo "  $$np rank(s):"; \
	    python3 tests/2d_inplane/lts_static_split_check.py check $$t $(LTS_SPLIT_PREFIX) || exit 1; \
	 done

# [cycle.pml] with [cycle.lts] rate2 (LTS with PML design rev1, piece P-b):
# the band's elements in one cluster, the band's memory stepped once per step
# of that cluster, the band elements' steps scaled by 1 / sqrt(kappa_e)
# (dt_rule = "element").  On 1 and 2 ranks the straight-fault band deck of the
# simple cases, with an event started at t = 0, runs LTS_PML_STEPS steps (20 s:
# the P wave is in the band from 15 s) with global stepping, with rate2 and
# the test key same_step = true (the LTS run's clusters at one step: global
# stepping to round-off) and with rate2 (LTS_PML_STEPS / 8 sync intervals).
# tests/2d_inplane/lts_pml_band_check.py compares them.  --max-steps keeps
# the QD regime after the phase short.  About 15 minutes.  The runs overwrite
# an earlier run's files.
LTS_PML_TMP = $(or $(TMPDIR),/tmp)/seas_test_lts_pml_band
LTS_PML_DECK = simple_cases/configs/D0_straight_hb_pml0p01.toml
LTS_PML_PREFIX = D0_straight_hb_pml0p01
LTS_PML_STEPS = 5336
test-lts-pml-band: seas_fd_cycle_driver
	@mkdir -p $(LTS_PML_TMP); \
	 python3 tests/2d_inplane/lts_pml_band_check.py decks $(LTS_PML_DECK) $(LTS_PML_TMP) || exit 1; \
	 for np in 1 2; do t=$(LTS_PML_TMP)/np$$np; mkdir -p $$t; \
	    for v in glob same lts; do n=$(LTS_PML_STEPS); [ $$v = lts ] && n=$$(($(LTS_PML_STEPS) / 8)); \
	       $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$np ./seas_fd_cycle_driver -d $(LTS_PML_TMP)/$$v.toml \
	          -o $$t/$$v --end-phase-after $$n --max-dyn-steps $$n --max-steps 3 > $$t/$$v.log 2>&1; rc=$$?; \
	       if [ $$rc -ne 0 ] && [ $$rc -ne 7 ]; then tail -15 $$t/$$v.log; \
	          echo "  FAILED: $$np rank(s), $$v: exit $$rc"; exit 1; fi; \
	    done; \
	    echo "  $$np rank(s):"; \
	    python3 tests/2d_inplane/lts_pml_band_check.py check $$t $(LTS_PML_PREFIX) || exit 1; \
	 done

# [cycle.quench] with [cycle.lts] rate2 (LTS with PML design rev1, piece P-c;
# NitscheLtsRun::SetQuench).  On 1 and 2 ranks the straight-fault deck of the
# simple cases without its band (the piece stands alone; LTS_QUENCH_BAND=band
# keeps the band), with an event started at t = 0 and a quench of d = 0.5 1/s
# below V_q = V_on, runs LTS_QUENCH_STEPS steps (24 s; the quench is on from
# about 19.6 s) with global stepping, with rate2 and same_step = true and
# with rate2.  tests/2d_inplane/lts_quench_check.py compares them.  About 15
# minutes.  The runs overwrite an earlier run's files.
LTS_QUENCH_TMP = $(or $(TMPDIR),/tmp)/seas_test_lts_quench
LTS_QUENCH_DECK = simple_cases/configs/D0_straight_hb_pml0p01.toml
LTS_QUENCH_PREFIX = D0_straight_hb_pml0p01
LTS_QUENCH_STEPS = 6400
LTS_QUENCH_BAND =
test-lts-quench: seas_fd_cycle_driver
	@mkdir -p $(LTS_QUENCH_TMP); \
	 python3 tests/2d_inplane/lts_quench_check.py decks $(LTS_QUENCH_DECK) $(LTS_QUENCH_TMP) $(LTS_QUENCH_BAND) || exit 1; \
	 for np in 1 2; do t=$(LTS_QUENCH_TMP)/np$$np; mkdir -p $$t; \
	    for v in glob same lts; do n=$(LTS_QUENCH_STEPS); [ $$v = lts ] && n=$$(($(LTS_QUENCH_STEPS) / 8)); \
	       $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$np ./seas_fd_cycle_driver -d $(LTS_QUENCH_TMP)/$$v.toml \
	          -o $$t/$$v --end-phase-after $$n --max-dyn-steps $$n --max-steps 3 > $$t/$$v.log 2>&1; rc=$$?; \
	       if [ $$rc -ne 0 ] && [ $$rc -ne 7 ]; then tail -15 $$t/$$v.log; \
	          echo "  FAILED: $$np rank(s), $$v: exit $$rc"; exit 1; fi; \
	    done; \
	    echo "  $$np rank(s):"; \
	    python3 tests/2d_inplane/lts_quench_check.py check $$t $(LTS_QUENCH_PREFIX) || exit 1; \
	 done

# [fault_geometry] vertex_normals = "averaged" with [cycle.lts] rate2 (LTS
# with PML design rev1, piece P-d).  On 1 and 2 ranks the D9f bend deck of the
# simple cases without its absorbing faces, with an event started at t = 0,
# runs LTS_AVG_STEPS steps (12 s: the rupture crosses the bend) with global
# stepping, with rate2 and same_step = true and with rate2.
# tests/2d_inplane/lts_avg_normals_check.py compares them.  About 10
# minutes.  The runs overwrite an earlier run's files.
LTS_AVG_TMP = $(or $(TMPDIR),/tmp)/seas_test_lts_avg_normals
LTS_AVG_DECK = simple_cases/configs/D9f_arc08L2p5_h250_avg_floor10.toml
LTS_AVG_PREFIX = D9f_arc08L2p5_h250_avg_floor10
LTS_AVG_STEPS = 6496
test-lts-avg-normals: seas_fd_cycle_driver
	@mkdir -p $(LTS_AVG_TMP); \
	 python3 tests/2d_inplane/lts_avg_normals_check.py decks $(LTS_AVG_DECK) $(LTS_AVG_TMP) || exit 1; \
	 for np in 1 2; do t=$(LTS_AVG_TMP)/np$$np; mkdir -p $$t; \
	    for v in glob same lts; do n=$(LTS_AVG_STEPS); [ $$v = lts ] && n=$$(($(LTS_AVG_STEPS) / 8)); \
	       $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$np ./seas_fd_cycle_driver -d $(LTS_AVG_TMP)/$$v.toml \
	          -o $$t/$$v --end-phase-after $$n --max-dyn-steps $$n --max-steps 3 > $$t/$$v.log 2>&1; rc=$$?; \
	       if [ $$rc -ne 0 ] && [ $$rc -ne 7 ]; then tail -15 $$t/$$v.log; \
	          echo "  FAILED: $$np rank(s), $$v: exit $$rc"; exit 1; fi; \
	    done; \
	    echo "  $$np rank(s):"; \
	    python3 tests/2d_inplane/lts_avg_normals_check.py check $$t $(LTS_AVG_PREFIX) || exit 1; \
	 done

# [cycle.opening] mode = "contact" with [cycle.lts] rate2 (LTS with PML design
# rev1).  On 1 and 2 ranks the D9f bend deck as in test-lts-avg-normals, with
# contact and numbers chosen so that one event opens the releasing side of the
# bend, runs LTS_CONTACT_STEPS steps (14.8 s) with global stepping, with rate2
# and same_step = true and with rate2.  The phase is ended in mid event, so the
# hand-back guard stops each run (exit 7): that is expected here.
# tests/2d_inplane/lts_contact_check.py compares them.  About 12 minutes.  The
# runs overwrite an earlier run's files.
LTS_CONTACT_TMP = $(or $(TMPDIR),/tmp)/seas_test_lts_contact
LTS_CONTACT_STEPS = 8000
test-lts-contact: seas_fd_cycle_driver
	@mkdir -p $(LTS_CONTACT_TMP); \
	 python3 tests/2d_inplane/lts_contact_check.py decks $(LTS_AVG_DECK) $(LTS_CONTACT_TMP) || exit 1; \
	 for np in 1 2; do t=$(LTS_CONTACT_TMP)/np$$np; mkdir -p $$t; \
	    for v in glob same lts; do n=$(LTS_CONTACT_STEPS); [ $$v = lts ] && n=$$(($(LTS_CONTACT_STEPS) / 8)); \
	       $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$np ./seas_fd_cycle_driver -d $(LTS_CONTACT_TMP)/$$v.toml \
	          -o $$t/$$v --end-phase-after $$n --max-dyn-steps $$n --max-steps 3 > $$t/$$v.log 2>&1; rc=$$?; \
	       if [ $$rc -ne 0 ] && [ $$rc -ne 7 ]; then tail -15 $$t/$$v.log; \
	          echo "  FAILED: $$np rank(s), $$v: exit $$rc"; exit 1; fi; \
	    done; \
	    echo "  $$np rank(s):"; \
	    python3 tests/2d_inplane/lts_contact_check.py check $$t $(LTS_AVG_PREFIX) || exit 1; \
	 done

# [cycle.finite_width] with [cycle.lts] rate2 (2.5D finite width with PML,
# quench and LTS, plan rev1, piece F-c).  On 1 and 2 ranks the straight-fault
# deck of the simple cases without its band (the piece stands alone;
# LTS_2D5_BAND=band keeps the band: pieces F-a and F-c together), with an
# event started at t = 0 and W = 15 km, runs LTS_2D5_STEPS steps (24 s) with
# global stepping, with rate2 and same_step = true, with rate2, and with
# global stepping without the 2.5D term (the control).
# tests/2d_inplane/lts_2d5_check.py compares them.  About 20 minutes.  The
# runs overwrite an earlier run's files.
LTS_2D5_TMP = $(or $(TMPDIR),/tmp)/seas_test_lts_2d5
LTS_2D5_STEPS = 6400
LTS_2D5_BAND =
test-lts-2d5: seas_fd_cycle_driver
	@mkdir -p $(LTS_2D5_TMP); \
	 python3 tests/2d_inplane/lts_2d5_check.py decks $(LTS_QUENCH_DECK) $(LTS_2D5_TMP) $(LTS_2D5_BAND) || exit 1; \
	 for np in 1 2; do t=$(LTS_2D5_TMP)/np$$np; mkdir -p $$t; \
	    for v in off glob same lts; do n=$(LTS_2D5_STEPS); [ $$v = lts ] && n=$$(($(LTS_2D5_STEPS) / 8)); \
	       $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$np ./seas_fd_cycle_driver -d $(LTS_2D5_TMP)/$$v.toml \
	          -o $$t/$$v --end-phase-after $$n --max-dyn-steps $$n --max-steps 3 > $$t/$$v.log 2>&1; rc=$$?; \
	       if [ $$rc -ne 0 ] && [ $$rc -ne 7 ]; then tail -15 $$t/$$v.log; \
	          echo "  FAILED: $$np rank(s), $$v: exit $$rc"; exit 1; fi; \
	    done; \
	    echo "  $$np rank(s):"; \
	    python3 tests/2d_inplane/lts_2d5_check.py check $$t $(LTS_QUENCH_PREFIX) || exit 1; \
	 done

# fd_cycle --early-sync (NitscheWaveOperator::SetEarlySync, DenseBlockStiffness::
# ExchangeHalo; speed-up options D0, 2026-11-11): the stiffness halo at the
# start of a dynamic step and the fault record's reduction at its end.  It
# moves MPI calls only, so the map-view tensor deck with its event started at
# t = 0 (the deck of test-handback-table), 200 dynamic steps, must write the
# same fault table, global file and on-fault station files byte for byte with
# the flag on and off, on 1, 2 and 4 ranks (each run stops at the guard, exit
# 7, as in test-handback-table).  The runs overwrite an earlier run's files.
EARLY_SYNC_TMP = $(or $(TMPDIR),/tmp)/seas_test_early_sync
test-early-sync: seas_fd_cycle_driver
	@mkdir -p $(EARLY_SYNC_TMP); \
	 python3 tests/2d_inplane/handback_table_check.py decks $(HANDBACK_DECK) $(EARLY_SYNC_TMP) || exit 1; \
	 for np in 1 2 4; do \
	    for v in off on; do fl=; [ $$v = on ] && fl=--early-sync; d=$(EARLY_SYNC_TMP)/np$${np}_$$v; \
	       $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$np ./seas_fd_cycle_driver -d $(EARLY_SYNC_TMP)/on.toml \
	          -o $$d --end-phase-after $(HANDBACK_STEPS) $$fl > $$d.log 2>&1; rc=$$?; \
	       if [ $$rc -ne 0 ] && [ $$rc -ne 7 ]; then tail -15 $$d.log; \
	          echo "  FAILED: $$np rank(s), early sync $$v: exit $$rc"; exit 1; fi; \
	    done; \
	    grep -q -- "--early-sync: stiffness halo" $(EARLY_SYNC_TMP)/np$${np}_on.log \
	       || { echo "  FAILED: $$np rank(s): the flag did not reach the operator"; exit 1; }; \
	    n=0; for f in $(EARLY_SYNC_TMP)/np$${np}_off/*_fault.txt $(EARLY_SYNC_TMP)/np$${np}_off/*_global.txt \
	                  $(EARLY_SYNC_TMP)/np$${np}_off/*_fltst_*.txt; do \
	       cmp -s $$f $(EARLY_SYNC_TMP)/np$${np}_on/$$(basename $$f) \
	          || { echo "  FAILED: $$np rank(s): $$(basename $$f) differs"; exit 1; }; n=$$((n + 1)); \
	    done; \
	    [ $$n -ge 3 ] || { echo "  FAILED: $$np rank(s): only $$n files compared"; exit 1; }; \
	    echo "  PASSED: $$np rank(s): $$n fault files byte-identical with --early-sync on and off"; \
	 done

# Option I, phase I3: the implicit QD fault solve in the cycle driver
# ([cycle.qd] integrator = "esdirk"; cycle/qd_phase.hpp).  The D9f contact
# deck with the key on, 80 QD steps: a run stopped at step 40 and restarted
# from its checkpoint repeats the uninterrupted run bit for bit, on 1 and 2
# ranks (tests/cycle/test_cycle_restart.sh; about 2 minutes per rank count).
# On 2 ranks the parallel direct solver builds a G whose last digits change
# from run to run, so there the three runs share one G file (GREEN_FILE=1):
# every fault, global, source and switch row is then identical, and the
# off-fault stations and the checkpoint's displacement, which come from the
# parallel solve itself, agree to VOLUME_TOL.
QD_ESDIRK_CYCLE_TMP ?= /tmp/seas_qd_esdirk_cycle
test-qd-esdirk-cycle: seas_fd_cycle_driver
	@for n in 1 2; do \
	   MPIRUN="$(MFEM_MPIEXEC)" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
	   GREEN_FILE=$$((n > 1)) VOLUME_TOL=1e-6 \
	   bash tests/cycle/test_cycle_restart.sh \
	      simple_cases/configs/D9f_arc08L2p5_h250_avg_floor10_contact_esdirk.toml \
	      $(QD_ESDIRK_CYCLE_TMP)_np$$n 40 $$n || exit 1; \
	done

# [cycle.relaxation] opening = "relaxed" with option I (Chunhui Zhao,
# 2026-11-13; cycle/normal_stress_relaxation.hpp): the deck refusals and the
# accepted deck (tests/2d_inplane/relaxed_opening_decks.sh), then the restart
# test of test-qd-esdirk-cycle on the D9f deck with the relaxation on (80 QD
# steps, stopped at 40 and restarted: bit for bit, 1 and 2 ranks).  The
# point test is T1e of test-qd-esdirk-points, the dynamic contact rule D4 of
# test-fault-contact-dyn.
QD_ESDIRK_RELAX_TMP ?= /tmp/seas_qd_esdirk_relax
test-qd-esdirk-relax: seas_fd_cycle_driver
	MPIEXEC="$(MFEM_MPIEXEC)" MPIEXEC_NP="$(MFEM_MPIEXEC_NP)" tests/2d_inplane/relaxed_opening_decks.sh
	@for n in 1 2; do \
	   MPIRUN="$(MFEM_MPIEXEC)" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
	   GREEN_FILE=$$((n > 1)) VOLUME_TOL=1e-6 \
	   bash tests/cycle/test_cycle_restart.sh \
	      simple_cases/configs/D9f_arc08L2p5_h250_avg_floor10_contact_esdirk_relax5k.toml \
	      $(QD_ESDIRK_RELAX_TMP)_np$$n 40 $$n || exit 1; \
	done

2d-inplane-all: test-fault-frame-field test-plate-loading test-okada-dislocation test-2d-fault-layout \
                test-2d-deck test-2d-nitsche test-2d-driver test-thermal-pressurization \
                test-normal-stress-hook test-2d-driver-tp test-2d-dislocation test-prestress-bp3 \
                test-2d-absorbing test-absorbing-all test-snapshots-volume test-tensor-prestress \
                test-handback-table test-tip-zone test-sigma-floor test-relaxation test-2d5-width test-2d5-dyn-driver \
                test-vertex-normals test-vertex-normal-kernels test-vertex-normal-cycle \
                test-fault-opening test-fault-contact test-fault-contact-dyn \
                test-fault-contact-qd test-exact-retry test-qd-esdirk-points test-qd-implicit-green \
                test-qd-esdirk-cycle test-qd-esdirk-relax test-par-dense-lu test-pml-partition test-early-sync test-lts-static-split test-pml-dt-rule test-lts-pml-band test-lts-quench test-lts-avg-normals test-lts-contact test-lts-2d5

.PHONY: test-fault-frame-field test-plate-loading test-okada-dislocation test-2d-fault-layout \
        test-2d-deck test-2d-nitsche test-2d-driver test-thermal-pressurization \
        test-normal-stress-hook test-2d-driver-tp test-2d-dislocation test-prestress-bp3 \
        test-2d-absorbing test-absorbing-all test-snapshots-volume test-tensor-prestress \
        test-handback-table test-tip-zone test-sigma-floor test-relaxation test-2d5-width test-2d5-dyn-driver \
        test-vertex-normals test-vertex-normal-kernels test-vertex-normal-cycle \
        test-fault-opening test-exact-retry test-qd-esdirk-points test-qd-implicit-green \
        test-qd-esdirk-cycle test-qd-esdirk-relax test-par-dense-lu test-pml-partition 2d-inplane-all

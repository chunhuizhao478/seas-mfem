# tests/ucerf3d/ucerf3d.mk - build file of the UCERF 3-D branch's tests
# (fully-dynamic-fault-opening-3d, plan rev5 Section 10).  It includes the PML
# test file (which includes the 2D in-plane file, cycle/cycle.mk and the
# shared Makefile unchanged), so every target there works through it too.
#
# Run from miniapps/seas:
#   make -f tests/ucerf3d/ucerf3d.mk test-qd-pcg          # Q1.1: U-Q1 (p = 1, 1 rank, amg and aux_h1) and U-Q2
#   make -f tests/ucerf3d/ucerf3d.mk test-vertex-normals-3d  # P1.1: U-N1, U-N2 (3-D vertex normals)
#   make -f tests/ucerf3d/ucerf3d.mk test-a-field-3d      # P1.2: U-A1 (3-D tip and depth a-field)
#   make -f tests/ucerf3d/ucerf3d.mk test-qd-pcg UQ1_ORDERS="1 2" UQ1_NPS="1 4"
#       (p = 2: MUMPS ran out of memory, error -13, on the 14 GB cloud host; a cluster node)

include tests/pml/pml.mk

TEST_QD_PCG_SRC = tests/ucerf3d/test_qd_pcg.cpp
TEST_QD_PCG_OBJ = $(TEST_QD_PCG_SRC:.cpp=.o)
$(TEST_QD_PCG_OBJ): %.o: $(SRC)%.cpp solver/qd_linear_solver.hpp $(CYCLE_HEADERS) \
                    $(CYCLE_STAMP) $(CYCLE_SHARED_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_qd_pcg: $(TEST_QD_PCG_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_QD_PCG_OBJ) $(CYCLE_LINK_OBJS) $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

UQ1_ORDERS ?= 1
UQ1_NPS ?= 1
UQ2_REFUSALS = 1:pcg_rel_tol needs [cycle.qd].solver = "pcg"|2:pcg_rel_tol must lie in (0, 1)|3:pcg_max_it must be >= 1|4:pcg_amg_strength must lie in (0, 1)|5:pcg_smoother must be "default" or "l1_jacobi"|6:solver must be "mumps", "mumps_blr" or "pcg"|7:function is a direct-solver build; use one or the other|8:needs a symmetric positive definite K|9:QD PCG: no convergence in 2 iterations|10:pcg_preconditioner must be "amg" or "aux_h1"|11:pcg_amg_strength needs [cycle.qd].pcg_preconditioner
test-qd-pcg: seas_test_qd_pcg
	@for p in $(UQ1_ORDERS); do for n in $(UQ1_NPS); do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_qd_pcg --order $$p || exit 1; \
	done; done
	@list='$(UQ2_REFUSALS)'; set -f; IFS='|'; for r in $$list; do n=$${r%%:*}; msg=$${r#*:}; \
	   out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_qd_pcg --refuse $$n 2>&1); rc=$$?; \
	   if [ $$rc -ne 0 ] && printf '%s' "$$out" | grep -qF "$$msg"; then \
	      echo "  PASSED: U-Q2 refusal $$n: exit $$rc, \"$$msg\""; \
	   else printf '%s\n' "$$out" | tail -15; echo "  FAILED: U-Q2 refusal $$n: exit $$rc, \"$$msg\" not found"; exit 1; fi; \
	done

.PHONY: test-qd-pcg

# P1.1 (plan rev5 U-N1, U-N2; Chunhui Zhao, 2026-11-25: "start P1"): the
# averaged vertex normals on 3-D fault surfaces on 1, 2 and 4 ranks, the
# VN3-DIGEST lines the same on every rank count, then the refusals.
TEST_VN3D_SRC = tests/ucerf3d/test_vertex_normals_3d.cpp
TEST_VN3D_OBJ = $(TEST_VN3D_SRC:.cpp=.o)
$(TEST_VN3D_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                  $(CYCLE_SHARED_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_vertex_normals_3d: $(TEST_VN3D_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_VN3D_OBJ) $(CYCLE_LINK_OBJS) $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

VN3D_TMP = $(or $(TMPDIR),/tmp)/seas_test_vertex_normals_3d
VN3D_REFUSALS = 1:junction|2:must be > 0
test-vertex-normals-3d: seas_test_vertex_normals_3d
	@mkdir -p $(VN3D_TMP)
	@set -e; for n in 1 2 4; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_vertex_normals_3d > $(VN3D_TMP)/np$$n.log 2>&1 \
	      || { cat $(VN3D_TMP)/np$$n.log; echo "  FAILED: $$n rank(s)"; exit 1; }; \
	   grep -E "PASSED|FAILED|tests passed|^  +[0-9]+ +[0-9.]+ " $(VN3D_TMP)/np$$n.log; \
	   grep "^VN3-DIGEST" $(VN3D_TMP)/np$$n.log > $(VN3D_TMP)/digest_np$$n.txt; done; \
	 d=$(VN3D_TMP)/digest_np; \
	 if cmp -s $${d}1.txt $${d}2.txt && cmp -s $${d}1.txt $${d}4.txt && [ "$$(wc -l < $${d}1.txt)" -eq 3 ]; \
	 then echo "  PASSED: U-N1d the vertex-normal tables are the same bit for bit on 1, 2 and 4 ranks"; \
	 else echo "  FAILED: U-N1d the VN3-DIGEST lines differ between rank counts"; exit 1; fi
	@list='$(VN3D_REFUSALS)'; set -f; IFS='|'; for r in $$list; do n=$${r%%:*}; msg=$${r#*:}; \
	   out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_vertex_normals_3d --refuse $$n 2>&1); rc=$$?; \
	   if [ $$rc -ne 0 ] && printf '%s' "$$out" | grep -qF "$$msg"; then \
	      echo "  PASSED: refusal $$n: exit $$rc, \"$$msg\""; \
	   else printf '%s\n' "$$out" | tail -15; echo "  FAILED: refusal $$n: exit $$rc, \"$$msg\" not found"; exit 1; fi; \
	done

.PHONY: test-vertex-normals-3d

# P1.2 (plan rev5 U-A1; Chunhui Zhao, 2026-11-25: "start P1"): the tip and
# depth a-field on 3-D fault surfaces on 1, 2 and 4 ranks, the AF3-DIGEST
# lines the same on every rank count, then the refusals.
TEST_AF3D_SRC = tests/ucerf3d/test_a_field_3d.cpp
TEST_AF3D_OBJ = $(TEST_AF3D_SRC:.cpp=.o)
$(TEST_AF3D_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) cycle/fault_surface_zones.hpp $(CYCLE_STAMP) \
                  $(CYCLE_SHARED_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_a_field_3d: $(TEST_AF3D_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_AF3D_OBJ) $(CYCLE_LINK_OBJS) $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

AF3D_TMP = $(or $(TMPDIR),/tmp)/seas_test_a_field_3d
AF3D_REFUSALS = 1:three or more fault faces|2:lies on no fault triangle
test-a-field-3d: seas_test_a_field_3d
	@mkdir -p $(AF3D_TMP)
	@set -e; for n in 1 2 4; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_a_field_3d > $(AF3D_TMP)/np$$n.log 2>&1 \
	      || { cat $(AF3D_TMP)/np$$n.log; echo "  FAILED: $$n rank(s)"; exit 1; }; \
	   grep -E "PASSED|FAILED|tests passed" $(AF3D_TMP)/np$$n.log; \
	   grep "^AF3-DIGEST" $(AF3D_TMP)/np$$n.log > $(AF3D_TMP)/digest_np$$n.txt; done; \
	 d=$(AF3D_TMP)/digest_np; \
	 if cmp -s $${d}1.txt $${d}2.txt && cmp -s $${d}1.txt $${d}4.txt && [ "$$(wc -l < $${d}1.txt)" -eq 4 ]; \
	 then echo "  PASSED: U-A1d the a samples are the same bit for bit on 1, 2 and 4 ranks"; \
	 else echo "  FAILED: U-A1d the AF3-DIGEST lines differ between rank counts"; exit 1; fi
	@list='$(AF3D_REFUSALS)'; set -f; IFS='|'; for r in $$list; do n=$${r%%:*}; msg=$${r#*:}; \
	   out=$$($(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_a_field_3d --refuse $$n 2>&1); rc=$$?; \
	   if [ $$rc -ne 0 ] && printf '%s' "$$out" | grep -qF "$$msg"; then \
	      echo "  PASSED: refusal $$n: exit $$rc, \"$$msg\""; \
	   else printf '%s\n' "$$out" | tail -15; echo "  FAILED: refusal $$n: exit $$rc, \"$$msg\" not found"; exit 1; fi; \
	done

.PHONY: test-a-field-3d

# P1.3 (plan rev5 U-C1 to U-C3; Chunhui Zhao, 2026-11-26: "proceed next"):
# unilateral contact in the 3-D Nitsche dynamic phase, on 1 and 2 ranks.
TEST_CT3D_SRC = tests/ucerf3d/test_contact_3d.cpp
TEST_CT3D_OBJ = $(TEST_CT3D_SRC:.cpp=.o)
$(TEST_CT3D_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) domain/nitsche_energy.hpp $(CYCLE_STAMP) \
                  $(CYCLE_SHARED_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_contact_3d: $(TEST_CT3D_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_CT3D_OBJ) $(CYCLE_LINK_OBJS) $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-contact-3d: seas_test_contact_3d
	@for n in 1 2; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_contact_3d || exit 1; \
	done

.PHONY: test-contact-3d

# P1.9 (plan rev5 U-G3): the device step's contact branch against the host
# (U-G3a U-C1's block, U-G3b U-C3's T-bend run, U-G3c that run with the
# relaxed contact condition, U-G3d with a reduction, a floor and the floor
# taper), debug device, 1 and 2 ranks.
test-device-contact-3d: seas_test_contact_3d
	@for n in 1 2; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_contact_3d --device cpu \
	      --write test_device_contact_3d_np$$n || exit 1; \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_contact_3d --device debug \
	      --compare test_device_contact_3d_np$$n || exit 1; \
	   rm -f test_device_contact_3d_np$$n.U-G3*; \
	done

.PHONY: test-device-contact-3d

# P1.4 (plan rev5 U-P1, U-P2): the CFS-PML band and local time stepping in
# 3-D.  U-P1 (--only P1C): a plane P wave down a uniform column, band on the
# two end walls, against the band's design law (about 8 min on 2 ranks).
# U-P2 (--only P2): an event on T-pml to 14 s with global stepping, rate-2
# clusters at one step and rate-2 clusters (about 12 min on 4 ranks).
TEST_PL3D_SRC = tests/ucerf3d/test_pml_lts_3d.cpp
TEST_PL3D_OBJ = $(TEST_PL3D_SRC:.cpp=.o)
$(TEST_PL3D_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                  $(CYCLE_SHARED_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_pml_lts_3d: $(TEST_PL3D_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_PL3D_OBJ) $(CYCLE_LINK_OBJS) $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-pml-lts-3d: seas_test_pml_lts_3d
	OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 2 ./seas_test_pml_lts_3d --only P1C
	OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 4 ./seas_test_pml_lts_3d --only P2 --sync 1480

.PHONY: test-pml-lts-3d

# P1.5 (plan rev5 U-S1): the constant stress tensor on a dipping fault
# (T-dip), on 1, 2 and 4 ranks.
TEST_TP3D_SRC = tests/ucerf3d/test_tensor_prestress_3d.cpp
TEST_TP3D_OBJ = $(TEST_TP3D_SRC:.cpp=.o)
$(TEST_TP3D_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                  $(CYCLE_SHARED_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_tensor_prestress_3d: $(TEST_TP3D_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_TP3D_OBJ) $(CYCLE_LINK_OBJS) $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-tensor-prestress-3d: seas_test_tensor_prestress_3d
	@for n in 1 2 4; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_tensor_prestress_3d || exit 1; \
	done

.PHONY: test-tensor-prestress-3d

# P1.6 (plan rev5 U-L1, U-L2): loading options A (walls and creep faces
# off the plane) and B (walls only) on T-bend: all checks on 1 rank; on 2
# and 4 ranks without the two static solves (the second direct solve on
# several ranks was killed for memory on the 15 GB test machine); then with
# the mesh cut along the fault and creep planes on 2 and 4 ranks, so that
# the creep faces are shared.
TEST_LD3D_SRC = tests/ucerf3d/test_loading_3d.cpp
TEST_LD3D_OBJ = $(TEST_LD3D_SRC:.cpp=.o)
$(TEST_LD3D_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                  $(CYCLE_SHARED_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_loading_3d: $(TEST_LD3D_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_LD3D_OBJ) $(CYCLE_LINK_OBJS) $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-loading-3d: seas_test_loading_3d
	$(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) 1 ./seas_test_loading_3d
	@for n in 2 4; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_loading_3d --no-solve || exit 1; \
	done
	@for n in 2 4; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_loading_3d --split --no-solve || exit 1; \
	done

.PHONY: test-loading-3d

# P1.7 (plan rev5 U-G1): the cycle's dynamic phase on the mfem device
# ([cycle].device_step) against the host: a switch-in, 300 steps and the
# hand-back state on T-plane, debug device, 1 and 2 ranks.
TEST_DC3D_SRC = tests/ucerf3d/test_device_cycle_3d.cpp
TEST_DC3D_OBJ = $(TEST_DC3D_SRC:.cpp=.o)
$(TEST_DC3D_OBJ): %.o: $(SRC)%.cpp $(CYCLE_HEADERS) $(CYCLE_STAMP) \
                  $(CYCLE_SHARED_HEADERS) $(MFEM_LIB_FILE) $(CONFIG_MK)
	@mkdir -p $(@D)
	$(MFEM_CXX) $(CYCLE_CXXFLAGS) -c $< -o $@

seas_test_device_cycle_3d: $(TEST_DC3D_OBJ) $(CYCLE_LINK_OBJS)
	$(MFEM_CXX) $(MFEM_LINK_FLAGS) -o $@ $(TEST_DC3D_OBJ) $(CYCLE_LINK_OBJS) $(CYCLE_LIBS)
	$(SEAS_POST_LINK_DEDUP_RPATH)

test-device-cycle-3d: seas_test_device_cycle_3d
	@for n in 1 2; do \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_device_cycle_3d --device cpu \
	      --write test_device_cycle_3d_np$$n.bin || exit 1; \
	   $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_device_cycle_3d --device debug \
	      --compare test_device_cycle_3d_np$$n.bin || exit 1; \
	   rm -f test_device_cycle_3d_np$$n.bin; \
	done

.PHONY: test-device-cycle-3d

# P1.8 (plan rev5 U-G2): the device step's quench, relaxation r and floor
# against the host, debug device, 1 and 2 ranks.
#   a  one step: quench on, r below (x < 0) and above (x >= 0) the floor
#   b  300 steps: the same, quench on after step 100 and off after step 200
#   c  300 steps: quench on from the start, r without a floor
TEST_DQ3D_A = --name U-G2a --steps 1 --min-peak 0 --quench-from 0 --relax --floor 24.5e6
TEST_DQ3D_B = --name U-G2b --steps 300 --quench-from 100 --quench-to 200 --relax --floor 24.5e6
TEST_DQ3D_C = --name U-G2c --steps 300 --quench-from 0 --relax
test-device-quench-3d: seas_test_device_cycle_3d
	@for n in 1 2; do \
	   for c in "$(TEST_DQ3D_A)" "$(TEST_DQ3D_B)" "$(TEST_DQ3D_C)"; do \
	      $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_device_cycle_3d $$c --device cpu \
	         --write test_device_quench_3d_np$$n.bin || exit 1; \
	      $(MFEM_MPIEXEC) $(MFEM_MPIEXEC_NP) $$n ./seas_test_device_cycle_3d $$c --device debug \
	         --compare test_device_quench_3d_np$$n.bin || exit 1; \
	      rm -f test_device_quench_3d_np$$n.bin; \
	   done; \
	done

.PHONY: test-device-quench-3d

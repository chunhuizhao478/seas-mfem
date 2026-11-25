# tests/ucerf3d/ucerf3d.mk - build file of the UCERF 3-D branch's tests
# (fully-dynamic-fault-opening-3d, plan rev5 Section 10).  It includes the PML
# test file (which includes the 2D in-plane file, cycle/cycle.mk and the
# shared Makefile unchanged), so every target there works through it too.
#
# Run from miniapps/seas:
#   make -f tests/ucerf3d/ucerf3d.mk test-qd-pcg          # Q1.1: U-Q1 (p = 1, 1 rank, amg and aux_h1) and U-Q2
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

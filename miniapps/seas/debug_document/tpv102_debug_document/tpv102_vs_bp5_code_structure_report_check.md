# Code Review: tpv102_vs_bp5_code_structure_report_2026-04-22.md

(Placed in `miniapps/seas/debug_document/tpv102_debug_document/` rather than
project root to avoid clobbering unrelated existing `REVIEW.md` files at
`/Users/chunhuizhao/projects/seas-mfem/REVIEW.md` and
`miniapps/seas/REVIEW.md`. This file is the `_check.md` companion to the
review subject and follows the existing directory convention.)

## Review Scope

- **Subject under review**: `tpv102_vs_bp5_code_structure_report_2026-04-22.md`
  (the proposal to reuse BP5 abstractions in the TPV102 dynamic path).
- **Review charter** (from `/code-review` invocation): "make sure your
  purpose of reuse CAN NOT modify existing BP5 functionality or cause break
  of BP5 run, if additional code is needed, then TPV102 should just create
  a separate file."
- **Authoritative rule** (`CLAUDE.md`, confirmed by user-memory
  `feedback_tpv102_bp5_no_shared_edit.md`): duplicate the BP5 utility into a
  new file under `dynamic/` with identical signature; do not edit BP5
  source; this eases the future quasi-dynamic + dynamic merge.
- **Files inspected to verify claims in the report**:
  - `miniapps/seas/common/fault_scatter.hpp` (FaultScatter / SharedFaultCommBlock)
  - `miniapps/seas/fault/face_quadrature.hpp` (FaceQuadrature)
  - `miniapps/seas/solver/time_stepper.hpp` (DormandPrinceRK45)
  - `miniapps/seas/fault/fault_basis.hpp` (FaultBasis / FaultBasisQPData)
  - `miniapps/seas/dynamic/wave_operator.{hpp,inl}` (fault-face flux call sites)
  - `miniapps/seas/drivers/tpv102_driver.cpp` (driver loop, nucleation)
  - `miniapps/seas/drivers/seas_driver.cpp` (BP5 reference loop)

## Findings

### [R-001] [CRITICAL] [report §B2, §C row 5] — FaceQuadrature reuse for nucleation requires modifying shared `fault/face_quadrature.hpp`

**Category:** DEVIATION (violates the no-BP5-modify rule)

**Description:**
Report §B2 proposes using `FaceQuadrature` (from `fault/face_quadrature.hpp`,
shared with BP5) as the fault-side L2 space for nucleation projection. The
report does not flag that `FaceQuadrature`'s quadrature rule is hard-wired
to `IntRules.Get(face_geom, 2 * vol_order + 1)` (face_quadrature.hpp:69-73),
which differs from the TPV102 wave operator's fault-face rule
`IntRules.Get(ftr->GetGeometryType(), 2 * order_)` (wave_operator.inl:428,
1823, 2060, etc.). The two rules are numerically distinct for every
practical order (e.g., `order=2` → 6 vs 12 QPs on a triangle); projecting
nucleation with one rule while evaluating flux at the other introduces
an aliasing term at every fault QP.

**Trigger:**
Any attempt to use `FaceQuadrature::ProjectQPToDOF` or the reference mass
matrix inverse from `FaceQuadrature` as the receptacle for nucleation while
leaving the wave operator's flux-integration rule unchanged.

**Actual behavior (of the report as written):**
Silently assumes that `FaceQuadrature`'s QP layout matches the wave
operator's QP layout. It doesn't. The reader of the report would land in
one of two failure modes:
  (a) Modify `FaceQuadrature` to accept an explicit quadrature order
      parameter — which edits a file shared with BP5 and violates the rule.
  (b) Change the wave operator's fault-face rule to match `FaceQuadrature`'s —
      which changes the accepted BP5 quadrature order on fault faces (BP5
      also consumes the wave operator's rule in shared unit tests such as
      `seas_test_face_quadrature`, `seas_test_fault_face_flux_frame_and_flux`).

**Expected behavior:**
Either of:
  - Explicitly state that the nucleation projection must live in a new
    `dynamic/dynamic_face_quadrature.hpp` (duplicating the small number of
    helpers needed — reference mass matrix, `e_q_`, `ir_` — with the
    `2 * order_` rule instead of `2 * vol_order + 1`).
  - Or state that nucleation will not be re-projected via FaceQuadrature
    and will remain QP-local (using the wave operator's rule directly),
    and list which FaceQuadrature-provided utility the dynamic path will
    duplicate.

**Suggested fix:** in the report, Part B §B2 and Part C row 5, add an
explicit "Duplicate into `dynamic/dynamic_face_quadrature.hpp` with
`2*order_` rule" paragraph, or replace the FaceQuadrature reuse with a
simpler "QP-local projection with the existing wave-operator rule."

```diff
 ### B2. `FaceQuadrature` as the fault-side L2 space
-`fault/face_quadrature.hpp` gives a proper H1 face basis + reference mass
-matrix inverse, already validated at p=0..4 on triangles and in use by
+`fault/face_quadrature.hpp` is shared with BP5 and hard-wires the
+quadrature rule to `2*vol_order+1`. The wave operator integrates fault
+faces at `2*order_`. To avoid modifying the BP5 header, duplicate the
+required pieces (e_q_, M_ref_inv_, ir_) into a new
+`dynamic/dynamic_face_quadrature.hpp` parameterised by the wave-operator
+rule. Do NOT add a "quadrature-order" parameter to the shared file.
+
+The duplicate gives a proper H1 face basis + reference mass matrix
+inverse, already validated at p=0..4 on triangles via the same
+pattern used by
 `elasticity_operator_traction.inl` for projecting QP traction to fault
 DOFs. Nucleation is by far the most natural use:
```

**Test case:**
```python
def test_R001_facequadrature_rule_mismatch():
    """The shared fault QP count from wave_operator.inl's 2*order rule
    must match the QP count used by any nucleation-L2 code path."""
    # For order=2 on a triangle:
    # IntRules.Get(TRIANGLE, 2*2)   -> 6 QPs  (wave operator rule)
    # IntRules.Get(TRIANGLE, 2*2+1) -> 6 QPs
    # For order=3:
    # IntRules.Get(TRIANGLE, 2*3)   -> 7 QPs
    # IntRules.Get(TRIANGLE, 2*3+1) -> 12 QPs   ← mismatch
    # Assert equal QP count for orders 1..5 between the rule the wave
    # operator actually uses (2*order) and whatever nucleation uses.
    for p in [1, 2, 3, 4, 5]:
        assert n_qp(2*p) == n_qp_nucleation(p), \
            f"order {p}: rule mismatch -> aliasing at every fault QP"
```

---

### [R-002] [CRITICAL] [report §B6, §C row 6] — DormandPrinceRK45 is not state-layout-agnostic; reusing it for TPV102 requires edits to `solver/time_stepper.hpp`

**Category:** DEVIATION (violates the no-BP5-modify rule)

**Description:**
Report §B6 says `DormandPrinceRK45` is a "drop-in replacement" for TPV102's
hand-rolled RK4. It is not. `time_stepper.hpp:629` hard-codes
`int state_per_node_ = 2`, and the stepper's V-guard logic at lines
253-269, 330, 354, 379, 404 interprets the state vector as
`[slip_c0, slip_c1, ..., slip_c{spn-2}, psi_i]` per node:

```cpp
int spn = state_per_node_;
int n_dofs = stage_k.Size() / spn;
for (int c = 0; c < spn - 1; c++)  // "slip rate components (not psi)"
```

For TPV102, the ODE state would be `[Q_bulk_{NUM_STATE*ndof_total}; slip_dip;
slip_strike; psi_fault]` — a composite vector with NO uniform "per-node"
layout. Feeding this into `DormandPrinceRK45` would:
  - Divide `stage_k.Size()` by `state_per_node_=3` and interpret each
    tuple of 3 floats as a fault-node (`slip_dip, slip_strike, psi`).
  - Scan wave-velocity components (in `Q_bulk[VX..VZ]` — floats on the
    order of `c_s ~ 3000 m/s` during rupture) with the V-guard, which
    triggers at `v_max > v_guard_factor * v_max_stage0`. For a typical
    stage0 ambient `V ~ 1e-12` and rupture onset pushing `VX ~ 1 m/s`,
    the guard rejects every step once nucleation starts.
  - Also applies the BP5 `SetDtMax(0.1 * seconds_per_year)` default
    (seas_driver.cpp:446), which is meaningless for a dynamic-rupture CFL
    around 2 ms.

**Trigger:**
The first TPV102 step in which any wave velocity exceeds
`v_guard_factor * (ambient noise)`, which happens immediately after
nucleation onset.

**Actual behavior (if implemented as the report suggests):**
Either (a) reject-cascade the step counter to its maximum, or (b) silently
corrupt the stepper's "per-node" indexing when applied to a composite
state, OR (c) the implementer is forced to add a dynamic-path branch into
`solver/time_stepper.hpp` (e.g., `SetPseudoStatePerNode(NUM_STATE + 3)` or
a V-guard-disabled mode) — editing a shared file. All three are violations
or silent bugs.

**Expected behavior:**
Explicitly state that `solver/time_stepper.hpp` cannot be reused as-is for
TPV102, and propose duplicating the RK4/RK45 stepping kernel into
`dynamic/dynamic_rk45.hpp` with a state-layout-free V-guard (or no V-guard
at all, since dynamic rupture has its own CFL-based stability criterion).

**Suggested fix:** in the report, Part B §B6 and Part C row 6:

```diff
 ### B6. `DormandPrinceRK45` (with adaptive dt and MPI error reduction)
-Already handles the MPI rank-consistent accept/reject
-(`MPI_Allreduce(err_norm, MPI_MAX)` — CLAUDE.md calls out this as critical),
-V-guard, PI controller, FSAL endpoint re-use. Drop-in replacement for
-TPV102's `for (step=0; step<nsteps; step++)` fixed-dt loop, yielding:
+BP5's `DormandPrinceRK45` (solver/time_stepper.hpp) hard-codes
+`state_per_node_` and a V-guard that scans "slip rate components" —
+both tied to the BP5 fault-state layout. It is NOT drop-in for a TPV102
+state vector that includes the bulk wave field. The dynamic path needs a
+new `dynamic/dynamic_rk45.hpp` that borrows the Butcher coefficients,
+PI controller, and MPI error reduction but omits the V-guard (dynamic
+rupture uses an explicit CFL criterion instead) and assumes a flat
+state. Do NOT add a "state layout" parameter to the shared BP5 stepper.
```

**Test case:**
```python
def test_R002_dp45_vguard_incompatible_with_bulk_state():
    """DormandPrinceRK45's V-guard interprets the state as
    [slip_c0..slip_c{spn-2}, psi] per node. Feeding a TPV102 composite
    state (bulk Q + fault slip + psi) triggers V-guard rejection
    on wave velocities. Document the incompatibility."""
    # Construct a trivial TDO with state = [Q_bulk ~ 1e-3 m/s velocities,
    #                                       slip_dip=0, slip_str=0, psi=0.4]
    # Set DP45 state_per_node=3, V-guard factor=100, ambient V_stage0 = 1e-12.
    # After one RK stage with any wave velocity > 100 * 1e-12 = 1e-10 m/s,
    # V-guard rejects. Assert that the stepper is not usable on this state.
    assert vguard_rejects_every_step_when_state_contains_bulk_velocities()
```

---

### [R-003] [MODERATE] [report §B3, §C row 2] — FaultScatter API assumes the BP5 owned-face / recv-face partition; using it without extension requires either modifying shared code or rebuilding the partition inside TPV102

**Category:** DEVIATION (partial violation / under-specified)

**Description:**
Report §B3 claims FaultScatter can be used in TPV102 "with FaultScatter in
the TPV102 path." Verified in `common/fault_scatter.hpp`:

- `SharedFaultCommBlock` carries `send_owned_faces` and `recv_local_faces`
  — indices into a per-rank owned-fault-face partition that BP5
  `ElasticityDomainOperator` builds (the "owned view" used by
  `RestrictToOwnedFault` / `ExpandOwnedToLocalFault`).
- `BeginScatter` packs a send buffer assuming the owned-data vector is
  indexed as `comps_per_dof * (owned_face * nbf_per_face + kk) + c`
  (fault_scatter.hpp:107-112). That layout is specific to the BP5
  ElasticityOperator's owned-fault numbering.
- TPV102 does NOT have an owned-fault-face partition in the BP5 sense.
  Both ranks run `FaultFaceFlux::Evaluate` on all shared-fault QPs;
  wave_operator.inl treats every rank as "owner" of the Elem1 side of
  each of its local shared faces. There is no BP5-style subset-and-scatter.

Consequence: "use FaultScatter directly" either (a) requires adding a
TPV102-specific owned-face partition builder inside `dynamic/` (fine, but
the report does not flag the ~100+ LOC of partition-construction logic),
or (b) requires extending FaultScatter to handle a different indexing
scheme (violates the rule).

Additionally, `fault_scatter.hpp:91` hard-codes `MPI_DOUBLE`. TPV102's
`MFEM_USE_SINGLE` builds use `float`. Not a BP5-modification issue, but
a separate portability concern that should be flagged alongside the B3
reuse claim.

**Trigger:**
Any attempt to construct a `FaultScatter` for dynamic-path shared-fault
communication using TPV102's existing shared-fault index space (which
doesn't match the `send_owned_faces` / `recv_local_faces` contract).

**Actual behavior (of the report as written):**
Glosses over the partition-rebuild cost. A naive implementer following
the report would either crash on an out-of-bounds index in
`send_bufs_[bi](j * block_size + kk * comps_per_dof + c)` or be forced
to edit `fault_scatter.hpp`.

**Expected behavior:**
State explicitly that TPV102 must build a new per-rank partition of its
fault QPs into owner-resolved "send" and "recv" sets before any
FaultScatter use, and place this partition-builder in
`dynamic/dynamic_fault_partition.hpp`. Either that, or duplicate a
trimmed-down scatter helper into `dynamic/dynamic_fault_scatter.hpp`.

**Suggested fix:** in the report, Part B §B3:

```diff
-### B3. `SharedFaultCommBlock` + `FaultScatter` (drop R-101 and `shared_fault_elem1_on_plus_`)
-
-`common/fault_scatter.hpp` resolves shared-fault cross-rank ownership via
-the vertex-triple key that is *already shared* (the shared_fault_key.hpp
-header explicitly plans this merge). With FaultScatter in the TPV102 path:
+### B3. Duplicate `SharedFaultCommBlock` + `FaultScatter` into `dynamic/`
+
+`common/fault_scatter.hpp` assumes BP5's owned-fault-face partition
+(send_owned_faces / recv_local_faces index into the BP5 ElasticityOperator
+owned view, with a hard-coded `MPI_DOUBLE` payload and
+`comps_per_dof * (face*nbf_per_face + kk) + c` packing). TPV102 has no
+such partition today (both ranks Evaluate on every shared QP).
+
+Do NOT reuse the shared header. Duplicate into
+`dynamic/dynamic_fault_scatter.hpp` with a TPV102-appropriate signature
+(per-shared-QP owner list, 9-component payload, MPITypeMap<real_t>), plus
+a new `dynamic/dynamic_fault_partition.hpp` that builds the owner-resolved
+partition using `shared_fault_key.hpp`'s FaceVertexKey (which is already
+intentionally bit-identical between dynamic/ and BP5's nested copy).
+
+With the duplicate in TPV102's path:
```

**Test case:**
```python
def test_R003_faultscatter_partition_contract():
    """FaultScatter's BeginScatter indexes send_bufs_ assuming
    comps_per_dof*(owned_face*nbf_per_face+kk)+c layout. TPV102's
    fault-QP layout is [interior_qps..., shared_qps...] with no
    owned-face partition. Document the contract mismatch."""
    # Construct minimal BP5-style blocks and a BP5-layout owned_data -> works
    # Construct TPV102-style fault-QP data with no owned partition -> OOB
    assert bp5_layout_scatter_ok()
    assert tpv102_style_scatter_without_new_partition_oob()
```

---

### [R-004] [MODERATE] [report §B1, §C row 1] — Canonical-rotation cache location is unspecified; a cache field inside `FaultBasisQPData` would modify BP5 source

**Category:** ASSUMPTION (ambiguous, risks BP5 modification)

**Description:**
Report §B1 proposes precomputing `T_can`, `Tinv_can` per QP and "storing
them inside WaveOperator" — but the obvious implementation choice is to
extend `FaultBasisQPData` (in `fault/fault_basis.hpp`, shared with BP5)
with two `DenseMatrix` members. That would modify a BP5 header and break
ABI/layout assumptions for existing BP5 tests (`seas_test_fault_basis`,
`seas_test_fault_basis_qp_orthonormality`, `seas_test_fault_basis_dip_strike_symmetry`).

The report does not:
  - Forbid the natural `FaultBasisQPData` extension.
  - Specify the storage schema (e.g., "`std::vector<DenseMatrix>
    T_can_per_qp_` as a private member of `WaveOperator`, indexed by
    `[basis_idx][q]`").
  - Address that `DenseMatrix` is 9×9 (NUM_STATE×NUM_STATE) and the cache
    size is `(N_interior_fault_qps + N_shared_fault_qps) * 2 * 9 * 9`
    doubles, which for a production-scale mesh is non-trivial memory
    (several hundred MB at 400 ranks × 1e5 fault QPs × 1296 bytes). That
    memory consideration may be fine but should be stated.

**Trigger:**
A reader of the report implements B1 by adding the most ergonomic
`DenseMatrix T_can; DenseMatrix Tinv_can;` fields to `FaultBasisQPData`.
BP5 tests fail at link or runtime when `FaultBasis::Compute` returns
default-constructed matrices that downstream BP5 code didn't expect.

**Expected behavior:**
State explicitly: cache lives in `WaveOperator` as a private member, NOT
in `FaultBasisQPData`. Give the exact field name and type.

**Suggested fix:** Part B §B1:

```diff
 ### B1. `FaultBasis` as the *sole* source of per-QP rotations (not a supplement)
 ...
-Change: precompute the canonical rotation matrices per QP at setup time
-(one `DenseMatrix T_can[nq_total]`, one `DenseMatrix Tinv_can[nq_total]`,
-stored inside `WaveOperator`). The 4 face-flux routines then read from the
-cache instead of rebuilding.
+Change: in `WaveOperator` (file `dynamic/wave_operator.hpp`, PRIVATE
+members), add `std::vector<DenseMatrix> cached_T_can_` and
+`std::vector<DenseMatrix> cached_Tinv_can_`, both sized
+`num_fault_faces * nbf_per_face` with indexing
+`[basis_idx * nbf_per_face + q]`. Populate at the end of the constructor
+right after `ComputeQPBasisShared`. Do NOT extend `FaultBasisQPData`
+(shared with BP5) with cache fields — that would alter the struct layout
+referenced by BP5's `seas_test_fault_basis*` unit tests and change the
+initialization contract.
```

**Test case:**
```python
def test_R004_cache_location_preserves_faultbasisqpdata_layout():
    """FaultBasisQPData must remain unchanged so BP5 tests still pass."""
    # Assert sizeof(FaultBasisQPData) matches pre-change snapshot
    # Run seas_test_fault_basis* under the new cache implementation
    assert sizeof_FaultBasisQPData_unchanged()
    assert bp5_fault_basis_tests_pass()
```

---

### [R-005] [MODERATE] [report §B4, §C row 3] — Duplicating `rate_state_fault.hpp` (~1,059 LOC) is correct under the rule but the report presents it as low-friction

**Category:** QUALITY / COST-TRANSPARENCY

**Description:**
Report §B4 proposes `dynamic/dynamic_rate_state_fault.hpp` as "duplicate
of rate_state_fault.hpp with ComputeRHS signature taking a Riemann trial
traction." `rate_state_fault.hpp` is 1,059 LOC. The duplicate will inherit:

- 2 constructors (BP2 + BP5) with diverging BP5-specific branches;
- `PreInit`, `Init`, `ComputeRHS` with BP2/BP5 compile-time dispatch;
- Monitoring/diagnostics (`monitor_interval_`, `monitor_stations_`);
- The `FRIC-GUARD` NaN detection block;
- The psi-space vs theta-space branch.

For TPV102, only a small subset is relevant (BP5-style vector path, psi
only, no BP2 branch, no per-node monitor, friction is already per-QP via
`FaultFaceFlux::ComputeTrialTraction`). The report doesn't indicate which
parts are kept, which are dropped, and where the duplicate will diverge
from the BP5 original over time.

Report §C row 3 rates this as "Medium-high" risk but does not distinguish
"correctness risk" (the actual algorithm change) from "maintenance risk"
(~1000 LOC of BP5-style code that must be kept in sync or allowed to
drift).

**Trigger:**
An implementer following the report literal-minded, duplicating all of
`rate_state_fault.hpp`. Months later, a BP5 fix lands in the original
(e.g., a numerical guard on sigma_n_eff) and silently does not propagate
to the dynamic duplicate.

**Expected behavior:**
The report should (a) narrow the duplicate to a minimal BP5-vector-psi
subset (target: ~300 LOC, not ~1000), (b) list which BP5 methods the
duplicate must NOT include, and (c) acknowledge the long-term
maintenance-drift cost explicitly as a downside of the rule.

**Suggested fix:** Part B §B4:

```diff
 ### B4. `RateStateFaultOperator<MeshType, 2>` as the dynamic friction operator
 ...
-Proposed refactor (follows the CLAUDE.md rule "duplicate into a new file
-under `dynamic/`, don't modify BP5 source"):
-
-```
-dynamic/dynamic_rate_state_fault.hpp    ← duplicate of rate_state_fault.hpp
-                                          with ComputeRHS signature taking a
-                                          Riemann trial traction instead of
-                                          tau_pre+elastic, and a pre-stress
-                                          discriminator for the total-Q path.
-```
+Proposed refactor (CLAUDE.md rule: no BP5 source edits, new file under
+`dynamic/`):
+
+  dynamic/dynamic_rate_state_fault.hpp    ~300 LOC, strictly BP5-vector-psi
+                                          subset of rate_state_fault.hpp.
+
+Explicitly OMIT from the duplicate:
+  - BP2 scalar constructor and `if constexpr (SlipComponents == 1)` branch;
+  - `monitor_interval_`, `monitor_stations_`, per-station traction print block;
+  - Theta-space (non-psi) branch (TPV102 is always psi);
+  - `PreInit`/`Init` — TPV102 state is constructed by `InitializeFaultDOFs`
+    / `InitializeStateTotal`, not by the Tandem 4-phase init;
+  - `SetSlip`/`GetSlip` — TPV102 carries slip directly in DOFData/state.
+
+Maintenance risk: any BP5-side fix to the friction solve, NaN guard, or
+sigma_n_eff branch must be manually mirrored into the dynamic duplicate.
+This is an accepted cost of the no-BP5-edit rule (per
+`feedback_tpv102_bp5_no_shared_edit.md`).
```

**Test case:**
```python
def test_R005_dynamic_rate_state_minimal_subset():
    """The dynamic duplicate should be <= ~400 LOC and must not contain
    BP2-scalar / monitor / theta-space code."""
    src = open("dynamic/dynamic_rate_state_fault.hpp").read()
    assert "SlipComponents == 1" not in src
    assert "monitor_stations_" not in src
    assert "use_psi_ = false" not in src
    assert lines(src) < 400
```

---

### [R-006] [MODERATE] [report §A5, Part D] — The `shared_fault_key.hpp`-merge claim presupposes a future BP5 edit the report does not acknowledge

**Category:** DEVIATION (the described end-state requires editing BP5)

**Description:**
Report §A5 quotes `shared_fault_key.hpp`'s header: "BP5's
`ElasticityOperator` currently owns an equivalent, bit-identical
`FaceVertexKey` nested struct and `MakeFaceKey` member function in
`miniapps/seas/domain/elasticity_operator.hpp` (circa lines 483-514). The
long-term plan is for BP5 to ALSO include this header and drop its nested
copies." The merge itself modifies `elasticity_operator.hpp` (a BP5 file).

Report §C row 2 ("move shared faults to `SharedFaultCommBlock` +
`FaultScatter`") implicitly depends on a canonical key that works across
both paths. If the dynamic path continues to use `dynamic/shared_fault_key.hpp`
while BP5 keeps its nested copy, that's fine. But the report's phrasing
"the vertex-triple key that is *already shared*" is misleading — it is
NOT shared today; there are two independent bit-identical copies that the
report implicitly plans to unify.

**Trigger:**
An implementer reads the report and performs the merge (removes BP5's
nested `FaceVertexKey` and includes `dynamic/shared_fault_key.hpp`),
thinking the report endorses that merge as part of B3. That edit
modifies BP5 source.

**Expected behavior:**
State explicitly that the "bit-identical copies are kept for now; the
header merge is NOT part of this plan and would be a separate,
BP5-modifying change to land only with explicit approval."

**Suggested fix:** Part D (and implicitly Part A5):

```diff
+### D0. Deferred: `shared_fault_key.hpp` BP5-merge
+`shared_fault_key.hpp` notes a future plan to have BP5 drop its nested
+`FaceVertexKey` copy and include the shared header. That merge modifies
+`domain/elasticity_operator.hpp` (BP5 source) and is explicitly NOT part
+of this proposal. For the duration of this migration the two
+bit-identical copies are kept; a refactor to unify them requires its
+own approval cycle.
```

**Test case:**
```python
def test_R006_shared_fault_key_merge_not_performed():
    """domain/elasticity_operator.hpp must still contain its nested
    FaceVertexKey copy after any of the B1..B6 landings."""
    assert "struct FaceVertexKey" in \
           open("domain/elasticity_operator.hpp").read()
```

---

### [R-007] [MODERATE] [report §C Row 2] — Deletion list includes BP5-shared constructs (indirect BP5-breakage risk)

**Category:** BUG (if followed literally)

**Description:**
Part C Row 2 says: "delete `shared_fault_elem1_on_plus_`, R-101 verifier,
`ComputeSharedFaceFluxRHS` core logic." All three live in
`dynamic/wave_operator.{hpp,inl}` — safe. But the row does NOT scope the
deletion carefully:

- `ComputeSharedFaceFluxRHS` also handles **non-fault shared faces**
  (wave_operator.inl:1718 "Non-fault shared face: standard welded Godunov
  flux."). Deleting the "core logic" wholesale removes the non-fault
  interior-face path through shared faces and breaks all non-fault DG
  communication across MPI seams — which is how the elastodynamics wave
  propagates across partition boundaries. The row as written wipes out
  wave propagation, not just the fault branch.
- `shared_fault_elem1_on_plus_` is only populated in the constructor
  conditional on `bc_.fault_attr > 0`. If deleted without guarding the
  read sites inside `ComputeSharedFaceFluxRHS` under the new code path,
  a BP5 test that happens to use the shared WaveOperator in a fixture
  (`seas_test_parallel_wave_operator`) would segfault on the first
  shared fault QP.

**Trigger:**
Implementer follows Row 2 literally: deletes the function body, deletes
the member. Any run with more than one rank and a shared non-fault face
loses cross-rank flux.

**Expected behavior:**
Scope the deletion to the fault-specific branches only, not the whole
function. State that the "non-fault shared face welded Godunov flux"
branch must be preserved verbatim.

**Suggested fix:** Part C Row 2:

```diff
-| 2 | B3 — move shared faults to `SharedFaultCommBlock` + `FaultScatter`;
-    delete `shared_fault_elem1_on_plus_`, R-101 verifier,
-    `ComputeSharedFaceFluxRHS` core logic |
+| 2 | B3 — duplicate SharedFaultCommBlock + FaultScatter into `dynamic/`
+    (see R-003). In `ComputeSharedFaceFluxRHS`, replace ONLY the
+    fault-branch body (wave_operator.inl:1540-1714). Preserve the
+    non-fault branch (wave_operator.inl:1716-1732, welded Godunov on
+    shared interior faces) verbatim — that's how wave energy crosses
+    partition seams. Delete `shared_fault_elem1_on_plus_` only AFTER
+    confirming it is read nowhere outside the fault branch.
+    The R-101 verifier can be removed only AFTER the new
+    single-owner-scatter pattern is validated end-to-end. |
```

**Test case:**
```python
def test_R007_nonfault_shared_flux_preserved():
    """After B3, MPI-seam wave propagation on non-fault shared faces
    must still produce the pre-change Godunov flux."""
    # Run seas_test_parallel_wave_operator on 4 ranks with a fixture that
    # has non-fault shared faces and non-zero initial Q far from the fault.
    # Assert rank-boundary wave amplitudes unchanged.
    assert parallel_wave_operator_test_passes()
```

---

### [R-008] [LOW] [report §A1] — "Intentional duplication" attribution is overstated

**Category:** QUALITY (correctness of the report's justification claim)

**Description:**
The report says (A1): 'The code itself documents this (wave_operator.inl:1769):
"The duplication (vs. refactoring ComputeFaceFluxRHS into a shared core)
is intentional..."' — but the cited comment is specifically about
**ADER vs RK4** duplication (the ADER variants are at 1769; the quoted
block is inside the ADER introduction). It says nothing about the
**interior vs shared** duplication, which is the root of pepper hypothesis
H1. The report conflates the two and attributes both to the "intentional"
comment.

**Trigger:**
A reviewer skeptical of the pepper-hypothesis rationale reads the
report, follows the line reference, sees the comment is about ADER-only
duplication, and discounts the rest of the argument.

**Expected behavior:**
Split the accounting: ADER vs RK4 is intentional and documented; interior
vs shared is undocumented and is the genuine H1/H2 driver.

**Suggested fix:** A1 opening:

```diff
 Line numbers: `ComputeFaceFluxRHS` (831–1322, 491 LOC), `ComputeSharedFaceFluxRHS`
 (1323–1773, 450 LOC), `ComputeADERFaceFluxRHS` (1774–2106, 332 LOC),
 `ComputeADERSharedFaceFluxRHS` (2107–2393, 286 LOC).

-The code itself documents this (wave_operator.inl:1769):
+The ADER/RK4 split is explicitly documented as intentional
+(wave_operator.inl:1769):
 > The duplication (vs. refactoring ComputeFaceFluxRHS into a shared core)
 > is intentional: `wave_operator.inl` is on the CLAUDE.md "Files Requiring
 > Extreme Care" list and the RK4 path must remain byte-identical.
+
+The interior/shared split has NO comparable justification in the code;
+it is an organic copy-paste that pre-dates R-801 and is the principal
+H1/H2 pepper vector identified by the companion report.
```

**Test case:** Not applicable (documentation correctness; non-executable).

---

### [R-009] [LOW] [report §A3] — The fluctuation-path `ApplyNucleation` alternative is not symmetry-equivalent to `ApplyNucleationTotal`

**Category:** POSSIBLE-BUG (the report recommends a simpler alternative that may not be numerically equivalent)

**Description:**
A3 says: "An alternative — keep nucleation in `DOFData.tau2_0` as the
fluctuation path does — exists in `tpv102_setup.hpp::ApplyNucleation` and
is simpler by a factor of ~5."

This is true *for the fluctuation Q path*. Under the v9.3.0 total-Q
migration (I-06), bulk Q carries the pre-stress tensor in global
coordinates; `EvaluateTotal` asserts `data.sigma_n0 == 0 && data.tau1_0
== 0 && data.tau2_0 == 0` (fault_face_flux.cpp:266-274). Writing
`data.tau2_0 = TPV102Params::tau_ini + dtau` under total-Q semantics
would trigger that `MFEM_ASSERT` in debug builds and silently double-count
pre-stress in release builds.

The report presents the switch as "just simpler" without noting that it
requires also rolling back the total-Q migration for nucleation — which
means either (a) reverting nucleation to fluctuation semantics while
other BCs stay total-Q (semantic split the I-06 migration explicitly
rejected), or (b) a bigger rework than A3 implies.

**Trigger:**
Implementer reads A3, replaces `ApplyNucleationTotal` with
`ApplyNucleation`, runs in debug mode, hits the EvaluateTotal assert.

**Expected behavior:**
A3 should note that simplifying nucleation to the DOFData path is
incompatible with the current I-06 total-Q invariant, and state which
additional change (revert nucleation-only to fluctuation, or extend
FaultFaceFlux with a nucleation-delta input parameter, etc.) goes
alongside the simplification.

**Suggested fix:** A3 closing paragraph:

```diff
-... The migration can be *scoped to the BC path* without dragging
-nucleation out of `DOFData`.
+... However, reverting nucleation to `DOFData.tau2_0` while leaving BCs
+on the total-Q path violates the `EvaluateTotal` invariant
+`data.tau2_0 == 0` (fault_face_flux.cpp:266-274) — the assert fires in
+debug and pre-stress is double-counted in release. The simplification
+therefore requires either adding a nucleation-delta parameter to
+EvaluateTotal (extending FaultFaceFlux signature) or restricting
+nucleation-delta into a small dedicated field on DOFData that
+EvaluateTotal explicitly consumes (e.g., `data.tau2_nuc`). Document
+this dependency when planning the simplification.
```

**Test case:**
```python
def test_R009_nucleation_tau2_0_violates_evaluate_total():
    """Writing tau2_0 = tau_ini + dtau under total-Q breaks the
    EvaluateTotal invariant."""
    dof = DOFData()
    dof.tau2_0 = 1e6
    # EvaluateTotal asserts tau2_0 == 0 in debug
    with pytest.raises(AssertionError, match="tau2_0"):
        fault_flux.EvaluateTotal(dof, Q_plus, Q_minus, ...)
```

---

## Summary

- **Critical issues:** 2 (R-001 FaceQuadrature rule mismatch, R-002 DP45
  not state-layout-agnostic — both would cause BP5 edits or silent
  corruption if the report is implemented literally)
- **Moderate issues:** 5 (R-003 FaultScatter API mismatch; R-004
  FaultBasisQPData cache risk; R-005 duplicate-LOC transparency;
  R-006 shared_fault_key merge ambiguity; R-007 deletion-scope bug in
  Part C Row 2)
- **Low issues:** 2 (R-008 over-attributed "intentional" quote; R-009
  A3 alternative breaks EvaluateTotal invariant)
- **Plan compliance** (to the `/code-review` charter "reuse CANNOT modify
  BP5 or break BP5 run, else duplicate"): **PARTIAL**. Of the 6 proposed
  reuses, three (B2 FaceQuadrature, B3 FaultScatter, B6 DP45) would
  require BP5-shared edits if taken literally. Two (B1 cache placement,
  B4 RSOP duplication) are under-specified on how to avoid BP5 edits.
  Only B5 (new top-level `SEASDynamicOperator`) is unambiguously safe
  as written.
- **Verdict:** **FAIL — must revise before executing.** The report as
  written would land in one of three failure modes: (a) BP5-source
  edits that violate `feedback_tpv102_bp5_no_shared_edit.md`, (b) BP5
  runtime breakage in shared unit tests (R-001, R-004, R-007), or (c)
  silent numerical corruption in the TPV102 runs the report is trying
  to stabilize (R-001 aliasing, R-002 V-guard cascades, R-009 pre-stress
  double-count).

## Unreviewed Areas

- I did not verify that the 6-step `SEASDynamicOperator::Mult` pseudocode
  in §B5 is physically well-posed for dynamic rupture (it borrows the
  BP5 quasi-dynamic structure but replaces the elliptic solve with a
  hyperbolic stage). That is a plan-level correctness question beyond
  the `/code-review` charter, which is scoped to "does the reuse break
  BP5?".
- I did not exhaustively verify that every private member of
  `ElasticityDomainOperator` consumed by `FaultGeometry` /
  `RateStateFaultOperator` has a dynamic-path equivalent. The report
  does not propose reusing `FaultGeometry` for the dynamic path (it is
  mentioned once in §B passim but not in the §C rollout table), so the
  BP5-impact question does not arise.
- I did not audit `dynamic/tpv102_setup_total.hpp` (697 LOC) for its own
  internal correctness; that is the scope of the original pepper-report
  hypotheses (H3, H4), not of this review.

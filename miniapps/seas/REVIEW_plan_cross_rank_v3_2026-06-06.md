# Code Review (v3 plan): General distributed-material (DMF) design

## Review Scope
- Plan: `document/bimaterial_dev/PLAN_cross_rank_material_exchange_2026-06-06.md` (v3, the DMF
  unified design).
- Prior rounds: `REVIEW_plan_cross_rank_material_exchange_2026-06-06.md` (R-001..R-007),
  `REVIEW_plan_cross_rank_round2_2026-06-06.md` (R-201..R-209).
- Code cross-checked: `dynamic/heterogeneous_material.{hpp,cpp}` (MaterialField modes,
  EvalAt/At, halfspace coeff), `dynamic/bimaterial_wave_operator.inl` (BuildGodunovFluxPool_,
  AssignFault eps-offset at :1028-1048, the GridFunction skip at :1034), `dynamic/wave_operator.inl:36`
  (state basis = `BasisType::GaussLobatto`).
- Domain: CLAUDE.md, the B0/B3 memory (halfspace `sign((x-x0).n)>=0`).

## Findings

### [R-301] CRITICAL — Phase 1: `ProjectCoefficient` onto an order-p GLL space pollutes the fault-face nodes ⇒ wrong per-side material for a jump-at-fault material (TPV6)

**Category:** BUG

**Description:**
The plan materializes a Coefficient into `mat_gf_` via `ProjectCoefficient` on an L2 space with
`order = state order p` and `BasisType::GaussLobatto`. GLL nodes INCLUDE element endpoints/face
nodes. For a halfspace material the coeff is evaluated AT the fault-face nodes (y=0), where it is
AMBIGUOUS — `sign((x-x0).n) >= 0` returns "near" (strong) at y=0 regardless of side. So a far-side
(weak) fault element gets STRONG material on its fault-face nodes; the element's nodal triple is
MIXED (interior=weak, fault-face=strong). `GetVectorValue` at any point near the fault then
INTERPOLATES that mix → the per-side material is polluted toward the wrong side. This is the exact
ambiguity the current code avoids with the eps-offset (bimaterial_wave_operator.inl:1028-1048 —
"perturb a face-QP IP a small eps along the fault normal" — i.e. evaluate strictly INSIDE the
element). The plan's projection THROWS AWAY the eps-offset disambiguation.

It also breaks the plan's own "byte-exact for halfspace / L2 projection exact" claim: a jump is
not representable in the order-p space, so the projection is a Gibbs-like mix at every element the
interface passes through (the fault elements AND the non-conforming welded-interface elements).

**Trigger:** any halfspace / jump-at-element-face material (TPV6/7, any bi-material fault) — i.e.
the primary target of this work.

**Actual behavior (planned):** fault-element nodal material is a strong/weak mix; per-side
impedances polluted.

**Expected behavior:** each fault element's per-side material is the pure material on its side
(what the eps-offset interior eval gives today).

**Suggested fix:** make the DMF an L2 **order-0** (piecewise-constant, one triple per element)
field, POPULATED by the existing per-element interior-point logic (the same centroid eval that
fills `per_elem_lmr_`, and the eps-offset side-pick for the fault) — NOT `ProjectCoefficient`.
Concretely: compute the per-element triple at a point STRICTLY INSIDE the element (centroid for
bulk; the eps-offset-selected element for the fault), and set the order-0 dof to it. Order-0 is
discontinuous (preserves the jump) AND unambiguous (no face node is ever sampled).
```diff
- mat_fec_ = make_unique<L2_FECollection>(order_, dim_, BasisType::GaussLobatto);  // order p
- mat_gf_->ProjectCoefficient(vc);   // pollutes fault-face nodes for a jump material
+ mat_fec_ = make_unique<L2_FECollection>(0, dim_);                                 // order 0
+ // fill per element from an INTERIOR point (reuse the per_elem_lmr_ centroid logic):
+ for (int e=0; e<ne_; ++e) { mat_gf_->SetElementDofValues(e, per_elem_lmr_[e]); }  // one triple/elem
```
If sub-element resolution at the fault is later required, use an interior-BIASED custom projection
(offset boundary nodes inward before coeff eval) — NOT plain `ProjectCoefficient`. Note in the plan.

**Test case:**
```cpp
// np=1 (serial) is enough to expose the pollution:
void test_R301_halfspace_fault_element_not_polluted() {
   // one fault element wholly on the FAR (weak) side, fault face at y=0.
   // DMF per-side material for that element == the WEAK triple (vp=3750...), NOT a
   // strong/weak mix.  Fails for order-p ProjectCoefficient (face nodes = strong);
   // passes for order-0 interior-point fill.
}
```

---

### [R-302] CRITICAL — Phase 4: fault per-side must use the SAME interior representation on BOTH sides; the plan mixes peer-from-DMF with local-eps-offset

**Category:** BUG / DEVIATION (the v1 R-001 trap, re-surfacing)

**Description:**
Phase 4 says the LOCAL side uses `MaterialAtLocal_(elem, eps_offset_ip)` and the PEER side uses
`MaterialAtNbr_(nbr_idx, eps_offset_ip)`. If `MaterialAtLocal_` reads an order-p `mat_gf_` at the
eps-offset (sub-element) while `MaterialAtNbr_` reads the peer's value, the two sides are NOT a
consistent representation (one interpolates the local order-p field at a sub-element point, the
other reads the peer's). For a SYMMETRIC depth profile this re-introduces asymmetry
(`Zp_plus != Zp_minus`) — the σ_n-leak regression round-1 R-001 was about. With the order-0 fix
(R-301) both sides read a per-element constant ⇒ consistent ⇒ symmetric. The plan must STATE that
both sides use the identical per-element representation, and DROP the "eps-offset IP" read (order-0
has one value per element; the eps-offset only SELECTS which element is `+`/`-`, it does not change
the value).

**Trigger:** symmetric depth-profile fault on a partition seam (TPV31 at np>1, no fault-locality).

**Suggested fix:** in Phase 4, after using the eps-offset only to pick the `+`/`-` ELEMENT, read
each side's per-element constant: `MaterialAtLocal_(Elem1)` and `MaterialAtNbr_(nbr_idx)` (no ip
needed for order-0). Assert symmetry in the test.
```diff
- LOCAL: MaterialAtLocal_(Elem1, eps_offset_ip1);  PEER: MaterialAtNbr_(nbr_idx, eps_offset_ip2);
+ // eps-offset selects the +/- element only; the value is that element's order-0 constant:
+ plusMat  = (e1_plus ? MaterialAtLocal_(Elem1) : MaterialAtNbr_(nbr_idx));
+ minusMat = (e1_plus ? MaterialAtNbr_(nbr_idx) : MaterialAtLocal_(Elem1));
```

**Test case:** `test_R302_shared_fault_symmetric_depthprofile` (the Phase-4 symmetric gate) — a
depth-profile fault on a seam ⇒ `|Zp_plus-Zp_minus| <= 1e-12*Zp_plus`. Fails if the two sides use
different (sub-element vs per-element) reads.

---

### [R-303] MODERATE — Phase 2: `per_elem_lmr_` ↔ DMF data-flow is circular

**Category:** DEVIATION (spec inconsistency)

**Description:**
Phase 2 says `per_elem_lmr_[e] = MaterialAtLocal_(e, centroid)` (reads the DMF), but Phase 1
populates the DMF from the per-element material. With the R-301 order-0 fix the DMF IS the
per-element material, so making `per_elem_lmr_` read the DMF is circular. The correct data flow:
compute the per-element triple ONCE (the existing centroid logic, extended to GridFunction), store
it in `per_elem_lmr_`, COPY it into the order-0 `mat_gf_`, `ExchangeFaceNbrData`; then only the
NEIGHBOUR and FAULT-peer reads use `MaterialAtNbr_`. `per_elem_lmr_` stays the source, not a
consumer.

**Suggested fix:** Phase 1 fills `per_elem_lmr_` (centroid) for ALL modes first, then sets
`mat_gf_` order-0 dofs from `per_elem_lmr_`. Phase 2 changes only `shared_face_neighbour_material_`
(→ `MaterialAtNbr_`) and removes the abort; `per_elem_lmr_`/`BuildGodunovFluxPool_` keep their
local centroid source (no DMF read).

**Test case:** `test_R303_per_elem_lmr_equals_dmf_local` — `mat_gf_` local value == `per_elem_lmr_`
to the bit (they must be the same source).

---

### [R-304] MODERATE — Plan: byte-exact claims are overstated; restate for the order-0 design

**Category:** DEVIATION (correctness of the stated contract)

**Description:**
The plan claims "byte-exact for Constant and piecewise-constant Coefficient ... L2 projection is
EXACT." With order-p `ProjectCoefficient` this is false (R-301). With the order-0 interior-point
fill (R-301 fix) the correct contract is:
- Constant: byte-exact (one constant everywhere).
- Halfspace (TPV6): byte-exact for the FLUX POOL + neighbour + per-side fault (uniform per
  element; centroid == eps-offset value for a piecewise-constant-per-element field).
- Depth profile (TPV31/TPV102): the FAULT per-side eta changes from "eps-offset at the fault-QP
  depth" (current) to "element-centroid depth" (order-0) ⇒ a SMALL re-baseline (NOT byte-exact),
  even though symmetry is preserved. The BULK already used the centroid ⇒ byte-exact there.
State this precisely so the implementer does not assert TPV31 fault byte-exactness and chase a
false failure.

**Suggested fix:** replace the Constraints "byte-exact / L2 projection exact" bullet with the
three-case statement above; mark the TPV31 fault-eta change as a documented re-baseline.

**Test case:** `test_R304_tpv6_fault_byteexact_vs_legacy` (halfspace per-side == legacy eps-offset
to 1e-12) + a documented (non-gating) TPV31 fault-eta diff note.

---

### [R-305] MODERATE — `--material-legacy-eval` fallback lands too late (Phase 5); byte-exact is unrecoverable during Phases 2-4

**Category:** EDGE_CASE (phasing)

**Description:**
The plan introduces the legacy-eval fallback only in Phase 5, but Phase 2 already flips the BULK to
the DMF (default) and Phase 4 the FAULT — so between Phase 2 and Phase 5 there is NO way to run the
byte-exact oracle (e.g. to re-baseline TPV31 or bisect a regression). The fallback is also the
oracle the Phase-5 `dmf_vs_legacy` gate needs.

**Suggested fix:** introduce `use_dmf_` + `--material-legacy-eval` in Phase 1 (default DMF=on for
the bulk only after Phase 2 wiring), so every intermediate phase can be A/B'd against legacy.
Move the `dmf_vs_legacy_constant_halfspace` gate to run as soon as Phase 2 lands.

**Test case:** `test_R305_dmf_vs_legacy_oracle` runnable after Phase 2 (bulk) and Phase 4 (fault).

---

### [R-306] MODERATE — GridFunction source → order-0 DMF population is unspecified (and `per_elem_lmr_` currently can't be built for GridFunction)

**Category:** ASSUMPTION

**Description:**
Two gaps for the CVM (GridFunction) case the plan calls primary:
1. `BuildGodunovFluxPool_` fills `per_elem_lmr_` via `EvalAt(center)`, which ABORTS for
   GridFunction (heterogeneous_material.hpp:206). So `per_elem_lmr_` cannot be built for a
   GridFunction material today — Phase 1/2 must add a GridFunction branch (use `At(elem, dof)` /
   `GetValue(elem, center)` to get the per-element value).
2. The plan's "direct dof copy" from the 3 scalar order-p source GFs into the order-0 vdim-3
   `mat_gf_` is a dimension/layout mismatch (order-p → order-0, scalar → vdim-3). It must REDUCE
   each source GF to one per-element value (centroid `GetValue` or element average), then write
   the order-0 dof.

**Suggested fix:** specify a single per-element reduction used for ALL modes — `MaterialPerElem_(e)`
→ {Constant: the constant; Coefficient: `coeff->Eval(T, center)`; GridFunction:
`gf->GetValue(e, center)`} — fill both `per_elem_lmr_` and the order-0 `mat_gf_` from it.

**Test case:** `test_R306_gridfunction_per_elem` — a GridFunction material ⇒ `per_elem_lmr_[e]` ==
the source GF value at element e's centroid (no abort), all components.

---

### [R-307] LOW — Phase 4 "eps-offset IP" wording is misleading under order-0

**Category:** QUALITY

**Description:**
Under the order-0 fix the eps-offset does NOT change the material value (one per element); it only
selects which element is `+`/`-`. The plan's repeated "MaterialAt*(.., eps_offset_ip)" implies the
ip matters. Clarify: the existing `interior_fault_elem1_on_plus_` / `shared_fault_elem1_on_plus_`
already encode the side; no per-QP material variation remains.

**Suggested fix:** reword Phase 4 to "use the precomputed `*_elem1_on_plus_` to pick `+`/`-`; read
each element's order-0 constant; ip is unused for the value."

---

### [R-308] LOW — DMF is built for homogeneous runs too; order-0 makes it cheap but the collective remains

**Category:** ASSUMPTION

**Description:**
Order-0 (1 dof/elem × 3) makes the storage negligible (the earlier order-p concern is resolved),
but `mat_gf_->ExchangeFaceNbrData()` adds one collective + a face-nbr buffer to EVERY run incl.
homogeneous BP5/TPV205, where material was a single scalar. Byte-exact for Constant holds; the only
cost is the (small) exchange. Acceptable, but the plan should note it (and that Constant could skip
the exchange as a future micro-opt — the legacy path already does).

**Suggested fix:** note the homogeneous cost; optionally skip the exchange when `n_shared==0` ... but
keep the collective on ALL ranks if ANY rank has shared faces (R-209) — so guard by a global
reduction, not a local `n_shared`.

---

## Summary
- Critical issues: 2 (R-301 ProjectCoefficient pollutes fault-face nodes → wrong per-side material;
  R-302 fault sides must use the identical per-element representation)
- Moderate issues: 4 (R-303 circular data flow; R-304 byte-exact claims; R-305 fallback phasing;
  R-306 GridFunction population)
- Low issues: 2 (R-307 eps-offset wording; R-308 homogeneous overhead)
- Plan compliance: N/A (plan). The DMF concept is sound and meets the "one general path"
  directive — BUT the chosen materialization (order-p `ProjectCoefficient`) is WRONG for the
  jump-at-fault materials that are the whole point. The fix (order-0, interior-point fill, same
  representation both sides) keeps the unified design and is correct.
- Verdict: **PASS WITH FIXES (plan must be revised before implementation)** — switch the DMF to
  L2 order-0 populated by the per-element interior-point logic (R-301/R-302/R-303/R-306), restate
  the byte-exact contract (R-304), and land the legacy oracle early (R-305).

## Note on what v3 got RIGHT (so the fix preserves it)
- One representation consumed uniformly (local + `ExchangeFaceNbrData`), no per-mode branching in
  the flux/fault hot path — keep this.
- L2 (discontinuous) so jumps are preserved — keep, but at order 0 (not p).
- GridFunction handled by the same path — keep (via the per-element reduction, R-306).
- The legacy direct-eval as a byte-exact oracle — keep, land it earlier (R-305).

## Unreviewed Areas
- The exact MFEM `GetVectorValue` / order-0 `GetFaceNbrElementVDofs` round-trip — gated by the
  Phase-1 np=2 test at implementation.
- Whether any production config actually supplies a GridFunction material today (grep found none in
  drivers/spatial) — CVM is future; the GridFunction path is built+gated but not yet exercised by a
  real config.

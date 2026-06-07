# Code Review (round 4): v4 plan — fresh adversarial pass for critical bugs

## Review Scope
- Plan: `document/bimaterial_dev/PLAN_cross_rank_material_exchange_2026-06-06.md` (v4, the uniform
  material-accessor design).
- Prior rounds: v3 review (`REVIEW_plan_cross_rank_v3_2026-06-06.md`, R-301..R-308) — v4 was the
  fix.  Round-2 (`_round2`, R-201..R-209) test-gating still applies.
- Code cross-checked: `bimaterial_wave_operator.inl` (offset_ip lambda :1047-1062, AssignFault
  interior loop, BuildGodunovFluxPool_), `heterogeneous_material.hpp` (modes, public coeff/GF
  members), `wave_operator.inl:36` (GLL basis).
- Verified v4 RESOLVED the v3 criticals: R-301 (no projection ⇒ no fault-face pollution), R-302
  (eps-offset both sides ⇒ exact symmetry on the asymmetric mesh), R-303 (no circular DMF), R-304
  (byte-exact restated), R-305 (no separate fallback needed — local eval unchanged), R-306
  (per-element read via the accessor, GridFunction no longer aborts).  Confirmed each is addressed
  in the v4 text.  Now hunting for NEW bugs introduced by v4.

## Findings

### [R-401] MODERATE [POSSIBLE] — Phase 4: `offset_ip` on the PEER (face-neighbour) transform `ftr->Elem2` needs a valid `Jacobian()` at an interior IP

**Category:** ASSUMPTION

**Description:**
`offset_ip` (inl:1047) calls `T.SetIntPoint(&ip); DenseMatrixInverse Jinv(T.Jacobian());` and uses
`T.GetDimension()`.  For an INTERIOR fault face, `T = ftr->Elem2` is a LOCAL element — Jacobian
valid.  v4 Phase 4 calls `offset_ip(*ftr->Elem2, ftr->GetElement2IntPoint(), in2)` where, for a
SHARED fault face, `ftr->Elem2` is the FACE-NEIGHBOUR element transform.  DG face integrators use
`Elem2` only at FACE QPs; v4 evaluates its `Jacobian()` at an eps-offset INTERIOR point.  This is a
new dependency — if the face-nbr `Elem2` transform is not fully populated for interior IPs (or its
Jacobian is degenerate), the peer eps-offset is garbage.  Likely OK (MFEM sets `Elem2` as a full
`IsoparametricTransformation` for the face-nbr element after `ExchangeFaceNbrData`), but it MUST be
verified at implementation, not assumed.

**Trigger:** any shared fault face at np>1 (the Phase-4 path).

**Suggested fix:** after `pmesh.ExchangeFaceNbrData()`, assert the peer transform is usable before
the loop, and rely on the gate below.  If `ftr->Elem2->Jacobian()` is invalid, compute the peer
eps-offset PHYSICAL point from the LOCAL side instead (`x_qp - eps*n_phys`) and evaluate the
Coefficient via a point-locating transform — but only if the assert fails.
```cpp
// in the shared-fault loop, before offset_ip on the peer:
MFEM_VERIFY(ftr->Elem2 != nullptr && ftr->Elem2->GetDimension() == dim_,
            "AssignFault(shared): peer transform not usable for eps-offset at sf=" << sf);
```

**Test case:** `test_R401_shared_fault_symmetric_depthprofile` (the Phase-4 symmetric gate, run
np=2 on an ASYMMETRIC partition) — if `offset_ip` on the peer is wrong, `Zp_plus != Zp_minus`
(or a crash); if correct, `|Zp_plus-Zp_minus| <= 1e-12*Zp_plus`.  This single gate both proves the
fix and exercises R-401.

---

### [R-402] MODERATE — Phase 1: `MaterialAtLocal_`/`MaterialAtNbr_` GridFunction path assumes SCALAR material GFs; not asserted

**Category:** ASSUMPTION

**Description:**
v4's GridFunction branch uses `lambda_gf->GetValue(elem, ip)` (a scalar read).  `MaterialField`
docs say the 3 GFs "share the SAME ParFiniteElementSpace used by the wave operator's state."  The
state space is vdim=NUM_STATE (9).  If a material GF is actually built on the vdim-9 state space,
`GetValue(elem, ip)` returns component 0 of a 9-vector, not the scalar material — silently wrong.
No GridFunction material is constructed anywhere today (grep found none in drivers/spatial), so
this is forward-looking, but the plan must PIN the contract: the material GFs are SCALAR (vdim=1)
on an L2 space of the state order.

**Trigger:** a future CVM GridFunction material (the stated goal).

**Suggested fix:** in `MakeGridFunction` / the accessor, assert `gf->FESpace()->GetVDim() == 1`
and the FE order == state order; document "material GFs are scalar L2(order p)".
```cpp
MFEM_VERIFY(lambda_gf->FESpace()->GetVDim() == 1,
            "MaterialField GridFunction must be scalar (vdim=1)");
```

**Test case:** `test_R402_gridfunction_scalar_space` — construct with a vdim=1 GF (passes) and a
vdim>1 GF (aborts).

---

### [R-403] LOW — Phase 1: `MaterialAtLocal_` is declared `const` but `mesh_.GetElementTransformation` / `Eval` need non-const

**Category:** QUALITY

**Description:**
The interior path uses `const_cast<MeshType &>(mesh_)` (inl:1067) precisely because the mesh
accessor is non-const.  A `const` `MaterialAtLocal_` calling `mesh_.GetElementTransformation(e)`
will not compile without the same const_cast.  Mirror the existing pattern.

**Suggested fix:** `auto &m = const_cast<MeshType&>(mesh_); ElementTransformation *T =
m.GetElementTransformation(elem);` inside the accessor (matches inl:1067).

---

### [R-404] LOW — accessors should assert `material_ != nullptr`

**Category:** EDGE_CASE

**Description:**
`AssignFault` currently guards `if (material_ == nullptr) return;` (inl:1029).  The new accessors,
called from the flux pool + neighbour + fault, should assert `material_` is set (the bimaterial
operator always has one, but a defensive assert prevents a null-deref if a future path omits it).

**Suggested fix:** `MFEM_ASSERT(material_ != nullptr, "MaterialAt*: material_ not set");` at the top
of each accessor.

---

### [R-405] LOW — Phase 1: 3 separate `ExchangeFaceNbrData` calls (one per material GF)

**Category:** QUALITY

**Description:**
v4 exchanges `lambda_gf`, `mu_gf`, `rho_gf` independently ⇒ 3 collectives at the ctor (GridFunction
mode only).  Correct but wasteful; a combined vdim-3 GF would be 1 collective.  Not a bug (one-time
ctor cost), note as a future micro-opt.  Keep all 3 on ALL ranks (R-209 collective symmetry).

---

## Summary
- Critical issues: 0 (v4 resolved the v3 criticals R-301/R-302; no new CRITICAL found).
- Moderate issues: 2 (R-401 peer-transform Jacobian — verify+gate; R-402 scalar-GF contract).
- Low issues: 3 (R-403 const_cast; R-404 null assert; R-405 combined exchange).
- Plan compliance: the v4 design is internally consistent and meets the user goal (one uniform
  workflow, GridFunction first-class, byte-exact for existing TPV, exact fault symmetry on any
  mesh).  The remaining items are implementation-time verifications, not design flaws.
- Verdict: **PASS WITH MINOR FIXES** — v4 is implementable.  Address R-401 (assert + the symmetric
  gate proves it) and R-402 (scalar-GF assert) during implementation; R-403/404/405 are mechanical.

## The one open decision for the user (carried from v4 R-7, not a bug)
v4 deliberately does NOT "materialize every material into one GridFunction" (the v3/earlier
directive) — it reads each representation natively through one accessor — because a jump-at-fault
Coefficient cannot be projected onto a GridFunction without either fault-face-node pollution
(R-301) or on-asymmetric-mesh asymmetry that seeds the TPV31 σ_n leak (R-302, and TPV31 mixed-flux
runs on the asymmetric 50 m mesh).  The user's GOAL (one general workflow, GridFunction/CVM
first-class, no per-problem branching) is met; only the leaf accessor knows the representation.
This is the one point to confirm before implementation.

## Unreviewed Areas
- The actual MFEM behaviour of `ftr->Elem2->Jacobian()` at an interior IP (R-401) — gated, not yet
  run.
- The vdim/order of a real CVM GridFunction (none exists today; R-402 pins the contract for when
  one does).

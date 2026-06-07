# Code Review: PLAN_bimaterial_central_contrast_guard_2026-06-06 (plan review)

## Review Scope
- Plan: `document/mixed_flux_dev/PLAN_bimaterial_central_contrast_guard_2026-06-06.md`
- Files reviewed: the plan itself, cross-checked against
  `dynamic/bimaterial_wave_operator.{hpp,inl}`, `dynamic/godunov_flux_bimaterial.{hpp,cpp}`,
  `dynamic/godunov_flux.hpp`, `dynamic/wave_operator.{hpp,inl}`,
  `drivers/spatial_dyn_driver.cpp`, `spatial/code/spatial_friction.{hpp,cpp}`.
- Domain context: `CLAUDE.md`, `miniapps/seas/CLAUDE.md`, the drdg3d comparison doc,
  this session's sigma_n-by-depth findings.
- Nature: this reviews a PLAN (no code yet). "Actual/Expected behavior" = what the plan
  specifies vs what it must specify; "Suggested fix" = an edit to the plan or a hard
  requirement the /code-fix or implementer must honor.

## Findings

### [R-001] MODERATE — Phase 1/4 mechanism test: self-side-only metric is not a valid dissipation measure across a contrast

**Category:** BUG

**Description:**
Phases 1 and 4 specify the mechanism test as
`dot(jump, F_upwind_self - F_central_self) > 0` and `-> 0.5*Zp*jump^2` in the
homogeneous limit. The bi-material upwind flux is **two-valued** at an impedance
contrast (`F_self != F_nbr`; see `bimaterial_wave_operator.inl:741-748` where side e1
and e2 use DIFFERENT matrix pairs `mat_e1_*` vs `mat_e2_*`), whereas the central flux
is single-valued. The true face energy-dissipation rate is the sum over BOTH sides
(`-(Q_L·F*_L·nL + Q_R·F*_R·nR)`), not the self-side difference alone. The self-side
quantity has no guaranteed sign when `A_self != A_nbr`, so the assertion `> 0` may fail
spuriously (false alarm) or, if it happens to pass, does not actually demonstrate that
central lacks the dissipation that upwind provides.

**Trigger:** Run the Phase-1/4 test with the TPV31 5 km contrast (Zs 7.99e6 vs 9.38e6).

**Actual behavior (as planned):** asserts a one-sided quantity that may be negative or
meaningless across the contrast.

**Expected behavior:** assert a genuine, sign-definite energy-dissipation measure that
is 0 for central and > 0 for upwind across the contrast.

**Suggested fix (edit the plan's test spec):**
```diff
- mechanism: across the TPV31 5 km contrast, for a normal-velocity jump,
-   `dot(jump, F_upwind_self - F_central_self) > 0` strictly (central discards
-   dissipation), and `-> 0.5*Zp*jump^2` in the homogeneous limit (within 1e-9 rel).
+ mechanism: use the FULL two-sided face dissipation.  Build both per-face upwind
+   matrix pairs (e1,e2 via BuildPerFaceFluxMatricesGlobal called twice, self/nbr
+   swapped) and the single central pair.  For a pure normal-velocity jump Q across
+   the contrast, form F_up_e1,F_up_e2 (ApplyPerFaceFlux) and F_ce, then
+     D = (Q_e1 - Q_e2) . ( 0.5*(F_up_e1 + F_up_e2) - F_ce )
+   Assert D > 0 strictly across the 5 km contrast (upwind dissipates the jump that
+   central does not) and D -> 0.5*Zp*jump^2 in the homogeneous limit (within 1e-9 rel,
+   where F_up_e1 == F_up_e2).  Equivalently, validate at the OPERATOR level: apply the
+   full BimaterialWaveOperator Mult to a normal-velocity perturbation localized at a
+   contrast face and show the discrete energy ||Q||_M^2 decreases with the face on
+   upwind and is ~conserved with it on central.
```

**Test case (intent):**
```cpp
// test_bimaterial_contrast_guard: mechanism
GodunovFlux soft(lam_s,mu_s,rho_s), hard(lam_h,mu_h,rho_h);  // 5km pair
real_t nor[3] = {0,0,1};                                     // z-normal layer face
DenseMatrix upL_e1,upN_e1,upL_e2,upN_e2,ceL,ceN;
BimaterialFlux::BuildPerFaceFluxMatricesGlobal(nor, soft, hard, upL_e1, upN_e1);
BimaterialFlux::BuildPerFaceFluxMatricesGlobal(nor, hard, soft, upL_e2, upN_e2);
BimaterialFlux::BuildPerFaceCentralMatricesGlobal(nor, soft, hard, ceL, ceN);
real_t Qself[NUM_STATE]={0}, Qnbr[NUM_STATE]={0};
Qself[VZ]=1.0;  // pure normal-velocity jump
real_t Fup1[NUM_STATE],Fup2[NUM_STATE],Fce[NUM_STATE];
BimaterialFlux::ApplyPerFaceFlux(upL_e1,upN_e1,Qself,Qnbr,Fup1);
BimaterialFlux::ApplyPerFaceFlux(upL_e2,upN_e2,Qself,Qnbr,Fup2);
BimaterialFlux::ApplyPerFaceFlux(ceL,ceN,Qself,Qnbr,Fce);
double D=0; for(int c=0;c<NUM_STATE;c++) D += (Qself[c]-Qnbr[c])*(0.5*(Fup1[c]+Fup2[c])-Fce[c]);
ASSERT_GT(D, 0.0);                       // central discards dissipation across contrast
```

---

### [R-002] MODERATE — Fallback-to-upwind relies on an unasserted invariant

**Category:** ASSUMPTION

**Description:**
The fix routes reclassified faces to upwind by erasing them from
`central_flux_face_set_`, so `InteriorFaceFlux_`/`SharedInteriorFaceFlux_`
(`bimaterial_wave_operator.inl:741-748,769-771`) read
`per_face_bimaterial_flux_[mesh_face]`. This is populated for all non-fault interior
faces by the ctor (`:50` -> `BuildPerFaceBimaterialFluxMatrices_`, `:307` assigns all
faces, loop builds all non-fault faces) — verified correct TODAY. But the plan never
asserts it. If the upwind build is later changed to skip central-set faces (a plausible
memory optimization), an erased face would dispatch to a default-constructed (0x0)
`DenseMatrix`, hitting the `ApplyPerFaceFlux` size `MFEM_VERIFY`
(`godunov_flux_bimaterial.cpp:398`) — or, worse, silent garbage in a release path.

**Trigger:** any future refactor that makes `per_face_bimaterial_flux_` sparse over the
corridor; latent today.

**Actual behavior (as planned):** silently assumes the upwind matrices exist.

**Expected behavior:** assert the fallback target exists before erasing.

**Suggested fix (add to Phase 2 requirements):**
```diff
  AFTER both loops: `for (int f : reclassified) central_flux_face_set_.erase(f);`
+ Before erasing each face, assert its upwind fallback exists:
+   MFEM_VERIFY(per_face_bimaterial_flux_[f][0][0].Height() == NUM_STATE,
+     "contrast-guard: reclassified face " << f << " has no bi-material upwind "
+     "matrices to fall back to (per_face_bimaterial_flux_ not populated).");
+ (For shared faces assert per_face_bimaterial_flux_[f][0][0] likewise, since the
+  shared upwind path reads side 0 only.)
```

**Test case (intent):**
```cpp
// build bimaterial op with a contrast corridor face, tol=0.05, SetMixedFluxMode(Adjacent)
// assert: every face erased from the central set has per_face_bimaterial_flux_[f][0][0]
//         sized NUM_STATE x NUM_STATE (fallback present) -> dispatch produces finite flux.
```

---

### [R-003] MODERATE — Shared-face neighbour-material convention is unresolved but is a correctness prerequisite (rank-consistency)

**Category:** ASSUMPTION

**Description:**
Phase 2 builds a temporary neighbour `GodunovFlux` from
`shared_face_neighbour_material_.at(mesh_face_idx)` to test the contrast on MPI shared
faces, but the storage convention (`{lambda,mu,rho}` vs `{vp,vs,rho}`) is only flagged
as risk R6, not resolved. `GodunovFlux`'s ctor is `(lambda, mu, rho)`
(`godunov_flux.hpp:43`). If the array actually holds `{vp,vs,rho}`, constructing
`GodunovFlux(arr[0],arr[1],arr[2])` yields nonsense impedances → wrong contrast → wrong
classification. Worse, the LOCAL side uses `FluxForElem_(e1)` (true material) while the
NEIGHBOUR side uses the (mis-converted) array → the two ranks owning the shared face can
reach DIFFERENT verdicts → one rank central, the other upwind on the same face →
**non-conservative cross-rank flux**. This is a correctness bug, not merely a risk.

**Trigger:** any np>1 TPV31 run with the guard enabled and a layer interface on a
partition seam.

**Actual behavior (as planned):** convention guessed at implementation time.

**Expected behavior:** convention read from `ExchangeBiMaterialNeighbours_` and
verified; both ranks classify identically.

**Suggested fix (promote R6 to a Phase-2 blocking step):**
```diff
+ Phase 2 step 0 (BLOCKING, do before coding the shared-face filter): read
+   ExchangeBiMaterialNeighbours_ end-to-end and document whether
+   shared_face_neighbour_material_ stores {lambda,mu,rho} or {vp,vs,rho}.  Build the
+   temporary neighbour GodunovFlux with the matching ctor args (convert vp/vs->lambda/mu
+   via mu=rho*vs^2, lambda=rho*vp^2-2mu if needed).
+ Add a parallel assertion test: the temporary neighbour GodunovFlux's GetZp()/GetZs()
+   equal the peer element's GetZp()/GetZs() to 1e-10 relative (proves the convention +
+   conversion are right, hence rank-consistent classification).
```

**Test case (intent):**
```cpp
// np=2, fixture with a known material on each rank straddling the seam:
//   neighbour-derived GodunovFlux.GetZp() ~= peer GodunovFlux.GetZp()  (rel 1e-10)
//   IsStrongContrast(local, neighbour, tol) identical on both ranks for the shared face
```

---

### [R-004] LOW — `all_continuous` mixed-flux mode is unaddressed by tests/acceptance

**Category:** EDGE_CASE

**Description:**
`mixed_flux` may be `"all_continuous"` (`spatial_friction.cpp:1153`), where
`central_flux_face_set_` covers ~95% of faces. The guard filters the set mode-agnostically,
so under all_continuous on a bi-material mesh it would reclassify EVERY layer-crossing
face — a large set, possibly not intended, and untested. The plan only specifies/derives
behavior and tests for `adjacent`.

**Trigger:** `--mixed-flux all_continuous` + bi-material + tol>=0.

**Suggested fix:**
```diff
+ Phase 2: state that the guard applies to ANY mixed_flux mode (it filters the set),
+   and add a Phase-2 acceptance check at all_continuous on the bi-material fixture
+   (all contrast faces reclassified; uniform faces retained).  If all_continuous is
+   out of scope, say so explicitly and assert mixed_flux != all_continuous when tol>=0.
```

**Test case (intent):** dispatch test with mode=AllContinuous, tol=0.05 -> only
contrast faces leave the central set.

---

### [R-005] LOW — Histogram double-counts shared faces under MPI_Reduce

**Category:** QUALITY

**Description:**
Phase 2 sums histogram bins across ranks via `MPI_Reduce`. A shared face is present in
`central_flux_face_set_` on BOTH owning ranks, so it is counted twice in the global
histogram and "reclassified N/M" count. Diagnostic-only (no correctness impact), but
the printed numbers mislead.

**Suggested fix:**
```diff
+ When tallying the histogram over shared faces, count each shared face on the LOWER-rank
+ owner only (or print "(shared faces double-counted)" next to the totals).
```

---

### [R-006] LOW — Disabled path (tol<0) must add ZERO new collective calls (deadlock safety)

**Category:** BUG

**Description:**
The plan gates the filter on `tol >= 0` but also adds an `MPI_Reduce` for the histogram.
Per memory `R-001` (no rank-conditional collective -> deadlock), the new collective MUST
be called by all ranks or none, and the byte-exact contract is for tol<0. The plan should
state that for tol<0 NO new code runs on ANY rank (no Reduce), and for tol>=0 the Reduce
is collective-uniform (tol is a config/CLI scalar, identical on all ranks, so the gate is
rank-uniform — but assert/ document that tol is broadcast/identical across ranks).

**Trigger:** tol set on some ranks only (shouldn't happen via config, but a CLI typo or
per-rank env could), or histogram Reduce called rank-conditionally.

**Suggested fix:**
```diff
+ Phase 2: the entire new block (filter + histogram + MPI_Reduce) is guarded by a SINGLE
+ `if (mixed_flux_contrast_tol_ >= 0)` evaluated identically on all ranks (tol comes from
+ the same config/CLI on every rank).  For tol<0 the function body is byte-identical to
+ today (no new statements, no collectives).  Add a debug assert that tol is identical
+ across ranks (MPI_Allreduce min==max) when >=0.
```

**Test case (intent):** np=2, tol=-1 -> output byte-identical to pre-change; np=2,
tol=0.05 with both ranks -> no hang, histogram prints once.

---

### [R-007] LOW — Contrast metric `max(cZp,cZs)`: document that it is intentionally conservative

**Category:** ASSUMPTION

**Description:**
The leak channel is the NORMAL (P) impedance Zp (sigma_n). The plan keys on
`max(cZp,cZs)`, which also fires on a shear-only contrast. For TPV31 cp and cs co-jump so
it is moot, but the choice is unstated. Using `max` is the safe (conservative) choice;
just record the rationale so a reviewer/implementer does not "simplify" it to Zp-only.

**Suggested fix:**
```diff
+ Note in the Key-math section: threshold on max(cZp,cZs) deliberately (catches any
+ impedance contrast that a non-dissipative central flux would mishandle, not only the
+ P channel); do NOT reduce to Zp-only.
```

---

## Summary
- Critical issues: 0
- Moderate issues: 3 (R-001 test validity, R-002 unasserted fallback invariant,
  R-003 shared-face convention / rank-consistency)
- Low issues: 4 (R-004 all_continuous, R-005 histogram double-count, R-006 disabled-path
  collective safety, R-007 metric rationale)
- Plan compliance: N/A (plan review — measured against the diagnosis + codebase facts)
- Verdict: PASS WITH FIXES — the approach is sound and the linchpin (upwind fallback
  exists for all corridor faces) is verified, but fix R-001/R-002/R-003 in the plan
  before implementation: the mechanism test must measure two-sided dissipation, the
  fallback invariant must be asserted, and the shared-face material convention must be
  resolved (it is a cross-rank correctness prerequisite, not just a risk).

## Unreviewed Areas
- The exact `ExchangeBiMaterialNeighbours_` storage convention (R-003) — not traced to
  ground truth in this pass; flagged as a blocking step for Phase 2.
- The Frontera A/B numerics (Phase 4) — cannot be evaluated locally (no full-mesh runs);
  acceptance thresholds in the plan are reasonable but unverified until the run returns.
- Whether the 1-element corridor ring is WIDE enough that reclassifying only in-corridor
  contrast faces removes the dominant seed — a physics-correctness uncertainty the plan
  correctly defers to the Frontera A/B; not decidable from the plan or unit tests.

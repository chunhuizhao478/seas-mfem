# Code Review: SAFS slip-law + strong-rate-weakening WIRING (2026-05-28)

## Review Scope
- Plan: `debug_document/spatial_dynamic_rupture_slip_law_srw_plan_2026-05-28.md`
  (+ prior plan review `REVIEW_slip_law_srw_plan_2026-05-28.md`)
- Files reviewed (the implemented diff):
  - `spatial/code/spatial_friction.hpp` (schema: `StateEvolutionKind`, `f_w`/`V_w`, `RateStatePerDOFParams::V_w`)
  - `spatial/code/spatial_friction.cpp` (parser + resolver `V_w` fill)
  - `dynamic/tpv104_substep_iterator.{hpp,cpp}` (additive nuc_callback overload)
  - `dynamic/friction_iterator.hpp` (`SlipLawSRWFrictionIterator`)
  - `dynamic/friction_iterator_factory.cpp` (dispatch + `V_w` injection)
  - `drivers/spatial_dyn_driver.cpp` (banners)
  - `Makefile` (driver + factory-test link)
  - `config/…_srw_Dc010_nuc8km_500m.toml`, both `…_8N_400r_dev_2hr[_restart]_safs.sbatch`
  - tests: `test_friction_iterator_factory.cpp`, `test_spatial_friction_config.cpp`
- Domain context: `CLAUDE.md` (Brent solver; friction extreme-care; np>1 RS path),
  `safs-nucleation-sigma-n-and-lengthscales` memory, the SRW physics in
  `friction/slip_law_srw_psi.hpp`.

Adversarial audit of the implemented wiring. The copied overload was diffed
line-for-line against the original `AdvanceWithSubStepStates` — the only
differences are the two intended deltas (nuc_callback, per-QP `d.b`) plus the
empty-callback and `d.b` finiteness guards; **no transcription error found there.**

## Findings

### [R-001] CRITICAL [friction_iterator_factory.cpp:64 MakeFrictionIterator] — SRW factory aborts on any MPI rank with zero fault DOFs

**Category:** BUG

**Description:**
The slip-SRW branch guards the resolved per-DOF `V_w` with
```cpp
MFEM_VERIFY(rs->V_w.Size() == rs->a.Size() && rs->V_w.Size() > 0, ...);
```
The `&& rs->V_w.Size() > 0` clause is wrong. On a parallel run, many ranks own
only far-field elements and have **zero fault DOFs** (`N == 0`):
`ResolveRateState` then returns `rs.a.Size() == 0` and (for SlipSRW)
`rs.V_w.SetSize(0)` → `rs.V_w.Size() == 0` (spatial_friction.cpp:1378-1382).
`MakeFrictionIterator` is called **unconditionally on every rank**
(`spatial_dyn_driver.cpp:2043` — the iterator is needed for `SetSubSteps`, which
is also unconditional at :2061). So on an `N==0` rank the guard evaluates
`0 == 0 && 0 > 0` → `false` → `MFEM_VERIFY` fails → `MPI_Abort` → the whole job
dies at construction.

The aging branch (`RateStateAgingFrictionIterator`) never touches `rs`, so it
runs fine at `N==0` — this is a **SRW-only regression** that the existing
8N/400r aging jobs do not have. It will abort the production 500 m run
immediately (a localized fault at 400 ranks guarantees many `N==0` ranks). It is
invisible to the np=1 unit tests (F6/F8 use `N==2`).

The correct invariant is `rs->V_w.Size() == rs->a.Size()` alone: that still
catches the "SlipSRW selected but `V_w` not filled" bug when `N > 0`
(`0 == N` is false → abort), while allowing the legitimate `N == 0` rank
(`0 == 0` → build an iterator with empty `V_w_`; `Advance` is never called on it
because `AdvanceADERWithSubStep_Spatial` guards the fault branch with
`if (n_total_fault_qps > 0)`).

**Trigger:** any `--mixed-flux`/serial run at `np > 1` (or even np=1 if a rank
geometry yields no fault QPs) with `state_evolution = "slip_srw"`; specifically
the 8N/400r dev sbatch.

**Actual behavior:** `MFEM_VERIFY` abort on every rank that owns no fault DOFs →
immediate job failure.

**Expected behavior:** ranks with no fault DOFs build the iterator and simply
skip the (skipped-anyway) fault `Advance`.

**Suggested fix:**
```diff
-               MFEM_VERIFY(rs->V_w.Size() == rs->a.Size() && rs->V_w.Size() > 0,
-                           "MakeFrictionIterator: rs->V_w (size " << rs->V_w.Size()
-                           << ") must be resolved to the per-DOF count (rs->a size "
-                           << rs->a.Size() << ") for slip_srw; check that "
-                           "ResolveRateState ran with state_evolution==SlipSRW.");
+               // NOTE: do NOT require Size() > 0 — ranks with no fault DOFs
+               // legitimately have N == 0 (rs.a.Size() == rs.V_w.Size() == 0).
+               // The "==" alone still catches a SlipSRW config whose V_w was
+               // left unfilled (Size()==0 while a.Size()==N>0).
+               MFEM_VERIFY(rs->V_w.Size() == rs->a.Size(),
+                           "MakeFrictionIterator: rs->V_w (size " << rs->V_w.Size()
+                           << ") must equal the per-DOF count (rs->a size "
+                           << rs->a.Size() << ") for slip_srw; ResolveRateState "
+                           "fills V_w to N (possibly 0) when state_evolution==SlipSRW.");
```
Also update test F8's intent comment (it relies on `N>0` + empty `V_w`, which
the `==` check still rejects) and add the `N==0` test below.

**Test case:**
```cpp
// tests/unit/test_friction_iterator_factory.cpp
static void F9_slip_srw_zero_fault_dofs_builds()
{
   // N == 0 rank: rs.a and rs.V_w both Size()==0 (ResolveRateState on an
   // empty fault).  Must BUILD, not abort (regression guard for R-001).
   FaultFaceFlux flux(kRho, kCp, kCs);
   const spatial::SpatialFrictionConfig cfg = MakeSlipSRWConfig(FrictionSolver::V0);
   spatial::RateStatePerDOFParams rs = MakeRsWithVw(/*n=*/0, /*fill_vw=*/true);
   std::unique_ptr<IFrictionIterator> fr = MakeFrictionIterator(cfg, flux, &rs);
   TEST_ASSERT(fr != nullptr && fr->WaveOpLaw() == FaultFrictionLaw::RateAndState,
               "slip_srw with N==0 fault DOFs builds (no abort) — R-001");
}
```
(`MakeRsWithVw(0, true)` must `SetSize(0)` every vector; with the current helper
`fill_vw=true` does `rs.V_w.SetSize(0); rs.V_w = 0.1;` which is a valid no-op on
an empty Vector.)

---

### [R-002] MODERATE [POSSIBLE] [config …_srw_Dc010_nuc8km_500m.toml + both sbatch headers] — the depth-profile "VW→VS arrest at 11 km" is an AGING-law property and does not transfer to SRW

**Category:** ASSUMPTION (physics / misleading documentation; possible stability impact)

**Description:**
The new SRW config and both sbatch headers carry over the aging-law claim that
the depth profile's VW→VS transition (a>b below ~11 km) makes "down-dip rupture
arrest in the VS region above the 16.6 km fault bottom." Under the slip law with
strong rate weakening the steady-state friction is
`f_ss(V) = f_w + (f_LV(V) − f_w)/(1 + (V/V_w)^8)^(1/8)`, which → `f_w = 0.2` at
`V ≫ V_w = 0.1 m/s` **regardless of the sign of (b−a)**. So at seismic slip
rates the nominal-VS region weakens to `f_w` just like the VW core; the
steady-state strengthening that produces aging-law arrest is bypassed once the
front drives `V > V_w`. Down-dip arrest under SRW is therefore stress-dependent
and not guaranteed by the a>b profile — the rupture may propagate further
down-dip than the aging config did, up to (or past) the reflection-contaminated
window. This is the documented np>1 RS "treat as PRELIMINARY" regime, now with a
weaker arrest mechanism.

This is not a code defect (the SRW formula is implemented correctly via
`UpdateStateAnalyticSlipLawSRW`), but the config/sbatch comments assert an arrest
behavior that the chosen law does not robustly provide — a reviewer or future
user will be misled, and the run's down-dip extent must be monitored rather than
assumed bounded.

**Trigger:** any SRW run that nucleates and propagates down-dip past ~11 km.

**Actual behavior:** comments claim aging-style VS arrest; SRW provides
f_w-weakening throughout at high V, so arrest is not assured.

**Expected behavior:** the config/sbatch should state that under SRW the VS
region does not guarantee arrest (f_ss→f_w at V≫V_w), and that down-dip extent is
to be monitored; keep the depth profile for the *direct-effect* (a) structure,
not as an arrest guarantee.

**Suggested fix (doc, in the config `[friction.rate_state.depth_profile]` comment
and both sbatch `Config:`/header blocks):**
```diff
-#    (a-b)(z) <- param_a_minus_b_500m_strongvw.csv (VW->VS transition ~11 km);
+#    (a-b)(z) <- param_a_minus_b_500m_strongvw.csv (VW->VS transition ~11 km).
+#    CAVEAT (SRW): unlike the aging law, the slip-law SRW steady state
+#    f_ss -> f_w at V >> V_w REGARDLESS of (b-a) sign, so the a>b region below
+#    ~11 km does NOT guarantee down-dip arrest the way it does for aging.
+#    Treat down-dip extent as something to MONITOR, not a built-in stop.
```

**Test case:** (physics caveat — covered by a documentation assertion, not a unit
test; the closest mechanical check is a steady-state evaluation)
```cpp
// tests/unit/test_slip_law_srw_brent.cpp (or a PsiSS unit check)
TEST(SlipLawSRW, FssApproachesFwAtHighVEvenWhenAExceedsB)
{
   // VS-region params: a > b (a-b < 0), as below 11 km in the profile.
   const real_t a = 0.026, b = 0.018, V0 = 1e-6, f0 = 0.6, fw = 0.2, Vw = 0.1;
   const real_t fss_hi = /* f_ss via PsiSS_SRW->f_ss at V = 10 m/s */;
   EXPECT_NEAR(fss_hi, fw, 0.05)
       << "SRW weakens to f_w at high V even for a>b (no aging-style VS arrest)";
}
```

---

### [R-003] LOW [tests/unit/test_friction_iterator_factory.cpp] — test gap let R-001 through; no N==0 factory case and no resolver-fills-V_w case

**Category:** QUALITY (missing coverage that masks a CRITICAL bug)

**Description:**
F6/F7/F8 exercise only `N>0` (F6 builds with `N=2`; F8 aborts with `N=2`,
unfilled `V_w`). None exercises `N==0`, which is exactly the production-MPI
condition R-001 fails on. Separately, no test asserts that `ResolveRateState`
actually *fills* `rs.V_w` to size `N` for SlipSRW (and leaves it empty for
aging) — the factory's correctness depends on that, but it is currently only
exercised with a hand-built `rs`, not the resolver's output.

**Trigger:** running the unit suite and believing SRW dispatch is covered.

**Actual behavior:** the suite is green while a CRITICAL `N==0` abort ships.

**Expected behavior:** add (a) the `F9` `N==0` factory test from R-001, and (b) a
resolver test asserting `rs.V_w.Size()==N` (filled `V_w_default`) for SlipSRW and
`rs.V_w.Size()==0` for aging.

**Suggested fix:** add `F9_slip_srw_zero_fault_dofs_builds()` (R-001) and a
resolver case in `test_spatial_friction_resolver.cpp`:
```cpp
// after resolving a slip_srw RateStateBlock over a small mesh fixture:
TEST_ASSERT(rs.V_w.Size() == rs.a.Size() && rs.V_w.Size() > 0,
            "SlipSRW: ResolveRateState fills per-DOF V_w to N");
for (int i = 0; i < rs.V_w.Size(); ++i)
   TEST_ASSERT(rs.V_w(i) == blk.V_w_default, "V_w filled with V_w_default");
// and for an aging block:
TEST_ASSERT(rs_aging.V_w.Size() == 0, "aging: V_w left empty");
```

---

## Resolved during review (no finding)
- **Overload transcription:** the SAFS callback overload
  (`tpv104_substep_iterator.cpp:697-893`) is a faithful copy of the original
  (`:509-684`) — validation block, Σ-deltaT check, memset, sub-step loop,
  slip accumulation, accumulator scaling, and `WriteBackState`-on-last all match;
  the only changes are `nuc_callback(t_sub_end, dt_sub)` (was
  `ApplyNucleationIncremental_TPV104`), `d.b` (was `state_evo_.GetB()`), and the
  two added guards. TPV104 byte-exact path untouched.
- **`SlipLawSRWPsi` ctor mapping** in the adapter
  (`a_default, b_default, V_0_default, f_0_default, f_w_default, V_w_default`)
  matches the ctor `(a, b, V0, f0, muW, V_w_default)` with `muW = f_w`; scalar
  a/b are unused placeholders (ψ uses per-QP d.a/d.b), only V0/f0/muW are read.
- **Resolver `V_w` fill** (spatial_friction.cpp:1378-1382) is after the per-DOF
  loop, gated on SlipSRW, `operator=(real_t)` fills all N; correct for both the
  ParMesh and serial `resolve_rs_impl<MeshT>` instantiations.
- **Driver banners** (`spatial_dyn_driver.cpp:817-826, 847-852, verify-dispatch
  is_srw`) are well-formed and compile; informational only.
- **Equilibrium-ψ seed** reused unchanged; primary branch is law-agnostic; the
  aging-specific fallback is unreachable at V_init=1e-12 (prior R-004, noted).

## Summary
- Critical issues: 1 (R-001)
- Moderate issues: 1 (R-002, POSSIBLE)
- Low issues: 1 (R-003)
- Plan compliance: FULL on the implemented scope; R-001 is an MPI edge case the
  plan's §4.5 V_w guard under-specified; R-002 is a physics caveat the plan §3
  flagged in spirit but the config/sbatch comments overstate.
- Verdict: **FAIL — must fix before submitting on Frontera.** R-001 will abort
  the 8N/400r production job at construction; it is a one-line fix + an N==0
  test. R-002 is a must-do documentation correction (and a run-monitoring note);
  R-003 closes the coverage gap that hid R-001.

## Unreviewed Areas
- The deferred plan §4.7 iterator-physics tests (per-QP-b ψ drive, SRW+Brent
  convergence, adapter↔bare parity) are still unimplemented — out of scope here,
  but R-002's test sketch overlaps the Brent/PsiSS test and should be folded in.
- End-to-end behavior on the actual 500 m mesh (nucleation → propagation →
  down-dip extent) is unverifiable locally (no local mesh runs) — must be watched
  in the first Frontera run, especially per R-002.
- `Makefile` link additions were validated only by a successful local link, not
  by a Frontera build.

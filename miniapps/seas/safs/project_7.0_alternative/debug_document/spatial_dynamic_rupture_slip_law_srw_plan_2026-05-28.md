# Plan — Slip-Law + Strong-Rate-Weakening (TPV104-style FL=103) for the SAFS spatial dynamic-rupture driver

**Date:** 2026-05-28
**Author:** implementation planning (`/code-plan`)
**Status:** READY — §9 decisions resolved + review findings R-001…R-005 folded in
(2026-05-28); cleared to implement. Review: `REVIEW_slip_law_srw_plan_2026-05-28.md`.
**Scope:** add a TPV104-style slip-law + strong-rate-weakening (SRW / "Fast
Velocity Weakening", SCEC FL=103) friction model to `seas_spatial_dyn_driver`,
then repoint the two dev sbatch jobs at it on the 500 m `_triq` mesh with
Dc=0.10, a=0.0127, b=0.0261.

Target jobs (to be re-pointed, not duplicated):
- `jobs/safs/spatial_dyn_ratestate_depthprofile_normalcap_pureupwind_triq_8N_400r_dev_2hr_safs.sbatch`
- `jobs/safs/spatial_dyn_ratestate_depthprofile_normalcap_pureupwind_triq_8N_400r_dev_2hr_restart_safs.sbatch`

---

## 1. Objective

Switch the SAFS dynamic-rupture friction from the **aging rate-and-state law**
(current) to the **slip law with strong rate weakening** (TPV104 FL=103), keeping
the nucleation-proven length scales (Dc=0.10 m, a=0.0127, b=0.0261, 8 km-diameter
gradual_overstress patch, 500 m `_triq` mesh). The dynamic stress drop will then
be driven by the strong-rate-weakening term (`f_w`, `V_w`) rather than by the
aging-law `(b−a)` slowness effect alone.

---

## 2. Decisive architecture findings (verified in code)

| Question | Finding | Consequence |
|---|---|---|
| What friction laws does the spatial driver support? | `FrictionLawKind { SlipWeakening, RateState }`; `RateState` == **aging** only, via `RateStateAgingFrictionIterator` wrapping `Tpv102SubStepIterator` + `AgingLawPsi` (`dynamic/friction_iterator.hpp:155`). | SRW is **not** a config flip — it needs new code. |
| Does the SRW law already exist? | Yes: `friction/slip_law_srw_psi.hpp` (`SlipLawSRWPsi`, `UpdateStateAnalyticSlipLawSRW`, `PsiSS_SRW`) — the FL=103 ψ-space law. Wired **only** into `seas_tpv104_driver`. | Reuse it; do not reimplement the physics. |
| Is the SRW sub-step iterator hardcoded to TPV104? | **No.** `Tpv104SubStepIterator(FaultFaceFlux&, const SlipLawSRWPsi&)` takes the friction scalars `a,b,V0,f0,muW` from the passed `SlipLawSRWPsi` (any values) and per-QP `V_w` as a per-call side-channel (`dynamic/tpv104_substep_iterator.hpp:64–198`). | SAFS values (Dc, a, b, f_w, V_w) flow through it unchanged. |
| Is the extension point anticipated? | **Yes.** `friction_iterator_factory.cpp:48–53` documents: *"The slip-law SRW variant lands in Phase 5, gated on a future `[friction.rate_state] state_evolution` field … branch on it"*. `rs` (per-DOF params) is already threaded into the factory but unused. | The design was scaffolded for exactly this. |
| Is the force solve law-specific? | **No.** Both aging and SRW route the per-QP friction solve through `FaultFaceFlux::ComputeStageState → FrictionSolver::Solve`; the regularized strength `τ = σ_n·a·asinh((V/2V0)·exp(ψ/a))` is identical. **Only the state-evolution (ψ_ss and dψ/dt) differs.** | No change to the friction solver. The aging RS path passes `Method::Brent` (CLAUDE.md); SRW must do the same. |
| Is the equilibrium-ψ seed law-specific? | **Mostly no.** `SeedEquilibriumPsi_RS` (`dynamic/spatial_setup.hpp:397`) primary branch seeds ψ from the asinh strength inversion at `V_init` — law-agnostic. Its `tau_eff≤0` fallback uses the AGING ψ_ss (R-004), unreachable at V_init=1e-12 (eta·V_init≈4.6e-6 Pa ≪ tau0). | Reused **unchanged**; do not raise V_init into the damping-dominated regime without an SRW ψ_ss fallback. |
| How does the spatial driver apply `gradual_overstress` nucleation? | Per-sub-step `nuc_callback(t_sub_end, dt_sub)` passed into `IFrictionIterator::Advance` (`spatial_dyn_driver.cpp:2073,503`); it writes `tau{1,2}_nuc` before each sub-step's friction solve. | The SRW iterator's `Advance` must honor this callback. |
| Does `Tpv104SubStepIterator` accept a nuc_callback? | **No** — it hard-codes `ApplyNucleationIncremental_TPV104` in the sub-step loop (`tpv104_substep_iterator.cpp:295,627`). `Tpv102SubStepIterator` **does** have a callback overload (`tpv102_substep_iterator.cpp:377–509`, the SAFS aging path uses it). | Add a parallel additive callback overload to the TPV104 iterator (see §4.4). |

**Net:** the difference between today (aging) and the target (SRW) is *only* the
state-evolution law object and a per-QP `V_w` field. The force solve, the
equilibrium seed, the ADER predictor, the mixed-flux/upwind dispatch, the
nucleation accumulator, the I/O, and the restart format are all unchanged.

---

## 3. Physics: the SRW slip law (what we are turning on)

Regularized strength (unchanged, law-agnostic):
```
τ = σ_n · a · asinh( (V / 2V0) · exp(ψ/a) )
```
State evolution (slip law, SRW steady state):
```
dψ/dt   = -(V/L) · (ψ - ψ_ss(V))            # slip law (L = Dc)
ψ_ss(V) = a · ln( (2V0/V) · sinh(f_ss(V)/a) )
f_ss(V) = f_w + (f_LV(V) - f_w) / (1 + (V/V_w)^8)^(1/8)
f_LV(V) = max(0, f0 - (b - a)·ln(V/V0))
```
- `f_w` = fully-weakened friction (TPV104: 0.2).
- `V_w` = weakening velocity (TPV104 VW core: 0.1 m/s).
- The analytic per-sub-step update is `UpdateStateAnalyticSlipLawSRW(...)`
  (`slip_law_srw_psi.hpp:92`), used exactly as the TPV104 driver uses it.

The aging-law length scales in the current config (L_b=2.51 km, L_inf=6.06 km,
8 km patch = 1.32× supercritical) governed nucleation under **aging**. Under the
**slip law**, spontaneous post-forcing propagation is governed differently
(slip-law nucleation has no finite L_∞ blow-up; the forced `gradual_overstress`
patch + SRW dynamic weakening drive the rupture). Keeping the proven a/b/Dc means
nucleation forcing is identical; the SRW term sets the dynamic stress drop. This
is consistent with the user's "Dc/a/b proven for nucleation; add SRW for
propagation" intent. **This physics shift is the main reason to validate the
first run as a stability/propagation smoke test, not a quantitative benchmark.**

---

## 4. Code change surface (file-by-file)

### 4.1 `spatial/code/spatial_friction.hpp` — config schema
Extend `RateStateBlock` (currently ends at the `depth_profile` member):
```cpp
// State-evolution selector: "aging" (default, unchanged) or "slip_srw".
std::string state_evolution = "aging";
// SRW (FL=103) globals — consumed only when state_evolution == "slip_srw".
real_t f_w_default   = 0.2;     // fully-weakened friction (TPV104 = 0.2)
real_t V_w_default   = 0.1;     // weakening velocity [m/s] (TPV104 core = 0.1)
```
Extend `RateStatePerDOFParams`:
```cpp
Vector V_w;   // per-DOF weakening velocity (SRW only; empty for aging)
```
Add an `enum class StateEvolutionKind { Aging, SlipSRW };` OR keep the string and
branch on it in the factory — recommend a small enum resolved at parse time and
stored on the block to avoid string compares on the hot path. (Decision: §9.)

### 4.2 `spatial/code/spatial_friction.cpp` — parser + resolver
- In the `[friction.rate_state]` parse block (near the existing `a_default` /
  `b_default` reads): parse `state_evolution`, `f_w`, `V_w`. Validate:
  - `state_evolution ∈ {"aging","slip_srw"}` (MFEM_ABORT otherwise, mirroring the
    `[meta].law` validator at line 650).
  - When `slip_srw`: `0 < f_w < f0`, `V_w > 0` (loud abort with the offending
    value).
- In `ResolveRateState` (both ParMesh and serial overloads): populate `out.V_w`
  to size `ndof`. For the constant-`V_w` case, fill with `blk.V_w_default`.
  (Spatially-varying `V_w` via `SpatialRule`/depth-profile is a later extension;
  flag in §9.)
- Aging path: leave `out.V_w` empty (size 0) so nothing downstream reads it.

### 4.3 `dynamic/friction_iterator.hpp` — new adapter
Add `SlipLawSRWFrictionIterator : public IFrictionIterator`, modeled on
`RateStateAgingFrictionIterator` but:
- Owns a `SlipLawSRWPsi state_evo_` constructed as
  `SlipLawSRWPsi(a=a_default /*placeholder*/, b=b_default /*placeholder*/, V0,
  f0, muW=f_w, V_w_default)`. **R-001/R-005:** in the SAFS path the iterator
  sources `a`, `b`, `L=Dc` PER-QP from `DOFData` (`d.a`, `d.b`, `d.Dc`), so the
  scalar `a` AND `b` passed here are **placeholders, unused** — only `V0`, `f0`,
  `muW` (the non-virtual `Get*()` getters) are read. Calling
  `state_evo_.SetProductionMode()` is a **harmless safety net** (the iterator
  evolves ψ via the free `UpdateStateAnalyticSlipLawSRW` + non-virtual getters,
  never the base virtuals production mode guards); set it anyway so a future edit
  that routes through a base virtual fails loudly rather than silently.
- Owns a `Tpv104SubStepIterator it_(flux, state_evo_)`.
- Holds a per-QP `std::vector<real_t> V_w_;` and a setter
  `void SetVw(std::vector<real_t> v) { V_w_ = std::move(v); }` (the
  `IFrictionIterator::Advance` signature carries no `V_w`, so the driver injects
  it once before the time loop).
- `Advance(...)` forwards to the new **callback overload** of
  `Tpv104SubStepIterator::AdvanceWithSubStepStates(dof_data, fault_coords, V_w_,
  Q_pw_plus, Q_pw_minus, dt, t_start, I_imp_plus, I_imp_minus,
  FrictionSolver::Method::Brent, nuc_callback)`.
  - **Pass `Method::Brent`** (CLAUDE.md: Brent, not Newton — matches the aging
    SAFS adapter, NOT the TPV104 default of `NewtonRaphsonStable`).
  - Reject an empty `nuc_callback` (mirror the `LswFrictionIterator` /
    `Tpv102` R-013 guard).
  - MFEM_VERIFY `V_w_.size() == dof_data.size()` (loud, before the call).
- Delete copy/move (it_ holds `const SlipLawSRWPsi&` → same dangling-ref footgun
  guarded in `RateStateAgingFrictionIterator`; member-init order: `state_evo_`
  before `it_`).
- `WaveOpLaw()` returns `FaultFrictionLaw::RateAndState` (same wave-op flux law
  as aging — the wave operator only distinguishes LSW vs RS).

### 4.4 `dynamic/tpv104_substep_iterator.{hpp,cpp}` — additive callback overload
Add an overload of `AdvanceWithSubStepStates` taking a trailing
`const std::function<void(real_t,real_t)>& nuc_callback`, copied from the
existing body with **exactly TWO changes** (R-001 — it is NOT "identical except
the nucleation line"):
1. the hard-coded `ApplyNucleationIncremental_TPV104(...)` call inside the
   sub-step loop (`tpv104_substep_iterator.cpp:627`) is replaced by
   `nuc_callback(t_sub_end, dt_sub)` — mirrors the proven
   `Tpv102SubStepIterator` R-013 pattern; AND
2. **the ψ-update sources `b` from the per-QP `d.b`, NOT the scalar
   `state_evo_.GetB()`** (the update is at `tpv104_substep_iterator.cpp:416`,
   and the SeisSol variant at `:658`):
   ```diff
   - state_evo_.GetB(),
   + d.b,   // per-QP depth-profile b (R-001; scalar GetB() would flatten b(z))
   ```
   **REQUIRED for decision #3** (depth-varying b(z) → VW→VS arrest at 11 km).
   `InitializeFaultDOFs_Spatial_RS` sets `d.b = rs.b(i)` per-DOF
   (`spatial_setup.hpp:341`). Without this the iterator would evolve every QP
   with one constant b — silently discarding the depth profile and mixing a
   scalar b with per-QP a in `f_LV = max(0, f0−(b−a)·ln(V/V0))`.
- Add a precondition in the overload:
  `MFEM_VERIFY(std::isfinite(d.b) && d.b > 0.0, ...)` (DOFData.b defaults to NaN).
- **Existing signatures untouched** → the TPV104 driver path keeps `GetB()` and
  its byte-exact regression is preserved (the new overload is only reachable from
  the SAFS adapter).
- Add the same empty-callback guard as Tpv102 (R-013).
- **RESOLVED (review):** the SRW sub-step friction pipeline reads **both**
  `tau1_nuc` and `tau2_nuc` — `ComputeStageState` folds both into the total and
  corrected traction (`fault_face_flux.cpp:183–184, 320–321`), so the SAFS
  `gradual_overstress` dip+strike injection works through this iterator unchanged
  (moot for the current config anyway: `delta_tau_dip_pa = 0.0`).

### 4.5 `dynamic/friction_iterator_factory.cpp` — dispatch
- Branch on `cfg.rate_state->state_evolution`:
  - `"aging"` → `RateStateAgingFrictionIterator` (unchanged).
  - `"slip_srw"` → build `SlipLawSRWFrictionIterator` from
    `*cfg.rate_state` (scalars) and **set its per-QP `V_w` from `rs->V_w`**
    (the `rs` argument is finally consumed — remove the `(void) rs;`).
    MFEM_VERIFY `rs != nullptr` and `rs->V_w.Size() == <ndof>` on this branch.
- Keep the `V_0_default == FrictionSolver::V0` guard for both RS branches.

### 4.6 `drivers/spatial_dyn_driver.cpp` — wiring
- After `MakeFrictionIterator(...)` (line 2043), the SRW path needs its per-QP
  `V_w` set. Two clean options (decision §9):
  - (a) Resolve `V_w` inside the factory from `rs` (preferred — keeps the driver
    law-agnostic; the factory already receives `rs`). Driver unchanged here.
  - (b) `dynamic_cast<SlipLawSRWFrictionIterator*>` in the driver and call
    `SetVw(...)`. Avoid — leaks the concrete type into the driver.
  → Recommend (a): the factory builds *and fully configures* the iterator.
- Banner (`spatial_dyn_driver.cpp:1606–1610`): extend the RS line to print
  "RateAndState (slip-law SRW)" / "SlipLawSRWFrictionIterator (Tpv104 SRW,
  Brent)" when `state_evolution == "slip_srw"`. Echo `f_w`, `V_w`.
- `SeedEquilibriumPsi_RS` call (line 1392): **unchanged** — reused as-is. R-004
  caveat: only the primary (`tau_eff>0`) branch is law-agnostic; the
  damping-dominated (`tau_eff≤0`) fallback seeds the AGING ψ_ss, but is
  unreachable at V_init=1e-12 (eta·V_init≈4.6e-6 Pa ≪ tau0). Do NOT raise V_init
  into the damping-dominated regime without adding an SRW ψ_ss fallback.
- `--print-derived` (`spatial_print_derived.cpp`): optionally add the SRW
  steady-state friction `f_ss(V_init)` / `f_w` / `V_w` to the derived dump.
  (Nice-to-have; not required for the run.)

### 4.7 Tests (mandatory — CLAUDE.md: friction is extreme-care)
- `tests/unit/test_friction_iterator_factory.cpp`: add a `state_evolution =
  "slip_srw"` config case → asserts the factory returns an iterator with
  `WaveOpLaw()==RateAndState` and (new) does not abort.
- New `tests/unit/test_slip_law_srw_iterator_parity.cpp` — **TWO tests (R-002)**:
  1. **Constant-b parity:** drive `SlipLawSRWFrictionIterator` and a bare
     `Tpv104SubStepIterator` with a CONSTANT `b` (so `d.b == GetB()`), a no-op
     nuc_callback, and **`Method::Brent` passed explicitly to BOTH sides**;
     assert bit-identical DOFData after one macro-step. Proves the adapter
     plumbing adds no perturbation. (With the R-001 fix the adapter reads `d.b`
     while the bare iterator reads `GetB()`, so this parity holds **only** for
     constant b — do NOT use a depth-varying fixture here, and do NOT let the
     bare iterator default to `NewtonRaphsonStable`.)
  2. **Per-QP b honored (catches R-001):** two QPs with identical `a, Dc, V, ψ,
     V_w` but DIFFERENT `d.b` (e.g. 0.0261 at the patch vs 0.005 near the VW→VS
     edge) must produce DIFFERENT ψ after one macro-step. If b were flattened to
     a scalar, ψ would be identical.
- `tests/unit/test_advance_interface_compiles.cpp`: extend so the new adapter is
  exercised through the `IFrictionIterator` interface.
- **SRW + Brent convergence (R-003):** `tests/unit/test_slip_law_srw_brent.cpp`
  — for representative `(ψ, σ_n, a)` points, build the stress from the analytic
  strength at a known `V*` and assert the Brent solve recovers `V*` to <1e-6
  relative. The TPV104 traces were validated with `NewtonRaphsonStable`; this
  guards the SAFS-mandated `Method::Brent` on the SRW envelope (CLAUDE.md).
- Run `make test` (baseline must stay green); run the TPV104 smoke/regression to
  confirm the additive overload left the oracle byte-exact.

---

## 5. Config + jobs changes

### 5.1 New config TOML
Create
`config/spatial_friction_rate_state_safs_projected_stress_srw_Dc010_nuc8km_500m.toml`
by copying the existing
`..._depthprofile_Dc010_nuc8km_500m.toml` and changing only:
- `[meta].description` → note "slip-law SRW (FL=103)".
- `[numerics].mixed_flux = "none"` (the existing file says `"adjacent"` but the
  sbatch overrides with `--mixed-flux none`; set the config self-consistent to
  the confirmed pure-upwind production stance — see the mixed-flux blow-up
  analysis 2026-05-26).
- `[friction.rate_state]`: add
  ```toml
  state_evolution = "slip_srw"
  f_w             = 0.2        # decision §9
  V_w             = 0.1        # m/s; decision §9 (constant vs spatially varying)
  ```
- Keep the depth-profile CSVs (`param_a_500m_strongvw.csv`,
  `param_a_minus_b_500m_strongvw.csv`) so a(patch)=0.0127, b(patch)=0.0261 and
  the VW→VS arrest at 11 km are preserved. (Decision §9: keep profile vs constant
  a/b.)

### 5.2 The two sbatch files (edit in place)
Both currently default to the **250 m** mesh + Dc=0.05 config + **presplit**.
Change in each:
- `CONFIG_TOML` default → the new SRW 500 m config (§5.1).
- `TRIQ_MESH` → `experimental_mesh_refinement/safs_fault_box_nwcut_500m_lcfar3000_z0embed_triq.msh`.
- `SAFS_MESH_MODE` default → **`serial`** (500 m `_triq` = 1.16 M tets, serial is
  OOM-safe; presplit not needed — drop the 250 m presplit-prefix default).
- CSV pre-flight list → `param_a_500m_strongvw.csv` / `param_a_minus_b_500m_strongvw.csv`.
- Banner / comment text: "slip-law SRW (FL=103); Dc=0.10, a=0.0127, b=0.0261,
  f_w=…, V_w=…; 500 m _triq mesh".
- Restart job: the 8N/400r rank count is unchanged, so the 400-file checkpoint
  check stays; update the stale `RESTART_FROM` default + the comment that warns
  about the 250 m/Dc=0.05 mismatch (now it must come from an SRW-500 m dev run).
- Keep dev queue 8N/400r/2 h, `--mixed-flux none`, the `sigma_n_strength_floor_pa`
  pre-flight grep (the new config keeps the 10 MPa floor).
- **Optional rename** (flag, do not do silently): the filenames say
  "ratestate_depthprofile"; an SRW-accurate name would be
  `spatial_dyn_srw_sliplaw_..._triq_500m_8N_400r_dev_2hr_safs.sbatch`. Decision §9.

---

## 6. Risk & extreme-care notes (CLAUDE.md compliance)

- **Friction is extreme-care.** The new code is *additive*: a new adapter, a new
  additive overload, new config fields. No existing friction signature, no
  `dieterich_ruina.hpp`, no `state_evolution.hpp`, no `fault_face_flux.cpp`, no
  `wave_operator.inl` is modified.
- **Brent, not Newton** — the adapter passes `Method::Brent` (CLAUDE.md), unlike
  the TPV104 default.
- **TPV104 byte-exact oracle** preserved: the new callback overload leaves the
  existing two `AdvanceWithSubStepStates` / `Advance` entry points unchanged; the
  parity unit test (§4.7) and the TPV104 regression guard this.
- **No reverts.** This does not touch any prior fix; it adds a parallel law.
- **Physics caveat:** switching aging→slip-law changes nucleation/propagation
  dynamics (no aging L_∞ instability). First run is a propagation/stability smoke
  test. The reflection-free window (~6.94 s, R-025) and the 100 s cap-monitor
  framing from the existing sbatch headers still apply.

---

## 7. Why not the alternatives

- **Switch the sbatch to `seas_tpv104_driver`:** infeasible. That driver hardcodes
  `TPV104Params` material (ρ=2670, c_s=3464), the TPV104 boxcar a(x,z)/V_w(x,z)
  geometry, a planar y=0 fault, and TPV104 nucleation — none compatible with the
  SAFS projected-stress tensor, curved CFM fault, 500 m mesh, or depth profile.
- **Write a fresh SAFS SRW sub-step iterator** (no touch to the TPV104 file):
  duplicates ~200 lines of proven friction-pipeline/ψ-update/imposed-state code
  that must be kept in sync. Higher long-term divergence risk than the additive
  overload. (Kept as a fallback if §4.4's `tau1_nuc` check reveals the TPV104
  pipeline is too strike-slip-specialized to generalize cleanly.)

---

## 8. Implementation order (for the implementer)

1. Schema: `RateStateBlock` + `RateStatePerDOFParams` fields (§4.1). Build.
2. Parser + resolver (§4.2). Add a parse unit test. Build + `make test`.
3. Additive callback overload on `Tpv104SubStepIterator` (§4.4) + verify
   `tau1_nuc`/`tau2_nuc` both read. Run TPV104 regression → must be byte-exact.
4. New adapter `SlipLawSRWFrictionIterator` (§4.3). Tests (§4.7): constant-b
   parity, per-QP-b (R-001 guard), and SRW+Brent convergence (R-003).
5. Factory dispatch + V_w wiring (§4.5). Factory unit test (§4.7).
6. Driver banner (§4.6). Build the driver.
7. New config TOML (§5.1).
8. Edit the two sbatch files (§5.2).
9. Full `make test`; dry-run the driver locally on a small mesh with the SRW
   config (`--tfinal` tiny) to confirm it nucleates and steps without abort.

---

## 9. DECISIONS — RESOLVED 2026-05-28

1. **`f_w` = 0.2** (TPV104 FL=103 default). Moderate dynamic stress drop.
2. **`V_w` = 0.1 m/s, CONSTANT** over the whole fault (no spatial-`V_w`
   resolution code for this first SRW run). `ResolveRateState` fills `rs.V_w`
   with the scalar `V_w_default`.
3. **Keep the depth-profile a(z)/b(z)** — `param_a_500m_strongvw.csv` /
   `param_a_minus_b_500m_strongvw.csv`: a(patch)=0.0127, b(patch)=0.0261, VW→VS
   arrest at 11 km. Unchanged from the existing Dc010 config.
4. **Reuse `Tpv104SubStepIterator` + additive nuc_callback overload** (§4.4).
   **R-001:** because decision #3 keeps depth-varying b(z), the overload MUST also
   source `b` from per-QP `d.b` (not the scalar `state_evo_.GetB()`) — see §4.4.
   The fresh-iterator fallback (§7) is unnecessary: the `tau1_nuc` concern is
   resolved (`ComputeStageState` folds both nuc components).
5. **Keep the existing sbatch filenames** — repurpose the two `*_ratestate_
   depthprofile_*` files in place (content + banners updated to SRW). No rename.
6. **`mixed_flux = "none"`** in the new config (matches the sbatch `--mixed-flux
   none` override and the confirmed pure-upwind production stance).

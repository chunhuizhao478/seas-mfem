# Implementation Plan: San Andreas Fault — Rate-and-State SeisSol -> ParaView Verification Activity (two groups)

> FOLDER NOTE: requested as `safs_seisol_v2_2_0_GroupActivity` (singular); saved in
> `safs_seisol_v2_2_0_GroupActivities` (plural), which already holds a full copy of the case
> (mesh, materials, parameter files). Rename the folder if you want the singular name; nothing
> in this plan depends on it.

## Overview

A two-group, hands-on activity on the San Andreas Fault (SAFS) using **rate-and-state friction
WITHOUT strong rate weakening**: **Group 1 = aging law (SeisSol FL=3)**, **Group 2 = slip law
(SeisSol FL=4)**. Both groups run the SAME three tasks: **[1]** with constant material AND with
the CVM material; **[2]** compute rupture speed (Vr/RT), seismic moment (M0), seismic potency
(P0) in ParaView and verify against SeisSol's own outputs; **[3]** build the surface PGV map and
correlate it with the fault rupture velocity. Pedagogy: do the **constant-material** case first
(so M0 = mu*P0 is exactly verifiable and participants learn to trust their ParaView pipeline
against SeisSol ground truth), THEN repeat on the **CVM** material.

This is a post-processing/pedagogy deliverable plus SeisSol input decks that MUST pass the Phase-6
smoke test (rupture nucleates AND propagates) before distribution -- plain RS without strong
weakening is prone to confined, non-propagating events; the decks ship with enlarged nucleation
and may still need the Risk-Assessment tuning. It does NOT
modify the solver. The friction model is the only physics change vs the existing `..._LSW` case:
LSW (FL=16) is replaced by rate-and-state (FL=3 / FL=4), derived from the existing SAFS RSSRW
(FL=103) setup with the strong-velocity-weakening fields removed.

## Established facts (verified in SeisSol `v1.3.1-1760-g49bdd63e4` at /Users/chunhuizhao/projects/SeisSol)

### Friction-law selection
- `DRParameters.h:33-35`: `RateAndStateAgingLaw = 3`, `RateAndStateSlipLaw = 4`,
  `RateAndStateFastVelocityWeakening = 103`. **Aging = FL 3, slip = FL 4.** FL=103 is the
  strong/fast velocity weakening we are REMOVING.
- `Factory.cpp:64,93`: FL=3 -> `AgingLaw<NoTP>`, FL=4 -> `SlipLaw<NoTP>` (no thermal pressurization).
  Both use the SAME initializer (`RateAndStateInitializer`) and the SAME parameters — the only
  difference is the state-evolution ODE. So one fault yaml serves both groups; only the par-file
  `FL` line changes.

### Rate-and-state parameters
- Global (in `&DynamicRupture`, required when RS), `DRParameters.cpp:100-104`:
  `rs_f0`, `rs_b`, `rs_sr0`, `rs_inisliprate1`, `rs_inisliprate2`.
- `rs_muw`/`mu_w` is read ONLY for FL=103 (`DRParameters.cpp:106-107`) — OMIT for FL=3/4.
- Per-fault-point (yaml), `RateAndStateInitializer.cpp:126-138`: `rs_a`, `rs_sl0` (alt `RS_sl0`);
  and `insertIfPresent("rs_b", ...)` => **spatial b(z) from the yaml IS honored for FL=3/4**
  (par-file `RS_b` is the fallback). `rs_srW` is read ONLY by the Fast-Velocity initializer — OMIT.
- Initial RS state is auto-computed from the initial slip rate + projected stress
  (`RateAndStateInitializer.cpp:102-118`); do NOT set state directly.

### Nucleation (friction-law-independent)
- `BaseDRInitializer.cpp:47-138`: nucleation stress is read generically and rotated to the fault
  CS. Traction-parameterized `Tnuc_s/Tnuc_d/Tnuc_n` (what SAFS uses) works for any FL, including
  FL=3/4. It is ramped in time by SeisSol's smoothStep over `[s_0, s_0+t_0]`.

### Quantities SeisSol auto-computes (the verification ground truth) — unchanged by FL
- Potency `P0 = integral(slip dA)` (`EnergyOutput.cpp:312-319`); Moment `M0 = integral(mu*slip dA)`
  (`:307-320`), `mu = 2*muP*muM/(muP+muM)`, `muSide = rho*Vs^2` => **homogeneous medium: M0 = 3.2e10 * P0 exactly.**
- `ASl` accumulated slip = `integral(|slip rate| dt)` (`FrictionSolverCommon.h:603`).
- `RT` = first time `|slip rate| > 1e-3 m/s` (`:480-483`); `Vr = 1/|grad RT|`
  (degree N-1 polynomial gradient, `ReceiverBasedOutput.cpp:463-523`; `Vr=0` where `RT==0`).
- `PSR` = running max of `|slip rate|` (`:509`).
- Energy CSV is long/tidy: `time,variable,simulation_index,measurement`; ground-truth rows are
  `variable == seismic_moment` and `variable == potency`.
- Free surface (`FreeSurfaceWriterExecutor.cpp:75`) writes `v1,v2,v3,u1,u2,u3`. **SeisSol does
  NOT output PGV** — PGV is a ParaView-computed product (task 3), not a verify-against-SeisSol item.

## Constraints
- **Per-job download < 1 GB** (QuakeWorx; no cluster shell — all post-processing in ParaView on a
  laptop). Each of the 4 runs is its own job and individually clears 1 GB (~0.47 GB double).
  VERIFY whether the 1 GB limit is per-job (assumed here) or a per-user/per-project quota; if the
  latter, delete each run's download before fetching the next, or coarsen intervals (see levers).
- **Mesh fixed:** `safs_mesh.puml.h5` — 121,316 fault facets (BC 3), 34,096 free-surface triangles
  (BC 1). Unchanged across all runs.
- **Constant material (exact M0):** `[rho=2670, mu=3.2e10, lambda=3.2e10]` (Poisson 0.25,
  Vs=3462 m/s, Vp~=5996 m/s) from `archive/safs_material.yaml`. CVM = `safs_material_cvm.yaml`
  (+ `safs_material_cvm.nc`, ASAGI).
- **Surface output MUST use the legacy compact path** (`surfacevtkorder = -1`): geometry once,
  one growing `.h5`, cell-centered. The heavy `surfacevtkorder = 2` path (used by the RSSRW deck)
  is multi-GB and is forbidden here.
- **Fault output legacy path**, `refinement = 0` (1 value/facet) to fit the budget (Phase 2).
- **No solver changes.** Only input decks + docs.
- **Plain-text math only in all docs** (no LaTeX); PDFs are generated with pandoc+xelatex and the
  source must stay ASCII (no Unicode math glyphs) so the PDF renders cleanly.

## Phase 1: Shared rate-and-state fault model (no strong rate weakening)

### Goal
A single fault yaml exists that drives BOTH the aging (FL=3) and slip (FL=4) laws on the SAFS
curved fault, plus the constant material file, ready to reference from the decks.

### Files to Create
- `safs_fault_rs.yaml` — RS fault model (stress + a(z),b(z) + Dc + Tnuc nucleation), NO `rs_srW`.
- `safs_material.yaml` — copy of `archive/safs_material.yaml` (homogeneous medium) into case root.

### Files to Reuse (unchanged)
- `safs_initial_stress.yaml` (regional effective stress tensor, P_p baked in), `safs_mesh.puml.h5`,
  `safs_material_cvm.yaml` + `safs_material_cvm.nc`.

### Detailed Requirements
1. `cp archive/safs_material.yaml safs_material.yaml` (verify values `rho 2670.0, mu 3.2e10, lambda 3.2e10`).
2. Create `safs_fault_rs.yaml` with EXACTLY this content (derived from `safs_seisol_v2_1_0_RSSRW/safs_fault.yaml`
   by removing `rs_srW` from the LuaMap; everything else — stress, a(z)/b(z), Dc, Tnuc — retained):

```yaml
!Switch
# SAFS rate-and-state fault model WITHOUT strong velocity weakening.
# Use with FL=3 (aging) OR FL=4 (slip): this yaml is identical for both laws;
# only parameters.par's FL line differs. Derived from safs_seisol_v2_1_0_RSSRW
# by REMOVING the strong-rate-weakening field rs_srW (and RS_muW from the deck).
# Frame: x=East, y=North, z=Up (UTM 11N). Compression NEGATIVE, right-lateral POSITIVE.

# --- Initial effective stress tensor (regional Hickman-Zoback, P_p=20 MPa baked in) ---
[s_xx, s_yy, s_zz, s_xy, s_yz, s_xz]: !Include safs_initial_stress.yaml

# --- Spatial RS parameters a(z), b(z): Allison & Dunham (2021) Fig 3a a(d);
#     a-b: VW (-0.004) for d<=8 km, ramp to neutral by 12 km, strong VS (+0.015) by 16 km
#     (deep arrest barrier). d = depth km = -z/1000.
[rs_a, rs_b]: !LuaMap
  returns: [rs_a, rs_b]
  function: |
    function f (x)
      local d = -x["z"] / 1000.0
      if (d < 0.0) then d = 0.0 end
      local a
      if (d <= 14.892395982783356) then
        a = 0.0062769230769231
            + (0.0255076923076923 - 0.0062769230769231) * d / 14.892395982783356
      else
        a = 0.0255076923076923
            + (0.14521923076923074 - 0.0255076923076923)
            * (d - 14.892395982783356) / (52.6829268292683 - 14.892395982783356)
      end
      local amb
      if (d <= 8.0) then
        amb = -0.004
      elseif (d <= 12.0) then
        amb = -0.004 + 0.004 * (d - 8.0) / 4.0
      else
        amb = 0.015 * math.min((d - 12.0) / 4.0, 1.0)
      end
      return { rs_a = a, rs_b = a - amb }
    end

# --- State-evolution distance Dc (rs_sl0) ---
[rs_sl0]: !ConstantMap
  map:
    rs_sl0: 0.10

# --- Nucleation: Gaussian strike overstress at the hypocenter, smoothStep over [0, t_0] ---
# Enlarged vs the SRW source (R 6000->8000, dtau 24->30 MPa): plain RS (FL=3/4) has a large
# nucleation length L_nuc ~ 9-12.5 km, so a bigger forced patch is needed to PROPAGATE.
# CONFIRM by the Phase-6 smoke test; tune further per Risk Assessment if it still arrests.
[Tnuc_s]: !LuaMap
  returns: [Tnuc_s]
  function: |
    function f (x)
      local dx = x["x"] - 606971.0
      local dy = x["y"] - 3707270.0
      local dz = x["z"] + 4965.62
      local r2 = dx*dx + dy*dy + dz*dz
      local R  = 8000.0
      return { Tnuc_s = 30.0e6 * math.exp(-r2 / (R*R)) }
    end

[Tnuc_n, Tnuc_d]: !ConstantMap
  map:
    Tnuc_n: 0.0
    Tnuc_d: 0.0
```

### Edge Cases to Handle
- **Spatial `rs_b` support.** First run log must show `RS parameter source ... b 0` (b from easi).
  If the build cannot read spatial `rs_b` for FL=3/4, b silently falls back to the constant
  `RS_b`; document and set `RS_b` to the hypocenter value 0.0168 (acceptable for the activity, not
  for production).
- **No `rs_srW`/`rs_muw` anywhere** for FL=3/4 (the law does not use them; `insertIfPresent`
  would read them but they would be inert — keep them absent to avoid confusion).

### Acceptance Criteria
- [ ] `safs_fault_rs.yaml` parses; SeisSol reads `rs_a`, `rs_b`, `rs_sl0`, `Tnuc_s` without error.
- [ ] `safs_material.yaml` present with the three constant values.
- [ ] Run log confirms RS (not LSW) and the b-source line.

### Dependencies
- Depends on: nothing. Required by: Phase 2.

## Phase 2: Parameter decks (aging FL=3 / slip FL=4) + combined output + size budget

### Goal
Four ready-to-submit decks (2 laws x 2 materials) exist, each emitting fault + free-surface +
energy output, each provably < 1 GB to download.

### Files to Create
- `parameters_aging.par`  — FL=3, defaults to CONSTANT material (Group 1, case 1).
- `parameters_slip.par`   — FL=4, defaults to CONSTANT material (Group 2, case 1).
- (CVM cases are the same two decks with one line changed — see requirement 3.)

### Detailed Requirements
1. Create `parameters_aging.par` with EXACTLY this content:

```fortran
! SAFS GROUP ACTIVITY — Group 1: rate-and-state AGING law (FL=3), NO strong rate weakening.
! Slip law (Group 2): copy to parameters_slip.par and set FL = 4 (only that line changes).
! Material: CONSTANT first (M0 = 3.2e10*P0 exactly verifiable); for the CVM case set
!   MaterialFileName = 'safs_material_cvm.yaml'.
&equations
MaterialFileName = 'safs_material.yaml'        ! CONSTANT. CVM case: 'safs_material_cvm.yaml'
Plasticity = 0
numflux = 'godunov'
numfluxnearfault = 'godunov'
/
&IniCondition
/
&DynamicRupture
FL = 3                                          ! 3 = aging (Group 1); 4 = slip (Group 2)
ModelFileName = 'safs_fault_rs.yaml'
! (GP-wise fault-parameter evaluation is the DEFAULT in SeisSol >= v1.3; no flag needed.
!  If the target QuakeWorx SeisSol is OLDER and requires it, add: GPwise = 1)
RS_f0  = 0.6
RS_sr0 = 1d-6
RS_b   = 0.0168                                 ! FALLBACK; spatial rs_b(z) from yaml (verify log "b 0")
RS_iniSlipRate1 = 1d-12                          ! V_init magnitude sets initial state
RS_iniSlipRate2 = 0.0
s_0 = 0.0
t_0 = 1.0                                        ! nucleation ramp duration [s]
XRef = 0.0
YRef = -1.0
ZRef = 0.0
refPointMethod = 1
OutputPointType = 4
/
&Elementwise
printtimeinterval_sec = 2.0                      ! fault snapshot cadence (STF + Vr/PSR maps)
! 1=SRs/SRd 7=Vr 8=ASl 9=PSR 10=RT selected:
OutputMask = 1 0 0 0 0 0 1 1 1 1 0 0
refinement_strategy = 2
refinement = 0                                   ! 1 value/facet. KEEP 0 for <1 GB.
/
&Pickpoint
/
&SourceType
/
&SpongeLayer
/
&MeshNml
MeshFile = 'safs_mesh.puml.h5'
meshgenerator = 'PUML'
/
&Discretization
CFL = 0.5
ClusteredLTS = 2
FixTimeStep = 0.1
/
&Output
OutputFile = 'output/safs'
Format = 10                                      ! no volume wavefield output
iOutputMask = 0 0 0 0 0 0 1 1 1
iPlasticityMask = 0 0 0 0 0 0 1
TimeInterval = 5.0                              ! moot (Format=10)
refinement = 1                                  ! VOLUME refinement; moot under Format=10. NOT the
                                                ! fault &Elementwise refinement (that is 0, above)
! --- Free surface output (task 3 PGV): LEGACY COMPACT path (geometry once) ---
SurfaceOutput = 1
SurfaceOutputRefinement = 0
SurfaceOutputInterval = 1.0                      ! v1/v2/v3 -> PGV via ParaView temporal max
surfacevtkorder = -1                             ! MANDATORY: -1 = compact. 2 would be multi-GB.
! --- Energy output (tasks 2: M0, P0 ground truth) ---
EnergyOutput = 1
EnergyTerminalOutput = 1
EnergyOutputInterval = 0.25
Checkpoint = 0
/
&AbortCriteria
EndTime = 100.0
/
&Analysis
/
&Debugging
/
```

2. Create `parameters_slip.par` = byte-identical to `parameters_aging.par` EXCEPT `FL = 4` and the
   header comment.
3. CVM cases: duplicate each deck (or change at submission) with
   `MaterialFileName = 'safs_material_cvm.yaml'`. Keep `safs_material_cvm.nc` in the run dir; the
   QuakeWorx SeisSol app must have ASAGI (it does — Kaikoura training case uses it).

### Download-size budget (per run, EndTime=100 s)

```
Fault &Elementwise (legacy XDMF, geometry once, refinement=0, 6 components):
  Nf=121,316; per snapshot = 6*121316*REAL ; snapshots = 100/2 + 1 = 51
    double: 6*121316*8 = 5.82 MB * 51 = 297 MB ; single: 149 MB ; + ~10 MB geometry
Free surface (legacy compact, refinement=0, 6 fields v1..u3, geometry once):
  Ns=34,096; per snapshot = 6*34096*REAL ; snapshots = 100/1 + 1 = 101
    double: 6*34096*8 = 1.64 MB * 101 = 165 MB ; single: 83 MB ; + ~1 MB geometry
Energy CSV (0.25 s, ~400 times * ~13 variables): < 1 MB
---------------------------------------------------------------------------
TOTAL per run  ~=  0.47 GB (double)   /   ~=  0.24 GB (single)     << 1 GB
```

**Answer: each run clears 1 GB with margin.** Levers (document in the worksheet):
- PGV looks clipped -> `SurfaceOutputInterval = 0.5` (surface ~330 MB; total ~0.64 GB double, still < 1 GB).
- Tighter quota -> `printtimeinterval_sec = 5.0` (fault ~122 MB; total ~0.29 GB double).
- NEVER set `surfacevtkorder >= 0` or fault `refinement = 1` on a double build (each blows past 1 GB).

### Acceptance Criteria
- [ ] All four decks parse; SeisSol selects FL=3 / FL=4 (log) and reads the constant or CVM material.
- [ ] `du -sh output/` for a finished run is 0.25-0.65 GB depending on precision/interval.
- [ ] Output contains `*-fault.xdmf` (arrays SRs,SRd,Vr,ASl,PSR,RT), `*-surface.xdmf`
      (v1,v2,v3,u1,u2,u3), and `*-energy.csv` (variables seismic_moment, potency).

### Dependencies
- Depends on: Phase 1. Required by: Phases 3-6.

## Phase 3: Constant-material ParaView verification session (tasks [1] learn + [2])

### Goal
A `WORKSHEET.md` that takes a participant from the downloaded constant-material run to five
quantities computed in ParaView, each checked against SeisSol ground truth — establishing trust in
the pipeline before the CVM case.

### Files to Create
- `WORKSHEET.md`

### ParaView feasibility (state up front)
- Scale: 121,316 fault triangles + 34,096 surface triangles x ~51-101 steps — small; runs on a
  laptop (a few hundred MB RAM). Open `output/safs-fault.xdmf` and `output/safs-surface.xdmf`
  (SeisSol XDMF, native ParaView, time series).
- Stock filters only: `Calculator`, `Cell Data to Point Data`, `Gradient`, `Integrate Variables`,
  `Plot Data Over Time`, `Temporal Statistics`, `Threshold`.
- Gotcha: fault/surface legacy output is CELL data. `Integrate Variables` on cell data does
  `sum(value*area)` (what we want). `Gradient` needs POINT data, so convert RT with
  `Cell Data to Point Data` first.

### The five verification tasks (identical pipeline for aging and slip; only the data differs)
1. **Potency P0** — last step; `Integrate Variables` on fault -> integrated `ASl` = integral(slip dA)
   [m^3]. Truth: energy CSV `variable==potency`, last time. Expect < ~1% (quadrature).
2. **Seismic moment M0** — `M0 = 3.2e10 * P0`; `Mw = (2/3)*log10(M0) - 6.07`. Truth: CSV
   `seismic_moment`, last time. Expect P0-band agreement (constant mu = harmonic mean exactly).
3. **Moment-rate / STF (CONSTANT material only for the `x 3.2e10` form)** — `Integrate Variables`
   -> `Plot Data Over Time` of `ASl` x 3.2e10 = M0(t); finite-difference for the rate. Truth: CSV
   `seismic_moment(t)` differenced at 0.25 s. Expect matching shape; ParaView coarser (2 s).
   ON CVM: the `x 3.2e10` factor is INVALID (mu varies); use SeisSol `seismic_moment(t)` directly,
   or the Phase-4 resampled-mu integral. Reusing `x 3.2e10` on CVM gives a wrong STF.
4. **Rupture velocity Vr** — `Threshold RT>0` -> `Cell Data to Point Data` -> `Gradient` of RT ->
   `Calculator 1/mag(RT_gradient)`. Truth: fault array `Vr`. Expect interior agreement; edge/kink
   outliers. Classify supershear where `Vr > Vs = 3462 m/s` (uniform here — constant medium).
5. **Peak slip rate PSR** — `Calculator sqrt(SRs^2+SRd^2)` -> `Temporal Statistics (Maximum)`.
   Truth: fault array `PSR`. Expect ParaView <= SeisSol (2 s undersampling) — the Nyquist lesson.

### Acceptance Criteria
- [ ] A ParaView novice reaches all five comparisons with only the listed filters.
- [ ] Each task names the ParaView source AND the SeisSol ground truth (array or CSV variable).
- [ ] Units stated everywhere (P0 m^3, M0 N*m, Vr m/s, PSR m/s).

### Dependencies
- Depends on: Phases 1-2 + one constant-material run per group. Required by: Phases 4-5.

## Phase 4: CVM-material run + material-effect comparison (task [1])

### Goal
Participants repeat the verified pipeline on the CVM run and quantify how heterogeneous material
changes the source quantities and rupture.

### Detailed Requirements (add a CVM section to `WORKSHEET.md`)
1. Re-run the same group's deck with `MaterialFileName = 'safs_material_cvm.yaml'`.
2. Recompute P0 (still exact: `Integrate ASl`), Vr, PSR, PGV exactly as in Phase 3.
3. **M0 caveat (teaching point):** with the CVM, `M0 = integral(mu(x)*slip dA)` and mu varies on
   the fault, so `M0 != mu_const * P0`. ParaView cannot form the exact M0 without mu on the fault
   (SeisSol does not output mu). Options to document:
   - (a) Compare ONLY P0 and the SeisSol-reported `seismic_moment`; explain the inequality.
   - (b) ADVANCED: load `safs_material_cvm.nc` in ParaView, `Resample With Dataset` mu onto the
     fault, `Calculator mu*ASl`, `Integrate Variables` -> approximate M0; compare to CSV.
   - The **STF/moment-rate is mu-weighted the same way**: do NOT reuse Phase-3's `x 3.2e10` on
     CVM; read the moment-rate from SeisSol `seismic_moment(t)` directly (or via option (b)).
4. Compare const vs CVM: P0, M0, peak Vr, supershear fraction, PGV pattern. Tabulate the deltas.

### Acceptance Criteria
- [ ] CVM P0 (ParaView) matches CSV potency < 1%.
- [ ] The const-vs-CVM comparison table is filled with both groups' numbers.
- [ ] The M0 inequality on CVM is explained (mu-weighting), not reported as an error.

### Dependencies
- Depends on: Phase 3 + one CVM run per group.

## Phase 5: PGV — rupture-velocity correlation (task [3])

### Goal
Each group produces a surface PGV map and correlates it with the fault rupture velocity, looking
for supershear -> PGV-high relationships.

### Detailed Requirements (add a task-3 section to `WORKSHEET.md`)
1. **PGV map** from the surface output: open `safs-surface.xdmf`; `Calculator vmag = sqrt(v1^2+v2^2+v3^2)`;
   `Temporal Statistics (Maximum)` -> `vmag_maximum` = PGV (m/s). (No SeisSol PGV ground truth —
   this is the analysis product; trust comes from the Phase-3 verification of the same temporal-max
   operation used for PSR.)
2. **Rupture velocity map** = the verified fault `Vr` (or the participant's `1/|grad RT|`).
3. **Overlay / correlate** (surface map vs fault-at-depth):
   - Visual: render the z=0 PGV surface and the fault colored by Vr in one view (map/top view).
   - Quantitative: extract along-strike profiles — `Plot Over Line` (or `Resample With Dataset`
     onto a polyline following the surface fault trace) for PGV, and the fault top-edge Vr — and
     overlay vs along-strike distance; compute the correlation coefficient.
   - Hypothesis to test: supershear segments (`Vr > 3462 m/s`) sit beneath PGV highs (Mach-cone
     enhancement). Contrast aging vs slip groups, and const vs CVM.

### Acceptance Criteria
- [ ] PGV map produced for both materials and both groups.
- [ ] An along-strike PGV-vs-Vr overlay plot exists with a stated correlation.
- [ ] A one-paragraph interpretation (does faster rupture => stronger surface shaking here?).

### Dependencies
- Depends on: Phase 3 (PGV uses the verified temporal-max) + the surface output from Phase 2.

## Phase 6: Facilitator answer key + mandatory smoke test

### Goal
A `ANSWER_KEY.md` with reference numbers, tolerances, and discrepancy explanations, AND a
documented pre-session smoke test that confirms each FL nucleates and produces a measurable event.

### Detailed Requirements
1. **Smoke test BEFORE distributing** (this is also the top-risk mitigation): run each
   `parameters_{aging,slip}.par` on the CONSTANT material to t~=20-30 s and confirm the energy CSV
   `seismic_moment` grows (rupture nucleated and propagated past the patch). If it stalls at the
   nucleation patch, apply mitigations from Risk Assessment and re-run.
2. Record per (law x material): final P0, seismic_moment, Mw; ParaView integrated-ASl vs CSV
   potency %; STF peak + time; Vr interior median + supershear fraction; PSR (ParaView vs SeisSol);
   PGV peak + the PGV-Vr correlation.
3. Tolerance table: P0 < 1%; M0 = P0 band (const); STF peak < 10%, peak time within 1 snapshot;
   Vr interior median < 10% (edge outliers ungraded); PSR ParaView <= SeisSol (teaching).
4. Explain each expected mismatch (quadrature; polynomial vs linear gradient; temporal sampling;
   ASl path-integral vs net slip; CVM mu-weighting of M0).

### Acceptance Criteria
- [ ] Smoke test shows nucleation + propagation for FL=3 AND FL=4 on the constant material.
- [ ] Every worksheet quantity has a reference value and tolerance.
- [ ] Each discrepancy has a cause, not just a number.

### Dependencies
- Depends on: Phases 1-5 + four official runs.

## Testing Strategy
- **Deck validity:** all four decks parse; logs show FL=3/FL=4, RS b-source, constant/CVM material.
- **Budget:** `du -sh output/` per run is 0.25-0.65 GB.
- **Self-consistency (core test):** P0(ParaView) ~= P0(CSV) and (constant) M0 = 3.2e10*P0 ~=
  seismic_moment(CSV) < 1%.
- **Nucleation/propagation:** energy CSV moment grows for both laws and both materials (smoke test).
- **Reproducibility:** a second participant on the same files gets the same numbers (deterministic).

## Risk Assessment
- **(Highest) Plain RS may under-propagate / fail to run away from the nucleation patch.** Without
  strong rate weakening, steady-state weakening is only `(a-b)*ln(V/V0)` (~5% of f0 from V0=1e-6 to
  ~1 m/s), and the RS nucleation length `L_nuc = mu*Dc/((b-a)*sn) ~= 9 km (CVM) / ~12.5 km (const)`
  is large vs the forced overstress footprint (~1.5 km). The forced patch WILL slip (so M0/P0/Vr/PGV
  are always non-zero and the verification works regardless), but it may not propagate far.
  Detect: `seismic_moment` plateaus right after `t_0`. Mitigations, in order: (a) enlarge the
  nucleation patch `R` 6000 -> 8000-9000 m so the forced region ~ L_nuc; (b) raise `Tnuc_s`
  24 -> 30 MPa; (c) widen/strengthen the VW band (make `amb` more negative or extend to d<=12 km);
  (d) accept a confined event for the verification activity. DECIDE via the Phase-6 smoke test.
- **Constant-material nucleation margin.** mu 23.4 -> 32 GPa raises L_nuc ~1.37x, so the constant
  case is the harder one to nucleate — run the smoke test on CONSTANT first.
- **Aging vs slip differ in rupture, by design.** Slip law weakens/localizes faster (often a
  slightly larger, sharper event) than aging (which heals). This is the inter-group comparison, not
  a bug — but it means the two groups' reference numbers differ; the ANSWER_KEY needs both.
- **Spatial `rs_b` not read on older builds** -> b falls back to constant `RS_b` (wrong away from
  hypocenter). Verify the log `b 0`; if `b 1`, document the limitation (still fine for the activity).
- **Surface output size trap.** `surfacevtkorder >= 0` or fault `refinement = 1` each break the
  1 GB budget on a double build — both are pinned in the decks; do not change.
- **ParaView gradient noise** at the curved/segmented SAF surface and fault edges — mask `RT>0` and
  compare interior medians, not per-cell extremes.
- **CVM M0 not exactly verifiable** in ParaView (mu varies, not output) — by design the constant
  case is the exact-verification anchor; CVM M0 is compared via the SeisSol value + the optional
  resample-mu extension.

## Deliverables checklist (in safs_seisol_v2_2_0_GroupActivities/)
- [ ] `safs_fault_rs.yaml`            (Phase 1)
- [ ] `safs_material.yaml`            (Phase 1, copied from archive/)
- [ ] `parameters_aging.par` (FL=3)   (Phase 2)
- [ ] `parameters_slip.par`  (FL=4)   (Phase 2)
- [ ] `WORKSHEET.md`                  (Phases 3-5)
- [ ] `ANSWER_KEY.md`                 (Phase 6, after 4 official runs + smoke test)
- [ ] `GROUP_ACTIVITY_PLAN.md` (this file) + `GROUP_ACTIVITY_PLAN.pdf`

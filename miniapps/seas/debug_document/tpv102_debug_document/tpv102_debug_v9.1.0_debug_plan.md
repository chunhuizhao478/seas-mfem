# TPV102 Debug v9.1.0 Plan

> **Author:** debugger agent. **Date:** 2026-04-20 (initial);
> **rev. 2026-04-20 post-review** to incorporate `REVIEW.md` findings
> R-001 – R-004. See §11 for the review-incorporation index.
> **Successor to** `tpv102_debug_v9.0.0_debug_plan.md` (Pelties-9 per-side
> flux patch) and `tpv102_debug_v9.0.0_seissol_flux_comparison.md`.
> **Artifacts this plan consumes:**
> - Station plots: `tpv102/plots_results_200m_p1_2.0s_400r_vmax_halt_job7667881/`
> - ParaView slip-rate screenshot: `~/Desktop/Screenshot 2026-04-20 at 7.55.14 PM.png`
> - Post-fix source (v9.0.0 Pelties-9 applied): `dynamic/wave_operator.inl`
>   lines 832–918 (interior branch) and 1240–1300 (shared branch).
> - Review of this plan + associated code:
>   `/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/REVIEW.md`.
>
> **Status:** investigation only. No source edits applied. No Frontera jobs
> launched. Waiting for user approval per `feedback_frontera_approval.md`
> before any `sbatch`.
>
> **Biggest revision after review.** R-V91-A's rank-1 hypothesis has
> changed: H-V91-A4 (QP values written at reference-triangle vertex
> positions) replaces H-V91-A1 (DG-jump artefact). The old A1/A2/A3 set
> is retained as secondary hypotheses but the reviewer's identification
> of the output-writer bug in `paraview_output.hpp:566–596` is now the
> primary lead. The Frontera 12 s re-run (§6) remains required but only
> for R-V91-B; R-V91-A is closable by the R-001 source patch alone.

---

## TODO — remaining actions to close R-V91

### Blocking path (must all pass)

- [ ] **YOU:** Read §2 (now updated with H-V91-A4), §3, and §11.
      Authorise the R-001 / R-002 / R-003 / R-004 source patches listed
      in §11 (they live in `paraview_output.hpp` and `fault_basis.hpp`).
      Per `CLAUDE.md` "Proposing Fixes", no source edit happens without
      this.
- [ ] **Debugger:** Apply R-001 / R-002 / R-003 / R-004 patches in a
      single C0 commit (see updated §8). R-001 subsumes R-002's
      loop-bound hardcoding.
- [ ] **Debugger:** Add the two diagnostic unit tests in §4 — now
      written against the post-R-001 CellData contract, not the
      pre-R-001 PointData contract — in a C1 commit.
- [ ] **Debugger:** Run `make test`; expect all new tests green and
      v9.0.0 regression suite still green.
- [ ] **Debugger:** Run the offline diagnostics in §5 against the existing
      job-7667881 VTU / station output. Post-R-001, §5.1 / §5.2 / §5.3
      become *supplementary* (verifying the fix closed R-V91-A); §5.4 /
      §5.5 remain primary for R-V91-B.
- [ ] **YOU:** Authorise the Phase-9.1 Frontera 12 s re-run (§6).
      Parameters drafted, sbatch not yet generated per
      `feedback_frontera_approval.md`. R-001 must be in-tree on the
      build before this sbatch is submitted so the 12 s fault-surface
      output is speckle-free from the start.
- [ ] **Debugger:** Generate
      `tpv102_200m_p1_12.0s_400rank_v91_run.sbatch` from the working
      module template once authorised; submit via `ibrun`.
- [ ] **Frontera:** Submit Phase 9.1; expect (a) no isolated-element speckles
      in the slip-rate VTU (R-V91-A sanity check post-R-001),
      (b) dip-channel deviation at `flt_0_7.5` within the SeisSol/PyLith
      envelope (< 2× reference noise floor) past t ≈ 5 s (R-V91-B pass).
- [ ] **Debugger:** Write `tpv102_debug_v9.1.0_fix.md` after each step,
      one entry per substantive action.
- [ ] **Debugger:** When R-V91 closes, move this plan's TODO checklist
      to DONE and update R-items in a new `tpv102_debug_v9.1.0_check.md`.

### Frozen / deferred

- Any attempt to "smooth" the fault-surface VTU output by welding adjacent
  triangles or switching to a continuous L2 projection. **Superseded by
  R-001** (cell-data rewrite in §11). Do not pursue a welding / L2
  projection variant — it would over-smooth compared to the face-averaged
  cell-data form prescribed by R-001.
- Any change to `FaultFaceFlux::Evaluate` eq.(7)–(12) or
  `GodunovFlux::Interior(n, Q, Q)` per-side call. The v9.0.0 Pelties-9
  patch has not been regressed — the two dip-deviation hypotheses in §3
  do not touch the imposed-state construction, and the reviewer's
  compliance check in `REVIEW.md` confirms the per-side flux in
  `wave_operator.inl:832-918` + `:1326-1340` agrees with
  `fault_face_flux.cpp`'s (7)/(9)/(11)–(12) pipeline.

---

## Dashboard — read me first

**Where we are.** v9.0.0 closed H-V9-J: the DG flux at fault faces now
evaluates `A·T·Q̃` per-side (Pelties 2012 eq. 9) instead of feeding both
imposed states into a welded-face Riemann solver. Frontera job 7667881
(200 m, P1, 400-rank, tfinal = 2.0 s) is the first post-fix production
run. Radiation is clearly now entering the bulk — V_strike peaks decay
(vs the v8.0.0 pin at 7.69871 m/s flat plateau), slip accumulates
off-axis, and off-hypocenter dip channels show the expected magnitudes
from elastic coupling. Two new issues surface from the post-fix plots:

- **R-V91-A — ParaView slip-rate surface shows isolated "speckle"
  elements** with visibly different values from their immediate
  neighbours, scattered across the rupture front (see screenshot).
  Unclear whether this is (a) a DG-visualisation artefact from the
  disconnected-triangle VTU writer, (b) a partition-boundary effect
  at the 400-rank decomposition, or (c) a DOFData mapping bug
  in the output writer.
- **R-V91-B — Hypocenter station flt_0_7.5 shows larger dip transients
  than SeisSol/PyLith reference.** Specifically at the hypocenter
  (x=0, down-dip=7.5 km), SEAS-MFEM's `slip_rate_dip` spikes to
  ~0.02 m/s and `τ_dip` spikes to −0.3 MPa at breakaway (t ≈ 1.5 s),
  whereas reference codes show ±0.003 m/s and ±0.05 MPa respectively.
  Off-axis stations (flt_9_7.5, flt_n9_7.5) agree with reference to
  within ~10%, so the effect is localised to the hypocenter / early
  rupture.

Running to tfinal = 12 s is required to see whether R-V91-B is a
transient-only deviation that decays (i.e., the hypocenter is slightly
more reactive to breakaway noise in MFEM but converges to reference) or
a persistent systematic offset (i.e., a coupling asymmetry between dip
and strike channels). At tfinal = 2 s we cannot distinguish these two.

**What you decide today.** Pick the first-cut hypothesis in §2 (R-V91-A)
and §3 (R-V91-B). No source is edited regardless of the choice — the §5
offline diagnostics and §4 unit tests can run locally.

### Phase progress

| Phase | Scope | Status | Detail |
|---|---|---|---|
| 0 | R-001 / R-002 / R-003 / R-004 source patches (commit C0) | **PROPOSED — awaiting user approval** | §11, §8 |
| 1 | Unit-test additions (output-consistency post-R-001, basis-symmetry post-R-003) | **PROPOSED** — §4 | run via `make test` on laptop |
| 2 | Laptop offline diagnostics against job 7667881 VTU + station output | **PROPOSED** — §5 | no build required; post-R-001 these become verification rather than diagnostic |
| 3 | Frontera 12 s re-run at same 200 m / 400 rank | **PROPOSED — awaiting user approval** | §6 |
| 4 | Cross-check against SeisSol 12 s reference dataset (paper dataset or new SeisSol run) | **OPEN** — §3.C | depends on data availability |

### Hypothesis status at a glance

| H-V91- | Status | Evidence | Test / diagnostic |
|---|---|---|---|
| **A4** | **RANK 1 for R-V91-A** (NEW, added post-review) | VTU writer at `paraview_output.hpp:566-596` places three per-QP values at reference-triangle **vertex** positions while the QPs themselves live at **interior** barycentric points `(1/6,1/6), (2/3,1/6), (1/6,2/3)` (MFEM `intrules.cpp:1271-1277` degree-2 rule, `// interior points`). Per-vertex rendering thus misinterprets interior samples as corner values. | R-001 source patch (see §11) — falsifies by fixing; §4.1 post-R-001 CellData continuity test |
| A1 | Demoted to RANK 2 | Any residual element-scale colour variation *after* R-001 is applied would be true DG-jump behaviour | §4.1 test still catches this as a secondary check |
| A2 | OPEN for R-V91-A (low priority after R-001) | 400-rank partition may interact with shared-fault DOF mapping | §5.2 per-rank speckle-map audit |
| A3 | OPEN for R-V91-A (low priority after R-001) | Output writer may read stale DOFData on orphan / shared-as-bdr faces | §5.3 fault-DOF-zero audit |
| **B1** | **RANK 1 for R-V91-B (until 12s run)** | FaultBasis dip/strike computed asymmetrically (`strike = up×n` with explicit normalize, `dip = strike×n` without explicit normalize). Review R-003 quantifies the drift at ~5 ULP — too small to produce 0.3 MPa spike alone, but orthonormality-enforcing patch is cheap and defensive | R-003 source patch; §4.2 unit test + §5.4 station-hierarchy audit |
| B2 | OPEN for R-V91-B | Godunov frame rotation uses `can_t1 = dip` first — ordering may perturb T_can numerical conditioning | §5.5 permutation diagnostic |
| B3 | OPEN for R-V91-B | SEAS-MFEM friction solver solves for V_abs only, then decomposes to V1/V2 scalar-multiplicatively; SeisSol uses Kaneko 2008 vector iteration | §3.4 SeisSol code audit (mostly DONE in §3.3 — algebraic equivalence confirmed); 12 s Frontera re-run remains the decisive test |

---

## 1. Symptoms — what the v9.0.0 post-fix run shows

**Run:** Frontera job 7667881 (sbatch submitted 2026-04-19 approx),
mesh = 200 m uniform, P1, 400 MPI ranks on 8 nodes, tfinal = 2.0 s,
v9.0.0 Pelties-9 per-side flux applied.
**Output directory on laptop:**
`tpv102/plots_results_200m_p1_2.0s_400r_vmax_halt_job7667881/`

### 1.1 Overall slip-rate signature (`tpv102_overview.png`)

V_strike peaks arrive at each station in order of rupture-front
arrival time: hypo (flt_0_7.5, yellow) at t ≈ 1.5 s with peak 4.1 m/s;
off-strike stations fire in sequence and the farthest (flt_12_12 green)
peaks at t ≈ 6.5 s with the largest peak 8.1 m/s. All peaks decay
monotonically after their arrival. **This is a working rupture.** The
R-1001 pin (v8.0.0) would have shown everything at 7.69871 m/s flat. So
the v9.0.0 Pelties-9 fix is producing radiation.

### 1.2 Hypocenter station `flt_0_7.5`

| Channel              | MFEM (black)      | DR/PyL reference (red/blue dashed) |
|----------------------|-------------------|-----------------------------------|
| V_strike peak        | 4.1 m/s           | 4.1 m/s (agreement to ~3%)        |
| V_strike decay       | matches reference | matches                           |
| Slip_strike final (t=2s) | 6 m            | 6 m (match)                       |
| Shear_strike         | 75 → 65 MPa       | match                             |
| Normal stress        | 120 MPa (flat)    | 120 MPa (flat)                    |
| log10(state)         | matches           | matches                           |
| **V_dip**            | spike ±0.02 m/s at t≈1.5 s, decays | noise ±0.003 m/s |
| **Slip_dip**         | transient −5e-4 m, recovers to 0 | noise ±1e-4 m |
| **Shear_dip**        | transient −0.3 MPa, recovers     | noise ±0.1 MPa |

The dip channels are **non-zero** (unlike v8.0.0's exactly-zero pin —
which was the H-V9-J smoking gun). They have the **wrong magnitude and
shape** compared to reference: MFEM shows a coherent breakaway spike
whereas reference shows continuous noise. At tfinal = 2 s the spike has
not yet fully decayed. **Open question:** does it decay to reference
envelope by t = 12 s, or does it persist / grow?

### 1.3 Off-hypocenter stations

| Station         | V_strike | Slip_strike | Slip_dip | Comment |
|-----------------|----------|-------------|----------|---------|
| `flt_0_3`       | peak 3.5 m/s at t≈3 s | 8 m | ≈ 0 | matches reference |
| `flt_0_12`      | peak 3.4 m/s at t≈3.5 s | 6 m | +8e-4 (MFEM) vs same ref | match |
| `flt_9_7.5`     | peak 5.2 m/s at t≈4.8 s | 6 m | −0.012 m | matches reference |
| `flt_n9_7.5`    | peak 5.5 m/s at t≈5 s | 7 m | +0.022 m | matches reference |

**Off-axis stations show reference-level dip-slip magnitudes** (the
non-zero "dip slip" off-axis is physical — it's geometry, rupture-front
radiation in 3D elastic medium). Only the on-axis hypocenter station
shows an unphysical spike. The spike is localised in both space and
time.

### 1.4 Slip-rate surface (ParaView screenshot)

Coloured rupture disc at t ≈ 1.5 s, centered at the hypocenter. Front
is a clear ring of high slip rate (orange/red). Inside the ring,
interior mostly blue (≈ 0) as expected (no back-slip). **But** scattered
within and around the ring are small isolated elements with visibly
different colours from their neighbours — some white/yellow hot spots
inside the ring interior; some small dark patches along the ring.

These are **element-scale speckles**, not a smooth wave pattern. They
do not correspond to any physical structure.

---

## 2. Issue R-V91-A — slip-rate surface speckles

### 2.1 Code path that writes the fault surface

`io/paraview_output.hpp::WriteFaultSurfaceVTU` (lines 517–722) writes
per-rank `fault_surface_r{rank}_c{cycle}.vtu` as an ASCII triangle mesh
with per-vertex data arrays. Key structural details:

- **Triangles are disconnected** (lines 567–577): every triangle has
  its own three vertices; adjacent triangles' vertices are duplicated
  in the vertex list and never welded. So within each triangle ParaView
  linearly interpolates the three per-vertex values, but **across adjacent
  triangles the interpolation is piecewise-independent**.
- **"Vertex" values are actually QP values at interior positions**
  — this is the post-review root cause (R-001 / H-V91-A4). Lines
  579–596 read `local_slip_rate(2*d)` etc. at
  `d = face_mesh_idx * nbf + k`. The DOF layout is QP-ordered
  (stride `nqp_per_face`, enforced by `WaveOperator::SetFaultDOFData`
  in `wave_operator.hpp:151`). For P1 triangle faces, the quadrature
  rule is `IntRules.Get(TRIANGLE, 2*order) = degree-2 rule`, which
  per MFEM `fem/intrules.cpp:1271-1277` is `AddTriPoints3(0, 1./6., 1./6.)`
  at the **interior** barycentric points
  `(1/6, 1/6), (2/3, 1/6), (1/6, 2/3)`. The VTU writer places these
  three interior-QP values at the three **corner** positions
  `(0, 0), (1, 0), (0, 1)` of the reference triangle (lines 566–577).
  ParaView then linearly interpolates from the "corners"; the
  resulting rendering samples the field at corner positions it never
  actually had, and — because adjacent triangles' interior QPs sit
  on opposite sides of the shared physical edge — the reconstruction
  can disagree across that edge by `O(h · |∇f|)` even when the
  underlying DG field is continuous. This produces the element-scale
  speckle in the screenshot.
- **R-002 latent bug:** the per-"vertex" loop is hardcoded
  `for (int k = 0; k < 3; k++)` (line 581), matching `nbf = 3` for
  P1 triangles by coincidence. For P2 triangles `nbf = 6` and QPs
  3-5 would be silently dropped. Fix subsumed by R-001 cell-data
  rewrite (the averaging loop iterates `k < nbf`).
- **Only `FTr != nullptr` faces are written** (line 564 early return):
  orphan faces (fi with `fault_face_elem1_[fi] == -1` in
  `InitFaultOutputBP5`) exist in the DOF vector but are skipped by the
  VTU writer.
- **Per rank writes its own VTU; rank 0 writes `.pvtu` index**
  (lines 690–721) listing every rank's `fault_surface_r{r}_c{c}.vtu`.
  ParaView loads the `.pvd` → all `.pvtu` → all `.vtu`; no further
  blending across ranks.

### 2.2 Why the L2-p0 volume-projection was deprecated

The docstring at lines 505–510 is explicit:

> "This replaces the L2-p0 volume projection which creates
> scattered-dot artifacts in ParaView."

So a prior version had scattered-dot artefacts from L2-p0 in the bulk
mesh. That version was replaced by this VTU-triangle writer. The
current speckles are **not** the same artefact as before — they are on
the fault surface mesh directly.

### 2.3 Hypotheses

#### H-V91-A4 — QP/vertex-position mismatch (RANK 1 post-review)

**Claim.** The speckles are caused by `WriteFaultSurfaceVTU`
(`paraview_output.hpp:566-596`) reading per-QP DOF values and writing
them as per-vertex PointData on triangles whose reference-frame
vertex positions are `(0,0), (1,0), (0,1)`, while the three QPs those
values came from live at `(1/6, 1/6), (2/3, 1/6), (1/6, 2/3)`.
ParaView's linear interpolation across the reference triangle thus
reconstructs the field at positions it was never sampled at, with
an error of order `h · |∇f|` that crosses element boundaries because
adjacent triangles' interior QPs are on opposite sides of the shared
edge. This is plan-level RANK 1 post-review.

**Evidence for.** Direct code reading: MFEM `fem/intrules.cpp:1271-1277`
degree-2 triangle rule uses `AddTriPoints3(0, 1./6., 1./6.)` at
interior barycentric coordinates (explicit `// interior points` comment).
`WaveOperator::SetFaultDOFData` (`wave_operator.hpp:151`) lays out
DOFs with stride `nqp_per_face`. `WriteFaultSurfaceVTU` pairs these
QP values with reference-frame corner positions without a
QP→vertex reconstruction. Reviewer's own analysis in `REVIEW.md` R-001
§Description.

**Evidence against.** None — the bug is a direct read of the code.

**Diagnostic to decide.** R-001 source patch: replace PointData with
CellData carrying the per-face QP-average. Applies at
`paraview_output.hpp:540-596` (accumulation loop) and `:690-716`
(PVTU block). Post-patch the speckles must disappear on existing
job-7667881 data when re-rendered. §4.1 unit test
`test_fault_surface_vtu_continuity` uses a known-smooth linear field
`f(x2,x3) = 0.5 + 0.1·x2 + 0.2·x3` and asserts (a) each cell value
equals the face-centroid average to 1e-10 relative, (b) adjacent-cell
value jump is ≤ 1e-12 for this linear field.

**If confirmed.** R-001 lands in commit C0 (§8). No 12 s Frontera
re-run is required to close R-V91-A — the fix is a pure output-path
change and can be validated on existing job-7667881 raw VTU by
re-writing with the patched writer. Speckles should vanish from
subsequent fault-surface plots.

#### H-V91-A1 — DG element-jump artefact (demoted to RANK 2)

**Claim.** After R-001 lands, any remaining element-scale colour
variation in the rupture-front region is genuine DG element-jump
behaviour from the piecewise-constant CellData form. This would not
be "speckle" any more — it would be a uniform ring of discontinuous
cell colours following the front.

**Evidence for (post-R-001).** P1 face-local DG basis admits genuine
jumps at element boundaries; a face-averaged CellData value smooths
within the face but preserves the inter-face jumps.

**Evidence against.** The speckles in the screenshot appear
*inside* the rupture disc interior where the field is smooth, not
along the front where DG jumps would concentrate — consistent with
H-V91-A4 (interior QPs rendered at vertex positions introduce error
proportional to the field gradient, which is nonzero even inside the
rupture zone from wave radiation), not with H-V91-A1.

**Diagnostic to decide.** §4.1 unit test above is also the A1
falsification test — CellData under a linear field must produce
zero inter-cell jumps up to 1e-12 by construction, so if the test
fails that points at a residual bug beyond A4.

**If residual issue confirmed.** Switch to a P1 nodal reconstruction
per-triangle (reviewer's R-001 Option 2: solve 3×3 linear map from
QPs to vertices). Larger surface change than R-001 Option 1;
deferred.

#### H-V91-A2 — Partition-boundary DOF mismatch

**Claim.** At 400-rank decomposition there are many shared-fault faces
between rank pairs. Per `wave_operator.inl:396–446`, each rank computes
`elem1_on_plus_` for shared faces based on a geometric projection rule;
both ranks then compute Pelties-9 per-side flux independently. If the
two ranks agree on `can_n` (canonical normal) and the per-DOF state
update is deterministic (same psi, same traction, same slip), the
per-rank `local_slip_rate` at shared-face DOFs should match bit-for-bit.
If there's a minor drift (rounding, reduction order), the two ranks'
per-vertex VTU values would differ slightly, but since ParaView doesn't
weld across VTU files, only the rendering would differ — visually
symmetric across the partition boundary but never producing an isolated
speckle.

**Evidence against being the root cause.** Speckles appear inside the
rupture interior, not preferentially along straight lines where rank
boundaries would be.

**Diagnostic to decide.** §5.2 — overlay the rank-decomposition
boundary on the speckle locations. If speckles correlate with rank
boundaries (e.g., all speckles lie on shared-face elements), H-V91-A2
fires. If they scatter randomly, H-V91-A2 is falsified.

#### H-V91-A3 — DOFData mapping read-stale

**Claim.** `WriteFaultSurfaceVTU` uses
`base = face_mesh_idx * nbf` (line 580) to index into
`local_slip_rate`. But `face_mesh_idx` is the face index *as passed
into process_face*, which for interior faces is `i` (the same as the
loop counter) and for shared faces is `n_int + i`. If the DOF vector
layout in the driver does NOT match this ordering — e.g., if the
driver packs DOFs by `interior_face_order` that differs from
`fault_interior_faces_`, or if it includes orphan faces that bump the
counter — then VTU reads misaligned DOF values.

**Evidence.** Lines 242–264 describe the orphan handling:
`fault_face_elem1_[i] = -1` is recorded for faces where the
interior-transformation lookup returns nullptr. The docstring says
"kept to preserve fi*nbf index alignment with fault_coords". So the
orphan DOFs ARE in the local_slip_rate vector at positions that the
VTU writer would read. But the VTU writer's `process_face` early-returns
on `!FTr`, so it wouldn't emit an orphan face's value. Still,
alignment must be preserved: if orphans are present, `fault_interior_faces_[i]`
must correspond to `i * nbf` in local_slip_rate. The driver enforces
this in `fault_geometry.hpp` / `fault_nodes.hpp`. Confirm.

**Diagnostic to decide.** §5.3 — for each face in the VTU, dump the
DOF indices it reads from and the corresponding coordinate from
`local_x2`/`local_x3`. Check that the VTU triangle's physical
coordinates match the `(x2, x3)` read from the DOF vector. Any mismatch
is an H-V91-A3 hit.

---

## 3. Issue R-V91-B — hypocenter dip-channel deviation

### 3.1 Where the dip/strike split happens in SEAS-MFEM

The canonical frame is set at every fault QP in
`dynamic/wave_operator.inl:769–788`:

```cpp
real_t can_n[3], can_t1[3], can_t2[3];
for (int d = 0; d < 3; d++) {
   can_n[d]  = qpd.sign_flipped ? -qpd.normal[d]   : qpd.normal[d];
   can_t1[d] = qpd.sign_flipped ? -qpd.tangent1[d] : qpd.tangent1[d];
   can_t2[d] = qpd.sign_flipped ? -qpd.tangent2[d] : qpd.tangent2[d];
}
GodunovFlux::BuildRotation       (can_n, can_t1, can_t2, T_can);
GodunovFlux::BuildRotationInverse(can_n, can_t1, can_t2, T_can_inv);
```

`qpd.tangent1` and `qpd.tangent2` come from `FaultBasis` (see
`fault/fault_basis.hpp:398–418`):

```cpp
// strike = normalize(up × n_ref)                  (first cross product)
// dip    = strike × n_ref                         (second cross product)
// tangent1 = dip, tangent2 = strike               (Tandem convention)
```

So for TPV102 with `ref_normal = (0, -1, 0)` and `up = (0, 0, 1)`:
- `strike  = (0,0,1) × (0,-1,0) = (1, 0, 0)`
- `dip     = (1,0,0) × (0,-1,0) = (0, 0, -1)`
- `tangent1 = (0, 0, -1) = dip (down)`
- `tangent2 = (1,  0, 0) = strike (+x)`

In `FaultFaceFlux::Evaluate` (`fault_face_flux.cpp:39–201`):
- `tau1_trial` is built from `Q_plus[VY] - Q_minus[VY]` and
  `Q_plus[SXY] + Q_minus[SXY]`. Fault-local `VY`/`SXY` correspond to
  `can_t1 = dip`.
- `tau2_trial` mirror on `VZ`/`SXZ` → `can_t2 = strike`.
- `V1 = V_abs * tau1_total / (strength + eta_s * V_abs)` → dip slip rate.
- `V2 = V_abs * tau2_total / (strength + eta_s * V_abs)` → strike slip rate.

**Key property.** For pure strike-slip TPV102, `tau1_0 = 0`,
`tau2_0 = tau_ini`. If the bulk trace is quiescent (`Q_plus = Q_minus = 0`),
then `tau1_trial = tau2_trial = 0` and `tau1_total = 0`, `tau2_total = tau_ini`
→ `V1 = V_abs * 0 / divisor = 0` exactly. Dip slip rate is zero as long
as dip trial traction is zero.

### 3.2 How the dip channel can pick up a transient at breakaway

At the moment of breakaway, the bulk trace is NOT quiescent — the
rupture-front elastic wave is already being radiated back to the fault.
For the hypocentre QP, the bulk state at t ≈ 1.5 s has:
- `Q_plus[VZ]` (strike velocity): O(m/s) — the fast slip wave
- `Q_plus[VY]` (dip velocity): O(mm/s) — radiation coupling
- `Q_plus[SXY]` (dip traction): O(10 kPa) — radiation coupling

Through eq. 7b–c these feed back into:
- `tau1_trial = eta_s * (VY^- - VY^+ + SXY^+/Zs^+ + SXY^-/Zs^-)` → O(10 kPa)
- `tau2_trial = eta_s * (VZ^- - VZ^+ + SXZ^+/Zs^+ + SXZ^-/Zs^-)` → O(10 MPa)

The friction solver then computes:
- `Theta = sqrt(tau1_total² + tau2_total²)` ≈ tau2_total (dominated by strike)
- `V1 = V_abs * tau1_total / divisor` — tiny dip slip rate tracking the tiny dip trial traction
- `V2 = V_abs * tau2_total / divisor` — large strike slip rate

So some dip slip rate is expected from bulk-radiation coupling.
**The question is whether MFEM's dip coupling is identical to SeisSol's
or has an asymmetry.**

### 3.3 SeisSol analogue — for direct comparison

`SeisSol/src/DynamicRupture/FrictionLaws/CpuImpl/RateAndState.h:200–256`
(`calcSlipRateAndTraction` — called from `updateFrictionAndSlip`):

```cpp
totalTraction1 = initialStressInFaultCS[3] + faultStresses.traction1  // XY
totalTraction2 = initialStressInFaultCS[5] + faultStresses.traction2  // XZ
divisor = strength + impAndEta.etaS * slipRateMagnitude
slipRate1 = slipRateMagnitude * totalTraction1 / divisor
slipRate2 = slipRateMagnitude * totalTraction2 / divisor
traction1 = faultStresses.traction1 - etaS * slipRate1
traction2 = faultStresses.traction2 - etaS * slipRate2
```

**Structurally identical** to SEAS-MFEM `fault_face_flux.cpp:141–156`.
Modulo notation, no difference. Note that SeisSol's traction1/traction2
are face-local XY/XZ in the fault-aligned frame — independent of which
cardinal direction is "dip" vs "strike"; that mapping is applied at
output time only.

So the vector-friction law does not distinguish dip from strike.

### 3.4 Where the asymmetry could enter

#### H-V91-B1 — FaultBasis dip computed with one more cross-product than strike

`fault/fault_basis.hpp:398–418`:

```cpp
// strike = normalize(up × n_ref)           // 1 cross product + normalize
// dip    = strike × n_ref                  // 1 more cross product
// (note: no separate normalize for dip — strike _|_ n_ref is already unit)
for (int d = 0; d < 3; d++) { tangent1[d] = dv[d]; tangent2[d] = s[d]; }
```

So `tangent1 = dip`, `tangent2 = strike`. The strike vector is
normalized explicitly; the dip is computed from `strike × n_ref` which
is unit *mathematically* but can drift by O(ε) in FP due to `strike`
normalization residuals propagating.

**Hypothesis.** For each face, the dip direction numerical error is
O(ε · |strike|). For a face off the fault centerline (where strike
length is O(1)), this is O(ε). But for faces near the fault edges or
near rank partition boundaries, numerical errors in `strike` can grow
by a few ulps, and dip inherits those errors + its own cross-product
precision. At the hypocentre, where multiple faces meet at the nucleation
centre, the per-face dip vectors may have micro-inconsistencies that
accumulate into a visible dip-channel residual during breakaway.

**Reviewer scaling estimate (R-003).** `|dip|` drifts by up to ~5 ULP
= O(10⁻¹⁵). At TPV102's 300 MPa stress scale this is ~10⁻¹³ Pa — far
below the 0.3 MPa hypocentre spike. **The reviewer explicitly
concludes R-003 cannot alone explain R-V91-B.** Fix is still landed
in commit C0 as a defensive, cheap orthonormality guarantee that
removes a latent asymmetry class (§11).

**Evidence for.** Dip spike is localised to the hypocentre; off-axis
stations (where individual face dip vectors are less sensitive to
breakaway amplification) match reference.

**Evidence against.** At double precision the effect is ~10⁻¹⁵, far
below a 0.3 MPa spike (which is 10⁻³ of the 300 MPa background stress
scale). More likely to be FP-rounding level, not physical level.
Reviewer agrees.

**Diagnostic.** §4.2 unit test — `test_fault_basis_dip_strike_symmetry`
perturbs the ref_normal and up and asserts dip and strike lengths are
both 1.0 to 10 ULP, and both orthogonal to the normal to 10 ULP.
Expected pre-fix: `|strike|` passes (explicit normalize), `|dip|`
fails by ~5 ULP. Expected post-R-003: both pass. Test verifies the
fix without claiming it explains R-V91-B.

#### H-V91-B2 — T_can rotation ordering (can_t1 first vs can_t2 first)

The `BuildRotation` call passes `(can_n, can_t1, can_t2)` in fixed
order. Internally, this builds a 9×9 Voigt rotation where `can_t1` maps
to the fault-local Y axis (for VY) and `can_t2` to the Z axis (for VZ).
By the assignment `tangent1 = dip`, SEAS-MFEM effectively puts dip on
VY and strike on VZ.

**Hypothesis.** There's nothing sacred about dip-on-Y and strike-on-Z.
In a pure strike-slip problem, if the Godunov flux Jacobian `A_x` has
any asymmetry between its Y-block and Z-block (it shouldn't for
isotropic elasticity), dip and strike would pick up different errors.

**Evidence for.** Concretely, `A_x^±` for 3D isotropic elastic wave
equation is block-diagonal with identical S-wave blocks for Y and Z
tangential components. So dip-on-Y and strike-on-Z should produce
identical numerical flux for identical trial traction. No asymmetry
expected from this source alone.

**Evidence against being the root cause.** Block symmetry is exact in
theory and tested in §3.1 of v9.0.0 (12/12 PASS).

**Diagnostic.** §5.5 permutation test — swap the tangent1/tangent2
assignment (or equivalently, rotate 90° around the normal and rerun
§3.1 test). Expected: zero change in flux magnitudes at QP level.
Any difference > ULP is an H-V91-B2 hit.

#### H-V91-B3 — SEAS-MFEM uses a different friction iteration scheme than SeisSol

From `dynamic/friction_solver.{hpp,cpp}` (already validated in v7.0.0
H-V7-friction-solver): SEAS-MFEM uses **Brent** on `|V|` in log10 space
for the rate-state friction equation (Eqs. 8 of the plan). From
SeisSol's `CpuImpl/RateAndState.h:150–197` (`updateStateVariableIterative`):
SeisSol uses **Kaneko 2008 iteration** — a fixed-point scheme that
alternates state-variable updates and slip-rate solves, averaging new
and old slip rates.

**Hypothesis.** The two schemes find the same `|V|` to solver tolerance.
But their intermediate iterates differ. If the plans' `dt` is adapted
on these iterates (SEAS-MFEM has RK4 explicit; SeisSol has ADER-LTS),
the effective sub-step precision could differ between codes at
breakaway transients, producing different dip/strike timing noise.

**Evidence for.** Breakaway is inherently transient. Differences between
iterative schemes are most visible there.

**Evidence against.** The friction solver acts on `|V|` and `Theta`
only, with zero dip/strike asymmetry in the inner iteration. Both
schemes should converge to the same `|V|`, then decompose to `V1`/`V2`
via the same ratio `totalTraction1 / divisor` vs `totalTraction2 / divisor`.
No asymmetry between dip and strike can be introduced by the friction
scheme alone.

**Diagnostic.** §3.C — audit SeisSol `updateStateVariableIterative`
and compare convergence criterion and tolerance to SEAS-MFEM's Brent
solver. Look specifically for (a) a stopping criterion that differs
between dip and strike components (there should be none), (b) any
Kaneko-style averaging that smooths the transient, and (c) the
post-convergence `calcSlipRateAndTraction` call — if SeisSol applies
any extra correction / filtering after convergence, SEAS-MFEM may miss
it.

### 3.5 Read-out: these three hypotheses have very different fix
surfaces

- **H-V91-B1** (basis asymmetry): fix is in `fault_basis.hpp` — add
  an explicit normalize on dip + tighten the cross-product convention.
  Fixes a pure-numerics issue; unlikely to move a 0.3 MPa spike.
- **H-V91-B2** (rotation ordering): fix is in `wave_operator.inl` —
  swap tangent assignment. Unlikely root cause; physical asymmetry is
  absent.
- **H-V91-B3** (friction iteration): fix is either in
  `friction_solver.cpp` (swap to Kaneko) or leave Brent alone —
  depends on whether the transient deviation persists to t = 12 s.

**The 12 s Frontera run in §6 is required to distinguish H-V91-B3 from
B1/B2.** If the dip spike decays below reference envelope by t = 5 s
(transient-only), the deviation is iteration-scheme sensitivity and
does not need a fix (the codes agree in the asymptotic limit). If it
persists, H-V91-B3 is the rank-1 hypothesis.

---

## 4. Unit tests to add (laptop-only, no Frontera)

### 4.1 `test_fault_surface_vtu_continuity` — post-R-001 CellData contract

**Purpose.** (a) Confirm the R-001 fix's CellData writer produces
correct face-averaged values for a known-smooth field. (b) Regress
against H-V91-A4 (QP-at-vertex mismatch) returning. (c) Catch
H-V91-A1 as a secondary check — under a linear field the
face-averaged CellData form must produce zero inter-cell jumps.

**Scaffold.** Under `tests/unit/test_fault_surface_vtu_continuity.cpp`.
Emulate `WriteFaultSurfaceVTU` (post-R-001 CellData form) on a small
4×4 fault-only triangulation. Fill `local_slip_rate` from a known-smooth
linear analytic field
`f(x2, x3) = 0.5 + 0.1 · x2 + 0.2 · x3` evaluated at each QP. Run the
VTU writer to a temp directory and parse the output.

**Acceptance targets.** Per reviewer R-001 test case:
- **Face-average correctness.** For each output triangle, the single
  CellData value must equal the face-centroid evaluation of `f`
  (= average of `f` at the three QPs by linearity of barycentric
  quadrature) to `1e-10` relative.
- **Inter-cell continuity under linear field.** For each pair of
  triangles sharing a physical edge, cell-value difference ≤ `1e-12`.
- **All-QPs-reach-output** (R-002 regression): with order=2 (nbf=6)
  and per-QP IDs `local_slip_rate[2*d] = d`, each face's cell value
  equals the mean of its six per-QP IDs (globally unique), so a
  pre-fix writer that reads only the first three QPs would return a
  distinguishable different value.

**Pre-fix behaviour (sanity).** Against the unpatched
`paraview_output.hpp`, the same test must fail the face-average
assertion by `O(h · ||∇f||)`; this is the regression gate for R-001.

**Makefile entry.** After `seas_test_fault_face_flux_frame`. Target
`test-fault-surface-vtu-continuity`.

### 4.2 `test_fault_basis_dip_strike_symmetry` — post-R-003 orthonormality contract

**Purpose.** (a) Confirm the R-003 fix (explicit normalize of dip in
`ComputeOrientedFrame`) produces an orthonormal basis to 10 ULP on
random perturbed normals. (b) Document the pre-fix ~5 ULP dip-length
drift as a regression gate.

**Scaffold.** Under `tests/unit/test_fault_basis_dip_strike_symmetry.cpp`.
Construct a minimal 2-element mesh with a single shared fault face in
the `(0, -1, 0)` plane. Compute `FaultBasis` with
`ref_normal = (0, -1, 0)`, `up = (0, 0, 1)`.

**Acceptance targets.** Per reviewer R-003 test case:
For 1000 random unit-vector perturbations of the normal
(`||delta|| < 0.1`, normalized), and the basis computed via
`FaultBasis::Compute`:
- `|strike| == 1.0` to 10 ULP (pre-fix: passes via explicit normalize).
- `|dip| == 1.0` to 10 ULP (pre-fix: fails by ~5 ULP; post-R-003: passes).
- `dot(dip, strike) == 0.0` to 10 ULP.
- `dot(dip, normal) == 0.0` to 10 ULP.
- `dot(strike, normal) == 0.0` to 10 ULP.

**Note.** The test deliberately does NOT assert that R-003 by itself
moves the R-V91-B hypocentre dip spike — the reviewer's scaling
estimate places that well below the observed 0.3 MPa. The test
guarantees orthonormality so downstream rotation matrices are
well-conditioned, eliminating this as a latent issue for any future
higher-order / smaller-dt configuration.

**Makefile entry.** After `seas_test_fault_basis`. Target
`test-fault-basis-dip-strike-symmetry`.

### 4.3 Regression preservation

All v9.0.0 tests (§3.1, §3.1b, §3.1c, §3.1d) must continue to pass.
The two new tests do not touch production code; run `make test` after
adding to confirm no regression.

---

## 5. Offline diagnostics against job 7667881 (no Frontera)

These run locally on the laptop against the already-downloaded output
directory.

**Status after review.** §5.1 / §5.2 / §5.3 were written as
diagnostics to identify the R-V91-A root cause. The review's R-001
identifies that root cause directly from source reading (H-V91-A4),
so these three diagnostics become **verification** rather than
discovery — they now confirm that the R-001 patch closes the speckle
artefact. §5.4 / §5.5 remain primary diagnostics for R-V91-B, which
is not addressed by any of the R-001–R-004 patches.

### 5.1 VTU speckle audit

**Input.**
`tpv102/plots_results_200m_p1_2.0s_400r_vmax_halt_job7667881/FaultSurface/fault_surface_r0_c10.vtu`
(or any mid-run cycle). Confirm this directory exists locally first
(the `plots_results_*` dir looks like a post-processed dataset; the
raw VTUs may live in a `run_*` dir).

**Procedure.**
1. Parse the VTU with any XML parser (Python/xml.etree in a small
   script under `tpv102/verification/`).
2. For each triangle i and each vertex v, record `(vertex_coord, value_at_vertex)`.
3. Identify pairs of triangles (i, j) that share a physical edge — two
   vertices of i are co-located (within 1e-10) with two vertices of j.
4. At each shared vertex, compute `|value_at_i - value_at_j|`. This is
   the DG jump.
5. Plot a histogram of DG jumps. Overlay on ParaView rendering —
   speckle locations correspond to top 5% of DG jumps.

**Expected outcome.**
- If speckles = top-of-histogram DG jumps → H-V91-A1 is the cause.
- If speckles have DG jumps in the bottom half of the histogram →
  speckles are NOT DG jumps; they are stale DOF values or a writer
  bug. H-V91-A3 prime.

### 5.2 Partition-boundary overlay

**Procedure.**
1. Parse each `fault_surface_r{r}_c{c}.vtu` per rank.
2. Record a "rank ownership" field for each triangle (which VTU file
   contains it).
3. Render in ParaView with rank-colouring alongside the slip-rate
   colouring. Compare speckle locations to rank-boundary locations.

**Expected outcome.**
- Speckles on rank-boundary elements → H-V91-A2 fires (partition
  mismatch).
- Speckles inside a single rank's domain → H-V91-A2 falsified.

### 5.3 DOFData-zero audit

**Purpose.** Check whether any DOFs in `local_slip_rate` at the time
of VTU write had not been updated in the previous RK step, i.e.
carried stale values.

**Procedure (requires DIAG build).** Add a compile-time diag block in
the driver that writes to `debug_dof_update_log_r{r}.txt` per rank:
per-DOF, last RK step at which the DOF was written. At VTU write time,
compare "last DOF update step" to current step; log any DOF that has
`step_diff > 1` (i.e., missed the last step's update).

**Acceptance.**
- Zero DOFs with `step_diff > 1` → H-V91-A3 falsified.
- Any non-zero → H-V91-A3 confirmed; bug in DOF update loop.

**Cost.** Requires rebuild with `-DSEAS_DIAG_FAULT_DOF_UPDATES=1`. Not
on the critical path if §5.1 resolves.

### 5.4 Station-hierarchy audit for dip

**Input.** The existing station plots at `flt_0_3`, `flt_0_7.5`,
`flt_0_12`, `flt_9_7.5`, `flt_n9_7.5`, `flt_12_3`, `flt_12_12`,
`flt_n12_3`, `flt_n12_12`. (All 10 station plots exist as PNGs.)

**Procedure.**
1. Extract the peak `|V_dip|` and peak `|Shear_dip|` at each station
   during the breakaway transient (t ≤ 2 s).
2. For each station, compute the ratio `(MFEM peak) / (DR peak)`.
3. Plot ratio vs `sqrt(x2² + (7.5 - x3)²)` (distance from hypocentre).

**Expected outcome.**
- Ratio monotonically decreases to 1.0 with distance → deviation is
  localised to hypocentre only, consistent with H-V91-B1 / B3.
- Ratio constant or increasing → deviation is NOT hypocentre-local;
  look for a global asymmetry.

### 5.5 Permutation diagnostic for T_can rotation

**Purpose.** Verify H-V91-B2 claim that tangent1/tangent2 ordering does
not affect the flux.

**Procedure (laptop unit test — extend §3.1 of v9.0.0).** At QP level,
for identical Q_plus/Q_minus/DOFData, compute F_h with:
- Case A: `BuildRotation(can_n, can_t1=dip, can_t2=strike)` (current).
- Case B: `BuildRotation(can_n, can_t1=strike, can_t2=dip)` (swapped).

After swapping inputs to `Evaluate` (tau1↔tau2 roles swapped), F_h in
both cases should be identical. Any difference > ULP = H-V91-B2 hit.

---

## 6. Frontera 12 s re-run specification (not yet launched)

### 6.1 Rationale

tfinal = 2.0 s of job 7667881 shows the breakaway transient and the
first ≈ 0.5 s of decay. To decide H-V91-B3 (friction-iteration
sensitivity) vs H-V91-B1/B2 (structural asymmetry) we need to see
whether the dip-channel transient in MFEM decays to reference envelope
asymptotically. SeisSol's dip noise at flt_0_7.5 in the reference
datasets continues past t = 10 s with ±0.003 m/s. MFEM needs the same
tfinal to be compared.

The strike channel is a secondary check — we expect visual agreement
with SeisSol to tfinal = 12 s since the v9.0.0 fix lines up the rupture
dynamics.

### 6.2 Proposed sbatch parameters

| Parameter | Value | Rationale |
|---|---|---|
| Mesh | `tpv102/mesh/tpv102_200m.msh` | Match job 7667881 |
| Order | P1 | Match job 7667881 |
| Ranks | 400 | Match job 7667881 |
| Nodes | 8 | Match job 7667881 |
| tfinal | **12.0 s** | v9.1.0 — extended from 2.0 |
| CFL | 0.5 | driver auto-CFL for explicit RK4 (per `feedback_explicit_cfl.md` don't set `--dt`) |
| Output dt | 0.1 s | 120 frames for the VTU series, roughly 600 MB of ParaView data |
| Station output | every step or 0.005 s | station time-series — small; existing driver flag `--station-dt` |
| Time limit | 4 hours | ballpark 6× the 2 s run |
| DIAG flags | None (v9.0.0 bisection is closed) | keep output small |

**Build flags.** No `SEAS_DIAG_*`. Match the working sbatch module list
from `feedback_sbatch_modules.md` — do NOT hand-prune `module load`
lines (fftw3, etc.).

**Halt condition.** User-configured V_max halt at 400 m/s (same as
job 7667881 — filename `vmax_halt`). Do NOT reduce this, we want to see
rupture saturation.

### 6.3 Acceptance after the 12 s run

#### Primary — R-V91 close criteria

- **R-V91-B-pass-1:** `flt_0_7.5` Slip_dip stays within ±2× SeisSol
  reference envelope (i.e., |MFEM slip_dip| ≤ 0.002 m) for t ≥ 5 s.
  Breakaway transient may still exceed this during 1.4 s ≤ t ≤ 2.5 s
  but must decay.
- **R-V91-B-pass-2:** `flt_0_7.5` V_dip RMS past t = 5 s matches
  SeisSol reference to within 2× (both should be ~0.003 m/s RMS).
- **R-V91-B-pass-3:** `flt_0_12`, `flt_0_3` dip channels agree with
  reference to within 2× (same envelope as SeisSol).

If all three pass → H-V91-B3 confirmed as the mechanism (transient
iteration sensitivity; both codes converge asymptotically). R-V91-B
closes with no further code change (R-003 has already landed in C0).

If any fail → H-V91-B1 or B2 become prime. Open a follow-on v9.1.1
plan and add a targeted unit test using the production data to bound
the asymmetry.

#### Secondary — VTU speckle visual (R-V91-A post-R-001 sanity)

With the R-001 CellData patch in the 12 s build, the fault-surface
VTU time series must be speckle-free from t = 0 onward. Any residual
speckle would fire H-V91-A1 (genuine DG jumps in the CellData form)
or indicate an incomplete R-001 patch. If speckles persist in the 12 s
output: re-open R-V91-A with a fresh hypothesis ladder; do not claim
R-V91-A closed.

---

## 7. Assumptions and caveats

1. The reference "DR" curves in the station plots are from the SeisSol
   DR benchmark dataset (`plots_results_*/tpv102_flt_*.png` — paper
   reference) — NOT from a fresh SeisSol run. If fidelity is doubted,
   run SeisSol locally or request the Christodoulou 2023 dataset. Out
   of scope for v9.1.0.
2. PyLith reference noise envelope is distinct from SeisSol's; we use
   SeisSol as primary reference since the v9.0.0 flux derivation is
   grounded in Pelties 2012 (SeisSol paper).
3. The 12 s run disk budget: 120 cycles × ~30 MB per VTU snapshot ≈
   3.6 GB for ParaView fault-surface data. Add station CSV (~50 MB).
   Confirm Frontera quota before launch.
4. The v9.0.0 Pelties-9 patch is assumed correct — no regression tests
   were added to the quasi-dynamic BP5 pathway that also uses
   `wave_operator.inl`. If v9.0.0 inadvertently broke BP5, it would
   appear only on a BP5 rerun. Out of scope for this plan; flagged as
   a suggestion for a separate BP5 regression sweep.

---

## 8. Commit plan — in execution order, once you approve

| # | Commit scope | Branch state after | CI expected |
|---|---|---|---|
| C0 | Source fixes from review: R-001 (`paraview_output.hpp` CellData rewrite, subsumes R-002 loop-bound), R-003 (`fault_basis.hpp` explicit dip normalize), R-004 (docstring update). All four findings in one commit so fault-surface + fault-basis code land together for a single CI round. | GREEN | passes |
| C1 | Add `tests/unit/test_fault_surface_vtu_continuity.cpp` + Makefile target + docstring. Written against the post-R-001 CellData contract; verifies C0's R-001 fix and regresses against H-V91-A4 regression. | GREEN | passes |
| C2 | Add `tests/unit/test_fault_basis_dip_strike_symmetry.cpp` + Makefile target + docstring. Written against the post-R-003 orthonormality contract. | GREEN | passes |
| C3 | Add `tpv102/verification/vtu_speckle_audit.py` (§5.1/5.2 offline diagnostic, pure post-processing). Post-R-001 this becomes a verification tool rather than a discovery diagnostic. | GREEN | N/A (script) |
| C4 | After §5 + §4 outcomes: either (a) CLOSE R-V91-A with a `tpv102_debug_v9.1.0_fix.md` documenting R-001 as the cause and fix, or (b) open v9.1.1 plan if a residual bug is found post-C0. | GREEN | passes |
| C5 | Generate `jobs/tpv102/tpv102_200m_p1_12.0s_400rank_v91.sbatch` — Frontera-only, not auto-submitted. | GREEN | N/A (sbatch) |
| C6 | After user-approved Frontera submission + results: CLOSE R-V91-B with fix doc, or open v9.1.1 plan. | GREEN | passes |

**C0 is the only commit with production-code changes.** C1/C2 add
tests; C3 adds post-processing; C4/C5/C6 are documentation and
Frontera artefacts. All seven commits individually reversible; C0
must land before C1/C2 so the tests exercise the fixed code.

---

## 9. Glossary / references

- **R-V91-A / R-V91-B:** regression IDs for the two issues raised in
  this plan (v9.1.0 issue labels).
- **R-001 … R-004:** review-finding IDs from `REVIEW.md` (2026-04-20).
  See §11 for mapping to plan sections.
- **H-V91-A4:** hypothesis IDs are plan-level; A4 was added post-review
  and is now RANK 1 for R-V91-A.
- **Pelties-9:** Pelties et al. 2012 eq. (9) — the per-side flux form
  SeisSol uses at fault faces. Applied in SEAS-MFEM as of v9.0.0.
- **BP5 convention:** tangent1 = dip, tangent2 = strike. Established in
  `CLAUDE.md` (seas miniapp) and `fault_basis.hpp` docstring.
- **`qpd.sign_flipped`:** per-QP flag set by FaultBasis when raw normal
  `n_raw` is anti-aligned with ref_normal; all basis vectors are negated
  so downstream code never applies a sign factor.
- **Station naming:** `flt_{x2_km}_{x3_km}` — `x2` along strike in km,
  `x3` down-dip depth in km. `flt_0_7.5` = hypocentre (0 km along
  strike, 7.5 km depth).
- **DR / PyLith reference:** `DRDG3D` (Pelties/SeisSol DR benchmark)
  and `PyLith` — the SCEC TPV102 community benchmark reference curves.
- **Screenshot referenced:** `~/Desktop/Screenshot 2026-04-20 at 7.55.14 PM.png`
  — ParaView rendering of fault_surface.pvd slip_rate_strike field at
  a mid-run cycle of job 7667881.
- **MFEM triangle degree-2 rule:** `IntRules.Get(TRIANGLE, 2)` →
  `AddTriPoints3(0, 1./6., 1./6.)` → the three barycentric points
  `(1/6, 1/6), (2/3, 1/6), (1/6, 2/3)` — all strictly interior to the
  reference triangle. Defined in `fem/intrules.cpp:1271-1277` with
  explicit `// interior points` comment. Central to H-V91-A4.

---

## 10. Unreviewed / carry-over items from REVIEW.md §Unreviewed Areas

Flagged by the reviewer as explicitly not covered in the 2026-04-20
pass, but worth tracking:

- **`friction_solver.hpp` inner Brent routine.** Out of scope
  for v9.1.0; behaviour is the open H-V91-B3 question. Existing tests
  `test_friction_solver.cpp` + `test_vector_friction.cpp` pin the
  bracket handling. No action this plan.
- **PML damping.** Applied after the fault-face loop in
  `wave_operator.inl:477-480`. Reviewer confirms a PML bug cannot
  produce hypocentre-localised dip noise. No action.
- **Driver RK4 loop `tpv102_driver.cpp:777-890`.** The prior
  2026-04-14 review found an output-consistency issue with `tau`/`V`
  fields at station-writer time (RK4-weighted V1/V2/tau_corr averaging
  at lines 858–874). That finding's suggested diff still stands and
  is independent of R-V91-A and R-V91-B. **This plan does not
  supersede that finding** — reviewer explicitly carries it open.
  Recommended action: address in a separate v9.1.2 ticket; do not
  batch with this v9.1.0 work.
- **`pseas.cpp` / `seas_pseas` / BP5 drivers.** Out of scope per plan
  §7 caveat (4). BP5 regression sweep against the v9.0.0 Pelties-9
  patch is proposed there but deferred.
- **`TPV102StationWriter::Write`** (`tpv102_setup.hpp:300-315`). Uses
  the same `DOFData.V1/V2/slip1/slip2/tau1_corr/tau2_corr` fields as
  the VTU writer and is covered by the 2026-04-14 output-consistency
  finding. No separate review.

---

## 11. Review incorporation index (REVIEW.md, 2026-04-20)

This plan was reviewed in
`/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/REVIEW.md` on
2026-04-20. Verdict: **PASS WITH FIXES**. Summary of findings and how
they have been incorporated into this plan.

### 11.1 Findings

| ID | Severity | File | Summary | Plan-level consequence | Source-level fix |
|---|---|---|---|---|---|
| **R-001** | **CRITICAL** | `paraview_output.hpp:566-596` (`WriteFaultSurfaceVTU`) | VTU writer places per-QP values at reference-triangle **vertex** positions while QPs live at interior barycentric `(1/6,1/6), (2/3,1/6), (1/6,2/3)`. Root cause of R-V91-A. | Added H-V91-A4 to §2.3 as RANK 1. Demoted H-V91-A1 to RANK 2. Updated §2.1 code-path analysis. Updated §4.1 test contract. Updated §5 preamble. Updated §6.3 secondary-acceptance. | Reviewer diff in `REVIEW.md` R-001 §Suggested fix — PointData→CellData at lines 540–596 and 696–716 of `paraview_output.hpp`. Averaging loop `k<nbf` (subsumes R-002). Land in commit **C0**. |
| **R-002** | MODERATE | `paraview_output.hpp:581` (same file) | Hardcoded `for (int k = 0; k < 3; k++)` silently drops QPs beyond 3 (order>1 or quad faces). Latent. | Documented in §2.1. §4.1 test includes a nbf=6 / order=2 sanity check. | Subsumed by R-001 averaging loop — no separate diff. |
| **R-003** | LOW | `fault_basis.hpp:413-418` (`ComputeOrientedFrame`) | `dip = strike × n_raw` not explicitly normalized; ~5 ULP drift from orthonormality. | Updated §3.4 H-V91-B1 with reviewer's scaling estimate (10⁻¹⁵ — too small to explain R-V91-B). Updated §4.2 test contract. | Reviewer diff in `REVIEW.md` R-003 §Suggested fix — insert `d_len = sqrt(...)`, `MFEM_VERIFY(d_len > 1e-12, ...)`, `dv[] *= 1.0/d_len` between the cross product and the tangent1/tangent2 assignment. Land in commit **C0**. |
| **R-004** | LOW | `paraview_output.hpp:505-510` + `:538` (docstring) | Docstring claims "smooth interpolation within each face in ParaView" but writer emits QP values at vertex positions. | Noted in §2.1 for future readers. | Reviewer diff in `REVIEW.md` R-004 §Suggested fix — rewrite docstring to describe post-R-001 CellData form. Land in commit **C0**. |

### 11.2 Plan-compliance verdict

Reviewer confirms:
- v9.0.0 Pelties-9 per-side flux in `wave_operator.inl:832-918`
  (interior) and `:1326-1340` (shared) is self-consistent with this
  plan's §3.1 math and with `fault_face_flux.cpp`'s (7)/(9)/(11)-(12)
  pipeline. **No new plan-deviation bugs in the solver.**
- Plan's §2 hypothesis list was incomplete — H-V91-A4 needed to be
  added before §5.1's audit ran. **Addressed in this revision.**
- R-V91-B is not addressable by any bug in the code as written; §6's
  12 s re-run is still required to distinguish H-V91-B3 from B1/B2.
  **Plan unchanged on this point.**

### 11.3 Revision log

- **2026-04-20 (initial):** Plan authored against station + screenshot
  data from job 7667881, with H-V91-A1/A2/A3 and H-V91-B1/B2/B3.
- **2026-04-20 (post-review, this revision):** H-V91-A4 added as
  RANK 1 for R-V91-A; §2.1 rewritten to cite MFEM integration-rule
  table; §4.1 and §4.2 test contracts rewritten to match reviewer's
  proposed post-fix assertions; §8 gains C0 (source-patch commit)
  before C1/C2 (test commits); §10 captures reviewer's unreviewed
  areas including the open 2026-04-14 tau/V-consistency finding;
  this §11 added as the review index.

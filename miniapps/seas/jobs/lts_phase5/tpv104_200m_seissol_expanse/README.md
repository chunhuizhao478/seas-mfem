# LTS Phase 5 — SeisSol side of the TPV104-200m head-to-head (Expanse)

The missing third number of the Phase-5 comparison: **SeisSol v1.1.3 o4 on the
exact same 200 m mesh** as the MFEM runs in `../tpv104_200m_lts_speed_expanse/`.

| leg | scheme | status |
|---|---|---|
| MFEM GTS | p3 / ADER-O4, global stepping | measured: **~0.27 sim-s/h** (job 52332649, canceled at t=1.52 s) |
| MFEM LTS | p3 / ADER-O4, clustered rate-2 (λ-wiggle, Nc=6) | job 52340525 |
| **SeisSol LTS** | o4 (=p3), clustered rate-2 (`ClusteredLTS=2`) | **this deck** |

## Provenance (verified 2026-07-21, workflow recon)

- `tpv104.puml.h5` (103 MiB, **not in git** — over GitHub's 100 MB limit, ships by scp)
  is pumgen's conversion (`pumgen -s msh2`) of a **byte-for-byte copy** of the MFEM mesh
  `miniapps/seas/tpv104/mesh/tpv104_200m.msh` (`.msh` md5 `e88ae223205530548b1d3438e6eb17c1`).
  Verified: 2,464,689 tets (= MFEM NE exactly), 401,211 vertices, BC decode
  {1=free: 25,142, 3=DR: 75,808 (÷2 = **37,904 fault triangles** = MFEM's fault-face count), 5=absorbing: 8,006}.
  puml md5 pinned in the sbatch: `5911b0b4ec29a8f61968d5b32ed0a11c`.
- `parameters.par` / `tpv104_fault.yaml` / `tpv104_material.yaml`: from the QuakeWorx
  tpv104 deck (`~/Downloads/tpv104/`). Full parameter-parity audit vs the MFEM deck:
  **every quantity matches** (a/b/Dc/f0/V0/f_w/V_w boxcars, σn=120 τ=40 MPa, λ/μ/ρ,
  compact-bell nucleation R=3000 m Δτ=45 MPa at (0,0,−7500), tanh transitions
  algebraically identical) — **except `t_0`, which was 0.5 s (stale) and is fixed to
  1.0 s here** (MFEM `T_nuc_s="1.0s"`; the fault.yaml header itself documents 1 s).

## Deck deltas vs the QuakeWorx original (all in `parameters.par` header too)

`t_0` 0.5→**1.0** · `EndTime` 12→**2.0** · fault output 0.5→**0.05 s** + `refinement` 1→**0** ·
surface output 0.5→**0.05 s** · `OutputPointType` 4→**5** + `&Pickpoint` filled (9 SCEC on-fault
stations) · off-fault receivers re-enabled (`pickdt=0.005`, 6 stations at z=−0.1 m) ·
`numfluxnearfault` line **removed** (needs ≥v1.2.0; binary is v1.1.3).
**Unchanged:** `CFL=0.5`, `ClusteredLTS=2` — each code runs its own production scheme.

## Fairness contract (incl. adversarial-verify disclosures, 2026-07-21)

- **Same mesh** (md5-pinned), **same order** (o4 = degree-3 basis, order factor 2·4−1 = 7 ==
  MFEM p3+ADER-4, factor 2·3+1 = 7; same h = 2×insphere formula in both codes), same physics
  window (2.0 s, t_0 = 1.0 s; all three legs), same node count (2 × 128-core), both fp64,
  comparable output cadence (0.05 s fault+surface), no volume output on either side.
- **Each code its own production time-stepping**: SeisSol `CFL=0.5` + auto-cluster LTS
  (v1.1.3 defaults, no wiggle); MFEM `cfl_dg_safety=1.0` + λ-wiggle Nc=6.
- **Courant number MATCHED at 0.5 on both codes.** The MFEM decks originally inherited
  TPV104's conservative `cfl=0.25` (half SeisSol's dt → a silent 2× step handicap); both MFEM
  decks are now `cfl=0.5`, so with the identical h formula + order factor the two codes take
  the SAME dt. Basis: SAFS production MFEM p3 upwind+ADER ran 0.5 at safety 1.0; the 2 s
  window incl. rupture is the stability gate. ⚠️ MFEM runs *before* the change (52332649 GTS,
  52340525 LTS, both cfl=0.25) keep a valid LTS-vs-GTS *ratio*; their absolute sim-s/h
  convert ×2 for cross-code comparison.
- **Known second-order asymmetries** (each ≲ a few %, disclosed not fixed):
  SeisSol's family config idles the reserved comm core (`SEISSOL_COMMTHREAD=0` with OMP=15
  → 120/128 cores busy, ~6% self-handicap — kept because it IS the proven production config);
  the MFEM binary carries Caliper instrumentation when built `USE_CALIPER=YES` (~1–3%, run
  `P5_CALI=0` for a clean A/B); MFEM writes a final V2 checkpoint (SeisSol Checkpoint=0);
  SeisSol's 9 on-fault ascii receivers sample every local step vs MFEM's per-sync stations
  (sub-1%, roughly offset by MFEM's 256-small-files-per-snapshot VTU pattern on Lustre).
- SeisSol runs 8 MPI × 15 OMP per node; MFEM 128 pure-MPI ranks per node. Per-node
  throughput (sim-s/wall-h on 2 nodes) is the comparison metric.

## Run it

Deployed like the other SeisSol decks — **rsync'd to `/expanse/lustre/projects/lbl107/czhao1/safs/`**
(alongside `v4_0_0_pref_case1_small_pgv` etc.), NOT run from the seas-mfem repo tree.
This repo folder is the source of truth; re-rsync after any edit.

```bash
# 1) sync the deck + the mesh (from the laptop):
rsync -avP /Users/chunhuizhao/projects/seas-mfem-spatial-dyn-driver/miniapps/seas/jobs/lts_phase5/tpv104_200m_seissol_expanse/ \
      czhao1@login.expanse.sdsc.edu:/expanse/lustre/projects/lbl107/czhao1/safs/tpv104_200m_seissol/
rsync -avP /Users/chunhuizhao/Downloads/tpv104/tpv104.puml.h5 \
      czhao1@login.expanse.sdsc.edu:/expanse/lustre/projects/lbl107/czhao1/safs/tpv104_200m_seissol/

# 2) on Expanse:
cd /expanse/lustre/projects/lbl107/czhao1/safs/tpv104_200m_seissol
sbatch run_tpv104_200m_seissol_expanse.sbatch        # account usc143 (the qstore binary's); or: sbatch -A lbl107 ...
tail -f tpv104-ss-o4.*.out
```

No build step: the proven v1.1.3 o4 binary at
`/expanse/projects/qstore/usc143/qwxdev/apps/expanse/rocky8.8/seisol_build/SeisSol/build-rome-o4-elastic/seissol-elastic-o4-f64`
is used as-is. The sbatch pre-flights the mesh md5, `t_0=1.0`, `EndTime=2.0`, and the
absence of `numfluxnearfault` before launching.

## Readout

1. **Speed**: the sbatch prints `SPEED METRIC: ... sim-s/wall-hour` at the end (init included;
   subtract the init wall time from the log for the stepping-only rate). Compare with the MFEM
   GTS (~0.27) and MFEM LTS numbers on the same 2 nodes.
2. **Receivers**: v1.1.3 **silently drops** unlocated receivers — the post-run block counts
   them; expect 6 off-fault + 9 on-fault.
3. **LTS clustering**: the log's cluster table shows what SeisSol's auto-clustering chose on
   this mesh (MFEM's report: RAW harmonic 3.63×, production λ=0.63 → true ≈3.68× — SeisSol's
   ideal is the same RAW number since the mesh is identical).
4. **Physics**: on-fault receiver traces + fault xdmf vs the MFEM stations/fault-VTU through
   t = 1.52 s (the GTS reference window) — same rupture, three codespaths.

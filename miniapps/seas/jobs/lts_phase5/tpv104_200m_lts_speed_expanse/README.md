# LTS Phase 5 — TPV104 200 m realized-speedup + benchmark deck (Expanse)

The **go/no-go measurement** for clustered LTS: does the machinery that was
built and byte-verified in Phases 0–4b actually deliver the projected speedup on
a production mesh, and does it still reproduce the SCEC TPV104 benchmark?

## What it runs

Two runs that differ in **exactly one config line** (`[numerics].lts`), asserted
by the sbatch before launch:

| Run | config | scheme |
|-----|--------|--------|
| GTS baseline | `tpv104_200m_lts_off.toml` | `lts="off"` — global time stepping |
| LTS | `tpv104_200m_lts_rate2.toml` | `lts="rate2"` — clustered rate-2 LTS (bulk + fault-half interleave) |

- **Mesh:** `tpv104/mesh/tpv104_200m.msh` (2.46 M tets, SCEC production resolution), via the in-folder symlink `mesh_tpv104_200m.msh`.
- **Numerics:** `order=3` / `ader_order=4`, pure upwind (`mixed_flux="none"`), `--cfl-dg-safety 1.0` (SeisSol-equivalent dt). Order-matched to SeisSol o4.
- **Perf levers ON for both** (`--deriv-cache --shared-ck-recursion`) — this is a throughput comparison, not a bit-exact gate.
- **`--lts-report` pre-step** (np=1) prints the predicted harmonic speedup and the hypothetical 256-rank partition imbalance before the timed runs.
- 256 ranks (2 nodes × 128), pure MPI. TPV104 has **no sidecars** (analytic material).
- ⚠️ The 200 m mesh (`tpv104/mesh/tpv104_200m.msh`, ~126 MB) is `*.msh` → **gitignored**; it does **not** arrive via `git clone/pull`. Stage it onto Expanse (scp/rsync into `tpv104/mesh/`) before submitting. The sbatch re-creates the in-folder symlink automatically and fails loudly if the mesh is absent.

## Prerequisite: build on Expanse (login-node compiles are forbidden)

```
bash <root>/build_expanse.sh                    # bit-exact / no Caliper
USE_CALIPER=YES bash <root>/build_expanse.sh    # for the per-component report
```
The binary must carry **LTS Phase 0–4b** (clustering + LTS-aware partition +
bulk rank-seams + fault-half interleave + the `[lts-partition]` imbalance
diagnostic) — i.e. built from the `safs-v4_0_0-alt-case1-mfem-speed` branch at
or after commit `6e3eb2b`.

## Submit

```
cd <root>/miniapps/seas/jobs/lts_phase5/tpv104_200m_lts_speed_expanse
sbatch run_tpv104_200m_lts_speed_expanse.sbatch
```

### Env knobs (all optional)

| var | default | meaning |
|-----|---------|---------|
| `P5_TFINAL_GTS` | `2.0s`  | GTS window — only a steady per-step rate is needed (the full 12 s GTS run is ~6× slower than LTS and unnecessary for the rate). |
| `P5_TFINAL_LTS` | `12.0s` | LTS window — full TPV104 rupture (benchmark stations + rupture speed) *and* steady rate. |
| `P5_RUN_GTS` / `P5_RUN_LTS` | `1` / `1` | toggle either run. |
| `P5_DERIV_OPT` | `--deriv-cache --shared-ck-recursion` | perf levers (`""` = off). |
| `P5_CALI` | `1` | Caliper per-component report (soft-degrades if the binary lacks libcaliper). |
| `P5_FAULT_INTERLEAVE` | `1` | `SEAS_LTS_FAULT_INTERLEAVE` for the RATE2 run (fault-half LTS on). |
| `P5_OUT_BASE` | `<folder>/runs` | output root (job-ID-keyed). |

## What to read out (the four comparisons)

1. **Realized LTS speedup** = `(LTS sim-s/wall-hr) / (GTS sim-s/wall-hr)` from the two run logs.
   - Acceptance gate: **≥15×** over the GTS-safety-1 baseline. Ideal: **38.7×**.
   - The `--lts-report` line gives the *predicted* harmonic speedup (≈6.0× element-update on the 200 m clustering) for cross-check.
2. **Realization attribution** — the `[lts-partition]` line (also logged live by the RATE2 run):
   - measured locally at np=256: **bulk per-cluster imbalance ~1.16×** (multi-constraint working), **fault imbalance ~2.34×** (friction-work gap, small on TPV104).
   - Caliper `Max ≫ Avg` per region confirms where any residual imbalance sits.
3. **vs SeisSol** — compare the LTS throughput to the **SeisSol o4 log on the same 200 m mesh**.
   - ⚠️ **A SeisSol TPV104-200m o4 timing log must be produced/located for this** — it is not in this repo. Point to it here once available.
4. **Benchmark physics** — diff the LTS on-fault stations (auto-written, `tag=tpv104`) against the **SCEC TPV104 reference** (rupture time, peak slip rate, final slip). This validates the fault-half interleave at full active rupture — the one piece only a production run can confirm.

## Expected outcome (Phase-0 projection)

Order-matched, the plan projects MFEM-LTS at **~1.2–1.5× faster than SeisSol o4**
on the same mesh (up to ~3× if a lower far-field order is later accepted). This
deck produces the single measured number that confirms or corrects that, and —
via `[lts-partition]` + Caliper — attributes any shortfall to load imbalance vs
communication so the next lever (e.g. the fault-work partition constraint for
SAFS) is unambiguous.

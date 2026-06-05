# Code Review: SAFS RS-SRW + LSW transfer + dev-queue jobs (2026-06-01)

## Review Scope
- Task: verify the SAFS problem-setup transfer from `seas-mfem-safs` (branch `safs`) into this repo is correct, and build dev-queue 2 h sbatch jobs for both physics (slip-weakening + rate-and-state).
- Files reviewed:
  - `safs/project_7.0_alternative/config/spatial_friction_rate_state_safs_projected_stress_srw_Dc010_nuc8km_500m.toml` (ported + 3 schema renames)
  - `safs/project_7.0_alternative/config/spatial_friction_slip_weakening_safs_projected_stress_resolution_Dc2.toml` (pre-existing)
  - `safs/project_7.0_alternative/friction/rate-and-state/param_a_500m_strongvw.csv`, `param_a_minus_b_500m_strongvw.csv` (ported)
  - `jobs/safs/spatial_dyn_ratestate_depthprofile_normalcap_pureupwind_triq_{10N_500r_normal_48hr,8N_400r_dev_2hr,8N_400r_dev_2hr_restart}_safs.sbatch` (converted to SRW)
  - `jobs/safs/spatial_dyn_slipweakening_normalcap_pureupwind_triq_{10N_500r_normal_48hr,8N_400r_dev_2hr,8N_400r_dev_2hr_restart}_safs.sbatch` (already correct)
- Domain context: `miniapps/seas/CLAUDE.md` (Gmsh v2.2, HDF5/ZFP I/O, no-local-mesh-runs), `spatial/code/spatial_friction.cpp` (config parser), project memory ([[safs-srw-port-from-safs-branch]], [[safs-vs-hrsref-config-api-divergence]]).

## What was verified correct
- **SRW config parses in this repo.** Every `[friction.rate_state]` key after the 3 renames maps to a key the parser reads; `eta="auto"` (parser line 52/694) and `[time].dt_initial="auto"` (parse_dt_initial, line 384) are accepted; the base-schema key set is byte-identical to an aging config this repo already parses (no unknown-key abort). `state_evolution="slip_law_srw"` is accepted (parser line 671); `f_w_default`/`V_w_default` are read (lines 663/664).
- **Floor pre-flight passes:** SRW config has `[friction].sigma_n_strength_floor_pa = 10.0e6` (the `^…sigma_n_strength_floor_pa…=` grep matches); LSW config has it at line 199.
- **All referenced inputs present:** SRW config, LSW config, both `_strongvw` CSVs, and the 500 m `_triq` `.msh` (62 MB, gitignored) are on disk.
- **Job↔config consistency:** all 3 RS jobs reference the SRW config + `_strongvw` CSVs with zero stray `depthprofile.toml`/`vwvs11km` references; all 3 SW jobs reference the LSW config.

## Findings

### [R-001] MODERATE [jobs/safs/*_ratestate_*_dev_2hr*.sbatch] — RS dev + dev-restart jobs were still on the AGING config (incomplete transfer) — FIXED THIS SESSION

**Category:** DEVIATION

**Description:**
The normal-queue RS job was converted to SRW, but its dev-queue siblings (`..._8N_400r_dev_2hr_safs.sbatch` and `..._8N_400r_dev_2hr_restart_safs.sbatch`) still defaulted to the aging config `spatial_friction_rate_state_safs_projected_stress_depthprofile.toml` + `param_a_vwvs11km.csv`. Since the user submits the **dev** job first to "check any issues," the smoke would have exercised the wrong physics (aging, not SRW), and the dev-restart job would have **resumed an SRW checkpoint with the aging config** — silently changing the friction law mid-stream.

**Trigger:** `sbatch …_dev_2hr_safs.sbatch` (smoke) or resuming via `…_dev_2hr_restart_safs.sbatch`.

**Actual behavior (pre-fix):** dev jobs ran aging-law RS, not SRW.

**Expected behavior:** dev jobs mirror the normal-queue SRW job exactly (config, CSVs, banner), differing only in queue/scale.

**Suggested fix (APPLIED):** same 5 edits as the normal job — `CONFIG_TOML` default → SRW config, CSV pre-flight → `_strongvw`, header + `EXPERIMENT`/`Config` banners → SRW. Verified: both dev jobs now reference the SRW config + `_strongvw` CSVs, `stray_aging=0`.

**Test case:**
```bash
test_R001_rs_dev_jobs_use_srw() {
  for f in jobs/safs/spatial_dyn_ratestate_depthprofile_normalcap_pureupwind_triq_8N_400r_dev_2hr*.sbatch; do
    grep -q 'srw_Dc010_nuc8km_500m.toml' "$f" || { echo "FAIL $f: not SRW"; exit 1; }
    grep -q 'param_a_500m_strongvw.csv'  "$f" || { echo "FAIL $f: not strongvw CSV"; exit 1; }
    grep -qE 'depthprofile\.toml|vwvs11km' "$f" && { echo "FAIL $f: stray aging ref"; exit 1; }
  done; echo PASS
}
```

---

### [R-002] MODERATE [POSSIBLE] [SRW physics parity across branches] — "same results" not verifiable locally; needs the dev smoke + verifier

**Category:** ASSUMPTION

**Description:**
The config + CSVs + mesh are faithful copies, but this repo's SRW *implementation* (`SlipLawSRWPsi` / the SRW friction iterator) is a **later lineage** than `seas-mfem-safs` — cf. this branch's commit `e0a99f0` "fix SRW friction-iterator factory abort on fault-less ranks (np>1)", which the `safs` branch does not have. Identical inputs do **not** guarantee bit-identical trajectories if the iterator changed. This cannot be checked locally (no local full-mesh runs; the binary isn't built here).

**Trigger:** Comparing this repo's SRW run against a prior `seas-mfem-safs` SRW run.

**Actual/Expected:** unknown until run; expected to match to the extent the SRW algorithm is unchanged (differences, if any, would be the np>1 fault-less-rank fix and any later SRW bugfix).

**Suggested fix (verification action, not a code edit):** run the **dev** RS job, confirm clean init (`--print-derived` gate passes) + bounded `V_max`, then run the on-fault verifier and compare station traces / nucleation timing to the source result:
```bash
SAFS_TFINAL=2s sbatch jobs/safs/spatial_dyn_ratestate_depthprofile_normalcap_pureupwind_triq_8N_400r_dev_2hr_safs.sbatch
# then: python3 safs/project_7.0_alternative/spatial/code/scripts/verify_spatial_dyn_smoke_safs.py \
#         --fault-vtkhdf <OUT>/fault.vtkhdf --toml-config <abs SRW toml> --log <log>
```

**Test case:** (system-level; run on Frontera) assert `V_max` finite and bounded through ≥ a few hundred steps, and `[derived]` RS gate does not abort.

---

### [R-003] LOW [POSSIBLE] [source SRW config] — source's plain `f_w`/`V_w` keys were likely silently ignored by a `_default`-reading parser; the port is more correct and physically equivalent

**Category:** ASSUMPTION

**Description:**
The source SRW config used top-level `f_w = 0.2` / `V_w = 0.1`. This repo's parser reads `f_w_default` / `V_w_default` and would **ignore** the un-suffixed keys, falling back to its built-in defaults (`f_w_default`=0.2, `V_w_default`=0.1 — parser lines 663/664). Those defaults coincide with the intended values, so the effective SRW physics is identical whether or not the rename is applied. The port makes the keys explicit (correct), and if the `safs`-branch parser had the same `_default` convention, the source itself was running on the (matching) defaults. No action needed; documented so a future value change (e.g. `f_w=0.3`) isn't silently dropped.

**Suggested fix:** none (already correct after the rename). Just be aware: always use the `_default` suffix for global SRW friction params in this repo.

---

### [R-004] LOW [jobs/safs/*_dev_2hr_safs.sbatch] — dev `tfinal=100s` will NOT finish in the 2 h wall, so the verifier + clean-exit path never run during a smoke

**Category:** QUALITY

**Description:**
The dev jobs inherit `SAFS_TFINAL=100s` (the normal-queue value). On the 2 h dev wall the run reaches only ~1.5–2 s of sim time and is then SIGTERM-killed → `RC!=0` → the post-run verifier block (gated on `RC==0`) is skipped. For a "submit and check any issues" smoke, a run that **finishes** is far more informative (clean shutdown, checkpoint flush, verifier).

**Suggested fix (recommend, not forced — preserves the established 100 s default):** submit the smoke with a short finishing `tfinal`:
```bash
SAFS_TFINAL=2s  sbatch jobs/safs/spatial_dyn_ratestate_depthprofile_normalcap_pureupwind_triq_8N_400r_dev_2hr_safs.sbatch
SAFS_TFINAL=2s  sbatch jobs/safs/spatial_dyn_slipweakening_normalcap_pureupwind_triq_8N_400r_dev_2hr_safs.sbatch
```
2 s clears the `--print-derived` gate, the gradual-overstress nucleation (`T_nuc=1 s`), and a post-nucleation advance, and exits cleanly so the verifier runs and the checkpoint/restart path is exercised.

**Test case:** (system-level) with `SAFS_TFINAL=2s`, assert the job reaches `driver_rc=0` and emits a `verify_summary_*.txt`.

---

### [R-005] LOW [jobs/safs/*.sbatch] — all SAFS jobs hard-require an `MFEM_USE_H5Z_ZFP=YES` build; without it they abort at pre-flight

**Category:** ASSUMPTION (pre-existing; flagged because the user is about to submit)

**Description:**
Every SAFS job greps `config/config.mk` for `MFEM_USE_H5Z_ZFP = YES` (and checks for `extern/h5z-zfp/install/plugin/libh5zzfp.so`) because they pass `--paraview-fault-zfp-tol 1.0e-12`. If `seas_spatial_dyn_driver` was built without H5Z_ZFP, the job exits 1 before running. This is intentional in this repo (HDF5+ZFP is the production fault-output path), but is a build prerequisite to confirm before the smoke.

**Suggested fix:** ensure the Frontera build has H5Z_ZFP (`bash build_frontera.sh` builds the plugin), **or** for a quick smoke on a non-ZFP build, override the fault output to VTU:
```bash
# (only if the build lacks H5Z_ZFP) — edit the smoke invocation to drop ZFP:
#   --paraview-fault-vtu   (instead of --paraview-fault-hdf5 --paraview-fault-zfp-tol …)
```

---

## Summary
- Critical issues: 0
- Moderate issues: 2 (R-001 incomplete dev transfer — **FIXED**; R-002 SRW cross-branch parity — verify via smoke)
- Low issues: 3 (R-003 source key-name latent no-op; R-004 dev tfinal for a finishing smoke; R-005 ZFP build dependency)
- Plan compliance: **FULL** after R-001 fix — both RS-SRW and LSW have consistent normal + dev + dev-restart jobs, all inputs present.
- Verdict: **PASS WITH FIXES** — R-001 applied this session; the transfer is internally consistent and parses. Remaining items are a verification action (R-002) and submit-time guidance (R-004/R-005), not code defects.

## Unreviewed Areas
- **Runtime numerical correctness** (does the SRW run reproduce the source's V_max(t)/nucleation): not checkable locally (no full-mesh runs); gated on the Frontera dev smoke (R-002).
- **HDF5 fault-output filename** (`fault.vtkhdf`) assumed by the verifier `--fault-vtkhdf` arg: pre-existing wiring, not re-derived from the driver source this round.
- **Volume/bulk ParaView output volume** on the 1.16 M-tet mesh under bare `--paraview`: pre-existing in all SAFS jobs; watch scratch usage during the smoke but not introduced by this transfer.
- **The LSW config's physics** (`law="slip_weakening"`, barrier nucleation): pre-existing in this repo, confirmed to parse + carry the floor; its numerical behaviour was not re-audited.

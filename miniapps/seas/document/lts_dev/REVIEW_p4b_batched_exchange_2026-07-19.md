# P4b — batched seam exchange: impl + review — 2026-07-19

Branch `safs-v4_0_0-alt-case1-mfem-speed` (local, unpushed).  Commits f351309
(I-batch), 06acc04 (D(k)-batch), 5580169 (test), 2deb2c4 (cleanup).

## What 4b (this increment) delivers
The COLLECTIVE-COUNT reduction of the D-7 3-buffer exchange: the per-cluster seam
exchanges, which were NUM_STATE per-component scalar `ExchangeFaceNbrData` calls in
4a, are now BATCHED (all-NUM_STATE) collectives via the vdim=NUM_STATE byNODES
`ghost_gf_full_state_`:
- seam time-integral `I` (Correct): NUM_STATE → **1** collective per correcting cluster.
- provider `D(k)` (Predict, c≥1): NUM_STATE·order → **order** collectives per cluster.

BYTE-IDENTICAL to 4a: `Ordering::byNODES` is component-major, so the component-major
`I`/`D(k)` copies straight into `GetData()`; the ghost is read via the layout-agnostic
`GetFaceNbrElementVDofs` vdof map (`nbr_all[nbr_vdofs[c*ndof2+i]]`, sign-decoded) —
the same per-(comp,dof) values as the 4a per-component `nbr_data`.  This mirrors the
already-validated `EvaluateBulkAtFaultQPsCanonical` batched exchange.  `ghost_dk_[c]`
is now `order` full-state Vectors (was `order·NUM_STATE`).  `BuildTickTable`
`n_collectives`: I term `num_state·|correct|` → `|correct|`; D(k) term
`num_state·ader_order` → `ader_order`.

## Byte-vs-4a verification (GREEN — all Δ identical to 4a)
| Gate | 4b result (== 4a) |
|---|---|
| diff-0 seam np=2 | Δenergy=4.441e-16, Δcen=0.0, Δmax=0.0 |
| diff-1 Nc=3 np=3 | 1.614e-16 / 1.780e-16 / 0.0 |
| fault interleave np=2 / np=10 | 1.666e-27 / 4.443e-27, PASS |
| SCEC stations (fault) | 9/9 byte-identical |
| per-sync collective counter | matches the REDUCED n_collectives |
| regression | predictor 31 · reorder 14 · scheduler 31 · layout 102 · clustering 4043 · partition 19 · tpv102-ADER 4 · seam np2/4 60/120 |

`test_lts_layout` T4 n_collectives updated (9,18,9,27 → 1,2,1,3) — the batched count
is unit-covered, not just MPI-gated.

## Adversarial review — PASS-WITH-FIXES (byte-identity confirmed + unit-covered)
Verified: the byNODES/vdof pack+read is character-for-character the validated
reference; `ndof2==ndof_per_el_` guarded; the D(k) cache is a genuine deep copy
(no aliasing when the same `q_gf_full` is reused for I then D(k)); `nbr_all` is
aliased once and never re-exchanged inside the face loop; the counter matches.
- **[LOW] FIXED** — the Stage-1 `lts_meta.num_state == NUM_STATE` pin is now
  vacuous (4b dropped both num_state factors from `n_collectives`); removed.
- **[LOW/perf] noted (not done)** — the batched pack zeros the whole `q_gf_full` +
  scans all `ne_` per exchange.  NOT a regression (4a did the same total work; 4b
  actually scans `ne_` 9× fewer times).  A dirty-list pack is a further
  constant-factor win, byte-identical — deferred.

## REMAINING 4b scope — the flux-PREMULTIPLIED per-face block (D-7 payload iii)
This increment is the BATCHED-buffer optimization (collective-count).  The deeper
EDGE innovation — sending the FLUX-PROJECTED payload (premultiplied by the
receiver's test functions) instead of the raw `I`/`D(k)` coefficients — is NOT
implemented.  It would: (a) cut the diff-1 `D(k)` bandwidth ~order× (send the
per-face flux block `n_seam_faces × NUM_STATE × nbf_face`, not the `order`-level
Taylor stack), and (b) halve the seam flux evaluation (one side computes + sends;
the other adds).  It is byte-comparable against this batched 4b and its benefit is
Expanse-scale.  RECOMMENDATION: measure the batched-4b comm on Expanse (Phase 5)
before the substantial premultiplied-block investment (measure-first; the plan
already measures Phase-5 on 4b, and the batched exchange may already suffice).

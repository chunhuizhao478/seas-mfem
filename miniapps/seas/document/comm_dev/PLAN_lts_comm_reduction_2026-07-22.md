# Implementation Plan: LTS communication reduction + overlap

**Date:** 2026-07-22 (v1, scoped from the evidence-round pre-scope) · **Branch:** `safs-v4_0_0-alt-case1-mfem-speed`
**Motivation:** `document/lts_dev/RESULTS_p5_tpv104_200m_2026-07-21.md` — exposed `MPI_Waitall` = 53 % of
MFEM-LTS wall (2341 of 4356 s/sim-s on TPV104-200m); 585 M P2P messages, 17.3 M waits, ~125
ghost-exchange rounds per sync at Nc=6.
**Companion:** `document/kernel_dev/PLAN_ader_kernel_efficiency_2026-07-21.md` (rev 3) — the composed
end state (both programs) is what approaches SeisSol; see its Summary composition table.

## Summary — read this first

**The problem.** Under LTS, more than half the wall-clock is spent *waiting for ghost exchanges
to complete*. The tick loop performs ~125 full-halo exchange rounds per sync interval — one per
cluster-correction plus one forecast round each — and every one blocks. The same solver under
global time stepping spends 0.00002 % there: this is entirely a structure-of-the-LTS-loop cost.

**The fix.** Three moves, in order of certainty: (1) **merge** — same-tick corrections touch
disjoint elements, so their exchanges superpose into one round per tick (125 → 32, the
theoretical minimum), *bitwise identical by construction*; (2) **overlap** — post the tick's
exchange early and complete it late, hiding the wire time behind the tick's ghost-independent
compute (friction, volume, interior faces), via a small custom exchanger built on MFEM's public
communication tables; (3) **shrink** — ship only seam-adjacent elements of active clusters
instead of the full halo (last, only if measurement says bytes matter).

**Expected outcome.** Exposed wait drops from 2341 to ~150–400 s/sim-s (53 % → 7–16 % of wall):
**~1.8–2.0× on the LTS leg by itself**, and the enabling half of the composed
"within ~2.4–3.9× of SeisSol" end state (the composed state is comm-bound, so this program is
not optional).

**Main risk.** An unknown fraction of the measured wait is **late-sender skew** (fault-rank
friction imbalance) rather than communication — skew is untouched by merging or overlap. Phase 0
decomposes wire-vs-skew before anything is built; if skew > ~40 %, the fault-weighted partition
(existing METIS ncon=2 v3 plan) is promoted to co-requisite and the payoff here is restated.

**What this does NOT do.** No kernel work (companion plan). No change to exchanged *values* or
read indices in Phases 1–2 — bitwise identity is the acceptance bar there. No GTS-path changes
except the explicitly pinned Phase-4 retrofits. No MFEM-core edits (the ldof tables are public).

## Glossary

| Label | Meaning |
|---|---|
| round | one `ExchangeFaceNbrData`-equivalent: pack → Isend/Irecv to ~34 neighbors → Waitall×2 |
| tick | one dt_base step of the LTS sync interval (32 ticks/sync at Nc=6) |
| I-exchange | per-correct exchange of the cluster's time-integral in `ghost_gf_full_state_` (vdim=9 byNODES) |
| forecast exchange | per-correct (c<Nc-1) exchange of the premultiplied seam forecast (P4b) |
| P-007 audit | driver assert: measured ghost-exchange count == `BuildTickTable.n_collectives`; a mismatch aborts (or hangs at np≥10) — must be updated in lockstep with any exchange change |
| R-004 hazard | byNODES ghost-unpack trap: every new vdim/packing has a different FaceNbrData interleaving; unpacks must go through `GetFaceNbrElementVDofs` maps, never slab formulas |
| skew | Waitall time caused by a late sender (load imbalance), not by wire/bandwidth |

## Exchange inventory (measured/verified call sites, Nc=6)

| site (`wave_operator.inl`) | trigger | per sync | payload | mergeable |
|---|---|---:|---|---|
| seam I-exchange (:2727, from correct step 2c :2482) | every Correct(c) | 63 | full halo vdim=9 (only cluster-c blocks nonzero) | → per-tick: 32 |
| forecast pre-pass (:2703, deep-copy :2706) | Correct(c), c<Nc-1 | 62 | full halo vdim=9 (only seam-coarse blocks nonzero) | → per-tick: 32; joinable with I via 2nd buffer/vdim=18 (never superposed — one ghost elem can carry both) |
| `PrepareSeamCoarseForecast` (:3000) | Predict(c≥1) | 0 | retain-only since P4b | already eliminated |
| substep fault exchange (:3712, driver :712) | GTS path only | 0 (LTS: interior-only under D-2) | vdim=9 | GTS retrofit (Phase 4) |
| shared-corrector 9-sequential (:6513-24) | GTS path only | 0 | 9 scalar halos | GTS retrofit: 9→1 (Phase 4, TPV byte pins) |
| driver V_max Allreduce + P-007 audit | per sync | 1–2 | scalar | keep (safety net) |

Theoretical minimum: **32/sync** (cluster 0 predicts *and* corrects every tick; its I is produced
in the same tick it is consumed — no cross-tick pipelining possible).

## Phases

**Phase 0 — Waitall decomposition + accounting baseline.** Timestamp post-vs-complete around
both exchange sites (Caliper or manual) → split the 2341 s/sim-s into wire vs skew, per
cluster-tick class; record bytes/round; verify the 125/sync count. *Gate:* decomposition sums to
the Phase-5 totals; no code-path change (byte-identity trivial). **This phase decides the skew
question and can ride in the same Expanse job as the kernel plan's B0.**

**Phase 1 — Per-tick merge (125 → 64).** Hoist I-pack and forecast-pack out of `Correct(c)`
into a per-tick pre-correct hook packing ALL correcting clusters (disjoint element sets) into
`ghost_gf_full_state_` + a second vdim=9 GF; the seam corrector consumes the tick buffers.
Validity: all predicts precede all corrects within a tick (`RunSyncInterval`); mode-0/1/2 reads
filter by `ghost_cluster_id_` so every correct reads exactly the bytes it reads today. *Gate:*
**bitwise-identical Q at every sync** vs HEAD on the LTS fixtures at np∈{1,2,4,10}; P-007 audit
updated in the same commit; `fc_all` deep copy retired only after the double-buffer lands
(R4F-003 lifetime rule).

**Phase 2 — Single-round tick exchange (64 → 32).** Post both buffers' Isend/Irecv together,
one wait phase (or one vdim=18 GF, I and forecast in separate component planes). *Gate:*
bitwise identity + a byNODES-unpack parity test cloned for the new packing (the R-004 hazard IS
this phase's risk); `n_collectives` = 1/tick.

**Phase 3 — Split-post overlap (`NbrExchangerSplit`).** New ~80-line file under `dynamic/`
replicating `pgridfunc.cpp:209-268` as Begin()/End() with persistent requests
(`MPI_Send_init/Recv_init` + `Startall`; safe because one round/tick ⇒ send buffer quiescent)
over the **public** `send_face_nbr_ldof`/`face_nbr_ldof` tables — no MFEM-core edit. Begin after
the tick's predict-phase bulk CK (earliest I_[0] exists); End just before the tick's first
seam-face sweep; the hiding window = all predict friction + correct volume/interior/fault work.
*Gate:* values and unpack indices unchanged ⇒ **bitwise identity mandatory**; np=10 no-hang;
`debug` device backend clean (Read()/Write() UseDevice gotcha); P-007 counts Begin() calls;
GPU-aware-MPI branch mirrored.

**Phase 4 — Sparse payloads + GTS retrofits + acceptance.** Per-(cluster, neighbor) ldof
sublists (bytes floor) ONLY if Phase 0 shows bandwidth (not skew) binding; GTS-path wins
(9-sequential → 1; pre-post the O substep exchanges) behind the TPV byte-exact pins. *Gate:*
LTS fixtures bitwise; TPV102/104/205 byte-exact for any GTS change; **acceptance: exposed
Waitall ≤ ~15–20 % of wall on the Phase-5 config** — else escalate the skew findings to the
fault-weighted-partition plan.

## Payoff arithmetic (from the pre-scope; assumptions stated there)

| lever | exposed Waitall (s/sim-s) | LTS wall | speedup |
|---|---:|---:|---:|
| today | 2341 | 4356 | 1.0× |
| merge → 64 | ~1200 | ~3280 | 1.35× |
| merge → 32 | ~600 | ~2680 | 1.65× |
| overlap alone | ~950–1400 | ~3030–3480 | 1.3–1.45× |
| **merge + overlap** | **~150–400** | **~2230–2480** | **~1.8–2.0×** |

Message shape: ~34 neighbors/rank, per-message O(10–100 KB) — not purely latency-bound; the 53 %
exposure is round-count × soft-synchronization skew, which is exactly what merging + overlap
attack. The skew caveat (Summary) bounds the downside.

## Risks

- **Skew misattribution** (Summary) — Phase-0 decomposition is the go/no-go.
- **R-004 unpack hazard** on every new vdim/packing — `GetFaceNbrElementVDofs`-only + parity
  test per layout.
- **Matched-collective lockstep** — `BuildTickTable.n_collectives`, the stepper, and the wave
  counter change in one commit or np≥10 hangs (R-1600 class); the P-007 audit stays, never
  relaxed.
- **Buffer lifetime** — persistent requests require the one-round-per-tick structure (Phase 3
  depends on Phase 2); `fc_all` copy removal requires the double-buffer.
- **Device path** — mirror the GPU-aware-MPI branch and the UseDevice side effects.
- Local fixtures cap at np≤10; production message shape is cluster-only; ask before submissions.

# TPV31 spatial — archived jobs

The active TPV31 spatial test matrix is the **four unified-bimaterial jobs** (one per
scheme × order), all wired with the full bimaterial stack:

| Active job | scheme | mesh | bimaterial fixes |
|---|---|---|---|
| `mixedflux_rk45/tpv31_p1_rk45_mixedflux_50m_normal_contrastguard.sbatch` | mixed flux p1 | original `tpv31_50m.msh` | volume (contrast guard) + fault (matrix) + cross-rank |
| `mixedflux_rk45/tpv31_p2_rk45_mixedflux_50m_normal_contrastguard.sbatch` | mixed flux p2 | original `tpv31_50m.msh` | volume + fault + cross-rank |
| `upwind_ader/tpv31_p1_pureupwind_50m_hybrid_normal.sbatch` | upwind p1 | symmetric `tpv31_50m_hybrid.msh` | fault (matrix) + cross-rank (no central flux ⇒ no volume guard) |
| `upwind_ader/tpv31_p2_pureupwind_50m_hybrid_normal.sbatch` | upwind p2 | symmetric `tpv31_50m_hybrid.msh` | fault + cross-rank |

Archived here (kept for reference, not part of the active matrix):

- `tpv31_p1_rk45_mixedflux_50m_normal.sbatch`, `tpv31_p2_..._normal.sbatch` — the **guard-OFF**
  mixed-flux arm (no `--mixed-flux-contrast-tol`); missing the bimaterial **volume** fix, so it
  exhibits the σ_n collapse at the 2400/5000/10000 m interfaces.  The A/B control for the guard.
- `tpv31_p1_pureupwind_50m_rvsmooth_split_normal.sbatch`, `tpv31_p2_..._rvsmooth_split_...` — the
  **rvsmooth-split** mesh variant of the upwind arm (volume-ratio regularization, not the y-mirror
  symmetric mesh).  Alternative near-fault-asymmetry control.

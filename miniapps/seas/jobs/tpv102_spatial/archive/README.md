# TPV102 spatial — archived jobs

**Note:** TPV102 is a HOMOGENEOUS benchmark (`[material_constant_fallback]`,
`interior_flux = "scalar"`).  It has no material contrast, so the bimaterial fixes
(volume contrast guard, per-side fault Riemann, cross-rank peer read) are inherently
N/A here — the cross-rank exchange is byte-exact for constant material.  TPV102 is kept
as the homogeneous comparison case, with the same 4-job layout as TPV6/TPV31:

| Active job | scheme | mesh |
|---|---|---|
| `mixedflux_rk45/tpv102_p1_rk45_mixedflux_10N_500r_normal.sbatch` | mixed flux p1 | original `tpv102_200m.msh` |
| `mixedflux_rk45/tpv102_p2_rk45_mixedflux_12N_600r_normal.sbatch` | mixed flux p2 | original `tpv102_200m.msh` |
| `upwind_ader/tpv102_p1_pureupwind_200m_hybrid_normal.sbatch` | upwind p1 | symmetric `tpv102_200m_hybrid.msh` |
| `upwind_ader/tpv102_p2_pureupwind_200m_hybrid_normal.sbatch` | upwind p2 | symmetric `tpv102_200m_hybrid.msh` |

Archived here (kept for reference, not part of the active matrix):

- `tpv102_p1_rk45_mixedflux_200m_xsym_normal.sbatch` — strike-symmetric (`xsym`) mesh variant of
  the p1 mixed-flux arm (dip-slip-drift A/B); the active matrix uses the original mesh.
- `tpv102_p1_rk45_mixedflux_8N_400r_dev.sbatch` — flex-queue dev/smoke run.
- `tpv102_p1_rk45_pureupwind_200m_flex.sbatch` — flex-queue pure-upwind dev run.
- `tpv102_p3_rk45_mixedflux_16N_800r_normal.sbatch` — p3 (order 3); the matrix only spans p1/p2.
- `tpv102_p1_pureupwind_200m_rvsmooth_normal.sbatch`, `tpv102_p2_..._rvsmooth_...` — rvsmooth mesh
  variant of the upwind arm (not the y-mirror symmetric/hybrid mesh).

The native-driver dev/bisect/init jobs live separately under `jobs/tpv102/` (untouched).

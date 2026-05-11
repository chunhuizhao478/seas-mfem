# Sensitivity sweep — CGAL corefine — SAFS 2000 m

Generated: 2026-05-08 16:54:07

## Per-config summary

| config | mesh_edge | min_edge | spacing | n_F | conformal | edge_min | edge_med | edge_max | q_min | q_med | elapsed | score |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| baseline | 1500 | 100 | 750 | 6098 | 7/7 | 0.0 | 1557.1 | 3341.4 | 0.4168 | 0.9161 | 0.4s | 0.872 |
| coarser | 2000 | 100 | 1000 | 4223 | 7/7 | 0.0 | 494.9 | 4358.1 | 0.4188 | 0.8943 | 0.3s | 0.658 |
| finer | 1000 | 100 | 500 | 11249 | 7/7 | 0.0 | 1187.9 | 2101.2 | 0.2609 | 0.9317 | 0.6s | 0.796 |
| tighter_floor | 1500 | 200 | 750 | 6049 | 7/7 | 0.0 | 1555.7 | 3341.4 | 0.3423 | 0.9159 | 0.4s | 0.857 |
| fine_floor | 1500 | 50 | 750 | 6098 | 7/7 | 0.0 | 1557.1 | 3341.4 | 0.4168 | 0.9161 | 0.4s | 0.872 |

**Best config:** `baseline` (`mesh_edge_size=1500m`, `min_edge=100m`, `polyline_spacing=750m`)

## Per-fault detail (best config: baseline)

| fault | V | F | edge_min | edge_med | edge_max | q_min | q_med | n_edges_below_100m |
|---|---|---|---|---|---|---|---|---|
| SAFS-SAFZ-COAV-Mission_Creek_fault_strand-CFM4_2000m | 868 | 1678 | 0.1 | 53.0 | 3102.1 | 0.5277 | 0.8639 | 1432 |
| SAFS-SAFZ-MJVS-San_Andreas_fault-CFM6_2000m | 396 | 658 | 0.1 | 1552.7 | 2924.9 | 0.5213 | 0.9186 | 291 |
| SAFS-SAFZ-MULT-Southern_San_Andreas_fault_and_Banning-C | 487 | 818 | 0.3 | 1561.6 | 2980.2 | 0.4168 | 0.9149 | 351 |
| SAFS-SAFZ-SBMT-Garnet_Hill_fault-CFM6_2000m | 507 | 877 | 0.3 | 1380.6 | 2986.2 | 0.5916 | 0.9104 | 378 |
| SAFS-SAFZ-SBMT-Mission_Creek_fault_strand-CFM4_2000m | 231 | 395 | 1038.1 | 1859.5 | 2950.1 | 0.5781 | 0.9172 | 0 |
| SAFS-SAFZ-SBMT-San_Andreas_fault-CFM6_2000m | 954 | 1672 | 0.0 | 1739.9 | 3341.4 | 0.5257 | 0.9247 | 450 |

## Score formula

```
score = 0.5 * (n_conformal_pairs / n_pairs)
      + 0.3 * max(0, 1 - |edge_med - mesh_edge_size| / mesh_edge_size)
      + 0.2 * q_overall_min
```

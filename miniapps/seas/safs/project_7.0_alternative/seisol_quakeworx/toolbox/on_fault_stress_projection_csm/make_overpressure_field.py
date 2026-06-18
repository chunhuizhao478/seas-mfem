#!/usr/bin/env python3
"""make_overpressure_field.py — build a recipe's fault-local overpressure .npz
(PLAN_gate_removal_csm_2026-06-17.md, Phase 3).

Projects the CSM C1 field at the recipe's closure k on the fault mesh (V3 PUML
by default), then calls fault_local_overpressure.build_overpressure_field to
make the 2D (s,z) graded DeltaPp field that lifts the gate cores to the S<1.7
floor without going supercritical, and saves it as an .npz consumed by
project_csm_stress_to_vtu.py / csm_stress_to_asagi.py via --overpressure-field.

Recipe B: --k-ratio 1.8 --mu-s 0.318 (mu_s = max mu_app(k=1.8) + 0.02).
Recipe C: --k-ratio 1.5 --mu-s 0.225.
"""
import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from gate_removal_sweep import (  # noqa: E402
    Field, GATES_KM, SEIS_ZKM, STRIKE_AZ, HYPO, DEFAULT_MESH, MU_D, S_TARGET,
    mu_s_for_k)
from fault_local_overpressure import (  # noqa: E402
    build_overpressure_field, save_overpressure_field)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--mesh", default=DEFAULT_MESH)
    ap.add_argument("--k-ratio", type=float, required=True)
    ap.add_argument("--mu-s", type=float, default=None,
                    help="static friction; default = max mu_app(k) + 0.02")
    ap.add_argument("--gates", default="SanGorgonio,g181,g221",
                    help="comma list of gate names from GATES_KM")
    ap.add_argument("--lambda-max", type=float, default=0.9)
    ap.add_argument("--target-mu", type=float, default=None,
                    help="aim point for the gate cores (>= floor, < mu_s); "
                    "default = floor. Aim above the floor to leave margin "
                    "against volume-grid blur after baking.")
    ap.add_argument("--recipe", default="B")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    field = Field(args.mesh)
    mu, tau, sn, active = field.project(args.k_ratio)
    mu_s = args.mu_s if args.mu_s is not None else mu_s_for_k(field, args.k_ratio)[0]
    bands = []
    for g in args.gates.split(","):
        if g not in GATES_KM:
            sys.exit(f"ERROR: unknown gate {g!r}; choices {list(GATES_KM)}")
        bands.append(GATES_KM[g])
    depth_m = np.maximum(-field.cent[:, 2], 0.0)
    op = build_overpressure_field(
        field.s_km, depth_m, sn, tau, field.Sv_total, field.Pp_h, field.seis,
        bands, mu_s, STRIKE_AZ, HYPO[:2], mu_d=MU_D, s_target=S_TARGET,
        lambda_max=args.lambda_max, seis_zkm=SEIS_ZKM,
        core_target_mu=args.target_mu,
        meta={"recipe": args.recipe, "k_ratio": float(args.k_ratio),
              "mesh": os.path.basename(args.mesh), "gates": args.gates})
    out = args.out or os.path.join(
        HERE, f"overpressure_{args.recipe}_2026-06-17.npz")
    save_overpressure_field(out, op)
    m = op["meta"]
    print(f"recipe {args.recipe}: k={args.k_ratio} mu_s={mu_s:.3f} "
          f"floor={m['floor']:.3f} gates={args.gates}")
    print(f"  grid {op['DPp'].shape} (s x depth); DeltaPp max {m['DPp_max_MPa']:.1f} "
          f"MPa, mean(nonzero) {m['DPp_mean_nonzero_MPa']:.1f} MPa; "
          f"under-raised heterogeneous cells {m['n_hetero_cells_under_raised']}")
    print(f"  wrote {out}")


if __name__ == "__main__":
    main()

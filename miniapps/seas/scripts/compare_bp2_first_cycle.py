#!/usr/bin/env python3
"""
BP2 First Earthquake Cycle: Comparison and Plotting Script

Compares simulation output against Erickson's finite-difference reference
solution (50m resolution) for the SCEC SEAS BP2-QD benchmark.

Usage:
    python3 compare_bp2_first_cycle.py \
        --sim-dir . --ref-dir bp2/benchmark_data_200m \
        --output-dir bp2_verification_plots
"""

import argparse
import os
import sys
import numpy as np

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("WARNING: matplotlib not available. Skipping plots.")


# Standard 12 SCEC probe depths [m]
PROBE_DEPTHS_M = [
    0.0, -2400.0, -4800.0, -7200.0, -9600.0,
    -12000.0, -14400.0, -16800.0, -19200.0,
    -24000.0, -28800.0, -36000.0
]

# Highlight depths for summary figures
HIGHLIGHT_DEPTHS_KM = [0.0, 7.2, 12.0, 24.0, 36.0]

SECONDS_PER_YEAR = 365.25 * 24 * 3600


def depth_to_km(depth_m):
    return abs(depth_m) / 1000.0


def make_sim_filename(prefix, depth_m):
    """Generate simulation output filename matching BenchmarkOutput::MakeFilename."""
    depth_km = depth_to_km(depth_m)
    if abs(depth_km - round(depth_km)) < 1e-6:
        depth_str = f"{int(round(depth_km))}"
    else:
        depth_str = f"{depth_km:.1f}"
    return f"{prefix}_z{depth_str}km.txt"


def make_ref_filename(depth_m):
    """Generate Erickson reference filename."""
    depth_km = depth_to_km(depth_m)
    if abs(depth_km - round(depth_km)) < 1e-6:
        depth_str = f"{int(round(depth_km))}"
    else:
        depth_str = f"{depth_km:.1f}"
    return f"bp2-qd-erickson-z{depth_str}km-res.txt"


def load_scec_file(filepath, t_max=None):
    """Load SCEC-format time series file (5 columns, # comment lines)."""
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if 'slip' in line and 'time' not in line.split()[0]:
                continue
            # Skip header line with column names
            if line.startswith('t ') or line.startswith('t\t'):
                continue
            try:
                vals = [float(x) for x in line.split()]
                if len(vals) >= 5:
                    if t_max is not None and vals[0] > t_max:
                        break
                    data.append(vals[:5])
            except ValueError:
                continue

    if not data:
        return None

    arr = np.array(data)
    return {
        'time': arr[:, 0],
        'slip': arr[:, 1],
        'log10_slip_rate': arr[:, 2],
        'shear_stress': arr[:, 3],
        'log10_state': arr[:, 4],
    }


def relative_l2_error(sim, ref):
    """Compute relative L2 error: ||sim - ref||_2 / ||ref||_2"""
    num = np.sqrt(np.sum((sim - ref)**2))
    den = np.sqrt(np.sum(ref**2))
    if den < 1e-30:
        return 0.0 if num < 1e-30 else 1e30
    return num / den


def interpolate_onto(x_ref, y_ref, x_target):
    """Linear interpolation of reference onto target time grid."""
    return np.interp(x_target, x_ref, y_ref)


def plot_probe_comparison(sim_data, ref_data, depth_km, errors, output_path):
    """Generate 4-panel comparison figure for a single probe depth."""
    if not HAS_MATPLOTLIB:
        return

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'BP2-QD Verification: z = {depth_km} km', fontsize=14)

    time_sim_yr = sim_data['time'] / SECONDS_PER_YEAR
    time_ref_yr = ref_data['time'] / SECONDS_PER_YEAR

    panels = [
        ('slip', 'Slip [m]', False),
        ('log10_slip_rate', 'log$_{10}$(Slip Rate) [m/s]', False),
        ('shear_stress', 'Shear Stress [MPa]', False),
        ('log10_state', 'log$_{10}$(State) [s]', False),
    ]
    error_keys = ['slip', 'slip_rate', 'stress', 'state']

    for idx, (key, ylabel, _) in enumerate(panels):
        ax = axes[idx // 2][idx % 2]
        ax.plot(time_ref_yr, ref_data[key], 'k--', linewidth=1.0,
                label='Reference (Erickson 50m)', alpha=0.8)
        ax.plot(time_sim_yr, sim_data[key], 'b-', linewidth=1.0,
                label='MFEM (200m)')
        ax.set_xlabel('Time [yr]')
        ax.set_ylabel(ylabel)
        err = errors.get(error_keys[idx], -1)
        if err >= 0:
            ax.set_title(f'{ylabel}  (L2 err: {err:.4f})')
        else:
            ax.set_title(ylabel)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_error_summary(all_errors, output_path):
    """Generate summary bar chart of L2 errors across all depths."""
    if not HAS_MATPLOTLIB:
        return

    depths = [e['depth_km'] for e in all_errors if e.get('has_data')]
    slip_errs = [e['slip'] for e in all_errors if e.get('has_data')]
    rate_errs = [e['slip_rate'] for e in all_errors if e.get('has_data')]
    stress_errs = [e['stress'] for e in all_errors if e.get('has_data')]
    state_errs = [e['state'] for e in all_errors if e.get('has_data')]

    if not depths:
        return

    x = np.arange(len(depths))
    width = 0.2

    fig, ax = plt.subplots(figsize=(14, 6))
    ax.bar(x - 1.5*width, slip_errs, width, label='Slip', color='#2196F3')
    ax.bar(x - 0.5*width, rate_errs, width, label='Slip Rate', color='#FF9800')
    ax.bar(x + 0.5*width, stress_errs, width, label='Stress', color='#4CAF50')
    ax.bar(x + 1.5*width, state_errs, width, label='State', color='#9C27B0')

    ax.set_xlabel('Depth [km]')
    ax.set_ylabel('Relative L2 Error')
    ax.set_title('BP2-QD First Cycle: Relative L2 Errors vs Erickson Reference')
    ax.set_xticks(x)
    ax.set_xticklabels([f'{d:.1f}' for d in depths], rotation=45)
    ax.legend()
    ax.grid(True, axis='y', alpha=0.3)
    ax.set_yscale('log')

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description='Compare BP2 simulation output against Erickson reference')
    parser.add_argument('--sim-dir', default='.', help='Simulation output directory')
    parser.add_argument('--ref-dir', default='bp2/benchmark_data_200m',
                        help='Reference data directory')
    parser.add_argument('--output-dir', default='bp2_verification_plots',
                        help='Plot output directory')
    parser.add_argument('--prefix', default='bp2_verify',
                        help='Simulation output file prefix')
    parser.add_argument('--t-max-yr', type=float, default=250.0,
                        help='Maximum time in years to compare')
    args = parser.parse_args()

    t_max_s = args.t_max_yr * SECONDS_PER_YEAR

    os.makedirs(args.output_dir, exist_ok=True)

    print(f"BP2 First Cycle Comparison")
    print(f"  Simulation dir: {args.sim_dir}")
    print(f"  Reference dir:  {args.ref_dir}")
    print(f"  Output dir:     {args.output_dir}")
    print(f"  t_max:          {args.t_max_yr} yr")
    print()

    all_errors = []

    # Header for summary table
    print(f"{'Depth(km)':>10s} {'Slip L2':>12s} {'Rate L2':>12s} "
          f"{'Stress L2':>12s} {'State L2':>12s} "
          f"{'Tnuc_sim(yr)':>14s} {'Tnuc_ref(yr)':>14s}")
    print('-' * 88)

    for depth_m in PROBE_DEPTHS_M:
        depth_km = depth_to_km(depth_m)
        sim_file = os.path.join(args.sim_dir,
                                make_sim_filename(args.prefix, depth_m))
        ref_file = os.path.join(args.ref_dir, make_ref_filename(depth_m))

        error_entry = {'depth_km': depth_km, 'has_data': False}

        # Load data
        sim_data = load_scec_file(sim_file, t_max_s) if os.path.exists(sim_file) else None
        ref_data = load_scec_file(ref_file, t_max_s) if os.path.exists(ref_file) else None

        if sim_data is None:
            print(f"{depth_km:10.1f} {'MISSING':>12s} (sim: {sim_file})")
            all_errors.append(error_entry)
            continue
        if ref_data is None:
            print(f"{depth_km:10.1f} {'MISSING':>12s} (ref: {ref_file})")
            all_errors.append(error_entry)
            continue

        # Interpolate reference onto simulation time grid
        ref_interp = {}
        for key in ['slip', 'log10_slip_rate', 'shear_stress', 'log10_state']:
            ref_interp[key] = interpolate_onto(
                ref_data['time'], ref_data[key], sim_data['time'])

        # Compute errors
        errors = {
            'slip': relative_l2_error(sim_data['slip'], ref_interp['slip']),
            'slip_rate': relative_l2_error(sim_data['log10_slip_rate'],
                                            ref_interp['log10_slip_rate']),
            'stress': relative_l2_error(sim_data['shear_stress'],
                                         ref_interp['shear_stress']),
            'state': relative_l2_error(sim_data['log10_state'],
                                        ref_interp['log10_state']),
        }

        # Nucleation time (first time log10(V) > -3)
        tnuc_sim_yr = -1.0
        tnuc_ref_yr = -1.0
        mask_sim = sim_data['log10_slip_rate'] > -3.0
        if np.any(mask_sim):
            idx = np.argmax(mask_sim)
            tnuc_sim_yr = sim_data['time'][idx] / SECONDS_PER_YEAR
        mask_ref = ref_data['log10_slip_rate'] > -3.0
        if np.any(mask_ref):
            idx = np.argmax(mask_ref)
            tnuc_ref_yr = ref_data['time'][idx] / SECONDS_PER_YEAR

        tnuc_sim_str = f"{tnuc_sim_yr:.1f}" if tnuc_sim_yr > 0 else "N/A"
        tnuc_ref_str = f"{tnuc_ref_yr:.1f}" if tnuc_ref_yr > 0 else "N/A"

        print(f"{depth_km:10.1f} {errors['slip']:12.4e} "
              f"{errors['slip_rate']:12.4e} {errors['stress']:12.4e} "
              f"{errors['state']:12.4e} {tnuc_sim_str:>14s} {tnuc_ref_str:>14s}")

        error_entry.update(errors)
        error_entry['has_data'] = True
        error_entry['tnuc_sim_yr'] = tnuc_sim_yr
        error_entry['tnuc_ref_yr'] = tnuc_ref_yr
        all_errors.append(error_entry)

        # Plot individual probe comparison
        if HAS_MATPLOTLIB:
            plot_path = os.path.join(
                args.output_dir, f'bp2_verify_z{depth_km:.1f}km.png')
            plot_probe_comparison(sim_data, ref_data, depth_km, errors,
                                  plot_path)

    # Summary plot
    if HAS_MATPLOTLIB:
        summary_path = os.path.join(args.output_dir,
                                     'bp2_verify_error_summary.png')
        plot_error_summary(all_errors, summary_path)
        print(f"\nPlots saved to: {args.output_dir}/")

    # Print overall assessment
    probes_with_data = [e for e in all_errors if e.get('has_data')]
    probes_with_eq = [e for e in probes_with_data
                      if e.get('tnuc_sim_yr', -1) > 0]

    print(f"\n{'=' * 40}")
    print(f"Probes with data:       {len(probes_with_data)} / {len(PROBE_DEPTHS_M)}")
    print(f"Probes detecting EQ:    {len(probes_with_eq)}")

    if probes_with_eq:
        for e in probes_with_eq:
            if e.get('tnuc_ref_yr', -1) > 0:
                rel_err = abs(e['tnuc_sim_yr'] - e['tnuc_ref_yr']) / e['tnuc_ref_yr']
                status = "OK" if rel_err < 0.10 else "WARNING"
                print(f"  z={e['depth_km']:.1f}km: Tnuc={e['tnuc_sim_yr']:.1f}yr "
                      f"(ref={e['tnuc_ref_yr']:.1f}yr, err={rel_err*100:.1f}%) "
                      f"[{status}]")

    has_nan = any(np.isnan(e.get('slip', 0)) or np.isinf(e.get('slip', 0))
                  for e in probes_with_data)
    print(f"NaN/Inf detected:       {'YES' if has_nan else 'NO'}")

    passed = len(probes_with_data) > 0 and len(probes_with_eq) > 0 and not has_nan
    print(f"Overall:                {'PASS' if passed else 'FAIL'}")
    print(f"{'=' * 40}")

    return 0 if passed else 1


if __name__ == '__main__':
    sys.exit(main())

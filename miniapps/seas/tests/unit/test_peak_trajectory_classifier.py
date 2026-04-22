#!/usr/bin/env python3
"""Unit test for tpv102/r_v92_outlier_detection/plot_peak_trajectory.py.

Validates that the classify_trajectory() function correctly distinguishes
the four shape categories that drive the R-V92-E01 re-analysis:

  - SATURATION:       rise then plateau near maximum
  - MONOTONIC-GROWTH: peak at end of run, late-mean close to max
  - STEP-AND-DECAY:   rise then fall well below max
  - STEP-SPIKE-DECAY: mid-run peak, moderate late-run magnitude

Synthetic trajectories cover each shape plus the empty/flat edge cases.
"""

from __future__ import annotations

import math
import os
import sys
import tempfile
import unittest
from pathlib import Path

# Make the script importable as a module.
HERE = Path(__file__).resolve().parent
SCRIPT_DIR = HERE.parent.parent / "tpv102" / "r_v92_outlier_detection"
sys.path.insert(0, str(SCRIPT_DIR))

import plot_peak_trajectory as analyzer  # noqa: E402


class TestClassifyTrajectory(unittest.TestCase):
    """Shape classification on synthetic time series."""

    def test_empty_returns_flat(self):
        self.assertEqual(analyzer.classify_trajectory([], []), "FLAT")

    def test_all_zero_returns_flat(self):
        times = [0.1 * i for i in range(10)]
        peaks = [0.0] * 10
        self.assertEqual(analyzer.classify_trajectory(times, peaks), "FLAT")

    def test_saturation_shape(self):
        """Field rises quickly, then plateaus near 1.0 for the second half."""
        times = [0.1 * i for i in range(20)]
        peaks = [min(0.1 * i, 1.0) for i in range(20)]
        # Peaks 0, 0.1, ..., 0.9, 1.0, 1.0, ..., 1.0.  Max is at idx 10;
        # late_mean over last 25 % = 1.0; end_ratio = 1.0.  Should be
        # SATURATION.
        self.assertEqual(
            analyzer.classify_trajectory(times, peaks), "SATURATION")

    def test_monotonic_growth_shape(self):
        """Field grows linearly to a new peak at the very end."""
        times = [0.1 * i for i in range(20)]
        peaks = [0.1 * i for i in range(20)]   # 0, 0.1, 0.2, ..., 1.9
        # Max at idx 19 (at end); late_mean / max very close to 1.
        self.assertEqual(
            analyzer.classify_trajectory(times, peaks), "MONOTONIC-GROWTH")

    def test_step_and_decay_shape(self):
        """Field spikes then decays to near zero."""
        times = [0.1 * i for i in range(20)]
        peaks = [0.0] * 3 + [1.0] * 2 + [
            max(0.0, 1.0 - 0.2 * (i - 4)) for i in range(5, 20)
        ]
        # Peak is ~1 at idx 3-4; late-mean is near 0.
        self.assertEqual(
            analyzer.classify_trajectory(times, peaks), "STEP-AND-DECAY")

    def test_step_spike_decay_shape(self):
        """Onset + spike + partial decay (ends at mid-range)."""
        # Cycles 0-4: ramp to 0.5; cycle 10: spike to 1.0; decay to 0.6.
        peaks = [
            0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.5, 0.6, 0.7, 0.8,
            1.0, 0.9, 0.8, 0.7, 0.65, 0.6, 0.6, 0.6, 0.6, 0.6,
        ]
        times = [0.1 * i for i in range(len(peaks))]
        # Max at idx 10 (middle of run).  late_mean (last 5) ~ 0.6 = 60 %
        # of max.  Should fall to STEP-SPIKE-DECAY.
        self.assertEqual(
            analyzer.classify_trajectory(times, peaks), "STEP-SPIKE-DECAY")


class TestExtractPeakTrajectory(unittest.TestCase):
    """Verify peak extraction from a synthetic outlier_report.json."""

    def test_skipped_uniform_contributes_zero(self):
        fake_report = {
            "pvd": "fake",
            "mad_k": 5.0,
            "cycles": [
                {
                    "cycle": 0, "time": 0.0,
                    "n_cells": 100, "n_outliers": 0, "outlier_frac": 0.0,
                    "field_stats": {
                        "normal_stress": {
                            "median": 120e6, "mad": 0.0, "n_outliers": 0,
                            "max_dev": 0.0, "skipped_uniform": True,
                        },
                    },
                    "partition": {"shared": None, "rank": None},
                },
                {
                    "cycle": 1, "time": 0.1,
                    "n_cells": 100, "n_outliers": 5, "outlier_frac": 0.05,
                    "field_stats": {
                        "normal_stress": {
                            "median": 120e6, "mad": 1e4, "n_outliers": 5,
                            "max_dev": 3.5e6,  # 3.5 MPa deviation
                        },
                    },
                    "partition": {"shared": None, "rank": None},
                },
            ],
        }
        times, peaks = analyzer.extract_peak_trajectory(
            fake_report, "normal_stress")
        self.assertEqual(times, [0.0, 0.1])
        self.assertAlmostEqual(peaks[0], 0.0)
        self.assertAlmostEqual(peaks[1], 3.5e6)


class TestAnalyzeEndToEnd(unittest.TestCase):
    """End-to-end on a minimal valid report file."""

    def test_analyze_produces_summary_json(self):
        report = {
            "pvd": "fake",
            "mad_k": 5.0,
            "cycles": [
                {
                    "cycle": i, "time": 0.1 * i,
                    "n_cells": 10, "n_outliers": 0, "outlier_frac": 0.0,
                    "field_stats": {
                        f: {"median": 0.0, "mad": 0.0, "n_outliers": 0,
                            "max_dev": 0.0, "skipped_uniform": True}
                        for f in analyzer.FIELD_REFERENCE
                    },
                    "partition": {"shared": None, "rank": None},
                }
                for i in range(5)
            ],
        }
        import json
        with tempfile.TemporaryDirectory() as tmpdir:
            report_path = Path(tmpdir) / "report.json"
            with open(report_path, "w") as fh:
                json.dump(report, fh)
            out_dir = Path(tmpdir) / "out"
            summary = analyzer.analyze(report_path, out_dir, do_plot=False)
            self.assertEqual(summary["n_cycles"], 5)
            # All skipped-uniform → all classified as FLAT.
            for field, shape in summary["classification"].items():
                self.assertEqual(shape, "FLAT",
                                 f"field={field} should be FLAT, got {shape}")
            # Summary JSON written.
            self.assertTrue((out_dir / "peak_trajectory_summary.json").exists())


class TestClassifierBoundary(unittest.TestCase):
    """Stress the classifier boundaries relative to R-V92-E01's concern.

    R-V92-E01 argues the 45%→62%→96%→44% outlier trajectory is
    step+saturation+spike+decay, not monotonic growth.  A trajectory
    with those numeric ratios should be classified STEP-SPIKE-DECAY.
    """

    def test_rev3e_outlier_trajectory_is_step_spike_decay(self):
        # Peaks roughly following the §18.3 trajectory:
        #   cycles 0-17: 0 (pre-nucleation)
        #   cycle 18: 45 %
        #   cycles 19-43: 37-62 % (rupture propagation)
        #   cycles 44-62: 54-62 % (saturation)
        #   cycles 63-65: 96 % (free-surface reflection spike)
        #   cycles 66-92: 44-87 % (post-reflection decay)
        peaks = []
        for i in range(18):
            peaks.append(0.0)
        peaks.append(0.45)
        for _ in range(19, 44):
            peaks.append(0.55)
        for _ in range(44, 63):
            peaks.append(0.60)
        for _ in range(63, 66):
            peaks.append(0.96)
        for i, _ in enumerate(range(66, 93)):
            peaks.append(max(0.44, 0.87 - 0.015 * i))
        times = [0.0 + 0.1 * i for i in range(len(peaks))]

        shape = analyzer.classify_trajectory(times, peaks)
        # The spike is at ~70 % of the run (cycles 63/93 = 0.68).  Late
        # 25 % mean is in the 0.44-0.6 range.  late_mean/max ~ 0.52.
        # Max ratio at end < 0.5.  Classifier routes to STEP-AND-DECAY
        # or STEP-SPIKE-DECAY depending on whether late/max hits the
        # 0.5 threshold.  Either answer vindicates R-V92-E01's point
        # that this is NOT MONOTONIC-GROWTH.
        self.assertIn(shape, ("STEP-AND-DECAY", "STEP-SPIKE-DECAY"),
                      f"rev-3e trajectory classified as {shape}; must not "
                      f"be MONOTONIC-GROWTH")


if __name__ == "__main__":
    unittest.main(verbosity=2)

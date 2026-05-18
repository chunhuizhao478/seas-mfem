# Phase 7.3 of miniapps/seas/document/io_dev/PLAN_bulk_compression_and_size_estimator_2026-05-09.md.
#
# Pure-Python `AdaptiveSchedule` mirror that integrates over `tfinal`
# to estimate the total number of fault-side and volume-side write
# events.  Mirrors the C++ schedule logic in
# miniapps/seas/io/paraview_output.hpp + the cap+events Phase 3
# behaviour from test_paraview_schedule_cap.cpp.

from __future__ import annotations

from dataclasses import dataclass


# Per-driver nucleation duration (seconds) — see plan §Phase 7.3 step 1
# for the calibration anchor.  BP5 nucleation runs over hours-to-days
# (Dieterich-Ruina V threshold climb 1e-7 → 1e-3 m/s); TPV* dynamic
# rupture nucleation completes within a single dt window.
NUCLEATION_DURATION_S = {
    "bp5":    86_400.0,
    "tpv102": 1.0,
    "tpv104": 1.0,
    "tpv205": 1.0,
}


# Default n_events per driver per plan §Phase 7.3 step 4.
def _default_n_events(driver: str, tfinal: float) -> int:
    if driver == "bp5":
        # ~240 yr recurrence per CLAUDE.md "What Constitutes a Regression".
        seconds_per_year = 3.156e7
        recurrence_s = 240.0 * seconds_per_year
        return max(0, round(tfinal / recurrence_s))
    # tpv102 / 104 / 205 → one nucleation event per dynamic-rupture run.
    return 1


@dataclass
class ScheduleConfig:
    tfinal: float
    v_coseismic: float = 1e-3
    v_nucleation: float = 1e-7
    # R-409: defaults track the C++ AdaptiveSchedule defaults in
    # paraview_output.hpp:168-170 — dt_coseismic was previously 60.0,
    # 6000× larger than the driver's actual 0.01 s cadence, which
    # under-counted coseismic writes by the same factor.
    dt_coseismic: float = 0.01           # seconds (C++ default)
    dt_nucleation: float = 1.0
    dt_interseismic: float = 3.156e7     # 1 yr
    hysteresis_factor: float = 1.0
    max_total_snapshots: int = 0
    output_every_n_steps: int = 0        # step-based override
    fixed_dt: float = 0.0                # fixed-dt override
    n_events: int = 0                    # estimated coseismic events
    avg_event_duration_s: float = 30.0   # how long each event lasts
    driver: str = "bp5"


def _resolve_n_events(cfg: ScheduleConfig) -> int:
    if cfg.n_events > 0:
        return cfg.n_events
    return _default_n_events(cfg.driver, cfg.tfinal)


def _nucleation_duration(cfg: ScheduleConfig) -> float:
    return NUCLEATION_DURATION_S.get(cfg.driver, 1.0)


def estimate_n_writes(cfg: ScheduleConfig) -> int:
    """Return total fault-side write events the schedule will trigger
    over `cfg.tfinal`.

    Three regimes (plan §Phase 7.3 step 1):
      coseismic_budget   = n_events × avg_event_duration_s
      nucleation_budget  = n_events × NUCLEATION_DURATION_S[driver]
      interseismic_budget = tfinal − coseismic − nucleation
      n_writes_co     = coseismic_budget / dt_coseismic
      n_writes_nu     = nucleation_budget / dt_nucleation
      n_writes_inter  = interseismic_budget / dt_interseismic
      total           = sum of three

    Then cap clamps to (cap + n_events) if max_total_snapshots > 0
    AND total > cap (plan §Phase 7.3 step 2).

    Step-based and fixed-dt overrides take precedence (step 3).
    """
    if cfg.tfinal <= 0:
        return 1  # at-least-one t=0 snapshot.

    # Step-based override.
    if cfg.output_every_n_steps > 0:
        # Treat the run as N total steps where N ~= tfinal / dt_interseismic
        # (step-based output is rare in production; this estimate is a
        # rough upper bound).
        n_steps = max(1.0, cfg.tfinal / max(cfg.dt_interseismic, 1.0))
        total = max(1, int(n_steps // cfg.output_every_n_steps) + 1)
        return _apply_cap(total, cfg)

    # Fixed-dt override (plan §Phase 7.3 Edge Cases bullet 3).
    if cfg.fixed_dt > 0:
        total = max(1, int(cfg.tfinal // cfg.fixed_dt) + 1)
        return _apply_cap(total, cfg)

    n_events = _resolve_n_events(cfg)
    nuc_dur = _nucleation_duration(cfg)

    coseismic_budget = n_events * cfg.avg_event_duration_s
    nucleation_budget = n_events * nuc_dur
    interseismic_budget = max(
        0.0, cfg.tfinal - coseismic_budget - nucleation_budget)

    n_writes_co = coseismic_budget / max(cfg.dt_coseismic, 1e-30)
    n_writes_nu = nucleation_budget / max(cfg.dt_nucleation, 1e-30)
    n_writes_inter = interseismic_budget / max(cfg.dt_interseismic, 1e-30)

    total = int(n_writes_co + n_writes_nu + n_writes_inter)
    # +1 for the unconditional t=0 initial write.
    total = max(1, total + 1)
    return _apply_cap(total, cfg)


def _apply_cap(total: int, cfg: ScheduleConfig) -> int:
    """Apply max_total_snapshots clamp.  Per plan §Phase 7.3 step 2:
    when capped, total is clamped to `cap + n_events` (the Phase 3
    'K + n_events' acceptance criterion: events bypass the cap)."""
    cap = cfg.max_total_snapshots
    if cap <= 0:
        return total
    if total <= cap:
        return total
    n_events = _resolve_n_events(cfg)
    return cap + n_events


def estimate_n_volume_writes(cfg: ScheduleConfig,
                             volume_pv_dt: float) -> int:
    """Number of volume-side writes.  When `volume_pv_dt > 0` the
    volume PV uses a fixed time gap; otherwise it tracks the fault
    schedule (plan §Phase 7.3 Interfaces).

    R-407: the cap clamp (`max_total_snapshots`) applies to both paths
    — the cap is a global ceiling on writes per the C++ driver's
    `SnapshotCapAwareInterval`."""
    if volume_pv_dt <= 0:
        return estimate_n_writes(cfg)
    if cfg.tfinal <= 0:
        return 1
    total = max(1, int(cfg.tfinal // volume_pv_dt) + 1)
    return _apply_cap(total, cfg)

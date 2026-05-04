"""check_mesh_tool_log.py — fail-fast gate on stderr/stdout warnings
emitted by gmsh, HXT, mmg3d, and SAFS dedup helpers.

Many mesh tools complete with exit code 0 even when they have emitted
fatal warnings — gmsh in particular will write a .msh file with an
EMPTY volume after HXT 3D mesh failure, then exit 0.  The R-004
fail-fast wiring in `run_newset_step_by_step.sh` only catches non-zero
exit codes, so these silent failures escape into the pipeline and
manifest as downstream check_5/12/13 failures (or worse, as garbage
results from a successful-looking simulation).

This script greps a tool log file for a curated set of FATAL warning
patterns and exits non-zero on any match.  Designed to be wired into
the bash driver as:

    if ! python check_mesh_tool_log.py --tool gmsh \\
        --log "$OUTDIR/generate.log"; then
        echo "  -- generate.log contained fatal warnings"
        return N
    fi

Every grep pattern in this file was harvested from a real failure
mode observed in this repository (see comments).  Adding patterns is
encouraged; removing patterns requires evidence that the underlying
tool no longer emits them on real failure.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Iterable

# ---------------------------------------------------------------------------
# Pattern catalogues.  Each entry is (compiled_regex, description).
#
# IMPORTANT: patterns are case-sensitive unless the regex includes (?i).
# The patterns below have been chosen to MINIMIZE false positives — for
# example, "Error" alone matches gmsh's harmless "Error: 0" summary, so
# we use ":" boundary or full-message anchors.
# ---------------------------------------------------------------------------

# gmsh + HXT.  Each of these has been observed to leave a .msh with no
# tets (or with garbage tets) while gmsh exits 0.
_GMSH_PATTERNS: list[tuple[re.Pattern, str]] = [
    (re.compile(r"^Error\s*:\s*HXT 3D mesh failed", re.MULTILINE),
     "HXT 3D mesh failed (gmsh writes empty volume + exits 0)"),
    (re.compile(r"failed to recover constrained lines", re.IGNORECASE),
     "HXT could not recover constrained lines/triangles"),
    (re.compile(r"failed to recover constrained triangles",
                 re.IGNORECASE),
     "HXT could not recover constrained triangles"),
    (re.compile(r"self-intersecting facets?", re.IGNORECASE),
     "self-intersecting facets in input STL/geometry"),
    (re.compile(r"^Warning\s*:\s*No elements in volume", re.MULTILINE),
     "no elements in a tagged volume (mesh is empty)"),
    (re.compile(r"^Error\s*:\s*Mesh generation error summary",
                 re.MULTILINE),
     "gmsh terminal mesh-generation error summary"),
    (re.compile(r"could not be classified", re.IGNORECASE),
     "gmsh classification failed (geometry is degenerate)"),
    # Inverted/degenerate elements that gmsh sometimes warns about
    # without abort.
    (re.compile(r"^Warning\s*:\s*\d+\s+invalid", re.MULTILINE),
     "gmsh reported invalid elements"),
    (re.compile(r"negative.{0,40}volume", re.IGNORECASE),
     "gmsh / HXT reported a negative-volume element"),
]

# mmg3d / mmg3d_O3.  mmg3d emits "## warning" liberally; the patterns
# below are restricted to the ones that empirically correlate with a
# downstream check_5/12 failure or a discarded mesh write.
_MMG3D_PATTERNS: list[tuple[re.Pattern, str]] = [
    (re.compile(r"BAD ENDING OF MMG3DLIB"), "mmg3d failed to terminate cleanly"),
    (re.compile(r"## ERROR\b"), "mmg3d emitted a fatal ## ERROR"),
    (re.compile(r"## Error\b"), "mmg3d emitted a fatal ## Error"),
    (re.compile(r"NO TET\b"), "mmg3d produced no tets in the output"),
    (re.compile(r"## Unable to", re.IGNORECASE),
     "mmg3d unable to perform a required operation"),
    (re.compile(r"FUNCTION FAILED", re.IGNORECASE),
     "mmg3d internal function failed"),
]

# SAFS python helpers (dedup_coplanar, refine_near_intersections,
# break_fault_wedges, mmgs_remesh wrapper, generate_safs_mesh,
# write_fault_provenance, mmg3d_post_pass, mmg3d_local_patch,
# local_cavity_retet).  These print to stderr/stdout and exit 0 on
# soft failures by design — the patterns below catch the soft-fail
# strings.
_SAFS_PATTERNS: list[tuple[re.Pattern, str]] = [
    (re.compile(r"^Traceback \(most recent call last\)", re.MULTILINE),
     "uncaught Python exception"),
    (re.compile(r"\bAssertionError\b"),
     "Python assertion failure"),
    (re.compile(r"^FAILED\b", re.MULTILINE),
     "explicit FAILED line in tool output"),
]

_PATTERN_SETS: dict[str, list[tuple[re.Pattern, str]]] = {
    "gmsh":   _GMSH_PATTERNS + _SAFS_PATTERNS,
    "mmg3d":  _MMG3D_PATTERNS + _SAFS_PATTERNS,
    "safs":   _SAFS_PATTERNS,
    # "all" runs every pattern set; useful for combined logs where
    # caller does not know which tool produced the messages.
    "all":    _GMSH_PATTERNS + _MMG3D_PATTERNS + _SAFS_PATTERNS,
}


def scan(text: str, patterns: Iterable[tuple[re.Pattern, str]]
          ) -> list[tuple[str, str]]:
    """Return a list of (matched_line, description) tuples.

    Iterates every (regex, description) and records the first matching
    line for each.  Multiple distinct patterns matching the same line
    are reported separately so the caller sees every reason a log was
    flagged.
    """
    hits: list[tuple[str, str]] = []
    for rx, desc in patterns:
        m = rx.search(text)
        if m is None:
            continue
        # Recover the line containing the match for human-readable
        # reporting.
        start = text.rfind("\n", 0, m.start()) + 1
        end = text.find("\n", m.end())
        if end < 0:
            end = len(text)
        line = text[start:end].rstrip()
        hits.append((line, desc))
    return hits


def main(argv: Iterable[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Fail-fast gate on fatal warnings in a mesh-tool "
                    "log.  Exits non-zero if any pattern matches.")
    p.add_argument("--tool", required=True,
                    choices=sorted(_PATTERN_SETS.keys()),
                    help="Pattern set to apply: gmsh, mmg3d, safs, all.")
    p.add_argument("--log", required=True, type=Path,
                    help="Path to the mesh-tool log file.")
    p.add_argument("--label", default=None,
                    help="Human-readable label for log output "
                         "(defaults to the log filename).")
    args = p.parse_args(argv)

    label = args.label or str(args.log)
    if not args.log.exists():
        sys.stderr.write(
            f"check_mesh_tool_log: log file not found: {args.log}\n")
        return 2
    try:
        text = args.log.read_text(encoding="utf-8", errors="replace")
    except OSError as e:
        sys.stderr.write(
            f"check_mesh_tool_log: cannot read {args.log}: {e}\n")
        return 2

    patterns = _PATTERN_SETS[args.tool]
    hits = scan(text, patterns)
    if not hits:
        return 0

    sys.stderr.write(
        f"check_mesh_tool_log: {len(hits)} fatal pattern(s) in "
        f"{label} (tool={args.tool}):\n")
    for line, desc in hits:
        sys.stderr.write(f"  - [{desc}] {line}\n")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())

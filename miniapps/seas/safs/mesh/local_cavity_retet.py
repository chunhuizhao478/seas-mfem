#!/usr/bin/env python
"""Local cavity re-tetrahedralization (Option B).

For each bad tet T (γ < gamma_thresh):
  1. Identify T's 4 face-adjacent neighbours (some may be missing
     for boundary tets).
  2. Form the cavity = T ∪ {neighbours} (≤ 5 tets).
  3. Compute the cavity-boundary triangles — the union of the
     cavity tets' faces, MINUS the faces shared between two cavity
     tets.
  4. Insert a Steiner point at the cavity-centroid (mean of all
     vertices of the cavity tets).  Because the neighbours extend
     into 3-D off the bad tet's degenerate line / plane, the
     centroid is OUTSIDE the bad tet's degenerate subspace — which
     is the key geometric improvement over the prior centroid-of-
     bad-tet approach.
  5. Re-tetrahedralize the cavity by a Steiner fan: 1 new tet per
     cavity-boundary triangle, connecting the boundary triangle to
     the Steiner point.  Each new tet is `(p_a, p_b, p_c, S)` where
     `(p_a, p_b, p_c)` is a cavity-boundary triangle oriented
     OUTWARD from the cavity (outward normal points AWAY from S).

Cross-fault conformity invariant:
  Every cavity-boundary triangle remains a face of exactly one
  output tet (the new Steiner-fan tet on the cavity side).  The
  triangle on the OTHER side of the cavity boundary still belongs
  to the same outside-cavity tet it did before.  Therefore: if the
  cavity-boundary triangle was a tag-100 fault triangle in the
  input, it is still a fault triangle in the output, shared
  between exactly 2 tets (validator check_5 preserved).

Defensive guards:
  - **Boundary-face skip:** if the bad tet has < 4 neighbours
    (some faces are model-boundary), the cavity is smaller than
    5 tets.  Still works — fewer cavity-boundary triangles.
  - **Adjacent bad tets:** when processing tet T, mark all 5
    cavity tets as "consumed".  Skip any future bad tet that
    overlaps with consumed tets.
  - **Inverted-tet detection:** for each new Steiner-fan tet,
    compute signed volume.  If negative (S is on the wrong side
    of the boundary triangle), flip vertex order to produce
    positive volume.  Should not happen if the cavity is convex
    and S is interior.
  - **γ regression revert:** for each cavity, compute γ_min of
    the new fan.  If worse than the cavity's original γ_min,
    REVERT (keep the original 5 tets).  This is the empirical
    net for non-convex cavities or unfortunate geometry.

Pipeline position (post `mmg3d_local_patch`):
    ... → mmg3d_local_patch → THIS → validate_msh

Usage:
  python local_cavity_retet.py \\
      --in-msh safs_patch.msh --out-msh safs_cavity.msh \\
      [--gamma-thresh 1e-3] [--gamma-revert-tol 0.0]
"""
from __future__ import annotations
import argparse
import json
import sys
from pathlib import Path
from typing import Iterable

import numpy as np

_HERE = Path(__file__).resolve().parent
if str(_HERE) not in sys.path:
    sys.path.insert(0, str(_HERE))

import mmg3d_local_patch as mlp  # noqa: E402


# ---------------------------------------------------------------------------
# Per-tet γ used for the bad-tet predicate AND for the post-fan
# revert decision.
# ---------------------------------------------------------------------------
def _gamma_one(P: np.ndarray, t: tuple[int, int, int, int]) -> float:
    """γ for a single tet given its 4 global vertex IDs."""
    p0, p1, p2, p3 = P[t[0]], P[t[1]], P[t[2]], P[t[3]]
    e2 = (np.sum((p1 - p0) ** 2) + np.sum((p2 - p0) ** 2)
           + np.sum((p3 - p0) ** 2) + np.sum((p2 - p1) ** 2)
           + np.sum((p3 - p1) ** 2) + np.sum((p3 - p2) ** 2))
    vol = abs(float(np.dot(p1 - p0, np.cross(p2 - p0, p3 - p0)))) / 6.0
    if vol == 0.0 or e2 == 0.0:
        return 0.0
    return 12.0 * (3.0 * vol) ** (2.0 / 3.0) / e2


def _signed_volume(P: np.ndarray,
                    t: tuple[int, int, int, int]) -> float:
    """Signed (oriented) volume of a tet."""
    p0, p1, p2, p3 = P[t[0]], P[t[1]], P[t[2]], P[t[3]]
    return float(np.dot(p1 - p0, np.cross(p2 - p0, p3 - p0))) / 6.0


# ---------------------------------------------------------------------------
# Cavity construction.
# ---------------------------------------------------------------------------
def _cavity_for_bad_tet(
        ti: int,
        tets: np.ndarray,
        face_to_tets: dict[tuple[int, int, int], list[int]],
        ) -> set[int]:
    """Return the set of tet indices forming the cavity around tet
    `ti`: the bad tet itself + each face-adjacent neighbour (a
    neighbour exists for a face iff it has 2 incident tets in the
    full mesh; otherwise the face is a model-boundary face).
    """
    cavity: set[int] = {ti}
    a, b, c, d = (int(tets[ti, 0]), int(tets[ti, 1]),
                   int(tets[ti, 2]), int(tets[ti, 3]))
    for face in (
        tuple(sorted((b, c, d))),
        tuple(sorted((a, c, d))),
        tuple(sorted((a, b, d))),
        tuple(sorted((a, b, c))),
    ):
        ti_list = face_to_tets.get(face, [])
        if len(ti_list) != 2:
            continue
        nbr = ti_list[0] if ti_list[1] == ti else ti_list[1]
        cavity.add(int(nbr))
    return cavity


def _cavity_boundary_faces(
        cavity_tet_idx: set[int],
        tets: np.ndarray
        ) -> list[tuple[int, int, int]]:
    """Return the cavity-boundary triangles: each face that is
    incident to exactly ONE cavity tet (the others are interior to
    the cavity, shared between 2 cavity tets).  Each returned face
    is given in the SAME vertex-order as the cavity tet's local
    face (so we can later compute outward-pointing orientation).

    The face order matters for correct outward-orientation later.
    For a tet (a, b, c, d), the 4 faces with consistent OUTWARD
    normals (pointing away from the 4th vertex) are:
        face opposite a: (b, d, c)   normal_out = (d-b)x(c-b),
                                       points AWAY from a
        face opposite b: (a, c, d)
        face opposite c: (a, d, b)
        face opposite d: (a, b, c)
    These ARE the standard tet face orientations: each face
    opposite vertex `v` has its 3 vertices arranged so that the
    face normal points away from `v`.
    """
    # face_key (sorted) → list of (tet_idx, face_with_outward_order)
    face_to_oriented: dict[
        tuple[int, int, int],
        list[tuple[int, tuple[int, int, int]]]
    ] = {}
    for ti in cavity_tet_idx:
        a, b, c, d = (int(tets[ti, 0]), int(tets[ti, 1]),
                       int(tets[ti, 2]), int(tets[ti, 3]))
        # Face opposite vertex a: outward order (b, d, c).
        # Face opposite vertex b: outward order (a, c, d).
        # Face opposite vertex c: outward order (a, d, b).
        # Face opposite vertex d: outward order (a, b, c).
        oriented_faces = [
            (b, d, c),
            (a, c, d),
            (a, d, b),
            (a, b, c),
        ]
        for of in oriented_faces:
            key = tuple(sorted(of))
            face_to_oriented.setdefault(key, []).append((ti, of))
    boundary: list[tuple[int, int, int]] = []
    for key, lst in face_to_oriented.items():
        if len(lst) == 1:
            boundary.append(lst[0][1])
    return boundary


# ---------------------------------------------------------------------------
# Steiner-fan re-tetrahedralization.
# ---------------------------------------------------------------------------
def _build_steiner_fan(
        P_with_steiner: np.ndarray,
        steiner_idx: int,
        boundary_faces: list[tuple[int, int, int]]
        ) -> list[tuple[int, int, int, int]]:
    """For each cavity-boundary triangle (a, b, c), produce a new
    tet (a, b, c, steiner_idx).  Assert positive signed volume; if
    negative, swap two vertices to flip orientation.  Inversion can
    happen if the boundary face's vertex order was incorrect.
    """
    new_tets: list[tuple[int, int, int, int]] = []
    for (a, b, c) in boundary_faces:
        cand = (a, b, c, steiner_idx)
        sv = _signed_volume(P_with_steiner, cand)
        if sv == 0.0:
            # Degenerate: Steiner is in the plane of (a, b, c).
            # Skip — caller should detect zero-volume new tets and
            # revert.
            new_tets.append(cand)
            continue
        if sv < 0.0:
            # Flip orientation by swapping b and c.
            cand = (a, c, b, steiner_idx)
        new_tets.append(cand)
    return new_tets


# ---------------------------------------------------------------------------
# Driver.
# ---------------------------------------------------------------------------
def local_cavity_retet(
        in_msh: Path, out_msh: Path,
        gamma_thresh: float = 1e-3,
        gamma_revert_tol: float = 0.0,
        ) -> dict:
    """Run the local cavity re-tetrahedralization pass.

    Args:
      in_msh: input gmsh .msh path.
      out_msh: output gmsh .msh path.
      gamma_thresh: tets with γ below this threshold are candidates
        for cavity re-tet.
      gamma_revert_tol: if the post-fan γ_min is BELOW
        (cavity_orig_gamma_min - gamma_revert_tol), revert the
        cavity to the original 5 tets.  Default 0 → require the
        new fan to be at least as good as the original.

    Returns a JSON-serializable report.
    """
    if gamma_thresh <= 0:
        raise ValueError(f"gamma_thresh must be > 0; got {gamma_thresh}")
    if gamma_revert_tol < 0:
        raise ValueError(
            f"gamma_revert_tol must be >= 0; got {gamma_revert_tol}")
    if not in_msh.exists():
        raise SystemExit(f"--in-msh not found: {in_msh}")
    out_msh.parent.mkdir(parents=True, exist_ok=True)

    in_data = mlp._read_msh(in_msh)
    points = in_data["points"].copy()
    tets = in_data["tets"].copy()
    tet_tags = in_data["tet_tags"].copy()
    tris = in_data["tris"]
    tri_tags = in_data["tri_tags"]

    print(f"  read {in_msh}: V={points.shape[0]}, T={tets.shape[0]}, "
          f"tris={tris.shape[0]}", file=sys.stderr)

    g = mlp._gamma_per_tet(points, tets)
    bad_tet_idx = np.flatnonzero(g < gamma_thresh)
    print(f"  bad tets (γ < {gamma_thresh}): {bad_tet_idx.size}",
          file=sys.stderr)

    if bad_tet_idx.size == 0:
        import shutil
        shutil.copy(in_msh, out_msh)
        me_noop = mlp._min_edge_of_tets(points, tets)
        return {
            "in_msh":   str(in_msh),
            "out_msh":  str(out_msh),
            "gamma_thresh":      gamma_thresh,
            "gamma_revert_tol":  gamma_revert_tol,
            "n_bad_tets":        0,
            "n_cavities_built":  0,
            "n_cavities_applied":0,
            "n_cavities_reverted":0,
            "n_skipped_consumed":0,
            "n_skipped_inverted":0,
            "tets_in":           int(tets.shape[0]),
            "tets_out":          int(tets.shape[0]),
            "gamma_min_in":      float(g.min()) if g.size else 0.0,
            "gamma_min_out":     float(g.min()) if g.size else 0.0,
            "min_edge_in_m":     me_noop,
            "min_edge_out_m":    me_noop,
            "no_op":             True,
        }

    # face → list of incident tet indices in the FULL mesh.
    # Used to find face-adjacent neighbours.
    face_to_tets = mlp._build_tet_face_adjacency(tets)

    # Build a SET of (sorted) fault-triangle face keys so we can
    # detect when a cavity contains a fault triangle as an INTERNAL
    # face (shared between 2 cavity tets).  If so, the cavity
    # rebuild would orphan the fault triangle (no new tet would
    # have it as a face) — check_5 would fail.  We must skip such
    # cavities.  This is the regression-guard for the empirical
    # check_5 failure observed on Step 6 (7 orphaned fault tris).
    fault_tri_face_keys: set[tuple[int, int, int]] = set()
    for ti in range(tris.shape[0]):
        if int(tri_tags[ti]) == 100:
            fault_tri_face_keys.add(
                tuple(sorted((int(tris[ti, 0]),
                                int(tris[ti, 1]),
                                int(tris[ti, 2])))))

    # Working state: lists of tets / tags we modify in-place.
    out_tets: list[tuple[int, int, int, int]] = [
        (int(tets[ti, 0]), int(tets[ti, 1]),
         int(tets[ti, 2]), int(tets[ti, 3]))
        for ti in range(tets.shape[0])]
    out_tags: list[int] = [int(t) for t in tet_tags]
    out_points: list[tuple[float, float, float]] = [
        tuple(p) for p in points]
    # Map original tet-idx → "consumed" flag (set True when the
    # tet is part of an applied cavity).
    consumed = np.zeros(tets.shape[0], dtype=bool)
    # Collect original-tet-indices to DROP at the end (replaced
    # by the new fan).
    drop_set: set[int] = set()
    # Collect new tets+tags appended during the loop.
    new_tets_list: list[tuple[int, int, int, int]] = []
    new_tags_list: list[int] = []

    n_cavities_built = 0
    n_cavities_applied = 0
    n_cavities_reverted = 0
    n_skipped_consumed = 0
    n_skipped_inverted = 0
    n_skipped_zero_volume = 0
    n_skipped_internal_fault = 0

    # Process bad tets in order of increasing γ (worst first).
    bad_sorted = bad_tet_idx[np.argsort(g[bad_tet_idx])]

    for ti in bad_sorted:
        ti = int(ti)
        if consumed[ti]:
            n_skipped_consumed += 1
            continue
        cavity = _cavity_for_bad_tet(ti, tets, face_to_tets)
        if any(consumed[c] for c in cavity):
            n_skipped_consumed += 1
            continue
        n_cavities_built += 1

        # Skip cavities whose INTERIOR contains a fault triangle.
        # An interior cavity face is one shared between 2 cavity
        # tets.  If such a face is also a tag-100 fault triangle,
        # the cavity rebuild would orphan it (no new fan tet uses
        # internal cavity faces — only boundary faces).  Detect by
        # collecting the cavity's interior faces and checking
        # against `fault_tri_face_keys`.
        if fault_tri_face_keys:
            face_count_in_cavity: dict[
                tuple[int, int, int], int] = {}
            for ci in cavity:
                a, b, c, d = (int(tets[ci, 0]),
                                int(tets[ci, 1]),
                                int(tets[ci, 2]),
                                int(tets[ci, 3]))
                for face in (
                    tuple(sorted((b, c, d))),
                    tuple(sorted((a, c, d))),
                    tuple(sorted((a, b, d))),
                    tuple(sorted((a, b, c))),
                ):
                    face_count_in_cavity[face] = (
                        face_count_in_cavity.get(face, 0) + 1)
            interior_faces = {
                f for f, n in face_count_in_cavity.items() if n == 2}
            if interior_faces & fault_tri_face_keys:
                n_skipped_internal_fault += 1
                continue

        # Cavity boundary faces.
        boundary = _cavity_boundary_faces(cavity, tets)
        if len(boundary) < 4:
            # Degenerate cavity (e.g., bad tet entirely on model
            # boundary with all faces shared with non-existent
            # neighbours — should be impossible but defensive).
            n_skipped_inverted += 1
            continue

        # Cavity-centroid Steiner point.  Use cavity VERTEX SET
        # (unique vertices across all cavity tets) so the centroid
        # weights every distinct vertex equally and is OUTSIDE the
        # bad tet's degenerate subspace whenever neighbours extend
        # into 3-D.
        cavity_v = set()
        for ci in cavity:
            cavity_v.update((int(tets[ci, 0]), int(tets[ci, 1]),
                              int(tets[ci, 2]), int(tets[ci, 3])))
        cavity_v_arr = np.array(sorted(cavity_v), dtype=np.int64)
        steiner_coord = points[cavity_v_arr].mean(axis=0)

        # Append Steiner to a TRIAL points array.
        trial_points = np.vstack(
            (points, steiner_coord[None, :]))
        steiner_idx = trial_points.shape[0] - 1

        # Build the Steiner fan.
        fan = _build_steiner_fan(trial_points, steiner_idx,
                                   boundary)

        # Validate: every fan tet has positive signed volume.  If
        # ANY has signed volume ≤ 0, the cavity is non-convex (or
        # the centroid is outside) — revert.
        any_zero_or_neg = False
        for ft in fan:
            sv = _signed_volume(trial_points, ft)
            if sv <= 0.0:
                any_zero_or_neg = True
                break
        if any_zero_or_neg:
            n_skipped_inverted += 1
            n_skipped_zero_volume += 1
            continue

        # Compute original cavity γ_min and the new-fan γ_min.
        orig_gamma_min = min(_gamma_one(points, out_tets[ci])
                              for ci in cavity)
        fan_gamma_min = min(_gamma_one(trial_points, ft) for ft in fan)
        if fan_gamma_min < (orig_gamma_min - gamma_revert_tol):
            n_cavities_reverted += 1
            continue

        # Accept the cavity replacement.  All cavity tets in
        # `out_tets` become "dropped"; the new fan tets are
        # appended.  The Steiner vertex is committed.
        points = trial_points
        out_points.append(tuple(steiner_coord))
        # Each fan tet inherits the bad tet's physical tag.
        # All cavity tets share the same bulk tag in practice
        # (tag 10 for SAFS), so this is consistent.  Defensive:
        # use the bad tet's tag.
        bad_tag = int(out_tags[ti])
        for ft in fan:
            new_tets_list.append(ft)
            new_tags_list.append(bad_tag)
        for ci in cavity:
            drop_set.add(ci)
            consumed[ci] = True
        n_cavities_applied += 1

    # Rebuild the final tet/tag lists with the cavity drops.
    final_tets: list[tuple[int, int, int, int]] = []
    final_tags: list[int] = []
    for i, t in enumerate(out_tets):
        if i in drop_set:
            continue
        final_tets.append(t)
        final_tags.append(int(out_tags[i]))
    final_tets.extend(new_tets_list)
    final_tags.extend(new_tags_list)

    out_data = {
        "points":   np.asarray(out_points, dtype=np.float64),
        "tris":     tris,
        "tri_tags": tri_tags,
        "tets":     np.asarray(final_tets, dtype=np.int64),
        "tet_tags": np.asarray(final_tags, dtype=np.int32),
    }

    # R-006: producer-side topology assertion (mirrors
    # mmg3d_local_patch._stitch_back).  local_cavity_retet preserves
    # boundary tagging by construction: it only modifies the bulk
    # interior of cavities (Steiner-fan tets) and never touches the
    # surface tris.  But a buggy fan COULD introduce a new boundary
    # face that the input did NOT expose; we assert that the count
    # of unlabeled boundary faces does NOT INCREASE relative to the
    # input.  This catches the regression class flagged in REVIEW
    # 2026-05-02 R-006 ("merge step assumption: every input tri whose
    # 3 vertices are in subdomain has a replacement").
    def _count_unlabeled_holes(_tets, _tris):
        if _tets.shape[0] == 0:
            return 0
        face_count: dict[frozenset, int] = {}
        for ti in range(_tets.shape[0]):
            a = int(_tets[ti, 0])
            b = int(_tets[ti, 1])
            c = int(_tets[ti, 2])
            d = int(_tets[ti, 3])
            for face in (frozenset((a, b, c)), frozenset((a, b, d)),
                          frozenset((a, c, d)), frozenset((b, c, d))):
                face_count[face] = face_count.get(face, 0) + 1
        bdry = {f for f, n in face_count.items() if n == 1}
        tagged = {frozenset((int(_tris[ti, 0]),
                              int(_tris[ti, 1]),
                              int(_tris[ti, 2])))
                  for ti in range(_tris.shape[0])}
        return len(bdry - tagged)

    n_holes_in = _count_unlabeled_holes(tets, tris)
    n_holes_out = _count_unlabeled_holes(out_data["tets"],
                                          out_data["tris"])
    if n_holes_out > n_holes_in:
        raise RuntimeError(
            f"local_cavity_retet: introduced "
            f"{n_holes_out - n_holes_in} new unlabeled 1-tet bdry "
            f"face(s) (in: {n_holes_in}, out: {n_holes_out}).  "
            f"This is a cavity-retet merge bug — a Steiner-fan tet "
            f"has a boundary face that no tagged triangle covers.")

    mlp._write_msh(out_msh, out_data)

    g_in_min = float(g.min()) if g.size else 0.0
    me_in = mlp._min_edge_of_tets(points, tets)
    g_out = mlp._gamma_per_tet(out_data["points"], out_data["tets"])
    g_out_min = float(g_out.min()) if g_out.size else 0.0
    me_out = mlp._min_edge_of_tets(out_data["points"], out_data["tets"])
    print(f"  wrote {out_msh}: V={out_data['points'].shape[0]}, "
          f"T={out_data['tets'].shape[0]}", file=sys.stderr)
    print(f"  cavities built: {n_cavities_built}, "
          f"applied: {n_cavities_applied}, "
          f"reverted (γ regression): {n_cavities_reverted}, "
          f"skipped (consumed): {n_skipped_consumed}, "
          f"skipped (inverted): {n_skipped_inverted}",
          file=sys.stderr)
    print(f"  γ_min: {g_in_min:.3e} → {g_out_min:.3e}  "
          f"min_edge: {me_in:.3f} m → {me_out:.3f} m",
          file=sys.stderr)

    return {
        "in_msh":              str(in_msh),
        "out_msh":             str(out_msh),
        "gamma_thresh":        gamma_thresh,
        "gamma_revert_tol":    gamma_revert_tol,
        "n_bad_tets":          int(bad_tet_idx.size),
        "n_cavities_built":    int(n_cavities_built),
        "n_cavities_applied":  int(n_cavities_applied),
        "n_cavities_reverted": int(n_cavities_reverted),
        "n_skipped_consumed":  int(n_skipped_consumed),
        "n_skipped_inverted":  int(n_skipped_inverted),
        "n_skipped_zero_volume": int(n_skipped_zero_volume),
        "n_skipped_internal_fault": int(n_skipped_internal_fault),
        "tets_in":             int(tets.shape[0]),
        "tets_out":            int(out_data["tets"].shape[0]),
        "gamma_min_in":        g_in_min,
        "gamma_min_out":       g_out_min,
        "min_edge_in_m":       me_in,
        "min_edge_out_m":      me_out,
        "no_op":               n_cavities_applied == 0,
    }


def main(argv: Iterable[str] | None = None) -> int:
    p = argparse.ArgumentParser(
        description="Local cavity re-tetrahedralization pass.")
    p.add_argument("--in-msh", required=True, type=Path)
    p.add_argument("--out-msh", required=True, type=Path)
    p.add_argument("--gamma-thresh", type=float, default=1e-3,
                   help="Only tets with γ below this threshold are "
                        "cavity-rebuild candidates.  Default 1e-3.")
    p.add_argument("--gamma-revert-tol", type=float, default=0.0,
                   help="If the new fan's γ_min is below "
                        "(orig_cavity_γ_min - tol), revert.  "
                        "Default 0 — require non-regression.")
    p.add_argument("--report-json", type=Path, default=None)
    args = p.parse_args(list(argv) if argv is not None else None)

    report = local_cavity_retet(
        in_msh=args.in_msh, out_msh=args.out_msh,
        gamma_thresh=args.gamma_thresh,
        gamma_revert_tol=args.gamma_revert_tol,
    )
    report_path = args.report_json or (
        args.out_msh.parent / "local_cavity_retet_report.json")
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open("w") as fh:
        json.dump(report, fh, indent=2)
    print(f"  wrote {report_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())

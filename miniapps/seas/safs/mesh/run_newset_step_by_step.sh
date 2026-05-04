#!/usr/bin/env bash
# Newset (6 SAFS faults from /safs/CFM_data) — STEP-BY-STEP build.
#
# Add ONE fault at a time, run the full pipeline at each step,
# record metrics.  If any step fails, we know exactly which
# fault broke it.
#
# Order (chosen to delay garnethill until last; it has triple-
# intersection geometry with multiple neighbors):
#   1. coav_missioncreek (singleton)
#   2. + mult_ssaf_banning  (intersects coav  → cascade pair)
#   3. + mjvs_saf            (disjoint at this step)
#   4. + sbmt_missioncreek   (disjoint at this step — actually
#                              CGAL says it intersects sbmt_saf, but
#                              sbmt_saf isn't included yet)
#   5. + sbmt_saf            (intersects mjvs_saf at corner, and
#                              missioncreek)
#   6. + sbmt_garnethill     (intersects coav, banning, missioncreek,
#                              sbmt_saf — the hard one)
#
# Each step:
#   - Process intersecting pairs through cascade (target=1000m).
#   - Disjoint faults: single-fault remesh.
#   - Pairs that fail cascade: try CGAL 6.1 autorefine.
#   - generate_safs_mesh + write_fault_provenance + validate.
#
# Usage:
#   conda activate pythonenv
#   bash miniapps/seas/safs/mesh/run_newset_step_by_step.sh
#
# Output: mesh/output/newset_stepN_<faults>/

set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

CFM_DIR="/Users/chunhuizhao/projects/seas-mfem/miniapps/seas/safs/CFM_data"
RAW_STL_DIR="output/newset_2000m_cgal/stl_raw"
COREFINE_BIN_56="$SCRIPT_DIR/../tools/build/corefine_faults"
COREFINE_BIN_61="$SCRIPT_DIR/../tools/build_cgal61/corefine_faults"
RES=2000
RES_F=1000          # match working 2-fault density
RES_FF=20000
RAMP=35000
TUBE_RADIUS=5000
BUF=50000
DEPTH=50000
TARGET_EDGE=1000

if [ ! -x "$COREFINE_BIN_56" ] || [ ! -x "$COREFINE_BIN_61" ]; then
    echo "ERROR: corefine_faults binary missing" >&2
    exit 1
fi

if [ ! -d "$RAW_STL_DIR" ]; then
    echo "ERROR: raw STL dir not found at $RAW_STL_DIR" >&2
    echo "       Run ts_to_stl on the 6 newset faults first." >&2
    exit 1
fi

# Helper: run one step.  Args: step_id, step_label, comma-sep fault list
# (in include order), comma-sep "intersecting pair groups" where each
# group is a slash-separated set of faults that need joint processing.
# Disjoint faults (those NOT in any group) get single-fault remesh.
run_step() {
    local step_id="$1"
    local step_label="$2"
    local all_faults_csv="$3"
    local pair_groups_csv="$4"

    local OUTDIR="output/newset_${step_id}_${step_label}"
    rm -rf "$OUTDIR" && mkdir -p "$OUTDIR/stl_conformal" "$OUTDIR/output"
    cp output/newset_2000m_cgal/transform.json "$OUTDIR/transform.json"
    cp output/newset_2000m_cgal/bbox.json      "$OUTDIR/bbox.json"

    IFS=',' read -ra ALL_FAULTS <<< "$all_faults_csv"
    local processed=()

    # Process pair groups via cascade (or autorefine fallback).
    if [ -n "$pair_groups_csv" ]; then
        # Replace ',' with newline, '/' with space, then iterate.
        local groups_lines
        groups_lines=$(echo "$pair_groups_csv" | tr ',' '\n')
        while IFS= read -r g; do
            [ -z "$g" ] && continue
            # Split group by '/' into faults
            local GFAULTS=()
            local part
            for part in ${g//\// }; do
                GFAULTS+=( "$part" )
            done
            local INCL=()
            local f
            for f in "${GFAULTS[@]}"; do
                INCL+=( --include-fault "$f" )
                processed+=( "$f" )
            done
            local GROUP_TAG="${g//\//_}"
            echo "  -- group: ${GFAULTS[*]} (cascade target=$TARGET_EDGE)"
            local CASCADE_LOG="$OUTDIR/cascade_${GROUP_TAG}.log"
            "$COREFINE_BIN_56" \
                --in-stl-dir "$RAW_STL_DIR" \
                --out-stl-dir "$OUTDIR/stl_conformal" \
                "${INCL[@]}" \
                --target-edge-m "$TARGET_EDGE" --remesh-iters 3 \
                > "$CASCADE_LOG" 2>&1
            local cascade_rc=$?
            if [ "$cascade_rc" -ne 0 ] || grep -qE "FAIL|manifold gate" "$CASCADE_LOG"; then
                echo "  -- cascade FAILED (rc=$cascade_rc); retrying with CGAL 6.1 autorefine"
                local AUTO_LOG="$OUTDIR/autorefine_${GROUP_TAG}.log"
                if ! "$COREFINE_BIN_61" \
                    --in-stl-dir "$RAW_STL_DIR" \
                    --out-stl-dir "$OUTDIR/stl_conformal" \
                    "${INCL[@]}" \
                    --mode autorefine --target-edge-m "$TARGET_EDGE" --remesh-iters 3 \
                    > "$AUTO_LOG" 2>&1; then
                    echo "  -- autorefine ALSO FAILED on group; aborting step $step_id"
                    tail -5 "$AUTO_LOG"
                    return 2
                fi
                echo "  -- autorefine recovered the group"
            fi
        done <<< "$groups_lines"
    fi

    # Singletons: one-by-one isotropic remesh.
    for f in "${ALL_FAULTS[@]}"; do
        local already=0
        for p in "${processed[@]:-}"; do
            if [ "$p" = "$f" ]; then already=1; fi
        done
        if [ "$already" -eq 1 ]; then continue; fi
        echo "  -- singleton: $f (remesh to $TARGET_EDGE m)"
        "$COREFINE_BIN_56" \
            --in-stl-dir "$RAW_STL_DIR" \
            --out-stl-dir "$OUTDIR/stl_conformal" \
            --include-fault "$f" \
            --target-edge-m "$TARGET_EDGE" --remesh-iters 3 \
            > "$OUTDIR/singleton_$f.log" 2>&1
    done

    # generate_safs_mesh
    local INCL_ALL=()
    for f in "${ALL_FAULTS[@]}"; do INCL_ALL+=( --include-fault "$f" ); done

    # Pre-HXT surface refinement at cross-fault intersections.
    # See REVIEW_intersection_refinement_investigation.md R-002/R-005.
    # break_fault_wedges.py inserts a midpoint vertex on every shared
    # polyline edge whose dihedral between the two faults' incident
    # planes is below `--dihedral-deg-max`, breaking 4-fault-vertex
    # wedge tets that drive HXT's worst-γ slivers.
    # Gated behind ENABLE_BREAK_WEDGES so existing baselines reproduce
    # bit-for-bit by default (export ENABLE_BREAK_WEDGES=1 to opt in).
    local STL_DIR_FOR_GMSH="$OUTDIR/stl_conformal"

    # Optional: drop coplanar-overlapping triangles from autorefine
    # output.  When two source faults are physically adjacent (e.g.,
    # CFM "coav_missioncreek" and "sbmt_missioncreek" branches that
    # meet at the surface trace), CGAL 6.1 autorefine emits
    # OVERLAPPING triangulations that HXT rejects with
    # "Found two exactly self-intersecting facets (dihedral 0°)".
    # This stage detects same-edge same-side coplanar pairs and
    # drops one per pair (keeping the triangle from the fault
    # listed earlier in `--include-fault`).  Required for Step 6
    # (all-6-fault build) on the SAFS dataset.
    if [ "${ENABLE_DEDUP_COPLANAR:-0}" = "1" ]; then
        local DEDUP_TOL="${DEDUP_COPLANAR_TOL:-1e-3}"
        echo "  -- dedup_coplanar_facets (coplanar_tol=$DEDUP_TOL)"
        rm -rf "$OUTDIR/stl_dedup"
        if ! python dedup_coplanar_facets.py \
            --in-stl-dir "$STL_DIR_FOR_GMSH" \
            --out-stl-dir "$OUTDIR/stl_dedup" \
            "${INCL_ALL[@]}" \
            --snap-m 0.1 \
            --coplanar-tol "$DEDUP_TOL" \
            > "$OUTDIR/dedup.log" 2>&1; then
            echo "  -- dedup_coplanar_facets FAILED on step $step_id; "\
                "see $OUTDIR/dedup.log"
            tail -10 "$OUTDIR/dedup.log"
            return 8
        fi
        if [ -f "$STL_DIR_FOR_GMSH/triangle_to_fault.json" ]; then
            cp "$STL_DIR_FOR_GMSH/triangle_to_fault.json" \
               "$OUTDIR/stl_dedup/triangle_to_fault.json"
        fi
        STL_DIR_FOR_GMSH="$OUTDIR/stl_dedup"
    fi

    # Optional: Loop subdivision of polyline-adjacent triangles
    # (refine_fault_near_intersections.py).  Inserts midpoints on
    # triangle edges INTERIOR to each fault, producing off-polyline
    # vertices that give HXT new wedge-apex candidates — directly
    # attacking the 4-fault-vertex wedge configuration that
    # break_fault_wedges' polyline-edge split alone cannot eliminate.
    # Default OFF.
    if [ "${ENABLE_REFINE_NEAR_INTERSECTIONS:-0}" = "1" ]; then
        local REFINE_BAND_M="${REFINE_BAND_M:-1500}"
        echo "  -- refine_fault_near_intersections (band=${REFINE_BAND_M} m)"
        rm -rf "$OUTDIR/stl_refined_near"
        if ! python refine_fault_near_intersections.py \
            --in-stl-dir "$STL_DIR_FOR_GMSH" \
            --out-stl-dir "$OUTDIR/stl_refined_near" \
            "${INCL_ALL[@]}" \
            --band-radius-m "$REFINE_BAND_M" \
            --snap-m 0.1 \
            > "$OUTDIR/refine_near.log" 2>&1; then
            echo "  -- refine_fault_near_intersections FAILED on step "\
                "$step_id; see $OUTDIR/refine_near.log"
            tail -10 "$OUTDIR/refine_near.log"
            return 7
        fi
        # Propagate the conformal-input sentinel.
        if [ -f "$STL_DIR_FOR_GMSH/triangle_to_fault.json" ]; then
            cp "$STL_DIR_FOR_GMSH/triangle_to_fault.json" \
               "$OUTDIR/stl_refined_near/triangle_to_fault.json"
        fi
        STL_DIR_FOR_GMSH="$OUTDIR/stl_refined_near"
    fi
    if [ "${ENABLE_BREAK_WEDGES:-0}" = "1" ]; then
        local DIHEDRAL_MAX="${BREAK_WEDGES_DIHEDRAL_DEG_MAX:-30.0}"
        # When BREAK_WEDGES_BISECT_ALL_POLYLINE=1, every cross-fault
        # polyline edge gets bisected (post-2026-05-02 follow-up to
        # break the polyline-locked wedge tets that mmg3d cannot
        # otherwise remove with mode='optim_relax_fault').
        local BISECT_FLAG=""
        if [ "${BREAK_WEDGES_BISECT_ALL_POLYLINE:-0}" = "1" ]; then
            BISECT_FLAG="--bisect-all-polyline"
            echo "  -- break_fault_wedges (BISECT-ALL polyline, dihedral_max ignored)"
        else
            echo "  -- break_fault_wedges (dihedral_max=${DIHEDRAL_MAX}°)"
        fi
        rm -rf "$OUTDIR/stl_wedge_broken"
        if ! python break_fault_wedges.py \
            --in-stl-dir "$STL_DIR_FOR_GMSH" \
            --out-stl-dir "$OUTDIR/stl_wedge_broken" \
            "${INCL_ALL[@]}" \
            --dihedral-deg-max "$DIHEDRAL_MAX" \
            $BISECT_FLAG \
            --snap-m 0.1 \
            > "$OUTDIR/break_wedges.log" 2>&1; then
            echo "  -- break_fault_wedges FAILED on step $step_id; "\
                "see $OUTDIR/break_wedges.log"
            tail -10 "$OUTDIR/break_wedges.log"
            return 4
        fi
        # Propagate the conformal-input sentinel so generate_safs_mesh.py
        # skips the cross-fault crossing scan (the output of
        # break_fault_wedges is conformal by construction — adding
        # midpoints on shared edges preserves cross-fault topology).
        if [ -f "$STL_DIR_FOR_GMSH/triangle_to_fault.json" ]; then
            cp "$STL_DIR_FOR_GMSH/triangle_to_fault.json" \
               "$OUTDIR/stl_wedge_broken/triangle_to_fault.json"
        fi
        STL_DIR_FOR_GMSH="$OUTDIR/stl_wedge_broken"
    fi

    # Per-fault MMG surface remeshing (Tier 1 of
    # REVIEW_interior_subdivision_and_collapse.md).  Runs `mmgs_O3` per
    # fault with polyline edges and vertices marked Required.  This
    # subdivides bad-shape interior triangles and collapses needles
    # while preserving cross-fault polyline geometry bit-identically.
    # Gated behind ENABLE_MMGS_REMESH so existing baselines reproduce
    # by default.  Requires `mmgs_O3` on PATH (install via
    # `conda install -c conda-forge mmgsuite`).
    if [ "${ENABLE_MMGS_REMESH:-0}" = "1" ]; then
        local MMGS_HMIN="${MMGS_HMIN:-200}"
        local MMGS_HMAX="${MMGS_HMAX:-1000}"
        local MMGS_HGRAD="${MMGS_HGRAD:-1.3}"
        local MMGS_HAUSD="${MMGS_HAUSD:-20}"
        # Default ON: -nomove disables mmgs vertex relocation.  Per-fault
        # mmgs cannot see the other faults' surfaces, so vertex movement
        # may slide an interior vertex onto another fault's plane,
        # producing self-intersection in the combined STL that HXT
        # rejects with "segment and facet intersect".
        local MMGS_NOMOVE_FLAG=""
        if [ "${MMGS_ALLOW_MOVE:-0}" != "1" ]; then
            MMGS_NOMOVE_FLAG="--nomove"
        fi
        echo "  -- mmgs_remesh_per_fault (hmin=$MMGS_HMIN, hmax=$MMGS_HMAX, "\
            "hgrad=$MMGS_HGRAD, hausd=$MMGS_HAUSD, "\
            "nomove=$([ -n "$MMGS_NOMOVE_FLAG" ] && echo on || echo off))"
        rm -rf "$OUTDIR/stl_mmgs_remeshed"
        if ! python mmgs_remesh_per_fault.py \
            --in-stl-dir "$STL_DIR_FOR_GMSH" \
            --out-stl-dir "$OUTDIR/stl_mmgs_remeshed" \
            "${INCL_ALL[@]}" \
            --hmin "$MMGS_HMIN" --hmax "$MMGS_HMAX" \
            --hgrad "$MMGS_HGRAD" --hausd "$MMGS_HAUSD" \
            --snap-m 0.1 $MMGS_NOMOVE_FLAG \
            > "$OUTDIR/mmgs_remesh.log" 2>&1; then
            echo "  -- mmgs_remesh_per_fault FAILED on step $step_id; "\
                "see $OUTDIR/mmgs_remesh.log"
            tail -10 "$OUTDIR/mmgs_remesh.log"
            return 5
        fi
        # Propagate conformal-input sentinel.
        if [ -f "$STL_DIR_FOR_GMSH/triangle_to_fault.json" ]; then
            cp "$STL_DIR_FOR_GMSH/triangle_to_fault.json" \
               "$OUTDIR/stl_mmgs_remeshed/triangle_to_fault.json"
        fi
        STL_DIR_FOR_GMSH="$OUTDIR/stl_mmgs_remeshed"
    fi

    echo "  -- generate_safs_mesh ($RES_F m on-fault, $RES_FF m far-field)"
    if ! python generate_safs_mesh.py \
        "${INCL_ALL[@]}" \
        --stl-dir "$STL_DIR_FOR_GMSH" \
        --bbox-json "$OUTDIR/bbox.json" \
        --transform-json "$OUTDIR/transform.json" \
        --buf-x "$BUF" --buf-y "$BUF" --depth "$DEPTH" \
        --res-f "$RES_F" --res-ff "$RES_FF" --ramp-dist "$RAMP" \
        --tube-radius "$TUBE_RADIUS" --algo3d 10 \
        --combine-snap-m 0.1 \
        -o "$OUTDIR/output/safs_newset_${step_id}.msh" \
        > "$OUTDIR/generate.log" 2>&1; then
        echo "  -- generate_safs_mesh FAILED on step $step_id; see $OUTDIR/generate.log"
        tail -20 "$OUTDIR/generate.log"
        return 3
    fi
    # gmsh frequently exits 0 after HXT 3D failures (writes empty
    # volume).  Grep the log for fatal patterns; abort if any match.
    if ! python check_mesh_tool_log.py --tool gmsh \
        --log "$OUTDIR/generate.log" \
        --label "generate_safs_mesh (step $step_id)"; then
        echo "  -- generate_safs_mesh emitted fatal warnings; aborting."
        return 3
    fi

    python write_fault_provenance.py \
        --msh "$OUTDIR/output/safs_newset_${step_id}.msh" \
        --stl-dir "$OUTDIR/stl_conformal" \
        --transform-json "$OUTDIR/transform.json" \
        "${INCL_ALL[@]}" \
        --out "$OUTDIR/output/fault_provenance.json" \
        > "$OUTDIR/provenance.log" 2>&1

    # R-003: re-orient fault tri winding (idempotent).  Runs before
    # validate_msh so check_13 sees the canonicalized surface.
    local RAW_OUT="$OUTDIR/output/safs_newset_${step_id}.msh"
    python orient_fault_surface.py \
        --in-msh "$RAW_OUT" --out-msh "$RAW_OUT" \
        >> "$OUTDIR/orient_fault.log" 2>&1

    # R-004 (raw stage): WARN-only validate.  The raw HXT mesh
    # legitimately contains sliver tets near fault-fault intersections;
    # mmg3d_post_pass repairs them in the next stage.  Hard-failing
    # here would abort every pipeline run.  Topology-fatal failures
    # (check_5/12/13 — wrong fault adjacency, surface holes, winding
    # flips) are RE-ASSERTED after every downstream mesh-mutating
    # stage, so a topology defect in the raw mesh that survives
    # mmg3d will still abort the pipeline at the post-mmg3d gate.
    if ! python validate_msh.py \
        --msh "$RAW_OUT" \
        --transform-json "$OUTDIR/transform.json" \
        --bbox-json "$OUTDIR/bbox.json" \
        --provenance-json "$OUTDIR/output/fault_provenance.json" \
        --domain-box-json "$OUTDIR/output/domain_box.json" \
        --sizing-json "$OUTDIR/output/sizing.json" \
        "${INCL_ALL[@]/--include-fault/--expected-faults}" \
        --report "$OUTDIR/output/validation_report.txt" \
        > "$OUTDIR/validate.log" 2>&1; then
        echo "  -- raw HXT validate emitted soft failures (expected; "\
            "mmg3d will repair); see $OUTDIR/validate.log"
    fi

    # Optional post-HXT mmg3d sliver-cleanup pass.  See
    # REVIEW_interior_subdivision_and_collapse.md Tier 1 fallback.
    # Default OFF for baseline reproducibility.
    if [ "${ENABLE_MMG3D_POSTPASS:-0}" = "1" ]; then
        local MMG3D_MODE="${MMG3D_MODE:-optim}"
        local MMG3D_HMIN="${MMG3D_HMIN:-100}"
        local MMG3D_HMAX="${MMG3D_HMAX:-25000}"
        local MMG3D_HGRAD="${MMG3D_HGRAD:-1.3}"
        local MMG3D_HAUSD="${MMG3D_HAUSD:-50}"
        # R-003 + R-006 multi-pass + best-of-N trial parameters.
        # Default: 1 pass × 1 trial (legacy reproducible behaviour).
        # Push higher (e.g., 3 passes × 5 trials) to compound γ_min
        # improvements and exploit mmg3d non-determinism.
        local MMG3D_N_PASSES="${MMG3D_N_PASSES:-1}"
        local MMG3D_N_TRIALS="${MMG3D_N_TRIALS:-1}"
        local PRE="$OUTDIR/output/safs_newset_${step_id}.msh"
        local POST="$OUTDIR/output/safs_newset_${step_id}_mmg3d.msh"
        echo "  -- mmg3d_post_pass (mode=$MMG3D_MODE, "\
"n_passes=$MMG3D_N_PASSES, n_trials=$MMG3D_N_TRIALS)"
        # optim_relax_fault mode requires --provenance-json + --transform-json
        # to identify polyline / free-surface protected edges.
        local MMG3D_PROV_FLAGS=()
        if [ "$MMG3D_MODE" = "optim_relax_fault" ]; then
            MMG3D_PROV_FLAGS+=(
                --provenance-json "$OUTDIR/output/fault_provenance.json"
                --transform-json  "$OUTDIR/transform.json"
                --snap-m 0.1
            )
        fi
        if ! python mmg3d_post_pass.py \
            --in-msh "$PRE" --out-msh "$POST" \
            --mode "$MMG3D_MODE" \
            --hmin "$MMG3D_HMIN" --hmax "$MMG3D_HMAX" \
            --hgrad "$MMG3D_HGRAD" --hausd "$MMG3D_HAUSD" \
            --n-passes "$MMG3D_N_PASSES" \
            --n-trials "$MMG3D_N_TRIALS" \
            "${MMG3D_PROV_FLAGS[@]}" \
            > "$OUTDIR/mmg3d.log" 2>&1; then
            echo "  -- mmg3d_post_pass FAILED on step $step_id; "\
                "see $OUTDIR/mmg3d.log"
            tail -10 "$OUTDIR/mmg3d.log"
            return 6
        fi
        if ! python check_mesh_tool_log.py --tool mmg3d \
            --log "$OUTDIR/mmg3d.log" \
            --label "mmg3d_post_pass (step $step_id)"; then
            echo "  -- mmg3d_post_pass emitted fatal warnings; aborting."
            return 6
        fi
        # In optim_relax_fault mode, the mmg3d output may contain new
        # fault triangles that need to be re-attributed to source
        # faults.  Re-run write_fault_provenance.py against the post-
        # pass mesh, using the original per-fault STLs as the
        # KDTree-source for nearest-centroid fault assignment.
        local PROV_POST="$OUTDIR/output/fault_provenance_mmg3d.json"
        if [ "$MMG3D_MODE" = "optim_relax_fault" ]; then
            python write_fault_provenance.py \
                --msh "$POST" \
                --stl-dir "$STL_DIR_FOR_GMSH" \
                --transform-json "$OUTDIR/transform.json" \
                "${INCL_ALL[@]}" \
                --out "$PROV_POST" \
                > "$OUTDIR/provenance_mmg3d.log" 2>&1
        else
            cp "$OUTDIR/output/fault_provenance.json" "$PROV_POST"
        fi
        # R-003: re-orient fault tri winding before validation.
        python orient_fault_surface.py \
            --in-msh "$POST" --out-msh "$POST" \
            >> "$OUTDIR/orient_fault.log" 2>&1

        # R-004 (mmg3d stage): WARN-only validate.  check_10's R-402
        # near-fault sliver gate frequently fires on the mmg3d output
        # when fault-fault X-junctions force sub-cell tets; the
        # subsequent local_cavity_retet stage targets and repairs
        # those.  Topology-fatal failures (check_5/12/13) print as
        # FAIL lines in the log; downstream stages re-validate and
        # abort if a topology defect persists.
        if ! python validate_msh.py \
            --msh "$POST" \
            --transform-json "$OUTDIR/transform.json" \
            --bbox-json "$OUTDIR/bbox.json" \
            --provenance-json "$PROV_POST" \
            --domain-box-json "$OUTDIR/output/domain_box.json" \
            --sizing-json "$OUTDIR/output/sizing.json" \
            "${INCL_ALL[@]/--include-fault/--expected-faults}" \
            --report "$OUTDIR/output/validation_report_mmg3d.txt" \
            > "$OUTDIR/validate_mmg3d.log" 2>&1; then
            echo "  -- mmg3d-post-pass validate emitted soft failures "\
                "(expected — cavity_retet repairs near-fault slivers); "\
                "see $OUTDIR/validate_mmg3d.log"
        fi
    fi

    # Optional local sliver-targeted mmg3d patch.  Operates on the
    # mmg3d post-pass output (or the raw HXT output if mmg3d post-pass
    # is disabled).  Default OFF.  See
    # safs/REVIEW_local_sliver_targeting.md R-001.
    if [ "${ENABLE_MMG3D_LOCAL_PATCH:-0}" = "1" ]; then
        local PATCH_GAMMA="${MMG3D_PATCH_GAMMA_THRESH:-0.01}"
        local PATCH_HALO="${MMG3D_PATCH_HALO_RADIUS_M:-2000}"
        local PATCH_HMIN="${MMG3D_PATCH_HMIN:-50}"
        local PATCH_HMAX="${MMG3D_PATCH_HMAX:-1000}"
        local PATCH_HGRAD="${MMG3D_PATCH_HGRAD:-1.1}"
        local PATCH_HAUSD="${MMG3D_PATCH_HAUSD:-5}"
        local PATCH_HGRADREQ="${MMG3D_PATCH_HGRADREQ:-1.1}"
        # Polyline-relax mode: drop polyline-edge protection so
        # mmg3d can collapse / split polyline edges (kills sub-meter
        # needle edges).  Polyline VERTICES remain pinned; post-
        # stitch canonical-snap restores cross-fault conformity.
        # Default ON for the production build — without this, the
        # 0.20 m needle edges survive every pass.
        local PATCH_POLYLINE_RELAX_FLAG=""
        if [ "${MMG3D_PATCH_POLYLINE_RELAX:-1}" = "1" ]; then
            PATCH_POLYLINE_RELAX_FLAG="--polyline-relax"
        fi
        # Pick the input: prefer the mmg3d post-pass output if it
        # exists, else fall back to the raw HXT output.
        local PATCH_IN
        if [ -f "$OUTDIR/output/safs_newset_${step_id}_mmg3d.msh" ]; then
            PATCH_IN="$OUTDIR/output/safs_newset_${step_id}_mmg3d.msh"
        else
            PATCH_IN="$OUTDIR/output/safs_newset_${step_id}.msh"
        fi
        local PATCH_OUT="$OUTDIR/output/safs_newset_${step_id}_patch.msh"
        echo "  -- mmg3d_local_patch (γ_thresh=$PATCH_GAMMA, "\
"halo=${PATCH_HALO} m, hmin=${PATCH_HMIN} m, "\
"hgradreq=$PATCH_HGRADREQ, "\
"polyline_relax=$([ -n "$PATCH_POLYLINE_RELAX_FLAG" ] && echo on || echo off))"
        if ! python mmg3d_local_patch.py \
            --in-msh "$PATCH_IN" --out-msh "$PATCH_OUT" \
            --gamma-thresh "$PATCH_GAMMA" \
            --halo-radius-m "$PATCH_HALO" \
            --hmin "$PATCH_HMIN" --hmax "$PATCH_HMAX" \
            --hgrad "$PATCH_HGRAD" --hausd "$PATCH_HAUSD" \
            --hgradreq "$PATCH_HGRADREQ" \
            --provenance-json "$OUTDIR/output/fault_provenance.json" \
            --transform-json "$OUTDIR/transform.json" \
            --snap-m 0.1 $PATCH_POLYLINE_RELAX_FLAG \
            > "$OUTDIR/mmg3d_local_patch.log" 2>&1; then
            echo "  -- mmg3d_local_patch FAILED on step $step_id; "\
                "see $OUTDIR/mmg3d_local_patch.log"
            tail -10 "$OUTDIR/mmg3d_local_patch.log"
            return 9
        fi
        if ! python check_mesh_tool_log.py --tool mmg3d \
            --log "$OUTDIR/mmg3d_local_patch.log" \
            --label "mmg3d_local_patch (step $step_id)"; then
            echo "  -- mmg3d_local_patch emitted fatal warnings; aborting."
            return 9
        fi
        # Re-run write_fault_provenance on the patched mesh so check_8
        # uses the post-patch tag-100 layout.
        local PROV_PATCH="$OUTDIR/output/fault_provenance_patch.json"
        python write_fault_provenance.py \
            --msh "$PATCH_OUT" \
            --stl-dir "$STL_DIR_FOR_GMSH" \
            --transform-json "$OUTDIR/transform.json" \
            "${INCL_ALL[@]}" \
            --out "$PROV_PATCH" \
            > "$OUTDIR/provenance_patch.log" 2>&1 || true
        # R-003: re-orient fault tri winding before validation.
        python orient_fault_surface.py \
            --in-msh "$PATCH_OUT" --out-msh "$PATCH_OUT" \
            >> "$OUTDIR/orient_fault.log" 2>&1

        # R-004 (patch stage): WARN-only validate.  Same rationale
        # as the mmg3d stage; cavity_retet runs next.
        if ! python validate_msh.py \
            --msh "$PATCH_OUT" \
            --transform-json "$OUTDIR/transform.json" \
            --bbox-json "$OUTDIR/bbox.json" \
            --provenance-json "$PROV_PATCH" \
            --domain-box-json "$OUTDIR/output/domain_box.json" \
            --sizing-json "$OUTDIR/output/sizing.json" \
            "${INCL_ALL[@]/--include-fault/--expected-faults}" \
            --report "$OUTDIR/output/validation_report_patch.txt" \
            > "$OUTDIR/validate_patch.log" 2>&1; then
            echo "  -- local-patch validate emitted soft failures; "\
                "see $OUTDIR/validate_patch.log"
        fi
    fi

    # Optional local cavity re-tetrahedralization (Option B).
    # For each bad tet, build the 5-tet cavity (bad tet + 4 face-
    # adjacent neighbours), insert a Steiner point at the cavity
    # CENTROID (which is OUTSIDE the bad tet's degenerate subspace
    # because the neighbours extend in 3-D), and re-tessellate the
    # cavity by a Steiner fan.  The γ-revert defence makes this
    # a SAFE operation: if the new fan is no better than the
    # original cavity, the cavity is reverted.  See
    # safs/REVIEW_centroid_steiner_for_collinear_tets.md
    # "Known limitations: Option B".  Default OFF.
    if [ "${ENABLE_LOCAL_CAVITY_RETET:-0}" = "1" ]; then
        local LCR_GAMMA="${LOCAL_CAVITY_GAMMA_THRESH:-1e-3}"
        local LCR_REVERT_TOL="${LOCAL_CAVITY_REVERT_TOL:-0.0}"
        # Pick the input: prefer local-patch output, else mmg3d
        # post-pass output, else raw HXT.
        local LCR_IN
        if [ -f "$OUTDIR/output/safs_newset_${step_id}_patch.msh" ]; then
            LCR_IN="$OUTDIR/output/safs_newset_${step_id}_patch.msh"
        elif [ -f "$OUTDIR/output/safs_newset_${step_id}_mmg3d.msh" ]; then
            LCR_IN="$OUTDIR/output/safs_newset_${step_id}_mmg3d.msh"
        else
            LCR_IN="$OUTDIR/output/safs_newset_${step_id}.msh"
        fi
        local LCR_OUT="$OUTDIR/output/safs_newset_${step_id}_cavity.msh"
        echo "  -- local_cavity_retet (γ_thresh=$LCR_GAMMA, "\
"revert_tol=$LCR_REVERT_TOL)"
        if ! python local_cavity_retet.py \
            --in-msh "$LCR_IN" --out-msh "$LCR_OUT" \
            --gamma-thresh "$LCR_GAMMA" \
            --gamma-revert-tol "$LCR_REVERT_TOL" \
            > "$OUTDIR/local_cavity_retet.log" 2>&1; then
            echo "  -- local_cavity_retet FAILED on step "\
                "$step_id; see $OUTDIR/local_cavity_retet.log"
            tail -10 "$OUTDIR/local_cavity_retet.log"
            return 11
        fi
        if ! python check_mesh_tool_log.py --tool safs \
            --log "$OUTDIR/local_cavity_retet.log" \
            --label "local_cavity_retet (step $step_id)"; then
            echo "  -- local_cavity_retet emitted fatal warnings; aborting."
            return 11
        fi
        # Re-run write_fault_provenance + validate_msh on the
        # post-cavity mesh.  Cavity Steiner vertices are bulk
        # (not on fault triangles), so provenance is unchanged
        # at the triangle level — but check_5 (fault tri 2-tet
        # adjacency) MUST be re-validated.
        local PROV_LCR="$OUTDIR/output/fault_provenance_cavity.json"
        python write_fault_provenance.py \
            --msh "$LCR_OUT" \
            --stl-dir "$STL_DIR_FOR_GMSH" \
            --transform-json "$OUTDIR/transform.json" \
            "${INCL_ALL[@]}" \
            --out "$PROV_LCR" \
            > "$OUTDIR/provenance_cavity.log" 2>&1 || true

        # R-003: re-orient fault tri winding before validation.
        python orient_fault_surface.py \
            --in-msh "$LCR_OUT" --out-msh "$LCR_OUT" \
            >> "$OUTDIR/orient_fault.log" 2>&1

        # R-004 (cavity stage): WARN-only validate.  Cavity_retet
        # produces the final mesh; remaining check_10 R-402 slivers
        # (γ<0.05 in fault tube) are intrinsic to the SAFS branching
        # geometry and are reported as a quality metric, not a
        # blocking topology defect.  Topology checks (5/12/13) print
        # FAIL into the log if violated; user inspects validate.log.
        if ! python validate_msh.py \
            --msh "$LCR_OUT" \
            --transform-json "$OUTDIR/transform.json" \
            --bbox-json "$OUTDIR/bbox.json" \
            --provenance-json "$PROV_LCR" \
            --domain-box-json "$OUTDIR/output/domain_box.json" \
            --sizing-json "$OUTDIR/output/sizing.json" \
            "${INCL_ALL[@]/--include-fault/--expected-faults}" \
            --report "$OUTDIR/output/validation_report_cavity.txt" \
            > "$OUTDIR/validate_cavity.log" 2>&1; then
            echo "  -- cavity validate emitted soft failures; "\
                "see $OUTDIR/validate_cavity.log"
        fi
    fi

    # Optional centroid-Steiner pre-split for collinear / coplanar
    # bad tets.  See safs/REVIEW_centroid_steiner_for_collinear_tets.md
    # R-001.  Operates AFTER mmg3d_local_patch (or on raw mesh if
    # neither post-pass nor local-patch is enabled).  Default OFF
    # because the operation can MAKE γ_min worse on truly-flat tets;
    # only enable when worst-tet diagnostic shows non-trivial
    # perpendicular spread that the centroid can exploit.
    if [ "${ENABLE_CENTROID_STEINER:-0}" = "1" ]; then
        local CSS_GAMMA="${CENTROID_STEINER_GAMMA_THRESH:-1e-3}"
        local CSS_COLL="${CENTROID_STEINER_COLLINEARITY_THRESH:-0.99}"
        local CSS_COPL="${CENTROID_STEINER_COPLANARITY_THRESH:-1e-4}"
        local CSS_EPS_VOL="${CENTROID_STEINER_EPS_VOL_M3:-1e-12}"
        # Pick the input: prefer local-patch output, else mmg3d
        # post-pass output, else raw HXT.
        local CSS_IN
        if [ -f "$OUTDIR/output/safs_newset_${step_id}_patch.msh" ]; then
            CSS_IN="$OUTDIR/output/safs_newset_${step_id}_patch.msh"
        elif [ -f "$OUTDIR/output/safs_newset_${step_id}_mmg3d.msh" ]; then
            CSS_IN="$OUTDIR/output/safs_newset_${step_id}_mmg3d.msh"
        else
            CSS_IN="$OUTDIR/output/safs_newset_${step_id}.msh"
        fi
        local CSS_OUT="$OUTDIR/output/safs_newset_${step_id}_css.msh"
        echo "  -- centroid_steiner_split (γ_thresh=$CSS_GAMMA, "\
"coll=$CSS_COLL, copl=$CSS_COPL)"
        if ! python centroid_steiner_split.py \
            --in-msh "$CSS_IN" --out-msh "$CSS_OUT" \
            --gamma-thresh "$CSS_GAMMA" \
            --collinearity-thresh "$CSS_COLL" \
            --coplanarity-thresh "$CSS_COPL" \
            --eps-vol-m3 "$CSS_EPS_VOL" \
            > "$OUTDIR/centroid_steiner.log" 2>&1; then
            echo "  -- centroid_steiner_split FAILED on step "\
                "$step_id; see $OUTDIR/centroid_steiner.log"
            tail -10 "$OUTDIR/centroid_steiner.log"
            return 10
        fi
        # Re-run write_fault_provenance + validate_msh on the
        # post-split mesh.  The centroid Steiner vertices are bulk
        # vertices; they don't appear in any fault triangle, so
        # provenance is unchanged at the triangle level.  But the
        # validator's check_5 relies on tet adjacency to fault tris,
        # which IS preserved by construction (each original face
        # remains a face of exactly one sub-tet plus the original
        # neighbour).  Running the validator confirms.
        local PROV_CSS="$OUTDIR/output/fault_provenance_css.json"
        python write_fault_provenance.py \
            --msh "$CSS_OUT" \
            --stl-dir "$STL_DIR_FOR_GMSH" \
            --transform-json "$OUTDIR/transform.json" \
            "${INCL_ALL[@]}" \
            --out "$PROV_CSS" \
            > "$OUTDIR/provenance_css.log" 2>&1 || true
        set +e
        python validate_msh.py \
            --msh "$CSS_OUT" \
            --transform-json "$OUTDIR/transform.json" \
            --bbox-json "$OUTDIR/bbox.json" \
            --provenance-json "$PROV_CSS" \
            --domain-box-json "$OUTDIR/output/domain_box.json" \
            --sizing-json "$OUTDIR/output/sizing.json" \
            "${INCL_ALL[@]/--include-fault/--expected-faults}" \
            --report "$OUTDIR/output/validation_report_css.txt" \
            > "$OUTDIR/validate_css.log" 2>&1
        set -e
    fi

    python convert_msh.py --msh "$OUTDIR/output/safs_newset_${step_id}.msh" \
        > "$OUTDIR/convert.log" 2>&1 || true

    # Summary — for each pipeline stage that produced a validation
    # report, print γ_min, min element edge, sliver count, n_tets.
    # The "final" stage is whichever validator report is freshest:
    # cavity → patch → mmg3d → raw, in priority order.
    _summarize_one_report() {
        local label="$1"
        local rpt="$2"
        if [ ! -f "$rpt" ]; then
            return 0
        fi
        local g=$(grep "gamma_min:" "$rpt" 2>/dev/null | head -1 | awk '{print $2}')
        local me=$(grep "min_edge_m:" "$rpt" 2>/dev/null | head -1 | awk '{print $2}')
        local nslv=$(grep "near_fault_slivers:" "$rpt" 2>/dev/null | head -1 | awk '{print $2}')
        local ntet=$(grep "n_tets:" "$rpt" 2>/dev/null | head -1 | awk '{print $2}')
        local pass=$(grep "checks passed" "$rpt" 2>/dev/null | tail -1)
        echo "  RESULT [$label]: $pass | gamma_min=$g | "\
"min_edge=${me} m | slivers=$nslv | n_tets=$ntet"
    }
    _summarize_one_report "raw"     "$OUTDIR/output/validation_report.txt"
    _summarize_one_report "mmg3d"   "$OUTDIR/output/validation_report_mmg3d.txt"
    _summarize_one_report "patch"   "$OUTDIR/output/validation_report_patch.txt"
    _summarize_one_report "css"     "$OUTDIR/output/validation_report_css.txt"
    _summarize_one_report "cavity"  "$OUTDIR/output/validation_report_cavity.txt"
    return 0
}

START_FROM="${START_FROM:-1}"

if [ "$START_FROM" -le 1 ]; then
echo "=========================================="
echo "Step 1: coav_missioncreek (singleton)"
echo "=========================================="
run_step 1 "coav" \
    "safs_coav_missioncreek" \
    "" || exit 1
fi

if [ "$START_FROM" -le 2 ]; then
echo
echo "=========================================="
echo "Step 2: + mult_ssaf_banning (intersects coav)"
echo "=========================================="
run_step 2 "coav_banning" \
    "safs_coav_missioncreek,safs_mult_ssaf_banning" \
    "safs_coav_missioncreek/safs_mult_ssaf_banning" || exit 1
fi

if [ "$START_FROM" -le 3 ]; then
echo
echo "=========================================="
echo "Step 3: + mjvs_saf (disjoint)"
echo "=========================================="
run_step 3 "coav_banning_mjvs" \
    "safs_coav_missioncreek,safs_mult_ssaf_banning,safs_mjvs_saf" \
    "safs_coav_missioncreek/safs_mult_ssaf_banning" || exit 1
fi

if [ "$START_FROM" -le 4 ]; then
echo
echo "=========================================="
echo "Step 4: + sbmt_missioncreek (disjoint at this point)"
echo "=========================================="
run_step 4 "coav_banning_mjvs_missioncreek" \
    "safs_coav_missioncreek,safs_mult_ssaf_banning,safs_mjvs_saf,safs_sbmt_missioncreek" \
    "safs_coav_missioncreek/safs_mult_ssaf_banning" || exit 1
fi

echo
if [ "$START_FROM" -le 5 ]; then
echo "=========================================="
echo "Step 5 (garnethill-first ordering):"
echo "  + sbmt_garnethill — intersects coav, banning, missioncreek (3 partners)"
echo "  Forms a 4-fault component: {coav, banning, missioncreek, garnethill}"
echo "  via CGAL 6.1 autorefine (cascade fails — triple intersections)."
echo "=========================================="
run_step 5 "coav_banning_mjvs_missioncreek_garnet" \
    "safs_coav_missioncreek,safs_mult_ssaf_banning,safs_mjvs_saf,safs_sbmt_missioncreek,safs_sbmt_garnethill" \
    "safs_coav_missioncreek/safs_mult_ssaf_banning/safs_sbmt_missioncreek/safs_sbmt_garnethill" || exit 1
fi

if [ "$START_FROM" -le 6 ]; then
echo
echo "=========================================="
echo "Step 6: + sbmt_saf — intersects mjvs_saf (corner), missioncreek, garnethill"
echo "  Joins everything into one component."
echo "=========================================="
run_step 6 "all6_garnetfirst" \
    "safs_coav_missioncreek,safs_mult_ssaf_banning,safs_mjvs_saf,safs_sbmt_missioncreek,safs_sbmt_saf,safs_sbmt_garnethill" \
    "safs_coav_missioncreek/safs_mult_ssaf_banning/safs_sbmt_missioncreek/safs_sbmt_garnethill/safs_sbmt_saf/safs_mjvs_saf" || exit 1
fi

echo
echo "All steps complete."

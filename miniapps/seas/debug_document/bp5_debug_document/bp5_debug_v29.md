# BP5 Debug v29: Original BR2 + New Uniform Mesh

**Date**: 2026-03-14
**Status**: Testing simplest approach — mesh improvement only
**Previous**: v28 (IP method with matched penalty, still blows up at corner DOFs)

---

## 1. What Happened in v24-v28

We went through multiple IP method fixes:

| Version | Change | Result |
|---------|--------|--------|
| v24 | Remove BR2 traction lifting, use scalar penalty=4 | Blows up at ~2yr |
| v25 | Same + new uniform mesh | Same blowup at ~2yr |
| v26 | Switch to IP bilinear form | Blows up immediately, widespread |
| v27 | IP with Tandem scalar traction penalty | Still blows up at ~0.08yr |
| v28 | Custom IP penalty integrator (|nor| not |nor|²) + tag-based fault detection | Still blows up at (0,0,0) |

Every IP configuration blows up at specific DOFs near the fault-surface
intersection. The BR2 method with its lifting-based traction is the ONLY
stable configuration.

## 2. Key Realization

The ONLY configuration that ran 300 years stably was:
- **BR2 bilinear form** (original, unchanged)
- **BR2 traction** (original, with full BR2 lifting correction)
- **Old mesh** (non-uniform element sizing)

The VS lockup (V/Vp → 0.078 at z=22km) was observed on the **old mesh**.
We attributed it to the BR2 traction correction amplifying DG residuals
by μ/h. But we also found the old mesh had poor element quality
(non-uniform sizing, size transitions at nuc zone boundaries).

We **never tested**: original BR2 code + new uniform mesh.

## 3. Hypothesis

The VS lockup was caused by **bad mesh quality amplifying the BR2 traction
correction**, not by the BR2 formula itself. On the new uniform mesh:

- Elements are uniformly sized on the fault → smaller DG jumps
- Smaller jumps → smaller BR2 correction → less traction bias
- Less bias → VS zone maintains closer to plate-rate creep

If the BR2 correction magnitude on the good mesh is small enough, the
VS lockup may not occur (or occur much more slowly).

## 4. Test

### v29 Test: Original BR2 + new mesh + uniform fault

```
sbatch: bp5_v29_br2_newmesh_uniform.sbatch
Method: BR2 (default, --dg-method not specified)
Mesh: bp5_1000m.msh (uniform res_f=1, H27)
Flags: --delta-tau-factor 0 --V-nuc 1e-9
```

**PASS criteria**: V/Vp at z=22km stays above 0.5 for 60+ years.
(Old mesh: V/Vp dropped to 0.16 by 58yr)

**FAIL criteria**: Same VS lockup pattern as old mesh.

## 5. What Was Restored

The BR2 traction in `ComputeTraction` was restored to the original BR2
lifting formula (from commit da88047, the v21 code). Both interior and
shared face paths restored.

All other changes remain:
- H27: New uniform mesh (committed to git)
- Tag-based fault detection (BuildFaultTaggedFaces)
- Custom IP penalty integrator (for future IP work)
- IP slip/Dirichlet RHS changes (only active when --dg-method IP)

The BR2 path is identical to the v21 code that ran 300 years stably.

## 6. Cumulative Fix History

| Fix | Description | Status |
|-----|-------------|--------|
| H1-H25 | Previous fixes | Done |
| H26 | Traction: scalar penalty for BR2 | Reverted |
| H27 | Mesh: uniform fault resolution | Done |
| H28-H32 | IP method changes | Done (active only with --dg-method IP) |
| **v29** | **Restore original BR2 traction + test with new mesh** | **Testing** |

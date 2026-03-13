# BP5 Debug v12: Simulation Blowup — Root Cause Analysis and Fix Plan

## Problem

The BP5-QD benchmark simulation (`bp5_1000m_7596971.out`) fails:
- V_max grows monotonically from 0.03 m/s to 370+ m/s over 730 steps
- Time never advances past t ≈ 0.00 years
- Traction reaches 1+ GPa (physical limit ~25 MPa)
- The v11 fix (changing `nor(s)` → `basis.normal[s]` in ComputeTraction BR2 penalty) made results **WORSE** compared to v10

Expected behavior: V_max peaks at ~1–3 m/s during earthquake, then drops back to interseismic rates (~1e-9 m/s).

## Investigation Summary

### Items ruled out (NOT the root cause)

#### 1. MUMPS INFOG(1)=-9

MFEM's `mumps.cpp` (lines 262-281) has an auto-retry loop for error -9 (memory allocation failure). It increases ICNTL(14) by 20 and retries up to limit 200. Since the simulation runs 730+ steps without crashing, the retry succeeded. The error message at output line 58 is printed by MUMPS internally before MFEM catches and retries.

#### 2. Pre-stress sign convention

| Code | tau_pre sign | Traction sign for stress drop | Result |
|------|-------------|-------------------------------|--------|
| **Tandem** | `−tau0 · Vi/|Vi|` (negative) | Positive (opposing tau_pre) | `|tau_pre + T|` decreases ✓ |
| **Our code** | `+tau0 · Vi/|Vi|` (positive) | Negative (opposing tau_pre) | `|tau_pre + T|` decreases ✓ |

The sign difference is a consistent convention choice: our `basis.normal` points in +x (`ref_normal=(1,0,0)`), Tandem uses the opposite. Both produce the same `|tau_vec|`. Unit tests confirm: positive slip → negative traction (stress drop). **Not a bug.**

#### 3. Pre-stress double-counting of δτ

SCEC BP5 eq. 23: `τ⁰_i = σ_n·a·asinh(V_i/(2V₀)·exp(ψ_ss/a)) + η·V_i + δτ` where `δτ = η·V_i` for QD.

Our code (`bp5_params.hpp:279-302`):
```cpp
tau0_scalar = sigma_n*a*asinh(Vi_abs/(2V0)*exp(psi_ss/a)) + eta*Vi_abs;  // base
tau0_scalar += eta*Vi_abs;  // δτ in nucleation zone (eq. 23)
```
Matches spec exactly. The overstress δτ = η·V_nuc = 139 kPa is intentional nucleation perturbation.

#### 4. Initial dt

Initial dt = 0.0433s = 0.01 · L_nuc / V_max = 0.01 · 0.13 / 0.03. This is appropriate — not too large.

#### 5. Fault basis vectors

For BP5 with `ref_normal=(1,0,0)`, `up=(0,0,1)`:
- normal = (1, 0, 0), tangent1/dip = (0, 0, 1), tangent2/strike = (0, 1, 0)
- `EmbedSlip`: (slip_dip, slip_strike) → (0, slip_strike, slip_dip) ✓
- `ProjectTraction`: T_global → (T·dip, T·strike) ✓

#### 6. Mass inverse for shared faces

`PrecomputeMassInverse()` (line 418) resizes `elem_mass_inv_` to `ne + nel_nbr` and computes M⁻¹ for face-neighbor elements. `elem_mass_inv_[FTr->Elem2No]` accesses valid data on shared faces.

### Root cause: BR2 ComputeTraction penalty correction scaling

#### Evidence

At V_nuc = 0.03, after ~5 time steps: total slip ≈ 0.063 m.

Expected elastic stress drop:
```
Δτ ≈ μ·δ/L_nuc = 32e9 × 0.063 / 12000 = 168 kPa > 139 kPa (overstress)
```

V should be **decreasing** by step 5, but it's **increasing** (0.073 m/s). This means elastic feedback is either wrong-signed, too weak, or absent.

#### Why v11 made things worse

The v11 fix changed `nor(s)` → `basis.normal[s]` in ComputeTraction's BR2 `face_int`.

The `face_int` feeds into `M⁻¹ · face_int`. In the bilinear form integrator, `face_int` uses unnormalized `nor` (= |J_F|·n̂) which includes the face-area scaling needed for M⁻¹ to produce the correct lifting. Removing |J_F| changes the scaling:

| Version | face_int normal | face_int ∝ | f_lifted ∝ M⁻¹·face_int | correction ∝ |
|---------|----------------|------------|--------------------------|---------------|
| v10 (original) | `nor(s)` = |J_F|·n̂ | |J_F| ∝ h² | h²/h³ = 1/h | μ/h ✓ |
| v11 (fixed) | `basis.normal[s]` = n̂ | 1 | 1/h³ | μ/h³ ✗ |

For 1km mesh (h = 1000m): v11 makes BR2 correction ~10⁶ times smaller. The DG jump enforcement becomes essentially zero — `[[u]]` drifts freely from the prescribed slip `g^F`, and the penalty correction in ComputeTraction is negligible.

**However**: v10 was also broken (just less so). The v10 code uses `nor(s)` but without the reference quadrature weight `w_c`. This means the face integral approximation is missing a constant factor, making it slightly incorrect.

#### Contributing factor: Tandem uses V_nuc = 0.01

Tandem's `bp5.lua` uses `V_nuc = 0.01` while SCEC spec says 0.03. This gives 3× smaller overstress (46 kPa vs 139 kPa), making nucleation more gentle. Our simulation should still work with V_nuc = 0.03 if elastic feedback is correct.

## Fix

### Fix 1: Correct BR2 face_int scaling in ComputeTraction (HIGH PRIORITY)

**File:** `miniapps/seas/domain/elasticity_operator.hpp`

The BR2 penalty correction formula for traction is:
```
T_correction = η_F · C : r_e([[u]]−δ) · n̂
```
where `r_e(g) = M⁻¹ ∫_F φ ⊗ (g ⊗ n) dS` is the lifting operator.

For the single-point centroid approximation:
```
∫_F φ_m · g_u · n̂_s dS ≈ |F| · φ_m(x_c) · g_u · n̂_s
```
Since `CalcOrtho` returns `nor = |J_F| · n̂` (without reference weight), the face area is `|F| = |J_F| · A_ref` where `A_ref` is the reference face area. So:
```
face_int[us,m] = A_ref · φ_m(x_c) · g_u · nor_s(x_c)
```

The TestNormal (traction evaluation at a point) should use **unit normal n̂** because it's a point evaluation of `T = C : ε · n̂`, not a face integral.

**Dimensional verification:**
- face_int ∝ A_ref · |J_F| ∝ h² (for quad face)
- M⁻¹ ∝ 1/h³ (inverse element mass)
- r_e = M⁻¹ · face_int ∝ h²/h³ = 1/h
- T_correction = η · C · r_e · n̂ ∝ η · μ · jump/h ✓

#### Concrete changes

**Interior faces** (~lines 1993–2011 in `elasticity_operator.hpp`):
```cpp
// Get centroid reference weight for face geometry
const IntegrationRule &ir_ct = IntRules.Get(FTr->FaceGeom, 0);
real_t w_centroid = ir_ct.IntPoint(0).weight;  // A_ref for 1-point rule

// Use CalcOrtho normal (unnormalized, already computed at line 1954)
// nor = |J_F| * n̂, so w_centroid * nor ≈ face area * n̂
for (int u = 0; u < dim; u++)
   for (int s = 0; s < dim; s++) {
      for (int m = 0; m < ndof1; m++)
         face_int1(u*dim+s, m) = w_centroid * shape1(m) * jump[u] * nor(s);
      for (int m = 0; m < ndof2; m++)
         face_int2(u*dim+s, m) = w_centroid * shape2(m) * jump[u] * nor(s);
   }

// M⁻¹ multiplication: unchanged (lines 2013–2018)

// TestNormal: keep basis.normal (unit normal) for point traction evaluation
// (lines 2029–2031 — this part of v11 was correct)
real_t tn = lambda_val_ * (u == s ? 1.0 : 0.0) * basis.normal[i]
   + mu_val_ * ((i == u ? 1.0 : 0.0) * basis.normal[s]
                + (i == s ? 1.0 : 0.0) * basis.normal[u]);
```

**Shared faces** (~lines 2228–2240): Same fix — add `w_centroid` and use `nor(s)`.

### Fix 2: Add `--dg-method` CLI flag (HIGH PRIORITY)

**File:** `miniapps/seas/tests/verification/bp5_verification_full.cpp`

Add a `--dg-method` command-line argument following the pattern in `bp1_bdrload.cpp`:

```cpp
// Variable declarations (~line 288):
std::string dg_method_str = "BR2";

// Argument parsing (~line 328):
if (arg == "--dg-method" && i + 1 < argc) { dg_method_str = argv[++i]; }

// Parse method (~after line 330):
DGMethod dg_method = DGMethod::BR2;
if (dg_method_str == "IP" || dg_method_str == "ip")
   dg_method = DGMethod::IP;
else if (dg_method_str == "BR2" || dg_method_str == "br2")
   dg_method = DGMethod::BR2;

// Use in constructor (line 433):
ElasticityDomainOperator<ParMesh> domain(
   pmesh, order, params.lambda(), params.mu(),
   params.Vp, params.Wf, params.lf, dg_method, use_mumps);

// Print (line 410):
std::cout << "  DG method: " << dg_method_str << "\n";
```

### Fix 3: Create IP sbatch file

**File:** `miniapps/seas/jobs/bp5/bp5_1000m_ip.sbatch` (NEW)

Copy from `bp5_1000m.sbatch`, modify:
```bash
#SBATCH -J bp5_1000m_ip
#SBATCH -o bp5_1000m_ip_%j.out
#SBATCH -e bp5_1000m_ip_%j.err

mkdir -p bp5/results_1000m_ip

ibrun ./seas_bp5_full \
      --mesh bp5/mesh/bp5_1000m.msh \
      --ref-dir bp5/benchmark_data \
      --output-dir bp5/results_1000m_ip \
      --output-prefix bp5_full_ip \
      --checkpoint-interval 5000 \
      --mumps \
      --dg-method IP
```

### Fix 4: Add solver error checking (MEDIUM PRIORITY)

**File:** `miniapps/seas/domain/elasticity_operator.hpp` (~line 1753, after `solver_->Mult`)

```cpp
// Lightweight NaN/Inf check on solver output
#ifndef NDEBUG
for (int i = 0; i < X_.Size(); i++) {
   MFEM_VERIFY(std::isfinite(X_(i)),
      "Domain solve produced NaN/Inf at DOF " << i);
}
#endif
```

## Files to Modify

1. **`miniapps/seas/domain/elasticity_operator.hpp`**
   - ComputeTraction interior faces (~lines 1993–2046): Fix BR2 face_int scaling
   - ComputeTraction shared faces (~lines 2228–2268): Same fix
   - Solve (~line 1753): Add NaN check

2. **`miniapps/seas/tests/verification/bp5_verification_full.cpp`**
   - Add `--dg-method` CLI argument

3. **`miniapps/seas/jobs/bp5/bp5_1000m_ip.sbatch`** (NEW)
   - IP method run configuration

4. **`miniapps/seas/tests/unit/test_elasticity_operator.cpp`**
   - Update BR2 traction tests to verify new scaling

## Verification

1. Run existing unit tests: `make seas_test_elasticity_operator && ./seas_test_elasticity_operator`
2. Submit both BP5 jobs on TACC:
   - `sbatch jobs/bp5/bp5_1000m.sbatch` — BR2 with fixed scaling
   - `sbatch jobs/bp5/bp5_1000m_ip.sbatch` — IP method
3. Compare V_max time series — both should show earthquake cycle behavior (V peaks, then decreases)
4. Verify patch test still passes (constant strain → exact traction recovery)

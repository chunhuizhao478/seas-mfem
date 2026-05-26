# Code Review: Phase R exact bi-material Riemann solver plan (2026-05-18, rev-4 — corrected)

## REVISION NOTE

This review supersedes the rev-3 draft (same file). Rev-3 carried R-001 as
CRITICAL on the assertion that the plan's
`Q* = (I − godunov) · Q_L + godunov · Q_R` does NOT collapse to standard
upwind in the homogeneous limit. **That finding was wrong**, based on a
sign error on my part when checking the eigenvalues of SeisSol's matR
columns. After re-verifying against SeisSol's actual code:

1. SeisSol's `getTransposedCoefficientMatrix` returns `Aᵀ` (transpose),
   not `A`. So `starMatrix(0)` in the yateto codegen is `Aᵀ`.
2. SeisSol's matR col 0 (built with **local** material and **positive**
   sqrt for `v_x`) corresponds to eigenvalue **−c_p^L** (LEFT-going from
   L's POV using L material), NOT +c_p^L as the plan's eigenvalue table
   labels. Verified analytically:

   ```
   A^x · (λ+2μ, λ, λ, 0, …, 0, +c_p, 0, 0)ᵀ
     = (−(λ+2μ)c_p, −λc_p, −λc_p, 0, …, 0, −(λ+2μ)/ρ, 0, 0)ᵀ
     = −c_p · (λ+2μ, λ, λ, 0, …, 0, +c_p, 0, 0)ᵀ        ⟹  eigvalue = −c_p
   ```

   Symmetrically, matR col 8 (built with neighbor material, **negative**
   sqrt for `v_x`) is the +c_p^R eigenvector under `A^R`.

3. Therefore chi = diag(1,1,1,0,…,0) selects the LEFT-going (negative-
   eigenvalue) subspace, godunov projects onto LEFT-going, and:

   ```
   A · godunov         = A⁻
   A · (I − godunov)   = A − A⁻ = A⁺
   F = A · Q* = A⁺ Q_L + A⁻ Q_R          ✓ STANDARD UPWIND
   ```

   2×2 acoustic toy problem (`A=[[0,−K],[−1/ρ,0]]`, c=√(K/ρ), Z=ρc):

   ```
   matR    = [[K, K], [+c, −c]]
   godunov = (1/2) · [[1, Z], [1/Z, 1]]
   A · godunov          = (1/2) · [[−c, −K], [−1/ρ, −c]] = A⁻     ✓
   A · (I − godunov)    = (1/2) · [[+c, −K], [−1/ρ, +c]] = A⁺     ✓
   ```

The plan's CORE FORMULA is correct. The bugs are in the eigenvalue
**labels** at PLAN lines 99-110 and in the muddled "homogeneous-limit
collapse" wording at lines 199-207. Both are MODERATE (documentation /
specification) and would not have caused the implementing agent to
produce wrong code if they copy the SeisSol matR construction verbatim.

R-001 is reformulated below as a MODERATE labeling bug; R-002 through
R-008 from the prior review are unchanged.

## Review Scope
- Plan: `miniapps/seas/safs/project_7.0_alternative/document/PLAN_phase_R_exact_bimaterial_riemann.md`
  (1,141 lines, 58 sections; rev-2 against SeisSol, then rev-3 cadence reordering)
- Reference cross-checked: SeisSol
  (`/Users/chunhuizhao/projects/SeisSol/src/Equations/elastic/Model/ElasticSetup.h:28-166` —
   `getTransposedCoefficientMatrix` and `getTransposedGodunovState`;
   `/Users/chunhuizhao/projects/SeisSol/src/Initializer/CellLocalMatrices.cpp:200-355` —
   `initializeCellMatrices`;
   `/Users/chunhuizhao/projects/SeisSol/codegen/kernels/aderdg.py:227-246, 368-385` —
   `computeFluxSolverLocal` and `localFlux` runtime kernel;
   `/Users/chunhuizhao/projects/SeisSol/src/Kernels/LinearCK/Local.cpp:130-180` —
   `nApNm1` application at runtime).
- Domain context: `miniapps/seas/CLAUDE.md` (sign conventions, byte-exact contract,
  "extreme care" file list);
  `miniapps/seas/dynamic/godunov_flux.{hpp,cpp}` (existing homogeneous flux to be
  preserved verbatim);
  `~/.claude/projects/-Users-chunhuizhao-projects-seas-mfem/memory/feedback_complete_sign_sites.md`.

## Findings

### [R-001] [MODERATE] [PLAN.md:99-110 eigenvalue table + PLAN.md:143-146 chi description] — Eigenvalue table has signs swapped; chi description calls left-going waves "local-outgoing"

**Category:** BUG (documentation; labels contradict the actual matR construction)

**Description:**
The plan's eigenvalue table at lines 99-110 reads (paraphrased):

```
| col | eigenvalue   | description                                | material |
|-----|--------------|--------------------------------------------|----------|
| 0   | +c_p^L       | local right-going P-wave                   | local    |
| 1   | +c_s^L       | local right-going S-wave, y-polarised      | local    |
| 2   | +c_s^L       | local right-going S-wave, z-polarised      | local    |
| 3   | 0            | (zero mode)                                | local    |
| 4   | 0            | (zero mode)                                | local    |
| 5   | 0            | (zero mode)                                | local    |
| 6   | −c_s^R       | neighbour left-going S-wave, z-polarised   | neighbour|
| 7   | −c_s^R       | neighbour left-going S-wave, y-polarised   | neighbour|
| 8   | −c_p^R       | neighbour left-going P-wave                | neighbour|
```

The actual matR construction (which the plan correctly reproduces from
SeisSol verbatim at lines 113-141, with cols 0/1/2 using `+sqrt(c)` for
the velocity entry and L material; cols 6/7/8 using `−sqrt(c)` for the
velocity entry and R material) yields the OPPOSITE eigenvalue signs:

- Cols 0, 1, 2: eigenvalues **−c_p^L, −c_s^L, −c_s^L** (LEFT-going, L material)
- Cols 6, 7, 8: eigenvalues **+c_s^R, +c_s^R, +c_p^R** (RIGHT-going, R material)

Verified analytically (see REVISION NOTE above and 2×2 toy problem).

In parallel, the plan's chi-description at lines 143-146:
> chi = diag(1, 1, 1, 0, 0, 0, 0, 0, 0)
> "selects the local-outgoing characteristic family (columns 0,1,2 of matR)"

is misleading: cols 0,1,2 are NOT "local-outgoing" in the standard sense
(waves leaving the L cell). They are LEFT-going waves with respect to the
face normal, i.e., waves emanating from the Riemann fan into L's interior.
In LeVeque §22 / Pelties (2012) terminology these are the "negative-speed
waves" or "left-going waves" — the ones we sum into Q* in
`Q* = Q_L + Σ_{s^p < 0} α^p r^p`.

**Trigger:** an implementer or reviewer who tries to verify the plan's
formula from the eigenvalue table will derive:
- chi picks cols 0,1,2 = "+c eigenvalues" (per table) = RIGHT-going / OUTGOING from L,
- godunov projects onto RIGHT-going (positive-eigenvalue) subspace,
- A · godunov = A⁺ (per the table's labels),
- A · (I − godunov) = A − A⁺ = A⁻,
- F = A · Q* = A⁻ Q_L + A⁺ Q_R — SWAPPED from standard upwind.

They will then conclude that the plan's math is wrong (which I did, in the
rev-3 draft of this review) and spend hours debugging a phantom bug.

**Actual behavior (if implemented per the plan's CODE, not its table):**
The implementer copies the matR construction from SeisSol verbatim,
copies the formula `Q* = (I − godunov) Q_L + godunov Q_R` verbatim, and
the code is CORRECT. R.1.T-1 byte-exact gate PASSES on first run.

**Actual behavior (if implemented per the plan's TABLE):**
A confused implementer who tries to re-derive matR from the eigenvalue
labels (e.g., uses `matR(6, 0) = -sqrt(...)` to force a "+c_p" eigenvector)
gets a matR with the WRONG columns. chi then selects the wrong subspace.
Standard upwind is reversed at the homogeneous limit. R.1.T-1 FAILS.

**Expected behavior:** the table labels match the actual eigenvalues of
the matR construction, and chi's role is described accurately.

**Suggested fix:**

```diff
@@ Eigenvalue / eigenvector table @@
-| col | eigenvalue | description                                | material |
-|-----|------------|--------------------------------------------|----------|
-| 0   | +c_p^L     | local right-going P-wave                   | local    |
-| 1   | +c_s^L     | local right-going S-wave, y-polarised      | local    |
-| 2   | +c_s^L     | local right-going S-wave, z-polarised      | local    |
-| 3   | 0          | (zero mode)                                | local    |
-| 4   | 0          | (zero mode)                                | local    |
-| 5   | 0          | (zero mode)                                | local    |
-| 6   | −c_s^R     | neighbour left-going S-wave, z-polarised   | neighbour|
-| 7   | −c_s^R     | neighbour left-going S-wave, y-polarised   | neighbour|
-| 8   | −c_p^R     | neighbour left-going P-wave                | neighbour|
+| col | eigenvalue | description                                | material |
+|-----|------------|--------------------------------------------|----------|
+| 0   | −c_p^L     | left-going P-wave (into L's interior)      | local    |
+| 1   | −c_s^L     | left-going S-wave, y-pol (into L)          | local    |
+| 2   | −c_s^L     | left-going S-wave, z-pol (into L)          | local    |
+| 3   | 0          | (zero mode)                                | local    |
+| 4   | 0          | (zero mode)                                | local    |
+| 5   | 0          | (zero mode)                                | local    |
+| 6   | +c_s^R     | right-going S-wave, z-pol (into R)         | neighbour|
+| 7   | +c_s^R     | right-going S-wave, y-pol (into R)         | neighbour|
+| 8   | +c_p^R     | right-going P-wave (into R's interior)     | neighbour|
+
+**Eigenvalue sign convention** — matR cols 0,1,2 have NEGATIVE
+eigenvalues despite using positive-`+sqrt` for the velocity entry; this
+is because the eigenvector structure
+`R = (λ+2μ, λ, λ, 0, …, 0, +c, 0, 0)ᵀ` satisfies `A·R = −c·R`, verified
+analytically (see REVIEW_phase_R_exact_bimaterial_riemann.md REVISION
+NOTE). Symmetrically, cols 6,7,8 with `−sqrt` velocity give POSITIVE
+eigenvalues under neighbor material. A standalone unit test
+(R-001 test below) MUST verify these eigenvalues numerically before
+anything downstream is built.
```

And at the chi-description (lines 143-146):

```diff
-    chi = diag(1, 1, 1, 0, 0, 0, 0, 0, 0)        ∈ ℝ^{9×9}
-
-This selects the local-outgoing characteristic family (columns 0,1,2
-of matR).
+    chi = diag(1, 1, 1, 0, 0, 0, 0, 0, 0)        ∈ ℝ^{9×9}
+
+This selects the LEFT-going characteristic family (cols 0,1,2 of
+matR), i.e., waves emanating from the Riemann fan into the L cell's
+interior.  In LeVeque (22.21) terms these are the "negative-speed
+waves" or "left-going waves" — the ones we sum into Q*:
+    Q* = Q_L + Σ_{s^p < 0} α^p r^p
+with `α = matR⁻¹ · (Q_R − Q_L)`.  Equivalently:
+    Q* = (I − godunov) · Q_L + godunov · Q_R
+where `godunov = matR · chi · matR⁻¹` is the projector onto the LEFT-
+going subspace.  In the homogeneous limit, A · godunov = A⁻, so
+F = A · Q* = A⁺ Q_L + A⁻ Q_R — the standard upwind flux.
```

**Test case:**

```cpp
// tests/unit/test_phaser_matr_eigenvalues.cpp (NEW, lands BEFORE R.1
// implementation begins to anchor the sign convention)
TEST(PhaseRMatR, R001_col_0_is_negative_cP_eigenvector_of_local_A)
{
    const real_t lambda = 2.0, mu = 1.0, rho = 1.0;
    const real_t lambda2mu = lambda + 2.0 * mu;
    const real_t cP = std::sqrt(lambda2mu / rho);

    // Build SeisSol's matR col 0 (with local material, +sqrt for v_x):
    mfem::Vector R(9); R = 0.0;
    R(0) = lambda2mu;
    R(1) = lambda;
    R(2) = lambda;
    R(6) = +cP;          // SeisSol convention: positive sqrt

    // Build A^x (local material, OUR convention NOT SeisSol's transpose):
    mfem::DenseMatrix A(9, 9); A = 0.0;
    A(0, 6) = -lambda2mu;
    A(1, 6) = -lambda;
    A(2, 6) = -lambda;
    A(3, 7) = -mu;
    A(5, 8) = -mu;
    A(6, 0) = -1.0 / rho;
    A(7, 3) = -1.0 / rho;
    A(8, 5) = -1.0 / rho;

    mfem::Vector AR(9);
    A.Mult(R, AR);

    // Verify A · R = -cP · R (NEGATIVE eigenvalue, NOT positive):
    for (int i = 0; i < 9; ++i)
    {
        EXPECT_NEAR(AR(i), -cP * R(i), 1e-12)
            << "matR col 0 eigenvalue mismatch at component " << i
            << " — expected -cP*R(i) = " << -cP * R(i)
            << ", got A*R(i) = " << AR(i)
            << ".  The plan's eigenvalue TABLE (PLAN.md:99-110) claims"
            << " col 0 has eigenvalue +cP, but the actual matR"
            << " construction (PLAN.md:113-141 + verified against SeisSol"
            << " ElasticSetup.h:101-105) gives -cP.";
    }
}

TEST(PhaseRMatR, R001_col_8_is_positive_cP_eigenvector_of_neighbor_A)
{
    // Symmetric to above: matR col 8 (neighbor material, -sqrt for v_x)
    // is the +cP_R eigenvector of A^R, NOT -cP_R as the table claims.
    const real_t lambda_R = 4.0, mu_R = 2.0, rho_R = 1.0;
    const real_t lambda2mu_R = lambda_R + 2.0 * mu_R;
    const real_t cP_R = std::sqrt(lambda2mu_R / rho_R);

    mfem::Vector R(9); R = 0.0;
    R(0) = lambda2mu_R;  R(1) = lambda_R;  R(2) = lambda_R;
    R(6) = -cP_R;        // SeisSol convention: NEGATIVE sqrt for cols 6/7/8

    mfem::DenseMatrix A_R(9, 9);  /* fill with neighbor material */
    /* ... */

    mfem::Vector AR(9); A_R.Mult(R, AR);
    for (int i = 0; i < 9; ++i)
        EXPECT_NEAR(AR(i), +cP_R * R(i), 1e-12);
}

TEST(PhaseRMatR, R001_homogeneous_limit_byte_exact_vs_GodunovFlux_Interior)
{
    // For homogeneous (lambda_L,mu_L,rho_L) = (lambda_R,mu_R,rho_R),
    // BimaterialFlux::ComputeFluxFaceLocal(...) must produce
    // bit-identical output to GodunovFlux::Interior(...) for any random
    // (Q_L, Q_R, normal).  This is the headline R.1.T-1 gate; running
    // it FIRST as a sanity check anchors the math sign convention.
    std::mt19937 rng(20260518);
    std::uniform_real_distribution<real_t> u(-1.0, 1.0);
    for (int trial = 0; trial < 100; ++trial)
    {
        const real_t lambda = 1.0 + u(rng), mu = 1.0 + u(rng), rho = 1.0 + u(rng);
        mfem::Vector Q_L(9), Q_R(9), nor(3);
        for (int i = 0; i < 9; ++i) { Q_L(i) = u(rng); Q_R(i) = u(rng); }
        nor(0) = u(rng); nor(1) = u(rng); nor(2) = u(rng);
        nor /= nor.Norml2();

        mfem::Vector F_bimat(9), F_homog(9);
        BimaterialFlux::ComputeFluxFaceLocal(
            lambda, mu, rho, lambda, mu, rho, Q_L, Q_R, nor, F_bimat);
        GodunovFlux flux(lambda, mu, rho);
        flux.Interior(Q_L, Q_R, nor, F_homog);

        for (int i = 0; i < 9; ++i)
            EXPECT_DOUBLE_EQ(F_bimat(i), F_homog(i))
                << "BIT-EXACT mismatch at trial " << trial << " comp " << i;
    }
}
```

---

### [R-002] [MODERATE] [PLAN.md:199-207 "Homogeneous-limit collapse"] — Wording conflates qGodLocal (projector) with per-side flux matrix (AplusT); implementer may drop the A multiplication

**Category:** BUG (specification — incorrect mathematical claim)

**Description:**
The plan's "Homogeneous-limit collapse" paragraph (lines 199-207) says:
> When `(λ_L, μ_L, ρ_L) == (λ_R, μ_R, ρ_R)`, the first 3 columns of matR
> equal the last 3 columns (up to sign), `godunov = matR · chi · matR⁻¹`
> collapses to the homogeneous upwind projector, and
> `qGodLocal · Q_L + qGodNeighbor · Q_R` reproduces the existing
> `GodunovFlux::Interior(...)` formula `A⁺ · Q_L + A⁻ · Q_R`
> bit-identically.

The claim `qGodLocal · Q_L + qGodNeighbor · Q_R = A⁺ Q_L + A⁻ Q_R` is
WRONG. qGodLocal is a 9×9 PROJECTOR (= `I − godunov^T` in SeisSol's
transposed form, or `I − godunov` in the row-vector form); it has no
Jacobian A inside it. So `qGodLocal · Q_L` is just a projected state
vector, not a flux. The actual flux from L's side is

```
F_local_facelocal = A · Q* = A · (I − godunov) · Q_L + A · godunov · Q_R
                          = A · qGodLocal · Q_L + A · qGodNeighbor · Q_R
                          ≠  qGodLocal · Q_L + qGodNeighbor · Q_R
```

SeisSol pre-folds the A multiplication into the precomputed `AplusT`
matrix at init time (CellLocalMatrices.cpp:255-265, codegen
aderdg.py:227-246):

```
AplusT_codegen = fluxScale · Tinv^T · QgodLocal · A^T · T^T
                              ^^^^^^^^^^^^^^^^^^^^^^^^
                              ↑ A is folded in HERE
```

The plan's section on `ComputeFluxFaceLocal` (lines 156-180) DOES
correctly write `F_local = A · Q*` (line 175). So the implementer's
runtime code path will be right.

But the "Homogeneous-limit collapse" paragraph (the section that the
implementer reads when checking their byte-exact gate) tells them they
should be comparing `qGodLocal · Q_L + qGodNeighbor · Q_R` against
`GodunovFlux::Interior`. They'll write the test that way, get junk
(because qGodLocal is a projector not a flux), and conclude the gate
is broken.

**Trigger:** the implementer writes R.1.T-1 byte-exact test referencing
the formula in the "Homogeneous-limit collapse" paragraph literally.

**Actual behavior (if test follows plan's wording):** the byte-exact
test computes `qGodLocal · Q_L + qGodNeighbor · Q_R` and compares it to
`GodunovFlux::Interior`'s flux output. The first is a projected state
(units: stress / velocity); the second is a flux (units: stress / time
/ velocity / time). Dimensional mismatch → test fails immediately or
with absurd numbers.

**Expected behavior:** the paragraph correctly identifies that the
flux comparison is `F = A · Q*` (or equivalently the SeisSol-style
precomputed flux `AplusT · Q_L + AminusT · Q_R`), not `qGodLocal · Q_L
+ qGodNeighbor · Q_R`.

**Suggested fix:**
```diff
-When `(λ_L, μ_L, ρ_L) == (λ_R, μ_R, ρ_R)`, the first 3 columns of matR
-equal the last 3 columns (up to sign), `godunov = matR · chi · matR⁻¹`
-collapses to the homogeneous upwind projector, and
-`qGodLocal · Q_L + qGodNeighbor · Q_R` reproduces the existing
-`GodunovFlux::Interior(...)` formula `A⁺ · Q_L + A⁻ · Q_R`
-bit-identically.
+When `(λ_L, μ_L, ρ_L) == (λ_R, μ_R, ρ_R)`, the matR construction
+yields a single-material eigenbasis: cols 0,1,2 are the LEFT-going
+(−c_p, −c_s, −c_s) eigenvectors and cols 6,7,8 are the RIGHT-going
+(+c_s, +c_s, +c_p) eigenvectors of the SAME A.  chi = diag(1,1,1,0,…,0)
+selects the LEFT-going subspace, so `godunov = matR · chi · matR⁻¹`
+is the projector onto the negative-eigenvalue subspace; in the
+homogeneous limit this gives `A · godunov = A⁻` (and therefore
+`A · (I − godunov) = A⁺`).
+
+The flux from L's side is `F_local = A · Q*` (per the
+`ComputeFluxFaceLocal` formula above).  Substituting:
+
+    F_local = A · [(I − godunov) Q_L + godunov Q_R]
+           = A · (I − godunov) · Q_L + A · godunov · Q_R
+           = A⁺ · Q_L + A⁻ · Q_R              ← STANDARD UPWIND
+
+which is bit-identical to `GodunovFlux::Interior(...)`'s output for
+the same `(Q_L, Q_R, normal, lambda, mu, rho)`.
+
+**Byte-exact gate must compare F_local (= A · Q*), NOT the bare
+projector `qGodLocal · Q_L + qGodNeighbor · Q_R`.**  qGodLocal alone
+is a 9×9 PROJECTOR (no Jacobian inside); it produces a STATE, not a
+FLUX.  SeisSol pre-folds the A multiplication into the precomputed
+per-face flux matrix `AplusT = T · A · qGodLocal · Tinv` at init time
+(CellLocalMatrices.cpp:255-265); our implementation does the same
+via `BuildPerFaceFluxMatricesGlobal`.
```

**Test case:**

```cpp
TEST(PhaseRFlux, R002_qGodLocal_alone_is_NOT_a_flux_dimensional_check)
{
    // Sanity test: qGodLocal · Q_L is dimensionally a STATE (units
    // match Q), not a FLUX (units match A · Q).  This test catches a
    // confused implementer who thinks qGodLocal IS the flux operator.
    const real_t lambda = 1.0, mu = 1.0, rho = 1.0;
    mfem::DenseMatrix qGodLocal(9, 9);
    BimaterialFlux::BuildQGodLocalFaceLocal_(
        lambda, mu, rho, lambda, mu, rho, qGodLocal);
    mfem::Vector Q_L(9);
    Q_L = 1.0;  // Unit input.
    mfem::Vector projected(9);
    qGodLocal.Mult(Q_L, projected);
    // For a unit state input on the homogeneous system, projected
    // should be a CONVEX COMBINATION (or related operation) on Q_L,
    // NOT a velocity-scaled flux of magnitude ~c.  A real flux for
    // a unit stress would be ~c_p times O(1), i.e., O(1000) for
    // crustal materials.  qGodLocal · Q_L should be O(1).
    for (int i = 0; i < 9; ++i)
        EXPECT_LT(std::abs(projected(i)), 10.0)
            << "qGodLocal · Q_L magnitude at comp " << i
            << " is " << projected(i)
            << "; expected O(1).  If this is O(c)~1000, the matrix"
            << " is being built as a FLUX instead of a PROJECTOR.";
}
```

---

### [R-003] [MODERATE] [plan §"Numerical / performance constraints" lines 285-290] — Memory budget arithmetic mixes global and per-rank quantities

**Category:** BUG (arithmetic error / specification ambiguity)

**Description:**
The constraint section claims:
> precompute `qGodLocal` (9×9 = 162 doubles) + `qGodNeighbor` (162 doubles) per
> HETEROGENEOUS face only. For SAFS 1000m mesh (~4M tets, ~16M faces, ≲ 5 %
> heterogeneous), this is ≲ 800 K faces × 324 doubles × 8 bytes = ≲ 2 GB per rank.

The arithmetic `800 K × 324 × 8 = 2.07 GB` correctly computes the GLOBAL memory
across all ranks combined. The "per rank" qualifier is incorrect: with `nprocs = 100`
the per-rank memory is ~20 MB, not 2 GB. Furthermore, the doubling note at line 777
("Memory cost doubles to 4 × 81 = 324 doubles per heterogeneous face") would have
us at 324 doubles ALREADY in the headline arithmetic, not 162 + 162. So the
"≲ 2 GB" number is the correct GLOBAL number for the doubled-POV storage at 5 %
heterogeneity.

Then R.3.T-5 (line 823) sets the budget as `< 4 GB per rank`. This is internally
inconsistent with the "≲ 2 GB per rank" headline claim — if 2 GB is the GLOBAL
budget, 4 GB per rank is 200× too generous; if 2 GB is per-rank, 4 GB per rank
is only 2× too generous.

The actual hazard: an SAFS production run with high heterogeneity (e.g., 30 %) at
`nprocs = 8` (small parallel run on a workstation) would consume
`~16 M × 4 × 0.3 × 324 × 8 = ~50 GB` GLOBALLY → `~6 GB per rank` at np=8. The
plan's R.3.T-5 budget would PASS the per-rank check while busting the
workstation's total RAM. The user would not get a warning.

**Trigger:** any heterogeneity ratio > 5 % combined with parallel partition
where per-rank memory budget is more permissive than the workstation's RAM
allowance.

**Actual behavior (if implemented per plan):** R.3.T-5 is satisfied; the
WaveOperator OOM-aborts at construction without a clear warning.

**Expected behavior:** the budget guard reports both the GLOBAL total AND
the per-rank maximum, and fails when EITHER exceeds a configurable bound.

**Suggested fix:**
```diff
-- **Memory bound**: precompute `qGodLocal` (9×9 = 162 doubles) +
-  `qGodNeighbor` (162 doubles) per HETEROGENEOUS face only.  For
-  SAFS 1000m mesh (~4M tets, ~16M faces, ≲ 5 % heterogeneous), this
-  is ≲ 800 K faces × 324 doubles × 8 bytes = ≲ 2 GB per rank.
-  Acceptable.  If the SAFS heterogeneity exceeds 20 % of faces, fall
-  back to per-QP runtime computation (acceptance gate R.3.T-5 below
-  catches this).  **Phase R MUST log the heterogeneous-face count
-  and the resulting per-face memory at WaveOperator ctor time**
-  so the budget can be monitored.
+- **Memory bound**: precompute `qGodLocal` (9×9 = 81 doubles) +
+  `qGodNeighbor` (81 doubles) per HETEROGENEOUS face per POV.  With
+  two POVs per face (see line 777), this is 324 doubles per
+  heterogeneous face.  For SAFS 1000m mesh (~4M tets, ~16M GLOBAL
+  interior faces, ≲ 5 % heterogeneous), TOTAL across all ranks is
+  ~800 K heterogeneous faces × 324 doubles × 8 bytes ≈ 2 GB GLOBAL,
+  i.e., ~20 MB per rank at np=100, ~250 MB per rank at np=8.
+  Acceptable on both.  **Phase R MUST log BOTH the global heterogeneous-
+  face count (gathered by `MPI_Allreduce(SUM)`) AND the per-rank
+  maximum (gathered by `MPI_Allreduce(MAX)`) at WaveOperator ctor
+  time** so the user sees both numbers and either can be the budget
+  trigger.
```

And update R.3.T-5:
```diff
-- [ ] **R.3.T-5 (memory budget)**: on the SAFS 1000m
-  `_lcfar3000.msh` + CVM-H sidecar, rank-0 log reports
-  heterogeneous-face count < 5 % of total interior faces AND
-  `per_face_bimaterial_flux_` total < 4 GB per rank.
+- [ ] **R.3.T-5 (memory budget)**: on the SAFS 1000m
+  `_lcfar3000.msh` + CVM-H sidecar, rank-0 log reports
+  heterogeneous-face count < 5 % of total interior faces AND
+  `per_face_bimaterial_flux_` GLOBAL total (`MPI_Allreduce(SUM)`)
+  < 4 GB AND per-rank maximum (`MPI_Allreduce(MAX)`) < 1 GB.
+  Both bounds must hold; the lower of the two is the binding gate.
```

**Test case:**
```cpp
TEST(PhaseRMemory, R003_per_rank_and_global_budget_both_reported)
{
   // Build a 2-rank partition with HIGHLY heterogeneous material
   // (every face heterogeneous), e.g. by giving each rank's local
   // elements a distinct (lambda, mu, rho).  After WaveOperator ctor,
   // rank 0 log must contain BOTH:
   //   - "global per_face_bimaterial_flux_ total: <X> bytes"
   //   - "per-rank max:                           <Y> bytes"
   // And the test ASSERTs Y > X / nprocs (per-rank max is at least
   // average; equal-partition lower bound).
   //
   // This catches the case where only the global number is reported
   // and a high-imbalance partition busts a single rank's RAM.
}
```

---

### [R-004] [MODERATE] [plan §Phase R.3 acceptance R.3.T-3] — Dispatch-counter test references R.4's mesh, but R.3 lands BEFORE R.4 per the commit cadence

**Category:** DEVIATION (broken dependency ordering)

**Description:**
R.3.T-3 (plan line 813-814) says:
> with `SEAS_DIAG_PHASER` build, on the R.4 two-layer mesh (`nz = 64`,
> `Z_R/Z_L = 2.0`), the count of `ApplyPerFaceFlux` calls per `Mult`
> equals the count of interface faces

But the cadence section (lines 41-74) places R.3 in Commit #2 and R.4 in Commit #3,
with the explicit rule "Commit #3 (R.4) cannot land until Commit #2 is green". R.3
therefore cannot test against an R.4 artifact that does not yet exist.

The Dependencies block of R.3 (line 832-834) partially acknowledges this:
> Depends on: Phase R.1, Phase R.2 (the two-layer fixture is used
> by R.3.T-3's dispatch counter check)

— pointing at R.2's `TwoLayerMaterial` fixture. But R.2 produces only the
Coefficient, not the mesh. The mesh is built inside `test_phaser_layered_p_wave.cpp`
(R.4's deliverable). Without a mesh, R.3.T-3 cannot count dispatch calls.

**Trigger:** implementer reaches Commit #2 (R.3), cannot run R.3.T-3 because no
heterogeneous mesh fixture exists yet.

**Actual behavior (if implemented per plan):** R.3.T-3 cannot be satisfied at
Commit #2 time. Implementer either (a) waits for R.4, breaking the cadence guard,
(b) creates an ad-hoc inline mesh inside R.3's commit and never reuses it, (c)
skips R.3.T-3.

**Expected behavior:** R.3.T-3 reuses ONLY R.1 and R.2 artifacts (both already
landed by Commit #1).

**Suggested fix:** Build the test mesh inside R.3's own test file (analogous to
the existing `test_phaseh_wave_operator_constant_parity.cpp` pattern, which builds
its own inline `MakeCartesian3D` mesh). The R.4 test then BUILDS ON the R.3
infrastructure for its own purposes.

```diff
-- [ ] **R.3.T-3 (heterogeneous dispatch counter)**: with
-  `SEAS_DIAG_PHASER` build, on the R.4 two-layer mesh (`nz = 64`,
-  `Z_R/Z_L = 2.0`), the count of `ApplyPerFaceFlux` calls per `Mult`
-  equals the count of interface faces (= `nx · ny ·
-  nbf_per_face = 16 · nbf_per_face`).  Proof that the dispatch
-  reaches the heterogeneous branch on the test fixture R.4 will use.
+- [ ] **R.3.T-3 (heterogeneous dispatch counter)**: with
+  `SEAS_DIAG_PHASER` build, build an inline `Mesh::MakeCartesian3D(4, 4, 8)`
+  mesh + `TwoLayerMaterial` (R.2 fixture) with `z_interface = 0.5`,
+  `Z_R/Z_L = 2.0`, INSIDE the R.3 test file.  Construct a
+  `WaveOperator<ParMesh>(...)` with the MaterialField from the fixture,
+  then assert the count of `ApplyPerFaceFlux` calls per `Mult` equals
+  the count of interface faces (= `nx · ny · nbf_per_face =
+  16 · nbf_per_face`).  Proof that the dispatch reaches the
+  heterogeneous branch.  R.4 reuses the same fixture + mesh pattern
+  for its full convergence sweep.
```

---

### [R-005] [MODERATE] [plan §Phase R.4 step 5] — Convergence-slope acceptance band excludes the theoretical optimum

**Category:** BUG (acceptance criterion fails on correct code)

**Description:**
R.4 step 5 + R.4.T-3 (lines 886, 915) require:
> log-log fit of L2 error vs `h` over `nz ∈ {32, 64, 128}` returns slope ∈
> `[0.8, 1.5]` (the FE order `p = 1` ± slack for discrete pulse-tracking error).

For DG on hyperbolic problems with polynomial order `p = 1`, the theoretical L2
convergence rate is `O(h^{p+1/2}) = O(h^{1.5})` (Cockburn-Shu 1991, the upper
bound for tightest DG). Many implementations achieve `O(h^{p+1}) = O(h^{2})`
on smooth solutions. The actual measured slope for a Gaussian pulse on a
2-layer medium typically lies in `[1.5, 2.0]`.

The plan's band `[0.8, 1.5]` ADMITS the suboptimal `1.0` rate (which a buggy
implementation might produce) and CAPS at exactly 1.5 — which would put any
slope at or above 1.5 (the THEORETICALLY CORRECT value) AT the boundary,
giving a 50 % chance of failing the gate to floating-point noise.

**Trigger:** a CORRECT implementation that achieves the theoretical
`O(h^{1.5})` or better rate.

**Actual behavior (if implemented per plan):** the acceptance gate may
spuriously fail on the upper bound; OR a buggy implementation with `O(h)`
convergence (instead of `O(h^{1.5})`) passes silently.

**Expected behavior:** the band should encompass the theoretical optimum
with reasonable slack, and EXCLUDE the suboptimal `O(h)` rate so a
first-order-only bug is caught.

**Suggested fix:**
```diff
-6. Convergence sweep: log-log fit L2 error vs `h = 1/nz`; assert
-   slope ∈ `[0.8, 1.5]` for `nz ∈ {32, 64, 128}` (`p = 1` ± slack).
+6. Convergence sweep: log-log fit L2 error vs `h = 1/nz`; assert
+   slope ∈ `[1.3, 2.5]` for `nz ∈ {32, 64, 128}`.  Theory: DG with
+   p=1 gives `O(h^{p+1/2}) = O(h^{1.5})` (Cockburn-Shu 1991) to
+   `O(h^{p+1}) = O(h^{2})` on smooth solutions; the lower bound
+   `1.3` admits modest pre-asymptotic deviation, the upper bound
+   `2.5` admits super-convergence on this particular fixture.
+   A measured slope < 1.0 indicates a first-order error (likely the
+   bi-material flux is degraded to a non-upwind first-order scheme);
+   fail in that case so the bug is caught.
```

```diff
-- [ ] **R.4.T-3**: log-log fit of L2 error vs `h = 1/nz` over
-  `nz ∈ {32, 64, 128}` returns slope ∈ `[0.8, 1.5]`.
+- [ ] **R.4.T-3**: log-log fit of L2 error vs `h = 1/nz` over
+  `nz ∈ {32, 64, 128}` returns slope ∈ `[1.3, 2.5]` (theoretical
+  optimum for DG-p1 on hyperbolic problems is `[1.5, 2.0]`).
```

---

### [R-006] [MODERATE] [plan §Phase R.1 step 4 + step 6] — Phantom "index transposition" warning misdirects implementer; misses real zero-mode-column scaling convention

**Category:** ASSUMPTION (specification gap that misdirects the implementer)

**Description:**
Phase R.1 §3 step 4 (line 437-441) warns:
> SeisSol uses a state-vector ordering of `(σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz, u, v, w)`
> that differs in the σ_xy/σ_yz/σ_xz position from our
> `(SXX, SYY, SZZ, SXY, SYZ, SXZ, VX, VY, VZ)` by one transposition.

The two orderings are IDENTICAL for SXY/SYZ/SXZ (both place σ_xy at index 3, σ_yz
at 4, σ_xz at 5). There is no transposition. The plan flags a phantom risk and may
cause the implementer to introduce one (silently swapping σ_yz ↔ σ_xz when porting
SeisSol's matR construction).

The REAL risk the plan misses: SeisSol's matR scales the zero-mode columns by
`(λ + 2μ)` for matrix conditioning at extreme impedance contrasts (see
`ElasticSetup.h:115-117`):
```
matR(4, 3) = local.lambda + 2 * local.mu;
matR(1, 4) = local.lambda + 2 * local.mu;
matR(2, 5) = local.lambda + 2 * local.mu;
```
The PHYSICAL zero-mode eigenvectors are unit vectors `e_yy`, `e_zz`, `e_yz`; SeisSol
stores `(λ + 2μ) · e_yy` etc. for matR conditioning. Naive porting that uses unit
zero-mode columns gives a DIFFERENT godunov matrix (off by a per-column scalar),
which is mathematically equivalent (the zero-mode subspace is annihilated by
A·(...)) but DIVERGES from SeisSol byte-for-byte. If the implementing agent
benchmarks their matR against SeisSol's `getTransposedGodunovState` output during
debugging, the byte-difference will be a red herring.

**Trigger:** implementer ports SeisSol's matR and (a) follows the plan's warning
to "swap σ_yz ↔ σ_xz", breaking the homogeneous-limit gate; or (b) replaces the
zero-mode scaling with unit vectors, finds matR doesn't match SeisSol's bytes
and spends hours debugging a non-bug.

**Actual behavior (if implemented per plan):**
- Variant (a): R.1.T-1 fails in mysterious ways.
- Variant (b): R.1.T-1 passes, but cross-validation against SeisSol's matrix
  bytes is misleading and debug time is wasted.

**Expected behavior:** the plan documents the actual mapping (no transposition)
AND flags the zero-mode scaling convention.

**Suggested fix:**
```diff
-   4. **WARNING — index ordering**: SeisSol uses a state-vector
-      ordering of `(σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz, u, v, w)`
-      that differs in the σ_xy/σ_yz/σ_xz position from our
-      `(SXX, SYY, SZZ, SXY, SYZ, SXZ, VX, VY, VZ)` by one
-      transposition.  Build `matR` against OUR state ordering, NOT
-      SeisSol's verbatim — the eigenvector content is the same, only
-      the row indices differ.  Cross-check by asserting against
-      `GodunovFlux::GetAx()` in the homogeneous limit (R.1.T-1).
+   4. **Index ordering** — our `(SXX, SYY, SZZ, SXY, SYZ, SXZ, VX, VY, VZ)`
+      matches SeisSol's `(σ_xx, σ_yy, σ_zz, σ_xy, σ_yz, σ_xz, u, v, w)`
+      VERBATIM for the WAVE-CARRYING rows.  No transposition needed
+      between the two.  Cross-check by asserting `GodunovFlux::Interior`
+      byte-equality in the homogeneous limit (R.1.T-1).
+
+      **WARNING — zero-mode column scaling**: SeisSol scales the
+      zero-mode columns by `(λ + 2μ)` for matR conditioning at extreme
+      impedance contrasts (see `ElasticSetup.h:115-117`):
+      ```
+      matR(4, 3) = lambda + 2*mu;   // σ_yz mode column
+      matR(1, 4) = lambda + 2*mu;   // σ_yy mode column
+      matR(2, 5) = lambda + 2*mu;   // σ_zz mode column
+      ```
+      The PHYSICAL zero-mode eigenvectors are unit vectors `e_yy`,
+      `e_zz`, `e_yz`; SeisSol stores `(λ + 2μ) · e_<...>` for
+      conditioning.  Either preserve SeisSol's scaling for byte-
+      exactness with their published matrices, OR use unit
+      eigenvectors (different godunov matrix but mathematically
+      equivalent because the zero-mode subspace is annihilated by
+      `A · (...)`).  Recommendation: USE SeisSol's scaling so a
+      side-by-side matR byte-equality check during debugging is
+      meaningful.
```

---

### [R-007] [LOW] [plan §"Acoustic special case"] — Acoustic threshold left unspecified; SeisSol uses `epsilon()`, plan uses `1e-12`

**Category:** ASSUMPTION (specification gap)

**Description:**
The plan §"Acoustic special case" (line 209-225) says:
> abort on `μ_L < ε` or `μ_R < ε` with a clear "SAFS acoustic regions out of scope"
> message

Then Phase R.1 step 2 line 421 hardcodes:
```cpp
MFEM_VERIFY(mu_self > 1e-12 && mu_nbr > 1e-12, ...);
```

SeisSol's `testIfAcoustic` (`src/Model/Common.h:33`):
```cpp
return std::abs(mu) <= std::numeric_limits<T>::epsilon();
```
For `double`, `epsilon() = 2.22e-16`, four orders of magnitude tighter than the plan's
`1e-12`. Discrepancy: a material with `mu = 1e-13` Pa would pass SeisSol's check
(treated as acoustic) but fail the plan's (rejected as having too-small μ).

This will not affect SAFS production (CVM-H has `μ_min > 1e9 Pa`), but a future
test on a synthetic low-shear material would behave inconsistently between the
two codes — and if the implementer later adds the acoustic branch (currently
out of scope), the threshold mismatch will be a silent bug.

**Suggested fix:**
```diff
-      MFEM_VERIFY(mu_self > 1e-12 && mu_nbr > 1e-12,
+      // Match SeisSol's testIfAcoustic threshold
+      // (src/Model/Common.h:33: |mu| <= epsilon()).  Anything below
+      // this is treated as acoustic and aborted (SAFS-out-of-scope).
+      MFEM_VERIFY(std::abs(mu_self) > std::numeric_limits<real_t>::epsilon()
+                  && std::abs(mu_nbr)  > std::numeric_limits<real_t>::epsilon(),
                   "BimaterialFlux::BuildGodunovStateFaceLocal: acoustic input "
```

---

### [R-008] [LOW] [plan §"Commit cadence" Commit #2 table + R.3.T-1] — Pre-Phase-R baseline references a hardcoded SHA (`b18bed5`)

**Category:** QUALITY (test reproducibility)

**Description:**
R.3.T-1 (line 803-805) and the Risk Assessment §6 (line 1056) both reference
"the pre-Phase-R baseline (Commit `b18bed5`)".

If Commit #2 (R.3 + Phase H Stage 2) lands AFTER any unrelated commit lands on
`feature/safs-quasi-dynamic` (e.g., a separate documentation commit, a bugfix to
an unrelated driver, a CI tweak), the baseline shifts but the SHA stays the same.
The acceptance criterion would then compare against an outdated baseline, possibly
masking a regression that landed in the intervening commits.

**Suggested fix:**
```diff
-- [ ] **R.3.T-1 — HEADLINE GATE: TPV / BP5 byte-exact regression.**
-  `make test-tpv104`, `make test-tpv205`, `make test-tpv102-local`
-  ...
-  bit-identical output to the pre-Phase-R baseline (Commit `b18bed5`).
+- [ ] **R.3.T-1 — HEADLINE GATE: TPV / BP5 byte-exact regression.**
+  `make test-tpv104`, `make test-tpv205`, `make test-tpv102-local`
+  ...
+  bit-identical output to the parent of Phase R Commit #2
+  (`git rev-parse HEAD^` AT THE TIME R.3 lands).  Capture the parent
+  SHA in the commit message as
+  `Pre-Phase-R baseline: <sha>` so the bisect surface is reproducible
+  even after the branch advances.  Do NOT hardcode the SHA in this
+  plan — it can stale.
```

---

### [R-009] [LOW] [plan §R.5.T-3] — Cross-reference acceptance pins on a line number that already shifts when the plan re-renders

**Category:** QUALITY (test fragility)

**Description:**
R.5.T-3 (line 980):
> `spatial_dynamic_rupture_plan.md:586` no longer carries the "future work"
> caveat; it cross-references this plan.

A line-number pin breaks the moment any earlier line is added/removed in
`spatial_dynamic_rupture_plan.md` — including the rev-2 cleanup that already
happened, which may have shifted line 586 to a different position. A future
maintainer expanding the parent plan's overview will shift line numbers and
break this acceptance gate without changing the substantive content the gate
checks.

**Suggested fix:**
```diff
-- [ ] **R.5.T-3**: `spatial_dynamic_rupture_plan.md:586` carries
-  the Phase R cross-reference (verify it matches this file's
-  path).
+- [ ] **R.5.T-3**: `spatial_dynamic_rupture_plan.md` no longer
+  contains the literal string "future-work item shared with that
+  plan" anywhere in the file (the legacy caveat).  Instead, a
+  grep for "PLAN_phase_R_exact_bimaterial_riemann.md" returns
+  AT LEAST ONE match in that file.  Verify via:
+  ```bash
+  ! grep -q 'future-work item shared' \
+         spatial_dynamic_rupture_plan.md \
+  && grep -q 'PLAN_phase_R_exact_bimaterial_riemann' \
+         spatial_dynamic_rupture_plan.md
+  ```
```

---

## Summary

- Critical issues: **0**  (rev-3's R-001 was a false alarm; see REVISION NOTE)
- Moderate issues: **6** (R-001 reformulated, R-002, R-003, R-004, R-005, R-006)
- Low issues: **3** (R-007, R-008, R-009)
- Plan compliance: N/A (this is a plan review, not a code review)
- Verdict: **PASS WITH FIXES — fix R-001 (eigenvalue table) and R-002
  (qGodLocal-vs-flux wording) before any implementation work begins.**

  These two are MODERATE-but-load-bearing: the plan's math FORMULAS are
  correct, but its documentation (eigenvalue table, homogeneous-limit
  paragraph) misrepresents what the formulas mean. An implementer
  following the SeisSol matR construction VERBATIM will produce
  correct code regardless. An implementer trying to verify the math
  from the table or to write the byte-exact gate from the prose will
  get junk and waste cycles. R-001 and R-002 fix the documentation
  so it matches the (correct) math; they do not change any code.

  R-003 through R-006 are correctness / specification issues that
  should be fixed before the implementing agent starts the corresponding
  phase. They do not block the start of R.1, but they do block clean
  completion of R.3, R.4, and the memory-budget guard.

  R-007, R-008, R-009 are quality issues that can land alongside the
  fixes above; they do not block any phase on their own.

  **The PHASE R PLAN IS SOUND** for the core math (Riemann projector,
  per-face precomputation pattern, SeisSol-style matR construction).
  The TPV/BP5 byte-exact regression headline gate (R.3.T-1) and the
  homogeneous-limit unit gate (R.1.T-1) WILL pass on correct
  implementations.

## Unreviewed Areas

- **`BuildPerFaceFluxMatricesGlobal` implementation details (plan lines 462-478)**:
  the plan's three-matrix product `T · A_self_facelocal · qGodLocal · T⁻¹`
  was not verified against SeisSol's `Tinv · QgodLocal · A · T` form on
  the FULL 9×9 elastic system (only on the 2×2 acoustic reduction).
  Per the corrected derivation in this review's REVISION NOTE, the
  two formulations should agree once the convention is clarified
  (SeisSol stores transposed forms for column-major BLAS; our plan
  uses row-major un-transposed forms). The implementing agent should
  ADD a unit test that builds both forms and asserts they agree on
  random inputs.

- **Cross-rank `ExchangeBiMaterialNeighbours_` MPI exchange (Phase H Stage 2)**:
  the plan defers Phase H Stage 2 to a separate sub-phase. The plan
  assumes R.4 lands in the same commit as Phase H Stage 2 but does
  not enumerate what Phase H Stage 2 must deliver (specifically the
  `MPI_Allgatherv` schema for sharing per-rank material parameters
  on shared faces). This is out of scope for THIS plan; flag for
  separate planning round.

- **`material_dedup.hpp` extraction** (plan line 537): the plan calls for
  extracting `GodunovFluxPool::round_sig` into a shared helper. I did not
  audit whether `GodunovFluxPool::Build`'s rounding rule is already
  refactorable without API churn. If it is not, the `IsHomogeneous` threshold
  drift risk (plan Risk Assessment #4) cannot be cleanly fixed via the
  shared-helper approach and a different mitigation is needed.

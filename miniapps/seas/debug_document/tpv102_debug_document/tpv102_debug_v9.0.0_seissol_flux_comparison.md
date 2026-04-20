# SeisSol vs SEAS-MFEM fault-face flux — formula comparison

> Companion to `tpv102_debug_v9.0.0_debug_plan.md`.  Read this **before**
> approving the Option A patch.  The goal is to let you verify that the
> proposed replacement formula (which is what SeisSol uses) is the right
> one, independently of my earlier narrative.
>
> Nothing in the code has been edited.  No tests have been run by the
> debugger yet.  This document summarises source reads only.
>
> Status: investigation complete, patch not yet drafted, waiting for your
> authorisation per the plan's TODO checklist.

---

## 0.  30-second read

SEAS-MFEM and SeisSol build **the same imposed state** at a fault face
(the friction-corrected "Godunov state").  They then feed it into the DG
bulk rhs in **two different ways**:

- SeisSol evaluates `A · T · Q_imposed` **once per side**, on each side
  using **that side's own imposed state**.  The plus element gets the
  plus-side's flux, the minus element gets the minus-side's flux, signs
  set by each element's outward normal.

- SEAS-MFEM stuffs **both imposed states into a Riemann solver**
  (`flux_.Interior(nor, Q_imp_plus_g, Q_imp_minus_g)`) to produce a
  **single** F_h, then assembles that single F_h into **both** elements
  with opposite signs — the standard welded-face DG pattern.

Under pure strike-slip the Riemann call returns exactly zero (algebraic
cancellation).  Even if it didn't, the "one F, two signs" assembly is
wrong at a fault because the two sides have genuinely different fluxes.
The SeisSol form fixes both issues at once.

---

## A.  Verification against the reference papers (added 2026-04-20, round 2)

User asked: does SeisSol actually follow its reference papers, or is our
current approach wrong, or is the proposed fix wrong?  I re-read both
benchmark_document PDFs in full with focus on the DG flux sections.
Verdict:

**SeisSol follows Pelties 2012 literally.  SEAS-MFEM's current flux is
inconsistent with Pelties 2012 — it uses Dumbser-Käser 2006 eq. (14), a
Riemann/upwind flux intended for welded faces, fed with imposed states.
The proposed Option A fix reproduces Pelties 2012 eq. (9)/(20) exactly.**

Primary sources:
- `benchmark_document/An arbitrary high-order discontinuous Galerkin
  method for elastic waves on unstructured meshes – II. The
  three-dimensional isotropic case.pdf`  (Dumbser & Käser 2006, GJI 167)
- `benchmark_document/Journal of Geophysical Research  Solid Earth -
  2012 - Pelties - ...with a.pdf`  (Pelties et al. 2012, JGR 117)

### A.1  Dumbser & Käser 2006 — the Riemann/upwind flux for welded faces

Paper eq. (14), p. 323:

    F_p^h  =  ½ T_pq (A_qr + |A_qr|) (T_rs)^{-1} Q̂_sl^(m)   Φ_l^(m)         ← self
           +  ½ T_pq (A_qr − |A_qr|) (T_rs)^{-1} Q̂_sl^(m_j) Φ_l^(m_j)       ← neighbour

This is **the upwind numerical flux**.  It is what SEAS-MFEM's
`flux_.Interior(nor, Q_self, Q_nbr)` implements literally
(`godunov_flux.cpp:305-349`).  Correct for **welded** faces — no fault.
`A` is the flux Jacobian (eq. (3), p. 321) of the self element using its
own material.  `T / T^{-1}` are the 9×9 rotations for stress+velocity.
`|A| = R |Λ| R^{-1}` is the absolute-value Jacobian from the characteristic
decomposition.

### A.2  Pelties 2012 — the fault-face flux is NOT a Riemann solver

The paper inherits the DG weak form and ADER time integration from
Dumbser-Käser 2006 (Pelties §3, "A detailed description of the adopted DG
scheme is given by Dumbser and Käser [2006]").  It then **replaces** the
flux formula at fault faces.  Key equations, p. 4:

Paper eq. (7)→(8) (DG weak form):

    ∫ Φ_k ∂Q_p/∂t dV  +  Σ_{faces} F_pk  −  ∫ (∂Φ_k/∂x A_pq + …) Q_q dV  =  0

Paper eq. (9) — **the flux across a face with normal aligned to x**:

    F_pk  =  A_pr ∫_t^{t+Δt} ∫_S Φ_k Q̃_r dS dt                             (9)

Paper §[12] verbatim:

> "where Q̃ stands for a suitable approximation of the unknowns on the
> fault and the integral covers the face S and a time interval of size
> Δt … At special boundaries, such as the free surface or faults, the
> values of Q̃ in (9) might be imposed in order to satisfy the physical
> boundary conditions. In the particular case of dynamic faults, we
> impose values derived from the Coulomb friction model (2)."

Paper eq. (9) is **one matrix A, one state Q̃, no split, no
Riemann**.  It is the unsplit flux evaluated at a single interface
state — not the Dumbser-Käser eq. (14) Riemann form.

Paper eq. (13), p. 4 — **the Godunov state** (used for welded / locked
faults):

    2σ_xx^G = (σ_xx^+ + σ_xx^-)  + ρ c_p (u^- − u^+)
    2σ_xy^G = (σ_xy^+ + σ_xy^-)  + (μ/c_s) (v^- − v^+)
    2σ_xz^G = (σ_xz^+ + σ_xz^-)  + (μ/c_s) (w^- − w^+)
    2u^G    = (u^+ + u^-)        + (1/(ρc_p)) (σ_xx^- − σ_xx^+)
    2v^G    = (v^+ + v^-)        + (c_s/μ) (σ_xy^- − σ_xy^+)
    2w^G    = (w^+ + w^-)        + (c_s/μ) (σ_xz^- − σ_xz^+)

This is the standard linear-Riemann exact solution at the interface.
For a welded face, Q̃ = Q^G and eq. (9) gives `F = A · Q^G`.  Classical
result: `½(A+|A|)Q^+ + ½(A−|A|)Q^- ≡ A · Q^G`, so for welded faces
Pelties eq. (9) and Dumbser-Käser eq. (14) are algebraically identical
— this is why SEAS-MFEM's welded-face flux is correct.

### A.3  Pelties 2012 — at an **active-slip** fault, Q̃ is SIDE-SPECIFIC

Paper §[15]: "When the fault is locked, the Godunov state is used in
place of Q̃ in the flux (9).  When active slip is expected, we have to
impose the shear stresses σ_xy and σ_xz on the fault according to the
Coulomb friction model (2) to obtain new traction values σ̃_xy,il and
σ̃_xz,il, which might be different from σ_xy^G and σ_xz^G.  Note that
fault-normal and fault-parallel components are uncoupled in equations
(13).  Because we ignore the possibility of fault opening, the Godunov
state values are assigned to the fault-normal component of stress and
velocity, σ̃_xx,il = σ_xx^G,il and ũ_il = u^G,il."

Paper eq. (15), p. 5 — **imposed velocities differ on each side**:

    ṽ^+_il = v^+_il + (c_s/μ) (σ̃_xy,il − σ_xy^+,il)
    ṽ^-_il = v^-_il − (c_s/μ) (σ̃_xy,il − σ_xy^-,il)
    w̃^+_il = w^+_il + (c_s/μ) (σ̃_xz,il − σ_xz^+,il)
    w̃^-_il = w^-_il − (c_s/μ) (σ̃_xz,il − σ_xz^-,il)

Paper §[16]: "These expressions are crucial for the understanding of
fault dynamics using fluxes, as they state that an imposed shear
traction instantly and locally generates an imposed velocity parallel
to the fault."

Because `c_s/μ = 1/(ρ c_s) = 1/Z_s = invZ_s`, these are **bit-for-bit
the same formulas** SEAS-MFEM computes in
`fault_face_flux.cpp:145-164`.  The imposed-state construction inside
SEAS-MFEM's `Evaluate` IS paper eq. (15).  (Proved by direct comparison
in §4 of this doc.)

Paper eq. (20), p. 5 — **the discrete flux computation**:

    F_pk  =  A_pr  Σ_i Σ_l  ω_i^S ω_l^T  Φ_k(ξ_i)  Q̃_{r,il}                (20)

Paper §[16] preceding eq. (20): "Using the shear stresses σ̃_xy,il and
σ̃_xz,il and the velocities from (15) all values of Q̃ at the interface
are known and the flux (9) can be computed with the discrete expression
[eq. (20)]."

**Eq. (20) is eq. (9) discretised.  One A, one Q̃.  No Riemann solver.**

Because Q̃ differs on the two sides (paper eq. (15): ṽ^+ ≠ ṽ^-), the
flux `F = A · Q̃` is evaluated **separately for each side**, using that
side's Q̃.  This is exactly what SeisSol does (`fluxSolver · QInterpolated`
per-side, `Neighbor.cpp:89-99`, verified in §3 of this doc).

### A.4  Direct verdict on the three possibilities

The user's question: "does SeisSol follow the paper, or is our current
approach wrong, or is the proposed fix wrong?"

1. **Does SeisSol follow the paper?**  YES.  `fluxSolver = A · T`
   (SeisSol `CellLocalMatrices.cpp:801-817`) combined with `I +=
   fluxSolver · QInterpolated` (SeisSol `dynamic_rupture.py:124-130`)
   is the literal discrete form of Pelties eq. (9) / eq. (20).

2. **Is SEAS-MFEM's current approach wrong?**  YES, with respect to
   Pelties 2012.  The code calls `flux_.Interior(nor, Q_imp_plus,
   Q_imp_minus)` which is Dumbser-Käser eq. (14) — the Riemann/upwind
   flux — fed with two different imposed states.  Pelties 2012 does
   **not** prescribe this form at a fault face.  It prescribes eq. (9):
   one A, one Q̃ per side.

3. **Is the proposed Option A fix wrong?**  NO.  Calling
   `flux_.Interior(nor, Q_imp_side, Q_imp_side)` — identical imposed
   state on both Riemann arguments — collapses to `T · (A_x^+ + A_x^-)
   · T^{-1} · Q_imp_side = T · A_x · T^{-1} · Q_imp_side = A · Q̃^side`
   in the global frame.  That is Pelties eq. (9) **verbatim**, just
   routed through SEAS-MFEM's existing `Interior` primitive because the
   identity `A^+ + A^- = A` holds for the characteristic decomposition.
   No new primitive is needed.

### A.5  Why SEAS-MFEM's current form cancels, restated in paper terms

The cancellation identified in plan §1 and this doc's §5 has a paper-
level interpretation.  For a pure-strike-slip fault, Pelties eq. (15)
makes `ṽ^+` and `ṽ^-` differ by exactly `−Δv` (the slip rate) in the
tangential components, while all stress components and the normal
velocity are equal.  If you take the two imposed states `Q̃^+` and `Q̃^-`
and plug them into Dumbser-Käser eq. (14) (i.e., the Riemann flux for
welded faces), you are effectively computing the Godunov state from
two already-friction-corrected half-states; the velocity jump exactly
offsets the tangential stress contribution in the upwind split, and
the flux collapses.  This is a **misuse** of the welded-face upwind
flux at a fault — it assumes the two half-states need a Riemann
projection, but they have already been constructed to BE the left/right
of the Riemann problem's solution.  Pelties 2012 avoids this precisely
by not invoking the Riemann solver: each side's flux is `A · Q̃^side`,
period.

### A.6  One subtlety worth stating

Pelties 2012 eq. (9) has a **single A per face**, i.e., the paper does
not distinguish `A_plus` vs `A_minus`.  The paper considers homogeneous
material at the face.  SeisSol's implementation generalises to
bimaterial by using each element's own star matrix on its own side;
for TPV102 (homogeneous) `A_plus = A_minus`.  SEAS-MFEM, which also uses
a single `GodunovFlux` instance per material, is consistent with the
paper for TPV102 and would need per-side material selection only for
future bimaterial problems.

### A.7  Conclusion for the user

The paper reading changes nothing about the proposed Option A patch —
it grounds it.  Before paper check, the patch was "follow SeisSol."
After paper check, the patch is "follow Pelties 2012 eq. (9) / eq. (20),
which is what SeisSol was built from."  Those are the same fix.

The converse — that SEAS-MFEM's current fault-face flux is a welded-face
upwind flux misapplied at a fault — is now documented as paper-grounded
evidence, not code interpretation.

---

## 1.  Cast of characters

**The fault-face quadrature loop (per-QP):**

```
                +--------------+
  Q_plus, Q_minus  →  |   Evaluate    |  →  Q_imp_plus, Q_imp_minus   (local frame)
  (bulk trace,       +--------------+
   local frame)           Eq. 7–12

  Q_imp_*_local  ─ T_can (local→global) ─→  Q_imp_*_g                (global frame)

                +--------------+
  Q_imp_*_g    →  |  Flux step    |  →  F_h  (global frame)
                +--------------+
                   (← the bug)

  F_h            ─ DG rhs assembly ──→   bulk rhs contribution
                                         (← second layer of bug)
```

**Step 1 "Evaluate" — correct in SEAS-MFEM.**  This is what `fault_face_flux.cpp`
computes: trial traction → friction solve → imposed states that satisfy
Eqs. 10–12.  Plan §A.2 / §3.1 confirms 12/12 PASS at unit level.  Mirrors
SeisSol `FrictionSolverCommon.h::postcomputeImposedStateFromNewStress`
(see §3 of this doc for term-by-term comparison).

**Step 2 "T_can" — correct in SEAS-MFEM.**  Standard rotation; axis-aligned
for TPV102's y=0 fault so the Voigt factor doesn't bite.

**Step 3 "Flux step" — WRONG in SEAS-MFEM, different from SeisSol.**  This
is where §3.1b (the one failing test) fires.

**Step 4 "DG rhs assembly" — also wrong for interior-fault faces.**
Secondary bug; see §6.

---

## 2.  What SEAS-MFEM does today at Step 3

Interior-fault branch, `miniapps/seas/dynamic/wave_operator.inl:811-816`:

```cpp
// Route +/− imposed states to self/nbr slots by elem1_on_plus
const real_t *Q_self_imp = elem1_on_plus ? Q_imp_plus_g  : Q_imp_minus_g;
const real_t *Q_nbr_imp  = elem1_on_plus ? Q_imp_minus_g : Q_imp_plus_g;
real_t F_h_total[NUM_STATE];
flux_.Interior(nor, Q_self_imp, Q_nbr_imp, F_h_total);   // ← ONE flux call, BOTH sides imposed
```

Shared-fault branch, `:1232`:

```cpp
flux_.Interior(can_n, Q_imp_plus_g, Q_imp_minus_g, F_h);  // same shape, same bug
```

`GodunovFlux::Interior` (`miniapps/seas/dynamic/godunov_flux.cpp:325-349`) is
a conventional Godunov / upwind Riemann solver:

```cpp
F = T · (A_x^+ · Q_self_rot + A_x^- · Q_nbr_rot)      // eq. ★
```

where `T` rotates from face-aligned back to global, and `A_x^±` are the
positive/negative characteristic projectors of the x-direction flux
Jacobian.  This is the right Riemann solver **for a welded DG face**.  At
a fault, both arguments come from the friction-corrected imposed state,
and eq. ★ degenerates — see §5.

**The DG rhs assembly for the interior-fault branch (`:857-887`):**

```cpp
rhs[Elem1] -= w · shape1 · F_h_total;    // line 861-862
rhs[Elem2] += w · shape2 · F_h_total;    // line 884-885
```

This "−F on Elem1, +F on Elem2" pattern assumes **one single-valued flux
on the face**, which is valid for a welded face but not for a fault (see
§6).

---

## 3.  What SeisSol does — verified against source

Main call site — `src/Kernels/LinearCK/Neighbor.cpp:89-99`:

```cpp
case FaceType::DynamicRupture: {
   dynamicRupture::kernel::nodalFlux drKrnl = m_drKrnlPrototype;
   drKrnl.fluxSolver    = cellDrMapping[face].fluxSolver;    // side-specific
   drKrnl.QInterpolated = cellDrMapping[face].godunov;       // side's imposed state
   drKrnl.Q             = data.get<LTS::Dofs>();             // self bulk (prefetch only)
   drKrnl.execute(cellDrMapping[face].side,
                  cellDrMapping[face].faceRelation);
}
```

Note: **no Riemann solver is called**.  `drKrnl` is a generated
matrix-vector kernel.

Kernel body — `codegen/kernels/dynamic_rupture.py:124-130`:

```python
I[k,p]  +=  V3mTo2nTWDivM[i,h][k,l]  *  QInterpolated[l,q]  *  fluxSolver[q,p]
```

In standard notation:

    I  +=  M_face_assembly  ·  Q_imposed  ·  fluxSolver

where `I` is the (time-integrated) bulk rhs for **that one element**, and
`M_face_assembly` folds (face→nodal projection) · (quadrature weights)
· (inverse mass).

Precomputed `fluxSolver` — `codegen/kernels/dynamic_rupture.py:87` and
`src/Initializer/CellLocalMatrices.cpp:801-817`:

```
fluxSolver_plus  = (−2 · S_plus  / (6 · V_plus ))  ·  A_plus   ·  T
fluxSolver_minus = (+2 · S_minus / (6 · V_minus))  ·  A_minus  ·  T
```

- `A_plus`, `A_minus` = normal-direction flux Jacobian built from the
  plus / minus element's **own material**.  For bimaterial faults these
  differ; for TPV102 (homogeneous) they are equal.
- `T` = global → face-aligned rotation.
- The scalar factors are SeisSol's ADER-specific geometric scale.  In
  SEAS-MFEM the equivalent scale is injected later via `w = ip.weight ·
  nor_len` in the assembly loop; the matrix part is just `A · T`.
- Opposite sign on plus vs minus encodes the outward normal of each
  element.

**So the contribution each element gets is:**

| element     | rhs update (schematic)                                |
|-------------|-------------------------------------------------------|
| plus-side   | `rhs_plus  += (−1) · (geom) · A_plus  · T · Q_imp_plus`  |
| minus-side  | `rhs_minus += (+1) · (geom) · A_minus · T · Q_imp_minus` |

Each side sees **its own** imposed state.  No L/R averaging.  No Riemann
solver.

---

## 4.  The imposed-state construction is already correct in SEAS-MFEM

Verified term-by-term — the only thing SEAS-MFEM does right here.

| quantity         | SEAS-MFEM `fault_face_flux.cpp:145-164`      | SeisSol `FrictionSolverCommon.h:357-362`  |
|------------------|----------------------------------------------|-------------------------------------------|
| `Q_imp_plus[v_n]`  | `v_n^+ + invZp_p · (σ_n_corr − σ_n^+)`       | `qIPlus[U] + invZp · (normalStress − qIPlus[N])`  |
| `Q_imp_plus[v_t1]` | `v_t1^+ + invZs_p · (τ1_corr − τ1^+)`        | `qIPlus[V] + invZs · (traction1 − qIPlus[T1])`    |
| `Q_imp_plus[v_t2]` | `v_t2^+ + invZs_p · (τ2_corr − τ2^+)`        | `qIPlus[W] + invZs · (traction2 − qIPlus[T2])`    |
| `Q_imp_plus[σ_n]`  | `σ_n_corr`                                   | `normalStress`                            |
| `Q_imp_plus[τ1]`   | `τ1_corr`                                    | `traction1`                               |
| `Q_imp_plus[τ2]`   | `τ2_corr`                                    | `traction2`                               |

Mirror on the minus side (flipped sign on the `invZ` term, from SeisSol
lines 367-372 in the same file; SEAS-MFEM lines 141-147).

Plan §3.1 12/12 PASS is consistent with this table.  **No change to
Evaluate is part of the Option A patch.**

---

## 5.  Why the Riemann call cancels to zero — the algebra in plain code

Setup.  TPV102 pure strike-slip at the post-breakaway fixture
(plan §3 and §A.3):

- canonical frame: `can_n = (0,−1,0)`, `can_t1 = (0,0,−1)` (dip),
  `can_t2 = (+1,0,0)` (strike).
- bulk Q at both sides = 0 (quiescent).
- friction produces `τ2_corr = −η_s · V` (Eq. 10), where V is the slip
  rate in the strike direction and η_s = Z_s/2 for homogeneous sides.
- `τ1_corr = 0`, `σ_n_corr = 0` (only strike channel is active).

Eqs. 11–12 under these assumptions give, in local-frame indices:

```
Q_imp_plus_local[v_t2]  = +V/2               Q_imp_minus_local[v_t2]  = −V/2
Q_imp_plus_local[τ_nt2] = τ2_corr = −η_s V   Q_imp_minus_local[τ_nt2] = τ2_corr = −η_s V
                                           (same value on both sides)
all other components                       = 0
```

Wait — it's actually `−V/2` on plus and `+V/2` on minus in the global
frame after rotation (plan §A.3 records `Q_imp_plus_g[VX] = −V/2`).  The
sign is set by the rotation, not the physics; what matters is:

- **The two sides agree on stress** (`τ2_corr` same on both).
- **The two sides disagree on velocity** by exactly `−V` along strike.

Now plug into the rotated `ApplySplitFlux`
(`godunov_flux.cpp:305-320`), restricted to the S-wave strike-tangential
sub-block (v_t2, τ_nt2):

```
F_rot[v_t2]  = A_x^+[v_t2, τ_nt2] · Q_L_rot[τ_nt2] + A_x^-[v_t2, τ_nt2] · Q_R_rot[τ_nt2]
             +  A_x^+[v_t2, v_t2]  · Q_L_rot[v_t2]  + A_x^-[v_t2, v_t2]  · Q_R_rot[v_t2]
```

With the characteristic entries of `A_x^±` for the S-wave sub-block:

```
A_x^+[v_t2, τ_nt2]   = −1/(2ρ)        A_x^+[v_t2, v_t2] = +c_s/2
A_x^-[v_t2, τ_nt2]   = −1/(2ρ)        A_x^-[v_t2, v_t2] = −c_s/2
```

Plugging Q_L = Q_imp_plus, Q_R = Q_imp_minus:

```
F_rot[v_t2]  = −1/(2ρ) · τ2_corr  − 1/(2ρ) · τ2_corr
             + c_s/2   · (V/2)    − c_s/2  · (−V/2)      // careful with sign of V on each side
             =  −τ2_corr/ρ  +  c_s · V / 2

           with τ2_corr = −η_s·V = −(Z_s/2)·V = −(ρ c_s / 2)·V:
             =  +c_s V / 2  −  c_s V / 2
             =  0          identically
```

**Zero.  Holds ∀ V, ρ, c_s.**  The plan §1 derivation in a slightly
different sign convention gives the same result.  `F_rot[τ_nt2]` cancels
likewise.  So `F = T · F_rot = 0` in the global frame too.

**Root cause in one line:** Eqs. 10–12 force the two imposed states to
satisfy `A_x^+ Q_imp_plus + A_x^- Q_imp_minus ≡ 0` in the tangential
sub-block.  That is the defining property of a friction-corrected
interface state — it is the **single value at the interface**, so the
Riemann projectors on either side agree and cancel.

---

## 6.  SeisSol's form doesn't cancel — why

SeisSol never asks "what is the upwind flux given Q_L, Q_R?"  It says:
"Given that the interface state is Q_imp, evaluate the conservation-form
flux `A · T · Q_imp` in each element's own frame."  Because:

- Plus-side contribution uses **only** Q_imp_plus (not Q_imp_minus),
  so the cancellation in §5 can't happen — the two halves of the
  cancellation come from different flux calls.
- The flux Jacobian A is the **full** normal-direction Jacobian
  (A = A^+ + A^-), not the one-sided projector.

Making this concrete on the strike sub-block: applying `A = A^+ + A^-`
to `Q_imp_plus_local` alone:

```
F_side[v_t2]  = A[v_t2, τ_nt2] · τ2_corr   +   A[v_t2, v_t2] · (+V/2)
             = −1/ρ · τ2_corr              +   0 · (V/2)           // A[v_t2, v_t2] = 0 (full Jacobian diagonal of S-wave sub-block is zero)
             = −τ2_corr / ρ
             = +η_s V / ρ
             = +c_s V / 2
```

Non-zero and physically sensible: O(c_s · V) per side, which for
TPV102 post-breakaway with V = O(1 m/s) and c_s = O(3 km/s) lands in the
O(10³ m²/s²) range.

(Back-of-envelope with the plan's fixture ρ≈2670, τ2_corr=−1.25 MPa:
F ≈ 1.25e6 / 2670 ≈ 468 m²/s², which is inside the plan's §3.1b
acceptance band `∈ (10, 1e4) m²/s²`.)

---

## 7.  The secondary assembly bug (layer 2)

Even if F_h ≠ 0, the interior-fault branch's existing DG rhs pattern
would still be wrong at a fault.  The current code does

```cpp
// wave_operator.inl:857-887
rhs[Elem1] -= w · shape1 · F_h_total;
rhs[Elem2] += w · shape2 · F_h_total;
```

This "−F on Elem1, +F on Elem2" is the standard **conservative** DG face
assembly: flux leaving Elem1 equals flux entering Elem2, because the two
elements share a single-valued interface flux.

At a fault, the flux on the two sides is **not** single-valued.  Plus
and minus elements both see a flux built from their own imposed state.
In SeisSol:

- plus element gets `A · T · Q_imp_plus`  (with plus-outward-normal sign)
- minus element gets `A · T · Q_imp_minus` (with minus-outward-normal sign)

These are **two independent fluxes**, not `±` of each other.  The signs
are independent — they come from each element's outward normal, not from
a conservation identity.

**Shared-fault branch already handles this correctly** via `continue`
after writing only Elem1 (`:1242`).  Each rank writes only its own
Elem1; the other rank writes its own Elem1 (= physical Elem2) with the
opposite sign on a separate flux evaluation.

**Interior-fault branch does NOT** — it writes both Elem1 and Elem2 from
one F_h.  So even in a hypothetical where Step 3 returned a correct
non-zero F_h_total from some averaging scheme, the assembly would still
mis-couple the two sides.

The Option A patch fixes both layers in one stroke by evaluating the flux
per-side and assembling per-side.

---

## 8.  Proposed replacement (pseudocode only — not applied)

### 8.1  Interior-fault branch — `wave_operator.inl:811-887`

Replace lines 811-816:

```cpp
const real_t *Q_self_imp = elem1_on_plus ? Q_imp_plus_g  : Q_imp_minus_g;
const real_t *Q_nbr_imp  = elem1_on_plus ? Q_imp_minus_g : Q_imp_plus_g;
real_t F_h_total[NUM_STATE];
flux_.Interior(nor, Q_self_imp, Q_nbr_imp, F_h_total);
```

with:

```cpp
// Two one-sided A · T · Q fluxes in the canonical-normal frame.
// Using flux_.Interior(n, Q, Q) gives T · (A_x^+ + A_x^-) · Q_rot = T · A · Q_rot,
// which is exactly SeisSol's `fluxSolver · QInterpolated` up to the
// geometric scale that SEAS-MFEM supplies through `w · nor_len` below.
real_t F_h_plus[NUM_STATE], F_h_minus[NUM_STATE];
flux_.Interior(can_n, Q_imp_plus_g,  Q_imp_plus_g,  F_h_plus);
flux_.Interior(can_n, Q_imp_minus_g, Q_imp_minus_g, F_h_minus);
```

Replace the Elem1/Elem2 assembly at lines 857-887 with per-side
assembly:

```cpp
// Elem1 sits on whichever canonical side elem1_on_plus indicates.
// Each side's rhs gets that side's own flux, signed by that side's
// outward normal (+can_n for plus, -can_n for minus).
if (elem1_on_plus) {
   for (int c = 0; c < NUM_STATE; c++) {
      for (int i = 0; i < ndof; i++) {
         rhs[c * ndof_total_ + dof_offset1 + i] -= w * shape1(i) * F_h_plus[c];
         rhs[c * ndof_total_ + dof_offset2 + i] += w * shape2(i) * F_h_minus[c];
      }
   }
} else {
   for (int c = 0; c < NUM_STATE; c++) {
      for (int i = 0; i < ndof; i++) {
         rhs[c * ndof_total_ + dof_offset1 + i] += w * shape1(i) * F_h_minus[c];
         rhs[c * ndof_total_ + dof_offset2 + i] -= w * shape2(i) * F_h_plus[c];
      }
   }
}
```

Using `can_n` (canonical, bit-identical across ranks) instead of MFEM's
rank-local `nor` matches the shared-fault branch's R-802 convention and
removes the need for the L/R swap previously at lines 811-814.

### 8.2  Shared-fault branch — `wave_operator.inl:1232-1244`

Replace:

```cpp
flux_.Interior(can_n, Q_imp_plus_g, Q_imp_minus_g, F_h);
const real_t accum_sign = elem1_on_plus ? +1.0 : -1.0;
for (int c = 0; c < NUM_STATE; c++)
   for (int i = 0; i < ndof; i++)
      rhs[c * ndof_total_ + dof_offset1 + i] -= accum_sign * w * shape1(i) * F_h[c];
```

with:

```cpp
// This rank owns Elem1 only; assemble that side's flux.
real_t F_h_side[NUM_STATE];
const real_t *Q_imp_side = elem1_on_plus ? Q_imp_plus_g : Q_imp_minus_g;
flux_.Interior(can_n, Q_imp_side, Q_imp_side, F_h_side);
const real_t assemble_sign = elem1_on_plus ? -1.0 : +1.0;
for (int c = 0; c < NUM_STATE; c++)
   for (int i = 0; i < ndof; i++)
      rhs[c * ndof_total_ + dof_offset1 + i] += assemble_sign * w * shape1(i) * F_h_side[c];
```

The paired rank (Elem1 on the other side) independently runs the same
block with `elem1_on_plus` flipped, assembling the other side's flux
with the opposite sign.  No ghost writes from either rank.

---

## 9.  Why the existing `flux_.Interior(n, Q, Q)` primitive is sufficient

No new flux routine is needed.  `ApplySplitFlux` at
`godunov_flux.cpp:305-320` already computes

    F_rot  =  A_x^+ · Q_L_rot  +  A_x^- · Q_R_rot

and the eigendecomposition gives `A_x = A_x^+ + A_x^-`, so calling
`Interior(n, Q, Q)` yields

    F  =  T · (A_x^+ + A_x^-) · Q_rot  =  T · A_x · Q_rot  =  A_n · Q

which is exactly SeisSol's `A · T · Q_imposed` operator, rotated to
global.  The only thing the Option A patch changes is **which Q we feed
in** and **how we assemble the result**.

(This is a happy accident of SEAS-MFEM's factored form.  It means the
patch touches only `wave_operator.inl`, not `godunov_flux.cpp`.)

---

## 10.  Assumptions and caveats worth flagging

These are not blockers for R-1001 — TPV102 doesn't exercise them — but
they would bite a future non-TPV102 problem.

1. **Bimaterial faults.**  SeisSol uses `A_plus` ≠ `A_minus` (each side's
   own material).  SEAS-MFEM has a single `GodunovFlux` instance per
   material, so the above works only for homogeneous faults.  If SEAS-MFEM
   ever adds a material jump across a fault, we will need per-side
   `GodunovFlux` selection.  Out of scope for v9.0.0.

2. **Voigt factor in `BuildRotation`.**  `godunov_flux.cpp:238-276`
   (BuildRotation) and `:196-235` (BuildRotationInverse) apply the Voigt
   symmetrisation on different pairs (`(i,j)` vs `(a,b)`).  For
   TPV102's axis-aligned canonical frame the off-diagonal entries never
   mix, so this never fires.  REVIEW already flagged this as
   "Unreviewed" (plan Appendix F); we will not open it in v9.0.0.

3. **Time-integrated vs instantaneous Q.**  SeisSol integrates Q over
   an ADER LTS step and weighs with time-quadrature `weight[o]`
   (FrictionSolverCommon.h:357).  SEAS-MFEM's RK45 uses instantaneous Q;
   the spatial formula is identical, the time-weight loop disappears.

4. **Sign of `F_h[VX_global]`.**  The §3.1b acceptance band `∈ (10,
   1e4) m²/s²` is a magnitude check; the plan does not prescribe a sign.
   Back-of-envelope (§6) predicts `|F_h[VX_global]| ≈ 468 m²/s²` for the
   post-breakaway fixture, comfortably in band.  Sign will be verifiable
   at first run; if it comes out with the opposite sign, it means one of
   the rotation matrices has a handedness we missed — a convention issue,
   not a physics issue — and we'd flip the outermost assembly sign.

---

## 11.  Regression criteria after the fix

From plan §3, reproduced here for easy reference:

- **§3.1b `test_fault_face_flux_frame_and_flux`**: currently `7/9`
  (A1-FAIL on `F_h[VX_global] = 0`).  Expected after fix: `9/9` with
  `|F_h[VX_global]| ∈ (10, 1e4) m²/s²`.
- **§3.1** (`test_fault_face_flux_frame`): must remain 12/12 PASS.
- **§3.1c** (`test_godunov_identity_normal_reversal`): must remain 5/5
  PASS at ULP.  The Godunov identity is a property of `Interior(n, L,
  R)` — unchanged by the patch.
- **§3.1d** (`test_fault_flux_interior_vs_shared_branch_equivalence`):
  currently 4/4 PASS.  The patch changes both branches **in the same
  way** (per-side flux + per-side assembly), so algebraic equivalence
  should continue to hold.  If this test breaks after the fix, it's a
  signal that the interior-branch and shared-branch per-side logic
  diverged — inspect the sign of `assemble_sign` vs `elem1_on_plus`.
- **All other `make test` unit tests**: must remain PASS.

Frontera gates (FP-5 at 1000 m / 100 rank, then Phase 5 at 200 m / 400
rank) come after the unit tests go green — see plan §TODO.

---

## 12.  Commit plan — in execution order, once you approve

Per plan §TODO and `CLAUDE.md` "Proposing Fixes":

1. **Commit A (failing-gate lock):** wire up
   `test_fault_face_flux_frame_and_flux` as it stands, so the failure
   mode F_h[VX_global] = 0 is in git history as a regression gate.
   Nothing else changes in this commit.
2. **Commit B (fix):** apply §8.1 and §8.2 to `wave_operator.inl` only.
   No changes to `fault_face_flux.cpp` or `godunov_flux.cpp`.
3. Re-run `make test`; expect §3.1b to flip 7/9 → 9/9 and all other tests
   to stay green.
4. Update `tpv102_debug_v9.0.0_fix.md` with one entry per step (per plan
   §TODO documentation rule).
5. Update `tpv102_debug_v9.0.0_check.md` with H-V9-J → Closed.
6. Handoff for FP-5 sbatch generation (awaiting your Frontera
   approval per `feedback_frontera_approval.md`).

---

## 14.  Detailed implementation plan (added 2026-04-20, round 3)

Concrete, actionable.  No source edits until you approve §13.

### 14.1  Commit sequence  (revised round 3 per R-F07)

**Default: single-commit strategy — tree stays GREEN at all times.**

| # | Commit scope | Branch state after commit | CI expected |
|---|---|---|---|
| C1 | Source fix + test updates + new unit tests + welded-fallback hard-abort, all together.  Source: `wave_operator.inl:811-887` (main path + fallback) and `:1232-1244`.  Tests: update §3.1b + §3.1d, add §3.1e + §3.1f.  `godunov_flux.hpp` header comment (R-F06).  `FaultFaceFlux::Evaluate` bimaterial `MFEM_VERIFY` (R-F08). | GREEN across `make test` | passes |
| C2 | Documentation: `tpv102_debug_v9.0.0_fix.md` step log entry; `tpv102_debug_v9.0.0_check.md` H-V9-J → Closed; this doc → mark §12/§14 done. | GREEN | passes |

Rationale for collapsing (v8's two-commit plan retired by R-F07): the
"failing-gate-in-git-history" signal is recoverable via `git log -p` on
the test file alone (the assertion change is visible in the diff), and
the collapsed form keeps the tree green for pull/bisect use by other
developers.  No loss of regression-gate evidence.

### 14.2  Source edit — interior-fault branch (`wave_operator.inl:811-887`)

**BEFORE** (the buggy site):

```cpp
// wave_operator.inl:811-816
const real_t *Q_self_imp = elem1_on_plus ? Q_imp_plus_g : Q_imp_minus_g;
const real_t *Q_nbr_imp  = elem1_on_plus ? Q_imp_minus_g : Q_imp_plus_g;
real_t F_h_total[NUM_STATE];
flux_.Interior(nor, Q_self_imp, Q_nbr_imp, F_h_total);
```

... plus the DIAG C-2 block at `:818-843` (prints `F_h_total`), DIAG C-3
block at `:845-879` (brackets the Elem1 accumulation, prints `F_h_total[VX]`),
Elem1 accumulation at `:857-864` (`rhs[Elem1] -= w·shape1·F_h_total`), and
Elem2 accumulation at `:880-887` (`rhs[Elem2] += w·shape2·F_h_total`).

**AFTER**:

```cpp
// Two one-sided A·T·Q fluxes per Pelties 2012 eq. (9).  Using
// flux_.Interior(can_n, Q, Q) exploits the identity A^+ + A^- = A
// (see godunov_flux.cpp:94-105 eigenvalue construction) so the global
// result is T · A_x · Tinv · Q = A_{can_n} · Q (the normal-direction
// flux Jacobian applied to a single imposed state).
real_t F_h_plus[NUM_STATE], F_h_minus[NUM_STATE];
flux_.Interior(can_n, Q_imp_plus_g,  Q_imp_plus_g,  F_h_plus);
flux_.Interior(can_n, Q_imp_minus_g, Q_imp_minus_g, F_h_minus);

// F_h_elem1 is the Elem1 side's own-side flux — used only by the DIAG
// prints below.  Name chosen (per R-F04 rename) to avoid confusion with
// the shared-fault branch's `Q_imp_side` (MPI-rank-local naming).
const real_t *F_h_elem1 = elem1_on_plus ? F_h_plus : F_h_minus;

#ifdef SEAS_DIAG_FAULT_FLUX
// C-2 FLUX — print Elem1's own-side flux (in global frame).
if (dof_idx >= 0 &&
    dof_idx < static_cast<int>(fault_dof_data_->size()) &&
    (*fault_dof_data_)[dof_idx].diag_print) {
   real_t F_v_mag = std::sqrt(
      F_h_elem1[VX]*F_h_elem1[VX] + F_h_elem1[VY]*F_h_elem1[VY]
    + F_h_elem1[VZ]*F_h_elem1[VZ]);
   real_t F_s_max = 0.0;
   for (int c = 0; c < 6; c++) {
      F_s_max = std::max(F_s_max, std::abs(F_h_elem1[c]));
   }
   std::fprintf(stderr,
      "[C-2 FLUX] rank=%d  dof=%d  side=%c  |F_v|=%.3e m2/s2  "
      "max|F_stress|=%.3e Pa*m/s  F_h[VX]=%+.3e  "
      "F_h[SXY]=%+.3e  F_h[VY]=%+.3e  F_h[SXZ]=%+.3e\n",
      g_seas_my_rank, dof_idx, elem1_on_plus ? '+' : '-',
      F_v_mag, F_s_max,
      F_h_elem1[VX], F_h_elem1[SXY], F_h_elem1[VY], F_h_elem1[SXZ]);
}
#endif

#ifdef SEAS_DIAG_FAULT_FLUX
const bool _c3_diag =
   (dof_idx >= 0 &&
    dof_idx < static_cast<int>(fault_dof_data_->size()) &&
    (*fault_dof_data_)[dof_idx].diag_print);
const int  _c3_probe = VX * ndof_total_ + dof_offset1 + 0;
const real_t _c3_pre = _c3_diag ? rhs[_c3_probe] : 0.0;
#endif

// Per-side DG assembly.  Each element gets its OWN imposed-state flux,
// signed by its outward normal relative to can_n.
//   Elem1 outward = +can_n  ⇒ rhs -= w·shape1·F_h_plus    (elem1_on_plus)
//   Elem1 outward = -can_n  ⇒ rhs += w·shape1·F_h_minus   (else)
// Elem2 is the opposite side; its outward normal flips the sign.
if (elem1_on_plus) {
   for (int c = 0; c < NUM_STATE; c++) {
      for (int i = 0; i < ndof; i++) {
         rhs[c * ndof_total_ + dof_offset1 + i] -= w * shape1(i) * F_h_plus[c];
         rhs[c * ndof_total_ + dof_offset2 + i] += w * shape2(i) * F_h_minus[c];
      }
   }
} else {
   for (int c = 0; c < NUM_STATE; c++) {
      for (int i = 0; i < ndof; i++) {
         rhs[c * ndof_total_ + dof_offset1 + i] += w * shape1(i) * F_h_minus[c];
         rhs[c * ndof_total_ + dof_offset2 + i] -= w * shape2(i) * F_h_plus[c];
      }
   }
}

#ifdef SEAS_DIAG_FAULT_FLUX
if (_c3_diag) {
   const real_t _c3_post = rhs[_c3_probe];
   std::fprintf(stderr,
      "[C-3 RHS]  rank=%d  dof=%d  side=%c  "
      "rhs[VX,elem1,dof0] pre=%+.3e post=%+.3e  "
      "delta=%+.3e  w=%.3e  shape1(0)=%.3e  F_h_elem1[VX]=%+.3e\n",
      g_seas_my_rank, dof_idx, elem1_on_plus ? '+' : '-',
      _c3_pre, _c3_post, _c3_post - _c3_pre, w, shape1(0),
      F_h_elem1[VX]);
}
#endif
```

Net change: ~45 lines replace the original ~30 lines at `:811-887`.

**Welded-fallback `else` branch** (`wave_operator.inl:889-905` in pre-fix
numbering) — the branch taken when `dof_idx < 0` or `qpd_ptr == nullptr`
at a fault face.  Pre-fix it silently falls back to the welded Riemann
flux `flux_.Interior(nor, Q_self, Q_nbr, F_h)`; per **R-F02 (CRITICAL)**
this would bypass Pelties eq. (9) for any fault QP that reaches it.
Post-fix: hard-abort so a ctor-population bug surfaces instead of
silently masking as under-radiation.

**BEFORE**:

```cpp
// wave_operator.inl:889-905 (pre-fix)
else {
   // Fault face but no DOFData mapping (or FaultBasis not populated) —
   // fall back to welded interior flux.
   flux_.Interior(nor, Q_self, Q_nbr, F_h);
   for (int c = 0; c < NUM_STATE; c++) {
      for (int i = 0; i < ndof; i++) {
         rhs[c * ndof_total_ + dof_offset1 + i] -= w * shape1(i) * F_h[c];
         rhs[c * ndof_total_ + dof_offset2 + i] += w * shape2(i) * F_h[c];
      }
   }
}
```

**AFTER**:

```cpp
else {
   // R-F02 (v9.0.0): welded-flux fallback at a fault face would bypass
   // Pelties-9 eq. (9) per-side form and silently under-radiate that
   // QP.  A missing DOFData/FaultBasis mapping at a face flagged as a
   // fault is a ctor-population bug — force it to surface rather than
   // mask as reduced bulk coupling.
   MFEM_ABORT("interior fault face f=" << f
              << " has no FaultBasis/DOFData mapping ("
              << "fb_idx=" << fb_idx
              << ", dof_idx=" << dof_idx
              << ", have_basis=" << have_basis << ").  "
              << "v9.0.0 Pelties-9 per-side flux requires a valid mapping; "
              << "welded-flux fallback is no longer physical.");
}
```

### 14.3  Source edit — shared-fault branch (`wave_operator.inl:1232-1244`)

**BEFORE**:

```cpp
// wave_operator.inl:1232-1241
flux_.Interior(can_n, Q_imp_plus_g, Q_imp_minus_g, F_h);
const real_t accum_sign = elem1_on_plus ? +1.0 : -1.0;
for (int c = 0; c < NUM_STATE; c++) {
   for (int i = 0; i < ndof; i++) {
      rhs[c * ndof_total_ + dof_offset1 + i] -=
         accum_sign * w * shape1(i) * F_h[c];
   }
}
continue;
```

**AFTER**:

```cpp
// One-sided A·T·Q per Pelties 2012 eq. (9).  This rank's Elem1 is on
// exactly one side of the fault; its own imposed state drives its rhs.
// The paired rank writes the opposite-side contribution separately.
real_t F_h_side[NUM_STATE];
const real_t *Q_imp_side = elem1_on_plus ? Q_imp_plus_g : Q_imp_minus_g;
flux_.Interior(can_n, Q_imp_side, Q_imp_side, F_h_side);

// Elem1 outward = +can_n (elem1_on_plus) ⇒ rhs -= w·shape1·F_h_side
// Elem1 outward = -can_n (else)          ⇒ rhs += w·shape1·F_h_side
const real_t assemble_sign = elem1_on_plus ? -1.0 : +1.0;
for (int c = 0; c < NUM_STATE; c++) {
   for (int i = 0; i < ndof; i++) {
      rhs[c * ndof_total_ + dof_offset1 + i] +=
         assemble_sign * w * shape1(i) * F_h_side[c];
   }
}
continue;
```

Net change: 10 lines replace 10 lines at `:1232-1241`.

**Shared-fault fallback `else` branch** (`wave_operator.inl:1246-1249`
in pre-fix numbering) — same R-F02 concern: at a shared fault face
where the DOFData lookup fails, the current code silently welds.
Post-fix: hard-abort for the same reasons as §14.2's fallback.

**BEFORE**:

```cpp
// wave_operator.inl:1246-1249 (pre-fix)
else {
   flux_.Interior(nor, Q_self, Q_nbr, F_h);
}
```

**AFTER**:

```cpp
else {
   // R-F02 (v9.0.0): see §14.2 for rationale.  Welded-flux fallback at
   // a SHARED fault face is equally unacceptable; abort instead.
   MFEM_ABORT("shared fault face sf=" << sf
              << " has no DOFData mapping (dof_idx=" << dof_idx
              << ").  v9.0.0 Pelties-9 per-side flux requires a valid "
              << "mapping; welded-flux fallback is no longer physical.");
}
```

(The outer `else` for NON-fault shared faces at `:1251-1254` is
unchanged — standard welded flux is still correct there.)

### 14.4  Build & run checklist  (revised round 3)

```bash
cd /Users/chunhuizhao/projects/seas-mfem/miniapps/seas

# Edit sources per §14.2 / §14.3 (main-path replacement + two hard-abort fallbacks).
# Per R-F06: add header comment to `dynamic/godunov_flux.hpp` near line 84:
#   // NOTE: BuildJacobian must remain public because
#   //   tests/unit/test_godunov_interior_equal_sides_identity.cpp (§15.4)
#   // reconstructs A_x externally to verify the identity underpinning
#   // the v9.0.0 Pelties-9 per-side fix (see wave_operator.inl:811-887).
#   void BuildJacobian(int dir, DenseMatrix &A) const;
#
# Per R-F08: add bimaterial guard in `dynamic/fault_face_flux.cpp::Evaluate`,
# near the top of the function (e.g. line 75, before Step 1):
#   // v9.0.0 Pelties-9 per-side flux (plan §10.1) uses a single
#   // GodunovFlux instance and therefore assumes A_plus == A_minus.
#   MFEM_VERIFY(data.Zp_plus == data.Zp_minus && data.Zs_plus == data.Zs_minus,
#               "Bimaterial fault face detected; v9.0.0 Pelties-9 fix assumes "
#               "homogeneous material.  Extend GodunovFlux to per-side A "
#               "before running this configuration.");
#
# Edit tests per §15.2 / §15.3.  Add new tests per §15.4 / §15.5.
# Register new tests in the Makefile (follow existing
# seas_test_fault_face_flux_frame_and_flux pattern).

conda activate mfem-dev
make clean
make -j8 seas_test_fault_face_flux seas_test_fault_face_flux_frame \
         seas_test_fault_face_flux_frame_and_flux \
         seas_test_godunov_identity_normal_reversal \
         seas_test_fault_flux_interior_vs_shared_branch_equivalence \
         seas_test_godunov_interior_equal_sides_identity \
         seas_test_fault_face_flux_per_side_assembly

# Run per-test (fast):
make test-fault-face-flux-frame
make test-fault-face-flux-frame-and-flux   # was 7/9, expect 9+/9+ after fix
make test-godunov-identity-normal-reversal
make test-fault-flux-interior-vs-shared-branch-equivalence
make test-godunov-interior-equal-sides-identity         # NEW (§15.4)
make test-fault-face-flux-per-side-assembly             # NEW (§15.5)

# Full unit-test regression:
make test

# R-F05: DIAG build smoke — verify SEAS_DIAG_FAULT_FLUX path still
# compiles post-rename and still emits [C-*] lines.  make test does NOT
# exercise this path by default; without this step, a stale F_h_total
# reference inside the ifdef block would ship unnoticed.
rm -f drivers/tpv102_driver.o dynamic/*.o
make seas_tpv102_driver SEAS_EXTRA_CPPFLAGS="-DSEAS_DIAG_FAULT_FLUX" \
    2>&1 | tee /tmp/post_fix_diag_build.log
grep -qiE 'error|undefined reference' /tmp/post_fix_diag_build.log \
    && { echo "R-F05 FAIL: DIAG build regressed"; exit 1; }
# Brief serial run — confirm [C-2 FLUX] fires with the new `side=` tag:
./seas_tpv102_driver --mesh tpv102/mesh/tpv102_1000m.msh --tfinal 0.001 \
    2>&1 | grep '\[C-2 FLUX\]' | head -1
# Expected: at least one [C-2 FLUX] ... side=+ or side=- ... line.

# Parallel smoke (must still PASS — no behaviour change for non-fault faces):
mpirun -np 4 seas_test_parallel_elasticity
mpirun -np 8 seas_test_bp5_parallel_smoke
```

**Pass criteria** (all must hold for R-1001 to be closable):
- §3.1 12/12 PASS (unchanged — tests Evaluate only).
- §3.1b 9/9 PASS with `|F_h_plus[VX]|, |F_h_minus[VX]| ∈ (10, 1e4)`.
- §3.1c 5/5 PASS at ULP (Godunov identity property unchanged).
- §3.1d 4/4 PASS after the updated per-side logic is mirrored in the test.
- §3.1e (new) primitive-identity test PASSES.
- §3.1f (new) per-side-flux-asymmetry test PASSES.
- All other `make test` targets stay green.

---

## 15.  Unit-test design (added 2026-04-20, round 3)

### 15.1  Existing test impact matrix

| Test | File | Exercises | Impact of fix |
|---|---|---|---|
| §3.1  | `test_fault_face_flux_frame.cpp` | Evaluate (Eqs. 7–12) in isolation | **None** — imposed-state construction unchanged. |
| §3.1b | `test_fault_face_flux_frame_and_flux.cpp` | Full chain to `flux_.Interior` | **Update required** (§15.2): call per-side. |
| §3.1c | `test_godunov_identity_normal_reversal.cpp` | `F(+n,L,R) + F(-n,R,L) = 0` | **None** — property of Interior, unchanged. |
| §3.1d | `test_fault_flux_interior_vs_shared_branch_equivalence.cpp` | Interior vs shared branch agreement on Elem1 | **Update required** (§15.3): both branches now per-side. |

### 15.2  §3.1b update — `test_fault_face_flux_frame_and_flux.cpp`

Replace the Step 3 block at `tests/unit/test_fault_face_flux_frame_and_flux.cpp:158-201` with:

```cpp
// ---- Step 3: per-side A·T·Q fluxes (post-fix, Pelties 2012 eq. 9) ----
GodunovFlux flux(TPV102Params::lambda, TPV102Params::mu, TPV102Params::rho);
real_t F_h_plus[NUM_STATE], F_h_minus[NUM_STATE];
flux.Interior(can_n, Q_imp_plus_g,  Q_imp_plus_g,  F_h_plus);
flux.Interior(can_n, Q_imp_minus_g, Q_imp_minus_g, F_h_minus);

auto print_flux = [](const char *label, const real_t *F) {
   std::cout << "  [Step 3] " << label
             << ": VX=" << F[VX] << " VY=" << F[VY] << " VZ=" << F[VZ]
             << " SXX=" << F[SXX] << " SXY=" << F[SXY] << " SXZ=" << F[SXZ] << "\n";
};
print_flux("F_h_plus",  F_h_plus);
print_flux("F_h_minus", F_h_minus);

// §3.1b regression gate: each side must produce a non-zero, O(100) m²/s²
// VX-component flux (Pelties eq. 9 + post-breakaway fixture).
real_t abs_F_plus_VX  = std::abs(F_h_plus[VX]);
real_t abs_F_minus_VX = std::abs(F_h_minus[VX]);
std::cout << "  DECISION: |F_h_plus[VX]|  = " << abs_F_plus_VX << " m²/s²\n";
std::cout << "  DECISION: |F_h_minus[VX]| = " << abs_F_minus_VX << " m²/s²\n";

TEST_ASSERT(abs_F_plus_VX  > 1e-10, "|F_h_plus[VX]|  > 1e-10 (no cancellation)");
TEST_ASSERT(abs_F_minus_VX > 1e-10, "|F_h_minus[VX]| > 1e-10 (no cancellation)");
TEST_ASSERT(abs_F_plus_VX  > 10.0 && abs_F_plus_VX  < 1.0e4,
            "|F_h_plus[VX]|  in (10, 1e4) m²/s²");
TEST_ASSERT(abs_F_minus_VX > 10.0 && abs_F_minus_VX < 1.0e4,
            "|F_h_minus[VX]| in (10, 1e4) m²/s²");

// Physics check: VX flux is stress-driven (A[VX,σ] ≠ 0, A[VX,VX] = 0).
// Plus/minus imposed states share the SAME corrected traction σ_corr
// (Pelties eq. 11d/12d) and the T_can rotation is identical for both
// sides, so Q_imp_{plus,minus}_g[SXY] are bit-identical, and
// F_{can_n}·Q_imp_plus_g and F_{can_n}·Q_imp_minus_g are ULP-identical
// in their VX row.  Per R-F03, a 5% tolerance was ~14 OOM too loose
// — tighten to 1e-10 so a subtle FP-ordering drift between sides
// surfaces as a test failure.
real_t rel_diff_VX = std::abs(F_h_plus[VX] - F_h_minus[VX]) /
                     std::max(abs_F_plus_VX, abs_F_minus_VX);
TEST_ASSERT(rel_diff_VX < 1.0e-10,
            "F_h_plus[VX] == F_h_minus[VX] to ULP (stress continuity)");

// Sign check — analytically, F_h_{plus,minus}[VX] for can_n=(0,-1,0)
// under right-lateral strike slip (V2 > 0, τ2_corr = -η_s·V2 < 0) is:
//   F_{can_n} = n·A = 0·A_x + (-1)·A_y + 0·A_z = -A_y
//   A_y[VX][SXY] = -1/ρ     (from BuildJacobian(dir=1), verified)
//   Q_imp_{plus,minus}_g[SXY] = +η_s·V2 (the rotation flips the sign
//     of the local τ_{n,t2}=-η_s·V2 because can_n is the -y axis).
//   F_h_{plus,minus}[VX] = (-A_y[VX][SXY]) · Q_imp_{plus,minus}_g[SXY]
//                        = -(-1/ρ) · (+η_s·V2)
//                        = +η_s·V2/ρ = +c_s·V2/2  (POSITIVE)
// A common derivation mistake (seen in round-3 REVIEW R-F01) is to
// apply A_y directly to Q_imp_g without the `n·A` sign flip, which
// gives the wrong sign.  The DG rhs update `rhs -= w·shape·F_h_plus`
// then produces rhs[Elem_plus][VX] < 0 — the sign asserted below and
// in §15.5 T-F.  See §18 for the full derivation and history.
TEST_ASSERT(F_h_plus[VX]  > 0.0,
            "F_h_plus[VX]  > 0 analytically (= +c_s·V2/2 for V2 > 0)");
TEST_ASSERT(F_h_minus[VX] > 0.0,
            "F_h_minus[VX] > 0 analytically (same value as plus)");

// Stress flux: the velocity jump across the fault is opposite on the
// two sides (ṽ^+ − v^+ = +δ, ṽ^- − v^- = −δ in Pelties eq. 15), and
// A[SXY, v] ≠ 0.  So F_h_plus[SXY] and F_h_minus[SXY] must be ≈
// opposite-signed (velocity-driven channel).
TEST_ASSERT(F_h_plus[SXY] * F_h_minus[SXY] < 0.0 ||
            (std::abs(F_h_plus[SXY]) < 1e-6 && std::abs(F_h_minus[SXY]) < 1e-6),
            "F_h_plus[SXY] and F_h_minus[SXY] have opposite signs (velocity-driven)");
```

Test count grows from 9 to ~12; currently-failing assertion `|F_h[VX]|
> 10` gets replaced by the stronger `|F_h_plus[VX]| > 10 ∧ |F_h_minus[VX]|
> 10`.

### 15.3  §3.1d update — `test_fault_flux_interior_vs_shared_branch_equivalence.cpp`

Rationale: both branches now use per-side logic, so the test must
compare per-side Elem1 contributions, not the old single-F_h form.

Replace the "interior branch replica" block around
`tests/unit/test_fault_flux_interior_vs_shared_branch_equivalence.cpp:148-155` with:

```cpp
// Post-fix interior-fault branch: per-side Interior call.
real_t F_h_plus[NUM_STATE], F_h_minus[NUM_STATE];
flux.Interior(can_n, Q_imp_plus_g,  Q_imp_plus_g,  F_h_plus);
flux.Interior(can_n, Q_imp_minus_g, Q_imp_minus_g, F_h_minus);

// Elem1 on plus:  rhs -= w·shape1·F_h_plus
// Elem1 on minus: rhs += w·shape1·F_h_minus
// Fold into a single "Elem1 contribution" expression:
const real_t *F_h_elem1 = elem1_on_plus ? F_h_plus : F_h_minus;
const real_t interior_sign = elem1_on_plus ? -1.0 : +1.0;
for (int c = 0; c < NUM_STATE; c++) {
   F_elem1_interior[c] = interior_sign * F_h_elem1[c];
}
```

Replace the "shared branch replica" block around `:167-172` with:

```cpp
// Post-fix shared-fault branch: one-sided Interior call per rank.
real_t F_h_shared[NUM_STATE];
const real_t *Q_imp_side = elem1_on_plus ? Q_imp_plus_g : Q_imp_minus_g;
flux.Interior(can_n, Q_imp_side, Q_imp_side, F_h_shared);

// rhs[Elem1] += assemble_sign·w·shape1·F_h_shared   (assemble_sign as §14.3)
const real_t shared_sign = elem1_on_plus ? -1.0 : +1.0;
for (int c = 0; c < NUM_STATE; c++) {
   F_elem1_shared[c] = shared_sign * F_h_shared[c];
}
```

Equivalence assertion on `F_elem1_interior == F_elem1_shared` stays
unchanged; both branches now produce identical Elem1 contribution
because they both evaluate `±A_{can_n} · Q_imp_side` for this rank's
elem1_on_plus.  Expected: 4/4 PASS at ULP (unchanged from pre-fix).

### 15.4  New test T-E — primitive identity  (file: `test_godunov_interior_equal_sides_identity.cpp`)

Purpose: verify the identity `flux_.Interior(n, Q, Q) ≡ A · T · Tinv ·
Q` (where A is the full flux Jacobian in the x-direction) for arbitrary
Q and n.  This is the identity the fix depends on.

```cpp
// tests/unit/test_godunov_interior_equal_sides_identity.cpp
//
// Verifies: flux_.Interior(n, Q, Q) = T · (A_x^+ + A_x^-) · Tinv · Q
//                                   = T · A_x · Tinv · Q
// for arbitrary unit n and arbitrary physically-sane Q.  This is the
// identity that the v9.0.0 Option A fix relies on.

#include "mfem.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../config/tpv102_params.hpp"

using namespace mfem;
using namespace mfem::seas;

static int nt = 0, np = 0, nf = 0;
#define T_NEAR(v, e, tol, msg) do { nt++; \
   if (std::abs((v)-(e)) < (tol)) { np++; std::cout<<"  PASS: "<<msg<<"\n"; } \
   else { nf++; std::cout<<"  FAIL ["<<__LINE__<<"]: "<<msg \
      <<" (got "<<(v)<<", expected "<<(e)<<")\n"; } } while (0)

// Build A_x (full Jacobian) via BuildJacobian, then compute expected
// F = T · A_x · Tinv · Q.
static void ExpectedFullFlux(const GodunovFlux &flux,
                             const real_t *n, const real_t *Q, real_t *F) {
   real_t t1[3], t2[3];
   GodunovFlux::BuildFrame(n, t1, t2);
   DenseMatrix T(NUM_STATE), Tinv(NUM_STATE), Ax(NUM_STATE);
   GodunovFlux::BuildRotation(n, t1, t2, T);
   GodunovFlux::BuildRotationInverse(n, t1, t2, Tinv);
   flux.BuildJacobian(0, Ax);   // full A_x, not A_x^±

   real_t Q_rot[NUM_STATE], tmp[NUM_STATE], F_rot[NUM_STATE];
   Tinv.Mult(Q, Q_rot);
   Ax.Mult(Q_rot, F_rot);
   T.Mult(F_rot, F);
}

int main() {
   GodunovFlux flux(TPV102Params::lambda, TPV102Params::mu, TPV102Params::rho);

   // 5 normals (covers x/y/z axes, up-fallback branch, oblique).
   real_t normals[5][3] = {
      { 1.0,  0.0,  0.0},
      { 0.0,  1.0,  0.0},
      { 0.0,  0.0,  1.0},   // triggers BuildFrame's up=(1,0,0) fallback
      { 1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0), 1.0/std::sqrt(3.0)},
      { 0.6, -0.8,  0.0}
   };
   real_t Q_samples[3][NUM_STATE] = {
      // stress-only (tests A[v, σ] rows)
      {1.0e6, 2.0e6, 3.0e6, 4.0e6, 5.0e6, 6.0e6, 0.0, 0.0, 0.0},
      // velocity-only (tests A[σ, v] rows)
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, -0.5, 0.3},
      // mixed physical (velocity O(1), stress O(MPa))
      {1e6, -2e6, 0.5e6, 3e5, -4e5, 6e5, 0.1, -0.2, 0.15}
   };

   int case_id = 0;
   for (int ni = 0; ni < 5; ni++) {
      for (int qi = 0; qi < 3; qi++) {
         real_t F_interior[NUM_STATE], F_expected[NUM_STATE];
         flux.Interior(normals[ni], Q_samples[qi], Q_samples[qi], F_interior);
         ExpectedFullFlux(flux, normals[ni], Q_samples[qi], F_expected);
         for (int c = 0; c < NUM_STATE; c++) {
            real_t tol = std::max(1.0e-10,
                                  1.0e-12 * std::abs(F_expected[c]));
            T_NEAR(F_interior[c], F_expected[c], tol,
                   "case " + std::to_string(case_id) + " c=" + std::to_string(c));
         }
         case_id++;
      }
   }

   std::cout << "Total=" << nt << " Passed=" << np << " Failed=" << nf << "\n";
   return (nf > 0) ? 1 : 0;
}
```

Expected: all ULP-level agreements (5 normals × 3 Q-samples × 9
components = 135 assertions).  Any failure implies the identity fails —
reject the fix until understood.

Note: requires `GodunovFlux::BuildJacobian` to be callable externally.
It currently is (public in `godunov_flux.hpp`).

### 15.5  New test T-F — per-side flux asymmetry  (file: `test_fault_face_flux_per_side_assembly.cpp`)

Purpose: verify the per-side assembly logic of §8.1 independent of
`wave_operator.inl` — i.e., check the algebra of the four
`(elem1_on_plus, c)` sign combinations using a minimal in-test rhs
scalar buffer.

```cpp
// tests/unit/test_fault_face_flux_per_side_assembly.cpp
//
// Reproduce §14.2's per-side DG assembly with a tiny rhs buffer and
// assert that (a) both Elem1 and Elem2 get non-zero contributions,
// (b) they are NOT ±-related (which was the welded-face assumption),
// (c) the signs match the outward-normal convention.

#include "mfem.hpp"
#include "../../dynamic/fault_face_flux.hpp"
#include "../../dynamic/godunov_flux.hpp"
#include "../../dynamic/wave_state.hpp"
#include "../../config/tpv102_params.hpp"

using namespace mfem; using namespace mfem::seas;
static int nt=0, np=0, nf=0;
#define T_OK(cond,msg) do{nt++; if(cond){np++;std::cout<<"  PASS: "<<msg<<"\n";} \
   else{nf++;std::cout<<"  FAIL ["<<__LINE__<<"]: "<<msg<<"\n";}}while(0)

static void RunCase(bool elem1_on_plus) {
   // Same TPV102 post-breakaway hypocentre fixture as §3.1b.
   real_t can_n[3]={0,-1,0}, can_t1[3]={0,0,-1}, can_t2[3]={1,0,0};
   DenseMatrix T_can(NUM_STATE);
   GodunovFlux::BuildRotation(can_n, can_t1, can_t2, T_can);

   DOFData d;
   d.Zp_plus=d.Zp_minus=TPV102Params::Zp; d.Zs_plus=d.Zs_minus=TPV102Params::Zs;
   d.eta_p=TPV102Params::Zp/2; d.eta_s=TPV102Params::Zs/2;
   d.sigma_n0=TPV102Params::sigma_n; d.tau1_0=0.0; d.tau2_0=81e6;
   d.a=TPV102Params::a_vw; d.Dc=TPV102Params::Dc;
   d.psi = d.a * std::log(2.0*TPV102Params::V0/1.0 *
                          std::sinh(81e6/(d.sigma_n0*d.a)));

   real_t Q_plus_local[NUM_STATE]={}, Q_minus_local[NUM_STATE]={};
   real_t Q_imp_plus_local[NUM_STATE], Q_imp_minus_local[NUM_STATE];
   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   ff.Evaluate(d, Q_plus_local, Q_minus_local,
               Q_imp_plus_local, Q_imp_minus_local);

   real_t Q_imp_plus_g[NUM_STATE], Q_imp_minus_g[NUM_STATE];
   T_can.Mult(Q_imp_plus_local,  Q_imp_plus_g);
   T_can.Mult(Q_imp_minus_local, Q_imp_minus_g);

   GodunovFlux flux(TPV102Params::lambda, TPV102Params::mu, TPV102Params::rho);
   real_t F_plus[NUM_STATE], F_minus[NUM_STATE];
   flux.Interior(can_n, Q_imp_plus_g,  Q_imp_plus_g,  F_plus);
   flux.Interior(can_n, Q_imp_minus_g, Q_imp_minus_g, F_minus);

   // Per-side assembly (§14.2 replica, one DOF per element for simplicity).
   real_t rhs_elem1[NUM_STATE]={}, rhs_elem2[NUM_STATE]={};
   const real_t w=1.0, shape1=1.0, shape2=1.0;
   if (elem1_on_plus) {
      for (int c = 0; c < NUM_STATE; c++) {
         rhs_elem1[c] -= w*shape1*F_plus[c];
         rhs_elem2[c] += w*shape2*F_minus[c];
      }
   } else {
      for (int c = 0; c < NUM_STATE; c++) {
         rhs_elem1[c] += w*shape1*F_minus[c];
         rhs_elem2[c] -= w*shape2*F_plus[c];
      }
   }

   const char *tag = elem1_on_plus ? "elem1_on_plus=true" : "elem1_on_plus=false";
   T_OK(std::abs(rhs_elem1[VX]) > 10.0, std::string("|rhs_Elem1[VX]| > 10  ") + tag);
   T_OK(std::abs(rhs_elem2[VX]) > 10.0, std::string("|rhs_Elem2[VX]| > 10  ") + tag);

   // Welded-face pattern would give rhs_elem2 = -rhs_elem1.  At a fault
   // they are NOT ±-related because Q_imp_plus ≠ Q_imp_minus in velocity.
   // SXY is the velocity-driven channel where the asymmetry lives most.
   real_t sum_SXY = rhs_elem1[SXY] + rhs_elem2[SXY];
   T_OK(std::abs(sum_SXY) > 1.0,
        std::string("rhs_Elem1[SXY] + rhs_Elem2[SXY] ≠ 0 (not welded)  ") + tag);

   // Outward-normal sign check — SIGN DERIVATION (important: see §18
   // for history; the round-3 REVIEW R-F01 proposed inverting these
   // asserts, incorrectly, by applying A_y instead of A_{can_n}).
   //
   // For right-lateral strike-slip TPV102 (V2 > 0, τ2_corr = -η_s·V2):
   //
   //   can_n = (0, -1, 0)  ⇒  A_{can_n} = n·A = -A_y
   //   A_y[VX][SXY] = -1/ρ         (verified: BuildJacobian(dir=1))
   //   Q_imp_{plus,minus}_g[SXY] = +η_s·V2  (common to both sides; the
   //                                          local τ_{n,t2} = -η_s·V2
   //                                          flips sign under rotation
   //                                          with can_n=-y, can_t2=+x).
   //
   //   F_h_{plus,minus}[VX]
   //     = A_{can_n}[VX][SXY] · Q_imp_side_g[SXY]
   //     = (-A_y[VX][SXY]) · Q_imp_side_g[SXY]
   //     = -(-1/ρ) · (+η_s·V2)
   //     = +η_s·V2/ρ  =  +c_s·V2/2    ⟹ POSITIVE.
   //
   // DG rhs assembly (plus side, outward normal = +can_n):
   //     rhs -= w·shape·F_h_plus  ⇒  rhs_plus[VX] < 0.
   // DG rhs assembly (minus side, outward normal = -can_n):
   //     rhs += w·shape·F_h_minus ⇒  rhs_minus[VX] > 0.
   //
   // The common sign-error is to compute F_h = A_y · Q directly,
   // omitting the `n·A` projection; this flips the sign and leads to
   // asserts that fail on correct code.  Do not fall into that trap.
   bool plus_side_rhs_sign = elem1_on_plus
                             ? (rhs_elem1[VX] < 0.0)
                             : (rhs_elem2[VX] < 0.0);
   bool minus_side_rhs_sign = elem1_on_plus
                              ? (rhs_elem2[VX] > 0.0)
                              : (rhs_elem1[VX] > 0.0);
   T_OK(plus_side_rhs_sign && minus_side_rhs_sign,
        std::string("rhs_plus[VX]<0 AND rhs_minus[VX]>0 (F_h>0, outward-normal) ")
        + tag);
}

int main() {
   RunCase(true);
   RunCase(false);
   std::cout << "Total=" << nt << " Passed=" << np << " Failed=" << nf << "\n";
   return (nf > 0) ? 1 : 0;
}
```

Expected: 8/8 PASS (4 assertions × 2 cases).  This is the test that
would have caught the old "one F_h, two signs" welded-face assembly
even if the flux cancellation hadn't been there.

### 15.6  Optional parallel extension (deferred)

Extend `tests/parallel/test_r101_shared_fault.cpp` to verify that the
two ranks sharing a fault face produce Elem1 contributions consistent
with each other under the per-side assembly — i.e., rank A's Elem1
(plus side) and rank B's Elem1 (minus side) jointly give `-F_plus` and
`+F_minus` respectively, and the ghost Q exchange is symmetric.

Not required for R-1001 close — existing R-101 verifier (Appendix E.1)
already checks DOFData cross-rank pairing at the level needed.
Defer to follow-up.

---

## 16.  Self-audit of the proposed patch (added 2026-04-20, round 3)

Bug-hunt on the pseudocode in §8 and the concrete code in §14 — before
the fix agent writes it.  Each finding has an ID the fix agent can cite.

### 16.1  Findings

**P-001 [CRITICAL resolved] — DIAG block variable rename.**

*Risk:* The DIAG C-2 / C-3 blocks at `wave_operator.inl:818-843` and
`:845-879` reference `F_h_total` (the pre-fix local name).  If the fix
replaces the local with `F_h_plus` / `F_h_minus` without updating the
DIAG blocks, the code fails to compile when `SEAS_DIAG_FAULT_FLUX` is
set.

*Resolution:* §14.2 replaces `F_h_total` references with `F_h_elem1`
(aliased from the elem1_on_plus ternary), and adds a `side=%c` flag
to each fprintf so the log captures which side was printed.  No
behaviour change when `SEAS_DIAG_FAULT_FLUX` is off.

*Test:* build both with and without `-DSEAS_DIAG_FAULT_FLUX=1`; both
must compile (not currently covered by `make test` — add
compile-smoke check to §14.4).

**P-002 [MODERATE] — `nor` vs `can_n` in the interior-fault branch.**

*Risk:* The pre-fix interior-fault site uses MFEM's `nor` (rank-local,
Elem1-outward).  §14.2 switches the `flux_.Interior` call to `can_n`
(canonical, bit-identical across ranks).  If any downstream code in
the same block references `nor`, behaviour changes subtly.

*Investigation:* grep inside `:811-887` — `nor` appears only in the
original `flux_.Interior(nor, ...)` call at line 816.  No other use.
Safe to switch.

*Test:* §15.5 T-F constructs both `elem1_on_plus` cases explicitly
using `can_n` and verifies the outward-normal sign convention through
`rhs_elem1` and `rhs_elem2`.

**P-003 [MODERATE] — `flux_.Interior(n, Q, Q)` identity relies on
A^+ + A^- = A exactly, not approximately.**

*Risk:* If the eigendecomposition in `godunov_flux.cpp:94-119` had a
bug that leaked non-zero off-diagonal entries into `Λ^+ + Λ^- − Λ`,
the identity would not hold to ULP, and per-side Interior calls would
accumulate rounding error.

*Investigation:* `Λ^+ = diag(cp, cs, cs, 0, 0, 0, 0, 0, 0)` and
`Λ^- = diag(0, 0, 0, 0, 0, 0, -cs, -cs, -cp)` by construction
(`godunov_flux.cpp:99-105`).  Their sum is the full eigenvalue matrix
`Λ = diag(cp, cs, cs, 0, 0, 0, -cs, -cs, -cp)`.  So `A^+ + A^- =
R(Λ^+ + Λ^-)R^{-1} = RΛR^{-1} = A` **exactly**, by construction.

*Test:* §15.4 T-E verifies the identity numerically at 135
assertions.  Any failure → reject fix.

**P-004 [LOW] — performance: two `Interior` calls per QP is 2× the
per-QP cost vs the pre-fix (non-working) single call.**

*Risk:* Doubles the fault-face flux cost.  For TPV102 at 200 m /
400-rank this is measurable but not dominant — fault faces are a
surface in a 3D volume, so O(h^2) vs the O(h^3) bulk work.

*Mitigation:* acceptable for v9.0.0.  Follow-up: add a
`GodunovFlux::ApplyFullFlux(n, Q, F)` method that does `T·A·T^{-1}·Q`
in one matvec (skipping the redundant `A_x^- · Q` since the two
halves sum to `A · Q`).  Defer — file as separate perf ticket.

**P-005 [LOW] — BuildFrame's internal tangent (t1, t2) differs from
FaultBasis's `can_t1` / `can_t2`.**

*Risk:* Inside `Interior`, BuildFrame generates its own orthogonal
(t1, t2) from the normal — these do not match `can_t1 = (0,0,-1)`,
`can_t2 = (1,0,0)` built in wave_operator.inl:758-760.  However, the
final F_h is in global coordinates (rotated back by `T = BuildRotation(...)`
inside `Interior`), so the choice of tangent is invisible to the caller.

*Investigation:* for `can_n = (0,-1,0)`, BuildFrame gives
t1 = (1, 0, 0), t2 = (0, 0, 1) — different from FaultBasis convention.
Verified by §15.4 T-E: the Interior output matches `T · A_x · T^{-1} · Q`
computed with BuildFrame's own frame.

*Consequence:* none for correctness.  Worth a one-line comment in the
patch explaining the internal frame is not the FaultBasis frame.

**P-006 [MODERATE — POSSIBLE] — per-side flux with homogeneous `A`
vs bimaterial.**

*Risk:* Pelties 2012 eq. (9) has a single `A` per face.  SeisSol
generalises to per-side A_plus / A_minus.  SEAS-MFEM's single
`GodunovFlux` instance (constructed from global ρ, λ, μ) means our
patch implicitly assumes homogeneous material across the fault.

*Consequence:* correct for TPV102 (homogeneous).  Would be wrong for
a bimaterial fault.  None of SEAS-MFEM's current benchmarks are
bimaterial.

*Mitigation:* add a runtime assertion in the ctor of `FaultFaceFlux`
that `data.Zp_plus == data.Zp_minus` (which the existing code already
implies but doesn't enforce).  Out of scope for R-1001 — file as
separate generalisation ticket.

**P-007 [LOW] — sign-sanity back-of-envelope value depends on actual
`tau2_corr` at fixture.**

*Risk:* §15.2 asserts `|F_h_plus[VX]| ∈ (10, 1e4)` based on the
estimate `τ/ρ ≈ 1.25e6/2670 ≈ 468`.  If TPV102 parameters (`rho`,
`mu`) are changed the fixture needs re-tuning.

*Mitigation:* the acceptance band `(10, 1e4)` is wide (three decades)
specifically to tolerate fixture drift.  Only fixture changes that
move `τ_corr` by >20× would fail.  Low risk.

### 16.2  Net verdict

- 1 CRITICAL finding, resolved by §14.2 DIAG rename.
- 2 MODERATE findings, resolved by investigation + tests.
- 4 LOW findings, mitigated or out of scope.
- **Plan compliance:** FULL — §14 / §15 implement precisely what §8
  prescribed, plus DIAG plumbing the §8 pseudocode didn't cover.
- **Verdict:** PASS WITH FIXES.  The fixes are already incorporated
  into §14 / §15 (no additional action needed beyond §14 / §15 as
  written).

---

## 18.  Round-3 review response (REVIEW.md findings R-F01…R-F08)

User requested a code-review on this plan on 2026-04-20.  Review saved
to `seas-mfem/REVIEW.md` with 8 findings (2 CRITICAL, 4 MODERATE, 2 LOW).
Dispositions and applied changes:

### 18.1  Accepted and applied in this revision

| ID | Severity | Finding | Resolution |
|----|----------|---------|------------|
| **R-F02** | CRITICAL | Welded-fallback `else` branch at `wave_operator.inl:889-905` (interior-fault) and `:1246-1249` (shared-fault) would silently use the welded Riemann flux at a fault QP whose DOFData lookup failed — bypassing Pelties-9. | §14.2 + §14.3 now specify `MFEM_ABORT` in both fallback branches.  A ctor-population bug surfaces instead of masking as reduced radiation. |
| **R-F03** | MODERATE | §15.2 tolerance `rel_diff_VX < 0.05` was ~14 OOM too loose; the two sides' `F_h[VX]` are bit-identical by stress continuity. | §15.2 tolerance tightened to `1.0e-10` with a comment explaining the exact-equality physics. |
| **R-F04** | MODERATE | `F_h_self` is ambiguous (DG-element vs MPI-rank sense); a clearer name is preferred. | Renamed `F_h_self` → `F_h_elem1` throughout §14.2. |
| **R-F05** | MODERATE | `make test` does not compile the `SEAS_DIAG_FAULT_FLUX` path, so a stale `F_h_total` reference inside an `#ifdef` block would ship unnoticed. | §14.4 now runs a post-fix DIAG-on build + a short driver run that greps for `[C-2 FLUX]` output as a smoke-test. |
| **R-F06** | MODERATE | §15.4 T-E silently depends on `GodunovFlux::BuildJacobian` staying public; a future refactor could break this without notice. | §14.4 now adds a header comment next to `BuildJacobian` declaration referencing the T-E test and the Pelties-9 fix. |
| **R-F07** | LOW | Two-commit "C1 RED, C2 GREEN" strategy leaves the tree broken for the duration of C2 review. | §14.1 revised to a single-commit strategy (everything GREEN at every push).  History of the regression-gate flip is still recoverable via `git log -p` on the test file. |
| **R-F08** | LOW | Plan's caveat §10.1 notes bimaterial is out of scope, but the patch doesn't enforce homogeneity at runtime. | §14.4 now adds an `MFEM_VERIFY(Zp_plus == Zp_minus && Zs_plus == Zs_minus, ...)` in `FaultFaceFlux::Evaluate`.  A future bimaterial run will abort with a pointer to this plan. |

### 18.2  Rejected with explanation

**R-F01 — sign of `rhs_plus[VX]` and `rhs_minus[VX]`**  ✗ **NOT APPLIED**

The reviewer asserts that the §15.5 T-F test's sign expectations are
inverted and would fail on correct code.  After re-deriving from first
principles, **the original plan is correct and the reviewer's proposed
inversion would make the test fail on correct code**.

Root cause of the review error — the trace at REVIEW.md line 59:

    F_h_plus[VX] = -(1/ρ) · Q_imp_plus_g[SXY] = -c_s·V/2

uses `A_y` directly as the flux Jacobian at the face.  But the flux at
a face with normal `can_n = (0, -1, 0)` is the **normal-direction**
Jacobian `A_{can_n} = n·A = -A_y` (there is exactly one contribution,
from `n_y = -1`).  Including the sign flip:

    A_{can_n}[VX][SXY]    = -A_y[VX][SXY] = -(-1/ρ) = +1/ρ
    Q_imp_plus_g[SXY]     = +η_s·V2       (rotation; common to both sides)
    F_h_plus[VX]          = A_{can_n}[VX][SXY] · Q_imp_plus_g[SXY]
                          = (+1/ρ) · (+η_s·V2)
                          = +c_s·V2/2    (POSITIVE, not negative)

Cross-check via `flux_.Interior(can_n, Q, Q)` in its internal BuildFrame
(t1 = (1,0,0) = global +x, t2 = (0,0,1) = global +z):

    Q_imp_plus_g rotates to:
      Q_rot[VY_local] = +x_component(Q_imp_plus_g) = Q_imp_plus_g[VX] = -V2/2
      Q_rot[SXY_local] = σ_{n,t1,local} = -η_s·V2   (global σ_xy projected)
    F_rot = A_x · Q_rot has only two non-zero entries:
      F_rot[VY_local] = A_x[VY][SXY] · Q_rot[SXY_local]
                      = (-1/ρ) · (-η_s·V2) = +c_s·V2/2
      F_rot[SXY_local] = A_x[SXY][VY] · Q_rot[VY_local]
                      = (-μ) · (-V2/2) = +μ·V2/2
    Rotate back: F_global[VX] = F_rot[VY_local] = +c_s·V2/2 ✓ (POSITIVE)

Both derivations give `F_h_plus[VX] = +c_s·V2/2 > 0`, which makes the
DG assembly produce `rhs_plus[VX] = -w·shape·F_h_plus[VX] < 0` — the
original plan assertion.  Applying the reviewer's inversion would
flip the test to `rhs_plus[VX] > 0` and the test would fail on the
correct code.

**Defensive actions taken in response to R-F01** (even though the fix
itself is rejected):

1. §15.2 (`test_fault_face_flux_frame_and_flux.cpp`) now asserts the
   sign of `F_h_plus[VX]` and `F_h_minus[VX]` explicitly (both > 0)
   with the derivation in-line.  A future reader will not have to
   reconstruct the `n·A` projection from scratch.
2. §15.5 (`test_fault_face_flux_per_side_assembly.cpp`) now includes a
   full sign-derivation comment block explicitly naming the
   `A_{can_n} = -A_y` projection and the common sign-error to avoid.
3. This §18 records the review disagreement and rationale so it's
   discoverable.

If the reviewer (or any future reader) still believes the signs are
inverted, the expected next step is to reproduce the 2-line derivation
above using `BuildJacobian(dir=1)` output, NOT to re-open without
showing an error in the derivation.

### 18.3  Unreviewed areas (carried forward from REVIEW.md §Unreviewed)

- Shared-fault MPI ghost-exchange post-fix: §15.6 (deferred parallel
  test) would close this.  Not a blocker.
- ADER time-weight factor: plan §10.3 is correct; spatial formula is
  independent of time integration.
- `BuildRotation` Voigt factor for oblique faults: §10.2 — out of scope
  for TPV102.
- 2× `Interior` cost at 400-rank production: §16.1 P-004 — acceptable
  for v9.0.0.
- Combined `SEAS_DIAG_FAULT_FLUX` + `SEAS_DIAG_GHOST_EXCHANGE` build:
  not currently exercised; note for a future DIAG smoke-test extension.

---

## 19.  What I need from you before executing

Read this document top-to-bottom.  Then answer:

- **Q0.** (Round 2.)  Does §A (reference-paper verification against
  Pelties 2012 eq. (9)/(15)/(20) and Dumbser-Käser 2006 eq. (14))
  convince you that the current code is inconsistent with the paper and
  that the proposed fix is paper-grounded — not just code-mimicry of
  SeisSol?
- **Q1.** Does §5 (the cancellation algebra) convince you that the
  Riemann call is the wrong primitive at a fault face?
- **Q2.** Does §3 + §6 + §9 (SeisSol's one-sided `A · T · Q_imp`, and
  why `flux_.Interior(n, Q, Q)` reproduces it exactly) convince you
  that Option A is a like-for-like port of SeisSol's formulation, not a
  new invention?
- **Q3.** Does §7 (the secondary assembly bug) convince you that the
  patch must change **both** the flux call and the rhs assembly
  pattern, not just the flux call?
- **Q4.** Are the caveats in §10 acceptable to leave as-is for v9.0.0?
- **Q5.** Approve to proceed with the (revised, round-3) commit plan in
  §12 (summary) / §14 (source edits + hard-abort fallbacks + header
  + bimaterial VERIFY + DIAG smoke) / §15 (test updates + new tests)?
- **Q6.** (Round 3.)  Does §16 (self-audit of the proposed patch)
  leave any residual concern that would block the patch?  If so,
  which P-xxx?
- **Q7.** (Round 3.)  Do you agree with the §18 dispositions of
  REVIEW.md findings R-F01…R-F08?  In particular, do you concur that
  R-F01 is a reviewer derivation error (explained in §18.2 with the
  `A_{can_n} = -A_y` cross-check) and should NOT be applied?  If you
  still believe R-F01 is correct, reply with the specific line in the
  §18.2 derivation you disagree with, rather than re-opening
  wholesale.

If yes to Q0–Q5 and no residual concern in Q6 and concurrence on Q7,
reply "approved — proceed with commit plan".  If no to any, say which
and why, and we iterate.

---

## References — file:line, verified 2026-04-20

**SEAS-MFEM:**
- `miniapps/seas/dynamic/fault_face_flux.cpp:39-65` — trial traction (Eq. 7)
- `miniapps/seas/dynamic/fault_face_flux.cpp:70-170` — Evaluate (Eqs. 8–12)
- `miniapps/seas/dynamic/fault_face_flux.cpp:145-164` — imposed-state construction (Eqs. 11–12)
- `miniapps/seas/dynamic/godunov_flux.cpp:125-188` — BuildJacobian (flux Jacobian A_x)
- `miniapps/seas/dynamic/godunov_flux.cpp:196-235` — BuildRotationInverse (global → face-aligned)
- `miniapps/seas/dynamic/godunov_flux.cpp:238-276` — BuildRotation (face-aligned → global)
- `miniapps/seas/dynamic/godunov_flux.cpp:305-320` — ApplySplitFlux (`A_x^+ Q_L + A_x^- Q_R`)
- `miniapps/seas/dynamic/godunov_flux.cpp:325-349` — Interior (full flux_.Interior call)
- `miniapps/seas/dynamic/wave_operator.inl:811-816` — interior-fault Step 3 (the bug site)
- `miniapps/seas/dynamic/wave_operator.inl:857-887` — interior-fault assembly (layer-2 bug)
- `miniapps/seas/dynamic/wave_operator.inl:1232-1244` — shared-fault Step 3 + assembly

**SeisSol (`/Users/chunhuizhao/projects/SeisSol/`):**
- `src/DynamicRupture/FrictionLaws/FrictionSolverCommon.h:143-222` — precomputeStressFromQInterpolated (friction solve)
- `src/DynamicRupture/FrictionLaws/FrictionSolverCommon.h:288-407` — postcomputeImposedStateFromNewStress (Eqs. 11–12)
- `src/Kernels/LinearCK/Neighbor.cpp:89-99` — DR case in the neighbouring-flux loop (main call site)
- `codegen/kernels/dynamic_rupture.py:87` — rotateFluxMatrix (`fluxSolver = fluxScale · A · T`)
- `codegen/kernels/dynamic_rupture.py:124-130` — nodalFluxGenerator (`I += M · Q · fluxSolver`)
- `src/Initializer/CellLocalMatrices.cpp:801-817` — fluxSolverPlus / fluxSolverMinus precompute

**Plan / context:**
- `debug_document/tpv102_debug_document/tpv102_debug_v9.0.0_debug_plan.md` §1, §2, §3, §A.3
- `debug_document/tpv102_debug_document/tpv102_debug_v9.0.0_fix.md` §1, §3

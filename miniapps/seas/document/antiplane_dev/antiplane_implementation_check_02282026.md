# Antiplane DG Implementation Check

**Date:** 2026-02-28
**Reference document:** *Discontinuous Galerkin For Sequences of Earthquakes and Aseismic Slip — MFEM Implementation* (Chunhui Zhao, February 2026)
**Scope:** Equation-by-equation comparison between the document's governing equations and the MFEM implementation in `miniapps/seas/`.

---

## 1. Governing Equations Reviewed

The document formulates the quasi-dynamic antiplane (Mode III) problem as a Poisson equation with a fault jump condition:

```
μΔu₃ = 0                              in Ω
σ₃₁n₁ + σ₃₂n₂ = T₃                   on Γᴺ
u₃ = u₃°                              on Γᴰ
R_{k3}(u₃⁺ − u₃⁻) = δₖ               on Γᶠ
```

This is discretized with a DG method (IP or BR2), producing:

- **Bilinear form** `B_h(u_h, v_h)` — assembled once and cached
- **Linear form** `L(v_h; δ)` — reassembled each time step from the current slip `δ`

The fault traction is extracted as the average normal flux and fed into the rate-and-state friction system.

---

## 2. Equation-by-Equation Comparison

### 2.1 Bilinear Form — Shear Modulus μ (Sections 1.5.1 & 1.5.2)

**Document (Eq. 28 BR2 / Eq. 36 IP):**
```
B_h(u,v) = ∫_Ω μ∇u·∇v dx − ∫_Γ μ([[u]]·{∇v·n} + {∇u·n}·[[v]]) ds + μ α(u,v)
```
μ multiplies every term.

**Implementation (`AssembleStiffness`, `domain/antiplane_operator.hpp:558`):**
```cpp
ConstantCoefficient one(1.0);
cached_a_->AddDomainIntegrator(new DiffusionIntegrator(one));
cached_a_->AddInteriorFaceIntegrator(
    new DGDiffusionIntegrator(one, sigma_, rep_kappa));
```

**Status — Design divergence (safe for constant μ).**
Both the stiffness matrix and the slip RHS use `μ = 1`. Since the linear system is `A·u = b` with the same implicit μ scaling on both sides, the solution `u_h` is identical to the system `(μA)·u = μb` when μ is constant. μ is then applied externally in `ComputeTraction`. This approach is correct for the current benchmarks (uniform shear modulus), but deviates from the document's formulation and would require a non-trivial refactor to support spatially varying μ.

---

### 2.2 Slip RHS — Shear Modulus μ (Document Eq. 31 / Eq. 39)

**Document (IP, L(v_h) fault terms, Eq. 39):**
```
L(v_h)|_{Γᶠ} = −Σ_{e∈Γᶠ} ∫_e {μ∇v_h·n} δ ds  +  Σ_{e∈Γᶠ} ∫_e μ(η_e/h_e) δ·[v_h] ds
```

**Implementation (`AssembleSlipContributionIP`, `antiplane_operator.hpp:1289-1296`):**
```cpp
elvec1(k) += sigma_ * dn1(k) * slip_imposed * w1;      // no μ
elvec1(k) += wq_penalty * slip_imposed * shape1(k);    // no μ
```

**Status — Design divergence (consistent with §2.1).**
The RHS is assembled without μ, consistent with the stiffness having μ = 1. Correct for constant μ; same caveat as §2.1.

---

### 2.3 Traction Extraction — Hardcoded x-Direction Normal (Algorithm Step 15)

**Document (Algorithm Step 15):**
```
τ_{qs,i} = { μ ∇u_h · n̂ }|_i
```
where `n̂` is the **unit fault normal** in the general sense.

**Implementation (`ComputeTraction`, `antiplane_operator.hpp:910-914`):**
```cpp
real_t avg_dudx = 0.5 * (grad1(0) + grad2(0));   // index 0 = x-component only
real_t tau_face = mu_ * avg_dudx;
```

**Status — Structural mismatch / hardcoded assumption.**
The implementation takes only `grad(0)` = ∂u/∂x, implicitly assuming the fault normal is exactly `n̂ = ê_x = (1, 0)`. The document's formula is `∇u · n̂` (dot product with the actual face normal vector). For the SEAS benchmarks with a perfectly vertical fault at x = 0, `ê_x` is the correct normal and the result is numerically exact. However, the code does not compute or use the face normal, making it incorrect for any fault geometry that is not perfectly axis-aligned.

The correct general implementation is:
```cpp
// Retrieve the face normal from the transformation
Vector nor(mesh_.Dimension());
CalcOrtho(FTr->Jacobian(), nor);
real_t nor_len = nor.Norml2();

// Dot gradient with unit normal from each side
real_t g1_dot_n = (grad1 * nor) / nor_len;   // ∂u⁺/∂n̂
real_t g2_dot_n = (grad2 * nor) / nor_len;   // ∂u⁻/∂n̂
real_t tau_face = mu_ * 0.5 * (g1_dot_n + g2_dot_n);
```

---

### 2.4 Traction Extraction — Single Midpoint vs. Face Quadrature (Algorithm Step 15)

**Document (Algorithm Step 15):** `τ_{qs,i}` is evaluated at fault node `i`. For a face-based DOF representation, this requires sampling the gradient at the representative face point.

**Implementation (`ComputeTraction`, `antiplane_operator.hpp:896-898`):**
```cpp
IntegrationPoint ip;
ip.x = 0.5;           // single midpoint — no quadrature loop
FTr->SetAllIntPoints(&ip);
```

**Status — Under-integrated for polynomial order p ≥ 2.**
The gradient within a DG element of order `p` is a polynomial of degree `p − 1`. For `p = 1`, the gradient is constant over the element, so midpoint evaluation is exact. For `p ≥ 2`, the gradient varies over the face and a single midpoint sample produces only an O(h) approximation to the face-averaged flux. The slip assembly in `AssembleSlipContributionIP/BR2` correctly uses a full quadrature loop over the face — `ComputeTraction` should use the same integration rule for consistency:

```cpp
int order = std::max(fe1->GetOrder(), fe2->GetOrder());
const IntegrationRule &ir = IntRules.Get(FTr->FaceGeom, 2*order + 1);
// ... quadrature loop over ir ...
```

---

### 2.5 Friction Coefficient f(V, θ) (Eq. 43)

**Document (Eq. 43):**
```
f(V,θ) = a · sinh⁻¹[ (V / 2V₀) · exp((f₀ + b · ln(V₀θ/Dᶜ)) / a) ]
```

**Implementation (`DieterichRuinaFriction::FrictionCoefficient`, `friction/dieterich_ruina.hpp:75-80`):**
```cpp
real_t log_arg  = cp_.V0 * theta / cp_.Dc;
real_t exp_arg  = (cp_.f0 + cp_.b * std::log(log_arg)) / a;
real_t sinh_arg = (V / (2.0 * cp_.V0)) * std::exp(exp_arg);
return a * std::asinh(sinh_arg);
```
**Status — Exact match.** ✓

---

### 2.6 State Evolution — Aging Law (Eq. 42)

**Document (Eq. 42):**
```
dθ/dt = 1 − V·θ / Dᶜ
```

**Implementation (`AgingLaw::Rate`, `friction/state_evolution.hpp:87-90`):**
```cpp
return 1.0 - V * theta / Dc;
```
**Status — Exact match.** ✓

---

### 2.7 Quasi-dynamic Stress Balance (Eq. 40–41)

**Document (Eq. 40–41):**
```
τ = F(V,θ) + ηV,     F(V,θ) = σₙ · f(V,θ),     η = μ / (2cₛ)
```

**Implementation (`SolveSlipRate`, `friction/dieterich_ruina.hpp:191-192`):**
```cpp
real_t F    = sigma_n * f + eta * V - tau;   // residual
real_t dF_dV = sigma_n * df_dV + eta;        // Jacobian
```
**Status — Exact match.** ✓

---

### 2.8 Newton-Raphson Convergence Criterion (Algorithm Steps 20–26)

**Document (Step 26):**
```
until  |V_new − V| / max(|V|, V_min) < 10⁻¹²   or   |F| < 10⁻¹² · τ
```

**Implementation (`dieterich_ruina.hpp:211-215`):**
```cpp
real_t rel_change = std::abs(V_new - V) / std::max(V, V_min_);
if (rel_change < tol || std::abs(F) < tol * tau) { break; }
```
**Status — Exact match.** ✓ (`tol = 1e-12`, `V_min_ = 1e-30`)

---

### 2.9 RK45 Error Norm (Algorithm Step 41)

**Document (Step 41):**
```
ε = max_i  |e_i| / (atol + rtol · |y_i^{n+1}|)     (weighted L-infinity)
```

**Implementation (`DormandPrinceRK45::Step`, `solver/time_stepper.hpp:364-375`):**
```cpp
real_t scale = atol_ + rtol_ * std::abs(y_tmp_(i));
real_t ei    = std::abs(err_(i)) / scale;
if (ei > err_norm) { err_norm = ei; }   // L-infinity over all DOFs
```
**Status — Exact match.** ✓
Note: `rtol_ = 1e-50` (effectively zero) to match Tandem/PETSc pure-absolute-tolerance mode.

---

### 2.10 RK45 Step Size Adaptation (Algorithm Step 43)

**Document (Step 43):**
```
Δt^{n+1} = min(Δt_max, max(Δt_min, γ · Δt^n · ε^{−1/5})),   γ = 0.9
```

**Implementation (`time_stepper.hpp:406-412`):**
```cpp
dt_new = safety_ * dt * std::pow(err_norm, -0.2);          // γ · dt · ε^{-1/5}
dt_new = std::max(dt_new, shrink_min_ * dt);               // per-step floor  (×0.1)
dt_new = std::min(dt_new, growth_max_ * dt);               // per-step ceiling (×10)
dt_new = std::max(dt_min_, std::min(dt_max_, dt_new));     // global bounds
```

**Status — Minor extension of the document formula.**
The per-step shrink/growth caps (`shrink_min_ = 0.1`, `growth_max_ = 10.0`) are PETSc defaults not present in the document. They add robustness without changing normal-operation behavior. Additionally, after a rejected step the implementation applies an extra `reject_safety_ = 0.5` (PETSc default), also not in the document.

---

### 2.11 Prestress τ₀ Computation (Algorithm Step 5)

**Document (Step 5):**
```
τ₀ = σₙ · aₘₐₓ · sinh⁻¹[(V_init/2V₀) · exp((f₀ + b·ln(V₀/V_init)) / aₘₐₓ)] + η · V_init
```

**Implementation (`BP2Params::tau0`, `config/bp2_params.hpp:172-178`):**
```cpp
real_t log_term = f0 + b * std::log(V0 / V_init);
real_t exp_term = std::exp(log_term / amax);
real_t arg      = (V_init / (2.0 * V0)) * exp_term;
real_t f_ss     = amax * std::asinh(arg);
return sigma_n * f_ss + eta() * V_init;
```
**Status — Exact match.** ✓

---

### 2.12 Initial State θ₀ (Algorithm Step 7)

**Document (Step 7):**
```
θ₀_i = (Dᶜ/V₀) · exp[(aᵢ/b)·ln(2V₀/V_init · sinh(f*_i/aᵢ)) − f₀/b]
    where  f*_i = (τ₀ + τ_{qs,i} − η_i · V_init) / σₙ
```

**Implementation (`DieterichRuinaFriction::InitialState`, `friction/dieterich_ruina.hpp:261-288`):**
```cpp
real_t tau_eff = tau0 - eta * V_init;         // (τ₀ + τ_qs,i) − η·V_init
real_t f       = tau_eff / sigma_n;           // = f*_i
// ...
real_t exp_arg = (a / cp_.b) * std::log(log_arg) - cp_.f0 / cp_.b;
real_t theta   = (cp_.Dc / cp_.V0) * std::exp(exp_arg);
```
Called from `RateStateFaultOperator::Init` with `tau = tau0_ + traction(i)`.
**Status — Exact match.** ✓

---

### 2.13 Seismogenic Zone Boundary Condition (Algorithm Steps 30–35)

**Document (Step 30–35):** For z_i < −W_f (below seismogenic zone):
```
k_s^{δ,i} = V_p   (prescribed plate rate)
k_s^{θ,i} = 0     (no state evolution)
```

**Implementation (`RateStateFaultOperator::ComputeRHS`, `fault/rate_state_fault.hpp:235-242`):**
```cpp
if (depths(i) < -params_.Wf)
{
   rate(i * StatePerNode + SlipIndex)  = params_.Vp;
   rate(i * StatePerNode + ThetaIndex) = 0.0;
   continue;
}
```
**Status — Exact match.** ✓

---

## 3. Summary Table

| # | Component | Document | Implementation | Status |
|---|-----------|----------|----------------|--------|
| 2.1 | Bilinear form μ scaling | μ in all terms | `one(1.0)` throughout | Design divergence — safe for const μ |
| 2.2 | Slip RHS μ scaling | μ in L(v_h) fault terms | no μ | Design divergence — consistent with 2.1 |
| **2.3** | **Traction normal direction** | **`{μ∇u·n̂}` with general n̂** | **`grad(0)` = ∂u/∂x only** | **Structural mismatch — hardcoded ê_x assumption** |
| **2.4** | **Traction quadrature** | **Face-node average flux** | **Single midpoint (ip.x=0.5)** | **Under-integrated for p ≥ 2** |
| 2.5 | Friction f(V,θ) | Eq. 43 | Exact | ✓ |
| 2.6 | Aging law dθ/dt | Eq. 42 | Exact | ✓ |
| 2.7 | Stress balance τ=F+ηV | Eq. 40–41 | Exact | ✓ |
| 2.8 | Newton convergence | Step 26 | Exact | ✓ |
| 2.9 | RK45 error norm | Step 41, L-inf | L-inf | ✓ |
| 2.10 | RK45 dt adaptation | Step 43, γ·dt·ε^{-1/5} | Same + PETSc per-step clips | Minor extension |
| 2.11 | Prestress τ₀ | Step 5 | Exact | ✓ |
| 2.12 | Initial θ₀ | Step 7 | Exact | ✓ |
| 2.13 | Below-Wf BC | Steps 30–35 | Exact | ✓ |

---

## 4. Prioritized Action Items

### Priority 1 — Structural (§2.3): Fix traction normal direction

The traction extraction should dot the gradient with the actual face normal rather than hardcoding the x-component. For the current benchmarks this is numerically correct (fault is exactly at x=0), but the code will silently produce wrong tractions for any non-axis-aligned fault geometry.

**Files:** `domain/antiplane_operator.hpp`, `ComputeTraction` (~line 879)

### Priority 2 — Accuracy (§2.4): Add quadrature loop to traction extraction

Replace the single midpoint evaluation with a proper integration rule consistent with the one used in `AssembleSlipContributionIP/BR2`. This is especially important for `order ≥ 2` runs where the gradient is not constant over the face.

**Files:** `domain/antiplane_operator.hpp`, `ComputeTraction` (~line 896)

### Priority 3 — Future-proofing (§2.1 & 2.2): Introduce μ coefficient

When non-uniform material properties are needed, the `DiffusionIntegrator` and `DGDiffusionIntegrator` in `AssembleStiffness`, and both `AssembleSlipContributionIP/BR2`, must use an actual `mu_` coefficient rather than `one(1.0)`.

**Files:** `domain/antiplane_operator.hpp`, `AssembleStiffness` (~line 553); `AssembleSlipContributionIP` (~line 1182); `AssembleSlipContributionBR2` (~line 1313)

---

## 5. Notes on Design Choices

- **μ = 1 in stiffness:** An intentional simplification valid for constant μ. The solution `u_h` is scale-invariant (dividing both sides of `Au = b` by μ gives the same u). Traction is correctly recovered as `μ · {∂u/∂x}` after the solve.
- **FSAL in RK45:** The FSAL (`k_[6]` reused as `k_[0]`) is implemented correctly with proper invalidation on step rejection.
- **Per-step dt clips:** `shrink_min = 0.1`, `growth_max = 10.0`, and `reject_safety = 0.5` are PETSc TSAdapt defaults adopted to match Tandem behavior; not specified in the document but well-established practice.

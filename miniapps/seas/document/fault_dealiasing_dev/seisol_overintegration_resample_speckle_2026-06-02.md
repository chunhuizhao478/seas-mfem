---
title: "Eliminating On-Fault Normal-Stress Speckle in Pure-Upwind DG Dynamic Rupture"
subtitle: "Why over-integration and slip-rate resampling (the SeisSol recipe) remove the speckle — mechanism, equations, and an implementation plan for the matrix/pure-upwind path (TPV31)"
author: "SEAS-MFEM dynamic-rupture team"
date: "2026-06-02"
---

# Purpose and scope

This document explains, from first principles and with equations, **why** the
two SeisSol techniques — **spatial over-integration of the fault flux** and
**slip-rate / state-variable resampling** — eliminate the on-fault normal-stress
"speckle" (and the secular $\sigma_n$ drift it drives) that our pure-upwind
discontinuous-Galerkin (DG) dynamic-rupture runs exhibit.

It is written to answer one question that has not yet been made concrete:

> *We can see that mixed flux removes the speckle and that SeisSol stays clean
> with pure upwind, but **how** do over-integration and resampling actually
> eliminate the speckle values we observe?*

The target application is **TPV31**, which is forced onto **pure upwind + ADER +
matrix (heterogeneous-material) Riemann** — the matrix Riemann solver is needed
for TPV31's depth-varying material, *not* a fault-normal material contrast (the
material is identical on both sides at each depth) — and can use *neither* mixed
flux *nor* the RK integrator. The two techniques here are **flux-agnostic and ADER-compatible**:
they change *how the nonlinear friction is represented*, not the flux or the
eigenvalues, so they apply unchanged to that path. The same fix benefits every
pure-upwind spatial-driver run (TPV102/104/205, SAFS).

This is an explanation + implementation plan for review. **No source code is
changed by this document.**

---

# 1. The symptom, and what we have already ruled out

On a planar fault with material identical on both sides (TPV31, TPV102/104/205),
the on-fault normal stress $\sigma_n$ is **analytically constant** for a
strike-slip rupture. In our pure-upwind runs it instead shows:

- a high-spatial-frequency **speckle** across the fault plane (the ParaView
  fault-surface field, amplitude comparable to $\sigma_{n0}$), and
- a secular **$\sigma_n$ collapse** at every station (TPV102: $120 \to 80,\,-30,\,75$ MPa;
  TPV31: $60.6 \to 23$ MPa, then blow-up).

Established by controlled experiment and code reading (see the companion
findings in `debug_document/tpv31_debug_document/`):

| observation | conclusion |
|---|---|
| homogeneous TPV102/104/205 (scalar Riemann) drift identically | **not** the heterogeneous-material / matrix path |
| mixed flux (touches only fault-*adjacent bulk* faces) removes it | the fault-face friction formula is fine; the seed is near-fault |
| serial $\equiv$ parallel to $10^{-14}$ (`test_rupture_tilted_fault_serial_vs_parallel`, $\theta=0$) | **intrinsic**, present in serial; **not** MPI / partition |
| planar fault $\Rightarrow |n\cdot\hat n_{\mathrm{ref}}|=1$, frame exact | **not** the orientation/frame bug (that is curvilinear-only) |
| p1 $\approx$ p2 (erratic, non-monotone); clean mixed-flux p1 $\approx$ p2 | **not** "just raise the order" — order does not converge it |

So the defect is **intrinsic to the discretization, present at p1 in serial, in
the common (scalar + matrix) pure-upwind path, and order-independent.** This
document identifies that defect as **aliasing of the nonlinear friction law**,
and explains why SeisSol's over-integration + resample remove it.

---

# 2. The one equation both codes share, and why $\sigma_n$ must be constant

At each fault evaluation point, in fault-aligned coordinates $(n, t_1, t_2)$,
both our code (`fault_face_flux.cpp:62`) and SeisSol
(`FrictionSolverCommon.h:182`) compute the identical Godunov interface state:

$$
\sigma_n^{*} \;=\; \eta_p\Big(\,\underbrace{v_n^{-}-v_n^{+}}_{[\![ v_n]\!]}
\;+\; \frac{\sigma_n^{+}}{Z_p^{+}}+\frac{\sigma_n^{-}}{Z_p^{-}}\Big),
\qquad
\eta_p=\frac{Z_p^{+}Z_p^{-}}{Z_p^{+}+Z_p^{-}},\quad Z_p=\rho c_p .
$$

$$
\tau_i^{*} \;=\; \eta_s\Big([\![ v_{t_i}]\!]+\frac{\tau_i^{+}}{Z_s^{+}}+\frac{\tau_i^{-}}{Z_s^{-}}\Big),
\quad i=1,2 .
$$

The friction law then solves for a slip rate $S$ from the shear balance and
returns the **corrected (radiated) shear traction**

$$
\tau_i^{\mathrm{corr}} \;=\; \tau_i^{*}-\eta_s\,S_i ,
$$

leaving the normal traction at its Godunov value. Writing
$\sigma_n^{\mathrm{tot}}=\sigma_{n0}+\sigma_n^{*}$ (perturbation about the
background), the normal-stress perturbation is

$$
\boxed{\;\delta\sigma_n \;=\; \eta_p\,[\![ v_n]\!] \;+\;
\eta_p\Big(\tfrac{\sigma_n^{+}}{Z_p^{+}}+\tfrac{\sigma_n^{-}}{Z_p^{-}}\Big)\;}
$$

For a planar strike-slip fault with symmetric material there is **no
normal-direction wave**: the physical $[\![ v_n]\!]=0$ and the normal-stress
fluctuation is zero, so $\delta\sigma_n\equiv 0$. **Any** spurious normal-velocity
jump is amplified by $\eta_p \approx Z_p/2 \approx 8\times10^{6}\ \mathrm{Pa\,s/m}$:

$$
\delta\sigma_n \;\approx\; \eta_p\,[\![ v_n]\!]
\quad\Longrightarrow\quad
[\![ v_n]\!] \approx \frac{40\ \mathrm{MPa}}{8\times10^{6}} \approx 5\ \mathrm{m/s}.
$$

The observed $\approx$40 MPa drift therefore corresponds to a **coherent**
$[\![ v_n]\!]\approx 5$ m/s — *of the same order as the slip rate itself*. This
is not faint round-off grit; it is a large, systematic leak of shear motion into
the normal channel. The rest of this document explains where that coherent
$[\![ v_n]\!]$ comes from.

---

# 3. Root cause: aliasing of the nonlinear friction

## 3.1 The friction law is the only nonlinearity, and it is not a polynomial

In a DG scheme the bulk operator for linear elasticity is **linear**: on affine
elements with exact quadrature, the nodal and modal formulations are *identical*
(a change of basis). The **only** nonlinearity in the entire dynamic-rupture
operator is the fault friction. At a fault point it returns

$$
S \;=\; G\big(\tau^{*},\,\sigma_n^{*},\,\psi,\,\dots\big),
\qquad
\tau^{\mathrm{corr}} \;=\; \tau^{*}-\eta_s S ,
$$

where $G$ is rate-and-state (an $\operatorname{asinh}$ form), or linear
slip-weakening with a **kink** at $\delta=d_c$, or a rupture front where $S$
jumps from $0$ to several m/s over one element. **None of these is a
polynomial.** Even when the inputs $\tau^{*}(x),\sigma_n^{*}(x),\psi(x)$ are
degree-$N$ polynomials on the face, the output $S(x)$ has spectral content at
**all** wavenumbers, including far above degree $N$.

## 3.2 What "aliasing" means here, precisely

A DG scheme can only carry the degree-$N$ projection of $S$. The $L^2$ projection
onto the degree-$N$ face polynomial space is

$$
\widehat S \;=\; \Pi_N S \;=\; M^{-1}\!\int_{F}\!\boldsymbol\phi\,S(x)\,dx ,
\qquad M_{ij}=\int_F \phi_i\phi_j\,dx ,
$$

with $\boldsymbol\phi=\{\phi_i\}$ the basis. The integral is evaluated by
**quadrature**,

$$
\int_F \phi_i\,S\,dx \;\approx\; \sum_{q} w_q\,\phi_i(x_q)\,S(x_q).
$$

Because $S$ is **not** a polynomial, *no finite quadrature is exact*. A rule that
is exact only to degree $2N$ (our minimal "mass-matrix order", `2*order`) cannot
distinguish the high-wavenumber content of $S$ (modes above $N$) from the
resolved low modes: that energy is **folded down — aliased — onto the
coefficients $\widehat S_0,\dots,\widehat S_N$.** This is aliasing error. It is
*not* reduced by mesh refinement of a fixed feature and it is *not* zero-mean in
the quantities that matter (Section 3.4).

## 3.3 Two channels: per-step flux aliasing, and secular state accumulation

The friction nonlinearity pollutes the solution through **two distinct
channels**, and the two SeisSol techniques each close one of them.

**(1) Per-step flux aliasing (instantaneous).** Each step the radiated traction
is assembled as the flux integral $\int_F\phi_i\,\tau^{\mathrm{corr}}\,dx$,
evaluated by quadrature. Because $\tau^{\mathrm{corr}}$ is non-polynomial
(Section 3.1), an under-integrated rule **folds its $>N$ content onto the
resolved coefficients $0..N$** (Section 3.2). The fold-over hits the **top
coefficient $c_N$ first** — it has the lowest fold-over threshold (Section 4.1) —
so the spurious energy appears as the most oscillatory representable mode,
degree $N$: *the per-element checkerboard you see as speckle.* **Over-integration
closes this channel:** with enough Gauss points the resolved coefficients are the
true $L^2$ values and no spurious top-mode energy is injected. This keeps the
solution's *genuine* degree-$N$ content — it removes the *spurious* part, it does
**not** delete the mode.

**(2) Secular state accumulation.** Slip and state are advanced every step by
integrating the *pointwise* friction output at the GPs. Those pointwise values
carry $>N$ content that the degree-$N$ DOF field cannot represent; left in, it
**accumulates in the slip/state field over many steps**, drifting it out of the
polynomial space and re-seeding speckle through the strength $G(\cdot)$.
**Resampling closes this channel:** it projects the per-step slip/state increment
back onto the degree-$N$ DOF space (Section 4.2), discarding only the
unrepresentable $>N$ part and keeping all modes $0..N$ (including mode $N$).

The two are **coupled, not independent**: the resample projector is the
*identity* whenever the GP count equals the DOF count (the minimal mass-matrix
rule), so it does nothing **until** the flux is over-integrated. Over-integration
is the prerequisite; resample is the secular cleanup layered on top.

## 3.4 The chain from friction aliasing to the $\sigma_n$ speckle and drift

Here is the causal chain the user asked to see spelled out:

$$
\underbrace{S(x),\ \tau^{\mathrm{corr}}(x)}_{\text{nonlinear friction}}
\;\xrightarrow[\text{no dealiasing}]{\ \Pi_N\ }\;
\underbrace{\text{spurious degree-}N\text{ top mode in }\tau^{\mathrm{corr}}}_{\text{checkerboard on the fault}}
\;\xrightarrow[\text{radiates}]{\ \text{bulk}\ }\;
\underbrace{\text{checkerboard wavefield, incl. } v_n}_{\text{spurious normal P-radiation}}
\;\xrightarrow{\ \eta_p\ }\;
\underbrace{\delta\sigma_n\ \text{speckle}}_{\times\,\eta_p\ \approx\ 8\times10^{6}} .
$$

1. The friction reaction $\tau^{\mathrm{corr}}=\tau^{*}-\eta_s S$ is the traction
   the fault **applies back to the bulk**. If $\tau^{\mathrm{corr}}$ carries the
   spurious top mode, the fault **radiates a checkerboard wavefield**.
2. A checkerboard shear source on an imperfect (non-mirror-symmetric) near-fault
   mesh does **not** radiate purely tangentially: it leaks a fault-normal
   particle-velocity component $v_n$ into the bulk next to the fault.
3. That $v_n$ returns to the fault as $[\![ v_n]\!]\neq0$ and is multiplied by
   $\eta_p\approx 8\times10^{6}$ in Eq. (box, Section 2) — turning a small slip-rate
   ringing into a **tens-of-MPa** $\sigma_n$ speckle.

**Why a zero-mean speckle becomes a one-way drift.** The friction *strength* uses
$|\sigma_n|$ (rate-state) or $\max(\sigma_n,0)$ (LSW). These are **even /
rectifying** functions. A symmetric $\pm\delta\sigma_n$ speckle therefore
produces an **asymmetric** change in strength, hence in slip, hence a **net
(secular) bias** in $\sigma_n$ — the observed collapse. For LSW, once a tensile
excursion crosses $\sigma_n\le0$ the strength drops to cohesion, the fault
free-slides, and the bias runs away (TPV31 blow-up). **Remove the speckle and
both the speckle and its rectified drift disappear.**

## 3.5 Why raising the order does not help (and can hurt)

Raising $N$ makes the rupture front *sharper relative to the element* and pushes
the aliased energy onto a higher-wavenumber top coefficient $c_N$. With no
degree-dependent dealiasing, more order folds *more* $>N$ energy onto the resolved
modes — exactly the erratic, non-converging p1$\to$p2 behaviour we measured.
SeisSol's cure is **order-aware**: the fault over-integration scales with the
order ($(N+2)^2$ GPs) and the resample projects onto the *current* degree-$N$
space, so aliasing is controlled at every order and high-order runs pay off; ours
has neither knob, so nothing tracks the order.

---

# 4. SeisSol's two mechanisms (with equations)

SeisSol uses pure-upwind Godunov flux **everywhere** (no mixed flux) and stays
clean because of *how it forms and represents the fault solution*:

## 4.1 Over-integration of the fault flux

SeisSol integrates the fault flux at $(N+2)^2$ Gauss points per triangular face
(Pelties et al. 2014, §3 — "these flux functions are integrated with a quadrature
based on $(N+2)^2$ GPs"), and assigns friction parameters per Gauss point. This
is **more** than the minimal mass-matrix rule. The extra points capture more of
the spectrum of the nonlinear integrand, so the aliasing of high modes onto the
resolved coefficients (Section 3.2) is strongly suppressed:

$$
\Big|\textstyle\sum_q w_q\phi_i(x_q)S(x_q)-\int_F\phi_i S\,dx\Big|
\ \xrightarrow[\ \text{more GPs}\ ]{}\ \text{small}.
$$

Our path uses `IntRules.Get(geom, 2*order)` (`spatial_dyn_driver.cpp:217`) —
the minimal mass-matrix order, **under-integrated for the nonlinear friction.**

### Worked example: how under-integration makes a high mode masquerade as a resolved one

This is the single most important fact, so here it is in numbers. Work on the
reference interval $[-1,1]$ with the Legendre basis of a degree-$N=2$ element,
$$P_0=1,\ P_1=x,\ P_2=\tfrac{3x^2-1}{2},\ P_3=\tfrac{5x^3-3x}{2},\
P_4=\tfrac{35x^4-30x^2+3}{8},$$
whose defining property is orthogonality — *different modes integrate to zero
against one another*:
$$\int_{-1}^{1} P_m(x)\,P_n(x)\,dx = 0,\qquad m\neq n .$$
A degree-2 element is therefore **blind** to any higher mode: the projection of
$P_3,P_4,\dots$ onto $\{P_0,P_1,P_2\}$ is exactly zero.

**What our actual rule is.** Our fault quadrature is `IntRules.Get(geom,2*order)`,
a rule exact to degree $2N$. For $N=2$ MFEM returns the **3-point Gauss rule**
(nodes $0,\pm\sqrt{3/5}$, weights $\tfrac89,\tfrac59,\tfrac59$; exact to degree
$5$). The contaminated coefficient that matters is the **top** one,
$c_2=\int_{-1}^1 P_2\,S\,dx$ — the fold-over reaches it first.

**The next mode up, $P_3$, does NOT alias under our rule.** With $S=P_3$ the
product $P_2P_3$ has degree $5$, which the 3-point rule integrates *exactly*, so
$$c_2^{\,(3\text{-pt})}=\int_{-1}^1 P_2\,P_3\,dx=0\qquad\checkmark$$
— the true value, by orthogonality. (So, unlike a naive degree-$2N{-}1$ rule, our
actual rule is already clean for $P_3$; the first mode that *bites* is one
higher.)

**Put the first biting mode into the slip rate:** $S(x)=P_4(x)$. The true top
coefficient is zero ($P_4\perp P_2$): $c_2=\int_{-1}^1 P_2\,P_4\,dx=0$.
**Now compute it with our 3-point rule.** The integrand $P_2P_4$ has degree
$6>5$, so the rule is *not* exact. At $x=\pm\sqrt{3/5}$ ($x^2=\tfrac35$):
$P_2=\tfrac25$, $P_4=-\tfrac{3}{10}$, so $P_2P_4=-\tfrac{3}{25}$; at $x=0$:
$P_2=-\tfrac12$, $P_4=\tfrac38$, so $P_2P_4=-\tfrac{3}{16}$. Hence
$$c_2^{\,(3\text{-pt})}=2\cdot\tfrac59\Big(\!-\tfrac{3}{25}\Big)+\tfrac89\Big(\!-\tfrac{3}{16}\Big)
=-\tfrac{2}{15}-\tfrac16=\boxed{-\tfrac{3}{10}}\;\neq\;0 .$$
The unresolved $P_4$ mode has **injected $-\tfrac{3}{10}$ into the resolved $P_2$
slot.** The 3-point rule cannot tell $P_4$ from $P_2$ at its three nodes and
charges $P_4$'s integral to $P_2$. *That nonzero $c_2$ is the spurious top-mode
wiggle — the speckle.*

**Over-integrate** with a 4-point Gauss rule (exact to degree $7\ge6$): the
degree-6 product is now integrated exactly, so
$$c_2^{\,(4\text{-pt})}=0 \qquad\checkmark$$
and $P_4$ again contributes nothing to the resolved modes. Over-integration did
**not smooth** $P_4$ ($P_4$ is untouched) — it stopped $P_4$ from *masquerading*
as $P_2$. In one line: $-\tfrac{3}{10}$ (contaminated) $\to 0$ (clean).

**The general rule.** A quadrature exact to degree $D$ evaluates
$\int_{-1}^{1}P_k\,S\,dx$ correctly *only* for the part of $S$ of degree
$\le D-k$; anything above that **folds onto** the coefficient $c_k$. The
$(N{+}1)$-point Gauss rule that `IntRules.Get(geom,2*order)` returns is exact to
degree $D=2N+1$ (Gauss over-delivers by one beyond the requested $2N$), so $c_N$
is correct up to degree $D-N=N+1$ and the first mode that folds onto it is
$\deg(P_N\,P_{N+2})=2N+2>2N+1$, i.e. **$P_{N+2}$** (for $N=2$: $P_4$, as above).
Over-integration raises $D$ and pushes the fold-over threshold above the modes
that bite. *(On simplex faces MFEM may deliver exactly degree $2N$ rather than
$2N+1$, in which case the first folding mode is $P_{N+1}$; tie the bound to the
rule actually returned for the face geometry in use.)*

**Why it is called "folding" (the sampling view).** On $M$ equally-spaced points
$x_j=2\pi j/M$,
$$e^{\,i(k+M)x_j}=e^{\,ikx_j}\,\underbrace{e^{\,i2\pi j}}_{=\,1}=e^{\,ikx_j},$$
so mode $k+M$ is *identical to mode $k$ at the samples*, and the measured
coefficient is the true one plus all of its aliases:
$$\widehat S_k^{\,(\text{sampled})}=\widehat S_k+\widehat S_{k+M}+\widehat S_{k-M}+\cdots$$
More points $M\Rightarrow$ the aliases $k\pm M$ sit farther away $\Rightarrow$ the
resolved band stays clean. ("A 1 kHz tone read as 100 Hz" is exactly
$e^{\,i(k+M)x_j}=e^{\,ikx_j}$.)

This $c_2=-\tfrac{3}{10}$ is precisely the spurious top-mode oscillation in the
flux that, via $\tau^{\mathrm{corr}}=\tau^{*}-\eta_s S$, radiates the checkerboard
and feeds $\delta\sigma_n=\eta_p\,(v_n^{-}-v_n^{+})$. **Over-integration removes
this contamination of the resolved modes — it is the per-step lever on the
speckle.** **Resampling** (next) plays the *complementary, secular* role: it
projects the *accumulated* slip/state increment back onto the degree-$N$ space so
the unrepresentable $>N$ content does not build up over many steps. Resampling
**keeps** the top mode $c_N$ (the genuine degree-$N$ content); it is the
over-integration, not the resample, that makes the per-step traction clean.

## 4.2 Resampling: keep the accumulated slip/state in the degree-$N$ space

After the pointwise friction solve, SeisSol **resamples** the per-step increment
of the accumulated state — the **slip-rate magnitude** for LSW
(`LinearSlipWeakening.cpp:18`), the **state-variable increment** $\Delta\psi$ for
rate-state (`FastVelocityWeakeningLaw.h:144`, `RateAndState.h:89`) — onto the DOF
polynomial space and re-evaluates at the (over-integrated) fault GPs
(matrix doc `LinearSlipWeakening.h:227`):

> "resampleMatrix first projects LocSR on the two-dimensional basis on the
> reference triangle with degree $\le$ ConvergenceOrder $-1$, and then evaluates
> the polynomial at the quadrature points."

**This is the full degree-$N$ space, not $N-1$.** In SeisSol the polynomial
degree equals ConvergenceOrder $-1$ (`FreeSurfaceIntegrator.h:42`:
`PolyDegree = ConvergenceOrder - 1`), so "degree $\le$ ConvergenceOrder $-1$"
means "degree $\le N$" — the *same* space as the DOFs. The comment states the
goal exactly: that "the state increment (slip) lies in the **same polynomial
space as the degrees of freedom**." Numerically, the generated resample matrices
have rank $(N+1)(N+2)/2$ — the full degree-$N$ DOF count — at every order
(verified for $N=1..6$).

Formally, with $V_{N}$ the Vandermonde from degree-$N$ modal coefficients to the
(over-integrated) GP values and $W=\mathrm{diag}(w_q)$ the quadrature weights,

$$
\boxed{\;R \;=\; V_{N}\,\big(V_{N}^{\mathsf T}W\,V_{N}\big)^{-1}V_{N}^{\mathsf T}W\;}
\qquad
\Delta \leftarrow R\,\Delta ,
$$

applied to the increment $\Delta$ (slip-rate magnitude, or $\Delta\psi$). $R$ is
the **$L^2$ projector onto the degree-$N$ DOF space, expressed in GP-value
space**: it **reproduces all modes $0..N$ exactly (the top mode $N$ included)**
and **discards only the $>N$ content** that the over-integrated GPs expose but the
DOF basis cannot represent. Two consequences:

- **$R$ is the identity when \#GP $=$ \#DOF** (the minimal mass-matrix rule, e.g.
  tensor Gauss on a quad face: $(N+1)^2$ GPs $=$ $(N+1)^2$ DOFs). It does nothing
  until the flux is over-integrated (Section 4.1) — the two techniques are
  **coupled**.
- **$R$ does not touch the radiated traction.** SeisSol builds
  $\tau^{\mathrm{corr}}=\tau^{*}-\eta_s S$ (and the directional slip) from the
  *un-resampled* friction-solve rate inside `calcSlipRateAndTraction`, *before*
  the resample; resample is applied only to the accumulated slip/state that drives
  the strength on the next step. So $R$ controls the **secular** out-of-space
  drift of the state, while **over-integration** controls the **per-step**
  radiated speckle (Sections 4.1 and 3.3).

**A caveat on naming.** SeisSol's `BiMaterialFault` specialization
(Prakash–Clifton) is a true **across-fault material contrast**, and there SeisSol
**disables** resampling because in that regime "resampling introduces artificial
oscillations" (`LinearSlipWeakening.h:261`). **TPV31 is *not* this case** — its
material is identical on both sides at each depth (our "matrix" Riemann handles
depth-heterogeneity, not a fault-normal jump), so resampling is appropriate and
should be **on** for TPV31/TPV102/104/205. The action item is the converse: on
any target with a genuine across-fault material contrast (some SAFS
configurations), resample must be **disabled**, matching SeisSol.

## 4.3 A worked 1-D illustration (what each technique actually does)

Take a single fault element as $[-1,1]$, degree $N=2$, Legendre basis
$P_0,P_1,P_2$. Let the friction response be a sharp rupture front,
$S(x)=\tfrac12\big(1+\tanh(x/\varepsilon)\big)$ with $\varepsilon=0.1$. As a
non-polynomial it has content at **all** degrees:
$S=a_0P_0+a_1P_1+a_2P_2+(\text{degree}>2\ \text{tail})$.

- **Under-integrated flux (minimal GPs):** the $>2$ tail folds onto the resolved
  coefficients (Section 4.1), worst onto the top coefficient $a_2$, so the
  per-step flux integral $\int_F\phi_i\,\tau^{\mathrm{corr}}$ is **wrong** —
  carrying a spurious, maximally-oscillatory degree-2 component. *Across many
  elements these spurious top modes tile into the checkerboard speckle.*
- **Over-integration (the per-step fix):** with enough GPs the resolved
  $a_0,a_1,a_2$ are the **true** $L^2$ coefficients — no folded-in energy. The
  radiated traction is the honest degree-2 projection of the friction, so it no
  longer carries a *spurious* top mode. (Its genuine $a_2$ — the true degree-2
  content — is **kept**; we remove the aliasing, not the mode.)
- **Resample $R$ (the secular fix):** the slip/state is advanced by integrating
  the *pointwise* friction at the over-integrated GPs, which still carries the
  $>2$ tail. $R$ projects that increment back onto degree 2
  ($R\widehat\Delta=a_0P_0+a_1P_1+a_2P_2$, top mode **kept**), discarding only the
  unrepresentable $>2$ tail. Without it the tail accumulates in the slip/state
  over many steps and re-seeds speckle through the strength $G(\cdot)$; with it
  the state stays a clean degree-2 field.

Over-integration makes the per-step radiated traction the true degree-2
projection (no spurious top mode) $\Rightarrow$ no checkerboard radiation
$\Rightarrow$ $[\![ v_n]\!]\to$ round-off $\Rightarrow$
$\delta\sigma_n=\eta_p[\![ v_n]\!]\to 0$; resample keeps the *accumulated* state in
the degree-2 space so the drift does not creep back. **That is how
over-integration (per step) + resample (secular) eliminate the speckle — and why
resample is a no-op without the over-integration underneath it.**

## 4.4 Where SeisSol applies friction in time (context, not the lever)

SeisSol solves friction at the $N$ Gauss–Legendre *time* nodes of one ADER step
and integrates with `timeWeights` to a single imposed state
(`FrictionSolverCommon.h:329`; Uphoff thesis Eq. 4.60), with one slip rate per
(space-GP, time-node) used for *both* the flux and the slip increment. Our
isolation run (pure-upwind + RK45, job 7762948) showed the **time integrator is
not the lever** — RK45 drifts the same as ADER. The spatial dealiasing (4.1, 4.2)
is the lever.

---

# 5. Our current code vs. the target (the exact gap)

| ingredient | SeisSol | SEAS-MFEM (pure-upwind ADER, spatial driver) |
|---|---|---|
| fault-flux quadrature | $(N+2)^2$ GPs (over-integrated) | `2*order` (mass-matrix minimal) — `spatial_dyn_driver.cpp:217` |
| slip / state dealiasing | $R$: project the per-step increment to degree $N$ (full DOF space) on the over-integrated GPs | **none** (`grep resample` = 0 hits) |
| resample target | LSW: slip-rate magnitude $\to$ accumulated slip; RS: state-variable increment $\Delta\psi$ | **none** |
| fault state storage | per-GP, kept in a degree-$N$ field via $R$ | per-GP, **independent** (`DOFData`, no projection) |

The two missing ingredients are exactly the dealiasing pair. Crucially, **both
are flux-agnostic and do not alter the eigenvalues**, so they are compatible with
ADER (no RK needed) and with the heterogeneous-material **matrix** Riemann path —
i.e. with TPV31's hard constraints. (Note this "matrix" path is for
depth-heterogeneity, *not* a fault-normal material contrast; see the
`BiMaterialFault` caveat in Section 4.2.)

---

# 6. Implementation plan (for `/code-implement` after review)

> Goal: add (A) over-integration of the fault flux and (B) slip/state resampling
> to the pure-upwind fault path — **in that order** (resample is a no-op without
> over-integration, Section 4.2) — keeping byte-exact behaviour when both knobs
> are off, and validate that they flatten $\sigma_n$ without mixed flux or RK.

## Phase 0 — Local reproduction harness (no production runs)
- Build a **serial**, homogeneous, **planar-fault** mini-rupture unit test large
  enough for a front to propagate (e.g. a structured slab of $\approx$20–40 fault
  faces; the existing 2-tet fixture is too small to develop the drift).
- Outputs: $\max_F|\sigma_n-\sigma_{n0}|(t)$ and $\max_F|[\![ v_n]\!]|(t)$.
- Acceptance: reproduces a growing $\sigma_n$ excursion at p1 in serial (confirms
  Sections 1–3 locally and gives an A/B testbed). Run a deliberately
  **asymmetric** vs **mirror-symmetric** mesh variant to confirm the near-fault
  mesh-asymmetry sensitivity.

## Phase 1 — Over-integration of the fault flux (the prerequisite)
- Add a separate, higher fault-quadrature degree (target $(N+2)^2$ GPs per face,
  SeisSol's choice, or a tunable `fault_overint_factor`) used for the friction
  solve **and** the flux-assembly integral $\int_F\phi_i F_h$, decoupled from the
  bulk mass-matrix rule.
- Gate behind `--fault-overint <k>` (default off $\Rightarrow$ falls back to
  `2*order` $\Rightarrow$ byte-exact with today).
- Assign friction parameters per over-integration GP. This is the **per-step
  speckle lever** (Section 4.1), and it must land **before** resample, because
  the resample projector is the identity unless \#GP $>$ \#DOF.

## Phase 2 — Resample operator $R$ on the fault face
- For the reference fault-face geometry, order $N$, and the **over-integrated** GP
  set from Phase 1, build $R$ (Eq. 4.2) once: GP-value $\to$ **degree-$N$** $L^2$
  projection $\to$ GP-value. Affine faces share one reference $R$ (the $|J_F|$
  scale cancels in the projector).
- New helper, e.g. `dynamic/fault_resample.hpp` (mirrors SeisSol's
  `init::resample`), with a unit test asserting:
  (i) $R$ reproduces **any degree-$N$ field exactly — including the pure
  degree-$N$ top mode** (NOT annihilates it);
  (ii) $R$ **removes degree-$>N$ content** sampled at the over-integration GPs
  (the part orthogonal to the degree-$N$ space goes to zero; the degree-$N$ best
  fit remains);
  (iii) $R$ is idempotent ($R^2=R$) with **rank $(N+1)(N+2)/2$** (full degree-$N$
  DOF count — matches SeisSol's resample matrices, verified $N=1..6$);
  (iv) at the minimal rule (\#GP $=$ \#DOF), $R=I$ (so Phase 2 alone is
  byte-exact).

## Phase 3 — Apply $R$ in the friction iterators
- In `friction_substep_iterator.{hpp,cpp}` (and the RK path `rk_time_stepper`):
  build $\tau^{\mathrm{corr}}$ and the directional slip from the **un-resampled**
  friction-solve rate (as today), then apply $R$ **only to the accumulated state
  increment**, per face — matching SeisSol (`calcSlipRateAndTraction` runs before
  the resample). **Do not rebuild $\tau^{\mathrm{corr}}$ from a resampled
  quantity.**
- **Resample target depends on the law** (Section 4.2): **LSW** $\to$ resample the
  slip-rate magnitude and integrate it into the accumulated slip (which drives the
  slip-weakening $\mu$); **rate-state** $\to$ resample the state-variable increment
  $\Delta\psi$.
- This needs a **per-face gather/scatter** of the per-QP `DOFData`; the QPs of one
  face form the polynomial $R$ acts on. The data are already laid out in per-face
  contiguous QP blocks of size `nbf = IntRules.Get(geom, 2*order).GetNPoints()`
  (`spatial_dyn_driver.cpp:217`); over-integration grows `nbf`. Assert
  QP-count-per-face $=$ rule size.
- Gate behind `--fault-resample` (default off). With Phase 1 also off, $R=I$ and
  the path is byte-exact with today.

```cpp
// Phase 3 sketch. tau_corr + directional slip come from the (un-resampled)
// friction-solve rate, exactly as today. Resample acts ONLY on the accumulated
// state increment, and only does anything when the GP set is over-integrated
// (Phase 1); off => byte-exact. DOFData::V1/V2 = friction-solve slip rate; the
// *_trial values live in the per-QP working state (fault_face_flux.hpp:153).
for (const Face &f : fault_faces) {
    for (int q = 0; q < nbf; ++q) {                  // nbf = QPs per face (grows w/ over-int)
        DOFData &d = dof[f.qp0 + q];
        d.slip1    += d.V1 * dt_sub;                 // un-resampled rate -> directional slip
        d.slip2    += d.V2 * dt_sub;
        d.tau1_corr = work[q].tau1_trial - d.eta_s * d.V1;  // un-resampled rate -> flux
        d.tau2_corr = work[q].tau2_trial - d.eta_s * d.V2;
    }
    if (fault_resample) {                            // dealias the ACCUMULATED state increment
        if (law == LSW) {                            // LSW: |V|*dt -> accumulated slip
            GatherFaceState(f, dSlipMag);
            Apply(R, out, dSlipMag);                 // out = R * dSlipMag (separate buffer)
            ScatterAddAccumSlip(f, out);
        } else {                                     // rate-state: state-variable increment
            GatherFaceState(f, dPsi);
            Apply(R, out, dPsi);
            ScatterAddState(f, out);
        }
    }
}
```

## Phase 4 — Wire into the spatial driver + matrix path
- Expose both knobs in `spatial_dyn_driver.cpp`; confirm they are honoured on the
  `interior_flux="matrix"` (`BimaterialWaveOperator`) path — these touch only the
  fault-flux assembly/iterator, not the interior Riemann, so no conflict with the
  matrix path or with ADER.
- Resample is **on** for TPV31/TPV102/104/205 (symmetric material). Add a guard
  that **disables resample on any genuine across-fault material contrast**
  (matching SeisSol's `BiMaterialFault`, Section 4.2).

## Acceptance criteria
- Both knobs off $\Rightarrow$ byte-exact vs current (regression contract
  preserved). `--fault-resample` alone (no over-integration) $\Rightarrow$ also
  byte-exact, because $R=I$ (Phase 2 test iv).
- On (Phase 0 harness, p1, serial, **`--fault-overint` on**, then
  `--fault-overint --fault-resample`): $\max_F|\sigma_n-\sigma_{n0}|$ drops by
  $\ge 10\times$ and $\max_F|[\![ v_n]\!]|\to$ round-off, with over-integration
  carrying the per-step reduction and resample suppressing the residual secular
  creep.
- On (p1 vs p2, both knobs on): the residual $\sigma_n$ error now **decreases**
  with order.

---

# 7. Validation plan

1. **Local A/B (Phase 0 harness):** off vs `--fault-overint` vs
   `--fault-overint --fault-resample`; report $\sigma_n$/$[\![ v_n]\!]$ curves and
   the p1$\to$p2 convergence.

```bash
# Phase 0 local A/B — serial, planar fault, no production run (compile + unit test)
make seas_test_fault_planar_serial
./seas_test_fault_planar_serial --order 1                                     # baseline: sigma_n drifts
./seas_test_fault_planar_serial --order 1 --fault-overint 2                   # per-step speckle down
./seas_test_fault_planar_serial --order 1 --fault-overint 2 --fault-resample  # + secular state cleanup
./seas_test_fault_planar_serial --order 2 --fault-overint 2 --fault-resample  # order now helps
# NB: --fault-resample ALONE (no over-integration) is a NO-OP: R = identity when
#     #GP == #DOF (Section 4.2), so it must be combined with --fault-overint.
```
2. **Frontera confirm (cheap, flex):** TPV102 pure-upwind ADER with the knobs on
   — predict $\sigma_n$ flattens toward 120 MPa like the mixed-flux+RK45 control
   (job 7760376) but on the **pure-upwind ADER** path.
3. **TPV31:** the decisive target — pure-upwind + ADER + matrix with the knobs on;
   predict constant $\sigma_n$, dip slip $\to0$, completion to 15 s.

---

# 8. Risks, uncertainties, and fallback

- **Sufficiency is not guaranteed.** Over-integration + resample remove the
  *source* of the speckle (the per-step flux aliasing and the secular
  out-of-space state accumulation, Section 3.3). If a *second* path also feeds
  $[\![ v_n]\!]$ — e.g. the bulk upwind dissipation on the fault-adjacent faces
  itself, the mechanism Zhang (2023) cures with mixed flux — then dealiasing
  alone may only *reduce*, not eliminate, the drift. The Phase 0 A/B is designed
  to detect this.
- **Fallback (still pure-upwind-compatible-ish):** a **Rusanov / Lax–Friedrichs
  near-fault flux** on the fault-adjacent faces of the matrix operator — central
  flux plus a *scalar* dissipation $\tfrac12\lambda_{\max}[\![ Q]\!]$. This is
  what SeisSol added as `fluxNearFault` (Oct 2024, `CellLocalMatrices.cpp:258`);
  unlike pure central it stays dissipative, so it is **ADER-stable** (no RK
  needed) and implementable on the matrix path. It changes the near-fault flux
  (a mild departure from strict pure upwind) but avoids the imaginary-axis
  instability of pure central + ADER.
- **Pure central + ADER is excluded:** non-dissipative $\Rightarrow$ imaginary-axis
  eigenvalues $\Rightarrow$ ADER amplifies them (the documented mixed-flux+ADER
  runaway). Pure central would require porting RK to the matrix path (two pieces
  of new C++); de-prioritised in favour of the dealiasing route.

---

# 9. References

- **Zhang, Liu & Chen (2023)**, *A Mixed-Flux-Based Nodal DG Method for 3-D
  Dynamic Rupture Modeling*, JGR Solid Earth 128, e2022JB025817. (SSOs from
  upwind on asymmetric fault-adjacent meshes; nodal-DG mixed-flux cure.)
- **Pelties, Gabriel & Ampuero (2014)**, *Verification of an ADER-DG method for
  complex dynamic rupture problems*, GMD 7, 847–866. ($(N+2)^2$ GP fault
  over-integration; selective upwind dissipation.)
- **Uphoff (2020)**, PhD thesis, TU München (fault Godunov state Eq. 4.51,
  imposed-state time integral Eq. 4.60 — referenced directly in SeisSol source).
- **Wollherr, Gabriel & Uphoff (2018)**, *Off-fault plasticity ... modal DG*,
  GJI 214, 1556–1584. ($L^2$ projection of fault state back to the polynomial
  space.)
- SeisSol source (read 2026-06-02): `FrictionSolverCommon.h:182,329`;
  `CpuImpl/LinearSlipWeakening.{h:227,261,cpp:18}` (LSW resamples the slip-rate
  magnitude; `BiMaterialFault` disables resample); `CpuImpl/RateAndState.h:89`
  and `CpuImpl/FastVelocityWeakeningLaw.h:144` (rate-state resamples the
  state-variable increment $\Delta\psi$, not the slip rate);
  `Solver/FreeSurfaceIntegrator.h:42` (`PolyDegree = ConvergenceOrder - 1`, i.e.
  polynomial degree $N=$ ConvergenceOrder$-1$ — so "degree $\le$
  ConvergenceOrder$-1$" $=$ degree $\le N$, the full DOF space);
  `codegen/matrices/dr_stroud_matrices_{2..7}.json` (resample matrices, rank
  $(N+1)(N+2)/2$ verified); `Kernels/LinearCK/Neighbor.cpp:71`;
  `Initializer/CellLocalMatrices.cpp:258`.
- SEAS-MFEM: `dynamic/fault_face_flux.cpp:62`; `dynamic/friction_substep_iterator.hpp:201`;
  `drivers/spatial_dyn_driver.cpp:217,811`; isolation job 7762948 (pure-upwind+RK45).

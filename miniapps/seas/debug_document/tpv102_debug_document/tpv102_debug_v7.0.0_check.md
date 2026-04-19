# Code Review v7.0.0: Fresh adversarial review of the v6 R-701 fix

**Versioning convention (introduced this round):** `v[a].[b].[c]`, where
`[a]` bumps on critical bugs, `[b]` on moderate bugs, `[c]` on low bugs.
Previous iteration was `v6` (implicit `v6.0.0`).  This review found a
critical bug, so the version bumps to `v7.0.0`.

**Baseline:** the v6 fix landed the R-701 refactor: `WaveOperator` now
owns a BP5 `FaultBasis` instance, `ComputeSharedFaceFluxRHS` uses the
canonical (ref-normal-aligned) frame reconstructed from BP5's
`FaultBasisQPData`, the R-501 owner-broadcast MPI exchange is deleted
(~270 lines), and `shared_face_peer_`/`SharedFacePeer` are removed.
Tests in `seas_test_r101_shared_fault` (2 ranks) all pass with
`max_rel_diff=0`.

**This review finds the test suite is insufficient to catch a latent
critical bug: the t1 axis has two incompatible conventions in the code
path now — `BuildFrame` (interior fault) gives t1 = strike; `FaultBasis`
(shared fault via R-701) gives t1 = dip.**  The two conventions
coexist in `ComputeFaceFluxRHS` (interior) and `ComputeSharedFaceFluxRHS`
(shared) of the same binary.  DOFData fields `V1, tau1_corr, slip1,
tau1_0` therefore have opposite physical meaning between interior-fault
QPs and shared-fault QPs on the same mesh, which breaks the TPV102
initialisation and output.

## Review Scope
- Plan/history:
  - `tpv102_debug_v6_check.md`, `tpv102_debug_v6_fix.md`.
  - Earlier rounds (`v1`..`v5`) for context on R-001/R-101/R-501 history.
- Files re-reviewed (full re-read, no relying on prior passes):
  - `miniapps/seas/dynamic/wave_operator.hpp`
  - `miniapps/seas/dynamic/wave_operator.inl`
  - `miniapps/seas/dynamic/godunov_flux.cpp` (`BuildFrame`)
  - `miniapps/seas/dynamic/fault_face_flux.cpp` (`Evaluate` / `ComputeTrialTraction`)
  - `miniapps/seas/dynamic/tpv102_setup.hpp` (`InitializeFaultDOFs`)
  - `miniapps/seas/fault/fault_basis.hpp` (public API, not modified)
  - `miniapps/seas/tests/parallel/test_r101_shared_fault.cpp`
- Domain context:
  - `miniapps/seas/CLAUDE.md` invariant: *"BuildFrame aligns t1 with the
    x-axis (along-strike), so tau1_0 = tau_ini is the along-strike
    pre-stress"* (applies to TPV102).
- Runtime verification:
  - `make seas_test_r101_shared_fault` — rebuilds clean.
  - `mpirun -np 2 ./seas_test_r101_shared_fault` — 10/10 PASS.
  - Cross-rank DOFData consistency holds (`max_rel_diff=0`).  **But**
    cross-rank consistency is necessary, not sufficient: both ranks
    can agree on the wrong physical direction and still pass.

## Findings

### [R-801] [CRITICAL] `dynamic/wave_operator.inl:ComputeSharedFaceFluxRHS` — `can_t1 / can_t2` from BP5 `FaultBasis` are `(dip, strike)`, but the rest of the TPV102 pipeline (driver init, interior-fault `ComputeFaceFluxRHS`, station output) assumes `(t1, t2) = (strike, +z)`; shared-fault DOFData is therefore driven as dip-slip while interior-fault DOFData is strike-slip on the same fault

**Category:** BUG — semantic inconsistency of the fault-local tangent axes across the two code paths in `WaveOperator`.

**Description:**
TPV102's initialisation and documented convention (`dynamic/tpv102_setup.hpp:63`,
`CLAUDE.md` "Sign Conventions") assumes that in the fault-local frame:
- `t1` points **along-strike** (= `+x` for TPV102's vertical y=0 fault).
- `V1, tau1_corr, slip1, tau1_0` are the strike-aligned scalar
  components.

This convention is enforced on interior fault faces by
`GodunovFlux::BuildFrame(nor, t1, t2)` at `wave_operator.inl:701`.  With
a typical CalcOrtho result `nor = (0, -1, 0)` on a TPV102 fault face,
`BuildFrame` gives:

| axis | BuildFrame result | physical meaning |
|---|---|---|
| `t1` | `up × nor = (0,0,1) × (0,-1,0) = (+1, 0, 0)` | +x = **strike** |
| `t2` | `nor × t1 = (0,-1,0) × (+1,0,0) = (0, 0, +1)` | +z (up)   |

The R-701 refactor for **shared** fault faces replaces this with BP5's
`FaultBasis::ComputeOrientedFrame` convention (`fault/fault_basis.hpp:54-56`):
*"tangent1 = dip, tangent2 = strike (Tandem convention)"*.  On the same
fault (`ref_normal = (0, -1, 0)`, `up = (0, 0, 1)`), canonical
`(pre-Step-5)` frame produces:

| axis | FaultBasis canonical result | physical meaning |
|---|---|---|
| `can_n`  | `ref-aligned (0, -1, 0)` | −y (fault normal)      |
| `can_t1` | `strike × n = (+1,0,0) × (0,-1,0) = (0, 0, -1)` | −z = **dip-down** |
| `can_t2` | `up × n = (0,0,1) × (0,-1,0) = (+1, 0, 0)` | +x = **strike** |

So on **shared** fault faces, `can_t1 = (0, 0, -1)` — **dip**, not
strike.  Since the code at `wave_operator.inl:1040-1050` rotates
`Q_self` into the canonical frame with this `(can_n, can_t1, can_t2)`
and passes the result to `FaultFaceFlux::Evaluate`, the quantity that
gets stored back as `DOFData.V1, tau1_corr, slip1` on these QPs is the
component along `can_t1` — i.e., the **dip** component, not strike.

The driver's `InitializeFaultDOFs` (`dynamic/tpv102_setup.hpp:67-86`)
writes the same initial condition to **every** fault DOF:

```cpp
d.tau1_0 = TPV102Params::tau_ini;   // 75e6 Pa  (driver says "strike")
d.V1     = TPV102Params::V_ini;     // 1e-12    (driver says "strike-slip")
d.tau2_0 = 0.0;
d.V2     = 0.0;
```

On interior fault DOFs, `(tau1_0, V1)` are **strike** components (75 MPa
along-strike, 1e-12 m/s strike-slip) — consistent with TPV102 spec.
On shared fault DOFs, the same values are interpreted by Evaluate as
**dip** components: 75 MPa of dip-down pre-stress, 1e-12 m/s dip-slip.
There is no strike pre-stress on shared fault DOFs at t=0.

As the rupture propagates from the nucleation zone (interior fault,
driven correctly as strike-slip) into any shared fault QP, the shared
QP tries to rupture in the wrong direction, the Riemann solver
produces nonsense (Brent's method may bracket the wrong side, or V_abs
comes out with the wrong magnitude because `tau1_total = tau1_0 +
tau1_trial` mixes strike-stress-in-dip-slot with dip-trial-from-flux),
and the simulation diverges quietly — but not in a way the existing
R-101 verifier catches, because **both ranks agree on the wrong
value**.

This is exactly the kind of failure mode the user's directive warned
about: *"later we need to wire quasi-dynamic and dynamic together"* —
they will not agree on fault kinematics if they use inconsistent
tangent conventions.

**Trigger:**
Any TPV102 run with at least one shared fault face.  (The R-501a
inline 2-tet test at `mpirun -np 2` already has this condition but
the test only checks cross-rank equality, not physical correctness.)

**Actual behavior:**
- Interior fault DOFs: `V1` = strike-slip rate, `tau1_corr` = corrected
  strike traction.  Consistent with driver init and TPV102 spec.
- Shared fault DOFs: `V1` = dip-slip rate, `tau1_corr` = corrected
  dip traction.  INCONSISTENT with driver init (which wrote strike
  magnitude there) and TPV102 spec (pure strike-slip).
- R-101 verifier: passes (both ranks agree).
- Station output: `V1` column is strike on some QPs, dip on others,
  depending on whether the station's nearest DOF is interior or shared.

**Expected behavior:**
Exactly one convention for `(t1, t2)` throughout every fault DOF and
every code path that writes/reads `DOFData.V1/V2/tau1_corr/tau2_corr/slip1/slip2`.

**Suggested fix — pick option A or option B:**

**Option A** (preferred — uses BP5 FaultBasis everywhere, aligns with
the "reduce code duplication" directive): switch the interior-fault
branch of `ComputeFaceFluxRHS` to also use `FaultBasis` (via the same
`(can_n, can_t1, can_t2)` reconstruction used in shared-fault), and
update `InitializeFaultDOFs` / `ApplyNucleation` to match BP5's
convention:

```diff
 // dynamic/tpv102_setup.hpp:InitializeFaultDOFs
-      d.tau1_0 = TPV102Params::tau_ini;
-      d.tau2_0 = 0.0;
-      // ...
-      d.V1 = TPV102Params::V_ini;  // pure mode-II initial slip
-      d.V2 = 0.0;
+      // BP5 convention: tangent1 = dip, tangent2 = strike.  TPV102 is
+      // pure strike-slip, so the along-strike pre-stress lives in
+      // tau2_0 and the along-strike initial slip rate in V2.
+      d.tau1_0 = 0.0;
+      d.tau2_0 = TPV102Params::tau_ini;    // strike pre-stress
+      // ...
+      d.V1 = 0.0;
+      d.V2 = TPV102Params::V_ini;          // strike initial slip rate

 // dynamic/tpv102_setup.hpp:ApplyNucleation
-      dof_data[i].tau1_0 = TPV102Params::tau_ini + dtau;
+      dof_data[i].tau2_0 = TPV102Params::tau_ini + dtau;

 // dynamic/wave_operator.inl:ComputeFaceFluxRHS fault branch — replace
 // BuildFrame with the FaultBasis canonical frame (same pattern as the
 // shared-fault branch):
-      real_t t1[3], t2[3];
-      GodunovFlux::BuildFrame(nor, t1, t2);
-      DenseMatrix T(NUM_STATE), Tinv(NUM_STATE);
-      GodunovFlux::BuildRotation(nor, t1, t2, T);
-      GodunovFlux::BuildRotationInverse(nor, t1, t2, Tinv);
+      // Look up this interior fault face in FaultBasis by face index.
+      int fb_idx = LookupInteriorFaultBasisIndex(f);   // new helper
+      MFEM_ASSERT(fb_idx >= 0, "interior fault face missing from FaultBasis");
+      const FaultBasisData &bd = fault_basis_->GetBasis(fb_idx);
+      MFEM_ASSERT(q < (int)bd.qp_data.size(),
+                  "FaultBasis::qp_data not populated for interior fault face");
+      const FaultBasisQPData &qpd = bd.qp_data[q];
+      real_t can_n[3], can_t1[3], can_t2[3];
+      for (int d = 0; d < 3; d++) {
+         can_n[d]  = qpd.sign_flipped ? -qpd.normal[d]   : qpd.normal[d];
+         can_t1[d] = qpd.sign_flipped ? -qpd.tangent1[d] : qpd.tangent1[d];
+         can_t2[d] = qpd.sign_flipped ? -qpd.tangent2[d] : qpd.tangent2[d];
+      }
+      DenseMatrix T(NUM_STATE), Tinv(NUM_STATE);
+      GodunovFlux::BuildRotation(can_n, can_t1, can_t2, T);
+      GodunovFlux::BuildRotationInverse(can_n, can_t1, can_t2, Tinv);
```

Plus update any consumer of `DOFData.V1/tau1_corr/slip1` — station writer,
ParaView output — to know that BP5's convention is in force (dip in
component 0, strike in component 1).  Searching the driver:

```
drivers/tpv102_driver.cpp:635-642: pv_local_slip_rate(2*i + 0) = d.V2;  // BP5 comp0 = dip
drivers/tpv102_driver.cpp:636:      pv_local_slip_rate(2*i + 1) = d.V1;  // BP5 comp1 = strike
```

The driver *already* swaps the components for ParaView output ("V1 →
strike channel, V2 → dip channel") because BP5's ParaView output uses
BP5 convention.  But the SOURCE of V1 (how it's populated during
Mult) is currently in two different conventions.  Option A makes the
source BP5-canonical everywhere and the driver's swap for ParaView
output is already correct.

**Option B** (preferred if Option A's driver rework is too invasive):
remap the BP5 FaultBasis output so shared-fault code sees
`(can_t1, can_t2) = (strike, dip-up)` matching BuildFrame:

```diff
 // dynamic/wave_operator.inl:ComputeSharedFaceFluxRHS (R-701 block):
-for (int d = 0; d < 3; d++) {
-   can_n[d]  = qpd.sign_flipped ? -qpd.normal[d]   : qpd.normal[d];
-   can_t1[d] = qpd.sign_flipped ? -qpd.tangent1[d] : qpd.tangent1[d];
-   can_t2[d] = qpd.sign_flipped ? -qpd.tangent2[d] : qpd.tangent2[d];
-}
+// R-801 fix: remap BP5's (dip, strike) → TPV102's (strike, +z) so
+// shared-fault DOFData has the same t1/t2 semantic as interior-fault
+// DOFData (which uses GodunovFlux::BuildFrame).
+// BP5 canonical: tangent1 = dip, tangent2 = strike.
+// TPV102 BuildFrame for vertical y=0 fault: t1 = strike, t2 = dip-up.
+// Map: tpv_t1 = bp5_t2 (strike), tpv_t2 = -bp5_t1 (dip-up = -dip-down).
+for (int d = 0; d < 3; d++) {
+   can_n[d]  = qpd.sign_flipped ? -qpd.normal[d]   : qpd.normal[d];
+   can_t1[d] = qpd.sign_flipped ? -qpd.tangent2[d] : qpd.tangent2[d];   // strike
+   can_t2[d] = qpd.sign_flipped ?  qpd.tangent1[d] : -qpd.tangent1[d];  // dip-up
+}
```

This option requires no driver changes, but adds a non-obvious remap
every time FaultBasis is consumed by TPV102 code.  Document it
prominently.

**Test case:**
```cpp
// tests/parallel/test_r101_shared_fault.cpp — new T-R801:
// Drive multi-stage RK4 on a nonzero Q with the TPV102 physical
// initial condition (tau_ini strike pre-stress), then verify that
// for a SHARED fault QP the resulting V1/tau1_corr/slip1 are in the
// STRIKE direction, not dip.
//
// Discriminator: for a pure-strike-slip TPV102 setup, |V1| should
// dominate |V2| on interior fault DOFs.  If the code correctly
// carries the strike convention through the shared-fault path,
// shared fault DOFs must also have |V1| >> |V2|.  Under the current
// R-701 code, shared fault DOFs come out with |V2| >> |V1| (swapped
// semantics).
static void TestR801_StrikeSlipConventionOnSharedFault(int rank,
                                                        MPI_Comm comm,
                                                        int nprocs)
{
   if (nprocs != 2) { /* SKIPPED */ return; }
   Mesh serial = BuildTwoTetSharedFaultMeshInline();
   int part[2] = {0, 1};
   ParMesh pmesh(comm, serial, part);
   BoundaryConfig bc; bc.natural_attrs={1}; bc.fault_attr=3; bc.absorbing_attrs={5};
   WaveOperator<ParMesh> wave(pmesh, 1, TPV102Params::lambda,
                              TPV102Params::mu, TPV102Params::rho, bc);

   // Populate fault DOFData with the standard TPV102 setup (pure
   // strike-slip, tau_ini along-strike).
   const int n_shared = wave.GetFaultSharedFaces().Size();
   int nqp = 0;
   if (n_shared > 0) {
      auto *ftr = pmesh.GetSharedFaceTransformations(wave.GetFaultSharedFaces()[0]);
      nqp = IntRules.Get(ftr->GetGeometryType(), 2).GetNPoints();
   }
   int max_nqp = nqp;
   MPI_Allreduce(&nqp, &max_nqp, 1, MPI_INT, MPI_MAX, comm);
   nqp = max_nqp;

   const int num_fault_total = n_shared * nqp;   // inline: no interior
   std::vector<Vector> fault_coords;
   BuildFaultCoords(pmesh, wave.GetFaultInteriorFaces(),
                    wave.GetFaultSharedFaces(), 1, nqp, fault_coords);
   std::vector<DOFData> dof_data;
   InitializeFaultDOFs(dof_data, num_fault_total, fault_coords);

   FaultFaceFlux ff(TPV102Params::rho, TPV102Params::cp, TPV102Params::cs);
   wave.SetFaultFlux(&ff);
   wave.SetFaultDOFData(&dof_data, nqp);

   // Drive Mult on a LARGE nonzero Q so flux-driven friction response is
   // non-trivial.  Scale chosen so V_abs from Brent solve is O(1 m/s).
   Vector Q(wave.Height()); Q = 0.0;
   for (int i = 0; i < Q.Size(); i++) { Q(i) = 1e6 * std::sin(0.1 * i); }
   Vector k(Q.Size());
   for (int s = 0; s < 4; s++) { wave.Mult(Q, k); }

   // Assertion: on a pure-strike-slip TPV102 setup, after Mult the
   // slip rate magnitude |V1| (supposed to be strike) must dominate
   // |V2| (supposed to be dip).  Current R-701 code fails this because
   // shared-face path has them swapped.
   for (int i = 0; i < num_fault_total; i++) {
      const DOFData &d = dof_data[i];
      real_t V1_mag = std::abs(d.V1);
      real_t V2_mag = std::abs(d.V2);
      // If V2 dominates V1 by >10x on a pure-strike-slip setup,
      // the convention is wrong.
      TEST_ASSERT(V1_mag > V2_mag,
                  "R-801: pure strike-slip expects |V1| > |V2| on "
                  "shared fault DOF " + std::to_string(i) + " (got |V1|="
                  + std::to_string(V1_mag) + ", |V2|=" + std::to_string(V2_mag) + ")");
   }
}
```

---

### [R-802] [CRITICAL] [POSSIBLE] `dynamic/wave_operator.inl:ComputeSharedFaceFluxRHS` — flux accumulation uses rank-local `nor` (from `CalcOrtho`) under the implicit assumption `nor_A = -nor_B`; the v6 fix report stated empirically `nor_A = nor_B` on the inline mesh; if that observation holds, the Godunov conservation identity `F(L, R, +n) = -F(R, L, -n)` does NOT compensate and both ranks add the SAME `F_h` to their Elem1 rhs — breaking momentum conservation across the shared fault

**Category:** BUG (POSSIBLE) — correctness of the flux accumulation

**Description:**
The v6 fix report `tpv102_debug_v6_fix.md:98-104` documents an
empirical observation:
> *"My added diagnostics showed both ranks receive the SAME
> `nor = (0, 1, 0)` and the SAME `sign_flipped = 1`, even though their
> local `Elem1` is on opposite physical sides."*

This motivated the geometric `elem1_on_plus_` swap flag.  Good — for
the Evaluate *inputs*, the swap + canonical frame reconstruction means
both ranks feed Evaluate with bit-identical arguments.

But the subsequent flux assembly at `wave_operator.inl:1093`:

```cpp
flux_.Interior(nor, Q_self_imp, Q_nbr_imp, F_h);
```

uses the rank-local `nor` — i.e., the raw CalcOrtho result on this
rank.  `Q_self_imp`, `Q_nbr_imp` are in the GLOBAL frame (rotated back
from canonical via `T_can`), bit-identical on both ranks after the
`elem1_on_plus` self/nbr swap.

So: on rank A, `F_h_A = flux_.Interior(nor_A, Q_imp_plus_g, Q_imp_minus_g)`.
On rank B with `elem1_on_plus = false` (opposite), the swap means
`(Q_self_imp, Q_nbr_imp) = (Q_imp_minus_g, Q_imp_plus_g)`, so
`F_h_B = flux_.Interior(nor_B, Q_imp_minus_g, Q_imp_plus_g)`.

The Godunov conservation identity is:
```
F(Q_L, Q_R, +n) = -F(Q_R, Q_L, -n)
```

If `nor_A = -nor_B` (the classical MFEM convention):
```
F_h_B = F(Q_imp_minus_g, Q_imp_plus_g, -n_A) = -F(Q_imp_plus_g, Q_imp_minus_g, +n_A) = -F_h_A.
```
Both ranks accumulate `rhs[Elem1] -= F_h * w`, which gives
`rhs[Elem1_A] -= F_h_A * w` and `rhs[Elem1_B] += F_h_A * w`.  Conservation ✓.

**But if the v6 empirical observation holds** and `nor_A = nor_B`:
```
F_h_B = F(Q_imp_minus_g, Q_imp_plus_g, +n_A) ≠ -F_h_A.
```

Using the splitting identity `F(L, R, +n) = A_n^+ L + A_n^- R`:
```
F_h_A = A_{n_A}^+ Q_imp_plus_g + A_{n_A}^- Q_imp_minus_g.
F_h_B = A_{n_A}^+ Q_imp_minus_g + A_{n_A}^- Q_imp_plus_g.
```
These two values are independent, not opposite.  Both ranks do
`rhs[Elem1] -= F_h * w`, giving TWO independent subtractions in the
rhs.  **Momentum is not conserved across the shared fault.**

Concretely, sum `rhs[Elem1_A] + rhs[Elem1_B] = -(F_h_A + F_h_B)·w ≠ 0`
when it should be (= 0 for conservation).

**The R-501a test does not catch this.**  The test only verifies
cross-rank DOFData equality (which depends on `Evaluate` inputs, not
on `F_h` accumulation).  It does not integrate multi-step dynamics
long enough to expose the conservation violation; after one Mult the
RHS bears the wrong flux, but since `Q` is not advanced between Mult
calls in the test, no drift compounds into observable output.

**Uncertainty:**  The v6 fix report says `nor_A = nor_B` on the inline
mesh but does not provide diagnostic evidence for the general TPV102
mesh.  MFEM's documented convention (`pmesh.hpp:592-597`) says Elem1
is the element whose outward normal matches the face normal — which
*should* produce `nor_A = -nor_B` on shared faces.  If the fix report's
empirical observation is specific to the inline-mesh topology (e.g.,
due to how `ParMesh` constructor partitions a hand-crafted mesh with
`int partition[2] = {0,1}`), then on real TPV102 meshes the classical
convention may hold and R-802 is a non-issue.  **This is why R-802 is
flagged POSSIBLE: it is unverified either way.**

**Trigger:**
Any TPV102 run with a shared fault face, if MFEM's shared-face
CalcOrtho actually returns identical (not opposite) normals on the
two ranks.

**Actual behavior:**
Both ranks' Elem1 rhs subtracts an independent `F_h` that are not
mutually opposite.  Sum is not zero.  Bulk momentum is not
conserved across the fault seam.

**Expected behavior:**
The flux accumulation must be conservation-compliant regardless of
MFEM's nor convention on the local rank.

**Suggested fix:**
Use the **canonical normal** `can_n` (which is identical on both ranks
by construction) for the flux computation, not the rank-local `nor`:

```diff
-flux_.Interior(nor, Q_self_imp, Q_nbr_imp, F_h);
+// R-802 fix: use canonical nor so both ranks compute identical F_h.
+// Rank A: F_h_A = F(Q_imp_plus_g, Q_imp_minus_g, can_n).
+// Rank B: F_h_B = F(Q_imp_minus_g, Q_imp_plus_g, can_n) due to
+//          elem1_on_plus swap.
+// The accumulation sign is then determined by elem1_on_plus:
+//   Elem1 on canonical-+ side: rhs[Elem1] -= F_h (flux leaves Elem1).
+//   Elem1 on canonical-- side: rhs[Elem1] += F_h (flux enters Elem1).
+flux_.Interior(can_n, Q_imp_plus_g, Q_imp_minus_g, F_h);
+const real_t accum_sign = elem1_on_plus ? +1.0 : -1.0;
+for (int c = 0; c < NUM_STATE; c++)
+   for (int i = 0; i < ndof; i++)
+      rhs[c*ndof_total_ + dof_offset1 + i] -= accum_sign * w * shape1(i) * F_h[c];
```

Conservation holds under this formulation regardless of whether MFEM
gives `nor_A = nor_B` or `nor_A = -nor_B` on shared faces.

**Test case:**
```cpp
// tests/parallel/test_r101_shared_fault.cpp — new T-R802:
// Conservation across the shared fault face: sum of per-rank
// rhs[Elem1]·w over the shared-face QPs must be ≈ 0.
//
// Preconditions: 2-tet inline, 2 ranks, Q = nonzero.
// After one Mult, gather rhs on both ranks and sum the element-wise
// contributions into a scalar.  On pre-R-802 code with nor_A = nor_B,
// the sum is O(|F_h|); with the fix it is zero to round-off.

static void TestR802_ConservationAcrossSharedFault(int rank, MPI_Comm comm,
                                                    int nprocs)
{
   // ... (as above, build the context, drive one Mult on nonzero Q)
   Vector k(Q.Size()); wave.Mult(Q, k);

   // Extract per-element sum of k contributions from rank-A's Elem1
   // and rank-B's Elem1 (both adjacent to the shared fault face).
   // For a tet element at the fault, the shared-face flux
   // contributes to the 3 face DOFs of that element.
   real_t local_sum = 0.0;
   // ... sum over the face DOFs of Elem1_local ...

   real_t global_sum = 0.0;
   MPI_Allreduce(&local_sum, &global_sum, 1, MPITypeMap<real_t>::mpi_type,
                 MPI_SUM, comm);

   // Tolerance: O(machine epsilon × |F_h|).  Pre-R-802 with the empirical
   // nor_A = nor_B, this sum is O(MPa × face-weight), orders of
   // magnitude above tolerance.
   TEST_ASSERT(std::abs(global_sum) < 1e-8,
               "R-802: flux-sum across shared fault must conserve "
               "(got " + std::to_string(global_sum) + ")");
}
```

---

## Summary
- Critical issues: **2** (R-801 — semantic t1 mismatch between interior
  and shared fault; R-802 POSSIBLE — conservation violation if MFEM
  gives `nor_A = nor_B` on shared faces)
- Moderate issues: **0**
- Low issues: **0**
- Plan compliance (against v6 directive "reuse BP5 FaultBasis, don't
  modify BP5"): **FULL on the reuse question; INCOMPLETE on the
  semantic-consistency question** — the reuse introduced a convention
  mismatch that the v6 tests do not cover.
- **Verdict: FAIL — must fix before proceeding to Frontera.**
  The R-701 refactor achieves cross-rank DOFData bit-identity, which
  was the v5 goal.  But it introduces a convention split between
  interior-fault and shared-fault code paths that corrupts the physical
  meaning of TPV102's strike-slip initial condition.  Before the next
  Frontera run:
    1. Fix R-801 (option A or option B).
    2. Decide R-802 empirically — add a diagnostic or test that
       verifies conservation across the shared fault.  Apply the
       R-802 fix if the conservation test fails.

## Version bump rationale
- Pre-existing: `v6` (implicit `v6.0.0`).
- This review: **2 CRITICAL** findings (R-801, R-802).  Bump `[a]` to 7.
- `[b]` and `[c]` reset to 0.
- **New version: `v7.0.0`**.  File: `tpv102_debug_v7.0.0_check.md`.

## Unreviewed Areas
- Whether the observed `nor_A = nor_B` on the inline 2-tet mesh (v6
  fix report) generalises to real TPV102 production meshes (1000 m,
  200 m) or is an inline-mesh artifact.  Requires a Frontera-side
  diagnostic run or a more comprehensive unit test.
- Interior-fault behaviour when CalcOrtho produces inconsistent normal
  orientations across different interior fault faces on the same rank.
  TPV102's Gmsh-generated mesh likely gives consistent orientation,
  but no test asserts this.  Not part of the R-701 regression; flagged
  for future work.
- BP5 non-regression (`seas_test_fault_basis`, `seas_test_elasticity_operator`)
  was reported by v6_fix as having pre-existing compile errors
  (stale 11-arg callers of `ProjectTractionToFaultDOFs`).  This review
  did not rebuild/test the BP5 suite to confirm the pre-existing-ness;
  those compile errors predate R-701 per the fix report.  Flag for
  separate cleanup.
- T-G1 (4-rank 1000 m driver run) from the v6 test plan was not
  executed in this review round — running it locally requires the
  conda mfem-dev env and a 1000 m mesh file.  Recommend running after
  R-801 fix to confirm the Frontera-shape test scenario passes.

# Conformal-class Higgs selector theorem

**Lane ID:** `CROSS-RFT-ALPHA-GRA-FA-CONFORMAL-CLASS-HIGGS-SELECTOR-V001`

**Official short name:** `CCHS`

**Date:** 2026-08-27

**Builder status:** `MUTABLE_PRESCREEN_READY__BUILDER_REPLAY_PASS__NOT_SEALED`

**Claim class:** exact conditional conformal-Laplacian uniqueness theorem;
exact characteristic-only no-go; exact conditional composition from AQ4DL
dimension four through the Visser-table conformal Higgs coupling to EY's
visible-sector induced Ricci coefficient

**Not claimed:** that causal order alone supplies full Weyl covariance; that
the Higgs equation physically exists before a volume representative is
selected; quantum Weyl invariance; an RG-fixed Higgs curvature coupling; a
complete ultraviolet spectrum; a total positive Newton coefficient; physical
record-to-metric soldering; gravity; or outcome selection

## 1. Exact question and answer

Suppose an E-EMERGENTSPACE construction reaches an intermediate continuum
stage at which complete causal order has fixed a smooth Lorentzian conformal
class `[g]`, but no physical volume observable has yet selected one metric
representative.  If a massless Higgs scalar equation must already be a
well-defined equation on `[g]`, does that requirement select the nonminimal
curvature coupling used by EY?

There are two different answers to two often-confused meanings of
"well-defined":

1. If only the characteristic cone must be independent of the representative,
   **no value of the curvature coupling is selected**.  Every operator in the
   family below has the same principal symbol and null characteristics.
2. If the complete weighted scalar equation must intertwine under every smooth
   positive Weyl rescaling, and no Weyl gauge field, dilaton, compensator, or
   higher-derivative replacement is admitted, then the curvature coupling is
   uniquely

   \[
    \xi_D={D-2\over4(D-1)}.
   \]

AQ4DL's conditional `D=4` output then gives the Visser-table value
`xi_H=1/6`.  In EY's declared visible-only proper-time census this yields

\[
 \operatorname{str}k_1^{\rm vis}=-{1\over2},\qquad
 C_{R,\rm vis}^{>}={\kappa_R^2-\mu^2\over64\pi^2}>0.
\]

The full-Weyl-intertwining premise is not supplied by the current
F3 + E-EMERGENTSPACE adoption, AQ4DL, or causal-cone recovery alone.  The
composition is therefore a sharp conditional selector and a concrete next
same-parent target, not an actual-world closure.

## 2. Frozen operator class

Let `D>2`, let `g` be one representative of a smooth Lorentzian conformal
class, and let

\[
 \widehat g=e^{2\sigma}g                                  \tag{FA01}
\]

for arbitrary smooth real `sigma`.  Consider one component of the unbroken
massless Higgs doublet, or the gauge-covariant kinetic Hessian acting on that
component, in the local natural second-order class

\[
 P_{g,\xi}\phi=-\nabla_g^2\phi+\xi R_g\phi.               \tag{FA02}
\]

Equation (FA02) uses the Visser table convention: its scalar heat-kernel
coefficient is `k_1=1/6-xi` in `D=4`, and the table-conformal value is
`xi=1/6`.  If a Lorentzian sign convention writes the wave equation with the
overall opposite sign, multiplying the complete equation by `-1` changes no
solution and does not change the table parameter.

For the complex Higgs doublet use
`H^dagger H=(1/2) sum_i phi_i^2`; hence
`xi R H^dagger H=(1/2)xi R sum_i phi_i^2`.  Each of the four real components
therefore carries the same table parameter `xi`, with no factor-of-two change.

The following restrictions are load bearing:

1. the principal part is the canonical two-derivative scalar Laplacian;
2. the only local dimension-two geometric zeroth-order term is `xi R`;
3. there is no mass term at this stage;
4. no independent Weyl connection, dilaton, compensator, density scale, or
   selected volume representative is supplied;
5. the internal Standard-Model connection, when retained, is a Weyl-weight
   zero one-form and does not supply an extra scale; and
6. "defined on `[g]`" means the full equation on a weighted scalar density,
   not merely equality of its characteristics.

Higher-derivative conformal operators, nonlocal kernels, and compensated Weyl
theories are lawful alternatives, but they are outside this theorem.

## 3. Full conformal-class premise

Freeze the prospective physical premise:

> **`PREVOLUME-WEYL-HIGGS`.**  Before any physical volume observable, clock
> normalization, or other scale datum selects a representative of `[g]`, the
> same parent already contains the massless Higgs kinetic equation (FA02) as a
> natural equation on the conformal class.  For every representative change
> `g -> exp(2 sigma) g`, a fixed scalar weight `w` maps solutions to solutions
> through one local intertwining law, with no compensator or additional
> background field.

"Before" here is logical construction/dependency order in the non-temporal
E-EMERGENTSPACE scale flow, not a dated cosmological epoch.

Operationally, this requires functions `w` and `w'`, independent of `sigma`
and of the solution, such that

\[
 P_{e^{2\sigma}g,\xi}\!left(e^{w\sigma}\phi\right)
 =e^{w'\sigma}P_{g,\xi}\phi                             \tag{FA03}
\]

for every smooth `sigma` and `phi`.  Equality only after choosing one special
representative, only on constant Weyl rescalings, or only at the level of the
null cone does not satisfy this premise.

## 4. Theorem CCHS-1 -- exact conformal-Laplacian uniqueness

Under the operator class of section 2 and `PREVOLUME-WEYL-HIGGS`,

\[
 \boxed{
 w=-{D-2\over2},\qquad
 w'=-{D+2\over2},\qquad
 \xi=\xi_D:={D-2\over4(D-1)}.}                    \tag{FA04}
\]

The exact intertwining law is

\[
 \boxed{
 P_{e^{2\sigma}g,\xi_D}
 \left(e^{-(D-2)\sigma/2}\phi\right)
 =e^{-(D+2)\sigma/2}P_{g,\xi_D}\phi.}             \tag{FA05}
\]

Thus the equation `P phi=0` depends only on `[g]` when `phi` is typed as the
corresponding weighted scalar density.

### Proof

In the convention of (FA02), direct transformation gives

\[
\begin{aligned}
 R_{e^{2\sigma}g}
 &=e^{-2\sigma}\left[
 R_g-2(D-1)\nabla^2\sigma
 -(D-1)(D-2)(\nabla\sigma)^2\right],               \tag{FA06}\\
 -\widehat\nabla^2(e^{w\sigma}\phi)
 &=e^{(w-2)\sigma}\left[
 -\nabla^2\phi-(2w+D-2)\nabla\sigma\!\cdot\!\nabla\phi
 -w(\nabla^2\sigma)\phi
 -w(w+D-2)(\nabla\sigma)^2\phi\right].            \tag{FA07}
\end{aligned}
\]

For (FA03) to hold for arbitrary `sigma` and `phi`, the independent
`grad(sigma).grad(phi)` coefficient must vanish:

\[
 2w+D-2=0
 \quad\Longrightarrow\quad
 w=-{D-2\over2}.                                   \tag{FA08}
\]

The independent `(nabla^2 sigma) phi` coefficient must also vanish:

\[
 -w-2\xi(D-1)=0
 \quad\Longrightarrow\quad
 \xi={D-2\over4(D-1)}.                             \tag{FA09}
\]

With (FA08)--(FA09), the remaining `(grad sigma)^2 phi` coefficient is

\[
 -w(w+D-2)-\xi(D-1)(D-2)=0,                        \tag{FA10}
\]

and the common prefactor in (FA07) gives `w'=w-2=-(D+2)/2`.
This proves existence and uniqueness inside the frozen operator class. QED.

The same calculation holds with a Weyl-weight-zero internal gauge connection
in place of the ordinary derivative because `sigma` is an internal singlet;
the curvature-coupling coefficient is unchanged.

## 5. Theorem CCHS-2 -- characteristic-only no-go

For every real `xi`, the principal symbol of (FA02) is

\[
 \operatorname{Prin}(P_{g,\xi})(x,k)=g^{\mu\nu}(x)k_\mu k_\nu. \tag{FA11}
\]

Under `g -> exp(2 sigma) g`, this symbol is multiplied by the nonzero factor
`exp(-2 sigma)`, so its zero set is unchanged.  Hence

\[
 \boxed{
 \text{common causal cone or characteristic propagation alone}
 \not\Longrightarrow \xi={D-2\over4(D-1)}.}        \tag{FA12}
\]

### Proof

The `xi R` term is zeroth order and does not enter the principal symbol.
Conformal rescaling multiplies the inverse metric by `exp(-2 sigma)`, which
does not change the characteristic zero set.  Therefore every `xi` has the
same conformal null cone, although only (FA04) satisfies the complete
intertwining law in the frozen class. QED.

This no-go prevents the F3 common-cone obligation, by itself, from being used
as a hidden proof of conformal Higgs coupling.

## 6. Corollary CCHS-3 -- conditional AQ4DL-to-EY composition

Assume one same parent satisfies all of the following:

1. AQ4DL's `QFRONT-DIM`, same-front Maxwell, and `MARGINAL-ALPHA` premises,
   so `D=q=4`;
2. a smooth causal-order stage genuinely fixes `[g]` before a volume
   representative is selected;
3. `PREVOLUME-WEYL-HIGGS` and the frozen operator class hold for the unbroken
   massless Higgs kinetic Hessian;
4. the selected representative and Visser/RIEHB curvature dictionary match
   the one later used by EY; and
5. EY's visible-only (`N_p=0`) zero-temperature proper-time shell,
   anomaly-free field census, thresholds, one-loop, and matching premises all
   hold.

Then CCHS-1 gives

\[
 \boxed{D=4\quad\Longrightarrow\quad\xi_H={1\over6}}             \tag{FA13}
\]

in EY's Visser-table convention.  EY then gives

\[
\begin{aligned}
 \operatorname{str}k_1^{\rm vis}
 &= {1\over6}-4\left({1\over6}\right)=-{1\over2},\\
 \boxed{
 C_{R,\rm vis}^{>}
 ={\kappa_R^2-\mu^2\over64\pi^2}>0.}              \tag{FA14}
\end{aligned}
\]

This closes EY's visible-only `xi_H` sign gate conditionally.  It does not
close the total-spectrum gate.  In particular, if `N_p` additional minimally
coupled real scalars are actually present in the same shell, EY gives

\[
 C_R^{\rm vis+p}
 ={3-N_p\over192\pi^2}(\kappa_R^2-\mu^2),          \tag{FA15}
\]

which is positive only for integer `N_p=0,1,2`, vanishes for `N_p=3`, and is
negative for `N_p>=4`.  Six PMMDC coordinates therefore cannot be silently
counted as six minimally coupled scalar determinants.  If separately proved
to be conformally coupled scalars, their individual `k_1` terms vanish at
this order; neither propagation nor that coupling is proved here.

## 7. Existing-premise audit

The exact conditional composition is not presently an unconditional program
result:

1. **F3 + E-EMERGENTSPACE adoption** requires an incidence-to-causal map, a
   common cone, continuum recovery, and a protected tensor sector.  It does
   not assert a smooth conformal-class-only stage or require Higgs dynamics to
   be fully Weyl-intertwined before scale/volume selection.
2. **AQ4DL** conditionally selects `D=q=4` from a classically marginal Maxwell
   sector.  It explicitly does not derive that sector, causal geometry, or a
   Higgs curvature coupling.
3. **The conditional G3 observation** that complete chronology fixes `[g]`
   requires a smooth causally regular Lorentzian continuum.  Even when those
   premises hold, the order theorem fixes kinematics; it does not require the
   Higgs equation to exist at the pre-volume stage.
4. **EO/EV/ET and the current gamma-soldering lanes** provide finite or
   conditional metric antecedents and refinement gates, not the physical
   `PREVOLUME-WEYL-HIGGS` law.

Therefore `PREVOLUME-WEYL-HIGGS` is a new, falsifiable same-parent
architecture/physics premise.  It must be derived from the record-conditioned
parent or tested as a prospective construction rule; it may not be inferred
from the word "causal."

## 8. Quantum, RG, and phase ceilings

1. **Classical versus quantum covariance.**  CCHS-1 is a classical operator
   theorem.  In four dimensions the renormalized Standard Model has running
   gauge, Yukawa, and scalar couplings and a trace/Weyl anomaly.  The complete
   quantum effective action therefore need not be a functional of `[g]`
   alone.
2. **RG custody of `xi_H`.**  Equation (FA13) fixes the curvature coefficient
   of the massless kinetic Hessian at the pre-volume matching stage.  An
   RG-improved or all-orders use of (FA14) additionally requires proof that
   the renormalized `xi_H` entering the declared shell remains `1/6`, or an
   explicit integration of its running and thresholds.  No beta function is
   assumed here.
3. **Strict one-loop scope.**  At strict one loop, inserting the frozen
   classical Hessian `xi_H=1/6` into EY gives (FA14).  Quantum corrections to
   that input, anomaly terms, and running belong to higher-order/matched
   custody, not to the uniqueness algebra.
4. **Representative selection is necessary.**  EY's physical crossover,
   proper-time window, measure, and induced Einstein--Hilbert term all use a
   selected metric representative and scale.  Exact Weyl covariance is only
   posited at the earlier conformal-class stage; it is not claimed to survive
   representative selection or the physical cutoff.
5. **Electroweak phase.**  The Higgs mass term, vacuum expectation value, and
   resolved broken-phase thresholds violate the massless premise.  The join
   applies only to EY's declared ultraviolet symmetric-gauge-basis shell.
6. **Anomaly does not change the classical selector into a contradiction.**
   It prevents promotion to an exact quantum Weyl symmetry.  Curvature-squared
   and logarithmic anomaly terms remain in RIEHB's remainder ledger and must
   be bounded on the Einstein band.

## 9. Controls and exact disposition

1. **Characteristic control.**  Use common null cones to select `xi`.
   Rejected by CCHS-2.
2. **Representative control.**  Choose a metric volume first and then call an
   equation on that chosen metric an equation intrinsically on `[g]`.
   Rejected.
3. **Compensator control.**  Add a Weyl gauge field or dilaton and retain the
   uniqueness conclusion.  Rejected; the operator class changed.
4. **Mass control.**  Keep a fixed nonzero Higgs mass at the pre-volume stage.
   Rejected; it supplies a scale and breaks (FA03).
5. **Higher-derivative control.**  Replace the canonical second-order
   principal part by another conformal operator.  Rejected; recompute the
   theorem in that class.
6. **Convention control.**  Insert `1/6` into EY without matching the Visser
   table curvature/Hessian convention.  Rejected.
7. **Quantum control.**  Promote classical conformal covariance to an exact
   anomaly-free quantum Standard Model.  Rejected.
8. **RG control.**  Hold `xi_H=1/6` across an arbitrary shell without a
   matching or running argument.  Rejected.
9. **Pair-mode control.**  Add six minimally coupled PMMDC scalar determinants
   while quoting the visible-only positive coefficient.  Rejected by (FA15).
10. **Gravity control.**  Relabel one positive visible one-loop coefficient as
    a complete Newton constant or record-origin gravity theorem.  Rejected.

**Disposition:**

`FULL_WEYL_INTERTWINING_OF_CANONICAL_MASSLESS_SCALAR_EQUATION_ON_CONFORMAL_CLASS_UNIQUELY_SELECTS_XI_D_(D_MINUS_2)_OVER_4_(D_MINUS_1)__CHARACTERISTIC_CONE_ALONE_SELECTS_NO_XI__AQ4DL_D4_PLUS_NEW_PREVOLUME_WEYL_HIGGS_PREMISE_GIVES_VISSER_TABLE_XI_H_ONE_OVER_6__EY_VISIBLE_ONLY_ONE_LOOP_COEFFICIENT_POSITIVE_DELTA_OVER_64_PI2__EXISTING_ADOPTED_PREMISES_DO_NOT_SUPPLY_PREVOLUME_WEYL_HIGGS__QUANTUM_ANOMALY_RG_THRESHOLDS_COMPLETE_SPECTRUM_SOLDERING_AND_GRAVITY_OPEN__MUTABLE_NOT_SEALED`

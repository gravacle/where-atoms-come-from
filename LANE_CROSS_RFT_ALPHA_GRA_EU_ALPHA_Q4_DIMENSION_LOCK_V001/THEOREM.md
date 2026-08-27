# Alpha--q4 dimension-lock selector theorem

**Lane ID:** `CROSS-RFT-ALPHA-GRA-EU-AQ4DL-V001`

**Official short name:** `AQ4DL`

**Date:** 2026-08-27

**Builder status:**
`MUTABLE_PRESCREEN_READY__BUILDER_REPLAY_PASS__NOT_SEALED`

**Claim class:** exact engineering-dimension theorem; exact conditional
same-front selector from a classically marginal visible Maxwell charge to a
four-operation count/contrast spacetime

**Not claimed:** that RFT alone selects four operations; that EO's supplied
q=4 antecedent is physically realized by nature; derivation of a Maxwell
sector from records; derivation or numerical calculation of alpha; exclusion
of dimensionful gauge theories or interacting fixed points in other
dimensions; gravity; Newton's constant; or outcome selection

## 1. Exact question

EO shows that a supplied operation-count front separates into one formation
depth and a contrast quotient.  For `q` complete reusable operations, the
same algebra has one depth coordinate and `q-1` independent contrasts.  This
does not select `q`.

The actual visible electromagnetic sector supplies a possible REQUIRE-side
selector.  The question is:

\[
 \begin{split}
 &\text{If the same emergent front supports one canonically normalized}\
 &\text{Maxwell sector whose inherited charge coupling is classically}\
 &\text{marginal and needs no compensating length/mass scale, must }q=4?
 \end{split}                                                   \tag{EU01}
\]

The answer is yes under the explicit same-front premises below.  The result
is a dimension selector, not a gravity theorem.

## 2. Count/contrast dimension premise

For an integer `q>=2`, let

\[
 P_q=I_q-{1\over q}\mathbf1\mathbf1^{\mathsf T},
 \qquad t=\mathbf1^{\mathsf T}m,
 \qquad \xi=P_qm.                                      \tag{EU02}
\]

Then

\[
 P_q^2=P_q,\qquad P_q\mathbf1=0,
 \qquad \operatorname{rank}P_q=q-1,                  \tag{EU03}
\]

and

\[
 \mathbb R^q=\operatorname{span}\{\mathbf1\}
 \oplus\mathbf1^\perp.                               \tag{EU04}
\]

Freeze the following physical premise:

> **`QFRONT-DIM`.**  In one admitted infrared window, formation depth is one
> nondegenerate physical time coordinate, the full contrast quotient
> `1^perp` supplies all and only the nondegenerate spatial tangent directions,
> and no additional, compactified, gapped, null, gauge, or quotient direction
> is counted as a physical tangent direction in that same window.

Under `QFRONT-DIM`, the emergent spacetime dimension is

\[
 \boxed{D=1+(q-1)=q.}                                \tag{EU05}
\]

Equation (EU05) is the general count/contrast dimension identity underlying
EO's `q=4` construction.  The sealed EO theorem realizes its algebra and
calibrated local metric conditionally at `q=4`; it does not prove
`QFRONT-DIM` for nature or select the value of `q`.

## 3. Same-front Maxwell marginality premise

On that same `D`-dimensional continuum front, use natural units only for
engineering dimensions: `[partial_mu]=1` and a local Lagrangian density has
mass dimension `D`.  Assume one visible compact `U(1)` eigenmode with a
canonically normalized two-derivative Maxwell term and at least one
canonically normalized charged matter field.  For a Dirac representative,

\[
 \mathcal L_{\rm vis}
 =-{1\over4}F_{\mu\nu}F^{\mu\nu}
 +\bar\psi i\gamma^\mu(\partial_\mu+i e_D A_\mu)\psi
 +\cdots.                                             \tag{EU06}
\]

The same conclusion follows with a canonical charged scalar.  Charge
normalization, the visible eigenmode, and the renormalization convention are
fixed before applying the theorem; without charged matter the normalization
of a free Abelian field does not define a physical alpha.

Freeze the selector premise:

> **`MARGINAL-ALPHA`.**  The physical inherited charge `e_D` is classically
> marginal, and `alpha=C e_D^2` with dimensionless convention factor `C`
> requires no additional length, mass, cutoff, compactification, density, or
> renormalization scale to cancel an engineering dimension.  Quantum running
> along the inherited trajectory is allowed; an engineering power of the RG
> scale is not being hidden in the definition of alpha.

This is stronger than merely saying that some scale-dependent symbol called
alpha can be made dimensionless.

## 4. Theorem AQ4DL-1 -- Maxwell engineering-dimension lock

Canonical normalization of (EU06) gives

\[
 [A_\mu]={D-2\over2},
 \qquad [\psi]={D-1\over2}.                         \tag{EU07}
\]

Because `partial_mu` and `e_D A_mu` occur in the same covariant derivative,

\[
 \boxed{[e_D]=1-[A_\mu]={4-D\over2}},
 \qquad [\alpha]=4-D.}                              \tag{EU08}
\]

Therefore

\[
 \boxed{
 \operatorname{MARGINAL\!\!\!-\!ALPHA}
 \Longrightarrow [e_D]=0
 \Longrightarrow D=4.}                             \tag{EU09}
\]

### Proof

The Maxwell kinetic term has dimension
`2+2[A_mu]=D`, which gives the first equation in (EU07).  The Dirac kinetic
term has dimension `1+2[psi]=D`, which gives the second.  Equality of the
dimensions of `partial_mu` and `e_D A_mu` gives (EU08).  A classically
marginal charge has `[e_D]=0`; hence `(4-D)/2=0` and `D=4`.  QED.

The normalization-equivalent convention

\[
 \mathcal L_{\rm vis}
 =-{1\over4g_D^2}\mathcal F_{\mu\nu}\mathcal F^{\mu\nu}
 +\bar\psi i\gamma^\mu(\partial_\mu+i\mathcal A_\mu)\psi,
 \qquad \mathcal A_\mu=g_DA_\mu                       \tag{EU10}
\]

has `[mathcal A_mu]=1` and

\[
 \boxed{[g_D^2]=4-D,}                                \tag{EU11}
\]

so the conclusion is convention independent.

## 5. Theorem AQ4DL-2 -- alpha selects q=4 inside the admitted class

If `QFRONT-DIM`, the same-front Maxwell premise, and `MARGINAL-ALPHA` all
hold, then (EU05) and (EU09) compose:

\[
 \boxed{
 \operatorname{QFRONT\!\!\!-\!DIM}
 \land\operatorname{MARGINAL\!\!\!-\!ALPHA}
 \Longrightarrow D=4
 \Longrightarrow q=4.}                             \tag{EU12}
\]

This is the precise ALLOW/REQUIRE statement:

- record/front algebra may **allow** other `q` values;
- a same-front canonically normalized visible Maxwell sector with a
  scale-free classically marginal inherited alpha **requires** `q=4` within
  this class.

It removes q=4 as a free numerological choice only after the electromagnetic
premises are earned on the same record front.  It neither produces the four
physical operations nor calculates the numerical value of alpha.

## 6. Dimensionless running couplings in `D!=4`

For `D!=4`, one may always define

\[
 \widehat e(\mu)=e_D\,\mu^{(D-4)/2},
 \qquad
 \widehat\alpha(\mu)=C_D e_D^2\mu^{D-4}.             \tag{EU13}
\]

These quantities are dimensionless, but their definition contains the
compensating scale `mu`.  Before loop or anomalous contributions,

\[
 \mu{d\widehat e\over d\mu}={D-4\over2}\widehat e,
 \qquad
 \mu{d\widehat\alpha\over d\mu}=(D-4)\widehat\alpha. \tag{EU14}
\]

An interacting fixed point can cancel the canonical term in a specific
`D!=4` theory.  That is a different premise: it does not make the canonical
Maxwell charge classically marginal, and it introduces a fixed-point dynamics
that this theorem does not exclude or derive.  Likewise, a lattice spacing,
compactification radius, density, mass, or cutoff can convert a dimensionful
`e_D` into a dimensionless ratio.  Such a compensator violates
`MARGINAL-ALPHA` rather than contradicting AQ4DL-1.

## 7. Same-world alpha inheritance boundary

The pinned SAI theorem says that, once one visible canonically normalized
`U(1)` parent and comparison context are physically fixed, compatible records
in that sector inherit one alpha RG trajectory and a genuinely different
alpha cannot be inserted as an ordinary same-sector subsystem.  SAI does not
derive the `U(1)` sector, its spacetime dimension, or classical marginality.

Consequently AQ4DL and SAI compose only as follows:

\[
 \begin{split}
 &\text{same record front earns `QFRONT-DIM` and a visible canonical Maxwell}\
 &\text{sector satisfying `MARGINAL-ALPHA`}\\
 &\qquad\Longrightarrow q=D=4,\\
 &\text{SAI same-sector ancestry}
 \qquad\Longrightarrow\text{that dimension lock and alpha trajectory are}\
 &\hspace{43mm}\text{shared across the admitted visible sector.}
 \end{split}                                                   \tag{EU15}
\]

This is a conditional actual-world selector.  It is not a derivation of the
Maxwell sector from RFT and not a proof about hypothetical sectors that fail
the premises.

## 8. Controls and falsifiers

1. **Scale-compensator control.**  In `D!=4`, call (EU13) an intrinsically
   scale-free alpha while omitting `mu`.  Dimensional analysis falsifies the
   claim.
2. **Free-photon control.**  Remove all charged matter and charge
   normalization.  Rescaling the free Abelian field makes `e_D` nonphysical;
   no alpha selector follows.
3. **Noncanonical control.**  Replace the Maxwell kinetic term with a
   higher-derivative, nonlocal, Lifshitz, or anomalously scaled operator.
   Equation (EU07) must be recomputed and AQ4DL does not apply unchanged.
4. **Fixed-point control.**  Use an interacting `D!=4` fixed point to cancel
   the canonical beta term.  This may define a dimensionless fixed-point
   coupling, but it does not satisfy the stated classical-marginality premise.
5. **Compactification control.**  Let `q>D` with extra contrast directions
   compactified, gapped, or gauge.  Then `QFRONT-DIM` fails; the theorem selects
   only the active infrared dimension, not the microscopic count rank.
6. **Circular-Maxwell control.**  Import a Maxwell theorem already assuming
   `3+1` dimensions to prove `D=4`.  AQ4DL instead uses the general-`D`
   engineering dimensions (EU07)--(EU08).
7. **Different-front control.**  Take alpha from an unrelated continuum model
   and q from the record front.  The same-front premise fails.
8. **Numerical-alpha control.**  Infer `alpha_obs` from `D=4`.  Marginality
   fixes only an engineering dimension, not a coupling value or boundary
   condition.
9. **Gravity control.**  Relabel the q=4 selector as a metric, curvature,
   protected spin-two pole, Einstein response, or gravity.  None follows here.

## 9. Exact disposition

This lane proves that a supplied q-operation count/contrast spacetime has
`D=q` under the explicit no-extra/no-collapse physical dimension premise, and
that a same-front canonical Maxwell charge has engineering dimension
`(4-D)/2`.  Requiring the inherited visible alpha to be classically marginal
without a compensating length/mass scale forces `D=4` and therefore `q=4`.
Dimensionless scale-dressed couplings and interacting fixed points in other
dimensions are explicitly outside that inference.  The result is a narrow
Alpha-to-q4 REQUIRE-side selector for the same-world gravity-origin chain, not
a derivation of electromagnetism, alpha's value, or gravity.

**Disposition:**

`GENERAL_Q_COUNT_CONTRAST_DIMENSION_D_EQUALS_Q_CONDITIONAL__CANONICAL_MAXWELL_CHARGE_DIMENSION_(4_MINUS_D)_OVER_2_EXACT__CLASSICALLY_MARGINAL_SCALE_FREE_INHERITED_ALPHA_REQUIRES_D4_AND_Q4__D_NOT_EQUAL_4_SCALE_DRESSED_AND_FIXED_POINT_COUPLINGS_NOT_EXCLUDED__SAME_FRONT_EM_AND_QFRONT_PREMISES_OPEN__ALPHA_VALUE_AND_GRAVITY_NOT_DERIVED`

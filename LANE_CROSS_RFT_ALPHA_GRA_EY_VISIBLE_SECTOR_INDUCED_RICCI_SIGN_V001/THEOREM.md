# Visible-sector induced Ricci-sign theorem

**Lane ID:** `CROSS-RFT-ALPHA-GRA-EY-VISIBLE-SECTOR-INDUCED-RICCI-SIGN-V001`

**Official short name:** `VSIRS`

**Date:** 2026-08-27

**Builder status:** `MUTABLE_PROOF_DRAFT__NOT_SEALED`

**Claim class:** exact one-loop proper-time coefficient census for the declared
unbroken visible Standard-Model field content; exact open sign condition with
nonminimal Higgs and additional pair-memory scalar modes; conditional RIEHB
coefficient calculation in a physically distinguished zero-bare-term
E-EMERGENTSPACE matching prescription

**Not claimed:** that the visible Standard Model is the complete ultraviolet
spectrum; that the record crossover equals the Planck scale; a
scheme-independent absolute Newton constant without a distinguished parent
matching rule; a small cosmological constant; physical record-to-metric
soldering; continuum emergence; gravity; or outcome selection

## 1. Exact question

RIEHB left the sign of the induced Ricci coefficient as a physical spectrum
gate.  That is essential for an arbitrary collection of fields, but our actual
visible sector is not an arbitrary collection.  Does its known spin census
already determine the sign of its one-loop contribution in the Sakharov
proper-time prescription?

Conditionally, yes, but the actual visible census is close to a cancellation
rather than sign-definite by spin count alone.  In the convention fixed below,
the exact heat-kernel supertrace is negative precisely on the open domain
\(\xi_H>1/24\), so the induced Einstein--Hilbert coefficient is positive
there.  Minimal Higgs coupling gives the opposite sign; table-conformal
coupling gives the positive sign.  Additional minimal bosonic pair-memory
modes tighten rather than relax the positivity bound.

This closes the sign of a concrete visible-field contribution conditional on
the owned \(\xi_H\) domain.  It closes the total RIEHB coefficient only if the
same parent proves that the declared shell census is complete and fixes the
zero-bare-term matching prescription.

## 2. Frozen coefficient and curvature convention

Use the one-loop proper-time convention summarized in Visser's induced-gravity
coefficient table; call it convention \(V\).  Visser writes the standard
Lorentzian gravity term as

\[
 -{1\over16\pi G}\int\sqrt{-g}\,R_V .
\]

RIEHB writes it as \(+(16\pi G)^{-1}\int\sqrt{-g}\,R_R\).  The frozen
dictionary is therefore

\[
 \boxed{R_R=-R_V,\qquad
 -{R_V\over16\pi G}=+{R_R\over16\pi G}.}           \tag{EY00}
\]

All Laplacian, Wick-rotation, and nonminimal-coupling signs are first evaluated
in convention \(V\), and only the completed invariant is then translated by
(EY00).  In particular, \(\xi_H\) in this lane means the coefficient in
Visser's scalar Hessian/table convention, defined by
\(k_1^V=1/6-\xi_H\).  A coefficient written in another Riemann or Hessian
convention must be mapped before insertion; reusing the symbol without that
map is forbidden.

For one complete species, including physical spin states and the required
gauge/ghost census, the convention-\(V\) Ricci heat-kernel supertrace
coefficients are

\[
 k_1(\hbox{real scalar})={1\over6}-\xi,qquad
 k_1(\hbox{Weyl spinor})=-{1\over6},qquad
 k_1(\hbox{massless vector})=-{2\over3}.            \tag{EY01}
\]

These are traced heat-kernel coefficients before the statistics weight.
The supertrace is
\[
 \operatorname{str}k_1
 :=\sum_A(-1)^{F_A}n_A k_{1,A},                    \tag{EY01a}
\]
so the Weyl table value \(-1/6\) acquires a second minus from
\((-1)^F\), while the bosonic massless-vector value \(-2/3\) does not.
Gauge/ghost accounting internal to the quoted vector coefficient must not be
applied a second time.  For the sharp proper-time window
\(s\in[\kappa_R^{-2},\mu^{-2}]\), in its massless asymptotic regime and with a
physically distinguished zero bare Ricci term,

\[
 {1\over G_{>}}
 =-{\operatorname{str}k_1\over2\pi}
   (\kappa_R^2-\mu^2),                              \tag{EY02}
\]

or equivalently

\[
 \boxed{
 C_R^{>}:={1\over16\pi G_{>}}
 =-{\operatorname{str}k_1\over32\pi^2}
   (\kappa_R^2-\mu^2).}                            \tag{EY03}
\]

The first equality is Visser's convention-\(V\) one-loop-dominance result; the
second follows from (EY00), so \(C_R^>\) is precisely RIEHB's positive-\(R_R\)
coefficient.  Equations (EY02)--(EY03) use \(\hbar=c=1\).  The proper-time
window is a covariant regulator representative, not an exact momentum
projector.  Mass thresholds replace the
common quadratic shell factor by the owned positive threshold integrals of
RIEHB.  Once different species have different thresholds, their terms must be
summed separately; the massless simplification below may not be used.

Primary coefficient source: Matt Visser,
[“Sakharov's induced gravity: a modern perspective”](https://arxiv.org/abs/gr-qc/0204062),
especially its gravity-action convention, equations (21) and (29), and its
low-spin \(k_1\) table.

## 3. Declared visible-sector census

Use the zero-temperature ultraviolet/symmetric-gauge-basis coefficient in a
proper-time window whose infrared end obeys
\(\mu\gg v_{\rm EW},m_i\) for every field treated as massless.  This is a
high-energy asymptotic statement, not an assumption of a thermal cosmological
electroweak-restored phase.  Take the minimal visible field content:

1. one complex Higgs doublet, hence four real scalar components with one common
   curvature coupling \(\xi_H\);
2. the complete three-generation anomaly-free Standard-Model chiral
   representation set: fifteen two-component Weyl fields per generation,
   hence forty-five Weyl species, without right-handed neutrinos; and
3. twelve massless gauge vectors: eight color, three weak, and one hypercharge.

The hypercharge vector \(B_\mu\) and neutral weak vector \(W^3_\mu\) are the
two gauge ancestors of the canonically normalized photon,
\(A_\mu=\cos\theta_W B_\mu+\sin\theta_W W^3_\mu\) in the displayed
convention.  Low-energy \(e\) and the running-alpha trajectory inherit both
\(g_Y\) and \(g_2\) through electroweak mixing.  Alpha does not change the
spin-count coefficient (EY01); it remains part of the same-parent field,
threshold, and complete-stress census.

The Weyl term is used only for the complete anomaly-free representation set,
not as forty-five unrelated chiral determinants.  The Standard-Model
hypercharges cancel the perturbative gauge and mixed gauge--gravitational
anomalies generation by generation, and the number of left-handed
\(SU(2)\) doublets is even, so the parity-even squared/Laplace-type determinant
census used by the quoted \(k_1\) coefficient is admitted.  If a proposed
extension is anomalous or cannot be reduced to the declared determinant
class, its Weyl coefficient may not be added by raw counting.

Adding three right-handed-neutrino Weyl modes, if they lie in the same shell,
makes the supertrace more positive and the induced-only Ricci sign less
favorable.  Additional fields are treated explicitly in section 5 rather than
silently excluded.  If the proper-time window crosses electroweak breaking or resolved
mass thresholds, the common factor below is invalid: species-specific
threshold integrals and the broken-phase vector, Goldstone, gauge-fixing, and
ghost ledger must replace it.

## 4. Theorem VSIRS-1 -- exact visible-sector sign

For the field content in section 3,

\[
\begin{aligned}
 \operatorname{str}k_1^{\rm vis}
 &=4\left({1\over6}-\xi_H\right)
   -45\left(-{1\over6}\right)
   +12\left(-{2\over3}\right)\\
 &=\boxed{{1\over6}-4\xi_H}.                       \tag{EY04}
\end{aligned}
\]

Therefore

\[
 \boxed{
 \xi_H>{1\over24}
 \quad\Longrightarrow\quad
 C_{R,\rm vis}^{>}
 ={24\xi_H-1\over192\pi^2}
 (\kappa_R^2-\mu^2)>0.}                            \tag{EY05}
\]

Table-minimal coupling \(\xi_H=0\) lies outside this domain and gives

\[
 \boxed{
 \operatorname{str}k_1^{\rm vis}={1\over6},\qquad
 C_{R,\rm vis}^{>}
 =-{1\over192\pi^2}(\kappa_R^2-\mu^2)<0.}         \tag{EY06}
\]

Table-conformal coupling \(\xi_H=1/6\) lies strictly inside the domain and
gives

\[
 \boxed{
 \operatorname{str}k_1^{\rm vis}=-{1\over2},\qquad
 C_{R,\rm vis}^{>}
 ={1\over64\pi^2}(\kappa_R^2-\mu^2)>0.}           \tag{EY06a}
\]

### Proof

Insert the complete species multiplicities into (EY01) and collect over a
common denominator to obtain (EY04).  Since the shell factor in (EY03) is
positive, \(C_R^{>}>0\) exactly when
\(\operatorname{str}k_1<0\).  Solving that inequality gives (EY05), and
the two displayed values give (EY06)--(EY06a). QED.

## 5. Theorem VSIRS-2 -- pair-memory and unknown-field margin

Let \(N_p\) additional minimally coupled real scalar modes belong to the same
fast shell and be counted nowhere else.  This includes the conservative
control in which the six PMMDC pair-memory deformation coordinates acquire
ordinary minimal scalar propagation.  Then

\[
 \operatorname{str}k_1^{\rm vis+p}
 ={1+N_p\over6}-4\xi_H,                           \tag{EY07}
\]

and

\[
 \boxed{
 24\xi_H>1+N_p
 \quad\Longrightarrow\quad
 C_{R,\rm vis+p}^{>}
 ={24\xi_H-1-N_p\over192\pi^2}
 (\kappa_R^2-\mu^2)>0.}                            \tag{EY08}
\]

For six minimal pair-memory modes and a minimally coupled Higgs,

\[
 \boxed{
 \operatorname{str}k_1^{\rm vis+p6}={7\over6},\qquad
 C_{R,\rm vis+p6}^{>}
 =-{7\over192\pi^2}(\kappa_R^2-\mu^2)<0.}        \tag{EY09}
\]

After the statistics weight, an added Weyl mode contributes \(+1/6\) to the
supertrace and can spoil positivity; an added standard massless vector
contributes \(-2/3\) and helps it.  Either extension is admitted only as part
of a lawful anomaly-free set in the same standard Laplace-type determinant
class.
An arbitrary additional scalar with nonminimal coupling contributes
\(1/6-\xi\), so the exact general condition is the corresponding total
\(\operatorname{str}k_1<0\); raw species count is not sufficient.

### Proof

Each added minimal real scalar contributes \(1/6\) to (EY04).  Substitution in
(EY03) gives (EY07)--(EY09). QED.

## 6. Conditional Newton-scale corollary

If the declared census is the complete fast shell, E-EMERGENTSPACE supplies a
physically distinguished matching rule with no pre-existing metric Ricci term,
all RIEHB measure/boundary/remainder conditions hold, and
\(24\xi_H>1+N_p\), then the visible-sector calculation predicts

\[
 \boxed{
 G_{\rm eff}(\mu)
 ={12\pi\over
 (24\xi_H-1-N_p)(\kappa_R^2-\mu^2)}}              \tag{EY10}
\]

in natural units on the admitted leading-derivative band.  Thus the remaining
dimensionful input is the independently derived record/geometry crossover
\(\kappa_R\), not accepted Newton \(G\).  Conversely, inserting measured
\(G\) to define \(\kappa_R\) is a fit, not an origin proof.

Equation (EY10) is not promoted for the actual world until the complete
spectrum, threshold, regulator, and matching custody pass.  It does establish
an exact calculation target: an independently measured or derived
\(\kappa_R\) and complete field census would determine \(G_{\rm eff}\) rather
than merely its sign.

## 7. RIEHB and record-ancestry interface

Under the \(24\xi_H>1+N_p\), complete-census, convention-dictionary, shell, and
distinguished-matching premises, VSIRS closes RIEHB's
`C_R^{eff}>0` gate for the declared parent.  It does not close:

1. same-parent record-query localization and physical soldering;
2. q4 front instantiation, shared-edge gluing, or continuum refinement;
3. the six-mode spatial deformation and initial-constraint gates;
4. complete stress/work/controller/boundary matching in a lineage
   intervention;
5. curvature-squared/nonlocal remainder bounds outside the declared band; or
6. the volume term and cosmological-constant problem.

The electromagnetic/alpha sector participates twice without double counting:
it helps organize actual record formation and matter, and its gauge field is
one member of the complete fast/slow action whose metric determinant and
stress contribute to gravitational response.  Alpha itself is not an extra
stress source and does not replace the complete action derivative.

## 8. Controls and exact disposition

1. **Incomplete-spectrum control.**  Call the visible contribution the total
   coefficient while unowned ultraviolet fields or matching terms remain.
   Rejected.
2. **Scheme control.**  Claim an absolute induced \(G\) while allowing an
   arbitrary bare/counterterm shift.  Rejected; the distinguished
   E-EMERGENTSPACE matching prescription is load bearing.
3. **Threshold control.**  Apply one common massless shell factor across
   resolved unequal masses.  Rejected.
4. **Gauge control.**  Count vector polarizations without the gauge/ghost
   coefficient in (EY01).  Rejected.
5. **Higgs control.**  Omit \(\xi_H\).  Rejected; (EY05) displays its exact
   open sign range.
6. **Pair-mode control.**  Add six scalar determinants before proving that the
   PMMDC deformation modes propagate as independent minimally coupled fields.
   Rejected; (EY09) is a conditional control and shows that such modes worsen
   the induced-only sign bound.
7. **Scale-fit control.**  solve (EY10) for \(\kappa_R\) using accepted \(G\)
   and call it a prediction.  Rejected.
8. **Lambda control.**  infer a small cosmological term from the positive
   Ricci coefficient.  Rejected.
9. **Gravity control.**  relabel one positive induced contribution as the
   final record-origin theorem without soldering, refinement, constraints,
   complete stress, and ancestry.  Rejected.
10. **Convention control.**  import the sign of Visser's \(R_V\) coefficient
    into RIEHB's \(R_R\) action without (EY00), or insert an unmapped
    nonminimal coupling.  Rejected.
11. **Electroweak-window control.**  call the ultraviolet gauge-basis census a
    thermal restoration theorem, or use its common massless factor across a
    window that resolves electroweak thresholds.  Rejected.
12. **Chiral control.**  sum isolated Weyl determinants without the complete
    anomaly-free representation and determinant-square premises.  Rejected.
13. **Statistics control.**  treat the Weyl table value \(-1/6\) as already
    fermion-weighted and then apply no explicit \((-1)^F\), or apply that
    weight twice.  Rejected; (EY01a) fixes the single statistics insertion.

**Disposition:**

`DECLARED_VISIBLE_GAUGE_BASIS_HAS_EXACT_SUPERTRACE_K1_ONE_OVER_6_MINUS_4_XI_H__POSITIVE_INDUCED_RICCI_COEFFICIENT_IFF_XI_H_GREATER_THAN_ONE_OVER_24_IN_THE_DECLARED_ZERO_PAIR_MODE_CENSUS__MINIMAL_HIGGS_SIGN_NEGATIVE__TABLE_CONFORMAL_HIGGS_SIGN_POSITIVE__SIX_MINIMAL_PAIR_MEMORY_SCALARS_TIGHTEN_BOUND_AND_GIVE_NEGATIVE_SIGN_AT_MINIMAL_HIGGS__CONDITIONAL_G_EFF_FORMULA_IN_TERMS_OF_INDEPENDENT_RECORD_CROSSOVER__COMPLETE_UV_SPECTRUM_THRESHOLDS_DISTINGUISHED_MATCHING_SOLDERING_REFINEMENT_CONSTRAINTS_STRESS_ANCESTRY_LAMBDA_AND_REAL_WORLD_GRAVITY_OPEN`

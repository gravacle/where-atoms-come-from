# Independent hostile audit: GL6BX V001

**Target:** `LANE_CROSS_RFT_GRA_GL6BX_EIGHTEEN_LINK_ACTIVE_COLLAR_FULL_PAIR_SOURCE_V001`  
**Audit:** `AUDIT_G_GL6BX_EIGHTEEN_LINK_ACTIVE_COLLAR_FULL_PAIR_SOURCE_V001`  
**Date:** 2026-09-02  
**Verdict:** `PASS_AFTER_PRESCREEN_AND_HBAR_REPAIR_RESEAL`  
**Target edits by auditor:** none

## 1. Verdict

The frozen GL6BX theorem is correct within its declared finite-collar,
source-before-projection, canonical two-state, and order-`h^6` scope.  The
eighteen-link geometry, complete source-free sixth-order matrix, nonzero but
commuting `T2` first vertex, source contacts, literal four-class cancellation,
transition normalization, CTP sign/scale separation, contact root exclusion,
and truncated stationary root exclusion all pass hostile reconstruction.

The decisive result is a narrowly scoped no-single-collar-tuning theorem:
neither the complete displayed contact polynomial nor the summed `h<=6`
stationary two-level surrogate reaches `h_E/2=h_T` at positive `r`.  It is not
a no-go for the untruncated parent, connected accumulation, or gravity.

## 2. Prescreen repairs and frozen custody

The hostile prescreen and units review found five issues before the final
target was frozen.

1. An early regression guard applied the pure-`T2` GL6BV onsite identity to
   `A1/E2` sources.  That is false because `A1/E2` act nontrivially inside
   `P`.  The author restricted the guard to `T2` while retaining the full
   defect-minus-locked source terms for every sector.
2. The first prose snapshot did not explicitly distinguish authenticated
   orientation classes from the equal one-per-class weights.  The repaired
   theorem and ledger state that this finite composition is prospectively
   imposed, not generated or weighted by the parent.
3. `EXACT_LEDGER.json` was initially outside manifest coverage and the packet
   verifier retained a stale `504/504` ledger token after the replay grew to
   `512/512`.
4. The first reseal retained a line-sensitive scope token that failed across
   a theorem line break.  The author repaired the verifier and resealed in
   content-manifest-seal order.
5. A later CTP review found that the normalized-time formula had restored
   `Delta t/hbar` in the sine but omitted the overall `1/hbar`; its small-time
   coefficient consequently also missed `1/hbar^2`.  The author restored the
   complete physical convention and performed a fresh full source replay and
   reseal.

The final audit pins the corrected snapshot only.  Its dependencies, all
fifteen payload bytes, manifest, and seal replay cleanly.  The auditor made no
target edit.

## 3. Geometry and two-word retained space

The audit rebuilt the period-four bipartite incidence graph, the canonical
alternating hexagon, and the deterministic degree-two `Q4` witness without
importing a GL6BX module.  Taking the union of all links incident on the six
cycle nodes gives exactly six cycle links and twelve spokes.  The six core
nodes have all four incident flips active; twelve boundary nodes have one
active link each.  Their frozen occupancies split six with one and six with
two.

An exhaustive `2^18` census gives exactly two zero-defect active words.  Their
symmetric difference is precisely the six cycle links.  Thus `dim P=2` is an
output of the declared active collar, not an imported isolated-hexagon
assumption.

## 4. Wave operator, Hermitian transform, and source derivatives

The author solves

\[
 C+D\chi=\chi(A+B\chi),\qquad
 H_B=A+B\chi,\qquad G=I+\chi^\dagger\chi,
 \qquad H_C=G^{1/2}H_BG^{-1/2}.                         \tag{A-BX01}
\]

The audit traced the coefficient recurrence and both source derivatives term
by term.  The first Bloch derivative retains the direct `C_a,D_a,A_a,B_a`
terms, both differentiated nonlinear factors, and the unknown `chi_a` at the
correct order.  The mixed derivative contains the two derivative orderings of
every product and all `chi_ab` insertions.  No fold is omitted.

For `S=G^{1/2}`, the implementation uses the noncommutative Sylvester rules

\[
 SS_a+S_aS=G_a,
 \qquad SS_{ab}+S_{ab}S=G_{ab}-S_aS_b-S_bS_a,           \tag{A-BX02}
\]

and differentiates `S^{-1}` and `S H_B S^{-1}` with every first/mixed product
term retained.  Exact Hermiticity checks at each order are therefore
substantive rather than a final symmetrization.

As an independent algorithmic cross-check, the audit used an
intermediate-normalized stationary-branch Rayleigh--Schrodinger recurrence,
not the author's Bloch/des-Cloizeaux code.  It gives

\[
 E_+^{(6)}=-{3511\over360},\qquad
 E_-^{(6)}={2159\over360},                               \tag{A-BX03}
\]

so their half-sum is `-169/90` and half-difference is `-63/8`.  Together with
the independently recovered lower orders,

\[
 H_C=(-9r^2-11r^4/4-169r^6/90)I-(63r^6/8)\sigma_x.      \tag{A-BX04}
\]

This explicitly guards the often-dropped sixth-order diagonal term as well as
the `GL6AO` off-diagonal coefficient.  The stationary gap is `63r^6/4`.

## 5. Source convention and nonzero `T2` vertex

For a nontrivial local `T1` source at sorted core node one, the independent
eigenbranch recurrence gives

\[
 V_T^{(4)}={79\over36}I,
 \qquad V_T^{(6)}={33307\over6480}I+{105\over8}\sigma_x, \tag{A-BX05}
\]

and the onsite second derivative

\[
 C_T^{(2)}=-4I,\quad C_T^{(4)}=-{187\over54}I,
 \quad C_T^{(6)}=-{184387\over19440}I-{165\over4}\sigma_x. \tag{A-BX06}
\]

The literal core-node/port transforms reproduce the reported pattern
`(-T3,+T1,+T2,+T2,+T1,-T3)`.  Every nonzero first vertex is a combination of
`I` and `sigma_x`; it commutes with (A-BX04) and has no retained-doublet
transition.  This is a stationary/contact `T2` ingredient, not a `T2` pole.

The `E2` source cannot use the pure-`T2` GL6BV shortcut.  At node zero, for
example, its locked-space first vertex already has order-zero matrix
`diag(4,0)`.  The repaired Bloch source equation keeps `D_a chi-chi A_a`, and
the independent order-`h^2` calculation uses
`delta M=M(defect)-M(locked)` in every denominator derivative.  A separate
full `2^18` sparse finite-difference/polar-canonicalization prescreen at
`r=0.15` matched the source-free, local `T1`, and inside-`P` local `E1` jets
through `h^6`, with the remaining discrepancy scaling as `h^8` and residuals
below `9e-8`.

## 6. Spatial moments and orientation multiplicity

The audit independently contracts the six source positions with

\[
 {1\over3}\left[\sum_i X_iX_i-
 {1\over2}(RC+CR)\right]                                \tag{A-BX07}
\]

for the normalized `a_*^2 k^2` coefficient.  At order `h^2`, the selected
`E/T` block is proportional to

\[
 \begin{pmatrix}1&1&2\\1&-1&0\end{pmatrix},             \tag{A-BX08}
\]

with scalar `-4`.  A literal permutation/translation census produces four
orientation classes, each of stabilizer six.  The selected block is invariant
under its own stabilizer, and one representative from each class sums to zero.
The cancellation is therefore checked on physical transforms, not inferred
from an abstract basis average.

The transition-gradient normalization is separately

\[
 R_{E,P}=8,\qquad R_{E,C}=8,\qquad
 R_{E,P+C}=32,
 \qquad R_{E,4\,selected\ classes}=128.                 \tag{A-BX09}
\]

The parent and child gradients are identical, so adding all six nodes doubles
the amplitude and quadruples the residue.  The final factor four is the class
multiplicity.  These are different factors; there is no second orientation
factor hidden in `32`.

The same independent order-`h^2` ledger gives

\[
 (h_A,h_E,h_T)_{k^0}=(-24,-72,-56),\qquad
 (h_A,h_E,h_T)_{k^2}=(4,4,4).                            \tag{A-BX10}
\]

Thus the leading contact solder defect is `4/2-4=-2`, and the nonzero `k^0`
contact cannot be canceled by the vanishing retained-doublet `k^0` residue.
No Ward identity follows.

Crucially, equal one-per-class multiplicity in (A-BX09) is a selected finite
composition.  GL6BX authenticates the four orientations but does not derive
their weights from the parent, refinement, stationarity, or accumulation.
Consequently its diagonal four-class result is not a macroscopic isotropy
theorem.

## 7. CTP separation

With the inherited physical half-susceptibility convention,

\[
 K^R(t)={i\over2\hbar}\theta(t)\langle[V(t),V(0)]\rangle,
 \qquad
 K^R_{\rm dyn}(t,k)={\theta(t)\over\hbar}R(k)
 \sin\!\left({\Delta t\over\hbar}\right).              \tag{A-BX11}
\]

This is one-half the standard raw susceptibility for a source Hamiltonian
`H+jV`.  At `k^0`, `R=0`.  At `k^2`, only `E2` is nonzero and begins at `128`;
together with `Delta=63U_dr^6/4`, this yields the small-time coefficient
`2016U_dr^6 t/hbar^2` and solder defect `1008U_dr^6 t/hbar^2`.  The regulated
zero-frequency normalization is `R/Delta`, so the stationary convention is
`F''=-2K^R(0)`.  The eliminated-`Q` energy Hessian instead appears in the
connected CTP functional as `-C delta(t)`.  The two supports must not be
equated, and restoring the microscopic high-gap response replaces the local
contact rather than adding a second copy.

## 8. Exact root exclusions

For the contact defect, the quadratic bracket in `x=r^2` has

\[
 \operatorname{disc}=-{236225414\over10125},\qquad
 x_{\rm vertex}={5521500\over2146938851},\qquad
 D_{\rm max}=-{4252057452\over2146938851}.               \tag{A-BX12}
\]

It is strictly negative for every positive `r`.

For the stationary surrogate, define the analytic scaled polynomial
`Q_6(x)=r^6D_6(r)` with `x=r^2`.  Clearing denominators and content gives the
ascending coefficients

```text
(-16796160000, 33592320000, 49921920000, 105059462400,
 -242165732400, -309862380600, -6525772091869).
```

The independent exact rational Sturm chain has degrees
`(6,5,4,3,2,1,0)`.  Its signs at `0+` and positive infinity are

```text
(-,+,+,-,+,+,+)   V=3
(-,-,+,+,-,-,+)   V=3.
```

Hence `Q_6` has zero positive roots.  Since `Q_6(0)<0`, the truncated
stationary defect is negative for every `r>0`, including the conservative
isolated-doublet domain `0<r<1/18`.

This is exact only for the summed `h<=6` two-level functional.  Dividing the
spectral term by a gap beginning at `r^6` means unknown order-eight transition
vertices and gap corrections already alter nonleading Laurent coefficients.
There is therefore no controlled full-parent threshold without a joint
transition-vertex/contact/gap expansion and a rigorous remainder bound.

The narrow future certification test is as follows.  Scale the complete
defect to an analytic function of `x`; isolate a candidate positive root by
exact Sturm arithmetic; keep its interval inside a proved spectral/source
domain; and bound the full remainder `R`.  Endpoint bounds
`|R(x_i)|<|Q(x_i)|` preserve a sign bracket.  The derivative condition
`sup|R'|<inf|Q'|` on the interval proves uniqueness, and
`|delta x|<=sup|R|/inf|Q'|` controls displacement.  Existing exact
coefficients and the Sturm implementation can be reused, but the higher
transition, contact, and gap jets must be added together.

## 9. Live-path consequence

GL6BX closes tuning of this one selected collar through the displayed order.
It does not close owner-correct connected sums, and it does not establish that
the parent supplies equal orientation weights.  The next lawful route is an
owner-correct connected accumulation of authenticated collars under one
microscopic source, stationary state, and boundary rule.  No Ricci tensor,
graviton, gravity law, or `G` is inserted or derived here.

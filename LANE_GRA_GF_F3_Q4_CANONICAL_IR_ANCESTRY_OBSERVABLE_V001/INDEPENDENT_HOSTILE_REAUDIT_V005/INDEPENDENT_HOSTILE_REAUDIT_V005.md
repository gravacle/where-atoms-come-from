# Independent hostile re-audit of repaired GF V005

**Target:** `LANE_GRA_GF_F3_Q4_CANONICAL_IR_ANCESTRY_OBSERVABLE_V001`  
**Target revision:** `V005`, root-repaired after the preserved V001--V004 rejections  
**Verdict:** **PASS -- NARROW DESIGN-CONTRACT SEAL**  
**Independent replay:** `270/270` audit checks passed  
**Author replay:** `104/104` reproduced

## Executive finding

No material defect remains in the frozen V005 prospective C1/G2 observable
contract.  The canonical-amplitude classifier is now disjoint and total over
the declared physical limit data: raw rank, availability of an applicable
same-parent canonical map, positive lower singular-value bound, finite upper
singular-value bound, and arbitrary raw scaling when no map exists.  Its exact
executable classifier agrees with that partition on the complete Boolean
truth table.

All three prior classifier defects remain repaired:

1. canonical amplitude and normalized ancestry are separate outputs;
2. every full-rank raw no-map case is INDETERMINATE, independent of raw
   scaling; and
3. an applicable map with unbounded `sigma_max(Z_Q^can)` is FAIL.

This is a seal of the **design contract only**.  V005 contains no matched-
family spectral payload and proves no positive G2, G3, gravity, or `G` result.

## Exhaustive amplitude decision

The frozen classifier is equivalent to

\[
\begin{array}{c|c|c}
\text{condition}&\text{amplitude label}&\text{reason}\\ \hline
\operatorname{rank}\mathcal N^{\rm raw}<2&\mathrm{FAIL}&
  \text{persistent source-channel rank loss}\\
\operatorname{rank}\mathcal N^{\rm raw}=2,\ D\ \text{absent}&
  \mathrm{INDETERMINATE}&\text{canonical normalization unavailable}\\
D\ \text{applicable},\ \liminf\sigma_{\min}>0,\
  \limsup\sigma_{\max}<\infty&\mathrm{PASS}&\text{finite nonzero vertex}\\
D\ \text{applicable and rank/lower loss or upper divergence}&
  \mathrm{FAIL}&\text{canonical amplitude criterion fails}.
\end{array}
\]

Independent enumeration covered all 16 combinations of raw-rank status, map
availability, lower-limit status, and upper-limit status.  Exactly one branch
fires in every case.  The exact author function was then tested on the full
32-row table obtained by independently varying the now-non-load-bearing raw
`controlled_power` flag; every output matched the specification.

## Exact regression families

### V004 upper-bound regression -- PASS

For

\[
 \mathcal N^{\rm raw}_{L}=I_2,\qquad
 D_L=\operatorname{diag}(L,1),\qquad
 Z^{\rm can}_{Q,L}=\operatorname{diag}(L,1),
\]

the vertex remains rank two with `sigma_min=1` while `sigma_max=L` is
unbounded.  V005 uniquely returns FAIL.

### Lower-bound loss -- PASS

For `D_L=diag(1/L,1)`, the finite-size vertex remains rank two but its lower
singular value tends to zero.  V005 uniquely returns FAIL.

### V003 no-map regression -- PASS

Full-rank raw numerators proportional to `L^{-1}I`, `I`, and `LI` all return
INDETERMINATE when no applicable same-parent `D` exists.  Thus neither
vanishing, constant, nor divergent raw scaling is confused with a canonical
PASS or physical FAIL.

### V002 ancestry separation -- PASS

For `Delta_L=1/L`, `R_L=(L/2)I`, and `S_L=L^2I`, one has `Z_can=I` while
`Omega_L=(1/(2L))I -> 0`.  V005 reports canonical-amplitude PASS, ancestry
FAIL, and overall G2 FAIL.  No component-label overlap remains.

## No-regression results

### Representation and scalar rejection -- PASS as design rules

The target requires a covariant unitary massless Poincare pole bundle with
little-group helicities `+2,-2`; a two-scalar doublet fails.  On the two-state
`{+2,-2}` fiber, the ISO(2) commutators permit no nonzero null-translation
matrix elements.  Equal slopes alone remain insufficient.  V005 supplies no
physical Poincare/helicity data.

### TT zero mode and ground query -- PASS

`Pi_TT(0)` remains forbidden and homogeneous zero-character data are
separately typed.  The exact nonzero-k projector is idempotent with rank two.
Every ground query uses rank one or the complete basis-invariant density
`rho0=P0/rank(P0)`, never a selected vector in a degenerate space.

### Affine binding and momentum registry -- PASS

The primitive-cell source is correctly called affine until one independently
sealed FD certificate binds a common physical scale.  Conjugate source and
operator conversions cancel exactly.  Independent `L=5,10` enumeration
reproduces every character, equal conjugate norms, the `q=3,4,8` rays,
`q_min=3`, the doubling injection, and new odd fine-cover characters.

### Ancestry and factorization -- PASS as design rules

The generalized ancestry spectrum is invariant under invertible source-basis
changes.  Exact branchwise pole amputation recovers the source-independent
vertex.  Factorization readiness remains explicitly distinct from G3
longitudinal decoupling and soft universality.

## Custody and prior rejections

All eleven frozen V005 target files, 25 dependencies, nine author-manifest
entries, the author seal, and the `104/104` author replay pass from their
correct roots.  All nine frozen V004 rejection artifacts are pinned; their
recursive custody preserves the V003, V002, and V001 rejection histories.
Nothing in this pass erases or rewrites those audits.

## Ceiling and disposition

The sealed contract still lacks the native all-`G_L` source ledger, matched
pole data, an earned canonical source map, the Poincare/helicity payload,
physical FD binding, stability data, and multipoint factorization data.  It
therefore establishes no numerical pole, positive G2, G3 soft theorem,
gravity result, or `G`.  GE was neither audited nor repinned.

`GF_V005_REAUDIT_PASS__NARROW_DESIGN_CONTRACT_SEALED__V001_TO_V004_REJECTIONS_PRESERVED__AMPLITUDE_CLASSIFIER_DISJOINT_AND_TOTAL__LOWER_AND_UPPER_LIMITS_EXHAUSTIVE__NO_D_RAW_SCALING_INDEPENDENT__V004_DIAG_L_1_REGRESSION_FAILS_AS_REQUIRED__AMPLITUDE_ANCESTRY_SEPARATE__HELICITY_SCALAR_K0_GROUND_AFFINE_REGISTRY_FACTORIZATION_GATES_PASS_AS_DESIGN_RULES__PHYSICS_UNEXECUTED__NO_POSITIVE_G2_NO_G3_NO_GRAVITY_NO_G__NO_GE_REPIN`

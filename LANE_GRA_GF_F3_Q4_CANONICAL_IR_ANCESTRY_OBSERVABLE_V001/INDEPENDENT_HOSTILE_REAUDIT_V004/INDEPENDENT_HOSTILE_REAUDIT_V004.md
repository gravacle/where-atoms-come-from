# Independent hostile re-audit of repaired GF V004

**Target:** `LANE_GRA_GF_F3_Q4_CANONICAL_IR_ANCESTRY_OBSERVABLE_V001`  
**Target revision:** `V004`, root-repaired after the preserved V001--V003 rejections  
**Verdict:** **REJECT -- REPAIR_REQUIRED**  
**Independent replay:** `208/208` audit checks passed  
**Author replay:** `102/102` reproduced

## Executive finding

V004 correctly repairs the V003 no-map hole.  On the explicitly scored
finite-upper-bound domain, all four requested cases are now disjoint and
total: raw rank loss is FAIL; full raw rank without an applicable `D` is
INDETERMINATE for every raw scaling; an applicable `D` with positive bounded
`Z_Q^can` is PASS; and canonical rank or lower-bound loss is FAIL.  Normalized
ancestry remains independent, and the V002 overlap regression remains
amplitude PASS, ancestry FAIL, overall G2 FAIL.

One adjacent load-bearing defect prevents a seal.  The machine contract makes
a finite `sigma_max(Z_Q^can)` necessary for PASS, but neither the theorem's
FAIL branch nor the executable classifier handles an unbounded upper singular
value.  The prose then has no branch, while the executable returns PASS.  This
is a positive-G2 false-positive path, not merely an unreported diagnostic.

## Exact counterexample

Let the full-rank raw numerator and a prospectively derived, applicable source
map be

\[
 \mathcal N^{\rm raw}_{C,L}=I_2,
 \qquad D_L=\begin{pmatrix}L&0\\0&1\end{pmatrix}.             \tag{A1}
\]

Then

\[
 \mathcal N^{\rm can}_{C,L}=D_LD_L^\dagger
       =\operatorname{diag}(L^2,1),
 \qquad Z^{\rm can}_{Q,L}=\operatorname{diag}(L,1).          \tag{A2}
\]

Thus `rank N_raw = rank N_can = 2`, `sigma_min=1`, but
`sigma_max=L -> infinity`.  Under the frozen V004 predicates:

- PASS is false because the contract explicitly requires finite `sigma_max`;
- FAIL is false because there is no raw/canonical rank loss and no vanishing
  canonical lower singular value; and
- INDETERMINATE is false because an applicable frozen `D` exists.

The declared branch vector is `(false,false,false)`.  The exact author
function has no upper-bound input and returns PASS from

```text
canonical_outcome(raw_rank=2, controlled_power=False,
                  map_frozen=True, lower_bound=1).
```

That contradicts the machine contract and would permit an unbounded source
vertex through a gate that expressly requires a bounded canonical amplitude.

## Minimal repair

Add the upper-bound result to the executable decision input and complete the
applicable-map FAIL branch:

1. PASS requires full raw/canonical rank, positive limiting lower singular
   value, and a uniform finite upper singular-value bound;
2. FAIL includes raw/canonical rank loss, lower-bound loss, **or upper-bound
   divergence** after an applicable `D`; and
3. full raw rank without an applicable `D` remains INDETERMINATE regardless of
   its raw scaling.

Add (A1)--(A2) as an exact regression and retain both earlier repaired
regressions.  This is the smallest change; it adds no record or gravity
machinery.

## Rechecked gates

### Requested four-case partition -- PASS on its bounded domain

All 16 combinations of raw rank `{1,2}`, map availability, canonical rank
`{1,2}`, and lower-bound status were exhaustively enumerated at finite upper
bound.  Exactly one branch fires for each.  Both controlled-power and finite
nonvanishing full-rank no-`D` inputs are uniquely INDETERMINATE.

### Ancestry separation -- PASS

For `Delta_L=1/L`, `R_L=(L/2)I`, and `S_L=L^2I`, V004 retains
`N_can=I`, `Z_can=I`, and `Omega_L=(1/(2L))I -> 0`: canonical amplitude PASS,
ancestry FAIL, overall G2 FAIL.

### Representation, zero mode, and ground query -- PASS as design rules

The target retains the unitary massless Poincare pole-bundle requirement,
little-group helicities `+2,-2`, scalar-doublet rejection, nonzero-character
TT domain, and basis-invariant complete-ground density.  Independent exact
tests again give a rank-two idempotent nonzero-k TT projector and reject
`Pi_TT(0)`.  No physical pole/helicity data are supplied.

### Affine binding, registry, ancestry, and factorization -- PASS

Physical naming still requires one independent FD binding.  Source/operator
unit conversions cancel exactly.  `L=5,10` registry enumeration reproduces all
characters, conjugates, frozen rays, minimum shell, and cover behavior.  The
generalized ancestry spectrum remains basis invariant and exact branchwise
pole amputation recovers the source-independent vertex.  Factorization remains
distinct from G3 longitudinal decoupling.

### Custody and ceiling -- PASS

All eleven frozen V004 target files, 25 dependencies, nine author-manifest
entries, the author seal, and author `102/102` replay pass.  All nine frozen
V003 rejection artifacts are pinned; their own custody preserves V002 and
V001.  The target still claims design only: no matched-family pole evidence,
Poincare/helicity result, positive G2, G3, gravity, or `G`.

## Disposition

`GF_V004_REAUDIT_REJECTED__V001_TO_V003_REJECTIONS_PRESERVED__V003_NO_D_TOTALITY_REPAIRED__BOUNDED_FOUR_CASE_PARTITION_PASS__AMPLITUDE_ANCESTRY_SEPARATION_PASS__FINITE_SIGMA_MAX_REQUIRED_BUT_NOT_IMPLEMENTED__UNBOUNDED_CANONICAL_VERTEX_COUNTEREXAMPLE_EXACT__AUTHOR_FALSE_PASS_MATERIAL__UPPER_DIVERGENCE_TO_FAIL_REPAIR_REQUIRED__NO_GE_REPIN`

GE must not be repinned to this rejected custody.

# Independent hostile re-audit of repaired GF V003

**Target:** `LANE_GRA_GF_F3_Q4_CANONICAL_IR_ANCESTRY_OBSERVABLE_V001`  
**Target revision:** `V003`, root-repaired after the preserved V001 and V002 rejections  
**Verdict:** **REJECT -- REPAIR_REQUIRED**  
**Independent replay:** `190/190` audit checks passed  
**Author replay:** `101/101` reproduced

## Executive finding

V003 correctly repairs the V002 overlap.  Canonical amplitude and normalized
ancestry are now independent outputs: the exact V002 family is amplitude
PASS, ancestry FAIL, and therefore overall G2 FAIL.  The amplitude branches
are mutually disjoint on every branch that V003 defines.

The advertised three-valued canonical-amplitude classifier is nevertheless
not total.  A full-rank raw numerator that does not vanish, when no applicable
same-parent canonical map `D` has been derived, satisfies none of the frozen
PASS, FAIL, or INDETERMINATE predicates.  The author verifier silently labels
that input FAIL through a fallback branch, contradicting the theorem and
machine contract.  This is a material prospective-decision defect, not a
positive-G2 loophole.

## Exact counterexample

On every nonzero member of any declared infrared sequence, take

\[
  \mathcal N^{\rm raw}_{C,L}(k_L)=I_2,                 \tag{A1}
\]

and suppose no applicable same-parent source map `D_L(k_L)` has been derived.
This is an exact, finite, nonvanishing, rank-two raw source.  V003 itself says
that `D=I` is only the raw reporting convention and is not proof that the raw
source is a canonical nonderivative field.  Therefore:

- PASS is false because no applicable frozen `D` exists;
- FAIL is false because the raw source has no rank loss and no applicable `D`
  produces canonical rank or lower-bound loss; and
- INDETERMINATE is false under the frozen words because (A1) does not
  "vanish with a controlled momentum power."

Thus the declared branch vector is `(false,false,false)`.  The exact author
function

```text
canonical_outcome(raw_rank=2, controlled_power=False,
                  map_frozen=False, lower_bound=1)
```

instead returns `FAIL` through its final fallback.  That value is not licensed
by the theorem or contract.  It also confuses an unavailable canonical source
normalization with a physical amplitude failure.

## Why the outer G2 conjunction does not cure it

The outer G2 rule remains conservative: without canonical-amplitude PASS the
packet cannot establish positive G2.  But an outer conjunction cannot make an
undefined component classifier total, nor reconcile the executable FAIL with
the prose's absence of a matching branch.  GF is a prospective observable
contract, so a submitted dataset must receive one consistent amplitude label
before that label is joined to ancestry, helicity, cone, and stability gates.

## Minimal repair

Make the no-map rule complete:

> Every full-rank raw numerator for which no applicable prospectively frozen
> same-parent `D` exists is canonical-amplitude INDETERMINATE.  Report its raw
> scaling without post-result reweighting.

Then the partition is total and disjoint:

1. persistent raw rank loss is FAIL;
2. full raw rank with no applicable `D` is INDETERMINATE;
3. full raw rank with applicable `D` and finite positive canonical singular-
   value bounds is PASS; and
4. full raw rank with applicable `D` but canonical rank/lower-bound loss is
   FAIL.

Normalized ancestry remains a separate output and the overall G2 rule remains
conjunctive.  Add (A1) as the exact regression alongside the already-correct
V002 overlap regression.

## Rechecked gates

### V002 overlap separation -- PASS

For `Delta_L=1/L`, `R_L=(L/2)I`, and `S_L=L^2 I`, one has
`N_can=I`, `Z_can=I`, and `Omega_L=(1/(2L))I -> 0`.  V003 now reports exactly
amplitude PASS, ancestry FAIL, overall G2 FAIL.  No amplitude-label overlap
remains.

### Wigner/Poincare and scalar rejection -- PASS as a design rule

The target retains the covariant unitary massless Poincare pole-bundle gate,
little-group weights `+2,-2`, and explicit scalar-doublet rejection.  The
ISO(2) commutators permit no nonzero null-translation matrix elements within
a two-state `{+2,-2}` fiber.  No physical Poincare/helicity data are supplied.

### TT zero mode and degenerate ground -- PASS

`Pi_TT(0)` remains forbidden; homogeneous zero-character data are separately
typed.  The nonzero projector is exactly idempotent with rank two.  Residues
and ancestry retain the basis-invariant complete-ground density
`rho0=P0/rank(P0)` rather than a selected vector.

### Affine binding and momentum registry -- PASS

The source remains affine until one independently sealed FD certificate binds
the common physical cell and length.  Operator/source conversions cancel
exactly.  Independent `L=5,10` enumeration reproduces every character, equal
conjugate norms, the `q=3,4,8` rays, `q_min=3`, the doubling injection, and new
odd fine-cover classes.

### Ancestry, factorization, and ceiling -- PASS as design rules

The generalized ancestry spectrum is basis invariant.  Exact branchwise pole
amputation recovers the source-independent vertex.  Factorization readiness
remains separate from G3 longitudinal decoupling.  V003 still claims design
only: no matched-family pole data, Poincare/helicity result, positive G2, G3,
gravity, or `G`.

## Custody and disposition

All eleven root-frozen V003 target files, 25 dependencies, nine author-
manifest entries, the author seal, and the `101/101` author replay pass from
their correct roots.  `PRIOR_REJECTIONS_CUSTODY.sha256` preserves all nine
artifacts in the sealed V002 rejection packet, including its embedded V001
rejection record.  Prior rejection bytes are preserved rather than replayed
against the later in-place author revision.

`GF_V003_REAUDIT_REJECTED__V001_AND_V002_REJECTIONS_PRESERVED__V002_OVERLAP_REPAIRED__AMPLITUDE_ANCESTRY_SEPARATION_PASS__AMPLITUDE_CLASSIFIER_DISJOINT_BUT_NONTOTAL__FULL_RANK_NONVANISHING_NO_D_COUNTEREXAMPLE_EXACT__AUTHOR_EXECUTABLE_PROSE_MISMATCH_MATERIAL__NO_D_FULL_RANK_TO_INDETERMINATE_REPAIR_REQUIRED__NO_GE_REPIN`

GE must not be repinned to this rejected custody.

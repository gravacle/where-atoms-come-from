# Independent hostile re-audit of repaired GF V002

**Target:** `LANE_GRA_GF_F3_Q4_CANONICAL_IR_ANCESTRY_OBSERVABLE_V001`  
**Target revision:** `V002`, repaired after the preserved V001 rejection  
**Verdict:** **REJECT -- REPAIR_REQUIRED**  
**Independent replay:** `153/153` audit checks passed  
**Author replay:** `94/94` reproduced

## Executive finding

V002 materially repairs all three defects that caused the V001 rejection:

1. it requires a covariant unitary massless Poincare pole bundle with
   little-group helicities `+2,-2` and rejects a two-scalar doublet;
2. it defines TT projection only for nonzero characters and types `k=0`
   homogeneous data separately; and
3. it replaces a selected ground vector by a rank-one projector or the
   normalized complete-ground density `rho0=P0/rank(P0)`.

Its affine/physical binding, momentum registry, normalized ancestry,
factorization-readiness boundary, custody, and strict design-only ceiling also
survive independent checks.

One new load-bearing decision-rule defect remains. The advertised
PASS/FAIL/INDETERMINATE **canonical-amplitude** classifier is not disjoint:

- PASS is triggered by a frozen applicable canonical map with full-rank,
  finite, strictly positive `Z_Q^can`; but
- FAIL is also triggered when the separate normalized-ancestry gate fails.

Those conditions can hold simultaneously.

## Exact counterexample

On a two-dimensional source/pole image, take the infrared family

\[
 \Delta_L={1\over L},\qquad
 R_L={L\over2}I_2,
 \qquad
 \mathcal N^{\rm can}_L=2\Delta_LR_L=I_2,
 \qquad Z^{\rm can}_{Q,L}=I_2.                    \tag{A1}
\]

Use a prospectively frozen applicable map `D=I`, and let the full source Gram
be

\[
 S_L=L^2I_2,
 \qquad
 \Omega_L=S_L^{-1/2}R_LS_L^{-1/2}
 ={1\over2L}I_2\longrightarrow0.                 \tag{A2}
\]

This family respects `0 <= R_L <= S_L` for the infrared sequence. Equation
(A1) has full raw/canonical rank, finite `sigma_max`, and
`liminf sigma_min(Z_Q^can)=1`, so the declared amplitude PASS condition is
true. Equation (A2) fails normalized ancestry, so the declared amplitude FAIL
condition is also true. It is not INDETERMINATE.

The outer positive-G2 conjunction still rejects this family because ancestry
is separately required. That prevents a false positive G2 promotion, but it
does not make the advertised three-valued amplitude output single-valued.
The ambiguity is especially material in a prospective decision contract: the
same submitted dataset can be reported as both canonical-amplitude PASS and
canonical-amplitude FAIL. `RESULT.json` compounds the inconsistency by
describing residue FAIL without the theorem/contract ancestry branch.

## Minimal repair

Remove `failed normalized ancestry` from the canonical-amplitude FAIL branch.
Keep ancestry as its already-defined independent ancestry/G2 gate. Then add
the exact regression (A1)--(A2) and require:

- canonical amplitude = PASS;
- normalized ancestry = FAIL; and
- overall positive G2 = FAIL.

This one separation makes the amplitude PASS/FAIL/INDETERMINATE classifier
disjoint without weakening the final G2 conjunction or changing any physics.
An alternative ordered-precedence rule would change the meaning of
"canonical amplitude" and must be stated explicitly across theorem,
contract, result, and verifier; it is not present in the frozen V002 bytes.

## Independent gate results

### Wigner/Poincare and two-scalar gate -- PASS as a design rule

The target requires a unitary massless Poincare representation, covariant pole
projectors/residues/vertices, and little-group rotation weights `+2,-2`.
On a two-state fiber, the ISO(2) relations
`[J,T_plus]=T_plus`, `[J,T_minus]=-T_minus` admit no nonzero translation-
generator matrix element between weights `+2` and `-2`; the null translations
therefore act trivially. A scalar doublet with weights `(0,0)` fails. This is
a correct Wigner-helicity gate, not merely a TT-rank test. No physical
Poincare/helicity data are supplied by V002.

### TT zero mode -- PASS

`Pi_TT(0)` is forbidden, `k=0` retains only separately typed unprojected
homogeneous data, and soft TT limits use punctured directional sequences. An
independent rational nonzero-k projector is idempotent with rank two; the same
construction rejects `k=0` by exact division-by-zero.

### Degenerate ground -- PASS

The normalized complete-ground query is invariant under an exact rational
ground-basis rotation. A selected ground vector changes its query and is
correctly forbidden.

### Affine versus physical binding -- PASS

The source conversion is correctly affine until an independently audited FD
cell/scale certificate binds `a_*`, `v3`, and the cell map physically for all
sectors. The conjugate source conversion exactly cancels the operator
conversion. The required physical certificate remains absent and is not
claimed.

### Momentum registry -- PASS

Independent enumeration at `L=5` and `L=10` reproduces all `L^3` characters,
equal conjugate norms, the three `q=3,4,8` rays, `q_min=3`, the doubling cover
map, and new odd fine-cover characters.

### Ancestry and factorization -- PASS as design rules

The generalized ancestry spectrum is invariant under an invertible source-
basis change. Exact branchwise `D_C (Z_Q^can)^+` amputation recovers the
source-independent vertex. Factorization readiness remains explicitly
separate from G3 longitudinal decoupling.

### Custody and ceiling -- PASS

All 11 frozen V002 target files, 25 dependencies, 9 author-manifest entries,
the author seal, and the `94/94` author replay pass from their correct roots.
The original V001 rejection is retained in `PRIOR_V001_REJECTION.md`. V002
still claims design only: no matched-family pole data, Poincare/helicity
result, positive G2, G3, gravity, or `G`.

## Disposition

`GF_V002_REAUDIT_REJECTED__V001_REJECTION_PRESERVED__WIGNER_HELICITY_GATE_PASS__SCALAR_DOUBLET_REJECTED__TT_K0_EXCLUSION_PASS__DEGENERATE_GROUND_PASS__AFFINE_PHYSICAL_BINDING_PASS__MOMENTUM_REGISTRY_PASS__ANCESTRY_AND_FACTORIZATION_DESIGN_PASS__CANONICAL_AMPLITUDE_CLASSIFIER_OVERLAP_MATERIAL__MINIMAL_SEPARATION_REPAIR_REQUIRED__NO_GE_REPIN`

GE must not be repinned to this rejected custody.

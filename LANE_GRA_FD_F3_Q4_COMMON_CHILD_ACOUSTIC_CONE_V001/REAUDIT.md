# Independent hostile re-audit -- F3/q4 common-child acoustic cone

**Lane:** `GRA-FD-F3-Q4-CCMAC-V001`  
**Re-audit date:** 2026-08-27  
**Corrected theorem SHA-256:**
`60d012766675c12e82dd1731e202a6c0ed48f24e2697f589b63eecc3cb650287`  
**Corrected self-audit SHA-256:**
`d9d91f84bdf2af3cbbe92c648c8066a29ac302255f48b1a938e9ece91d83c3ed`  
**Corrected verifier SHA-256:**
`16d72697c3475110aae714c453c74fb45e043f326ed5006dd9c7064f50631bb8`

**Disposition:**
`PASS__EXACT_FINITE_Q4_COMMON_CHILD_INCIDENCE__EXACT_PROSPECTIVE_DETUNED_F3_FORM_SCHUR_SIBLING_KERNEL__EXACT_A3_BULK_SYMBOL_COVOLUME_AND_REFINING_AFFINE_ATLAS__EXACT_SUPPLIED_MASSLESS_ACOUSTIC_ACTION_TO_SMOOTH_IR_PRINCIPAL_CONE__RAW_APPEND_CONE_REMAINS_POLYHEDRAL__ALL_PHYSICAL_LIFTS_PHASE_COMMON_PROBES_CONSTRAINTS_STRESS_RGRL_B_AND_GRAVITY_REMAIN_OPEN`

## 1. Re-audit result

The corrected verifier passes `40/40`.  I independently rechecked the five
failure classes in the first audit.  All five are materially repaired, and no
new algebraic, custody, refinement, sector, or physical-binding defect was
introduced.

The original `AUDIT.md` remains the historical verdict on the earlier theorem
hash.  This `REAUDIT.md` supersedes that verdict only for the corrected hashes
listed above.

## 2. Closure of the five required corrections

### 2.1 Detuning and current-parent custody -- closed

`DETUNED-Q4-CARRIER-LIFT` now separately supplies:

- coexisting `S_N` and `S_(N+1)` carrier modes;
- the q4 append support;
- the inherited F3-form scalar transfer amplitude `t`;
- the new positive child/parent detuning `Delta`; and
- preparation, maintenance, source, work, controller, support, boundary,
  failure, and quarantine ownership.

The theorem explicitly says that current F3 supplies only the scalar transfer
form, not the carrier realization or detuning.  Section 3 and the final
summary consistently call the construction prospective and conditional on
the complete lift.  The prior source-custody defect is closed.

### 2.2 Exact low block and perturbative remainder -- closed

With `K_N=B_N^dagger B_N`, the corrected theorem defines the exact negative
spectral branch by

\[
 f_{\Delta,t}(K_N)
 ={\Delta I-\sqrt{\Delta^2I+4t^2K_N}\over2}.
\]

This is the correct parent-space spectral representation of the low branch of
the Hermitian bipartite block.  The frozen domain is explicit and the
remainder is stated in operator norm.  For each scalar eigenvalue `lambda` of
`K_N`, put `y=t^2 lambda/Delta^2`.  Taylor's theorem applied to

\[
 g(y)={1-\sqrt{1+4y}\over2},
 \qquad g'''(y)=-12(1+4y)^{-5/2},
\]

gives

\[
 |g(y)+y-y^2|\le2y^3.
\]

Functional calculus therefore gives exactly the displayed bound

\[
 \|R_N\|\le {2|t|^6\over\Delta^5}\|K_N\|^3.
\]

The Schur sign, denominator, and operator order remain correct.

### 2.3 One-common-child comparator -- closed

The three-mode block is now explicitly a separately truncated comparator with
incidence row `B=(1,1)`.  Its low bright energy is compared only with the
per-common-child matrix `-(t^2/Delta)[[1,1],[1,1]]`.  The text explicitly
denies that this is the full two-parent principal block of `K_N`.  The prior
finite-block mismatch is closed.

### 2.4 Massless action and refining physical speed -- closed

`MASSLESS-A3-COLLECTIVE-ACTION` now explicitly supplies all load-bearing
ingredients:

- a real collective field;
- exact constant-shift symmetry;
- absence of a mass term;
- the complete nearest-root quadratic action; and
- positive `chi,kappa`.

The theorem states that the Schur calculation earns only the `A3` support
algebra and does not derive the field, symmetry, coefficients, or phase from
current F3.  The acoustic dispersion and smooth infrared principal cone then
follow exactly from the supplied action.

The continuum composition now also freezes the required finite-speed scaling

\[
 {\kappa(a_*)\over\chi(a_*)}
 ={3c_{\rm phys}^2\over16a_*^2},
\]

and correctly changes the remainder to
`O(c_phys^2 a_*^2 |k|^4)`.  It expressly states that fixed `kappa/chi` would
give `c_* -> 0`.  The mathematical affine refinement is no longer being used
to imply a nondegenerate physical-time continuum by itself.

### 2.5 Pair-field type join and sector criterion -- closed

The theorem now separates PMMDC's finite statistical tangent from a
propagating q4/F3 field and introduces `Q4-PAIR-FIELD-LIFT` as an explicit,
unproved same-parent interface.  That interface owns the `S4` type join,
compact support, gluing, propagation, constitutive action, retained-sector
choice, and all physical ports.

For edge order `(12,13,14,23,24,34)`, the displayed opposite-edge involution
and projectors satisfy

\[
 P_A+P_E+P_T=I_6,
 \qquad (\operatorname{rank}P_A,
         \operatorname{rank}P_E,
         \operatorname{rank}P_T)=(1,2,3),
\]

and are mutually orthogonal idempotents.  The conventional irrep label `E`
is restored.  For every retained sector the theorem requires
`chi_X>0,kappa_X>0`, and (FD32) now equates speeds only over the retained set
`R`.  Projected-out sectors correctly impose no physical speed condition.

## 3. Independent replay of the unchanged mathematical core

The following results remain exact.

1. `B_N^dagger B_N=4I+A_N` globally on every finite q4 front.
2. The parent Schur equation has the displayed negative sign and
   `t^2/(Delta-z)` coefficient.
3. The tetrahedral roots have second moment `(16/3)I_V`.
4. The infinite-interior symbols obey `|b(k)|^2=4+A(k)` and the scalar band is
   even and quadratic at the symmetric `k=0` point.
5. The primitive `A3` covolume is `16a_*^3/(3sqrt(3))`; the tetrahedron volume,
   `N^3/6` bulk ratio, and exact stars-and-bars boundary terms are correct.
6. Under `a_* -> 0`, `N -> infinity`, and `Na_* -> L`, compact interior sets
   recede from the boundary in lattice units as `min_a m_a -> infinity`.
7. The supplied scalar wave action has
   `omega^2=(2kappa/chi) sum[1-cos]` and the infrared principal polynomial
   `-omega^2+c_*^2|k|^2`.
8. This collective characteristic cone is a different mathematical object
   from the four-ray raw append cone; the former does not change or contradict
   the latter.

## 4. Physical claim ceiling

The packet now passes because it proves a narrow conditional theorem, not
because the open physics has been converted into notation.  It still does
**not** prove:

- that current F3 or nature instantiates either lift;
- that current F3 generates the detuning or massless collective action;
- a common physical clock, probe cone, or absolute cell/volume binding;
- a protected helicity-two pole or constraint/Ward reduction;
- a lineage-conditioned retarded metric kernel or universal stress vertex;
- RGRL-B, Einstein dynamics, numerical `G`, or gravity.

Those ceilings are explicit in the theorem, self-audit, and verifier phrases.
The passing `40/40` executable replay checks the finite/algebraic claims; it is
not empirical or physical promotion of the open antecedents.

## 5. Final verdict

**PASS.**  No material defect remains in the corrected packet's declared
claim.  Its earned result is an exact finite common-child/Schur construction,
an exact `A3` bulk and refinement calculation, and an exact implication from
a prospectively supplied massless collective action to a smooth infrared
Lorentz principal cone.  The same-parent derivation of the physical lifts,
collective phase, tensor constraints, stress response, and gravity remains
the next scientific problem.


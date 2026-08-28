# Independent hostile audit — FV Coulomb-DPAR projected source rank

**Audit date:** 2026-08-27  
**Audited lane:** `LANE_GRA_FV_F3_Q4_COULOMB_DPAR_PROJECTED_SOURCE_RANK_V001`  
**Disposition:** `PASS_AFTER_FV_PURE_PREMISE_AND_BYTE_HYGIENE_REPAIR`  
**Independent replay:** `83/83` exact checks passed  
**Builder replay after repair:** `89/89` exact checks passed

## 1. Bottom line

The repaired FV theorem survives independent hostile audit.  Under FU
`S1`--`S9` **and** the newly explicit `S10 / FV-PURE` complete-source premise,
the ideal-Coulomb DPAR source has exact formal nonidentity projected-operator
rank six on the covering-matched q4 diamond-ice family.  The independent
calculation reproduced

\[
 \overline T_{C,n}(\lambda)=
 \frac{21}{8}(8-15\lambda)I-
 \frac{63}{8}(2-5\lambda)D_d,
\]

all `720` H6 histories, all `3600` proper-prefix denominators, the complete
`4 x 64` environment-independent Hermitian result, direct `E` rank two,
ring `A1+T2` rank four, and the actual projected-operator witness

\[
 \det W=-\frac{4678629417}{256}\ne0.
\]

This is an off-shell, formal effective-operator result.  It is not a
retarded/CTP response theorem, Ward identity, propagating tensor pole,
gravity theorem, `RGRL-B` derivation, or calculation of `G`.

## 2. Defects found before acceptance

### `M1` — material premise-inheritance defect, repaired

The first frozen FV core claimed its pure ideal-Coulomb result under FU
`S1`--`S9` alone.  That implication did not follow.  FU explicitly permits
residual local mutual kernels and ranked covariant remainders, requiring only
that the net `E` coefficient not cancel.  Such permitted terms can change the
H6 `A1/T2` coefficients, the cancellation slopes, and possibly the rank.

The repaired core now introduces `S10 / FV-PURE`.  It states that the complete
nonidentity first derivative before Feshbach is exactly the ideal-Coulomb
pair source plus the unchanged FS edge source; residual local terms must
vanish or coalesce into the same normalized radial law, and cross-node,
boundary, controller, and other remainder derivatives must be absent or
identity/reference terms.  The theorem and result now state explicitly that
`FV-PURE` is stronger than FU `S1`--`S9`.  This repairs the implication without
pretending to derive the physical solder.

### `N1` — source-byte hygiene defect, repaired

The first theorem contained a literal carriage-return byte inside the TeX
summation subscript for `\rm flippable`.  The repaired theorem contains no
forbidden ASCII control bytes.  The core was rebuilt and resealed after the
repair; this audit pins only the repaired hashes.

No material mathematical defect remains in the repaired claim.

## 3. Independent calculation

The audit implementation does not import or execute the builder verifier.
It independently:

1. rebuilt the normalized tetrahedral edge dyads `D_a` and the six normalized
   root dyads `R_ab` in `(xx,yy,zz,2xy,2xz,2yz)` coordinates;
2. rebuilt every permutation of the six H6 flips and recomputed every proper
   prefix gap from degree changes;
3. recomputed the pair-source derivative of each virtual energy relative to
   the appropriate initial ice endpoint;
4. differentiated the six flip numerators and all five resolvents in each
   history;
5. evaluated both endpoint-referenced directions and formed the explicitly
   frozen Hermitian forward/reverse average;
6. exhausted all four missing-label orientations and every one of the `64`
   local external-link assignments;
7. independently rebuilt the periodic `G_5` graph, its `500` elementary
   hexagons, and its absence of four-cycles;
8. used a separate integral-flow implementation to construct four actual
   global degree-two ice states and their switched ice endpoints; and
9. evaluated two global diagonal-difference functionals and four distinct
   global off-diagonal ring functionals before taking rank and determinant.

The independent path census is

| sorted proper-prefix gaps | multiplicity |
|---|---:|
| `(2,2,2,2,2)` | 96 |
| `(2,2,2,2,4)` | 144 |
| `(2,2,2,4,4)` | 216 |
| `(2,2,4,4,4)` | 192 |
| `(2,2,4,4,6)` | 72 |

and

\[
 \sum_{\pi\in S_6}\prod_{r=1}^{5}\delta(S_r)^{-1}=\frac{63}{8}.
\]

The endpoint-dependent residue of a one-way calculation is diagonal,
traceless `E`; it changes sign under endpoint reversal and cancels in the
Hermitian average for all `256` audited local environments.  No environment
average was used.

At `\lambda=-1/2`, the four actual ring matrix-element rows are

\[
\begin{array}{c|rrrrrr}
d&xx&yy&zz&2xy&2xz&2yz\\ \hline
0&231/8&231/8&231/8&-189/8&-189/8&-189/8\\
1&231/8&231/8&231/8& 189/8& 189/8&-189/8\\
2&231/8&231/8&231/8& 189/8&-189/8& 189/8\\
3&231/8&231/8&231/8&-189/8& 189/8& 189/8
\end{array}
\]

and span `A1+T2` with rank four.  The two normalized diagonal differences are

\[
 (-1,1,0,0,0,0),\qquad(-1,0,1,0,0,0),
\]

and span `E`.  These are operator evaluations, not a rank of six coefficient
vectors: the diagonal rows are expectation differences on three distinct
global direction-pair coverings, while the ring rows are four distinct
off-diagonal matrix elements between global ice states.  Every identity is
annihilated by all six functionals.

## 4. Hostile boundary tests

- **Coefficient rank versus operator rank:** passed using actual global
  projected matrix elements and a nonzero exact determinant.
- **Missing-label independence:** passed; `G_5` realizes `125` hexagons of
  each missing-label orientation, and four explicit completed ring entries
  are distinct.
- **Endpoint energy subtraction:** passed independently for both initial
  endpoints.  Reusing the forward endpoint would fail the Hermiticity test.
- **One-sided residue:** not promoted.  Its `E` contamination cancels only in
  the stated forward/reverse convention.
- **Identities:** rejected by diagonal differences and off-diagonal
  functionals; no scalar was counted as rank.
- **Lower folds:** `G_5` has girth six.  Orders two and four cannot connect
  different ice states, so lower kernels and their folds are diagonal and
  cannot enter the selected H6 ring entries.
- **Through H8:** an ice-to-ice odd edge support is Eulerian.  At girth six
  and support at most eight, the only nonempty connected possibilities are a
  hexagon or octagon; repeated links dress H6 and empty support is diagonal.
  Thus the six-witness determinant has
  `det W(h)=(-4678629417/256)h^24+O(h^26)` after the declared common factors.
  This is a formal-order statement, not a finite-coupling convergence bound.
- **Special slopes:** independently reproduced.  `lambda=2/5` removes ring
  `T2`; `lambda=3/5` removes ring `A1`; `lambda=0` removes direct `E`; the
  Coulomb value `-1/2` avoids all three.
- **Sign and normalization:** four identical physical ring factors
  `-h^6/U_d^5` leave the normalized determinant sign unchanged; the two
  direct rows carry `U_d`.  The first possible H8 replacement changes
  `h^24` to `h^26`.
- **Physical solder:** not rederived.  `S1`--`S9` and `FV-PURE` remain explicit
  antecedents.  Failure of any antecedent requires a fresh complete-source
  calculation.
- **Effective-operator scope:** the endpoint-referenced Hermitian average is
  the frozen off-shell convention being audited.  No representation-
  independent spectral or causal response is inferred from its static rank.

## 5. Publication-safe statement

On the covering-matched q4 diamond-ice family, conditional on FU `S1`--`S9`
and the explicit `S10 / FV-PURE` complete-source identification, the
source-before-Feshbach ideal-Coulomb deformation produces six independent
nonidentity projected source operators as a formal expansion through H8:
direct diagonal `E` rank two plus Hermitian H6 ring `A1+T2` rank four.  The
result closes the projected-source rank prerequisite for this completed
branch; the next scientific step is a dynamical retarded/CTP spectral test,
not promotion of static rank to gravity.

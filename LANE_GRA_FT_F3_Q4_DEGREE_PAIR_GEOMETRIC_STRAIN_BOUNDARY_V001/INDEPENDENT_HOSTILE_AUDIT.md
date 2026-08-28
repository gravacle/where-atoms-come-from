# Independent hostile audit

**Lane:** `GRA-FT-F3-Q4-DPGSB-V001`

**Audit date:** 2026-08-27

**Disposition:**
`PASS_AFTER_MICROSCOPIC_VS_PROJECTED_RANK__DPAR_SUFFICIENCY__A1_REPLACEMENT__AND_COMPATIBLE_H6_WITNESS_REPAIRS__DPAR_UNADOPTED__CTP_WARD_TENSOR_GRAVITY_G_OPEN`

## 1. Hostile findings and repairs

The core same-source-off nonuniqueness theorem survives, but the submitted
packet required four material scope repairs.

1. **Microscopic rank was being conflated with projected rank.** The six
   normalized sibling-root dyads and six independent microscopic pair Pauli
   strings give an exact six-coordinate source before Feshbach reduction. On
   the ice fiber, however, the direct pair-operator image is only `A1+E`,
   rank three. The inherited edge/Feshbach image has only an upper rank bound
   in this chain. All exact rank-six language is now explicitly microscopic;
   neither projected nor CTP rank six is asserted.
2. **Absolute minimality was not proved.** `DPAR` is one compact sufficient
   constitutive law. After the coframe, affine map, source convention, and
   exact `S4` covariance are fixed, its linear content is one common nonzero
   slope rather than six fitted weights. That does not prove it is the unique
   or logically necessary physical closure. The words “single minimal
   premise” were removed. `DPAR` remains neither inherited nor adopted.
3. **The scalar accounting was ambiguous.** The DPAR root family carries an
   `A1` component as well as `E`. It therefore replaces the FS scalar
   deformation of the degree term; it is not added on top of it. The one-edge
   source remains unchanged. Adopting DPAR would prospectively revise the
   FQ17a/FS degree-source rule while preserving `H[0]`, not derive the root
   weights from the unchanged additive rule.
4. **The local commutator needed compatible global states.** A rank-two table
   of local ring parities alone does not prove a nonzero matrix element of the
   flippability-projected global H6 operator. The audit independently built a
   periodic degree-two ice state for both alternating orientations of every
   one of the twelve `G_5` hexagons through the reference vertex. The exact
   H6 sign and normalization then hold on actual finite matrix elements, and
   distinct ring symmetric differences prevent cross-ring cancellation.

The audit also made the DPAR derivative explicit. For real symmetric `j`,

\[
F(j)=I-{j\over2}+O(j^2),\qquad
{ |F(j)r_{ab}|^2\over |r_{ab}|^2}
=1-j:\widehat R_{ab}+O(j^2).
\]

Thus a real differentiable law with `g'(1)=lambda` gives

\[
{\partial H\over\partial j_{ij}}\bigg|_0
=-{U_d\lambda\over2}\sum_{a<b}\widehat R_{ab}^{ij}P_{ab},
\qquad
Q^{ij}=-2{\partial H\over\partial j_{ij}}\bigg|_0
=U_d\lambda\sum_{a<b}\widehat R_{ab}^{ij}P_{ab}.
\]

This independently fixes the sign, the factor two, and the normalized root
tensor. Real `g`, real `j`, and Hermitian `P_ab` preserve Hermiticity.

## 2. Independent exact reconstruction

The hostile executable imports no builder code and uses only Python standard
library exact arithmetic. It verifies:

- the frozen dependency bytes and appended-byte tamper rejection;
- rank four for the four normalized tetrahedral edge dyads, their exact
  diagonal-traceless two-dimensional `E` null, and rank four after all
  additive pair weights are included;
- rank six for the six normalized sibling-root dyads;
- all 24 tetrahedral label permutations as exact orthogonal actions, including
  covariance of every normalized root dyad;
- `(d-2)^2=I+(1/2) sum_(a<b) Z_a Z_b` on all sixteen local basis states and
  full-Hilbert independence of the six pair Pauli strings;
- ice pair rank three, centered ice pair rank two, root-pair type `A1+E`,
  exact `T2` annihilation, and scalar value `-2I`;
- identical source-off identity/pair coefficients for FS and DPAR, microscopic
  derivative ranks four and six, and direct projected root-pair rank three;
- every coordinate and pair in the DPAR differentiation, including off-
  diagonal contraction normalization, sign, and `Q=-2 dH/dj` factor;
- a nonzero quadratic-contact Hessian with zero source-off gradient, which
  prevents linear source rank from being promoted to full CTP rank;
- the 250-vertex, 500-edge periodic `G_5` graph, its twelve hexagons through
  one vertex, and all six local ring-edge types;
- 24 exact bipartite `b`-matching constructions extending the twelve rings
  and two alternating orientations to global degree-two ice states; and
- exact nonzero H6 matrix-element channels of rank two on `E`, with the
  `+2 J_6 U_d lambda` commutator factor and no cross-ring cancellation.

## 3. Result that survives

The source-off reduced parent plus its existing pair observables does not
determine its geometric strain derivative. The frozen FS query and the DPAR
root query are two Hermitian, source-before-Feshbach families with exactly the
same `H[0]` and different microscopic derivative ranks. That is a genuine
nonuniqueness theorem.

Conditional on separately deriving, adopting, or independently calibrating
DPAR with nonzero slope, the microscopic source has exact rank six and its
two new directions use already existing centered ice-pair `E` operators. On
the explicitly constructed compatible finite ice states, the inherited H6
ring dynamics acts nontrivially on both local `E` directions at formal leading
order. This is a conditional local source-and-commutator result, not an
inheritance theorem for DPAR.

## 4. Ceilings retained

The audit does not establish exact projected source rank six, a stationary
state, a uniform or `k=0` response, a nondegenerate state-dependent CTP
kernel, Ward constraints, a tensor pole, helicity two, RGRL-B, gravity, or
`G`. Generated H6/H8 source derivatives and prospective contacts must remain
under source-before-Feshbach custody. The nonzero block-local formal
order-six coefficient cannot be relabeled as a collective gravity mode.

The independent executable completes **646/646** checks. Its
full replay is frozen in `INDEPENDENT_HOSTILE_VERIFICATION.txt`.

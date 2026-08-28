# Independent hostile audit — q4 Clifford collective cone

**Lane:** `GRA-FC-Q4CCC-V001`  
**Date:** 2026-08-27  
**Audited:** `THEOREM.md` and
`verify_q4_clifford_collective_cone.py`

**Disposition:**
`FAIL_AS_A_CURRENT_SAME_PARENT_Q4_F3_RESULT__PASS_AS_AN_EXACT_PROSPECTIVE_CONSTANT_COEFFICIENT_A3_STENCIL_THEOREM_AFTER_CUSTODY_NARROWING`

## 1. Replay and independent algebra

The supplied verifier passes `32/32`.  Independent calculation confirms the
mathematical core:

- the tetrahedral frame has Gram entries `1,-1/3` and second moment
  `(4/3) I`;
- the six roots have second moment `(16/3) I`;
- in the displayed Cartesian frame, the six-dyad coefficient matrix has
  determinant `-524288/729`, hence exact rank six;
- the primitive-root Gram determinant is `256/27`, giving covolume
  `16/(3 sqrt(3))`;
- (FC07) has zero gradient and quadratic coefficient
  `(16/3)t a_*^2`;
- for (FC10), assuming the declared shift representation,
  `(T-T^dagger)` is anti-Hermitian, commutes with the internal Clifford
  matrix, and the prefactor makes the term Hermitian;
- `c_0=3/16` makes the principal map the identity;
- (FC14)--(FC16) and the `O(a_*^2 |k|^3)` expansion are correct; and
- the six weight variations span `Sym^2(V)` and give
  `delta(B^T B)=2 delta B` at `B=I`.

The central result is therefore real: **if** one supplies the constant-
coefficient Clifford root stencil on a coherent `A3` carrier lattice, its
infrared principal polynomial is the Weyl/Lorentz polynomial while the raw
q4 append cone remains polyhedral.

## 2. Material defect 1 — (FC07) is not the existing EO scalar SWAP

Lines 28--30, 98--127, 335, 368, and the final disposition attribute (FC07)
to the already displayed/current scalar sibling-SWAP dynamics.  That custody
is false.

EO's `FOCST-3a` supplies one **prospective local common-child** Hamiltonian on
`L,C,R,G`.  It does not supply the direct, translation-invariant root-graph
hopping whose symbol is (FC07).  EO explicitly says that current F3 has not
earned the q4 endpoint-to-carrier lift and that a globally glued signed walk
remains open (`FOCST` lines 398--411).  A global product or sum of overlapping
common-child swaps is not proved to equal (FC07), and a scheduled discrete-
time swap network need not have the same one-band dispersion.

Accordingly, the exact obstruction is:

> the prospectively supplied **even scalar nearest-root stencil (FC07)** has a
> quadratic band edge.

It is not an obstruction proved for the current EO/F3 dynamics.  This defect
requires claim narrowing; the algebra need not change.

## 3. Material defect 2 — the coherent carrier lift is new, not inherited

The q4 merger earns an abstract count front, sibling diamonds, and an `A3`
adjacency relation.  It does not earn a globally coexisting orthogonal
carrier-position factor with unitary/partial shifts `T_alpha`.  In the bounded
q4 witness, the reachable code obeys `m=m(w)`.  Changing `m` while acting as
the identity on its complete word register generally leaves that reachable
code.  EO avoids this by requiring a distinct probe-occupancy factor over
already realized coexisting endpoints and labels that requirement
prospective.

Thus lines 51--61, 131, 145--162, and 338--340 overstate what is inherited.
The theorem must add an explicit `COHERENT-A3-CARRIER-LIFT` premise: a separate
carrier-position Hilbert factor over coexisting endpoint modes, compatible
shift maps, one common structural `Z`, and complete route/controller/boundary
ports.  With that premise, (FC10) adds no new **combinatorial** graph, but it
does add a previously unearned coherent physical realization of that graph.

This also means the claim class cannot yet be “same-parent” in the sense of an
already derived F3 successor.  It is an exact conditional construction whose
same-parent derivation remains the next physics problem, as line 162 correctly
acknowledges.

## 4. Material defect 3 — the bounded witness does not by itself supply the
continuum sequence

The inherited q4 witness is finite at each fixed cap `R`.  Condition (FC19)
cannot hold in one such model.  The local limit additionally requires a
prospectively supplied family

`a_* -> 0`, `R(a_*) -> infinity`, `sum_a m_a(a_*) <= R(a_*)`,

with compatible interior carrier stencils and boundary completions.  No
uniform resource bound or physical scalable parent follows from the bounded
witness.

Moreover, lines 249--251 are literally false as written: a finite compact
neighborhood cannot contain the complete infinite lattice `a_* Lambda`.  The
correct local exhaustion statement is that, for every compact `K subset V`,
the points of the translated full lattice lying in `K` eventually coincide
with points available in the translated finite simplex patch.  Under that
correction and the explicit refining family, the local differential-operator
limit is sound.

## 5. Material ceiling — chirality and the inherited `S4` parent

The caveat at lines 139--143 is mathematically correct: a two-component Weyl
fiber realizes the oriented spin action, while the squared cone is invariant
under full `S4`.  It does not, however, make (FC10) itself a linearly unitary
full-`S4` operator.  The q4 parent uses complete-port `S4` covariance.  A
lawful symmetry-preserving successor must either supply the doubled Pin
carrier or explicitly declare and physically own a chirality/orientation
selection and the resulting reduction of microscopic symmetry.

The occupied two-content qutrit sector of F3 supplies a two-dimensional
vector space, but not this spin representation or its type join.  Its current
transfer preserves the record-content label and its parent has a declared
content-exchange covariance.  Clifford hopping generally mixes those
contents.  Reusing that sector therefore requires a proof that record-value
custody and content symmetry are preserved; otherwise a separate coin/spinor
factor is owed.  The theorem's admission that the vertex is not derived
prevents a gravity overclaim, but “already present ... carrier” is still too
strong.

## 6. Volume statement

(FC21) is the exact primitive `A3` lattice covolume.  (FC22) is the exact
coordinate four-parallelotope covolume after an independent time-length
calibration; under the lock it agrees with `16 a_*^4/(3 sqrt(3))`.  It is not
EO's top spatial tetrahedron volume, and it is `24` times EO's causal
four-simplex volume.  Those are different cells, not an algebraic conflict.

Calibration alone does not physically bind or assign ownership of this
abstract lattice cell.  Until the separate physical cell-binding/census gate
is supplied, the result should be called a **calibrated coordinate
covolume**, not a physically instantiated event-cell volume.

## 7. Lesser claim defects and verifier ceiling

- “Minimal” and “smallest” are not proved relative to a declared operator
  class.  The two-band Clifford construction is economical, but minimality
  should be removed or formally scoped.
- A finite truncated root stencil can already be made a Hermitian finite
  matrix by including each retained edge with its adjoint.  A terminal or
  quarantine completion is a physical complete-port obligation, not required
  merely for algebraic self-adjointness.
- The phrase “physical distance” before `a_*` is physically calibrated should
  read “relational chart distance.”
- The verifier checks the finite algebra and the presence of ceiling phrases.
  It does not test source custody, coherent support, the bounded-to-continuum
  family, `S4` implementation, record-content preservation, or physical cell
  binding.  Its `32/32` result must not be used as a full promotion gate.

## 8. Promotion decision

Do **not** promote the packet with its present “current scalar dynamics,”
“inherited carrier motion,” or “same-parent” wording.  After the narrow
corrections above, promote the following result:

`EXACT_EVEN_SCALAR_A3_STENCIL_QUADRATIC_NO_GO__EXACT_PROSPECTIVE_CLIFFORD_A3_STENCIL_WITH_WEYL_LORENTZ_IR_CONE__EXACT_ROOT_COVOLUME_AND_SIX_WEIGHT_COMETRIC_RANK__COHERENT_A3_CARRIER_LIFT_REFINING_PARENT_FAMILY_CONTENT_SPIN_JOIN_FULL_S4_VERTEX_AND_PHYSICAL_CELL_BINDING_OPEN__NO_CURRENT_F3_DERIVATION_AND_NO_GRAVITY`

# Independent hostile audit -- q4/F3 RGRL-B constraint-origin screen

**Lane:** `GRA-FP-F3-Q4-RCOS-V001`

**Audit date:** 2026-08-27

**Disposition:** `PASS_WITH_STRICT_FINITE_AND_ROTOR_COMPLETION_CEILINGS`

## 1. Independence and result

I did not import, call, or accept the builder verifier as evidence.  I rebuilt
the finite periodic diamond incidence matrices, the hard-core ice function
algebra, the `S4` characters, the closed-ring charge changes, the one-link
microscopic commutator, the rotor phase-space count, and the Maxwell/linearized
ADM principal symbols in a separate executable.

The central conclusion survives hostile review.  The currently owned q4/F3
ice branch contains one local scalar compact-`U(1)` Gauss species.  Its local
occupation-pair identities are reducible hard-core consequences of that
scalar law, not an independently inherited scalar-plus-vector first-class
packet.  Matching an `S4` module dimension or the final count of two physical
polarizations does not turn this spin-one architecture into RGRL-B.

No material theorem correction was required.  This audit makes one point
especially explicit: saying that a pair relation is the zero operator is a
statement about the projected algebra.  The decisive independence test is
that, in the hard-core parent algebra, the relations factor through the one
scalar ice generator; they therefore do not supply four new gauge directions.

## 2. Finite incidence and inherited Gauss algebra

For each independently generated periodic diamond quotient with linear size
`L=2,3,4`, the audit finds

\[
 |V|=2L^3,\qquad |E|=4L^3,\qquad
 \operatorname{rank}B=|V|-1.
\]

Breadth-first traversal separately confirms connectedness.  Every oriented
incidence column sums to zero, and exact rational elimination finds no second
row dependency.  Thus the sole linear dependency is the global one.  This is
the finite realization of one scalar equation per vertex, not four constraint
species per vertex.

For alternating cycles of lengths six and eight, the endpoint occupation
change lies in `ker B`.  A one-link change does not.  This rederives the
operator statements without using the theorem's ring calculation:

- the leading hexagon ring commutes with every inherited Gauss generator;
- the order-eight alternating octagon does likewise;
- a dressed hexagon has the same closed-cycle endpoint, while its diagonal
  coefficient is a function of commuting occupations, so it also preserves
  Gauss; and
- the unprojected `-h sum_e X_e` term fails the test because the exact
  two-state matrix commutator `[E_e,X_e]` is nonzero.

The exact symmetry therefore belongs to the projected effective branch.  At
finite `U_d`, the supplied microscopic parent admits virtual charge defects
and does not own an exact microscopic tensor gauge algebra.

## 3. Why four pair equations are not four gauge generators

Write `S=s_1+s_2+s_3+s_4` and `j_ab=s_as_b`, with `s_a^2=1`.  The audit checks
on all sixteen hard-core configurations, not just the six ice states, that

\[
 2(j_{12}-j_{34})
 =S(s_1+s_2-s_3-s_4),
\]

with the two analogous opposite-pair formulas, and

\[
 2\left(2+\sum_{a<b}j_{ab}\right)=S^2.
\]

Consequently all four affine pair relations vanish when `S=0`, but three are
the scalar constraint multiplied by state-dependent diagonal factors and the
fourth has no linear part on the constraint surface.  Their formal coefficient
vectors have rank four only when six pair symbols are incorrectly treated as
independent coordinates.  On the actual ice fiber the pair-image rank is
three (`A1+E`), its centered tangent rank is two (`E`), and all four relation
operators act as zero.

Pair observables commute with the inherited Gauss algebra because they are
functions of electric occupations.  Gauge invariance under an existing orbit
does not manufacture a new orbit.  To reinterpret the six pair symbols as a
new rank-two canonical field would require independent conjugates, a
symplectic form, four nontrivial generators, Hamiltonian preservation, and
closure.  None is present in this finite parent.  The audit therefore accepts
the theorem's distinction between algebraic identities and first-class
constraints.

## 4. `S4`, parity, and continuum spin

Independent character reconstruction gives

\[
 \mathbb R^{\Omega_2}=A_1\oplus E\oplus T_2,
 \qquad
 \operatorname{Sym}^2(T_2)=A_1\oplus E\oplus T_2.
\]

That equality is only an `S4` restriction.  Complement sends every one-link
function to its negative and leaves every pair function invariant.  The
available odd one-link `T2` is therefore the vector/Maxwell channel.  The
even local pair image contains `A1+E`, not the missing tensor `T2`.  In the
continuum little group the physical characters are `2 cos(theta)` for
helicity one and `2 cos(2 theta)` for helicity two; at `theta=pi/2` they are
`0` and `-2`.  Equal tetrahedral irrep labels do not establish equal `O(3)`
spin, parity, pole, or canonical type.

The audit also confirms the fixed-momentum rank trap.  One scalar amplitude
and its three spatial derivatives have rank one.  Independently specifiable
lapse-plus-shift data have rank four.  A derivative list cannot be counted as
three new gauge functions.

## 5. Rotor and linearized-ADM count

The standard compact-rotor completion gives

\[
 \dim\Gamma_{\rm phys}=2\bigl(|E|-(|V|-1)\bigr).
\]

On the periodic diamond quotients this is twice `2L^3+1` configurations.  The
extra zero-momentum harmonic configuration is a global finite-volume mode;
it is not another local generator.  This count is not asserted to be the
dimension of the finite six-state hard-core fiber.

At generic nonzero integer momentum the independently reconstructed symbols
have ranks

\[
 \operatorname{rank}(k_iE^i)=1,
 \qquad
 \operatorname{rank}(k_j\pi^{ij})=3,
 \qquad
 \operatorname{rank}
 (k_ik_jh^{ij}-|k|^2h^i{}_i)=1.
\]

The spatial symmetric-tensor gauge map has rank three, the scalar-curvature
row annihilates it, and transverse plus trace conditions leave a two-
dimensional TT space.  Thus Maxwell and linearized ADM both end with four
physical phase dimensions, but they start with different canonical fields
and quotient different first-class algebras.  Equal output count is not
equivalence.

## 6. Custody, negative tests, and ceilings

The independent replay validates all eleven frozen dependency digests and
requires every exact lane payload member in `MANIFEST.sha256`.  The outer
`SEAL.sha256` authenticates the manifest bytes.  Negative tests detect a
changed theorem byte, a changed dependency byte, a changed manifest byte,
and omission of a required payload member.

Run:

```bash
python3 LANE_GRA_FP_F3_Q4_RGRLB_CONSTRAINT_ORIGIN_SCREEN_V001/independent_hostile_audit.py
```

Expected result: `SUMMARY 114/114 independent hostile checks passed`.

The accepted result has four strict ceilings:

1. It is a finite/local algebra theorem, not a thermodynamic no-go.  A new
   collective redundancy at a critical limit remains logically open.
2. The Dirac count belongs to the stated rotor completion, not to a literal
   count of the finite hard-core ice Hilbert space.
3. The Maxwell continuum typing is imported from the separately pinned
   `MAXWELL-IR` result.
4. Neither equal physical polarization count nor local `S4` character proves
   a tensor pole, helicity two, a common cone, nonlinear constraint closure,
   RGRL-B, Einstein dynamics, or numerical `G`.

## 7. Sharp next lawful test

The next no-lab test is the fixed-parent connected TT four-point/spectral
calculation already identified by the order-eight lane.  Any isolated even
tensor candidate must pass two tests together:

1. nonzero thermodynamic TT residue on a common linear cone; and
2. four independently derived null directions of its effective action with
   rank-three vector plus independent scalar symbols and constraint-preserving
   closure.

A tensor-looking level without that null structure is not RGRL-B.  If the
fixed link/ring parent fails either test, the minimal lawful successor is a
genuinely collective loop/surface variable derived from the same parent,
with its conjugate and gauge algebra obtained from link commutators.  Merely
renaming a local pair or adding a rescue interaction is not warranted.

**Final audit result:**
`PASS_WITH_STRICT_FINITE_AND_ROTOR_COMPLETION_CEILINGS`.

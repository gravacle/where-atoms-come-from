# Independent Hostile Audit — GL6CP Orientation Accumulation

## Disposition: REPAIR_REQUIRED

The central representation-theory result survives independent
reconstruction, but the packet cannot be promoted in its frozen form.  Its
orientation-moment antecedent is not mathematically defined on the quotient
space it names, and the target is not packaged or sealed.

## Results independently recovered

Using the `S4` conjugacy classes rather than the target enumeration, the
standard three-dimensional character gives

\[
 \chi_W=(6,2,2,0,0),\qquad
 \chi_{\operatorname{Sym}^2W}=(21,5,5,0,1),
\]

and hence

\[
 \dim[\operatorname{Sym}^2W]^{S_4}=3,qquad
 \dim[\operatorname{Sym}^2W\otimes W]^{S_4}=9.
\]

An independent infinitesimal `SO(3)` calculation, in a different six-pair
coordinate order and at two modular primes, gives constraint ranks `19` and
`122`.  Four explicit quadratic invariants lie in the common kernel, so the
rotational dimensions are exactly `2` and `4`.  Because the rotational
subspace is contained in the cubic subspace, Haar projection restricted to
the nine-dimensional cubic space has rank four and removes five quadratic
anisotropies.

Restricting the four rotational scalars to
`t=(h_yz,h_zx,h_xy)` independently gives

\[
 0,\quad2|k|^2I,\quad0,\quad(|k|^2I-D)+O.
\]

Thus the projected `T2` block obeys `B+C=0`.  The conversion to
`c_eff+d_eff=kappa_eff/2` is valid only under the target's stated additional
condition that the accumulated response admits the GL6CO pullback
parameterization; the target does preserve that qualification.

The degree-six requirement is also correct: a quadratic response contains
four coframe factors from its two symmetric-tensor legs and two from its
momentum leg.  Componentwise moment error bounded by `epsilon` therefore
gives

\[
 \|\Delta C\|_{\max}\le\epsilon\|C\|_1
\]

with coefficient one.  Nonnegative, owner-correct occurrence weights and
typewise control are genuine antecedents, not conclusions of group theory.

Finally, for

\[
 A_{ij}=\tfrac12(E_{ii}-E_{jj}),\qquad
 B_{ij}=E_{ij}+E_{ji},
\]

the three commutators are the three elementary skew generators and have rank
three.  The ordered group-commutator formula is correct through second order.
Capacity is not dynamics.  This proves algebraic capacity only; the target
correctly leaves physical transport, inverse legs, retained ordering, and
orientation dynamics open.

## Material defect CP-A1 — quotient coordinates

The theorem declares `Q_C in SO(3)/Gamma` and makes `mu_(tau,L)` a measure on
that quotient.  CP10 then integrates raw coordinate monomials
`prod Q_ij`.  Those monomials are not functions on `SO(3)/Gamma`.

There is an explicit kill test.  The identity and any nontrivial proper cubic
element `gamma` represent the same identity coset.  Choosing a quarter-turn
with `gamma_00=0` gives

\[
 (I_{00})^6=1,\qquad(\gamma_{00})^6=0.
\]

Thus `epsilon_(tau,L)^(6)` changes with the representative even though the
contracted response (CP07) does not.  CP10 is formally undefined as written.

Two narrow repairs preserve the substantive theorem:

1. define a same-parent authenticated lift `Q_C in SO(3)`—including the
   ordering that fixes its representative—and state that the raw-moment
   criterion is sufficient but lift-dependent; or
2. retain quotient-valued coframes and replace every raw monomial by its
   right-`Gamma` average.  Those averaged monomials are quotient functions,
   Haar averages are unchanged, and the same `l1` bound follows because the
   intrinsic tensor is `Gamma` invariant.

The second form is intrinsic and is the recommended repair.

## Improper port actions are not a defect

The target's underlying `S4` port action contains determinant-minus-one
matrices.  Independent enumeration verifies that

\[
 R\longmapsto\widetilde R=\det(R)R
\]

is a faithful 24-element proper cubic group.  Since both the pair field and
the inversion-even momentum input carry rank two,
`rho_W(tilde R)=rho_W(R)`.  The target's proper-cubic replacement is therefore
valid, though the repaired theorem should state this construction explicitly.

## Custody defect CP-A2

The frozen target contains six files but no `verify_packet.py`,
`MANIFEST.sha256`, or `SEAL.sha256`.  Its README nevertheless instructs the
reader to run `verify_packet.py`.  The audit pins all six target bytes, but
the author packet itself must be packaged and sealed after CP-A1 is repaired.

## Scope assessment

The claim ceiling is otherwise sound.  The target does not treat a port
relabeling as a new orientation, does not average `T2` alone, does not infer
Haar weights from dense support, and does not claim that F3 already realizes
the required orientation distribution.  No phase, masslessness, Ricci,
gravity, or `G` conclusion follows at this gate.

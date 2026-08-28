# Independent hostile audit -- complete reduced q4 source-rank audit

**Lane:** `GRA-FS-F3-Q4-CSRAV-V001`

**Audit date:** 2026-08-27

**Disposition:**
`PASS_AFTER_PROSPECTIVE_WEIGHT_OPERATOR_RANK_AND_NONUNIFORM_SOURCE_REPAIRS__REDUCED_QUERY_MICRO_RANK4__EFFECTIVE_RANK_AT_MOST4__UNREDUCED_PHYSICAL_COMPLETION_OPEN`

## 1. Executive verdict and material repairs

The branch-relative obstruction survives hostile reconstruction.  On the one
frozen reduced CW/FM parent and the one prospectively frozen additive FQ17a
query, the microscopic source has exact operator rank four, `A1+T2`, and the
two diagonal-traceless `E` directions are exact nulls.  A fixed Feshbach map
cannot restore those missing first derivatives; the through-order-eight
effective image therefore has rank **at most** four.  This is enough to fail
the frozen six-off-shell-direction prerequisite without running the expensive
CTP spectrum.

Four repairs were required before that statement was auditable.

1. The original packet said FQ17a *fixed* the degree-square tensor.  FQ17a
   instead requires occurrence multiplicities and node weights to be frozen
   prospectively.  The corrected theorem now declares occurrence one on each
   of the degree square's four actual incident supports.  This gives
   `W_v=sum_a D_a=(4/3)I` as one lawful additive query; it is not a unique
   tensor derived from the source-free Hamiltonian.
2. Rank four of four tensor weights is only an upper bound on an
   operator-valued source.  The corrected proof now supplies the lower bound:
   occupation-basis matrix elements differing on only one link isolate the
   four distinct `X_e` operators.  The theorem also freezes `h != 0`, without
   which the claimed exact lower rank need not hold.
3. The displayed word formula

   \[
   M_w=\sum_aN_aD_a-(m-1)W_v
   \]

   is exact under one common uniform source.  For a block-resolved source,
   virtual-gap derivatives carry history-dependent energy fractions and are
   not generally that integer formula.  The corrected statement uses only
   what remains exact: every such derivative is a linear combination of the
   same `D_a` and scalar `W_v` tensors.
4. Omitted storage, carrier, formation, feedback, boundary-exchange,
   controller, and port terms are now described as definition-level
   exclusions of the selected reduced Hamiltonian.  They are not zero-weight
   assignments and are not claimed physically absent from unreduced BS04.

The corrected builder replay passes `75/75`; the independent hostile replay
passes `129/129`.

## 2. Periodic family, cover, and no-wrap proof

I rebuilt `G_5` and `G_10` independently from the four shifts
`{0,e1,e2,e3}`.  They have respectively `2L^3` vertices and `4L^3` edges,
are connected, bipartite, closed, and four-regular.  `G_5` has no simple
four-cycle, exactly `4*5^3` simple hexagons, and every hexagon uses three q4
labels twice.

Every `2p`-step same-part displacement is a sum of `p` forward shifts minus
`p` reverse shifts.  For `p<=4`, each coordinate has magnitude at most four,
so it cannot be a nonzero multiple of five.  Exhaustion of all label balances
through eight steps and direct lifting of every simple `G_5` six- and
eight-cycle both confirm zero winding.  Thus all cycles through order eight
are inherited local cycles, not quotient artifacts.  Reduction modulo five
has eight-point fibers from `G_10`, and every neighborhood maps bijectively;
it is the claimed eight-sheeted graph cover.  The same coordinate proof
extends to `L_r=5*2^r`.

This is a covering-matched mathematical family.  It remains neither a nested
Hilbert-space inclusion nor a physical record lineage, and periodic physical
identification is not derived from F3.

## 3. Complete term census versus prospective source choice

CW and FM explicitly select

\[
H=U_d\sum_v(d_v-2)^2-h\sum_eX_e,
\qquad E_R=0,quad U_d>0,quad h\ne0.
\]

Relative to this displayed reduced Hamiltonian, flips and degree squares are
the complete nonzero microscopic term list; the formal onsite class has zero
coefficient.  The broader BS04 dependency separately contains carrier,
formation, feedback, and port sectors.  Their absence from the selected
equation is legitimate branch reduction, but says nothing about their
physical absence.  Reattaching any such sector requires its own prospectively
frozen tensor and reopens the source-rank census.

The degree square is genuinely supported on its four incident links.  FQ17a
therefore permits the prospectively declared occurrence-one additive weight
`W_v=sum_aD_a`.  The audit does not promote that permission to uniqueness.
An independently derived nonadditive node tensor, root/cross-dyad query,
rotating coframe, storage tensor, or port tensor defines a different source
packet and is outside this theorem.

## 4. Exact microscopic operator rank

Using symmetric coordinates `(xx,yy,zz,xy,xz,yz)`, with the conventional
factor two in off-diagonal contractions, the four tetrahedral dyads have
exact rational rank four.  Their sum is scalar, and their simultaneous kernel
is exactly

\[
\mathcal N_E
=\{\operatorname{diag}(x,y,z):x+y+z=0\},
\qquad \dim\mathcal N_E=2.
\]

I then constructed the four one-link Pauli flips on a four-link occupation
block and the actual diagonal `(d-2)^2` operator.  The four flips are linearly
independent, and the six source-coordinate operators built from flips plus
the scalar degree term have exact rank four.  Both displayed `E` contractions
vanish as operator identities.  This validates exact microscopic rank four,
not merely rank four of a list of coefficient tensors.  At `h=0`, the exact
lower bound is withdrawn but the `E` null and at-most-four ceiling remain.

## 5. Projection, identities, folds, and nonuniform sources

Independent dual-number differentiation reproduces the uniform-source word
tensor `M_w` for every four-label occurrence composition at orders two, four,
six, and eight.  Every row lies in `span{D_a}` and annihilates both `E`
generators.

That finite word census is not itself custody for every block-resolved fold.
The load-bearing proof is the fixed-map chain rule.  The source is inserted
before Feshbach reduction; every microscopic first derivative in an `E`
direction is zero.  Differentiating an independent rational Schur complement
confirms that the effective first derivative is then zero as well.  Resolvent
insertions, scalar identities, energy derivatives, and self-consistency folds
are all outputs of the same fixed map, so none acquires an independently
assignable post-projection source tensor.  An explicit two-block replay also
confirms four internal directions per block and both local `E` nulls when the
source is nonuniform.

The microscopic map has exact rank four.  Projection can lower it.  The
theorem therefore correctly claims only
`rank(Q_eff)<=4`, not exact projected rank four.  Referring to the union of
microscopic and projected inventories as rank four does not change that
ceiling, because the exact microscopic component is retained in that union.

## 6. Quadratic contact boundary

A prospectively frozen `O(j^2)` contact can have a nonzero `EE` or mixed
`E-A1` Hessian.  I constructed such a Hessian explicitly.  Its gradient at
`j=0` is nevertheless zero, so it cannot create the missing source-off linear
`Q_E`, a spectral pole of that absent linear operator, or its canonical
commutator moment.  The theorem does not claim the contact Hessian or every
higher CTP derivative vanishes.

## 7. Exact scientific disposition

The earned result is narrow and decisive:

\[
\boxed{
\text{frozen reduced CW/FM parent + frozen additive query}
\longrightarrow \operatorname{rank}Q_{\rm micro}=4
\longrightarrow \ker Q\supset E
\longrightarrow \operatorname{rank}Q_{\rm eff}\le4.}
\]

Therefore the current reduced Q4-BLOCK additive strain query cannot satisfy
its own six-direction antecedent.  A full six-channel CTP calculation cannot
repair that missing linear source rank.  This is not a no-go theorem for a
complete physical BS parent, a root/cross-dyad or rotating-frame query, a
different collective variable, a tensor phase, RGRL-B, gravity, or `G`.

## 8. Reproducibility and custody

Run:

```bash
python3 -B LANE_GRA_FS_F3_Q4_COMPLETE_SOURCE_RANK_AUDIT_V001/independent_hostile_audit.py
```

The independent executable imports no builder code.  It verifies all eleven
dependency hashes, reconstructs the graph and algebra, checks the six-member
base payload manifest and builder seal, and applies appended-byte negative
tests.  `AUDIT_MANIFEST.sha256` and `AUDIT_SEAL.sha256` separately seal the
hostile audit bytes and transcript without recursive self-hashing.

# Record-conditioned tetrahedral holonomy witness theorem

**Lane ID:** `CROSS-RFT-GRA-ER-RECORD-TETRAHEDRAL-HOLONOMY-WITNESS-V001`

**Official short name:** `RTHW`

**Date:** 2026-08-27

**Builder status:**
`SOURCE_FROZEN_PENDING_INDEPENDENT_HOSTILE_AUDIT__UPSTREAM_EP_EQ_SEALED`

**Claim class:** exact finite compact-unitary two-probe common-frame witness;
exact identical signed-gate-inventory curved/flat comparison; exact
tetrahedral shared-face-set compatibility and natural-coframe torsion no-go;
conditional post-formation qualified-record routing composition

**Not claimed:** autonomous derivation of the routing law; an EP35 complete
stress match; lineage-specific causality; pointwise shared-face gluing; a
selected discrete Levi--Civita connection; torsion freedom; a refinement
limit; a massless spin-two pole; universal stress coupling; Einstein gravity;
Newton's constant; or outcome selection

## 1. Exact question and boundary

EQ supplies an exact bounded four-record-stream `Q4` parent once its supplied
record predicates and finite physical antecedent are satisfied.  EP proves the
exact typing of post-formation common-frame holonomy and keeps physical link
realization, face gluing, torsion, connection selection, refinement, stress
matching, and gravity separate.

This lane asks the next smaller question:

\[
 \text{Can one qualified record control an exactly unitary, two-probe}
 \quad\text{common Lorentz connection with nonzero diamond holonomy?} \tag{ER01}
\]

The answer is yes in one explicit finite model.  The construction uses only
compact spatial rotations, so it avoids the noncompact-spin-similarity versus
positive-Hilbert-unitary obstruction.  It also exposes two remaining physical
ceilings exactly: complete stress matching is not earned, and the natural
constant tetrahedral coframe is not torsion free.

The record is read only after the bounded Boolean cell has formed.  Neither the
record-controlled routing operation nor either probe is inserted into Q4's
complete-output equality.  Doing so would force the two route outputs to agree
and eliminate the curvature witness.

## 2. Exact tetrahedral rotation packet

Use the four standard tetrahedral unit vectors

\[
\begin{aligned}
 n_1&=(1,1,1)/\sqrt3,& n_2&=(1,-1,-1)/\sqrt3,\\
 n_3&=(-1,1,-1)/\sqrt3,& n_4&=(-1,-1,1)/\sqrt3.
\end{aligned}                                                   \tag{ER02}
\]

Then `n_i dot n_i=1` and `n_i dot n_j=-1/3` for `i != j`.  The null rays
`ell_i=(1,n_i)` therefore have tetrahedral Gram data

\[
 K_{ii}=0,\qquad K_{ij}=-{4\over3}\quad(i\ne j)                 \tag{ER03}
\]

in `eta=diag(-1,1,1,1)`.

Let `R_1` implement the even permutation `(234)` and let `R_2` implement
`(134)`:

\[
 R_1=
 \begin{pmatrix}0&0&1\\1&0&0\\0&1&0\end{pmatrix},
 \qquad
 R_2=
 \begin{pmatrix}0&-1&0\\0&0&1\\-1&0&0\end{pmatrix}.          \tag{ER04}
\]

They obey

\[
 R_i^{\mathsf T}R_i=I_3,\qquad \det R_i=1,\qquad
 R_i^3=I_3,\qquad \operatorname{tr}R_i=0.                      \tag{ER05}
\]

`R_1` fixes `n_1` and cycles the other three rays.  `R_2` fixes `n_2` and
cycles the other three rays.  Thus every `R_i` or `R_i^{-1}` used on a
direction-`i` link preserves the corresponding three-ray shared-face set.

Embed the rotations in the proper orthochronous Lorentz group by

\[
 \Lambda_i=\operatorname{diag}(1,R_i),\qquad
 \Lambda_i^{-1}=\operatorname{diag}(1,R_i^{-1}).                \tag{ER06}
\]

### Theorem RTHW-1 -- exact face-set-preserving Lorentz links

Every map in (ER06) lies in `SO^+(1,3)`, preserves the complete Gram packet
(ER03), fixes the link's indexed null ray, and permutes the other three rays
within the intrinsic shared-face set.

### Proof

Equation (ER05) proves proper spatial orthogonality, and the block embedding
proves Lorentz isometry, determinant one, and time orientation.  Direct action
of (ER04) on (ER02) gives `(234)` and `(134)`.  Each permutation fixes its
indexed ray and preserves all inner products. QED.

This is shared **face-set** compatibility for a symmetric cell.  It neither
chooses a pointwise face identification nor proves a unique metric-compatible
torsion-free connection.

## 3. Lawful compact unitary realization on two probes

Let bold `i,j,k` denote quaternion units.  Choose the unit quaternions

\[
 q_1={1+\mathbf i+\mathbf j+\mathbf k\over2},\qquad
 q_2={1-\mathbf i+\mathbf j+\mathbf k\over2}.                  \tag{ER07}
\]

Quaternion conjugation on imaginary vectors projects `q_1` to `R_1` and
`q_2` to `R_2`.  In a positive two-level probe block the corresponding spin
unitary is

\[
 U(q)=q_0 I-i(q_x\sigma_x+q_y\sigma_y+q_z\sigma_z),            \tag{ER08}
\]

and satisfies

\[
 U(q)U(q)^\dagger=I,\qquad
 U(q)\sigma_jU(q)^\dagger=R(q)^k{}_j\sigma_k.                  \tag{ER09}
\]

Use two independently prepared and independently read probe blocks `s=1,2`,
each receiving the same `U(q)` on every admitted link.  Their extracted
coefficient maps are therefore the same `Lambda` in (ER06).  This is a lawful
compact-unitary realization of EP's common-link gate, not merely an algebraic
noncompact spin representative.

The spin lift of a Lorentz rotation has a central sign.  Fix that sign before
scoring.  On the direction-1 bottom link use `q_1`; on its top link use
`-q_1^{-1}`.  This same signed choice is used in both lifecycle arms.  The sign
is not discarded: it is chosen so that the Lorentz-flat arm below is also flat
on the raw admitted spin-probe routes rather than hiding a central `-I` loop.

## 4. One exact record-controlled diamond

Take count directions `1,2`, base vertex `v`, and common child
`c=v+e_1+e_2`.  Freeze the four link positions before scoring.  The extracted
Lorentz maps are:

| link | KEEP | BREAK |
|---|---:|---:|
| `T_1(v)` | `Lambda_1` | `Lambda_1` |
| `T_2(v)` | `Lambda_2` | `Lambda_2^{-1}` |
| `T_1(v+e_2)` | `Lambda_1^{-1}` | `Lambda_1^{-1}` |
| `T_2(v+e_1)` | `Lambda_2^{-1}` | `Lambda_2` |

Thus KEEP and BREAK have the exact same Lorentz gate inventory

\[
 \{\Lambda_1,\Lambda_1^{-1},\Lambda_2,\Lambda_2^{-1}\}         \tag{ER10}
\]

and the same per-direction multisets.  BREAK merely exchanges the bottom and
top placements of the two direction-2 modules.

Use EP's child-based route convention:

\[
 P_{21}=T_1(v+e_2)T_2(v),\qquad
 P_{12}=T_2(v+e_1)T_1(v),\qquad
 H=P_{21}P_{12}^{-1}.                                          \tag{ER11}
\]

On the spatial block, KEEP gives

\[
 P_{21}^K=R_1^{-1}R_2=(124),\qquad
 P_{12}^K=R_2^{-1}R_1=(142),                                  \tag{ER12}
\]

and hence

\[
 \boxed{H_K=(142)\ne I},                                      \tag{ER13}
\]

with matrix

\[
 H_K=
 \begin{pmatrix}0&-1&0\\0&0&-1\\1&0&0\end{pmatrix}.        \tag{ER14}
\]

BREAK instead gives

\[
 P_{21}^B=R_1^{-1}R_2^{-1}
 =R_2R_1=P_{12}^B
 =\operatorname{diag}(-1,1,-1),                               \tag{ER15}
\]

so

\[
 \boxed{H_B=I}.                                                \tag{ER16}
\]

The selected spin lifts close the central-sign loophole.  With
`T_1(v)=q_1`, `T_1(v+e_2)=-q_1^{-1}`, and the corresponding `q_2` or
`q_2^{-1}` assignments from the table,

\[
 P_{21}^{B,\rm spin}=\mathbf j=P_{12}^{B,\rm spin},\qquad
 H_B^{\rm spin}=1.                                             \tag{ER17}
\]

The KEEP spin holonomy is noncentral and projects to (ER14).  Both independent
probe blocks therefore return the same nontrivial Lorentz conjugacy class in
KEEP and the identity in BREAK.

### Theorem RTHW-2 -- exact inventory-matched common-frame holonomy witness

The table above is a finite positive-Hilbert-unitary, two-probe realization in
which two arms have identical signed gate inventory and identical
per-direction gate multisets, yet

\[
 [H_K]\ne[H_B].                                                \tag{ER18}
\]

The difference is the relational placement and noncommutative ordering of the
same modules, not a difference in gate count or one-link spectrum.

### Proof

Equations (ER07)--(ER09) prove lawful common-probe implementation.  Exact
matrix or `A_4` multiplication gives (ER12)--(ER16).  Quaternion multiplication
with the frozen lift signs gives (ER17).  The identity is not conjugate to the
nonidentity element (ER14). QED.

## 5. Exact record-conditioned routing law and its ceiling

Let `r in {K,B}` be a prospectively qualified, retained record sector with an
authenticated custody ledger.  After formation, let it control one finite
program-register swap:

\[
 U_{\rm route}
 =|K\rangle\!\langle K|\otimes I
  +|B\rangle\!\langle B|\otimes\operatorname{SWAP}_{2,\rm bottom/top}. \tag{ER19}
\]

The two direction-2 program registers initially contain the pair
`(q_2,q_2^{-1})`; the controlled swap selects the KEEP or BREAK placement in
the table.  `U_route` is exactly unitary, uses the same two program modules,
and leaves the direction-1 pair unchanged.

### Theorem RTHW-3 -- exact conditional record-to-holonomy composition

Conditional on the supplied qualified-record sector and the declared routing
law (ER19), the record value deterministically controls the common two-probe
holonomy:

\[
 r=K\Longrightarrow[H]=[(142)],\qquad
 r=B\Longrightarrow[H]=[I].                                   \tag{ER20}
\]

This is the first exact compact-unitary record-conditioned common-frame
curvature witness in the current program.  It proves that a formed physical
memory variable **can** be made constitutive of a finite common Lorentz
connection.  It does not derive (ER19) from record formation alone and does
not yet identify lineage specifically in an EP35 intervention.

For a lineage-specific claim, the actual KEEP and BREAK prepared states must
also show EP36's prospectively qualified nonzero lineage-sector distribution
contrast with authenticated retention/quarantine handoff.  Certification must
be nondemolition/noncreating or use a causally separate randomized
subensemble.  The KEEP and BREAK adjoints must agree on the complete EP35
stress/work match algebra, and every other future-active difference must pass
the exhaustive unmatched-variable census.  This lane supplies none of those
stronger physical matches.

No step chooses a quantum outcome.  The theorem starts after `r` is a formed,
qualified record and asks what a declared physical routing law does with it.

## 6. What is and is not stress matched

Every one-link spatial rotation in (ER10) has spectrum

\[
 \{1,e^{2\pi i/3},e^{-2\pi i/3}\},                             \tag{ER21}
\]

and both arms use the exact same signed spin-unitary multiset.  Consequently
every additive one-link action that is only a class function of the individual
rotation is exactly equal between the arms.

That is not a complete physical stress match.  The bottom/top direction-2
placements, controller-record correlations, route products, and plaquette
holonomy deliberately differ.  A nonadditive plaquette energy, probe
back-reaction, pulse orientation, controller memory, work history, or boundary
response can therefore differ.  Until those ports satisfy EP35--EP36 and the
exhaustive census, (ER18) is an **inventory-matched intervention witness**, not
a lineage-isolated gravity experiment.

## 7. Shared-face pass and exact torsion failure

Because every cell carries (ER03), the intrinsic shared-face Gram condition
passes.  More strongly, every direction-`i` link in the table fixes `n_i` and
permutes the complementary three-ray face set.  This earns a symmetric
face-automorphism packet.

It does not earn EP's discrete torsion-free closure.  Use the natural constant
coframe `ell_i(v)=ell_i` and evaluate EP21 at the base for directions `1,2`.
In KEEP,

\[
 T_1(v)^{-1}\ell_2=\ell_4,\qquad
 T_2(v)^{-1}\ell_1=\ell_4,                                   \tag{ER22}
\]

so the two sides are `ell_1+ell_4` and `ell_2+ell_4`, which are unequal.  In
BREAK,

\[
 T_1(v)^{-1}\ell_2=\ell_4,\qquad
 T_2(v)^{-1}\ell_1=\ell_3,                                   \tag{ER23}
\]

so the two sides are `ell_1+ell_4` and `ell_2+ell_3`, also unequal.

### Theorem RTHW-4 -- natural-coframe Levi--Civita no-go

The exact finite witness preserves the intrinsic symmetric face data and its
face label set, but it is not the torsion-free connection of the inherited
constant tetrahedral coframe.  Therefore (ER13) is curvature of the supplied
common Lorentz connection, not demonstrated Levi--Civita/Regge curvature.

An alternate cell-dependent coframe, face identification, or connection rule
may repair torsion, but it must be constructed and must pass EP's independent
connection-selection/uniqueness gate.  It cannot be inferred from local Gram
symmetry or from (ER18).

## 8. Shortest nonduplicate route toward gravity

The exact progress is now

\[
 \boxed{
 \text{bounded qualified record}
 \longrightarrow\text{unitary common-frame link placement}
 \longrightarrow\text{nontrivial discrete Lorentz holonomy}.} \tag{ER24}
\]

The shortest remaining route is not more record machinery.  It is one
same-parent successor satisfying all of the following:

1. a face-compatible coframe and a separately selected metric-compatible
   torsion-free connection;
2. a controlled small-cell family with `H=I+O(area)`, fixed logarithm branch,
   shape regularity, and a common-probe refinement limit;
3. an EP35--EP36 stress-matched, state-qualified lineage intervention showing
   that the curvature response changes while the complete stress ledger does
   not; and
4. the already-owned common massless tensor, universal-stress, protected-pole,
   and RIEHB premises.

At that point EP supplies the continuum common-metric join and RIEHB supplies
the nonlinear leading-derivative Einstein--Hilbert back-reaction.  This lane
does not rebuild either theorem.

The natural infinitesimal screening calculation is a same-inventory family
`R_i(epsilon)=exp(epsilon J_i)`, for which a route commutator begins as

\[
 H(\epsilon)=I+\epsilon^2[J_1,J_2]+O(\epsilon^3).               \tag{ER25}
\]

Equation (ER25) is a target calculation, not a continuum claim here.  Its
coframes, faces, torsion, stress ledger, and area calibration must be supplied
before `log H/area` can be called Riemann curvature.

## 9. Controls and falsifiers

1. **Inside-Q4 control.**  Put the probes or routing programs into Q4's
   complete merger output.  Exact merger equality forces the two route outputs
   to agree; any retained nontrivial `H` then contradicts the parent.
2. **One-probe control.**  Remove either independent probe.  A route-dependent
   internal unitary is then not common geometry.
3. **Spin-center control.**  Use the naive unsigned lifts.  BREAK has Lorentz
   identity but a central spin loop `-I`; unless separately retained, that is a
   hidden fiber confound.  The signed choice leading to (ER17) removes it.
4. **Inventory control.**  The four signed link modules must be the same in both
   arms.  Replacing a module after observing the score voids the witness.
5. **Face control.**  Use a direction-`i` gate that moves `n_i`.  The exact
   shared-face-set pass is then lost even if holonomy remains nontrivial.
6. **Torsion control.**  Equations (ER22)--(ER23) must not be relabeled as a
   Levi--Civita pass.
7. **Stress control.**  Matching counts, spectra, or additive class actions
   does not match controller, work, plaquette, boundary, or total stress.
8. **Certification control.**  A measurement that creates the later response
   cannot certify record-conditioned causality on the same scored run.
9. **Finite-cell control.**  A 120-degree finite holonomy without refinement is
   not `log H/area`, a smooth curvature tensor, or gravity.

## 10. Exact disposition

Conditioned on the supplied bounded Q4 record parent and the post-formation
routing law, this lane proves an exact positive-Hilbert-unitary two-probe model
in which the same signed tetrahedral rotation modules yield nonidentity common
Lorentz holonomy in one record sector and identity holonomy in the other.  It
also proves that the natural constant coframe is not torsion free.

**Disposition:**

`COMPACT_UNITARY_TWO_PROBE_TETRAHEDRAL_LINKS_EXACT__IDENTICAL_SIGNED_GATE_INVENTORY_AND_PER_DIRECTION_MULTISET_EXACT__QUALIFIED_RECORD_CONTROLLED_ROUTING_GIVES_NONTRIVIAL_VERSUS_IDENTITY_COMMON_LORENTZ_HOLONOMY_EXACT_CONDITIONAL__INTRINSIC_FACE_AND_FACE_SET_COMPATIBILITY_EXACT__NATURAL_CONSTANT_COFRAME_TORSION_FAILS_EXACT__COMPLETE_STRESS_MATCH_LINEAGE_ISOLATION_CONNECTION_SELECTION_REFINEMENT_AND_GRAVITY_OPEN`

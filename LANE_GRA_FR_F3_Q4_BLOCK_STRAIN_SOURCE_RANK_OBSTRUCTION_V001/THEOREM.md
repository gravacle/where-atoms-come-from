# Frozen q4 block-strain source rank obstruction

**Lane ID:** `GRA-FR-F3-Q4-BSSRO-V001`

**Short name:** `BSSRO`

**Date:** 2026-08-27

**Claim class:** exact microscopic q4 edge-dyad rank and `S4` type theorem;
exact additive multi-edge linear-source closure; exact preservation of the
linear source null at `j=0` under fixed Feshbach reduction even with general
`O(j^2)` contacts; all-order null preservation for the invariant-contact
subclass; exact CTP spectral/commutator boundary; exact additive-edge-source
pass-condition obstruction; explicit complete-source ceiling

**Status:**
`FOUR_TETRAHEDRAL_MICROSCOPIC_EDGE_DYADS_SPAN_A1_PLUS_T2_RANK4__DIAGONAL_TRACELESS_E_IS_EXACT_NULL2__FQ17A_ADDITIVE_MULTI_EDGE_LINEAR_WEIGHTS_RETAIN_NULL2__FIXED_FESHBACH_PRESERVES_THE_LINEAR_E_NULL_AT_SOURCE_OFF_EVEN_WITH_GENERAL_OJ2_CONTACTS__INVARIANT_CONTACTS_PRESERVE_THE_ALL_ORDER_NULL__GENERAL_E_SEAGULL_MAY_CHANGE_CONTACT_HESSIANS_BUT_NOT_THE_LINEAR_CONJUGATE_OR_ITS_SPECTRAL_POLE__ADDITIVE_EDGE_SOURCE_SUBCLASS_FAILS_SIX_OFFSHELL_DIRECTION_PREREQUISITE__COMPLETE_NODE_PORT_BOUNDARY_CONTROLLER_LINEAR_WEIGHT_INVENTORY_OPEN__CROSS_DYAD_ROOT_SOURCE_AND_OTHER_COLLECTIVE_VARIABLES_OPEN`

**Not claimed:** that the fixed F3/q4 parent has no collective tensor phase;
that the six `A3` root dyads fail to span; that every possible nonlocal or
cross-dyad query has rank four; that a thermodynamic composite cannot acquire
a tensor pole or emergent null algebra; that the separately frozen
node/port/boundary/controller linear weights required by BS20 and FQ all lie
in the edge-dyad span; that the complete Q4-BLOCK-STRAIN-CTP source has already
been ranked; that adopted RGRL-B is false; or that gravity or `G` has been
derived or excluded.

## 1. Exact question and frozen dependencies

`FQ` froze one successor, `Q4-BLOCK-STRAIN-CTP`, and required it to pass
three outputs together: six independent off-shell strain directions, a
nondegenerate conjugate/response packet, and the independent vector-plus-
scalar null architecture.  It also froze the source-before-Feshbach rule.
For a microscopic term `H_xi` with actual one-edge support occurrences,

\[
 m_\xi^{ij}=\sum_{e\in\operatorname{supp}_1(\xi)}
 \nu_{\xi e}\,\Delta_eX^i\Delta_eX^j .            \tag{FR01}
\]

This packet asks the smallest prior question: does the already frozen
**edge-supported additive component** interrogate six independent symmetric
directions? No Hamiltonian, pair attraction, tensor projector, cross-dyad
term, blocked root-edge source, or fitted contact is added.

BS20 and FQ also require onsite, node, port, boundary, and controller linear
weights to be frozen with the complete term list. Their tensor weights are
not fixed by (FR01), and this packet does not invent them. Consequently this
is an exact obstruction for the FQ17a additive edge-supported source
subclass, not yet a rank theorem for the complete BS20 source.

The load-bearing dependency bytes are:

| role | dependency | SHA-256 |
|---|---|---|
| declared BS20/CTP parent | `LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001/MICRO_ACTION.md` | `00eba581b90fb9f0b25e3fad1362b055049824897433ba06e356ab9b1f6c76ec` |
| q4 root/dyad comparison | `LANE_GRA_FC_F3_Q4_CLIFFORD_COLLECTIVE_CONE_V001/THEOREM.md` | `28b6319e3187337da8ebef2212b030ff6e5b9f8168d9844ae172d94f3e0641a6` |
| diamond affine edge typing | `LANE_GRA_FE_F3_Q4_DIAMOND_ICE_CARRIER_JOIN_V001/THEOREM.md` | `4cc63e3e5853b4250a2a5b78256d41b83b195cf527901819a62a04ef53f8d932` |
| frozen projected parent | `LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/THEOREM.md` | `5e68e4a8c62ad89cff309781a5cb54e071092e82594e60c04e7992414bc18dbe` |
| active pair-response inventory | `LANE_GRA_FJ_F3_Q4_AUTHENTICATED_LINK_PAIR_RESPONSE_V001/THEOREM.md` | `05f4a619a6f80aa40c48570ab4035ab874426502a31468a08f435e66610bd769` |
| ice observable typing | `LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/THEOREM.md` | `cd9d6c2ad704136b8fb89dea54fdcfa3fc57d393bf37762ddd5bfe5dc8bbab98` |
| through-order-eight Feshbach parent | `LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/THEOREM.md` | `78f0687c9f597c96d235224dad45c204d12d7e6c973c270e3305a901efc75b25` |
| finite TT precursor | `LANE_GRA_FO_F3_Q4_FINITE_TT_FOUR_POINT_V001/THEOREM.md` | `44fc28edc9820d2b4ea67cef9f83beef60e53bfe6407b582ffbeccfe42f756c5` |
| precursor audit | `LANE_GRA_FO_F3_Q4_FINITE_TT_FOUR_POINT_V001/INDEPENDENT_AUDIT.md` | `84d8c02c3e560198f6a9fae04f5ee81bc72c98354e02f1f58d506a8c3171c453` |
| frozen successor specification | `LANE_GRA_FQ_F3_Q4_COLLECTIVE_METRIC_ORIGIN_SCREEN_V001/THEOREM.md` | `07445c035ed4c5167a5a20280c4db69a5101eeb71831cdeb126b29702d04b69d` |
| successor hostile correction | `LANE_GRA_FQ_F3_Q4_COLLECTIVE_METRIC_ORIGIN_SCREEN_V001/INDEPENDENT_HOSTILE_AUDIT.md` | `91aa35170432684a47278e46ee2b9d56658a43acc8acbbb84480d047cdbe6dcf` |

## 2. Four microscopic append-edge dyads have rank four

Use the exact Cartesian realization of the tetrahedral q4 append vectors,

\[
 n_1={1\over\sqrt3}(1,1,1),\quad
 n_2={1\over\sqrt3}(1,-1,-1),\quad
 n_3={1\over\sqrt3}(-1,1,-1),\quad
 n_4={1\over\sqrt3}(-1,-1,1),                    \tag{FR02}
\]

and put `D_a=n_a n_a^T`.  Their Frobenius Gram matrix is

\[
 \langle D_a,D_b\rangle=(n_a\!\cdot n_b)^2
 =\begin{cases}1&a=b,\\1/9&a\ne b.\end{cases}     \tag{FR03}
\]

It has eigenvalue `4/3` on the uniform vector and `8/9` on the three
contrasts.  Hence

\[
 \boxed{\operatorname{rank}\operatorname{span}\{D_a\}=4,
 \qquad \operatorname{span}\{D_a\}\cong A_1\oplus T_2.}     \tag{FR04}
\]

The `T2` name is not inferred from dimension.  Permuting the four append
labels gives the four-point permutation module; after removing its uniform
line its characters on
`1,(12),(12)(34),(123),(1234)` are `(3,1,-1,0,-1)`, exactly the frozen
tetrahedral `T2` convention.

For a symmetric source `j`, its microscopic contraction is

\[
 x_a(j)=j:D_a=n_a^Tjn_a.                           \tag{FR05}
\]

The simultaneous kernel of all four rows is exactly

\[
 \boxed{\mathcal N_E={\operatorname{diag}(x,y,z):x+y+z=0\},
 \qquad \dim\mathcal N_E=2.}                      \tag{FR06}
\]

Indeed, the four sign patterns in (FR02) force every off-diagonal entry to
zero and then force the trace to zero.  A basis is
`diag(1,-1,0)` and `diag(1,1,-2)`.  Under the same frozen `S4` convention this
is the missing `E` summand of
`Sym^2(V)=A1+E+T2`.

### Theorem `BSSRO-1` -- microscopic source-rank theorem

The frozen one-edge BS20/FQ17a component has four, not six, independent local
symmetric directions. Its exact two-dimensional null is the diagonal-
traceless `E` sector. Separately assigned non-edge linear weights are not part
of this rank statement.

## 3. Additive multi-edge weights cannot restore `E`

Every weight allowed by (FR01) is

\[
 m_\xi=\sum_aN_{\xi a}D_a,
 \qquad N_{\xi a}=\sum_{e:a(e)=a}\nu_{\xi e}.       \tag{FR07}
\]

Therefore `delta j_E:m_xi=0` for every `delta j_E in N_E`, independently of
the number of edges, flip-word length, block, boundary location, or repeated
occurrence.  In particular, the BS06 vertex degree term has four incident
directions and

\[
 \sum_{a=1}^4D_a={4\over3}I,                       \tag{FR08}
\]

so its additive source is purely `A1`.  A hexagon or octagon word merely
changes the nonnegative integers in (FR07).

This result survives nonuniform sources block by block.  At every block the
local contraction factors through the quotient

\[
 \operatorname{Sym}^2(V)\longrightarrow
 \operatorname{Sym}^2(V)/\mathcal N_E
 \cong A_1\oplus T_2.                              \tag{FR09}
\]

Fourier transformation, blocking, and a uniform q4 affine coframe mix blocks
but cannot increase the internal tensor rank.

### Theorem `BSSRO-2` -- additive closure

Every edge-supported microscopic term admitted by the frozen FQ17a additive
rule is exactly blind to two local `E` source directions. No post-projection
loop weighting may alter this statement without changing the query rule.

## 4. Source-before-Feshbach preserves the linear null exactly

For the additive edge-supported source subclass, separate its linear source
from the contact slot allowed by BS20a:

\[
 H[j]=\widehat H(\{x_{\xi}(j)\})+R[j],
 \qquad x_\xi(j)=j(\beta_\xi):m_\xi,
 \qquad R[j]=O(j^2).                              \tag{FR10}
\]

For any block-local `delta j_E in N_E`, the invariant part obeys the exact
finite-shift identity

\[
 \widehat H[x(j+\delta j_E)]=\widehat H[x(j)].      \tag{FR11}
\]

The complete BS20a contact `R` need not obey (FR11). What follows without any
additional contact assumption is the weaker but load-bearing operator
identity

\[
 D_EH\big|_{j=0}=0,                                \tag{FR11a}
\]

because `R=O(j^2)`. Thus an allowed quadratic contact can change a second or
higher source derivative, but it cannot manufacture a linear `E`-conjugate
operator at the source-off point.

Equations (FR10)--(FR11a) do not include an independently assigned linear
onsite/node/port/boundary/controller weight outside the four-dyad span. Such a
term would add a new first derivative and must be included in the complete
source-rank census before the full Q4-BLOCK source is scored.

With fixed `P`, `Q=1-P`, the exact energy-dependent Feshbach operator is

\[
 H_{\rm eff}(E,j)=PHP+PHQ(E-QHQ)^{-1}QHP.          \tag{FR12}
\]

The chain rule and (FR11a) prove

\[
 D_EH_{\rm eff}(E,j)\big|_{j=0}=0.                 \tag{FR13}
\]

The same first-derivative statement holds for the fixed formal
self-consistent branch, every proper resolvent insertion, every energy
derivative and fold, and the complete through-order-eight series. If, in
addition, `R[j]` itself factors only through the invariants `x_xi`, then the
stronger finite-shift identity and all-order consequence hold:

\[
 H_{\rm eff}(E,j+\delta j_E)=H_{\rm eff}(E,j),
 \qquad D_ED_{A_2}\cdots D_{A_r}H_{\rm eff}=0.      \tag{FR14}
\]

Products of microscopic sources can generate higher-rank coefficient tensors
such as `D_a tensor D_b`; within this invariant-contact subclass they still
vanish when any external leg is in `N_E`. General BS20a contacts are treated
separately below.

The effective source-conjugate operator

\[
 \mathcal Q_{\rm eff}^{ij}=-2{\partial H_{\rm eff}\over\partial j_{ij}}
 \bigg|_{j=0}                                      \tag{FR15}
\]

therefore obeys

\[
 \boxed{\delta j_{E,ij}\mathcal Q_{\rm eff}^{ij}=0
 \quad\text{for both independent }E\text{ directions}.}      \tag{FR16}
\]

This is an operator identity before a state, temperature, volume limit, or
response prescription is chosen.

### Theorem `BSSRO-3` -- linear Feshbach rank monotonicity for the frozen source

For fixed projectors and branch, the first derivative of the Feshbach map for
the additive edge-supported source at source off factors through (FR09).
General `O(j^2)` contacts do not alter that derivative, so this component's
linear source rank cannot exceed four. The stronger all-order factorization
holds only for contacts that themselves factor through the four invariants.
No rank is assigned here to still-unfrozen non-edge linear weights.

## 5. CTP derivatives, contacts, and commutator moments

For the invariant-contact subclass, the complete CTP functional obeys

\[
 W[j_++\delta j_{E,+},j_-+\delta j_{E,-}]=W[j_+,j_-]. \tag{FR17}
\]

For general BS20a contacts, (FR17) need not hold. Nevertheless the source-off
linear operator in (FR16) is zero, so the noncontact retarded commutator
kernel of two linear source operators has the factorized form

\[
 \chi=C^T\widehat\chi C,
 \qquad \operatorname{rank}\chi\le4,
 \qquad C\delta j_E=0,                             \tag{FR18}
\]

where `C:j -> (j:D_1,...,j:D_4)`.  Likewise, for

\[
 M_n^{A B}=\left\langle\left[
 (\operatorname{ad}_{H_{\rm eff}})^n\mathcal Q^A_{\rm eff},
 \mathcal Q^B_{\rm eff}\right]\right\rangle,       \tag{FR19}
\]

equation (FR16) makes every moment with an `E`-polarized leg exactly zero.
No state choice can repair an operator that vanishes before expectation.

The `O(j^2)` contact slot in BS20a requires a precise separation.

1. **Invariant contacts.** Resolvent derivatives, folds, normalization
   identities, and microscopic contacts built only from the invariants
   `x_xi` inherit (FR17) and vanish with an `E` leg.
2. **General `E`-dependent contacts.** A prospectively frozen and genuinely
   derived term such as `j_E^2 C_E` may be a legitimate BS20a contact and can
   have a nonzero instantaneous Hessian. It need not be called a new query.
   Because it is `O(j^2)`, however, its first derivative at source off is
   zero. It supplies no linear `Q_E`, no `Q_E` spectral residue, no commutator
   moment with an `E`-polarized linear leg, and no canonical conjugate for
   that missing direction. A contact chosen only after seeing the result
   would instead violate the prospective freeze.

An instantaneous Hessian must not be counted as a dynamical off-shell field
or as a rank-two Ward generator.  The two source-blind directions occur
**before** any gauge quotient; they are missing configuration directions, not
four successful first-class constraints.

### Theorem `BSSRO-4` -- contact and response boundary

Contacts that factor through the four additive invariants retain the `E`
nulls at every order. A general allowed quadratic `E` contact can change the
source Hessian, so the theorem does not set every `E`-leg CTP derivative to
zero. It still cannot satisfy the missing linear-conjugate,
commutator-moment, or `Q_E` pole requirements at source off.

## 6. Why the six `A3` root dyads do not rescue this source

The already proved six sibling roots are

\[
 \alpha_{ab}=n_b-n_a,
 \qquad 1\le a<b\le4,                              \tag{FR20}
\]

and their dyads span `Sym^2(V)` with rank six.  But

\[
 \alpha_{ab}\alpha_{ab}^T
 =D_a+D_b-n_an_b^T-n_bn_a^T.                      \tag{FR21}
\]

FQ17a assigns the two-edge microscopic support the **additive** weight
`D_a+D_b`; it does not contain the two cross terms in (FR21).  The sibling
root is a blocked two-step displacement between same-sublattice sites,
whereas `D_a,D_b` are the four microscopic append-edge dyads in the fixed
diamond parent.

Weighting a generated ring or two-step transition by (FR21) only after the
Feshbach reduction is exactly the post hoc operation which FQ's hostile audit
forbids.  A cross-dyad or explicitly blocked root-edge source remains a
possible successor if derived and frozen prospectively, but it is a different
query.  Its existence does not change the rank of the source scored here.

## 7. Main theorem and disposition

### Theorem `BSSRO-5` -- frozen-source pass-condition obstruction

The current additive edge-supported BS20/FQ17a q4 source component has exact
local linear rank four, `A1+T2`, and exact `E` nullity two. Its additive linear
microscopic operator, the first derivative of its fixed Feshbach effective
action, linear effective operators, noncontact retarded correlators, and
nested commutator moments factor through that four-dimensional quotient.
General `O(j^2)` contacts can alter higher source derivatives but not this
linear rank. Therefore an additive-edge-only source instantiation cannot
satisfy `FQ`'s prerequisite of **six independent off-shell strain-coordinate
directions before quotient**:

\[
 \boxed{
 \operatorname{rank}D_jH_{\rm eff}\le4<6
 \quad\Longrightarrow\quad
 \texttt{the additive-edge-only Q4 source subclass fails}.}       \tag{FR22}
\]

This is a useful early falsifier, not a failure of the entire gravity program
and not yet closure of the complete Q4-BLOCK source. The latter still requires
a prospective census of every onsite/node/port/boundary/controller linear
weight. It closes only the additive edge-supported microscopic block-strain
subclass. The interacting H6/H8 parent, the connected four-one-link route, a
genuinely collective loop/surface variable, a prospectively derived cross-
dyad or blocked-root source, thermodynamic tensor emergence, adopted RGRL-B,
gravity, and `G` all retain their previous status.

The scientific consequence is narrow and decisive: do not spend resources
computing a six-channel CTP spectrum from the additive edge weights alone.
The immediate successor is `Q4-COMPLETE-SOURCE-RANK-AUDIT`: freeze the
complete same-parent microscopic term list and derive every non-edge linear
weight, then rank the combined source before Feshbach. A six-channel CTP
spectrum is warranted only if that combined prospective source earns the
missing two `E` directions without relabeling a post-Feshbach transition or
adding a source-off interaction.

# Complete linear-source rank audit for the reduced F3/q4 pure-ice branch

**Lane ID:** `GRA-FS-F3-Q4-CSRAV-V001`

**Short name:** `CSRAV`

**Date:** 2026-08-27

**Claim class:** exact covering-matched periodic q4-diamond family; exact
complete term census and one prospectively frozen additive source query for
the reduced CW/FM pure-incidence parent;
exact microscopic linear-source rank and through-order-eight projected
rank ceiling/null;
exact `O(j^2)` contact boundary; explicit unreduced-BS underdetermination

**Status:**
`COVERING_MATCHED_PERIODIC_Q4_DIAMOND_FAMILY_FROZEN__COMPLETE_REDUCED_CW_FM_TERM_CENSUS__ONE_LINK_FLIPS_RANK4_A1_PLUS_T2__DEGREE_NODE_WEIGHT_SCALAR__ER0_ONSITE_INACTIVE__CARRIER_STORAGE_FORMATION_FEEDBACK_BOUNDARY_CONTROLLER_AND_PORT_SECTORS_EXCLUDED_BY_SELECTED_PARENT__SOURCE_BEFORE_FESHBACH_H6_H8_FOLDS_AND_IDENTITIES_RETAIN_E_NULL2__GENERAL_OJ2_CONTACT_MAY_HAVE_E_HESSIAN_BUT_NO_LINEAR_QE__REDUCED_BRANCH_FAILS_SIX_DIRECTION_PREREQUISITE__UNREDUCED_BS_PHYSICAL_COMPLETION_REMAINS_UNDERDETERMINED`

**Not claimed:** that the complete unreduced BS04 Hamiltonian has rank four;
that periodicity alone removes a controller or port; that F3 derives the
periodic physical identification; that a separately derived rotating
coframe, root/cross-dyad, surface, storage, carrier, feedback, node, boundary,
controller, or port weight lies in the four-dyad span; that an `O(j^2)`
contact has zero Hessian; that every possible collective variable has rank
four; a tensor pole, RGRL-B, gravity, or `G`.

## 1. Exact question and byte custody

The hostile audit of `FR` proved that the FQ17a additive edge-supported
source has rank four, but correctly refused to promote that result to the
complete BS20 source.  BS20 and FQ require every onsite, node, boundary,
controller, and port linear weight to be inventoried rather than silently
assigned zero.  This lane performs that census for exactly one already owned
branch: the reduced, closed, periodic CW/FM pure-incidence parent.  It does
not call that branch the complete physical BS parent.

The load-bearing dependency bytes are:

| role | dependency | SHA-256 |
|---|---|---|
| BS20 source and full-parent term contract | `LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001/MICRO_ACTION.md` | `00eba581b90fb9f0b25e3fad1362b055049824897433ba06e356ab9b1f6c76ec` |
| q4 diamond family | `LANE_GRA_FE_F3_Q4_DIAMOND_ICE_CARRIER_JOIN_V001/THEOREM.md` | `4cc63e3e5853b4250a2a5b78256d41b83b195cf527901819a62a04ef53f8d932` |
| q4 family audit | `LANE_GRA_FE_F3_Q4_DIAMOND_ICE_CARRIER_JOIN_V001/INDEPENDENT_AUDIT.md` | `9d7ef0419b3022dba0db1add7a46d145ebe4b6ec035f73a9b760e63b978d1b2b` |
| reduced degree-two parent | `LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/THEOREM.md` | `5e68e4a8c62ad89cff309781a5cb54e071092e82594e60c04e7992414bc18dbe` |
| CW hostile audit | `LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/INDEPENDENT_HOSTILE_AUDIT.md` | `a91caa20d16b0a1194333f9b51d96546a4ea24d55e23bf1f04c7d249641af8db` |
| through-order-eight parent | `LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/THEOREM.md` | `78f0687c9f597c96d235224dad45c204d12d7e6c973c270e3305a901efc75b25` |
| FM audit | `LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/INDEPENDENT_AUDIT.md` | `53893c7198241f0f8f6aa766f3595fb75b83d208581833c32656b28d7c7f02b9` |
| frozen CTP successor | `LANE_GRA_FQ_F3_Q4_COLLECTIVE_METRIC_ORIGIN_SCREEN_V001/THEOREM.md` | `07445c035ed4c5167a5a20280c4db69a5101eeb71831cdeb126b29702d04b69d` |
| FQ hostile audit | `LANE_GRA_FQ_F3_Q4_COLLECTIVE_METRIC_ORIGIN_SCREEN_V001/INDEPENDENT_HOSTILE_AUDIT.md` | `91aa35170432684a47278e46ee2b9d56658a43acc8acbbb84480d047cdbe6dcf` |
| additive-source obstruction | `LANE_GRA_FR_F3_Q4_BLOCK_STRAIN_SOURCE_RANK_OBSTRUCTION_V001/THEOREM.md` | `62c7aaee9433a9ffa970ff6e38bac5585200cf40d6fca2cb70477e7e1e7524eb` |
| FR hostile audit and scope repair | `LANE_GRA_FR_F3_Q4_BLOCK_STRAIN_SOURCE_RANK_OBSTRUCTION_V001/INDEPENDENT_HOSTILE_AUDIT.md` | `d2da0796cfec7cff8f1d7da5c9bc449d38acdbae089dd9778fb5f19cb6e42b88` |

## 2. One frozen covering-matched periodic family

Let

\[
 L_r=5\,2^r,\qquad r=0,1,2,\ldots .                 \tag{FS01}
\]

For every `L=L_r`, take cells `x in (Z/LZ)^3`, two vertices `A_x,B_x`,
and four q4-labelled links

\[
 e_{x,a}:A_x\longrightarrow B_{x+s_a},\qquad
 (s_1,s_2,s_3,s_4)=(0,e_1,e_2,e_3).                \tag{FS02}
\]

This is the periodic quotient already proved lawful in `FE`, with the lower
bound strengthened from `L>=4` to `L>=5`.  It has

\[
 |V(G_L)|=2L^3,\qquad |E(G_L)|=4L^3,                \tag{FS03}
\]

is simple, connected, bipartite, boundaryless, and coordination four.  A
quotient-induced noncontractible closed path needs at least `L` same-part
two-step moves and therefore at least `2L>=10` graph edges.  Hence every
simple cycle through length eight is inherited from the infinite diamond
graph: the H6 hexagons are plaquette complete, and every H8 octagon is an
ordinary local FM octagon rather than a new wrapping artifact.

Reduction modulo `L_r` makes `G_(L_(r+1)) -> G_(L_r)` an eight-sheeted graph
cover.  Thus this is one **covering-matched family of finite quotients** with
a common q4 label set and affine coframe.  It is not an inclusion of finite
Hilbert spaces, a nested physical record lineage, or a derivation of periodic
identification from F3.  The distinction is part of the theorem.

### Theorem `CSRAV-1` -- support admissibility

The family (FS01)--(FS02) is a single matched finite-quotient family on which
the CW H6 and FM H8 support premises hold exactly at every size.  It has no
geometric boundary.  This topological fact alone says nothing about whether
the broader BS parent has external ports.

## 3. Frozen reduced parent and prospective term decomposition

On each `G_L`, freeze exactly the FM symmetric-detuning slice

\[
 H=H_0+V_X,\qquad
 H_0=U_d\sum_v(d_v-2)^2,\qquad
 V_X=-h\sum_eX_e,\qquad U_d>0,                    \tag{FS04}
\]

Freeze `h != 0`.  This nonzero-flip condition is required only for the exact
microscopic lower rank bound; the `E` null and projected rank ceiling remain
valid at `h=0` as well.

with the formal CW onsite coefficient `E_R` fixed to zero.  This equation is
the complete selected Hamiltonian.  It is a reduced pure-incidence branch of
BS04, not a claim that the other BS sectors vanish in nature.

Freeze the displayed BS06 decomposition rather than refactoring it after a
response is seen:

1. one Hermitian flip term `F_(x,a)=-hX_(x,a)` for every link;
2. one formal onsite term `O_(x,a)=E_R n_(x,a)=0` for every link; and
3. one unsplit degree-square term
   `V_v=U_d(d_v-2)^2` for every vertex.

There is no standalone microscopic identity term in this decomposition.  The
constant contained algebraically inside a degree square inherits the source
factor of that complete square; it is not given an independent linear
weight.  Expanding the square into edge monomials would give an equivalent
rank bound only if all pieces retain their FQ17a additive custody.  Assigning
the constant a new tensor after the split would define a different query.

Use the tetrahedral coframe

\[
 n_1={1\over\sqrt3}(1,1,1),\quad
 n_2={1\over\sqrt3}(1,-1,-1),\quad
 n_3={1\over\sqrt3}(-1,1,-1),\quad
 n_4={1\over\sqrt3}(-1,-1,1),                    \tag{FS05}
\]

and `D_a=n_an_a^T`.  Reverse bonds at a `B` vertex have the same dyads.  Assign
each flip and formal onsite term its one-edge weight `D_a`.  The degree square
is a multi-edge incidence term with exactly the four incident directions.
This lane prospectively freezes its FQ17a occurrence multiplicities to one
for each actual incident one-link support.  That is a lawful additive query
choice, not a unique tensor forced by the source-free CW/FM Hamiltonian.  It
gives

\[
 W_v=\sum_{a=1}^4D_a={4\over3}I.                  \tag{FS06}
\]

For a spatially resolved source, take one primitive-cell block per `x`,
assign outgoing `A_x` link terms to block `x`, and assign the two displayed
node squares at `A_x,B_x` to block `x`.  This fixes the microscopic block map
before response.  Other translation-covariant regroupings cannot enlarge the
internal span because every attached tensor is still a `D_a` or `W_v`.

With `R[j]=O(j^2)` the frozen source-deformed parent is therefore

\[
\begin{split}
 H[j]={}&\sum_{x,a}\left(1-{1\over2}j(x):D_a\right)F_{x,a}\\
 &+\sum_{x,s=A,B}\left(1-{1\over2}j(x):W_v\right)V_{s,x}
 +R[j],                                             \tag{FS07}
\end{split}
\]

while the zero onsite terms contribute no operator to the linear conjugate.
Equation (FS07) is a source query: `H[j=0]=H` exactly.

## 4. Complete term and weight census

The following table distinguishes absence proved by the selected dependency
from a weight that was merely inconvenient.

| class requested by BS20/FQ | selected CW/FM status | frozen linear weight | reason |
|---|---|---|---|
| one-link flip `-hX_e` | present, nonzero | `D_a` | explicit `V_X` in (FS04) |
| onsite detuning `E_R n_e` | formal class present, coefficient exactly zero | `D_a`, multiplying the zero operator | FM freezes `E_R=0` |
| vertex degree square | present, nonzero | prospectively frozen `W_v=sum_a D_a=(4/3)I` | explicit `H_0`; this query assigns occurrence one to each of its four actual incident supports under FQ17a |
| storage occupation/copying | excluded from selected branch | none in this parent | (FS04) is the complete reduced Hamiltonian |
| carrier transfer/current | excluded from selected branch | none in this parent | same |
| formation/writer/reservoir | excluded from selected branch | none in this parent | same |
| incidence feedback | excluded from selected branch | none in this parent | same |
| independent content-node term | excluded from selected branch | none in this parent | only degree nodes occur in (FS04) |
| geometric boundary | absent | none | every `G_L` is periodic and closed |
| boundary exchange | excluded from selected branch | none in this parent | periodicity alone would not prove this; (FS04) does |
| controller/clock | excluded from selected branch | none in this parent | (FS04), not topology |
| source/reader/support port | excluded from selected branch | none in this parent | (FS04), not topology |
| H6 ring terms | generated, retained | source-before-Feshbach output in `span{D_a}` | CW plus (FS07) |
| H8 dressed-H6 and octagon terms | generated, retained | source-before-Feshbach output in `span{D_a}` | FM plus (FS07) |
| projected scalar identities/folds | generated, retained | source-before-Feshbach output in `span{D_a}` | never dropped |
| general contact | retained | no linear weight because `R=O(j^2)` | Hessian may be nonzero |

This table does **not** assign weights to the omitted BS sectors.  It proves
that they are absent from this dependency branch.  Reattaching any of them
creates a larger physical parent whose complete source is underdetermined
until its linear tensor is frozen.

### Theorem `CSRAV-2` -- completeness relative to the selected parent

The rows above exhaust every microscopic term and every projected term
through order eight in the reduced CW/FM branch, together with every linear
weight in the prospectively frozen query (FS07).  No boundary, controller,
port, storage, carrier, formation,
feedback, or extra node contribution has been set to zero merely to obtain
the rank result; those terms are not in the selected Hamiltonian.  This is a
definition-level exclusion by the displayed reduced parent, not a derivation
that those sectors vanish in a physical completion.  The theorem makes no
corresponding absence statement for unreduced BS04.

## 5. Exact complete microscopic source rank

The Gram matrix of the four edge dyads is

\[
 \langle D_a,D_b\rangle=
 \begin{cases}1,&a=b,\\1/9,&a\ne b,\end{cases}      \tag{FS08}
\]

so their span has rank four and tetrahedral type `A1+T2`.  The simultaneous
source kernel is

\[
 \mathcal N_E={\operatorname{diag}(x,y,z):x+y+z=0\},
 \qquad\dim\mathcal N_E=2,                         \tag{FS09}
\]

with basis `diag(1,-1,0)` and `diag(1,1,-2)`.  Since `W_v` is their sum, the
degree-node terms add no direction.  The zero onsite coefficient adds no
operator.  The four nonzero flip terms also attain the upper bound as an
**operator** source: in the link-occupation basis, the matrix element between
two configurations differing only on link `(x,a)` isolates `X_(x,a)` and
hence isolates `j:D_a`.  With `h != 0`, the four distinct link flips are
linearly independent.  Therefore

\[
 \boxed{\operatorname{rank}\mathcal Q_{\rm micro}=4=A_1+T_2,
 \qquad\ker\mathcal Q_{\rm micro}=E.}             \tag{FS10}
\]

This is the rank of the complete frozen **linear** source for the reduced
parent, block by block.  It is not a gauge quotient and the two blind
directions are not thereby promoted to constraints.

## 6. H6, H8, identities, folds, and resolvents

FQ requires source deformation before projection.  Apply fixed-`P_2`
Feshbach reduction to (FS07), retaining every virtual history, resolvent,
fold, and scalar identity through order eight.  For a length-`m` flip word
under one common spatially uniform source
with label occurrence counts `N_a`, `sum_a N_a=m`, its source-off logarithmic
linear tensor is, up to the overall operator coefficient,

\[
 M_w=\sum_aN_aD_a-(m-1)W_v.                       \tag{FS11}
\]

The first term comes from the `m` flip numerators.  At `E_R=0` every virtual
gap is proportional to the uniformly source-scaled degree penalty; the
second term comes from the `m-1` resolvents.  Energy derivatives and
self-consistency folds retain the same total homogeneity
`h^m/U_d^(m-1)`.  Equation (FS11) is not asserted block by block for a
nonuniform source: a local derivative of a virtual gap then carries
history-dependent scalar energy fractions.  Each such derivative is still a
linear combination of the same `D_a` and scalar `W_v` factors, which is the
only fact needed for the internal-rank conclusion.

Every scalar diagonal history has a weight of the form (FS11) and is kept.
For an elementary H6 ring missing label `d`,

\[
 M_{6,d}=2\sum_{a\ne d}D_a-5W_v=-3W_v-2D_d.       \tag{FS12}
\]

An H8 octagon has `sum_aN_a=8`; a dressed H6 word has the corresponding eight
occurrences, including its repeated link.  In both cases (FS11) applies to
every constituent history.  Summing histories and folds can change operator
coefficients but cannot leave `span{D_a}`.

Equivalently, the fixed Feshbach map obeys the chain rule: its first
derivative is a linear map of the microscopic first derivatives.  Because
both `E` contractions vanish before projection, they vanish for the complete
through-order-eight effective operator:

\[
 \delta j_E:\mathcal Q_{\rm eff}=0,
 \qquad
 \operatorname{rank}\mathcal Q_{\rm eff}\le4.     \tag{FS13}
\]

The microscopic flips retain rank four, while projection can only lower the
rank.  For the complete joint microscopic-plus-projected inventory the exact
rank is therefore four.  No generated ring or identity is an independently
assignable post-Feshbach source weight.

### Theorem `CSRAV-3` -- projected rank obstruction

All H6 rings, H8 dressed rings and octagons, scalar identities, resolvents,
and folds generated from the frozen source-deformed parent retain the two
`E` null directions.  The complete microscopic source has exact rank four;
its through-order-eight effective image has rank at most four and therefore
fails FQ's six-independent-off-shell-direction prerequisite before the
expensive CTP spectrum is computed.

## 7. The `O(j^2)` contact caveat

A general prospectively frozen BS20 contact may contain an `E`-dependent
quadratic seagull.  It can give a nonzero `EE` or mixed `E-A1` Hessian and can
alter higher CTP derivatives.  What `R[j]=O(j^2)` cannot do is supply a
source-off linear operator:

\[
 {\partial R\over\partial j_{ij}}\bigg|_{j=0}=0.  \tag{FS14}
\]

Thus it supplies neither a missing linear `Q_E`, nor a `Q_E` spectral pole,
nor a canonical commutator moment built from that absent linear conjugate.
This statement does not set the seagull Hessian to zero.  A contact that is
linear in `E` would be a new linear source datum and is not an `O(j^2)`
contact.

## 8. Exact disposition and next boundary

For the one frozen covering-matched periodic CW/FM pure-ice branch,

\[
 \boxed{
 \text{complete linear source}
 \longrightarrow A_1+T_2\text{ rank }4
 \longrightarrow E\text{ null }2
 \longrightarrow \text{FQ six-direction prerequisite fails}.} \tag{FS15}
\]

This closes the complete-source census that `FR` left open **for this reduced
branch**.  It does not close the source rank of an unreduced physical BS
completion.  If storage, carriers, record feedback, boundary exchange,
controllers, or ports are rejoined, their independently frozen linear
weights are the minimal missing data.  A source with a root/cross-dyad or an
independently rotating coframe is likewise a different prospective query,
not a repair of (FS07).

The result is proof-relevant even though negative: the frozen reduced
Q4-BLOCK-STRAIN-CTP route cannot meet its own six-coordinate antecedent, so
running its full spectral kernel would not rescue that source.  The theorem
does not exclude collective tensor physics under a differently derived
complete parent, and it makes no gravity claim.

# F3/q4 carrier-lift derivability and same-incidence no-go theorem

**Lane ID:** `GRA-FF-F3-Q4-CLDNG-V001`

**Short name:** `CLDNG`

**Date:** 2026-08-27

**Claim class:** exact factor/reachability audit; exact restricted one-carrier
F3 generator; exact current-parent detuning obstruction; exact same-incidence
incompatibility between full q4 carrier hopping and `d_*=2` diamond ice;
minimal physical-antecedent isolation

**Builder status:** `READY_FOR_INDEPENDENT_HOSTILE_AUDIT`

**Not claimed:** that either missing lift is adopted; an autonomous graph
selector; a new `K_eT_e` interaction; a second support field; a graph reward;
a source-off detuning; a thermodynamic phase; visible electromagnetism; a
tensor mode; gravity; or `G`

## 1. Frozen question and dependencies

`BQ4RSW` supplies a bounded q4 record-stream architecture with count fronts

\[
 S_N=\{m\in\mathbb N_0^4:\mathbf1^{\mathsf T}m=N\},
 \qquad |S_N|={N+3\choose3},                         \tag{FF01}
\]

complete retained word/provenance custody, and the append maps

\[
 V_{a,N}=C_{a,N}\otimes W_{a,N}\otimes J_N.          \tag{FF02}
\]

`Q4DICJ` proves that the translation-invariant completion of the incidence
between `S_N` and `S_(N+1)` is exactly the diamond net.  `CCMAC` then proves
the Schur kernel under a prospective `DETUNED-Q4-CARRIER-LIFT`, while
`Q4DICJ` states the distinct `Q4-CARRIER/EDGE-LIFT` needed to put the F3
`d_*=2` link parent on that support.

The unchanged F3 seed supplies:

- one qutrit `psi_v=span{|B>,|0>,|1>}` at every eligible vertex;
- one binary incidence link `n_e` at every eligible edge;
- content-blind carrier hopping `-t sum_e n_e T_e^psi`;
- uniform carrier onsite energy `epsilon_psi sum_v q_v^psi`;
- link detuning, degree energy, record-occupation feedback, and current-square
  feedback; and
- a symbolic complete-port slot whose concrete matrices and ownership remain
  a completion contract.

The question here is not whether one can write either lift as a new model.
It is whether the already declared q4 and F3 factors and operators derive it
without a new graph reward, unowned hardware, or an unlisted interaction.

## 2. Exact q4 factor and reachability boundary

At depth `N`, the q4 active count factor is

\[
 \mathcal Q_N=\operatorname{span}\{|m\rangle:m\in S_N\}.       \tag{FF03}
\]

It is the fixed-total subspace of four physical counter streams.  The
complete reachable code is

\[
 \mathcal H_N^{\rm code}
 =\operatorname{span}\{|m(w)\rangle_Q|w\rangle_Z|p_N\rangle_P:
 w\in\{1,2,3,4\}^N\}.                              \tag{FF04}
\]

Different words with the same count remain orthogonal in `Z`.  KEEP makes
`Z` reference-stably blind to the descended active front; it does not erase
or identify those complete histories.

By contrast, on a supplied finite F3 vertex set `V`, one fixed-content,
one-carrier sector is

\[
 \mathcal H_{\psi,1}^{(x)}
 =\operatorname{span}\{|v,x\rangle:v\in V\}
 \cong\ell^2(V),\qquad x\in\{0,1\}.                \tag{FF05}
\]

This sector is part of the tensor product of one qutrit per coexisting
vertex.  A basis isometry

\[
 \iota_N:\mathcal Q_N\longrightarrow
 \mathcal H_{\psi,1}^{(x)},
 \qquad |m\rangle\mapsto|v_m,x\rangle              \tag{FF06}
\]

exists whenever enough F3 vertex factors have first been allocated and
labelled.  Equation (FF06) is a dimension-level encoding, not an existing
q4/F3 interaction.

### Theorem `CLDNG-1` -- q4 append is not frozen-history carrier hopping

No existing BQ4 append generator is the history-blind Hermitian q4-incidence
carrier transfer used by `CCMAC`.

#### Proof

On every complete codeword, (FF02) acts as

\[
 |m(w)\rangle_Q|w\rangle_Z|p_N\rangle_P
 \longmapsto
 |m(w)+e_a\rangle_Q|wa\rangle_Z|p_{N+1}\rangle_P. \tag{FF07}
\]

It necessarily writes the append label into a fresh word/provenance slot and
allocates new scaffolding.  The prospective carrier operator instead has the
form

\[
 H_{B_N}=-t\sum_{m\in S_N}\sum_{a=1}^4
 (|m+e_a\rangle\langle m|+|m\rangle\langle m+e_a|)
 \otimes I_Z,                                      \tag{FF08}
\]

on already coexisting modes and a frozen structural-history factor.  The
word matrix element in (FF07) is `|wa><w|`, whereas (FF08) is identity on one
common `Z`.  The former is a directed formation write between depth-specific
factors; the latter is a reversible post-formation carrier hop.  They are not
the same operator or channel. QED.

In particular, the two complete routes `wab` and `wba` remain orthogonal even
when their counts coincide.  Tracing their blind history can descend to one
count front, but it does not create coherent carrier amplitudes between
coexisting sites.  A separate history-blind carrier factor or a complete-port
recombiner remains necessary.

There is also a literal node-count mismatch in the microscopic F3 seed, which
uses the same number of qutrit nodes in consecutive composition layers:

\[
 |S_{N+1}|-|S_N|={N+3\choose2}>0.                  \tag{FF09}
\]

One can pad the smaller layer to `|S_(N+1)|`, but the extra qutrits and their
all-to-all possible links must then be quarantined.  Padding is a viable
construction technique, not a derivation of the quarantine from BQ4.

## 3. Exact positive core already present in F3

The missing lift must not obscure what F3 already supplies.  Let
`G=(V_+ sqcup V_-,E)` be one prospectively fixed simple bipartite eligible
graph.  Freeze a source-off sector with:

1. one conserved F3 carrier of one fixed content `x`;
2. the formation/copy couplings off;
3. a fixed incidence word `n`, with every incidence-changing BS06/FPMH pulse
   off during this carrier hold (or an independently owned exact pin);
4. fixed record-storage eigenvalues and a covariantly matched port block on
   which every carrier-relevant onsite/port term is off or a common scalar;
   and
5. no carrier source or sink during the hold.

On a link `e={u,v}`, literal BS07 algebra gives

\[
 (J_e^\psi)^2=(T_e^\psi)^2
 =q_u^\psi+q_v^\psi-2q_u^\psi q_v^\psi.           \tag{FF10}
\]

In the one-carrier sector the collision term vanishes.  If
`G_n=(V,{e:n_e=1})` has adjacency `A_n` and degree matrix `D_n`, every
incidence, degree, record-rebate, and fixed-port term independent of carrier
position contributes one common scalar `C_(n,r,port)`.  Therefore the exact
carrier generator is

\[
 \boxed{
 H_{\psi,1}^{(x)}
 =C_{n,r,\rm port}I+\epsilon_\psi I
 +\lambda_JD_n-tA_n.}                              \tag{FF11}
\]

The other carrier content gives a unitarily identical block.

### Theorem `CLDNG-2` -- exact restricted q4 carrier generator

If the q4 incidence graph has been prospectively supplied as `G_elig` and
the saturated word `n_e=1` is fixed on every q4 edge, then

\[
 A_n=\begin{pmatrix}0&B_N^\dagger\\B_N&0\end{pmatrix}.          \tag{FF12}
\]

On a regular degree-`d` realization, (FF11) becomes

\[
 \boxed{
 H_{\psi,1}^{(x)}
 =(C+\epsilon_\psi+d\lambda_J)I
 -t\begin{pmatrix}0&B_N^\dagger\\B_N&0\end{pmatrix}.}         \tag{FF13}
\]

For the deep/periodic q4 diamond support, `d=4`.  After subtracting the common
scalar, the unchanged F3 parent therefore contains exactly the required
content-blind scalar incidence transfer.  No new hopping coefficient, graph
reward, or fitted sibling edge is required.

#### Proof

Equation (FF10) makes the current-square sum count the number of active edges
incident to the occupied carrier site, which is the diagonal matrix `D_n`.
BS09 gives `-tA_n`.  Every remaining frozen term is common in the declared
block.  Regularity gives `D_n=dI`, and bipartiteness gives (FF12). QED.

This is an exact **restricted generator**, not the complete lift.  The graph,
saturated incidence word, node-label map, preparation, hold, and ports remain
prospectively supplied.

## 4. The support solder is absent from the current composition

The microscopic F3 seed introduces a link qubit for every possible
adjacent-layer arrow.  Its graph-general successors call the subset of link
registers which exist `G_elig` and state explicitly that `G_elig` is supplied.
The q4 architecture, meanwhile, carries counts and operation histories but no
F3 node qutrit or F3 link qubit indexed by every pair `(m,m+e_a)`.

No declared cross-architecture term:

1. implements the address isometry (FF06);
2. binds a q4 append key `(m,a)` to one F3 link identity;
3. prepares or conserves the q4 sparse eligible-link mask while quarantining
   every nonedge and padded node; or
4. transfers a q4 count label into a post-formation carrier location while
   acting as identity on retained history and complete references.

Preparing an F3 incidence eigenword as part of its initial-state contract is
allowed, but is preparation, not derivation.  With the inherited BS06 link
flip active on a possible nonedge, a zero occupation is not an invariant
support mask.  PESC can use a separately formed pair-memory field `K_e` to
gate successor incidence, but BQ4 does not itself produce one retained `K_e`
factor per q4 append edge.  Importing the PESC programming mission is an
additional support-formation antecedent with its own writer, route, owner,
work, failure, and boundary census.

Accordingly q4 proves the graph relation, and F3 owns link/carrier dynamics on
a supplied eligible graph, but the current composition contains no physical
solder between the two statements.

## 5. Exact current-parent detuning obstruction

`CCMAC` requires, after a common scalar is removed,

\[
 H_{\rm FD}=
 \begin{pmatrix}0&-tB_N^\dagger\\-tB_N&\Delta_\chi I\end{pmatrix},
 \qquad \Delta_\chi>0,                              \tag{FF14}
\]

where `Delta_chi` is a child-carrier onsite offset.  It is not the BS06 link
coefficient conventionally also named `Delta`.

### Theorem `CLDNG-3` -- no positive child offset from the source-off F3 bulk

On the regular q4/diamond block, every currently declared source-off F3 bulk
term gives zero child/parent carrier detuning.

#### Proof

In the fixed-incidence one-carrier block:

- `epsilon_psi` is the same at every vertex;
- the BS06 link detuning and degree penalty act on `n`, not on carrier
  location, and are one scalar after `n` is fixed;
- the BS11 record-occupation rebate also acts on `n` and fixed storage, not on
  carrier location;
- the BS11 current-square term is `lambda_JD_n`; on the degree-four diamond
  block it is `4lambda_J I` on both bipartition classes; and
- source-off BS10 copy couplings are off, while the symbolic BS12 port slot
  has no frozen layer-staggered matrix.

Thus (FF13), after its common scalar is removed, has zero diagonal on both
parts.  It equals (FF14) only at `Delta_chi=0`, which does not support the
low-parent Schur reduction. QED.

On a regular periodic diamond hold this obstruction also survives exact
covariant elimination of auxiliary factors.  Let `R_+-` be the graph
automorphism exchanging the two FCC cosets and simultaneously exchange every
matched storage and port factor.  The uniform source-off parent commutes with
`R_+-`.  Any Feshbach projector chosen covariantly with this exchange produces
an exact effective operator which also commutes with `R_+-`.  After removal
of its common scalar, a child-only term is proportional to

\[
 \Pi_- -\Pi_+,
 \qquad
 R_{+-}(\Pi_- -\Pi_+)R_{+-}^\dagger
 =-(\Pi_- -\Pi_+),                                 \tag{FF14a}
\]

and is therefore forbidden in that symmetric block.  A conditioned
sublattice-asymmetric state, schedule, or port can break this conclusion, but
then that asymmetry and its complete ownership are precisely an additional
detuning antecedent.

The finite nonnegative slab does not repair the result.  Every parent has
degree four, while a child `c` has

\[
 d_c=|\{a:c_a>0\}|\le4.                            \tag{FF15}
\]

The current-square relative shift is

\[
 \lambda_J(d_c-4)\le0,                             \tag{FF16}
\]

For `lambda_J>0` this is boundary-dependent, nonuniform, and exactly zero in
the deep interior; for `lambda_J=0` it is identically zero.  It is never a
uniform positive child gap.  The fact that an inner microscopic layer can
participate in two scheduled slabs also does not change the one-slab bulk
operator into (FF14); using asymmetric exposure as a detuning would require a
separately frozen controller/boundary schedule and its work ledger.

A concrete BS12 completion could prospectively supply

\[
 H_{\rm stagger}=\Delta_\chi\sum_{c\in V_-}q_c^\psi,
 \qquad\Delta_\chi>0,                              \tag{FF17}
\]

but its source, maintenance, switching, work, recoil, clock, failure,
boundary, and reset ownership would have to be frozen.  Equation (FF17) is a
minimal possible antecedent, not an existing-parent theorem and not adopted
here.  It is uniform on one bipartition and therefore is not a graph-shape
reward.

## 6. Exact same-`n` incompatibility of the FD and FE uses

The two proposed lifts use different roles of F3 incidence.

1. To obtain the full q4 carrier matrix `B_N` from unchanged BS09, every q4
   eligible edge must have `n_e=1`, because its hopping coefficient is
   `-tn_eT_e`.  On the regular q4/diamond support this gives
   `d_v(n)=4`.
2. The `Q4-CARRIER/EDGE-LIFT` used by the diamond-ice theorem keeps the same
   `n_e` binary field dynamical and restricts its low manifold to
   `d_v(n)=2`.

### Theorem `CLDNG-4` -- one incidence word cannot realize both exact blocks

For `t!=0` on any nonempty regular q4/diamond support, no binary word `n`
simultaneously gives the full BS09 q4-incidence carrier transfer and belongs
to the `d_*=2` ice manifold.

#### Proof

Equality of the BS09 hopping operator with the full q4 adjacency requires
`n_e=1` on each eligible edge, hence degree four at each vertex.  Ice requires
degree two at each vertex.  Since `4!=2`, the sectors are disjoint. QED.

In a quantum superposition of ice words, `n_e` remains an operator.  Replacing
it by an expectation value does not give the exact tensor factor
`B_N \otimes I_n`; BS09 instead entangles carrier motion with the instantaneous
ice configuration.  Nor may the saturated and ice free energies be added:
they are different slices of one field.

An exact calculation of the coupled `n+psi` dynamics could still discover a
different infrared carrier kernel.  `CLDNG-4` rules out identification with
the bare simultaneous FD-full-support and FE-ice blocks; it does not rule out
a separately derived collective effective law.

One could postulate a direct eligibility-gated carrier term `K_eT_e` or a
second saturated support field while leaving `n` in the ice sector.  PESC has
already proved that either is a genuinely new kinetic ingredient.  This
packet records that option only as the boundary of the no-go.  It does not
adopt, recommend, or install it.

The two conditional calculations remain separately lawful:

- FD may use a supplied saturated q4 carrier support plus an owned detuning;
  or
- FE may use the supplied q4 eligible skeleton with dynamical `d_*=2` link
  ice.

They are not yet one simultaneous unchanged-parent phase.

## 7. Complete-port obstruction and the minimal missing antecedents

BQ4's ideal logical KEEP/BREAK circuit explicitly says that physical heat,
recoil, supply, clock, route, and control matching remains to be instantiated.
F3's BS12 symbol is likewise a completion slot, not proof that concrete port
matrices and ownership exist.  Tensoring the two parents without a cross term
retains all ports but leaves the theories dynamically independent.  Adding an
unlisted cross term would require a new ownership audit.

The narrowest common missing interface is therefore:

**`Q4-SUPPORT-SOLDER`.**  One finite, port-complete physical construction
must:

1. allocate and prospectively identify one coexisting F3 vertex mode for
   every q4 front label in the chosen two-front slab;
2. retain and quarantine all padded/off-code vertices and links;
3. bind each append key `(m,a)` to exactly one F3 eligible-link identity and
   exclude every nonedge without a fitted graph-energy reward;
4. prepare the declared carrier/incidence sector from the q4-labelled source
   and prove its reachability;
5. act as identity on the retained/future-blind q4 structural history during
   the post-formation carrier window; and
6. own every source, transducer, route, work, heat, recoil, support, clock,
   boundary, invalid, failure, quarantine, reset, and reference port.

This antecedent can be pursued as a finite programmed construction.  It need
not autonomously select diamond and it need not add a graph reward.  It is
nevertheless not contained in the present BQ4 or F3 operator census.

The lift requirements then separate exactly:

\[
\begin{array}{c|c}
\text{target}&\text{minimal still-supplied antecedent}\\ \hline
\text{FE q4 eligible edge lift}&\text{`Q4-SUPPORT-SOLDER`}\\
\text{FD detuned carrier lift}&\text{`Q4-SUPPORT-SOLDER` + owned (FF17)}\\
\text{simultaneous FD full hopping + FE ice}&
\text{an additional distinct kinetic/support field, not adopted}
\end{array}                                         \tag{FF18}
\]

The positive result (FF13) means the missing physics is not another carrier
hopping law for either separate conditional lane.  It is the physical solder,
the FD offset, and—only if simultaneity is demanded—the separation of the two
incompatible roles currently assigned to `n`.

## 8. Result and disposition

The existing theories nearly meet, but they do not yet compose automatically:

\[
 \boxed{
 \text{q4 earns the diamond incidence relation}
 \quad+\quad
 \text{F3 earns scalar carrier/link dynamics on supplied support},
 }
 \tag{FF19}
\]

while no current operator physically solders the q4 labels to F3 factors, no
source-off bulk term gives the positive FD child gap, and the same F3
incidence word cannot be both saturated and two-in/two-out.

**Disposition:**

`EXACT_Q4_FACTOR_AND_APPEND_VERSUS_CARRIER_HOP_SEPARATION__EXACT_RESTRICTED_F3_ONE_CARRIER_GENERATOR_ALREADY_SUPPLIES_SCALAR_TRANSFER__Q4_SUPPORT_SOLDER_AND_POSITIVE_CHILD_DETUNING_NOT_DERIVED__FULL_BS09_Q4_HOPPING_REQUIRES_D4_WHILE_DIAMOND_ICE_REQUIRES_D2_ON_THE_SAME_N_FIELD__SIMULTANEOUS_EXACT_BLOCKS_INCOMPATIBLE__MINIMAL_PHYSICAL_ANTECEDENTS_ISOLATED__NO_NEW_KT_TERM_SECOND_FIELD_GRAPH_REWARD_PHASE_OR_GRAVITY_ADOPTED`

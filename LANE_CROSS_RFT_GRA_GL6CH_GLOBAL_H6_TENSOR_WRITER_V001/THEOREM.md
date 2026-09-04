# GL6CH — GLOBAL SIXTH-ORDER TENSOR-WRITER THEOREM

## Status and exact scope

This packet isolates the first configuration-changing pair-memory source
vertex of the degree-two F3 parent.  It starts from the sealed `GL6AO`
classification of sixth-order locked transitions and the sealed `GL6BX`
source-before-projection convention, but independently re-enumerates all
`6!=720` histories and all local source components.

The result is an exact, owner-once, arbitrary-locked-state operator theorem
on the declared girth-six `Q4`/infinite diamond-incidence parent.
It establishes a candidate-field-dependent future writer: a physical
six-cycle transition amplitude changes when a local `T2` pair field is
changed.  It does **not** establish that the candidate field is dynamically
made by prior records, select a phase, derive a continuum, produce a pole or
causal cone, identify a Ricci operator, prove gravity, or calculate `G`.

## 1. Parent and source chart

On the declared period-four incidence quotient, and equivalently for the
finite-support linked term in the infinite incidence parent, use

\[
 D=\sum_v(k_v-2)^2,\qquad W=-\sum_eX_e,\qquad
 H(j)=U_dD+hW+\sum_v j_v^TM_v .                         \tag{CH01}
\]

At a four-port constraint node, in pair order

\[
 {cal P}=(01,02,03,12,13,23),\qquad
 (M_v)_{ab}=Z_{v,a}Z_{v,b}.                              \tag{CH02}
\]

The pure pair-tensor subspace is spanned by the three orthogonal vectors

\[
 t_1=e_{01}-e_{23},\quad t_2=e_{02}-e_{13},\quad
 t_3=e_{03}-e_{12},\qquad t_i^Tt_j=2\delta_{ij}.          \tag{CH03}
\]

Write `P_T` for its Euclidean projector and restrict the source in the main
theorem to `j_v=P_Tj_v`.  Every degree-two locked word has

\[
 M_{01}=M_{23},\quad M_{02}=M_{13},\quad M_{03}=M_{12},
 \qquad P_TM_v=0.                                        \tag{CH04}
\]

Thus a pure-`T2` source vanishes in the locked subspace itself.  It probes the
virtual formation histories rather than inserting a locked-space field by
hand.

## 2. Geometric tensor attached to a cycle vertex

Let an alternating six-cycle `c` use local ports `{a,b}` at vertex `v`, and
let `{c,d}` be the complementary port pair.  Define

\[
 \boxed{\Theta_{v,c}=e_{ab}-e_{cd}} .                    \tag{CH05}
\]

Then

\[
 P_T\Theta_{v,c}=\Theta_{v,c},\qquad
 \|\Theta_{v,c}\|^2=2,                                  \tag{CH06}
\]

and it is orthogonal to the `A1` trace and both `E2` directions.  This tensor
is fixed by the port incidence of the cycle; it does not depend on which of
the two alternating locked words is used.

## 3. Exact 720-history calculation

Number the cycle links modulo six and let their initial flip charges alternate
as `(+1,-1,+1,-1,+1,-1)`, or its negative.  For a proper nonempty prefix
`S` of a flip order, the intermediate defect energy is exactly

\[
 E(S)=|\partial_cS|>0,                                   \tag{CH07}
\]

the number of cycle vertices incident on exactly one selected link.  The
source-free Q-only contribution of an order `pi` is

\[
 w_\pi=-\prod_{r=1}^{5}{1\over E(S_r)}.                  \tag{CH08}
\]

All 720 orders give nine energy profiles, and exact summation gives

\[
 \sum_{\pi\in S_6}w_\pi=-{63\over8}.                    \tag{CH09}
\]

For completeness, the canonical Hermitian first-source derivative can be
performed before restricting to `T2`.  If `m_v^0`, `m_v^1`, and `m_v(S_r)`
are the initial, final, and intermediate six-pair words, endpoint
symmetrization gives

\[
 g_v=\sum_\pi\left(\prod_{r=1}^5{1\over E(S_r)}\right)
 \sum_{r=1}^5 {m_v(S_r)-\tfrac12(m_v^0+m_v^1)\over E(S_r)}.
                                                               \tag{CH10}
\]

The executable exhausts both alternating directions, all six vertices, all
six local cycle-port pairs, both assignments of the adjacent cycle edges,
and both exterior locked signs: 288 local contexts.  In every context,

\[
 \boxed{g_v={105\over8}e_{ab}},\qquad
 \boxed{P_Tg_v={105\over16}\Theta_{v,c}}.                \tag{CH11}
\]

Equivalently, for `j_v=Theta_{v,c}`, its intermediate score is exactly `2`
when `v` is defective and `0` otherwise, and

\[
 \partial_{j_v=\Theta_{v,c}}
 \langle s\mathbin\triangle c|H_{\rm eff}|s\rangle
 ={105\over8}{h^6\over U_d^6}.                           \tag{CH12}
\]

Equation (CH11) corrects a tempting but false full-vector shortcut.  The
expression `(105/32)(Delta M^a+Delta M^b)` has the same `A1` and `T2`
projections, but differs in `E2`.  It omits canonical endpoint
symmetrization.  Direct enumeration and the independently sealed `GL6BX`
canonical implementation both return `(105/8)e_ab`.

The `A1` derivative is universally `105/8` per cycle vertex.  It rescales the
ring amplitude but is scalar in pair space.  The non-scalar result needed
here is the `T2` projection in (CH11).

## 4. Why no lower-order fold creates the result

For the diagonal owner-once histories, the replay independently obtains, at
every locked node,

\[
 V_v^{(2)}=-M_v,
 \qquad
 V_v^{(4)}=-{4\over9}{\bf1}_6-{37\over12}M_v.             \tag{CH13}
\]

The second identity is checked on all six degree-two central words and all
`3^4` compatible simple radius-one neighborhoods, for 486 cases.  From
(CH04),

\[
 P_TV_v^{(0)}=P_TV_v^{(2)}=P_TV_v^{(4)}=0.               \tag{CH14}
\]

The local identity alone would not exclude an alternating four-cycle in an
arbitrary simple 4-regular bipartite graph.  The theorem is narrower: its
executable globally checks that the declared `Q4` incidence has zero
four-cycles, while `GL6AO` supplies the same girth-six condition for the
declared infinite parent.  Hence there is no configuration-changing locked
operator at order two or four in this parent.

Moreover, no proper prefix of the six distinct cycle flips returns to the
locked sector.  A diagonal source does not change configuration.  Therefore
every folded term through order six either factors through a nonexistent
lower-order locked transition or multiplies a lower pure-`T2` first vertex
that is zero.  Canonical metric corrections cannot change the off-diagonal
answer at this first configuration-changing order.  Equations (CH09)--(CH12)
are consequently the complete order-six off-diagonal first-`T2` source term,
not a selected-collar residue.

## 5. Global arbitrary-locked-state operator

For each undirected elementary six-cycle let

\[
 T_c=P\left(\prod_{e\in c}X_e\right)P.                   \tag{CH15}
\]

It annihilates a nonalternating locked word and exchanges the two alternating
words.  It is real and self-adjoint.  Every order-six transition between
distinct locked configurations is of this form by `GL6AO`, while the replay
checks the coefficient in both directions.

Hence the complete off-diagonal part through first order in a pure-`T2`
source is

\[
\boxed{
 [H_{\rm eff}(j_T)]_{\rm off}^{(6)}=
 -{63\over8}{h^6\over U_d^5}\sum_cT_c
 +{105\over16}{h^6\over U_d^6}
  \sum_cT_c\sum_{v\in c}j_v^T\Theta_{v,c}
 }+O\!\left({h^6\|j_T\|^2\over U_d^7},
             {h^8\over U_d^7},
             {h^8\|j_T\|\over U_d^8}\right).            \tag{CH16}
\]

The cycles are counted once.  On `Q4` there are 256 such cycles, every link
belongs to six, and the four three-port orientation classes contain 64
cycles each.  Equation (CH16) is valid on every locked basis word of the
declared parent, not only on the deterministic `GL6CC` witness or the
`GL6BX` collar.

## 6. Orientation-balanced tensor access

An elementary hexagon uses three of the four port orientations.  If `d` is
the missing port, its sum of six vertex tensors is

\[
 u_d=2\sum_{\{a,b\}\subset\{0,1,2,3\}\setminus\{d\}}
       (e_{ab}-e_{\overline{ab}}).                        \tag{CH17}
\]

In `(t_1,t_2,t_3)` coordinates, with missing-port order `d=3,2,1,0`,

\[
 u_d=2(1,1,-1),\quad2(1,-1,1),\quad
     2(-1,1,1),\quad2(-1,-1,-1).                         \tag{CH18}
\]

The exact tetrahedral identities are

\[
 \sum_du_d=0,\qquad \|u_d\|^2=24,\qquad
 u_d^Tu_{d'}=-8\ (d\ne d'),\qquad
 \boxed{\sum_du_du_d^T=32P_T}.                           \tag{CH19}
\]

Thus an orientation-balanced dense family has rank-three uniform `T2`
writer access.  The zero vector sum in (CH19) is not loss of access: the four
orientations multiply distinct cycle toggles, while their quadratic Gram is
strictly supported on all of `T2`.

The locked diagonal record image lies in `A1+E2`, whereas (CH16)--(CH19)
supply all three `T2` writer directions.  Composed algebraically, the
diagonal record vertex and the off-diagonal writer therefore access the five
non-trace directions `E2+T2`.  This is a source-access statement, not yet a
metric, spacetime, or gravity statement.

## 7. Physical meaning and remaining gate

Equation (CH16) is the first exact global mechanism in this parent by which a
candidate pair field changes what record-forming transition is written next.
It is configuration-changing, state-independent in coefficient, finite
range, and survives the owner-once sum that cancels the order-four `T2`
residues.  It therefore supplies a lawful microscopic feedback channel that
an accumulation theory may identify with a field made by retained records.

That identification has not been made by this theorem.  The next gate is to
place (CH16), the diagonal record vertex, and the global contact in one
stationary accumulated state; derive the candidate field from the same
physical degrees of freedom; and test the resulting inverse response for a
common cone and the required long-wavelength tensor form.  No graviton or
preinserted Einstein kernel is used here.

`PASS__GL6CH_GLOBAL_H6_TENSOR_WRITER__EXACT_720_HISTORY_MINUS_63_OVER_8__CANONICAL_FULL_GRADIENT_105_OVER_8_EAB__PURE_T2_PROJECTION_105_OVER_16_THETA__ARBITRARY_LOCKED_STATE_OWNER_ONCE_OPERATOR__LOWER_H0_H2_H4_T_FIRST_VERTICES_ZERO__FOUR_ORIENTATION_TETRAHEDRON_GRAM_32_PT_RANK3__CANDIDATE_FIELD_DEPENDENT_FUTURE_WRITER__NO_PHASE_RICCI_GRAVITY_OR_G`

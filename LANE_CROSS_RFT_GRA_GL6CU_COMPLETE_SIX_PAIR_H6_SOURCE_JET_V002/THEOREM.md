# GL6CU V002 — COMPLETE SIX-PAIR SOURCE-JET THEOREM THROUGH `h^6`

## Status and claim class

This packet proves an exact **effective-Hamiltonian operator theorem**.  It
computes the complete raw six-pair first-source and mixed-second-source jet
of the selected diagonal branch and every incident canonical cycle writer
through sixth hopping order.  This source-before-projection result is formed
before any `A1+E2+T2` projection.  The proof is structural on every finite
simple degree-four bipartite parent of girth at least six and is replayed
literally on the inherited period-four `Q4` F3 regulator.

It is not a stationary connected response, a Schur/Legendre quotient, or a
one-particle-irreducible (`1PI`) kernel.  In particular a raw source
direction `delta j` belongs to the dual probe space and is **not** the
field-space Ward deformation

\[
 \delta C=D_C^{-1}(k\odot\xi).
\]

No such identification is made here.

V002 preserves the V001 operator calculation and repairs only its frozen
packet verifier, geometric-support versus computed-Jet census, compact-ledger
description, and CTP branch-sign custody.  The sealed V001 target and its
`REPAIR_REQUIRED` audit remain immutable dependencies.

## 1. Parent and source convention

At every degree-four constraint node `v`, order the six pair observables as

\[
 M_v=(M_{01},M_{02},M_{03},M_{12},M_{13},M_{23}),
 \qquad M_{ab}=z_a z_b .                              \tag{CU01}
\]

For arbitrary, mutually independent, nonuniform sources use

\[
 H(j)=U_dD+hW+\sum_{v,ab}j_{v,ab}M_{v,ab},
 \qquad \eta_{v,ab}=j_{v,ab}/U_d .                    \tag{CU02}
\]

The locked projector `P`, the occupation basis, `z`, and `M` are fixed and
source-independent.  On the selected `Q4` branch the bare `P` source retains
all 768 raw coordinates and has zero source Hessian.  For a downstream
normalization `H=H0-JY`, this packet's plus-source convention means
`Y_A=-V_A` when `V_A=partial_A H_can`; the second derivative must be tracked
independently rather than assigned the same single-derivative sign.

The retained space is one selected locked branch for the diagonal inventory
and the local two-endpoint `P` span for each active incident cycle.  This is
the complete through-`h^6` operator row incident to that branch, not a
materialization of every row of the exponentially large global locked
manifold and not a stationary superposition on it.

The returned `Jet` of a coefficient `K(eta)` is

\[
 \left(K(0),\;\partial_AK(0),\;
       \partial_A\partial_BK(0)\right),               \tag{CU03}
\]

where the last entry is the actual mixed derivative, not a Taylor
coefficient divided by two.  Thus an order-`h^(2m)` first derivative carries
the physical scale `h^(2m)/U_d^(2m)`, while its energy Hessian carries
`h^(2m)/U_d^(2m+1)`.

## 2. Owner-once diagonal theorem

For a distinct flipped-edge support `S`, solve the exact
intermediate-normalized Rayleigh recurrence with every defect denominator
kept as the affine Jet

\[
 \Delta_T(\eta)=D(T)+\sum_A\eta_A
       \bigl(M_A(T)-M_A(\varnothing)\bigr).            \tag{CU04}
\]

Define its linked owner by subset Möbius inversion,

\[
 K^{\rm own}_{S,n}=K^{\rm raw}_{S,n}
   -\sum_{\varnothing\ne T\subsetneq S}K^{\rm own}_{T,n}. \tag{CU05}
\]

Because a length-`2m` return word contains at most `m` distinct edges, the
complete diagonal support through `h^6` has at most three distinct edges.
Simplicity, bipartiteness, and girth at least six leave exactly these
connected supports:

\[
\begin{array}{c|c}
 h^2 & \text{one edge}\\
 h^4 & \text{one edge, one two-edge wedge}\\
 h^6 & \text{one edge, one wedge, one three-edge star,
                     one length-three path}.
\end{array}                                            \tag{CU06}
\]

Disconnected supports factor into independent components and their linked
Jets vanish.  Equations (CU04)--(CU06) retain the endpoint locked-space
source, all denominator derivatives (`b'` and `d'` in earlier notation),
and all repeated-flip/fold contributions.  No projected-sector darkness is
used in the derivation.

## 3. Canonical alternating-cycle theorem

At sixth order the only configuration-changing locked-to-locked history is
an alternating elementary six-cycle.  For each such undirected owner take
the two endpoint states as `P` and the other 62 cycle subsets as `Q`.  The
engine solves the exact source-raised graph equation for the Bloch map
`chi`, constructs

\[
 H_B=A+B\chi,\qquad S=P+\chi^\dagger\chi,              \tag{CU07}
\]

and returns the des-Cloizeaux canonical operator

\[
 H_{\rm can}=S^{1/2}H_BS^{-1/2}.                       \tag{CU08}
\]

All quantities in (CU07)--(CU08) are sparse exact rational Jets and are
retained through source degree two and hopping degree six.  This includes
the endpoint `P` source, denominator derivatives, matching, Bloch metric,
square-root, inverse-square-root, and canonical fold terms.  Hermiticity is
verified at every hopping order.

Equivalently, differentiating the complete canonical map `F` includes
`D^2F[G_A,G_B]+DF[G_AB]`.  For (CU02) the microscopic source embedding is
affine, so bare `G_AB=0`; induced denominator, matching, fold, and canonical
second derivatives are nevertheless nonzero and retained.  A later
nonlinear physical-source embedding with its own nonzero `G_AB` is outside
this packet.

For every active cycle `c`, the source-free and first-source parts reproduce
the frozen result

\[
 (H_{\rm can}^{(6)})_{01}\big|_{\eta=0}=-{63\over8},
 \qquad
 \partial_{\eta_{v,ab}}(H_{\rm can}^{(6)})_{01}
 =\begin{cases}
 105/8,&v\in c\text{ and }(a,b)\text{ are its two ports at }v,\\
 0,&\text{otherwise}.
 \end{cases}                                           \tag{CU09}
\]

The mixed second-source writer Jet is nonzero and is frozen exactly by the
local stencil ledger, rather than inferred from (CU09).

## 4. Literal inherited `Q4` census

The selected inherited background has 128 constraint nodes, 256 links, and
64 active alternating cycle transitions.  The graph-wide geometric-support
census is

\[
 256_{\rm link}+768_{\rm wedge}+512_{\rm star}
 +2304_{\rm path}+256_{\rm six\ cycle}=4096.           \tag{CU10}
\]

These are 4096 **enumerated supports**, not 4096 computed selected-row
operator Jets.  The executable computes

\[
 3840_{\rm diagonal}+64_{\rm active\ incident\ cycles}
 =3904.                                                   \tag{CU10a}
\]

The other 192 graph-wide cycle supports are enumerated and classified as
inactive on that row; they are not canonically differentiated.  One active
undirected cycle owns both Hermitian matrix directions.  The executable
preserves every raw source coordinate and classifies the computed exact local
stencils without identifying translated or rotated sources.  No full global
locked-manifold operator is claimed.

The compact ledger commits 1654 reconstructed local stencil classes across
eight family/order records.  It stores counts, aggregate class/stencil hashes,
at most four representative metadata records, and one literal
`first_exact_stencil` per family/order record.  The sealed executable—not the
JSON text alone—reconstructs all 1654 classes and verifies those aggregate
hashes.

The aggregate diagonal source-off values are

\[
 K^{(2)}=-{256\over2},\qquad
 K^{(4)}=-{7\,256\over24},\qquad
 K^{(6)}=-{893\,256\over1080},                         \tag{CU11}
\]

matching the frozen `GL6AO` values.  Raw contractions reproduce the frozen
`GL6CF` order-two symbol, the `GL6BV` tensor contact, the `GL6CG` pointwise
order-two and order-four first-source operators, the `GL6CN` diagonal
order-six tensor darkness, and the `GL6CH` cycle writer.

## 5. Hessian versus branchwise connected contact

The displayed second derivative is the **energy Hessian**

\[
 K''_{AB}=\partial_A\partial_B H_{\rm can}.             \tag{CU12}
\]

The branch convention is defined here rather than inferred from an unpinned
label.  Let

\[
 Z[j^+,j^-]=\operatorname{Tr}\!\left(
 T e^{-i\int H(j^+)dt}\,\rho\,
 \bar T e^{+i\int H(j^-)dt}\right),
 \qquad W[j^+,j^-]=-i\log Z[j^+,j^-],                    \tag{CU13}
\]

and write `sigma=+1` on the forward branch and `sigma=-1` on the backward
branch.  Direct differentiation gives

\[
 {\delta W\over\delta j_A^\sigma(t)}
 =-\sigma\langle H_{,A}^\sigma(t)\rangle,
\]

so the same-branch coincident connected contact operator contributed by the
local source-second term is

\[
 \boxed{K^{{\rm CTP,contact},\sigma\sigma}_{AB}(t,s)
   =-\sigma K''_{AB}\,\delta(t-s),\qquad
 K^{{\rm direct},+-}_{AB}=K^{{\rm direct},-+}_{AB}=0.} \tag{CU14}
\]

Thus the forward `++` direct contact is `-K'' delta(t-s)` and the backward
`--` direct contact is `+K'' delta(t-s)`.  The pinned `GL6BV` theorem,
manifest, and seal contain the inherited forward-branch convention and are a
regression check; (CU13)--(CU14) supply the defining custody inside this
packet.  The Hamiltonian convention remains `H=H0+jM`.  Under the downstream
normalization `H=H0-JY` with the same source coordinate, `Y_A=-V_A` for
`V_A=partial_A H_can`.

Equation (CU14) does not include the two-first-vertex spectral term, state
variation, constraints, matching to a retained field, or a Schur/Legendre
quotient.  Those belong to the later complete physical response.

## 6. Disposition

What is proved is the complete owner-once selected-branch and incident-cycle
canonical Hamiltonian source Jet through `h^6`, with arbitrary nonuniform
six-pair sources, on the stated F3 parent class and with a literal full
selected-`Q4` replay.  The packet closes the
previous source-second inventory ambiguity at the operator layer.

What remains open is the stationary connected response, the same-state
physical `1PI`/quotient kernel, the map from probe space to retained-field
space, the parent longitudinal Ward null, Ricci/Einstein completion,
accumulation/refinement, gravity, and `G`.

`PASS__GL6CU_V002_COMPLETE_SELECTED_BRANCH_RAW_SIX_PAIR_CANONICAL_HAMILTONIAN_SOURCE_JET_THROUGH_H6__4096_GEOMETRIC_SUPPORTS_ENUMERATED__3904_SELECTED_ROW_NONBARE_JETS_COMPUTED__1654_STENCIL_CLASSES_HASH_COMMITTED_AND_EXECUTABLY_RECONSTRUCTED__BRANCHWISE_CTP_SIGN_DEFINED__CONNECTED_1PI_WARD_GRAVITY_G_OPEN`

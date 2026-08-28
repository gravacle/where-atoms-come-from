# Ice-projected q4/F3 hybrid tensor-response theorem

**Lane ID:** `GRA-FK-F3-Q4-IHTR-V001`

**Short name:** `IHTR`

**Date:** 2026-08-27

**Claim class:** exact local diamond-ice representation theorem; exact
ice-projected one-link/pair rank theorem; exact hybrid spatial-tensor
representation-isomorphism candidate conditional on an independently owned
scalar and nonzero sector normalizations; exact symmetric-ice Fisher-query
calculation; exact commutator and compressed linked-ring response under the
inherited sixth-order F3 ring Hamiltonian

**Not claimed:** that all six unprojected pair operators remain independent
inside ice; that the independently owned scalar has already been physically
bound to the same parent; that the frozen PMMDC intertwiners are already the
physical metric response of the actual ice query; arbitrary PMMDC state
preparation; that an active correlation is automatically a qualified record;
an autonomous/scalable q4 support; an all-orders thermodynamic theorem; a
massless tensor pole, helicity two, universal stress coupling, RGRL-B,
Einstein dynamics, gravity, or `G`

## 1. Exact question and inherited domain

`FE` identifies the q4 append-incidence bulk with the diamond net and imports
the `d_*=2` ice manifold under the explicit carrier/edge-lift antecedent.
`FH` physically binds a supplied finite q4 edge list to authenticated F3/PESC
link factors, while retaining the distinction between support memory `K_e`
and active link occupation `n_e`.  `FJ` proves that before an ice projection,
four active link factors carry six independent Walsh-pair operators and an
inherited finite response.  `PMMDC` proves the representation identity

\[
 \operatorname{Sym}^2(V)=A_1\oplus E\oplus T_2,
 \qquad V=\mathbf1^\perp\subset\mathbb R^4,        \tag{FK01}
\]

and supplies exact equivariant maps from tetrahedral one-body and pair
directions to those tensor sectors.  `CW` gives the inherited leading
diamond-ice Hamiltonian

\[
 H_{\rm ring}^{(6)}=E_{\rm scalar}P_2
 -J_6\sum_{C\in\mathcal P_6}B_C,
 \qquad
 J_6={63h^6\over8U_d^5}>0                         \tag{FK02}
\]

at symmetric detuning on supplied plaquette-complete coordination-four
diamond support.  The nonzero-response statements below use `h\ne0`, so
`J_6>0`; at `h=0` the exact representation and Fisher statements survive but
the ring response vanishes.

The support typing is load-bearing.  `FH` authenticates the q4 append links
of a finite raw slab and proves only local interior diamond neighborhoods;
its global degree-two ice sector is empty.  Therefore every **global** use of
(FK02) additionally assumes the compatible coordination-four boundary or
periodic completion, or the controlled infinite-support definition, already
required by `FE` and `CW`.  This packet does not derive that domain or claim
that any added boundary links inherit q4 append authentication.  The local
ice identities need only an interior
degree-four vertex with a compatible ice state, while the linked-ring
statements need an admitted interior hexagon and a compatible global ice
state in which that hexagon is flippable.

The narrow questions are:

1. What becomes of the one-link and pair observable representations after
   imposing the local two-in/two-out constraint exactly?
2. Can the surviving modes still span one spatial symmetric tensor without
   adding an interaction?
3. Does the inherited ring law make those modes respond across vertices?

The answers are exact but typed: one-link contrasts supply `T2`; pair
variations supply `E`, not six independent modes; an independently owned
scalar can supply `A1`; and those three representations admit an exact
candidate isomorphism to `Sym^2(V)`. In the simplest symmetric ice Fisher
query, however, the pair `E` is a genuine first-order metric tangent while
the one-link `T2` derivative vanishes. Inherited ring flips still give
nonzero `E` and `T2` operator response. This is a finite representation,
query-ceiling, and dynamics closure, not gravity.

## 2. Local ice fiber and exact `S4` decomposition

At one interior degree-four diamond vertex `v`, order the incident physical
links by their q4 append labels and define

\[
 s_a:=Z_{e_a}=1-2n_{e_a},\qquad
 j_{ab}:=s_as_b,\qquad 1\le a<b\le4.              \tag{FK03}
\]

The local ice fiber is

\[
 \Omega_{2,v}=\{s\in\{-1,+1\}^4:\sum_as_a=0\},
 \qquad |\Omega_{2,v}|=6.                         \tag{FK04}
\]

Every element has two plus and two minus entries.  Let `P_(2,v)` denote its
projector.  All identities below are operator identities after compression
by `P_(2,v)` and therefore also hold on the global ice manifold.

### Theorem `IHTR-1` -- exact collapse and decomposition

On `Omega_(2,v)`,

\[
 \sum_as_a=0,\qquad
 j_{ab}=j_{cd}\quad\text{when }\{a,b,c,d\}=\{1,2,3,4\},
 \qquad
 \sum_{a<b}j_{ab}=-2I.                            \tag{FK05}
\]

Consequently:

\[
 \operatorname{span}\{s_1,s_2,s_3,s_4\}\cong T_2,
 \qquad \dim=3,                                  \tag{FK06}
\]

while

\[
 \operatorname{span}\{j_{ab}:a<b\}\cong A_1\oplus E,
 \qquad \dim=3.                                  \tag{FK07}
\]

The pair `T2` sector present in the unprojected six-edge label module is the
kernel of the ice restriction: its three opposite-edge differences vanish.
The pair `A1` sector is the fixed constant in (FK05), so the space of
**normalized pair-state variations** is only the two-dimensional `E` sector.

The six-dimensional local permutation module, equivalently the complete
space of diagonal functions on `Omega_(2,v)`, is

\[
 \mathbb R^{\Omega_{2,v}}
 =\operatorname{span}\{I,s_a,j_{ab}\}
 \cong A_1\oplus E\oplus T_2.                    \tag{FK08}
\]

#### Proof

The first identity is the constraint.  Its square gives

\[
 0=\left(\sum_as_a\right)^2=4I+2\sum_{a<b}j_{ab},
\]

which proves the last identity in (FK05).  Also
`s_1s_2s_3s_4=+1` because exactly two signs are negative.  Hence
`j_ab j_cd=1` for complementary pairs, and since each pair observable has
eigenvalues `+/-1`, `j_ab=j_cd`.

The four one-link functions have the sole linear relation `sum s_a=0`; their
span is the standard tetrahedral representation, called `T2` in the frozen
PMMDC convention.  For pairs put

\[
 p_1=j_{12}=j_{34},\qquad
 p_2=j_{13}=j_{24},\qquad
 p_3=j_{14}=j_{23}.                               \tag{FK09}
\]

Then `p_1+p_2+p_3=-I`.  Their uniform direction is `A1`, and the two
independent centered differences are `E`.  The local six-state permutation
character on the conjugacy classes
`1,(12),(12)(34),(123),(1234)` is

\[
 (6,2,2,0,0)
 =(1,1,1,1,1)+(2,0,2,-1,0)+(3,1,-1,0,-1),        \tag{FK10}
\]

which is exactly `A1+E+T2`.  Odd one-link functions and even pair functions
have zero intersection, so their dimensions `3+3` exhaust the six-state
diagonal module.  QED.

This is a material correction to any attempt to carry FJ's unprojected
six-pair independence unchanged into the ice phase.  FJ remains exact on the
unconstrained four-link factor; its pair-only six-mode statement does not
survive `P_2`.

## 3. Exact representation candidate and symmetric-query ceiling

Define the centered pair coordinates

\[
 q_i=p_i+{1\over3}I,\qquad q_1+q_2+q_3=0.         \tag{FK11}
\]

The five nonconstant statistics `s_a` modulo their one relation and `q_i`
modulo their one relation form a basis of all normalized local probability
tangents.  Equivalently, the ice-restricted exponential family

\[
 p_{\theta,\kappa}(s)
 ={1\over Z}\exp\!\left(\sum_a\theta_as_a+
                         \sum_{i=1}^3\kappa_iq_i\right),
 \quad \sum_a\theta_a=\sum_i\kappa_i=0,          \tag{FK12}
\]

is a minimal saturated five-parameter family on the six ice states.  This is
a mathematical state family; no preparation theorem is implied.

Let

\[
 \mathbf1=(1,1,1,1)^T,\quad
 P=I_4-{1\over4}\mathbf1\mathbf1^T,\quad
 v_a=Pe_a,\quad V=\mathbf1^\perp.                 \tag{FK13}
\]

Use the two frozen PMMDC equivariant maps

\[
 \begin{aligned}
 L(x)&=P\operatorname{diag}(x)P\big|_V,
     &&x\in V,\\
 M(y)&=\sum_{a<b}y_{ab}(v_av_b^T+v_bv_a^T),
     &&y\in\mathbb R^{\mathcal E_4}.             \tag{FK14}
 \end{aligned}
\]

`L` is an isomorphism from the one-link `T2` onto the `T2` summand of
`Sym^2(V)`.  Restrict `M` to the ice-surviving pair subspace

\[
 E_{\rm pair}=\{y:y_{12}=y_{34},\ y_{13}=y_{24},\
 y_{14}=y_{23},\ \sum_{a<b}y_{ab}=0\}.            \tag{FK15}
\]

Because the full PMMDC edge map is equivariant and invertible, this
restriction is an isomorphism onto the unique `E` summand. These maps are
representation intertwiners. Their use here does **not** identify either map
with the derivative of a physical information metric.

### Theorem `IHTR-2A` -- hybrid representation-rank closure

If one separately owns a scalar `rho` which changes the common information-
metric accumulation or scale, then

\[
 \boxed{
 \mathcal T_{\boldsymbol c}(\rho,y,x)
 =c_A\rho I_V+c_EM(y)+c_TL(x):
 A_1\oplus E\oplus T_2\overset{\cong}{\longrightarrow}
 \operatorname{Sym}^2(V),\qquad c_Ac_Ec_T\ne0.} \tag{FK16}
\]

Thus the exact representation answer is **yes**:

\[
 \boxed{\text{one-link }T_2+\text{pair }E+
 \text{independently owned scalar }A_1
 =\operatorname{Sym}^2(V).}                       \tag{FK17}
\]

No new microscopic interaction is required for this algebraic span. The
three summands occur with multiplicity one, have dimensions `3,2,1`, and are
inequivalent, so the direct sum is injective and six-dimensional. Equation
(FK16) is an exact **representation-isomorphism candidate**, not a completed
physical metric solder. Choosing `c_A=c_E=c_T=1` would only choose units for
three abstract intertwiners. The physical sector normalizations, their
relative signs, and their common calibration remain to be derived.

The word **independently** is load-bearing.  A uniform pair coupling inside
the ice family adds

\[
 \lambda\sum_{a<b}j_{ab}=-2\lambda I             \tag{FK18}
\]

to the exponent and is removed by normalization.  It creates no state
tangent and no physical scale response.  The required `A1` can be supplied,
for example, by an actually owned copy/accumulation scale of the type
isolated in PMMDC, but this theorem neither derives that scalar nor proves
that it shares the pair/one-link parent, clock, ports, and calibration.

### Theorem `IHTR-2B` -- exact symmetric-ice Fisher result

The representation candidate can be tested against the most direct local
query already present in (FK12). Regard `s` as the localization statistic in
`V` and define its Fisher covariance

\[
 F_{\rm ice}(\theta,\kappa)
 =\operatorname{Cov}_{p_{\theta,\kappa}}(s)\big|_V. \tag{FK18a}
\]

At the permutation-symmetric uniform ice point,

\[
 \boxed{F_{\rm ice}(0,0)={4\over3}I_V.}           \tag{FK18b}
\]

For `x in V`, use the one-link score `xi_x(s)=sum_a x_a s_a`. For
`y in E_pair`, use the edge-pair score
`eta_y(s)=sum_(a<b)y_ab s_a s_b`. Exact differentiation gives

\[
 \boxed{
 D_xF_{\rm ice}(0,0)=0,\qquad
 D_yF_{\rm ice}(0,0)={8\over3}M(y).}              \tag{FK18c}
\]

The second coefficient uses exactly the six-edge source convention in the
definition of `eta_y`; using one coefficient per opposite matching rescales
it by one half. The identities follow because the uniform ice measure is
invariant under `s -> -s`: `ss^T` is even and the one-link score is odd. A
direct sum over the six states yields the pair formula. Its image has rank
two and is exactly `E`.

Therefore the simplest symmetric ice query has an `A1` baseline and an exact
first-order pair-`E` metric tangent, but it has a **first-order `T2` no-go**.
The mathematical five-parameter state family remains saturated, and the
one-link `T2` operators remain real and dynamically active; neither fact
makes `L(x)` the first derivative of this Fisher metric. A nonzero physical
`T2` metric response requires a derived nonsymmetric background, a
higher-order solder, or another already-owned query mechanism. Which one is
real, and the three sector calibrations in (FK16), remain open.

## 4. Exact inherited ring commutator algebra

For an elementary diamond hexagon `C`, let

\[
 U_C=\prod_{e\in C}X_e,
 \qquad B_C=P_CU_C=U_CP_C,                         \tag{FK19}
\]

where `P_C` is the diagonal projector onto the two alternating occupations
of that hexagon.  For any diagonal Walsh observable

\[
 W_A=\prod_{e\in A}Z_e,                            \tag{FK20}
\]

the flippability projector commutes with `W_A`, while

\[
 B_CW_A=(-1)^{|A\cap C|}W_AB_C.                   \tag{FK21}
\]

### Theorem `IHTR-3` -- exact one-link and pair dynamics

Under (FK02),

\[
 \boxed{
 [H_{\rm ring}^{(6)},W_A]
 =2J_6\sum_{C:\ |A\cap C|\ \mathrm{odd}}W_AB_C.} \tag{FK22}
\]

At one vertex traversed by a ring `C`, call the two ring edges `p,q` and the
two external incident edges `r,s`.  Then:

- `s_p` and `s_q` are nonconserved under `B_C`, while `s_r,s_s` commute with
  that ring term;
- the four crossing pairs `j_pr,j_ps,j_qr,j_qs` are nonconserved;
- `j_pq` and `j_rs` commute with that ring term.

In a flippable local configuration, `s_p=-s_q` and `s_r=-s_s`.  The four
crossing pairs are therefore one signed `E` direction.  Varying the local
ring-edge pair over the six choices spans all of `E`.  The corresponding
one-link differences are tetrahedral roots `e_p-e_q` and span all of `T2`.
Hence no surviving normalized sector is frozen by representation alone.

For two diagonal Walsh observables `A` and `D` (or two components having
definite ring parity), a second exact identity is

\[
 \boxed{
 [[H_{\rm ring}^{(6)},A],D]
 =-4J_6\sum_{C:\ A,D\ \mathrm{both\ odd\ on}\ C}ADB_C.} \tag{FK23}
\]

Choose an admitted hexagon `C` and a compatible ice basis state `|n>` in
which `C` is alternating, and choose crossing-pair observables at two
different vertices of `C`.  Both are odd on that hexagon, so the
`|n>`-to-`|n\triangle C>` matrix element of its term in (FK23) has a
nonzero coefficient.  Distinct hexagons connect ice basis states by distinct
symmetric differences, so that matrix element cannot be canceled by another
ring.  Under this explicit flippability premise, the full inherited
sixth-order Hamiltonian therefore has a nonzero cross-vertex
operator-response channel.

Equations (FK22)--(FK23) are exact for the complete sixth-order ring sum.
They do not assume an isolated plaquette or insert a pair-pair interaction.

## 5. Exact compressed linked-ring response

Fix one ice configuration `|n>` in which `C` is alternating and define

\[
 |\bar n\rangle=U_C|n\rangle,\qquad
 Q_C=|n\rangle\langle n|+|\bar n\rangle\langle\bar n|. \tag{FK24}
\]

In the ordered two-state basis, distinct plaquette flips leave this
two-dimensional compressed block, while `C` exchanges its two basis states.
The scalar-diagonal theorem in `CW` therefore gives

\[
 Q_CH_{\rm ring}^{(6)}Q_C
 =E_CI-J_6\sigma_x.                               \tag{FK25}
\]

Every one-link or pair observable odd on `C` compresses to

\[
 Q_COQ_C=o_O\sigma_z,
 \qquad o_O\in\{-1,+1\},                          \tag{FK26}
\]

while an even observable compresses to a scalar.  In particular,

\[
 [E_CI-J_6\sigma_x,o_O\sigma_z]
 =2iJ_6o_O\sigma_y.                               \tag{FK27}
\]

For complex energy `z` away from the two poles, define

\[
 R_E(z)={1\over z-E}-{1\over z+E}
 ={2E\over z^2-E^2}.                              \tag{FK28}
\]

The exact zero-temperature spectral response of the **compressed linked-ring
Hamiltonian** is

\[
 \boxed{
 \chi^{R,Q_C}_{OD}(z)=o_Oo_D R_{2J_6}(z).}        \tag{FK29}
\]

For two crossing-pair observables at different vertices of `C`, this is an
exact nonzero cross-vertex kernel.  At `z=i\kappa`, `\kappa>0`,

\[
 R_{2J_6}(i\kappa)
 =-{4J_6\over\kappa^2+4J_6^2}\ne0.               \tag{FK30}
\]

The compression `Q_C H Q_C` is exact, but `Q_C` need not be invariant under
the full plaquette sum because other rings may carry either basis state out
of the block.  Therefore (FK29) is not asserted to be the full many-ring
resolvent.  The full-system statement earned without that qualification is
the nonzero operator response (FK23).

The poles in (FK29) lie at a finite ring gap `2J_6`.  They are not a massless
tensor pole.  The separately imported thermodynamic `U(1)` evidence concerns
a spin-one Maxwell phase and cannot be relabeled as helicity-two gravity.

## 6. What this closes

The exact revised chain is

\[
 \boxed{
 \begin{gathered}
 \text{authenticated interior q4/F3 diamond links}+
 \text{supplied compatible }d_*=2\text{ global domain}\\
 \longrightarrow\text{local ice module }A_1\oplus E\oplus T_2\\
 \longrightarrow
 \bigl(\text{one-link }T_2\bigr)\oplus
 \bigl(\text{pair }E\bigr)\oplus
 \bigl(\text{independent scalar }A_1\bigr)
 \overset{\rm representation}{\cong}\operatorname{Sym}^2(V)\\
 \longrightarrow
 F_{\rm ice}(0,0)={4\over3}I_V,\quad
 D_EF_{\rm ice}={8\over3}M(E),\quad D_{T_2}F_{\rm ice}=0\\
 \longrightarrow\text{exact }E/T_2\text{ ring nonconservation and
 cross-vertex response}.
 \end{gathered}}                                  \tag{FK31}
\]

This closes the following parts of the prior open gates:

1. **`Q4-PAIR-SOLDER`, ice-projected representation and symmetric-query
   parts.** The pair-only six-mode solder does not survive ice, but its exact
   `E` part does. The symmetric ice Fisher query realizes that `E` tangent
   exactly. Together with one-link `T2` and a separately owned scalar, the
   abstract local spatial-tensor representation rank is complete; the
   physical metric rank is not, because the same query has zero first-order
   `T2` response.
2. **`PAIR-FIELD-DYNAMICS`, finite projected part.**  The inherited sixth-
   order ring law makes the `E` pair sector and `T2` one-link sector
   noncommuting and supplies an exact common-ring cross-vertex channel.  No
   `j-j` rescue interaction is needed.
3. **PMMDC rank interpretation.** Inside ice, the normalized state-space
   coordinates are five hybrid modes rather than six independent pair
   couplings. Representation availability is not yet physical metric
   availability: scalar ownership and the nonzero `T2` metric solder are
   distinct remaining questions.

## 7. What remains for a massless tensor or gravity claim

The following are still open and are not machinery embellishments; they are
the physics gates separating (FK31) from gravity:

1. **Physical metric solder and calibration.** Derive a nonzero `T2` metric
   response despite the symmetric-query no-go; bind an independently variable
   `A1` accumulation/scale mode to the same physical episode, query, clock,
   and ports; and derive the three sector normalizations and relative signs.
2. **State preparation and record qualification.**  Demonstrate an open,
   calibrated preparation/readout neighborhood for the five normalized ice
   modes; prove that nominated retained observables pass formation,
   distinguishability, retention, and lineage BREAK rather than merely living
   on authenticated support.
3. **Collective tensor spectrum.**  Derive the long-distance response matrix
   of the hybrid six-mode tensor, its constraints and gauge redundancy, and a
   protected linear **helicity-two** pole with nonzero residue.  A finite
   `2J_6` ring pole and the spin-one ice photon do not satisfy this gate.
4. **Gluing, refinement, and support.**  Establish compatible shared-edge
   transport, a shape-regular continuum limit, autonomous/scalable support,
   and volume-uniform control of the `O(h^8)` parent remainder.
5. **Gravity endpoint.**  Establish universal stress coupling, RGRL-B,
   nonlinear back-reaction/Einstein response, and the calibrated source law
   needed for `G`.

**Disposition:**

`EXACT_LOCAL_ICE_MODULE_A1_PLUS_E_PLUS_T2__ONE_LINK_T2_EXACT__PAIR_SPAN_A1_PLUS_E_BUT_NORMALIZED_PAIR_TANGENT_E_ONLY__HYBRID_T2_PLUS_E_PLUS_INDEPENDENT_A1_IS_AN_EXACT_REPRESENTATION_ISOMORPHISM_CANDIDATE__SYMMETRIC_ICE_FISHER_BASELINE_A1_AND_PAIR_E_TANGENT_EXACT_BUT_FIRST_ORDER_T2_ZERO__INHERITED_SIXTH_ORDER_RING_GIVES_EXACT_E_T2_COMMUTATORS_AND_NONZERO_COMMON_RING_CROSS_VERTEX_RESPONSE__PAIR_ONLY_SIX_MODE_SOLDER_DOES_NOT_SURVIVE_ICE__PHYSICAL_METRIC_SOLDER_SECTOR_CALIBRATION_SCALAR_OWNERSHIP_PREPARATION_COLLECTIVE_HELICITY2_GLUE_ALL_ORDERS_STRESS_RGRLB_GRAVITY_AND_G_OPEN`

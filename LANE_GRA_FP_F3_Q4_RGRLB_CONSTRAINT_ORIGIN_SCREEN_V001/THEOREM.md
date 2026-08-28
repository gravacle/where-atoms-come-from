# q4/F3 ice RGRL-B constraint-origin screen

**Lane ID:** `GRA-FP-F3-Q4-RCOS-V001`

**Short name:** `RCOS`

**Date:** 2026-08-27

**Claim class:** exact finite projected-ice constraint algebra; exact local
`S4` representation and relation census; exact finite incidence rank; exact
continuum principal-symbol type comparison conditional on the already supplied
Maxwell infrared identification

**Status:**
`PROJECTED_ICE_HAS_ONE_SCALAR_U1_GAUSS_SPECIES__PAIR_RELATIONS_ARE_ALGEBRAIC_NOT_NEW_FIRST_CLASS_GENERATORS__S4_MODULE_MATCH_IS_NOT_CONTINUUM_SPIN2_EQUIVALENCE__RGRLB_NOT_MICROSCOPICALLY_DERIVED_FROM_CURRENT_Q4_ICE_BRANCH__FINITE_SCREEN_NOT_THERMODYNAMIC_NO_GO`

**Not claimed:** an all-orders volume-uniform F3 phase theorem; an exact gauge
symmetry of the unprojected finite-`U_d` microscopic Hamiltonian; an autonomous
rank-two field; an emergent thermodynamic tensor constraint; a tensor pole;
helicity two; RGRL-B from below; Einstein dynamics; gravity; or numerical `G`

## 1. Question, frozen parent, and proof boundary

The current no-lab gravity route has two facts which must not be conflated.
The q4, `d_*=2` F3 branch supplies an exact degree-two ice projection and,
through orders six and eight, only closed-loop transitions plus scalar
diagonal terms.  Its supplied infrared limit is compact Maxwell.  Separately,
the adopted real-world gravity theorem requires RGRL-B: six physical spatial
metric-deformation fields, a right inverse onto the metric tangent, their
stationarity equations, and the complete independent Ward/constraint packet
needed to reduce that off-shell rank-two field to two physical helicities.

This packet asks the narrow inherited-algebra question:

> Does the already fixed q4/F3 ice parent contain the scalar-plus-vector
> first-class constraint architecture of RGRL-B, or only the one-scalar Gauss
> architecture of compact `U(1)`?

No field, interaction, metric, grid, Ward identity, or fitted coefficient is
added.  The load-bearing inputs are frozen byte-for-byte:

| role | dependency | SHA-256 |
|---|---|---|
| declared microscopic parent | `LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001/MICRO_ACTION.md` | `00eba581b90fb9f0b25e3fad1362b055049824897433ba06e356ab9b1f6c76ec` |
| exact projected ice/ring parent | `LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/THEOREM.md` | `5e68e4a8c62ad89cff309781a5cb54e071092e82594e60c04e7992414bc18dbe` |
| audit of projected parent | `LANE_CROSS_ALPHA_GRA_CW_F3_PURE_KINETIC_U1_SUPPORT_SCREEN_V001/INDEPENDENT_HOSTILE_AUDIT.md` | `a91caa20d16b0a1194333f9b51d96546a4ea24d55e23bf1f04c7d249641af8db` |
| exact local representation input | `LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/THEOREM.md` | `cd9d6c2ad704136b8fb89dea54fdcfa3fc57d393bf37762ddd5bfe5dc8bbab98` |
| audit of local input | `LANE_GRA_FK_F3_Q4_ICE_HYBRID_TENSOR_RESPONSE_V001/INDEPENDENT_AUDIT.md` | `c52eab9d701d1c6e82f1d7ec395841f4d2810e96cccbc3e2504760b6742e81e4` |
| supplied Maxwell-IR/type input | `LANE_GRA_FL_F3_Q4_MAXWELL_COMPOSITE_POLE_SCREEN_V001/THEOREM.md` | `98e2b3bc7a1c998d7839dc1a6b435cc1c8ed6d5a622ba45f63571be9ef646452` |
| audit of Maxwell/type input | `LANE_GRA_FL_F3_Q4_MAXWELL_COMPOSITE_POLE_SCREEN_V001/INDEPENDENT_AUDIT.md` | `327bf6a4476c4c6382757dc156a96c6032233d34c25c1f7935e2582acf6c607a` |
| exact order-eight current-parent input | `LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/THEOREM.md` | `78f0687c9f597c96d235224dad45c204d12d7e6c973c270e3305a901efc75b25` |
| audit of order-eight input | `LANE_GRA_FM_F3_Q4_INHERITED_TT_KERNEL_BOUNDARY_V001/INDEPENDENT_AUDIT.md` | `53893c7198241f0f8f6aa766f3595fb75b83d208581833c32656b28d7c7f02b9` |
| adopted real-world target | `FINAL_GRAVITY_REAL_WORLD_THEOREM_V001.md` | `1caabded24b861932b319ed715556a5d4123b2cff5ea3004676e12c4c76de155` |
| explicit six-equation/four-constraint target | `FINAL_GRAVITY_SAME_WORLD_COMPOSITION_GATE_MATRIX_V001.md` | `8bde0b0a8bff8561b929d8a32412d6e3a51d810c89f16f163a1f35488a3772e6` |

The finite theorems below apply to any connected supplied coordination-four
bipartite graph on which the degree-two projection and named closed-loop
effective operators are defined.  The continuum interpretation is conditional
on the already supplied `MAXWELL-IR` premise.  A finite graph cannot prove or
exclude a distinct emergent tensor redundancy in a thermodynamic limit.

## 2. Exact inherited Gauss algebra

Let `G=(V=A cup B,E)` be connected and orient every edge from `A` to `B`.
Let `B_{ve}=+1` at the `A` endpoint, `-1` at the `B` endpoint, and zero
otherwise.  With occupation `n_e in {0,1}`, put

\[
 E_e=n_e-\frac12,
 \qquad
 G_v=(BE)_v=\eta_v(d_v-2),
 \qquad
 \eta_v=\begin{cases}+1&v\in A,\\-1&v\in B.\end{cases}       \tag{FP01}
\]

The inherited degree-two sector is exactly

\[
 \Omega_2(G)=\{n:G_v(n)=0\ \text{for every }v\}.             \tag{FP02}
\]

Because all `E_e` commute,

\[
 [G_v,G_w]=0.                                                \tag{FP03}
\]

For a connected graph the oriented incidence matrix has

\[
 \operatorname{rank}B=|V|-1,
 \qquad \sum_vG_v=0,                                        \tag{FP04}
\]

and (FP04) is the sole linear dependence.  In the unsigned degree variables
`Q_v=d_v-2`, the same dependency is

\[
 \sum_{v\in A}Q_v-\sum_{v\in B}Q_v=0.                       \tag{FP05}
\]

Let `T_e^+` and `T_e^-` change `E_e` by `+1` and `-1`.  Then

\[
 [G_v,T_e^\pm]=\pm B_{ve}T_e^\pm.                            \tag{FP06}
\]

A directed closed-loop operator contains one entering and one leaving ladder
operator at every visited vertex, so

\[
 [G_v,R_C]=[G_v,R_C^\dagger]=0.                              \tag{FP07}
\]

The exact sixth-order hexagon term and the exact eighth-order octagon and
dressed-hexagon terms are all of this closed-loop type; their diagonal pieces
are scalar.  Hence the controlled effective Hamiltonian obeys

\[
 [G_v,H_{\rm eff}^{(6)}+H_{\rm eff}^{(8)}]=0.                 \tag{FP08}
\]

More generally, for the exact Feshbach operator restricted to `P_2`,
`G_vP_2=P_2G_v=0`, so any `P_2H_effP_2` preserves (FP02).  Constraint
preservation therefore returns `dot G_v=0`; it does not generate a second
independent constraint species.

This statement must be kept separate from the unprojected parent.  Its
single-link tunnelling includes `-h sum_e X_e`, and (FP06) gives

\[
 [G_v,-h\sum_eX_e]\ne0.                                     \tag{FP09}
\]

At finite `U_d`, exact eigenstates can contain virtual degree defects.  Thus
the one-scalar Gauss law is exact in the projected effective theory and is an
emergent low-energy constraint of the unprojected F3 parent; the full
microscopic action has not supplied an exact tensor gauge algebra.

### Theorem `RCOS-1` -- inherited finite constraint algebra

On every connected supplied graph, the current controlled q4/F3 ice branch
has one Abelian scalar incidence-constraint species per vertex, with the one
global dependency (FP04).  Its inherited closed-loop dynamics preserves that
algebra through the exact owned order eight.  Neither preservation nor the
Dirac consistency step produces an independent vector constraint or another
scalar constraint.

## 3. Exact degree-of-freedom comparison

The following finite canonical count belongs to the standard compact-rotor
completion of the link algebra, not literally to the dimension of the
hard-core six-state Hilbert fiber.  For a connected graph, put
`r=|V|-1`.  The link rotor has `2|E|` phase-space dimensions.  Reduction by
the `r` first-class Gauss generators gives

\[
 \dim\Gamma_{\rm phys}=2(|E|-r).                             \tag{FP10}
\]

On a periodic diamond graph with `N_c` primitive cells,

\[
 |V|=2N_c,\qquad |E|=4N_c,\qquad r=2N_c-1,                  \tag{FP11}
\]

so the configuration count is `2N_c+1`: two transverse configurations at
each nonzero momentum and three global harmonic configurations at zero
momentum.  Winding sectors and the extra periodic zero-mode count are global
features, not extra local first-class constraints.

In the continuum Maxwell phase, excluding temporal multipliers,

\[
 (A_i,E^i):6\text{ phase dimensions}
 \xrightarrow{\ \partial_iE^i=0\ \text{first class}\ }
 4\text{ phase dimensions}=2\text{ configurations}.         \tag{FP12}
\]

For a symmetric spatial tensor,

\[
 (h_{ij},\pi^{ij}):12\text{ phase dimensions}
 \xrightarrow{\ \mathcal H_i=0,\ \mathcal H_\perp=0\ }
 4\text{ phase dimensions}=2\text{ configurations},         \tag{FP13}
\]

because the three vector and one scalar constraints are independent and
first class.  The equal final number `2` is not an equivalence: Maxwell starts
with three configuration components and quotients one scalar gauge function;
the rank-two packet starts with six and quotients an independent vector plus
scalar packet.  Nonlinear gravity further requires closure of the
hypersurface-deformation algebra and propagation of all four constraints.

### Theorem `RCOS-2` -- the equal-polarization false positive

The q4/F3 projected constraint count is exactly the compact-`U(1)` count.
It cannot be promoted to the RGRL-B count merely because both Maxwell and a
massless rank-two theory have two propagating polarizations.

## 4. Local `S4` census and the continuum type mismatch

At one coordination-four vertex the four link variables form

\[
 \mathbb R^4_{\rm link}=A_1\oplus T_2.                       \tag{FP14}
\]

The scalar ice equation removes `A1`, leaving the three-dimensional one-link
`T2`.  With `s_a=1-2n_a` and `j_{ab}=s_as_b`, exact ice projection gives

\[
 \begin{aligned}
 j_{12}&=j_{34},&j_{13}&=j_{24},&j_{14}&=j_{23},\\
 \sum_{a<b}j_{ab}&=-2.
 \end{aligned}                                               \tag{FP15}
\]

Consequently

\[
 \operatorname{span}\{j_{ab}\}=A_1\oplus E,
 \qquad \delta j\text{ at fixed normalization}=E,           \tag{FP16}
\]

while all diagonal functions on the six ice states obey

\[
 \mathbb R^{\Omega_2}=A_1\oplus E\oplus T_2.                \tag{FP17}
\]

Abstractly, for the tetrahedral standard vector module `V`,

\[
 \operatorname{Sym}^2(V)=A_1\oplus E\oplus T_2.             \tag{FP18}
\]

Equations (FP17)-(FP18) are a real finite-group isomorphism, but they are not
a physical tensor solder.  Under continuum rotations:

- the actual one-link `T2` is the complement-odd polar-vector `ell=1` flux;
- a symmetric tensor is `Sym2(ell=1)=ell=0 plus ell=2`;
- the even `ell=2` sector restricts to `E plus T2` under `S4`; and
- ice projection retains only its local pair `E`, while its even pair `T2`
  is identically killed by (FP15).

The two `T2` copies are isomorphic after restriction to `S4`, but have
different continuum rotation characters.  Around a propagation axis,

\[
 \chi_{\rm photon}(\theta)=2\cos\theta,
 \qquad
 \chi_{\rm TT}(\theta)=2\cos2\theta.                         \tag{FP19}
\]

At `theta=pi/2` these are `0` and `-2`.  Complement parity, continuum spin,
spectral support, and canonical type all agree that the one-link `T2` cannot
be renamed as the missing tensor `T2`.

There is another finite-group trap.  A scalar plus vector constraint packet
has local `S4` type `A1 plus T2`, and four samples or the first spatial jet of
one scalar can also be arranged as `A1 plus T2`.  At fixed nonzero Fourier
momentum, however, one scalar amplitude `chi(k)` and its gradient
`ik_i chi(k)` have symbol rank one.  An independently specifiable lapse plus
shift packet has rank four.  Spatial derivatives of one scalar gauge function
are not three new gauge functions.

### Theorem `RCOS-3` -- representation/type boundary

The local q4 ice function algebra has enough entries to imitate the
six-dimensional `S4` character of a symmetric tensor.  Its physical odd and
even subspaces do not have the `O(3)`, parity, canonical, constraint, or
spectral type of one rank-two field.  Therefore local `S4` rank is not a
microscopic RGRL-B realization.

## 5. Why the four pair identities are not four first-class constraints

The four relations in (FP15) are especially tempting because their formal
coefficient space has `A1 plus T2`, the same `S4` character as one scalar plus
one vector.  They nevertheless do not generate gauge redundancy.

First, every `j_ab` is a function of electric/occupation variables and is
invariant under the already inherited gauge action:

\[
 [G[\lambda],j_{ab}]=0.                                      \tag{FP20}
\]

This makes `j_ab` an observable of the existing gauge algebra, not a generator
of another gauge orbit.

Second, after restriction to the six-state ice Hilbert space, every relation
operator in (FP15) is literally the zero operator.  Its commutator and
Hamiltonian vector field are zero.  It cannot remove a canonical pair from a
six-component tensor phase space because that tensor phase space and its six
conjugate momenta were never derived.

Third, products and nonlinear functions of the existing `G_v` do not add
independent gauge directions.  On the constraint surface, the differential of
a smooth function of `G` is a linear combination of the `dG_v`; a product with
no linear part has vanishing differential there.  Closed-ring interactions
change gauge-invariant response functions but do not enlarge the generator
span in (FP03)-(FP04).

If the six pair symbols were promoted by hand to independent canonical
coordinates, a symplectic form, conjugate variables, gauge action, and
first-class closure would all be new inputs.  The current parent supplies none
of them, so no first-class classification can be inherited from the four
algebraic identities alone.

### Theorem `RCOS-4` -- bilinear constraint no-go on the owned fiber

No local pair or bilinear composite of the present q4 ice variables inherits
an independent first-class constraint from the current finite algebra.  The
four pair relations are exact kinematic consequences of the one scalar ice
law and `s_a^2=1`; they are not a scalar-plus-vector gauge packet.

## 6. Principal-symbol Ward discriminator

At a nonzero continuum momentum `k`, Maxwell has the scalar constraint symbol

\[
 k_iE^i=0,\qquad \operatorname{rank}=1.                       \tag{FP21}
\]

For a symmetric spatial tensor the momentum-constraint symbol

\[
 (C\pi)^i=k_j\pi^{ij}                                        \tag{FP22}
\]

has rank three, while the independent scalar linear-curvature row

\[
 C_\perp(h)=k_ik_jh^{ij}-|k|^2h^i{}_i                       \tag{FP23}
\]

has rank one.  The spatial gauge map

\[
 (D\xi)_{ij}=k_i\xi_j+k_j\xi_i                              \tag{FP24}
\]

has rank three and is annihilated by (FP23).  Transversality plus trace
removes four of the six symmetric components, leaving a two-dimensional TT
space.

The corresponding off-shell identities are likewise different in type.  A
`U(1)` Ward identity is generated by one scalar gauge function.  A covariant
rank-two response has one vector identity, whose spatial canonical split is
the independent vector-plus-scalar constraint packet.  Neither (FP21) nor
arbitrary derivatives of it supply (FP22)-(FP24).

The executable replay verifies these ranks exactly at generic integer
`k=(1,2,3)`.  Rank is stable on the open set `k ne 0`; the special zero mode
is separately identified as a global finite-volume issue.

## 7. Main theorem and promotion ceiling

### Theorem `RCOS-5` -- current-parent RGRL-B constraint-origin screen

Under the frozen q4 degree-two projection and owned closed-loop effective
dynamics, the only inherited local first-class architecture is the
one-scalar Abelian Gauss redundancy of compact `U(1)`.  The exact local pair
identities do not add an independent first-class generator.  The same local
`S4` character as a symmetric tensor or as a scalar-plus-vector packet is a
finite-group coincidence defeated by parity, continuum rotation type,
canonical rank, and the principal-symbol test.

Therefore the current controlled q4/F3 ice branch has **not** derived from
below the six-field plus independent scalar/vector Ward-constraint content of
RGRL-B.  This does not retract the adopted RGRL-B real-world theorem premise;
it locates its present microscopic derivation boundary exactly.

What is proved, imported, and open is:

| statement | status |
|---|---|
| six local ice states and pair identities | exact finite theorem |
| one scalar Abelian incidence algebra and sole global dependency | exact on every connected supplied finite graph |
| closed hexagon/octagon preservation through owned order eight | exact on the admitted finite support/order |
| unprojected single-link tunnelling breaks exact Gauss | exact microscopic commutator |
| local `S4` module/type mismatch | exact finite theorem |
| Maxwell one-scalar versus rank-two scalar-plus-vector symbol ranks | exact algebraic comparison |
| Maxwell continuum interpretation | imported `MAXWELL-IR` premise |
| thermodynamic survival or birth of another redundancy | open |
| tensor pole, common cone, nonlinear closure, RGRL-B from below | open |

This is not an all-architectures or thermodynamic no-go.  It rules out the
specific claim that already available local q4 ice links/pairs secretly
contain the required RGRL-B constraint packet.

## 8. Minimal lawful successor

The next calculation should not add record machinery or a rescue interaction.
It should jointly test the already fixed
`H_eff^(6)+H_eff^(8)` parent in the same thermodynamic analysis:

1. compute the commonly normalized connected even/TT four-point function and
   channel-two-particle-irreducible kernel;
2. solve the finite-volume spectrum or Bethe-Salpeter problem and test whether
   a nonzero-residue common-cone tensor pole survives;
3. simultaneously extract the null directions of its effective action and
   require an independently derived rank-three vector symbol plus an
   independent scalar constraint, with constraint-preserving closure; and
4. reject the local composite route if the pole dissolves into the two-photon
   continuum or the required null/constraint ranks fail.

If that fixed-parent calculation fails, the minimal distinct same-parent
successor is a genuinely collective loop/surface relational variable with an
inherited conjugate, not another local pair relabeling.  Its six-component
off-shell field, vector-plus-scalar constraints, closure, and Ward identity
must all be derived before it can instantiate RGRL-B.  No such architecture is
inserted here.

## 9. Disposition

\[
 \boxed{
 \begin{gathered}
 \text{q4/F3 degree-two projection}
 \longrightarrow \text{one scalar ice Gauss law}
 \longrightarrow \text{compact spin-one }U(1),\\
 \text{local pairs}
 \longrightarrow \text{gauge-invariant }A_1\oplus E
 \not\longrightarrow \text{new first-class constraints},\\
 \text{six local }S_4\text{ slots}
 \not\longrightarrow \text{RGRL-B rank-two Ward packet}.
 \end{gathered}}
\]

This narrows the no-lab gravity program to physics in the fixed collective
dynamics: a tensor pole and its protecting constraint algebra must emerge
together.  Neither can be supplied by renaming the Maxwell photon or the local
pair identities.


# Minimal translation-owning recoil parent for the encoded F3 hold

**Lane:** `LANE_GRA_GD_F3_Q4_TRANSLATION_OWNING_RECOIL_PARENT_V001`  
**Short name:** `TORP`  
**Date:** 2026-08-28  
**Plan gate:** `B1`  
**Claim class:** exact autonomous link-support/reservoir recoil completion;
exact factor-edge-local translation generator on one auxiliary common
mechanical torus and equal-and-opposite momentum ledger;
exact full-code reduction to GA/FV/FY nonidentity physics modulo one common
reference scalar; exact code-independent same-mechanical-state hold no-go for
nonzero recoil

**Not claimed:** a derived recoil scale or placement, a physical strain solder,
a complete spacetime source, `T^{0j}`, a stress Ward identity, a protected
tensor pole, a common cone, gravity, or `G`.

## 1. Frozen B1 question

GA supplied a charge-conserving link/reservoir encoder and an internal scalar
`U(1)` current, but did not locate a spatial current or own mechanical recoil.
GB then proved that finite translation of a fixed numerical diamond support
does not supply local physical momentum.  B1 therefore asks for the smallest
autonomous completion in which the already-required link support, reservoir,
controller, and boundary factors own translations and every active exchange
has an equal-and-opposite impulse, while the complete encoded hold still
reduces to GA/FV/FY.

The construction below adds no tensor field and no new record mechanism.  It
dresses the already-present GA charge transfer by the minimal relative
translation operator required to put its recoil somewhere.

## 2. Translation-owning factors and autonomous Hamiltonian class

Give every admitted factor `a` in the declared support/reservoir/controller/
boundary ledger the same auxiliary mechanical torus

\[
 \mathbb T^3=\mathbb R^3/\Lambda,
 \qquad \mathcal H_a^{\rm mech}=L^2(\mathbb T^3),
 \qquad P_a=-i\hbar\nabla_{r_a}.                    \tag{GD01}
\]

The exponentials of `r_a` and the self-adjoint `P_a` are the exact Weyl
translation pair; no globally single-valued angle operator is needed.  For
each link support `L_e` and its GA reservoir `R_e`, choose
`kappa_e in 2 Lambda^*`, so `+/-hbar kappa_e/2` are allowed torus momenta,
and define

\[
 U_e=\exp[i\kappa_e\!\cdot(r_{R_e}-r_{L_e})].       \tag{GD02}
\]

It obeys

\[
 [P_{L_e}^i,U_e]=-\hbar\kappa_e^iU_e,
 \qquad [P_{R_e}^i,U_e]=+\hbar\kappa_e^iU_e,
 \qquad[P_{L_e}+P_{R_e},U_e]=0.                    \tag{GD03}
\]

Write `A_e=sigma_e^+ T_(e,-)` for GA's forward charge transfer.  Replace only
the already-present inherited flip by

\[
 \mathsf X_e=A_eU_e+A_e^\dagger U_e^\dagger .      \tag{GD04}
\]

This operator is Hermitian, conserves GA total charge, and commutes with the
pair total momentum.  Its coefficient is the inherited `h_e`; no rescue
interaction is added.

The source-off autonomous class consists of the lifted F3 Hamiltonian, the
mechanical kinetic terms, and only relative-coordinate exchanges or
potentials:

\[
 H_{\rm GD}[j]
 =\iota_{\rm rec}(H_{\rm F3}[j])
 +\sum_a{P_a^2\over2M_a}
 +H_{\rm ctrl}+H_\partial
 +\Pi_{C,{\rm on}}H_{R\partial}.                  \tag{GD05}
\]

Here `iota_rec` is the ordered, termwise encoded substitution
`Z_e -> Z_e tensor I`, `X_e -> mathsf X_e`.  It is an exact representation
of the inherited algebra **after restriction to the image of `W_E`**; no
global Pauli-algebra identity is claimed on charge sectors outside that
image.  `H_ctrl` and `H_partial` own the admitted controller and exterior
boundary factors.  Every active term in `H_(R partial)` must exchange charge
and recoil with an explicit
boundary factor through a relative Weyl operator of the form (GD02).  During
the GA response hold the controller is in an invariant `Pi_(C,off)` sector,
so the outer port is off.  The one fixed Hamiltonian is autonomous; the hold
is a declared invariant controller sector, not a time-dependent c-number
switch.  Explicitly, `Pi_(C,on)Pi_(C,off)=0`,
`[H_ctrl,Pi_(C,off)]=0`, and `H_(R partial)` acts trivially on the controller,
so `[Pi_(C,on),H_(R partial)]=0`; the displayed port term is Hermitian and the
off sector invariant.

All coefficients and pair potentials in this admitted class depend only on
relative coordinates.  Hence simultaneous translation of every included
factor is an exact symmetry.  This defines the conditional controller/port
class; it does not construct a nontrivial controller transition or energy-
current packet.  A classical clamp or omitted exterior is not in the class,
because it would be an unowned momentum port.

## 3. Exact charge-recoil encoder

For clarity take equal recoil masses `M_e` on `L_e,R_e`; the unequal-mass
condition is given below.  Choose `p_e/hbar in Lambda^*`.  Since
`kappa_e in 2 Lambda^*`, every one-factor momentum displayed below belongs to
the torus spectrum.  In the fixed pair center-momentum sector `2p_e`, put

\[
\begin{aligned}
 |\chi_{e,-}\rangle
 &=|p_e+\tfrac12\hbar\kappa_e\rangle_L
   |p_e-\tfrac12\hbar\kappa_e\rangle_R,\\
 |\chi_{e,+}\rangle
 &=|p_e-\tfrac12\hbar\kappa_e\rangle_L
   |p_e+\tfrac12\hbar\kappa_e\rangle_R .          \tag{GD06}
\end{aligned}
\]

Extend GA's isometry to

\[
\begin{aligned}
 W_e|-\rangle&=|-\rangle_L|+\rangle_R|\chi_{e,-}\rangle,\\
 W_e|+\rangle&=|+\rangle_L|-\rangle_R|\chi_{e,+}\rangle .  \tag{GD07}
\end{aligned}
\]

Because `U_e|chi_(e,-)>=|chi_(e,+)>`, direct multiplication gives

\[
 \boxed{
 \mathsf X_eW_e=W_eX_e,
 \quad(q_*Z_e+Q_{R_e})W_e=0,
 \quad(P_{L_e}+P_{R_e})W_e=2p_eW_e.}               \tag{GD08}
\]

The recoil kinetic energy also restricts to one scalar on both codewords:

\[
 \left({P_{L_e}^2\over2M_e}+{P_{R_e}^2\over2M_e}\right)W_e
 =\left({|p_e|^2\over M_e}
       +{\hbar^2|\kappa_e|^2\over4M_e}\right)W_e. \tag{GD09}
\]

For unequal masses, write the two momenta as the mass-weighted center piece
plus relative momentum `q`.  A kick sends `q->q-hbar kappa_e`; equality of
the two kinetic energies requires
`q dot (hbar kappa_e)=hbar^2|kappa_e|^2/2`.  The symmetric half-kick above is
the minimal equal-mass realization.

Tensoring (GD07) over all links gives `W_E`.  Since (GD08)--(GD09) hold for
each link separately, `W_E` covers the complete inherited link Hilbert,
including every ice `P` and off-ice/virtual `Q` configuration.  It is not a
ground-subspace-only encoding.

### Theorem `TORP-1` -- exact full-code translation-owning lift

On any finite F3 link set and fixed mechanical center-momentum sectors,
`W_E` intertwines the complete inherited `Z,X` algebra with the charge- and
momentum-conserving recoil algebra.  The only source-off mechanical addition
on the full code is

\[
 \epsilon_{\rm rec}I
 =\sum_e\left({|p_e|^2\over M_e}
       +{\hbar^2|\kappa_e|^2\over4M_e}\right)I.    \tag{GD10}
\]

Thus a nonzero local recoil can coexist exactly with the inherited link
dynamics; it need not change the FY Hamiltonian.

## 4. Exact factor-edge momentum and charge ledger

Let `Y_e=A_eU_e` and `H_(flip,e)=-h_e(Y_e+Y_e^dagger)`.  Equations
(GD03)--(GD04) give

\[
 J^P_{L_e\to R_e}
 =ih_e\kappa_e(Y_e-Y_e^\dagger),                  \tag{GD11}
\]

and the exact operator identities

\[
 \boxed{
 \dot P_{L_e}+J^P_{L_e\to R_e}=0,
 \qquad
 \dot P_{R_e}-J^P_{L_e\to R_e}=0.}               \tag{GD12}
\]

The GA charge current on the same exchange is

\[
 I^Q_{L_e\to R_e}
 =-{2ih_eq_*\over\hbar}(Y_e-Y_e^\dagger),
 \qquad
 J^P_{L_e\to R_e}
 =-{\hbar\kappa_e\over2q_*}I^Q_{L_e\to R_e}.     \tag{GD13}
\]

The second equality is a property of this minimal recoil completion.  Its
free vector `kappa_e` is not a derived universal charge-to-momentum law and
must not be called `T^{0j}`.

More generally, every active factor-pair term

\[
 H_{ab}=-g_{ab}(B_{ab}U_{ab}+B_{ab}^\dagger U_{ab}^\dagger),
 \qquad U_{ab}=e^{i\kappa_{ab}\cdot(r_b-r_a)},      \tag{GD14}
\]

where the internal transfer is mechanical-translation neutral,
`[P_a^i,B_(ab)]=[P_b^i,B_(ab)]=0`, has one antisymmetric momentum current with
`dot P_a+J_(a->b)=0`, `dot P_b-J_(a->b)=0`.  A relative potential contributes
its ordinary equal-and-opposite gradient force to the same ledger.  If
`B_(ab)` itself carries mechanical coordinate or momentum dependence, its
additional commutators must be retained and (GD14)'s bare Weyl-current formula
does not apply.  Summing
over the complete support/reservoir/controller/boundary factor graph yields

\[
 \boxed{{d\over dt}\sum_aP_a=0}                   \tag{GD15}
\]

with every internal factor edge counted once.  This is the auxiliary
factor-edge translation and recoil ownership absent in GB's fixed-support
parent.  Because `kappa_e` and the exact momentum-sector states are not yet
bound to GC's diamond positions, (GD11) is not a physical diamond-space
current.

### Theorem `TORP-2` -- exact closed factor-edge momentum ledger

The declared autonomous relative-coordinate class owns every impulse it
contains.  Conditionally, an active outer port is in this class only when its
exterior charge and recoil factor appear in (GD15).  A reservoir-only outer operation still
leaves the GA code, exactly as GA proved; momentum closure does not make it a
code-preserving hold operation.

## 5. Preservation of FV/FY on the encoded hold

Freeze the **`GD-TRANSLATION-OWNING-FULL-CODE-HOLD`**:

1. every link/reservoir pair is prepared by (GD07) in a fixed center-momentum
   sector;
2. the controller is in an invariant port-off sector during response;
3. boundary and controller momenta are in fixed kinetic-energy sectors;
4. `kappa_e`, the recoil masses, the isometry `W_E`, and all added mechanical
   terms are independent of the frozen spatial query `j_ij`, except for at
   most one declared common full-code identity/reference source; and
5. every other support/reservoir/controller/boundary term restricts to one
   common scalar on the full `P+Q` code, with first `j_ij` derivative absent
   or the same declared identity.

Let `W_hold` denote `W_E` tensored with the declared controller-off and fixed
boundary/controller momentum states.  Then for one scalar `c_GD[j]`,
including (GD10) and the fixed controller/boundary kinetic reference,

\[
 H_{\rm GD}[j]W_{\rm hold}
 =W_{\rm hold}(H_{\rm F3}[j]+c_{\rm GD}[j]I).     \tag{GD16}
\]

Consequently the consistent Feshbach reference shift gives

\[
 H_{{\rm eff},\rm GD}[j;z+c_{\rm GD}[j]]W_{P,\rm hold}
 =W_{P,\rm hold}(H_{{\rm eff},\rm F3}[j;z]
                  +c_{\rm GD}[j]I_P).             \tag{GD17}
\]

Every nonidentity FV/FY coefficient, virtual gap, H2/H4/H6 fold and ring
matrix element, commutator, rank, native phase, and connected response is
unchanged.  A source-dependent kick, mass, support potential, controller, or
port is outside this hold and requires a new source/rank calculation.  In
particular, source independence is not enough if an added term restricts to
a nonidentity operator, preserving GA's repaired hold condition.

### Theorem `TORP-3` -- exact B1 reduction

The minimal recoil parent passes B1's **algebraic existence** condition on the
declared hold: one autonomous translation-invariant Hamiltonian owns the
nonzero flip recoil and every included boundary impulse, yet reduces on the
full encoded hold to GA/FV/FY nonidentity physics modulo one common reference
identity.  This is existence on an auxiliary common mechanical torus, not a
claim that the factor edges have already been placed in physical diamond
space.

This is not B2.  The inherited `j_ij` remains a frozen source query; this lane
does not derive `T^{00}`, `T^{0i}`, `T^{ij}`, source work, contacts, or the
native divergence from one spacetime source family.

## 6. Sharp minimality and physical-solder boundary

Suppose instead that the encoder factorized as logical code times the **same**
normalized mechanical state `|chi>` for both GA codewords while
`kappa_e!=0`.  Exact restriction to `X_e` would
require `U_e|chi>=e^{i phi}|chi>`.  In the momentum basis, `U_e` is a bilateral
shift along an infinite nonzero reciprocal-lattice orbit.  Its coefficients
would obey `|c_(p+n hbar kappa)|=|c_p|` for every integer `n`; square
summability then forces every coefficient on the orbit to vanish.  Hence no
nonzero normalizable eigenvector exists.  If `kappa_e=0`, the same-state hold is
possible but (GD11) is zero and no recoil is owned.

### Theorem `TORP-4` -- code-independent same-mechanical-state hold no-go

Within the minimal torus/Weyl class, a nonzero exact recoil is incompatible
with a code-independent logical-times-one-mechanical-state hold.  This does
**not** exclude link/reservoir product momentum states: each `chi_(e,+/-)` in
(GD06) is such a product.  It proves only that the two logical codewords
cannot reuse the same mechanical state.  The half-kick correlation in
(GD06)--(GD07), or an equivalent code-dependent recoil correlation, is
necessary rather than decorative machinery.

Three physical quantities remain deliberately unselected:

- the magnitude/direction of `kappa_e` and its binding to the GC `A3`
  support geometry;
- a normalizable localized-support state or controlled wavepacket limit,
  rather than the exact momentum-sector witness used here; and
- the position/source dependence that turns inherited `j_ij` into one
  autonomously owned stress insertion with controller work and contacts.

Those are the bounded input to B2.  The exact momentum code proves B1
existence; it does not manufacture the missing stress Ward packet.

## 7. Disposition

`EXACT_MINIMAL_RELATIVE_WEYL_RECOIL_DRESSING_OF_THE_EXISTING_GA_CHARGE_FLIP__EXACT_NORMALIZABLE_HALF_KICK_CORRELATED_CODE__DRESSED_FLIP_INTERTWINES_WITH_X_ON_THE_COMPLETE_P_PLUS_Q_LINK_HILBERT__RECOIL_KINETIC_ENERGY_IS_ONE_FULL_CODE_SCALAR__EVERY_ACTIVE_RELATIVE_EXCHANGE_OWNS_EQUAL_AND_OPPOSITE_MOMENTUM__EXPLICIT_BOUNDARY_FACTOR_CLOSES_OUTER_LEDGER_BUT_ACTIVE_RESERVOIR_ONLY_PORT_BREAKS_HOLD__FV_FY_NONIDENTITY_SOURCE_AND_H6_RESPONSE_PRESERVED_MODULO_ONE_COMMON_REFERENCE_SHIFT__NONZERO_RECOIL_CODE_INDEPENDENT_SAME_MECHANICAL_STATE_HOLD_NO_GO__KICK_SCALE_GEOMETRIC_SOLDER_SPACETIME_SOURCE_CONTACTS_STRESS_WARD_TENSOR_CONE_GRAVITY_AND_G_OPEN`

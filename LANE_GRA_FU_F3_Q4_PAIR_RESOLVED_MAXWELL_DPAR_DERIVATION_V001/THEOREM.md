# Pair-resolved Maxwell realization of the q4 degree-pair affine response

**Lane ID:** `GRA-FU-F3-Q4-PRMDD-V001`

**Short name:** `PRMDD`

**Date:** 2026-08-27

**Claim class:** exact tetrahedral lumped-capacitance first-variation
obstruction; exact four-terminal elastance matching theorem; exact central-
kernel derivation of `DPAR`; exact ideal-Coulomb slope and normalization;
exact conditional composition with the prior FT microscopic rank theorem;
explicit same-parent, port, cancellation, and circularity boundary

**Status:**
`LUMPED_TETRAHEDRAL_NODE_CHARGING_IS_A1_ONLY_AND_E_NULL__PAIR_RESOLVED_FOUR_TERMINAL_ELASTANCE_REPRODUCES_THE_UD_DEGREE_PAIR_OPERATOR_AT_SOURCE_OFF__CENTRAL_MUTUAL_KERNEL_DERIVES_DPAR__IDEAL_FIXED_COUPLING_COULOMB_GIVES_LAMBDA_MINUS_ONE_HALF_AND_CONDITIONAL_UD_ALPHA_LENGTH_NORMALIZATION__FT_MICROSCOPIC_RANK6_FOLLOWS_ONLY_AFTER_COMPLETE_PHYSICAL_SOLDER_AND_NONCANCELLATION__UD_REMAINS_A_GAUSS_CHARGE_PENALTY_NOT_IR_MAXWELL_STIFFNESS__VISIBLE_U1_AND_GRAVITY_NOT_DERIVED`

**Not claimed:** that the current reduced F3 parent already contains a
capacitance matrix, physical charge, physical length, or Maxwell Green kernel;
that a bare inherited `X_a` is a gauge-invariant physical-charge flip;
that the q4 emergent compact `U(1)` is visible electromagnetism; that `U_d` is
the transverse infrared Maxwell electric stiffness; that the ideal Coulomb
kernel is the complete QED potential; that alpha is calculated; that a
projected or CTP source has rank six; a Ward packet, tensor pole, RGRL-B,
gravity, or `G`.

## 1. Exact question and dependency custody

`FT` proved that the q4 degree square already contains six physical pair
operators and that a pair-separation law called `DPAR` would supply the two
linear `E` strain directions missing from the additive edge source.  It also
proved that the source-off Hamiltonian cannot choose its own derivative:
`DPAR` was neither inherited nor adopted.

This lane asks the next narrower physics question:

> Can an explicit electromagnetic or capacitive realization of the existing
> degree penalty derive the `DPAR` derivative, and which apparently natural
> realization provably cannot?

The load-bearing dependency bytes are frozen in `DEPENDENCIES.sha256`.  Their
roles and ceilings are:

| role | dependency | inherited boundary |
|---|---|---|
| microscopic incidence parent | `LANE_GRA_BS_F3_QIRN_MICRO_ACTION_V001/MICRO_ACTION.md` | `U_d(d-d_*)^2` is a parametric, distance-free degree penalty; the current and port terms contain no frozen mutual-capacitance law |
| q4 coframe and physical interface | `LANE_GRA_FE_F3_Q4_DIAMOND_ICE_CARRIER_JOIN_V001/THEOREM.md` | tetrahedral roots are exact, while physical coexistence, absolute length, and visible charge remain open |
| Gauss interpretation | `LANE_GRA_CL_F3_DEGREE_LOCK_GAUGE_PHASE_BRIDGE_V001/THEOREM.md` | `U_d` exactly penalizes discrete divergence charge |
| stiffness separation | `LANE_CROSS_ALPHA_GRA_F3_U1_STIFFNESS_IDENTIFIABILITY_V001/THEOREM.md` | `U_d` is not thereby transverse Maxwell electric stiffness |
| Coulomb-side and sixth-order controls | `LANE_CROSS_ALPHA_GRA_F3_COULOMB_TANGENT_V001/THEOREM.md`; `LANE_CROSS_ALPHA_GRA_F3_DIAMOND_SIXTH_ORDER_V001/THEOREM.md` | supplied-support spin-one phase comparisons do not identify visible EM or physical geometry |
| existing carrier current | `LANE_GRA_CS_F3_CARRIER_RESPONSE_SUPPORT_SELECTION_V001/THEOREM.md` | the inherited `J^2` term is edge-local and supplies no sibling mutual-current geometry |
| visible-sector boundaries | `LANE_GRA_CA_F3_UNIFORM_RESOURCE_MATCHING_LOCALITY_OBSTRUCTION_V001/JOINED_MULTIFORCE_PARENT.md`; `LANE_RFT_ALPHA_GAUGE_NORMALIZATION_V001/THEOREM.md`; `EU`; `EY` | alpha/Maxwell use needs a fixed physical sector and same-front solder; fixed-background or induced-curvature endpoints cannot be imported to manufacture the pregeometric source |
| conditional strain theorem | `LANE_GRA_FT_F3_Q4_DEGREE_PAIR_GEOMETRIC_STRAIN_BOUNDARY_V001/THEOREM.md` and hostile audit | `DPAR` with nonzero slope gives exact microscopic rank six, but only rank-three `A1+E` after direct ice restriction; CTP and gravity remain open |

The present completion is therefore **prospective**.  It may physically
realize the already displayed Gauss-charge penalty without identifying that
penalty with the distinct infrared Maxwell stiffness.

## 2. Exact q4 algebra

Use the regular tetrahedral unit vectors

\[
 n_1={1\over\sqrt3}(1,1,1),\quad
 n_2={1\over\sqrt3}(1,-1,-1),\quad
 n_3={1\over\sqrt3}(-1,1,-1),\quad
 n_4={1\over\sqrt3}(-1,-1,1).                    \tag{FU01}
\]

Let the four commuting link involutions be `Z_a`, put

\[
 n_a^{\rm occ}={1-Z_a\over2},\qquad
 d=\sum_an_a^{\rm occ},\qquad P_{ab}=Z_aZ_b.       \tag{FU02}
\]

Then

\[
 d-2=-{1\over2}\sum_aZ_a,
\qquad
 \boxed{U_d(d-2)^2=U_dI+{U_d\over2}\sum_{a<b}P_{ab}.} \tag{FU03}
\]

For a declared positive scale `a_*`, define sibling terminal positions and
roots

\[
 x_a=a_*n_a,\qquad r_{ab}=x_b-x_a,\qquad
 r_0^2=|r_{ab}|^2={8a_*^2\over3},\qquad
 \widehat R_{ab}={r_{ab}r_{ab}^{\mathsf T}\over r_0^2}. \tag{FU04}
\]

Equation (FU04) is initially a relational coframe and scale declaration.  It
is not an absolute physical length until assumption `S2` below is discharged.

## 3. Lumped tetrahedral charging cannot supply `E`

Consider the strongest symmetric single-node charging realization,

\[
 H_{\rm lump}[j]=A[j](d-2)^2+B[j]I,                \tag{FU05}
\]

where, for example, `A[j]=q_L^2/(2C_\Sigma[j])`, and where `A` and `B` are
real differentiable scalar functions at `j=0`.  Tetrahedral covariance means

\[
 A[O_\pi jO_\pi^{\mathsf T}]=A[j],\qquad
 B[O_\pi jO_\pi^{\mathsf T}]=B[j]                 \tag{FU06}
\]

for every orthogonal tetrahedral label permutation `O_pi`.

### Lemma `PRMDD-L1` -- unique invariant linear functional

Every tetrahedrally invariant linear functional on real symmetric `3 x 3`
matrices is proportional to `tr j`.

#### Proof

Identify a linear functional with one symmetric gradient tensor `T` under the
Frobenius pairing.  Equation (FU06) requires
`O_pi^T T O_pi=T` for all tetrahedral permutations.  The natural three-
dimensional tetrahedral representation is irreducible, so its real symmetric
commutant is the scalar line.  Equivalently, direct use of the sign vectors in
(FU01) first kills every off-diagonal entry of `T` and then equates its three
diagonal entries.  Hence `T=cI`, and the functional is `c tr j`. QED.

It follows that

\[
 D_jH_{\rm lump}|_0
 =c_A(\operatorname{tr}j)(d-2)^2+c_B(\operatorname{tr}j)I. \tag{FU07}
\]

Both diagonal-traceless `E` strains therefore annihilate the entire first
variation.  On the ice fiber `d=2`, even the non-reference term vanishes
identically for every value of the scalar capacitance.

### Theorem `PRMDD-1` -- lumped-capacitance obstruction

A tetrahedrally symmetric single total-node capacitance, or a shared-node law
proportional to `(sum_a I_a)^2`, has an `A1`-only first variation and exact
`E` nullity two.  Geometry dependence of its one scalar coefficient cannot
close the FT strain-source gap.  This is a no-go for the lumped realization,
not for a resolved mutual-capacitance matrix.

The inherited BS11 current square does not evade this result.  It is a sum of
edge-local `n_e(J^\psi_{uv})^2` terms; its exact reduction contains endpoint
occupation products, not sibling signed currents `I_aI_b`, and neither the
current amplitude nor its coefficient is a function of (FU04).

## 4. Pair-resolved four-terminal elastance

Freeze four coexisting terminal modes and a complete physical reference or
ground.  Let `C(F)` be the resulting real, nonsingular four-terminal
capacitance matrix and let

\[
 {\cal E}(F)=C(F)^{-1}                              \tag{FU08}
\]

on the full four-terminal charge domain.  This full inverse is load bearing
for the operator equality (FU14) on all sixteen q4 states and for the
unprojected FT rank composition (FU29).

If an ungrounded capacitance matrix has a common-potential zero mode, its
inverse on the neutral quotient may instead be used **only after restriction
to** `sum_a Z_a=0`, in particular on the six ice states.  That quotient does
not define the off-ice defect energies, cannot establish (FU14) as a full
operator identity, and cannot discharge (FU29).  Silently replacing the full
inverse by a quotient inverse is forbidden.

Solder the link involutions prospectively to signed terminal charges

\[
 q_a=q_*Z_a,\qquad q_*\ne0,                        \tag{FU09}
\]

and define the resolved electrostatic energy

This charge solder changes the status of the inherited link flip.  With
`X_a=\sigma_a^+ + \sigma_a^-` and
`[Z_a,\sigma_a^\pm]=\pm2\sigma_a^\pm`, one has

\[
 [q_*Z_a,X_a]=2q_*(\sigma_a^+-\sigma_a^-)\ne0.    \tag{FU09a}
\]

Thus a bare inherited `X_a` cannot be promoted to visible charged dynamics.
A charge-conserving completion must own a reference/reservoir charge `Q_R`
and transfer operators `T_{a,-}`, `T_{a,+}=T_{a,-}^dagger` with

\[
 [Q_R,T_{a,\pm}]=\pm2q_*T_{a,\pm},\qquad
 \widetilde X_a=\sigma_a^+T_{a,-}+\sigma_a^-T_{a,+},
 \qquad[Q_{\rm tot},\widetilde X_a]=0,             \tag{FU09b}
\]

where `Q_tot=q_* sum_a Z_a+Q_R`.  Exact composition with the inherited
action additionally requires a fixed-total-charge encoded subspace on which
`\widetilde X_a` is unitarily equivalent to the inherited `X_a`, with every
reservoir state, transfer current, source-work term, and port retained in the
ledger.  Otherwise the charged completion is a larger dynamical parent and
requires a fresh source and rank audit.  If `Z_a` is kept only as an internal
pseudocharge, neither the visible-`U(1)` nor the alpha normalization below is
licensed.

Now define the resolved electrostatic energy

\[
 H_C(F)={q_*^2\over2}Z^{\mathsf T}{\cal E}(F)Z
       +E_{\rm ref}(F)I.                           \tag{FU10}
\]

Because `Z_a^2=I`, the diagonal elastances contribute only an identity:

\[
 H_C(F)=
 \left[{q_*^2\over2}\sum_a{\cal E}_{aa}(F)+E_{\rm ref}(F)\right]I
 +q_*^2\sum_{a<b}{\cal E}_{ab}(F)P_{ab}.           \tag{FU11}
\]

At the regular tetrahedral point, exact `S4` covariance gives

\[
 {\cal E}_{aa}(I)=e_s,
 \qquad {\cal E}_{ab}(I)=e_m\quad(a\ne b).         \tag{FU12}
\]

The full grounded matrix is positive definite precisely when
`e_s-e_m>0` on the three contrast modes and `e_s+3e_m>0` on the common mode.
Positive definiteness can therefore coexist with the required positive
mutual coefficient.

Impose the source-off matching conditions

\[
 \boxed{q_*^2e_m={U_d\over2}},\qquad
 E_{\rm ref}(I)=U_d-2q_*^2e_s.                    \tag{FU13}
\]

Substitution in (FU11) gives the exact operator equality

\[
 \boxed{H_C(I)=U_dI+{U_d\over2}\sum_{a<b}P_{ab}
                   =U_d(d-2)^2.}                 \tag{FU14}
\]

The second equation in (FU13) is an energy-reference assignment; only the
first fixes a nonidentity observable coefficient.  A different diagonal
self-energy changes the reference, not the six pair operators.

### Theorem `PRMDD-2` -- exact source-off elastance realization

A tetrahedrally symmetric, fully grounded pair-resolved four-terminal
elastance with (FU09) and (FU13) realizes the existing q4 degree square exactly
at source off on the full q4 Hilbert space.  It adds no new incidence-sector
source-off operator.  It is nevertheless a new **physical completion** of the
previously parametric F3 action because F3 did not derive the terminal charge
solder, the reference conductor, `C(F)`, or the charge-conserving dynamical
lift (FU09b).  A neutral-quotient construction proves only the ice-restricted
version and is not PRMDD-2.

## 5. A central mutual kernel derives `DPAR`

Assume now that the complete off-diagonal elastance in the scored local block
has one real differentiable central kernel,

\[
 {\cal E}_{ab}(F)=V(|Fr_{ab}|),\qquad a\ne b,
 \qquad V(r_0)\ne0.                                \tag{FU15}
\]

Define

\[
 x_{ab}(F)={|Fr_{ab}|^2\over r_0^2},\qquad
 g(x)={V(r_0\sqrt{x})\over V(r_0)}.                \tag{FU16}
\]

Then `g(1)=1`.  For a real symmetric strain source with

\[
 F(j)=I-{j\over2}+O(j^2),                          \tag{FU17}
\]

one has

\[
 x_{ab}(F(j))=1-j:\widehat R_{ab}+O(j^2).          \tag{FU18}
\]

The ordinary chain rule gives

\[
 \boxed{\lambda:=g'(1)={r_0V'(r_0)\over2V(r_0)}.} \tag{FU19}
\]

Using (FU13), the nonidentity part of (FU10) becomes

\[
 {U_d\over2}\sum_{a<b}
 g\!\left({|Fr_{ab}|^2\over r_0^2}\right)P_{ab}, \tag{FU20}
\]

which is exactly the `DPAR` normal form of FT.  Differentiating the physical
family, rather than assigning a tensor after projection, yields

\[
 H_C[j]=H_C[0]-{U_d\lambda\over2}
 j_{ij}\sum_{a<b}\widehat R_{ab}^{ij}P_{ab}+O(j^2)
 +H_{\rm id}[j],                                   \tag{FU21}
\]

and under FT's convention `Q^{ij}=-2 partial H/partial j_ij|_0`,

\[
 \boxed{Q_{\rm pair}^{ij}=U_d\lambda
 \sum_{a<b}\widehat R_{ab}^{ij}P_{ab}.}           \tag{FU22}
\]

Real `V`, real `F`, and symmetric real elastance preserve Hermiticity.  The
common central law is `S4` covariant because label permutations orthogonally
permute the six roots.

### Theorem `PRMDD-3` -- physical derivation of the DPAR component

Under the pair-resolved solder (FU08)--(FU15), the affine derivative of the
same source-off degree energy is exactly `DPAR`, with slope (FU19).  The two
missing `E` directions are present iff the **complete net** `E` slope remains
nonzero after the terms in Section 7 are included.  A stationary kernel with
`V'(r_0)=0` supplies only a quadratic contact at leading order and fails the
linear FT gate.

## 6. Ideal Coulomb specialization and alpha normalization

For a homogeneous three-dimensional vacuum or medium with a coupling held
fixed at one declared matching scale, take the ideal mutual elastance

\[
 V_C(r)={1\over4\pi\epsilon_0\epsilon_r r}.        \tag{FU23}
\]

Equations (FU16) and (FU19) give

\[
 \boxed{g_C(x)=x^{-1/2},\qquad \lambda_C=-{1\over2}.} \tag{FU24}
\]

Let `q_*=kappa e(mu_0)` and use the rationalized SI convention

\[
 \alpha(\mu_0)={e(\mu_0)^2\over4\pi\epsilon_0\hbar c}. \tag{FU25}
\]

The nonidentity source-off matching condition in (FU13) is then

\[
 {U_d\over2}={q_*^2\over4\pi\epsilon_0\epsilon_r r_0}
 ={\kappa^2\alpha(\mu_0)\hbar c\over\epsilon_r r_0},
\]

or

\[
 \boxed{U_d={2\kappa^2\alpha(\mu_0)\hbar c
                    \over\epsilon_r r_0},
 \qquad r_0=a_*\sqrt{8/3}.}                       \tag{FU26}
\]

In rationalized natural units `hbar=c=epsilon_0=1`, this is the identical
statement `e^2=4pi alpha` and `U_d=2 kappa^2 alpha/(epsilon_r r_0)`.

Equation (FU26) is **not inherited F3** and is not a calculation of alpha.  It
is the coefficient match inside this prospective physical completion.  It
uses an independently normalized charge, physical length, dielectric and
Maxwell kernel.  Conversely, solving it for alpha without independently
owning those quantities would be circular parameter fitting.

The fixed-coupling qualifier is load bearing.  For a complete radial kernel
`V(r)=A(r)/r`,

\[
 \lambda=-{1\over2}+{r_0\over2}{A'(r_0)\over A(r_0)}. \tag{FU27}
\]

If `A(r)` is represented by a running `alpha(mu)` with `mu proportional 1/r`
and `beta_alpha=mu d alpha/d mu`, the extra term is
`-beta_alpha/(2alpha)`.  Dielectric dispersion, screening, conductor shape,
and radiative corrections belong in the complete `V`; they may not be hidden
while retaining the exact `-1/2` label.

### Theorem `PRMDD-4` -- ideal-Coulomb conditional result

The ideal fixed-coupling Coulomb completion derives the nonzero DPAR slope
`-1/2` and the normalization (FU26).  This realizes a possible microscopic
origin of the **Gauss-charge penalty**.  It does not identify `U_d` with
`U_E^IR`, `K_B^IR`, the bare ring coefficient, photon velocity, or the alpha
invariant extracted from matched infrared stiffnesses.  It is fully
compatible with the U1SI, Coulomb-tangent, and sixth-order proof ceilings.

## 7. Complete-source, boundary, common-mode, and cancellation conditions

A local point-pair formula is not automatically a complete physical action.
Write the full source derivative in the local pair sector as

\[
 D_jK_{ab}|_0=-{U_d\lambda\over2}
 (j:\widehat R_{ab})+D_jK_{ab}^{\rm rem}|_0.       \tag{FU28}
\]

All of the following are load bearing.

1. **Diagonal/self terms.**  Since `Z_a^2=I`, they are identity operators.
   Their source derivatives and work still remain in the complete ledger, but
   they cannot supply or cancel a nonidentity pair `E` operator.
2. **Lumped/common mode.**  A term `c[j](sum_a Z_a)^2` is the degree square
   itself.  Under exact tetrahedral covariance its independent coefficient
   derivative is `A1` only by PRMDD-1.  It cannot cancel the centered `E`
   response.
3. **Other local mutual kernels.**  Under exact `S4` covariance their `E` part
   is another scalar multiple of the unique `E` intertwiner.  Define
   `lambda_E^net=lambda+lambda_E^rem`.  Exact rank closure requires
   `lambda_E^net != 0`; equality is a real physical cancellation and fails the
   theorem.  Exact identification of the *whole* pair source with (FU20)
   additionally requires the remainder either vanish or share the same common
   normalized radial law.  An extra `A1` identity or common-mode response may
   coexist but must not be double counted as a second degree energy.
4. **Broken tetrahedral boundaries.**  A support, ground, dielectric, reader,
   or controller that breaks `S4` produces a general six-by-six derivative.
   Then one must rank the measured full map; a common `lambda` may not be
   asserted from the central subterm alone.
5. **Nonlocal/cross-node terms.**  Coulomb or circuit elimination generally
   produces `Z_{v,a}Z_{w,b}` operators as well.  On the full product Hilbert
   space these Walsh strings are independent of the local `P_ab`.  They cannot
   be dropped or absorbed into (FU03).  For exact composition with FT they
   must be shielded to zero, reduce to declared identity/reference terms, or
   already belong to the prospectively frozen larger parent whose full source
   is re-ranked.
6. **Capacitance gauge and ports.**  The full theorem requires a physical
   ground/reference and a nonsingular four-terminal matrix.  A neutral
   quotient is admissible only for the explicitly ice-restricted result and
   cannot support the off-ice defect spectrum or unprojected rank.  Source
   work, dielectric, recoil/support, boundary fields and any dissipative
   environment must be owned once.
7. **Charge-conserving dynamics.**  If (FU09) is a physical charge solder,
   every inherited flip must be replaced by or exactly encoded as the dressed
   transfer (FU09b).  The compensating reservoir/current, its work, and its
   ports are part of the same parent.  Bare `X_a` is not gauge-invariant
   charged dynamics.  If the dressed parent is not exactly equivalent on a
   frozen encoded subspace, its complete source must be re-ranked.
8. **Order of reductions.**  The complete `F -> H_full(F)` family is frozen
   first.  Any exact electrostatic elimination is performed with its source
   dependence retained, and only then is the fixed incidence Feshbach map
   applied.  Assigning (FU22) after inspecting a projected response is
   forbidden.

These conditions distinguish a derivation from a suggestive circuit analogy.

## 8. Conditional composition with FT

Freeze the following physical-solder packet:

- **`S1 — coexistence`:** `Q4-CARRIER/EDGE-LIFT` or an equivalent complete
  construction makes the four q4 sibling factors simultaneous physical modes;
- **`S2 — length`:** `r_0` and `F` are operational relational shape/length
  variables earned independently of the sought gravitational endpoint;
- **`S3 — charge and transfer`:** (FU09) maps the same link involutions to one
  canonically normalized physical charge sector at a declared scale, and
  (FU09b) supplies a conserved, gauge-invariant compensating transfer for
  every flip.  On the frozen scored subspace the dressed flips must be exactly
  equivalent to the inherited `X_a`; bare `X_a` is insufficient;
- **`S4 — complete grounded field`:** one nonsingular, port-complete grounded
  capacitance/Green-kernel parent owns every reference charge, compensating
  reservoir/current, boundary, common mode, source, support and work channel;
- **`S5 — source-off match`:** (FU13)--(FU14) hold without an omitted
  nonidentity source-off interaction on the full off-ice as well as ice q4
  state space; a neutral quotient does not discharge this premise;
- **`S6 — constitutive law`:** (FU15) or a fully ranked covariant generalization
  is derived before response scoring;
- **`S7 — noncancellation`:** the complete pair-sector `E` coefficient obeys
  `lambda_E^net != 0`;
- **`S8 — source ordering`:** the physical deformation is inserted before the
  fixed Feshbach reduction, with the FR/FS one-edge source left unchanged; and
- **`S9 — noncircular ancestry`:** charge, length and kernel are not inferred
  from the same emergent Maxwell, metric, Ricci, or gravity endpoint that the
  calculation is intended to explain.

Under `S1`--`S9`, PRMDD-3 derives the DPAR `E` component rather than adopting
it.  FT's exact theorem then applies:

\[
 \boxed{
 \operatorname{rank}D_jH_{\rm micro}|_{j=0}=6
 \quad\text{before Feshbach}.}                    \tag{FU29}
\]

The unchanged one-edge dyads supply `A1+T2`; the nonzero resolved pair
derivative supplies the missing `E`.  The source-off Hamiltonian in the scored
incidence sector remains (FU03).  Direct restriction of the pair map to local
ice still has only `A1+E` rank three, and the complete projected edge/Feshbach
rank is not upgraded by notation.

### Theorem `PRMDD-5` -- conditional microscopic rank closure

`S1`--`S9` plus the already sealed FT theorem imply exact microscopic
source-before-Feshbach rank six.  If any physical solder is absent, if the
full `E` derivative cancels, or if nonlocal source-off operators enlarge the
parent without a fresh rank audit, the implication does not fire.

## 9. Circularity ceiling and disposition

There are three different Maxwell statements, and this lane keeps them
separate:

1. a prospective electrostatic/capacitive **physical realization** of the
   microscopic Gauss-charge penalty `U_d`;
2. the conditional spin-one compact-`U(1)` Coulomb phase generated from the
   F3 degree lock and ring dynamics; and
3. the actual visible electromagnetic sector with its inherited running
   alpha and complete stress contribution.

Using statement 2 or a fixed-background version of statement 3 to create the
physical length, charge, compensating charge-transfer dynamics, or kernel
needed to derive statement 1 would close a logical circle.  EU's same-front
premises and circular-Maxwell control, EY's
open physical record-to-metric solder, and the CA post-emergence matching rule
all forbid that promotion.  EY is also explicit that alpha is not a second
stress source; (FU26) is a local coefficient match, not an extra gravitational
source.

The exact lawful result is therefore

\[
 \boxed{
 \begin{gathered}
 \text{tetrahedral lumped charging}\Longrightarrow A_1\text{ only and }E\text{ null},\\
 \text{pair-resolved central elastance}+\text{complete independent solder}
 \Longrightarrow \text{derived DPAR},\\
 \text{ideal fixed-coupling Coulomb}\Longrightarrow
 \lambda=-1/2\text{ and conditional (FU26)},\\
 \text{derived noncancelled DPAR}+\text{FT}
 \Longrightarrow \text{exact microscopic rank six},\\
 \text{no independent solder}\Longrightarrow
 \text{conditional physical completion only, not inherited F3}.
 \end{gathered}}                                    \tag{FU30}
\]

Here `complete independent solder` includes a full grounded elastance on the
off-ice charge sectors and a charge-conserving dressed realization of every
link flip.  A neutral quotient and bare `X_a` can establish neither the full
operator identity nor the unprojected rank claim.  This advances the
microscopic geometry-source problem without promoting a tensor pole, Einstein
endpoint, gravity, or `G`.

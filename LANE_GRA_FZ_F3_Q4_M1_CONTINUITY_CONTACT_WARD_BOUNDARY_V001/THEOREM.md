# F3/q4 nonzero-momentum continuity, contact, and Ward boundary

**Lane:** `LANE_GRA_FZ_F3_Q4_M1_CONTINUITY_CONTACT_WARD_BOUNDARY_V001`  
**Short name:** `CCWB`  
**Date:** 2026-08-28  
**Claim class:** exact finite projected-charge theorem; exact cyclotomic
supplied-embedding longitudinal diagnostic; exact rational TT quotient; sharp
current/contact dependency theorem

**Not claimed:** a continuum derivative, a complete electromagnetic Ward
identity, a diffeomorphism Ward identity, a protected tensor gauge mode, a
massless pole, RGRL-B, gravity, Newton's law, or `G`.

## 1. Frozen question

FY constructed the complete native-support spatial source through H6 at the
nonzero quotient momentum `m=1`.  It also found `k_i Q^{ij} != 0`, while
correctly declining to call that a Ward failure because no temporal density,
current, or contact term had been constructed.

This lane asks the next typed question without changing the Hamiltonian:

> Does the exact FY source already belong to a discrete charge-continuity and
> contact-complete Ward packet, and what transverse response follows from the
> source that actually exists?

The answer has three parts.  The projected incidence charge has an exact but
trivial continuity law.  The spatial source is exactly nontransverse under
FY's supplied continuum-embedding contraction and has no inherited
temporal/current/contact partner in the frozen parent.  Its kinematic TT
quotient is nevertheless exact and two-dimensional.  Because the physical
discrete divergence has not been derived, the embedding contraction is a
diagnostic, not the left side of a proved Ward identity.

The load-bearing parents are FU through FY.  In particular, FU already proved
that a physical charge interpretation requires the dressed transfer and
reservoir of FU09b; `FV-PURE` and FY froze a six-component spatial-strain
source after projection, not the missing temporal and contact packet.

## 2. Exact projected incidence-charge continuity

On the selected FO component let

\[
 G_v=d_v-2,
 \qquad P_{\rm ice}=\sum_{n\in\mathcal C_{180}}|n\rangle\langle n|.
 \tag{FZ01}
\]

Every one of the 180 basis states has degree two at every one of the sixty
vertices.  The H6 ring generator connects 420 undirected pairs of such states
and never leaves that component.  Therefore

\[
 \boxed{P_{\rm ice}G_vP_{\rm ice}=0,
 \qquad [H_6,P_{\rm ice}G_vP_{\rm ice}]=0
 \quad\text{for all }v.}
 \tag{FZ02}
\]

The projected incidence-charge continuity equation is consequently the exact
identity `0=0`.  It admits a zero-current representation, but it does not
select one: divergence-free or circulating projected currents could coexist.
No bond-current allocation is constructed in this parent.

### Theorem `CCWB-1` -- trivial-charge boundary

The FY spatial source cannot acquire its missing Ward partner from the
projected Gauss charge.  Equation (FZ02) is a genuine conservation result, but
it supplies no nonzero density and no derived current operator with which to
form the missing spatial stress balance.  A separately derived divergence-free
current would not change the density identity, but it is not furnished by
(FZ02).  Off-ice virtual charge sectors were used to derive the H2/H4/H6
coefficients; their physical reservoir currents and source contacts were not
retained as operators in the 180-state response space.

This statement concerns the incidence Gauss charge.  A canonically normalized
visible electromagnetic charge is even more demanding: FU09b requires a
reference/reservoir charge and dressed transfer before the inherited flip can
be interpreted as charged dynamics.

## 3. Exact supplied-embedding longitudinal diagnostic

The frozen shortest reciprocal representative is

\[
 q={1\over30}(1,5,-11),
 \qquad k={\pi\sqrt3\over60}(7,15,-17).
 \tag{FZ03}
\]

For tensor coordinates `(xx,yy,zz,2xy,2xz,2yz)`, twice the three
longitudinal covectors are

\[
 \begin{split}
 \ell_x^{(2)}&=(14,0,0,15,-17,0),\\
 \ell_y^{(2)}&=(0,30,0,7,0,-17),\\
 \ell_z^{(2)}&=(0,0,-34,0,7,15).
 \end{split}                                           \tag{FZ04}
\]

Reducing the native pair-source phase polynomial modulo `Phi_240` gives a
nonzero exact remainder.  One frozen witness is

\[
\begin{split}
 -14&+28z^2+14z^8+14z^{10}+14z^{24}-14z^{26}
 -28z^{34}\\
 &-14z^{42}-14z^{48}-14z^{56}+28z^{58}
 \not\equiv0\pmod{\Phi_{240}} .                    \tag{FZ05}
\end{split}
\]

Thus FY's supplied embedding contraction `k_i Q_pair^{ij}(m=1)` is exactly
nonzero over `Q(zeta_240)`; this is not a floating residual and is not yet a
derived discrete divergence.

FY proved

\[
 Q_{\rm complete}^{ij}(1;x)=
 \rho f_E(x)D_1^{ij}+R_1^{ij},
 \quad
 f_E=1-x^2-{37\over12}x^4-{16247\over900}x^6,       \tag{FZ06}
\]

where `D_1` is diagonal in the FO configuration basis and `R_1` is strictly
off diagonal.  At the two frozen samples,

\[
\begin{array}{c|c|c}
x&f_E(x)&\rho f_E(x)\\ \hline
2/5&2415673/3515625&2415673/113400\\
1/2&15853/57600&31706/14175
\end{array}                                          \tag{FZ07}
\]

so the nonzero diagonal embedding-longitudinal witness cannot be cancelled by
the ring source.  Independently, one complete 720+720 Hermitian ring entry
also has a nonzero exact embedding-longitudinal cyclotomic remainder.  Hence
the complete source remains nontransverse under this supplied contraction even
at a root of `f_E`; the ring entry cannot be cancelled by a diagonal term.

### Theorem `CCWB-2` -- supplied-embedding zero-slot diagnostic

The complete through-H6 `m=1` source fails `k_i Q^{ij}=0` exactly under FY's
supplied continuum-embedding contraction.  Therefore a spatial-only
completion with zero temporal/current and contact slots fails **under that
contraction**, and it will also fail for any later derived discrete divergence
proved proportional or equivalent to it on the frozen source space.

This does **not** rule out a complete Ward identity, and it does not yet rule
out zero temporal/contact slots for an unknown physical `Delta_m`.  Equation
(FZ03) is the supplied diamond embedding of the cyclic graph momentum.  A true
discrete continuity law must derive its bond divergence symbol from the same
physical current allocation; a future `Delta_m` need not equal `k_i` until
that equivalence is proved.  Replacing it silently by the continuum
contraction in (FZ03) is not licensed.

## 4. What a complete continuity/contact identity requires

Let `Delta_m` be the derived discrete divergence, let `P_m^j` be a temporal
momentum/current density, and put

\[
 L_m^j=(\Delta_m)_i Q_m^{ij}.
 \tag{FZ08}
\]

The operator identity required before a Ward claim is

\[
 \boxed{{i\over\hbar}[H,P_m^j]+L_m^j=0.}            \tag{FZ09}
\]

For a finite Hermitian `H`, an algebraic solution of (FZ09) exists iff

\[
 \boxed{\Pi_E L_m^j\Pi_E=0\quad\text{for every energy }E,} \tag{FZ10}
\]

where `Pi_E` is the complete projector onto the possibly degenerate energy
eigenspace.  This follows because the range of `ad_H` is the Hilbert--Schmidt
orthogonal complement of the block-diagonal commutant.  When (FZ10) holds,
one possible energy-basis solution is

\[
 (P_m^j)_{ab}={i\hbar(L_m^j)_{ab}\over E_a-E_b}
 \quad(E_a\ne E_b),                                \tag{FZ11}
\]

plus an arbitrary operator commuting with `H`.  Equation (FZ11) is an
algebraic inverse-Liouvillian construction.  It need not be local, descended
from a microscopic current, source-covariant, or port complete, and therefore
cannot by itself be promoted to physical continuity.

With

\[
 \chi^R_{AB}(t)=-{i\over\hbar}\theta(t)
 \langle[A(t),B]\rangle
\]

and the `e^{i omega t}` Fourier convention, (FZ09) implies

\[
 (\Delta_m)_i\chi^R_{Q^{ij}B}(\omega)
 -i\omega\chi^R_{P^jB}(\omega)
 =-{i\over\hbar}\langle[P_m^j,B]\rangle .          \tag{FZ12}
\]

The right side is the equal-time contact.  If `B` or `P` depends on the
external source, the corresponding second source derivative/seagull must be
added as well.  FY contains only `dH/dj_ij`; it contains neither the
`j_00/j_0i` derivatives nor the second derivatives needed to evaluate
(FZ12).

Two conservation types must not be conflated.  FU09b's dressed transfer is
the necessary microscopic `U(1)` charge-current completion.  Equations
(FZ08)--(FZ12), with `P^j=T^{0j}`, are the separate stress/translation Ward
packet.  The charge current contributes to the complete electromagnetic and
reservoir stress ledger, but charge conservation alone does not derive
energy--momentum conservation or its metric contacts.

### Theorem `CCWB-3` -- sharp dependency theorem

The frozen FY packet neither proves nor falsifies (FZ09)--(FZ12), because the
physical `Delta_m` and the objects on their temporal and contact sides are not
defined by that parent.  CCWB-2 prevents filling the missing slots with zero
only after adopting the supplied embedding contraction, or after proving the
physical divergence equivalent to it.  For an otherwise unknown physical
`Delta_m`, even the zero-current/contact question remains undecided.  Any
algebraic inverse-Liouvillian solution must remain a diagnostic until it is
derived from a local, port-complete source family.

## 5. Exact transverse response that does follow

Let `r=(7,15,-17)`, `r^2=563`, and

\[
 P_{ij}=\delta_{ij}-{r_ir_j\over563}.
 \tag{FZ13}
\]

For a symmetric spatial tensor `A`, define

\[
 \Lambda_{\rm TT}(A)=PAP-{1\over2}P\,\operatorname{tr}(PA).
 \tag{FZ14}
\]

Exact rational arithmetic gives

\[
 P^2=P,\quad Pr=0,\quad\operatorname{tr}P=2,
 \quad\Lambda_{\rm TT}^2=\Lambda_{\rm TT},
 \quad\operatorname{rank}\Lambda_{\rm TT}=2.       \tag{FZ15}
\]

Every image is exactly transverse and traceless.  FY's plus/cross sources
span this same two-dimensional kinematic image.  At `x=2/5` and `x=1/2`, the
complete source has TT ground-image rank two with positive Gram eigenvalues;
all responding finite-graph poles remain at positive gap.

### Theorem `CCWB-4` -- exact TT quotient, finite response only

The FY source has a mathematically exact two-component TT quotient and excites
both components at the two frozen samples.  This is the strongest transverse
response presently licensed.  It is a kinematic projection relative to the
supplied embedding direction, not a Ward-derived constraint, not a helicity
theorem, and not a massless graviton.

## 6. Minimal physical completion

No new record mechanism is required.  The minimum missing physics is the
already identified charge/current completion carried through the source
reduction:

1. instantiate FU09b's `Q_R`, `T_(a,+/-)`, and
   `X_tilde_a=sigma_a^+T_(a,-)+sigma_a^-T_(a,+)` so that
   `[Q_tot,X_tilde_a]=0`;
2. allocate the compensating reservoir and port currents on the same native
   A/B/link supports and derive the finite bond-divergence `Delta_m`;
3. freeze one complete source family
   `H[j_00,j_0i,j_ij; ports]`, so `T^00`, `T^0i`, `T^ij`, boundary exchange,
   and all second-derivative contacts come from the same Hamiltonian; and
4. perform Feshbach projection only after those source and contact terms are
   retained, then test (FZ09), (FZ10), and (FZ12).

The first missing current is therefore not arbitrary.  For the inherited
flip block `H_flip=-h X_tilde_a`, FU09b's commutators give exactly

\[
 \dot q_a={2ihq_*\over\hbar}
 (\sigma_a^+T_{a,-}-\sigma_a^-T_{a,+}),
 \qquad \dot Q_R=-\dot q_a .                       \tag{FZ16}
\]

Equation (FZ16) is the minimal conserved transfer-current seed.  What remains
is to locate it on the native physical supports, retain its source dependence
and ports through the six-step virtual histories, and derive rather than
assign the corresponding divergence and contacts.  A complete `T^{0j}` must
then be derived from the stress of that transfer, reservoir, field, support,
and ports; it is not identified with the scalar charge current itself.

This is a current/contact completion of the existing physical solder, not an
expansion of RFT.  If the completion changes the encoded source beyond
`FV-PURE`, the source rank and FY response must be recomputed rather than
assumed.

## 7. Disposition

`EXACT_PROJECTED_GAUSS_CHARGE_CONTINUITY_IS_TRIVIAL__EXACT_COMPLETE_M1_SOURCE_IS_LONGITUDINAL_UNDER_THE_SUPPLIED_EMBEDDING_CONTRACTION__ZERO_TEMPORAL_CONTACT_SLOTS_FAIL_ONLY_FOR_THAT_CONTRACTION_OR_A_PROVED_EQUIVALENT__PHYSICAL_DISCRETE_DIVERGENCE_AND_FULL_WARD_STATUS_REMAIN_UNDECIDED__EXACT_RATIONAL_EMBEDDING_TT_PROJECTOR_HAS_RANK_TWO_AND_BOTH_FINITE_RESPONSE_IMAGES_ARE_ACTIVE__MINIMAL_FU09B_CURRENT_PORT_SOURCE_COMPLETION_IS_IDENTIFIED__NO_NEW_RECORD_MACHINERY_WARD_GRAVITY_OR_G_PROMOTED`

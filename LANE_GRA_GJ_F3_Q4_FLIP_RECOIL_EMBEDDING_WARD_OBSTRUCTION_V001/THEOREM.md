# Flip-recoil supplied-embedding stress-Ward obstruction

**Lane:** `LANE_GRA_GJ_F3_Q4_FLIP_RECOIL_EMBEDDING_WARD_OBSTRUCTION_V001`  
**Short name:** `FRSEWO`  
**Date:** 2026-08-29  
**Plan gate:** `RF3a`  
**Claim class:** exact finite conditional no-go for the direct GD encoded
momentum density under FY/FZ's supplied embedding contraction; exact type
separation of momentum density and recoil flux; exact contact and boundary
ceiling

**Not claimed:** a native diamond-space divergence, a no-go for every
translation-owning F3 completion, a complete
`H[h_00,h_0i,h_ij;ports]`, a Feshbach-dressed momentum density, continuum
stress conservation, a massless mode, gravity, Newton's law, or `G`.

## 1. Frozen question

GD supplied the first exact translation-owning completion of the GA/FU09b
charge transfer.  Each encoded link/reservoir pair owns equal-and-opposite
recoil, and on GD's source-independent full-code hold the inherited FV/FY
nonidentity Hamiltonian and spatial source remain unchanged modulo one common
reference scalar.

GD explicitly did **not** derive `T^{0j}`, the placement of its auxiliary
factor edge in physical diamond space, a native divergence, or a spacetime
source family.  FZ independently proved that FY's complete spatial source
has a nonzero diagonal longitudinal component under FY's supplied continuum-
embedding contraction, while preserving the ceiling that the physical native
divergence `Delta_m` is still undefined.

This lane asks the smallest typed composition question:

> If the only momentum density added to the frozen FY response is the direct
> GD encoded link/reservoir momentum, can it close the source-off stress Ward
> identity under the supplied embedding contraction?

The answer is no.  The result is conditional on the supplied embedding
contraction, or on a later native divergence proved proportional or equivalent
to it on the scored source.  It is not an unconditional no-go for an unknown
physical `Delta_m`.

More exactly, the excluded assignment is the bare, directly projected,
scalar-weighted GD `P_L/P_R` assignment under the frozen FY spatial source.
No stronger class of momentum-density construction is excluded.

## 2. Direct GD momentum density is configuration diagonal

For one encoded link/reservoir pair, GD's two codewords have momenta

\[
 \begin{aligned}
 P_{L_e}^jW_e&=W_e\left(p_e^jI-
               {\hbar\kappa_e^j\over2}Z_e\right),\\
 P_{R_e}^jW_e&=W_e\left(p_e^jI+
               {\hbar\kappa_e^j\over2}Z_e\right).
 \end{aligned}                                      \tag{GJ01}
\]

This is just GD06--GD09 written as an operator on the encoded logical link.
It holds for every kick magnitude and direction.  Give the two factors any
prospectively frozen local or Fourier weights `a_e(m),b_e(m)`.  The direct
encoded momentum density is then

\[
 T^{0j}_{{\rm GD,dir},m}
 =\sum_e\left[a_e(m)P_{L_e}^j+b_e(m)P_{R_e}^j\right]
 =c_m^jI+\sum_ed_{e,m}^jZ_e .                       \tag{GJ02}
\]

It is diagonal in the incidence-configuration basis.  Multiplying its terms
by arbitrary support phases, moving a reservoir between declared support
labels, or relabeling factor edges changes `c,d` but not this operator type.
The ice projector is itself diagonal in that basis, so the direct projected
density remains diagonal.

For every finite operator `H` and every operator `D` diagonal in a basis,

\[
 [H,D]_{nn}=H_{nn}D_{nn}-D_{nn}H_{nn}=0.            \tag{GJ03}
\]

Consequently

\[
 \boxed{
 \operatorname{diag}_{\rm config}
 \left({i\over\hbar}
 [H_6,T^{0j}_{{\rm GD,dir},m}]\right)=0 .}          \tag{GJ04}
\]

No spectral assumption, perturbative denominator, kick calibration, or
continuum approximation enters (GJ04).

## 3. The exact supplied-embedding diagonal does not vanish

FZ's supplied `m=1` embedding has `k` proportional to `(7,15,-17)`.  Its
exact cyclotomic calculation gives a diagonal pair-source longitudinal
witness

\[
\begin{split}
 w(z)={}&-14+28z^2+14z^8+14z^{10}+14z^{24}-14z^{26}
 -28z^{34}\\
 &-14z^{42}-14z^{48}-14z^{56}+28z^{58}
 \not\equiv0\pmod{\Phi_{240}} .                    \tag{GJ05}
\end{split}
\]

Because `deg w=58<deg Phi_240=64`, its displayed nonzero remainder is already
decisive.  FY's complete H2/H4/H6 diagonal lift multiplies the same native
pair ledger by

\[
 f_E(x)=1-x^2-{37\over12}x^4-{16247\over900}x^6 .  \tag{GJ06}
\]

At the first frozen sample,

\[
 x={2\over5},\qquad \rho={15625\over504},\qquad
 \rho f_E={2415673\over113400}\ne0.                \tag{GJ07}
\]

Thus the complete diagonal entry of the supplied-embedding contraction

\[
 L_{{\rm emb},m}^{j}:=(\Delta_{{\rm emb},m})_i
                       T_m^{ij}                    \tag{GJ08}
\]

is nonzero.  The H6 ring source is strictly off diagonal and cannot cancel a
diagonal matrix entry.

### Theorem `FRSEWO-1` -- direct-density conditional obstruction

On GD's source-independent encoded hold, take the temporal momentum density
to consist only of (GJ02), retain FY's complete spatial source, and use FZ's
supplied embedding contraction, or a native divergence already proved
equivalent to it on that source.  Then

\[
 \boxed{
 {i\over\hbar}[H_6,T^{0j}_{{\rm GD,dir},m}]
 +L_{{\rm emb},m}^{j}\ne0 .}                       \tag{GJ09}
\]

#### Proof

The first term has zero configuration diagonal by (GJ04).  The second has
the nonzero exact diagonal entry (GJ05)--(GJ07).  Equality of two finite
operators requires equality of every diagonal entry.  Therefore the sum
cannot vanish.  QED.

This proves that adding GD's direct half-kick momentum labels to FY does not
complete RF3.  It does not prove that the full F3 parent lacks a lawful
momentum density.

## 4. The off-diagonal recoil current does not evade the theorem

The strongest existing apparent escape is GD's exact recoil current.  Put

\[
 Y_e=A_eU_e,qquad
 H_{{\rm flip},e}=-h_e(Y_e+Y_e^\dagger),
\]

so GD derives

\[
 J^P_{L_e\to R_e}
 =ih_e\kappa_e(Y_e-Y_e^\dagger),                  \tag{GJ10}
\]

and

\[
 \dot P_{L_e}+J^P_{L_e\to R_e}=0,qquad
 \dot P_{R_e}-J^P_{L_e\to R_e}=0.                \tag{GJ11}
\]

Equation (GJ10) is off diagonal.  If it were deliberately substituted for a
momentum **density**, its commutator with `H_flip` could indeed have a
diagonal component.  The verifier checks this hostile countercase exactly.
It does not evade FRSEWO-1 for three independent reasons.

1. GD derives (GJ10) as the oriented flux in (GJ11), not as `T^{0j}`.  It is
   already on the divergence/current side of its two-factor balance and
   cannot be reassigned to the density side without changing the physical
   source definition.
2. A single link transfer leaves the local degree-two ice fiber.  Therefore
   its direct projection obeys
   `P_ice J^P_(L_e->R_e) P_ice=0`.  A nonzero effective ice operator would
   have to be derived by carrying an `h_0i` source and this current through
   the same virtual histories; GD did not do that.
3. GD's exact incidence is the auxiliary link-support/reservoir graph.  It
   closes (GJ11), but `kappa_e`, reservoir placement, and the map of that
   factor edge into FY's A/B/link-midpoint diamond supports are unselected.
   Hence its factor-edge divergence is not FZ's supplied embedding
   contraction and is not a derived native diamond `Delta_m`.

### Theorem `FRSEWO-2` -- typed current boundary

The existing off-diagonal recoil current closes the GD auxiliary L/R impulse
ledger exactly, but supplies no existing nonzero projected `T^{0j}` and no
typed identity with FY's spatial source.  It cannot be reassigned by notation.
A Feshbach-dressed momentum density derived from a complete `h_0i` source
remains open and is not excluded by FRSEWO-1.

## 5. Native divergence, local indexing, and boundaries

The diagonal argument survives arbitrary scalar Fourier phases and local
support regroupings in (GJ02).  It therefore is not an artifact of choosing
one numerical reservoir phase.

It is nevertheless load-bearing that FZ has not derived the native physical
`Delta_m`.  A different lawful divergence could annihilate the diagonal FY
source or combine it with support/field stresses not present in `FV-PURE`.
FRSEWO-1 becomes a native no-go only after that divergence is proved
proportional or equivalent to the supplied contraction on the scored source.

Likewise, an active boundary contribution would change the operator balance
to

\[
 {i\over\hbar}[H,T^{0j}_m]+(\Delta_m)_iT^{ij}_m
 +B_m^j=0.                                         \tag{GJ12}
\]

GD's response hold puts the outer port in an invariant off sector.  A new
active boundary with a diagonal `B_m^j` could alter (GJ09), but it changes the
parent and must own its charge, momentum, work, placement, and source
derivatives.  FRSEWO makes no no-go claim for such a completed parent.

## 6. Why contacts cannot repair the source-off identity

Write an allowed spatial contact as `R(j)=O(j^2)`.  Then

\[
 {\partial R\over\partial j_{ij}}\bigg|_{j=0}=0,
 \qquad
 {\partial^2R\over\partial j_{ij}\partial j_{kl}}
 \bigg|_{j=0}\ \text{may be nonzero}.              \tag{GJ13}
\]

The first equation means that `R` changes neither the source-off spatial
operator nor (GJ09).  The second says that it can contribute a seagull to a
differentiated response identity.  Similarly, a mixed term `h_0i j_kl C`
has a nonzero mixed contact but contributes neither direct source operator
when both sources vanish.

Equal-time contacts in FZ12 arise after differentiating a time-ordered or
retarded correlator.  They enforce the response Ward identity once the
operator continuity equation exists; they do not alter the source-off
operator equation from which that response identity follows.

### Theorem `FRSEWO-3` -- contact ceiling

No prospectively frozen `O(j^2)` or mixed `h_0j` contact can cancel the first-
source, source-off diagonal mismatch in (GJ09).  A term with a nonzero first
derivative, a nonzero background source, or an active boundary is a different
first-source parent and requires a new complete rank and Ward calculation.

## 7. Feshbach and positive-construction interface

FY is an effective ice-space calculation.  A complete physical observable
need not equal the direct projection of its microscopic density: source
insertions in the `P/Q` blocks and derivatives of the resolvent can generate
additional effective operators.  FRSEWO-1 therefore rules out only the
literal direct GD-density completion.  It does not rule out the calculation
FZ actually requested.

The recommended next positive construction must freeze before projection one
family

\[
 H[h_{00},h_{0i},h_{ij};\mathrm{ports}],           \tag{GJ14}
\]

with all of the following in the same parent:

1. FU09b/GA charge-conserving transfers;
2. GD recoil plus the momentum owner of the FU pair interaction itself;
3. terminal, reservoir, field/support, controller, and boundary energy and
   momentum densities;
4. the full strain dependence rather than only its first derivative;
5. all second source derivatives and equal-time contacts; and
6. a native factor-to-diamond support map from which `Delta_m` is derived.

Only then should the complete density, current, stress, port, and contact
packet be carried through the fixed ice Feshbach reduction.  That calculation
may create a lawful off-diagonal effective `T^{0j}`, modify FY's spatial
source, derive a different native divergence, or fail.  Every outcome must be
recomputed rather than inferred from this no-go.

This interface is a constructive recommendation, not a logical-necessity
theorem.  FRSEWO-1 alone does not prove that pair-field/support momentum is
the only repair.  Dynamical position-weighted localization, interaction
contributions to `T^{0j}`, a modified spatial source, or another lawful
source-complete parent remain untested.

## 8. Disposition

`EXACT_BARE_DIRECTLY_PROJECTED_SCALAR_WEIGHTED_GD_PL_PR_DENSITY_IS_CONFIGURATION_DIAGONAL__ITS_COMMUTATOR_HAS_ZERO_CONFIGURATION_DIAGONAL__FROZEN_FY_FZ_SPATIAL_SOURCE_HAS_A_NONZERO_EXACT_DIAGONAL_LONGITUDINAL_ENTRY_UNDER_THE_SUPPLIED_EMBEDDING_CONTRACTION__THAT_BARE_DIRECT_DENSITY_THEREFORE_FAILS_THE_CONDITIONAL_SOURCE_OFF_WARD_TEST__OFF_DIAGONAL_GD_RECOIL_CURRENT_CLOSES_ONLY_ITS_AUXILIARY_FACTOR_EDGE_FLUX_LEDGER_AND_CANNOT_BE_REASSIGNED_AS_T0J__O_J2_AND_MIXED_CONTACTS_CANNOT_REPAIR_THE_FIRST_SOURCE_OFF_OPERATOR_IDENTITY__NATIVE_DELTA_POSITION_WEIGHTED_LOCALIZATION_INTERACTION_T0J_FESHBACH_DRESSED_MOMENTUM_MODIFIED_SPATIAL_SOURCE_ACTIVE_BOUNDARY_CONTINUUM_GRAVITY_AND_G_OPEN`

This result eliminates one insufficient completion.  It does not refute F3,
RGRL, Gravity Formation Theory, or gravity emergence.

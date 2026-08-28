# Finite homogeneous response of the FV rank-witness source

**Lane ID:** `GRA-FW-F3-Q4-CHPR-V001`  
**Short name:** `CHPR`  
**Date:** 2026-08-27  
**Claim class:** exact finite-sector source construction; exact operator and
commutator ranks; deterministic complete finite-matrix ground-state Lehmann
response; exact algebraic pole and residue certificate; strict
finite-sector/zero-momentum boundary

**Status:**
`FV_WITNESS_TWO_PIECE_SOURCE_RESTRICTS_TO_RANK5_MOD_IDENTITY_ON_THE_FO_180_STATE_COMPONENT__HOMOGENEOUS_SOURCE_TO_COMMUTATOR_RANK3__GROUND_STATE_RETARDED_AND_FIRST_MOMENT_RANK2__TWO_EXACT_RANK1_POLES__GENERATED_DIAGONAL_Q2_Q4_Q6_DERIVATIVES_NOT_INCLUDED__FINITE_COMPONENT_ZERO_MOMENTUM_SUBOPERATOR_RESULT_NOT_THE_COMPLETE_FIXED_ORDER_CTP_RESPONSE`

**Not claimed:** that the six FV operators lose rank on the full projected
family; that the complete generated diagonal derivatives and folds through
H6 have been included; that the displayed rank two is an upper bound on the
complete fixed-order response; that a spatially resolved source has been
differentiated through the microscopic parent; that the finite component
realizes the thermodynamic state; that an accidental finite-sector
conservation law is a Ward identity; or that a massless tensor pole, RGRL-B,
gravity, or `G` has been derived or excluded.

## 1. Exact question and conditional custody

FV proved that the ideal-Coulomb source-before-Feshbach derivative has six
independent nonidentity projected operators on the covering-matched q4
diamond-ice **family**:

\[
 \operatorname{rank}_{\rm nonid}D_jH_{\rm eff}^{(\le6)}
 =2_E+1_{A_1}+3_{T_2}=6.                         \tag{FW01}
\]

FV correctly withheld a retarded or CTP rank.  Its complete fixed-`P2`
source through leading H6 is

\[
 Q_{\rm eff}^{(\le6)}=Q_{\rm pair}^{(0)}
 +Q_{\rm diag}^{(2)}+Q_{\rm diag}^{(4)}+Q_{\rm diag}^{(6)}
 +Q_{\rm ring}^{(6)}+Q_{\rm id}^{(\le6)}.        \tag{FW01a}
\]

FV computed the exact rank-closing direct and ring entries, but it did not
compute the generated diagonal operators in (FW01a), because their values
were unnecessary for its matrix-element rank witnesses.  Those omitted
diagonal operators can matter in a retarded calculation.

This lane therefore asks a deliberately smaller dynamical question without
changing the Hamiltonian:

> What response do FV's two explicit rank-closing pieces have in the exact FO
> 180-state translation-closed ring component and its unique sector ground
> state?

Define the load-bearing scope marker

\[
 \boxed{Q_{\rm FV-WITNESS}:=Q_{\rm pair}^{(0)}
       +Q_{{\rm ring,irr}}^{(6)}.}                \tag{FW01b}
\]

`FV-WITNESS` is a source-suboperator diagnostic.  It omits
`Q_diag^(2,4,6)`, every differentiated fold that contributes there, and
identity/reference terms.  It is not called the complete FV source, the
complete through-H6 source, or the FQ CTP packet.  Its response rank is not an
upper bound on any of those complete objects.

The construction is explicitly conditional on all FU/FV physical-solder
premises `S1`--`S10`, including the `FV-PURE` premise that excludes an
unfrozen residual nonidentity affine kernel at the scored order.  This lane
does not independently derive that physical completion.  It pins the final
FV core bytes and the independently audited FO finite sector in
`DEPENDENCIES.sha256`.

This is a lawful first finite response screen because it composes two already
earned results: FV supplies two exact source suboperators, while FO supplies a
fully enumerated source-off Hamiltonian, translation action, unique sector
ground state, and complete eigensystem.  No interaction, projector, kinetic
term, or fitted response weight is added.  The complete generated-diagonal
calculation is the immediate successor, not an optional refinement.

## 2. Frozen source and scale

Use the FO sector Hamiltonian in units `J6=1`,

\[
 H=-\sum_C B_C,
 \qquad E_0=-2(1+\sqrt2),                         \tag{FW02}
\]

on the 180 states generated from FO's fixed winding seed.  It contains 420
undirected ring transitions and is closed under the cyclic translation
`x -> x+1`.

For the direct pair source retain the FV Coulomb slope `lambda=-1/2`.  The
source tensor in `U_d` units is

\[
 Q_{\rm pair}^{ij}/U_d=-{1\over2}
 \sum_v\sum_{a<b}\widehat R_{ab}^{ij}P_{v,ab}.    \tag{FW03}
\]

For a flippable hexagon with missing q4 label `d`, FV's complete
forward/reverse differentiation gives, after converting its FV11 coefficient
to `J6` units,

\[
 Q_{{\rm ring},d}^{ij}/J_6
 =-{31\over6}\delta^{ij}+{9\over2}D_d^{ij}.       \tag{FW04}
\]

Equation (FW04) includes the differentiated six numerators and all
endpoint-referenced resolvents already certified by FV.  It is not a
post-Feshbach hand weight.

Define the physical scale ratio

\[
 \rho={U_d\over J_6}>0.                           \tag{FW05}
\]

All formulas below use `J6=1`; therefore the witness operator is
`Q=rho Q_pair+Q_ring`.  No conclusion depends on the artificial numerical
choice `rho=1`.  The verifier separately replays `rho=1` and `rho=2` and
checks the analytic direct/ring factorization.  This scale convention does
not restore the omitted generated diagonal terms.

Use the Frobenius-orthonormal tensor basis

\[
\begin{aligned}
 A&={xx+yy+zz\over\sqrt3},\\
 E_1&={xx-yy\over\sqrt2},\qquad
 E_2={xx+yy-2zz\over\sqrt6},\\
 T_{xy}&={2xy\over\sqrt2},\qquad
 T_{xz}={2xz\over\sqrt2},\qquad
 T_{yz}={2yz\over\sqrt2}.
\end{aligned}                                    \tag{FW06}
\]

The direct term occupies only `A1+E`; the ring term occupies only `A1+T2`.
Every homogeneous source commutes with cyclic translation, so this packet is
strictly a momentum-zero calculation.

## 3. Finite-component operator identities

The six family operators do not remain six modulo identity after restriction
to this one small winding component.  Exact enumeration gives

\[
 \boxed{Q_A={60\rho\,\mathbf1+11H\over\sqrt3}.}  \tag{FW07}
\]

Thus the uniform `A1` source is nonidentity but conserved: it is an energy
rescaling plus an identity reference.  In the `E` block define

\[
 Q_{E_c}={\sqrt3\over2}Q_{E_1}+{1\over2}Q_{E_2},
 \qquad
 Q_{E_a}=-{1\over2}Q_{E_1}+{\sqrt3\over2}Q_{E_2}. \tag{FW08}
\]

They correspond respectively to the normalized diagonal tensors
`diag(2,-1,-1)/sqrt(6)` and `diag(0,1,-1)/sqrt(2)`.  On the entire selected
component,

\[
 \boxed{Q_{E_c}=16\sqrt6\,\rho\,\mathbf1,}        \tag{FW09}
\]

while `Q_Ea` is nonconserved.  Finally,

\[
 [H,Q_{T_{yz}}]=0,                               \tag{FW10}
\]

and its unique-ground-state eigenvalue is `3 sqrt(2)-6`.

Consequently the six homogeneous sources have exact operator rank five after
the identity is quotiented, not six.  This does **not contradict FV**.  FV's
rank-six theorem uses local diagonal differences and four distinct global
ring matrix-element witnesses on the covering-matched family.  FW restricts
all sources to one uniform query on one 180-state connected component; one
of FV's `E` directions is constant there.

The axis singled out in (FW09)--(FW10) belongs to the chosen quotient and
winding component.  It is not an `S4`-invariant thermodynamic statement.

## 4. State-independent commutator rank

Before choosing a density matrix, form the Hilbert-Schmidt Gram matrix

\[
 C_{AB}=\operatorname{Tr}\!\left([H,Q_A]^\dagger[H,Q_B]\right). \tag{FW11}
\]

In the order `(A,E1,E2,Txy,Txz,Tyz)`, exact finite enumeration gives

\[
 C=
 \begin{pmatrix}
 0&0&0&0&0&0\\
 0&960\rho^2&-960\sqrt3\rho^2&0&0&0\\
 0&-960\sqrt3\rho^2&2880\rho^2&0&0&0\\
 0&0&0&25920&0&0\\
 0&0&0&0&25920&0\\
 0&0&0&0&0&0
 \end{pmatrix}.                                  \tag{FW12}
\]

Therefore

\[
 \boxed{\operatorname{rank}(Q\mapsto[H,Q])=3}    \tag{FW13}
\]

for every `rho != 0`.  The nonconserved subspace is spanned by
`(E_a,Txy,Txz)`; `(A,E_c,Tyz)` is conserved on this component.  This
source-to-commutator rank is an operator statement.  It still does not equal
a state-specific retarded rank.

## 5. Exact ground-state spectral response

Let `|0>` be the unique FO sector ground state.  For positive excitation
energy `Delta_n`, define the residue vectors

\[
 v_{n,A}=\langle n|Q_A|0\rangle,
 \qquad R(\Delta)=\sum_{n:\Delta_n=\Delta}v_nv_n^{\mathsf T}. \tag{FW14}
\]

The complete 180-state Lehmann calculation has support at only two gaps:

\[
 \Delta_1=2+2\sqrt2,\qquad
 \Delta_2=4+2\sqrt2.                             \tag{FW15}
\]

Both excited energies lie in the zero cyclic-momentum block.  In the basis
(FW06), their exact rank-one residue vectors can be chosen as

\[
\begin{aligned}
 r_1&=\left(0,{\rho\over\sqrt2},
 -\rho\sqrt{3\over2},-{3\over\sqrt2},
 -{3\over\sqrt2},0\right),\\
 r_2&=\left(0,0,0,{3\over\sqrt2},
 -{3\over\sqrt2},0\right),
\end{aligned}                                    \tag{FW16}
\]

so that

\[
 R_1=r_1r_1^{\mathsf T},\quad
 R_2=r_2r_2^{\mathsf T},\quad
 \|r_1\|^2=2\rho^2+9,\quad \|r_2\|^2=9,
 \quad r_1\cdot r_2=0.                           \tag{FW17}
\]

The algebraic statement is not inferred from a numerical pattern.  Uniform
translation reduces the calculation exactly to the six-orbit block

\[
 H_0=\begin{pmatrix}
0&-1&-1&-1&-1&-2\\
-1&0&-1&-1&0&-1\\
-1&-1&0&0&-1&-1\\
-1&-1&0&0&-1&-1\\
-1&0&-1&-1&0&-1\\
-2&-1&-1&-1&-1&0
\end{pmatrix},                                   \tag{FW17a}
\]

with normalized ground vector

\[
 g_0=\left({1\over2},{1\over2\sqrt2},{1\over2\sqrt2},
 {1\over2\sqrt2},{1\over2\sqrt2},{1\over2}\right).          \tag{FW17b}
\]

The exact projectors onto excited energies zero and two are the rational
polynomials

\[
 P_0={(H_0^2+4H_0-4I)(H_0-2I)\over8},\qquad
 P_2={(H_0^2+4H_0-4I)H_0\over16}.               \tag{FW17c}
\]

The uniform source blocks lie on the exact denominator-twelve lattice.
Direct multiplication of those blocks, (FW17b), and (FW17c) gives (FW16),
and the centered source vector has zero remainder outside `P0+P2`.  This is
the load-bearing algebraic proof that no third pole is present.  Complete
180-state diagonalization and degenerate-eigenspace residue sums are an
independent numerical replay.

Hence the exact algebraic two-pole retarded form of `FV-WITNESS` is

\[
 \boxed{
 \chi^R(\omega)=
 R_1\!\left({1\over\omega-\Delta_1+i0}
            -{1\over\omega+\Delta_1+i0}\right)
 +R_2\!\left({1\over\omega-\Delta_2+i0}
            -{1\over\omega+\Delta_2+i0}\right).} \tag{FW18}
\]

It follows that

\[
 \boxed{\operatorname{rank}_{\rm ground}\chi^R=2.} \tag{FW19}
\]

The positive static Kubo matrix

\[
 K={2R_1\over\Delta_1}+{2R_2\over\Delta_2}       \tag{FW20}
\]

has the two nonzero eigenvalues

\[
 (2\rho^2+9)(\sqrt2-1),
 \qquad {9\over2}(2-\sqrt2).                    \tag{FW21}
\]

Thus neither static coefficient rank nor operator commutator rank may be
reported as the ground-state retarded rank.

## 6. First moment and the nonconserved dark direction

With the FQ convention

\[
 M_n^{AB}=\langle0|[(\operatorname{ad}_H)^nQ_A,Q_B]|0\rangle, \tag{FW22}
\]

the equal-time moment is `M0=0` and

\[
 \boxed{M_1=-2(\Delta_1R_1+\Delta_2R_2),
 \qquad \operatorname{rank}M_1=2.}               \tag{FW23}
\]

Therefore `n_star=1`, but its common-normalization rank is two rather than
six.  This `FV-WITNESS` packet has a degenerate canonical moment at
homogeneous momentum.  It does **not** establish that the complete source in
(FW01a) fails the `Q4-BLOCK-STRAIN-CTP` condition, because its generated
diagonal terms have not yet been included.

The four-dimensional ground-response nullspace is not four conservation
laws.  Three nulls are the component-conserved directions `(A,E_c,Tyz)`.
The fourth is

\[
 Q_{\rm dark}=3Q_{E_a}-\sqrt2\rho\,
 {Q_{T_{xy}}+Q_{T_{xz}}\over\sqrt2}.             \tag{FW24}
\]

After subtracting its ground expectation,
`Q_dark|0>=0`, while `[H,Q_dark] != 0`.  It is a state-specific dark
direction, **not a Ward identity**, constraint, or gauge generator.

## 7. Disposition and exact ceiling

This result supplies the first exact composition

\[
 \boxed{
 \text{FV rank-six family source}
 \supset Q_{\rm FV-WITNESS}
 \longrightarrow
 \text{FO finite homogeneous response}
 =\text{operator rank }5
 \to\operatorname{rank}\operatorname{ad}_H 3
 \to\operatorname{rank}\chi^R_{|0\rangle}2.}     \tag{FW25}
\]

That hierarchy is the physics result for the two rank-closing pieces.  It
proves by construction why an off-shell projected operator rank cannot be
promoted directly to a dynamical metric field.  It also identifies the
homogeneous `A1` witness direction as an energy rescaling.  The latter is a
conservation channel in this finite calculation, not yet an emergent
gravitational constraint.

This is a **finite-sector, homogeneous-source result**.  It does not close the
full F3 collective route because:

1. the FO component is one small winding sector, not the complete ice Hilbert
   space or a nested thermodynamic family;
2. FV derived a uniform affine derivative, not a spatially resolved
   source-before-Feshbach density with a frozen term-ownership convention;
3. the source is at `k=0`, where energy conservation necessarily suppresses
   the `A1` response;
4. `Q_diag^(2,4,6)` and its differentiated folds are omitted, so even the
   complete fixed-order H6 source response remains to be computed;
5. only leading `H6` source-off dynamics is included, while a complete FQ
   test requires the fixed through-`H8` parent, contacts, ports, density
   family, and blocked CTP action; and
6. no ungauge-fixed Ward/constraint calculation is performed.

The immediate lawful no-laboratory calculation is to compute every generated
diagonal derivative and fold in `Q_diag^(2,4,6)` on the same FO component and
rerun the operator/commutator/Lehmann hierarchy for the complete fixed-order
source.  Only after that should the program derive and freeze a local block
source, carry it through the complete finite parent, and compute nonzero-
momentum response on a prospectively matched family.

Nothing here establishes or excludes a thermodynamic massless tensor,
RGRL-B, gravity, or a numerical value of `G`.

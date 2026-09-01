# Authenticated locked-E source overlap and retained-K loop-selection obstruction

**Lane:** `LANE_CROSS_RFT_GRA_GL6AQ_AUTHENTICATED_E_LOOP_SELECTION_OBSTRUCTION_V001`  
**Short name:** `GL6AQ V001`  
**Date:** 2026-08-31  
**Status:** isolated author theorem with exact replay; independent hostile review
required before freeze  
**Inputs:** the exact frozen `GL6AM` snapshot sealed through its distinct hostile
audit packet, and the exact author/audit-sealed `GL6AN` snapshot

**Claim class:** exact authenticated pair-source overlap with the local locked
`E` coordinate; exact zero direct/linear selection rule for a one-cell
retained-`K` transverse source; exact six-support gating of the sealed
alternating-hexagon `E` displacement; sharp obstruction to promoting that
finite linked coefficient to a nonzero stationary bulk defect contrast.

**Not claimed:** preparation or selection of a locked bulk state; an infinite
locked-sector projector; completeness of the order-six effective Hamiltonian;
strictly positive retarded or dissipative weight; a gap or pole; physical
momentum, length, speed, or cone; stress, Ricci/Einstein response, gravity, or
`G`.

## 1. Exact authenticated source overlap with locked E

At one original parent use the six pair order

\[
 {cal P}=(12,13,14,23,24,34)
\]

and the unsigned vertex/pair incidence matrix

\[
 R_{a,A}={\bf1}_{a\in A}.
\]

The sealed `GL6AN` theorem proves that the two-dimensional space

\[
 {cal E}:=\ker R                                                   \tag{AQ01}
\]

is exactly the local `E` irrep.  For `c in E`, define the authenticated
`GL6AM` pair source/read

\[
 O_x(c):=\sum_{A\in{cal P}}c_A M_{x,A},\qquad
 M_{x,ab}=Z_{(x,a)}Z_{(x,b)}.                                    \tag{AQ02}
\]

It is inserted by the exact finite pulse

\[
 \exp\!\left({i\eta\over2}O_x(c)\right)                          \tag{AQ03}
\]

and is a finite real linear combination of complete pair reads.  Thus no new
source or read operator is postulated.

Let `P_x^(2)` be the ordinary local projector onto `k_x=2`.  For each of the
six locked basis configurations, let `M(z)` be its six pair values.  The
affine locked identity is

\[
 R M(z)=(-1,-1,-1,-1)^T.                                        \tag{AQ04}
\]

Consequently every locked displacement lies in `E`.  More strongly, exact
enumeration of the six local locked configurations gives, for every real
`c in E`,

\[
 {1\over6}\sum_{k(z)=2}\big(c^TM(z)\big)^2
 ={8\over3}\,c^Tc.                                               \tag{AQ05}
\]

The mean vanishes because `c` is orthogonal to the uniform pair vector.
Therefore

\[
 c\ne0\quad\Longrightarrow\quad
 P_x^{(2)}O_x(c)P_x^{(2)}\ne0.                                  \tag{AQ06}
\]

Equations (AQ03)--(AQ06) are an exact nonzero authenticated operator overlap
with the locked `E` coordinate.  They do not select a state.

## 2. What stationary pulse/read positivity does and does not add

For any chosen joint stationary/translation/`S4`-invariant `GL6AM` state
`omega`, contract its positive pair correlation measure with `c`:

\[
 \mu_c({\cal B})=\sum_{A,B}c_A\mu_{AB}({\cal B})c_B\ge0,
 \qquad
 \mu_c(\mathbb R)=
 \omega\!\left((O_x(c)-\omega(O_x(c)))^2\right).                \tag{AQ07}
\]

This is positivity, not strict positivity.  It supplies no nonzero retarded
commutator and no positive dissipative measure without additional state
input.

The normalized product trace `tr_0` makes the distinction exact.  It is a
lawful joint stationary/translation/`S4`-invariant state of the quasi-local
spin algebra.  Pauli orthogonality gives

\[
 \mu_c^{\rm tr}(\mathbb R)=\operatorname{tr}_0(O_x(c)^2)=c^Tc>0
 \quad(c\ne0),                                                   \tag{AQ08}
\]

but tracial cyclicity gives

\[
 \operatorname{tr}_0([\tau_t(A),B])=0                           \tag{AQ09}
\]

for all local `A,B`.  Hence its `GL6AM` retarded kernel is zero even though
its `E` correlation measure has strictly positive total mass.  Positive
stationary homogeneous pulse/read correlation mass is not a proof of
response, a pole, or the order-six locked mechanism; no KMS or passivity
premise is being inserted.

## 3. Exact zero direct and one-cell linear K-to-E projection

On the all-formed branch a finite retained word changes the transverse
coefficient through the sealed `GL6AM` defect

\[
 V_\kappa=h\sum_{e\in D}(1-\kappa_e)X_e,
 \qquad \kappa_e\in\{0,1\}.                                    \tag{AQ10}
\]

In the positive strong-lock regulator of `GL6AN`, one flip creates two unit
degree defects.  Therefore

\[
 P_{\cal Q}X_eP_{\cal Q}=0,
 \qquad P_{\cal Q}V_\kappa P_{\cal Q}=0.                        \tag{AQ11}
\]

The direct Pauli Hilbert--Schmidt overlap also vanishes:

\[
 \operatorname{tr}_0(O_x(c)X_e)=0.                              \tag{AQ12}
\]

There is a complementary one-cell symmetry selection rule.  At the
`S4`-fixed origin, the four port sources `(X_1,...,X_4)` carry the permutation
representation `A1+T2`; it contains no `E`.  The six pair reads carry
`A1+E+T2`.  Thus any `S4`-equivariant one-cell cross kernel `K_{A a}(t)` obeys

\[
 P_E K(t)=0.                                                     \tag{AQ13}
\]

In particular, a stationary homogeneous linear-response reinterpretation of
a one-cell transverse `K` source cannot drive the pair `E` channel.  A
generic multi-cell defect word is not `S4`-closed, so (AQ13) is not extended
to it; `GL6AM` assigns such a word no homogeneous sector split in the first
place.

## 4. The first sealed configuration-changing K coupling is order six

Take the exact alternating hexagon `C` and locked configurations `|i>` and
`|f>` constructed in sealed `GL6AN`.  At each of its three parent nodes, the
loop exchanges one occupied and one unoccupied incident loop link while the
other two incident links remain one occupied and one unoccupied.  If

\[
 \delta M_x:=M_x(f)-M_x(i),
\]

then direct local enumeration gives

\[
 R\delta M_x=0,
 \qquad \delta M_x\ne0,
 \qquad \|\delta M_x\|^2=16.                                   \tag{AQ14}
\]

Thus `delta M_x` is exactly `E`, and the authenticated read with
`c=delta M_x` distinguishes the two locked configurations by

\[
 \langle f|O_x(\delta M_x)|f\rangle
 -\langle i|O_x(\delta M_x)|i\rangle=16.                        \tag{AQ15}
\]

The sealed `GL6AN` path sum enumerates all `6!` orderings and gives
`-63/8`.  If the six transverse terms on `C` carry their physical retained
coefficients `kappa_e`, every leading path flips every edge exactly once.
The energy denominators are unchanged, so the exact linked matrix element is

\[
 \boxed{
 \langle f|H_{\mathrm{eff},{\cal Q}}^{(6)}(\kappa)|i\rangle_C
 =-{63\over8}{h^6\over U_d^5}\prod_{e\in C}\kappa_e.}          \tag{AQ16}
\]

For the binary retained word, removing any one of the six supports kills
this particular leading loop entry.  Equations (AQ14)--(AQ16) prove a
nonzero, genuinely nonlinear retained-`K` gate on an `E`-changing collective
operation.  This is not a nonzero linear projection: (AQ11)--(AQ13) remain
exact.  No claim is made that all lower-order diagonal shifts for a
nonuniform defect word are absent or configuration independent; those terms
are outside the sealed homogeneous scalar census used here.

## 5. Sharp stationary-bulk obstruction

The linked coefficient (AQ16) cannot be inserted into the stationary
`GL6AM` response functional from the sealed premises.

1. `GL6AN` proves (AQ16) on one finite period-four locked projector and only
   as a linked strong-lock coefficient.  It does not construct an infinite
   locked projector, the complete order-six effective generator, or a
   stationary locked bulk state.
2. `GL6AM` makes `gamma^kappa` a lawful bulk defect dynamics, but the original
   homogeneous state is generally not stationary for it.  Its matched
   contrast is nonequilibrium and receives neither the positive stationary
   measure (AQ07) nor an `A1/E/T2` scalar decomposition.
3. The product trace is invariant under every automorphism of this UHF spin
   algebra.  Hence, for every local read and every pair of finite words,

   \[
   \Delta_{\kappa|\kappa'}^{\operatorname{tr}_0}(A,t)
   =\operatorname{tr}_0(\gamma_t^\kappa(A)
    -\gamma_t^{\kappa'}(A))=0.                                 \tag{AQ17}
   \]

Equation (AQ17) refutes a **universal** nonzero stationary defect-contrast
claim, including one inferred merely from positive `E` correlation mass.
The sealed inputs neither prove nor refute the **existence** of some other
selected stationary locked bulk state with a nonzero `K`-conditioned `E`
contrast.  That existential statement remains open because the required
state and thermodynamic effective dynamics have not been supplied.

## 6. Exact conclusion and ceiling

The strongest lawful composition is therefore

\[
 \boxed{\begin{gathered}
 \text{authenticated pair pulse/read}\ \longrightarrow\
 \text{nonzero local locked-}E\text{ operator overlap};\\
 \text{retained-}K\text{ source}\ \longrightarrow\
 \text{zero direct and one-cell linear }E\text{ projection};\\
 \text{six retained loop supports}\ \longrightarrow\
 -{63\over8}{h^6\over U_d^5}\prod_{e\in C}\kappa_e
 \text{ on an }E\text{-changing finite linked move};\\
 \text{sealed stationary bulk composition}\ \not\Longrightarrow\
 \text{nonzero defect }E\text{ contrast}.
 \end{gathered}}                                                \tag{AQ18}
\]

Closing the last implication requires an independently justified bulk-state
preparation/selection, nonzero-variance and nonzero-response hypotheses, and
a controlled thermodynamic order-six effective-dynamics construction.  No
state, pole, physical momentum, cone, stress/Ricci/Einstein law, gravity, or
`G` is inferred here.

`PASS__AUTHENTICATED_PAIR_SOURCE_HAS_EXACT_NONZERO_LOCAL_LOCKED_E_OVERLAP__STATIONARY_CORRELATION_MEASURE_POSITIVE_BUT_NOT_STRICT_RESPONSE__TRANSVERSE_RETAINED_K_SOURCE_HAS_ZERO_DIRECT_AND_ONE_CELL_LINEAR_E_PROJECTION__SIX_RETAINED_SUPPORTS_GATE_THE_E_CHANGING_HEXAGON_BY_MINUS63_OVER8_TIMES_THEIR_PRODUCT__UNIVERSAL_NONZERO_STATIONARY_DEFECT_CONTRAST_REFUTED_BY_TRACE__EXISTENTIAL_SELECTED_LOCKED_BULK_RESPONSE_OPEN__NO_STATE_POLE_MOMENTUM_CONE_RICCI_GRAVITY_OR_G`

# q4 normalized operator-neutral CTP kernel theorem

**Lane:** `LANE_CROSS_RFT_GRA_GL6W_Q4_NORMALIZED_OPERATOR_NEUTRAL_CTP_KERNEL_V001`  
**Short name:** `GL6W V001`  
**Date:** 2026-08-31  
**Status:** author packet frozen after independent hostile pre-freeze review;
post-freeze custody audit pending  
**Claim class:** exact formed-KEEP six-source active-system CTP response on the literal
GL6T/GL6U star; closed all-time factorized kernel; exact finite interacting
spectral compiler; exact direct affine-system seagull census and full-apparatus
Schur boundary

**Not claimed:** an unconditional average over every formation branch; that
pair observables are themselves records; a continuum metric source; a common
clock/matter/EM localization map; intercell gluing; a causal cone;
Ward/Bianchi closure; an infrared pole or operator; Ricci response, gravity,
or `G`; a completed controller/work/boundary/reference CTP or causal Schur
kernel.

## 1. Same-parent source and CTP convention

On the formed-KEEP branch of audited GL6T/GL6U, the active star Hamiltonian is

\[
 H_K=-h\sum_{a=1}^4X_a+\Delta d
 +U_d\left[(d-d_\star)^2+
       \sum_{a=1}^4(n_a-d_\star)^2\right],
 \quad d=\sum_an_a,\quad n_a={1-Z_a\over2}.          \tag{W01}
\]

The state at the response entrance is

\[
 \rho_\tau=e^{-iH_K\tau/\hbar}|0^4\rangle\langle0^4|
             e^{iH_K\tau/\hbar}.                    \tag{W02}
\]

GL6V supplies six exact pulse sources for

\[
 M_A=Z_aZ_b,\qquad A=(ab)\in{\cal E}_4.             \tag{W03}
\]

Freeze one common registered source clock `T_*` and define

\[
 E_\star={\hbar\over T_\star},\qquad Q_A=E_\star M_A,
 \qquad J_A(t)=j_AT_\star f_A(t).                    \tag{W04}
\]

Here `f_A` is the GL6V inverse-time profile with integral one.  Then its exact
source Hamiltonian is in the BS21 half-source convention,

\[
 H_{\rm src}(t)=-{1\over2}\sum_AJ_A(t)Q_A
 =-{\hbar\over2}\sum_Aj_Af_A(t)M_A.                 \tag{W05}
\]

`E_*` is measured probe-clock normalization, not a fitted bulk coupling or a
candidate gravitational constant.  A different common `T_*` is a declared
source-coordinate change and must be carried through any later same-query
metric calibration.

On the two CTP branches use

\[
 H_\pm=H_K-{1\over2}\sum_AJ_{\pm,A}Q_A,\qquad
 J_c={J_++J_-\over2},\quad J_\Delta=J_+-J_-.         \tag{W06}
\]

Equation (W06) is insertion-functional notation for GL6V's exact finite
pulse list, with parent evolution off or refocused during each source
sandwich and source-off evolution between insertions.  It does not promote
that scheduled dilation to arbitrary simultaneous smooth evolution under
`H_K+H_src(t)`.

With `W=-i hbar log Z`, define

\[
 q_B(t)=2{\delta W\over\delta J_{\Delta,B}(t)}
       =\langle Q_B(t)\rangle,\qquad
 {\cal G}^R_{BA}(t,s)=
 2{\delta^2W\over\delta J_{\Delta,B}(t)\delta J_{c,A}(s)}.
                                                               \tag{W07}
\]

Direct variation of the unitary gives the exact physical convention

\[
 \boxed{{\cal G}^R_{BA}(t,s)=
 {i\over2\hbar}\Theta(t-s)
 \langle[Q_B(t),Q_A(s)]\rangle_{\rho_\tau}.}         \tag{W08}
\]

This sign is fixed by (W05).  Relative to the GL6T diagnostic
`chi^R_MM=-(i/hbar)Theta<[M_B,M_A]>`,

\[
 \boxed{{\cal G}^R=-{E_\star^2\over2}\chi^R_{MM}.} \tag{W09}
\]

Thus a source-sign convention cannot be changed downstream without changing
the reported response.

Define the symmetric connected physical noise

\[
 {\cal N}_{BA}(t,s)={1\over2}
 \langle\{\delta Q_B(t),\delta Q_A(s)\}\rangle.     \tag{W10}
\]

Then

\[
 2W_{J_cJ_\Delta}={\cal G}^A=({\cal G}^R)^T,\qquad
 W_{J_\Delta J_\Delta}={i\over4\hbar}{\cal N},\qquad
 W_{J_cJ_c}=0.                                      \tag{W11}
\]

Equations (W08)--(W11), together with the mean, reconstruct the full
quadratic **active-system pulse-insertion** functional in this declared
six-source sector.

## 2. Closed kernel for the factorized comparator

Set `U_d=0` and put

\[
 \epsilon=\sqrt{\Delta^2+4h^2},\qquad
 \sigma={2h\over\epsilon},\qquad c={\Delta\over\epsilon},
 \qquad\omega={\epsilon\over\hbar}.                 \tag{W12}
\]

For `0<=s<t` after the response entrance, define

\[
\begin{aligned}
 z_t&=c^2+\sigma^2\cos[\omega(\tau+t)],
 &z_s&=c^2+\sigma^2\cos[\omega(\tau+s)],\\
 u&=c^2+\sigma^2\cos[\omega(t-s)],\\
 v&=\sigma^2c\{\sin[\omega(\tau+t)]
       -\sin[\omega(\tau+s)]-\sin[\omega(t-s)]\}.
\end{aligned}                                        \tag{W13}
\]

The exact one-link two-time correlator is

\[
 \boxed{\langle Z(t)Z(s)\rangle=u+iv.}              \tag{W14}
\]

Let `A_L` be the adjacency matrix of the line graph `L(K4)`.  Pair
factorization now gives the full active-system physical retarded kernel

\[
 \boxed{{\cal G}^R(t,s)=
 -{E_\star^2\over\hbar}\Theta(t-s)
 \left[2uv\,I_6+vz_tz_s\,A_L\right].}              \tag{W15}
\]

Its `S4` eigenchannels are

\[
\begin{aligned}
 {\cal G}^R_{A_1}&=-{2E_\star^2v\over\hbar}
                  (u+2z_tz_s),\\
 {\cal G}^R_{E_2}&=-{2E_\star^2v\over\hbar}
                  (u-z_tz_s),\\
 {\cal G}^R_{T_2}&=-{2E_\star^2uv\over\hbar}.
\end{aligned}                                        \tag{W16}
\]

The exact symmetric connected kernel is

\[
 \boxed{{\cal N}(t,s)=E_\star^2
 [d_N(t,s)I_6+a_N(t,s)A_L],}                        \tag{W17}
\]

where

\[
 d_N=u^2-v^2-z_t^2z_s^2,\qquad
 a_N=z_tz_s(u-z_tz_s).                              \tag{W18}
\]

At equal times its sector eigenvalues are

\[
 E_\star^2(1-z^2)(1+5z^2),\quad
 E_\star^2(1-z^2)^2,\quad
 E_\star^2(1-z^4)                                  \tag{W19}
\]

on `(A1,E2,T2)`, respectively, and are nonnegative.  The advanced kernel is
the time-reversed transpose of (W15).

At the response entrance,

\[
 \boxed{
 \partial_t{\cal G}^R(t,0)|_{0^+}
 =-{E_\star^2\over2\hbar^2}D(\tau),}              \tag{W20}
\]

where `D` is exactly the audited GL6T double-commutator matrix.  Thus the
closed two-time calculation reproduces, rather than replaces, the entrance
result with the source sign and units fixed.

## 3. Exact interacting finite spectral compiler

For every admitted `U_d>=0`, equations (W01)--(W03) define a finite exact
compiler.  Write

\[
 U(t)=e^{-iH_Kt/\hbar},\quad
 M_A(t)=U(t)^\dagger M_AU(t),\quad
 \mu_A(t)=\operatorname{Tr}[\rho_\tau M_A(t)],      \tag{W21}
\]

\[
 C_{BA}(t,s)=\operatorname{Tr}[\rho_\tau M_B(t)M_A(s)].
                                                               \tag{W22}
\]

The complete interacting active-system kernels are exactly

\[
 \boxed{
 {\cal G}^R_{BA}={iE_\star^2\over2\hbar}\Theta(t-s)
 [C_{BA}(t,s)-C_{AB}(s,t)],}                       \tag{W23}
\]

\[
 \boxed{
 {\cal N}_{BA}=E_\star^2\left[
 {C_{BA}(t,s)+C_{AB}(s,t)\over2}-\mu_B(t)\mu_A(s)
 \right].}                                         \tag{W24}
\]

No continuum or Gaussian approximation occurs in (W21)--(W24).  The
permutation-invariant Hamiltonian admits the exact Schur--Weyl reduction

\[
 J=2\ (5\times5,A_1),\qquad
 J=1\ (3\times3,T_2\text{ multiplicity }3),\qquad
 J=0\ (1\times1,E_2\text{ multiplicity }2),         \tag{W25}
\]

After subtracting the common blank scalar, with `m=-J,...,J`, `r=2-m`,

\[
 (H_J)_{mm}=r\Delta+U_dr(r+1-4d_\star),\qquad
 (H_J)_{m,m+1}=-h\sqrt{J(J+1)-m(m+1)}.             \tag{W26}
\]

Equivalently, in increasing occupation number,

\[
 \mathsf H_{J=2}=\begin{pmatrix}
 \delta_0&-2h&0&0&0\\
 -2h&\delta_1&-\sqrt6h&0&0\\
 0&-\sqrt6h&\delta_2&-\sqrt6h&0\\
 0&0&-\sqrt6h&\delta_3&-2h\\
 0&0&0&-2h&\delta_4
 \end{pmatrix},\quad
 \mathsf H_{J=1}=\begin{pmatrix}
 \delta_1&-\sqrt2h&0\\
 -\sqrt2h&\delta_2&-\sqrt2h\\
 0&-\sqrt2h&\delta_3
 \end{pmatrix},\quad \mathsf H_{J=0}=(\delta_2),   \tag{W26a}
\]

where `delta_r` is the diagonal entry in (W26).  The pair-module projectors
are fixed with no fitted weights by

\[
 P_{A_1}={A_L(A_L+2I)\over24},\qquad
 P_{E_2}={A_L(A_L-4I)\over12},\qquad
 P_{T_2}=-{(A_L-4I)(A_L+2I)\over8}.                \tag{W26b}
\]

The blank enters only the `J=2,m=2` state, while one pair insertion can visit
all three blocks.  Finite diagonalization of (W26) plus the literal Walsh
matrices (W03) therefore computes every entry of (W23)--(W24) at arbitrary
times with no new parameter.

Let `A_O` connect opposite edges of `K4`.  Every resulting `S4`-covariant
six-pair kernel has the unique form

\[
 K=k_dI_6+k_aA_L+k_bA_O,                            \tag{W27}
\]

and hence

\[
 \boxed{
 K_{A_1}=k_d+4k_a+k_b,\quad
 K_{E_2}=k_d-2k_a+k_b,\quad
 K_{T_2}=k_d-k_b.}                                 \tag{W28}
\]

The interacting entrance slope is again (W20), now with audited GL6U
`D=-8hxI_6-4hyA_L`.  In particular, its inherited interaction-owned piece

\[
 y-xz^2=-{16\over3}h^3U_d(\tau/\hbar)^4
       +O((\tau/\hbar)^6)                          \tag{W29}
\]

enters the normalized retarded kernel.  This proves an actual interlink
response contribution; it does not prove a collective pole or continuum
stiffness.

## 4. Reference, direct contact, and full-apparatus Schur boundary

On matched formed BREAK, the transverse term is absent, the blank is an exact
eigenstate, and every `M_A` commutes with the diagonal Hamiltonian.  Therefore

\[
 {\cal G}^{R,{\rm BREAK}}={\cal N}^{\rm BREAK}=0.   \tag{W30}
\]

The KEEP-minus-BREAK contrast is thus exactly the KEEP kernel above.

The GL6V direct source is affine, so its direct **system** seagull is zero.
Its garbage-free query dilation is exactly uncomputed before source-off
response, and the terminal read occurs only after the last insertion.  This
proves that no query-bit parity copy contaminates (W23)--(W24).

It does **not** by itself prove that controller, clock, work,
switching-boundary, reference, source-setting, or full formation-instrument
response blocks vanish.  GL6V treats the schedule as a retained external
control boundary.  Those owners require their own physical CTP ledger before
a complete causal Schur quotient can be taken.  Therefore (W23)--(W24) are
the normalized active-system block that a future complete kernel must embed;
they are not silently called the complete physical cell/Schur kernel.

This is a complete kernel only for the declared formed-KEEP six-source
section.  It is not the complete cell response to all matter, EM, clock,
boundary, or geometric sources, and it does not discard the other formation
statuses from the full experimental instrument.

## 5. Operator-neutral next gate

GL6W closes

\[
 \boxed{
 \text{authenticated links}\to
 \text{physical six-source/read dilation}\to
 \text{normalized active-system two-time CTP response}\to
 \text{exact }A_1\oplus E_2\oplus T_2\text{ operator data}.} \tag{W31}
\]

Active-system overlap work may proceed in parallel, but no complete gluing or
infrared promotion is earned until the missing apparatus ledger and causal
Schur quotient in Section 4 are closed.  After that gate, the proof must show
whether the same physical query localizes records, clocks, matter, and EM into
one refining causal structure.  Conservation/Ward identities and the
infrared operator content must then be derived from the glued response.  Only
after those tests may
the response be projected as

\[
 H_2=\kappa_RM_R^\perp+H_{\rm res},                 \tag{W32}
\]

with Ricci/Einstein earned only if the complete residual and required
identities are controlled.  Equation (W32) is a downstream test, not an
assumption in this packet.

`PASS__EXACT_GL6V_NORMALIZED_SIX_SOURCE_CTP_CONVENTION__CLOSED_ALL_TIME_Ud0_RETARDED_ADVANCED_NOISE_KERNEL__ENTRANCE_MATCH_TO_GL6T_WITH_PHYSICAL_SIGN_AND_UNITS__EXACT_Ud_POSITIVE_FINITE_SPECTRAL_COMPILER__GL6U_INTERLINK_DEFECT_PROPAGATES_INTO_KERNEL__MATCHED_BREAK_ZERO__DIRECT_SYSTEM_SEAGULL_ZERO_ONLY__FULL_APPARATUS_SCHUR_OPEN__NO_METRIC_RICCI_GRAVITY_OR_G`

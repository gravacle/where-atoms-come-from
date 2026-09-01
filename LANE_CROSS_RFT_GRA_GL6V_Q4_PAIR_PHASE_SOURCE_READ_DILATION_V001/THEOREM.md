# q4 pair phase-source and complete-read dilation theorem

**Lane:** `LANE_CROSS_RFT_GRA_GL6V_Q4_PAIR_PHASE_SOURCE_READ_DILATION_V001`  
**Short name:** `GL6V V001`  
**Date:** 2026-08-31  
**Status:** author packet frozen after independent hostile pre-freeze review;
post-freeze custody audit pending  
**Claim class:** exact finite selected-apparatus composition; six independent
source-normalized pair pulse coordinates on the literal GL6T/GL6U links;
exact source-ancilla return; complete four-link/pair terminal read; exact CTP
equal-source normalization and direct affine-source contact census

**Not claimed:** that the probe is an autonomous F3 bulk law; that a pair
observable is a record; a continuous background metric source; a completed
interacting two-time CTP/contact/Schur kernel; same-query physical metric
calibration; cell gluing, common cone, Ward/Bianchi closure, infrared operator,
Ricci response, gravity, or `G`.

## 1. Exact remaining source problem

Audited GL6T and frozen GL6U use the six pair operators

\[
 M_A=Z_aZ_b,\qquad A=(ab)\in{\cal E}_4.             \tag{V01}
\]

Their entrance double commutators are exact, but GL6T correctly calls the
pair sources formal: no physical unit or complete branch schedule had been
attached to (V01).

Hostile-clean GL4P already supplies the necessary finite query circuit.  On
the exact `N=0` GL6T/GL6U star write `n_a=(1-Z_a)/2`.  Give each pair one
fresh ready bit `Q_A`, with `Z_Q|0>=|0>`, and define the garbage-free
star-local realization

\[
 \widetilde V_6:=
 \prod_{a<b}{\rm CNOT}_{n_a\to Q_{ab}}
              {\rm CNOT}_{n_b\to Q_{ab}}
 \equiv V_6^{(n,Q)}\otimes I_{K,G,{\rm source},{\rm spec}},             \tag{V02a}
\]

where the order is immaterial.  Each direct gate has the exact bounded form

\[
 {\rm CNOT}_{n\to Q}
 =\exp[-i\pi n(1-X_Q)/2].                              \tag{V02b}
\]

It has no transport, counter, work, or scratch register on which a link value
can escape.  The externally registered switching schedule, program, clock,
source-setting metadata, status, and ideal-model failure flag are event-blind
and are acted on by the displayed identity.  The failure flag has one
deterministic zero-failure outcome in this exact selected apparatus.  Thus
the full scheduled dilation is a deterministic unitary instrument, not a
success-filtered branch.

Equivalently on `(n,Q)`,

\[
 V_6=\sum_{q\in\{\pm1\}^4}\Pi_q\otimes
 \prod_{a<b}X_{Q_{ab}}^{(1-q_aq_b)/2},\qquad
 \Pi_q=\prod_{a=1}^4{I+q_aZ_a\over2}.              \tag{V02}
\]

GL4P proves that (V02) is unitary, involutive, and nondemolition for the
complete commuting `Z/pair` algebra.  Equations (V02a)--(V02b) select its
narrow garbage-free direct-gate realization on the literal links.  They are
a finite controlled-apparatus schedule, not a new autonomous F3 bulk term.

## 2. Six normalized source coordinates

After applying (V02), act on the six query bits with the registered phase
operation

\[
 P_Q(j)=\exp\!\left({i\over2}\sum_{A\in{\cal E}_4}j_AZ_{Q_A}\right),
 \qquad j_A\in\mathbb R,                            \tag{V03}
\]

and exactly uncompute the query and every possible event-bearing apparatus
factor:

\[
 U_{\rm src}(j)=\widetilde V_6^\dagger
 (P_Q(j)\otimes I_D)\widetilde V_6.                \tag{V04}
\]

The `j_A` are prospectively set external probe coordinates, not event labels,
record values, fitted bulk couplings, or new gravitational constants.  The
selected apparatus contains six independently registered phase controllers.
Their programs, pulse profiles, clocks, source-setting metadata,
boundary/status outputs, and deterministic ideal-model failure flag are
retained in the complete output and are event-blind.  During each
copy/phase/uncompute window, the system response Hamiltonian is switched off
or exactly refocused as in GL4P's controlled-apparatus schedule; outside the
window all probe couplings are off.

For a physical time profile, let `f_A(t)` be a frozen real registered profile
with units `1/time` and

\[
 \int_{I_A}f_A(t)\,dt=1.                            \tag{V05}
\]

The phase controller has Hamiltonian

\[
 H_{Q,\rm src}(t)=-{\hbar\over2}
 \sum_Aj_Af_A(t)Z_{Q_A}.                            \tag{V06}
\]

Thus `j_A` is a dimensionless pulse area.  Equivalently, for a square pulse
of registered duration `T_A`, its source energy is `E_A=hbar/T_A` and
`H=-(E_A/2)j_A Z_Q` on that interval.  This fixes the source units without
using an observed gravitational quantity.

## 3. Theorem `GL6V-1` -- exact phase-sandwich source

On the ready-ancilla subspace,

\[
 \boxed{
 U_{\rm src}(j)\bigl(|\psi\rangle\otimes|0^6\rangle_Q
                         \otimes|D_0\rangle\bigr)
 =\exp\!\left({i\over2}\sum_AM_Aj_A\right)|\psi\rangle
  \otimes|0^6\rangle_Q\otimes|D_0\rangle.}         \tag{V07}
\]

Hence (V04) is the exact finite dilation of the impulsive system source

\[
 \boxed{
 H_{\rm src}(t)=-{\hbar\over2}\sum_Aj_Af_A(t)M_A,} \tag{V08}
\]

with no escaped parity copy and no residual source term after the pulse.
Every source direction is independent: at the origin,

\[
 {\partial U_{\rm src}\over\partial j_A}\bigg|_0
 ={i\over2}M_A,\qquad
 2^{-4}\operatorname{Tr}(M_AM_B)=\delta_{AB}.      \tag{V09}
\]

### Proof

In a simultaneous `Z` eigenstate `|q>`, each pair of direct CNOTs in (V02a)
writes the pair bit
`(1-q_aq_b)/2`, so `Z_(Q_ab)` has eigenvalue `q_aq_b=M_ab`.  Equation (V03)
therefore contributes the phase
`exp[(i/2) sum_A j_A q_aq_b]`.  The reversed gate list returns all six
ancillas to ready without changing that phase.  Because (V02a) acts as
identity on `D`, every other retained apparatus factor remains in the same
common state `|D_0>`; there is no escaped parity copy.  Linearity proves
(V07).  Differentiation
gives (V09); Walsh-pair orthogonality on all sixteen four-bit strings gives
the trace identity.  QED.

The construction is exact for one pulse or any finite time-ordered list of
such pulses separated by the inherited system evolution.  Functional
derivatives at specified insertion times are therefore physical pulse
derivatives.  No claim is made that one finite sandwich equals simultaneous
continuous evolution under `H_parent+H_src` for an arbitrary profile without
the declared pulse/refocusing schedule.

## 4. CTP branch convention and direct contact

On the two CTP branches use separately frozen pulse areas `j_(+,A)` and
`j_(-,A)`, and define

\[
 j_{c,A}={j_{+,A}+j_{-,A}\over2},\qquad
 j_{\Delta,A}=j_{+,A}-j_{-,A}.                    \tag{V10}
\]

The finite closed-time functional uses the complete deterministic forward unitary
`U[j_+]`, the complete backward adjoint `U[j_-]^dagger`, and traces every
system and retained apparatus output.  For equal branch sources,

\[
 \boxed{Z[j,j]=\operatorname{Tr}(U[j]\rho U[j]^\dagger)=1.} \tag{V11}
\]

There is no success conditioning in (V11): every ideal-model failure outcome
other than the retained zero-failure flag has zero support.  Were a
nonunitary success branch selected instead, its equal-source trace would be
its success probability rather than one.  Equation (V11) is an exact
equal-source normalization, not an assumption that equal
CTP sources generate a nonzero response direction.  Query and source
profiles remain distinct when the later response kernel is contracted.

The physical source Hamiltonian (V08) is affine in every `j_A`, so its direct
system-source seagull is exactly

\[
 \boxed{
 {\partial^2H_{\rm src}\over\partial j_A\partial j_B}=0.} \tag{V12}
\]

Equation (V12) does not declare the complete physical metric/contact ledger
zero.  Time-ordering contacts, controller/work response, boundary terms,
source-profile limits, and any later geometric attachment remain explicit
owners in the full CTP calculation.

## 5. Theorem `GL6V-2` -- complete terminal read

After all response insertions, the terminal instrument retains the complete
sixteen-outcome computational query of the four literal link bits

\[
 q=(q_1,q_2,q_3,q_4)\in\{\pm1\}^4,                \tag{V13}
\]

together with `K/G`, controller, clock, work, status, failure, source-history,
boundary, and reference outputs.  Every pair result is then the deterministic
function

\[
 m_{ab}=q_aq_b.                                    \tag{V14}
\]

Alternatively, a terminal application of (V02) may retain all sixty-four
formal pair-bit outcomes, of which eight have nonzero support.  That parity-
only read has the exact global-sign degeneracy `q~ -q`; it is sufficient for
the pair algebra but is not substituted for (V13) in any future same-query
localization claim.

The phase-sandwich source itself leaves its six ancillas blank and unescaped,
so it does not dephase the response-ready links.  A terminal read is performed
only after the final CTP insertion.  Equations (V13)--(V14) are complete and
noncreating; impossible/failure outcomes remain registered rather than being
postselected away.

## 6. What is closed and the operator-neutral next gate

This packet earns the finite same-link chain

\[
 \boxed{
 \text{four authenticated support-gated links}
 \to\text{six independently normalized pair-source pulses}
 \to\text{source-off inherited response}
 \to\text{complete four-link and six-pair read}.}  \tag{V15}
\]

The phase controller is probe apparatus, not a force added to the F3 parent.
It supplies the missing operational derivatives with fixed units.  It neither
chooses the bulk values of `(h,Delta,U_d,d_star)` nor derives a gravitational
constant.

The immediate next theorem is to calculate, on this exact source/read parent,
the full retarded, advanced, symmetric/noise, direct-contact, controller,
boundary, reference, and nuisance response, then take the complete causal
Schur quotient.  Only afterward may cells be glued, a common cone and
Ward/Bianchi identities be tested, and the resulting infrared operator be
decomposed.  No Ricci template is inserted at the microscopic or q4 stage.

EW's exact `R^6 -> Sym^2(V)` Fisher isomorphism remains algebraically reusable,
but it is not yet a physical metric calibration for this coherent GL6T/GL6U
query law.  Clocks, matter, and electromagnetic probes must still be shown to
localize through the same complete query and fix its physical scale.

`PASS__HOSTILE_CLEAN_GL4P_SIX_PARITY_COPY_ATTACHED_TO_LITERAL_GL6T_GL6U_LINKS__EXACT_COPY_PHASE_UNCOPY_SOURCE_IDENTITY__SIX_INDEPENDENT_DIMENSIONLESS_HBAR_CLOCK_NORMALIZED_PAIR_PULSE_COORDINATES__SOURCE_ANCILLAS_RETURN_BLANK__CTP_EQUAL_SOURCE_NORMALIZATION__DIRECT_AFFINE_SOURCE_SEAGULL_ZERO_ONLY__COMPLETE_SIXTEEN_OUTCOME_LINK_AND_PAIR_READ__NO_BULK_FORCE_METRIC_RICCI_GRAVITY_OR_G`

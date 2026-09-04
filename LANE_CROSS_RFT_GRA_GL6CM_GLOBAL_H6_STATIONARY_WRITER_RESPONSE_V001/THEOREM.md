# GL6CM — GLOBAL H6 STATIONARY WRITER-RESPONSE THEOREM

## Status and exact scope

This packet converts the exact `GL6CH` source-dependent ring writer into a
same-state stationary response on any finite connected component of the
degree-two `H6` parent.  It proves reciprocity, positivity, the exact common
amplitude null, and a sharp kernel criterion.  Composing with `GL6CJ` gives a
parent-derived pair-source response; `GL6CK` supplies its first strict
two-overlap witness.

The result concerns the **spectral order-six block**.  It does not include the
source-before-Feshbach second derivative (contact), select a thermodynamic
phase, prove locality of the complete bulk kernel, identify a metric, derive
Ricci or Einstein dynamics, prove gravity, or calculate `G`.

## 1. Finite stationary component

Let `C` be a nontrivial connected component of a finite locked flip graph and
let `T_c` denote every distinct elementary cycle-toggle operator active on
that component.  After removal of the common scalar, the order-six
source-free locked-sector Hamiltonian is

\[
 H_0=-J\sum_cT_c,\qquad
 J={63\over8}{h^6\over U_d^5}>0.                       \tag{CM01}
\]

Perron--Frobenius gives a unique normalized component ground state `|0>`,
ground energy `E_0`, and a strictly positive finite-component gap.  Put

\[
 Q=1-|0\rangle\langle0|,\qquad
 R=Q(H_0-E_0)^{-1}Q.                                  \tag{CM02}
\]

The operator `R` is strictly positive on `Q`.

For a pure pair-tensor source, `GL6CH` gives the complete off-diagonal
source-linear order-six writer block

\[
 V_{\rm wr}^{(6)}(j)=\lambda_T\sum_c x_c(j)T_c .       \tag{CM03}
\]

\[
 \lambda_T={105\over16}{h^6\over U_d^6},\qquad
 x_c(j)=\sum_{v\in c}j_v^T\Theta_{v,c}.                \tag{CM04}
\]

Here `j` has energy units and `lambda_T` is dimensionless.  For a scalar
source coordinate `s`, write `j=s\widehat j` and

\[
 B_{\widehat j}:=\sum_cx_c(\widehat j)T_c,\qquad
 O_{\widehat j}:=\lambda_TB_{\widehat j}.              \tag{CM05}
\]

The theorem below deliberately uses only this writer block.  In particular,
`GL6CH` does not classify a possible diagonal first-source operator at order
six, and the full Feshbach Hamiltonian also has source-second contact terms.

## 2. Exact stationary writer response

Define the explicitly writer-only linear family

\[
 \widetilde H_{\widehat j}(s)=H_0+sO_{\widehat j},       \tag{CM05a}
\]

and let `E_wr(s)` be its analytic ground branch.  Ordinary nondegenerate
perturbation theory, applied to the same `H_0`, `|0>`, and writer operator,
gives

\[
 {d^2E_{\rm wr}(s)\over ds^2}\bigg|_{s=0}
 =-2\langle0|O_{\widehat j}R O_{\widehat j}|0\rangle. \tag{CM06}
\]

This is, by definition, the two-writer-vertex spectral contribution.  It is
not the second derivative of the full source-dependent Feshbach Hamiltonian:
a possible diagonal order-six first-source vertex, the order-`h^2`
source-second contact, and higher source-second vertices are outside this
calculation.

More generally define the bilinear spectral response

\[
 \boxed{
 {\cal K}_{\rm spec}(j,k)
 =2\lambda_T^2\operatorname{Re}
 \langle0|B_jR B_k|0\rangle .}                         \tag{CM07}
\]

Then

\[
 {\cal K}_{\rm spec}(j,j)\ge0,\qquad
 {\cal K}_{\rm spec}(j,k)={\cal K}_{\rm spec}(k,j),   \tag{CM08}
\]

and the kernel is exact:

\[
 \boxed{
 {\cal K}_{\rm spec}(j,j)=0
 \iff QB_j|0\rangle=0
 \iff B_j|0\rangle\in\operatorname{span}\{|0\rangle\}.} \tag{CM09}
\]

Thus positivity is not a mechanical stiffness postulate.  It is the spectral
theorem applied to the pair-source-sensitive transition operator generated
by the same microscopic parent.

In cycle-amplitude coordinates, (CM07) is the positive-semidefinite matrix

\[
 K_{cc'}=2\lambda_T^2\operatorname{Re}
 \langle0|T_cR T_{c'}|0\rangle.                        \tag{CM10}
\]

For a common amplitude variation `y_c=q` for every active cycle,

\[
 \sum_cqT_c=-{q\over J}H_0,                            \tag{CM11}
\]

so its action on `|0>` has no excited-state component.  Consequently

\[
 \boxed{K\mathbf1=0.}                                  \tag{CM12}
\]

This is an exact global common-rescaling null on every finite component.  It
does not imply that the kernel is only one-dimensional; (CM09) is the exact
test for additional dark combinations.

## 3. Pullback to the one microscopic pair-source probe

Let `W` be the finite-range writer map of `GL6CH/GL6CJ`,

\[
 (Wj)_c=x_c(j)=\sum_{v\in c}j_v^T\Theta_{v,c}.          \tag{CM13}
\]

The stationary pair-source block is the exact pullback

\[
 \boxed{K_T^{\rm spec}=W^TKW.}                         \tag{CM14}
\]

It is reciprocal and positive semidefinite without an orbit average or a
postulated material coefficient.  At each Q4 node the local writer has rank
three and normal `8P_T`, but this local rank does not by itself prove that the
global pullback has no phase-dependent kernel.  Equation (CM09) supplies the
missing global test.

Equations (CM07)--(CM14) use the same stationary state and the same parent
for every source leg.  They therefore replace the earlier nonstationary
fixed-background row as the lawful spectral assembly.  The full physical
Hessian must still add the owner-once contact
`<partial_j partial_j H_eff>` in this same state.

## 4. First strict accumulation witness

For one isolated ring the component is `K2`.  Every amplitude perturbation is
proportional to the source-free Hamiltonian, so (CM09) gives zero response.

For the `GL6CK` two-ring overlap star, let `w_0,w_1` be the two amplitude
derivatives in this writer-only linear family.  Exact diagonalization gives

\[
 -E_0''(0)={\sqrt2\over4J}(w_0-w_1)^2.                 \tag{CM15}
\]

This is a positive **difference law**: the common direction is null and the
relative direction is strict.  At one shared endpoint choose
`d j/ds=Theta_0`; the two geometric cycle tensors are orthogonal, so

\[
 w_0={105\over8}{h^6\over U_d^6},\qquad w_1=0,          \tag{CM16}
\]

and

\[
 \boxed{
 -E_0''(0)={175\sqrt2\over32}{h^6\over U_d^7}>0.}      \tag{CM17}
\]

An isolated writer therefore only changes a clock/amplitude, whereas the
first concrete overlap supplies a relative mode and a nonzero stationary
response.  This is the exact local sign and mechanism required by the
overlap-quadratic-form route.  It is not permission to add one copy of
(CM15) for every dense-parent edge: larger linked clusters can generate
cross terms and must be evaluated owner-once.

## 5. Effective-coordinate meaning

On any finite source subspace on which a chosen response block `K=-E''` is
strictly positive, the source-to-expectation map is locally invertible.  To
fix the sign convention, put `W=-E`, define the expectation coordinate by
`phi=dW/ds=-dE/ds`, and use
`Gamma(phi)=s phi-W(s)=E(s)+s phi`.  Its Legendre Hessian obeys

\[
 \Gamma''={\cal K}^{-1}>0.                              \tag{CM18}
\]

Thus a later quantity called field stiffness is not an extra microscopic
spring: it is the inverse of the derived stationary susceptibility.  For the
literal one-direction witness (CM17), the spectral-only inverse curvature is

\[
 ({\cal K}_{\rm spec})^{-1}
 ={16\sqrt2\over175}{U_d^7\over h^6}.                  \tag{CM19}
\]

Equations (CM18)--(CM19) do not yet make the probe source autonomous.  A
source-free composite-field action requires the full contact-plus-spectral
Hessian, a controlled stationary/thermodynamic domain, and real-time
response.  In particular, (CM19) is neither a Ricci coefficient nor `G`.

## 6. What has and has not closed

The microscopic chain now contains, in one parent:

\[
 \text{locked pair structure}
 \longrightarrow \text{six-direction dressed source access}
 \longrightarrow \text{future ring writer}
 \longrightarrow \text{strict overlap response}.       \tag{CM20}
\]

The remaining calculation is no longer whether a microscopic writer-sector
feedback channel exists.  It is whether the complete dense stationary
response—including record authentication, its contact, all linked clusters,
refinement, and causal continuation—has only the gravitational
long-wavelength structure.

`PASS__GL6CM_FINITE_COMPONENT_SAME_STATE_H6_WRITER_RESPONSE__RECIPROCAL_POSITIVE_SEMIDEFINITE_SPECTRAL_KERNEL__EXACT_QB_GROUND_KERNEL_CRITERION__COMMON_RING_RESCALING_NULL__SAME_PARENT_PAIR_SOURCE_PULLBACK_WTKW__ISOLATED_RING_ZERO__TWO_OVERLAP_STRICT_DIFFERENCE_RESPONSE_175SQRT2_OVER32_H6_UD7__INVERSE_RESPONSE_NOT_IMPORTED_STIFFNESS__RECORD_AUTH_CONTACT_BULK_PHASE_REALTIME_RICCI_GRAVITY_G_OPEN`

# GL6CQ — SAME-STATE STATIONARY-RESPONSE MOMENT SUM-RULE THEOREM

## Status and scope

This packet replaces the undetermined coefficients `kappa,b,c,d` of the
GL6CO stationary cycle symbol by exact zeroth- and second-moment projections
of one connected cycle susceptibility.  It also gives the `A1` coefficient
`alpha` and the mixed `A1-T2` coefficient `eta`.  Combining those moments
with the GL6BV contact probability `p`, evaluated in the **same stationary
state**, converts CO29--CO30 into directly measurable real-space sum rules.

The theorem concerns the explicitly defined contact-plus-two-writer spectral
sector.  It does not prove that a state satisfies either sum rule, select a
phase, establish convergence at criticality, complete every source-second
vertex, perform response-to-1PI inversion, derive Ricci or Einstein dynamics,
prove gravity, or calculate `G`.

## 1. One state and one bare cycle susceptibility

On a finite periodic connected component, let `|0>` be the unique
translation-, inversion-, and `S4`-invariant stationary ground state of the
same locked cycle Hamiltonian

\[
 H_0=-J\sum_cT_c,
 \qquad J={63\over8}{h^6\over U_d^5},                 \tag{CQ01}
\]

and put

\[
 Q_0=1-|0\rangle\langle0|,
 \qquad {\cal R}=Q_0(H_0-E_0)^{-1}Q_0.               \tag{CQ02}
\]

For a cycle of orientation `d` in cell `R`, define the **bare connected**
static susceptibility

\[
 \boxed{
 K^{\rm bare}_{dd'}(R)=
 2\,\operatorname{Re}\langle0|T_{0d}{\cal R}T_{Rd'}|0\rangle .} \tag{CQ03}
\]

It contains neither `lambda_T^2` nor `mu^2`.  Equivalently,

\[
 K^{\rm bare}_{dd'}(R)=2\int_0^\infty d\tau\,
 \operatorname{Re}\langle\delta T_{0d}(\tau)
                         \delta T_{Rd'}(0)\rangle_0,   \tag{CQ04}
\]

whenever the integral converges.  Thus the inputs below are connected
Euclidean correlation observables, not fitted material constants.

For the bulk formulas, take a sequence of such periodic components and
assume that the connected correlator has a translation-covariant
thermodynamic limit.  Equations below are first exact finite-volume spectral
identities; their real-space moment form refers to that limit.  This avoids
treating a wrapped finite-torus displacement as a physical second moment.

Let `chi_d` be the inherited GL6CL center of orientation `d`, let `L_R` be
the inherited Bravais displacement, and define the centered separation

\[
 X_{R;dd'}=L_R+\chi_{d'}-\chi_d.                       \tag{CQ05}
\]

The centered Fourier symbol is

\[
 \widehat K_{dd'}(k)=\sum_R
 e^{ik\cdot X_{R;dd'}}K^{\rm bare}_{dd'}(R).           \tag{CQ06}
\]

These inherited coordinates remain a calculation chart; no physical length
calibration is asserted.

## 2. Convergent real-space moments

Assume entrywise absolute second-moment convergence,

\[
 \sum_R(1+|X_{R;dd'}|^2)
       |K^{\rm bare}_{dd'}(R)|<\infty .                \tag{CQ07}
\]

Define

\[
 Z_{dd'}=\sum_RK^{\rm bare}_{dd'}(R),
 \qquad
 M_{dd'}^{mn}=\sum_RX_{R;dd'}^mX_{R;dd'}^n
 K^{\rm bare}_{dd'}(R).                               \tag{CQ08}
\]

Reciprocity and inversion cancel the odd moment, so

\[
 \boxed{\widehat K(k)=Z-{1\over2}k_mk_nM^{mn}+o(|k|^2).} \tag{CQ09}
\]

Finite fourth absolute moment upgrades the remainder to `O(|k|^4)`; an
exponential moment gives an analytic symbol near zero.

Use the orthogonal cycle solder

\[
 u_d={1\over2},\qquad Q_{di}={1\over2}(T_d)_i,
 \qquad S=(u\;Q),                                     \tag{CQ10}
\]

and denote the three columns of `Q` by `Q_i`.  Project the moments as

\[
\begin{aligned}
 M_{AA}^{mn}&=u^TM^{mn}u,\\
 M_{AT,i}^{mn}&=u^TM^{mn}Q_i,\\
 M_{TT,ij}^{mn}&=Q_i^TM^{mn}Q_j.                       \tag{CQ11}
\end{aligned}
\]

The common-amplitude spectral null gives `Z u=0`.  `S4` invariance then
implies `Q^T ZQ=kappa I`.  Define the observable contractions

\[
\begin{aligned}
 Z_T&=\sum_iQ_i^TZQ_i,\\
 M_\perp&=\sum_{i\ne m}M_{TT,ii}^{mm},\\
 M_\parallel&=\sum_iM_{TT,ii}^{ii},\\
 M_\times&=\sum_{i\ne j}M_{TT,ij}^{ij}.               \tag{CQ12}
\end{aligned}
\]

## 3. Exact coefficient projections

For the GL6CO invariant form

\[
 S^T\widehat K(k)S=
 \begin{pmatrix}
 \alpha |k|^2&\eta(k_yk_z,k_zk_x,k_xk_y)\\
 \eta(k_yk_z,k_zk_x,k_xk_y)^T&
 \kappa I+b|k|^2I+cD(k)+dO(k)
 \end{pmatrix}+o(|k|^2),                              \tag{CQ13}
\]

the coefficients are not free placeholders.  They are exactly

\[
\boxed{
 \kappa={Z_T\over3},\qquad
 \alpha=-{1\over6}\sum_mM_{AA}^{mm},}                 \tag{CQ14}
\]

\[
\boxed{
 \eta=-{1\over3}\left(
 M_{AT,x}^{yz}+M_{AT,y}^{zx}+M_{AT,z}^{xy}\right),}    \tag{CQ15}
\]

and

\[
\boxed{
 b=-{M_\perp\over12},\qquad
 c=-{M_\parallel\over6}+{M_\perp\over12},\qquad
 d=-{M_\times\over6}.}                               \tag{CQ16}
\]

Equivalently,

\[
 Z_T=3\kappa,\quad M_\perp=-12b,\quad
 M_\parallel=-6(b+c),\quad M_\times=-6d.             \tag{CQ17}
\]

The mixed coefficient `eta` and scalar coefficient `alpha` do not enter the
common tensor pullback through order `k^2`: the zero-mode writer lands in
cycle `T2`, while its first cycle-`A1` component is already quadratic.
They are nevertheless included because they are observable moments needed
for a later full `A1+E2+T2` response audit.

## 4. Normalization custody: the writer appears once

GL6CM uses

\[
 \lambda_T={105\over16}{h^6\over U_d^6}
\]

with the cycle score `j.Theta`.  GL6CL uses

\[
 \mu={105\over8}{h^6\over U_d^6}=2\lambda_T           \tag{CQ18}
\]

with canonical pair coordinates.  For a pure `T2` source,
`j.Theta_ab=2j_ab`, so

\[
 \lambda_T(j\cdot\Theta_{ab})=\mu j_{ab}.             \tag{CQ19}
\]

These are the same vertex in two bases, not two factors to multiply.
Because (CQ03) is bare, the spectral block in the GL6CL common coordinate
`j_+=(j_P+j_C)/2` is

\[
 \boxed{{\cal H}^{H6}_{T,j_+}(k)=
 \mu^2B_T(k)^*\widehat K(k)B_T(k).}                   \tag{CQ20}
\]

whereas the orthonormal common source used by the GL6BV contact is
`widehat j_+=(j_P+j_C)/sqrt(2)=sqrt(2)j_+`.  Therefore

\[
 \boxed{\widehat{\cal H}^{H6}_T(k)=
 {\mu^2\over2}B_T(k)^*\widehat K(k)B_T(k).}           \tag{CQ20a}
\]

The half is a coordinate normalization, not another writer factor; `mu^2`
is still applied exactly once.  Independent exact composition in the GL6CL
coordinate reproduces

\[
 {1\over\mu^2}{\cal H}^{H6}_{T,j_+}=
 8\kappa I+(-2\kappa+8b)|k|^2I
 +(-16\kappa+8c)D+(12\kappa+8d)O+o(|k|^2).           \tag{CQ21}
\]

## 5. The same-state contact observable

In the same state `|0>`, let

\[
 p=\langle0|\Pi^{(v,w;a)}_{\rm same}|0\rangle
  =\Pr_0[\sigma_{z_v}(a)=\sigma_{z_w}(a)]             \tag{CQ22}
\]

for an edge of port `a`.  Translation and `S4` invariance make it independent
of the selected edge.  This is the equal-time expectation of the local
same-partner projector, not a value imported from another orbit mixture.

With

\[
 g_{\rm ct}={h^2\over4U_d^3},                          \tag{CQ23}
\]

the GL6BV common contact contributes the quadratic coefficients

\[
 A_{\rm ct}={4\over3}(1-4p),\qquad
 B_{\rm ct}=0,
 \qquad C_{\rm ct}={8\over3}(2p-1).                   \tag{CQ24}
\]

No `mu` multiplies this contact.  It is the adiabatically eliminated form of
the high-gap one-link Kubo history; adding that high-gap history again would
double count it.

## 6. Exact observable sum rules

For the defined same-state contact-plus-two-writer spectral Hessian, the
necessary-and-sufficient `T2` **quadratic-gradient** symmetric-tensor
extension condition is

\[
 {\mu^2\over2}[-4\kappa+8(c+d)]
 +g_{\rm ct}{8\over3}(2p-1)=0.                         \tag{CQ25}
\]

Substituting (CQ14)--(CQ17) and multiplying by `3/2` gives the exact
real-space observable form

\[
\boxed{
 {\mu^2\over2}[-2Z_T+M_\perp-2M_\parallel-2M_\times]
 +4g_{\rm ct}(2p-1)=0.}                               \tag{CQ26}
\]

The additional held-out reference-shape condition is

\[
 {\mu^2\over2}[-2\kappa+8b]
 +g_{\rm ct}{4\over3}(1-4p)=0,                        \tag{CQ27}
\]

or equivalently

\[
\boxed{
 -{\mu^2\over2}[Z_T+M_\perp]+2g_{\rm ct}(1-4p)=0.}  \tag{CQ28}
\]

Thus the earlier coefficient question has become a concrete measurement:
integrate the connected cycle correlations with weights `1` and `X_mX_n`,
measure `p` in that same state, and test (CQ26).  Testing the stronger
reference shape requires both (CQ26) and (CQ28).

Neither test removes or interprets the zero-momentum tensor term.  In
particular, these equations do not prove background stationarity,
masslessness, or a gauge null.

Equations (CQ26)--(CQ28) assert exact equivalences for the displayed
two-block Hessian; they do not assert that nature or a selected phase makes
their left-hand sides vanish.  If a later complete source-first calculation
finds an additional source-second block at the same order, its independently
measured moment projection must be added rather than hidden inside these
terms.

## 7. Thermodynamic and critical limits

On every finite connected component the resolvent and sums are finite, but a
physical infinite-volume coefficient requires a controlled sequence of
stationary states and convergence of (CQ07).  A spectral gap plus a proved
clustering bound is sufficient; finite-volume Perron--Frobenius uniqueness
alone is not.

If the zeroth moment converges but the second moment diverges, `kappa` may
exist while `b,c,d` do not.  At a critical point the symbol can instead have
nonanalytic behavior such as `|k|^sigma` or `k^2 log|k|`.  In that case the
analytic CO quadratic classification and (CQ26)--(CQ28) are not the correct
test; the nonanalytic leading kernel must be measured and matched directly.
Moment divergence is therefore a physical criticality diagnostic, not a
license to assign infinite coefficients.

## 8. Exact disposition

Established:

1. every GL6CO stationary coefficient through `k^2` is an explicit moment of
   the bare connected cycle susceptibility;
2. the `CM` and `CL` writer conventions agree through `mu=2 lambda_T`, and
   the writer is applied exactly once;
3. the GL6BV contact and cycle moments are lawfully combined in one state;
4. CO29 and CO30 are now the observable real-space sum rules (CQ26) and
   (CQ28); and
5. their precise convergence and critical-nonanalyticity boundary is known.

Still open are evaluation of the moments in a selected thermodynamic state,
proof that either sum rule vanishes, completion of any further source-second
blocks, full `E2-T2/E2-E2` response, response-to-1PI inversion, record
authentication, refinement, causal continuation, Ricci/Einstein dynamics,
gravity, and `G`.

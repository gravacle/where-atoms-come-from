# Direct authenticated-pair finite-mission collar theorem

**Lane:** `LANE_CROSS_RFT_GRA_GL6BA_AUTHENTICATED_PAIR_FINITE_MISSION_COLLAR_V001`  
**Short name:** `GL6BA V001`  
**Date:** 2026-08-31  
**Status:** author completion candidate; distinct hostile audit required before
promotion  
**Claim class:** exact finite-collar reduction of the full same-parent F3
dynamics for an authenticated binary pair mission; exact `A_3` shell
refinement of the inherited boundary-Cauchy estimate; state-independent
all-finite-duration error; direct treatment of the admitted `R=2,5/2`
members without the `GL6AY` high-ratio hypothesis

**Disposition:**
`PASS__FULL_F3_NOT_EFFECTIVE_LOCKED_GENERATOR__EXACT_INTERACTION_CLUSTER_COLLAR__A3_SHELL_10L2_PLUS2__UNIFORM_OUTSIDE_VOLUME_OPERATOR_TAIL__AUTHENTICATED_BINARY_PAIR_DTV_TAIL__ALL_FINITE_SIGMA_AND_ALL_FINITE_R__R2_AND_R5_OVER2_DIRECTLY_LICENSED__ONLY_MEMBER_CLOCK_AND_COLLAR_STATE_PAYLOAD_OPEN__NO_GRAVITON_RICCI_GRAVITY_OR_G`

**Not claimed:** selection of a physical value of `R`; a numerical mission
duration; that the collar bound is numerically sharp; one fixed collar valid
for infinite time; total variation control of all authentication/status
flags; an exact locked phase; a selected state; gravity or `G`.

## 1. The shortest lawful branch after `GL6AZ`

`GL6AZ` proves that the admitted ratios `R=2` and `R=5/2` lie outside the
very conservative sufficient domain of the `GL6AY` prethermal normal-form
theorem.  That is a boundary of that proof route, not a failure of either F3
member.  The present packet takes the other branch named by `GL6AZ`: calculate
the finite authenticated mission directly under the full F3 Hamiltonian.

On the selected ideal all-formed/`MATCH` active factor, the exact parent is

\[
 H=U_dN_{\rm def}-h\sum_pX_p,
 \qquad N_{\rm def}=\sum_v(k_v-2)^2,
 \qquad U_d,h>0.                                             \tag{BA01}
\]

Put

\[
 R={U_d\over h},\qquad
 s={h(t-t_F)\over\hbar},\qquad
 \sigma_{\rm obs}={h(t_Q-t_F)\over\hbar}.                   \tag{BA02}
\]

On the homogeneous authenticated `A_3` incidence parent, expansion of
`N_def` gives, up to a scalar that drops out of dynamics, the dimensionless
quasi-local interaction

\[
 \boxed{
 \overline H_R={H\over h}
 =\sum_{p\in\mathbb L}(-X_p-6Rn_p)
  +2R\sum_{\{p,q\}\in E_\infty}n_pn_q,
 \qquad n_p={1-Z_p\over2}.}                                  \tag{BA03}
\]

Every link-graph site has exactly six partners.  The onsite product dynamics
in (BA03) is treated exactly and never enlarges support.  In its interaction
picture each pair term retains two-link support and norm at most `2R`; hence
the inherited dimensionless dressing-complete influence rate is

\[
 \lambda_R=48R.                                               \tag{BA04}
\]

No locked-space projection, effective Hamiltonian, large-`R` expansion, or
prethermal hypothesis is used below.

## 2. Exact causal collars and their size

Let the authenticated pair label be `beta=(0,{a,b})` and set

\[
 M_\beta=Z_{(0,a)}Z_{(0,b)},\qquad
 X_\beta=\{(0,a),(0,b)\}.                                    \tag{BA05}
\]

For an integer `m>=0`, the exact interaction-cluster collar is

\[
 C_m(\beta)=\{p\in\mathbb L:
       d_L(p,X_\beta)\le m\},                                \tag{BA06}
\]

where `d_L` is distance in the physical link-interaction graph.  A
time-ordered interaction-picture word with at most `m` pair-interaction
insertions cannot leave (BA06): each insertion meets the already reached
support and crosses at most one link-graph edge.  Thus every Dyson/Krylov
cluster through pair-interaction order `m` is exactly reproduced on
`C_m(beta)`.  The statement treats all onsite rotations exactly; `m` is not
silently renamed a physical time.

For one finite authenticated FPSS realization, use the centered cell collar

\[
 B_L=\{x\in A_3:d_{A_3}(x,0)\le L\},\qquad
 \Lambda_L=B_L\times\{1,2,3,4\}.                             \tag{BA07}
\]

A link-graph step is either same-parent (`P`) or same-child (`C`).  A maximal
consecutive `C` run telescopes:

\[
 (x,a_0)\xrightarrow C\cdots\xrightarrow C
 (x+e_{a_0}-e_{a_j},a_j).
\]

It changes the cell by at most one `A_3` step, and two distinct `C` runs are
separated by at least one `P` step.  Hence

\[
 d_L(p,q)\ge2d_{A_3}(\operatorname{parent}p,
                     \operatorname{parent}q)-1
 \quad\hbox{for distinct parent cells},                      \tag{BA08}
\]

and therefore

\[
 C_m(\beta)\subseteq\Lambda_{\lceil m/2\rceil}.              \tag{BA08a}
\]

Alternating `C,P,C,...,C` geodesics attain equality.  Thus
`B_ceil(m/2)` is the smallest centered **cell ball** guaranteed by geometry
to contain every pair-interaction Dyson/Krylov cluster through order `m`.
This does not claim that every observable or parameter choice saturates
every allowed cluster.

For the ordinary-time Taylor series there is a further exact F3 filtration.
Initially `M_beta` is diagonal, every pair term is diagonal, and a newly
reached link enters diagonally.  Before each pair edge can extend support, a
transverse onsite `X` commutator is required at its inside endpoint.  An
omitted interaction across `B_L` therefore first becomes possible only at
ordinary nested-commutator order `4L+2`.  Consequently `B_L` exactly
reproduces all ordinary-time Taylor coefficients through order `4L+1`; for a
requested order `m>=1`, the smallest geometry/commutator-licensed centered
cell radius is `max(0,ceil((m-1)/4))`.  This coefficient statement is not used
to turn a finite Taylor polynomial into an uncontrolled finite-time answer.
Alternating transverse onsite and diagonal pair commutators along an
attaining geodesic gives a nonzero formal boundary-reaching word at order
`4L+2`, so no smaller radius is guaranteed by the F3 term algebra alone.

The exact `A_3` shell and ball counts are

\[
 |B_r\setminus B_{r-1}|=10r^2+2\quad(r\ge1),                 \tag{BA09}
\]

\[
 \boxed{|B_L|={10L^3+15L^2+11L+3\over3},\qquad
 |\Lambda_L|=4|B_L|.}                                       \tag{BA10}
\]

Every coordinate of a point in `B_L` lies in `[-L,L]`.  Therefore the
translation

\[
 x\longmapsto(L+1,L+1,L+1,L+1)+x,\qquad N=4(L+1),           \tag{BA10a}
\]

embeds the complete collar strictly inside one finite FPSS slab while
preserving every port and literal shared-child relation.  More generally,
every finite inherited open exterior used below has such a finite
authenticated FPSS ancestry/embedding.  The induced open cut at
`partial(Lambda_L)` is a mathematical restriction used to certify a
calculation; it is not asserted to be an independently performed physical
switch or a separately authenticated mission.

To prove (BA09), classify a nonzero `x in A_3` by `p` positive and `q`
negative coordinates.  Its distance is the common positive/negative sum
`r`, and the exact census is

\[
 \sum_{\substack{p,q\ge1\\p+q\le4}}
 {4\choose p}{4-p\choose q}
 {r-1\choose p-1}{r-1\choose q-1}=10r^2+2.                  \tag{BA11}
\]

Thus the physical finite calculation is explicitly sized before any matrix,
Krylov, or tensor-network evaluation is attempted.

## 3. Uniform finite-parent collar error

Let `mathfrak F_L` be the family of complete finite all-formed/`MATCH` FPSS
missions whose authenticated translation places `Lambda_L` strictly in the
interior of the mission's full active-link graph `Omega`.  On that graph the
exact finite Hamiltonian (BA01) is the restriction of (BA03): the inherited
onsite coefficient is `-6R`, and every inherited parent/shared-child pair
present in the complete FPSS graph is retained.  Let `tau_s^(R,Omega)` be
this finite-dimensional physical dynamics.  Let
`tau_s^(R,L)` retain precisely the onsite terms on `Lambda_L` and the pair
terms whose two endpoints lie in `Lambda_L`, extended by the identity on the
rest of `Omega`.  Equivalently one may retain the product onsite dynamics on
all of `Omega`; outside onsite factors commute with the collared observable
and give the same answer.  No wraparound or fitted boundary term is added.
The induced cut at `partial(Lambda_L)` is only a proof device; `Omega`, not
the cut, is the complete authenticated finite mission.

`GL6AK` separately proves that the net of these finite inherited open
dynamics has a boundary-independent quasi-local limit, denoted `tau_s^R`.
That limit will be used only as a mathematical corollary after the finite
authenticated statement has been proved.

Only shared-child pair terms cross the cell-collar boundary.  Their exact
number is

\[
 \boxed{C_L=36L^2+36L+12=12(3L^2+3L+1).}                    \tag{BA12}
\]

Indeed, fix one of the twelve directed `A_3` displacements `e_a-e_b`.  A cell
`x` exits `B_L` along that displacement precisely when it lies on the radius
`L` shell with `x_a>=0` and `x_b<=0`.  The bivariate positive/negative-sum
generating function for those cells is

\[
 { (1-zw)^2\over(1-z)^3(1-w)^3}.
\]

Its `z^Lw^L` coefficient is

\[
 {L+2\choose2}^2-2{L+1\choose2}^2+{L\choose2}^2
 =3L^2+3L+1.                                                 \tag{BA13}
\]

Multiplication by twelve proves (BA12), with no coarse volume-shell sum.

Split the complete finite `Omega` interaction into the collared Hamiltonian,
the exterior
Hamiltonian, and the `C_L` physical cross-boundary pair terms.  Interior and
exterior parts act on disjoint link factors, so the decoupled evolution of
`M_beta` is exactly `tau_s^(R,L)(M_beta)`.  For a crossing pair
`Z={p,q}`, with `p` inside and `q` outside,

\[
 [2Rn_pn_q,\tau_u^{(R,L)}(M_\beta)]
 =2R[n_p,\tau_u^{(R,L)}(M_\beta)]n_q.                        \tag{BA14}
\]

Thus only the inside endpoint enters.  The inherited arbitrary-support
commutator theorem gives the following bound.  An induced-collar link
distance can only exceed its infinite-graph value, so writing the latter as
`d_L` is conservative:

\[
 \|[n_p,\tau_u^{(R,L)}(M_\beta)]\|
 \le2\sum_{z\in X_\beta}T_{d_L(z,p)}(48R|u|),
 \qquad T_d(x)=\sum_{k=d}^{\infty}{x^k\over k!}.              \tag{BA15}
\]

Every crossing edge appends an outside endpoint at cell radius `L+1`.
Equation (BA08) applied to that endpoint and removal of the last link step
give `d_L(z,p)>=2L` for its inside endpoint.  The exact boundary-only Duhamel
estimate is therefore

\[
\boxed{
\begin{aligned}
 &\sup_{\Omega\in\mathfrak F_L}
\|\tau_s^{(R,\Omega)}(M_\beta)
       -\tau_s^{(R,L)}(M_\beta)\|\\
 &\quad\le {1\over12}
 \sum_{Z=(p,q)\in\partial\Lambda_L}
 \sum_{z\in X_\beta}
 T_{d_L(z,p)+1}(48R|s|)\\
 &\quad\le
 2(3L^2+3L+1)T_{2L+1}(48R|s|)
 =:\varepsilon_L^{\partial}(R,s).
 \end{aligned}}                                               \tag{BA16}
\]

The coefficient `1/12` is not fitted: Duhamel supplies `2R`, the commutator
bound supplies `2`, and integration supplies `1/(48R)`.  The coarse second
line then uses two initial links and the exact `C_L` boundary census.
Equation (BA16) is uniform in every complete authenticated finite exterior
in `mathfrak F_L`.  The same operator estimate also holds for an arbitrary
mathematical induced open `Omega` containing `Lambda_L`: such a cut contains
at most the `C_L` crossing terms summed in (BA16), but the cut itself is not
thereby promoted to an authenticated physical mission.

Taking the `GL6AK` norm limit proves the identical estimate with
`tau_s^(R,Omega)` replaced by the mathematical quasi-local dynamics
`tau_s^R`.  Combining either version with the trivial operator ceiling gives

\[
\boxed{
 \sup_{\Omega\in\mathfrak F_L}
\|\tau_s^{(R,\Omega)}(M_\beta)
       -\tau_s^{(R,L)}(M_\beta)\|
\le\min\{2,\varepsilon_L^{\partial}(R,s)\}.}                \tag{BA17}
\]

For `R=0`, distinct links factorize and the collar error is exactly zero.
For every finite `R>0` and finite `s`, the factorial tail dominates the
quadratic boundary count, so

\[
 \lim_{L\to\infty}\varepsilon_L^{\partial}(R,s)=0           \tag{BA18}
\]

uniformly for `s` in compact intervals.  This is an all-finite-duration
certificate, not a claim that one finite collar controls infinite time.

## 4. Certified collar radius

For a prospectively declared binary-probability tolerance
`delta in (0,1)`, define

\[
 \boxed{
L_{\rm cert}(R,|s|,\delta)=
 \min\left\{L\ge0:
(3L^2+3L+1)T_{2L+1}(48R|s|)\le\delta\right\}.}             \tag{BA19}
\]

Equation (BA18) proves that (BA19) is finite for every finite input triple.
`L_cert` is the smallest collar licensed by this explicit conservative
bound; it need not be the smallest collar attainable after exact
cancellation or a sharper model-specific estimate.

There is also a closed analytic sizing envelope.  For every `mu>0`,

\[
 T_d(x)\le e^{xe^\mu-\mu d}.                                 \tag{BA19a}
\]

If `0<x<d=2L+1`, choosing `mu=ln(d/x)` gives

\[
 (3L^2+3L+1)T_{2L+1}(x)
 \le(3L^2+3L+1)\left({ex\over2L+1}\right)^{2L+1}.           \tag{BA19b}
\]

This proves directly that a finite certified collar exists with radius
scaling no worse than linearly in `R|s|` plus the accuracy demand.  It is an
upper-envelope statement, not a measured propagation speed.

At `x=0`, the exact tail is `T_d(0)=0` for every `d>=1`; no logarithmic
optimization is invoked.

## 5. The finite authenticated binary result

Fix one of the complete finite authenticated exterior missions
`Omega in mathfrak F_L` above.  Let
`omega_Omega` be its common postformation state on the selected active factor,
and let `omega_L` be the exact reduction of that same state to `Lambda_L`.
Immediately after the system sampling endpoint `t_Q`, perform the same
complete flag-retaining terminal read with the response off or refocused.  Its
binary pair marginal has projectors

\[
 P_\pm={I\pm M_\beta\over2}.                                 \tag{BA20}
\]

Define the finite-parent and collared binary results by

\[
 p_\pm^\Omega={1\over2}
 [1\pm\omega_\Omega(\tau_{\sigma_{\rm obs}}^{(R,\Omega)}(M_\beta))],\qquad
p_\pm^{(L)}={1\over2}
[1\pm\omega_L(\tau_{\sigma_{\rm obs}}^{(R,L)}(M_\beta))]. \tag{BA21}
\]

Binary total variation and (BA16)--(BA17) give the direct finite-mission
theorem

\[
\boxed{
 D_{\rm TV}(p^\Omega,p^{(L)})
\le\min\{1,\,(3L^2+3L+1)
T_{2L+1}(48R|\sigma_{\rm obs}|)\}.}                        \tag{BA22}
\]

All terminal authentication, route, clock, status, and failure outputs remain
retained.  Equation (BA22) controls only the selected-factor binary pair
marginal; it neither postselects `MATCH` nor bounds the joint distribution of
all flags.

If `omega` is a quasi-local state compatible with an exhaustion of these
finite parents, define `p^infinity` from `omega(tau^R(M_beta))`.  The `GL6AK`
norm limit and (BA22) give the same right-hand side for
`D_TV(p^infinity,p^(L))`.  This is a mathematical completion/extension of the
finite authenticated theorem, not a claim that one infinite record or one
infinite authenticated mission was physically realized.

The finite quantity on the right of (BA21) can be evaluated by exact finite
exponentiation, Krylov propagation, or a separately controlled tensor method
on the `4|B_L|` active qubits.  The collar theorem controls the omitted
physical exterior independently of the numerical method used inside.

## 6. Direct status of the admitted moderate ratios

No condition in (BA16)--(BA22) requires `R` to be large.  For the two already
admitted H6 member witnesses,

\[
 \begin{array}{c|c|c}
 h/U_d&R=U_d/h&48R|\sigma_{\rm obs}|\\ \hline
 1/2&2&96|\sigma_{\rm obs}|\\
 2/5&5/2&120|\sigma_{\rm obs}|
 \end{array}                                                  \tag{BA23}
\]

and (BA22) applies without extrapolating `GL6AY`.  This proves that both
members are directly treatable for every finite mission duration.  It does
not select either member as nature's value or assert that the resulting
collar is computationally small.

The shortest physical payload remains the one isolated by `GL6AZ`:

\[
 R={\Delta_{\rm def}\over2A_X},\qquad
 \sigma_{\rm obs}={A_X(t_Q-t_F)\over\hbar},                  \tag{BA24}
\]

with `A_X=h`, `Delta_def=2U_d`, a formation-completion/sampling timestamp
pair on the same clock, and custody that (BA01) governs that interval.  To
produce a numerical pair probability rather than only a state-independent
exterior error, the selected mission must additionally provide the reduced
postformation state `omega_L`; the prepared-blank mission is one already
defined special case.  The tolerance `delta` is a prospective numerical
accuracy specification, not a new physical constant.

## 7. Strict ceilings

1. `GL6BA` bypasses the unusable `GL6AY` sufficient regime; it does not repair
   or extrapolate that normal-form theorem.
2. The comparison is full F3 dynamics versus a spatial collar of the same
   full F3 dynamics, not full dynamics versus a locked/effective generator.
3. `L_cert` is sufficient and explicit, not asserted numerically optimal.
4. A fixed collar is not uniform for unbounded time; every declared finite
   mission has a finite certified collar.
5. No numerical `R`, `sigma_obs`, or postformation state is manufactured.
6. Only the authenticated binary pair marginal is bounded, not the full flag
   distribution.
7. No selected GNS state, pole, continuum, graviton, Ricci form, Einstein
   equation, gravity identification, or Newton constant is claimed.

No graviton, Ricci target, Einstein equation, gravity identification, or `G`
is used as a premise or inferred as a conclusion.

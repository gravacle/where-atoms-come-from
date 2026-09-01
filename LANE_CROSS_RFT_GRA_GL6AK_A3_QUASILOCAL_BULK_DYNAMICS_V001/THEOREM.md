# A3 quasi-local bulk-dynamics and stationary spectral-measure theorem

**Lane:** `LANE_CROSS_RFT_GRA_GL6AK_A3_QUASILOCAL_BULK_DYNAMICS_V001`  
**Short name:** `GL6AK V001`  
**Date:** 2026-08-31  
**Status:** author frozen after independent hostile pre-freeze review;
post-freeze custody/replay audit required  
**Claim class:** exact infinite locally finite completion of the authenticated
`A3` shared-child atlas; exact thermodynamic dynamics for the inherited
all-formed F3 interaction; explicit exhaustion-boundary error series; exact
translation and `S4` covariance; existence of joint stationary bulk states;
positive matrix-valued stationary spectral measures for the six local pair
observables

**Not claimed:** one infinite terminal record query; autonomous selection of
the all-formed route or a stationary state; convergence of finite-volume
ground or Gibbs states; uniqueness, clustering, a gap, a pole, hydrodynamics,
or a propagating mode; physical length, momentum, or speed; a Lorentz/common
physical cone; a continuum, Ward/Bianchi closure, Ricci/Einstein form,
gravity, or `G`.

## 1. The exact question and selected parent

`GL6AA` supplies, on every finite FPSS slab, independently queried parent
cell IDs, four port labels, literal shared children, the displacement
`e_a-e_b`, and its complete MATCH/MISMATCH/BLANK/COLLISION instrument.
`GL6AI` supplies the inherited all-formed link Hamiltonian in the exact form

\[
 H_{N}=\sum_e[-hX_e+\varepsilon_\star n_e]
       +\sum_{\{e,f\}\in E(L_N)}2U_dn_en_f+C_N,
 \qquad
 \varepsilon_\star=\Delta+2U_d(1-2d_\star),                 \tag{AK01}
\]

and proves a size-independent quasi-local influence tail with

\[
 J=2|U_d|,\qquad \Delta_L=6,\qquad
 \lambda_{\rm F3}={4J\Delta_L\over\hbar}
                  ={48|U_d|\over\hbar}.                      \tag{AK02}
\]

The present packet asks whether those finite authenticated neighborhoods
define one boundary-independent bulk dynamics and whether that dynamics has
stationary response measures before any pole, metric, or Ricci form is
assumed.  The answer is yes for the selected homogeneous all-formed member.
The selection of that member is an input; the limit, covariance, stationary
state existence, and spectral-measure statements are derived below.

## 2. Infinite authenticated incidence pattern

Let

\[
 A_3:=\{x\in\mathbb Z^4:\mathbf1^Tx=0\},\qquad
 \mathbb L:=A_3\times\{1,2,3,4\}.                            \tag{AK03}
\]

The site `(x,a)` is the active link leaving parent cell `x` through port
`a`.  Define two kinds of unordered pair-interaction edge:

\[
 \begin{aligned}
 I(x;a,b)&=\{(x,a),(x,b)\}, &&a<b,\\
 S(x;a,b)&=\{(x,a),(x+e_a-e_b,b)\}, &&a<b.
 \end{aligned}                                                \tag{AK04}
\]

Write `E_infty` for all edges in (AK04).  Every site `(x,a)` has precisely
the three other ports in its parent cell and the three ports

\[
 (x+e_a-e_b,b),\qquad b\ne a,                                \tag{AK05}
\]

with the same literal child.  Hence

\[
 \boxed{\deg_{E_\infty}(x,a)=3+3=6.}                          \tag{AK06}
\]

This is the homogeneous interior value of the exact `GL6AI` census, not an
inserted regular grid.

The infinite incidence pattern has finite authenticated ancestry in the
following precise sense.  For any finite set \(F\subset A_3\)—enlarged first
to include any finite collar required by the calculation—choose

\[
 m_i\ge 1-\min_{x\in F}x_i,
 \qquad m_i\ge1,
 \qquad N=\sum_i m_i.                                        \tag{AK07}
\]

Then \(x\mapsto m+x\) embeds that finite collared set into the strict interior
of \(S_N\).  It sends (AK04) to the literal FPSS parent and
shared-child incidences and preserves every label and displacement.  Thus
each finite local calculation can be placed inside one finite `GL6AA` MATCH
mission.  Equation (AK07) does **not** turn the infinite net into one infinite
record or certify a simultaneous terminal query at infinitely many sites.

## 3. Quasi-local algebra and inherited interaction

For finite \(F\subset\mathbb L\), let

\[
 \mathfrak A_F=\bigotimes_{p\in F}M_2(\mathbb C),\qquad
 \mathfrak A=\overline{\bigcup_{F\Subset\mathbb L}\mathfrak A_F}. \tag{AK08}
\]

On every site put the inherited one-link term, and on every edge put the
inherited degree-law pair term,

\[
 h_p=-hX_p+\varepsilon_\star n_p,
 \qquad
 \Phi_Z=2U_dn_pn_q\quad(Z=\{p,q\}\in E_\infty).              \tag{AK09}
\]

There is no new coupling or free weight in (AK09).  The all-formed route
means the inherited transverse gate is `beta_p=1` at every displayed site.
Route, controller, failure, and atlas-MISMATCH branches remain in the finite
complete instruments; (AK09) is the normalized active factor of the selected
all-formed MATCH member, not a success-filtered assertion about every branch.

The onsite terms define the exact product automorphism `alpha_t^on`.  In its
interaction picture,

\[
 \Phi_Z(t)=\alpha_{-t}^{\rm on}(\Phi_Z),\qquad
 \|\Phi_Z(t)\|\le J=2|U_d|,                                  \tag{AK10}
\]

with unchanged two-site support.  This treats the onsite dynamics exactly;
it is not discarded from the parent.

## 4. Explicit finite-volume boundary-error series

Use the authenticated cell distance

\[
 d_{A_3}(x,y)={1\over2}\|x-y\|_1.                            \tag{AK11}
\]

Let `B_R={x:d_A3(x,0)<=R}`, `Lambda_R=B_R x {1,2,3,4}`, and let
`tau_t^(R)` be the interaction-picture finite-volume dynamics containing
exactly the inherited pair terms whose two endpoints lie in `Lambda_R`,
combined with the global product onsite dynamics.  No artificial wraparound
or boundary coupling is added.

Let \(A\in\mathfrak A_X\), where \(X\subset\Lambda_{r_X}\), and write
\(|X|\) for its number of
active-link sites.  For a pair edge `Z`, define

\[
 r(Z)=\min_{p\in Z}d_{A_3}(\operatorname{parent}(p),0).        \tag{AK12}
\]

A cell belongs to six within-cell pair terms and twelve shared-child pair
terms.  Assign each pair term to an endpoint cell of smallest radius, with a
fixed tie rule.  Therefore the number of pair terms with `r(Z)=r` obeys the
uniform coarse bound

\[
 N_r\le18\,|B_r\setminus B_{r-1}|
     \le18(2r+1)^3.                                           \tag{AK13}
\]

Here \(B_{-1}:=\varnothing\), so (AK13) also covers the origin shell.

The last inequality follows because \(d_{A_3}(x,0)\le r\) implies every
coordinate \(x_i\in[-r,r]\), while the fourth coordinate is fixed by the
first three.

For completeness, the arbitrary-support extension used here is not inferred
by decomposing an arbitrary many-site operator.  Rerun the positive
Duhamel--Jacobi recursion of `GL6AI` with initial support indicator
\(\mathbf1_X\).  Linearity of that positive recursion and its one-site bound
give, for finite supports \(X,Y\),

\[
 \|[\tau_u(A_X),B_Y]\|
 \le2\|A_X\|\|B_Y\|
 \sum_{x\in X}\sum_{y\in Y}
 T_{d_L(x,y)}(\lambda_{\rm F3}|u|).                           \tag{AK13a}
\]

This is a support union bound on the recursion, not a norm-unsafe expansion
of \(A_X\) into one-site operators.  Apply (AK13a) to the two-site operator
\(B_Y=\Phi_Z(s)\), use \(\|\Phi_Z(s)\|\le J\), and use
\(d_L(x,y)\ge r(Z)-r_X\) for both endpoints.  For every finite intermediate
volume and every \(Z\) with \(r(Z)\ge r_X\), this yields

\[
 \|[\tau_u(A),\Phi_Z(s)]\|
 \le4J\|A\|\,|X|\,
 T_{r(Z)-r_X}(\lambda_{\rm F3}|u|),                           \tag{AK14}
\]

where

\[
 T_d(z)=\sum_{k=d}^{\infty}{z^k\over k!}.                    \tag{AK15}
\]

The factor four in (AK14) is the two-endpoint sum applied to the
dressing-complete `GL6AI` one-site bound.  If `S>R>r_X`, the exact Duhamel
identity, (AK13), and (AK14) yield

\[
 \begin{split}
 \|\tau_t^{(S)}(A)-\tau_t^{(R)}(A)\|
 &\le {72J\|A\||X|\over\hbar}
 \sum_{r=R}^{\infty}(2r+1)^3
 \int_0^{|t|}T_{r-r_X}(\lambda_{\rm F3}u)\,du.              \tag{AK16}
 \end{split}
\]

For \(U_d\ne0\), \(\lambda_{\rm F3}=24J/\hbar\), so the promised explicit
series is

\[
 \boxed{
 \mathcal E_R(A,t):=
 3\|A\||X|\sum_{r=R}^{\infty}(2r+1)^3
 T_{r-r_X+1}(\lambda_{\rm F3}|t|),
 \quad
 \|\tau_t^{(S)}(A)-\tau_t^{(R)}(A)\|\le\mathcal E_R(A,t).
 }                                                            \tag{AK17}
\]

For `U_d=0`, distinct sites factorize and the finite-volume dynamics of `A`
is exactly volume independent.  For every finite `t_0`, the factorial tail
in (AK17) beats the cubic shell factor uniformly on `|t|<=t_0`; hence

\[
 \lim_{R\to\infty}\sup_{|t|\le t_0}\mathcal E_R(A,t)=0.       \tag{AK18}
\]

It follows that

\[
 \boxed{\tau_t(A):=\lim_{R\to\infty}\tau_t^{(R)}(A)}         \tag{AK19}
\]

exists in norm, uniformly on compact time intervals, extends by density to
a strongly continuous automorphism group of `A`, and is independent of the
chosen locally complete open exhaustion.  More explicitly, for any two
exhaustions made only from restrictions of (AK09), replace `R` in (AK17) by
the largest interaction-radius about `X` that both exhaustions contain; that
radius tends to infinity and gives the same limit.  Arbitrarily strong or
new boundary laws are outside this statement.

## 5. Exact bulk covariance

For \(z\in A_3\) and \(\sigma\in S_4\), define

\[
 \theta_z(x,a)=(x+z,a),\qquad
 \rho_\sigma(x,a)=(\sigma x,\sigma a).                        \tag{AK20}
\]

Both maps preserve (AK04) and every coefficient in (AK09).  Finite-volume
covariance plus the unique limit (AK19) therefore gives

\[
 \boxed{
 \tau_t\theta_z=\theta_z\tau_t,qquad
 \tau_t\rho_\sigma=\rho_\sigma\tau_t.
 }                                                            \tag{AK21}
\]

This is exact covariance under authenticated displacement relabelings and
the four operation labels.  It is not Poincaré covariance and does not assign
a physical length or angle to the atlas.

## 6. Existence of a joint stationary bulk state

No finite-volume state limit is assumed.  Start with any state `omega_0` on
`A` and take the continuous-time averages

\[
 \omega_T(A)={1\over T}\int_0^T\omega_0(\tau_t(A))\,dt.       \tag{AK22}
\]

Weak-* compactness of the state space supplies a cluster point
`omega_stat`.  For fixed `s`,

\[
 |\omega_T(\tau_s(A))-\omega_T(A)|
 \le {2|s|\|A\|\over T},                                     \tag{AK23}
\]

so every such cluster point is stationary.

The group \(A_3\) is isomorphic to \(\mathbb Z^3\).  One explicit Følner
sequence is

\[
 F_R=\{(k_1,k_2,k_3,-k_1-k_2-k_3):|k_i|\le R\}.              \tag{AK24}
\]

Translation-average `omega_stat` over `F_R` and take another weak-* cluster
point.  Since (AK21) makes translation and time evolution commute, the new
state is both stationary and translation invariant.  Finally average it over
the finite group `S4`.  The result, denoted `bar omega`, obeys

\[
 \boxed{
 \bar\omega\tau_t=\bar\omega,qquad
 \bar\omega\theta_z=\bar\omega,qquad
 \bar\omega\rho_\sigma=\bar\omega.
 }                                                            \tag{AK25}
\]

Equation (AK25) proves existence, not uniqueness, purity, a ground-state or
KMS property, clustering, or operational preparation of `bar omega`.

## 7. Positive six-channel stationary spectral measure

At the reference cell define the six local pair observables

\[
 M_{ab}=Z_{(0,a)}Z_{(0,b)},\qquad 1\le a<b\le4,               \tag{AK26}
\]

and center them in the state (AK25):

\[
 \widehat M_A=M_A-\bar\omega(M_A)\mathbf1,
 \qquad A\in\binom{\{1,2,3,4\}}2.                            \tag{AK27}
\]

In the GNS representation of `bar omega`, stationarity implements `tau_t` by
a strongly continuous unitary group

\[
 U(t)=e^{itL},\qquad
 U(t)\pi(B)\Omega=\pi(\tau_t(B))\Omega.                       \tag{AK28}
\]

Let `P_L` be the projection-valued spectral measure of the self-adjoint
Liouvillian \(L\) and \(\psi_A=\pi(\widehat M_A)\Omega\).  Because
\(U(t)=e^{itL}\), the spectral coordinate \(\nu\) below has units of angular
frequency.  It is distinct from the energy coordinate used in `GL6AJ`; the
conversion is \(E=\hbar\nu\).  Then

\[
 \boxed{
 \mu_{AB}(\mathcal B)
 :=\langle\psi_A,P_L(\mathcal B)\psi_B\rangle
 }                                                            \tag{AK29}
\]

defines a finite positive matrix-valued Borel measure: for every
\(c\in\mathbb C^6\),

\[
 c^\dagger\mu(\mathcal B)c
 =\|P_L(\mathcal B)\sum_Ac_A\psi_A\|^2\ge0.                  \tag{AK30}
\]

The stationary correlation and retarded pair response are therefore

\[
 F_{AB}(t):=\bar\omega(\widehat M_A\tau_t(\widehat M_B))
 =\int_{\mathbb R}e^{it\nu}\mu_{AB}(d\nu),                  \tag{AK31}
\]

\[
 \boxed{
 \chi^R_{AB}(t)
 =-{i\over\hbar}\Theta(t)
 \left[F_{AB}(-t)-F_{BA}(t)\right].
 }                                                            \tag{AK32}
\]

`S4` invariance makes `mu` commute with the six-pair permutation
representation.  Since that representation is multiplicity-free,

\[
 \mathbb C^6\cong A_1\oplus E\oplus T_2,
 \qquad
 \boxed{
 \mu=\mu_{A_1}P_{A_1}+\mu_EP_E+\mu_{T_2}P_{T_2},
 }                                                            \tag{AK33}
\]

where each coefficient is a positive scalar measure.  This is an exact
stationary operator decomposition.  It does not say that any sector is
gapless, pole-dominated, propagating, or gravitational.

## 8. Displacement characters are not yet physical momentum

Let `varphi(k_1,k_2,k_3)=(k_1,k_2,k_3,-k_1-k_2-k_3)`.  Every
\(\vartheta\in\mathbb T^3\) defines the character

\[
 \chi_\vartheta(\varphi(k))=e^{i\vartheta\cdot k}.             \tag{AK34}
\]

For a finitely supported envelope `f:A3->C`, the bulk-aligned observable

\[
 M_A(f,\vartheta)
 =\sum_{x\in A_3}f(x)\chi_\vartheta(x)\theta_x(M_A)           \tag{AK35}
\]

is a legitimate local observable and therefore has the positive measure
(AK29).  The label `vartheta` is only a character of the authenticated
displacement group.  It is **not** physical momentum.  An unsmeared plane
wave, a thermodynamic spectral density at fixed character, and a
long-wavelength limit require additional summability, normalization,
refinement, and physical calibration theorems.

## 9. Exact result and next gate

The derived chain is

\[
 \boxed{\begin{gathered}
 \text{finite independently authenticated F3 neighborhoods}
 \longrightarrow \text{one locally finite }A_3\text{ interaction}\
 \longrightarrow \text{explicit boundary-Cauchy series}
 \longrightarrow \text{unique quasi-local bulk dynamics}\
 \longrightarrow \text{exact }A_3\rtimes S_4\text{ covariance}
 \longrightarrow \text{joint stationary states and positive}\
 \text{six-channel }A_1/E/T_2\text{ spectral measures}.
 \end{gathered}}                                               \tag{AK36}
\]

This closes the mathematical passage from the finite transient parent to a
stationary bulk response object without importing a graviton, pole, metric,
or Ricci ansatz.  The shortest remaining physics gate is to give a stationary
state an authenticated bulk-windowed source/read mission and determine,
under controlled enlargement/refinement, whether its character-resolved
`A1`, `E`, and `T2` measures show a stable collective infrared law.  Only the
resulting complete operator may later be decomposed as
`H2=kappa_R M_R+H_res` through the conditional `GL6L` bridge.

`PASS__FINITE_AUTHENTICATED_NEIGHBORHOODS_EMBED_EXACTLY_IN_F3_SLABS__INFINITE_A3_SITE_AND_SHARED_CHILD_INCIDENCE_DEGREE6__NO_NEW_COUPLING__EXPLICIT_BOUNDARY_ERROR_SERIES__NORM_CAUCHY_THERMODYNAMIC_DYNAMICS__OPEN_EXHAUSTION_INDEPENDENCE__EXACT_A3_TRANSLATION_AND_S4_COVARIANCE__JOINT_STATIONARY_STATE_EXISTS_WITHOUT_FINITE_STATE_LIMIT__POSITIVE_SIX_CHANNEL_LIOUVILLIAN_SPECTRAL_MEASURE__A1_E_T2_DECOMPOSITION__CHARACTER_NOT_PHYSICAL_MOMENTUM__NO_POLE_CONE_RICCI_GRAVITY_OR_G__AUTHOR_FROZEN_POSTFREEZE_AUDIT_REQUIRED`

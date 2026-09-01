# Finite authenticated bulk pulse/defect response theorem

**Lane:** `LANE_CROSS_RFT_GRA_GL6AM_AUTHENTICATED_BULK_RESPONSE_FUNCTIONAL_V001`  
**Short name:** `GL6AM V001`  
**Date:** 2026-08-31  
**Status:** author draft after hostile self-scope; independent hostile review
required before freeze  
**Claim class:** exact finite-window operational binding of the authenticated
atlas, formation-support, pair-source, and complete-read missions to the
sealed `GL6AK` quasi-local dynamics; explicit boundary-independent pulse-word
and local-defect limits; stationary finite-window retarded and positive
Liouvillian spectral functionals; inherited factorial relational influence
tail; exact covariance and restricted `A1/E/T2` sector statement

**Not claimed:** one infinite record/query; autonomous selection of the
homogeneous all-formed member or of a stationary state; convergence of finite
ground/Gibbs/apparatus states; promotion of the finite `GL6AF/GL6AG`
coefficients to stationary bulk coefficients; positive dissipative spectral
weight; detailed balance or fluctuation--dissipation; a gap, pole, mode,
physical momentum, length, or speed; strict microcausality or a Lorentz cone;
stress, Ward/Bianchi closure, Ricci/Einstein form, gravity, or `G`.

## 1. Narrow question and answer

`GL6AA` makes any finite collection of parent IDs, port labels, literal
shared children, and terminal link values jointly queryable without letting
the source/controller read the expected atlas.  `GL6AF` and `GL6AG` make a
finite formation word a retained physical `K`-support branch and use the
complete pair read.  Their inherited `GL6V/GL6W` apparatus supplies the exact
pair pulse

\[
 \exp\!\left({\mathrm i\over2}\sum_Aj_AM_A\right),
 \qquad M_A=Z_aZ_b,                                           \tag{AM01}
\]

with all source ancillas returned blank and all terminal outcomes retained.
`GL6AK` separately supplies the infinite all-formed `A_3` dynamics and joint
stationary invariant states.  The missing question is whether a **finite**
authenticated intervention/read mission has one open-boundary-independent
bulk value.

The answer is yes, in two precisely different senses.

1. A finite sequence of the pulse insertions (AM01), followed by a finite
   read, converges in operator norm with an explicit sum of the sealed
   `GL6AK` boundary tails.  In any chosen `GL6AK` stationary state its first
   derivative is therefore a lawful thermodynamic-limit retarded response,
   and every finite source/read window has a positive matrix-valued
   Liouvillian correlation measure.
2. Changing finitely many all-formed transverse supports to an authenticated
   `GL6AF/GL6AG` word is a bounded local perturbation of the generator.  Its
   cocycle dynamics and matched finite-read contrasts also have
   boundary-independent limits.  Such a defect contrast is generally
   nonequilibrium and is **not** assigned the stationary positivity or the
   homogeneous `A1/E/T2` decomposition of the unperturbed state.

This is a thermodynamic completion of every fixed finite mission, not a limit
in which the mission itself becomes infinite.

## 2. Authenticated finite source and read operators

Use the sealed bulk site set and dynamics

\[
 \mathbb L=A_3\times\{1,2,3,4\},\qquad
 \tau_t:\mathfrak A\longrightarrow\mathfrak A.                \tag{AM02}
\]

For a cell `x`, a pair `A={a,b}`, and
`alpha=(x,A)`, put

\[
 M_\alpha=M_{x,ab}:=Z_{(x,a)}Z_{(x,b)},\qquad
 Q_\alpha:=E_\star M_\alpha,qquad E_\star={\hbar\over T_\star}.
                                                                    \tag{AM03}
\]

`T_star` is the registered probe clock of `GL6W`; it is not a bulk coupling.
For a finite set `I` of such labels and real settings `j_alpha`, the exact
finite source unitary is

\[
 V_I(\mathbf j)=
 \exp\!\left({\mathrm i\over2}
              \sum_{\alpha\in I}j_\alpha M_\alpha\right).       \tag{AM04}
\]

All displayed pair operators commute, including across cells, so (AM04) is
also any scheduled product of the corresponding sealed phase sandwiches.
The response Hamiltonian is off or refocused during each sandwich and is
source-off between sandwiches, exactly as in `GL6V/GL6W`; no arbitrary
simultaneous smooth drive is inserted here.

For every fixed finite `I`, `GL6AK` equation (AK07) embeds a collar containing
all source and read supports in the strict interior of one finite FPSS slab.
`GL6AA` then appends actual ID, port, shared-child, clock, and complete
terminal records.  Its ID/label banks are held as the identity during the
response and are spectators of (AM02)--(AM04).  Hence they operationally
index the local algebra without changing its dynamics.  Every
`MATCH/MISMATCH/BLANK/COLLISION/FAILURE` output remains in that finite
instrument; the selected ideal `MATCH` active factor is not declared an
unconditional law.

## 3. Explicit pulse-word boundary limit

Let `tau_t^(R)` be the locally complete open restriction used in `GL6AK`, and
let

\[
 \mathcal E_R(A,t)=3\|A\|\,|X_A|
 \sum_{r=R}^{\infty}(2r+1)^3
 T_{r-r_A+1}(\lambda_{\rm F3}|t|),
 \qquad \lambda_{\rm F3}={48|U_d|\over\hbar},                  \tag{AM05}
\]

for local `A` supported on `X_A` inside cell radius `r_A`, with `R>r_A`.
For `U_d=0`, set `mathcal E_R=0`, as in the sealed theorem.

Choose a fixed finite ordered pulse list
`0<=s_1<=...<=s_n<=t` and finite unitaries `V_k` of the form (AM04).  Define

\[
 W_R=\tau_{s_n}^{(R)}(V_n)\cdots\tau_{s_1}^{(R)}(V_1),\qquad
 \mathcal O_R(B)=W_R^*\tau_t^{(R)}(B)W_R,                      \tag{AM06}
\]

where `B` is any fixed local read; the authenticated pair read has
`B=M_beta` or a finite real linear combination of such operators.  For
`S>R` large enough to contain every bare support, telescoping products and
unitary invariance give the exact comparison

\[
 \boxed{
 \|\mathcal O_S(B)-\mathcal O_R(B)\|
 \le \mathcal E_R(B,t)
   +2\|B\|\sum_{k=1}^n\mathcal E_R(V_k,s_k).}                  \tag{AM07}
\]

Indeed,
`||W_S-W_R|| <= sum_k ||tau_(s_k)^(S)(V_k)-tau_(s_k)^(R)(V_k)||`,
and each term is bounded by (AM05).  The right side tends to zero uniformly
when the finite list of times stays in a compact interval.  Thus

\[
 \boxed{
 \mathcal O(B):=W^*\tau_t(B)W
 =\lim_{R\to\infty}\mathcal O_R(B)}                            \tag{AM08}
\]

exists in norm and is independent of every locally complete open exhaustion
allowed by `GL6AK`.  The same product estimate gives a unique two-branch
pulse-insertion functional

\[
 \mathcal Z_\omega[\mathbf j_+,\mathbf j_-]
 :=\omega\!\left(W(\mathbf j_-)^*W(\mathbf j_+)\right),        \tag{AM09}
\]

and a unique sourced read
`omega(mathcal O(B))` for every state `omega` on `mathfrak A`.  No
finite-volume state limit is needed for this operator statement.  If a net
of finite-state extensions converges weak-* to `omega`, (AM07) also gives the
corresponding expectation limit; `GL6AK` does not assert that ground or Gibbs
states supply such a net.

## 4. Stationary retarded response and relational causal tail

Let `omega` be any one of the nonempty family of joint
time/translation/`S4`-invariant states whose existence is proved in `GL6AK`.
This is a conditional choice, not a selected vacuum or equilibrium phase.
Differentiating the exact source pulse with the sealed `GL6W` half-source
sign convention gives

\[
 \boxed{
 \mathcal G^R_{\beta\alpha}(t)
 ={\mathrm iE_\star^2\over2\hbar}\Theta(t)\,
 \omega\!\left([\tau_t(M_\beta),M_\alpha]\right).}            \tag{AM10}
\]

Equivalently, `mathcal G^R=-(E_star^2/2)chi^R_MM` with the `GL6W`
indexing.  Embed each finite-open observable in `mathfrak A` and define
`mathcal G^(R)` by the same formula with `tau_t` replaced by `tau_t^(R)` and
with `omega` restricted to that embedded algebra.  Applying (AM05) to
`M_beta` gives the explicit boundary comparison

\[
 \left|\mathcal G^{R,(S)}_{\beta\alpha}(t)
       -\mathcal G^{R,(R)}_{\beta\alpha}(t)\right|
 \le {E_\star^2\over\hbar}\,\mathcal E_R(M_\beta,t),          \tag{AM11}
\]

so (AM10) is a thermodynamic-limit response and not merely a formal GNS
commutator.

Let `X_alpha,X_beta` be the two-link supports and let `d_L` be the infinite
authenticated link-graph distance.  Passing the sealed `GL6AI` arbitrary-
support commutator bound through the `GL6AK` norm limit yields

\[
 \boxed{
 |\mathcal G^R_{\beta\alpha}(t)|
 \le {E_\star^2\over\hbar}\Theta(t)
 \sum_{p\in X_\beta}\sum_{q\in X_\alpha}
 T_{d_L(p,q)}(\lambda_{\rm F3}t).}                            \tag{AM12}
\]

The same formula with the operator norms and finite support sums retained
holds for arbitrary finite source/read combinations.  Equation (AM10)
vanishes exactly for negative response time because the authenticated read
is after the source.  Equation (AM12) is factorial quasi-local suppression,
not zero outside a spatial cone: at every `t>0` its tail may be nonzero.

## 5. Positive finite-window spectral functional

Center every finite-window pair operator,

\[
 \widehat M_\alpha=M_\alpha-\omega(M_\alpha)\mathbf1.
\]

In the GNS representation of `omega`, let
`U(t)=exp(i t L)` implement `tau_t`, let `P_L` be the spectral measure of
`L`, and put `psi_alpha=pi(widehat M_alpha)Omega`.  For every finite
source/read window `I`,

\[
 \boxed{
 \mu_{\alpha\beta}(\mathcal B)
 :=\langle\psi_\alpha,P_L(\mathcal B)\psi_\beta\rangle,
 \qquad \alpha,\beta\in I}                                   \tag{AM13}
\]

is a finite positive matrix-valued Borel measure:

\[
 \sum_{\alpha,\beta}\overline c_\alpha
 \mu_{\alpha\beta}(\mathcal B)c_\beta
 =\left\|P_L(\mathcal B)\sum_\alpha c_\alpha\psi_\alpha
  \right\|^2\ge0.                                             \tag{AM14}
\]

With

\[
 F_{\alpha\beta}(t)
 :=\omega(\widehat M_\alpha\tau_t(\widehat M_\beta))
 =\int_{\mathbb R}e^{\mathrm it\nu}\mu_{\alpha\beta}(d\nu), \tag{AM15}
\]

the authenticated retarded kernel is the exact spectral functional

\[
 \boxed{
 \mathcal G^R_{\beta\alpha}(t)
 ={\mathrm iE_\star^2\over2\hbar}\Theta(t)
 \left[
  \int e^{-\mathrm it\nu}\mu_{\beta\alpha}(d\nu)
 -\int e^{+\mathrm it\nu}\mu_{\alpha\beta}(d\nu)
 \right].}                                                     \tag{AM16}
\]

Thus the correlation/noise spectral measure is positive for every finite
window.  At a merely stationary state its support may lie on both signs of
frequency.  The commutator measure entering (AM16) is a signed difference of
positive measures and need not be positive at positive frequency.  No
detailed balance, passivity, KMS relation, fluctuation--dissipation theorem,
gap, or pole follows.

Translation and label invariance give only the exact covariance laws

\[
 \mu_{(x,A),(y,B)}=\mu_{(0,A),(y-x,B)},\qquad
 \mu_{(\sigma x,\sigma A),(\sigma y,\sigma B)}
 =\mu_{(x,A),(y,B)}.                                           \tag{AM17}
\]

For the six pairs in one cell, or for the six operators

\[
 M_A(f)=\sum_xf(x)M_{x,A}                                     \tag{AM18}
\]

with one finite real envelope satisfying `f(sigma x)=f(x)` for every
`sigma in S4`, the six-channel measure commutes with the pair permutation
representation.  The sealed multiplicity-free decomposition therefore gives

\[
 \boxed{
 \mu^{(f)}=\mu_{A_1}^{(f)}P_{A_1}
           +\mu_E^{(f)}P_E+\mu_{T_2}^{(f)}P_{T_2},}             \tag{AM19}
\]

with positive scalar measures.  An arbitrary envelope, a partial formation
word, or a localized defect need not close under `S4`; for it (AM17), not
(AM19), is the lawful statement.  Finite character-modulated envelopes are
included through their real and imaginary Hermitian parts, but their
character is not physical momentum and no fixed-character infinite-window
density is proved.

## 6. Authenticated finite formation words are local defects

The all-formed `GL6AK` onsite term is `-hX_p+epsilon_star n_p`.  Let
`D` be a finite subset of `mathbb L` and let an authenticated formation branch have
`kappa_p in {0,1}` on `D`, with the selected all-formed support retained
outside `D`.  Its active generator differs from the bulk generator by the
bounded local self-adjoint operator

\[
 \boxed{
 V_\boldsymbol\kappa
 =h\sum_{p\in D}(1-\kappa_p)X_p.}                              \tag{AM20}
\]

Adding (AM20) cancels precisely the transverse term when `kappa_p=0` and
changes neither `epsilon_star` nor any pair interaction.  This is the exact
`P^K` branch restriction used structurally by `GL6AF/GL6AG`, now placed in a
finite collared bulk mission; it is not a continuously fitted coupling.

Let the finite-volume interaction-picture cocycle solve

\[
 \mathrm i\hbar\,{dW_{R,\boldsymbol\kappa}(t)\over dt}
 =\tau_t^{(R)}(V_\boldsymbol\kappa)
  W_{R,\boldsymbol\kappa}(t),\qquad W_{R,\boldsymbol\kappa}(0)=\mathbf1,
                                                                    \tag{AM21}
\]

and define

\[
 \gamma_t^{(R),\boldsymbol\kappa}(A)
 =W_{R,\boldsymbol\kappa}(t)^*\tau_t^{(R)}(A)
  W_{R,\boldsymbol\kappa}(t).                                  \tag{AM22}
\]

Unitary Duhamel comparison and (AM05) give

\[
 \|W_{S,\boldsymbol\kappa}(t)-W_{R,\boldsymbol\kappa}(t)\|
 \le {1\over\hbar}\int_0^{|t|}
       \mathcal E_R(V_\boldsymbol\kappa,u)\,du,                \tag{AM23}
\]

\[
 \boxed{
 \|\gamma_t^{(S),\boldsymbol\kappa}(A)
       -\gamma_t^{(R),\boldsymbol\kappa}(A)\|
 \le \mathcal E_R(A,t)
 +{2\|A\|\over\hbar}\int_0^{|t|}
       \mathcal E_R(V_\boldsymbol\kappa,u)\,du.}              \tag{AM24}
\]

Consequently `gamma_t^kappa` is a unique boundary-independent strongly
continuous cocycle-perturbed automorphism group.  For two finite
authenticated words the matched read contrast

\[
 \Delta_{\boldsymbol\kappa|\boldsymbol\kappa'}^\omega(A,t)
 :=\omega\!\left(\gamma_t^{\boldsymbol\kappa}(A)
                 -\gamma_t^{\boldsymbol\kappa'}(A)\right)      \tag{AM25}
\]

is therefore a lawful thermodynamic-limit number for every state `omega`.

The finite-route `GL6AI` estimate is uniform in the complete `K` word, so a
single-site Duhamel comparison may be telescoped across every site on which
two finite words differ.  For a read `B_Y` supported on `Y`, its infinite-
volume matched contrast obeys the state-independent relational tail

\[
 \boxed{
 |\Delta_{\boldsymbol\kappa|\boldsymbol\kappa'}^\omega(B_Y,t)|
 \le {2h\|B_Y\|\over\hbar}
 \sum_{p\in D}|\kappa_p-\kappa'_p|
 \sum_{q\in Y}\int_0^{|t|}
 T_{d_L(p,q)}(\lambda_{\rm F3}u)\,du.}                         \tag{AM25a}
\]

For `U_d=0`, every remote cross-link term in (AM25a) is exactly zero, as in
`GL6AI`; an overlapping onsite read is not claimed to vanish.  For nonzero
`U_d`, (AM25a) is again a factorial quasi-local envelope, not strict spatial
support.

The `GL6AK` invariant state is generally not stationary under
`gamma^kappa`.  Each defect dynamics separately has stationary states by
weak-* time averaging, but neither authentication nor (AM20) chooses one or
identifies states between two defect words.  Hence (AM25) has no derived
positive stationary spectral measure, detailed balance, or homogeneous
sector split.  Only a symmetry that stabilizes the complete defect word may
be retained.

## 7. Exact scope of the bridge

The earned chain is

\[
 \boxed{\begin{gathered}
 \text{finite AA-authenticated cells/ports}
 \longrightarrow\text{finite physical pair pulse or retained }K\text{ word}
 \\
 \longrightarrow\text{sealed AK boundary-Cauchy dynamics}
 \longrightarrow\text{boundary-independent finite read and factorial
 relational tail};\\
 \text{stationary homogeneous pulse branch}
 \longrightarrow\text{retarded factorial tail + positive finite-window }
 \mu\\
 \longrightarrow A_1\oplus E\oplus T_2
 \text{ only for an }S_4\text{-closed window}.
 \end{gathered}}                                                \tag{AM26}
\]

What is **not** transported is equally exact.  `GL6AF`'s `U_d=0`, prepared-
blank entrance slope and `GL6AG`'s finite `N=1`, prepared-blank order-twelve
neighbor coefficient concern different states and finite generators.  They
authenticate the source/read and retained-branch ancestry used here, but
their numerical coefficients are not stationary bulk coefficients and are
not inserted into (AM10)--(AM25).

The next nonformal physics gate is therefore not another existence theorem.
It is an independently specified preparation or selection of a bulk state,
plus a controlled sequence of finite authenticated envelopes with the
summability/refinement and physical calibration needed to test an infrared
law.  Until then there is no physical momentum, pole, cone, Ricci response,
gravity, or `G`.

`PASS__FINITE_AA_AUTHENTICATED_SOURCE_READ_WINDOWS_ONLY__GL6V_W_PAIR_PULSES_BIND_TO_AK_DYNAMICS__EXPLICIT_PULSE_WORD_BOUNDARY_TAIL__THERMODYNAMIC_RETARDED_FUNCTIONAL_FOR_EVERY_CHOSEN_AK_INVARIANT_STATE__INHERITED_FACTORIAL_RELATIONAL_CAUSAL_TAIL_NOT_STRICT_CONE__FINITE_WINDOW_CORRELATION_MEASURE_POSITIVE__COMMUTATOR_DISSIPATIVE_SIGN_NOT_POSITIVE_WITHOUT_PASSIVITY_OR_KMS__A1_E_T2_ONLY_FOR_S4_CLOSED_WINDOWS__FINITE_AF_AG_K_WORD_IS_BOUNDED_LOCAL_DEFECT_WITH_EXPLICIT_COCYCLE_BOUNDARY_LIMIT__DEFECT_CONTRAST_NOT_STATIONARY_SPECTRAL__AF_AG_FINITE_COEFFICIENTS_NOT_PROMOTED__NO_INFINITE_QUERY_STATE_SELECTION_PHYSICAL_MOMENTUM_POLE_CONE_RICCI_GRAVITY_OR_G__INDEPENDENT_HOSTILE_REVIEW_REQUIRED`

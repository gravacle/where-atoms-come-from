# S4 categorical gamma--QFI ancestry theorem

**Lane ID:** `CROSS-RFT-GRA-EU-S4-CATEGORICAL-GAMMA-QFI-ANCESTRY-V001`

**Official short name:** `SCGQA`

**Date:** 2026-08-27

**Builder status:** `SOURCE_FROZEN_PENDING_INDEPENDENT_HOSTILE_AUDIT`

**Claim class:** exact finite-dimensional same-parent model theorem; exact
classical-query/state-fidelity type join; exact SLD-QFI and product-count
sufficiency theorem; exact conditional coordinate/scale match to EO's
tetrahedral contrast Gram

**Not claimed:** that the exponential family by itself forms a record; that
its parameter is physical position; that information geometry is spacetime;
that actual record descendants are independent products; that the actual
world instantiates this family; a physical soldering law, selected
Levi--Civita transport, continuum refinement, curvature response, stress
coupling, or gravity

## 1. Exact question and inherited boundary

GSGB separated the classical complete-query fidelity `gamma_Q` from the
squared Uhlmann state fidelity `gamma_state`, because they need not agree for
an arbitrary read.  It also left open whether the record contrast and the
rank-three QFI coframe candidate can belong to one state family rather than to
direct-product spectators.

This lane closes that type-and-parent gap in one exact model.  It supplies a
commuting four-outcome family for which:

1. the complete basis read is sufficient and saturates fidelity data
   processing exactly;
2. `gamma_Q=gamma_state` on the same physical family;
3. the SLD QFI is a positive rank-three tensor on
   `V=1^perp` and is exactly `S_4`-isotropic at the symmetric point;
4. independent product marks make gamma multiply and QFI add; and
5. the q4 count is a sufficient statistic whose fixed-depth contrasts are
   exactly of `A_3` type.

One independently qualified positive-margin record episode is then bound to
this family.  The record predicate is supplied by URFT formation, retention,
query, and lineage custody; it is not inferred from a softmax distribution.

The theorem is an intermediate ancestry lemma.  The real-world physical
soldering, connection, refinement, complete-response, and stress gates remain
explicit in section 9.

## 2. One exact `S_4`-covariant commuting family

Let

\[
 \mathbf1=(1,1,1,1)^{\mathsf T},\qquad
 P=I_4-\frac14\mathbf1\mathbf1^{\mathsf T},\qquad
 V=\mathbf1^\perp=P\mathbb R^4.                    \tag{EU01}
\]

Fix one dimensionless `lambda>0`.  For `theta in V`, define

\[
 p_a(\theta)=\frac{e^{\lambda\theta_a}}{Z(\theta)},
 \qquad
 Z(\theta)=\sum_{b=1}^4e^{\lambda\theta_b},          \tag{EU02}
\]

and on one four-level carrier

\[
 \rho(\theta)=\sum_{a=1}^4p_a(\theta)|a\rangle
 \langle a|.                                        \tag{EU03}
\]

Every state is full rank.  Restriction to `V` removes the common-shift gauge;
for `lambda>0`, the family is identifiable because

\[
 \log\frac{p_a}{p_b}=\lambda(\theta_a-\theta_b).     \tag{EU04}
\]

It covers the entire interior of the four-outcome simplex: for any strictly
positive `p`, its unique inverse in `V` is

\[
 \theta_a=\frac1\lambda\left(\log p_a-
 \frac14\sum_b\log p_b\right).                     \tag{EU04a}
\]

For `pi in S_4`, let `U_pi|a>=|pi(a)>` and let `Pi` be its real permutation
matrix.  Then

\[
 Z(\Pi\theta)=Z(\theta),\qquad
 \rho(\Pi\theta)=U_\pi\rho(\theta)U_\pi^\dagger.    \tag{EU05}
\]

Thus the whole family, not just four labels at one point, is exactly
`S_4`-covariant.  The unique `S_4`-fixed parameter in `V` is `theta=0`, where
`p_a=1/4`.

## 3. Theorem SCGQA-1 -- exact complete-query/state-gamma join

Use squared classical fidelity

\[
 \gamma_Q[p,q]
 =\left(\sum_{a=1}^4\sqrt{p_aq_a}\right)^2           \tag{EU06}
\]

and squared Uhlmann fidelity

\[
 \gamma_{\rm state}(\rho,\sigma)
 =\left(\operatorname{Tr}\sqrt{\sqrt\rho\,\sigma
 \sqrt\rho}\right)^2.                              \tag{EU07}
\]

The complete orthogonal basis query

\[
 M_a=|a\rangle\langle a|                            \tag{EU08}
\]

has law `p(theta)` and is sufficient for the family: its complete outcome law
determines the whole diagonal state.  For every `theta,theta' in V`,

\[
 \boxed{
 \gamma_Q[p(\theta),p(\theta')]
 =\gamma_{\rm state}[\rho(\theta),\rho(\theta')]
 =\frac{\left(\sum_a
 e^{\lambda(\theta_a+\theta'_a)/2}\right)^2}
 {Z(\theta)Z(\theta')}.}                            \tag{EU09}
\]

The value is one exactly when `theta=theta'`, and lies strictly between zero
and one otherwise.

### Proof

The matrices in (EU03) commute.  Therefore the positive square root inside
(EU07) is diagonal with entries `sqrt(p_a(theta)p_a(theta'))`; taking the
trace and squaring gives (EU06) and (EU09).  Equality to one in
Cauchy--Schwarz requires the normalized square-root probability vectors to be
equal, hence `p(theta)=p(theta')`, and identifiability on `V` then gives
`theta=theta'`.  Full support excludes zero fidelity. QED.

The data-processing direction is important.  For any quantum-to-classical
read `M`, fidelity monotonicity says

\[
 \gamma_{\rm state}(\rho,\sigma)
 \leq \gamma_Q[M(\rho),M(\sigma)].                  \tag{EU10}
\]

Reads can destroy distinguishability and thereby increase fidelity.  The
basis read (EU08) saturates (EU10); the equality is not assumed for an
arbitrary read.

## 4. Theorem SCGQA-2 -- exact SLD QFI and local gamma Hessian

In the ambient four-coordinate notation, the symmetric logarithmic
derivatives are

\[
 L_i(\theta)=\lambda\sum_{a=1}^4
 (\delta_{ai}-p_i)|a\rangle\langle a|,              \tag{EU11}
\]

and the SLD quantum Fisher tensor is

\[
 \boxed{
 {\cal F}_{ij}(\theta)
 =\lambda^2\bigl[p_i\delta_{ij}-p_ip_j\bigr]
 =\lambda^2[\operatorname{diag}p-pp^{\mathsf T}]_{ij}.}
                                                               \tag{EU12}
\]

It has null vector `1` and is positive definite after restriction to `V`.
At the symmetric point,

\[
 \boxed{{\cal F}(0)=\frac{\lambda^2}{4}P.}          \tag{EU13}
\]

For `dtheta in V`, the same squared fidelity in (EU09) obeys

\[
 -\log\gamma_{\rm state}(\theta,\theta+d\theta)
 =\frac14{\cal F}_{ij}(\theta)d\theta^id\theta^j
 +o(\|d\theta\|^2),                                \tag{EU14}
\]

or, with the Hessian taken in the second-endpoint displacement at
coincidence,

\[
 {\cal F}_{ij}=-2\,\partial_i\partial_j
 \log\gamma_{\rm state}\big|_{d\theta=0}.          \tag{EU15}
\]

### Proof

Differentiating (EU02) gives

\[
 \partial_i\log p_a=\lambda(\delta_{ai}-p_i).
                                                               \tag{EU16}
\]

Because `rho` is diagonal and full rank, (EU11) solves
`partial_i rho=(rho L_i+L_i rho)/2`.  Taking
`Tr[rho(L_iL_j+L_jL_i)/2]` gives (EU12).  For any real `x`,

\[
 x^{\mathsf T}{\cal F}x
 =\lambda^2\left[\sum_ap_ax_a^2-
 (\sum_ap_ax_a)^2\right],                           \tag{EU17}
\]

the positive weighted variance.  It vanishes only for a constant `x`; on
`V`, that means `x=0`.  Substitution of `p_a=1/4` gives (EU13).  Expanding
the exact expression (EU09) at coincidence gives (EU14)--(EU15). QED.

The basis query also saturates the quantum Cramer--Rao information bound on
this commuting family: its classical score is (EU16), so its classical
Fisher tensor is exactly (EU12).

## 5. Theorem SCGQA-3 -- product marks, sufficient q4 counts, and `A_3`

For `N>=1`, take an actual tensor-product family

\[
 \rho_N(\theta)=\rho(\theta)^{\otimes N}.            \tag{EU18}
\]

A complete word read has law

\[
 P_\theta(w)=\prod_{r=1}^Np_{w_r}(\theta),
 \qquad w\in\{1,2,3,4\}^N.                         \tag{EU19}
\]

Let `m_a(w)` count letter `a`, so `sum_a m_a=N`.  The count query has the
multinomial law

\[
 P_\theta(m)=\frac{N!}{\prod_am_a!}
 \prod_{a=1}^4p_a(\theta)^{m_a}.                    \tag{EU20}
\]

The likelihood in (EU19) depends on `w` only through `m`, so the q4 count is
a sufficient statistic.  More strongly, its coarse-graining loses no
fidelity:

\[
 \boxed{
 \begin{aligned}
 \gamma_Q[P_\theta(m),P_{\theta'}(m)]
 &=\gamma_Q[P_\theta(w),P_{\theta'}(w)]\\
 &=\gamma_{\rm state}[\rho_N(\theta),\rho_N(\theta')]\\
 &=\gamma_{\rm state}[\rho(\theta),\rho(\theta')]^N.
 \end{aligned}}                                     \tag{EU21}
\]

The product QFI and gamma information are

\[
 \boxed{{\cal F}_N=N{\cal F},\qquad
 I_{\gamma,N}:=-\log\gamma_N=N I_{\gamma,1}.}       \tag{EU22}
\]

On the fixed-depth count slice, define the observed count contrast

\[
 \xi_{\rm obs}=Pm=m-\frac N4\mathbf1.               \tag{EU23}
\]

For any two admitted counts `m,m'`,

\[
 m'-m\in A_3:=\{z\in\mathbb Z^4:\mathbf1^{\mathsf T}z=0\}.
                                                               \tag{EU24}
\]

The elementary count exchanges are the roots `e_b-e_a`, exactly the
fixed-depth contrast type used by EO.  Conversely, every `z in A_3` is the
difference of two nonnegative q4 counts at one common sufficiently large
depth, so the difference lattice across fixed-depth slices is exactly `A_3`.

### Proof

The factorization of (EU19) proves sufficiency.  Before squaring, the count
Bhattacharyya coefficient is

\[
 \begin{aligned}
 \sum_{|m|=N}\sqrt{P_\theta(m)P_{\theta'}(m)}
 &=\sum_{|m|=N}\frac{N!}{\prod_am_a!}
 \prod_a\bigl(\sqrt{p_ap'_a}\bigr)^{m_a}\\
 &=\left(\sum_a\sqrt{p_ap'_a}\right)^N,
 \end{aligned}                                      \tag{EU25}
\]

by the multinomial theorem.  Squaring gives (EU21), including equality with
product-state fidelity.  Product SLDs add over factors and cross terms have
zero mean, giving `F_N=N F`; logarithms give the second equality in (EU22).
Fixed total count gives (EU24). QED.

### Corollary SCGQA-3a -- bounded classical metric uniqueness

On the interior of the classical four-outcome probability simplex, impose the
standard regular classical statistical-category premise: a local Riemannian
metric must be invariant under congruent Markov embeddings and their
sufficient-statistic inverses.  Chentsov's theorem then makes the
Fisher--Rao metric unique up to one positive overall scalar.  Pulling that
metric back through (EU02) gives exactly (EU12), and the sufficient count map
of (EU20) preserves it.

This corollary is bounded in two ways.  It does not say the physical metric
must satisfy that statistical-category premise; `E-EMERGENTSPACE` or an
equivalent localization theorem must establish that.  It also makes no claim
that the SLD metric is the unique quantum monotone metric.  Here the quantum
and classical tensors agree only because the family is commuting and the
complete basis query saturates data processing.

Equation (EU22) applies to actual conditional tensor-product marks.  Mere
record redundancy does not imply that factorization; correlated descendants
require their joint state and joint fidelity.

## 6. Theorem SCGQA-4 -- binding one qualified record episode

Supply one authenticated physical record episode `R` with alternatives
`r in {K,B}` that independently passes the adopted URFT formation, retention,
lineage, and complete-query predicates with positive margin.  Require, as an
explicit same-parent model premise, that its alternative carrier states are

\[
 \rho_K=\rho(\theta_K),\qquad
 \rho_B=\rho(\theta_B),\qquad
 \theta_K\ne\theta_B,                               \tag{EU26}
\]

and that (EU08) is its authenticated complete query.

The binding includes a complete-port census: every alternative-dependent
query degree of freedom in the claimed episode is contained in (EU26).  Any
omitted factor must be identical in `K,B`, retained with the untouched
reference, and descend reference-stably; it cannot carry an unqueried
alternative-dependent spectator contrast.

Then the query gamma and state gamma of that same record are exactly one
object:

\[
 \boxed{
 \operatorname{REC}(R)\Longrightarrow
 \gamma_Q[p_K,p_B]
 =\gamma_{\rm state}(\rho_K,\rho_B)<1,
 \qquad I_\gamma(R)>0.}                             \tag{EU27}
\]

If `N` qualified conditionally independent descendants of that same supplied
family are physically retained, (EU21)--(EU22) give their exact count query,
gamma accumulation, and QFI addition.

### Proof

URFT qualifies `R`; it is not invoked to derive (EU02).  The explicit
same-parent binding (EU26) and complete read (EU08) invoke SCGQA-1 on the very
states that carry the supplied record alternatives.  Distinct parameters give
strictly subunit fidelity.  The conditional product conclusion is SCGQA-3.
QED.

This closes the direct type ambiguity for the supplied episode.  It does not
prove that every record has four outcomes, is commuting, is iid, or belongs to
this exponential family.

## 7. Exact coordinate and scale bind to EO

The parameter `theta` and the observed count contrast `xi_obs=Pm` are
different typed objects.  Statistics alone does not license
`theta=xi_obs`, and neither object is physical distance before soldering.

### 7.1 A declared half-contrast parameter chart

Let EO's unit append contrast be

\[
 v_a=Pe_a,\qquad
 v_a^{\mathsf T}v_b=\delta_{ab}-\frac14.             \tag{EU28}
\]

If a physical coordinate binding is prospectively supplied as

\[
 \theta(\xi)=\frac12\xi,\qquad
 d\theta_a=\frac12v_a,                              \tag{EU29}
\]

then, for one mark at the symmetric point and
`s_theta=ell_gamma^2 F(0)`, the four spatial edge Gram entries are

\[
 d\theta_a^{\mathsf T}s_\theta d\theta_b
 =\frac{\ell_\gamma^2\lambda^2}{16}
 \left(\delta_{ab}-\frac14\right).                  \tag{EU30}
\]

Therefore EO's target `a^2` on the diagonal and `-a^2/3` off diagonal is
obtained exactly by the one frozen scale lock

\[
 \boxed{\ell_\gamma^2\lambda^2=\frac{64a^2}{3}.}   \tag{EU31}
\]

Adding EO's independently calibrated common clock contribution `-a^2` then
reproduces its four-null-edge Lorentz Gram.  Equation (EU31) is a coordinate
and unit conversion under the supplied chart (EU29), not a prediction of
`a`, `lambda`, or physical spacetime.

For a general supplied chart `theta=kappa xi` and `N` product marks, the
pulled-back uniform QFI is

\[
 {\cal F}^{(\xi)}_N
 =\frac{N\kappa^2\lambda^2}{4}P,                    \tag{EU32}
\]

and the exact EO scale condition is

\[
 \boxed{\ell_\gamma^2N\kappa^2\lambda^2
 =\frac{16a^2}{3}.}                                 \tag{EU33}
\]

The choice `N=1,kappa=1/2` is precisely (EU31).  A direct identity chart
`theta=xi` instead gives `ell_gamma^2 lambda^2=16a^2/3`; confusing the two
charts creates a factor-of-four error.

### 7.2 The locally inferred mean-count chart

The model itself supplies only the expectation map

\[
 \bar\xi(\theta)=\mathbb E_\theta[Pm]
 =N\left(p(\theta)-\frac14\mathbf1\right).          \tag{EU34}
\]

At `theta=0`,

\[
 d\bar\xi=\frac{N\lambda}{4}d\theta,
 \qquad
 {\cal F}^{(\bar\xi)}_N=\frac4N P.                 \tag{EU35}
\]

If that mean-count coordinate, rather than (EU29), is physically soldered to
EO's unit count contrast, the corresponding scale lock is

\[
 \boxed{\ell_\gamma^2=\frac{Na^2}{3}.}             \tag{EU36}
\]

Equations (EU31), (EU33), and (EU36) are mutually consistent pullbacks in
different declared coordinates.  None chooses the physical chart.

### 7.3 Cross-lane length-convention ledger

This lane uses the sealed GSGB convention

\[
 s_{\rm GSGB}=\ell_F^2{\cal F}.                    \tag{EU36a}
\]

The separate ET packet writes the same candidate metric as

\[
 q_{\rm ET}=\frac{\ell_B^2}{4}Q.                   \tag{EU36b}
\]

When `Q=F` in the same coordinates and both symbols denote the same physical
metric, the conversion is

\[
 \boxed{\ell_B=2\ell_F.}                           \tag{EU36c}
\]

Every `ell_gamma` in (EU30)--(EU36) means `ell_F` in (EU36a).  It must not be
compared numerically with ET's `ell_B` without (EU36c).  For example, the
half-contrast lock (EU31) becomes
`ell_B^2 lambda^2=256a^2/3`, while the direct identity-chart lock becomes
`ell_B^2 lambda^2=64a^2/3`.  These are the same metrics under different
length-symbol conventions, not different physics.

## 8. What this lane closes

This lane establishes one exact same-parent ancestry segment:

\[
 \boxed{
 \begin{aligned}
 &\text{supplied qualified record }R
 \longrightarrow \rho(\theta_r)\\
 &\longrightarrow
 \gamma_Q=\gamma_{\rm state}<1\\
 &\longrightarrow
 {\cal F}=\lambda^2(\operatorname{diag}p-pp^{\mathsf T})\\
 &\xrightarrow{\theta=0}
 {\cal F}=\frac{\lambda^2}{4}P\\
 &\xrightarrow{\text{product marks}}
 \{\gamma_N=\gamma_1^N,\ {\cal F}_N=N{\cal F},
 \ m\text{ sufficient},\ \Delta m\in A_3\}.
 \end{aligned}}                                     \tag{EU37}
\]

Unlike GSGB's abstract separation, the classical record gamma, state gamma,
QFI, and q4 count now live in one explicitly constructed family.  No separate
spectator QFI factor is used.

## 9. What remains for the real-world gravity theorem

The following are still physics gates, not missing algebra inside SCGQA:

1. **Physical instantiation and lineage.**  Show that actual qualified record
   descendants realize the common categorical family or an experimentally
   validated generalization, with KEEP/BREAK custody and without a spectator
   direct product.
2. **Physical chart (`G-SOLDER`).**  Establish which, if any, map from record
   state parameters or count statistics to relational displacement is the
   one read by clocks, matter, electromagnetism, and both frame probes.
3. **Shared metric and transport.**  Show that neighboring cells agree on
   shared-face lengths and that actual probe transport is the selected,
   metric-compatible torsion-free connection of that same coframe.
4. **Controlled refinement.**  Supply a noncollapsed shape-regular scale
   family whose near-identity holonomies converge to the Levi--Civita
   curvature of the common metric.
5. **Same-sector response and complete stress.**  Bind that metric sector to
   the authenticated tensor response, complete conserved stress, protected
   spin-two, coefficient, remainder, variation-span, and causal-custody gates
   already isolated by the gravity program.

Only after these joins may the existing conditional Einstein--Hilbert and
back-reaction theorem be invoked as a real-world gravity result.  In
particular,

\[
 \boxed{\text{QFI information metric}\ne
 \text{physical spacetime metric without `G-SOLDER`.}}        \tag{EU38}
\]

## 10. Exact theorem status

`S4_CATEGORICAL_SAME_PARENT_GAMMA_Q_STATE_FIDELITY_JOIN_EXACT__SLD_QFI_RANK3_AND_S4_ISOTROPY_EXACT__IID_PRODUCT_GAMMA_MULTIPLICATION_QFI_ADDITION_AND_Q4_COUNT_SUFFICIENCY_EXACT__ONE_SUPPLIED_QUALIFIED_RECORD_EPISODE_BOUND_MODEL_CONDITIONALLY__EO_COORDINATE_SCALE_PULLBACK_EXACT__ACTUAL_WORLD_PHYSICAL_SOLDERING_LEVI_CIVITA_REFINEMENT_COMPLETE_STRESS_AND_GRAVITY_OPEN`

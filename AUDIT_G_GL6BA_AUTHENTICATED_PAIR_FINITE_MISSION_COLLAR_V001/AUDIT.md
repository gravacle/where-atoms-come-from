# Distinct hostile audit — GL6BA authenticated-pair finite-mission collar

**Target:** `LANE_CROSS_RFT_GRA_GL6BA_AUTHENTICATED_PAIR_FINITE_MISSION_COLLAR_V001/`  
**Frozen theorem SHA-256:** `d7ce0a7527a68f49e6ea2ee8edbb400a142fbb49297d8fe99cae78ffa0154ab0`  
**Frozen author-manifest SHA-256:** `6e14332230f713d51e393a5889fe78964fe0e63588b4b841533fa6af7ef19103`  
**Frozen author-seal-file SHA-256:** `34f29b3c03d53c4dbc9736d1bf7a7785e0a49a2aad04299b69a3804290c5971e`  
**Disposition:** `PASS__FULL_F3_FINITE_FPSS_COLLAR_EXACT__FINITE_BOUNDARY_COEFFICIENTS_CROSSING_CENSUS_PORT_DISTANCE_DUHAMEL_AND_TAYLOR_CONSTANTS_EXACT__AUTHENTICATED_BINARY_MARGINAL_ONLY__ALL_FINITE_R_AND_TIME__R2_AND_R5_OVER2_LICENSED_WITHOUT_ADHH__NO_GRAVITY_PROMOTION`

## 1. Custody and defect history

All twelve author files are byte-pinned in `AUDITED_TARGETS.sha256`.  The
independent hostile replay imports no author code.  It rebuilds the finite
FPSS incidences, `A_3` balls, link graph, boundary count, Taylor filtration,
state reduction, and binary total variation directly with standard-library
arithmetic.  Normal and optimized Python runs execute the same checks.

One scope defect was found and repaired while the author packet was still
mutable.  An earlier formulation allowed a merely induced mathematical open
set `Omega` to stand as the authenticated mission.  Embedding an induced cut
authenticates its internal relations but not arbitrary omitted boundary
couplings.  The frozen theorem instead takes its primary supremum over
`mathfrak F_L`, the complete finite all-formed/`MATCH` FPSS exteriors in which
the collar lies strictly inside the full authenticated active graph.  It
correctly demotes arbitrary induced opens to mathematical corollaries with at
most the displayed crossing terms.  This audit passes only those repaired,
sealed bytes.  No author byte was edited by this audit.

## 2. Exact finite-parent Hamiltonian, including FPSS boundaries

The selected exact parent is

\[
 H=U_d\sum_v(k_v-2)^2-h\sum_pX_p,
 \qquad R=U_d/h>0.
\]

The dangerous boundary question is whether a finite FPSS child star with
fewer than four incident links changes the displayed `-6R` onsite
coefficient.  It does not.  For every structural star `I_v`, of any admitted
size,

\[
 \left(\sum_{e\in I_v}n_e-2\right)^2
 =4-3\sum_{e\in I_v}n_e
  +2\sum_{\{e,f\}\subset I_v}n_en_f.                 \tag{A1}
\]

Every active link has exactly two original endpoints, including at the finite
FPSS boundary.  Therefore every link receives `-3R-3R=-6R`, and every
surviving pair has one physical vertex owner and coefficient `2R`.  The exact
finite dimensionless interaction is consequently

\[
 \overline H_{R,\Omega}
 =\sum_{p\in\Omega}(-X_p-6Rn_p)
  +2R\sum_{\{p,q\}\in E_\Omega}n_pn_q+C_\Omega.       \tag{A2}
\]

For `Lambda_L=B_L x {1,2,3,4}` inside a complete authenticated `Omega`, split

\[
 \overline H_{R,\Omega}=H_L+H_{\rm ext}+V_\partial+C_\Omega,
 \qquad
 V_\partial=2R\sum_{Z=\{p,q\}\in\partial\Lambda_L}n_pn_q. \tag{A3}
\]

The interior and exterior Hamiltonians have disjoint link supports.  Hence
their decoupled evolution of the central pair observable is exactly the
collar evolution.  No norm of the extensive exterior, fitted boundary term,
locked-space projection, effective Hamiltonian, or ADHH normal form enters.

## 3. Exact authenticated `A_3` collar and crossing census

The bulk link sites and interaction edges are

\[
 \mathbb L=A_3\times\{1,2,3,4\},\qquad
 A_3=\{x\in\mathbb Z^4:\mathbf1^Tx=0\},
\]

\[
 I(x;a,b)=\{(x,a),(x,b)\},\qquad
 S(x;a,b)=\{(x,a),(x+e_a-e_b,b)\}.
\]

Every link site has three same-parent and three same-child neighbors.  With

\[
 d_{A_3}(x,y)=\tfrac12\|x-y\|_1,
 \qquad B_L=\{x:d_{A_3}(x,0)\le L\},
\]

independent sign-composition gives

\[
 |B_r\setminus B_{r-1}|=10r^2+2,
 \qquad
 |B_L|={10L^3+15L^2+11L+3\over3}.                  \tag{A4}
\]

For one directed root `e_a-e_b`, an outward boundary cell has `x_a>=0` and
`x_b<=0`.  Its positive/negative-sum generating function is

\[
 { (1-zw)^2\over(1-z)^3(1-w)^3},
\]

whose diagonal coefficient is

\[
 {L+2\choose2}^2-2{L+1\choose2}^2+{L\choose2}^2
 =3L^2+3L+1.                                          \tag{A5}
\]

There are twelve directed roots.  Orienting every cut edge from its unique
inside endpoint counts each physical shared-child pair exactly once, so

\[
 \boxed{C_L=12(3L^2+3L+1).}                         \tag{A6}
\]

The translation by `(L+1,L+1,L+1,L+1)` places every collar parent at positive
FPSS coordinates.  Each outward neighbor is still nonnegative and both links
in (A6) terminate at the same literal child.  The complete finite exterior,
not the induced cut, carries the physical authentication.

## 4. Port-aware distance and interaction-cluster collar

A link-graph step is same-parent (`P`) or same-child (`C`).  A maximal
consecutive `C` run telescopes:

\[
 (x,a_0)\longrightarrow
 (x+e_{a_0}-e_{a_j},a_j).
\]

It changes cell distance by at most one.  If a path has `r` nonempty `C`
runs, it has at least `r` `C` steps and `r-1` separating `P` steps.  Therefore

\[
 \boxed{d_L(p,q)\ge
 2d_{A_3}(\operatorname{parent}p,\operatorname{parent}q)-1.} \tag{A7}
\]

A crossing edge ends in cell radius `L+1`.  Applying (A7) and deleting its
last link step gives, for either central support link `z`,

\[
 \boxed{d_L(z,p_{\rm inside})\ge2L.}                 \tag{A8}
\]

The replay finds attaining endpoints at every tested radius.  The same run
argument proves that every cluster with at most `m` interaction-picture pair
insertions lies in cell radius `ceil(m/2)`.  This is a support statement, not
a physical-time truncation.

## 5. Boundary-only Duhamel bound and every constant

The exact inherited arbitrary-support commutator estimate on the induced
collar is

\[
 \|[n_p,\tau_u^{(R,L)}(M_\beta)]\|
 \le2\sum_{z\in X_\beta}T_{d_L(z,p)}(48R|u|),
 \qquad |X_\beta|=2.                                  \tag{A9}
\]

For one crossing pair, Duhamel contributes `2R`; (A9) contributes `2`; and

\[
 \int_0^{|s|}T_d(48Ru)\,du
 ={1\over48R}T_{d+1}(48R|s|).                         \tag{A10}
\]

Thus the exact coefficient per crossing pair and initial support link is

\[
 (2R)(2){1\over48R}={1\over12}.                       \tag{A11}
\]

Using (A6), (A8), and the two initial links gives

\[
\begin{aligned}
 \sup_{\Omega\in\mathfrak F_L}
 \|\tau_s^{(R,\Omega)}(M_\beta)-\tau_s^{(R,L)}(M_\beta)\|
 &\le {1\over12}\sum_{Z\in\partial\Lambda_L}
       \sum_{z\in X_\beta}T_{d_L(z,p)+1}(48R|s|)\\
 &\le\boxed{2(3L^2+3L+1)T_{2L+1}(48R|s|)}.          \tag{A12}
\end{aligned}
\]

There is no missing endpoint factor, extra extensive-volume factor, or
off-by-one in the tail.  The trivial operator ceiling is two.  At `x=0`
because `s=0`, the two evolutions agree identically and `T_d(0)=0` for
`d>=1`.  At `R=0`, treated separately so that no division by `R` is used,
the links factorize and the collar error is again exactly zero.  No
logarithmic tail optimization is invoked.  At `L=0`, the exact census is
twelve and (A12) correctly starts with `T_1`.

For every fixed finite `R` and finite `s`, the factorial tail beats the
quadratic boundary count.  Hence a finite certified collar exists.  This does
not make one fixed collar uniform for unbounded time or prove that the
certified matrix is computationally small.

## 6. Ordinary Taylor filtration

This is separate from the interaction-picture finite-time estimate.  The
initial `M_beta` is diagonal.  Every pair term is diagonal, so the first
commutator that reaches a new link introduces that link diagonally.  Before a
second pair edge can continue along the path, a transverse onsite `X`
commutator must activate the new inside endpoint.  Diagonal onsite terms
cannot shorten this sequence.

The outside endpoint of a cut edge has link distance

\[
 D\ge2L+1.
\]

A nonzero exterior word therefore needs at least `D` pair commutators and
`D` transverse commutators:

\[
 k_{\rm first}\ge2D\ge4L+2.                           \tag{A13}
\]

Consequently the full and collared ordinary-time Taylor coefficients match
through order `4L+1`.  The replay instantiates the real two-root `M_beta` and
alternates transverse and pair commutators along a literal `C,P,...,C`
attaining geodesic.  This gives a nonzero formal ordered word at order
`4L+2`, so no smaller radius is guaranteed by the exact F3 term algebra.  The
formal-word statement is not promoted to a claim that the fully summed
coefficient is nonzero on every finite graph; such a claim would require a
separate cancellation analysis.  The independent replay also computes
aggregate full-versus-cut nested commutators on path reductions and finds
their first difference at exactly order `2D`.

This coefficient theorem does not promote a finite Taylor polynomial to an
uncontrolled finite-time approximation; (A12) is the finite-time control.

## 7. Exact state reduction, binary total variation, and no postselection

Fix one complete finite authenticated `Omega in mathfrak F_L` and its common
postformation state `omega_Omega`.  The collar state must be the exact
reduction of that same state.  For every collar operator `A`, including an
entangled exterior state,

\[
 \omega_\Omega(A\otimes I_{\rm ext})=\omega_L(A).       \tag{A14}
\]

The complete terminal link read has all sixteen outcomes
`q in {+1,-1}^4`, together with every retained route, clock, status,
authentication, and failure output.  The pair map `m_ab=q_aq_b` partitions
the sixteen outcomes into eight plus and eight minus outcomes and discards
nothing.  Its projectors are

\[
 P_\pm={I\pm M_\beta\over2}.
\]

For the two binary distributions produced from the common state and its exact
reduction,

\[
\begin{aligned}
 D_{\rm TV}(p^\Omega,p^{(L)})
 &=|p_+^\Omega-p_+^{(L)}|\\
 &=\tfrac12|\omega_\Omega(\tau^{(R,\Omega)}M_\beta
                      -\tau^{(R,L)}M_\beta)|\\
 &\le\boxed{\min\{1,(3L^2+3L+1)
 T_{2L+1}(48R|\sigma_{\rm obs}|)\}}.                 \tag{A15}
\end{aligned}
\]

The `1/2` is exact and the probability cap is one.  No read outcome or flag
is postselected: the pair marginal sums every retained flag value.  The
all-formed/`MATCH` active member is a prospectively selected deterministic
premise, not an observed-success conditioning step.  Equation (A15) controls
only that selected-factor binary pair marginal.  It does not bound total
variation of the joint flag/output distribution or an arbitrary random route
mixture.

## 8. Finite complete exteriors and the infinite ceiling

The primary physical statement is finite: `Omega` is one complete
authenticated FPSS mission, its full Hamiltonian retains every inherited
term, and the induced collar is only the comparison device.  A mathematical
induced open set has no more than the `C_L` crossing terms, so the same norm
bound holds, but that fact alone does not authenticate its omitted boundary.

Passing the uniform finite estimate through the independently sealed `GL6AK`
norm limit gives the same quasi-local bound.  This is a mathematical
completion.  It does not assert one infinite authenticated record, one
infinite terminal query, or an operationally prepared infinite state.

## 9. Moderate ratios and exact ceiling

No step above assumes large `R`, expands in `1/R`, or uses ADHH/`GL6AY`.
Within the inherited positive-coupling domain it is valid for every finite
`R` and finite mission duration.  The two admitted parameter assignments give

\[
 R=2:\quad48R|\sigma_{\rm obs}|=96|\sigma_{\rm obs}|,
\]

\[
 R=5/2:\quad48R|\sigma_{\rm obs}|=120|\sigma_{\rm obs}|.
\]

They are therefore directly licensed for a finite-collar calculation under
the full F3 Hamiltonian.  This does not validate the H6 approximation at
moderate coupling, select either ratio as nature's value, supply the mission
duration, select the postformation state, or guarantee practical numerical
cost.

The theorem proves no phase, record lifetime, full retained-output bound,
stationary state, pole, continuum, physical cone, graviton, Ricci tensor,
Einstein equation, gravity identification, or Newton constant `G`.  No
conventional gravity premise appears in the derivation.

Within these exact ceilings, every hostile attack survived.

**Hostile verdict: PASS.**

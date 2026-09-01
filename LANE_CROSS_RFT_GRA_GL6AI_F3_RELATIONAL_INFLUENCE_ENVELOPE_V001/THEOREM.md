# F3 authenticated relational-influence envelope theorem

**Lane:** `LANE_CROSS_RFT_GRA_GL6AI_F3_RELATIONAL_INFLUENCE_ENVELOPE_V001`  
**Short name:** `GL6AI V001`  
**Date:** 2026-08-31  
**Status:** author frozen after independent hostile pre-freeze review and
hostile-audited `GL6AH` dependency pin; post-freeze custody/replay audit
required  
**Claim class:** exact interaction-picture decomposition of the inherited
finite FPSS/F3 link Hamiltonian; exact uniform active-link degree census;
finite-graph commutator and matched source-branch Duhamel tails; exact
link-to-authenticated-cell distance comparison; uniform-in-`N` exponential
relational-influence envelope

**Not claimed:** exact finite-speed support or strict microcausality; a Lorentz
cone; a common physical cone
shared with matter or light; a stationary or long-wavelength propagating mode;
physical length or speed calibration; a continuum limit; Ward/Bianchi closure;
Ricci or Einstein form; gravity; or `G`.

## 1. Exact question and parent

The hostile-audited FPSS shared-child theorem `GL6Y` fixes, for every finite
`N`, the active append-edge set

\[
 E_N=\{e=(m,m+e_a):m\in S_N,\ a=1,\ldots,4\}                 \tag{AI01}
\]

and the inherited formed-KEEP link Hamiltonian

\[
 H_N=-h\sum_{e\in E_N}X_e+\Delta\sum_{e\in E_N}n_e
 +U_d\sum_{v\in S_N\sqcup S_{N+1}}(d_v-d_\star)^2,
 \qquad n_e={1-Z_e\over2}.                                  \tag{AI02}
\]

`GL6AA` independently authenticates the finite parent-cell atlas, including
literal shared children and the cell distance

\[
 d_{\mathcal G_N}(m,n)={1\over2}\|m-n\|_1.                  \tag{AI03}
\]

The immediately preceding `GL6AH` packet supplies a nonzero matched
record-formation signal over one direct shared-child connector in the full six
pair-coordinate read.  The present packet asks a different question: does the
same complete finite F3 generator impose a size-independent upper envelope on
how fast any local influence can spread across the authenticated atlas?

The answer below is an upper-bound theorem on the actual F3 parent.  It does
not infer a gravitational cone from the nonzero direct-edge coefficient.

## 2. Exact onsite-plus-pair decomposition

Fix one complete route sector and write
`beta_e=P_e^K in {0,1}` for the conserved KEEP gate on link `e`.  Formation,
writer, source, and query pulses are off during the interval.  Expanding each
degree square and using `n_e^2=n_e` gives, up to the retained
branch-independent guard scalar,

\[
 \boxed{
 H_{N,\beta}=H_{\rm on,\beta}+H_{\rm pair}+C_N,
 }
                                                                  \tag{AI04}
\]

\[
 H_{\rm on,\beta}
 =\sum_{e\in E_N}\left[-h\beta_eX_e+
 \varepsilon_\star n_e\right],
 \qquad
 \varepsilon_\star=\Delta+2U_d(1-2d_\star),                    \tag{AI05}
\]

\[
 \boxed{
 H_{\rm pair}=\sum_{\{e,f\}\in E(L_N)}2U_dn_en_f.
 }
                                                                  \tag{AI06}
\]

Here `L_N` is the line graph of the literal FPSS parent-child incidence graph:
two active link registers are adjacent exactly when their physical links share
a parent or a child.  The incidence graph is simple, so a pair is counted once.
Every term in (AI05) is one-link local.  Consequently
`exp(-it H_on,beta/hbar)` is a product of one-link unitaries and cannot enlarge
operator support.  In its interaction picture,

\[
 \Phi_{ef}^{(\beta)}(t)
 =2U_d\,n_e^{(\beta)}(t)n_f^{(\beta)}(t),
 \qquad
 \|\Phi_{ef}^{(\beta)}(t)\|\le J,\quad J:=2|U_d|,               \tag{AI07}
\]

with the same two-link support for every time and every route word.  Thus `h`,
`Delta`, and `d_star` rotate local operator axes but do not enter the
support-growth constant.  No term or retained port has been deleted.

## 3. Exact active-link degree census

Let `e=(m,c)` with `c=m+e_a`.  It has exactly three other links at parent `m`.
At child `c`, the incident active links are `(c-e_b,c)` for precisely the
indices with `c_b>0`.  If

\[
 q(c):=|\{b:c_b>0\}|,
\]

then the child contributes `q(c)-1` further partners, disjoint from the three
parent partners.  Hence

\[
 \boxed{
 \deg_{L_N}(e)=3+[q(c)-1]=q(c)+2\le6.
 }
                                                                  \tag{AI08}
\]

The upper bound is attained whenever all four child coordinates are positive.
It is independent of `N`; write `Delta_L=6` for the uniform ceiling.

## 4. A finite-graph commutator lemma with all dressing retained

For clarity, this section proves a conservative bound rather than silently
using a path-only constant that discards interaction dressing.  Define the
nonnegative influence matrix on active links by

\[
 \mathsf J_{ef}=\begin{cases}
 J,&e\ne f,\ e\sim f,\\
 J\deg_{L_N}(e),&e=f,\\
 0,&\text{otherwise}.
 \end{cases}                                                     \tag{AI09}
\]

Its diagonal is the exact bookkeeping term required by the Duhamel-Jacobi
recursion for interactions that touch the already reached support without
advancing its endpoint.  It is not a Hamiltonian self-coupling.  Equation
(AI08) gives

\[
 \|\mathsf J\|_\infty\le2J\Delta_L.                            \tag{AI10}
\]

Let `A_e` and `B_f` be operators on individual active links.  Applying the
interaction-picture Duhamel identity, the Jacobi identity, and
`||[R,S]||<=2||R||||S||` gives the componentwise Volterra inequality

\[
 C_{ef}(t)\le2\|B_f\|\delta_{ef}
 +{2\over\hbar}\sum_g\mathsf J_{eg}
 \int_0^{|t|}C_{gf}(s)\,ds,                                    \tag{AI11}
\]

where `C_ef` is the supremum of the commutator norm over unit-norm `A_e`.
Picard iteration of this positive inequality yields, for `e!=f`,

\[
 \boxed{
 \|[\tau_t^{(N,\beta)}(A_e),B_f]\|
 \le2\|A_e\|\|B_f\|
 \left[e^{2|t|\mathsf J/\hbar}\right]_{ef}.
 }
                                                                  \tag{AI12}
\]

Every nonzero off-diagonal factor of `mathsf J` crosses one edge of `L_N`.
Diagonal factors do not change the endpoint.  Therefore a matrix power
connecting `e` to `f` contains at least

\[
 d_L(e,f):=\operatorname{dist}_{L_N}(e,f)                        \tag{AI13}
\]

off-diagonal factors and in particular has total power at least `d_L(e,f)`.
Using (AI10), define

\[
 \boxed{
 \lambda_{\rm F3}:={4J\Delta_L\over\hbar}
 ={48|U_d|\over\hbar},
 \qquad
 T_d(x):=\sum_{r=d}^\infty{x^r\over r!}.
 }
                                                                  \tag{AI14}
\]

Then

\[
 \boxed{
 \|[\tau_t^{(N,\beta)}(A_e),B_f]\|
 \le2\|A_e\|\|B_f\|T_{d_L(e,f)}(\lambda_{\rm F3}|t|).
 }
                                                                  \tag{AI15}
\]

For an output supported on a finite link set `Y`, the same recursion and the
triangle inequality give the right side of (AI15) summed over `f in Y`.

The factor `48`, rather than the tempting path-only `24`, is deliberate:
`24|U_d|/hbar=2J Delta_L/hbar` counts advancing off-diagonal walks but does not
by itself bound arbitrary insertions that dress an already reached support.
No sharper constant is claimed without a separate exact no-dressing lemma.

## 5. Retained source-`K` Duhamel bound

Choose one literal source link `s`.  Compare two complete route sectors that
have the same normalized active-link Hilbert factor, common blank active-link
input, and identical `beta_e` for every `e!=s`, while `beta_s=1` in the first
sector and `beta_s=0` in the second.  This is exactly the retained source-`K`
comparator used by `GL6AH`.  It does not identify the `beta_s=0` ancestry with
both a sham/KEEP mission and a formed/BREAK mission; those are distinct full
instrument branches even when their normalized active-link generator agrees.
On the common active-link factor the two Hamiltonians differ by the physical
F3 term

\[
 V_s=-hX_s,qquad\|V_s\|=h.                                    \tag{AI16}
\]

All `K/G`, formation/sham, controller, failure, and terminal-query outcomes
remain in the complete instrument.  Define the signed, branch-normalized
active-output contrast

\[
 \Delta_s\langle B_Y(t)\rangle
 :=\langle B_Y(t)\rangle_{\beta_s=1,\beta_{-s}}
  -\langle B_Y(t)\rangle_{\beta_s=0,\beta_{-s}}.                \tag{AI16a}
\]

This is a fixed-sector dynamical contrast on the common active factor, not
success filtering and not a formal deletion switch.  A terminal observable
that directly reads the orthogonal `K/G` route record is intentionally not the
remote active-output observable bounded here.

For the common blank active-link input and any output `B_Y` supported on the
finite link set `Y`, the exact Duhamel identity plus (AI15) gives

\[
 \boxed{
 |\Delta_s\langle B_Y(t)\rangle|
 \le {2h\|B_Y\|\over\hbar}\sum_{f\in Y}\int_0^{|t|}
 T_{d_L(s,f)}(\lambda_{\rm F3}u)\,du.
 }
                                                                  \tag{AI17}
\]

For `U_d!=0`,

\[
 \int_0^{|t|}T_d(\lambda_{\rm F3}u)\,du
 ={1\over\lambda_{\rm F3}}
 T_{d+1}(\lambda_{\rm F3}|t|).                                \tag{AI18}
\]

If `U_d=0`, distinct links factorize and every cross-link contrast in (AI17)
is exactly zero.  GL6Y's separate all-BREAK comparator sets every `beta_e=0`;
its blank active block is an exact eigenstate of the remaining diagonal
Hamiltonian.  That separate statement is not substituted for the single-link
source-`K` comparator in (AI16a)--(AI17).

## 6. Authenticated cell-distance descent

Each active link register belongs to one parent cell `C_m`.  A link-graph edge
is of exactly one of two types:

1. a same-parent step, which stays in `C_m`; or
2. a same-child step from `C_m` to an adjacent authenticated cell `C_n`.

Project any link path from `e in C_m` to `f in C_n` to its parent-cell labels
and delete repeated labels.  The result is a walk in `G_N`; hence

\[
 \boxed{
 d_L(e,f)\ge d_{\mathcal G_N}(m,n)
 ={1\over2}\|m-n\|_1.
 }
                                                                  \tag{AI19}
\]

This statement is physical on the selected `GL6AA` `MATCH` atlas because the
shared child and both endpoint IDs are independently queried.  Every atlas
`MISMATCH`, `BLANK`, `COLLISION`, and failure result remains in the complete
terminal census and lies outside the asserted matched-cell identification.

Equations (AI15)--(AI19) therefore give one bound uniform in slab size `N` for
every fixed finite source and output cell block.  They do not assign a length
to a cell step.

## 7. Exponential relational-influence envelope

For every `mu>0`, the elementary marked-tail estimate

\[
 T_d(x)\le e^{xe^\mu-\mu d}                                    \tag{AI20}
\]

turns (AI15) into

\[
 \boxed{
 \|[\tau_t(A_e),B_f]\|
 \le2\|A_e\|\|B_f\|
 e^{-\mu[d_{\mathcal G_N}(m,n)-v_\mu|t|]},
 \qquad
 v_\mu={\lambda_{\rm F3}e^\mu\over\mu}.
 }
                                                                  \tag{AI21}
\]

If `s in C_m`, define the source-to-output cell separation

\[
 d_{\rm cell}(s,Y):=
 \min_{f\in Y}d_{\mathcal G_N}(m,\operatorname{parent}(f)).      \tag{AI21a}
\]

Then (AI17) implies

\[
 \boxed{
 |\Delta_s\langle B_Y(t)\rangle|
 \le {2h\|B_Y\||Y||\,|t|\over\hbar}
 e^{-\mu[d_{\rm cell}(s,Y)-v_\mu|t|]}.
 }
                                                                  \tag{AI22}
\]

At `mu=1`, the certified relational-envelope velocity is

\[
 \boxed{v_1={48e|U_d|\over\hbar}}
 \quad\text{cell steps per unit parent time, as an upper envelope only}. \tag{AI23}
\]

It is not a measured signal velocity.  `GL6AH` supplies a nonzero direct-edge
formed contrast, while (AI21)--(AI23) supply the uniform upper envelope.  The
combination earns an exponentially quasi-local **relational record-channel
influence envelope** on the selected authenticated atlas; it does not earn
exact finite-speed support or a stationary bulk mode.  At every nonzero time
the analytic tail in (AI15) is generally nonzero outside any nominal front.

## 8. Exact ceiling and next gate

The exact chain is

\[
 \boxed{\begin{gathered}
 \text{qualified F3 link-record route sectors}
 \to\text{ inherited }2U_dn_en_f\text{ connectors}\\
 \to\text{ uniform link degree }\le6
 \to\text{ exact finite-graph commutator/Duhamel tails}\\
 \to\text{ authenticated cell-distance bound}
 \to\text{ uniform exponential relational-influence envelope}.
 \end{gathered}}                                                   \tag{AI24}
\]

The result is a transient finite-parent quasi-locality bound in the record
channel.  It is
not a Lorentz cone because no common stationary dispersion, physical clock and
length calibration, or equality with independently sourced matter/light
responses has been proved.  The shortest next physics gate is therefore to
test the same homogeneous formed parent on a stationary/KMS or otherwise
controlled collective background for a stable long-wavelength propagating
mode.  Only after gluing, common-cone, Ward/Bianchi, and infrared operator
decomposition may any complete residual be compared with Ricci through the
conditional `GL6L` bridge.

`PASS__ACTUAL_F3_ONSITE_PLUS_PAIR_DECOMPOSITION_EXACT__ACTIVE_LINK_DEGREE_LE6_EXACT__CONSERVATIVE_DRESSING_COMPLETE_COMMUTATOR_MATRIX_BOUND_LAMBDA48UD_OVER_HBAR__NORMALIZED_SOURCE_K_BETA1_MINUS_BETA0_DUHAMEL_TAIL__AUTHENTICATED_CELL_DISTANCE_DESCENT__UNIFORM_N_EXPONENTIALLY_QUASILOCAL_RELATIONAL_INFLUENCE_ENVELOPE__DIRECT_NONZERO_EDGE_INPUT_FROM_HOSTILE_AUDITED_GL6AH__ANALYTIC_TAIL_NOT_EXACT_FINITE_SPEED__NOT_A_LORENTZ_OR_COMMON_PHYSICAL_CONE__NO_STATIONARY_MODE_RICCI_GRAVITY_OR_G__AUTHOR_FROZEN_POSTFREEZE_AUDIT_REQUIRED`

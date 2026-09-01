# Fixed-program formation-pattern `E`-coordinate source theorem

**Lane:** `LANE_CROSS_RFT_GRA_GL6AF_FORMATION_PATTERN_E2_SOURCE_V001`  
**Short name:** `GL6AF V001`  
**Date:** 2026-08-31  
**Status:** author draft; independent hostile audit required before freeze  
**Claim class:** exact branchwise response of the complete GL6T `F/S`
instrument for every retained formation pattern; exact two-record threshold
for the first nonzero response bilinear restricted to the fixed GL6AB `E`
plane; exact normalized GL6W entrance slope; fixed-program equal-count
orientation witness

**Not claimed:** that a pair observable is itself a record; an unconditional
postselected experiment; a collective interaction or stiffness; a multi-cell
bulk source; stress-energy; a continuum cone; Ricci response; gravity; or
`G`.

## 1. The source question left open by GL6AD

GL6AD proves that the earliest finite multi-cell `E2` opening has a cubic
piece proportional to an endpoint port-multiplicity imbalance.  On the
supplied FPSS boundary that imbalance may be a program artifact.  It cannot
be called a physical lineage source until changing authenticated formation
history, while leaving the program and response law fixed, changes a tensor
response.

That source test is already possible in the exact `N=0` F3 parent.  It does
not require a material scaffold or a gravitational ansatz.  The answer is:

\[
 \boxed{0\hbox{ or }1\text{ formed lineages}\Rightarrow E^TDE=0,
 \qquad 2\text{ formed lineages}\Rightarrow
 \operatorname{rank}(E^TDE)=1.}                   \tag{AF01}
\]

The result is source-side only.  It does not replace the multi-cell
propagation calculation.

## 2. One common physical parent, all branches retained

Use the exact four-link FPSS star and the complete tensor-product FPMH
formation instrument of GL6T.  On each physical link retain the formation
alternative, route, `K/G` support, controller, failure, and complete terminal
query.  The common response Hamiltonian on the `U_d=0` comparator is

\[
 H=\sum_{a=0}^3[-hP_a^KX_a+\Delta n_a],
 \qquad h,\Delta>0.                                \tag{AF02}
\]

No event label is inserted into (AF02).  On a fixed retained `F/S` and KEEP
branch define the measured support eigenvalues

\[
 \kappa_a\in\{0,1\},\qquad P_a^K=\kappa_a.         \tag{AF03}
\]

`kappa_a=1` means that link's independently authenticated formation was
routed to KEEP; `kappa_a=0` means its retained sham alternative.  GL6T's
terminal query proves the `F/S` contrast for each `a` with every other event
vector fixed.  Hence (AF03) is a physical branch restriction of one operator,
not a freely assigned response weight.  Every complementary branch remains
in the device census; no successful branch is renormalized into an
unconditional law.

More explicitly, the retained-pattern projector is

\[
 \Pi_{\boldsymbol\kappa}
 =\prod_{a=0}^3(P_a^K)^{\kappa_a}
                 (I-P_a^K)^{1-\kappa_a}.          \tag{AF03a}
\]

The support projectors commute with the response Hamiltonian and pair
queries,

\[
 [\Pi_{\boldsymbol\kappa},H]
 =[\Pi_{\boldsymbol\kappa},M_A]=0.                \tag{AF03b}
\]

They also commute with the GL6V source/read dilation because that dilation
acts on the active-link and probe registers and is the identity on the
`K/G/history` registers.  Thus

\[
 [\Pi_{\boldsymbol\kappa},U_{\rm src}(j)]=0,       \tag{AF03c}
\]

the probability of a retained `kappa` branch is independent of the source
amplitude, and GL6W source differentiation may be taken inside the normalized
fixed branch.  This is a branchwise response theorem, not an unconditional
success-filtered force law.

All active link qubits start blank.  After the common source-off prewait
`tau`, a formed link has the exact GL6T Bloch coordinates

\[
 x={2h\Delta\over\epsilon^2}(1-\cos\vartheta),
 \qquad
 z={\Delta^2+4h^2\cos\vartheta\over\epsilon^2},
 \qquad
 \epsilon=\sqrt{\Delta^2+4h^2},\quad
 \vartheta={\epsilon\tau\over\hbar},              \tag{AF04}
\]

while a sham link remains blank.  Thus

\[
 \widetilde x_a=\kappa_ax,
 \qquad z_a=1-\kappa_a(1-z).                       \tag{AF05}
\]

The fixed program, couplings, clock interval, source/read normalization, and
raw Hilbert space are identical for all `kappa`.  Only the retained physical
formation pattern differs.

## 3. Exact response for every formation pattern

Order the six pair observables as

\[
 M_{ab}=Z_aZ_b,qquad
 (ab)=(01),(02),(03),(12),(13),(23).               \tag{AF06}
\]

At the response entrance define, exactly as in GL6T,

\[
 D^{(\boldsymbol\kappa)}_{BA}
 =\operatorname{Tr}\rho_{\boldsymbol\kappa}(\tau)
 [[H,M_B],M_A].                                    \tag{AF07}
\]

Direct Pauli commutation and product-state evaluation give the complete
branchwise matrix:

\[
 \boxed{
 D_{BA}^{(\boldsymbol\kappa)}=
 \begin{cases}
 -4h(\widetilde x_a+\widetilde x_b),
   &A=B=\{a,b\},\\[2pt]
 -4h\widetilde x_s z_u z_v,
   &A=\{s,u\},\ B=\{s,v\},\ u\ne v,\\[2pt]
 0,&A\cap B=\varnothing.
 \end{cases}}                                      \tag{AF08}
\]

The second line is symmetric under `A<->B`.  Equation (AF08) is exact for
all sixteen retained patterns, not a small-time expansion.  It reduces to
the audited GL6T matrix for `kappa=(1,1,1,1)` and to zero when every
`kappa_a=0`.

## 4. Exact accumulation threshold on the fixed `E` plane

Use the GL6AB nonorthogonal program-frame embedding

\[
 E=\begin{pmatrix}
 1&1\\-1&0\\0&-1\\0&-1\\-1&0\\1&1
 \end{pmatrix},
 \qquad
 w_{ab}:=E_{(ab),:}.                               \tag{AF09}
\]

For every retained pattern define the restricted response bilinear

\[
 B_E^{(\boldsymbol\kappa)}:=
 E^TD^{(\boldsymbol\kappa)}E.                     \tag{AF09a}
\]

If the two coordinates of this fixed plane are changed by an invertible
matrix `S`, then `B_E` changes by congruence, `B_E -> S^T B_E S`; its rank is
therefore basis independent.  A partial formation pattern generally breaks
`S4`, however, so `D^(kappa)` need not preserve the `E` plane and can mix it
with its complement.  The theorem determines the exact `E x E` restriction;
it does not call that restriction an invariant `E2` Schur block for a
broken-symmetry branch.

For no formed link or exactly one formed link, (AF08) obeys

\[
 \boxed{E^TD^{(\boldsymbol\kappa)}E=0.}            \tag{AF10}
\]

Now let exactly the two links `a,b` be formed.  All other links remain in
their retained sham branches.  Exact contraction of (AF08) gives

\[
 \boxed{
 E^TD^{(ab)}E
 =-16hx(1-z)\,w_{ab}^{T}w_{ab}.}                  \tag{AF11}
\]

Here `w_ab` is a row covector, so `w_ab^T w_ab` is a `2 x 2` rank-one
matrix.  Under

\[
 0<\vartheta<2\pi,                                 \tag{AF12}
\]

GL6T proves `x>0` and `-1<z<1`; consequently the coefficient in (AF11) is
strictly nonzero.  The first nonzero direction in the fixed `E` restriction
therefore opens at two formed lineages, not at zero or one.

For this declared finite response target, (AF10) is the exact `REQUIRE`
side—one qualified lineage is insufficient—and (AF11) supplies the exact
`ALLOW` witness—one qualified pair is sufficient.  This is not promoted to
a universal two-record rule for arbitrary observables.

The six equal-count patterns do not insert six weights.  They are six
retained branches of the same `S4`-covariant instrument.  Their three
distinct covector lines are

\[
 (1,1),\qquad(-1,0),\qquad(0,-1),                 \tag{AF13}
\]

and span the fixed two-coordinate `E` plane.  Thus changing *which* two
lineages form rotates the opened rank-one restricted tensor direction while
preserving the number of formed records,
the response parameters, and the branch energy spectrum up to the inherited
port permutation.

For completeness, three formed lineages and four formed lineages have
rank-two `E x E` restrictions on the open domain (AF12).  For the representative
`kappa=(1,1,1,0)`,

\[
 E^TD^{(1110)}E
 =-4hx(3-2z-z^2)
 \begin{pmatrix}2&1\\1&2\end{pmatrix},           \tag{AF14}
\]

while for all four formed lineages

\[
 E^TD^{(1111)}E
 =-16hx(1-z^2)
 \begin{pmatrix}2&1\\1&2\end{pmatrix}.          \tag{AF15}
\]

Equations (AF10)--(AF15) are a response-rank threshold, not a claim of a
thermodynamic phase transition.

## 5. Physical normalized entrance slope

GL6W closes the six physical source/read dilations and proves

\[
 \left.\partial_t{\cal G}^R(t,0)\right|_{0^+}
 =-{E_\star^2\over2\hbar^2}D(\tau).               \tag{AF16}
\]

Therefore the two-record branch has the exact normalized fixed-plane slope

\[
 \boxed{
 E^T\left.\partial_t{\cal G}^{R,(ab)}(t,0)
 \right|_{0^+}E
 ={8E_\star^2hx(1-z)\over\hbar^2}
 w_{ab}^{T}w_{ab}.}                               \tag{AF17}
\]

Its sign and normalization are inherited from the same physical CTP
convention; they are not fitted.  The direct system seagull remains zero as
in GL6V/GL6W.

## 6. What this establishes and what remains

The exact source-side chain is

\[
 \boxed{
 \text{retained formation alternatives}
 \to\text{physical }K\text{-support pattern}
 \to\text{source-off pair response}
 \to\text{two-lineage fixed-}E\text{ opening}
 \to\text{pattern-selected tensor direction}.}    \tag{AF18}
\]

This proves a qualified formation-to-response chain: authenticated formation
with no-bypass retention supplies the physical `K`-support pattern, and that
pattern changes the finite response.  Equal record count with different
qualified placements gives permuted restricted tensor directions in one
fixed theory.

The Hamiltonian reads the physical support projector `P^K`, not the semantic
fact that a later query certifies `REC`.  An otherwise supplied `K=1` support
sector would reproduce the same response.  GL6AF therefore does **not** prove
an authentication-sensitive force.  It proves the narrower same-parent
ancestry chain

\[
 \text{qualified formation}\to\text{retained physical }K
 \to\text{fixed-plane response},                 \tag{AF18a}
\]

with record authentication establishing that the displayed `K` support can
be physically produced and retained by the F3 lifecycle.

The ceiling is equally important.  At `U_d=0` the link dynamics factorize;
AF11 is pair-coordinate response, not collective stiffness.  GL6U separately
proves that restoring the inherited degree interaction produces a genuine
inter-link factorization defect.  GL6AF does not yet solder the fixed-cell
source to GL6AD's multi-cell port imbalance, prove propagation around a
homogeneous interior, derive a common cone or conservation law, or identify
an infrared gravitational operator.  Those are the next operator-neutral
gates.

`PASS__ONE_COMMON_F3_PARENT_ALL_FS_BRANCHES_RETAINED__EXACT_ALL_16_FORMATION_PATTERN_RESPONSE__ZERO_OR_ONE_RECORD_FIXED_E_RESTRICTION_NULL__TWO_RECORD_FIRST_RANK1_FIXED_E_OPENING__BROKEN_S4_NO_INVARIANT_BLOCK_CLAIM__BRANCH_PROJECTOR_SOURCE_INDEPENDENT__PHYSICAL_K_NOT_SEMANTIC_REC_DYNAMICS__EQUAL_COUNT_PATTERN_SELECTS_COVECTOR_DIRECTION__THREE_AND_FOUR_RECORD_FIXED_E_RANK2__GL6W_NORMALIZED_CTP_ENTRANCE_SLOPE__PAIR_OBSERVABLE_NOT_CALLED_RECORD__NO_COLLECTIVE_STIFFNESS_BULK_STRESS_RICCI_GRAVITY_OR_G_CLAIM`

# Strict-lock T2 source-before-Feshbach contact and response theorem

**Lane:** `LANE_CROSS_RFT_GRA_GL6BV_STRICT_LOCK_T2_SOURCE_FESHBACH_CONTACT_V001`  
**Short name:** `GL6BV V001`  
**Date:** 2026-09-02  
**Status:** author theorem with exact replay; independent hostile audit required  
**Claim class:** exact local defect-frame algebra; exact static order-`h^2`
source-before-Feshbach functional in the pure local-pair `T2` source chart;
exact first and second source derivatives; exact arbitrary-locked-state onsite
and edge contacts; exact translation/`S4` expectation parametrization; leading
gapped microscopic retarded response; bounded period-four witness; active-
support scope test against the selected isolated-hexagon `E2` control

**Not claimed:** a complete six-pair source functional; a dynamically selected
finite-coupling state or thermodynamic phase; a low-energy `T2` pole; equality
of a static contact and a retarded coefficient; physical momentum, distance,
or metric normalization; `SO(3)`; a Ward identity; a common cone; a
Ricci/Einstein law; gravity; a graviton; or `G`.

## 1. Result first

At every degree-four constraint node use pair order

\[
 {\cal P}=(01,02,03,12,13,23),\qquad M_{ab}=Z_aZ_b.      \tag{BV01}
\]

Let `O` exchange opposite pairs and put

\[
 P_T={I-O\over2}.                                         \tag{BV02}
\]

For a degree-one or degree-three configuration, its pair vector is one of
four vectors `tau_a`.  They are pure `T2` and obey the exact tight-frame
identities

\[
 P_T\tau_a=\tau_a,\qquad
 \tau_a^T\tau_b=8\delta_{ab}-2,\qquad
 \sum_a\tau_a=0,\qquad
 \boxed{\sum_a\tau_a\tau_a^T=8P_T}.                    \tag{BV03}
\]

Starting from any local degree-two configuration and flipping each of its
four incident links once produces the multiset
`{tau_0,tau_1,tau_2,tau_3}`.  This fact is statewise; no average or selected
phase is used.

Now insert a genuine microscopic pair source before projection,

\[
 H(j)=U_d\sum_vq_v^2-h\sum_eX_e+\sum_v j_v^TM_v,
 \qquad P_Tj_v=j_v.                                      \tag{BV04}
\]

The restriction `P_T j=j` is part of the theorem.  It does not replace the
operator by a projected surrogate: the source still couples to the original
six pair observables.  Since `P M_T P=0`, the order-`h^2` Feshbach energy is
exact in the static source within its denominator domain:

\[
 \boxed{
 F^{(2)}(j)
 =-h^2\sum_{\eta\in\Omega}\sum_{e=\{v,w\}}
 { |\eta\rangle\langle\eta|\over
  2U_d+j_v^Tt_{v,e}(\eta)+j_w^Tt_{w,e}(\eta)} .}          \tag{BV05}
\]

Here `Omega` is the locked basis and `t_(v,e)(eta)` is the actual defect pair
vector at endpoint `v` after flipping `e`.  Every physical edge owns one and
only one denominator.

At source off, the complete local first derivative cancels whenever all four
incident one-link terms are active.  The Hessian is

\[
 \boxed{
 F^{(2)}_{j_vj_v}(0)=-{2h^2\over U_d^3}P_T,}              \tag{BV06}
\]

independent of the locked configuration, while for adjacent nodes
`e={v,w}`,

\[
 \boxed{
 F^{(2)}_{j_vj_w}(0;\eta)
 =-{h^2\over4U_d^3}
 t_{v,e}(\eta)t_{w,e}(\eta)^T .}                         \tag{BV07}
\]

The edge block is rank one and state-dependent; its reverse block is its
transpose.  Non-neighbor blocks vanish at this order.  Thus the anticipated
onsite coefficient `-2h^2/U_d^3 P_T` is confirmed exactly for `H+j.M` in the
pure-`T2` chart, but it is not the whole spatial kernel.

The same eliminated defects carry a leading microscopic retarded response at
gap `2U_d`.  In the inherited `i/2` commutator convention its onsite part is

\[
 K^R_{vv,T}(t)
 ={2h^2\over U_d^2}\,\theta(t)\sin(2U_dt)P_T+cdots .     \tag{BV08}
\]

It is gapped.  Equation (BV06) is its adiabatically eliminated static
curvature, not evidence for a low-energy `T2` mode.

## 2. Local defect algebra

Choose the representative with one exceptional negative port and define

\[
 (\tau_a)_{bc}=(-1)^{{\bf1}_{a\in\{b,c\}}}.              \tag{BV09}
\]

Explicitly,

\[
\begin{aligned}
 \tau_0&=(-1,-1,-1,+1,+1,+1),\\
 \tau_1&=(-1,+1,+1,-1,-1,+1),\\
 \tau_2&=(+1,-1,+1,-1,+1,-1),\\
 \tau_3&=(+1,+1,-1,+1,-1,-1).
\end{aligned}                                             \tag{BV10}
\]

A global sign reversal of all four `Z` values does not change a pair vector,
so (BV10) represents both degree one and degree three.  Their four-spin
product is `-1`; consequently opposite pair values have opposite sign and
`O tau_a=-tau_a`.  Their component sum is zero.  This proves that they are
pure `T2`.  Direct inner products give the Gram matrix in (BV03), and the
last frame identity follows because the four-vector Gram has nonzero
eigenvalue eight on its three-dimensional centered subspace.

For a locked spin vector `z`, let `sigma_z(a)` be the unique other port with
the same sign as port `a`.  Flipping port `a` leaves one exceptional sign at
`sigma_z(a)`, up to the irrelevant global reversal.  Therefore

\[
 t_a(z)=M(z\text{ with }z_a\mapsto-z_a)
       =\tau_{\sigma_z(a)}.                               \tag{BV11}
\]

The map `sigma_z` is the union of the two same-sign pairs and hence a
permutation of the four labels.  Equations (BV03) follow at every locked
node in every configuration.

## 3. Exact source convention and Feshbach derivation

Work first on the simple period-four regulator of `GL6AN`, or on a finite
collar for which the four incident links and both endpoints of every retained
edge are fixed.  Let `P` project onto `q_v=0` at every constraint node and
`Q=1-P`.  The pair source is diagonal and preserves these fixed blocks.  With

\[
 A=PHP,\quad B=PHQ,\quad C=QHP,\quad D=QHQ,\quad
 R=(z-D)^{-1},                                            \tag{BV12}
\]

the Feshbach operator is `A+BRC`.  This is precisely the derivative order of
`GL6BR`: source first, then resolvent differentiation, then any CTP or
spatial reduction.

For the source in (BV04), every locked node has zero `T2` pair value, so

\[
 P\sum_vj_v^TM_vP=0.                                     \tag{BV13}
\]

One link flip creates one unit defect at each endpoint, of source-off cost
`2U_d`.  All other nodes remain locked and have zero source energy.  The
intermediate-state cost is therefore

\[
 \Delta_e(j;\eta)=2U_d+j_v^Tt_{v,e}(\eta)
                         +j_w^Tt_{w,e}(\eta).             \tag{BV14}
\]

On the simple inherited incidence, two flips return to the lock only by
flipping the same edge twice.  Evaluating (BV12) at `z=0` gives (BV05).
Because (BV13) removes the order-one source shift of every locked state, the
source-dependent spectral value begins at order `h^2`; putting it back into
the resolvent changes the answer only at order `h^4`.  Thus each diagonal
entry of (BV05) is also the physical energy of its analytic locked branch
through order `h^2`, not merely an arbitrary fixed-`z` derivative.

The source domain is the connected neighborhood in which every displayed
denominator stays nonzero.  A sufficient local bound is
`max_v ||j_v|| < U_d/sqrt(6)`.  No claim is made through a defect-level
crossing.

Differentiating before setting the source to zero gives

\[
 {\partial F^{(2)}\over\partial j_v}
 =h^2\sum_{\eta,e\ni v}{t_{v,e}(\eta)\over\Delta_e^2}
       |\eta\rangle\langle\eta|,                         \tag{BV15}
\]

and

\[
 {\partial^2F^{(2)}\over\partial j_v\partial j_w}
 =-2h^2\sum_{\eta,e\supset\{v,w\}}
 {t_{v,e}(\eta)t_{w,e}(\eta)^T\over\Delta_e^3}
 |\eta\rangle\langle\eta|.                             \tag{BV16}
\]

For `v=w`, the sum is over its four incident edges.  Equations (BV03),
(BV15), and (BV16) prove

\[
 F^{(2)}_{j_v}(0)=0,qquad
 F^{(2)}_{j_vj_v}(0)=-{h^2\over4U_d^3}\,8P_T,            \tag{BV17}
\]

which is (BV06).  For distinct endpoints, one edge gives (BV07).

It is useful to collect the result without suppressing ownership.  Define

\[
 ({\cal A}_\eta j)_e
 =t_{v,e}(\eta)^Tj_v+t_{w,e}(\eta)^Tj_w .                 \tag{BV18}
\]

Then the complete real-space Hessian is

\[
 \boxed{F^{(2)\prime\prime}(0;\eta)
 =-{h^2\over4U_d^3}{\cal A}_\eta^*{\cal A}_\eta.}        \tag{BV19}
\]

It is negative semidefinite, has the onsite and edge blocks above, and has
no other blocks.  Equation (BV19) is also the bounded full-supercell Bloch
factorization on any periodic locked witness after inserting the appropriate
edge phases.

At `j=0` the locked manifold is degenerate.  For a generic multisite source,
the physical ground energy is the minimum of the branch functions in
(BV05); because (BV07) depends on the branch, that minimum can fail to have a
single bilinear Hessian at the origin.  Equations (BV19) and (BV29) are
therefore locked-branch or fixed-density Hessians, not a claimed unique
ground-energy Hessian.  The one-site onsite block (BV06) is universal across
four-incident-edge active branches.

The conclusion changes if the source is not restricted to `T2`.  A general
six-pair source acts nontrivially in `P`; the physical denominator then uses
the defect-minus-locked pair vector, and the spectral derivative must carry
the `P`-space source energy.  Extending (BV05) unchanged to `A1` or `E` would
be an error.

## 4. Contact versus microscopic retarded response

In the `GL6BR` CTP convention, eliminating `Q` turns (BV16) into the local-in-
time effective contact

\[
 W^{\rm ct}_{vw}(t,s)
 =-\langle F^{(2)}_{j_vj_w}(0)\rangle\,\delta(t-s).       \tag{BV20}
\]

Thus the sign in the connected generating functional is opposite the energy
Hessian sign.  At this order the effective first-source vertex (BV15)
vanishes identically at source off, so there is no additional low-energy
two-first-vertex Kubo term from `F^(2)`.

The same fact can be checked before elimination.  To first order in `h/U_d`,
a locked state contains the virtual component

\[
 |\eta\rangle_{\rm dressed}
 =|\eta\rangle+{h\over2U_d}\sum_e|\eta^e\rangle+cdots . \tag{BV21}
\]

The `T2` pair read has transition vector
`(h/(2U_d))t_(v,e)` to the defect state `eta^e`.  Therefore, for a stationary
locked density at leading order and

\[
 K^R={i\over2}\theta(t)
       \langle[M_T(t),M_T(0)^T]\rangle,                  \tag{BV22}
\]

the leading gapped kernel is

\[
 K_T^R(t)
 ={h^2\over4U_d^2}\theta(t)\sin(2U_dt)
   \langle{\cal A}_\eta^*{\cal A}_\eta\rangle+\cdots .   \tag{BV23}
\]

Its onsite specialization is (BV08), and its nearest-neighbor block has the
same state-dependent outer product as (BV07).  At zero frequency the stated
commutator normalization gives the exact leading-order relation

\[
 F^{(2)\prime\prime}(0)=-2K_T^R(\omega=0).               \tag{BV24}
\]

One may use the full microscopic Kubo history (BV23), or the adiabatically
eliminated contact (BV20).  Adding both as independent effects would double
count the same virtual defects.  Neither representation contains a
low-energy `T2` pole at this order.

## 5. Translation- and S4-invariant stationary expectation

No stationary phase is selected by the inputs.  Nevertheless, every
translation- and `S4`-invariant locked density has an exact one-parameter
form for the equal-time edge expectation.

For an edge of port `a`, neither endpoint defect label equals `a`.  Let

\[
 p=\Pr_\rho[\sigma_{z_v}(a)=\sigma_{z_w}(a)].             \tag{BV25}
\]

Translation and `S4` invariance make `p` independent of the edge and port.
The port stabilizer has only the equal-label and unequal-label orbits.  Put

\[
 Q_a={\tau_a\tau_a^T\over6},\qquad
 \beta=4p-{4\over3}.                                     \tag{BV26}
\]

Then

\[
 C_a:=\langle t_{v,e}t_{w,e}^T\rangle_\rho
 =\boxed{\beta P_T+\left({2\over3}-\beta\right)Q_a}.      \tag{BV27}
\]

This formula applies to the diagonal of a quantum density as well as to a
classical distribution, because the contact operators are diagonal in the
locked basis.  It exposes the state dependence in one scalar rather than
silently averaging it away.

Let `z_a=exp(i theta_a)` be the four edge-character phases in the symmetric
`A3` gauge, `sum_a theta_a=0`, and define

\[
 C(\theta)=\sum_az_aC_a,qquad
 I_2(\theta)=\sum_a\theta_a^2.                            \tag{BV28}
\]

The expectation of the edge Gram in parent/child sublattice order is

\[
 \boxed{
 G_\rho(\theta)=
 \begin{pmatrix}8P_T&C(\theta)\\C(\theta)^*&8P_T\end{pmatrix},
 \qquad
 \langle F^{(2)\prime\prime}\rangle
 =-{h^2\over4U_d^3}G_\rho.}                              \tag{BV29}
\]

At zero character,

\[
 C(0)=\gamma P_T,qquad
 \gamma={8\over3}(4p-1),                                 \tag{BV30}
\]

so the common and relative parent/child blocks have Gram coefficients
`8+gamma` and `8-gamma`.  They are not fixed by symmetry alone.

The full small-character expansion is

\[
\begin{aligned}
 C(\theta)={}&\gamma P_T
 +i\left({2\over3}-\beta\right)\sum_a\theta_aQ_a\\
 &-{1\over2}\left[
   \beta I_2(\theta)P_T
  +\left({2\over3}-\beta\right)\sum_a\theta_a^2Q_a
  \right]+O(|\theta|^3).                                 \tag{BV31}
\end{aligned}
\]

The second line is the parity-even `k^2` coefficient.  It is separate from
the `k^0` contact (BV30), and both are separate from the frequency-dependent
retarded kernel (BV23).  The direction `theta` is a dimensionless
translation character, not calibrated physical momentum.

## 6. Smallest authenticated periodic hexagon-support witness

The smallest sealed periodic regulator known here that has no wrapped
four-cycle and contains an alternating hexagon is the `Q4` witness of
`GL6AN/GL6AR`.  It has 64 cells, 128 constraint nodes, and 256 links.  The
exact deterministic background contains 218 links whose two endpoint
same-sign partners agree.  Averaging this one witness over its 64 translations
and 24 port permutations gives

\[
 p_{Q4}={218\over256}={109\over128}.                       \tag{BV32}
\]

This orbit mixture is a lawful order-`h^2` regulator: the source-off
effective Hamiltonian is scalar at that order.  It is not a parent-selected
phase, and stationarity under the full order-six hexagon Hamiltonian is not
asserted.

For this bounded witness,

\[
 \beta={199\over96},\qquad
 {2\over3}-\beta=-{45\over32},\qquad
 \gamma={77\over12}.                                     \tag{BV33}
\]

Thus the exact `k^0` common/relative Gram coefficients are

\[
 G_+(0)={173\over12}P_T,qquad
 G_-(0)={19\over12}P_T.                                  \tag{BV34}
\]

Let `Pi_+` and `Pi_-` be the fixed common and relative parent/child
sublattice projectors.  For the parity-even diagonal projections onto these
fixed subspaces, through quadratic character order,

\[
\begin{aligned}
 G_{+,{\rm diag}}^{\rm even}(\theta)={}&{173\over12}P_T
 -{199\over192}I_2(\theta)P_T
 +{45\over64}\sum_a\theta_a^2Q_a+O(|\theta|^4),\\
 G_{-,{\rm diag}}^{\rm even}(\theta)={}&{19\over12}P_T
 +{199\over192}I_2(\theta)P_T
 -{45\over64}\sum_a\theta_a^2Q_a+O(|\theta|^4).
                                                               \tag{BV35}
\end{aligned}
\]

Multiplication by `-h^2/(4U_d^3)` gives these projected energy-Hessian
blocks.  They are not eigenbranch dispersions or a Schur-reduced `k^2`
coefficient: the imaginary term linear in `theta` in (BV31) mixes the fixed
common and relative subspaces and can feed a quadratic term into either
operation.  The complete symbol is (BV29)--(BV31).  The verifier also
constructs the unaveraged `Q4` edge factorization `A_eta^* A_eta`; no large
matrix or phase selection is inferred from the compact orbit result.

## 7. Isolated-hexagon active-support obstruction

The sealed strict-lock packet `GL6BW` uses the fixed-boundary alternating-
hexagon doublet

\[
 H_c=-J\sigma_x,qquad
 J={63\over8}{h^6\over U_d^5},\qquad
 \rho_+=|+\rangle\langle+|,                              \tag{BV36}
\]

and obtains, after its explicitly chosen four-orientation completion,

\[
 K_E^R(t)=32a_*^2\sin(2Jt)P_E.                            \tag{BV37}
\]

The two locked basis configurations can also be embedded in the full `Q4`
parent.  In that separate embedding all four incident one-link terms are
active at every constraint node, so both configurations obey (BV06).  This
is the algebraic check performed by the verifier.

It does **not** attach (BV06) to `rho_+`.  The imported isolated control has a
six-link dynamic core and an eighteen-link *projected support*.  A projected
support statement is not a declaration that every one of those eighteen
microscopic transverse links is active in source-before-Feshbach histories.
If only the six core links are active, each cycle node has two, not four,
incident virtual flips.  Its local expressions are instead

\[
 F^{(2)}_{j_v}(0)={h^2\over4U_d^2}\sum_{e\in c:\,e\ni v}t_{v,e},
 \qquad
 F^{(2)}_{j_vj_v}(0)=-{h^2\over4U_d^3}
 \sum_{e\in c:\,e\ni v}t_{v,e}t_{v,e}^T,                \tag{BV38}
\]

where each sum has only two terms.  At a cycle node the two defect vectors
are distinct, their sum is nonzero, and their frame has rank two.  The full
hexagon toggle swaps the unordered two-vector set, so this partial anisotropic
contact is the same in the two `K2` basis configurations; it is nevertheless
not `8P_T`, and the four-edge conclusion (BV17) does not follow.

Declaring all eighteen projected-support links active might restore the four-
incident-edge identity at the six cycle nodes, but it defines a new
microscopic collar.  Its locked subspace need not remain the `K2` doublet,
its order-four/order-six source-free effective Hamiltonian has not been
replayed, and `rho_+` has not been proved stationary for it.  Therefore the
current packets supply same-microscopic-parent and pair-source-chart
ingredients, but no one completed source-first `E/T` functional and no one
common stationary state.

The sign convention in `GL6BW` is `H-J.M`; identifying `J=-j` makes the
linear vertices agree, and the second derivatives are unchanged.

There is consequently no Ricci result.  The available projected `E` control
has spectral scale `2J=O(h^6/U_d^5)`, whereas the full-parent `T` result
(BV23) has spectral scale `2U_d` and becomes an
instantaneous contact after elimination.  The four-orientation sum in
`GL6BW` is also a chosen finite composition, not a homogeneous parent state.
Most importantly, `h_E/2=h_T` is a condition on like coefficients of one
completed, identically normalized, soldered causal kernel.  A low-energy
retarded `E` coefficient cannot be equated to a static or high-gap `T`
coefficient.

## 8. Exact next gate and ceiling

The next bounded calculation is now identified.  Declare the complete
eighteen-link active microscopic collar and its fixed exterior, enumerate its
actual locked subspace, then insert the full six-pair source before the
canonical effective-Hamiltonian construction through orders `h^2`, `h^4`,
and `h^6`.  Every active incident flip, direct history, and fold must be
included.  The resulting `E-E`, `E-T`, and `T-T` contacts and first-source
vertices must then enter one branchwise CTP response.  Only that calculation
can decide whether a stationary low-energy `E/T` comparison exists on the
collar at the same temporal/spatial coefficient, orientation normalization,
and solder.

No large `Q4` census is needed before that test.  Nothing in (BV01)--(BV37)
adds `SO(3)`, a graviton, Einstein equations, a fitted contact, an imported
material law, or a Ward identity.

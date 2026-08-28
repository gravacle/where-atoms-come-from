# F3/q4 common-child mediation and collective-cone theorem

**Lane ID:** `GRA-FD-F3-Q4-CCMAC-V001`

**Short name:** `CCMAC`

**Date:** 2026-08-27

**Claim class:** exact finite q4 incidence/Schur-complement theorem; exact
relational `A3` bulk symbol and volume theorem; exact conditional collective
acoustic-cone theorem; exact separation of raw-front and collective cones

**Status:**
`EXACT_Q4_COMMON_CHILD_INCIDENCE__EXACT_PROSPECTIVE_DETUNED_F3_FORM_SCHUR_SIBLING_KERNEL__EXACT_A3_BULK_SECOND_MOMENT_COVOLUME_AND_REFINING_AFFINE_ATLAS__EXACT_SUPPLIED_MASSLESS_ACOUSTIC_ACTION_TO_SMOOTH_IR_PRINCIPAL_CONE__CARRIER_DETUNING_AND_PAIR_FIELD_LIFTS_COLLECTIVE_PHASE_COMMON_PROBES_CONSTRAINTS_STRESS_RGRL_B_AND_GRAVITY_OPEN`

**Not claimed:** that F3 already physically instantiates all q4 count modes as
coexisting carrier sites; that the bounded q4 witness supplies a global
coherent walk; that the collective phase exists; that controller time is
already physical proper time; a protected helicity-two pole, common matter
cone, diffeomorphism constraint, Einstein dynamics, or a numerical `G`.

## 1. Why this calculation is the next F3 test

The raw q4 append cone has four null extreme rays and remains polyhedral.  That
is an exact obstruction if **raw append reachability itself** is required to
be the final Lorentz cone.  Gravity Formation Theory instead proposes a
collective phase on record-built relational support.  A crystalline solid is
the elementary analogy: finitely many microscopic bond directions can support
an infrared acoustic equation whose principal characteristic cone is smooth.

The relevant question is therefore narrower:

> Does q4 common-child incidence, when physically lifted into one completely
> owned carrier parent with a supplied detuning, generate an isotropic
> quadratic sibling kernel; and would a separately supplied massless
> collective action on that support have a smooth `3+1` principal cone
> without inventing a direction net or a background grid?

The answer is **yes at the exact algebraic and explicitly conditional
collective levels**.  The calculation does not derive the detuning, the
massless coordinate, or its stiffness from current F3.  It identifies those
objects as the remaining physical interface rather than hiding them in the
word "phase."

## 2. Finite q4 fronts and the common-child incidence operator

For `N >= 0`, let

\[
 S_N=\{m\in\mathbb N_0^4:\mathbf1^{\mathsf T}m=N\}.
 \tag{FD01}
\]

Let `B_N : ell^2(S_N) -> ell^2(S_(N+1))` be the unweighted append-incidence
operator

\[
 (B_N)_{c m}=1
 \quad\Longleftrightarrow\quad
 c=m+e_a\text{ for one }a\in\{1,2,3,4\}.
 \tag{FD02}
\]

Every parent has exactly four children.  Define the sibling graph `G_N` on
`S_N` by joining two distinct fronts precisely when they share a child, and
let `A_N` be its adjacency matrix.

### Theorem `CCMAC-1` -- exact global common-child kernel

\[
 \boxed{B_N^\dagger B_N=4I_{S_N}+A_N.}
 \tag{FD03}
\]

Moreover, for distinct `m,m' in S_N`,

\[
 (A_N)_{m'm}=1
 \quad\Longleftrightarrow\quad
 m'-m=e_b-e_a\quad(a\ne b)
 \tag{FD04}
\]

with both endpoints nonnegative.  In the all-positive interior, every vertex
has the twelve oriented `A3` root neighbors.

#### Proof

The diagonal entry `(B^dagger B)_(mm)` counts the four children `m+e_a`.
For `m != m'`, an off-diagonal entry counts their common children.  Equality
`m+e_a=m'+e_b` is equivalent to (FD04).  A distinct pair can share at most
one child, so every off-diagonal entry is either zero or one and equals the
sibling adjacency.  QED.

Equation (FD03) is global finite algebra, not yet a physical carrier graph.
The bounded q4 witness realizes the count-front Hilbert spaces and complete
history custody, while the earlier FOCST theorem supplies only a prospective
local endpoint-to-carrier SWAP.  `CCMAC` therefore introduces the following
explicit interface rather than silently promoting those modes:

**`DETUNED-Q4-CARRIER-LIFT`.**  One F3-form parent physically realizes the declared
`S_N` and `S_(N+1)` front labels as coexisting carrier modes, acts identically
on the retained/future-blind history factor, and uses the q4 append incidence
as the support of the already declared content-symmetric F3 transfer term of
amplitude `t`.  In addition, it supplies a positive source-off child/parent
detuning `Delta`; the preparation and maintenance of that offset and every
associated source, work, controller, support, boundary, failure, and
quarantine port are explicitly owned.  Current F3 supplies the scalar
transfer form, not this detuning or carrier realization.  The complete lift
is a construction target, not a consequence of recordhood.

## 3. A prospectively detuned F3-form transfer gives an exact Schur kernel

Under `DETUNED-Q4-CARRIER-LIFT`, use the F3 scalar transfer amplitude `t` on
the q4 append support and the separately supplied positive source-off
detuning `Delta` of child modes relative to parent modes.  The finite
Hermitian block Hamiltonian is

\[
 H_N=
 \begin{pmatrix}
  0&-tB_N^\dagger\\
  -tB_N&\Delta I
 \end{pmatrix}.
 \tag{FD05}
\]

The two record contents remain an inert identity multiplicity, exactly as in
the F3 carrier term; no content outcome is selected.

For spectral parameter `z != Delta`, eliminating the child amplitude gives
the exact parent Schur equation

\[
 \boxed{
 \left[-zI-{t^2\over\Delta-z}
 (4I+A_N)\right]\psi_N=0.}
 \tag{FD06}
\]

Put `K_N=B_N^\dagger B_N=4I+A_N`.  The exact negative/low spectral branch,
represented on the parent space by functional calculus, is

\[
 H_{N,\mathrm{low}}
 :=f_{\Delta,t}(K_N)
 ={\Delta I-\sqrt{\Delta^2I+4t^2K_N}\over2}.
 \tag{FD07a}
\]

On the frozen smallness domain

\[
 {4t^2\|K_N\|\over\Delta^2}\le {1\over2},
 \tag{FD07b}
\]

Taylor functional calculus gives the operator-norm statement

\[
 H_{N,\mathrm{low}}
 =-{t^2\over\Delta}(4I+A_N)
 +{t^4\over\Delta^3}(4I+A_N)^2+R_N,
 \qquad
 \|R_N\|\le
 {2|t|^6\over\Delta^5}\|4I+A_N\|^3.
 \tag{FD07}
\]

Thus common-child mediation generates the **global scalar sibling hopping
kernel** without inserting sibling edges individually.  It is still
conditional on `DETUNED-Q4-CARRIER-LIFT`; (FD05) does not claim that the
current F3 microscopic parent has already realized the carriers or detuning.

As a separately truncated **one-common-child comparator**, retain only one
sibling pair and its shared child.  Its incidence row is `B=(1,1)` and the
three-mode block has one dark endpoint state of energy zero and bright
eigenvalues

\[
 E_\pm={\Delta\pm\sqrt{\Delta^2+8t^2}\over2}.
 \tag{FD08}
\]

The low bright energy is `-2t^2/Delta+O(t^4/Delta^3)`, exactly matching the
per-common-child contribution
`-(t^2/Delta)[[1,1],[1,1]]`.  It is not the full two-parent principal block
of `K_N`, because the full q4 incidence retains each parent's other children.

## 4. Exact `A3` bulk symbol and the scalar first-order obstruction

Put

\[
 V=\mathbf1^\perp\subset\mathbb R^4,
 \qquad P=I-\tfrac14\mathbf1\mathbf1^{\mathsf T},
 \qquad n_a={2\over\sqrt3}Pe_a,
 \tag{FD09}
\]

and define the six unoriented roots

\[
 \alpha_{ab}=n_b-n_a,\qquad 1\le a<b\le4.
 \tag{FD10}
\]

They obey the exact tight-frame identity

\[
 \boxed{
 \sum_{a<b}\alpha_{ab}\alpha_{ab}^{\mathsf T}
 ={16\over3}I_V.}
 \tag{FD11}
\]

At lattice scale `a_*`, the translation-invariant interior symbol of the
sibling adjacency is

\[
 A(k)=2\sum_{a<b}\cos(a_* k\!\cdot\!\alpha_{ab}),
 \tag{FD12}
\]

while the parent-to-child incidence symbol obeys

\[
 |b(k)|^2
 =\left|\sum_a e^{ia_*k\cdot n_a}\right|^2
 =4+A(k).
 \tag{FD13}
\]

Therefore

\[
 A(k)=12-{16\over3}a_*^2|k|^2+O(a_*^4|k|^4),
 \tag{FD14}
\]

and the scalar one-carrier band obtained from (FD06) is even in `k`; its
gradient at `k=0` is exactly zero.  The current scalar transfer by itself is
not a first-order Weyl/Dirac cone.  Calling its quadratic band a Lorentz cone
would be false.

## 5. Theorem `CCMAC-2` -- a collective phase has a smooth principal cone

Now impose one explicit dynamical antecedent instead of changing the graph.

**`MASSLESS-A3-COLLECTIVE-ACTION`.**  On the lifted support, a controlled
block supplies one real collective coordinate `phi_m(tau)`, an exact shift
symmetry `phi -> phi+constant` that forbids its mass term, and the complete
quadratic nearest-root action

\[
 S^{(2)}={1\over2}\int d\tau\left[
 \chi\sum_m\dot\phi_m^2
 -\kappa\sum_{m}\sum_{a<b}
 (\phi_{m+e_b-e_a}-\phi_m)^2\right],
 \qquad \chi,\kappa>0.
 \tag{FD15}
\]

Each undirected root edge is counted once by the selected six positive roots.
The time `tau` is the F3 controller/response time.  It becomes physical proper
time only after the common-clock and common-probe gates pass.

This is a separately supplied action antecedent.  The common-child Schur
calculation proves the `A3` support algebra but does **not** derive `phi`, its
shift symmetry, `chi`, `kappa`, or the absence of other quadratic terms from
current F3.

The exact bulk dispersion is

\[
 \boxed{
 \omega^2(k)={2\kappa\over\chi}
 \sum_{a<b}[1-\cos(a_*k\cdot\alpha_{ab})].}
 \tag{FD16}
\]

Using (FD11),

\[
 \boxed{
 \omega^2(k)=c_*^2|k|^2+O(a_*^4|k|^4),
 \qquad
 c_*^2={16\kappa a_*^2\over3\chi}.}
 \tag{FD17}
\]

Hence the infrared principal polynomial is

\[
 \boxed{-\omega^2+c_*^2|k|^2,}
 \tag{FD18}
\]

with the smooth `3+1` Lorentz characteristic cone.  The tetrahedral/fcc
anisotropy first enters beyond the quadratic principal term.

#### Proof

Fourier transformation of the sibling difference operator gives (FD16).
Taylor expansion and (FD11) give (FD17).  Positivity of `chi,kappa` gives one
hyperbolic time direction and three positive spatial principal directions,
which proves (FD18).  QED.

This result **does not contradict** the fixed-finite-direction raw-front
obstruction.  The two theorems concern different objects:

\[
 \begin{array}{c|c}
 \text{raw append convex cone}&\text{collective response characteristic}\
 \hline
 \text{four microscopic extreme rays}&
 \text{zero set of an infrared quadratic principal symbol}\\
 \text{remains polyhedral}&
 \text{is smooth by the isotropic root second moment}
 \end{array}
 \tag{FD19}
\]

Thus an `AFR` direction net is sufficient for a raw-front construction but is
not necessary for this collective-wave route.  What must be derived here is
the stable phase and its positive `chi,kappa`, not infinitely many
microscopic bond directions.

For a dynamical continuum composed with the geometric refinement
`a_* -> 0`, a finite nonzero physical speed is an additional scale binding:

\[
 {\kappa(a_*)\over\chi(a_*)}
 ={3c_{\rm phys}^2\over16a_*^2}.
 \tag{FD18a}
\]

Under this binding,

\[
 \omega^2=c_{\rm phys}^2|k|^2
 +O(c_{\rm phys}^2a_*^2|k|^4).
 \tag{FD18b}
\]

If `kappa/chi` is instead held fixed, `c_* -> 0`; geometric refinement alone
does not produce a nondegenerate physical-time continuum.  Equation (FD18a),
or an equivalent controller-time rescaling, remains an open dynamics/clock
calibration rather than an F3 consequence.

## 6. Exact refining affine gluing and geometric volume

Define the relational position map on one fixed-depth front by

\[
 X_{a_*}(m)=a_*\sum_am_an_a.
 \tag{FD20}
\]

Sibling differences are the translated roots `a_* alpha_ab`; their transition
maps are base-point independent and every commuting rhombus closes.  The
infinite difference lattice has primitive basis, for example,

\[
 a_*\alpha_{14},\quad a_*\alpha_{24},\quad
 a_*\alpha_{34},
 \tag{FD21}
\]

with exact cell covolume

\[
 \boxed{v_3={16\over3\sqrt3}a_*^3.}
 \tag{FD22}
\]

The convex hull of `X_(a_*)(S_N)` is the regular tetrahedral domain

\[
 \mathcal T_{N,a_*}=\operatorname{conv}\{Na_*n_1,\ldots,Na_*n_4\}
 \tag{FD23}
\]

of volume

\[
 \operatorname{Vol}_3(\mathcal T_{N,a_*})
 ={8(Na_*)^3\over9\sqrt3},
 \qquad
 {\operatorname{Vol}_3(\mathcal T_{N,a_*})\over v_3}
 ={N^3\over6}.
 \tag{FD24}
\]

This agrees with the leading term of the exact front count

\[
 |S_N|={N+3\choose3}
 ={N^3\over6}+N^2+{11N\over6}+1,                  
 \tag{FD25}
\]

with the lower powers constituting the finite boundary census.

For a refinement sequence

\[
 a_*\to0,\qquad N\to\infty,\qquad Na_*\to L>0,
 \tag{FD26}
\]

the translated primitive cells give a shape-regular affine atlas on every
compact subset of the interior of
`conv{L n_1,...,L n_4}`.  The transition maps are translations, their Jacobian
is the identity, and every compact interior set has distance from the boundary
that diverges in **lattice units**, equivalently `min_a m_a -> infinity` along
its corresponding point sequence.  Thus the relational q4 bulk has an exact
refining **mathematical three-manifold** realization; no Cartesian grid is
inserted.

If `tau` is independently bound to one common physical clock and the
collective cone (FD18) is shared by the admitted probes, the corresponding
metric is

\[
 ds^2=-c_*^2d\tau^2+dX^2.
 \tag{FD27}
\]

A controller-time cell `Delta tau` then has coordinate four-volume

\[
 v_4=c_*\Delta\tau\,v_3.
 \tag{FD28}
\]

Equations (FD22)--(FD28) are geometric calibration formulas.  Calling them
actual physical volume still requires the independent RGRL-A clock, probe,
and absolute-volume binding; neither event count nor gamma fixes `a_*`.

## 7. Six pair-memory fields and the exact RGRL-B screening condition

PMMDC's finite statistical tangent is the six-dimensional edge representation

\[
 \mathscr E=A_1\oplus E\oplus T_2,
 \qquad \dim=(1,2,3),
 \tag{FD29}
\]

and its exact Jacobian `D : mathscr E -> Sym^2(V)` is invertible.  This finite
statistical result is not yet a propagating q4/F3 field.  Its use here
therefore requires a second explicit interface:

**`Q4-PAIR-FIELD-LIFT`.**  In one same-parent realization, the four binary
port labels of the PMMDC family are type-identified, equivariantly under
`S4`, with the four reusable q4 operation labels.  Its six tangent variables
`J_(ab)` become compactly supported fields `j_(ab)(m,tau)` on the same
carrier support; their gluing, propagation, constitutive action, retained
sector, and every source/work/controller/boundary/failure/quarantine port are
owned in that parent.  PMMDC and current F3 do not derive this lift.

For completeness, order the six edges as
`(12,13,14,23,24,34)`.  Let `O` exchange opposite edges
`12<->34`, `13<->24`, and `14<->23`, and let `u=(1,...,1)^T`.  The exact
orthogonal sector projectors are

\[
 P_A={uu^{\mathsf T}\over6},\qquad
 P_T={I-O\over2},\qquad
 P_E={I+O\over2}-P_A,
 \tag{FD29a}
\]

with ranks `1,3,2`, respectively, and sum `I_6`.  In the prospective
factorized quadratic collective class under `Q4-PAIR-FIELD-LIFT`, let

\[
 \Chi=\chi_AP_A+\chi_EP_E+\chi_TP_T,
 \quad
 \Kappa=\kappa_AP_A+\kappa_EP_E+\kappa_TP_T,
 \quad
 M^2=m_A^2P_A+m_E^2P_E+m_T^2P_T.
 \tag{FD30}
\]

The three exact sector dispersions are

\[
 \omega_X^2(k)
 ={m_X^2\over\chi_X}
 +{\kappa_X\over\chi_X}\lambda_{A_3}(k),
 \quad
 \lambda_{A_3}(k)=2\sum_{a<b}[1-\cos(a_*k\cdot\alpha_{ab})],
 \tag{FD31}
\]

for every retained `X in R subseteq {A,E,T}`.  Require
`chi_X>0` and `kappa_X>0` for every retained propagating sector.  Then all
retained pair modes share one healthy gapless principal cone in this class
exactly when

\[
 \boxed{
 m_X=0\quad(X\in R),
 \qquad
 {\kappa_X\over\chi_X}=r>0\quad(X\in R).}
 \tag{FD32}
\]

If all three sectors are retained, this is precisely
`kappa_A/chi_A=kappa_E/chi_E=kappa_T/chi_T`.  A projected-out sector imposes
no speed condition on the physical retained modes.

Finite `S4` symmetry alone does not force (FD32).  Nor does it remove four of
the six symmetric-tensor components.  A successful F3 derivation must produce
the common ratio, the constraint/Ward packet leaving the healthy helicity-two
sector, and the complete universal-stress vertex.  More general
orientation-coupled derivative kernels require a correspondingly stronger
principal-symbol audit; (FD30) is a minimal screen, not a classification of
all `S4` actions.

## 8. What has advanced and what remains

This theorem advances the F3 route in four proof-relevant ways.

1. It globally composes q4 common-child incidence into the exact sibling
   kernel `(4I+A_N)` without drawing sibling edges by hand.
2. It proves that an F3-form scalar transfer supplies that kernel through an
   exact Schur complement once `DETUNED-Q4-CARRIER-LIFT` holds; the newly
   supplied detuning and its physical ownership remain explicit.
3. It proves that the fixed q4 root set already has the isotropic second moment
   needed for a smooth collective `3+1` characteristic cone **if** the
   complete massless action is supplied.  Raw direction densification is not
   required on this route.
4. It proves the refining affine bulk atlas and the exact spatial cell-volume
   calibration, while preserving the physical scale/clock ceiling.

The concentrated open physics is now:

\[
 \boxed{
 \begin{gathered}
 \text{derive `DETUNED-Q4-CARRIER-LIFT` in one complete F3 parent;}\\
 \text{derive `Q4-PAIR-FIELD-LIFT` and a stable lineage-conditioned
 massless collective phase; calculate }
 (\Chi,\Kappa,M^2);\\
 \text{derive (or refute) the common-cone condition (FD32), the physical
 constraint/Ward reduction,}\\
 \text{the common clock/probe/volume binding, full matched lineage ports,
 and the RGRL-B stress vertex.}
 \end{gathered}}
 \tag{FD33}
\]

The theorem does not adopt a successor action and does not rescue a failed
parent by adding arbitrary graph machinery.  It turns the existing q4/F3
pieces into one sharply falsifiable dynamics calculation.

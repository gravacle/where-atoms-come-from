# Spatial-deformation and constraint-propagation closure theorem

**Lane ID:** `CROSS-RFT-GRA-EX-SPATIAL-DEFORMATION-CONSTRAINT-PROPAGATION-V001`

**Official short name:** `SDCP`

**Date:** 2026-08-27

**Builder status:** `MUTABLE_HOSTILE_REPAIR__NOT_SEALED`

**Claim class:** exact nonlinear residual-closure theorem on a Gaussian-normal
patch; exact replacement of full ten-component composite-metric tangent span by
six spatial metric deformations plus diffeomorphism Ward identity and initial
constraint custody; conditional composition with the q4 pair-memory Fisher
deformation family and RIEHB

**Not claimed:** physical realization of the pair-memory family; derivation of
the common Lorentzian metric, Gaussian foliation, initial constraints, complete
matter equations, record-to-metric soldering, continuum refinement, a positive
induced Ricci coefficient, gravity in nature, Newton's constant, outcome
selection, singularity resolution, or a black-hole S-matrix

## 1. Exact question

RIEHB gives the exact composite-field equation

\[
 (Dg)^*{\cal E}_g+F_{\rm explicit}=0.              \tag{EX01}
\]

It used dense range of the full metric deformation map as one sufficient way
to infer \({\cal E}_g=0\).  The q4 record construction naturally offers six
pair-memory deformations of a three-dimensional Fisher metric, not ten
independent components of a four-dimensional metric.  Does that necessarily
leave four missing Einstein equations?

No.  For a diffeomorphism-invariant complete effective action, six independent
spatial metric variations give the six evolution equations.  The four normal
projections are constraints.  The nonlinear Ward identity propagates those
constraints homogeneously, so prospectively authenticated zero constraint data
on one Cauchy slice give the remaining four equations everywhere in the
Gaussian-normal development.

The exact replacement for unrestricted dense \(Dg\) is therefore

\[
 \boxed{
 \text{full six-dimensional spatial tangent span}
 +\text{complete Ward identity}
 +\text{zero initial constraint custody}.}         \tag{EX02}
\]

This is a reduction of the gravity gate, not a claim that current record
physics already satisfies it.

## 2. Typed premises

All premises refer to one same-parent retained effective theory.

### SDCP-P1 -- common smooth Lorentzian patch

There is an admitted globally hyperbolic neighborhood on which the already
earned common metric can be written in Gaussian-normal coordinates

\[
 ds^2=-d\tau^2+h_{ij}(\tau,x)dx^idx^j,              \tag{EX03}
\]

with \(h\) positive definite and at least \(C^2\).  This is a gauge choice on
an already earned Lorentzian geometry, not a microscopic background input.
The theorem is local up to the first caustic of the chosen normal congruence;
overlapping admitted patches require compatible constraint custody.

### SDCP-P2 -- complete covariant residual and Ward identity

After every retained matter, electromagnetic, record, writer, reservoir,
support, read, boundary, ghost, and constraint variable has been varied, and
its residual has either been solved or retained explicitly in one enlarged
Noether system, define the symmetric metric residual by

\[
 \delta\Gamma
 ={1\over2}\int d^4x\sqrt{-g}\,
 {\cal E}^{\mu\nu}\delta g_{\mu\nu}.               \tag{EX04}
\]

Before any on-shell reduction, diffeomorphism covariance supplies the complete
Noether identity

\[
 \nabla_\mu {\cal E}^{\mu}{}_{\nu}
 +\sum_B {\cal N}^{B}_{\nu}[{\cal E}_B]
 +{\cal G}_\nu+{\cal B}_\nu=0,                    \tag{EX04a}
\]

where every nonmetric Euler--Lagrange residual is displayed as
\({\cal E}_B\), \({\cal N}^{B}_{\nu}\) is its convention-fixed Noether
differential operator after integration by parts, and
\({\cal G}_\nu,{\cal B}_\nu\) retain all gauge-fixing,
ghost, regulator, anomaly, and boundary terms.  This premise requires those
nonmetric equations actually to be solved and the gauge/regulator/boundary
terms actually to vanish (or, in a separately declared enlarged theorem, to
enter a proved homogeneous coupled propagation system).  Only then does
(EX04a) reduce to the complete on-shell identity

\[
 \boxed{\nabla_\mu{\cal E}^{\mu\nu}=0.}             \tag{EX05}
\]

Merely listing another retained field in an enlarged residual is not enough:
its residual must vanish or its source term must remain.  Boundary anomalies,
explicit backgrounds, and unmatched controller forces are absent or are
retained in (EX04a) and solved by the enlarged theorem.

### SDCP-P3 -- six same-parent spatial deformation channels

There are six retained collective spacetime fields
\(J_A(\tau,x)\), \(A=1,\ldots,6\), admitting arbitrary compactly supported
variations in the declared patch.  Their only relevant local dependence in
the declared sector is through the spatial metric.  Any apparent explicit
force must vanish, or factor through other retained Euler--Lagrange residuals
that have separately been solved to zero.  At every point in the admitted
patch,

\[
 M_{Aij}:={\partial h_{ij}\over\partial J_A}         \tag{EX06}
\]

spans \(\operatorname{Sym}^2(T^*\Sigma)\):

\[
 \boxed{\operatorname{rank}\{M_{Aij}\}=6.}          \tag{EX07}
\]

For cellwise or nonlocal collective variables, (EX07) is replaced by
injectivity of the adjoint deformation map on the spatial residual in a
declared function space, with a controlled refinement whose admissible
variations approximate arbitrary compactly supported spatial metric
variations.  A finite six-parameter vector and pointwise matrix rank may not
be substituted for that functional-range statement.

### SDCP-P4 -- prospective initial constraint custody

On one Cauchy surface \(\Sigma_0\) fixed before the scored evolution, the four
normal residuals vanish:

\[
 \boxed{{\cal E}^{00}|_{\Sigma_0}=0,\qquad
        {\cal E}^{0i}|_{\Sigma_0}=0.}               \tag{EX08}
\]

These conditions must follow from the same-parent initial/boundary equations,
admissibility law, or independently checked constraints.  Observing a later
Einstein solution and then selecting (EX08) is not custody.

The residual components and coefficients must also lie in a uniqueness class
for the homogeneous propagation equations below.  For the displayed classical
proof it is sufficient that \({\cal E}^{0\nu}\) be continuous and
\(K^i{}_j\) locally integrable in \(\tau\) along each normal curve.  A weak or
higher-derivative realization must state and prove its corresponding
distributional uniqueness theorem.

## 3. Theorem SDCP-1 -- six pair deformations close the spatial residual

Assume SDCP-P1--P3 and stationarity under every compactly supported
\(\delta J_A\).  Then

\[
 \boxed{{\cal E}^{ij}=0}                             \tag{EX09}
\]

throughout the admitted patch.

### Proof

At fixed Gaussian lapse and shift, \(\delta g_{00}=\delta g_{0i}=0\) and
\(\delta g_{ij}=M_{Aij}\delta J_A\).  The chain rule, arbitrary compactly
supported \(\delta J_A\), and the
zero/factorized-explicit-force premise give

\[
 0={\delta\Gamma\over\delta J_A}
 ={1\over2}\sqrt{-g}\,{\cal E}^{ij}M_{Aij}          \tag{EX10}
\]

in the local realization, with the corresponding adjoint equation in the
functional realization.  Since the six \(M_A\) span the six-dimensional space
of symmetric spatial tensors, their annihilator is zero.  Hence
\({\cal E}^{ij}=0\).  QED.

The conclusion is nonlinear: no expansion \(g=\eta+h\) has been used.

## 4. Theorem SDCP-2 -- Ward propagation supplies the four constraints

Assume SDCP-P1--P4 and (EX09).  Then

\[
 \boxed{{\cal E}^{\mu\nu}=0}                        \tag{EX11}
\]

throughout the Gaussian-normal development of \(\Sigma_0\).

### Proof

Let \(h=\det h_{ij}\) and

\[
 K^i{}_j={1\over2}h^{ik}\partial_\tau h_{kj}.       \tag{EX12}
\]

For a symmetric contravariant tensor,

\[
 \nabla_\mu{\cal E}^{\mu\nu}
 ={1\over\sqrt h}\partial_\mu
   (\sqrt h\,{\cal E}^{\mu\nu})
 +\Gamma^\nu_{\mu\lambda}{\cal E}^{\mu\lambda}.
                                                               \tag{EX13}
\]

In Gaussian-normal coordinates,
\(\Gamma^j{}_{0i}=K^j{}_i\),
\(\Gamma^j{}_{00}=\Gamma^0{}_{00}=\Gamma^0{}_{0i}=0\).
Using (EX09), the spatial components \(\nu=j\) of (EX05) reduce exactly to

\[
 \boxed{
 \partial_\tau(\sqrt h\,{\cal E}^{0j})
 =-2\sqrt h\,K^j{}_i{\cal E}^{0i}.}                 \tag{EX14}
\]

At each spatial point this is a homogeneous linear ordinary differential
system.  Its unique solution with the initial data (EX08) is
\({\cal E}^{0i}=0\).

The normal component \(\nu=0\) then reduces to

\[
 \partial_\tau(\sqrt h\,{\cal E}^{00})
 +\partial_i(\sqrt h\,{\cal E}^{i0})=0.            \tag{EX15}
\]

The second term vanishes because \({\cal E}^{i0}=0\), so the initial value in
(EX08) gives \({\cal E}^{00}=0\).  Together with (EX09), this proves
(EX11). QED.

Equation (EX14) is the exact nonlinear reason that four additional microscopic
metric slots are not automatically required.  It does not remove the physical
obligation to earn (EX05) and (EX08).

## 5. Corollary SDCP-3 -- reduced RIEHB closure

Suppose the same parent also satisfies RIEHB's common-metric, complete
Wilsonian census, positive matched Ricci coefficient, controlled derivative
remainder, and explicit-force premises.  Write its complete residual as

\[
 {\cal E}_{\mu\nu}
 =G_{\mu\nu}+\Lambda_{\rm eff}g_{\mu\nu}
 -8\pi G_{\rm eff}T^{\rm complete}_{\mu\nu}
 -\Delta^{\rm rem}_{\mu\nu},                       \tag{EX16}
\]

in units with the declared speed factors absorbed, and require
\(G_{\rm eff}>0\).  Under SDCP-P1--P4,

\[
 \boxed{
 G_{\mu\nu}+\Lambda_{\rm eff}g_{\mu\nu}
 =8\pi G_{\rm eff}T^{\rm complete}_{\mu\nu}
 +\Delta^{\rm rem}_{\mu\nu}.}                     \tag{EX17}
\]

Thus RIEHB's unrestricted dense-\(Dg\) condition can be replaced, on this
foliated same-parent domain, by the weaker injectivity condition

\[
 \boxed{
 \ker((D_\Sigma g)^*)
 \cap\{\text{Ward-compatible residuals with zero initial constraints}\}
 =\{0\}.}                                           \tag{EX18}
\]

SDCP-1--2 prove (EX18) when the six spatial fields have the required local or
functional range.  A positive leading Einstein--Hilbert term supplies the
Fierz--Pauli/helicity-two kinetic operator about an admitted solution.  An
actual protected retarded massless pole still requires low-energy dominance,
no cancellation by the controlled remainder, hyperbolicity, causal boundary
conditions, and suitable local or asymptotic background control.

## 6. Pair-memory Fisher composition

The pending EW lane supplies the focused candidate for SDCP-P3.  On four
binary record ports \(s_a=\pm1\), use tetrahedral contrast vectors
\(v_a\in V=\mathbf1^\perp\) and

\[
 p_{J,\theta}(s)
 ={1\over Z(J,\theta)}
 \exp\!\left[\theta\mathbin\cdot
       \sum_av_as_a+\sum_{a<b}J_{ab}s_as_b\right].  \tag{EX19}
\]

At \(\theta=0,J=0\), the localization Fisher tensor for \(\theta\) is
isotropic.  The six pair couplings obey

\[
 {\partial {\cal F}_{\theta}\over\partial J_{ab}}
 \bigg|_{0}=v_a\odot v_b,\qquad
 v_a\odot v_b:=v_av_b^{\mathsf T}+v_bv_a^{\mathsf T}. \tag{EX20}
\]

The six tetrahedral symmetric products form a basis of
\(\operatorname{Sym}^2(V)\).  Therefore, if the same complete record query is
physically soldered pointwise by

\[
 h=\ell_F^2{\cal F}_\theta,                         \tag{EX21}
\]

the inverse-function theorem supplies the six-dimensional matrix rank at one
cell on an open neighborhood of the symmetric point.  Global-flip symmetry
keeps every one-port marginal uniform, so these are genuine
pair-memory/correlation deformations rather than hidden single-port loading
changes.

Equations (EX19)--(EX21) are a same-parent mathematical deformation mechanism.
They satisfy the algebraic part of SDCP-P3 only after EW's exact factors,
complete-query type join, positivity domain, and record qualification survive
independent audit.  To satisfy the actual functional premise, the parent must
also supply one six-field copy \(J_A(\tau,x)\), or a cellwise refinement whose
admissible variations become dense in the declared compactly supported
function space.  For every such variation it must recompute all ten signed
simplex intervals and the lapse/null Gram data, preserve Lorentz signature,
nondegeneracy, the independently fixed clock bind, and exact shared-edge
gluing.  The finite EW inverse-function theorem alone proves none of this local
field lift.  It also does not prove that the actual world uses this family or
that its pair variables are physically soldered to the common metric.

## 7. Controls and countermodels

1. **Rank control.**  If the six \(M_A\) have rank below six, a nonzero spatial
   residual in their annihilator satisfies every collective equation.
2. **Single-family control.**  A q4 label count or one three-parameter softmax
   family does not imply (EX07).  Dimension of a localization tangent and
   dimension of its metric-deformation space are different questions.
3. **Ward control.**  Without (EX05), setting \({\cal E}^{ij}=0\) places no
   condition on \({\cal E}^{0\mu}\).
4. **Initial-constraint control.**  Even with (EX05), a nonzero constant
   \(\sqrt h{\cal E}^{00}\) in a static patch survives when
   \({\cal E}^{0i}=0\).  Equation (EX08) is load bearing.
5. **Matter-shell control.**  A metric Ward identity with off-shell matter
   source terms is not (EX05).  The complete enlarged residual must be solved.
6. **Explicit-force control.**  A direct \(J_A\) force can cancel the pullback
   metric residual in (EX10).  It must vanish; a factorized force is removable
   only when every factor residual has separately been solved to zero.
7. **Gauge control.**  Gaussian-normal coordinates are chosen only after a
   common Lorentzian metric is earned.  Using them to define the microscopic
   record family would violate E-EMERGENTSPACE.
8. **Patch control.**  Caustics or unmatched patch boundaries terminate the
   local theorem; they cannot be crossed by notation.
9. **Induced-action control.**  SDCP closes a variational rank/constraint gate.
   It does not prove \(C_R^{\rm eff}>0\), common-cone field content, or remainder
   control.
10. **Ancestry control.**  A pair-memory family and a metric in a direct product
    can satisfy their separate equations.  Same-parent soldering and a matched
    lineage intervention remain necessary for a record-origin claim.
11. **Functional-lift control.**  Six global parameters can annihilate six
    integrated projections while leaving a nonzero local spatial residual.
    Cellwise/local fields or the declared adjoint-injectivity and refinement
    theorem are load bearing.

## 8. Exact disposition

This lane proves that the final induced-gravity theorem does not need ten
independent record deformation channels.  Six full-rank spatial metric
deformations, complete diffeomorphism Ward custody, and zero initial
constraints imply the complete nonlinear metric equation.  The q4
pair-memory Fisher family is the first exact candidate for those six channels;
its physical soldering and actual-world realization remain open.

**Disposition:**

`FULL_SIX_SPATIAL_DEFORMATION_RANK_IMPLIES_ALL_SPATIAL_METRIC_RESIDUALS_ZERO__COMPLETE_NONLINEAR_WARD_IDENTITY_PROPAGATES_FOUR_ZERO_INITIAL_CONSTRAINTS__FULL_METRIC_EQUATION_FOLLOWS_WITHOUT_TEN_INDEPENDENT_MICROSCOPIC_METRIC_SLOTS__PAIR_MEMORY_FISHER_FAMILY_IS_CANDIDATE_SAME_PARENT_SPATIAL_TANGENT_COMPLETION__PHYSICAL_SOLDERING_INITIAL_CONSTRAINT_CUSTODY_POSITIVE_INDUCED_COEFFICIENT_REFINEMENT_ANCESTRY_AND_REAL_WORLD_INSTANTIATION_OPEN`

# Record-conditioned induced Einstein--Hilbert back-reaction theorem

**Lane:** `LANE_CROSS_RFT_MGFT_INDUCED_EH_BACKREACTION_V001`

**Short name:** `RIEHB`

**Date:** 2026-08-26

**Claim class:** exact conditional finite-band Wilsonian heat-kernel theorem;
exact composite-metric variational-span gate; full nonlinear
leading-derivative metric action on an already-earned smooth domain

**Not claimed:** derivation of a manifold, dimension, Lorentzian signature,
common cone, a metric from gamma alone, positivity or numerical calculation of
Newton's constant from current F3 data, cancellation of the cosmological term,
an all-scale exact action, a UV completion, singularity resolution, or a proof
that nature realizes the parent

## 1. Physical question

The record program does not expect microscopic records themselves to be tiny
pieces of classical gravity.  The working sequence is instead

\[
 \text{record-conditioned relational phase}
 \longrightarrow g_{\mu\nu}^{\rm eff}
 \longrightarrow \Gamma_{\rm eff}[g,\text{slow fields}]
 \longrightarrow \text{metric back-reaction}.
 \tag{IB01}
\]

Suppose the first arrow has been earned.  Does ordinary collective quantum
field physics supply a concrete mechanism for the second and third arrows?

Yes, conditionally.  Once fast physical modes propagate on one common smooth
Lorentzian metric, their covariant functional determinants admit the Ricci-
scalar invariant and generically generate it.  Species and nonminimal-coupling
contributions can cancel, so a nonzero positive coefficient remains an
explicit gate.  When present, that term is the full nonlinear
Einstein--Hilbert functional, not merely its linearization.  The full Einstein
equation, however, follows from the underlying collective equations only if
the emergent metric variations have complete physical tangent span and no
independent collective force remains.

## 2. Frozen parent `IBP`

Everything below belongs to one same-parent coarse-graining family.  No metric
or action from a separate model may be joined after seeing the result.

### IBP.1 -- earned common geometry

On a relational domain `D`, prior record-conditioned dynamics has earned a
smooth four-dimensional Lorentzian metric `g` and a common low-energy cone for
the admitted matter, electromagnetic, record, clock, and probe modes.  A
Euclidean continuation or an equivalent causal determinant prescription is
available for the local coefficient calculation.  Boundary conditions and
edge modes are fixed prospectively.

This is the output sought from the active support/locality/cone lane.  It is a
premise here, not a result of the determinant.

### IBP.2 -- Wilsonian fast/slow split

At physical coarse-graining scales `mu<kappa_R`, split the complete same-parent
field content into a fast shell `chi_>` and retained slow modes `chi_<`.  The
split is covariant on the common metric, and every determinant, ghost,
constraint, boundary, and measure contribution is counted exactly once.  A
sharp proper-time representative uses only
`s in [kappa_R^-2,mu^-2]`; a smooth covariant profile is also admissible if its
moments are owned.  The complementary infrared/nonlocal contribution belongs
to `Gamma_<` or to one separately owned remainder and is not integrated again.
The fast Hessians
are Laplace type after the declared gauge/constraint treatment,

\[
 P_A(g)=-\nabla_g^2+{\cal E}_A(g)+m_A^2,
 \tag{IB02}
\]

or lie in a separately proved generalized heat-kernel class.  Statistics and
multiplicity are carried by signed weights `sigma_A`.

The physical ultraviolet end of the emergent phase is `kappa_R`; it is the
scale where this metric description or the record-conditioned collective
phase ceases to apply.  It is not silently called the Planck scale.

### IBP.3 -- covariant derivative expansion

For curvatures, masses/potentials relevant to the expansion, and retained
momenta small compared with `mu`, the fast-shell functional determinant admits
the local asymptotic expansion

\[
 \operatorname{Tr}e^{-sP_A(g)}
 \sim {1\over(4\pi s)^2}\int_D\!d^4x\sqrt{|g|}\,
 \operatorname{tr}\left[a_{0,A}+s a_{1,A}+s^2a_{2,A}+\cdots\right].
 \tag{IB03}
\]

General covariance makes `a_0` a scalar constant, `a_1` a linear combination
of `R` and mass/potential scalars, and `a_2` a sum of curvature-squared,
field-strength, potential, and total-derivative invariants.

### IBP.4 -- physical coefficient, origin bookkeeping, and retained source

Keep the pre-existing and fast-mode contributions separate,

\[
 C_R^{\rm eff}(\mu)=C_R^{\rm pre}(\mu)+C_R^{>}(\mu).
 \tag{IB04a}
\]

After all fast species, measures, and matching terms assigned to this
Wilsonian parent are included, the determinant contribution `C_R^>(mu)` is
finite.  It can vanish by cancellation.  For the strict **induced-origin**
conclusion require a physically distinguished microscopic
regulator/matching prescription, fixed by the parent rather than by a later
renormalization convention, in which `C_R^pre=0` by an exact microscopic
condition or symmetry and `C_R^> != 0`.  Without that extra physical
prescription, finite local terms can be shifted between preterm, matching, and
determinant, so this lane proves only an induced contribution/renormalization
and the matched total.  For healthy back-reaction require the total coefficient
`C_R^eff` of

\[
 \int_D d^4x\sqrt{-g}\,R
\]

to be positive in the Lorentzian convention

\[
 C_R^{\rm eff}={1\over16\pi G_{\rm eff}}>0.
 \tag{IB04}
\]

The retained slow action defines one complete total stress tensor

\[
 T_{\mu\nu}^{<}:=-{2\over\sqrt{-g}}
 {\delta\Gamma_{<}\over\delta g^{\mu\nu}},
 \tag{IB05}
\]

including matter, EM, record, writer, support, reservoir, work, read, and
boundary terms at the claimed scale exactly once.

`IBP.4` is an explicit origin, sign, and normalization gate.  The heat-kernel
operator basis alone does not guarantee a nonzero positive coefficient for an
arbitrary spectrum.  A freely adjustable unowned counterterm would destroy a
numerical prediction.

## 3. Theorem RIEHB-1 -- the finite fast shell generates the local gravity basis

Define the regulated one-loop difference relative to a fixed reference metric
`g_0` by

\[
 \Delta\Gamma_>[g;g_0]
 =-{\hbar\over2}\sum_A\sigma_A
 \int_{\kappa_R^{-2}}^{\mu^{-2}}{ds\over s}
 \operatorname{Tr}\!\left(e^{-sP_A(g)}-e^{-sP_A(g_0)}\right),
 \tag{IB06}
\]

with the overall signed convention absorbed consistently into `sigma_A`.
Substitution of (IB03) gives

\[
 \Delta\Gamma_>
 =\int_D\!d^4x\sqrt{-g}\left[
 C_0^{>}(\mu)+C_R^{>}(\mu)R+
 C_{R^2}^{>}(\mu)R^2+C_{C^2}^{>}(\mu)C_{\alpha\beta\gamma\delta}^2
 +\cdots\right]-(g\to g_0).
 \tag{IB07}
\]

For the sharp massless shell the two leading proper-time integrals are exactly

\[
 \int_{\kappa_R^{-2}}^{\mu^{-2}} ds\,s^{-3}
 ={\kappa_R^4-\mu^4\over2},
 \qquad
 \int_{\kappa_R^{-2}}^{\mu^{-2}} ds\,s^{-2}
 =\kappa_R^2-\mu^2.
 \tag{IB08}
\]

With mass, the second expression becomes

\[
 I_1(m,\kappa_R,\mu)=
 \int_{\kappa_R^{-2}}^{\mu^{-2}} ds\,s^{-2}e^{-m^2s}
 =\kappa_R^2-\mu^2-m^2\log(\kappa_R^2/\mu^2)
 +O\!\left(m^4\mu^{-2}\right)
 \tag{IB09}
\]

in the declared `m^2/mu^2 -> 0` light-threshold asymptotic convention; for
general mass the integral on the first line is retained exactly.  Hence
`C_R^>` is a
calculable signed spectral/threshold sum once the complete fast spectrum,
nonminimal couplings, physical cutoff/crossover, measure, and matching
terms are owned.

Under `IBP.4`, include any owned pre/matching volume term and write

\[
 C_0^{\rm eff}:=C_0^{\rm pre}+C_0^{>}
 =-{\Lambda_{\rm eff}\over8\pi G_{\rm eff}}.
\]

Then the leading matched two-derivative metric action is

\[
 \boxed{
 \Gamma_{\rm grav}^{(0,2)}[g]
 ={1\over16\pi G_{\rm eff}}
 \int_D d^4x\sqrt{-g}\,(R-2\Lambda_{\rm eff}).}
 \tag{IB10}
\]

Equation (IB10) is nonlinear in `g` to all orders.  It is only an expansion in
derivatives/curvature and loops.  At momenta `k` where

\[
 \left|{C_{R^2}k^2\over C_R^{\rm eff}}\right|\ll1,
 \qquad
 \left|{C_{C^2}k^2\over C_R^{\rm eff}}\right|\ll1,
 \tag{IB11}
\]

the Einstein--Hilbert term is the leading metric dynamics.

### Proof

Insert (IB03) into (IB06), integrate term by term on the declared asymptotic
domain, and group the resulting diffeomorphism scalars by derivative order.
The `a_0` integral gives (IB08)'s quartic volume coefficient; the curvature
part of `a_1` gives the quadratic Ricci-scalar coefficient; `a_2` gives
logarithmic curvature-squared terms.  No expansion of `g=eta+h` was used, so
the scalar `sqrt(-g)R` is the full nonlinear functional.  Positivity and the
definition (IB04) then give (IB10), while power counting gives (IB11).  QED.

## 4. Theorem RIEHB-2 -- back-reaction and one common leading coefficient

If `g` is an independently variable effective field and boundary variations
are fixed or cancelled by the admitted boundary action, stationarity of

\[
 \Gamma_{\rm eff}=\Gamma_{\rm grav}^{(0,2)}+\Gamma_<+\Gamma_{\rm rem}
\]

gives

\[
 \boxed{
 G_{\mu\nu}+\Lambda_{\rm eff}g_{\mu\nu}
 =8\pi G_{\rm eff}T_{\mu\nu}^{<}
 +\Delta_{\mu\nu}^{\rm rem}.}
 \tag{IB12}
\]

Here `Delta_rem` is the separately bounded variation of curvature-squared,
nonlocal, higher-loop, boundary, and matching terms.  One coefficient
multiplies the complete variational total stress because all admitted slow
modes enter through the same metric variation (IB05), not because recordhood
gives each species a separately tuned charge.  Diffeomorphism invariance gives
the corresponding Ward/Bianchi compatibility for the complete on-shell
parent.  This leading common coefficient does not by itself prove the weak or
strong equivalence principle, exclude nonminimal curvature couplings, or
remove forces carried by additional retained backgrounds.

Equation (IB12) is the concrete macroscopic feedback loop:

\[
 \text{fast collective fluctuations create metric stiffness}
 \longrightarrow T^< \text{ deforms }g
 \longrightarrow g \text{ changes future matter/EM/record propagation}.
 \tag{IB13}
\]

It does not say that gamma itself is a stress tensor or that information
curvature is algebraically spacetime curvature.

## 5. Theorem RIEHB-3 -- exact composite-metric rank gate

In a genuine emergent theory the metric may not be independently variable.
Let the retained collective backgrounds be `Phi^A`, with

\[
 g_{\mu\nu}=g_{\mu\nu}(\Phi),
 \qquad
 Dg:\delta\Phi\mapsto\delta g.
 \tag{IB14}
\]

The exact chain rule is

\[
 0={\delta\Gamma\over\delta\Phi^A}
 =(Dg)^*_{A}{}^{\mu\nu}{\cal E}_{\mu\nu}+F_A^{\rm explicit},
 \qquad
 F_A^{\rm explicit}:=
 \left.{\delta\Gamma\over\delta\Phi^A}\right|_g,
 \tag{IB15a}
\]

where

\[
 {\mathcal E}_{\mu\nu}:={1\over\sqrt{-g}}
 {\delta\Gamma\over\delta g^{\mu\nu}}.
\]

An owned but nonzero `F_explicit` cannot be dropped.  For the pure metric
equation require it to vanish, or prove that all apparent explicit dependence
factors through an enlarged effective-field vector whose residuals are
included in `E`.  Under that condition (IB15a) reduces to

\[
 (Dg)^*{\cal E}=0.                                  \tag{IB15}
\]

Therefore:

1. `(Dg)^* E=0` does **not** generally imply `E=0`.
2. The exact Hilbert-space condition is `ker((Dg)^*)={0}`, equivalently that
   the range of `Dg` is dense in the complete physical metric-variation space
   (symmetric variations modulo gauge/boundary nulls).  Then (IB15) implies
   the full effective metric equation `E=0`.
3. Surjectivity of `Dg` is sufficient.  In a finite-mode truncation it is
   exactly full row rank.  If the infinite-dimensional range is closed,
   density and surjectivity coincide.
4. If the range is not dense, its nonzero orthogonal complement supplies a
   residual invisible to every allowed collective variation.  The parent
   earns only projected metric dynamics.  Mere nonsurjectivity is not enough
   for this conclusion when the range is dense but nonclosed.

### Proof

Equation (IB15a) is the functional chain rule.  After the explicit-force gate,
(IB15) is its pure pullback.  In finite mode truncation it is the elementary
linear-algebra statement `D^T E=0`; full row rank makes `D^T` injective.  In a
Hilbert realization, `ker(D^*)` is the orthogonal complement of the closure of
`range(D)`, proving the dense-range criterion and the non-dense countermodel.
The declared Banach realization must supply the corresponding separating-dual
condition rather than import Hilbert orthogonality by name.  QED.

This rank gate is the exact repair to the tempting claim that the determinant
alone derives full Einstein dynamics for any emergent metric.  It also gives a
focused target for F3: derive enough independent collective deformation
channels to span the physical metric response, rather than adding six tensor
bits to every microscopic record.

## 6. Record and gamma ancestry

The theorem composes with RFT/MGFT only if all of the following are same-parent:

1. retained records causally contribute to forming or maintaining the
   relational phase that earns `g`;
2. the complete common-cone field content in `IBP.2` is the content of that
   same phase;
3. `kappa_R` is its measured crossover/breakdown scale;
4. KEEP versus whole-lineage BREAK, with energy/stress/EM/work/boundary ports
   matched, changes at least one of the phase metric, the fast spectral sum,
   or the induced response coefficient;
5. the tangent-span map and explicit-force gate in `RIEHB-3` belong to the same
   collective parent.

Gamma can diagnose formation and, in qualified equilibrium windows, help
reconstruct response spectra.  It does not replace the determinant's complex
phase, the common-cone premise, the tangent-span/explicit-force gates, or the
sign and scale matching in (IB04).

## 7. Exact advance and remaining frontier

This theorem closes a previously missing **mechanism class** from earned
geometry to gravitational back-reaction:

\[
 \boxed{
 \text{common emergent Lorentzian metric}
 +\text{covariant fast field spectrum}
 +C_R^{\rm eff}>0
 +F_{\Phi}^{\rm explicit}=0
 +\ker(Dg)^*=\{0\}
 \Longrightarrow
 \text{full nonlinear EH back-reaction at leading derivative}.}
 \tag{IB16}
\]

It shifts the concentrated gravity-origin work to measurable/derivable
objects: the record-conditioned stable local phase, its common Lorentzian
metric, the complete fast spectrum and crossover, the sign/normalization of
`C_R`, the independent collective-force residual, and the physical tangent
span of its deformation map.  It does not close those objects by assertion.

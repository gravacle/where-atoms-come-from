# Protected massless spin-2 completion audit

Statements are labeled **EXACT**, **CONDITIONAL**, **EXCLUSION**, or **OPEN**. The
alternatives and falsifiers were sealed before this calculation. This lane neither
selects a theory nor assigns empirical or proof weight.

## 1. Increment beyond the existing threshold result

GRA-AS established a coherent generic threshold

\[
 \det[I-\chi^R K^R]=0,
 \qquad
 \Xi=\lambda_{\max}(\chi^R K^R)=1
\]

in its restricted static form. A protected `1/k^2` response can yield a `1/r`
potential in three spatial dimensions. That result did not determine the pole's
Lorentz representation, numerator, source, residue, gauge identity, or nonlinear
vertices.

The missing implication is therefore

```text
generic protected massless pole
    != protected massless helicity-2 graviton.
```

This lane tests what additional structure would make the implication valid.

## 2. What representations can record order supply?

Work locally in a record rest frame only for representation counting. At nonzero
spatial momentum, rotations decompose candidate order variables as follows.

| candidate | linear SO(3) content | helicity-2 available linearly? |
|---|---|---|
| scalar `X` | one scalar | no |
| four-vector `V_mu` | two scalar pieces plus two transverse-vector components | no |
| symmetric `Q_mu_nu` | four scalar, four vector and two tensor components | yes, but not isolated or protected |

For `Q_mu_nu`, the ten components split schematically into

\[
 Q_{00}:1_S,
 \quad Q_{0i}:1_S+2_V,
 \quad Q_{ij}:2_S+2_V+2_T.                              \tag{AT01}
\]

**EXACT representation exclusion.** A scalar or one vector cannot be the linear
carrier of a helicity-2 pole on an isotropic background. A bilinear such as
`V_mu V_nu` is a composite operator and belongs to `AT-A2`, not to a linear vector
mode. A nonzero vector expectation also selects a preferred direction.

An isotropic symmetric background has the form

\[
 \bar Q_{\mu\nu}=q_g\eta_{\mu\nu}
 +q_u{u_\mu u_\nu\over c^2}.                            \tag{AT02}
\]

`q_u != 0` preserves spatial rotations in the record frame but breaks boost
invariance. The tensor fluctuation in (AT01) may still exist, but `c_T=c`, a common
matter cone, and Lorentzian soft-theorem premises then require an explicit infrared
restoration/decoupling result.

**EXACT limitation.** Merely possessing two spatial transverse-traceless components
does not make `Q_mu_nu` a massless graviton. A generic symmetric-tensor action has
additional scalar/vector components and no identity forcing them to be constraints.
The tensor order parameter must acquire a gauge redundancy or an equivalent complete
degenerate constraint algebra.

## 3. The required linear massless spin-2 structure

Let `h_mu_nu` be a dimensionless symmetric fluctuation. In a locally Lorentzian,
local, two-derivative domain, the massless Fierz-Pauli kinetic operator may be written
compactly as

\[
 S_2={Z_2\over2}\int d^4x\,h_{\mu\nu}
       {\cal E}^{\mu\nu\rho\sigma}h_{\rho\sigma}
 -{1\over2c}\int d^4x\,h_{\mu\nu}T^{\mu\nu},           \tag{AT03}
\]

where normalization factors may be moved between `h`, `Z_2`, and the coupling, but
the physical pole residue may not. The operator obeys

\[
 \partial_\mu {\cal E}^{\mu\nu\rho\sigma}h_{\rho\sigma}=0,
 \qquad
 \delta h_{\mu\nu}=\partial_\mu\xi_\nu+\partial_\nu\xi_\mu. \tag{AT04}
\]

Gauge invariance of the source term requires

\[
 \partial_\mu T^{\mu\nu}=0.                             \tag{AT05}
\]

Four first-class gauge generators and their constraints reduce ten configuration
components to two propagating helicities. The positive physical spin-2 residue
requires the correctly signed `Z_2`; a mass term or generic potential breaks (AT04)
and changes the degree count.

**CONDITIONAL uniqueness.** Under local Lorentz invariance, a single symmetric
field, at most two derivatives, and exactly two massless helicities, the healthy
quadratic theory is in the massless Fierz-Pauli/linearized-Einstein family up to
normalization, field redefinition, boundary terms, and gauge fixing. This statement
does not apply to nonlocal, Lorentz-breaking, higher-derivative, or multi-field
degenerate systems.

The static Newtonian constraint is not itself a radiative scalar. Conserved sources
excite constrained `h_00`/longitudinal components while the on-shell radiation
contains the two helicity-2 states. A TT projection of a static monopole would reject
GR itself and is not used here.

## 4. Why the long-range coupling must be universal

Under a Lorentz-invariant S-matrix with a massless spin-2 pole, the leading soft
emission factor from external leg `i` has the form

\[
 {\cal M}_{n+1}\longrightarrow
 \left[\sum_i\eta_i g_i
 {p_i^\mu p_i^\nu\over p_i\!\cdot q}\right]
 \epsilon_{\mu\nu}(q)\,{\cal M}_n,                     \tag{AT06}
\]

where `eta_i` distinguishes incoming and outgoing momenta. The unphysical shift

\[
 \epsilon_{\mu\nu}\rightarrow\epsilon_{\mu\nu}
 +q_\mu\zeta_\nu+q_\nu\zeta_\mu                       \tag{AT07}
\]

decouples only if

\[
 \sum_i\eta_i g_i p_i^\nu=0.                            \tag{AT08}
\]

Ordinary momentum conservation supplies `sum eta_i p_i=0`. For arbitrary processes
and species, (AT08) follows from it only when all nonzero long-range charges share
one coupling `g_i=g_*`, including the spin-2 field's own energy-momentum.

**CONDITIONAL result.** This is the Weinberg-style equivalence/universal-coupling
constraint. It is stronger than composition-independent static acceleration: it
requires the same normalized coupling to matter, antimatter, radiation, binding
energy, clocks, probes, and the gravitational field itself. It assumes the stated
soft Lorentzian particle/S-matrix regime; it does not prove that the regime emerges
from records.

## 5. Nonlinear bootstrap and the assumption-scoped IR family

Once `h_mu_nu` couples to stress, its own stress changes the source. Iterating the
same coupling while preserving the gauge identity deforms both the action and gauge
transformations. Under the standard assumptions

```text
one massless spin-2 field;
locality and Lorentz/Poincare invariance;
unitary positive kinetic sign;
at most two derivatives in the leading equations;
consistent coupling to the complete conserved stress tensor;
no independent long-range spin-2 charge,
```

the deformation closes into the Einstein-Hilbert family with a possible
cosmological term, up to field redefinitions and boundary terms:

\[
 S_{IR}={c^3\over16\pi G_*}\int d^4x\sqrt{-g}(R-2\Lambda)
       +S_{all}[g,\Psi].                                 \tag{AT09}
\]

Four-dimensional Lovelock uniqueness independently says that a symmetric,
divergence-free natural metric tensor made from the metric and at most its first two
derivatives lies in the `G_mu_nu + Lambda g_mu_nu` family under the theorem's
naturality, regularity, and quasilinearity assumptions.

**CONDITIONAL unique-family result.** The endpoint is not an arbitrary landscape:
within the frozen leading-order assumptions, the admissible family is GR plus
`Lambda`, normalization, and diffeomorphism-invariant higher-scale corrections. This
is an IR-family restriction, not a record derivation and not canonical selection.

The result does not exclude scalar-tensor theories, higher-curvature EFT terms,
nonlocal actions, Lorentz-breaking phases, extra dimensions, or other additional
fields. Those evade at least one uniqueness premise and must separately pass the
extra-mode, local-GR, stability, and common-metric gates. Local finite-order `C^2`
terms generically expose a wrong-sign spin-2 pole if treated as a fundamental theory;
as EFT operators their extra pole must stay outside the admitted cutoff.

For several initially independent massless Fierz-Pauli fields with positive internal
metric, locality and at most two derivatives, the consistent-deformation results
exclude Yang-Mills-like cross-interactions; the allowed deformation is a sum of
separate Einstein-Hilbert sectors under those assumptions. Thus a record tensor plus
an already-existing metric cannot simply be declared two interacting universal
gravitons. The physical long-range matter coupling must reduce to one admitted
metric combination or enter a theory outside the theorem's domain with its complete
constraint audit.

## 6. Ward/Bianchi closure is an off-shell requirement

For a diffeomorphism-invariant effective action `Gamma[g,Psi_R,Psi_known]`, the
Noether identity has the schematic form

\[
 2\nabla_\mu\!\left({1\over\sqrt{-g}}
 {\delta\Gamma\over\delta g_{\mu\nu}}\right)
 +\sum_A E_A\,\nabla^\nu\Psi_A
 +\text{representation terms}=0.                        \tag{AT10}
\]

The representation terms are essential when record order is vector or tensor.
Only after every exchanging field satisfies its equation does (AT10) reduce to the
on-shell Bianchi/source relation. Imposing `nabla T=0` beside an otherwise incomplete
filter or tensor equation is not an off-shell construction.

The complete source must be obtained once from metric variation:

```text
T_total = T_known + T_EM + T_record + T_reservoir
          + T_interaction + T_support/boundary,
```

with the partition adjusted to avoid overlap. If a record contribution is moved to
the geometric Euler tensor, it may not also remain on the source side. For an open
or dissipative record sector, a retarded in-in influence action or local enlarged
system must carry the corresponding noise/supply and Ward identity; a static kernel
does not provide them.

## 7. Protection through continuing record accumulation

Let a tensor inverse propagator near the threshold contain

\[
 D_2^{-1}(\omega,k;X)=Z_T(X)
 [-(\omega+i0)^2/c_T^2(X)+k^2+m_T^2(X)] +\cdots.         \tag{AT11}
\]

A generic critical transition gives `m_T^2(X_c)=0` at one point. If records continue
to accumulate, `X` moves away and the range becomes finite. Long-range gravity in
the whole post-threshold phase requires

\[
 m_T^2(X)=0\quad\text{for every admitted }X\ge X_c,       \tag{AT12}
\]

protected by an exact/emergent gauge redundancy, or a separately derived stable
self-organized critical mechanism. Landau criticality alone does not supply (AT12).

Suppressing the positive overall SI conversion carried by `cal Z_T`, the minimally
admitted tensor action on an isotropic background is

\[
 S_T^{(2)}={1\over8}\int dt\,d^3x\,a^3 {\cal Z}_T
 \left[\dot h_{ij}^{TT}\dot h_{ij}^{TT}
 -{c_T^2\over a^2}(\partial_k h_{ij}^{TT})^2
 -c_T^2m_T^2 h_{ij}^{TT}h_{ij}^{TT}\right].              \tag{AT13}
\]

Here `[m_T]=m^-1`. The endpoint requires `cal Z_T>0`, `c_T^2>0`, `c_T=c`,
`m_T=0`, exactly two tensor helicities, no hidden scalar/vector ghost, and the same
source normalization in static and radiative sectors. A persistent record-rest-frame
operator can change the time and spatial coefficients independently, so `c_T=c` is
not automatic.

## 8. Can a threshold switch gravity on?

Write the SI curvature coefficient as

\[
 S_R=C_R(X)\int d^4x\sqrt{-g}\,R,
 \qquad [C_R]={\rm kg\,s^{-1}},                           \tag{AT14}
\]

with

\[
 C_R={c^3\over16\pi G_{eff}}.                           \tag{AT15}
\]

This exposes three distinct mechanisms.

### 8.1 Vanishing kinetic coefficient

If `C_R -> 0`, then `G_eff -> infinity`; the tensor kinetic term vanishes. For a
canonically normalized fluctuation `h_c ~ sqrt(C_R) h`, the Einstein cubic vertex
scales schematically as

\[
 {1\over\sqrt{C_R}}h_c(\partial h_c)^2.                  \tag{AT16}
\]

The effective strong-coupling scale collapses as `C_R -> 0`. Therefore setting a
factor `F(X)R` to zero below the threshold does **not** represent weak or absent
gravity in the same controlled metric EFT.

**EXCLUSION.** The naive `F(X_c)=0` switch-on is not an admissible controlled
gravity-off interpretation without a separate ultraviolet/transition completion.
This is not a no-go against all phase emergence; it is a no-go against reading a
zero Einstein kinetic coefficient as a benign off state.

### 8.2 Residue/coupling onset

Weak decoupling is `G_eff -> 0`, equivalently `C_R -> infinity`. A positive massless
composite pole may instead develop a residue

\[
 D_2^R(k;X)\supset {R_2(X)P^{(2)}\over k^2+i0},
 \qquad R_2(X<X_c)=0,\quad R_2(X>X_c)>0.                 \tag{AT17}
\]

This is kinematically coherent if `R_2` has positive spectral sign and universal
source coupling. In metric variables the corresponding infinite-to-finite stiffness
is nonanalytic at the off phase. It therefore needs a microscopic composite or
phase-transition derivation; it cannot be inferred from the scalar threshold alone.

### 8.3 New EFT only above the phase boundary

The metric and its gauge redundancy may be valid collective variables only for
`X>X_c`, with no continuation of (AT14) below the boundary. That avoids interpreting
`C_R=0` as the off state, but it moves the burden to deriving the new Hilbert space,
causal matching, initial state, universal source map, and finite positive residue at
the boundary.

A scalar record order multiplying `R`, such as `F(X)R`, assumes a metric and
diffeomorphism structure already exist. If `X` is dynamical, it also produces scalar-
metric mixing and generally an extra polarization unless stabilized, screened, or
decoupled. Such a branch can modulate gravity but is not by itself a spin-2 origin.

## 9. Composite emergence and the Weinberg-Witten boundary

Under the Weinberg-Witten premises, a theory with a Lorentz-covariant conserved
stress tensor cannot contain a massless spin greater than one state that carries
nonzero four-momentum in the theorem's particle sense. This constrains both composite
and elementary candidates.

**CONDITIONAL exclusion.** `AT-A2` is excluded if its microscopic theory has all of
the theorem's premises while the proposed composite graviton carries the corresponding
energy-momentum. It is not an unconditional ban on emergent gravity.

Lawful evasion classes include, without selection:

- microscopic Lorentz invariance is absent and only emerges in the infrared;
- no gauge-invariant Lorentz-covariant local microscopic stress tensor exists;
- spacetime/locality itself is emergent or the map is nonlocal;
- the gravitational energy is necessarily gauge-dependent/pseudotensorial; or
- the relevant state is not an ordinary particle within the theorem's domain.

Each evasion creates a positive obligation: show infrared Lorentz/common-cone
recovery, a complete conserved operational energy ledger, locality/causality at the
claimed scale, and the protected two-helicity pole. Merely naming an evasion earns no
gravity credit.

## 10. O0 versus O2

| question | O0 source-side | O2 metric-side |
|---|---|---|
| Where is spin-2? | In the unchanged Einstein-Hilbert metric | In the metric action/projector to be selected |
| Does the record sector create the graviton? | no; it supplies `T_R` to a pre-existing graviton | not if a nonzero EH term is retained; possibly only in a new composite/activated phase |
| Positive feature | two helicities, `c_T=c`, and GR nonlinear vertices are inherited | could in principle bind record order to the geometric response |
| Missing object | off-shell record/filter action or entropy-current completion and record-to-`T_R` map | covariant projector/action, gauge origin, healthy constraints, causal state, and record-to-O2 map |
| Primary hazard | ordinary stored energy is simply reclassified as new source | double counting, extra modes, wrong residue, zero-kinetic strong coupling |

**EXACT comparison.** O0 currently satisfies the spin-2 endpoint only by assumption:
it retains Einstein-Hilbert gravity. It is a viable record-source theory branch but
cannot prove gravity emergence. The presently written O2 architecture also contains
`S_EH`, so it modifies a pre-existing graviton unless a distinct zero-EH microscopic
phase and controlled transition are derived. A source/geometry rewrite is not new
physics unless the variational ownership and observables differ without double count.

The most economical surviving threshold interpretation is therefore `AT-A0` or the
modulation form of `AT-A3`: record accumulation selects or tunes a phase, while the
protected metric gauge sector is separately present. That economy is not canonical
selection and does not satisfy the desired record-origin proof.

The genuinely emergent alternatives `AT-A1`, `AT-A2`, and the new-EFT form of
`AT-A3` remain logically open but lack their defining construction.

## 11. Minimal SI residue and source map

Use `x^0=ct`, dimensionless `g_mu_nu`, and `d^4x` in `m^4`. Define the complete
material stress by

\[
 T_{\mu\nu}=-{2c\over\sqrt{-g}}
 {\delta S_{all}\over\delta g^{\mu\nu}},
 \qquad [T_{\mu\nu}]={\rm J\,m^{-3}}.                   \tag{AT18}
\]

Variation of (AT09) gives

\[
 G_{\mu\nu}+\Lambda g_{\mu\nu}
 ={8\pi G_*\over c^4}T_{\mu\nu}.                       \tag{AT19}
\]

The dimensions close:

```text
[c^3/G_*] = kg s^-1,
[d^4x R] = m^2,
[S] = kg m^2 s^-1 = J s,
[8 pi G_* T/c^4] = m^-2.
```

In the weak static limit,

\[
 \nabla^2\Phi=4\pi G_*\rho_{total},
 \qquad Z_\Phi={1\over4\pi G_*}={4C_R\over c^3},        \tag{AT20}
\]

with `[Z_Phi]=kg s^2 m^-3`. The same `G_*` must normalize light/clock response,
tidal multipoles, tensor radiation, and nonlinear self-coupling. It is an output
only after the pole residue and source transduction are independently calibrated;
inserting measured Newtonian `G` into (AT09) proves nothing about record origin.

## 12. Exact bounded result

The lane establishes four assumption-scoped results:

1. **Representation exclusion:** scalar and single-vector record order cannot
   linearly supply helicity-2; a symmetric tensor is necessary but not sufficient.
2. **IR unique family:** once a positive massless helicity-2 pole, Lorentzian soft
   regime, universal coupling, locality, and two-derivative consistent self-coupling
   are admitted, the leading one-metric endpoint lies in the Einstein-Hilbert plus
   cosmological-term family.
3. **Switch-on exclusion:** a vanishing Einstein kinetic coefficient is a
   strong-coupling limit, not a controlled gravity-off phase.
4. **Composite conditional no-go:** a composite graviton must explicitly evade at
   least one Weinberg-Witten premise and then recover the missing infrared properties.

What remains open is the central microscopic arrow:

```text
record accumulation / threshold
    -> exact or emergent diffeomorphism redundancy
    -> positive protected tensor pole and finite SI residue
    -> universal complete-stress coupling
    -> controlled nonlinear metric phase.
```

No parent supplies that arrow. The exact disposition is

```text
ASSUMPTION_SCOPED_UNIQUE_IR_GR_FAMILY_IDENTIFIED;
RECORD_TO_GAUGE_ORIGIN_AND_CONTROLLED_THRESHOLD_ONSET_NOT_DERIVED.
```

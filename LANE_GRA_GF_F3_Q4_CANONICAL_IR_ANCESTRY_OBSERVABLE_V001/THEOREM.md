# Canonical infrared ancestry-observable contract for G2

**Lane:** `LANE_GRA_GF_F3_Q4_CANONICAL_IR_ANCESTRY_OBSERVABLE_V001`  
**Packet revision:** `V005`, repaired after the independent V004 rejection  
**Short name:** `CIAO`  
**Date:** 2026-08-28  
**Plan gate:** prospective `C1`, serving direct-shortcut `G2`  
**Claim class:** exact observable and decision contract; no new spectral data

**Not claimed:** completion of the native source on every `G_L`, a matched
finite-size response, a massless pole, Ward-qualified helicity two,
factorization, a common physical cone, positive `G2`, a soft theorem,
gravity, or `G`.

## 1. Frozen question and ancestry invariant

Phase A supplies one fixed record-parent tensor source, exact finite-family
geometry, local coefficient ancestry, and finite nonzero sampled response.
It does not supply the thermodynamic pole.  The direct shortcut requires the
infrared object to remain this same source descendant rather than a later
inserted tensor field.

For `G_L`, `L=5*2^r`, define the complete ground projector `P_(0,L)` and
candidate projector or controlled cluster `P_C,L(k)`.  Before response is
inspected, the ground query must be fixed as either (i) a rank-one `P_(0,L)`,
or (ii) the basis-invariant density `rho_(0,L)=P_(0,L)/rank(P_(0,L))` in one
exactly typed ground/flux/boundary sector.  In the second case all correlators
below use the same `rho_0`; choosing an individual vector in the degenerate
ground space is forbidden.  The load-bearing ground-to-pole map to be tested
is

\[
 P_{C,L}(k)Q^{TT\dagger}_{L,\rm vol}(k)P_{0,L}
 =Z^{{\rm can}}_{Q,L}(k)h^{TT\dagger}_{L,\rm can}(k)P_{0,L}
   +\mathcal R_L(k)P_{0,L},                        \tag{GF01}
\]

where the field has a positive unit time-kinetic convention, `Z_Q^can` is
measured rather than normalized away, and the remainder is analytic or
subleading in the same declared limit.  For a mixed `rho_0`, (GF01) is read as
an operator map on `ran(P_0)`, while residue and ancestry Grams use `rho_0`
and are not formed from an arbitrary purification.  This packet fixes every
object in (GF01) before its response is calculated.

### Theorem `CIAO-0` -- no-substitution design rule

A G2 result is admissible only if its Hamiltonian, source derivative,
physical-volume conversion, clock, scale, momentum class, TT projector,
spectral cluster, and canonical source vertex match this frozen contract.
Changing any of those after examining a spectrum creates a different test
and cannot prove Phase-A ancestry.

## 2. Affine-volume source normalization and physical binding

GC gives the affine primitive-cell volume

\[
 v_3={16\over3\sqrt3}a_*^3,
 \qquad V_L=L^3v_3.                                \tag{GF02}
\]

Until one independently custodied FD cell/scale certificate binds `a_*` to a
physical rod and the primitive cell to a physical measure, `v_3` and `V_L`
are **affine** volumes.  After that one global binding, the same quantities
are physical volumes.  No tensor, EM, or matter sector may receive its own
length fit.
The operators `q_x` are integrated-cell operators; (GF05) would not be the
normalization of a pointwise density with a separately chosen cell weight.

Let `q_x^{ij}([n])` be the **complete** integrated-cell source derivative,
including pair, diagonal-fold, hopping-numerator, and ring supports:

\[
 q_x^{ij}([n])=-2{\partial H_L[j]\over
                  \partial j_{ij,x}(-[n])}\bigg|_{j=0}.      \tag{GF03}
\]

This definition is prospective; the all-`G_L` native ledger is still a
required calculation.  With `N_L=L^3`, GC's count-normalized mode is

\[
 Q^{ij}_{L,\rm cnt}([n])={1\over\sqrt{2N_L}}
 \sum_x e^{2\pi i\,n\cdot x/L}q_x^{ij}([n]).       \tag{GF04}
\]

The unique conversion to an integrated-cell, affine-volume-normalized
Fourier mode is therefore

\[
 \boxed{
 Q^{ij}_{L,\rm aff}([n])
 =\sqrt{2\over v_3}\,Q^{ij}_{L,\rm cnt}([n])
 ={1\over\sqrt{V_L}}
 \sum_x e^{2\pi i\,n\cdot x/L}q_x^{ij}([n]).}     \tag{GF05}
\]

The conjugate source must transform oppositely,

\[
 j_{\rm aff}=\sqrt{v_3/2}\,j_{\rm cnt},
 \qquad j_{\rm aff}:Q_{\rm aff}=j_{\rm cnt}:Q_{\rm cnt}.   \tag{GF06}
\]

After the independent binding certificate is supplied, the same paired
objects are denoted `Q_vol` and `j_vol`; no further multiplicative conversion
is allowed.  Thus (GF05) is fixed by the source pairing and cell measure, not
by the observed residue.  It must not be confused with a wavefunction-residue
fit.

### Theorem `CIAO-1` -- exact count-to-affine-volume conversion

Given GC's cell volume and FY's count convention, (GF05)--(GF06) are forced.
Any other multiplicative conversion changes the declared affine source or
its conjugate query.  Calling it physical additionally requires the single
independently custodied FD binding certificate.  A later momentum-dependent
deconvolution is admissible only after the corresponding derivative map has
been derived from the same parent and frozen in a successor contract before
response inspection.

## 3. One parent clock, scale, and momentum family

For `U_d>0` and `h!=0`, use the inherited ring scale as the sole clock,

\[
 x={h\over U_d},\qquad
 E_*=J_6={63\over8}U_dx^6>0,
 \qquad t_*={\hbar\over J_6}.                      \tag{GF07}
\]

Every tensor, EM, and admitted matter-probe energy is an eigenvalue or pole
of the same parent `H_L` and is reported as `delta=Delta/J_6`.  Every momentum
uses the same `a_*` and is reported as `p=a_*k`.  No sectorwise time or length
rescaling is permitted.

For every class `[n] in (Z/LZ)^3`, let `hat n` be the lexicographically first
global minimizer in its class of

\[
 q(n)=4\sum_i n_i^2-\left(\sum_i n_i\right)^2.
\]

Then

\[
 k_{L,[n]}={2\pi\over L}A^{-\mathsf T}\widehat n,
 \qquad
 |k_{L,[n]}|^2={3\pi^2\over4L^2a_*^2}q(\widehat n).          \tag{GF08}
\]

The raw source registry is defined on the complete set of `L^3` characters,
including conjugate class `[-n]`.  TT projection, pole clusters, residues,
and ancestry are defined only on the `L^3-1` nonzero characters.  At
`[n]=[0]`, the packet stores the unprojected homogeneous symmetric-tensor
source and exact scalar/point-group sector data; it does not define
`Pi_TT(0)`.  Every `k->0` TT limit is a directional limit along a declared
nonzero sequence.  Conjugacy is a character-level statement; at a
Brillouin-zone boundary the lexicographically selected minimum representative
of `[-n]` need not be the literal negative of the selected representative of
`[n]`, although their reciprocal norms agree.  The mandatory low-momentum
screen freezes the three
inequivalent rays and their conjugates:

\[
 n^{(1)}=(1,0,0),\quad q=3;\qquad
 n^{(2)}=(1,1,0),\quad q=4;\qquad
 n^{(3)}=(1,-1,0),\quad q=8.                      \tag{GF09}
\]

The first gives `|k_min|=3pi/(2La_*)`.  A base character pulls back as
`[n] -> [2n]` on `G_(2L)` and preserves its physical `k`; each cover's new
odd lowest modes are new response points, not inherited spectra.

For `k!=0`, apply GC's exact six-component symmetric-tensor projector

\[
 \Pi^{TT}_{ij,kl}(k)={1\over2}
 (P_{ik}P_{jl}+P_{il}P_{jk}-P_{ij}P_{kl}),
 \qquad P_{ij}=\delta_{ij}-{k_ik_j\over |k|^2},   \tag{GF10}
\]

directly to the six-component source.  Pass/fail quantities use the image
projector and its eigenvalues, not a freely rotated plus/cross basis.  A
deterministic projected Gram--Schmidt basis may be stored for arrays, but it
has no physical claim role.

### Theorem `CIAO-2` -- complete kinematic registry

Equation (GF08) uniquely registers every finite momentum and its conjugate.
Equation (GF10) assigns a two-dimensional TT quotient only to nonzero
characters.  The homogeneous character is separately typed as above.  The
three rays test distinct reciprocal norms.  This registry is exact kinematics
only and cannot supply a gap, pole, representation, or cone.

## 4. Basis-independent pole and coalescing-cluster observable

At fixed momentum and exact conserved quantum numbers, let `P_E` denote the
**complete** projector onto an energy eigenspace.  Prospectively select a
candidate cluster `C_L(k)` and define

\[
 P_{C,L}(k)=\sum_{E\in C_L(k)}P_E.                 \tag{GF11}
\]

The selection rule may use only frozen sector labels and an energy window or
branch-continuation rule declared without reference to individual vectors.
Exact degeneracy is retained whole.  If two finite-volume polarizations are
split, both are retained and must obey

\[
 {\operatorname{diam}\{\Delta_s:s\in C_L\}\over |k|}\to0.  \tag{GF12}
\]

An accidental degeneracy with extra states must be split by exact conserved
quantum numbers or included.  A positive G2 result requires exactly two
physical states in the typed pole sector; source rank two alone cannot hide
an extra propagating mode.

Let `rho_(0,L)` denote the frozen rank-one ground projector or the normalized
complete-ground density fixed in Section 1.  For the two-dimensional TT image
at `k!=0`, define the positive raw-source spectral matrices

\[
 R^{{\rm raw},(s)}_{AB,L}(k)=
 {\rm Tr}\!\left[\rho_{0,L}Q^{TT}_{A,L,\rm vol}(k)P_s
 Q^{TT\dagger}_{B,L,\rm vol}(k)\right]_c,         \tag{GF13}
\]

\[
 R^{\rm raw}_{C,L}=\sum_{s\in C_L}R^{{\rm raw},(s)}_L,
 \qquad
 \mathcal N^{\rm raw}_{C,L}
 =\sum_{s\in C_L}2\Delta_{s,L}R^{{\rm raw},(s)}_L.         \tag{GF14}
\]

The factor `2Delta` belongs specifically to GC's frozen even-Euclidean
quadratic-denominator convention.  The raw numerator may vanish with a
derived momentum power; a constant is not assumed.

On the TT image with its inherited Frobenius inner product, the raw source
vertex is the unique positive square root

\[
 Z^{\rm raw}_{Q,L}=(\mathcal N^{\rm raw}_{C,L})^{1/2}.      \tag{GF15}
\]

Its two eigenvalues (equivalently, its singular values) are invariant under
orthonormal TT reporting rotations and must be reported.  The matrix entries
of a square root are not claimed invariant under an arbitrary nonorthogonal
coordinate change.  Equation (GF15) does not authorize replacing `Q_vol` by
`N^{-1/2}Q_vol`, which would force unit residue and erase the test.

The canonical pole field has positive unit time-kinetic term.  A derived
canonical source vertex may be claimed only if a same-parent source map
`D_L(k)` has been derived and frozen before response inspection.  It gives

\[
 \mathcal N^{\rm can}_{C,L}
 =D_L\mathcal N^{\rm raw}_{C,L}D_L^\dagger,
 \qquad
 Z^{\rm can}_{Q,L}=(\mathcal N^{\rm can}_{C,L})^{1/2}.       \tag{GF15a}
\]

`D=I` is the raw reporting convention, not by itself a proof that a
derivative source is a nonderivative canonical field.  Removing a momentum
power requires the separately derived and prospectively frozen map.

The amplitude decision is three-valued:

1. **PASS:** an applicable frozen `D_L(k)` exists and both singular values of
   `Z_Q^can` have finite strictly positive limiting lower bounds on every
   declared ray;
2. **FAIL:** the complete raw source is identically or persistently rank
   deficient, or an applicable frozen map produces canonical rank loss, a
   vanishing canonical lower singular value, or an unbounded canonical upper
   singular value; and
3. **INDETERMINATE:** the raw numerator remains full rank but no applicable
   same-parent canonical map has yet been derived, whether the raw vertex is
   nonvanishing or carries a controlled momentum power.

INDETERMINATE is neither positive G2 nor a refutation of a derivative-source
pole.  Exact positive spectral matrices cannot have a negative eigenvalue;
a reported negative value beyond certified numerical error is instead an
input/covariance failure.

This three-valued output classifies **canonical amplitude only**.  Normalized
ancestry is the independent gate in Section 5.  In particular, a family may
have canonical-amplitude PASS and ancestry FAIL; its overall G2 result is
then FAIL.  The two labels must not be collapsed into one classifier.

### Theorem `CIAO-3` -- degeneracy-safe raw and canonical source vertices

The spectral sums (GF13)--(GF14) are invariant under arbitrary basis rotations
within exact excited eigenspaces, under reordering of a split cluster, and
under arbitrary basis rotations in a degenerate ground space when the frozen
normalized complete-ground density is used.
Under an orthonormal TT reporting rotation their coordinate matrices,
including the positive square root (GF15), transform by conjugation, so their
ranks and spectra are invariant.  They distinguish physical pole-space rank,
source-visible rank, sign, and amplitude.  They do not manufacture nonzero
amplitude.

## 5. Normalized Phase-A ancestry

Using the same frozen `rho_0` as in (GF13), let

\[
 S_{AB,L}={\rm Tr}[\rho_{0,L}\delta Q^{TT}_{A,L,\rm vol}
                     \delta Q^{TT\dagger}_{B,L,\rm vol}],
 \qquad
 R_{AB,L}={\rm Tr}[\rho_{0,L}\delta Q^{TT}_{A,L,\rm vol}P_{C,L}
                     \delta Q^{TT\dagger}_{B,L,\rm vol}].  \tag{GF16}
\]

Here `delta Q=Q-Tr(rho_0 Q)` removes the disconnected ground component (and
is automatic at a nonzero translation character when `rho_0` is translation
invariant).  For rank-one `rho_0=|0><0|`, (GF16) reduces to the source-vector
Gram used in V001.

When `S_L` is positive definite, define the normalized ancestry matrix

\[
 \Omega_L=S_L^{-1/2}R_LS_L^{-1/2}.                \tag{GF17}
\]

Equivalently, its eigenvalues are the generalized eigenvalues of
`R_L v=eta S_L v`, avoiding any basis choice.  Since `0<=P_C<=I`,
`0<=Omega<=I`.  The same-parent two-polarization ancestry gate is

\[
 \boxed{\liminf_{r\to\infty}\lambda_{\min}(\Omega_{L_r}(k_r))>0}              \tag{GF18}
\]

on each declared IR ray, after only the source renormalization already
frozen in this contract.  A singular `S` or `eta_min->0` fails this complete
Phase-A source channel.  It does not exclude a wholly source-invisible tensor
sector, which could not establish this shortcut's ancestry in any event.

### Theorem `CIAO-4` -- scale-independent ancestry discriminator

The generalized spectrum of `(R,S)` is invariant under every invertible
change of TT reporting basis and every common nonzero source-unit conversion.
Thus (GF18) tests projection of the Phase-A source state into the candidate
infrared sector rather than an arbitrary amplitude convention.

## 6. Helicity-two representation, spectral, and common-cone discriminator

Spatial TT rank two is a kinematic source quotient, not a helicity theorem.
A positive G2 result must additionally construct, on the same nonzero IR
sequence and pole bundle, a unitary massless Poincare representation in the
common limiting cone.  Its translations must have the measured pole
four-momentum; rotations must obey

\[
 U(R)P_C(k)U(R)^\dagger=P_C(Rk),                  \tag{GF18a}
\]

with covariant pole residues and amputated vertices.  On the complexified
two-state fiber, rotations by `theta` about `k` must have eigenvalues

\[
 \boxed{\{e^{+2i\theta},e^{-2i\theta}\}},          \tag{GF18b}
\]

equivalently little-group generator spectrum `{+2,-2}`.  Helicity zero,
`+/-1`, a reducible scalar doublet, or an arbitrary momentum-dependent
identification of two scalar poles with a TT reporting basis fails G2.  The
finite lattice need not have exact continuous rotations, but the continuum
pole projectors, residues, and multipoint amplitudes must converge uniformly
to (GF18a)--(GF18b); equality of dispersive slopes is insufficient.

This representation gate does not assume the G3 Ward result.  Longitudinal
or gauge-null replacement in physical amplitudes and soft universality remain
owned by G3.

For every branch `alpha` in the typed two-state cluster and every frozen ray,
a positive spectral result requires

\[
 {\Delta_{\alpha,L}(k)\over\hbar|k|}\to c_T>0,
 \qquad {\Delta_{\alpha,L}-\Delta_{\beta,L}\over\hbar|k|}\to0,               \tag{GF19}
\]

with one isotropic `c_T`, `rank P_C=rank R_C=2`, the helicity representation
(GF18a)--(GF18b), nonvanishing normalized ancestry, and no negative-norm or
extra longitudinal long-range mode.  The cluster must lie below the multiparticle threshold by
a controlled margin, or remain an exact symmetry-protected delta function at
threshold.  It must also survive a prospectively frozen neighborhood of
admitted higher inherited operators with volume-uniform bounds.

All speeds are compared using the same parent units.  For
`X in {TT,EM,matter probe}`, define

\[
 \widehat c_X^2=\lim_{k\to0}
 { [E_X(k)^2-E_X(0)^2]/J_6^2\over |k|^2a_*^2}.    \tag{GF20}
\]

For a gapless branch this reduces to
`(Delta_X/J6)^2/(|k|a_*)^2`.  A common candidate Lorentz cone requires one
positive isotropic value of (GF20) for every admitted sector and ray, with no
sectorwise clock or length fit.  Equality of slopes is necessary but is not,
by itself, a proof of the Poincare representation required above.

The exact failures are:

- a nonzero limiting tensor mass, `z!=1`, or directional limiting velocity;
- continuum-only weight, dissolving residue, or loss of an isolated or
  protected cluster;
- physical/projected rank other than two, a negative spectral sign, or an
  extra long-range mode;
- failure of continuum Poincare covariance or little-group helicities
  `{+2,-2}`;
- `eta_min->0`, canonical amplitude FAIL, or canonical amplitude
  INDETERMINATE at the attempted positive-G2 promotion;
- incompatible tensor/EM/matter values of (GF20); or
- loss of the candidate under admitted inherited operators.

A controlled raw numerator power without a derived canonical map is exactly
the INDETERMINATE outcome of Section 4.

## 7. Factorization interface required by G3

G2 must deliver one of two interfaces for the **same** projector and source.

1. **Asymptotic/LSZ route.**  The two isolated stable states admit a
   volume-uniform infinite-volume limit.  Mixed correlators have the same
   `Z_Q^can`, and pole amputation gives finite source-independent vertices.
2. **Retarded/1PI route.**  In a common soft neighborhood,

   \[
   G^R_{Q\mathcal O\cdots}
   =Z_Q^{\rm can}D_C^{-1}\Gamma_{h\mathcal O\cdots}
    +G^R_{\rm analytic},                           \tag{GF21}
   \]

   and `D_C (Z_Q^can)^+ G^R` has a finite,
   source-independent, volume-uniform limit.  Here `+` is the Moore--Penrose
   inverse on the rank-two pole image.  Before exact coalescence, `D_C` is the
   diagonal branch-denominator matrix (or the equivalent inverse kernel on
   the full cluster projector), never a scalar average of split gaps.

At least one route must pass.  A branch cut that replaces the pole,
nonfactorizing residues, incompatible amputated limits from two source
probes, or the absence of a uniform soft neighborhood fails factorization.
Analytic contacts may alter the regular term but cannot replace a missing
pole residue.

This is factorization **readiness** for the helicity-typed pole, not the G3
longitudinal-decoupling or soft-universality theorem.

## 8. Frozen pass/fail packet and present disposition

An executable G2 packet must contain:

1. the complete native `q_x^{ij}([n])` ledger on every scored `G_L`;
2. a rank-one ground projector or the fixed normalized complete-ground
   density, plus exact flux/boundary sector and shared `H_L,J6,a_*` metadata;
3. full TT correlators or spectral projectors and covariance;
4. every raw momentum class, the three nonzero rays and conjugates, and
   separately typed homogeneous data with no `Pi_TT(0)`;
5. cluster gaps, full projectors, residue matrices, continuum thresholds,
   and ancestry Grams;
6. pole-bundle rotation/boost covariance and little-group data sufficient to
   prove helicities `{+2,-2}` and reject a scalar doublet;
7. matched EM and matter-probe dispersions under the same clock;
8. the higher-inherited-operator deformation ledger;
9. multipoint data for at least one factorization route; and
10. an independently audited and sealed FD cell/scale certificate binding the
    one `a_*`, `v_3`, and cell map physically for every admitted sector.

### Theorem `CIAO-5` -- prospective C1 contract

Equations (GF02)--(GF21), including the representation and three-valued
amplitude gates, form one source-bound and basis-independent decision
procedure for G2.  It cannot be rescued by a
post-result normalization, a fitted interaction, an inserted tensor field,
or sector-specific clocks.

The present V005 author packet preserves the V002--V003 repairs and closes
the missing upper-bound branch found by the independent V004 re-audit.  It
awaits a new independent seal and has
not assembled the all-`G_L` native ledger, produced matched-family pole data,
constructed the Poincare representation, derived a canonical source map, or
supplied the physical FD binding certificate.  Consequently it proves no
positive G2 statement.

## 9. Disposition

`V005_REPAIRED_AFTER_V004_REJECTION__AMPLITUDE_CLASSIFIER_DISJOINT_TOTAL_AND_UPPER_BOUNDED__AMPLITUDE_AND_ANCESTRY_SEPARATE__PROSPECTIVE_C1_G2_OBSERVABLE_CONTRACT__Q_AFF_EQUALS_SQRT_2_OVER_V3_TIMES_Q_CNT__PHYSICAL_VOLUME_REQUIRES_INDEPENDENT_FD_BINDING__TT_ONLY_AT_NONZERO_K__K0_HOMOGENEOUS_DATA_SEPARATELY_TYPED__GROUND_QUERY_BASIS_SAFE__RAW_CANONICAL_RESIDUE_PASS_FAIL_INDETERMINATE__POINCARE_LITTLE_GROUP_HELICITIES_PLUS_MINUS_TWO_REQUIRED__TWO_SCALAR_COUNTEREXAMPLE_REJECTED__COMMON_CONE_NOT_USED_AS_BOOST_PROOF__FACTORIZATION_INTERFACE_EXACTLY_TYPED__DESIGN_ONLY__NO_NUMERICAL_POLE_NO_POSITIVE_G2_NO_G3_NO_GRAVITY_NO_G`

# SPAG--EIR consistency and empirical path audit

**Audit ID:** `GRAVITY-SPAG-EIR-PATH-V001`

**Date:** 2026-08-27

**Scope:** one physics consistency check plus a bounded screen of primary
experimental sources.  This memo does not amend RGRL, RIEC, GFT, or the frozen
SPAG protocol.

**Disposition:**
`REAL_SIGNAL_CLASSIFICATION_DEFECT_FOUND__NO_PUBLIC_DATASET_IDENTIFIES_THE_SPAG_PRIMARY_ESTIMAND__NIST_BIPM_IS_THE_CLOSEST_ADAPTABLE_FORCE_PLATFORM__A_NEW_EIGHT_CELL_DISCOVERY_RUN_CAN_BE_EXECUTED_AFTER_THE_SOURCE_CLASSIFICATION_IS_FROZEN`

## 1. Result

No public acquisition found here can estimate the registered SPAG coefficient

\[
 \beta_{TM}={1\over8}\sum_{t,d,m\in\{-1,+1\}}tm\,q_{t,d,m}.
 \tag{PA01}
\]

The reason is identification, not merely sensitivity.  None of the public
gravity datasets contains independently randomized target-lineage and
dummy-lineage redistributions.  Page--Geilker contains the nearest historical
branch variable, but its decision and mass placement are deterministically
locked.  Existing data can calibrate ordinary gravity, noise, transfer, and
systematics; they cannot be relabelled after the fact as KEEP/BREAK data.

The screen also exposed a real internal physics constraint.  A nonzero SPAG
response cannot be attributed to an RGRL constitutive coordinate **alone** if
the two arms have identical complete stress, identical remainder, and identical
gravitational initial/boundary data while both satisfy the same leading
Einstein equation.  Such a signal must enter the physical equations as one of:

1. a lineage-dependent contribution to the complete source;
2. a lineage-dependent remainder or extra leading field;
3. different gravitational initial, incoming, or boundary data; or
4. an unmatched ordinary collateral channel.

This is a classification defect in the present SPAG interpretation, not a
reason to abandon the experiment.  The experiment becomes sharper when it is
registered as a test for an **ancestry-conditioned departure from the matched
EIR solution**, with the possible source buckets fixed before data exposure.

## 2. Exact EIR consistency lemma

Let `K` and `B` denote KEEP and BREAK.  The closed GFT endpoint equation is

\[
 G_{\mu\nu}[g_a]+\Lambda g_{a\mu\nu}
 ={8\pi G_{\rm eff}\over c^4}T^{\rm complete}_{a\mu\nu}
 +\Delta^{\rm rem}_{a\mu\nu},
 \qquad a\in\{K,B\}.                              \tag{PA02}
\]

In the weak-field SPAG domain, subtraction and linearization about the common
background give

\[
 {\cal E}^{(1)}\delta g
 ={8\pi G_{\rm eff}\over c^4}\delta T^{\rm complete}
 +\delta\Delta^{\rm rem},                         \tag{PA03}
\]

where `delta` means KEEP minus BREAK.  The retarded solution projected onto the
registered response coordinate is therefore

\[
 \delta q_+
 =\mathcal P_+\mathcal G_R
 \left[{8\pi G_{\rm eff}\over c^4}\delta T^{\rm complete}
       +\delta\Delta^{\rm rem}\right]
 +\delta q_{+,\rm hom}+\delta q_{+,\rm coll}.     \tag{PA04}
\]

All common terms proportional to `delta g`--including ordinary induced matter
response and any admitted linearized remainder response--are understood to be
included in the full linearized operator.  The right-side deltas in
(PA03)--(PA04)
denote independent arm differences after that common response is owned.

Here `G_R` is the retarded leading-Einstein Green operator,
`delta q_hom` owns different initial/incoming/boundary data, and
`delta q_coll` owns ordinary experimental mismatch.

Consequently, under

\[
 \delta T^{\rm complete}=0,\qquad
 \delta\Delta^{\rm rem}=0,\qquad
 \delta q_{+,\rm hom}=0,\qquad
 \delta q_{+,\rm coll}=0,                         \tag{PA05}
\]

well-posedness gives

\[
 \boxed{\delta g=0\ \text{up to gauge},\qquad \beta_{TM}=0.} \tag{PA06}
\]

The same conclusion holds nonlinearly on any locally unique Cauchy-development
domain with identical complete source and Cauchy/boundary data.  It does not
require global uniqueness of every Einstein solution.

RGRL-B's `J_A` fields can provide a right-inverse coordinate system for the
metric variation tangent without supplying a new term on the right side of
(PA03).  Thus

\[
 \delta g=(D_Jg)\,\delta J
\]

is a constitutive/kinematic tangent statement.  It is not, by itself, an
on-shell forcing law.  If a KEEP/BREAK intervention moves the stationary
solution, the term that moves it must be owned in (PA04).

## 3. Exact signal classification

| Observed source of `beta_TM` | EIR/GFT classification | What a positive result would mean |
|---|---|---|
| `delta T_lineage != 0`, where `T_lineage` is the metric variation of an independently frozen record/lineage term in `Gamma_<^retained` | part of `T_complete`, counted exactly once | leading-Einstein-compatible response to a newly identified physical source term; not a source-free history effect |
| lineage changes `C_R`, `G_eff`, or another unsuppressed kinetic coefficient | extra leading constitutive field unless derived as one constant on the connected phase | departure from EIR-4's one-metric constant-coefficient endpoint; it narrows or falsifies that endpoint packet |
| `delta Delta_rem != 0` | owned `Gamma_rem` correction/extra bounded mode | RGRL-related only after an independent type-join; not the leading Einstein response |
| `delta q_hom != 0` | changed initial, incoming, asymptotic, or apparatus boundary data | history-conditioned gravitational boundary response; not a local source coefficient |
| `delta q_coll != 0` | ordinary force, stress, EM, heat, support, controller, termination, geometry, or analysis leakage | no lineage result |
| nonzero residual with none of the above independently identified | unclassified ancestry-correlated anomaly | discovery result only; it cannot yet be called RGRL-C or Gravity Formation confirmation |

A lineage-dependent `G_eff(L)` acting on the ordinary mass geometry is not a
loophole.  It makes the mass-odd SPAG interaction nonzero, but it also changes
the endpoint action coefficient.  Unless a new dynamical coefficient field and
its Ward/stationarity equations are supplied, that is outside the fixed
constant-coefficient EIR theorem.

The adopted RGRL-C language matches "complete stress" while requiring a
nonzero matched lineage-to-metric column.  That can be read consistently only
as an off-shell tangent/constitutive statement, or with one of the source
buckets above left unmatched by definition.  It cannot simultaneously mean a
new on-shell force at fixed (PA05).

## 4. Consequence for the registered SPAG claim

SPAG presently matches **ordinary** measured stress and collateral in several
operational passages, but its theory basis inherits RGRL-C's stronger
"complete stress matched" phrase.  Before construction, the program should
freeze one of two scientifically clean interpretations:

### Interpretation A -- EIR-consistency test

Treat the entire lineage redistribution as a physical-state intervention and
model every admitted lineage term in `Gamma_<^retained`.  SPAG then tests a
frozen `delta T_lineage` and its predicted retarded metric response.  A pass is
a `LINEAGE_DEPENDENT_COMPLETE_SOURCE_COLUMN_PASS`.

This is the only route by which a persistent signal can remain inside the
leading Einstein equation without changing its coefficient or boundary data.
It requires a quantitative source functional before confirmatory Run B; a
post-hoc fitted residual is insufficient.

### Interpretation B -- discovery of physics outside the leading endpoint

Match every term admitted to `T_complete` and search for a residual.  A
positive result is then classified as remainder/extra-mode, boundary-history,
ordinary collateral, or unclassified anomaly.  It is a valuable falsification
test of the joint `RGRL + EIR` package, but it cannot be advertised as direct
confirmation of the already assumed leading Einstein composition.

Under either interpretation, a null is scientifically useful as a bound on the
registered ancestry-conditioned column.  Without a derived coefficient floor,
it is not a universal refutation of RGRL.

## 5. Public experiment and data screen

### 5.1 Page--Geilker: closest historical logic, exact rank failure

Page and Geilker used a radioactive-decay/Geiger decision to choose the
macroscopic source-mass placement and measured the corresponding torsion
response.  It is the strongest existing branch-following analogue.  But its
decision variable and mass placement obey `X=M` on every run.  The design has
no off-diagonal support at fixed mass geometry, so a decision/record
coefficient and the ordinary mass coefficient are exactly inseparable.  The
published result therefore cannot estimate (PA01).

Primary source: D. N. Page and C. D. Geilker,
["Indirect Evidence for Quantum Gravity," *Physical Review Letters* 47, 979
(1981)](https://doi.org/10.1103/PhysRevLett.47.979).

### 5.2 NIST/BIPM 2026: closest adaptable force apparatus

The replicated BIPM torsion balance at NIST is the strongest directly
adaptable platform found in this screen.  It already has:

- four source masses on a stepper-motor carousel;
- two sign-reversing source angles separated by `37.67 degrees`;
- grounded-vacuum electrostatic shielding;
- independent free-deflection and electrostatic-servo readout methods;
- complete dimensional metrology and extensive environmental/systematic
  modelling; and
- a large calibrated peak-to-peak signal: `31.1979 nN m` for four copper
  source masses of `11.19168 kg` each.

Its published configuration-averaged Type-A torque uncertainties are
`0.0002--0.0006 nN m`, or `2e-13--6e-13 N m`.  Those numbers demonstrate an
instrument scale, not a prospective SPAG detection limit; the complete
ancestry-run covariance and new route systematics would have to be measured.

The acquisition contains no target/dummy lineage factors, authenticated
KEEP/BREAK routing, or scored sham network.  The paper exposes fitted torque
summaries and uncertainty budgets but no public event-level lineage dataset.
It therefore cannot test (PA01) retrospectively.

Primary/official sources: S. Schlamminger et al.,
["Redetermination of the gravitational constant with the BIPM torsion balance
at NIST," *Metrologia* 63, 025012
(2026)](https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=961075), and the
[NIST publication record](https://www.nist.gov/publications/redetermination-gravitational-constant-bipm-torsion-balance-nist).

### 5.3 Westphal: closest small-source torsion endpoint

Westphal et al. measured the field of a sub-100-mg gold source with a shielded
torsion pendulum and spatial source modulation.  It is an excellent small
source and electrostatic-systematics precedent, but the source motion is not a
closed record-lineage factorial.  No `L_T` or `L_D` exists in the acquisition.

Primary source: T. Westphal et al.,
["Measurement of gravitational coupling between millimetre-sized masses,"
*Nature* 591, 225--228
(2021)](https://www.nature.com/articles/s41586-021-03250-7).

### 5.4 Fuchs: best public force-response traces, no lineage factor

Fuchs et al. measured a kilogram-source gravitational signal with a `0.43 mg`
levitated detector near `26.7 Hz` and deposited lock-in time traces.  The
deposit is useful for detector noise, lock-in estimators, and source-position
response development.  The source is a continuously driven wheel, not a
writer-off retained mass record, and neither ancestry factor is present.
Moreover, the deposit does not by itself contain the full raw source-motion,
SI-transfer, calibration, and covariance chain needed to turn a post-hoc
history label into a physical source effect.

Primary/data sources: T. M. Fuchs et al.,
["Measuring gravity with milligram levitated masses," *Science Advances* 10,
eadk2949 (2024)](https://pmc.ncbi.nlm.nih.gov/articles/PMC10889343/), and
[Zenodo 10300430](https://zenodo.org/records/10300430).

### 5.5 Panda: strongest optional atom-probe component, not joined data

Panda et al. measured a miniature source's attraction with a lattice atom
interferometer,

\[
 a_{\rm mass}=33.3\pm5.6_{\rm stat}\pm2.7_{\rm syst}\ {m nm\,s^{-2}}.
\]

This is a strong precedent for a later compositionally different probe.  Its
public deposit contains presented-data products while the analysis code is
request-only; more importantly, the acquisition has no SPAG ancestry factors
and no common physical parent with any torsion experiment.  It cannot be
software-joined to create a common-freefall SPAG result.

Primary/data sources: C. D. Panda et al.,
["Measuring gravitational attraction with a lattice atom interferometer,"
*Nature* 631, 515--520
(2024)](https://www.nature.com/articles/s41586-024-07561-3), and
[Zenodo 10995225](https://doi.org/10.5281/zenodo.10995225).

### 5.6 Yan: low-frequency response precedent, no external source contrast

Yan et al. operated an optomechanical torsion pendulum with a `0.6 mHz`
eigenfrequency for three months and reported `0.3 microrad/sqrt(Hz)` sensitivity
near `2.5 mHz`.  It constrains apparatus design for low-frequency anomalous
response searches.  It does not contain a moved external source, a retained
mass record, or either SPAG lineage factor.

Primary source: T. Yan et al.,
["First result for testing semiclassical gravity effect with a torsion
balance," *Physical Review D* 111, 082007
(2025)](https://doi.org/10.1103/PhysRevD.111.082007).

### 5.7 Preparation-history theory is not an existing outcome

Fedida and Kent show that distinguishing different preparations with the same
reduced density operator is nonstandard physical content: ordinary mixture
equivalence forbids it, while some nonlinear or semiclassical gravity models
violate it.  This is closely analogous to why a same-final-state SPAG response
must appear somewhere in the physical state/action rather than remain a bare
historical label.  Their result is theoretical and supplies no SPAG outcome.

Primary source: S. Fedida and A. Kent,
["Mixture equivalence principles and postquantum theories of gravity,"
*Physical Review D* 111
(2025)](https://doi.org/10.1103/pttr-6kj7).

## 6. Smallest executable experiment

The smallest **scientifically interpretable** new acquisition is a
single-probe, discovery/bound retrofit of a high-signal motorized Cavendish
platform.  The NIST/BIPM architecture is the best concrete template; a new
Page--Geilker-scale apparatus could implement the same logic with lower
metrological maturity.

The minimum acquisition retains all eight registered cells:

\[
 M\times L_T\times L_D
 \in\{-1,+1\}^3.                                  \tag{PA07}
\]

It adds only the missing causal controls to an existing source-position
experiment:

1. Randomize the common mass setting `nu=M` before either original record is
   formed.
2. Prepare one authenticated descendant command buffer and one ancestry-
   disjoint buffer with identical registered motion waveforms in each target
   and dummy family.
3. Route the selected target waveform to the source carousel and the selected
   dummy waveform to a mechanically matched dummy; route both complements to
   measured isolated terminations.
4. In every arm execute the same reset, carousel trajectory, passive final
   hold, route cut, settling interval, and terminal source metrology.
5. Populate balanced eight-cell superblocks and separately powered sham-route
   cells in which controller/selector records change while ancestry allocation
   does not.
6. Read one torsion coordinate first.  The NIST free-deflection and servo
   methods are valuable internal cross-checks but use the same pendulum and do
   not constitute SPAG's independent common-freefall probe.
7. Freeze at least two hold-time strata and both physical source orientations.
   A decaying hold-time effect is classified first as transient/homogeneous or
   collateral; a persistent mass-odd effect proceeds to the source-bucket
   analysis in section 3.
8. Run in discovery/bound mode.  No RGRL-derived coefficient floor presently
   exists, so a confirmatory null can bound only this apparatus column.

Before the first scored gravity exposure, the protocol must also freeze which
of `delta T_lineage`, `delta Delta_rem`, or gravitational boundary data is the
candidate theory channel.  If none is predicted, the maximum honest positive
verdict is

`REPRODUCED_ANCESTRY_CORRELATED_GRAVITY_RESIDUAL__SOURCE_BUCKET_UNRESOLVED`.

Only after a positive torsion result and an acquisition-disjoint replication
should a Panda-class atom channel be added for the common-freefall test.

## 7. Program consequence

The immediate physics path is therefore:

\[
 \boxed{
 \text{derive/freeze the lineage source bucket}
 \ \longrightarrow\ 
 \text{retrofit the eight-cell torsion acquisition}
 \ \longrightarrow\ 
 \text{held-out replication}
 \ \longrightarrow\ 
 \text{optional atom common-response test}.}
\]

Public data already establish that the detector and source-motion components
exist.  They do not establish or bound the registered SPAG ancestry
coefficient.  The only theory repair needed before execution is physical, not
mechanical: specify where a nonzero lineage-conditioned metric solution lives
in the complete EIR equation.

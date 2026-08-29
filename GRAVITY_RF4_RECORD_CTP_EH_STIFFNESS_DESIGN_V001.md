# RF4 record-conditioned CTP / Einstein--Hilbert stiffness design

**Design ID:** `GFT-RF4-RCEHS-V001`

**Date:** 2026-08-29

**Status:** `PROOF_DESIGN__NO_GRAVITON_PREMISE__AWAITS_RF3_COMPLETE_SOURCE`

**Purpose:** isolate the shortest rigorous theorem that turns the qualified
record-derived metric and a complete conserved RF3 source into a positive
Einstein--Hilbert stiffness.  This is a downstream gravity calculation, not a
new record mechanism.

## 1. Result first

The strict microscopic RF4 calculation is not a search for a graviton and is
not a second reconstruction of the metric.  It is the extraction of one local
two-derivative coefficient from the **complete same-parent metric effective
action**:

\[
 \boxed{
 C_R^{\rm eff}
 =\mathscr P_R\!\left[\Gamma^{(2)}_{\Delta c,\,\rm complete}\right]>0.}
 \tag{RF4D01}
\]

Here `Delta/c` are the closed-time-path difference/average metric variables,
and `P_R` is a prospectively normalized projection onto the quadratic symbol
of `int sqrt(-g) R`.  The kernel must contain the microscopic collective
action, every integrated fast field, all constraints, measures, reservoirs,
writers, supports, ports, boundary terms, second source derivatives, and
matching terms exactly once.

If the same parent has already earned a local four-dimensional, parity-even,
diffeomorphism-covariant, metric-only conservative action through two
derivatives, then (RF4D01) fixes the entire nonlinear leading action:

\[
 \boxed{
 \Gamma_{\rm grav}^{(0,2)}[g]
 =C_R^{\rm eff}\int d^4x\sqrt{-g}
       (R-2\Lambda_{\rm eff}),\qquad C_R^{\rm eff}>0.}
 \tag{RF4D02}
\]

No microscopic particle, helicity-two pole, soft theorem, or graviton is a
premise.  A linear tensor-particle description can be recovered only after
(RF4D02) is established.

There is also a shorter actual-world branch.  Once RF3 and the leading-action
classification have independently identified the same RGRL metric, the
complete conserved source, and the Einstein--Hilbert **form**, the observed
attractive Newtonian response may calibrate the positive matched total
coefficient.  That branch proves the form and identifies the realized
response as gravity without claiming a parameter-free microscopic derivation
of the coefficient.  Section 10 states its noncircularity conditions.

## 2. Why one coefficient is enough, and why its sign is not automatic

On a connected four-dimensional phase, a local scalar action of one metric
with at most two derivatives has, modulo a boundary term, the form

\[
 S_{(0,2)}[g]
 =\int d^4x\sqrt{-g}\,[C_0+C_RR].                 \tag{RF4D03}
\]

Equivalently, its natural divergence-free metric Euler tensor is a linear
combination of `g_mn` and `G_mn`.  This is the RIEHB/Lovelock uniqueness step.
Full nonlinear covariance therefore makes the coefficient measured in any
one non-null curvature direction the common coefficient of the complete
`sqrt(-g)R` functional.

Covariance does **not** fix that coefficient's sign or make it nonzero.  For
any otherwise lawful completion one may add

\[
 \lambda\int d^4x\sqrt{-g}R                       \tag{RF4D04}
\]

with either sign while preserving record formation, metric kinematics, and
the diffeomorphism Ward identity.  Likewise, a healthy unitary fast field can
make a negative, zero, or positive contribution depending on its spin,
curvature coupling, contacts, and matching convention.  Thus the following
implication is false:

\[
 \text{qualified records + conserved source + passivity}
 \ \not\Longrightarrow\ C_R^{\rm eff}>0.          \tag{RF4D05}
\]

The sign is a genuine datum rather than a consequence of the structural
premises.  It must be obtained either from the complete microscopic
calculation above or from an independent same-metric, same-source endpoint
calibration.  It is not a missing record postulate.

## 3. Correct closed-time-path object

For the complete RF3 source family define

\[
 Z[g^+,g^-]
 =\operatorname{Tr}\!\left(U[g^+]\rho\,U[g^-]^\dagger\right),
 \qquad W=-i\hbar\log Z,                          \tag{RF4D06}
\]

and use

\[
 g_c={g^++g^-\over2},\qquad g_\Delta=g^+-g^-.    \tag{RF4D07}
\]

The physical retarded inverse-response/equation kernel is the mixed object

\[
 \mathcal K^R_{\mu\nu,\rho\sigma}
 =\left.{\delta^2\Gamma_{\rm CTP}
  \over\delta g_\Delta^{\mu\nu}\,\delta g_c^{\rho\sigma}}
 \right|_{g_\Delta=0}.                            \tag{RF4D08}
\]

It is not `Gamma_CTP[g,g]`, which vanishes by closed-time-path unitarity, and
it is not the connected susceptibility unless the required 1PI inversion and
collective-field terms have been performed.  A real local conservative term
appears as

\[
 \Gamma_{\rm CTP}^{\rm loc}[g^+,g^-]
 =S_{\rm loc}[g^+]-S_{\rm loc}[g^-].              \tag{RF4D09}
\]

Dissipative, noise, memory, and other branch-mixing terms remain in their
owned nonlocal/remainder sector.  Equations (RF4D08)--(RF4D09) prevent a
common but invalid promotion of a retarded correlator, Euclidean free-energy
curvature, or influence gamma directly into an Einstein--Hilbert action.

Because `C=<Y>` is an expectation coordinate, its collective action is obtained
by the typed CTP Legendre transform from the natural source `J`; `W_JJ` is a
connected susceptibility and `Gamma_CC` is its 1PI inverse only on a fixed
invertible physical quotient.  If other retained collective fields `chi`
mix with `C`, the metric block used in (RF4D08) is the full Schur complement

\[
 \boxed{
 K_{CC}^{\rm dr}
 =\Gamma_{CC}-\Gamma_{C\chi}
  \Gamma_{\chi\chi}^{-1}\Gamma_{\chi C},}         \tag{RF4D09a}
\]

with gauge, constraint, zero-mode, boundary, and causal-inverse conventions
fixed before inversion.  Equation (RF4D09a), not an entrywise reciprocal of
the FY response, owns feedback from every retained field.  A genuinely slow
or gapless `chi` that cannot be lawfully eliminated remains an explicit field
and takes the parent outside the metric-only premise below.

## 4. Complete source/contact formula

Let `j_A` be the complete metric source coordinates and

\[
 H_A:={\partial H\over\partial j_A}=-{1\over2}Q_A,
 \qquad H_{AB}:={\partial^2H\over\partial j_A\partial j_B}.
 \tag{RF4D10}
\]

For a finite thermal parent, the exact zero-frequency Euclidean Hessian is

\[
 \boxed{
 {\partial^2F\over\partial j_A\partial j_B}
 =\langle H_{AB}\rangle_\beta
 -\int_0^\beta d\tau\,
  \langle\delta H_A(-i\tau)\,\delta H_B(0)\rangle_\beta.}
 \tag{RF4D11}
\]

The first term is the second-source-derivative contact/seagull.  The second is
the noncommuting Kubo--Mori stress covariance.  The real retarded kernel is
the analytically continued CTP counterpart with its equal-time contacts.
Neither term may be omitted.  In particular, positivity of the covariance
does not fix the sign of their difference.  Equation (RF4D11) supplies the
complete connected/source curvature entering the Legendre and Schur steps;
it is not by itself the final collective 1PI stiffness.

The same formula is evaluated after the declared fast/slow split.  A local
preterm, fast determinant, constraint/measure term, and matching term are
summed before the sign is scored:

\[
 C_R^{\rm eff}
 =C_R^{\rm collective}+C_R^{>}+C_R^{\rm constraint/measure}
  +C_R^{\rm boundary/match}+C_R^{\rm pre}.         \tag{RF4D12}
\]

Strict induced origin additionally requires a physically distinguished
same-parent prescription fixing `C_R^pre=0`.  That stronger statement is not
needed to identify a positive gravitational response; a positive matched
total is sufficient.

## 5. Record-coordinate pullback and coefficient extractor

Let `C_e=<Y_e>` be the six qualified observable pair-memory coordinates, and
let

\[
 D_C:\mathbb R^6\overset\cong\longrightarrow
       \operatorname{Sym}^2(V)                    \tag{RF4D13}
\]

be the exact PMICS/PMSR metric tangent after the physical RGRL/F3 solder and
scale calibration.  For a nonzero spatial covector `k`, choose a pair vector
`t` such that

\[
 h_t=D_Ct,\qquad Ph_tP\ne0,                       \tag{RF4D14}
\]

and preferably choose one of PMICS's two transverse trace-free tidal
directions.  Let `K_CC` be the complete 1PI `Delta/c` kernel pulled back to
the pair-memory fields.  Define the reference kernel `K_R` by the second
variation of

\[
 I_R[g]=\int\sqrt{-g}R                            \tag{RF4D15}
\]

in exactly the same metric, Fourier, Wick/retarded, boundary, and symmetric-
tensor conventions.  On a near-flat admitted patch, subtract the
zero-derivative contribution and set

\[
 \boxed{
 C_R(t,k)=
 {\langle t,[K_{CC}(0,k)-K_{CC}(0,0)]t\rangle
  \over
  \langle D_Ct,[K_R(0,k)-K_R(0,0)]D_Ct\rangle}.}  \tag{RF4D16}
\]

The denominator is nonzero for (RF4D14).  The RIEHB Lorentzian convention is
used to orient its sign, so `C_R>0` is exactly positive physical metric
stiffness rather than an untracked Wick-rotation convention.

For a finite cell/refinement family use `C_R^(N)(t,k_N)`.  The pass condition
is a common limit with a prospective positive lower bound

\[
 \sup_{t,\hat k}\left|C_R^{(N)}(t,k_N)-C_R^{\rm eff}\right|
 \le\epsilon_N,
 \qquad \epsilon_N\to0,
 \qquad C_R^{\rm eff}\ge c_*>0,                  \tag{RF4D17}
\]

and a nonempty band on which every four-derivative, nonlocal, anisotropic,
boundary, and extra-field remainder is smaller than the leading `k^2` term.
The supremum is taken on the normalized physical curvature quotient, not on
the three pure-gradient directions.

Full nonlinear covariance makes one nonzero normalized direction algebraically
sufficient to fix `C_R`.  The executable audit nevertheless scores both PMICS
tidal directions, the trace/constraint direction, several momentum
orientations, and successive refinements.  Disagreement is a falsifier of the
claimed common metric/covariant limit rather than permission to average the
coefficients.

## 6. Theorem `RCEHS` -- record-conditioned Einstein--Hilbert stiffness

Assume one same-parent family satisfies:

1. **qualified metric ancestry:** observable retained pair-memory fields have
   the authenticated formation/retention/distinguishability/lineage custody
   of URFT and the physical metric/curvature realization of RF1--RF2;
2. **complete RF3 source:** `H[g_00,g_0i,g_ij;ports]` supplies all ten source
   components, second derivatives, currents, work, constraints, boundaries,
   and contacts before projection, and its complete CTP functional satisfies
   the off-shell Ward identity and the controlled common-cone/refinement
   limit;
3. **leading action class:** the conservative zero/two-derivative sector is
   local, four-dimensional, parity even, diffeomorphism covariant, and
   metric-only after lawful algebraic elimination; every retained additional
   field and every branch-mixing/nonlocal term is separately owned and
   bounded;
4. **unique matching:** the complete 1PI and fast/slow contribution census is
   finite and counted exactly once in one prospectively fixed physical
   matching prescription; and
5. **positive coefficient calculation:** (RF4D16)--(RF4D17) hold with
   `C_R^eff>=c_*>0` and a controlled derivative remainder on a nonempty band.

Then the conservative leading metric action is exactly (RF4D02), modulo the
declared boundary term.  Its coefficient is common to all admitted metric
directions and all complete-stress sectors.  Together with EX's six spatial
stationarity equations, complete on-shell Ward identity, and independently
owned initial constraints, the full metric residual is the Einstein residual
plus the bounded remainder.

### Proof

Premises 1--2 identify the pair-memory deformation with the physical metric
variation and make its complete `Delta/c` Hessian a same-parent source
derivative rather than a joined correlator.  Premise 3 and the local
two-derivative invariant classification give (RF4D03).  The quotient
isomorphism of PMICS makes the denominator in (RF4D16) nonzero on the scored
tidal directions, so premises 4--5 fix the unique Ricci coefficient and its
positive sign.  The common limit and remainder bound give (RF4D02) on the
declared band.  EX then turns the six spatial equations into all ten metric
equations under its independently stated constraint premises.  No particle
premise enters. QED.

## 7. What is already proved

1. **PMICS/RF1:** the observable pair-memory metric tangent has exactly three
   nonzero-momentum intrinsic-curvature quotient directions and three pure
   spatial-gradient directions.  This supplies the lawful probes `t` in
   (RF4D14).
2. **PMSR/RF2, finite commuting ceiling:** observable pair expectations, the
   Fisher metric, and the physical pair strain source are exact mixed
   derivatives of one finite Gibbs generating functional under the complete
   DPAR premise.  The noncommuting full-F3 lift remains open.
3. **FY/FZ:** the finite projected F3 parent has a six-channel spatial source,
   a two-dimensional active TT quotient, and exact finite gapped response.
   It lacks the temporal/current/contact-complete Ward family required by
   premise 2.
4. **CTP/KMS lanes:** a complete influence-gamma functional plus KMS and
   causal data can reconstruct the dissipative retarded kernel, but only up
   to the independently owned contact polynomial.  This supports the response
   method and proves why gamma or passivity alone cannot close RF4.
5. **EP:** tetrahedral symmetry permits independent `A1`, `E`, and `T2` CTP
   kernels.  It does not force isotropy, a tensor pole, or an EH coefficient.
6. **RIEHB:** on an earned common metric, a complete covariant fast shell with
   positive matched `C_R` yields the full nonlinear EH action and
   back-reaction.  RCEHS is the shorter direct-CTP coefficient route to that
   same gate; it does not repeal RIEHB's remainder and matching guards.
7. **EX:** six local spatial metric equations plus the complete on-shell Ward
   identity and prospectively zero initial constraints give the full metric
   equation.
8. **EY/FA:** the declared visible Standard-Model shell has a positive partial
   Ricci contribution at the table-conformal Higgs value.  This is not the
   complete total in (RF4D12).
9. **FF:** a q4 basis isometry is not a physical q4/F3 support solder, and the
   old same-incidence carrier/ice construction is incompatible.  RCEHS may be
   applied only after RF2--RF3 establish one actual same-parent metric/source
   family; it cannot join a q4 metric to an unrelated F3 response by matching
   dimensions.

## 8. Exact remaining calculation

After RF3 closes, RF4 needs one bounded calculation, not a new architecture:

1. form the complete finite `H_N[g^+,g^-;ports]` and its exact first and second
   metric derivatives;
2. construct the connected `Delta/c` source Hessian with (RF4D11), perform the
   typed `J <-> C` Legendre inversion, and form the complete Schur complement
   (RF4D09a), retaining the collective-field action and declared fast/slow
   matching;
3. pull it back and push it forward through the audited `C <-> g` isomorphism;
4. evaluate (RF4D16) on the exact PMICS curvature quotient for several
   momenta/orientations and on at least three refinements;
5. prove a uniform common-coefficient bound, `C_R>=c_*>0`, and an explicit
   `O(k^4/kappa_R^2)` or stronger remainder bound;
6. replay the result independently from source matrices, not from a copied
   coefficient table; and
7. only then compose with EX/RIEHB and promote RF4.

For a finite Hamiltonian this can be implemented with exact symbolic source
derivatives, exact or interval-certified diagonalization/Lehmann sums, and
rigorous finite-size/remainder bounds.  A laboratory is not logically needed
for this **model-conditional** theorem.  A laboratory or qualifying external
data are needed to establish that nature realizes the selected parent and its
matching/crossover packet.

## 9. Hidden conventional assumptions that must remain visible

1. A conserved source does not imply positive stiffness.
2. CTP unitarity makes the equal-branch action vanish; the mixed `Delta/c`
   kernel, not `Gamma[g,g]`, carries the equation of motion.
3. A connected retarded susceptibility is not automatically the 1PI inverse
   metric kernel.
4. Legendre inversion must be performed on the physical gauge/constraint
   quotient and with every mixed retained-field block included.
5. Spectral/passivity positivity fixes dissipative residue signs, not the
   real local `k^2` contact coefficient.
6. A Euclidean sign cannot be imported into the Lorentzian `C_R` convention
   without a fixed analytic-continuation dictionary.
7. Locality, a smooth manifold, a common cone, and metric-only leading order
   are macroscopic theorem premises to be earned by RF3/refinement, not
   consequences of the word “record.”
8. Integrating out a gapless retained field can produce a nonlocal metric
   functional or an extra leading mode; it cannot be called algebraic
   elimination.
9. A one-momentum finite-graph `k^2` fit is not a continuum coefficient.
10. A visible-sector heat-kernel subtotal is not the complete matched total.
11. `C_R^pre=0` is a stronger strict-origin/matching claim; positive total
    stiffness is enough for the gravity-identification theorem.

## 10. Short actual-world endpoint-matching branch

The first-principles coefficient calculation in sections 4--8 is sufficient
but is not logically necessary for identifying the realized actual-world
response as gravity.  A shorter branch is available if RF3 and the
leading-action theorem close first.

### Theorem `RCEHS-END` -- positive total stiffness by independent Newton matching

Assume, independently of a numerical fit to Newton's law:

1. RF1--RF3 identify one qualified record-conditioned RGRL metric `g_R`, its
   physical scale, and one complete variational stress
   `T_mn^complete` in the same parent;
2. the complete off-shell Ward identity, locality, four-dimensional
   metric-only zero/two-derivative class, and remainder bounds give the form

   \[
    \Gamma_{\rm grav}^{(0,2)}
    =\int\sqrt{-g_R}\,[C_R^{\rm eff}R+C_0^{\rm eff}],
                                                        \tag{RF4D18}
   \]

   with the sign and magnitude not yet fixed;
3. source energy/mass is calibrated inertially or nongravitationally from the
   already-complete stress ledger, rather than defined from the gravitational
   response being scored;
4. rods, clocks, causal propagation, and the weak-field potential/acceleration
   are referred to the already-identified `g_R`, rather than to a metric whose
   scale was chosen to make Newton's law hold;
5. the scored static response is attractive and has the operational Newton
   normalization `G_obs>0` on a nonempty band; and
6. no unsuppressed scalar, vector, second-metric, affine, nonlocal, boundary,
   or other leading channel contributes to that same static normalization;
   all permitted corrections satisfy the frozen remainder bound.

Then the weak-field reduction of (RF4D18) and the independent endpoint
measurement give, in natural units,

\[
 \boxed{
 C_R^{\rm eff}={1\over16\pi G_{\rm obs}}>0.}       \tag{RF4D19}
\]

and `Lambda_eff=-C_0^eff/(2C_R^eff)` when the volume coefficient is retained.

With `x^0=ct`, the corresponding SI action coefficient is

\[
 \boxed{
 C_{R,\rm SI}^{\rm eff}={c^3\over16\pi G_{\rm obs}}>0.}       \tag{RF4D20}
\]

This closes the positive **total** RF4 coefficient for gravity identification.
It does not decompose that total into microscopic, fast-shell, pre-existing,
measure, boundary, or matching pieces.

### Proof

The leading-action classification supplies the EH family before the endpoint
normalization is inspected.  Its linear static equation for the already
normalized complete stress has one free coefficient, `C_R^eff`.  Premises
3--5 compare that prediction with an independently normalized attractive
source-to-metric response and therefore fix the coefficient by (RF4D19).
Premise 6 prevents a scalar/extra-field admixture from making the measured
Newton coefficient a different combination of couplings.  Positivity of
`G_obs` fixes the healthy sign in the frozen RIEHB convention. QED.

### Three distinct claims

This branch separates three results that should not be made to stand or fall
together.

1. **Form and record-conditioned mechanism.**  RF1--RF3 plus the
   local/covariant metric-only classification prove that qualified record
   lineage constitutes the metric tangent and complete source whose leading
   response has the EH form.  EX supplies the full equation from the six
   spatial equations, Ward identity, and owned constraints.  This is the
   theoretical gravity-formation result.
2. **Empirical coefficient matching.**  The observed same-metric,
   same-source attractive Newtonian response fixes the one remaining total
   coefficient to (RF4D19).  This is an actual-world calibration, analogous
   in logical role to inheriting the observed electromagnetic coupling on the
   already-identified electromagnetic sector.  It is not a weakness or an
   unproved sign assumption.
3. **Parameter-free or strict induced-origin derivation.**  Deriving
   `G_obs` numerically from F3 parameters, proving the complete microscopic
   split (RF4D12), and physically fixing `C_R^pre=0` remain deeper results.
   They are not required to establish that the record-conditioned response
   has the form, positive sign, and observed strength of gravity.

### Circularity audit

The endpoint branch is noncircular only if all of the following separations
are respected.

- The EH form must follow from the complete Ward/locality/metric-only theorem,
  not be assumed because the data were reduced with general relativity.
- The source mass/energy normalization must come from the complete independent
  stress ledger, not from `F=Gm_1m_2/r^2` itself.
- The RGRL metric and its scale must be identified by the record/causal-volume,
  clock, and common-probe construction before `G_obs` is fitted.
- The endpoint observable should be the underlying calibrated weak-field
  response, or an accepted `G_obs` value used explicitly as an empirical
  input.  A tabulated value whose extraction assumed the target model is a
  calibration within that model, not independent evidence for the EH form.
- The measured coefficient is the matched total.  EY/FA, a fast determinant,
  a collective contribution, or a putative bare term may not be added to it
  again.
- Extra leading fields must be excluded or jointly fitted before identifying
  Newton's coefficient with `1/(16 pi C_R)`.  Otherwise `G_obs` can be a
  mixture of the metric stiffness and additional charges.

Under these guards there is no same-metric/source circularity: the theory
first identifies what metric and source are being related and restricts the
law to a one-coefficient EH family; the observation then supplies that
coefficient.  What remains open is why the microscopic parent realizes that
number, not whether the resulting positive response is gravity.

## 11. Decision boundary

No theory decision is needed if the microscopic coefficient is positive, if
the guarded endpoint branch closes, or if a specific omitted
contact/remainder defect is found and repaired within the same parent.  A
genuine theory decision arises only if both lawful closure branches fail--for
example, the complete audited same-parent coefficient is zero or negative,
the actual endpoint cannot be joined to that same metric/source without an
unsuppressed extra field, and two inequivalent non-ad-hoc physical parent
repairs remain after the source, matching, and refinement census is exhausted.

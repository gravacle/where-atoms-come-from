# Record-mediated localization derivability-boundary theorem

**Lane ID:** `CROSS-RFT-GRA-EZ-RECORD-MEDIATED-LOCALIZATION-BOUNDARY-V001`

**Official short name:** `RMLB`

**Date:** 2026-08-27

**Status:** `BOUNDED_PROOF_DRAFT__INDEPENDENT_AUDIT_PENDING`

**Claim class:** exact implication from AURFT/U-DCL plus terminal-record
localization; exact Fisher data-processing corollary; exact finite quantum
countermodel to a universal classical localization query; exact noisy-record
countermodel to query/QFI saturation; exact metric-spectator countermodel to
proper-metric soldering

**Not claimed:** that the actual world uses either countermodel; that records
are irrelevant to localization; that CQLTS is false; that the EU/FERS model is
not a valid sufficient class; or that a classical information metric is
spacetime without a physical soldering law

## 1. Exact question and operational premise

AURFT proves universal record coverage inside the adopted U-DCL physical axiom
system.  U-DCL is deliberately metric-free: it supplies a complete acyclic
incidence, one apparatus-inclusive joint separator, exact complete
instruments, and incidence-respecting future composition for each admitted
record.  It does not name a Lorentz metric.

Add the proposed operational fact:

> **TRL -- terminal-record localization.** Every scored localization,
> displacement, clock, matter-probe, or electromagnetic-probe measurement in
> the admitted finite mission terminates in an independently bona-fide record,
> and the reported score is obtained from that record's registered query by a
> prospectively declared, parameter-independent read kernel.

This is the formal reading of “terminates in a bona-fide record **of the
localization score**.”  If “terminates in” means only that some unrelated
record occurs later in time, not even record mediation follows; temporal
succession alone supplies no kernel in (EZ01).

For a measurement setting/instrument `m`, let `x` denote the physical
configuration parameter, `z` the fixed context, `Q_m` the complete raw
registered query outcome for that one instrument, and `Y_m` its reported
localization score.

The question is which of the following follow from AURFT/U-DCL plus TRL:

1. every scored localization observation is record-mediated;
2. one classical record query is sufficient for **all** localization
   measurements;
3. its Fisher tensor equals the state QFI and is the physical metric; and
4. that metric equals the causal/proper-time metric followed by matter and EM.

## 2. Theorem RMLB-1 -- exact per-measurement record mediation

Assume AURFT/U-DCL and TRL.  Then, for every admitted scored localization
measurement `m`, there is one U-DCL physical witness for its terminal record
and a Markov kernel `K_m` such that

\[
 \boxed{
 \Pr(Y_m\in B\mid x,z,m)
 =\int K_m(B\mid q)\,p_m(dq\mid x,z)}              \tag{EZ01}
\]

for every measurable score set `B`.  The distribution `p_m` is the registered
query marginal of the same apparatus-inclusive separator and complete
instrument supplied by that record's DCL witness.  If the registered outcome
is itself the score, `K_m` is the identity kernel.

### Proof

TRL supplies an independently admitted record `r_m`, its registered query
`Q_m`, and the declared read kernel `K_m`.  U-DCL supplies
`DCL_phys(r_m)`.  D2 supplies one positive normalized joint separator including
the apparatus, history, carrier, environment, and query degrees; D3 supplies
the complete instrument producing `Q_m`; and D4 prevents an undeclared direct
source-label bypass.  Marginalizing the same joint object gives `p_m` and
applying the declared read kernel gives (EZ01).  No outcome-conditioned
separator or cross-instrument branch mixing is used. QED.

Equation (EZ01) is stronger than the bare statement that a record happens
after a measurement: it gives exact same-parent causal/statistical custody for
the score.  It is also strictly weaker than FERS-P2.  DCL **completeness** means
that no outcome or physical port of this chosen instrument is omitted.
Statistical **sufficiency** in FERS means that one and the same classical query
Blackwell-dominates every admitted localization instrument.  The first does
not imply the second.

There is one additional conclusion inside a **single fixed DCL mission**.  If
several localization settings belong to the complete post-frontier instrument
family of the same witness, D2--D4 give one common pre-query statistical
experiment:

\[
 \boxed{
 p_m(dq\mid x,z)
 = {\cal J}^{b}_m[\Omega^b_{x,z}](dq)
 \quad\text{for every admitted }m.}                \tag{EZ01a}
\]

For `b=K`, `${\cal J}_m` is a stochastic instrument acting on the classical
separator law; for `b=W`, it is a CP instrument acting on the common C*-state.
This is a genuine common **separator experiment**, and is the strongest
same-mission sufficiency supplied directly by U-DCL.  It is not yet one
physically readable classical record query.  The classical separator can
contain inaccessible crossing variables, while the quantum separator can
require mutually incompatible reads.  Different finite missions also need not
share one separator family.

## 3. Corollary RMLB-2 -- the information that follows automatically

Suppose `p_m(q|x)` is differentiable, dominated, and has parameter-independent
support on the declared tangent.  Its classical Fisher tensor is

\[
 F^{Q_m}_{ij}(x)
 =\int \partial_i\log p_m(q|x)\,
        \partial_j\log p_m(q|x)\,p_m(dq|x).         \tag{EZ02}
\]

The Markov data-processing theorem applied to (EZ01) gives

\[
 \boxed{F^{Y_m}\preceq F^{Q_m}.}                   \tag{EZ03}
\]

In the W branch, let the complete pre-query separator states, with every
parameter-dependent apparatus degree included, form a differentiable family
`rho_x`.  Require the frozen query instrument to induce one
parameter-independent POVM `M_m(dq)` such that
`p_m(dq|x)=Tr[rho_x M_m(dq)]`.  The Braunstein--Caves bound then gives,
tangent by tangent,

\[
 \boxed{F^{Q_m}\preceq F^{\rm SLD}_{\rho}.}        \tag{EZ04}
\]

Thus terminal recordhood supplies a legitimate per-instrument Fisher tensor
whenever the ordinary regularity conditions hold, and U-DCL supplies its
physical custody.  Neither inequality is generally an equality.  AURFT/U-DCL
does not impose Chentsov naturality across the category of localization
experiments, one common normalization, a sufficient recovery map, or QFI
saturation.

## 4. Countermodel RMLB-C1 -- terminal records do not imply one complete query

Consider a finite W-branch world with a two-dimensional carrier whose admitted
source family is tomographically complete (for example, all qubit density
operators, or any four density operators whose real linear span is
`Herm(2)`).  After the U-DCL frontier a
freely chosen local setting `m in {Z,X}` selects one of the sharp instruments

\[
 P^Z_\pm={I\pm\sigma_z\over2},\qquad
 P^X_\pm={I\pm\sigma_x\over2}.                     \tag{EZ05}
\]

Each outcome is copied into an orthogonal stable pointer bit and queried later.
The source-to-pointer incidence is a finite DAG; the carrier, setting path,
apparatus, pointer, failure outcome, and query are all included in one joint
state and complete CP instrument.  The future instrument receives the freely
chosen setting on its declared local wire and the source state only through
the carrier.  Hence every pointer record has a direct `DCL_phys` witness and
TRL holds.  A toy world in which all records are generated by these finite
acyclic pointer instruments satisfies U-DCL.

Suppose nevertheless that one classical standard-Borel query POVM `M(dk)` were
sufficient for both sharp readings on this tomographically complete family.
Then measurable stochastic kernels `a(i|k)` and `b(j|k)` reproduce the two
outcome probabilities for every admitted state.  Tomographic completeness
upgrades those equalities of probabilities to the operator identities

\[
 P^Z_i=\int a(i|k)M(dk),\qquad
 P^X_j=\int b(j|k)M(dk).                           \tag{EZ06}
\]

Define

\[
 G_{ij}:=\int a(i|k)b(j|k)M(dk).                  \tag{EZ07}
\]

`G` is a joint POVM with marginals `P^Z` and `P^X`, and
`0 <= G_ij <= P^Z_i` as well as `0 <= G_ij <= P^X_j`.  A positive operator
bounded by a projection is supported inside that projection.  Hence the two
sharp marginals give

\[
 P^Z_iP^X_j=G_{ij}=P^X_jP^Z_i.                    \tag{EZ08}
\]

But

\[
 [P^Z_+,P^X_+]={i\over2}\sigma_y\ne0,             \tag{EZ09}
\]

a contradiction.  Hence all localization readings can terminate in bona-fide
U-DCL records while no single classical terminal query Blackwell-dominates
them.  The common quantum carrier state does not repair the conclusion: it is
not a broadcast classical record query, and extracting the two sharp reads
requires incompatible instruments.

This countermodel is decisive against

\[
 \text{AURFT/U-DCL + TRL}
 \Longrightarrow\text{FERS-P2 complete-query sufficiency}.   \tag{EZ10}
\]

The missing query law is independent of metric matching.  For example, retain
the full tomographically complete source family above as the domain on which a
putative universal query must work, and use the real pure-state curve
`|psi_x>=cos(x/2)|0>+sin(x/2)|1>`, `0<x<pi`, only as its metric-matched slice.
The recorded Z query has constant Fisher information one on that slice and may
be prospectively soldered to a constant one-dimensional proper metric.  The
separately recorded sharp X instrument is still incompatible with Z.

The incompatibility already follows even if one tests the putative common
query only on this real curve.  If its Z postprocessing has positive effect
`A` and reproduces `P^Z_+` for every `x`, the Bloch expansion gives
`A=P^Z_+ + c\sigma_y`, because the real curve spans the `I`, `\sigma_x`, and
`\sigma_z` expectation directions.  But
`det(P^Z_+ + c\sigma_y)=-c^2`; positivity forces `c=0`, so `A=P^Z_+`.
The same argument for the X postprocessing gives `B=P^X_+`.  A common query
would therefore jointly measure the two incompatible sharp effects, which is
impossible by (EZ07)--(EZ09).  Thus soldering one query does not create a
universal classical sufficient query.

## 5. Countermodel RMLB-C2 -- record mediation does not saturate QFI

Let

\[
 \rho_t={I+t\sigma_z\over2},\qquad -1<t<1,         \tag{EZ11}
\]

and query it with the informative noisy binary POVM

\[
 E_y={I+y\eta\sigma_z\over2},qquad
 y\in\{-1,+1\},\quad 0<\eta<1.                    \tag{EZ12}
\]

Copy `y` into the same kind of stable finite pointer record.  This again obeys
U-DCL and TRL.  Direct calculation gives

\[
 p_t(y)={1+y\eta t\over2},qquad
 F^{Q}(t)={\eta^2\over1-\eta^2t^2},qquad
 F^{\rm SLD}(t)={1\over1-t^2}.                    \tag{EZ13}
\]

Therefore

\[
\boxed{
 F^{\rm SLD}(t)-F^Q(t)
 ={1-\eta^2\over(1-t^2)(1-\eta^2t^2)}>0.}         \tag{EZ14}
\]

The squared-fidelity separation is also exact.  For `t=0`, `t'=u` with
`0<|u|<1`, the commuting state family and its noisy recorded query give

\[
 \gamma_{\rm state}(0,u)
 ={1+\sqrt{1-u^2}\over2},\qquad
 \gamma_Q(0,u)
 ={1+\sqrt{1-\eta^2u^2}\over2}
 >\gamma_{\rm state}(0,u).                         \tag{EZ14a}
\]

The record is bona fide and informative, the instrument is complete, and the
query is exactly registered, yet it loses state information.  Hence neither
`F_cl=F_SLD` nor query-fidelity/state-fidelity equality follows from U-DCL or
terminal recordhood.  EU proves those equalities for its commuting sufficient
family; actual-world use requires an instantiation or saturation premise.

## 6. Countermodel RMLB-C3 -- information geometry is not proper geometry

Take a finite K-branch localization device with physical coordinate
`x in (-1,1)` and complete recorded bit

\[
 p_x(+)=\frac{1+x}{2},\qquad
 p_x(-)=\frac{1-x}{2}.                              \tag{EZ15}
\]

Its Fisher metric is

\[
 F(x)={1\over1-x^2}.                                \tag{EZ16}
\]

Place the finite apparatus and ordinary matter/EM probes on the Lorentzian
domain

\[
 ds^2=-d\tau^2+L^2dx^2+dy^2+dz^2.                 \tag{EZ17}
\]

Let every device and probe score be written to a finite stable pointer record.
All missions can have complete acyclic K witnesses, so AURFT/U-DCL and TRL
hold.  Within the binary device, (EZ15) is already the complete raw query.  But
there is no one constant information-to-length scale `ell_F` for which

\[
 L^2=\ell_F^2F(x)                                  \tag{EZ18}
\]

on any nontrivial open `x` interval.

At a declared finite resolution, every admitted clock/matter/EM localization
score in this toy world may be a fixed copy or stochastic postprocessing of
the same raw bit.  The bit is then Blackwell-complete for that whole admitted
score family, while (EZ18) still fails.  Thus universal score sufficiency does
not by itself solder the statistical metric to the proper metric.

More generally, the same DCL incidence and record channel can be placed in
direct product with different smooth Lorentzian spectator metrics, including
nonconstant conformal rescalings that preserve the causal order seen by DCL
while changing proper lengths and times.  Matter and EM may consistently
follow the selected spectator metric.  U-DCL's truth value and (EZ15)--(EZ16)
do not change, because U-DCL constrains incidence and composability rather than
a line element.

Thus even an informative complete terminal record query does not make its
Fisher geometry the physical proper metric.  A same-parent soldering law is
logically independent.

## 7. Strongest derivable chain and current coarse-interface residual

The exact implication boundary is

\[
\boxed{
 \begin{aligned}
 &\text{AURFT/U-DCL + TRL}\\
 &\quad\Longrightarrow
 \left\{\begin{array}{l}
   \text{per-instrument record factorization (EZ01)},\\
   \text{one common typed separator experiment inside each fixed DCL mission
   (EZ01a)}
 \end{array}\right.\\
 &\quad\xRightarrow{\ \text{ordinary smoothness}\ }
   \text{per-instrument Fisher tensor and DPI (EZ02)--(EZ04)},
 \end{aligned}}                                    \tag{EZ19}
\]

but the same antecedent does not imply a single sufficient classical query,
QFI saturation, Chentsov selection as the operational metric, or equality to
the causal/proper metric.

Consequently the original `RW-L` premise can be narrowed, but not eliminated.
Record mediation itself is no longer residual.  The shortest currently
exhibited noncircular coarse-interface sufficient addition is one two-clause
law:

> **RLS -- record-localization soldering law.** On one admitted same-parent
> infrared domain: (A) one authenticated terminal record query is
> Blackwell-sufficient for every localization read made by clocks, matter, EM,
> and independent probes, and in the W branch it is sufficient on the declared
> state family in the stronger sense needed to preserve squared fidelity and
> saturate SLD QFI; and (B) the operational infinitesimal proper
> metric belongs to the full FERS-P3/P4 naturality class: regularity under
> congruent Markov embeddings and sufficient inverses, stochastic contraction,
> product additivity, and one common normalization.  One prospectively fixed
> scale/time bind supplies the physical units, and FERS-P6 requires actual
> transported frames to implement the proper time-oriented length-preserving
> face maps.

Clause A is not hidden inside the word `complete`; RMLB-C1 and C2 refute that
shortcut.  Clause B is not hidden inside the word `record`; RMLB-C3 refutes
that shortcut.  Under clause A, the actual query/state gamma and Fisher/QFI
types join.  Under clause B, Chentsov plus FERS selects

\[
 q=\ell_F^2F,
\]

and FERS edge/face uniqueness identifies physical transport.  If the actual
EU family is independently instantiated, its exact saturation theorem can
supply the W-branch part of clause A; it still does not supply clause B.

RLS is therefore the shortest currently exhibited **coarse-interface
sufficient residual within the current EU/FERS/induced-action route**.  C1/C2
and C3 prove that its query/type-join and metric/transport interfaces are
logically independent.  They do not prove that each clause's exact subpremise
list is unique, or that no different gravity route can use a different law.
RLS is one falsifiable physical law, not a new record-formation axiom and not a
definition of spacetime by information.

## 8. Disposition

`AURFT_PLUS_TERMINAL_RECORD_LOCALIZATION_PROVES_EXACT_PER_INSTRUMENT_RECORD_MEDIATION_COMMON_TYPED_SEPARATOR_EXPERIMENT_WITHIN_EACH_FIXED_DCL_MISSION_AND_FISHER_DATA_PROCESSING__DOES_NOT_PROVE_ONE_PHYSICALLY_READABLE_BLACKWELL_COMPLETE_CLASSICAL_QUERY__DOES_NOT_PROVE_QUERY_STATE_FIDELITY_OR_QFI_SATURATION__DOES_NOT_PROVE_INFORMATION_METRIC_EQUALS_CAUSAL_PROPER_METRIC__RECORD_LOCALIZATION_SOLDERING_RLS_IS_THE_SHORTEST_CURRENTLY_EXHIBITED_COARSE_INTERFACE_TWO_CLAUSE_SUFFICIENT_RESIDUAL_WITHIN_THE_CURRENT_ROUTE`

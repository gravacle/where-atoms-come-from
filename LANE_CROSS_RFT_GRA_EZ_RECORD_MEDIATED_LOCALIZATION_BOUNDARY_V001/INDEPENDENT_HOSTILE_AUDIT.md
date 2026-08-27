# Independent hostile audit: RMLB

**Lane:** `CROSS-RFT-GRA-EZ-RECORD-MEDIATED-LOCALIZATION-BOUNDARY-V001`

**Date:** 2026-08-27

**Audit mode:** independent no-edit review of the stabilized theorem and
verifier; this artifact was written by the hostile auditor after the builder
bytes stabilized

**Audited theorem SHA-256:**
`c6095e1374268c9afecd8179cf622c612c2d121a9eab43b8886f9cca3ff86ff4`

**Audited verifier SHA-256:**
`4a6788c25b34d986d2ad4e8a5ed604589b2f2716f34fcaaa6cd070d7427d2702`

## Verdict

**ACCEPT / CLEAN at the declared derivability-boundary claim class.**

The theorem correctly derives per-instrument record mediation and, inside one
fixed DCL mission, a common typed pre-query separator experiment. It correctly
shows that these conclusions do not imply one physically readable classical
query sufficient for every localization instrument, query/state fidelity or
Fisher/QFI saturation, or equality between an information metric and the
causal/proper metric.

This audit does not promote RLS from a residual physical law to a consequence
of AURFT/U-DCL, does not establish that the actual world satisfies RLS, and
does not establish gravity.

## Hostile tests

### 1. TRL, per-instrument factorization, and the common-separator ceiling

The implication to (EZ01) uses the strong and explicitly declared meaning of
TRL: a bona-fide terminal record **of the localization score**, together with
a prospectively fixed, parameter-independent read kernel. A merely later,
unrelated record would not imply the factorization, and the theorem says so.

The DCL witness supplies same-parent custody for the registered query and its
complete apparatus-inclusive instrument. Applying the declared Markov kernel
therefore gives the score distribution without outcome-conditioned separator
selection or a hidden source-label bypass. The theorem also keeps two notions
properly distinct: DCL completeness closes the ports of a chosen instrument;
Blackwell sufficiency would require one classical query to dominate every
admitted localization instrument.

For settings in one complete post-frontier family of one fixed DCL witness,
(EZ01a) legitimately identifies one common typed separator experiment. In the
K branch the downstream maps are stochastic instruments; in the W branch they
are CP instruments. This does not create a readable classical query: a
classical separator can include inaccessible variables, and a quantum
separator can require incompatible measurements. No cross-mission common
separator is claimed.

### 2. Fisher data processing and the Braunstein--Caves bound

Under the stated domination, differentiability, fixed-support, and
parameter-independent-kernel assumptions, Markov data processing has the
correct direction,

\[
 F^{Y_m}\preceq F^{Q_m}.
\]

In the W branch the theorem places every parameter-dependent apparatus degree
inside the pre-query state and requires the frozen query to be represented by
one parameter-independent POVM. The Braunstein--Caves inequality then has the
correct direction,

\[
 F^{Q_m}\preceq F^{\rm SLD}_{\rho}.
\]

Neither inequality supplies equality, a recovery map, Chentsov naturality, or
a shared normalization. The theorem does not infer any of those extra
conditions from terminal recordhood.

### 3. Universal-query countermodel

The qubit source family is tomographically complete. If one standard-Borel
query POVM and two stochastic postprocessings reproduced the sharp Z and X
statistics, tomography would upgrade the statistical equalities to effect
identities. The operator-valued integral

\[
 G_{ij}=\int a(i\mid k)b(j\mid k)M(dk)
\]

is then a joint POVM. Its two sharp marginals imply
`0 <= G_ij <= P^Z_i` and `0 <= G_ij <= P^X_j`; projection support gives

\[
 P^Z_iP^X_j=G_{ij}=P^X_jP^Z_i.
\]

This contradicts
`[P^Z_+,P^X_+]=(i/2)sigma_y != 0`. The standard-Borel outcome space therefore
does not evade the ordinary sharp-observable incompatibility theorem.

The metric-matched one-parameter slice is also sound. Equality of a candidate
Z effect with `P^Z_+` on the real pure-state curve leaves only an undetected
term `c sigma_y`; positivity forces it away because
`det(P^Z_+ + c sigma_y)=-c^2`. The same argument holds for X. Thus even on
that slice a common query would jointly measure the two incompatible sharp
effects. The recorded Z Fisher information can equal a constant
one-dimensional proper metric while universal classical query sufficiency
still fails.

### 4. Noisy-record QFI and fidelity separation

For the declared commuting qubit family and noisy binary POVM, direct
calculation gives

\[
 F^Q(t)={\eta^2\over1-\eta^2t^2},\qquad
 F^{\rm SLD}(t)={1\over1-t^2},
\]

and hence the strictly positive gap

\[
 F^{\rm SLD}(t)-F^Q(t)
 ={1-\eta^2\over(1-t^2)(1-\eta^2t^2)}.
\]

The fidelity convention is consistently squared fidelity. At `t=0` and
`t'=u`, the state and recorded-query values are respectively

\[
 \gamma_{\rm state}={1+\sqrt{1-u^2}\over2},\qquad
 \gamma_Q={1+\sqrt{1-\eta^2u^2}\over2},
\]

so `gamma_Q > gamma_state` for `0<eta<1` and nonzero `u`. This is the correct
data-processing direction and proves that a complete bona-fide record can be
informative without saturating state QFI or preserving state fidelity.

### 5. Metric-spectator countermodel

For the complete recorded bit `p_x(+)= (1+x)/2`, the replayed Fisher metric is
`F(x)=1/(1-x^2)`. It cannot be related by one constant scale to the constant
proper component `L^2` on a nontrivial open interval. At fixed finite
resolution the same bit may nevertheless Blackwell-dominate every admitted
score.

Taking the record channel in direct product with different smooth Lorentzian
spectator metrics preserves the DCL incidence and query statistics while
changing proper lengths and times. Positive conformal rescalings can preserve
causal order while changing the line element. This is a decisive countermodel
to deriving metric soldering from metric-free record incidence alone.

### 6. Scope and minimality of the residual law

The two RLS interfaces are independently necessary within the exhibited
coarse-interface route: the incompatible/noisy quantum examples block the
query/type-join shortcut, while the spectator-metric example blocks the
metric/transport shortcut. The theorem includes the full stated FERS-P3/P4
naturality class, a prospectively fixed physical scale/time bind, and the
FERS-P6 requirement that actual transported frames realize the proper
time-oriented length-preserving face maps.

The minimality language is correctly bounded. RLS is the shortest **currently
exhibited coarse-interface sufficient residual within the current
EU/FERS/induced-action route**. The theorem does not claim global logical
minimality, uniqueness of every subpremise, or exclusion of a different route
to gravity.

## Replay and custody

Fresh replay produced:

```text
PASS 49/49
PER_INSTRUMENT_RECORD_MEDIATION_EXACT
INCOMPATIBLE_RECORDED_QUERIES_REFUTE_UNIVERSAL_CLASSICAL_SUFFICIENCY
NOISY_BONA_FIDE_RECORD_STRICTLY_FAILS_QFI_SATURATION
METRIC_SPECTATOR_REFUTES_PROPER_GEOMETRY_DERIVATION
RLS_COARSE_INTERFACE_TWO_CLAUSE_RESIDUAL_ISOLATED
```

The duplicate lane-code collision was safely resolved by the EZ rename; the
audited theorem uses EZ consistently. No theorem or verifier byte was edited
by this audit. No dependency ledger, verification transcript, result, or
manifest existed in the draft lane at audit time; those remain builder custody
before any source freeze.

**Independent disposition:**

`CLEAN_EXACT_PER_INSTRUMENT_RECORD_MEDIATION_AND_FIXED_MISSION_COMMON_TYPED_SEPARATOR__UNIVERSAL_READABLE_CLASSICAL_QUERY_QFI_FIDELITY_AND_PROPER_METRIC_SOLDERING_DO_NOT_FOLLOW__RLS_IS_ONLY_THE_SHORTEST_CURRENTLY_EXHIBITED_COARSE_INTERFACE_RESIDUAL_WITHIN_THE_CURRENT_ROUTE__ACTUAL_WORLD_INSTANTIATION_AND_GRAVITY_REMAIN_OPEN`

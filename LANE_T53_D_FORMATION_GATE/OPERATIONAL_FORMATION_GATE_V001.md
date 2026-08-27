# Operational formation gate — V001

**Status:** FROZEN, DEFAULT-REFUTED, NOT YET APPLIED  
**Theory authority:** `PRINCIPAL_DECISION.md`  
**Unit of judgment:** one prospectively closed claim about one target state on one bona-fide actual physical record surface  
**Positive outcomes:** `PRIMARY_FORMATION_CONFIRMED` or `DERIVATIVE_FORMATION_CONFIRMED`

## 1. Purpose and boundary

This gate asks whether a declared physical target becomes a record of a declared earlier event. It does not ask merely whether a material responds to a drive. It applies the same constitutive predicate to engineered or externally driven formation and to autonomous or endogenous formation. The gate begins refuted: no positive outcome is available unless every mandatory clause below is satisfied by actual physical data and then independently reproduced.

A positive result is local to the declared surface, carrier, event class, environment, mission time, and operational range. It is not a universal-law verdict. Many local confirmations cannot by accumulation establish that the terms cover every bona-fide surface. A new bona-fide surface that requires an undeclared constitutive term or science-code branch reopens the general theory.

The one dynamical object is

`C_t = C_psi(X, W_t; p_C)`, `R_(t+) ~ K_theta(R_t, C_t, B_t, G_t; p_S)`, and `Y_t ~ M_phi(R_t, Q_t; p_I)`,

with a physical closure gate `G_t` and a provenance graph `G_X`. `C_psi` is the causal transduction or coupling path by which event/writer forcing reaches the target. `G_t` marks the measurable surface event at which incorporation closes and storage begins. `p_C`, `p_S`, and `p_I` are prospectively measurable coupling, surface, and instrument primitives, not post-hoc mechanism labels. `W` may be an engineered or exogenous interaction or an event-bounded autonomous or endogenous interaction. Persistence is the post-closure evolution `K_theta(R_t, OFF, H, CLOSED; p_S)`. Direct formation, derivative copying, read-only access, rewriting, erasure, transient response, and no change are typed causal paths through this same object. A mere retrospective taxonomy does not satisfy the gate.

`R` is representation-independent: it may be discrete, continuous, count-valued, field-valued, distributed, or mixed. “Definite physical target” means an actual per-specimen state or trajectory resolved within prospective uncertainty, not a mandatory two-state bit or point microstate. The information statistic, causal contrasts, provenance, hold test, and read test must operate on the declared representation without changing the constitutive predicate.

## 2. Closed claim tuple

Before outcome-bearing acquisition, the registration must close the tuple

`Omega = (S, T, X, P, D_causal, A, I_f, C_psi, G_t, H, tau_m, Q, Z, J, U, K_theta, M_phi, G_X, F, V_URM)`,

where:

- `S` is the actual physical surface and sample population;
- `T` is the target record-bearing physical coordinate or operational state class;
- `X` is the earlier event or event label;
- `P` is the provenance hypothesis, `PRIMARY` or `DERIVATIVE(parent_record_id)`;
- `D_causal` is the closed causal design, either `RANDOMIZED_INTERVENTION` for a controllable coupling or `PROSPECTIVE_NATURAL_CONTRAST` for an uncontrollable event;
- `A` is the design-specific allocation to `FORMATION` or `SHAM`, randomized for a controllable intervention and prospectively matched and response-blind for a natural contrast;
- `I_f` is the bounded formation interaction;
- `C_psi` is the measurable and falsifiable causal transduction or coupling path from `X` and `W` into the target dynamics;
- `G_t` is the measurable and falsifiable surface closure gate separating formation from subsequent storage;
- `H` is the common writer-off hold or storage condition;
- `tau_m` is the declared mission time measured from certified interaction cessation;
- `Q` is the read intervention and its back-action assay;
- `Z` is the event-level readout and registered physical-state interpretation;
- `J` is the prospectively selected information or proper-prediction statistic relating `Z` to `X`;
- `U` is the calibration-derived uncertainty and detection-floor system;
- `K_theta` is the frozen causal state-transition law, including writer-off hold dynamics, derived from declared measurable primitives;
- `M_phi` is the frozen readout law relating the physical state and read intervention to observed data;
- `G_X` is the event-information provenance graph, including parent-record edges and rewrite or erasure transitions;
- `F` is the complete prospective falsifier and outcome map; and
- `V_URM` is the immutable released URM version and official input contract.

No universally fixed numerical tolerance is introduced here. Each mission time, operating range, calibration uncertainty, detection floor, equivalence margin, read-disturbance bound, sample design, and decision threshold must be justified from the surface and instrument and frozen before outcome-bearing data are visible. A negative conclusion is admissible only when the registered sensitivity could have detected the claimed effect at its declared scale; otherwise the result is `UNSCOREABLE`.

## 3. Normative clauses

The following one-line clauses are normative and are mirrored byte-for-byte as strings in `OPERATIONAL_FORMATION_GATE_V001.json`.

- **FG-001**: The claim tuple, admissible operating range, population, exclusions, analysis, uncertainty system, mission time, falsifiers, and immutable URM release SHALL be content-hashed and time-fixed before outcome-bearing acquisition or, for a response-blind archival holdout, before any numerical access or unblinding.
- **FG-002**: The scored evidence SHALL be raw observations from an identified bona-fide actual physical surface; simulation, synthetic data, model output, and software self-tests SHALL NOT satisfy this clause.
- **FG-003**: The causal signal from the declared event X SHALL precede the target update it explains, and X or its independently measured forcing SHALL precede the physical closure gate and READ; labels, timing, balance or sampling law, and blind provenance SHALL be recorded independently of the target readout.
- **FG-004**: The causal design SHALL be fixed as RANDOMIZED_INTERVENTION when the formation coupling is controllable or PROSPECTIVE_NATURAL_CONTRAST when the event is uncontrollable, and no other design SHALL support a positive classification.
- **FG-005**: RANDOMIZED_INTERVENTION SHALL use concealed FORMATION versus SHAM assignment differing only in the declared coupling, while PROSPECTIVE_NATURAL_CONTRAST SHALL use independently measured forcing, prospectively matched negative time or surface controls, blind response analysis, and explicit alternative-cause falsifiers; allocation, attrition, exclusions, and deviations SHALL be fully reported in either design.
- **FG-006**: Before the formation interval, the declared target coordinate SHALL carry no resolvable information about X under the prospective statistic and calibration-derived absence floor, with sensitivity sufficient to detect an effect of the post-hold claimed scale.
- **FG-007**: The candidate carrier SHALL be a definite physical target on every scored trial, resolved as a registered operational state or state-equivalence class within prospective instrument uncertainty; an ensemble-average response without event-level target states SHALL NOT satisfy this clause.
- **FG-008**: The start and cessation of the formation interaction SHALL be detected by a prospectively registered physical boundary assay, and no formation input or causally equivalent surrogate SHALL remain active during the common hold.
- **FG-009**: FORMATION and SHAM targets SHALL undergo the same declared hold or storage environment after interaction cessation, and retained information SHALL be scored at the prospectively declared mission time tau_m rather than only while or immediately after the interaction is active.
- **FG-010**: The registered read Q SHALL recover information about X from the target without replaying, recreating, or supplying the formation interaction, and its measured disturbance and error SHALL remain within prospectively declared instrument-specific bounds; a destructive read is admissible if it does not recreate formation and the pre-read target state remains identifiable.
- **FG-011**: Relative to both BEFORE and SHAM under the frozen analysis, FORMATION SHALL cause a calibration-resolved increase in event-level accessible information about X that survives to tau_m, with the registered uncertainty interval clearing the prospective information floor and every named causal falsifier remaining untriggered.
- **FG-012**: A positive target directly coupled to X SHALL be classified PRIMARY_FORMATION_CONFIRMED, while a positive target written through an identified parent record SHALL be classified DERIVATIVE_FORMATION_CONFIRMED and SHALL carry an unbroken machine-readable lineage edge to that parent.
- **FG-013**: Primary versus derivative provenance and engineered versus autonomous interaction origin SHALL be reported as independent axes, and ambiguity between primary and derivative provenance SHALL block both positive classifications.
- **FG-014**: Calibration, raw-data custody, missingness, exclusions, multiplicity, stopping, randomization checks, contamination controls, read back-action, and all surface-specific uncertainty floors SHALL follow their prospective rules without outcome-conditioned repair.
- **FG-015**: At least one independent external team SHALL repeat the complete physical protocol on its own actual specimens or instances with its own apparatus and newly acquired raw data, without private project assistance, copied outcome data, or custom science logic.
- **FG-016**: The originating and independent raw bundles SHALL each run through the same immutable public URM formation contract without source modification or surface-specific science code; only validated declarative instrument mappings admitted by that release are allowed.
- **FG-017**: PRIMARY_FORMATION_CONFIRMED and DERIVATIVE_FORMATION_CONFIRMED SHALL remain unavailable unless FG-001 through FG-016 and FG-019 through FG-026 all pass in both the originating experiment and the independent reproduction wherever the clause applies.
- **FG-018**: Fragment redundancy SHALL be measured and reported when available as record strength, reach, or robustness, but SHALL NOT be treated as a constitutive requirement and SHALL NOT rescue failure of any mandatory formation clause.
- **FG-019**: Primary formation, derivative copying, autonomous or endogenous formation, engineered or exogenous formation, writer-off persistence, read-only or preexisting information, rewriting, erasure, transient response, and no change SHALL be represented as typed causal paths through one K_theta and M_phi with one event-provenance graph, not as separately scoped theories.
- **FG-020**: K_theta, its writer-off hold dynamics, and M_phi SHALL be derived or predicted prospectively from declared measurable physical primitives and SHALL predict held-out event-level outcome distributions and mode-specific falsifiers; retrospective trajectory labels or outcome-fitted science branches SHALL NOT satisfy this clause.
- **FG-021**: R SHALL admit discrete, continuous, count-valued, field-valued, distributed, or mixed physical representations, and definite state SHALL mean an actual per-specimen state or trajectory resolved within prospective uncertainty rather than an assumed binary bit or point microstate.
- **FG-022**: Every claimed active or inactive mode SHALL prospectively derive a measurable and falsifiable causal transduction C_psi and surface closure gate G_t from physical primitives, binding X through C_psi into K_theta under G_t and then setting the formation input OFF for post-closure hold dynamics.
- **FG-023**: READ_ONLY SHALL have a frozen three-way predicate that confirms recovery of preexisting target information by M_phi without causal information creation or formation replay, refutes that branch under adequate sensitivity when recovery fails or Q creates or supplies formation, and otherwise returns READ_ONLY_UNSCOREABLE.
- **FG-024**: PERSISTENCE SHALL have a frozen three-way predicate that confirms survival of closure-time target information under K_theta with C_psi OFF through tau_m, refutes the registered retention prediction under adequate sensitivity, and otherwise returns PERSISTENCE_UNSCOREABLE.
- **FG-025**: REWRITE SHALL have a frozen three-way predicate that confirms a causal K_theta transition from prospectively identified old content or lineage to prospectively identified new content or lineage that persists through its declared post-rewrite hold, refutes that prediction under adequate sensitivity, and otherwise returns REWRITE_UNSCOREABLE.
- **FG-026**: ERASURE SHALL have a frozen three-way predicate that confirms causal loss of prospectively identified accessible target information under K_theta through its declared post-erasure hold, refutes that prediction under adequate sensitivity, and otherwise returns ERASURE_UNSCOREABLE; read failure or hidden information alone SHALL NOT confirm erasure.

## 4. Exact positive predicate

For a closed claim `Omega`, let `g_k(Omega,r)` be the Boolean score of clause `FG-k` at reproduction site `r`, using only the frozen rule and that site's actual raw bundle. Let `r=0` be the originating site and let `R_ext` contain at least one eligible independent external site. The default-refuted positive predicate is

`FORMATION(Omega) = [AND over k in {1..16,19..26} of g_k(Omega,0)] AND [exists r in R_ext: AND over every applicable k in {1..16,19..26} of g_k(Omega,r)] AND [lineage_class(Omega) is unambiguous]`.

`FG-017` enforces the experimental and reproduction conjunction, `FG-018` is the redundancy guardrail, `FG-019` through `FG-022` enforce the single predictive, representation-independent theory with physical coupling and closure terms, and `FG-023` through `FG-026` close the read-only, persistence, rewrite, and erasure branches. A missing, indeterminate, unrun, non-independent, or integrity-failed clause is not a pass. Statistical or calibration uncertainty is propagated according to the frozen surface-specific rule; it is never resolved in favor of formation by default.

If `FORMATION(C)` is true, the provenance axis selects exactly one positive outcome:

- `PRIMARY_FORMATION_CONFIRMED` if the bounded interaction was directly causally coupled to `X` and did not obtain `X` through an already formed parent record.
- `DERIVATIVE_FORMATION_CONFIRMED` if the target was previously uninformative and an identified parent record causally transferred information about `X` into it; the lineage certificate is mandatory.

## 5. Exhaustive outcome map and precedence

The scorer applies the following precedence so that one closed claim receives exactly one outcome. Run-level findings may be retained as evidence, but neither positive label is emitted until the independent-reproduction clauses close.

1. **`UNSCOREABLE`** — required prospective fields, actual-data identity, integrity, sensitivity, intervention validity, uncertainty closure, public-URM execution, or eligible independent reproduction are absent or invalid, so the constitutive predicate cannot be decided.
2. **`MIXTURE_OR_UNRESOLVED`** — valid evidence reveals an unregistered mixture, incompatible site results, unresolved causal alternatives, or ambiguous primary-versus-derivative provenance that the frozen rule cannot separate.
3. **`READ_OR_PREEXISTING`** — the declared target was already informative about `X` before the alleged formation interval, or the alleged interaction only read, selected, exposed, or relabeled existing target information without causally creating previously absent information there.
4. **`TRANSIENT_RESPONSE`** — a causal target response or information increase exists during or immediately after the formation interval but fails the frozen retention decision at `tau_m` under common hold.
5. **`NO_RECORD`** — with adequate registered sensitivity, the formation contrast creates no calibration-resolved event-level information about `X` in the definite target state, even if a non-informative physical response occurs.
6. **`PRIMARY_FORMATION_CONFIRMED`** — every applicable clause passes at the originating and eligible independent sites and the unambiguous lineage is direct from `X` to the target.
7. **`DERIVATIVE_FORMATION_CONFIRMED`** — every applicable clause passes at the originating and eligible independent sites, the target lacked prior information, and the unambiguous lineage identifies a causal parent-record-to-target transfer.

Copying is therefore not grouped with preexistence. A genuine causal copy is derivative formation at the target. Conversely, lineage preservation prevents a derivative record from masquerading as a primary one.

## 6. Mode-specific predicates under the same law

These are companion judgments, not separately scoped theories. Each uses the same frozen `C_psi`, `G_t`, `K_theta`, `M_phi`, measurable primitives, prospective uncertainty system, actual-data rule, and independent-reproduction rule. For each mode, missing integrity, calibration, causal contrast, physical boundary, adequate sensitivity, public-URM execution, or eligible reproduction returns that mode's `UNSCOREABLE` outcome before either other outcome is available.

### Read-only

Let `I_pre` be accessible information in the target immediately before `Q`, let `I_Q` be the information recovered by `Q`, and let the frozen causal assays test whether `Q` injects `X`, replays `C_psi`, or drives a formation transition.

- `READ_ONLY_CONFIRMED` iff `I_pre` is present above its prospective floor, `M_phi` predicts `I_Q` within uncertainty, `Q` recovers that information without a calibration-resolved causal increase in target information, no formation input is supplied, and read disturbance remains within its prospective bound.
- `READ_ONLY_REFUTED` iff sensitivity is adequate and either the registered recovery fails or `Q` measurably creates, supplies, or replays formation rather than only reading the existing target.
- `READ_ONLY_UNSCOREABLE` iff the pre-read state, `Q` coupling, disturbance, sensitivity, public-URM execution, or independent reproduction is unresolved or invalid.

### Persistence

Let `I_close` be target information at the measured `G_t` closure boundary and `I_tau` the information after common hold to `tau_m` with `C_psi=OFF`.

- `PERSISTENCE_CONFIRMED` iff `I_close` is present, the physical assays certify closure and `C_psi=OFF`, and `I_tau` satisfies the frozen `K_theta` retention prediction within prospective uncertainty at the declared mission time.
- `PERSISTENCE_REFUTED` iff those causal and boundary assays are valid, sensitivity is adequate, and `I_tau` violates the frozen retention prediction or falls on its prospectively declared non-retained side.
- `PERSISTENCE_UNSCOREABLE` iff closure, input-off status, common hold, sensitivity, public-URM execution, or independent reproduction is unresolved or invalid.

### Rewrite

Let the target initially carry prospectively identified content or lineage `L_old`, and let the declared rewrite coupling carry `L_new` through `C_psi` under `G_t`.

- `REWRITE_CONFIRMED` iff the target's causal transition from `L_old` to the prospectively predicted `L_new` is resolved relative to SHAM, matches `K_theta`, preserves the new provenance edge, and persists through the declared post-rewrite hold with formation input off.
- `REWRITE_REFUTED` iff causal controls and sensitivity are adequate and the registered transition or post-rewrite retention prediction fails.
- `REWRITE_UNSCOREABLE` iff old/new content, lineage, coupling, closure, hold, sensitivity, public-URM execution, or independent reproduction is unresolved or invalid.

### Erasure

Let `I_old` be prospectively identified accessible target information before the declared erase coupling and `I_post` its value after the declared post-erasure hold.

- `ERASURE_CONFIRMED` iff the erase coupling causes the prospectively predicted loss of `I_old` relative to SHAM under `K_theta`, the loss survives the post-erasure hold, and independent alternate reads and carrier assays exclude mere read failure or reversible hiding within their prospective sensitivity.
- `ERASURE_REFUTED` iff causal controls and sensitivity are adequate and the target retains recoverable `I_old` contrary to the frozen erasure prediction.
- `ERASURE_UNSCOREABLE` iff loss cannot be separated from instrument failure, inaccessible-but-retained information, closure or hold failure, inadequate sensitivity, public-URM failure, or missing independent reproduction.

## 7. Required evidence bundle

Each originating or reproduction bundle must be closed, content-addressed, and include:

1. claim tuple and frozen protocol instance conforming to `FORMATION_PROTOCOL_V001.schema.json`;
2. physical surface, specimen or instance, apparatus, operator, location, time, environment, and calibration identities;
3. raw BEFORE, interaction-boundary, FORMATION/SHAM, COMMON-HOLD, mission-time READ, read-back-action, and calibration observations;
4. event labels and timestamps held independently from target readout until the declared blind is opened;
5. for `RANDOMIZED_INTERVENTION`, a randomization commitment, assignment log, block structure, attrition, exclusions, and deviations; for `PROSPECTIVE_NATURAL_CONTRAST`, independently measured forcing, the frozen matching rule and negative time or surface controls, response-blind access records, alternative-cause falsifiers, attrition, exclusions, and deviations;
6. writer-off or interaction-cessation certificate and common-hold equivalence assay;
7. prospective uncertainty budgets, detection and absence floors, power or sensitivity justification, and all decision thresholds with units;
8. event-level target states, readouts, frozen information statistic, uncertainty intervals, null and sham contrasts, and falsifier vector;
9. read-disturbance assay and evidence that READ did not recreate formation;
10. provenance certificate, including source-event linkage for primary claims or parent-record lineage for derivative claims;
11. fragment-redundancy measurements when acquired, marked non-constitutive;
12. complete raw-to-URM transformation manifest, immutable URM release hash, machine result, logs, and environment lock;
13. custody hashes for every raw and derived object and an explicit inventory of missing or excluded objects; and
14. for an independent site, signed declarations of independence, independently measured specimens, instances, or archives, its own apparatus, newly measured raw data, public-only assistance, and absence of custom science logic; a unique natural event need not be recreated when independent material carrying that event is remeasured.
15. frozen `C_psi`, `G_t`, `K_theta`, and `M_phi` specifications, measurable-primitive inputs, event-provenance graph, held-out predictive distributions, and mode-specific falsifiers for formation, copying, persistence, read-only behavior, rewriting, erasure, transient response, and no change.

## 8. Actual-data and reproduction rules

“Actual data” means measurements made on the identified physical target by a real instrument under the registered protocol or retained in a response-blind archival holdout. Data need not be newly generated: an archival bundle is admissible only when the protocol, analysis, uncertainty system, controls, and immutable URM release were frozen before the scoring team received numerical access or unblinded response labels, and custody proves that blindness. A deposit whose outcome-bearing responses were already numerically accessed by that team can motivate terms or test ingestion but cannot satisfy this prospective gate. Simulations and synthetic controls can validate software but carry zero weight on the physical predicate.

The independent team must be capable of running the public instructions without private coaching from the program and must acquire a genuinely independent physical measurement bundle. For a repeatable event this may recreate the event. For a unique natural event it may instead newly measure independent specimens or archives that carry the same event; reusing the originating raw measurements is not reproduction. It may perform ordinary apparatus setup and instantiate the contract's declared units, calibration, uncertainty floors, and instrument mapping. It may not change the constitutive predicate, analysis logic, outcome map, or URM source. Its mission time must equal the frozen claim mission time; its apparatus-specific uncertainty floors must be independently fixed from blind calibration and must possess adequate sensitivity for the claim.

If the official public URM cannot accept the independent bundle through its declared input contract, the surface claim is `UNSCOREABLE` and the general any-surface ambition remains open. Adding a new constitutive term or surface-specific science branch after seeing the result does not repair that run; it refutes or reopens the broader theory and requires a new prospective version.

Rewriting and erasure are not relabeled formation. They are transitions of an already record-bearing target under the same `K_theta`: rewriting changes the prospectively declared content or provenance, while erasure causally removes accessible information under a declared post-intervention hold. Read-only access is the `Q` branch of `M_phi` with no new target information. A surface on which a branch is physically unavailable must receive a prospective negative prediction from the same primitives; the branch may not simply be omitted after observation.

## 9. Relation to the program's completion condition

This gate can establish reproducible formation only for its closed claim scope. The program completes only when real physicists can use the released URM to reproduce the record-formation and gravity-emergence results on any and every bona-fide real record surface under a constructive universal coverage argument. This document supplies the formation predicate; it does not assert that completion condition has been met.

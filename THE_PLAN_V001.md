# THE PLAN — to a PROVEN process

**Call it with `./ledger/plan.py`. `ledger/plan.tsv` is the source of truth; this file is
generated. Tasks append, IDs never renumber, and every task carries a DONE_WHEN so "done" is
checkable rather than a judgement call.**

## THE DEFINITION OF DONE

> **A full demonstration of a PROVEN process that explains how EM, gravity, and alpha form
> THE WORLD'S RECORDS — with EVERY PROOF GROUNDED IN AT LEAST ONE NAMED PHYSICAL RECORD
> WITH REAL NUMBERS — and the model, so that the proof can be checked by anyone AGAINST
> THOSE RECORDS.**
>
> — the principal, 2026-08-19; grounding clause added 2026-08-20 after ~20 real records
> were tested against the five clauses and ZERO satisfied them (H-3).

> **PROOF means: a test ANY PHYSICIST ANYWHERE can run against THEIR OWN real-world data
> and confirm works as asserted. Otherwise we are not proving anything other than that we
> can construct equations.** — the principal, 2026-08-20

**Nothing below counts as finished until every clause of that is true. The grounding clause
converts H-3, H-5 and H-6 from open hazards into standing acceptance criteria.**

**22 of 33 tasks done.**  ▶ marks the CRITICAL PATH — the chain the claim rests on.

| DONE | DOING | TODO | BLOCKED |
|---|---|---|---|
| 22 | 1 | 6 | 4 |

---

## PHASE A  MAKE IT LEGIBLE  —  3/3 done

| # | task | status | DONE WHEN | depends |
|---|---|:---:|---|---|
| ▶ **T-1** | **THE GLOSSARY. Map every program term to its standard name and owner — record = logical operator, writer = conjugate logical, channel = OURS coined 2026-08-19, Gamma = placeholder for (H_1,<,>). Mark each term BORROWED or OURS** | **DONE** | GLOSSARY.md exists, every term in the anchor and the live claim documents appears in it, and each is marked BORROWED (with owner) or OURS | — |
| ▶ **T-3** | **FIX O-19. record_model.py builds the *-algebra when the commutant needs only the generators. This is why it caps out at dim 32** | **DONE** | the model computes the commutant from generators, validate_model.py still passes 12/12, and a dim-256 carrier completes | — |
| **T-2** | ONE CURRENT STATEMENT OF THE PROCESS. The register is 6000+ lines and append-only; THE_CLAIM predates most of the work. Write the process as it now stands, in mapped terms | **DONE** | PROCESS_V002.md states the process end to end, every sentence cites a ledger row, and no withdrawn claim appears | T-1 |

## PHASE B  MAKE THE MODEL THE PROCESS MODEL  —  4/4 done

| # | task | status | DONE WHEN | depends |
|---|---|:---:|---|---|
| ▶ **T-4** | **ADD THE ENVIRONMENT. The model has no bath, no coupling, no chi. Add them so formation is computable from (H, {L_k}, coupling)** | **DONE** | model.formation(coupling) returns chi and reproduces F-20's 0.97527192 and O-18's three-row table | T-3 |
| ▶ **T-5** | **ADD THE CHANNEL CRITERION. G-16 — chi > 0 iff the coupling pairs oddly with the record's conjugate — is nowhere in the model or MODEL.md** | **DONE** | model.channel(coupling, record) returns the pairing and predicts chi>0 correctly on all cycles of the 2x2 torus | T-4 |
| ▶ **T-6** | **MULTI-RECORD AND DEPENDENCIES. Extend independence() to report which records share a channel and which formation events interact** | **DONE** | the model returns, for a family, which members can be formed without disturbing the others — measured, not asserted | T-5 |
| **T-7** | VALIDATE THE FORMATION HALF. Every formation result in the register re-derived through the model rather than a one-off lane script | **DONE** | validate_model.py covers F-20, F-21, F-23, G-16, C-17, C-18 and passes | T-6 |

## PHASE C  BREAK THE SINGLE CARRIER  —  4/5 done

| # | task | status | DONE WHEN | depends |
|---|---|:---:|---|---|
| ▶ **T-8** | **CHOOSE THE SECOND CARRIER FAMILY. PF-3B failed because the bouquet has cosystole 1. Find a non-manifold with BOTH distances >= 2, or state why none is reachable** | **DONE** | a named second carrier with both distances >= 2 is built and sealed, or O-14 is closed with a proof that none exists at reachable dimension | T-3 |
| ▶ **T-13** | **O-15. The bath is 3 qubits and chi(t) recurs. Show the arrow survives a bath large enough not to recur, and that redundancy evens out** | **DONE** | chi(t) monotone to a plateau on a larger bath, with fragment redundancy measured | T-3 |
| ▶ **T-9** | **RE-RUN EVERY LOAD-BEARING ROW ON TWO CARRIERS -- REOPENED 2026-08-20 BY EXTERNAL REVIEW. It was marked DONE while satisfying NEITHER branch of its own DONE_WHEN: the string SINGLE-CARRIER appears ZERO times in ledger/status_ledger.tsv, and only about 24 of 162 rows cite two structurally different carriers. REPLICATE.md's claim that single-carrier rows are marked was FALSE of the shipped ledger and has been corrected in place** | TODO | every PROVED row either cites two structurally different carriers or is explicitly marked SINGLE-CARRIER in the ledger | T-8 |
| **T-11** | PF-3B. Clause (v) under DEF-A on the second carrier — does admissible still give zero flippers? | **DONE** | the second carrier reproduces 0 admissible flippers against a nonzero any-unitary control | T-8 |
| **T-10** | PARAMETER INDEPENDENCE. Show no conclusion moves with beta, lambda, bath size or lattice size | **DONE** | a sweep table per PROVED row showing the conclusion invariant, with the noise floor printed beside every fit | T-9 |

## PHASE D  CLOSE THE OPEN PHYSICS  —  3/3 done

| # | task | status | DONE WHEN | depends |
|---|---|:---:|---|---|
| **T-12** | PF-4. The necessary condition. C-4 is PARTIAL: a sufficient family and a proven non-CSS escape, no established necessary condition | **DONE** | either a condition covering both is stated and tested, or it is proved none narrower than the clauses exists | — |
| **T-14** | ALPHA'S ROLE, RESTATED AFTER PF-6. Theorem D's eps^d was rescoped to GENERIC single-site perturbations; the Z-only case gives eps^n*. Restate alpha's role in mapped terms | **DONE** | the alpha row states which perturbation class it holds for, with both exponents measured | T-1 |
| **T-22** | GIVE ALPHA REAL CONTENT (H-2). Nothing in the corpus is dimensionful and the fine structure constant is connected to no number. Either connect it or state plainly that 'alpha' names a generic coupling here | **DONE** | either a registered result in which alpha's VALUE matters, or an explicit statement in the anchor that alpha stands for a generic coupling strength | — |

## PHASE E  THE CLAIM: DOES GRAVITY EMERGE?  —  8/13 done

| # | task | status | DONE WHEN | depends |
|---|---|:---:|---|---|
| ▶ **T-20** | **DOES CLASSICAL GRAVITY EMERGE FROM RECORDS? Build a COLLECTIVE quantity from MANY records and ask whether it obeys anything Einstein-constraint-like — without importing the constraint (D-1 binds: the classical form may be recovered, never assumed or used as the test)** | **DONE** | a many-record collective quantity is defined, computed, and either shown to satisfy a constraint of the right type or shown not to — with the classical form never used as a criterion | T-6, T-9, T-13 |
| **T-21** | REVISIT X-2 AT THE COLLECTIVE LEVEL. The record-level reclassification stands; the emergence question it foreclosed does not | TODO | X-2 carries an explicit note distinguishing the record-level question (category error, closed) from the collective one (open, = T-20) | T-20 |
| ▶ **T-23** | **IS THE CORRELATION ENERGY A SOURCE? O-47 exhibited an interaction energy between records that coexists with FREE single-record writes and costs no clause -- the first structure with a source's shape not paid for by breaking something. Test it against all five standards** | **DONE** | the correlation energy is measured against extensive, additive, not-a-count, sign-definite and power-law-falloff; every verdict marked EXACT or TREND; and every effect labelled INSERTED or INDUCED | T-20 |
| ▶ **T-24** | **CLAUSE (v) ON THE CARRIERS WHERE THE ENERGY LIVES. Protection was never tested in O-42, O-44 or O-47, nor anywhere in O-36 on the non-abelian carrier. Without it these are records in four clauses, not five** | TODO | protection is checked on the correlation carriers and on D(D_4), with the contractible-region size swept, or it is stated plainly that these are four-clause records and every claim resting on them is scoped that way | T-23 |
| **T-25** | SETTLE WHAT ADMISSIBLE MEANS. DEF-A -- a unitary with [U,H] = 0 -- is marked PROVISIONAL on one carrier family, its physical-channel disjunct is UNTESTED, and clause (v) is stated with the same word so protection moves with it. O-44 showed it was briefly the sole support of the program's strongest negative result | TODO | either a definition is adopted and clauses (iv) and (v) are re-run under it on two carriers, or the framework states plainly that DEF-A is a CHOICE and names every row that turns on it | — |
| **T-26** | RE-AUDIT EVERY REGISTERED NULL AGAINST THE VENUE RULE. Four measurements in one session were taken where the effect could not appear -- abelian carriers for a frame rotation, one-qubit bath sites for a capacity question, gauge-invariant states for transport, a non-degenerate H for a record. D-22 says a permutation-symmetric carrier has no geometry to detect at all | TODO | every registered null either cites the carrier's demonstrated capacity to show the effect -- a positive control in the same table -- or is rescoped to absence-within-range or absence-in-this-venue | — |
| **T-27** | BRING PROCESS_V002 AND THE ANCHOR CURRENT. One session withdrew C-37 and C-55, rescoped C-29, C-30, C-32, C-33, C-36 and C-43, failed O-42, and reopened then reclosed O-4. The process document and the anchor sentence predate all of it | TODO | PROCESS_V002.md and the anchor cite no WITHDRAWN or FAILED row, every sentence traces to a live ledger row, and the anchor's claim about gravity matches what the ledger now says | T-23 |
| ▶ **T-28** | **IMPLEMENT THE AMENDED DEFINITION IN THE MODEL. records(t_m) = slow Lindblad eigenmodes with \|Re lambda\| <= 1/t_m, and ADMISSIBLE(beta, W) as a smallest-bath-plus-one-ladder-site search in O-44's style. This is the one-move model change the census names** | **DONE** | the model returns records under (i')-(v') with declared tolerances, and REPRODUCES the exact records on the three registered carriers in the t_m -> infinity, W -> 0 limit | — |
| ▶ **T-29** | **THE FIRST PROOF GROUNDED IN A WORLD RECORD: the HDD CoCrPt grain as an (H,{L_k}) two-state macrospin with Arrhenius jumps -- K_uV = 61 kT, tau = 2e17 s, splitting 0.007-7 kT, write head 5e-19 J. WARNING from the census: C-14's count law collapses to zero records under ANY well asymmetry, so the width-clustered quasi-degenerate multiplicity fix must be adopted FIRST or the run provably returns zero** | **DONE** | the model takes the grain, returns its record under the amended clauses, and yields a number a physicist can check against published retention data | T-28 |
| ▶ **T-30** | **DOES THE GRAVITY-STRENGTH EXCLUSION SURVIVE AT FINITE TOLERANCE? Re-run C-52's trace-ratio computation on a METASTABLE record at declared tolerance. The Z[i] quantisation was proved over exact signed Paulis at tolerance zero; exclusions proved at tolerance zero do not automatically survive at finite tolerance. This single computation decides whether the exclusion was a theorem about the WORLD or about the corner** | **DONE** | C-52 and C-53 restated at finite tolerance on a metastable carrier, with the verdict stated either way | T-28 |
| **T-31** | RE-RUN THE DEGENERACY-DEPENDENT GRAVITY CHAIN WITH ASYMMETRIC WELLS. C-14's count law k = min_E v_2(m_E) and the W-42/W-43/W-44/W-51 record-to-geometry chain depend on engineered degeneracy; under any well asymmetry the multiplicities collapse. The chain survives only if the record COUNT is insensitive to well-depth asymmetry -- an explicit re-check, not a presumption | **DONE** | the count law and the geometry chain are re-run with asymmetric wells and either survive with a stated tolerance or are withdrawn | T-28 |
| ▶ **T-32** | **IS THE CONFIGURATION ENERGY OF REAL RECORDS A SOURCE? Test the five standards on written magnetic media rather than on the toy chain, with a demagnetised control** | **DONE** | the five standards are scored on a real record with a control that would have failed them, and the verdict stated either way | T-30 |
| ▶ **T-33** | **DOES THE CALCULATION AGREE WITH ANY RECORD SURFACE? Verify the record laws across structurally different mechanisms, with surfaces where they must NOT apply as controls** | **DONE** | two laws verified on at least four mechanisms at machine precision, with at least one control where the model declines | T-29 |

## PHASE F  WRITE AND CHECK THE PROOF  —  0/4 done

| # | task | status | DONE WHEN | depends |
|---|---|:---:|---|---|
| ▶ **T-15** | **THE THREE ROLES AS ONE STATEMENT. EM supplies the complex and the record; Gamma supplies space, writer, protection and channel; alpha sets the cost** | **BLOCKED** | A-PR restated with every clause citing a two-carrier row, and no clause resting on a WITHDRAWN or PARTIAL row | T-9, T-14, T-20, T-23 |
| ▶ **T-16** | **WRITE THE PROOF. End to end: what is claimed, the theorems, the evidence, and what each rests on** | **BLOCKED** | PROOF_V002.md, every step citing a ledger row, no step resting on a single carrier unless marked | T-15, T-23 |
| **T-17** | STATE THE LIMITS IN THE PROOF ITSELF. PF-6 failed; X-4, T-VI.3, T-VI.4 are BLOCKED; empirical contact is zero | **BLOCKED** | the limits section names every BLOCKED row and says plainly that this is a mathematical result with no empirical contact | T-16 |
| ▶ **T-18** | **EXTERNAL CHECK. Someone outside the program reproduces the result from REPLICATE.md alone** | **BLOCKED** | reproduce.sh passes end to end on a clean clone, and every headline number in PROOF_V002 is traceable to a sealed lane output | T-16 |

## PHASE G  CONTINUOUS  —  0/1 done

| # | task | status | DONE WHEN | depends |
|---|---|:---:|---|---|
| **T-19** | THE GATE STAYS GREEN. on_finding.sh passes after every finding, for the life of the plan | **DOING** | the gate has passed on every commit since it was built | — |

---

## STATUS VOCABULARY

| | |
|---|---|
| `DONE` | the completion criterion in DONE_WHEN is met and verified |
| `DOING` | in progress |
| `TODO` | not started, not blocked |
| `BLOCKED` | cannot start until its dependency clears |
| `DROPPED` | deliberately abandoned; the reason is in the register |

---

> **The plan is not the goal. The program is. If a task turns out to be the wrong task,
> it is DROPPED with the reason in the register, and a better one is appended.**

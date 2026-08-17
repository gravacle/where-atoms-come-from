# W-17 — RULING ON THE FRAME CHALLENGE INSTRUMENT (DRAFT, for the principal's seal)

**Subject:** `FRAME_CHALLENGE_V001.md` (81 lines, 2026-08-17, PROPOSED NOT ADOPTED), validated by
retrodiction over four route decisions R1-R4 under reading cutoffs, plus a cutoff audit and an
instrument audit.

> ## VERDICT
> **THE INSTRUMENT WORKS, AT TWO FIFTHS OF ITS STATED SIZE. ADOPT `F1` (AMENDED) AND `F5`. DELETE
> `F2` AND `F3`. DEMOTE `F4` TO A GREP GATE. DO NOT ADOPT THE FIVE-TEST VERSION AS WRITTEN — IT
> CONVICTS THE PARTITION OF THE INTEGERS INTO EVENS AND ODDS, AND IT CONVICTS THE ONE BINARY W-15
> ITSELF CERTIFIES AS FORCED WITH A PROOF UNDER IT.**
>
> One certified retrodiction: **`F5` at R1**, and it is strengthened here to stand without the
> carrier the cutoff audit objected to. One hard **SOUND** (`F2` at R2) plus one in-scope sound
> steelman (`F1` at R4). The 19-of-20 fire rate reported by the instrument audit is **not** evidence
> the routes were broken, and this ruling settles that by running the control the round omitted.

---

## 0. LINEAGE AND STANDING — DECLARED FIRST, BECAUSE IT DISCOUNTS EVERYTHING BELOW

I am **Claude Opus 5**. So are all four frame lanes, the cutoff auditor, the instrument auditor, and
every register row from **W-07** onward (`REGISTER_V001.md:1112`, `:1245`). This ruling is
**adversarial checking inside one lineage**, which is precisely the grade `CUSTODY_V001.md:44-53`
assigns and precisely the grade W-03's still-unfired reopen conditions (`REGISTER_V001.md:87`,
`:158`) exist to lift. **This ruling is layer fourteen of one block. Discount it as one block.**

My own working files: `LANE_W17_RULING/` — `sound_frame_control.py`, `sound_frame_control.OUT.txt`,
`CONFOUNDS.txt`, sealed in `SEALS.sha256` (3/3 verify). All four frame-lane seals verify:
**11/11, 8/8, 7/7, 2/2.**

---

## 1. THE CONTROL THE ROUND DID NOT RUN — AND IT CHANGES THE READING OF THE WHOLE ROUND

The instrument audit named the missing experiment and did not run it. **It is run here.** The
question is not "did the tests fire" but "do they fire on frames that are *sound*". Two controls:

* **CONTROL 1 — PARITY.** Arms: `n` even / `n` odd. Route form: *take the even branch or the odd
  branch?* An unimpeachable partition of `Z`. Anything that fires here is a fault machine.
* **CONTROL 2 — `diagonal / non-diagonal` on `U(3)`.** **W-15's own named paradigm** of a binary the
  mathematics handed over: *"`diagonal / non-diagonal` was forced and has a proof under it"*
  (`REGISTER_V001.md:1509-1511`).

### RESULTS (`LANE_W17_RULING/sound_frame_control.OUT.txt`, sealed `a25e0070…`)

| frame | F1 (letter) | F2 | F3 | F4 | F5 | fired |
|---|---|---|---|---|---|---|
| **PARITY** | **FIRE** | **FIRE** | ok | ok | ok | **2 / 5** |
| **DIAGONAL / NON-DIAGONAL** | **FIRE** | **FIRE** | **FIRE** | *indeterminate* | ok | **3.5 / 5** |

The numbers behind each firing:

* **`F2` on PARITY = distance 0.** The instrument lists *change of coordinates* among its candidate
  maps. `f(n) = n+1` carries the odd arm **exactly onto** the even arm: symmetric difference `2` of
  `100000` (boundary term only), density `0.00002000 -> 0`. **`F2` returns "the arms are the same
  object under a map" for evens versus odds.**
* **`F2` on `diagonal / non-diagonal` = distance 0.000e+00, on a positive-dimensional set.** For the
  3-cycle `P` and **any** diagonal unitary `D`, `(P.D)^3` is diagonal. Over `20000` draws: minimum
  off-diagonal norm of `U = P.D` itself `1.732051` (solidly in arm B), maximum off-diagonal norm of
  `U^3` **`0.000e+00`**. A 3-real-parameter slab of arm B lands inside arm A under a listed map.
  **This is the corpus's own instance**: `||T_F^3 - M_dF|| = 6.18e-16`, `||T_C^3 - M_c|| = 4.34e-16`
  (`REGISTER_V001.md:1519-1534`).
* **`F3` on `diagonal / non-diagonal` = measure-zero cell.** Haar on `U(3)`, `20000` draws:
  `mu(diagonal) = 0.000000` (0 draws below `1e-9`, minimum observed off-diagonal norm
  `2.466693e-01`), `mu(non-diagonal) = 1.000000`. `F3`'s stated rule — *"an empty cell, **a
  measure-zero cell**, or an overlap means the predicate does not partition"* — **fires on a
  predicate that partitions `U(3)` exactly, overlap 0, union 1.**
* **`F1` on both = an embedding always exists.** See §5.
* **`F4` on `diagonal / non-diagonal` = indeterminate, by construction.** Sentence 1, *"a basis is
  designated"*, is true in-corpus (S1's vertex classes). Sentence 2, *"diagonality is a property of
  the operator alone"*, is false — `2000` of `2000` Haar unitaries are diagonal in their own
  eigenbasis, so under a free basis **arm B is empty**. Two defensible sentences, opposite verdicts.
  **`F4`'s output is chosen by the lane, not forced by the frame.**

**Now the same two frames, under the amended instrument recommended in §4** (`F1` with the exclusion
clause; `F2` deleted; `F3` demoted to a reporting line; `F4` demoted to a grep gate; `F5` unchanged):

| frame | F1 (amended) | F5 | fired |
|---|---|---|---|
| PARITY | ok | ok | **0 / 2** |
| DIAGONAL / NON-DIAGONAL | ok | ok | **0 / 2** |

**RULING ON THE FIRE RATE.** 19-of-20 across R1-R4 was reported as though it were evidence about the
routes. It is not, on its own: the same instrument fires 2-to-3.5 times out of five on frames nobody
disputes. **What survives the control is exactly the two tests that clear both controls — and those
are the same two tests that produced 4 of the round's 5 genuinely-new findings.** That coincidence is
the strongest single result in this validation.

*Confounds on the controls, recorded not answered: neither control is a route decision
(`FRAME_CHALLENGE_V001.md:14-20` restricts the instrument to route decisions); I did not construct a
sound ROUTE decision, so the controls bound the false-positive rate on sound partitions, not on
sound routes. Full ledger `LANE_W17_RULING/CONFOUNDS.txt` C1-C7.*

---

## 2. DID IT CATCH A FRAME FAILURE FROM BEHIND THE CUTOFF?

### R1 — **YES. ONE CATCH, CERTIFIED, AND IT IS `F5`.**

**The test that fired:** `F5`, THE NULL OPTION, on the sentence at
`FOUNDING_DESIGN_V001.md:63-65` — *"An inductive limit of finite objects is not finite — which is
precisely how it escapes recurrence."* That sentence is the founding design's **single stated escape**
from its own §4 obstruction, and `FOUNDING_DESIGN_V001.md:99-102` sends stage S3 out *"with the
inductive-limit template as the starting point."* **It is the sentence that licensed the program's
central stage.**

**The number:** the escape inference is invalid, because non-finiteness is **necessary** for
non-recurrence and not **sufficient**. R1's lane measured this on K1's electric spectrum — perfect
revival `|A(2*pi)| = 1.000000000000000` at truncations `M = 1, 5, 20, 60, 120` (dimensions 9 to
58081).

**The cutoff audit objected that this is K1-dependent, and K1 postdates the founding decision. THE
OBJECTION IS ANSWERED HERE.** The inference dies on *any* commensurate infinite spectrum, with no
complex, no gauge group and no fibre (`sound_frame_control.OUT.txt`, final block):

```
    dim =       10   E_n = n   |A(t=2*pi)| = 1.000000000000000
    dim =      100   E_n = n   |A(t=2*pi)| = 1.000000000000000
    dim =    10000   E_n = n   |A(t=2*pi)| = 1.000000000000000
    dim =  1000000   E_n = n   |A(t=2*pi)| = 1.000000000000000
```

**Two lines of arithmetic, requiring nothing but the sentence being tested.** This is now
cutoff-legal in the strongest available sense: it needs no corpus at all.

**And it is new against the entire corpus.** Machine-checked over all 1696 lines of
`REGISTER_V001.md`, W-01 through W-16: occurrences of **`inductive limit`: 0.** Occurrences of
**`quasi-local`: 0.** Sixteen rows, eleven of them adversarial, **and not one examined the escape
route that licensed S3.** Checked again across the twelve sealed artifacts (S1, S2, S2-audit, S3,
S3-audit, S4, W07, W08, W10, W11, W13, PUBLICATION): **one** occurrence, at
`S2_FORMATION_CONDITION_ON_K1_AUDIT_V001.md:376-377`, and it **cites the template as a reference for
what "subsystem" means** — it does not test whether the escape follows. **Total: one mention, zero
tests, in the entire corpus.** W-08 does declare *"THE FOUNDING OBSTRUCTION IS FALSE AS AN INFERENCE"*
(`REGISTER_V001.md:769`), but for a **different defect** — the single-cell `|Z_k|` versus the product
`|Omega_N|`. R1's `F5` attacks the **escape**, not the obstruction. Neither is a restatement of the
other.

**Was it leaked?** No. The instrument's two worked examples cover `F2` (`:41-44`) and `F4`
(`:54-56`). **`F5` has no worked example anywhere in the instrument** — §1's `F5` entry is four lines
of bare test statement. The session auto-memory leak carried *"forced crossing does not exist"*,
which is the program's terminal conclusion and not the invalidity of the escape inference. **This
catch is the only one in the round that no leak vector reaches.**

**SCORED: the instrument could have caught it at the time, `F5` is the test that fires, and the
number is `1.000000000000000`.**

### R3 — **NO. NOT A CERTIFIED CATCH, AND I SCORE IT AS A MISS RATHER THAN LEAVING IT OPEN.**

R3's central finding — *"'transport convention' misnames the operative variable, which is the
clock"* — is within one refinement of the **ERRATUM AGAINST W-11** (`REGISTER_V001.md:1519`),
which names the **advance lattice** `L_F Z x L_C Z`. But:

1. **The direction was handed to R3 by its own grep window.** The cutoff audit established that
   `grep -n "^#" REGISTER_V001.md` returned line 1519 — *"ERRATUM AGAINST W-11 — **ITS HEADLINE
   NAMES THE WRONG OBJECT**, AND TWO LANES SAID SO AND WERE OVERRIDDEN"* — i.e. R3's conclusion, in
   one line, in the output R3 read. R3 disclosed seeing the headings and did **not** identify that
   this one carried its finding.
2. **`F2`'s number and its intended reading were pre-supplied** by `FRAME_CHALLENGE_V001.md:41-44`.
3. **The uncontaminated residue restates W-10.** R3's `F4` value-spread (`0.405295883` nats over five
   published rows; two rows sharing every published invariant) is W-10's own defect **N-1**.
4. **`F3`'s load-bearing leg is a method the corpus has already convicted.** *"0 of 3000 Haar draws"*
   is exhaustion over a sample, and the **ERRATUM AGAINST W-14** (`REGISTER_V001.md:1422`, `:1455-1462`)
   rules that sampling is *structurally blind to measure-zero sets, and the interesting sets in this
   problem are measure-zero*. My Control 2 makes the same point from the other side: on
   `diagonal / non-diagonal` the entire content of arm A sits at Haar measure 0.

**SCORED: R3 produced no catch this round can certify.** Note that this is a judgement about
*evidence*, not about R3's competence — its lane is careful, its seals verify, and its disclosure was
better than its self-assessment.

---

## 3. DID IT REPORT SOUND ANYWHERE? — YES, TWICE, AND IT MATTERS

**An instrument that never clears a frame is worthless. This one clears frames.**

* **HARD SOUND: `F2` at R2.** The two arms of *STOP versus STOP-FALLS-REBUILD* genuinely differ on
  `pi_build`, and the difference has checkable content: the wedge-growth arithmetic reproduces
  exactly — `k = (5^n - 1)/4 = [1, 6, 31, 156, 781]`, `V = [5, 25, 125, 625, 3125]`, every
  consecutive pair dividing, forced because `5 = 1 (mod 4)`. **Independently reproduced here.**
* **IN-SCOPE SOUND: `F1` at R4.** Inside the transport convention the connection/state binary *is* an
  exhaustive partition: `rank(state) = 3`, `rank(connection) = 2`, `rank(both) = 5 = 3 + 2`, by SVD.
  R4 ran the steelman, reported it, and then showed the partition fails only once W-11 makes the
  convention movable. **That is the instrument behaving exactly as designed.**
* **AND THE AMENDED INSTRUMENT CLEARS BOTH CONTROLS, 0 fired of 4 slots** (§1).

**SCORED: specificity is demonstrated for the amended two-test instrument. It is NOT demonstrated for
the five-test instrument as written.**

---

## 4. WHICH TESTS EARNED THEIR PLACE — RECOMMENDED DELETIONS

Classification of the round's 19 fired slots (instrument audit, spot-checked here): **8 restate
register text already on record · 6 attach a new number to an already-registered fact · 5 genuinely
new.** The five genuinely-new are **R1/F1, R1/F5, R2/F4, R2/F5, R4/F1** — i.e. **`F1` twice, `F5`
twice, `F4` once. `F2`: zero. `F3`: zero.**

### KEEP — `F5`, AND RUN IT FIRST

Grade: **strongest of the five.** It produced the round's only certified retrodiction (§2), it
produced two of the five genuinely-new findings, it cleared both controls, and it is the cheapest to
run — its deliverable is a **witness for each arm**, or a demonstration that no witness exists. It is
also the one test the instrument text says nothing prescriptive about, which is the pattern §4.5
below turns into a rule. **Move it to position one:** if neither branch obtains, the embedding
question does not arise.

### KEEP — `F1`, AMENDED (see §5 for the amendment)

Grade: **load-bearing and most abusable.** It generated the re-posed question in **4 of 4** routes
and produced the round's one register-reopening result (R4/F1-D: a non-uniform charge moves `lambda`
by `0.110749834` nats at identical `pi` and identical `(u,v)`, with a uniform-charge control at
`1.6e-07` — this fires **W-14's own REOPENS clause verbatim**, `REGISTER_V001.md:1417`). It also
fires on both controls **as written**, and must not be adopted without the exclusion clause.

### DELETE — `F2` (DEGENERACY)

Four independent reasons, any two of which would be enough:

1. **Zero new content in the round.** All three `F2` numbers were forced by the arms' own
   definitions: R1's `dist(union_n A(K_n), A_inf) = 0` **is the definition** of a C\*-inductive limit
   (a control that could not have failed, and the lane banked it instead of voiding it); R3's
   `T^3 = M_gamma` is **printed** at `S3_THE_CROSSING_AUDIT_V001.md:182` and has been in the register
   since W-04 at `REGISTER_V001.md:338`; R4's is N3 read forward.
2. **It fires on both controls, with distance exactly 0 on parity.**
3. **It cannot distinguish a symmetry from a quotient.** A map exchanging the arms is, in a sound
   binary, *evidence the partition is clean* — `n -> n+1` on parity is the archetype. A map
   *collapsing* one arm into the other is evidence of degeneracy. The corpus's real case is the
   second kind (`M_gamma` is a **power** of `T`, so the corpus's operator is a **quotient** of
   COR-F's), and `F2` as written returns the same number for both.
4. **It does not prevent its own founding failure mode.** `F2` exists because W-11 leg A *"read the
   degeneracy number backwards"*, and `F2`'s remedy is to **mandate the direction of reading**
   (*"report the number AS A FRAME RESULT"*, `:38`). Mandating a direction is the same error with the
   sign flipped, and this round shows the cost: the mandated direction was leaked to R1 and R3 and
   both lanes said their `F2` reading was not independently arrived at.

**SALVAGE — one clause, folded into `F1`:** *if a listed map carries one arm into the other, state
whether it is a **symmetry** (arm-exchanging; evidence the partition is sound) or a **quotient**
(collapsing; evidence it is not), and give the fibre dimension.* R3 in fact computed the right number
for this without needing `F2`: the fibre of `S -> S^L` over `M_gamma` has **real dimension 6**, so
`M_gamma` is a function of `T` and `T` is not a function of `M_gamma`.

### DELETE — `F3` (CARVING), AS A TEST

**Its stated rule is false.** *"An empty cell, a measure-zero cell, or an overlap means the predicate
does not partition"* — the middle clause convicts `diagonal / non-diagonal` at
`mu(diagonal) = 0.000000` over 20000 Haar draws, on a predicate that partitions `U(3)` exactly. **A
measure-zero cell means the dichotomy is non-generic, not that it fails to partition** — and in this
program the non-generic cell is repeatedly where the content is: **W-13** rules N1 *true at the
corpus's own published connection* and *false on a comeager set*, i.e. the corpus's headline lives on
a measure-zero set; the **ERRATUM AGAINST W-14** rules sampling structurally blind to exactly those
sets.

**SALVAGE, two pieces:**
* Keep *"measure both cells under a stated measure and report the two numbers"* as a **reporting
  line with no verdict attached.** Reporting is cheap and was informative (R4's three measures, and
  its finding that the answer changes with the measure, is good work).
* R2's `mu_reopen` leg — *both arms fire, discharge, or name **0 of the 11** live registered reopen
  conditions* — is the one strong product of `F3` this round, and **it is not a measure argument at
  all.** It is a relevance count against the register's own live questions. **Move it into `F1` as a
  mandatory line: for each arm, how many of the register's standing reopen conditions does it
  touch?**

### DEMOTE — `F4` (PRESUPPOSITION), TO A GREP GATE

* **4 of 4 fired; 3 of 4 restate register text already on record** — R1's *"the physical algebra does
  not factorize"* is `REGISTER_V001.md:318` verbatim; R4's `grep -ic constraint = 0` restates
  `REGISTER_V001.md:397` (*"no backreaction, no Gauss law, no constraint, no equation of motion"*).
  I reproduced R4's byte count exactly: **427,872 bytes** over its ten sealed artifacts, `constraint`
  **0**, `gauss law` **0**.
* **Its one genuinely-new product — R2's custody-§8 consequence (both arms live in rows carrying no
  reopen condition, hence permanently closed by the program's own rule) — was obtained by counting
  pointers, i.e. in gate mode, not in measurement mode.**
* **Its verdict is set by lane discretion.** Control 2, measured: two defensible presupposition
  sentences, opposite verdicts.

**AS A GATE IT IS VALUABLE AND NEARLY FREE:** before testing the presupposition as a claim, **grep it
against the register.** If the register already asserts its negation, the finding is **procedural** —
*a route was taken over a contradiction already on the record* — and belongs to **W-16 gate 4**, not
to the frame measurement. That is a real and serious finding; it just is not a discovery about the
frame, and scoring it as one inflates the instrument.

### 4.5 THE PATTERN THAT EXPLAINS THE WHOLE ROUND

| test | does the instrument pre-announce a worked case? | genuinely-new content produced |
|---|---|---|
| `F2` | **YES** — W-11 leg A, `4.64e-15`, `:41-44` | **0 of 4** |
| `F4` | **YES** — R1's factorization, `:54-56` | 1 of 4, and that one in gate mode |
| `F1` | partial — positive exemplars at `:35` | 2 of 4 |
| `F3` | no | 0 of 4 |
| `F5` | **no** | **2 of 4, including the round's only certified retrodiction** |

**The two tests whose answers the instrument handed the lanes before they started produced one new
finding between them. The test the instrument said nothing about produced the round's best result.**
Deleting `F2` and demoting `F4` is therefore not only a cost cut — **it removes the two channels
through which the registrar's own answers reached the frame lanes.**

---

## 5. IS `F1` TRIVIALLY SATISFIABLE? — **YES, AS WRITTEN. PROVED WITH A NUMBER. FIXED WITH ONE CLAUSE.**

The auditor was told to attack this and was right to. `F1` says: *"Find the parameter space
containing both arms. **If one exists, the binary is a SAMPLE of a space and not a partition of
it.**"* **Two arms always lie in some space.** Measured:

* **PARITY:** `Z/2` embeds in the family `{Z/m}`. Off-axis points exist at every `m != 2`. **An
  embedding exists → by the letter, `F1` fires on evens versus odds.**
* **`diagonal / non-diagonal`:** embed in `U(3)` (real dimension 9) with coordinate `delta` = the
  off-diagonal Frobenius norm on `[0, 2.449490]`. Interior points exist, **and the decision-relevant
  functional differs there** — phase-blindness spread of `Z = <s, U s>` at fixed moduli:
  `0.000000` at `delta = 0`, `0.001798` at `delta = 0.141186`, `0.103101` at `delta = 1.000000`,
  `0.325605` at `delta = 1.414214`. **So this frame fires even under R4's stronger `(a)+(b)`
  discriminator.**

**THE MISSING CLAUSE IS EXCLUSION.** In both controls, every interior point **lies inside one of the
two arms**. The coordinate *refines an arm*; it does not exhibit a third option. Measured: integers
in neither parity arm **0 of 100000**; `U(3)` points in neither cell **0 of 20000**.

### THE AMENDED `F1` — three clauses, all of which must be met for the test to fire

1. **EXCLUSION.** Exhibit a point of the embedding space that lies in **NEITHER named arm**, and give
   the count or the construction. *(Parity: 0. `diagonal/non-diagonal`: 0. R4: the direction
   `th_11 - th_10 - th_01`, residual/response `0.808` as run, `0.846` on the auditor's recomputation.
   R2: **15 of 17** pointer-verified build targets claimed by neither arm. R1: the interior of the
   `L` axis, where a finite carrier is durable to horizon at mu-fraction `1.000`.)*
2. **RELEVANCE.** Measure the decision-relevant quantity at that point and show it differs. For a
   **route** decision the quantity is *what gets built next* and may be categorical — but the lane
   must **declare which quantity it is using**. *(This is where R2's `F1` is weak — `|D| = 5440` is a
   cardinality over a finite option set with no metric — and where it is nonetheless admissible,
   because 15 of 17 targets are genuinely unclaimed.)*
3. **PROVENANCE.** Tag each coordinate **DERIVED** (with the citation that forces it) or **COINED**
   (a lane word). **W-16 gate 1 already demands this tagging for options; extend it to coordinates.**
   *(R1's `s` = spectral incommensurability: COINED, and R1 said so. R3's stride family: COINED, and
   R3 said so. R4's spectral space: DERIVED from `REGISTER_V001.md:577`, `:1171-1172`.)*

Under clause 1 alone, both controls clear. **Under all three, `F1` retains every one of its
productive results this round and loses both false positives.**

---

## 6. WHAT DOES IT COST, AND IS THAT PROPORTIONATE?

**Measured cost of this round:** 4 frame lanes producing **32 sealed files**, **108,367 bytes of
script** and **65,898 bytes of captured output**, plus 2 audit lanes = **6 agent-lanes for 4 route
decisions**, i.e. **1.5 lanes per route decision**.

**Against the program's own baseline that is cheap.** W-08 spent **14 agents** on *one* question
(`REGISTER_V001.md:772-773`); W-11 ran **twelve** `LANE_W11_*` directories; W-10 ran **eleven**. At
1.5 lanes per route the instrument costs roughly **one tenth of a normal round**.

**After the recommended deletions it costs less than half of that** — two tests plus one grep, and
the deleted pair (`F2`, `F3`) are the two that require the most numerical work for the least return.
A frame lane becomes **one lane**.

**PROPORTIONATE — TO ROUTE DECISIONS ONLY, AND `FRAME_CHALLENGE_V001.md:14-22` ALREADY SAYS SO, FOR
THE RIGHT REASON:** *"Running it on every question would make it ceremony. Ceremony is how the
isolation ledger died."* **Endorsed without amendment.** Four route decisions in sixteen register
rows is one frame lane per four rows.

**THE REAL COST THIS ROUND WAS NOT THE LANES.** Three of four lanes spent a substantial fraction of
their output on **contamination bookkeeping** — R1's `cutoff_held` field alone runs longer than its
`F1` — because the brief ordered them to read a post-cutoff instrument that names their answers.
**Fix the brief (§7, amendment A) and the cost falls again.**

---

## 7. WHAT READS TWO WAYS — EXACTLY ONE ITEM, AND THE TIC IS DECLARED

**W-15 measured this registrar returning "reads two ways, scored neither" at 5.85x the predecessor's
rate and called it *"a template, not a judgement — it wears the form of restraint while performing a
conversion."* I am therefore spending exactly one refusal in this ruling, and I am naming what I
refused to refuse.**

**Scored rather than deferred:** the 19/20 fire rate (settled by running the control, §1) · whether
R3's catch counts (**scored: no**, §2) · whether `F1` is trivially satisfiable (**scored: yes**, §5)
· whether the instrument reports SOUND (**scored: yes, twice**, §3) · which tests to delete
(**scored: two deletions and one demotion**, §4). **Each of those could have been written as "reads
two ways." None of them is.**

### THE ONE ITEM: **CAN A FRAME LANE'S BRIEF BE KEPT CLEAN BY *THIS* REGISTRAR AT ALL?**

* **Reading 1 — it is an artifact of this round.** The leak vectors are enumerable and every one is
  in a file the registrar controls: `FRAME_CHALLENGE_V001.md:41-44` (`F2`'s number **and** its
  direction), `:54-56` (`F4`'s answer for R1), `:35-36` (`F1`'s shape for R4, `F1`'s method exemplar
  for R2), `:14-20` (the route list itself). Strip section 1's worked examples and section 0's route
  list into a registrar-only appendix and the channel closes.
* **Reading 2 — it is structural.** The registrar writes the brief; the lanes are the registrar's own
  model; and **the leak vector that fired earliest for R1 was the session auto-memory line**
  *"forced crossing does not exist; STOP-FALLS-REBUILD awaits the principal"* — present at turn zero,
  before any file was opened, reachable by **no brief edit whatsoever**.

**NOT SCORED — and the reason is not that the question is hard, it is that the deciding evidence does
not exist.** The experiment that would settle it is **already registered and has never been run**:
the **lineage-independent lane** W-03 specified, standing as a reopen condition on two rows
(`REGISTER_V001.md:87`, `:158`), unfired for anything after W-06. `FRAME_CHALLENGE_V001.md:76-78` and
W-16's *"WHAT THIS DOES NOT FIX"* both name this and neither can decide it either. **Run one frame
lane on a different model, on a route already framed by this one, and compare. Until then this stays
open, and it is the only thing in this ruling that does.**

---

## 8. THE INSTRUMENT AS IT SHOULD BE ADOPTED

> ### THE FRAME CHALLENGE — V002 (recommended)
>
> **RUNS ON ROUTE DECISIONS ONLY** (`V001` §0, unchanged and endorsed).
>
> **`F5` — THE NULL OPTION. RUN FIRST.** Exhibit a witness for each arm, or show that no witness
> exists. *"Neither" is an admissible answer requiring its own positive evidence, never a residue.*
>
> **`F1` — EMBED. RUN SECOND, AND ONLY IF `F5` LEAVES A CHOICE STANDING.** Give the space, its
> dimension, and the coordinates of each arm. **The test fires only if all three clauses are met:**
> **(1) EXCLUSION** — a point in NEITHER arm is exhibited, with a count or a construction;
> **(2) RELEVANCE** — the decision-relevant quantity is named and measured there and differs
> (categorical is allowed for routes, but must be declared);
> **(3) PROVENANCE** — every coordinate tagged DERIVED (with citation) or COINED.
> **Two mandatory reporting lines, no verdict attached:** the measure of each cell under a stated
> measure *(salvaged from `F3`)*; and, for each arm, how many of the register's standing reopen
> conditions it touches *(salvaged from R2's `mu_reopen`)*. **One mandatory clause when a map relates
> the arms** *(salvaged from `F2`)*: is it a **symmetry** or a **quotient**, and what is the fibre
> dimension?
>
> **THE PRESUPPOSITION GATE (was `F4`) — a grep at registration, not a test.** Write the sentence the
> question assumes, then grep it against `REGISTER_V001.md` **before** testing it. If the register
> already asserts its negation, register the **procedural** finding — *a route was closed over a
> contradiction on the record* — under **W-16 gate 4**.
>
> **THE BRIEF.** The frame lane's brief carries the question, the three `F1` clauses, and `F5`.
> **It carries no worked examples, no observed-failure annotations, and no route list.** Those move
> to a registrar-only appendix. **Measured cost of shipping them in V001: `F2`'s and `F4`'s entire
> output for this round.**
>
> **NOT CLAIMED EXHAUSTIVE.** Two tests is itself a coined enumeration; per W-15's discriminator it
> is a place to look for a third.

---

## 9. WHAT WOULD OVERTURN THIS RULING

**REOPENS IF:**

* a **sound ROUTE decision** — not a sound partition — is constructed and the amended `F1`+`F5` fires
  on it *(my controls bound the false-positive rate on partitions only; `CONFOUNDS.txt` C3)*; **or**
* a **lineage-independent** frame lane, run on a route this registrar already framed, returns a
  materially different enumeration *(§7; `REGISTER_V001.md:87`, `:158`)*; **or**
* `F2` or `F3`, run under a brief carrying **no worked example**, produces a genuinely-new finding on
  any route — which would show this round measured the brief and not the tests; **or**
* R1's `F5` catch is found to be pre-empted somewhere in the corpus *(checked and closed: register
  **0** occurrences of `inductive limit`, **0** of `quasi-local`; twelve sealed artifacts, **1**
  occurrence, a reference and not a test — `S2_FORMATION_CONDITION_ON_K1_AUDIT_V001.md:376-377`)* —
  reopens only if the pre-emption is found in the workflow journals or lane directories, which I did
  not scan.

**GRADE OF THIS RULING:** adversarially-checked inside one lineage, never independently corroborated
(`CUSTODY_V001.md:44-53`). **Layer fourteen. Discount it as one block with W-07 through W-16.**

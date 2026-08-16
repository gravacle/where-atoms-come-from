# W-07 — THE RECURRENCE OBSTRUCTION, ISOLATED — V001 — 2026-08-16

**The lineage-independent lane W-03 specified and nobody ran.** Corpus authored by Claude Fable 5
(trailer on all 9 commits). This lane: Claude Opus 5. Custody §4's shared-lineage caveat does not
apply between this lane and the corpus it audits. It is still ONE lane, and §7 discounts it.

Conventions, seeds, grid, dressing and isolation ledger: `LANE_W07_RECURRENCE_ISOLATION/PUBLISHED_CONVENTIONS.txt`.
All figures below regenerate from that directory; `SEALS.sha256` covers 15 files.

---

## 0. WHAT WAS READ BEFORE ANYTHING WAS COMMISSIONED

Per the process rule the program paid for. Sealed corrections read in full before any computation:
**S3 audit COR-A … COR-L** (`S3_THE_CROSSING_AUDIT_V001.md:789-800`), including **COR-F** at `:794`
and its §2.4 body at `:160-209`; **S2 audit COR-A … COR-H** (`S2_FORMATION_CONDITION_ON_K1_AUDIT_V001.md:611-648`)
with its CHOICE LEDGER and FLAG BLOCK G1-G6; **S1 in full** — the artifact no adversary ever read;
the **ERRATUM AGAINST W-02**; and W-01 → W-06 in order.

Two of those corrections are load-bearing below and neither was carried into the row it bears on:
**COR-E** (`:793`) — *"`sup|Z_k|` is a window artefact stated as an equality … the true supremum is 1.
Label it a lower bound."* **COR-F** (`:794`) — transport is not diagonal in general.

---

## 1. THE FINDING IN ONE SENTENCE

**S1's published connection has finite order 4 in `U(1)`, and every recurrence figure in this corpus
was measured there** — so the one obstruction the register calls undented was tested at the single
point of its parameter space where it cannot fail, and off that point it does not occur at all
within the ranges the corpus itself uses.

---

## 2. WHAT W-06 GOT RIGHT, REPRODUCED INDEPENDENTLY

Reproduced first, because a lane that only attacks is not a check.

**(a) The dressed observable is exactly gauge-invariant under the FULL action.** Building
`t_v = W(tree path v0→v)^{-1} s_v` from S1's own edge transports and moving **both** `a_e` and `s_v`
per S1:63 — the half S3 did not implement — gives

```
max | A(g·s, g·a) − A(s,a) |  =  3.600e-16     2000 gauge transforms x 4 observables
```
W-06 registered `4.45e-16`. Same order, independently derived. **CONFIRMED.**

**(b) The dressed algebra sees what S3's test cannot.** On S1's published connection, same branches:

```
max_v | |x_v|^2 − |y_v|^2 |   (every diagonal / S3-gauge-invariant observable) = 6.939e-18
| A_23[x] − A_23[y] |          (the dressed observable)                        = 0.384349931183
```
**S3:187-188's *"indistinguishable by EVERY gauge-invariant carrier observable"* is FALSE, and W-06's
restoration of S3-0's broad form from VACUOUS to FALSE stands.** This lane does not dent it.

**(c) W-06's recurrence figure reproduces exactly.** Same carrier, same connection, generic ready
state, dressed observable `A_23`, `k ≤ 4000`:

```
cells with separation < 1e-9 :  1000 of 4000        min over k <= 4000 : 6.729e-19
```
W-06 registered *"falls below 1e-9 on 1000 of 4000 cells"* and *"returns to 2.221e-16."* **Both
reproduced.** That agreement is what licenses treating the reconstruction as W-06's object.

---

## 3. THE FINDING — `ord(ρ) = 4`, AND IT IS THE WHOLE OF THE EFFECT

`W_F = −1` and `W_C = −i` are **S1 §6's own published instance**, reproduced exactly from S1's bytes.

```
ord(W_F) = ord(−1) = 2      ord(W_C) = ord(−i) = 4      W_F = W_C^2
< W_F , W_C >  =  < −i >  =  Z_4                     order 4
```

**S1's published connection is of finite order 4.** Its orbit in `T^2` is four points — not a dense
line, and not even the rational subtorus the W-02 erratum caught at S3/S4's headline `f=2.0, c=1.1`.
The corpus has exactly two distinguished connections and **both are arithmetically degenerate, in
two different ways, and the record names neither property.**

The dressed separation at cell `k` is `amp · |ρ^k − 1|` with `ρ` the branch ratio. `1000 = 4000/4`
is `ord(ρ)`.

### THE ISOLATION — held fixed: carrier, ready state, observable, dressing, `k`-range, code path

| connection | `arg(ρ)/2π` | `ord(ρ)` | min over `k ≤ 4000` | cells `< 1e-9` |
|---|---|---|---|---|
| **S1 PUBLISHED** `(π/3 ×3, π/2 ×3)` | `−0.250000000` | **4** | `6.729e-19` | **1000 of 4000** |
| GENERIC `√2 / √3` | `−0.317837245` | ∞ | `1.567e-04` | **0 of 4000** |
| RANDOM seed 1 | `−0.077362414` | ∞ | `1.722e-05` | **0 of 4000** |
| RANDOM seed 2 | `−0.046247977` | ∞ | `3.095e-04` | **0 of 4000** |
| RANDOM seed 3 | `+0.014316520` | ∞ | `3.057e-04` | **0 of 4000** |

**Five for five.** One variable moved. **The operative variable is `ord(ρ)` — the order of the branch
ratio in `U(1)` — finite versus infinite.** Not "the connection", which is too coarse to be a finding.

*Recorded rather than silently fixed:* leg B's first generic connection was built from `φ` and `φ²`,
and `φ² = φ+1` forces `W_F = W_C` exactly. My own confound, caught on re-read, redone at leg E with
independent irrationals. The verdict did not move; the ledger should show it anyway.

### THE SCALING — `attained` versus `approached` are different obstructions

| `K` | PUBLISHED `ord(ρ)=4`: min | exact zeros | GENERIC: min | zeros |
|---|---|---|---|---|
| `1e3` | `0` | `250` | `2.847e-03` | `0` |
| `1e4` | `0` | `2500` | `4.154e-04` | `0` |
| `1e5` | `0` | `25000` | `3.745e-05` | `0` |
| `1e6` | `0` | `250000` | `3.377e-06` | `0` |
| `1e7` | `0` | `2500000` | `3.043e-07` | `0` |

On the published connection the dressed record is annihilated **to exact zero on exactly `K/4` cells,
at every `K`, forever.** Off it the separation is **never zero at any `k`**, and the worst near-return
floor falls like `~2π/K`.

---

## 4. THE SAME CUT, ONE LEVEL DOWN — THE FOUNDING OBSTRUCTION

`FOUNDING_DESIGN_V001.md §4` promotes *"a finite discrete spectrum is recurrent"* to **the obstruction
the construction must answer at the start.** W-01 measures it (`0.0247` at `n=42`, recurring to
`0.99994`); COR-H corrects the figure; COR-E corrects it again and labels it a lower bound. On which
connection was it measured? On S1's.

`Z_k = P0·conj(W_F)^k W_C^k + PF·conj(W_F)^k + PC·W_C^k` (S2 Theorem A's three unit-modulus
coefficients). Same isolation, `k ≤ 10^6`:

| connection | `max|Z_k|`, `p` published | cells `> 1−1e-12` | `max|Z_k|`, `p` generic | cells `> 1−1e-12` |
|---|---|---|---|---|
| **S1 PUBLISHED** | `1.000000000000000` | **500000** | `1.000000000000000` | **250000** |
| GENERIC | `0.999999999998574` | **0** | `0.999999999996579` | **0** |
| S3/S4 HEADLINE `f=2.0,c=1.1` | `0.999999999997325` | **0** | `0.999999998553642` | **0** |

**`sup_k |Z_k| = 1` is ATTAINED — exactly, periodically, on a positive fraction of all cells — on
S1's published connection, and only APPROACHED everywhere else.** COR-E saw half of this and labelled
the number a lower bound; nothing in the record distinguishes attained from approached.

They are not the same obstruction. **ATTAINED** means the write is exactly undone, forever, on 25%
of cells; nothing downstream can recover it. **APPROACHED** means the obstruction is quantitative —
a Diophantine question about how fast `|ρ^k − 1|` dips as a function of `k`, against how fast the
record grows. That question has a real answer and **nobody has asked it.**

---

## 5. THE CUSTODY GAP — W-06 HAS NO ARTIFACT AND NO CODE

`w07_f_custody.sh`, run from the repo root:

```
4.45e-16        -> REGISTER_V001.md          3*sqrt(3)     -> REGISTER_V001.md
2.221e-16       -> REGISTER_V001.md          8.42e-04      -> REGISTER_V001.md
4.11e-01        -> REGISTER_V001.md          1000 of 4000  -> REGISTER_V001.md
0.5196          -> (nothing, anywhere)
```

W-01 and W-02 left sealed build+audit artifacts. W-03 and W-05 left ten lane directories of `.py`
and `.OUT.txt`. **W-06 left nothing.** Every number it produced — the dressing, the `3√3/10`
separation, the recurrence figure, the 52-subgroup sweep, the Bell `R_N` transplant
(`1.000000000000`, `4.11e-01`, `8.42e-04`), the Schmidt `[1,0]`, `dim A^G = 5`, the wedge sequence —
exists in exactly one file, the register row asserting it.

**W-06 convicted IMP-1 of entering a governing clause with *"no digest, no file:line, no named
ruling — which under custody §1 makes it flagged, not inherited."* By that test, applied to itself,
W-06's own computational content is flagged, not inherited.** The disposition of record rests on it.

Two readings, and this lane does not distinguish them: the lanes ran and the registrar failed to
seal and commit; or less than that. §2 above is evidence for the first — two of W-06's registered
figures reproduce exactly from an independent reconstruction, which is not what one expects of
numbers that were not computed.

**`3√3/10 = 0.51961524` is the exception: it does not reproduce.** Under the dressed observable that
reproduces W-06's other two figures, the separation factor on S1's published connection is
`|(−1)^k − (−i)^k| ∈ {√2, 2, √2, 0}`; a factor of `√3` requires an element of order 3, and
`⟨W_F, W_C⟩ = Z_4` contains none. Either a different observable, a different connection, or a
different normalisation produced it. Unresolvable while the code does not exist.

---

## 6. THE 52-PARTITION SWEEP IS A VACUOUS CONTROL

W-06 offers, as leg (b) of its refutation of the registrar's gauge-group hypothesis: *"sweeping all
52 partition subgroups of `U(1)^5`, exactly ONE yields S1 §4's invariant parameter count 2."*
Re-run: `Bell(5) = 52` swept, distribution `{2:1, 3:10, 4:25, 5:15, 6:1}`, the unique winner the
discrete partition. **The figure is right.**

**And it could not have come out otherwise.** A `k`-block subgroup has constants acting trivially,
so it acts through at most `k−1` parameters; `rank ≤ k−1`; `invariants = 6 − rank ≥ 7 − k`. Reaching
`2` needs `rank 4` needs `k = 5`. One line, no sweep.

**"Could not have failed" voids a CONTROL** — the program's own rule, and this is a control, not a
theorem. **Leg (b) is void.** W-06's leg (a) — that `U(1)^V` is the group of `U(1)`-valued 0-cochains
of the complex S1 §5 already deploys in `f = da`, hence intrinsic — is an argument, not a sweep, and
**it survives.** The gauge-group hypothesis stays refuted; one of its two legs does not.

---

## 7. WHAT THIS LANE DOES NOT ESTABLISH

Stated before the summary, not after it.

1. **It does not show a durable record exists.** `sup_k|Z_k| = 1` still holds generically as a
   supremum. At `k ≤ 10^6` the branch comparison returns to within `1.4e-12` of complete
   indistinguishability. For any physical reading of "record", that may erase it just as thoroughly
   as exact recurrence. **What changed is that the question is now quantitative and was previously
   treated as settled.**
2. **The reading is two-way.** Either the founding obstruction is an artefact of a degenerate
   connection and the program's starting wall is softer than recorded; **or** it survives generically
   in quantitative form and the only defect is the register's phrasing. **This lane does not
   distinguish them and does not score itself as having.**
3. **The dressed observable is a reconstruction**, inferred from two registered figures, not read
   from W-06's code, which is not on disk.
4. **This is one lane, and it is layer six.** Five consecutive layers have each been caught by the
   next and the rate has not fallen. W-06 wrote *"discount this layer too."* **Discount this one.**
   Its most likely failure mode is the corpus's recurring one — misnaming the operative variable.
   I have named it `ord(ρ)`, finite versus infinite. If that is the wrong name, everything in §3
   and §4 is a confounded headline of exactly the type now five times of record.

---

## 8. WHAT FOLLOWS, AND WHAT IS NOT MINE TO DECIDE

**The three decisions remain with the principal and nothing here takes them:** accept
STOP-FALLS-REBUILD or hold at STOP; publish the Mahler-measure note; rename the project.

What this lane does change is the **shape** of the first one. STOP-FALLS-REBUILD was issued with the
recurrence obstruction standing as *"the one thing in this corpus that no demolition and no
restoration has dented."* That sentence is measured at `ord(ρ) = 4`. The rebuild's cost line —
*"the untouched recurrence obstruction may kill it"* — is a statement about a measure-zero set of
connections.

**The cheap decisive question this opens, which is not a rebuild and does not presume one:** the
Mahler measure `λ = m(p00 + p10 x + p01 y + p11 xy)` is the entropy of an algebraic `Z^2`-action, and
the near-return depth of `|ρ^k − 1|` is governed by the same arithmetic. **Does the record's growth
rate outrun the near-return floor, or not?** Both sides are already computed objects in this corpus.
That is one calculation, on N1 — the result the record already proposes for publication — and it
decides whether the founding obstruction survives off the degenerate point.

---

## 9. CUSTODY

Built under `CUSTODY_V001.md`. Lane code, conventions, seeds and outputs sealed at
`LANE_W07_RECURRENCE_ISOLATION/SEALS.sha256`, 15 files. Every figure on this page regenerates from
that directory. **GRADE: INDEPENDENTLY-DERIVED (Opus 5 against a Fable 5 corpus), NOT
ADVERSARIALLY-CHECKED — no refuter has read it.** That is the inverse of every other grade in this
corpus, and it is the weaker half that is missing here.

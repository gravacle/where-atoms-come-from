# C-72 — THE RULING BRIEF — 2026-08-20

**For the principal, before ruling.** This document exists because writing `PROOF_V002` found a defect
in one of the program's two `PROVED` rows, the registrar applied your bar and moved the row, and no
pre-registered rule governs that move — so the ratification is yours, exactly as the `T-48` judge's
`N2` scoping was. Nothing here asks you to take a number on trust; every figure names the command.

**Current state on disk:** `C-72` is `PARTIAL`, `BLOCKED_BY T-50`. `C-71` is `PROVED`. The program
reports **one** `PROVED` row. This is reversible with one `status.py set` and a register line.

---

## 1. WHAT THE ROW CLAIMED

`C-72`'s registered `GROUNDED` cell, as it stood this morning, verbatim:

> QUANTITY: `|sum| / sum|.|` of the accumulated per-record quantity of a medium, dimensionless.
> **PREDICTION: 1 for every data pattern on an occupancy-encoded surface** (one-sign carrier
> injection); 1 under DC saturation and screening toward zero for random or DC-free data on an
> orientation-encoded surface; null (within the surface's declared unwritten tolerance) when unwritten.
> **FALSIFIER: an occupancy-encoded page whose ratio depends on the data pattern** or goes positive;
> an orientation-encoded DC-saturated track that screens, or an AC-erased track that does not.

The bar requires four things of a `PROVED` row: the quantity in units, the predicted value, what
would falsify it, and `RECORDS VERIFIED` on two structurally different real surfaces. This row named
all four. Requirement (4) is the only one `ledger/status.py` enforces; (1)–(3) are enforced by the
registrar reading the row, which is how this sat unexamined.

---

## 2. THE DEFECT

**The `1` was an algebraic identity, not a measurement.**

`ProjectModel.formation_occupancy` scores a programmed cell at `−q` and **an erased cell at exactly
zero**. A cell scored zero contributes nothing to the numerator and nothing to the denominator, so

```
|Σ| / Σ|·|  =  |−qS| / (qS)  =  1
```

for **every** 0/1 page whatsoever, for every `q`, at every `N`. There is no branch on which it could
have come out otherwise. That is `D-10`: a check that cannot fail.

**The sealed lane said so itself, and the row was written anyway.** `LANE_T34_NAND/t34_nand.py:63-68`:
the ratio is 1 *identically*, and the draws beneath it are *a consistency check of that premise, not
independent evidence*. The lane was honest; the ledger row promoted the tautology to a prediction.

**And the same model layer contradicts the premise twelve lines away.** `model/geometry.py:924-942`
returns, from one seed-11 draw sequence, both the 1000 written 0/1 pages **and** a `±5 e` over-erase
residual field for erased cells. `model/project_model.py:202-218` scores the written pages with erased
cells at exactly zero (line 212) and the residual separately (line 215). **The two models of the same
cell were never combined.**

Combine them — the model's own residual on the model's own sealed pages:

```
min  ρ = 0.9681335687351514
mean ρ = 0.9738459653763885
max  ρ = 0.9797987440417643
pages at exactly 1.0: 0        programmed fraction f ∈ [0.444, 0.552]
```

**ρ depends on the data pattern.** That is `C-72`'s **own registered falsifier**, fired by `C-72`'s own
model. The row is falsified by its own author.

### Reproduced four times, independently

| who | result |
|---|---|
| registrar (me), from `geometry.py` directly | `0.9681335687351514` / `0.9797987440417643` |
| erratum lane, own script under `ERRATUM_C72/` | match, bit for bit |
| refuter A (computation), rebuilt before reading the lane | match; also in electron units, max diff `4.44e-16` |
| refuter B (instrument), own path | match; also under flipped sign convention |

Run it yourself:

```bash
cd "/Users/bgm/MB Work/where-atoms-come-from/model" && python3 -c "import numpy as np, geometry as GE; p=GE.occupancy_patterns(); q=GE.N_E*GE.E_CHARGE; err=p['unwritten_e']*GE.E_CHARGE; rs=[abs(np.where(w==1,-q,err).sum())/np.abs(np.where(w==1,-q,err)).sum() for w in p['written']]; print(min(rs), max(rs))"
```

---

## 3. WHAT SURVIVES, AND IT IS NOT NOTHING

The substance of `C-72` — occupancy accumulates, orientation screens — **holds, and is better
understood than before.**

**The floor.** In signed electron units, a programmed cell at `−N_E` (one sign, injection-constitutive),
an erased cell at `r_i` of either sign with `|r_i| ≤ Δ`, programmed fraction `f`:

```
ρ  ≥  RHO(f)  =  ( N_E·f − Δ(1−f) ) / ( N_E·f + Δ(1−f) ),      valid for f > Δ/(N_E+Δ)
```

`RHO(1) = 1`, `RHO(0.5) = 0.904762`, `RHO(0.444) = 0.882153` (the leanest sealed page), and the bound
is **void below `f = 0.047619`**, where a nearly-erased page may screen arbitrarily.

**Confirmed and TIGHT**, not merely valid. Refuter A brute-forced the `(Σr, Σ|r|)` polytope at eleven
values of `f` and recovered the closed form to `1e-9`, with the minimiser always at every residual
pinned to `+Δ`; 20 000 adversarial edge pages produced zero violations; one cell at `Δ + 1e-9` breaks
it. The bound is therefore exactly a restatement of *"no cell exceeds the declared tolerance."*

**The discriminator, and this is the part worth having.** The floor **contains no `N`**. An orientation
surface carrying real data screens as `N^{-1/2}`. So the two encodings differ **in kind**, and their
separation **grows**:

```
N = 10³   occupancy 0.965292   orientation 0.091817    11×
N = 10⁴   occupancy 0.971562   orientation 0.030901    31×
N = 10⁵   occupancy 0.972795   orientation 0.008651   112×
```

Fitted log-log slope `+0.511`, against an all-programmed control at `0.000` and a sampling-noise
surrogate at `−0.0002`. The correction costs about `0.012` of one order out of roughly `2.7`.

**What died is a tautology. What lived is the only half that ever had a failing branch.**

---

## 4. THE PROPOSED REPAIR, AND WHY IT DID NOT SURVIVE

The erratum lane proposed a rescoped `PREDICTION` and `FALSIFIER` built on the floor: read the page
cell-by-cell on a flash tester, set `q_i = C_fg(V_t,i − V_t,neutral)`, and require `ρ ≥ RHO(f)` from
the part's own datasheet `Q_p` and `Δ`.

**Two independent refuters, with their own machinery, both returned `REFUTED`.** Thirty-two residual
defects between them. They fall into three groups, and each group alone is disqualifying.

### (a) The measurement is not baseline-robust, and the baseline does not exist

`ρ` is scale-free in `C_fg` — it cancels — but it is **not translation-free**. On one fixed histogram,
sweeping the neutral reference across the part's own read window moves `ρ` **from `0.001348` to
`1.000000`** and flips the sign of `Σq_i`.

- A neutral reference wrong by **`−4.777 e` out of 100** — 4.78% of the program window, about `0.19 V`
  on a `4 V` window — puts a **perfect page below its own floor** and fires the proposed falsifier.
- Wrong by `+5 e` the other way returns **`ρ = 1.000000` exactly** — restoring the very tautology the
  erratum exists to destroy.
- `V_t,neutral` is **on no datasheet**, and the proposed text explicitly **forbids** the one reference
  a `V_t` histogram actually supplies: the erased population's own median.

The whole discriminating range of `ρ` lives inside a `±Δ` window around the true neutral. **The
measurement requires knowing, to better than `Δ`, the very quantity `Δ` bounds.** That is circular.

### (b) The falsifier clauses fire on correct physics

- **No `f`-guard on the `N`-decay clause.** A fixed-length record written into growing pages gives
  fitted slopes `−0.572` (1000 cells), `−0.889` (100), `−0.878` (10). Every one is a healthy
  occupancy page with every programmed cell one-signed, and every one trips the clause.
- **Nominal `Q_p`.** Actual programmed charge 5% below the nominal used in the floor, with a
  physically-signed in-spec residual: `ρ = 0.900000` against a floor of `0.904762`. Below. Fires.
- **Self-contradiction.** The prediction claims scope down to `f = 0.0476`; the separation clause
  asserts a gap the floor only delivers above `f = 0.0776`. Both were stated as in scope.

### (c) The orientation half rests on an identification that is false

`ρ = M_r/M_s` was asserted, with real CoCrPt squareness put in the `0.85–0.95` class. But `ρ` is blind
to easy-axis tilt and squareness is not: with **zero** grains reversed and 30° dispersion,
`S = 0.872536` while `ρ = 0.999153`. Three-dimensional random easy axes give `S = 0.5000` with no
grain reversed at all — the proposed relation would report 25% reversed. Separately, `ε` has **no
operational definition**: a reader may set `ε := (1−S)/2` and the clause can never fire. And the
`0.85–0.95` class has **zero grounding in this repository** — `grep` returns no squareness figure
anywhere, and `CITATIONS.md` pins three anchors, none of them this.

---

## 5. WHY `C-71` SURVIVES AND `C-72` DOES NOT

This is the crux, and it is structural rather than a matter of degree.

**`C-71` claims a SIGN.** Its measurement is *the shift of a programmed page's `V_t` distribution
against the same part's own erased population* — **a within-part comparison**. It needs no absolute
reference, because both populations sit on the same instrument in the same read. A baseline offset
moves both distributions together and cancels.

**`C-72` claims a RATIO of signed quantities.** `|Σq_i|/Σ|q_i|` requires each `q_i` to be signed
against an absolute neutral, and there is no within-part comparison that supplies one.

I checked `C-71`'s own low-`f` limit rather than assuming it: over 2000 draws per `f` on 1000-cell
pages under the model's own residual, the net-positive page fraction is `0.000` at `f = 0.005, 0.01,
0.02, 0.05, 0.1, 0.5`; at `f = 0.005` the page sum is `−499 e`, an order inside the declared `5000 e`
unwritten tolerance — so such a page reads as *unwritten* and the sign claim does not apply to it.
`C-71`'s `PREDICTION` now carries that guard explicitly, and `5000 e / 100 e` is fifty cells, so the
guard binds below about `f = 0.05` on a 1000-cell page. **It is a real restriction, not decoration.**

---

## 6. THE GATE — INCLUDING WHERE I OVERSTATED IT

`model/validate_geometry.py`'s `C-71 occupancy: written ratio == 1 for EVERY pattern` compared against
a literal and could not fail. **It is removed**, and six checks stand in its place: the floor, the
discriminator at matched `N`, the widening with `N`, and three `D-15` controls in different
configurations under the identical treatment. `36/36`, chain `PASS`.

**My commit message said "six checks with real failing branches." That was too strong.** I ran the
mutation myself — occupancy's write made two-signed, i.e. occupancy behaving like orientation:

```
                 BASELINE    MUTATED
  floor            PASS        FAIL
  discriminator    PASS        FAIL
  widens           PASS        PASS   ← the structural claim's own check does not flip
  control A        PASS        FAIL
  control B        PASS        FAIL
  control C        PASS        PASS
```

**Four of six flip, not six.** And the one that does *not* flip is `widens` — the check that gates the
structural claim the whole surviving argument rests on. Refuter B measured it failing on **9 of 50
seeds** under that mutation: about 18% power. Separately, `control A`'s equality half cannot fail
under tolerance variation at all — I verified it passes at `Δ = 1, 5, 20, 50, 90, 99`, because
`RHO(1) ≡ 1` and the all-programmed page has zero erased cells for the treatment to act on.

The replacement is a genuine improvement on an identity that could never fail. It is **not** yet a
gate I would call clean, and the document's own scope cells now say so. The register needs one
correction line to match; I have not written it pending your ruling.

---

## 7. THE OPTIONS, AND WHAT EACH COSTS

| | what it does | what it costs |
|---|---|---|
| **A. Ratify `PARTIAL`** *(what is on disk)* | `C-72` waits on `T-50` for a prediction that survives refutation. Program reports **one** `PROVED` row | The headline number halves. `T-18`'s external check now has one row to offer, not two |
| **B. `PROVED`, rescoped to the sign law** | Rewrite `C-72` down to the sign statement plus the matched-`N` discriminator, keep the status | `C-72` becomes close to a restatement of `C-71` on a second encoding. The distinct content — the ratio law — is exactly what broke, so the row keeps its grade by giving up its subject |
| **C. `PROVED` unchanged, register the defect** | Status stands; the finding is recorded | The row would assert a prediction the register calls an artifact. This is the failure mode the bar exists to prevent, and I would argue against it |
| **D. Something else you see** | — | — |

**What I would rule, and the argument against it.**

**A.** Your bar requires a predicted value and a falsifier. `C-72`'s are known-defective — the
prediction is an identity, the falsifier is triggered by the row's own model — and the replacement did
not survive two refuters. A row in that state is not a row a physicist anywhere can run. `PARTIAL` is
the vocabulary's own word for it: one direction established, the other open and named in `BLOCKED_BY`.

**The honest argument against A** is that `PARTIAL` understates what is actually known. The
discriminator is computed, adversarially confirmed, structural rather than numerical, and now gated.
Someone reading `PARTIAL` will not guess that the thing which broke was the tautology and the thing
which survived was the measurement. If that reading matters more than the grade's precision, **B** is
defensible — provided the rescoped row says plainly that it no longer claims a ratio law.

**What does not change under any option:** the defect is registered, `T-50` carries the repair, the
artifact gate is gone, and `C-71` is unaffected.

---

## 8. WHERE TO READ THE PRIMARY MATERIAL

- `PROOF_V002.md` §3 — the FORMATION layer, `P-FORM-9` (the floor), `P-FORM-12` (the bar)
- `REGISTER_V001.md` — the `T-16` entry, section *THE ERRATUM*
- `ledger/status_ledger.tsv` — `C-71` and `C-72` `GROUNDED` cells as they now read
- `LANE_T34_NAND/t34_nand.py:63-68` — the lane's own statement that the ratio is 1 identically
- `model/validate_geometry.py` — the six replacement checks
- scratchpad `ERRATUM.json`, `ERRATUM_REFUTED.json` — the proposal and the two refutations in full

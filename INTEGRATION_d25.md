# INTEGRATION_d25 — T-55: the D-25 gate made unbypassable (family d25, 2026-08-21)

Deliverables of this family: `model/checks_d25.py` (the check block; it is also this family's
module — no separate layer module, per the T-55 brief) and this document.  The later T-54
core repair also removed `checks_writing.py`'s direct `RecordSurface` import/construction;
that probe now enters through `URM.surface()` before deliberately removing provenance.
No edit to project_model.py, any validate_*.py, geometry.py, grounded.py, record_model.py,
any lane, ledger, register, manifest, or sealed document is made here.  The registrar
applies the remaining edits below.

## 1. What this closes (the P-DEF-7 finding)

PROOF_V002 P-DEF-7: the D-25 guard (`URM.surface`) is real and tested, and bypassed by the
program's own validators — `model/validate_project.py` constructs `RecordSurface` directly, so
no proof number passed through the guard; and D-25 stage (2) (source-pinned external anchors
as GATED checks) was met only by a substring test (`"Weller" in ok.provenance`).

`model/checks_d25.py` supplies:

- **The construction scan** — an AST scan of the model tree that FAILS if any file except
  `project_model.py` binds, imports, or calls the name `RecordSurface` (tool-refusal form:
  the name may not appear in code at all outside the gate's home file; aliased imports are
  caught; a file that does not parse is itself an offense). Chained into the validator, any
  FUTURE file that constructs a surface around the gate fails the suite — the gate becomes
  unbypassable going forward, not merely present.
- **D-25 stage (2)'s first real gate** — the pinned measured azobenzene half-life beside the
  model's computed lifetime at the pinned dH, with stated tolerance and semantics, plus a
  power control (§4).

## 2. Current offenders (actual assembled-tree scan after the core repair)

Read-only run on 2026-08-21 (`python3 -B model/checks_d25.py`, 0.28 s) reports exactly:

```
OFFENDERS validate_project.py:4 (import); validate_project.py:20 (reference);
validate_project.py:53 (reference); validate_project.py:54 (reference);
validate_project.py:44 (reference)
```

Normalized in source order, those are:

```
validate_project.py:4  (import)     from project_model import RecordSurface, ProjectModel
validate_project.py:20 (reference)  grain = RecordSurface("CoCrPt grain", ...)
validate_project.py:44 (reference)  s = RecordSurface(nm, mech, dE, Eb, T, f0)
validate_project.py:53 (reference)  zirc = RecordSurface("Zircon U-238", ...)
validate_project.py:54 (reference)  cmb = RecordSurface("CMB photon", ...)
```

One file, five sites. No other file under `model/` binds the name.  In particular,
`checks_writing.py` is clean after its URM-routed refusal-probe repair. (The sealed lanes also
construct directly — LANE_T41_EXTERNAL among them — but they are sealed, predate the gate,
and lie outside the scanned tree by design: the scan root is `model/`.)

## 3. The exact minimal edits to `model/validate_project.py` (verbatim old → new)

The registrar patch has eight minimal semantic edits: five remove the scanned import/
constructor sites, two bind existing display rows to existing provenance-registry keys,
and one chains the ten checks.  Each `old` line occurs exactly once in the assembled tree.

**Edit 1 — line 4 (the import):**
```
old: from project_model import RecordSurface, ProjectModel
new: from project_model import ProjectModel, URM
```
(The later `from project_model import URM` at line 70 becomes redundant but is harmless;
removing it is optional.)

**Edit 2 — line 20 (the T-29 grain; provenance from the registry key "CoCrPt grain" →
Weller & Moser 1999):**
```
old: grain = RecordSurface("CoCrPt grain", "magnetic anisotropy", 3.0 * G.KB * 300, 2.0e5 * 1.26e-24, 300.0, 1e9)
new: grain = URM.surface("CoCrPt grain", "magnetic anisotropy", 3.0 * G.KB * 300, 2.0e5 * 1.26e-24, 300.0, 1e9)
```

**Edit 3 — line 36 (bind the T-33 row to its registry key; the name is used only for
construction, no check detail prints it):**
```
old: SURFACES = [("CoCrPt HDD grain", "magnetic anisotropy", 3.0 * G.KB * 300, 60.8 * G.KB * 300, 300.0, 1e9),
new: SURFACES = [("CoCrPt grain", "magnetic anisotropy", 3.0 * G.KB * 300, 60.8 * G.KB * 300, 300.0, 1e9),
```

**Edit 4 — line 40 (same, azobenzene):**
```
old:             ("Azobenzene cis/trans", "photoisomerisation", 0.60 * eV, 1.05 * eV, 300.0, 1e13),
new:             ("Azobenzene", "photoisomerisation", 0.60 * eV, 1.05 * eV, 300.0, 1e13),
```
(Alternative the registrar may prefer: keep both display names and add alias entries
"CoCrPt HDD grain" / "Azobenzene cis/trans" to `PROVENANCE` — that touches
project_model.py, which this family may not edit, so the rename is what was tested.)

**Edit 5 — line 44 (the six T-33 surfaces through the gate; all six names are now
registry keys):**
```
old:     s = RecordSurface(nm, mech, dE, Eb, T, f0)
new:     s = URM.surface(nm, mech, dE, Eb, T, f0)
```

**Edit 6 — line 53 (the zircon control; explicit provenance — it is a control surface, and
its provenance says so):**
```
old: zirc = RecordSurface("Zircon U-238", "nuclear decay", 4.27e6 * eV, 4.27e6 * eV, 300.0, 1e21)
new: zirc = URM.surface("Zircon U-238", "nuclear decay", 4.27e6 * eV, 4.27e6 * eV, 300.0, 1e21, provenance="control surface, census GR1 entry 3 (zircon U-Pb): decay is temperature-independent; the model must DECLINE")
```

**Edit 7 — line 54 (the CMB control, same pattern; census GR1 entry 4):**
```
old: cmb = RecordSurface("CMB photon", "free flight", 0.0, 0.0, 2.7, 0.0, thermal=False)
new: cmb = URM.surface("CMB photon", "free flight", 0.0, 0.0, 2.7, 0.0, thermal=False, provenance="control surface, census GR1 entry 4 (CMB photon polarisation): no bath, free flight; the model must DECLINE")
```

**Edit 8 — chain-in, immediately after line 82 (the registry check), before
`print("=" * 78)`:**
```
old: check("D-25 registry supplies pinned provenance", "Weller" in ok.provenance, ok.provenance[:60])
new: check("D-25 registry supplies pinned provenance", "Weller" in ok.provenance, ok.provenance[:60])
# D-25 (T-55): construction scan + stage-(2) external-anchor gates — see model/checks_d25.py
from checks_d25 import run_d25_checks
run_d25_checks(check)
```
(The three new lines are inserted at column 0, no indentation.)

## 4. Where the checks chain in, and what each gates

`run_d25_checks(check, root=None)` takes the host validator's `check(name, cond, detail)`
(the validate_geometry.py idiom) and, optionally, an override scan root (default: the
directory of checks_d25.py — the model tree; the override exists so a candidate tree can be
audited before integration). Chained per Edit 8 it runs inside `validate_project.py`, which
`validate_geometry.py` already chains, so one `python3 model/validate_geometry.py` exercises
everything. The ten checks:

1. `D-25 scan: no RecordSurface binding outside project_model.py in the model tree` — the
   unbypassability gate. FAILS pre-integration (lists §2's offenders), PASSES after Edits 1–7.
2. `D-25 scan control: planted direct construction is detected` — D-15 positive control
   beside the scan's zero (synthetic tree, planted offender flagged, clean file unflagged).
3. `D-25 scan probe: aliased import (as RS) is detected` — definition-not-shortcut probe:
   the scan keys on the imported name, not call text.
4. `D-25 refusal probe: whitespace-only provenance is refused` — API fidelity beyond the
   validator's gated blank case.
5. `D-25 refusal probe: corner tier refuses lowercase 'def-a'` — the DEF-A self-declaration
   is exact.
6. `D-25 positive control beside the refusals: declared surfaces construct through the gate`
   — the gate returns on declared surfaces (world via registry, corner via DEF-A).
7. `D-25 anchor (envelope)` — measured mean lifetime t1/2/ln2 = 1.731e5 s inside the model
   envelope [9.221e1, 9.221e5] s at pinned dH = 0.915 eV over f0 ∈ [1e9, 1e13] /s.
   Stated semantics: four-decade-envelope containment — weak contact, no prediction.
8. `D-25 anchor (Eyring point)` — model t1/2 = 4.170e4 s at the pinned (dH, dS) vs measured
   1.2e5 s: −0.46 decades, inside the stated ONE-DECADE tolerance. Stated semantics:
   order-of-magnitude consistency with both parameters supplied by the pinned source; the
   residual factor 2.9 is that source's own internal spread; per the C-69 demotion note the
   closed-form agreement is an identity of the construction — the model's own content is the
   Liouvillian-mode readback.
9. `D-25 anchor power control` — the pre-correction midpoint convention fails the same gate
   at −5.37 decades (the 5–6 order contradiction the C-69 demotion note records): the gate
   can reject a wrong convention.
10. `D-25 anchor cross-check` — the gated lifetime equals the independent closed form to
    2.4e-16 relative (D-8: the decision-path number is computed through the mode extraction).

This replaces the substring test as D-25 stage (2)'s first real gate. The substring check
itself (line 82) is KEPT — it gates the registry-lookup mechanics, P-DEF-7's narration names
it, and deleting a named gate would stale the proof; its role is demoted from stage-(2)
stand-in to registry-mechanics check.

Pinned anchors carried as literals (D-8: they ARE the anchors, stated as such in the file
header): t1/2 = 1.2e5 s (1.4 d, benzene, 35 C), dH = 0.915 eV, dS = −50.2 J/mol/K — ACS
Cent. Sci. PMC9951306 via LANE_T41_EXTERNAL/CITATIONS.md item 1; dE = 0.60 eV is the corpus's
azobenzene value (enters only the negligible reverse rate). Every model-side number on a
decision path is computed through `URM.surface` + `ProjectModel.lifetime`.

## 5. ProjectModel method signatures to add: NONE, and the observation-entry story

This family adds no layer method: the D-25 gate already lives in `URM.surface`, and its laws
are ProjectModel's existing ones. The family's additions are gates.

How a NEW observation of this family's kind enters the URM:

- **A new record surface** enters ONLY through `URM.surface(name, mechanism, dE, E_b, T, f0,
  thermal=, provenance=, tier=)` — either register its pinned source in
  `project_model.PROVENANCE` under its name, or pass `provenance=` explicitly (controls state
  they are controls); corner carriers self-declare `provenance="DEF-A"`. Constructing
  `RecordSurface` directly anywhere else under `model/` now fails the validator suite: the
  scan (check 1) is the enforcement.
- **A new external measured number** enters as a new gated anchor block in
  `model/checks_d25.py` following the azobenzene pattern: (i) pin value, conditions, and
  source in LANE_T41_EXTERNAL/CITATIONS.md (or a successor citations lane); (ii) carry the
  pinned values as declared anchor literals in the file header; (iii) compute the model-side
  number through `URM.surface` + the layer method — never write it down; (iv) state the
  tolerance and the semantics (what the comparison does and does not establish) in the
  check's detail string; (v) put a power control beside it — a wrong convention or mutated
  input that the same gate measurably rejects (the azobenzene midpoint history is the
  template); (vi) gate any zero with a positive control (D-15).

## 6. Verification record

Pre-integration audit that exposed the bypass (preserved as history):

- `python3 -B model/checks_d25.py`: **9 PASS, 1 FAIL, exit 1, 0.28 s**.  The sole
  failure is the intended construction scan and lists exactly the five sites in §2;
  all nine refusal/anchor/control checks pass.
- `python3 -B model/checks_writing.py`: **57 PASS, 0 FAIL, exit 0, 5.37 s**; its
  D-25 refusal probe uses `URM.surface()` and contributes no scan offender.

Final integrated record after the eight registrar edits:

- `python3 checks_d25.py` standalone: **10 PASS, 0 FAIL, exit 0** (< 5 s).
- `python3 validate_project.py` (patched): **24 PASS, 0 FAIL, exit 0** — the live file's
  14 gates (the 11 T-46-era gates plus the three D-25 guard gates that landed with T-46)
  plus the 10 d25 gates (the run is dominated by the pre-existing T-28 [[6,4,2]]
  eigendecomposition, as before; the d25 block adds < 5 s).

Edits 1–7 are now landed and edit 8 chains the clean scan plus all nine
positive/refusal/anchor controls into `validate_project.py`. The live standalone D-25
block is 10/10; the earlier 9/1 result above is the dated pre-integration state that made
the five-site bypass visible.

Post-integration count: `validate_project.py` prints **24 PASS** (14 + 10).
`validate_geometry.py` now narrates that count correctly, and `validate_urm.py` chains the
project/D-25 validator after its 176 family gates and the 33 geometry gates.

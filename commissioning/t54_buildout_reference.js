export const meta = {
  name: 'T54-T55-urm-buildout',
  description: 'T-54 + T-55: home the unhomed machinery into the URM (T-46 pattern) and make D-25 unbypassable — the URM as the framework new observations enter',
  phases: [
    { title: 'Build', detail: 'five builders in parallel, standalone modules + integration specs, no shared-file edits' },
    { title: 'Verify', detail: 'two adversarial verifiers, T-46 standard: re-derive from sealed lanes, API-fidelity probes' },
  ],
}

const REPO = '/Users/bgm/MB Work/where-atoms-come-from'

const PREAMBLE = `
REPO: "${REPO}"   — THE PATH CONTAINS A SPACE. Always quote it in every shell command.

You are building out the UNIVERSAL RECORD MODEL (task T-54/T-55). THE PRINCIPAL'S DIRECTIVE,
2026-08-21, governing this work: **the URM is the world model — the framework new observations are
added INTO.** Not a museum of past results: an extensible object where a new record surface enters
through the provenance gate, a new law enters as a layer method with a validator gate, and a new
external number enters as a gated comparison. Build accordingly.

READ FIRST: model/project_model.py (the whole file — the layer architecture and D-25 gate),
model/geometry.py's header (the T-46 fold-in pattern you are replicating: ported machinery, each
function's docstring naming its claim row, sealed source, and owners), model/VERIFY_T46 (what the
verifier will hold you to), and your family's sealed lanes named below. Then your ledger rows.

FILE DISCIPLINE — COLLISION SAFETY IS ABSOLUTE: you write ONLY (a) one NEW module file under model/
named for your family, (b) one NEW check-block file model/checks_<family>.py containing check()
calls in the validate_geometry.py idiom (def run_<family>_checks(check): ...), and (c) one
INTEGRATION_<family>.md at the repo root of your module documenting: the ProjectModel method
signatures to add (signature + docstring text, claim row, sealed source, owners), where your checks
chain in, and the observation-entry story (how a NEW observation of this family's kind would enter
through your functions). NEVER edit project_model.py, validate_*.py, geometry.py, grounded.py,
record_model.py, any lane, the ledger, the register, or any document. NEVER run git or reproduce.sh.
The registrar integrates.

STANDARDS (the T-46 bar, and the verifier re-derives everything):
 - Ported machinery reproduces SEALED numbers exactly — your check block must gate them, plus at
   least one API-fidelity probe BEYOND the gated range (definition-not-shortcut).
 - D-8: no literal expected value on any decision path where a computed comparison is possible;
   where a sealed number IS the anchor, gate it as the sealed anchor it is, stated as such.
 - D-15: every zero gated with a positive control beside it. D-25: world-tier inputs through
   provenance. Owners named. Exact arithmetic where the lanes used it (F_2, Fractions).
 - Your module must be importable standalone AND cheap: total check runtime under ~120s.
`

phase('Build')

const FAMILIES = [
  { key: 'arrow', brief: `THE ARROW LAYER (the family §3 of the proof says is owed a LAYER).
Rows: F-17 (chi threshold at weight d: 0.00000000 at weight 1 for all 24 observables, 0.11448276 at
weight 2), F-18 (weight-1 entangles I(S:B)=0.04549256 yet transfers ZERO record bits), F-19
(system-only operations cannot reduce I(S:B), invariance 3.686e-14), F-21 (fragment redundancy
0.789/0.048/0.679 under weight-d, exactly zero under weight-1). Sealed lanes: LANE_F1_ARROW,
LANE_PF2_DYNAMICAL, and the T-9 battery entries for F-20/F-21. Existing machinery:
record_model.RecordModel.formation (with fragment=), and RecordModel.redundancy — which no validator
calls today; wire it. Module: model/arrow.py. The observation-entry story: how a new bath/fragment
observation would be scored through these functions.` },
  { key: 'countlaw', brief: `THE COUNT LAW HOMED (C-86; C-76's kill and what survived it).
The census law k(t_m): both widths from clause (ii') applied to the record's own mode with zero
adjustable content — delta_pop(t_m) = kT*ln(expm1(B/kT - ln f0 t_m)), delta_coh = hbar/t_m — and the
dated staircase with drop times t*_i = f0^-1 exp((B_i-dE_i)/kT)/(1+e^{-dE_i/kT}). Sealed lanes:
LANE_T47_A_WIDTH..T47_D_REGISTER; model/count_law.py (the corner k = min_E v2(m_E), already model-side
— reference it, do not duplicate). Module: model/countlaw.py with a ProjectModel-ready census method
taking a LIST of RecordSurface objects and t_m, returning k and the drop schedule. Gate the sealed
staircase numbers + the T-31 asymmetric-well behavior + the C-76 kill as a control (the chosen-width
form must NOT be reachable through your API). Observation-entry: how a real device census would run
through it (this is the machinery the Saira/Woodside grounding lanes will call).` },
  { key: 'classes', brief: `THE REACHABLE CLASSES HOMED (C-87 — appears nowhere in model/ today).
A Gamma-priced coupling reaches exactly three classes: exponential below mu_c, the earned dimension's
own critical Green's class at mu_c (the critical identity IS the venue's discrete Poisson equation),
divergent above. Sealed lanes: LANE_T44_A.. (the kernel/corner machinery), LANE_T44_B_WORLD (the D=3
member: mu_c = 1/6 exact three ways, exponent bracket containing 1, G in [1.5039,1.5544], d*G(d) onto
[0.476369,0.487321]). Module: model/classes.py: the signed/priced kernel on a declared venue graph,
mu_c located by the resolvent route (exact rationals), the class verdict per mu as a computed
boolean triple. Gate the sealed D=3 numbers + the cross-dimension discriminator (D=1 linear / D=2
log / D=3 power on the same instrument) + one API-fidelity probe off the sealed grid.
Observation-entry: a new venue graph enters as an adjacency structure with declared provenance.` },
  { key: 'writing', brief: `THE WRITING TIER HOMED (C-91 — the E1/E2/E3 ensembles; the family C-93's
responsive-venue computation needs). Sealed lanes: LANE_T48_A_DERIVATION (unitarity forces
conservation for every weighting; conserving <=> critical, det(I-K)=0 with exact nonzero controls),
LANE_T48_B_CORNER (link-uniformity earned from the writer algebra, 1/deg = mu_c in-lane),
LANE_T48_C_WORLD (E1 transport critical at every dE; E2 conserving, uniform at dE=0; E3 the model's
own erase channel, never critical, gap ln(mu_c/mu) = ln(1 + e^{dE/kT}/l) with f0, E_b dropping out).
Module: model/writing.py: the three ensemble constructors on a declared venue + surface constants,
the conservation/criticality verdicts as computed booleans, the closed-form gap CHECKED against the
computed one (never sourced from it). Gate the sealed identities + controls (CTRL-BIAS-LINK,
CTRL-LEAK). Observation-entry: this is where 'does one written record shift dE for an adjacent
write?' (C-93's named next step) will run — state the entry point for that computation explicitly
in the integration doc, but do NOT run it (it is T-49-adjacent physics, commissioned separately).` },
  { key: 'd25', brief: `T-55 — THE D-25 GATE MADE UNBYPASSABLE. The proof's P-DEF-7 finding: the
guard is real, tested, and bypassed by the program's own validators (validate_project.py line ~20
constructs RecordSurface directly), so no proof number passed through it. Your deliverables differ
from the other families': (a) model/checks_d25.py — a check block that greps the model tree and
FAILS if any file outside project_model.py constructs RecordSurface directly (the tool-refusal form
of the rule; list current offenders in your integration doc); (b) INTEGRATION_d25.md — the exact
minimal edits to each offending validator to route through URM.surface() (verbatim old/new lines —
the registrar applies them); (c) in the same check block, ONE gated external-anchor comparison: a
pinned measured value (from LANE_T41_EXTERNAL/CITATIONS.md — e.g. the azobenzene t1/2 = 1.2e5 s at
35 C in benzene) placed beside the model's computed lifetime envelope for the pinned dH = 0.915 eV,
with a STATED tolerance and semantics, replacing the substring test as D-25 stage (2)'s first real
gate. State honestly in the check's detail string what the comparison does and does not establish
(the azobenzene convention history is the cautionary tale — read the C-69 demotion note first).` },
]

const built = await parallel(FAMILIES.map((f) => () => agent(
  `${PREAMBLE}\n\nYOUR FAMILY: ${f.key}\n\n${f.brief}\n\nBuild, test standalone (run your module and your check block; every check PASS), and return: the file paths, the check names with their measured values, the integration spec summary, and anything that surprised you.`,
  { label: `build:${f.key}`, phase: 'Build' }
)))

phase('Verify')

const verdicts = await parallel([
  () => agent(`${PREAMBLE}
YOU ARE VERIFIER ONE (families: arrow, countlaw, d25), default REFUTED, the VERIFY_T46 standard:
(1) run every check block standalone and confirm PASS counts; (2) re-derive the gated numbers FROM
THE SEALED LANE OUTPUTS DIRECTLY (grep the sealed .txt, then call the new API fresh — never read the
number from the check source); (3) API-fidelity probes on untested cases — definition, not shortcut;
(4) D-8/D-15 scan of every check; (5) for d25: confirm the tree-grep check actually fails on the
current offenders and that the integration doc's edits are minimal and correct.
BUILDER RETURNS: ${JSON.stringify(built.filter(Boolean).slice(0, 3)).slice(0, 9000)}
Write your working to "${REPO}/LANE_T54_VERIFY/" part A. Return verdict per family with your numbers.`, { label: 'verify:A', phase: 'Verify' }),
  () => agent(`${PREAMBLE}
YOU ARE VERIFIER TWO (families: classes, writing), default REFUTED, the VERIFY_T46 standard — same
five duties as your sibling, on the two heaviest families: re-derive the D=3 member numbers and the
E1/E2/E3 identities from LANE_T44_B_WORLD and LANE_T48_* sealed outputs directly, then through the
new API fresh; attack the exact-rational claims (are they actually exact on the measurement path?);
probe one venue OFF the sealed set per family.
BUILDER RETURNS: ${JSON.stringify(built.filter(Boolean).slice(3, 5)).slice(0, 9000)}
Write your working to "${REPO}/LANE_T54_VERIFY/" part B. Return verdict per family with your numbers.`, { label: 'verify:B', phase: 'Verify' }),
])

return { built: built.filter(Boolean).map((b) => (typeof b === 'string' ? b.slice(0, 1200) : b)), verdicts: verdicts.filter(Boolean).map((v) => (typeof v === 'string' ? v.slice(0, 2500) : v)) }

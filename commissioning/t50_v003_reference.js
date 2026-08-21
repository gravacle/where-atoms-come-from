export const meta = {
  name: 'T50-v003-terminal',
  description: 'T-50 DESIGN ONE V003: the identity-level repairs, under the pre-registered terminal boundary — this round ends the question either way',
  phases: [
    { title: 'Build', detail: 'V003 with the three judged repairs, each with its measured failing branch' },
    { title: 'Refute', detail: 'both refuter roles, strongest prior kills first' },
    { title: 'Judge', detail: 'survival, or the pre-registered fallback rescope — never no-runnable-form' },
  ],
}

const REPO = '/Users/bgm/MB Work/where-atoms-come-from'

const PREAMBLE = `
REPO: "${REPO}"   — THE PATH CONTAINS A SPACE. Always quote it in every shell command.

T-50 DESIGN ONE V003 — THE TERMINAL ROUND, under a boundary REGISTERED IN THE PROGRAM'S RECORD
before this workflow launched (REGISTER_V001.md, the V002 entry): if V003 falls to another
D1-density-shape-family kill, the outcome is the FALLBACK RESCOPE — the falsifier-grade 2-SE claim
(never fired falsely in any round) becomes the registered prediction, point band scoped or
withdrawn — and NEVER "no runnable form." This round ends the T-50 observable question either way.

READ FIRST, in full: "${REPO}/LANE_T50/V002/JUDGMENT_V002.txt" (especially §4, the named repairs,
each carrying its D-8 obligation), both V002 refuter verdicts (LANE_T50/V002/VERIFY_A/VERDICT_A2.txt
and VERIFY_B/), v2_design.md and the sealed v2_pipeline.py, and the register's V002 entry.

WRITE ONLY in "${REPO}/LANE_T50/V003/" (create it) and scratchpad. Never edit anything else; never
run git / reproduce.sh / r3.sh. Binding: the principal's ruling and three constraints (the T-50 plan
row), the accumulation directive, D-8, D-15, D-1.
`

phase('Build')

const design = await agent(`${PREAMBLE}
YOUR TASK — BUILD V003 from the sealed V002 pipeline with EXACTLY the three judged repairs, each
with its measured failing branch or the condition is not registered:
 (1) C1 closed at the IDENTITY level: bound the RELATIVE trend |d log10 medf / d log10 N| <= tol_rel
     over THE FIT'S OWN admissible window (guard point onward) — subsuming the scale gap and the
     window gap. tol_rel co-tuned over the FULL admissible density range INCLUDING the
     low-density-skew diagonal and the island literals; the sparse-pattern discrete-median bias at
     N=16-32 re-measured; the band arithmetic stated correctly (tol_rel bounds the exponent shift
     identically and must sit inside the 0.1 half-band with named margin). Failing branch: every
     cascade and island mask of BOTH V002 kills, at measured rates.
 (2) The point certificate counts surrogates outside BOTH band edges for the point sentence (lower
     edge alone still governs SEAM/falsifier arming). Failing branch: the K1/K2 kill reads must fail
     certification or void, at measured rate.
 (3) The orientation data-balance condition registered on the reader's own write record (threshold
     on the declared pattern's one-way excess, co-tuned), wired into the fire_c path AND the
     guard-scope table — clause (c) must never again fire on data the prediction makes no claim
     about. Failing branch measured.
Change NOTHING else — V003 is V002 plus these three; every other check, control and certificate
carries over byte-comparable. Deliver v3_design.md (the deltas stated as deltas), v3_pipeline.py ->
v3_run.txt, v3_mutations.py -> v3_mutations.txt (the V002 suite PLUS the K1/K2/A-K1 kill masks as
permanent members, every designated cell at its rate). Return the summary, the three failing-branch
rates, and honest risks.`, { label: 'build:V003', phase: 'Build' })

phase('Refute')

const refs = await parallel([
  () => agent(`${PREAMBLE}
REFUTER A — COMPUTATION, default REFUTED. THE V003: ${JSON.stringify(design).slice(0, 9000)}
Your prior kills first: the D1 density-shape family in EVERY moment and window — new mask shapes of
your own construction beyond the suite's; then gain/offset/quantifier attacks. THE BOUNDARY MATTERS
FOR YOUR REPORT: distinguish explicitly whether any kill you land is D1-FAMILY (density-shape
through the C1 condition) or a NEW family — the judge's branch turns on exactly that. Work in
LANE_T50/V003/VERIFY_A/. Return verdict, defects (family-tagged), your numbers.`, { label: 'refuteA:V003', phase: 'Refute' }),
  () => agent(`${PREAMBLE}
REFUTER B — INSTRUMENT, default REFUTED. THE V003: ${JSON.stringify(design).slice(0, 9000)}
The bar unchanged: a stranger with a flash tester and a polar Kerr microscope runs this from the
protocol alone. Your prior kills first (the certificate's edges, drift, fixed pattern, gain, the
blank/railed/full pages, every guard's honest-choice sweep, every falsifier hunted for an
ordinary-physics trigger) — and tag every kill D1-FAMILY or NEW. Work in LANE_T50/V003/VERIFY_B/.
Return verdict, defects (family-tagged), your numbers.`, { label: 'refuteB:V003', phase: 'Refute' }),
])

phase('Judge')

const judgment = await agent(`${PREAMBLE}
THE JUDGE, adding no measurement, scoring under the REGISTERED boundary. INPUTS:
DESIGN: ${JSON.stringify(design).slice(0, 6000)}
REFUTATIONS: ${JSON.stringify(refs.filter(Boolean)).slice(0, 12000)}
Rule: (i) V003 SURVIVES both refuters -> draft the register entry, the C-72 GROUNDED replacement
(V003's prediction + falsifier), the mutation-suite permanence note, and the flagged-for-principal
question of C-72's return to PROVED (both DONE_WHEN halves would now be met — drafted, never ruled).
(ii) V003 falls to a D1-FAMILY kill -> THE BOUNDARY FIRES: draft the fallback rescope — the
falsifier-grade 2-SE claim as the registered prediction, point band scoped or withdrawn, with the
exact GROUNDED text — and the register entry saying the boundary fired as registered.
(iii) V003 falls ONLY to a genuinely NEW family -> the boundary does not cover it; rule honestly
between one further named repair and escalation to the principal, with the new family's repair
prospects stated. Write LANE_T50/V003/JUDGMENT_V003.txt. Return verdict, branch, drafts.`, { label: 'judge:V003', phase: 'Judge' })

return { design: typeof design === 'string' ? design.slice(0, 1500) : design, verdicts: refs.filter(Boolean).map((r) => (typeof r === 'string' ? r.slice(0, 400) : r)), judgment: typeof judgment === 'string' ? judgment.slice(0, 3000) : judgment }

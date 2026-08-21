export const meta = {
  name: 'LANE-T51A-second-lump',
  description: 'LANE_T51_A: the Second Lump first computation — connectivity gate, F(D) on toric 3x3, adversarial verification, judged against the pre-registered V1-V5',
  phases: [
    { title: 'Gate', detail: 'redefine w_enc over connected enclosing strings; resolve the 3x2 counterexample' },
    { title: 'Measure', detail: 'F(D) sector-exact on toric 3x3 with all five controls' },
    { title: 'Verify', detail: 'adversarial rebuild with independent machinery, fresh points off the grid' },
    { title: 'Judge', detail: 'score V1-V5 verbatim; name the next step per branch' },
  ],
}

const REPO = '/Users/bgm/MB Work/where-atoms-come-from'

const PREAMBLE = `
REPO: "${REPO}"   — THE PATH CONTAINS A SPACE. Always quote it in every shell command.

You are an agent of the "Where Atoms Come From" program building LANE_T51_A — the first computation
of the field-side instrument THE SECOND LUMP. THE COMMISSION IS FIXED AND PRE-REGISTERED in
"${REPO}/FIELD_INSTRUMENT_V001.md" §1 — read it in full FIRST, then the full design at
"${REPO}/LANE_T51_IDENT/T51_DESIGN_the-second-lump.json", then BOTH Second Lump critiques in
"${REPO}/LANE_T51_IDENT/CRITIQUES.jsonl" (their repairs are BINDING — chains fail by being
under-read), then CORE_FRAMEWORK_V001.md (the anchor), then the C-92 row:
    awk -F'\\t' '$1=="C-92"' "${REPO}/ledger/status_ledger.tsv"

WRITE ONLY inside "${REPO}/LANE_T51_A/" (create it) and your scratchpad. NEVER edit the ledger, the
register, any document, model/, replicate/, or any other lane. NEVER run git. NEVER run
replicate/reproduce.sh or r3.sh.

BINDING DISCIPLINE:
 - D-1 ABSOLUTE: no classical gravitational form is required of, or tested against, anything here.
   The shape of F vs separation is an OUTPUT. The word Newton must not appear in this lane.
 - THE PRINCIPAL'S DIRECTIVE (2026-08-20, binding, quoted in C-92): "the mechanism for accumulation
   is whatever it proves to be." No outcome of any computation in this lane is a failure against an
   imported standard; every outcome registers as the surface's own law. No "kill" framing anywhere.
 - D-8: verdicts are COMPUTED BOOLEANS with both branches reachable; no literal expected values on
   any decision path; every fit carries a noise floor beside it.
 - D-15: every reported zero has a positive control beside it in a DIFFERENT configuration.
   Construction certificates (algebraic identities) are labeled certificates, never controls.
 - D-24: separations are stated ONLY in earned quantities (d_gen, connected w_enc); no lattice
   coordinate enters any claim.
 - Exact where feasible: F_2 arithmetic for all coset/connectivity computations; sector-exact
   diagonalization (declare float tolerances for eigenvalues); state every tolerance before use.
 - Owners named for every borrowed idea. Errors go in the lane's D24_AUDIT.txt, never narrated.
`

phase('Gate')

const gate = await agent(
  `${PREAMBLE}

YOUR TASK — THE CONNECTIVITY GATE (must land before any V3 scoring; the audit's computed
counterexample is why): the sign-constrained coset minimum w_enc admits DISCONNECTED representatives
that the dynamics cancels — at 3x2 the coset minimum is 3 while the measured onset converged to 4.0.

PRIOR ART TO CHECK AGAINST, NOT TO COPY: the fourth-angle rigidity lane (C-93) independently
computed connected w_enc at the 3x2 calibration AND four-plus 3x3 placements, with onset orders
(LANE_T51_IDENT/FOURTH_rigidity.py and FOURTH_rigidity_OUT.txt, plus the critic's rebuild). It found
a SECOND counterexample at the commissioned 3x3 venue itself: old-spec 4, connected 5, measured
onset 4.99. Your independent computation must be built first and compared after — agreement is a
cross-check, disagreement is a finding to report, never to reconcile silently.

Build, in "${REPO}/LANE_T51_A/":
 1. g1_connected_wenc.py — the redefined quantity: w_enc_conn = the minimal weight of a CONNECTED
    admissible enclosing string (connected in the venue's edge-adjacency graph restricted to the
    string's support), computed by exact F_2 machinery. Reuse geometry.py's coset method as the
    starting point (import it; do not copy-paste silently — name what you reuse). Enumerate or
    branch-and-bound over the coset with a connectivity filter; PROVE minimality at these sizes
    (exhaustive within a declared weight bound, bound stated and justified).
 2. Verify the counterexample resolves: at the 3x2 calibration placement, print the disconnected
    coset minimum (expected from the audit: 3), the connected minimum, and the audit's measured
    onset 4.0 beside them — the connected minimum must be checked against the onset by a COMPUTED
    comparison, not asserted. If the connected minimum does NOT match the measured onset, that is a
    FINDING — report it plainly; do not tune the definition until it matches (D-8 in definition
    space). Also compute both quantities for every placement the 3x3 run will use.
 3. Sealed output g1_connected_wenc.txt with every number, the enumeration bound, and the
    connectivity definition stated precisely enough that a stranger reimplements it.

Return: the connected w_enc values for all placements (3x2 calibration + both 3x3 separations), the
counterexample resolution verdict, and anything that surprised you.`,
  { label: 'gate:connected-wenc', phase: 'Gate' }
)

phase('Measure')

const measure = await agent(
  `${PREAMBLE}

THE GATE RESULT (input, not to be re-run — its lane files are on disk in LANE_T51_A/):
${JSON.stringify(gate).slice(0, 6000)}

YOUR TASK — THE MEASUREMENT, exactly as commissioned in FIELD_INSTRUMENT_V001.md §1:

Venue: toric 3x3 torus (18 edges), sector-exact blocks via the diagonal conserved algebra
(plaquettes, source values, Z-winding pair); the 3x2 run as calibration (reproduce the design
exploration's numbers first — F(contact) = -1.99e-3, F(far) = +1.94e-4 at lam = 0.05 — as a
calibration CHECK, stated as reproduction, before any new number).

Probe: star-hole pair at adjacent vertices, connector weight 1. Source: plaquette-hole pair at the
two distinct earned separations 3x3 affords, each stated as (d_gen, w_enc_conn) from the gate.
Mediator: V = lam * sum_e Z_e, lam in {0.02, 0.05, 0.10}.

Compute per (placement, lam): Delta(b=+1), Delta(b=-1), F = the difference; onset order of F in lam
(log-log fit WITH noise floor printed beside it).

CONTROLS, all five, two-way, in the same tables:
 C1 F written vs unwritten (source-absent).
 C2 Gamma-equivalent placement swap at equal earned separation.
 C3 WINDING-SECTOR SWEEP — compute F in every winding sector; MANDATORY before any sign of F is
    attributed to content. Report the decomposition plainly.
 C4 onset-order bracket against connected w_enc at both placements.
 C5 back-action: source-sector energies with probe present/absent, against a DECLARED tolerance.
CERTIFICATES (labeled as such, never counted as controls): unwritten-equals-source-absent
(conserving-quadrature identity); X-quadrature probe quiescence (commutation identity).

Scripts a1_*.py .. a5_*.py (or fewer, your structure), each with a sealed .txt. Hilbert-space sizes
and block dimensions PRINTED and checked against the declared ~2^11 bound before any eigh runs. If a
computation exceeds memory/time at 3x3, SAY SO and deliver the largest exact venue that fits — never
silently sample.

Return: the full F table, all control outcomes, the calibration verdict, and every tolerance used.`,
  { label: 'measure:F-of-D', phase: 'Measure' }
)

phase('Verify')

const verify = await agent(
  `${PREAMBLE}

YOUR TASK — ADVERSARIAL VERIFICATION, default REFUTED. A builder produced LANE_T51_A's gate and
measurement (files on disk in "${REPO}/LANE_T51_A/"; its returns below). Rebuild the load-bearing
computations with INDEPENDENT MACHINERY — your own state construction, your own sector reduction
(or full-venue where it fits), fresh lambda values OFF the declared grid, at least one placement
re-derived from scratch — and attack:
 A. the sector exactness claim (is the source bit truly an exact quantum number under this mediator?
    prove or refute by commutation computation, not citation);
 B. the connectivity gate (re-derive w_enc_conn independently for one placement; try to construct a
    lighter connected enclosing string than the builder's minimum);
 C. the winding-sector attribution (C3): can you make the sign of F flip by moving winding sector
    alone? If yes, is any sign claim in the builder's tables thereby unlicensed?
 D. the onset fits (refit with your own estimator and noise floor; check the bracket verdicts);
 E. the D-1/directive scan: any sentence in the lane's outputs that requires a shape, imports a
    standard, or frames an outcome as failure.
GATE + MEASUREMENT RETURNS:
${JSON.stringify({ gate, measure }).slice(0, 9000)}

Write your rebuild into "${REPO}/LANE_T51_A/VERIFY/" with sealed outputs. Return: verdict
NOT_REFUTED or REFUTED per item A-E with your own numbers beside the builder's, and any finding.`,
  { label: 'verify:independent', phase: 'Verify' }
)

phase('Judge')

const judgment = await agent(
  `${PREAMBLE}

YOUR TASK — THE JUDGE. You add NO measurement. Score the PRE-REGISTERED rule from
FIELD_INSTRUMENT_V001.md §1 VERBATIM — read it from the document, quote each of V1-V5, and score
each as a computed boolean from the lane's sealed numbers (builder + verifier). The rule was
registered before the numbers existed; you apply it on its own stated conditions, exactly as the
T-44/T-48 judges did.

INPUTS: everything in "${REPO}/LANE_T51_A/" including VERIFY/; the returns:
${JSON.stringify({ gate, measure, verify }).slice(0, 9000)}

 - If the verifier REFUTED anything load-bearing, the affected verdicts are UNSCOREABLE — say so;
   never score around a refutation.
 - Score V1-V5. State which branch of the pre-registered rule fires (ALL PASS -> the earned-sense
   field-side statement, with the next step being the two-source composition MEASUREMENT at 4x4 —
   framed per the principal's directive as measuring what composition IS, every outcome registering;
   V1 FALSE -> the mediator-family sweep next; partial -> name exactly what blocks scoring).
 - Draft the register entry for the landing: what was computed, what the verdicts are, the
   corrections trail (from D24_AUDIT files), the next step named. The registrar lands it; you draft.
 - Write your judgment to "${REPO}/LANE_T51_A/JUDGMENT.txt".

Return: the scored rule, the branch, the draft register entry, and the ranked next step.`,
  { label: 'judge', phase: 'Judge' }
)

return { gate, measure, verify, judgment }

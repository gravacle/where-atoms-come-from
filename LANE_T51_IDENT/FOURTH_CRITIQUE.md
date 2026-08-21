# T-51 FOURTH DESIGN — COMBINED CRITIQUE (D-1/D-24 auditor + measurability skeptic) — 2026-08-21

**Design critiqued:** `FOURTH_DESIGN.md` (THE RIGID RULER — the rigidity theorem), with
`FOURTH_rigidity.py` / `FOURTH_rigidity_OUT.txt`. Posture: default REFUTED; every load-bearing
number re-run from independently written machinery (own F2 kit, own enumeration, own Hamiltonian
build, own partial trace, own cut-rank — only the sealed coordinate conventions matched so
placements name the same physics). Verification scripts: scratchpad
`fourth_verify.py`, `fourth_probe_extra.py`, plus an inline winding-restricted sweep
(session scratchpad `/private/tmp/claude-501/-Users-bgm-MB-Work/04dac55c-b963-46a5-9ea4-74afb63af51b/scratchpad/`).

## VERDICT: **RIGIDITY_CONFIRMED** — with the repairs below binding before registration.

Every computed exhibit reproduces; the three proof mechanisms are sound at their stated tier;
one precision hole in R-III(b) was found, computed, and closed in repaired form; two wording
overclaims need rescoping. No unrepaired D-1 or D-24 violation. No repair touches a number.

---

## 1. INDEPENDENT VERIFICATION — every load-bearing computation reproduces

| computation | design's value | my independent value | status |
|---|---|---|---|
| d_W class minima, L=2,3,4 (per-pair, own affine solve, standard Z̄ loops) | {L, L, 2L} | {2,2,4}, {3,3,6}, {4,4,8} | EQUAL |
| translation invariance d_W(s⊕c,s'⊕c)=d_W(s,s'), all 64 triples × 3 sizes | holds | holds | EQUAL |
| 3×2 class histograms, both placements (brute enumeration over all 2^12 Z-vectors) | as printed | bit-identical, class size 128 | EQUAL |
| (w_direct, w_enc-old, w_enc-conn), 4 placements | (1,3,4) (1,3,3) (1,4,5) (1,4,4) | same; conn minima also confirmed by hand (path-through-detector arithmetic) | EQUAL |
| F at reference λ, 4 placements | +1.130978e−4 / −8.292711e−4 / −2.605e−5 / +1.063e−4 | +1.130978e−4 / −8.292711e−4 / −2.604969e−5 / +1.063428e−4 | EQUAL (6–7 sig figs) |
| F onset first rungs | 3.995 / 2.980 / 4.988 / 3.980 | 3.995 / 2.980 / 4.988 / 3.980 (full ladders identical to 3 decimals) | EQUAL |
| Δ onset ≈ 1 = w_direct, both sectors, all placements; witness ≥ 0.990 | holds | holds | EQUAL |
| winding annex, 4 sectors at λ=0.05 | +1.131e−4 / −1.864e−5 / −1.229e−4 / +1.865e−5 | same to 4 sig figs | EQUAL |
| region tier: trace-dists, δS at λ=0 and 0.05 | 2.9e−16 / 1.000; 0 / 0; 8.86e−4, +4.56e−6 / 1.000, +1.80e−3 | same (own partial trace) | EQUAL |
| cut-ranks (R_away, R_cross) | (7, 5) | (7, 5) (own GF(2) kernel route) | EQUAL |
| mixed-span (Y-dressed) minima vs pure-Z | equal | equal | EQUAL |
| committed OUT vs fresh rerun of the lane script | 22/22 | byte-identical | EQUAL |

**Two adversarial extensions the design did not run, both strengthening it:**

- **Off-script placement** (probe stars (0,0),(0,1) vertical; source ((1,1),(1,2)) vertical, 3×3):
  fresh static computation gives old-spec 4 / connected 5; measured F onset first rung **4.987**.
  The onset–connected lock is not placement-tuned; five placements now confirm it, two of them
  from machinery that shares nothing with the design's.
- **Swept-region enumeration at 3×2** (all 495 four-edge regions, λ=0): region entropy is exactly
  unmoved on **every** region (worst |δS| = 4.4e−16 bits) — R-III(c) exhaustively confirmed at
  this venue — and the response-iff-unavoidable statement holds with **zero violations** once the
  representative class is winding-preserving (finding F-1 below). This discharges the design's
  named next step (b) at the 3×2 venue, in the repaired form.

## 2. FINDINGS

**F-1 — R-III(b) as written is false; the repaired statement is exhaustively true (REPAIR, computed).**
The statement "a region's reduced state responds to a distant write iff … every admissible
representative of the content's writer crosses the region oddly" omits the winding-sector
qualifier. Computed counterexample: at 3×2, region R = {h(0,0), h(1,0), v(2,0), v(2,1)} is
avoided by the weight-2 representative X on {v(0,1), v(1,1)} (verified: it flips exactly the two
source plaquettes), yet the region responds at trace-distance 1.000 — because that representative
crosses Z̄₁ once and therefore transports the b=+1 ground state into the **wrong winding sector**.
45 of 495 four-edge regions violate the unqualified iff; with the representative class restricted
to **winding-sector-preserving** members (class size 32 of 128), violations drop to **0/495**.
This is the same declared winding-class convention the judgment already binds to w_enc; R-III(b)
must carry it explicitly. The design's honesty clause ("mechanism-plus-instance, not a swept
enumeration") anticipated a gap; this is the gap, and it is now closed at 3×2 in repaired form.

**F-2 — The δS-crossing exhibit is a source-quantity; rescope it (REPAIR).** Verified: R_cross =
{h(2,1), h(2,0), v(2,1), v(0,1)} is **exactly the source plaquette (2,1)'s own 4-edge support**
and contains the minimal writer edge v(2,1). Consequences: (i) R-III(b)'s word "distant" is wrong
for the exhibited crossing instance — the write acts inside the region; (ii) §2 finding 3's
δS = +1.80e−3 bits "on the crossing region" is the source's own dressing entropy — a source-side
quantity — so §5.2's blanket sentence ("amplitude differentials … at places where the content is
not") is false as applied to that exhibit. Only the away value (+4.6e−6 bits) is entitled to
field-side language, and it carries no same-table noise floor. Repair: rescope finding 3 as a
source-adjacent dressing response (or exhibit a crossing region disjoint from the source support,
which needs a larger venue), and attach a measured floor to the away-δS before calling it a
response. The theorem is untouched; F, read at the probe, remains the clean field-side quantity.

**F-3 — One-sentence overclaim in §5.6 (REPAIR).** "Proving its reading is the *only* place a
field can live on this surface" — the theorem excludes the **integer tier**; it does not prove
uniqueness of the Second Lump's functional among real-valued responses (the program's own Hanging
Clock bath channel and Cert-Shadow's dynamical arm are other real-valued readings). Repair: "the
only tier", not "the only place".

**F-4 — R-I scope note (REPAIR, one clause).** The sketch quantifies over "admissible writers";
the argument and C-78's sealed metric live on the Pauli/z=0-reduced tier (the z=0 reduction was
verified exhaustively in the sealed lane at L=2,3). O-4 admissibility ([U,H]=0) admits non-Pauli
unitaries whose label action need not be a translation. The theorem is true of **the metric as
sealed**; say so, and leave general-unitary weight as the O-4-shaped scope boundary it already is.

**F-5 — Minor overstatement.** §1 R-II claims "full … weight histograms at … two 3×3 placements"
as computed exhibits; the OUT prints 3×3 minima only. Print them or write "minima exhibited".

**F-6 — Estimator declaration (presentational).** The Δ gates score the last rung, the F gates the
first; both defensible (first rung is the principled λ→0 estimator) but the round's D-8 standard
asks one declared estimator with its measured noise floor in the same table. At 3×3 the first-rung
F values sit ~2.8e−11 against an eigh floor ~1e−13 — safe by two orders, and my independent
rebuild reproduces every rung to 3 decimals, so stability is confirmed; the repair is to declare
and tabulate, not to recompute.

**Pressed and clean — import hunt (D-1).** No classical form is required or tested anywhere; grep
finds no banned comparison (Newton unnamed; no inverse-square, 1/r, G, geodesic; the only "metric"
is the program's own earned d_W). Every gate tests surface-computed integers. The negative branch
registers per C-92's directive.

**Pressed and clean — earned concepts (D-24).** d_W (C-78, sealed), w_min/w_enc conventions (C-80
sealed + the judgment's connected gate — implemented exactly as the judgment's acceptance
instance demands and reproducing it), cut-rank (C-81, sealed), enclosure/holonomy (the judged
first lane's vocabulary), "bias is mass" (verbatim in C-91's row), E1/E2/E3 (C-91's world-tier
ensembles), clause (iv′) free-energy floor (the anchor). d_gen descriptive only, per the judgment.
Lattice coordinates appear as construction labels only. **One earned-concept gap at the world
tier:** "adjacent write" in §3's named next step has no earned separation on the E1/E2/E3 venues —
the same DECLARED/UNEARNED scope the judgment pinned to the Clock's chain medium. The world-tier
computation must carry that scope tag from birth and include a separation sweep before any
field-side (as opposed to metric-response) language attaches; otherwise its positive branch reads
as write-disturb at contact — the census's program-disturb/ATI analogues the design itself cites.

**Pressed and clean — source-quantity disguise on the main readings.** F is the probe's own law
differential, read at the probe: not a source quantity. The rigidity integers are properties of
pairs, classes, and regions, never of the source's value. The two places source-side quantities
creep in are F-2's δS-crossing exhibit and, prospectively, the world-tier adjacency computation —
both repairable by scoping, neither load-bearing for the theorem.

**Certificate hygiene (D-15) — clean.** Certificates labeled and never counted (R-I invariance,
R-II histogram equality, cut-rank state-freeness); the could-fail members are real: the mixed span
could have beaten pure-Z and did not; the onsets could have tracked old-spec w_enc where
old ≠ connected and demonstrably did not; both region directions gated.

## 3. WHERE THIS SITS RELATIVE TO THE JUDGED SEQUENCE (recommendation)

This design completes the owed round without disturbing the ranking, and its static content should
not wait for the sequence. Treat it the way the judge treated Cert-Shadow's static arm: the
rigidity theorem (with F-1's winding qualifier and F-4's sealed-metric scope) is registrable now,
sharing a register entry or sitting adjacent to that arm as the geometry-currency member of the
same three-currency family mirroring C-80's zero-beyond-contact — and one piece is operationally
urgent rather than merely registrable: the 3×3 connected-gate instance (old-spec 4, connected 5,
onset 4.99 at probe (0,0),(1,0), source ((1,1),(2,1)), now confirmed by three independent
machineries and by a fifth off-script placement) must reach LANE_T51_A **before it fires**, since
V3's bracket at that placement is 5, not 4. Nothing here rivals the Second Lump — the onset-lock
verification strengthens its mechanism claim at every placement tested — and nothing displaces the
Clock or Cert-Shadow, whose channels the theorem does not touch (F-3's repair keeps that boundary
honest). The responsive-venue statement is the right handoff of the geometry-response question to
the world tier, subject to the adjacency-scope repair above; its two named next steps are
lane-sized, and my 3×2 sweep has already discharged next step (b) at the calibration venue —
the 3×3 sweep remains. Commissioning order stays: Second Lump first; this design's theorem and
V3 instance enter the record now.

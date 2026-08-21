# COMMISSIONING REFERENCES — dead workflow scripts preserved at the 2026-08-21 cutover

These are the WORKFLOW SCRIPTS of runs the session limit killed mid-flight, preserved verbatim so
the successor can re-commission without reconstructing briefs. **Workflow resume is same-session
only — these CANNOT be resumed; re-commission fresh, with the sealed on-disk artifacts as INPUT.**

- `t51a_verify_judge_reference.js` — LANE_T51_A. Gate + Measure phases COMPLETED and their sealed
  artifacts are in LANE_T51_A/ (ALL-PASS = TRUE on V1–V5, **UNVERIFIED** — see HANDOFF_2026-08-21.md
  §3.1). Re-commission ONLY the Verify and Judge phases, feeding the sealed lane as input.
- `t54_buildout_reference.js` — the URM build-out. Four families built (arrow, countlaw, writing,
  d25 — files in model/ + INTEGRATION_*.md); the `classes` FAMILIES entry is the brief for the one
  unbuilt family; the two Verify phase prompts are the verification standard.
- `t50_v003_reference.js` — T-50 V003, the terminal round. Nothing ran; the whole script is the
  commission, under the register's pre-registered boundary.

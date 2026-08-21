"""T-50 DESIGN TWO -- ADVERSARIAL MUTATION SUITE (the DONE_WHEN(b) bar).

The principal's ruling: the structural-claim check must fail at a rate of 1.000 under an
adversarial mutation suite -- not the 0.18 the old widens check measured -- and the suite
covers at least (m1) the two-signed write, (m2) a non-zero-mean residual, and (m3) a
programmed density falling with N.

WHAT "FAIL" MEANS PER MUTATION, stated before running, because the three mutations are
not the same kind of corruption:

  m1 TWO-SIGNED WRITE   corrupts the CLAIM ITSELF (occupancy no longer accumulates).
                        The claim check -- "an in-scope occupancy read ACCUMULATES" --
                        must FAIL (verdict SCREENS or INDETERMINATE, never ACCUMULATES).

  m2 NON-ZERO-MEAN      does NOT corrupt the claim: a biased in-spec residual is healthy
     RESIDUAL (mu=3e)   physics (refuter A, residual defect 8: real over-erase IS biased;
                        a claim check firing here would repeat the refuted floor's defect
                        of firing on correct physics).  What must fail at 1.000 is the
                        SENTINEL: the absolute-zero-referenced (V_t,neutral-style) verdict,
                        which false-fires under the bias -- measured beside the claim
                        check, which must be INVARIANT (rate 0.000 of moving).

  m3 DENSITY FALLING    corrupts the MEASUREMENT'S PRECONDITION.  The pipeline must
     WITH N             REFUSE at 1.000 (any reported verdict at all counts as the gate
                        failing to catch the mutation).

  m4 COMMON-MODE 0.5e   (the T-50 row's own attack, added)  same contract as m2: claim
                        check invariant, sentinel false-fires at 1.000.

  m5 DC-BIASED DATA     (added)  data at 55% ones sold as 'random' on an orientation
     SOLD AS RANDOM     track: the balance guard must route it OUT of the screening
                        clause at 1.000 (dc_loaded flag), else the screening clause
                        would misapply.

NO CHECK BELOW COMPARES AGAINST A LITERAL: every check is a verdict/guard produced by
the full pipeline on a fresh read, and every check's failing branch is exercised by at
least one suite member (D-8).  The baseline row is the suite's D-15 control in a
different configuration: same checks, unmutated reads, expected pass.
"""

import sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from d2_observable import (read_occupancy, read_occupancy_unwritten, read_orientation,
                           read_orientation_erased, estimate, estimate_naive)

N = 1 << 15
NMIN, NMAX = 128, 2048
SEEDS = 50
OUT = []

def emit(s=""):
    print(s)
    OUT.append(s)

def occ_pair(seed, f=0.5, **kw):
    rng = np.random.default_rng(seed)
    v_w, _ = read_occupancy(N, f, rng, **kw)
    kw_null = {k: kw[k] for k in ("mu", "offset", "gain") if k in kw}
    v_e = read_occupancy_unwritten(N, rng, **kw_null)
    return v_w, v_e

emit("=" * 100)
emit("T-50 DESIGN TWO -- ADVERSARIAL MUTATION SUITE -- 50 seeds per row")
emit("=" * 100)

# ---------------------------------------------------------------- baseline (D-15 control)
acc = 0
for s in range(SEEDS):
    r = estimate(*occ_pair(50000 + s), "occupancy", NMIN, NMAX)
    if r["verdict"] == "ACCUMULATES":
        acc += 1
emit(f"\nBASELINE (unmutated occupancy, the D-15 control in a different configuration):")
emit(f"  claim check ACCUMULATES: {acc}/{SEEDS} pass  (expected: pass; this row is what "
     "the mutations are measured against)")
base_pass = acc

# ---------------------------------------------------------------- m1 two-signed write
fail = 0
verds = {}
for s in range(SEEDS):
    r = estimate(*occ_pair(51000 + s, two_signed=True), "occupancy", NMIN, NMAX)
    v = r["verdict"]
    verds[v] = verds.get(v, 0) + 1
    if v != "ACCUMULATES":
        fail += 1
emit(f"\nm1 TWO-SIGNED WRITE: claim check must FAIL (not ACCUMULATES)")
emit(f"  fail rate: {fail}/{SEEDS} = {fail/SEEDS:.3f}   verdicts: {verds}")
m1_rate = fail / SEEDS

# ---------------------------------------------------------------- m2 non-zero-mean residual
claim_moved = 0
sentinel_fired = 0
for s in range(SEEDS):
    r0 = estimate(*occ_pair(52000 + s), "occupancy", NMIN, NMAX)
    rm = estimate(*occ_pair(52000 + s, mu=3.0), "occupancy", NMIN, NMAX)
    if rm["verdict"] != r0["verdict"] or \
            abs(rm.get("D", np.nan) - r0.get("D", np.nan)) > 0.05:
        claim_moved += 1
    nv = estimate_naive(*occ_pair(52000 + s, mu=3.0), NMIN, NMAX)
    if not (nv["D"] > 0.5):          # the naive differential collapses under the bias
        sentinel_fired += 1
emit(f"\nm2 NON-ZERO-MEAN RESIDUAL (mu = 3 e, in-spec, the physically-signed case):")
emit(f"  claim check moved:      {claim_moved}/{SEEDS} = {claim_moved/SEEDS:.3f}   "
     "(contract: 0.000 -- firing here would repeat refuted defect A8/B2, the floor "
     "violated by healthy biased physics)")
emit(f"  sentinel false-fired:   {sentinel_fired}/{SEEDS} = {sentinel_fired/SEEDS:.3f}   "
     "(contract: 1.000 -- the absolute-zero reference IS corrupted by the bias; "
     "this is the measured reason the self-reference is load-bearing)")
m2_claim, m2_sent = claim_moved / SEEDS, sentinel_fired / SEEDS

# ---------------------------------------------------------------- m3 density falling with N
caught = 0
for s in range(SEEDS):
    rng = np.random.default_rng(53000 + s)
    v_w, _ = read_occupancy(N, 0.5, rng, fixed_record=1024)
    v_e = read_occupancy_unwritten(N, rng)
    r = estimate(v_w, v_e, "occupancy", NMIN, NMAX)
    if r["verdict"] == "REFUSED":
        caught += 1
emit(f"\nm3 PROGRAMMED DENSITY FALLING WITH N (fixed 1024-cell record, growing block):")
emit(f"  pipeline REFUSED: {caught}/{SEEDS} = {caught/SEEDS:.3f}   (contract: 1.000; any "
     "reported verdict is the gate failing to catch the mutation)")
m3_rate = caught / SEEDS

# ---------------------------------------------------------------- m4 the row's 0.5e offset
claim_moved4 = 0
sentinel_fired4 = 0
for s in range(SEEDS):
    r0 = estimate(*occ_pair(54000 + s), "occupancy", NMIN, NMAX)
    rm = estimate(*occ_pair(54000 + s, offset=0.5), "occupancy", NMIN, NMAX)
    if rm["verdict"] != r0["verdict"] or \
            abs(rm.get("D", np.nan) - r0.get("D", np.nan)) > 0.05:
        claim_moved4 += 1
    nv = estimate_naive(*occ_pair(54000 + s, offset=0.5), NMIN, NMAX)
    if not (nv["D"] > 0.5):
        sentinel_fired4 += 1
emit(f"\nm4 COMMON-MODE OFFSET 0.5 e/cell (the T-50 row's computed attack):")
emit(f"  claim check moved:      {claim_moved4}/{SEEDS} = {claim_moved4/SEEDS:.3f}   (contract: 0.000)")
emit(f"  sentinel false-fired:   {sentinel_fired4}/{SEEDS} = {sentinel_fired4/SEEDS:.3f}   (contract: 1.000)")
m4_claim, m4_sent = claim_moved4 / SEEDS, sentinel_fired4 / SEEDS

# ---------------------------------------------------------------- m5 DC-biased data as random
routed = 0
for s in range(SEEDS):
    rng = np.random.default_rng(55000 + s)
    ss = (rng.random(N) < 0.55).astype(float) * 2 - 1
    v_w = 1.0 * ss + rng.normal(0.0, 0.1, N)
    v_e = read_orientation_erased(N, rng)
    r = estimate(v_w, v_e, "orientation", NMIN, NMAX)
    if r.get("guards", {}).get("dc_loaded", False):
        routed += 1
emit(f"\nm5 DC-BIASED DATA (55% ones) SOLD AS RANDOM on an orientation track:")
emit(f"  balance guard routed to DC clause: {routed}/{SEEDS} = {routed/SEEDS:.3f}   "
     "(contract: 1.000; the screening clause must never be applied to DC-loaded data)")
m5_rate = routed / SEEDS

# ---------------------------------------------------------------- summary
emit("\n" + "=" * 100)
emit("SUITE SUMMARY (contract -> measured)")
emit(f"  baseline claim-check pass ................ {base_pass}/{SEEDS}")
emit(f"  m1 two-signed write, claim FAILS ......... 1.000 -> {m1_rate:.3f}")
emit(f"  m2 biased residual, claim INVARIANT ...... 0.000 -> {m2_claim:.3f}")
emit(f"  m2 biased residual, sentinel FIRES ....... 1.000 -> {m2_sent:.3f}")
emit(f"  m3 falling density, pipeline REFUSES ..... 1.000 -> {m3_rate:.3f}")
emit(f"  m4 0.5e common mode, claim INVARIANT ..... 0.000 -> {m4_claim:.3f}")
emit(f"  m4 0.5e common mode, sentinel FIRES ...... 1.000 -> {m4_sent:.3f}")
emit(f"  m5 DC-loaded data, guard ROUTES .......... 1.000 -> {m5_rate:.3f}")
ok = (base_pass >= SEEDS - 2 and m1_rate == 1.0 and m2_claim == 0.0 and m2_sent == 1.0
      and m3_rate == 1.0 and m4_claim == 0.0 and m4_sent == 1.0 and m5_rate == 1.0)
emit(f"  BOOL SUITE (all contracts met): {ok}")
emit("=" * 100)

with open(os.path.join(HERE, "d2_mutation_suite.txt"), "w") as fh:
    fh.write("\n".join(OUT) + "\n")

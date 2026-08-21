"""REFUTER B / ATTACK 1 -- THE BLANK CONTROL PAGE.

A physicist's first control: run the sealed pipeline on an all-erased (never-written or
freshly erased) page presented as the WRITTEN sector.  The PREDICTION's own scope bound
is n_min >= 4(1-f)/f, unbounded as f -> 0, and the FALSIFIER text claims the mostly-0xFF
page REFUSES ("out of scope, measured at f = 0.005").  The question: what does the SEALED
ESTIMATOR -- which the design declares to be the protocol -- return at f = 0 exactly?

The code path under attack, d2_observable._core:
    if encoding == "occupancy" and 0.0 < f_hat < 1.0:   # scope guard
    if encoding == "occupancy" and 0.0 < f_hat < 1.0:   # density guard
Both guards are SKIPPED when f_hat == 0.0 exactly.  A blank page whose every cell sits
inside the 6*e_sd erased band measures f_hat = 0.0 exactly and bypasses the scope bound.
"""

import sys, os
import numpy as np

LANE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, LANE)
from d2_observable import (read_occupancy_unwritten, estimate, DELTA)

N = 1 << 15
NMIN, NMAX = 128, 2048
OUT = []

def emit(s=""):
    print(s)
    OUT.append(s)

emit("=" * 100)
emit("REFUTER B / ATTACK 1 -- BLANK PAGE PRESENTED AS THE WRITTEN SECTOR")
emit("=" * 100)

emit("\n1a. MODEL RESIDUAL (uniform{-5..5}, the design's own read model), 20 seeds")
verdicts = {}
for s in range(20):
    rng = np.random.default_rng(600000 + s)
    v_w = read_occupancy_unwritten(N, rng)      # a blank page, read as 'written'
    v_e = read_occupancy_unwritten(N, rng)      # the unwritten reference range
    r = estimate(v_w, v_e, "occupancy", NMIN, NMAX)
    v = r["verdict"]
    verdicts[v] = verdicts.get(v, 0) + 1
    if s < 3:
        emit(f"  seed {s}: verdict={v}  f_hat={r.get('f', float('nan')):.6f}  "
             f"D={r.get('D', float('nan')):+.4f}  SE={r.get('se', float('nan')):.4f}  "
             f"guards evaluated: scope&density {'SKIPPED (f==0)' if r.get('f')==0.0 else 'ran'}")
emit(f"  verdicts over 20 blank pages: {verdicts}")

emit("\n1b. GAUSSIAN RESIDUAL (a real erased distribution has tails), sd = 3 e, 20 seeds")
verdicts_g = {}
fhats = []
for s in range(20):
    rng = np.random.default_rng(610000 + s)
    v_w = rng.normal(0.0, 3.0, N)
    v_e = rng.normal(0.0, 3.0, N)
    r = estimate(v_w, v_e, "occupancy", NMIN, NMAX)
    v = r["verdict"]
    verdicts_g[v] = verdicts_g.get(v, 0) + 1
    fhats.append(r.get("f", float("nan")))
emit(f"  verdicts over 20 Gaussian blank pages: {verdicts_g}")
emit(f"  f_hat values: min {np.nanmin(fhats):.6f}  max {np.nanmax(fhats):.6f}  "
     f"(6-sigma band on 32768 cells leaves essentially no outliers)")

emit("\n1c. THE CONTRAST CASE: the same page with a HANDFUL of outlier cells (f > 0 strictly)")
rng = np.random.default_rng(620000)
v_w = read_occupancy_unwritten(N, rng).astype(float)
v_w[:3] = 100.0                                  # three programmed cells in 32768
v_e = read_occupancy_unwritten(N, rng)
r = estimate(v_w, v_e, "occupancy", NMIN, NMAX)
emit(f"  3 programmed cells (f_hat={r.get('f', float('nan')):.6f}): verdict={r['verdict']}")
if r["verdict"] == "REFUSED":
    emit(f"    reason: {r['reason'][:110]}")

emit("\nREADING. The falsifier's clause (i) fires on 'an occupancy-encoded surface whose read")
emit("passes every guard IN THE READ ITSELF ... and returns verdict SCREENS'.  At f_hat = 0.0")
emit("exactly -- the ordinary blank control page -- the sealed estimator SKIPS the scope and")
emit("density guards (both are conditioned on 0.0 < f_hat < 1.0), evaluates every remaining")
emit("guard, and the guards PASS.  Whatever verdict prints above at 1a/1b is the verdict a")
emit("stranger's run reports on a blank page.  The design's own falsifier text asserts the")
emit("mostly-0xFF page REFUSES; f = 0.005 does, f = 0.0 does not.")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "vb1_blank_page.txt"), "w") as fh:
    fh.write("\n".join(OUT) + "\n")

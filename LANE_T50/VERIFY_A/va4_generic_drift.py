"""REFUTER A -- va4: generic (one-sided) written-leg drift vs the centred case, and
the null-leg control.  Output: va4_generic_drift.txt.  Finding: one-sided drift is
misrouted to the DC clause (ACCUMULATES, dc_loaded=True) at every amplitude tested;
mean-centred drift (what an AC-coupled read channel produces) fires clause (ii)
(ACCUMULATES, dc_loaded=False, all guards passing) -- see va2/va3; null-leg drift is
properly REFUSED.  The guards protect only the leg they watch."""
import sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
from d2_observable import estimate

rows = []


def say(s):
    print(s)
    rows.append(s)


say("va4 -- GENERIC (ONE-SIDED) DRIFT: the two branches")
N = 1 << 15
for d in (0.3, 0.6, 1.0):
    r = np.random.default_rng(410 + int(d * 10))
    s_ = (r.integers(0, 2, N) * 2 - 1).astype(float)
    v_w = s_ + r.normal(0, 0.1, N) + d * np.linspace(0, 1, N)
    v_e = (r.integers(0, 2, N) * 2 - 1).astype(float) + r.normal(0, 0.1, N)
    res = estimate(v_w, v_e, "orientation")
    tag = ("REFUSED " + res["reason"].split(":")[0]) if res.get("refused") else \
        f"D={res['D']:+.4f} {res['verdict']} dc_loaded={res['guards'].get('dc_loaded')}"
    say(f"  one-sided ramp 0..{d:.1f}: {tag}")
say("  (an AC-coupled / DC-servo'd read channel nulls the track mean in hardware,")
say("   turning generic drift into the centred case that fires clause (ii))")
for d in (0.3, 1.0):
    r = np.random.default_rng(430 + int(d * 10))
    s_ = (r.integers(0, 2, N) * 2 - 1).astype(float)
    v_w = s_ + r.normal(0, 0.1, N)
    v_e = (r.integers(0, 2, N) * 2 - 1).astype(float) + r.normal(0, 0.1, N) \
        + d * np.linspace(0, 1, N)
    res = estimate(v_w, v_e, "orientation")
    tag = ("REFUSED " + res["reason"].split(":")[0]) if res.get("refused") else \
        f"D={res['D']:+.4f} {res['verdict']}"
    say(f"  null-leg-only drift 0..{d:.1f}: {tag}   (the guarded leg: refusal expected)")

with open(os.path.join(HERE, "va4_generic_drift.txt"), "w") as fh:
    fh.write("\n".join(rows) + "\n")

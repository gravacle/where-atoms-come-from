"""REFUTER B / ATTACK 2 ROBUSTNESS -- drift SHAPE and track LENGTH.

The linear ramp is the most ordinary instrument drift there is (tip height, temperature).
A ramp is NOT zero-mean about the track unless centred; a centred ramp (drift crossing
zero mid-scan) has zero mean and defeats the dc-balance guard exactly like the sine.
Sweep shape {centred ramp, 1-period sine} x amplitude x track length.
Reference region flat (its own part of the frame).  Healthy verdict: SCREENS.
"""

import sys, os
import numpy as np

LANE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, LANE)
from d2_observable import read_orientation, read_orientation_erased, estimate

NMIN, NMAX = 128, 2048
OUT = []

def emit(s=""):
    print(s)
    OUT.append(s)

def row(label, r):
    if r.get("refused"):
        emit(f"  {label:<46s} REFUSED: {r['reason'][:55]}")
    else:
        emit(f"  {label:<46s} beta_w={r['beta_w']:+.3f} D={r['D']:+.4f} SE={r['se']:.4f} "
             f"dc={r['guards'].get('dc_loaded')}  --> {r['verdict']}")

emit("=" * 100)
emit("REFUTER B / DRIFT SHAPES x LENGTHS (balanced data, flat reference region)")
emit("=" * 100)
for NO in (1 << 15, 1 << 17):
    emit(f"\n  track length N = {NO}")
    for shape in ("ramp", "sine"):
        for A in (0.08, 0.125, 0.20):
            rng = np.random.default_rng(900000 + NO % 997 + int(A * 1000))
            v_w = read_orientation(NO, rng, kind="random", m=1.0, read_sd=0.1)
            v_e = read_orientation_erased(NO, rng, m=1.0, read_sd=0.1)
            i = np.arange(NO)
            if shape == "ramp":
                v_w = v_w + A * (2.0 * i / NO - 1.0)      # centred linear drift, zero mean
            else:
                v_w = v_w + A * np.sin(2 * np.pi * i / NO)
            r = estimate(v_w, v_e, "orientation", NMIN, NMAX)
            row(f"{shape}, A={A:.3f}", r)

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "vb5_drift_shapes.txt"), "w") as fh:
    fh.write("\n".join(OUT) + "\n")

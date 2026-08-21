"""REFUTER B / ATTACK 2 REFINED -- the interleave guard sees only the null leg's
HALF-MEAN DIFFERENCE.  A background even-symmetric about the region midpoint (any even
harmonic: bow that rises and falls, scan returning to focus) has EQUAL half-means and
sails through; a background confined to the written track's part of the frame never
touches the guard at all.  Both are ordinary states of one real scan frame.

Sweep: written track carries a one-period (odd) or two-period (even) zero-mean bow of
amplitude A_w; the AC-erased reference region carries an even-symmetric bow of amplitude
A_e in {0, 0.05, 0.10}.  DC-balanced random data everywhere -- the healthy verdict is
SCREENS.  Every ACCUMULATES with all guards green is falsifier clause (ii) firing on
ordinary physics.
"""

import sys, os
import numpy as np

LANE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, LANE)
from d2_observable import read_orientation, read_orientation_erased, estimate

N = 1 << 17
NMIN, NMAX = 128, 2048
OUT = []

def emit(s=""):
    print(s)
    OUT.append(s)

def pair(seed, A_w, A_e, periods_w=1.0):
    rng = np.random.default_rng(seed)
    v_w = read_orientation(N, rng, kind="random", m=1.0, read_sd=0.1)
    v_e = read_orientation_erased(N, rng, m=1.0, read_sd=0.1)
    i = np.arange(N)
    v_w = v_w + A_w * np.sin(2 * np.pi * periods_w * i / N)
    # even-symmetric bow on the reference region: equal half-means by construction
    v_e = v_e + A_e * np.sin(2 * np.pi * 2.0 * i / N)
    return v_w, v_e

def row(label, r):
    g = r.get("guards", {})
    if r.get("refused"):
        emit(f"  {label:<48s} REFUSED: {r['reason'][:60]}")
    else:
        emit(f"  {label:<48s} beta_w={r['beta_w']:+.3f} beta_e={r['beta_e']:+.3f} "
             f"D={r['D']:+.4f} SE={r['se']:.4f} dc={g.get('dc_loaded')} "
             f"null={g.get('null_leg_ok')} intl={g.get('interleave_ok')} "
             f"res={g.get('ref_resolution_ok')}  --> {r['verdict']}")

emit("=" * 112)
emit("REFUTER B / ATTACK 2 REFINED -- WRITTEN-TRACK BOW vs A HALF-BALANCED (EVEN) REFERENCE BOW")
emit("=" * 112)

emit("\n3a. Reference region FLAT (bow lives on the written track's part of the frame)")
for A_w in (0.05, 0.10, 0.125, 0.15, 0.20):
    row(f"A_w={A_w:.3f} (1 period), A_e=0", estimate(*pair(730000 + int(A_w * 1000), A_w, 0.0),
                                                     "orientation", NMIN, NMAX))

emit("\n3b. Reference region carries an EVEN bow (equal half-means -- interleave guard blind)")
for A_w, A_e in ((0.10, 0.05), (0.15, 0.05), (0.15, 0.10), (0.20, 0.10)):
    row(f"A_w={A_w:.3f}, A_e={A_e:.2f} (even)",
        estimate(*pair(740000 + int(A_w * 1000) + int(A_e * 100), A_w, A_e),
                 "orientation", NMIN, NMAX))

emit("\n3c. Written bow also EVEN (2 periods; nothing changes -- no guard sees the written leg)")
for A_w in (0.125, 0.15):
    row(f"A_w={A_w:.3f} (2 periods), A_e=0",
        estimate(*pair(750000 + int(A_w * 1000), A_w, 0.0, periods_w=2.0),
                 "orientation", NMIN, NMAX))

emit("\n3d. RATE at the working point A_w=0.15 (1 period), A_e=0.05 (even), 50 seeds")
verds = {}
ds = []
for s in range(50):
    r = estimate(*pair(760000 + s, 0.15, 0.05), "orientation", NMIN, NMAX)
    verds[r["verdict"]] = verds.get(r["verdict"], 0) + 1
    if not r.get("refused"):
        ds.append(r["D"])
emit(f"  verdicts: {verds}")
if ds:
    ds = np.array(ds)
    emit(f"  D over non-refused reads: mean {ds.mean():+.4f}, min {ds.min():+.4f}, "
         f"max {ds.max():+.4f}")

emit("\n3e. And the SAME even-null trick at smaller written bows -- where does the flip begin?")
for A_w in (0.04, 0.06, 0.08):
    row(f"A_w={A_w:.3f} (1 period), A_e=0.02 (even)",
        estimate(*pair(770000 + int(A_w * 1000), A_w, 0.02), "orientation", NMIN, NMAX))

emit("\nREADING.  The guards watch the null leg's half-mean difference and its slope, and the")
emit("written leg's MEAN (dc balance).  A zero-mean smooth background on the WRITTEN track is")
emit("watched by nothing.  Whatever ACCUMULATES verdicts print above were produced by a")
emit("DC-balanced data track (true state: screening) read with an ordinary scan bow.")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "vb3_drift_refined.txt"), "w") as fh:
    fh.write("\n".join(OUT) + "\n")

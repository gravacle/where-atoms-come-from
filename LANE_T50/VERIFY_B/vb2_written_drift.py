"""REFUTER B / ATTACK 2 -- ZERO-MEAN BACKGROUND DRIFT ON THE WRITTEN TRACK.

The falsifier text claims drifting reads are covered by 'interleave / null-leg guards'.
Both of those guards look ONLY at the reference/null population (the AC-erased region).
Nothing in the pipeline inspects the WRITTEN track's own smooth background.  A scanning
read (MFM / scanning Hall / spin-stand with imperfect leveling) puts a smooth additive
background on the map; scan bow and tip-height drift are ZERO-MEAN after leveling and
their amplitude differs from region to region of the same frame.

The attack: DC-balanced random data track (the SCREENS case, healthy physics), plus a
zero-mean sinusoidal background A*sin(2*pi*i/N) on the written region.  The erased
region carries its own background with amplitude A_e <= A (region-to-region asymmetry
of the same scan).  Sweep A and the asymmetry; record every guard and the verdict.

Numbers: block-sum drift power ~ n^2 A^2/2 vs signal power ~ n*(m^2+sd^2); parity at
n* = 2(1+sd^2)/A^2.  A = 0.125 (12.5%% of grain amplitude) puts n* at the ladder foot.
The dc-balance guard sees mean(u_w) ~ 0 (the background is zero-mean by construction),
the interleave and null-leg guards see only the erased region.  Predicted hole: false
ACCUMULATES with every guard green -- falsifier clause (ii) fired by an ordinary scan.
"""

import sys, os
import numpy as np

LANE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, LANE)
from d2_observable import read_orientation, read_orientation_erased, estimate

N = 1 << 17          # the design's own V4 track length
NMIN, NMAX = 128, 2048
OUT = []

def emit(s=""):
    print(s)
    OUT.append(s)

def drifted_pair(seed, A_w, A_e, periods_w=1.0, periods_e=1.0):
    rng = np.random.default_rng(seed)
    v_w = read_orientation(N, rng, kind="random", m=1.0, read_sd=0.1)
    v_e = read_orientation_erased(N, rng, m=1.0, read_sd=0.1)
    i = np.arange(N)
    v_w = v_w + A_w * np.sin(2 * np.pi * periods_w * i / N)   # zero-mean background
    v_e = v_e + A_e * np.sin(2 * np.pi * periods_e * i / N)
    return v_w, v_e

def row(label, r):
    g = r.get("guards", {})
    if r.get("refused"):
        emit(f"  {label:<52s} REFUSED: {r['reason'][:70]}")
    else:
        emit(f"  {label:<52s} beta_w={r['beta_w']:+.3f} beta_e={r['beta_e']:+.3f} "
             f"D={r['D']:+.4f} SE={r['se']:.4f} dc_loaded={g.get('dc_loaded')} "
             f"null_ok={g.get('null_leg_ok')} intl_ok={g.get('interleave_ok')}  "
             f"--> {r['verdict']}")

emit("=" * 110)
emit("REFUTER B / ATTACK 2 -- ZERO-MEAN WRITTEN-TRACK BACKGROUND, DC-BALANCED DATA (healthy = SCREENS)")
emit("=" * 110)

emit("\n2a. BASELINE (no background): the healthy verdict")
row("A_w=0.00, A_e=0.00", estimate(*drifted_pair(700001, 0.0, 0.0), "orientation", NMIN, NMAX))

emit("\n2b. SYMMETRIC background (same amplitude both regions -- the best case for the design)")
for A in (0.05, 0.10, 0.15, 0.25):
    row(f"A_w={A:.2f}, A_e={A:.2f}",
        estimate(*drifted_pair(700010 + int(A * 100), A, A), "orientation", NMIN, NMAX))

emit("\n2c. ASYMMETRIC background (written region carries more of the frame's bow -- ordinary)")
for A_w, A_e in ((0.10, 0.05), (0.125, 0.05), (0.15, 0.05), (0.15, 0.075), (0.20, 0.10),
                 (0.25, 0.10)):
    row(f"A_w={A_w:.3f}, A_e={A_e:.3f}",
        estimate(*drifted_pair(700100 + int(A_w * 1000), A_w, A_e), "orientation", NMIN, NMAX))

emit("\n2d. HOW OFTEN, over 50 seeds, at the ordinary operating point A_w=0.15, A_e=0.05")
verds = {}
for s in range(50):
    r = estimate(*drifted_pair(710000 + s, 0.15, 0.05), "orientation", NMIN, NMAX)
    verds[r["verdict"]] = verds.get(r["verdict"], 0) + 1
emit(f"  verdicts: {verds}")
emit("  every ACCUMULATES above is falsifier clause (ii) firing on a DC-balanced, correctly")
emit("  screening track whose only sin is a 15% zero-mean scan background -- with every")
emit("  guard the pipeline owns showing green.")

emit("\n2e. CONTROL: the same background amplitudes on the NULL leg alone (the guarded direction)")
for A_e in (0.15, 0.25):
    row(f"A_w=0.000, A_e={A_e:.2f}",
        estimate(*drifted_pair(720000 + int(A_e * 100), 0.0, A_e), "orientation", NMIN, NMAX))
emit("  (when the drift sits on the region the guards actually watch, the pipeline behaves;")
emit("   the written track is the unwatched direction)")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "vb2_written_drift.txt"), "w") as fh:
    fh.write("\n".join(OUT) + "\n")

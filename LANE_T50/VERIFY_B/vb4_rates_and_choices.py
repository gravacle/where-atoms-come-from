"""REFUTER B / ATTACKS 3-5.

4a. RATE of the written-bow false ACCUMULATES (clause ii) at the ordinary operating
    point: A_w = 0.125, reference region flat, 50 seeds.
4b. ATTACK 3 -- DIFFERENTIAL INTER-REGION OFFSET (orientation).  The reference is a
    DIFFERENT region of the scan; a DC background difference delta between the two
    regions loads mean(u_w) directly.  The dc-balance guard fires at z > 4, i.e. at
    delta = 4*sigma/sqrt(N) -- 1.1% of grain amplitude at N = 2^17, SHRINKING with N.
    Consequence measured here: an honestly DC-balanced track is routed to the DC clause
    (dc_loaded=True), where the design predicts D = 1 -- and no falsifier clause exists.
    Measure the misrouting threshold and what D the misrouted track actually reports.
4c. ATTACK 4 -- READER FREE CHOICES, occupancy: read quantization (the V_t sweep step),
    program-disturb shift of the written sector's erased band, sector length, n_max.
    Honest choices must not flip a verdict; refusals are acceptable outcomes.
4d. ATTACK 5 -- bit-run structure on orientation (real bits span multiple grains/read
    samples): runs of length L in {8, 32, 64, 128} at the sealed ladder.
"""

import sys, os
import numpy as np

LANE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, LANE)
from d2_observable import (read_occupancy, read_occupancy_unwritten, read_orientation,
                           read_orientation_erased, estimate, N_E, DELTA)

NMIN, NMAX = 128, 2048
OUT = []

def emit(s=""):
    print(s)
    OUT.append(s)

def row(label, r):
    g = r.get("guards", {})
    if r.get("refused"):
        emit(f"  {label:<52s} REFUSED: {r['reason'][:60]}")
    else:
        emit(f"  {label:<52s} beta_w={r['beta_w']:+.3f} beta_e={r['beta_e']:+.3f} "
             f"D={r['D']:+.4f} SE={r['se']:.4f} dc={g.get('dc_loaded')}  --> {r['verdict']}")

emit("=" * 112)
emit("REFUTER B / ATTACKS 3-5 -- RATES, MISROUTING, AND THE READER'S FREE CHOICES")
emit("=" * 112)

# ---------------------------------------------------------------- 4a rate of clause-(ii) false fire
emit("\n4a. RATE: written bow A_w=0.125 (1 period), reference region flat, 50 seeds, N=2^17")
NO = 1 << 17
verds = {}
ds = []
for s in range(50):
    rng = np.random.default_rng(800000 + s)
    v_w = read_orientation(NO, rng, kind="random", m=1.0, read_sd=0.1)
    v_e = read_orientation_erased(NO, rng, m=1.0, read_sd=0.1)
    v_w = v_w + 0.125 * np.sin(2 * np.pi * np.arange(NO) / NO)
    r = estimate(v_w, v_e, "orientation", NMIN, NMAX)
    verds[r["verdict"]] = verds.get(r["verdict"], 0) + 1
    if not r.get("refused"):
        ds.append(r["D"])
ds = np.array(ds)
emit(f"  verdicts: {verds}")
emit(f"  D: mean {ds.mean():+.4f}  min {ds.min():+.4f}  max {ds.max():+.4f}")
emit("  clause (ii) fire rate = ACCUMULATES / 50 on a healthy DC-balanced (screening) track")

# ---------------------------------------------------------------- 4b differential offset misrouting
emit("\n4b. DIFFERENTIAL INTER-REGION OFFSET delta on an honestly DC-balanced track (N=2^17)")
emit("    guard threshold: z=4 at delta = 4*sigma/sqrt(N) = "
     f"{4 * 1.005 / np.sqrt(NO):.4f} of grain amplitude")
for delta in (0.005, 0.010, 0.012, 0.020, 0.050):
    rng = np.random.default_rng(810000 + int(delta * 10000))
    v_w = read_orientation(NO, rng, kind="random", m=1.0, read_sd=0.1) + delta
    v_e = read_orientation_erased(NO, rng, m=1.0, read_sd=0.1)
    r = estimate(v_w, v_e, "orientation", NMIN, NMAX)
    row(f"delta={delta:.3f} (background mismatch, balanced data)", r)
emit("  every dc=True row above is an honestly balanced track ROUTED OUT of the screening")
emit("  clause into the DC clause, whose prediction (D=1) it then contradicts (D~0) --")
emit("  and the falsifier set contains NO clause for the DC route.  Also note the guard")
emit("  threshold TIGHTENS as 1/sqrt(N): the better the read, the more likely the misroute.")

# ---------------------------------------------------------------- 4c occupancy free choices
emit("\n4c. OCCUPANCY -- READER FREE CHOICES (honest choices, f=0.5 page, N=2^15)")
N = 1 << 15

emit("  -- quantization: V_t sweep step q, in units of the 10 e erased spread (e_sd~3 e)")
for q in (1.0, 3.0, 6.0, 12.0):
    rng = np.random.default_rng(820000 + int(q * 10))
    v_w, _ = read_occupancy(N, 0.5, rng)
    v_e = read_occupancy_unwritten(N, rng)
    vq_w = np.round(v_w / q) * q
    vq_e = np.round(v_e / q) * q
    r = estimate(vq_w, vq_e, "occupancy", NMIN, NMAX)
    row(f"sweep step q={q:.0f} e", r)

emit("  -- floor-clipped erased read (sweep cannot go below the erased nominal): v = max(v, 0)")
rng = np.random.default_rng(821000)
v_w, _ = read_occupancy(N, 0.5, rng)
v_e = read_occupancy_unwritten(N, rng)
r = estimate(np.maximum(v_w, 0.0), np.maximum(v_e, 0.0), "occupancy", NMIN, NMAX)
row("half-censored erased distribution", r)

emit("  -- program disturb: written sector's erased cells shifted up by d e (unwritten range clean)")
for d in (3.0, 9.0, 30.0):
    rng = np.random.default_rng(822000 + int(d * 10))
    data = (rng.random(N) < 0.5).astype(int)
    resid = rng.integers(-DELTA, DELTA + 1, N).astype(float) + d      # disturbed erased cells
    prog = N_E + rng.normal(0.0, 2.0, N)
    v_w = np.where(data == 1, prog, resid)
    v_e = read_occupancy_unwritten(N, rng)                            # clean unwritten range
    r = estimate(v_w, v_e, "occupancy", NMIN, NMAX)
    row(f"disturb shift d={d:.0f} e", r)

emit("  -- sector length and ladder top (the reader's own sizes)")
for (Nn, nmx) in ((1 << 13, 512), (1 << 15, 2048), (1 << 17, 8192)):
    rng = np.random.default_rng(823000 + Nn % 997)
    v_w, _ = read_occupancy(Nn, 0.5, rng)
    v_e = read_occupancy_unwritten(Nn, rng)
    r = estimate(v_w, v_e, "occupancy", 128, nmx)
    row(f"N={Nn}, ladder 128..{nmx}", r)

# ---------------------------------------------------------------- 4d bit-run structure
emit("\n4d. ORIENTATION -- REAL BITS SPAN RUNS OF GRAINS/SAMPLES (balanced data, runs of L)")
for L in (8, 32, 64, 128):
    rng = np.random.default_rng(830000 + L)
    nbits = NO // L + 1
    bits = (rng.integers(0, 2, nbits) * 2 - 1).astype(float)
    s = np.repeat(bits, L)[:NO]
    v_w = 1.0 * s + rng.normal(0.0, 0.1, NO)
    v_e = read_orientation_erased(NO, rng, m=1.0, read_sd=0.1)
    r = estimate(v_w, v_e, "orientation", NMIN, NMAX)
    row(f"run length L={L} (bit = {L} read samples)", r)
emit("  (the sealed ladder starts at 128 samples; the verdict must hold for honest")
emit("   oversampling choices L << n_min and degrade to refusal/indeterminate, never to a")
emit("   false ACCUMULATES, as L approaches n_min)")

with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "vb4_rates_and_choices.txt"), "w") as fh:
    fh.write("\n".join(OUT) + "\n")

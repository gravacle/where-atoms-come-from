"""REFUTER A -- va1: REBUILD the design's headline model-side numbers with my own
estimator (va_estimator.py), my own read models, my own seeds; then cross-check my
estimator against the sealed pipeline on identical arrays; then the sealed seed-11
continuity from model/geometry.py directly."""
import sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))                       # LANE_T50 (sealed pipeline)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(HERE)), "model"))
from va_estimator import run, core, occ, ori
import d2_observable as D2

N = 1 << 15
rows = []


def say(s):
    print(s)
    rows.append(s)


say("=" * 96)
say("va1 -- REBUILD OF THE DESIGN'S OPERATING POINTS (independent code, independent seeds)")
say("=" * 96)

# ---- R1: occupancy ensemble, my code, my seeds (design claims mean 1.009 sd 0.100,
#          ACCUMULATES with 0 false SCREENS over 200)
ds, ses, verds = [], [], {}
for s in range(200):
    rng = np.random.default_rng(777000 + s)
    v_w, v_e = occ(N, 0.5, rng)
    r = run(v_w, v_e, "occupancy")
    verds[r["verdict"]] = verds.get(r["verdict"], 0) + 1
    if not r["refused"]:
        ds.append(r["D"]); ses.append(r["se"])
ds, ses = np.array(ds), np.array(ses)
say(f"R1 occupancy f=0.5, 200 fresh reads: D mean {ds.mean():.4f} sd {ds.std():.4f} "
    f"(design: 1.009 / 0.100); mean jk SE {np.nanmean(ses):.4f} calib "
    f"{ds.std()/np.nanmean(ses):.2f} (design 0.90); verdicts {verds}")

# ---- R1b: occupancy across f, my seeds (design: all f in scope ACCUMULATES)
for f in (0.01, 0.05, 0.1, 0.3, 0.7, 0.9, 1.0):
    rng = np.random.default_rng(int(1000 * f) + 88)
    v_w, v_e = occ(N, f, rng)
    r = run(v_w, v_e, "occupancy")
    say(f"R1b f={f:.2f}: " + (f"REFUSED ({r['why']})" if r["refused"] else
        f"beta_w={r['beta_w']:.4f} beta_e={r['beta_e']:.4f} D={r['D']:+.4f} "
        f"SE={r['se']:.4f} {r['verdict']}"))

# ---- R2: orientation ensemble (design: mean 0.018 sd 0.161, 0 false ACCUMULATES,
#          calib 1.00)
ds, ses, verds, nacc = [], [], {}, 0
for s in range(200):
    rng = np.random.default_rng(888000 + s)
    v_w, v_e = ori(N, rng, kind="random")
    r = run(v_w, v_e, "orientation")
    verds[r["verdict"]] = verds.get(r["verdict"], 0) + 1
    if not r["refused"]:
        ds.append(r["D"]); ses.append(r["se"])
        if r["verdict"] == "ACCUMULATES":
            nacc += 1
ds, ses = np.array(ds), np.array(ses)
say(f"R2 orientation random, 200 reads: D mean {ds.mean():.4f} sd {ds.std():.4f} "
    f"(design: 0.018 / 0.161); calib {ds.std()/np.nanmean(ses):.2f} (design 1.00); "
    f"false ACCUMULATES {nacc}/200; verdicts {verds}")

# ---- R2b: DC-saturated (design: D = 1.020 ACCUMULATES, dc_loaded True)
rng = np.random.default_rng(999)
v_w, v_e = ori(1 << 17, rng, kind="dc")
r = run(v_w, v_e, "orientation")
say(f"R2b DC-saturated 2^17: D={r['D']:+.4f} {r['verdict']} dc_loaded={r['g'].get('dc_loaded')}")

# ---- R3: cross-check my estimator against the sealed pipeline on IDENTICAL arrays
say("\nR3 CROSS-CHECK (same arrays through both implementations):")
worst = 0.0
for s in range(20):
    rng = np.random.default_rng(661000 + s)
    v_w, v_e = occ(N, 0.5, rng)
    a = run(v_w, v_e, "occupancy")
    b = D2.estimate(v_w, v_e, "occupancy")
    dd = abs(a["D"] - b["D"])
    worst = max(worst, dd)
    if a["verdict"] != b["verdict"]:
        say(f"  seed {s}: VERDICT MISMATCH {a['verdict']} vs {b['verdict']}")
for s in range(20):
    rng = np.random.default_rng(662000 + s)
    v_w, v_e = ori(N, rng, kind="random")
    a = run(v_w, v_e, "orientation")
    b = D2.estimate(v_w, v_e, "orientation")
    if not (a["refused"] or b["refused"]):
        worst = max(worst, abs(a["D"] - b["D"]))
    if a["verdict"] != b["verdict"]:
        say(f"  ori seed {s}: VERDICT MISMATCH {a['verdict']} vs {b['verdict']}")
say(f"  max |D_mine - D_sealed| over 40 paired reads: {worst:.2e}")

# ---- R4: sealed seed-11 continuity, my estimator, geometry.py read-only
say("\nR4 SEALED SEED-11 CONTINUITY (design: 993/1000 D>0.5, mean 0.957, 7 refusals):")
import geometry as GE
pats = GE.occupancy_patterns()
err = pats["unwritten_e"].astype(float)
rngc = np.random.default_rng(424242)          # my own null draws, declared
ds, nref = [], 0
for w in pats["written"]:
    v_w = np.where(w == 1, 100.0, err)
    v_e = rngc.integers(-5, 6, 1000).astype(float)
    r = core(v_w, v_e, "occupancy", 4, 64)
    if r["refused"]:
        nref += 1
    else:
        ds.append(r["D"])
ds = np.array(ds)
say(f"  my numbers: {int((ds > 0.5).sum())}/1000 pages D>0.5, mean {ds.mean():.4f}, "
    f"min {ds.min():.4f}, max {ds.max():.4f}, refusals {nref}")

# ---- R5: affine invariance is EXACT (my code, algebraic claim of the design)
say("\nR5 AFFINE INVARIANCE (v -> g(v+b), both sectors), my draws:")
rng = np.random.default_rng(555)
v_w, v_e = occ(N, 0.5, rng)
r0 = run(v_w, v_e, "occupancy")
worst = 0.0
for g_, b_ in ((1, 0.5), (1, -50), (2, 7), (0.5, -0.5), (3.7, 123.4)):
    r = run(g_ * (v_w + b_), g_ * (v_e + b_), "occupancy")
    worst = max(worst, abs(r["D"] - r0["D"]))
say(f"  max |dD| over 5 affine maps: {worst:.2e}  (design claims exact cancellation)")

# ---- R6: coverage of the stated PREDICTION 'D = 1 within 2 jackknife SE' per read
say("\nR6 PER-READ COVERAGE OF THE STATED PREDICTION (not the falsifier):")
viol_o, tot_o = 0, 0
for s in range(200):
    rng = np.random.default_rng(777000 + s)
    v_w, v_e = occ(N, 0.5, rng)
    r = run(v_w, v_e, "occupancy")
    if not r["refused"] and np.isfinite(r["se"]):
        tot_o += 1
        if abs(r["D"] - 1.0) > 2 * r["se"]:
            viol_o += 1
viol_r, tot_r = 0, 0
for s in range(200):
    rng = np.random.default_rng(888000 + s)
    v_w, v_e = ori(N, rng, kind="random")
    r = run(v_w, v_e, "orientation")
    if not r["refused"] and np.isfinite(r["se"]):
        tot_r += 1
        if abs(r["D"] - 0.0) > 2 * r["se"]:
            viol_r += 1
say(f"  occupancy: |D-1| > 2SE on {viol_o}/{tot_o} honest in-scope reads "
    f"({100*viol_o/max(tot_o,1):.1f}%)")
say(f"  orientation: |D-0| > 2SE on {viol_r}/{tot_r} honest in-scope reads "
    f"({100*viol_r/max(tot_r,1):.1f}%)")
say("  the falsifier's own bar (SCREENS/ACCUMULATES against the 0.5 midpoint) is what")
say("  protects the row; the prediction SENTENCE 'D = 1 within 2 jackknife SE' for EVERY")
say("  in-scope read is a per-read coverage claim and the rate above is its violation rate")

with open(os.path.join(HERE, "va1_rebuild.txt"), "w") as fh:
    fh.write("\n".join(rows) + "\n")

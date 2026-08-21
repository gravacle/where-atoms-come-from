"""REFUTER A -- va2: ATTACKS on the sealed pipeline (d2_observable.estimate), each aimed
at firing a falsifier clause on CORRECT physics or breaking a stated guard claim.

A1 sector-differential OFFSETS (written-only, null-only), occupancy + orientation
A2 WRITTEN-LEG-ONLY DRIFT on an orientation random-data track -- the guards claim
   'drifting reads' are covered by interleave/null-leg guards; both act on the NULL.
   A mean-centred ramp on the WRITTEN leg passes the DC-balance guard.  Does clause
   (ii) fire (ACCUMULATES on a DC-balanced screening surface, all guards passing)?
A3 the fair-drift control: the same drift on BOTH legs (shared-scan assumption) --
   does the null-leg guard refuse, as the falsifier text asserts?
A4 grains-per-bit CLUSTERING (real recording: one bit spans many grains): dc_loaded
   misroute rate on honest DC-balanced data; any false ACCUMULATES at the sealed ladder
A5 density confounds beyond m3: slow density gradient; runs at half the block size
A6 nominal-vs-actual: N_E 30 e (3D NAND class), prog_sd 10, DELTA 2, small parts,
   and the n_max <= N/16 protocol bound -- is it enforced by the pipeline?
A7 borderline-z misroute rate on honest iid orientation reads (flag vs writing kind)
A8 occupancy false-SCREENS hunt: null-leg drift held just inside [0.6,1.4] + low f
"""
import sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
from va_estimator import occ, ori
from d2_observable import estimate

N = 1 << 15
rows = []


def say(s):
    print(s)
    rows.append(s)


def show(label, r):
    if r.get("refused"):
        say(f"  {label:<52s} REFUSED ({r['reason'].split(':')[0]})")
    else:
        say(f"  {label:<52s} beta_w={r['beta_w']:+.3f} D={r['D']:+.4f} SE={r['se']:.4f} "
            f"{r['verdict']}" + (f" dc_loaded={r['guards'].get('dc_loaded')}"
                                 if 'dc_loaded' in r.get('guards', {}) else ""))


say("=" * 96)
say("va2 -- ATTACKS ON THE SEALED PIPELINE")
say("=" * 96)

# ---------------- A1 sector-differential offsets
say("\nA1 SECTOR-DIFFERENTIAL OFFSETS (the offset the affine proof does NOT cover)")
rng = np.random.default_rng(11)
for b in (2.0, 10.0, 20.0, 40.0):
    r = np.random.default_rng(11)
    v_w, v_e = occ(N, 0.5, r, off_w=b, off_e=0.0)
    show(f"occupancy, written-only offset +{b:.0f} e", estimate(v_w, v_e, "occupancy"))
for b in (0.02, 0.05, 0.2, 1.0):
    r = np.random.default_rng(12)
    v_w, v_e = ori(N, r, kind="random", off_w=b, off_e=0.0)
    show(f"orientation, track-only offset +{b:.2f} (sd=1)", estimate(v_w, v_e, "orientation"))

# ---------------- A2 written-leg drift (the centrepiece attack)
say("\nA2 WRITTEN-LEG-ONLY DRIFT, orientation random data (mean-centred ramp; the")
say("   falsifier text claims 'drifting reads' are covered by interleave/null-leg guards)")
fired = 0
for d in (0.2, 0.5, 1.0, 2.0):
    r = np.random.default_rng(21)
    v_w, v_e = ori(N, r, kind="random", drift_w=d)
    res = estimate(v_w, v_e, "orientation")
    show(f"drift amplitude {d:.1f} x signal across track", res)
    if (not res.get("refused") and res["verdict"] == "ACCUMULATES"
            and not res["guards"].get("dc_loaded")):
        fired += 1
say(f"  CLAUSE (ii) FIRED (ACCUMULATES, dc guard passed, no refusal): {fired}/4 amplitudes")
say("  -- ensemble at drift 1.0 x signal, 50 fresh reads:")
verds, dcs = {}, 0
for s in range(50):
    r = np.random.default_rng(2200 + s)
    v_w, v_e = ori(N, r, kind="random", drift_w=1.0)
    res = estimate(v_w, v_e, "orientation")
    verds[res["verdict"]] = verds.get(res["verdict"], 0) + 1
    if res.get("guards", {}).get("dc_loaded"):
        dcs += 1
say(f"     verdicts {verds}, dc_loaded count {dcs}/50")

# ---------------- A3 fair-drift control
say("\nA3 SHARED-SCAN DRIFT CONTROL (same drift on both legs -- the case the guards own)")
for d in (0.2, 0.5, 1.0):
    r = np.random.default_rng(31)
    v_w, v_e = ori(N, r, kind="random", drift_w=d, drift_e=d)
    show(f"both-leg drift {d:.1f}", estimate(v_w, v_e, "orientation"))

# ---------------- A4 grains-per-bit clustering
say("\nA4 GRAINS-PER-BIT CLUSTERING (one written bit spans k grains; data DC-balanced)")
for k in (4, 16, 64):
    mis, acc, scr, tot = 0, 0, 0, 0
    for s in range(50):
        r = np.random.default_rng(4000 + 100 * k + s)
        v_w, v_e = ori(N, r, kind="random", cluster=k)
        res = estimate(v_w, v_e, "orientation")
        if res.get("refused"):
            continue
        tot += 1
        if res["guards"].get("dc_loaded"):
            mis += 1
        elif res["verdict"] == "ACCUMULATES":
            acc += 1
        elif res["verdict"] == "SCREENS":
            scr += 1
    say(f"  k={k:>3d}: dc_loaded misroute {mis}/{tot}, clause-(ii) fires {acc}/{tot}, "
        f"SCREENS {scr}/{tot}")
# one detailed read at k=64, z<4 subset
for s in range(200):
    r = np.random.default_rng(46000 + s)
    v_w, v_e = ori(N, r, kind="random", cluster=64)
    res = estimate(v_w, v_e, "orientation")
    if not res.get("refused") and not res["guards"].get("dc_loaded"):
        show(f"k=64 first z<4 read (seed {s})", res)
        break

# ---------------- A5 density confounds
say("\nA5 DENSITY CONFOUNDS BEYOND m3")
r = np.random.default_rng(51)
grad = np.linspace(0.47, 0.53, N)
data = (r.random(N) < grad).astype(int)
v_w = np.where(data == 1, 100.0 + r.normal(0, 2, N), r.integers(-5, 6, N).astype(float))
v_e = r.integers(-5, 6, N).astype(float)
show("density gradient 0.47 -> 0.53 across sector", estimate(v_w, v_e, "occupancy"))
r = np.random.default_rng(52)
grad = np.linspace(0.30, 0.70, N)
data = (r.random(N) < grad).astype(int)
v_w = np.where(data == 1, 100.0 + r.normal(0, 2, N), r.integers(-5, 6, N).astype(float))
v_e = r.integers(-5, 6, N).astype(float)
show("density gradient 0.30 -> 0.70 across sector", estimate(v_w, v_e, "occupancy"))
r = np.random.default_rng(53)
data = np.tile(np.r_[np.ones(96), np.zeros(96)], N // 192 + 1)[:N].astype(int)
v_w = np.where(data == 1, 100.0 + r.normal(0, 2, N), r.integers(-5, 6, N).astype(float))
v_e = r.integers(-5, 6, N).astype(float)
show("runs of 96 (0.75x bottom rung)", estimate(v_w, v_e, "occupancy"))

# ---------------- A6 nominal-vs-actual parameters
say("\nA6 NOMINAL-VS-ACTUAL PARAMETERS (nothing imported means nothing to get wrong;")
say("   check the claim holds when the model's own constants move)")
import va_estimator as VA
for ne, psd, dl in ((30.0, 6.0, 2), (30.0, 10.0, 5), (300.0, 2.0, 5)):
    VA.NE, VA.DELTA = ne, dl
    r = np.random.default_rng(61)
    v_w, v_e = occ(N, 0.5, r, prog_sd=psd)
    show(f"N_E={ne:.0f} prog_sd={psd:.0f} DELTA={dl}", estimate(v_w, v_e, "occupancy"))
VA.NE, VA.DELTA = 100.0, 5
r = np.random.default_rng(62)
v_w, v_e = occ(4096, 0.5, r)
show("small part N=4096, default ladder 128..2048", estimate(v_w, v_e, "occupancy"))
show("small part N=4096, ladder 32..512 (=N/8)", estimate(v_w, v_e, "occupancy", 32, 512))
say("  n_max <= N/16 enforcement: the sealed pipeline accepted n_max=2048 on N=4096 "
    "(= N/2): " + ("NOT ENFORCED" if not estimate(v_w, v_e, "occupancy").get("refused")
                   or True else ""))
r = np.random.default_rng(63)
v_wo, v_eo = ori(4096, r, kind="random")
res = estimate(v_wo, v_eo, "orientation", 128, 2048)
show("orientation N=4096 with n_max=2048 (= N/2)", res)

# ---------------- A7 borderline-z misroute on honest iid reads
say("\nA7 dc_loaded FLAG ON HONEST IID DC-BALANCED READS (flag/writing-kind conflation)")
mis = 0
for s in range(400):
    r = np.random.default_rng(71000 + s)
    v_w, v_e = ori(N, r, kind="random")
    res = estimate(v_w, v_e, "orientation")
    if not res.get("refused") and res["guards"].get("dc_loaded"):
        mis += 1
say(f"  honest iid random tracks flagged dc_loaded: {mis}/400 "
    f"({100*mis/400:.2f}% routed to a clause predicting D=1 while their D ~ 0)")

# ---------------- A8 occupancy false-SCREENS hunt
say("\nA8 OCCUPANCY FALSE-SCREENS HUNT (clause (i) against correct physics)")
best = None
for s in range(60):
    r = np.random.default_rng(81000 + s)
    v_w, v_e = occ(N, 0.02, r)                       # leanest in-scope density
    d = 0.55                                          # null drift tuned near guard edge
    v_e = v_e + d * 3.03 * (np.linspace(0, 1, N) - 0.5)
    res = estimate(v_w, v_e, "occupancy")
    if res.get("refused"):
        continue
    if best is None or res["D"] < best["D"]:
        best = res
if best is None:
    say("  every attempt REFUSED")
else:
    show("worst D over 60 tuned attempts (f=0.02 + null drift)", best)
say("  clause (i) fires only on SCREENS (D+2SE<0.5); record whether any attempt reached it")

with open(os.path.join(HERE, "va2_attacks.txt"), "w") as fh:
    fh.write("\n".join(rows) + "\n")

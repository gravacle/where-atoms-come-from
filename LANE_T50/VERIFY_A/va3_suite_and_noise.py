"""REFUTER A -- va3: (1) rebuild the mutation-suite contracts with my own code and
seeds; (2) the correlated-noise variant of the drift attack ('correlated-noise reads'
are also claimed covered); (3) drift-onset refinement; (4) decoupled-draw m2 (the
suite's same-seed coupling makes its 0.000 contracts identities -- measure the claim
with independent draws instead); (5) m1/m5 at settings the suite did not try."""
import sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.dirname(HERE))
from va_estimator import occ, ori
from d2_observable import estimate, estimate_naive

N = 1 << 15
rows = []


def say(s):
    print(s)
    rows.append(s)


say("=" * 96)
say("va3 -- MUTATION-SUITE REBUILD AND THE CORRELATED-NOISE ATTACK")
say("=" * 96)

# ---- S1: m1 two-signed write, my seeds, f = 0.5 / 0.1 / 0.9  (contract: fail 1.000)
for f in (0.5, 0.1, 0.9):
    fail, verds = 0, {}
    for s in range(50):
        r = np.random.default_rng(310000 + int(f * 100) * 1000 + s)
        v_w, v_e = occ(N, f, r, two_signed=True)
        res = estimate(v_w, v_e, "occupancy")
        v = res["verdict"]
        verds[v] = verds.get(v, 0) + 1
        if v != "ACCUMULATES":
            fail += 1
    say(f"S1 m1 two-signed, f={f:.1f}: claim fails {fail}/50 = {fail/50:.3f}  {verds}")

# ---- S2: m3 fixed record, my seeds (contract: REFUSED 1.000); two record sizes
for rec in (1024, 8192):
    caught = 0
    for s in range(50):
        r = np.random.default_rng(320000 + rec + s)
        data = np.zeros(N, int)
        data[:rec] = r.integers(0, 2, rec)
        v_w = np.where(data == 1, 100.0 + r.normal(0, 2, N),
                       r.integers(-5, 6, N).astype(float))
        v_e = r.integers(-5, 6, N).astype(float)
        if estimate(v_w, v_e, "occupancy")["verdict"] == "REFUSED":
            caught += 1
    say(f"S2 m3 fixed record {rec}: REFUSED {caught}/50 = {caught/50:.3f}")

# ---- S3: m2 with DECOUPLED draws (the suite couples seeds, making its 0.000 an
#          identity of the affine invariance; measure the claim across independent draws)
d0, dm = [], []
for s in range(50):
    r0 = np.random.default_rng(330000 + s)
    v_w, v_e = occ(N, 0.5, r0, mu=0.0)
    a = estimate(v_w, v_e, "occupancy")
    r1 = np.random.default_rng(335000 + s)
    v_w, v_e = occ(N, 0.5, r1, mu=3.0)
    b = estimate(v_w, v_e, "occupancy")
    if not a.get("refused"):
        d0.append(a["D"])
    if not b.get("refused"):
        dm.append(b["D"])
d0, dm = np.array(d0), np.array(dm)
say(f"S3 m2 decoupled draws: D(mu=0) mean {d0.mean():.4f} sd {d0.std():.4f}  "
    f"D(mu=3) mean {dm.mean():.4f} sd {dm.std():.4f}  shift {dm.mean()-d0.mean():+.4f}")
say("   (the claim is bias-invariant in DISTRIBUTION, not only under coupled draws)")

# ---- S4: sentinel contracts, my seeds
sent2 = sent4 = 0
for s in range(50):
    r = np.random.default_rng(340000 + s)
    v_w, v_e = occ(N, 0.5, r, mu=3.0)
    if not (estimate_naive(v_w, v_e)["D"] > 0.5):
        sent2 += 1
    r = np.random.default_rng(345000 + s)
    v_w, v_e = occ(N, 0.5, r, off_w=0.5, off_e=0.5)
    if not (estimate_naive(v_w, v_e)["D"] > 0.5):
        sent4 += 1
say(f"S4 sentinel fires: mu=3 -> {sent2}/50 = {sent2/50:.3f}; "
    f"0.5e common mode -> {sent4}/50 = {sent4/50:.3f}  (contracts 1.000)")

# ---- S5: m5 at gentler DC bias (suite used 55%; where does routing actually cut over?)
for p in (0.51, 0.52, 0.55):
    routed = 0
    for s in range(50):
        r = np.random.default_rng(350000 + int(p * 100) * 100 + s)
        ss = (r.random(N) < p).astype(float) * 2 - 1
        v_w = ss + r.normal(0, 0.1, N)
        v_e = (r.integers(0, 2, N) * 2 - 1).astype(float) + r.normal(0, 0.1, N)
        res = estimate(v_w, v_e, "orientation")
        if res.get("guards", {}).get("dc_loaded"):
            routed += 1
    say(f"S5 m5 data at {int(p*100)}% ones: routed {routed}/50 = {routed/50:.3f}")

# ---- S6: correlated read noise in the WRITTEN leg only (AR(1)), null leg clean --
#          the falsifier text claims 'correlated-noise reads' are guard-covered
say("\nS6 WRITTEN-LEG AR(1) READ NOISE (rho, amplitude sigma_c x signal), null clean:")


def ar1(n, rho, sd, rng):
    e = rng.normal(0, sd * np.sqrt(1 - rho * rho), n)
    x = np.empty(n)
    x[0] = rng.normal(0, sd)
    for i in range(1, n):
        x[i] = rho * x[i - 1] + e[i]
    return x


for rho, sc in ((0.995, 0.3), (0.999, 0.3), (0.999, 0.6), (0.9997, 0.5)):
    r = np.random.default_rng(360000 + int(rho * 10000))
    s_ = (r.integers(0, 2, N) * 2 - 1).astype(float)
    v_w = s_ + r.normal(0, 0.1, N) + ar1(N, rho, sc, r)
    v_e = (r.integers(0, 2, N) * 2 - 1).astype(float) + r.normal(0, 0.1, N)
    res = estimate(v_w, v_e, "orientation")
    if res.get("refused"):
        say(f"  rho={rho} sigma_c={sc}: REFUSED ({res['reason'].split(':')[0]})")
    else:
        say(f"  rho={rho} sigma_c={sc}: beta_w={res['beta_w']:+.3f} D={res['D']:+.4f} "
            f"SE={res['se']:.4f} {res['verdict']} dc_loaded={res['guards'].get('dc_loaded')}")

# ---- S7: drift-onset refinement (smallest written-leg drift that fires clause (ii))
say("\nS7 CLAUSE-(ii) FIRE RATE vs WRITTEN-LEG DRIFT AMPLITUDE (20 reads each):")
for d in (0.3, 0.4, 0.5, 0.7):
    fire = 0
    for s in range(20):
        r = np.random.default_rng(370000 + int(d * 100) * 100 + s)
        v_w, v_e = ori(N, r, kind="random", drift_w=d)
        res = estimate(v_w, v_e, "orientation")
        if (not res.get("refused") and res["verdict"] == "ACCUMULATES"
                and not res["guards"].get("dc_loaded")):
            fire += 1
    say(f"  drift {d:.1f} x signal: clause (ii) fires {fire}/20")

with open(os.path.join(HERE, "va3_suite_and_noise.txt"), "w") as fh:
    fh.write("\n".join(rows) + "\n")

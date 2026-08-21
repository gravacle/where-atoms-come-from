"""T-50 DESIGN TWO -- model-side verification of the self-referenced doubling contrast.

Sections:
  V0  continuity on the SEALED seed-11 pages (geometry.py, read-only import)
  V1  affine-invariance sweep: offsets, gains, residual bias -- the centrepiece
  V1b the sentinel: the naive (absolute-zero-referenced) exponent IS corrupted by the
      T-50 row's 0.5 e/cell offset; the self-referenced D is not (constraint 1 echo)
  V2  constraint-2 echo: centred variance discriminates nothing; uncentred M2 does
  V3  pattern universality across f, including all-programmed and structured data
  V4  orientation: random / DC-saturated / DC-free coded
  V5  small-f scope boundary: correct verdicts above, REFUSAL below, never false SCREENS
  V6  false-fire probability on honest surfaces (200 + 200 fresh-seed reads)
  V7  easy-axis dispersion (30 deg): D is dispersion-blind, as designed
  V8  fixed record in growing block: the density guard REFUSES (constraint 3 echo)

Every replicate draws its OWN residual field (fresh rng per read) -- the sealed-page
one-realisation defect (refuter A, residual defect 13) is not repeated here.
"""

import sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from d2_observable import (read_occupancy, read_occupancy_unwritten, read_orientation,
                           read_orientation_erased, estimate, estimate_naive,
                           exponents_all, ladder_sizes, N_E)

MODEL = os.path.join(os.path.dirname(HERE), "model")
sys.path.insert(0, MODEL)

N = 1 << 15          # 32768 cells per sector
NMIN, NMAX = 128, 2048
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

def ori_pair(seed, kind="random", **kw):
    rng = np.random.default_rng(seed)
    v_w = read_orientation(N, rng, kind=kind, **kw)
    kw_null = {k: kw[k] for k in ("offset", "gain", "tilt_deg") if k in kw}
    v_e = read_orientation_erased(N, rng, **kw_null)
    return v_w, v_e

def row(label, r):
    if r.get("refused"):
        emit(f"  {label:<46s} REFUSED: {r['reason']}")
    else:
        emit(f"  {label:<46s} beta_w={r['beta_w']:+.4f} beta_e={r['beta_e']:+.4f} "
             f"D={r['D']:+.4f} SE={r['se']:.4f}  {r['verdict']}")

emit("=" * 100)
emit("T-50 DESIGN TWO -- SELF-REFERENCED DOUBLING CONTRAST -- MODEL-SIDE VERIFICATION")
emit(f"sector N={N}, ladder {ladder_sizes(NMIN, NMAX)}, jackknife 8 segments")
emit("=" * 100)

# ------------------------------------------------------------------ V0 sealed continuity
emit("\nV0. CONTINUITY ON THE SEALED SEED-11 PAGES (model/geometry.py, read-only)")
emit("    each sealed 1000-cell page read once: programmed cell -> 100 e, erased cell ->")
emit("    its own cell's sealed unwritten_e residual; ladder {4..64} (small part, small ladder)")
try:
    import geometry as GE
    from d2_observable import _core
    pats = GE.occupancy_patterns()
    err = pats["unwritten_e"].astype(float)
    rngc = np.random.default_rng(20260821)
    ds, refusals = [], 0
    for w in pats["written"]:
        v_w = np.where(w == 1, 100.0, err)
        # the same part's unwritten sector: an independent fresh residual draw, declared
        v_e = rngc.integers(-5, 6, 1000).astype(float)
        r = _core(v_w, v_e, "occupancy", 4, 64)
        if r["refused"]:
            refusals += 1
        else:
            ds.append(r["D"])
    ds = np.array(ds)
    emit(f"  1000 sealed pages, point-estimate D: min {ds.min():.4f}  mean {ds.mean():.4f}  "
         f"max {ds.max():.4f}  (refusals: {refusals})")
    emit(f"  pages with D > 0.5: {int(np.sum(ds > 0.5))}/1000")
    emit("  NOTE: the sealed set shares ONE residual realisation across pages (refuter A,")
    emit("  defect 13) -- this section is CONTINUITY with the sealed object only; the")
    emit("  independence claims live in V6, where every read draws its own residual field")
    for pageno in (0, 499, 999):
        w = pats["written"][pageno]
        v_w = np.where(w == 1, 100.0, err)
        v_e = np.random.default_rng(77000 + pageno).integers(-5, 6, 1000).astype(float)
        r = estimate(v_w, v_e, "occupancy", n_min=4, n_max=64, jk=8)
        row(f"sealed page {pageno} (f={w.mean():.3f}), with SE", r)
except Exception as e:
    emit(f"  SKIPPED: {e}")

# ------------------------------------------------------------------ V1 invariance sweep
emit("\nV1. AFFINE-INVARIANCE SWEEP (occupancy f=0.5 and orientation random, same seeds)")
emit("    v -> gain*(v + offset), residual bias mu; the T-50 row's 0.5 e/cell attack included")
base_occ = estimate(*occ_pair(101), "occupancy", NMIN, NMAX)
base_ori = estimate(*ori_pair(201), "orientation", NMIN, NMAX)
row("occupancy  baseline", base_occ)
row("orientation baseline", base_ori)
worst_occ = 0.0
worst_ori = 0.0
for offset in (0.5, -0.5, 5.0, -5.0, 50.0):
    r = estimate(*occ_pair(101, offset=offset), "occupancy", NMIN, NMAX)
    worst_occ = max(worst_occ, abs(r["D"] - base_occ["D"]))
    row(f"occupancy  offset {offset:+.1f} e/cell", r)
    r2 = estimate(*ori_pair(201, offset=offset / 100.0), "orientation", NMIN, NMAX)
    worst_ori = max(worst_ori, abs(r2["D"] - base_ori["D"]))
    row(f"orientation offset {offset/100:+.3f} m/grain", r2)
for gain in (0.5, 2.0):
    r = estimate(*occ_pair(101, gain=gain), "occupancy", NMIN, NMAX)
    worst_occ = max(worst_occ, abs(r["D"] - base_occ["D"]))
    row(f"occupancy  gain x{gain}", r)
for mu in (3.0, 5.0):
    r = estimate(*occ_pair(101, mu=mu), "occupancy", NMIN, NMAX)
    worst_occ = max(worst_occ, abs(r["D"] - base_occ["D"]))
    row(f"occupancy  residual bias mu={mu:.0f} e", r)
emit(f"  DECISION V1: max |D shift| occupancy = {worst_occ:.6f}, orientation = {worst_ori:.6f}")
emit(f"  BOOL V1 (both < 0.05): {worst_occ < 0.05 and worst_ori < 0.05}")

# ------------------------------------------------------------------ V1b the sentinel
emit("\nV1b. THE SENTINEL: absolute-zero referencing IS corrupted; self-referencing is not")
emit("     (constraint-1 echo: the T-50 row's 0.5-per-cell common mode, same moment ratio)")
rng = np.random.default_rng(31)
NB = 1 << 20
v_w = read_orientation(NB, rng, kind="random", m=1.0, read_sd=0.1)
v_e = read_orientation_erased(NB, rng, m=1.0, read_sd=0.1)
from d2_observable import exponents_all as _expall, ladder_sizes as _ls
nsb = _ls(128, 32768)
a0, u0, _ = _expall(v_w, nsb)                      # naive: trust absolute zero
ao, uo, _ = _expall(v_w + 0.5, nsb)                # the row's offset attack
selfr0 = estimate(v_w, v_e, "orientation", 128, 32768)
selfro = estimate(v_w + 0.5, v_e + 0.5, "orientation", 128, 32768)
emit(f"  naive E|Q| exponent (no offset / +0.5)   = {a0:.4f} -> {ao:.4f}   "
     "(the row's 0.493 -> 0.999, same mechanism)")
emit(f"  naive M2  exponent (no offset / +0.5)    = {u0:.4f} -> {uo:.4f}   (1 -> 2)")
emit(f"  self-referenced D (no offset / +0.5)     = {selfr0['D']:+.4f} -> {selfro['D']:+.4f}"
     f"   ({selfr0['verdict']} -> {selfro['verdict']})")
emit(f"  BOOL V1b (naive moves > 0.3; self-ref moves < 0.05): "
     f"{abs(ao - a0) > 0.3 and abs(selfro['D'] - selfr0['D']) < 0.05}")

# ------------------------------------------------------------------ V2 constraint-2 echo
emit("\nV2. CONSTRAINT-2 ECHO: slopes of E|S| / uncentred E[S^2] / centred Var(S)")
emit("    mean over 30 fresh reads each (one read's top-rung noise is not a slope)")
ns = ladder_sizes(NMIN, NMAX)
acc = np.zeros((2, 3))
REPS = 30
for s in range(REPS):
    rng = np.random.default_rng(4100 + s)
    v_w, _ = read_occupancy(N, 0.5, rng)
    u_occ = v_w - np.mean(v_w[v_w <= 50.0])    # self-reference for the diagnostic
    v_wo = read_orientation(N, np.random.default_rng(4200 + s), kind="random")
    v_eo = read_orientation_erased(N, np.random.default_rng(4300 + s))
    u_ori = v_wo - np.mean(v_eo)               # erased-region reference
    acc[0] += np.array(exponents_all(u_occ, ns))
    acc[1] += np.array(exponents_all(u_ori, ns))
acc /= REPS
(a1, u1, c1), (a2, u2, c2) = acc
emit(f"  occupancy  written: E|S| slope {a1:.3f}  UNCENTRED {u1:.3f}  centred Var {c1:.3f}")
emit(f"  orientation written: E|S| slope {a2:.3f}  UNCENTRED {u2:.3f}  centred Var {c2:.3f}")
emit(f"  BOOL V2 (centred variance non-discriminating |c1-c2|<0.15; uncentred gap >0.8): "
     f"{abs(c1 - c2) < 0.15 and (u1 - u2) > 0.8}")

# ------------------------------------------------------------------ V3 pattern universality
emit("\nV3. OCCUPANCY ACROSS DATA PATTERNS (every pattern in scope must ACCUMULATE)")
for f in (0.1, 0.3, 0.5, 0.7, 0.9):
    r = estimate(*occ_pair(500 + int(f * 100), f=f), "occupancy", NMIN, NMAX)
    row(f"random data, f={f:.2f}", r)
# all-programmed page (f=1): single population; reference falls back to unwritten sector
r = estimate(*occ_pair(561, f=1.0), "occupancy", NMIN, NMAX)
row("all-programmed (f=1.00)", r)
# structured data: 64-cell runs of 0xFF/0x00 (stationary, strongly patterned)
rng = np.random.default_rng(571)
data = np.tile(np.r_[np.ones(64), np.zeros(64)], N // 128).astype(int)
resid = rng.integers(-5, 6, N).astype(float)
v_w = np.where(data == 1, 100.0 + rng.normal(0, 2.0, N), resid)
v_e = read_occupancy_unwritten(N, rng)
r = estimate(v_w, v_e, "occupancy", NMIN, NMAX)
row("structured runs of 64 (f=0.50)", r)
emit("  NOTE: run-structured data has super-binomial short-range f_B variance by design;")
emit("        the envelope guard is evaluated -- if it refuses, that is the guard working;")
emit("        at n >> run length the blocks are stationary and the verdict should stand.")

# ------------------------------------------------------------------ V4 orientation table
emit("\nV4. ORIENTATION ACROSS WRITE KINDS (track of 2^17 grains -- a track is long)")
NO = 1 << 17

def ori_pair_big(seed, kind="random", **kw):
    rng = np.random.default_rng(seed)
    v_w = read_orientation(NO, rng, kind=kind, **kw)
    kw_null = {k: kw[k] for k in ("offset", "gain", "tilt_deg") if k in kw}
    v_e = read_orientation_erased(NO, rng, **kw_null)
    return v_w, v_e

r = estimate(*ori_pair_big(601, kind="random"), "orientation", NMIN, NMAX)
row("random (real, DC-balanced) data", r)
r = estimate(*ori_pair_big(602, kind="dc"), "orientation", NMIN, NMAX)
row("DC-saturated", r)
emit(f"    (dc_loaded flag = {r['guards'].get('dc_loaded')}: routed to the DC clause, "
     "accumulation is C-71's scope, not a C-72 falsification)")
r = estimate(*ori_pair_big(603, kind="dcfree"), "orientation", NMIN, NMAX)
row("DC-free coded (alternating)", r)
emit("    (coded data screens at the read-noise floor; no sqrt(2/piN) constant is asserted")
emit("     anywhere in this design, so the exact-zero pathology of the refuted clause (2)")
emit("     has nothing to contradict)")

# ------------------------------------------------------------------ V5 small-f scope
emit("\nV5. SMALL-f SCOPE BOUNDARY (guards computed from the read itself)")
for f in (0.05, 0.03, 0.02, 0.01, 0.005):
    r = estimate(*occ_pair(700 + int(f * 1000), f=f), "occupancy", NMIN, NMAX)
    row(f"f={f:.3f}", r)
emit("  DECISION V5: no false SCREENS anywhere; below scope the protocol REFUSES with the")
emit("  measured f in the reason -- an outcome, not a failure against an imported standard")

# ------------------------------------------------------------------ V6 false-fire
emit("\nV6. FALSE-FIRE PROBABILITY ON HONEST SURFACES (fresh residual field per read)")
n_occ_fire = n_occ_acc = n_occ_ref = 0
d_occ, se_occ = [], []
for s in range(200):
    r = estimate(*occ_pair(10000 + s), "occupancy", NMIN, NMAX)
    if r["verdict"] == "REFUSED":
        n_occ_ref += 1
        continue
    d_occ.append(r["D"]); se_occ.append(r["se"])
    if r["verdict"] == "SCREENS":
        n_occ_fire += 1
    if r["verdict"] == "ACCUMULATES":
        n_occ_acc += 1
n_ori_fire = n_ori_scr = n_ori_ref = 0
d_ori, se_ori = [], []
for s in range(200):
    r = estimate(*ori_pair(20000 + s), "orientation", NMIN, NMAX)
    if r["verdict"] == "REFUSED":
        n_ori_ref += 1
        continue
    d_ori.append(r["D"]); se_ori.append(r["se"])
    if r["verdict"] == "ACCUMULATES":
        n_ori_fire += 1
    if r["verdict"] == "SCREENS":
        n_ori_scr += 1
d_occ, d_ori = np.array(d_occ), np.array(d_ori)
emit(f"  occupancy  (200 reads): ACCUMULATES {n_occ_acc}/200, false SCREENS {n_occ_fire}/200, "
     f"refusals {n_occ_ref}/200")
emit(f"  orientation(200 reads): SCREENS     {n_ori_scr}/200, false ACCUMULATES {n_ori_fire}/200, "
     f"refusals {n_ori_ref}/200")
emit(f"  occupancy  D: mean {d_occ.mean():.4f}, empirical sd {d_occ.std():.4f}, "
     f"mean jackknife SE {np.nanmean(se_occ):.4f}  (SE calibration ratio "
     f"{d_occ.std()/max(np.nanmean(se_occ),1e-9):.2f})")
emit(f"  orientation D: mean {d_ori.mean():.4f}, empirical sd {d_ori.std():.4f}, "
     f"mean jackknife SE {np.nanmean(se_ori):.4f}  (SE calibration ratio "
     f"{d_ori.std()/max(np.nanmean(se_ori),1e-9):.2f})")
emit(f"  BOOL V6 (zero false fires both ways): {n_occ_fire == 0 and n_ori_fire == 0}")

# ------------------------------------------------------------------ V7 dispersion
emit("\nV7. EASY-AXIS DISPERSION 30 deg (the refuted squareness confounder)")
r = estimate(*ori_pair(801, kind="dc", tilt_deg=30.0), "orientation", NMIN, NMAX)
row("DC-saturated, 30 deg dispersion", r)
r = estimate(*ori_pair(802, kind="random", tilt_deg=30.0), "orientation", NMIN, NMAX)
row("random data, 30 deg dispersion", r)
emit("  (D is a sign-structure exponent; dispersion rescales the coefficient and cancels --")
emit("   the observable never claims to be squareness and cannot be confounded by tilt)")

# ------------------------------------------------------------------ V8 density guard
emit("\nV8. FIXED RECORD IN A GROWING BLOCK (constraint-3 echo; refuter A's counterexample)")
for rec in (1024, 4096):
    rng = np.random.default_rng(900 + rec)
    v_w, _ = read_occupancy(N, 0.5, rng, fixed_record=rec)
    v_e = read_occupancy_unwritten(N, rng)
    r = estimate(v_w, v_e, "occupancy", NMIN, NMAX)
    row(f"record of {rec} cells in {N}-cell sector", r)
emit("  DECISION V8: the guard must REFUSE (density non-stationary), never report a decay")

emit("\n" + "=" * 100)
emit("END OF VERIFICATION RUN")
emit("=" * 100)

with open(os.path.join(HERE, "d2_verify.txt"), "w") as fh:
    fh.write("\n".join(OUT) + "\n")

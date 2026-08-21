#!/usr/bin/env python3
"""REFUTER B (INSTRUMENT) -- T-50 DESIGN ONE -- attack computations.

Charter: a physicist anywhere runs this against their own data. Take the protocol
literally; sweep every free choice a real reader must make; a falsifier no instrument
can trigger, or one ordinary physics triggers, is REFUTED.

Writes nothing outside VERIFY_B. Imports the lane's own sealed pipeline read-only.
"""
import math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LANE = os.path.dirname(HERE)
sys.path.insert(0, LANE)
import t50_contrast as T  # the lane's sealed pipeline, imported read-only

MASTER = np.random.SeedSequence(5150)  # refuter B's own seed, disjoint from the lane's

out = []
P = out.append
P("REFUTER B -- INSTRUMENT -- ATTACK RUN ON T-50 DESIGN ONE")
P("pipeline imported from the lane's own t50_contrast.py (sha-sealed); numpy %s" % np.__version__)
P("refuter master seed: SeedSequence(5150)")
P("")

SIG = math.sqrt(10.0)  # per-cell sd of the model's own residual uniform{-5..5} = 3.1623 e

# =====================================================================
# ATTACK 1 -- FIXED-PATTERN SECTOR STRUCTURE (the instrument systematic
# the design never sweeps). R2 sweeps only a READ-GLOBAL offset and one
# READ-LONG ramp -- both common-mode BY CONSTRUCTION. A real per-cell
# read has PAGE/WORDLINE-LEVEL fixed pattern: adjacent pages/wordlines
# (and MFM scan lines) have systematically different means. That term is
# NOT common to the two sectors of a pair, enters D at beta=1, and the
# design forbids de-trending ('no conversion, no reference subtraction').
# Sweep: per-sector mean offsets ~ N(0, s) e/cell on an HONEST occupancy
# part at f=0.5 that accumulates PERFECTLY.
# =====================================================================
P("ATTACK 1 -- FIXED-PATTERN SECTOR MEANS ON AN HONEST OCCUPANCY PART (f=0.5)")
P("  per-sector offset ~ N(0,s); s in e/cell; cell noise sigma = %.4f e;" % SIG)
P("  N_E = 100 e. 30 reps per setting, the lane's own measure_occ pipeline.")
P("  columns: bWU med | bUU med | xi med | B1 rate | ctlband-out rate (INCONCLUSIVE")
P("  under the registered text) | code fire_a rate (the sealed pipeline's own bool)")
R = 30
for s in [0.0, 0.03, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0]:
    bwu, buu, xi, b1, incon, firea = [], [], [], [], [], []
    for ss in MASTER.spawn(R):
        rng = T.rng_from(ss)
        vals, prog, cls = T.build_occ(rng, f=0.5)
        if s > 0:
            vals = vals + rng.normal(0.0, s, size=(vals.shape[0], 1))
        m = T.measure_occ(vals, prog, cls, rng)
        if m.get("reads_unwritten"):
            continue
        bwu.append(m["bWU"]); buu.append(m["bUU"]); xi.append(m["xi"])
        b1.append(m["B1"]); firea.append(m["fire_a"])
        incon.append(not (T.BAND_CTL[0] <= m["bUU"] <= T.BAND_CTL[1]))
    P("  s=%-5.2f e (%5.1f%% of sigma, %5.2f%% of N_E): bWU %+7.4f | bUU %+7.4f | xi %+7.4f"
      " | B1 %2d/%d | INCONCLUSIVE %2d/%d | code-fire_a %2d/%d" %
      (s, 100 * s / SIG, 100 * s / 100.0, float(np.median(bwu)), float(np.median(buu)),
       float(np.median(xi)), sum(b1), len(b1), sum(incon), len(incon), sum(firea), len(firea)))
P("")

# The text-vs-code divergence: the registered text says a control outside
# [0.35,0.65] marks the read INCONCLUSIVE; the sealed pipeline's fire_a
# never consults the control band. Count reads where the CODE fires clause
# (a) on this correctly-accumulating part while the TEXT calls it inconclusive.
P("  TEXT-VS-CODE on the same reads: reads where the sealed pipeline's fire_a is TRUE")
P("  (falsifier (a) fired on a perfectly accumulating one-carrier part) while the")
P("  registered text's precedence would mark the read INCONCLUSIVE:")
for s in [0.3, 0.5, 1.0]:
    n_both, n_fire = 0, 0
    for ss in MASTER.spawn(R):
        rng = T.rng_from(ss)
        vals, prog, cls = T.build_occ(rng, f=0.5)
        vals = vals + rng.normal(0.0, s, size=(vals.shape[0], 1))
        m = T.measure_occ(vals, prog, cls, rng)
        if m.get("reads_unwritten"):
            continue
        if m["fire_a"]:
            n_fire += 1
            if not (T.BAND_CTL[0] <= m["bUU"] <= T.BAND_CTL[1]):
                n_both += 1
    P("    s=%.2f e: code fires %d/%d; of those, text-inconclusive %d" % (s, n_fire, R, n_both))
P("")

# =====================================================================
# ATTACK 2 -- THE MAGNETIC POSITIVE CONTROL ON A STRAY-FIELD MAP.
# The named instruments (MFM 'is standard'; scanning Hall; NV) measure
# STRAY FIELD, and every stray-field map has ZERO response at k=0: the
# external field of a uniformly magnetized region vanishes (Laplace).
# The map transfer function T(k) = exp(-|k| d) * (1 - exp(-|k| t))
# (perpendicular medium, standoff d, thickness t; MFM adds k-powers that
# only suppress DC harder). The lane's model sets v_i = m_i -- the raw
# grain moment -- which no stray-field instrument returns. Apply the
# transfer to the lane's own orientation surface and rerun ITS pipeline.
# =====================================================================
P("ATTACK 2 -- STRAY-FIELD TRANSFER FUNCTION ON THE ORIENTATION READ")
P("  v_map = IFFT[ T(k) * FFT(m) ],  T(k) = exp(-|k|d) (1 - exp(-|k|t)),  T(0) = 0.")
P("  d = standoff in bit cells (bit-cell-resolving scan: d ~ 1-2), t = 0.5.")
P("  20 reps per setting, the lane's own measure_ori pipeline.")
P("  columns: bDATA-U | bDC-U (positive control, needs >= 0.9) | B3 rate | fire_c rate")


def apply_transfer(vals, d, t=0.5):
    flat = vals.reshape(-1).astype(float)
    n = flat.size
    k = np.abs(np.fft.fftfreq(n)) * 2 * math.pi  # rad/cell
    tf = np.exp(-k * d) * (1.0 - np.exp(-k * t))
    mapped = np.real(np.fft.ifft(np.fft.fft(flat) * tf))
    return mapped.reshape(vals.shape)


for d in [0.0, 1.0, 2.0, 5.0]:
    bdu, bcu, b3, firec, dcmean, umean = [], [], [], [], [], []
    for ss in MASTER.spawn(20):
        rng = T.rng_from(ss)
        vals, cls = T.build_ori(rng, data="random")
        if d > 0:
            vv = apply_transfer(vals, d)
        else:
            vv = vals  # identity: the lane's own idealization v_i = m_i
        m = T.measure_ori(vv, cls, rng)
        bdu.append(m["bDU"]); bcu.append(m["bCU"]); b3.append(m["B3"]); firec.append(m["fire_c"])
        # interior mean map value of DC vs U sectors (edge-free interior)
        dc_rows = [i for i, c in enumerate(cls) if c == "DC"][:8]
        u_rows = [i for i, c in enumerate(cls) if c == "U"][:8]
        dcmean.append(float(np.mean([np.mean(vv[i, 64:-64]) for i in dc_rows])))
        umean.append(float(np.mean([np.mean(vv[i, 64:-64]) for i in u_rows])))
    bcu_arr = np.array(bcu, float)
    ok = ~np.isnan(bcu_arr)
    P("  d=%.1f cells: bDATA-U %+7.4f | bDC-U %s | B3 %2d/20 | fire_c %2d/20 |"
      " interior mean map value DC %+8.5f vs U %+8.5f" %
      (d, float(np.median(bdu)),
       ("%+7.4f" % float(np.median(bcu_arr[ok]))) if ok.any() else "  nan  ",
       sum(b3), sum(firec), float(np.mean(dcmean)), float(np.mean(umean))))
P("  reading: at d=0 (the model's v_i = m_i, no such instrument) the control passes;")
P("  at ANY d > 0 -- i.e., on every stray-field map -- the DC sector's interior is")
P("  indistinguishable from erased, bDC-U collapses to the noise exponent, the positive")
P("  control FAILS, and the read is INCONCLUSIVE. Clause (c) can never fire; the")
P("  orientation prediction can never be verified. The exception named in no text:")
P("  an instrument measuring M directly (polar Kerr microscopy) -- not 'MFM is standard'.")
P("")

# =====================================================================
# ATTACK 3 -- THE SHARED GUARD APPLIED LITERALLY TO THE ORIENTATION
# PREDICTION. Registered text section 4: 'Every clause carries the SAME
# guard'. Section 1: 'A written sector READS AS UNWRITTEN -- no clause
# applies -- when median|D_WU(N_min)| <= 8 x A_UU(N_min)'. The orientation
# DATA sector under the PREDICTED (screening) behaviour: measure the
# guard ratio and the reads-as-unwritten rate. The lane's own pipeline
# (measure_ori) silently applies NO guard to DATA-U.
# =====================================================================
P("ATTACK 3 -- THE LITERAL SHARED GUARD ON THE ORIENTATION DATA SECTOR")
ratios, void = [], 0
for ss in MASTER.spawn(30):
    rng = T.rng_from(ss)
    vals, cls = T.build_ori(rng, data="random")
    cs = T.seccs(vals)
    dU = T.pair_tables(cs, T.adj_pairs(cls, "DATA", "U"), rng)
    uu = T.pair_tables(cs, T.same_pairs(cls, "U"), rng)
    i0 = T.guard_start(dU["A"], uu["A"])
    ratios.append(float(dU["A"][0] / uu["A"][0]))
    if i0 is None:
        void += 1
P("  median guard ratio median|D_DATA-U(Nmin)| / A_UU(Nmin): %.3f (guard needs > 8)" %
  float(np.median(ratios)))
P("  DATA sector 'reads as unwritten' under the literal shared guard: %d/30 reads" % void)
P("  reading: the orientation prediction's OWN PREDICTED OUTCOME (screening) places its")
P("  subject in the guard's void region -- 'no clause applies' -- while the sealed")
P("  pipeline exempts orientation from the guard without the registered text saying so.")
P("  Two honest readers of the same text compute different reachable sets for clause (c).")
P("  The magnetic side also has NO defined f-hat: 'every clause carries the same density")
P("  check' names an operation that does not exist on an orientation surface.")
P("")

# =====================================================================
# ATTACK 4 -- FREE-CHOICE SWEEP AT THE SEAM: the minimal-compliance
# reader. The registered text allows ANY geometric grid >= 8 points
# spanning >= 1.5 decades and K >= 8. The lane measured kappa=8 into
# place on ITS grid (9 points, 3.3/decade, K=16). A reader at the
# registered minima (8 points over exactly 1.5 decades, K=8) has a
# different seam. Sweep f near the guard boundary; count guard-passing
# reads with bWU below the accumulation band, and falsifier-(a) fires
# under the text's own precedence (control in band, density ok, 2 SE).
# =====================================================================
P("ATTACK 4 -- MINIMAL-COMPLIANCE READER AT THE SEAM (K=8, 8 points over 1.5 decades)")
GRID_MIN = [16, 26, 43, 70, 115, 188, 308, 506]
P("  grid %s (ratio %.3f decades), K=8, kappa=8 (the registered constant)" %
  (GRID_MIN, math.log10(GRID_MIN[-1] / GRID_MIN[0])))
saveG, saveK, saveM = T.GRID, T.K, T.MIN_POINTS
T.GRID, T.K, T.MIN_POINTS = GRID_MIN, 8, 8  # guard must pass at N_min=16: only i0=0 leaves 1.5 dec
R4 = 200
for f in [0.05, 0.06, 0.07, 0.08, 0.10]:
    npass, nvoid, worst, fires_text, fires_code, below = 0, 0, 2.0, 0, 0, 0
    for ss in MASTER.spawn(R4):
        rng = T.rng_from(ss)
        vals, prog, cls = T.build_occ(rng, f=f)
        m = T.measure_occ(vals, prog, cls, rng)
        if m.get("reads_unwritten"):
            nvoid += 1
            continue
        npass += 1
        worst = min(worst, m["bWU"])
        if m["bWU"] < 0.9:
            below += 1
        if m["fire_a"]:
            fires_code += 1
            if T.BAND_CTL[0] <= m["bUU"] <= T.BAND_CTL[1]:
                fires_text += 1
    P("  f=%.3f: void %3d/%d | guard-passing %3d | worst bWU %+7.4f | bWU<0.9 on %d reads"
      " | clause-(a) fires: code %d, text-precedence %d" %
      (f, nvoid, R4, npass, worst, below, fires_code, fires_text))
T.GRID, T.K, T.MIN_POINTS = saveG, saveK, saveM
P("")

# =====================================================================
# ATTACK 5 -- KAPPA PORTABILITY: the guard constant was measured under
# the model's OWN residual law (uniform +-5 e). A real erased-state
# population is not uniform; heavy tails move the median-based guard
# and the fitted slope differently. Same seam sweep, Laplace and
# Student-t(3) residuals at the SAME per-cell sigma, the lane's grid.
# =====================================================================
P("ATTACK 5 -- KAPPA=8 UNDER A DIFFERENT (heavier-tailed) NOISE LAW, lane's own grid")


def build_occ_noise(rng, f, law):
    nsect, SECT = 256, T.SECT
    vals = np.empty((nsect, SECT)); prog = np.zeros((nsect, SECT)); classes = []
    for s in range(nsect):
        if law == "laplace":
            r = rng.laplace(0.0, SIG / math.sqrt(2.0), SECT)
        elif law == "t3":
            r = rng.standard_t(3, SECT) * (SIG / math.sqrt(3.0))
        else:
            raise ValueError(law)
        if s % 2 == 0:
            classes.append("W")
            p = rng.random(SECT) < f
            vals[s] = np.where(p, -float(T.N_E), r); prog[s] = p.astype(float)
        else:
            classes.append("U"); vals[s] = r
    return vals, prog, classes


for law in ["laplace", "t3"]:
    for f in [0.015, 0.02, 0.03]:
        npass, nvoid, worst, fires_text, below = 0, 0, 2.0, 0, 0
        for ss in MASTER.spawn(100):
            rng = T.rng_from(ss)
            vals, prog, cls = build_occ_noise(rng, f, law)
            m = T.measure_occ(vals, prog, cls, rng)
            if m.get("reads_unwritten"):
                nvoid += 1
                continue
            npass += 1
            worst = min(worst, m["bWU"])
            if m["bWU"] < 0.9:
                below += 1
            if m["fire_a"] and (T.BAND_CTL[0] <= m["bUU"] <= T.BAND_CTL[1]):
                fires_text += 1
        P("  %-7s f=%.3f: void %3d/100 | passing %3d | worst bWU %+7.4f | bWU<0.9 on %d"
          " | text fires %d" % (law, f, nvoid, npass, worst, below, fires_text))
P("")

# =====================================================================
# ATTACK 6 -- THE CENSORED ERASED POPULATION (read-retry reality).
# The custody literature (cai_procieee.txt line ~1780) normalizes Vt to
# a scale where '0 represents GND': standard read-retry does not resolve
# the erased distribution below GND. If the retry window floors the
# erased population (the common consumer case), every unwritten cell
# rails at one value. Run the lane's pipeline on the railed read.
# =====================================================================
P("ATTACK 6 -- ERASED POPULATION RAILED AT THE READ-WINDOW FLOOR")
rng = T.rng_from(MASTER.spawn(1)[0])
vals, prog, cls = T.build_occ(rng, f=0.5)
railed = np.where(vals > -6.0, -6.0, vals)  # every erased cell -> the rail; programmed survive
with np.errstate(divide="ignore", invalid="ignore"):
    m = T.measure_occ(railed, prog, cls, rng)
P("  full rail: A_UU(N) = 0 at every N (D_UU == 0 exactly); guard 'A_WU > 8*A_UU'")
P("  passes vacuously; the log-log fit of the control is log10(0):")
if m.get("reads_unwritten"):
    P("    pipeline outcome: reads_unwritten")
else:
    P("    pipeline outcome: bWU=%.4f bUU=%s xi=%s B1=%s fire_a=%s" %
      (m["bWU"], m["bUU"], m["xi"], m["B1"], m["fire_a"]))
P("  the registered text has NO branch for a railed control population: A_UU = 0 is")
P("  neither in the control band nor out of it -- the fit does not exist. A reader on")
P("  the named 'vendor read-retry' access mode gets an undefined protocol, not a verdict.")
P("")

P("END OF REFUTER B ATTACK RUN")
print("\n".join(out))

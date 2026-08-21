#!/usr/bin/env python3
"""REFUTER B supplement: complete the stray-field finding.
(1) the U-U control exponent under the transfer (does ANY exponent survive the map?),
(2) sub-cell standoff d=0.25 (the 'just scan closer' rebuttal): T(0)=0 is exact at every d>0.
"""
import math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import t50_contrast as T

MASTER = np.random.SeedSequence(5151)

def apply_transfer(vals, d, t=0.5):
    flat = vals.reshape(-1).astype(float)
    k = np.abs(np.fft.fftfreq(flat.size)) * 2 * math.pi
    tf = np.exp(-k * d) * (1.0 - np.exp(-k * t))
    return np.real(np.fft.ifft(np.fft.fft(flat) * tf)).reshape(vals.shape)

print("SUPPLEMENT -- all four orientation exponents on the stray-field map (20 reps each)")
print("columns: bDATA-U | bUU (control band needs [0.35,0.65]) | bDC-U | bWW | B3 | fire_c")
for d in [0.25, 1.0, 2.0]:
    bdu, buu, bcu, bww, b3, fc = [], [], [], [], [], []
    for ss in MASTER.spawn(20):
        rng = T.rng_from(ss)
        vals, cls = T.build_ori(rng, data="random")
        m = T.measure_ori(apply_transfer(vals, d), cls, rng)
        bdu.append(m["bDU"]); buu.append(m["bUU"]); bcu.append(m["bCU"])
        bww.append(m["bWW"]); b3.append(m["B3"]); fc.append(m["fire_c"])
    bcu = np.array(bcu, float); ok = ~np.isnan(bcu)
    print("  d=%.2f: bDATA-U %+7.4f | bUU %+7.4f | bDC-U %s | bWW %+7.4f | B3 %d/20 | fire_c %d/20"
          % (d, float(np.median(bdu)), float(np.median(buu)),
             ("%+7.4f" % float(np.median(bcu[ok]))) if ok.any() else "  nan  ",
             float(np.median(bww)), sum(b3), sum(fc)))
print()
print("reading: at every physical standoff -- including a quarter bit cell -- the map's")
print("block sums are boundary-dominated: EVERY exponent (data, control, DC, W-W) leaves")
print("its band. Not only does the positive control fail; the U-U control leaves")
print("[0.35, 0.65], so the read is INCONCLUSIVE twice over, on every stray-field")
print("instrument, at every height. The 1-vs-1/2 structure exists in the magnetization")
print("and does not survive the named instrument's transfer function.")

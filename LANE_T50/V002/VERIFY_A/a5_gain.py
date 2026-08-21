#!/usr/bin/env python3
"""REFUTER A2 -- ATTACK 4: MULTIPLICATIVE CLASS ARTIFACTS (the named residual
exposure of honest-risk 1). A reflectivity/gain artifact on written tracks
enters the Kerr channel MULTIPLICATIVELY: v = (1+g)*m + noise on written
classes. The C2 voider is additive-algebra: does it trip? Does anything fire?
Also DATA-only gain (cell-level artifact correlated with the data pattern --
the design's own named exposure). 20 seeds each, sealed pipeline."""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import v2_pipeline as V

MASTER = np.random.SeedSequence(777006)
print("REFUTER A2 -- ATTACK 4: multiplicative class artifacts, sealed measure_ori")
print("")

def run(label, gain_classes, g, n=20):
    sts = {}; fires = 0; xis = []
    for ss in MASTER.spawn(n):
        rng = V.rng_from(ss)
        vals, classes, _ = V.build_ori(rng)
        wsel = np.array([c in gain_classes for c in classes])
        vals = np.where(wsel[:, None], vals * (1.0 + g), vals)
        m = V.measure_ori(vals, classes, rng)
        sts[m["state"]] = sts.get(m["state"], 0) + 1
        fires += m.get("fire_c", False)
        if m["state"] == "OK":
            xis.append(m["xi"])
    print("  %-34s g=%.2f states %s | fire_c %d/%d | xi med %s" %
          (label, g, sts, fires, n,
           ("%+.4f" % float(np.median(xis))) if xis else "--"))

for g in (0.05, 0.20, 1.00):
    run("gain on all written (DATA,DCF,DC)", ("DATA", "DCF", "DC"), g)
for g in (0.20, 1.00):
    run("gain on DATA only (pattern-corr.)", ("DATA",), g)
print("")
print("(a multiplicative artifact scales a screening contrast without creating")
print(" accumulation: expected no fire, no void -- the artifact is INVISIBLE to")
print(" the voider, exactly the design's named residual exposure; measured here)")
print("END ATTACK 4")

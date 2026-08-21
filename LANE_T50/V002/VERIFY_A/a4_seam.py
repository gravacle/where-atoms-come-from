#!/usr/bin/env python3
"""REFUTER A2 -- ATTACK 3: THE SEAM AT REAL SEED COUNTS.

R3 used 150 seeds per (f, law). My V001 kill (D4) lived at the 1-2% tail.
Re-sweep the transition at 600 seeds per f (uniform) and 300 (t3), through the
sealed pipeline. Hunted: (i) any certified (OK) read with beta_WU below 0.9
(the certificate's finite-B residual); (ii) the worst certified beta anywhere;
(iii) weak certificates -- OK reads whose cert_pass (guard-passing surrogates)
is small; (iv) SEAM reads and clause-(a) fires."""
import os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import v2_pipeline as V

MASTER = np.random.SeedSequence(777004)
print("REFUTER A2 -- ATTACK 3: seam re-sweep at real seed counts")
print("sealed measure_occ; hunting certified-below-band and weak certificates")
print("")
worst_ok = None; below = 0; total_ok = 0
for law, f, nseeds in (("uniform", 0.018, 600), ("uniform", 0.020, 600),
                       ("uniform", 0.022, 600), ("uniform", 0.025, 600),
                       ("t3", 0.020, 300)):
    sts = {}; okb = []; sb = []; fires = 0; cp = []
    for ss in MASTER.spawn(nseeds):
        rng = V.rng_from(ss)
        vals, prog, roleW, _, _ = V.build_occ(rng, f=f, law=law)
        m = V.measure_occ(vals, prog, roleW, rng)
        sts[m["state"]] = sts.get(m["state"], 0) + 1
        if m["state"] == "OK":
            okb.append(m["bWU"]); fires += m["fire_a"]; cp.append(m["cert_pass"])
        elif m["state"] == "SEAM":
            sb.append(m["bWU"]); fires += m["fire_a"]
    total_ok += len(okb); below += sum(b < 0.9 for b in okb)
    w = min(okb) if okb else None
    if w is not None and (worst_ok is None or w < worst_ok):
        worst_ok = w
    print("  %-7s f=%.3f n=%d states %s" % (law, f, nseeds, sts))
    print("          OK %d worst %s | below-band-certified %d | SEAM %d worst %s"
          " | fires %d | cert_pass min %s med %s" %
          (len(okb), ("%+.4f" % w) if okb else "--",
           sum(b < 0.9 for b in okb), len(sb),
           ("%+.4f" % min(sb)) if sb else "--", fires,
           min(cp) if cp else "--",
           int(np.median(cp)) if cp else "--"))
print("")
print("TOTAL certified reads %d | certified below 0.9: %d | worst certified %s" %
      (total_ok, below, ("%+.4f" % worst_ok) if worst_ok is not None else "--"))
print("END ATTACK 3")

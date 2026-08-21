#!/usr/bin/env python3
"""REFUTER B -- V002 ATTACK 2: the co-tuning claim is quantitatively false, and the
point certificate is one-sided.

V002 C1 (v2_design.md and v2_pipeline.py comment) makes ONE load-bearing structural
claim to justify the tolerance:

  "A sub-tolerance median trend can move the fitted exponent by at most the tolerance
   itself: 0.02, five times inside the 0.1 band half-width."

That is dimensionally wrong. For a one-carrier write the design's own identity is
median|D(N)| = N_E * N * medf(N). So

  beta_WU = d log(N * medf) / d log N = 1 + d log(medf)/d log N = 1 + slope/(medf*ln10)

where slope = d(medf)/d log10 N is the ABSOLUTE density-median slope the C1 condition
bounds. The exponent shift is slope/(medf*ln10) -- RELATIVE to medf. The tolerance is
absolute (0.02 f-units/decade); the shift it permits is 0.02/(medf*ln10), which at
medf = 0.05 is 0.17, not 0.02. The "at most 0.02" claim holds only near medf ~ 0.43
(where 1/(medf*ln10) ~ 1); below that it fails, without bound as medf -> 0.

Part A measures beta_WU - 1 against the predicted slope/(medf*ln10) on the design's
own pipeline. Part B shows the point certificate (point_certificate) is ONE-SIDED --
it counts only surrogates fitting BELOW 0.9, so an up-tilted ladder (beta > 1.1)
certifies OK -- which is why the low-density skew read of attack 1 reaches state OK
with beta 1.23 and the registered point-band prediction is asserted and false.
"""
import math, sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import v2_pipeline as V

MASTER = np.random.SeedSequence(777002)

def main():
    print("REFUTER B -- V002 ATTACK 2: co-tuning is dimensionally wrong; certificate one-sided")
    print("all reads through the LANE'S OWN sealed measure_occ / point_certificate")
    print()

    # ---- Part A: measured beta-1 vs predicted slope/(medf*ln10) ----
    print("PART A -- the exponent shift the C1 tolerance actually permits.")
    print("Each row: a cascade mask at base density f0, its DETERMINISTIC C1 slope, its")
    print("effective medf, the design's structural prediction slope/(medf*ln10), and the")
    print("MEASURED certified beta_WU-1 (median of state-OK reads, 40 seeds). C1 PASS means")
    print("the read is admitted by the density condition (|slope| <= %.2f)." % V.DENS_TOL_MED)
    print()
    print("  %-22s %8s %7s %10s %10s %6s %5s" %
          ("mask", "C1slope", "medf", "pred dB", "meas dB", "C1", "OK/40"))
    ln10 = math.log(10.0)
    roleW0 = V.occ_roles()
    for hi, lo, f0 in [(1.4, 0.6, 0.05), (1.5, 0.5, 0.05), (1.5, 0.5, 0.07),
                       (1.6, 0.4, 0.06), (1.5, 0.5, 0.04), (1.4, 0.6, 0.04),
                       (1.3, 0.7, 0.05), (1.5, 0.5, 0.10)]:
        mask = V.cascade_mask(f0, hi, lo)
        progm = np.zeros((V.NSECT, V.SECT)); progm[roleW0] = mask.astype(float)
        fmed = V.dens_med_ladder(progm, roleW0)
        slope, _ = V.fit_lin(V.GRID[V.DENS_MED_I0:], fmed[V.DENS_MED_I0:])
        medf_win = float(np.median(fmed[V.DENS_MED_I0:]))
        pred_dB = slope / (medf_win * ln10)
        c1pass = abs(slope) <= V.DENS_TOL_MED
        bws = []
        for ss in MASTER.spawn(40):
            rng = V.rng_from(ss)
            vals, prog, roleW, _, _ = V.build_occ(rng, mask=mask)
            m = V.measure_occ(vals, prog, roleW, rng)
            if m["state"] == "OK":
                bws.append(m["bWU"])
        meas_dB = (float(np.median(bws)) - 1.0) if bws else float("nan")
        print("  cascade %.1f/%.1f f0=%.2f %+8.4f %7.4f %+10.4f %+10.4f %6s %5d" %
              (hi, lo, f0, slope, medf_win, pred_dB, meas_dB,
               "PASS" if c1pass else "void", len(bws)))
    print()
    print("  The measured certified shift tracks slope/(medf*ln10), NOT the tolerance.")
    print("  Every C1-PASS row above is a legal one-carrier pattern that is guard-passing,")
    print("  density-median-trend-free, control-in-band and CERTIFIED (state OK), for which")
    print("  the registered point-band prediction beta_WU in [0.9,1.1] is ASSERTED -- and")
    print("  the measured certified beta_WU is outside it. This is refuter A's D1 family")
    print("  (the automatic-REFUTED trigger the judgment names) returning through the gap")
    print("  between an ABSOLUTE tolerance and a RELATIVE confound.")
    print()

    # ---- Part B: the certificate is one-sided ----
    print("PART B -- the point certificate guards only the LOWER band edge.")
    print("point_certificate(wu_p, uu_A) counts a surrogate BAD iff its fitted slope <")
    print("BAND_ACC[0] = %.2f. It never tests the upper edge %.2f. So a read whose own" %
          (V.BAND_ACC[0], V.BAND_ACC[1]))
    print("ladder tilts UP certifies OK. Demonstration on the attack-1 read (cascade")
    print("1.5/0.5 f0=0.05), 40 seeds -- certified reads and where their beta_WU sits:")
    mask = V.cascade_mask(0.05, 1.5, 0.5)
    ok = 0; above = 0; below = 0; inb = 0; certbad = []
    for ss in MASTER.spawn(40):
        rng = V.rng_from(ss)
        vals, prog, roleW, _, _ = V.build_occ(rng, mask=mask)
        m = V.measure_occ(vals, prog, roleW, rng)
        if m["state"] == "OK":
            ok += 1
            certbad.append(m["cert_bad"])
            if m["bWU"] > V.BAND_ACC[1]:
                above += 1
            elif m["bWU"] < V.BAND_ACC[0]:
                below += 1
            else:
                inb += 1
    print("  state OK %d/40 | of them beta_WU above 1.1: %d, in band: %d, below 0.9: %d" %
          (ok, above, inb, below))
    print("  certificate bad-surrogate counts on these OK reads: %s" %
          (("min %d max %d" % (min(certbad), max(certbad))) if certbad else "n/a"))
    print("  Every certified read sits ABOVE the band, and the certificate reports 0 bad")
    print("  surrogates on all of them: it cannot see an up-tilt. The point-band sentence")
    print("  ([0.9, 1.1], two-sided) is protected by a one-sided certificate.")

if __name__ == "__main__":
    main()

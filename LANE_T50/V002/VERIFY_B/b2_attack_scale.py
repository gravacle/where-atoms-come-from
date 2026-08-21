#!/usr/bin/env python3
"""REFUTER B (instrument) -- V002 ATTACK 1: THE C1 TOLERANCE DOES NOT SCALE.

The registered density condition (C1) bounds the ABSOLUTE trend of the median
per-block programmed fraction: |slope of median f-hat on log10 N| <= 0.02 (f units
per decade). But the exponent corruption is RELATIVE: for a one-carrier write
median|D| = N_E * N * medf(N), so

    delta beta  ~=  slope / (medf * ln 10).

The tolerance is absolute in f while the confound is relative -- so the SAME skew
structure that V002 voids at f0 = 0.5 (cascade slopes 0.036..0.088, all
VOID_DENSITY_MEDIAN at 1.000) sails UNDER the tolerance when the declared pattern's
base density is scaled down (slope scales with f0, the relative trend -- and the
exponent corruption -- does not).

Attack: the lane's OWN cascade/two-level mask builders (declared literal masks,
identical machinery, mask seed 20260821) at low base density, run through the
lane's OWN sealed measure_occ. For each read: state, density-median slope,
beta_WU, certificate outcome. Then the two registered kill criteria:
  (i)  the POINT-BAND sentence: a read in state OK (certified) with beta_WU
       outside [0.9, 1.1] on a legal one-carrier pattern = the registered
       prediction false on correct physics. NOTE the certificate is ONE-SIDED
       (point_certificate counts only surrogates fitting BELOW 0.9): an up-tilted
       ladder certifies.
  (ii) CLAUSE (b) as registered ("two declared data patterns on the SAME part,
       both reads state OK, beta_WU differ by more than 0.2 beyond 2 SE of the
       difference"): the mask read against a uniform f=0.5 read on the same part
       geometry = a falsifier firing on two legal one-carrier patterns.
"""
import math, sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import v2_pipeline as V

SEEDS = 30
MASTER = np.random.SeedSequence(777001)

def run_mask(mask, seeds=SEEDS):
    rows = []
    for ss in MASTER.spawn(seeds):
        rng = V.rng_from(ss)
        vals, prog, roleW, _, _ = V.build_occ(rng, mask=mask)
        m = V.measure_occ(vals, prog, roleW, rng)
        rows.append(m)
    return rows

def run_uniform(f, seeds=SEEDS):
    rows = []
    for ss in MASTER.spawn(seeds):
        rng = V.rng_from(ss)
        vals, prog, roleW, _, _ = V.build_occ(rng, f=f)
        m = V.measure_occ(vals, prog, roleW, rng)
        rows.append(m)
    return rows

def main():
    print("REFUTER B -- V002 ATTACK 1: C1 tolerance is absolute, the confound is relative")
    print("masks built by the LANE'S OWN builders (mask seed %d); measured by the lane's" % V.MASK_SEED)
    print("own sealed measure_occ; %d seeds per configuration; attack master seed 777001" % SEEDS)
    print("tolerance DENS_TOL_MED = %.2f (absolute f units per decade)" % V.DENS_TOL_MED)
    print()

    # deterministic mask slopes first (the C1 condition object itself)
    print("MASK LADDER -- deterministic density-median slope by base density f0:")
    print("  (V002's sealed R2 measured these same shapes at f0=0.5: slopes")
    print("   0.0361..0.0878, all VOID at 1.000. Watch the slope scale with f0.)")
    configs = []
    for hi, lo in ((1.5, 0.5), (1.6, 0.4), (1.7, 0.3)):
        for f0 in (0.5, 0.20, 0.10, 0.05):
            mask = V.cascade_mask(f0, hi, lo)
            progm = np.zeros((V.NSECT, V.SECT))
            roleW0 = V.occ_roles()
            progm[roleW0] = mask.astype(float)
            fmed = V.dens_med_ladder(progm, roleW0)
            ds, _ = V.fit_lin(V.GRID[V.DENS_MED_I0:], fmed[V.DENS_MED_I0:])
            passes = abs(ds) <= V.DENS_TOL_MED
            print("  cascade %.2f/%.2f f0=%.2f  f=%.4f  slope %+.4f  %s" %
                  (hi, lo, f0, mask.mean(), ds, "PASSES C1" if passes else "voided"))
            if passes and f0 < 0.5:
                configs.append(("cascade %.2f/%.2f f0=%.2f" % (hi, lo, f0), mask))
    print()

    # uniform partner (the ordinary declared pattern on the same part)
    uni = run_uniform(0.5)
    uni_ok = [m for m in uni if m["state"] == "OK"]
    print("UNIFORM f=0.50 PARTNER: OK %d/%d  beta_WU median %+.4f" %
          (len(uni_ok), SEEDS, float(np.median([m["bWU"] for m in uni_ok]))))
    print()

    print("THE GUARD-PASSING SKEW FAMILY -- full sealed pipeline per read:")
    kills_point = 0
    kills_b = 0
    for name, mask in configs:
        rows = run_mask(mask)
        states = [m["state"] for m in rows]
        ok = [m for m in rows if m["state"] == "OK"]
        seam = [m for m in rows if m["state"] == "SEAM"]
        print("  %-26s states %s" % (name, V._statecount(states)))
        if ok:
            bws = np.array([m["bWU"] for m in ok])
            slopes = [m["med_slope"] for m in ok]
            over = [m for m in ok if m["bWU"] > V.BAND_ACC[1]]
            under = [m for m in ok if m["bWU"] < V.BAND_ACC[0]]
            print("    OK (CERTIFIED) %d/%d  beta_WU median %+.4f  IQR [%+.4f, %+.4f]  max %+.4f" %
                  (len(ok), SEEDS, float(np.median(bws)),
                   float(np.percentile(bws, 25)), float(np.percentile(bws, 75)),
                   float(bws.max())))
            print("    density-median slope of these OK reads: median %+.4f (tolerance %.2f)" %
                  (float(np.median(slopes)), V.DENS_TOL_MED))
            print("    POINT-BAND VIOLATIONS among CERTIFIED reads: above 1.1: %d/%d, below 0.9: %d/%d" %
                  (len(over), len(ok), len(under), len(ok)))
            kills_point += len(over) + len(under)
            # clause (b) per the registered rule: pair each OK mask read with each OK uniform read
            fires = 0; total = 0; margins = []
            for mm in ok:
                for mu in uni_ok:
                    d = abs(mm["bWU"] - mu["bWU"])
                    se = math.hypot(mm["seWU"], mu["seWU"])
                    total += 1
                    margins.append(d - 2 * se - 0.2)
                    if d - 2 * se > 0.2:
                        fires += 1
            print("    CLAUSE (b) vs uniform f=0.5 (both OK, registered rule): FIRES %d/%d pairings"
                  "  fire margin median %+.4f" % (fires, total, float(np.median(margins))))
            kills_b += fires
        if seam:
            print("    SEAM %d  beta_WU median %+.4f" %
                  (len(seam), float(np.median([m["bWU"] for m in seam]))))
    print()
    print("SUMMARY: certified point-band violations %d; clause-(b) fires on legal" % kills_point)
    print("one-carrier pattern pairs %d. Certificate one-sidedness check: point_certificate" % kills_b)
    print("counts a surrogate bad only if its fitted slope < %.2f (the lower edge);" % V.BAND_ACC[0])
    print("an up-tilted ladder certifies OK by construction.")

if __name__ == "__main__":
    main()

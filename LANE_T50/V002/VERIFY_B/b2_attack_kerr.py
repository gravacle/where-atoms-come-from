#!/usr/bin/env python3
"""REFUTER B -- V002 ATTACK 3: the Kerr raw-value law (honest risk #1), and a reseed
robustness check of attack 1.

The design names polar Kerr microscopy and concedes (risk #1) that a MULTIPLICATIVE
state-correlated artifact could evade the additive DCF voider. Real magneto-optic
imaging carries exactly such terms: the measured analyzer-difference signal is
theta_K(M) plus a magnetization-INDEPENDENT reflectivity/ellipticity background that
differs between written and erased domain structure, and the Kerr rotation itself can
be weighted by local reflectivity that depends on the written pattern.

I test the concrete physical channel the design's own additive voider does NOT model:
a per-cell reflectivity weight correlated with the written BIT, measured = M*(1+g*b)
where b in {0,1} marks a domain in the "1" state. For a DATA block this adds
g*sum(b_i) ~ g*N*p1, a term LINEAR in N -- exactly the accumulation signature clause
(c) tests. Question: does the DCF voider catch it (INCONCLUSIVE, safe) or does clause
(c) FIRE on correctly-screening data (a false fire, a kill)?

Then Part C reseeds attack 1 under a fresh master seed to show the point-band
violation is not a seed artifact.
"""
import math, sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import v2_pipeline as V

MASTER = np.random.SeedSequence(777003)

def build_ori_kerr_mult(rng, gain, data_p1=0.5, nsect=V.NSECT_ORI):
    """Same as V.build_ori (kerr), but with a per-cell multiplicative reflectivity
    weight (1 + gain*b_i) where b_i = 1 for a cell in the '1'/up domain -- a
    magnetization-magnitude-independent, PATTERN-correlated Kerr weighting. Everything
    else identical to the lane's kerr model (grain average is trivial here since we
    weight the resolved cell; PSF, read noise, edge margin all follow the lane)."""
    classes = [V.ORI_CYCLE[s % 6] for s in range(nsect)]
    M = np.empty((nsect, V.SECT)); B = np.empty((nsect, V.SECT))
    for s, c in enumerate(classes):
        if c == "DATA":
            b = (rng.random(V.SECT) < data_p1).astype(float)
        elif c == "DCF":
            b = np.tile([1.0, 0.0], V.SECT // 2)
        elif c == "DC":
            b = np.ones(V.SECT)
        else:  # U: AC-erased, half up half down at grain level -> cell mean ~0, b~0.5
            b = (rng.binomial(V.GRAINS, 0.5, V.SECT) / V.GRAINS)
        m = 2.0 * b - 1.0
        M[s] = m; B[s] = b
    meas = M * (1.0 + gain * B)          # the multiplicative pattern-correlated weight
    flat = np.convolve(meas.reshape(-1), V._psf_kernel(), mode="same")
    vals = flat.reshape(nsect, V.SECT) + rng.normal(0.0, V.READ_NOISE, (nsect, V.SECT))
    return vals, classes

def main():
    print("REFUTER B -- V002 ATTACK 3: Kerr multiplicative pattern-correlated weight")
    print("measured = M*(1 + g*b), b = up-domain indicator; read otherwise the lane's own")
    print("kerr model. Screening DATA (p1=0.5). Does the DCF voider catch it, or does")
    print("clause (c) FIRE on correct screening physics? 30 seeds per level.")
    print()
    print("  gain g | states | fire_c | beta_DATA-U med | beta_DCF-U med | beta_DC-U med")
    for g in (0.0, 0.05, 0.10, 0.20, 0.40):
        sts = []; fc = 0; bdu = []; bfu = []; bcu = []
        for ss in MASTER.spawn(30):
            rng = V.rng_from(ss)
            vals, classes = build_ori_kerr_mult(rng, g)
            m = V.measure_ori(vals, classes, rng)
            sts.append(m["state"]); fc += m.get("fire_c", False)
            if "bDU" in m: bdu.append(m["bDU"])
            if "bFU" in m: bfu.append(m["bFU"])
            if m.get("bCU") is not None and isinstance(m.get("bCU"), float) \
               and not math.isnan(m.get("bCU", float("nan"))):
                bcu.append(m["bCU"])
        def med(x): return ("%+.4f" % float(np.median(x))) if x else "  --  "
        print("  g=%.2f  %s  fire_c %d/30  %s  %s  %s" %
              (g, V._statecount(sts), fc, med(bdu), med(bfu), med(bcu)))
    print()
    print("  Reading: if the DCF voider (INCONCLUSIVE_CROSSTALK) catches it, the design is")
    print("  safe on this vector -- a named non-verdict. If clause (c) fires with state OK,")
    print("  it is a false fire on correct screening physics (the risk-#1 exposure realised).")
    print()

    # ---- Part C: reseed the attack-1 point-band violation ----
    print("PART C -- reseed of attack 1 (cascade 1.5/0.5 f0=0.05), fresh master seed,")
    print("60 seeds: the point-band violation is not a seed artifact.")
    mask = V.cascade_mask(0.05, 1.5, 0.5)
    ok = above = 0; bws = []
    for ss in MASTER.spawn(60):
        rng = V.rng_from(ss)
        vals, prog, roleW, _, _ = V.build_occ(rng, mask=mask)
        m = V.measure_occ(vals, prog, roleW, rng)
        if m["state"] == "OK":
            ok += 1; bws.append(m["bWU"])
            if m["bWU"] > V.BAND_ACC[1]:
                above += 1
    print("  state OK %d/60 | certified beta_WU above 1.1: %d/%d | median %+.4f | max %+.4f" %
          (ok, above, ok, float(np.median(bws)) if bws else float("nan"),
           float(np.max(bws)) if bws else float("nan")))

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""REFUTER B -- V002 ATTACK 1, AIRTIGHT: trace one low-density skew read through EVERY
registered gate, proving it satisfies the full point-band quantifier while beta_WU sits
outside [0.9, 1.1]. Pre-empts any 'it did not really pass gate X' rebuttal.

Also confirms the same shape at f0=0.5 (which V002's R2 tested) correctly VOIDS -- so
the design verified the high-density corner and the low-density uniform corner but
never the low-density SKEW corner where D1 survives.
"""
import math, sys, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import v2_pipeline as V

def trace(mask, label, seed):
    rng = V.rng_from(np.random.SeedSequence(seed))
    vals, prog, roleW, _, _ = V.build_occ(rng, mask=mask)
    # replay measure_occ's gate order with prints
    v = V.detrend_addr(vals, roleW)
    cs = V.seccs(v)
    uidx = list(np.where(~roleW)[0])
    upool, A_uu = V.uu_pool(cs, V.same_role_pairs(uidx), rng)
    print("  [%s] seed %d" % (label, seed))
    print("    GATE 1 railed?  A_UU(16) = %.4f  -> %s" %
          (A_uu[0], "RAILED" if A_uu[0] == 0 else "pass"))
    pairs = V.adj_wu_pairs(roleW)
    print("    pairs available: %d (need >= %d)" % (len(pairs), V.K_MIN))
    wu_p, A_wu, _ = V.pair_pool(cs, pairs, rng)
    i0 = V.guard_start(A_wu, A_uu)
    print("    GATE 2 void guard: first N_min with A_WU > %.0f*A_UU -> i0=%s (N=%s)  "
          "margin=%.2f  -> %s" %
          (V.KAPPA_VOID, i0, V.GRID[i0] if i0 is not None else None,
           (A_wu[i0] / A_uu[i0]) if i0 is not None else float("nan"),
           "READS_UNWRITTEN" if i0 is None else "pass"))
    if i0 is None:
        print("    -> state READS_UNWRITTEN (this seed); try another")
        return None
    fmed = V.dens_med_ladder(prog, roleW)
    dslope, _ = V.fit_lin(V.GRID[V.DENS_MED_I0:], fmed[V.DENS_MED_I0:])
    print("    GATE 3 density-median: slope on N>=64 = %+.4f  (tolerance %.2f)  -> %s" %
          (dslope, V.DENS_TOL_MED,
           "VOID_DENSITY_MEDIAN" if abs(dslope) > V.DENS_TOL_MED else "PASS"))
    if abs(dslope) > V.DENS_TOL_MED:
        print("    -> state VOID_DENSITY_MEDIAN")
        return None
    Ns = V.GRID[i0:]
    bWU, seWU = V.fit_loglog(Ns, A_wu[i0:])
    bUU, seUU = V.fit_loglog(Ns, A_uu[i0:])
    xi = bWU - bUU
    print("    GATE 4 control band: beta_UU = %+.4f in [%.2f, %.2f]?  -> %s" %
          (bUU, V.BAND_CTL[0], V.BAND_CTL[1],
           "pass" if V.BAND_CTL[0] <= bUU <= V.BAND_CTL[1] else "INCONCLUSIVE_CONTROL"))
    if not (V.BAND_CTL[0] <= bUU <= V.BAND_CTL[1]):
        print("    -> state INCONCLUSIVE_CONTROL")
        return None
    cp, cb = V.point_certificate(wu_p, A_uu, rng)
    state = "OK" if (cb == 0 and cp > 0) else "SEAM"
    print("    GATE 5 certificate: %d guard-passing surrogates, %d below-band  -> state %s"
          % (cp, cb, state))
    print("    ===> beta_WU = %+.4f  (registered point-band prediction: [%.2f, %.2f])" %
          (bWU, V.BAND_ACC[0], V.BAND_ACC[1]))
    print("         xi = %+.4f (>= %.2f: %s)" % (xi, V.XI_MIN, xi >= V.XI_MIN))
    inband = V.BAND_ACC[0] <= bWU <= V.BAND_ACC[1]
    print("         POINT-BAND PREDICTION HOLDS: %s  <-- %s" %
          (inband, "OK" if inband else "FALSE on this legal, certified, guard-passing "
                                       "one-carrier pattern"))
    return bWU

def main():
    print("REFUTER B -- V002: full gate trace of a low-density skew read")
    print()
    print("A. cascade 1.5/0.5 at f0=0.05 (C1 slope +0.0194, UNDER the 0.02 tolerance):")
    mask_lo = V.cascade_mask(0.05, 1.5, 0.5)
    got = None
    for sd in (101, 102, 103, 104, 105):
        b = trace(mask_lo, "cascade 1.5/0.5 f0=0.05", sd)
        print()
        if b is not None and got is None:
            got = b
    print("B. THE SAME SHAPE at f0=0.5 (what V002's R2 tested) -- correctly voided:")
    mask_hi = V.cascade_mask(0.5, 1.5, 0.5)
    trace(mask_hi, "cascade 1.5/0.5 f0=0.50", 201)
    print()
    print("The design verified skew masks VOID at f0=0.5 (R2) and low-density UNIFORM")
    print("patterns PASS (R2 honest table, f=0.05 uniform slope 0.0019). It never tested")
    print("the low-density SKEW corner, where the median trend falls under the absolute")
    print("0.02 tolerance while the relative confound (slope/(medf*ln10)) does not -- and")
    print("the point-band prediction is asserted and false there.")

if __name__ == "__main__":
    main()

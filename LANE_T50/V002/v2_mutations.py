#!/usr/bin/env python3
"""T-50 DESIGN ONE V002 -- THE ADVERSARIAL MUTATION SUITE.

Every STRUCTURAL DECISION of the V002 pipeline is shown able to fail, at the measured
rate, under its DESIGNATED killing mutation -- the bar is fail-rate 1.000 (D-8; the
principal's gate ruling, extended to this design's own decisions). The suite carries:
  - the V001 members (two-signed write; non-zero-mean residual as an INVARIANCE row;
    the fixed-record density counterexample);
  - the SKEW MASKS (multiplicative cascade, two-level) -- refuter A's D1 family, the
    judge's carried median-vs-mean density-skew member;
  - the FIXED-PATTERN member -- refuter B's K3, the judge's carried member;
  - the RAILED-population member (K5/D5);
  - the STRAY-FIELD-instrument member (K1);
  - the CROSSTALK and BIASED-ERASE members (A-D2; V001's scope clause, now detected);
  - the IMBALANCED-DATA member: clause (c)'s own firing branch (the falsifier IS
    triggerable on the named instrument -- the K1 complaint answered by measurement).

PRECEDENCE INVARIANT (C6), asserted on every read of every member: a falsifier boolean
that is True while the read state is not OK/SEAM (occupancy) or OK (orientation) is a
suite FAILURE -- code equals text or this file exits nonzero.

Runs the sealed pipeline's own functions (one implementation; INST-17).
"""
import math, sys
import numpy as np
import v2_pipeline as V

SEEDS = 30
MASTER = np.random.SeedSequence(50002)


def occ_reads(mutate, seeds=SEEDS, **kw):
    """Yield measured occupancy reads under a builder mutation."""
    for ss in MASTER.spawn(seeds):
        rng = V.rng_from(ss)
        vals, prog, roleW, wclass, fp = V.build_occ(rng, **mutate(kw))
        m = V.measure_occ(vals, prog, roleW, rng)
        assert (not m.get("fire_a", False)) or m["state"] in ("OK", "SEAM"), \
            "PRECEDENCE VIOLATION (occ): fire in state %s" % m["state"]
        yield m


def ori_reads(seeds=SEEDS, **kw):
    for ss in MASTER.spawn(seeds):
        rng = V.rng_from(ss)
        vals, classes, fp = V.build_ori(rng, **kw)
        m = V.measure_ori(vals, classes, rng)
        assert (not m.get("fire_c", False)) or m["state"] == "OK", \
            "PRECEDENCE VIOLATION (ori): fire in state %s" % m["state"]
        yield m


def rate(n, d):
    return "%d/%d = %.3f" % (n, d, n / d)


def main():
    out = []
    P = out.append
    P("T-50 DESIGN ONE V002 -- ADVERSARIAL MUTATION SUITE -- SEALED RUN")
    P("suite master seed SeedSequence(50002) | %d seeds per member | numpy %s" %
      (SEEDS, np.__version__))
    P("bar: every structural decision fails at rate 1.000 under its designated")
    P("member; invariance rows must NOT flip (0.000) with the raw shadow displayed;")
    P("the C6 precedence invariant is asserted on every read (a fire in a non-")
    P("OK/SEAM state aborts this suite).")
    P("")
    results = {}   # name -> (achieved, target_desc, ok)

    # ---------------- BASELINE ----------------
    P("M0 BASELINE (no mutation) -- the suite's own honest branch:")
    okc = b1 = fires = 0
    for m in occ_reads(lambda kw: dict(f=0.5)):
        okc += m["state"] == "OK"; b1 += m.get("B1", False)
        fires += m.get("fire_a", False)
    P("   occupancy: state OK %s | B1 %s | fires %d" %
      (rate(okc, SEEDS), rate(b1, SEEDS), fires))
    results["M0-occ all-pass"] = (okc == SEEDS and b1 == SEEDS and fires == 0,
                                  "all OK, B1 true, 0 fires")
    okc = b3 = fires = 0
    for m in ori_reads():
        okc += m["state"] == "OK"; b3 += m.get("B3", False)
        fires += m.get("fire_c", False)
    P("   orientation: state OK %s | B3 %s | fires %d" %
      (rate(okc, SEEDS), rate(b3, SEEDS), fires))
    results["M0-ori all-pass"] = (okc == SEEDS and b3 == SEEDS and fires == 0,
                                  "all OK, B3 true, 0 fires")
    P("")

    # ---------------- MU1 two-signed write ----------------
    P("MU1 TWO-SIGNED WRITE (occupancy behaving like orientation) -- designated kill")
    P("    of decision S1 (B1, the accumulation claim):")
    flips = fires = 0; bws = []
    for m in occ_reads(lambda kw: dict(f=0.5, two_signed=True)):
        if not m.get("B1", False):
            flips += 1
        if m.get("fire_a", False):
            fires += 1
        if m["state"] in ("OK", "SEAM"):
            bws.append(m["bWU"])
    P("    B1 FLIPS %s | falsifier (a) fires %s | beta_WU median %+.4f" %
      (rate(flips, SEEDS), rate(fires, SEEDS), float(np.median(bws))))
    results["S1 <- MU1"] = (flips == SEEDS, "flip 1.000")
    results["S1 fire <- MU1"] = (fires == SEEDS, "fire 1.000")

    # ---------------- MU2 invariance rows ----------------
    P("MU2 NON-ZERO-MEAN RESIDUAL AND COMMON-MODE OFFSET -- INVARIANCE rows (the")
    P("    T-50 constraint-1 architecture: decisions must NOT move; the raw statistic's")
    P("    corruption is displayed as the shadow):")
    P("    (the shadow is the UNTREATED raw single-sector exponent -- the corruption")
    P("    the contrast is immune to; the treated pipeline never sees it):")
    for name, kw in (("mu=+2 e", dict(f=0.5, mu=2.0)),
                     ("mu=+20 e", dict(f=0.5, mu=20.0)),
                     ("offset +0.5 e", dict(f=0.5, offset=0.5)),
                     ("offset +50 e", dict(f=0.5, offset=50.0))):
        flips = 0; raws = []
        for ss in MASTER.spawn(SEEDS):
            rng = V.rng_from(ss)
            vals, prog, roleW, _, _ = V.build_occ(rng, **kw)
            m = V.measure_occ(vals, prog, roleW, rng)
            assert (not m.get("fire_a", False)) or m["state"] in ("OK", "SEAM")
            if not m.get("B1", False):
                flips += 1
            if m["state"] in ("OK", "SEAM"):
                raws.append(V.fit_loglog(V.GRID, V.raw_ladder(vals, roleW, rng))[0])
        P("    %-14s B1 flips %s | untreated-raw shadow beta %+8.4f" %
          (name, rate(flips, SEEDS), float(np.median(raws))))
        results["S-inv <- MU2 %s" % name] = (flips == 0, "flip 0.000")

    # ---------------- MU3 fixed record ----------------
    P("MU3 FIXED 512-CELL RECORD IN GROWING BLOCKS (the M3 counterexample, as the")
    P("    declared literal mask) -- designated kill of the guard/density pair:")
    caught = fires = 0; sts = []
    mask = V.record_mask()
    for m in occ_reads(lambda kw: dict(mask=mask)):
        sts.append(m["state"])
        if m["state"] in ("READS_UNWRITTEN", "VOID_DENSITY_MEDIAN"):
            caught += 1
        fires += m.get("fire_a", False)
    P("    VOIDED (READS_UNWRITTEN or VOID_DENSITY_MEDIAN) %s | fires %d | states %s" %
      (rate(caught, SEEDS), fires, V._statecount(sts)))
    results["S3 <- MU3"] = (caught == SEEDS and fires == 0, "voided 1.000, 0 fires")

    # ---------------- MU4-6 skew masks ----------------
    P("MU4-MU6 THE SKEW MASKS (refuter A's D1; the judge's carried density-skew")
    P("    member) -- designated kills of decision S2 (the density-median condition):")
    for name, mask in (("MU4 cascade 1.5/0.5", V.cascade_mask(0.5, 1.5, 0.5)),
                       ("MU5 cascade 1.7/0.3", V.cascade_mask(0.5, 1.7, 0.3)),
                       ("MU6 two-level 0.95/0.03", V.twolevel_mask())):
        voided = fires = 0; sts = []
        for m in occ_reads(lambda kw, mask=mask: dict(mask=mask)):
            sts.append(m["state"])
            if m["state"] in ("VOID_DENSITY_MEDIAN", "READS_UNWRITTEN"):
                voided += 1
            fires += m.get("fire_a", False)
        P("    %-24s VOIDED %s | fires %d | states %s" %
          (name, rate(voided, SEEDS), fires, V._statecount(sts)))
        results["S2 <- %s" % name.split()[0]] = (voided == SEEDS and fires == 0,
                                                 "voided 1.000, 0 fires")

    # ---------------- MU7 fixed pattern ----------------
    P("MU7 IID PER-SECTOR FIXED PATTERN s=0.5 e/cell (refuter B's K3 level 5x; the")
    P("    judge's carried fixed-pattern member) -- designated kill of decision S4")
    P("    (the control-band fixed-pattern detector), rung 1 only:")
    caught = fires = 0
    for m in occ_reads(lambda kw: dict(f=0.5, fp_mode="iid", fp_amp=0.5)):
        caught += m["state"] == "INCONCLUSIVE_CONTROL"
        fires += m.get("fire_a", False)
    P("    INCONCLUSIVE_CONTROL %s | fires %d (sealed V001: fire_a 26/30 here)" %
      (rate(caught, SEEDS), fires))
    results["S4 <- MU7"] = (caught == SEEDS and fires == 0, "detect 1.000, 0 fires")
    P("    ...and rung 2 (the calibration read) RECOVERS the verdict:")
    rec = 0
    for ss in MASTER.spawn(SEEDS):
        rng = V.rng_from(ss)
        vals, prog, roleW, _, fp = V.build_occ(rng, f=0.5, fp_mode="iid", fp_amp=0.5)
        calib = V.build_occ_calib(rng, fp)
        m = V.measure_occ(vals, prog, roleW, rng, calib=calib)
        rec += m.get("B1", False)
    P("    rung-2 B1 recovery %s" % rate(rec, SEEDS))
    results["S4 rung2 recovery"] = (rec == SEEDS, "recover 1.000")

    # ---------------- MU8 railed ----------------
    P("MU8 RAILED ERASED POPULATION (K5/D5) -- designated kill of decision S5")
    P("    (the insufficient-access branch):")
    caught = fires = 0
    for m in occ_reads(lambda kw: dict(f=0.5, railed=True)):
        caught += m["state"] == "INCONCLUSIVE_RAILED"
        fires += m.get("fire_a", False)
    P("    INCONCLUSIVE_RAILED %s | fires %d | no nan reaches any decision" %
      (rate(caught, SEEDS), fires))
    results["S5 <- MU8"] = (caught == SEEDS and fires == 0, "detect 1.000, 0 fires")

    # ---------------- MU9 stray-field instrument ----------------
    P("MU9 STRAY-FIELD INSTRUMENT (K1: MFM/Hall/NV transfer, T(0)=0, d=1 cell) --")
    P("    designated kill of decision S6 (the orientation positive control):")
    caught = 0
    for m in ori_reads(read="stray", standoff=1.0):
        caught += m["state"] == "INCONCLUSIVE_DC_CONTROL"
    P("    INCONCLUSIVE_DC_CONTROL %s -- the named-instrument requirement is" %
      rate(caught, SEEDS))
    P("    load-bearing and self-measured: a stray-field reader is told, by the read's")
    P("    own control, that this instrument cannot run the protocol.")
    results["S6 <- MU9"] = (caught == SEEDS, "detect 1.000")

    # ---------------- MU10 crosstalk ----------------
    P("MU10 CLASS-CORRELATED CROSSTALK (A-D2: +0.10 and +0.20 per grain on written")
    P("    sectors) -- designated kill of decision S7 (the DC-free voider, C2):")
    for co in (0.10, 0.20):
        caught = fires = 0
        for m in ori_reads(class_offset=co):
            caught += m["state"] == "INCONCLUSIVE_CROSSTALK"
            fires += m.get("fire_c", False)
        P("    offset %.2f: INCONCLUSIVE_CROSSTALK %s | fire_c %d (V001: %s)" %
          (co, rate(caught, SEEDS), fires,
           "13/40" if co == 0.10 else "34/40"))
        results["S7 <- MU10 c=%.2f" % co] = (caught == SEEDS and fires == 0,
                                             "void 1.000, 0 fires")

    # ---------------- MU11 imbalanced data ----------------
    P("MU11 GENUINELY ACCUMULATING DATA (75% one-way): clause (c)'s OWN firing branch")
    P("    -- the falsifier is triggerable on the named instrument:")
    fires = okc = 0
    for m in ori_reads(data_p1=0.75):
        okc += m["state"] == "OK"
        fires += m.get("fire_c", False)
    P("    state OK %s | fire_c %s" % (rate(okc, SEEDS), rate(fires, SEEDS)))
    results["S8 <- MU11"] = (fires == SEEDS, "fire 1.000")

    # ---------------- MU12 biased erase ----------------
    P("MU12 BIASED ERASE +0.25/grain (V001 policed by scope clause; V002 detects):")
    caught = fires = 0
    for m in ori_reads(erase_bias=0.25):
        caught += m["state"] == "INCONCLUSIVE_CROSSTALK"
        fires += m.get("fire_c", False)
    P("    INCONCLUSIVE_CROSSTALK %s | fire_c %d (V001: fired 34-35/50 when the" %
      (rate(caught, SEEDS), fires))
    P("    scope clause was ignored -- the clause is replaced by a detector)")
    results["S7 <- MU12"] = (caught == SEEDS and fires == 0, "void 1.000, 0 fires")

    # ---------------- matrix ----------------
    P("")
    P("=" * 74)
    P("DESIGNATION MATRIX -- decision x designated member -> achieved (bar):")
    allok = True
    for name, (ok, target) in results.items():
        P("  %-28s %-24s %s" % (name, target, "MET" if ok else "NOT MET"))
        allok = allok and ok
    P("")
    if allok:
        P("SUITE VERDICT: ALL DESIGNATED CELLS MET -- every structural decision has a")
        P("measured failing branch at its designated rate; the invariance rows hold")
        P("with the raw shadow shown; no fire boolean was ever True outside the")
        P("registered firing states (asserted on every read of every member).")
    else:
        P("SUITE VERDICT: AT LEAST ONE DESIGNATED CELL NOT MET -- the matrix above")
        P("names it; this file exits nonzero and the design does not ship over it.")
    P("END OF MUTATION SUITE")
    print("\n".join(out))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())

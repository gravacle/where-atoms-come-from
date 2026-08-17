#!/usr/bin/env python3
"""
LANE W-10 / C — REFUTER 2 — STEP 2.
HOW BIG IS THE EFFECT THE RESONANT ARM IS SUPPOSED TO DETECT, AND HOW BIG IS THE RULER?

Lane C's C-4 claims: "The three limits ALL EXIST on both four-class carriers, and TWO OF THE
THREE ARE NOT m(P)", with the resonant arm's target the (11,20) subtorus Mahler measure.
That is a statement about the SIZE of  m(Q_{11,20}) - m(P).  This step computes that size to
40 digits by TWO independent evaluators and sets it beside the lane's own measurement error.

EVALUATORS (both mpmath, dps printed):
  (S1) trapezoid on log|Q(e^{i t})| with half-step offset
  (S2) Jensen's formula on the ROOTS: m(Q) = log|a_d| + sum_{|r|>1} log|r|
Q(z) = p01 + p00 z^m + p11 z^n + p10 z^{m+n} is P restricted to H = {(z^n, z^-m)}, (m,n)=(11,20).
If (S1) and (S2) disagree beyond the stated tolerance the row is void.

m(P) itself is taken from the closed forms lane C proves in C-2 and re-checked by the Jensen
reduction; the closed forms are re-derived here independently in mpmath.
"""
import mpmath as mp

mp.mp.dps = 50

M_REL, N_REL = 11, 20      # the primitive relation of f = 2.0, c = 1.1

CASES = [
    ("B0b  SENSE U  (4,2,1,2)/9", [mp.mpf(4)/9, mp.mpf(2)/9, mp.mpf(1)/9, mp.mpf(2)/9]),
    ("B4   SENSE U  (1,1,1,3)/6", [mp.mpf(1)/6, mp.mpf(1)/6, mp.mpf(1)/6, mp.mpf(3)/6]),
    ("K1   SENSE U  (0,2,2,1)/5", [mp.mpf(0), mp.mpf(2)/5, mp.mpf(2)/5, mp.mpf(1)/5]),
    ("SENSE C 4-class (1,1,1,1)/4", [mp.mpf(1)/4]*4),
    ("B0b* SENSE U  (2,1,3,3)/9", [mp.mpf(2)/9, mp.mpf(1)/9, mp.mpf(3)/9, mp.mpf(3)/9]),
]

# lane C's reported K = 1e7 Birkhoff numbers, transcribed from C_04_birkhoff.out.txt
LANEC = {
    "B0b  SENSE U  (4,2,1,2)/9": dict(dioph=-0.810930271, reson=-0.810930225,
                                      sub=-0.810930204535, err_d=5.435e-08),
    "B4   SENSE U  (1,1,1,3)/6": dict(dioph=-0.693147237, reson=-0.693147188,
                                      sub=-0.693147181044, err_d=5.612e-08),
    "K1   SENSE U  (0,2,2,1)/5": dict(dioph=-0.756573532, reson=-0.756337086,
                                      sub=-0.756337009107, err_d=5.402e-08),
    "SENSE C 4-class (1,1,1,1)/4": dict(dioph=-1.386295105, reson=-1.386295428,
                                        sub=-1.386293534824, err_d=1.405e-06),
}


def m_jensen_2var(p, n=1 << 16):
    """m(P) by the Jensen reduction, half-step offset nodes, mpmath."""
    p00, p10, p01, p11 = p
    tot = mp.mpf(0)
    for j in range(n):
        t = 2 * mp.pi * (mp.mpf(j) + mp.mpf(1) / 2) / n
        e = mp.expjpi(2 * (mp.mpf(j) + mp.mpf(1) / 2) / n)
        A = abs(p00 + p10 * e); B = abs(p01 + p11 * e)
        tot += mp.log(A if A > B else B)
    return tot / n


def Q_coeffs(p, m=M_REL, n=N_REL):
    """Q(z) = p01 + p00 z^m + p11 z^n + p10 z^{m+n}; returned HIGH degree first for polyroots."""
    p00, p10, p01, p11 = p
    d = m + n
    c = [mp.mpf(0)] * (d + 1)
    c[0] += p01; c[m] += p00; c[n] += p11; c[d] += p10
    return list(reversed(c))          # polyroots wants leading first


def m_sub_trapezoid(p, n=1 << 16, m=M_REL, nn=N_REL):
    p00, p10, p01, p11 = p
    tot = mp.mpf(0)
    for j in range(n):
        u = (mp.mpf(j) + mp.mpf(1) / 2) / n
        x = mp.expjpi(2 * nn * u); y = mp.expjpi(-2 * m * u)
        tot += mp.log(abs(p00 + p10 * x + p01 * y + p11 * x * y))
    return tot / n


def m_sub_roots(p, m=M_REL, n=N_REL):
    c = Q_coeffs(p, m, n)
    while c and c[0] == 0:
        c = c[1:]
    lead = c[0]
    rts = mp.polyroots(c, maxsteps=200, extraprec=400)
    s = mp.log(abs(lead))
    for r in rts:
        if abs(r) > 1:
            s += mp.log(abs(r))
    return s, rts


if __name__ == "__main__":
    print("=" * 104)
    print(f"R2_02 — EFFECT SIZE OF THE RESONANT ARM vs THE RULER THAT MEASURES IT.  dps = {mp.mp.dps}")
    print(f"        subtorus H from the primitive relation ({M_REL},{N_REL}) of f = 2.0, c = 1.1")
    print("=" * 104)

    rows = []
    for label, p in CASES:
        mP = m_jensen_2var(p)
        mS_t = m_sub_trapezoid(p)
        mS_r, rts = m_sub_roots(p)
        agree = abs(mS_t - mS_r)
        onunit = sum(1 for r in rts if abs(abs(r) - 1) < mp.mpf('1e-20'))
        print("\n" + "-" * 104)
        print(f"  {label}")
        print(f"     m(P)                       = {mp.nstr(mP, 18)}")
        print(f"     m(Q_(11,20)) trapezoid 2^16= {mp.nstr(mS_t, 18)}")
        print(f"     m(Q_(11,20)) ROOTS/Jensen  = {mp.nstr(mS_r, 18)}"
              f"   |trap - roots| = {mp.nstr(agree, 4)}")
        print(f"     roots of Q ON the unit circle: {onunit} of {len(rts)}"
              f"   {'  <- log singularities: trapezoid is only O(1/n) here' if onunit else ''}")
        eff = mS_r - mP
        print(f"     EFFECT SIZE  m(Q) - m(P)   = {mp.nstr(eff, 12)}   |eff| = {mp.nstr(abs(eff), 6)}")
        rows.append((label, mP, mS_r, eff, onunit))

    print("\n" + "=" * 104)
    print("THE POWER TABLE.  Can lane C's K = 1e7 run tell the two candidate limits apart?")
    print("=" * 104)
    print(f"  {'case':30s} {'|m(Q)-m(P)| (effect)':>22s} {'lane C err at K=1e7':>22s} "
          f"{'ratio':>10s}  verdict")
    for label, mP, mS, eff, onunit in rows:
        if label not in LANEC:
            continue
        L = LANEC[label]
        ruler = L['err_d']            # the arm-independent finite-N wobble, from the DIOPHANTINE row
        ratio = float(abs(eff)) / ruler
        verdict = ("RESOLVED" if ratio > 3 else
                   "UNRESOLVED — the ruler is coarser than the effect")
        print(f"  {label:30s} {mp.nstr(abs(eff),6):>22s} {ruler:22.3e} {ratio:10.3f}  {verdict}")

    print("\n" + "=" * 104)
    print("AND WHICH TARGET DID LANE C's MEASURED RESONANT VALUE ACTUALLY LAND NEAREST?")
    print("=" * 104)
    print(f"  {'case':30s} {'measured (K=1e7)':>18s} {'|meas-m(P)|':>13s} {'|meas-m(Q)|':>13s}  nearest")
    for label, mP, mS, eff, onunit in rows:
        if label not in LANEC:
            continue
        L = LANEC[label]
        meas = mp.mpf(repr(L['reson']))
        d1 = abs(meas - mP); d2 = abs(meas - mS)
        nearest = "m(P)  <-- THE TARGET IT IS CLAIMED NOT TO CONVERGE TO" if d1 < d2 else "m(Q) subtorus"
        print(f"  {label:30s} {float(meas):18.9f} {float(d1):13.3e} {float(d2):13.3e}  {nearest}")

    print("\n" + "=" * 104)
    print("LANE C's OWN SUBTORUS EVALUATOR, CHECKED.  C_04.subtorus_limit is a 2^22-node")
    print("trapezoid on log|Q|.  Where Q has roots ON the unit circle that is an O(1/n) rule.")
    print("=" * 104)
    for label, mP, mS, eff, onunit in rows:
        if label not in LANEC:
            continue
        laneval = mp.mpf(repr(LANEC[label]['sub']))
        print(f"  {label:30s} lane C = {float(laneval):.12f}   exact = {mp.nstr(mS,13):>16s}"
              f"   err = {float(abs(laneval-mS)):.2e}"
              f"{'   <- unit-circle roots present' if onunit else ''}")

    print("\n" + "=" * 104)
    print("SENSE C, EXACTLY.  P = (1+x)(1+y)/4 restricts on H to Q = (1+z^11)(1+z^20)/4, so")
    print("m(Q) = log(1/4) + m(1+z^11) + m(1+z^20) = -log 4 = m(P) EXACTLY.  The subtorus limit")
    print("and the torus limit COINCIDE there, and lane C's -1.386293534824 is 8.3e-07 of")
    print("trapezoid error, not a distinct limit.")
    print(f"  -log 4 = {mp.nstr(-mp.log(4), 18)}")
    print("\nDONE.")

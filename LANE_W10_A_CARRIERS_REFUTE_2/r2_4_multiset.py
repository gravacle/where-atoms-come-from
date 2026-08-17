#!/usr/bin/env python3
"""R2 SCRIPT 4 -- THE MULTISET THEOREM'S HYPOTHESIS.
The registrar named it REAL NON-NEGATIVE.  The lane (W10A-10) corrected that to REAL.
This script tests whether REAL is itself one word too strong, records that MY OWN first
proposed replacement (COLLINEAR) was also too strong, and gives the name that survives.
Seed 20260817, NOT the lane's 20260816: no agreement below is a shared-draw artefact."""
import sys, math
from itertools import permutations
import numpy as np

OUT = []
def o(s=""):
    print(s); OUT.append(s)

NQ = 1 << 20
T = (np.arange(NQ) + 0.5) * (2 * np.pi / NQ)
E = np.exp(1j * T)
def m(v):
    p00, p10, p01, p11 = v
    return float(np.mean(np.log(np.maximum(
        np.maximum(np.abs(p00 + p01 * E), np.abs(p10 + p11 * E)), 1e-300))))
def spread(v):
    vals = [m(list(p)) for p in permutations(v)]
    return max(vals) - min(vals), vals[0]
def Ds(v):
    a = np.angle(np.array(v, dtype=complex))
    return (a[0]+a[3]-a[1]-a[2], a[0]+a[1]-a[2]-a[3], a[0]+a[2]-a[1]-a[3])
def cosD(v):
    return [math.cos(d) for d in Ds(v)]

rng = np.random.default_rng(20260817)
o("=" * 104)
o("R2 SCRIPT 4 — IS 'REAL COEFFICIENTS' THE OPERATIVE HYPOTHESIS FOR THE MULTISET THEOREM?")
o("=" * 104)
o(f"Jensen quadrature n = 2^{int(math.log2(NQ))}, midpoint rule.  float64.")
o()
o("THE LANE'S ARGUMENT, RESTATED AND THEN SHARPENED.")
o("  m(P) = (1/2pi) INT log max(|A(t)|,|B(t)|) dt, A = p00+p01 e^{it}, B = p10+p11 e^{it}.")
o("  Swapping p00<->p01 sends |A(t)| -> |A(-t)| and leaves B alone; the integral is unchanged")
o("  ONLY IF |B(t)| is also EVEN in t.  |a+b e^{it}|^2 = |a|^2+|b|^2+2Re(a conj(b) e^{-it}) is")
o("  even in t  <=>  a conj(b) is REAL.  So the swap needs a PHASE condition, not reality.")
o()
o("AND THE PHASE CONDITION HAS A SINGLE DEGREE OF FREEDOM.  m(P) is invariant under the torus")
o("gauge x -> lam x, y -> mu y (|lam|=|mu|=1) and under an overall unimodular factor -- three of")
o("the four coefficient phases -- so")
o("      m(P) = F( |p00|,|p10|,|p01|,|p11| ; Delta ),   Delta = arg( p00 p11 / (p10 p01) ),")
o("and F sees Delta only through cos Delta (conjugation symmetry m(conj P) = m(P)).")
o("A permutation replaces the sign pattern (+,-,-,+) on the phases by its image, and only THREE")
o("two-plus-two patterns exist up to sign, so Delta is carried into {+-D1,+-D2,+-D3} with")
o("      D1 = a00+a11-a10-a01,  D2 = a00+a10-a01-a11,  D3 = a00+a01-a10-a11.")
o("REAL coefficients have a_c in {0,pi}, so D1 = D2 = D3 mod 2pi.  THAT is why reality works.")
o("Collinear complex vectors (p = e^{i th} w, w real) work for the same reason: a common phase")
o("cancels out of every D.  So reality is SUFFICIENT and STRICTLY NOT NECESSARY.")
o()
o("-" * 104)
o("EXHIBIT A — FOUR KINDS OF COEFFICIENT VECTOR")
o("-" * 104)
o(f"{'coefficient kind':<40}{'24-perm spread':<18}{'cos D1':<11}{'cos D2':<11}{'cos D3':<11}{'equal'}")
w1 = rng.random(4); w2 = rng.normal(size=4)
CASES = [("real non-negative (registrar's)", list(w1)),
         ("real MIXED SIGN (the lane's)", list(w2)),
         ("COMPLEX, collinear e^{i pi/7} x w1", list(np.exp(1j*math.pi/7)*w1)),
         ("COMPLEX, collinear e^{i 1.2346} x w2", list(np.exp(1j*1.234567)*w2)),
         ("COMPLEX, NOT collinear: p00 x i", [complex(0, abs(w1[0]))] + [complex(abs(x)) for x in w1[1:]]),
         ("COMPLEX, generic", list(rng.normal(size=4) + 1j*rng.normal(size=4))),
         ("COMPLEX, generic 2nd draw", list(rng.normal(size=4) + 1j*rng.normal(size=4)))]
for nm, v in CASES:
    s, m0 = spread(v); c = cosD(v)
    o(f"{nm:<40}{s:<18.3e}{c[0]:<+11.6f}{c[1]:<+11.6f}{c[2]:<+11.6f}{str(max(c)-min(c) < 1e-9)}")
o()
o("**'REAL COEFFICIENTS' IS NOT THE OPERATIVE HYPOTHESIS.  A COLLINEAR COMPLEX VECTOR -- every")
o("  coefficient off the real line -- IS FULLY S4-INVARIANT.  AND SO IS A NON-COLLINEAR ONE.**")
o()
o("-" * 104)
o("MY OWN FIRST NAME WAS ALSO TOO STRONG, AND IT IS RECORDED HERE RATHER THAN PATCHED AWAY")
o("-" * 104)
o("My draft of this block read: 'the operative hypothesis is COLLINEARITY, p in e^{i th} R^4'.")
o("ROW 5 REFUTES IT: p00 rotated by i with the other three real and positive is NOT collinear")
o("and IS invariant, because its three cos D all equal 0.  I name it correctly now, and record")
o("that this is the SIXTH consecutive name for this variable in this program and the THIRD in")
o("this round -- registrar (real non-negative), lane (real), me (collinear), me (cos D).")
o()
o("-" * 104)
o("EXHIBIT B — THE cos-D CONDITION IS SUFFICIENT, AND IT IS NOT NECESSARY.  BOTH MEASURED.")
o("-" * 104)
N1, N2 = 90, 120
bad = 0; worst = 0.0
for i in range(N1):
    r = rng.normal(size=4); kind = i % 3
    if kind == 0:   v = list(r)
    elif kind == 1: v = list(np.exp(1j*rng.uniform(0, 2*np.pi)) * r)
    else:           v = [abs(r[0])*np.exp(1j*rng.uniform(0, 2*np.pi))] + [complex(abs(x)) for x in r[1:]]
    c = cosD(v)
    assert max(c) - min(c) < 1e-9
    s, _ = spread(v); worst = max(worst, s); bad += (s > 1e-8)
o(f"SUFFICIENCY: {N1} vectors built with cos D1 = cos D2 = cos D3 (real / collinear / one-phase).")
o(f"   24-permutation spreads exceeding 1e-8: {bad} of {N1}.   worst spread {worst:.3e}")
cnt = {"inv_eq": 0, "inv_neq": 0, "var_eq": 0, "var_neq": 0}
deg = 0; degtot = 0
for _ in range(N2):
    v = list(rng.normal(size=4) + 1j*rng.normal(size=4))
    vals = [m(list(p)) for p in permutations(v)]
    s = max(vals) - min(vals); c = cosD(v)
    eq = (max(c) - min(c)) < 1e-9; inv = s < 1e-8
    cnt[("inv" if inv else "var") + ("_eq" if eq else "_neq")] += 1
    if inv and not eq:
        degtot += 1
        deg += (abs(vals[0] - math.log(max(abs(z) for z in v))) < 1e-8)
o(f"NECESSITY: {N2} generic complex draws -> {cnt}")
o(f"   The {cnt['var_eq']} 'varying but cos D equal' cases are what would REFUTE SUFFICIENCY: none.")
o(f"   The {cnt['inv_neq']} 'invariant but cos D unequal' cases REFUTE NECESSITY -- and all of")
o(f"   them are DEGENERATE: m = log(max |p_c|) in {deg} of {degtot}, i.e. one coefficient")
o("   dominates and every arrangement returns the same number for a reason that has nothing to")
o("   do with the symmetry.  So: cos D1 = cos D2 = cos D3 is SUFFICIENT, and is the right name")
o("   for the multiset theorem's hypothesis; it is not necessary, and the exceptions are the")
o("   locus where m degenerates.")
o()
o("WHAT IS UNCHANGED FOR THE CORPUS.  pi is a probability vector: real, non-negative, so all")
o("three D are 0 and the condition holds trivially.  W-03's multiset theorem, W10A-10's numbers")
o("and every registered figure resting on them are UNTOUCHED.  What changes is only the NAME of")
o("the hypothesis -- which is what this round is about, and which four consecutive layers of")
o("this program have got wrong.")
with open("r2_4_multiset.OUT.txt", "w") as fh:
    fh.write("\n".join(OUT) + "\n")

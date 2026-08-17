# W10-D REFUTER 1 -- LEG R3.  D-13 NAMES THE OPERATIVE VARIABLE ONE HYPOTHESIS TOO STRONG --
# THE EXACT DEFECT IT CONVICTS THE REGISTRAR OF, IN THE SAME ROW.
#
#   registrar's name : REAL NON-NEGATIVITY
#   lane D's name    : REALITY            ("the registrar's name is one hypothesis too strong")
#   MY FIRST NAME    : COLLINEARITY IN C  -- ALSO WRONG, AND RECORDED BELOW RATHER THAN PATCHED
#   the actual name  : HOW MANY OF THE FOUR **JENSEN-ADJACENT COEFFICIENT PAIRS** HAVE
#                      conj(p_i) p_j REAL.  4 or 2 -> 24.  1 -> 16.  0 -> 8.
#
# The four adjacent pairs are the two Jensen groupings' branches:
#     pairing in x : {p00,p01} and {p10,p11}      pairing in y : {p00,p10} and {p01,p11}
# |a + b e^{it}| = |b + a e^{it}| pointwise iff conj(a) b is REAL -- lane D's own stated
# mechanism.  It is a condition on PAIRS, and it is neither a condition on the field of the
# coefficients nor on their collective geometry.
#
# Consequence for D-13 as filed:
#   (i)  "REALITY" is sufficient, not operative: vectors with ZERO real entries give |G| = 24.
#   (ii) "COMPLEX coefficients drop |G| to exactly 8" is FALSE AS STATED: rotating exactly ONE
#        coefficient of a real vector by an arbitrary phase leaves |G| = 24, and there is a
#        third regime, |G| = 16, that D-13's dichotomy has no room for.
#
# float64.  Integrator: Jensen-in-x on the full circle, trapezoid, n = 2^20 (lane D used 2^18;
# every count below is re-checked at 2^18 and 2^21 and does not move).
import numpy as np
from itertools import permutations

rng = np.random.default_rng(20260816)
PAIRS = [(0, 2), (1, 3), (0, 1), (2, 3)]     # index order (00,10,01,11)
PAIRNAME = ["{00,01}", "{10,11}", "{00,10}", "{01,11}"]


def mahler(p, n=1 << 20):
    t = (np.arange(n)+0.5)*2*np.pi/n
    e = np.exp(1j*t)
    return np.log(np.maximum(np.abs(p[0] + p[2]*e), np.abs(p[1] + p[3]*e))).mean()


def Gset(p, n=1 << 20, tol=1e-9):
    b = mahler(p, n)
    return [q for q in permutations(range(4)) if abs(mahler(np.asarray(p)[list(q)], n) - b) < tol]


def n_adj_real(p, tol=1e-12):
    return sum(1 for i, j in PAIRS if abs(np.imag(np.conj(p[i])*p[j])) < tol)


def collinear(p, tol=1e-9):
    return max(abs(np.imag(np.conj(p[i])*p[j])) for i in range(4) for j in range(4)) < tol


def crossings(p):
    """how many of the 24 permutations give a genuinely CROSSING pair of Jensen branches.
       0 means m = log(max modulus) by domination for every arrangement -- a DEGENERATE probe
       on which |G| = 24 carries no information."""
    n = 0
    t = (np.arange(4096)+.5)*2*np.pi/4096
    e = np.exp(1j*t)
    for q in permutations(range(4)):
        v = np.asarray(p)[list(q)]
        A, B = np.abs(v[0]+v[2]*e), np.abs(v[1]+v[3]*e)
        if (A > B).any() and (B > A).any():
            n += 1
    return n


print("="*100)
print("== R3-0  MY OWN CONFOUND FIRST, RECORDED NOT PATCHED ==")
print("="*100)
print("  My first probe for 'is reality the operative variable?' was B0b's weight vector with one")
print("  entry rotated: (4/9, 2/9, 1/9, (2/9)e^{0.3i}).  It returned |G| = 24 and I nearly filed")
print("  that as the counterexample.  IT IS DEGENERATE: B0b's Jensen branches do not cross in ANY")
print("  of the 24 arrangements, so m = log(max modulus) by domination and |G| = 24 could not have")
print("  come out otherwise.  Every probe below prints its crossing count, and a probe with 0")
print("  crossings is not used as evidence.")
bad = np.array([4/9, 2/9, 1/9, (2/9)*np.exp(0.3j)])
print(f"    the degenerate probe: crossings = {crossings(bad)} of 24, |G| = {len(Gset(bad))}  -- DISCARDED")
good = np.array([0.55, 0.20, 0.40, 0.25])
print(f"    the probe base used instead, r = {good}: crossings = {crossings(good*np.array([1,-1,1,1]))} of 24")
print("  NOTE: lane D's own non-negative arms B0b and B4 are BOTH in the degenerate class")
print(f"    B0b (4/9,2/9,1/9,2/9): crossings = {crossings(np.array([4/9,2/9,1/9,2/9]))} of 24")
print(f"    B4  (1/6,1/6,1/6,3/6): crossings = {crossings(np.array([1/6,1/6,1/6,3/6]))} of 24")
print("  -- so two of lane D's seven arms in leg 2D could not have shown |G| < 24 either.  Its")
print("  non-degenerate arms (its 'generic' and its two negative-entry vectors) carry the row.")

print("\n"+"="*100)
print("== R3-A  ARM DIFF -- THE PROBES, AND THE ONE THING THAT MOVES ==")
print("="*100)
print("  Moving variable: WHICH ADJACENT COEFFICIENT PAIRS HAVE A REAL PRODUCT.  Moduli, the")
print("  integrator, the node count, the 24 permutations and the 1e-9 tolerance are fixed.")
r = np.array([0.55, 0.20, 0.40, 0.25])
CASES = [("C1 real, all + (non-negative)   ", r.astype(complex)),
         ("C2 real, one NEGATIVE entry     ", (r*np.array([1, -1, 1, 1])).astype(complex)),
         ("C3 e^{i0.9} * C2  (collinear)   ", np.exp(0.9j)*(r*np.array([1, -1, 1, 1])))]
for pos in range(4):
    v = (r*np.array([1, -1, 1, 1])).astype(complex)
    v[pos] *= np.exp(0.7j)
    CASES.append((f"C4 C2 with entry {pos} rotated 0.7 ", v))
v = (r*np.array([1, -1, 1, 1])).astype(complex)
v[0] *= np.exp(0.3j); v[1] *= np.exp(1.1j)
CASES.append(("C5 two entries rotated (adj-real=1)", v))
v = (r*np.array([1, -1, 1, 1])).astype(complex)
v[1] *= np.exp(0.4j); v[2] *= np.exp(0.9j); v[3] *= np.exp(1.3j)
CASES.append(("C6 torus shift (adj-real=0)     ", v))
CASES.append(("C7 lane D's complex generic arm ", np.array([0.4+0.1j, 0.2-0.3j, 0.3+0.25j, 0.1+0.4j])))
CASES.append(("C8 lane D's complex phases arm  ", np.array([1.0, np.exp(0.7j), np.exp(2.1j), np.exp(-1.3j)])))
for i in range(3):
    CASES.append((f"C9 random complex {i+1}           ", rng.normal(size=4)+1j*rng.normal(size=4)))
for lab, v in CASES:
    print(f"    {lab}  {np.array2string(np.asarray(v), precision=3)}")
seen = set(np.array2string(np.asarray(v), precision=14) for _, v in CASES)
print(f"    pairwise-distinct arms: {len(seen)} of {len(CASES)}")

print("\n"+"="*100)
print("== R3-B  |G| AGAINST THE ADJACENT-PAIR COUNT ==")
print("="*100)
print(f"  {'case':36s} {'#real entries':>13s} {'collinear':>10s} {'adj-real':>9s} {'|G|':>4s} "
      f"{'cross/24':>9s} {'adj pairs that are real':>26s}")
for lab, v in CASES:
    v = np.asarray(v)
    G = Gset(v)
    nr = int(sum(abs(np.imag(z)) < 1e-12 for z in v))
    which = ",".join(PAIRNAME[i] for i in range(4)
                     if abs(np.imag(np.conj(v[PAIRS[i][0]])*v[PAIRS[i][1]])) < 1e-12) or "none"
    print(f"  {lab:36s} {nr:13d} {str(collinear(v)):>10s} {n_adj_real(v):9d} {len(G):4d} "
          f"{crossings(v):9d} {which:>26s}")
print("\n  READ, AGAINST D-13 AS FILED:")
print("   * C3 has ZERO real entries and |G| = 24.  REALITY IS NOT THE OPERATIVE VARIABLE --")
print("     lane D's correction of the registrar is RIGHT, and its replacement name is wrong too.")
print("   * C4 (four rows) are COMPLEX, NOT collinear, non-degenerate, and give |G| = 24.  So my")
print("     own first replacement, COLLINEARITY, is ALSO one hypothesis too strong.  Recorded.")
print("   * C5 gives |G| = 16 -- a THIRD regime.  D-13's stated dichotomy 'real -> 24, complex ->")
print("     exactly 8' has no room for it, and 16 is not even a possible subgroup order in S4")
print("     (16 does not divide 24): G(p) = {sigma : m(sigma p) = m(p)} IS NOT A GROUP, because")
print("     permutation does not commute with the phase freedom.  The register, W-03, lane D and")
print("     the registrar all write '|G|' for it.")
print("   * C6-C9 have adj-real = 0 and give exactly 8 = the Newton-polygon D4.")
print("   THE NAME OF RECORD: |G| is decided by HOW MANY OF THE FOUR JENSEN-ADJACENT PAIRS")
print("   {00,01} {10,11} {00,10} {01,11} HAVE conj(p_i) p_j REAL.  Non-negativity gives 4;")
print("   reality gives 4; collinearity gives 4; but 2 also gives 24, 1 gives 16, 0 gives 8.")

print("\n  STABILITY OF THE THREE COUNTS UNDER THE INTEGRATOR (n = 2^18 / 2^20 / 2^21):")
for lab, v in [CASES[2], CASES[3], CASES[-5], CASES[-4]]:
    v = np.asarray(v)
    print(f"    {lab}  |G| = " + " / ".join(str(len(Gset(v, n=1 << e))) for e in (18, 20, 21)))

print("\n"+"="*100)
print("== R3-C  THE D4 THAT IS ALWAYS PRESENT, IDENTIFIED FROM THE SUBSTITUTIONS ==")
print("="*100)


def closure(gens):
    S = {(0, 1, 2, 3)}
    fr = [(0, 1, 2, 3)]
    while fr:
        x = fr.pop()
        for g in gens:
            y = tuple(x[g[i]] for i in range(4))
            if y not in S:
                S.add(y); fr.append(y)
    return S


D4 = closure([(1, 0, 3, 2), (2, 3, 0, 1), (0, 2, 1, 3)])
print("  x -> 1/x : (00 10)(01 11);  y -> 1/y : (00 01)(10 11);  x <-> y : (10 01).")
print(f"  |<those three>| = {len(D4)} -- present for EVERY coefficient vector, complex or not.")
print(f"  D4 + the four adjacent-pair transpositions closes to |S4| = "
      f"{len(closure([(1,0,3,2),(2,3,0,1),(0,2,1,3),(2,1,0,3),(0,3,2,1),(1,0,2,3),(0,1,3,2)]))}")
for lab, v in [CASES[-1], CASES[2]]:
    v = np.asarray(v)
    G = set(Gset(v))
    print(f"  {lab}  |G| = {len(G):2d}  equals D4 exactly? {G == D4}")

print("\n"+"="*100)
print("== R3-D  D-14's COUNT, RE-RUN OFF LANE D's RESONANT TEST POINT ==")
print("="*100)
print("  Lane D measured '2 of 24 permutations preserve |Z_k| pointwise' at (f,c) = (1.3, 2.0),")
print("  which leg R1 shows is EXACTLY RESONANT (20f - 13c = 0).  Re-run at generic connections,")
print("  generic weights, longer k-range:")
CLSE = ((0, 0), (1, 0), (0, 1), (1, 1))


def Zk(p, f, c, k):
    return sum(p[i]*np.exp(1j*k*(-CLSE[i][0]*f + CLSE[i][1]*c)) for i in range(4))


kk = np.arange(1, 41)
pts = [("lane D's (1.3, 2.0)  RESONANT ", (1.3, 2.0)),
       ("2pi(sqrt2-1), 2pi(sqrt3-1)    ", (2*np.pi*(np.sqrt(2)-1), 2*np.pi*(np.sqrt(3)-1))),
       ("(e, pi/e)                     ", (np.e, np.pi/np.e))]
for i in range(3):
    pts.append((f"rng connection {i+1}               ", tuple(rng.uniform(-np.pi, np.pi, 2))))
pgen = np.array([0.31, 0.17, 0.29, 0.23])
print(f"  generic weights p = {pgen}, k = 1..40")
for lab, (fv, cv) in pts:
    base = np.abs(Zk(pgen, fv, cv, kk))
    surv = [q for q in permutations(range(4))
            if np.abs(np.abs(Zk(pgen[list(q)], fv, cv, kk)) - base).max() < 1e-12]
    print(f"    {lab}  |Z_k|-preserving perms = {len(surv)} of 24  -> {surv}")
print("  D-14's COUNT SURVIVES -- identity and W-03's involution (3,2,1,0), 2 of 24, at every")
print("  connection tested including non-resonant ones.  Only the point it was measured at was")
print("  mislabelled 'generic'.")
print("\n  ITS NOVELTY CLAIM DOES NOT SURVIVE.  Both halves are already in the register:")
print("    REGISTER_V001.md:196-197  'lambda_B is a function of the MULTISET of the four class")
print("                               weights -- 24 of 24 permutations invariant'  [LAMBDA LEVEL]")
print("    REGISTER_V001.md:198-200  'multiplying Z_k by conj(u)^k conj(v)^k leaves |Z_k| fixed")
print("                               ... an exact symmetry at EVERY connection'   [|Z_k| LEVEL]")
print("  W-03 states the 24 at the lambda level and the involution at the |Z_k| level, in the same")
print("  row.  What is new in D-14 is only the SHARPNESS -- that no THIRD permutation survives.")
print("  D-14's self-flag ('If it is already somewhere in the corpus, this lane has committed the")
print("  corpus's own defect of not reading its own record') is the correct reading and it fires.")

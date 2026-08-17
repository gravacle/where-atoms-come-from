# LANE W11-R-CROSS  LEG F -- WHAT THE SURVIVING STATEMENT ACTUALLY SAYS, AND WHAT IT DOES NOT.
#
# F1.  The steelman's "CLOCK-FREE THEOREM OF RECORD" -- "Z(a,b) = <T_F^a s, T_C^b s> is a function
#      of pi alone EXACTLY on the lattice L_F.Z x L_C.Z ... No convention appears in that sentence"
#      -- HAS A CONVENTION IN IT: T.  The "exactly" half (falsity off the lattice) is a property of
#      COR-F's tick and of no other admissible root.  Under the uniform root the same object is
#      pi-only at EVERY (a,b), so the lattice does not characterise anything about the corpus.
# F2.  The correctly scoped survivor, checked: at joint loop closure the observable is pi-only for
#      EVERY root, because U^L = M.  It is a statement about WHEN you look, not about WHAT moves.
# F3.  And "carrier-independent" is false in the sense Reading A asserts, by the register's own
#      W-09: which classes the incidence occupies is carrier data and it changes the answer.
#      Reproduced here in my own code on K1 (three classes) and B0b (four classes).
# F4.  On B0b the corpus's own clock compares 4k edges against 3k edges.  "k circuits vs k circuits"
#      is a synchronisation choice, not a carrier-supplied one.
import numpy as np
import xlib as X

rng = np.random.default_rng(20260817)
CASES = [("K1 ", X.K1_LOOP_F, X.K1_LOOP_C, 5, np.array([1.0, 0.37, 0.91, 2 ** 0.5, 0.23, 1.77]),
          np.array([0.40, 0.15, 0.15, 0.15, 0.15])),
         ("B0b", X.B0B_LOOP_F, X.B0B_LOOP_C, 9, rng.uniform(0, 2 * np.pi, 18),
          np.array([.10, .12, .09, .14, .11, .11, .11, .11, .11]))]

print("== F1  THE 'CLOCK-FREE' THEOREM IS ROOT-SCOPED: THE LATTICE IS A PROPERTY OF T, NOT OF THE")
print("       OBSERVABLE.  Sweep (a,b) independently, once per transport ==")
NMAX = 13
print(f"  {'carrier':<6}{'transport':<22}{'pi-only (a,b) pairs, a,b < 13':>32}{'predicted lattice count':>26}")
for nm, lf, lc, NV, aa, w in CASES:
    LF, LC = len(lf), len(lc)
    w = w / w.sum()
    states = [np.sqrt(w) + 0j] + X.random_pi_identical(rng, lf, lc, NV, w, k=4)
    assert X.arms_differ(*states), "STATE ARMS BYTE-IDENTICAL -- leg void"
    ops = [("COR-F edge tick T", X.T_edge(lf, aa, NV), X.T_edge(lc, aa, NV)),
           ("uniform root D", X.D_uniform(lf, aa, NV), X.D_uniform(lc, aa, NV))]
    for rn, uF, uC in ops:
        hits = 0
        for A in range(1, NMAX):
            for B in range(1, NMAX):
                v = [abs(X.Z(uF, uC, s, A, B)) for s in states]
                hits += (max(v) - min(v) < 1e-12)
        print(f"  {nm:<6}{rn:<22}{hits:>32}{((NMAX-1)//LF)*((NMAX-1)//LC):>26}")
print("  -> for T the count is the lattice count; for D it is ALL 144 pairs.  The sentence")
print("     'invisible EXACTLY at joint loop closure' is true of COR-F's tick and false of the")
print("     other admissible root, so it names a convention after all -- the one the steelman's")
print("     own leg 5 rejects.  It cannot be offered as the convention-free statement of record.")

print("\n== F2  THE CORRECTLY SCOPED SURVIVOR: AT JOINT LOOP CLOSURE, EVERY ROOT IS THE CORPUS'S ==")
print("       (this is leg B1 restated as the positive claim, with the quantifier where it belongs)")
for nm, lf, lc, NV, aa, w in CASES:
    LF, LC = len(lf), len(lc)
    MF, MC = X.M_circuit(lf, aa, NV), X.M_circuit(lc, aa, NV)
    w = w / w.sum()
    s = np.sqrt(w) * np.exp(1j * rng.uniform(0, 2 * np.pi, NV))
    worst = 0.0
    for _ in range(150):
        uF = X.random_root(lf, aa, NV, rng, "generic")
        uC = X.random_root(lc, aa, NV, rng, "generic")
        for A in range(1, 4):
            for B in range(1, 4):
                worst = max(worst, abs(X.Z(uF, uC, s, LF * A, LC * B) - X.Z(MF, MC, s, A, B)))
    print(f"  {nm}  max over 150 roots and 9 lattice points | Z[U] - Z[M_gamma] | = {worst:.2e}")
print("  -> 'the incidence is invisible at joint loop closure' is TRUE, CONVENTION-FREE, and an")
print("     IDENTITY: U^L = M is the defining equation of every candidate.  What it is NOT is a")
print("     fact about transport, about the carrier, or about the physics.")

print("\n== F3  AND WHICH CLASSES THE INCIDENCE OCCUPIES IS CARRIER DATA (W-09), REPRODUCED HERE ==")
def fires(p, f, c):
    pts = []
    if p[0] > 1e-12: pts.append(1 + 0j)
    if p[1] > 1e-12: pts.append(np.exp(-1j * f))
    if p[2] > 1e-12: pts.append(np.exp(1j * c))
    if p[3] > 1e-12: pts.append(np.exp(1j * (c - f)))
    ws = np.array([x for x in (p[0], p[1], p[2], p[3]) if x > 1e-12])
    # 0 in convex hull of the occupied unit-modulus characters?
    ang = np.sort(np.angle(np.array(pts)))
    if len(ang) < 2: return False
    gaps = np.diff(np.concatenate([ang, [ang[0] + 2 * np.pi]]))
    return bool(np.max(gaps) <= np.pi + 1e-12)
N = 200000
fc = rng.uniform(-np.pi, np.pi, size=(N, 2))
for nm, lf, lc, NV, aa, w in CASES:
    p = X.pi_of(np.sqrt(w / w.sum()) + 0j, lf, lc, NV)
    occ = ",".join(t for t, x in zip(("00", "10", "01", "11"), p) if x > 1e-12)
    hit = 0
    flip = 0
    for i in range(N):
        f, c = fc[i]
        a1 = fires(p, f, c)
        hit += a1
        flip += (a1 != fires(p, -f, c))
    print(f"  {nm} occupied classes {{{occ}}}   firing region = {hit/N:.4f}   "
          f"verdict changes under f -> -f at {flip}/{N} points")
print("  -> exactly W-09's 1/4 (three classes, curvature-aware) and 1/2 (four classes,")
print("     curvature-BLIND).  The functional's FORM is carrier-free; its ANSWER is not.")
print("     Reading A's own slogan -- 'the physics does not depend on the complex' -- is false,")
print("     and the steelman concedes it.  What is left to decide is only the STATUS of the form.")

print("\n== F4  THE CORPUS'S CLOCK ON B0b COMPARES DIFFERENT NUMBERS OF EDGES ==")
lf, lc = X.B0B_LOOP_F, X.B0B_LOOP_C
print(f"  |gamma_F| = {len(lf)}, |gamma_C| = {len(lc)}.  At circuit k the two branches have taken")
print(f"  {len(lf)}k and {len(lc)}k edge steps.  Equal only at k = 0.  lcm = {np.lcm(len(lf),len(lc))}:")
print(f"  at n = 12 edge steps BOTH branches are at loop closure, but at circuit counts 3 and 4.")
print("  So 'k circuits against k circuits' is itself a synchronisation stipulation, no more")
print("  carrier-supplied than 'n edges against n edges'.  The steelman says this and it is right.")

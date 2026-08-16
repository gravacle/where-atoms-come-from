#!/usr/bin/env python3
"""
LANE R-GROUP REFUTER — part 3.

Fixes two defects I found in my own part 2 and then makes both demonstrations
EXACT (finite orbits, no estimator error at all).

DEFECT I FOUND IN MY OWN RUN (recorded, not silently fixed):
  part 2 chose SU(2) eigen-angles alpha=0.9, gamma=1.3. 0.9/1.3 = 9/13 EXACTLY,
  so 13*alpha - 9*gamma = 0 and the pair orbit lies on a SUBTORUS. My part-2 line
  "-> H = T^2" was FALSE — the same error S3 section 6(f) made and the erratum
  against W-02 corrects. Part 3 reports the subtorus case CORRECTLY (as a 1-D
  Mahler measure over H) and adds a genuinely generic pair.
"""
import numpy as np
from fractions import Fraction

NV = 5
A = np.array([1, 1, 1, 0, 0])
B = np.array([1, 0, 0, 1, 1])
p_S3 = np.array([0.4, 0.15, 0.15, 0.15, 0.15])
P0, Q0, R0 = 0.4, 0.3, 0.3          # class weights (1,1), (1,0), (0,1)


def banner(t):
    print("\n" + "=" * 78); print(t); print("=" * 78)


def su2(ax, ang):
    ax = np.asarray(ax, float); ax = ax / np.linalg.norm(ax)
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]])
    sz = np.array([[1, 0], [0, -1]], complex)
    n = ax[0] * sx + ax[1] * sy + ax[2] * sz
    return np.cos(ang / 2) * np.eye(2) - 1j * np.sin(ang / 2) * n


def Z_direct(WF, WC, S, ks):
    d = S.shape[1]
    out = np.empty(len(ks), complex)
    for i, k in enumerate(ks):
        FF = np.linalg.matrix_power(WF, int(k))
        CC = np.linalg.matrix_power(WC, int(k))
        tot = 0j
        for v in range(NV):
            L = FF @ S[v] if A[v] else S[v]
            R = CC @ S[v] if B[v] else S[v]
            tot += np.vdot(L, R)
        out[i] = tot
    return out


def char_coeffs(WF, WC, S, alpha, gamma):
    eF, PF = np.linalg.eig(WF)
    eC, QC = np.linalg.eig(WC)
    aF, aC = np.angle(eF), np.angle(eC)
    co = {}
    def add(k, val): co[k] = co.get(k, 0j) + val
    for v in range(NV):
        if A[v] and B[v]:
            x = PF.conj().T @ S[v]; y = QC.conj().T @ S[v]; Rm = PF.conj().T @ QC
            for j in range(2):
                for l in range(2):
                    add((-int(round(aF[j] / alpha)), int(round(aC[l] / gamma))),
                        np.conj(x[j]) * Rm[j, l] * y[l])
        elif A[v]:
            x = PF.conj().T @ S[v]
            for j in range(2): add((-int(round(aF[j] / alpha)), 0), abs(x[j]) ** 2)
        elif B[v]:
            y = QC.conj().T @ S[v]
            for l in range(2): add((0, int(round(aC[l] / gamma))), abs(y[l]) ** 2)
        else:
            add((0, 0), np.vdot(S[v], S[v]))
    return co


# =============================================================================
banner("A. SU(2) — lambda_B IS A MAHLER MEASURE OVER L-PERP, VERBATIM AS S4 SAYS")

rng = np.random.default_rng(2718281)            # SEED PUBLISHED
S = rng.normal(size=(NV, 2)) + 1j * rng.normal(size=(NV, 2))
for v in range(NV):
    S[v] *= np.sqrt(p_S3[v]) / np.linalg.norm(S[v])

for label, (alpha, gamma) in [("RESONANT  13a = 9g", (0.9, 1.3)),
                              ("GENERIC   a=0.9 g=sqrt(2)", (0.9, np.sqrt(2.0)))]:
    WF = su2([0, 0, 1], 2 * alpha)
    WC = su2([1, 0.4, 0.2], 2 * gamma)
    co = char_coeffs(WF, WC, S, alpha, gamma)
    ks = np.arange(0, 2001)
    dev = np.max(np.abs(Z_direct(WF, WC, S, ks) -
                        sum(c * np.exp(1j * ks * (m * alpha + n * gamma))
                            for (m, n), c in co.items())))
    # the relation lattice of the EIGEN-ANGLES
    rel = [(m, n) for m in range(-40, 41) for n in range(-40, 41)
           if (m, n) != (0, 0) and
           abs((m * alpha + n * gamma) / (2 * np.pi) -
               round((m * alpha + n * gamma) / (2 * np.pi))) < 1e-11]
    prim = min(rel, key=lambda t: abs(t[0]) + abs(t[1])) if rel else None
    print("\n  %s" % label)
    print("    ||[W_F,W_C]|| = %.6f   (non-abelian)" % np.linalg.norm(WF @ WC - WC @ WF))
    print("    character-sum validation vs genuine matrix powers, k<=2000: %.3e" % dev)
    print("    relation lattice of (alpha,gamma): primitive generator", prim,
          "-> rank L =", 0 if prim is None else 1)

    # direct schedule-B average
    kk = np.arange(1, 20_000_001)
    Zb = sum(c * np.exp(1j * kk * (m * alpha + n * gamma)) for (m, n), c in co.items())
    lam_dir = float(np.mean(np.log(np.abs(Zb))))
    print("    lambda_B direct, schedule B, N = 2e7        = %.9f" % lam_dir)

    if prim is None:
        # H = T^2: 2-variable Mahler measure, Jensen exact in y
        nx = 400_000
        t = (np.arange(nx) + 0.5) / nx * 2 * np.pi
        x = np.exp(1j * t)
        C = {n: np.zeros(nx, complex) for n in (-1, 0, 1)}
        for (m, n), c in co.items(): C[n] += c * x ** m
        out = np.empty(nx)
        for i in range(nx):
            r = np.roots([C[1][i], C[0][i], C[-1][i]])
            out[i] = np.log(abs(C[1][i])) + np.sum(np.log(np.maximum(1.0, np.abs(r))))
        lam_m = float(np.mean(out))
        print("    m(P) on H = T^2 (Jensen, 4e5 nodes)         = %.9f" % lam_m)
    else:
        # H = the circle { (e^{i n0 t}, e^{-i m0 t}) } annihilated by prim=(m0,n0)
        m0, n0 = prim
        g = np.gcd(abs(m0), abs(n0))
        e1, e2 = n0 // g, -m0 // g            # direction of L-perp
        nt = 4_000_000
        t = (np.arange(nt) + 0.5) / nt * 2 * np.pi
        Pv = sum(c * np.exp(1j * t * (m * e1 + n * e2)) for (m, n), c in co.items())
        lam_m = float(np.mean(np.log(np.abs(Pv))))
        print("    H = circle t -> (e^{i%dt}, e^{i%dt});  int_H log|P| dHaar = %.9f"
              % (e1, e2, lam_m))
    print("    |direct - Haar/Mahler| = %.2e" % abs(lam_dir - lam_m))

print("""
  READ THIS. At SU(2) — genuinely non-abelian, d = 2 — S4 section 2.2's structural
  theorem holds WORD FOR WORD:  lambda_B = int_{L-perp} log|P| d(Haar).  What changes
  is only the NEWTON POLYGON (unit square {0,1}^2 -> 3x3 square {-1,0,1}^2) and the
  fact that the coefficients are COMPLEX and depend on the state's DIRECTION in the
  fibre, not just on |s(v)|^2.  The character lattice and the Mahler measure DO NOT
  DIE AT d>1 UNDER A POWER SCHEDULE.""")


# =============================================================================
banner("B. THE SCHEDULE WALL, EXACTLY — d=1, U(1), UNIT CHARGE, RANK-2 LATTICE")

print("Both holonomies roots of unity -> S4's TIER 1: the orbit is FINITE and every")
print("average below is an EXACT finite sum. No estimator, no grid, no seed.\n")


def Zval(u, v, i, j):
    return P0 * u ** i * v ** j + Q0 * u ** i + R0 * v ** j


for (Nf, Nc) in [(3, 3), (4, 4), (3, 5), (6, 6)]:
    u = np.exp(-2j * np.pi / Nf)          # u = conj(W_F), W_F a primitive Nf-th root
    v = np.exp(+2j * np.pi / Nc)          # v = W_C
    L = int(np.lcm(Nf, Nc))
    # schedule B, k_n = n : orbit is the DIAGONAL {(u^k, v^k)}, period L
    lam_B = float(np.mean([np.log(abs(Zval(u, v, k, k))) for k in range(L)]))
    # word schedule: cell (i,j) runs any word with i F-letters and j C-letters.
    # At d=1 the operators commute (part 1: ||[M_F,M_C]|| = 0.000e+00), so the cell
    # value is Zval(u,v,i,j). Orbit is the FULL PRODUCT Z_Nf x Z_Nc.
    lam_W = float(np.mean([[np.log(abs(Zval(u, v, i, j))) for j in range(Nc)]
                           for i in range(Nf)]))
    print("  W_F = zeta_%d, W_C = zeta_%d :" % (Nf, Nc))
    print("     schedule B  (diagonal orbit, %2d points)  lambda = %.12f" % (L, lam_B))
    print("     WORD sched. (product orbit, %2d points)   lambda = %.12f"
          % (Nf * Nc, lam_W))
    print("     DIFFERENCE = %.12f" % (lam_B - lam_W))

print("""
  EXACT, closed-form, zero numerical uncertainty. Same carrier K1. Same U(1).
  Same unit charge. Same fibre dimension ONE. Same class weights (0.4,0.3,0.3).
  Same connection. Same relation lattice L. ONLY THE SCHEDULE VARIES — and
  lambda moves.  S4's "lambda_B is a function of L alone" is a theorem about
  the DIAGONAL orbit that the power schedule happens to trace; it is not a
  statement about the group, and it fails inside the box the claim calls intact.""")


# =============================================================================
banner("C. AND THE WORD SCHEDULE NEEDS NO NON-ABELIANNESS TO DO THIS")

print("At d=1 the two loop operators M_F, M_C are diagonal, hence commute, hence")
print("every word equals M_F^i M_C^j. Verified in part 1: ||[M_F,M_C]|| = 0.000e+00,")
print("and 3/3 random length-12 words matched their abelianisation to 4e-16.")
print("So the word-schedule closure H is ABELIAN, there IS a character lattice,")
print("there IS a Mahler measure — and lambda still moves. The claim's stated")
print("MECHANISM for the schedule route ('a word schedule would generate a")
print("genuinely non-abelian closure') is FALSE at d=1, and the route works anyway.")

print("\nWhat a word schedule actually changes, stated correctly:")
print("  power schedule: H = closure of the CYCLIC group <(u,v)>  -> MONOTHETIC,")
print("     H = L-perp, and lambda_B = m over L-perp.")
print("  word schedule:  H = closure of <u> x <v>  -> a PRODUCT, generally BIGGER")
print("     than L-perp. Abelian either way. The monothetic-ness, not the")
print("     abelian-ness, is what the schedule destroys.")

print("\nDONE.")

#!/usr/bin/env python3
"""
LANE R-GROUP REFUTER — part 2. Tightening the two loose numerics and adding
the continuity scan that locates the wall exactly.

Same conventions as r_group_refuter.py (carrier, d1, d2, Z_k, schedule B).
"""
import numpy as np

NV = 5
A = np.array([1, 1, 1, 0, 0])
B = np.array([1, 0, 0, 1, 1])
p_S3 = np.array([0.4, 0.15, 0.15, 0.15, 0.15])
f0, c0 = 2.0, 1.1


def banner(t):
    print("\n" + "=" * 78); print(t); print("=" * 78)


# ---------------------------------------------------------------- Mahler, exact-in-y
def mahler_via_jensen(coeffs, nx=200_000):
    """m(P), P(x,y) = sum_{m,n} c_{mn} x^m y^n with n in {-1,0,1}.
    Inner integral over y done EXACTLY by Jensen (roots of a quadratic)."""
    ms = sorted(set(m for (m, n) in coeffs))
    t = (np.arange(nx) + 0.5) / nx * 2 * np.pi
    x = np.exp(1j * t)
    C = {n: np.zeros(nx, dtype=complex) for n in (-1, 0, 1)}
    for (m, n), c in coeffs.items():
        C[n] += c * x ** m
    # P = C[-1]/y + C[0] + C[1] y ; multiply by y (|y|=1, no contribution):
    #   Q(y) = C[1] y^2 + C[0] y + C[-1]
    a2, a1, a0 = C[1], C[0], C[-1]
    out = np.zeros(nx)
    for i in range(nx):
        co = [a2[i], a1[i], a0[i]]
        while len(co) > 1 and abs(co[0]) < 1e-300:
            co = co[1:]
        if len(co) == 1:
            out[i] = np.log(abs(co[0])) if abs(co[0]) > 0 else -np.inf
            continue
        r = np.roots(co)
        out[i] = np.log(abs(co[0])) + np.sum(np.log(np.maximum(1.0, np.abs(r))))
    return float(np.mean(out)), ms


def su2(ax, ang):
    ax = np.asarray(ax, float); ax = ax / np.linalg.norm(ax)
    sx = np.array([[0, 1], [1, 0]], complex)
    sy = np.array([[0, -1j], [1j, 0]])
    sz = np.array([[1, 0], [0, -1]], complex)
    n = ax[0] * sx + ax[1] * sy + ax[2] * sz
    return np.cos(ang / 2) * np.eye(2) - 1j * np.sin(ang / 2) * n


def Z_direct(WF, WC, S, ks):
    """Z_k by GENUINE matrix powers — the reference implementation."""
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


def Z_from_chars(coeffs, alpha, gamma, ks):
    ks = np.asarray(ks, dtype=np.float64)
    Z = np.zeros(len(ks), complex)
    for (m, n), c in coeffs.items():
        Z += c * np.exp(1j * ks * (m * alpha + n * gamma))
    return Z


# =============================================================================
banner("4'. SU(2) — DIRECT LAMBDA_B vs MAHLER MEASURE, BOTH TIGHTENED")

alpha, gamma = 0.9, 1.3
WFs = su2([0, 0, 1], 2 * alpha)
WCs = su2([1, 0.4, 0.2], 2 * gamma)
rng4 = np.random.default_rng(2718281)          # SEED PUBLISHED (same as part 1)
Ssu = rng4.normal(size=(NV, 2)) + 1j * rng4.normal(size=(NV, 2))
for v in range(NV):
    Ssu[v] *= np.sqrt(p_S3[v]) / np.linalg.norm(Ssu[v])

eF, PF = np.linalg.eig(WFs)
eC, QC = np.linalg.eig(WCs)
aF, aC = np.angle(eF), np.angle(eC)
coeffs = {}
def add(k, val): coeffs[k] = coeffs.get(k, 0j) + val
for v in range(NV):
    if A[v] and B[v]:
        x = PF.conj().T @ Ssu[v]; y = QC.conj().T @ Ssu[v]; R = PF.conj().T @ QC
        for j in range(2):
            for l in range(2):
                add((-int(round(aF[j] / alpha)), int(round(aC[l] / gamma))),
                    np.conj(x[j]) * R[j, l] * y[l])
    elif A[v]:
        x = PF.conj().T @ Ssu[v]
        for j in range(2): add((-int(round(aF[j] / alpha)), 0), abs(x[j]) ** 2)
    elif B[v]:
        y = QC.conj().T @ Ssu[v]
        for l in range(2): add((0, int(round(aC[l] / gamma))), abs(y[l]) ** 2)
    else:
        add((0, 0), np.vdot(Ssu[v], Ssu[v]))

chk = np.arange(0, 2001)
Zref = Z_direct(WFs, WCs, Ssu, chk)
Zchr = Z_from_chars(coeffs, alpha, gamma, chk)
print("  VALIDATION: max |Z_k(genuine matrix powers) - character sum|, k<=2000 = %.3e"
      % np.max(np.abs(Zref - Zchr)))
for N in (60_000, 300_000, 2_000_000, 20_000_000):
    Z = Z_from_chars(coeffs, alpha, gamma, np.arange(1, N + 1))
    print("  lambda_B, schedule B, N=%10d  =  %.9f" % (N, np.mean(np.log(np.abs(Z)))))
m_su, ms = mahler_via_jensen(coeffs, nx=200_000)
print("  m(P) Jensen-exact-in-y, 200000 x-nodes  =  %.9f" % m_su)
print("  Newton polygon x-exponents:", ms, " y-exponents: {-1,0,1}")
print("  -> the 3x3 square {-1,0,1}^2, not the unit square {0,1}^2.")

# is (alpha,gamma,2pi) rationally independent?  (so H = T^2)
rel = [(m, n) for m in range(-30, 31) for n in range(-30, 31)
       if (m, n) != (0, 0) and abs((m * alpha + n * gamma) / (2 * np.pi) -
                                   round((m * alpha + n * gamma) / (2 * np.pi))) < 1e-12]
print("  relations (m,n) with m*alpha+n*gamma in 2piZ, |m|,|n|<=30 :", rel, "-> H = T^2")


# =============================================================================
banner("5b'. THE SCHEDULE WALL AT d=1 — PROVABLY EQUIDISTRIBUTING WORD SCHEDULE")

u0, v0_ = np.exp(-1j * f0), np.exp(1j * c0)

def Z_two_index(u, v, aa, bb, p):
    return (p[0] * u ** aa * v ** bb + (p[1] + p[2]) * u ** aa + (p[3] + p[4]) * v ** bb)

# schedule B (the corpus): a_n = b_n = n
n = np.arange(1, 4_000_001, dtype=np.int64)
lam_B = float(np.mean(np.log(np.abs(Z_two_index(u0, v0_, n, n, p_S3)))))

# WORD schedule: cell indexed by (i,j) over a 2000x2000 square; the cell runs any
# word with i letters F and j letters C. At d=1 every such word equals M_F^i M_C^j
# (proved in part 1: the operators commute), so the cell value is Z_{(i,j)}.
G = 2000
I, J = np.meshgrid(np.arange(1, G + 1, dtype=np.int64),
                   np.arange(1, G + 1, dtype=np.int64), indexing="ij")
lam_W = float(np.mean(np.log(np.abs(Z_two_index(u0, v0_, I.ravel(), J.ravel(), p_S3)))))

# reference: generic 2-torus Mahler measure, Jensen in y
tt = (np.arange(4_000_000) + 0.5) / 4_000_000 * 2 * np.pi
xx = np.exp(1j * tt)
m_gen = float(np.mean(np.log(np.maximum(np.abs(0.4 + 0.3 * xx), 0.3))))

print("  connection f=2.0 c=1.1 (RESONANT: -11f+20c = 0). SAME carrier, SAME class")
print("  weights (0.4,0.3,0.3), SAME group U(1), SAME unit charge, SAME d=1.")
print("  ONLY THE SCHEDULE IS VARIED.\n")
print("  schedule B  (k_n = n, powers)   lambda = %.9f   [record subtorus -0.767014993]" % lam_B)
print("  WORD schedule (i,j) 2000x2000   lambda = %.9f" % lam_W)
print("  m(0.4+0.3x+0.3y)  (full torus)  lambda = %.9f   [record -0.767507880]" % m_gen)
print("  |word - full-torus Mahler| = %.2e" % abs(lam_W - m_gen))
print("  |schedule B - word|        = %.6f   <-- THE SCHEDULE MOVES LAMBDA" % abs(lam_B - lam_W))

# convergence of the word schedule
for G2 in (200, 600, 1200, 2000):
    I2, J2 = np.meshgrid(np.arange(1, G2 + 1, dtype=np.int64),
                         np.arange(1, G2 + 1, dtype=np.int64), indexing="ij")
    lw = float(np.mean(np.log(np.abs(Z_two_index(u0, v0_, I2.ravel(), J2.ravel(), p_S3)))))
    print("     G=%4d  lambda_word = %.9f   (dev from full-torus Mahler %.2e)"
          % (G2, lw, abs(lw - m_gen)))


# =============================================================================
banner("3'. WHERE THE WALL ACTUALLY IS — CONTINUITY SCAN IN NON-CENTRALITY")

print("Family: W_F = diag(e^{if}, e^{i(f+delta)}), W_C = diag(e^{ic}, e^{i(c+delta)}).")
print("delta = 0 is CENTRAL (scalar). delta != 0 is abelian, d=2, NON-central.")
print("Statistic: spread of |Z_1| over 12 states with IDENTICAL class weights.\n")
rngd = np.random.default_rng(8675309)          # SEED PUBLISHED
states = []
for _ in range(12):
    T = rngd.normal(size=(NV, 2)) + 1j * rngd.normal(size=(NV, 2))
    for v in range(NV):
        T[v] *= np.sqrt(p_S3[v]) / np.linalg.norm(T[v])
    states.append(T)

def Z1(WF, WC, S):
    tot = 0j
    for v in range(NV):
        L = WF @ S[v] if A[v] else S[v]
        R = WC @ S[v] if B[v] else S[v]
        tot += np.vdot(L, R)
    return tot

print("   delta        ||[W_F,W_C]||    dist(W_F, centre)     spread |Z_1|")
for delta in (0.0, 1e-6, 1e-4, 1e-2, 0.1, 0.5, 1.0):
    WF = np.diag([np.exp(1j * f0), np.exp(1j * (f0 + delta))])
    WC = np.diag([np.exp(1j * c0), np.exp(1j * (c0 + delta))])
    vals = np.array([abs(Z1(WF, WC, T)) for T in states])
    dc = np.linalg.norm(WF - (np.trace(WF) / 2) * np.eye(2))
    print("  %8.1e     %.3e        %.6e        %.6e"
          % (delta, np.linalg.norm(WF @ WC - WC @ WF), dc, vals.max() - vals.min()))

print("\nSame scan, NON-ABELIAN family: W_F = su2(z, 2f), W_C = su2(axis, 2c) rotated")
print("towards the centre by theta (theta=0 -> scalar +-I only at special angles);")
print("here we interpolate the SU(2) pair towards a scalar by shrinking the angle:")
print("   ang          ||[W_F,W_C]||    dist(W_F, centre)     spread |Z_1|")
for ang in (0.0, 1e-6, 1e-4, 1e-2, 0.1, 0.5, 1.0):
    WF = su2([0, 0, 1], 2 * ang)
    WC = su2([1, 0.4, 0.2], 2 * ang)
    vals = np.array([abs(Z1(WF, WC, T)) for T in states])
    dc = np.linalg.norm(WF - (np.trace(WF) / 2) * np.eye(2))
    print("  %8.1e     %.3e        %.6e        %.6e"
          % (ang, np.linalg.norm(WF @ WC - WC @ WF), dc, vals.max() - vals.min()))

print("\nThe statistic degrades with DISTANCE FROM THE CENTRE, not with the commutator")
print("and not with d. At ||[W_F,W_C]|| = 0 exactly (abelian, d=2) the spread is")
print("LARGE whenever W_F is off-centre. The wall is CENTRALITY.")

print("\nDONE.")

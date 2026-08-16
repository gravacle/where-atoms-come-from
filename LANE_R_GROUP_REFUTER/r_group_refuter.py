#!/usr/bin/env python3
"""
LANE R-GROUP REFUTER — "Where Atoms Come From"
Target claim: box_wall — "the wall is d=1, not U(1); abelian is not a wall; the
untouched wall is the schedule".

Written from scratch. No code reused from any axis lane.
Everything printed. Seeds and grids published inline.

CONVENTIONS PUBLISHED HERE (S4's defect of record was not publishing these).

Vertex order:  (v0, v1, v2, v3, v4)  -> indices 0..4
Edge order:    (e1..e6) as S1 section 1:
   e1: v0->v1   e2: v1->v2   e3: v2->v0
   e4: v0->v3   e5: v3->v4   e6: v4->v0
Face order:    (F,) attached along e1.e2.e3
d1[edge, vertex] = +1 at target, -1 at source        (6 x 5)
d2[face, edge]   = +1 with the traversal orientation  (1 x 6)

Loop F = {v0,v1,v2}   Loop C = {v0,v3,v4}
class(v) = (a_v, b_v) = ([v in loop F], [v in loop C])
  v0=(1,1)  v1=(1,0)  v2=(1,0)  v3=(0,1)  v4=(0,1)

Transport operator, W-01 / S4 section 2 literal form:
  (M_gamma s)(v) = W(gamma) s(v)  if v on gamma, else s(v)
Cell value:
  Z_k = < M_F^k s , M_C^k s >  =  sum_v  s(v)^dag (W_F^dag)^{k a_v} (W_C)^{k b_v} s(v)
At d=1 with p_v=|s_v|^2, u=conj(W_F), v=W_C:  Z_k = sum_v u^{k a_v} v^{k b_v} p_v.
Schedule B: k_n = n.   lambda_B = lim (1/N) sum_{n<=N} log|Z_n|.
"""
import numpy as np

np.set_printoptions(precision=12, suppress=False, linewidth=140)

# ----------------------------------------------------------------------------
# 0. THE CARRIER, PUBLISHED
# ----------------------------------------------------------------------------
NV, NE, NF = 5, 6, 1
EDGES = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)]  # (source, target)
d1 = np.zeros((NE, NV), dtype=int)
for i, (s, t) in enumerate(EDGES):
    d1[i, s] -= 1
    d1[i, t] += 1
d2 = np.zeros((NF, NE), dtype=int)
d2[0, 0] = d2[0, 1] = d2[0, 2] = 1          # F attached along e1.e2.e3

LOOP_F = [0, 1, 2]
LOOP_C = [0, 3, 4]
A = np.array([1 if v in LOOP_F else 0 for v in range(NV)])   # a_v
B = np.array([1 if v in LOOP_C else 0 for v in range(NV)])   # b_v


def banner(t):
    print("\n" + "=" * 78)
    print(t)
    print("=" * 78)


banner("0. CARRIER K1 — INCIDENCE MATRICES PUBLISHED")
print("d1 (edge x vertex), rows e1..e6, cols v0..v4:")
print(d1)
print("d2 (face x edge), row F, cols e1..e6:")
print(d2)
print("d2 @ d1  (must be the zero 1x5 vector):", (d2 @ d1).ravel())
print("V,E,F =", NV, NE, NF, "  chi =", NV - NE + NF)
print("a_v =", A, "  b_v =", B)


# ----------------------------------------------------------------------------
# GENERAL TRANSPORT — ANY FIBRE DIMENSION, ANY UNITARY GROUP
# ----------------------------------------------------------------------------
def Zk_general(WF, WC, S, ks):
    """S: (5,d) complex, S[v] = s(v).  WF,WC: (d,d) unitary. Returns Z_k array.
    Uses the W-01 literal operator (M_gamma s)(v) = W s(v) on the loop."""
    d = S.shape[1]
    out = np.zeros(len(ks), dtype=complex)
    for i, k in enumerate(ks):
        FF = np.linalg.matrix_power(WF, int(k))
        CC = np.linalg.matrix_power(WC, int(k))
        tot = 0.0 + 0.0j
        for v in range(NV):
            L = FF if A[v] else np.eye(d, dtype=complex)
            R = CC if B[v] else np.eye(d, dtype=complex)
            tot += np.vdot(L @ S[v], R @ S[v])     # <M_F^k s, M_C^k s> at v
        out[i] = tot
    return out


def Zk_general_covariant(WF, WC, S, ks, GF, GC):
    """Claim's C-1 variant: based holonomy conjugated to each vertex,
    W^(v) = g_v W g_v^{-1}, with per-vertex gauge frames GF[v], GC[v]."""
    d = S.shape[1]
    out = np.zeros(len(ks), dtype=complex)
    for i, k in enumerate(ks):
        tot = 0.0 + 0.0j
        for v in range(NV):
            if A[v]:
                Wv = GF[v] @ WF @ GF[v].conj().T
                L = np.linalg.matrix_power(Wv, int(k)) @ S[v]
            else:
                L = S[v]
            if B[v]:
                Wv = GC[v] @ WC @ GC[v].conj().T
                R = np.linalg.matrix_power(Wv, int(k)) @ S[v]
            else:
                R = S[v]
            tot += np.vdot(L, R)
        out[i] = tot
    return out


def Zk_d1(u, v, p, ks):
    """Closed form at d=1 (S4 section 2)."""
    ks = np.asarray(ks)
    return np.array([sum(p[w] * u ** (k * A[w]) * v ** (k * B[w]) for w in range(NV))
                     for k in ks], dtype=complex)


def lambda_direct(Zfun, N):
    ks = np.arange(1, N + 1)
    Z = Zfun(ks)
    return float(np.mean(np.log(np.abs(Z))))


def mahler_2var(coeff_dict, ngrid=4000, seed=None):
    """m(P) for P(x,y)=sum c_{mn} x^m y^n on the 2-torus, tensor midpoint grid."""
    t = (np.arange(ngrid) + 0.5) / ngrid * 2 * np.pi
    X = np.exp(1j * t)[:, None]
    Y = np.exp(1j * t)[None, :]
    P = np.zeros((ngrid, ngrid), dtype=complex)
    for (m, n), c in coeff_dict.items():
        P += c * (X ** m) * (Y ** n)
    return float(np.mean(np.log(np.abs(P))))


def mahler_jensen_3term(p00, p10, p01, p11, n=2_000_000):
    """m(p00 + p10 x + p01 y + p11 xy) by Jensen in y:
       inner integral = log max(|p00+p10 x|, |p01+p11 x|)."""
    t = (np.arange(n) + 0.5) / n * 2 * np.pi
    x = np.exp(1j * t)
    return float(np.mean(np.log(np.maximum(np.abs(p00 + p10 * x),
                                           np.abs(p01 + p11 * x)))))


# ----------------------------------------------------------------------------
# 1. BASELINE — my code must reproduce the corpus at d=1
# ----------------------------------------------------------------------------
banner("1. BASELINE AT d=1 — does my implementation match the record?")

p_S3 = np.array([0.4, 0.15, 0.15, 0.15, 0.15])
f0, c0 = 2.0, 1.1                      # S3's headline connection (resonant)
u0, v0_ = np.exp(-1j * f0), np.exp(1j * c0)

Zclosed = Zk_d1(u0, v0_, p_S3, range(0, 201))
S1fib = np.zeros((NV, 1), dtype=complex)
rng = np.random.default_rng(20260816)          # SEED PUBLISHED
phases = rng.uniform(0, 2 * np.pi, NV)         # arbitrary vertex phases
S1fib[:, 0] = np.sqrt(p_S3) * np.exp(1j * phases)
Zmat = Zk_general(np.array([[np.exp(1j * f0)]]), np.array([[np.exp(1j * c0)]]),
                  S1fib, range(0, 201))
print("closed form vs direct matrix action, max |diff| over k<=200 : %.3e"
      % np.max(np.abs(Zclosed - Zmat)))

Zk400 = Zk_d1(u0, v0_, p_S3, range(1, 401))
print("min |Z_k|, k<=400  = %.6f at k=%d   (record: 0.024654 at k=42)"
      % (np.min(np.abs(Zk400)), 1 + int(np.argmin(np.abs(Zk400)))))
Zk4000 = Zk_d1(u0, v0_, p_S3, range(1, 4001))
print("sup |Z_k|, k<=4000 = %.6f at k=%d   (record: 0.999941 at k=377)"
      % (np.max(np.abs(Zk4000)), 1 + int(np.argmax(np.abs(Zk4000)))))
print("#{|Z_k|>0.99, k<=4000} = %d               (record: 37)"
      % int(np.sum(np.abs(Zk4000) > 0.99)))

m_generic = mahler_jensen_3term(0.4, 0.3, 0.3, 0.0)
print("m(0.4 + 0.3x + 0.3y) by Jensen, n=2e6      = %.9f   (record: -0.767507880)"
      % m_generic)
lamB_res = lambda_direct(lambda ks: Zk_d1(u0, v0_, p_S3, ks), 4_000_000)
print("lambda_B direct, f=2.0 c=1.1, N=4e6        = %.9f   (record: -0.767014993 SUBTORUS)"
      % lamB_res)
print("relation check -11f + 20c = %.3e  -> RESONANT, orbit on a subtorus" %
      (-11 * f0 + 20 * c0))
fg, cg = 1.0, np.sqrt(2.0)
lamB_gen = lambda_direct(lambda ks: Zk_d1(np.exp(-1j * fg), np.exp(1j * cg), p_S3, ks),
                         4_000_000)
print("lambda_B direct, f=1.0 c=sqrt2, N=4e6      = %.9f   (generic torus value)"
      % lamB_gen)


# ----------------------------------------------------------------------------
# 2. THE MISSING CONTROL — d = 2 WITH CENTRAL (SCALAR) HOLONOMY
#    This is the control that isolates FIBRE DIMENSION from NON-CENTRALITY.
#    The claim never ran it.
# ----------------------------------------------------------------------------
banner("2. THE CONTROL THAT ISOLATES d — FIBRE DIMENSION 2..7, SCALAR HOLONOMY")

print("For each d: W_F = e^{if} I_d, W_C = e^{ic} I_d, RANDOM state directions in C^d.")
print("If d were the wall, these must differ from the d=1 answer. They do not.\n")
rng2 = np.random.default_rng(11235)            # SEED PUBLISHED
for d in (1, 2, 3, 7):
    WF = np.exp(1j * f0) * np.eye(d, dtype=complex)
    WC = np.exp(1j * c0) * np.eye(d, dtype=complex)
    S = rng2.normal(size=(NV, d)) + 1j * rng2.normal(size=(NV, d))
    # normalise each fibre to the prescribed class weight p_S3
    for v in range(NV):
        S[v] *= np.sqrt(p_S3[v]) / np.linalg.norm(S[v])
    Zg = Zk_general(WF, WC, S, range(0, 201))
    print("  d=%d : max |Z_k(d) - Z_k(d=1 closed form)| over k<=200 = %.3e"
          % (d, np.max(np.abs(Zg - Zclosed))))

print("\nSix states with IDENTICAL class weights but different directions, d=2 SCALAR:")
spread = []
for trial in range(6):
    d = 2
    S = rng2.normal(size=(NV, d)) + 1j * rng2.normal(size=(NV, d))
    for v in range(NV):
        S[v] *= np.sqrt(p_S3[v]) / np.linalg.norm(S[v])
    spread.append(abs(Zk_general(np.exp(1j * f0) * np.eye(2),
                                 np.exp(1j * c0) * np.eye(2), S, [1])[0]))
spread = np.array(spread)
print("  |Z_1| over six states = ", spread)
print("  SPREAD = %.6e     (W-03 records 0.4247 under SU(2))" % (spread.max() - spread.min()))
print("  -> at d=2 with CENTRAL holonomy the cancellation of section 2 HOLDS EXACTLY.")


# ----------------------------------------------------------------------------
# 3. d = 2 ABELIAN, NON-CENTRAL — is this a NEW phenomenon, or charge in disguise?
# ----------------------------------------------------------------------------
banner("3. d=2 ABELIAN NON-CENTRAL — EXACT d=1 SURROGATE")

q1, q2 = 1, 2                                    # the two U(1) weights (charges)
WF2 = np.diag([np.exp(1j * q1 * f0), np.exp(1j * q2 * f0)])
WC2 = np.diag([np.exp(1j * q1 * c0), np.exp(1j * q2 * c0)])
print("W_F = diag(e^{i q f}), W_C = diag(e^{i q c}), q=(%d,%d): commutator norm = %.3e"
      % (q1, q2, np.linalg.norm(WF2 @ WC2 - WC2 @ WF2)))

rng3 = np.random.default_rng(31415)              # SEED PUBLISHED
S = rng3.normal(size=(NV, 2)) + 1j * rng3.normal(size=(NV, 2))
for v in range(NV):
    S[v] *= np.sqrt(p_S3[v]) / np.linalg.norm(S[v])
Zab = Zk_general(WF2, WC2, S, range(0, 401))

# the surrogate: TWO DISJOINT d=1 copies of K1, charges q1 and q2, weights = |S[v,j]|^2
w1 = np.abs(S[:, 0]) ** 2
w2 = np.abs(S[:, 1]) ** 2
Zsur = (Zk_d1(np.exp(-1j * q1 * f0), np.exp(1j * q1 * c0), w1, range(0, 401))
        + Zk_d1(np.exp(-1j * q2 * f0), np.exp(1j * q2 * c0), w2, range(0, 401)))
print("max |Z_k(d=2 abelian) - Z_k(d=1 two-copy charged surrogate)|, k<=400 = %.3e"
      % np.max(np.abs(Zab - Zsur)))
print("-> the d=2 abelian system IS a d=1 system on the disjoint carrier K1 (+) K1")
print("   with charges (%d,%d). Nothing about dimension two is doing any work." % (q1, q2))

print("\nSix states, identical class weights, d=2 ABELIAN NON-CENTRAL:")
sp = []
for trial in range(6):
    T = rng3.normal(size=(NV, 2)) + 1j * rng3.normal(size=(NV, 2))
    for v in range(NV):
        T[v] *= np.sqrt(p_S3[v]) / np.linalg.norm(T[v])
    sp.append(abs(Zk_general(WF2, WC2, T, [1])[0]))
sp = np.array(sp)
print("  |Z_1| =", sp, "\n  SPREAD = %.6f  -> class-weight statistic FAILS, abelian, d=2" %
      (sp.max() - sp.min()))
print("  and it fails through the SAME object as W-03's charge run: the character set.")

# is the character set still the four corners of the unit square?
print("\n  characters present, d=1 charge 1 : (m,n) in {(0,0),(1,0),(0,1),(1,1)}  4 corners")
print("  characters present, d=2 q=(1,2)  : (m,n) in {(0,0),(1,0),(2,0),(0,1),(0,2),"
      "(1,1),(2,2)}  -> NOT a unit square (S4-1's hypothesis gone), exactly as the "
      "charge run of record.")


# ----------------------------------------------------------------------------
# 4. d = 2 NON-ABELIAN SU(2) — DOES THE CHARACTER LATTICE / MAHLER MEASURE DIE?
# ----------------------------------------------------------------------------
banner("4. SU(2), POWER SCHEDULE — THE CHARACTER LATTICE AND MAHLER MEASURE SURVIVE")


def su2(ax, ang):
    ax = np.asarray(ax, dtype=float)
    ax = ax / np.linalg.norm(ax)
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]])
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    n = ax[0] * sx + ax[1] * sy + ax[2] * sz
    return np.cos(ang / 2) * np.eye(2) - 1j * np.sin(ang / 2) * n


alpha, gamma = 0.9, 1.3                       # HALF-angles enter as +-alpha, +-gamma
WFs = su2([0, 0, 1], 2 * alpha)               # eigenvalues e^{+-i alpha}
WCs = su2([1, 0.4, 0.2], 2 * gamma)           # eigenvalues e^{+-i gamma}
print("commutator norm ||[W_F,W_C]|| = %.6f  (genuinely non-abelian)"
      % np.linalg.norm(WFs @ WCs - WCs @ WFs))

rng4 = np.random.default_rng(2718281)         # SEED PUBLISHED
Ssu = rng4.normal(size=(NV, 2)) + 1j * rng4.normal(size=(NV, 2))
for v in range(NV):
    Ssu[v] *= np.sqrt(p_S3[v]) / np.linalg.norm(Ssu[v])

# --- reproduce the record's breakage: identical class weights, different |Z_1|
sp = []
for trial in range(6):
    T = rng4.normal(size=(NV, 2)) + 1j * rng4.normal(size=(NV, 2))
    for v in range(NV):
        T[v] *= np.sqrt(p_S3[v]) / np.linalg.norm(T[v])
    sp.append(abs(Zk_general(WFs, WCs, T, [1])[0]))
sp = np.array(sp)
print("six states, identical class weights, SU(2): |Z_1| spread = %.4f  "
      "(record phenomenon reproduced)" % (sp.max() - sp.min()))

# --- NOW: the character decomposition the claim says does not exist
eF, PF = np.linalg.eig(WFs)
eC, QC = np.linalg.eig(WCs)
aF = np.angle(eF)          # (+-alpha)
aC = np.angle(eC)          # (+-gamma)
print("eigen-angles of W_F:", aF, "   of W_C:", aC)

coeffs = {}


def add(key, val):
    coeffs[key] = coeffs.get(key, 0j) + val


# class (0,0): none on K1 (all five vertices lie on a loop) -> constant term absent
for v in range(NV):
    if A[v] and B[v]:                      # v0
        x = PF.conj().T @ Ssu[v]
        y = QC.conj().T @ Ssu[v]
        R = PF.conj().T @ QC
        for j in range(2):
            for l in range(2):
                add((-int(round(aF[j] / alpha)), int(round(aC[l] / gamma))),
                    np.conj(x[j]) * R[j, l] * y[l])
    elif A[v]:
        x = PF.conj().T @ Ssu[v]
        for j in range(2):
            add((-int(round(aF[j] / alpha)), 0), np.abs(x[j]) ** 2)
    elif B[v]:
        y = QC.conj().T @ Ssu[v]
        for l in range(2):
            add((0, int(round(aC[l] / gamma))), np.abs(y[l]) ** 2)
    else:
        add((0, 0), np.vdot(Ssu[v], Ssu[v]))

print("\nNewton polygon of the SU(2) cell value (exponents (m,n) of e^{i k (m*alpha + n*gamma)}):")
for kk in sorted(coeffs):
    print("   (%2d,%2d)  c = %+.9f %+.9fi" % (kk[0], kk[1], coeffs[kk].real, coeffs[kk].imag))

ks = np.arange(0, 301)
Zsu = Zk_general(WFs, WCs, Ssu, ks)
Zchar = np.zeros(len(ks), dtype=complex)
for (m, n), c in coeffs.items():
    Zchar += c * np.exp(1j * ks * (m * alpha + n * gamma))
print("\nmax |Z_k(SU(2), direct) - character sum| over k<=300 = %.3e" %
      np.max(np.abs(Zsu - Zchar)))

lam_su_direct = lambda_direct(lambda kk: Zk_general(WFs, WCs, Ssu, kk), 60000)
# Mahler: substitute x=e^{i alpha}, y=e^{i gamma}; generic (alpha,gamma) -> H = T^2
lam_su_mahler = mahler_2var(coeffs, ngrid=3000)
print("lambda_B direct (schedule B, N=60000)      = %.9f" % lam_su_direct)
print("lambda_B as Mahler measure m(P) on T^2     = %.9f   (3000^2 grid)" % lam_su_mahler)
print("|difference| = %.3e" % abs(lam_su_direct - lam_su_mahler))
print("-> at SU(2), under a POWER schedule, lambda_B IS a Mahler measure over the")
print("   closure H of the pair orbit. The Newton polygon grows from the UNIT SQUARE")
print("   {0,1}^2 to the 3x3 square {-1,0,1}^2 and the coefficients become COMPLEX")
print("   and state-dependent. The MECHANISM is untouched; only the STATISTIC dies.")


# ----------------------------------------------------------------------------
# 5. THE SCHEDULE — DOES IT REQUIRE NON-ABELIAN? NO. IT BITES AT d=1.
# ----------------------------------------------------------------------------
banner("5. WORD SCHEDULES AT d=1 — words abelianize, and lambda MOVES ANYWAY")

# 5a. at d=1 the two loop operators commute: every word = M_F^a M_C^b
MF1 = np.diag([np.exp(1j * f0) if A[v] else 1.0 for v in range(NV)])
MC1 = np.diag([np.exp(1j * c0) if B[v] else 1.0 for v in range(NV)])
print("d=1: ||[M_F, M_C]|| = %.3e   -> words abelianize EXACTLY" %
      np.linalg.norm(MF1 @ MC1 - MC1 @ MF1))
MFs = np.eye(NV * 2, dtype=complex)   # placeholder, replaced below for SU(2)

# SU(2) analogue of the same two operators, to contrast
def big(W, mask, d):
    M = np.zeros((NV * d, NV * d), dtype=complex)
    for v in range(NV):
        M[v * d:(v + 1) * d, v * d:(v + 1) * d] = W if mask[v] else np.eye(d)
    return M


MFb, MCb = big(WFs, A, 2), big(WCs, B, 2)
print("SU(2): ||[M_F, M_C]|| = %.6f   -> words do NOT abelianize" %
      np.linalg.norm(MFb @ MCb - MCb @ MFb))

rngw = np.random.default_rng(161803)          # SEED PUBLISHED
for trial in range(3):
    w = rngw.integers(0, 2, size=12)          # a word: 0 = F letter, 1 = C letter
    Md1 = np.eye(NV, dtype=complex)
    Md2 = np.eye(NV * 2, dtype=complex)
    for letter in w:
        Md1 = (MF1 if letter == 0 else MC1) @ Md1
        Md2 = (MFb if letter == 0 else MCb) @ Md2
    a, b = int(np.sum(w == 0)), int(np.sum(w == 1))
    ab1 = np.linalg.matrix_power(MF1, a) @ np.linalg.matrix_power(MC1, b)
    ab2 = np.linalg.matrix_power(MFb, a) @ np.linalg.matrix_power(MCb, b)
    print("  word %s -> (a,b)=(%d,%d) : d=1 dev %.3e | SU(2) dev %.3e"
          % ("".join("FC"[x] for x in w), a, b,
             np.linalg.norm(Md1 - ab1), np.linalg.norm(Md2 - ab2)))

# 5b. THE TEST. Same carrier, same connection (RESONANT), same class weights,
#     same group U(1), same charge 1, same d=1. ONLY THE SCHEDULE CHANGES.
banner("5b. ONE OBJECT VARIED: THE SCHEDULE. d=1, U(1), UNIT CHARGE, RESONANT CONNECTION")


def lambda_two_index(u, v, p, aa, bb):
    Z = np.array([sum(p[w] * u ** (aa[i] * A[w]) * v ** (bb[i] * B[w])
                      for w in range(NV)) for i in range(len(aa))], dtype=complex)
    return float(np.mean(np.log(np.abs(Z)))), Z


Nw = 400_000
# canonical clock, as the corpus: a_n = b_n = n
aa = np.arange(1, Nw + 1)
lam_pow, _ = lambda_two_index(u0, v0_, p_S3, aa, aa)
# the claim's own proposed WORD schedule, at d=1: cell n runs a random word of
# length n; its abelianization is (a_n, b_n) with a_n + b_n = n.
rngs = np.random.default_rng(577215664)       # SEED PUBLISHED
aw = np.zeros(Nw, dtype=np.int64)
bw = np.zeros(Nw, dtype=np.int64)
for n in range(1, Nw + 1):
    a = int(rngs.binomial(n, 0.5))            # #F letters in a length-n fair word
    aw[n - 1], bw[n - 1] = a, n - a
lam_word, _ = lambda_two_index(u0, v0_, p_S3, aw, bw)

print("connection f=2.0, c=1.1  (-11f+20c = 0, RESONANT: orbit on a SUBTORUS)")
print("  schedule B, k_n = n            lambda = %.9f   (record subtorus  -0.767014993)"
      % lam_pow)
print("  WORD schedule, |w_n| = n       lambda = %.9f   (generic torus    -0.767507880)"
      % lam_word)
print("  m(0.4+0.3x+0.3y) Jensen        lambda = %.9f" % m_generic)
print("  |word - generic Mahler| = %.2e     |power - word| = %.6f"
      % (abs(lam_word - m_generic), abs(lam_pow - lam_word)))
print("\n  The word schedule at d=1 has an ABELIAN closure (the operators commute),")
print("  there IS a character lattice, there IS a Mahler measure — and lambda MOVES,")
print("  off the subtorus value onto the full-torus value. The schedule wall does NOT")
print("  need non-abelianness. It bites inside the box the claim says is intact.")

# 5c. and it is not a resonance-only effect: independent-count schedules change H
banner("5c. WHY: the pair-orbit stops being MONOTHETIC")
print("power schedule:  H = closure{(u^n, v^n)}  = L^perp,  L = {(m,n): u^m v^n = 1}")
print("word schedule:   H = closure{(u^a_n, v^b_n)} with a_n, b_n asymptotically")
print("                 independent  ->  H = closure<u> x closure<v> = T^2 here.")
print("Rank of the relation lattice L at f=2.0,c=1.1:")
best = None
for m in range(-40, 41):
    for n in range(-40, 41):
        if (m, n) == (0, 0):
            continue
        r = (-m * f0 + n * c0) / (2 * np.pi)
        if abs(r - round(r)) < 1e-12:
            if best is None or abs(m) + abs(n) < abs(best[0]) + abs(best[1]):
                best = (m, n)
print("  smallest relation (m,n) with W_F^{-m} W_C^{n} = 1 :", best)
print("  -> rank L = 1 under the power schedule (a circle), rank 0 under the word")
print("     schedule (all of T^2). SAME connection. SAME L. DIFFERENT lambda.")
print("  So 'lambda_B is a function of L alone' is a statement about the SCHEDULE,")
print("  not about the group, and it already fails at d=1, U(1), unit charge.")


# ----------------------------------------------------------------------------
# 6. THE GAUGE COUNT — is d>1 'inside S1's own definition'?
# ----------------------------------------------------------------------------
banner("6. S1 SECTION 4's SATURATED INVARIANT COUNT AT d>1")

print("S1 section 4: 6 edge params - 4 effective gauge params = 2 invariants,")
print("and 'exactly two invariants exist, so the count is saturated'.")
print("At fibre dimension d the same count reads:")
for d in (1, 2, 3):
    edge = NE * d * d                # dim U(d)^6
    gauge = NV * d * d - 1           # dim U(d)^5 minus the global stabiliser
    print("   d=%d :  edges %3d  -  effective gauge %3d  =  %3d gauge-invariant reals"
          % (d, edge, gauge, edge - gauge))
print("Numerical check at d=2 (rank of the gauge-orbit tangent space at a random point):")
rngg = np.random.default_rng(999983)          # SEED PUBLISHED


def rand_u(d, rg):
    Z = rg.normal(size=(d, d)) + 1j * rg.normal(size=(d, d))
    Q, R = np.linalg.qr(Z)
    return Q @ np.diag(np.diag(R) / np.abs(np.diag(R)))


for d in (1, 2):
    Us = [rand_u(d, rngg) for _ in range(NE)]
    # tangent to the gauge orbit: dU_e = i(X_target U_e - U_e X_source), X_v hermitian
    basis = []
    herm = []
    for i in range(d):
        E = np.zeros((d, d), dtype=complex); E[i, i] = 1; herm.append(E)
    for i in range(d):
        for j in range(i + 1, d):
            E = np.zeros((d, d), dtype=complex); E[i, j] = 1; E[j, i] = 1; herm.append(E)
            E = np.zeros((d, d), dtype=complex); E[i, j] = -1j; E[j, i] = 1j; herm.append(E)
    for v in range(NV):
        for H in herm:
            vec = []
            for e, (s, t) in enumerate(EDGES):
                dU = np.zeros((d, d), dtype=complex)
                if t == v:
                    dU += 1j * H @ Us[e]
                if s == v:
                    dU -= 1j * Us[e] @ H
                # express in the tangent space of U(d) at U_e:  dU = i K U_e, K herm
                K = -1j * dU @ Us[e].conj().T
                vec.extend([K[a, b].real for a in range(d) for b in range(d)])
                vec.extend([K[a, b].imag for a in range(d) for b in range(d)])
            basis.append(vec)
    M = np.array(basis)
    r = np.linalg.matrix_rank(M, tol=1e-8)
    print("   d=%d : rank of gauge action = %d,  dim U(d)^6 = %d,  invariants = %d"
          % (d, r, NE * d * d, NE * d * d - r))
print("At d=2 the invariant count is 5, not 2, and W_F, W_C do not exhaust it")
print("(tr(W_F W_C) and the eigenbasis angle are extra). S1 section 4's SATURATION —")
print("its stated reason for believing nothing is hidden — is FALSE at d>1.")


# ----------------------------------------------------------------------------
# 7. ROBUSTNESS: does any of this depend on the claim's C-1 operator choice?
# ----------------------------------------------------------------------------
banner("7. ROBUSTNESS UNDER THE CLAIM'S C-1 (gauge-covariant based holonomy)")

GF = [rand_u(2, rngg) for _ in range(NV)]
GC = [rand_u(2, rngg) for _ in range(NV)]
Zc_scalar = Zk_general_covariant(np.exp(1j * f0) * np.eye(2),
                                 np.exp(1j * c0) * np.eye(2), Ssu, range(0, 201), GF, GC)
print("d=2 SCALAR under C-1: max |Z_k - d=1 closed form| = %.3e"
      % np.max(np.abs(Zc_scalar -
                      Zk_d1(u0, v0_, np.array([np.vdot(Ssu[v], Ssu[v]).real
                                               for v in range(NV)]), range(0, 201)))))
print("-> the d=2-scalar counterexample is INDEPENDENT of C-1: it holds under both")
print("   operator conventions, because a scalar commutes with every gauge frame.")

print("\nDONE.")

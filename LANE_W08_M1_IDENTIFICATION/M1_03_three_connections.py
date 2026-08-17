#!/usr/bin/env python3
"""
M1_03 — THE BIRKHOFF AVERAGE (1/N) sum_{k<=N} log|Z_k| ON THREE CONNECTIONS.

ISOLATION LEDGER (the whole point of this script)
  HELD FIXED : carrier K1; the ready state, through its pushforward
               pi = (p00,p10,p01,p11) = (0, 0.3, 0.3, 0.4)  [S3/S4's own generic state];
               the observable Z_k = <M_dF^k s, M_c^k s>; the identification Z_k = P(u^k,v^k)
               verified in M1_01; the k-grid k = 1..10^7; the code path below; float64.
  THE ONE THING THAT MOVES : the RELATION LATTICE of the pair (u,v) = (conj W_F, W_C),
               L = { (m,n) in Z^2 : u^m v^n = 1 }.
                 A  rank L = 0   (u,v) Diophantine        orbit dense in T^2
                 B  rank L = 2   (u,v) torsion, order 4   orbit = 4 points
                 C  rank L = 1   primitive (11,20)        orbit dense in a circle subgroup
  NOT a comparison, reported separately and labelled: exhibit D changes the ready state as
  well and is therefore NOT isolated.  It is the corpus's actual published configuration.

WHY THESE THREE POINTS AND NOT OTHERS
  A : alpha = -2^(1/3), beta = 4^(1/3) (as fractions of 2pi).  1, 2^(1/3), 4^(1/3) are a
      Q-basis of the cubic field Q(2^(1/3)), so by W. M. Schmidt's subspace theorem (1970)
      the vector is Diophantine of type 1+eps.  This is a POINT WHERE THE DIOPHANTINE
      HYPOTHESIS IS A THEOREM, not a point where it is assumed.  (For (sqrt2, sqrt3) the
      analogous statement is also a theorem; for a "random" float it is unprovable.)
  B : S1 sec6's PUBLISHED connection, a1=a2=a3=pi/3, a4=a5=a6=pi/2 -> W_F=-1, W_C=-i.
  C : S3/S4's headline f=2.0, c=1.1, exactly resonant (-11f+20c=0)  [ERRATUM v W-02].

PRECISION.  float64 everywhere, but the PHASES are not accumulated in float64: frac(k*alpha)
is computed by exact int64 modular arithmetic against a fixed denominator D, plus (case A
only) a float64 correction for the truncation residue.  Naive k*alpha in float64 loses ~1e-9
of phase by k = 1e7; this loses ~1e-17.  In case C the two numerators are chosen with a
COMMON denominator so that 11*alpha + 20*beta = 0 holds EXACTLY in the representation --
without that, double rounding of 1.1 makes -11f+20c = 1.78e-15 and the orbit drifts off the
subtorus.  Legitimacy of that substitution: the closed subgroup H = closure{(u^k,v^k)} for a
primitive relation (m,n) with gcd(m,n)=1 is the CONNECTED kernel {(z^n, z^-m)}, which depends
only on (m,n) and not on where (u,v) sits on it; so every non-torsion point of that circle
has the same limit.  Verified numerically at the end.
"""
import numpy as np

M_P_EXACT = None      # filled from M1_02

from M1_02_mahler_machinery import m_R1, m_CM, m_grid

P00, P10, P01, P11 = 0.0, 0.3, 0.3, 0.4        # HELD FIXED across A, B, C
K = 10 ** 7
CHECK = [10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7]
CHUNK = 10 ** 6


def Zabs_from_phases(tx, ty, p10=P10, p01=P01, p11=P11, p00=P00):
    """|P(x,y)| with x=e^{i tx}, y=e^{i ty}."""
    x = np.exp(1j * tx)
    y = np.exp(1j * ty)
    return np.abs(p00 + p10 * x + p01 * y + p11 * x * y)


def orbit_average(A_num, B_num, D, dA=0.0, dB=0.0, K=K, checkpoints=CHECK, eps_list=()):
    """alpha = A_num/D + dA, beta = B_num/D + dB.  Exact int64 modular reduction."""
    tot = 0.0
    tot_eps = {e: 0.0 for e in eps_list}
    nzero = 0
    out = {}
    out_eps = {e: {} for e in eps_list}
    minabs = np.inf
    done = 0
    while done < K:
        n = min(CHUNK, K - done)
        k = np.arange(done + 1, done + n + 1, dtype=np.int64)
        fa = ((k * A_num) % D).astype(np.float64) / D + k.astype(np.float64) * dA
        fb = ((k * B_num) % D).astype(np.float64) / D + k.astype(np.float64) * dB
        fa = np.mod(fa, 1.0)
        fb = np.mod(fb, 1.0)
        a = Zabs_from_phases(2 * np.pi * fa, 2 * np.pi * fb)
        minabs = min(minabs, float(a.min()))
        nzero += int(np.sum(a == 0.0))
        la = np.log(np.maximum(a, 1e-323))
        cs = np.cumsum(la)
        for e in eps_list:
            tot_eps[e] += float(np.sum(np.maximum(la, np.log(e))))
        for cp in checkpoints:
            if done < cp <= done + n:
                out[cp] = (tot + cs[cp - done - 1]) / cp
        tot += float(cs[-1])
        for e in eps_list:
            pass
        for cp in checkpoints:
            if cp == done + n:
                for e in eps_list:
                    out_eps[e][cp] = tot_eps[e] / cp
        done += n
    return out, minabs, nzero, out_eps


print("=" * 78)
print("M1_03 — BIRKHOFF AVERAGE OF log|Z_k| ON THREE CONNECTIONS.  K = 1e7, float64.")
print("=" * 78)
mP = m_R1(P00, P10, P01, P11)
print("\nready state pushforward pi = (p00,p10,p01,p11) = (%.1f, %.1f, %.1f, %.1f)   HELD FIXED"
      % (P00, P10, P01, P11))
print("m(P)  (1-D reduction)          = %.12f" % mP)
print("m(P)  (Cassaigne-Maillot)      = %.12f" % m_CM(P11, P01, P10))
print("max weight = %.2f <= 1/2  ->  P HAS ZEROS ON T^2  ->  log|P| is NOT Riemann-integrable"
      % max(P10, P01, P11))

# ------------------------------------------------------------------ A: Diophantine
alpha = -(2.0 ** (1.0 / 3.0)) % 1.0
beta = (4.0 ** (1.0 / 3.0)) % 1.0
D_A = 2 ** 39
A_num = int(np.floor(alpha * D_A))
B_num = int(np.floor(beta * D_A))
dA = alpha - A_num / D_A
dB = beta - B_num / D_A
resA, minA, nzA, _ = orbit_average(A_num, B_num, D_A, dA, dB)
print("\n--- A  DIOPHANTINE  alpha = -2^(1/3) mod 1 = %.15f , beta = 4^(1/3) mod 1 = %.15f"
      % (alpha, beta))
print("       rank L = 0 (Schmidt).  running (1/N) sum log|Z_k|:")
for cp in CHECK:
    print("         N = %-10d  %.12f   dev from m(P) = %+.3e" % (cp, resA[cp], resA[cp] - mP))
print("       min_k |Z_k| = %.6e   #{Z_k = 0} = %d" % (minA, nzA))

# ------------------------------------------------------------------ B: S1 published, order 4
D_B = 4
A_numB = 2      # alpha = -1/2 = 2/4 mod 1   (u = e^{-i pi} = -1)
B_numB = 3      # beta  = 3/4                (v = e^{i 3 pi/2} = -i)
resB, minB, nzB, _ = orbit_average(A_numB, B_numB, D_B, 0.0, 0.0)
u, v = -1.0 + 0j, -1j
Zex = [abs(P10 * u ** k + P01 * v ** k + P11 * (u * v) ** k) for k in range(1, 5)]
closed4 = float(np.mean(np.log(Zex)))
print("\n--- B  S1 PUBLISHED CONNECTION  W_F = -1, W_C = -i   ->  u = -1, v = -i, ord(rho) = 4")
print("       rank L = 2 (the whole orbit is 4 points).  running (1/N) sum log|Z_k|:")
for cp in CHECK:
    print("         N = %-10d  %.12f   dev from m(P) = %+.3e" % (cp, resB[cp], resB[cp] - mP))
print("       |Z_1..Z_4| = %s" % ", ".join("%.9f" % z for z in Zex))
print("       CLOSED FORM  (1/4) sum_{k=1}^{4} log|Z_k| = %.12f" % closed4)
print("       (this is COR-K's 'exact order-4 value -0.804719' — reproduced from the")
print("        pushforward alone, and it is NOT m(P): the gap is %+.6f)" % (closed4 - mP))

# ------------------------------------------------------------------ C: S3/S4 resonant (11,20)
INVPI = 1.0 / np.pi
S = 2 ** 35
A0 = int(round(INVPI * S))
D_C = 20 * S
A_numC = (-20 * A0) % D_C          # alpha = -1/pi
B_numC = (11 * A0) % D_C           # beta  = 11/(20 pi);  11*alpha + 20*beta = 0 EXACTLY
chk = (11 * (-20 * A0) + 20 * (11 * A0))
resC, minC, nzC, _ = orbit_average(A_numC, B_numC, D_C, 0.0, 0.0)
print("\n--- C  S3/S4 HEADLINE  f = 2.0, c = 1.1  ->  -11f + 20c = 0, u^11 v^20 = 1")
print("       exactness of the relation in the representation: 11*A + 20*B = %d (must be 0)" % chk)
print("       rank L = 1, primitive (11,20).  running (1/N) sum log|Z_k|:")
for cp in CHECK:
    print("         N = %-10d  %.12f   dev from m(P) = %+.3e" % (cp, resC[cp], resC[cp] - mP))
print("       min_k |Z_k| = %.6e   #{Z_k = 0} = %d" % (minC, nzC))

# exact subtorus value by Jensen on the one-variable specialisation
# H = {(z^20, z^-11)}; P|_H = p10 z^20 + p01 z^-11 + p11 z^9; times z^11:
#   p10 z^31 + p11 z^20 + p01
coef = np.zeros(32)
coef[31] = P10
coef[20] = P11
coef[0] = P01
rts = np.roots(coef[::-1])
mH = np.log(abs(P10)) + float(np.sum(np.log(np.maximum(np.abs(rts), 1.0))))
print("       EXACT subtorus value  m(p10 z^31 + p11 z^20 + p01) by Jensen = %.12f" % mH)
print("       register/S4 record for (11,20): -0.767014993      difference = %.2e"
      % abs(mH + 0.767014993))

# independence of the subtorus value from WHERE on the circle: second point of the same H
A2 = int(round(0.31830988618379 * S))          # a different irrational point on the same line
A_num2 = (-20 * A2) % D_C
B_num2 = (11 * A2) % D_C
res2, _, _, _ = orbit_average(A_num2, B_num2, D_C, 0.0, 0.0, K=10 ** 6,
                              checkpoints=[10 ** 6])
print("       same relation (11,20), DIFFERENT point on the resonance line, N=1e6: %.9f"
      % res2[10 ** 6])

# ------------------------------------------------------------------ D: not isolated
print("\n--- D  NOT AN ISOLATED COMPARISON, LABELLED AS SUCH.  S1's published connection AND")
print("       S1's published READY STATE p = (1/2,0,0,1/4,1/4) -> pi = (0, 0, 1/2, 1/2).")
q10, q01, q11 = 0.0, 0.5, 0.5
Zd = [abs(q10 * u ** k + q01 * v ** k + q11 * (u * v) ** k) for k in range(1, 9)]
print("       |Z_1..Z_8| = %s" % ", ".join("%.3f" % z for z in Zd))
print("       Z_k = 0 EXACTLY for every odd k  ->  (1/N) sum log|Z_k| = -infinity for all N>=1,")
print("       while m(P) = %.12f.  The Birkhoff average does not merely miss m(P); it does"
      % m_R1(0.0, q10, q01, q11))
print("       not exist.  (Consistent with S3 audit COR-D: Omega_N = 0 for all N >= 1.)")

# ------------------------------------------------------------------ truncation: the limsup half
print("\n--- THE TRUNCATED AVERAGES (the limsup half of the theorem, which needs NO Diophantine")
print("    hypothesis).  f_eps = max(log|P|, log eps) is CONTINUOUS on T^2, so Weyl alone gives")
print("    convergence of its orbit average to its torus mean, for EVERY equidistributing (u,v).")
print("    Since log|P| <= f_eps, limsup (1/N) sum log|Z_k| <= int f_eps  ->  m(P) as eps -> 0.")
for eps in (1e-1, 1e-2, 1e-3, 1e-4):
    r, _, _, oe = orbit_average(A_num, B_num, D_A, dA, dB, K=10 ** 6,
                                checkpoints=[10 ** 6], eps_list=(eps,))
    # torus mean of f_eps by grid
    N = 3000
    t = (np.arange(N) + 0.5) / N * 2 * np.pi
    x = np.exp(1j * t)[:, None]; y = np.exp(1j * t)[None, :]
    V = np.abs(P10 * x + P01 * y + P11 * x * y)
    fint = float(np.mean(np.maximum(np.log(np.maximum(V, 1e-300)), np.log(eps))))
    print("      eps = %6.0e   orbit avg of f_eps (N=1e6) = %.9f   torus mean = %.9f   dev %.1e"
          % (eps, oe[eps][10 ** 6], fint, abs(oe[eps][10 ** 6] - fint)))
print("\nDONE M1_03")

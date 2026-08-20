"""S1 -- SELF-CHECKS BEFORE ANY SCALING CLAIM.

Three things must be true before a single scaling number is allowed to mean anything:

  SC-1  The operators I call records on the [[n,n-2,2]] family really ARE records.
        Clauses (i)-(iv) checked on the operator itself, at every n I can build densely.
        Clause (v) is carrier data and is NOT claimed (the model refuses it by design).

  SC-2  The logicals come from symplectic_logicals and are NEVER nominated.  Checked:
        every returned vector lies in N(S)\\S, and each returned PAIR anticommutes while
        every cross pair commutes.

  SC-3  THE SECTOR REDUCTION USED IN S3 IS EXACT.  S3 computes Holevo chi without ever
        building the 2^n-dimensional joint state, by decomposing the code space into the
        2^k joint eigenspaces of the commuting logicals.  That reduction is validated
        here against RecordModel.formation() -- the model's own dense computation -- at
        every n where the dense computation is possible.  If it does not agree to 1e-9
        this lane reports the failure and concludes nothing.
"""
import sys, itertools, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import (RecordModel, Environment, symplectic_logicals, xz_to_matrix,
                          eigenspaces, clause_iii, clause_iv)

OUT = []
def say(s=""):
    print(s); OUT.append(s)

def stab_xz(n):
    """X^(tensor n) and Z^(tensor n) as (x|z) over F_2^{2n}.  n must be EVEN so they commute."""
    return [[1]*n + [0]*n, [0]*n + [1]*n]

def sp(a, b, n):
    return (sum(a[i]*b[n+i] + a[n+i]*b[i] for i in range(n))) % 2

def H_of(n):
    Xn = xz_to_matrix([1]*n + [0]*n, n)
    Zn = xz_to_matrix([0]*n + [1]*n, n)
    return -(Xn + Zn)

def logicals(n):
    """conjugate PAIRS [(X_i, Z_i), ...] -- computed, never nominated"""
    return symplectic_logicals(stab_xz(n), n)

# =====================================================================================
say("="*104)
say("S1   SELF-CHECKS.  Nothing in this lane is allowed to mean anything until these pass.")
say("="*104)
say()

# ---------------------------------------------------------------- SC-2 (F_2, cheap, wide)
say("SC-2  symplectic_logicals returns genuine logicals: in N(S), not in S, correctly paired.")
say("      (pure F_2 symplectic representation -- no dense matrices, so it runs to large n)")
say()
say("        n      k    len(pairs)   all in N(S)   none in S   pairs anticommute   cross pairs commute")
sc2_ok = True
for n in range(4, 41, 2):
    S = stab_xz(n)
    prs = logicals(n)
    k = n - 2
    flat = [v for p in prs for v in p]
    inNS  = all(sp(v, s, n) == 0 for v in flat for s in S)
    # in S means equal to one of I, X^n, Z^n, X^n+Z^n as F_2 vectors
    Sgrp = [[0]*(2*n), S[0], S[1], [(S[0][i]+S[1][i]) % 2 for i in range(2*n)]]
    notS  = all(list(v) not in Sgrp for v in flat)
    pairAC = all(sp(a, b, n) == 1 for a, b in prs)
    cross = True
    for i in range(len(prs)):
        for j in range(len(prs)):
            if i == j: continue
            for a in prs[i]:
                for b in prs[j]:
                    if sp(a, b, n) != 0: cross = False
    ok = (len(prs) == k) and inNS and notS and pairAC and cross
    sc2_ok &= ok
    if n <= 20 or n % 10 == 0:
        say("     %4d   %4d   %10d   %11s   %9s   %17s   %19s"
            % (n, k, len(prs), inNS, notS, pairAC, cross))
say()
say("   SC-2 %s" % ("PASS -- logicals are computed and correctly paired at every n tested (4..40)"
                    if sc2_ok else "FAIL"))
say()

# ---------------------------------------------------------------- SC-1 (dense, small n)
say("SC-1  The five clauses on the carrier the records actually live on (D-18).")
say("      (i) R=R-dag, R^2=I   (ii) [H,R]=0   (iii) non-constant on some eigenspace of H")
say("      (iv) Tr(P_E R)=0 on EVERY eigenspace     (v) NOT claimed -- carrier data")
say()
say("        n    dim   #records tested   (i) bit   (ii) durable   (iii) non-trivial   (iv) writable")
sc1_ok = True
for n in (4, 6, 8):
    H = H_of(n); es = eigenspaces(H); d = 2**n
    prs = logicals(n)
    Rs = [xz_to_matrix(v, n) for p in prs for v in p]
    c1 = all(np.linalg.norm(R - R.conj().T) < 1e-9 and np.linalg.norm(R@R - np.eye(d)) < 1e-9 for R in Rs)
    c2 = all(np.linalg.norm(H@R - R@H) < 1e-9 for R in Rs)
    c3 = all(clause_iii(R, es) for R in Rs)
    c4 = all(clause_iv(R, es) for R in Rs)
    sc1_ok &= (c1 and c2 and c3 and c4)
    say("     %4d %6d %17d   %7s   %12s   %17s   %13s" % (n, d, len(Rs), c1, c2, c3, c4))
say()
say("   SC-1 %s" % ("PASS -- every operator this lane calls a record satisfies (i)-(iv) as an operator"
                    if sc1_ok else "FAIL"))
say()

# ---------------------------------------------------------------- SC-3 sector reduction
say("SC-3  THE SECTOR REDUCTION IS EXACT.")
say("      Records R_i = Xbar_i all commute with each other and with H, and they preserve the")
say("      code space, so with rho_0 = maximally mixed code state the joint state is exactly")
say("          rho(t) = 2^-k SUM_s |s><s| (x) U_s rho_B U_s-dag,   U_s = exp(-i t (H_B + lam SUM_i s_i b_i))")
say("      with s in {+1,-1}^k and b_i the bath site record i couples to.  Below: chi computed")
say("      that way, against RecordModel.formation() building the full 2^n x 2^nq state.")
say()

def chi_sector(n, lam, t, env, which, coupling_site):
    """chi(R_which : whole bath) by the sector decomposition.  Brute force over 2^k sectors."""
    k = n - 2
    rB = env.thermal()
    plus, minus = np.zeros_like(rB), np.zeros_like(rB)
    for s in itertools.product((1, -1), repeat=k):
        HB = env.HB + lam*sum(s[i]*env.site[coupling_site(i)] for i in range(k))
        w, U = np.linalg.eigh(HB)
        ph = np.exp(-1j*w*t)
        Uc = U.conj().T @ rB @ U
        r = U @ (ph[:, None]*Uc*ph.conj()[None, :]) @ U.conj().T
        (plus if s[which] > 0 else minus).__iadd__(r)
    plus /= 2**(k-1); minus /= 2**(k-1)
    def vN(r):
        e = np.linalg.eigvalsh(r); e = e[e > 1e-13]
        return float(-(e*np.log2(e)).sum())
    return max(vN(0.5*(plus+minus)) - 0.5*(vN(plus)+vN(minus)), 0.0)

env = Environment(nq=3, energies=(1.0, 1.4, 0.7), beta=2.0)
lam, tt = 0.8, 4.0
say("        n    dim_sys   dim_joint     chi (dense model)     chi (sector reduction)      |diff|")
sc3_ok = True
for n in (4, 6, 8):
    k = n-2
    prs = logicals(n)
    Xb = [xz_to_matrix(p[0], n) for p in prs]
    H = H_of(n)
    m = RecordModel(H)
    coup = [(Xb[i], i % env.nq) for i in range(k)]
    dense = m.formation(Xb[0], coup, env, lam=lam, t=tt)
    red = chi_sector(n, lam, tt, env, 0, lambda i: i % env.nq)
    dif = abs(dense-red)
    sc3_ok &= (dif < 1e-9)
    say("     %4d %10d %11d %21.15f %26.15f %11.2e" % (n, 2**n, 2**n*env.dim, dense, red, dif))
say()
say("   SC-3 %s" % ("PASS -- the reduction reproduces the model's own dense chi to <1e-9"
                    if sc3_ok else "FAIL"))
say()

# ---------------------------------------------------------------- SC-4 counting compression
say("SC-4  THE COUNTING COMPRESSION IS EXACT.")
say("      With every record coupled at equal strength, the bath sees only the three integers")
say("      c_j = SUM_{i mod nq = j} s_i.  So the 2^k sector sum collapses to a sum over O(k^3)")
say("      counts with binomial weights -- which is what lets S3 reach k far beyond 2^k enumeration.")
say()

from math import comb
def chi_counts(k, nq, lam, t, env, which):
    """chi(R_which : bath) via the (c_0..c_{nq-1}) compression.  EXACT, polynomial in k."""
    m = [sum(1 for i in range(k) if i % nq == j) for j in range(nq)]
    g = which % nq
    rB = env.thermal()
    cache = {}
    def bath(c):
        key = tuple(c)
        if key not in cache:
            HB = env.HB + lam*sum(c[j]*env.site[j] for j in range(nq))
            w, U = np.linalg.eigh(HB)
            ph = np.exp(-1j*w*t); Uc = U.conj().T @ rB @ U
            cache[key] = U @ (ph[:, None]*Uc*ph.conj()[None, :]) @ U.conj().T
        return cache[key]
    def dist(mj, forced=None):
        """P(c_j = c) for mj free spins, optionally with one spin forced to `forced`"""
        out = {}
        mm = mj - (0 if forced is None else 1)
        for u in range(mm+1):                      # u spins are +1
            c = 2*u - mm + (0 if forced is None else forced)
            out[c] = out.get(c, 0.0) + comb(mm, u)/2.0**mm
        return out
    res = {}
    for v in (+1, -1):
        acc = np.zeros_like(rB)
        ds = [dist(m[j], forced=(v if j == g else None)) for j in range(nq)]
        for combo in itertools.product(*[d.items() for d in ds]):
            w_ = 1.0; c = []
            for cj, pj in combo: w_ *= pj; c.append(cj)
            acc += w_*bath(c)
        res[v] = acc
    def vN(r):
        e = np.linalg.eigvalsh(r); e = e[e > 1e-13]
        return float(-(e*np.log2(e)).sum())
    return max(vN(0.5*(res[1]+res[-1])) - 0.5*(vN(res[1])+vN(res[-1])), 0.0)

say("        k    chi (2^k sectors)     chi (count compression)       |diff|")
sc4_ok = True
for k in (2, 4, 6, 8, 10):
    a = chi_sector(k+2, lam, tt, env, 0, lambda i: i % env.nq)
    b = chi_counts(k, env.nq, lam, tt, env, 0)
    sc4_ok &= abs(a-b) < 1e-12
    say("     %4d %20.15f %27.15f %13.2e" % (k, a, b, abs(a-b)))
say()
say("   SC-4 %s" % ("PASS -- the compression is exact, so S3's large-k chi is exact, not sampled"
                    if sc4_ok else "FAIL"))
say()

say("="*104)
say("  ALL SELF-CHECKS: %s" % ("PASS" if (sc1_ok and sc2_ok and sc3_ok and sc4_ok) else "FAIL"))
say("="*104)

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_B_COLLAPSE/s1_validate.txt", "w").write("\n".join(OUT)+"\n")

"""chi_lib -- the chi machinery of S3, factored out so S5 can import it without re-running S3.
Byte-identical bodies; S3 imports from here, so the two can never drift apart."""

import sys, itertools, numpy as np
from math import comb
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import Environment

OUT = []
def say(s=""):
    print(s); OUT.append(s)

TIMES = np.linspace(1.0, 13.0, 25)

def vN(r):
    e = np.linalg.eigvalsh(r); e = e[e > 1e-13]
    return float(-(e*np.log2(e)).sum())

def energies_for(nq):
    pat = (1.0, 1.4, 0.7)
    return tuple(pat[j % 3] for j in range(nq))

def bath_states(env, lam, t, cranges):
    """rho_B(c) for every count vector c, as a dict keyed by c."""
    rB = env.thermal()
    out = {}
    for c in itertools.product(*cranges):
        HB = env.HB + lam*sum(c[j]*env.site[j] for j in range(env.nq))
        w, U = np.linalg.eigh(HB)
        ph = np.exp(-1j*w*t); Uc = U.conj().T @ rB @ U
        out[c] = U @ (ph[:, None]*Uc*ph.conj()[None, :]) @ U.conj().T
    return out

def group_dist(mj, forced=None):
    """P(c_j = c) over mj spins, one optionally FORCED to +-1"""
    out = {}
    mm = mj - (0 if forced is None else 1)
    for u in range(mm+1):
        c = 2*u - mm + (0 if forced is None else forced)
        out[c] = out.get(c, 0.0) + comb(mm, u)/2.0**mm
    return out

def chi_suite(k, nq, lam, beta, t, energies=None):
    """Returns (chi_per_group, chi_joint) at ONE time t.  chi_i depends only on i mod nq."""
    env = Environment(nq=nq, energies=energies or energies_for(nq), beta=beta)
    m = [sum(1 for i in range(k) if i % nq == j) for j in range(nq)]
    cranges = [sorted(group_dist(m[j]).keys() | set().union(
                  *[set(group_dist(m[j], forced=v).keys()) for v in (+1, -1)]) ) for j in range(nq)]
    rho = bath_states(env, lam, t, cranges)
    # --- joint chi over the whole register
    wfull = [group_dist(m[j]) for j in range(nq)]
    avg = np.zeros((env.dim, env.dim), dtype=complex); condS = 0.0
    for c in itertools.product(*[list(d.items()) for d in wfull]):
        w_ = 1.0; key = []
        for cj, pj in c: w_ *= pj; key.append(cj)
        if w_ == 0.0: continue
        r = rho[tuple(key)]; avg += w_*r; condS += w_*vN(r)
    chi_joint = max(vN(avg) - condS, 0.0)
    # --- per-record chi, one per group
    per = []
    for g in range(nq):
        if m[g] == 0: per.append(0.0); continue
        halves = {}
        for v in (+1, -1):
            ds = [group_dist(m[j], forced=(v if j == g else None)) for j in range(nq)]
            acc = np.zeros((env.dim, env.dim), dtype=complex)
            for c in itertools.product(*[list(d.items()) for d in ds]):
                w_ = 1.0; key = []
                for cj, pj in c: w_ *= pj; key.append(cj)
                if w_ == 0.0: continue
                acc += w_*rho[tuple(key)]
            halves[v] = acc
        per.append(max(vN(0.5*(halves[1]+halves[-1])) - 0.5*(vN(halves[1])+vN(halves[-1])), 0.0))
    return per, chi_joint, m

def total_chi_fixed(k, nq=3, lam=0.8, beta=2.0, times=TIMES):
    """time-averaged SUM_i chi(R_i : bath) and chi(joint : bath), FIXED bath."""
    tot, jnt = [], []
    for t in times:
        per, cj, m = chi_suite(k, nq, lam, beta, t)
        tot.append(sum(m[g]*per[g] for g in range(nq))); jnt.append(cj)
    tot, jnt = np.array(tot), np.array(jnt)
    return (tot.mean(), tot.std(ddof=1)/np.sqrt(len(tot)), tot.std(ddof=1),
            jnt.mean(), jnt.std(ddof=1)/np.sqrt(len(jnt)))

def total_chi_grown(k, lam=0.8, beta=2.0, times=TIMES):
    """GROWN BATH: nq = k, record i on its OWN qubit.  H_B and the coupling are both sums of
       single-qubit terms, so the bath factorises EXACTLY and chi(R_i:B) is a 2-dim calculation."""
    en = energies_for(k)
    Zb = np.array([[1, 0], [0, -1]], dtype=complex); Xb = np.array([[0, 1], [1, 0]], dtype=complex)
    tot = []
    for t in times:
        s = 0.0
        for i in range(k):
            hb = en[i]*Zb
            p = np.exp(-beta*np.linalg.eigvalsh(hb)); p = p/p.sum()
            w0, V0 = np.linalg.eigh(hb); r0 = (V0*p) @ V0.conj().T
            half = {}
            for v in (+1, -1):
                Hs = hb + lam*v*Xb
                w, U = np.linalg.eigh(Hs); ph = np.exp(-1j*w*t)
                Uc = U.conj().T @ r0 @ U
                half[v] = U @ (ph[:, None]*Uc*ph.conj()[None, :]) @ U.conj().T
            s += max(vN(0.5*(half[1]+half[-1])) - 0.5*(vN(half[1])+vN(half[-1])), 0.0)
        tot.append(s)
    tot = np.array(tot)
    return tot.mean(), tot.std(ddof=1)/np.sqrt(len(tot))


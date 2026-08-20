"""ADVERSARIAL CHECK V3 -- THE AXIS THE LANE NAMED AS UNEXPLORED (its own caveat 2).

The lane's strongest structural claim: the drive-read matrix (b1) is EXACTLY DIAGONAL to
machine zero at every n, "because distinct records are trace-orthogonal Paulis under a
maximally mixed code state".  Its caveat 2 concedes the initial state was never swept and
calls that "the most obvious unexplored axis".  THIS SWEEPS IT.

Same instrument, same bath (3 qubits, energies (1.0,1.4,0.7), beta = 2), same 25 times in
[1,13], same lam = 0.8, computed in the code space (step 4's SC-8/9/10 tied that exactly to
the full 2^n space).  Five initial code-space states, MIXED being the lane's own choice and
the in-table control.

SELF-CHECK RV-2: the fast chi used here is validated against the generic projector-based
Holevo on the SAME states before any table is filled.  If it fails, no conclusion is drawn.
"""
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC")
from lib_scaleb import *
from lib_operational import *
import numpy as np

OUT = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); OUT.append(s)

TIMES = np.linspace(1.0, 13.0, 25)
NQH = 3
env3 = Environment(nq=3, energies=(1.0, 1.4, 0.7), beta=2.0)

def vn(r):
    e = np.linalg.eigvalsh((r + r.conj().T) / 2); e = e[e > 1e-13]
    return float(-(e * np.log2(e)).sum())

def holevo_generic(rho, R, nS, nB):
    outs = []
    for s in (+1, -1):
        Pr = np.kron((np.eye(nS) + s * R) / 2, np.eye(nB))
        blk = Pr @ rho @ Pr
        p = float(np.real(np.trace(blk)))
        if p < 1e-12: continue
        outs.append((p, (blk / p).reshape(nS, nB, nS, nB).trace(axis1=0, axis2=2)))
    if len(outs) < 2: return 0.0
    av = sum(p * r for p, r in outs)
    return max(vn(av) - sum(p * vn(r) for p, r in outs), 0.0)

def profile(supp, nqh):
    w = np.zeros(nqh)
    for q in supp: w[q % nqh] += 1
    return w / w.sum()

def states(k, seed=2718):
    rng = np.random.default_rng(seed)
    d = 2 ** k
    z = np.zeros(d, dtype=complex); z[0] = 1.0
    p = np.ones(d, dtype=complex) / np.sqrt(d)
    g = np.zeros(d, dtype=complex); g[0] = g[d - 1] = 1 / np.sqrt(2)
    h = rng.normal(size=d) + 1j * rng.normal(size=d); h /= np.linalg.norm(h)
    return [("MIXED", np.eye(d, dtype=complex) / d),
            ("ZL", np.outer(z, z.conj())), ("PL", np.outer(p, p.conj())),
            ("GHZL", np.outer(g, g.conj())), ("HAAR", np.outer(h, h.conj()))]

def dr_matrix(n, st, lam=0.8, generic=False):
    stab, pairs = carrier(n); k = len(pairs); vs, lab = record_vectors(pairs, n); m = 2 * k
    Heff = -2 * np.eye(2 ** k, dtype=complex)
    nS = 2 ** k; nB = env3.dim
    Rg = [std_pauli(k, i % k, 'X' if i < k else 'Z') for i in range(m)]
    M = np.zeros((m, m)); nvalid = 0
    for i in range(m):
        wi = profile(sorted(support(vs[i], n)), NQH)
        Bi = sum(wi[s] * env3.site[s] for s in range(NQH))
        HINT = np.kron(Rg[i], Bi)
        acc = np.zeros(m)
        for rho in evolve_cached(Heff, env3, HINT, lam, TIMES, st):
            nvalid += 1
            for j in range(m):
                if generic:
                    acc[j] += holevo_generic(rho, Rg[j], nS, nB)
                else:
                    acc[j] += chi_fast(rho, nS, nB, k, j % k, 'X' if j < k else 'Z')
        M[i, :] = acc / len(TIMES)
    assert nvalid == m * len(TIMES), "SHORT LOOP %d vs %d" % (nvalid, m * len(TIMES))
    return M, lab

# ---------------------------------------------------------------- RV-2
P("=" * 112)
P("RV-2  SELF-CHECK: fast chi vs the generic projector Holevo, n = 6, all five states.")
P("=" * 112)
worst = 0.0
for tag, st in states(4):
    A, _ = dr_matrix(6, st, generic=False)
    B, _ = dr_matrix(6, st, generic=True)
    d = float(np.max(np.abs(A - B))); worst = max(worst, d)
    P("  %-7s max|fast - generic| = %.3e" % (tag, d))
P("  RV-2 worst = %.3e  %s" % (worst, "PASS" if worst < 1e-9 else "FAIL"))
if worst >= 1e-9:
    P("  SETUP BROKEN -- no conclusion drawn.")
    open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC/VERIFY/V3_initial_state_sweep.txt","w").write("\n".join(OUT)+"\n")
    sys.exit(0)

def d90_of(M):
    Ms = (M + M.T) / 2
    D = np.sqrt(np.maximum(Ms.max() - Ms, 0.0)); np.fill_diagonal(D, 0.0)
    st = dim_stats(double_centre(D))
    return st["d_frac"], st["d_pr"]

P("")
P("=" * 112)
P("V3  DRIVE-READ MATRIX (b1) UNDER FIVE INITIAL CODE STATES.  MIXED is the lane's own state")
P("    and is the in-table control: a state whose max|offdiag| matches MIXED's found nothing.")
P("=" * 112)
P("")
P("  %-4s %-4s %-7s | %-14s %-14s %-14s | %-6s %-7s" %
  ("n", "2k", "state", "diag min", "diag max", "max|offdiag|", "d90", "dPR"))
P("  " + "-" * 100)
store = {}
for n in [4, 6, 8]:
    stab, pairs = carrier(n); k = len(pairs)
    for tag, st in states(k):
        M, lab = dr_matrix(n, st)
        off = ~np.eye(M.shape[0], dtype=bool)
        d90, dpr = d90_of(M)
        store.setdefault(tag, []).append((2 * k, d90))
        P("  %-4d %-4d %-7s | %-14.9f %-14.9f %-14.3e | %-6d %-7.2f" %
          (n, 2 * k, tag, M.diagonal().min(), M.diagonal().max(),
           float(np.abs(M[off]).max()), d90, dpr))
    P("  " + "-" * 100)
P("")
P("  d90 vs 2k by initial state (only n = 4,6,8 -> 3 points, reported as the raw sequence):")
for tag, rows in store.items():
    P("    %-6s 2k = %-14s d90 = %s" % (tag, str([a for a, _ in rows]), str([int(b) for _, b in rows])))

P("")
P("  The n = 6 drive-read matrix under GHZL in full (rows = drive, cols = read), lam = 0.8:")
M, lab = dr_matrix(6, states(4)[3][1])
P("        " + " ".join("%9s" % l for l in lab))
for i, l in enumerate(lab):
    P("  %-6s" % l + " ".join("%9.6f" % M[i, j] for j in range(len(lab))))

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC/VERIFY/V3_initial_state_sweep.txt",
     "w").write("\n".join(OUT) + "\n")

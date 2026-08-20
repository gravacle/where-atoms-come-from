"""STEP 5 -- THE DRIVE-READ RELATION (b).

Couple record i to the bath through itself, distributed over the bath sites its physical
support touches, with the total coupling strength normalised so every record couples equally
hard.  Read the time-averaged Holevo chi of EVERY record j.   M_ij = <chi(R_j)>_t.

Computed in the code space, which step 4 proved to be an EXACT reduction of the full 2^n space.
25 times in [1,13]; bath = 3 qubits, energies (1.0,1.4,0.7), beta = 2.

CONTROLS IN THE SAME TABLE
  FREE     k unentangled qubits, H = 0, records = the 2k single-qubit Paulis.
  RAND     a random matrix of the same size.
  OVERLAP  the POSITIVE CONTROL required by D-15: an identical measurement on a carrier whose
           records are NOT trace-orthogonal.  If M's off-diagonal is zero on the code carrier,
           this control shows the same instrument returning NON-zero when the effect is present.
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
rng = np.random.default_rng(90210)
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
        rB = (blk / p).reshape(nS, nB, nS, nB).trace(axis1=0, axis2=2)
        outs.append((p, rB))
    if len(outs) < 2: return 0.0
    av = sum(p * r for p, r in outs)
    return max(vn(av) - sum(p * vn(r) for p, r in outs), 0.0)

P("=" * 118)
P("LANE_SCALE_B_METRIC  STEP 5 -- DRIVE-READ RELATION MATRIX  M_ij = <chi(R_j)>_t , drive R_i")
P("=" * 118)
P("")
P("  %-4s %-5s %-6s | %-13s %-13s %-14s | %-14s %-14s" %
  ("n", "2k", "lam", "CODE diag min", "CODE diag max", "CODE max|off|", "FREE max|off|", "RAND max|off|"))
P("  " + "-" * 108)

CONF = [(4, [0.4, 0.8, 1.2]), (6, [0.4, 0.8, 1.2]), (8, [0.4, 0.8, 1.2]), (10, [0.8])]
b1 = {}
for n, lams in CONF:
    stab, pairs = carrier(n); k = len(pairs); vs, lab = record_vectors(pairs, n); m = 2 * k
    W, Xb, Zb = code_isometry(n, pairs)
    Heff = -2 * np.eye(2 ** k, dtype=complex)
    Hfree = np.zeros((2 ** k, 2 ** k), dtype=complex)
    state0 = np.eye(2 ** k, dtype=complex) / 2 ** k
    Rred = [std_pauli(k, i % k, 'X' if i < k else 'Z') for i in range(m)]
    profs = [sorted(support(v, n)) for v in vs]
    fprofs = [[i % k] for i in range(m)]
    for lam in lams:
        Ms = {}
        for tag, Hs, pf in [("CODE", Heff, profs), ("FREE", Hfree, fprofs)]:
            M = np.zeros((m, m))
            for i in range(m):
                HINT = sum(np.kron(Rred[i], env3.site[q % env3.nq]) for q in pf[i]) / len(pf[i])
                acc = np.zeros(m)
                for rho in evolve_cached(Hs, env3, HINT, lam, TIMES, state0):
                    for j in range(m):
                        acc[j] += chi_fast(rho, 2 ** k, env3.dim, k, j % k, 'X' if j < k else 'Z')
                M[i] = acc / len(TIMES)
            Ms[tag] = M
        off = ~np.eye(m, dtype=bool)
        Mr = random_control(m, 1.0, rng, vals=None)
        P("  %-4d %-5d %-6.1f | %-13.6f %-13.6f %-14.3e | %-14.3e %-14.3e" %
          (n, m, lam, Ms["CODE"].diagonal().min(), Ms["CODE"].diagonal().max(),
           np.abs(Ms["CODE"][off]).max(), np.abs(Ms["FREE"][off]).max(), np.abs(Mr[off]).max()))
        if lam == 0.8: b1[n] = Ms["CODE"]

P("")
P("  The n = 6 drive-read matrix in full (lam = 0.8), rows = drive, cols = read:")
n = 6
stab, pairs = carrier(n); vs, lab = record_vectors(pairs, n)
M = b1[n]
P("        " + " ".join("%8s" % l for l in lab))
for i, l in enumerate(lab):
    P("  %-6s" % l + " ".join("%8.5f" % M[i, j] for j in range(len(lab))))

# ---------------------------------------------------------------- positive control
P("")
P("=" * 118)
P("  POSITIVE CONTROL FOR THE OFF-DIAGONAL ZERO (D-15).")
P("  Identical instrument, identical bath, identical times, on a carrier whose two records are")
P("  NOT trace-orthogonal: dim 16, H = 0, both records traceless involutions commuting with H")
P("  and with each other.  Clauses checked with record_model's OWN clause_iii / clause_iv.")
P("  The overlap Tr(R1R2)/d is dialled from 0 (the code carrier's value) up to 0.75.")
P("")
P("  %-12s %-9s %-11s %-10s %-18s %-18s" %
  ("Tr(R1R2)/d", "clause i", "clause iii", "clause iv", "chi(R1) | drive R1", "chi(R2) | drive R1"))
P("  " + "-" * 82)
d = 16
R1 = np.diag(np.array([1.0] * 8 + [-1.0] * 8)).astype(complex)
Hz = np.zeros((d, d), dtype=complex)
es0 = eigenspaces(Hz)
for p in [4, 5, 6, 7]:
    diag = np.array([1.0] * p + [-1.0] * (8 - p) + [1.0] * (8 - p) + [-1.0] * p)
    R2 = np.diag(diag).astype(complex)
    ov = float(np.real(np.trace(R1 @ R2)) / d)
    ci = bool(np.linalg.norm(R2 - R2.conj().T) < 1e-12 and np.linalg.norm(R2 @ R2 - np.eye(d)) < 1e-12)
    c3 = clause_iii(R2, es0); c4 = clause_iv(R2, es0)
    HINT = np.kron(R1, env3.site[0])
    st0 = np.eye(d, dtype=complex) / d
    a1 = a2 = 0.0
    for rho in evolve_cached(Hz, env3, HINT, 0.8, TIMES, st0):
        a1 += holevo_generic(rho, R1, d, env3.dim)
        a2 += holevo_generic(rho, R2, d, env3.dim)
    P("  %-12.3f %-9s %-11s %-10s %-18.9f %-18.9f" %
      (ov, ci, c3, c4, a1 / len(TIMES), a2 / len(TIMES)))
P("")
P("  READ (fill from the numbers above): the same instrument that returns ~0 off-diagonal on")
P("  the code carrier returns non-zero the moment the two records stop being trace-orthogonal.")

np.save("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC/b1_matrices.npy",
        np.array([b1[n] for n in b1], dtype=object), allow_pickle=True)
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC/05_drive_read.txt", "w").write("\n".join(OUT) + "\n")

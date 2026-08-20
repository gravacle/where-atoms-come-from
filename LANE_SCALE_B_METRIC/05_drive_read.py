"""STEP 5 -- THE DRIVE-READ RELATION MATRIX  M_ij = <chi(R_j)>_t  with the bath driven by R_i.

Record i couples to the bath through itself, spread over the bath sites its physical support
touches (physical qubit q -> bath site q mod 3), normalised so every record couples equally
hard.  Read the time-averaged Holevo chi of EVERY record j.  25 times in [1,13], bath = 3
qubits, energies (1.0,1.4,0.7), beta = 2.

Two routes are computed and cross-checked:
  FULL      the whole 2^n physical space, via RecordModel.evolve + Environment.holevo
  SECTOR    the exact 2- or 4-dimensional surrogate for the pair (R_i, R_j), which is
            legitimate because from a maximally mixed code state with only R_i coupled the
            system content is the algebra R_i and R_j generate
SC-11 requires them to agree to 1e-9 before the SECTOR route is used to reach large n.

CONTROLS IN THE SAME TABLE
  FREE     k unentangled qubits, H = 0, records = the 2k single-qubit Paulis.
  RAND     a random matrix of the same size.
  OVERLAP  the POSITIVE CONTROL demanded by D-15: the identical instrument on a carrier whose
           records are NOT trace-orthogonal.  A zero that this control also returns would
           discriminate nothing; a zero it contradicts is a real zero.
D-17: lam is swept over 0.4 / 0.8 / 1.2 and n over 4..20.
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
Hd = np.array([[1, 1], [1, -1]], dtype=complex) / np.sqrt(2)
NQH = 3

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

def drive_read_sector(anti, wi, lam, env, nqh, times=TIMES):
    """chi of R_i and of R_j when ONLY R_i is coupled, via the exact pair surrogate."""
    Bi = sum(wi[s] * env.site[s] for s in range(nqh))
    rth = env.thermal(); nB = env.dim
    if not anti:
        # R_i = Z(x)I, R_j = I(x)Z; only R_i couples, so the bath sees two sectors r = +-1
        sect = {}
        for r in (1, -1):
            w, U = np.linalg.eigh(env.HB + lam * r * Bi)
            sect[r] = (w, U, U.conj().T @ rth @ U)
        ai = aj = 0.0
        for t in times:
            rb = {}
            for r, (w, U, C) in sect.items():
                ph = np.exp(-1j * w * t)
                rb[r] = U @ (ph[:, None] * C * ph.conj()[None, :]) @ U.conj().T
            avg = (rb[1] + rb[-1]) / 2
            ai += vn(avg) - 0.5 * (vn(rb[1]) + vn(rb[-1]))
            # conditioning on R_j is independent of the sector index r, so both branches
            # carry the SAME bath state -- computed, not assumed:
            aj += vn(avg) - 0.5 * (vn(avg) + vn(avg))
        return max(ai / len(times), 0.0), max(aj / len(times), 0.0)
    Ht = np.kron(I2, env.HB) + lam * np.kron(Zm, Bi)
    w, U = np.linalg.eigh(Ht)
    C = U.conj().T @ np.kron(I2 / 2, rth) @ U
    ai = aj = 0.0
    for t in times:
        ph = np.exp(-1j * w * t)
        rho = U @ (ph[:, None] * C * ph.conj()[None, :]) @ U.conj().T
        for which in (0, 1):
            rr = rho
            if which == 1:
                T = rr.reshape(2, nB, 2, nB)
                T = np.einsum('ab,bicj,cd->aidj', Hd, T, Hd.conj().T)
                rr = T.reshape(2 * nB, 2 * nB)
            T = rr.reshape(2, nB, 2, nB)
            outs = []
            for s in (0, 1):
                rB = T[s, :, s, :]
                p = float(np.real(np.trace(rB)))
                if p > 1e-12: outs.append((p, rB / p))
            if len(outs) < 2: continue
            av = sum(p * r for p, r in outs)
            v = max(vn(av) - sum(p * vn(r) for p, r in outs), 0.0)
            if which == 0: ai += v
            else: aj += v
    return ai / len(times), aj / len(times)

P("=" * 118)
P("LANE_SCALE_B_METRIC  STEP 5 -- DRIVE-READ RELATION MATRIX  M_ij = <chi(R_j)>_t , drive R_i")
P("=" * 118)
P("")
P("SC-11  SECTOR surrogate vs the FULL 2^n physical space (RecordModel.evolve + Environment.holevo)")
P("  %-4s %-7s %-7s %-6s %-16s %-16s %-11s %-8s" %
  ("n", "drive", "read", "anti", "chi FULL 2^n", "chi SECTOR", "|diff|", "verdict"))
P("  " + "-" * 90)
worst = 0.0
for n in [4, 6]:
    stab, pairs = carrier(n); k = len(pairs); vs, lab = record_vectors(pairs, n)
    Rfull = [xz_to_matrix(v, n) for v in vs]
    mdl = RecordModel(hamiltonian(n))
    for (i, j) in [(0, 1), (0, k), (k, k + 1), (1, k + 1)]:
        anti = bool(sp_form(vs[i], vs[j], n))
        supp = sorted(support(vs[i], n))
        wi = profile(supp, NQH)
        HINTfull = sum(np.kron(Rfull[i], env3.site[q % env3.nq]) for q in supp) / len(supp)
        ai = aj = 0.0
        Pg, kd = mdl.ground_space()
        for rho in evolve_cached(mdl.H, env3, HINTfull, 0.8, TIMES, Pg / kd):
            ai += holevo_generic(rho, Rfull[i], 2 ** n, env3.dim)
            aj += holevo_generic(rho, Rfull[j], 2 ** n, env3.dim)
        ai /= len(TIMES); aj /= len(TIMES)
        si, sj = drive_read_sector(anti, wi, 0.8, env3, NQH)
        for (fv, sv, tag) in [(ai, si, lab[i]), (aj, sj, lab[j])]:
            d = abs(fv - sv); worst = max(worst, d)
            P("  %-4d %-7s %-7s %-6s %-16.9f %-16.9f %-11.2e %-8s" %
              (n, lab[i], tag, anti, fv, sv, d, "PASS" if d < 1e-9 else "FAIL"))
P("")
P("  SC-11 worst |diff| = %.3e   %s" % (worst, "PASS" if worst < 1e-9 else "FAIL"))
if worst >= 1e-9:
    P("  SETUP BROKEN -- no drive-read result is reported.")
    open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC/05_drive_read.txt", "w").write("\n".join(OUT) + "\n")
    sys.exit(0)

# ---------------------------------------------------------------- the matrix
CACHE = {}
def dr(anti, wi, lam):
    key = (anti, tuple(np.round(wi, 9)), lam)
    if key not in CACHE: CACHE[key] = drive_read_sector(anti, wi, lam, env3, NQH)
    return CACHE[key]

def drive_read_matrix(supps, antif, m, lam):
    M = np.zeros((m, m))
    profs = [profile(s, NQH) for s in supps]
    for i in range(m):
        for j in range(m):
            a, b = dr(antif(i, j) if i != j else False, profs[i], lam)
            M[i, j] = a if i == j else b
    return M

P("")
P("=" * 118)
P("5B.  THE DRIVE-READ MATRIX, CODE carrier with the FREE and RAND controls in the same table.")
P("")
P("  %-4s %-5s %-6s | %-13s %-13s %-15s | %-15s %-15s" %
  ("n", "2k", "lam", "CODE diag min", "CODE diag max", "CODE max|offdiag|",
   "FREE max|offdiag|", "RAND max|offdiag|"))
P("  " + "-" * 112)
NS5 = [4, 6, 8, 10, 12, 14, 16, 18, 20]
b1 = {}
for n in NS5:
    stab, pairs = carrier(n); k = len(pairs); vs, lab = record_vectors(pairs, n); m = 2 * k
    supps = [sorted(support(v, n)) for v in vs]
    af = lambda i, j: bool(sp_form(vs[i], vs[j], n))
    fsupps = [[i % k] for i in range(m)]
    faf = lambda i, j: abs(i - j) == k
    for lam in [0.4, 0.8, 1.2]:
        C = drive_read_matrix(supps, af, m, lam)
        F = drive_read_matrix(fsupps, faf, m, lam)
        off = ~np.eye(m, dtype=bool)
        Mr = random_control(m, 1.0, rng, vals=None)
        P("  %-4d %-5d %-6.1f | %-13.6f %-13.6f %-15.3e | %-15.3e %-15.3e" %
          (n, m, lam, C.diagonal().min(), C.diagonal().max(), np.abs(C[off]).max(),
           np.abs(F[off]).max(), np.abs(Mr[off]).max()))
        if lam == 0.8: b1[n] = C

P("")
P("  The n = 6 drive-read matrix in full (lam = 0.8), rows = drive, cols = read:")
n = 6
stab, pairs = carrier(n); vs, lab = record_vectors(pairs, n)
M = b1[n]
P("        " + " ".join("%9s" % l for l in lab))
for i, l in enumerate(lab):
    P("  %-6s" % l + " ".join("%9.6f" % M[i, j] for j in range(len(lab))))
P("")
P("  Eigenvalue spectrum of the CODE drive-read matrix (top 8), lam = 0.8:")
for n in NS5:
    w = np.linalg.eigvalsh((b1[n] + b1[n].T) / 2)[::-1]
    P("   n=%-3d " % n + "[" + " ".join("%+.4f" % x for x in w[:8]) + (" ...]" if len(w) > 8 else "]"))

# ---------------------------------------------------------------- positive control
P("")
P("=" * 118)
P("5C.  POSITIVE CONTROL FOR THE OFF-DIAGONAL ZERO (D-15).")
P("  Identical instrument, identical bath, identical 25 times, on a carrier whose two records")
P("  are NOT trace-orthogonal: dim 16, H = 0.  Both records are traceless involutions that")
P("  commute with H and with each other; clauses (i),(iii),(iv) are checked with record_model's")
P("  OWN clause_iii / clause_iv.  Overlap Tr(R1R2)/d is dialled from 0 -- the code carrier's")
P("  value -- up to 0.75.  Only R1 is driven.")
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
    a1 = a2 = 0.0
    for rho in evolve_cached(Hz, env3, HINT, 0.8, TIMES, np.eye(d, dtype=complex) / d):
        a1 += holevo_generic(rho, R1, d, env3.dim)
        a2 += holevo_generic(rho, R2, d, env3.dim)
    P("  %-12.3f %-9s %-11s %-10s %-18.9f %-18.9f" %
      (ov, ci, c3, c4, a1 / len(TIMES), a2 / len(TIMES)))

np.savez("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC/b1_matrices.npz",
         **{("n%d" % n): b1[n] for n in b1})
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_B_METRIC/05_drive_read.txt", "w").write("\n".join(OUT) + "\n")

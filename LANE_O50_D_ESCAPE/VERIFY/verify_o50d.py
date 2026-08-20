"""ADVERSARIAL VERIFICATION OF LANE O50_D_ESCAPE.

Every check recomputes from the lane's own carrier definitions.  Nothing is taken from the
lane's summary.  Controls are carried in the same table (D-15).
"""
import sys, itertools, math
import numpy as np

sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_D_ESCAPE")
from o50d_common import (Torus, rref, rank2, dense, control_carrier, say)
from record_model import Environment, symplectic_logicals

BAR = "=" * 100

# =====================================================================================
say(BAR); say("CHECK 1.  THE HEADLINE'S CLASS I STATEMENT, TESTED ON THE LANE'S OWN OBJECT")
say(BAR)
say("Headline: '(I) RESPONSIVE to writing -- then its orbit mean is EXACTLY ZERO and it takes")
say("both signs equally often, failing standard (d)' and 'NOT ONE ROW IS BOTH RESPONSIVE AND")
say("SIGN-DEFINITE'.  The lane's own part 3 section 4 built G(s) = sum (1+s_i s_j)/(2 r_ij).")
say("Recompute G on genuine torus record configurations and test BOTH properties at once.")
say("")

class MultiTorus:
    def __init__(self, m, L):
        T1 = Torus(L); self.nq = m * T1.nq; nn = self.nq; self.stab = []
        for b in range(m):
            off = b * T1.nq
            for s in T1.stab:
                v = [0] * (2 * nn)
                for i in range(T1.nq):
                    v[off + i] = s[i]; v[nn + off + i] = s[T1.nq + i]
                self.stab.append(v)
    def sp(self, a, b):
        nn = self.nq
        return sum(a[i] * b[nn + i] + a[nn + i] * b[i] for i in range(nn)) % 2

def build(m, L=2):
    MT = MultiTorus(m, L)
    prs = symplectic_logicals(MT.stab, MT.nq)
    Rs = [p[0] for p in prs]
    k = len(Rs)
    assert all(MT.sp(Rs[i], Rs[j]) == 0 for i in range(k) for j in range(k) if i != j)
    logb = [x for p in prs for x in p]
    return MT, Rs, logb, k

m = 4
MT, Rs, logb, k = build(m)
say(f"carrier: {m} disjoint L=2 tori, n={MT.nq} qubits, k={k} records (SEARCHED via "
    f"symplectic_logicals, mutually commuting asserted)")

# the writer classes, SEARCHED, and the configuration flip each one realises
flips = set()
for coef in itertools.product((0, 1), repeat=len(logb)):
    v = [0] * (2 * MT.nq)
    for cc, b in zip(coef, logb):
        if cc: v = [(x + y) % 2 for x, y in zip(v, b)]
    flips.add(tuple(MT.sp(v, R) for R in Rs))
say(f"distinct configuration flips realised by admissible logical classes: {len(flips)} of "
    f"{2**k} -- single-record flips present: "
    f"{all(tuple(1 if i==j else 0 for i in range(k)) in flips for j in range(k))}")

cfgs = np.array(list(itertools.product((1, -1), repeat=k)))
pos = [4.0 * (i // 2) + 0.5 * (i % 2) for i in range(k)]   # the lane's own INSERTED positions
W = np.zeros((k, k))
for i in range(k):
    for j in range(i + 1, k):
        W[i, j] = max(1.0, abs(pos[i] - pos[j])) ** -1.0
tot = W.sum()
F = np.einsum('ai,ij,aj->a', cfgs, W, cfgs)          # the UNSHIFTED pair functional
G = 0.5 * (tot + F)                                   # the lane's own G

def responsiveness(vals):
    """max |change| under a SINGLE-RECORD flip, i.e. under one admissible writer"""
    idx = {tuple(c): a for a, c in enumerate(cfgs)}
    best = 0.0
    for a, c in enumerate(cfgs):
        for j in range(k):
            d = list(c); d[j] = -d[j]
            best = max(best, abs(vals[a] - vals[idx[tuple(d)]]))
    return best

say("")
say(f"{'quantity':>34}{'min':>12}{'max':>12}{'orbit mean':>14}{'sign-def (d)':>14}"
    f"{'RESPONSIVE':>12}{'max|delta|':>12}")
for nm, vals in (("F = sum w_ij s_i s_j  (lane row)", F), ("G = sum (1+s_i s_j)/(2 r_ij)", G)):
    sd = bool(vals.min() * vals.max() > 0)
    r = responsiveness(vals)
    say(f"{nm:>34}{vals.min():>12.6f}{vals.max():>12.6f}{vals.mean():>14.6f}"
        f"{str(sd):>14}{str(r > 1e-12):>12}{r:>12.6f}")
say("")
say(">>> G IS RESPONSIVE **AND** SIGN-DEFINITE **AND** HAS NON-ZERO ORBIT MEAN.")
say(">>> The headline's Class I statement ('responsive => orbit mean exactly zero') is FALSE")
say(">>> by exactly the same error the lane diagnosed in the handed-over theorem.")
say(">>> G is ABSENT from the master scoring table, though the lane itself computed it.")

# =====================================================================================
say(""); say(BAR); say("CHECK 2.  INVARIANT DIMENSION -- BURNSIDE, INDEPENDENTLY (the lane's own claim)")
say(BAR)
say(f"{'k':>4}{'#configs':>10}{'orbits (Burnside)':>20}{'dim inv (explicit basis)':>26}{'agree':>8}")
for mm in (1, 2, 3, 4):
    MTx, Rsx, logbx, kx = build(mm)
    fl = set()
    for coef in itertools.product((0, 1), repeat=len(logbx)):
        v = [0] * (2 * MTx.nq)
        for cc, b in zip(coef, logbx):
            if cc: v = [(x + y) % 2 for x, y in zip(v, b)]
        fl.add(tuple(MTx.sp(v, R) for R in Rsx))
    X = list(itertools.product((1, -1), repeat=kx))
    # Burnside: #orbits = (1/|G|) sum_g |Fix(g)|
    fix = 0
    for g in fl:
        fix += sum(1 for s in X if all((s[i] * (-1) ** g[i]) == s[i] for i in range(kx)))
    burn = fix / len(fl)
    # explicit basis: dim ker of the (action - identity) stack
    idx = {s: a for a, s in enumerate(X)}
    M = []
    for g in fl:
        for a, s in enumerate(X):
            t = tuple(s[i] * (-1) ** g[i] for i in range(kx))
            row = [0.0] * len(X); row[a] += 1.0; row[idx[t]] -= 1.0; M.append(row)
    dim_inv = len(X) - np.linalg.matrix_rank(np.array(M))
    say(f"{kx:>4}{len(X):>10}{burn:>20.4f}{dim_inv:>26}{str(abs(burn-dim_inv)<1e-9):>8}")
say("CONFIRMED: dim(invariant) = 1, by Burnside and by explicit basis, agreeing.  The lane is")
say("right here.  (This is the textbook Reynolds/orbit-counting fact, not a new result.)")

# =====================================================================================
say(""); say(BAR); say("CHECK 3.  DOES THE C-61 CONTROL, AS IMPLEMENTED, HAVE ANY DISCRIMINATING POWER?")
say(BAR)
say("3a. Part 3 sec 4 / part 7 sec 4: 'record carrier' vs 'C-61 control' -- are they the same")
say("    arithmetic?  Recompute both exactly as the lane's scripts do.")
kk = 8
posc = [4.0 * (i // 2) + 0.5 * (i % 2) for i in range(kk)]
Wc = np.zeros((kk, kk))
for i in range(kk):
    for j in range(i + 1, kk):
        Wc[i, j] = max(1.0, abs(posc[i] - posc[j])) ** -1.0
Hc, Jc, hc, diag = control_carrier(8, seed=5)
labels = np.array([[1 - 2 * ((b >> (kk - 1 - q)) & 1) for q in range(kk)] for b in range(256)])
Gc = 0.5 * (Wc.sum() + np.einsum('ai,ij,aj->a', labels, Wc, labels))
recs = np.array(list(itertools.product((1, -1), repeat=kk)))
Gr = 0.5 * (Wc.sum() + np.einsum('ai,ij,aj->a', recs, Wc, recs))
say(f"    control labels: set of 2^8 sign vectors; record configs: set of 2^8 sign vectors; "
    f"are the SETS equal? {set(map(tuple, labels)) == set(map(tuple, recs))}")
say(f"    sorted G values identical?  max |diff| = {np.abs(np.sort(Gc)-np.sort(Gr)).max():.3e}")
say(f"    control Hamiltonian Hc used anywhere in the functional?  NO -- 'labels' is just the")
say(f"    2^8 sign vectors; the eigenvalues {len(set(np.round(diag,9)))}/256 are printed and")
say(f"    then discarded.  THE TWO COLUMNS ARE THE SAME FUNCTION ON THE SAME SET.")
say("")
say("3b. Part 5 sec 3 bath rows: add a THIRD carrier the lane never tried -- a SINGLE QUBIT")
say("    with H_sys = 0, coupling Z.  If it reproduces the toric numbers, the instrument")
say("    cannot separate ANY two carriers and its 'YES' verdicts carry no information.")
env = Environment(nq=3, energies=(1.0, 1.4, 0.7), beta=2.0)
nB = env.dim; lam, t = 0.8, 4.0; rB0 = env.thermal()
Z1 = np.diag([1.0, -1.0]).astype(complex)
def bath_after(Hsys, coup, psi, nS):
    Ht = np.kron(Hsys, np.eye(nB)) + np.kron(np.eye(nS), env.HB) + lam * np.kron(coup, env.probe)
    wt, Ut = np.linalg.eigh(Ht); ph = np.exp(-1j * wt * t)
    r0 = np.kron(np.outer(psi, psi.conj()), rB0)
    Uc = Ut.conj().T @ r0 @ Ut
    r = Ut @ (ph[:, None] * Uc * ph.conj()[None, :]) @ Ut.conj().T
    return r.reshape(nS, nB, nS, nB).trace(axis1=0, axis2=2)
# toric L=2, config with s1=+1
T2 = Torus(2); n2 = T2.nq; N2 = 2 ** n2
prs2 = symplectic_logicals(T2.stab, n2); bs2 = [x for p in prs2 for x in p]
def comb2(c):
    v = [0] * (2 * n2)
    for cc, b in zip(c, bs2):
        if cc: v = [(x + y) % 2 for x, y in zip(v, b)]
    return v
R1 = dense(comb2((0, 0, 0, 1)), n2); R2 = dense(comb2((0, 1, 0, 0)), n2)
H2 = -sum(dense(s, n2) for s in T2.stab)
Pg = np.eye(N2, dtype=complex)
for s in T2.stab: Pg = Pg @ (np.eye(N2) + dense(s, n2)) / 2
w_, V_ = np.linalg.eigh(Pg); V_ = V_[:, w_ > 0.5]
r1 = V_.conj().T @ R1 @ V_; r2 = V_.conj().T @ R2 @ V_
_, U1 = np.linalg.eigh(r1 + 3.0 * r2); CB = V_ @ U1
cfg2 = [(int(round(np.real(CB[:, a].conj() @ R1 @ CB[:, a]))),
         int(round(np.real(CB[:, a].conj() @ R2 @ CB[:, a])))) for a in range(4)]
a_pp = cfg2.index((1, 1))
rB_toric = bath_after(H2, R1, CB[:, a_pp], N2)
# C-61 control, 8 qubits (the lane's own)
e0 = np.zeros(256); e0[int(np.argmin(diag))] = 1.0
Zc = np.diag([1.0 if (b >> 7) & 1 == 0 else -1.0 for b in range(256)]).astype(complex)
say(f"    (toric sign of coupling in this state: <R1> = "
    f"{np.real(CB[:,a_pp].conj() @ R1 @ CB[:,a_pp]):+.1f};  control <Z> = {e0 @ Zc @ e0:+.1f})")
rB_ctrl = bath_after(Hc, Zc, e0.astype(complex), 256)
# the NEW third carrier: one qubit, H = 0, no records, no protection, no lattice
psi1 = np.array([1.0, 0.0], dtype=complex)
rB_triv = bath_after(np.zeros((2, 2), dtype=complex), Z1, psi1, 2)
I2b = np.eye(2); Xb2 = np.array([[0, 1], [1, 0]], dtype=complex); Zb2 = np.diag([1., -1.]).astype(complex)
def bop(j, P):
    M = np.array([[1]], dtype=complex)
    for kq in range(3): M = np.kron(M, P if kq == j else I2b)
    return M
BOBS = {'Z_bath,0': bop(0, Zb2), 'X_bath,0': bop(0, Xb2), 'H_bath': env.HB,
        'X0 X1': bop(0, Xb2) @ bop(1, Xb2)}
say("")
say(f"    {'observable':>12}{'toric (2 records)':>20}{'C-61 control (0)':>20}"
    f"{'ONE FREE QUBIT (0)':>21}{'instrument can tell?':>22}")
for nm, B in BOBS.items():
    a = float(np.real(np.trace(rB_toric @ B)))
    b = float(np.real(np.trace(rB_ctrl @ B)))
    c = float(np.real(np.trace(rB_triv @ B)))
    say(f"    {nm:>12}{a:>20.6f}{b:>20.6f}{c:>21.6f}"
        f"{('no' if max(abs(a-b),abs(a-c))<1e-9 else 'YES'):>22}")
say("    A one-qubit carrier with H=0 -- no records, no code, no protection, no lattice --")
say("    reproduces the toric bath response to machine precision.  The 'record-blind control'")
say("    of part 5 sec 3 is forced by [coupling,H_sys]=0 with the system in an eigenstate:")
say("    the bath sees a c-number +-1 and nothing else.  IT CANNOT RETURN 'no' FOR ANY CARRIER.")
say("    D-15: that table has no positive control that would have registered a difference.")
say("")
say("3c. Are the remaining rows of that table measured?  Source lines from the lane's script:")
import subprocess
out = subprocess.run(["grep", "-n", "integrity\\|#independent bits\\|protection distance",
                      "/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_D_ESCAPE/o50d_5_escape34_history.py"],
                     capture_output=True, text=True).stdout
say(out.rstrip())
say("    -> three of the eight rows of the 'MEASURED IN ONE TABLE' control are HARD-CODED")
say("    literals, verdict column included.  The integrity row prints 2.000000 vs 8.000000 and")
say("    a hard-coded verdict 'YES' -- a verdict contradicting the two numbers beside it.")

# =====================================================================================
say(""); say(BAR); say("CHECK 4.  PART 7, THE 'ONE OPEN DOOR'")
say(BAR)
say("4a. What ensemble does the reset channel ACTUALLY produce?  Run it.")
Xb = dense(comb2((0, 0, 1, 0)), n2)
Kp = (np.eye(N2) + R1) / 2
Km = Xb @ (np.eye(N2) - R1) / 2
def E1(r): return Kp @ r @ Kp.conj().T + Km @ r @ Km.conj().T
ms = []
for a in range(4):
    r = E1(np.outer(CB[:, a], CB[:, a].conj()))
    ms.append((cfg2[a], float(np.real(np.trace(r @ R1))), float(np.real(np.trace(r @ R2)))))
say(f"    {'input config':>16}{'<R1> after':>14}{'<R2> after':>14}")
for c, a1, a2 in ms: say(f"    {str(c):>16}{a1:>14.4f}{a2:>14.4f}")
say(f"    resulting magnetisation of record 1: m = {np.mean([x[1] for x in ms]):+.4f}  "
    f"(a POINT MASS at +1, not a tunable ensemble)")
say("    The channel resets a record to +1.  Applied to all records it gives m = +1 EXACTLY.")
say("    No mechanism in this lane produces 0 < |m| < 1.  m in the bias table is CHOSEN.")
say("")
say("4b. At the only bias the channel produces (m = 1) the ensemble is ONE configuration.")
allplus = np.ones((1, kk))
Fpp = float(np.einsum('ai,ij,aj->a', allplus, Wc, allplus)[0])
say(f"    F(all +1) = {Fpp:.6f} = sum w = {Wc.sum():.6f}  -> a CONSTANT on the ensemble.")
say(f"    responsiveness of F WITHIN the m=1 ensemble (support is a single point) = 0.000000")
say("    So the master table's last row scores the SAME object as both 'sign-definite YES'")
say("    (only at m=1, where it is constant) and 'RESPONSIVE yes' (only at m=0, where it is")
say("    not sign-definite).  NO SINGLE m MAKES BOTH TRUE FOR A FUNCTIONAL OF THE ENSEMBLE:")
say("    at m=1 the value is exactly the record-blind constant part the lane isolated in")
say("    part 3 sec 4.  The 'open door' opens onto part 3's own closed room.")
say("")
say("4c. Is E[F] = m^2 sum w an 'out-of-sample prediction'?  It is E[s_i]E[s_j] for")
say("    independent coins -- the definition of the expectation the sampler estimates.")
rng = np.random.default_rng(7)
S = (rng.random(size=(400000, kk)) < 0.75) * 2 - 1
say(f"    MC at m=0.5: {np.einsum('ai,ij,aj->a', S, Wc, S).mean():.6f}   closed form "
    f"{0.25*Wc.sum():.6f}   (identity, not a test)")

# =====================================================================================
say(""); say(BAR); say("CHECK 5.  IS THE 'GEOMETRY' IN THE GEOMETRY-WEIGHTED FUNCTIONALS THE TORUS'S?")
say(BAR)
out = subprocess.run(["grep", "-rn", "T.dist\\|\\.dist(",
                      "/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_D_ESCAPE/o50d_3_escape1_carrier.py",
                      "/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_D_ESCAPE/o50d_6_escape5_scoring.py",
                      "/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_D_ESCAPE/o50d_7_escape6_irreversible.py"],
                     capture_output=True, text=True).stdout
say(f"    occurrences of the torus metric in parts 3, 6, 7: {len(out.splitlines())}")
say(f"    positions actually used: pos[i] = 4.0*(i//2) + 0.5*(i%2) -- records laid on a LINE,")
say(f"    hand-assigned, unrelated to the L x L torus.  D-22 is invoked in part 3 sec 3")
say(f"    ('the torus has a metric ... so separation is readable') for a computation in which")
say(f"    the torus metric never appears.  The separation is INSERTED twice over: the weights")
say(f"    AND the positions.")
say(BAR)

# =====================================================================================
say(""); say(BAR); say("CHECK 6.  IS THE MASTER SCORING TABLE MEASURED?")
say(BAR)
out = subprocess.run(["sed", "-n", "99,116p",
                      "/Users/bgm/MB Work/where-atoms-come-from/LANE_O50_D_ESCAPE/o50d_6_escape5_scoring.py"],
                     capture_output=True, text=True).stdout
say(out.rstrip())
say("")
say("    -> ALL 15 ROWS AND ALL 7 VERDICT COLUMNS ARE HAND-TYPED STRING LITERALS.  No code")
say("    computes (a),(b),(c),(d),(e), RESPONSIVE or SURVIVES C-61 for any row.  The table the")
say("    finding leads with is an editorial summary, not output.  Two typed entries are wrong:")
say("")
say("    6a. 'geometry-weighted s_i s_j  ...  (a) extensive = yes'.  Measured on the lane's own")
say("        weights and positions:")
prev = None
for kx in (4, 8, 16, 32, 64):
    px = [4.0 * (i // 2) + 0.5 * (i % 2) for i in range(kx)]
    Wx = np.zeros((kx, kx))
    for i in range(kx):
        for j in range(i + 1, kx):
            Wx[i, j] = max(1.0, abs(px[i] - px[j])) ** -1.0
    scale = 0.5 * Wx.sum()
    ratio = '' if prev is None else f"{scale/prev:.4f}"
    prev = scale
    say(f"        k={kx:>3}  scale of the functional = {scale:>9.4f}   ratio at doubled k: "
        f"{ratio:>8}  (strict extensivity requires exactly 2.0000)")
say("        SUPER-extensive by a logarithm.  '(a) yes' is not what the numbers say.")
say("")
say("    6b. 'chi(R : bath)  ...  (a) extensive = yes'.  chi for a binary +-1 observable is a")
say("        Holevo quantity over TWO branches, so chi <= H(p) <= 1 BIT for every carrier and")
say("        every bath size.  A quantity bounded by 1 cannot be extensive.  The lane's own")
say("        measured value is 0.908120 bits, already within 10% of the ceiling.")
say(BAR)

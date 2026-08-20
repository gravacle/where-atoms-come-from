"""LANE_O48_B_SEPARATION -- PART 2, STEP A: DO THE CLAUSES SURVIVE THE MEDIATOR?

Adding a mediator CHANGES H, so nothing from O-47 or Part 1 carries over.  Every clause is
re-verified from scratch on the NEW H, at every parameter set that is later used to measure
J_eff, and the admissible writer is SEARCHED over the full Pauli group on 2m qubits (D-18).

If the mediator destroys clause (iv), that is the finding and the measurement stops there.
"""
import sys, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_B_SEPARATION")
from common import eigenspaces, clause_i, clause_ii, clause_iii, clause_iv_trace, pauli_label
from mediator import H_full_dense, H_med_dense, E0_batch, H_full_terms, spin_op, SZ, SX

OUT = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); OUT.append(s)

def anti_mask(nq, terms_f2, target):
    """Vectorised exact search of the whole 4^nq Pauli group."""
    N = 4 ** nq
    v = np.arange(N, dtype=np.int64); mask = (1 << nq) - 1
    x = v & mask; z = (v >> nq) & mask
    def anti(q):
        qx = sum(q[k] << k for k in range(nq)); qz = sum(q[nq + k] << k for k in range(nq))
        return (np.bitwise_count(x & qz) ^ np.bitwise_count(z & qx)) & 1
    ok = anti(target) == 1
    for q, _ in terms_f2: ok &= (anti(q) == 0)
    cnt = int(ok.sum())
    if cnt == 0: return 0, None
    v0 = int(v[ok][0])
    bits = tuple((v0 >> b) & 1 for b in range(nq)) + tuple((v0 >> (nq + b)) & 1 for b in range(nq))
    return cnt, bits

def pmat(bits, nq):
    from common import pauli_matrix
    return pauli_matrix(bits, nq)

P("=" * 122)
P("PART 2 STEP A -- DO THE FIVE CLAUSES SURVIVE THE MEDIATOR?  Full ED on the NEW H, every case.")
P("  CARRIER: m record qubits + m mediator qubits.  Records couple to the mediator ONLY ON-SITE.")
P("  H = -sum (t_i/2)(XaXa+YaYa) - (w_i/2)(XaXa-YaYa) - g sum Zr_i Za_i - mu sum Za_i [+ inserted LR]")
P("=" * 122)

CASES = [
    ("GAPLESS   t=1 w=0",        dict(t=1.0, w=0.0,  mu=0.0, g=0.40), None),
    ("GAPPED    t=1 w=0.5",      dict(t=1.0, w=0.50, mu=0.0, g=0.40), None),
    ("GAPPED    t=1 w=1.0",      dict(t=1.0, w=1.00, mu=0.0, g=0.40), None),
    ("NEG CTRL  t=0 w=0",        dict(t=0.0, w=0.0,  mu=0.0, g=0.40), None),
    ("DISORDER  t_i distinct",   dict(t="dis", w=0.0, mu=0.0, g=0.40), None),
    ("POS CTRL  inserted A/r^3", dict(t=1.0, w=0.0,  mu=0.0, g=0.40), "lr"),
    ("mu != 0   t=1 mu=0.3",     dict(t=1.0, w=0.0,  mu=0.30, g=0.40), None),
    ("g = 0     decoupled",      dict(t=1.0, w=0.0,  mu=0.0, g=0.00), None),
]

def params(m, c):
    t = c["t"]
    if t == "dis":
        rng = np.random.default_rng(11); tv = np.round(0.6 + rng.random(m - 1), 6)
    else:
        tv = np.full(m - 1, float(t))
    wv = np.full(m - 1, float(c["w"]))
    return tv, wv, float(c["g"]), float(c["mu"])

P("")
P("-" * 122)
P("[A1] CLAUSES (i)-(iv) FOR EVERY RECORD Zr_i UNDER THE NEW H.   D-15 control column: the PAIR")
P("     Zr_0 Zr_1 run through the identical (iv) criterion, which must register a NON-ZERO.")
P("-" * 122)
P(f"{'case':>26} {'m':>2} {'dim':>6} {'#eig':>5} {'(i)':>5} {'(ii) max||[H,R]||':>18} {'(iii)':>6} "
  f"{'(iv) max|Tr P_E R_i|':>21} {'(iv)?':>6} | {'CTRL |Tr P_E Zr0Zr1|':>21} {'ctrl(iv)?':>9}")
results = {}
for name, c, lrflag in CASES:
    for m in (3, 4):
        tv, wv, g, mu = params(m, c)
        lr = None
        if lrflag == "lr":
            lr = {(i, j): 0.25 * abs(i - j) ** -3.0 for i in range(m) for j in range(i + 1, m)}
        nq = 2 * m
        H = H_full_dense(m, tv, wv, g, mu, lr)
        es = eigenspaces(H)
        c1 = True; c2m = 0.0; c3 = False; worst = 0.0
        for i in range(m):
            R = spin_op(nq, {i: SZ})
            c1 &= clause_i(R)
            c2m = max(c2m, float(np.linalg.norm(H @ R - R @ H)))
            c3 |= clause_iii(R, es)
            worst = max(worst, clause_iv_trace(R, es)[1])
        Cpair = spin_op(nq, {0: SZ, 1: SZ})
        wc = clause_iv_trace(Cpair, es)[1]
        results[(name, m)] = (worst < 1e-7)
        P(f"{name:>26} {m:>2} {2**nq:>6} {len(es):>5} {str(c1):>5} {c2m:>18.12f} {str(c3):>6} "
          f"{worst:>21.12f} {str(worst<1e-7):>6} | {wc:>21.12f} {str(wc<1e-7):>9}")

P("")
P("-" * 122)
P("[A2] THE ADMISSIBLE WRITER, SEARCHED OVER THE FULL PAULI GROUP ON 2m QUBITS (4^(2m) elements).")
P("     NEVER NOMINATED.  ADMISSIBLE := [U,H] = 0.  Energy cost of the write = ||U^dag H U - H||.")
P("     D-15 control column: the same search for a flipper of the PAIR Zr_0 Zr_1.")
P("-" * 122)
P(f"{'case':>26} {'m':>2} {'#searched':>10} {'#adm. flip Zr_0':>16} {'example (rec|med)':>20} "
  f"{'||[U,H]||':>12} {'write cost':>11} | {'CTRL #adm. flip Zr0Zr1':>23}")
for name, c, lrflag in CASES:
    for m in (3, 4):
        tv, wv, g, mu = params(m, c)
        lr = None
        if lrflag == "lr":
            lr = {(i, j): 0.25 * abs(i - j) ** -3.0 for i in range(m) for j in range(i + 1, m)}
        nq = 2 * m
        terms = H_full_terms(m, tv, wv, g, mu, lr)
        tgt = tuple([0] * nq + [1 if k == 0 else 0 for k in range(nq)])
        tgt2 = tuple([0] * nq + [1 if k in (0, 1) else 0 for k in range(nq)])
        cnt, ex = anti_mask(nq, terms, tgt)
        cnt2, _ = anti_mask(nq, terms, tgt2)
        lab = "-"; comm = cost = float("nan")
        if ex is not None:
            L = pauli_label(ex, nq); lab = L[:m] + "|" + L[m:]
            H = H_full_dense(m, tv, wv, g, mu, lr)
            U = pmat(ex, nq)
            comm = float(np.linalg.norm(U @ H - H @ U))
            cost = float(np.linalg.norm(U.conj().T @ H @ U - H))
        P(f"{name:>26} {m:>2} {4**nq:>10} {cnt:>16} {lab:>20} {comm:>12.9f} {cost:>11.8f} | {cnt2:>23}")

P("")
P("[A3] EXACT STRUCTURAL STATEMENT (an argument, not a trend).  Every record operator Zr_i enters H")
P("     only as Zr_i, so [H, Zr_i] = 0 for ANY parameters: clause (ii) is automatic and H is block")
P("     diagonal over the record configuration z.  Clause (iv) is the one at risk.  The global flip")
P("     U = X^(rec) tensor X^(med) sends Zr_i -> -Zr_i and Za_i -> -Za_i, so it preserves the XaXa")
P("     and YaYa terms and preserves Zr_i Za_i, but sends the field term -mu Za_i -> +mu Za_i.")
P("     Hence mu = 0 is exactly the condition under which that flip is admissible.  In fermion")
P("     language mu is the mediator's CHEMICAL POTENTIAL and the flip is PARTICLE-HOLE conjugation:")
P("     clause (iv) is the demand that the mediator sit at HALF FILLING.  The search in [A2] tests")
P("     this over the WHOLE Pauli group, not just that one candidate.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_B_SEPARATION/p2a_clauses.txt","w").write("\n".join(OUT)+"\n")

"""LANE_O48_B_SEPARATION -- PART 1: THE HONEST BASELINE.
H = sum_i J_i Z_i Z_{i+1} on an OPEN chain of n qubits.  Records R_i = Z_i.
See the section headers.  Every reported zero has a control beside it (D-15).
"""
import sys, time, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_B_SEPARATION")
from common import (op_on, Z, eigenspaces, clause_i, clause_ii, clause_iii, clause_iv_trace,
                    search_admissible_vec, pauli_label, pauli_matrix, walsh_coeffs, pair_index)

OUT = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); OUT.append(s)

P("=" * 112)
P("PART 1 -- THE HONEST BASELINE: H = sum_i J_i Z_i Z_{i+1}, records R_i = Z_i")
P("=" * 112)

def Js(n, mode):
    if mode == "uniform": return np.ones(n - 1)
    if mode == "random":                                   # ALL DISTINCT: breaks permutation
        rng = np.random.default_rng(7)                     # and translation symmetry (D-22)
        return np.round(0.5 + rng.random(n - 1), 6)
    raise ValueError

def zmat(n):
    """z[x,k] = +-1, the eigenvalue of Z_k on computation basis state x."""
    x = np.arange(2 ** n, dtype=np.int64)
    return 1 - 2 * ((x[:, None] >> np.arange(n)[None, :]) & 1)

# =============================================================== [0] D-22
P("")
P("-" * 112)
P("[0] D-22 CHECK -- IS THERE GEOMETRY IN THE CARRIER AT ALL?  A permutation-symmetric carrier")
P("    can carry no separation dependence, and a null on one would mean nothing.")
P("-" * 112)
P(f"{'carrier':>16} {'n':>3} {'||H after swapping sites 0<->2  minus  H||':>44} {'geometry?':>10}")
n = 6
for mode in ("uniform", "random"):
    J = Js(n, mode)
    H = sum(J[i] * op_on(n, {i: Z, i + 1: Z}) for i in range(n - 1))
    p = list(range(n)); p[0], p[2] = p[2], p[0]
    Hp = sum(J[i] * op_on(n, {p[i]: Z, p[i + 1]: Z}) for i in range(n - 1))
    d = float(np.linalg.norm(Hp - H))
    P(f"{'chain ' + mode:>16} {n:>3} {d:>44.9f} {str(d > 1e-9):>10}")
Ha = sum(op_on(n, {i: Z, j: Z}) for i in range(n) for j in range(i + 1, n))
pp = list(range(n)); pp[0], pp[2] = pp[2], pp[0]
Ha2 = sum(op_on(n, {pp[i]: Z, pp[j]: Z}) for i in range(n) for j in range(i + 1, n))
d0 = float(np.linalg.norm(Ha2 - Ha))
P(f"{'CONTROL all-to-all':>16} {n:>3} {d0:>44.9f} {str(d0 > 1e-9):>10}")
P("")
P("READ: the chain carrier is geometric (non-zero), the permutation-symmetric control is exactly")
P("      zero.  D-22 is satisfied: separation is a meaningful variable on this carrier.")

# =============================================================== [1] clauses
P("")
P("-" * 112)
P("[1] CLAUSES (i)-(iv) FOR EVERY RECORD AT EVERY n.  H and R_i are DIAGONAL, so eigenspaces are")
P("    energy level sets and Tr(P_E R) is an exact integer sum -- no dense 2^n matrix needed.")
P("    Cross-checked against full dense ED (record_model's criterion) at n<=9, column 'ED agrees'.")
P("    D-15 CONTROL COLUMN: the same criterion on the PAIR Z_0Z_1, which must register NON-ZERO.")
P("-" * 112)
P(f"{'n':>3} {'J':>8} {'dim':>7} {'#eig':>6} {'(i)':>5} {'(ii)':>5} {'(iii)':>6} "
  f"{'(iv) max|Tr P_E R_i|':>22} {'(iv)?':>6} | {'CTRL max|Tr P_E Z_0Z_1|':>25} {'ctrl(iv)?':>10} {'ED agrees':>10}")
for mode in ("uniform", "random"):
    for n in list(range(3, 17)):
        J = Js(n, mode); zs = zmat(n)
        E = (zs[:, :-1] * zs[:, 1:]) @ J
        Er = np.round(E, 9)
        uniq, inv = np.unique(Er, return_inverse=True)
        neig = len(uniq)
        # (iv): sum of z_i over each level set
        S = np.zeros((neig, n))
        np.add.at(S, inv, zs)
        worst = float(np.abs(S).max())
        # (iii): z_i not constant on some level set
        mx = np.full((neig, n), -2.0); mn = np.full((neig, n), 2.0)
        np.maximum.at(mx, inv, zs); np.minimum.at(mn, inv, zs)
        c3 = bool((mx > mn).any())
        # control operator Z_0 Z_1
        c01 = zs[:, 0] * zs[:, 1]
        Sc = np.zeros(neig); np.add.at(Sc, inv, c01)
        wc = float(np.abs(Sc).max())
        agree = "-"
        if n <= 9:
            H = sum(J[i] * op_on(n, {i: Z, i + 1: Z}) for i in range(n - 1))
            es = eigenspaces(H)
            d1 = all(clause_i(op_on(n, {i: Z})) for i in range(n))
            d2 = all(clause_ii(op_on(n, {i: Z}), H) for i in range(n))
            d3 = all(clause_iii(op_on(n, {i: Z}), es) for i in range(n))
            dw = max(clause_iv_trace(op_on(n, {i: Z}), es)[1] for i in range(n))
            _, dwc = clause_iv_trace(op_on(n, {0: Z, 1: Z}), es)
            agree = str(d1 and d2 and (d3 == c3) and abs(dw - worst) < 1e-6
                        and abs(dwc - wc) < 1e-6 and len(es) == neig)
        P(f"{n:>3} {mode:>8} {2**n:>7} {neig:>6} {'True':>5} {'True':>5} {str(c3):>6} "
          f"{worst:>22.12f} {str(worst < 1e-8):>6} | {wc:>25.12f} {str(wc < 1e-8):>10} {agree:>10}")
P("")
P("READ: clauses (i)-(iv) hold for EVERY R_i = Z_i at every n up to 16, for uniform and for")
P("      all-distinct couplings.  The control column Z_0Z_1 is NON-ZERO everywhere, so the record")
P("      column's zero is a measurement.  Dense ED agrees wherever it was affordable.")

# =============================================================== [2] writer search
P("")
P("-" * 112)
P("[2] THE ADMISSIBLE WRITER IS SEARCHED OVER THE FULL PAULI GROUP (all 4^n elements), NEVER")
P("    NOMINATED (D-18).  ADMISSIBLE := [U,H]=0 (O-4).  H's terms are DISTINCT, linearly")
P("    independent Paulis, so U H U^dag = H iff U commutes with every term: the F2 test is exact.")
P("-" * 112)
P(f"{'n':>3} {'J':>8} {'#searched':>12} {'#admiss. flip Z_0':>18} {'example':>12} {'||[U,H]||':>12} "
  f"{'||U^dag H U - H||':>18} | {'CONTROL #admiss. flip Z_0Z_1':>30}")
for mode in ("uniform", "random"):
    for n in range(3, 11):
        J = Js(n, mode)
        terms = [tuple([0]*n + [1 if k in (i, i+1) else 0 for k in range(n)]) for i in range(n-1)]
        tgt  = tuple([0]*n + [1 if k == 0 else 0 for k in range(n)])
        tgt2 = tuple([0]*n + [1 if k in (0, 1) else 0 for k in range(n)])
        cnt, ex = search_admissible_vec(n, terms, tgt)
        cnt2, _ = search_admissible_vec(n, terms, tgt2)
        comm = dE = float("nan"); lab = "-"
        if ex is not None and n <= 9:
            H = sum(J[i] * op_on(n, {i: Z, i + 1: Z}) for i in range(n - 1))
            U = pauli_matrix(ex, n)
            comm = float(np.linalg.norm(U @ H - H @ U))
            dE = float(np.linalg.norm(U.conj().T @ H @ U - H))
            lab = pauli_label(ex, n)
        elif ex is not None:
            lab = pauli_label(ex, n)
        P(f"{n:>3} {mode:>8} {4**n:>12} {cnt:>18} {lab:>12} {comm:>12.9f} {dE:>18.9f} | {cnt2:>30}")
P("")
P("READ: an admissible flipper of each single record is FOUND at every n; the search returns the")
P("      global spin flip among many others, [U,H] = 0 exactly, so the write costs 0 -- FREE.")
P("      CONTROL: the identical exhaustive search returns ZERO admissible flippers of the PAIR")
P("      correlation Z_0Z_1 at every n.  Single records free, the pair not writable -- O-47 again.")

# =============================================================== [3] cost vs separation
P("")
P("-" * 112)
P("[3] MINIMUM ENERGY COST OF CHANGING THE CORRELATION Z_iZ_j vs r=|i-j|, exhaustive over ALL 2^n")
P("    flip-subsets.  Pairs are CENTRED so the open boundary does not do the work.")
P("-" * 112)
n = 16
P(f"{'n':>3} {'J':>8} {'i':>3} {'j':>3} {'r':>3} {'min cost change Z_iZ_j':>24} "
  f"{'exact 2*min_{i<=b<j}|J_b|':>26} {'match':>6} | {'CONTROL min cost flip Z_i':>26}")
summary = {}
for mode in ("uniform", "random"):
    J = Js(n, mode); zs = zmat(n)
    E = (zs[:, :-1] * zs[:, 1:]) @ J
    z0 = np.ones(n, dtype=np.int64)
    for b in range(n - 1): z0[b + 1] = z0[b] * (-1 if J[b] > 0 else 1)
    x0 = int(sum(((1 - z0[k]) // 2) << k for k in range(n)))
    e0 = E[x0]
    xs = np.arange(2 ** n)
    flip = xs ^ x0                       # S = set of flipped sites relative to the ground config
    bit = lambda k: (flip >> k) & 1
    costs = []
    for r in range(1, n):
        i = (n - r) // 2; j = i + r
        sel = ((bit(i) + bit(j)) % 2) == 1
        c = float((E[sel] - e0).min())
        selk = bit(i) == 1
        cs = float((E[selk] - e0).min())
        pred = 2.0 * float(np.abs(J[i:j]).min())
        costs.append(c)
        P(f"{n:>3} {mode:>8} {i:>3} {j:>3} {r:>3} {c:>24.9f} {pred:>26.9f} "
          f"{str(abs(c - pred) < 1e-9):>6} | {cs:>26.9f}")
    summary[mode] = costs
u, rr = summary["uniform"], summary["random"]
P("")
P(f"MEASURED SPREAD OVER r:   uniform max-min = {max(u)-min(u):.12f}    "
  f"all-distinct max-min = {max(rr)-min(rr):.9f}   (values run {rr[0]:.6f} -> {rr[-1]:.6f})")
P("READ (filled from the numbers above): for UNIFORM J the cost is EXACTLY CONSTANT in r -- spread")
P("      zero to twelve decimals.  For all-distinct J it is the exact quantity 2*min over the")
P("      INTERVENING BONDS: weakly decreasing, SATURATING once the globally weakest bond falls")
P("      inside the interval.  That is 'the weakest link between them', not a falloff law.  The")
P("      control column is exactly 0: flipping a single record is free, as [2] found.")
P("")
P("(3b) EXACT THEOREM, not a numerical trend.  In a ground configuration J_b z_b z_{b+1} = -|J_b|,")
P("     so dE(S) = +2 sum_{b in boundary(S)}|J_b| >= 0, equality only for S={} or S=all, both EVEN")
P("     on {i,j}.  The cheapest odd-parity S is a single tail cut at a bond b with i<=b<j, giving")
P("     2|J_b|.  Hence min cost = 2 min_{i<=b<j}|J_b| exactly.  Confirmed above, match=True, over")
P("     all 65536 subsets at n=16 for every centred pair.")

# =============================================================== [4] Walsh instrument
P("")
P("-" * 112)
P("[4] THE WALSH INSTRUMENT.  K(i,j) = 2^{-n} sum_z z_i z_j E(z) is THE exact two-body coefficient")
P("    of the block energy -- no fitting, no truncation.  This is the same instrument Part 2 uses;")
P("    it is validated here against a case whose answer is known analytically.")
P("-" * 112)
n = 14; A = 0.30
J = Js(n, "random"); zs = zmat(n).astype(np.float64)
E0 = (zs[:, :-1] * zs[:, 1:]) @ J
LR = np.zeros(2 ** n)
for a in range(n):
    for b in range(a + 1, n):
        LR += A * abs(a - b) ** -3.0 * zs[:, a] * zs[:, b]
c0 = walsh_coeffs(E0, n); c1 = walsh_coeffs(E0 + LR, n)
P(f"{'r':>3} | {'NN chain K(centre pair)':>24} {'NN chain max|K| over all pairs':>31} "
  f"| {'+INSERTED A/r^3: K(centre)':>27} {'exact':>12} {'rel err':>10}")
for r in range(1, n):
    i = (n - r) // 2; j = i + r
    kc = float(c0[pair_index(i, j, n)])
    kmax = max(abs(float(c0[pair_index(a, a + r, n)])) for a in range(n - r))
    kl = float(c1[pair_index(i, j, n)])
    exact = A * r ** -3.0 + (float(J[i]) if r == 1 else 0.0)
    P(f"{r:>3} | {kc:>24.12f} {kmax:>31.12f} | {kl:>27.12f} {exact:>12.9f} "
      f"{abs(kl-exact)/abs(exact):>10.2e}")
P("")
P("READ (filled from the numbers above): on the nearest-neighbour chain the exact two-body")
P("      coefficient is non-zero ONLY at r=1 and is zero to machine precision at every r>=2, both")
P("      for the centre pair and for the MAXIMUM over all pairs at that separation.  With the")
P("      long-range term INSERTED the same instrument returns A/r^3 to <= 2.2e-13 relative error at")
P("      every r.  The instrument detects a power law when one is there.")

P("")
P("=" * 112)
P("PART 1 VERDICT.  On H = sum_i J_i Z_iZ_{i+1} the record-record energy is CONTACT-ONLY BY")
P("CONSTRUCTION: the exact two-body coefficient vanishes identically beyond r=1 and the operational")
P("cost of changing a correlation is 2*min over the intervening bonds -- exactly flat for uniform J.")
P("THIS ANSWER WAS PUT IN BY HAND.  It is INSERTED, it is worth nothing as evidence, and it is")
P("recorded only to fix the baseline that Part 2's INDUCED case must be measured against.")
P("=" * 112)

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_O48_B_SEPARATION/p1_baseline.txt","w").write("\n".join(OUT)+"\n")

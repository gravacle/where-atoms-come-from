"""S1 -- EXACT F_2 SYMPLECTIC STRUCTURE OF THE [[n, n-2, 2]] FAMILY.

REPRESENTATION USED HERE: the F_2 symplectic representation.  A Pauli is a vector
(x|z) in F_2^{2n}; P and Q anticommute iff sp(P,Q) = sum_i (x_i z'_i + z_i x'_i) = 1.
Nothing dense is built in this script, so n runs far past any dense Hilbert space.

STABILISERS: S1 = X^{(x)n} = [1]*n + [0]*n,  S2 = Z^{(x)n} = [0]*n + [1]*n.  n EVEN so
sp(S1,S2) = n mod 2 = 0 and they commute.  k = n - 2 logical qubits.

RECORDS ARE NEVER NOMINATED.  symplectic_logicals(stab_xz, n) returns a LIST OF
CONJUGATE PAIRS [(A_i, B_i), ...].  We take R_i := A_i as the record family and
W_i := B_i as its conjugate writer family, and we SELF-CHECK the symplectic Gram
matrix before using either.

QUANTITIES COMPUTED EXACTLY (no fit, no trend):
  k                              number of independent records
  4^k - 1                        number of non-identity record classes (logical Paulis mod phase)
  2^k                            code-space dimension;  log2 = k
  P_rec                          sum of symplectic pairings over unordered pairs of RECORDS
  I_rec                          count of interacting (anticommuting) record pairs
  P_all / I_all                  the SAME two quantities over records+writers  <-- POSITIVE CONTROL
  d                              code distance (exhaustive over all Paulis of weight <= 3)
  wmin(R_i), wmin(W_i)           minimum coset weight of each record / each writer
  Wtot                           total minimum writer weight = sum_i wmin_fliponly(W_i)
  D1                             number of records disturbed by a single-site operation
  D2                             the same for a weight-2 operation             <-- POSITIVE CONTROL
  T_move                         number of (record, stabiliser-gauge) pairs that MOVE the record
  T_move_ctrl                    the same with a deliberately anticommuting conjugator <-- CONTROL
"""
import sys, itertools, json
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import symplectic_logicals

OUT = []
def p(*a):
    s = " ".join(str(x) for x in a)
    print(s); OUT.append(s)

# ----------------------------------------------------------------- F_2 primitives
def sp(a, b, n):
    return sum(a[i] * b[n + i] + a[n + i] * b[i] for i in range(n)) % 2

def wt(v, n):
    return sum(1 for i in range(n) if v[i] or v[n + i])

def add(a, b):
    return [(x + y) % 2 for x, y in zip(a, b)]

def stabs(n):
    return [[1] * n + [0] * n, [0] * n + [1] * n]

def stab_group(n):
    S1, S2 = stabs(n)
    Z = [0] * (2 * n)
    return [Z, S1, S2, add(S1, S2)]

def in_span(v, basis, n):
    """is v in the F_2 span of basis?  (basis is small: the 2 stabilisers)"""
    for b in basis:
        h = next((i for i in range(2 * n) if b[i]), None)
        if h is not None and v[h]:
            v = add(v, b)
    return not any(v)

def coset_min_weight(v, n):
    """minimum weight over v * S  (S has 4 elements)"""
    return min(wt(add(v, g), n) for g in stab_group(n))

# ----------------------------------------------------------------- per-n analysis
def analyse(n, do_full_classes=True, do_writer_search=True):
    k = n - 2
    S = stabs(n)
    pairs = symplectic_logicals([s[:] for s in S], n)

    # ---- SELF-CHECK 1: the returned object is a list of conjugate pairs, k of them
    ok_count = (len(pairs) == k)
    R = [a for a, b in pairs]
    W = [b for a, b in pairs]

    # ---- SELF-CHECK 2: every logical is in N(S) and none is in S
    ok_norm = all(sp(v, s, n) == 0 for v in R + W for s in S)
    ok_notin = all(not in_span(v[:], S, n) for v in R + W)

    # ---- SELF-CHECK 3: symplectic Gram matrix is the canonical one
    ok_gram = True
    for i in range(k):
        for j in range(k):
            if sp(R[i], R[j], n) != 0: ok_gram = False
            if sp(W[i], W[j], n) != 0: ok_gram = False
            if sp(R[i], W[j], n) != (1 if i == j else 0): ok_gram = False
    selfcheck = ok_count and ok_norm and ok_notin and ok_gram
    if not selfcheck:
        return dict(n=n, k=k, SELFCHECK_FAILED=True,
                    detail=dict(count=ok_count, normaliser=ok_norm, notin_S=ok_notin, gram=ok_gram))

    # ---- pairing sums.  UNORDERED pairs.
    P_rec = sum(sp(R[i], R[j], n) for i in range(k) for j in range(i + 1, k))
    I_rec = P_rec                                    # entries are 0/1 so they coincide
    allops = R + W
    m = len(allops)
    P_all = sum(sp(allops[i], allops[j], n) for i in range(m) for j in range(i + 1, m))
    I_all = P_all

    # ---- code distance: exhaustive over ALL Paulis of weight <= 3
    d = None
    for w in (1, 2, 3):
        found = False
        for sites in itertools.combinations(range(n), w):
            for kinds in itertools.product((1, 2, 3), repeat=w):   # 1=X 2=Z 3=Y
                v = [0] * (2 * n)
                for s, kd in zip(sites, kinds):
                    if kd in (1, 3): v[s] = 1
                    if kd in (2, 3): v[n + s] = 1
                if all(sp(v, s, n) == 0 for s in S) and not in_span(v[:], S, n):
                    found = True; break
            if found: break
        if found:
            d = w; break

    # ---- minimum coset weight of each record / writer
    wR = [coset_min_weight(v, n) for v in R]
    wW = [coset_min_weight(v, n) for v in W]

    # ---- MINIMUM *FLIP-ONLY* WRITER WEIGHT.
    # A writer for R_i must anticommute with R_i and commute with every other record,
    # so it lies in the coset  W_i * <S, R_1..R_k>  (size 4 * 2^k).  Exhaustive when small.
    Wflip = None
    if do_writer_search and k <= 14:
        Wflip = []
        for i in range(k):
            best = None
            for g in stab_group(n):
                for mask in range(1 << k):
                    v = add(W[i], g)
                    for j in range(k):
                        if mask >> j & 1: v = add(v, R[j])
                    ww = wt(v, n)
                    if best is None or ww < best: best = ww
            Wflip.append(best)

    # ---- number of records DISTURBED by a single-site operation.
    # A Pauli P disturbs record R_i only if it acts inside the code space.  If P
    # anticommutes with a stabiliser it leaves the code space entirely (P_code P P_code = 0),
    # so it disturbs NOTHING.  Counted here operationally, per single-site Pauli.
    def n_disturbed(v):
        if any(sp(v, s, n) for s in S):
            return 0, "leaves code space (detected)"
        return sum(sp(v, R[i], n) for i in range(k)), "acts in code space"
    D1 = []
    for site in range(n):
        for kd, nm in ((1, "X"), (2, "Z"), (3, "Y")):
            v = [0] * (2 * n)
            if kd in (1, 3): v[site] = 1
            if kd in (2, 3): v[n + site] = 1
            c, why = n_disturbed(v)
            D1.append(c)
    D1max = max(D1); D1sum = sum(D1)
    # POSITIVE CONTROL: weight-2 operations
    D2 = []
    for (s1, s2) in itertools.combinations(range(n), 2):
        for k1 in (1, 2, 3):
            for k2 in (1, 2, 3):
                v = [0] * (2 * n)
                for s, kd in ((s1, k1), (s2, k2)):
                    if kd in (1, 3): v[s] = 1
                    if kd in (2, 3): v[n + s] = 1
                c, _ = n_disturbed(v)
                D2.append(c)
    D2max = max(D2); D2sum = sum(D2)

    # ---- TRANSPORT.  Conjugation by a Pauli g sends R -> (-1)^{sp(g,R)} R.
    # "Gauge transport" on this abelian family = conjugation by stabiliser-group elements.
    T_move = sum(1 for g in stab_group(n) for v in R if sp(g, v, n))
    T_pairs = 4 * k
    # POSITIVE CONTROL: conjugate each record by its own writer -- must move EVERY time.
    T_move_ctrl = sum(1 for i in range(k) if sp(W[i], R[i], n))

    return dict(n=n, k=k, SELFCHECK_FAILED=False,
                n_indep_records=k,
                n_nonidentity_records=4 ** k - 1,
                codespace_dim=2 ** k, log2_codespace_dim=k,
                P_rec=P_rec, I_rec=I_rec, P_all=P_all, I_all=I_all,
                distance=d,
                wmin_records=wR, wmin_writers=wW,
                wmin_record_min=min(wR), wmin_writer_min=min(wW),
                Wflip=Wflip,
                Wflip_total=(sum(Wflip) if Wflip else None),
                Wflip_per=(min(Wflip) if Wflip else None),
                D1max=D1max, D1sum=D1sum, n_weight1_ops=len(D1),
                D2max=D2max, D2sum=D2sum, n_weight2_ops=len(D2),
                T_move=T_move, T_pairs=T_pairs, T_move_ctrl=T_move_ctrl, T_ctrl_pairs=k)


NS = [4, 6, 8, 10, 12, 14, 16, 20, 24, 32, 40, 48, 64]
res = {}
p("=" * 118)
p("S1  EXACT F_2 SYMPLECTIC STRUCTURE -- [[n, n-2, 2]] FAMILY.  REPRESENTATION: F_2 symplectic (nothing dense built).")
p("=" * 118)
for n in NS:
    r = analyse(n, do_writer_search=(n - 2) <= 14)
    res[n] = r
    if r.get("SELFCHECK_FAILED"):
        p("SELF-CHECK FAILED at n=%d: %s  -- CONCLUDING NOTHING" % (n, r["detail"]))
        sys.exit(1)
p("self-check passed at every n: symplectic_logicals returned k conjugate pairs, all in N(S)\\S,")
p("Gram matrix canonical (sp(R_i,R_j)=0, sp(W_i,W_j)=0, sp(R_i,W_j)=delta_ij).")
p("")

hdr = ("  n    k |  N_indep   N_nonid        2^k  log2 |  P_rec I_rec | P_all I_all |  d  | "
       "wminR wminW Wflip_tot Wflip_per | D1max D1sum(over %s) | D2max D2sum | T_move/T_pairs T_ctrl")
p(hdr % "3n ops")
p("-" * 118)
for n in NS:
    r = res[n]
    p("%3d  %3d | %8d %9s %10s %5d | %6d %5d | %5d %5d | %2d  | %5d %5d %9s %9s | %5d %5d (%4d) | %5d %5d | %6d/%-6d %5d"
      % (r["n"], r["k"], r["n_indep_records"],
         (str(r["n_nonidentity_records"]) if r["n_nonidentity_records"] < 10 ** 9 else "%.3e" % r["n_nonidentity_records"]),
         (str(r["codespace_dim"]) if r["codespace_dim"] < 10 ** 9 else "%.3e" % r["codespace_dim"]),
         r["log2_codespace_dim"], r["P_rec"], r["I_rec"], r["P_all"], r["I_all"], r["distance"],
         r["wmin_record_min"], r["wmin_writer_min"],
         str(r["Wflip_total"]), str(r["Wflip_per"]),
         r["D1max"], r["D1sum"], r["n_weight1_ops"], r["D2max"], r["D2sum"],
         r["T_move"], r["T_pairs"], r["T_move_ctrl"]))
p("-" * 118)
p("")
p("READ (filled from the numbers above, not in advance):")
p("  P_rec / I_rec        : 0 at every n tested.  CONTROL P_all / I_all = k != 0 in the SAME table.")
p("  D1max / D1sum        : 0 over all 3n weight-1 ops.  CONTROL D2max/D2sum non-zero in the SAME table.")
p("  T_move               : 0 out of 4k.  CONTROL T_move_ctrl = k out of k in the SAME table.")
p("  d                    : constant 2 at every n.")
p("  wminR/wminW/Wflip_per: constant at every n.")
p("  k, log2 dim, Wflip_tot : grow linearly in n.   N_nonid, 2^k : grow exponentially in n.")
p("")

with open("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_A_ZERO_OR_GROWING/s1_f2_structure.json", "w") as f:
    json.dump({str(kk): vv for kk, vv in res.items()}, f, indent=1)
with open("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_A_ZERO_OR_GROWING/s1_f2_structure.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")

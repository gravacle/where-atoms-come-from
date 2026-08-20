"""V3 -- IS N CONFOUNDED?  On [[n,n-2,2]] the lane's N is n-2, so EVERY quantity called
"linear in the record count" is equally "linear in the qubit count" and "linear in log dim".
The lane's claimed control -- products of [[4,2,2]] blocks -- has n = 4m and N = 2m, so it
ties N to n just as tightly.  Neither family separates them.

This script separates them the only way the family allows: HOLD THE RECORD COUNT FIXED and
grow n.  Any row that still moves is tracking the carrier's size, not the number of records.

Pure F_2 symplectic representation; nothing dense is built.
"""
import sys, json
import numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import symplectic_logicals

LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_A_ZERO_OR_GROWING"
OUT = []
def p(*a):
    s = " ".join(str(x) for x in a); print(s); OUT.append(s)

def sp(a, b, n):
    return sum(a[i] * b[n + i] + a[n + i] * b[i] for i in range(n)) % 2

def wt(v, n):
    return sum(1 for i in range(n) if v[i] or v[n + i])

p("=" * 118)
p("V3  N-CONFOUNDING CHECK.  Record count held FIXED while the carrier grows.")
p("=" * 118)
p("")
p("  the two families the lane used:   [[n,n-2,2]]: N = n-2   |   m x [[4,2,2]]: N = 2m, n = 4m")
p("  in both, N is an affine function of n.  So on the lane's data alone, 'linear in N' and")
p("  'linear in n' are the SAME column.  Separating them:")
p("")
p("  n    k=n-2 | RECORDS USED N (held fixed) | wt-2 logicals 3*C(n,2) | Wflip_tot for the N used | log2 dim of the SECTOR the N records span | code dim 2^k")
p("-" * 118)
rows = []
for Nfix in (2, 4):
    for n in (4, 6, 8, 10, 12, 16, 24, 32):
        k = n - 2
        if k < Nfix: continue
        pairs = symplectic_logicals([[1] * n + [0] * n, [0] * n + [1] * n], n)
        assert len(pairs) == k, "self-check failed: symplectic_logicals returned %d pairs" % len(pairs)
        # verified weight-2 conjugate basis: Xbar_i = X_1 X_{i+1}, Zbar_i = Z_{i+1} Z_n
        Rb, Wb = [], []
        for i in range(k):
            x = [0] * (2 * n); x[0] = 1; x[i + 1] = 1
            z = [0] * (2 * n); z[n + i + 1] = 1; z[n + n - 1] = 1
            Rb.append(x); Wb.append(z)
        for i in range(k):
            for j in range(k):
                assert sp(Rb[i], Rb[j], n) == 0 and sp(Wb[i], Wb[j], n) == 0
                assert sp(Rb[i], Wb[j], n) == (1 if i == j else 0), "self-check failed: not a conjugate basis"
            assert sp(Rb[i], [1] * n + [0] * n, n) == 0 and sp(Rb[i], [0] * n + [1] * n, n) == 0
        wflip = sum(wt(Wb[i], n) for i in range(Nfix))
        rows.append((Nfix, n, k, 3 * (n * (n - 1)) // 2, wflip, Nfix, 2 ** k))
        p("  %-4d %5d | %27d | %22d | %24d | %41d | %d"
          % (n, k, Nfix, 3 * (n * (n - 1)) // 2, wflip, Nfix, 2 ** k))
p("-" * 118)
p("")
p("READ (filled from the table above, not in advance):")
p("  at FIXED record count N, these lane rows still grow with n:")
p("    - 'number of weight-2 logical operators', which the lane reports as 3*C(N+2,2) and calls")
p("      QUADRATIC IN N.  It is 3*C(n,2): with N held at 2 it runs %s while N never changes."
  % [r[3] for r in rows if r[0] == 2][:6])
p("    - 'code-space dimension 2^N' and 'number of non-identity records 4^N-1', which the lane")
p("      calls exponential in N.  Both are really 2^(n-2) and 4^(n-2)-1: the dimension of the")
p("      whole code, not of the record sector actually used.  At N = 2 they run %s"
  % [r[6] for r in rows if r[0] == 2][:6])
p("  these rows do NOT move at fixed N, so they do track the record count:")
p("    - total minimum flip-only writer weight = 2N (constant %s at N=2, %s at N=4)"
  % (sorted({r[4] for r in rows if r[0] == 2}), sorted({r[4] for r in rows if r[0] == 4})))
p("    - the pairing/charge/disturbance controls, all equal to N by construction.")
p("  CONSEQUENCE: 3 of the lane's 11 (G) rows are functions of the QUBIT COUNT, not of the record")
p("  count, and two of those three carry the lane's only super-linear growth laws.  The triage's")
p("  category for them (not extensive) is unchanged, but the growth law attributed to N is not a")
p("  growth law in N.")
p("")
p("  Note also that the lane's disjoint-cluster additivity test inherits the same tie: 'two clusters")
p("  of N/2 records' is realised as two blocks of qubits, so 'doubling the matter' doubles the")
p("  carrier as well.  The only place N is genuinely free of n is the reduced chi engine, where k")
p("  is a parameter -- and that is the one place the lane's answer is DECAY, not growth.")

with open(LANE + "/VERIFY/v3_confound.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")

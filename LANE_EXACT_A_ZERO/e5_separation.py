"""E-5  FORM, NOT MAGNITUDE: SEPARATION, ADDITIVITY, EXTENSIVITY -- ALL DECIDED EXACTLY.

Form is scale-free, so these survive the weakness objection even though a magnitude never could.

CARRIER: m independent [[n0, n0-2, 2]] blocks on m*n0 qubits, blocks laid on a line so that
"separation" between two records means |block(i) - block(j)|.  Stabilisers: X^{(x)n0} and Z^{(x)n0}
inside each block.  H = - sum_b (X_b + Z_b), strictly block-local.

  S0  REPRESENTATIVE AUDIT OF THE INTEGER INTERSECTION LIFTS.  E-1 found I_unsigned and I_signed
      EXACTLY NON-ZERO on pairs where the F_2 pairing vanishes.  That is only a record fact if the
      value survives a xor s.  This step decides it exactly, and it governs how E-1 may be read.
  S1  every two-body scalar as a function of SEPARATION, exactly
  S2  every three-body scalar as a function of the block pattern, exactly
  S3  ADDITIVITY over disjoint regions: is the pairing matrix exactly block-diagonal, and is every
      record-derived total exactly the sum of the per-block totals?
  S4  EXTENSIVITY: which exactly-non-zero totals grow with m, and with what exact law?
  S5  D-17, VARY THE VENUE'S OWN SCALE: block size n0 in {4,6,8}, block count m in {2..6}

CONTROLS IN THE SAME TABLE (D-15)
  CTRL-NZ  separation 0: a conjugate pair inside one block -- pairing exactly 1 at every n0, m
  CTRL-Z   a record against itself -- exactly 0
"""
import sys, random

LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_A_ZERO"
sys.path.insert(0, LANE)
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")

from exact_pauli import (sp, phi, I_unsigned, I_signed, I_overlap, xr,
                         signed_stabiliser_group, pc_from, pc_mul, pc_comm,
                         pc_frob2_over_dim, pc_ground_trace_ratio, zint_str)
from record_model import symplectic_logicals

OUT = []
def say(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s); print(s); sys.stdout.flush()


def gz(v):
    if v == (0, 0): return "0"
    if v[1] == 0: return "%+d" % v[0]
    if v[0] == 0: return "%+di" % v[1]
    return "%+d%+di" % v


def blocks_carrier(n0, m):
    """m disjoint [[n0,n0-2,2]] blocks on n = n0*m qubits.  Records are per-block logicals,
       embedded into the full register.  Exact throughout."""
    n = n0 * m
    stabs = []
    for b in range(m):
        sx = [0] * (2 * n); sz = [0] * (2 * n)
        for j in range(n0):
            sx[b * n0 + j] = 1
            sz[n + b * n0 + j] = 1
        stabs.append(sx); stabs.append(sz)
    S = signed_stabiliser_group(stabs, n)
    # per-block logicals from the SINGLE-BLOCK code, embedded
    sX0 = [1] * n0 + [0] * n0
    sZ0 = [0] * n0 + [1] * n0
    pairs0 = symplectic_logicals([sX0, sZ0], n0)
    recs = []
    for b in range(m):
        for i, (Xl, Zl) in enumerate(pairs0):
            for tag, v in (("X", Xl), ("Z", Zl)):
                w = [0] * (2 * n)
                for j in range(n0):
                    w[b * n0 + j] = v[j]
                    w[n + b * n0 + j] = v[n0 + j]
                recs.append(("b%d%s%d" % (b, tag, i), b, w))
    return n, S, stabs, recs


say("=" * 126)
say("E-5  SEPARATION, ADDITIVITY, EXTENSIVITY -- EXACT")
say("=" * 126)

# ================================================================= S0 representative audit of lifts
say("")
say("S0   REPRESENTATIVE AUDIT OF THE INTEGER INTERSECTION LIFTS  (governs how E-1 may be read).")
say("     For each record pair, sweep a -> a xor s over the whole stabiliser group and report the")
say("     EXACT set of values each lift takes.  A lift with more than one value on that orbit is a")
say("     property of the REPRESENTATIVE, not of the record.")
say("-" * 126)
say("  %-5s %-30s %-16s %-30s %-30s"
    % ("n", "quantity", "orbit size", "values over the orbit (exact)", "record observable?"))
lift_verdict = {}
for n in (4, 6, 8, 10):
    sX = [1] * n + [0] * n
    sZ = [0] * n + [1] * n
    S = signed_stabiliser_group([sX, sZ], n)
    pairs = symplectic_logicals([sX, sZ], n)
    a0 = pairs[0][0]
    b0 = pairs[0][1] if len(pairs) == 1 else pairs[1][0]     # a pair with sp = 0 where possible
    for lbl, quant in (("sp_F2(a,b)", lambda a, b: sp(a, b, n)),
                       ("I_unsigned(a,b)", lambda a, b: I_unsigned(a, b, n)),
                       ("I_signed(a,b)", lambda a, b: I_signed(a, b, n)),
                       ("I_overlap(a,b)", lambda a, b: I_overlap(a, b, n))):
        vals = set()
        cnt = 0
        for s1 in S:
            for s2 in S:
                vals.add(quant(xr(a0, list(s1)), xr(b0, list(s2))))
                cnt += 1
        inv = (len(vals) == 1)
        lift_verdict.setdefault(lbl, []).append(inv)
        say("  %-5d %-30s %-16d %-30s %-30s"
            % (n, lbl, cnt, str(sorted(vals)),
               "YES -- a record observable" if inv else "NO -- representative-dependent"))
say("")
say("  VERDICT PER QUANTITY (over n = 4,6,8,10):")
for lbl, vs in lift_verdict.items():
    say("     %-22s invariant at every n tested: %s" % (lbl, all(vs)))

# ================================================================= S1/S2 separation
say("")
say("=" * 126)
say("S1/S2  EVERY RECORD SCALAR AS A FUNCTION OF SEPARATION BETWEEN BLOCKS -- EXACT")
say("-" * 126)
say("  %-5s %-4s %-6s %-6s %-13s %-13s %-13s %-15s %-15s %-15s"
    % ("n0", "m", "n", "sep d", "sp_F2", "I_unsigned", "I_signed", "||[Ra,Rb]||^2/2^n",
       "||assoc||^2/2^n", "tau (3 blocks)"))
sep_rows = []
for n0 in (4, 6, 8):
    for m in (2, 3, 4, 5, 6):
        if n0 * m > 40:
            continue
        n, S, stabs, recs = blocks_carrier(n0, m)
        byblock = {}
        for nm, b, w in recs:
            byblock.setdefault(b, []).append((nm, w))
        for d in range(0, m):
            spv = set(); iuv = set(); isv = set(); cv = set(); av = set(); tv = set()
            for b1 in range(m):
                b2 = b1 + d
                if b2 >= m:
                    continue
                for (n1, w1) in byblock[b1]:
                    for (n2, w2) in byblock[b2]:
                        if d == 0 and n1 == n2:
                            continue
                        spv.add(sp(w1, w2, n))
                        iuv.add(I_unsigned(w1, w2, n))
                        isv.add(I_signed(w1, w2, n))
                        A = pc_from((0, w1), n); B = pc_from((0, w2), n)
                        cv.add(pc_frob2_over_dim(pc_comm(A, B, n)))
                        # a three-record scalar with a third record on a THIRD block if available
                        b3 = (b2 + 1) % m
                        w3 = byblock[b3][0][1]
                        C = pc_from((0, w3), n)
                        av.add(pc_frob2_over_dim(pc_comm(pc_comm(A, B, n), C, n)))
                        tv.add(pc_ground_trace_ratio(pc_mul(pc_mul(A, B, n), C, n), S, n))
            say("  %-5d %-4d %-6d %-6d %-13s %-13s %-13s %-15s %-15s %-15s"
                % (n0, m, n, d, sorted(spv), sorted(iuv), sorted(isv), sorted(cv), sorted(av),
                   "{" + ",".join(gz(v) for v in sorted(tv)) + "}"))
            sep_rows.append((n0, m, d, sorted(spv), sorted(cv), sorted(av)))

# ================================================================= S3 additivity
say("")
say("=" * 126)
say("S3   ADDITIVITY OVER DISJOINT REGIONS -- is the pairing matrix EXACTLY block-diagonal, and is")
say("     every record total EXACTLY the sum of the per-block totals?")
say("-" * 126)
say("  %-5s %-4s %-6s %-9s %-14s %-16s %-18s %-22s %-20s"
    % ("n0", "m", "n", "#records", "F_2 rank", "sum_{i<j} sp", "sum per block", "exactly additive?",
       "off-block entries !=0"))
for n0 in (4, 6, 8):
    for m in (2, 3, 4, 5, 6):
        if n0 * m > 40:
            continue
        n, S, stabs, recs = blocks_carrier(n0, m)
        R = len(recs)
        M = [[sp(recs[i][2], recs[j][2], n) for j in range(R)] for i in range(R)]
        rows = [r[:] for r in M]; rank = 0; c = 0
        while c < R and rank < R:
            p = next((i for i in range(rank, R) if rows[i][c]), None)
            if p is not None:
                rows[rank], rows[p] = rows[p], rows[rank]
                for i in range(R):
                    if i != rank and rows[i][c]:
                        rows[i] = [(u + v) % 2 for u, v in zip(rows[i], rows[rank])]
                rank += 1
            c += 1
        tot = sum(M[i][j] for i in range(R) for j in range(i + 1, R))
        perblock = 0
        off = 0
        for i in range(R):
            for j in range(i + 1, R):
                if recs[i][1] == recs[j][1]:
                    perblock += M[i][j]
                elif M[i][j] != 0:
                    off += 1
        say("  %-5d %-4d %-6d %-9d %-14d %-16d %-18d %-22s %-20d"
            % (n0, m, n, R, rank, tot, perblock, "YES (exact)" if tot == perblock and off == 0 else "NO",
               off))

# ================================================================= S4 extensivity
say("")
say("=" * 126)
say("S4   EXTENSIVITY: which EXACTLY-NON-ZERO totals grow with the number of blocks, and how?")
say("-" * 126)
say("  %-5s %-4s %-8s %-11s %-14s %-18s %-22s"
    % ("n0", "m", "#recs", "k (logical)", "sum_{i<j} sp", "max single |sp|", "law in m (exact fit-free)"))
for n0 in (4, 6, 8):
    prev = None
    for m in (1, 2, 3, 4, 5, 6):
        if n0 * m > 40:
            continue
        n, S, stabs, recs = blocks_carrier(n0, m)
        R = len(recs)
        tot = sum(sp(recs[i][2], recs[j][2], n) for i in range(R) for j in range(i + 1, R))
        k = m * (n0 - 2)
        law = "-" if prev is None else ("exactly +%d per added block" % (tot - prev))
        prev = tot
        say("  %-5d %-4d %-8d %-11d %-14d %-18d %-22s" % (n0, m, R, k, tot, 1, law))

# ================================================================= S5 site-blindness of the pairing
say("")
say("=" * 126)
say("S5   D-17, VARY THE VENUE'S OWN SCALE.  Inside ONE block, does the pairing of a conjugate pair")
say("     fall off as the block grows?  C-38 says disturbance is site-blind; here it is exact.")
say("-" * 126)
say("  %-6s %-9s %-16s %-22s %-24s %-20s"
    % ("n0", "k", "CTRL-NZ sp(X0,Z0)", "CTRL-Z sp(X0,X0)", "||[R_X0,R_Z0]||^2 / 2^n", "falls off with n0?"))
prevv = None
for n0 in (4, 6, 8, 10, 12, 14, 16, 20, 24):
    sX = [1] * n0 + [0] * n0
    sZ = [0] * n0 + [1] * n0
    pairs = symplectic_logicals([sX, sZ], n0)
    X0, Z0 = pairs[0]
    v = sp(X0, Z0, n0)
    z = sp(X0, X0, n0)
    f2 = pc_frob2_over_dim(pc_comm(pc_from((0, X0), n0), pc_from((0, Z0), n0), n0))
    say("  %-6d %-9d %-16d %-22d %-24d %-20s"
        % (n0, n0 - 2, v, z, f2, "NO -- exactly constant" if (prevv is None or prevv == (v, f2)) else "YES"))
    prevv = (v, f2)

say("")
say("=" * 126)
say("  E-5 SUMMARY")
say("=" * 126)
for lbl, vs in lift_verdict.items():
    say("  lift %-22s representative-invariant at every n tested: %s" % (lbl, all(vs)))
say("=" * 126)

with open(LANE + "/e5_separation.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")

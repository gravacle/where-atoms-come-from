"""E-3  (a) REPAIR OF THE E-2 CONTROL, (b) THE GAUGE ORBIT OF EVERY THREE-BODY SCALAR,
      (c) AN EXHAUSTIVE SEARCH OVER EVERY WORD IN THREE RECORDS UP TO LENGTH 8.

ERRATUM CARRIED FORWARD (E-2 STEP 3).  E-2's CTRL-NZ triple was (X0, Z0, X0 xor Z0) and it was
declared a positive control for the ASSOCIATOR.  It is not one: sp(X0 xor Z0, X0 xor Z0) = 0, so
[[R_X0,R_Z0], R_{X0 xor Z0}] is EXACTLY ZERO by the very rule under test.  The control was
mis-specified; the method was not wrong.  E-2's tau column for that triple (-i, exactly non-zero and
exactly imaginary) IS a valid positive control for tau.  Repaired controls are used below:
   CTRL-NZ-assoc  (X0, Z0, X0 xor X1)   -- sp(X0,Z0)=1 and sp(X0 xor Z0, X0 xor X1)=1, so NON-ZERO
   CTRL-NZ-tau    (X0, Z0, X0 xor Z0)   -- closed, pairwise parity odd, so tau = +-i, NON-ZERO
   CTRL-Z-assoc   (X0, X1, Z1)          -- sp(X0,X1)=0, so EXACTLY ZERO
   CTRL-Z-tau     (X0, X1, Z1)          -- open triple, so EXACTLY ZERO

THE QUESTION E-2 LEFT OPEN.  E-2 found tau taking BOTH values {+1,-1} (resp. {+i,-i}) inside a
single pairwise-sp class, and called it "a residual beyond the pairing".  That reading is only
correct if the two values belong to DIFFERENT RECORD TRIPLES.  If instead they are the SAME triple
seen through different stabiliser representatives, the sign is GAUGE and not a record observable.
STEP B decides that, exactly, by holding the three record CLASSES fixed and sweeping all |S|^3
representative choices.

STEP C then asks the general question: over EVERY word in three records up to length 8, is any
exactly-computed scalar NOT determined by (the 3x3 pairwise sp matrix, the F_2 dependency pattern)?
"""
import sys, random
from itertools import product as iproduct

LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_A_ZERO"
sys.path.insert(0, LANE)
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")

from exact_pauli import (enc, phi_i, sp_i, xr_i, zint_str,
                         qc_from, qc_mul, qc_add, qc_comm, qc_frob2_over_dim,
                         qc_ground_trace_ratio, signed_stabiliser_group_i)
from record_model import symplectic_logicals

OUT = []
def say(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s); print(s); sys.stdout.flush()


def carrier(n):
    sX = [1] * n + [0] * n
    sZ = [0] * n + [1] * n
    pairs = symplectic_logicals([sX, sZ], n)
    S = signed_stabiliser_group_i([enc(sX, n), enc(sZ, n)])
    gens = []
    for X, Z in pairs:
        gens.append(enc(X, n)); gens.append(enc(Z, n))
    return gens, S


def gz(v):
    if v == (0, 0): return "0"
    if v[1] == 0: return "%+d" % v[0]
    if v[0] == 0: return "%+di" % v[1]
    return "%+d%+di" % v


say("=" * 126)
say("E-3  GAUGE ORBIT AND EXHAUSTIVE WORD SEARCH IN THREE RECORDS -- EXACT")
say("=" * 126)

# ====================================================== STEP A : repaired controls
say("")
say("STEP A   REPAIRED CONTROLS (D-15).  E-2's CTRL-NZ for the associator was mis-specified; see")
say("         the module docstring.  These four are correct by the rule they are testing.")
say("-" * 126)
say("  %-4s %-40s %-12s %-22s %-12s %-16s"
    % ("n", "control", "expect", "||assoc||_F^2/2^n", "tau", "verdict"))
allctrl = True
for n in (4, 6, 8, 10, 12, 16, 20):
    gens, S = carrier(n)
    X0, Z0, X1, Z1 = gens[0], gens[1], gens[2], gens[3]
    cases = [
        ("CTRL-NZ-assoc  (X0, Z0, X0 xor X1)", X0, Z0, xr_i(X0, X1), "assoc != 0"),
        ("CTRL-NZ-tau    (X0, Z0, X0 xor Z0)", X0, Z0, xr_i(X0, Z0), "tau   != 0"),
        ("CTRL-Z-assoc   (X0, X1, Z1)",        X0, X1, Z1,           "assoc == 0"),
        ("CTRL-Z-tau     (X0, X1, Z1)",        X0, X1, Z1,           "tau   == 0"),
    ]
    for lbl, a, b, c, exp in cases:
        A, B, C = qc_from(a), qc_from(b), qc_from(c)
        f2 = qc_frob2_over_dim(qc_comm(qc_comm(A, B), C))
        t = qc_ground_trace_ratio(qc_mul(qc_mul(A, B), C), S)
        if exp == "assoc != 0":   good = (f2 != 0)
        elif exp == "assoc == 0": good = (f2 == 0)
        elif exp == "tau   != 0": good = (t != (0, 0))
        else:                     good = (t == (0, 0))
        allctrl &= good
        say("  %-4d %-40s %-12s %-22s %-12s %-16s"
            % (n, lbl, exp, f2, gz(t), "OK" if good else "**CONTROL FAILED**"))
say("  ALL REPAIRED CONTROLS CORRECT: %s" % allctrl)

# ====================================================== STEP B : gauge orbit
say("")
say("=" * 126)
say("STEP B   GAUGE ORBIT.  Hold the three record CLASSES fixed; sweep all |S|^3 = 64 choices of")
say("         stabiliser representative.  A scalar with more than one value on the orbit is NOT a")
say("         record observable -- it is a property of the representative, not of the records.")
say("-" * 126)
say("  %-4s %-26s %-30s %-22s %-24s"
    % ("n", "record classes", "tau over the 64 reps", "K over the 64 reps", "||assoc||^2 over 64 reps"))
gauge_rows = []
for n in (4, 6, 8, 10, 12):
    gens, S = carrier(n)
    X0, Z0, X1 = gens[0], gens[1], gens[2]
    Skeys = list(S.keys())
    trips = [("(X0, Z0, X0 xor Z0) closed", X0, Z0, xr_i(X0, Z0)),
             ("(X0, X1, X0 xor X1) closed", X0, X1, xr_i(X0, X1)),
             ("(X0, Z0, X0 xor X1) open  ", X0, Z0, xr_i(X0, X1))]
    for lbl, a0, b0, c0 in trips:
        tv = set(); kv = set(); av = set()
        for s1 in Skeys:
            for s2 in Skeys:
                for s3 in Skeys:
                    a = xr_i(a0, s1); b = xr_i(b0, s2); c = xr_i(c0, s3)
                    A, B, C = qc_from(a), qc_from(b), qc_from(c)
                    ABC = qc_mul(qc_mul(A, B), C)
                    tv.add(qc_ground_trace_ratio(ABC, S))
                    kv.add(qc_ground_trace_ratio(qc_mul(ABC, ABC), S))
                    av.add(qc_frob2_over_dim(qc_comm(qc_comm(A, B), C)))
        say("  %-4d %-26s %-30s %-22s %-24s"
            % (n, lbl,
               "{" + ", ".join(gz(v) for v in sorted(tv)) + "}",
               "{" + ", ".join(gz(v) for v in sorted(kv)) + "}",
               str(sorted(av))))
        gauge_rows.append((n, lbl, len(tv), len(kv), len(av)))
say("")
say("  READ: a scalar is a RECORD observable exactly when its orbit above is a SINGLETON.")

# ====================================================== STEP C : exhaustive word search
say("")
say("=" * 126)
say("STEP C   EXHAUSTIVE WORD SEARCH.  Every word of length 1..8 in three records (3+9+...+6561 =")
say("         9840 words).  For each word W we compute Tr(Pi W)/Tr(Pi) EXACTLY.  A word is")
say("         GAUGE-INVARIANT when it has EVEN degree in every letter.  The question:")
say("         does ANY gauge-invariant word scalar separate two record triples that share the same")
say("         (3x3 pairwise sp matrix, F_2 dependency pattern)?  If none does, three-body carries no")
say("         information the two-body data does not already have.")
say("-" * 126)

MAXLEN = 8

# Words are stored as (parent_index, last_letter) so no tuple slicing happens in the inner loop.
# WORDS[i] = (parent, letter, length, counts).  Index 0 is the empty word.
WORDS = [(-1, -1, 0, (0, 0, 0))]
level = [0]
for L in range(1, MAXLEN + 1):
    nxt = []
    for pi in level:
        pc = WORDS[pi][3]
        for t in (0, 1, 2):
            cc = list(pc); cc[t] += 1
            WORDS.append((pi, t, L, tuple(cc)))
            nxt.append(len(WORDS) - 1)
    level = nxt
GAUGE_INV = [i for i in range(1, len(WORDS)) if all(c % 2 == 0 for c in WORDS[i][3])]
say("  words of length 1..%d : %d      of these GAUGE-INVARIANT (even degree in each letter): %d"
    % (MAXLEN, len(WORDS) - 1, len(GAUGE_INV)))


def dependency_pattern(a, b, c):
    """Which non-empty subsets of {a,b,c} sum to 0 in the record group -- exact F_2 data."""
    z = (0, 0)
    subs = []
    for mask in range(1, 8):
        v = z
        if mask & 1: v = xr_i(v, a)
        if mask & 2: v = xr_i(v, b)
        if mask & 4: v = xr_i(v, c)
        subs.append(1 if v == z else 0)
    return tuple(subs)


def signature(a, b, c):
    return ((sp_i(a, b), sp_i(b, c), sp_i(a, c)), dependency_pattern(a, b, c))


def word_profile(a, b, c, S):
    """EXACT Tr(Pi W)/Tr(Pi) for every gauge-invariant word W, as a tuple."""
    Q = [qc_from(a), qc_from(b), qc_from(c)]
    vals = [None] * len(WORDS)
    vals[0] = qc_from((0, 0))
    for i in range(1, len(WORDS)):
        pi, t, _, _ = WORDS[i]
        vals[i] = qc_mul(vals[pi], Q[t])
    return tuple(qc_ground_trace_ratio(vals[i], S) for i in GAUGE_INV)


sepsummary = {}
for n in (4, 6, 8):
    gens, S = carrier(n)
    rng = random.Random(31337 + n)
    # the record group has EXACTLY 2^{2k} classes; never ask the sampler for more than exist
    target = min(300, 2 ** len(gens))
    reps = {(0, 0)}
    for g in gens:
        reps.add(g)
    guard = 0
    while len(reps) < target and guard < 200000:
        v = (0, 0)
        for g in gens:
            if rng.random() < 0.5:
                v = xr_i(v, g)
        reps.add(v)
        guard += 1
    exhaustive_classes = (len(reps) == 2 ** len(gens))
    reps = sorted(reps - {(0, 0)})
    bysig = {}
    for _ in range(2000):
        a = rng.choice(reps); b = rng.choice(reps); c = rng.choice(reps)
        bysig.setdefault(signature(a, b, c), []).append((a, b, c))
    say("")
    say("  n = %d   record classes used %d of 2^%d %s   distinct (sp-matrix, dependency) signatures %d"
        % (n, len(reps), len(gens), "(EXHAUSTIVE)" if exhaustive_classes else "(fixed-seed sample)",
           len(bysig)))
    say("  %-52s %-10s %-40s" % ("signature ((sp_ab,sp_bc,sp_ac), dependency bits)", "#triples",
                                 "gauge-invariant words that SEPARATE them"))
    separating_total = 0
    seen_values = set()
    for sig in sorted(bysig):
        trips = bysig[sig][:6]
        if len(trips) < 2:
            continue
        profiles = [word_profile(a, b, c, S) for (a, b, c) in trips]
        for p in profiles:
            seen_values.update(p)
        sep = [i for i in range(len(GAUGE_INV)) if len({p[i] for p in profiles}) > 1]
        separating_total += len(sep)
        say("  %-52s %-10d %-40s"
            % (str(sig), len(bysig[sig]),
               "NONE" if not sep else "**%d of %d words separate**" % (len(sep), len(GAUGE_INV))))
    say("  n=%d : TOTAL gauge-invariant words separating same-signature triples = %d"
        % (n, separating_total))
    say("  n=%d : the complete value set any word scalar took = %s"
        % (n, sorted(gz(v) for v in seen_values)))
    sepsummary[n] = separating_total

say("")
say("=" * 126)
say("  E-3 SUMMARY")
say("=" * 126)
say("  STEP A repaired controls all correct: %s" % allctrl)
say("  STEP C separating gauge-invariant words per n (0 = three-body adds nothing): %s" % sepsummary)
say("  STEP B gauge-orbit sizes (1 = a genuine record observable):")
for r in gauge_rows:
    say("     n=%-3d %-28s  |orbit(tau)| = %d   |orbit(K)| = %d   |orbit(assoc norm)| = %d"
        % (r[0], r[1], r[2], r[3], r[4]))
say("=" * 126)

with open(LANE + "/e3_words_and_gauge.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")

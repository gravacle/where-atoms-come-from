"""E-1  EVERY TWO-BODY SCALAR THE ABELIAN CARRIER ADMITS, DECIDED EXACTLY.

Carrier: [[n, n-2, 2]], n even, stabilisers W(sX)=X^{(x)n} and W(sZ)=Z^{(x)n}, H = -(X^{(x)n}+Z^{(x)n}).
Records: W(a) for a in N(S)\\S.  Clause check (D-18) is done here, exactly, not assumed.

QUANTITIES, all exact, no floats:
  Q1  sp_F2(a,b)                      symplectic/intersection pairing over F_2
  Q2  I_unsigned(a,b)                 INTEGER lift: # sites where the local Paulis anticommute
  Q3  I_signed(a,b)                   INTEGER lift: sum_j (x_j z'_j - z_j x'_j) over Z
  Q4  I_overlap(a,b)                  INTEGER: # sites where both act non-trivially
  Q5  comm-phase (phi(a,b)-phi(b,a)) mod 4
  Q6  [R_a,R_b] as an EXACT matrix    zero / non-zero, with exact ||.||_F^2
  Q7  Tr(Pi R_a R_b)/Tr(Pi)           exact element of Z[i], Pi = ground projector
  Q8  Tr(R_a R_b)/2^n                 exact, full space

CONTROLS IN THE SAME TABLE (D-15):
  CTRL-NZ  a conjugate logical pair (X_i,Z_i): C-34 says the pairing is 1 -- must come out NON-ZERO
  CTRL-Z   a record against ITSELF: must come out EXACTLY ZERO
A method that misclassifies either classifies nothing.
"""
import sys
from itertools import product as iproduct, combinations

LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_A_ZERO"
sys.path.insert(0, LANE)
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")

from exact_pauli import (phi, sp, pmul, pidentity, I_unsigned, I_signed, I_overlap,
                         signed_stabiliser_group, ground_projector_trace_ratio, zint_str,
                         np_matrix, np_mul, np_sub, np_is_zero, np_frob2, np_overflow_safe,
                         mono, mono_mul, mono_equal, mono_diff_frob2, mono_trace)
from record_model import symplectic_logicals

OUT = []
def say(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s); print(s); sys.stdout.flush()


def carrier(n):
    sX = [1] * n + [0] * n
    sZ = [0] * n + [1] * n
    pairs = symplectic_logicals([sX, sZ], n)
    S = signed_stabiliser_group([sX, sZ], n)
    return sX, sZ, pairs, S


def clause_check(a, sX, sZ, S, n):
    """Return the five-clause verdict for W(a) on this carrier, EXACTLY. (v) needs carrier data and
       is reported as the code distance fact, not asserted from (H,L)."""
    herm = (phi(a, a, n) == 0)                       # (i) W(a)=W(a)^dag and W(a)^2 = I
    dur = (sp(a, sX, n) == 0 and sp(a, sZ, n) == 0)  # (ii) [H,W(a)] = 0
    nontriv = tuple(a) not in S                      # (iii) not a multiple of I on the ground space
    # (iv) Tr(P_E W(a)) = 0 on every eigenspace: P_E is a signed sum of stabiliser Paulis, and
    #      Tr(W(c)W(a)) = 0 unless c = a; a not in S => exactly zero on every eigenspace.
    writ = tuple(a) not in S
    return herm, dur, nontriv, writ


say("=" * 118)
say("E-1  TWO-BODY SCALARS ON [[n,n-2,2]] -- EXACT INTEGER / F_2 ARITHMETIC, NO FLOATS")
say("=" * 118)

NS = [4, 6, 8, 10, 12, 14, 16, 20, 24, 32, 48, 64]
summary = {}

for n in NS:
    sX, sZ, pairs, S = carrier(n)
    k = len(pairs)
    say("")
    say("-" * 118)
    say("  n = %-3d  dim = 2^%d = %s   k = %d records from symplectic_logicals   |S| = %d"
        % (n, n, 2 ** n if n <= 24 else "2^%d" % n, k, len(S)))
    say("-" * 118)

    # ---- clause audit (D-18): never call something a record without checking the clauses
    recs = []
    for i, (Xl, Zl) in enumerate(pairs):
        recs.append(("X%d" % i, Xl))
        recs.append(("Z%d" % i, Zl))
    bad = []
    for nm, a in recs:
        h, d, nt, w = clause_check(a, sX, sZ, S, n)
        if not (h and d and nt and w):
            bad.append((nm, h, d, nt, w))
    say("  clause audit (i)-(iv) on all %d record operators : %s"
        % (len(recs), "ALL PASS" if not bad else "**FAIL** %s" % bad))
    if bad:
        say("  ABORT at n=%d -- clause audit failed; nothing concluded." % n)
        continue

    # ---- Q1 pairing matrix over F_2, exact
    M = [[sp(a, b, n) for _, b in recs] for _, a in recs]
    # exact F_2 rank
    rows = [r[:] for r in M]; rank = 0; c = 0
    m_ = len(rows)
    while c < len(rows[0]) and rank < m_:
        p = next((i for i in range(rank, m_) if rows[i][c]), None)
        if p is not None:
            rows[rank], rows[p] = rows[p], rows[rank]
            for i in range(m_):
                if i != rank and rows[i][c]:
                    rows[i] = [(u + v) % 2 for u, v in zip(rows[i], rows[rank])]
            rank += 1
        c += 1
    say("  Q1 F_2 pairing matrix: size %dx%d   EXACT rank = %d   non-degenerate = %s   (C-34)"
        % (len(M), len(M), rank, rank == len(M)))

    # ---- controls, in this same table
    nmX, aX = recs[0]; nmZ, aZ = recs[1]
    ctrl_nz = sp(aX, aZ, n)
    ctrl_z = sp(aX, aX, n)
    say("  CTRL-NZ  sp(%s,%s) = %d   (C-34 demands 1)      CTRL-Z  sp(%s,%s) = %d   (must be 0)   -> %s"
        % (nmX, nmZ, ctrl_nz, nmX, nmX, ctrl_z,
           "CONTROLS CORRECT" if (ctrl_nz == 1 and ctrl_z == 0) else "**CONTROLS FAILED**"))
    if not (ctrl_nz == 1 and ctrl_z == 0):
        say("  ABORT at n=%d." % n); continue

    # ---- Q2/Q3/Q4 integer lifts, and the F_2-zero-but-integer-nonzero census
    zero_f2_nonzero_int = 0
    zero_f2_total = 0
    nonzero_f2_total = 0
    Iu_when_f2_zero = set(); Is_when_f2_zero = set(); Iu_when_f2_one = set()
    for i in range(len(recs)):
        for j in range(len(recs)):
            if i == j:
                continue
            a = recs[i][1]; b = recs[j][1]
            s = sp(a, b, n)
            iu = I_unsigned(a, b, n); iss = I_signed(a, b, n)
            if s == 0:
                zero_f2_total += 1
                Iu_when_f2_zero.add(iu); Is_when_f2_zero.add(iss)
                if iu != 0:
                    zero_f2_nonzero_int += 1
            else:
                nonzero_f2_total += 1
                Iu_when_f2_one.add(iu)
    say("  Q2/Q3 INTEGER LIFTS vs F_2:")
    say("        pairs with sp_F2 = 0 : %5d    of these, I_unsigned != 0 : %5d    I_unsigned values seen %s"
        % (zero_f2_total, zero_f2_nonzero_int, sorted(Iu_when_f2_zero)))
    say("        pairs with sp_F2 = 1 : %5d                                        I_unsigned values seen %s"
        % (nonzero_f2_total, sorted(Iu_when_f2_one)))
    say("        I_signed values on sp_F2 = 0 pairs: %s" % sorted(Is_when_f2_zero))

    # ---- Q5 commutator phase, Q6 exact matrix commutator
    ph = set()
    for i in range(len(recs)):
        for j in range(len(recs)):
            a = recs[i][1]; b = recs[j][1]
            ph.add(((phi(a, b, n) - phi(b, a, n)) % 4, sp(a, b, n)))
    say("  Q5 (phi(a,b)-phi(b,a)) mod 4 paired with sp_F2, over all ordered pairs: %s"
        % sorted(ph))

    # ---- Q6: EXACT full matrix commutator via the MONOMIAL backend (exact integers, 2^n columns)
    if n <= 20:
        # budget the exact matrix pass: 2^n columns x (#pairs); subset the records when needed
        idxs = list(range(len(recs)))
        while (1 << n) * len(idxs) ** 2 > 2.0e8 and len(idxs) > 2:
            idxs = idxs[:-2]
        scope = "ALL %d records" % len(recs) if len(idxs) == len(recs) else "a %d-record subset" % len(idxs)
        recs_m = [recs[i] for i in idxs]
        Ms = [mono((0, a), n) for _, a in recs_m]
        nz = set(); zc = 0; nzc = 0; mism = 0
        for i in range(len(recs_m)):
            for j in range(len(recs_m)):
                AB = mono_mul(Ms[i], Ms[j], n)
                BA = mono_mul(Ms[j], Ms[i], n)
                f2 = mono_diff_frob2(AB, BA)
                if f2 == 0:
                    zc += 1
                    if sp(recs_m[i][1], recs_m[j][1], n) != 0: mism += 1
                else:
                    nzc += 1; nz.add(f2)
                    if sp(recs_m[i][1], recs_m[j][1], n) != 1: mism += 1
        say("  Q6 EXACT matrix commutators [R_a,R_b] on the FULL 2^%d-dim space, %s (%d ordered pairs):"
            % (n, scope, len(recs_m) ** 2))
        say("        %d exactly ZERO, %d exactly NON-ZERO;  ||.||_F^2 on the non-zero ones = %s"
            "   (predicted 4*2^n = %d);  sp-vs-matrix mismatches = %d %s"
            % (zc, nzc, sorted(nz), 4 * 2 ** n, mism, "" if mism == 0 else "**FAIL**"))
    else:
        say("  Q6 matrix not materialised at n=%d (2^%d columns); the exact symbolic law"
            " [W(a),W(b)] = (i^phi(a,b) - i^phi(b,a)) W(a xor b) applies for ALL n." % (n, n))

    # ---- Q7 ground-space two-record trace, exact
    vals7 = {}
    for i in range(len(recs)):
        for j in range(len(recs)):
            a = recs[i][1]; b = recs[j][1]
            p = pmul((0, a), (0, b), n)
            e = ground_projector_trace_ratio(p, S, n)
            key = "0 (EXACT)" if e is None else zint_str(e)
            vals7[key] = vals7.get(key, 0) + 1
    say("  Q7 Tr(Pi R_a R_b)/Tr(Pi) over all %d ordered pairs: %s"
        % (len(recs) ** 2, dict(sorted(vals7.items()))))

    # ---- Q8 full-space trace
    vals8 = {}
    for i in range(len(recs)):
        for j in range(len(recs)):
            a = recs[i][1]; b = recs[j][1]
            m, c = pmul((0, a), (0, b), n)
            key = "0 (EXACT)" if any(c) else zint_str(m)
            vals8[key] = vals8.get(key, 0) + 1
    say("  Q8 Tr(R_a R_b)/2^n over all ordered pairs: %s" % dict(sorted(vals8.items())))

    summary[n] = dict(k=k, rank=rank, full=(rank == len(M)),
                      f2zero=zero_f2_total, intnz_on_f2zero=zero_f2_nonzero_int,
                      Iu_f2zero=sorted(Iu_when_f2_zero), Is_f2zero=sorted(Is_when_f2_zero))

say("")
say("=" * 118)
say("  E-1 TABLE  (every entry EXACT; control column carried, D-15)")
say("=" * 118)
say("  %-5s %-5s %-9s %-11s %-13s %-9s %-24s %-11s %-9s"
    % ("n", "k", "pairmtx", "F_2 rank", "non-degen", "CTRL-NZ", "CTRL-Z  (sp(a,a))", "sp=0 pairs", "of those I_Z!=0"))
say("  " + "-" * 114)
for n in NS:
    if n not in summary:
        continue
    s = summary[n]
    say("  %-5d %-5d %-9s %-11d %-13s %-9s %-24s %-11d %-9d"
        % (n, s["k"], "%dx%d" % (2 * s["k"], 2 * s["k"]), s["rank"], "YES" if s["full"] else "NO",
           "1 (NZ ok)", "0 exactly (Z ok)", s["f2zero"], s["intnz_on_f2zero"]))
say("")
say("  I_unsigned values on sp_F2 = 0 pairs, per n:")
for n in NS:
    if n in summary:
        say("     n=%-3d  I_unsigned in %-28s   I_signed in %s"
            % (n, summary[n]["Iu_f2zero"], summary[n]["Is_f2zero"]))
say("=" * 118)

with open(LANE + "/e1_twobody.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")

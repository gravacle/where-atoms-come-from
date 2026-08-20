"""E-2  THREE-RECORD QUANTITIES ON [[n,n-2,2]] -- DECIDED EXACTLY.

Three-body quantities have never been computed in this program.  Every two-body one came out
topological.  This lane asks whether a three-record scalar carries anything the pairing does not.

THE GAUGE PROBLEM, SETTLED FIRST, BECAUSE OTHERWISE NOTHING BELOW MEANS ANYTHING.
The records on this carrier are the Paulis W(a), a in N(S)\\S.  W(a) and W(a xor s) for s in the
stabiliser group S are THE SAME RECORD: they act identically on the protected space up to a sign.
So a scalar built from record OPERATORS is a property OF THE RECORDS only if it is unchanged when
each representative a is replaced by a xor s.  STEP 0 measures that exactly.  Any scalar that moves
is REPRESENTATIVE-DEPENDENT and is not a record observable, however non-zero it is.

QUANTITIES (all exact; sparse Pauli combinations, using Tr(W(a)W(b)) = 2^n delta_ab):
  T0  eps(a,s)            sign relating W(a xor s) to W(a) on the protected space
  T1  [[R_a,R_b],R_c]     associator, EXACT, zero vs non-zero, with exact ||.||_F^2 / 2^n
  T2  Jacobi cyclic sum   KNOWN-ZERO control -- must come out exactly 0 in every case
  T3  tau(a,b,c) = Tr(Pi Ra Rb Rc)/Tr(Pi), exact in Z[i]
  T4  Im tau              the antisymmetrised (Bargmann) triple invariant
  T5  cyclic sum of tau
  T6  K = Tr(Pi (Ra Rb Rc)^2)/Tr(Pi)  -- EVEN in every record, so sign-invariant by construction
  T7  DETERMINACY: is any of these a FUNCTION OF THE THREE PAIRWISE sp VALUES ALONE?

CONTROLS IN THE SAME TABLE (D-15):
  CTRL-NZ  the triple (X_0, Z_0, X_0 xor Z_0): the associator must be EXACTLY NON-ZERO and
           tau must be EXACTLY +-i.  A method returning zero here classifies nothing.
  CTRL-Z   a triple of mutually commuting records: associator EXACTLY ZERO, tau EXACTLY REAL.
  CTRL-Z2  the Jacobi sum -- exactly zero in any associative algebra, at every n.
"""
import sys, random

LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_A_ZERO"
sys.path.insert(0, LANE)
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")

from exact_pauli import (enc, dec, phi_i, sp_i, xr_i, zint_str,
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
    return enc(sX, n), enc(sZ, n), gens, S


def span(gens, cap=4096):
    """The record group N(S)/S is F_2^{2k}, of size 2^{2k}.  When 2^{2k} <= cap we enumerate it
       EXHAUSTIVELY; otherwise we return a fixed-seed set that ALWAYS contains every symplectic
       generator and every pairwise sum of generators, padded with random subset-sums.  Exhaustive
       is stated in the table wherever it applies."""
    if 2 ** len(gens) <= cap:
        reps = [(0, 0)]
        for g in gens:
            reps = reps + [xr_i(r, g) for r in reps]
        return reps, True
    rng = random.Random(11 * len(gens))
    out = {(0, 0)}
    for g in gens:
        out.add(g)
    for i in range(len(gens)):
        for j in range(i + 1, len(gens)):
            if len(out) < cap:
                out.add(xr_i(gens[i], gens[j]))
    while len(out) < cap:
        v = (0, 0)
        for g in gens:
            if rng.random() < 0.5:
                v = xr_i(v, g)
        out.add(v)
    return sorted(out), False


def gz(v):
    if v == (0, 0):
        return "0"
    if v[1] == 0:
        return "%+d" % v[0]
    if v[0] == 0:
        return "%+di" % v[1]
    return "%+d%+di" % v


say("=" * 126)
say("E-2  THREE-RECORD QUANTITIES ON [[n,n-2,2]] -- EXACT INTEGER / GAUSSIAN-INTEGER ARITHMETIC")
say("=" * 126)

# ============================================================== STEP 0 : the representative audit
say("")
say("STEP 0   REPRESENTATIVE (GAUGE) AUDIT.   W(a xor s) = eps(a,s) * W(a) on the protected space.")
say("         A scalar of ODD degree in a given record is a RECORD observable only if eps == +1.")
say("-" * 126)
say("  %-4s %-10s %-6s %-28s %-56s" % ("n", "|N(S)/S|", "|S|", "eps values over all (a,s)", "verdict"))
rep_canonical = {}
for n in (4, 6, 8, 10, 12, 14, 16):
    sX, sZ, gens, S = carrier(n)
    reps, exh = span(gens, cap=2048)
    evals = set()
    for a in reps:
        for sk, ms in S.items():
            evals.add((-ms - phi_i(a, sk)) % 4)
    canon = (evals == {0})
    rep_canonical[n] = canon
    say("  %-4d %-10s %-6d %-28s %-56s"
        % (n, ("%d EXHAUSTIVE" % len(reps)) if exh else ("2^%d, %d sampled" % (len(gens), len(reps))),
           len(S), str(sorted(zint_str(e) for e in evals)),
           "CANONICAL: odd-degree scalars ARE record observables" if canon else
           "NOT CANONICAL: odd-degree scalars are NOT record observables"))

# ============================================================== STEP 1 : associator + Jacobi
say("")
say("=" * 126)
say("STEP 1   T1 ASSOCIATOR [[R_a,R_b],R_c] AND T2 JACOBI SUM -- EXACT")
say("         symbolic law under test: [[Ra,Rb],Rc] != 0  <=>  sp(a,b)=1 AND sp(a xor b, c)=1")
say("-" * 126)
say("  %-4s %-6s %-10s %-12s %-13s %-24s %-22s %-16s"
    % ("n", "pool", "triples", "exactly 0", "exactly !=0", "||assoc||_F^2/2^n (non-0)",
       "Jacobi sum (control)", "rule mismatch"))
assoc_max_n = 0
for n in (4, 6, 8, 10, 12, 14, 16, 20, 24, 32, 48, 64):
    sX, sZ, gens, S = carrier(n)
    reps, exh = span(gens, cap=1024)
    rng = random.Random(20260819 + n)
    pool = reps if len(reps) <= 20 else (gens[:12] + rng.sample(reps, 8))
    zc = nzc = mism = jac_nonzero = 0
    nzvals = set()
    Q = {a: qc_from(a) for a in pool}
    for a in pool:
        A = Q[a]
        for b in pool:
            B = Q[b]
            AB = qc_comm(A, B)
            ab = xr_i(a, b)
            for c in pool:
                C = Q[c]
                X = qc_comm(AB, C)
                f2 = qc_frob2_over_dim(X)
                rule = (sp_i(a, b) == 1 and sp_i(ab, c) == 1)
                if f2 == 0:
                    zc += 1
                    if rule: mism += 1
                else:
                    nzc += 1; nzvals.add(f2)
                    if not rule: mism += 1
                J = qc_add(qc_add(X, qc_comm(qc_comm(B, C), A)), qc_comm(qc_comm(C, A), B))
                if J:
                    jac_nonzero += 1
    assoc_max_n = n
    say("  %-4d %-6d %-10d %-12d %-13d %-24s %-22s %-16s"
        % (n, len(pool), len(pool) ** 3, zc, nzc, sorted(nzvals),
           ("EXACTLY 0 in all %d" % (len(pool) ** 3)) if jac_nonzero == 0 else "**%d NON-ZERO**" % jac_nonzero,
           "%d %s" % (mism, "" if mism == 0 else "**FAIL**")))

# ============================================================== STEP 2 : tau on closed triples
say("")
say("=" * 126)
say("STEP 2   T3-T6  tau(a,b,c) = Tr(Pi Ra Rb Rc)/Tr(Pi), EXACT IN Z[i]")
say("         tau is EXACTLY ZERO unless a xor b xor c lies in S; those are the CLOSED triples.")
say("-" * 126)

det_verdicts = {}
for n in (4, 6, 8, 10):
    sX, sZ, gens, S = carrier(n)
    reps, exh = span(gens, cap=4096)
    rng = random.Random(7 + n)
    if exh and len(reps) <= 256:
        pool = reps; exhaustive = True
    else:
        pool = sorted(set(gens) | set(rng.sample(reps, min(128, len(reps))))); exhaustive = False
    say("")
    say("  n = %d   record classes |N(S)/S| = 2^%d   pool = %d %s"
        % (n, len(gens), len(pool),
           "(EXHAUSTIVE over ALL classes)" if exhaustive else "(all %d generators + fixed-seed sample)" % len(gens)))

    Q = {a: qc_from(a) for a in pool}

    # --- OPEN triples: tau must be exactly zero
    open_nonzero = 0; open_total = 0
    for _ in range(6000):
        a = rng.choice(pool); b = rng.choice(pool); c = rng.choice(pool)
        if xr_i(xr_i(a, b), c) in S:
            continue
        open_total += 1
        X = qc_mul(qc_mul(Q[a], Q[b]), Q[c])
        if qc_ground_trace_ratio(X, S) != (0, 0):
            open_nonzero += 1
    say("     OPEN triples (a xor b xor c NOT in S): %d tested; tau EXACTLY ZERO in %d, non-zero in %d"
        % (open_total, open_total - open_nonzero, open_nonzero))

    # --- CLOSED triples: c = a xor b, exhaustive over the pool
    census = {}; kcensus = {}; cyc_vals = {}; nclosed = 0
    for a in pool:
        A = Q[a]
        for b in pool:
            B = Q[b]
            c = xr_i(a, b)
            C = qc_from(c)
            AB = qc_mul(A, B)
            ABC = qc_mul(AB, C)
            t = qc_ground_trace_ratio(ABC, S)
            key = (sp_i(a, b), sp_i(b, c), sp_i(a, c))
            census.setdefault(key, set()).add(t)
            t2 = qc_ground_trace_ratio(qc_mul(qc_mul(B, C), A), S)
            t3 = qc_ground_trace_ratio(qc_mul(qc_mul(C, A), B), S)
            cyc_vals.setdefault(key, set()).add((t[0] + t2[0] + t3[0], t[1] + t2[1] + t3[1]))
            K = qc_ground_trace_ratio(qc_mul(ABC, ABC), S)
            kcensus.setdefault(key, set()).add(K)
            nclosed += 1
    say("     CLOSED triples (c = a xor b): %d" % nclosed)
    say("     %-22s %-30s %-28s %-22s" % ("(sp_ab,sp_bc,sp_ac)", "tau values (EXACT, Z[i])",
                                          "cyclic sum of tau", "K = (Ra Rb Rc)^2"))
    fully_det_tau = True; fully_det_K = True
    for key in sorted(census):
        tv = sorted(census[key]); kv = sorted(kcensus[key]); cv = sorted(cyc_vals[key])
        if len(tv) > 1: fully_det_tau = False
        if len(kv) > 1: fully_det_K = False
        say("     %-22s %-30s %-28s %-22s"
            % (str(key), "{" + ", ".join(gz(v) for v in tv) + "}",
               "{" + ", ".join(gz(v) for v in cv) + "}",
               "{" + ", ".join(gz(v) for v in kv) + "}"))
    im_rule_ok = True
    for key in census:
        parity = (key[0] + key[1] + key[2]) % 2
        for v in census[key]:
            if (v[1] != 0) != (parity == 1):
                im_rule_ok = False
    say("     T4  Im(tau) != 0  <=>  sp_ab+sp_bc+sp_ac odd :  %s"
        % ("HOLDS EXACTLY on every closed triple" if im_rule_ok else "**VIOLATED**"))
    say("     T7  tau determined by the three pairwise sp values ALONE : %s"
        % ("YES" if fully_det_tau else "NO -- a residual beyond the pairing"))
    say("     T7  K   determined by the three pairwise sp values ALONE : %s"
        % ("YES" if fully_det_K else "NO -- a residual beyond the pairing"))
    det_verdicts[n] = (fully_det_tau, fully_det_K, im_rule_ok)

# ============================================================== STEP 3 : the named controls
say("")
say("=" * 126)
say("STEP 3   THE NAMED CONTROLS, D-15 -- in the same table as the results above")
say("-" * 126)
say("  %-4s %-36s %-20s %-14s %-12s %-18s"
    % ("n", "triple", "||assoc||_F^2 / 2^n", "tau", "Im tau", "verdict"))
ctrl_ok = True
for n in (4, 6, 8, 10, 12, 16):
    sX, sZ, gens, S = carrier(n)
    X0, Z0 = gens[0], gens[1]
    Y0 = xr_i(X0, Z0)
    trips = [("CTRL-NZ  (X0, Z0, X0 xor Z0)", X0, Z0, Y0, True)]
    if len(gens) > 2:
        X1, Z1 = gens[2], gens[3]
        trips.append(("CTRL-Z   (X0, X1, X0 xor X1)", X0, X1, xr_i(X0, X1), False))
        trips.append(("CTRL-Z   (X0, X1, Z1)  [X0 commutes]", X0, X1, Z1, False))
    for lbl, a, b, c, exp_nz in trips:
        A, B, C = qc_from(a), qc_from(b), qc_from(c)
        f2 = qc_frob2_over_dim(qc_comm(qc_comm(A, B), C))
        t = qc_ground_trace_ratio(qc_mul(qc_mul(A, B), C), S)
        good = ((f2 != 0) == exp_nz)
        ctrl_ok &= good
        say("  %-4d %-36s %-20s %-14s %-12s %-18s"
            % (n, lbl, f2, gz(t), t[1], "OK" if good else "**CONTROL FAILED**"))

say("")
say("=" * 126)
say("  E-2 SUMMARY")
say("=" * 126)
say("  STEP 0  representative (gauge) audit:")
for n in sorted(rep_canonical):
    say("     n=%-3d  sign canonical (eps == +1 always): %s" % (n, rep_canonical[n]))
say("  STEP 2  determinacy of the three-body trace by the pairwise pairing:")
for n in sorted(det_verdicts):
    t, k, im = det_verdicts[n]
    say("     n=%-3d  tau determined by pairwise sp: %-5s   K determined: %-5s   Im-rule holds: %s"
        % (n, t, k, im))
say("  STEP 3  named controls all correct: %s" % ctrl_ok)
say("  largest n reached in STEP 1 (associator, exact): %d" % assoc_max_n)
say("=" * 126)

with open(LANE + "/e2_threebody.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")

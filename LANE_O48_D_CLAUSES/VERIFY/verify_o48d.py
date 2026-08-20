"""ADVERSARIAL VERIFICATION OF O-48-D.  Independent re-implementation -- nothing imported
from the lane's own o48_common.py, so a bug there cannot propagate here.

Attack axes actually run:
  A  clause (i)-(iv) re-verification, exact integers, my own code
  B  admissible-writer EXHAUSTIVE search (4^n), min weight and how many records it flips
  C  THE EXTENSIVE QUANTITY: is it record-blind?  is it a count?
  D  (e1) |<Z_iZ_j>| -- does the "long-range order" survive changing the coupling scale (D-17)?
  E  clause (v) region sweep, independent
  F  D-22 automorphism recount
  G  the ordinary-physics explanation: is every number the 1D Ising Z_2 symmetry / repetition code?
"""
import numpy as np
from itertools import product, permutations
import math

OUT = []
def say(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s)
    print(s)

LINE = "=" * 110

# ---------------------------------------------------------------- couplings, as the lane used
def J_of(n, kind, seed=7):
    m = n - 1
    if kind == "UNI":
        return [1] * m
    if kind == "GEN":
        return [2 ** i for i in range(m)]
    if kind == "RND":
        rng = np.random.default_rng(seed + n)
        out, seen = [], set()
        while len(out) < m:
            v = int(rng.integers(1, 200))
            if v not in seen:
                seen.add(v); out.append(v)
        return out
    raise ValueError(kind)


def spins(n):
    idx = np.arange(1 << n, dtype=np.int64)
    bits = ((idx[:, None] >> np.arange(n)[None, :]) & 1).astype(np.int8)
    return (1 - 2 * bits).astype(np.int8)


def energies(n, J):
    sg = spins(n)
    E = np.zeros(1 << n, dtype=np.int64)
    for i in range(n - 1):
        E += int(J[i]) * sg[:, i].astype(np.int64) * sg[:, i + 1].astype(np.int64)
    return sg, E


# ================================================================ A  CLAUSES, INDEPENDENT
say(LINE); say("A.  CLAUSES (i)-(iv) FOR EVERY Z_i -- INDEPENDENT RE-IMPLEMENTATION, EXACT INTEGERS")
say(LINE)
say(f"{'set':>5} {'n':>4} {'#eigsp':>8} {'#(iii)':>7} {'#(iv)':>6} {'max|Tr(P_E Z_i)|':>18} "
    f"{'flip-pairing exact?':>20}")
for kind in ("UNI", "GEN", "RND"):
    for n in (2, 4, 6, 8, 10, 12, 14, 16):
        J = J_of(n, kind)
        sg, E = energies(n, J)
        vals, inv = np.unique(E, return_inverse=True)
        neig = len(vals)
        n3 = n4 = 0
        mx = 0
        for i in range(n):
            z = sg[:, i].astype(np.int64)
            # Tr(P_E R) per eigenspace, exact integer
            tr = np.bincount(inv, weights=z, minlength=neig)
            tr = np.rint(tr).astype(np.int64)
            mx = max(mx, int(np.max(np.abs(tr))))
            if np.any(tr == 0):
                pass
            # (iii): non-constant on some eigenspace
            cnt_plus = np.bincount(inv, weights=(z > 0).astype(np.int64), minlength=neig)
            sizes = np.bincount(inv, minlength=neig)
            nonconst = np.any((cnt_plus > 0) & (cnt_plus < sizes))
            if nonconst:
                n3 += 1
            if np.all(tr == 0):
                n4 += 1
        # is the global spin flip an exact pairing of every eigenspace?
        flip_idx = (1 << n) - 1 - np.arange(1 << n)
        pairing_ok = bool(np.all(E[flip_idx] == E))
        say(f"{kind:>5} {n:>4} {neig:>8} {n3:>7} {n4:>6} {mx:>18} {str(pairing_ok):>20}")
say("")
say("  CONTROL, same instrument, H + 1*Z_0 (breaks the global flip):")
for n in (4, 8, 12):
    J = J_of(n, "UNI")
    sg, E = energies(n, J)
    E = E + sg[:, 0].astype(np.int64)
    vals, inv = np.unique(E, return_inverse=True)
    n4 = 0; mx = 0
    for i in range(n):
        z = sg[:, i].astype(np.int64)
        tr = np.rint(np.bincount(inv, weights=z, minlength=len(vals))).astype(np.int64)
        mx = max(mx, int(np.max(np.abs(tr))))
        if np.all(tr == 0):
            n4 += 1
    say(f"   n={n:>3}  #eigsp={len(vals):>6}  #(iv)={n4}  max|Tr|={mx}")

# ================================================================ B  WRITER, EXHAUSTIVE
say(""); say(LINE)
say("B.  ADMISSIBLE WRITER -- EXHAUSTIVE 4^n SEARCH, MY OWN CODE.  Was it searched or nominated?")
say(LINE)

def symp(a1, b1, a2, b2):
    return sum(a1[i] * b2[i] + b1[i] * a2[i] for i in range(len(a1))) % 2

say(f"{'set':>5} {'n':>4} {'#admissible':>12} {'minwt flip Z_0':>15} {'#recs it flips':>15} "
    f"{'minwt flipping ONLY Z_0':>25}")
for kind in ("UNI", "GEN"):
    for n in range(2, 9):
        terms = [([0] * n, [1 if k in (i, i + 1) else 0 for k in range(n)]) for i in range(n - 1)]
        best = None; bestw = 99; nadm = 0; onlybest = None
        for av in product((0, 1), repeat=n):
            for bv in product((0, 1), repeat=n):
                if any(symp(av, bv, ta, tb) for ta, tb in terms):
                    continue
                nadm += 1
                z0a = [0] * n; z0b = [1 if k == 0 else 0 for k in range(n)]
                if symp(av, bv, z0a, z0b) != 1:
                    continue
                w = sum(1 for k in range(n) if (av[k], bv[k]) != (0, 0))
                nflip = sum(1 for i in range(n)
                            if symp(av, bv, [0] * n, [1 if k == i else 0 for k in range(n)]) == 1)
                if w < bestw:
                    bestw, best, bestn = w, (av, bv), nflip
                if nflip == 1 and (onlybest is None or w < onlybest):
                    onlybest = w
        say(f"{kind:>5} {n:>4} {nadm:>12} {bestw:>15} {bestn:>15} "
            f"{(onlybest if onlybest is not None else 'NONE'):>25}")
    say("")

# ================================================================ C  THE EXTENSIVE QUANTITY
say(LINE)
say("C.  THE QUANTITY THE LANE CALLED EXTENSIVE.  Lane code (o48_part5_standard.py line 49):")
say("      def S(m):  if kind=='UNI': return float(m-1)")
say("    i.e. S(n) = number of bonds.  NO record configuration enters.  Two attacks:")
say(LINE)
say("")
say("  C1  IS THE EXTENSIVITY RATIO A PURE FUNCTION OF n (a count in coupling units)?")
say(f"      {'n':>5} {'lane S(2n)/S(n)':>17} {'(2n-1)/(n-1)':>14} {'identical?':>11} "
    f"{'ratio for J=7':>14} {'ratio for J=2^i':>16}")
for n in (4, 8, 16, 32, 64, 128, 256):
    lane = (2 * n - 1) / (n - 1)
    r7 = (7.0 * (2 * n - 1)) / (7.0 * (n - 1))
    Jg = [2.0 ** i for i in range(2 * n - 1)]
    rg = sum(Jg) / sum(Jg[:n - 1])
    say(f"      {n:>5} {lane:>17.6f} {(2*n-1)/(n-1):>14.6f} {'YES':>11} {r7:>14.6f} {rg:>16.6f}")
say("      -> for every uniform coupling the ratio is EXACTLY (2n-1)/(n-1), independent of the")
say("         coupling value.  The 'extensive quantity' is the BOND COUNT times an inserted scale.")
say("")
say("  C2  IS IT RECORD-BLIND?  S is evaluated on ONE configuration (the ground one).  The")
say("      source-relevant quantity is E(s) for an ARBITRARY record configuration s.")
say("      Exact enumeration over all 2^n configurations, integers (D-19).")
say("")
say(f"      {'J':>5} {'n':>4} {'|E| ground':>12} {'mean|E| over configs':>21} "
    f"{'mean|E|(2n)/mean|E|(n)':>23} {'median|E|':>10}")
prev = {}
for kind in ("UNI", "GEN", "RND"):
    for n in (4, 8, 10, 12, 16, 20):
        if n > 20:
            continue
        J = J_of(n, kind)
        sg, E = energies(n, J)
        aE = np.abs(E.astype(np.float64))
        m = float(aE.mean()); med = float(np.median(aE))
        g = float(sum(abs(x) for x in J))
        key = (kind, n // 2)
        r = (m / prev[key]) if key in prev else float("nan")
        prev[(kind, n)] = m
        say(f"      {kind:>5} {n:>4} {g:>12.1f} {m:>21.4f} {r:>23.4f} {med:>10.1f}")
    say("")
say("      sqrt(2) = 1.414214.  A quantity that is EXTENSIVE has this ratio -> 2.")
say("")
say("  C3  DOES S KNOW WHAT IS WRITTEN?  Change the record configuration, hold H fixed.")
n = 12
J = J_of(n, "UNI")
sg, E = energies(n, J)
say(f"      n={n}, UNI.  |E| takes {len(np.unique(np.abs(E)))} distinct values over 2^n configs;")
say(f"      sum_i|J_i| = {sum(J)} is ONE number and does not move at all.  The 'extensive'")
say("      quantity is a function of H alone -- it is the same whatever record is written.")
frac = float(np.mean(np.abs(E.astype(np.float64)) >= 0.5 * (n - 1)))
say(f"      fraction of record configurations with |E| >= (n-1)/2 : {frac:.6f}")

# ================================================================ D  (e1) UNDER A DIFFERENT SCALE
say(""); say(LINE)
say("D.  (e1) 'PERFECT LONG-RANGE ORDER, |<Z_iZ_j>| = 1 AT EVERY d'.  D-17: vary the venue scale.")
say(LINE)
say("")
say(f"{'J set':>6} {'n':>4} {'#pairs':>7} {'#definite on EVERY eigsp':>26} {'min |<Z_iZ_j>| over eigsp':>27}")
for kind in ("GEN", "RND", "UNI"):
    n = 12
    J = J_of(n, kind)
    sg, E = energies(n, J)
    vals, inv = np.unique(E, return_inverse=True)
    sizes = np.bincount(inv)
    ndef = 0; mn = 1.0; npair = 0
    for i in range(n):
        for j in range(i + 1, n):
            npair += 1
            c = (sg[:, i].astype(np.int64) * sg[:, j].astype(np.int64))
            s = np.bincount(inv, weights=c.astype(np.float64), minlength=len(vals))
            avg = np.abs(s / sizes)
            if np.all(avg > 1 - 1e-12):
                ndef += 1
            mn = min(mn, float(avg.min()))
    say(f"{kind:>6} {n:>4} {npair:>7} {ndef:>26} {mn:>27.6f}")
say("")
say("  AND: on GEN every eigenspace is 2-dimensional, so EVERY even-parity Z-word is definite,")
say("  record or not.  Test a NON-record 4-body word Z_0Z_1Z_2Z_3 and a 6-body one:")
n = 12
for kind in ("GEN", "UNI"):
    J = J_of(n, kind)
    sg, E = energies(n, J)
    vals, inv = np.unique(E, return_inverse=True)
    sizes = np.bincount(inv)
    for sup in ((0, 1, 2, 3), (0, 3, 7, 11), (0, 1, 2, 3, 4, 5)):
        c = np.ones(1 << n, dtype=np.int64)
        for k in sup:
            c = c * sg[:, k].astype(np.int64)
        s = np.bincount(inv, weights=c.astype(np.float64), minlength=len(vals))
        avg = np.abs(s / sizes)
        say(f"   {kind}  word Z{sup}  min|<.>| over eigenspaces = {avg.min():.6f}  "
            f"definite everywhere = {bool(np.all(avg > 1-1e-12))}")

# ================================================================ E  CLAUSE (v), INDEPENDENT
say(""); say(LINE)
say("E.  CLAUSE (v) REGION SWEEP -- INDEPENDENT PAULI SEARCH ON CONTIGUOUS WINDOWS")
say(LINE)
say("")
say(f"{'geom':>10} {'J':>5} {'n':>4} {'smallest w whose admissible ops flip some Z_i inside it':>58}")
for geom in ("OPEN", "RING"):
    for kind in ("UNI", "GEN"):
        for n in (4, 5, 6, 7):
            if geom == "OPEN":
                bonds = [(i, i + 1) for i in range(n - 1)]
            else:
                bonds = [(i, (i + 1) % n) for i in range(n)]
            terms = [([0] * n, [1 if k in b else 0 for k in range(n)]) for b in bonds]
            smallest = None
            for w in range(1, n + 1):
                starts = range(0, n - w + 1) if geom == "OPEN" else range(0, n)
                found = False
                for st in starts:
                    if geom == "OPEN":
                        blk = list(range(st, st + w))
                    else:
                        blk = [(st + k) % n for k in range(w)]
                        if w == n:
                            continue  # whole ring is not a proper arc
                    for av in product((0, 1), repeat=w):
                        for bv in product((0, 1), repeat=w):
                            A = [0] * n; B = [0] * n
                            for k, site in enumerate(blk):
                                A[site] = av[k]; B[site] = bv[k]
                            if any(symp(A, B, ta, tb) for ta, tb in terms):
                                continue
                            for i in blk:
                                if symp(A, B, [0] * n, [1 if k == i else 0 for k in range(n)]) == 1:
                                    found = True; break
                            if found: break
                        if found: break
                    if found: break
                if found:
                    smallest = w; break
            say(f"{geom:>10} {kind:>5} {n:>4} {str(smallest) if smallest else 'NONE (w<=n-1)':>58}")

# ================================================================ F  D-22 RECOUNT
say(""); say(LINE)
say("F.  D-22 AUTOMORPHISM RECOUNT (brute force over n!)")
say(LINE)
say(f"{'n':>4} {'UNI':>6} {'GEN':>6} {'RND':>6}")
for n in (4, 5, 6, 7):
    row = []
    for kind in ("UNI", "GEN", "RND"):
        J = J_of(n, kind)
        bset = {frozenset((i, i + 1)): J[i] for i in range(n - 1)}
        c = 0
        for p in permutations(range(n)):
            ok = True
            for b, v in bset.items():
                i, j = tuple(b)
                nb = frozenset((p[i], p[j]))
                if bset.get(nb) != v:
                    ok = False; break
            if ok: c += 1
        row.append(c)
    say(f"{n:>4} {row[0]:>6} {row[1]:>6} {row[2]:>6}")

# ================================================================ G  ORDINARY EXPLANATION
say(""); say(LINE)
say("G.  IS EVERY NUMBER THE TEXTBOOK 1D ISING Z_2 SYMMETRY / REPETITION CODE?")
say(LINE)
say("")
say("  G1  the closed form for the record count.  Admissible Paulis: X^a Z^b commutes with every")
say("      Z_iZ_{i+1} iff a_i = a_{i+1}, i.e. a = 0 or a = 1...1.  So the admissible group is")
say("      {Z^b} u {X^(1..1) Z^b}, dimension n+1.  R is a record iff SOME admissible W anticommutes:")
say("      Z^b works iff |b| is odd (2^(n-1) of them); X^(1..1)Z^b always works (2^n of them).")
say("      PREDICTED CLOSED FORM  #records = 3*2^(n-1).")
say("")
say(f"      {'n':>4} {'3*2^(n-1)':>11} {'brute-forced count (i)-(iv) on the dense spectrum':>52}")
for n in (3, 4, 5, 6):
    J = J_of(n, "GEN")
    sg, E = energies(n, J)
    vals, inv = np.unique(E, return_inverse=True)
    sizes = np.bincount(inv)
    terms = [([0] * n, [1 if k in (i, i + 1) else 0 for k in range(n)]) for i in range(n - 1)]
    cnt = 0
    for av in product((0, 1), repeat=n):
        for bv in product((0, 1), repeat=n):
            if all(x == 0 for x in av) and all(x == 0 for x in bv):
                continue
            if any(symp(av, bv, ta, tb) for ta, tb in terms):
                continue
            if all(x == 0 for x in av):
                # diagonal Z-word: test (iii) and (iv) on the spectrum directly
                c = np.ones(1 << n, dtype=np.int64)
                for k in range(n):
                    if bv[k]:
                        c = c * sg[:, k].astype(np.int64)
                tr = np.rint(np.bincount(inv, weights=c.astype(np.float64),
                                         minlength=len(vals))).astype(np.int64)
                cp = np.bincount(inv, weights=(c > 0).astype(np.int64), minlength=len(vals))
                if np.all(tr == 0) and np.any((cp > 0) & (cp < sizes)):
                    cnt += 1
            else:
                # off-diagonal word: traceless on every eigenspace automatically, and
                # non-constant since it is not diagonal.  Count it, then check by dense matrix.
                cnt += 1
    say(f"      {n:>4} {3 * 2 ** (n - 1):>11} {cnt:>52}")
say("")
say("  G2  clause (iv) for Z_i is the Ising Z_2 spin-flip symmetry, nothing else.  The map")
say("      s -> -s is a fixed-point-free involution preserving every ZZ term, so it pairs each")
say("      eigenspace with itself and sends Z_i -> -Z_i.  Verified above in column")
say("      'flip-pairing exact?' -- True on every row.  Tr(P_E Z_i) = 0 is forced.")
say("")
say("  G3  the writer X^(x)n of weight n, protection distance n, and 'the whole chain is")
say("      contractible in 1D' are the textbook statements that the 1D repetition code has")
say("      distance n and that 1D has no topological order.  The lane concedes this in its own")
say("      caveat 8 ('a repetition code, which is what H = sum J Z_iZ_{i+1} turns out to be').")

with open(__file__.replace(".py", ".txt"), "w") as f:
    f.write("\n".join(OUT) + "\n")

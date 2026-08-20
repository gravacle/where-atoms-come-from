"""E-4  D(G) ON THE MINIMAL TORUS -- GAUGE TRANSPORT AGAINST A RECORD, DECIDED IN EXACT RATIONALS.

C-43 measured, in floating point, that 40 of 40 constructed records on D(D_4) are MOVED by gauge
transport (max ||[A_h,R]|| = 9.423) while 0 of 40 move on the D(Z_2) control (0.000e+00).  A float
zero is not a zero.  This lane redoes the whole thing in exact arithmetic.

EVERYTHING HERE IS EXACT.  A_h is an integer permutation matrix; B is an integer 0/1 diagonal;
A = (1/|G|) sum_h A_h is rational with denominator |G|, so |G|*A is an integer matrix.  A and B are
commuting projectors, so H = -(A+B) has AT MOST three eigenvalues and its spectral projectors are
    P_{-2} = A B        P_{-1} = A(I-B) + (I-A)B        P_0 = (I-A)(I-B)
which are EXACT RATIONALS with denominator |G|.  No eigen-solver, no tolerance, no sampling.

WHAT IS DECIDED
  D1  A_h == I exactly?                      -- an EXACT ARGUMENT for every abelian G, at every |G|
  D2  dim E, exactly, per eigenspace         -- and its PARITY, which gates clause (iv)
  D3  dim of the H-commutant   = sum_E (dim E)^2
      dim of the H+transport commutant = sum_E (1/|G|) sum_h |chi_E(h)|^2
      the GAP is exactly zero  <=>  every record automatically commutes with all A_h
  D4  an EXACT RATIONAL WITNESS RECORD R (Hermitian, R^2 = I, [H,R] = 0, Tr(P_E R) = 0 on every E)
      with [A_h, R] computed EXACTLY, reporting the exact rational ||[A_h,R]||_F^2

CONTROL IN THE SAME TABLE (D-15)
  D(Z_2), D(Z_4), D(Z_2xZ_2): abelian; the SAME construction must return EXACTLY ZERO.
  D(D_4), D(Q_8): non-abelian; the SAME construction must return EXACTLY NON-ZERO.
  D(Z_3): |G| not a power of two; C-41 says no record exists at all -- the parity column must say so.
A method that does not split these classifies nothing.
"""
import sys
from fractions import Fraction as F

LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_EXACT_A_ZERO"
OUT = []
def say(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s); print(s); sys.stdout.flush()


# ------------------------------------------------------------------ groups (exact, tiny)
def group_Zn(n):
    return list(range(n)), (lambda a, b: (a + b) % n), "Z_%d" % n


def group_ZmxZn(m, n):
    return ([(a, b) for a in range(m) for b in range(n)],
            (lambda x, y: ((x[0] + y[0]) % m, (x[1] + y[1]) % n)), "Z_%dxZ_%d" % (m, n))


def group_dihedral(k):
    els = [(a, b) for a in range(k) for b in range(2)]
    def mul(x, y):
        a1, b1 = x; a2, b2 = y
        return ((a1 + (a2 if b1 == 0 else -a2)) % k, (b1 + b2) % 2)
    return els, mul, "D_%d" % k


def group_Q8():
    els = [(s, u) for s in (0, 1) for u in range(4)]
    base = {(0, 0): (0, 0), (0, 1): (0, 1), (0, 2): (0, 2), (0, 3): (0, 3),
            (1, 0): (0, 1), (2, 0): (0, 2), (3, 0): (0, 3),
            (1, 1): (1, 0), (2, 2): (1, 0), (3, 3): (1, 0),
            (1, 2): (0, 3), (2, 1): (1, 3),
            (2, 3): (0, 1), (3, 2): (1, 1),
            (3, 1): (0, 2), (1, 3): (1, 2)}
    def mul(x, y):
        s1, u1 = x; s2, u2 = y
        s3, u3 = base[(u1, u2)]
        return ((s1 + s2 + s3) % 2, u3)
    return els, mul, "Q_8"


# ------------------------------------------------------------------ exact matrix helpers
def zeros(N):
    return [[0] * N for _ in range(N)]


def ident(N):
    M = zeros(N)
    for i in range(N):
        M[i][i] = 1
    return M


def mmul(A, B, N=None):
    Bt = list(zip(*B))
    return [[sum(a * b for a, b in zip(row, col)) for col in Bt] for row in A]


def madd(A, B, N=None):
    return [[x + y for x, y in zip(ra, rb)] for ra, rb in zip(A, B)]


def msub(A, B, N=None):
    return [[x - y for x, y in zip(ra, rb)] for ra, rb in zip(A, B)]


def mscale(A, s, N=None):
    return [[x * s for x in ra] for ra in A]


def mtrace(A, N=None):
    return sum(A[i][i] for i in range(len(A)))


def mzero(A):
    return all(v == 0 for row in A for v in row)


def mfrob2(A):
    return sum(v * v for row in A for v in row)


# ------------------------------------------------------------------ the carrier
def build(els, mul, name):
    n = len(els)
    idx = {g: i for i, g in enumerate(els)}
    e = next(g for g in els if mul(g, g) == g)
    inv = {}
    for g in els:
        inv[g] = next(h for h in els if mul(g, h) == e)
    N = n * n

    perms = []
    for h in els:
        sigma = [0] * N
        for i1, g1 in enumerate(els):
            for i2, g2 in enumerate(els):
                a = mul(mul(h, g1), inv[h])
                b = mul(mul(h, g2), inv[h])
                sigma[i1 * n + i2] = idx[a] * n + idx[b]      # column i1*n+i2 -> row sigma[...]
        perms.append(sigma)

    # NA = |G| * A  as an EXACT INTEGER matrix
    NA = zeros(N)
    for sigma in perms:
        for c in range(N):
            NA[sigma[c]][c] += 1

    Bd = [0] * N
    for i1, g1 in enumerate(els):
        for i2, g2 in enumerate(els):
            comm = mul(mul(g1, g2), mul(inv[g1], inv[g2]))
            Bd[i1 * n + i2] = 1 if comm == e else 0
    Bm = zeros(N)
    for i in range(N):
        Bm[i][i] = Bd[i]
    return dict(name=name, n=n, N=N, els=els, mul=mul, idx=idx, inv=inv, e=e,
                perms=perms, NA=NA, B=Bm, Bd=Bd)


def perm_matrix(sigma, N):
    M = zeros(N)
    for c in range(N):
        M[sigma[c]][c] = 1
    return M


# ------------------------------------------------------------------ exact rational linear algebra
def rref_fraction(rows, ncols):
    rows = [list(r) for r in rows]
    piv = []
    r = 0
    for c in range(ncols):
        p = next((i for i in range(r, len(rows)) if rows[i][c] != 0), None)
        if p is None:
            continue
        rows[r], rows[p] = rows[p], rows[r]
        pv = rows[r][c]
        rows[r] = [x / pv for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][c] != 0:
                f = rows[i][c]
                rows[i] = [x - f * y for x, y in zip(rows[i], rows[r])]
        piv.append(c)
        r += 1
        if r == len(rows):
            break
    return rows[:r], piv


def range_basis(Pnum, den, N):
    """EXACT rational basis of the range of the projector P = Pnum/den, taken from its columns."""
    cols = [[F(Pnum[i][c], den) for i in range(N)] for c in range(N)]
    basis = []
    # exact Gaussian elimination on the column set
    pivots = []
    reduced = []
    for col in cols:
        v = list(col)
        for (pi, pv) in zip(pivots, reduced):
            if v[pi] != 0:
                f = v[pi]
                v = [a - f * b for a, b in zip(v, pv)]
        nz = next((i for i in range(N) if v[i] != 0), None)
        if nz is None:
            continue
        v = [a / v[nz] for a in v]
        pivots.append(nz); reduced.append(v); basis.append(list(col))
    return basis


def gram_schmidt_exact(vs, N):
    """Exact rational Gram-Schmidt WITHOUT normalisation -- stays inside Q, no square roots."""
    out = []
    for v in vs:
        w = list(v)
        for u in out:
            num = sum(a * b for a, b in zip(w, u))
            if num != 0:
                den = sum(a * a for a in u)
                f = num / den
                w = [a - f * b for a, b in zip(w, u)]
        if any(x != 0 for x in w):
            out.append(w)
    return out


def build_witness_record(P_list, N):
    """R = sum_E (P_E - 2 Q_E), Q_E an exact rational orthogonal projector onto HALF of E.
       Hermitian and R^2 = I by construction; Tr(P_E R) = 0 requires dim E even."""
    R = [[F(0) for _ in range(N)] for _ in range(N)]
    ok = True
    for (Pnum, den, dimE) in P_list:
        if dimE == 0:
            continue
        if dimE % 2 != 0:
            ok = False
            continue
        basis = range_basis(Pnum, den, N)
        og = gram_schmidt_exact(basis, N)
        assert len(og) == dimE, "range basis has the wrong dimension: %d vs %d" % (len(og), dimE)
        # P_E - 2 Q_E  with Q_E = sum over the first half
        for i in range(N):
            for j in range(N):
                R[i][j] += F(Pnum[i][j], den)
        for u in og[:dimE // 2]:
            d = sum(a * a for a in u)
            for i in range(N):
                if u[i] == 0:
                    continue
                ui = u[i]
                for j in range(N):
                    if u[j] != 0:
                        R[i][j] -= 2 * ui * u[j] / d
    return R, ok


# ------------------------------------------------------------------ main
say("=" * 128)
say("E-4  D(G) GAUGE TRANSPORT AGAINST A RECORD -- EXACT RATIONAL ARITHMETIC (fractions.Fraction)")
say("=" * 128)

CARRIERS = [group_Zn(2), group_Zn(3), group_Zn(4), group_ZmxZn(2, 2),
            group_dihedral(4), group_Q8(), group_Zn(8), group_Zn(5)]

say("")
say("STEP 1   D1/D2  -- is A_h the identity?  exact eigenspace dimensions and their parity")
say("-" * 128)
say("  %-12s %-6s %-6s %-9s %-34s %-30s %-18s"
    % ("carrier", "|G|", "dim", "abelian", "exact eigenvalue:dim(E)", "all dim E even? (clause iv)",
       "A_h == I exactly?"))
rows = []
for els, mul, name in CARRIERS:
    C = build(els, mul, name)
    N = C["N"]; n = C["n"]
    NA = C["NA"]; Bm = C["B"]
    I = ident(N)
    nI = mscale(I, n, N)
    # exact spectral projectors, numerators over denominator n
    IB = msub(I, Bm, N)
    nIA = msub(nI, NA, N)
    P2num = mmul(NA, Bm, N)                       # A B          (eig -2)
    P1num = madd(mmul(NA, IB, N), mmul(nIA, Bm, N), N)   # A(I-B)+(I-A)B  (eig -1)
    P0num = mmul(nIA, IB, N)                      # (I-A)(I-B)   (eig 0)
    # exact self-checks
    tot = madd(madd(P2num, P1num, N), P0num, N)
    resolution_ok = mzero(msub(tot, nI, N))
    ab_ok = mzero(msub(mmul(NA, Bm, N), mmul(Bm, NA, N)))
    idem_ok = mzero(msub(mmul(NA, NA, N), mscale(NA, n, N)))
    dims = []
    for lbl, Pn in (("-2", P2num), ("-1", P1num), ("0", P0num)):
        t = F(mtrace(Pn, N), n)
        assert t.denominator == 1, "eigenspace dimension is not an integer"
        dims.append((lbl, int(t)))
    abelian = all(mul(g, h) == mul(h, g) for g in els for h in els)
    Ah_is_I = all(all(s[c] == c for c in range(N)) for s in C["perms"])
    all_even = all(d % 2 == 0 for _, d in dims)
    say("  %-12s %-6d %-6d %-9s %-34s %-30s %-18s"
        % ("D(%s)" % name, n, N, "yes" if abelian else "NO",
           " ".join("%s:%d" % t for t in dims),
           "YES" if all_even else "NO -- clause (iv) impossible on an odd E",
           "YES (exact)" if Ah_is_I else "no"))
    assert resolution_ok and ab_ok and idem_ok, "exact projector self-check failed for %s" % name
    rows.append(dict(C=C, name=name, n=n, N=N, abelian=abelian, Ah_is_I=Ah_is_I,
                     dims=dims, all_even=all_even,
                     P=[(P2num, n, dims[0][1]), (P1num, n, dims[1][1]), (P0num, n, dims[2][1])]))
say("  exact self-checks on every carrier: sum of projectors == I, [A,B] == 0, A^2 == A : ALL PASS")

say("")
say("=" * 128)
say("STEP 2   D3  -- EXACT COMMUTANT DIMENSIONS.  chi_E(h) = Tr(P_E A_h) is an exact rational;")
say("         dim End_G(E) = (1/|G|) sum_h chi_E(h) chi_E(h^-1) must be a positive INTEGER.")
say("         GAP = dim(commutant of H) - dim(commutant of H and all A_h).  GAP = 0 means every")
say("         record commutes with transport automatically; GAP > 0 means a generic record is MOVED.")
say("-" * 128)
say("  %-12s %-9s %-26s %-26s %-14s %-26s"
    % ("carrier", "abelian", "dim comm(H) = sum dimE^2", "dim comm(H,A_h) = sum m_i^2",
       "GAP (exact)", "verdict"))
for r in rows:
    C = r["C"]; N = r["N"]; n = r["n"]
    els = C["els"]; inv = C["inv"]; idx = C["idx"]
    total_H = sum(d * d for _, d in r["dims"])
    total_G = 0
    for (Pnum, den, dimE) in r["P"]:
        if dimE == 0:
            continue
        chi = {}
        for hi, h in enumerate(els):
            sigma = C["perms"][hi]
            # Tr(P A_h) = sum_i P[sigma(i)][i]
            chi[h] = F(sum(Pnum[sigma[i]][i] for i in range(N)), den)
        s = F(0)
        for h in els:
            s += chi[h] * chi[inv[h]]
        s = s / n
        assert s.denominator == 1 and s >= 1, "dim End_G(E) is not a positive integer: %s" % s
        total_G += int(s)
    gap = total_H - total_G
    say("  %-12s %-9s %-26d %-26d %-14d %-26s"
        % ("D(%s)" % r["name"], "yes" if r["abelian"] else "NO", total_H, total_G, gap,
           "transport acts trivially" if gap == 0 else "a generic record IS MOVED"))
    r["gap"] = gap

say("")
say("=" * 128)
say("STEP 3   D4  -- AN EXACT RATIONAL WITNESS RECORD, AND [A_h, R] IN EXACT RATIONALS.")
say("         R = sum_E (P_E - 2 Q_E), Q_E an exact rational orthogonal projector onto half of E.")
say("         Hermitian and R^2 = I by construction; every check below is an EXACT equality test.")
say("-" * 128)
say("  %-12s %-7s %-9s %-9s %-11s %-13s %-30s %-16s"
    % ("carrier", "built?", "R=R^T", "R^2=I", "[H,R]=0", "Tr(P_E R)=0", "max ||[A_h,R]||_F^2 (EXACT)",
       "classification"))
for r in rows:
    if r["N"] > 64:
        say("  %-12s %-7s %-9s %-9s %-11s %-13s %-30s %-16s"
            % ("D(%s)" % r["name"], "skipped", "-", "-", "-", "-",
               "dim %d too large for exact Fraction linear algebra" % r["N"], "-"))
        continue
    C = r["C"]; N = r["N"]; n = r["n"]
    if not r["all_even"]:
        say("  %-12s %-7s %-9s %-9s %-11s %-13s %-30s %-16s"
            % ("D(%s)" % r["name"], "NO", "-", "-", "-", "impossible",
               "no record exists: an eigenspace has ODD dimension", "(Z) vacuous"))
        continue
    R, ok = build_witness_record(r["P"], N)
    herm = all(R[i][j] == R[j][i] for i in range(N) for j in range(N))
    R2 = mmul(R, R, N)
    Iq = [[F(1) if i == j else F(0) for j in range(N)] for i in range(N)]
    sq = mzero(msub(R2, Iq, N))
    # [H,R] = 0  <=>  [A,R] = 0 and [B,R] = 0
    NAq = [[F(v) for v in row] for row in C["NA"]]
    Bq = [[F(v) for v in row] for row in C["B"]]
    commH = mzero(msub(mmul(NAq, R, N), mmul(R, NAq, N))) and \
            mzero(msub(mmul(Bq, R, N), mmul(R, Bq, N)))
    trs = []
    for (Pnum, den, dimE) in r["P"]:
        if dimE == 0:
            trs.append(F(0)); continue
        Pq = [[F(Pnum[i][j], den) for j in range(N)] for i in range(N)]
        trs.append(mtrace(mmul(Pq, R, N), N))
    tr_ok = all(t == 0 for t in trs)
    worst = F(0)
    for sigma in C["perms"]:
        Ah = [[F(v) for v in row] for row in perm_matrix(sigma, N)]
        Cm = msub(mmul(Ah, R, N), mmul(R, Ah, N))
        f2 = sum(v * v for row in Cm for v in row)
        if f2 > worst:
            worst = f2
    cls = "(NZ) exactly non-zero" if worst != 0 else "(Z) EXACTLY ZERO"
    say("  %-12s %-7s %-9s %-9s %-11s %-13s %-30s %-16s"
        % ("D(%s)" % r["name"], "YES", herm, sq, commH, tr_ok, str(worst), cls))
    r["worst"] = worst

say("")
say("=" * 128)
say("  E-4 SUMMARY -- classification with the control carried in the same table")
say("=" * 128)
say("  %-12s %-9s %-16s %-14s %-34s" % ("carrier", "abelian", "A_h == I?", "GAP", "||[A_h,R]||_F^2 verdict"))
for r in rows:
    say("  %-12s %-9s %-16s %-14s %-34s"
        % ("D(%s)" % r["name"], "yes" if r["abelian"] else "NO",
           "YES" if r["Ah_is_I"] else "no", r.get("gap", "-"),
           ("EXACTLY 0" if r.get("worst") == 0 else str(r.get("worst")))
           if "worst" in r else ("no record: odd eigenspace" if not r["all_even"] else "not built")))
say("=" * 128)

with open(LANE + "/e4_dg_exact.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")

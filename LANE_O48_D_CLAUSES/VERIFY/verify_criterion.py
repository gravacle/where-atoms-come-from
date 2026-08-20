"""ATTACK ON THE HEADLINE'S STRONGEST POSITIVE CLAIM:
   "The true requirement is A PAIR OF ANTICOMMUTING SYMMETRIES OF H, one of which is the record --
    this predicted the record count exactly on all 19 carriers tested ... 19 carriers, 19
    agreements, 0 exceptions."

PREDICTED (lane) = #{Pauli R != I : [R,H]=0 and SOME PAULI W has [W,H]=0 and WR=-RW}
MEASURED (lane)  = #{Pauli R : clauses (i)-(iv), (iv) via Tr(P_E R)=0 on every eigenspace}

The lane tested 19 HAND-PICKED carriers.  Here the same two counts are run on a RANDOM SAMPLE
of Pauli Hamiltonians -- commuting and non-commuting -- to see whether the agreement is a
theorem or a property of the sample.  Any row with PREDICTED != MEASURED is an exception.
"""
import numpy as np
from itertools import product

OUT = []
def say(*a):
    s = " ".join(str(x) for x in a); OUT.append(s); print(s)

I2 = np.eye(2, dtype=complex)
Xm = np.array([[0, 1], [1, 0]], dtype=complex)
Zm = np.array([[1, 0], [0, -1]], dtype=complex)
Ym = 1j * Xm @ Zm

def pmat(a, b):
    M = np.array([[1]], dtype=complex)
    for i in range(len(a)):
        x, z = a[i], b[i]
        P = I2 if (x, z) == (0, 0) else (Xm if (x, z) == (1, 0) else (Zm if (x, z) == (0, 1) else Ym))
        M = np.kron(M, P)
    return M

def symp(a1, b1, a2, b2):
    return sum(a1[i] * b2[i] + b1[i] * a2[i] for i in range(len(a1))) % 2

def analyse(n, terms, tol=1e-8):
    """terms = [(a,b,J)].  Returns (predicted, measured, n_commuting, exceptions list)."""
    H = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for a, b, J in terms:
        H = H + J * pmat(a, b)
    w, V = np.linalg.eigh(H)
    # group eigenvalues
    groups, cur = [], [0]
    for k in range(1, len(w)):
        if abs(w[k] - w[k - 1]) < tol:
            cur.append(k)
        else:
            groups.append(cur); cur = [k]
    groups.append(cur)
    words = [(a, b) for a in product((0, 1), repeat=n) for b in product((0, 1), repeat=n)]
    adm = [(a, b) for a, b in words
           if not any(symp(a, b, ta, tb) for ta, tb, _ in terms)]
    admset = [(a, b) for a, b in adm if not (all(x == 0 for x in a) and all(x == 0 for x in b))]
    pred = 0; meas = 0; exc = []
    for a, b in admset:
        has_w = any(symp(a, b, wa, wb) == 1 for wa, wb in adm)
        R = pmat(a, b)
        Rp = V.conj().T @ R @ V
        ok4 = True; nonconst = False
        for g in groups:
            blk = Rp[np.ix_(g, g)]
            tr = np.trace(blk).real
            if abs(tr) > 1e-6:
                ok4 = False
            ev = np.linalg.eigvalsh((blk + blk.conj().T) / 2)
            if len(g) > 1 and (ev.max() - ev.min()) > 1e-6:
                nonconst = True
            if len(g) == 1:
                pass
        m = ok4 and nonconst
        if has_w: pred += 1
        if m: meas += 1
        if has_w != m:
            exc.append(((a, b), has_w, m))
    return pred, meas, len(admset), exc


say("=" * 100)
say("RANDOM-SAMPLE TEST OF THE 'ANTICOMMUTING SYMMETRY PAIR' CRITERION")
say("=" * 100)
say("")
say(f"{'family':>26} {'n':>3} {'#terms':>7} {'#adm R':>7} {'PREDICTED':>10} {'MEASURED':>9} {'agree?':>7}")

rng = np.random.default_rng(11)
tot = 0; bad = 0; badrows = []
for trial in range(120):
    n = int(rng.integers(2, 5))
    k = int(rng.integers(1, 5))
    seen = set(); terms = []
    while len(terms) < k:
        a = tuple(int(x) for x in rng.integers(0, 2, n))
        b = tuple(int(x) for x in rng.integers(0, 2, n))
        if (a, b) in seen or (all(x == 0 for x in a) and all(x == 0 for x in b)):
            continue
        seen.add((a, b))
        terms.append((a, b, int(rng.integers(1, 6))))
    comm = all(symp(t1[0], t1[1], t2[0], t2[1]) == 0
               for i, t1 in enumerate(terms) for t2 in terms[i + 1:])
    p, m, na, exc = analyse(n, terms)
    tot += 1
    if p != m:
        bad += 1
        badrows.append((n, terms, p, m, exc[:4]))
    if trial < 14 or p != m:
        say(f"{('random ' + ('COMM' if comm else 'NONCOMM')):>26} {n:>3} {k:>7} {na:>7} "
            f"{p:>10} {m:>9} {str(p == m):>7}")

say("")
say(f"  RANDOM SAMPLE: {tot} Pauli Hamiltonians, n = 2..4, 1..4 terms, commuting and not.")
say(f"  EXCEPTIONS (PREDICTED != MEASURED): {bad}")

say("")
say("  NAMED CARRIERS INCLUDING ONES THE LANE DID NOT TEST:")
named = []
# AF Ising triangle -- the model that broke a previous headline in this program
named.append(("AF Ising triangle Z_iZ_j", 3,
              [((0,0,0),(1,1,0),1), ((0,0,0),(0,1,1),1), ((0,0,0),(1,0,1),1)]))
named.append(("frustrated square ZZ", 4,
              [((0,0,0,0),(1,1,0,0),1), ((0,0,0,0),(0,1,1,0),1),
               ((0,0,0,0),(0,0,1,1),1), ((0,0,0,0),(1,0,0,1),1)]))
named.append(("Heisenberg XXX 2-site", 2,
              [((1,1),(0,0),1), ((1,1),(1,1),1), ((0,0),(1,1),1)]))
named.append(("Heisenberg XXX 3-site chain", 3,
              [((1,1,0),(0,0,0),1), ((1,1,0),(1,1,0),1), ((0,0,0),(1,1,0),1),
               ((0,1,1),(0,0,0),1), ((0,1,1),(0,1,1),1), ((0,0,0),(0,1,1),1)]))
named.append(("transverse-field Ising n=3", 3,
              [((0,0,0),(1,1,0),1), ((0,0,0),(0,1,1),1),
               ((1,0,0),(0,0,0),1), ((0,1,0),(0,0,0),1), ((0,0,1),(0,0,0),1)]))
named.append(("4-qubit ZZ ring (GEN J)", 4,
              [((0,0,0,0),(1,1,0,0),1), ((0,0,0,0),(0,1,1,0),2),
               ((0,0,0,0),(0,0,1,1),4), ((0,0,0,0),(1,0,0,1),8)]))
for nm, n, terms in named:
    p, m, na, exc = analyse(n, terms)
    say(f"{nm:>32} {n:>3} {len(terms):>7} {na:>7} {p:>10} {m:>9} {str(p == m):>7}")
    if p != m:
        for e in exc[:5]:
            say(f"      EXCEPTION word a={e[0][0]} b={e[0][1]}  pauli-witness={e[1]}  clauses={e[2]}")

with open(__file__.replace(".py", ".txt"), "w") as f:
    f.write("\n".join(OUT) + "\n")

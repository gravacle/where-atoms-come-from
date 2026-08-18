"""LANE O-2/O-3, PART 3.  WHICH CLAUSE DIES WHEN A DEGENERACY SPLITS AT FIRST ORDER?

ROW O-3 registered this SIDE-STEP, at LOW cost:
      "prove any first-order-split degeneracy cannot carry a noise-surviving record."
This script tests whether that is right.  It is NOT.  The side-step is FALSE AS WRITTEN, and this
script exhibits two explicit counterexamples in which a degeneracy splits at FIRST order and a record
survives ALL FIVE CLAUSES.  It then states the repair, which is still low cost but is a DIFFERENT
argument from the one registered.

THE THREE THINGS THAT CAN HAPPEN, stated exactly.  H(eps) = H0 + eps V; R a record for H0, P the
projector on the degenerate eigenspace E0; Phi(V) = PVP - (tr PVP/n)P.

  (A)  IF [V,R] != 0:  CLAUSE (ii) DIES, and it dies LINEARLY --
              ||[H(eps),R]|| = eps ||[V,R]||
       for every eps, WHATEVER the splitting exponent.  Note the consequence: on a distance-d code a
       weight-1 perturbation splits the level only at order eps^d, but it violates (ii) at order eps.
       (ii) and the splitting exponent are DIFFERENT clocks.  This is worth stating because the
       program has been using the splitting exponent as if it measured the record's survival.

  (B)  IF Phi(V) HAS NON-DEGENERATE SPECTRUM (the multiplet is FULLY resolved at first order):
       CLAUSE (iii) DIES, and no deformation of R rescues it.  Every eigenspace descended from E0 is
       one-dimensional, so any R obeying (i)+(ii) is constant on each of them, i.e. R is a function
       of H(eps).  The record space is EMPTY, not merely "the old R broke".

  (C)  IF Phi(V) != 0 BUT COMMUTES WITH R AND IS DEGENERATE ON R's BLOCKS:
       NOTHING DIES.  The degeneracy splits at first order and R still satisfies (i)-(v).
       This is the counterexample to the registered side-step.

  AND CLAUSE (v) DOES NOT DEPEND ON eps OR ON V AT ALL.  It is a property of the pair (E0, R) and the
  set of contractible operations.  So "the splitting kills the protection" is a category error: the
  splitting cannot kill (v), because (v) never mentioned the perturbation.  What is TRUE is an
  EQUIVALENCE at fixed E0, proved in Part 2 --
        some contractible V has Phi(V) != 0   <=>   the Knill-Laflamme condition fails at
        contractible weight   <=>   a contractible operator acts on E0 as a non-trivial logical.
  That equivalence, not the implication registered in O-3, is what the exhaustiveness argument needs.

A NOTE ON CLAUSE (v) AS IT WAS PREVIOUSLY MEASURED (LANE_P1_DEFINITION).  There, (v) was scored as
      min over local L of || L^dag R L + R ||
with the local set L = {Z on one link}.  That set consists entirely of operators that COMMUTE with a
Z-type record, so the test could not have failed.  Swept over ALL single-qubit operators the same
bare operator identity FAILS even for the toric code: X on any link of the record's cycle sends
Z_cycle -> -Z_cycle exactly.  What survives is the CODE-SPACE form of (v): no contractible operator
ACTS ON E0 as a logical operation anticommuting with R.  This script reports both, and flags that (v)
is not well posed until obstruction O-4 ("admissible" is undefined) is settled.
"""
import itertools
import numpy as np

rng = np.random.default_rng(90210)
I2 = np.eye(2, dtype=complex)
SX = np.array([[0, 1], [1, 0]], dtype=complex)
SY = np.array([[0, -1j], [1j, 0]], dtype=complex)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
PAULI = {'I': I2, 'X': SX, 'Y': SY, 'Z': SZ}


def pmat(s):
    M = np.array([[1.0 + 0j]])
    for ch in s:
        M = np.kron(M, PAULI[ch])
    return M


def weight(s):
    return sum(1 for ch in s if ch != 'I')


def symp(a, b):
    c = 0
    for x, y in zip(a, b):
        if x != 'I' and y != 'I' and x != y:
            c ^= 1
    return c


def pmul(a, b):
    out = []
    for x, y in zip(a, b):
        if x == 'I':
            out.append(y)
        elif y == 'I':
            out.append(x)
        elif x == y:
            out.append('I')
        else:
            out.append({'XY': 'Z', 'YX': 'Z', 'XZ': 'Y', 'ZX': 'Y', 'YZ': 'X', 'ZY': 'X'}[x + y])
    return ''.join(out)


def stab_group(gens, n):
    S = {'I' * n}
    frontier = ['I' * n]
    while frontier:
        nxt = []
        for a in frontier:
            for g in gens:
                b = pmul(a, g)
                if b not in S:
                    S.add(b); nxt.append(b)
        frontier = nxt
    return S


def code_H(gens, n):
    N = 2 ** n
    H = np.zeros((N, N), complex)
    for g in gens:
        H -= pmat(g)
    H = (H + H.conj().T) / 2
    ev, U = np.linalg.eigh(H)
    tol = 1e-8 * max(1.0, abs(ev).max())
    sel = np.abs(ev - ev[0]) < tol
    Q = U[:, sel]
    return H, Q, Q @ Q.conj().T


def phi(V, Q):
    B = Q.conj().T @ V @ Q
    m = B.shape[0]
    return B - np.trace(B) / m * np.eye(m)


def eigenspaces(H, tol=1e-9):
    ev, U = np.linalg.eigh(H)
    out = []
    i = 0
    while i < len(ev):
        j = i
        while j + 1 < len(ev) and abs(ev[j + 1] - ev[i]) < tol * max(1.0, abs(ev[i])):
            j += 1
        out.append((ev[i], U[:, i:j + 1]))
        i = j + 1
    return out


def clause_report(R, H, Wop, contractible, name, code_space_v=None, verbose=True):
    """the five clauses, using the LANE_P1_DEFINITION measures."""
    N = R.shape[0]
    c1 = max(np.linalg.norm(R - R.conj().T), np.linalg.norm(R @ R - np.eye(N)))
    c2 = np.linalg.norm(H @ R - R @ H)
    nontriv = 0.0
    for _, Uk in eigenspaces(H):
        b = Uk.conj().T @ R @ Uk
        if b.shape[0] > 1:
            nontriv = max(nontriv, float(np.linalg.norm(b - np.trace(b) / b.shape[0] * np.eye(b.shape[0]))))
    c4 = np.linalg.norm(Wop.conj().T @ R @ Wop + R)
    bare_v = min(float(np.linalg.norm(L.conj().T @ R @ L + R)) for L in contractible)
    if verbose:
        print(f"    {name}")
        print(f"      (i)   bit         max(||R-R+||,||R^2-I||)  = {c1:.2e}   {'PASS' if c1 < 1e-9 else 'FAIL'}")
        print(f"      (ii)  durable     ||[H,R]||                = {c2:.3e}   {'PASS' if c2 < 1e-9 else 'FAIL'}")
        print(f"      (iii) nontrivial  max deviation from scalar= {nontriv:.4f}    {'PASS' if nontriv > 1e-6 else 'FAIL'}")
        print(f"      (iv)  writable    ||W+RW + R||             = {c4:.2e}   {'PASS' if c4 < 1e-9 else 'FAIL'}")
        print(f"      (v)   protected   BARE operator form, min over contractible L of ||L+RL+R|| "
              f"= {bare_v:.3f}  {'PASS' if bare_v > 1e-6 else 'FAIL'}")
        if code_space_v is not None:
            print(f"      (v')  protected   CODE-SPACE form: is there a contractible operator acting on E0")
            print(f"                        as a logical anticommuting with R?  {code_space_v}")
    return dict(c1=c1, c2=c2, iii=nontriv, c4=c4, v=bare_v)


def code_space_protection(gens, n, S, Rstring, wmax):
    """search for a Pauli of weight <= wmax that acts on the code space as a logical ANTICOMMUTING
    with the record R (given as a Pauli string).  Returns the witness or None."""
    for w in range(1, wmax + 1):
        for supp in itertools.combinations(range(n), w):
            for letters in itertools.product('XYZ', repeat=w):
                s = ['I'] * n
                for k, ch in zip(supp, letters):
                    s[k] = ch
                s = ''.join(s)
                if all(symp(s, g) == 0 for g in gens) and s not in S and symp(s, Rstring) == 1:
                    return s
    return None


print("=" * 104)
print("LANE O-3.  WHICH CLAUSE DIES WHEN A FIRST-ORDER-SPLIT DEGENERACY IS PERTURBED?")
print("=" * 104)

# ------------------------------------------------------------------ the carrier: toric code 2x2
h = lambda i, j: (j % 2) * 2 + (i % 2)
v = lambda i, j: 4 + (j % 2) * 2 + (i % 2)
gens = []
for j in range(2):
    for i in range(2):
        s = ['I'] * 8
        for L in [h(i, j), h(i - 1, j), v(i, j), v(i, j - 1)]:
            s[L] = 'X' if s[L] == 'I' else 'I'
        gens.append(''.join(s))
        s = ['I'] * 8
        for L in [h(i, j), v(i + 1, j), h(i, j + 1), v(i, j)]:
            s[L] = 'Z' if s[L] == 'I' else 'I'
        gens.append(''.join(s))


def independent(gs, n):
    rows, keep = [], []
    for g in gs:
        vec = 0
        for k, ch in enumerate(g):
            if ch in 'XY':
                vec |= 1 << (2 * k)
            if ch in 'ZY':
                vec |= 1 << (2 * k + 1)
        cur = vec
        for r in rows:
            p = r.bit_length() - 1
            if cur >> p & 1:
                cur ^= r
        if cur:
            rows.append(cur); rows.sort(reverse=True); keep.append(g)
    return keep


gens = independent(gens, 8)
n = 8
S = stab_group(gens, n)
H0, Q, P = code_H(gens, n)
m = Q.shape[1]
print(f"\nCARRIER: toric code on the 2x2 torus.  n={n} qubits, {len(gens)} independent stabilisers,")
print(f"         code space dimension {m} = 2^{n - len(gens)}, distance d = 2 (verified in Part 2).")

# find two commuting logical Z-type operators of weight 2, and their conjugate X-type partners
logicals = []
for w in (2,):
    for supp in itertools.combinations(range(n), w):
        for letters in itertools.product('XYZ', repeat=w):
            s = ['I'] * n
            for k, ch in zip(supp, letters):
                s[k] = ch
            s = ''.join(s)
            if all(symp(s, g) == 0 for g in gens) and s not in S:
                logicals.append(s)
# search for a genuine symplectic basis Z1,X1,Z2,X2 among the weight-2 logicals:
#   {Z1,X1} anticommute, {Z2,X2} anticommute, and the two pairs commute with each other.
Z1 = X1 = Z2 = X2 = None
for a in logicals:
    for b in logicals:
        if symp(a, b) != 1:
            continue
        for c in logicals:
            if symp(c, a) or symp(c, b):
                continue
            for dd in logicals:
                if symp(dd, a) or symp(dd, b) or symp(dd, c) != 1:
                    continue
                Z1, X1, Z2, X2 = a, b, c, dd
                break
            if Z2 is not None and X2 is not None:
                break
        if X2 is not None:
            break
    if X2 is not None:
        break
assert X2 is not None, "no symplectic basis of weight-2 logicals found"
print(f"         logical operators found: Zbar1 = {Z1}   Xbar1 = {X1}   Zbar2 = {Z2}   Xbar2 = {X2}")
print(f"         commutation check: [Zbar1,Zbar2] = {'0' if symp(Z1, Z2) == 0 else 'NOT 0'},  "
      f"{{Zbar2,Xbar2}} = {'anticommute' if symp(Z2, X2) == 1 else 'COMMUTE (BAD)'}")

contract1 = [pmat(s) for s in
             [''.join('XYZ'[c] if k == q else 'I' for k in range(n)) for q in range(n) for c in range(3)]]

print("\nSELF-CHECK -- all five clauses at eps = 0, on the unperturbed carrier.")
Rmat = pmat(Z2)
Wmat = pmat(X2)
wit = code_space_protection(gens, n, S, Z2, 1)
_ = clause_report(Rmat, H0, Wmat, contract1, "R = Zbar2, W = Xbar2, eps = 0",
                  code_space_v=("NO -- none exists at weight 1" if wit is None else f"YES: {wit}"))
print("      NOTE the (v) row: the BARE operator form already FAILS, because X on a link of the")
print("      record's cycle flips it exactly.  The P1 lane scored (v) as PASS only because its local")
print("      set was Z-type, which commutes with a Z-type record by construction.  This is an")
print("      ERRATUM-GRADE observation about how (v) has been measured, and it is obstruction O-4.")

# ============================================================ CASE A
print("\n" + "-" * 104)
print("CASE A.  CONTRACTIBLE (weight-1) PERTURBATION.  Splitting is second order (d=2).  What dies?")
print("-" * 104)
VA = sum(rng.normal() * pmat(s) for s in
         [''.join('XYZ'[c] if k == q else 'I' for k in range(n)) for q in range(n) for c in range(3)])
VA = (VA + VA.conj().T) / 2
VA = VA / np.linalg.norm(VA) * np.linalg.norm(H0)
print(f"    ||Phi(VA)|| on E0 = {np.linalg.norm(phi(VA, Q)):.3e}   (zero: no first-order splitting, d=2)")
comm = np.linalg.norm(VA @ Rmat - Rmat @ VA)
print(f"    ||[VA,R]|| = {comm:.4f}")
def e0_cluster(He, P0, mm):
    """the mm eigenvalues descended from E0, by overlap. Tolerance-free reporting."""
    ev, U = np.linalg.eigh(He)
    ov = np.einsum('ij,jk,ki->i', U.conj().T, P0, U).real
    idx = np.argsort(-ov)[:mm]
    return np.sort(ev[idx])


def min_gap_whole(He):
    ev = np.linalg.eigvalsh(He)
    d = np.diff(np.sort(ev))
    return float(d.min())


print(f"\n    {'eps':>8s} {'splitting':>14s} {'||[H(eps),R]||':>16s} {'eps*||[VA,R]||':>16s} "
      f"{'ratio':>9s} {'min gap in E0 cluster':>23s} {'min gap, WHOLE spectrum':>25s}")
for eps in (1e-4, 1e-3, 1e-2, 1e-1):
    He = H0 + eps * VA
    e = e0_cluster(He, P, m)
    spl = e[-1] - e[0]
    c2 = np.linalg.norm(He @ Rmat - Rmat @ He)
    print(f"    {eps:8.1e} {spl:14.4e} {c2:16.6e} {eps * comm:16.6e} "
          f"{c2 / (eps * comm):9.6f} {float(np.diff(e).min()):23.4e} {min_gap_whole(He):25.4e}")
print("    => CLAUSE (ii) FAILS LINEARLY IN eps (ratio 1.000000), while the LEVEL splitting is")
print("       quadratic.  The record's durability and the level's splitting run on different clocks.")
He = H0 + 1e-3 * VA
print()
_ = clause_report(Rmat, He, Wmat, contract1, "R = Zbar2 (unchanged), H = H0 + 1e-3 VA")
print(f"      the WHOLE spectrum of H(eps) is non-degenerate (min gap over all 256 levels = "
      f"{min_gap_whole(He):.3e}),")
print(f"      so EVERY eigenspace is one-dimensional and NO R obeying (i)+(ii) can satisfy (iii).")
print(f"      => under a GENERIC contractible perturbation both failures occur at once: (ii) for the")
print(f"         original R, (iii) for every deformed one.  Nothing rescues the record.")

# ============================================================ CASE B
print("\n" + "-" * 104)
print("CASE B.  A PERTURBATION THAT FULLY RESOLVES THE MULTIPLET AT FIRST ORDER.  What dies?")
print("-" * 104)
VB = 1.0 * pmat(Z1) + 0.37 * pmat(X1) + 0.61 * pmat(Z2) + 0.23 * pmat(X2)
VB = (VB + VB.conj().T) / 2
ph = phi(VB, Q)
evph = np.linalg.eigvalsh(ph)
print(f"    Phi(VB) eigenvalues on E0 = {np.round(evph, 6)}")
print(f"    non-degenerate? {'YES' if np.min(np.diff(np.sort(evph))) > 1e-8 else 'NO'}  "
      f"(min gap {np.min(np.diff(np.sort(evph))):.3e})")
for eps in (1e-4, 1e-3, 1e-2):
    He = H0 + eps * VB
    e = e0_cluster(He, P, m)
    print(f"    eps = {eps:.1e}:  E0-cluster eigenvalues = {np.round(e, 8)}   "
          f"min gap within the cluster = {float(np.diff(e).min()):.4e}")
print("    => every eigenspace descended from E0 is ONE-DIMENSIONAL.  Any R with [H(eps),R]=0 is")
print("       diagonal in that basis, hence CONSTANT on each of them: on the former E0, R is a")
print("       function of H(eps).  CLAUSE (iii) DIES THERE FOR EVERY R -- no deformation rescues it.")
print("       THIS is the only case in which the registered O-3 side-step is correct.")
print("    HONEST CAVEAT: clause (iii) quantifies over SOME eigenspace of H, not over E0, and this")
print(f"       V_B is built from logical operators only, so the EXCITED syndrome sectors are")
print(f"       untouched and remain degenerate (largest remaining eigenspace has dimension "
      f"{max(Uk.shape[1] for _, Uk in eigenspaces(H0 + 1e-3 * VB))}).")
print("       A record could in principle live up there.  CASE A's generic contractible perturbation")
print("       is the one that removes every degeneracy in the whole spectrum at once.")

# ============================================================ CASE C
print("\n" + "-" * 104)
print("CASE C.  COUNTEREXAMPLE 1 TO THE REGISTERED SIDE-STEP.  A perturbation that splits the")
print("         degeneracy AT FIRST ORDER while a record survives ALL FIVE CLAUSES.")
print("-" * 104)
VC = pmat(Z1)                      # a logical operator: weight 2 = d, splits at first order
print(f"    V = Zbar1 = {Z1}  (weight {weight(Z1)} = d).  Record R = Zbar2 = {Z2}, writer W = Xbar2.")
print(f"    ||Phi(V)|| on E0 = {np.linalg.norm(phi(VC, Q)):.4f}   -> FIRST-ORDER splitting")
for eps in (1e-4, 1e-3, 1e-2):
    He = H0 + eps * VC
    ev, U = np.linalg.eigh(He)
    ov = np.einsum('ij,jk,ki->i', U.conj().T, P, U).real
    idx = np.argsort(-ov)[:m]
    e = np.sort(ev[idx])
    print(f"      eps = {eps:.1e}:  splitting = {e[-1] - e[0]:.6e}   ratio to 2*eps = "
          f"{(e[-1] - e[0]) / (2 * eps):.6f}  <- LINEAR, first order")
He = H0 + 1e-3 * VC
print()
_ = clause_report(Rmat, He, Wmat, contract1, "R = Zbar2 under H = H0 + 1e-3 Zbar1",
                  code_space_v=("NO -- none exists at weight 1" if wit is None else f"YES: {wit}"))
print("    => (i),(ii),(iii),(iv) ALL PASS at a FIRST-ORDER-SPLIT degeneracy.  The splitting operator")
print("       COMMUTES with the record, so it resolves a DIFFERENT logical direction and leaves R's")
print("       own two-fold degeneracy intact.  THE REGISTERED SIDE-STEP IS FALSE AS WRITTEN.")

# ============================================================ CASE D
print("\n" + "-" * 104)
print("CASE D.  COUNTEREXAMPLE 2, AND THE STRONGER ONE: the splitting perturbation is CONTRACTIBLE")
print("         (weight 1) and still a record survives.  Carrier = [[3,1,1]] (x) [[5,1,3]].")
print("-" * 104)
gensD = ["ZZIIIIII", "IZZIIIII",
         "IIIXZZXI", "IIIIXZZX", "IIIXIXZZ", "IIIZXIXZ"]
nD = 8
SD = stab_group(gensD, nD)
HD, QD, PD = code_H(gensD, nD)
mD = QD.shape[1]
ZA = "ZIIIIIII"                    # logical Z of the fragile [[3,1,1]] block, WEIGHT 1
ZB = "IIIZZZZZ"                    # logical Z of the protected [[5,1,3]] block
XB = "IIIXXXXX"                    # logical X of the protected block
print(f"    n = {nD}, code space dim {mD} = 2^{nD - len(gensD)}.  Block A = [[3,1,1]] on qubits 1-3 (d=1),")
print(f"    Block B = [[5,1,3]] on qubits 4-8 (d=3).  OVERALL code distance = 1.")
for s, lab in [(ZA, "Zbar_A"), (ZB, "Zbar_B"), (XB, "Xbar_B")]:
    ok = all(symp(s, g) == 0 for g in gensD) and s not in SD
    print(f"      {lab} = {s}  weight {weight(s)}   is a non-trivial logical: {ok}")
print(f"      [Zbar_A, Zbar_B] = {'0' if symp(ZA, ZB) == 0 else 'NOT 0'}   "
      f"{{Zbar_B, Xbar_B}} = {'anticommute' if symp(ZB, XB) == 1 else 'COMMUTE (BAD)'}")
VD = pmat(ZA)
print(f"\n    ||Phi(V)|| for V = Zbar_A (WEIGHT 1, CONTRACTIBLE) = {np.linalg.norm(phi(VD, QD)):.4f}"
      f"   -> FIRST-ORDER splitting by a contractible perturbation")
for eps in (1e-4, 1e-3, 1e-2):
    He = HD + eps * VD
    ev, U = np.linalg.eigh(He)
    ov = np.einsum('ij,jk,ki->i', U.conj().T, PD, U).real
    idx = np.argsort(-ov)[:mD]
    e = np.sort(ev[idx])
    print(f"      eps = {eps:.1e}:  splitting = {e[-1] - e[0]:.6e}   ratio to 2*eps = "
          f"{(e[-1] - e[0]) / (2 * eps):.6f}  <- LINEAR")
contractD = [pmat(''.join('XYZ'[c] if k == q else 'I' for k in range(nD)))
             for q in range(nD) for c in range(3)]
witD = code_space_protection(gensD, nD, SD, ZB, 2)
He = HD + 1e-3 * VD
print()
_ = clause_report(pmat(ZB), He, pmat(XB), contractD,
                  "R = Zbar_B under H = H0 + 1e-3 Zbar_A  (contractible splitter)",
                  code_space_v=("NO -- no operator of weight <= 2 acts as a logical anticommuting "
                                "with R" if witD is None else f"YES: {witD}"))
print("    => a CONTRACTIBLE, weight-1 perturbation splits the degeneracy at FIRST order, and the")
print("       record in the protected sector survives (i),(ii),(iii),(iv) and the code-space form of")
print("       (v).  So even the REPAIRED side-step -- 'first-order split by a contractible")
print("       perturbation' -- does not by itself kill the record.  What matters is whether the")
print("       splitting resolves THE RECORD'S OWN logical direction.")

# ============================================================ THE EQUIVALENCE THAT DOES HOLD
print("\n" + "-" * 104)
print("THE STATEMENT THAT IS TRUE, AND IS WHAT EXHAUSTIVENESS ACTUALLY NEEDS.")
print("-" * 104)
print("""    Define, for a record R on a code space E0,
        d_R  :=  the minimum weight of an operator that acts on E0 as a logical operation
                 ANTICOMMUTING with R.
    THEOREM.  the following are equivalent:
        (1) clause (v) holds for R against operations of contractible weight <= w;
        (2) d_R > w;
        (3) no perturbation of weight <= w resolves R's own +-1 eigenspaces at first order,
            i.e. Phi(V) commutes with R|E0 for every V of weight <= w.
    (1)<=>(2) is the definition.  (2)<=>(3) is Knill-Laflamme applied to the R-graded code:
    a weight-<=w operator with Phi(V) not commuting with R|E0 IS an undetectable logical
    anticommuting with R, and conversely.""")
print("\n    MEASURED CHECK of (2)<=>(3) on all four carriers used above:")
for lab, gg, nn, SS, QQ, Rstr, wmax in [
        ("toric 2x2, R = Zbar2", gens, n, S, Q, Z2, 2),
        ("[[3,1,1]](x)[[5,1,3]], R = Zbar_B", gensD, nD, SD, QD, ZB, 3),
        ("[[3,1,1]](x)[[5,1,3]], R = Zbar_A", gensD, nD, SD, QD, ZA, 3)]:
    Rm = pmat(Rstr)
    R0 = QQ.conj().T @ Rm @ QQ
    for w in range(1, wmax + 1):
        worst = 0.0
        for supp in itertools.combinations(range(nn), w):
            for letters in itertools.product('XYZ', repeat=w):
                s = ['I'] * nn
                for k, ch in zip(supp, letters):
                    s[k] = ch
                Vm = pmat(''.join(s))
                A = phi(Vm, QQ)
                worst = max(worst, float(np.linalg.norm(A @ R0 - R0 @ A)))
        witness = code_space_protection(gg, nn, SS, Rstr, w)
        agree = (worst > 1e-9) == (witness is not None)
        print(f"      {lab:<36s} w={w}:  max ||[Phi(V), R|E0]|| = {worst:.3e}   "
              f"d_R <= w witness: {witness if witness else 'none'}   "
              f"{'AGREE' if agree else 'DISAGREE'}")

print("\n" + "=" * 104)
print("VERDICT ON ROW O-3")
print("=" * 104)
print("""  THE REGISTERED SIDE-STEP IS REFUTED AS WRITTEN.  "Any first-order-split degeneracy cannot carry
  a noise-surviving record" is false: CASE C splits at first order with all five clauses intact, and
  CASE D does it with a CONTRACTIBLE splitter.

  WHICH CLAUSE DIES, precisely:
     - keep R fixed and [V,R] != 0  ->  CLAUSE (ii), and it dies LINEARLY in eps even when the level
       splitting is of order eps^d.  (ii) and the splitting exponent are different clocks.
     - deform R to restore (ii), with the multiplet FULLY resolved  ->  CLAUSE (iii), for every R.
     - CLAUSE (v) dies from NEITHER.  It does not depend on eps or on V; it is a property of the pair
       (E0, R).  "The splitting destroys the protection" is a category error.

  THE REPAIR, still low cost but a DIFFERENT argument: replace the splitting-order criterion by the
  R-graded distance d_R, and use the equivalence (1)<=>(2)<=>(3) above.  Exhaustiveness then reads:
     a record survives contractible noise of weight <= w  IFF  d_R > w
  which is Knill-Laflamme, is an EQUIVALENCE rather than an implication, and does not mention the
  perturbation order at all.  COST ASSESSMENT: the registered "LOW" is right about the effort and
  wrong about the route.""")

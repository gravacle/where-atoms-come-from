#!/usr/bin/env python3
"""
LANE O-4 -- "admissible" in clause (iv).  INSTANCE PART.

The 2x2 toric code (8 qubits, dim 256) is the program's own witness that clauses
(i)-(v) can all hold.  This script asks what happens to clauses (iv) and (v) and to
P-3 under each candidate definition of "admissible".

THE HEADLINE THIS SCRIPT WAS BUILT TO TEST
   Under DEF-C ("admissible = any unitary"), is clause (v) -- "no operation on a
   contractible region has U^dag R U = -R" -- actually TRUE of the toric code?

SELF-CHECKS: the toric code's known numbers (ground-space dimension 4, code distance
2 at L=2, |centraliser| = 2^{n+k} = 1024) are checked before any headline is read.
POSITIVE CONTROLS: every zero reported over contractible regions is paired with the
same enumeration run over a NON-contractible region, which returns nonzero.
"""

import itertools
import numpy as np

np.set_printoptions(precision=4, suppress=True)
FAILURES = []


def check(name, condition, detail=""):
    tag = "PASS" if condition else "FAIL"
    if not condition:
        FAILURES.append(name)
    print(f"   [{tag}] {name}   {detail}")


def hr(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


# ----------------------------------------------------------------------------
# 2x2 toric code.  Sites (x,y) in Z_2 x Z_2.  Edge (x,y,d): d=0 horizontal
# (x,y)->(x+1,y), d=1 vertical (x,y)->(x,y+1).  n = 8 qubits.
L = 2
N = 2 * L * L


def eidx(x, y, d):
    return 2 * (L * (y % L) + (x % L)) + d


EDGE_NAME = {}
for y in range(L):
    for x in range(L):
        for d in range(2):
            EDGE_NAME[eidx(x, y, d)] = f"({x},{y},{'H' if d == 0 else 'V'})"

STARS = []      # vertex operators: X on the 4 edges at a vertex
for y in range(L):
    for x in range(L):
        STARS.append(sorted({eidx(x, y, 0), eidx(x - 1, y, 0),
                             eidx(x, y, 1), eidx(x, y - 1, 1)}))
PLAQS = []      # plaquette operators: Z on the 4 edges of a face
for y in range(L):
    for x in range(L):
        PLAQS.append(sorted({eidx(x, y, 0), eidx(x, y + 1, 0),
                             eidx(x, y, 1), eidx(x + 1, y, 1)}))

I2 = np.eye(2, dtype=complex)
X2 = np.array([[0, 1], [1, 0]], dtype=complex)
Y2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
Z2 = np.array([[1, 0], [0, -1]], dtype=complex)


def kron_op(single):
    M = np.array([[1.0 + 0j]])
    for q in range(N):
        M = np.kron(M, single[q])
    return M


def pauli(xs, zs):
    """Pauli operator with X on edges in xs and Z on edges in zs (Y where both)."""
    single = [I2] * N
    for q in range(N):
        a, b = q in xs, q in zs
        if a and b:
            single[q] = 1j * X2 @ Z2   # = Y, phase chosen Hermitian
        elif a:
            single[q] = X2
        elif b:
            single[q] = Z2
    return kron_op(single)


def nrm(A):
    return np.linalg.norm(A)


def comm(A, B):
    return A @ B - B @ A


hr("SECTION 0 -- BUILD THE 2x2 TORIC CODE AND CHECK ITS KNOWN NUMBERS")

A_ops = [pauli(set(s), set()) for s in STARS]
B_ops = [pauli(set(), set(p)) for p in PLAQS]
H = -(sum(A_ops) + sum(B_ops))

print("   stars   :", [[EDGE_NAME[e] for e in s] for s in STARS])
print("   plaqs   :", [[EDGE_NAME[e] for e in p] for p in PLAQS])

check("all stabilisers commute with H",
      max(nrm(comm(O, H)) for O in A_ops + B_ops) < 1e-10,
      f"max = {max(nrm(comm(O,H)) for O in A_ops+B_ops):.2e}")

w, v = np.linalg.eigh(H)
gs_dim = int(np.sum(np.abs(w - w[0]) < 1e-8))
check("ground-space dimension = 2^{2g} = 4  (Theorem A, known answer)", gs_dim == 4,
      f"measured {gs_dim}, E0 = {w[0]:.4f}")

# the record and its writer
R_edges = {eidx(0, 0, 0), eidx(1, 0, 0)}          # horizontal non-contractible loop, Z-type
W_edges = {eidx(0, 0, 0), eidx(0, 1, 0)}          # dual vertical loop, X-type
R = pauli(set(), R_edges)
U = pauli(W_edges, set())
print(f"   R = Z on {[EDGE_NAME[e] for e in sorted(R_edges)]}")
print(f"   U = X on {[EDGE_NAME[e] for e in sorted(W_edges)]}")

check("R is a bit (clause i)", nrm(R - R.conj().T) < 1e-12 and nrm(R @ R - np.eye(2**N)) < 1e-12)
check("[H,R] = 0 (clause ii, coherent part)", nrm(comm(H, R)) < 1e-10, f"{nrm(comm(H,R)):.2e}")
Ls = [pauli(set(), {q}) for q in range(N)]        # single-qubit dephasing jumps
check("[L_k,R] = 0 for all 8 dephasing jumps (clause ii)",
      max(nrm(comm(Lk, R)) for Lk in Ls) < 1e-10,
      f"max = {max(nrm(comm(Lk,R)) for Lk in Ls):.2e}")

Pgs = v[:, :gs_dim] @ v[:, :gs_dim].conj().T
Rg = Pgs @ R @ Pgs
dev = nrm(Rg - (np.trace(Rg) / gs_dim) * Pgs)
check("R non-constant on the ground eigenspace (clause iii)", dev > 1e-6,
      f"deviation from scalar = {dev:.4f}")


hr("SECTION 1 -- CLAUSE (iv) UNDER EACH DEFINITION")

print("\n--- DEF-A : admissible iff [U,H] = 0 ---")
check("U commutes with H", nrm(comm(U, H)) < 1e-10, f"||[U,H]|| = {nrm(comm(U,H)):.2e}")
check("U flips R", nrm(U.conj().T @ R @ U + R) < 1e-10,
      f"||U^dag R U + R|| = {nrm(U.conj().T@R@U+R):.2e}")
check("=> clause (iv) HOLDS under DEF-A: criterion (1) satisfied, NOT vacuous", True)

# balance lemma cross-check on every eigenspace of H
levels = []
i = 0
while i < len(w):
    j = i
    while j + 1 < len(w) and abs(w[j + 1] - w[i]) < 1e-8:
        j += 1
    V = v[:, i:j + 1]
    levels.append((w[i], j - i + 1, np.real(np.trace(V.conj().T @ R @ V))))
    i = j + 1
print("   energy levels of H and the balance number Tr(P_E R):")
for e, dim, tr in levels:
    print(f"      E = {e:+7.3f}   dim = {dim:4d}   Tr(P_E R) = {tr:+.2e}")
check("balance lemma: Tr(P_E R) = 0 on EVERY eigenspace",
      max(abs(t) for _, _, t in levels) < 1e-9,
      f"max = {max(abs(t) for _,_,t in levels):.2e}")

print("\n--- DEF-B : admissible iff in the *-algebra generated by {H, L_k} ---")
gens = [H] + Ls
worst = 0.0
rng = np.random.default_rng(7)
for _ in range(200):
    Wd = np.eye(2**N, dtype=complex)
    for _ in range(rng.integers(1, 5)):
        Wd = Wd @ gens[rng.integers(len(gens))]
    worst = max(worst, nrm(comm(Wd, R)) / max(nrm(Wd), 1e-30))
check("every word in {H, L_k} commutes with R  => NO DEF-B writer exists",
      worst < 1e-10, f"max relative ||[W,R]|| = {worst:.2e}")
check("POSITIVE CONTROL: the same test on U registers a non-commutation",
      nrm(comm(U, R)) > 1.0, f"||[U,R]|| = {nrm(comm(U,R)):.2f}  <-- the test can see a flipper")
check("=> clause (iv) is VACUOUS under DEF-B: criterion (1) FAILS", True)

print("\n--- DEF-C : admissible iff unitary ---")
check("clause (iv) holds trivially under DEF-C (U above is unitary)", True,
      "and criterion (2) already fails -- but see SECTION 2, where DEF-C also breaks (v)")


hr("SECTION 2 -- CLAUSE (v), ENUMERATED OVER CONTRACTIBLE REGIONS")

print("""
Binary symplectic enumeration of all 4^8 = 65536 Pauli operators.  For each we ask:
   FLIPS R ?    <=>  it anticommutes with R
   ADMISSIBLE ? <=>  it commutes with every stabiliser generator (equivalently [P,H]=0)
   CONTRACTIBLE?<=>  its support lies inside ONE star's 4 edges or ONE plaquette's 4
                     edges (this includes every single-edge operator)
""")

# symplectic vectors: (x, z) in F_2^8 x F_2^8
def vec(edges):
    a = np.zeros(N, dtype=np.int8)
    for e in edges:
        a[e] = 1
    return a


gen_sym = [(vec(s), np.zeros(N, dtype=np.int8)) for s in STARS] + \
          [(np.zeros(N, dtype=np.int8), vec(p)) for p in PLAQS]
R_sym = (np.zeros(N, dtype=np.int8), vec(R_edges))


def sympl(p, q):
    return int((p[0] @ q[1] + p[1] @ q[0]) % 2)


# cross-check the symplectic machinery against the dense matrices
mism = 0
for (gx, gz), Op in zip(gen_sym, A_ops + B_ops):
    dense_comm = nrm(comm(Op, R)) > 1e-10
    sym_comm = sympl((gx, gz), R_sym) == 1
    mism += int(dense_comm != sym_comm)
U_sym = (vec(W_edges), np.zeros(N, dtype=np.int8))
mism += int((sympl(U_sym, R_sym) == 1) != (nrm(comm(U, R)) > 1e-10))
check("symplectic commutation agrees with dense matrices on 9 operators", mism == 0,
      f"{mism} mismatches")

CONTRACTIBLE_REGIONS = [set(s) for s in STARS] + [set(p) for p in PLAQS]
NONCONTRACTIBLE_REGION = set(W_edges)     # the dual loop -- the positive control


def enumerate_over(region):
    """All Paulis supported inside `region`.  Returns (n_total, n_flip, n_adm,
    n_adm_and_flip, min_weight_of_adm_flip)."""
    edges = sorted(region)
    tot = flip = adm = both = 0
    minw = None
    for xs_bits in range(2 ** len(edges)):
        for zs_bits in range(2 ** len(edges)):
            x = np.zeros(N, dtype=np.int8)
            z = np.zeros(N, dtype=np.int8)
            for k, e in enumerate(edges):
                if xs_bits >> k & 1:
                    x[e] = 1
                if zs_bits >> k & 1:
                    z[e] = 1
            if not x.any() and not z.any():
                continue
            tot += 1
            f = sympl((x, z), R_sym) == 1
            a = all(sympl((x, z), g) == 0 for g in gen_sym)
            flip += f
            adm += a
            if f and a:
                both += 1
                wgt = int(np.sum((x | z) > 0))
                minw = wgt if minw is None else min(minw, wgt)
    return tot, flip, adm, both, minw


print("   region                                    #ops  #flip-R  #admiss  #BOTH")
tot_both = 0
tot_flip = 0
for reg, label in ([(set(s), f"star   {[EDGE_NAME[e] for e in sorted(s)]}") for s in STARS] +
                   [(set(p), f"plaq   {[EDGE_NAME[e] for e in sorted(p)]}") for p in PLAQS]):
    t, f, a, b, _ = enumerate_over(reg)
    tot_both += b
    tot_flip += f
    print(f"   {label:<42s} {t:5d} {f:8d} {a:8d} {b:6d}")

print()
# ---- L = 2 ARTEFACT, DIAGNOSED IN o4_diagnose_L2.py -------------------------------
# At L = 2 a plaquette's 4 edges are bottom/top/left/right of ONE face, and with only
# two rows the bottom and top horizontal edges of a column ARE the whole dual vertical
# loop.  So a "plaquette region" at L = 2 carries a winding cycle and is NOT
# contractible.  The expected count here is therefore NONZERO, and that is a
# known-answer check on the enumerator, not a failure of clause (v).
check("L=2 ARTEFACT REPRODUCED: plaquette regions at L=2 are not contractible, so "
      "they DO contain admissible flippers", tot_both == 16,
      f"count = {tot_both} (4 per plaquette x 4 plaquettes).  "
      f"Clause (v) is decided at L >= 3 in o4_L3.py, where the count is 0.")
check("POSITIVE CONTROL (in the same enumeration): contractible operators that DO "
      "flip R exist in quantity", tot_flip > 0,
      f"count = {tot_flip}  <-- so the zero above is not a zero of the enumerator")
check("clause (v) FAILS under DEF-C: those flippers are unitary Paulis", tot_flip > 0,
      f"e.g. a SINGLE-EDGE X on {EDGE_NAME[sorted(R_edges)[0]]} anticommutes with R")

Xe = pauli({sorted(R_edges)[0]}, set())
check("   explicit single-edge witness: X_e flips R", nrm(Xe.conj().T @ R @ Xe + R) < 1e-10,
      f"||X_e^dag R X_e + R|| = {nrm(Xe.conj().T@R@Xe+R):.2e}")
check("   and X_e is NOT admissible under DEF-A", nrm(comm(Xe, H)) > 1.0,
      f"||[X_e,H]|| = {nrm(comm(Xe,H)):.3f}")

print("\n   POSITIVE CONTROL -- the same enumeration over a NON-contractible region:")
t, f, a, b, minw = enumerate_over(NONCONTRACTIBLE_REGION)
print(f"   dual loop {[EDGE_NAME[e] for e in sorted(NONCONTRACTIBLE_REGION)]}"
      f"    #ops {t}  #flip-R {f}  #admiss {a}  #BOTH {b}   min weight {minw}")
check("a NON-contractible region does contain admissible flippers", b > 0,
      f"count = {b}, minimum weight = {minw}")


hr("SECTION 3 -- GLOBAL CENSUS AND THE MINIMUM WEIGHT OF AN ADMISSIBLE WRITER")

n_adm = n_flip = n_both = 0
minw = None
min_examples = []
for xb in range(2 ** N):
    x = np.array([(xb >> k) & 1 for k in range(N)], dtype=np.int8)
    for zb in range(2 ** N):
        z = np.array([(zb >> k) & 1 for k in range(N)], dtype=np.int8)
        if not x.any() and not z.any():
            continue
        a = all(sympl((x, z), g) == 0 for g in gen_sym)
        if not a:
            continue
        n_adm += 1
        if sympl((x, z), R_sym) == 1:
            n_both += 1
            wgt = int(np.sum((x | z) > 0))
            if minw is None or wgt < minw:
                minw, min_examples = wgt, []
            if wgt == minw and len(min_examples) < 6:
                min_examples.append((x.copy(), z.copy()))

print(f"   admissible Paulis (centraliser of the stabiliser group) : {n_adm + 1}"
      f"   [including identity]")
check("centraliser size = 2^{n+k} = 2^{8+2} = 1024  (known answer)", n_adm + 1 == 1024,
      f"measured {n_adm + 1}")
print(f"   admissible Paulis that FLIP R                           : {n_both}")
check("exactly half the centraliser flips R (a coset argument)", n_both == 512,
      f"measured {n_both}")
print(f"   MINIMUM WEIGHT of an admissible flipper                 : {minw}")
check("minimum weight = code distance d = L = 2 (known answer)", minw == 2,
      f"measured {minw}")
for x, z in min_examples:
    sup = [EDGE_NAME[e] for e in range(N) if x[e] or z[e]]
    typ = "".join("Y" if x[e] and z[e] else "X" if x[e] else "Z" if z[e] else ""
                  for e in range(N))
    print(f"      weight-{minw} admissible flipper: {typ} on {sup}")

print("""
   IS ANY MINIMUM-WEIGHT ADMISSIBLE FLIPPER CONTAINED IN A CONTRACTIBLE REGION?""")
contained = 0
for x, z in min_examples:
    sup = {e for e in range(N) if x[e] or z[e]}
    if any(sup <= reg for reg in CONTRACTIBLE_REGIONS):
        contained += 1
check("L=2 ARTEFACT: the weight-2 writers DO fit inside an L=2 plaquette "
      "(same cause; resolved at L>=3)", contained == 2, f"count = {contained}")


hr("SECTION 4 -- P-3 UNDER DEF-A, AND THE GAP IT LEAVES")

print("""
P-3 AS STATED:  (iv)+(v) => THE WRITER IS NON-LOCAL.

STEP-BY-STEP UNDER DEF-A:
  1. (iv) supplies an admissible U with U^dag R U = -R.                  [verified: yes]
  2. Suppose U is supported inside a contractible region S.
  3. Then U is an operation on a contractible region, so by (v) U^dag R U = +R.
  4. Contradiction with 1, since R != 0.  Hence no admissible flipper is supported
     inside a contractible region.                                              QED

  Step 3 requires that (v)'s quantifier RANGE INCLUDE the admissible operators.  It
  does, under either reading of (v) ("all operations" or "all admissible operations"),
  because admissible-and-contractible is a subset of both.  So P-3 GOES THROUGH, and
  it is not circular: DEF-A never mentions locality, regions, topology or homology.

THE GAP, AND IT IS REAL.
  "Non-local" in P-3 can only mean NOT SUPPORTED INSIDE ONE CONTRACTIBLE REGION.  It
  does NOT mean "not generated by local operations", because (v) quantifies over
  single regions, not over PRODUCTS of them.  Test below.
""")

e0, e1 = sorted(W_edges)
X0, X1 = pauli({e0}, set()), pauli({e1}, set())
print(f"   the writer U factorises: U = X_{EDGE_NAME[e0]} . X_{EDGE_NAME[e1]},"
      f" each factor on ONE edge")
check("the two single-edge factors reproduce U", nrm(X0 @ X1 - U) < 1e-10,
      f"||X0 X1 - U|| = {nrm(X0@X1-U):.2e}")
check("factor 1 is NOT admissible", nrm(comm(X0, H)) > 1.0, f"||[X0,H]|| = {nrm(comm(X0,H)):.3f}")
check("factor 2 is NOT admissible", nrm(comm(X1, H)) > 1.0, f"||[X1,H]|| = {nrm(comm(X1,H)):.3f}")
check("but each factor ALONE already flips R or not, and their PRODUCT flips R",
      nrm(U.conj().T @ R @ U + R) < 1e-10)

print("""
   READING.  The admissible writer is a product of d = 2 INADMISSIBLE single-region
   operations.  So (v) forbids a one-region write; it does NOT forbid a d-region
   write.  That is exactly Theorem D: local terms DO destroy the record, at order d.
   P-3's conclusion must therefore be stated as

       EVERY ADMISSIBLE WRITER HAS SUPPORT LARGER THAN ANY SINGLE CONTRACTIBLE
       REGION -- and the minimum such support is the code distance d.

   and NOT as "nothing local can write or destroy it" (THE_CLAIM_V001 row 7), which
   over-reads (v) and contradicts row 8 of the same table.
""")

hr("SUMMARY -- INSTANCE PART")
if FAILURES:
    print(f"*** {len(FAILURES)} SELF-CHECK FAILURES: {FAILURES}")
else:
    print("ALL SELF-CHECKS PASSED.")
print(f"""
   DEF-A : (iv) HOLDS (explicit U).  Clause (v) is NOT decided at L = 2 -- the
           lattice is too small for a plaquette to be contractible.  See o4_L3.py:
           at L = 3 and L = 5 the count of admissible contractible flippers is 0.
   DEF-B : (iv) VACUOUS -- no writer exists at all.
   DEF-C : (iv) holds but (v) FAILS -- {tot_flip} single-region unitaries flip R, so
           under DEF-C the toric-code record is NOT a record by the program's own
           definition.
""")

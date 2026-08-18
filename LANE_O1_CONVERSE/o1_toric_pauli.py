#!/usr/bin/env python3
"""
LANE O-1 -- THE CONVERSE ON THE PROGRAM'S OWN CARRIER.

The abstract counterexamples in o1_gap_a / o1_gap_b show the converse is false in general.
This file asks the sharper question: does the PROGRAM'S OWN FLAGSHIP CARRIER -- the toric
code of PROOF_V001 Theorems A-D -- satisfy clause (ii) once realistic local noise is
present?

  RESULT.  It does not.  With single-edge X and Z jump operators the commutant of
  alg{I, H, L_k, L_k^dag} is the SCALARS, so by clause (ii) the toric code has NO RECORD
  AT ALL -- while Theorem C's projector-level protection is untouched, exactly zero.

  So the exact commutator in clause (ii) and the projector condition in Theorem C are NOT
  the same condition, and the program has been reading a theorem about the second as
  support for the first.  ||[L,R]|| = 2.000 and ||P [L,R] P|| = 0.000 in the SAME system.

METHOD.  All jump operators here are Paulis, so the commutant of the generated algebra is
the span of the Paulis commuting with every generator, and the whole computation is exact
F_2 symplectic linear algebra -- no floating point, no tolerance.  The dense 256x256
cross-check at L=2 confirms it against explicit matrices.
"""

import itertools
import numpy as np

np.set_printoptions(precision=6, suppress=True, linewidth=170)


# ------------------------------------------------------------------ F_2 linear algebra
def rref2(M):
    M = M.copy() % 2
    rows, cols = M.shape
    piv, r = [], 0
    for c in range(cols):
        p = None
        for i in range(r, rows):
            if M[i, c]:
                p = i
                break
        if p is None:
            continue
        M[[r, p]] = M[[p, r]]
        for i in range(rows):
            if i != r and M[i, c]:
                M[i] ^= M[r]
        piv.append(c)
        r += 1
        if r == rows:
            break
    return M[:r], piv


def rank2(M):
    if M.size == 0:
        return 0
    return rref2(M)[0].shape[0]


def nullspace2(M):
    """Basis of {v : M v = 0 mod 2}, as rows."""
    rows, cols = M.shape
    R, piv = rref2(M)
    free = [c for c in range(cols) if c not in piv]
    basis = []
    for f in free:
        v = np.zeros(cols, dtype=np.uint8)
        v[f] = 1
        for i, c in enumerate(piv):
            v[c] = R[i, f]
        basis.append(v)
    return np.array(basis, dtype=np.uint8) if basis else np.zeros((0, cols), dtype=np.uint8)


def in_span2(v, B):
    if B.shape[0] == 0:
        return not v.any()
    return rank2(np.vstack([B, v])) == rank2(B)


# ------------------------------------------------------------------ toric code
class Toric:
    """L x L torus.  Edge index: h(x,y) = 2*(y*L+x), v(x,y) = 2*(y*L+x)+1."""

    def __init__(self, L):
        self.L = L
        self.n = 2 * L * L

    def h(self, x, y):
        return 2 * ((y % self.L) * self.L + (x % self.L))

    def v(self, x, y):
        return 2 * ((y % self.L) * self.L + (x % self.L)) + 1

    def vertex_ops(self):
        """A_v = X on the four edges incident to v."""
        out = []
        for y in range(self.L):
            for x in range(self.L):
                xs = np.zeros(self.n, dtype=np.uint8)
                for e in (self.h(x, y), self.h(x - 1, y), self.v(x, y), self.v(x, y - 1)):
                    xs[e] ^= 1
                out.append(np.concatenate([xs, np.zeros(self.n, dtype=np.uint8)]))
        return np.array(out)

    def plaquette_ops(self):
        """B_p = Z on the four edges bounding p."""
        out = []
        for y in range(self.L):
            for x in range(self.L):
                zs = np.zeros(self.n, dtype=np.uint8)
                for e in (self.h(x, y), self.h(x, y + 1), self.v(x, y), self.v(x + 1, y)):
                    zs[e] ^= 1
                out.append(np.concatenate([np.zeros(self.n, dtype=np.uint8), zs]))
        return np.array(out)

    def stabilizers(self):
        return np.vstack([self.vertex_ops(), self.plaquette_ops()])

    def single_edge(self, e, kind):
        p = np.zeros(2 * self.n, dtype=np.uint8)
        if kind == "X":
            p[e] = 1
        else:
            p[self.n + e] = 1
        return p


def symplectic_matrix(gens, n):
    """Rows: for each generator g, the linear functional v -> omega(g,v) over F_2."""
    rows = []
    for g in gens:
        gx, gz = g[:n], g[n:]
        rows.append(np.concatenate([gz, gx]) % 2)   # omega((x|z),(a|b)) = x.b + z.a
    return np.array(rows, dtype=np.uint8) if rows else np.zeros((0, 2 * n), dtype=np.uint8)


def commutant_and_records(T, noise):
    """Return (nu, n_logical) : nullity of the symplectic system, and the number of
    independent Pauli RECORDS (commutant elements acting non-trivially on the ground space).
    """
    n = T.n
    S = T.stabilizers()
    gens = np.vstack([S, noise]) if len(noise) else S
    M = symplectic_matrix(gens, n)
    W = nullspace2(M)                       # Paulis commuting with H and every L_k
    nu = W.shape[0]
    Sred, _ = rref2(S)                      # the stabilizer subgroup
    # records = W modulo (W intersect stabilizer group)
    inter = 0
    if nu:
        stack = np.vstack([Sred, W]) if Sred.shape[0] else W
        inter = rank2(Sred) + nu - rank2(stack)
    return nu, nu - inter


def show(t):
    print("\n" + "=" * 98)
    print(t)
    print("=" * 98)


# ------------------------------------------------------------------ dense cross-check
I2 = np.eye(2, dtype=complex)
Xp = np.array([[0, 1], [1, 0]], dtype=complex)
Zp = np.array([[1, 0], [0, -1]], dtype=complex)


def pauli_dense(vec, n):
    xs, zs = vec[:n], vec[n:]
    op = np.array([[1.0 + 0j]])
    for i in range(n):
        m = np.eye(2, dtype=complex)
        if xs[i]:
            m = m @ Xp
        if zs[i]:
            m = m @ Zp
        op = np.kron(op, m)
    return op


def main():
    print("O-1 ON THE PROGRAM'S OWN CARRIER -- TORIC CODE + PAULI NOISE (exact F_2)")

    # ---------------------------------------------------------------- self-checks
    show("SELF-CHECKS -- every expected value known independently of this code")
    ok = True

    def chk(label, got, want):
        nonlocal ok
        g = (got == want)
        ok = ok and g
        print(f"  [{'PASS' if g else 'FAIL'}] {label:64s} got {got:5d}  want {want:5d}")

    for L in (2, 3, 4):
        T = Toric(L)
        chk(f"L={L}: independent stabilizer generators = 2L^2 - 2",
            rank2(T.stabilizers()), 2 * L * L - 2)
        nu, k = commutant_and_records(T, [])
        chk(f"L={L}: centralizer dimension = 2n - (2L^2-2)", nu, 2 * T.n - (2 * L * L - 2))
        chk(f"L={L}: independent logical Paulis with no noise = 4 (2 logical qubits)", k, 4)
        # every stabilizer commutes with every other
        S = T.stabilizers()
        M = symplectic_matrix(S, T.n)
        bad = int(((M @ S.T) % 2).sum())
        chk(f"L={L}: stabilizer generators mutually commute (violations)", bad, 0)

    # dense check at L=2: ground space dimension 4 = 2^{2g}
    T2 = Toric(2)
    S2 = T2.stabilizers()
    H = np.zeros((2 ** T2.n, 2 ** T2.n), dtype=complex)
    for s in S2:
        H -= pauli_dense(s, T2.n)
    w = np.linalg.eigvalsh(H)
    gdim = int(np.sum(np.abs(w - w[0]) < 1e-8))
    chk("L=2 dense: ground-space dimension (Theorem A predicts 2^{2g} = 4)", gdim, 4)
    print(f"\n  SELF-CHECKS: {'ALL PASS' if ok else 'SOME FAILED -- STOP'}")
    if not ok:
        raise SystemExit(1)

    # ---------------------------------------------------------------- the null
    show("EXPERIMENT 1 -- CLAUSE (ii) UNDER LOCAL NOISE.  Records surviving each noise set.")
    print(f"  {'L':>2} {'n=edges':>8}  {'noise model':38s} {'dim commutant':>14} {'records':>8}")
    for L in (2, 3, 4):
        T = Toric(L)
        models = [
            ("none (positive control)", []),
            ("dephasing: Z on every edge", [T.single_edge(e, "Z") for e in range(T.n)]),
            ("bit flip: X on every edge", [T.single_edge(e, "X") for e in range(T.n)]),
            ("X and Z on ONE edge", [T.single_edge(0, "X"), T.single_edge(0, "Z")]),
            ("X and Z on EVERY edge", [T.single_edge(e, k) for e in range(T.n) for k in "XZ"]),
        ]
        for name, noise in models:
            nu, k = commutant_and_records(T, noise)
            print(f"  {L:2d} {T.n:8d}  {name:38s} {2**nu if nu < 60 else -1:>14} {k:8d}")
        print()
    print("  READING: the toric code satisfies clause (ii) under PURE DEPHASING (2 records")
    print("  survive: the Z-type logicals) and under PURE BIT-FLIP (2 records: the X-type),")
    print("  but under BOTH -- i.e. under generic single-qubit noise -- the commutant of the")
    print("  *-algebra is the scalars and the toric code HAS NO RECORD by clause (ii).")

    # ---------------------------------------------------------------- how much noise
    show("EXPERIMENT 2 -- HOW MUCH LOCAL NOISE DOES IT TAKE?  Full (X,Z) noise on a random\n"
         "                fraction f of edges; records surviving, averaged over 200 draws.")
    rng = np.random.default_rng(11)
    for L in (3, 4):
        T = Toric(L)
        print(f"\n  L = {L}, n = {T.n} edges, code distance d = L = {L}")
        print(f"  {'noisy edges':>12} {'f':>6}  {'mean records':>13}  {'P(some record)':>15}")
        for ne in range(0, T.n + 1, max(1, T.n // 10)):
            tot, alive = 0, 0
            for _ in range(200):
                idx = rng.choice(T.n, size=ne, replace=False)
                noise = [T.single_edge(int(e), kk) for e in idx for kk in "XZ"]
                nu, k = commutant_and_records(T, noise)
                tot += k
                alive += (k > 0)
            print(f"  {ne:12d} {ne/T.n:6.2f}  {tot/200:13.3f}  {alive/200:15.3f}")
    print("\n  Records die well before the noise covers the lattice: a set of noisy edges")
    print("  meeting every non-contractible cycle AND every non-contractible cocycle is")
    print("  enough, and such a set has size O(L), not O(L^2).")

    # ---------------------------------------------------------------- the two conditions
    show("EXPERIMENT 3 -- THE DECISIVE COMPARISON.  Clause (ii) and Theorem C are different\n"
         "                conditions, and they disagree on the program's own carrier.")
    T = Toric(2)
    n = T.n
    S = T.stabilizers()
    Hd = np.zeros((2 ** n, 2 ** n), dtype=complex)
    for s in S:
        Hd -= pauli_dense(s, n)
    wv, vv = np.linalg.eigh(Hd)
    gs = vv[:, np.abs(wv - wv[0]) < 1e-8]
    P = gs @ gs.conj().T
    print(f"  L=2 toric code, {n} qubits, dim {2**n}.  Ground space dimension {gs.shape[1]}.")

    # a Z-type logical: a non-contractible cycle of Z's -- the horizontal loop y=0
    zc = np.zeros(2 * n, dtype=np.uint8)
    for x in range(T.L):
        zc[n + T.h(x, 0)] = 1
    R = pauli_dense(zc, n)
    print(f"  R = Z-loop on the y=0 horizontal cycle.  R=R^dag: "
          f"{np.linalg.norm(R - R.conj().T):.2e}   R^2=I: {np.linalg.norm(R@R - np.eye(2**n)):.2e}")
    print(f"  [H,R] = {np.linalg.norm(Hd @ R - R @ Hd):.2e}")
    onP = gs.conj().T @ R @ gs
    print(f"  (iii) R is non-constant on the ground eigenspace: "
          f"||R|_gs - (tr/4) I|| = "
          f"{np.linalg.norm(onP - np.trace(onP)/4*np.eye(4)):.4f}   eigenvalues "
          f"{np.round(np.linalg.eigvalsh(onP),6)}")

    print("  (all norms below are SPECTRAL norms, so a full violation reads 2.0000)")
    print(f"\n  {'local jump operator L':28s} {'||[L,R]||':>12}  {'||P [L,R] P||':>15}  verdict")
    for e in range(4):
        for kind in "XZ":
            Lp = pauli_dense(T.single_edge(e, kind), n)
            c = np.linalg.norm(Lp @ R - R @ Lp, 2)     # SPECTRAL norm
            pc = np.linalg.norm(P @ (Lp @ R - R @ Lp) @ P, 2)
            verdict = "clause (ii) OK" if c < 1e-9 else \
                      ("CLAUSE (ii) FAILS, Thm C holds" if pc < 1e-9 else "both fail")
            print(f"  {kind}_e  on edge {e:2d}{'':14s} {c:12.4f}  {pc:15.3e}  {verdict}")

    print("\n  THE NUMBER: for the single-edge X jump operators lying ON the record's cycle,")
    print("  ||[L,R]|| = 2.0000 exactly -- clause (ii) is violated at O(1), not at O(eps) --")
    print("  while ||P [L,R] P|| = 0 to machine precision, which is Theorem C.")
    print("\n  CONSEQUENCE FOR O-1.  Clause (ii) is an operator identity; Theorem C is a")
    print("  statement about the compression to the record space.  Every 'durability' result")
    print("  the program has PROVED is of the second kind.  The converse of the degeneracy")
    print("  proposition cannot be closed while clause (ii) is stated as an exact operator")
    print("  commutator, because the program's own carrier fails it.")


if __name__ == "__main__":
    main()

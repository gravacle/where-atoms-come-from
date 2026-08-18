#!/usr/bin/env python3
"""
LANE O-1 -- the one step of Theorem O1-A that is quoted from the literature rather than
proved here, checked numerically so it is not taken on trust.

STEP QUOTED:  for a von Neumann algebra A on a finite-dimensional space and a projection
P in A, the compression map pi : A' -> B(PH), pi(X) = X|_{PH}, has image equal to the
commutant of P A P inside B(PH).

CONSEQUENCE USED:  a record non-constant on the eigenspace E_lam exists  <=>  the
commutant of P A P in B(E_lam) is non-trivial  <=>  P A P is a proper subalgebra of
B(E_lam).

Checked below on random systems with structured degeneracies: dim pi(A') is compared
against dim (PAP)' computed independently.  Both sides are computed by different code
paths (one compresses the commutant, the other takes the commutant of the compression).
"""

import numpy as np
from o1_core import (commutant_basis, star_closure, generated_algebra, eigenspaces, span_dim)

np.set_printoptions(precision=6, suppress=True, linewidth=170)


def main():
    rng = np.random.default_rng(31337)
    print("CHECK: dim pi(A')  ==  dim (P A P)'   for P an eigenprojector of H (P is in A)\n")
    print(f"  {'D':>3} {'m':>3}  {'dim A':>6} {'dim A_prime':>12} {'dim pi(A_prime)':>16} "
          f"{'dim (PAP)_prime':>16}  {'dim PAP':>8} {'m^2':>5}  ok")
    allok = True
    trials = 0
    nontrivial = 0
    for D, blocks in ((4, [2, 2]), (6, [3, 3]), (6, [2, 4]), (6, [2, 2, 2]),
                      (8, [4, 4]), (8, [2, 3, 3]), (6, [6]), (8, [8])):
        for t in range(4):
            # H constant on a random grouping of the blocks; L block-diagonal generic
            H = np.zeros((D, D), dtype=complex)
            L = np.zeros((D, D), dtype=complex)
            o = 0
            energies = [int(rng.integers(0, 2)) for _ in blocks]
            for i, s in enumerate(blocks):
                H[o:o + s, o:o + s] = energies[i] * np.eye(s)
                B = rng.normal(size=(s, s)) + 1j * rng.normal(size=(s, s))
                L[o:o + s, o:o + s] = B + 5.0 * (i + 1) * np.eye(s)
                o += s
            A = generated_algebra([H, L], D)
            Ap = commutant_basis(star_closure([H, L]), D)[0]
            lam, cols = eigenspaces(H)[0]
            m = cols.shape[1]
            if m < 2:
                continue
            piAp = [cols.conj().T @ X @ cols for X in Ap]
            PAP = [cols.conj().T @ X @ cols for X in A]
            dim_pi = span_dim(piAp)
            # commutant of PAP inside B(E)
            dim_comm = commutant_basis(star_closure(PAP), m)[1]
            dim_PAP = span_dim(PAP)
            ok = (dim_pi == dim_comm)
            allok = allok and ok
            trials += 1
            nontrivial += int(dim_PAP < m * m)
            print(f"  {D:3d} {m:3d}  {len(A):6d} {len(Ap):12d} {dim_pi:16d} {dim_comm:16d}  "
                  f"{dim_PAP:8d} {m*m:5d}  {'PASS' if ok else 'FAIL'}")
    print(f"\n  {trials} systems, of which {nontrivial} have a reducible compression and "
          f"{trials - nontrivial} an irreducible one")
    print(f"  (so the check is exercised on BOTH outcomes, not only on the easy side)")
    print(f"\n  QUOTED STEP: {'VERIFIED on every system tested' if allok else 'FAILED'}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
ROW O-10, PART 2 -- WHY THE THRESHOLD IS d, AND WHAT THE WITNESS ACTUALLY IS.

Part 1 (o10_threshold.py) measured threshold == d on six codes at d = 1, 2, 3.
Part 2 does two things:

(A) STATES AND CHECKS THE GENERAL PROOF, so the claim is not left at d <= 3.

  THEOREM.  Let S be a stabiliser group on n qubits with -I not in S,
  rank_{F2}(S) = n-k, k >= 1, and code projector Pc = (1/|S|) sum_{s in S} s.
  For any Pauli E:
    (a) if E anticommutes with some s in S,  then Pc E Pc = 0            (proportional, c=0)
    (b) if E = lambda*s for some s in S,     then Pc E Pc = lambda Pc    (proportional)
    (c) otherwise (E in N(S) but E not in <phases, S>),
        Pc E Pc is NOT proportional to Pc.
  Hence  min{ wt(E) : Pc E Pc not prop. Pc }  =  min{ wt(E) : E in N(S)\<i,S> }  =  d.

  PROOF of (a): pick s in S with Es = -sE.  Then Pc E Pc = Pc E (s Pc) = -Pc s E Pc
                = -(Pc s) E Pc = -Pc E Pc, so Pc E Pc = 0.   [uses s Pc = Pc s = Pc]
  PROOF of (b): E Pc = lambda s Pc = lambda Pc, so Pc E Pc = lambda Pc.
  PROOF of (c): E commutes with all of S, so E Pc = Pc E and Pc E Pc = E Pc.
                Suppose E Pc = c Pc.  E is unitary, so Pc E^dag E Pc = |c|^2 Pc,
                i.e. Pc = |c|^2 Pc; Pc != 0 since k >= 1, so |c| = 1.
                Then S' = < S, c^{-1} E > has -I not in S' and, because the symplectic
                vector of E lies OUTSIDE the F2 row space of S, rank(S') = n-k+1.
                Its +1 eigenspace therefore has dimension 2^{k-1}.  But c^{-1}E acts as
                +1 on the whole 2^k-dimensional code space of S, which is contained in
                that +1 eigenspace: 2^k <= 2^{k-1}.  Contradiction.  []

  The proof is n-, k-, d-, CSS- and degeneracy-independent.  What Part 1 rules out is
  an IMPLEMENTATION artefact, and what it exhibits is that DEGENERACY does not lower the
  threshold (Shor carries weight-2 stabiliser elements while d = 3).

  Checked numerically here, on every code:
    C-a  every weight<d Pauli that anticommutes with some stabiliser has Pc E Pc == 0 exactly
    C-b  every weight<d Pauli in the stabiliser group has Pc E Pc == lambda Pc, |lambda| == 1
    C-c  the rank-halving step: rank(projector of <S, c^{-1}E>) == 2^{k-1} for the witness E

(B) IDENTIFIES THE WITNESS AS A RECORD OPERATOR.
  For the weight-d witness E, restrict E to the code space in an orthonormal basis:
  R = V^dag E V, V an isometry onto the code space.  We check against CORE_FRAMEWORK
  clause (i) BIT:   R = R^dag  and  R^2 = I,  and report the spectrum.
  A balanced +-1 split (tr R = 0) is what "distinguishes SAME-energy states" means at
  clause (iii).  This also explains the observed residual exactly:
      X(E) = ||Pc E Pc - c Pc||_F^2 = tr(R^dag R) - |tr R|^2 / 2^k = 2^k - |tr R|^2 / 2^k
  so X == 2^k is the arithmetic signature of a BALANCED bit.  Predicted before reading:
  [[4,2,2]] and [[8,2,2]] -> 4 ; the k=1 codes -> 2.
"""

import itertools
import numpy as np

from o10_threshold import (parse_pauli, pauli_str, weight, code_projector,
                           apply_left, toric_2x2, phase_vector)

TOL = 1e-10


def pauli_matrix(a, b, n):
    """HERMITIAN Pauli string.  X(a)Z(b) is XZ = -iY on shared positions, so multiply by
    i^popcount(a&b) to recover the genuine Hermitian Pauli.  Phase is irrelevant to Part 1's
    proportionality test but essential here, where we check R = R^dag."""
    dim = 1 << n
    M = np.zeros((dim, dim), dtype=np.complex128)
    idx = np.arange(dim, dtype=np.int64)
    ph = phase_vector(b, dim)
    M[idx ^ a, idx] = ph
    return (1j ** bin(a & b).count('1')) * M


def isometry(P, kdim):
    w, V = np.linalg.eigh(P)
    order = np.argsort(-w)
    w = w[order]
    V = V[:, order]
    assert np.allclose(w[:kdim], 1.0, atol=1e-9), w[:kdim]
    assert np.allclose(w[kdim:], 0.0, atol=1e-9), w[kdim:kdim + 3]
    return V[:, :kdim]


def rank_of_extended_group(gens, a, b, c, n):
    """rank check for proof step (c): projector of <S, c^{-1}E> must have rank 2^{k-1}."""
    dim = 1 << n
    P = np.eye(dim, dtype=np.complex128)
    for (ga, gb) in gens:
        P = 0.5 * (P + apply_left(ga, gb, P))
    E = pauli_matrix(a, b, n) / c
    P = 0.5 * (P + E @ P)
    return int(np.linalg.matrix_rank(P, tol=1e-9))


def analyse(name, gens_str, n, kbits, d):
    print('=' * 78)
    print('CODE %s   n=%d k=%d d=%d' % (name, n, kbits, d))
    print('=' * 78)
    gens = [parse_pauli(s) for s in gens_str]
    P = code_projector(gens, n)
    dim = 1 << n
    kdim = 2 ** kbits
    V = isometry(P, kdim)

    # --- symplectic membership test for stabiliser group
    basis, pivots = [], []
    for (ga, gb) in gens:
        cur = (ga << n) | gb
        for p, bv in zip(pivots, basis):
            if cur >> p & 1:
                cur ^= bv
        if cur:
            basis.append(cur)
            pivots.append(cur.bit_length() - 1)

    def in_S(a, b):
        cur = (a << n) | b
        for p, bv in zip(pivots, basis):
            if cur >> p & 1:
                cur ^= bv
        return cur == 0

    def commutes(a, b):
        return all((bin(a & gb).count('1') + bin(ga & b).count('1')) % 2 == 0
                   for (ga, gb) in gens)

    # --- C-a / C-b : classify every Pauli of weight < d
    worst_a = 0.0
    worst_b = 0.0
    n_a = n_b = 0
    lam_bad = 0
    for w in range(1, d):
        for support in itertools.combinations(range(n), w):
            for choice in itertools.product((1, 2, 3), repeat=w):
                a = b = 0
                for pos, ch in zip(support, choice):
                    bit = 1 << (n - 1 - pos)
                    if ch & 1:
                        a |= bit
                    if ch & 2:
                        b |= bit
                E = pauli_matrix(a, b, n)
                A = P @ E @ P
                if not commutes(a, b):
                    n_a += 1
                    worst_a = max(worst_a, np.max(np.abs(A)))
                else:
                    assert in_S(a, b), 'weight<d element of N(S)\\S would contradict d'
                    n_b += 1
                    lam = np.trace(A) / kdim
                    worst_b = max(worst_b, np.max(np.abs(A - lam * P)))
                    if abs(abs(lam) - 1.0) > 1e-9:
                        lam_bad += 1
    print('  C-a  %5d Paulis of weight<d anticommute with S : max|PcEPc| = %.3e -> %s'
          % (n_a, worst_a, 'PASS' if worst_a < TOL else 'FAIL'))
    print('  C-b  %5d Paulis of weight<d lie IN S           : max|PcEPc - lam Pc| = %.3e, |lam|!=1 count %d -> %s'
          % (n_b, worst_b, lam_bad, 'PASS' if (worst_b < TOL and lam_bad == 0) else 'FAIL'))
    if n_b == 0:
        print('       (no weight<d stabiliser elements: this code is NON-DEGENERATE below d)')

    # --- the weight-d witnesses
    wits = []
    for support in itertools.combinations(range(n), d):
        for choice in itertools.product((1, 2, 3), repeat=d):
            a = b = 0
            for pos, ch in zip(support, choice):
                bit = 1 << (n - 1 - pos)
                if ch & 1:
                    a |= bit
                if ch & 2:
                    b |= bit
            if commutes(a, b) and not in_S(a, b):
                wits.append((a, b))
    print('  weight-d elements of N(S)\\S found: %d' % len(wits))

    # analyse the first few, and the rank-halving step on the first
    shown = 0
    all_bit = True
    all_balanced = True
    for (a, b) in wits:
        E = pauli_matrix(a, b, n)
        R = V.conj().T @ E @ V
        herm = np.max(np.abs(R - R.conj().T))
        sq = np.max(np.abs(R @ R - np.eye(kdim)))
        ev = np.linalg.eigvalsh((R + R.conj().T) / 2) if herm < 1e-9 else None
        trR = np.trace(R)
        X = np.trace(R.conj().T @ R).real - abs(trR) ** 2 / kdim
        if herm >= 1e-9 or sq >= 1e-9:
            all_bit = False
        if abs(trR) > 1e-9:
            all_balanced = False
        if shown < 3:
            print('    witness %-12s  R=R^dag dev %.2e  R^2=I dev %.2e  tr R = %+.3f%+.3fi  spec %s  X = %.6f'
                  % (pauli_str(a, b, n), herm, sq, trR.real, trR.imag,
                     np.array2string(np.round(ev, 6), precision=3) if ev is not None else 'n/a', X))
            shown += 1
    print('  BIT clause (i) on ALL %d weight-d witnesses: R=R^dag and R^2=I -> %s'
          % (len(wits), 'PASS' if all_bit else 'FAIL'))
    print('  balanced (tr R = 0) on ALL witnesses -> %s ; predicted X = 2^k = %d'
          % ('PASS' if all_balanced else 'FAIL', kdim))

    # --- C-c rank halving
    a, b = wits[0]
    E = pauli_matrix(a, b, n)
    R = V.conj().T @ E @ V
    c = np.trace(R) / kdim
    # E does not act as a scalar, so pick c = +1 branch of the proof: use eigenvalue +1
    rk = rank_of_extended_group(gens, a, b, 1.0, n)
    print('  C-c  rank(projector of <S, E>) = %d ; 2^(k-1) = %d -> %s'
          % (rk, kdim // 2, 'PASS' if rk == kdim // 2 else 'FAIL'))
    print('       (the proof\'s contradiction made concrete: adjoining the witness HALVES the space,')
    print('        so the witness cannot have been acting as a scalar on the full code space)')
    print()
    return dict(name=name, all_bit=all_bit, all_balanced=all_balanced,
                ca=(worst_a < TOL), cb=(worst_b < TOL and lam_bad == 0),
                cc=(rk == kdim // 2), nwit=len(wits), degenerate=(n_b > 0))


def main():
    codes = [
        ('[[3,1,1]] repetition (CONTROL)', ['ZZI', 'IZZ'], 3, 1, 1),
        ('[[4,2,2]]', ['XXXX', 'ZZZZ'], 4, 2, 2),
        None,  # toric filled in below
        ('[[5,1,3]] perfect (NON-CSS)', ['XZZXI', 'IXZZX', 'XIXZZ', 'ZXIXZ'], 5, 1, 3),
        ('[[7,1,3]] Steane (CSS)', ['IIIXXXX', 'IXXIIXX', 'XIXIXIX',
                                    'IIIZZZZ', 'IZZIIZZ', 'ZIZIZIZ'], 7, 1, 3),
        ('[[9,1,3]] Shor (DEGENERATE)', ['ZZIIIIIII', 'IZZIIIIII',
                                         'IIIZZIIII', 'IIIIZZIII',
                                         'IIIIIIZZI', 'IIIIIIIZZ',
                                         'XXXXXXIII', 'IIIXXXXXX'], 9, 1, 3),
    ]
    tg, tn = toric_2x2()
    codes[2] = ('[[8,2,2]] toric 2x2', tg, tn, 2, 2)

    res = [analyse(*c) for c in codes]

    print('=' * 78)
    print('O-10 PART 2 SUMMARY')
    print('=' * 78)
    hdr = '%-34s %5s %5s %5s %8s %10s %6s' % ('code', 'C-a', 'C-b', 'C-c', 'BIT(i)', 'balanced', '#wit')
    print(hdr)
    print('-' * len(hdr))
    ok = True
    for r in res:
        ok = ok and r['ca'] and r['cb'] and r['cc'] and r['all_bit'] and r['all_balanced']
        print('%-34s %5s %5s %5s %8s %10s %6d'
              % (r['name'], 'PASS' if r['ca'] else 'FAIL', 'PASS' if r['cb'] else 'FAIL',
                 'PASS' if r['cc'] else 'FAIL', 'PASS' if r['all_bit'] else 'FAIL',
                 'PASS' if r['all_balanced'] else 'FAIL', r['nwit']))
    print()
    print('  degenerate-below-d codes in the set: %s'
          % ', '.join(r['name'] for r in res if r['degenerate']))
    print()
    print('VERDICT PART 2: %s' % ('every step of the general proof verified numerically on every code; '
                                  'each weight-d witness IS a balanced record operator (R=R^dag, R^2=I, tr R=0)'
                                  if ok else 'A STEP FAILED -- SEE TABLE'))


if __name__ == '__main__':
    main()

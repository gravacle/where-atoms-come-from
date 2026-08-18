#!/usr/bin/env python3
"""
ROW O-10, PART 3 -- THE THRESHOLD IS ABOUT COUPLINGS, NOT ABOUT PAULIS.

F-13 states the RFP needs "a coupling of weight >= d".  Parts 1 and 2 measured and proved
the threshold for PAULI operators.  A coupling in a Hamiltonian or a Lindblad operator is
a general operator, not a Pauli.  This part closes that gap.

COROLLARY.  Let O be ANY operator (not necessarily Hermitian, not necessarily a Pauli)
whose Pauli expansion O = sum_E c_E E contains only terms of weight < d.  Then
Pc O Pc = lambda Pc for a scalar lambda.
  PROOF.  Linearity: Pc O Pc = sum_E c_E (Pc E Pc).  By the theorem of Part 2, each
  Pc E Pc with wt(E) < d equals lambda_E Pc (with lambda_E = 0 when E anticommutes with
  some stabiliser).  So Pc O Pc = (sum_E c_E lambda_E) Pc.   []
In particular EVERY operator supported on fewer than d qubits qualifies, since its Pauli
expansion lives entirely on that support.  So no coupling of weight < d can distinguish
code states -- the threshold is a property of the coupling's WEIGHT, not of its being Pauli.

TEST HERE: random *dense* Hermitian and random non-Hermitian operators supported on
(d-1) qubits, and random Pauli-sparse operators with all terms of weight < d.
  SELF-CHECK: residual ||Pc O Pc - lambda Pc||_F must be 0 to machine precision.
  POSITIVE CONTROL: the identical random construction on d qubits (one more qubit)
  must produce a NONZERO residual for at least one draw -- otherwise the null above is
  uncontrolled.  We report how many of the draws fire.
"""

import itertools
import numpy as np

from o10_threshold import parse_pauli, pauli_str, code_projector, toric_2x2
from o10_proof_check import pauli_matrix, toric_2x2 as _t  # noqa: F401

RNG = np.random.default_rng(20260818)


def embed(op_small, support, n):
    """Embed a 2^m x 2^m operator acting on qubits `support` into the full 2^n space."""
    m = len(support)
    dim = 1 << n
    # build via Pauli expansion on the support -- exact and avoids index gymnastics
    paulis1 = {
        'I': np.eye(2, dtype=np.complex128),
        'X': np.array([[0, 1], [1, 0]], dtype=np.complex128),
        'Y': np.array([[0, -1j], [1j, 0]], dtype=np.complex128),
        'Z': np.array([[1, 0], [0, -1]], dtype=np.complex128),
    }
    out = np.zeros((dim, dim), dtype=np.complex128)
    for labels in itertools.product('IXYZ', repeat=m):
        B = paulis1[labels[0]]
        for L in labels[1:]:
            B = np.kron(B, paulis1[L])
        coeff = np.trace(B.conj().T @ op_small) / (2 ** m)
        if abs(coeff) < 1e-15:
            continue
        s = ['I'] * n
        for q, L in zip(support, labels):
            s[q] = L
        a, b = parse_pauli(''.join(s))
        out += coeff * pauli_matrix(a, b, n)
    return out


def resid(P, O, kdim):
    A = P @ O @ P
    lam = np.trace(A) / kdim
    return np.max(np.abs(A - lam * P))


def run(name, gens_str, n, kbits, d, ndraws=40):
    print('=' * 78)
    print('CODE %s   n=%d k=%d d=%d' % (name, n, kbits, d))
    print('=' * 78)
    gens = [parse_pauli(s) for s in gens_str]
    P = code_projector(gens, n)
    kdim = 2 ** kbits
    supports_lo = list(itertools.combinations(range(n), d - 1))
    supports_hi = list(itertools.combinations(range(n), d))

    def draw(m, herm):
        M = RNG.normal(size=(1 << m, 1 << m)) + 1j * RNG.normal(size=(1 << m, 1 << m))
        return (M + M.conj().T) / 2 if herm else M

    worst_lo = 0.0
    n_lo = 0
    if d - 1 >= 1:
        for herm in (True, False):
            for _ in range(ndraws):
                sup = supports_lo[RNG.integers(len(supports_lo))]
                O = embed(draw(d - 1, herm), sup, n)
                worst_lo = max(worst_lo, resid(P, O, kdim))
                n_lo += 1
        print('  weight-(d-1)=%d random couplings, %d draws (Hermitian and non-Hermitian):'
              % (d - 1, n_lo))
        print('    max ||Pc O Pc - lambda Pc||_max = %.3e  -> %s'
              % (worst_lo, 'PASS (all proportional)' if worst_lo < 1e-9 else 'FAIL'))
    else:
        print('  d=1: no weight-(d-1) operators other than multiples of I; step vacuous.')

    # POSITIVE CONTROL at weight d
    fired = 0
    worst_hi = 0.0
    for herm in (True, False):
        for _ in range(ndraws):
            sup = supports_hi[RNG.integers(len(supports_hi))]
            O = embed(draw(d, herm), sup, n)
            r = resid(P, O, kdim)
            worst_hi = max(worst_hi, r)
            if r > 1e-9:
                fired += 1
    print('  POSITIVE CONTROL, weight-d=%d random couplings, %d draws: %d FIRED, max residual %.6f -> %s'
          % (d, 2 * ndraws, fired, worst_hi, 'PASS' if fired > 0 else 'FAIL (uncontrolled null above)'))

    # also: a Pauli-sparse random operator with mixed supports, all weights < d
    if d - 1 >= 1:
        allP = []
        for w in range(1, d):
            for sup in itertools.combinations(range(n), w):
                for ch in itertools.product('XYZ', repeat=w):
                    s = ['I'] * n
                    for q, L in zip(sup, ch):
                        s[q] = L
                    allP.append(''.join(s))
        worst_mix = 0.0
        for _ in range(20):
            O = np.zeros((1 << n, 1 << n), dtype=np.complex128)
            picks = RNG.choice(len(allP), size=min(12, len(allP)), replace=False)
            for i in picks:
                a, b = parse_pauli(allP[i])
                O = O + (RNG.normal() + 1j * RNG.normal()) * pauli_matrix(a, b, n)
            worst_mix = max(worst_mix, resid(P, O, kdim))
        print('  mixed-support Pauli-sparse couplings, all terms weight<d, 20 draws:')
        print('    max residual = %.3e -> %s' % (worst_mix, 'PASS' if worst_mix < 1e-9 else 'FAIL'))
    else:
        worst_mix = 0.0
    print()
    return dict(name=name, lo=worst_lo, fired=fired, hi=worst_hi, mix=worst_mix, d=d)


def main():
    tg, tn = toric_2x2()
    codes = [
        ('[[3,1,1]] repetition (CONTROL)', ['ZZI', 'IZZ'], 3, 1, 1),
        ('[[4,2,2]]', ['XXXX', 'ZZZZ'], 4, 2, 2),
        ('[[8,2,2]] toric 2x2', tg, tn, 2, 2),
        ('[[5,1,3]] perfect (NON-CSS)', ['XZZXI', 'IXZZX', 'XIXZZ', 'ZXIXZ'], 5, 1, 3),
        ('[[7,1,3]] Steane (CSS)', ['IIIXXXX', 'IXXIIXX', 'XIXIXIX',
                                    'IIIZZZZ', 'IZZIIZZ', 'ZIZIZIZ'], 7, 1, 3),
        ('[[9,1,3]] Shor (DEGENERATE)', ['ZZIIIIIII', 'IZZIIIIII',
                                         'IIIZZIIII', 'IIIIZZIII',
                                         'IIIIIIZZI', 'IIIIIIIZZ',
                                         'XXXXXXIII', 'IIIXXXXXX'], 9, 1, 3),
    ]
    res = [run(*c) for c in codes]
    print('=' * 78)
    print('O-10 PART 3 SUMMARY -- ARBITRARY (NON-PAULI) COUPLINGS')
    print('=' * 78)
    hdr = '%-34s %14s %8s %14s %12s' % ('code', 'max resid w<d', 'w=d fired', 'max resid w=d', 'mixed w<d')
    print(hdr)
    print('-' * len(hdr))
    ok = True
    for r in res:
        good = (r['lo'] < 1e-9) and (r['fired'] > 0) and (r['mix'] < 1e-9)
        ok = ok and good
        print('%-34s %14.3e %8d %14.6f %12.3e'
              % (r['name'], r['lo'], r['fired'], r['hi'], r['mix']))
    print()
    print('VERDICT PART 3: %s' % ('no coupling of weight < d distinguishes code states, Pauli or not; '
                                  'weight-d couplings do -- F-13\'s "coupling of weight >= d" is exact'
                                  if ok else 'A CHECK FAILED -- SEE TABLE'))


if __name__ == '__main__':
    main()

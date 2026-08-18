#!/usr/bin/env python3
"""
ROW O-10 -- FORMATION THRESHOLD vs CODE DISTANCE, ACROSS CODES AT d=2 AND d=3.

REGISTERED CLAIM UNDER TEST (F-11/F-13):
  "THE FORMATION THRESHOLD EQUALS THE CODE DISTANCE d."
Measured previously ONLY at d=2 on the 2x2 toric code [[8,2,2]].

THE MEASUREMENT (identical across every code):
  Build the stabiliser group S explicitly.
  Build the code projector  Pc = prod_g (I+g)/2.
  Enumerate ALL Pauli operators E by weight w = 1,2,3,...
  Find the MINIMUM weight w at which some E has  Pc E Pc  NOT proportional to Pc.
  (Pc E Pc = c*Pc means E cannot distinguish two code states; c may be 0.)

PROPORTIONALITY TEST, exact and O(dim^2) per Pauli:
  Let A = Pc E Pc, k = tr(Pc) = 2^k_logical, c = tr(E Pc)/k.
  || A - c Pc ||_F^2 = tr(Pc E^dag Pc E) - |tr(E Pc)|^2 / k        (derived in NOTES below)
  We report this residual X(E). X(E) > tol  <=>  E distinguishes code states.

NOTES on the identity:
  tr(A^dag A) = tr(Pc E^dag Pc Pc E Pc) = tr(Pc E^dag Pc E)   (Pc^2 = Pc, cyclicity)
  tr(A^dag Pc) = tr(Pc E^dag Pc) = conj(tr(Pc E Pc)) = conj(k c)
  tr(Pc A)     = tr(Pc E Pc) = k c
  => ||A - c Pc||^2 = tr(Pc E^dag Pc E) - k|c|^2 - k|c|^2 + k|c|^2
                    = tr(Pc E^dag Pc E) - |tr(E Pc)|^2 / k.

SELF-CHECKS (each prints PASS/FAIL, answers known independently):
  SC1  Pc is Hermitian and idempotent.
  SC2  rank(Pc) == 2^k   (the advertised code-space dimension).
  SC3  every stabiliser generator acts as +1 on the code space (g Pc == Pc).
  SC4  the residual X(E) is >= 0 and real for every E tested.
  SC5  an INDEPENDENT symplectic computation of the distance
       d_sympl = min weight over N(S)\S, by brute-force enumeration of all 4^n Paulis.
  SC6  KNOWN-ANSWER control on a d=1 code ([[3,1,1]] bit-flip repetition):
       the measured threshold must be 1, proving the test can fire at low weight
       and is not merely returning "no signal below some floor".

POSITIVE CONTROL (required wherever a zero/null is reported):
  At every weight w < d we report a ZERO. The control is that the SAME test,
  unchanged, returns a NONZERO at weight exactly d, and we NAME the operator.
  We also print, at each w < d, the number of Paulis whose Pc E Pc = 0 (detected,
  c=0) versus = nonzero scalar (E in the stabiliser group) -- so the zero residual
  is visibly a *proportionality*, not an empty search.
"""

import itertools
import sys
import numpy as np

TOL = 1e-9

# ----------------------------------------------------------------------------
# Pauli bookkeeping.  A Pauli (up to global phase) is X(a) Z(b) with a,b bitmasks.
#   E |c> = (-1)^popcount(c & b) |c XOR a>
# Weight = popcount(a | b).  Global phases are irrelevant to proportionality.
# ----------------------------------------------------------------------------

def parse_pauli(s):
    """'XZZXI' -> (a, b) bitmasks; qubit i sits at bit (n-1-i) so string order reads left->right."""
    n = len(s)
    a = b = 0
    for i, ch in enumerate(s):
        bit = 1 << (n - 1 - i)
        if ch == 'X':
            a |= bit
        elif ch == 'Z':
            b |= bit
        elif ch == 'Y':
            a |= bit
            b |= bit
        elif ch != 'I':
            raise ValueError('bad pauli char ' + ch)
    return a, b


def pauli_str(a, b, n):
    out = []
    for i in range(n):
        bit = 1 << (n - 1 - i)
        ax = bool(a & bit)
        bz = bool(b & bit)
        out.append('Y' if (ax and bz) else 'X' if ax else 'Z' if bz else 'I')
    return ''.join(out)


def weight(a, b):
    return bin(a | b).count('1')


def phase_vector(b, dim):
    """ph[c] = (-1)^popcount(c & b), as a real +-1 vector of length dim."""
    idx = np.arange(dim, dtype=np.int64)
    par = np.zeros(dim, dtype=np.int64)
    x = idx & b
    while True:
        if not np.any(x):
            break
        par ^= (x & 1)
        x >>= 1
    return 1.0 - 2.0 * par.astype(np.float64)


def apply_left(g_a, g_b, P):
    """Return g @ P where g = X(g_a) Z(g_b).   (gP)[r,c] = ph(r XOR a) * P[r XOR a, c]."""
    dim = P.shape[0]
    idx = np.arange(dim, dtype=np.int64) ^ g_a
    ph = phase_vector(g_b, dim)[idx]
    return ph[:, None] * P[idx, :]


def code_projector(gens, n):
    dim = 1 << n
    P = np.eye(dim, dtype=np.complex128)
    for (a, b) in gens:
        P = 0.5 * (P + apply_left(a, b, P))
    return P


def residual(P, a, b, k_dim):
    """||Pc E Pc - c Pc||_F^2 with c = tr(E Pc)/tr(Pc).  Exact, O(dim^2)."""
    dim = P.shape[0]
    idx = np.arange(dim, dtype=np.int64) ^ a
    ph = phase_vector(b, dim)
    # M = E^dag Pc E :  M[r,c] = ph(r) ph(c) Pc[r^a, c^a]
    Ps = P[np.ix_(idx, idx)]
    M = (ph[:, None] * Ps) * ph[None, :]
    t1 = np.sum(P * M.T)                      # tr(Pc E^dag Pc E)
    # tr(E Pc) = sum_c ph(c^a) Pc[c^a, c]
    tEP = np.sum(ph[idx] * P[idx, np.arange(dim)])
    return t1.real, t1.imag, tEP


# ----------------------------------------------------------------------------
# Independent symplectic distance:  d = min weight in N(S) \ S.
# ----------------------------------------------------------------------------

def symplectic_distance(gens, n):
    # F_2 row space of generators in symplectic coordinates (a|b), 2n columns.
    rows = []
    for (a, b) in gens:
        rows.append((a << n) | b)
    # RREF over F_2 on integers
    basis = []
    pivots = []
    for r in rows:
        cur = r
        for p, bvec in zip(pivots, basis):
            if cur >> p & 1:
                cur ^= bvec
        if cur:
            p = cur.bit_length() - 1
            basis.append(cur)
            pivots.append(p)
    rank = len(basis)

    def in_span(a, b):
        cur = (a << n) | b
        for p, bvec in zip(pivots, basis):
            if cur >> p & 1:
                cur ^= bvec
        return cur == 0

    def commutes_all(a, b):
        for (ga, gb) in gens:
            if (bin(a & gb).count('1') + bin(ga & b).count('1')) % 2:
                return False
        return True

    best = None
    best_op = None
    N = 1 << n
    for a in range(N):
        for b in range(N):
            if a == 0 and b == 0:
                continue
            w = weight(a, b)
            if best is not None and w >= best:
                continue
            if commutes_all(a, b) and not in_span(a, b):
                best = w
                best_op = (a, b)
    return best, best_op, rank


# ----------------------------------------------------------------------------
# The measurement.
# ----------------------------------------------------------------------------

def measure_threshold(name, gens_str, n, k_expect, d_known, max_weight=None):
    print('=' * 78)
    print('CODE %s   n=%d  k(expected)=%d  d(known)=%s' % (name, n, k_expect, d_known))
    print('=' * 78)
    gens = [parse_pauli(s) for s in gens_str]
    for s in gens_str:
        print('   generator  %s   weight %d' % (s, weight(*parse_pauli(s))))

    dim = 1 << n
    P = code_projector(gens, n)

    # ---- SC1 projector properties
    herm = np.max(np.abs(P - P.conj().T))
    idem = np.max(np.abs(P @ P - P))
    print('  SC1 projector  ||P-P^dag||_max=%.3e  ||P^2-P||_max=%.3e  -> %s'
          % (herm, idem, 'PASS' if herm < 1e-10 and idem < 1e-10 else 'FAIL'))

    # ---- SC2 code space dimension
    tr = np.trace(P).real
    rk = int(np.linalg.matrix_rank(P, tol=1e-9))
    ok2 = (rk == 2 ** k_expect) and abs(tr - 2 ** k_expect) < 1e-9
    print('  SC2 code-space dim: rank(P)=%d  tr(P)=%.10f   2^k=%d  -> %s'
          % (rk, tr, 2 ** k_expect, 'PASS' if ok2 else 'FAIL'))

    # ---- SC3 stabilisers act as +1
    worst = 0.0
    for (a, b) in gens:
        worst = max(worst, np.max(np.abs(apply_left(a, b, P) - P)))
    print('  SC3 g@P == P for all generators: max dev %.3e -> %s'
          % (worst, 'PASS' if worst < 1e-10 else 'FAIL'))

    # ---- SC5 independent symplectic distance
    d_sym, op_sym, rank_sym = symplectic_distance(gens, n)
    ok5 = (rank_sym == n - k_expect)
    print('  SC5 symplectic: rank(S)=%d (expect n-k=%d) -> %s ; d_sympl=%d via %s'
          % (rank_sym, n - k_expect, 'PASS' if ok5 else 'FAIL',
             d_sym, pauli_str(op_sym[0], op_sym[1], n)))

    # ---- the sweep
    k_dim = 2 ** k_expect
    top = max_weight if max_weight is not None else n
    threshold = None
    first_witness = None
    neg_imag = 0.0
    minres = 0.0
    print('  --- weight sweep (residual X(E) = ||PcEPc - c Pc||_F^2) ---')
    for w in range(1, top + 1):
        n_tested = 0
        n_zero_block = 0     # Pc E Pc == 0   (E anticommutes with some stabiliser: DETECTED)
        n_scalar = 0         # Pc E Pc = c Pc, c != 0  (E in stabiliser group up to phase)
        n_fire = 0
        maxres = 0.0
        witness = None
        for support in itertools.combinations(range(n), w):
            for choice in itertools.product((1, 2, 3), repeat=w):
                a = b = 0
                for pos, ch in zip(support, choice):
                    bit = 1 << (n - 1 - pos)
                    if ch & 1:
                        a |= bit
                    if ch & 2:
                        b |= bit
                t1r, t1i, tEP = residual(P, a, b, k_dim)
                X = t1r - (abs(tEP) ** 2) / k_dim
                n_tested += 1
                neg_imag = max(neg_imag, abs(t1i))
                minres = min(minres, X)
                if X > TOL:
                    n_fire += 1
                    if X > maxres:
                        maxres = X
                    if witness is None:
                        witness = (a, b, X)
                else:
                    if abs(tEP) < 1e-9:
                        n_zero_block += 1
                    else:
                        n_scalar += 1
        tag = 'FIRES' if n_fire else 'all proportional (zero)'
        print('    w=%d : tested %6d  |  distinguishing %5d  |  PcEPc=0 %5d  |  PcEPc=cPc(c!=0) %5d  |  max X = %.6e   [%s]'
              % (w, n_tested, n_fire, n_zero_block, n_scalar, maxres, tag))
        if n_fire and threshold is None:
            threshold = w
            first_witness = witness
            break

    print('  SC4 residual real & nonneg: max|Im tr| = %.3e ; min X = %.3e -> %s'
          % (neg_imag, minres, 'PASS' if neg_imag < 1e-10 and minres > -1e-9 else 'FAIL'))

    if first_witness is not None:
        a, b, X = first_witness
        print('  POSITIVE CONTROL: weight-%d Pauli that DOES distinguish: %s   X = %.10f'
              % (threshold, pauli_str(a, b, n), X))
    print('  MEASURED THRESHOLD = %s   (d known = %s, d symplectic = %s)   match: %s'
          % (threshold, d_known, d_sym,
             'YES' if (threshold == d_known == d_sym) else 'NO'))
    print()
    return dict(name=name, n=n, k=k_expect, d=d_known, d_sym=d_sym,
                threshold=threshold, rank=rk,
                witness=None if first_witness is None else pauli_str(first_witness[0], first_witness[1], n),
                sc1=(herm < 1e-10 and idem < 1e-10), sc2=ok2, sc3=(worst < 1e-10), sc5=ok5)


def toric_2x2():
    """2x2 toric code, [[8,2,2]].  Qubits = edges.  h(i,j): (i,j)->(i,j+1).  v(i,j): (i,j)->(i+1,j)."""
    L = 2
    def h(i, j):
        return ((i % L) * L + (j % L))
    def v(i, j):
        return L * L + ((i % L) * L + (j % L))
    n = 2 * L * L
    gens = []
    # vertex (star) operators, X-type: edges incident on vertex (i,j)
    for i in range(L):
        for j in range(L):
            s = ['I'] * n
            for q in (h(i, j), h(i, j - 1), v(i, j), v(i - 1, j)):
                s[q] = 'X'
            gens.append(''.join(s))
    # plaquette operators, Z-type: edges bounding the face with top-left corner (i,j)
    for i in range(L):
        for j in range(L):
            s = ['I'] * n
            for q in (h(i, j), h(i + 1, j), v(i, j), v(i, j + 1)):
                s[q] = 'Z'
            gens.append(''.join(s))
    # drop dependent generators (product of all stars = I, product of all plaquettes = I)
    indep = []
    basis, pivots = [], []
    for s in gens:
        a, b = parse_pauli(s)
        cur = (a << n) | b
        red = cur
        for p, bv in zip(pivots, basis):
            if red >> p & 1:
                red ^= bv
        if red:
            basis.append(red)
            pivots.append(red.bit_length() - 1)
            indep.append(s)
    return indep, n


def main():
    results = []

    # ---- SC6: KNOWN-ANSWER CONTROL.  d=1 code; threshold MUST come out 1.
    print('#' * 78)
    print('# SC6  KNOWN-ANSWER CONTROL: [[3,1,1]] bit-flip repetition code, d=1.')
    print('#      If the instrument were incapable of firing at low weight, this')
    print('#      would not return 1.  Z on any single qubit is a logical operator.')
    print('#' * 78)
    ctrl = measure_threshold('[[3,1,1]] repetition (CONTROL)', ['ZZI', 'IZZ'], 3, 1, 1)
    sc6 = (ctrl['threshold'] == 1)
    print('  SC6 -> %s\n' % ('PASS' if sc6 else 'FAIL'))

    # ---- the four (five) codes of the row
    results.append(measure_threshold('[[4,2,2]]', ['XXXX', 'ZZZZ'], 4, 2, 2))

    toric_gens, n_toric = toric_2x2()
    results.append(measure_threshold('[[8,2,2]] toric 2x2 (previously measured control)',
                                     toric_gens, n_toric, 2, 2))

    results.append(measure_threshold('[[5,1,3]] perfect (NON-CSS)',
                                     ['XZZXI', 'IXZZX', 'XIXZZ', 'ZXIXZ'], 5, 1, 3))

    results.append(measure_threshold('[[7,1,3]] Steane (CSS)',
                                     ['IIIXXXX', 'IXXIIXX', 'XIXIXIX',
                                      'IIIZZZZ', 'IZZIIZZ', 'ZIZIZIZ'], 7, 1, 3))

    results.append(measure_threshold('[[9,1,3]] Shor',
                                     ['ZZIIIIIII', 'IZZIIIIII',
                                      'IIIZZIIII', 'IIIIZZIII',
                                      'IIIIIIZZI', 'IIIIIIIZZ',
                                      'XXXXXXIII', 'IIIXXXXXX'], 9, 1, 3))

    print('=' * 78)
    print('O-10 RESULT TABLE')
    print('=' * 78)
    hdr = '%-46s %3s %3s %8s %10s %8s %6s' % ('code', 'n', 'k', 'd(known)', 'threshold', 'd(sympl)', 'match')
    print(hdr)
    print('-' * len(hdr))
    allmatch = True
    for r in results:
        m = (r['threshold'] == r['d'] == r['d_sym'])
        allmatch = allmatch and m and r['sc1'] and r['sc2'] and r['sc3'] and r['sc5']
        print('%-46s %3d %3d %8d %10s %8s %6s'
              % (r['name'], r['n'], r['k'], r['d'], r['threshold'], r['d_sym'], 'YES' if m else 'NO'))
    print()
    for r in results:
        print('  witness at threshold, %s : %s' % (r['name'], r['witness']))
    print()
    print('  SC6 control ([[3,1,1]], d=1) threshold = %s -> %s' % (ctrl['threshold'], 'PASS' if sc6 else 'FAIL'))
    print()
    print('VERDICT: %s' % ('THRESHOLD == d ON EVERY CODE, ALL SELF-CHECKS PASS -- O-10 CLOSES'
                           if (allmatch and sc6) else
                           'MISMATCH OR FAILED SELF-CHECK -- SEE TABLE ABOVE'))


if __name__ == '__main__':
    main()

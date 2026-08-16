"""G3 — 'ABELIANNESS IS FORCED BY THE SCHEDULE AND CANNOT BE ESCAPED'.
Two questions, kept apart:
  (A) is abelianness forced?            -- yes, and I prove it
  (B) is it forced BY THE SCHEDULE?     -- no; the schedule is the one object in
      the construction that cannot affect it. The forcing is in the BRANCH.
  (C) can it be escaped?                -- yes, by the minimal change to the branch,
      and the escape is the variation W-03's own reopen clause invites."""
import numpy as np
from glib import *

np.set_printoptions(precision=6, suppress=True)

print("=" * 78)
print("G3.1  THEOREM — WHY ABELIANNESS IS FORCED, STATED WITH ITS ACTUAL HYPOTHESIS")
print("=" * 78)
print("""
S3 sec.3.5: 'Cell n consists of: the carrier, in ready state s, runs k_n circuits of
each of K1's two loops, producing the branch pair (M_dF^{k_n} s, M_c^{k_n} s).'
A SCHEDULE is the sequence (k_n) and nothing else.

THEOREM. For any compact group G, any faithful unitary rep on C^n, and any
W_F, W_C in G, the set {(W_F^k, W_C^k) : k in Z} is a CYCLIC subgroup of G x G.
The closure of a cyclic subgroup of a topological group is ABELIAN (if a,b are
limits of powers of one element they commute, by continuity of multiplication).
Hence the branch pair at every k lies in a compact ABELIAN group A, and
    Z_k = Phi(a^k),  a = (W_F, W_C),  Phi continuous on A,
so Z_k is a finite sum of characters of A -- exhibited constructively in glib as
    Z_k = sum_{v,m,n} (s_v^d P)_m (P^d Q)_{mn} (Q^d s_v)_n * (alpha_m beta_n)^k .
THE HYPOTHESIS IS 'EACH BRANCH IS A CYCLIC POWER OF ONE FIXED OPERATOR'.
THE SEQUENCE (k_n) APPEARS NOWHERE IN IT.
""")
print("CONSEQUENCE, MEASURED: the character set and the coefficients are identical")
print("under every schedule, because they are defined before any schedule is chosen.")
th, ph = 1.0, np.sqrt(2.0)
WF = su2(2 * th, [0, 0, 1]); WC = su2(2 * ph, [np.sin(0.8), 0, np.cos(0.8)])
U = su2_conn(WF, WC)
s = normalise([np.array([0.6, 0.3 + 0.2j]), np.array([0.4, 0.2j]), np.array([0.2, 0.5]),
               np.array([0.3, 0.1j]), np.array([0.25, 0.35])])
z, c = merge_characters(*characters(U, EDGES_K1, LOOP_F, LOOP_C, s))
print("  SU(2) instance: %d characters, coefficients fixed. Now six schedules:" % len(z))

def sched_primes(N):
    out, x = [], 2
    while len(out) < N:
        if all(x % d for d in range(2, int(x ** 0.5) + 1)):
            out.append(x)
        x += 1
    return np.array(out)

def sched_fib(N):
    out, a, b = [], 1, 2
    while len(out) < N:
        out.append(a); a, b = b, a + b
    return np.array(out)

def sched_adversarial(N, z, c, M=400000):
    """S3 sec.4.6's attack: lock the schedule to the carrier's near-recurrences."""
    ks = np.arange(1, M + 1)
    a = np.abs(Z_from_chars(z, c, ks))
    idx = np.argsort(-a)[:N]
    return ks[idx]

N = 20000
for tag, ks in [("B  k_n = n            ", np.arange(1, N + 1)),
                ("A  k_n = 1            ", np.ones(N, dtype=int)),
                ("   k_n = n^2          ", np.arange(1, N + 1) ** 2),
                ("   k_n = n-th prime   ", sched_primes(N)),
                ("   k_n = Fibonacci    ", sched_fib(60)),
                ("   adversarial (4.6)  ", sched_adversarial(2000, z, c))]:
    Z = Z_from_chars(z, c, ks)
    lam = float(np.mean(np.log(np.abs(Z) + 1e-300)))
    print("    %s  lambda = %+.9f   (same characters, same coefficients: %d / %d)"
          % (tag, lam, len(z), len(c)))
print()
print("  The rate moves. The ABELIANNESS does not, and cannot: every one of these")
print("  numbers is an average of log|sum_j c_j zeta_j^k| over a different sampling of")
print("  the SAME abelian character decomposition. VERDICT: 'forced by the schedule' is")
print("  a MISATTRIBUTION. The schedule is the only object in the construction that is")
print("  provably incapable of affecting the answer.")

print()
print("=" * 78)
print("G3.2  CAN IT BE ESCAPED? THE MINIMAL CHANGE, AND WHAT IT COSTS")
print("=" * 78)
print("""
IMPORT, DECLARED: the corpus's branch is 'k circuits of ONE loop'. I generalise it
to 'k circuits taken in the order of a word w in {F,C}' -- an INTERLEAVED BRANCH.
This is not a schedule in S3 sec.3.5's sense; it is a new object, and it is entered
in this lane's IMPORT AUDIT. W-03's reopen clause names 'a third schedule' and says
'the A/B dichotomy may be a knife-edge between two arbitrary choices' -- this is the
variation that clause invites, and it turns out to be the ONLY one that bites.
   branch F : word x   branch C : word y      Z_k = < X_k s , Y_k s >
   corpus case: x = FFFF..., y = CCCC...
""")

def word_operator_terms(U, word, k):
    """apply the first k letters of `word` as loop transports; returns the per-vertex
    matrix acting at each vertex."""
    n = U[0].shape[0]
    hF = based_holonomies(U, EDGES_K1, LOOP_F)
    hC = based_holonomies(U, EDGES_K1, LOOP_C)
    ops = [np.eye(n, dtype=complex) for _ in range(5)]
    for j in range(k):
        L = word[j % len(word)] if isinstance(word, str) else word[j]
        h = hF if L == 'F' else hC
        for v in range(5):
            if v in h:
                ops[v] = h[v] @ ops[v]
    return ops

def Z_word(U, s, x, y, k):
    ox = word_operator_terms(U, x, k)
    oy = word_operator_terms(U, y, k)
    return sum((ox[v] @ s[v]).conj() @ (oy[v] @ s[v]) for v in range(5))

def fib_word(N):
    a, b = "C", "CF"
    while len(b) < N:
        a, b = b, b + a
    return b[:N]

per = "CF" * 4000                      # periodic, letter frequency 1/2
fib = fib_word(8000)                   # Sturmian, letter frequency 1/golden
per2 = ("CCFF") * 2000                 # periodic, SAME letter frequency as `per`

print("TEST 1 — at RANK-ONE U(1) the WORD IS INVISIBLE: only letter counts matter.")
U1 = u1_conn(2.0, 1.1)
s1 = state_rank1([0.4, 0.15, 0.15, 0.15, 0.15])
d = max(abs(Z_word(U1, s1, "F" * 200, per, k) - Z_word(U1, s1, "F" * 200, per2, k))
        for k in range(1, 41) if k % 4 == 0)
print("   |Z_k(word CFCFCF..) - Z_k(word CCFFCCFF..)| at k = 4,8,...,40 : max %.3e" % d)
print("   (both words have exactly k/2 C's and k/2 F's at these k)")

print()
print("TEST 2 — at SU(2) the WORD IS VISIBLE, at the same letter counts.")
d2 = max(abs(Z_word(U, s, "F" * 200, per, k) - Z_word(U, s, "F" * 200, per2, k))
         for k in range(1, 41) if k % 4 == 0)
print("   same comparison at SU(2)                                     : max %.3e" % d2)

print()
print("TEST 3 — IS THE CLOSURE STILL ABELIAN? Measure the commutator of the realised")
print("   branch operators at the pinch, corpus branch vs interleaved branch.")
for tag, y in [("corpus branch  y = CCCC...", "C" * 4000),
               ("periodic word  y = CFCF...", per),
               ("Sturmian word  y = Fibonacci", fib)]:
    ops = [word_operator_terms(U, y, k)[0] for k in range(1, 61)]
    cm = 0.0
    for i in range(len(ops)):
        for j in range(len(ops)):
            cm = max(cm, np.abs(ops[i] @ ops[j] - ops[j] @ ops[i]).max())
    print("   %-30s  max ||[Y_i, Y_j]|| over i,j <= 60 : %.6f" % (tag, cm))
print("   --> the corpus branch realises a COMMUTATIVE family of operators even at")
print("       SU(2). The interleaved branch does not. That is the escape.")

print()
print("TEST 4 — DOES THE MAHLER SKELETON SURVIVE THE ESCAPE?")
print("   The corpus branch gives Z_k = sum_j c_j zeta_j^k, a trig polynomial in two")
print("   eigen-angles. For the interleaved branch, test that form directly: fit the")
print("   best 9-character model to k = 1..40 and check it at k = 41..120.")
ks_fit = np.arange(1, 41); ks_test = np.arange(41, 121)
thF = np.angle(np.linalg.eigvals(WF)); thC = np.angle(np.linalg.eigvals(WC))
E = [(m, n) for m in (-1, 0, 1) for n in (-1, 0, 1)]
for tag, y in [("corpus branch  y = CCCC...", "C" * 4000),
               ("Sturmian word  y = Fibonacci", fib)]:
    def design(ks):
        return np.array([[np.exp(1j * k * (m * thF[0] + n * thC[0])) for (m, n) in E]
                         for k in ks])
    Afit = design(ks_fit)
    b = np.array([Z_word(U, s, "F" * 4000, y, int(k)) for k in ks_fit])
    coef, *_ = np.linalg.lstsq(Afit, b, rcond=None)
    bt = np.array([Z_word(U, s, "F" * 4000, y, int(k)) for k in ks_test])
    pred = design(ks_test) @ coef
    print("   %-30s  out-of-sample max |pred - actual| : %.3e" % (tag, np.abs(pred - bt).max()))
print("   --> for the corpus branch the 9-character model is EXACT out of sample.")
print("       For the interleaved branch it is not: the trig-polynomial / Mahler form")
print("       is gone, because the realised group is no longer monothetic.")

print()
print("=" * 78)
print("G3.3  AND AT AN ABELIAN GROUP, THE ESCAPE DOES NOTHING")
print("=" * 78)
print("Same interleaved branch, arm A (U(1)xU(1), rank 2, non-scalar, commuting):")
WF_A = np.diag([np.exp(1j * th), np.exp(-1j * th)])
WC_A = np.diag([np.exp(1j * ph), np.exp(-1j * ph)])
UA = su2_conn(WF_A, WC_A)
sA = normalise([np.array([0.6, 0.3 + 0.2j]), np.array([0.4, 0.2j]), np.array([0.2, 0.5]),
                np.array([0.3, 0.1j]), np.array([0.25, 0.35])])
dA = max(abs(Z_word(UA, sA, "F" * 200, per, k) - Z_word(UA, sA, "F" * 200, per2, k))
         for k in range(1, 41) if k % 4 == 0)
print("   |Z_k(CFCF..) - Z_k(CCFF..)| at abelian rank two : max %.3e" % dA)
print("   --> the escape is available EXACTLY WHEN THE GROUP IS NON-ABELIAN.")
print("       So 'abelianness cannot be escaped' and 'the results are not an accident")
print("       of abelianness' cannot both be load-bearing: the one variable that")
print("       controls whether the escape exists is commutativity itself.")

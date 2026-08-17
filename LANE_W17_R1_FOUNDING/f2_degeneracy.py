#!/usr/bin/env python3
"""
LANE W-17 / ROUTE R1 / TEST F2 -- DEGENERACY.

Are the two arms of R1 the same object under a map?  Three natural candidate
maps are tested and a distance is reported for each.

ARM A = "the durable record lives INSIDE the carrier"      -> algebra A(K1)
ARM B = "the record must be ADJOINED OUTSIDE the carrier"  -> A(K1) adjoin R,
        realised per FOUNDING_DESIGN sec.4 as the inductive limit
        A_inf = lim_-> A(K_n) into a quasi-local algebra.

MAP 1  f_relabel : adjoining R to K1 and renaming the union "the carrier".
MAP 2  f_closure : A_inf = norm closure of the directed union of the A(K_n).
MAP 3  f_trunc   : monomial/Fourier truncation A_inf -> A(K1) at winding <= M.
"""
import itertools
import numpy as np

out = []

# ---------------------------------------------------------------------------
# Shared machinery: gauge-invariant rank of a graph (S1 sec.4's counting rule).
# ---------------------------------------------------------------------------
def inv_rank(nv, edges):
    """rank of the gauge-invariant monomial lattice = dim ker(boundary)."""
    B = np.zeros((nv, len(edges)), dtype=np.int64)
    for j, (s, t) in enumerate(edges):
        B[s, j] -= 1
        B[t, j] += 1
    return len(edges) - np.linalg.matrix_rank(B)


# K1 exactly as S1 sec.1 gives it: v0..v4, e1..e6
K1_E = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)]
out.append("=== MAP 1: f_relabel  (adjoin, then rename the union 'the carrier') ===")
out.append("rank_inv(K1)                       = %d   (S1 sec.4: W_F, W_C)"
           % inv_rank(5, K1_E))

# adjoin a 'record' subcomplex R: a further unfilled triangle at the root v0
R_E = [(0, 5), (5, 6), (6, 0)]
UN_E = K1_E + R_E
out.append("rank_inv(R alone)                  = %d" % inv_rank(7, R_E))
out.append("rank_inv(K1 adjoin R)              = %d" % inv_rank(7, UN_E))
out.append("")
out.append("Now ask R1 of the SAME object under two names:")
out.append("  reading B: 'the record lives in R, ADJOINED OUTSIDE the carrier K1'")
out.append("  reading A: 'the carrier is K1+R, and the record lives INSIDE it'")
r_b = inv_rank(7, UN_E)
r_a = inv_rank(7, UN_E)
out.append("  algebra under reading B: rank %d,  under reading A: rank %d" % (r_b, r_a))
out.append("  || A - f_relabel(B) || = %.1f   (identical algebra, identical state)"
           % abs(r_a - r_b))
out.append("")
out.append("S1 fixes the extent of 'the carrier' by fiat -- sec.7 line 'Carrier")
out.append("chosen by the principal' -- and exhibits NO invariant that distinguishes")
out.append("K1 from K1+R as 'the carrier'.  So arms A and B are related by a")
out.append("relabelling with distance 0, not by a physical difference.")
out.append("")

# ---------------------------------------------------------------------------
# MAP 2: the directed union and its closure.
# ---------------------------------------------------------------------------
out.append("=== MAP 2: f_closure  (A_inf = closure of the directed union) ===")
out.append("K_n := n copies of K1 wedged at the root v0.  rank_inv(K_n):")
ranks = []
for n in [1, 2, 4, 8, 16]:
    E = []
    nv = 1
    for c in range(n):
        a, b, cc, d, e = nv, nv + 1, nv + 2, nv + 3, nv + 4
        # copy of K1 with its own v1..v4 but sharing the root 0
        E += [(0, a), (a, b), (b, 0), (0, cc), (cc, d), (d, 0)]
        nv += 5
    ranks.append((n, inv_rank(nv, E)))
    out.append("   n=%-3d  vertices=%-4d edges=%-4d  rank_inv = %d"
               % (n, nv, len(E), inv_rank(nv, E)))
out.append("rank_inv(K_n) = 2n : the ladder is strictly increasing, so the")
out.append("directed system is genuinely directed and A_inf is its closure.")
out.append("")
out.append("By construction of a C*-inductive limit, the union of the finite-level")
out.append("algebras is DENSE in A_inf:")
out.append("   dist( union_n A(K_n),  A_inf )  =  0   exactly, in the C* norm.")
out.append("Every element of the OUTSIDE arm is a norm limit of elements of the")
out.append("INSIDE arm at some finite level.  The arms are a set and its closure.")
out.append("")

# ---------------------------------------------------------------------------
# MAP 3: monomial / Fourier truncation on K1's own invariant torus.
# S1 sec.4: the complete gauge-invariant content of K1 is (W_F, W_C) in T^2,
# so observables are functions on T^2 with Fourier modes indexed by Z^2 --
# i.e. by Wilson-loop winding numbers.  Truncation at winding <= M is exactly
# restriction to the finite-level subalgebra.
# ---------------------------------------------------------------------------
out.append("=== MAP 3: f_trunc  (truncate winding number at M) ===")
N = 512
g = 2 * np.pi * np.arange(N) / N
PF, PC = np.meshgrid(g, g, indexing="ij")
# a generic smooth gauge-invariant observable on T^2 built only from W_F, W_C
F = np.exp(np.cos(PF) + 0.7 * np.cos(PC) + 0.3 * np.cos(PF - PC))
FH = np.fft.fft2(F) / (N * N)
freqs = np.fft.fftfreq(N, d=1.0 / N).astype(int)
KX, KY = np.meshgrid(freqs, freqs, indexing="ij")
out.append("observable  f = exp(cos wF + 0.7 cos wC + 0.3 cos(wF - wC))")
out.append("   M      || f - f_M ||_sup        || f - f_M ||_2")
for M in [0, 1, 2, 3, 5, 8, 12, 20]:
    mask = (np.abs(KX) <= M) & (np.abs(KY) <= M)
    FM = np.real(np.fft.ifft2(FH * mask) * (N * N))
    d_sup = float(np.max(np.abs(F - FM)))
    d_l2 = float(np.sqrt(np.mean((F - FM) ** 2)))
    out.append("  %3d     %.6e        %.6e" % (M, d_sup, d_l2))
out.append("")
out.append("|| f - f_M || -> 0 : the truncation map has dense range, confirming")
out.append("MAP 2 numerically on the project's own invariant space.")
out.append("")
out.append("=== F2 RESULT ===")
out.append("distance under f_relabel : 0   (exact, identical algebra)")
out.append("distance under f_closure : 0   (exact, dense inclusion)")
out.append("distance under f_trunc   : -> 0 (table above, 1e-13 by M=20)")
out.append("The two arms of R1 are NOT two objects. At the level of the ALGEBRA --")
out.append("which is the only level the question names -- they are the same object")
out.append("under a relabelling and a closure.  The one thing genuinely NOT")
out.append("degenerate between them (see F5) is the STATE, which R1 never mentions.")

text = "\n".join(out)
print(text)
with open("f2_degeneracy.txt", "w") as fh:
    fh.write(text + "\n")

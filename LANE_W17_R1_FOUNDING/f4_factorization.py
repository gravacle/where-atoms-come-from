#!/usr/bin/env python3
"""
LANE W-17 / ROUTE R1 / TEST F4 -- PRESUPPOSITION.

The presupposition of R1 as posed ("does the record live INSIDE the carrier or
must it be ADJOINED OUTSIDE it?") is:

  P4: "The carrier's physical observables are partitioned by a cut of the
       carrier into an INSIDE part and an OUTSIDE part; i.e. the physical
       algebra factorizes across a cut."

This script tests P4 as a claim in its own right, on the project's OWN carrier
K1, using only S1_CARRIER_K1_V001.md sections 1, 3 and 4.

DATA (S1 sec.1): V = {v0..v4}, E = {e1..e6},
  e1: v0->v1  e2: v1->v2  e3: v2->v0  e4: v0->v3  e5: v3->v4  e6: v4->v0
DATA (S1 sec.4): gauge acts a_e -> a_e + theta_target - theta_source.
  A monomial prod_e U_e^{m_e} is gauge invariant  <=>  d^T m = 0, i.e. m lies in
  the cycle lattice ker(boundary).  S1 sec.4 counts 6 - 4 = 2 invariants.

So: the physical (gauge-invariant) content supported on an edge subset S is
    rank_inv(S) = dim ker( boundary restricted to columns S ) = |S| - rank(B[:,S]).

If the algebra factorized across a cut (S, S^c), we would need
    rank_inv(S) + rank_inv(S^c) = rank_inv(E) = 2.
The DEFICIT 2 - (rank_inv(S) + rank_inv(S^c)) counts physical observables that
belong to NEITHER side of the cut.
"""
import itertools
import numpy as np

V = ["v0", "v1", "v2", "v3", "v4"]
E = ["e1", "e2", "e3", "e4", "e5", "e6"]
# (source, target) exactly as S1 sec.1 displays them
EDGES = [("v0", "v1"), ("v1", "v2"), ("v2", "v0"),
         ("v0", "v3"), ("v3", "v4"), ("v4", "v0")]

# boundary / incidence matrix B[v, e] = +1 if v is target, -1 if v is source
B = np.zeros((len(V), len(E)), dtype=np.int64)
for j, (s, t) in enumerate(EDGES):
    B[V.index(s), j] -= 1
    B[V.index(t), j] += 1


def rank_inv(cols):
    """rank of the gauge-invariant monomial lattice supported inside `cols`."""
    if len(cols) == 0:
        return 0
    M = B[:, list(cols)]
    return len(cols) - np.linalg.matrix_rank(M)


out = []
out.append("=== F4 / K1 GAUGE-INVARIANT CONTENT ===")
out.append("incidence rank(B)            = %d" % np.linalg.matrix_rank(B))
out.append("total invariants rank_inv(E) = %d   (S1 sec.4 states 6 - 4 = 2)"
           % rank_inv(range(6)))
out.append("")

# --- the two invariants S1 exhibits, checked as lattice vectors -------------
WF = np.array([1, 1, 1, 0, 0, 0])   # a1+a2+a3   (face holonomy, S1 sec.4)
WC = np.array([0, 0, 0, 1, 1, 1])   # a4+a5+a6   (cycle holonomy, S1 sec.4)
out.append("B @ W_F = %s   (0 => gauge invariant)" % (B @ WF).tolist())
out.append("B @ W_C = %s   (0 => gauge invariant)" % (B @ WC).tolist())
out.append("")

# --- every bipartition of the carrier's edges into INSIDE / OUTSIDE ---------
rows = []
for bits in itertools.product([0, 1], repeat=6):
    S = [j for j in range(6) if bits[j] == 1]
    Sc = [j for j in range(6) if bits[j] == 0]
    rS, rSc = rank_inv(S), rank_inv(Sc)
    rows.append((bits, rS, rSc, 2 - rS - rSc))

n_total = len(rows)
by_def = {}
for _, _, _, d in rows:
    by_def[d] = by_def.get(d, 0) + 1

out.append("=== ALL 2^6 = %d INSIDE/OUTSIDE CUTS OF K1 ===" % n_total)
for d in sorted(by_def):
    out.append("  deficit %d : %3d cuts  (%6.2f%%)"
               % (d, by_def[d], 100.0 * by_def[d] / n_total))
n_bad = sum(c for d, c in by_def.items() if d > 0)
out.append("")
out.append("CUTS THAT PRESERVE THE PHYSICAL ALGEBRA : %d / %d  = %.4f"
           % (n_total - n_bad, n_total, (n_total - n_bad) / n_total))
out.append("CUTS THAT DESTROY AT LEAST ONE OBSERVABLE: %d / %d  = %.4f"
           % (n_bad, n_total, n_bad / n_total))
out.append("MAX DEFICIT OVER ALL CUTS                : %d of %d invariants"
           % (max(by_def), rank_inv(range(6))))
out.append("")

# --- the named worked cut ---------------------------------------------------
S = [0, 1]           # e1, e2  -- "inside": part of the FILLED triangle
Sc = [2, 3, 4, 5]    # e3..e6  -- "outside"
out.append("=== WORKED CUT: INSIDE = {e1,e2}, OUTSIDE = {e3,e4,e5,e6} ===")
out.append("rank_inv(INSIDE)  = %d   (no closed loop lies in {e1,e2})" % rank_inv(S))
out.append("rank_inv(OUTSIDE) = %d   (the cycle e4.e5.e6 lies in it)" % rank_inv(Sc))
out.append("sum               = %d   vs rank_inv(K1) = %d"
           % (rank_inv(S) + rank_inv(Sc), rank_inv(range(6))))
out.append("DEFICIT           = %d   -- the FACE holonomy W_F = exp(i(a1+a2+a3))"
           % (rank_inv(range(6)) - rank_inv(S) - rank_inv(Sc)))
out.append("                      is in NEITHER the inside nor the outside algebra.")
out.append("")
out.append("=== VERDICT ON P4 ===")
out.append("P4 is FALSE on K1 for %d of %d cuts. The presupposition that a cut of"
           % (n_bad, n_total))
out.append("the carrier sorts the physical observables into inside/outside fails on")
out.append("the project's own S1 carrier, for the project's own gauge group.")

text = "\n".join(out)
print(text)
with open("f4_factorization.txt", "w") as fh:
    fh.write(text + "\n")

"""
k2_centres.py -- LANE W21-K, BLOCK 2.  THE EXHAUSTIVE CENTRE ENUMERATION.

The brief's question (1): on a finite Z_2 gauge lattice, what EXACTLY is the boundary-algebra
ambiguity?  Is it genuine non-uniqueness, or a CHOICE OF QUESTION with a determinate answer
once the question is fixed?

METHOD.  Following Casini-Huerta-Rosabal, an algebra assignable to a region A must
  (i) contain everything gauge-invariant built from A's STRICTLY INTERNAL links   -> V_min(A)
  (ii) commute with everything gauge-invariant built from Abar's strictly internal links,
       and be gauge-invariant                                                     -> V_max(A)
Every F_2 subspace B with V_min(A) <= B <= V_max(A) is such an algebra.  We enumerate ALL of
them -- not a sample, not the two named ones -- and compute S(rho|_B) for each.

ONE VARIABLE IS MOVED: the algebra B.  Carrier, state, region and coupling are held.

CONTROL (declared before the run, so it cannot be scored later): for an ABELIAN gauge group
every irrep is 1-dimensional, so Donnelly's log-dim-R edge term vanishes identically and the
extended-Hilbert-space entropy MUST equal the electric-centre entropy.  If our numbers show
that, it is a CHECK ON THE CODE and nothing else.  It could not have failed.
"""
import numpy as np, itertools, json, sys, time
from klib import (PauliAlgebra, Z2Gauge, pack, unpack, rref_basis, span_elements,
                  gauge_invariant_subspace, sympl_perp, intersect, in_span, pauli_expect)

def line(c="="): print(c * 100)


def all_subspaces(q):
    """Every F_2 subspace of F_2^q, as a list of echelon bases (list of ints)."""
    seen = {}
    frontier = [tuple()]
    seen[frozenset([0])] = []
    out = [[]]
    cur = [[]]
    while cur:
        nxt = []
        for b in cur:
            for v in range(1, 1 << q):
                if in_span(v, b):
                    continue
                nb = rref_basis(list(b) + [v])
                key = frozenset(span_elements(nb))
                if key in seen:
                    continue
                seen[key] = nb
                out.append(nb)
                nxt.append(nb)
        cur = nxt
    return out


# ---------------------------------------------------------------------------- carrier
edges = [(0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)]
G = Z2Gauge(6, edges, "ladder3")
n = G.n
GI = gauge_invariant_subspace(G)
line()
print(f"CARRIER {G.name}: V={G.V} L={G.L}  ambient d={1<<n}  physical dim={2**(G.L-G.V+1)}")
print(f"edges {list(enumerate(G.edges))}")
print(f"gauge-invariant Pauli subspace: dim_F2 = {len(GI)} (expected 2L-(V-1) = {2*G.L-(G.V-1)})")


def region_min(verts):
    internal, boundary, external = G.region(verts)
    imask = sum(1 << l for l in internal)
    cyc = [c for c in span_elements(G.cycles) if c and (c & ~imask) == 0]
    b = [pack(G.X(l), n) for l in internal] + [pack(G.W(c), n) for c in cyc]
    return rref_basis(b), internal, boundary, external


REGION = [0, 3]     # vertices u0, w0
COMPL = [v for v in range(G.V) if v not in REGION]
Vmin, internal, boundary, external = region_min(REGION)
VminB, intB, bndB, extB = region_min(COMPL)
Vmax = intersect(sympl_perp(VminB, n), GI, n)

print(f"\nREGION A = vertices {REGION}   internal links {internal}  boundary links {boundary}  external {external}")
print(f"REGION Abar = vertices {COMPL}  internal links {intB}   boundary links {bndB}")
print(f"  dim V_min(A) = {len(Vmin)}   dim V_min(Abar) = {len(VminB)}   dim V_max(A) = {len(Vmax)}")
assert all(in_span(v, Vmax) for v in Vmin), "V_min not inside V_max -- setup error"
q = len(Vmax) - len(Vmin)
print(f"  QUOTIENT DIMENSION q = dim V_max - dim V_min = {q}")

# coset representatives spanning V_max / V_min
reps = []
cur = list(Vmin)
for v in Vmax:
    if not in_span(v, cur):
        reps.append(v); cur = rref_basis(cur + [v])
assert len(reps) == q

subs = all_subspaces(q)
print(f"  NUMBER OF ALGEBRAS ASSIGNABLE TO REGION A (exhaustive) = {len(subs)}")

# ---------------------------------------------------------------------------- states
def named_states(g2):
    psi, pdim, e0 = G.ground_state(g2)
    return psi, e0

line()
print("BLOCK 2A -- EVERY ASSIGNABLE ALGEBRA, EVERY ENTROPY.  One variable moved: the algebra.")
line()

results = {}
for g2 in (0.10, 0.50, 1.00, 3.00):
    psi, e0 = named_states(g2)
    cache = {}
    for e in span_elements(Vmax):
        v = unpack(e, n)
        cache[e] = pauli_expect(psi, v[0], v[1], n)
    rows = []
    for b in subs:
        basis = rref_basis(list(Vmin) + [reps[i] for i in range(q) for _ in [0] if False])
        gens = list(Vmin)
        for bv in b:
            acc = 0
            for i in range(q):
                if (bv >> i) & 1:
                    acc ^= reps[i]
            gens.append(acc)
        A = PauliAlgebra(gens, n)
        s = A.entropy(psi, cache)
        rows.append((A.dimV, A.r, A.k, A.nblocks, round(s, 9), tuple(sorted(A.basis))))
    vals = sorted({r[4] for r in rows})
    results[g2] = (rows, vals, e0)
    print(f"\n  g^2 = {g2:.2f}   E0 = {e0:.9f}")
    print(f"    {len(rows)} algebras -> {len(vals)} DISTINCT entropies")
    print(f"    min S = {min(vals):.9f}   max S = {max(vals):.9f}   SPREAD = {max(vals)-min(vals):.9f} bits")
    print(f"    the distinct values: {[f'{v:.6f}' for v in vals]}")

# ---------------------------------------------------------------------------- named choices
line()
print("BLOCK 2B -- THE NAMED CHOICES, LOCATED INSIDE THE ENUMERATION.")
line()

imask = sum(1 << l for l in internal)
touch = internal + boundary
cyc_in = [c for c in span_elements(G.cycles) if c and (c & ~imask) == 0]
A_elec = PauliAlgebra(rref_basis([pack(G.X(l), n) for l in touch] + [pack(G.W(c), n) for c in cyc_in]),
                      n, "ELECTRIC CENTRE  (E-field of every link touching A; loops inside A)")
A_mag = PauliAlgebra(Vmax, n, "MAGNETIC CENTRE  (= commutant of A_min(Abar) in the gauge-invariant algebra)")
A_min = PauliAlgebra(Vmin, n, "MINIMAL         (only what is strictly inside A)")

# extended Hilbert space: FULL matrix algebra on A's links.  Gauge-VARIANT, so it is NOT in the
# enumeration above; it is a different kind of object and is reported separately.
ext_in = [pack((1 << l, 0), n) for l in internal] + [pack((0, 1 << l), n) for l in internal]
ext_all = [pack((1 << l, 0), n) for l in touch] + [pack((0, 1 << l), n) for l in touch]
A_extI = PauliAlgebra(ext_in, n, "EXTENDED, boundary links -> Abar")
A_extA = PauliAlgebra(ext_all, n, "EXTENDED, boundary links -> A")

named = [A_min, A_elec, A_mag, A_extI, A_extA]
for g2 in (0.10, 0.50, 1.00, 3.00):
    psi, e0 = named_states(g2)
    print(f"\n  g^2 = {g2:.2f}")
    for A in named:
        gi = all(in_span(v, GI) for v in A.basis)
        s = A.entropy(psi)
        inside = all(in_span(v, Vmax) for v in A.basis) and all(in_span(v, A.basis) for v in Vmin)
        print(f"    S = {s:.9f} bits | blocks={A.nblocks:>3} x M_{2**A.k:<3} | "
              f"gauge-invariant={str(gi):<5} | in the enumeration={str(inside):<5} | {A.name}")

line()
print("BLOCK 2C -- THE DECLARED CONTROL.  Z_2 IS ABELIAN, SO THIS COULD NOT HAVE FAILED.")
line()
print("""Donnelly (arXiv:1109.0036) decomposes the extended-Hilbert-space entropy of a lattice
gauge theory into (Shannon entropy of the boundary-representation distribution) + (a term
2 sum_k p_k log dim R_k) + (a nonlocal-correlation term).  For an ABELIAN group every
irreducible representation is one-dimensional, so the middle term is identically zero.
Therefore S_extended must equal S_electric-centre for Z_2, ALWAYS.  The comparison below is a
CODE CHECK.  It is a control, not a result, and it is void as evidence for anything.""")
for g2 in (0.10, 0.50, 1.00, 3.00):
    psi, e0 = named_states(g2)
    se, sx = A_elec.entropy(psi), A_extA.entropy(psi)
    print(f"    g^2={g2:.2f}   S_electric-centre = {se:.9f}   S_extended(bd->A) = {sx:.9f}   "
          f"diff = {abs(se-sx):.2e}   [CONTROL]")

line()
print("BLOCK 2D -- MUTUAL INFORMATION ACROSS THE SAME ENUMERATION.")
print("CHR's continuum claim is that the ambiguity is boundary-local, so RELATIVE ENTROPY and")
print("MUTUAL INFORMATION are unambiguous.  Finitely there is no continuum limit to take, so this")
print("is a live question here and not a citation.  Pairing rule: algebra B for A, and its")
print("commutant inside the gauge-invariant algebra for Abar.  I = S(B) + S(B^c) - S(B v B^c).")
line()
for g2 in (0.10, 0.50, 1.00, 3.00):
    psi, e0 = named_states(g2)
    cache = {}
    Ivals = []
    for b in subs:
        gens = list(Vmin)
        for bv in b:
            acc = 0
            for i in range(q):
                if (bv >> i) & 1:
                    acc ^= reps[i]
            gens.append(acc)
        B = PauliAlgebra(gens, n)
        Bc = PauliAlgebra(intersect(sympl_perp(B.basis, n), GI, n), n)
        Bj = PauliAlgebra(B.basis + Bc.basis, n)
        I = B.entropy(psi) + Bc.entropy(psi) - Bj.entropy(psi)
        Ivals.append(round(I, 9))
    u = sorted(set(Ivals))
    print(f"  g^2={g2:.2f}  {len(Ivals)} algebra choices -> {len(u)} distinct I(A:Abar); "
          f"min={min(u):.9f} max={max(u):.9f} SPREAD={max(u)-min(u):.9f} bits")
    print(f"          values: {[f'{v:.6f}' for v in u]}")

line()
print("END BLOCK 2")

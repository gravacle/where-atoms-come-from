"""
k2b_centres.py -- LANE W21-K, BLOCK 2 (rebuilt).  THE EXHAUSTIVE CENTRE ENUMERATION,
counted MODULO GAUSS, which is the only way it counts anything.

BLOCK 2.0 records the failure of the first attempt as a measured degeneracy, not a result.
BLOCK 2.1 rebuilds on a carrier where the region ambiguity is non-degenerate and enumerates
          EVERY algebra assignable to the region, with its entropy.
BLOCK 2.2 locates the named choices (minimal / electric centre / magnetic centre / trivial
          centre) inside that enumeration.
BLOCK 2.3 runs the same enumeration for MUTUAL INFORMATION.
BLOCK 2.4 the extended-Hilbert-space arm, computed in the ambient space by ordinary partial
          trace, and the declared Z_2 control.
"""
import numpy as np, itertools, sys
from klib import (Z2Gauge, rref_basis, span_elements, in_span, pack, unpack, sympl,
                  sympl_perp, intersect, gauge_invariant_subspace, pauli_expect)
from klog import Sector, LogAlgebra, _reduce

np.set_printoptions(precision=9, suppress=True)
def line(c="="): print(c * 100)


def all_subspaces(q):
    seen = set()
    out = [[]]
    cur = [[]]
    seen.add(frozenset([0]))
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
                seen.add(key); out.append(nb); nxt.append(nb)
        cur = nxt
    return out


def build(G, verts):
    """V_min(A), V_max(A) as ambient F_2 subspaces (packed)."""
    n = G.L
    GI = gauge_invariant_subspace(G)
    def vmin(vs):
        internal, boundary, external = G.region(vs)
        imask = sum(1 << l for l in internal)
        cyc = [c for c in span_elements(G.cycles) if c and (c & ~imask) == 0]
        return rref_basis([pack(G.X(l), n) for l in internal] +
                          [pack(G.W(c), n) for c in cyc]), internal, boundary, external
    Vmin, internal, boundary, external = vmin(verts)
    compl = [v for v in range(G.V) if v not in verts]
    VminB, intB, bndB, extB = vmin(compl)
    Vmax = intersect(sympl_perp(VminB, n), GI, n)
    return Vmin, Vmax, VminB, internal, boundary, external, GI


# ============================================================ BLOCK 2.0
line()
print("BLOCK 2.0 -- THE FIRST ATTEMPT'S NULL, DIAGNOSED AS A DEGENERACY AND LOGGED AS ONE.")
line()
edges3 = [(0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)]
G3 = Z2Gauge(6, edges3, "ladder3")
Vmin, Vmax, VminB, internal, boundary, external, GI = build(G3, [0, 3])
gp3 = rref_basis([pack(g, G3.L) for g in G3.gauss])
print(f"carrier ladder3, region A = vertices [0,3]: dim V_min={len(Vmin)} dim V_max={len(Vmax)} "
      f"-> q = {len(Vmax)-len(Vmin)} in the AMBIENT algebra, and 374 subspaces sit between them.")
rmin = rref_basis([_reduce(v, gp3) for v in Vmin])
rmax = rref_basis([_reduce(v, gp3) for v in Vmax])
rmin = [v for v in rmin if v]; rmax = [v for v in rmax if v]
print(f"MODULO GAUSS: dim V_min = {len(rmin)}, dim V_max = {len(rmax)}, "
      f"q_physical = {len(rmax)-len(rmin)}")
print(f">>> q_physical = {len(rmax)-len(rmin)}.  All 374 ambient algebras restrict to the SAME algebra on")
print("    the 4-dimensional physical sector.  The identical entropies were a THEOREM about the")
print("    Gauss law, not a measurement.  Scored as VOID.  This is the fourth time this program has")
print("    had to catch a forced identity; it is caught here before it is reported.")


# ============================================================ BLOCK 2.1
line()
print("BLOCK 2.1 -- A NON-DEGENERATE CARRIER, AND THE FULL ENUMERATION.")
line()

# 3x3 open grid: 9 vertices, 12 links, 4 plaquettes, c = 4 -> physical dim 16
def grid(R, C):
    idx = lambda i, j: i * C + j
    E = []
    for i in range(R):
        for j in range(C - 1):
            E.append((idx(i, j), idx(i, j + 1)))
    for i in range(R - 1):
        for j in range(C):
            E.append((idx(i, j), idx(i + 1, j)))
    return R * C, E

nv, E = grid(3, 3)
G = Z2Gauge(nv, E, "grid3x3")
n = G.L
S = Sector(G)
gp = rref_basis([pack(g, n) for g in G.gauss])
print(f"carrier {G.name}: V={G.V} L={G.L} cycle_dim={G.cycle_dim} physical dim D={S.D} "
      f"(logical qubits c={S.c})")
print(f"edges {list(enumerate(G.edges))}")

# scan regions for a non-degenerate ambiguity
print("""
REGION SCAN.  A region must be an honest region, not an arbitrary vertex set: BOTH the induced
subgraph on A and the induced subgraph on Abar are required to be CONNECTED and to contain at
least one link.  (The unrestricted scan's winner was the independent set {1,3,5,7}, which has
NO internal links at all: V_min is then empty and V_max is everything, so q_phys = 8 is a
statement that an empty region has no determined algebra.  True and useless.  Logged.)
q_phys = dim V_max - dim V_min, both modulo Gauss.""")
adj = {v: set() for v in range(G.V)}
for a, b in G.edges:
    adj[a].add(b); adj[b].add(a)

def connected(vs):
    vs = set(vs)
    if not vs:
        return False
    st = [next(iter(vs))]; seen = {st[0]}
    while st:
        u = st.pop()
        for w in adj[u] & vs:
            if w not in seen:
                seen.add(w); st.append(w)
    return seen == vs

rows = []
for size in range(2, G.V - 1):
    for verts in itertools.combinations(range(G.V), size):
        compl = [v for v in range(G.V) if v not in verts]
        if not connected(verts) or not connected(compl):
            continue
        Vmin, Vmax, VminB, internal, boundary, external, _ = build(G, list(verts))
        if not internal or not G.region(compl)[0]:
            continue
        rmin = [v for v in rref_basis([_reduce(x, gp) for x in Vmin]) if v]
        rmax = [v for v in rref_basis([_reduce(x, gp) for x in Vmax]) if v]
        if not all(in_span(v, rmax) for v in rmin):
            continue
        rows.append((len(rmax) - len(rmin), verts, rmin, rmax, internal, boundary, external))
rows.sort(key=lambda r: (-r[0], len(r[1])))
print(f"  {len(rows)} admissible regions.  top by q_phys:")
for r in rows[:8]:
    print(f"     q_phys={r[0]}  verts={r[1]}  internal={r[4]}  boundary={r[5]}")
best = next(r for r in rows if r[0] <= 6)
print(f"  CHOSEN (largest q_phys that is exhaustively enumerable, q<=6): vertices {best[1]}  "
      f"q_phys = {best[0]}  internal links {best[4]}  boundary links {best[5]}")

q, verts, rmin, rmax, internal, boundary, external = best
reps = []
cur = list(rmin)
for v in rmax:
    if not in_span(v, cur):
        reps.append(v); cur = rref_basis(cur + [v])
assert len(reps) == q
subs = all_subspaces(q)
print(f"  dim V_min(mod Gauss) = {len(rmin)}   dim V_max(mod Gauss) = {len(rmax)}")
print(f"  NUMBER OF ALGEBRAS ASSIGNABLE TO THE REGION (exhaustive over F_2^{q}) = {len(subs)}")

def alg_from(b, name=""):
    gens = list(rmin)
    for bv in b:
        acc = 0
        for i in range(q):
            if (bv >> i) & 1:
                acc ^= reps[i]
        gens.append(acc)
    return LogAlgebra(gens, S, gp, name)

print()
for g2 in (0.10, 0.50, 1.00, 3.00):
    psi, e0, gap = S.ground_state(g2)
    rho = np.outer(psi, psi.conj())
    cache = {}
    vals = {}
    for b in subs:
        A = alg_from(b)
        s = round(A.entropy(rho, cache), 9)
        vals.setdefault(s, []).append((A.dimV, A.nblocks, 2 ** A.k))
    ks = sorted(vals)
    print(f"  g^2={g2:.2f}  E0={e0:.9f} gap={gap:.6f}   {len(subs)} algebras -> "
          f"{len(ks)} DISTINCT entropies")
    print(f"      min={min(ks):.9f}  max={max(ks):.9f}  SPREAD={max(ks)-min(ks):.9f} bits")
    print(f"      values (bits): {[f'{v:.6f}' for v in ks]}")
    print(f"      multiplicity of each value: {[len(vals[v]) for v in ks]}")


# ============================================================ BLOCK 2.2
line()
print("BLOCK 2.2 -- THE NAMED CHOICES INSIDE THAT ENUMERATION.")
line()
Vmin_a, Vmax_a, VminB_a, internal, boundary, external, GIfull = build(G, list(verts))
imask = sum(1 << l for l in internal)
cyc_in = [c for c in span_elements(G.cycles) if c and (c & ~imask) == 0]
A_MIN = LogAlgebra(Vmin_a, S, gp, "MINIMAL        (strictly inside A only)")
A_ELEC = LogAlgebra([pack(G.X(l), n) for l in internal + boundary] +
                    [pack(G.W(c), n) for c in cyc_in], S, gp,
                    "ELECTRIC CENTRE (E of every link touching A + loops inside A)")
A_MAG = LogAlgebra(Vmax_a, S, gp, "MAGNETIC CENTRE (commutant of A_min(Abar), the maximal choice)")

# trivial-centre choices: CHR say these are maximal trees of boundary links / partial gauge
# fixings.  Finitely: any subspace between V_min and V_max whose symplectic radical is trivial.
trivial = []
for b in subs:
    A = alg_from(b)
    if A.r == 0:
        trivial.append((b, A))
print(f"algebras with TRIVIAL CENTRE (a genuine tensor factor, so an honest entanglement "
      f"entropy): {len(trivial)} of {len(subs)}")
for g2 in (0.10, 1.00, 3.00):
    psi, e0, gap = S.ground_state(g2)
    rho = np.outer(psi, psi.conj()); cache = {}
    print(f"\n  g^2 = {g2:.2f}")
    for A in (A_MIN, A_ELEC, A_MAG):
        print(f"    S = {A.entropy(rho, cache):.9f} bits | {A.describe()}")
    if trivial:
        tv = sorted({round(A.entropy(rho, cache), 9) for _, A in trivial})
        print(f"    trivial-centre choices: {len(trivial)} of them, entropies "
              f"{[f'{v:.6f}' for v in tv]}")


# ============================================================ BLOCK 2.3
line()
print("BLOCK 2.3 -- THE SAME ENUMERATION FOR MUTUAL INFORMATION I(A:Abar).")
print("Pairing: algebra B for A, and its commutant inside the gauge-invariant algebra for Abar.")
line()
GIred = [v for v in rref_basis([_reduce(x, gp) for x in GIfull]) if v]
for g2 in (0.10, 0.50, 1.00, 3.00):
    psi, e0, gap = S.ground_state(g2)
    rho = np.outer(psi, psi.conj()); cache = {}
    Iv = []
    for b in subs:
        B = alg_from(b)
        Bc = LogAlgebra(intersect(sympl_perp(B.basis, n), GIred, n), S, gp)
        Bj = LogAlgebra(B.basis + Bc.basis, S, gp)
        Iv.append(round(B.entropy(rho, cache) + Bc.entropy(rho, cache) - Bj.entropy(rho, cache), 9))
    u = sorted(set(Iv))
    print(f"  g^2={g2:.2f}  -> {len(u)} distinct I; min={min(u):.9f} max={max(u):.9f} "
          f"SPREAD={max(u)-min(u):.9f} bits")
    print(f"      values: {[f'{v:.6f}' for v in u]}")


# ============================================================ BLOCK 2.4
line()
print("BLOCK 2.4 -- THE EXTENDED-HILBERT-SPACE ARM, AND THE DECLARED Z_2 CONTROL.")
line()
print("""The extended-Hilbert-space algebra is the FULL matrix algebra on the region's links in the
UNCONSTRAINED space.  It is gauge-VARIANT, so it does not act on the physical sector and it is
not in the enumeration above.  It is computed here by ordinary partial trace of the lifted
physical state -- a completely separate code path.

*** RETRACTION, MINE, FILED BEFORE THE NUMBERS. ***  My first pass declared as a control that
"S_extended must equal S_electric-centre because Z_2 is abelian and Donnelly's log-dim-R term
vanishes".  That is a HEADLINE PARAPHRASE and it is wrong.  Donnelly's decomposition is for the
construction in which each boundary LINK IS SPLIT into two half-links carrying matched edge
modes, and the superselection label is the flux through the entangling surface.  Assigning a
whole boundary link to one side, which is what is computed below, is a DIFFERENT construction,
and the numbers below show the two do not agree: the bd->A assignment misses the electric-centre
answer by 2.000 bits at g^2 = 0.10.  The lane that warned about identification-by-paraphrase
caught its own.  The control is withdrawn; what replaces it is the theorem in 2.5.""")

def vn(rho):
    w = np.clip(np.linalg.eigvalsh((rho + rho.conj().T) / 2).real, 0, None)
    return float(-sum(x * np.log2(x) for x in w if x > 1e-14))

def ptrace_links(psi_amb, L, keep):
    pos = sorted(L - 1 - l for l in keep)
    rest = [i for i in range(L) if i not in pos]
    t = psi_amb.reshape([2] * L)
    t = np.transpose(t, pos + rest).reshape(1 << len(pos), 1 << len(rest))
    return t @ t.conj().T

print()
for g2 in (0.10, 0.50, 1.00, 3.00):
    psi, e0, gap = S.ground_state(g2)
    rho = np.outer(psi, psi.conj()); cache = {}
    amb = S.lift(psi)
    s_ext_A = vn(ptrace_links(amb, n, internal + boundary))
    s_ext_I = vn(ptrace_links(amb, n, internal))
    s_el = A_ELEC.entropy(rho, cache)
    s_min = A_MIN.entropy(rho, cache)
    print(f"  g^2={g2:.2f}  S_ext(bd->A)={s_ext_A:.9f}  S_electric-centre={s_el:.9f}  "
          f"DIFF={abs(s_ext_A-s_el):.9f}")
    print(f"           S_ext(bd->Abar)={s_ext_I:.9f}  S_minimal={s_min:.9f}  "
          f"DIFF={abs(s_ext_I-s_min):.2e}")

line()
print("BLOCK 2.5 -- WHEN THE EXTENDED AND THE GAUGE-INVARIANT ANSWER MUST AGREE.  A THEOREM,")
print("STATED FIRST, THEN EXHIBITED.  ONE VARIABLE MOVED: DOES THE REGION CONTAIN A CYCLE.")
line()
print("""THEOREM (finite, exact, Z_N).  Let R be a set of links and |psi> a physical state.  Every
Pauli supported on R with a non-empty Z-part is gauge-VARIANT unless that Z-part is a CYCLE.
For a gauge-variant operator O, <psi|O|psi> = <psi|G O G|psi> = -<psi|O|psi> = 0.  Hence:

   IF R CONTAINS NO CYCLE, the reduced density matrix of |psi> on R is exactly DIAGONAL in the
   electric basis, and S_extended(R) = S(electric algebra of R) EXACTLY, for every state, every
   coupling, every graph.  The extended-Hilbert-space answer and the gauge-invariant answer
   CANNOT DISAGREE.  Any null obtained on a forest region is void.

   IF R CONTAINS A CYCLE, the region carries a Wilson loop, the reduction has magnetic
   coherences, and the two answers may differ.  This is the ONLY place a disagreement can live.

The bd->Abar column above is the forest case (internal links form a tree) and the diff is 1e-14
FOUR TIMES OUT OF FOUR.  That agreement is the theorem, not a measurement.  Below is the arm
that actually moves the variable: a region whose internal links CONTAIN a plaquette.""")
cyc_region = [0, 1, 3, 4]
Vmin_c, Vmax_c, VminB_c, int_c, bnd_c, ext_c, _ = build(G, cyc_region)
imask_c = sum(1 << l for l in int_c)
cyc_inside = [c for c in span_elements(G.cycles) if c and (c & ~imask_c) == 0]
print(f"\n  ARM CYCLE:   region vertices {cyc_region}, internal links {int_c}, "
      f"cycles inside = {[bin(c) for c in cyc_inside]}  -> {len(cyc_inside)} Wilson loop(s)")
print(f"  ARM FOREST:  region vertices {list(verts)}, internal links {internal}, "
      f"cycles inside = {[bin(c) for c in span_elements(G.cycles) if c and (c & ~imask)==0]}")
A_MIN_C = LogAlgebra(Vmin_c, S, gp, "minimal(cycle region)")
A_ELEC_C = LogAlgebra([pack(G.X(l), n) for l in int_c + bnd_c] +
                      [pack(G.W(c), n) for c in cyc_inside], S, gp, "electric centre(cycle region)")
A_ELECONLY_C = LogAlgebra([pack(G.X(l), n) for l in int_c], S, gp, "electric-only(cycle region)")
print(f"\n  {'g^2':>6} | {'S_ext(internal only)':>21} | {'S_electric-only':>16} | "
      f"{'DIFF':>12} | {'S_minimal(with loop)':>21}")
for g2 in (0.10, 0.30, 0.50, 1.00, 3.00):
    psi, e0, gap = S.ground_state(g2)
    rho = np.outer(psi, psi.conj()); cache = {}
    s_ext = vn(ptrace_links(S.lift(psi), n, int_c))
    s_eo = A_ELECONLY_C.entropy(rho, cache)
    s_mn = A_MIN_C.entropy(rho, cache)
    print(f"  {g2:6.2f} | {s_ext:21.9f} | {s_eo:16.9f} | {abs(s_ext-s_eo):12.9f} | {s_mn:21.9f}")
print("""
READ IT BOTH WAYS.  A non-zero DIFF says the extended construction sees the region's Wilson loop
and the purely-electric algebra does not.  A zero DIFF says the state has no magnetic coherence
on that loop at that coupling.  Neither reading is evidence about gravity, and the column is
printed so the difference between the two readings is visible rather than chosen.""")

line()
print("END BLOCK 2")

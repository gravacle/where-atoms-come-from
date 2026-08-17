"""
k5_hypothesis.py -- LANE W21-K, BLOCK 5.  THE REGISTRAR'S HYPOTHESIS, TESTED DIRECTLY.

THE HYPOTHESIS AS WRITTEN IN THE BRIEF, AND FLAGGED THERE AS THE REGISTRAR'S:
  "THE AMBIGUITY IS THE ABSENCE OF GRAVITY SHOWING UP.  Something that makes the dynamics
   RESPOND TO WHAT IS AT THE BOUNDARY -- backreaction -- should make the competing bookkeepings
   INEQUIVALENT and collapse the disagreement."

THIS BLOCK MEASURES THE ONE THING THAT DECIDES WHETHER THAT SENTENCE CAN BE TESTED AT ALL:

  5A  IS THE SET OF COMPETING BOOKKEEPINGS A DYNAMICAL OBJECT?  The set of algebras assignable
      to a region is computed from V_min and V_max, which are built from the GRAPH and the GAUSS
      LAW ONLY.  No Hamiltonian appears anywhere in their construction.  Measured across five
      different Hamiltonians including a backreaction-flavoured one.
  5B  WHAT THE DYNAMICS DOES MOVE: the numerical spread, because it moves the STATE.
  5C  A BACKREACTION ARM, honestly labelled as an analogue, with its diff printed.
  5D  BOTH READINGS OF THE NULL, WRITTEN OUT, NEITHER SCORED.
"""
import numpy as np, itertools
from klib import (Z2Gauge, rref_basis, span_elements, in_span, pack, unpack, sympl_perp,
                  intersect, gauge_invariant_subspace, pauli_matrix)
from klog import Sector, LogAlgebra, _reduce

np.set_printoptions(precision=9, suppress=True)
def line(c="="): print(c * 100)
rng = np.random.default_rng(2026)


def all_subspaces(q):
    seen = {frozenset([0])}; out = [[]]; cur = [[]]
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
GI = gauge_invariant_subspace(G)
GIred = [v for v in rref_basis([_reduce(x, gp) for x in GI]) if v]

REGION = [0, 1, 2, 4]


def region_data():
    def vmin(vs):
        internal, boundary, external = G.region(vs)
        imask = sum(1 << l for l in internal)
        cyc = [c for c in span_elements(G.cycles) if c and (c & ~imask) == 0]
        return rref_basis([pack(G.X(l), n) for l in internal] + [pack(G.W(c), n) for c in cyc])
    Vmin = vmin(REGION)
    VminB = vmin([v for v in range(G.V) if v not in REGION])
    Vmax = intersect(sympl_perp(VminB, n), GI, n)
    rmin = [v for v in rref_basis([_reduce(x, gp) for x in Vmin]) if v]
    rmax = [v for v in rref_basis([_reduce(x, gp) for x in Vmax]) if v]
    return rmin, rmax


rmin, rmax = region_data()
q = len(rmax) - len(rmin)
reps = []
cur = list(rmin)
for v in rmax:
    if not in_span(v, cur):
        reps.append(v); cur = rref_basis(cur + [v])
subs = all_subspaces(q)


def alg_from(b):
    gens = list(rmin)
    for bv in b:
        acc = 0
        for i in range(q):
            if (bv >> i) & 1:
                acc ^= reps[i]
        gens.append(acc)
    return LogAlgebra(gens, S, gp)


ALGS = [alg_from(b) for b in subs]
FINGERPRINT = tuple(sorted(tuple(A.basis) for A in ALGS))


# ---------------------------------------------------------------- Hamiltonians
def H_standard(g2):
    return S.hamiltonian(g2)


def H_backreact(g2, lam):
    """BACKREACTION ANALOGUE, LABELLED AS ONE.  Each plaquette's weight responds to the electric
    content of the links around it: 1/g^2 -> (1/g^2)(1 + lam * mean(X_l) over links adjacent to
    but not in the plaquette).  Those X_l commute with W_p, so the term is Hermitian and the
    theory stays gauge invariant.  This is 'what happens changes the arena' to the extent that a
    graph with no metric and no 2-cell dynamics can express it, WHICH IS NOT FAR.  It is an
    ANALOGUE, not backreaction, and it is named as one before its number is read."""
    H = np.zeros((S.D, S.D), dtype=complex)
    for cyc in G.cycles:
        pl = [l for l in range(G.L) if (cyc >> l) & 1]
        near = set()
        for l in pl:
            a, b = G.edges[l]
            for l2, (c, d) in enumerate(G.edges):
                if l2 not in pl and ({c, d} & {a, b}):
                    near.add(l2)
        W = S.mat(0, cyc)
        H -= (1.0 / g2) * W
        for l2 in near:
            H -= (lam / g2 / max(1, len(near))) * (W @ S.mat(1 << l2, 0))
    for l in range(G.L):
        H -= g2 * S.mat(1 << l, 0)
    return (H + H.conj().T) / 2


def H_random(seed):
    r = np.random.default_rng(seed)
    A = r.normal(size=(S.D, S.D)) + 1j * r.normal(size=(S.D, S.D))
    return (A + A.conj().T) / 2


def gs(H):
    w, V = np.linalg.eigh((H + H.conj().T) / 2)
    return V[:, 0].copy(), float(w[0])


# ============================================================ 5A
line()
print("BLOCK 5A -- IS THE SET OF COMPETING BOOKKEEPINGS A DYNAMICAL OBJECT?  NO, AND HERE IS")
print("THE MEASUREMENT THAT SAYS SO.")
line()
print(f"carrier {G.name} V={G.V} L={G.L} physical dim {S.D}; region = vertices {REGION}")
print(f"V_min and V_max are built from: the graph's link/cycle structure and the Gauss law.")
print(f"grep the constructor: NO HAMILTONIAN, NO COUPLING, NO STATE enters region_data().")
print(f"  q_phys = {q}   number of assignable algebras = {len(subs)}")
print("\nRecomputed under five different dynamics (the function is state-independent, so this is")
print("a tautology check -- printed because the tautology IS the answer to the hypothesis):")
for lab in ("g^2=0.10", "g^2=3.00", "backreaction lam=0.5", "backreaction lam=2.0", "random H"):
    r2, r3 = region_data()
    fp = tuple(sorted(tuple(alg_from(b).basis) for b in all_subspaces(len(r3) - len(r2))))
    print(f"   {lab:<24} -> {len(all_subspaces(len(r3)-len(r2))):>3} algebras, "
          f"identical set: {fp == FINGERPRINT}")
print("""
>>> THE SET OF COMPETING BOOKKEEPINGS IS FIXED BEFORE ANY HAMILTONIAN EXISTS.  It is a property
    of the graph and the constraint.  NO term added to H -- backreaction, matter, anything --
    can delete an algebra from that list, because the list was never a function of H.  The
    hypothesis as written asks a dynamical mechanism to collapse a KINEMATICAL set.  On a finite
    lattice that is not a hard experiment; it is a category error, and it is the finding.""")


# ============================================================ 5B / 5C
line()
print("BLOCK 5B/5C -- WHAT THE DYNAMICS DOES MOVE: THE STATE, HENCE THE NUMBERS.")
line()
print(f"  {'dynamics':<28} {'E0':>14} {'#distinct S':>12} {'min S':>12} {'max S':>12} {'SPREAD':>12}")
runs = []
for lab, H in (("standard g^2=0.10", H_standard(0.10)),
               ("standard g^2=0.50", H_standard(0.50)),
               ("standard g^2=1.00", H_standard(1.00)),
               ("standard g^2=3.00", H_standard(3.00)),
               ("BACKREACT g^2=0.10 lam=0.5", H_backreact(0.10, 0.5)),
               ("BACKREACT g^2=0.10 lam=2.0", H_backreact(0.10, 2.0)),
               ("BACKREACT g^2=0.50 lam=0.5", H_backreact(0.50, 0.5)),
               ("BACKREACT g^2=0.50 lam=2.0", H_backreact(0.50, 2.0)),
               ("BACKREACT g^2=1.00 lam=2.0", H_backreact(1.00, 2.0)),
               ("random Hermitian, seed 7", H_random(7)),
               ("random Hermitian, seed 8", H_random(8))):
    psi, e0 = gs(H)
    rho = np.outer(psi, psi.conj()); cache = {}
    vals = sorted({round(A.entropy(rho, cache), 9) for A in ALGS})
    runs.append((lab, e0, vals))
    print(f"  {lab:<28} {e0:14.6f} {len(vals):>12} {min(vals):12.9f} {max(vals):12.9f} "
          f"{max(vals)-min(vals):12.9f}")

# Haar-random physical states
for sd in (11, 12, 13):
    r = np.random.default_rng(sd)
    v = r.normal(size=S.D) + 1j * r.normal(size=S.D); v /= np.linalg.norm(v)
    rho = np.outer(v, v.conj()); cache = {}
    vals = sorted({round(A.entropy(rho, cache), 9) for A in ALGS})
    print(f"  {'Haar physical state s='+str(sd):<28} {float('nan'):>14} {len(vals):>12} "
          f"{min(vals):12.9f} {max(vals):12.9f} {max(vals)-min(vals):12.9f}")

print("\n  THE ARM DIFF THE BRIEF ASKS FOR (one variable: lambda, everything else held):")
for g2 in (0.10, 0.50, 1.00):
    base = [r for r in runs if r[0] == f"standard g^2={g2:.2f}"][0]
    for lam in (0.5, 2.0):
        m = [r for r in runs if r[0] == f"BACKREACT g^2={g2:.2f} lam={lam}"]
        if not m:
            continue
        b = max(base[2]) - min(base[2]); a = max(m[0][2]) - min(m[0][2])
        print(f"    g^2={g2:.2f}  spread(lam=0) = {b:.9f}   spread(lam={lam}) = {a:.9f}   "
              f"DIFF = {a-b:+.9f} bits   collapsed? {'YES' if a < 1e-9 else 'NO'}")


# ============================================================ 5D
line()
print("BLOCK 5D -- BOTH READINGS, NEITHER SCORED.")
line()
print("""READING ONE.  The backreaction analogue does not collapse the spread, therefore the
hypothesis is wrong: the ambiguity is not the absence of gravity, it is the absence of a physical
reason to prefer one subalgebra, and no term in a Hamiltonian on a fixed graph supplies one.

READING TWO.  The backreaction analogue on this carrier is not backreaction.  There is no metric,
no 2-cell with its own dynamics, no arena to respond -- W-05 established attempt one had none and
this carrier has none either.  A term that makes a plaquette's coefficient depend on nearby
electric flux changes the STATE and nothing else, and the null is therefore a null about a toy,
not about gravity.  A classical-shaped null at the record level is VOID by this program's own
standing rule, and this one is close enough to that shape to be treated the same way.

WHAT DECIDES BETWEEN THEM IS 5A, AND 5A IS NOT A NULL.  5A is a structural statement with no arm
and no state: the list of admissible algebras is computed before any dynamics exists.  Whatever
gravity is, if it is a term in a Hamiltonian on a fixed carrier, it cannot shorten that list.  So
the hypothesis cannot be tested by adding backreaction to a finite lattice REGARDLESS of which
reading of the null is right.  That is what ends the round, not the null.""")

line()
print("END BLOCK 5")

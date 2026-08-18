"""O-8 PART A: are NON-CSS stabiliser codes inside G-12's class?
A CSS code IS a length-2 F_2 chain complex.  So the question 'is candidate C covered by G-12?'
is exactly 'is C a CSS code in SOME single-qubit-Clifford basis?'
Single-qubit Cliffords act on the symplectic pair (x,z) as GL(2,F_2) = S_3 (6 elements);
phases are irrelevant to CSS-ness.  We search that group EXHAUSTIVELY where feasible."""
import sys, itertools, random
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O78_FRACTON_NONABELIAN")
from gf2 import *
from models import steane, perfect5, toric2d, chamon

def hr(t=""): print("\n" + "="*78); print(t); print("="*78)

# GL(2,F_2) elements as (a,b,c,d): (x,z) -> (ax+bz, cx+dz). Index 0 = Hadamard, 2 = identity.
GL2 = []
for a in range(2):
    for b in range(2):
        for c in range(2):
            for d in range(2):
                if (a*d ^ b*c) == 1: GL2.append((a,b,c,d))
assert len(GL2) == 6

def apply_lc(gens, n, choice):
    """choice[i] indexes GL2 for qubit i.  (x,z) -> (a x + b z, c x + d z)."""
    out = []
    for g in gens:
        X = xpart(g,n); Z = zpart(g,n)
        nX = 0; nZ = 0
        for i in range(n):
            xi = (X >> i) & 1; zi = (Z >> i) & 1
            if not (xi or zi): continue
            a,b,c,d = GL2[choice[i]]
            if (a*xi ^ b*zi): nX |= (1 << i)
            if (c*xi ^ d*zi): nZ |= (1 << i)
        out.append(mk(nX,nZ,n))
    return out

def is_css_now(gens, n):
    return css_split(gens, n)[3]

def exhaustive_lc_search(name, gens, n, expect):
    total = 6**n
    print(f"{name}: exhaustive search over all {total} single-qubit-Clifford bases (6^{n})")
    found = None
    cnt = 0
    for choice in itertools.product(range(6), repeat=n):
        cnt += 1
        if is_css_now(apply_lc(gens, n, choice), n):
            found = choice; break
    verdict = "CSS FORM EXISTS" if found else "NO CSS FORM IN ANY LOCAL-CLIFFORD BASIS"
    print(f"    searched {cnt}/{total}.  RESULT: {verdict}")
    if found: print(f"    witness basis = {found}")
    ok = (bool(found) == expect)
    print(f"    SELF-CHECK (expected {'a CSS form' if expect else 'none'}): {'PASS' if ok else 'FAIL'}")
    return bool(found)

hr("PART A0 -- POSITIVE CONTROLS for the local-Clifford search")
print("""A search that can only ever answer NO is worthless.  Control 1 hides a CSS code
behind a random local Clifford: the search MUST find it again.""")
def code422():
    n = 4
    v = 0b1111
    return [mk(v,0,n), mk(0,v,n)], n
rng = random.Random(11)
cg4, cn4 = code422()
sc4 = tuple(rng.randrange(6) for _ in range(cn4))
print(f"  [[4,2,2]] scrambled by {sc4}: CSS in that basis? {is_css_now(apply_lc(cg4,cn4,sc4), cn4)}")
exhaustive_lc_search("  CONTROL scrambled [[4,2,2]]", apply_lc(cg4,cn4,sc4), cn4, expect=True)
sg, sn = steane()
scramble = tuple(rng.randrange(6) for _ in range(sn))
scrambled = apply_lc(sg, sn, scramble)
print(f"  Steane scrambled by local Clifford {scramble}: CSS in that basis? {is_css_now(scrambled, sn)}")
print(f"  k preserved under local Clifford? {code_k(scrambled,sn)[0]} == {code_k(sg,sn)[0]} : "
      f"{'PASS' if code_k(scrambled,sn)[0]==code_k(sg,sn)[0] else 'FAIL'}")
exhaustive_lc_search("  CONTROL scrambled Steane", scrambled, sn, expect=True)

print("\n  Control 2: the XZZX code (Bonilla Ataides et al. 2021) -- non-CSS as written,")
print("  known to be Hadamard-equivalent to a CSS surface code on alternate qubits.")
def xzzx_ring(L):
    """XZZX toric code on an L x L torus: qubits on VERTICES, stabiliser X Z X Z on each plaquette."""
    n = L*L
    idx = lambda x,y: (x%L)*L + (y%L)
    g = []
    for x in range(L):
        for y in range(L):
            X = (1<<idx(x,y)) | (1<<idx(x+1,y+1))
            Z = (1<<idx(x+1,y)) | (1<<idx(x,y+1))
            g.append(mk(X,Z,n))
    return g, n
xg, xn = xzzx_ring(4)
print(f"  XZZX L=4: n={xn} k={code_k(xg,xn)[0]} commute_violations={all_commute(xg,xn)} "
      f"CSS as written = {is_css_now(xg,xn)}")
IDENT = GL2.index((1,0,0,1)); HAD = GL2.index((0,1,1,0))
print(f"  (GL2 index: identity={IDENT}, Hadamard={HAD})")
had_choice = tuple((HAD if (i//4 + i%4)%2 else IDENT) for i in range(xn))   # H on one sublattice
print(f"  after Hadamard on the odd sublattice: CSS = {is_css_now(apply_lc(xg,xn,had_choice), xn)}")
print("  -> the search target is reachable; a NO answer elsewhere is therefore informative.")

hr("PART A1 -- THE PERFECT [[5,1,3]] CODE")
pg, pn = perfect5()
print(f"  n={pn} k={code_k(pg,pn)[0]} d=3 (exact, PART 0 of o7) commute_violations={all_commute(pg,pn)}")
found5 = exhaustive_lc_search("  [[5,1,3]]", pg, pn, expect=False)

hr("PART A2 -- DOES THE TWO-NUMBER STRUCTURE OF G-12 EVEN EXIST HERE?")
print("""G-12 requires TWO numbers: k-systole (min weight of a Z-logical) and k-cosystole
(min weight of an X-logical).  Those are only defined if the logical group SPLITS into
pure-X and pure-Z sectors.  Count pure-X and pure-Z logical representatives:""")
def sector_count(name, gens, n):
    _, pivS = rank_gf2(gens)
    N = normaliser_basis(gens, n)
    pureX = [v for v in N if zpart(v,n)==0]
    pureZ = [v for v in N if xpart(v,n)==0]
    # dimension of the pure-X part of S-perp
    rowsZ = []
    for bit in range(n):
        row = 0
        for i,v in enumerate(N):
            if (zpart(v,n)>>bit)&1: row |= (1<<i)
        rowsZ.append(row)
    cx = len(nullspace_basis(rowsZ, len(N)))
    rowsX = []
    for bit in range(n):
        row = 0
        for i,v in enumerate(N):
            if (xpart(v,n)>>bit)&1: row |= (1<<i)
        rowsX.append(row)
    cz = len(nullspace_basis(rowsX, len(N)))
    # how many of those are NON-trivial logicals (outside S)?
    def nontrivial(coefs):
        best = None
        for c in coefs:
            v = 0; cc=c; i=0
            while cc:
                if cc&1: v ^= N[i]
                cc >>= 1; i += 1
            if v and not in_span(pivS, v):
                w = pweight(v,n)
                if best is None or w < best: best = w
        return best
    bx = nontrivial(nullspace_basis(rowsZ, len(N)))
    bz = nontrivial(nullspace_basis(rowsX, len(N)))
    print(f"  {name:26s} dim(S-perp)={len(N):3d}  pure-X subspace dim={cx:3d}  pure-Z dim={cz:3d}  "
          f"min wt pure-X logical={bx}  pure-Z={bz}")
sector_count("2D toric L=3 (CSS)", *toric2d(3))
sector_count("Steane (CSS)", *steane())
sector_count("[[5,1,3]] (non-CSS)", *perfect5())
print("""
  READ THIS ROW CAREFULLY.  For [[5,1,3]] the pure-X and pure-Z logical sectors are
  1-dimensional (X^5 and Z^5), giving 'systole' = 'cosystole' = 5.  But the EXACT
  distance is 3, achieved only by a MIXED-type logical (e.g. Y-containing).  So on a
  non-CSS code the systole/cosystole pair does not compute the protection at all:
  min(systole, cosystole) = 5 while d = 3.""")

hr("PART A3 -- CHAMON MODEL under sublattice-uniform local Cliffords")
print("""Chamon is translation invariant, so a CSS-ifying Clifford would plausibly be
sublattice-uniform.  We search all assignments constant on the 2-, 4- and 8-site
sublattices of the cubic lattice (6^2 = 36, 6^4 = 1296, 6^8 = 1679616 bases).""")
for L in (4,):
    cg, cn = chamon(L)
    for mod, label in ((2,"parity (x+y+z) mod 2"), (4,"(x+2y+4z) mod 4"), (8,"(x%2,y%2,z%2) octant")):
        if mod == 2:   cls = lambda x,y,z: (x+y+z) % 2
        elif mod == 4: cls = lambda x,y,z: (x + 2*y + 4*z) % 4
        else:          cls = lambda x,y,z: (x%2) + 2*(y%2) + 4*(z%2)
        lab = [0]*cn
        for x in range(L):
            for y in range(L):
                for z in range(L):
                    lab[((x%L)*L+(y%L))*L+(z%L)] = cls(x,y,z)
        nc = max(lab)+1
        hit = None; tried = 0
        space = 6**nc
        if space <= 2000:
            it = itertools.product(range(6), repeat=nc)
        else:
            rr = random.Random(7)
            it = (tuple(rr.randrange(6) for _ in range(nc)) for _ in range(4000))
        for combo in it:
            tried += 1
            ch = [combo[lab[i]] for i in range(cn)]
            if is_css_now(apply_lc(cg, cn, ch), cn): hit = combo; break
        mode = "EXHAUSTIVE" if 6**nc <= 2000 else "SAMPLED"
        print(f"  L={L} classes={nc:2d} ({label:22s}) {mode} searched {tried:7d}/{6**nc:8d}  -> "
              f"{'CSS FOUND '+str(hit) if hit else 'NO CSS FORM in this family'}")
print("\nDONE.")

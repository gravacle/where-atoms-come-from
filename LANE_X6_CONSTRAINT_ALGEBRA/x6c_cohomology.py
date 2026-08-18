"""X6c.  WHAT THE RECORD SPACE *IS*, AS A QUOTIENT OF THE CONSTRAINT STRUCTURE.

The lane asks whether the record-level constraint is the SAME constraint as gravity's.  Before that
can be answered the record space has to be identified structurally, not by analogy.  This script
does exactly one thing: it exhibits the record space as an explicit quotient of two GF(2) vector
spaces built out of the constraint itself.

    C_2  --d2-->  C_1  --d1-->  C_0
    plaquettes    links         vertices

    Z_1 = ker d1   the CYCLE space         (claim: = the Gauss-law / physical sector)
    B_1 = im  d2   the BOUNDARY space      (claim: = the contractible-loop / no-record sector)
    H_1 = Z_1/B_1  the record space

Everything below is exact GF(2) linear algebra plus exhaustive enumeration.  No coupling constant
appears anywhere in this file, and no conclusion is drawn from any spectrum: the objects measured
are dimensions of kernels and images, which are scale-free properties of the algebra.

TOY_SEPARATION.  The lattice is a finite exact model of the constraint algebra.  The quantity of
interest (dim H_1) is reported at three sizes precisely so that its INDEPENDENCE of size can be
read off; any quantity that did move with size would be out of scope as a lattice artefact.
"""
import itertools, numpy as np

# ----------------------------------------------------------------------------------------------
# GF(2) linear algebra.  Rows are Python ints used as bitmasks over the link index set.
# ----------------------------------------------------------------------------------------------
def gf2_rowreduce(rows):
    """Return (list of pivot rows in echelon form, list of pivot bit positions)."""
    basis = []      # echelon rows
    piv = []        # pivot bit position of each echelon row
    for r in rows:
        cur = r
        for b, p in zip(basis, piv):
            if (cur >> p) & 1:
                cur ^= b
        if cur:
            p = cur.bit_length() - 1
            basis.append(cur); piv.append(p)
            # keep it sorted by pivot descending for determinism
            order = sorted(range(len(basis)), key=lambda k: -piv[k])
            basis = [basis[k] for k in order]; piv = [piv[k] for k in order]
    return basis, piv

def gf2_rank(rows):
    return len(gf2_rowreduce(rows)[0])

def gf2_reduce_vec(v, basis, piv):
    """Reduce vector v against an echelon basis; returns the residue (0 iff v in span)."""
    cur = v
    for b, p in zip(basis, piv):
        if (cur >> p) & 1:
            cur ^= b
    return cur

def gf2_kernel_basis(rows, ncols):
    """Basis (as bitmasks over ncols) of the kernel of the map x -> (rows . x), rows = list of
    bitmask rows of an (m x ncols) matrix.  Standard augmented row reduction."""
    # augment each COLUMN-tracking identity: work on the transpose.
    # Build m x ncols matrix as list of rows; solve M x = 0.
    M = list(rows)
    # echelon on columns 0..ncols-1
    piv_col_of_row = {}
    r = 0
    Mw = M[:]
    used = []
    for c in range(ncols):
        pr = None
        for i in range(r, len(Mw)):
            if (Mw[i] >> c) & 1:
                pr = i; break
        if pr is None:
            continue
        Mw[r], Mw[pr] = Mw[pr], Mw[r]
        for i in range(len(Mw)):
            if i != r and ((Mw[i] >> c) & 1):
                Mw[i] ^= Mw[r]
        piv_col_of_row[r] = c
        used.append(c)
        r += 1
        if r == len(Mw):
            break
    free = [c for c in range(ncols) if c not in used]
    ker = []
    for f in free:
        v = 1 << f
        for rr, c in piv_col_of_row.items():
            if (Mw[rr] >> f) & 1:
                v |= (1 << c)
        ker.append(v)
    return ker

def bits(mask, n):
    return tuple((mask >> k) & 1 for k in range(n))

def mask_of(tup):
    m = 0
    for k, v in enumerate(tup):
        if v: m |= (1 << k)
    return m

# ----------------------------------------------------------------------------------------------
# CANONICAL CARRIER (identical to LANE_G1 / LANE_G2 so the numbers are commensurable)
# ----------------------------------------------------------------------------------------------
def torus(nx, ny):
    vid = lambda i, j: (j % ny) * nx + (i % nx)
    E = []; ind = {}
    for j in range(ny):
        for i in range(nx): ind[('h', i, j)] = len(E); E.append((vid(i, j), vid(i + 1, j)))
    for j in range(ny):
        for i in range(nx): ind[('v', i, j)] = len(E); E.append((vid(i, j), vid(i, j + 1)))
    PL = [[ind[('h', i, j)], ind[('v', (i + 1) % nx, j)], ind[('h', i, (j + 1) % ny)], ind[('v', i, j)]]
          for j in range(ny) for i in range(nx)]
    return nx * ny, E, PL, ind

def open_grid(nx, ny):
    """Planar (genus 0, disk) nx x ny grid of plaquettes -- NO wrap-around.  Used only as an
    independently-known control: a disk has H_1 = 0."""
    VX, VY = nx + 1, ny + 1
    vid = lambda i, j: j * VX + i
    E = []; ind = {}
    for j in range(VY):
        for i in range(nx): ind[('h', i, j)] = len(E); E.append((vid(i, j), vid(i + 1, j)))
    for j in range(ny):
        for i in range(VX): ind[('v', i, j)] = len(E); E.append((vid(i, j), vid(i, j + 1)))
    PL = [[ind[('h', i, j)], ind[('v', i + 1, j)], ind[('h', i, j + 1)], ind[('v', i, j)]]
          for j in range(ny) for i in range(nx)]
    return VX * VY, E, PL, ind

def d1_rows(NV, E):
    """Boundary map d1 : C_1 -> C_0 over GF(2), one bitmask row per VERTEX over the link index set.
    Row v has bit k set iff link k is incident to v (an odd number of times)."""
    L = len(E)
    rows = []
    for v in range(NV):
        m = 0
        for k, (a, b) in enumerate(E):
            c = (1 if a == v else 0) + (1 if b == v else 0)
            if c % 2: m |= (1 << k)
        rows.append(m)
    return rows

def d2_rows(PL):
    """Plaquette map d2 : C_2 -> C_1 over GF(2), one bitmask row per PLAQUETTE."""
    return [mask_of_links(p) for p in PL]

def mask_of_links(links):
    m = 0
    for k in links: m ^= (1 << k)
    return m

# ----------------------------------------------------------------------------------------------
# THE MEASUREMENT
# ----------------------------------------------------------------------------------------------
def analyse(name, NV, E, PL, exhaustive=True):
    L = len(E); NP = len(PL)
    D1 = d1_rows(NV, E)
    D2 = d2_rows(PL)

    rank_d1 = gf2_rank(D1)
    dim_Z1 = L - rank_d1
    dim_B1 = gf2_rank(D2)

    Zb = gf2_kernel_basis(D1, L)                 # basis of Z_1
    Bb, Bp = gf2_rowreduce(D2)                   # echelon basis of B_1

    # d1 . d2 = 0 : every plaquette really is a cycle, so B_1 subset Z_1.
    d2d1_zero = True
    for p in D2:
        for rv in D1:
            if bin(p & rv).count('1') % 2 != 0:
                d2d1_zero = False
    dim_H1 = dim_Z1 - dim_B1

    print(f"\n=== {name}   links L={L}  vertices V={NV}  plaquettes F={NP} ===")
    print(f"  rank(d1)              = {rank_d1}      (connected graph predicts V-1 = {NV-1})")
    print(f"  dim Z_1 = ker(d1)     = {dim_Z1}      (graph formula L-V+1 = {L-NV+1})")
    print(f"  dim B_1 = im(d2)      = {dim_B1}      (torus predicts F-1 = {NP-1})")
    print(f"  d1 o d2 = 0 (B_1 subset Z_1)? {d2d1_zero}")
    print(f"  dim H_1 = dim Z_1 - dim B_1 = {dim_Z1} - {dim_B1} = {dim_H1}")
    print(f"  |Z_1| = 2^{dim_Z1} = {2**dim_Z1}   |B_1| = 2^{dim_B1} = {2**dim_B1}   "
          f"|H_1| = 2^{dim_H1} = {2**dim_H1}")

    # ---- explicit coset enumeration of Z_1 / B_1 ------------------------------------------
    Zelems = []
    for coeffs in itertools.product(range(2), repeat=dim_Z1):
        v = 0
        for c, b in zip(coeffs, Zb):
            if c: v ^= b
        Zelems.append(v)
    classes = {}
    for v in Zelems:
        r = gf2_reduce_vec(v, Bb, Bp)            # canonical coset representative
        classes.setdefault(r, []).append(v)
    sizes = sorted({len(c) for c in classes.values()})
    print(f"  explicit cosets of B_1 in Z_1: {len(classes)} classes, class sizes {sizes}")

    # ---- GAUSS LAW == Z_1 : exhaustive membership test -------------------------------------
    match_rate = None; nchecked = None; nphys = None
    if exhaustive:
        Zset = set(Zelems)
        agree = 0; total = 0; nphys = 0
        for s in itertools.product(range(2), repeat=L):
            # CANONICAL Gauss law, copied verbatim from the program construction
            g = all((sum(s[k] for k, (a, b) in enumerate(E) if a == v)
                     - sum(s[k] for k, (a, b) in enumerate(E) if b == v)) % 2 == 0
                    for v in range(NV))
            inZ = mask_of(s) in Zset
            if g: nphys += 1
            if g == inZ: agree += 1
            total += 1
        match_rate = agree / total
        nchecked = total
        print(f"  Gauss-law sector size D = {nphys}   (2^dim Z_1 = {2**dim_Z1})  "
              f"{'MATCH' if nphys == 2**dim_Z1 else 'MISMATCH'}")
        print(f"  membership test over ALL {total} link configurations: "
              f"'satisfies Gauss law' == 'lies in ker(d1)' agreement = {agree}/{total} "
              f"= {100.0*match_rate:.6f}%")

    return dict(name=name, L=L, NV=NV, NP=NP, rank_d1=rank_d1, dim_Z1=dim_Z1, dim_B1=dim_B1,
                dim_H1=dim_H1, nclasses=len(classes), class_sizes=sizes,
                match_rate=match_rate, nchecked=nchecked, nphys=nphys, d2d1_zero=d2d1_zero)

# ----------------------------------------------------------------------------------------------
# SELF-CHECKS AGAINST INDEPENDENTLY KNOWN ANSWERS
# ----------------------------------------------------------------------------------------------
def selfcheck():
    print("=" * 96)
    print("SELF-CHECK  (cases whose answers are known independently of this script)")
    print("=" * 96)
    ok = True

    # (a) GF(2) rank / kernel machinery on standard graphs: dim Z_1 = |E| - |V| + 1 (connected).
    graphs = {
        "path P_4 (tree)":      (4, [(0,1),(1,2),(2,3)],                       0),
        "cycle C_5":            (5, [(0,1),(1,2),(2,3),(3,4),(4,0)],           1),
        "complete K_4":         (4, [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)],     3),
        "theta graph":          (3, [(0,1),(0,1),(0,2),(1,2)],                 2),
    }
    for gname, (nv, ed, known) in graphs.items():
        got = len(ed) - gf2_rank(d1_rows(nv, ed))
        good = (got == known)
        ok &= good
        print(f"  [{'PASS' if good else 'FAIL'}] dim Z_1({gname}) = {got}, known = {known}")

    # (b) genus-0 control: a planar 2x2 grid of plaquettes is a DISK, so H_1 must be 0.
    NV, E, PL, _ = open_grid(2, 2)
    D1 = d1_rows(NV, E); D2 = d2_rows(PL)
    z = len(E) - gf2_rank(D1); b = gf2_rank(D2); h = z - b
    good = (h == 0)
    ok &= good
    print(f"  [{'PASS' if good else 'FAIL'}] planar 2x2 grid (disk, genus 0): "
          f"dim Z_1={z}, dim B_1={b}, dim H_1={h}, known = 0")

    # (c) Euler characteristic of the torus triangulation: V - E + F = 0.
    for (nx, ny) in [(2,2),(2,3),(3,3)]:
        NV, E, PL, _ = torus(nx, ny)
        chi = NV - len(E) + len(PL)
        good = (chi == 0)
        ok &= good
        print(f"  [{'PASS' if good else 'FAIL'}] Euler char torus {nx}x{ny}: "
              f"V-E+F = {NV}-{len(E)}+{len(PL)} = {chi}, known = 0")

    print(f"  SELF-CHECK BLOCK 1: {'PASS' if ok else 'FAIL'}")
    return ok

# ----------------------------------------------------------------------------------------------
def main():
    ok = selfcheck()
    if not ok:
        print("\nABORT: self-check block 1 failed; no headline number is reported.")
        return

    print("\n" + "=" * 96)
    print("MEASUREMENT: Z_1, B_1, H_1 ON Z2 TORI")
    print("=" * 96)
    res = []
    for (nx, ny) in [(2,2),(2,3),(3,3)]:
        NV, E, PL, ind = torus(nx, ny)
        r = analyse(f"torus {nx}x{ny}", NV, E, PL, exhaustive=True)
        r['nx'], r['ny'] = nx, ny
        # explicit non-contractible representatives: the two winding Wilson loops
        hloop = mask_of_links([ind[('h', i, 0)] for i in range(nx)])   # winds in x
        vloop = mask_of_links([ind[('v', 0, j)] for j in range(ny)])   # winds in y
        D1 = d1_rows(NV, E); D2 = d2_rows(PL)
        Bb, Bp = gf2_rowreduce(D2)
        def is_cycle(m): return all(bin(m & rv).count('1') % 2 == 0 for rv in D1)
        hc, vc = is_cycle(hloop), is_cycle(vloop)
        rh = gf2_reduce_vec(hloop, Bb, Bp)
        rv_ = gf2_reduce_vec(vloop, Bb, Bp)
        rhv = gf2_reduce_vec(hloop ^ vloop, Bb, Bp)
        print(f"  winding loops: x-loop in Z_1? {hc}   y-loop in Z_1? {vc}")
        print(f"  nontrivial in H_1 (residue != 0)? x-loop {rh != 0}, y-loop {rv_ != 0}, "
              f"x+y {rhv != 0}   -> {1 + (rh!=0) + (rv_!=0) + (rhv!=0)} distinct classes exhibited "
              f"(|H_1| = {2**r['dim_H1']})")
        r['loops_ok'] = hc and vc and rh != 0 and rv_ != 0 and rhv != 0
        res.append(r)

    # ------------------------------------------------------------------------------------------
    print("\n" + "=" * 96)
    print("5.  SIDE BY SIDE: DOES dim H_1 TRACK AREA, OR ONLY TOPOLOGY?")
    print("=" * 96)
    print(f"  {'lattice':<12}{'area(F)':>9}{'L':>5}{'V':>5}{'dim Z_1':>9}{'dim B_1':>9}{'dim H_1':>9}"
          f"{'|H_1|':>8}{'match rate':>13}")
    for r in res:
        print(f"  {r['name']:<12}{r['NP']:>9}{r['L']:>5}{r['NV']:>5}{r['dim_Z1']:>9}{r['dim_B1']:>9}"
              f"{r['dim_H1']:>9}{2**r['dim_H1']:>8}{100.0*r['match_rate']:>12.6f}%")
    areas = [r['NP'] for r in res]
    print(f"\n  area F ranges over {areas} (a factor of {max(areas)/min(areas):.2f})")
    print(f"  dim Z_1 = {[r['dim_Z1'] for r in res]}   -> MOVES with area (= F + 1)")
    print(f"  dim B_1 = {[r['dim_B1'] for r in res]}   -> MOVES with area (= F - 1)")
    print(f"  dim H_1 = {[r['dim_H1'] for r in res]}   -> CONSTANT")
    const = len({r['dim_H1'] for r in res}) == 1
    print(f"  dim H_1 depends on AREA? {'NO' if const else 'YES'}    "
          f"depends on TOPOLOGY only? {'YES' if const else 'NO'}")

    # ------------------------------------------------------------------------------------------
    print("\n" + "=" * 96)
    print("FINAL SELF-CHECK")
    print("=" * 96)
    all_ok = True
    for r in res:
        c1 = (r['dim_H1'] == 2)
        c2 = (r['dim_Z1'] - r['dim_B1'] == r['dim_H1'])
        c3 = (r['match_rate'] == 1.0)
        c4 = (r['nphys'] == 2 ** r['dim_Z1'])
        c5 = (r['nclasses'] == 2 ** r['dim_H1'])
        c6 = (r['class_sizes'] == [2 ** r['dim_B1']])
        c7 = r['d2d1_zero']
        c8 = r['loops_ok']
        all_ok &= (c1 and c2 and c3 and c4 and c5 and c6 and c7 and c8)
        print(f"  {r['name']}:")
        print(f"    [{'PASS' if c1 else 'FAIL'}] dim H_1 == 2 (= 2g, g=1)                 got {r['dim_H1']}")
        print(f"    [{'PASS' if c2 else 'FAIL'}] dim Z_1 - dim B_1 == dim H_1             "
              f"{r['dim_Z1']} - {r['dim_B1']} = {r['dim_Z1']-r['dim_B1']}")
        print(f"    [{'PASS' if c3 else 'FAIL'}] Gauss law <=> ker(d1), exhaustive        "
              f"rate {r['match_rate']}")
        print(f"    [{'PASS' if c4 else 'FAIL'}] |Gauss sector| == 2^dim Z_1              "
              f"{r['nphys']} vs {2**r['dim_Z1']}")
        print(f"    [{'PASS' if c5 else 'FAIL'}] #cosets == 2^dim H_1                     "
              f"{r['nclasses']} vs {2**r['dim_H1']}")
        print(f"    [{'PASS' if c6 else 'FAIL'}] every coset has 2^dim B_1 elements       "
              f"{r['class_sizes']} vs [{2**r['dim_B1']}]")
        print(f"    [{'PASS' if c7 else 'FAIL'}] d1 o d2 == 0")
        print(f"    [{'PASS' if c8 else 'FAIL'}] both winding loops are nontrivial in H_1")
    print(f"\n  OVERALL SELF-CHECK: {'PASS' if (all_ok and ok) else 'FAIL'}")

    # ------------------------------------------------------------------------------------------
    print("\n" + "=" * 96)
    print("6.  WHAT THE RECORD SPACE IS, IN TERMS OF THE CONSTRAINT")
    print("=" * 96)
    print("""  The Gauss-law constraint is not a filter applied to the record space; it IS the
  definition of the cycle space, and the record space is what is left of it after the
  plaquettes are divided out:

      record space  =  H_1(T^2; GF(2))  =  Z_1 / B_1
                    =  ker(Gauss law on link configurations) / im(plaquette map)
                    =  {configurations the constraint ADMITS} / {configurations the constraint
                       admits that are generated by contractible plaquette loops}

  ONE SENTENCE: the record space is the quotient of the constraint's KERNEL (the Gauss-law-
  satisfying link configurations, = the cycle space Z_1) by the constraint's own IMAGE from
  one level up (the span of the plaquettes, = the boundary space B_1), i.e. exactly the
  degree-1 GF(2) homology of the surface -- so a record is a constraint-admissible
  configuration MODULO everything the constraint can generate locally.""")

if __name__ == "__main__":
    main()

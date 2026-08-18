# LANE X6-a — EXACTLY WHAT ALGEBRA DO THE GAUSS CONSTRAINTS CLOSE INTO?
#
# The lattice is a FINITE EXACT MODEL OF AN ALGEBRA, not a stand-in for spacetime. Everything
# reported here is algebraic structure (generators, closure, rank, structure constants vs
# structure functions) -- scale-free properties that survive finiteness.
#
# G_v = prod_{k incident to v} sigma^x_k, built on the FULL 2^L link Hilbert space.
#
# TWO REPRESENTATIONS ARE BUILT AND CROSS-CHECKED, because "the computational basis" is a
# convention and the algebra must not depend on it:
#   REP A ("flip rep")  -- computational basis = sigma^z eigenbasis. G_v is a PERMUTATION matrix
#                          (it flips the links incident on v). Commutation here is a NON-TRIVIAL
#                          fact: permutation matrices generically do not commute.
#   REP B ("electric basis") -- computational basis = sigma^x eigenbasis. G_v is DIAGONAL with
#                          entries (-1)^(sum of incident bits). This is the basis in which the
#                          program's canonical physical-sector construction is literally the
#                          Gauss law, so the canonical snippet is run verbatim here.
# Commutation in REP B is trivially true (diagonal matrices commute); it is reported only as a
# consistency check on REP A, never as the evidence.
#
# Everything is exact: permutations are integers, matrices are int8, no floating-point tolerance
# is used anywhere a zero is claimed (max|entry| == 0 exactly).

import itertools
import numpy as np

rng = np.random.default_rng(20260818)

def rule(t):
    print()
    print("=" * 92)
    print(t)
    print("=" * 92)

def sub(t):
    print()
    print("-- " + t)


# ---------------------------------------------------------------------------------------------
# CARRIER
# ---------------------------------------------------------------------------------------------
def build_torus(nx, ny):
    """nx x ny periodic square lattice. Link ('h',i,j): (i,j)->(i+1,j). Link ('v',i,j): (i,j)->(i,j+1)."""
    ind = {}
    E = []
    def V(i, j):
        return (i % nx) + nx * (j % ny)
    k = 0
    for j in range(ny):
        for i in range(nx):
            ind[('h', i, j)] = k; E.append((V(i, j), V(i + 1, j))); k += 1
            ind[('v', i, j)] = k; E.append((V(i, j), V(i, j + 1))); k += 1
    NV = nx * ny
    L = len(E)
    # CANONICAL TORUS PLAQUETTES (program convention, note the wrap-around)
    PL = [[ind[('h', i, j)], ind[('v', (i + 1) % nx, j)], ind[('h', i, (j + 1) % ny)], ind[('v', i, j)]]
          for j in range(ny) for i in range(nx)]
    return ind, E, PL, NV, L


def incidence_masks(E, NV, L):
    """mask_v = bitmask of links k such that sigma^x_k appears an ODD number of times in G_v.
       (A self-loop would contribute sigma^x twice = identity; handled correctly by the parity.)"""
    masks = []
    for v in range(NV):
        m = 0
        for k, (a, b) in enumerate(E):
            cnt = (1 if a == v else 0) + (1 if b == v else 0)
            if cnt % 2 == 1:
                m |= (1 << k)
        masks.append(m)
    return masks


# ---------------------------------------------------------------------------------------------
# GF(2) LINEAR ALGEBRA ON THE GENERATORS  (exact, integer)
# ---------------------------------------------------------------------------------------------
def gf2_pivots(masks):
    """Row-reduce the vertex-link incidence matrix over GF(2), tracking which generator
       combination produced each pivot row. Returns dict pivot_bit -> (row_value, combo)."""
    piv = {}
    for u, val in enumerate(masks):
        v, c = val, (1 << u)
        for p in sorted(piv, reverse=True):
            if (v >> p) & 1:
                v ^= piv[p][0]; c ^= piv[p][1]
        if v:
            piv[v.bit_length() - 1] = (v, c)
    return piv


def gf2_solve(piv, m):
    """Return a combo c (bitmask over generators) with XOR_{u in c} mask_u == m, or None."""
    v, c = m, 0
    for p in sorted(piv, reverse=True):
        if (v >> p) & 1:
            v ^= piv[p][0]; c ^= piv[p][1]
    return None if v else c


def popcount(x):
    return bin(x).count("1")


# ---------------------------------------------------------------------------------------------
# SELF-CHECKS THAT MUST BE ABLE TO FAIL  (positive controls)
# ---------------------------------------------------------------------------------------------
def positive_controls(L):
    """The commutator machinery must report NON-ZERO when the answer is known to be non-zero,
       and the GF(2) rank machinery must report NV-c, not a hardcoded NV-1."""
    ok = True
    sub("POSITIVE CONTROL 1 — the norm machinery can see a non-zero commutator")
    Xp = np.array([[0, 1], [1, 0]], dtype=float)
    Zp = np.array([[1, 0], [0, -1]], dtype=float)
    I2 = np.eye(2)
    def kr(ops):
        M = np.array([[1.0]])
        for o in ops:
            M = np.kron(M, o)
        return M
    X0 = kr([Xp] + [I2] * (L - 1))
    Z0 = kr([Zp] + [I2] * (L - 1))
    C = X0 @ Z0 - Z0 @ X0
    got = float(np.linalg.norm(C, 'fro'))
    exp = 2.0 * np.sqrt(2.0) * (2.0 ** ((L - 1) / 2.0))
    good = abs(got - exp) < 1e-9
    ok &= good
    print(f"   || [sigma^x_0 , sigma^z_0] ||_F  on 2^{L} =  {got:.10f}   (independent closed form "
          f"2*sqrt(2)*2^((L-1)/2) = {exp:.10f})   {'PASS' if good else 'FAIL'}")

    sub("POSITIVE CONTROL 2 — the GF(2) rank machinery is not hardcoding NV-1")
    # Two DISJOINT edges: NV=4, L=2, c=2 components. Incidence rank over GF(2) = NV - c = 2,
    # so the generated group is Z2^2 and prod_v G_v = I is NOT the only relation.
    E2 = [(0, 1), (2, 3)]
    m2 = incidence_masks(E2, 4, 2)
    r2 = len(gf2_pivots(m2))
    good2 = (r2 == 2)
    ok &= good2
    print(f"   two disjoint edges (NV=4, L=2, 2 components): GF(2) rank = {r2}   "
          f"(known independently: NV - #components = 4 - 2 = 2)   {'PASS' if good2 else 'FAIL'}")
    print(f"   -> the code returns NV-2 here, so a report of NV-1 on the torus is a measurement, "
          f"not a hardcode.")
    return ok


# ---------------------------------------------------------------------------------------------
# THE MAIN BATTERY
# ---------------------------------------------------------------------------------------------
def run(nx, ny, do_dense_spectral, n_states):
    ind, E, PL, NV, L = build_torus(nx, ny)
    dim = 2 ** L
    masks = incidence_masks(E, NV, L)
    idx = np.arange(dim, dtype=np.int64)

    rule(f"CARRIER:  {nx} x {ny} TORUS   —   NV = {NV} vertices,  L = {L} links,  "
         f"{len(PL)} plaquettes,  full Hilbert space dim 2^{L} = {dim}")
    print(f"   links E = {E}")
    print(f"   G_v flip masks (bit k set <=> sigma^x_k appears in G_v):")
    for v in range(NV):
        print(f"      G_{v}: links {[k for k in range(L) if (masks[v]>>k)&1]}   mask = {masks[v]:0{L}b}")

    # ------------------------------------------------------------------ SELF-CHECK: D
    rule("SELF-CHECK (required):  D  ==  2^(L - NV + 1)")
    st = [s for s in itertools.product(range(2), repeat=L)
          if all((sum(s[k] for k, (a, b) in enumerate(E) if a == v)
                  - sum(s[k] for k, (a, b) in enumerate(E) if b == v)) % 2 == 0 for v in range(NV))]
    idxmap = {s: i for i, s in enumerate(st)}
    D = len(st)
    Dexp = 2 ** (L - NV + 1)
    p1 = (D == Dexp)
    print(f"   canonical physical-sector construction (program snippet, verbatim):  D = {D}")
    print(f"   independently known torus value 2^(L-NV+1) = 2^({L}-{NV}+1) = {Dexp}")
    print(f"   ROUTE 1 (REP B, G_v diagonal, physical = Gauss-satisfying configurations):  "
          f"{'PASS' if p1 else 'FAIL'}")

    # independent route: dimension of the common +1 eigenspace of the FLIP rep = # of group orbits
    piv = gf2_pivots(masks)
    rank = len(piv)
    group_masks = set()
    for S in range(1 << NV):
        m = 0
        for u in range(NV):
            if (S >> u) & 1:
                m ^= masks[u]
        group_masks.add(m)
    gorder = len(group_masks)
    seen = np.zeros(dim, dtype=bool)
    orbits = 0
    gm = sorted(group_masks)
    for s in range(dim):
        if not seen[s]:
            orbits += 1
            for m in gm:
                seen[s ^ m] = True
    p2 = (orbits == Dexp)
    print(f"   ROUTE 2 (REP A, G_v permutations, physical = common +1 eigenspace = # group orbits, "
          f"the action is free so each orbit gives exactly one invariant vector):  dim = {orbits}   "
          f"{'PASS' if p2 else 'FAIL'}")
    selfcheck = p1 and p2
    print(f"   SELF-CHECK: {'PASS' if selfcheck else 'FAIL'}")
    if not selfcheck:
        raise SystemExit("SELF-CHECK FAILED — no headline number will be reported.")

    # ------------------------------------------------------------------ EXPLICIT MATRICES
    rule("EXPLICIT MATRICES ON THE FULL 2^L SPACE  —  built TWICE, by independent routes")
    def perm_matrix(mask):
        P = np.zeros((dim, dim), dtype=np.int8)
        P[idx ^ mask, idx] = 1
        return P
    Xp8 = np.array([[0, 1], [1, 0]], dtype=np.int8)
    I28 = np.eye(2, dtype=np.int8)
    agree = True
    for v in range(NV):
        M = np.array([[1]], dtype=np.int8)
        for k in range(L - 1, -1, -1):                    # link 0 = least significant bit
            M = np.kron(M, Xp8 if (masks[v] >> k) & 1 else I28)
        Pv = perm_matrix(masks[v])
        same = bool(np.array_equal(M, Pv))
        agree &= same
        del M
    print(f"   REP A  G_v built as an explicit Kronecker product of {L} Pauli factors "
          f"(sigma^x on incident links, identity elsewhere), for all {NV} vertices,")
    print(f"   compared entrywise against the combinatorial bit-flip permutation matrix: "
          f"IDENTICAL for every v = {agree}")
    diagB = [np.array([(-1) ** popcount(s & masks[v]) for s in range(dim)], dtype=np.int8)
             for v in range(NV)]
    print(f"   REP B  G_v = diag((-1)^(# incident links excited)) on 2^{L}: built for all {NV} "
          f"vertices.")

    # ------------------------------------------------------------------ 1. COMMUTATORS
    rule("1.  COMMUTATOR NORMS  || [G_v , G_w] ||  FOR EVERY PAIR")
    print("   REP A (permutation rep — the non-trivial one). Products formed by exact permutation")
    print("   composition, then materialised as explicit int8 matrices; the commutator is an exact")
    print("   integer matrix, so max|entry| == 0 means EVERY norm is exactly 0.")
    print()
    print(f"   {'pair':>10} {'max|entry|':>12} {'Frobenius':>14} {'spectral':>14}")
    maxfro = 0.0
    maxspec = 0.0
    maxent = 0
    for v in range(NV):
        for w in range(v + 1, NV):
            Pvw = perm_matrix(masks[v] ^ masks[w])          # G_v G_w
            Pwv = perm_matrix(masks[w] ^ masks[v])          # G_w G_v
            C = Pvw.astype(np.int16) - Pwv.astype(np.int16)
            me = int(np.abs(C).max())
            fro = float(np.sqrt(float((C.astype(np.int64) ** 2).sum())))
            if me == 0:
                spec = 0.0                                   # C is the exact zero matrix
                spec_s = "0.0 (exact)"
            elif do_dense_spectral:
                spec = float(np.linalg.norm(C.astype(float), 2)); spec_s = f"{spec:.10f}"
            else:
                spec = float('nan'); spec_s = "not computed"
            maxent = max(maxent, me); maxfro = max(maxfro, fro)
            if spec == spec:
                maxspec = max(maxspec, spec)
            print(f"   {'[G%d,G%d]'%(v,w):>10} {me:>12d} {fro:>14.10f} {spec_s:>14}")
            del Pvw, Pwv, C
    print()
    print(f"   MAXIMUM OVER ALL {NV*(NV-1)//2} PAIRS:  max|entry| = {maxent}   "
          f"max ||[G_v,G_w]||_F = {maxfro:.10f}   max ||[G_v,G_w]||_2 = {maxspec:.10f}")
    # rep B cross-check
    mb = 0
    for v in range(NV):
        for w in range(v + 1, NV):
            mb = max(mb, int(np.abs(diagB[v] * diagB[w] - diagB[w] * diagB[v]).max()))
    print(f"   REP B cross-check (diagonal rep, commutation is automatic there): "
          f"max|entry| over all pairs = {mb}")

    # ------------------------------------------------------------------ 2. INVOLUTION
    rule("2.  G_v^2  —  IS EACH GENERATOR AN INVOLUTION?")
    allinv = True
    for v in range(NV):
        sq_mask = masks[v] ^ masks[v]
        Psq = perm_matrix(sq_mask)
        isI = bool(np.array_equal(Psq, np.eye(dim, dtype=np.int8)))
        devB = int(np.abs(diagB[v] * diagB[v] - 1).max())
        allinv &= (isI and devB == 0)
        print(f"   G_{v}^2 = I ?   REP A: max|G_v^2 - I| = {int(np.abs(Psq - np.eye(dim,dtype=np.int8)).max())}"
              f"   REP B: max|G_v^2 - I| = {devB}   -> {isI and devB==0}")
        del Psq
    print(f"   EVERY generator is an involution: {allinv}   (order of each G_v = 2)")

    # ------------------------------------------------------------------ 3. GLOBAL RELATION
    rule("3.  THE PRODUCT OF ALL G_v OVER THE LATTICE")
    total = 0
    for v in range(NV):
        total ^= masks[v]
    Ptot = perm_matrix(total)
    prod_isI = bool(np.array_equal(Ptot, np.eye(dim, dtype=np.int8)))
    dB = 1
    for v in range(NV):
        dB = dB * diagB[v]
    devB = int(np.abs(dB - 1).max())
    print(f"   REP A: prod_v G_v is the flip by mask {total:0{L}b};  max|prod_v G_v - I| = "
          f"{int(np.abs(Ptot - np.eye(dim,dtype=np.int8)).max())}  ->  identity = {prod_isI}")
    print(f"   REP B: max|prod_v G_v - I| = {devB}")
    print(f"   REASON (structural, not numeric): every link has exactly two endpoints, so each")
    print(f"   sigma^x_k occurs exactly twice in prod_v G_v and cancels. This is a scale-free")
    print(f"   property of ANY graph, independent of lattice size.")
    print(f"   => THE {NV} GENERATORS ARE NOT INDEPENDENT: there is exactly one global relation.")
    print(f"   EXACT RANK of the generated group (GF(2) rank of the vertex-link incidence matrix)")
    print(f"   = {rank}   vs NV = {NV}   -> corank {NV - rank}")
    # exhibit the full relation lattice: kernel of c -> XOR_u c_u mask_u
    kernel = []
    for S in range(1 << NV):
        m = 0
        for u in range(NV):
            if (S >> u) & 1:
                m ^= masks[u]
        if m == 0:
            kernel.append(S)
    print(f"   COMPLETE relation subgroup {{ S subset of vertices : prod_{{v in S}} G_v = I }} has "
          f"{len(kernel)} elements: {[format(S,'0%db'%NV) for S in kernel]}")
    print(f"   (i.e. the ONLY relations are the trivial one and the product of ALL vertices)")
    del Ptot

    # ------------------------------------------------------------------ 4. THE GROUP
    rule("4.  THE GROUP GENERATED BY {G_v}")
    k = int(round(np.log2(gorder)))
    print(f"   |group| enumerated explicitly over all 2^NV = {1<<NV} words: {gorder} distinct operators")
    print(f"   |group| = 2^k with k = {k};   NV = {NV};   k = NV - {NV - k}")
    print(f"   GF(2) rank of the incidence matrix = {rank}  (agrees with k: {rank == k})")
    print(f"   Every element is an involution and all elements commute (measured above), and the")
    print(f"   group is generated by {NV} involutions modulo exactly {len(kernel)-1} nontrivial relation.")
    print()
    print(f"   >>> THE GROUP IS  Z2^(NV-1) = Z2^{NV-1}  of order {2**(NV-1)}  <<<   "
          f"(measured order {gorder}, match = {gorder == 2**(NV-1)})")
    print(f"   It is ABELIAN and ELEMENTARY ABELIAN (exponent 2). As an F2-vector space it is the")
    print(f"   row space of the incidence matrix; its dimension NV-1 is the graph-theoretic")
    print(f"   statement 'connected graph', which is scale-free.")

    # ------------------------------------------------------------------ 5. STRUCTURE CONSTANTS
    rule("5.  STRUCTURE CONSTANTS vs STRUCTURE FUNCTIONS  —  THE DECISIVE TEST")
    print("   THE QUESTION: G_v G_w expanded in the generator basis,")
    print("       G_v G_w  =  prod_u G_u^{c_u(v,w)}   with c_u in F2,")
    print("   are the c_u the SAME for every state, or do they depend on the field configuration?")
    print("   (A structure-FUNCTION algebra -- e.g. the ADM hypersurface-deformation algebra, whose")
    print("    bracket of two Hamiltonian constraints returns the momentum constraint with a")
    print("    coefficient built from the inverse spatial metric -- has coefficients that change")
    print("    from field configuration to field configuration. That is what we are looking for.)")
    print()
    print("   METHOD (state-resolved, not assertion): for each ordered pair (v,w) and each sampled")
    print("   state |psi>, compute y = G_v G_w |psi> by applying the two operators in sequence, then")
    print("   test EVERY one of the 2^NV generator words g against y. The set of words matching on")
    print("   THAT state is read off from the state itself; the coefficients c_u are then whatever")
    print("   the state says they are. We compare those read-off coefficients across states.")
    print()

    # ---- word table
    words = []
    for S in range(1 << NV):
        m = 0
        for u in range(NV):
            if (S >> u) & 1:
                m ^= masks[u]
        words.append((S, m))

    def coeff_readoff(PSI, label, check_coset=False):
        """PSI: (nstates, dim) real array. Returns per-pair list of per-state canonical coefficient
           vectors, read off from the state by exhaustive matching.
           If check_coset, also verify the matching word set on EVERY state is exactly the
           predicted coset { e_v + e_w + K : K in the relation subgroup }."""
        nst = PSI.shape[0]
        rows = []
        maxvar = 0
        coset_ok = True
        for v in range(NV):
            for w in range(NV):
                if v == w:
                    continue
                tmp = PSI[:, idx ^ masks[w]]          # G_w |psi>
                y = tmp[:, idx ^ masks[v]]            # G_v G_w |psi>
                match = np.zeros((nst, 1 << NV), dtype=bool)
                for (S, m) in words:
                    z = PSI[:, idx ^ m]
                    match[:, S] = np.all(z == y, axis=1)
                # canonical read-off: smallest matching word, as a coefficient vector over F2
                cvec = np.zeros((nst, NV), dtype=np.int8)
                nomatch = 0
                for i in range(nst):
                    ms = np.nonzero(match[i])[0]
                    if len(ms) == 0:
                        nomatch += 1
                        continue
                    S = int(ms.min())
                    for u in range(NV):
                        cvec[i, u] = (S >> u) & 1
                var = int(np.abs(cvec - cvec[0]).max()) if nst > 0 else 0
                nmatch = match.sum(axis=1)
                if check_coset:
                    S0 = (1 << v) | (1 << w)
                    pred = np.zeros(1 << NV, dtype=bool)
                    for K in kernel:
                        pred[S0 ^ K] = True
                    for i in range(nst):
                        if not np.array_equal(match[i], pred):
                            coset_ok = False
                rows.append((v, w, cvec[0].copy(), var, int(nmatch.min()), int(nmatch.max()), nomatch))
                maxvar = max(maxvar, var)
        return rows, maxvar, coset_ok

    # ---- (a) generic states on the FULL space: the sharp test
    sub(f"(a) {n_states} RANDOM GENERIC STATES ON THE FULL 2^{L} SPACE  (the sharp test: on a")
    print("    generic state, a word matches only if it IS the right group element)")
    PSI = rng.normal(size=(n_states, dim))
    rowsA, varA, cosetA = coeff_readoff(PSI, "generic", check_coset=True)
    print(f"   {'pair (v,w)':>12}  {'c read off from the states':<28} {'state-to-state var':>19} "
          f"{'#matching words':>17}")
    for (v, w, c0, var, nmn, nmx, nm) in rowsA:
        print(f"   {'(%d,%d)'%(v,w):>12}  {''.join(str(int(x)) for x in c0):<28} {var:>19d} "
              f"{('%d..%d'%(nmn,nmx)):>17}")
    print(f"   MAXIMUM STATE-TO-STATE VARIATION IN THE STRUCTURE COEFFICIENTS = {varA}")
    print(f"   Reading the table: the coefficient string is indexed by vertex u = 0..{NV-1}. The")
    print(f"   matching set on every state has exactly {rowsA[0][4]} elements = |relation subgroup| = "
          f"{len(kernel)},")
    print(f"   namely e_v + e_w AND its complement e_v + e_w + (1,1,...,1) -- the two are the same")
    print(f"   operator because prod_v G_v = I. The canonical read-off prints whichever of the two")
    print(f"   has the smaller integer code, which is why e.g. (0,{NV-1}) prints as the complement.")
    print(f"   EXPLICIT COSET CHECK: on EVERY sampled state the matching word set equals exactly")
    print(f"   {{ e_v + e_w + K : K in the relation subgroup }} for EVERY ordered pair:  {cosetA}")
    print(f"   -> the composition law is not merely constant-looking, it is the predicted coset on")
    print(f"      each individual state, with no state-dependent members.")

    # ---- (b) random states inside the physical sector
    sub(f"(b) {n_states} RANDOM PHYSICAL STATES")
    print("    REP A physical sector = the common +1 eigenspace of all G_v = gauge-invariant states,")
    print(f"    dimension {orbits}. Built as normalised group-orbit sums with random coefficients.")
    reps = []
    seen2 = np.zeros(dim, dtype=bool)
    for s in range(dim):
        if not seen2[s]:
            reps.append(s)
            for m in gm:
                seen2[s ^ m] = True
    B = np.zeros((len(reps), dim))
    for i, s in enumerate(reps):
        for m in gm:
            B[i, s ^ m] = 1.0
    coef = rng.normal(size=(n_states, len(reps)))
    PHYS = coef @ B
    # verify they really are physical
    worst = 0.0
    for v in range(NV):
        worst = max(worst, float(np.abs(PHYS[:, idx ^ masks[v]] - PHYS).max()))
    print(f"    verification that the sampled states are physical: max_v max|G_v psi - psi| = {worst:.3e}")
    rowsB, varB, _ = coeff_readoff(PHYS, "physical")
    nmn = min(r[4] for r in rowsB); nmx = max(r[5] for r in rowsB)
    print(f"    MAXIMUM STATE-TO-STATE VARIATION IN THE STRUCTURE COEFFICIENTS = {varB}")
    print(f"    #matching words on physical states ranges {nmn}..{nmx} -- on a gauge-INVARIANT state")
    print(f"    every group element acts as the identity, so all {1<<NV} words match. That is")
    print(f"    degenerate BY CONSTRUCTION and is why (a) is the sharp test; the point here is that")
    print(f"    the matching set is still the SAME set for every physical state (variation {varB}).")

    # ---- (c) REP B: random superpositions over the canonical physical configurations
    sub(f"(c) {n_states} RANDOM SUPERPOSITIONS OVER THE CANONICAL PHYSICAL CONFIGURATIONS (REP B)")
    PB = np.zeros((n_states, dim))
    phys_idx = np.array([sum(s[k] << k for k in range(L)) for s in st], dtype=np.int64)
    PB[:, phys_idx] = rng.normal(size=(n_states, len(phys_idx)))
    worddiag = []
    for (S, m) in words:
        g = np.ones(dim, dtype=np.int8)
        for u in range(NV):
            if (S >> u) & 1:
                g = g * diagB[u]
        worddiag.append(g)
    varC = 0
    for v in range(NV):
        for w in range(NV):
            if v == w:
                continue
            y = diagB[v] * (diagB[w] * PB)
            matchmat = np.zeros((n_states, 1 << NV), dtype=bool)
            for (S, m) in words:
                matchmat[:, S] = np.all(worddiag[S] * PB == y, axis=1)
            cv = np.zeros((n_states, NV), dtype=np.int8)
            for i in range(n_states):
                nz = np.nonzero(matchmat[i])[0]
                S = int(nz.min())
                for u in range(NV):
                    cv[i, u] = (S >> u) & 1
            varC = max(varC, int(np.abs(cv - cv[0]).max()))
    print(f"    MAXIMUM STATE-TO-STATE VARIATION IN THE STRUCTURE COEFFICIENTS = {varC}")
    print(f"    (note: REP B states are supported on the physical configurations, where every G_v")
    print(f"     acts as +1; as in (b) the matching set is maximal but IDENTICAL across states.)")

    # ---- (d) operator-level identity: the coefficients are exact as MATRICES
    sub("(d) OPERATOR-LEVEL IDENTITY (state-independence in its strongest form)")
    worstop = 0
    for v in range(NV):
        for w in range(NV):
            if v == w:
                continue
            c = gf2_solve(piv, masks[v] ^ masks[w])
            m = 0
            for u in range(NV):
                if (c >> u) & 1:
                    m ^= masks[u]
            A1 = perm_matrix(masks[v] ^ masks[w])
            A2 = perm_matrix(m)
            worstop = max(worstop, int(np.abs(A1.astype(np.int16) - A2.astype(np.int16)).max()))
            del A1, A2
    print(f"    max over all {NV*(NV-1)} ordered pairs of  max| G_v G_w - prod_u G_u^{{c_u}} |  =  {worstop}")
    print(f"    The expansion holds as an identity of MATRICES on the whole 2^{L} space, so the")
    print(f"    coefficients cannot depend on the state: there is no state left for them to depend on.")

    maxvar_all = max(varA, varB, varC)
    print()
    print(f"   >>> MAXIMUM STATE-TO-STATE VARIATION IN THE STRUCTURE COEFFICIENTS, ALL SAMPLES: "
          f"{maxvar_all}  <<<")
    print(f"   >>> THE COEFFICIENTS ARE STRUCTURE CONSTANTS, NOT STRUCTURE FUNCTIONS. <<<")

    return dict(nx=nx, ny=ny, NV=NV, L=L, D=D, Dexp=Dexp, selfcheck=selfcheck,
                maxent=maxent, maxfro=maxfro, maxspec=maxspec, allinv=allinv,
                prod_is_I=prod_isI, rank=rank, gorder=gorder, k=k, nrel=len(kernel),
                varA=varA, varB=varB, varC=varC, worstop=worstop, maxvar=maxvar_all,
                coset=cosetA)


# ---------------------------------------------------------------------------------------------
rule("X6-a  —  WHAT ALGEBRA DO THE GAUSS CONSTRAINTS CLOSE INTO?")
print("Z2 lattice gauge theory, G_v = prod_{k incident v} sigma^x_k on the FULL link Hilbert space.")
print("All zeros below are EXACT integer zeros (max|entry| == 0), not small floats.")

rule("MACHINERY CONTROLS  —  checks that are able to FAIL")
c1 = positive_controls(8)
print(f"   CONTROLS: {'PASS' if c1 else 'FAIL'}")
if not c1:
    raise SystemExit("CONTROL FAILED — no headline number will be reported.")

r22 = run(2, 2, do_dense_spectral=True, n_states=48)
r23 = run(2, 3, do_dense_spectral=False, n_states=48)

rule("SUMMARY — BOTH LATTICES")
print(f"   {'lattice':>10} {'NV':>4} {'L':>4} {'D':>6} {'2^(L-NV+1)':>12} {'selfchk':>9} "
      f"{'max||[G,G]||':>14} {'G^2=I':>7} {'prod=I':>8} {'rank':>6} {'|grp|':>7} {'k':>4} "
      f"{'maxvar':>8} {'coset':>7}")
for r in (r22, r23):
    print(f"   {'%dx%d'%(r['nx'],r['ny']):>10} {r['NV']:>4} {r['L']:>4} {r['D']:>6} {r['Dexp']:>12} "
          f"{'PASS' if r['selfcheck'] else 'FAIL':>9} {r['maxfro']:>14.6f} "
          f"{str(r['allinv']):>7} {str(r['prod_is_I']):>8} {r['rank']:>6} {r['gorder']:>7} "
          f"{r['k']:>4} {r['maxvar']:>8} {str(r['coset']):>7}")
print()
print("   VERDICT (algebraic structure only; scale-free, no coupling anywhere in this lane):")
print("   The Gauss constraints close into an ELEMENTARY ABELIAN 2-GROUP  Z2^(NV-1):")
print("     * all generators commute exactly            -> the algebra is ABELIAN")
print("     * every generator squares to the identity   -> exponent 2")
print("     * prod over ALL vertices = identity         -> exactly one global relation, rank NV-1")
print("     * composition coefficients are integers in F2 that do not move from state to state")
print("       -> STRUCTURE CONSTANTS. Measured state-to-state variation: exactly 0.")
print()
print("   NOT ESTABLISHED HERE: anything about the ADM/hypersurface-deformation algebra itself.")
print("   No gravitational constraint algebra was computed in this script. The comparison that")
print("   the lane's falsifier turns on is made elsewhere; this file supplies only the Z2 side.")

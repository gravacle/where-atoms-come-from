# pin_carrier.py -- LANE W20_PRE.  PRE-REGISTRATION PHASE.  NO STATE IS BUILT ANYWHERE IN THIS FILE.
#
# Everything here is INCIDENCE COMBINATORICS on tri_chain12 plus finite-dimensional OPERATOR ALGEBRA.
# No Hamiltonian is diagonalised, no vector of length 2^n is ever allocated, no ground state, no
# Haar state, no coupling is evaluated.  The point of the phase is to fix every convention and to
# publish the vacuity map BEFORE any state can embarrass a lane.
#
# CONVENTIONS INHERITED VERBATIM FROM W-19 (W19_RULING/rule_verify.py, header):
#   Z_2 on each link.  Z|j> = (-1)^j |j>,  X|j> = |j+1 mod 2>.
#   G_v = prod_{l incident to v} X_l.   PHYSICAL sector: G_v = eta_v, eta_v in {+1,-1}.
#   W_z = prod_{l in z} Z_l for z in the GF(2) cycle space.
#   H = -(1/g2) sum_{p in PLAQ} W_p - g2 sum_l X_l.
#
# GF(2) LINEAR ALGEBRA is done on 12-bit masks over the edge set; vertex subsets on 8-bit masks.
import itertools, sys
from collections import deque

LOG = []
def P(*a):
    s = " ".join(str(x) for x in a); LOG.append(s); print(s, flush=True)

def rule(t=""):
    P("\n" + "=" * 100); P(t); P("=" * 100)

# ---------------------------------------------------------------- GF(2) helpers
def rank2(vecs):
    b = []
    for v in vecs:
        for x in b: v = min(v, v ^ x)
        if v: b = sorted(b + [v], reverse=True)
    return len(b), b

def span2(basis):
    out = {0}
    for b in basis:
        out |= {x ^ b for x in out}
    return sorted(out)

def bits(m, n):
    return [i for i in range(n) if m >> i & 1]

def pop(m): return bin(m).count("1")

# ---------------------------------------------------------------- THE CARRIER
def tri_chain12():
    return 8, [(0, 7), (0, 1), (0, 2), (1, 2), (1, 3), (2, 3),
               (3, 4), (4, 5), (4, 6), (5, 6), (5, 7), (6, 7)]

V, E = tri_chain12()
L = len(E)

rule("BLOCK 1 -- THE CARRIER, PUBLISHED")
P("carrier name        : tri_chain12   (inherited from W-19 sec 2a / sec 8.1; edge list byte-identical)")
P("V = %d   L = %d" % (V, L))
P("EDGE LIST, link index first (this is the published carrier):")
for i, (a, b) in enumerate(E):
    P("   link %2d = (%d,%d)" % (i, a, b))
deg = [0] * V
for (a, b) in E: deg[a] += 1; deg[b] += 1
P("degrees             : %s   -> cubic: %s" % (deg, all(d == 3 for d in deg)))
P("simple              : %s" % (len(set(tuple(sorted(e)) for e in E)) == L and all(a != b for a, b in E)))

adj = [[] for _ in range(V)]
for i, (a, b) in enumerate(E):
    adj[a].append((b, i)); adj[b].append((a, i))

def connected(edgeset_mask, verts):
    vs = [v for v in range(V) if verts >> v & 1]
    if not vs: return True
    seen = {vs[0]}; dq = deque([vs[0]])
    while dq:
        x = dq.popleft()
        for (y, i) in adj[x]:
            if (edgeset_mask >> i & 1) and (verts >> y & 1) and y not in seen:
                seen.add(y); dq.append(y)
    return len(seen) == len(vs)

P("connected           : %s" % connected((1 << L) - 1, (1 << V) - 1))

# ---------------------------------------------------------------- INCIDENCE: the two spaces
rule("BLOCK 2 -- THE INCIDENCE, AND THE DERIVATION OF 'SURFACE'")
P("The carrier is a 1-complex.  Its ONLY given structure is the incidence map")
P("     d : C_1 -> C_0 ,  d(edge) = its two endpoints        (over GF(2))")
P("and its transpose, the coboundary")
P("     delta : C^0 -> C^1 ,  delta(v) = star(v) = the edges at v .")
P("Two subspaces of C_1 are therefore FORCED, and nothing else is:")
P("   Z_1 = ker d      the CYCLE space      (Z-strings that are gauge invariant)")
P("   B^1 = im delta   the COCYCLE / CUT space (X-strings that are FIXED TO A SCALAR by Gauss)")

star = [0] * V
for i, (a, b) in enumerate(E):
    star[a] |= 1 << i; star[b] |= 1 << i
cocyc_basis_r, cocyc_basis = rank2(star)
COCYC = span2(cocyc_basis)

def dpartial(mask):
    r = 0
    for i in bits(mask, L):
        a, b = E[i]; r ^= (1 << a) ^ (1 << b)
    return r

CYC = [m for m in range(1 << L) if dpartial(m) == 0]
cyc_rank, cyc_basis = rank2(CYC)

P("")
P("dim Z_1 (cycle space)   = L - V + 1 = %d      |Z_1| = %d" % (cyc_rank, len(CYC)))
P("dim B^1 (cocycle space) = V - 1     = %d      |B^1| = %d   (nonzero: %d)"
  % (cocyc_basis_r, len(COCYC), len(COCYC) - 1))
P("dim Z_1 + dim B^1 = %d = L  : each is the ANNIHILATOR of the other under the GF(2)"
  % (cyc_rank + cocyc_basis_r))
P("intersection-form pairing <x,y> = |x AND y| mod 2 .  Verified exhaustively below.")
bad = [(z, c) for z in CYC for c in COCYC if pop(z & c) % 2]
P("pairs (z in Z_1, c in B^1) with |z AND c| ODD : %d   of %d checked"
  % (len(bad), len(CYC) * len(COCYC)))
bicyc = sorted(set(CYC) & set(COCYC))
P("BUT THEY ARE NOT COMPLEMENTARY AS SUBSPACES.  Over GF(2) a space can meet its own annihilator.")
P("Z_1 cap B^1 (the BICYCLE space) has dimension %d, %d elements:" % (rank2(bicyc)[0], len(bicyc)))
for x in bicyc:
    if x: P("   bicycle: links %s (weight %d) -- simultaneously a Wilson loop and a Gauss surface"
            % (bits(x, L), pop(x)))
P("CORRECTION TO THE OBVIOUS READING: 'the magnetic and electric surfaces are disjoint objects' is")
P("FALSE on this carrier.  %d supports are both.  What IS true is the commutation statement above:"
  % (len(bicyc) - 1))
P("X^{delta A} and Z^z always commute, because the PAIRING vanishes even where the SUPPORTS coincide.")
P("")
P(">>> DERIVED READING E  (ELECTRIC SURFACE).  A surface is a nonzero element of B^1, i.e. an EDGE")
P("    CUT delta(A) indexed by a VERTEX SET A.  It is forced by the Gauss law as an OPERATOR")
P("    IDENTITY:  prod_{v in A} G_v = X^{delta(A)} = prod_{v in A} eta_v , a SCALAR.")
P("    Count of distinct nonzero electric surfaces on this carrier: %d" % (len(COCYC) - 1))
P("    (delta has kernel {empty, V}, so A and its complement name the same surface: 2^8/2 - 1 = 127.)")
P(">>> DERIVED READING M  (MAGNETIC SURFACE).  A surface is a nonzero element of Z_1, i.e. a CYCLE.")
P("    It is forced by gauge invariance of Z-strings: Z^z commutes with every G_v iff z in ker d.")
P("    Count of distinct nonzero magnetic surfaces: %d" % (len(CYC) - 1))
P(">>> BOTH READINGS SURVIVE DERIVATION AND BOTH ARE CARRIED, as the brief instructs.  They are")
P("    MUTUAL ANNIHILATORS (not complements -- see the bicycle space), hence X^{delta A} and Z^z")
P("    ALWAYS COMMUTE (0 odd pairs above), while their SUPPORTS may coincide.")
P(">>> READING P (PLAQUETTE SET) DOES NOT SURVIVE.  A graph has no 2-cells.  A plaquette set is a")
P("    CHOICE of cycle basis, not an incidence datum -- see BLOCK 7, where the choice is 4-fold")
P("    degenerate on this very carrier.  It is COINED, and it is pinned in BLOCK 7 rather than")
P("    derived here.")

# ---------------------------------------------------------------- THE SYSTEM REGION
rule("BLOCK 3 -- THE SYSTEM REGION S, AND THE COLLAPSE OF READING R (REGION BOUNDARY)")
S_LINKS = [1, 2, 3]
S_MASK = sum(1 << i for i in S_LINKS)
A_SET = {0, 1, 2}
A_MASK = sum(1 << v for v in A_SET)
P("S (declared, inherited from W-19 sec 8.1) = links %s = %s"
  % (S_LINKS, [E[i] for i in S_LINKS]))
P("S is a cycle (triangle) on vertices {0,1,2}:  d(S) = %d  (0 means closed)" % dpartial(S_MASK))
induced = [i for i, (a, b) in enumerate(E) if a in A_SET and b in A_SET]
P("A = %s .  Links with BOTH endpoints in A: %s" % (sorted(A_SET), induced))
P("=> S IS VERTEX-INDUCED:  S = E[A].  This is what makes the third reading collapse.")
SIGMA_MASK = star[0] ^ star[1] ^ star[2]
SIGMA = bits(SIGMA_MASK, L)
P("delta(A) = star(0) XOR star(1) XOR star(2) = links %s = %s" % (SIGMA, [E[i] for i in SIGMA]))
cross = [i for i, (a, b) in enumerate(E) if (a in A_SET) != (b in A_SET)]
P("links with EXACTLY ONE endpoint in A (the 'region boundary' reading R): %s" % cross)
P("=> READING R == READING E EXACTLY, because S is vertex-induced.  Reading R is not independent")
P("   and is not carried separately.  SURFACE OF S := Sigma = delta({0,1,2}) = %s , |Sigma| = %d."
  % (SIGMA, len(SIGMA)))
P("boundary vertices of S (incident to both S and E\\S): %s"
  % sorted({v for i in cross for v in E[i] if v in A_SET}))

ENV = [i for i in range(L) if i not in S_LINKS]
ENV_MASK = sum(1 << i for i in ENV)
P("ENVIRONMENT E_env = %s   (|E_env| = %d)" % (ENV, len(ENV)))

# ---------------------------------------------------------------- GAUSS RELATIONS AT S
rule("BLOCK 4 -- THE GAUSS RELATIONS AT S.  THIS IS THE VACUITY MAP'S FOUNDATION.")
for v in sorted(A_SET):
    P("G_%d = X^{star(%d)} = X on links %s  =  eta_%d" % (v, v, bits(star[v], L), v))
P("product over A: X on links %s = eta_0*eta_1*eta_2   <-- the surface identity"
  % bits(SIGMA_MASK, L))
P("")
P("Rewritten so the forcing is visible (each vertex of A has star = 1 surface link + 2 S links):")
for v in sorted(A_SET):
    ins = sorted(set(bits(star[v], L)) & set(S_LINKS)); outs = sorted(set(bits(star[v], L)) - set(S_LINKS))
    P("   X_%d X_%d  =  eta_%d * X_%d        (S-pair %s forced by the single surface link %s)"
      % (ins[0], ins[1], v, outs[0], ins, outs))
P("")
P("=> THE THREE PAIRWISE X-PRODUCTS INSIDE S ARE EACH FORCED BY ONE SURFACE LINK.")
P("   Only 2 of the 3 are independent (their product is the surface identity), so the Gauss law")
P("   pins EXACTLY 2 BITS of S's electric algebra, and pins them TO THE SURFACE.")
P("   What is NOT forced: any single X_i on S, and the whole Z sector.")

# ---------------------------------------------------------------- THE FRAGMENTS
rule("BLOCK 5 -- THE NAIVE PARTITION, AND THE MEASUREMENT THAT REFUTES IT BEFORE IT IS USED")
P("This block evaluates W-19 sec 8.1's suggested partition (sizes 2,2,2,3, chosen BY SIZE) and")
P("the naive stratification rule 'a fragment with no surface link is Gauss-free'.  BOTH FAIL.")
P("The partition actually adopted is declared in pin_arms.py BLOCK 3.")
FRAGS = {"F1": [0, 4], "F2": [5, 6], "F3": [7, 8, 9], "F4": [10, 11]}
order = ["F1", "F2", "F3", "F4"]
allf = sorted(sum(FRAGS.values(), []))
P("F1 = %s = %s" % (FRAGS["F1"], [E[i] for i in FRAGS["F1"]]))
P("F2 = %s = %s" % (FRAGS["F2"], [E[i] for i in FRAGS["F2"]]))
P("F3 = %s = %s" % (FRAGS["F3"], [E[i] for i in FRAGS["F3"]]))
P("F4 = %s = %s" % (FRAGS["F4"], [E[i] for i in FRAGS["F4"]]))
P("pairwise disjoint   : %s" % (len(allf) == len(set(allf))))
P("union == E_env      : %s   (sizes %s, total %d)"
  % (allf == ENV, [len(FRAGS[k]) for k in order], len(allf)))
P("size cap |F| <= floor(|E_env|/2) = %d : max declared size = %d  -> satisfied"
  % (len(ENV) // 2, max(len(FRAGS[k]) for k in order)))

def forced_bits(fmask):
    """GF(2) dim of { delta(A') restricted to S : delta(A') subset of S union F }."""
    got = []
    for am in range(1 << V):
        c = 0
        for v in bits(am, V): c ^= star[v]
        if c & ~(S_MASK | fmask): continue
        got.append(c & S_MASK)
    return rank2(got)[0], sorted(set(got))

def cyc_in(mask):
    return rank2([z for z in CYC if not (z & ~mask)])[0]

P("")
P("PER-FRAGMENT STRATIFICATION  --  computed from the incidence alone, before any state:")
P("%-4s %-12s %-8s %-14s %-16s %-10s" % ("frag", "links", "|F cap Sigma|", "electric bits", "electric bits", "magnetic"))
P("%-4s %-12s %-8s %-14s %-16s %-10s" % ("", "", "", "of S FORCED", "of S FREE (of 3)", "bits in F"))
P("-" * 84)
for k in order:
    fm = sum(1 << i for i in FRAGS[k])
    nf, _ = forced_bits(fm)
    P("%-4s %-12s %-8d %-14d %-16d %-10d"
      % (k, str(FRAGS[k]), pop(fm & SIGMA_MASK), nf, 3 - nf, cyc_in(fm)))
P("-" * 84)
allm = sum(1 << i for i in ENV)
P("%-4s %-12s %-8d %-14d %-16d %-10d" % ("ALL", "E_env", pop(allm & SIGMA_MASK),
                                          forced_bits(allm)[0], 3 - forced_bits(allm)[0], cyc_in(allm)))
P("")
P(">>> READ COLUMN 3 AGAINST COLUMN 4.  THE NAIVE RULE IS REFUTED BY ITS OWN TABLE:")
P("    F3 = [7,8,9] and F4 = [10,11] contain ZERO surface links and STILL each force one bit of")
P("    S's electric algebra.  The witnesses, printed rather than asserted:")
for k in ["F3", "F4"]:
    fm = sum(1 << i for i in FRAGS[k])
    for am in range(1 << V):
        c = 0
        for v in bits(am, V): c ^= star[v]
        if (c & ~(S_MASK | fm)) or not (c & S_MASK): continue
        if pop(am) <= V // 2:
            P("       %s : delta(%s) = links %s  =>  X on %s = eta * X on %s"
              % (k, bits(am, V), bits(c, L), bits(c & S_MASK, L), bits(c & ~S_MASK, L)))
            break
P("    Sigma = [0,4,5] is NOT the only Gauss surface enclosing S's electric pairs.  delta({0,7}) =")
P("    [1,2,10,11] and delta({0,5,6,7}) = [1,2,7,8] enclose the SAME pair X_1X_2 while avoiding")
P("    Sigma entirely.  THE SURFACE THAT FORCES A GIVEN RECORD BIT IS A COSET, NOT A SET.")
P("")
P(">>> CONSEQUENCE FOR THE DESIGN, AND IT IS WHY THIS PHASE EXISTS: the vacuity of a fragment is")
P("    NOT |F cap Sigma|.  It is D(F), a GF(2) rank over all 2^8 vertex subsets.  A lane that")
P("    stratified by surface-link count would have called F3 and F4 'live', measured a plateau")
P("    that the Gauss law had already guaranteed one bit of, and reported it.")
P(">>> W-19's SUGGESTED PARTITION HAS D = 2,2,1,1 : NO FRAGMENT IS FREE.  IT IS NOT ADOPTED.")

# ---------------------------------------------------------------- MAGNETIC CONTENT
rule("BLOCK 6 -- MAGNETIC CONTENT: WHAT CAN CARRY A WILSON LOOP, AND WHAT CANNOT")
tris = []
for c in CYC:
    if pop(c) == 3: tris.append(c)
P("all triangles (weight-3 cycles): %d" % len(tris))
for t in tris:
    vs = sorted({v for i in bits(t, L) for v in E[i]})
    P("   links %s  on vertices %s" % (bits(t, L), vs))
P("")
P("cycle space restricted to the ENVIRONMENT E_env = %s : dimension %d" % (ENV, cyc_in(ENV_MASK)))
envc = sorted([z for z in CYC if z and not (z & ~ENV_MASK)], key=pop)
for z in envc:
    P("   env cycle: links %s (weight %d)" % (bits(z, L), pop(z)))
P("=> the three environment cycles PAIRWISE INTERSECT (all contain link 9 or overlap on 7,8):")
for i in range(len(envc)):
    for j in range(i + 1, len(envc)):
        P("   |%s AND %s| = %d" % (bits(envc[i], L), bits(envc[j], L), pop(envc[i] & envc[j])))
P(">>> DERIVED, AND IT FORCES THE DESIGN: AT MOST ONE PAIRWISE-DISJOINT ENVIRONMENT FRAGMENT CAN")
P("    CARRY A WILSON LOOP.  There is no partition of E_env into disjoint fragments with two")
P("    magnetic fragments.  F3 = [7,8,9] is that fragment, and it is the SHAPE TWIN of S:")
P("    S = triangle on {0,1,2}, F3 = triangle on {4,5,6}, both weight 3, disjoint, and F3 shares")
P("    NO surface link with S.")
P("")
P("W_S := Z_1 Z_2 Z_3, the Wilson loop of S.  Is it determined by any environment operator?")
P("   Gauss forces X-strings on B^1 to scalars.  Z_1 (cycles) is the ORTHOGONAL COMPLEMENT of B^1")
P("   (BLOCK 2: zero odd pairs).  NOTHING in the Gauss law constrains any Z-string.")
P("   env cycle space = span%s ;  is S's cycle %s a member?  %s"
  % ([bits(z, L) for z in envc], bits(S_MASK, L), S_MASK in span2(rank2(envc)[1])))
P("   => NO.  W_S is not equal to any environment Wilson loop, and no Gauss constraint touches it.")
P(">>> THE MAGNETIC CHANNEL ON S CANNOT BE GAUSS-FORCED.  That is the derived reason the falsifier")
P("    is LIVE there and VACUOUS in the electric channel on F1/F2.")

# ---------------------------------------------------------------- PLAQUETTES
rule("BLOCK 7 -- THE PLAQUETTE SET IS COINED, IT IS 4-FOLD DEGENERATE, AND IT IS PINNED HERE")
byw = {}
for c in CYC:
    if c: byw.setdefault(pop(c), []).append(c)
P("cycle weight spectrum: %s" % {w: len(byw[w]) for w in sorted(byw)})
# matroid greedy, but enumerate ALL minimum-weight cycle bases to expose the tie
def greedy_all():
    order_c = sorted([c for c in CYC if c], key=lambda c: (pop(c), c))
    red = []; chosen = []
    for c in order_c:
        v = c
        for b in red: v = min(v, v ^ b)
        if v:
            red = sorted(red + [v], reverse=True); chosen.append(c)
            if len(chosen) == cyc_rank: break
    return chosen
G = greedy_all()
P("matroid-greedy min-weight cycle basis (weight, then edge-index tuple, ascending):")
for c in G: P("   plaquette links %s (weight %d)" % (bits(c, L), pop(c)))
P("total basis weight = %d" % sum(pop(c) for c in G))
# count all min-weight bases
minw = sum(pop(c) for c in G)
cnt = 0; examples = []
for comb in itertools.combinations([c for c in CYC if c], cyc_rank):
    if sum(pop(c) for c in comb) != minw: continue
    if rank2(list(comb))[0] == cyc_rank:
        cnt += 1
        if len(examples) < 8: examples.append(comb)
P("number of DISTINCT minimum-weight cycle bases (total weight %d): %d" % (minw, cnt))
P("they differ ONLY in the 5th generator; the four triangles are in every one of them:")
fifth = sorted({tuple(sorted(bits(c, L))) for comb in examples for c in comb if pop(c) > 3})
for f in fifth: P("   candidate 5th plaquette: links %s (weight %d)" % (list(f), len(f)))
P("")
P(">>> PINNED, BY AN EXPLICIT RULE STATED BEFORE ANY STATE:")
P("    PLAQ = the 4 triangles, plus the lexicographically smallest edge-index tuple among the")
P("    minimum-weight 5th generators.  That is:")
for c in G: P("       %s" % bits(c, L))
P("    THE TIE IS %d-FOLD AND IT IS RECORDED, NOT HIDDEN.  Changing sum_p W_p changes H, so the" % cnt)
P("    plaquette set is NOT a convention -- it is part of the dynamics.  ARM D2 varies it.")

# ---------------------------------------------------------------- SURFACES, ENUMERATED
rule("BLOCK 8 -- THE 127 ELECTRIC SURFACES, AND WHERE Sigma SITS AMONG THEM")
sizes = {}
for am in range(1, 1 << V):
    if am >> 0 & 1 and am == (1 << V) - 1: pass
    c = 0
    for v in bits(am, V): c ^= star[v]
    if c: sizes.setdefault(pop(c), set()).add(c)
P("distinct nonzero cocycles by size |delta(A)|:")
for w in sorted(sizes): P("   |Sigma| = %2d : %3d surfaces" % (w, len(sizes[w])))
P("total distinct nonzero surfaces = %d" % sum(len(s) for s in sizes.values()))
P("minimum surface size = %d .  NOT a vertex star (the carrier is cubic, every star has 3 links):"
  % min(sizes))
for c in sizes[min(sizes)]:
    P("   the unique 2-link surface is %s -- tri_chain12 is 2-EDGE-CONNECTED, NOT 3-edge-connected."
      % bits(c, L))
    P("   Links 0 and 6 are a bottleneck separating {0,1,2,3} from {4,5,6,7}, and link 0 is IN Sigma.")
    P("   CONFOUND, RECORDED NOT FIXED: the declared surface shares a link with the carrier's")
    P("   global bottleneck, so 'surface' and 'bottleneck' are not independent on this carrier.")
P("Sigma = delta({0,1,2}) has |Sigma| = %d -- it is a MINIMUM-SIZE surface, and it is a GENUINE cut"
  % pop(SIGMA_MASK))
P("   (both sides connected: A side %s, complement side %s)"
  % (connected(S_MASK, A_MASK), connected(ENV_MASK & ~SIGMA_MASK, ((1 << V) - 1) ^ A_MASK)))
P("surfaces of size 3 whose enclosed region induces a TRIANGLE (i.e. candidate system regions):")
for am in range(1, 1 << V):
    c = 0
    for v in bits(am, V): c ^= star[v]
    if pop(c) != 3: continue
    ind = [i for i, (a, b) in enumerate(E) if (am >> a & 1) and (am >> b & 1)]
    if len(ind) == 3 and dpartial(sum(1 << i for i in ind)) == 0:
        P("   A = %s  ->  S = links %s ,  Sigma = links %s" % (bits(am, V), ind, bits(c, L)))

# ---------------------------------------------------------------- CHARGE SECTORS
rule("BLOCK 9 -- THE CHARGE SECTOR: WHAT IS ADMISSIBLE, AND THE FORMATION-ROUTE PAIRS")
P("prod_{v in V} G_v = X^{delta(V)} = X^0 = I  =>  prod_v eta_v = +1 : an EVEN number of charges.")
n_sectors = 2 ** (V - 1)
P("admissible charge sectors: 2^(V-1) = %d" % n_sectors)
P("The electric flux through Sigma is FIXED BY THE CHARGE SECTOR ALONE:")
P("     X^Sigma = eta_0 * eta_1 * eta_2  =  the enclosed charge parity of A = {0,1,2}.")
P("")
P("FORMATION-ROUTE PAIRS.  'Same surface, same flux through it, different formation.'")
P("Sigma is held byte-identical in all of them; only WHERE the charge sits moves.")
routes = [
    ("C0  vacuum",              [],       "reference"),
    ("C1  charges at {0,4}",    [0, 4],   "one charge INSIDE A, one outside -> flux(Sigma) = -1"),
    ("C2  charges at {1,5}",    [1, 5],   "one charge INSIDE A, one outside -> flux(Sigma) = -1"),
    ("C3  charges at {4,5}",    [4, 5],   "BOTH outside A                   -> flux(Sigma) = +1"),
    ("C4  charges at {0,1}",    [0, 1],   "BOTH inside A                    -> flux(Sigma) = +1"),
]
P("%-22s %-10s %-12s %s" % ("route", "eta = -1 at", "flux(Sigma)", "note"))
for nm, ch, note in routes:
    fl = 1
    for v in ch:
        if v in A_SET: fl = -fl
    P("%-22s %-10s %-12s %s" % (nm, str(ch), "%+d" % fl, note))
P("")
P(">>> THE TWO FORMATION-ROUTE CONTRASTS, BOTH WITH THE SURFACE AND ITS FLUX HELD FIXED:")
P("    PAIR ALPHA : C1 vs C2 .  Same Sigma, flux(Sigma) = -1 in both, charge placed at a DIFFERENT")
P("                 vertex of A (0 vs 1) and a different vertex outside (4 vs 5).")
P("    PAIR BETA  : C0 vs C3 .  Same Sigma, flux(Sigma) = +1 in both, but C3 has a real charge pair")
P("                 sitting entirely OUTSIDE A.  The surface and its flux are identical; the")
P("                 formation of the configuration is not.")
P("    PAIR GAMMA : C0 vs C4 .  Same Sigma, flux(Sigma) = +1 in both, charge pair entirely INSIDE A.")
P("  If the record content on S is the same across a pair, FORMATION is NARRATION on that pair and")
P("  only the surface + its flux matter.  If it differs, FORMATION is CAUSAL_EARNED and the number")
P("  is the difference.  NEITHER OUTCOME IS SCORED AS CONFIRMATION OF THE PREFERRED HYPOTHESIS.")

# ---------------------------------------------------------------- COMPOSITES
rule("BLOCK 10 -- JOIN vs UNION: THE GAP, MEASURED FROM THE INCIDENCE")
for a, b in [("F3", "F4"), ("F1", "F2"), ("F2", "F3")]:
    ma = sum(1 << i for i in FRAGS[a]); mb = sum(1 << i for i in FRAGS[b])
    P("%s cup %s : magnetic bits  %d (in %s) + %d (in %s) = %d   vs   %d in the union"
      % (a, b, cyc_in(ma), a, cyc_in(mb), b, cyc_in(ma) + cyc_in(mb), cyc_in(ma | mb)))
P(">>> F3 cup F4 gains ONE FULL BIT that neither fragment has: the Wilson loop [9,10,11] straddles")
P("    them.  So JOIN != UNION-OF-ALGEBRAS by exactly 1 bit on this partition.  AXIS 5 is pinned")
P("    by NOT USING COMPOSITES IN R_delta at all; any composite reported is reported separately")
P("    with the join algebra named and this 1-bit gap quoted.")

# ---------------------------------------------------------------- WRITE
rule("BLOCK 11 -- WHAT THIS FILE DID NOT DO")
P("No Hilbert space was constructed.  No Hamiltonian was diagonalised.  No coupling was evaluated.")
P("No state, ground or random, exists anywhere in this lane.  Every number above is incidence")
P("combinatorics on an 8-vertex 12-edge graph, exhaustively enumerated (2^12 = 4096 edge masks,")
P("2^8 = 256 vertex masks).  Nothing here can be embarrassed by a state because nothing here")
P("depends on one -- which is the entire point of running this phase first.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_W20_PRE/OUT_pin_carrier.txt", "w").write("\n".join(LOG) + "\n")

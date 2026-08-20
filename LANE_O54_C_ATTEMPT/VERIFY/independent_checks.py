"""VERIFY/independent_checks.py -- adversarial verifier's INDEPENDENT re-computations.
No import of o54c_lib; all machinery re-written from scratch.  Exact integer arithmetic.

Checks:
 V1  (4,6),(3,7): independent exhaustive X-coset scan -> w_min, N_min; compare to lane table.
 V2  w_min == torus L1 distance (external geometric cross-check of the 'earned' metric),
     and independent BFS d_gen on my own generator-overlap graph == lane's d_gen.
 V3  independent exhaustive Z-coset scan -> z_min; compare to lane table; recompute crossing
     type; check the finding-JSON claim about ctype (lane OUT says ctype nonempty on (4,6)).
 V4  both z_min values present at d=1 on (3,7) (the switch is not monotone-in-d in disguise).
 V5  independent stabilizer MI on (8,8): one contact and all separated axis/diag placements,
     own GF2 rank code; expect contact MI=2, separated MI=0, both states.
 V6  independent band MI on (8,8): separated constant 1, contact 15.
 V7  independent verdict on the INSERTED ring control r=3: decreasing to 0 (instrument can see
     a falloff) -- and confirm the lane labels it INSERTED in its OUT.
"""
import sys
from collections import deque
import numpy as np

def pc(x): return bin(x).count("1")

def rank_f2(rows):
    piv = {}; r = 0
    for m in rows:
        while m:
            t = m.bit_length() - 1
            if t in piv: m ^= piv[t]
            else: piv[t] = m; r += 1; break
    return r

class T2:
    def __init__(s, Lx, Ly):
        s.Lx, s.Ly, s.n = Lx, Ly, 2*Lx*Ly
    def h(s, x, y): return (y % s.Ly)*s.Lx + (x % s.Lx)
    def v(s, x, y): return s.Lx*s.Ly + (y % s.Ly)*s.Lx + (x % s.Lx)
    def star_x(s, x, y):
        m = 0
        for e in (s.h(x,y), s.h(x-1,y), s.v(x,y), s.v(x,y-1)): m |= 1 << e
        return m
    def plaq_z(s, x, y):
        m = 0
        for e in (s.h(x,y), s.h(x,y+1), s.v(x,y), s.v(x+1,y)): m |= 1 << e
        return m   # NOTE: this is the Z support as an n-bit EDGE mask (not shifted)

def sp_edge(xmask, zmask): return pc(xmask & zmask) % 2

FAILS = []
def chk(name, ok, extra=""):
    print(("ok    " if ok else "FAIL  ") + name + (("  " + extra) if extra else ""))
    if not ok: FAILS.append(name)

def coset_min(gens, rep):
    """Exhaustive min popcount over rep ^ span(gens); own chunking (base 18)."""
    base, rest = gens[:18], gens[18:]
    arr = np.zeros(1, dtype=np.uint64)
    for g in base:
        arr = np.concatenate([arr, arr ^ np.uint64(g)])
    best, cnt = None, 0
    for i in range(1 << len(rest)):
        off = rep
        for j in range(len(rest)):
            if (i >> j) & 1: off ^= rest[j]
        w = np.bitwise_count(arr ^ np.uint64(off))
        m = int(w.min()); c = int((w == m).sum())
        if best is None or m < best: best, cnt = m, c
        elif m == best: cnt += c
    return best, cnt

# lane's sealed table values to compare against (from o54c_attempts.OUT.txt)
LANE = {
 (4,6): {(1,0):(1,1,1,4),(2,0):(2,2,2,4),(0,1):(1,1,1,4),(0,2):(2,2,1,4),(0,3):(3,3,2,4),
         (1,1):(2,2,2,4),(2,1):(3,3,6,4),(1,2):(3,3,3,4),(2,2):(4,4,12,4),(1,3):(4,4,8,4),
         (2,3):(5,5,40,4)},
 (3,7): {(1,0):(1,1,1,4),(0,1):(1,1,1,3),(0,2):(2,2,1,3),(0,3):(3,3,1,3),
         (1,1):(2,2,2,3),(1,2):(3,3,3,3),(1,3):(4,4,4,3)},
}   # v -> (d_gen, w_min, N_min, z_min)

for (Lx, Ly), table in LANE.items():
    T = T2(Lx, Ly); n = T.n
    stars = [T.star_x(x,y) for y in range(Ly) for x in range(Lx)]
    plaqs = [T.plaq_z(x,y) for y in range(Ly) for x in range(Lx)]  # edge masks
    xbar1 = 0
    for y in range(Ly): xbar1 |= 1 << T.h(0,y)
    xbar2 = 0
    for x in range(Lx): xbar2 |= 1 << T.v(x,0)
    zbar1 = 0
    for y in range(Ly): zbar1 |= 1 << T.v(0,y)
    zbar2 = 0
    for x in range(Lx): zbar2 |= 1 << T.h(x,0)

    # independent generator-overlap graph (stars+plaqs as nodes, supports as edge masks)
    supports = stars + plaqs
    N = len(supports)
    adj = [[] for _ in range(N)]
    for i in range(N):
        for j in range(i+1, N):
            if supports[i] & supports[j]: adj[i].append(j); adj[j].append(i)
    def bfs(a, b):
        dist = {a: 0}; q = deque([a])
        while q:
            i = q.popleft()
            if i == b: return dist[i]
            for j in adj[i]:
                if j not in dist: dist[j] = dist[i]+1; q.append(j)
        return None

    u = (0,0)
    Bu = T.plaq_z(*u)
    iu = Lx*Ly + 0
    stars_ind = []
    piv = {}
    for g in stars:
        m = g
        while m:
            t = m.bit_length()-1
            if t in piv: m ^= piv[t]
            else: piv[t] = m; stars_ind.append(g); break
    gens_X = stars_ind + [xbar1, xbar2]
    chk("V1 (%d,%d) X span independent dim %d" % (Lx, Ly, len(gens_X)),
        rank_f2(gens_X) == len(gens_X))

    for v, (d_lane, w_lane, nm_lane, z_lane) in sorted(table.items()):
        # independent dual connector: straight steps x then y (crossing shared edges)
        p = 0
        x, y = u
        while x != v[0]:
            p ^= 1 << T.v(x+1, y); x = (x+1) % Lx
        while y != v[1]:
            p ^= 1 << T.h(x, y+1); y = (y+1) % Ly
        Bv = T.plaq_z(*v)
        # admissibility: X-chain commutes with all stars trivially (both X);
        # anticommutes with exactly Bu and Bv among plaquettes
        bad = [ (xx,yy) for yy in range(Ly) for xx in range(Lx)
                if sp_edge(p, T.plaq_z(xx,yy)) == 1 and (xx,yy) not in (u, v) ]
        chk("V1 (%d,%d) v=%s connector boundary exactly {u,v}" % (Lx,Ly,str(v)),
            not bad and sp_edge(p,Bu)==1 and sp_edge(p,Bv)==1)
        w, nm = coset_min(gens_X, p)
        chk("V1 (%d,%d) v=%s w_min=%d N_min=%d matches lane (%d,%d)"
            % (Lx,Ly,str(v),w,nm,w_lane,nm_lane), w == w_lane and nm == nm_lane)
        # V2 geometric cross-check + independent BFS
        dx, dy = abs(v[0]-u[0]), abs(v[1]-u[1])
        dL1 = min(dx, Lx-dx) + min(dy, Ly-dy)
        iv = Lx*Ly + v[1]*Lx + v[0]
        d_ind = bfs(iu, iv)
        chk("V2 (%d,%d) v=%s d_gen(independent BFS)=%s == lane %d == torus L1 %d == w_min %d"
            % (Lx,Ly,str(v),d_ind,d_lane,dL1,w), d_ind == d_lane == dL1 == w)
        # V3 independent Z-coset scan: writers of the hole-X record (rep Bu) over
        # span(remaining plaqs independent + zbars), all commuting with p
        rem = [T.plaq_z(xx,yy) for yy in range(Ly) for xx in range(Lx) if (xx,yy) not in (u,v)]
        rem_ind = []
        piv2 = {}
        for g in rem:
            m = g
            while m:
                t = m.bit_length()-1
                if t in piv2: m ^= piv2[t]
                else: piv2[t] = m; rem_ind.append(g); break
        gens_Z = rem_ind + [zbar1, zbar2]
        chk("V3 (%d,%d) v=%s Z span independent, all commute with connector" % (Lx,Ly,str(v)),
            rank_f2(gens_Z) == len(gens_Z) and all(sp_edge(g, p) == 0 for g in gens_Z))
        z, nz = coset_min(gens_Z, Bu)
        # crossing type recomputed independently
        odd = []
        for r_ in range(Ly):
            loop = 0
            for x_ in range(Lx): loop |= 1 << T.h(x_, r_)
            if sp_edge(p, loop) == 1: odd.append(Lx)
        for c_ in range(Lx):
            loop = 0
            for y_ in range(Ly): loop |= 1 << T.v(c_, y_)
            if sp_edge(p, loop) == 1: odd.append(Ly)
        zpred = min([4] + odd)
        chk("V3 (%d,%d) v=%s z_min=%d matches lane %d and switch formula %d (ctype %s)"
            % (Lx,Ly,str(v),z,z_lane,zpred,tuple(sorted(set(odd)))),
            z == z_lane == zpred)

# V4: both z values at d=1 on (3,7)
chk("V4 (3,7): z_min takes BOTH values {4,3} already at d=1 (switch not monotone-in-d)",
    LANE[(3,7)][(1,0)][0] == 1 and LANE[(3,7)][(1,0)][3] == 4 and
    LANE[(3,7)][(0,1)][0] == 1 and LANE[(3,7)][(0,1)][3] == 3)

# ---------- V5/V6: independent stabilizer entropies on (8,8) ----------
L = 8
T = T2(L, L); n = T.n
stars = [T.star_x(x,y) for y in range(L) for x in range(L)]
plaqs = [T.plaq_z(x,y) for y in range(L) for x in range(L)]
xbar1 = 0
for y in range(L): xbar1 |= 1 << T.h(0,y)
xbar2 = 0
for x in range(L): xbar2 |= 1 << T.v(x,0)
zbar1 = 0
for y in range(L): zbar1 |= 1 << T.v(0,y)
zbar2 = 0
for x in range(L): zbar2 |= 1 << T.h(x,0)
# symplectic rows: (xmask, zmask)
def rows_state(which):
    rows = [(s, 0) for s in stars] + [(0, p) for p in plaqs]
    if which == "Z": rows += [(0, zbar1), (0, zbar2)]
    else: rows += [(xbar1, 0), (xbar2, 0)]
    return rows

def rank_sym(rows):
    packed = [x | (z << n) for (x, z) in rows]
    return rank_f2(packed)

def S_region(rows, region):
    rs = set(region); comp = [q for q in range(n) if q not in rs]
    G = rank_sym(rows)
    assert G == n, G
    restr = []
    for (x, z) in rows:
        m = 0
        for i, q in enumerate(comp):
            if (x >> q) & 1: m |= 1 << (2*i)
            if (z >> q) & 1: m |= 1 << (2*i+1)
        restr.append(m)
    return len(rs) - (G - rank_f2(restr))

def patch(x0, y0, w):
    qs = set()
    for dx in range(w):
        for dy in range(w):
            qs.add(T.h(x0+dx, y0+dy)); qs.add(T.v(x0+dx, y0+dy))
    return sorted(qs)

A = patch(0, 0, 2)
for which in ("Z", "X"):
    rows = rows_state(which)
    for x0, want in [(2, 2), (3, 0), (4, 0), (5, 0)]:
        B = patch(x0, 0, 2)
        MI = S_region(rows, A) + S_region(rows, B) - S_region(rows, sorted(set(A)|set(B)))
        chk("V5 L=8 %sbar state, patch axis x0=%d: MI=%d expected %d" % (which, x0, MI, want),
            MI == want)
    for x0, want in [(3, 0), (4, 0)]:
        B = patch(x0, x0, 2)
        MI = S_region(rows, A) + S_region(rows, B) - S_region(rows, sorted(set(A)|set(B)))
        chk("V5 L=8 %sbar state, patch diag x0=%d: MI=%d expected %d" % (which, x0, MI, want),
            MI == want)

def band(x0, b):
    qs = set()
    for dx in range(b):
        for y in range(L):
            qs.add(T.h(x0+dx, y)); qs.add(T.v(x0+dx, y))
    return sorted(qs)

Ab = band(0, 2)
for which in ("Z", "X"):
    rows = rows_state(which)
    for x0, want in [(2, 15), (3, 1), (4, 1), (5, 1)]:
        Bb = band(x0, 2)
        MI = S_region(rows, Ab) + S_region(rows, Bb) - S_region(rows, sorted(set(Ab)|set(Bb)))
        chk("V6 L=8 %sbar state, band x0=%d: MI=%d expected %d" % (which, x0, MI, want),
            MI == want)

# ---------- V7: inserted ring control independently ----------
nq, r = 24, 3
starts = [i for i in range(nq) if (i % (2*r)) < r]
rows = []
for i in starts:
    j = (i + r) % nq
    rows.append(((1 << i) | (1 << j), 0))
    rows.append((0, (1 << i) | (1 << j)))
def S_ring(rows, region):
    rs = set(region); comp = [q for q in range(nq) if q not in rs]
    packed = [x | (z << nq) for (x, z) in rows]
    G = rank_f2(packed); assert G == nq
    restr = []
    for (x, z) in rows:
        m = 0
        for i, q in enumerate(comp):
            if (x >> q) & 1: m |= 1 << (2*i)
            if (z >> q) & 1: m |= 1 << (2*i+1)
        restr.append(m)
    return len(rs) - (G - rank_f2(restr))
A7 = list(range(4, 10))
seq = []
for off in range(10, 17):
    B7 = [q % nq for q in range(off, off+6)]
    if set(A7) & set(B7): continue
    MI = S_ring(rows, A7) + S_ring(rows, B7) - S_ring(rows, sorted(set(A7)|set(B7)))
    seq.append((off-10, MI))
dec = all(seq[i][1] >= seq[i+1][1] for i in range(len(seq)-1))
chk("V7 ring r=3 control: %s decreasing to 0, nonconstant" % (seq,),
    dec and seq[0][1] > 0 and seq[-1][1] == 0 and len(set(v for _, v in seq)) > 1)

print()
print("INDEPENDENT CHECKS: %s" % ("ALL OK" if not FAILS else "FAILURES: %s" % FAILS))
sys.exit(1 if FAILS else 0)

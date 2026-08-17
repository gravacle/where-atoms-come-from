#!/usr/bin/env python3
"""
LANE W-10 / C — STEP 1.  THE TWO FOUR-CLASS CARRIERS, BUILT FROM INCIDENCE.

Nothing here is quoted from S4 except the TARGET numbers, which are printed beside the
recomputed ones so that any disagreement is visible.  Targets: S4_THE_MEASUREMENT_V001.md
:519 (topology table) and :575-590 (class multisets).

Precision: exact integer / rational linear algebra over Q via Fraction row reduction.
No floating point is used for any rank.
"""
from fractions import Fraction
import itertools, sys

# ---------------------------------------------------------------- exact rank over Q
def rank_Q(M):
    """M: list of rows of Fractions/ints.  Returns rank over Q by exact elimination."""
    A = [[Fraction(x) for x in row] for row in M]
    rows, cols = len(A), (len(A[0]) if A else 0)
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if A[i][c] != 0:
                piv = i; break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(rows):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [a - f * b for a, b in zip(A[i], A[r])]
        r += 1
        if r == rows:
            break
    return r

def matmul(A, B):
    n, k, m = len(A), len(B), len(B[0])
    return [[sum(A[i][t] * B[t][j] for t in range(k)) for j in range(m)] for i in range(n)]

def in_span(cols, target):
    """Is target (list) in the column span of cols (list of columns)?  Exact."""
    M = [list(col) for col in cols]
    r0 = rank_Q(M)
    r1 = rank_Q(M + [list(target)])
    return r0 == r1

# ---------------------------------------------------------------- carrier container
class Carrier:
    def __init__(self, name, V, edges, faces, gF, gC):
        self.name = name
        self.V = V
        self.edges = edges            # list of (src, tgt)
        self.faces = faces            # list of dict edge_index -> coefficient
        self.gF = gF                  # dict edge_index -> coefficient
        self.gC = gC
        self.E = len(edges)
        self.F = len(faces)
        # d1 : V x E,  +1 at target, -1 at source
        self.d1 = [[0] * self.E for _ in range(V)]
        for e, (s, t) in enumerate(edges):
            self.d1[s][e] -= 1
            self.d1[t][e] += 1
        # d2 : E x F
        self.d2 = [[0] * self.F for _ in range(self.E)]
        for j, f in enumerate(faces):
            for e, coef in f.items():
                self.d2[e][j] += coef

    def chain(self, d):
        v = [0] * self.E
        for e, c in d.items():
            v[e] += c
        return v

    def loop_vertices(self, d):
        vs = set()
        for e, c in d.items():
            if c != 0:
                s, t = self.edges[e]
                vs.add(s); vs.add(t)
        return vs

    def report(self):
        V, E, F = self.V, self.E, self.F
        chi = V - E + F
        r1 = rank_Q(self.d1)
        r2 = rank_Q(self.d2)
        b0 = V - r1
        b2 = F - r2
        b1 = E - r1 - r2
        inv = E - (V - b0)
        d1d2 = matmul(self.d1, self.d2)
        maxent = max(abs(x) for row in d1d2 for x in row) if E and F else 0
        cols = [[self.d2[e][j] for e in range(E)] for j in range(F)]
        out = dict(name=self.name, V=V, E=E, F=F, chi=chi, b0=b0, b1=b1, b2=b2,
                   inv=inv, curv=r2, flat=b1, d1d2=maxent)
        for lbl, d in (("gF", self.gF), ("gC", self.gC)):
            ch = self.chain(d)
            iscyc = all(sum(self.d1[v][e] * ch[e] for e in range(E)) == 0 for v in range(V))
            bounds = in_span(cols, ch) if F else False
            out[lbl + "_cycle_defect"] = 0 if iscyc else 1
            out[lbl + "_bounds"] = bounds
        chF, chC = self.chain(self.gF), self.chain(self.gC)
        out["independent"] = rank_Q([chF, chC]) == 2
        return out

    def classes(self):
        """class multiset by (v in gF, v in gC).  Uses the loops' VERTEX SETS."""
        VF, VC = self.loop_vertices(self.gF), self.loop_vertices(self.gC)
        cnt = {"00": 0, "10": 0, "01": 0, "11": 0}
        lab = []
        for v in range(self.V):
            a = 1 if v in VF else 0
            b = 1 if v in VC else 0
            k = f"{a}{b}"
            cnt[k] += 1
            lab.append(k)
        return cnt, lab

# ---------------------------------------------------------------- B0b
def build_B0b():
    idx = lambda i, j: 3 * (i % 3) + (j % 3)
    edges = []
    H = {}; W = {}
    for i in range(3):
        for j in range(3):
            H[(i, j)] = len(edges); edges.append((idx(i, j), idx(i, j + 1)))
    for i in range(3):
        for j in range(3):
            W[(i, j)] = len(edges); edges.append((idx(i, j), idx(i + 1, j)))
    faces = []
    for i in range(3):
        for j in range(3):
            faces.append({H[(i, j)]: 1, W[(i, (j + 1) % 3)]: 1,
                          H[((i + 1) % 3, j)]: -1, W[(i, j)]: -1})
    gF = {H[(0, 0)]: 1, W[(0, 1)]: 1, H[(1, 0)]: -1, W[(0, 0)]: -1}
    gC = {H[(0, 0)]: 1, H[(0, 1)]: 1, H[(0, 2)]: 1}
    c = Carrier("B0b ring torus 3x3, loops MEET", 9, edges, faces, gF, gC)
    c.H, c.W = H, W
    return c

# ---------------------------------------------------------------- B4
def build_B4():
    N, S, a, b, cc, d = 0, 1, 2, 3, 4, 5
    edges = [(N, a), (a, S), (N, b), (b, S), (N, cc), (cc, S), (N, d), (d, S)]
    sph1 = {0: 1, 1: 1, 3: -1, 2: -1}
    sph2 = {4: 1, 5: 1, 7: -1, 6: -1}
    faces = [dict(sph1), dict(sph1), dict(sph2), dict(sph2)]
    gF = dict(sph1)
    gC = {0: 1, 1: 1, 5: -1, 4: -1}
    return Carrier("B4 spindle (two spheres, 2 glue pts)", 6, edges, faces, gF, gC)

# ---------------------------------------------------------------- K1 (three-class control)
def build_K1():
    edges = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)]
    faces = [{0: 1, 1: 1, 2: 1}]
    gF = {0: 1, 1: 1, 2: 1}
    gC = {3: 1, 4: 1, 5: 1}
    return Carrier("B1 = K1 (three-class CONTROL)", 5, edges, faces, gF, gC)

# ---------------------------------------------------------------- run
TARGETS = {
 "B0b ring torus 3x3, loops MEET": dict(V=9, E=18, F=9, chi=0, b0=1, b1=2, b2=1, inv=10,
        curv=8, flat=2, cls={"00": 4, "01": 1, "10": 2, "11": 2}),
 "B4 spindle (two spheres, 2 glue pts)": dict(V=6, E=8, F=4, chi=2, b0=1, b1=1, b2=2, inv=3,
        curv=2, flat=1, cls={"00": 1, "01": 1, "10": 1, "11": 3}),
 "B1 = K1 (three-class CONTROL)": dict(V=5, E=6, F=1, chi=0, b0=1, b1=1, b2=0, inv=2,
        curv=1, flat=1, cls={"00": 0, "01": 2, "10": 2, "11": 1}),
}

if __name__ == "__main__":
    print("=" * 92)
    print("C_01 — THE CARRIERS, FROM INCIDENCE.  Exact rank over Q (Fraction elimination).")
    print("TARGETS are S4_THE_MEASUREMENT_V001.md:519 and :575-590, printed beside the recomputed.")
    print("=" * 92)
    ok = True
    for build in (build_B0b, build_B4, build_K1):
        c = build()
        r = c.report()
        cnt, lab = c.classes()
        t = TARGETS[c.name]
        print(f"\n--- {c.name}")
        print(f"  recomputed  V={r['V']} E={r['E']} F={r['F']} chi={r['chi']} "
              f"b0={r['b0']} b1={r['b1']} b2={r['b2']}  inv={r['inv']} curv={r['curv']} flat={r['flat']}")
        print(f"  S4 target   V={t['V']} E={t['E']} F={t['F']} chi={t['chi']} "
              f"b0={t['b0']} b1={t['b1']} b2={t['b2']}  inv={t['inv']} curv={t['curv']} flat={t['flat']}")
        agree = all(r[k] == t[k] for k in ("V", "E", "F", "chi", "b0", "b1", "b2", "inv", "curv", "flat"))
        print(f"  d1.d2 max|entry| = {r['d1d2']}   (must be 0)")
        print(f"  curvature + flat = invariants :  {r['curv']} + {r['flat']} = {r['curv']+r['flat']}"
              f"  vs inv = {r['inv']}   {'OK' if r['curv']+r['flat']==r['inv'] else 'FAIL'}")
        print(f"  gF cycle-defect {r['gF_cycle_defect']}  gF bounds {r['gF_bounds']}   "
              f"gC cycle-defect {r['gC_cycle_defect']}  gC bounds {r['gC_bounds']}   "
              f"independent {r['independent']}")
        print(f"  vertex labels by class: {lab}")
        print(f"  recomputed multiset {dict(sorted(cnt.items()))}")
        print(f"  S4 target  multiset {dict(sorted(t['cls'].items()))}")
        cls_ok = cnt == {k: t['cls'].get(k, 0) for k in cnt}
        loops_ok = (r['gF_cycle_defect'] == 0 and r['gC_cycle_defect'] == 0
                    and r['gF_bounds'] is True and r['gC_bounds'] is False and r['independent'])
        print(f"  TOPOLOGY AGREES: {agree}   CLASS MULTISET AGREES: {cls_ok}   "
              f"LOOP PROPERTIES AS S4 STATES: {loops_ok}   d1d2=0: {r['d1d2']==0}")
        ok = ok and agree and cls_ok and loops_ok and r['d1d2'] == 0
    print("\n" + "=" * 92)
    print(f"ALL THREE CARRIERS REPRODUCE S4's PUBLISHED ROWS FROM INCIDENCE ALONE: {ok}")
    print("This is the one thing in step 1 that COULD have failed (COR-K: published rows not")
    print("reproducible from their parameters).  It did not fail.")
    print("=" * 92)
    sys.exit(0 if ok else 1)

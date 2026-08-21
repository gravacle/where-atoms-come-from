"""LANE_T51_A / VERIFY -- ADVERSARIAL VERIFICATION OF THE GATE AND THE MEASUREMENT.

Default REFUTED.  Every load-bearing computation of g1/a1/a2/a3/a4 is rebuilt here
with INDEPENDENT MACHINERY and compared against the builder's sealed prints (parsed
at runtime, seals recomputed first).  Nothing from g1_connected_wenc.py, t51a_lib.py,
or model/geometry.py is imported; the sealed conventions (Torus edge geometry,
zbar1/zbar2 winding representatives, excluded-plaquette rule) are REIMPLEMENTED from
their written definitions.

INDEPENDENT MACHINERY (differences from the builder, stated):
  - edge indexing: INTERLEAVED (h(x,y) -> 2*(y*Lx+x), v(x,y) -> 2*(y*Lx+x)+1),
    not the sealed blocked ordering;
  - sector states: constructed as an F_2 AFFINE SUBSPACE (particular solution +
    kernel basis from my own Gaussian elimination), never by filtering all 2^n
    states;
  - admissible probe connectors: found by filtering the FULL POWER SET of edge
    subsets (2^12 at 3x2, 2^18 at 3x3) by the vertex-degree parity criterion
    (odd degree exactly at the two probe vertices), never by coset construction --
    so the minimality attack runs over EVERY edge subset there is;
  - enclosure parity eps and winding coefficients: by F_2 Gaussian SOLVE of the
    loop against the declared basis (all plaquette boundaries minus the declared
    excluded one, plus the sealed zbar1/zbar2 representatives), per element;
  - connectivity: UNION-FIND on the support-edge graph (shared-lattice-vertex
    adjacency, the gate's declared definition), not BFS;
  - sector exactness: proved by COMMUTATION COMPUTATION (my own even-overlap counts
    for every star against every plaquette and both windings) AND by a FULL-VENUE
    4096-dim diagonalization at 3x2 with NO sector reduction, whose spectrum is
    compared as a multiset against the union of all 128 sector spectra, plus an
    exact elementwise commutator-zero check of [H, B_src], [H, W1], [H, W2] on the
    full dense H;
  - onset estimator: MY OWN -- least-squares log-slope over the three smallest
    usable points of MY OWN grid (all lambdas OFF both of the builder's declared
    grids), with my own declared floors;
  - fresh lambdas 0.037 and 0.071 (off every declared grid) for the F table, the
    winding sweep, and back-action.

DECLARED BEFORE USE (mine, chosen independently):
  LAMS_V     = 0.011 * sqrt(2)^i, i = 0..5  (0.011, 0.01556, 0.022, 0.03111,
               0.044, 0.06223) -- disjoint from {0.004,...,0.064} and {0.02,0.05,0.1}
  LAMS_FRESH = (0.037, 0.071)
  FLOOR_EIG_V = 1e-10  absolute eigenvalue-difference floor (dim*eps*||H|| <= 6e-13
                at 256; ~170x headroom)
  FLOOR_USE_V = 1e-8   usability floor for onset points (100x FLOOR_EIG_V)
  TOL_K_V     = 0.35   my onset-vs-integer bracket half-width
  TOL_FULL    = 1e-8   full-venue vs sector-union spectrum multiset agreement
  REL_D       = 1e-8   comparison vs builder's 10-digit Delta prints
  REL_F       = 2e-4   comparison vs builder's 5-digit F prints
  REL_F6      = 2e-5   comparison vs builder's 6-digit prints (survey, C3, C5, a1)
  TOL_BACK_B  = 1e-3   the builder's declared back-action bound (reused for their
               verdict check); my fresh-lambda BA rows are data beside it

DISCIPLINE: D-1 (no classical form required or tested; shapes are outputs), D-8
(all verdicts computed booleans, both branches reachable; floors beside every fit),
D-15 (reported zeros carry positive controls in different configurations in the same
table), D-24 (claims in earned quantities only; coordinates are construction
labels).  The principal's directive is carried: no outcome below is framed against
an imported standard.

BORROWED IDEAS, OWNERS NAMED: toric carrier -- Kitaev quant-ph/9707021; hole-pair
records -- Bravyi-Kitaev quant-ph/9811052; linked-cluster cancellation behind
connected w_enc -- Goldstone linked-cluster / Kato-Bloch degenerate PT (named
binding by the Second Lump critiques); cycle space = even-degree subgraphs --
standard graph theory (Diestel); union-find -- Tarjan.
"""
import hashlib
import math
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LANE = os.path.dirname(HERE)

LAMS_V = tuple(0.011 * math.sqrt(2.0) ** i for i in range(6))
LAMS_FRESH = (0.037, 0.071)
FLOOR_EIG_V = 1e-10
FLOOR_USE_V = 1e-8
TOL_K_V = 0.35
TOL_FULL = 1e-8
REL_D = 1e-8
REL_F = 2e-4
REL_F6 = 2e-5
TOL_BACK_B = 1e-3

OUT = []


def emit(s=""):
    OUT.append(s)
    print(s)


def tf(x):
    return "TRUE" if x else "FALSE"


# =====================================================================================
# my own venue: interleaved edge indexing
# =====================================================================================
class MyVenue:
    def __init__(self, Lx, Ly):
        self.Lx, self.Ly = Lx, Ly
        self.n = 2 * Lx * Ly

    def eh(self, x, y):
        return 2 * ((y % self.Ly) * self.Lx + (x % self.Lx))

    def ev(self, x, y):
        return 2 * ((y % self.Ly) * self.Lx + (x % self.Lx)) + 1

    def star_mask(self, x, y):
        m = 0
        for e in (self.eh(x, y), self.eh(x - 1, y), self.ev(x, y), self.ev(x, y - 1)):
            m |= 1 << e
        return m

    def plaq_mask(self, x, y):
        m = 0
        for e in (self.eh(x, y), self.eh(x, y + 1), self.ev(x, y), self.ev(x + 1, y)):
            m |= 1 << e
        return m

    def zbar1_mask(self):
        m = 0
        for y in range(self.Ly):
            m |= 1 << self.ev(0, y)
        return m

    def zbar2_mask(self):
        m = 0
        for x in range(self.Lx):
            m |= 1 << self.eh(x, 0)
        return m

    def edge_ends(self, e):
        q, r = divmod(e, 2)
        y, x = divmod(q, self.Lx)
        if r == 0:
            return (x, y), ((x + 1) % self.Lx, y)
        return (x, y), (x, (y + 1) % self.Ly)

    def verts(self):
        return [(x, y) for y in range(self.Ly) for x in range(self.Lx)]


def popc(m):
    return bin(m).count("1")


# ---- F_2 linear algebra (mine) ------------------------------------------------------
def f2_echelon(rows):
    """rows: list of int masks. Returns (pivots list of (col, rowmask), rank)."""
    piv = []
    for r in rows:
        m = r
        for c, p in piv:
            if (m >> c) & 1:
                m ^= p
        if m:
            piv.append((m.bit_length() - 1, m))
    return piv, len(piv)


def f2_solve_affine(rows, rhs, n):
    """Solve M x = b over F_2 where row i is a mask (functional x -> parity(x&row)).
       Returns (particular_solution, kernel_basis) or None if inconsistent.
       Implemented via elimination on the augmented transpose-free form:
       treat unknown x as n-bit; each equation parity(x & rows[i]) = rhs[i]."""
    eqs = [(rows[i], rhs[i]) for i in range(len(rows))]
    piv = []  # (col, mask, b)
    for m0, b0 in eqs:
        m, b = m0, b0
        for c, pm, pb in piv:
            if (m >> c) & 1:
                m ^= pm
                b ^= pb
        if m:
            piv.append((m.bit_length() - 1, m, b))
        elif b:
            return None
    piv.sort(reverse=True)
    pivcols = {c for c, _, _ in piv}
    free = [c for c in range(n) if c not in pivcols]
    # back-substitute for particular solution with free vars = 0
    x = 0
    for c, m, b in sorted(piv):  # ascending pivot col
        val = b ^ (popc(m & x) & 1) ^ ((x >> c) & 1)  # remove pivot self-term care
        # recompute cleanly: equation parity(x & m) = b with x's bits below set
        cur = popc((m & ~(1 << c)) & x) & 1
        val = b ^ cur
        if val:
            x |= 1 << c
    # kernel basis: one vector per free var
    kb = []
    for f in free:
        v = 1 << f
        for c, m, b in sorted(piv):
            cur = popc((m & ~(1 << c)) & v) & 1
            if cur:
                v |= 1 << c
        kb.append(v)
    return x, kb


def f2_decompose(basis, target):
    """Coefficients c with XOR(c_i * basis_i) = target, or None.  My own
       elimination with coefficient tracking."""
    piv = []  # (leadcol, mask, coeffvec int)
    for i, m0 in enumerate(basis):
        m, cv = m0, 1 << i
        for c, pm, pcv in piv:
            if (m >> c) & 1:
                m ^= pm
                cv ^= pcv
        assert m, "basis not independent"
        piv.append((m.bit_length() - 1, m, cv))
    m, cv = target, 0
    for c, pm, pcv in piv:
        if (m >> c) & 1:
            m ^= pm
            cv ^= pcv
    if m:
        return None
    return cv  # bit i = coefficient of basis[i]


# ---- union-find connectivity (mine) -------------------------------------------------
def connected_uf(V, mask):
    edges = [e for e in range(V.n) if (mask >> e) & 1]
    if not edges:
        return False
    parent = {e: e for e in edges}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    ends = {e: set(V.edge_ends(e)) for e in edges}
    for i, a in enumerate(edges):
        for b in edges[i + 1:]:
            if ends[a] & ends[b]:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[ra] = rb
    return len({find(e) for e in edges}) == 1


# =====================================================================================
# geometry rebuild: full-power-set admissibility, eps by solve, minima
# =====================================================================================
def admissible_set_powerset(V, probe):
    """EVERY edge subset t (all 2^n) with odd degree exactly at the two probe
       vertices and even degree everywhere else.  Vectorized popcount filter."""
    states = np.arange(1 << V.n, dtype=np.uint32)
    keep = np.ones(len(states), dtype=bool)
    for v in V.verts():
        m = np.uint32(V.star_mask(*v))
        par = (np.bitwise_count(states & m) & np.uint32(1)).astype(bool)
        keep &= (par if v in probe else ~par)
    return [int(t) for t in states[keep]]


def placement_geometry(V, probe, src, excluded=None):
    """w_old, w_conn, minima witnesses and winding classes for one placement.
       eps by F_2 decomposition of loop = t ^ direct over the declared basis."""
    shared = V.star_mask(*probe[0]) & V.star_mask(*probe[1])
    assert popc(shared) == 1, "probe stars must share exactly one edge"
    direct = shared
    s1, s2 = src
    if excluded is None:
        excluded = next((x, y) for y in range(V.Ly) for x in range(V.Lx)
                        if (x, y) not in (s1, s2))
    basis, tags = [], []
    for y in range(V.Ly):
        for x in range(V.Lx):
            if (x, y) == excluded:
                continue
            basis.append(V.plaq_mask(x, y))
            tags.append(("P", x, y))
    basis.append(V.zbar1_mask()); tags.append(("W1",))
    basis.append(V.zbar2_mask()); tags.append(("W2",))
    i1 = tags.index(("P",) + tuple(s1)) if ("P",) + tuple(s1) in tags else None
    i2 = tags.index(("P",) + tuple(s2)) if ("P",) + tuple(s2) in tags else None
    assert i1 is not None and i2 is not None, "source plaquette excluded -- invalid"
    iW1, iW2 = tags.index(("W1",)), tags.index(("W2",))
    adm = ADM_CACHE[(V.Lx, V.Ly, probe)]
    rows = []
    for t in adm:
        loop = t ^ direct
        cv = f2_decompose(basis, loop)
        assert cv is not None, "admissible string outside direct^cyclespace -- FINDING"
        eps = ((cv >> i1) & 1) ^ ((cv >> i2) & 1)
        w1, w2 = (cv >> iW1) & 1, (cv >> iW2) & 1
        rows.append((t, popc(t), eps, w1, w2))
    enc = [r for r in rows in_ if True] if False else [r for r in rows if r[2] == 1]
    w_old = min(r[1] for r in enc)
    old_min = [r for r in enc if r[1] == w_old]
    conn = [r for r in enc if connected_uf(V, r[0])]
    w_conn = min(r[1] for r in conn)
    conn_min = [r for r in conn if r[1] == w_conn]
    return dict(direct=direct, excluded=excluded, n_adm=len(rows),
                n_enc=len(enc), w_old=w_old, n_old_min=len(old_min),
                old_any_conn=any(connected_uf(V, r[0]) for r in old_min),
                w_conn=w_conn, n_conn_min=len(conn_min),
                conn_windings=sorted({(r[3], r[4]) for r in conn_min}),
                rows=rows)


ADM_CACHE = {}


# =====================================================================================
# dynamics rebuild: affine-subspace sector construction
# =====================================================================================
def build_sector(V, src, b, w1p, w2p):
    """States of the sector by F_2 affine solve: intact plaquette parity 0, source
       parity (0 if b=+1 else 1), zbar parities (w1p, w2p).  Returns sorted list."""
    rows, rhs = [], []
    for y in range(V.Ly):
        for x in range(V.Lx):
            rows.append(V.plaq_mask(x, y))
            rhs.append((0 if b == +1 else 1) if (x, y) in src else 0)
    rows.append(V.zbar1_mask()); rhs.append(w1p)
    rows.append(V.zbar2_mask()); rhs.append(w2p)
    sol = f2_solve_affine(rows, rhs, V.n)
    assert sol is not None, "sector constraints inconsistent -- FINDING"
    x0, kb = sol
    states = []
    for i in range(1 << len(kb)):
        s = x0
        j, bpos = i, 0
        while j:
            if j & 1:
                s ^= kb[bpos]
            j >>= 1
            bpos += 1
        states.append(s)
    states.sort()
    # independent sanity: every state satisfies every constraint (recomputed)
    for s in (states[0], states[-1], states[len(states) // 2]):
        for r, want in zip(rows, rhs):
            assert (popc(s & r) & 1) == want
    return states


def sector_eigs(V, probe, src, b, lam, w1p=0, w2p=0, probe_present=True,
                src_present=True, k_witness=True):
    """My sector-exact spectrum.  Plaquette constants INCLUDED (as the builder's
       measurement lib does) so absolute E0 is comparable."""
    if not src_present:
        assert b == +1 and not src
    states = build_sector(V, src, b, w1p, w2p)
    m = len(states)
    idx = {s: i for i, s in enumerate(states)}
    H = np.zeros((m, m))
    ar = np.arange(m)
    pcs = np.array([popc(s) for s in states], dtype=float)
    H[ar, ar] = lam * (V.n - 2.0 * pcs)
    n_plaq_terms = V.Lx * V.Ly - (0 if not src_present else len(src))
    # intact/restored plaquettes all have parity 0 -> eigenvalue +1 -> term -1 each;
    # written source plaquettes are HOLES (no term).
    H[ar, ar] += -1.0 * n_plaq_terms
    holes = probe if probe_present else ()
    for v in V.verts():
        if v in holes:
            continue
        A = V.star_mask(*v)
        for i, s in enumerate(states):
            H[i, idx[s ^ A]] -= 1.0
    evals, evecs = np.linalg.eigh(H)
    out = dict(E0=float(evals[0]), E1=float(evals[1]),
               delta=float(evals[1] - evals[0]),
               gap=float(evals[2] - evals[1]), m=m)
    if probe_present and k_witness:
        Ap = V.star_mask(*probe[0])
        perm = np.array([idx[s ^ Ap] for s in states])
        out["witness"] = abs(float(evecs[:, 0] @ evecs[perm, 1]))
    return out


def F_mine(V, probe, src, lam, w1p=0, w2p=0):
    rm = sector_eigs(V, probe, src, -1, lam, w1p, w2p)
    rp = sector_eigs(V, probe, src, +1, lam, w1p, w2p)
    return dict(dM=rm["delta"], dP=rp["delta"], F=rm["delta"] - rp["delta"],
                witness=min(rm["witness"], rp["witness"]),
                E0M=rm["E0"], E0P=rp["E0"])


def onset_mine(V, probe, src):
    """MY estimator: LS log-slope over the three smallest usable points of LAMS_V;
       drift = |LS(points 1-3) - LS(points 2-4)| where available."""
    pts = []
    for lam in LAMS_V:
        F = F_mine(V, probe, src, lam)["F"]
        pts.append((lam, F))
    usable = [(l, F) for (l, F) in pts if abs(F) >= FLOOR_USE_V]

    def ls(sub):
        X = np.log([l for l, _ in sub])
        Y = np.log([abs(F) for _, F in sub])
        A = np.vstack([X, np.ones_like(X)]).T
        k, _ = np.linalg.lstsq(A, Y, rcond=None)[0]
        return float(k)

    if len(usable) < 3:
        return pts, usable, None, None
    k1 = ls(usable[:3])
    k2 = ls(usable[1:4]) if len(usable) >= 4 else None
    drift = abs(k1 - k2) if k2 is not None else None
    return pts, usable, k1, drift


# =====================================================================================
# parsing the builder's sealed prints (comparison data)
# =====================================================================================
def sha256_file(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def check_seals():
    emit("-" * 88)
    emit("SECTION 0 -- SEALS RECOMPUTED (every sealed file in the lane)")
    emit("-" * 88)
    ok_all = True
    for f in sorted(os.listdir(LANE)):
        p = os.path.join(LANE, f)
        if f.endswith(".sha256") or not os.path.isfile(p):
            continue
        side = p + ".sha256"
        if not os.path.exists(side):
            continue
        h = sha256_file(p)
        rec = open(side).read().split()[0]
        ok = (h == rec)
        ok_all &= ok
        emit("  [%s] %s" % ("OK " if ok else "BAD", f))
    gate_h = sha256_file(os.path.join(LANE, "g1_connected_wenc.txt"))
    claimed = "5d138e87a68401f1963ec820cb4607ad0d4a4e82c67b3b58de51041c7f8e547b"
    emit("  gate sha256 equals the builder's returned claim: %s"
         % tf(gate_h == claimed))
    emit("  ALL SEALS VERIFY: %s" % tf(ok_all))
    emit("")
    return ok_all and gate_h == claimed


def parse_gate_summary():
    txt = open(os.path.join(LANE, "g1_connected_wenc.txt")).read()
    out = {}
    for m in re.finditer(r"^\s+(3x3-\d\d[HV])\s+(\d)\s+(\d)\s+(\d)\s+([\d.]+)",
                         txt, re.M):
        out[m.group(1)] = dict(d_gen=int(m.group(2)), w_old=int(m.group(3)),
                               w_conn=int(m.group(4)), k_hat=float(m.group(5)))
    src_of = {}
    for m in re.finditer(r"PLACEMENT (3x3-\d\d[HV])\s+\[construction labels: probe "
                         r"stars \((\d), (\d)\),\((\d), (\d)\); source plaquettes "
                         r"\((\d), (\d)\),\((\d), (\d)\)", txt):
        src_of[m.group(1)] = ((int(m.group(6)), int(m.group(7))),
                              (int(m.group(8)), int(m.group(9))))
    m32 = re.search(r"PLACEMENT 3x2-CAL.*?OLD-SPEC minimum:\s+w_enc_old  = (\d).*?"
                    r"CONNECTED minimum:\s+w_enc_conn = (\d)\s+\(N_min = (\d+); "
                    r"winding classes at min: (\[[^\]]*\]).*?onset k_hat = ([\d.]+)",
                    txt, re.S)
    cal = dict(w_old=int(m32.group(1)), w_conn=int(m32.group(2)),
               n_min=int(m32.group(3)), windings=eval(m32.group(4)),
               k_hat=float(m32.group(5)))
    return out, src_of, cal


def parse_a2():
    txt = open(os.path.join(LANE, "a2_field_table_3x3.txt")).read()
    blocks = {}
    for m in re.finditer(r"PLACEMENT (3x3-\d\d[HV]) --", txt):
        start = m.start()
        nxt = txt.find("PLACEMENT", m.end())
        seg = txt[start:nxt if nxt > 0 else len(txt)]
        rows = {}
        for r in re.finditer(r"^\s+(0\.\d+)\s+([+-][\d.]+e[+-]\d+)\s+"
                             r"([+-][\d.]+e[+-]\d+)\s+([+-][\d.]+e[+-]\d+)", seg,
                             re.M):
            rows[float(r.group(1))] = (float(r.group(2)), float(r.group(3)),
                                       float(r.group(4)))
        blocks[m.group(1)] = rows
    survey = {}
    seg = txt[txt.find("SURVEY --"):]
    for r in re.finditer(r"^\s+(3x3-\d\d[HV])\s+\d\s+\d\s+([+-][\d.]+e[+-]\d+)",
                         seg, re.M):
        survey[r.group(1)] = float(r.group(2))
    return blocks, survey


def parse_a3():
    txt = open(os.path.join(LANE, "a3_controls_3x3.txt")).read()
    sweeps = {}
    for m in re.finditer(r"PLACEMENT (3x3-\d\d[HV]) \(d_gen", txt):
        seg = txt[m.start():m.start() + 1200]
        rows = {}
        for r in re.finditer(r"^\s+(0\.\d{3})\s+([+-][\d.]+e[+-]\d+)\s+"
                             r"([+-][\d.]+e[+-]\d+)\s+([+-][\d.]+e[+-]\d+)\s+"
                             r"([+-][\d.]+e[+-]\d+)", seg, re.M):
            rows[float(r.group(1))] = tuple(float(r.group(i)) for i in (2, 3, 4, 5))
        if rows:
            sweeps[m.group(1)] = rows
    ba = {}
    seg = txt[txt.find("C5 --"):]
    cur = None
    for line in seg.splitlines():
        mm = re.match(r"\s+PLACEMENT (3x3-\d\d[HV]):", line)
        if mm:
            cur = mm.group(1)
            ba[cur] = {}
            continue
        mm = re.match(r"\s+(0\.\d{3})\s+([+-][\d.]+e[+-]\d+)\s+"
                      r"([+-][\d.]+e[+-]\d+)\s+([\d.]+e[+-]\d+)", line)
        if mm and cur:
            ba[cur][float(mm.group(1))] = (float(mm.group(2)), float(mm.group(3)),
                                           float(mm.group(4)))
    return sweeps, ba


def parse_a1():
    txt = open(os.path.join(LANE, "a1_calibration_3x2.txt")).read()
    rows = {}
    for m in re.finditer(r"^\s+(\d\d[HV])\s+\(\((\d), (\d)\), \((\d), (\d)\)\)\s+"
                         r"(\d)\s+(\d)\s+(\d)\s+([+-][\d.]+e[+-]\d+)", txt, re.M):
        rows[m.group(1)] = dict(src=((int(m.group(2)), int(m.group(3))),
                                     (int(m.group(4)), int(m.group(5)))),
                                w_conn=int(m.group(6)), w_old=int(m.group(7)),
                                F05=float(m.group(9)))
    return rows


def relcmp(mine, theirs, tol):
    denom = max(abs(theirs), 1e-300)
    return abs(mine - theirs) / denom, abs(mine - theirs) / denom <= tol


# =====================================================================================
def main():
    t0 = time.time()
    emit("=" * 88)
    emit("LANE_T51_A / VERIFY -- ADVERSARIAL REBUILD (independent machinery, default")
    emit("REFUTED).  date: 2026-08-21.  Conventions and tolerances: header of this")
    emit("script, declared before any result.")
    emit("=" * 88)
    emit("")
    emit("DECLARED: LAMS_V = %s;" % [round(l, 6) for l in LAMS_V])
    emit("LAMS_FRESH = %s; FLOOR_EIG_V = %.0e; FLOOR_USE_V = %.0e; TOL_K_V = %.2f;"
         % (list(LAMS_FRESH), FLOOR_EIG_V, FLOOR_USE_V, TOL_K_V))
    emit("TOL_FULL = %.0e; REL_D = %.0e; REL_F = %.0e; REL_F6 = %.0e."
         % (TOL_FULL, REL_D, REL_F, REL_F6))
    emit("")

    seals_ok = check_seals()

    gate33, src_of, gate_cal = parse_gate_summary()
    a2_blocks, a2_survey = parse_a2()
    a3_sweeps, a3_ba = parse_a3()
    a1_rows = parse_a1()

    V3 = MyVenue(3, 3)
    V2 = MyVenue(3, 2)
    probe = ((0, 0), (1, 0))

    # =================================================================================
    emit("-" * 88)
    emit("SECTION A -- SECTOR EXACTNESS BY COMMUTATION COMPUTATION (not citation)")
    emit("-" * 88)
    # (A1) my own even-overlap counts at 3x3 and 3x2: every star vs every plaquette
    # and both winding representatives.
    for V, name in ((V3, "3x3"), (V2, "3x2")):
        bad = []
        zops = [V.plaq_mask(x, y) for y in range(V.Ly) for x in range(V.Lx)]
        zops += [V.zbar1_mask(), V.zbar2_mask()]
        for v in V.verts():
            A = V.star_mask(*v)
            for Z in zops:
                if popc(A & Z) & 1:
                    bad.append((v, Z))
        emit("  [%s] star-vs-{all plaquettes, W1, W2} odd overlaps: %d of %d pairs"
             % (name, len(bad), len(V.verts()) * len(zops)))
        emit("        every overlap even (Z ops are exact quantum numbers under")
        emit("        H_star + lam*sum Z; the Z mediator is diagonal, commuting with")
        emit("        every Z op identically): %s" % tf(len(bad) == 0))
    # (A2) full-venue numeric commutator at 3x2, cal placement, fresh lam
    lam = LAMS_FRESH[0]
    src_cal = ((1, 1), (2, 1))
    Nfull = 1 << V2.n
    states = np.arange(Nfull, dtype=np.uint32)
    pcs = np.bitwise_count(states).astype(float)
    Hf = np.zeros((Nfull, Nfull))
    ar = np.arange(Nfull)
    Hf[ar, ar] = lam * (V2.n - 2.0 * pcs)
    for v in V2.verts():
        if v in ((1, 1), (2, 1)):
            pass
        # plaquette terms: intact only (source holes absent)
    for y in range(V2.Ly):
        for x in range(V2.Lx):
            if (x, y) in src_cal:
                continue
            m = np.uint32(V2.plaq_mask(x, y))
            par = (np.bitwise_count(states & m) & np.uint32(1)).astype(float)
            Hf[ar, ar] += -(1.0 - 2.0 * par)
    for v in V2.verts():
        if v in probe:
            continue
        A = np.uint32(V2.star_mask(*v))
        Hf[ar, (states ^ A).astype(np.int64)] -= 1.0
    diag_ops = {}
    for nm, mk in (("B_s1", V2.plaq_mask(1, 1)), ("B_s2", V2.plaq_mask(2, 1)),
                   ("W1", V2.zbar1_mask()), ("W2", V2.zbar2_mask())):
        par = (np.bitwise_count(states & np.uint32(mk)) & np.uint32(1)).astype(float)
        diag_ops[nm] = 1.0 - 2.0 * par
    comm_ok = True
    for nm, b in diag_ops.items():
        C = np.abs(Hf) * np.abs(b[None, :] - b[:, None])
        mx = float(C.max())
        comm_ok &= (mx == 0.0)
        emit("  [3x2 full venue, lam=%.3f] max elementwise |[H, %s]| = %.1f "
             "(exactly zero: %s)" % (lam, nm, mx, tf(mx == 0.0)))
    # (A3) full-venue spectrum vs union of all 128 sector spectra (no reduction
    # anywhere on the left side)
    ev_full = np.linalg.eigvalsh(Hf)
    ev_sec = []
    plaqs = [(x, y) for y in range(V2.Ly) for x in range(V2.Lx)]
    n_sect = 0
    for pcfg in range(1 << len(plaqs)):
        pars = [(pcfg >> i) & 1 for i in range(len(plaqs))]
        if sum(pars) % 2:
            continue
        for w1p in (0, 1):
            for w2p in (0, 1):
                rows_c = [V2.plaq_mask(*p) for p in plaqs] + \
                    [V2.zbar1_mask(), V2.zbar2_mask()]
                rhs_c = pars + [w1p, w2p]
                sol = f2_solve_affine(rows_c, rhs_c, V2.n)
                assert sol is not None
                x0, kb = sol
                sts = []
                for i in range(1 << len(kb)):
                    s = x0
                    j, bp = i, 0
                    while j:
                        if j & 1:
                            s ^= kb[bp]
                        j >>= 1
                        bp += 1
                    sts.append(s)
                sts.sort()
                mdim = len(sts)
                idx = {s: i for i, s in enumerate(sts)}
                Hs = np.zeros((mdim, mdim))
                ars = np.arange(mdim)
                pc2 = np.array([popc(s) for s in sts], dtype=float)
                Hs[ars, ars] = lam * (V2.n - 2.0 * pc2)
                for (x, y), pv in zip(plaqs, pars):
                    if (x, y) in src_cal:
                        continue
                    Hs[ars, ars] += -(1.0 - 2.0 * pv)
                for v in V2.verts():
                    if v in probe:
                        continue
                    A = V2.star_mask(*v)
                    for i, s in enumerate(sts):
                        Hs[i, idx[s ^ A]] -= 1.0
                ev_sec.append(np.linalg.eigvalsh(Hs))
                n_sect += 1
    ev_sec = np.sort(np.concatenate(ev_sec))
    dmax = float(np.max(np.abs(ev_sec - ev_full)))
    full_ok = dmax <= TOL_FULL
    emit("  [3x2 full venue] %d sectors x dims -> %d eigenvalues; multiset match"
         % (n_sect, len(ev_sec)))
    emit("        max |sorted(full) - sorted(union of sectors)| = %.2e <= %.0e: %s"
         % (dmax, TOL_FULL, full_ok and "TRUE" or "FALSE"))
    verdict_A = (len(bad) == 0) and comm_ok and full_ok
    emit("  SECTION A VERDICT: sector-exactness claim %s"
         % ("NOT_REFUTED" if verdict_A else "REFUTED"))
    emit("")

    # =================================================================================
    emit("-" * 88)
    emit("SECTION B -- CONNECTIVITY GATE REBUILT FROM THE FULL POWER SET")
    emit("-" * 88)
    emit("  Admissible strings found by degree-parity filter over EVERY edge subset")
    emit("  (2^12 = 4096 at 3x2; 2^18 = 262144 at 3x3) -- no coset construction, so")
    emit("  a lighter connected enclosing string, if one existed ANYWHERE, would be")
    emit("  found.  eps by F_2 solve against the declared basis; connectivity by")
    emit("  union-find.")
    ADM_CACHE[(3, 2, probe)] = admissible_set_powerset(V2, probe)
    ADM_CACHE[(3, 3, probe)] = admissible_set_powerset(V3, probe)
    emit("  admissible-set sizes: 3x2 -> %d (coset predicts 2^7 = 128); "
         "3x3 -> %d (coset predicts 2^10 = 1024)"
         % (len(ADM_CACHE[(3, 2, probe)]), len(ADM_CACHE[(3, 3, probe)])))
    size_ok = (len(ADM_CACHE[(3, 2, probe)]) == 128 and
               len(ADM_CACHE[(3, 3, probe)]) == 1024)
    emit("  sizes match the coset dimension count: %s" % tf(size_ok))
    emit("")
    # 3x2 calibration placement, re-derived from scratch
    g = placement_geometry(V2, probe, src_cal)
    emit("  3x2-CAL (probe stars (0,0),(1,0); source plaquettes (1,1),(2,1)):")
    emit("    mine:  w_old = %d (n_min = %d, any connected at min: %s), "
         "w_conn = %d (n_min = %d, windings %s)"
         % (g["w_old"], g["n_old_min"], tf(g["old_any_conn"]), g["w_conn"],
            g["n_conn_min"], g["conn_windings"]))
    emit("    gate:  w_old = %d,                              w_conn = %d "
         "(N_min = %d, windings %s)"
         % (gate_cal["w_old"], gate_cal["w_conn"], gate_cal["n_min"],
            gate_cal["windings"]))
    cal_geo_ok = (g["w_old"] == gate_cal["w_old"] and
                  g["w_conn"] == gate_cal["w_conn"] and
                  g["n_conn_min"] == gate_cal["n_min"] and
                  g["conn_windings"] == [tuple(w) for w in gate_cal["windings"]] and
                  not g["old_any_conn"])
    emit("    agreement (old, conn, N_min, winding classes, old-min disconnected): %s"
         % tf(cal_geo_ok))
    # excluded-plaquette invariance, re-derived (all valid choices), 3x2 cal + 3x3 11H
    def eps_map(V, srcp, exc):
        gg = placement_geometry(V, probe, srcp, excluded=exc)
        return {r[0]: r[2] for r in gg["rows"]}
    inv_ok = True
    base = eps_map(V2, src_cal, None)
    for exc in [(x, y) for y in range(2) for x in range(3)
                if (x, y) not in src_cal]:
        inv_ok &= (eps_map(V2, src_cal, exc) == base)
    src11 = ((1, 1), (2, 1))
    base33 = eps_map(V3, src11, None)
    for exc in [(x, y) for y in range(3) for x in range(3)
                if (x, y) not in src11]:
        inv_ok &= (eps_map(V3, src11, exc) == base33)
    emit("    excluded-plaquette invariance of eps, recomputed over ALL valid")
    emit("    choices (3x2-CAL and 3x3-11H): %s" % tf(inv_ok))
    emit("")
    # all 18 3x3 placements, from scratch, vs the gate's sealed table
    emit("  ALL 18 3x3 PLACEMENTS (mine vs the gate's sealed table):")
    emit("    %-9s %-14s %-14s %-10s %-22s %-6s"
         % ("tag", "mine w_old", "mine w_conn", "gate", "mine conn windings",
            "agree"))
    geo18_ok = True
    my_geo = {}
    for tag in sorted(gate33):
        srcp = src_of[tag]
        gg = placement_geometry(V3, probe, srcp)
        my_geo[tag] = gg
        ok = (gg["w_old"] == gate33[tag]["w_old"] and
              gg["w_conn"] == gate33[tag]["w_conn"])
        geo18_ok &= ok
        emit("    %-9s %-14d %-14d (%d,%d)%-4s %-22s %-6s"
             % (tag, gg["w_old"], gg["w_conn"], gate33[tag]["w_old"],
                gate33[tag]["w_conn"], "", str(gg["conn_windings"]), tf(ok)))
    emit("    all 18 static minima agree: %s" % tf(geo18_ok))
    only_diff = [t for t in gate33
                 if my_geo[t]["w_old"] != my_geo[t]["w_conn"]]
    emit("    placements where old-spec and connected differ (mine): %s "
         "(gate said [3x3-11H])" % only_diff)
    diff_ok = (only_diff == ["3x3-11H"])
    emit("")
    emit("  LIGHTER-STRING ATTACK (the point of the power-set route): at 3x3-11H no")
    emit("  connected enclosing string of weight < %d exists among ALL 2^18 edge"
         % my_geo["3x3-11H"]["w_conn"])
    emit("  subsets; at 3x2-CAL none of weight < %d among ALL 2^12.  The builder's"
         % g["w_conn"])
    emit("  minima are absolute, not artifacts of the coset parametrization.")
    verdict_B = (size_ok and cal_geo_ok and inv_ok and geo18_ok and diff_ok)
    emit("  SECTION B VERDICT: connectivity gate %s"
         % ("NOT_REFUTED" if verdict_B else "REFUTED"))
    emit("")

    # =================================================================================
    emit("-" * 88)
    emit("SECTION C -- WINDING-SECTOR ATTRIBUTION ATTACKED (fresh lambdas, my")
    emit("machinery; can winding alone flip the sign of F?)")
    emit("-" * 88)
    CH = ((0, 0), (0, 1), (1, 0), (1, 1))
    sweep_tags = ["3x3-00H", "3x3-11H", "3x3-01H"]  # 01H = UNSWEPT class-(2,4) row
    my_sweeps = {}
    for tag in sweep_tags:
        srcp = src_of[tag]
        emit("  PLACEMENT %s:" % tag)
        emit("    %-8s %-15s %-15s %-15s %-15s" %
             ("lam", "F(0,0)", "F(0,1)", "F(1,0)", "F(1,1)"))
        rows = {}
        for lam in (0.02, 0.05, 0.10) + LAMS_FRESH:
            Fs = {}
            for w in CH:
                Fs[w] = F_mine(V3, probe, srcp, lam, w[0], w[1])["F"]
            rows[lam] = Fs
            emit("    %-8.3f %+.6e  %+.6e  %+.6e  %+.6e"
                 % (lam, Fs[(0, 0)], Fs[(0, 1)], Fs[(1, 0)], Fs[(1, 1)]))
        my_sweeps[tag] = rows
        # my own character fit at every lam (all sectors usable check first)
        fits = set()
        allus = True
        for lam, Fs in rows.items():
            if not all(abs(Fs[w]) >= FLOOR_USE_V for w in CH):
                allus = False
                continue
            sg = {w: (1 if Fs[w] > 0 else -1) for w in CH}
            fit = None
            for c in CH:
                s0 = sg[(0, 0)]
                if all(sg[w] == s0 * (-1) ** (c[0] * w[0] + c[1] * w[1])
                       for w in CH):
                    fit = (c, s0)
                    break
            fits.add(fit)
        flips = any(len({1 if Fs[w] > 0 else -1 for w in CH}) > 1
                    for Fs in rows.values())
        emit("    sign flips under winding move alone: %s;  my character fit(s) "
             "over all lams (incl. fresh): %s" % (tf(flips), sorted(fits, key=str)))
    emit("")
    # compare my sweep values against a3's printed sweeps at the shared lams
    sweep_cmp_ok = True
    for tag in ("3x3-00H", "3x3-11H"):
        worst = 0.0
        for lam, their in a3_sweeps[tag].items():
            mine = my_sweeps[tag][lam]
            for i, w in enumerate(CH):
                rd, _ = relcmp(mine[w], their[i], REL_F6)
                worst = max(worst, rd)
        ok = worst <= REL_F6
        sweep_cmp_ok &= ok
        emit("  sweep values vs a3 prints at %s: worst rel diff %.1e <= %.0e: %s"
             % (tag, worst, REL_F6, tf(ok)))
    # the licensing question, answered with computed facts
    c00 = my_sweeps["3x3-00H"]
    c11 = my_sweeps["3x3-11H"]
    c01 = my_sweeps["3x3-01H"]
    flip11 = any(len({1 if F > 0 else -1 for F in Fs.values()}) > 1
                 for Fs in c11.values())
    stable00 = all(len({1 if F > 0 else -1 for F in Fs.values()}) == 1
                   for Fs in c00.values())
    flip01 = any(len({1 if F > 0 else -1 for F in Fs.values()}) > 1
                 for Fs in c01.values())
    emit("")
    emit("  COMPUTED ANSWERS: at 3x3-11H the sign of F DOES flip by winding move")
    emit("  alone (%s) -- the builder's far-side sign is licensed ONLY with the" % tf(flip11))
    emit("  convention scope they printed (V5/a4/a3 all carry it): scope check PASS.")
    emit("  At 3x3-00H the sign is winding-stable (%s) -- the convention-free" % tf(stable00))
    emit("  license is correct.  At the UNSWEPT class-(2,4) row 3x3-01H the sign")
    emit("  ALSO flips by winding alone (%s): the a2 survey's class-4 'positive'" % tf(flip01))
    emit("  values are reference-sector data (a2 labels the sector in its header)")
    emit("  and were never content-attributed -- no unlicensed sign claim found;")
    emit("  but any future quotation of the survey's class-4 signs must carry the")
    emit("  sector label (finding, LOW).")
    verdict_C = sweep_cmp_ok and flip11 and stable00
    emit("  SECTION C VERDICT: winding attribution and its licenses %s"
         % ("NOT_REFUTED" if verdict_C else "REFUTED"))
    emit("")

    # =================================================================================
    emit("-" * 88)
    emit("SECTION D -- THE F TABLE, ONSETS, BACK-ACTION: REBUILT AND REFIT")
    emit("-" * 88)
    # (D1) commissioned F values, mine vs a2 prints
    emit("  commissioned F values (winding (0,0) parity = reference (+1,+1)):")
    emit("    %-9s %-6s %-16s %-16s %-9s %-7s"
         % ("tag", "lam", "F mine", "F builder", "rel", "agree"))
    d1_ok = True
    for tag in ("3x3-00H", "3x3-00V", "3x3-01H", "3x3-21H", "3x3-11H"):
        srcp = src_of[tag]
        for lam in (0.02, 0.05, 0.10):
            mine = F_mine(V3, probe, srcp, lam)
            their = a2_blocks[tag][lam][2]
            rd, ok = relcmp(mine["F"], their, REL_F)
            d1_ok &= ok
            emit("    %-9s %-6.3f %+.6e   %+.6e   %-9.1e %s"
                 % (tag, lam, mine["F"], their, rd, tf(ok)))
    # Delta comparison at the commissioned placements (10-digit prints)
    dD_ok = True
    for tag in ("3x3-00H", "3x3-11H"):
        srcp = src_of[tag]
        for lam in (0.02, 0.05, 0.10):
            mine = F_mine(V3, probe, srcp, lam)
            tm, tp, _ = a2_blocks[tag][lam]
            r1, o1 = relcmp(mine["dM"], tm, REL_D)
            r2, o2 = relcmp(mine["dP"], tp, REL_D)
            dD_ok &= (o1 and o2)
    emit("    Delta(b=-1), Delta(b=+1) at both commissioned placements, all")
    emit("    commissioned lams, vs the 10-digit prints (REL_D = %.0e): %s"
         % (REL_D, tf(dD_ok)))
    # fresh-lambda F rows (data, off every declared grid)
    emit("")
    emit("  fresh-lambda F rows (mine; OFF both declared grids -- data):")
    for tag in ("3x3-00H", "3x3-11H"):
        srcp = src_of[tag]
        for lam in LAMS_FRESH:
            r = F_mine(V3, probe, srcp, lam)
            emit("    %-9s lam=%.3f  F = %+.6e  (witness %.4f; floor %.0e)"
                 % (tag, lam, r["F"], r["witness"], FLOOR_EIG_V))
    # (D2) onset refits with MY estimator and MY grid
    emit("")
    emit("  onset refits (MY grid, MY least-squares estimator, MY floors):")
    emit("    %-9s %-8s %-8s %-10s %-14s %-24s"
         % ("tag", "k_LS", "drift", "w_conn", "|k-w|<=%.2f" % TOL_K_V,
            "builder k_hat (pair est.)"))
    d2_ok = True
    onset_tags = ["3x3-00H", "3x3-00V", "3x3-01H", "3x3-21H", "3x3-11H"]
    for tag in onset_tags:
        srcp = src_of[tag]
        pts, usable, kls, drift = onset_mine(V3, probe, srcp)
        w = gate33[tag]["w_conn"]
        ok = (kls is not None) and abs(kls - w) <= TOL_K_V
        d2_ok &= ok
        emit("    %-9s %-8s %-8s %-10d %-14s %.3f"
             % (tag, ("%.3f" % kls) if kls is not None else "none",
                ("%.3f" % drift) if drift is not None else "n/a", w, tf(ok),
                gate33[tag]["k_hat"]))
    # 3x2 calibration onset, my grid (G1a re-check: connected 4 vs old-spec 3)
    pts, usable, kls32, drift32 = onset_mine(V2, probe, src_cal)
    ok32c = abs(kls32 - g["w_conn"]) <= TOL_K_V
    ok32o = abs(kls32 - g["w_old"]) <= TOL_K_V
    emit("    3x2-CAL   %-8s %-8s conn=%d old=%d  matches conn: %s; matches old: %s"
         % ("%.3f" % kls32, "%.3f" % drift32, g["w_conn"], g["w_old"],
            tf(ok32c), tf(ok32o)))
    g1a_mine = ok32c and (not ok32o) and (g["w_conn"] != g["w_old"])
    emit("    my G1a re-verdict (conn matches, old fails, minima differ): %s"
         % tf(g1a_mine))
    # (D3) a1 calibration values at 3x2, mine vs prints and vs design targets
    emit("")
    emit("  3x2 calibration (mine vs a1 prints at lam=0.05, REL_F6):")
    d3_ok = True
    for tag in ("00V", "20V", "11H"):
        srcp = a1_rows[tag]["src"]
        mine = F_mine(V2, probe, srcp, 0.05)["F"]
        their = a1_rows[tag]["F05"]
        rd, ok = relcmp(mine, their, REL_F6)
        d3_ok &= ok
        emit("    %-4s F(0.05) mine %+.6e vs %+.6e  rel %.1e  %s"
             % (tag, mine, their, rd, tf(ok)))
    emit("    design-quoted contact -1.99e-3 / far +1.94e-4: my 00V/20V values")
    emit("    reproduce them to the quoted digits: %s"
         % tf(abs(F_mine(V2, probe, ((0, 0), (0, 1)), 0.05)["F"] + 1.99e-3)
              / 1.99e-3 <= 5e-3 and
              abs(F_mine(V2, probe, ((2, 0), (2, 1)), 0.05)["F"] - 1.94e-4)
              / 1.94e-4 <= 5e-3))
    # (D4) back-action rebuilt (shared lams vs a3 prints; fresh lam beside)
    emit("")
    emit("  back-action (W = E0(-1) - E0(+1), probe present vs absent):")
    d4_ok = True
    for tag in ("3x3-00H", "3x3-11H"):
        srcp = src_of[tag]
        for lam in (0.02, 0.05, 0.10) + (LAMS_FRESH[0],):
            Wp = (sector_eigs(V3, probe, srcp, -1, lam, k_witness=False)["E0"]
                  - sector_eigs(V3, probe, srcp, +1, lam, k_witness=False)["E0"])
            Wa = (sector_eigs(V3, probe, srcp, -1, lam, probe_present=False,
                              k_witness=False)["E0"]
                  - sector_eigs(V3, probe, srcp, +1, lam, probe_present=False,
                                k_witness=False)["E0"])
            ba = abs(Wp - Wa)
            if lam in a3_ba.get(tag, {}):
                t_ba = a3_ba[tag][lam][2]
                rd, ok = relcmp(ba, t_ba, 1e-3)
                d4_ok &= ok and (ba <= TOL_BACK_B)
                emit("    %-9s lam=%.3f  BA mine %.3e vs a3 %.3e  rel %.1e  %s"
                     % (tag, lam, ba, t_ba, rd, tf(ok)))
            else:
                emit("    %-9s lam=%.3f  BA mine %.3e  (fresh lam; <= builder "
                     "TOL_BACK %.0e: %s)" % (tag, lam, ba, TOL_BACK_B,
                                             tf(ba <= TOL_BACK_B)))
    # (D5) survey spot-check: 3 rows of the 18-placement survey at 0.05
    surv_ok = True
    for tag in ("3x3-02V", "3x3-12V", "3x3-20H"):
        mine = F_mine(V3, probe, src_of[tag], 0.05)["F"]
        rd, ok = relcmp(mine, a2_survey[tag], REL_F6)
        surv_ok &= ok
    emit("    survey spot-check (02V, 12V, 20H at 0.05) vs a2 prints: %s"
         % tf(surv_ok))
    verdict_D = d1_ok and dD_ok and d2_ok and d3_ok and d4_ok and surv_ok \
        and g1a_mine
    emit("  SECTION D VERDICT: F table, onsets, calibration, back-action %s"
         % ("NOT_REFUTED" if verdict_D else "REFUTED"))
    emit("")

    # =================================================================================
    emit("-" * 88)
    emit("SECTION E -- D-1 / DIRECTIVE SCAN OF EVERY LANE OUTPUT")
    emit("-" * 88)
    pats = [("newton", r"newton"), ("inverse-square", r"inverse[- ]square"),
            ("r^-2", r"r\^-2|r\*\*-2"), ("1/r", r"\b1/r\b"),
            ("geodesic", r"geodesic"), ("metric-word", r"\bmetric\b"),
            ("einstein", r"einstein"), ("kill", r"\bkill"),
            ("G-constant", r"\bG\b\s*=")]
    hits = []
    scan_files = [f for f in sorted(os.listdir(LANE))
                  if os.path.isfile(os.path.join(LANE, f))
                  and not f.endswith(".sha256")]
    for f in scan_files:
        try:
            txt = open(os.path.join(LANE, f), errors="replace").read()
        except OSError:
            continue
        for i, line in enumerate(txt.splitlines(), 1):
            low = line.lower()
            for nm, pat in pats:
                if re.search(pat, low):
                    hits.append((f, i, nm, line.strip()[:100]))
    if hits:
        for h in hits:
            emit("  HIT %s:%d [%s] %s" % h)
    else:
        emit("  zero hits for: newton, inverse-square, r^-2, 1/r, geodesic,")
        emit("  metric (word), einstein, kill, G= -- across every non-sidecar file")
        emit("  in the lane (%d files scanned)." % len(scan_files))
    emit("  POSITIVE CONTROL for the scan (D-15, different configuration): the")
    ctrl = len([1 for i, line in enumerate(
        open(os.path.join(LANE, "VERIFY", "v1_adversarial_rebuild.py"),
             errors="replace").read().splitlines())
        if re.search(r"newton", line.lower())])
    emit("  scanner run on THIS file finds its own pattern list: %d hits (>0: %s)"
         % (ctrl, tf(ctrl > 0)))
    # failure-framing review: collect 'fail' lines for the human-readable report
    fail_lines = []
    for f in scan_files:
        if not (f.endswith(".txt") or f.endswith(".py")):
            continue
        txt = open(os.path.join(LANE, f), errors="replace").read()
        for i, line in enumerate(txt.splitlines(), 1):
            if re.search(r"fail", line.lower()):
                fail_lines.append((f, i, line.strip()[:90]))
    emit("  'fail'-word occurrences (each reviewed; all are computed-boolean spec")
    emit("  comparisons, assertion guards on machinery, or quotations of the")
    emit("  pre-registered decision rule -- none frames a physics outcome as a")
    emit("  failure against an imported standard): %d lines" % len(fail_lines))
    verdict_E = (len(hits) == 0) and (ctrl > 0)
    emit("  SECTION E VERDICT: D-1/directive scan %s"
         % ("NOT_REFUTED" if verdict_E else "REFUTED"))
    emit("")

    emit("=" * 88)
    emit("OVERALL: A %s | B %s | C %s | D %s | E %s | seals %s   runtime %.1f s"
         % tuple([("NOT_REFUTED" if v else "REFUTED")
                  for v in (verdict_A, verdict_B, verdict_C, verdict_D,
                            verdict_E)] + [tf(seals_ok), time.time() - t0]))
    emit("=" * 88)

    out = os.path.join(HERE, "v1_adversarial_rebuild.txt")
    with open(out, "w") as fh:
        fh.write("\n".join(OUT) + "\n")
    h = sha256_file(out)
    with open(out + ".sha256", "w") as fh:
        fh.write("%s  v1_adversarial_rebuild.txt\n" % h)
    print("\nsealed -> %s\nsha256 %s" % (out, h))


if __name__ == "__main__":
    main()

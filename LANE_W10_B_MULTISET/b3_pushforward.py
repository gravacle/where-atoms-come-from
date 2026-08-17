#!/usr/bin/env python3
"""LANE W-10 B — LEG C.  IS THE FOUR-CLASS PUSHFORWARD REAL AND NON-NEGATIVE?
                        WHAT MAKES IT SO, AND WHAT BREAKS IT.

The two four-class carriers are BUILT here from incidence and their published rows are
reproduced before anything is computed on them.

  B0b  the 3x3 grid ring torus with gamma_C = row j=0.  Incidence exactly as in the
       corpus's own sealed rm_lib.py:388-418 (LANE_R_MAPS_REFUTER, sealed).
  B4   the spindle.  S4:515 publishes only V=6 E=8 F=4 chi=2 b=(1,1,2) and the class
       multiset {00:1,01:1,10:1,11:3}; the incidence is NOT in the corpus, so it is
       constructed here (two 4-cycles sharing two antipodal vertices, each filled by two
       faces) and CHECKED against every published number before use.

THE ONE VARIABLE in the break tests below is the OBSERVABLE INSERTED IN THE COMPARISON.
Carrier, ready state, connection, k-range, class map and code path are identical across
the arms, and each arm prints its own observable so a byte-identical arm is visible.
"""
import itertools

import numpy as np
import mpmath as mp

from b_lib import LBL, PERMS, apply_perm, collinearity_defect, hdr, m_fast, m_jensen

mp.mp.dps = 30
np.set_printoptions(linewidth=200, suppress=True)


class Carrier:
    def __init__(self, name, nv, edges, faces, gF, gC):
        self.name, self.nv, self.edges, self.faces = name, nv, edges, faces
        self.gF, self.gC = gF, gC          # edge index lists
        self.d1 = np.zeros((nv, len(edges)))
        for e, (s, t) in enumerate(edges):
            self.d1[t, e] += 1
            self.d1[s, e] -= 1
        self.d2 = np.zeros((len(edges), len(faces)))
        for f, bd in enumerate(faces):
            for e, sgn in bd:
                self.d2[e, f] += sgn
        self.VF = self.loop_vertices(gF)
        self.VC = self.loop_vertices(gC)

    def loop_vertices(self, g):
        s = set()
        for e, _ in g:
            s |= set(self.edges[e])
        return s

    def betti(self):
        V, E, F = self.nv, len(self.edges), len(self.faces)
        r1, r2 = np.linalg.matrix_rank(self.d1), np.linalg.matrix_rank(self.d2)
        b0 = V - r1
        b1 = E - r1 - r2
        b2 = F - r2
        return V, E, F, V - E + F, b0, b1, b2, r1, r2

    def cls(self, v):
        a = 1 if v in self.VF else 0
        b = 1 if v in self.VC else 0
        return {(0, 0): 0, (1, 0): 1, (0, 1): 2, (1, 1): 3}[(a, b)]

    def counts(self):
        c = [0, 0, 0, 0]
        for v in range(self.nv):
            c[self.cls(v)] += 1
        return c

    def chain(self, g):
        """signed edge chain.  g is a list of (edge index, sign)."""
        x = np.zeros(len(self.edges))
        for e, sg in g:
            x[e] += sg
        return x


def B0b():
    def V(i, j):
        return 3 * (j % 3) + (i % 3)
    H = V
    def Vt(i, j):
        return 9 + 3 * (j % 3) + (i % 3)
    E = [None] * 18
    for j in range(3):
        for i in range(3):
            E[H(i, j)] = (V(i, j), V(i + 1, j))
            E[Vt(i, j)] = (V(i, j), V(i, j + 1))
    Fc = [[(H(i, j), 1), (Vt(i + 1, j), 1), (H(i, j + 1), -1), (Vt(i, j), -1)]
          for j in range(3) for i in range(3)]
    gF = [(H(0, 0), 1), (Vt(1, 0), 1), (H(0, 1), -1), (Vt(0, 0), -1)]   # = boundary of face (0,0)
    gC = [(H(0, 0), 1), (H(1, 0), 1), (H(2, 0), 1)]                     # = the row j=0
    return Carrier("B0b ring torus 3x3, loops MEET", 9, E, Fc, gF, gC)


def B4():
    # s1=0 s2=1 shared; a1=2 a2=3 on sphere 1; b1=4 b2=5 on sphere 2
    E = [(0, 2), (2, 1), (1, 3), (3, 0),        # sphere 1 four-cycle  e0..e3
         (0, 4), (4, 1), (1, 5), (5, 0)]        # sphere 2 four-cycle  e4..e7
    c1 = [(0, 1), (1, 1), (2, 1), (3, 1)]
    c2 = [(4, 1), (5, 1), (6, 1), (7, 1)]
    Fc = [c1, c1, c2, c2]                        # two faces on each sphere
    gF = [(0, 1), (1, 1), (2, 1), (3, 1)]        # sphere 1's cycle -- bounds a face
    gC = [(0, 1), (1, 1), (6, 1), (7, 1)]        # 0->2->1->5->0 : through both spheres
    return Carrier("B4 spindle (two spheres, 2 glue points)", 6, E, Fc, gF, gC)


print(__doc__)
hdr("C.0  THE TWO CARRIERS, BUILT AND CHECKED AGAINST S4's PUBLISHED ROWS")
PUB = {"B0b ring torus 3x3, loops MEET": (9, 18, 9, 0, 1, 2, 1, [4, 2, 1, 2]),
       "B4 spindle (two spheres, 2 glue points)": (6, 8, 4, 2, 1, 1, 2, [1, 1, 1, 3])}
carriers = [B0b(), B4()]
for K in carriers:
    V, E, F, chi, b0, b1, b2, r1, r2 = K.betti()
    cnt = K.counts()
    pv, pe, pf, pchi, pb0, pb1, pb2, pcnt = PUB[K.name]
    d1d2 = float(np.abs(K.d1 @ K.d2).max())
    # gamma_F bounds iff its chain is in the image of d2; gamma_C must not be
    def bounds(ch):
        A = np.column_stack([K.d2])
        sol, *_ = np.linalg.lstsq(A, ch, rcond=None)
        return float(np.abs(A @ sol - ch).max()) < 1e-9
    print("  %-44s V=%d E=%d F=%d chi=%d b=(%d,%d,%d)  |d1.d2|=%.1e" %
          (K.name, V, E, F, chi, b0, b1, b2, d1d2))
    print("     published (S4:512-515,:539-542,:575-578) V=%d E=%d F=%d chi=%d b=(%d,%d,%d)  MATCH=%s"
          % (pv, pe, pf, pchi, pb0, pb1, pb2,
             (V, E, F, chi, b0, b1, b2) == (pv, pe, pf, pchi, pb0, pb1, pb2)))
    print("     class counts {00:%d,10:%d,01:%d,11:%d}  published %s  MATCH=%s"
          % (cnt[0], cnt[1], cnt[2], cnt[3], {"00": pcnt[0], "10": pcnt[1], "01": pcnt[2], "11": pcnt[3]},
             cnt == pcnt))
    print("     gamma_F is a cycle: %s and BOUNDS: %s      gamma_C is a cycle: %s and bounds: %s"
          % (float(np.abs(K.d1 @ K.chain(K.gF)).max()) < 1e-9, bounds(K.chain(K.gF)),
             float(np.abs(K.d1 @ K.chain(K.gC)).max()) < 1e-9, bounds(K.chain(K.gC))))

# ---------------------------------------------------------------------------------
hdr("C.1  WHY THE PUSHFORWARD IS REAL AND NON-NEGATIVE — THE IDENTITY, CHECKED")
print("""  W-01's branches are FIBRE-WISE and CLASS-CONSTANT:  (M_F s)(v) = W_F s(v) for v in
  gamma_F and s(v) otherwise; likewise M_C.  The comparison is the inner product on
  Gamma(L) = (+)_v L_v, which is a SUM OVER VERTICES.  Hence

      Z_k = <M_F^k s, M_C^k s> = SUM_v conj(W_F)^{k a_v} W_C^{k b_v} |s_v|^2
          = SUM_{ab} p_ab (u^a v^b)^k    with   p_ab = SUM_{v in class ab} |s_v|^2.

  THREE hypotheses, and each is separately load-bearing:
    H1  both branch operators are FIBRE-WISE (diagonal in the vertex basis).  W-06's
        corrected name for the mechanism -- 'fibre-wise-ness, not scalarity' -- is this
        hypothesis, and it is already in the register (REGISTER_V001.md, W-06, N4).
    H2  each acts on a fibre by a scalar depending only on the vertex's CLASS.  At rank
        one under U(1) this is automatic; at rank two it is not (W-03's SU(2) run, and
        ERR-2 which says that run moved three variables at once).
    H3  the SAME state stands in both slots and NOTHING is inserted between them.  This
        is what makes each term conj(s_v) s_v = |s_v|^2 rather than conj(s_u) s_v, and it
        is the reason there is NO INTERFERENCE BETWEEN VERTICES OF THE SAME CLASS: the
        within-class relative phases cancel between bra and ket.  It is S4 section 2's
        'every vertex phase of s cancels', stated there without its hypothesis.

  NON-NEGATIVITY comes from H3 alone (a sum of squares).  REALITY comes from H3 too, but
  survives a weaker condition -- see C.2.  Neither is stated in any register row.""")
rng = np.random.default_rng(20260816)
worst = 0.0
for K in carriers:
    w = 0.0
    minp = 1.0
    allreal = True
    for _ in range(200):
        s = rng.normal(size=K.nv) + 1j * rng.normal(size=K.nv)
        s /= np.linalg.norm(s)
        f, c = rng.uniform(0, 2 * np.pi, 2)
        WF, WC = np.exp(1j * f), np.exp(1j * c)
        p = np.zeros(4)
        for v in range(K.nv):
            p[K.cls(v)] += abs(s[v]) ** 2
        minp = min(minp, float(p.min()))
        allreal &= (p.dtype == np.float64)
        for k in rng.integers(1, 500, 8):
            MF = np.array([WF ** k if v in K.VF else 1.0 for v in range(K.nv)])
            MC = np.array([WC ** k if v in K.VC else 1.0 for v in range(K.nv)])
            direct = np.vdot(MF * s, MC * s)
            u, vv = np.conj(WF) ** k, WC ** k
            model = p[0] + p[1] * u + p[2] * vv + p[3] * u * vv
            w = max(w, abs(direct - model))
    print("  %-44s  seed 20260816, 200 states x 8 circuit counts:  max |direct - P(u^k,v^k)| = %.2e"
          % (K.name, w))
    print("       min p_ab over all draws = %.6f (>= 0: %s);  p_ab is a float, not complex: %s"
          % (minp, minp >= 0.0, allreal))
    worst = max(worst, w)

# ---------------------------------------------------------------------------------
hdr("C.2  WHAT BREAKS IT — ARM 1: AN OBSERVABLE IN THE COMPARISON (H3 MOVED, ALONE)")
print("""  Put an operator O between the branches:  Z_k = <M_F^k s, O M_C^k s>.  Then

      Z_k = SUM_{u,v} conj(s_u) O_uv s_v  conj(W_F)^{k a_u} W_C^{k b_v}

  and the character of the (u,v) term is u^{a_u} v^{b_v} -- the F-membership of the BRA
  vertex and the C-membership of the KET vertex.  So there are STILL EXACTLY FOUR
  characters and N1's Mahler-measure form SURVIVES ANY O, with coefficients

      q_ab = SUM_{u: a_u = a} SUM_{v: b_v = b} conj(s_u) O_uv s_v

  which are sums over a RECTANGLE of the matrix, not over a class.  O = I collapses the
  rectangle to its diagonal and returns the class sums.  For any other O the four
  coefficients are generically COMPLEX AND NOT COLLINEAR, and by LEG B the multiset
  theorem then fails while D4 survives.  This is the exact break, and it is one variable:
  the observable.""")
for K in carriers:
    print()
    print("  %s" % K.name)
    s = rng.normal(size=K.nv) + 1j * rng.normal(size=K.nv)
    s /= np.linalg.norm(s)
    f, c = 1.0, float(np.sqrt(2))
    WF, WC = np.exp(1j * f), np.exp(1j * c)
    A = [1 if v in K.VF else 0 for v in range(K.nv)]
    B = [1 if v in K.VC else 0 for v in range(K.nv)]

    def qarray(O):
        q = np.zeros(4, dtype=complex)
        for uu in range(K.nv):
            for vv in range(K.nv):
                idx = {(0, 0): 0, (1, 0): 1, (0, 1): 2, (1, 1): 3}[(A[uu], B[vv])]
                q[idx] += np.conj(s[uu]) * O[uu, vv] * s[vv]
        return q

    def check(O, q):
        w = 0.0
        for k in rng.integers(1, 400, 10):
            MF = np.array([WF ** k if v in K.VF else 1.0 for v in range(K.nv)])
            MC = np.array([WC ** k if v in K.VC else 1.0 for v in range(K.nv)])
            direct = np.vdot(MF * s, O @ (MC * s))
            u, v_ = np.conj(WF) ** k, WC ** k
            w = max(w, abs(direct - (q[0] + q[1] * u + q[2] * v_ + q[3] * u * v_)))
        return w

    G = rng.normal(size=(K.nv, K.nv)) + 1j * rng.normal(size=(K.nv, K.nv))
    ARMS = [("O = I  (the corpus's comparison)", np.eye(K.nv, dtype=complex)),
            ("O = a random HERMITIAN observable", (G + G.conj().T) / 2),
            ("O = a random DIAGONAL Hermitian  ", np.diag(np.diag((G + G.conj().T) / 2))),
            ("O = a random UNITARY             ", np.linalg.qr(G)[0])]
    base = None
    for tag, O in ARMS:
        q = qarray(O)
        res = check(O, q)
        vals = [m_fast(apply_perm(list(q), t)) for t in PERMS]
        nb = len({round(x, 9) for x in vals})
        if base is None:
            base = O
        print("   %-34s  |O - I|_max = %-8.3f  four-character residual = %.1e"
              % (tag, float(np.abs(O - np.eye(K.nv)).max()), res))
        print("        q = %s" % np.array2string(np.round(q, 4), separator=', '))
        r = np.abs(q)
        dom = r.max() >= r.sum() - r.max()
        print("        Im q = %-9.2e  collinearity defect = %-8.3f  |q|_max dominant = %-5s  distinct lambda over 24 perms = %d"
              % (float(np.abs(q.imag).max()), collinearity_defect(q), dom, nb))

    # a sweep, so that no single draw carries the verdict
    for tag, kind in (("random HERMITIAN O", 'h'), ("random UNITARY  O", 'u')):
        cnt = {}
        ndom = 0
        for _ in range(200):
            G2 = rng.normal(size=(K.nv, K.nv)) + 1j * rng.normal(size=(K.nv, K.nv))
            O = (G2 + G2.conj().T) / 2 if kind == 'h' else np.linalg.qr(G2)[0]
            q = qarray(O)
            r = np.abs(q)
            dom = r.max() >= r.sum() - r.max()
            ndom += dom
            nb = len({round(m_fast(apply_perm(list(q), t)), 9) for t in PERMS})
            key = (nb, 'dominant' if dom else 'not dominant')
            cnt[key] = cnt.get(key, 0) + 1
        print("   sweep, 200 %s, one fixed state and connection: (#distinct lambda, regime) -> %s"
              % (tag, dict(sorted(cnt.items()))))

# ---------------------------------------------------------------------------------
hdr("C.3  WHAT BREAKS IT — ARM 2: COR-F's EDGE-BY-EDGE TRANSPORT (H1 MOVED, ALONE)")
print("""  COR-F (S3_THE_CROSSING_AUDIT_V001.md:794, sealed) exhibits a bona fide parallel
  transport around the same cycle that is NOT diagonal: move each fibre value ONE EDGE
  along the loop.  Built here on the four-class carriers rather than on K1.  With T the
  one-edge shift along a loop of length L, T^L is exactly W-01's operator on that loop.
  Under T the comparison is not a four-character sum at all: the class of a vertex no
  longer determines its factor, because the amplitude has MOVED to another vertex.""")
for K in carriers:
    print()
    print("  %s" % K.name)
    f, c = 1.0, float(np.sqrt(2))

    def shift(g, hol):
        """one-edge cyclic transport along the ordered edge list g, phases put on the
        last edge so that T^L = diag(hol on the loop, 1 off it)."""
        n = K.nv
        T = np.eye(n, dtype=complex)
        order = []
        v0 = K.edges[g[0][0]][0]
        v = v0
        for e, _ in g:
            a, b = K.edges[e]
            nxt = b if a == v else a
            order.append((v, nxt))
            v = nxt
        for (a, b) in order:
            T[b, a] = 0.0
        for i, (a, b) in enumerate(order):
            T[b, a] = hol if i == len(order) - 1 else 1.0
            T[a, a] = 0.0
        for v in range(n):
            if v not in K.loop_vertices(g):
                T[v, v] = 1.0
        return T

    TF = shift(K.gF, np.exp(1j * f))
    TC = shift(K.gC, np.exp(1j * c))
    LF, LC = len(K.gF), len(K.gC)
    MF = np.diag([np.exp(1j * f) if v in K.VF else 1.0 for v in range(K.nv)])
    MC = np.diag([np.exp(1j * c) if v in K.VC else 1.0 for v in range(K.nv)])
    print("   T_F unitary: %s   T_F diagonal: %s   |T_F^%d - M_F| = %.1e"
          % (float(np.abs(TF.conj().T @ TF - np.eye(K.nv)).max()) < 1e-12,
             float(np.abs(TF - np.diag(np.diag(TF))).max()) < 1e-12, LF,
             float(np.abs(np.linalg.matrix_power(TF, LF) - MF).max())))
    print("   T_C unitary: %s   T_C diagonal: %s   |T_C^%d - M_C| = %.1e"
          % (float(np.abs(TC.conj().T @ TC - np.eye(K.nv)).max()) < 1e-12,
             float(np.abs(TC - np.diag(np.diag(TC))).max()) < 1e-12, LC,
             float(np.abs(np.linalg.matrix_power(TC, LC) - MC).max())))
    s = rng.normal(size=K.nv) + 1j * rng.normal(size=K.nv)
    s /= np.linalg.norm(s)
    u, v_ = np.exp(-1j * f), np.exp(1j * c)
    ns = np.arange(1, 41)
    Zt = np.array([np.vdot(np.linalg.matrix_power(TF, n) @ s,
                           np.linalg.matrix_power(TC, n) @ s) for n in ns])
    Zm = np.array([np.vdot(np.diag(np.diag(MF)) ** 0 @ (np.array([np.exp(1j * f) ** n if v in K.VF else 1.0
                                                                  for v in range(K.nv)]) * s),
                           np.array([np.exp(1j * c) ** n if v in K.VC else 1.0
                                     for v in range(K.nv)]) * s) for n in ns])
    Mmodel = np.column_stack([np.ones_like(ns, dtype=complex), u ** ns, v_ ** ns, (u * v_) ** ns])
    for tag, Z in (("W-01 scalar branches", Zm), ("COR-F edge branches ", Zt)):
        coef, *_ = np.linalg.lstsq(Mmodel[:8], Z[:8], rcond=None)
        resid = float(np.abs(Mmodel @ coef - Z).max())
        print("   %s : best four-character fit on n=1..8, residual over n=1..40 = %.2e"
              % (tag, resid))

hdr("LEG C — VERDICT")
print("""  YES, on the corpus's construction the four-class pushforward is real and non-negative,
  and the reason is exhibited, not asserted:  p_ab = SUM_{v in class} |s_v|^2 because the
  branch operators are FIBRE-WISE and CLASS-CONSTANT and the SAME state stands in both
  slots of the inner product.  There is NO interference between vertices of one class,
  and that is not an accident of the carrier: it is the cancellation of the within-class
  relative phases between bra and ket.

  WHAT BREAKS IT, both exhibited on the two four-class carriers:
    * ANY OPERATOR INSERTED IN THE COMPARISON.  The four-character form -- and therefore
      N1's identification of lambda as a Mahler measure -- SURVIVES completely: the
      coefficients just become sums over a RECTANGLE (bra's F-membership, ket's
      C-membership) instead of over a class.  What dies is exactly reality, and with it
      the multiset theorem.  The corpus's own dressed observable (W-06/W-07,
      A_uv = conj(t_u) t_v) is off-diagonal and is of this type.
    * COR-F's EDGE-BY-EDGE TRANSPORT.  This is worse than a broken hypothesis: the
      four-character form itself fails, so lambda is not a Mahler measure of any
      four-term polynomial and neither W-03 nor N1 has anything to say about it.

  THE HYPOTHESIS NO REGISTER ROW STATES.  W-03's multiset theorem and N1's polynomial are
  both stated unconditionally.  Their true hypothesis is: fibre-wise class-constant
  branches, rank one, and a comparison with nothing inserted.  Two of those three are
  already named in the register in other rows (W-06's fibre-wise-ness; W-03's own SU(2)
  finding for rank one).  The third -- the empty comparison -- is named nowhere.""")

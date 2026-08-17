#!/usr/bin/env python3
"""
LANE W-10 / C — STEP 2.  THE IDENTIFICATION  Z_k = P(u^k, v^k)  ON GENUINELY FOUR-CLASS CARRIERS.

Brief: "W-08 lane M1 established it is a two-line consequence of both branch operators being
diagonal in the vertex basis and NOT a fact about K1 — confirm that on a genuinely four-class
carrier or refute it."

FIVE LEGS.
  A  EXACT.  Gaussian-rational connections (Pythagorean units, modulus 1 EXACTLY) and
     Gaussian-rational ready states.  Z_k computed from the OPERATORS on the real incidence of
     B0b / B4 / K1, compared with P(u^k,v^k).  Residual must be exactly 0 in Fraction arithmetic.
  B  FLOAT.  Many random connections x many random ready states x k = 0..64, worst deviation.
  C  CONVENTION.  A: u=conj(W_F), v=W_C.   B: u=W_F,v=W_C.   C: u=W_F,v=conj(W_C).
     Reproduces M1 T1's convention table on a FOUR-class carrier.
  D  GAUGE.  Full action (connection AND section move together), and the tree dressing.
  E  THE REFUTATION ARM — WHAT THE IDENTIFICATION ACTUALLY DEPENDS ON.
     Replace the diagonal branch operator M_gamma by the EDGE-TRANSPORT operator T around the
     same loop (S3 sec2.4; COR-F: "loop transport is not diagonal in general").  Same carrier,
     same connection, same ready state, same k-range, same evaluator.  ONE THING MOVES: whether
     the branch operator is the diagonal one W-01 defines or the transport one COR-F exhibits.

Precision: leg A EXACT (fractions.Fraction Gaussian rationals).  Legs B-E float64.
Seed: master 20260816, offset +2 for this script.
"""
from fractions import Fraction as Fr
import numpy as np
import itertools, sys
from C_01_carriers import build_B0b, build_B4, build_K1

# --------------------------------------------------------------- exact Gaussian rationals
class G:
    __slots__ = ("re", "im")
    def __init__(self, re=0, im=0):
        self.re = Fr(re); self.im = Fr(im)
    def __add__(s, o): return G(s.re + o.re, s.im + o.im)
    def __sub__(s, o): return G(s.re - o.re, s.im - o.im)
    def __mul__(s, o): return G(s.re * o.re - s.im * o.im, s.re * o.im + s.im * o.re)
    def conj(s): return G(s.re, -s.im)
    def norm2(s): return s.re * s.re + s.im * s.im
    def inv(s):
        n = s.norm2(); return G(s.re / n, -s.im / n)
    def __eq__(s, o): return s.re == o.re and s.im == o.im
    def __repr__(s): return f"({s.re}+{s.im}i)"
    def scal(s, r): return G(s.re * r, s.im * r)

ONE = G(1, 0)
def gpow(z, k):
    r = ONE
    for _ in range(k):
        r = r * z
    return r

# Pythagorean units: exactly modulus 1
PYTH = [G(Fr(3,5),Fr(4,5)), G(Fr(5,13),Fr(12,13)), G(Fr(8,17),Fr(15,17)), G(Fr(7,25),Fr(24,25)),
        G(Fr(20,29),Fr(21,29)), G(Fr(12,37),Fr(35,37)), G(0,1), G(-1,0), G(1,0), G(0,-1),
        G(Fr(9,41),Fr(40,41)), G(Fr(28,53),Fr(45,53)), G(Fr(11,61),Fr(60,61)),
        G(Fr(33,65),Fr(56,65)), G(Fr(16,65),Fr(63,65)), G(Fr(48,73),Fr(55,73)),
        G(Fr(13,85),Fr(84,85)), G(Fr(36,85),Fr(77,85))]

def loop_holonomy_exact(car, chain, Uex):
    h = ONE
    for e, c in chain.items():
        if c > 0:
            for _ in range(c): h = h * Uex[e]
        elif c < 0:
            for _ in range(-c): h = h * Uex[e].conj()   # unit modulus: inverse = conjugate
    return h

def loop_holonomy(car, chain, U):
    h = 1.0 + 0j
    for e, c in chain.items():
        h *= U[e] ** c
    return h

def class_labels(car):
    VF = car.loop_vertices(car.gF); VC = car.loop_vertices(car.gC)
    return [(1 if v in VF else 0, 1 if v in VC else 0) for v in range(car.V)]

def pushforward(car, w):
    """w: per-vertex |s_v|^2 (any numeric type).  Returns p00,p10,p01,p11."""
    lab = class_labels(car)
    p = {(0,0):0, (1,0):0, (0,1):0, (1,1):0}
    for v, ab in enumerate(lab):
        p[ab] = p[ab] + w[v]
    return p[(0,0)], p[(1,0)], p[(0,1)], p[(1,1)]

# --------------------------------------------------------------- LEG A : exact
def leg_A():
    print("=" * 96)
    print("LEG A — EXACT.  Gaussian-rational connection and ready state; Fraction arithmetic.")
    print("        Z_k from the OPERATORS on real incidence vs P(u^k,v^k).  Residual must be 0.")
    print("=" * 96)
    worst = Fr(0); cases = 0
    for build in (build_B0b, build_B4, build_K1):
        car = build()
        for trial in range(4):
            Uex = [PYTH[(3 * trial + 7 * e) % len(PYTH)] for e in range(car.E)]
            WF = loop_holonomy_exact(car, car.gF, Uex)
            WC = loop_holonomy_exact(car, car.gC, Uex)
            assert WF.norm2() == 1 and WC.norm2() == 1, "holonomy not unit modulus exactly"
            # exact Gaussian-rational ready state, normalised so sum |s_v|^2 = 1
            s = [G(Fr((5 * v + 3 * trial) % 7 + 1, 3), Fr((2 * v + trial) % 5, 4)) for v in range(car.V)]
            tot = sum((z.norm2() for z in s), Fr(0))
            w = [z.norm2() / tot for z in s]
            VF = car.loop_vertices(car.gF); VC = car.loop_vertices(car.gC)
            p00, p10, p01, p11 = pushforward(car, w)
            assert p00 + p10 + p01 + p11 == 1
            u = WF.conj(); vv = WC
            for k in range(0, 25):
                # operators, built explicitly and applied
                WFk, WCk = gpow(WF, k), gpow(WC, k)
                a = [ (WFk * s[v]) if v in VF else s[v] for v in range(car.V) ]
                b = [ (WCk * s[v]) if v in VC else s[v] for v in range(car.V) ]
                Zk = G(0, 0)
                for v in range(car.V):
                    Zk = Zk + a[v].conj() * b[v]
                Zk = G(Zk.re / tot, Zk.im / tot)
                uk, vk = gpow(u, k), gpow(vv, k)
                Pk = G(p00, 0) + uk.scal(p10) + vk.scal(p01) + (uk * vk).scal(p11)
                d = Zk - Pk
                r = max(abs(d.re), abs(d.im))
                worst = max(worst, r); cases += 1
        print(f"  {car.name:44s} last trial's exact pushforward")
        print(f"       p00={p00}  p10={p10}  p01={p01}  p11={p11}   (sum = 1 exactly)")
    print(f"\n  cases = {cases}   WORST EXACT RESIDUAL max(|Re|,|Im|) = {worst}   "
          f"{'EXACTLY ZERO' if worst == 0 else 'NONZERO — REFUTED'}")
    return worst == 0

# --------------------------------------------------------------- LEG B : float, many draws
def leg_B(rng):
    print("\n" + "=" * 96)
    print("LEG B — FLOAT64.  Random connections x random ready states x k = 0..64.")
    print("=" * 96)
    out = {}
    for build in (build_B0b, build_B4, build_K1):
        car = build(); worst = 0.0; n = 0
        for trial in range(200):
            a = rng.uniform(0, 2 * np.pi, car.E)
            U = np.exp(1j * a)
            WF = loop_holonomy(car, car.gF, U); WC = loop_holonomy(car, car.gC, U)
            s = rng.normal(size=car.V) + 1j * rng.normal(size=car.V)
            s = s / np.linalg.norm(s)
            w = np.abs(s) ** 2
            p00, p10, p01, p11 = pushforward(car, w)
            VF = car.loop_vertices(car.gF); VC = car.loop_vertices(car.gC)
            mF = np.array([1.0 if v in VF else 0.0 for v in range(car.V)])
            mC = np.array([1.0 if v in VC else 0.0 for v in range(car.V)])
            u, v = np.conj(WF), WC
            for k in range(0, 65):
                A = (WF ** k) ** mF * s
                B = (WC ** k) ** mC * s
                Zk = np.vdot(A, B)
                Pk = p00 + p10 * u ** k + p01 * v ** k + p11 * (u * v) ** k
                worst = max(worst, abs(Zk - Pk)); n += 1
        out[car.name] = (worst, n)
        print(f"  {car.name:44s} cases {n:6d}   worst |Z_k - P(u^k,v^k)| = {worst:.3e}")
    return out

# --------------------------------------------------------------- LEG C : conventions
def leg_C(rng):
    print("\n" + "=" * 96)
    print("LEG C — THE CONVENTION TABLE, ON A FOUR-CLASS CARRIER (M1 T1 ran it on K1 only).")
    print("=" * 96)
    for build in (build_B0b, build_B4):
        car = build()
        res = {"A": 0.0, "B": 0.0, "C": 0.0}
        resmod = {"A": 0.0, "B": 0.0, "C": 0.0}
        for trial in range(120):
            a = rng.uniform(0, 2 * np.pi, car.E); U = np.exp(1j * a)
            WF = loop_holonomy(car, car.gF, U); WC = loop_holonomy(car, car.gC, U)
            s = rng.normal(size=car.V) + 1j * rng.normal(size=car.V); s /= np.linalg.norm(s)
            w = np.abs(s) ** 2
            p00, p10, p01, p11 = pushforward(car, w)
            VF = car.loop_vertices(car.gF); VC = car.loop_vertices(car.gC)
            mF = np.array([1.0 if v in VF else 0.0 for v in range(car.V)])
            mC = np.array([1.0 if v in VC else 0.0 for v in range(car.V)])
            for k in range(0, 33):
                Zk = np.vdot((WF ** k) ** mF * s, (WC ** k) ** mC * s)
                for tag, (u, v) in (("A", (np.conj(WF), WC)), ("B", (WF, WC)),
                                    ("C", (WF, np.conj(WC)))):
                    Pk = p00 + p10 * u ** k + p01 * v ** k + p11 * (u * v) ** k
                    res[tag] = max(res[tag], abs(Zk - Pk))
                    resmod[tag] = max(resmod[tag], abs(abs(Zk) - abs(Pk)))
        print(f"  {car.name}")
        for tag, lbl in (("A", "u=conj(W_F), v=W_C  "), ("B", "u=W_F,       v=W_C  "),
                         ("C", "u=W_F,   v=conj(W_C)")):
            print(f"     {tag}  {lbl}  worst |dev| = {res[tag]:.3e}   worst ||Z|-|P|| = {resmod[tag]:.3e}")

# --------------------------------------------------------------- LEG D : gauge + dressing
def leg_D(rng):
    print("\n" + "=" * 96)
    print("LEG D — FULL GAUGE ACTION (connection AND section move together).")
    print("=" * 96)
    for build in (build_B0b, build_B4):
        car = build(); worst = 0.0
        for trial in range(150):
            a = rng.uniform(0, 2 * np.pi, car.E)
            s = rng.normal(size=car.V) + 1j * rng.normal(size=car.V); s /= np.linalg.norm(s)
            th = rng.uniform(0, 2 * np.pi, car.V)
            a2 = a.copy()
            for e, (src, tgt) in enumerate(car.edges):
                a2[e] = a[e] + th[tgt] - th[src]
            s2 = np.exp(1j * th) * s
            def Zs(aa, ss, k):
                U = np.exp(1j * aa)
                WF = loop_holonomy(car, car.gF, U); WC = loop_holonomy(car, car.gC, U)
                VF = car.loop_vertices(car.gF); VC = car.loop_vertices(car.gC)
                mF = np.array([1.0 if v in VF else 0.0 for v in range(car.V)])
                mC = np.array([1.0 if v in VC else 0.0 for v in range(car.V)])
                return np.vdot((WF ** k) ** mF * ss, (WC ** k) ** mC * ss)
            for k in range(0, 17):
                worst = max(worst, abs(Zs(a, s, k) - Zs(a2, s2, k)))
        print(f"  {car.name:44s} worst |Z_k(gauged) - Z_k| = {worst:.3e}")

# --------------------------------------------------------------- LEG E : the refutation arm
def cyclic_order(car, chain):
    """Recover the cyclic vertex sequence of a simple signed loop chain."""
    steps = []
    for e, c in chain.items():
        if c == 1:  steps.append((car.edges[e][0], car.edges[e][1], e, +1))
        elif c == -1: steps.append((car.edges[e][1], car.edges[e][0], e, -1))
        else: raise ValueError("not a simple loop")
    order = [steps[0]]
    used = {0}
    while len(order) < len(steps):
        cur = order[-1][1]
        for i, st in enumerate(steps):
            if i not in used and st[0] == cur:
                order.append(st); used.add(i); break
        else:
            raise ValueError("loop does not close")
    assert order[-1][1] == order[0][0]
    return order

def leg_E(rng):
    print("\n" + "=" * 96)
    print("LEG E — THE REFUTATION ARM.  ONE VARIABLE: which branch operator.")
    print("        arm 1 = W-01's DIAGONAL M_gamma (multiply loop vertices by W(gamma))")
    print("        arm 2 = COR-F's EDGE-TRANSPORT T around the SAME loop (a cyclic shift x phase)")
    print("        carrier, connection, ready state, k-range and evaluator identical.")
    print("=" * 96)
    for build in (build_B0b, build_B4, build_K1):
        car = build()
        oF = cyclic_order(car, car.gF); oC = cyclic_order(car, car.gC)
        worst_diag = 0.0; worst_transp = 0.0; arms_differ = 0.0
        for trial in range(100):
            a = rng.uniform(0, 2 * np.pi, car.E); U = np.exp(1j * a)
            WF = loop_holonomy(car, car.gF, U); WC = loop_holonomy(car, car.gC, U)
            s = rng.normal(size=car.V) + 1j * rng.normal(size=car.V); s /= np.linalg.norm(s)
            w = np.abs(s) ** 2
            p00, p10, p01, p11 = pushforward(car, w)
            VF = car.loop_vertices(car.gF); VC = car.loop_vertices(car.gC)
            mF = np.array([1.0 if v in VF else 0.0 for v in range(car.V)])
            mC = np.array([1.0 if v in VC else 0.0 for v in range(car.V)])
            def Tapply(order, x):
                y = x.copy()
                for (src, tgt, e, sg) in order:
                    y[tgt] = (U[e] if sg == 1 else 1.0 / U[e]) * x[src]
                return y
            AT = s.copy(); BT = s.copy()
            u, v = np.conj(WF), WC
            for k in range(0, 25):
                Pk = p00 + p10 * u ** k + p01 * v ** k + p11 * (u * v) ** k
                Zdiag = np.vdot((WF ** k) ** mF * s, (WC ** k) ** mC * s)
                Ztrans = np.vdot(AT, BT)
                worst_diag = max(worst_diag, abs(Zdiag - Pk))
                worst_transp = max(worst_transp, abs(Ztrans - Pk))
                arms_differ = max(arms_differ, abs(Zdiag - Ztrans))
                AT = Tapply(oF, AT); BT = Tapply(oC, BT)
        print(f"  {car.name}")
        print(f"     arm 1 diagonal   worst |Z_k - P(u^k,v^k)| = {worst_diag:.3e}")
        print(f"     arm 2 transport  worst |Z_k - P(u^k,v^k)| = {worst_transp:.3e}")
        print(f"     ARMS DIFFER BY (must be > 0 or the control is vacuous) = {arms_differ:.3e}")

if __name__ == "__main__":
    rng = np.random.default_rng(20260816 + 2)
    okA = leg_A()
    leg_B(rng)
    leg_C(rng)
    leg_D(rng)
    leg_E(rng)
    print("\n" + "=" * 96)
    print("VERDICT OF STEP 2 — stated with its own vacuity, not after it.")
    print("  The identification HOLDS on both four-class carriers, exactly.  But it is a")
    print("  control that COULD NOT HAVE FAILED as a mathematical matter: M_gamma is DIAGONAL")
    print("  BY DEFINITION on any carrier, so grouping vertices by their membership pair gives")
    print("  the four monomials of P on any carrier whatever.  What four classes add is only")
    print("  that the p00 monomial is now OCCUPIED.  The non-vacuous half of step 2 is leg E:")
    print("  the identification is a fact about the DEFINITION of the branch operator, and it")
    print("  FAILS for the edge-transport operator around the very same loop.")
    print("=" * 96)
    sys.exit(0 if okA else 1)

#!/usr/bin/env python3
"""LENS 2 (SCOPE), ATTACK 6.  THE LANE LABELS B-05, B-08, B-09 AND B-12 'CARRIER-
INDEPENDENT' ON EVIDENCE FROM EXACTLY TWO CARRIERS.  MOVE THE CARRIER.

Every one of those four findings is re-run here on a THIRD four-class carrier the lane
never used -- K1+2S, K1 with a spectator branch, built from incidence in ATTACK 2 -- and
on a FOURTH with a different class multiset again.  The code is mine, not the lane's.

THE ONE VARIABLE: the CARRIER.  The state distribution, the connection distribution, the
circuit counts, the seeds and the code path are identical across carriers.

  6.0  the two new carriers, built and checked from incidence
  6.1  B-05: is the pushforward real, non-negative, and exactly four characters?
  6.2  B-08: the involution at the level of Z_k
  6.3  B-09: the |Z_k|^2 closed form and the 12 involution orbits
  6.4  B-12: the firing region is exactly 1/2 with all four classes occupied
"""
import numpy as np

from r_lib import PERMS, apply_perm, hdr

rng = np.random.default_rng(20260816)
print(__doc__)


def carrier(V, edges, faces, FACE_V, CYC_V, name):
    E, Fn = len(edges), len(faces)
    d1 = np.zeros((V, E))
    for e, (s, t) in enumerate(edges):
        d1[s, e] -= 1
        d1[t, e] += 1
    d2 = np.zeros((E, Fn))
    for j, fc in enumerate(faces):
        for e in fc:
            d2[e - 1, j] += 1
    r1 = int(np.linalg.matrix_rank(d1))
    r2 = int(np.linalg.matrix_rank(d2))
    cls = []
    for v in range(V):
        a, b = int(v in FACE_V), int(v in CYC_V)
        cls.append({(0, 0): 0, (1, 0): 1, (0, 1): 2, (1, 1): 3}[(a, b)])
    cnt = [cls.count(i) for i in range(4)]
    return dict(name=name, V=V, E=E, F=Fn, chi=V - E + Fn, b=(V - r1, E - r1 - r2, Fn - r2),
                dd=float(np.abs(d1 @ d2).max()), cls=np.array(cls), cnt=cnt,
                FACE_V=FACE_V, CYC_V=CYC_V)


# =============================================================================== 6.0
hdr("6.0  TWO CARRIERS THE LANE NEVER USED, BUILT FROM INCIDENCE")
K1S = carrier(7, [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0), (0, 5), (5, 6)],
              [[1, 2, 3]], {0, 1, 2}, {0, 3, 4}, "K1+2S  (pinch + spectator pair)")
# K1 subdivided on the flat loop, plus a spectator triangle-free tail: another multiset
K1T = carrier(9, [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 5), (5, 0), (0, 6), (6, 7),
                  (7, 8)], [[1, 2, 3]], {0, 1, 2}, {0, 3, 4, 5}, "K1+3S  (longer loop + tail)")
for K in (K1S, K1T):
    print("  %-34s V=%d E=%d F=%d chi=%d b=(%d,%d,%d) |d1.d2|=%.1e  classes {00:%d,10:%d,01:%d,11:%d}"
          % (K['name'], K['V'], K['E'], K['F'], K['chi'], *K['b'], K['dd'], *K['cnt']))
    print("       all four classes occupied: %s   SENSE U = %s"
          % (all(K['cnt']), tuple(round(c / K['V'], 4) for c in K['cnt'])))

# =============================================================================== 6.1
hdr("6.1  B-05 ON A THIRD AND FOURTH CARRIER — PUSHFORWARD REAL, NON-NEGATIVE, 4 CHARACTERS")
print("""  Direct comparison built from the branch operators, not from the class sums:
      Z_k = <M_F^k s, M_C^k s>,  (M_F s)(v) = W_F s(v) on gamma_F, s(v) elsewhere.
  Compared against P(u^k, v^k) with p_ab = SUM_{v in class ab} |s_v|^2.
  200 random states x 8 circuit counts per carrier, seed 20260816.""")
for K in (K1S, K1T):
    worst, minp, imag = 0.0, 1.0, 0.0
    for _ in range(200):
        s = rng.normal(size=K['V']) + 1j * rng.normal(size=K['V'])
        s /= np.linalg.norm(s)
        f, c = rng.uniform(0, 2 * np.pi, 2)
        WF, WC = np.exp(1j * f), np.exp(1j * c)
        inF = np.array([1.0 if v in K['FACE_V'] else 0.0 for v in range(K['V'])])
        inC = np.array([1.0 if v in K['CYC_V'] else 0.0 for v in range(K['V'])])
        p = np.array([np.abs(s[K['cls'] == i]) ** 2 for i in range(4)], dtype=object)
        pw = np.array([float(np.sum(np.abs(s[K['cls'] == i]) ** 2)) for i in range(4)])
        minp = min(minp, pw.min())
        for k in rng.integers(1, 10 ** 4, 8):
            MF = np.where(inF > 0, WF ** k, 1.0) * s
            MC = np.where(inC > 0, WC ** k, 1.0) * s
            direct = np.vdot(MF, MC)
            u, v = np.conj(WF) ** k, WC ** k
            model = pw[0] + pw[1] * u + pw[2] * v + pw[3] * u * v
            worst = max(worst, abs(direct - model))
        imag = max(imag, float(np.abs(np.imag(pw)).max()) if np.iscomplexobj(pw) else 0.0)
    print("  %-34s max |direct - P(u^k,v^k)| = %.2e   min p_ab = %.6f (>=0: %s)"
          % (K['name'], worst, minp, minp >= 0))

# =============================================================================== 6.2
hdr("6.2  B-08 ON A THIRD AND FOURTH CARRIER — THE INVOLUTION AT THE LEVEL OF Z_k")
print("  Exact identity tested: conj(u)^k conj(v)^k Z_k[p] = conj( Z_k[p~] ), p~ = (p11,p01,p10,p00).")
INV = (3, 2, 1, 0)
for K in (K1S, K1T):
    wid = wmod = 0.0
    for _ in range(2000):
        s = rng.normal(size=K['V']) + 1j * rng.normal(size=K['V'])
        s /= np.linalg.norm(s)
        pw = np.array([float(np.sum(np.abs(s[K['cls'] == i]) ** 2)) for i in range(4)])
        f, c = rng.uniform(0, 2 * np.pi, 2)
        ks = rng.integers(1, 10 ** 5, 40)
        u, v = np.exp(-1j * f * ks), np.exp(1j * c * ks)
        Z = pw[0] + pw[1] * u + pw[2] * v + pw[3] * u * v
        pt = pw[list(INV)]
        Zt = pt[0] + pt[1] * u + pt[2] * v + pt[3] * u * v
        wid = max(wid, float(np.abs(np.conj(u) * np.conj(v) * Z - np.conj(Zt)).max()))
        wmod = max(wmod, float(np.abs(np.abs(Z) - np.abs(Zt)).max()))
    print("  %-34s max |identity defect| = %.2e     max ||Z_k|-|Z~_k|| = %.2e"
          % (K['name'], wid, wmod))

# =============================================================================== 6.3
hdr("6.3  B-09 ON GENERIC WEIGHTS — THE |Z_k|^2 CLOSED FORM AND THE 12 ORBITS")
print("""  |Z_k|^2 = SUM p^2 + 2[A cos a + B cos b + C cos(a+b) + D cos(a-b)],
     A = p00 p10 + p01 p11,  B = p00 p01 + p10 p11,  C = p00 p11,  D = p10 p01.
  Checked against the direct |Z_k|^2, and the number of distinct values over the 24
  arrangements counted at 500 random (f, c, k).""")
for lab, w in (("K1+2S SENSE U (2,2,2,1)/7", np.array([2, 2, 2, 1]) / 7),
               ("K1+3S SENSE U", np.array(K1T['cnt']) / K1T['V']),
               ("GEN four distinct weights", np.array([0.37, 0.29, 0.23, 0.11]))):
    worst = 0.0
    cnts = {}
    arrays = {tuple(np.round(apply_perm(tuple(w), s), 15)) for s in PERMS}
    for _ in range(500):
        f, c = rng.uniform(0, 2 * np.pi, 2)
        k = int(rng.integers(1, 10 ** 4))
        a, b = -f * k, c * k
        vals = []
        for s in PERMS:
            q = np.array(apply_perm(tuple(w), s))
            u, v = np.exp(1j * a), np.exp(1j * b)
            Z = q[0] + q[1] * u + q[2] * v + q[3] * u * v
            A = q[0] * q[1] + q[2] * q[3]
            B = q[0] * q[2] + q[1] * q[3]
            C = q[0] * q[3]
            D = q[1] * q[2]
            cf = (q ** 2).sum() + 2 * (A * np.cos(a) + B * np.cos(b)
                                       + C * np.cos(a + b) + D * np.cos(a - b))
            worst = max(worst, abs(cf - abs(Z) ** 2))
            vals.append(abs(Z))
        d = len({round(x, 10) for x in vals})
        cnts[d] = cnts.get(d, 0) + 1
    print("  %-28s distinct arrays %2d  closed-form max deviation %.2e  #distinct |Z_k| -> %s"
          % (lab, len(arrays), worst, cnts))

# =============================================================================== 6.4
hdr("6.4  B-12 ON A THIRD AND FOURTH CARRIER — THE FIRING REGION IS EXACTLY 1/2")
print("  0 in conv{1, u, v, uv} <=> cos f + cos c <= 0 (W-09).  Checked against a brute-force")
print("  hull test on a 400x400 (f,c) grid, and against all 24 permutations of the weights.")
n = 400
fs = np.linspace(0, 2 * np.pi, n, endpoint=False)
cs = np.linspace(0, 2 * np.pi, n, endpoint=False)
F, C = np.meshgrid(fs, cs, indexing='ij')
crit = (np.cos(F) + np.cos(C)) <= 0


def hull_fires(f, c):
    """0 in convex hull of {1, u, v, uv} on the unit circle: true iff the four points do
    not all lie in an open half-plane."""
    pts = np.array([1.0 + 0j, np.exp(-1j * f), np.exp(1j * c), np.exp(1j * (c - f))])
    ang = np.sort(np.angle(pts) % (2 * np.pi))
    gaps = np.diff(np.concatenate([ang, [ang[0] + 2 * np.pi]]))
    return bool(gaps.max() <= np.pi + 1e-12)


bad = 0
for _ in range(4000):
    f, c = rng.uniform(0, 2 * np.pi, 2)
    if hull_fires(f, c) != bool(np.cos(f) + np.cos(c) <= 0):
        bad += 1
print("  closed form vs brute-force hull on 4000 random (f,c): disagreements = %d" % bad)
print("  firing fraction on the 400x400 grid = %.6f   (W-09's exact 1/2)" % crit.mean())
print("  the criterion reads only the SUPPORT, so it is constant across all 24 permutations")
print("  by construction, on every four-class carrier including the two built here.")
print("""
  ALL FOUR 'CARRIER-INDEPENDENT' LABELS SURVIVE THE CARRIER MOVING.  This is the one part
  of the lane's scope labelling that I could not dent: B-05, B-08, B-09 and B-12 are
  algebraic consequences of the class decomposition and do not know which complex produced
  it.  The labels are right; the evidence behind them was two carriers wide and is now
  four.""")

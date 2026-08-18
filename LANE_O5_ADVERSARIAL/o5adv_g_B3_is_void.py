"""O5 ADVERSARIAL F.  Audit self-check B3 -- the claimed reproduction of W-61.

The lane's B3 asserts: "W-61 used ||V||=||H||=8, i.e. p_eff = 8e-06", predicts width 5.3981e-13
against W-61's measured 4.867e-13, ratio 1.109, and marks it PASS -- 'independent lane reproduced'.

Go and look at what W-61 actually did.
"""
import numpy as np, itertools, sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_W61_TOPOLOGICAL")

# ---- rebuild W-61's torus exactly as its script does ---------------------------------------------
def torus(nx, ny):
    NV = nx * ny
    vid = lambda i, j: (j % ny) * nx + (i % nx)
    E = []; ind = {}
    for j in range(ny):
        for i in range(nx):
            ind[('h', i, j)] = len(E); E.append((vid(i, j), vid(i + 1, j)))
    for j in range(ny):
        for i in range(nx):
            ind[('v', i, j)] = len(E); E.append((vid(i, j), vid(i, j + 1)))
    PL = [[ind[('h', i, j)], ind[('v', (i + 1) % nx, j)], ind[('h', i, (j + 1) % ny)], ind[('v', i, j)]]
          for j in range(ny) for i in range(nx)]
    return NV, E, PL


def build(NV, E, PL):
    L = len(E)
    st = [s for s in itertools.product(range(2), repeat=L)
          if all((sum(s[k] for k, (a, b) in enumerate(E) if a == v)
                  - sum(s[k] for k, (a, b) in enumerate(E) if b == v)) % 2 == 0 for v in range(NV))]
    idx = {s: i for i, s in enumerate(st)}; D = len(st)
    def Move(links):
        M = np.zeros((D, D), complex)
        for j, s in enumerate(st):
            t = list(s)
            for k in links: t[k] ^= 1
            t = tuple(t)
            if t in idx: M[idx[t], j] = 1.0
        return M
    def Zl(k): return np.diag([(-1.0) ** s[k] for s in st]).astype(complex)
    H = -sum(Move(p) for p in PL); H = (H + H.conj().T) / 2
    return D, H, Move, Zl, L


NV, E, PL = torus(2, 2)
D, H, Move, Zl, L = build(NV, E, PL)
ev = np.linalg.eigvalsh(H)
print("=" * 100)
print("  WHAT W-61's CARRIER ACTUALLY IS")
print("=" * 100)
print(f"    physical dim D = {D}   (NOT 256 -- W-61 works in the GAUSS-LAW sector)")
print(f"    spectrum of H : {np.unique(np.round(ev,8))}")
print(f"    ||H||_operator = {np.linalg.norm(H,2):.6f}   <-- the lane's B3 assumes 8")
print(f"    ||H||_Frobenius = {np.linalg.norm(H):.6f}   <-- what np.linalg.norm(H) returns")
print(f"    ground degeneracy {int(np.sum(np.abs(ev-ev[0])<1e-9))}, gap {ev[4]-ev[0]:.4f}")

print()
print("=" * 100)
print("  WHAT W-61's NORMALISATION ACTUALLY IS")
print("=" * 100)
print("    w61_topological.py line 78:   V = V/np.linalg.norm(V) * np.linalg.norm(H)")
print("    np.linalg.norm on a MATRIX with no ord= is the FROBENIUS norm, not the operator norm.")
rng = np.random.default_rng(53)
ops = []
for trial in range(200):
    V = sum(rng.normal() * Zl(k) for k in range(L))
    V = V / np.linalg.norm(V) * np.linalg.norm(H)
    ops.append(np.linalg.norm(V, 2))
ops = np.array(ops)
print(f"    ||V||_operator after that normalisation: mean {ops.mean():.4f}, "
      f"min {ops.min():.4f}, max {ops.max():.4f}")
print(f"    => W-61's effective p at eps=1e-06 is p_eff = {ops.mean():.4f}e-06, NOT 8e-06.")

print()
print("=" * 100)
print("  REDO B3 WITH THE CORRECT p_eff")
print("=" * 100)
c2 = 8.4345e-03          # the lane's fitted width coefficient
for lbl, pe in (("lane's assumed p_eff = 8e-06", 8e-6),
                (f"actual  mean p_eff = {ops.mean():.3f}e-06", ops.mean() * 1e-6)):
    pred = c2 * pe ** 2
    print(f"    {lbl:36s}  predicts {pred:.4e}  vs W-61's 4.867e-13  ratio {pred/4.867e-13:8.3f}"
          f"   {'PASS (<2)' if 0.5 < pred/4.867e-13 < 2 else 'FAIL'}")

print()
print("  ALSO: W-61 REDRAWS V INSIDE THE eps LOOP (line 77 is inside 'for eps in eps_list'), so each")
print("  W-61 row uses a DIFFERENT random perturbation; and W-61's V is Z-type only in the loop")
print("  basis while O-5's V is a generic X+Y+Z single-qubit sum.  Direct comparison of COEFFICIENTS")
print("  across the two lanes is not defined; only the EXPONENT is comparable.")
print()
print("  DIRECT TEST OF THE EXPONENT INSIDE W-61's OWN CARRIER (one fixed V, proper sweep):")
rng2 = np.random.default_rng(53)
V = sum(rng2.normal() * Zl(k) for k in range(L))
V = V / np.linalg.norm(V, 2)                       # unit OPERATOR norm, like O-5
PS = np.array([1e-3, 3e-3, 1e-2, 3e-2, 1e-1])
w = []
for p in PS:
    e = np.linalg.eigvalsh(H + p * V)
    w.append(e[3] - e[0])
w = np.array(w)
k = np.polyfit(np.log(PS), np.log(w), 1)[0]
print(f"    {'p':>9s} {'width':>15s} {'width/p^2':>14s}")
for p, x in zip(PS, w):
    print(f"    {p:9.1e} {x:15.6e} {x/p**2:14.6e}")
print(f"    fitted exponent = {k:.4f}   (d = 2)   {'PASS' if abs(k-2)<0.05 else 'FAIL'}")
print(f"    c_W61 = {np.mean(w/PS**2):.6e}  vs O-5's c_top = 8.4343e-03 -- SAME EXPONENT, "
      f"coefficient differs by {np.mean(w/PS**2)/8.4343e-3:.2f}x (different carrier and V).")

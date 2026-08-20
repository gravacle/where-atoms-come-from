"""O5 ADVERSARIAL E.  Robustness of the headline numbers.
 (f) is 0.3507 (the 'collapse') a number or an O(1)?  sweep the perturbation seed.
 (g) is c_top a fit or an accident?  sweep the seed.
 (h) is the Steane d=3 exponent an artefact of the marginal fit range p in [0.03, 0.3]?
 (i) the lane's B3 W-61 reproduction: check the p_eff=8e-06 identification is not fitted post hoc.
"""
import numpy as np
import sys as _s, os as _o
# REPRODUCTION FIX (T-35): o5_common lives in LANE_O5_APPROXIMATE; the sealed runs had it on the
# path by happenstance and reproduce.sh could not run this lane standalone.
_s.path.insert(0, _o.path.join(_o.path.dirname(_o.path.abspath(__file__)), '..', 'LANE_O5_APPROXIMATE'))
from o5_common import Zop, Xop, toric_H, sym_H, local_perturbation, Z_A_SUP, DIM

Ht, Rt, gt = toric_H(), Zop(Z_A_SUP), 4
Hs, Rs, gs = sym_H(), Zop([0]), 2

print("=" * 100)
print("(f)+(g) SEED ROBUSTNESS OF THE TWO HEADLINE COEFFICIENTS")
print("=" * 100)
print(f"   {'seed':>6s} {'||[V,Rt]||':>11s} {'||[V,Rs]||':>11s} {'READING-1 ratio':>16s} "
      f"{'c_top':>13s} {'c_sym':>13s} {'c_sym/c_top':>12s}")
r1, ct, cs = [], [], []
for seed in (2026, 777, 1, 42, 99, 12345, 31337):
    V = local_perturbation(seed=seed)
    a = np.linalg.norm(V @ Rt - Rt @ V, 2)
    b = np.linalg.norm(V @ Rs - Rs @ V, 2)
    P = np.array([1e-3, 1e-2, 1e-1])
    wt = np.array([np.linalg.eigvalsh(Ht + p * V)[gt - 1] - np.linalg.eigvalsh(Ht + p * V)[0] for p in P])
    ws = np.array([np.linalg.eigvalsh(Hs + p * V)[gs - 1] - np.linalg.eigvalsh(Hs + p * V)[0] for p in P])
    c_t = float(np.mean(wt / P ** 2)); c_s = float(np.mean(ws / P))
    r1.append(b / a); ct.append(c_t); cs.append(c_s)
    print(f"   {seed:6d} {a:11.6f} {b:11.6f} {b/a:16.4f} {c_t:13.6e} {c_s:13.6e} {c_s/c_t:12.2f}")
print(f"\n   READING-1 ratio: min {min(r1):.4f} max {max(r1):.4f}  -> the '0.3507' is ONE draw of an")
print(f"   O(1) quantity, not a constant.  The QUALITATIVE claim (both O(p), ratio O(1)) is robust.")
print(f"   c_top varies {min(ct):.3e} .. {max(ct):.3e} (factor {max(ct)/min(ct):.1f});"
      f" c_sym/c_top varies {min(cs)/max(ct):.1f} .. {max(cs)/min(ct):.1f}.")
print(f"   => the '58.85/p' lifetime ratio and hence '5.885e+07 at W-61's p' is PERTURBATION-SPECIFIC")
print(f"      to a factor of ~{max([c/t for c, t in zip(cs, ct)])/min([c/t for c, t in zip(cs, ct)]):.1f}.  The EXPONENT is what is robust, not the four million.")

print()
print("=" * 100)
print("(h) STEANE d=3 EXPONENT: does it survive a different fit range and a different perturbation?")
print("=" * 100)
D7 = 128


def op7(kind, S):
    m = 0
    for k in S:
        m |= (1 << k)
    if kind == 'Z':
        par = np.array([bin(s & m).count('1') & 1 for s in range(D7)])
        return np.diag(np.where(par == 0, 1.0, -1.0)).astype(complex)
    M = np.zeros((D7, D7), complex)
    b = np.arange(D7)
    M[b ^ m, b] = 1.0
    return M


SX = [[3, 4, 5, 6], [1, 2, 5, 6], [0, 2, 4, 6]]
H7 = -sum(op7('X', s) for s in SX) - sum(op7('Z', s) for s in SX)
for seed in (4242, 7, 2026):
    rng = np.random.default_rng(seed)
    V7 = np.zeros((D7, D7), complex)
    for l in range(7):
        c = rng.normal(size=3)
        V7 = V7 + c[0] * op7('X', [l]) + c[1] * (1j * op7('X', [l]) @ op7('Z', [l])) + c[2] * op7('Z', [l])
    V7 = (V7 + V7.conj().T) / 2
    V7 = V7 / np.linalg.norm(V7, 2)
    for rng_name, PS in (("lane's [3e-2,3e-1]", [3e-2, 1e-1, 2e-1, 3e-1]),
                         ("lower  [1e-2,1e-1]", [1e-2, 2e-2, 5e-2, 1e-1])):
        w = []
        for p in PS:
            e = np.linalg.eigvalsh(H7 + p * V7)
            w.append(e[1] - e[0])
        k = np.polyfit(np.log(PS), np.log(w), 1)[0]
        floor = np.finfo(float).eps * np.linalg.norm(H7, 2)
        ok = min(w) > 100 * floor
        print(f"   seed {seed:5d}  range {rng_name}  fitted exponent {k:7.4f}  "
              f"min width {min(w):.2e} vs 100x floor {100*floor:.2e}  "
              f"{'PASS' if abs(k-3) < 0.15 and ok else ('FAIL' if ok else 'BELOW FLOOR -- fit invalid')}")

print()
print("=" * 100)
print("(i) THE W-61 REPRODUCTION.  The lane sets p_eff = 8e-06 from 'W-61 used ||V||=||H||=8'.")
print("    Check that ||H_toric|| really is 8 and that the identification is not a free parameter")
print("    tuned to make the ratio 1.109.")
print("=" * 100)
print(f"    ||H_toric||_op = {np.linalg.norm(Ht,2):.6f}   (so ||V||=||H|| => p_eff = 8*eps)  ")
c2 = 8.4345e-03
for scale, lbl in ((1.0, "p_eff = eps        (||V||=1)"), (8.0, "p_eff = 8*eps      (lane's choice)"),
                   (4.0, "p_eff = 4*eps      (alt)")):
    pred = c2 * (scale * 1e-6) ** 2
    print(f"    {lbl:34s} predicts {pred:.4e}   vs W-61's 4.867e-13   ratio {pred/4.867e-13:8.3f}")
print("""    Only the 8x identification lands within a factor 2.  It is a ONE-PARAMETER choice justified
    by an external convention in another lane, and it is the ONLY choice that makes B3 pass.  B3 is
    therefore a consistency check with one adjustable input, not an independent reproduction.""")

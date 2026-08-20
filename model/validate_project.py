"""EVERY SEALED HEADLINE NUMBER, REPRODUCED THROUGH THE ONE ENTRY POINT."""
import numpy as np, sys
import grounded as G
from project_model import RecordSurface, ProjectModel
I2 = np.eye(2); X = np.array([[0, 1], [1, 0]], dtype=complex); Z = np.array([[1, 0], [0, -1]], dtype=complex)
def word(n, s):
    M = np.array([[1]], dtype=complex)
    for c in s: M = np.kron(M, {'I': I2, 'X': X, 'Z': Z}[c])
    return M
eV = 1.602176634e-19
M = ProjectModel(); n_pass = 0; n_fail = 0
def check(name, cond, detail=""):
    global n_pass, n_fail
    if cond: n_pass += 1; print(f"  PASS  {name}  {detail}")
    else:    n_fail += 1; print(f"  FAIL  {name}  {detail}")

print("VALIDATE THE FULL-PROJECT MODEL")
print("=" * 78)
# 1. T-29: the grain's lifetime, sealed 5.6279e16 s
grain = RecordSurface("CoCrPt grain", "magnetic anisotropy", 3.0 * G.KB * 300, 2.0e5 * 1.26e-24, 300.0, 1e9)
tau = M.lifetime(grain)
check("T-29 grain tau", abs(tau - 5.6279e16) / 5.6279e16 < 1e-3, f"tau = {tau:.4e} s (sealed 5.6279e+16)")
# 2. T-29: the five amended clauses pass with the sealed margins
c = M.clauses(grain, t_m=10 * 3.156e7, C_v=1e-3, W_write=5e-19)
check("T-29 clauses (ii')-(v') all pass",
      c['ii']['durable'] and c['iii']['passes'] and c['iv']['passes'] and c['v']['passes'],
      f"margins: ii {1/(10*3.156e7)/c['ii']['rate']:.1e}x, iii {c['iii']['margin']:.1e}x, "
      f"iv {c['iv']['over_landauer']:.0f}x Landauer, v E_b/kT = {c['v']['E_b_over_kT']:.1f}")
# 3. T-33: both laws on the six mechanisms, machine precision
SURFACES = [("CoCrPt HDD grain", "magnetic anisotropy", 3.0 * G.KB * 300, 60.8 * G.KB * 300, 300.0, 1e9),
            ("NAND floating gate", "trapped charge", 0.30 * eV, 1.60 * eV, 300.0, 1e13),
            ("DNA base tautomer", "chemical bond", 0.35 * eV, 1.30 * eV, 310.0, 1e13),
            ("Fe(II) spin crossover", "spin transition", 2.10 * G.KB * 300, 0.85 * eV, 300.0, 1e12),
            ("Azobenzene cis/trans", "photoisomerisation", 0.60 * eV, 1.05 * eV, 300.0, 1e13),
            ("Alanine enantiomer", "parity violation", 1.0e-13 * eV, 1.40 * eV, 300.0, 1e13)]
worst1 = worst2 = 0.0
for nm, mech, dE, Eb, T, f0 in SURFACES:
    s = RecordSurface(nm, mech, dE, Eb, T, f0)
    v = M.steady_value(s); pred = np.tanh(dE / (2 * G.KB * T))
    worst1 = max(worst1, abs(v - pred))
    kT = G.KB * T
    gd = f0 * np.exp(-(Eb + dE / 2) / kT); gu = f0 * np.exp(-(Eb - dE / 2) / kT)
    worst2 = max(worst2, abs(M.lifetime(s) - 1 / (gu + gd)) / (1 / (gu + gd)))
check("T-33 law 1 on six mechanisms", worst1 < 2e-15, f"worst abs err {worst1:.2e}")
check("T-33 law 2 on six mechanisms", worst2 < 1e-9, f"worst rel err {worst2:.2e}")
# 4. T-33: the model DECLINES on the controls
zirc = RecordSurface("Zircon U-238", "nuclear decay", 4.27e6 * eV, 4.27e6 * eV, 300.0, 1e21)
cmb = RecordSurface("CMB photon", "free flight", 0.0, 0.0, 2.7, 0.0, thermal=False)
check("T-33 declines on zircon", M.lifetime(zirc) is None, "rates unrepresentable -> None")
check("T-33 declines on the CMB photon", M.lifetime(cmb) is None, "no bath -> None")
# 5. T-28: the DEF-A corner, sealed dims 96 / 1536 / 24
for name, H, want in (("[[4,2,2]]", -(word(4, 'XXXX') + word(4, 'ZZZZ')), 96),
                      ("[[6,4,2]]", -(word(6, 'X' * 6) + word(6, 'Z' * 6)), 1536),
                      ("3q Ising", -(word(3, 'ZZI') + word(3, 'IZZ')), 24)):
    r = M.corner(H)
    check(f"T-28 corner {name}", r['slow_dim'] == r['commutant_dim'] == want,
          f"slow {r['slow_dim']} = commutant {r['commutant_dim']} = sealed {want}")
# 6. T-34/T-32 machinery: written vs unwritten ratios
cfg_w = M.configuration(1.602e-17, np.ones(1000))
cfg_u = M.configuration(1.602e-17, np.zeros(1000))
check("formation: written ratio = 1", cfg_w['ratio'] == 1.0, f"ratio {cfg_w['ratio']}")
check("formation: unwritten is null", cfg_u['ratio'] is None and cfg_u['sum'] == 0.0, "sum 0, ratio undefined")
print("=" * 78)
print(f"  {n_pass} PASS, {n_fail} FAIL")
sys.exit(1 if n_fail else 0)

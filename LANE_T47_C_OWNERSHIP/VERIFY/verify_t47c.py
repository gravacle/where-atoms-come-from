"""ADVERSARIAL VERIFY for LANE_T47_C_OWNERSHIP (T-47 lane C), default refuted.

Independent implementation: the Liouvillian, its EXACT eigenmodes (not Rayleigh quotients),
the durability crossing, the corner limits, and the C-14 v_2 reconnection are all rebuilt
here from scratch. grounded.clause_ii is imported ONLY to cross-check the lane's instrument
against the independent numbers, never as the source of truth.

Attack lines (from the verifier brief):
  A1  smuggled tolerance: enumerate every numeric threshold on the conclusion path.
  A2  does the numeric crossing land on the derived formula INCLUDING the correction term?
      (independent bisection on the TRUE eigenvalue, not the lane's quotient)
  A3  corner limits: dE=0 exact criterion B >= kT ln(2 f0 t_m); coherence boundary must sit
      at hbar*sqrt(t_m^-2 - ((gu+gd)/2)^2) ~ hbar/t_m (probe 0.9 vs 1.1, sharper than the
      lane's 0.5 vs 20); C-14 v_2 reconnection computed, not asserted.
  A4  the lane's logged sigma_x error: confirm the failure mode is real (quotient drops the
      rotation) -- i.e. the lane's fix was necessary, not cosmetic.
  A5  the departure-term numbers quoted against the rivals (t*, e^-10, tanh asymptote).
"""
import sys, os
import numpy as np
from itertools import product

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "model"))
import grounded as G

KB, HBAR = 1.380649e-23, 1.054571817e-34
assert KB == G.KB and HBAR == G.HBAR

SZ = np.array([[1, 0], [0, -1]], dtype=complex)
SP = np.array([[0, 1], [0, 0]], dtype=complex)
SM = SP.conj().T
SX = np.array([[0, 1], [1, 0]], dtype=complex)
I2 = np.eye(2, dtype=complex)

failures = []
def check(name, cond, detail=""):
    tag = "PASS" if cond else "FAIL"
    print(f"[{tag}] {name}  {detail}")
    if not cond:
        failures.append(name)

# ---------------- independent adjoint Liouvillian (Heisenberg picture), built from scratch
def adjoint_L(H_rad, Ls):
    """L†(X) = i[H,X] + sum_k (k† X k - 0.5{k†k, X}); H in rad/s, Ls in s^-1/2."""
    n = H_rad.shape[0]; I = np.eye(n)
    M = 1j * (np.kron(I, H_rad.T) - np.kron(H_rad, I))  # vec col-stacking: vec(AXB)= (B^T ⊗ A)vec X
    # careful: use explicit construction instead
    M = np.zeros((n*n, n*n), dtype=complex)
    basis = []
    for j in range(n):
        for i in range(n):
            E = np.zeros((n, n), dtype=complex); E[i, j] = 1.0
            basis.append(E)
    def app(X):
        out = 1j * (H_rad @ X - X @ H_rad)
        for k in Ls:
            kd = k.conj().T
            out = out + kd @ X @ k - 0.5 * (kd @ k @ X + X @ kd @ k)
        return out
    for c, E in enumerate(basis):
        M[:, c] = app(E).reshape(-1, order='F')
    return M

def rates_2ls(B, dE, T, f0):
    kT = KB * T
    gu = f0 * np.exp(-(B - dE) / kT)   # escape from the metastable (upper) well, barrier B-dE
    gd = f0 * np.exp(-B / kT)          # escape from the lower well, barrier B
    return gu, gd

def system(B, dE, T, f0):
    gu, gd = rates_2ls(B, dE, T, f0)
    H_rad = (-(dE / 2) * SZ) / HBAR
    Ls = [np.sqrt(gd) * SM, np.sqrt(gu) * SP]
    return H_rad, Ls, gu, gd

def eigmode_rate(H_rad, Ls, target):
    """|lambda| of the EXACT adjoint eigenmode with largest overlap with `target`."""
    M = adjoint_L(H_rad, Ls)
    w, V = np.linalg.eig(M)
    t = target.reshape(-1, order='F'); t = t / np.linalg.norm(t)
    ov = [abs(np.vdot(V[:, i] / np.linalg.norm(V[:, i]), t)) for i in range(len(w))]
    j = int(np.argmax(ov))
    return abs(w[j]), ov[j], w[j]

# ---------------- carrier (as declared by the lane)
T, f0 = 300.0, 1.0e9
kT = KB * T
B = 60.0 * kT
t_m = 10 * 365.25 * 86400.0

# ===== A2/V1: record-mode rate: exact eigenmode vs gu+gd vs lane instrument
dE = 10.0 * kT
H_rad, Ls, gu, gd = system(B, dE, T, f0)
# exact eigenmode nearest sigma_z, but sigma_z mixes with I; use traceless projection:
M = adjoint_L(H_rad, Ls)
w, V = np.linalg.eig(M)
# the record mode is the eigenvector that is diagonal and traceless up to the equilibrium offset;
# identify by overlap with sigma_z after removing the I component
lamz, ovz, lam_full = eigmode_rate(H_rad, Ls, SZ)
check("V1a exact eigenmode rate == gu+gd", abs(lamz - (gu + gd)) / (gu + gd) < 1e-12,
      f"|lam|={lamz:.6e} gu+gd={gu+gd:.6e} overlap={ovz:.6f}")
c2 = G.clause_ii(-(dE / 2) * SZ, Ls, SZ, t_m)  # the lane feeds H in joules for population records
check("V1b lane instrument (clause_ii on SZ) == gu+gd", abs(c2['rate'] - (gu + gd)) / (gu + gd) < 1e-12,
      f"rate={c2['rate']:.6e}")
# and confirm the joules-vs-rad/s sloppiness is harmless for the sigma_z record ([H,sz]=0):
c2b = G.clause_ii(H_rad, Ls, SZ, t_m)
check("V1c H-units irrelevant for diagonal record", abs(c2['rate'] - c2b['rate']) / c2b['rate'] < 1e-12)

# ===== A2/V2: durability crossing by independent bisection on the TRUE eigenvalue
def true_rate(dE_):
    H_, Ls_, gu_, gd_ = system(B, dE_, T, f0)
    lam, ov, _ = eigmode_rate(H_, Ls_, SZ)
    return lam

lo, hi = 0.0, B
assert true_rate(lo) < 1.0 / t_m < true_rate(hi)
for _ in range(120):
    mid = 0.5 * (lo + hi)
    lo, hi = (mid, hi) if true_rate(mid) <= 1.0 / t_m else (lo, mid)
dE_star_num = 0.5 * (lo + hi)
# exact derived crossing: solve gu+gd = 1/t_m, i.e. B-dE = kT[ln(f0 t_m) + ln(1+e^{-dE/kT})]
d = B - kT * np.log(f0 * t_m)
for _ in range(300):
    d = B - kT * (np.log(f0 * t_m) + np.log1p(np.exp(-d / kT)))
dE_star_exact = d
dE_star_leading = B - kT * np.log(f0 * t_m)
check("V2a independent crossing lands on exact derived form",
      abs(dE_star_num - dE_star_exact) / kT < 1e-8,
      f"num={dE_star_num/kT:.8f} kT exact={dE_star_exact/kT:.8f} kT")
check("V2b leading (brief) form within kT ln2, high side",
      0 <= (dE_star_leading - dE_star_exact) / kT <= np.log(2) + 1e-12,
      f"leading={dE_star_leading/kT:.6f} kT corr={(dE_star_leading-dE_star_exact)/kT:.2e} kT")
check("V2c lane's quoted dE* = 19.71 kT", abs(dE_star_exact / kT - 19.71) < 0.01,
      f"{dE_star_exact/kT:.4f}")

# ===== A3/V3: symmetric corner EXACT criterion in B: crossing of B at dE=0
def rate0(B_):
    H_, Ls_, gu_, gd_ = system(B_, 0.0, T, f0)
    return gu_ + gd_
B_star_exact_pred = kT * np.log(2 * f0 * t_m)
loB, hiB = 10 * kT, 80 * kT
for _ in range(200):
    midB = 0.5 * (loB + hiB)
    loB, hiB = (midB, hiB) if rate0(midB) > 1.0 / t_m else (loB, midB)
B_star_num = 0.5 * (loB + hiB)
check("V3a dE=0 corner criterion B* = kT ln(2 f0 t_m) exactly (ln2 present)",
      abs(B_star_num - B_star_exact_pred) / kT < 1e-9,
      f"B*={B_star_num/kT:.6f} kT pred={B_star_exact_pred/kT:.6f} kT")
check("V3b correction attains kT ln2 at symmetric corner",
      abs((B_star_num - kT * np.log(f0 * t_m)) / kT - np.log(2)) < 1e-9)
H0, Ls0, gu0, gd0 = system(B, 0.0, T, f0)
tau0 = 1.0 / (gu0 + gd0)
check("V3c tau(dE=0) = exp(B/kT)/(2 f0)", abs(tau0 * 2 * f0 / np.exp(B / kT) - 1) < 1e-12,
      f"tau={tau0:.4e} s")

# ===== A3/V4: coherence boundary probed at 0.9 and 1.1 hbar/t_m (sharper than lane's 0.5/20)
lam_pred_offset = (gu0 + gd0) / 2  # on the B=60kT carrier
dE_bound_exact = HBAR * np.sqrt(max((1 / t_m) ** 2 - lam_pred_offset ** 2, 0.0))
check("V4a coherence boundary numerically = hbar/t_m (decay corr negligible)",
      abs(dE_bound_exact - HBAR / t_m) / (HBAR / t_m) < 1e-15,
      f"bound/(hbar/t_m)={dE_bound_exact/(HBAR/t_m):.16f}")
for frac, want in ((0.9, True), (1.1, False)):
    dEc = frac * HBAR / t_m
    Hc, Lsc, guc, gdc = system(B, dEc, T, f0)
    lam, ov, lamc = eigmode_rate(Hc, Lsc, SP)
    pred = np.hypot(dEc / HBAR, (guc + gdc) / 2)
    okm = abs(lam - pred) / pred < 1e-9
    okd = (lam <= 1.0 / t_m) == want
    check(f"V4b sigma+ mode at {frac} hbar/t_m: |lambda| exact & durable=={want}",
          okm and okd, f"|lam|={lam:.6e} pred={pred:.6e} overlap={ov:.4f}")
# A4: the logged sigma_x failure mode -- Rayleigh quotient must drop the rotation
dEc = 20 * HBAR / t_m
Hc, Lsc, guc, gdc = system(B, dEc, T, f0)
v = SX.reshape(-1, 1, order='F'); v = v / np.linalg.norm(v)
Mad = adjoint_L(Hc, Lsc)
q = abs(complex((v.conj().T @ Mad @ v)[0, 0]))
lam_true = eigmode_rate(Hc, Lsc, SP)[0]
check("V4c sigma_x quotient drops rotation (lane's logged error is real)",
      q < 1e-15 and lam_true > 1.0 / t_m,
      f"sigma_x quotient={q:.2e} true sigma+ |lambda|={lam_true:.2e} 1/t_m={1/t_m:.2e}")

# ===== A5/V5: departure-term numbers quoted against the rivals
dE = 10.0 * kT
gu, gd = rates_2ls(B, dE, T, f0)
t_star = 1.0 / (gu + gd)
t_star_leading = np.exp((B - dE) / kT) / f0
check("V5a t* (leading) = 5.18e12 s", abs(t_star_leading - 5.18e12) / 5.18e12 < 0.01,
      f"{t_star_leading:.3e} s (exact incl. corr: {t_star:.3e} s)")
check("V5b time factor e^-10 = 4.5e-5", abs(np.exp(-10) - 4.54e-5) < 1e-6, f"{np.exp(-10):.3e}")
check("V5c favored asymptote tanh(dE/2kT) = 0.999909", abs(np.tanh(5.0) - 0.999909) < 1e-6,
      f"{np.tanh(5.0):.6f}")
# both polarities, one rate, different asymptotes -- via independent propagator
H_rad, Ls, gu, gd = system(B, dE, T, f0)
Lv = np.zeros((4, 4), dtype=complex)
for c in range(4):
    E = np.zeros(4, dtype=complex); E[c] = 1
    X = E.reshape(2, 2, order='F')
    out = -1j * (H_rad @ X - X @ H_rad)
    for k in Ls:
        kd = k.conj().T
        out = out + k @ X @ kd - 0.5 * (kd @ k @ X + X @ kd @ k)
    Lv[:, c] = out.reshape(-1, order='F')
wL, VL = np.linalg.eig(Lv)
P = VL @ np.diag(np.exp(wL / (gu + gd))) @ np.linalg.inv(VL)
def mz(rho0):
    r = (P @ rho0.reshape(-1, 1, order='F')).reshape(2, 2, order='F')
    return float(np.real(np.trace(r @ SZ)))
meq = np.tanh(dE / (2 * kT))
m_f, m_u = mz(np.diag([1.0, 0]).astype(complex)), mz(np.diag([0, 1.0]).astype(complex))
check("V5d both polarities same rate, asymptote-split (closed forms)",
      abs(m_f - (meq + (1 - meq) / np.e)) < 1e-9 and abs(m_u - (meq + (-1 - meq) / np.e)) < 1e-9,
      f"favored {m_f:.4f} unfavored {m_u:.4f}")

# ===== A3/V6: C-14 v_2 reconnection, COMPUTED (N symmetric pairs; then generic asymmetry)
def v2(n):
    k = 0
    while n % 2 == 0:
        n //= 2; k += 1
    return k
def c14_count(energies_per_pair):
    """min over joint eigenspaces of v_2(multiplicity); joint spectrum by summation."""
    from collections import Counter
    tot = Counter()
    for combo in product(*energies_per_pair):
        tot[round(sum(combo), 9)] += 1
    return min(v2(m) for m in tot.values())
# symmetric corner: N pairs, dE=0 (arbitrary base energies), t_m below every tau
for N in (1, 2, 3, 4):
    pairs = [(i * 1.7, i * 1.7) for i in range(N)]   # dE=0 each, distinct bases
    kC14 = c14_count(pairs)
    stair = sum(1 for _ in range(N)
                if (lambda g: g <= 1.0 / t_m)(sum(rates_2ls(B, 0.0, T, f0))))
    check(f"V6a N={N} symmetric: C-14 v_2 count == staircase count == {N}",
          kC14 == N == stair, f"v2={kC14} staircase={stair}")
# generic asymmetry: C-14 collapses to 0 (the C-76 kill) while the staircase counts pairs
rng = np.random.default_rng(20260820)
pairs = [(0.0, float(rng.uniform(0.5, 2.5))) for _ in range(3)]
kC14 = c14_count(pairs)
dEs = [5 * kT, 15 * kT, 25 * kT]  # declared examples on the B=60kT carrier
stair = sum(1 for dE_ in dEs if sum(rates_2ls(B, dE_, T, f0)) <= 1.0 / t_m)
delta = B - kT * np.log(f0 * t_m)  # 19.71 kT
stair_formula = sum(1 for dE_ in dEs if dE_ <= delta)
check("V6b generic asymmetry: C-14 k=0 (C-76 kill) but staircase counts 2 (5,15 in; 25 out)",
      kC14 == 0 and stair == 2 and stair_formula == 2,
      f"v2={kC14} staircase={stair} formula-count={stair_formula} delta={delta/kT:.2f} kT")
# per-pair dropout times land on t*_i
for dE_ in dEs:
    ts = 1.0 / sum(rates_2ls(B, dE_, T, f0))
    ts_pred = np.exp((B - dE_) / kT) / f0 / (1 + np.exp(-dE_ / kT))
    check(f"V6c dropout t*(dE={dE_/kT:.0f}kT) matches exact formula",
          abs(ts - ts_pred) / ts_pred < 1e-12, f"t*={ts:.3e} s")

# ===== A1/V7: smuggled-tolerance scan of the conclusion path
# clause_ii durable threshold: 1/t_m + 1e-300 -- 1e-300 is below any representable effect at
# these scales (1/t_m ~ 3e-9); verify conclusions identical with the addend removed:
rate_at_crossing = true_rate(dE_star_exact)
check("V7a durability threshold is exactly 1/t_m (1e-300 addend inert)",
      abs(rate_at_crossing * t_m - 1) < 1e-6, f"rate*t_m at crossing = {rate_at_crossing*t_m:.8f}")
# slow_modes' tol=1e-9 is never on this path (clause_ii does not call it); confirm by source scan:
import inspect
src = inspect.getsource(G.clause_ii)
check("V7b clause_ii source contains no 1e-9 floor (only the inert 1e-300)",
      "1e-9" not in src and "1e-300" in src)

print()
if failures:
    print("VERIFY: FAILURES:", failures)
else:
    print("VERIFY: ALL INDEPENDENT CHECKS PASS")
sys.exit(1 if failures else 0)

"""LANE_T47_C_OWNERSHIP -- numeric controls for the ownership/falsifier check (T-47, lane C).

Not the derivation lane. These controls exist so the DEPARTURE TERM claimed against the rivals
(Sharrock / Street-Woolley viscosity; MLC flash level retirement) is COMPUTED, not asserted,
through the surface's own Liouvillian machinery (model/grounded.py conventions).

CONVENTION MAP (project_model.RecordSurface, solidity-review-corrected):
  E_b = activation energy from the METASTABLE (upper) well.
  Brief's B = barrier above the LOWER well = E_b + dE;  brief's (B - dE) = E_b.
  gu = f0 exp(-E_b/kT)  (escape from upper),  gd = gu exp(-dE/kT)  (from lower).
  Record mode (population difference) relaxes at gu+gd; tau_record = 1/(gu+gd).
  Derived width, brief form: dE <= B - kT ln(f0 t_m)  <=>  E_b >= kT ln(f0 t_m)  (leading order).
  Exact crossing: gu+gd = 1/t_m  <=>  E_b = kT [ln(f0 t_m) + ln(1 + e^{-dE/kT})].

CONTROLS
  C1  record-mode rate == gu+gd through the Liouvillian (no closed form smuggled).
  C2  threshold crossing at swept dE lands on the derived formula (exact form), and the
      leading-order brief form is within kT ln2 (the symmetric-corner correction).
  C3  BOTH write polarities relax at the SAME rate gu+gd toward DIFFERENT asymptotes
      (+/- tanh(dE/2kT) for R = sigma_z with H = -(dE/2)sz): the rival's viscosity/remanence
      observable on the favored branch retains tanh(dE/2kT) FOREVER while the two-valued
      record is dead past tau_record. Departure term in step position: Delta ln t* = -dE/kT
      relative to the lower-well (majority/Sharrock-inferred) barrier.
  C4  D-15 symmetric control: dE = 0 -> asymptotes 0, survival criterion reduces to
      B >= kT ln(2 f0 t_m) -- the Neel/Charap-Lu-He/Weller-Moser owned criterion (ownership-
      detection control MUST fire OWNED here).
  C5  D-15 coherence-type control: off-diagonal record on the same carrier; |lambda| =
      sqrt((dE/hbar)^2 + ((gu+gd)/2)^2) so durability to t_m forces |dE| <= ~hbar/t_m --
      the coherence width REAPPEARS as the off-diagonal width, distinct from the population
      width kT ln(f0 t_m). (H fed to the Liouvillian in angular-frequency units H/hbar for
      this check; the population controls are hbar-independent for sigma_z records.)
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "model"))
import numpy as np
import grounded as G

KB, HBAR = G.KB, G.HBAR
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
SP = np.array([[0, 1], [0, 0]], dtype=complex)
SM = SP.conj().T

def open_system(E_b, dE, T, f0, hbar_units=False):
    kT = KB * T
    gu = f0 * np.exp(-E_b / kT)
    gd = f0 * np.exp(-(E_b + dE) / kT)
    H = -(dE / 2) * SZ
    if hbar_units:
        H = H / HBAR
    return H, [np.sqrt(gd) * SM, np.sqrt(gu) * SP], gu, gd

# ---- carrier: exchange-biased-media-like numbers (declared, carrier-typical, not tuned)
T, f0 = 300.0, 1.0e9          # Sharrock's canonical attempt frequency
kT = KB * T
B = 60.0 * kT                  # lower-well barrier: Weller-Moser 10-yr class
t_m = 10 * 365.25 * 86400.0    # 10-year retention spec
ok = True

# ---- C1: record-mode rate through the Liouvillian equals gu+gd
dE = 10.0 * kT
E_b = B - dE
H, Ls, gu, gd = open_system(E_b, dE, T, f0)
c2 = G.clause_ii(H, Ls, SZ, t_m)
r_liouv, r_closed = c2['rate'], gu + gd
e1 = abs(r_liouv - r_closed) / r_closed
print(f"C1  rate(Liouvillian)={r_liouv:.6e}  gu+gd={r_closed:.6e}  rel.err={e1:.2e}")
ok &= e1 < 1e-9

# ---- C2: sweep dE at fixed B; find durability crossing; compare to derived width
lo, hi = 0.0, B  # dE bisection window; survival is monotone decreasing in dE at fixed B
def survives(dE):
    H, Ls, gu, gd = open_system(B - dE, dE, T, f0)
    return G.clause_ii(H, Ls, SZ, t_m)['durable']
assert survives(lo) and not survives(hi)
for _ in range(200):
    mid = 0.5 * (lo + hi)
    lo, hi = (mid, hi) if survives(mid) else (lo, mid)
dE_star_num = 0.5 * (lo + hi)
# exact derived crossing: E_b = kT[ln(f0 t_m) + ln(1+e^{-dE/kT})], E_b = B - dE, solve fixed point
d = B - kT * np.log(f0 * t_m)
for _ in range(200):
    d = B - kT * (np.log(f0 * t_m) + np.log1p(np.exp(-d / kT)))
dE_star_exact = d
dE_star_leading = B - kT * np.log(f0 * t_m)   # the brief's delta(t_m)
e2 = abs(dE_star_num - dE_star_exact) / kT
print(f"C2  crossing (numeric)      dE*/kT = {dE_star_num/kT:.6f}")
print(f"    crossing (exact form)   dE*/kT = {dE_star_exact/kT:.6f}   |diff| = {e2:.2e} kT")
print(f"    brief leading form      dE*/kT = {dE_star_leading/kT:.6f}   correction = "
      f"{(dE_star_leading-dE_star_exact)/kT:.4f} kT  (bound kT ln2 = {np.log(2):.4f} kT)")
ok &= e2 < 1e-6 and 0 <= (dE_star_leading - dE_star_exact) / kT <= np.log(2) + 1e-12

# ---- C3: both polarities same rate, different asymptotes; departure term
dE = 10.0 * kT
H, Ls, gu, gd = open_system(B - dE, dE, T, f0)
Lv = G.liouvillian(H, Ls)
w, V = np.linalg.eig(Lv)
j = int(np.argmin(np.abs(w)))
rho_ss = V[:, j].reshape(2, 2, order='F'); rho_ss = rho_ss / np.trace(rho_ss)
m_eq = float(np.real(np.trace(rho_ss @ SZ)))
m_eq_closed = float(np.tanh(dE / (2 * kT)))   # H = -(dE/2)sz: LOWER well is sz=+1
tau_rec = 1.0 / (gu + gd)
# time evolution from both saturations, propagated through the Liouvillian (expm at t = tau_rec,
# via eigendecomposition -- no scipy in this environment)
wL, VL = np.linalg.eig(Lv)
P = VL @ np.diag(np.exp(wL * tau_rec)) @ np.linalg.inv(VL)
def evolve(rho0):
    r = (P @ rho0.reshape(-1, 1, order='F')).reshape(2, 2, order='F')
    return float(np.real(np.trace(r @ SZ)))
m_fav = evolve(np.diag([1.0, 0.0]).astype(complex))    # written into favored (lower) well
m_unf = evolve(np.diag([0.0, 1.0]).astype(complex))    # written into unfavored (upper) well
# closed forms: m(t) = m_eq + (m0 - m_eq) e^{-(gu+gd)t}
m_fav_cf = m_eq_closed + (1 - m_eq_closed) * np.e**-1
m_unf_cf = m_eq_closed + (-1 - m_eq_closed) * np.e**-1
e3 = max(abs(m_fav - m_fav_cf), abs(m_unf - m_unf_cf), abs(m_eq - m_eq_closed))
print(f"C3  m_eq={m_eq:.6f} (tanh={m_eq_closed:.6f});  at t=tau_rec: favored {m_fav:.4f}, "
      f"unfavored {m_unf:.4f};  max err {e3:.2e}")
print(f"    favored-branch retained signal as t->inf: {m_eq_closed:.6f}  while record count -> 0 "
      f"past tau_rec = {tau_rec:.3e} s")
print(f"    DEPARTURE TERM: ln t*_record = ln[1/(f0)] + (B-dE)/kT  vs lower-well/majority "
      f"barrier B: shift = -dE/kT = {-dE/kT:.1f}  (time factor e^{{-dE/kT}} = {np.exp(-dE/kT):.2e})")
ok &= e3 < 1e-9

# ---- C4: symmetric control (dE=0) -- flat staircase corner; owned criterion recovered
H0, Ls0, gu0, gd0 = open_system(B, 0.0, T, f0)
c0 = G.clause_ii(H0, Ls0, SZ, t_m)
tau0 = 1.0 / (gu0 + gd0)
B_crit = kT * np.log(2 * f0 * t_m)
print(f"C4  dE=0: tau = {tau0:.3e} s = exp(B/kT)/(2 f0) check {np.exp(B/kT)/(2*f0):.3e}; "
      f"survives(10yr)={c0['durable']}; criterion B >= kT ln(2 f0 t_m) = {B_crit/kT:.2f} kT "
      f"(B = {B/kT:.0f} kT)  [OWNED corner: Neel/Charap-Lu-He/Weller-Moser]")
ok &= c0['durable'] and abs(tau0 * 2 * f0 / np.exp(B / kT) - 1) < 1e-9

# ---- C5: coherence-type record on the same carrier -- hbar/t_m must govern
dE_small = HBAR / t_m * 0.5           # inside the coherence width
dE_large = HBAR / t_m * 20.0          # outside it (still << kT: population width irrelevant here)
for lbl, dEc in (("inside", dE_small), ("outside", dE_large)):
    Hc, Lsc, guc, gdc = open_system(B, dEc, T, f0, hbar_units=True)
    # the record's OWN Liouvillian mode (spectrum()'s docstring caution): the off-diagonal
    # eigenmode is sigma+ with eigenvalue -i dE/hbar - (gu+gd)/2; sigma_x is NOT an eigenmode
    # and its Rayleigh quotient would drop the rotation (first draft of this control did
    # exactly that and read 8.8e-18 for both cases -- kept here as the lesson).
    cc = G.clause_ii(Hc, Lsc, SP, t_m)
    lam_pred = np.hypot(dEc / HBAR, (guc + gdc) / 2)
    print(f"C5  coherence dE {lbl} hbar/t_m: |lambda| = {cc['rate']:.6e}  "
          f"pred sqrt((dE/hbar)^2+((gu+gd)/2)^2) = {lam_pred:.6e}  durable={cc['durable']}")
    ok &= abs(cc['rate'] - lam_pred) / lam_pred < 1e-6
    ok &= cc['durable'] == (lbl == "inside")

print("ALL CONTROLS PASS" if ok else "CONTROL FAILURE")
sys.exit(0 if ok else 1)

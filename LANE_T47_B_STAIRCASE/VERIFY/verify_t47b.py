"""ADVERSARIAL VERIFY of LANE_T47_B_STAIRCASE (T-47 creator step).

Independent re-derivation and re-computation. All instruments in this file are written from
scratch (own GKSL vectorisation, own mode selection by decay-fit AND by eigen-expansion, own
Davies-Metropolis build, own basin lumping). grounded.py / project_model.py are imported ONLY
where the lane's claim is ABOUT those instruments (caveat 3a/3b), i.e. to test the claim.

Attack list:
  V1  the width: re-derive delta(t_m) from clause (ii') symbolically; scan the lane script for
      any chosen tolerance entering the COUNT (the C-76 kill).
  V2  numeric crossing: independent bisection on an EXPM DECAY-FIT instrument (not eigenvalues)
      must land on delta_exact including the -1 correction term; candidate offsets re-derived.
  V3  coherence control: independent bisection of the coherence-mode crossing vs hbar/t_m.
  V4  ensembles: recompute tau_i (own eig instrument), t*_exact, the 1201-point staircase both
      ways, the candidate's 82-point discrepancy, on BOTH declared ensembles.
  V5  symmetric control: threshold = e^(B/kT)/(2 f0); flatness.
  V6  corner limit: reproduce the common-B clustered-v2 table rows (agree AND disagree rows) and
      the chaining counterexample with own clustering code.
  V7  T-31 carrier: full independent rebuild (H, Davies-Metropolis, 256-dim Liouvillian, basin
      lumping); reproduce dE_i, tau_i, staircase [2,1,0], v2-exact = 0, generator classicality.
  V8  departure term: independent closed-form m(t) (rate-equation solution, no Liouvillian) vs
      the lane's exact-propagator numbers; departure sums.
  V9  instrument-scope claims (caveat 3a/3b): test grounded.clause_ii on the 16-dim carrier and
      the H-units convention consequence.
"""
import sys, os
import numpy as np

def expm(M):
    """Own matrix exponential: scaling-and-squaring on the Taylor series (no scipy here)."""
    M = np.asarray(M, complex)
    nrm = np.linalg.norm(M, 1)
    s = max(0, int(np.ceil(np.log2(max(nrm, 1e-300)))) + 4) if nrm > 0 else 0
    A = M / (2 ** s)
    E = np.eye(M.shape[0], dtype=complex); term = np.eye(M.shape[0], dtype=complex)
    for k in range(1, 30):
        term = term @ A / k
        E = E + term
    for _ in range(s):
        E = E @ E
    return E

KB = 1.380649e-23
HBAR = 1.054571817e-34
EV = 1.602176634e-19

def say(*a):
    print(*a); sys.stdout.flush()

FAILS = []
def vcheck(name, ok, detail=""):
    if not ok:
        FAILS.append(name)
    say(f"  [{'OK' if ok else 'REFUTE'}] {name}" + (f"  ({detail})" if detail else ""))

# ---------------------------------------------------------------- own instruments
def gksl(H, Ls):
    """Column-stacked GKSL generator, written independently."""
    H = np.asarray(H, complex); n = H.shape[0]; I = np.eye(n)
    L = -1j * (np.kron(I, H) - np.kron(H.T, I))
    for A in Ls:
        A = np.asarray(A, complex); AdA = A.conj().T @ A
        L += np.kron(A.conj(), A) - 0.5 * (np.kron(I, AdA) + np.kron(AdA.T, I))
    return L

SZ = np.array([[1, 0], [0, -1]], complex)
SX = np.array([[0, 1], [1, 0]], complex)
SP = np.array([[0, 1], [0, 0]], complex)   # |lower><upper| ... convention as in project_model
SM = SP.conj().T

def two_level(B, dE, T, f0):
    """Own build of the two-state surface with the SAME declared convention:
       B = barrier above the LOWER well; escape from upper = f0 e^-((B-dE)/kT),
       escape from lower = f0 e^-(B/kT). H in ANGULAR FREQUENCY units."""
    kT = KB * T
    gu = f0 * np.exp(-(B - dE) / kT)   # from metastable (upper) well
    gl = f0 * np.exp(-B / kT)          # from lower well
    H = -(dE / 2) * SZ / HBAR
    return H, [np.sqrt(gl) * SM, np.sqrt(gu) * SP], gu, gl

def mode_rate_expansion(H, Ls, R):
    """Own eigen-expansion mode selector (independent code, same physics)."""
    Lad = gksl(H, Ls).conj().T
    w, U = np.linalg.eig(Lad)
    v = np.asarray(R, complex).reshape(-1, 1, order='F'); v = v / np.linalg.norm(v)
    c = np.linalg.solve(U, v).ravel()
    wt = np.abs(c) * np.linalg.norm(U, axis=0)
    wt[int(np.argmin(np.abs(w)))] = -1.0
    return w[int(np.argmax(wt))]

def tau_pop_decayfit(B, dE, T, f0, t_probe):
    """FULLY independent instrument: exact propagator at two times, fit the exponential of
       m(t) - m_inf. No eigen-mode selection at all."""
    H, Ls, gu, gl = two_level(B, dE, T, f0)
    L = gksl(H, Ls)
    rho0 = np.diag([0.0, 1.0]).astype(complex)   # written into the UPPER (metastable) well
    v0 = rho0.reshape(-1, 1, order='F')
    wL, UL = np.linalg.eig(L)
    c0 = np.linalg.solve(UL, v0)
    def m_at(t):
        rt = (UL @ (np.exp(wL * t) * c0.ravel())[:, None]).reshape(2, 2, order='F')
        return float(np.real(np.trace(rt @ SZ)))
    # three equally spaced probes eliminate m_eq entirely:
    # (m1-m2)/(m2-m3) = e^(dt/tau) for a single exponential. (A two-point fit against the
    # ANALYTIC tanh equilibrium failed in verify: once the decay finishes to within 1e-16 the
    # residual made the ratio ~1 and read 'immortal' — my instrument's bug, logged here.)
    dt = 0.5 * t_probe
    m1, m2, m3 = m_at(dt), m_at(2 * dt), m_at(3 * dt)
    a, b = m1 - m2, m2 - m3
    if not (np.isfinite(a) and np.isfinite(b)) or a == 0 or b == 0 or a * b < 0:
        return 0.0        # decay finished before the first probe: dead at this probe scale
    r = a / b
    if r <= 1.0:
        return np.inf     # no resolvable decay over the window: alive at this probe scale
    return dt / np.log(r)

def tau_pop_eig(B, dE, T, f0):
    H, Ls, gu, gl = two_level(B, dE, T, f0)
    return 1.0 / abs(mode_rate_expansion(H, Ls, SZ))

# ---------------------------------------------------------------- the law, re-derived
def delta_exact(B, T, f0, t_m):
    kT = KB * T
    arg = np.exp(B / kT) / (f0 * t_m) - 1.0
    return -np.inf if arg <= 0 else kT * np.log(arg)

def ln_tstar_exact(B, dE, T, f0):
    kT = KB * T
    return (B - dE) / kT - np.log(f0) - np.log1p(np.exp(-dE / kT))

def ln_tstar_cand(B, dE, T, f0):
    return (B - dE) / (KB * T) - np.log(f0)

say("=" * 100)
say("VERIFY T-47 LANE B — independent rerun")
say("=" * 100)

# ================================================================ V1: derivation audit
say("")
say("V1. THE WIDTH — algebra re-derived; tolerance scan")
# clause (ii') registered in grounded.py: durable iff |lambda| <= 1/t_m.  The record observable's
# decaying population mode has lambda = -(g_u+g_l) EXACTLY (checked below at machine precision):
B0, T0, f00 = 1.2 * EV, 350.0, 1e9
worst_sum = 0.0
for dE_ev in (0.0, 0.05, 0.2, 0.4):
    dE = dE_ev * EV
    H, Ls, gu, gl = two_level(B0, dE, T0, f00)
    lam = mode_rate_expansion(H, Ls, SZ)
    worst_sum = max(worst_sum, abs(abs(lam) - (gu + gl)) / (gu + gl))
vcheck("V1a: record-mode rate = g_u+g_l exactly (own instrument)", worst_sum <= 1e-9,
       f"worst rel {worst_sum:.1e}")
# g_u+g_l <= 1/t_m  <=>  f0 t_m e^(-B/kT)(e^(dE/kT)+1) <= 1  <=>  dE <= kT ln(e^(B/kT)/(f0 t_m)-1)
kT = KB * T0
for t_m in (1e2, 1e6):
    d = delta_exact(B0, T0, f00, t_m)
    lhs = f00 * np.exp(-B0 / kT) * (np.exp(d / kT) + 1.0)
    vcheck(f"V1b: delta({t_m:.0e}) sits exactly on g_u+g_l = 1/t_m", abs(lhs * t_m - 1) <= 1e-12,
           f"resid {abs(lhs*t_m-1):.1e}")
say("  tolerance scan of the lane script: grep for chosen widths / floors")
os.system('grep -n "1e-" "%s" | grep -v "^.*#" | head -30' %
          (os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "t47b_staircase.py")))

# ================================================================ V2: the crossing, own instrument
say("")
say("V2. NUMERIC CROSSING on a decay-fit instrument (no eigen-selection), incl. correction term")
lane_A = {1e2: 0.436077474379, 1e4: 0.297180902789, 1e6: 0.158128505564, 1e7: 0.087211155930}
lane_cand_off = {1e2: 0.000001, 1e4: 0.000076, 1e6: 0.007605, 1e7: 0.077913}
worst = 0.0
for t_m, d_lane in lane_A.items():
    d_mine = delta_exact(B0, T0, f00, t_m)
    lo, hi = 0.0, B0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        tau = tau_pop_decayfit(B0, mid, T0, f00, t_m)
        if tau >= t_m: lo = mid
        else: hi = mid
    dE_star = 0.5 * (lo + hi)
    rel = abs(dE_star - d_mine) / d_mine
    worst = max(worst, rel)
    cand = B0 - kT * np.log(f00 * t_m)
    off = (cand - d_mine) / (kT * np.log(2.0))
    say(f"    t_m={t_m:.0e}: delta_mine={d_mine/EV:.12f} eV (lane {d_lane:.12f}); "
        f"decay-fit crossing rel err {rel:.1e}; cand offset {off:.6f} (lane {lane_cand_off[t_m]:.6f})")
    vcheck(f"V2: table value reproduced at t_m={t_m:.0e}", abs(d_mine / EV - d_lane) <= 1e-9)
    vcheck(f"V2: cand offset reproduced at t_m={t_m:.0e}", abs(off - lane_cand_off[t_m]) <= 1e-4)
vcheck("V2: independent decay-fit crossing lands on delta_exact (<=1e-6 rel)", worst <= 1e-6,
       f"worst {worst:.1e}")
# does the crossing distinguish exact from candidate at t_m=1e7?  gap = 0.0779*kT*ln2
gap = lane_cand_off[1e7] * kT * np.log(2.0) / EV
say(f"    exact-vs-candidate gap at 1e7 s = {gap:.6f} eV vs crossing precision ~1e-8 eV: "
    f"the instrument RESOLVES the correction term" if worst * 0.087 < gap else "NOT RESOLVED")
vcheck("V2: correction term resolved by the instrument", worst * lane_A[1e7] < gap)

# ================================================================ V3: coherence control
say("")
say("V3. COHERENCE-TYPE RECORD — hbar/t_m from the same clause (own bisection)")
B1, T1, f1, t_m1 = 1.4 * EV, 350.0, 1e9, 3.156e7
H, Ls, gu, gl = two_level(B1, 0.0, T1, f1)
g_tot = gu + gl
lo, hi = 0.0, 1e-38
for _ in range(200):
    mid = 0.5 * (lo + hi)
    Hm, Lm, _, _ = two_level(B1, mid, T1, f1)
    lam = mode_rate_expansion(Hm, Lm, SX)
    if abs(lam) <= 1.0 / t_m1: lo = mid
    else: hi = mid
dE_coh = 0.5 * (lo + hi)
w_coh = HBAR / t_m1
w_ex = HBAR * np.sqrt(1.0 / t_m1**2 - (g_tot / 2) ** 2)
say(f"    g_tot={g_tot:.4e} (lane 1.386e-11); dE*={dE_coh:.9e} (lane 3.341482230e-42); "
    f"hbar/t_m={w_coh:.9e}")
vcheck("V3: coherence crossing = hbar sqrt(1/t^2-(g/2)^2), rel<=1e-9",
       abs(dE_coh - w_ex) / w_ex <= 1e-9, f"rel {abs(dE_coh-w_ex)/w_ex:.1e}")
vcheck("V3: lane's 2.39e-8 rel offset from hbar/t_m reproduced",
       abs(abs(dE_coh - w_coh) / w_coh - 2.39e-8) <= 1e-9)
d_pop = delta_exact(B1, T1, f1, t_m1)
vcheck("V3: ratio of the two widths ~1.22e22 reproduced",
       abs(d_pop / w_coh / 1.22e22 - 1) < 0.01, f"ratio {d_pop/w_coh:.3e}")

# ================================================================ V4: ensembles
say("")
say("V4. ENSEMBLES — k(t_m) two ways, own code, both declared grids")
def verify_ensemble(name, Bs, dEs, T, f0, lane_taus):
    N = len(Bs)
    taus = np.array([tau_pop_eig(Bs[i], dEs[i], T, f0) for i in range(N)])
    lt_ex = np.array([ln_tstar_exact(Bs[i], dEs[i], T, f0) for i in range(N)])
    lt_cd = np.array([ln_tstar_cand(Bs[i], dEs[i], T, f0) for i in range(N)])
    worst_step = float(np.max(np.abs(np.log(taus) - lt_ex)))
    rel_lane = float(np.max(np.abs(taus / lane_taus - 1)))
    tg = np.logspace(0, 12, 1201); lg = np.log(tg)
    k1 = (taus[None, :] >= tg[:, None]).sum(axis=1)
    k2 = (lt_ex[None, :] >= lg[:, None]).sum(axis=1)
    kc = (lt_cd[None, :] >= lg[:, None]).sum(axis=1)
    d12 = int(np.max(np.abs(k1 - k2))); ndc = int((np.abs(k1 - kc) > 0).sum())
    say(f"    {name}: max|k_direct-k_exact|={d12} over 1201; candidate differs at {ndc} pts "
        f"(lane: 0 and 82); step law worst |dln|={worst_step:.1e}; taus vs lane rel {rel_lane:.1e}")
    vcheck(f"V4[{name}]: two routes agree everywhere", d12 == 0)
    vcheck(f"V4[{name}]: candidate discrepancy count = 82", ndc == 82, f"got {ndc}")
    vcheck(f"V4[{name}]: steps land on t*_exact <=1e-10", worst_step <= 1e-10)
    vcheck(f"V4[{name}]: lane's tau table reproduced <=1e-6 rel", rel_lane <= 1e-6)
    vcheck(f"V4[{name}]: staircase monotone decreasing", bool(np.all(np.diff(k1) <= 0)))
Bs_u = np.array([0.95 + 0.40 * i / 11.0 for i in range(12)]) * EV
dEs_u = np.array([0.18 * ((7 * i + 3) % 12) / 11.0 for i in range(12)]) * EV
lane_tau_u = np.array([7.845972e+03, 6.995896e+02, 3.315836e+04, 8.896269e+05, 1.302688e+05,
                       5.009554e+06, 4.979281e+05, 2.265904e+07, 1.884988e+06, 9.156432e+07,
                       3.026220e+09, 3.536011e+08])
verify_ensemble("UNIFORM", Bs_u, dEs_u, 350.0, 1e9, lane_tau_u)
zs = np.array([-1.5, -1.1, -0.7, -0.35, 0.0, 0.0, 0.35, 0.7, 1.1, 1.5])
zp = np.array([0.7, -0.35, 1.5, 0.0, -1.5, 1.1, -0.7, 0.35, -1.1, 0.0])
Bs_l = 1.10 * np.exp(0.12 * zs) * EV
dEs_l = 0.05 * np.exp(0.55 * zp) * EV
lane_tau_l = np.array([1.366468e+03, 1.542332e+04, 8.137589e+03, 2.466660e+05, 2.251496e+06,
                       3.165801e+05, 8.073097e+06, 1.994613e+07, 3.417463e+08, 1.470174e+09])
verify_ensemble("LOGNORMAL", Bs_l, dEs_l, 350.0, 1e9, lane_tau_l)

# ================================================================ V5: symmetric control
say("")
say("V5. SYMMETRIC CONTROL")
tau_s = tau_pop_eig(1.0 * EV, 0.0, 350.0, 1e9)
t_half = np.exp(1.0 * EV / (KB * 350.0)) / (2e9)
say(f"    tau(dE=0) = {tau_s:.6e} s; e^(B/kT)/(2 f0) = {t_half:.6e} s (lane 1.254112e5)")
vcheck("V5: threshold = e^(B/kT)/(2f0) to 1e-12", abs(tau_s - t_half) / t_half <= 1e-12)
vcheck("V5: lane's 1.254112e5 reproduced", abs(tau_s / 1.254112e5 - 1) <= 1e-6)

# ================================================================ V6: corner limit / proxy
say("")
say("V6. CORNER LIMIT — clustered v2 proxy, own clustering code")
def v2(m):
    k = 0
    while m > 0 and m % 2 == 0: m //= 2; k += 1
    return k
def clustered_mults(levels, width):
    lv = np.sort(np.asarray(levels)); out = [[lv[0], 1]]
    for x in lv[1:]:
        if x - out[-1][0] <= width: out[-1][1] += 1; out[-1][0] = x
        else: out.append([x, 1])
    return [m for _, m in out]
def joint_levels(dEs):
    lv = np.zeros(1)
    for d in dEs:
        lv = np.concatenate([lv, lv + d])
    return lv
m_sym = clustered_mults(joint_levels(np.zeros(8)), 0.0)
vcheck("V6: symmetric corner min v2 = 8", min(v2(m) for m in m_sym) == 8, f"mults {m_sym}")
Bc, Tc, f0c = 1.2 * EV, 350.0, 1e9
dEs_c = np.array([0.000, 0.005, 0.31, 0.45, 0.60, 0.75]) * EV
taus_c = np.array([tau_pop_eig(Bc, d, Tc, f0c) for d in dEs_c])
rows = {1e1: (4, 6), 1e4: (2, 2), 1e7: (2, 2), 1e8: (0, 1)}
for t_m, (kd_lane, kv_lane) in rows.items():
    d = delta_exact(Bc, Tc, f0c, t_m)
    kd = int((taus_c >= t_m).sum())
    ml = clustered_mults(joint_levels(dEs_c), max(d, 0.0))
    kv = min(v2(m) for m in ml)
    vcheck(f"V6: common-B row t_m={t_m:.0e} -> (k={kd_lane}, v2={kv_lane})",
           kd == kd_lane and kv == kv_lane, f"got ({kd},{kv})")
t_x = 1e4
d = delta_exact(Bc, Tc, f0c, t_x)
dE_pair = np.array([0.90 * d, 1.05 * d])
taus_p = np.array([tau_pop_eig(Bc, x, Tc, f0c) for x in dE_pair])
kd = int((taus_p >= t_x).sum())
ml = clustered_mults(joint_levels(dE_pair), d)
kv = min(v2(m) for m in ml)
say(f"    chaining: taus = {taus_p[0]:.3e}, {taus_p[1]:.3e} (lane 2.678e4, 6.110e3); "
    f"k={kd}, proxy v2={kv}, mults {ml}")
vcheck("V6: chaining counterexample (proxy 2 vs k 1)", kv == 2 and kd == 1)
# is tau(1.05 delta) < t_m genuinely (the dead record)?
vcheck("V6: the 1.05-delta record is genuinely dead at t_m", taus_p[1] < t_x,
       f"tau {taus_p[1]:.3e} vs {t_x:.0e}")

# ================================================================ V7: T-31 carrier, full rebuild
say("")
say("V7. T-31 GENERIC-ASYMMETRY CARRIER — independent rebuild")
I2 = np.eye(2, dtype=complex)
def word(s):
    M = np.array([[1.0]], complex)
    for c in s: M = np.kron(M, {'I': I2, 'X': SX, 'Z': SZ}[c])
    return M
Hsym = -(word('XXXX') + word('ZZZZ'))
def generic(eps):
    return eps * sum((1.0 + 0.6180339887 * j) *
                     word(''.join('Z' if i == j else 'I' for i in range(4))) for j in range(4))
COUPLE = [word(''.join(p if i == j else 'I' for i in range(4))) for j in range(4) for p in 'XZ']
Z1, Z2 = word('ZZII'), word('ZIZI')
kT_c, f0_c = 0.20, 1.0
def davies(H):
    E, V = np.linalg.eigh(H)
    Ls = []
    for A in COUPLE:
        At = V.conj().T @ A @ V
        groups = {}
        for m in range(16):
            for n in range(16):
                if abs(At[m, n]) < 1e-12: continue
                groups.setdefault(round(float(E[m] - E[n]), 9),
                                  np.zeros((16, 16), complex))[m, n] = At[m, n]
        for wk, M in groups.items():
            Ls.append(np.sqrt(f0_c * np.exp(-max(wk, 0.0) / kT_c)) * (V @ M @ V.conj().T))
    return E, V, Ls
BASINS = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
def basin(H, Ls):
    n = 16
    E, V = np.linalg.eigh(H)
    gib = V @ np.diag(np.exp(-(E - E.min()) / kT_c)) @ V.conj().T
    Pi = {(s1, s2): 0.25 * ((np.eye(n) + s1 * Z1) @ (np.eye(n) + s2 * Z2))
          for s1 in (1, -1) for s2 in (1, -1)}
    L = gksl(H, Ls)
    w, U = np.linalg.eig(L)
    idx = np.argsort(np.abs(np.real(w)))[:4]
    lam = np.real(w[idx])
    Y = np.zeros((4, 4))
    for col, j in enumerate(idx):
        Mj = U[:, j].reshape(n, n, order='F')
        for i, b in enumerate(BASINS):
            Y[i, col] = float(np.real(np.trace(Pi[b] @ Mj)))
    Q = Y @ np.diag(lam) @ np.linalg.inv(Y)
    gen_dev = max(float(np.max(np.abs(Q.sum(axis=0)))),
                  float(-min(0.0, np.min(Q - np.diag(np.diag(Q))))))
    slow = np.sort(np.abs(lam))
    pb = np.array([float(np.real(np.trace(Pi[b] @ gib))) for b in BASINS]); pb /= pb.sum()
    taus = []
    for rec in (0, 1):
        rate = 0.0
        for s in (1, -1):
            jdx = [j for j, b in enumerate(BASINS) if b[rec] == s]
            wts = pb[jdx] / pb[jdx].sum()
            outr = [sum(Q[i, j] for i, b in enumerate(BASINS) if b[rec] != s) for j in jdx]
            rate += float(np.dot(wts, outr))
        taus.append(1.0 / rate)
    return Q, taus, gen_dev / max(slow[-1], 1e-300)
def shell(E, V):
    order = np.argsort(E); bot, sad = order[:4], E[order[4]]
    out = []
    for R in (Z1, Z2):
        Elo = {1: np.inf, -1: np.inf}
        for nn in bot:
            sval = float(np.real(V[:, nn].conj() @ (R @ V[:, nn])))
            s = int(round(sval))
            assert abs(sval - s) < 1e-6, f"basin label not sharp: {sval}"
            Elo[s] = min(Elo[s], E[nn])
        out.append((abs(Elo[1] - Elo[-1]), sad - max(Elo[1], Elo[-1])))
    return out
lane_F = {0.00: ((0.0, 0.0), (5.5066e3, 5.5066e3)),
          0.05: ((0.0698, 0.0641), (5.5008e3, 5.5147e3)),
          0.16: ((0.5684, 0.5124), (3.6284e3, 4.1648e3))}
taus16 = None
for eps, ((d1l, d2l), (t1l, t2l)) in lane_F.items():
    H = Hsym + generic(eps)
    E, V, Ls = davies(H)
    Q, taus, qdev = basin(H, Ls)
    sd = shell(E, V) if eps > 0 else [(0.0, np.sort(E)[4] - np.sort(E)[0])] * 2
    mults = clustered_mults(np.linalg.eigvalsh(H), 1e-9)
    kv = min(v2(m) for m in mults)
    say(f"    eps={eps:.2f}: dE=({sd[0][0]:.4f},{sd[1][0]:.4f}) lane ({d1l},{d2l}); "
        f"tau=({taus[0]:.4e},{taus[1]:.4e}) lane ({t1l:.4e},{t2l:.4e}); qdev={qdev:.1e}; v2={kv}")
    vcheck(f"V7: eps={eps} dE reproduced <=1e-3", abs(sd[0][0] - d1l) <= 1e-3 and abs(sd[1][0] - d2l) <= 1e-3)
    vcheck(f"V7: eps={eps} taus reproduced <=1e-3 rel",
           abs(taus[0] / t1l - 1) <= 1e-3 and abs(taus[1] / t2l - 1) <= 1e-3)
    vcheck(f"V7: eps={eps} Q classical to 1e-6 of slow scale", qdev <= 1e-6, f"{qdev:.1e}")
    vcheck(f"V7: eps={eps} exact-mult v2 = {2 if eps == 0 else 0}", kv == (2 if eps == 0 else 0))
    if eps == 0.16: taus16 = taus
tg = np.logspace(np.log10(min(taus16)) - 2, np.log10(max(taus16)) + 2, 401)
kc = np.array([sum(1 for t in taus16 if t >= x) for x in tg])
vcheck("V7: staircase at eps=0.16 takes values [2,1,0]",
       sorted(set(kc.tolist()), reverse=True) == [2, 1, 0])
# f0_eff re-derivation at eps=0: ln f0_eff = B/kT - ln tau - ln 2
H0 = Hsym; E0 = np.linalg.eigvalsh(H0)
_, taus0, _ = basin(H0, davies(H0)[2])
lnf0 = (np.sort(E0)[4] - np.sort(E0)[0]) / kT_c - np.log(taus0[0]) - np.log(2.0)
say(f"    ln f0_eff = {lnf0:.4f} (lane 0.6931 = ln 2); f0_eff = {np.exp(lnf0):.4f}")
vcheck("V7: f0_eff = 2.000 reproduced (measured at symmetric point only)",
       abs(np.exp(lnf0) - 2.0) <= 2e-3)

# ================================================================ V8: departure term
say("")
say("V8. DEPARTURE TERM — closed-form m(t) (no Liouvillian at all)")
Bg, Tg, f0g = 1.1 * EV, 350.0, 1e9
kTg = KB * Tg
xs = np.array([0.2, 0.5, 0.8, 1.2, 1.8, 2.5, 3.2, 4.0, 5.0, 6.5])
t_probe = 1e6
lane_mlo = [0.752394, 0.759500, 0.768654, 0.784738, 0.818536, 0.870783, 0.923612, 0.964039,
            0.986614, 0.996998]
lane_mup = [-0.697573, -0.603483, -0.485130, -0.285304, 0.097794, 0.574181, 0.873992, 0.963400,
            0.986614, 0.996998]
k_probe, n_pop, dep = 0, 0, 0.0
worst_m = 0.0
for i, x in enumerate(xs):
    dE = x * kTg
    gu = f0g * np.exp(-(Bg - dE) / kTg); gl = f0g * np.exp(-Bg / kTg)
    tau = 1.0 / (gu + gl)
    meq = np.tanh(dE / (2 * kTg))
    mlo = meq + (1 - meq) * np.exp(-t_probe / tau)
    mup = meq + (-1 - meq) * np.exp(-t_probe / tau)
    worst_m = max(worst_m, abs(mlo - lane_mlo[i]), abs(mup - lane_mup[i]))
    alive = tau >= t_probe
    k_probe += alive
    n_pop += (mlo > 0)
    if not alive: dep += meq
say(f"    closed-form m(t) vs lane table worst |dm| = {worst_m:.2e}")
vcheck("V8: lane's exact-propagator m(t) = rate-equation closed form <=1e-5", worst_m <= 1e-5)
vcheck("V8: k=4, Sharrock survivors 10, departure 6 grains", k_probe == 4 and n_pop == 10)
vcheck("V8: departure sum m_eq = 5.4339", abs(dep - 5.4339) <= 1e-3, f"{dep:.4f}")
for s, lane_d in ((0.5, 4.8388), (1.0, 6.6955), (2.0, 8.0865), (4.0, 9.0453)):
    dd = float(np.tanh(s * xs / 2.0).sum())
    vcheck(f"V8: departure at scale {s} = {lane_d}", abs(dd - lane_d) <= 1e-3, f"{dd:.4f}")

# ================================================================ V9: instrument-scope claims
say("")
say("V9. THE LANE'S REGISTER-WORTHY INSTRUMENT CLAIMS (importing the model's own instruments)")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), "model"))
import grounded as GR
from project_model import RecordSurface, ProjectModel
# (3a) grounded.liouvillian takes H in angular frequency; project_model passes joules.
#      Consequence claim: harmless for population records ([H,R]=0), load-bearing for coherences.
s = RecordSurface("v", "t", 0.1 * EV, 1.0 * EV, 350.0, 1e9)
H_j, Ls_j, R_j = s.open_system()
r_joule = GR.clause_ii(H_j, Ls_j, R_j, 1.0)['rate']          # H in J (the model's own path)
r_omega = GR.clause_ii(H_j / HBAR, Ls_j, R_j, 1.0)['rate']   # H in 1/s (the lane's path)
vcheck("V9a: population-mode rate identical under both H units (harmless as claimed)",
       abs(r_joule - r_omega) / r_omega <= 1e-12, f"rel {abs(r_joule-r_omega)/r_omega:.1e}")
r_coh_j = abs(mode_rate_expansion(H_j, Ls_j, SX))
r_coh_o = abs(mode_rate_expansion(H_j / HBAR, Ls_j, SX))
say(f"    coherence-mode |lambda|: H-in-joules {r_coh_j:.3e} vs H-in-omega {r_coh_o:.3e} "
    f"(ratio {r_coh_j/r_coh_o:.3e}) — load-bearing exactly as the caveat says")
vcheck("V9a: coherence mode differs wildly under the units confusion (claim confirmed)",
       r_coh_j / r_coh_o < 1e-10 or r_coh_j / r_coh_o > 1e10)
# (3b) clause_ii Rayleigh quotient on the 16-dim carrier reads ~0.497 vs true ~3.6e3
H16 = Hsym + generic(0.16)
q = GR.clause_ii(H16, davies(H16)[2], Z1, 1.0)['rate']
say(f"    grounded.clause_ii on 16-dim carrier: 1/rate = {1.0/q:.3f} (lane 0.497; true tau_1 = {taus16[0]:.3e})")
vcheck("V9b: clause_ii misread on multi-shell carrier reproduced",
       abs(1.0 / q - 0.497) <= 5e-3 and taus16[0] / (1.0 / q) > 1e3)

# ================================================================ verdict data
say("")
say("=" * 100)
say(f"VERIFY: {len(FAILS)} REFUTE lines" + ("" if not FAILS else " -> " + "; ".join(FAILS)))
say("=" * 100)
sys.exit(1 if FAILS else 0)

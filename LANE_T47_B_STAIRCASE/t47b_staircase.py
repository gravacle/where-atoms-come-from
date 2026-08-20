"""T-47 / LANE B — THE COUNT LAW k(t_m) ON MULTI-RECORD CARRIERS.

THE CANDIDATE (from the brief): a record is two-valued at retention t_m only if BOTH values are
durable, tau_record = 1/(g_u+g_l) >= t_m, giving the width delta(t_m) = B - kT ln(f0 t_m) and the
staircase k(t_m) = #{i : dE_i <= B_i - kT ln(f0 t_m)}.

THE COMPLETION (derived here, then verified against the model's own instrument): the candidate
keeps only the faster escape. The record's own Liouvillian mode decays at g_u + g_l EXACTLY
(both directions), so the survival condition  f0 t_m e^(-B/kT) (e^(dE/kT)+1) <= 1  solves in
closed form with NO free parameter:

    EXACT WIDTH      delta(t_m)  =  kT ln( e^(B/kT) / (f0 t_m)  -  1 )
    EXACT STEP       t_m*(i)     =  f0^-1 e^((B_i-dE_i)/kT) / (1 + e^(-dE_i/kT))

The candidate's width is the dE >> kT asymptote; its error is kT ln(1+e^(-dE/kT)) <= kT ln 2,
i.e. the candidate's step times are HIGH by a factor <= 2, worst at dE = 0 (where the corner
bound f0^-1 e^(B/kT) is high by exactly 2 — the population difference relaxes at g_u+g_l, the
same factor grounded.py's spectrum() docstring already warns about).

BOTH WIDTHS FROM ONE CLAUSE: clause (ii') is |lambda| <= 1/t_m on the record's own mode.
  population mode      lambda = -(g_u+g_l)                    -> the Arrhenius width above
  coherence  mode      lambda = -(g_u+g_l)/2 +- i dE/hbar     -> |lambda|<=1/t_m needs
                       dE <= hbar sqrt(1/t_m^2 - ((g_u+g_l)/2)^2) ~ hbar/t_m
The coherence width hbar/t_m REAPPEARS for off-diagonal records from the imaginary part of the
same eigenvalue; C-76's failure was using it for a population record.

Every quantity in this file is the surface's own (B_i, dE_i, T, f0 declared per carrier; t_m the
retention spec) or a DERIVED width; the only tolerances are numerical instrument floors, declared
where they appear (D-24). No files outside this lane. Returns are data."""
import sys, os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'model'))
import grounded as G
from project_model import RecordSurface

KB, HBAR = G.KB, G.HBAR
EV = 1.602176634e-19
def say(*a): print(*a); sys.stdout.flush()
CHECKS = []
def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok)))
    say(f"    [{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))

# ------------------------------------------------------------------ the derived law, log-space
def ln_tstar_exact(B, dE, T, f0):
    """ln of the exact step time from clause (ii') on the population mode."""
    kT = KB * T
    return (B - dE) / kT - np.log(f0) - np.log1p(np.exp(-dE / kT))

def ln_tstar_cand(B, dE, T, f0):
    """ln of the CANDIDATE step time (keeps only the faster escape)."""
    return (B - dE) / (KB * T) - np.log(f0)

def delta_exact(B, T, f0, t_m):
    """The derived width, closed form. Negative -> no record with this B survives t_m."""
    kT = KB * T
    arg = np.exp(B / kT) / (f0 * t_m) - 1.0
    return -np.inf if arg <= 0 else kT * np.log(arg)

def delta_cand(B, T, f0, t_m):
    return B - KB * T * np.log(f0 * t_m)

# ------------------------------------------------------------------ the instrument (no formula)
def record_mode_rate(H_omega, Ls, R):
    """The record's own Liouvillian mode, by full eigendecomposition of the adjoint generator.
       H_omega in ANGULAR-FREQUENCY units (H/hbar) so populations and coherences carry the same
       time unit. Decomposes vec(R) in the eigenbasis of L-adjoint, excludes the ONE stationary
       mode (the eigenvalue of smallest modulus — the GKSL steady state's conserved partner), and
       returns the eigenvalue of the dominant remaining component. No scale-relative floor: an
       earlier draft excluded '|lambda| < 1e-14 * spectral scale', and the coherence frequency
       dE/hbar ~ 1e15 rad/s silently swallowed the slow population mode — the same class of error
       C-75/C-76 record. Selection is by expansion weight, never by magnitude."""
    Lad = G.liouvillian(H_omega, Ls).conj().T
    w, U = np.linalg.eig(Lad)
    v = np.asarray(R, dtype=complex).reshape(-1, 1, order='F')
    v = v / np.linalg.norm(v)
    c = np.linalg.solve(U, v).ravel()
    weight = np.abs(c) * np.linalg.norm(U, axis=0)
    j0 = int(np.argmin(np.abs(w)))                 # the stationary mode, excluded by identity
    weight[j0] = -1.0
    return w[int(np.argmax(weight))]

def surface_2level(B, dE, T, f0):
    """The model's own two-state surface. Convention map: the brief's B is the barrier above the
       LOWER well; project_model's E_b is the activation from the METASTABLE (upper) well, so
       E_b = B - dE. Returns (H_omega, Ls, R)."""
    s = RecordSurface("lane", "t47b", dE, B - dE, T, f0)
    H, Ls, R = s.open_system()
    return H / HBAR, Ls, R   # instrument note: grounded.liouvillian needs H in 1/s units

def tau_direct(B, dE, T, f0):
    """Route (1): the record-mode lifetime from the eigenvalue instrument. No formula."""
    H_omega, Ls, R = surface_2level(B, dE, T, f0)
    lam = record_mode_rate(H_omega, Ls, R)
    q = G.clause_ii(H_omega, Ls, R, 1.0)['rate']            # the model's own clause instrument
    assert abs(abs(lam) - q) <= 1e-9 * max(q, 1e-300), (lam, q)
    return 1.0 / abs(lam)

say("=" * 100)
say("T-47 LANE B    THE COUNT LAW k(t_m) ON MULTI-RECORD CARRIERS")
say("=" * 100)

# ================================================================== A. THE WIDTH, DERIVED & HIT
say("")
say("A.  THE DERIVED WIDTH — clause (ii') on the record's own population mode, closed form")
say("    survive t_m  <=>  g_u+g_l <= 1/t_m  <=>  f0 t_m e^(-B/kT)(e^(dE/kT)+1) <= 1")
say("    <=>  dE <= delta(t_m) = kT ln( e^(B/kT)/(f0 t_m) - 1 )        [exact, no free parameter]")
say("    candidate delta_cand = B - kT ln(f0 t_m) is the dE>>kT asymptote; error <= kT ln 2.")
say("")
B0, T0, f0_0 = 1.2 * EV, 350.0, 1e9
say(f"    surface: B = 1.2 eV, T = 350 K, f0 = 1e9 Hz   (kT = {KB*T0/EV:.6f} eV)")
say(f"    {'t_m (s)':>10} {'delta_exact (eV)':>17} {'dE* bisected on tau(dE)=t_m (eV)':>34} "
    f"{'rel err':>10} {'delta_cand (eV)':>16} {'cand offset/kTln2':>18}")
worst_cross = 0.0
for t_m in (1e2, 1e4, 1e6, 1e7):
    dlt = delta_exact(B0, T0, f0_0, t_m)
    lo, hi = 0.0, B0                       # tau(dE) is monotone decreasing; bisect the crossing
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if tau_direct(B0, mid, T0, f0_0) >= t_m: lo = mid
        else: hi = mid
    dE_star = 0.5 * (lo + hi)
    rel = abs(dE_star - dlt) / dlt
    worst_cross = max(worst_cross, rel)
    dc = delta_cand(B0, T0, f0_0, t_m)
    off = (dc - dlt) / (KB * T0 * np.log(2.0))
    say(f"    {t_m:>10.0e} {dlt/EV:>17.12f} {dE_star/EV:>34.12f} {rel:>10.1e} {dc/EV:>16.12f} {off:>18.6f}")
check("A: numeric crossing lands on delta_exact, all t_m, rel err <= 1e-9", worst_cross <= 1e-9,
      f"worst {worst_cross:.2e}")
say("    the candidate width overshoots by kT ln(1+e^(-delta/kT)) — its 'cand offset' column tends")
say("    to 1 (=kT ln2) as delta -> 0 and to 0 as delta >> kT: the asymptote, quantified.")

# ================================================================== B. BOTH WIDTHS, ONE CLAUSE
say("")
say("B.  BOTH WIDTHS FROM THE ONE CLAUSE (D-15 control: a coherence-type record)")
B1, T1, f1 = 1.4 * EV, 350.0, 1e9
t_m1 = 3.156e7                                     # one year
H_omega, Ls, _ = surface_2level(B1, 0.0, T1, f1)   # rates at dE ~ 1e-42 J are dE-independent
g_tot = G.clause_ii(H_omega, Ls, np.diag([1.0, -1.0]).astype(complex), 1.0)['rate']
lo, hi = 0.0, 1e-38                                 # bisect |lambda_coh(dE)| = 1/t_m
SX = np.array([[0, 1], [1, 0]], dtype=complex)
for _ in range(200):
    mid = 0.5 * (lo + hi)
    Ho, Lz, _ = surface_2level(B1, mid, T1, f1)
    lam = record_mode_rate(Ho, Lz, SX)
    if abs(lam) <= 1.0 / t_m1: lo = mid
    else: hi = mid
dE_coh = 0.5 * (lo + hi)
w_coh = HBAR / t_m1
w_exact = HBAR * np.sqrt((1.0 / t_m1) ** 2 - (g_tot / 2.0) ** 2)
say(f"    surface B = 1.4 eV, T = 350 K, f0 = 1e9, t_m = 1 yr;  g_u+g_l = {g_tot:.3e} 1/s")
say(f"    coherence-mode crossing, bisected on the eigenvalue:  dE* = {dE_coh:.9e} J")
say(f"    hbar/t_m = {w_coh:.9e} J     exact hbar*sqrt(1/t_m^2-((gu+gl)/2)^2) = {w_exact:.9e} J")
check("B: coherence width = hbar/t_m reappears (<=1e-6 rel, decay corr negligible here)",
      abs(dE_coh - w_coh) / w_coh <= 1e-6, f"rel {abs(dE_coh-w_coh)/w_coh:.2e}")
d_pop = delta_exact(B1, T1, f1, t_m1)
say(f"    THE TWO WIDTHS OF THE SAME SURFACE AT THE SAME t_m, from the same |lambda| <= 1/t_m:")
say(f"      population (diagonal) record:    delta = {d_pop/EV:.4f} eV = {d_pop:.3e} J")
say(f"      coherence (off-diagonal) record: delta = {w_coh:.3e} J")
say(f"      ratio {d_pop/w_coh:.2e} — C-76 quantified: the coherence width applied to a population")
say(f"      record undercounts by ~{np.log10(d_pop/w_coh):.0f} orders of magnitude; k = 0 was that error.")

# ================================================================== C. ENSEMBLES, k TWO WAYS
say("")
say("C.  k(t_m) TWO WAYS ON DECLARED MULTI-RECORD SURFACES")
say("    route 1: eigenvalue instrument per record (tau from the Liouvillian, no formula)")
say("    route 2: the derived staircase (exact form); the candidate form tabulated beside it")

def ensemble_report(name, Bs, dEs, T, f0, tgrid):
    N = len(Bs)
    say("")
    say(f"    -- {name}:  N = {N}, T = {T} K, f0 = {f0:.0e} Hz")
    say(f"    {'i':>3} {'B_i (eV)':>10} {'dE_i (eV)':>10} {'tau_i eig (s)':>14} "
        f"{'t*_exact (s)':>14} {'|d ln|':>9} {'t*_cand/t*_exact':>17}")
    taus, lt_ex, lt_cd = [], [], []
    worst = 0.0
    for i in range(N):
        tau = tau_direct(Bs[i], dEs[i], T, f0)
        le, lc = ln_tstar_exact(Bs[i], dEs[i], T, f0), ln_tstar_cand(Bs[i], dEs[i], T, f0)
        taus.append(tau); lt_ex.append(le); lt_cd.append(lc)
        dln = abs(np.log(tau) - le); worst = max(worst, dln)
        say(f"    {i:>3} {Bs[i]/EV:>10.4f} {dEs[i]/EV:>10.4f} {tau:>14.6e} "
            f"{np.exp(le):>14.6e} {dln:>9.1e} {np.exp(lc-le):>17.6f}")
    taus, lt_ex, lt_cd = np.array(taus), np.array(lt_ex), np.array(lt_cd)
    lg = np.log(tgrid)
    k1 = (taus[None, :] >= tgrid[:, None]).sum(axis=1)
    k2 = (lt_ex[None, :] >= lg[:, None]).sum(axis=1)
    kc = (lt_cd[None, :] >= lg[:, None]).sum(axis=1)
    d12, dc = np.abs(k1 - k2), np.abs(k1 - kc)
    say(f"    k_direct vs k_staircase(exact):    max discrepancy over {len(tgrid)} t_m = {d12.max()}"
        f"   (grid points differing: {(d12>0).sum()})")
    say(f"    k_direct vs k_staircase(candidate): max discrepancy = {dc.max()}, at "
        f"{(dc>0).sum()} grid points — the candidate holds each record a factor (1+e^(-dE/kT)) too long")
    say(f"    step-location law: max |ln tau_eig - ln t*_exact| over records = {worst:.2e}")
    check(f"C[{name}]: two routes agree at every t_m", d12.max() == 0)
    check(f"C[{name}]: steps land on t* = f0^-1 e^((B-dE)/kT)/(1+e^(-dE/kT)), <=1e-10",
          worst <= 1e-10, f"worst {worst:.2e}")
    say(f"    staircase is decreasing in ln t_m: {bool(np.all(np.diff(k1) <= 0))}; "
        f"k({tgrid[0]:.0e} s) = {k1[0]}, k({tgrid[-1]:.0e} s) = {k1[-1]}")
    return taus

tgrid = np.logspace(0, 12, 1201)
# UNIFORM: declared grids, no sampling
Bs_u  = np.array([0.95 + 0.40 * i / 11.0 for i in range(12)]) * EV
dEs_u = np.array([0.18 * ((7 * i + 3) % 12) / 11.0 for i in range(12)]) * EV
tau_u = ensemble_report("UNIFORM (declared grid)", Bs_u, dEs_u, 350.0, 1e9, tgrid)
# LOGNORMAL-LIKE DISCRETE: declared z-grid, B = 1.1 e^(0.12 z), dE = 0.05 e^(0.55 z') eV
zs  = np.array([-1.5, -1.1, -0.7, -0.35, 0.0, 0.0, 0.35, 0.7, 1.1, 1.5])
zp  = np.array([0.7, -0.35, 1.5, 0.0, -1.5, 1.1, -0.7, 0.35, -1.1, 0.0])
Bs_l  = 1.10 * np.exp(0.12 * zs) * EV
dEs_l = 0.05 * np.exp(0.55 * zp) * EV
tau_l = ensemble_report("LOGNORMAL-LIKE DISCRETE (declared z-grid)", Bs_l, dEs_l, 350.0, 1e9, tgrid)

# ================================================================== D. D-15 SYMMETRIC CONTROL
say("")
say("D.  D-15 CONTROL — SYMMETRIC CARRIER: k(t_m) must be FLAT to a single threshold")
Bs_s, dEs_s = np.full(8, 1.0 * EV), np.zeros(8)
taus_s = np.array([tau_direct(Bs_s[i], 0.0, 350.0, 1e9) for i in range(8)])
t_half = np.exp(1.0 * EV / (KB * 350.0)) / (2 * 1e9)
t_cand = np.exp(1.0 * EV / (KB * 350.0)) / 1e9
tg = np.logspace(0, 8, 801)
ks = (taus_s[None, :] >= tg[:, None]).sum(axis=1)
say(f"    8 records, B = 1.0 eV, dE = 0: tau_i (eig) all equal: "
    f"{taus_s.min():.6e} .. {taus_s.max():.6e} s")
say(f"    k = 8 on every grid t_m < tau, k = 0 above: values taken = {sorted(set(ks.tolist()))}")
flat = np.all((ks == 8) | (ks == 0)) and ks[0] == 8 and ks[-1] == 0
check("D: symmetric carrier flat at 8, single step to 0", flat)
say(f"    threshold measured {taus_s[0]:.6e} s = e^(B/kT)/(2 f0) = {t_half:.6e} s;")
say(f"    the candidate corner 'any t_m below f0^-1 e^(B/kT)' = {t_cand:.6e} s is HIGH by exactly 2:")
say(f"    at dE = 0 both wells escape, tau = 1/(2 f0 e^(-B/kT)). Corner limit holds, corrected by 2.")
check("D: dE=0 threshold = e^(B/kT)/(2f0) to 1e-12", abs(taus_s[0] - t_half) / t_half <= 1e-12)

# ================================================================== E. C-14 / C-76 RECONCILIATION
say("")
say("E.  THE CORNER LIMIT — mapping the staircase back to C-14's v_2 arithmetic (closes C-76)")
def v2(m):
    k = 0
    while m % 2 == 0 and m > 0: m //= 2; k += 1
    return k
def clustered_mults(levels, width):
    lv = np.sort(np.asarray(levels)); out = [[lv[0], 1]]
    for x in lv[1:]:
        if x - out[-1][0] <= width: out[-1][1] += 1; out[-1][0] = x   # chain rule, as in T-31
        else: out.append([x, 1])
    return [m for _, m in out]
def joint_levels(dEs):
    lv = [0.0]
    for d in dEs: lv = [x for x in lv] + [x + d for x in lv]
    return np.array(lv)

# symmetric corner: the arithmetic reconnects exactly
m_sym = clustered_mults(joint_levels(np.zeros(8)), 0.0)
say("    SYMMETRIC N=8 corner: joint spectrum = 2^8 coincident levels -> "
    + (f"multiplicities {m_sym}, min v_2 = {min(v2(m) for m in m_sym)} = k_direct = 8"
       if min(v2(m) for m in m_sym) == 8 else f"MISMATCH: {m_sym}"))
check("E: symmetric corner v_2 arithmetic = staircase count", min(v2(m) for m in m_sym) == 8)

# common-B carrier: clustered-v2 with the DERIVED width vs the record-mode count, across t_m
Bc, Tc, f0c = 1.2 * EV, 350.0, 1e9
dEs_c = np.array([0.000, 0.005, 0.31, 0.45, 0.60, 0.75]) * EV
taus_c = np.array([tau_direct(Bc, d, Tc, f0c) for d in dEs_c])
say(f"    COMMON-B CARRIER (B = 1.2 eV): dE_i/eV = {np.round(dEs_c/EV,3).tolist()}")
say(f"    {'t_m (s)':>10} {'delta(t_m) eV':>14} {'k_direct':>9} {'min v_2 clustered':>18}  multiplicities")
agree = mism = 0
for t_m in np.logspace(1, 8, 29):
    dlt = delta_exact(Bc, Tc, f0c, t_m)
    kd = int((taus_c >= t_m).sum())
    ml = clustered_mults(joint_levels(dEs_c), max(dlt, 0.0))
    kv = min(v2(m) for m in ml)
    if kv == kd: agree += 1
    else: mism += 1
    if t_m in (1e1, 1e4, 1e7) or kv != kd:
        say(f"    {t_m:>10.1e} {dlt/EV:>14.4f} {kd:>9} {kv:>18}  {ml}")
say(f"    agreement over 29 t_m: {agree} agree, {mism} disagree. WHERE it agrees: live splittings")
say(f"    well inside delta, dead ones well outside (the [4,56,4]-type rows: min v_2 = 2 = k).")
say(f"    WHERE it disagrees: (a) large delta CHAIN-MERGES the whole spectrum (the [64] rows read")
say(f"    6 while k = 4 or 3); (b) at delta <= 0 an EXACTLY degenerate pair still doubles every")
say(f"    level (the [2,2,...] row reads 1 while k = 0: degeneracy without durability — C-76's own")
say(f"    lesson). The arithmetic reconnects EXACTLY at the symmetric corner and in the separated")
say(f"    regime; clustered-v_2 is a spectral PROXY, and the law must be registered on record")
say(f"    modes, where both routes agree at every t_m (section C).")
# the honest boundary of the proxy: chaining pulls a DEAD record's level into a live cluster
t_x = 1e4
dlt = delta_exact(Bc, Tc, f0c, t_x)
dE_pair = np.array([0.90 * dlt, 1.05 * dlt])
taus_p = np.array([tau_direct(Bc, d, Tc, f0c) for d in dE_pair])
kd = int((taus_p >= t_x).sum())
ml = clustered_mults(joint_levels(dE_pair), dlt)
kv = min(v2(m) for m in ml)
say(f"    PROXY BOUNDARY (chaining): B = 1.2 eV, t_m = 1e4 s, delta = {dlt/EV:.4f} eV,")
say(f"      dE = (0.90, 1.05) delta -> tau = ({taus_p[0]:.3e}, {taus_p[1]:.3e}) s;")
say(f"      record-mode count k = {kd}; chain-clustered multiplicities {ml} -> min v_2 = {kv}.")
say(f"      The spectral proxy overcounts ({kv} vs {kd}): consecutive gaps 0.90d, 0.15d, 0.90d all")
say(f"      <= delta chain a dead record's level into the cluster. THE LAW LIVES ON RECORD MODES;")
say(f"      clustered-v_2 is a corner proxy, exact at dE = 0, not the law.")
check("E: chaining counterexample exhibits proxy overcount (v2=2 vs k=1)", kv == 2 and kd == 1)

# ================================================================== F. THE T-31 CARRIER
say("")
say("F.  THE T-31 GENERIC-ASYMMETRY CARRIER — the kill becomes a staircase")
say("    [[4,2,2]] Hsym = -(XXXX+ZZZZ) + generic Z-fields eps*(1+0.618 j). Carrier units DECLARED:")
say("    hbar = 1, f0 = 1 (Metropolis attempt rate), kT = 0.20 x the carrier's own unit; bath =")
say("    Davies-Metropolis single-qubit X_j, Z_j couplings, detailed balance; the numerical Bohr-")
say("    frequency grouping floor is 1e-9 (instrument floor, declared).")
I2 = np.eye(2, dtype=complex)
PX = np.array([[0, 1], [1, 0]], dtype=complex)
PZ = np.array([[1, 0], [0, -1]], dtype=complex)
def word(s):
    M = np.array([[1.0]], dtype=complex)
    for c in s: M = np.kron(M, {'I': I2, 'X': PX, 'Z': PZ}[c])
    return M
Hsym = -(word('XXXX') + word('ZZZZ'))
def generic(eps):
    return eps * sum((1.0 + 0.6180339887 * j) *
                     word(''.join('Z' if i == j else 'I' for i in range(4))) for j in range(4))
COUPLE = [word(''.join(p if i == j else 'I' for i in range(4))) for j in range(4) for p in 'XZ']
ZBAR = [word('ZZII'), word('ZIZI')]
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
                key = round(float(E[m] - E[n]), 9)
                groups.setdefault(key, np.zeros((16, 16), complex))[m, n] = At[m, n]
        for wk, Mat in groups.items():
            gam = f0_c * np.exp(-max(wk, 0.0) / kT_c)
            Ls.append(np.sqrt(gam) * (V @ Mat @ V.conj().T))
    return E, V, Ls
def shell_data(E, V):
    """(B_i, dE_i) read off the carrier's own spectrum: bottom shell labeled by the exact
       conserved Zbar_i; saddle = the lowest excited-shell level."""
    order = np.argsort(E)
    bot, sad = order[:4], E[order[4]]
    out = []
    for R in ZBAR:
        Elo = {+1: np.inf, -1: np.inf}
        for n in bot:
            s = int(round(float(np.real(V[:, n].conj() @ (R @ V[:, n])))))
            Elo[s] = min(Elo[s], E[n])
        dE = abs(Elo[+1] - Elo[-1]); Eup = max(Elo[+1], Elo[-1])
        out.append(dict(dE=dE, B_minus_dE=sad - Eup))
    return out

BASINS = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
def basin_instrument(H, Ls):
    """The carrier's own slow sector, no formula and no rate theory: basins = the joint
       (Zbar_1, Zbar_2) sectors; the 4-state generator Q is read off the FULL 256-dim
       Liouvillian's own 4 slowest eigenmodes through the basin readout Y[i,j] =
       tr(Pi_i unvec(u_j)):  Q = Y diag(lambda_slow) Y^-1 — exact at slow times by construction.
       (An instantaneous-outflow tr(Pi_b' L[T_b]) draft overcounted every rate by ~2 — the
       transition-state recrossing factor: a walker excited to the mid shell falls back to either
       side. Reading Q off the generator's own slow modes removes it; reported, not patched.)
       VALIDATED below: Q must come out a CLASSICAL MARKOV GENERATOR (columns sum to 0,
       off-diagonal rates >= 0) — that the slow sector is a classical 4-well chain is then
       measured, not assumed. Per-record lifetime = the record's own mode of its lumped
       two-state chain, lumping weights the carrier's own conditional equilibrium."""
    n = H.shape[0]
    Ev, Vv = np.linalg.eigh(H)
    gib = Vv @ np.diag(np.exp(-(Ev - Ev.min()) / kT_c)) @ Vv.conj().T
    Pi = {(s1, s2): 0.25 * ((np.eye(n) + s1 * ZBAR[0]) @ (np.eye(n) + s2 * ZBAR[1]))
          for s1 in (1, -1) for s2 in (1, -1)}
    L = G.liouvillian(H, Ls)
    w, U = np.linalg.eig(L)
    idx = np.argsort(np.abs(np.real(w)))[:4]
    lam = np.real(w[idx])
    assert np.max(np.abs(np.imag(w[idx]))) <= 1e-12 * max(1.0, np.max(np.abs(lam)))
    Y = np.zeros((4, 4))
    for col, j in enumerate(idx):
        Mj = U[:, j].reshape(n, n, order='F')
        for i, b in enumerate(BASINS):
            Y[i, col] = float(np.real(np.trace(Pi[b] @ Mj)))
    Q = Y @ np.diag(lam) @ np.linalg.inv(Y)
    gen_dev = max(float(np.max(np.abs(Q.sum(axis=0)))),
                  float(-min(0.0, np.min(Q - np.diag(np.diag(Q))))))
    pb = np.array([float(np.real(np.trace(Pi[b] @ gib))) for b in BASINS])
    pb = pb / pb.sum()
    taus = []
    for rec in (0, 1):
        rate = 0.0
        for s in (1, -1):                       # escape from value s, spectator at cond. equil.
            jdx = [j for j, b in enumerate(BASINS) if b[rec] == s]
            wts = pb[jdx] / pb[jdx].sum()
            out_rate = [sum(Q[i, j] for i, b in enumerate(BASINS) if b[rec] != s) for j in jdx]
            rate += float(np.dot(wts, out_rate))
        taus.append(1.0 / rate)
    return Q, taus, gen_dev, np.sort(np.abs(lam))

say("")
say("    per-record mode: the basin-lumped two-state chain (weights = the carrier's own")
say("    conditional Gibbs; hop rates read off the full 256-dim generator, no formula).")
say(f"    {'eps':>6} {'dE_1':>8} {'dE_2':>8} {'tau_1 basin':>12} {'tau_2 basin':>12} "
    f"{'tau_1 pred':>12} {'tau_2 pred':>12} {'|d ln| 1':>9} {'|d ln| 2':>9} {'v2 exact k':>11}")
lnf0_eff = None
rows_f = []
worst_qval = 0.0
for eps in (0.0, 0.05, 0.10, 0.16):
    H = Hsym + generic(eps)
    E, V, Ls = davies(H)
    Q, taus, gen_dev, slow_rates = basin_instrument(H, Ls)
    qdev = gen_dev / max(slow_rates[-1], 1e-300)   # generator violation, relative to slow scale
    worst_qval = max(worst_qval, qdev)
    mults = clustered_mults(np.linalg.eigvalsh(H), 1e-9)
    kv2 = min(v2(m) for m in mults)
    if eps == 0.0:
        sd = [dict(dE=0.0, B_minus_dE=np.sort(E)[4] - np.sort(E)[0])] * 2
        # f0_eff: the carrier's own attempt rate, MEASURED once at the symmetric point (a carrier
        # constant like f0 itself; not fitted to any eps > 0 prediction)
        lnf0_eff = sd[0]['B_minus_dE'] / kT_c - np.log(taus[0]) - np.log(2.0)
    else:
        sd = shell_data(E, V)
    preds = [np.exp(s['B_minus_dE'] / kT_c - lnf0_eff - np.log1p(np.exp(-s['dE'] / kT_c)))
             for s in sd]
    dl = [abs(np.log(taus[i] / preds[i])) for i in range(2)]
    rows_f.append(dict(eps=eps, taus=taus, preds=preds, dl=dl, kv2=kv2, sd=sd, qdev=qdev))
    say(f"    {eps:>6.2f} {sd[0]['dE']:>8.4f} {sd[1]['dE']:>8.4f} {taus[0]:>12.4e} {taus[1]:>12.4e} "
        f"{preds[0]:>12.4e} {preds[1]:>12.4e} {dl[0]:>9.3f} {dl[1]:>9.3f} {kv2:>11}")
say(f"    (f0_eff measured once at eps = 0: ln f0_eff = {lnf0_eff:.4f}; B_i - dE_i and dE_i read off")
say(f"    the carrier's own spectrum — saddle = lowest excited level, wells = labeled bottom shell)")
say(f"    basin-generator validation: Q carries the full generator's own 4 slowest eigenvalues by")
say(f"    construction, and comes out a CLASSICAL MARKOV GENERATOR — column sums / negative")
say(f"    off-diagonals deviate by at most {worst_qval:.2e} of the slow scale, over all eps. The")
say(f"    slow sector IS a classical 4-well chain: measured, not assumed.")
check("F: slow sector is a classical 4-state generator (<= 1e-6 of slow scale)",
      worst_qval <= 1e-6, f"worst {worst_qval:.2e}")
r = rows_f[-1]
t_lo, t_hi = min(r['taus']), max(r['taus'])
tgc = np.logspace(np.log10(t_lo) - 2, np.log10(t_hi) + 2, 401)
kc = np.array([int(sum(1 for t in r['taus'] if t >= x)) for x in tgc])
say(f"    STAIRCASE at eps = 0.16: k = 2 for t_m < {t_lo:.4e}, k = 1 to {t_hi:.4e}, k = 0 above")
say(f"      (values taken on the grid: {sorted(set(kc.tolist()), reverse=True)}); exact-multiplicity")
say(f"      v_2 count at the same eps = {r['kv2']} at EVERY t_m — T-31's kill row — while the")
say(f"      record-mode count steps 2 -> 1 -> 0. The binary 0 was the corner arithmetic, not the count.")
check("F: generic-asymmetry carrier gives k=2 below computed threshold, steps down",
      list(sorted(set(kc.tolist()), reverse=True)) == [2, 1, 0] and r['kv2'] == 0)
check("F: step order tracks the derived exponent (smaller B-dE dies first)",
      (r['sd'][0]['B_minus_dE'] < r['sd'][1]['B_minus_dE']) == (r['taus'][0] < r['taus'][1]))
worst_f = max(max(x['dl']) for x in rows_f if x['eps'] > 0)
say(f"    formula-vs-instrument residual on the abstract carrier: max |d ln tau| = {worst_f:.3f}")
say(f"    over an exponent (B-dE)/kT of ~{min(s['B_minus_dE'] for s in rows_f[-1]['sd'])/kT_c:.1f}-10.0;")
say(f"    the residual is the multi-pathway prefactor the two-well formula compresses into one f0")
say(f"    (reported, not tuned away). The LAW's content here is the staircase and its exponent.")
check("F: exponent carries the step (residual << step separation ln(t2*/t1*)+2 ln-units)",
      worst_f < abs(np.log(r['taus'][1] / r['taus'][0])) + 2.0)
say("    INSTRUMENT-SCOPE FINDING (register-worthy): grounded.clause_ii's Rayleigh quotient on this")
q_flat = G.clause_ii(Hsym + generic(0.16), davies(Hsym + generic(0.16))[2], ZBAR[0], 1.0)['rate']
say(f"    16-dim carrier reads 1/rate = {1.0/q_flat:.3f} (vs the true {r['taus'][0]:.3e}): the trace")
say(f"    inner product weights all shells uniformly and the excited-shell flips swamp the record's")
say(f"    slow mode. On multi-shell carriers the record's own mode must be taken in the basin-lumped")
say(f"    slow sector — same failure class as the C-69/spectrum() 2x-mode lesson, one level up.")

# ================================================================== G. THE DEPARTURE TERM
say("")
say("G.  THE DEPARTURE TERM — Sharrock-style population counting vs k(t_m), same carrier")
say("    Sharrock / Street-Woolley / Neel decay counting integrates surviving MAGNETISATION;")
say("    k(t_m) counts FUNCTIONING TWO-VALUED RECORDS. On an asymmetric carrier they differ.")
Bg, Tg, f0g = 1.1 * EV, 350.0, 1e9
kTg = KB * Tg
xs = np.array([0.2, 0.5, 0.8, 1.2, 1.8, 2.5, 3.2, 4.0, 5.0, 6.5])   # dE_i/kT, declared
dEs_g = xs * kTg
t_probe = 1e6
say(f"    carrier: 10 grains, B = 1.1 eV, T = 350 K, f0 = 1e9; dE_i/kT declared = {xs.tolist()}")
say(f"    probe t = {t_probe:.0e} s.  m_i(t) from the Liouvillian propagator e^(Lt), exact.")
say(f"    {'i':>3} {'dE/kT':>6} {'tau_i (s)':>12} {'record@t':>9} {'m(t) lower-written':>19} "
    f"{'m(t) upper-written':>19} {'m_eq':>7}")
def m_of_t(B, dE, T, f0, t, written):
    H_omega, Ls, R = surface_2level(B, dE, T, f0)
    L = G.liouvillian(H_omega, Ls)
    w, U = np.linalg.eig(L)
    rho0 = np.diag([1.0, 0.0]) if written == 'lower' else np.diag([0.0, 1.0])
    c = np.linalg.solve(U, rho0.reshape(-1, 1, order='F'))
    rt = (U @ (np.exp(w[:, None] * t) * c)).reshape(2, 2, order='F')
    return float(np.real(np.trace(rt @ R)))
dep_grains = []
for i, dE in enumerate(dEs_g):
    tau = tau_direct(Bg, dE, Tg, f0g)
    alive = tau >= t_probe
    mlo = m_of_t(Bg, dE, Tg, f0g, t_probe, 'lower')
    mup = m_of_t(Bg, dE, Tg, f0g, t_probe, 'upper')
    meq = np.tanh(dE / (2 * kTg))
    if not alive: dep_grains.append((i, mlo, meq))
    say(f"    {i:>3} {xs[i]:>6.1f} {tau:>12.4e} {'ALIVE' if alive else 'dead':>9} "
        f"{mlo:>19.6f} {mup:>19.6f} {meq:>7.4f}")
k_probe = sum(1 for dE in dEs_g if tau_direct(Bg, dE, Tg, f0g) >= t_probe)
n_pop = sum(1 for dE in dEs_g if m_of_t(Bg, dE, Tg, f0g, t_probe, 'lower') > 0)
say(f"    at t = {t_probe:.0e} s: functioning records k = {k_probe};  Sharrock-counted surviving")
say(f"    moments (lower-written, sign criterion) = {n_pop} of 10.  THE EXACT CONFIGURATIONS THAT")
say(f"    DIFFER: grains {[g[0] for g in dep_grains]} — upper value decayed (tau_i < t), lower value")
say(f"    persists at m ~ m_eq; magnetisation counted, record not. Departure = {n_pop - k_probe} grains;")
say(f"    magnetisation-weighted departure = sum over dead of m_eq = "
    f"{sum(g[2] for g in dep_grains):.4f} (of N = 10).")
check("G: departure configurations exist and are exactly the dead-but-magnetised grains",
      n_pop - k_probe == len(dep_grains) and len(dep_grains) > 0)
say("")
say("    departure vs the asymmetry distribution (same grains, dE scaled by s; t at which all dead):")
say(f"    {'scale s':>8} {'<dE>/kT':>8} {'k(t->inf)':>10} {'Sharrock survivors':>19} "
    f"{'departure sum m_eq':>19}")
for s in (0.5, 1.0, 2.0, 4.0):
    meqs = np.tanh(s * xs / 2.0)
    say(f"    {s:>8.1f} {s*xs.mean():>8.2f} {0:>10} {10:>19} {meqs.sum():>19.4f}")
say("    the departure GROWS with asymmetry (-> N as dE/kT -> inf) and vanishes on the symmetric")
say("    carrier (m_eq = 0): Sharrock counting and k(t_m) coincide only at dE = 0 — the corner.")
say("    Sharrock's own object M(t) = sum m_i(t) NEVER counts two-valuedness: an exchange-biased")
say("    medium DC-written along bias holds M ~ sum tanh(dE_i/2kT) FOREVER while k(t_m) -> 0.")

# ================================================================== VERDICT DATA
say("")
say("=" * 100)
npass = sum(1 for _, ok in CHECKS if ok); nfail = len(CHECKS) - npass
say(f"  {npass} PASS, {nfail} FAIL over {len(CHECKS)} checks")
say("=" * 100)
sys.exit(1 if nfail else 0)

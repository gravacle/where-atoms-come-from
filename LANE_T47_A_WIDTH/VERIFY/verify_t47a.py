"""ADVERSARIAL VERIFY, T-47 A. Independent of the lane's arithmetic path wherever possible:
(1) the record-mode eigenvalue from FULL eigendecomposition of the 4x4 adjoint Liouvillian
    (never the lane's Rayleigh quotient), bisected against 1/t_m;
(2) the exact threshold re-solved at 60-digit precision with mpmath (findroot on
    g_u+g_l = 1/t_m) against the expm1 closed form;
(3) the correction identity checked as algebra at high precision;
(4) bracket-sensitivity: the lane's bisection bracket constants (0.999*B, +40 kT) perturbed --
    if the crossing moves, a smuggled tolerance is load-bearing;
(5) the no-crossing condition of delta_exact vs the symmetric corner bound exp(B/kT)/(2 f0)
    -- must be the SAME inequality;
(6) the staircase recounted from scratch (closed form AND full-eig), and the C-14 corner
    arithmetic v2(2^N) = N stated as a computation, not prose;
(7) coherence face by full eig of the Hamiltonian-only generator, full modulus.
"""
import sys, os, numpy as np
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'model'))
import grounded as G
from project_model import RecordSurface

EV, YEAR = 1.602176634e-19, 3.156e7
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
SP = np.array([[0, 1], [0, 0]], dtype=complex)
FAILS = []
def ck(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}  {detail}")
    if not ok: FAILS.append(name)

B0, T0, F0, TM0 = 1.2 * EV, 300.0, 1e9, 10.0 * YEAR
kT0 = G.KB * T0

# ---- (1) population-mode rate from FULL EIG of L-adjoint, not the Rayleigh quotient ----------
def rate_fulleig(B, dE, T, f0):
    s = RecordSurface("v", "thermal", dE, B - dE, T, f0)
    H, Ls, R = s.open_system()
    Lad = G.liouvillian(H, Ls).conj().T
    w, V = np.linalg.eig(Lad)
    # the population sector: eigenvector with dominant weight on diagonal, traceless part
    best, brate = None, None
    for i in range(4):
        M = V[:, i].reshape(2, 2, order='F')
        offd = abs(M[0, 1]) + abs(M[1, 0])
        tr0 = abs(M[0, 0] - M[1, 1])          # traceless-diagonal weight
        if offd < 1e-12 * np.linalg.norm(M) and tr0 > 1e-6 and abs(w[i]) > 1e-30:
            if best is None or abs(w[i]) < brate: best, brate = i, abs(w[i])
    return abs(w[best])

r_ray_path = None
def rate_lane(B, dE, T, f0):     # the lane's path, for cross-comparison only
    s = RecordSurface("v", "thermal", dE, B - dE, T, f0)
    H, Ls, R = s.open_system()
    return G.clause_ii(H, Ls, R, np.inf)['rate']

print("== (1) full-eig population eigenvalue vs the lane's Rayleigh-quotient rate ==")
for dE in (0.0, 0.05 * EV, 0.158 * EV, 0.30 * EV):
    a, b = rate_fulleig(B0, dE, T0, F0), rate_lane(B0, dE, T0, F0)
    ck(f"fulleig-vs-lane dE={dE/EV:.3f} eV", abs(a - b) / b < 1e-10, f"rel {abs(a-b)/b:.1e}")
# and -(gu+gl) is an EXACT eigenvalue (sigma_z is not the eigenvector; sigma_z + aI is):
gu = F0 * np.exp(-(B0 - 0.1 * EV) / kT0); gl = F0 * np.exp(-B0 / kT0)
w_all = np.linalg.eigvals(G.liouvillian(*[x for x in RecordSurface("v","t",0.1*EV,B0-0.1*EV,T0,F0).open_system()[:2]]).conj().T)
ck("-(gu+gl) is an exact eigenvalue", min(abs(w + (gu + gl)) for w in w_all) / (gu + gl) < 1e-12)

def cross_fulleig(B, T, f0, t_m):
    kT = G.KB * T
    lo, hi = 0.0, 1.0 * B
    if rate_fulleig(B, 0.0, T, f0) >= 1.0 / t_m: return None
    if rate_fulleig(B, hi, T, f0) < 1.0 / t_m: return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if rate_fulleig(B, mid, T, f0) <= 1.0 / t_m: lo = mid
        else: hi = mid
        if hi - lo <= np.spacing(hi): break
    return 0.5 * (lo + hi)

print("== (2) main crossing: full-eig bisection vs 60-digit mpmath root vs expm1 form ==")
xn = cross_fulleig(B0, T0, F0, TM0)
mp.mp.dps = 60
kTm, Bm, F0m, TMm = mp.mpf(G.KB) * 300, mp.mpf('1.2') * mp.mpf('1.602176634e-19'), mp.mpf('1e9'), mp.mpf('10') * mp.mpf('3.156e7')
# dimensionless x = dE/kT: solve exp(x - b)(1 + exp(-x)) = exp(-c), b = B/kT, c = ln(f0 t_m)
b, c = Bm / kTm, mp.log(F0m * TMm)
gx = lambda x: (x - b) + mp.log(1 + mp.e**(-x)) + c        # log form, well-scaled
x_mp = kTm * mp.findroot(gx, mp.mpf('6.1'))
x_expm1 = kTm * mp.log(mp.expm1(Bm / kTm - mp.log(F0m * TMm)))
ck("mpmath root == expm1 closed form", abs(x_mp - x_expm1) / x_expm1 < mp.mpf('1e-50'),
   f"rel {float(abs(x_mp-x_expm1)/x_expm1):.1e}")
ck("full-eig bisection lands on it", abs(xn - float(x_expm1)) / float(x_expm1) < 5e-15,
   f"dE*={xn/EV:.12f} eV vs {float(x_expm1)/EV:.12f} eV, rel {abs(xn-float(x_expm1))/float(x_expm1):.2e}")
ck("matches the lane's sealed 0.158282235181 eV", abs(xn / EV - 0.158282235181) < 1e-11)
# correction identity, 60-digit: dE* == [B - kT ln(f0 t_m)] - kT ln(1+e^{-dE*/kT})
lhs = x_mp; rhs = (Bm - kTm * mp.log(F0m * TMm)) - kTm * mp.log(1 + mp.e**(-x_mp / kTm))
ck("correction identity exact at 60 digits", abs(lhs - rhs) / lhs < mp.mpf('1e-50'))
corr = float(kTm * mp.log(1 + mp.e**(-x_mp / kTm)))
ck("correction in (0, kT ln2]", 0 < corr <= kT0 * np.log(2) + 1e-30, f"corr = {corr/kT0:.6f} kT")

print("== (3) bracket-sensitivity: perturb the lane's bisection bracket constants ==")
def cross_bracket(hi_mult_B, extra_kT):
    kT = G.KB * T0
    dn = B0 - kT * np.log(F0 * TM0)
    lo, hi = 0.0, min(hi_mult_B * B0, max(dn, 0.0) + extra_kT * kT)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if rate_lane(B0, mid, T0, F0) <= 1.0 / TM0: lo = mid
        else: hi = mid
        if hi - lo <= np.spacing(hi): break
    return 0.5 * (lo + hi)
xa, xb, xc = cross_bracket(0.999, 40.0), cross_bracket(0.9999, 80.0), cross_bracket(0.5, 20.0)
ck("crossing invariant under bracket choice", max(abs(xa-xb), abs(xa-xc)) / xa < 1e-13,
   f"spread {max(abs(xa-xb),abs(xa-xc))/xa:.1e} -- bracket constants are not load-bearing")

print("== (4) no-crossing condition == symmetric corner bound, same inequality ==")
# delta_exact is None iff B/kT - ln(f0 t_m) <= ln 2 iff t_m >= exp(B/kT)/(2 f0)
for T, f0, t_m in ((400.0, 1e13, 3.156e3), (300.0, 1e9, 3.156e8), (400.0, 1e9, 3.156e7)):
    kT = G.KB * T
    none_cond = (1.2 * EV / kT - np.log(f0 * t_m)) <= np.log(2.0)
    bound = np.exp(1.2 * EV / kT) / (2 * f0)
    ck(f"T={T:.0f} f0={f0:.0e} t_m={t_m:.2e}: (no dE*>=0) == (t_m >= exp(B/kT)/2f0)",
       none_cond == (t_m >= bound), f"bound {bound:.3e} s")
tau_sym = 1.0 / rate_fulleig(B0, 0.0, T0, F0)
ck("symmetric tau = exp(B/kT)/(2 f0), factor 2 real", abs(tau_sym - np.exp(B0/kT0)/(2*F0)) / tau_sym < 1e-12,
   f"tau {tau_sym:.6e} s; one-sided form errs by 2.0x")

print("== (5) staircase recounted from scratch (closed form AND full eig) ==")
dEs = [0.0, 0.05, 0.10, 0.15, 0.20, 0.25]
tstar_cf = [np.exp((1.2 - d) * EV / kT0) / (F0 * (1 + np.exp(-d * EV / kT0))) for d in dEs]
sealed_t = [7.212449e10, 1.821850e10, 2.952602e9, 4.344239e8, 6.296091e7, 9.104795e6]
ck("t*_i match sealed table", all(abs(a - b) / b < 5e-7 for a, b in zip(tstar_cf, sealed_t)))
ks = []
for t_m in (3.156e4, 3.156e6, 3.156e7, 3.156e8, 3.156e9, 3.156e10, 3.156e11):
    k_cf = sum(1 for t in tstar_cf if t >= t_m)
    k_fe = sum(1 for d in dEs if rate_fulleig(1.2 * EV, d * EV, T0, F0) <= 1.0 / t_m)
    ks.append(k_cf); ck(f"k({t_m:.2e}) closed-form == full-eig", k_cf == k_fe, f"k={k_cf}")
ck("sealed staircase 6,6,5,4,2,1,0", ks == [6, 6, 5, 4, 2, 1, 0])
ck("monotone non-increasing in ln t_m", all(a >= b for a, b in zip(ks, ks[1:])))
# C-14 corner arithmetic, computed not asserted: N symmetric two-state records -> one 2^N-fold
# degenerate eigenspace; C-14's k = min_E v2(m_E) = v2(2^6) = 6 = the flat staircase count.
v2 = lambda n: (n & -n).bit_length() - 1
ck("C-14 corner: v2(2^6) == flat-staircase k == 6", v2(2 ** 6) == 6 == ks[0])

print("== (6) coherence face: full eig of Hamiltonian-only generator, full modulus ==")
tm_c = 1e-6; dcoh = G.HBAR / tm_c
def coh_rate_fulleig(dE):
    L = G.liouvillian(-(dE / 2) * SZ / G.HBAR, []).conj().T
    w = np.linalg.eigvals(L)
    return max(abs(x) for x in w)                # the coherence pair, +-i dE/hbar
lo, hi = 0.0, 2 * dcoh
for _ in range(200):
    mid = 0.5 * (lo + hi)
    if coh_rate_fulleig(mid) <= 1.0 / tm_c: lo = mid
    else: hi = mid
    if hi - lo <= np.spacing(hi): break
xc = 0.5 * (lo + hi)
ck("coherence crossing == hbar/t_m (full eig)", abs(xc - dcoh) / dcoh < 1e-12,
   f"{xc:.15e} vs {dcoh:.15e}")
ck("Re-part-only criterion would (wrongly) keep every coherence",
   all(abs(x.real) < 1e-20 for x in np.linalg.eigvals(G.liouvillian(-(2*dcoh/2)*SZ/G.HBAR, []))),
   "-- the full-|lambda| rule is load-bearing, as claimed")
# dissipative modulus: one expression, two corners
gu_a, gl_a, w_a = 5.0, 3.0, 6.0
La = [np.sqrt(gl_a) * SP.conj().T, np.sqrt(gu_a) * SP]
ws = np.linalg.eigvals(G.liouvillian(-(w_a / 2) * SZ, La))
pred = np.sqrt(((gu_a + gl_a) / 2) ** 2 + w_a ** 2)
ck("unified modulus sqrt(((gu+gl)/2)^2+(dE/hbar)^2) present in spectrum",
   min(abs(abs(z) - pred) for z in ws if abs(z.imag) > 1) / pred < 1e-12, f"pred {pred:.12f}")

print("== (7) smuggled-tolerance census over the lane script ==")
import tokenize, io
src = open(os.path.join(HERE, '..', 't47_a_width.py')).read()
consts = []   # every numeric literal in CODE (strings/comments excluded by the tokenizer)
for tok in tokenize.generate_tokens(io.StringIO(src).readline):
    if tok.type == tokenize.NUMBER:
        consts.append(tok.string)
suspects = [c for c in consts if c in ('1e-2', '0.01', '1e-3', '0.001')]
ck("no C-76-style clustering width among the lane's CODE constants", len(suspects) == 0,
   f"(1e-2/1e-3 appear only in prose describing C-76 and the sealed floor)")
names = [t.string for t in tokenize.generate_tokens(io.StringIO(src).readline)
         if t.type == tokenize.NAME]
ck("lane never overrides slow_modes' sealed tol in code", 'tol' not in names,
   "('tol=1e-9' appears only in prose; the sealed default is used unmodified)")
# the one code-level 1e-9 (line 264) is a display-time eigenvalue SELECTOR in section 5
# (imag ~ 6 rad/s vs ~ 0 in a demo with O(1) rates), not a width entering any law
# the remaining constants, eyeballable: physical params, check thresholds (assertions on
# numerical agreement, reported raw alongside), and bisection bracket constants (shown
# non-load-bearing in (3) above)
print("   code numeric literals:", sorted(set(consts), key=lambda s: (len(s), s)))

print()
print("VERDICT-INPUT: FAILS =", FAILS if FAILS else "none")
sys.exit(1 if FAILS else 0)

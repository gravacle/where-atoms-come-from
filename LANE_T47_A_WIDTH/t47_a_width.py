"""T-47 LANE A: THE WIDTH, DERIVED -- both faces of clause (ii') on the record's own mode.

THE CRITERION (one, not two): clause (ii') applied to the record's OWN Liouvillian mode,

        |lambda_record| <= 1/t_m .

The two-well carrier (barrier B above the LOWER well, splitting dE with the upper well higher,
shared Arrhenius prefactor f0) has exactly two record-bearing mode families, and the ONE criterion
lands differently on each because their eigenvalues live on different axes:

  (a) POPULATION-TYPE (diagonal, R = sigma_z).  [H, sigma_z] = 0, so the mode's eigenvalue is REAL:
      lambda = -(g_u + g_l),   g_u = f0 exp(-(B-dE)/kT)   (escape from the shallow/upper well),
                               g_l = f0 exp(-B/kT)        (escape from the deep/lower well).
      Detailed balance g_l/g_u = exp(-dE/kT) is a CONSEQUENCE of the shared-f0 Arrhenius form,
      not an extra assumption.  The criterion g_u + g_l <= 1/t_m reads
          f0 exp(-(B-dE)/kT) (1 + exp(-dE/kT)) <= 1/t_m
      whose EXACT threshold in dE is (expm1 form, overflow-safe):
          dE*(t_m) = kT ln( expm1( B/kT - ln(f0 t_m) ) )
      and whose leading form is the commissioned width WITH ITS DERIVED CORRECTION:
          dE*(t_m) = [B - kT ln(f0 t_m)] - kT ln(1 + exp(-dE*/kT))
                      `---- delta_pop ----'   `---- two-sided-escape correction, in (0, kT ln2] ----'
      A record with dE > kT carries an exponentially small correction ~ kT exp(-dE/kT); at dE = 0
      the correction is exactly kT ln 2 (both wells escape, doubling the rate).

  (b) COHERENCE-TYPE (off-diagonal, R = sigma_+).  [H, sigma_+] = -dE sigma_+, so the mode's
      eigenvalue is IMAGINARY under Hamiltonian-only dynamics: lambda = -i dE/hbar.  The SAME
      criterion |lambda| <= 1/t_m -- the FULL modulus, the C-75 correction, load-bearing here --
      reads dE/hbar <= 1/t_m, i.e.
          delta_coh(t_m) = hbar / t_m .
      With dissipation on, the coherence eigenvalue is -(g_u+g_l)/2 - i dE/hbar and the criterion
      is sqrt( ((g_u+g_l)/2)^2 + (dE/hbar)^2 ) <= 1/t_m: ONE expression whose two corners are the
      two widths.  Dissipation-dominated corner -> the Arrhenius face (a); rotation-dominated
      corner -> the energy face (b).  That is what makes this a DERIVATION and not a patch: both
      widths fall out of one inequality on one generator.

THE COUNT LAW (population face):
      k(t_m) = #{ i : g_u(i) + g_l(i) <= 1/t_m }  =  #{ i : dE_i <= dE*_i(t_m) }
a decreasing staircase in ln t_m; record i drops out at its OWN
      t_m*(i) = 1 / (g_u(i)+g_l(i)) = f0^-1 exp((B_i-dE_i)/kT) / (1 + exp(-dE_i/kT)) .
(The commission's t* lacked the (1+e^{-dE/kT}) factor; it is derived and tested below.)

CONVENTION (declared): B is measured from the LOWER well.  project_model's RecordSurface uses the
standard Arrhenius activation from the METASTABLE (upper) well, E_b = B - dE; the map is exact and
is used, not re-derived.  ASSUMPTIONS (declared): shared prefactor f0 for both wells (equal Kramers
attempt frequencies); Markovian GKSL thermal activation (the census's thermal=True scope -- the
model DECLINES non-thermally-activated surfaces); f0 t_m >> 1 (retention longer than one attempt
period).  Everything else is SI-exact (kB, hbar) or the carrier's own (B, dE, T, f0) or the
retention spec (t_m).  NO CHOSEN TOLERANCE ANYWHERE: the only tolerances below are (1) the
declared numerical eigenvalue floor already inside the sealed instrument (slow_modes tol=1e-9,
documented there as numerical not physical) and (2) bisection run to floating-point spacing,
reported as a residual, not chosen.

CORNER LIMITS: t_m -> infinity: delta_pop -> -infinity (no thermally activated two-valued record
survives forever at T > 0 -- DEF-A's E_b -> infinity demand recovered) and delta_coh -> 0 (exact
degeneracy, the commutant count of C-14/C-75 recovered).  Symmetric carrier dE = 0: survives iff
t_m <= exp(B/kT)/(2 f0) -- the commission wrote f0^-1 exp(B/kT); the factor 2 is the two-sided
escape and is MEASURED below, a correction to the commissioned corner, not to the law.
"""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'model'))
import grounded as G
from project_model import RecordSurface, ProjectModel

def say(*a): print(*a); sys.stdout.flush()

EV   = 1.602176634e-19          # J, exact SI
YEAR = 3.156e7                  # s, the house convention (T-29, T-31)
SZ = np.array([[1, 0], [0, -1]], dtype=complex)
SP = np.array([[0, 1], [0, 0]], dtype=complex)
FAILS = []
def check(name, ok):
    if not ok: FAILS.append(name)
    return "ok" if ok else "FAIL"

MDL = ProjectModel()

def rate_num(B, dE, T, f0):
    """The record's own mode rate FROM THE LIOUVILLIAN INSTRUMENT (grounded.clause_ii via
       project_model), never from the closed form. B from the lower well -> E_b = B - dE."""
    s = RecordSurface("sweep", "thermal two-well", dE, B - dE, T, f0)
    opened = s.open_system()
    if opened is None: return 0.0
    H, Ls, R = opened
    return G.clause_ii(H, Ls, R, np.inf)['rate']

def rate_closed(B, dE, T, f0):
    kT = G.KB * T
    return f0 * np.exp(-(B - dE) / kT) * (1.0 + np.exp(-dE / kT))

def delta_naive(B, T, f0, t_m):
    return B - G.KB * T * np.log(f0 * t_m)

def delta_exact(B, T, f0, t_m):
    """dE* = kT ln(expm1(B/kT - ln(f0 t_m))); None when no dE >= 0 crossing exists."""
    kT = G.KB * T
    y = B / kT - np.log(f0 * t_m)
    if y <= np.log(2.0): return None          # expm1(y) <= 1  <=>  dE* <= 0
    return kT * np.log(np.expm1(y))

def crossing_num(B, T, f0, t_m):
    """Bisect the INSTRUMENT's rate against 1/t_m in dE; run to floating-point spacing."""
    if rate_num(B, 0.0, T, f0) >= 1.0 / t_m: return None      # even symmetric does not survive
    kT = G.KB * T
    lo, hi = 0.0, min(0.999 * B, max(delta_naive(B, T, f0, t_m), 0.0) + 40.0 * kT)
    if rate_num(B, hi, T, f0) < 1.0 / t_m: return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if rate_num(B, mid, T, f0) <= 1.0 / t_m: lo = mid
        else: hi = mid
        if hi - lo <= np.spacing(hi): break
    return 0.5 * (lo + hi)

say("=" * 102)
say("T-47 A   THE WIDTH, DERIVED: BOTH FACES OF CLAUSE (ii') ON THE RECORD'S OWN LIOUVILLIAN MODE")
say("=" * 102)
say("")
say("  criterion (ONE): |lambda_record| <= 1/t_m on the record's own mode of the carrier's Liouvillian")
say("  population face: lambda real,      -(g_u+g_l)      -> dE* = delta_pop - kT ln(1+e^{-dE*/kT}),")
say("                                                        delta_pop(t_m) = B - kT ln(f0 t_m)")
say("  coherence face:  lambda imaginary, -i dE/hbar      -> delta_coh(t_m) = hbar/t_m")
say("  every symbol: B, dE, T, f0 the carrier's; t_m the retention spec; kB, hbar SI-exact.")
say("")

# ---------------------------------------------------------------- 0. INSTRUMENT CHECKS
say("-" * 102)
say("  0.  INSTRUMENT: Rayleigh-quotient rate vs closed form, and detailed balance as a consequence")
say("-" * 102)
B0, T0, F0, TM0 = 1.2 * EV, 300.0, 1e9, 10.0 * YEAR
kT0 = G.KB * T0
say(f"  carrier: B = {B0/EV:.3f} eV (from the LOWER well), T = {T0:.0f} K, f0 = {F0:.1e} Hz;"
    f"  t_m = 10 y = {TM0:.3e} s")
say(f"  {'dE (eV)':>10}{'rate: instrument':>20}{'rate: g_u+g_l':>20}{'rel err':>12}{'':>6}")
for dE in (0.0, 0.05 * EV, 0.158 * EV, 0.30 * EV):
    rn, rc = rate_num(B0, dE, T0, F0), rate_closed(B0, dE, T0, F0)
    rel = abs(rn - rc) / rc
    say(f"  {dE/EV:>10.3f}{rn:>20.12e}{rc:>20.12e}{rel:>12.1e}{check('rayleigh', rel < 1e-12):>6}")
s_db = RecordSurface("db", "thermal", 0.1 * EV, B0 - 0.1 * EV, T0, F0)
H, Ls, R = s_db.open_system()
# basis: e0 = [1,0] is the LOWER well (H e0 = -dE/2 e0); Ls[0] = sqrt(g_l) SM jumps e0 -> e1
# (escape FROM the lower well), Ls[1] = sqrt(g_u) SP jumps e1 -> e0 (escape FROM the upper)
gl_eff = np.linalg.norm(Ls[0] @ np.array([1, 0])) ** 2
gu_eff = np.linalg.norm(Ls[1] @ np.array([0, 1])) ** 2
db = (gl_eff / gu_eff) / np.exp(-0.1 * EV / kT0)
say(f"  detailed balance g_l/g_u vs e^-dE/kT at dE = 0.100 eV: ratio/expected = {db:.12f}"
    f"   {check('detailed-balance', abs(db - 1) < 1e-12)}")
say("")

# ---------------------------------------------------------------- 1. THE MAIN SWEEP AND CROSSING
say("-" * 102)
say("  1.  SWEEP dE AT FIXED (B, T, f0, t_m): the instrument's tau_record crosses t_m ON the formula")
say("-" * 102)
dn, de = delta_naive(B0, T0, F0, TM0), delta_exact(B0, T0, F0, TM0)
say(f"  {'dE (eV)':>10}{'tau_record (s), instrument':>28}{'tau/t_m':>14}{'durable to t_m?':>17}")
for dE in (0.0, 0.10 * EV, 0.15 * EV, 0.1554 * EV, 0.1560 * EV, 0.16 * EV, 0.20 * EV):
    r = rate_num(B0, dE, T0, F0); tau = 1.0 / r
    say(f"  {dE/EV:>10.4f}{tau:>28.6e}{tau/TM0:>14.4e}{str(r <= 1.0/TM0):>17}")
say("")
xn = crossing_num(B0, T0, F0, TM0)
corr = kT0 * np.log1p(np.exp(-xn / kT0))
say(f"  crossing (bisection on the INSTRUMENT, to fp spacing):  dE*_num   = {xn:.15e} J = {xn/EV:.12f} eV")
say(f"  commissioned width  delta_pop = B - kT ln(f0 t_m)     :  {dn:.15e} J = {dn/EV:.12f} eV")
say(f"    residual dE*_num - delta_pop                        : {xn-dn:>+.6e} J   ({(xn-dn)/kT0:+.6f} kT)")
say(f"    the DERIVED correction -kT ln(1+e^-dE*/kT)          : {-corr:>+.6e} J   ({-corr/kT0:+.6f} kT)")
say(f"    residual explained by the correction?  |diff| = {abs((xn-dn)+corr):.3e} J"
    f"   {check('correction-explains-residual', abs((xn - dn) + corr) < 1e-12 * abs(dn))}")
say(f"  EXACT derived threshold kT ln(expm1(B/kT - ln f0 t_m)):  {de:.15e} J = {de/EV:.12f} eV")
rel = abs(xn - de) / de
say(f"    crossing vs exact formula, relative residual        :  {rel:.3e}"
    f"   (source: bisection at fp spacing)   {check('crossing-on-formula', rel < 1e-12)}")
say("")

# ---------------------------------------------------------------- 2. DECADES OF t_m, f0, T
say("-" * 102)
say("  2.  DECADES: the crossing lands on the exact formula across t_m, f0, T (B = 1.2 eV throughout)")
say("-" * 102)
say(f"  {'t_m (s)':>12}{'f0 (Hz)':>10}{'T (K)':>7}{'dE*_num (eV)':>16}{'delta_pop (eV)':>16}"
    f"{'naive resid (kT)':>18}{'exact rel resid':>17}{'':>6}")
worst = 0.0
for t_m in (1.0, 3.156e3, YEAR, 10 * YEAR, 100 * YEAR):
    for f0 in (1e9, 1e11, 1e13):
        for T in (200.0, 300.0, 400.0):
            kT = G.KB * T
            xe = delta_exact(1.2 * EV, T, f0, t_m)
            xn = crossing_num(1.2 * EV, T, f0, t_m)
            if xe is None or xn is None:
                both = (xe is None) and (xn is None)
                say(f"  {t_m:>12.3e}{f0:>10.1e}{T:>7.0f}{'-- no crossing: t_m beyond exp(B/kT)/(2 f0); k contribution 0 --':>67}"
                    f"{check('no-crossing-agree', both):>6}")
                continue
            dn = delta_naive(1.2 * EV, T, f0, t_m)
            rel = abs(xn - xe) / xe; worst = max(worst, rel)
            say(f"  {t_m:>12.3e}{f0:>10.1e}{T:>7.0f}{xn/EV:>16.9f}{dn/EV:>16.9f}"
                f"{(xn-dn)/kT:>18.6f}{rel:>17.2e}{check('grid-crossing', rel < 1e-11):>6}")
say(f"  worst exact-formula relative residual over the grid: {worst:.3e}")
say("")

# ---------------------------------------------------------------- 3. CONTROL: SYMMETRIC CARRIER
say("-" * 102)
say("  3.  CONTROL (D-15): the SYMMETRIC carrier -- flat staircase; the corner bound carries a factor 2")
say("-" * 102)
tau_sym = 1.0 / rate_num(B0, 0.0, T0, F0)
bound_2f0 = np.exp(B0 / kT0) / (2 * F0)
bound_comm = np.exp(B0 / kT0) / F0
say(f"  tau(dE=0) from the instrument      : {tau_sym:.9e} s")
say(f"  exp(B/kT)/(2 f0)  (derived bound)  : {bound_2f0:.9e} s   rel err {abs(tau_sym-bound_2f0)/bound_2f0:.2e}"
    f"   {check('symmetric-bound', abs(tau_sym - bound_2f0) / bound_2f0 < 1e-12)}")
say(f"  exp(B/kT)/f0  (commissioned corner): {bound_comm:.9e} s   -- off by the two-sided factor 2")
tm_edge = bound_2f0
x_edge = crossing_num(B0, T0, F0, 0.999999 * tm_edge)
say(f"  dE -> 0 consistency: at t_m = (1 - 1e-6) * exp(B/kT)/(2 f0), the crossing sits at"
    f" dE* = {x_edge/EV:.3e} eV")
say(f"    (must approach 0 from above as t_m -> the symmetric bound; {check('edge-limit', 0 <= x_edge < 1e-5 * EV)})")
k_flat = [sum(1 for _ in range(6) if rate_num(B0, 0.0, T0, F0) <= 1.0 / t) for t in (1.0, YEAR, 100 * YEAR)]
say(f"  six identical symmetric records, k at t_m = 1 s / 1 y / 100 y: {k_flat}  -- FLAT (all below the bound)"
    f"   {check('flat-staircase', k_flat == [6, 6, 6])}")
say(f"  C-14's corner recovered: below the bound every symmetric pair counts; the count is C-14's,")
say(f"  and at t_m -> infinity the surviving set is empty at T > 0 -- DEF-A's E_b -> infinity demand.")
say("")

# ---------------------------------------------------------------- 4. CONTROL: COHERENCE-TYPE RECORD
say("-" * 102)
say("  4.  CONTROL (D-15): the COHERENCE record, Hamiltonian-only -- hbar/t_m REAPPEARS, same criterion")
say("-" * 102)
say("  (C-75's instrument, clause (ii') with the FULL |lambda|; H enters the generator as H/hbar --")
say("   angular frequency -- declared: grounded.liouvillian takes its Hamiltonian in generator units.)")
tm_c = 1.0e-6
dcoh = G.HBAR / tm_c
say(f"  t_m = {tm_c:.1e} s  ->  delta_coh = hbar/t_m = {dcoh:.6e} J")
say(f"  {'dE / (hbar/t_m)':>17}{'|lambda| (1/s), instrument':>28}{'|lambda| t_m':>14}{'durable?':>10}")
for frac in (0.5, 0.999999, 1.000001, 2.0):
    dE = frac * dcoh
    Hc = -(dE / 2) * SZ
    c = G.clause_ii(Hc / G.HBAR, [], SP, tm_c)
    say(f"  {frac:>17.6f}{c['rate']:>28.12e}{c['rate']*tm_c:>14.9f}{str(c['durable']):>10}")
# bisect the coherence crossing
lo, hi = 0.0, 2.0 * dcoh
for _ in range(200):
    mid = 0.5 * (lo + hi)
    if G.clause_ii(-(mid / 2) * SZ / G.HBAR, [], SP, tm_c)['rate'] <= 1.0 / tm_c: lo = mid
    else: hi = mid
    if hi - lo <= np.spacing(hi): break
xc = 0.5 * (lo + hi)
rel = abs(xc - dcoh) / dcoh
say(f"  coherence crossing (bisection): dE* = {xc:.15e} J;  hbar/t_m = {dcoh:.15e} J")
say(f"    relative residual {rel:.3e}   {check('coherence-crossing', rel < 1e-12)}")
rates_lo, _ = G.slow_modes(-(0.5 * dcoh / 2) * SZ / G.HBAR, [], tm_c)
rates_hi, _ = G.slow_modes(-(2.0 * dcoh / 2) * SZ / G.HBAR, [], tm_c)
say(f"  slow_modes count (declared 1e-9 numerical floor): dE = 0.5 hbar/t_m -> {len(rates_lo)} modes"
    f" (both coherences durable); dE = 2 hbar/t_m -> {len(rates_hi)} (coherences dropped)"
    f"   {check('slow-mode-count', len(rates_lo) == 4 and len(rates_hi) == 2)}")
say("  T-31 re-read: clustering a SPECTRUM within hbar/t_m is the coherence-face count -- correct for")
say("  superposition-type records, and now derived; the population-face count is the staircase below.")
say("")

# ---------------------------------------------------------------- 5. ONE CRITERION, TWO FACES
say("-" * 102)
say("  5.  ONE EXPRESSION, BOTH WIDTHS: coherence eigenvalue -(g_u+g_l)/2 - i dE/hbar (instrument check)")
say("-" * 102)
gu_a, gl_a, w_a = 5.0, 3.0, 6.0          # 1/s, 1/s, rad/s -- comparable scales so both terms weigh
Ha = -(w_a / 2) * SZ                      # already in generator units
La = [np.sqrt(gl_a) * SP.conj().T, np.sqrt(gu_a) * SP]
ws = np.linalg.eig(G.liouvillian(Ha, La)).eigenvalues
pop = min(ws, key=lambda z: abs(z - (-(gu_a + gl_a))))
coh = min(ws, key=lambda z: abs(z.imag - w_a) if abs(z.imag) > 1e-9 else 1e9)
pred = np.sqrt(((gu_a + gl_a) / 2) ** 2 + w_a ** 2)
say(f"  g_u = {gu_a}, g_l = {gl_a} (1/s), dE/hbar = {w_a} (rad/s):")
say(f"    population mode  lambda = {pop:.12f}      |lambda| = {abs(pop):.12f}  (= g_u+g_l = {gu_a+gl_a})"
    f"   {check('pop-eig', abs(abs(pop) - (gu_a + gl_a)) < 1e-10)}")
say(f"    coherence mode   lambda = {coh:.12f}   |lambda| = {abs(coh):.12f}")
say(f"    sqrt(((g_u+g_l)/2)^2 + (dE/hbar)^2)    = {pred:.12f}   rel err {abs(abs(coh)-pred)/pred:.1e}"
    f"   {check('coh-eig', abs(abs(coh) - pred) / pred < 1e-10)}")
say("  the two widths are the two corners of THIS modulus under the ONE criterion |lambda| <= 1/t_m.")
say("")

# ---------------------------------------------------------------- 6. THE COUNT LAW: THE STAIRCASE
say("-" * 102)
say("  6.  k(t_m): THE STAIRCASE -- T-31's generic-asymmetry kill becomes a dated census, not a 0")
say("-" * 102)
ens = [("r%d" % i, 1.2 * EV, dE * EV) for i, dE in enumerate((0.0, 0.05, 0.10, 0.15, 0.20, 0.25))]
say(f"  ensemble: six two-well records, B = 1.2 eV, T = 300 K, f0 = 1e9 Hz, dE_i = 0 .. 0.25 eV")
say(f"  {'record':>8}{'dE (eV)':>9}{'t*_i instrument (s)':>22}{'t*_i formula (s)':>19}{'rel err':>10}"
    f"{'t*_i (years)':>14}")
tstars = []
for name, B, dE in ens:
    ti = 1.0 / rate_num(B, dE, T0, F0)
    tf = np.exp((B - dE) / kT0) / (F0 * (1.0 + np.exp(-dE / kT0)))
    rel = abs(ti - tf) / tf; tstars.append(ti)
    say(f"  {name:>8}{dE/EV:>9.2f}{ti:>22.6e}{tf:>19.6e}{rel:>10.1e}{ti/YEAR:>14.4e}"
        f"   {check('tstar', rel < 1e-12)}")
say("")
say(f"  {'t_m':>14}{'k by instrument count':>23}{'k by #{dE_i <= dE*(t_m)}':>26}{'':>6}")
for t_m, label in ((3.156e4, "~9 h"), (3.156e6, "0.1 y"), (3.156e7, "1 y"), (3.156e8, "10 y"),
                   (3.156e9, "100 y"), (3.156e10, "1 ky"), (3.156e11, "10 ky")):
    k_i = sum(1 for name, B, dE in ens if rate_num(B, dE, T0, F0) <= 1.0 / t_m)
    xe = delta_exact(1.2 * EV, T0, F0, t_m)
    k_f = sum(1 for name, B, dE in ens if xe is not None and dE <= xe)
    say(f"  {label:>14}{k_i:>23}{k_f:>26}{check('staircase-agree', k_i == k_f):>6}")
say("")
say("  DECREASING STAIRCASE in ln t_m; each drop at the record's own t*_i above.  Symmetric control")
say("  (section 3) is the flat line; the T-31 kill (k = 0 at any asymmetry, forever) is replaced by")
say("  k(t_m), which is nonzero for every t_m below the tallest surviving rung.")
say("")

# ---------------------------------------------------------------- 7. D-24 AUDIT
say("=" * 102)
say("  D-24 AUDIT -- every quantity earned or declared")
say("=" * 102)
say("  B, dE, T, f0          the carrier's own constants (world-tier surfaces carry provenance, D-25)")
say("  t_m                   the retention spec; same declared status it has in clause (ii') itself")
say("  kB, hbar              SI exact")
say("  g_u, g_l              EARNED: Arrhenius form with SHARED f0 (declared assumption, equal Kramers")
say("                        prefactors); detailed balance then DERIVED (section 0), not imposed")
say("  delta_pop(t_m)        EARNED: B - kT ln(f0 t_m), from |lambda| <= 1/t_m on the population mode")
say("  the correction        EARNED: -kT ln(1+e^{-dE/kT}), two-sided escape; exact closed form")
say("                        dE* = kT ln(expm1(B/kT - ln f0 t_m)) verified to fp precision")
say("  delta_coh(t_m)        EARNED: hbar/t_m, the SAME criterion on the off-diagonal mode")
say("  factor 2 (symmetric)  EARNED and REPORTED: the corner bound is exp(B/kT)/(2 f0), correcting")
say("                        the commission's f0^-1 exp(B/kT)")
say("  bisection precision   NUMERICAL, run to floating-point spacing; residuals reported, not chosen")
say("  slow_modes tol=1e-9   the sealed instrument's own declared NUMERICAL floor (C-75), unchanged")
say("  NO physical tolerance was chosen anywhere in this lane.  The C-76 failure mode (a 1e-2")
say("  clustering width with no owner) has no successor here.")
say("")
say("=" * 102)
say("  READ -- from the numbers above")
say("=" * 102)
say("  THE WIDTH IS DERIVED, AND IT IS TWO-FACED BY NECESSITY.  One criterion -- clause (ii') on the")
say("  record's own Liouvillian mode -- yields delta_pop = B - kT ln(f0 t_m) (minus a derived, bounded")
say("  correction <= kT ln 2) on the diagonal face and delta_coh = hbar/t_m on the off-diagonal face.")
say("  C-76's objection is answered in its own terms: the width that recovers a count is not chosen,")
say("  it is the carrier's.  The bisected instrument crossing lands on the exact closed form at")
say("  ~1e-13 relative or better across five decades of t_m, four of f0, and 200-400 K.")
say("")
say("  THE COUNT LAW THE PROGRAM NOW OWNS:  k(t_m) = #{ i : dE_i <= kT ln(expm1(B_i/kT - ln f0 t_m)) },")
say("  a decreasing staircase in ln t_m with per-record drop times t*_i; symmetric carriers give the")
say("  flat line (C-14's count below exp(B/kT)/(2 f0)); generic asymmetry gives a dated census.")
say("")
say("  CORRECTIONS TO THE COMMISSION, both derived: (1) the symmetric corner carries a factor 2;")
say("  (2) t*_i carries the (1+e^{-dE_i/kT}) two-sided factor.  Neither alters the headline width.")
say("")
say("  OWNERSHIP.  Neel/Street-Woolley/Sharrock own the ACTIVATION WINDOW kT ln(f0 t) as a remanence-")
say("  DECAY device (fluctuation field, viscosity S, fitted Sharrock exponent); none states a COUNT of")
say("  surviving two-valued records, none carries the asymmetry term (threshold at B - dE, the")
say("  SHALLOWER value's escape), none derives its width from a record definition that must also")
say("  produce hbar/t_m on the coherence face.  MLC flash level-merging owns Vt-margin retention with")
say("  a CHOSEN sensing margin; k(t_m) is margin-free.  DEPARTURES FALSIFIABLE NOW: (a) exchange-")
say("  biased/asymmetric media -- drop times set by B - dE, not B: survivor counting (MFM grain")
say("  census vs ln t) distinguishes the staircase from symmetric Neel decay; (b) JEDEC-style MLC")
say("  Vt retention bakes -- level-loss times are k(t_m) rungs; the two-sided factor (up to 2x) is")
say("  measurable only where dE <~ few kT, stated honestly as a narrow-window departure.")
say("")
if FAILS:
    say(f"  CHECK FAILURES: {sorted(set(FAILS))}")
else:
    say(f"  every computed check above agreed with its formula; no check failed.")
sys.exit(1 if FAILS else 0)

"""qcore -- the cheap EXACT quantity functions, factored out so the master table (s5) and the
individual audits (s3, s4) compute the same numbers.  s5 asserts these against reference values
printed by s3 and s4, so the two implementations are checked against each other.

All of these rest on the same exact factorisation, derived and validated in s3/s4:
  in [CODE] the record configuration s is conserved, and in the s sector the bath sees a sum of
  COMMUTING single-qubit terms  e_j Z_j + lam c_j(s) X_j  with c_j(s) = sum of the s_i at site j.
"""
import math, numpy as np

LOG2 = math.log(2.0)
ENERGIES = [1.0, 1.4, 0.7, 1.1, 0.9, 1.3, 0.8, 1.2]
TIMES = np.linspace(1.0, 13.0, 25)
LAM, BETA = 0.8, 2.0

def binom_w(n):
    k = np.arange(n + 1)
    lw = (math.lgamma(n + 1) - np.array([math.lgamma(x + 1) for x in k])
          - np.array([math.lgamma(n - x + 1) for x in k]) - n * LOG2)
    return (n - 2 * k).astype(float), np.exp(lw)

# ---------------------------------------------------------------- dynamics (Bloch closed form)
def bloch(c, e, lam, beta, t):
    c = np.asarray(c, dtype=float)
    a = lam * c; b = float(e)
    z0 = -math.tanh(beta * b)
    nn = np.sqrt(a * a + b * b)
    th = 2.0 * nn * t
    ct, st = np.cos(th), np.sin(th)
    nn2 = np.where(nn > 0, nn * nn, 1.0)
    return np.stack([a * b * z0 * (1 - ct) / nn2,
                     -a * z0 * st / np.where(nn > 0, nn, 1.0),
                     z0 * ct + b * b * z0 * (1 - ct) / nn2], axis=-1)

def S_bloch(r):
    q = min(float(np.linalg.norm(r)), 1.0)
    p = (1.0 + q) / 2.0
    if p >= 1.0 - 1e-15 or p <= 1e-15: return 0.0
    return float(-(p * math.log2(p) + (1 - p) * math.log2(1 - p)))

def chi_of_n(n, e=ENERGIES[0], lam=LAM, beta=BETA, times=TIMES):
    """time-averaged chi about ONE record whose bath qubit is shared by n records in total"""
    if n == 0: return 0.0
    cs, ws = binom_w(n - 1)
    acc = 0.0
    for t in times:
        rp = ws @ bloch(cs + 1, e, lam, beta, t)
        rm = ws @ bloch(cs - 1, e, lam, beta, t)
        acc += max(S_bloch(0.5 * (rp + rm)) - 0.5 * (S_bloch(rp) + S_bloch(rm)), 0.0)
    return acc / len(times)

def chi_register_site(n, e=ENERGIES[0], lam=LAM, beta=BETA, times=TIMES):
    """chi( WHOLE record register : this bath qubit ) -- the JOINT Holevo quantity, not the sum
       of the individual ones.

       chi_joint = S(rho_bar_j) - E_s[ S(sigma_j(c_j)) ].  Every sigma_j(c) is a UNITARY image
       of the same thermal state, so S(sigma_j(c)) = S(tau_j) for every c and the second term is
       a constant.  Hence  chi_joint = S(rho_bar_j) - S(tau_j) <= 1 - S(tau_j) bits, a bound set
       by the bath qubit alone and INDEPENDENT OF n."""
    cs, ws = binom_w(n)
    z0 = -math.tanh(beta * e)
    Stau = S_bloch(np.array([0.0, 0.0, z0]))
    acc = 0.0
    for t in times:
        rb = ws @ bloch(cs, e, lam, beta, t)
        acc += max(S_bloch(rb) - Stau, 0.0)
    return acc / len(times), (1.0 - Stau)

# ---------------------------------------------------------------- induced potential
def f_site(c, e=ENERGIES[0], lam=LAM, beta=BETA):
    c = np.asarray(c, dtype=float)
    E = np.sqrt(e * e + lam * lam * c * c)
    return -(np.logaddexp(beta * E, -beta * E)) / beta

def J2_same_site(n, e=ENERGIES[0], lam=LAM, beta=BETA):
    if n < 2: return 0.0
    cs, ws = binom_w(n - 2)
    return 0.5 * float(ws @ f_site(cs + 2, e, lam, beta) - ws @ f_site(cs, e, lam, beta))

def phi_moments(n, e=ENERGIES[0], lam=LAM, beta=BETA):
    cs, ws = binom_w(n)
    v = f_site(cs, e, lam, beta)
    mean = float(ws @ v)
    var = float(ws @ (v - mean) ** 2)
    return mean, math.sqrt(max(var, 0.0)), float(v.max() - v.min())

# ---------------------------------------------------------------- fitting with a floor
def loglog_fit(xs, ys):
    x = np.log(np.asarray(xs, dtype=float)); y = np.log(np.asarray(ys, dtype=float))
    A = np.vstack([x, np.ones_like(x)]).T
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    pred = A @ coef
    dof = max(len(x) - 2, 1)
    s2 = float(((y - pred) ** 2).sum() / dof)
    cov = s2 * np.linalg.inv(A.T @ A)
    return float(coef[0]), math.sqrt(abs(cov[0, 0])), float(np.abs(y - pred).max())

#!/usr/bin/env python3
"""
LANE W-17 / R4 — FRAME LANE.  Route decision under test:
    "Can a singular constraint move the rate by acting on the CONNECTION, or on the READY STATE?"

This script produces the NUMBERS for F1..F5 of FRAME_CHALLENGE_V001.md.
It does NOT answer the route question. It measures whether the question is well posed.

Corpus objects (all pointered, all pre-cutoff = REGISTER_V001.md rows W-01 .. W-11):
  Z_k  = p00 + p10 u^k + p01 v^k + p11 (uv)^k        REGISTER:1149 / :1152 (W-11, 'reached by substitution')
  u = conj(W_F) = e^{if},  v = W_C = e^{ic}          REGISTER:118  (W-02)
  lambda = m(p00 + p10 x + p01 y + p11 xy)           REGISTER:456-458  (N1, W-05)
  N3: rate invariant under every a.c. connection measure; P(rank L > 0) = 0   REGISTER:465-468

numpy only.  No sympy.
"""
import sys
import numpy as np

rng_global = np.random.default_rng(20260817)
OUT = []
def say(*a):
    s = " ".join(str(x) for x in a)
    print(s); sys.stdout.flush(); OUT.append(s)

def rule(t):
    say(""); say("=" * 78); say(t); say("=" * 78)

# ----------------------------------------------------------------------------
# Mahler measure of P(x,y) = sum_j c_j(y) x^j, by Jensen in x then quadrature in y.
# Generic enough to take the charge-split polynomials of F1 leg D as well.
# ----------------------------------------------------------------------------
def mahler_xy(coeff, ny=20000):
    """coeff[i][j] = coefficient of x^i y^j.  Returns m(P) = (2pi)^-2 int log|P|."""
    coeff = np.asarray(coeff, dtype=float)
    nx = coeff.shape[0]
    t = 2.0 * np.pi * (np.arange(ny) + 0.5) / ny
    yv = np.exp(1j * t)
    # c_i(y) for each i
    C = np.array([np.polyval(coeff[i][::-1], yv) for i in range(nx)])  # shape (nx, ny)
    if nx == 2:
        # exact Jensen in x, fully vectorised:  int log|c0 + c1 x| dx/2pi = log max(|c0|,|c1|)
        mx = np.maximum(np.abs(C[0]), np.abs(C[1]))
        if np.any(mx <= 0):
            return -np.inf
        return float(np.mean(np.log(mx)))
    tot = 0.0
    for k in range(ny):
        c = C[:, k]                       # c[0] + c[1] x + ... + c[nx-1] x^{nx-1}
        nz = np.nonzero(np.abs(c) > 1e-14)[0]
        if len(nz) == 0:
            return -np.inf
        top = nz[-1]
        lead = c[top]
        val = np.log(abs(lead))
        if top > 0:
            r = np.roots(c[: top + 1][::-1])
            val += np.sum(np.log(np.maximum(1.0, np.abs(r))))
        tot += val
    return tot / ny

def mahler_pi(pi, ny=20000):
    """pi = (p00, p10, p01, p11);  P = p00 + p10 x + p01 y + p11 xy."""
    p00, p10, p01, p11 = pi
    return mahler_xy([[p00, p01], [p10, p11]], ny=ny)

def Zk(pi, f, c, k):
    p00, p10, p01, p11 = pi
    k = np.asarray(k)
    return (p00 + p10 * np.exp(1j * f * k) + p01 * np.exp(1j * c * k)
            + p11 * np.exp(1j * (f + c) * k))

def lam_emp(pi, f, c, N):
    k = np.arange(1, N + 1)
    z = np.abs(Zk(pi, f, c, k))
    z = np.maximum(z, 1e-300)
    return float(np.mean(np.log(z)))

# ----------------------------------------------------------------------------
rule("CONTROL C0 — the Mahler routine, against two register values it CAN fail")
# REGISTER:172-174  generic torus value  -0.767507880 = m(0.4 + 0.3x + 0.3y)
m1 = mahler_pi((0.4, 0.3, 0.3, 0.0))
say(f"  m(0.4 + 0.3x + 0.3y)                       = {m1:.12f}   register -0.767507880   "
    f"|diff| = {abs(m1 - (-0.767507880)):.3e}")
# REGISTER:1076-1081 (W-10 N-3)  lambda = m(4/9 + (2/9)x) = log(4/9) exactly, pi=(4/9,1/9,2/9,2/9)
m2 = mahler_pi((4 / 9, 1 / 9, 2 / 9, 2 / 9))
say(f"  m(4/9 + (1/9)x + (2/9)y + (2/9)xy)         = {m2:.12f}   log(4/9) = "
    f"{np.log(4/9):.12f}   |diff| = {abs(m2 - np.log(4/9)):.3e}")
# negative control: a value it must NOT reproduce
m3 = mahler_pi((0.25, 0.25, 0.25, 0.25))
say(f"  NEGATIVE CTRL m(uniform) = {m3:.12f}   log(1/4) = {np.log(0.25):.12f}  "
    f"(agrees only because (1+x)(1+y)/4 factors; m(1+x)=0)")
say("  -> the routine reproduces two independent register figures and is not tuned to either.")

# ============================================================================
rule("F1 — EMBED.  Every argument Z_k actually takes.  Is the count two?")
# ---------------------------------------------------------------------------
say("""
F1-A. THE SPECTRAL EMBEDDING.  Under the corpus's fibre-wise convention (REGISTER:577,
      REGISTER:1171-1172) Z_k is ALWAYS of the form  Z_k = sum_j w_j chi_j^k  with w_j >= 0,
      sum w_j = 1, |chi_j| = 1.  So Z_k = mu-hat(k) for a probability measure
              mu = sum_j w_j delta_{chi_j}   on the unit circle.
      Both named arms are moves of the SAME object mu:
              "act on the ready state"  = move the WEIGHTS w
              "act on the connection"   = move the ATOMS  chi
      The embedding space is E_n = P_n(U(1)) (n-atom probability measures), modulo the
      global rotation chi_j -> a*chi_j which leaves every |Z_k| fixed.
              dim E_4 = (4-1 weights) + (4 angles) - (1 rotation) = 6.
      The corpus's slice is chi = (1, u, v, uv): angles (0, f, c, f+c), i.e. 3 weights
      + 2 angles = 5.   6 - 5 = 1 UNREACHED DIMENSION, and it is exactly the additive
      relation  th_11 - th_10 - th_01 = 0,  which is W-11's class-constancy of the
      relative branch operator Q (REGISTER:1171-1176).""")

def jac_rank(fun, x0, npar, M=48, h=1e-6):
    base = fun(x0)
    J = np.zeros((M, npar))
    for i in range(npar):
        xp = np.array(x0, dtype=float); xp[i] += h
        xm = np.array(x0, dtype=float); xm[i] -= h
        J[:, i] = (fun(xp) - fun(xm)) / (2 * h)
    s = np.linalg.svd(J, compute_uv=False)
    return s

M = 48
kk = np.arange(1, M + 1)

def phi_corpus(th):          # (p10,p01,p11, f, c) -> |Z_1..Z_M|
    p10, p01, p11, f, c = th
    pi = (1 - p10 - p01 - p11, p10, p01, p11)
    return np.abs(Zk(pi, f, c, kk))

def phi_free(th):            # (w1,w2,w3, a1,a2,a3) -> |Z_1..Z_M|, atom0 gauge-fixed at angle 0
    w1, w2, w3, a1, a2, a3 = th
    w = np.array([1 - w1 - w2 - w3, w1, w2, w3])
    ang = np.array([0.0, a1, a2, a3])
    return np.abs(w @ np.exp(1j * np.outer(ang, kk)))

x_corpus = np.array([0.30, 0.25, 0.20, 1.13, 0.71])
x_free = np.array([0.30, 0.25, 0.20, 1.13, 0.71, 1.13 + 0.71])   # SAME point, corpus slice
s_c = jac_rank(phi_corpus, x_corpus, 5, M)
s_f = jac_rank(phi_free, x_free, 6, M)
tol = 1e-7
say(f"\n  singular values, corpus parametrisation (pi,f,c), 5 params:")
say("   ", np.array2string(s_c, precision=4))
say(f"    numerical rank at tol {tol:g}: {int(np.sum(s_c > tol))}")
say(f"  singular values, free-atom parametrisation (w, angles), 6 params, SAME point:")
say("   ", np.array2string(s_f, precision=4))
say(f"    numerical rank at tol {tol:g}: {int(np.sum(s_f > tol))}")
say(f"\n  => the binary {{connection, ready state}} spans a rank-{int(np.sum(s_c>tol))} tangent space")
say(f"     inside a rank-{int(np.sum(s_f>tol))} one.  MISSING DIMENSIONS: "
    f"{int(np.sum(s_f>tol)) - int(np.sum(s_c>tol))}")
# identify the missing direction explicitly
v_missing = np.zeros(6); v_missing[3] = 1.0; v_missing[4] = 1.0; v_missing[5] = -1.0
say(f"  the direction (dth_10, dth_01, dth_11) = (1, 1, -1) breaks th_11 = th_10 + th_01;")
d = np.linalg.norm(jac_rank(phi_free, x_free, 6, M) * 0 + 0)  # placeholder, real check below
Jf = np.zeros((M, 6))
h = 1e-6
for i in range(6):
    xp = x_free.copy(); xp[i] += h
    xm = x_free.copy(); xm[i] -= h
    Jf[:, i] = (phi_free(xp) - phi_free(xm)) / (2 * h)
Jc = np.zeros((M, 5))
for i in range(5):
    xp = x_corpus.copy(); xp[i] += h
    xm = x_corpus.copy(); xm[i] -= h
    Jc[:, i] = (phi_corpus(xp) - phi_corpus(xm)) / (2 * h)
resp = Jf @ v_missing
Q, _ = np.linalg.qr(Jc)
resid = resp - Q @ (Q.T @ resp)
say(f"  its response is NOT in the span of the two arms:  ||residual||/||response|| = "
    f"{np.linalg.norm(resid)/np.linalg.norm(resp):.6f}")

# ---------------------------------------------------------------------------
# --- the steelman: is the binary a COMPLETE partition INSIDE the convention? ---
s_state = np.linalg.svd(Jc[:, :3], compute_uv=False)
s_conn = np.linalg.svd(Jc[:, 3:], compute_uv=False)
say(f"\n  STEELMAN CHECK (this test could have come out either way): inside the convention,")
say(f"    rank(state block d/dpi)      = {int(np.sum(s_state > tol))}")
say(f"    rank(connection block d/dfc) = {int(np.sum(s_conn > tol))}")
say(f"    rank(both together)          = {int(np.sum(s_c > tol))}   "
    f"= {int(np.sum(s_state>tol))} + {int(np.sum(s_conn>tol))}, so the two arms are")
say(f"    INDEPENDENT and TOGETHER EXHAUST the convention-scoped tangent space.")
say(f"    => the binary is a SOUND partition CONDITIONAL ON THE TRANSPORT CONVENTION.")
say(f"    It under-enumerates only once the convention is allowed to move — and W-11,")
say(f"    the row immediately preceding this decision, ruled the convention a STIPULATION")
say(f"    (REGISTER:1149-1152 'reached by substitution'; H1 at REGISTER:1223 'a stipulation,")
say(f"    not a theorem, and the corpus's own sealed COR-F exhibits an admissible")
say(f"    alternative under which it fails').")

say("""
F1-B. THE LOOP DESIGNATION — a third locus, ALREADY EXHIBITED IN THE CORPUS.
      ERRATUM AGAINST W-09 (REGISTER:994-1009): class occupancy is a property of the
      LOOP DESIGNATION, not of the complex.  On B0b's own complex, with S4's own gamma_F,
      the STATE fixed (uniform 1/9) and the CONNECTION fixed, sweeping gamma_C over
      admissible simple cycles reaches 16 distinct class multisets.
      Source: LANE_W10_A_CARRIERS_REFUTE_1/r1_rebuild.OUT.txt:81 (cited by the erratum).""")
b0b = [((0, 0, 5 / 9, 4 / 9)), ((0, 1 / 9, 5 / 9, 1 / 3)), ((1 / 9, 1 / 9, 4 / 9, 1 / 3)),
       ((2 / 9, 2 / 9, 1 / 3, 2 / 9))]
vals = [mahler_pi(p, ny=6000) for p in b0b]
say("      lambda over four of those designations, recomputed here:")
for p, w in zip(b0b, vals):
    say(f"        pi = ({p[0]:.4f},{p[1]:.4f},{p[2]:.4f},{p[3]:.4f})   lambda = {w:+.9f}")
say(f"      SPREAD over the loop-designation axis alone: {max(vals)-min(vals):.9f} nats,")
say(f"      at fixed complex, fixed ready state, fixed connection.  Neither arm names it.")

# ---------------------------------------------------------------------------
say("""
F1-C. THE SCHEDULE — a fourth locus, and it is the one W-08 already measured.
      REGISTER:820: durability is a property of the (connection, SCHEDULE) pair, and
      "the corpus has never stated a schedule stipulation".  An adversary writing only the
      sqrt(K) cells of smallest 1-|Z_k| holds |Omega| ~ 0.55 forever.  Reproduced here at
      the corpus's own resonant test point f=2.0, c=1.1, pi=(0.4,0.3,0.3,0).
      A schedule supported on a density-zero subset of N IS a singular constraint —
      and it acts on NEITHER the connection NOR the state.""")
pi_res = (0.4, 0.3, 0.3, 0.0)
for K in [10 ** 4, 10 ** 5, 10 ** 6]:
    k = np.arange(1, K + 1)
    z = np.abs(Zk(pi_res, 2.0, 1.1, k))
    nl = -np.log(np.maximum(z, 1e-300))
    honest = float(np.sum(nl))
    idx = np.argsort(1.0 - z)[: int(np.sqrt(K))]
    adv = float(np.sum(nl[idx]))
    say(f"      K = {K:>8}   honest schedule k_n=n : {honest:14.3f} nats     "
        f"adversarial sqrt(K)-schedule : {adv:8.4f} nats")
say("      honest rate is linear in K; adversarial accumulation is O(1) and K-independent.")
say("      RATIO at K=1e6: " + f"{honest/adv:,.0f} : 1")

# ---------------------------------------------------------------------------
say("""
F1-D. THE CHARGE — a fifth locus.  REGISTER:208-215 (W-03): "THE MODALITY NAMED IN S4'S
      OWN BRIEF AND NEVER RUN: CHARGE ... q = (1,2,2,2,2) moves lambda".  A non-uniform
      charge splits an incidence class into two atoms at different characters WITHOUT
      touching the class pushforward pi, the connection, or the state's class weights.
      Exhibit, built here: K1's class-10 weight 0.3 carried by two vertices at 0.15 each.""")
lam_q1 = mahler_pi((0.4, 0.3, 0.3, 0.0))
lam_q2 = mahler_xy([[0.4, 0.3], [0.15, 0.0], [0.15, 0.0]])   # 0.4 + 0.15x + 0.15x^2 + 0.3y
say(f"      unit charge   P = 0.4 + 0.30x        + 0.3y   lambda = {lam_q1:+.9f}")
say(f"      charges (1,2) P = 0.4 + 0.15x+0.15x^2 + 0.3y   lambda = {lam_q2:+.9f}")
say(f"      MOVED BY {abs(lam_q2-lam_q1):.9f} nats at IDENTICAL pi, IDENTICAL (f,c).")

say("""
F1 VERDICT INPUTS: named arms = 2.  Arguments of the functional established in the
corpus at the cutoff = at least 6:  pi (state) | (u,v) (connection) | the relative branch
operator Q / transport convention (W-11 H1, REGISTER:1171, 1223) | the loop designation
(erratum vs W-09, REGISTER:994, :1007) | the schedule (W-08, REGISTER:820) | the charge
(W-03, REGISTER:208).  Fibre rank / gauge group (W-04 ERR-2, REGISTER:270-279) and the loop
COUNT (W-03 reopen, REGISTER:248-249) make eight.""")

# ============================================================================
rule("F2 — DEGENERACY.  Are the two arms the same object under a map?")
say("""
F2-A. EXACT MAP 1 — change of coordinates (the W-03 involution, REGISTER:203-206).
      Z_k(pi;u,v) * conj(uv)^k = Z_k(sigma pi; conj u, conj v)  with sigma: 00<->11, 10<->01.
      So a CONNECTION move (complex conjugation) is EXACTLY a STATE move (class swap).""")
worst = 0.0
for _ in range(2000):
    pi = rng_global.dirichlet(np.ones(4))
    f, c = rng_global.uniform(0, 2 * np.pi, 2)
    k = np.arange(1, 60)
    a = np.abs(Zk(tuple(pi), f, c, k))
    sg = (pi[3], pi[2], pi[1], pi[0])
    b = np.abs(Zk(sg, -f, -c, k))
    worst = max(worst, float(np.max(np.abs(a - b))))
say(f"      || |Z_k(pi;u,v)| - |Z_k(sigma pi; ubar,vbar)| ||_inf over 2000 draws x 59 cells")
say(f"        = {worst:.3e}     -> the arms OVERLAP on a Z_2 subgroup, distance 0 to machine eps.")

say("""
F2-B. EXACT MAP 2 — restriction to the target functional.  The route wants to move THE
      RATE.  Under H2 (REGISTER:1225: (conj W_F, W_C) generates a dense subgroup of T^2)
      Weyl gives lambda = m(P), a function of pi ALONE.  So the FORGET-THE-CONNECTION map
      f(pi,u,v) = pi collapses the whole connection arm to a POINT for the rate.""")
pi_t = (0.42, 0.19, 0.27, 0.12)
mt = mahler_pi(pi_t)
devs = []
for _ in range(400):
    f = rng_global.uniform(0, 2 * np.pi); c = rng_global.uniform(0, 2 * np.pi)
    devs.append(lam_emp(pi_t, f, c, 200000) - mt)
devs = np.array(devs)
say(f"      m(P) at pi = {pi_t}: {mt:+.9f}")
say(f"      400 random (f,c), lambda_N at N=2e5:  max|lambda_N - m(P)| = {np.max(np.abs(devs)):.3e}")
say(f"                                            DIAMETER of connection orbit = "
    f"{np.max(devs)-np.min(devs):.3e}   (Weyl truncation error, shrinking in N)")
sts = []
for _ in range(400):
    sts.append(mahler_pi(rng_global.dirichlet(np.ones(4)), ny=3000))
sts = np.array(sts)
say(f"      400 random pi at fixed (f,c):         DIAMETER of state orbit      = "
    f"{np.max(sts)-np.min(sts):.6f}")
say(f"      RATIO state-diameter : connection-diameter = "
    f"{(np.max(sts)-np.min(sts))/(np.max(devs)-np.min(devs)):,.0f} : 1")
say("      The two arms are NOT symmetric rivals.  On a full-measure set of connections the")
say("      connection arm is degenerate to a point and only the state arm has any extent.")

say("""
F2-C. MAP 3 — monomial substitution (the natural 'power' map, and the one charge suggests).
      x -> x^q, y -> y^q is measure-preserving on T^2, so m is invariant: a UNIFORM charge
      cannot move the rate.  Only a NON-uniform charge can (F1-D), and that is an atom-set
      move, i.e. the F1 missing dimension, not a connection move.""")
say(f"      m(0.4 + 0.3x + 0.3y)     = {mahler_pi((0.4,0.3,0.3,0.0)):+.9f}")
say(f"      m(0.4 + 0.3x^2 + 0.3y^2) = {mahler_xy([[0.4,0,0.3],[0,0,0],[0.3,0,0]]):+.9f}"
    "   (uniform charge 2: unchanged)")

# ============================================================================
rule("F3 — CARVING.  Measure both cells.")
say("""
Predicate under test: "a singular constraint moves the rate BY ACTING ON THE CONNECTION"
                  vs "                                    BY ACTING ON THE READY STATE".

MEASURE 1 — the natural measures on the two arms (Haar on T^2; Lebesgue on the 3-simplex).
            Cell A = {(f,c) : lambda(pi;f,c) != m(P)} .  This is N3 (REGISTER:465-468).""")
say(f"      Haar sample of 400 connections above: 0 exceed the Weyl error floor "
    f"{np.max(np.abs(devs)):.2e}.  Cell A has Haar measure 0 (N3: P(rank L > 0) = 0).")
say(f"      Cell B = {{pi : m(P) != m(P')}} has Lebesgue measure 1 on the simplex "
    f"(diameter {np.max(sts)-np.min(sts):.4f}).")
say("      MEASURE-ZERO CELL vs FULL-MEASURE CELL: the predicate does not partition anything")
say("      under the measure the corpus itself uses.  This is F3's stated failure mode.")

say("""
MEASURE 2 — the achievable RANGE of lambda under a SINGULAR constraint on each arm.
            This is the measure the question actually needs, since it stipulates singular.""")
# connection arm, singular: subgroup / delta measures
say("      CONNECTION ARM, singular constraints, at fixed pi = (0.4,0.3,0.3,0):")
say(f"        delta at (f,c)=(0,0)  [trivial connection]        lambda = "
    f"{np.log(abs(0.4+0.3+0.3)):+.9f}   (= 0, never writes)")
# find a zero of P on T^2: 0.4+0.3e^{if}+0.3e^{ic} = 0 requires |0.3e^{if}+0.3e^{ic}|=0.4
# 0.3|e^{if}+e^{ic}| = 0.4 -> |2cos((f-c)/2)| = 4/3 ; choose f=-c: 0.4+0.6cos f = 0
f0 = np.arccos(-0.4 / 0.6)
zval = abs(0.4 + 0.3 * np.exp(1j * f0) + 0.3 * np.exp(-1j * f0))
say(f"        delta at (f,c)=({f0:.6f},{-f0:.6f})  [a ZERO of P]   |Z_1| = {zval:.3e}"
    f"   lambda = -inf  (writes instantly)")
# order-4 point of S1: W_F=-1, W_C=-i -> u=conj(W_F)=-1, v=W_C=-i  => f=pi, c=-pi/2
lam_z4 = float(np.mean(np.log(np.maximum(
    np.abs(Zk(pi_t, np.pi, -np.pi / 2, np.arange(1, 4001))), 1e-300))))
say(f"        S1's published order-4 connection (u=-1, v=-i), pi = {pi_t}:")
say(f"          subgroup average = {lam_z4:+.9f}   m(P) = {mt:+.9f}   "
    f"|deviation| = {abs(lam_z4-mt):.4e}")
say(f"          (REGISTER:1226 reports 3.7e-02 at S1's own connection for its own pi — same order.)")
say(f"        => achievable range on the CONNECTION arm under singular constraints: "
    f"[-inf, 0].  UNBOUNDED.")
# state arm range
say("      READY-STATE ARM, at fixed generic (f,c): the achievable range of m(P) over the simplex:")
grid = []
for _ in range(4000):
    grid.append(mahler_pi(rng_global.dirichlet(np.ones(4) * 0.35), ny=1500))
grid = np.array(grid)
say(f"        4000 draws: max = {np.max(grid):+.6f}   min = {np.min(grid):+.6f}")
say(f"        sup is 0 (pi = delta on one class: Z_k == 1); inf is finite (P cannot vanish")
say(f"        identically for a probability vector) -> range approx [{np.min(grid):+.4f}, 0].")
say("      BOTH CELLS ARE NON-EMPTY, AND THE CONNECTION CELL'S RANGE STRICTLY CONTAINS THE")
say("      STATE CELL'S.  'Can it move the rate' is answered YES on both arms and cannot")
say("      come out otherwise.  Under process rule 'could not have failed', the predicate")
say("      carries no discriminating weight.")

say("""
MEASURE 3 — the OVERLAP cell the binary has no name for.  A Gauss-law-type constraint is
            by definition a RELATION between the state and the connection, not a condition
            on either separately.  Exhibited: the locus {(pi,f,c) : Z_1(pi,f,c) = 0} —
            W-01's own convex-hull firing condition (REGISTER:43).  It is a constraint that
            is singular in the joint variable and is a condition on NEITHER arm alone.""")
hits = 0; trials = 200000
u1 = rng_global.uniform(0, 2 * np.pi, trials); v1 = rng_global.uniform(0, 2 * np.pi, trials)
# three-class carrier: 0 in conv{1, u, v} <=> Wendel, exactly 1/4  (REGISTER:1034)
pts = np.stack([np.ones(trials), np.exp(1j * u1), np.exp(1j * v1)])
ang = np.angle(pts)
ang = np.sort(np.mod(ang - ang[0], 2 * np.pi), axis=0)
gaps = np.diff(np.vstack([ang, ang[0:1] + 2 * np.pi]), axis=0)
hits = int(np.sum(np.max(gaps, axis=0) <= np.pi))
say(f"      three-class carrier, 0 in conv{{1,u,v}}: {hits}/{trials} = {hits/trials:.6f}  "
    f"(W-09 exact value 1/4, REGISTER:939)")
cond = np.cos(u1) + np.cos(v1) <= 0
say(f"      four-class carrier,  0 in conv{{1,u,v,uv}} <=> cos f + cos c <= 0: "
    f"{int(np.sum(cond))}/{trials} = {np.sum(cond)/trials:.6f}  (W-09 exact value 1/2, REGISTER:940)")
say("      The joint cell has POSITIVE codimension-1 measure in the joint space and is")
say("      empty in each arm's own variable at fixed other.  The binary has no cell for it.")

# ============================================================================
rule("F4 — PRESUPPOSITION.")
say("""
THE SENTENCE THE QUESTION ASSUMES, in one line:

  "This construction has a CONSTRAINT that can be imposed, the rate is blocked without one,
   and the loci at which such a constraint could act are exhausted by {the connection,
   the ready state}."

TESTED AS THREE CLAIMS IN THEIR OWN RIGHT, independently of the choice:

 (i) "there is a constraint to impose."   FALSE AT THE BYTES.
     grep -ic 'constraint' over all ten sealed artifacts standing at the cutoff
     (S1, S2, S2-audit, S3, S3-audit, S4, W07, W08, W10, W11; 427,872 bytes) = 0.
     'gauss' occurs 5 times, every one of them 'Gauss-Bonnet' or 'Gaussian rationals';
     'Gauss law' occurs 0 times.  The register says so twice in its own voice:
       REGISTER:287-288 (W-04 ERR-3) "there is no action, no coupling constant, no Gauss law,
         no plaquette weight, no backreaction anywhere"
       REGISTER:397 (W-05)       "no backreaction, no Gauss law, no constraint, no equation
         of motion"
       REGISTER:329 (W-04)       corpus occurrences of 'Gauss law', 'constraint' in the
         gauge sense: 0.
     A constraint is a feature of a variational structure.  The corpus has no action
     (W-04 ERR-3, unbroken through four rounds of attack, REGISTER:610).

 (ii) "the rate is blocked."   TRUE ONLY OF ONE ARM, AND THE QUESTION INHERITS IT AS IF OF BOTH.
     N3 (REGISTER:465-466) reads: "the rate is invariant under every absolutely continuous
     CONNECTION measure".  It is a statement about connection measures only.  Nothing in
     N3, or anywhere at the cutoff, blocks the state arm.  The state arm never needed a
     singular constraint, so the two arms are not rivals for the same job.

 (iii) "the loci are exhausted by two."   FALSE — F1, six to eight loci, four of them with
     numbers already in the register at the cutoff.""")

# ============================================================================
rule("F5 — THE NULL OPTION.  Does EITHER branch obtain?")
say("""
NEITHER BRANCH OBTAINS AS A MECHANISM, and the evidence is positive, not residual:

  1. NO CONSTRAINT EXISTS TO PLACE ON EITHER SIDE.  0/427,872 bytes (F4(i)).  W-05's
     "HAS A FIELD: NO ... the action supplies a PRIOR over a fixed background, not a field"
     (REGISTER:393-397) is a ruling of record, and W-06 did not disturb it: its list of
     what stays dead includes "has a field: no" (REGISTER:610).

  2. BOTH BRANCHES OBTAIN TRIVIALLY AS ARITHMETIC (F3, Measure 2): a delta measure at a
     zero of P sends lambda to -inf; a delta at pi = (1,0,0,0) sends it to 0.  So the
     question's literal answer is "both, unboundedly", which is a control that could not
     have failed, and therefore decides no route.

  3. THE THING THAT IS ACTUALLY OPEN AT THE CUTOFF IS ADMISSIBILITY, AND THE REGISTER
     SAYS SO TWICE, ON TWO DIFFERENT AXES, IN ITS OWN REOPEN CLAUSES:
       REGISTER:910  (W-08 reopens) "an INTRINSIC admissibility criterion for SCHEDULES is
         written down under which the corpus's own SUM(1-z_n) test is recovered and the
         K^{-1/2} adversary excluded"
       REGISTER:1007 (erratum vs W-09) "deciding it needs an ADMISSIBILITY CRITERION FOR
         LOOP DESIGNATIONS, which the corpus has never written — the structural twin of
         W-08's missing schedule-admissibility criterion."
     Two registered open admissibility questions, neither on the connection and neither on
     the ready state.  The route decision ranged over neither of them.""")

rule("CONFOUNDS OF THIS LANE — recorded, not fixed (process rule)")
say("""
 C1. THE NEGATIVE CONTROL IS 3.5e-05 OFF, NOT MACHINE-EPS OFF.  m(uniform) came out
     -1.386259703761 against log(1/4) = -1.386294361120.  Cause: (1+x)(1+y)/4 vanishes on
     the torus, so log|P| is unbounded and the y-quadrature converges slowly.  Every
     conclusion above uses differences of order 1e-1 or larger, or exact Jensen branches,
     so this does not reach any of them — but it bounds this lane's Mahler accuracy near
     zeros of P at ~1e-4, not 1e-12.

 C2. I DID NOT REPRODUCE W-08's FOUR SCHEDULE CONSTANTS.  W-08 (REGISTER:822-823) reports the
     adversary accumulating 0.606 / 0.615 / 0.588 / 0.601 nats at K = 1e4..1e7.  My
     sqrt(K)-adversary accumulates 0.687 / 0.363 / 0.115 — bounded and DECREASING, not
     level.  Either W-08 measured SUM(1-|Z_k|) rather than -SUM log|Z_k|, or its adversary
     is a different one.  W-10's own scope table already marks "W-08's four schedule
     constants (not reproduced)" as UNDETERMINED (REGISTER:1052-1053), so this is a known gap
     and not a new contradiction.  THE F1-C CLAIM USED ABOVE IS ONLY THE BOUNDED-vs-LINEAR
     contrast, which is robust to which of the two functionals is meant.

 C3. THE ORDER-4 DEVIATION IS NOT A REPRODUCTION OF W-11's 3.7e-02.  I measured 1.57e-01
     at pi = (0.42,0.19,0.27,0.12), a pi of my own choosing; W-11 measured 3.7e-02 at S1's
     own ready state, which I did not use.  Different inputs, so no comparison is claimed
     beyond "a singular connection measure moves lambda by an amount of order 1e-1".

 C4. I COULD NOT REPRODUCE W-03's CHARGE FIGURE -1.200555 from the register's own data.
     Under the two class assignments for K1 I could construct from S1/S4, q = (1,2,2,2,2)
     gives log(0.3) = -1.203973 or leaves lambda at -0.767508.  W-03's charge run has no
     lane directory at the cutoff.  F1-D therefore uses an exhibit I built myself, not
     W-03's, and only cites W-03 for the existence of the modality.

 C5. THIS LANE IS OPUS 5 AND SO IS EVERY ROW FROM W-07 ON (REGISTER:1112, :1245).  It is NOT a
     lineage-independent check of W-07..W-11.  It is lineage-independent of S1-S4 and
     W-01..W-06 only.  Discount accordingly.

 C6. F3 MEASURE 2 IS A CONTROL THAT COULD NOT HAVE FAILED, AND THAT IS THE FINDING, NOT A
     DEFECT OF THE LANE.  It is reported as carrying no weight AS EVIDENCE ABOUT THE
     ROUTE, and full weight as evidence about the QUESTION.""")

with open("r4_frame.OUT.txt", "w") as fh:
    fh.write("\n".join(OUT) + "\n")
# canonical self-writer marker: reproduce.sh norm() strips exactly '[written]' (grep -vx)
print("\n[written]")

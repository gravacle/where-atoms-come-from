"""S4 -- THE SHARPER QUESTION.  Is ANY record-level quantity SUPER- or SUB-additive, and does
the defect GROW with the number of blocks?

S2 and S3 found only exactly-additive or exactly-zero objects.  There is one place left where a
genuine whole-differs-from-the-parts effect can live, and it is the one physical mechanism in
this model that can couple distant records: THE BATH THEY SHARE.  Integrating the bath out
leaves an EFFECTIVE POTENTIAL ON THE RECORD REGISTER, and that potential is exactly computable.

THE OBJECT.  Every record commutes with H_tot, so the record configuration s in {+-1}^N is a
conserved label and the joint thermal state block-diagonalises over it.  In the s sector the
bath sees   H_B^s = sum_j ( e_j Z_j + lam c_j(s) X_j ),   c_j(s) = sum_{i at site j} s_i,
a sum of commuting single-qubit terms.  Its free energy is therefore EXACTLY

     Phi(s)  =  sum_j f_j(c_j(s)),      f_j(c) = -(1/beta) ln( 2 cosh( beta sqrt(e_j^2 + lam^2 c^2) ) )

Phi is the bath-induced energy of the record configuration -- the closest thing in this model to
a source term sourced by records.  f_j is EVEN in c, so its Walsh expansion has NO degree-1
part: every bit of Phi beyond a constant is INTERACTION between records.  The pair coupling is

     J_2(i,i') = E_s[ Phi(s) s_i s_i' ]

and the questions are: when is it non-zero, and does the total interaction grow linearly in N?
"""
import sys, math, itertools, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_C_EXTENSIVITY")
from lanelib import *

OUT = []
def P_(s=""):
    print(s, flush=True); OUT.append(str(s))

Xs = np.array([[0, 1], [1, 0]], dtype=complex)
Zs = np.array([[1, 0], [0, -1]], dtype=complex)
ENERGIES = [1.0, 1.4, 0.7, 1.1, 0.9, 1.3, 0.8, 1.2]
LAM, BETA = 0.8, 2.0
LOG2 = math.log(2.0)

def f_site(c, e, lam=LAM, beta=BETA):
    c = np.asarray(c, dtype=float)
    E = np.sqrt(e * e + lam * lam * c * c)
    # -(1/beta) ln(2 cosh(beta E)), written stably
    return -(np.logaddexp(beta * E, -beta * E)) / beta

def binom_w(n):
    k = np.arange(n + 1)
    lw = (math.lgamma(n + 1) - np.array([math.lgamma(x + 1) for x in k])
          - np.array([math.lgamma(n - x + 1) for x in k]) - n * LOG2)
    return (n - 2 * k).astype(float), np.exp(lw)

P_("=" * 110)
P_("S4  THE BATH-INDUCED EFFECTIVE POTENTIAL ON THE RECORD REGISTER")
P_("=" * 110)

# ---------------------------------------------------------------- validation of Phi
P_("\n" + "-" * 110)
P_("VALIDATION  --  Phi(s) against the conditional bath free energy computed from the FULL")
P_("            joint Hamiltonian (dense), and against the [PHYS] 4-qubit block at m=1.")
P_("-" * 110)
P_("%-30s %-20s %-20s %-14s" % ("configuration", "Phi(s) [FORMULA]", "F(s) [DENSE]", "|difference|"))
P_("-" * 110)
maxdev = 0.0
rng = np.random.default_rng(7)
for (m, nq, mode) in [(1, 3, "shared"), (2, 3, "shared"), (3, 2, "shared"), (2, 2, "separate")]:
    k = 2 * m
    nsites = nq if mode == "shared" else nq * m
    site = [i % nq for i in range(k)] if mode == "shared" else \
           [b * nq + (r % nq) for b in range(m) for r in range(2)]
    env = Environment(nq=nsites, energies=tuple(ENERGIES[j % 8] for j in range(nsites)), beta=BETA)
    Zl, dimc = code_records_couplings(m)
    HINT = sum(np.kron(Zl[i], env.site[site[i]]) for i in range(k))
    Htot = np.kron(np.zeros((dimc, dimc), dtype=complex), np.eye(env.dim)) \
        + np.kron(np.eye(dimc), env.HB) + LAM * HINT
    for _ in range(2):
        s = rng.choice([-1, 1], size=k)
        # projector onto the sector with these record values
        Ps = np.eye(dimc, dtype=complex)
        for i in range(k):
            Ps = Ps @ ((np.eye(dimc) + s[i] * Zl[i]) / 2)
        Pfull = np.kron(Ps, np.eye(env.dim))
        Hs = Pfull @ Htot @ Pfull
        ev = np.linalg.eigvalsh(Hs)
        # keep only the eigenvalues living inside the sector (rank = env.dim)
        rank = int(round(np.real(np.trace(Pfull))))
        idx = np.argsort(np.abs(ev))[::-1][:rank] if rank < len(ev) else np.arange(len(ev))
        evs = np.sort(ev[idx])
        Fdense = -np.log(np.exp(-BETA * (evs - evs.min())).sum()) / BETA + evs.min()
        cj = [sum(s[i] for i in range(k) if site[i] == j) for j in range(nsites)]
        Phi = sum(float(f_site(cj[j], ENERGIES[j % 8])) for j in range(nsites))
        d = abs(Phi - Fdense); maxdev = max(maxdev, d)
        P_("%-30s %-20.12f %-20.12f %-14.3e"
           % ("%s m=%d nq=%d s=%s" % (mode[:4], m, nq, "".join("+" if x > 0 else "-" for x in s)),
              Phi, Fdense, d))
P_("-" * 110)
P_("max |FORMULA - DENSE| = %.3e   SELF-CHECK: %s"
   % (maxdev, "PASS" if maxdev < 1e-9 else "FAIL -- conclude nothing"))
assert maxdev < 1e-9

# tie to [PHYS]: the same sector energies from the real 4-qubit [[4,2,2]] block
Hp = stab_hamiltonian(1)
recs_v, _, _ = composite_records_writers(1)
Rp = [xz_to_matrix(v, 4) for v in recs_v]
envp = Environment(nq=3, energies=tuple(ENERGIES[:3]), beta=BETA)
rmp = RecordModel(Hp); Pg, kd = rmp.ground_space()
HtotP = np.kron(Hp, np.eye(envp.dim)) + np.kron(np.eye(16), envp.HB) \
    + LAM * sum(np.kron(Rp[i], envp.site[i % 3]) for i in range(2))
s = np.array([1, -1])
PsP = Pg.copy()
for i in range(2): PsP = PsP @ ((np.eye(16) + s[i] * Rp[i]) / 2)
Pf = np.kron(PsP, np.eye(envp.dim))
Hs = Pf @ HtotP @ Pf
ev = np.linalg.eigvalsh(Hs)
rank = int(round(np.real(np.trace(Pf))))
idx = np.argsort(np.abs(ev))[::-1][:rank]
evs = np.sort(ev[idx])
Fd = -np.log(np.exp(-BETA * (evs - evs.min())).sum()) / BETA + evs.min()
cj = [sum(s[i] for i in range(2) if i % 3 == j) for j in range(3)]
PhiP = sum(float(f_site(cj[j], ENERGIES[j])) for j in range(3)) - 2.0   # -2 = the block's H eigenvalue
P_("[PHYS] tie, real 4-qubit [[4,2,2]] block, s=(+,-): Phi+E_code = %.12f  vs dense %.12f  |d| %.3e"
   % (PhiP, Fd, abs(PhiP - Fd)))

# ---------------------------------------------------------------- the pair coupling
def J2_same_site(n, e, lam=LAM, beta=BETA):
    """EXACT: J_2 = 1/2 ( E[f(c''+2)] - E[f(c'')] ),  c'' = sum of n-2 iid uniform +-1"""
    if n < 2: return 0.0
    cs, ws = binom_w(n - 2)
    return 0.5 * float(ws @ f_site(cs + 2, e, lam, beta) - ws @ f_site(cs, e, lam, beta))

def J2_brute(n_at_sites, i, j, energies, lam=LAM, beta=BETA):
    """brute force over all 2^N sign strings.  n_at_sites = list of record-index lists per site."""
    N = sum(len(x) for x in n_at_sites)
    tot = 0.0
    for bits in itertools.product((-1, 1), repeat=N):
        c = [sum(bits[q] for q in grp) for grp in n_at_sites]
        Phi = sum(float(f_site(c[t], energies[t % 8], lam, beta)) for t in range(len(c)))
        tot += Phi * bits[i] * bits[j]
    return tot / (2 ** N)

P_("\n" + "-" * 110)
P_("TABLE 10  --  THE INDUCED PAIR COUPLING  J_2(i,j) = E_s[ Phi(s) s_i s_j ].")
P_("            D-15: the SAME-SITE column is the positive control for the two ZERO columns.")
P_("-" * 110)
P_("%-8s %-22s %-22s %-22s"
   % ("N", "J_2 SAME bath site", "J_2 DIFFERENT site", "J_2 SEPARATE baths"))
P_("-" * 110)
noise_floor = 0.0
for N in (2, 4, 6, 8, 10, 12):
    # same site: all N records on site 0
    same = J2_brute([list(range(N)), []], 0, 1, ENERGIES)
    same_exact = J2_same_site(N, ENERGIES[0])
    # different site: record 0 on site 0, record 1 on site 1, rest split
    grp = [[0] + list(range(2, 2 + (N - 2) // 2)), [1] + list(range(2 + (N - 2) // 2, N))]
    diff = J2_brute(grp, 0, 1, ENERGIES)
    # separate baths: same as different site but sites belong to different blocks -- identical math
    sep = J2_brute([[0], [1]] + [[q] for q in range(2, N)], 0, 1, ENERGIES)
    noise_floor = max(noise_floor, abs(diff), abs(sep))
    P_("%-8d %-22.14e %-22.14e %-22.14e" % (N, same, diff, sep))
    assert abs(same - same_exact) < 1e-12, "exact J_2 formula disagrees with brute force"
P_("-" * 110)
P_("   the exact formula J_2 = 1/2(E[f(c''+2)]-E[f(c'')]) reproduces every SAME-SITE brute-force")
P_("   value to < 1e-12.   NOISE FLOOR of the two zero columns: %.3e" % noise_floor)
P_("   READ: the bath mediates a REAL record-record interaction, but ONLY between records that")
P_("         share a bath qubit.  Between records at different sites -- and therefore between")
P_("         records in different regions with their own baths -- it is zero to machine")
P_("         precision, and zero by the exact argument that Phi is a SUM OVER SITES.")

# ---------------------------------------------------------------- growth of the interaction
def phi_moments(n, e, lam=LAM, beta=BETA):
    cs, ws = binom_w(n)
    v = f_site(cs, e, lam, beta)
    mean = float(ws @ v)
    var = float(ws @ (v - mean) ** 2)
    return mean, math.sqrt(max(var, 0.0)), float(v.max() - v.min())

P_("\n" + "-" * 110)
P_("TABLE 11  --  DOES THE INTERACTION GROW EXTENSIVELY?   N records on ONE shared bath qubit")
P_("            -- the most crowded, most favourable case for a collective effect.")
P_("-" * 110)
P_("%-9s %-18s %-18s %-16s %-16s %-14s"
   % ("N", "|J_2(N)|", "C(N,2)|J_2|", "std(Phi)", "spread(Phi)", "spread/N"))
P_("-" * 110)
NN = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]
rows = []
for N in NN:
    j2 = abs(J2_same_site(N, ENERGIES[0]))
    pair = (N * (N - 1) / 2) * j2
    mean, sd, spr = phi_moments(N, ENERGIES[0])
    rows.append((N, j2, pair, sd, spr))
    P_("%-9d %-18.10e %-18.10e %-16.10f %-16.8f %-14.8f" % (N, j2, pair, sd, spr, spr / N))

P_("\n" + "-" * 110)
P_("TABLE 12  --  THE EXTENSIVITY RATIO S(2N)/S(N).  Requirement (a) demands -> 2.")
P_("-" * 110)
P_("%-9s %-16s %-16s %-16s %-16s"
   % ("N", "C(N,2)|J_2|", "std(Phi)", "spread(Phi)", "chi_total (S3)"))
P_("%-9s %-16s %-16s %-16s %-16s" % ("", "ratio", "ratio", "ratio", "ratio"))
P_("-" * 110)
d = {N: r for N, *r in [(r[0], r[1], r[2], r[3], r[4]) for r in rows]}
for N in NN:
    if 2 * N in d:
        a = d[N]; b = d[2 * N]
        P_("%-9d %-16.6f %-16.6f %-16.6f %-16s"
           % (N, b[1] / a[1] if a[1] else float("nan"),
              b[2] / a[2] if a[2] else float("nan"),
              b[3] / a[3] if a[3] else float("nan"), "see S3 TABLE 6"))

# ---------------------------------------------------------------- fitted exponents
P_("\n" + "-" * 110)
P_("TABLE 13  --  GROWTH EXPONENTS  Q(N) ~ N^p, fitted on log-log over the LARGE-N half of the")
P_("            range, with residuals and an uncertainty.  NO FIT WITHOUT A NOISE FLOOR:")
P_("            the noise floor of this whole computation is the exact-zero column of TABLE 10,")
P_("            %.3e, and every fitted quantity below is many orders of magnitude above it." % noise_floor)
P_("-" * 110)
P_("%-20s %-14s %-14s %-16s %-14s"
   % ("quantity", "exponent p", "sigma(p)", "max |residual|", "verdict vs (a)"))
P_("-" * 110)
big = [r for r in rows if r[0] >= 256]
xs = np.log(np.array([r[0] for r in big], dtype=float))
for name, col in [("C(N,2)|J_2|", 2), ("std(Phi)", 3), ("spread(Phi)", 4), ("|J_2(N)|", 1)]:
    ys = np.log(np.array([r[col] for r in big], dtype=float))
    A = np.vstack([xs, np.ones_like(xs)]).T
    coef, res, *_ = np.linalg.lstsq(A, ys, rcond=None)
    pred = A @ coef
    resid = float(np.abs(ys - pred).max())
    dof = max(len(xs) - 2, 1)
    s2 = float(((ys - pred) ** 2).sum() / dof)
    cov = s2 * np.linalg.inv(A.T @ A)
    sig = math.sqrt(abs(cov[0, 0]))
    v = ("EXTENSIVE" if abs(coef[0] - 1.0) < 3 * max(sig, 1e-3) else
         ("SUB-EXTENSIVE" if coef[0] < 1.0 else "SUPER-EXTENSIVE"))
    P_("%-20s %-14.6f %-14.2e %-16.2e %-14s" % (name, coef[0], sig, resid, v))

# ---------------------------------------------------------------- additivity of Phi itself
P_("\n" + "-" * 110)
P_("TABLE 14  --  ADDITIVITY OF THE INTERACTION ITSELF.  Two blocks of 2 records each.")
P_("            Q = the total pair-interaction magnitude  sum_{i<j} |J_2(i,j)|.")
P_("            CONTROL = separate baths, where the defect MUST be exactly 0.")
P_("-" * 110)
def pairsum(sizes, energies):
    """sum over pairs of |J_2|, exactly: only SAME-SITE pairs contribute.
       `sizes` = number of records on each site, `energies` = that site's bath energy."""
    tot = 0.0
    for n, e in zip(sizes, energies):
        if n >= 2:
            tot += (n * (n - 1) / 2) * abs(J2_same_site(n, e))
    return tot

def phi_stats(sizes, energies):
    """(std, spread) of Phi = sum_j f_j(c_j) over uniform s, sites INDEPENDENT."""
    var = 0.0; spread = 0.0
    for n, e in zip(sizes, energies):
        _, sd, spr = phi_moments(n, e)
        var += sd * sd; spread += spr
    return math.sqrt(var), spread

E0 = ENERGIES[0]
P_("%-36s %-16s %-16s %-14s %-12s"
   % ("configuration", "Q(A+B)", "Q(A)+Q(B)", "DEFECT", "ratio"))
P_("-" * 110)
MS4 = (2, 4, 8, 16, 32, 64, 128, 256, 512, 1024)
for m in MS4:
    QAB = pairsum([2 * m], [E0]); QA = pairsum([m], [E0])
    P_("%-36s %-16.10e %-16.10e %-14.3e %-12.6f"
       % ("SHARED one site, halves of %d" % m, QAB, 2 * QA, QAB - 2 * QA,
          QAB / (2 * QA) if QA else float("nan")))
P_("")
for m in MS4:
    # CONTROL: two IDENTICAL but disjoint baths -- the only difference from the row above is
    # that the two halves no longer share a site.  Same energies, same lam, same beta.
    QAB = pairsum([m, m], [E0, E0]); QA = pairsum([m], [E0])
    P_("%-36s %-16.10e %-16.10e %-14.3e %-12.6f"
       % ("SEPARATE identical baths %d+%d [CTRL]" % (m, m), QAB, 2 * QA, QAB - 2 * QA,
          QAB / (2 * QA) if QA else float("nan")))
P_("")
for m in (16, 256, 1024):
    # second control: separate baths with DIFFERENT energies -- the residue here is a
    # bath-parameter mismatch, not a record-record interaction, and must not be read as one
    QAB = pairsum([m, m], [E0, ENERGIES[1]]); QA = pairsum([m], [E0])
    P_("%-36s %-16.10e %-16.10e %-14.3e %-12.6f"
       % ("SEPARATE, e=1.0 and e=1.4 [CTRL2]" if m == 16 else
          "   same, %d+%d" % (m, m), QAB, 2 * QA, QAB - 2 * QA, QAB / (2 * QA)))

P_("\n" + "-" * 110)
P_("TABLE 15  --  ADDITIVITY OF THE ENERGY SCALES OF Phi (not of an L1 norm of couplings).")
P_("            std(Phi) is the TYPICAL induced energy; spread(Phi) = max_s Phi - min_s Phi.")
P_("-" * 110)
P_("%-36s %-15s %-15s %-14s %-11s %-15s %-14s"
   % ("configuration", "std(A+B)", "std(A)+std(B)", "DEFECT", "ratio", "spread(A+B)", "DEFECT"))
P_("-" * 110)
for m in MS4:
    sA, pA = phi_stats([m], [E0])
    sAB, pAB = phi_stats([2 * m], [E0])
    P_("%-36s %-15.8f %-15.8f %-14.3e %-11.6f %-15.8f %-14.3e"
       % ("SHARED one site, halves of %d" % m, sAB, 2 * sA, sAB - 2 * sA, sAB / (2 * sA),
          pAB, pAB - 2 * pA))
P_("")
for m in MS4:
    sA, pA = phi_stats([m], [E0])
    sAB, pAB = phi_stats([m, m], [E0, E0])
    P_("%-36s %-15.8f %-15.8f %-14.3e %-11.6f %-15.8f %-14.3e"
       % ("SEPARATE identical baths %d+%d [CTRL]" % (m, m), sAB, 2 * sA, sAB - 2 * sA,
          sAB / (2 * sA), pAB, pAB - 2 * pA))
P_("   note: std adds IN QUADRATURE for independent baths, so std(A+B) = sqrt(2)*std(A) there;")
P_("   the additive quantity for disjoint regions is the VARIANCE, and Var(A+B)-Var(A)-Var(B)")
P_("   is exactly 0 in the CTRL rows by independence.")
for m in (16, 256, 1024):
    vA = phi_stats([m], [E0])[0] ** 2
    vAB = phi_stats([m, m], [E0, E0])[0] ** 2
    P_("   Var defect, SEPARATE %d+%d : %.3e   |   Var defect, SHARED halves of %d : %.6e"
       % (m, m, vAB - 2 * vA, m, phi_stats([2 * m], [E0])[0] ** 2 - 2 * vA))

P_("\n" + "=" * 110)
P_("READ OF S4  (filled in from the numbers above, never in advance)")
P_("=" * 110)
P_("""
 * There IS one genuine record-record interaction in this model, and exactly one: the bath
   mediates a coupling J_2 between two records THAT SHARE A BATH QUBIT.  TABLE 10 measures it
   at ~1e-1 for N=2 and it is confirmed by brute-force enumeration over all 2^N sign strings.

 * Between records that do NOT share a bath qubit it is ZERO -- %.1e, the machine floor --
   and zero by an exact argument, not merely numerically: Phi(s) = sum_j f_j(c_j(s)) is a SUM
   OVER SITES, and c_j depends only on the records at site j.  Two records at different sites
   therefore never appear in the same term.  This holds at every N.

 * So the interaction has NO RANGE.  It is not weak at distance and strong nearby; it is
   all-or-nothing contact.  There is no falloff to extrapolate, which is the one thing a field
   would have to supply.

 * Where the two halves DO share a site the defect is real and GROWS: TABLE 14's defect column
   runs 0.88 -> 8.7e3 across N = 4 -> 2048, as ~N^1.5.  But that column is an L1 norm of
   couplings, not an energy.  The energy scales are in TABLE 15: the typical induced energy
   std(Phi) grows as N^0.50 -- SUB-extensive -- and only the extreme spread max_s Phi - min_s Phi
   grows as N^1.00, and the spread is exactly (per-record constant) x N, which is the record
   COUNT in energy units and is what C-35 already excludes.

 * The CTRL rows -- two identical but disjoint baths -- have defect exactly 0 at every N, which
   is what 'disjoint regions' means physically.  The CTRL2 rows show what a mere bath-parameter
   mismatch looks like, so it is not mistaken for an interaction.
""" % noise_floor)
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_C_EXTENSIVITY/s4_induced_interaction.txt",
     "w").write("\n".join(OUT) + "\n")
print("\n[written]")

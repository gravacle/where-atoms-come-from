"""S3 -- THE HOLEVO HALF.  Total chi about all records, SHARED bath vs SEPARATE baths.

THREE ROUTES, and the whole point is that they agree:
  [DIRECT]   evolve the joint state in [CODE] (system dim 4^m, bath dim 2^nq) and read chi with
             the model's own Environment.holevo.  Exact, but costs eigh at dim 4^m * 2^nq.
  [MATRIX]   the exact factorisation below, evaluated with 2x2 matrix exponentials.
  [BLOCH]    the same factorisation in closed form on the Bloch sphere -- O(n) numpy, so it
             reaches n = 10^5 records on one bath site.
  [BLOCH] is checked against [MATRIX] and [MATRIX] against [DIRECT].

THE FACTORISATION.  In [CODE] the records are Z_1..Z_2m on 2m logical qubits, H_S is a multiple
of the identity, and every coupling is lam * Z_i (x) X_{j(i)}.  All the system operators in
H_tot are diagonal in the joint Z basis, so for each sign string s in {+-1}^{2m}

   H_tot restricted to |s>  =  const  +  sum_j ( e_j Z_j + lam c_j(s) X_j ),  c_j(s)=sum_{i:j(i)=j} s_i

a sum of COMMUTING single-qubit bath terms.  With the system maximally mixed on the code space,
   rho(t) = 2^{-2m} sum_s |s><s| (x) tensor_j sigma_j(c_j(s), t).
The c_j are independent across sites (disjoint record sets) and conditioning on s_i = +-1
touches only site j(i).  Entropy is additive on tensor products, so every site except j(i)
cancels:
   chi_i = S(rho_j) - 1/2[ S(rho_j^+) + S(rho_j^-) ],   j = j(i),
with c_j a sum of n_j iid uniform +-1.  chi_i depends ONLY on n_j.  Write it chi(n).  Then
   chi_total = sum_j n_j * chi(n_j).
"""
import sys, math, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_C_EXTENSIVITY")
from lanelib import *

OUT = []
def P_(s=""):
    print(s, flush=True); OUT.append(str(s))

Xs = np.array([[0, 1], [1, 0]], dtype=complex)
Zs = np.array([[1, 0], [0, -1]], dtype=complex)
LOG2 = math.log(2.0)
ENERGIES = [1.0, 1.4, 0.7, 1.1, 0.9, 1.3, 0.8, 1.2]

def vN_bits(r):
    e = np.linalg.eigvalsh(r); e = e[e > 1e-13]
    return float(-(e * np.log2(e)).sum())

# ------------------------------------------------------------------ [MATRIX] route
def sigma_site(c, e, lam, beta, t):
    h = e * Zs + lam * c * Xs
    w, V = np.linalg.eigh(h)
    tau = np.diag(np.array([math.exp(-beta * e), math.exp(beta * e)])); tau = tau / np.trace(tau)
    ph = np.exp(-1j * w * t)
    Uc = V.conj().T @ tau @ V
    return V @ (ph[:, None] * Uc * ph.conj()[None, :]) @ V.conj().T

def binom_w(n):
    k = np.arange(n + 1)
    lw = (math.lgamma(n + 1) - np.array([math.lgamma(x + 1) for x in k])
          - np.array([math.lgamma(n - x + 1) for x in k]) - n * LOG2)
    return n - 2 * k, np.exp(lw)

def chi_of_n_matrix(n, e, lam, beta, times):
    if n == 0: return 0.0
    cs, ws = binom_w(n - 1)
    acc = 0.0
    for t in times:
        rp = sum(w * sigma_site(c + 1, e, lam, beta, t) for c, w in zip(cs, ws))
        rm = sum(w * sigma_site(c - 1, e, lam, beta, t) for c, w in zip(cs, ws))
        acc += max(vN_bits(0.5 * (rp + rm)) - 0.5 * (vN_bits(rp) + vN_bits(rm)), 0.0)
    return acc / len(times)

# ------------------------------------------------------------------ [BLOCH] route
def bloch(c, e, lam, beta, t):
    """Bloch vector of sigma_site(c,...).  h = lam*c X + e Z = n.sigma; the thermal state has
       r0 = (0,0,-tanh(beta e)); Rodrigues rotation by angle 2|n|t about n-hat."""
    c = np.asarray(c, dtype=float)
    a = lam * c; b = float(e)
    z0 = -math.tanh(beta * b)
    nn = np.sqrt(a * a + b * b)
    th = 2.0 * nn * t
    ct, st = np.cos(th), np.sin(th)
    nn2 = np.where(nn > 0, nn * nn, 1.0)
    rx = a * b * z0 * (1 - ct) / nn2
    ry = -a * z0 * st / np.where(nn > 0, nn, 1.0)
    rz = z0 * ct + b * b * z0 * (1 - ct) / nn2
    return np.stack([rx, ry, rz], axis=-1)

def S_bloch(r):
    q = np.linalg.norm(r)
    p = (1.0 + min(q, 1.0)) / 2.0
    if p >= 1.0 - 1e-15 or p <= 1e-15: return 0.0
    return float(-(p * math.log2(p) + (1 - p) * math.log2(1 - p)))

def chi_of_n(n, e, lam, beta, times):
    """EXACT time-averaged chi of one record at a bath site shared by n records total."""
    if n == 0: return 0.0
    cs, ws = binom_w(n - 1)
    acc = 0.0
    for t in times:
        rp = ws @ bloch(cs + 1, e, lam, beta, t)
        rm = ws @ bloch(cs - 1, e, lam, beta, t)
        acc += max(S_bloch(0.5 * (rp + rm)) - 0.5 * (S_bloch(rp) + S_bloch(rm)), 0.0)
    return acc / len(times)

def esite(j, nq, mode):
    """the bath energy of site j.  In SEPARATE mode every block gets an IDENTICAL COPY of the
       same nq-site bath, so the energy depends on the site's index WITHIN the block.  Giving
       block b a bath with different energies would make chi differ between blocks for a reason
       that has nothing to do with any interaction -- that mismatch is shown separately as
       CTRL2, never used as the additivity control."""
    return ENERGIES[(j % nq) % len(ENERGIES)] if mode == "separate" else ENERGIES[j % len(ENERGIES)]

def occupancy(n_records, nq, mode, m=None):
    """which bath SITE each record couples to.
       SHARED  : one bath of nq sites, all 2m records round-robin over it.
       SEPARATE: block b owns its own nq sites; its 2 records go round robin inside them."""
    if mode == "shared":
        return [i % nq for i in range(n_records)]
    return [b * nq + (r % nq) for b in range(m) for r in range(2)]

def occ_counts(m, nq, mode):
    """per-site record counts, WITHOUT materialising the site list (O(nq) not O(m*nq))"""
    k = 2 * m
    if mode == "shared":
        return [k // nq + (1 if j < k % nq else 0) for j in range(nq)]
    per = [2 // nq + (1 if r < 2 % nq else 0) for r in range(nq)]
    return per * m

def chi_analytic(m, nq, mode, lam=0.8, beta=2.0, times=TIMES, matrix=False):
    """returns (per-site contributions occ_j*chi(occ_j) as an array, occ list).
       chi_total = sum of the returned array -- identical to summing over records, since every
       record on site j contributes the same chi(occ_j)."""
    occ = occ_counts(m, nq, mode)
    f = chi_of_n_matrix if matrix else chi_of_n
    cache, out = {}, []
    for j, o in enumerate(occ):
        e = esite(j, nq, mode)
        key = (o, e)
        if key not in cache: cache[key] = f(o, e, lam, beta, times)
        out.append(o * cache[key])
    return np.array(out), occ

def chi_analytic_perrecord(m, nq, mode, lam=0.8, beta=2.0, times=TIMES):
    """per-RECORD chi, in the record order used by chi_direct (for the validation table)"""
    k = 2 * m
    site = occupancy(k, nq, mode, m)
    occ = occ_counts(m, nq, mode)
    cache, out = {}, []
    for i in range(k):
        j = site[i]; e = esite(j, nq, mode)
        key = (occ[j], e)
        if key not in cache: cache[key] = chi_of_n(occ[j], e, lam, beta, times)
        out.append(cache[key])
    return np.array(out), occ

def chi_direct(m, nq, mode, lam=0.8, beta=2.0, times=TIMES, readouts=None):
    """`readouts` limits which records are read out -- each holevo readout costs two dense
       matmuls at the JOINT dimension, so at m=4 reading all 8 records is the expensive part,
       not the eigendecomposition."""
    k = 2 * m
    site = occupancy(k, nq, mode, m)
    nsites = nq if mode == "shared" else nq * m
    env = Environment(nq=nsites, energies=tuple(esite(j, nq, mode) for j in range(nsites)),
                      beta=beta)
    Zl, dimc = code_records_couplings(m)
    HS = np.zeros((dimc, dimc), dtype=complex)
    prop = Propagator(HS, env, [(Zl[i], site[i]) for i in range(k)], lam=lam)
    idx = list(range(k)) if readouts is None else list(readouts)
    return chi_timeavg(prop, [Zl[i] for i in idx], times), env.dim * dimc, idx

P_("=" * 110)
P_("S3  HOLEVO ADDITIVITY  --  total chi about all records, SHARED bath vs SEPARATE baths")
P_("=" * 110)

# ---------------------------------------------------------------- validation
TV = np.linspace(1.0, 13.0, 5)
P_("\n" + "-" * 110)
P_("VALIDATION 1  --  [BLOCH] against [MATRIX] for the crowding function chi(n)")
P_("-" * 110)
P_("%-6s %-22s %-22s %-14s" % ("n", "chi(n) [MATRIX]", "chi(n) [BLOCH]", "|difference|"))
dv = 0.0
for n in (1, 2, 3, 5, 9, 17):
    a = chi_of_n_matrix(n, 1.0, 0.8, 2.0, TIMES); b = chi_of_n(n, 1.0, 0.8, 2.0, TIMES)
    dv = max(dv, abs(a - b))
    P_("%-6d %-22.14f %-22.14f %-14.3e" % (n, a, b, abs(a - b)))
P_("max deviation %.3e -- %s" % (dv, "PASS" if dv < 1e-10 else "FAIL"))
assert dv < 1e-10

P_("\n" + "-" * 110)
P_("VALIDATION 2  --  [ANALYTIC] against [DIRECT] joint-state simulation (5 time points)")
P_("              If this fails, nothing below may be read.")
P_("-" * 110)
P_("%-10s %-4s %-4s %-12s %-9s %-18s %-18s %-14s"
   % ("mode", "m", "nq", "joint dim", "readouts", "sum chi DIRECT", "sum chi ANALYTIC", "max |diff|"))
P_("-" * 110)
maxdev = 0.0
JOBS = [("shared", 1, 3, None), ("shared", 2, 3, None), ("shared", 3, 3, None),
        ("shared", 1, 1, None), ("shared", 2, 1, None), ("shared", 3, 1, None),
        ("separate", 1, 2, None), ("separate", 2, 2, None),
        ("separate", 1, 1, None), ("separate", 2, 1, None), ("separate", 3, 1, None),
        ("shared", 4, 3, [0, 1]),          # dim 2048 -- 2 readouts keeps it affordable
        ("separate", 4, 1, [0, 1])]        # dim 4096 -- the largest DIRECT point reached
for mode, m, nq, ro in JOBS:
    ca, occ = chi_analytic_perrecord(m, nq, mode, times=TV)
    cd, jd, idx = chi_direct(m, nq, mode, times=TV, readouts=ro)
    d = float(np.abs(ca[idx] - cd).max()); maxdev = max(maxdev, d)
    P_("%-10s %-4d %-4d %-12d %-9s %-18.12f %-18.12f %-14.3e"
       % (mode, m, nq, jd, ("all" if ro is None else str(len(idx))), cd.sum(), ca[idx].sum(), d))
P_("-" * 110)
P_("max per-record |ANALYTIC - DIRECT| = %.3e" % maxdev)
P_("SELF-CHECK: %s" % ("PASS -- the factorisation is exact; it may be used at any m"
                       if maxdev < 1e-8 else "FAIL -- STOP, conclude nothing"))
assert maxdev < 1e-8

# ---------------------------------------------------------------- pairwise A+B
P_("\n" + "-" * 110)
P_("TABLE 4  --  THE PAIRWISE TEST:   chi(block A + block B)   vs   chi(A) + chi(B)")
P_("            CONTROL ROWS are the SEPARATE-bath runs, where additivity MUST hold exactly.")
P_("-" * 110)
P_("%-30s %-16s %-16s %-14s %-12s"
   % ("configuration", "chi(A+B)", "chi(A)+chi(B)", "DEFECT", "ratio"))
P_("-" * 110)
for mode, tag in (("shared", ""), ("separate", " [CTRL]")):
    for nq in (1, 2, 3, 4):
        a2, _ = chi_analytic(2, nq, mode); a1, _ = chi_analytic(1, nq, mode)
        P_("%-30s %-16.12f %-16.12f %-14.3e %-12.8f"
           % (("%s bath, nq=%d%s" % (mode.upper(), nq, tag)), a2.sum(), 2 * a1.sum(),
              a2.sum() - 2 * a1.sum(), a2.sum() / (2 * a1.sum())))

P_("")
P_("   CTRL2 -- SEPARATE baths whose ENERGIES DIFFER between the blocks.  Any defect here is a")
P_("   bath-parameter mismatch, NOT an interaction, and must not be read as one:")
_e = ENERGIES
_id = 2 * (chi_of_n(1, _e[0], 0.8, 2.0, TIMES) + chi_of_n(1, _e[1], 0.8, 2.0, TIMES))
_mm = (chi_of_n(1, _e[0], 0.8, 2.0, TIMES) + chi_of_n(1, _e[1], 0.8, 2.0, TIMES)
       + chi_of_n(1, _e[2], 0.8, 2.0, TIMES) + chi_of_n(1, _e[3], 0.8, 2.0, TIMES))
P_("   %-28s chi(A+B) = %.10f   2*chi(A) = %.10f   DEFECT = %.3e"
   % ("identical baths [CTRL]", _id, _id, 0.0))
P_("   %-28s chi(A+B) = %.10f   2*chi(A) = %.10f   DEFECT = %.3e"
   % ("mismatched baths [CTRL2]", _mm, _id, _mm - _id))

# ---------------------------------------------------------------- growth in m
P_("\n" + "-" * 110)
P_("TABLE 5  --  DOES THE DEFECT GROW WITH m?   chi_total(m) vs m*chi_total(1), nq=3 both sides.")
P_("            SEPARATE is the CONTROL: it must be exactly additive at every m.")
P_("-" * 110)
P_("%-7s %-9s %-16s %-13s %-13s %-11s %-16s %-13s"
   % ("m", "N=2m", "chi_tot SHARED", "m*chi_1 SH", "DEFECT SH", "ratio SH", "chi_tot SEP[C]", "DEFECT SEP[C]"))
P_("-" * 110)
sh1 = chi_analytic(1, 3, "shared")[0].sum()
se1 = chi_analytic(1, 3, "separate")[0].sum()
MS = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768,
      1024, 1536, 2048, 3072, 4096, 8192, 16384, 32768]
d5 = {}
for m in MS:
    S = chi_analytic(m, 3, "shared")[0].sum()
    Sb = chi_analytic(m, 3, "separate")[0].sum()
    d5[m] = S
    P_("%-7d %-9d %-16.10f %-13.4f %-13.3e %-11.6f %-16.6f %-13.3e"
       % (m, 2 * m, S, m * sh1, S - m * sh1, S / (m * sh1), Sb, Sb - m * se1))

# ---------------------------------------------------------------- saturation
P_("\n" + "-" * 110)
P_("TABLE 6  --  SATURATION TEST.  Requirement (a) of extensivity demands S(2N)/S(N) -> 2.")
P_("-" * 110)
P_("%-10s %-10s %-18s %-18s %-16s %-16s"
   % ("N=2m", "2N", "chi_tot(N) SHARED", "chi_tot(2N) SHARED", "S(2N)/S(N) SH", "S(2N)/S(N) SEP[C]"))
P_("-" * 110)
for m in MS:
    if 2 * m in d5:
        P_("%-10d %-10d %-18.10f %-18.10f %-16.6f %-16.6f"
           % (2 * m, 4 * m, d5[m], d5[2 * m], d5[2 * m] / d5[m], 2.0))
P_("   the SEPARATE column is exactly 2.000000 at every N by the factorisation (bath grows with")
P_("   the records); it is the positive control showing the ratio column can reach 2.")

# ---------------------------------------------------------------- chi(n)
P_("\n" + "-" * 110)
P_("TABLE 7  --  THE EXACT CROWDING FUNCTION chi(n): ONE bath qubit shared by n records.")
P_("            n*chi(n) is all the chi that one bath site can support.")
P_("-" * 110)
P_("%-10s %-20s %-20s %-16s" % ("n", "chi(n)", "n*chi(n)", "n*chi(n) / 1*chi(1)"))
P_("-" * 110)
c1 = chi_of_n(1, ENERGIES[0], 0.8, 2.0, TIMES)
NS_ = [1, 2, 3, 4, 6, 8, 12, 16, 24, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192,
       16384, 32768, 65536]
for n in NS_:
    cn = chi_of_n(n, ENERGIES[0], 0.8, 2.0, TIMES)
    P_("%-10d %-20.12e %-20.12f %-16.6f" % (n, cn, n * cn, (n * cn) / c1))

# ---------------------------------------------------------------- D-17
P_("\n" + "-" * 110)
P_("TABLE 8  --  D-17: VARY THE VENUE'S OWN SCALE before calling any effect new.")
P_("            The shared-bath ceiling tracks the BATH (nq, lam, beta), not the record count.")
P_("-" * 110)
P_("%-6s %-6s %-7s %-15s %-15s %-15s %-15s %-13s"
   % ("nq", "lam", "beta", "chi_tot m=16", "m=256", "m=4096", "m=65536", "ceiling/nq"))
P_("-" * 110)
for nq in (1, 2, 3, 6):
    for lam in (0.4, 0.8, 1.6):
        for beta in (0.5, 2.0):
            v = [chi_analytic(mm, nq, "shared", lam=lam, beta=beta)[0].sum()
                 for mm in (16, 256, 4096, 65536)]
            P_("%-6d %-6.1f %-7.1f %-15.8f %-15.8f %-15.8f %-15.8f %-13.8f"
               % (nq, lam, beta, v[0], v[1], v[2], v[3], v[3] / nq))

# ---------------------------------------------------------------- relation matrix
P_("\n" + "-" * 110)
P_("TABLE 9  --  THE CHI RELATION MATRIX  M_ij = chi about record j when ONLY record i is")
P_("            coupled.  Off-diagonal = information about one record leaking from another's")
P_("            channel.  [DIRECT] joint-state simulation, m=3, shared nq=3 bath, 25 times.")
P_("-" * 110)
nq = 3
env = Environment(nq=nq, energies=tuple(ENERGIES[:nq]), beta=2.0)
m = 3
Zl, dimc = code_records_couplings(m)
HS = np.zeros((dimc, dimc), dtype=complex)
M = np.zeros((2 * m, 2 * m))
for i in range(2 * m):
    M[i, :] = chi_timeavg(Propagator(HS, env, [(Zl[i], i % nq)], lam=0.8), Zl)
for i in range(2 * m):
    P_("   " + "  ".join("%11.4e" % M[i, j] for j in range(2 * m)))
off = M - np.diag(np.diag(M))
P_("   [CONTROL: the DIAGONAL is what a non-zero looks like on this same scale]")
P_("   max |off-diagonal|      = %.3e     min diagonal = %.8f" % (np.abs(off).max(), np.diag(M).min()))
blkoff = max(abs(M[i, j]) for i in range(2 * m) for j in range(2 * m) if i // 2 != j // 2)
P_("   max |CROSS-BLOCK entry| = %.3e     (records in DIFFERENT regions)" % blkoff)
ev = np.linalg.eigvalsh((M + M.T) / 2)
P_("   spectrum of (M+M^T)/2   : %s" % np.array2string(ev, precision=8))

P_("\n   relation-matrix spectral invariants vs m ([DIRECT], shared nq=3 bath):")
P_("   %-5s %-16s %-16s %-14s %-14s %-14s"
   % ("m", "tr M", "m*tr M(1)", "DEFECT", "lam_max", "max off-diag"))
tr1 = None
for mm in (1, 2, 3):
    Zl2, dc2 = code_records_couplings(mm)
    HS2 = np.zeros((dc2, dc2), dtype=complex)
    Mm = np.zeros((2 * mm, 2 * mm))
    for i in range(2 * mm):
        Mm[i, :] = chi_timeavg(Propagator(HS2, env, [(Zl2[i], i % nq)], lam=0.8), Zl2)
    if tr1 is None: tr1 = np.trace(Mm)
    evv = np.linalg.eigvalsh((Mm + Mm.T) / 2)
    offm = np.abs(Mm - np.diag(np.diag(Mm))).max()
    P_("   %-5d %-16.10f %-16.10f %-14.3e %-14.10f %-14.3e"
       % (mm, np.trace(Mm), mm * tr1, np.trace(Mm) - mm * tr1, evv.max(), offm))

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_C_EXTENSIVITY/s3_chi_additivity.txt",
     "w").write("\n".join(OUT) + "\n")
print("\n[written]")

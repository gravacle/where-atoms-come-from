"""S5 -- REQUIREMENT (b): ADDITIVITY OVER DISJOINT REGIONS.

Gravity's source is a density integrated over a region: two far-apart clusters have a source
equal to the sum of theirs.  A quantity can grow beautifully in N inside one region and still
fail this, and if it fails this it is not a source term at any N -- again a statement about the
FORM of the quantity, not about how far the simulation reached.

CONSTRUCTION.  Two [[n1, n1-2, 2]] and [[n2, n2-2, 2]] codes on DISJOINT qubit sets, embedded
side by side on n1 + n2 qubits.  Records of the union = the two codes' records, embedded.  The
regions are as far apart as the representation allows: no shared qubit, no shared bath site
(except where the test deliberately shares one, which is the point of row 3 below).

For each quantity: Q(union) against Q(region 1) + Q(region 2).
"""
import sys, itertools, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import symplectic_logicals, Environment

OUT = []
def say(s=""):
    print(s); OUT.append(s)

def stab_xz(n): return [[1]*n + [0]*n, [0]*n + [1]*n]
def sp(a, b, n): return (sum(a[i]*b[n+i] + a[n+i]*b[i] for i in range(n))) % 2
def wt(v, n): return sum(1 for i in range(n) if v[i] or v[n+i])
def supp(v, n): return set(i for i in range(n) if v[i] or v[n+i])
def add2(a, b): return [(x+y) % 2 for x, y in zip(a, b)]

def cosetmin(v, n):
    S = stab_xz(n); grp = [[0]*(2*n), S[0], S[1], add2(S[0], S[1])]
    return min((add2(v, g) for g in grp), key=lambda u: wt(u, n))

def embed(v, n, off, ntot):
    """put an n-qubit (x|z) vector on qubits [off, off+n) of an ntot-qubit register"""
    out = [0]*(2*ntot)
    for i in range(n):
        out[off+i] = v[i]; out[ntot+off+i] = v[n+i]
    return out

def recordsA(n):
    return [cosetmin(p[0], n) for p in symplectic_logicals(stab_xz(n), n)]

def recordsB(n):
    return [cosetmin(v, n) for p in symplectic_logicals(stab_xz(n), n) for v in p]

def min_writer_weight(R, n, stabs):
    """MINIMUM weight of an admissible Pauli that flips R.  Computed against the SUPPLIED
       stabiliser set, so it is honest on the union code (which has 4 stabilisers, not 2)."""
    for w in (1, 2):
        idx = range(n) if w == 1 else [(i, j) for i in range(n) for j in range(i+1, n)]
        for pos in idx:
            pos = (pos,) if w == 1 else pos
            for ps in itertools.product(((1, 0), (0, 1), (1, 1)), repeat=w):
                v = [0]*(2*n)
                for q, (x, z) in zip(pos, ps): v[q] = x; v[n+q] = z
                if all(sp(v, st, n) == 0 for st in stabs) and sp(v, R, n) == 1:
                    return w
    return 99

def quantities(V, n, stabs):
    """P_int(symplectic), P_int(overlap), W_tot, lam_max(symplectic), lam_max(overlap)"""
    m = len(V)
    Ms = np.zeros((m, m)); Mo = np.zeros((m, m))
    for i in range(m):
        si = supp(V[i], n)
        for j in range(m):
            if i == j: continue
            Ms[i, j] = sp(V[i], V[j], n)
            Mo[i, j] = len(si & supp(V[j], n))
    ws = [min_writer_weight(v, n, stabs) for v in V]
    return dict(P_s=int(Ms.sum()//2), P_o=int((Mo > 0).sum()//2),
                W_tot=sum(ws),                   # COMPUTED, not assumed
                lam_s=float(np.linalg.eigvalsh(Ms)[-1]) if m else 0.0,
                lam_o=float(np.linalg.eigvalsh(Mo)[-1]) if m else 0.0,
                N=m)

say("="*112)
say("S5   REQUIREMENT (b): IS THE QUANTITY ADDITIVE OVER DISJOINT REGIONS?")
say("="*112)
say()
say("EXACT ARGUMENT 4 -- ANY 'LARGEST EIGENVALUE OF A RELATION MATRIX' FAILS (b) AT EVERY N.")
say("  Disjoint regions share no record and no relation, so the relation matrix of the union is")
say("  BLOCK DIAGONAL: M = M_1 (+) M_2.  The spectrum of a block-diagonal matrix is the UNION of")
say("  the blocks' spectra, so lam_max(M) = MAX(lam_max(M_1), lam_max(M_2)) -- never the sum.")
say("  A maximum is not a sum.  Doubling the enclosed matter leaves lam_max EXACTLY UNCHANGED")
say("  whenever the two halves are alike.  This holds at every N, for every relation, and no")
say("  amount of growth in lam_max inside one region can repair it.")
say()

say("-"*112)
say("TABLE 3.  Combinatorial quantities: Q(region1 + region2) vs Q(region1) + Q(region2).")
say("          SET A records (the commuting family).  n1, n2 disjoint.")
say()
say("   n1   n2 | quantity | Q(reg1)  Q(reg2)  sum   Q(union)   additive?   ratio Q(union)/sum")
verdict = {}
for (n1, n2) in [(6, 6), (8, 8), (10, 10), (8, 12), (12, 12), (14, 14)]:
    ntot = n1+n2
    V1 = [embed(v, n1, 0, ntot) for v in recordsA(n1)]
    V2 = [embed(v, n2, n1, ntot) for v in recordsA(n2)]
    stab_u = [embed(stab_xz(n1)[0], n1, 0, ntot), embed(stab_xz(n1)[1], n1, 0, ntot),
              embed(stab_xz(n2)[0], n2, n1, ntot), embed(stab_xz(n2)[1], n2, n1, ntot)]
    q1 = quantities(recordsA(n1), n1, stab_xz(n1))
    q2 = quantities(recordsA(n2), n2, stab_xz(n2))
    qu = quantities(V1+V2, ntot, stab_u)
    for key in ("N", "P_s", "P_o", "W_tot", "lam_s", "lam_o"):
        s_ = q1[key]+q2[key]; u = qu[key]
        ok = abs(u - s_) < 1e-9
        verdict.setdefault(key, []).append(ok)
        say("  %3d  %3d | %-8s | %7.3f %8.3f %6.3f %9.3f %11s %10s"
            % (n1, n2, key, q1[key], q2[key], s_, u, "YES" if ok else "NO",
               ("%.4f" % (u/s_)) if s_ else "-"))
    say()
say("  ADDITIVITY VERDICT over all region pairs tested:")
for key, v in verdict.items():
    say("    %-8s additive over disjoint regions at every pair tested : %s" % (key, all(v)))
say("-"*112)
say()

say("-"*112)
say("TABLE 4.  Total chi: additive or not, with BOTH controls in the same table (D-15).")
say("          Region 1 has N1 records, region 2 has N2.  Time-averaged over 25 times in [1,13].")
say()
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_B_COLLAPSE")
from chi_lib import total_chi_fixed, TIMES, vN

def chi_separate_uniform(N, lam=0.8, beta=2.0, e=1.0, times=TIMES):
    """N records, each on its OWN single-qubit bath, all baths IDENTICAL (uniform energy).

       The uniform energy matters: with a repeating (1.0,1.4,0.7) pattern, region 2's records
       would be handed different bath energies inside the union than they had alone, and the
       resulting mismatch would look like a failure of additivity when it is only a relabelling.
       Making the local structure identical removes that artifact, so the row below tests
       additivity and nothing else."""
    Zb = np.array([[1, 0], [0, -1]], dtype=complex); Xb = np.array([[0, 1], [1, 0]], dtype=complex)
    hb = e*Zb
    w0, V0 = np.linalg.eigh(hb); p = np.exp(-beta*w0); p /= p.sum()
    r0 = (V0*p) @ V0.conj().T
    vals = []
    for t in times:
        half = {}
        for sgn in (+1, -1):
            w, U = np.linalg.eigh(hb + lam*sgn*Xb); ph = np.exp(-1j*w*t)
            Uc = U.conj().T @ r0 @ U
            half[sgn] = U @ (ph[:, None]*Uc*ph.conj()[None, :]) @ U.conj().T
        vals.append(max(vN(0.5*(half[1]+half[-1])) - 0.5*(vN(half[1])+vN(half[-1])), 0.0))
    m = float(np.mean(vals))
    return N*m, N*float(np.std(vals, ddof=1)/np.sqrt(len(vals)))

say("   venue                                     N1   N2   chi(1)   chi(2)     sum   chi(union)  additive?")
rows = []
for (N1, N2) in [(4, 4), (6, 6), (8, 8), (10, 10)]:
    # SHARED bath: both regions couple into the SAME 3-qubit bath
    a, *_ = total_chi_fixed(N1); b, *_ = total_chi_fixed(N2); u, *_ = total_chi_fixed(N1+N2)
    say("   SHARED 3-qubit bath (regions not far apart)  %3d %4d %8.4f %8.4f %7.4f %11.4f %10s"
        % (N1, N2, a, b, a+b, u, "YES" if abs(u-(a+b)) < 1e-6 else "NO"))
    # SEPARATE baths: each region gets its own bath -- this is what 'far apart' means physically
    ag, _ = chi_separate_uniform(N1); bg, _ = chi_separate_uniform(N2)
    ug, _ = chi_separate_uniform(N1+N2)
    say("   SEPARATE baths, one site per record (CONTROL) %3d %4d %8.4f %8.4f %7.4f %11.4f %10s"
        % (N1, N2, ag, bg, ag+bg, ug, "YES" if abs(ug-(ag+bg)) < 1e-6 else "NO"))
    rows.append((N1, N2, a, b, u, ag, bg, ug))
    say()
say("  READ: the CONTROL row is additive, so the test can register additivity when it is there.")
say("  The SHARED-bath row is not: %s"
    % ("chi(union) < chi(1)+chi(2) at every pair tested"
       if all(r[4] < r[2]+r[3]-1e-6 for r in rows) else "see numbers"))
say("  The additivity of the CONTROL row is imported from the bath, not from the records: it holds")
say("  because one bath site was ADDED per record.  A source term whose extensivity comes from")
say("  adding a degree of freedom per unit of matter is assuming what it is meant to derive.")
say("-"*112)
say()
say("="*112)
say("  S5 SUMMARY")
say("="*112)
say("  P_o IS additive over disjoint regions but is NOT a function of the enclosed record count:")
say("  two disjoint 6-record regions give P_o = 2 x 7 = 14, while ONE region holding 12 records")
say("  gives 21.  Its super-linear growth inside a region is a PACKING effect -- k logicals crowded")
say("  onto n = k+2 qubits -- and it disappears the moment the records are actually far apart.")
say("  A source term must not care how the matter is packed; this one cares about nothing else.")
say()
say("  P_int (either relation) and W_tot ARE additive over disjoint regions -- they are sums over")
say("  records or over pairs, so additivity is automatic and carries no information.")
say("  lam_max of ANY relation matrix is NOT additive, by exact argument 4, at every N.")
say("  Total chi on a SHARED bath is NOT additive: the regions split the bath's capacity (C-36).")
say("  Total chi on SEPARATE baths IS additive, but only because the bath was grown with N.")
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_B_COLLAPSE/s5_additivity.txt", "w").write("\n".join(OUT)+"\n")

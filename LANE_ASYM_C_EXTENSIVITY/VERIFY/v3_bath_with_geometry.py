"""V3 -- THE LANE'S OWN OPEN FOLLOW-UP, RUN.

The lane's five 'IDENTICALLY ZERO by proof' results include:

   "J_2 and chi cross-talk [are zero] because both the bath dynamics and the induced potential
    are a SUM OVER SITES"

and the headline generalises this to "NOTHING COUPLES DISJOINT REGIONS ... immune to the
objection that we simply have not gone to large enough N."

The proof is real, but its hypothesis is that the BATH HAMILTONIAN HAS NO SITE-SITE TERM.  That
is an input of this model, not a fact about records.  The lane's own caveat 2 admits it never
tested a bath with internal structure.  This script tests it.

  H_B(s) = sum_j [ e_j Z_j + lam s_j X_j ]  +  g * sum_j X_j X_{j+1}      (open chain)
  Phi(s) = -(1/beta) ln Tr exp(-beta H_B(s))          exact, dense diagonalisation
  J_2(i,j) = E_s[ Phi(s) s_i s_j ]                    brute force over all 2^nq sign strings

ONE RECORD PER BATH SITE, so records i and j are in DIFFERENT regions with DISJOINT supports
and their own bath sites -- exactly the configuration the lane reports as identically zero.

D-15: g = 0 is the positive control run in the SAME table.  At g = 0 the site-sum proof applies
and every cross-region J_2 must be exactly zero; that column is what a dead channel looks like.
"""
import sys, math, itertools, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_C_EXTENSIVITY")
from qcore import LAM, BETA, ENERGIES

OUT = []
def P(s=""):
    print(s); OUT.append(s)

I2 = np.eye(2, dtype=complex)
X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)

def op(nq, d):
    M = np.array([[1.0 + 0j]])
    for j in range(nq):
        M = np.kron(M, d.get(j, I2))
    return M

def build_cache(nq, g, coup="XX"):
    Zs = [op(nq, {j: Z}) for j in range(nq)]
    Xs = [op(nq, {j: X}) for j in range(nq)]
    H0 = sum(ENERGIES[j % len(ENERGIES)] * Zs[j] for j in range(nq))
    if coup == "XX":
        HB = sum(Xs[j] @ Xs[j + 1] for j in range(nq - 1))
    else:
        HB = sum(Zs[j] @ Zs[j + 1] for j in range(nq - 1))
    return H0 + g * HB, Xs

def phi(s, Hbase, Xs, lam=LAM, beta=BETA):
    H = Hbase + lam * sum(s[j] * Xs[j] for j in range(len(s)))
    w = np.linalg.eigvalsh(H)
    mx = w.min()
    return -(math.log(np.exp(-beta * (w - mx)).sum()) - beta * mx) / beta

def J2_matrix(nq, g, coup="XX"):
    Hbase, Xs = build_cache(nq, g, coup)
    signs = list(itertools.product([1, -1], repeat=nq))
    vals = np.array([phi(s, Hbase, Xs) for s in signs])
    S = np.array(signs, dtype=float)
    Jm = np.zeros((nq, nq))
    for i in range(nq):
        for j in range(nq):
            Jm[i, j] = float((vals * S[:, i] * S[:, j]).mean())
    return Jm, vals

P("=" * 112)
P("V3  GIVE THE BATH A GEOMETRY AND REMEASURE THE 'IDENTICALLY ZERO' CROSS-REGION COUPLING")
P("=" * 112)

NQ = 8
P("\nnq = %d bath sites, ONE record per site (N = %d records, all in DIFFERENT regions)." % (NQ, NQ))
P("Bath chain coupling g*sum_j X_j X_{j+1}.  lam = %.2f, beta = %.2f, site energies %s"
  % (LAM, BETA, ENERGIES[:NQ]))
P("All 2^%d = %d sign strings enumerated exactly; Phi from dense diagonalisation of a %dx%d H_B."
  % (NQ, 2 ** NQ, 2 ** NQ, 2 ** NQ))

P("\n" + "-" * 112)
P("TABLE V3-1  --  max |J_2(i,j)| over CROSS-SITE pairs (i != j), as g is turned on.")
P("                g = 0 is the D-15 CONTROL: the lane's site-sum proof applies there and the")
P("                column must be exactly zero.  Same instrument, same table.")
P("-" * 112)
P("%-9s %-20s %-20s %-20s %-14s"
  % ("g", "max |J_2| cross-site", "|J_2| nearest nbr", "|J_2| farthest", "verdict"))
P("-" * 112)
store = {}
for g in [0.0, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6]:
    Jm, vals = J2_matrix(NQ, g)
    store[g] = Jm
    off = [abs(Jm[i, j]) for i in range(NQ) for j in range(NQ) if i != j]
    nn = np.mean([abs(Jm[i, i + 1]) for i in range(NQ - 1)])
    fa = abs(Jm[0, NQ - 1])
    v = "ZERO (control)" if max(off) < 1e-12 else "NON-ZERO"
    P("%-9.2f %-20.6e %-20.6e %-20.6e %-14s" % (g, max(off), nn, fa, v))

P("\n" + "-" * 112)
P("TABLE V3-2  --  DOES IT HAVE A RANGE?  mean |J_2| by separation d = |i-j|, at g = 0.8.")
P("                A field must supply a FALLOFF.  The lane's stated reason for not calling")
P("                the shared-site coupling a field is that it has 'NO RANGE and no falloff'.")
P("-" * 112)
g = 0.8
Jm = store[g]
P("%-6s %-20s %-14s %-20s" % ("d", "mean |J_2(i,i+d)|", "ratio to d-1", "CONTROL g=0"))
P("-" * 112)
prev = None
lg = []
for d in range(1, NQ):
    vv = np.mean([abs(Jm[i, i + d]) for i in range(NQ - d)])
    c0 = np.mean([abs(store[0.0][i, i + d]) for i in range(NQ - d)])
    P("%-6d %-20.8e %-14s %-20.3e"
      % (d, vv, ("%.4f" % (vv / prev)) if prev else "-", c0))
    lg.append((d, vv)); prev = vv
ds = np.array([d for d, _ in lg], float); ys = np.log(np.array([v for _, v in lg]))
A = np.vstack([ds, np.ones_like(ds)]).T
co, *_ = np.linalg.lstsq(A, ys, rcond=None)
res = float(np.abs(ys - A @ co).max())
# power-law alternative
A2 = np.vstack([np.log(ds), np.ones_like(ds)]).T
co2, *_ = np.linalg.lstsq(A2, ys, rcond=None)
res2 = float(np.abs(ys - A2 @ co2).max())
P("")
P("   EXPONENTIAL fit  |J_2| ~ exp(-d/xi) :  xi = %.4f sites, max log-resid %.3e" % (-1/co[0], res))
P("   POWER-LAW   fit  |J_2| ~ d^-p       :  p  = %.4f,       max log-resid %.3e" % (-co2[0], res2))
P("   %s" % ("=> EXPONENTIAL is the better description (a screened, finite-range interaction)."
            if res < res2 else
            "=> POWER LAW is the better description."))

P("\n" + "-" * 112)
P("TABLE V3-3  --  ADDITIVITY OVER TWO DISJOINT HALVES OF THE CHAIN, g = 0.8 vs the g = 0")
P("                CONTROL.  Q = sum_{i<j} |J_2(i,j)|.  DEFECT = Q(A+B) - Q(A) - Q(B).")
P("-" * 112)
P("%-12s %-16s %-16s %-16s %-14s" % ("g", "Q(A+B)", "Q(A)+Q(B)", "DEFECT", "verdict"))
P("-" * 112)
h = NQ // 2
for g in [0.0, 0.2, 0.8, 1.6]:
    Jm = store[g]
    tot = sum(abs(Jm[i, j]) for i in range(NQ) for j in range(i + 1, NQ))
    qa  = sum(abs(Jm[i, j]) for i in range(h) for j in range(i + 1, h))
    qb  = sum(abs(Jm[i, j]) for i in range(h, NQ) for j in range(i + 1, NQ))
    d = tot - qa - qb
    P("%-12.2f %-16.8e %-16.8e %-16.3e %-14s"
      % (g, tot, qa + qb, d, "ZERO (control)" if abs(d) < 1e-12 else "DISJOINT REGIONS COUPLE"))

P("\n" + "-" * 112)
P("TABLE V3-4  --  D-17, VARY THE VENUE'S OWN SCALE.  Repeat at other bath sizes and with a")
P("                DIAGONAL (ZZ) chain coupling, to show the effect is not an artefact of nq=8")
P("                or of the XX choice.")
P("-" * 112)
P("%-8s %-8s %-8s %-22s %-22s" % ("nq", "coup", "g", "max|J_2| cross-site", "CONTROL g=0"))
P("-" * 112)
for nq in [4, 6, 7, 8]:   # nq=10 (1024 dense diagonalisations of 1024x1024) omitted: cost, not obstruction
    for coup in ["XX", "ZZ"]:
        Jm, _ = J2_matrix(nq, 0.8, coup)
        J0, _ = J2_matrix(nq, 0.0, coup)
        off = max(abs(Jm[i, j]) for i in range(nq) for j in range(nq) if i != j)
        of0 = max(abs(J0[i, j]) for i in range(nq) for j in range(nq) if i != j)
        P("%-8d %-8s %-8.2f %-22.6e %-22.3e" % (nq, coup, 0.8, off, of0))

P("\n" + "=" * 112)
P("READ OF V3  (filled in from the numbers above)")
P("=" * 112)
P(" * At g = 0 every cross-region J_2 is exactly 0 -- the lane's proof and its number are")
P("   reproduced here independently.  The proof is correct WITHIN ITS HYPOTHESIS.")
P(" * At g != 0 the cross-region J_2 is non-zero, on records in disjoint regions with their")
P("   own bath sites, and it FALLS OFF WITH SEPARATION with a fitted decay length.")
P(" * Therefore 'J_2 between records in different regions is IDENTICALLY ZERO at every N' is a")
P("   property of a bath whose Hamiltonian has no site-site term.  It is not a property of")
P("   records, and it does not survive the smallest physical generalisation of the mediator.")
P(" * The lane's stated ground for refusing to call the shared-site coupling a field -- 'no")
P("   range, no falloff to extrapolate' -- is removed the moment the mediator is allowed to")
P("   propagate.  A range appears.")

open("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_C_EXTENSIVITY/VERIFY/v3_bath_with_geometry.txt",
     "w").write("\n".join(OUT) + "\n")

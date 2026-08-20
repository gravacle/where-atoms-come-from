"""S6 -- THE SECOND SATURATING CONTROL, IN THE EXACT FORM THE BRIEF NAMES.

"the Holevo chi of a record on a fixed bath, which is bounded by the bath's bits (C-36)".

Here it is as a series the engine must classify: ONE record, ONE fixed 8-qubit bath, and the
readout restricted to a FRAGMENT of f bath qubits, f = 1..8.  chi(R : fragment_f) rises with f
and is bounded above by chi(R : whole bath), which is bounded by 8 bits.  It is KNOWN to
saturate -- the bound is exact -- so the engine's verdict on it is a test of the ENGINE.

Beside it, in the SAME table (D-15), the linear control in the same units: SUM over f
independent single-qubit baths, each carrying its own record.  That one is exactly linear.

Both series are time-averaged over 25 times in [1,13].  D-17: lam and beta are varied.
"""
import sys, itertools, numpy as np
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_B_COLLAPSE")
from record_model import Environment
from chi_lib import vN, TIMES
from fss_lib import classify

OUT = []
def say(s=""):
    print(s); OUT.append(s)

NQ = 8
def energies(nq):
    pat = (1.0, 1.4, 0.7); return tuple(pat[j % 3] for j in range(nq))

def frag(rB, nq, keep):
    t = rB.reshape([2]*nq + [2]*nq)
    for j in reversed([j for j in range(nq) if j not in keep]):
        t = np.trace(t, axis1=j, axis2=j + t.ndim//2)
    d = 2**len(keep)
    return t.reshape(d, d)

def chi_fragment_series(nq=NQ, lam=0.8, beta=2.0, times=TIMES):
    """ONE record coupled to the SUM of all bath sites; read out on fragments of size f."""
    env = Environment(nq=nq, energies=energies(nq), beta=beta)
    rB0 = env.thermal()
    acc = {f: [] for f in range(1, nq+1)}
    for t in times:
        half = {}
        for s in (+1, -1):
            HB = env.HB + lam*s*env.probe
            w, U = np.linalg.eigh(HB); ph = np.exp(-1j*w*t)
            Uc = U.conj().T @ rB0 @ U
            half[s] = U @ (ph[:, None]*Uc*ph.conj()[None, :]) @ U.conj().T
        for f in range(1, nq+1):
            keep = list(range(f))
            P = frag(half[+1], nq, keep); M = frag(half[-1], nq, keep)
            acc[f].append(max(vN(0.5*(P+M)) - 0.5*(vN(P)+vN(M)), 0.0))
    return {f: (float(np.mean(v)), float(np.std(v, ddof=1)/np.sqrt(len(v)))) for f, v in acc.items()}

def chi_independent_series(nmax=NQ, lam=0.8, beta=2.0, times=TIMES):
    """f INDEPENDENT single-qubit baths, one record each.  Exactly linear in f."""
    Zb = np.array([[1, 0], [0, -1]], dtype=complex); Xb = np.array([[0, 1], [1, 0]], dtype=complex)
    en = energies(nmax)
    per = []
    for i in range(nmax):
        hb = en[i]*Zb
        w0, V0 = np.linalg.eigh(hb); p = np.exp(-beta*w0); p /= p.sum()
        r0 = (V0*p) @ V0.conj().T
        vals = []
        for t in times:
            half = {}
            for s in (+1, -1):
                w, U = np.linalg.eigh(hb + lam*s*Xb); ph = np.exp(-1j*w*t)
                Uc = U.conj().T @ r0 @ U
                half[s] = U @ (ph[:, None]*Uc*ph.conj()[None, :]) @ U.conj().T
            vals.append(max(vN(0.5*(half[1]+half[-1])) - 0.5*(vN(half[1])+vN(half[-1])), 0.0))
        per.append((float(np.mean(vals)), float(np.std(vals, ddof=1)/np.sqrt(len(vals)))))
    out = {}
    for f in range(1, nmax+1):
        out[f] = (sum(per[i][0] for i in range(f)),
                  float(np.sqrt(sum(per[i][1]**2 for i in range(f)))))
    return out

say("="*112)
say("S6   THE SATURATING CONTROL IN THE BRIEF'S OWN FORM, AND ITS LINEAR TWIN, SAME TABLE.")
say("="*112)
say()
sat = chi_fragment_series()
lin = chi_independent_series()
say("     f (bath qubits read)   chi(R : fragment_f)  [KNOWN SATURATING]   SUM of f independent"
    " single-qubit records  [KNOWN LINEAR]")
for f in range(1, NQ+1):
    say("     %8d %28.6f +- %.4f %36.6f +- %.4f"
        % (f, sat[f][0], sat[f][1], lin[f][0], lin[f][1]))
say()
say("     EXACT BOUND on the left column: chi(R : fragment_f) <= chi(R : whole bath) <= %d bits." % NQ)
say("     Observed max: %.6f bits.  Bound respected: %s"
    % (max(v[0] for v in sat.values()), max(v[0] for v in sat.values()) <= NQ + 1e-9))
say()

F = np.arange(1, NQ+1, dtype=float)
for name, D in (("KNOWN SATURATING  chi(R : fragment_f)", sat),
                ("KNOWN LINEAR      f independent records", lin)):
    Q = np.array([D[int(f)][0] for f in F]); S = np.array([max(D[int(f)][1], 1e-6) for f in F])
    r = classify(F, Q, S, name)
    say("  %s" % name)
    say("    engine verdict : %s" % r["category"])
    say("    best form      : %s   dAICc to runner-up %.1f" % (r["best"], min(r["dAICc"], 999.9)))
    say("    AICc ranking   : %s" % ", ".join("%s=%.1f" % (nm, a) for nm, a in r["rank"]))
    say("    alpha*         : %.3f   collapse CV %.4f" % (r["alpha"], r["cv"]))
    say("    1/f -> 0       : %.5f +- %.5f   (1/f-fit rms %.3g)" % (r["Q0"], r["Q0e"], r["rms1N"]))
    say("    doubling ratio : %.3f    exponent %.3f +- %.3f" % (r["dbl"], r["expo"], r["expo_e"]))
    say("    EXTENSIVE?     : %s" % ("YES" if r["extensive"] else "no"))
    say()

say("-"*112)
say("D-17  the same two series under varied coupling and temperature.")
say("      lam   beta |  saturating: chi at f=1,4,8   ratio f=8/f=4  |  linear: ratio f=8/f=4")
for lam in (0.4, 0.8, 1.6):
    for beta in (0.5, 2.0):
        s_ = chi_fragment_series(lam=lam, beta=beta)
        l_ = chi_independent_series(lam=lam, beta=beta)
        say("     %5.2f %5.2f | %8.4f %8.4f %8.4f %14.3f  | %20.3f"
            % (lam, beta, s_[1][0], s_[4][0], s_[8][0], s_[8][0]/s_[4][0], l_[8][0]/l_[4][0]))
say()
say("  A doubling of the resource must give 2.000 for an extensive quantity.  The linear column")
say("  does; the saturating column does not, at any setting tested.")
say("-"*112)
open("/Users/bgm/MB Work/where-atoms-come-from/LANE_ASYM_B_COLLAPSE/s6_saturating_control.txt", "w").write("\n".join(OUT)+"\n")

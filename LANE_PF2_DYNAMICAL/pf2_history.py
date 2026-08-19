"""PF-2: the arrow as a HISTORY, and REDUNDANCY.  Pre-registered at commit bb3d372.

System: toric code 2x2 (dim 256).  Bath: 3 qubits (dim 8), so it has FRAGMENTS.
Coupling lambda * A (x) sum_j X_j.  Initial state: PRODUCT -- system maximally mixed on the
code space, bath thermal.  chi(Zbar:B) is therefore EXACTLY 0 at t=0.

[Zbar,H_tot] = 0 with A = Zbar, so <Zbar> is a constant of motion. Forming a record is not
the record changing value; it is the environment coming to hold information about it."""
import sys, numpy as np
def say(*a): print(*a); sys.stdout.flush()
exec(open('/Users/bgm/MB Work/where-atoms-come-from/LANE_F7_OCCUPANCY/f7_davies.py').read().split('say("="*104); say("0.')[0])
nS = 2**L; nq = 3; nB = 2**nq; beta = 2.0
I2b = np.eye(2); Xb = np.array([[0,1],[1,0]], dtype=complex); Zb = np.array([[1,0],[0,-1]], dtype=complex)
def bop(j, P):
    M = np.array([[1]], dtype=complex)
    for k in range(nq): M = np.kron(M, P if k == j else I2b)
    return M
om = np.array([1.0, 1.4, 0.7])
HB = sum(om[j] * bop(j, Zb) for j in range(nq))
COUP_B = sum(bop(j, Xb) for j in range(nq))
E0, V0 = np.linalg.eigh(H0); gs = int(np.sum(np.abs(E0 - E0[0]) < 1e-9))
Pg = V0[:, :gs] @ V0[:, :gs].conj().T
Ze = op({ind[('h',0,0)]: Z}, L)

def vN(r):
    e = np.linalg.eigvalsh(r); e = e[e > 1e-13]; return float(-(e * np.log2(e)).sum())
def trace_out_system(r):
    return r.reshape(nS, nB, nS, nB).trace(axis1=0, axis2=2)
def frag(rB, keep):
    """partial trace of the 3-qubit bath down to the qubits in `keep`."""
    t = rB.reshape([2]*nq + [2]*nq)
    for j in reversed([j for j in range(nq) if j not in keep]):
        t = np.trace(t, axis1=j, axis2=j + t.ndim // 2)
    d = 2**len(keep); return t.reshape(d, d)
def chi(r, O, keep=None):
    """Holevo information about the +-1 observable O, held by the bath (or a fragment)."""
    out = []
    for s in (+1, -1):
        P = np.kron((np.eye(nS) + s * O) / 2, np.eye(nB))
        blk = P @ r @ P; p = np.real(np.trace(blk))
        if p < 1e-12: continue
        rB = trace_out_system(blk / p)
        out.append((p, frag(rB, keep) if keep is not None else rB))
    if len(out) < 2: return 0.0
    av = sum(p * rb for p, rb in out)
    return max(vN(av) - sum(p * vN(rb) for p, rb in out), 0.0)

def run(A, lam, times):
    Ht = np.kron(H0, np.eye(nB)) + np.kron(np.eye(nS), HB) + lam * np.kron(A, COUP_B)
    w, U = np.linalg.eigh(Ht)
    pB = np.exp(-beta * np.linalg.eigvalsh(HB)); pB /= pB.sum()
    wB, VB = np.linalg.eigh(HB); rB0 = (VB * pB) @ VB.conj().T
    r0 = np.kron(Pg / gs, rB0)
    Uc = U.conj().T @ r0 @ U
    out = []
    for t in times:
        ph = np.exp(-1j * w * t)
        r = U @ (ph[:, None] * Uc * ph.conj()[None, :]) @ U.conj().T
        out.append(r)
    return out

say("="*100); say("PF-2  THE ARROW AS A HISTORY, AND REDUNDANCY"); say("="*100)
say(f"  system dim {nS}, ground degeneracy {gs};  bath = {nq} qubits (dim {nB});  beta = {beta}")
say(f"  SELF-CHECK ||[Zbar,H_S]|| = {np.linalg.norm(Zbar@H0-H0@Zbar):.1e}   "
    f"[Zbar, coupling] = {np.linalg.norm(Zbar@Zbar-Zbar@Zbar):.1e} (A = Zbar commutes with itself)")

times = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0]
say("\n1-2.  THE HISTORY:  chi(Zbar:B)(t) from a PRODUCT state, with <Zbar> conserved")
say(f"  {'t':>7}{'<Zbar>':>12}{'chi(Zbar:B) bits':>20}{'I(S:B) bits':>15}")
rs = run(Zbar, 0.8, times)
for t, r in zip(times, rs):
    rS = r.reshape(nS,nB,nS,nB).trace(axis1=1,axis2=3); rB = trace_out_system(r)
    say(f"  {t:>7.2f}{np.real(np.trace(rS@Zbar)):>12.6f}{chi(r,Zbar):>20.8f}"
        f"{vN(rS)+vN(rB)-vN(r):>15.8f}")

say("\n  CONTROL -- reverse the FULL evolution. Unitary dynamics is exactly reversible, so if the")
say("  instrument reported irreversibility in the CLOSED system it would be manufacturing it.")
back = run(Zbar, 0.8, [4.0, 0.0, -4.0])
say(f"    chi at t=+4.0 : {chi(back[0],Zbar):.8f}")
say(f"    chi at t= 0.0 : {chi(back[1],Zbar):.8f}   <- must be 0")
say(f"    chi at t=-4.0 : {chi(back[2],Zbar):.8f}   <- must equal t=+4 by time-reversal symmetry")

say("\n3-4.  REDUNDANCY (O-13) -- does EACH bath fragment hold a copy?  t = 4.0")
say(f"  {'coupling':<26}{'weight':>8}{'chi whole bath':>17}{'frag {0}':>11}{'frag {1}':>11}{'frag {2}':>11}")
for nm, A, wt in (("Zbar (logical)", Zbar, 2), ("Z_e  (single site)", Ze, 1)):
    r = run(A, 0.8, [4.0])[0]
    cs = [chi(r, Zbar, keep=[j]) for j in range(nq)]
    say(f"  {nm:<26}{wt:>8}{chi(r,Zbar):>17.8f}{cs[0]:>11.6f}{cs[1]:>11.6f}{cs[2]:>11.6f}")

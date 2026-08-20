"""V5b -- corrected independent transport-fixed-record construction (character/isotypic route,
a different code path from the lane's gauge_record, which uses no characters).  My first pass
(v5_fixed_record.py) split raw rank inside an isotypic block, which is wrong whenever d_rho>1;
it failed its own transport test on exactly D_8 and SD_16, the carriers with a 2-dimensional
irrep.  That failure is the self-check working, and is reported rather than hidden."""
import sys, numpy as np
LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_A_TRANSPORT"
sys.path.insert(0, LANE)
import glib
from carriers import census, minimal_torus, eigblocks, perm_apply, generic_record, subset_sums
def say(*a): print(*a); sys.stdout.flush()
say("="*120)
say("V5b  CORRECTED INDEPENDENT TRANSPORT-FIXED RECORD (character/isotypic route)")
say("="*120)
say(f"  {'carrier':<12}{'dim':>6}{'blocks used':>13}{'||R^2-I||':>11}{'||[H,R]||':>11}{'max|Tr P_E R|':>15}"
    f"{'non-scalar':>12}{'max||[A_h,R]||':>16}{'CONTROL generic':>17}")
for name in ("D_4","Q_8","D_8","Q_16","SD_16","M_4(2)","Pauli16"):
    G = next(g for g in glib.ladder(16) if g.name == name)
    H, perms, D = minimal_torus(G); blocks = eigblocks(H)
    cl, chi, d, cls_of = G.chars()
    rng = np.random.default_rng(11)
    M = rng.normal(size=(D,D)); M = (M+M.T)/2
    Mb = sum(M[np.ix_(p,p)] for p in perms)/len(perms)
    R = np.zeros((D,D)); nb = 0
    for ev, Q in blocks:
        P = Q @ Q.conj().T; mdim = Q.shape[1]; mins = []
        for r in range(len(cl)):
            A = np.zeros((D,D), dtype=complex)
            for g in range(G.n): A[perms[g], np.arange(D)] += np.conj(chi[r][cls_of[g]])
            A *= d[r]/G.n
            Pi = P @ A @ P; Pi = (Pi + Pi.conj().T)/2
            w,V = np.linalg.eigh(Pi); K = V[:, np.abs(w-1) < 1e-7]
            if not K.shape[1]: continue
            S = K.conj().T @ Mb @ K; S = (S+S.conj().T)/2
            ww, VV = np.linalg.eigh(S); i = 0
            while i < len(ww):
                j = i
                while j+1 < len(ww) and abs(ww[j+1]-ww[i]) < 1e-7: j += 1
                mins.append(K @ VV[:, i:j+1]); i = j+1
        assert sum(k.shape[1] for k in mins) == mdim, f"{name}: minimal ranks != dim E"
        sel = subset_sums([(k.shape[1],1) for k in mins], mdim//2)
        assert sel is not None, f"{name}: no balanced subset on eigenspace"
        nb += len(mins)
        Pr = sum(mins[i] @ mins[i].conj().T for i,t in enumerate(sel) if t)
        R += np.real(2*Pr - P)
    inv=np.linalg.norm(R@R-np.eye(D)); com=np.linalg.norm(H@R-R@H)
    trm=0.0; ns=False
    for ev,Q in blocks:
        S=Q.conj().T@R@Q; trm=max(trm,abs(complex(np.trace(S))))
        if np.linalg.norm(S-(np.trace(S)/S.shape[0])*np.eye(S.shape[0]))>1e-7: ns=True
    mv=max(np.linalg.norm(perm_apply(p,R)-R) for p in perms)
    Rg=generic_record(blocks,np.random.default_rng(3),D)
    mvg=max(np.linalg.norm(perm_apply(p,Rg)-Rg) for p in perms)
    say(f"  {name:<12}{D:>6}{nb:>13}{inv:>11.1e}{com:>11.1e}{trm:>15.1e}{str(ns):>12}{mv:>16.2e}{mvg:>17.4f}")

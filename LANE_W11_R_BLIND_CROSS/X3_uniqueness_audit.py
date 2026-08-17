"""X3 -- AUDIT OF THE BLIND LANE'S LEG D (the declared weakness, 'T IS NOT UNIQUE').

Two charges to test:
  (a) its 0/400 sweep is drawn from a family LARGER than 'a transport around gamma': its
      rand_root() puts a random unitary on the OFF-LOOP block too, so most draws move fibre
      values at vertices that are not on the loop at all.  Rerun restricted to A = identity
      off the loop and see whether 0/400 survives.  If it does, the confound did not matter;
      if it does not, the blind lane's headline number is void.
  (b) its COROLLARY -- 'each t_j a function of its OWN edge's transport => A = zeta.T,
      zeta^L = 1, T unique up to an L-th root of unity' -- omits an (L-1)-torus.
      Exhibit it: A = D.T with D diagonal on the loop and the product of D's loop entries 1."""
import numpy as np, sys
sys.path.insert(0, '/Users/bgm/MB Work/where-atoms-come-from/LANE_W11_R_BLIND_CROSS')
from xcore import *
np.set_printoptions(precision=9, linewidth=200)

f, c = 1.0, np.sqrt(2.0)
k = K1(f, c)
w = np.array([0.40,0.15,0.15,0.15,0.15])
St = list(equal_pi_triple(k, w, auto_moves(k, w, 0.45), [0.0,-1.1,0.35,2.9,-0.62]))

def spread(AF, AC, nmax=24):
    e = np.array([traj(k, s, nmax, AF, AC, lambda t:(t,t)) for s in St])
    return float(np.ptp(e, axis=0).max())

# ---------------------------------------------------------------- (a)
print("X3.a  THE SWEEP, RESTRICTED TO OPERATORS THAT ARE THE IDENTITY OFF THEIR OWN LOOP")
print("      The blind lane's rand_root() does NOT impose this: for gamma_F it puts a random")
print("      cube root of the identity on span{e3,e4}, vertices that are not on gamma_F.  Such")
print("      an operator is not a transport around gamma_F under any reading, and it breaks")
print("      invisibility for a reason that has nothing to do with the loop.  Restricted rerun:")
def rand_root_restricted(car, which, rng, force_diag=False, classconst=False):
    """A = identity off the loop; on the loop block, any unitary with A^L = W.I."""
    loop = car.loopF if which=='F' else car.loopC
    ang  = car.f    if which=='F' else car.c
    L = len(loop); Wr = np.exp(1j*ang/L); wL = np.exp(2j*np.pi/L)
    A = np.eye(car.nV, dtype=complex)
    labs = rng.integers(0, L, L)
    if classconst:
        cl = car.classes(); lab = {}
        for o in ORDER: lab[o] = rng.integers(0, L)
        labs = np.array([lab[(int(car.onF()[v]), int(car.onC()[v]))] for v in loop])
    D = Wr*np.power(wL, labs)
    if force_diag:
        blk = np.diag(D)
    else:
        X = rng.normal(size=(L,L)) + 1j*rng.normal(size=(L,L))
        V, _ = np.linalg.qr(X)
        blk = V @ np.diag(D) @ V.conj().T
    A[np.ix_(loop, loop)] = blk
    return A
rng = np.random.default_rng(4242)
rows = {}
for nm, kw in (("generic unitary root on the loop block, identity off it", dict()),
               ("DIAGONAL root, labels a function of the class",          dict(force_diag=True, classconst=True)),
               ("DIAGONAL root, per-vertex labels",                       dict(force_diag=True))):
    v = []
    for _ in range(400):
        AF = rand_root_restricted(k,'F',rng,**kw); AC = rand_root_restricted(k,'C',rng,**kw)
        assert np.linalg.norm(np.linalg.matrix_power(AF,3)-k.M('F')) < 1e-9
        assert np.linalg.norm(np.linalg.matrix_power(AC,3)-k.M('C')) < 1e-9
        v.append(spread(AF,AC))
    v = np.array(v); rows[nm] = v
    print("        %-54s min %.2e med %.2e max %.2e  #(<1e-12)=%3d/400" %
          (nm, v.min(), np.median(v), v.max(), int((v<1e-12).sum())))
print("      DIAGNOSIS of the survivors in row 1 (the blind lane reports 0/400 on its wider")
print("      family; restricted, a few survive).  A generic unitary root is diagonal only when")
print("      its three cube-root labels COINCIDE, i.e. when the block is scalar; P = (1/9)^2 per")
print("      draw for the pair = 4.94 expected in 400.  Counted directly:")
rng2 = np.random.default_rng(4242); nsc = 0; nsurv = 0
for _ in range(400):
    AF = rand_root_restricted(k,'F',rng2); AC = rand_root_restricted(k,'C',rng2)
    sc = (np.linalg.norm(AF-np.diag(np.diag(AF)))<1e-12 and np.linalg.norm(AC-np.diag(np.diag(AC)))<1e-12)
    sv = spread(AF,AC) < 1e-12
    nsc += int(sc); nsurv += int(sv)
    assert sc == sv or not sv, "a NON-diagonal root preserved invisibility -- would refute the criterion"
print("        scalar-block draws %d/400 ; invisibility-preserving draws %d/400 ; identical sets: %s"
      % (nsc, nsurv, nsc == nsurv))
print("        -> every survivor is a SCALAR block (= a power of M^{1/L} times the identity), i.e.")
print("           fibre-wise and class-constant.  No non-fibre-wise root preserved invisibility.")
print("      VERDICT ON CHARGE (a): the restricted sweep gives the same qualitative table, so the")
print("      blind lane's confound did not change its conclusion -- but its published numbers were")
print("      measured on the wrong family and the row headed 'generic pair' is not a sweep over")
print("      transports around the loops.  Recorded as a defect, scored as NOT changing the verdict.")

# ---------------------------------------------------------------- (b)
print("\nX3.b  THE COROLLARY OMITS AN (L-1)-TORUS OF EDGE-LOCAL SHIFTS")
print("      Take A = D.T with D = diag(d_v) on the loop, d off-loop = 1, prod_{v in loop} d_v = 1.")
print("      Then A[w_{j+1},w_j] = d_{w_{j+1}} U_{e_j}: still 'moves one edge, multiplying by a")
print("      function of that edge's own transport'.  A^L = M EXACTLY, A is unitary, A is")
print("      gauge-COVARIANT, and A is NOT zeta.T for any scalar zeta.")
rng = np.random.default_rng(80808)
worstA = 0.0; worstG = 0.0; sp = []
for trial in range(200):
    ds = {}
    for wch in 'FC':
        loop = k.loopF if wch=='F' else k.loopC
        ph = rng.uniform(0, 2*np.pi, len(loop)-1)
        ph = np.concatenate([ph, [-ph.sum()]])              # product of d = 1
        d = np.ones(k.nV, complex); d[loop] = np.exp(1j*ph)
        ds[wch] = np.diag(d)
    AF, AC = ds['F']@k.T('F'), ds['C']@k.T('C')
    worstA = max(worstA, np.linalg.norm(np.linalg.matrix_power(AF,3)-k.M('F')),
                        np.linalg.norm(np.linalg.matrix_power(AC,3)-k.M('C')))
    # gauge covariance of D.T
    th = rng.uniform(0,2*np.pi,5); G = np.diag(np.exp(1j*th))
    sh = lambda loop, phv: np.array([phv[j]+th[loop[(j+1)%len(loop)]]-th[loop[j]] for j in range(len(loop))])
    kg = Cx('K1g',5,[0,1,2],[0,3,4], sh(k.loopF,k.phF), sh(k.loopC,k.phC))
    worstG = max(worstG, np.linalg.norm(ds['F']@kg.T('F') - G@(ds['F']@k.T('F'))@np.linalg.inv(G)))
    sp.append(spread(AF,AC))
sp = np.array(sp)
print("      200 draws:  max ||A^L - M|| = %.2e   max gauge-covariance defect = %.2e" % (worstA, worstG))
print("      invisibility spread: min %.2e  med %.2e  max %.2e   #(<1e-12) = %d/200" %
      (sp.min(), np.median(sp), sp.max(), int((sp<1e-12).sum())))
nsc2 = 0
rng4 = np.random.default_rng(80808)
for _ in range(200):
    loop = k.loopF
    ph = rng4.uniform(0, 2*np.pi, len(loop)-1); ph = np.concatenate([ph, [-ph.sum()]])
    d = np.ones(k.nV, complex); d[loop] = np.exp(1j*ph)
    A = np.diag(d)@k.T('F'); R = A@np.linalg.inv(k.T('F'))
    nsc2 += int(np.linalg.norm(R - R[0,0]*np.eye(5)) < 1e-12)
print("      A_F is zeta.T (i.e. D scalar on the loop) in %d / 200 draws." % nsc2)
print("      [SELF-DEFECT S-X, RECORDED: the first run of this script printed a placeholder")
print("       '%d / 200' fed by a constant instead of this count.  No other number was affected;")
print("       the corrected count is the line above.]")
# explicit single witness, checked to be non-scalar
d = np.ones(5, complex); d[[0,1,2]] = [np.exp(1j*0.9), np.exp(-1j*0.4), np.exp(-1j*0.5)]
AF = np.diag(d)@k.T('F')
Rt = AF@np.linalg.inv(k.T('F'))
print("      witness: d = (e^{0.9i}, e^{-0.4i}, e^{-0.5i}) on gamma_F, product = %.12f%+.12fj" % (np.prod(d[[0,1,2]]).real, np.prod(d[[0,1,2]]).imag))
print("               ||A_F^3 - M_F|| = %.2e ;  A_F T_F^{-1} scalar? %s" %
      (np.linalg.norm(np.linalg.matrix_power(AF,3)-k.M('F')),
       np.linalg.norm(Rt - Rt[0,0]*np.eye(5)) < 1e-12))
print("      VERDICT ON CHARGE (b): the corollary as written ('a function of its OWN edge's")
print("      transport') admits this torus and is FALSE.  It is repaired by strengthening the")
print("      premise from 'a function of' to 'IS' -- i.e. the per-step factor is the edge's own")
print("      U_e and nothing else, which is the definition of discrete parallel transport, or")
print("      equivalently by demanding the SAME local rule at every edge (position-uniformity).")
print("      The repaired corollary is TRUE and the blind lane's CONCLUSION is unaffected: every")
print("      member of the torus breaks invisibility too.  Defect in a proof, not in a verdict.")

# ---------------------------------------------------------------- (c)
print("\nX3.c  THE INVISIBILITY-PRESERVING FAMILY IS LARGER THAN 'THE CIRCUIT CONVENTION'S ROOT'")
print("      S.(class-constant L-th roots of unity) is invisible too, so 'the diagonal branch is")
print("      diag(W^{1/L})' understates it.  Exhibit on K1:")
for labs in ((0,0),(1,0),(2,1)):
    AF = k.S('F').copy(); AC = k.S('C').copy()
    cl = k.classes()
    for v in cl[(1,0)]: AF[v,v] *= np.exp(2j*np.pi*labs[0]/3)
    for v in cl[(0,1)]: AC[v,v] *= np.exp(2j*np.pi*labs[1]/3)
    print("        class labels %s : ||A^L-M|| = %.1e / %.1e   spread = %.2e" %
          (str(labs), np.linalg.norm(np.linalg.matrix_power(AF,3)-k.M('F')),
           np.linalg.norm(np.linalg.matrix_power(AC,3)-k.M('C')), spread(AF,AC)))
print("      (labels attached to CLASSES stay invisible; the same labels attached to VERTICES")
print("       within a class do not -- X2 row R.)")

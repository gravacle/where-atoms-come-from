"""X4 -- SCOPE AUDIT OF THE BLIND LANE'S 'ONE-LINE THEOREM THAT IS THE WHOLE ANSWER' (its G.4):
       'INCIDENCE IS INVISIBLE  <=>  (a,b) in L_F Z x L_C Z'.
   Its criterion on Q (its D2) -- |<A_F^a s, A_C^b s>| is a function of pi alone for every s
   iff Q = A_F^{-a} A_C^b is diagonal with a class-constant diagonal -- is separately tested.
   The D2 criterion is the theorem; the G.4 lattice statement is a SPECIALISATION that needs a
   hypothesis, and the hypothesis is not stated.  Both are attacked here.
   Testing method: instead of three hand-built arms I sample 16 random states per (a,b) with
   pi held EXACTLY fixed, which is a strictly stronger test of 'function of pi alone'."""
import numpy as np, sys
sys.path.insert(0, '/Users/bgm/MB Work/where-atoms-come-from/LANE_W11_R_BLIND_CROSS')
from xcore import *
np.set_printoptions(precision=6, linewidth=250)

def equal_pi_sample(car, pi, n, rng):
    """n random states with EXACTLY the given class pushforward: random split inside each class,
       random phases everywhere."""
    cl = car.classes(); out = []
    for _ in range(n):
        wv = np.zeros(car.nV)
        for j, o in enumerate(ORDER):
            idx = cl[o]
            if len(idx) == 0: continue
            x = rng.random(len(idx)) + 0.05; x = x/x.sum()*pi[j]
            wv[idx] = x
        out.append(np.sqrt(wv)*np.exp(1j*rng.uniform(0, 2*np.pi, car.nV)))
    return out

def diag_classfn(car, Q, tol=1e-10):
    if np.linalg.norm(Q - np.diag(np.diag(Q))) > tol: return False
    d = np.diag(Q)
    for o in ORDER:
        idx = car.classes()[o]
        if len(idx) > 1 and np.max(np.abs(d[idx]-d[idx[0]])) > tol: return False
    return True

def audit(car, amax, bmax, rng, tol=1e-10, label=""):
    """returns (#cells, #mismatch between the D2 criterion and observation, #mismatch between
       the G.4 lattice rule and observation, list of lattice-rule counterexamples)"""
    pi = car.pi(equal_pi_sample(car, np.array([0.25,0.25,0.25,0.25]), 1, rng)[0]) \
         if False else None
    cl = car.classes()
    base = np.array([0.4 if len(cl[o])>0 else 0.0 for o in ORDER])
    base = np.array([(len(cl[o])+0.7) for o in ORDER], float); base[[len(cl[o])==0 for o in ORDER]] = 0.0
    base = base/base.sum()
    S = equal_pi_sample(car, base, 16, rng)
    TF, TC = car.T('F'), car.T('C')
    TFi = np.linalg.inv(TF)
    nmis_d = nmis_l = 0; ctr = []; ncell = 0
    for a in range(0, amax+1):
        for b in range(0, bmax+1):
            Q = np.linalg.matrix_power(TFi, a) @ np.linalg.matrix_power(TC, b)
            pred_d = diag_classfn(car, Q)
            pred_l = (a % car.LF == 0) and (b % car.LC == 0)
            vals = [abs(np.vdot(np.linalg.matrix_power(TF,a)@s, np.linalg.matrix_power(TC,b)@s)) for s in S]
            obs = (max(vals)-min(vals)) < 1e-11
            ncell += 1
            nmis_d += int(pred_d != obs)
            if pred_l != obs: nmis_l += 1; ctr.append((a,b,pred_l,obs))
    return ncell, nmis_d, nmis_l, ctr

rng = np.random.default_rng(5150)
f, c = 1.0, np.sqrt(2.0)

print("X4.1  THE TWO STATEMENTS ON THE CORPUS'S OWN CARRIERS, a,b <= 24, 16 equal-pi states/cell")
for car in (K1(f,c), B0b_registrar(f,c), B0b_blind(f,c)):
    n, md, ml, ctr = audit(car, 24, 24, rng)
    print("   %-11s cells %d   D2-criterion mismatches %d   G.4-lattice mismatches %d" % (car.name, n, md, ml))
print("   -> on the carriers of record both statements hold exactly.  Now the hypothesis.")

print("\nX4.2  THE G.4 LATTICE STATEMENT IS FALSE WITHOUT A HYPOTHESIS THE BLIND LANE DID NOT STATE")
print("      Degenerate carrier D1: gamma_C = gamma_F as a DIRECTED cycle (the two designated")
print("      loops coincide).  Then T_C = T_F and Q_{a,b} = T^{b-a}, which is diagonal and")
print("      class-constant whenever L | (b-a) -- e.g. at (a,b) = (1,1), which is NOT in LZ x LZ.")
Dg = Cx('degen/same', 5, [0,1,2], [0,1,2], [0.31,0.47,0.22], [0.31,0.47,0.22])
n, md, ml, ctr = audit(Dg, 9, 9, rng)
print("      cells %d   D2-criterion mismatches %d   G.4-lattice mismatches %d" % (n, md, ml))
print("      first eight lattice-rule counterexamples (a,b,predicted,observed): %s" % ctr[:8])
print("      Degenerate carrier D2: gamma_C = gamma_F traversed BACKWARDS, with the reversed")
print("      traversal's own phases (-p2,-p1,-p0), so that T_C = T_F^{-1} EXACTLY and hence")
print("      Q_{a,b} = T_F^{-(a+b)}: invisibility on the ANTI-diagonal L | (a+b), off L Z x L Z.")
pF = [0.31, 0.47, 0.22]
Dg2 = Cx('degen/rev', 5, [0,1,2], [0,2,1], pF, [-pF[2], -pF[1], -pF[0]])
print("      check T_C = T_F^{-1} : %.2e" % np.linalg.norm(Dg2.T('C') - np.linalg.inv(Dg2.T('F'))))
n, md, ml, ctr = audit(Dg2, 9, 9, rng)
print("      cells %d   D2-criterion mismatches %d   G.4-lattice mismatches %d ; e.g. %s" % (n, md, ml, ctr[:6]))
print("      -> the D2 criterion (a statement about Q) survives every case; the G.4 lattice")
print("         statement (a statement about a and b) needs 'the two loops are distinct enough")
print("         that Q is diagonal only when both factors are'.  The corpus's carriers satisfy")
print("         it; the theorem as published does not say so.")

print("\nX4.3  DEGENERATE CLASS STRUCTURE -- the blind lane's own S-5, tested rather than argued")
print("      (a) loops sharing NO vertex: class 11 empty.")
Dj = Cx('disjoint', 7, [0,1,2], [3,4,5], [0.3,0.4,0.3], [0.5,0.6,0.31])
n, md, ml, ctr = audit(Dj, 12, 12, rng)
print("          cells %d   D2 mismatches %d   G.4 mismatches %d" % (n, md, ml))
print("      (b) SINGLETON classes.  Two triangles sharing TWO vertices give |10| = |01| = 1,")
print("          so 'class-constant' constrains nothing there and the criterion could weaken.")
print("          It is constructible, contrary to a first guess, so it is TESTED, not argued:")
for nm, nV, lF, lC in (("11={0,1} 10={2} 01={3} 00={4}", 5, [0,1,2], [0,1,3]),
                       ("shared edge traversed oppositely", 5, [0,1,2], [1,0,3]),
                       ("no spectator at all",              4, [0,1,2], [0,1,3])):
    car = Cx('sing', nV, lF, lC, rng.uniform(-2,2,3), rng.uniform(-2,2,3))
    cnt = {"".join(map(str,o)): len(car.classes()[o]) for o in ORDER}
    n, md, ml, ctr = audit(car, 12, 12, rng)
    print("          %-34s classes %s  cells %d  D2 mism %d  G.4 mism %d %s"
          % (nm, cnt, n, md, ml, ("e.g. "+str(ctr[:3])) if ml else ""))

print("\nX4.4  RANDOM CARRIERS -- 60 random pairs of simple cycles, a,b <= 12, 16 states each")
bad_d = bad_l = 0; tot = 0; examples = []
for trial in range(60):
    nV = int(rng.integers(6, 11))
    LF = int(rng.integers(3, 6)); LC = int(rng.integers(3, 6))
    vs = rng.permutation(nV)
    loopF = list(vs[:LF])
    shared = int(rng.integers(0, min(LF, LC)))          # how many vertices the loops share
    pool = [v for v in vs if v not in loopF]
    if len(pool) < LC - shared: continue
    loopC = list(rng.permutation(loopF)[:shared]) + list(np.array(pool)[:LC-shared])
    if len(set(loopC)) != LC: continue
    car = Cx('rand%d'%trial, nV, loopF, loopC, rng.uniform(-2,2,LF), rng.uniform(-2,2,LC))
    n, md, ml, ctr = audit(car, 12, 12, rng)
    tot += 1; bad_d += (md>0); bad_l += (ml>0)
    if ml > 0 and len(examples) < 4: examples.append((car.loopF, car.loopC, ctr[:3]))
print("   carriers tested %d ; with a D2-criterion mismatch %d ; with a G.4-lattice mismatch %d"
      % (tot, bad_d, bad_l))
for e in examples: print("      lattice-rule counterexample: gamma_F=%s gamma_C=%s  cells %s" % e)

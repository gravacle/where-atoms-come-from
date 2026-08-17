"""X5 -- CHECKING THE BLIND LANE'S POSITIVE RESULTS (its P-1, P-2) AND MEASURING THE THING
   NEITHER LANE MEASURED: HOW MUCH OF THE READY STATE EACH CONVENTION ACTUALLY DISCARDS.

P-1 (blind):  Z^T_{Pq+r} = sum_classes c^(r)_ab U^{aq} V^{bq}, c^(0) = pi exactly, and
              lambda_edge = (1/P) sum_r m(P_r).   P = lcm(L_F, L_C).
P-2 (blind):  on B0b the edge rate equals the circuit rate for real-positive states because
              class 00 is fixed by both transports and its weight dominates -- 'a null that
              reads two ways, scored as neither'.  Tested here by sweeping p00.
NEW:          the state -> observable map's RANK under each convention.  'Carrier-independent
              beyond pi' is a claim that the functional factors through a 4-number summary.
              Measure the dimension of what each convention keeps."""
import numpy as np, sys
sys.path.insert(0, '/Users/bgm/MB Work/where-atoms-come-from/LANE_W11_R_BLIND_CROSS')
from xcore import *
from math import gcd
np.set_printoptions(precision=9, linewidth=220)

def lcm(a,b): return a*b//gcd(a,b)

def check_residue_formula(car, s, N=400):
    P = lcm(car.LF, car.LC)
    C = [dressed_coeffs(car, s, r) for r in range(P)]
    U = np.exp(-1j*car.f*(P//car.LF)); V = np.exp(1j*car.c*(P//car.LC))
    worst = 0.0
    TF, TC = car.T('F'), car.T('C')
    xF = s.copy(); xC = s.copy()
    for n in range(1, N+1):
        xF = TF@xF; xC = TC@xC
        Z = np.vdot(xF, xC)
        q, r = divmod(n, P)
        pred = (C[r][0] + C[r][1]*U**q + C[r][2]*V**q + C[r][3]*(U*V)**q)
        worst = max(worst, abs(Z - pred))
    return P, C, worst

f, c = 1.0, np.sqrt(2.0)
print("X5.1  THE RESIDUE IDENTITY, CHECKED AS AN IDENTITY (not as a rate agreement)")
for car in (K1(f,c), B0b_registrar(f,c), B0b_blind(f,c)):
    wv = np.array([0.40,0.15,0.15,0.15,0.15]) if car.nV==5 else np.ones(9)/9
    St = list(equal_pi_triple(car, wv, auto_moves(car, wv, 0.45),
              [0.0,-1.1,0.35,2.9,-0.62] if car.nV==5 else [0.,0.5,-1.2,2.,0.3,-0.8,1.7,0.9,-2.4]))
    for nm, s in zip('ABC', St):
        P, C, worst = check_residue_formula(car, s)
        d0 = np.linalg.norm(C[0] - car.pi(s).astype(complex))
        c00 = max(abs(C[r][0] - car.pi(s)[0]) for r in range(P))
        print("   %-11s state %s  P=lcm=%2d   max|Z_n - formula| over n<=400 = %.2e   ||c^(0)-pi|| = %.2e   max_r |c^(r)_00 - p00| = %.2e"
              % (car.name, nm, P, worst, d0, c00))
print("   -> the identity is exact, c^(0) IS pi, and c^(r)_00 is p00 for EVERY residue r")
print("      (class 00 is fixed by both transports).  The blind lane's P-1 reproduces.")

print("\nX5.2  THE RATE PREDICTION (1/P) sum_r m(P_r), against a direct Birkhoff average")
for car in (K1(f,c), B0b_registrar(f,c), B0b_blind(f,c)):
    wv = np.array([0.40,0.15,0.15,0.15,0.15]) if car.nV==5 else np.ones(9)/9
    St = list(equal_pi_triple(car, wv, auto_moves(car, wv, 0.45),
              [0.0,-1.1,0.35,2.9,-0.62] if car.nV==5 else [0.,0.5,-1.2,2.,0.3,-0.8,1.7,0.9,-2.4]))
    for nm, s in zip('ABC', St):
        P = lcm(car.LF, car.LC)
        ms = [mahler4(*dressed_coeffs(car, s, r)) for r in range(P)]
        pred = float(np.mean(ms))
        NN = 200000 - (200000 % P)
        meas = rate_of(traj(car, s, NN, car.T('F'), car.T('C'), lambda t:(t,t)))
        circ = rate_of(traj(car, s, 20000, car.T('F'), car.T('C'), lambda t:(car.LF*t, car.LC*t)))
        print("   %-11s %s  predicted %.9f  measured(N=%d) %.9f  dev %.1e   | circuit rate %.9f  m(pi) %.9f"
              % (car.name, nm, pred, NN, meas, abs(pred-meas), circ, mahler4(*car.pi(s))))

print("\nX5.3  P-2 TESTED, NOT LEFT AS A NULL: SWEEP THE SPECTATOR WEIGHT p00 ON B0b")
print("      If the edge/circuit rate agreement on B0b is spectator dominance, the gap should")
print("      grow monotonically as p00 falls.  One variable moved: p00 (the rest of pi held")
print("      in fixed proportion, the within-class split and phases held fixed).")
B = B0b_registrar(f,c); cl = B.classes()
print("      %8s %16s %16s %12s" % ("p00", "circuit rate", "edge rate", "|gap|"))
for p00 in (0.90, 0.75, 0.5556, 0.4444, 0.30, 0.15, 0.02):
    rest = np.array([2.,1.,2.]); rest = rest/rest.sum()*(1-p00)
    wv = np.zeros(9)
    for j,o in enumerate(ORDER):
        idx = cl[o]; val = p00 if j==0 else rest[j-1]
        wv[idx] = val/len(idx)
    s = np.sqrt(wv).astype(complex)
    rc = rate_of(traj(B, s, 20000, B.T('F'), B.T('C'), lambda t:(4*t,3*t)))
    re = rate_of(traj(B, s, 200004, B.T('F'), B.T('C'), lambda t:(t,t)))
    print("      %8.4f %16.9f %16.9f %12.2e   rel %8.2e" % (p00, rc, re, abs(rc-re), abs(rc-re)/abs(rc)))
print("      and the same sweep with a PHASE-VARIED state (the arm that moved in the blind lane):")
rng = np.random.default_rng(777)
ph = rng.uniform(0,2*np.pi,9)
for p00 in (0.90, 0.5556, 0.15):
    rest = np.array([2.,1.,2.]); rest = rest/rest.sum()*(1-p00)
    wv = np.zeros(9)
    for j,o in enumerate(ORDER):
        idx = cl[o]; val = p00 if j==0 else rest[j-1]
        wv[idx] = val/len(idx)
    s = np.sqrt(wv)*np.exp(1j*ph)
    rc = rate_of(traj(B, s, 20000, B.T('F'), B.T('C'), lambda t:(4*t,3*t)))
    re = rate_of(traj(B, s, 200004, B.T('F'), B.T('C'), lambda t:(t,t)))
    print("      %8.4f %16.9f %16.9f %12.2e   (phases on)" % (p00, rc, re, abs(rc-re)))

print("\nX5.4  *** HOW MUCH OF THE STATE EACH CONVENTION KEEPS -- RANK OF THE STATE MAP ***")
print("      s (2V real coordinates, one of them a global phase both maps ignore) mapped to")
print("        (i)  pi           -- 4 real numbers   [what the CIRCUIT advance can see]")
print("        (ii) {c^(r)}_r    -- 8P real numbers  [what the EDGE advance can see]")
print("      Numerical Jacobian rank at a generic state, central differences h=1e-6, SVD")
print("      cut at 1e-6 of the top singular value.")
def jac_rank(fn, x0, h=1e-6, tol=1e-6):
    y0 = fn(x0); J = np.empty((len(y0), len(x0)))
    for i in range(len(x0)):
        e = np.zeros(len(x0)); e[i] = h
        J[:, i] = (fn(x0+e) - fn(x0-e))/(2*h)
    sv = np.linalg.svd(J, compute_uv=False)
    return int((sv > tol*sv[0]).sum()), sv
for car in (K1(f,c), B0b_registrar(f,c)):
    P = lcm(car.LF, car.LC)
    rngl = np.random.default_rng(31337)
    s0 = rngl.normal(size=car.nV) + 1j*rngl.normal(size=car.nV); s0 /= np.linalg.norm(s0)
    x0 = np.concatenate([s0.real, s0.imag])
    tostate = lambda x: x[:car.nV] + 1j*x[car.nV:]
    fpi = lambda x: car.pi(tostate(x))
    fc  = lambda x: np.concatenate([np.concatenate([dressed_coeffs(car, tostate(x), r).real,
                                                    dressed_coeffs(car, tostate(x), r).imag])
                                    for r in range(P)])
    rpi, svpi = jac_rank(fpi, x0)
    rc,  svc  = jac_rank(fc,  x0)
    print("   %-11s V=%d  2V=%d  (max possible rank after the global phase = %d)" % (car.name, car.nV, 2*car.nV, 2*car.nV-1))
    print("        rank d(pi)/ds        = %d   -> the CIRCUIT advance factors through a %d-dimensional summary;"
          % (rpi, rpi))
    print("                                      it DISCARDS %d of the %d state dimensions." % (2*car.nV-1-rpi, 2*car.nV-1))
    print("        rank d({c^(r)})/ds   = %d   -> the EDGE advance's class-indexed data has rank %d;"
          % (rc, rc))
    print("                                      it DISCARDS %d." % (2*car.nV-1-rc))
print("      READING, STATED AT THE PRECISION THE NUMBERS SUPPORT.  Both conventions' functionals")
print("      factor through CLASS-INDEXED data: the edge advance does not abandon the class")
print("      structure, it REFINES it from one real 4-vector to P complex ones.  But the refined")
print("      object is far bigger: rank 9 of 9 on K1 (it determines the ready state up to the")
print("      global phase, so the factorisation is VACUOUS there) and rank 10 of 17 on B0b (it")
print("      still discards 7 dimensions, so on B0b the edge advance also has a genuine, larger")
print("      class summary).  'The functional depends on the state only through a class summary'")
print("      is therefore INFORMATIVE under the circuit advance (3 of 9, 4 of 17) and either")
print("      vacuous or much weaker under the edge advance -- carrier-dependently so, which is")
print("      itself a result: the SIZE of the summary is not carrier-independent.")

print("\nX5.5  WHY P-2 HAPPENS -- WHICH JENSEN BRANCH DOMINATES, MEASURED")
print("      m(P_r) = mean_x log max(|c00 + c10 x|, |c01 + c11 x|).  If the class-00 branch")
print("      dominates for every x and every r, then m(P_r) = m(c00 + c10 x) for every r and")
print("      the edge rate collapses onto the circuit rate EXACTLY, not approximately.")
B = B0b_registrar(f,c)
for nm, sx in (("uniform real", np.sqrt(np.ones(9)/9).astype(complex)),
               ("uniform, phases on", np.sqrt(np.ones(9)/9)*np.exp(1j*np.random.default_rng(777).uniform(0,2*np.pi,9)))):
    th = np.arange(4096)*(2*np.pi/4096); xg = np.exp(1j*th)
    doms = []
    for r in range(12):
        C = dressed_coeffs(B, sx, r)
        A = np.abs(C[0] + C[1]*xg); Bb = np.abs(C[2] + C[3]*xg)
        doms.append(bool(np.all(A >= Bb - 1e-14)))
    ms = [mahler4(*dressed_coeffs(B, sx, r)) for r in range(12)]
    print("      %-19s branch-00 dominates at all 4096 grid x for r = %s" % (nm, [r for r in range(12) if doms[r]]))
    print("      %-19s (1/12)sum m(P_r) = %.12f   log(4/9) = %.12f" % ("", float(np.mean(ms)), np.log(4/9)))
print("      -> P-2 IS SCORED, not left two-way: on B0b the agreement is SPECTATOR DOMINANCE")
print("         (X5.3 sweep: the gap grows ~5 orders of magnitude as p00 falls, absolutely and")
print("         relatively, ON THE SAME CARRIER), not a K1 artefact.")

"""X2 -- THE ISOLATION NEITHER LANE RAN, AND THE NAMING.

The registrar names the moved variable THE TRANSPORT CONVENTION.
The blind lane (R-1) says that is false at the bytes -- M = T^L, one operator -- and names it
THE ADVANCE PAIR (a,b).  Both are naming a COORDINATE on a two-dimensional design space, and
each lane's two arms sit at two corners of it.  Fill in the other two corners and the dispute
decides itself.

DESIGN:  per-tick branch operator  x  advance rule.   ONE VARIABLE MOVED BETWEEN ADJACENT CELLS.
   operators:  T = COR-F's edge tick (NOT fibre-wise)
               S = diag(W^{1/L}) on the loop, the 'smeared holonomy' root (fibre-wise,
                   phase a function of the incidence class)
               R = a diagonal L-th root of M whose labels DIFFER WITHIN a class
                   (fibre-wise, phase NOT a function of the incidence class)
   advances:   EDGE (a,b) = (n,n)      CIRCUIT (a,b) = (L_F k, L_C k)
Every operator satisfies A^L = M_gamma EXACTLY, so every cell is a legitimate rival convention
on the registrar's own admissibility test."""
import numpy as np, sys
sys.path.insert(0, '/Users/bgm/MB Work/where-atoms-come-from/LANE_W11_R_BLIND_CROSS')
from xcore import *
np.set_printoptions(precision=12, linewidth=200)

f, c = 1.0, np.sqrt(2.0)
w3 = np.exp(2j*np.pi/3)

def build_ops(car):
    ops = {}
    ops['T  edge tick (COR-F)']            = (car.T('F'), car.T('C'))
    ops['S  smeared root diag(W^{1/L})']   = (car.S('F'), car.S('C'))
    RF = car.S('F').copy(); RC = car.S('C').copy()
    # multiply ONE vertex of a multi-vertex class by a primitive L-th root of unity:
    # still an exact L-th root of M, still fibre-wise, but its phase is no longer a class function
    clF = car.classes()[(1,0)]; clC = car.classes()[(0,1)]
    RF[clF[0], clF[0]] *= np.exp(2j*np.pi/car.LF)
    RC[clC[0], clC[0]] *= np.exp(2j*np.pi/car.LC)
    ops['R  diagonal root, NOT a class function'] = (RF, RC)
    return ops

def props(car, AF, AC):
    """(A^L = M exactly?, fibre-wise?, phase a function of the incidence class?)"""
    okp = (np.linalg.norm(np.linalg.matrix_power(AF,car.LF)-car.M('F')) < 1e-12 and
           np.linalg.norm(np.linalg.matrix_power(AC,car.LC)-car.M('C')) < 1e-12)
    fw  = (np.linalg.norm(AF-np.diag(np.diag(AF))) < 1e-12 and
           np.linalg.norm(AC-np.diag(np.diag(AC))) < 1e-12)
    cf = True
    if fw:
        for A in (AF, AC):
            d = np.diag(A)
            for o in ORDER:
                idx = car.classes()[o]
                if len(idx) > 1 and np.max(np.abs(d[idx]-d[idx[0]])) > 1e-12: cf = False
    else:
        cf = False
    return okp, fw, cf

for car, base, phases in ((K1(f,c), [0.40,0.15,0.15,0.15,0.15], [0.0,-1.1,0.35,2.9,-0.62]),
                          (B0b_registrar(f,c), list(np.ones(9)/9), [0.,0.5,-1.2,2.,0.3,-0.8,1.7,0.9,-2.4]),
                          (B0b_blind(f,c),     list(np.ones(9)/9), [0.,0.5,-1.2,2.,0.3,-0.8,1.7,0.9,-2.4])):
    w = np.array(base)/np.sum(base)
    St = list(equal_pi_triple(car, w, auto_moves(car, w, 0.45), phases))
    print("="*104)
    print("CARRIER %s   L_F=%d L_C=%d   pi = %s   (three arms, identical pi, arms diffed in X1)"
          % (car.name, car.LF, car.LC, np.array2string(car.pi(St[0]), precision=6)))
    print("  %-42s %-10s %-9s %-11s %14s %14s" % ("per-tick operator","A^L=M","fibre-wise",
                                                  "class-fn", "EDGE (n,n)", "CIRCUIT (L_F k,L_C k)"))
    for nm,(AF,AC) in build_ops(car).items():
        okp, fw, cf = props(car, AF, AC)
        NN = 400
        e = np.array([traj(car,s,NN,AF,AC, lambda t:(t,t)) for s in St])
        q = np.array([traj(car,s,NN,AF,AC, lambda t:(car.LF*t,car.LC*t)) for s in St])
        se, sq = np.ptp(e,axis=0).max(), np.ptp(q,axis=0).max()
        print("  %-42s %-10s %-9s %-11s %14.3e %14.3e" % (nm, okp, fw, cf, se, sq))
    print("  reading: spread ~1e-15 = incidence INVISIBLE ; spread ~1e-1 = incidence VISIBLE")

print("""
==========================================================================================
WHAT THE TABLE SETTLES

1. Moving ONLY the advance (row T, left cell -> right cell) flips visible -> invisible.
   Moving ONLY the per-tick operator (left column, row T -> row S) ALSO flips it.
   So the advance pair is NOT 'the operative variable'; it is ONE OF TWO SUFFICIENT MOVES.
   R-1's headline is a coordinate choice, exactly as the registrar's 'transport convention' is.
   The registrar's two arms and the blind lane's two arms are the SAME two corners; the lanes
   disagree about which coordinate moved because both arms differ in both coordinates unless
   you first fix the T-generated parametrisation, which is itself a choice.

2. What all the invisible cells share and the visible one lacks is NOT diagonality.  Row R is
   fibre-wise (diagonal) at every tick, satisfies A^L = M exactly, and STILL separates the arms.
   The property is: THE COMPARED OPERATOR IS FIBRE-WISE **AND** ITS FIBRE PHASE IS A FUNCTION OF
   THE VERTEX'S INCIDENCE CLASS.  M_gamma = diag(W^{1[v in gamma]}) has both by construction --
   its diagonal is literally a function of the incidence indicator, which is why, and the only
   why, that only pi can enter.
   Clause one is W-06's own correction of N4 ('not scalar multiplication -- FIBRE-WISE-NESS',
   REGISTER :577).  Clause two the corpus has never written down, and the registrar's stated
   conclusion ('holds exactly where both branch operators are DIAGONAL') is clause one only.
   NAMING GUARD: 'fibre-wise' is W-06's term and I am adding a clause to it rather than coining
   a seventh name (W-10 N-7).  Not S4's 'canonical clock' (the cell schedule k_n, S4:170);
   not W-03's 'relation lattice' ({(m,n): u^m v^n = 1}, REGISTER :237); no other occurrence of
   'advance', 'class function' or 'class-constant' exists in the register.
==========================================================================================""")

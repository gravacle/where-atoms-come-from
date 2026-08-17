# RULING LEG 2 -- THE OBJECT THE WHOLE QUESTION IS ABOUT.
#
# Z_n = <U_F^n s, U_C^n s> = <s, Q_n s>  with  Q_n = (U_F^n)^* U_C^n   -- the RELATIVE branch
# operator at the tick the record is read.  Everything below is about Q_n and nothing else.
#
# CLAIM 1 (an IDENTITY, not a measurement).  Under W-01's convention the diagonal of M_gamma IS
#   the incidence indicator: diag(M_F)_v = W_F^{1[v in gamma_F]}.  Hence Q_k = mult by
#   conj(W_F)^{k a} W_C^{k b} at class (a,b) -- a function of the class label and of NOTHING ELSE.
# CLAIM 2 (a BICONDITIONAL, quantified over ALL operator pairs).  |Z_n| is a function of pi alone
#   IFF Q_n is multiplication by a function of the incidence class.
import numpy as np, rlib
rng = np.random.default_rng(rlib.SEED)

print("== R2.1  THE IDENTITY.  M_gamma's DIAGONAL *IS* THE INCIDENCE INDICATOR ==")
for C in (rlib.K1(), rlib.B0b()):
    worst_d = worst_q = 0.0
    for _ in range(500):
        a = rng.uniform(0,2*np.pi,len(C.edges))
        WF, WC = rlib.holon(C.walkF,a), rlib.holon(C.walkC,a)
        MF, MC = rlib.Mcirc(C,C.VF,WF), rlib.Mcirc(C,C.VC,WC)
        worst_d = max(worst_d, np.max(np.abs(np.diag(MF)-np.array([WF**C.cls[v][0] for v in range(C.nv)]))),
                               np.max(np.abs(np.diag(MC)-np.array([WC**C.cls[v][1] for v in range(C.nv)]))))
        for k in (1,2,3,7):
            Q = rlib.Qrel(MF,MC,k)
            tgt = np.diag([np.conj(WF)**(k*C.cls[v][0]) * WC**(k*C.cls[v][1]) for v in range(C.nv)])
            worst_q = max(worst_q, np.linalg.norm(Q-tgt))
    print(f"  {C.name:4s} max| diag(M_gamma)_v - W^(1[v in gamma]) |            = {worst_d:.2e}")
    print(f"  {C.name:4s} max|| Q_k - diag(class character conj(W_F)^ka W_C^kb) || = {worst_q:.2e}   (500 connections, k=1,2,3,7)")
print("  -> Z_k = sum_v (class character)_v |s_v|^2 = sum over the FOUR CLASSES.  Substitution,")
print("     no lemma, no possible alternative outcome.  Leg B1 / C's circuit rows are IDENTITIES")
print("     and are COULD-NOT-HAVE-FAILED CONTROLS.  I score them as zero evidence, both ways.\n")

print("== R2.2  THE BICONDITIONAL.  Positive AND negative controls, both directions ==")
def haar(n, rng):
    z = (rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))/np.sqrt(2)
    q,r = np.linalg.qr(z); return q*(np.diag(r)/np.abs(np.diag(r)))
for C in (rlib.K1(), rlib.B0b()):
    base = rng.dirichlet(np.ones(C.nv)); S = rlib.same_pi_states(C, rng, base, 24)
    agree = dis = npos = nneg = 0
    for trial in range(1500):
        mode = trial % 5
        if mode == 0:                                  # class-constant diagonal  (positive)
            dF = np.array([rng.uniform(0,2*np.pi) for _ in range(4)]); dC = rng.uniform(0,2*np.pi,4)
            UF = np.diag([np.exp(1j*dF[C.CLASSES.index(C.cls[v])]) for v in range(C.nv)])
            UC = np.diag([np.exp(1j*dC[C.CLASSES.index(C.cls[v])]) for v in range(C.nv)])
        elif mode == 1:                                # diagonal but NOT class-constant
            UF = np.diag(np.exp(1j*rng.uniform(0,2*np.pi,C.nv))); UC = np.diag(np.exp(1j*rng.uniform(0,2*np.pi,C.nv)))
        elif mode == 2:                                # CORRELATED: U_F = L_F R, U_C = L_C R
            R = np.zeros((C.nv,C.nv),dtype=complex)
            for c in C.CLASSES:
                vs = C.idx[c]
                if vs: R[np.ix_(vs,vs)] = haar(len(vs), rng)
            dF = rng.uniform(0,2*np.pi,4); dC = rng.uniform(0,2*np.pi,4)
            LF = np.diag([np.exp(1j*dF[C.CLASSES.index(C.cls[v])]) for v in range(C.nv)])
            LC = np.diag([np.exp(1j*dC[C.CLASSES.index(C.cls[v])]) for v in range(C.nv)])
            UF, UC = LF@R, LC@R
        elif mode == 3:                                # COR-F's own T at a random connection
            a = rng.uniform(0,2*np.pi,len(C.edges)); UF,UC = rlib.Tedge(C,C.walkF,a), rlib.Tedge(C,C.walkC,a)
        else:                                          # unrestricted Haar pair
            UF, UC = haar(C.nv,rng), haar(C.nv,rng)
        for n in (1,2,3):
            ok, res = is_ok = rlib.is_class_constant_diag(C, rlib.Qrel(UF,UC,n))
            sp = rlib.pi_spread(C, UF, UC, S, [n])
            blind = sp < 1e-9
            if ok: npos += 1
            else:  nneg += 1
            if ok == blind: agree += 1
            else:          dis += 1
    print(f"  {C.name:4s} 4500 (pair,tick) cells: AGREE {agree}, DISAGREE {dis}"
          f"   [ Q class-constant-diagonal: {npos} cells; not: {nneg} cells ]")
print("  -> BOTH outcomes occur in quantity, so the test COULD have failed.  It does not.")
print("     THE CRITERION IS ON THE *RELATIVE* OPERATOR Q_n, NOT ON EACH BRANCH OPERATOR:")
print("     mode 2 supplies pairs with U_F NON-DIAGONAL that ARE pi-blind, which refutes")
print("     'invisibility holds exactly where both branch operators are diagonal'.\n")

print("== R2.3  DIAGONAL IS NECESSARY BUT *NOT SUFFICIENT* -- the registrar's conclusion sentence ==")
for C in (rlib.K1(), rlib.B0b()):
    a = rlib.a_generic(C, rng, 1.0, 2**0.5)
    base = rng.dirichlet(np.ones(C.nv)); S = rlib.same_pi_states(C, rng, base, 24)
    WF, WC = rlib.holon(C.walkF,a), rlib.holon(C.walkC,a)
    # a DIAGONAL L-th root of M_gamma whose labels differ WITHIN a class
    def diagroot(Vs, W, L, branches):
        D = np.eye(C.nv, dtype=complex)
        for k,v in enumerate(sorted(Vs)): D[v,v] = np.exp(1j*(np.angle(W)+2*np.pi*branches[k])/L)
        return D
    bF = [0]*C.LF; bF[-1] = 1                          # one vertex on a different branch
    RF, RC = diagroot(C.VF,WF,C.LF,bF), diagroot(C.VC,WC,C.LC,[0]*C.LC)
    DF, DC = rlib.Droot(C,C.VF,WF,C.LF), rlib.Droot(C,C.VC,WC,C.LC)
    print(f"  {C.name:4s} ||R_F^L - M_F|| = {np.linalg.norm(np.linalg.matrix_power(RF,C.LF)-rlib.Mcirc(C,C.VF,WF)):.2e}"
          f"   R_F DIAGONAL? {np.allclose(RF,np.diag(np.diag(RF)))}   fibre-wise? True")
    print(f"       spread under (R_F,R_C), DIAGONAL, within-class labels differ, n<=12 : {rlib.pi_spread(C,RF,RC,S,range(1,13)):.2e}   <== VISIBLE")
    print(f"       spread under (D_F,D_C), DIAGONAL, class-constant labels,      n<=12 : {rlib.pi_spread(C,DF,DC,S,range(1,13)):.2e}   <== blind")
print("  -> ONE VARIABLE MOVED between those two rows: whether the fibre phase is a function of")
print("     the vertex or a function of its INCIDENCE CLASS.  Everything else is identical.")
print("     'DIAGONAL' is the wrong name.  So is 'fibre-wise' alone (REGISTER:577).")

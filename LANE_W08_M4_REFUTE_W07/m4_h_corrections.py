# W-08 / M4 leg H — WHICH PUBLISHED NUMBERS OF W-07 ARE WRONG, AND WHAT THEY SHOULD READ.
# Recompute W-07 sec3's isolation table with the ratio that actually governs its own observable
# (rho_true = W_F^dF W_C^-dC, m4_b), on W-07's own connections, state, observable and k-range.
#
# ISOLATION LEDGER.  Held fixed: everything W-07 held fixed, plus the connections themselves
# (reconstructed from W-07's own leg E source lines).  Moved: the ratio formula ALONE.
import numpy as np
FACE_V={0,1,2}; CYC_V={0,3,4}; TREE={1:(0,),2:(0,1),3:(3,),4:(3,4)}
def dress(s,a):
    u=np.exp(1j*np.asarray(a)); t=np.array(s,dtype=complex)
    for v,p in TREE.items():
        w=1.0+0j
        for e in p: w*=u[e]
        t[v]=s[v]/w
    return t
rng=np.random.default_rng(20260816)
s=rng.normal(size=5)+1j*rng.normal(size=5); s/=np.linalg.norm(s)     # W-07 leg E's state
K=4000; k=np.arange(1,K+1)
cases=[("S1 PUBLISHED (pi/3 x3, pi/2 x3)", np.array([np.pi/3]*3+[np.pi/2]*3)),
       ("GENERIC sqrt2/sqrt3            ", np.array([2*np.pi*(2**0.5)/3]*3+[2*np.pi*(3**0.5)/3]*3)),
       ("RANDOM seed 1                  ", np.random.default_rng(1).uniform(0,2*np.pi,6)),
       ("RANDOM seed 2                  ", np.random.default_rng(2).uniform(0,2*np.pi,6)),
       ("RANDOM seed 3                  ", np.random.default_rng(3).uniform(0,2*np.pi,6))]
u,v=2,3; dF=(v in FACE_V)-(u in FACE_V); dC=(v in CYC_V)-(u in CYC_V)
print("  W-07 sec3 table, as published (from leg E, rho_E = conj(W_F)^dF W_C^-dC)")
print("  vs the same table computed with rho_true = W_F^dF W_C^-dC on the same inputs.\n")
print(f"  {'connection':<33} {'arg(rho_E)/2pi':>15} {'arg(rho_true)/2pi':>18} {'W-07 min':>11} {'TRUE min':>11} {'cells<1e-9':>11}")
pub=[("-0.250000000","6.729e-19","1000"),("-0.317837245","1.567e-04","0"),("-0.077362414","1.722e-05","0"),
     ("-0.046247977","3.095e-04","0"),("+0.014316520","3.057e-04","0")]
for (tag,a),(argp,minp,cellp) in zip(cases,pub):
    WF=np.exp(1j*sum(a[:3])); WC=np.exp(1j*sum(a[3:]))
    rho_t=WF**dF*WC**(-dC); rho_e=np.conj(WF)**dF*WC**(-dC)
    t=dress(s,a); amp=abs(np.conj(t[u])*t[v])
    D=amp*np.abs(WF**(dF*k)-WC**(dC*k))
    print(f"  {tag:<33} {np.angle(rho_e)/(2*np.pi):>15.9f} {np.angle(rho_t)/(2*np.pi):>18.9f} {minp:>11} {D.min():>11.3e} {int((D<1e-9).sum()):>11}")
print()
print("  The 'arg(rho)/2pi' column of W-07 sec3 is the argument of a DIFFERENT group element on")
print("  four of five rows.  The 'min over k<=4000' column is a min of a different sequence on")
print("  those four rows.  The 'cells < 1e-9' column is unaffected (both elements are irrational")
print("  there), so the VERDICT of the table survives its own arithmetic — but the table's own")
print("  numbers are not reproducible from the quantity the page says it is tabulating.")
print("  This is COR-K's defect class, of record against S3 sec6(f), recurring at W-07 sec3.\n")
print("  And row 1's 'min = 6.729e-19': the exact value is 0 (m4_g G3a).  W-07's own leg C2")
print("  prints 0 for the same quantity.  The page publishes the float residue in sec3 and the")
print("  exact 0 in the sec3 scaling table, three lines apart, without noting they are the same k.")

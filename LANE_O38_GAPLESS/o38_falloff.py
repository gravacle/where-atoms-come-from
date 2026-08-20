"""O-38: DOES THE EXTENSIVE QUANTITY HAVE A RANGE, AND IS ITS FALLOFF A POWER LAW?

C-47 established the first record-level quantity meeting gravity's source standard: the bath-induced
effective potential Phi(s) = -(1/beta) ln Tr exp(-beta H_B(s)) is strictly extensive at a FIXED
environment, exactly additive over disjoint regions, and not a count.

But it is CONTACT-OR-NOTHING on every venue tested -- exactly zero between records that do not share
a bath site -- and where a mediator was allowed to propagate, the only falloff ever measured was
EXPONENTIAL, xi = 1.1250 sites, on a GAPPED bath. A screened interaction with a fixed screening
length is ruled out at any N by its form, however large N gets. GRAVITY NEEDS A POWER LAW.

THE POINT: this is a BATH-HAMILTONIAN CHOICE, not a record-count problem. A gapped mediator has
exponentially decaying correlations by construction; a GAPLESS one has power-law correlations by
construction. Nobody has run it.

  H_B = sum_j e_j Z_j + g sum_j X_j X_{j+1}      transverse-field Ising chain
  critical at |e| = |g| -- GAPLESS, correlations power-law
  gapped at |e| >> |g|  -- the D-15 CONTROL, which must reproduce an exponential falloff

Two lumps of records at bath-site separation r. The measured quantity is the CROSS-LUMP part of the
induced potential, U(r) = spread(Phi_AB) - spread(Phi_A) - spread(Phi_B). The variable is r, named
first per the relevance test.

DEAD-CONFIGURATION CONTROL: with one record per site Phi is even in each s_i and the whole effect is
identically zero. Records must CROWD onto a lump's sites for U to be non-zero at all. That control is
run and reported, because a null without it means nothing."""
import sys, itertools, numpy as np
def say(*a): print(*a); sys.stdout.flush()
I2=np.eye(2); X=np.array([[0,1],[1,0]],dtype=float); Z=np.array([[1,0],[0,-1]],dtype=float)
def op(nq,j,P):
    M=np.array([[1.0]])
    for i in range(nq): M=np.kron(M,P if i==j else I2)
    return M
def bath_H(nq,e,g):
    """transverse-field Ising chain: sum_j e Z_j + g sum_j X_j X_{j+1} (open chain)"""
    H=sum(e*op(nq,j,Z) for j in range(nq))
    H=H+sum(g*(op(nq,j,X)@op(nq,j+1,X)) for j in range(nq-1))
    return H
def phi(nq,e,g,lam,beta,coupling,s):
    """free energy of the bath given record signs s. coupling: list of (record_index, site)."""
    H=bath_H(nq,e,g).copy()
    for (ri,site) in coupling: H=H+lam*s[ri]*op(nq,site,Z)
    w=np.linalg.eigvalsh(H)
    m=w.min()
    return -(1.0/beta)*(np.log(np.exp(-beta*(w-m)).sum())-beta*m)
def spread(nq,e,g,lam,beta,coupling,nrec):
    vals=[phi(nq,e,g,lam,beta,coupling,s) for s in itertools.product((1,-1),repeat=nrec)]
    return max(vals)-min(vals)
say("="*104); say("O-38   IS THE FALLOFF A POWER LAW ON A GAPLESS MEDIATOR?"); say("="*104)
NQ=10; BETA=2.0; LAM=0.8; NA=NB=2      # 2 records per lump, both CROWDED onto one site each
say(f"  bath: open transverse-field Ising chain, {NQ} sites, beta={BETA}, lam={LAM}")
say(f"  two lumps of {NA} records each, CROWDED (both records of a lump on the SAME site)")
say(f"  U(r) = spread(Phi_AB) - spread(Phi_A) - spread(Phi_B)   -- the cross-lump induced energy")
say("")
for tag,e,g in (("GAPPED   e=1.0 g=0.1  (CONTROL: must be exponential)",1.0,0.1),
                ("CRITICAL e=1.0 g=1.0  (GAPLESS: power-law correlations)",1.0,1.0),
                ("ORDERED  e=0.1 g=1.0  (CONTROL: gapped the other side)",0.1,1.0)):
    say("-"*104); say(f"  {tag}"); say("-"*104)
    say(f"    {'r':>4}{'siteA':>7}{'siteB':>7}{'U(r)':>16}{'|U| vs U(1)':>15}")
    rows=[]
    for r in range(1,NQ-1):
        sa=0; sb=r
        if sb>=NQ: break
        cAB=[(0,sa),(1,sa),(2,sb),(3,sb)]
        cA =[(0,sa),(1,sa)]
        cB =[(0,sb),(1,sb)]
        U=spread(NQ,e,g,LAM,BETA,cAB,4)-spread(NQ,e,g,LAM,BETA,cA,2)-spread(NQ,e,g,LAM,BETA,cB,2)
        rows.append((r,U))
        say(f"    {r:>4}{sa:>7}{sb:>7}{U:>16.9f}{(abs(U)/abs(rows[0][1]) if abs(rows[0][1])>1e-14 else float('nan')):>15.6f}")
    good=[(r,abs(U)) for r,U in rows if abs(U)>1e-11]
    if len(good)<4:
        say(f"    -> only {len(good)} points above the 1e-11 floor; NO falloff can be fitted"); continue
    r=np.array([x[0] for x in good],dtype=float); u=np.array([x[1] for x in good])
    # power law: log|U| vs log r     exponential: log|U| vs r
    pw=np.polyfit(np.log(r),np.log(u),1); rp=np.log(u)-np.polyval(pw,np.log(r))
    ex=np.polyfit(r,np.log(u),1);        re=np.log(u)-np.polyval(ex,r)
    say(f"    power law   |U| ~ r^{pw[0]:+.4f}      rms residual {np.sqrt((rp**2).mean()):.5f}")
    say(f"    exponential |U| ~ exp({ex[0]:+.4f} r)  xi = {(-1/ex[0] if ex[0]<0 else float('nan')):.4f} sites"
        f"   rms residual {np.sqrt((re**2).mean()):.5f}")
    say(f"    -> {'POWER LAW fits better' if np.sqrt((rp**2).mean())<np.sqrt((re**2).mean()) else 'EXPONENTIAL fits better -- screened, ruled out at any N by form'}")
say("")
say("-"*104); say("  DEAD-CONFIGURATION CONTROL: one record per site, NOT crowded"); say("-"*104)
say("  Phi is then even in each s_i and the cross term must vanish identically.")
say(f"    {'r':>4}{'U(r) uncrowded':>18}")
for r in (1,2,3,4):
    cAB=[(0,0),(1,1),(2,r+1),(3,r+2)]
    cA=[(0,0),(1,1)]; cB=[(0,r+1),(1,r+2)]
    if r+2>=NQ: break
    U=spread(NQ,1.0,1.0,LAM,BETA,cAB,4)-spread(NQ,1.0,1.0,LAM,BETA,cA,2)-spread(NQ,1.0,1.0,LAM,BETA,cB,2)
    say(f"    {r:>4}{U:>18.9f}")
say("    -> if these are NOT ~0, the crowding hypothesis behind C-47 is wrong.")
say("")
say("="*104); say("  READ -- from the numbers above"); say("="*104)

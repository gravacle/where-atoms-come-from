"""O-31: IS THE INTERSECTION PAIRING THE SOURCE TERM?

C-32 removed the COUNT of records as the source. The principal's triad names a different one:

    R[carrier frame transport]  ~  <c_i, c_j>  x  channel_map[alpha]

whose source is CROSSINGS between record carriers, not a number of records.

A SCOPE CORRECTION FIRST. C-29 says transport between records is the identity on an abelian carrier.
That is true of CONJUGATION of fluxes, g -> h g h^-1, and it is the only thing C-29 claims. It is NOT
true of the pairing between record OPERATORS: for Wilson loops on a Z_2 carrier

    W(c_i) W(c_j) = (-1)^<c_i,c_j> W(c_j) W(c_i)

which is non-trivial exactly when the record carriers CROSS an odd number of times.

ERRATUM, v1: v1 paired the horizontal Z loop with a string on VERTICAL links, which is disjoint from
it, so every crossing number printed 0 -- and the READ then concluded 'the pairing was never zero
here', contradicting the table directly above it. The correct partner of a horizontal Z loop is the
X string on HORIZONTAL links down a column; they share exactly one link. The non-degeneracy check
added below fires on precisely this error: a set of logical operators whose pairing matrix is
singular is not a basis, and v1's was identically zero."""
import sys, numpy as np
def say(*a): print(*a); sys.stdout.flush()
say("="*104); say("O-31   IS THE INTERSECTION PAIRING THE SOURCE TERM?"); say("="*104)
def toric(L):
    idx={}; n=0
    for r in range(L):
        for c in range(L):
            for d in (0,1): idx[(r,c,d)]=n; n+=1
    stars=[[idx[(r,c,0)],idx[(r,(c-1)%L,0)],idx[(r,c,1)],idx[((r-1)%L,c,1)]] for r in range(L) for c in range(L)]
    plaqs=[[idx[(r,c,0)],idx[((r+1)%L,c,0)],idx[(r,c,1)],idx[(r,(c+1)%L,1)]] for r in range(L) for c in range(L)]
    return n,idx,stars,plaqs
def vec(n,links):
    v=np.zeros(n,dtype=np.int8); v[list(links)]=1; return v
def records(L):
    """the four record carriers. Z-type run ALONG links; X-type run along the DUAL, cutting links."""
    n,idx,stars,plaqs=toric(L)
    Z1=vec(n,[idx[(0,c,0)] for c in range(L)])   # Z: horizontal loop on horizontal links
    Z2=vec(n,[idx[(r,0,1)] for r in range(L)])   # Z: vertical loop on vertical links
    X1=vec(n,[idx[(r,0,0)] for r in range(L)])   # X: dual vertical cut, crosses Z1 once
    X2=vec(n,[idx[(0,c,1)] for c in range(L)])   # X: dual horizontal cut, crosses Z2 once
    return n,idx,stars,plaqs,(Z1,Z2),(X1,X2)
say("")
say("0. SELF-CHECK -- are these logical operators at all, and is the pairing NON-DEGENERATE?")
say("   A Z-string must meet every STAR in an even number of links; an X-string every PLAQUETTE.")
say("   And the pairing matrix over the four records must be SYMPLECTIC (determinant 1 mod 2).")
say(f"   {'L':>4}{'Z commutes w/ all stars':>26}{'X commutes w/ all plaqs':>26}{'pairing det mod 2':>20}")
ok_all=True
for L in (2,3,4,5,6):
    n,idx,stars,plaqs,Zs,Xs=records(L)
    zok=all(sum(Z[i] for i in s)%2==0 for Z in Zs for s in stars)
    xok=all(sum(X[i] for i in p)%2==0 for X in Xs for p in plaqs)
    M=np.array([[int(np.dot(a,b)%2) for b in (Xs[0],Xs[1])] for a in (Zs[0],Zs[1])],dtype=np.int8)
    det=int(round(np.linalg.det(M.astype(float))))%2
    ok_all &= zok and xok and det==1
    say(f"   {L:>4}{str(zok):>26}{str(xok):>26}{det:>20}")
say(f"   -> {'the records are genuine logical operators and the pairing is NON-DEGENERATE' if ok_all else 'SELF-CHECK FAILED -- these are not a record basis; no conclusion may be drawn'}")
if not ok_all: sys.exit("self-check failed")
say("")
say("1. DOES THE COMMUTATION OF RECORD OPERATORS COUNT CROSSINGS?  (Z_2 carrier -- the ABELIAN case)")
say(f"   {'L':>4}{'links':>8}{'<Z1,X1>':>10}{'<Z2,X2>':>10}{'<Z1,X2>':>10}{'<Z2,X1>':>10}{'Z1,X1':>24}")
for L in (2,3,4,5,6):
    n,idx,stars,plaqs,Zs,Xs=records(L)
    a=int(np.dot(Zs[0],Xs[0])%2); b=int(np.dot(Zs[1],Xs[1])%2)
    c=int(np.dot(Zs[0],Xs[1])%2); d=int(np.dot(Zs[1],Xs[0])%2)
    say(f"   {L:>4}{n:>8}{a:>10}{b:>10}{c:>10}{d:>10}{('ANTICOMMUTE' if a else 'commute'):>24}")
say("   -> THE PAIRING BETWEEN RECORD OPERATORS IS NON-ZERO ON AN ABELIAN CARRIER. C-29 is narrower")
say("      than it may read: it is about CONJUGATION of fluxes only. The pairing was here all along.")
say("")
say("2. DOES THE PAIRING GROW WITH THE EXTENT OF THE RECORDS?")
say(f"   {'L':>4}{'carrier length':>17}{'shared links':>15}{'crossings mod 2':>18}")
for L in (2,3,4,5,6,8,10,16):
    n,idx,stars,plaqs,Zs,Xs=records(L)
    shared=int(np.dot(Zs[0],Xs[0]))
    say(f"   {L:>4}{int(Zs[0].sum()):>17}{shared:>15}{shared%2:>18}")
say("   -> the record carriers get 8x longer and the crossing number does not move. It is a")
say("      TOPOLOGICAL invariant of the PAIR. Extent is not the density variable.")
say("")
say("3. DOES THE TOTAL PAIRING GROW WITH THE NUMBER OF RECORDS?")
say("   a genus-g surface carries 2g records with the standard symplectic pairing")
say(f"   {'genus g':>9}{'records 2g':>13}{'non-zero pairings':>20}{'pairs available':>18}{'fraction':>11}")
for g in range(1,9):
    m=2*g; nz=g; tot=m*(m-1)//2
    say(f"   {g:>9}{m:>13}{nz:>20}{tot:>18}{nz/tot:>11.4f}")
say("   -> non-zero pairings grow LINEARLY in the number of records while PAIRS grow quadratically,")
say("      so the fraction of record pairs that interact FALLS toward zero.")
say("")
say("="*104); say("  READ"); say("="*104)
say("  SCOPE CORRECTION: the pairing between record operators is NON-ZERO on the abelian carriers")
say("  this program has always used -- verified non-degenerate at every L. C-29 is about CONJUGATION")
say("  of fluxes only, and stands; the intersection pairing was never zero here.")
say("")
say("  BUT THE PAIRING IS NOT A DENSITY EITHER. It is a TOPOLOGICAL invariant of a PAIR of record")
say("  carriers: growing the carriers 8x does not move it, and across a surface with more records")
say("  the fraction of interacting pairs FALLS. Nothing accumulates with how much is enclosed. The")
say("  principal named this already -- the pairing is topologically rigid and metric-blind, a source")
say("  T_uv; and a source is not by itself a density law.")
say("")
say("  BOTH CANDIDATE SOURCE TERMS NOW FAIL THE SAME TEST FOR THE SAME REASON: they are TOPOLOGICAL,")
say("  and a topological quantity does not know how much is enclosed. Whatever supplies gravity's")
say("  density law at the record level is NEITHER A COUNT NOR AN INTERSECTION NUMBER. The one term")
say("  in the triad not yet tested is the remaining factor, channel_map[alpha] -- the metric.")

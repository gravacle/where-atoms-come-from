"""W-39b.  DO RECORDS CROWD EACH OTHER OUT, OR DID THE BATH JUST GET WEAKER?

W-39: protection margin falls 43.8x -> 2.0x -> 1.1x as 1, 2, 3 records are packed onto the carrier.
CONFOUND: the protecting bath necessarily SHRINKS as records are added (8 -> 5 -> 2 links), because
it must avoid every boundary. A weaker environment decoheres less, so the margin could fall for that
reason alone and have nothing to do with crowding.

SEPARATE THEM. Hold the TOTAL dissipation fixed by scaling gamma with 1/|bath|, so every packing is
compared at the same environmental strength. If the margin still collapses, records genuinely crowd
each other. If it recovers, W-39's reading is an artifact of bath size and must be withdrawn.

SECOND CONTROL: give each k the SAME SIZED bath (2 links) whether or not it needs to be small, so
bath size is literally constant across k.
"""
import itertools, numpy as np
exec(open('w39_capacity.py').read().split('print("W-39  HOW MANY')[0])

def rates(links,g2=0.01,gam=0.5):
    H=-MAG-g2*ELEC
    M=-1j*(np.kron(H,Id)-np.kron(Id,H.T))
    for X in [Zop(st,[k],N) for k in links]: M+=gam*(np.kron(X,X.conj())-np.kron(Id,Id))
    w,U=np.linalg.eig(M.conj().T)
    rate=-np.conj(w).real; U=U/np.linalg.norm(U,axis=0)
    out={}
    for S,(sup,O) in REG.items():
        v=(O/np.linalg.norm(O)).reshape(-1)
        ov=np.abs(U.conj().T@v); ov=ov/max(ov.sum(),1e-30)
        out[S]=float((ov*rate).sum())
    return out

def packings():
    got={}
    for k in range(1,5):
        for combo in itertools.combinations(list(REG),k):
            vecs=[sum(1<<i for i in S) for S in combo]; basis=[]; rank=0
            for v in vecs:
                for b in basis: v=min(v,v^b)
                if v: basis.append(v); basis.sort(reverse=True); rank+=1
            if rank<k: continue
            union=set().union(*[REG[S][0] for S in combo])
            free=[x for x in range(L) if x not in union]
            if free and k not in got: got[k]=(combo,free)
    return got

PK=packings()
print("W-39b  DOES THE MARGIN COLLAPSE SURVIVE A CONSTANT-STRENGTH ENVIRONMENT?")
print()
print("  RUN 1 -- gamma fixed at 0.5 (this is W-39, repeated). Bath size varies with k.")
print(f"  {'k':>2s} {'|bath|':>7s} {'gamma':>7s} {'gamma*|bath|':>13s} {'worst member':>13s} {'best non-member':>16s} {'margin':>8s}")
print("  "+"-"*76)
for k in sorted(PK):
    combo,free=PK[k]
    r=rates(free,gam=0.5)
    wm=max(r[S] for S in combo); bn=min(v for S,v in r.items() if S not in combo)
    print(f"  {k:2d} {len(free):7d} {0.5:7.3f} {0.5*len(free):13.2f} {wm:13.4e} {bn:16.4e} {bn/max(wm,1e-30):7.1f}x")

print()
print("  RUN 2 -- TOTAL dissipation held FIXED: gamma = 4.0/|bath|, so gamma*|bath| = 4.0 always.")
print(f"  {'k':>2s} {'|bath|':>7s} {'gamma':>7s} {'gamma*|bath|':>13s} {'worst member':>13s} {'best non-member':>16s} {'margin':>8s}")
print("  "+"-"*76)
for k in sorted(PK):
    combo,free=PK[k]
    g=4.0/len(free)
    r=rates(free,gam=g)
    wm=max(r[S] for S in combo); bn=min(v for S,v in r.items() if S not in combo)
    print(f"  {k:2d} {len(free):7d} {g:7.3f} {g*len(free):13.2f} {wm:13.4e} {bn:16.4e} {bn/max(wm,1e-30):7.1f}x")

print()
print("  RUN 3 -- bath size held LITERALLY CONSTANT at 2 links for every k (a subset of the free set).")
print(f"  {'k':>2s} {'bath':>10s} {'worst member':>13s} {'best non-member':>16s} {'margin':>8s}")
print("  "+"-"*56)
for k in sorted(PK):
    combo,free=PK[k]
    if len(free)<2: print(f"  {k:2d}   (fewer than 2 free links)"); continue
    b=free[:2]
    r=rates(b,gam=0.5)
    wm=max(r[S] for S in combo); bn=min(v for S,v in r.items() if S not in combo)
    print(f"  {k:2d} {str(b):>10s} {wm:13.4e} {bn:16.4e} {bn/max(wm,1e-30):7.1f}x")

print()
print("  READING: if RUN 2 and RUN 3 still show the margin falling with k, the collapse is CROWDING.")
print("  If the margin is restored once strength or size is held fixed, W-39's reading was the bath.")

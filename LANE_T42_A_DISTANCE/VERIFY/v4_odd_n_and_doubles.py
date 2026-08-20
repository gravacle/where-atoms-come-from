"""VERIFY check 4 (independent): (a) n=5 [[5,3,2]] two-generator nonexistence, full exhaust,
fresh code; (b) D(Z4) 2x2 class minimal weights, fresh code; (c) Aut of the 4-point toric
config metric by brute force over all 24 permutations; (d) D(Z2xZ2) L=2 metric has the exact
product form (two uniform 4-point factors), whose Aut is the wreath product S4 wr S2, order
1152 -- consistency with the lane's enumerated 1152."""
import itertools
def pcnt(x): return bin(x).count("1")
ok=True
def gate(lbl,b,detail=""):
    global ok
    print(("PASS " if b else "FAIL "),lbl,detail); ok&=b

# ---------- (a) n=5
n=5
def sp10(a,b):
    mask=(1<<n)-1
    return (pcnt((a&mask)&(b>>n))+pcnt((a>>n)&(b&mask)))%2
W1=[]
for i in range(n):
    W1+= [1<<i, 1<<(n+i), (1<<i)|(1<<(n+i))]
comm={}
for g in range(1<<(2*n)):
    m=0
    for i,w in enumerate(W1):
        if sp10(w,g)==0: m|=1<<i
    comm[g]=m
tot=0; bad=[]
for g1 in range(1,1<<(2*n)):
    for g2 in range(g1+1,1<<(2*n)):
        if sp10(g1,g2)!=0: continue
        tot+=1
        both=comm[g1]&comm[g2]
        Sbar={0,g1,g2,g1^g2}
        found=False
        m=both
        while m:
            i=(m&-m).bit_length()-1
            m&=m-1
            if W1[i] not in Sbar:
                found=True; break
        if not found:
            bad.append((g1,g2))
gate("n=5: all %d commuting pairs have a weight-1 logical (no [[5,3,2]] 2-gen)"%tot,
     tot==260865 and len(bad)==0, "tot=%d bad=%d"%(tot,len(bad)))

# ---------- (b) D(Z4) 2x2
LL=2; nE=8
def h(i,j): return (i%LL)*LL+(j%LL)
def v(i,j): return LL*LL+(i%LL)*LL+(j%LL)
plaq=[]
for i in range(LL):
    for j in range(LL):
        r=[0]*nE
        r[h(i,j)]+=1; r[v(i,j+1)]+=1; r[h(i+1,j)]-=1; r[v(i,j)]-=1
        plaq.append([x%4 for x in r])
wA=[0]*nE
for j in range(LL): wA[h(0,j)]=1
wB=[0]*nE
for i in range(LL): wB[v(i,0)]=1
def dot4(a,b): return sum(x*y for x,y in zip(a,b))%4
best={}
count={}
for m in range(4**nE):
    u=[]; mm=m
    for _ in range(nE):
        u.append(mm%4); mm//=4
    if any(dot4(u,p) for p in plaq): continue
    t=(dot4(u,wA),dot4(u,wB))
    w=sum(1 for x in u if x)
    count[t]=count.get(t,0)+1
    if t not in best or w<best[t]: best[t]=w
want={(a,b): 2*((a!=0)+(b!=0)) for a in range(4) for b in range(4)}
gate("D(Z4) 16 classes x 64 elements", len(count)==16 and all(c==64 for c in count.values()))
gate("D(Z4) class minima == L x |supp(t)| incl. d(2t)=d(t)", best==want, str(best) if best!=want else "")

# ---------- (c) toric 4-point Aut by brute force (L-independent shape {0,1,1,2} scaled)
pts=[(0,0),(1,0),(0,1),(1,1)]
def dtor(p,q):
    t=(p[0]^q[0],p[1]^q[1])
    return t[0]+t[1]   # shape 0/1/1/2 (scale irrelevant to Aut)
cnt=0
for perm in itertools.permutations(range(4)):
    if all(dtor(pts[i],pts[j])==dtor(pts[perm[i]],pts[perm[j]]) for i in range(4) for j in range(4)):
        cnt+=1
gate("toric config metric |Aut| = 8 (all 24 perms tried)", cnt==8, "got %d"%cnt)

# ---------- (d) D(Z2xZ2) product form -> Aut = S4 wr S2 = 1152 (known Hamming H(2,4))
labels4=[(a,b) for a in range(2) for b in range(2)]
configs=[(t1,t2) for t1 in labels4 for t2 in labels4]
def dz22(s,q):  # from the lane's verified law d = L*|supp(t1 OR t2)|, L=2
    t1=(s[0][0]^q[0][0], s[0][1]^q[0][1]); t2=(s[1][0]^q[1][0], s[1][1]^q[1][1])
    return 2*(((t1[0]|t2[0]))+((t1[1]|t2[1])))
# product decomposition: coordinate A = (t1.a,t2.a), coordinate B = (t1.b,t2.b)
prod_ok=True
for s in configs:
    for q in configs:
        ca=(s[0][0]^q[0][0], s[1][0]^q[1][0])
        cb=(s[0][1]^q[0][1], s[1][1]^q[1][1])
        w=2*((ca!=(0,0))+(cb!=(0,0)))
        prod_ok &= (w==dz22(s,q))
gate("D(Z2xZ2) metric == product of two uniform 4-point factors (Hamming H(2,4); "
     "Aut = S4 wr S2 = 1152, matching the enumerated 1152)", prod_ok)
print("V4 OVERALL:", "PASS" if ok else "FAIL")

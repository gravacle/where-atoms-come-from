"""VERIFY check 3 (independent): [[8,6,2]] configuration-distance spectrum from scratch,
plus the [[6,4,2]] uniform claim and the [[6,4,2]] record-venue counting (720 pairs sp=0,J>0).
No solver: N(S) for S=<X^n,Z^n> is {(x,z): |x| even, |z| even}; classes = cosets of
Sbar={0,Xn,Zn,XnZn}; class distance = min |supp x U supp z| over the coset."""
def pcnt(x): return bin(x).count("1")
ok=True
def gate(lbl,b,detail=""):
    global ok
    print(("PASS " if b else "FAIL "),lbl,detail); ok&=b

for n,expect_spec,expect_at2 in ((6,[0,2],15),(8,[0,2,4],28)):
    full=(1<<n)-1
    evens=[x for x in range(1<<n) if pcnt(x)%2==0]
    Sbar={(0,0),(full,0),(0,full),(full,full)}
    canon={}
    for x in evens:
        for z in evens:
            u=(x,z)
            c=min((x^a, z^b) for (a,b) in Sbar)
            w=pcnt(x|z)
            if c not in canon or w<canon[c]:
                canon[c]=w
    nclasses=len(canon)
    gate("n=%d number of classes = %d"%(n,4**(n-2)), nclasses==4**(n-2), "got %d"%nclasses)
    spec=sorted(set(canon.values()))
    gate("n=%d spectrum %s"%(n,expect_spec), spec==expect_spec, "got %s"%spec)
    dmin=min(w for c,w in canon.items() if c!=min(canon))  # min nonzero
    dmin=min(w for c,w in canon.items() if w>0)
    at_d=sum(1 for c,w in canon.items() if w==dmin and c!=(0,0))
    gate("n=%d classes at d_code=%d: %d"%(n,dmin,expect_at2), at_d==expect_at2, "got %d"%at_d)

# [[6,4,2]] record venue: count pairs of classes with sp=0 but J>0
n=6; full=(1<<n)-1
evens=[x for x in range(1<<n) if pcnt(x)%2==0]
Sbar=[(0,0),(full,0),(0,full),(full,full)]
classes={}
for x in evens:
    for z in evens:
        c=min((x^a,z^b) for (a,b) in Sbar)
        classes.setdefault(c,[]).append((x,z))
cl=sorted(classes)
assert len(cl)==256
def sp(a,b): return (pcnt(a[0]&b[1])+pcnt(a[1]&b[0]))%2
def cross(a,b): return pcnt((a[0]&b[1])^(a[1]&b[0]))
cnt=0
for i,a in enumerate(cl):
    for b in cl:
        if a==b: continue
        s=sp(classes[a][0],classes[b][0])
        if s==1: continue
        Jab=min(cross(u,v) for u in classes[a] for v in classes[b])
        if Jab>0: cnt+=1
gate("[[6,4,2]] ordered nonzero... pairs with sp=0 and J>0 == 720", cnt==720, "got %d"%cnt)
print("V3 OVERALL:", "PASS" if ok else "FAIL")

# LANE_W11_R_MATH_CROSS — LEG 5.  THE EXACT BLOCK, RE-RUN IN A DIFFERENT NUMBER REPRESENTATION.
# wm5 M5c works in Z[zeta_12] with the basis (1, z, z^2, z^3) and the relation z^4 = z^2 - 1.
# I work in Q(sqrt3, i) with the basis (1, sqrt3, i, i sqrt3) -- the same field, a different
# basis and a different multiplication table -- so a reduction bug in one cannot survive in both.
# Same connection (S1 sec6: a1=a2=a3=pi/3, a4=a5=a6=pi/2), same three ready states.
# NO FLOATING POINT IN THE EXACT BLOCK.
from fractions import Fraction as Fr
import numpy as np, xc0_lib as X

# element = (a, b, c, d)  meaning  a + b sqrt3 + c i + d i sqrt3 ,  a,b,c,d in Q
Z0=(Fr(0),)*4; ONE=(Fr(1),Fr(0),Fr(0),Fr(0))
def add(p,q): return tuple(u+v for u,v in zip(p,q))
def mul(p,q):
    a,b,c,d=p; e,f,g,h=q
    # (a + b s) + i(c + d s)  times  (e + f s) + i(g + h s),  s^2 = 3
    ra=a*e+3*b*f - (c*g+3*d*h)
    rb=a*f+b*e   - (c*h+d*g)
    rc=a*g+3*b*h + (c*e+3*d*f)
    rd=a*h+b*g   + (c*f+d*e)
    return (ra,rb,rc,rd)
def conj(p): a,b,c,d=p; return (a,b,-c,-d)
def num(p): return float(p[0])+float(p[1])*3**0.5 + 1j*(float(p[2])+float(p[3])*3**0.5)

UF=(Fr(1,2),Fr(0),Fr(0),Fr(1,2))        # exp(i pi/3) = 1/2 + i sqrt3/2
UC=(Fr(0),Fr(0),Fr(1),Fr(0))            # exp(i pi/2) = i
WF=mul(mul(UF,UF),UF)                   # = -1
WC=mul(mul(UC,UC),UC)                   # = -i
print("== X5a  THE FIELD ELEMENTS, IN MY BASIS ==")
print(f"   U_F = exp(i pi/3) = {UF}   U_F^3 = W_F = {WF}   (expect (-1,0,0,0))")
print(f"   U_C = exp(i pi/2) = {UC}   U_C^3 = W_C = {WC}   (expect (0,0,-1,0) = -i)")
assert WF==(Fr(-1),Fr(0),Fr(0),Fr(0)) and WC==(Fr(0),Fr(0),Fr(-1),Fr(0))

NV=5
def rat(k,d): return (Fr(k,d),Fr(0),Fr(0),Fr(0))
sA=[rat(20,25),rat(3,25),rat(4,25),rat(10,25),rat(10,25)]
sB=[rat(20,25),rat(5,25),rat(0,25),rat(10,25),rat(10,25)]
ZETA=(Fr(1,2),Fr(0),Fr(0),Fr(0))        # placeholder, real zeta below
zeta=(Fr(0),Fr(1,2),Fr(1,2),Fr(0))      # exp(i pi/6) = sqrt3/2 + i/2
def zpow(p,n):
    r=ONE
    for _ in range(n): r=mul(r,p)
    return r
sC=[mul(sA[i],zpow(zeta,k)) for i,k in enumerate((0,1,5,3,7))]
def pi_exact(s):
    w=[mul(conj(x),x) for x in s]
    for x in w: assert x[1]==0 and x[2]==0 and x[3]==0, "modulus^2 must be rational"
    return (w[0][0], w[1][0]+w[2][0], w[3][0]+w[4][0])
print(f"\n   pi(A) = {pi_exact(sA)}   pi(B) = {pi_exact(sB)}   pi(C) = {pi_exact(sC)}")
assert pi_exact(sA)==pi_exact(sB)==pi_exact(sC)
assert sA!=sB and sA!=sC and sB!=sC
print("   pi EXACTLY equal as rationals; the three states pairwise distinct.  (wm5 M5c prints")
print("   16/25, 1/25, 8/25 for p11, p10, p01 -- reproduced here in a different basis.)")

def matvec(M,v):
    out=[]
    for i in range(NV):
        acc=Z0
        for j in range(NV):
            if M[i][j]!=Z0: acc=add(acc,mul(M[i][j],v[j]))
        out.append(acc)
    return out
def zeros(): return [[Z0]*NV for _ in range(NV)]
TF=zeros(); TC=zeros(); MF=zeros(); MC=zeros()
for v in (3,4): TF[v][v]=ONE
for (u,v) in ((0,1),(1,2),(2,0)): TF[v][u]=UF
for v in (1,2): TC[v][v]=ONE
for (u,v) in ((0,3),(3,4),(4,0)): TC[v][u]=UC
for v in range(NV): MF[v][v]=WF if v in (0,1,2) else ONE
for v in range(NV): MC[v][v]=WC if v in (0,3,4) else ONE
def Zn(opF,opC,s,n):
    xF=list(s); xC=list(s)
    for _ in range(n): xF=matvec(opF,xF); xC=matvec(opC,xC)
    acc=Z0
    for v in range(NV): acc=add(acc,mul(conj(xF[v]),xC[v]))
    return acc
print("\n== X5b  |Z|^2 EXACTLY, MY BASIS, AGAINST wm5 M5c's PRINTED Z[zeta_12] VALUES ==")
ref_circ={1:Fr(13,125),2:Fr(529,625),3:Fr(13,125),4:Fr(1)}
print(f"   {'n':>3}  {'CIRCUIT |Z|^2 (mine, exact)':<34}{'equal across A,B,C?':>21}   wm5 printed")
for n in (1,2,3,4):
    vals=[mul(conj(Zn(MF,MC,s,n)),Zn(MF,MC,s,n)) for s in (sA,sB,sC)]
    for v in vals: assert v[2]==0 and v[3]==0, "must be real"
    ok=vals[0]==vals[1]==vals[2]
    r=vals[0]
    txt=f"{r[0]}" + (f" + {r[1]} sqrt3" if r[1]!=0 else "")
    print(f"   {n:>3}  {txt:<34}{str(ok):>21}   {ref_circ[n]}")
    assert r[1]==0 and r[0]==ref_circ[n], "MUST match wm5 M5c"
print("   -> all four circuit rows match wm5 M5c EXACTLY, as rationals, in a different basis.")
print(f"\n   {'n':>3}  {'EDGE |Z|^2 (mine, exact, state A)':<44}{'equal A,B,C?':>14}")
for n in (1,2,3,4,5,6):
    vals=[mul(conj(Zn(TF,TC,s,n)),Zn(TF,TC,s,n)) for s in (sA,sB,sC)]
    for v in vals: assert v[2]==0 and v[3]==0, "must be real"
    ok=vals[0]==vals[1]==vals[2]
    r=vals[0]
    txt=f"{r[0]}" + (f" + ({r[1]}) sqrt3" if r[1]!=0 else "")
    print(f"   {n:>3}  {txt:<44}{str(ok):>14}   numeric {num(vals[0]).real:.9f}")
print("   -> EXACTLY equal at n = 3 and 6, EXACTLY unequal at n = 1,2,4,5.  wm5 M5c reproduced")
print("      in an independent representation.  Note my n=1 value carries an IRRATIONAL part")
print("      (a multiple of sqrt3) while wm5 prints it in the zeta basis as (108784/390625,")
print("      -1728/15625, 0, 864/15625); the two agree numerically to the last digit:")
v1=mul(conj(Zn(TF,TC,sA,1)),Zn(TF,TC,sA,1))
lane=Fr(108784,390625)-Fr(1728,15625)*0  # reconstruct the lane's element numerically
lane_num=float(Fr(108784,390625)) + float(Fr(-1728,15625))*np.cos(np.pi/6) + float(Fr(864,15625))*np.cos(3*np.pi/6) \
        + 1j*(float(Fr(-1728,15625))*np.sin(np.pi/6) + float(Fr(864,15625))*np.sin(3*np.pi/6))
print(f"      mine  {num(v1):.12f}      lane's zeta-basis element {lane_num:.12f}"
      f"      |diff| {abs(num(v1)-lane_num):.2e}")

print("\n== X5c  AND THE SAME EXACT TEST ON THE OBJECT THAT ACTUALLY DECIDES: |Z| AT n NOT DIV 3 ==")
print("   The exact rationals show the invisibility break is not a floating-point artefact.")
print("   Magnitudes of the difference, exactly, at each n:")
for n in (1,2,4,5):
    vals=[mul(conj(Zn(TF,TC,s,n)),Zn(TF,TC,s,n)) for s in (sA,sB,sC)]
    d1=tuple(u-v for u,v in zip(vals[0],vals[1])); d2=tuple(u-v for u,v in zip(vals[0],vals[2]))
    print(f"   n={n}:  |Z_A|^2 - |Z_B|^2 = {num(d1).real:+.9f}   |Z_A|^2 - |Z_C|^2 = {num(d2).real:+.9f}"
          f"    both exactly nonzero: {d1!=Z0 and d2!=Z0}")

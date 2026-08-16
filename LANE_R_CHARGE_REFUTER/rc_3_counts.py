"""Supplement: the 243-charge enumeration at a quadrature accurate enough for the tolerance."""
import itertools, numpy as np, rc_lib as R
Ecls = R.class_vectors(5, R.K1_EDGES, R.K1_GAMMA_F, R.K1_GAMMA_C)
p_S3 = np.array([0.4,0.15,0.15,0.15,0.15])
def Ech(q): return (Ecls.T*np.asarray(q,dtype=np.int64)).T
base = R.lambda_B_generic(Ecls,p_S3,Nx=32768)
same=0; rank={0:0,1:0,2:0}; vals={}
for q in itertools.product(range(3),repeat=5):
    E=Ech(q); b=R.delta_lattice(E,p_S3); rank[len(b)]+=1
    lam=R.lambda_B_generic(E,p_S3,Nx=8192)
    if abs(lam-base)<1e-7: same+=1
    vals[round(lam,8)]=vals.get(round(lam,8),0)+1
print(f"base (unit charge) = {base:.12f}")
print(f"243 assignments q in {{0,1,2}}^5:  rank Delta = 2:{rank[2]}  1:{rank[1]}  0:{rank[0]}")
print(f"lambda_B^gen EQUAL to the unit-charge value (tol 1e-7, Nx=8192): {same} of 243")
print(f"distinct values: {len(vals)}")
# class-homogeneous locus, exhaustive over {1,2,3}^3 and with zeros allowed
ch=0; chsame=0
for q0,q1,q3 in itertools.product(range(0,4),repeat=3):
    q=[q0,q1,q1,q3,q3]; ch+=1
    lam=R.lambda_B_generic(Ech(q),p_S3,Nx=8192)
    if abs(lam-base)<1e-7: chsame+=1
print(f"class-homogeneous charges (q0,q1,q3) in {{0..3}}^3: {ch} total, "
      f"{chsame} give EXACTLY the unit-charge rate")

pres=[]
for q in itertools.product(range(3),repeat=5):
    lam=R.lambda_B_generic(Ech(q),p_S3,Nx=8192)
    if abs(lam-base)<1e-7: pres.append(q)
homog=[q for q in pres if q[1]==q[2] and q[3]==q[4]]
print(f"of the {len(pres)} rate-preserving charges, class-homogeneous (q1=q2 and q3=q4): {len(homog)}")
print("preserving set:", pres)

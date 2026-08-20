"""O-50-B shared machinery.  Nothing here nominates anything; every writer is SEARCHED for
   or DECIDED by an exact criterion, per D-18."""
import numpy as np, itertools, sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import RecordModel, symplectic_logicals, xz_to_matrix, eigenspaces

TOL = 1e-9
I2 = np.eye(2); X = np.array([[0,1],[1,0]],dtype=complex)
Z = np.array([[1,0],[0,-1]],dtype=complex); Y = 1j*X@Z

def kron(*ms):
    M = np.array([[1]],dtype=complex)
    for m in ms: M = np.kron(M,m)
    return M

def pauli_on(n, spec):
    """spec: dict qubit->'X'/'Y'/'Z'"""
    d = {'X':X,'Y':Y,'Z':Z,'I':I2}
    return kron(*[d.get(spec.get(i,'I')) for i in range(n)])

# ---------------------------------------------------------------- EXACT commutant (avoids D-21)
def exact_commutant(gens, tol=1e-8):
    """Basis of {M : [M,g]=0 for all g}, by an EXACT nullspace of the stacked commutator map.
       NOT sampling based.  record_model.commutant() is sampling based and returns a SHORT
       basis (D-21); this routine is used instead wherever a commutant dimension matters."""
    n = gens[0].shape[0]
    rows = []
    for g in gens:
        # vec(gM - Mg) = (I kron g - g^T kron I) vec(M)
        rows.append(np.kron(np.eye(n), g) - np.kron(g.T, np.eye(n)))
    A = np.vstack(rows)
    u,s,vh = np.linalg.svd(A)
    ns = vh[np.sum(s > tol*max(1.0,s[0] if len(s) else 1.0)):].conj().T if len(s) else np.eye(n*n)
    k = ns.shape[1]
    return [ns[:,j].reshape(n,n) for j in range(k)]

# ---------------------------------------------------------------- clause checks (D-18)
def clause_i(R):   return np.linalg.norm(R-R.conj().T)<1e-8 and np.linalg.norm(R@R-np.eye(R.shape[0]))<1e-8
def clause_ii(R,H,Ls=()):
    ok = np.linalg.norm(R@H-H@R)<1e-8
    for L in Ls: ok = ok and np.linalg.norm(R@L-L@R)<1e-8
    return ok
def clause_iii_(R,es):
    for _,P,m in es:
        M=P@R@P
        if np.linalg.norm(M-(np.trace(M)/m)*P)>1e-8: return True
    return False
def clause_iv_(R,es): return all(abs(np.trace(P@R))<1e-7 for _,P,_ in es)
def imbalance(R,es):  return max(abs(np.trace(P@R)) for _,P,_ in es)

# ---------------------------------------------------------------- the dimension table
def joint_blocks(H, family, tol=1e-7):
    """{(energy_index, sign_tuple): orthonormal column block}.  All members of `family`
       must commute with H and with each other."""
    es = eigenspaces(H)
    out = {}
    for ei,(_,PE,m) in enumerate(es):
        w,V = np.linalg.eigh(PE); C = V[:, w>0.5]
        groups = {(): C}
        for R in family:
            ng = {}
            for lab,Cc in groups.items():
                Rs = Cc.conj().T@R@Cc
                ws,Vs = np.linalg.eigh(Rs)
                for s in (+1,-1):
                    idx=[i for i in range(len(ws)) if (ws[i]>0)==(s>0)]
                    if idx: ng[lab+(s,)] = Cc@Vs[:,idx]
            groups = ng
        for lab,Cc in groups.items(): out[(ei,lab)] = Cc
    return es, out

def dim_table(H, family):
    es, blocks = joint_blocks(H, family)
    k = len(family)
    cfg = list(itertools.product((1,-1),repeat=k))
    d = {}
    for ei in range(len(es)):
        for s in cfg: d[(ei,s)] = 0
    for (ei,lab),C in blocks.items(): d[(ei,lab)] = C.shape[1]
    return es, blocks, d, cfg

# ---------------------------------------------------------------- EXACT writer criterion
def realisable_flips(d, cfg, nE):
    """LEMMA (proved in t2): an ADMISSIBLE U with U*R_i U = eps_i R_i for all i exists
       IFF dim V_(E,sigma) = dim V_(E,eps.sigma) for every energy shell E and every sigma.
       Returns the full group of realisable flip patterns eps in (Z_2)^k."""
    k = len(cfg[0])
    good = []
    for m in range(1<<k):
        eps = tuple((m>>i)&1 for i in range(k))
        ok = True
        for ei in range(nE):
            for s in cfg:
                t = tuple(-x if e else x for e,x in zip(eps,s))
                if d[(ei,s)] != d[(ei,t)]: ok=False; break
            if not ok: break
        if ok: good.append(eps)
    return good

def build_flip_unitary(blocks, eps, nE, cfg):
    """the block permutation implementing eps.  Manifestly unitary, manifestly [U,H]=0
       (it never leaves an energy shell).  Built, then VERIFIED, never assumed."""
    n = None
    for C in blocks.values(): n = C.shape[0]; break
    U = np.zeros((n,n),dtype=complex); seen=set()
    for ei in range(nE):
        for s in cfg:
            if (ei,s) in seen: continue
            t = tuple(-x if e else x for e,x in zip(eps,s))
            Cs = blocks.get((ei,s)); Ct = blocks.get((ei,t))
            if Cs is None and Ct is None: seen.add((ei,s)); seen.add((ei,t)); continue
            if Cs is None or Ct is None or Cs.shape[1]!=Ct.shape[1]: return None
            if s==t: U += Cs@Cs.conj().T
            else:    U += Cs@Ct.conj().T + Ct@Cs.conj().T
            seen.add((ei,s)); seen.add((ei,t))
    return U

def orbits_of(G, cfg):
    seen=set(); orbs=[]
    for s in cfg:
        if s in seen: continue
        o = sorted({tuple(-x if e else x for e,x in zip(g,s)) for g in G})
        for x in o: seen.add(x)
        orbs.append(o)
    return orbs

def invariant_characters(G, k):
    out=[]
    for m in range(1<<k):
        S = tuple(i for i in range(k) if (m>>i)&1)
        if all(sum(g[i] for i in S)%2==0 for g in G): out.append(S)
    return out

def f2rank(G,k):
    basis=[]
    for g in G:
        v=list(g)
        for b in basis:
            h=next((i for i in range(k) if b[i]),None)
            if h is not None and v[h]: v=[x^y for x,y in zip(v,b)]
        if any(v): basis.append(v)
    return len(basis)

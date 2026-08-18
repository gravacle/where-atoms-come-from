"""Stabiliser-model builders shared by the O-7 and O-8 scripts."""
import sys
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_O78_FRACTON_NONABELIAN")
from gf2 import *

# ---------------------------------------------------------------- builders
def toric2d(L):
    """qubits on edges of L x L torus. n = 2L^2. k = 2, d = L."""
    n = 2*L*L
    q = lambda x,y,d: 2*((x%L)*L + (y%L)) + d
    g = []
    for x in range(L):
        for y in range(L):
            # vertex X-star
            xs = [q(x,y,0), q(x-1,y,0), q(x,y,1), q(x,y-1,1)]
            v = 0
            for i in xs: v ^= (1 << i)
            g.append(mk(v,0,n))
            # plaquette Z
            zs = [q(x,y,0), q(x,y+1,0), q(x,y,1), q(x+1,y,1)]
            v = 0
            for i in zs: v ^= (1 << i)
            g.append(mk(0,v,n))
    return g, n

def toric3d(L):
    """qubits on edges of L^3 torus. n = 3L^3. k = 3 (H_1(T^3) = F_2^3)."""
    n = 3*L*L*L
    idx = lambda x,y,z: ((x%L)*L + (y%L))*L + (z%L)
    q = lambda x,y,z,d: 3*idx(x,y,z) + d
    e = [(1,0,0),(0,1,0),(0,0,1)]
    g = []
    for x in range(L):
        for y in range(L):
            for z in range(L):
                # vertex X-star: 6 edges
                v = 0
                for d in range(3):
                    v ^= (1 << q(x,y,z,d))
                    v ^= (1 << q(x-e[d][0], y-e[d][1], z-e[d][2], d))
                g.append(mk(v,0,n))
                # plaquettes: 3 per site
                for (d1,d2) in [(0,1),(1,2),(2,0)]:
                    a = e[d1]; b = e[d2]
                    v = 0
                    v ^= (1 << q(x,y,z,d1))
                    v ^= (1 << q(x+b[0],y+b[1],z+b[2],d1))
                    v ^= (1 << q(x,y,z,d2))
                    v ^= (1 << q(x+a[0],y+a[1],z+a[2],d2))
                    g.append(mk(0,v,n))
    return g, n

def xcube(L):
    """X-cube (Vijay-Haah-Fu). qubits on edges of L^3 torus, n = 3L^3.
       cube operator = X on the 12 edges; vertex cross = Z on 4 coplanar edges."""
    n = 3*L*L*L
    idx = lambda x,y,z: ((x%L)*L + (y%L))*L + (z%L)
    q = lambda x,y,z,d: 3*idx(x,y,z) + d
    e = [(1,0,0),(0,1,0),(0,0,1)]
    g = []
    for x in range(L):
        for y in range(L):
            for z in range(L):
                # cube with low corner (x,y,z): 12 edges, X-type
                v = 0
                for d in range(3):
                    u,w = [k for k in range(3) if k != d]
                    for a in (0,1):
                        for b in (0,1):
                            p = [x,y,z]
                            p[u] += a; p[w] += b
                            v ^= (1 << q(p[0],p[1],p[2],d))
                g.append(mk(v,0,n))
                # vertex crosses, Z-type, planes xy, yz, zx (third is the product of the other two)
                for (d1,d2) in [(0,1),(1,2)]:
                    v = 0
                    for d in (d1,d2):
                        v ^= (1 << q(x,y,z,d))
                        v ^= (1 << q(x-e[d][0], y-e[d][1], z-e[d][2], d))
                    g.append(mk(0,v,n))
    return g, n

def checkerboard(L):
    """Checkerboard model. qubits on SITES of L^3 torus (L even). n = L^3.
       On each cube of one colour: X on all 8 corners AND Z on all 8 corners."""
    assert L % 2 == 0
    n = L*L*L
    idx = lambda x,y,z: ((x%L)*L + (y%L))*L + (z%L)
    g = []
    for x in range(L):
        for y in range(L):
            for z in range(L):
                if (x+y+z) % 2: continue          # one colour of cubes only
                v = 0
                for a in (0,1):
                    for b in (0,1):
                        for c in (0,1):
                            v ^= (1 << idx(x+a,y+b,z+c))
                g.append(mk(v,0,n)); g.append(mk(0,v,n))
    return g, n

# Haah cubic code 1 as a polynomial (Laurent) pair over F_2[x,y,z]/(x^L-1,...)
HAAH_A = [(0,0,0),(1,0,0),(0,1,0),(0,0,1)]           # 1 + x + y + z
HAAH_B = [(0,0,0),(1,1,0),(0,1,1),(1,0,1)]           # 1 + xy + yz + zx

def haah(L):
    """Haah cubic code 1: 2 qubits per site of L^3 torus, n = 2L^3.
       sigma_Z = (a, b);  sigma_X = (bbar, abar)  -> CSS condition holds identically."""
    n = 2*L*L*L
    idx = lambda x,y,z: ((x%L)*L + (y%L))*L + (z%L)
    q = lambda x,y,z,s: 2*idx(x,y,z) + s
    g = []
    for x in range(L):
        for y in range(L):
            for z in range(L):
                vz = 0
                for (a,b,c) in HAAH_A: vz ^= (1 << q(x+a,y+b,z+c,0))
                for (a,b,c) in HAAH_B: vz ^= (1 << q(x+a,y+b,z+c,1))
                g.append(mk(0,vz,n))
                vx = 0
                for (a,b,c) in HAAH_B: vx ^= (1 << q(x-a,y-b,z-c,0))
                for (a,b,c) in HAAH_A: vx ^= (1 << q(x-a,y-b,z-c,1))
                g.append(mk(vx,0,n))
    return g, n

def chamon(L, sublattice=False):
    """Chamon model: 1 qubit per site of L^3 torus.
       O_r = X_{r+ex} X_{r-ex} Y_{r+ey} Y_{r-ey} Z_{r+ez} Z_{r-ez}.  NON-CSS."""
    n = L*L*L
    idx = lambda x,y,z: ((x%L)*L + (y%L))*L + (z%L)
    g = []
    for x in range(L):
        for y in range(L):
            for z in range(L):
                if sublattice and (x+y+z) % 2: continue
                X = 0; Z = 0
                for s in (1,-1):
                    X ^= (1 << idx(x+s,y,z))                       # X
                    X ^= (1 << idx(x,y+s,z)); Z ^= (1 << idx(x,y+s,z))   # Y = X and Z
                    Z ^= (1 << idx(x,y,z+s))                       # Z
                g.append(mk(X,Z,n))
    return g, n

def steane():
    n = 7
    Hs = ["0001111","0110011","1010101"]
    g = []
    for h in Hs:
        v = 0
        for i,c in enumerate(h):
            if c=="1": v |= (1<<i)
        g.append(mk(v,0,n)); g.append(mk(0,v,n))
    return g, n

def perfect5():
    """[[5,1,3]] : XZZXI and cyclic shifts."""
    n = 5
    base = "XZZXI"
    g = []
    for s in range(4):
        w = base[-s:] + base[:-s] if s else base
        X = 0; Z = 0
        for i,c in enumerate(w):
            if c in "XY": X |= (1<<i)
            if c in "ZY": Z |= (1<<i)
        g.append(mk(X,Z,n))
    return g, n


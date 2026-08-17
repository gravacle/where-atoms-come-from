# carriers.py -- the LADDER OF CARRIERS for LANE W19-A.
# Each entry: name -> (V, edges).  Edges are directed (tail, head); multi-edges allowed.
from collections import deque

def theta():                       # T1's carrier. 2 vertices, 3 parallel links.
    return 2, [(0,1),(0,1),(0,1)]

def theta_sub(k):                  # theta with each link subdivided k times (degree-2 vertices).
    edges=[]; V=2
    for _ in range(3):
        prev=0
        for j in range(k):
            edges.append((prev,V)); prev=V; V+=1
        edges.append((prev,1))
    return V, edges

def ladder(nsq):                   # chain of nsq squares. V=2(nsq+1), L=3nsq+1, C=nsq.
    n=nsq+1; V=2*n; edges=[]
    for i in range(n-1): edges.append((i, i+1))          # bottom rail
    for i in range(n-1): edges.append((n+i, n+i+1))      # top rail
    for i in range(n):   edges.append((i, n+i))          # rungs
    return V, edges

def cube():                        # Q3, cubic, girth 4. V=8 L=12 C=5
    V=8; edges=[]
    for x in range(8):
        for bit in range(3):
            y = x ^ (1<<bit)
            if y > x: edges.append((x,y))
    return V, edges

def grid(nx, ny, periodic):        # square lattice, nx*ny vertices
    V=nx*ny; idx=lambda i,j: i*ny+j; edges=[]
    for i in range(nx):
        for j in range(ny):
            if periodic or i+1<nx: edges.append((idx(i,j), idx((i+1)%nx, j)))
            if periodic or j+1<ny: edges.append((idx(i,j), idx(i, (j+1)%ny)))
    return V, edges

def petersen():                    # cubic, girth 5, the (3,5)-cage. V=10 L=15 C=6
    V=10; edges=[]
    for i in range(5): edges.append((i,(i+1)%5))
    for i in range(5): edges.append((5+i, 5+((i+2)%5)))
    for i in range(5): edges.append((i, 5+i))
    return V, edges

def lcf(n, shifts):
    edges=[(i,(i+1)%n) for i in range(n)]
    seen=set(frozenset(e) for e in edges)
    for i in range(n):
        j=(i+shifts[i%len(shifts)])%n
        k=frozenset((i,j))
        if k not in seen: seen.add(k); edges.append((i,j))
    return n, edges

def heawood():                     # LCF [5,-5]^7. cubic, girth 6, the (3,6)-cage.
    return lcf(14,[5,-5])          # == hexagonal (honeycomb) tiling of the torus, 7 hexagons

def mobius_kantor():               # LCF [5,-5]^8. cubic, girth 6. V=16 L=24 -- OVER THE CEILING.
    return lcf(16,[5,-5])

def pappus():                      # LCF [5,7,-7,7,-7,-5]^3. cubic girth 6. V=18 -- over ceiling.
    return lcf(18,[5,7,-7,7,-7,-5])

LADDER = [
  ("theta_3link",        theta()),
  ("theta_subdiv1",      theta_sub(1)),
  ("theta_subdiv2",      theta_sub(2)),
  ("torus_2x2",          grid(2,2,True)),
  ("ladder_2sq",         ladder(2)),
  ("ladder_3sq",         ladder(3)),
  ("ladder_5sq",         ladder(5)),
  ("grid_3x3_open",      grid(3,3,False)),
  ("cube_Q3",            cube()),
  ("torus_3x3",          grid(3,3,True)),
  ("petersen",           petersen()),
  ("heawood_honeycomb7", heawood()),
]

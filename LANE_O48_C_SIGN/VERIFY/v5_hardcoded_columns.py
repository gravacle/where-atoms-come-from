"""
VERIFY-5.  TWO PRINTED COLUMNS IN THE LANE'S TABLES ARE CONSTANTS, NOT MEASUREMENTS.

(A) s1_clauses.py, clause_table():
        c_i  = True
        if n <= 5:  ...dense check...
        c_ii = True
        if n <= 5:  ...dense check...
    The reported table runs n = 3,4,6,8,10,12 for 5 families = 30 rows, and prints "True True"
    in all 30.  Only the n=3 and n=4 rows (10 of 30) are ever tested.  The other 20 print the
    initialiser.  The finding's own text says clauses (i)-(iv) were "verified EXPLICITLY on the
    carrier at every n reported, never assumed".

    Demonstration: feed the SAME code path an operator that manifestly FAILS clause (ii)
    (R = X_0, which anticommutes with the bond terms touching site 0) and watch the column.

(B) s1_clauses.py, writer table:  `dE = 0 if cnt else None`  -- the energy-change column is a
    literal, not a measurement.  (Measured independently in v4_ordinary.py: it is 0, but it is 0
    for every admissible W and every H, by [W,H]=0.)
"""
import itertools
import numpy as np

I2=np.eye(2,dtype=complex); X2=np.array([[0,1],[1,0]],dtype=complex)
Z2=np.array([[1,0],[0,-1]],dtype=complex)
def kron(o):
    r=np.array([[1.0+0j]])
    for a in o: r=np.kron(r,a)
    return r
def dense_H(J,n):
    H=np.zeros((2**n,2**n),dtype=complex)
    for i,j in enumerate(J):
        ops=[I2]*n; ops[i]=Z2; ops[i+1]=Z2; H+=j*kron(ops)
    return H

def lane_clause_ii_column(J,n,op_at_site_0):
    """EXACTLY the lane's control flow for the (ii) column."""
    c_ii = True                      # <-- initialiser
    if n <= 5:                       # <-- gate
        Hd = dense_H(J,n)
        ops=[I2]*n; ops[0]=op_at_site_0; O=kron(ops)
        c_ii &= np.linalg.norm(Hd@O - O@Hd) < 1e-9
    return c_ii

def truth_clause_ii(J,n,op_at_site_0):
    Hd = dense_H(J,n)
    ops=[I2]*n; ops[0]=op_at_site_0; O=kron(ops)
    return bool(np.linalg.norm(Hd@O - O@Hd) < 1e-9)

print("="*104)
print("VERIFY-5   THE (i) AND (ii) COLUMNS ARE DEFAULT-TRUE ABOVE n = 5")
print("="*104)
print()
print(f"  {'operator':<12} {'n':>3} {'lane prints (ii)':>18} {'truth':>8} {'agree?':>8}   note")
for op,name in ((Z2,"R = Z_0"),(X2,"R = X_0")):
    for n in (3,4,5,6,8,10,12):
        J=[i+1 for i in range(n-1)]
        if n>10 and name=="R = X_0": pass
        lane=lane_clause_ii_column(J,n,op)
        tru=truth_clause_ii(J,n,op) if n<=12 else None
        note = "" if lane==tru else "  <-- LANE PRINTS True FOR A FALSE CLAUSE"
        print(f"  {name:<12} {n:>3} {str(lane):>18} {str(tru):>8} {str(lane==tru):>8}{note}")
    print()
print("  READ: above n = 5 the lane's (ii) column returns True for ANY operator, including one")
print("  that does not commute with H at all.  In the reported 30-row clause table, 20 rows")
print("  (n = 6, 8, 10, 12 x 5 families) carry that untested True for clauses (i) and (ii).")
print("  Clauses (iii) and (iv) ARE computed at every n; (i) and (ii) are not.")

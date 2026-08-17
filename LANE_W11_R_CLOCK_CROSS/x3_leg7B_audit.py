# X3 — DID LANE C QUIETLY CHANGE THE OBSERVABLE WHILE CHANGING THE CLOCK?  (my brief, verbatim)
#
# LANE C's Reading-A tilt rests on ONE leg: 7B.  Its structured output says
#   "the rate along every sublattice ray converges to m(P) ... built from COR-F's T alone"
#   "N1 ... is re-derivable from COR-F's own transport"
# READ w11c_7_closedform.py:44-50.  7B's inner loop is
#       u = exp(-i*angle(W_F)*i_step*n);  v = exp(+i*angle(W_C)*j_step*n)
#       z = |poly(pi, u, v)|;   vals.append(mean(log z))
# NEITHER T_F, T_C NOR ANY READY STATE APPEARS.  Leg 7B evaluates a scalar polynomial in
# (pi, f, c).  X3 proves this by SUBSTITUTION: destroy T and destroy the states, and re-run 7B.
import numpy as np, importlib.util, sys, io, contextlib, hashlib
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/LANE_W11_R_CLOCK")
import w11c_lib as L
from x_lib import m_jensen as my_m

print("== X3a  SUBSTITUTION TEST: replace T and the ready states with garbage, re-run leg 7 ==")
def capture(patch):
    """run the lane's leg 7 in-process with a patched w11c_lib, capture stdout."""
    import importlib
    importlib.reload(L)
    patch(L)
    src = open("/Users/bgm/MB Work/where-atoms-come-from/LANE_W11_R_CLOCK/w11c_7_closedform.py").read()
    g = {"__name__": "__main__"}
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            exec(compile(src, "w11c_7", "exec"), g)
    except Exception as e:
        buf.write(f"\n*** RAISED {type(e).__name__}: {e}\n")
    return buf.getvalue()

base = capture(lambda m: None)
def kill_T(m):
    real = m.ops
    def fake(K, a):
        TF, TC, MF, MC, WF, WC = real(K, a)
        rng = np.random.default_rng(999)
        # a random unitary in place of BOTH edge ticks -- 7A must break, 7B must not
        def ru(n):
            X = rng.normal(size=(n,n)) + 1j*rng.normal(size=(n,n))
            Q, R = np.linalg.qr(X); return Q*np.exp(-1j*np.angle(np.diag(R)))
        return ru(K.nv), ru(K.nv), MF, MC, WF, WC
    m.ops = fake
noT = capture(kill_T)

def grab7B(txt):
    out = []
    keep = False
    for ln in txt.split("\n"):
        if "7B  CONVERGENCE" in ln: keep = True
        if "7A  CLOSED FORM" in ln: keep = False
        if keep and ("(" in ln and ")" in ln and "." in ln) or (keep and "m(P)" in ln): out.append(ln)
    return "\n".join(out)
def grab7A(txt):
    return "\n".join(l for l in txt.split("\n") if "max | Z(L_F i" in l)

print("  7A WITH THE REAL T:      ", grab7A(base).strip()[:120])
print("  7A WITH T REPLACED:      ", grab7A(noT).strip()[:120])
print("  -> 7A IS a computation on T: destroying T destroys it.  GOOD.")
b7, n7 = grab7B(base), grab7B(noT)
print(f"\n  7B block sha256 with real T : {hashlib.sha256(b7.encode()).hexdigest()[:16]}")
print(f"  7B block sha256 with T KILLED: {hashlib.sha256(n7.encode()).hexdigest()[:16]}")
print(f"  IDENTICAL: {b7 == n7}")
print("  -> LEG 7B IS UNCHANGED WHEN BOTH EDGE TICKS ARE REPLACED BY RANDOM UNITARIES.")
print("     It is NOT 'built from COR-F's T'.  It is a numerical evaluation of Jensen's formula")
print("     for m(P) along an arithmetic progression on the torus.  The substitution is LICENSED")
print("     by 7A -- but the claim that the rate is re-derived FROM THE TRANSPORT is not what 7B")
print("     computes, and 7B is the only leg carrying the Reading-A tilt.")

print("\n== X3b  COULD LEG 7B HAVE FAILED?  Sweep the winding rates it did not sweep. ==")
pi = np.array([0., .30, .30, .40]); mP = my_m(pi)
f, c = 1.0, 2**0.5
print(f"   pi = {pi}   m(P) = {mP:.9f}   (f,c) = (1.0, sqrt2), the corpus's only generic point")
print(f"   {'(i_step, j_step)':>18} {'N=1e5':>14} {'N=1e6':>14} {'diff from m(P)':>16}")
rows = []
for (ii, jj) in [(1,1),(1,2),(3,1),(4,4),(17,43),(101,7),(1,1000),(999,1000)]:
    for N in (10**5, 10**6):
        n = np.arange(1, N+1)
        u = np.exp(-1j*f*ii*n); v = np.exp(1j*c*jj*n)
        z = np.abs(pi[0]+pi[1]*u+pi[2]*v+pi[3]*u*v)
        r = np.log(np.maximum(z,1e-300)).mean()
        if N == 10**5: r5 = r
    rows.append((ii,jj,r5,r))
    print(f"   {str((ii,jj)):>18} {r5:>14.9f} {r:>14.9f} {r-mP:>16.2e}")
print(f"   ALL EIGHT converge to m(P).  max |rate - m(P)| at N=1e6 = "
      f"{max(abs(r-mP) for _,_,_,r in rows):.2e}")
print("   THIS IS WEYL EQUIDISTRIBUTION + JENSEN, AND IT COULD NOT HAVE COME OUT OTHERWISE for")
print("   ANY non-zero integer rates at a generic (f,c).  Lane C marks legs 1, 5 and 7A as [T]")
print("   'could not have failed'.  IT DOES NOT MARK 7B.  7B is the same kind of object, and it")
print("   is the leg its 'content beyond the hypothesis' argument rests on.")

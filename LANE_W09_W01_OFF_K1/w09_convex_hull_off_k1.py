# W-09 leg A — W-01's convex-hull criterion, run in the register's own words on carriers the
# corpus ALREADY OWNS. Class multisets quoted verbatim from S4_THE_MEASUREMENT_V001.md:575-590.
#
# ISOLATION LEDGER. Held fixed: the hull test, the (f,c) grid, the seed, the evaluator, the
# firing predicate. Moved: WHICH CHARACTERS ARE OCCUPIED, and nothing else. That is a property
# of the carrier's incidence alone -- which vertices lie in gamma_F, in gamma_C, in both, in neither.
#
# W-01 (REGISTER_V001.md:43): the overlap "vanishes iff 0 lies in the convex hull of three
# unit-modulus coefficients". Read as the EXISTENCE statement it is: the set of connections at
# which SOME ready state on that carrier fires. (FIRE => HULL is unconditional; HULL => FIRE is
# the direction M3's refuter showed fails pointwise -- quantifier restored here.)
import numpy as np
rng = np.random.default_rng(20260816)
N = 200000

CARRIERS = [                      # name, class multiset as S4 publishes it, S4 line
 ("B1  K1 (as handed)",           {"01":2, "10":2, "11":1},            "S4:578"),
 ("B2  K1 both filled",           {"01":2, "10":2, "11":1},            "S4:581"),
 ("B1s K1 subdivided",            {"01":5, "10":5, "11":1},            "S4:583"),
 ("B3  horn torus",               {"01":2, "10":2, "11":1},            "S4:577"),
 ("B1p K1-bridged",               {"01":3, "10":3},                    "S4:582"),
 ("B1q K1-bridged + SPECTATOR",   {"00":1, "01":3, "10":3},            "S4:582"),
 ("B0b ring torus, loops meet",   {"00":4, "01":1, "10":2, "11":2},    "S4:575"),
 ("B4  spindle",                  {"00":1, "01":1, "10":1, "11":3},    "S4:579"),
]

def chars(occ, f, c):
    """unit-modulus character of each OCCUPIED class. u = conj(W_F) = e^{-if}, v = W_C = e^{ic}."""
    m = {"00": np.ones_like(f, dtype=complex), "10": np.exp(-1j*f),
         "01": np.exp(1j*c),                   "11": np.exp(1j*(c-f))}
    return [m[k] for k in occ]

def zero_in_hull(pts):
    """0 in conv{unit-modulus points} <=> they do not all lie in an open half-plane."""
    A = np.sort(np.angle(np.stack(pts, axis=0)), axis=0)
    g = np.diff(np.concatenate([A, A[:1] + 2*np.pi], axis=0), axis=0)
    return g.max(axis=0) <= np.pi + 1e-12

f = rng.uniform(-np.pi, np.pi, N); c = rng.uniform(-np.pi, np.pi, N)
closed = (np.cos(f) + np.cos(c) <= 0)

print("== W-01's FIRING REGION, ONE VARIABLE MOVING: WHICH CHARACTERS THE INCIDENCE OCCUPIES ==")
print(f"  {'carrier':<30} {'#chars':>6} {'characters':<14} {'firing region':>14} "
      f"{'cos f+cos c<=0?':>16} {'f->-f flips':>12}")
rows = []
for name, mult, src in CARRIERS:
    occ = sorted(mult)
    fire  = zero_in_hull(chars(occ, f, c))
    fireN = zero_in_hull(chars(occ, -f, c))          # ONE variable: sign of the CURVATURE only
    agree = int((fire == closed).sum())
    lbl = " ".join({"00":"1","10":"u","01":"v","11":"uv"}[k] for k in occ)
    print(f"  {name:<30} {len(occ):>6} {lbl:<14} {fire.mean():>14.6f} "
          f"{agree:>10}/{N} {int((fire!=fireN).sum()):>12}")
    rows.append((name, len(occ), fire.mean(), agree, int((fire != fireN).sum())))

print()
print("  READING THE COLUMNS:")
print("  * firing region 1/4 on every carrier whose occupied characters are {u,v,uv} -- K1's set.")
print("  * firing region 1/2 the moment class 00 is occupied AND uv is too: the criterion")
print("    acquires the closed form cos f + cos c <= 0 and agrees with it on every draw.")
print("  * B1p ({u,v}) NEVER fires: two unit points can only hull 0 if antipodal, measure zero.")
print("  * B1q is the sharp case. It HAS a spectator (p00>0) but NO vertex in both loops (p11=0),")
print("    so its three characters are {1,u,v} -- a DIFFERENT three from K1's {u,v,uv}.")
print("  * LAST COLUMN IS W-01's ADVERTISED VIRTUE: 'it distinguishes curvature from flat")
print("    holonomy, which K1 exists to separate.' Sending f -> -f alone (curvature reversed,")
print("    flat holonomy untouched) must change the verdict if the criterion sees curvature.")

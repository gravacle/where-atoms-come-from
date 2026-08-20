"""V4 -- DOES THE FINDING'S OWN TABLE CONTAIN A PAIR THAT CONTRADICTS ITS CONCLUSION?

The finding concludes: "The variation tracks the transport group G/Z(G) and the carrier's
size, never the record count."

TEST 1  MATCHED PAIRS.  Search all 35 carriers for pairs that hold FIXED everything the
        conclusion names as the determinant (|G|, Hilbert dimension, |G/Z|) and differ in the
        variable the conclusion says is inert (the record count).  If such a pair exists AND
        phi_fix differs on it, the conclusion is contradicted by the table above it.
        The mirror search is run in the same table (D-15): pairs holding the record count and
        the dimension fixed while |G/Z| differs.

TEST 2  WHICH INVARIANT IS phi_fix REALLY CONSTANT ON?  I print, for every carrier, the
        NORMALISED transport profile  h -> ( |C_G(h)|/|G| , k(C_G(h))|C_G(h)|/|G|^2 )  as a
        multiset with weights.  All of chi_{E}(h) are built from exactly these two numbers and
        scale as |G|^2, so phi_fix = sum_E <chi_E,chi_E> / sum_E dim_E^2 is a function of this
        normalised profile ALONE.  Carriers sharing the profile are isoclinic (Hall 1940: the
        commutator map G/Z x G/Z -> G' is the isoclinism invariant, and the class equation is
        determined by it up to the |Z| factor).  If the phi-classes coincide exactly with the
        profile-classes, then "phi_fix constant while the record count x16" is the ORDINARY
        statement "phi_fix is an isoclinism invariant and the record count scales with |Z|".

TEST 3  IS phi_fix INTRINSICALLY BLIND TO m?  Counterfactual: hold the excited sectors exactly
        as measured and vary m alone.  If dphi/dm is large and positive, then phi_fix is NOT
        record-count-independent as a function -- it is record-count-independent only because
        m cannot be varied independently on this carrier family.  That is a CONFOUND, not a
        dissociation.
"""
import sys, numpy as np
from collections import defaultdict
LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_SCALE_A_TRANSPORT"
sys.path.insert(0, LANE)
import glib
from carriers import census, isotypic, phi
def say(*a): print(*a); sys.stdout.flush()

R = []
for G in glib.ladder(64):
    ce = census(G); p = phi(isotypic(G, ce), ce['dims'])[2]
    prof = defaultdict(int)
    for h in range(G.n):
        c = int(ce['cent_sz'][h])
        C = G.centralizer(h); seen = set(); kk = 0
        for g in C:
            if g in seen: continue
            orb = {G.conj(x, g) for x in C}
            seen |= orb; kk += 1
        prof[(round(c / G.n, 12), round(kk * c / (G.n ** 2), 12))] += 1
    tot = sum(prof.values())
    nprof = tuple(sorted((k, round(v / tot, 12)) for k, v in prof.items()))
    R.append(dict(name=G.name, n=G.n, dim=G.n**2, m=ce['dims'][-2], GZ=G.n//ce['Z'],
                  ab=G.abelian, phi=p, prof=nprof, dims=[ce['dims'][v] for v in (-2,-1,0)],
                  chis=ce['chis'], nn=G.n))

say("="*128)
say("V4   MATCHED PAIRS, THE INVARIANT phi_fix IS REALLY CONSTANT ON, AND THE COUNTERFACTUAL")
say("="*128)
say("")
say("  TEST 1a  PAIRS WITH |G|, dim AND |G/Z| ALL EQUAL BUT DIFFERENT RECORD COUNT")
say(f"  {'carrier A':<14}{'carrier B':<14}{'|G|':>5}{'dim':>6}{'|G/Z|':>7}{'records A':>11}{'records B':>11}"
    f"{'phi A':>9}{'phi B':>9}{'phi differs?':>14}")
found1 = 0
for i in range(len(R)):
    for j in range(i+1, len(R)):
        a, b = R[i], R[j]
        if a['n']==b['n'] and a['dim']==b['dim'] and a['GZ']==b['GZ'] and a['m']!=b['m']:
            found1 += 1
            say(f"  {a['name']:<14}{b['name']:<14}{a['n']:>5}{a['dim']:>6}{a['GZ']:>7}{a['m']:>11}{b['m']:>11}"
                f"{a['phi']:>9.4f}{b['phi']:>9.4f}{str(abs(a['phi']-b['phi'])>1e-9):>14}")
if not found1: say("    (none)")
say("")
say("  TEST 1b  MIRROR, SAME TABLE: PAIRS WITH dim AND RECORD COUNT EQUAL BUT DIFFERENT |G/Z|")
say(f"  {'carrier A':<14}{'carrier B':<14}{'|G|':>5}{'dim':>6}{'records':>9}{'|G/Z| A':>9}{'|G/Z| B':>9}"
    f"{'phi A':>9}{'phi B':>9}{'phi differs?':>14}")
found2 = 0
for i in range(len(R)):
    for j in range(i+1, len(R)):
        a, b = R[i], R[j]
        if a['dim']==b['dim'] and a['m']==b['m'] and a['GZ']!=b['GZ']:
            found2 += 1
            say(f"  {a['name']:<14}{b['name']:<14}{a['n']:>5}{a['dim']:>6}{a['m']:>9}{a['GZ']:>9}{b['GZ']:>9}"
                f"{a['phi']:>9.4f}{b['phi']:>9.4f}{str(abs(a['phi']-b['phi'])>1e-9):>14}")
if not found2: say("    (none)")
say(f"")
say(f"    pairs found: TEST 1a = {found1}   TEST 1b = {found2}")

say("")
say("  TEST 2  PROFILE CLASSES vs phi CLASSES")
byprof = defaultdict(list); byphi = defaultdict(list)
for r in R:
    byprof[r['prof']].append(r['name']); byphi[round(r['phi'], 10)].append(r['name'])
say(f"  {'phi_fix':>10}   {'record counts in the class':<44}{'#profile classes inside':>26}")
for p in sorted(byphi, reverse=True):
    names = byphi[p]
    ms = ", ".join(str(next(r['m'] for r in R if r['name']==nm)) for nm in names)
    profs = {next(r['prof'] for r in R if r['name']==nm) for nm in names}
    say(f"  {p:>10.6f}   {ms[:44]:<44}{len(profs):>26}")
say("")
say(f"  distinct phi values: {len(byphi)}    distinct normalised transport profiles: {len(byprof)}")
say(f"  profile class == phi class on every carrier: "
    f"{all(len({round(next(r['phi'] for r in R if r['name']==nm),10) for nm in v})==1 for v in byprof.values()) and len(byprof)==len(byphi)}")
say("    -> phi_fix is a function of the NORMALISED transport profile alone; the record count")
say("       inside a profile class is |Z| x (a family constant), which is why it can grow x16")
say("       with phi_fix bit-identical.  That is an isoclinism identity, not a measurement.")

say("")
say("  TEST 3  COUNTERFACTUAL dphi/dm: hold the excited sectors as measured, vary m alone.")
say(f"  {'carrier':<13}{'m measured':>12}{'phi at m':>10}{'phi at m/2':>12}{'phi at 2m':>11}{'phi at 4m':>11}{'phi as m->inf':>15}")
for nm in ("D_4", "D_8", "D_16", "ES_2^(1+4)", "D_32"):
    r = next(x for x in R if x['name'] == nm)
    d = r['dims']; n = r['nn']
    f1 = float(np.sum(r['chis'][-1]**2)/n); f0 = float(np.sum(r['chis'][0]**2)/n)
    def ph(mm): return (mm*mm + f1 + f0)/(mm*mm + d[1]**2 + d[2]**2)
    m = d[0]
    say(f"  {nm:<13}{m:>12}{ph(m):>10.4f}{ph(m/2):>12.4f}{ph(2*m):>11.4f}{ph(4*m):>11.4f}{1.0:>15.4f}")
say("    -> phi_fix rises monotonically toward 1 with the record count when the record count is")
say("       varied INDEPENDENTLY.  It is flat in the lane's data only because m, dim(E(-1)),")
say("       dim(E(0)) all scale together by the same factor on this carrier family.")

"""T-7: does the MODEL reproduce every registered FORMATION result?

validate_model.py covers the existence half -- records constructed from (H,{L_k}) alone.
This covers the formation half: every number below was produced by a one-off lane script, and
is re-derived here THROUGH THE MODEL. Where a lane and the model disagree, one of them is wrong
and we want to know which."""
import sys, numpy as np
sys.path.insert(0, '/Users/bgm/MB Work/where-atoms-come-from/model')
from record_model import RecordModel, Environment
def say(*a): print(*a); sys.stdout.flush()
g = {}
exec(open('/Users/bgm/MB Work/where-atoms-come-from/LANE_F7_OCCUPANCY/f7_davies.py')
     .read().split('say("="*104); say("0.')[0], g)
H0, Zbar, Zbar2 = g['H0'], g['Zbar'], g['Zbar2']
Z, X, op, L, PLAQ, ind = g['Z'], g['X'], g['op'], g['L'], g['PLAQ'], g['ind']
Y = 1j * (X @ Z)
m = RecordModel(H0, []); env = Environment()
P = F = 0
def chk(lbl, got, tgt, tol=1e-7):
    global P, F
    ok = abs(got - tgt) < tol; P += ok; F += (not ok)
    say(f"  [{'PASS' if ok else 'FAIL'}] {lbl:<52} {got:>13.8f}  target {tgt:.8f}")
def chkb(lbl, got, tgt):
    global P, F
    ok = (got == tgt); P += ok; F += (not ok)
    say(f"  [{'PASS' if ok else 'FAIL'}] {lbl:<52} {str(got):>13}  target {tgt}")

say("=" * 96); say("VALIDATE THE FORMATION HALF -- every lane result, through the model"); say("=" * 96)

say("\nF-20  the history: chi grows from a product state while the record's value is fixed")
chk("chi at t = 1.0, coupling Zbar", m.formation(Zbar, Zbar, env, lam=0.8, t=1.0), 0.97527192)
chk("chi at t = 0   (product state, must be exactly 0)", m.formation(Zbar, Zbar, env, lam=0.8, t=0.0), 0.0)

say("\nF-21  redundancy: each bath fragment independently holds the record")
for j, tgt in ((0, 0.789366), (1, 0.048377), (2, 0.678602)):
    chk(f"fragment {{{j}}} at t = 4.0", m.formation(Zbar, Zbar, env, lam=0.8, t=4.0, fragment=[j]), tgt, tol=1e-5)
chk("weight-1 coupling, whole bath (must be 0)",
    m.formation(Zbar, op({ind[('h',0,0)]: Z}, L), env, lam=0.8, t=4.0), 0.0)

say("\nC-18  gauge invariance AND locality forbid formation")
BP = [op({l: Z for l in p}, L) for p in PLAQ]
chk("gauge-invariant local (plaquettes)",
    m.formation(Zbar, [(BP[i], i) for i in range(len(BP))], env, lam=0.8, t=4.0), 0.0)
chk("NOT gauge-invariant, local (sum Z_l)",
    m.formation(Zbar, [(op({l: Z}, L), l) for l in range(L)], env, lam=0.8, t=4.0), 0.21703158)
chk("gauge-invariant, NON-local (Zbar)", m.formation(Zbar, Zbar, env, lam=0.8, t=4.0), 0.90811968)

say("\nG-16  the channel criterion decides formation, and its OLD form does not")
chkb("channel(Zbar, Zbar) opens", m.channel(Zbar, Zbar)['opens_channel'], True)
chkb("channel(Zbar, Zbar2) does NOT", m.channel(Zbar, Zbar2)['opens_channel'], False)
chkb("channel(Zbar, Zbar*Zbar2) does NOT -- the old form said it did",
     m.channel(Zbar, Zbar @ Zbar2)['opens_channel'], False)
chk("  and chi for that coupling really is 0", m.formation(Zbar, Zbar @ Zbar2, env, lam=0.8, t=4.0), 0.0)

say("\nC-17  clause (ii) is a condition on the ENVIRONMENT: generic noise admits no record")
for lbl, Ls, want in (("single-site Z noise: records possible", [op({l: Z}, L) for l in range(L)], True),
                      ("single-site X,Y,Z noise: scalars only",
                       [op({l: Pm}, L) for l in range(L) for Pm in (X, Y, Z)], False)):
    mm = RecordModel(H0, Ls)
    chkb(lbl, len(mm.projs) > 1, want)

say("\nC-19  records compose: forming one leaves the other untouched")
res = m.formation_independence([Zbar, Zbar2], [Zbar, Zbar2], env, lam=0.8, t=4.0)
chkb("both couplings form independently", all(r['independent'] for r in res), True)
chk("  bath learns nothing about the other record",
    max(r['learned'][1 - r['targets'][0]] for r in res), 0.0)

say("\n" + "=" * 96)
say(f"  {P} PASS, {F} FAIL")
say("  Every number above was first produced by a one-off lane script and is re-derived here")
say("  through the model. The model is now the process model, not a parallel implementation.")
sys.exit(1 if F else 0)

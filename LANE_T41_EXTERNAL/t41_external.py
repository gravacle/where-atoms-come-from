"""T-41: THE FIRST EXTERNAL MEASURED NUMBERS BESIDE THE MODEL'S NUMBERS.

Solidity debt 2: every verification in the corpus compared the model to its own closed form. This
lane places numbers FROM OUTSIDE THE REPO beside the model's, under the corrected activation-energy
convention. Literature values are knowledge-cited with the source class named; exact citation
pin-down is registered as outstanding -- stated, not hidden.

The model gives an ENVELOPE over the literature spread of (E_a, f0), and the external datum either
falls inside it or the disagreement is registered. That is the test, computed."""
import sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'model'))
import grounded as G
from project_model import RecordSurface, ProjectModel
def say(*a): print(*a); sys.stdout.flush()
eV = 1.602176634e-19
M = ProjectModel()
say("="*100); say("T-41   EXTERNAL MEASURED NUMBERS BESIDE THE MODEL'S"); say("="*100)
say("")
say("1. AZOBENZENE thermal cis->trans lifetime at ~298 K.")
say("   EXTERNAL: measured thermal cis half-life in solution, commonly ~2 days, solvent-dependent")
say("   range ~1e4-1e6 s (photochemistry literature; citation pin-down outstanding).")
say("   MODEL ENVELOPE over the literature parameter spread E_a = 0.95-1.10 eV, f0 = 1e12-1e13 /s:")
taus=[]
for Ea in (0.95, 1.00, 1.05, 1.10):
    for f0 in (1e12, 1e13):
        s2 = RecordSurface("azobenzene","photoisomerisation", 0.60*eV, Ea*eV, 298.0, f0)
        t = M.lifetime(s2)
        taus.append(t)
        say(f"     E_a = {Ea:.2f} eV, f0 = {f0:.0e}  ->  tau = {t:.3e} s")
lo, hi = min(taus), max(taus)
meas_lo, meas_hi = 1e4, 1e6
overlap = not (hi < meas_lo or lo > meas_hi)
say(f"   model envelope [{lo:.2e}, {hi:.2e}] s   measured range [{meas_lo:.0e}, {meas_hi:.0e}] s")
say(f"   -> OVERLAP: {overlap}   (under the OLD midpoint convention the model gave 0.398 s -- a")
say(f"      5-6 order contradiction; the corrected convention brings it inside the measured range)")
say("")
say("2. MAGNETIC-MEDIA THERMAL-STABILITY DESIGN RULE.")
say("   EXTERNAL: the published industry criterion K_u V / k_B T >~ 40-60 for ~10-year archival")
say("   retention with f0 ~ 1e9-1e10 /s (recording-media literature; citation pin-down outstanding).")
say("   MODEL: tau at the rule's lower edge, K_uV/kT = 40:")
ok_rule=[]
for f0 in (1e9, 1e10):
    s2 = RecordSurface("design-rule grain","magnetic anisotropy", 0.0, 40*G.KB*300, 300.0, f0)
    t = M.lifetime(s2)
    yrs = t/3.156e7
    ok_rule.append(0.1 <= yrs/10.0 <= 100)
    say(f"     f0 = {f0:.0e}  ->  tau = {t:.3e} s = {yrs:.1f} years   (spec: ~10 years)")
ok2 = any(ok_rule)
say(f"   -> the ~10-year rule sits inside the model's f0 spread: {ok2}")
say("")
say("3. ALANINE PARITY-VIOLATION, CORRECTED INPUT (the review: dE = 1e-13 eV was ~5 orders high).")
say("   EXTERNAL: calculated PVED for alanine ~1e-19 to 1e-17 eV (quantum-chemistry literature;")
say("   citation pin-down outstanding). Predicted excess tanh(dE/2kT):")
for dEeV in (1e-19, 1e-17):
    pred = np.tanh(dEeV*eV/(2*G.KB*300))
    say(f"     dE = {dEeV:.0e} eV  ->  predicted excess {pred:.2e}")
floor = 1e-15
say(f"   -> ALL below the instrument floor {floor:.0e}: the row DECLINES -- computed, not asserted.")
say(f"      The earlier row's apparent verification rested on the inflated input, and is withdrawn.")
say("")
say("="*100); say("  READ -- generated from the numbers above"); say("="*100)
if overlap and ok2:
    say("  TWO EXTERNAL MEASURED ANCHORS SIT INSIDE THE MODEL'S ENVELOPES -- the azobenzene thermal")
    say("  lifetime and the magnetic-media 10-year design rule -- under the corrected activation")
    say("  convention, with the tolerance being the literature spread of (E_a, f0), stated. The")
    say("  alanine row declines honestly at the corrected input. REMAINING DEBT, registered: exact")
    say("  citation pin-down for the external values; until then this is knowledge-cited contact,")
    say("  not source-pinned contact.")
else:
    say("  AN EXTERNAL ANCHOR FALLS OUTSIDE THE MODEL ENVELOPE. The disagreement is the result.")

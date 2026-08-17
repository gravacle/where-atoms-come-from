"""
r5 -- REPRODUCTION OF THE TARGET LANE, BEFORE ANY ATTACK IS BELIEVED.
Verify the target's seals, re-run all five of its scripts in a clean process, and diff the output
against its sealed OUT_*.txt.  A refutation of a lane whose arithmetic I could not reproduce would
be worthless, and so would a refutation that quietly used different numbers.
"""
import hashlib, os, subprocess, sys

TGT = "/Users/bgm/MB Work/where-atoms-come-from/LANE_W19_B_CRITERION"
SCRIPTS = ["b0_selftest", "b1_discrimination", "b2_sbs_vs_rdelta",
           "b3_gauge_choices", "b4_threshold"]

print("=" * 104)
print("r5  REPRODUCTION OF LANE W19-B")
print("=" * 104)
print(f"target directory: {TGT}\n")

print("SEALS (target's own SEALS.sha256, verified file by file):")
seals = {}
for line in open(os.path.join(TGT, "SEALS.sha256")):
    h, f = line.split()
    seals[f] = h
bad = 0
for f, h in sorted(seals.items()):
    got = hashlib.sha256(open(os.path.join(TGT, f), "rb").read()).hexdigest()
    ok = (got == h)
    bad += (not ok)
    print(f"   [{'ok ' if ok else 'FAIL'}] {f:<28} {got[:16]}...")
print(f"   {len(seals)} files, {bad} mismatches\n")

print("RE-RUN AND DIFF (each script executed in a fresh python3 process):")
allok = True
for s in SCRIPTS:
    r = subprocess.run([sys.executable, f"{s}.py"], cwd=TGT, capture_output=True, text=True)
    sealed = open(os.path.join(TGT, f"OUT_{s}.txt")).read()
    same = (r.stdout == sealed)
    allok &= same and r.returncode == 0
    print(f"   [{'ok ' if same else 'DIFF'}] {s:<22} exit={r.returncode}  "
          f"sealed {len(sealed)} chars, rerun {len(r.stdout)} chars  "
          f"{'BYTE-IDENTICAL' if same else 'DIFFERS'}")
print(f"\nREPRODUCTION {'COMPLETE -- all five scripts byte-identical' if allok else 'FAILED'}")
print("""
Therefore every number this refutation lane disputes is a number I can regenerate, and every
number it AGREES with is agreed on the target's own output rather than on a paraphrase of it.""")

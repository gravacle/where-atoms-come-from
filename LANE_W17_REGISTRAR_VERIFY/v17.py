# Registrar's verification of W-17's one certified retrodiction, and of its custody claim.
import numpy as np, subprocess, pathlib
print("== V1  THE FOUNDING ESCAPE INFERENCE, TESTED ==")
print("  FOUNDING_DESIGN_V001.md:63-65 : \"An inductive limit of finite objects is not finite --")
print("  which is precisely how it escapes recurrence.\"  That sentence licensed stage S3.")
print("  It infers NON-RECURRENCE from NON-FINITENESS. Non-finiteness is NECESSARY, not SUFFICIENT.")
print()
print("  Counterexample needing no complex, no gauge group, no fibre -- just a commensurate spectrum:")
print(f"  {'dim':>10}  {'spectrum':>10}  {'|A(t = 2 pi)|':>20}")
for d in (10, 100, 10_000, 1_000_000):
    E=np.arange(1,d+1,dtype=float)                       # E_n = n, commensurate
    A=np.exp(-1j*E*2*np.pi).mean()
    print(f"  {d:>10}  {'E_n = n':>10}  {abs(A):>20.15f}")
print("  -> perfect revival at t = 2 pi at EVERY dimension including the limit. An infinite")
print("     commensurate spectrum recurs exactly. THE ESCAPE INFERENCE IS INVALID.")
print()
print("  and the contrast, so the test could have failed -- an INCOMMENSURATE spectrum:")
rng=np.random.default_rng(20260824)
for d in (10, 100, 10_000):
    E=np.sort(rng.uniform(0,1,d))
    t=np.linspace(1,1e4,20001)
    A=np.abs(np.exp(-1j*np.outer(t,E)).mean(axis=1))
    print(f"  {d:>10}  {'random':>10}  max over t in [1,1e4] = {A.max():.6f}")
print()
print("== V2  WAS THE ESCAPE EVER TESTED IN THIS PROGRAM? ==")
reg=pathlib.Path("../REGISTER_V001.md").read_text()
print(f"  REGISTER_V001.md is {len(reg.splitlines())} lines, W-01 through W-16.")
for term in ("inductive limit","quasi-local","UHF","Glimm"):
    print(f"    occurrences of {term!r:<18} in the register: {reg.lower().count(term.lower())}")
print("  -> the escape route that licensed the program's central stage was never examined by any")
print("     of sixteen register rows, eleven of them adversarial.")

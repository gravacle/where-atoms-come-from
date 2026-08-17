# W-15 — measure the registrar's own analytic shape against the predecessor's, on the same file.
# Raised by the principal, 2026-08-17, from a prose tic: "You're locked to an x/y way of thinking
# about almost everything." Tested rather than conceded. Split at W-07, the Fable 5 / Opus 5 boundary.
import re, pathlib
t=pathlib.Path("../REGISTER_V001.md").read_text() if pathlib.Path("../REGISTER_V001.md").exists() \
  else pathlib.Path("REGISTER_V001.md").read_text()
i=t.index("## W-07"); early,mine=t[:i],t[i:]
we,wm=len(early.split()),len(mine.split())
PATS={"'rather than'":r"rather than",
      "'X, not Y'":r"\bnot\b[^.\n]{0,40}\bbut\b|\bnever a\b|\b, not \b",
      "'versus' / 'vs'":r"\bversus\b|\bvs\b",
      "biconditional ('iff')":r"\biff\b|if and only if|BICONDITIONAL",
      "'either ... or'":r"\beither\b[^.\n]{0,60}\bor\b",
      "'two ways' / 'two readings'":r"two ways|two readings|READING A|READING B"}
print(f"  {'construction':<30}{'Fable 5 /1k':>13}{'Opus 5 /1k':>13}{'ratio':>9}")
for n,p in PATS.items():
    a=len(re.findall(p,early,re.I)); b=len(re.findall(p,mine,re.I))
    ra,rb=1000*a/we,1000*b/wm
    print(f"  {n:<30}{ra:>13.2f}{rb:>13.2f}{(rb/ra if ra else float('inf')):>9.2f}")
print(f"\n  words: Fable 5 {we}, Opus 5 {wm}")
print("\n  READING: the PHRASE-level tic is inherited house style -- 'rather than' 1.50 vs 1.53,")
print("  'X, not Y' 3.60 vs 2.51 (the registrar is LOWER). What is distinctively the registrar's")
print("  is the STRUCTURAL habit: biconditionals at 2.5x and two-way readings at ~6x.")
print("  The defect is not talking in x/y. It is CONVERTING questions into x/y, because a binary")
print("  is decidable by a computation and an open question is not.")

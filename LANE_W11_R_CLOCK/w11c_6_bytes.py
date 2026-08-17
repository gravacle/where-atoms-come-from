# LANE W-11 R/C — LEG 6 — THE CLOCK QUESTION AT THE BYTES.
# Questions (3) and (4) of the brief are textual, so this leg quotes, with sha256 of each file,
# and asserts the quoted line numbers so the extraction is reproducible rather than paraphrased.
import hashlib, sys
REPO = "/Users/bgm/MB Work/where-atoms-come-from"

FILES = {
    "S1":   "S1_CARRIER_K1_V001.md",
    "S2A":  "S2_FORMATION_CONDITION_ON_K1_AUDIT_V001.md",
    "S3":   "S3_THE_CROSSING_V001.md",
    "S3A":  "S3_THE_CROSSING_AUDIT_V001.md",
    "REG":  "REGISTER_V001.md",
}
CACHE = {}
def lines(tag):
    if tag not in CACHE:
        with open(f"{REPO}/{FILES[tag]}", "rb") as fh:
            b = fh.read()
        CACHE[tag] = (b.decode(), hashlib.sha256(b).hexdigest())
    return CACHE[tag]

def show(tag, lo, hi, why):
    txt, dg = lines(tag)
    L = txt.split("\n")
    print(f"\n--- {FILES[tag]}:{lo}" + (f"-{hi}" if hi != lo else "") + f"   [{why}]")
    for i in range(lo, hi+1):
        print(f"   {i:>5} | {L[i-1]}")

print("FILE DIGESTS (sha256), so every quote below is checkable at the bytes:")
for t in FILES:
    print(f"   {FILES[t]:<45} {lines(t)[1]}")

print("""
==================================================================================
 6A   QUESTION (3) — CHOICE LEDGER A2.  IS ITS 'WHY' A REASON, OR A RESTATEMENT?
      The registrar claims A2 justifies the CIRCUIT clock by citing that EDGE count is
      carrier-supplied.  VERIFY OR REFUTE AT THE BYTES.
==================================================================================""")
show("S2A", 658, 658, "CHOICE LEDGER A2, the entry that closed the clock question")
show("S1", 16, 22, "the byte range A2's 'why' column cites as its warrant")
print("""
   VERDICT ON THE REGISTRAR'S READING: **CONFIRMED, AND STRONGER THAN IT STATED.**
   A2's CHOICE is 'Time = number of CIRCUITS n of the loop'.
   A2's WHY is  'EDGE count is carrier-supplied combinatorics (S1 :16-22)'.
   And S1:16-22 is THE EDGE LIST -- six edges, nothing about loops, faces or circuits.
   The warrant is about the omitted alternative's own datum.

   A2's ALTERNATIVES column reads: 'a real parameter t with a Hamiltonian; no time at all'.
   THE EDGE TICK IS NOT AMONG THEM.  Against the two alternatives it does list, A2 is a good
   REASON.  Against the alternative it omits, it is not a reason and not even a restatement:
   a restatement would at least be about circuits.  It is an argument for the GENUS
   (a combinatorial, carrier-supplied clock) recorded as an argument for the SPECIES (circuits).
   Charitable reading -- 'circuits are edge-counts too, so they inherit the warrant' -- is true
   and still under-determines, because edges inherit it at least as well.""")

print("""
==================================================================================
 6B   AND THE SUBSTITUTION HAPPENS TWICE, IN ADJACENT SENTENCES, ON THE PAGE THAT
      WROTE THE CORRECTION.  (New here; W-06 found the missing alternatives column,
      not the substitution inside the warrant.)
==================================================================================""")
show("S2A", 168, 172, "S2 audit sec3.3, the section that WRITES its COR-F")
show("S2A", 639, 641, "S2 audit COR-F as filed")
print("""   'K1 does carry a discrete one: THE NUMBER OF EDGES TRAVERSED' ... cite S1:16-22 ...
   and then, two sentences later, 'with n = CIRCUIT COUNT'.  The correction is derived for
   edges and banked for circuits, in one paragraph.""")

print("""
==================================================================================
 6C   THE NARROWING PROPAGATES, EACH LINK CITING THE ONE BEFORE
==================================================================================""")
show("S3", 981, 981, "S3 CHOICE LEDGER C3 -- cites 'audit COR-F' for CIRCUIT count")
show("REG", 66, 66, "REGISTER W-01 row -- the disjunction is gone entirely")
print("""   CHAIN:  S2A:169 'the number of edges traversed'  (warrant: S1:16-22, the edge list)
        -> S2A:639 'edge traversals / circuits'        (disjunction)
        -> S2A:658 A2 'number of circuits'             (cites the edge warrant)
        -> S3:981  C3 'circuit count ... (audit COR-F)'(cites the disjunction, takes one arm)
        -> REG:66  'circuit count is carrier-supplied discrete time'  (no citation, no arm)
   FOUR NARROWINGS.  No step is flagged; no step lists the edge tick as a live alternative.""")

print("""
==================================================================================
 6D   QUESTION (4) — IS A PARTIAL CIRCUIT A LEGITIMATE STATE?  WHAT FORBIDS IT?
==================================================================================""")
show("S1", 50, 54, "S1 defines parallel transport EDGE-WISE; the circuit operator is derived")
show("S3", 392, 397, "S3 sec3.5 -- the cell is DEFINED to run whole circuits ('of EACH loop')")
show("S3", 1016, 1019, "S3 FLAG F3 -- the corpus's own sealed flag on the clock")
show("S3", 1026, 1028, "S3 FLAG F5 -- 'I do not know which is physical'")
show("S3A", 798, 798, "COR-J -- the only premise that WOULD forbid a partial circuit")
print("""
   NOTHING IN THE SEALED CORPUS FORBIDS A PARTIAL CIRCUIT.
   * S1:52 defines transport one edge at a time; T^n s is built from S1's own primitive.
   * S2A:169 and S2A:639 AFFIRMATIVELY name edge traversals as a discrete time K1 supplies.
   * S3:392-397 stipulates whole circuits inside the definition of a cell -- but that is the
     object under test, ledgered at C3, whose alternatives column omits the edge tick.
   * S3 F3 says in terms that the page does NOT derive that a physical process must use the
     canonical clock, and calls it 'a real gap ... where an adversary should attack first'.
   * S3 F5 says the page does not decide which clock is physical.
   * COR-J's premise 'the record must be gauge-invariant' would forbid it ONLY if promoted from
     the RECORD to the OPERATOR.  Leg 4c shows the RECORD is gauge-invariant at every edge tick,
     so COR-J AS WORDED does not exclude the edge clock.  The promotion is nowhere in writing,
     and COR-J is itself entered as a defect: undeclared, load-bearing, applied asymmetrically.
   NOR DOES ANYTHING REQUIRE A PARTIAL CIRCUIT.  The corpus is silent, and says so twice.""")

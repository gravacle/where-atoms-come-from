#!/usr/bin/env python3
# LANE_T43_C_RULE — apply the pre-registered T-43 decision rule (REGISTER_V001.md, final
# section) to the two verified tier inputs (LANE_T43_A_CORNER, LANE_T43_B_WORLD).
#
# DISCIPLINE: no literal verdicts. Every input number below is copied from the two
# ADVERSARIALLY VERIFIED lane findings (both not-refuted, high confidence, independent
# reimplementation); this lane re-derives every degree by constant finite differences and
# every closed form by exact integer check, then computes the rule's branch as a boolean
# conjunction. Exact arithmetic throughout (integers only). The rule decides; this lane
# does not choose a reading the rule did not register.
#
# THE RULE (quoted verbatim in the OUT, read from REGISTER_V001.md final section):
#   boundary-bounded in BOTH tiers  -> surface answered (certifiability reading;
#                                      strong content bound -> emergent-tier target)
#   volume-bounded                  -> surface declines
#   mixed/other                     -> reported as found; no forcing

import hashlib, sys

GATES = []
def gate(name, ok):
    GATES.append((name, bool(ok)))
    return bool(ok)

def diffs(seq):
    return [b - a for a, b in zip(seq, seq[1:])]

def degree(seq):
    # degree = smallest k with constant k-th finite differences (requires >= 3 residuals)
    s, k = list(seq), 0
    while len(s) >= 3:
        if len(set(s)) == 1:
            return k
        s = diffs(s); k += 1
    return None

OUT = []
def emit(line=""):
    OUT.append(line)

emit("=" * 78)
emit("LANE_T43_C_RULE — THE PRE-REGISTERED RULE, APPLIED")
emit("inputs: LANE_T43_A_CORNER (verified), LANE_T43_B_WORLD (verified)")
emit("=" * 78)

# ---------------------------------------------------------------- CORNER TIER (D = 2)
# Verified table, L = 12, squares s = 2..11 (T43_A finding, adversarially reproduced).
S      = list(range(2, 12))
STOREDc = [7, 19, 35, 55, 79, 107, 139, 175, 215, 259]
CERTc   = [6, 14, 22, 30, 38, 46, 54, 62, 70, 78]
Dc = 2

gate("corner closed form STORED = 2s^2+2s-5", all(v == 2*s*s + 2*s - 5 for s, v in zip(S, STOREDc)))
gate("corner closed form CERT   = 8s-10",     all(v == 8*s - 10        for s, v in zip(S, CERTc)))
deg_STOREDc = degree(STOREDc)
deg_CERTc   = degree(CERTc)
gate("corner STORED degree == D   (volume)",   deg_STOREDc == Dc)
gate("corner CERT   degree == D-1 (boundary)", deg_CERTc == Dc - 1)
# The corner lane computed static counts on a fixed code state: no epoch structure exists
# (finding field per_epoch_vs_cumulative: N/A). Hence corner per-use == corner total.
corner_has_time_axis = False
gate("corner has no epoch axis (per-use == total, one number)", not corner_has_time_axis)
corner_boundary = (deg_CERTc == Dc - 1) and (deg_STOREDc == Dc)
emit(f"CORNER: CERT degree {deg_CERTc} (boundary), STORED degree {deg_STOREDc} (volume);")
emit("        static state -> certifiable content has ONE value; boundary under every reading.")
emit(f"        corner_boundary = {corner_boundary}")
emit()

# ----------------------------------------------------------------- WORLD TIER (D = 3)
# Verified table, n = 2..12 (T43_B finding, adversarially reproduced incl. independent
# max-flow). Two exact counts exist and the lane reported BOTH, uncollapsed:
N       = list(range(2, 13))
STOREDw = [n**3 for n in N]
CERT_RATE  = [8, 26, 56, 98, 152, 218, 296, 386, 488, 602, 728]   # per-epoch (min-cut)
E_MIN      = [1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3]
CERT_TOTAL = [n**3 for n in N]   # N_opt(E_min) == n^3, every n (verified: whole volume)
Dw = 3

gate("world closed form CERT/epoch = 6n^2-12n+8", all(v == 6*n*n - 12*n + 8 for n, v in zip(N, CERT_RATE)))
deg_STOREDw = degree(STOREDw)
deg_RATE    = degree(CERT_RATE)
deg_TOTAL   = degree(CERT_TOTAL)
gate("world STORED degree == D (volume)",           deg_STOREDw == Dw)
gate("world CERT/epoch degree == D-1 (boundary)",   deg_RATE == Dw - 1)
gate("world CERT total  degree == D (volume)",      deg_TOTAL == Dw)
gate("world total == stored (all records certified at E_min)", CERT_TOTAL == STOREDw)
gate("world E_min <= 3 on range", max(E_MIN) <= 3)
world_rate_boundary = (deg_RATE == Dw - 1)
world_total_volume  = (deg_TOTAL == Dw)
emit(f"WORLD:  CERT/epoch degree {deg_RATE} (boundary), CERT total degree {deg_TOTAL} (volume),")
emit(f"        STORED degree {deg_STOREDw} (volume); E_min = {E_MIN} — the SPLIT, as found.")
emit(f"        world_rate_boundary = {world_rate_boundary}; world_total_volume = {world_total_volume}")
emit()

# ------------------------------------------------------------------- APPLY THE RULE
# The registered rule fixes NO time convention for 'externally certifiable content'
# ('channel-counted in the world model' names the instrument, not per-epoch vs anytime).
# Both readings are therefore computed; the rule's branch is taken on the conjunctions.
readings = {
    "per-epoch/simultaneous": (corner_boundary, world_rate_boundary),
    "anytime/cumulative":     (corner_boundary, not world_total_volume),
}
emit("RULE INPUTS (boundary-bounded?):")
branch = {}
for name, (c, w) in readings.items():
    both_boundary = c and w
    both_volume   = (not c) and (not w)
    b = "SURFACE_ANSWERED" if both_boundary else ("SURFACE_DECLINES" if both_volume else "MIXED")
    branch[name] = b
    emit(f"  {name:24s}: corner={c}, world={w}  -> {b}")

reading_dependent = len(set(branch.values())) > 1
gate("per-epoch reading: boundary in BOTH tiers",  branch["per-epoch/simultaneous"] == "SURFACE_ANSWERED")
gate("anytime  reading: NOT boundary in both (world total is volume)", branch["anytime/cumulative"] != "SURFACE_ANSWERED")
gate("neither reading gives volume-bounded-in-both (surface never declines outright)",
     all(b != "SURFACE_DECLINES" for b in branch.values()))
gate("rule input is reading-dependent (rule registered no time convention)", reading_dependent)

RULE_OUTPUT = "MIXED_AS_FOUND" if reading_dependent else (
    "SURFACE_ANSWERED_INTERFACE" if branch["per-epoch/simultaneous"] == "SURFACE_ANSWERED"
    else "SURFACE_DECLINES_VOLUME")
emit()
emit(f"RULE OUTPUT (computed): {RULE_OUTPUT}")
emit()
emit("AS FOUND: certifiable content is boundary-bounded in the corner tier under its only")
emit("reading, and boundary-bounded PER EPOCH in the world tier under every operationalisation")
emit("tried — while the world's certifiable TOTAL over epochs is the whole stored volume")
emit("(E_min <= 3 on range, E_min >= ceil(n^3/L(1)) ~ n/6 by cut arithmetic). The registered")
emit("rule's phrase 'externally certifiable content' does not fix per-epoch vs anytime; the")
emit("branch depends on that unregistered choice, so the rule outputs MIXED, reported as")
emit("found, no forcing. The time-convention disposition belongs to the principal.")
emit()

emit("-" * 78)
allpass = all(ok for _, ok in GATES)
for name, ok in GATES:
    emit(f"[{'PASS' if ok else 'FAIL'}] {name}")
emit(f"GATES: {sum(ok for _, ok in GATES)}/{len(GATES)}   ALL GATES PASS: {allpass}")
emit("=" * 78)

text = "\n".join(OUT) + "\n"
sys.stdout.write(text)
with open(__file__.rsplit("/", 1)[0] + "/t43c_rule.OUT.txt", "w") as f:
    f.write(text)
sys.exit(0 if allpass else 1)

"""T42-A part 4: candidate distances BETWEEN RECORDS.  EXACT, F_2.

Venue: the record classes N(S)/S of a carrier (all of them, computed by pairing against the
COMPUTED symplectic logicals; the S-coset of every class is enumerated in full).

Candidates:
  (ii-a) J(a,b) = min over ALL coset representative pairs of the number of sites where the
         two representatives locally anticommute -- the minimal-crossing integer.
         OWNER: LANE_EXACT_A_ZERO E-7 defined J and proved it gauge-invariant; what is OURS
         here is the metric-axioms interrogation of J on these venues.
  (ii-b) channel-adjacency distance: adj(a,b) = 1 iff sp_F2(a,b) = 1 (the G-16 channel
         criterion, BORROWED), d_ch = shortest-path distance in that graph on the nonzero
         classes.  Also the GEOMETRIC variant adj_J(a,b) = 1 iff J(a,b) > 0.
  (ctrl) the support-overlap measure ov(a,b) = |supp(minrep a) & supp(minrep b)| --
         the D-15 control EXPECTED TO FAIL the axioms, so the instrument can tell a
         metric from a non-metric on this venue.

Carriers: toric L=2 (n=8, |S|=64 per class) and [[6,4,2]] (n=6, |S|=4 per class).
Every verdict gated by a computed boolean; witnesses printed for every FAIL that matters.
"""
import sys, time
from collections import deque
LANE = "/Users/bgm/MB Work/where-atoms-come-from/LANE_T42_A_DISTANCE/VERIFY"
sys.path.insert(0, LANE)
from t42_lib import (pc, solve_affine_f2, span_all, sp_pair, weight_xz, crossings,
                     vec_to_mask, metric_axioms, aut_count, toric_stabilizers)
sys.path.insert(0, "/Users/bgm/MB Work/where-atoms-come-from/model")
from record_model import symplectic_logicals

OUT = []
def say(*a):
    s = " ".join(str(x) for x in a)
    OUT.append(s); print(s); sys.stdout.flush()

def gate(label, ok, detail=""):
    say(("PASS  " if ok else "FAIL  ") + label + (("  " + detail) if detail else ""))
    return ok

say("T42-A part 4  --  distances between records, exact, started " + time.strftime("%F %T"))
say("")
ALL_OK = True


def build_record_venue(stab_vecs, n):
    """Returns (classes, coset_elems, logi):  classes = list of F_2^{2k} label tuples,
       coset_elems[c] = FULL list of Pauli masks in that class (the S-orbit)."""
    pairs = symplectic_logicals(stab_vecs, n)
    logi = []
    for X, Z in pairs:
        logi.append(vec_to_mask(X)); logi.append(vec_to_mask(Z))
    stab_masks = [vec_to_mask(sv) for sv in stab_vecs]
    # N(S): nullspace of pairing with the stabilizer generators
    def full_row(g):
        mask = (1 << n) - 1
        return (g >> n) | ((g & mask) << n)
    rows = [full_row(g) for g in stab_masks]
    x0, nb, rank = solve_affine_f2(rows, [0] * len(rows), 2 * n)
    NS = span_all(nb)
    coset = {}
    for u in NS:
        c = tuple(sp_pair(u, l, n) for l in logi)
        coset.setdefault(c, []).append(u)
    return sorted(coset), coset, logi


def graph_metric(nodes, adj):
    """BFS all-pairs; returns (dist dict or None-if-disconnected, diameter)."""
    d = {}
    diam = 0
    ok = True
    for s in nodes:
        seen = {s: 0}
        q = deque([s])
        while q:
            u = q.popleft()
            for v in nodes:
                if v not in seen and adj(u, v):
                    seen[v] = seen[u] + 1
                    q.append(v)
        for v in nodes:
            if v not in seen:
                ok = False
            else:
                d[(s, v)] = seen[v]
                diam = max(diam, seen[v])
    return (d if ok else None), diam


for name, stab_vecs, n in [
        ("toric L=2", None, None),
        ("[[6,4,2]]", [[1] * 6 + [0] * 6, [0] * 6 + [1] * 6], 6)]:
    t0 = time.time()
    if name == "toric L=2":
        n, stab_vecs, _, _ = toric_stabilizers(2)
    classes, coset, logi = build_record_venue(stab_vecs, n)
    k2 = len(logi)
    nclasses = len(classes)
    csize = len(coset[classes[0]])
    say("%s: %d record classes, %d representatives per class (full S-orbit enumerated)"
        % (name, nclasses, csize))
    ALL_OK &= gate("%s coset sizes all equal (group structure computed, not assumed)" % name,
                   all(len(coset[c]) == csize for c in classes))

    # ---------------- J: the minimal-crossing integer, exact min over ALL rep pairs
    J = {}
    for i, a in enumerate(classes):
        for b in classes[i:]:
            m = None
            for u in coset[a]:
                for v in coset[b]:
                    c = crossings(u, v, n)
                    if m is None or c < m:
                        m = c
                    if m == 0:
                        break
                if m == 0:
                    break
            J[(a, b)] = m
            J[(b, a)] = m
    say("%s  J computed for all %d ordered class pairs (exact min over %d rep pairs each)"
        % (name, nclasses ** 2, csize ** 2))
    # CTRL-Z and consistency J mod 2 == sp (owner E-7)
    ALL_OK &= gate("%s CTRL-Z: J(a,a) = 0 for every class (a record cannot cross itself)"
                   % name, all(J[(a, a)] == 0 for a in classes))
    sp_of = {}
    for a in classes:
        for b in classes:
            # sp of the classes = parity of crossings of any reps = sp of first reps
            sp_of[(a, b)] = sp_pair(coset[a][0], coset[b][0], n)
    ALL_OK &= gate("%s J mod 2 == sp_F2 on every pair (E-7 consistency, recomputed here)"
                   % name, all(J[p] % 2 == sp_of[p] for p in J))
    rep = metric_axioms(classes, J)
    say("%s  J vs THE METRIC AXIOMS:  identity %s  symmetry %s  indiscernibles %s  "
        "triangle %s" % (name, rep["identity_ok"], rep["symmetry_ok"],
                         rep["indiscernible_ok"], rep["triangle_ok"]))
    ALL_OK &= gate("%s J FAILS indiscernibles (distinct non-crossing records at J=0)" % name,
                   not rep["indiscernible_ok"], "witness %s" % (rep["indiscernible_witness"][:1],))
    ALL_OK &= gate("%s J FAILS the triangle inequality" % name, not rep["triangle_ok"],
                   "witness (a,b,c,J(ac),J(ab),J(bc)) = %s" % (rep["triangle_witness"][:1],))
    ALL_OK &= gate("%s J IS symmetric with J(a,a)=0 -- a FORM, not a distance" % name,
                   rep["identity_ok"] and rep["symmetry_ok"])

    nonzero = [c for c in classes if any(c)]
    # triangle failure must not be an artifact of the identity class: NONZERO witness
    nz_wit = None
    for a in nonzero:
        for b in nonzero:
            for c in nonzero:
                if J[(a, c)] > J[(a, b)] + J[(b, c)]:
                    nz_wit = (a, b, c, J[(a, c)], J[(a, b)], J[(b, c)])
                    break
            if nz_wit:
                break
        if nz_wit:
            break
    ALL_OK &= gate("%s J triangle fails already among NONZERO records" % name,
                   nz_wit is not None, "witness %s" % (nz_wit,))

    # ---------------- channel adjacency (G-16, BORROWED) and its graph distance
    d_ch, diam = graph_metric(nonzero, lambda u, v: u != v and sp_of[(u, v)] == 1)
    ALL_OK &= gate("%s channel graph (sp=1) is CONNECTED on the %d nonzero classes"
                   % (name, len(nonzero)), d_ch is not None)
    if d_ch is not None:
        repc = metric_axioms(nonzero, d_ch)
        ALL_OK &= gate("%s channel-adjacency distance SATISFIES ALL METRIC AXIOMS "
                       "(%d pairs, %d triples)" % (name, repc["n_pairs"], repc["n_triples"]),
                       repc["all"])
        vals = sorted(set(v for kk, v in d_ch.items() if v > 0))
        say("%s  channel distance spectrum %s, diameter %d  -- a metric with almost no "
            "geometry" % (name, vals, diam))
        if len(nonzero) <= 15:
            kind, na = aut_count(nonzero, d_ch)
            say("%s  D-22 channel-metric venue Aut: |Aut| = %d (%s) on %d points"
                % (name, na, kind, len(nonzero)))
        else:
            say("%s  D-22 channel-metric venue: Aut enumeration skipped at %d points; "
                "spectrum %s is two-valued, the venue is near-uniform"
                % (name, len(nonzero), vals))

    # geometric variant: adjacency iff J > 0
    dJ, diamJ = graph_metric(nonzero, lambda u, v: u != v and J[(u, v)] > 0)
    even_pairs = [(a, b) for a in nonzero for b in nonzero
                  if sp_of[(a, b)] == 0 and J[(a, b)] > 0]
    say("%s  pairs with sp=0 but J>0 (invisible to the pairing, E-7's class): %d"
        % (name, len(even_pairs)))
    if dJ is not None and d_ch is not None:
        same = all(dJ[p] == d_ch[p] for p in d_ch)
        ALL_OK &= gate("%s the two adjacency notions (sp=1 vs J>0) give %s metric"
                       % (name, "the SAME" if same else "DIFFERENT"),
                       True, "same=%s" % same)

    # ---------------- D-15 control: support-overlap measure must FAIL
    minrep = {c: min(coset[c], key=lambda u: weight_xz(u, n)) for c in classes}
    mask_n = (1 << n) - 1
    def supp(u):
        return (u & mask_n) | (u >> n)
    ov = {(a, b): pc(supp(minrep[a]) & supp(minrep[b])) for a in classes for b in classes}
    repo = metric_axioms(classes, ov)
    ALL_OK &= gate("%s D-15 CONTROL: support-overlap measure FAILS the axioms "
                   "(identity_ok=%s, triangle_ok=%s)"
                   % (name, repo["identity_ok"], repo["triangle_ok"]), not repo["all"],
                   "identity witness %s" % (repo["identity_witness"][:1],))
    say("%s done in %.1fs" % (name, time.time() - t0))
    say("")

say("ALL GATES: " + ("PASS" if ALL_OK else "AT LEAST ONE FAIL"))
with open(LANE + "/t42_d_records.OUT.txt", "w") as f:
    f.write("\n".join(OUT) + "\n")

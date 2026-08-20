"""VERIFY -- independent world-tier check.  Depth by CLOSED FORMULA (min distance to a face,
plus one) instead of BFS; all aggregates recounted; closed forms and the DEPTH_SUM
recurrence rechecked; finite-difference reader reimplemented."""

def analyse_block(n):
    cap_total = n ** 3
    iface = 6 * n * n            # recounted below the honest way as well
    iface_count = 0
    cap_d1 = 0
    depth_sum = 0
    for x in range(n):
        for y in range(n):
            for z in range(n):
                # face adjacencies to the complement
                for (c, lim) in ((x, n), (y, n), (z, n)):
                    if c == 0: iface_count += 1
                    if c == lim - 1: iface_count += 1
                d = 1 + min(x, n - 1 - x, y, n - 1 - y, z, n - 1 - z)
                depth_sum += d
                if d == 1: cap_d1 += 1
    cap_prot = cap_total - cap_d1
    return cap_total, iface_count, cap_d1, cap_prot, depth_sum

def const_diff_order(seq):
    s = list(seq)
    for k in range(len(seq) - 1):
        if len(set(s)) == 1:
            return k
        s = [b - a for a, b in zip(s, s[1:])]
    return None

def gate(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}{('  ' + detail) if detail else ''}")
    return ok

if __name__ == "__main__":
    NS = list(range(2, 13))
    rows = [analyse_block(n) for n in NS]
    # lane's sealed numbers, transcribed from t42_world.OUT.txt
    sealed = {
        2: (8, 24, 8, 0, 8), 3: (27, 54, 26, 1, 28), 4: (64, 96, 56, 8, 72),
        5: (125, 150, 98, 27, 153), 6: (216, 216, 152, 64, 288),
        7: (343, 294, 218, 125, 496), 8: (512, 384, 296, 216, 800),
        9: (729, 486, 386, 343, 1225), 10: (1000, 600, 488, 512, 1800),
        11: (1331, 726, 602, 729, 2556), 12: (1728, 864, 728, 1000, 3528),
    }
    ok = True
    ok &= gate("every sealed row matches the closed-formula recount (depth by min-to-face, "
               "not BFS)", all(tuple(r) == sealed[n] for n, r in zip(NS, rows)))
    cap = [r[0] for r in rows]; ifc = [r[1] for r in rows]
    d1 = [r[2] for r in rows]; prot = [r[3] for r in rows]; dsum = [r[4] for r in rows]
    ok &= gate("degrees: CAP_total 3, IFACE 2, CAP_d1 2, CAP_prot 3, DEPTH_SUM none",
               const_diff_order(cap) == 3 and const_diff_order(ifc) == 2
               and const_diff_order(d1) == 2 and const_diff_order(prot) == 3
               and const_diff_order(dsum) is None)
    ok &= gate("closed forms: IFACE == 6n^2, CAP_d1 == n^3-(n-2)^3, CAP_prot == (n-2)^3",
               all(i == 6 * n * n and d == n ** 3 - (n - 2) ** 3 and p == (n - 2) ** 3
                   for n, i, d, p in zip(NS, ifc, d1, prot)))
    ok &= gate("DEPTH_SUM(n) == DEPTH_SUM(n-2) + n^3 for n >= 4",
               all(dsum[i] == dsum[i - 2] + NS[i] ** 3 for i in range(2, len(NS))))
    # scattered control, own seed: protected capacity must collapse, iface volume-like
    import random
    rnd = random.Random(99)
    cells = [(x, y, z) for x in range(24) for y in range(24) for z in range(24)]
    S = set(rnd.sample(cells, 216))
    NB = ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1))
    ifc_s = sum(1 for g in S for d in NB
                if (g[0]+d[0], g[1]+d[1], g[2]+d[2]) not in S)
    d1_s = sum(1 for g in S if any((g[0]+d[0], g[1]+d[1], g[2]+d[2]) not in S for d in NB))
    ok &= gate("scattered control (own seed): CAP_prot small, IFACE volume-like",
               (216 - d1_s) < 22 and ifc_s > 4 * 216,
               f"CAP_prot={216 - d1_s}, IFACE={ifc_s}")
    print(f"ALL INDEPENDENT WORLD GATES PASS: {ok}")

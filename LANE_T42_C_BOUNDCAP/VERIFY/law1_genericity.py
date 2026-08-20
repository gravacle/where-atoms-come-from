"""VERIFY -- ATTACK on the ownership claims around LAW-1 and the 'two independent routes'.

CLAIM UNDER TEST 1 (algebraic): the lane's two cut-rank routes are NOT independent.
  IR2  = r_S - r_in - r_out  with  r_in = (m - rk_st_Rc - 1) + (m - rk_pl_Rc - 1),
                                   r_out = (m - rk_st_R - 1) + (m - rk_pl_R - 1),
                                   r_S = 2(m - 1)
  IR2b = (rk_st_R - (m - rk_st_Rc - 1)) + (rk_pl_R - (m - rk_pl_Rc - 1))
  Substituting: IR2 - IR2b = 2(m-1) - r_out - rk_st_R - rk_pl_R = 0 IDENTICALLY.
  Both are functions of the same four ranks (rk_st_R, rk_st_Rc, rk_pl_R, rk_pl_Rc);
  the gate 'IR2 == IR2b' can never fail, whatever the carrier or the code state.
  Likewise SYN == IR2 + r_in is the same substitution.  (The cap_record direct-vs-formula
  check IS a genuine cross-check; these two are not.)

CLAIM UNDER TEST 2 (numerical): LAW-1  IR2 + 2 r_in + cap_record == 2|R|  is an identity
of ANY CSS stabiliser group and ANY region -- rank-nullity plus the symplectic pairing --
not a law of the record surface.  We test it on RANDOM CSS codes (random H; X-checks =
rows of H, Z-checks = kernel basis of H) with RANDOM regions.  If it never fails there,
it carries no venue information and 'OURS as a discovered law of this venue' is an
overclaim (it is still a correct and useful bookkeeping identity).
"""
import random

def rank_and_basis(vs):
    piv = {}
    for v in vs:
        while v:
            t = v.bit_length() - 1
            if t in piv:
                v ^= piv[t]
            else:
                piv[t] = v
                break
    return len(piv), list(piv.values())

def span_dim_inside(vs, mask):
    d, _ = rank_and_basis(vs)
    dp, _ = rank_and_basis([v & ~mask for v in vs])
    return d - dp

def kernel_vectors(rows, mask, n):
    cols = [c for c in range(n) if (mask >> c) & 1]
    piv = {}
    out = []
    for idx, c in enumerate(cols):
        v = 0
        for i, rw in enumerate(rows):
            if (rw >> c) & 1:
                v |= 1 << i
        tag = 1 << idx
        while v:
            t = v.bit_length() - 1
            if t in piv:
                pv, pt = piv[t]
                v ^= pv
                tag ^= pt
            else:
                piv[t] = (v, tag)
                break
        if v == 0:
            m = 0
            for j, cc in enumerate(cols):
                if (tag >> j) & 1:
                    m |= 1 << cc
            out.append(m)
    return out

def kernel_basis_of_matrix(rows, n):
    """basis of the right kernel of the matrix whose rows are `rows` (bitmask columns)."""
    return kernel_vectors(rows, (1 << n) - 1, n)

def gate(label, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}{('  ' + detail) if detail else ''}")
    return ok

if __name__ == "__main__":
    rnd = random.Random(2026)
    trials = 0
    law1_holds = True
    routes_agree = True
    for t in range(200):
        n = rnd.randrange(6, 16)
        kx = rnd.randrange(1, n - 2)
        Xchecks = []
        for _ in range(kx):
            v = rnd.getrandbits(n)
            if v: Xchecks.append(v)
        if not Xchecks: continue
        # Z checks: any subset of the kernel of the X-check matrix (guarantees CSS commutation)
        kerb = kernel_basis_of_matrix(Xchecks, n)
        if not kerb: continue
        Zchecks = [b for b in kerb if rnd.random() < 0.7] or kerb[:1]
        for _ in range(5):
            R = rnd.getrandbits(n)
            Rc = ((1 << n) - 1) ^ R
            AQ = bin(R).count("1")
            r_in = span_dim_inside(Xchecks, R) + span_dim_inside(Zchecks, R)
            r_out = span_dim_inside(Xchecks, Rc) + span_dim_inside(Zchecks, Rc)
            r_S = rank_and_basis(Xchecks)[0] + rank_and_basis(Zchecks)[0]
            IR2 = r_S - r_in - r_out
            SYN = rank_and_basis([v & R for v in Xchecks])[0] \
                + rank_and_basis([v & R for v in Zchecks])[0]
            IR2b = SYN - r_in
            # cap_record exactly as the lane defines it (direct span-membership route)
            kerX = kernel_vectors(Zchecks, R, n)   # X-ops in R commuting with all Z-checks
            kerZ = kernel_vectors(Xchecks, R, n)
            _, bX = rank_and_basis(Xchecks)
            _, bZ = rank_and_basis(Zchecks)
            capX = rank_and_basis(bX + kerX)[0] - len(bX)
            capZ = rank_and_basis(bZ + kerZ)[0] - len(bZ)
            cap_record = capX + capZ
            trials += 1
            if IR2 + 2 * r_in + cap_record != 2 * AQ:
                law1_holds = False
            if IR2 != IR2b:
                routes_agree = False
    ok = True
    ok &= gate(f"LAW-1 holds on EVERY random CSS code + random region tried "
               f"({trials} trials, codes far from toric): it is a GENERIC stabiliser "
               f"identity, not a venue law", law1_holds and trials > 500,
               f"trials={trials}")
    ok &= gate("IR2 == IR2b on every trial as well -- consistent with the algebraic "
               "reduction in the header: the 'two independent routes' share the same four "
               "ranks and the gate cannot fail", routes_agree)
    print(f"ATTACK GATES PASS: {ok}")

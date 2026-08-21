"""T-50 part (b) MUTATION HARNESS for the hardened C-71/C-72 contrast block.

Applies each mutation of the adversarial suite to the block's model generators and
ASSERTS the designated response of every check on every seed:

  required by the ruling            realised here
  -------------------------------   -----------------------------------------------
  the two-signed write              M1: occupancy write emits +-N_E per programmed
                                        cell (occupancy behaving like orientation)
  a non-zero-mean in-spec residual  M2a-i:  +0.5 e/cell common-mode offset -- the
                                            registrar's own constraint-1 case;
                                    M2a-ii: erased residual U{3..5} e (mean +4, every
                                            cell inside the declared +-5 e), common;
                                    M2b:    the same in-spec shifted residual applied
                                            SECTOR-DIFFERENTIALLY (+-4 e by family) --
                                            the within-part null broken from in spec
  programmed density falling with N M3: a fixed 32-cell record in a growing block
                                        (refuter A's counterexample, verbatim)
  plus, to give every remaining check its measured failing branch:
                                    M4: one-way (DC) data write on orientation
                                    M5: AC-erased DC positive control
                                    M6: carrier sign split across the write

THE ASSERTIONS (any violation exits nonzero):
  * unmutated baseline M0: all 8 checks PASS on all 50 seeds;
  * every (mutation x check) cell is DETERMINISTIC across 50 seeds: fail rate exactly
    1.000 where designated, exactly 0.000 elsewhere -- no 0.18-power cells anywhere;
  * every check has at least one suite mutation that fails it at rate 1.000;
  * constraint 1 measured, both halves: under M2a (common-mode, in-spec) every decision
    check is INVARIANT (PASS 50/50) while the RAW single-sector SHADOW exponent -- which
    is on no decision path -- flips at rate 1.000; under M0 the shadow does not flip;
  * install simulation: the block exec'd inside a validate_geometry.py-style namespace
    (check() predefined, as the registrar will install it) reproduces the module-mode
    run bit-for-bit at the declared seed.

Supplementary MEASURED-ONLY rows (rates reported, nothing asserted -- outcomes are
measured, never failed against an imported standard): M7 in-spec read-long drift;
M8 orientation sector-differential in-spec shift (its bite is bounded by the declared
read tolerance and crosses inside the grid; the measured partial rate is the honest
outcome, and asserting it would tune the gate to the mutation).
"""
import sys, os, time
import collections
import importlib.util
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
BLOCK = os.path.join(HERE, "c72_check_block.py")
spec = importlib.util.spec_from_file_location("c72_check_block", BLOCK)
blk = importlib.util.module_from_spec(spec)
spec.loader.exec_module(blk)

N_SEEDS = 50
CHECK_IDS = ["S1", "S2", "S3", "S4", "S5", "S6", "B2", "A2"]

SUITE = [
    # (row id, description, mutation, designated FAIL set, assert shadow flip?)
    ("M0", "baseline (no mutation)", None, set(), False),
    ("M1", "two-signed write (occupancy behaves like orientation)",
     {"kind": "two_signed_write"}, {"S1", "S4", "S5", "A2"}, False),
    ("M2a-i", "non-zero-mean in-spec residual: +0.5 e/cell COMMON-MODE offset "
              "(the registrar's constraint-1 case; +FRAC/grain on orientation)",
     {"kind": "offset_common", "v_e": 0.5, "v_g": blk.FRAC}, set(), True),
    ("M2a-ii", "non-zero-mean in-spec residual: erased cells U{3..5} e, COMMON "
               "(mean +4 e, every cell inside the declared +-5 e)",
     {"kind": "residual_mean_inspec"}, set(), True),
    ("M2b", "non-zero-mean in-spec residual, SECTOR-DIFFERENTIAL (+-4 e by family; "
            "the within-part null premise broken from inside spec)",
     {"kind": "residual_diff_inspec"}, {"B2", "S5"}, False),
    ("M3", "programmed density falling with N (fixed 32-cell record in a growing "
           "block -- refuter A's counterexample)",
     {"kind": "density_falling"}, {"S1", "S4", "S6"}, False),
    ("M4", "one-way (DC) data write on the orientation surface",
     {"kind": "one_way_data"}, {"S2", "S4"}, False),
    ("M5", "AC-erased DC positive control (two-signed write on the DC sectors)",
     {"kind": "ac_erased_dc"}, {"S3"}, False),
    ("M6", "carrier sign split across the write (half the sectors hole-signed)",
     {"kind": "carrier_split"}, {"S5"}, False),
]

SUPPLEMENTARY = [
    ("M7", "in-spec read-long drift (-Delta -> +Delta ramp across the whole read)",
     {"kind": "drift_inspec"}),
    ("M8", "orientation sector-differential in-spec shift (+-FRAC by family)",
     {"kind": "ori_diff_inspec"}),
]

KILL_TABLE = {   # check -> designated killing mutations (each must reach rate 1.000)
    "S1": ["M1", "M3"], "S2": ["M4"], "S3": ["M5"], "S4": ["M1", "M3", "M4"],
    "S5": ["M1", "M2b", "M6"], "S6": ["M3"], "B2": ["M2b"], "A2": ["M1"],
}


def tag_of(name):
    for cid in CHECK_IDS:
        if f" {cid} " in name:
            return cid
    raise AssertionError(f"unrecognised check name: {name}")


def run_once(row_idx, seed_idx, mutation):
    rec = []

    def ck(name, cond, detail=""):
        rec.append((tag_of(name), bool(cond)))

    st = blk.run_c72_contrast_checks(
        ck, seed=[blk.DECLARED_SEED, row_idx, seed_idx],
        mutation=mutation, verbose=False)
    outcome = dict(rec)
    assert list(outcome) == CHECK_IDS, f"check set changed: {list(outcome)}"
    shadow_flip = (st["xi_shadow"] is not None and st["b_uu_null"] is not None
                   and st["xi_shadow"] > st["b_uu_null"] / 2.0)
    return outcome, shadow_flip, st


def main():
    t0 = time.time()
    print("T-50 HARDENED GATE -- ADVERSARIAL MUTATION MATRIX -- 2026-08-21")
    print(f"block: LANE_T50/GATE/c72_check_block.py   declared seed {blk.DECLARED_SEED}")
    print(f"grid {blk.GRID}  K_PAIRS {blk.K_PAIRS}  DELTA_E {blk.DELTA_E}  "
          f"FRAC {blk.FRAC}  DENSITY_TOL {blk.DENSITY_TOL}  seeds per row {N_SEEDS}")
    print("=" * 100)

    print("\nBASELINE, DECLARED SEED (the block's own full printout):")
    print("-" * 100)
    tally = {"pass": 0, "fail": 0}

    def ck_v(name, cond, detail=""):
        tally["pass" if cond else "fail"] += 1
        print(f"  {'PASS' if cond else 'FAIL'}  {name}  {detail}")

    st_declared = blk.run_c72_contrast_checks(ck_v, verbose=True)
    assert tally == {"pass": 8, "fail": 0}, f"declared-seed baseline not all-PASS: {tally}"
    print(f"  -> declared-seed baseline: {tally['pass']} PASS, {tally['fail']} FAIL")

    print("\nINSTALL SIMULATION (block exec'd with a validate_geometry-style check() "
          "predefined, as the registrar installs it):")
    print("-" * 100)
    sim_rec = []

    def ck_sim(name, cond, detail=""):
        sim_rec.append((name, bool(cond)))

    ns = {"check": ck_sim, "__file__": BLOCK}
    exec(compile(open(BLOCK).read(), BLOCK, "exec"), ns)
    mod_rec = []

    def ck_mod(name, cond, detail=""):
        mod_rec.append((name, bool(cond)))

    blk.run_c72_contrast_checks(ck_mod, verbose=False)
    assert [r for r in sim_rec] == [r for r in mod_rec], "install-sim mismatch"
    print(f"  installed-mode run == module-mode run at declared seed: "
          f"{len(sim_rec)} checks, identical names and verdicts  -> OK")

    # ---------------- the matrix ----------------
    fail_counts = {}
    shadow_flips = {}
    margins = collections.defaultdict(list)
    for ri, (rid, desc, mut, expect_fail, expect_shadow) in enumerate(SUITE):
        fc = collections.Counter()
        sf = 0
        for s in range(N_SEEDS):
            outcome, shadow_flip, st = run_once(ri, s, mut)
            for cid, ok in outcome.items():
                if not ok:
                    fc[cid] += 1
            sf += int(shadow_flip)
            if rid == "M0":
                margins["S6 baseline |dens_slope|"].append(abs(st["dens_slope"]))
                margins["S1 baseline margin"].append(
                    st["xi_occ"] - st["b_uu_null"] / 2)
                margins["SHADOW baseline xi"].append(st["xi_shadow"])
            if rid == "M3":
                margins["S6 under M3 |dens_slope|"].append(abs(st["dens_slope"]))
            if rid == "M1":
                margins["S1 under M1 margin"].append(
                    st["xi_occ"] - st["b_uu_null"] / 2)
            if rid in ("M2a-i", "M2a-ii"):
                margins[f"SHADOW under {rid} xi"].append(st["xi_shadow"])
        fail_counts[rid] = fc
        shadow_flips[rid] = sf

    print("\nMUTATION x CHECK MATRIX (fail count / %d seeds; '*' = designated FAIL):"
          % N_SEEDS)
    print("-" * 100)
    hdr = f"{'row':8s}" + "".join(f"{c:>7s}" for c in CHECK_IDS) + f"{'SHADOW':>9s}"
    print(hdr)
    ok_matrix = True
    for rid, desc, mut, expect_fail, expect_shadow in SUITE:
        cells = []
        for cid in CHECK_IDS:
            n = fail_counts[rid][cid]
            mark = "*" if cid in expect_fail else " "
            cells.append(f"{n:>5d}{mark} ")
            want = N_SEEDS if cid in expect_fail else 0
            if n != want:
                ok_matrix = False
                print(f"  !! NON-DETERMINISTIC CELL {rid} x {cid}: {n}/{N_SEEDS} "
                      f"(designated {want})")
        srow = f"{shadow_flips[rid]:>7d}{'*' if expect_shadow else ' '}"
        print(f"{rid:8s}" + "".join(cells) + srow + f"   {desc}")
        if expect_shadow and shadow_flips[rid] != N_SEEDS:
            ok_matrix = False
            print(f"  !! SHADOW did not flip on every seed under {rid}: "
                  f"{shadow_flips[rid]}/{N_SEEDS}")
    if shadow_flips["M0"] != 0:
        ok_matrix = False
        print(f"  !! SHADOW flipped under baseline: {shadow_flips['M0']}/{N_SEEDS}")
    assert ok_matrix, "matrix is not deterministic at the designated verdicts"
    print(f"  every cell deterministic across {N_SEEDS} seeds: designated FAILs at "
          f"rate 1.000, everything else at rate 0.000  -> OK")

    print("\nFAIL-RATE PER CHECK UNDER ITS DESIGNATED KILLING MUTATION(S):")
    print("-" * 100)
    for cid in CHECK_IDS:
        rates = []
        for rid in KILL_TABLE[cid]:
            r = fail_counts[rid][cid] / N_SEEDS
            rates.append(f"{rid}:{r:.3f}")
            assert fail_counts[rid][cid] == N_SEEDS, f"{cid} under {rid} not 1.000"
        print(f"  {cid}: {'  '.join(rates)}   -> rate 1.000")
    print("  (the old widens check measured 0.18 under the two-signed write; "
          "S1, its successor, measures 1.000)")

    print("\nMARGINS (the decision constants' measured placement):")
    print("-" * 100)
    for k in sorted(margins):
        v = margins[k]
        print(f"  {k:32s} min {min(v):+.4f}  max {max(v):+.4f}")
    b_lo = max(margins["S6 baseline |dens_slope|"])
    m3_hi = min(margins["S6 under M3 |dens_slope|"])
    print(f"  DENSITY_TOL = {blk.DENSITY_TOL}: sits {blk.DENSITY_TOL/b_lo:.1f}x above "
          f"the worst baseline trend ({b_lo:.4f}) and {m3_hi/blk.DENSITY_TOL:.1f}x "
          f"below the smallest M3 trend ({m3_hi:.4f}) -- measured into place")

    # ---------------- supplementary measured-only rows ----------------
    print("\nSUPPLEMENTARY MEASURED-ONLY ROWS (20 seeds; rates reported, "
          "nothing asserted):")
    print("-" * 100)
    for ri, (rid, desc, mut) in enumerate(SUPPLEMENTARY):
        fc = collections.Counter()
        sf = 0
        for s in range(20):
            outcome, shadow_flip, st = run_once(100 + ri, s, mut)
            for cid, ok in outcome.items():
                if not ok:
                    fc[cid] += 1
            sf += int(shadow_flip)
        cells = "  ".join(f"{cid}:{fc[cid]}/20" for cid in CHECK_IDS if fc[cid])
        print(f"  {rid}: {desc}")
        print(f"      fails: {cells if cells else 'none'}   shadow flips: {sf}/20")

    print("=" * 100)
    print(f"VERDICT: baseline all-PASS ({N_SEEDS} seeds + declared seed); matrix "
          f"deterministic; every structural check fails at rate 1.000 under the "
          f"suite; constraint-1 invariance measured with the shadow flipping at "
          f"1.000 under both in-spec non-zero-mean rows; install simulation exact.")
    print(f"elapsed {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()

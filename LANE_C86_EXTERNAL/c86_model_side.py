"""C-86 EXTERNAL-DATA RUN -- MODEL SIDE.  THE COMPARISON, BUILT FROM THE CLOSED FORMS.

This lane runs the MODEL side of the C-86 promotion gate against the two pinned arms
(C86_NAND_PINNED_SOURCES.md, MAGNETIC_ARM_CITATIONS_V001.md).  The closed forms are the
instrument-verified ones from LANE_T47_B_STAIRCASE (instrument-vs-formula <= 6.0e-15 relative
per the ledger row; NOT re-derived here, only evaluated):

    EXACT WIDTH      delta(t_m)  =  kT ln( e^(B/kT) / (f0 t_m) - 1 )
    EXACT STEP       t*_i        =  f0^-1 e^((B_i - dE_i)/kT) / (1 + e^(-dE_i/kT))
    STAIRCASE        k(t_m)      =  #{ i : dE_i <= delta_i(t_m) }  =  #{ i : g_u+g_l <= 1/t_m }
    SYMMETRIC CORNER t*(dE=0)    =  e^(B/kT) / (2 f0)
    BOTH-VALUES FACT (1 + e^(-dE/kT))  in (1, 2]   -- the record dies at the SUM of escapes
    DEPARTURE TERM   M(t) = sum_surviving m_i  +  sum_dead tanh(dE_i / 2kT)
    CLOCK PER RECORD t*_i(T_u)/t*_i(T_b) = e^((B-dE_i)(1/kT_u - 1/kT_b))
                       * (1+e^(-dE_i/kT_b))/(1+e^(-dE_i/kT_u))
                     -> the record's own effective activation energy is B_i - dE_i (level-
                        dependent), plus a bounded correction of at most a factor 2 in time.

DISCIPLINE.  D-8: no fit is performed anywhere in this file; every model number comes from a
closed form evaluated at PINNED constants whose provenance is stated at the point of use;
unpinned constants are never given a value -- they are either swept transparently as a declared
family or the row is scored NOT-COMPARABLE.  D-15: the cycling table (Table 4 of the same
paper) is recomputed as the control separating the retention drive from the cycling drive.
Noise floors are declared before any statement that depends on them.  A mismatched semantic is
scored NOT-COMPARABLE, never bent into agreement (the azobenzene lesson); nulls and mismatches
are recorded with their two-way readings.

THE MAPPING (stated once, used throughout the NAND arm).  Each programmed cell is one record;
its two values are {the written state, the equilibrated state}; B_i = the detrapping barrier
(pinned class value 1.0--1.1 eV); dE_i = alpha * |v_i - v_eq| with v the normalized-Vt level
and alpha the normalized-unit -> eV calibration, VENDOR-WITHHELD (pin 7's own statement), so
alpha is unpinned; f0 is unpinned (no public source pins a NAND attempt frequency).  Any
alpha > 0 preserves order, so ORDINAL predictions are calibration-free; ABSOLUTE positions are
not.  Record death = the written value no longer durable; the observable proxy in pin 7 is
mean drift and adjacent-gap closure (margin-free by that dataset's semantics).

Pins are cited as N-<row> (NAND arm file) and M-<row> (magnetic arm file)."""
import math, sys
def say(*a): print(*a); sys.stdout.flush()
K_EV = 8.617333262e-5          # eV/K, CODATA, exact convention as in c86_nand_arm.py
YR   = 8760.0 * 3600.0         # 1 year = 8760 h (the convention of pin N-6's table)

say("=" * 98)
say("C-86 MODEL SIDE -- k(t_m) AND THE CLOSED FORMS AGAINST THE PINNED ARMS")
say("=" * 98)

# ------------------------------------------------------------------ closed forms (T-47B)
def t_star(B, dE, T, f0):
    kT = K_EV * T
    return math.exp((B - dE) / kT) / (f0 * (1.0 + math.exp(-dE / kT)))
def t_corner(B, T, f0):        # dE = 0 symmetric corner
    return math.exp(B / (K_EV * T)) / (2.0 * f0)
def clock_factor(Ea, Tu, Tb):  # single-Ea bake clock: t_use / t_bake
    return math.exp((Ea / K_EV) * (1.0 / Tu - 1.0 / Tb))

# ================================================================= 1. NAND ORDINAL TESTS
say("")
say("1. NAND ORDINAL PREDICTIONS (calibration-free under the stated mapping) vs pin N-7")
say("   Data: Cai et al., Proc. IEEE 105(9) 2017, Table 5 (transcribed in c86_nand_arm.py;")
say("   TLC, room temperature -- companion study pins 20 C -- 2,000 P/E, normalized-Vt units).")
say("   NOISE FLOOR: table quantization 0.1 units; no statement below 0.1 units is made.")
states = ["ER", "P1", "P2", "P3", "P4", "P5", "P6", "P7"]
mean = {
    "1 day":    [-92.7, 66.6, 128.1, 191.9, 254.9, 318.3, 384.3, 448.1],
    "1 week":   [-86.7, 67.5, 128.1, 191.4, 253.8, 316.5, 381.8, 444.9],
    "1 month":  [-84.4, 68.6, 128.7, 191.6, 253.5, 315.8, 380.9, 443.6],
    "3 months": [-75.6, 72.8, 131.6, 193.3, 254.3, 315.7, 380.2, 442.2],
    "1 year":   [-69.4, 76.6, 134.2, 195.2, 255.3, 316.0, 379.6, 440.8],
}
D = [mean["1 year"][i] - mean["1 day"][i] for i in range(8)]
say("")
say("   MODEL P-ORD-1: signed drift dm_i/dt has sign toward equilibrium and magnitude")
say("   monotone in |dE_i| (both the displacement and the rate g_u+g_l grow with dE), so the")
say("   signed 1d->1y drift must be STRICTLY DECREASING in level with EXACTLY ONE sign change.")
say("   measured drifts: " + "  ".join(f"{states[i]}:{D[i]:+.1f}" for i in range(8)))
steps = [D[i] - D[i + 1] for i in range(7)]
mono = all(s > 0 for s in steps)
signch = sum(1 for i in range(7) if D[i] > 0 >= D[i + 1] or D[i] >= 0 > D[i + 1])
say(f"   decreasing steps: {sum(s>0 for s in steps)}/7, min step {min(steps):.1f} units"
    f" (= {min(steps)/0.1:.0f}x the 0.1 floor); sign changes: {signch} (between P4 and P5)")
say(f"   -> P-ORD-1 {'HOLDS' if mono and signch==1 else 'FAILS'} on pin N-7, above floor.")
say("")
say("   MODEL P-ORD-2: the first-closing pair is the pair containing the record farthest from")
say("   equilibrium (largest |dE|): ER-P1 overall; among programmed pairs, the pair farthest")
say("   from the zero crossing on the fast side: P1-P2.")
gaps = {}
for i in range(7):
    g0 = mean["1 day"][i + 1] - mean["1 day"][i]
    g1 = mean["1 year"][i + 1] - mean["1 year"][i]
    gaps[f"{states[i]}-{states[i+1]}"] = g1 - g0
say("   measured d(gap)/yr: " + "  ".join(f"{k}:{v:+.1f}" for k, v in gaps.items()))
fast_all = min(gaps, key=gaps.get)
fast_prog = min((k for k in gaps if not k.startswith("ER")), key=gaps.get)
say(f"   fastest-closing overall: {fast_all}; fastest programmed: {fast_prog}"
    f"  -> P-ORD-2 {'HOLDS' if (fast_all=='ER-P1' and fast_prog=='P1-P2') else 'FAILS'}.")
say("")
say("   HONESTY LINE: P-ORD-1/2 are shared with the ENTIRE activated-relaxation class the")
say("   C-86 row concedes by name (Neel/Street-Woolley/Sharrock; MLC margin engineering).")
say("   Their agreement supports the model's INPUTS; it does not discriminate the owned law.")

# ---------------------------------------------------------------- small-t shape, no fit
say("")
say("2. SMALL-t SHAPE (pin N-7 five-point series; pin N-10 pins log-t from 7 min on 3D)")
say("   MODEL: a single record approaches equilibrium as 1 - e^(-t/tau) (one rate); a state's")
say("   mean is an ensemble over ~1e6+ cell-records; with dispersed B_i the ensemble sum gives")
say("   log-t over a window -- but that form is CONCEDED (Street-Woolley viscosity class) and")
say("   the dispersion width is unpinned.  Here: slopes d(drift)/d(ln t), displayed, NO FIT.")
tdays = {"1 day": 1.0, "1 week": 7.0, "1 month": 30.0, "3 months": 91.0, "1 year": 365.0}
say("   (time labels quantized to days as printed; '3 months' taken as 91 d -- label floor)")
for s_i, s in [(0, "ER"), (1, "P1"), (7, "P7")]:
    ts = list(tdays)
    sl = []
    for a, b in zip(ts, ts[1:]):
        dd = mean[b][s_i] - mean[a][s_i]
        dl = math.log(tdays[b] / tdays[a])
        sl.append(dd / dl)
    say(f"   {s}: interval slopes (units per ln t): " + ", ".join(f"{x:+.2f}" for x in sl))
say("   slopes are non-monotone at the 5-point resolution: the series separates NEITHER the")
say("   single-exponential approach NOR pure log-t.  Shape row -> not deciding at this floor.")

# ---------------------------------------------------------------- D-15 control
say("")
say("3. D-15 CONTROL (pin N-7 paper, Table 4: cycling 0 -> 3,000 P/E at fixed 1-day retention)")
cyc0  = [-110.0, 65.9, 127.4, 191.6, 254.9, 318.4, 384.8, 448.3]
cyc3k = [-84.1, 68.3, 128.2, 193.1, 255.7, 319.2, 385.4, 449.1]
dc = [cyc3k[i] - cyc0[i] for i in range(8)]
say("   0->3k P/E drift: " + "  ".join(f"{states[i]}:{dc[i]:+.1f}" for i in range(8)))
say("   ALL programmed states move UP, no sign change, not monotone-decreasing: the cycling")
say("   drive fails P-ORD-1's signature, so section 1's structure is retention's own.")

# ================================================================= 4. NAND ABSOLUTE k(t_m)
say("")
say("4. THE ABSOLUTE STAIRCASE k(t_m) FOR THIS DEVICE -- the executable corner, as a family")
say("   Pinned: B = 1.0 eV (pin N-5, CLASS) .. 1.1 eV (pins N-2/N-4/N-6, the JEDEC-class and")
say("   Micron-primary detrapping value); T = 293.15 K (20 C, pin N-8's explicit statement,")
say("   pin N-7 'room temperature').  UNPINNED: f0 (no public NAND pin exists) and alpha (the")
say("   normalized-unit -> eV calibration; withheld by pin N-7's own text).  dE_i > 0 only")
say("   moves a record EARLIER (factor (1+e^(dE_i/kT))/2 > 1 in rate), so the dE = 0 corner")
say("   t* = e^(B/kT)/(2 f0) is the LATEST possible first drop: k(t_m) = 8 requires")
say("   t_m <= t*(corner).  The family, swept transparently (D-8: a sweep, not a fit):")
T_RT = 293.15
kT_RT = K_EV * T_RT
say(f"   kT = {kT_RT:.6f} eV;  e^(B/kT) = {math.exp(1.0/kT_RT):.3e} (B=1.0), "
    f"{math.exp(1.1/kT_RT):.3e} (B=1.1)")
say(f"   {'f0 (Hz)':>10} | {'t*corner B=1.0 eV':>20} | {'t*corner B=1.1 eV':>20}")
def human(t):
    for u, s in (("s", 1.0), ("h", 3600.0), ("d", 86400.0), ("yr", YR)):
        pass
    if t < 3600: return f"{t:.2e} s"
    if t < 86400: return f"{t/3600:.1f} h"
    if t < YR: return f"{t/86400:.1f} d"
    return f"{t/YR:.1f} yr"
for f0 in (1e9, 1e10, 1e11, 1e12, 1e13):
    a, b = t_corner(1.0, T_RT, f0), t_corner(1.1, T_RT, f0)
    say(f"   {f0:>10.0e} | {human(a):>20} | {human(b):>20}")
say("")
say("   MEASURED (pin N-7, margin-free): k(1 yr, RT, 2k P/E) = 8 of 8 state-classes alive --")
say("   no adjacent pair merges inside the public window.  Model + pinned B then EXCLUDES:")
for B in (1.0, 1.1):
    f0max = math.exp(B / kT_RT) / (2.0 * YR)
    say(f"     B = {B:.1f} eV: survival of the dE=0 corner to 1 yr requires f0 <= {f0max:.2e} Hz"
        f"  (tighter by (1+e^(dE/kT))/2 for every dE_i > 0)")
say("   The textbook attempt-frequency decade 1e12-1e13 Hz is EXCLUDED at both pinned B; the")
say("   surviving region is physically ordinary (<= ~1e11 Hz at B = 1.1).  This is a")
say("   CONSTRAINT the data places on the model's unpinned constant, not an agreement -- with")
say("   f0 and alpha both unpinned, the absolute staircase has no public-data comparison.")

# ================================================================= 5. THE BAKE CLOCK
say("")
say("5. THE BAKE-CLOCK MAPPING AND ITS BAND (the semantic gate for every bake-positioned row)")
say("   Pinned mapping (pins N-2, N-6): t_use = t_bake * exp(Ea/k (1/T_use - 1/T_bake)),")
say("   Ea = 1.100 +- 0.001 eV (table inversion) with the CLASS pin 1.0 eV (N-5) as the low")
say("   edge of the pinned band.  Band of the mapping across the pinned Ea spread:")
Tu = 298.15
say(f"   {'T_bake':>8} | {'accel @1.0 eV':>14} | {'accel @1.1 eV':>14} | {'band ratio':>10}")
for TbC in (60, 80, 100, 120):
    Tb = TbC + 273.15
    a0, a1 = clock_factor(1.0, Tu, Tb), clock_factor(1.1, Tu, Tb)
    say(f"   {TbC:>6} C | {a0:>14.3e} | {a1:>14.3e} | {a1/a0:>10.2f}x")
say("   MODEL'S OWN CLOCK: per record the effective Ea is B_i - dE_i (level-dependent), plus")
say("   a two-sided correction bounded by a factor 2 in time.  A single-Ea clock is exact")
say("   ONLY if every record shares B - dE; with dE_i > 0 spanning levels, the model says the")
say("   1.1 eV bake under-delivers use-aging for EVERY displaced record (its true Ea is")
say("   smaller), i.e. bake-certified time OVERSTATES delivered aging -> bakes UNDERESTIMATE")
say("   real retention errors, worst for the largest-dE (outer) levels.")
say("   MEASURED (pin N-9): bakes underestimate the real 366-day errors, and per-read-voltage")
say("   error ordering CHANGES with temperature.  SIGN AGREES.  ATTRIBUTION CONFOUNDED: the")
say("   sources themselves (pins N-3, N-8, N-11) blame two mechanisms with different Ea, a")
say("   channel the single-escape model does not contain; and the pinned 1.1 eV was itself")
say("   measured from Delta-VT transients, which under the model mixes B - dE_i values.  So")
say("   the sign-level agreement is a SHARED NULL of the single-Ea clock -- it cannot")
say("   discriminate C-86 from mechanism mixing.  Bake-positioned ABSOLUTE times carry a")
say("   multiplicative band >= 1.5-2.6x (Ea spread) x <= 2x (two-sided correction) BEFORE the")
say("   unbounded measured mechanism systematic: absolute bake rows are NOT-COMPARABLE.")

# ================================================================= 6. MAGNETIC ARM
say("")
say("6. MAGNETIC ARM -- model numbers where the pins admit any")
say("")
say("6a. Pin M-1 (Wernsdorfer N=1 Co particle).  MODEL: the record's own Liouvillian mode is")
say("    ONE exponential at every (T, H) -- measured: exponential at every (T, H) pinned.")
say("    FORM AGREES; but this is the field-tilted Neel-Brown fragment CONCEDED BY NAME, and")
say("    the pinned absolutes (E0 = 214,000 K, tau0 = 3e-9 s) are OUTPUTS of fitting that same")
say("    form -- any absolute 'comparison' is circular by construction.  No census at N = 1.")
say("")
say("6b. Pin M-2 (Krause Fe/W(110) per-island ensemble).  Pinned: E_b = E0 + e_DW*N,")
say("    E0 = 61+-5 meV, e_DW = 7.5+-0.4 meV/row, worked T = 53.6 K; nu0 spans 1e13-1e16 Hz")
say("    WITHIN one ensemble.  The shared-f0 closed form is scope-excluded by C-86's own")
say("    register note; this data CONFIRMS the exclusion and the stress is computable:")
kT_K = K_EV * 53.6
per_row = 7.5e-3 / kT_K
spread = math.log(1e16 / 1e13)
say(f"    Delta E_b per atomic row / kT = {per_row:.2f}; ln(nu0 spread) = {spread:.2f}")
say(f"    -> a forced single f0 can scramble the death ORDER across up to "
    f"{spread/per_row:.1f} adjacent-N islands,")
say(f"    and mispositions t*_i by up to e^{spread:.1f} ~ 1e3 in time whichever f0 is chosen")
say(f"    (half-spread floor: >= e^{spread/2:.2f} = {math.exp(spread/2):.0f}x).  Real-media")
say("    staircases must run through the record-mode instrument (per-record f0_i), exactly as")
say("    the register's scope note declares.  No census was measured -> owned law untouched.")
say("")
say("6c. Pin M-3 (SMTJ, Funatsu 2022; raw dwell times public, Zenodo 6767828).  The one")
say("    public class carrying BOTH escape rates AND occupancy on one record.  MODEL NUMBERS")
say("    (parameter-free, the wholly-owned both-values factor at the step):")
say(f"    {'dE/kT':>8} | {'1+e^(-dE/kT)':>13} | {'tanh(dE/2kT)':>13}")
for x in (0.0, 1.0, 2.0, 3.0, 4.0):
    say(f"    {x:>8.1f} | {1.0+math.exp(-x):>13.4f} | {math.tanh(x/2.0):>13.4f}")
say("    At the balanced field the record lifetime is HALF the per-branch dwell time -- the")
say("    factor-2 corner that separates C-86's t* from the favored-branch (Sharrock-class)")
say("    convention, decaying to 1 as e^(-dE/kT).  The paper verifies exponential/Poisson")
say("    intervals (form AGREES -- conceded fragment) but its Delta values ASSUME tau0 = 1 ns,")
say("    so comparing t* built from those Deltas against the same dwell data is CIRCULAR (the")
say("    azobenzene trap, rate-level form).  The non-circular test -- occupancy and summed")
say("    rate from the RAW dwell sequences with tau0 marginalized -- needs the Zenodo download,")
say("    which this run did not perform (file-download gate); it stays the named next step.")
say("    THIS RUN INGESTS NO SMTJ NUMBERS -> row scored NOT-COMPARABLE (this run).")
say("")
say("6d. Pins M-4/M-5 (ASI PEEM census movies / SQUID remanence twin).  The ONLY public")
say("    count-grade time series.  MODEL PREREQUISITE: dE_i is the neighbors' dipolar field")
say("    and changes as neighbors flip; C-86's staircase is stated for fixed (B_i, dE_i), so")
say("    no number may be compared before the state-dependent-dE extension exists.  The")
say("    moment-side twin (M-5) is on DIFFERENT samples: cross-sample comparison is")
say("    semantics-broken.  NOT-COMPARABLE (extension required / cross-sample).")
say("")
say("6e. Pins M-6/M-7/M-8 (York-protocol 'grain census' / MFM heated-stage / multi-barrier")
say("    particle).  The commissioned census observable does not exist in these classes: the")
say("    field's census is a FITTED volume distribution (never an integer count), MFM signal")
say("    is a moment proxy, and real grains are not guaranteed two-state (M-8) -- the")
say("    two-valuedness check is a precondition the model imposes on any future census run.")
say("    NOT-COMPARABLE (absent observable / proxy semantics / scope hazard).")
say("")
say("6f. THE DEPARTURE TERM sum_dead tanh(dE_i/2kT) -- the wholly-owned discriminator")
say("    (remanence persists while records die).  It needs census AND remanence on ONE sample")
say("    with per-record constants.  NO pinned dataset on either arm carries both.  The term")
say("    is evaluated against nothing in this run; no hypothetical value is printed (D-8).")

# ================================================================= 7. THE COMPARISON TABLE
say("")
say("=" * 98)
say("7. THE COMPARISON TABLE -- measured beside model, tolerance, verdict, why")
say("=" * 98)
rows = [
 ("C1", "N-7 Cai'17 Table 5: signed drift order",
  "7/7 steps decreasing, 1 sign change (P4|P5); min step 2.4u",
  "strictly decreasing, exactly 1 crossing (P-ORD-1)",
  "floor 0.1u; min step 24x floor", "AGREES (ordinal)",
  "calibration-free under the mapping; but shared with the conceded activated-relaxation"
  " class -- supports inputs, does not discriminate the owned law"),
 ("C2", "N-7: first-closing pair",
  "ER-P1 fastest (-13.3u/yr); P1-P2 fastest programmed (-3.9u/yr)",
  "pair holding the extreme record: ER-P1; then P1-P2 (P-ORD-2)",
  "floor 0.1u", "AGREES (ordinal)",
  "same non-discrimination caveat as C1"),
 ("C3", "N-7: absolute first drop t*_1",
  "no pair merges within 1 yr at RT/2k P/E (k = 8 of 8 alive)",
  "t*_1 = f0^-1 e^((B-dE_max)/kT)/(1+e^(-dE_max/kT)); needs f0, alpha",
  "f0, alpha vendor-withheld", "NOT-COMPARABLE (absolute)",
  "constraint recorded instead: with pinned B, survival to 1 yr forces f0 <= 2.5e9 Hz"
  " (B=1.0) .. 1.3e11 Hz (B=1.1) at the dE=0 corner; 1e12-1e13 Hz excluded"),
 ("C4", "N-8 HPCA'15 MLC: down-side order",
  "dP3 > dP2 > dP1, P1 ~constant (figure-resolution)",
  "drift magnitude monotone in level on the down side",
  "figure-resolution", "AGREES (ordinal)",
  "same class caveat; the paper itself rejects bake acceleration on record"),
 ("C5", "N-10 SIGMETRICS'18: small-t form",
  "V = A log t + B from 7 min; RBER decade in ~3 h (3D CT)",
  "single record: exponential approach; ensemble with dispersed B: log-t (conceded form,"
  " dispersion unpinned)",
  "anonymized axes", "NOT-COMPARABLE (owned-law)",
  "the log-t form is Street-Woolley-class conceded territory; matching it tests nothing"
  " wholly owned; N-7's 5-point series separates neither form (sec 2)"),
 ("C6", "N-9 HotStorage'21: bake vs real year",
  "bakes UNDERESTIMATE 366-day errors; per-voltage ordering changes with T",
  "single-Ea clock must fail; direction: bake under-delivers aging for every dE_i > 0",
  "figure-resolution", "AGREES (sign only)",
  "shared null of the single-Ea clock -- mechanism mixing (N-3/N-8/N-11) predicts failure"
  " too; level-resolved reordering not extractable at figure resolution"),
 ("C7", "N-11 Malavena fresh-3D bake: level structure",
  "Delta-VT NONMONOTONIC in level (depassivation + detrapping, opposite signs)",
  "single-channel fixed-(B,dE) predicts monotone",
  "figure-resolution", "NOT-COMPARABLE (premise broken)",
  "the source itself identifies a second, chemically distinct channel: the fixed-constants"
  " premise fails on fresh 3D; REFUTES any naive single-channel application -- the channel"
  " check is now a mandatory precondition. Post-cycled data (same pin): monotone, agrees"
  " ordinally at figure resolution"),
 ("C8", "N-1/N-12 JEDEC & Mielke class: retirement",
  "UBER/FFR at vendor margins, ECC-downstream",
  "margin-free integer census k(t_m)",
  "--", "NOT-COMPARABLE (semantics)",
  "margin-laden drive-level observables are not a census; per-level structure unpublished"),
 ("C9", "N-2/N-5/N-6: the bake clock itself",
  "Ea = 1.100+-0.001 eV (table inversion) vs 1.0 eV class pin",
  "per-record effective Ea = B_i - dE_i, level-dependent; correction <= 2x",
  "band 1.5x (60 C) - 2.6x (120 C)", "NOT-COMPARABLE (absolute gate)",
  "the clock is one contested convention; every bake-positioned absolute inherits >= 3-5x"
  " multiplicative band before the unbounded mechanism systematic"),
 ("C10", "M-1 Wernsdorfer N=1: escape form",
  "P(t) exponential at every (T,H); telegraph residences exponential",
  "record mode = ONE exponential",
  "no bars pinned on tau0", "AGREES (form; conceded)",
  "field-tilted Neel-Brown is conceded by name; absolutes (E0, tau0) are outputs of the"
  " same form -- circular; no census at N=1"),
 ("C11", "M-2 Krause ensemble: shared-f0 form",
  "nu0 spans 1e13-1e16 Hz within ONE ensemble",
  "closed form assumes one f0 (declared scope); stress computed: order scrambled up to"
  " ~4 adjacent-N islands, positions off up to ~1e3x",
  "E0 +-5 meV, e_DW +-0.4 meV/row", "NOT-COMPARABLE (scope, confirmed)",
  "the data confirms the register's own scope exclusion; record-mode instrument route"
  " required for real media; no census measured"),
 ("C12", "M-3 SMTJ Zenodo: both-values factor & occupancy",
  "tau_P, tau_AP, occupancy public (raw); exponential/Poisson verified",
  "1/t* = g_u + g_l (factor (1+e^(-dE/kT)) in (1,2]); <m> = tanh(dE/2kT) -- tabulated",
  "tau0 = 1 ns ASSUMED in source", "NOT-COMPARABLE (this run)",
  "the ONLY located executable test of an owned corner (the both-values factor); needs the"
  " raw download with tau0 marginalized -- not performed this run (download gate); paper-"
  "level Delta comparison would be circular"),
 ("C13", "M-4/M-5 ASI census movies / remanence twin",
  "integer vertex counts vs t (PEEM); remanence on different samples",
  "staircase stated for FIXED (B_i, dE_i); dipolar dE_i is state-dependent",
  "--", "NOT-COMPARABLE (extension req.)",
  "only public count-grade series; comparable only after a state-dependent-dE extension,"
  " else any null reads two ways; cross-sample moment comparison semantics-broken"),
 ("C14", "M-6/M-7/M-8 exchange-bias census / MFM / multi-barrier",
  "fitted volume distributions; moment proxies; >2-state reversal",
  "integer survivor count of verified two-state records",
  "--", "NOT-COMPARABLE (absent/proxy/scope)",
  "the commissioned census observable does not exist publicly in this medium; producing it"
  " is lab work; two-valuedness must be verified per record first"),
 ("C15", "BOTH ARMS: departure term sum_dead tanh(dE_i/2kT)",
  "no dataset carries census AND remanence on one sample",
  "M(t) = sum_surv m_i + sum_dead tanh(dE_i/2kT)",
  "--", "NOT-COMPARABLE (no data)",
  "the wholly-owned discriminator is touched by NOTHING located; this is the promotion"
  " gate's operative gap"),
]
for r in rows:
    say("")
    say(f"  {r[0]:>4}  {r[1]}")
    say(f"        measured : {r[2]}")
    say(f"        model    : {r[3]}")
    say(f"        tolerance: {r[4]}")
    say(f"        VERDICT  : {r[5]}")
    say(f"        why      : {r[6]}")

# ================================================================= 8. VERDICT LINES
say("")
say("=" * 98)
say("8. MODEL-SIDE VERDICT (the comparison comes out as it comes out)")
say("=" * 98)
say("  * Score: 5 AGREES (C1, C2, C4, C6, C10) -- every one either ordinal/sign-level or a")
say("    conceded fragment; 10 NOT-COMPARABLE; 0 DISAGREES.  No agreement row discriminates")
say("    the wholly-owned content (integer census, both-values factor, departure term) from")
say("    the fragments C-86 concedes by name.  The owned content is touched by NO located")
say("    public dataset.")
say("  * The run's real products: (i) the f0 exclusion (pinned B + measured 1-yr survival")
say("    excludes f0 >= ~1e11-1e12 Hz at RT -- a falsifiable constraint any future calibrated")
say("    dataset must respect); (ii) the quantified bake-clock band (1.5-2.6x from the Ea")
say("    pins, <= 2x two-sided, unbounded mechanism systematic) gating all absolute bake")
say("    comparisons; (iii) data-confirmation of the shared-f0 scope exclusion (order")
say("    scrambling up to ~4 neighbors, ~1e3x positions); (iv) the mandatory channel check")
say("    (C7): fixed-single-channel premises FAIL on fresh 3D by the sources' own mechanism")
say("    decomposition.")
say("  * PROMOTION: this run does NOT support FORMAL -> PROVED.  The strengthened bar asked")
say("    the external data to reach the owned law; the located data reaches only conceded")
say("    territory and constraints.  Both readings remain open: law-untested (the data never")
say("    enters the census regime) vs law-supported (no located datum contradicts it beyond")
say("    the scope-gated C7/C11 rows).")
say("  * NEXT STEPS (named): (1) the Zenodo 6767828 grounding run -- occupancy vs summed-rate")
say("    with tau0 marginalized: the only located executable test of an owned corner;")
say("    (2) the state-dependent-dE extension, then the ASI PEEM re-analysis (census and")
say("    moment from the same frames -- the only public route to the departure term);")
say("    (3) the in-house read-retry census (pin N-7/N-9 methodology, retail chips) for")
say("    margin-free absolute staircase positions and the f0-exclusion check.")
sys.exit(0)

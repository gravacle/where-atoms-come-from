"""C-86 EXTERNAL-DATA RUN, NAND/MLC ARM -- ARITHMETIC ON PINNED PUBLIC NUMBERS ONLY.

No fit is performed anywhere in this file (D-8).  Every number is transcribed from a pinned
source listed in C86_NAND_PINNED_SOURCES.md; the script only differences, ratios, and inverts
published tables, and states the noise floor of each arithmetic before using it.

Three instruments:
  1. INDUSTRY-Ea INVERSION (control): the bake-equivalence table published by Li et al.
     (HotStorage 2021, slide 5) is inverted to the activation energy it silently assumes.
     If the industry table and JESD218B Annex A (which states 1.1 eV in text) disagree,
     the JEDEC-class 'activation energy' pin is not a single number and must be reported so.
  2. LEVEL-RESOLVED DRIFT from Cai/Ghose/Haratsch/Luo/Mutlu (Proc. IEEE 105(9) 2017,
     Table 5): per-state mean/sigma of TLC threshold-voltage distributions at 5 retention
     times, room temperature, 2,000 P/E.  Computed margin-free: adjacent-gap trajectories,
     gap-in-sigma-units trajectories, and whether any pair approaches merger inside the
     public window.  NOISE FLOOR: table quantized at 0.1 normalized units; per-state sigma
     9.1-12.8 units; any statement below 0.1 units is unreadable and none is made.
  3. CONTROL (D-15): Table 4 of the same paper (P/E cycling at fixed 1-day retention) --
     cycling moves programmed states UP, retention moves the high states DOWN, so the
     instrument distinguishes the two drives before any retention statement is trusted.
"""
import math, sys
def say(*a): print(*a); sys.stdout.flush()
K_EV = 8.617333262e-5  # eV/K

say("=" * 98)
say("C-86 NAND/MLC ARM -- WHAT THE PUBLIC NUMBERS SAY, COMPUTED, WITH FLOORS")
say("=" * 98)

# ---------------------------------------------------------------- 1. industry-Ea inversion
say("")
say("1. THE BAKE-EQUIVALENCE TABLE THE INDUSTRY ACTUALLY USES, INVERTED TO ITS Ea")
say("   Source: Li/Ye/Kuo/Xue, HotStorage 2021 slides p.5 (public PDF): 1 year at 25 C is")
say("   declared equivalent to the following bake hours:")
T0 = 298.15
t0_h = 8760.0
table = [(60, 97.65), (80, 11.16), (100, 1.61), (120, 0.28)]
say(f"   {'T_bake (C)':>12}{'hours':>10}{'implied Ea (eV)':>18}")
eas = []
for TC, th in table:
    T = TC + 273.15
    ea = math.log(t0_h / th) * K_EV / (1.0 / T0 - 1.0 / T)
    eas.append(ea)
    say(f"   {TC:>12}{th:>10}{ea:>18.3f}")
say(f"   spread {min(eas):.3f}-{max(eas):.3f} eV; JESD218B 6.1.3 states its Annex A stress-time")
say("   calculation assumes 1.1 eV IN TEXT.  The published table and the standard agree to the")
say("   0.28-hour quantization of the table's last row: ONE number, 1.1 eV, is what the")
say("   JEDEC-class world runs on -- and JESD218B 3.19 itself says the OTHER retention")
say("   mechanism (leak through the transfer dielectric) 'can be weakly accelerated or even")
say("   decelerated by high temperature', i.e., the standard's own text denies that a single")
say("   Ea covers retention.  The 1.1 eV is the DETRAPPING channel only.")

# ------------------------------------------------- 2. level-resolved drift, Cai Table 5
say("")
say("2. LEVEL-RESOLVED RETENTION DRIFT (Cai et al., Proc. IEEE 105(9) 2017, Table 5; TLC,")
say("   room temperature, 2,000 P/E; normalized-Vt units, nominal max = 512, GND = 0;")
say("   absolute volts are vendor-proprietary BY THE PAPER'S OWN STATEMENT)")
states = ["ER", "P1", "P2", "P3", "P4", "P5", "P6", "P7"]
times  = ["1 day", "1 week", "1 month", "3 months", "1 year"]
mean = {
    "1 day":    [-92.7, 66.6, 128.1, 191.9, 254.9, 318.3, 384.3, 448.1],
    "1 week":   [-86.7, 67.5, 128.1, 191.4, 253.8, 316.5, 381.8, 444.9],
    "1 month":  [-84.4, 68.6, 128.7, 191.6, 253.5, 315.8, 380.9, 443.6],
    "3 months": [-75.6, 72.8, 131.6, 193.3, 254.3, 315.7, 380.2, 442.2],
    "1 year":   [-69.4, 76.6, 134.2, 195.2, 255.3, 316.0, 379.6, 440.8],
}
# Table 5 sigma rows: the public PDF text-extraction truncates after P5 at "10.8"; the
# 1-day row is complete and P6/P7 1-year values are NOT visible in the extraction.
# Only the rows actually read are used; nothing is interpolated.
sigma_1day = [48.2, 9.7, 9.7, 9.4, 9.3, 9.1, 9.3, 8.5]   # ER..P7, complete row as printed
sigma_1yr_partial = {"P1": 12.8, "P2": 12.4, "P3": 12.0, "P4": 12.0}  # visible cells only
say("")
say("   NOISE FLOOR: quantization 0.1 units; drift statements below 0.1 units are not made.")
say("   Per-state sigma at 1 day: " + ", ".join(f"{s}={v}" for s, v in zip(states, sigma_1day)))
say("")
say("   (a) per-state mean drift, 1 day -> 1 year (units):")
d = [mean["1 year"][i] - mean["1 day"][i] for i in range(8)]
say("       " + "  ".join(f"{states[i]}:{d[i]:+.1f}" for i in range(8)))
say("       sign structure: ER,P1,P2,P3 UP; P4 ~flat (+0.4); P5,P6,P7 DOWN, magnitude")
say("       increasing with level (P7 largest down-drift).  The paper's own reading: TAT")
say("       (field-driven, only downward) dominates high states; hole loss lifts low states.")
say("")
say("   (b) adjacent-gap trajectories, margin-free (gap = mean[i+1]-mean[i], units):")
say(f"       {'pair':>8}" + "".join(f"{t:>11}" for t in times) + f"{'d(gap)/yr':>12}")
for i in range(7):
    row = [mean[t][i + 1] - mean[t][i] for t in times]
    say(f"       {states[i]+'-'+states[i+1]:>8}" + "".join(f"{g:>11.1f}" for g in row)
        + f"{row[-1]-row[0]:>+12.1f}")
say("       every programmed-pair gap SHRINKS; the fastest-closing programmed pair over the")
say("       full year is P1-P2; ER-P1 closes fastest of all (ER rises 23.3).")
say("")
say("   (c) resolvability in sigma units (gap / (sigma_i+sigma_j)), where sigma is public:")
for pair, (i, j) in {"P1-P2": (1, 2), "P2-P3": (2, 3), "P3-P4": (3, 4)}.items():
    g0 = mean["1 day"][j] - mean["1 day"][i]
    g1 = mean["1 year"][j] - mean["1 year"][i]
    s0 = sigma_1day[i] + sigma_1day[j]
    r0 = g0 / s0
    out = f"       {pair}: 1-day gap {g0:.1f} = {r0:.1f} combined-sigma"
    a, b = states[i], states[j]
    if a in sigma_1yr_partial and b in sigma_1yr_partial:
        s1 = sigma_1yr_partial[a] + sigma_1yr_partial[b]
        out += f"; 1-year gap {g1:.1f} = {g1/s1:.1f} combined-sigma"
    say(out)
say("       NO adjacent pair approaches merger inside the public window: the first")
say("       retirement event -- the first step of C-86's staircase -- lies BEYOND one year")
say("       at room temperature and 2,000 P/E for this chip.  Public data shows the")
say("       APPROACH to the first drop, never a drop.")

# ---------------------------------------------------------------- 3. control (D-15)
say("")
say("3. CONTROL: the same paper's Table 4 (cycling at fixed 1-day retention, 0 -> 3,000 P/E)")
cyc0 = [-110.0, 65.9, 127.4, 191.6, 254.9, 318.4, 384.8, 448.3]
cyc3k = [-84.1, 68.3, 128.2, 193.1, 255.7, 319.2, 385.4, 449.1]
dc = [cyc3k[i] - cyc0[i] for i in range(8)]
say("       0->3k P/E drift: " + "  ".join(f"{states[i]}:{dc[i]:+.1f}" for i in range(8)))
say("       cycling moves EVERY programmed state UP (more injection through degraded oxide);")
say("       retention moves the high states DOWN.  Opposite signatures on P5-P7: the")
say("       instrument separates the two drives, so 2(a)'s sign structure is retention's own.")
say("")
say("=" * 98)
say("VERDICT LINES (data side only; the model-side comparison is a later lane's job):")
say("  * Level-resolved, margin-free PUBLIC data exists (per-state mean/sigma vs time) but in")
say("    NORMALIZED units at ROOM temperature, and the staircase's first drop is outside the")
say("    published window.  JEDEC-class BAKE data with level-resolved retirement positions is")
say("    not public at any granularity found; drive-level bake results are UBER/margin-defined.")
say("  * The absolute calibration the parameter-free drop-time formula needs (B_i, dE_i in eV,")
say("    f0) is vendor-proprietary in every located source.  Order/shape tests are executable;")
say("    absolute t*_i positions are not, on public data as found.")
sys.exit(0)

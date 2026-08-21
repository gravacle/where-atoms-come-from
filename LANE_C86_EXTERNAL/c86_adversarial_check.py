"""C-86 EXTERNAL-DATA RUN -- ADVERSARIAL CHECK (default REFUTED), 2026-08-21.

Independent re-fetch and re-derivation pass over the two pinned arms
(C86_NAND_PINNED_SOURCES.md, MAGNETIC_ARM_CITATIONS_V001.md) and the model-side comparison
(c86_model_side.py).  Sources RE-FETCHED FIRST-HAND this run, from outside the repository:

  R1. Cai/Ghose/Haratsch/Luo/Mutlu, Proc. IEEE 105(9) 2017 -- arXiv:1706.08642 PDF downloaded,
      Appendix Tables 4 and 5 extracted from page 34 of the PDF (pypdf), full text searched.
  R2. Cai/Luo/Haratsch/Mai/Mutlu, HPCA 2015 -- ETH open PDF downloaded, text searched.
  R3. Li/Ye/Kuo/Xue, HotStorage 2021 slides -- hotstorage.org/2021/2021-slides/
      How_the_Common-SLIDES-LIQiao.pdf downloaded, all 18 slides extracted.
  R4. Malavena et al., IEEE JEDS 11 (2023) -- re.public.polimi.it open PDF downloaded.
  R5. JEDEC JC-64.8 chair's public deck (Alvin Cox, jedec.org) -- downloaded; carries the
      JESD218 class table and the 1.1 eV / second-mechanism language first-hand from JEDEC.
  R6. Wernsdorfer review arXiv:cond-mat/0101104 -- PDF downloaded, Sec. 3.2-3.3 read.
  R7. Hamburg SPM Triannual Report 2008-2010 (Krause et al. section) -- PDF downloaded.
  R8. Funatsu et al., Nat. Commun. 13 (2022) -- nature.com HTML fetched (open access);
      Zenodo record 6767828 landing page fetched.
  R9. Farhan PRL 111, 057204 existence/abstract via APS/PubMed (23952441).

No fit is performed; every number below is a transcription from a re-fetched source or a
closed-form evaluation at stated constants.  D-8/D-15 as in the sibling instruments.
"""
import math, sys
def say(*a): print(*a); sys.stdout.flush()
K_EV = 8.617333262e-5
YR = 8760.0 * 3600.0
say("=" * 98)
say("C-86 ADVERSARIAL CHECK -- RE-FETCH, RE-DERIVE, ATTACK; DEFAULT REFUTED")
say("=" * 98)

# ----------------------------------------------------- 1. pin N-7 transcription, cell-by-cell
say("")
say("1. PIN N-7 (Cai'17 Table 5) -- LANE TRANSCRIPTION vs THE PDF RE-FETCHED THIS RUN")
states = ["ER", "P1", "P2", "P3", "P4", "P5", "P6", "P7"]
# means as transcribed in c86_nand_arm.py / c86_model_side.py:
lane_mean = {
    "1 day":    [-92.7, 66.6, 128.1, 191.9, 254.9, 318.3, 384.3, 448.1],
    "1 week":   [-86.7, 67.5, 128.1, 191.4, 253.8, 316.5, 381.8, 444.9],
    "1 month":  [-84.4, 68.6, 128.7, 191.6, 253.5, 315.8, 380.9, 443.6],
    "3 months": [-75.6, 72.8, 131.6, 193.3, 254.3, 315.7, 380.2, 442.2],
    "1 year":   [-69.4, 76.6, 134.2, 195.2, 255.3, 316.0, 379.6, 440.8],
}
# means as read from the re-fetched PDF (arXiv:1706.08642, p.34, Table 5 top block):
pdf_mean = {
    "1 day":    [-92.7, 66.6, 128.1, 191.9, 254.9, 318.3, 384.3, 448.1],
    "1 week":   [-86.7, 67.5, 128.1, 191.4, 253.8, 316.5, 381.8, 444.9],
    "1 month":  [-84.4, 68.6, 128.7, 191.6, 253.5, 315.8, 380.9, 443.6],
    "3 months": [-75.6, 72.8, 131.6, 193.3, 254.3, 315.7, 380.2, 442.2],
    "1 year":   [-69.4, 76.6, 134.2, 195.2, 255.3, 316.0, 379.6, 440.8],
}
bad = sum(1 for t in lane_mean for i in range(8) if lane_mean[t][i] != pdf_mean[t][i])
say(f"   means: {40-bad}/40 cells IDENTICAL to the source -- transcription CONFIRMED.")
say("   sigma rows: the PDF's Table 5 bottom block is COMPLETE (no truncation):")
pdf_sigma = {
    "1 day":    [48.2, 9.7, 9.7, 9.4, 9.3, 9.1, 9.5, 9.1],
    "1 week":   [46.4, 10.7, 10.8, 10.5, 10.6, 10.3, 10.6, 10.6],
    "1 month":  [46.8, 11.3, 11.2, 11.0, 10.9, 10.8, 11.2, 11.1],
    "3 months": [45.9, 12.0, 11.8, 11.5, 11.4, 11.4, 11.7, 11.7],
    "1 year":   [45.9, 12.8, 12.4, 12.0, 12.0, 11.9, 12.3, 12.4],
}
lane_sigma_1day = [48.2, 9.7, 9.7, 9.4, 9.3, 9.1, 9.3, 8.5]   # as in c86_nand_arm.py
say("   TRANSCRIPTION ERROR FOUND (lane sigma_1day, P6/P7): lane has 9.3, 8.5;")
say(f"   the PDF has {pdf_sigma['1 day'][6]}, {pdf_sigma['1 day'][7]}.  The wrong pair equals the tail of Table 4's")
say("   P/E=0 sigma row (…, 9.3, 8.5): a row-slip in a truncated text extraction.  The claim")
say("   that 1-year P6/P7 sigmas are 'NOT visible' is also false: they are published (12.3, 12.4).")
say("   IMPACT AUDIT: resolvability lines used only P1..P4 sigmas (correctly transcribed);")
say("   NO verdict moves.  The pinned range 'sigma 8.5-12.8' is wrong: programmed-state range")
say("   is 9.1-12.8.  Bonus: with correct values, Table 5's 1-day row equals Table 4's 2,000")
say("   P/E row in BOTH mean and sigma -- the 2,000 P/E identification gets STRONGER.")
g0 = pdf_mean["1 day"][2] - pdf_mean["1 day"][1]; s0 = pdf_sigma["1 day"][1] + pdf_sigma["1 day"][2]
g1 = pdf_mean["1 year"][2] - pdf_mean["1 year"][1]; s1 = pdf_sigma["1 year"][1] + pdf_sigma["1 year"][2]
say(f"   re-check P1-P2 resolvability: {g0/s0:.2f} -> {g1/s1:.2f} combined-sigma (lane said 3.2 -> 2.3): OK")
D = [pdf_mean["1 year"][i] - pdf_mean["1 day"][i] for i in range(8)]
steps = [D[i] - D[i+1] for i in range(7)]
say(f"   re-check P-ORD-1 on source values: decreasing steps {sum(s>0 for s in steps)}/7,")
say(f"   min step {min(steps):.1f} u, sign changes {sum(1 for i in range(7) if D[i]>=0>D[i+1])}: HOLDS as published.")

# ----------------------------------------------------- 2. semantics attacks on pin N-7
say("")
say("2. PIN N-7 SEMANTICS -- what the source does and does not say (re-fetched text)")
say("   CONFIRMED VERBATIM: 'absolute threshold voltage values are proprietary information to")
say("   flash vendors', normalized scale nominal max 512, 0 = GND.")
say("   NOT IN THE SOURCE: the Proc. IEEE 2017 paper NOWHERE states the retention temperature")
say("   ('room tem', 'C', 'Celsius', 'ambient': zero hits outside unrelated context).  The")
say("   'room temperature' semantics is imported from HPCA 2015 -- a DIFFERENT paper on")
say("   DIFFERENT chips (2y-nm MLC, not this TLC).  The pinned file discloses the import but")
say("   the dataset's own temperature is UNPINNED; every T-dependent number downstream")
say("   inherits this.  Also unpinned: whether T was controlled over the year (HotStorage")
say("   states 'air condition at 25 C'; Cai'17 states nothing).")

# ----------------------------------------------------- 3. the bake clock, re-derived
say("")
say("3. PIN N-6 / N-2 (bake clock) -- INDEPENDENT RE-DERIVATION from the re-fetched slides")
say("   Slide 5 (verbatim table): 60C 97.65 h | 80C 11.16 h | 100C 1.61 h | 120C 0.28 h,")
say("   titled 'Baking hours to obtain equivalent retention effects of 25C for 1 year'.")
T0 = 298.15
for TC, th in [(60, 97.65), (80, 11.16), (100, 1.61), (120, 0.28)]:
    T = TC + 273.15
    ea = math.log(8760.0/th) * K_EV / (1.0/T0 - 1.0/T)
    say(f"     {TC} C -> implied Ea = {ea:.4f} eV")
say("   CONFIRMS the lane's inversion (1.100 +- 0.001 eV).  JEDEC-side language CONFIRMED")
say("   first-hand from the JC-64.8 chair's deck (jedec.org): class table Client 30C/1yr/")
say("   FFR 3%/UBER 1e-15, Enterprise 40C/3mo/3%/1e-16 -- and, verbatim: retention 'is to be")
say("   verified both for a temperature-accelerated mechanism (1.1eV) and a non-temperature-")
say("   accelerated mechanism' -- the two-mechanism caveat is JEDEC's own, as pinned (N-3).")
say("   Malavena JEDS 2023 re-fetched: 'a single activation energy EA nearly equal to 1.1 eV'")
say("   and 'normalized to the same arbitrary constant' CONFIRMED VERBATIM; the abstract adds")
say("   that detrapping and depassivation have 'similar activation energy and comparable")
say("   magnitude on fresh devices' -- STRENGTHENS row C7's mechanism confounding (an")
say("   Arrhenius merge at 1.1 eV cannot by itself isolate one channel).")

# ----------------------------------------------------- 4. ATTACK: the f0 exclusion (C3/product i)
say("")
say("4. ATTACK SUSTAINED -- the f0 'exclusion' rests on TWO free identifications + unpinned T")
say("   (a) MAPPING SWITCH.  The model-side mapping paragraph declares 'each programmed cell")
say("       is one record'; section 4 then reads k(1yr) = 8 of 8 from class-level pair")
say("       non-merging.  Under the file's own cell-level mapping the census is over ~1e6")
say("       records and the SOURCE ITSELF says cells die inside the window: Cai'17 (re-fetch),")
say("       'the number of data retention and read disturb errors remains low at LOW retention")
say("       age' -- i.e. grows with age; pin N-9 measures per-voltage RBER rising over 366")
say("       days at 25 C.  Cell-level record deaths INSIDE the public year are a measured")
say("       fact; 'no drop within 1 yr' is true only of the 8-class re-mapping.")
say("   (b) Ea IDENTIFICATION.  Section 4 sets B = the pinned 1.1 eV; section 5 of the SAME")
say("       file says the pinned 1.1 eV 'was itself measured from Delta-VT transients, which")
say("       under the model mixes B - dE_i values'.  Both cannot hold.  If 1.1 eV = mixture")
say("       of (B - dE_i), then B = 1.1 + dE and the corner bound relaxes by e^(dE/kT):")
kT = K_EV * 293.15
for dE in (0.0, 0.05, 0.10, 0.20):
    f0max = math.exp((1.1 + dE)/kT) / (2.0 * YR)
    say(f"         dE = {dE:.2f} eV -> f0max = {f0max:.2e} Hz")
say("       The 'textbook decade 1e12-1e13 Hz excluded' claim dissolves for dE >~ 0.06-0.12 eV,")
say("       and dE in eV is exactly the unpinned calibration (alpha) the file elsewhere refuses")
say("       to value.  The free choice was resolved toward the stronger-sounding product.")
say("   (c) TEMPERATURE.  T = 293.15 K is imported (sec. 2 above).  Sensitivity, B = 1.1 eV:")
for TC in (20, 25, 30):
    kTx = K_EV * (TC + 273.15)
    say(f"         T = {TC} C -> f0max = {math.exp(1.1/kTx)/(2.0*YR):.2e} Hz")
say("   VERDICT: the C3 row's NOT-COMPARABLE stands, but the recorded 'f0 exclusion' is")
say("   REFUTED AS A PRODUCT.  Honest restatement: CONDITIONAL on (8-class mapping) AND")
say("   (Ea = B) AND (T = 20 C), survival of class corners to 1 yr needs f0 <= ~1e11 Hz.")
say("   Same defect infects section 5's 'under-delivers for EVERY displaced record': under")
say("   the mixture reading, records with dE below the mixture mean have true Ea > 1.1 eV")
say("   and the bake OVER-delivers for them.  C6 stays a shared null -- 'sign AGREES' is")
say("   conditional on Ea = B and should carry that flag.")

# ----------------------------------------------------- 5. P-ORD provenance
say("")
say("5. ATTACK PARTIAL -- P-ORD-1/2 provenance (rows C1/C2)")
say("   P-ORD-1's 'exactly one sign change' is guaranteed for ANY drift-toward-equilibrium")
say("   field once v_eq is free -- and v_eq is located (between P4 and P5) from the SAME")
say("   table.  Real content: magnitude monotone in |dE| -- which 7/7 decreasing steps do")
say("   test, above floor.  P-ORD-2's 'fastest programmed pair = P1-P2' needs CONVEXITY of")
say("   drift vs level on the up side, which the model file asserts ('farthest from the zero")
say("   crossing on the fast side') but never derives; measured up-side closure rates are")
say("   3.9 / 2.8 / 2.9 u -- P2-P3 vs P3-P4 differ by 0.1 u = the quantization floor, so the")
say("   data itself cannot certify convexity.  AGREES verdicts stand as data statements; as")
say("   MODEL predictions they are under-derived, and both files already concede they are")
say("   shared with Neel/Street-Woolley/Sharrock-class relaxation -- non-discriminating.")

# ----------------------------------------------------- 6. magnetic arm re-fetches
say("")
say("6. MAGNETIC ARM RE-FETCHES")
say("   M-1 Wernsdorfer: CONFIRMED VERBATIM (tau0 ~= 3e-9 s, E0 = 214 000 K, Hsw0 = 143.05 mT,")
say("   activated volume ~ (25 nm)^3 ~= SEM volume, P(t) exponential at every (T,H)).")
say("   M-2 Krause/Triannual: CONFIRMED VERBATIM ('nu0 is on the order of 10^13 and 10^16 Hz'")
say("   across the islands of one ensemble; E0 = 61+-5 meV; e_DW = 7.5+-0.4 meV/AR;")
say("   K = 0.55+-0.03 meV/atom; w = 2.15+-0.35 nm).  Re-derived stress numbers:")
kT536 = K_EV * 53.6
spread = math.log(1e16/1e13)
say(f"     dE_b per row / kT(53.6 K) = {7.5e-3/kT536:.2f}; ln nu0-spread = {spread:.2f};")
say(f"     order-scramble reach = {spread/(7.5e-3/kT536):.1f} rows; position floor e^{spread/2:.2f} = {math.exp(spread/2):.0f}x: all CONFIRMED.")
say("   M-3 Funatsu: physics pins CONFIRMED VERBATIM (D = 34 nm, TMR 73%, RA 5.5 Ohm um2,")
say("   tau_P/tau_AP separate vs H_z, exponential/Poisson verified, shortest relaxation")
say("   ~0.3 ms, 'attempt frequency tau0 of 1 ns was assumed').")
say("   *** AVAILABILITY CLAIM REFUTED ***: Zenodo 6767828 is RESTRICTED -- landing page:")
say("   'This datasets will be Open Access after the paper is published', 0 downloads, files")
say("   hidden behind login + access request (checked 2026-08-21, paper published 2022).")
say("   AND the paper's own statement scopes the deposit as 'data that support the plots',")
say("   not raw dwell-time sequences.  The arm's 'raw data public' / 'raw dwell times:")
say("   Zenodo 6767828' is wrong on access TODAY and unverified on content depth.  Named")
say("   next step (1) (the Zenodo grounding run) is NOT executable as written: it needs an")
say("   access request AND a content check.  Row C12's NOT-COMPARABLE (this run) stands.")
say("   M-4 Farhan PRL 111, 057204: existence, PEEM method, string->domain regimes confirmed")
say("   (APS/PubMed 23952441).  M-5..M-8 are CLASS pins at stated depth -- not re-fetched.")

# ----------------------------------------------------- 7. verdict table
say("")
say("=" * 98)
say("7. ADVERSARIAL VERDICT PER COMPARISON ROW (lane verdict -> adversarial disposition)")
say("=" * 98)
rows = [
 ("C1",  "AGREES (ordinal)",        "UPHELD; data re-verified cell-for-cell; non-discriminating (conceded class)"),
 ("C2",  "AGREES (ordinal)",        "UPHELD as data; model prediction under-derived (convexity unstated, v_eq from same table)"),
 ("C3",  "NOT-COMPARABLE (abs)",    "UPHELD; but the recorded f0-exclusion product REFUTED AS STATED (sec. 4 above)"),
 ("C4",  "AGREES (ordinal)",        "UPHELD; HPCA'15 re-fetched, 20 C + 'almost constant' P1 + anti-bake quote verbatim"),
 ("C5",  "NOT-COMPARABLE (owned)",  "UPHELD"),
 ("C6",  "AGREES (sign only)",      "UPHELD as shared null; 'EVERY displaced record' derivation wrong under mixture reading"),
 ("C7",  "NOT-COMPARABLE (premise)","UPHELD and STRENGTHENED (Malavena abstract: similar Ea, comparable magnitude)"),
 ("C8",  "NOT-COMPARABLE (sem.)",   "UPHELD; JEDEC class table re-verified from JEDEC's own deck"),
 ("C9",  "NOT-COMPARABLE (gate)",   "UPHELD; inversion re-derived independently, 1.100+-0.001 eV confirmed"),
 ("C10", "AGREES (form; conceded)", "UPHELD; Wernsdorfer numbers verbatim; circularity note correct"),
 ("C11", "NOT-COMPARABLE (scope)",  "UPHELD; nu0 spread verbatim in source; stress arithmetic re-derived"),
 ("C12", "NOT-COMPARABLE (run)",    "UPHELD; but 'raw data public' REFUTED -- Zenodo restricted, content depth unverified"),
 ("C13", "NOT-COMPARABLE (ext.)",   "UPHELD; Farhan class existence confirmed"),
 ("C14", "NOT-COMPARABLE (absent)", "UPHELD at CLASS depth (not re-fetched)"),
 ("C15", "NOT-COMPARABLE (none)",   "UPHELD; the departure term is touched by nothing located -- the operative gap"),
]
for r in rows:
    say(f"   {r[0]:>4}  lane: {r[1]:<26} adversarial: {r[2]}")
say("")
say("PROMOTION QUESTION (the strengthened bar):")
say("  * Rows satisfying 'measured external number beside the model's, agreeing within stated")
say("    tolerance' ON OWNED CONTENT: NONE.  The 5 AGREES rows are ordinal/sign/form-level and")
say("    each sits in territory the C-86 row concedes by name.")
say("  * Rows FALSIFYING C-86: NONE.  C7 breaks a mapping premise (second channel on fresh 3D),")
say("    which is scope, not law failure; per the two-way-null rule it cannot be read as a loss.")
say("  * Honest verdict: NOT-YET-COMPARABLE.  Missing data, named: (i) level-resolved margin-")
say("    free bake retirement times + absolute calibration (B_i, dE_i in eV, f0, volts) --")
say("    vendor-held per the sources' own statements; (ii) integer survivor census + remanence")
say("    on ONE magnetic sample with per-record constants -- not published anywhere located,")
say("    lab work; (iii) the SMTJ raw dwell data -- RESTRICTED on Zenodo, access request +")
say("    content check required before next-step (1) can run.")
say("  * FORMAL stands; FORMAL -> PROVED is NOT supported by this run.  Register corrections")
say("    owed: the sigma transcription slip (sec. 1), the unpinned N-7 temperature (sec. 2),")
say("    the conditionalized f0 constraint (sec. 4), the Zenodo access status (sec. 6).")
sys.exit(0)
